#!/usr/bin/python

# See below for name and description
# Copyright (C) 2020 Richard J. Edwards <dr.r.edwards@icloud.com>
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
Module:       synbad
Description:  Synteny-based scaffolding assessment and adjustment
Version:      0.12.2
Last Edit:    02/09/23
GitHub:       https://github.com/slimsuite/synbad
Copyright (C) 2020  Richard J. Edwards - See source code for GNU License Notice

Function:
    SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
    between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

    Synbad will use or create:

    1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
    spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
    correspond to misassemblies.

    2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

    Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
    by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
    should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

    Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
    sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
    new fields with 5', 3' or Both.

    The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
    corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
    If the local hits are against two different sequences from the other genome, the two sequence names will be
    entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
    small gaps), `SynSpan` or `GapSpan` will have a value of zero.

    Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values, and subsequent local hit
    analysis to fix/update ratings:

    * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
    * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
    * `Div` = `Divergence` = Translocation or Breakpoint gap with flanking collinear hits.
    * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
    * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
    * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
    * `Inv` = `Inversion` = Flanking hits are on alternative strands.
    * `InvBrk` = `Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.
    * `InvDupFix` = `Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)
    * `InvFix` = `Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)
    * `Long` = `Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).
    * `Null` = No mapping between genomes for that gap.
    * `Span` = `Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
    * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).
    * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
    * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.

Dependencies:
    SynBad needs Minimap2, samtools and kat installed for full functionality.
    For `gapass` gap mode, Flye also needs to be installed.
    To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genome1=FILE    : Genome assembly used as the query in the GABLAM/Mashmap searches []
    genome2=FILE    : Genome assembly used as the searchdb in the GABLAM/Mashmap searches []
    basefile=X      : Prefix for output files [synbad]
    gablam=X        : Optional prefix for GABLAM/Mashmap search [defaults to $BASEFILE.map]
    mapper=X        : Whether to use minimap2, busco or mashmap (dev only) for all-by-all mapping [minimap2]
    gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
    maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
    maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
    synreadspan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]
    checkflanks=LIST: List of lengths flanking gaps that must also be spanned by reads [0,100,1000]
    spannedflank=INT: Required flanking distance for synreadspan "Spanned" rating [0]
    maxoverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
    chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    ### ~ Correction and Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    correct=LIST    : List of edit types to try to fix in the assembly (invert/extract/relocate/break/join; T/True=all) [invert]
    fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic if no corrections made [False]
    fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
    minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
    minctglen=INT   : Extract any contigs below a minimum length threshold [500]
    minbadctg=INT   : Extract any contigs with bad flanking ratings below a minimum length threshold [5000]
    minscafflen=INT : Remove any scaffolds (inc. detached/extracted contigs) below minimum length threshold [500]
    gapsize=INT     : Size of gaps to add when relocating assembling chunks [500]
    rejoin=T/F      : Whether to rejoin original gaps that end up split into termini but not too short [True]
    ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    masked1=FASFILE : Optional masked fasta file for assembly comparison [$BASEFILE1.masked.fasta]
    bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    busco1=FILE     : Optional BUSCO full results file for genome 1 []
    genomesize1=INT : Haploid genome 1 size (bp) [0]
    scdepth1=NUM    : Single copy ("diploid") read depth for genome 1. If zero, will use SC BUSCO mode [0]
    masked2=FASFILE : Optional masked fasta file for assembly comparison [$BASEFILE2.masked.fasta]
    bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
    reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    busco2=FILE     : Optional BUSCO full results file for genome 2 []
    genomesize2=INT : Haploid genome 2 size (bp) [0]
    scdepth2=NUM    : Single copy ("diploid") read depth for genome 2. If zero, will use SC BUSCO mode [0]
    mapflanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
    mapflanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
    fullmap=T/F     : Whether to abort if not all flanks can be mapped [True]
    ### ~ HiC Gap Flank options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hicbam1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    hicbam2=FILE    : Optional BAM file of HiC reads mapped onto assembly 2 [$BASEFILE2.HiC.bam]
    gapflanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
    pureflanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
    hicscore=X      : HiC scoring mode (pairs/score/wtscore) [wtscore]
    hicmin=X        : Min. number of HiC read pairs for a "best" HiC pairing ruling [3]
    hicmode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
    hicdir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
    hicdir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
    ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    newacc1=X       : Scaffold name prefix for updated Genome 1 output [None]
    newacc2=X       : Scaffold name prefix for updated Genome 2 output [None]
    bestpair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
    update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]
    hidegaps=LIST   : List of SynBad gap types to "hide" in final outputs. Will need to be revealed again later []
    force=T/F       : Whether to force regeneration of SynBad results tables [False]
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_lrbridge, rje_mashmap, rje_rmd, rje_seqlist, rje_sequence
import diploidocus, gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version without fragment=T implementation.
    # 0.1.1 - Added minlocid=PERC and minloclen=INT.
    # 0.1.2 - Added additional translocation skipping for SynTrans rating.
    # 0.1.3 - Modified the SynBad classification text.
    # 0.1.4 - Modified code to be able to run without long read mapping. Added dochtml output.
    # 0.2.0 - Added fragment=T output.
    # 0.3.0 - Added chromosome scaffold Translocation restriction.
    # 0.4.0 - Added an Duplication rating in place of Breakpoint for overlapping flanking hits; added top sequence pairs.
    # 0.5.0 - Added HiFi read type. Changed default to gapmode=gapspan.
    # 0.5.1 - Added some extra bug-fixes for running Diploidocus checkpos.
    # 0.6.0 - Added additional gap-spanning classes and qry-hit pair compression. Added dev reworking on main compression.
    # 0.6.1 - Tidy up of code and transition of reworked code to main run (no longer dev=T only).
    # 0.7.0 - Added code for fixing Inversions and well-supported translocations. Updated defaults.
    # 0.8.0 - Added hicbam=FILE inputs and initial HiC assessment. Fixed inversion bug.
    # 0.8.1 - Adding hicbest and summary table output. Fixed some calculation and re-run bugs.
    # 0.8.2 - Separated and tidied HiC processing from contig flank/end processing. Fixed summary. Added correct=T/F option.
    # 0.8.3 - Fixed HIC flank mapping error.
    # 0.8.4 - Added simple duplicity analysis with KAT kmers and Diploidocus CNV. Added extract and relocate edits.
    # 0.8.5 - Replace swap edits with break and join edits -> only join if pair are both scaffold ends. Add rejoin=T/F. Fixed major gap update bug.
    # 0.8.6 - Small bug fixes for partial input.
    # 0.9.0 - Added hidegaps=LIST option for hiding gaps. Add MashMap in place of GABLAM (dev=True). Fixed naming clashes.
    # 0.10.0- Added mapper=busco option to use BUSCO genes in place of GABLAM.
    # 0.10.1- Fixed end of sequence gap bug for contig/flank generation.
    # 0.10.2- Fixed bug with new filenaming for Diploidocus wrapping of DepthKopy. (May need better fix.)
    # 0.11.0- Added alternative masked input files for the actual pairwise synteny comparisons.
    # 0.11.1- Py3 bug fixes.
    # 0.12.0- Added output of a QC map in Telociraptor Format.
    # 0.12.1- Fixed the implementation of checkflanks=LIST and spannedflank=INT settings.
    # 0.12.2- Fixed the correct=LIST bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add chr1 and chr2 chromosome prefix identifers (PAFScaff prefixes) to restrict Translocations to chromosomes.
    # [ ] : Add GapGFF function from rje_genomics with additional Gap rating information. (Or BED?)
    # [X] : Add fullcollapse=T/F option to collapse entire (pairwise) matching regions between gaps.
    # [Y] : Add tidying of gap ratings and suggested modifications to scaffolds.
    # [ ] : Make sure that Frag gap type is being handled appropriately.
    # [?] : Add runmode with different options.
    # [Y] : Add hicbam file loading and flank/join assessements.
    # [Y] : Add modified fasta output following edits.
    # [ ] : Build LocusFixer chassis and use for re-assembling duplication gap regions.
    # [Y] : Execute inversions and output updated fasta file.
    # [ ] : Empirically establish the best values for maxoverlap, maxsynskip and maxsynspan.
    # [ ] : Option to add DepthCharge gaps where there might be misassemblies.
    # [Y] : Option to add KAT kmer analysis of regions to help assess duplicity. (Can use with depth or alone.)
    # [Y] : Standardise the fields (CamelCase)
    # [?] : Switch early from Qry/Hit to SeqName/Hit and change output file names? (Or just do this at end.)
    # [ ] : Add final output where Qry and Hit are replaced for Genome1/2 output.
    # [ ] : Rationalise the final output to the really useful stuff.
    # [?] : Add minimum HiC read count and/or HiC score filter?
    # [Y] : Add hicbest table with the top pairs and mark as single or reciprocal.
    # [Y] : Separate the flank/contig generation from the HiC processing.
    # [Y] : Add HiCNone to summary table.
    # [Y] : Add kmerreads1/2 input and KAT read kmer frequency analysis.
    # [ ] : Rationalise so that: remake=T remakes qry, hit, gap and corrections tables and anything with gap rating flanks.
    # [ ] : - force remakes everything except re-running GABLAM and re-parsing BAM files.
    # [ ] : - fullforce for force-running GABLAM and BAM parsing etc.
    # [ ] : - update should update the qry/hit and gaps/corrections table but load the synbad (corrected) maps for more corrections.
    # [ ] : Need to make sure lack of HiC data will not cause SynBad to crash. Add some internal bool settings.
    # [ ] : Add option to add the tighter assembly versus assembly minimap2 settings and/or other mappings.
    # [ ] : Replace read mapping and SC depth analysis with rje_readcore.py. (Maybe via Diploidocus first.)
    # [ ] : Add parsing of Juicer merge_dups.txt for HiC mapping.
    # [ ] : Add running and parsing of MashMap to replace Minimap2 for initial mapping (2 x map filter)
    # [ ] : Add running of DepthKopy for read depth profiles in place of KAT?
    # [ ] : Add generation of ChromSyn output for (a) whole assemblies, and (b) chromzoom above X bp (1 Mbp?)
    # [ ] : Rather than gap masking, add option to output different gap sizes - can be used with the mingap=INT setting later?
    # [ ] : Add chromsyn output to SynBad.
    # [ ] : Consider adding TEL and CEN repeats to the synteny tables.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SynBad', '0.12.2', 'September 2024 ', '2020')
    description = 'Synteny-based scaffolding assessment and adjustment'
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
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
#i# Descriptions of SynBad gap types.
gapdesc = {'Aln':'`Aligned` = Gap is found in the middle of a local alignment to the Hit',
           'Syn':'`Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).',
           'Ins':'`Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.',
           'Brk':'`Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.',
           'Dup':'`Duplication` = Overlapping flanking hits on the same strand.',
           'Inv':'`Inversion` = Flanking hits are on alternative strands.',
           'Tran':'`Translocation` = `SynSpan` indicates matches are on different scaffolds.',
           'Frag':'`Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.',
           'Term':'`Terminal` = Gap is between a local alignment and the end of the query sequence.',
           'Span':'`Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.',
           'HiC':'`HiC-best` = Any translocation gaps that have reciprocal best-rated HiC support.',
           'DupHiC':'`HiC-best Duplication` = Any duplication gaps that have reciprocal best-rated HiC support.',
           'Null':'No mapping between genomes for that gap.',
           'Div':'`Divergence` = Translocation or Breakpoint gap with flanking collinear hits.',
           'InvBrk':'`Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.',
           'InvFix':'`Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)',
           'InvDupFix':'`Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)',
           'Long':'`Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).'}
#########################################################################################################################
#i# Setup lists of gap types for use in different situations
#i# For now, InvFix are not in goodgaps to avoid nested inversions etc. Recommend re-running on fixed fasta.
puregaps = ['Syn', 'Aln', 'Span', 'Long','HiC','DupHiC']
goodgaps = ['Syn', 'Aln', 'Span', 'Div', 'Ins', 'Long','HiC','DupHiC']
skipgaps = ['Syn', 'Aln', 'Span', 'Long', 'Div', 'Dup', 'Ins','DupHiC','HiC']  #X# 'InvFix', 'InvDupFix',
fraggaps = ['Brk', 'Inv', 'InvBrk', 'Frag', 'Tran']
hicgaps = ['Brk', 'InvBrk', 'Frag', 'Tran', 'Null']
#########################################################################################################################
## SYNBAD Database Table Details
#########################################################################################################################
#i# Genertic database field formatting dictionary for use with all tables except where noted
dbformats = {}
for field in ['SeqLen','Start','End','Span0','AlnNum','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd']: dbformats[field] = 'int'
for field in ['GapLen','MaxFlank5','MaxFlank3','Span0','Span100','Span1000','Span5000','GapSpan','HitStart','HitEnd','Non','Alt','Trans','Inv','Dup']: dbformats[field] = 'int'
for field in ['Score','CtgLen','ID1','ID2','Pairs','SeqBP','ReadBP','ModeX','katMed','kat5','kat3']: dbformats[field] = 'int'
for field in ['WTScore','Score','HiCScore','MeanX','CN','CIsyst','CIrand','katMean','CN5','CN3']: dbformats[field] = 'num'
#########################################################################################################################
#i# Database table fields and keys
dbfields = {'blocks':'Qry,QryStart,QryEnd,Hit,HitStart,HitEnd,Strand,AlnNum,Length,Identity,Non,QryGap,SynType,Flank5,Flank3'.split(','),
            'contigs':'SeqName,Start,End,CtgLen,Name,Flank5,Flank3,SynBad'.split(','),
            'corrections':'SeqName,Start,End,Flank1,Flank2,Edit,Details,SynBad'.split(','),
            'flanks':'SeqName,Start,End,Name,Score,Strand,Pair'.split(','),
            'flanks.checkcnv':'SeqName,Start,End,Name,Score,Strand,Pair,SeqBP,ReadBP,MeanX,ModeX,CN,CIsyst,CIrand'.split(','),
            'full':'Qry,QryStart,QryEnd,Hit,HitStart,HitEnd,Strand,AlnNum,Length,Identity,Non,QryGap,HitGap,Qry5,Qry3,Hit5,Hit3'.split(','),
            'gap':'SeqName,Start,End,SeqLen,GapLen,GapName,Span0,Span100,Span1000,Span5000,GapSpan,SynSpan,SynBad,GapFlank5,GapFlank3,HiCScore,BestFlank5,BestFlank3'.split(','),
            'hicbest':'Flank1,Flank2,ID1,ID2,Pairs,Type,Score,WTScore,Best'.split(','),
            'hicpairs':'Flank1,Flank2,ID1,ID2,Pairs,Type,Score,WTScore'.split(','),
            'main':'Qry,QryStart,QryEnd,Hit,HitStart,HitEnd,Strand,Length,Identity,Non,Alt,Trans,Inv,Dup,QryGap,HitGap,Hit5,Hit3'.split(','),
            'summary':'File,SeqNum,CtgNum,Good,Fix,Dupl,Bad,HiCBest,HiCPart,HiCWeak,HiCPoor,HiCNone,AltBest,AltPart,AltWeak,HiCScore,Aln,Brk,Div,Dup,DupHiC,Frag,HiC,Ins,Inv,InvBrk,InvDupFix,InvFix,Long,Null,Span,Syn,Term,Tran'.split(',')}
dbkeys = {'blocks':'Qry,QryStart,QryEnd'.split(','),
          'contigs':'SeqName,Start,End'.split(','),
          'corrections':'SeqName,Start,End,Edit'.split(','),
          'flanks':'SeqName,Start,End'.split(','),
          'flanks.checkcnv':'SeqName,Start,End'.split(','),
          'full':'Qry,QryStart,QryEnd'.split(','),
          'gap':'SeqName,Start,End'.split(','),
          'hicbest':'Flank1,Flank2'.split(','),
          'hicpairs':'Flank1,Flank2'.split(','),
          'main':'Qry,QryStart,QryEnd'.split(','),
          'summary':['File']}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SynBad Class                                                                                               #
#########################################################################################################################
class SynBad(rje_obj.RJE_Object):
    '''
    SynBad Class. Author: Rich Edwards (2020).

    Str:str
    - BAM1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    - BAM2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    - BUSCO1=FILE     : Optional BUSCO full results file for genome 1 []
    - BUSCO2=FILE     : Optional BUSCO full results file for genome 2 []
    - Chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    - Chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    - GABLAM=X        : Optional prefix for GABLAM search [defaults to $BASEFILE.map]
    - GapMode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    - Genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    - Genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    - HiCBAM1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    - HiCBAM2=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    - MapFlanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
    - MapFlanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
    - Masked1=FASFILE : Alternative masked input file for pairwise comparison []
    - Masked2=FASFILE : Alternative masked input file for pairwise comparison []
    - HiCScore=X      : HiC scoring mode (pairs/score/wtscore) [wtscore]
    - HiCMode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
    - HiCDir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
    - HiCDir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
    - NewAcc1=X       : Scaffold name prefix for updated Genome 1 output [None]
    - NewAcc2=X       : Scaffold name prefix for updated Genome 2 output [None]
    - PAF1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    - PAF2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]

    Bool:boolean
    - BestPair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
    - Correct=T/F     : Whether to try to fix errors in the assembly [True]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
    - FullMap=T/F     : Whether to abort if not all flanks can be mapped [True]
    - Mapper=X        : Whether to use mashmap or minimap2 for all-by-all mapping [minimap2]
    - PureFlanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
    - Rejoin=T/F      : Whether to rejoin original gaps that end up split into termini but not too short [True]
    - Update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]

    Int:integer
    - GapFlanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
    - GapSize=INT     : Size of gaps to add when relocating assembling chunks [500]
    - GenomeSize1=INT : Haploid genome 2 size (bp) [0]
    - GenomeSize2=INT : Haploid genome 2 size (bp) [0]
    - HiCMin=X        : Min. number of HiC read pairs for a "best" HiC pairing ruling [3]
    - MaxOverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
    - MaxSynSkip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
    - MaxSynSpan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
    - MinBadCtg=INT   : Extract any contigs with bad flanking ratings below a minimum length threshold [5000]
    - MinCtgLen=INT   : Extract any contigs below a minimum length threshold [500]
    - MinScaffLen=INT : Remove any scaffolds (inc. detached/extracted contigs) below minimum length threshold [500]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    - MinReadSpan=INT     : Min number of Span0 reads in gaps table to prevent fragmentation [1]
    - SpannedFlank=INT : Required flanking distance for synreadspan "Spanned" rating [0]
    - SynReadSpan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]

    Num:float
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
    - SCDepth1=NUM    : Optional single copy read depth for genome 1 []
    - SCDepth2=NUM    : Optional single copy read depth for genome 2 []

    File:file handles with matching str filenames

    List:list
    - CheckFlanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000]
    - Correct=LIST    : List of edit types to try to fix in the assembly (invert/extract/relocate; T/True=all) [invert]
    - FragTypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
    - HideGaps=LIST   : List of SynBad gap types to "hide" in final outputs. Will need to be revealed again later []
    - Reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - Reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - ReadType2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - qry = Assembly Map list for Genome 1
    - hit = Assembly Map list for Genome 2

    Dict:dictionary
    - FlankMap = {qry:{GapFlank:OldFlank},hit:{GapFlank:OldFlank}}
    - HiCOffProb = {qh:BadHiCProb}

    Obj:RJE_Objects
    - DB = Database object
    - Genome1 = Genome 1 SeqList Object
    - Genome2 = Genome 2 SeqList Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM1','BAM2','BUSCO1','BUSCO2','Chr1','Chr2','GABLAM','GapMode','Genome1','Genome2','GenomeSize1','GenomeSize2','HiCBAM1','HiCBAM2',
                        'MapFlanks1','MapFlanks2','Masked1','Masked2','HiCScore','HiCMode','HiCDir1','HiCDir2','NewAcc1','NewAcc2','PAF1','PAF2','Mapper']
        self.boollist = ['BestPair','Correct','DocHTML','Fragment','FullMap','Rejoin','Update']
        self.intlist = ['GapFlanks','GapSize','GenomeSize1','GenomeSize2','HiCMin','MaxOverlap','MaxSynSkip','MaxSynSpan','MinBadCtg','MinCtgLen','MinScaffLen','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan']
        self.numlist = ['MinLocID','SCDepth1','SCDepth2']
        self.filelist = []
        self.listlist = ['Correct','FragTypes','HideGaps','Reads1','Reads2','ReadType1','ReadType2','qry','hit']
        self.dictlist = ['FlankMap','HiCOffProb']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GapMode':'gapspan','HiCScore':'wtscore','HiCMode':'synbad','Mapper':'minimap2'})
        self.setBool({'BestPair':False,'Correct':True,'Fragment':False,'PureFlanks':True,'FullMap':True,'Rejoin':True,'Update':False})
        self.setInt({'GapFlanks':10000,'GapSize':500,'HiCMin':3,'MinBadCtg':5000,'MinCtgLen':500,'MinScaffLen':500,'MaxOverlap':500,'MaxSynSkip':4,'MaxSynSpan':25000,'MinLocLen':1000,'MinReadSpan':1,'SpannedFlank':0,'SynReadSpan':5})
        self.setNum({'MinLocID':50.0,'SCDepth1':0.0,'SCDepth2':0.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['Correct'] = ['invert'] #,'extract','relocate','break','join']
        self.list['CheckFlanks'] = [0,100,1000]
        self.list['FragTypes'] = fraggaps
        self.dict['FlankMap'] = {'qry':{}, 'hit':{}}
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
                self._cmdReadList(cmd,'str',['Chr1','Chr2','GABLAM','GapMode','HiCScore','HiCMode','GenomeSize1','GenomeSize2','NewAcc1','NewAcc2','Mapper'])   # Normal strings
                self._cmdReadList(cmd,'path',['HiCDir1','HiCDir2'])  # String representing directory path
                self._cmdReadList(cmd,'file',['BAM1','BAM2','BUSCO1','BUSCO2','Genome1','Genome2','HiCBAM1','HiCBAM2','MapFlanks1','MapFlanks2','Masked1','Masked2','PAF1','PAF2'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BestPair','DocHTML','Fragment','FullMap','PureFlanks','Rejoin','Update'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['GapFlanks','GapSize','HiCMin','MinBadCtg','MinCtgLen','MinScaffLen','MaxOverlap','MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan'])   # Integers
                self._cmdReadList(cmd,'perc',['MinLocID']) # Percentage
                self._cmdReadList(cmd,'num',['SCDepth1','SCDepth2']) # Percentage
                self._cmdRead(cmd,type='num',att='SCDepth1',arg='scdep1')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='num',att='SCDepth2',arg='scdep2')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='GenomeSize1',arg='gensize1')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='GenomeSize2',arg='gensize2')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['FragTypes','HideGaps','ReadType1','ReadType2'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['CheckFlanks'])  # List of integers (split on commas or file lines)
                self._cmdReadList(cmd,'lclist',['Correct'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Reads1','Reads2']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStrLC('HiCScore') == 'pairs': self.setStr({'HiCScore':'Pairs'})
        elif self.getStrLC('HiCScore') == 'score': self.setStr({'HiCScore':'Score'})
        elif self.getStrLC('HiCScore') == 'wtscore': self.setStr({'HiCScore':'WTScore'})
        else:
            self.printLog('#HIC','HiCScore={0} not recognised: set to WTScore')
            self.setStr({'HiCScore':'WTScore'})
        for att in ['GenomeSize1','GenomeSize2']:
            if self.getStrLC(att):
                try: self.setInt({att:rje_seqlist.bpFromStr(self.getStrLC(att))})
                except:
                    self.errorLog('Problem with {0}={1}. Setting {0}=0'.format(att.lower(),self.getStrLC(att)))
                    self.setInt({att:0})
        self.setBool({'Correct':self.list['Correct'][0] not in ['','none','f','false']})
        if 'inv' in self.list['Correct']: self.list['Correct'].append('invert')
        self.setStr({'Mapper':self.getStrLC('Mapper')})
        if self.getStrLC('Mapper') not in ['minimap2','blastn','mashmap','busco']:
            self.warnLog('Mapper="{0}" not recognised: defaulting to mapper=minimap2'.format(self.getStrLC('Mapper')))
            self.setStr({'Mapper': 'minimap2'})
        # Adjust the flanking list
        if 0 not in self.list['CheckFlanks']:
            self.list['CheckFlanks'].append(0)
        if self.getInt('SpannedFlank') not in self.list['CheckFlanks']:
            self.list['CheckFlanks'].append(self.getInt('SpannedFlank'))
            self.printLog('#SPAN','Added spannedflank={0} length to checkflanks=LIST'.format(self.getInt('SpannedFlank')))
        self.list['CheckFlanks'].sort()
        while self.list['CheckFlanks'][0] < 0: self.list['CheckFlanks'] = self.list['CheckFlanks'][1:]
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # SynBad:  Synteny-based scaffolding adjustment

        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

        SynBad will use or create:

        1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
        spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
        correspond to misassemblies.

        2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

        Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
        by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
        should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

        Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
        sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
        new fields with the number of flanking 5', 3' gaps.

        The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
        corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
        If the local hits are against two different sequences from the other genome, the two sequence names will be
        entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
        small gaps), `SynSpan` or `GapSpan` will have a value of zero.

        Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:

        * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
        * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).
        * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
        * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
        * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
        * `Inv` = `Inversion` = Flanking hits are on alternative strands.
        * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
        * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
        * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.
        * `Span` = `Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
        * `Null` = No mapping between genomes for that gap.

        If `chr1=X` and/or `chr2=X` chromosome scaffold prefixes are provided then `Translocation` will be restricted to
        matches between two different chromosome scaffolds. Matches including one or more non-chromosome scaffolds will
        be classed as `Fragmentation`.

        Following gap fixes for inversions and translocations/breakpoints (below), following

        After the initial mapping and GABLAM unique hit reduction, SynBad will further compress the Qry-Hit pairs by
        combining collinear hits between Qry-Hit pairs into larger blocks. Merged hits must be within `maxsynspan=INT`
        bp, overlap by no more than `maxoverlap=INT` bp, and have no intervening assembly gaps. An initial pass compares
        adjacent hits in position order. This is followed by iterative compressions where up to `maxsynskip=INT`
        intervening hits will be incorporated, providing there are no assembly gaps between the merged hits. This will
        start with the longest local hit, scanning downstream and then upstream for hits to merge, with increasing numbers
        of intervening hits up to `maxsynskip=INT`. Once no more merging is found for that hit, the next longest hit will
        be considered. If any hits are merged during this process, it will be repeated until convergence is reached.

        During merging, `Length` and `Identity` statistics will be summed, and additional statistics will be added:

        * `Non` = The unaligned bases in the query between the GABLAM local hits.
        * `Alt` = The number of bases (`Length`) of hits to different hit scaffolds
        * `Trans` = The number of bases (`Length`) of hits to the same scaffold but outside the merged hits.
        * `Inv` = The number of bases (`Length`) of intermediate hits that are in the right place but inverted strand.
        * `Dup` = The number of overlapping bases in the query scaffolds between adjacent hits.

        NOTE: This will be done for both genomes. For Genome 2, the "hits" are the `Qry` scaffolds, and all assessments
        of collinearity will be performed on the `Hit` scaffolds.

        Next, SynBad will assess inversions and apparent translocations for flanking collinear hits and instigate a
        number of updated ratings or suggested fixes. The following new SynBad gap ratings may be added:

        * `Div` = `Divergence` = Translocation or Breakpoint gap with flanking collinear hits.
        * `InvBrk` = `Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.
        * `InvFix` = `Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)
        * `InvDupFix` = `Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)
        * `Long` = `Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).

        If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
        `minreadspan=INT` reads span the gap.

        A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
        where possible. [GapSpanner](https://github.com/slimsuite/gapspanner) can also be used to gap-fill spanned gaps
        prior to running SynBad.

        ## Dependencies

        SynBad needs Minimap2 installed. For `gapass` gap mode, Flye also needs to be installed. For KAT kmer assessments
        of flanks, KAT must be installed (or pre-run on the genomes).
        To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.


        ## Commandline options

        ```
        ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
        genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
        basefile=X      : Prefix for output files [synbad]
        gablam=X        : Optional prefix for GABLAM search [defaults to $BASEFILE.map]
        gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
        minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
        maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
        maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
        synreadspan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]
        checkflanks=LIST: List of lengths flanking gaps that must also be spanned by reads [0,100,1000]
        spannedflank=INT: Required flanking distance for synreadspan "Spanned" rating [0]
        maxoverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
        chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
        chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
        ### ~ Correction and Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        correct=LIST    : List of edit types to try to fix in the assembly (invert/extract/relocate/break/join; T/True=all) [True]
        fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic if no corrections made [False]
        fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
        minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
        minctglen=INT   : Extract any contigs below a minimum length threshold [500]
        minbadctg=INT   : Extract any contigs with bad flanking ratings below a minimum length threshold [5000]
        minscafflen=INT : Remove any scaffolds (inc. detached/extracted contigs) below minimum length threshold [500]
        gapsize=INT     : Size of gaps to add when relocating assembling chunks [500]
        rejoin=T/F      : Whether to rejoin original gaps that end up split into termini but not too short [True]
        ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        masked1=FASFILE : Optional masked fasta file for assembly comparison [$BASEFILE1.masked.fasta]
        bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
        paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
        reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        busco1=FILE     : Optional BUSCO full results file for genome 1 []
        scdep1=NUM      : Optional single copy read depth for genome 1 []
        masked2=FASFILE : Optional masked fasta file for assembly comparison [$BASEFILE2.masked.fasta]
        bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
        paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
        reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        busco2=FILE     : Optional BUSCO full results file for genome 2 []
        scdep2=NUM      : Optional single copy read depth for genome 2 []
        mapflanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
        mapflanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
        fullmap=T/F     : Whether to abort if not all flanks can be mapped [True]
        ### ~ HiC Gap Flank options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        hicbam1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
        hicbam2=FILE    : Optional BAM file of HiC reads mapped onto assembly 2 [$BASEFILE2.HiC.bam]
        gapflanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
        pureflanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
        hicscore=X      : HiC scoring mode (pairs/score/wtscore) [wtscore]
        hicmode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
        hicmin=X        : Min. number of HiC read pairs for a "best" HiC pairing ruling [3]
        hicdir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
        hicdir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
        ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        newacc1=X       : Scaffold name prefix for updated Genome 1 output [None]
        newacc2=X       : Scaffold name prefix for updated Genome 2 output [None]
        bestpair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
        update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]
        force=T/F       : Whether to force regeneration of SynBad results tables [False]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        ## SynBad workflow

        _Details of the SynBad workflow will be added in a future release._

        **NOTE:** SynBad is still under development. Features and details will continue to be updated.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return self.docHTML()
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.synBad()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if self.getStrLC('GapMode') not in ['gapspan','gapass','gapfill']:
                self.printLog('#MODE','GapMode "{0}" not recognised. Setting gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                self.setStr({'GapMode':'gapspan'})
            if self.getStrLC('GapMode') not in ['gapspan','gapass']:
                self.printLog('#MODE','GapMode "{0}" no longer supported. Please run Diploidocus separately or use gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                if rje.yesNo('Set gapmode=gapspan?'):
                    self.setStr({'GapMode':'gapspan'})
                    self.printLog('#MODE','Setting gapmode=gapspan.')
                else: return False
            ### ~ [2] Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('Genome1')): raise IOError('Genome1 file "{0}" not found!'.format(self.getStr('Genome1')))
            if not rje.exists(self.getStr('Genome2')): raise IOError('Genome2 file "{0}" not found!'.format(self.getStr('Genome2')))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def docHTML(self):  ### Generate the Diploidocus Rmd and HTML documents.                                        # v0.1.0
        '''Generate the Diploidocus Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            rtxt += rje.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2020 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown SynBad documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### Loading and formatting data methods                                                                     #
#########################################################################################################################
    def update(self): return self.force() or self.getBool('Update')
#########################################################################################################################
    def addTable(self,tclass,ttype='',expect=False,replace=False,make=True): ### Loads/creates and returns tclass.ttype table, else None/Error
        '''
        Loads/creates and returns tclass.ttype table, else None/Error.
        >> tclass:str = 'qry/hit/main'
        >> ttype:str = table subtype (e.g. gap or corrections)
        >> expect:bool [False] = whether to raise an error if table not found (True) or create empty table (False)
        >> replace:bool [False] = whether to delete existing table if found (True) and then load/create new one (False).
        >> make:bool [True] = whether to create an empty table if missing (True) or return None (False)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tclass = tclass.lower()
            tname = '{0}.{1}'.format(tclass,ttype.lower())
            if not ttype: tname = tclass
            db = self.db()
            dkey = ttype
            if not ttype: dkey = tclass
            if dkey in ['qry','hit']: dkey = 'main'
            table = self.dbTable(tclass,ttype,expect=False)
            if table and replace: db.deleteTable(table)
            elif table: return table
            mainkeys = dbkeys[dkey]
            if ttype in ['full','','blocks'] and tclass == 'hit':
                mainkeys = ['Hit','HitStart','HitEnd']
            table = db.addTable(mainkeys=mainkeys, name=tname, ignore=[], expect=expect)
            if table:
                table.dataFormat(dbformats)
            if not table and make:
                table = db.addEmptyTable(tname,dbfields[dkey],dbkeys[dkey],log=True)
            #self.deBug(', '.join(db.tableNames()))
            return table
        except:
            raise
#########################################################################################################################
    def dbTable(self,tclass,ttype='',expect=True,load=True):     ### Returns tclass.ttype table, else None/Error
        '''
        Returns tclass.ttype table, else None/Error.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tname = '{0}.{1}'.format(tclass.lower(),ttype.lower())
            if not ttype: tname = tclass.lower()
            table = self.db(tname)
            if table: return table
            elif expect:
                table = self.addTable(tclass,ttype,expect=True,replace=False,make=False)
                if table: return table
                raise ValueError('Cannot find table: "{0}"'.format(tname))
            return None
        except:
            raise
#########################################################################################################################
    def seqObjSetup(self):  ### Loads the two genomes into sequence list objects
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:  ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList1' not in self.obj or not self.obj['SeqList1']:
                base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
                cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base1)]
                seqcmd = self.cmd_list + cmd1 + ['summarise=F', 'gapstats=F', 'raw=F', 'dna=T']
                self.obj['SeqList1'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList2' not in self.obj or not self.obj['SeqList2']:
                base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
                cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base2)]
                seqcmd = self.cmd_list + cmd2 + ['summarise=F', 'gapstats=F', 'raw=F', 'dna=T']
                self.obj['SeqList2'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            return
        except:
            self.errorLog('%s.seqObjSetup error' % self.prog())
#########################################################################################################################
    ### <4> ### Main SynBad Workflow                                                                                    #
#########################################################################################################################
    #i# SynBad processing workflow:
    #i# [1] Load genomes = self.diploidocusGapRun()
    #i# [2] Extract gaps = self.diploidocusGapRun()
    #i# [3] GABLAM Search of Qry versus Hit = self.runGABLAM()
    #i# [4] Use the genome mapping to rate all the gaps = self.synBadMapping(cdb1,cdb2)
    #i# [5] Generate flanks and contigs and output assembly maps = self.contigMapping()
    #i# [6] Generate (or map) flanks and map on the HiC data if available = self.synBadHiCMap()
    #i# [7] Add Gaps to Qry/Hit table and compress regions = self.synBadCompress()
    #i# [8] Analyse gaps for reclassifications and possible fixes = self.tidyTables()
    #i# [9] Save primary results files pre-correction = self.saveTables()
    #i# [10] Summarise gap ratings = self.synBadSummarise()
    #i# [11] Fix assembly where possible = self.mapCorrections()
    #i# [12] Fragment assembly on dodgy gaps (optional) = self.mapFragment()
#########################################################################################################################
    def synBad(self):   ### Main SynBad run method
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#
            #!# Add possibility to run without a second genome - pure read-spanning and HiC verification

            ### ~ [1] Run GapSpanner on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Replace this to use GapSpanner code rather than Diploidocus #!#
            self.headLog('GAPSPANNER GAP ANALYSIS', line='=')
            self.infoLog('SynBad uses GapSpanner to extract a table of assembly gaps from each input assembly.')
            self.infoLog('If long-read sequencing data is provided, gaps will also be assessed for spanning read support.')
            (cdb1,cdb2) = self.runGapSpanner()
            if not cdb1 or not cdb2: return False

            ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GENOME VS GENOME GABLAM/MASHMAP SEARCH',line='=')
            self.infoLog('GABLAM or MashMap is used to search Genome1 against Genome2. Local alignments are filtered by length and identity.')
            self.infoLog('For each genome, local hits are reduced to unique (non-overlapping) alignments with the other genome.')
            self.infoLog('This unique hit reduction is performed based on number of identical aligned bases.')
            self.infoLog('Partially overlapping local hits are trimmed.')
            if not self.runGABLAM(): return False

            ### ~ [3] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP MAPPING',line='=')
            infotxt = ['SynBad processes the GABLAM tables, adding the GapSpanner gaps.',
                       'Local alignments between genomes are used to establish gap classes based on hit distance and orientation.',
                       'For each genome, the local hit tables are sorted and `QryGap` and `HitGap` fields added.',
                       'Any local alignments with flanking gaps are flagged and the gap tables updated with `GapSpan` and `SynSpan` fields.',
                       'See docs for details of SynBad classification.']
            for txt in infotxt: self.infoLog(txt)
            if not self.synBadMapping(cdb1,cdb2): return False

            ### ~ [4] Generate contigs, CtgEnds and initial assembly map based on raw SynBad ~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('CONTIG ENDS AND ASSEMBLY MAPPING',line='=')
            infotxt = ['Terminal regions from Contigs are extracted as `Flank` sequences for assembly maps and HiC mapping.',
                       'This generates the *.flanks.*, *.contigs.* and *.map.* outputs.']
            for txt in infotxt: self.infoLog(txt)
            if not self.contigMapping(): return False
            if not self.mapNewFlanks(): return False
            if self.debugging() and not self.dataCheck(): return False

            ### ~ [5] HiC read mapping and flank pairing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('HIC READ PAIR FLANK MAPPING AND FLANK PAIR SCORING',line='=')
            infotxt = ['Mapped HiC read pairs are reduced to Contig End flank regions and then filtered to pairs in different Flank regions.',
                       'Pairs of Flank regions (with a focus on Gap Flanks) are compared and scored for shared read pairs.',
                       'This generates the *.hicpairs.tdt and *.hicbest.tdt output.']
            for txt in infotxt: self.infoLog(txt)
            if not self.synBadHiCMap(): return False
            if self.debugging() and not self.dataCheck(): return False

            ### ~ [6] Add Gaps to Qry/Hit table and compress regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD HIT COMPRESSION',line='=')
            infotxt = ['Local hits are compressed to ungapped collinear regions using the MaxSynSpan, MaxOverlap and MaxSynSkip settings.',
                       'Gaps are also added to the main *.qry.tdt and *.hit.tdt tables.',
                       'This generates the *.full.tdt output.']
            for txt in infotxt: self.infoLog(txt)
            if not self.synBadCompress(): return False
            if self.debugging() and not self.dataCheck(): return False

            ### ~ [7] Analyse gaps for reclassifications and possible fixes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP & SCAFFOLD TIDY PROCESSING',line='=')
            infotxt = ['Compressed local hits are scanned for gaps with additional support and/or regions between pairs of gaps that can be fixed.',
                       'Tran/Brk gaps with sufficient flanking support are re-classified as "Div" (Divergent).',
                       'Brk gaps with adjacent flanks exceeding maxsynspan bp are re-classified as "Long".',
                       'Brk/InvBrk/Frag/Tran/Null gaps with Best HiC pair support are re-classified as "HiC".',
                       'Dup gaps with Best HiC pair support are re-classified as "DupHiC".',
                       'Inversions and Duplicate Inversions are identified and marked for fixing if possible',
                       'Hit regions are further compressed to synteny blocks using the updated Gap classification.']
            for txt in infotxt: self.infoLog(txt)
            if not self.tidyTables(): return False
            if self.debugging() and not self.dataCheck(): return False

            ### ~ [8] Output and summarise primary results pre-correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD OUTPUT',line='=')
            if not self.saveTables(backup=False): return False
            if not self.synBadSummarise(): return False
            for qh in ['qry','hit']:
                if not self.saveTelociraptorMaps(qh,mapname='gap'): return False

            ### ~ [9] Update assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD ASSEMBLY CORRECTION',line='=')
            # This generates the *.synbad.txt and *.synbad.fasta corrected assembly outputs.
            if self.getBool('Correct'):
                if not self.mapCorrections(): return False
                if self.debugging() and not self.dataCheck(): return False
            else:
                self.printLog('#EDIT','No assembly correction edits (correct=F).')
                self.infoLog('See *.corrections.tdt tables for recommended edits.')

            ### ~ [10] Fragment Gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP FRAGMENTATION',line='=')
            # If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
            # `minreadspan=INT` reads span the gap.
            if self.getBool('Fragment'):
                return self.mapFragment()
            else:
                self.printLog('#FRAG','No fragmentation (fragment=F).')

            return True
        except: self.errorLog('%s.synBad error' % self.prog()); return False
#########################################################################################################################
    ### <5> ### GapSpanner Gap Generation and Read Spanning                                                             #
#########################################################################################################################
    ### ~ [1] Run GapSpanner on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    def runGapSpanner(self):    ### Runs GapSpanner (Diploidocus) and returns the gap tables.
        '''
        Runs GapSpanner (Diploidocus) and returns the gap tables.
        :return: (cdb1,cdb2)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            dcmd = ['bam=', 'paf=', 'reads=', 'readtype=ont']
            wanted = ['qry.gap.tdt','hit.gap.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad *.gap.tdt files found: skipping GapSpanner extraction (force=F)')
                cdb1 = db.addTable('%s.qry.gap.tdt' % db.baseFile(), mainkeys=['SeqName', 'Start', 'End'], name='qry.gap', ignore=[], expect=True)
                cdb1.dataFormat(dbformats)
                cdb2 = db.addTable('%s.hit.gap.tdt' % db.baseFile(), mainkeys=['SeqName', 'Start', 'End'], name='hit.gap', ignore=[], expect=True)
                cdb2.dataFormat(dbformats)
                return (cdb1,cdb2)
            # if not self.force():
            #     self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
            #     self.setBool({'Force':True})

            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#
            #!# Update this to cycle qry and hit as for other methods

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # SynBad will use or create:
            #
            # 1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
            # spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
            # correspond to misassemblies.
            #
            # ==> tigersnake.v2.7.pafscaff.checkpos.tdt <==
            # #       seqname start   end     seqlen  gaplen  gap     MaxFlank5       MaxFlank3       Span0   Span100 Span1000        Span5000
            # 209     NSCUCHR1.01     57638   57737   341225390       100     NSCUCHR1.01.57638-57737 57637   341167653       1       1       1       0
            # 2       NSCUCHR1.01     103648  104147  341225390       500     NSCUCHR1.01.103648-104147       103647  341121243       14      14      12      5

            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base1)]
            rundip = False
            if self.getStrLC('BAM1'): cmd1.append('bam={0}'.format(self.getStr('BAM1')))
            if not self.getStrLC('PAF1') and rje.exists('{0}.paf'.format(base1)): self.setStr({'PAF1':'{0}.paf'.format(base1)})
            if self.getStrLC('PAF1'): cmd1.append('paf={0}'.format(self.getStr('PAF1'))); rundip = True
            if self.list['Reads1']: cmd1.append('reads={0}'.format(','.join(self.list['Reads1']))); rundip = True
            if self.list['ReadType1']: cmd1.append('readtype={0}'.format(','.join(self.list['ReadType1'])))
            dip1 = diploidocus.Diploidocus(self.log, self.cmd_list + dcmd + cmd1)
            cdb1 = None
            if rundip:
                dip1.run()
                cdb1 = dip1.db('checkpos')
                if cdb1:
                    cdb1.baseFile(basefile)
                    cdb1.setStr({'Name': 'qry.gap'})
                    self.db().list['Tables'].append(cdb1)
                else:
                    cdb1 = db.addTable('%s.checkpos.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
            if not cdb1:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base1))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base1, log=self.log):
                    seqcmd = self.cmd_list + cmd1 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList1'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                cdb1 = db.addTable('%s.gaps.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
                cdb1.addField('Span0', evalue=0)
            cdb1.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb1.addField('GapSpan', evalue='.')
            cdb1.addField('SynSpan', evalue='.')
            cdb1.addField('SynBad', evalue='Null')
            cdb1.index('seqname')
            if '#' in cdb1.fields(): cdb1.dropField('#')
            self.updateGapFields(cdb1)

            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base2)]
            rundip = False
            if self.getStrLC('BAM2'): cmd2.append('bam={0}'.format(self.getStr('BAM2')))
            if self.getStrLC('PAF2'): cmd2.append('paf={0}'.format(self.getStr('PAF2'))); rundip = True
            if self.list['Reads2']: cmd2.append('reads={0}'.format(','.join(self.list['Reads2']))); rundip = True
            if self.list['ReadType2']: cmd2.append('readtype={0}'.format(','.join(self.list['ReadType2'])))
            dip2 = diploidocus.Diploidocus(self.log, self.cmd_list + dcmd + cmd2)
            cdb2 = None
            if rundip:
                dip2.run()
                cdb2 = dip2.db('checkpos')
                if cdb2:
                    cdb2.baseFile(basefile)
                    cdb2.setStr({'Name': 'hit.gap'})
                    self.db().list['Tables'].append(cdb2)
                else:
                    cdb2 = db.addTable('%s.checkpos.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hit.gap', ignore=[], expect=True)
            if not cdb2:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base2))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base2, log=self.log):
                    seqcmd = self.cmd_list + cmd2 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList2'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                cdb2 = db.addTable('%s.gaps.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hit.gap', ignore=[], expect=True)
                cdb2.addField('Span0', evalue=0)
            cdb2.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb2.addField('GapSpan', evalue='.')
            cdb2.addField('SynSpan', evalue='.')
            cdb2.addField('SynBad', evalue='Null')
            cdb2.index('seqname')
            if '#' in cdb2.fields(): cdb2.dropField('#')
            self.updateGapFields(cdb2)

            return(cdb1,cdb2)
        except: self.errorLog('%s.runGapSpanner error' % self.prog()); raise
#########################################################################################################################
    def updateGapFields(self,gapdb):    ### Converts lower case fields to CamelCase for SynBad output
        '''
        Converts lower case fields to CamelCase for SynBad output.
        '''
        fieldconversion = {'seqname':'SeqName','start':'Start','end':'End','seqlen':'SeqLen','gaplen':'GapLen','gap':'GapName'}
        for field in ['seqname','start','end','seqlen','gaplen','gap']:
            if field in gapdb.fields(): gapdb.renameField(field,fieldconversion[field])
#########################################################################################################################
    ### <6> ### GABLAM Assembly Mapping                                                                                 #
#########################################################################################################################
    ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    def runGABLAM(self):   ### Runs GABLAM of two assemblies
        '''
        Runs GABLAM of two assemblies and loads results into database tables.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            wanted = ['qry.tdt','hit.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt and hit.tdt files found: skipping GABLAM (force=F)')
                return True
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#
            ## ~ [0a] Alternative masked fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fasta1 = self.getStr('Genome1')
            fasta2 = self.getStr('Genome2')
            masked1 = self.getStr('Masked1')
            masked2 = self.getStr('Masked2')
            if self.getStrLC('Masked1'):
                rje.checkForFiles([masked1],log=self.log,ioerror=True)
                fasta1 = masked1
            elif rje.checkForFiles([rje.baseFile(fasta1)+'.masked.fasta'],log=self.log,ioerror=False):
                fasta1 = rje.baseFile(fasta1)+'.masked.fasta'
            if self.getStrLC('Masked2'):
                rje.checkForFiles([masked2],log=self.log,ioerror=True)
                fasta2 = masked2
            elif rje.checkForFiles([rje.baseFile(fasta2)+'.masked.fasta'],log=self.log,ioerror=False):
                fasta2 = rje.baseFile(fasta2)+'.masked.fasta'

            ### ~ [1] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Mapper') == 'mashmap':
                if self.dev(): return self.runMashMap()
                else:
                    self.warnLog('Reverted to mapper=minimap2 (dev=F)')
                    self.setStr({'Mapper':'minimap2'})
            gabbase = basefile + '.map'
            if self.getStrLC('GABLAM'): gabbase = self.getStr('GABLAM')
            self.printLog('#GABLAM','GABLAM output basefile: {0}'.format(gabbase))
            # 2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.
            # ==> tiger.v.najna.XXXunique.tdt <==
            # Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd
            # In each case, Qry is Genome1 and Hit is Genome2
            quniq = '{0}.qryunique.tdt'.format(gabbase)
            huniq = '{0}.hitunique.tdt'.format(gabbase)
            ## ~ [2a] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.force() or not rje.exists(quniq) or not rje.exists(huniq):
                if self.getStrLC('Mapper') == 'busco':
                    self.buscoMap()
                else:
                    gabcmd = ['seqin={0}'.format(fasta1),'searchdb={0}'.format(fasta2),'mapper=minimap','minlocid=0','minloclen={0}'.format(self.getInt('MinLocLen')),'basefile={0}'.format(gabbase)]
                    gabobj = gablam.GABLAM(self.log,self.cmd_list+gabcmd)
                    gabobj.gablam()
            ## ~ [2b] Load tables and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                ufile = '{0}.{1}unique.tdt'.format(gabbase,qh.lower())
                udb = db.addTable(ufile,mainkeys=['Qry','Hit','AlnNum'],name=qh.lower(),ignore=[],expect=True)
                #udb = db.addTable(ufile,mainkeys=[qh,'{0}Start'.format(qh),'{0}End'.format(qh)],name=qh.lower(),ignore=[],expect=True)
                udb.dataFormat({'AlnNum':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                lenx = 0; idx = 0
                for entry in list(udb.entries()):
                    if entry['Length'] < self.getInt('MinLocLen'): udb.dropEntry(entry); lenx += 1
                    elif (100.0 * entry['Identity'] / entry['Length']) < self.getNum('MinLocID'): udb.dropEntry(entry); idx += 1
                self.printLog('#MINCUT','Dropped %s entries < %s bp and %s < %.1f%% identity' % (rje.iStr(lenx),rje.iStr(self.getInt('MinLocLen')),rje.iStr(idx),self.getNum('MinLocID')))
                udb.addField('HitStart',evalue=0)
                udb.addField('HitEnd',evalue=0)
                udb.addField('Strand',evalue='+')
                for entry in udb.entries():
                    if entry['SbjStart'] > entry['SbjEnd']:
                        entry['HitStart'] = entry['SbjEnd']
                        entry['HitEnd'] = entry['SbjStart']
                        entry['Strand'] = '-'
                    else:
                        entry['HitStart'] = entry['SbjStart']
                        entry['HitEnd'] = entry['SbjEnd']
                udb.newKey([qh,'{0}Start'.format(qh),'{0}End'.format(qh)])
                udb.setFields(['Qry','QryStart','QryEnd','Hit','HitStart','HitEnd','Strand','AlnNum','Length','Identity'])
                udb.addField('Non',evalue=0)
                udb.addField('QryGap',evalue='')
                udb.addField('HitGap',evalue='')
                udb.addField('Qry5',evalue='')
                udb.addField('Qry3',evalue='')
                udb.addField('Hit5',evalue='')
                udb.addField('Hit3',evalue='')

            return True
        except: self.errorLog('%s.runGABLAM error' % self.prog()); return False
#########################################################################################################################
    def buscoMap(self): ### Uses BUSCO mapping to generate mock gablam tables
        '''
        Loads BUSCO tables and uses to make synteny blocks, based on PAFScaff code. In each case, the Complete BUSCO
        genes from the "Complete" set can be mapped onto "Complete", "Fragmented" or "Duplicated" genes in the other
        assembly.

        Output tables are qryunique and hitunique local hits tables mimicking a GABLAM run using Minimap2:

        # Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd

        In each case, Qry is Genome1 and Hit is Genome2.
        :return:
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Compile BUSCO results')
            db = self.db()
            #i# v5 is the default:
            v5head = ['BuscoID','Status','Contig','Start','End','Strand','Score','Length','OrthoDBurl','Description']
            #!# Future updates could expand to non-BUSCO input and have custom headers.
            #!# Need ID, SeqName, Start, End, Strand
            qbusco = {'qry':'BUSCO1','hit':'BUSCO2'}
            hbusco = {'qry': 'BUSCO2', 'hit': 'BUSCO1'}
            basefile = self.baseFile()
            gabbase = basefile + '.map'
            if self.getStrLC('GABLAM'): gabbase = self.getStr('GABLAM')
            quniq = '{0}.qryunique.tdt'.format(gabbase)
            huniq = '{0}.hitunique.tdt'.format(gabbase)

            ### ~ [2] Load and filter data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Load tables: BUSCO1 and BUSCO2
            for btype in ['BUSCO1', 'BUSCO2']:
                tabhead = v5head
                bfile = self.getStr(btype)
                if not rje.exists(bfile):
                    raise IOError('{0} file "{1}" not found (pafin=busco)'.format(btype, bfile))
                fdb = db.addTable(bfile, mainkeys='auto', headers=tabhead, expect=True, name=btype)
                fdb.dropEntriesDirect('Status', ['Missing'], inverse=False)
                fdb.newKey(['Contig', 'Start', 'End'])
                fdb.keepFields(['BuscoID', 'Status', 'Contig', 'Start', 'End', 'Strand', 'Length'])
                fdb.dataFormat({'Start': 'int', 'End': 'int', 'Length': 'int'})
                fdb.index('BuscoID')
                if self.v() > 0 or self.debugging():
                    fdb.indexReport('Contig')
                else:
                    fdb.index('Contig')

            ### ~ [2] Load and filter data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                #i# Copy and filter to Complete in the focal assembly
                # i# Reduce to matching BuscoIDs
                qdb = db.copyTable(self.db(qbusco[qh]),'qbusco',add=False)
                hdb = db.copyTable(self.db(hbusco[qh]),'hbusco',add=False)
                #i# Reduced to Complete genes in Genome2 to avoid multiple mapping of query regions
                hdb.dropEntriesDirect('Status',['Complete'],inverse=True)
                qdb.dropEntriesDirect('BuscoID',hdb.orderedDataList('BuscoID'),inverse=True)
                hdb.dropEntriesDirect('BuscoID',qdb.orderedDataList('BuscoID'),inverse=True)
                ### ~ [3] Generate paf table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                #pfields = rje.split('# Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length BuscoID')
                pfields = rje.split('Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd Strand')
                tname = '{0}unique'.format(qh)
                pafdb = db.addEmptyTable(fields=pfields, keys=['AlnNum'], name=tname)

                ## ~ [3a] Synteny block method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Add counters for making synteny blocks
                i = 1
                for entry in qdb.entrySort():
                    entry['AlnNum'] = i
                    i += 1
                i = 1
                for entry in hdb.entrySort():
                    entry['AlnNum'] = i
                    i += 1
                #i# Generate synteny blocks
                pentry = {}; prevq = -1; prevh = -1
                for entry in qdb.entrySort():
                    #if entry['BuscoID'] not in fdb.index('BuscoID'): continue
                    for hentry in hdb.indexEntries('BuscoID',entry['BuscoID']):
                        #i# Set strand
                        strand = '+'
                        if entry['Strand'] != hentry['Strand']: strand = '-'
                        #i# New block, or continue?
                        newblock = abs(entry['AlnNum'] - prevq) != 1 or abs(hentry['AlnNum'] - prevh) != 1
                        prevq = entry['AlnNum']; prevh = hentry['AlnNum']
                        if not newblock and pentry and (entry['Contig'] != pentry['Qry'] or hentry['Contig'] != pentry['Hit'] or strand != pentry['Strand']):
                            newblock = True
                        if newblock:
                            if pentry: pafdb.addEntry(pentry)
                            #?# Add some kind of expectation of overlap by chance?
                            #?# Add Busco scores to the BitScore field?
                            pentry = {'AlnNum':pafdb.entryNum()+1,'Positives':1,'BitScore':0,'Expect':0,
                                      'Qry':entry['Contig'],'QryStart':entry['Start'],'QryEnd':entry['End'],'Strand':strand,
                                      'Hit':hentry['Contig'],'SbjStart':hentry['Start'],'SbjEnd':hentry['End'],
                                      'Identity':min(entry['Length'],hentry['Length']),
                                      'Length':max(entry['Length'],hentry['Length'])}
                        else:
                            pentry['QryEnd'] = entry['End']
                            if strand == '+': pentry['SbjEnd'] = hentry['End']
                            else: pentry['SbjStart'] = hentry['Start']
                            pentry['Identity'] += min(entry['Length'],hentry['Length'])
                            pentry['Length'] += max(entry['Length'],hentry['Length'])
                            pentry['Positives'] += 1
                        self.debug(pafdb.entrySummary(pentry,collapse=True))
                if pentry: pafdb.addEntry(pentry)
                self.printLog('#BUSCO','{0} BUSCO Complete synteny blocks added.'.format(rje.iStr(pafdb.entryNum())))
                if qh == 'hit':
                    for entry in pafdb.entries():
                        eswap = [entry['Qry'], entry['QryStart'], entry['QryEnd'], entry['Hit'], entry['SbjStart'], entry['SbjEnd']]
                        [entry['Hit'], entry['SbjStart'], entry['SbjEnd'], entry['Qry'], entry['QryStart'], entry['QryEnd']] = eswap
                for entry in pafdb.entries():
                    if entry['Strand'] == '-': [entry['SbjStart'], entry['SbjEnd']] = [entry['SbjEnd'], entry['SbjStart']]
                pafdb.dropField('Strand')
                if qh == 'qry':
                    pafdb.newKey(['Qry','QryStart','QryEnd'])
                    pafdb.saveToFile(quniq)
                else:
                    pafdb.newKey(['Hit','SbjStart','SbjEnd'])
                    pafdb.saveToFile(huniq)

            return True
        except: self.errorLog('%s.buscoMap error' % self.prog()); return False
#########################################################################################################################
    def runMashMap(self):   ### Runs Reciprocal MashMap mappings of two
        '''
        Runs Reciprocal MashMap mappings of two assemblies and loads results into database tables.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            gabbase = basefile + '.map'
            if self.getStrLC('GABLAM'): gabbase = self.getStr('GABLAM')
            self.printLog('#MASH','MashMap output basefile: {0}'.format(gabbase))
            ### ~ [1] Run MashMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Run MashMap, load results and convert:
            #i# udb.setFields(['Qry','QryStart','QryEnd','Hit','HitStart','HitEnd','Strand','AlnNum','Length','Identity'])
            # In each case, Qry is Genome1 and Hit is Genome2
            quniq = '{0}.qryunique.tdt'.format(gabbase)
            huniq = '{0}.hitunique.tdt'.format(gabbase)
            for qh in ('Qry','Hit'):
                ufile = '{0}.{1}unique.tdt'.format(gabbase,qh.lower())
                if self.force() or not rje.exists(ufile):
                    mashcmd = ['seqin={0}'.format(self.getStr('Genome1')), 'searchdb={0}'.format(self.getStr('Genome2'))]
                    if qh == 'Hit':
                        mashcmd = ['seqin={0}'.format(self.getStr('Genome2')), 'searchdb={0}'.format(self.getStr('Genome1'))]
                    mashcmd = mashcmd + ['basefile={0}.{1}'.format(gabbase,qh.lower()),'mashfilt=map']
                    mashobj = rje_mashmap.MashMap(self.log, self.cmd_list + mashcmd)
                    mashout = mashobj.run()
                    mashobj.obj['DB'] = db
                    mashobj.obj['DB'].info['Delimit'] = '\t'
                    udb = mashobj.parseOut(mashout,qh.lower())
                    #!# NOTE: These hits are not unique by query, despite map filter #!#
                    udb.info['Delimit'] = '\t'
                    if qh == 'Hit':
                        ex = 0.0; etot = udb.entryNum()
                        for entry in udb.entries():
                            self.progLog('Reformatting mashmap output: {0:.2f}%'.format(ex/etot)); ex += 100.0
                            (entry['Qry'],entry['QryStart'],entry['QryEnd'],entry['Hit'],entry['HitStart'],entry['HitEnd']) = (entry['Hit'],entry['HitStart'],entry['HitEnd'],entry['Qry'],entry['QryStart'],entry['QryEnd'])
                    udb.saveToFile(ufile)
                else:
                    mhead = ['Qry','QryLen','QryStart','QryEnd','Strand','Hit','HitLen','HitStart','HitEnd','PercID']
                    mkeys = ['Qry','QryStart','QryEnd','Strand','Hit','HitStart','HitEnd']
                    udb = db.addTable(ufile,mainkeys=mkeys,datakeys='All',delimit=' ',ignore=['#'],name=qh.lower(),expect=True,replace=True)
                    udb.dataFormat({'QryLen':'int','QryStart':'int','QryEnd':'int','HitLen':'int','HitStart':'int','HitEnd':'int','PercID':'num','AlnNum':'int','Length':'int','Identity':'int'})
                #i# Reformat for SynBad
                udb.newKey([qh,'{0}Start'.format(qh),'{0}End'.format(qh)])
                udb.setFields(['Qry','QryStart','QryEnd','Hit','HitStart','HitEnd','Strand','AlnNum','Length','Identity'])
                ## Add extra fields for SynBad
                udb.addField('Non',evalue=0)
                udb.addField('QryGap',evalue='')
                udb.addField('HitGap',evalue='')
                udb.addField('Qry5',evalue='')
                udb.addField('Qry3',evalue='')
                udb.addField('Hit5',evalue='')
                udb.addField('Hit3',evalue='')
            rje.checkForFiles([quniq,huniq])
            return True
        except: self.errorLog('%s.runMashMap error' % self.prog()); return False
#########################################################################################################################
    ### <7> ### Main SynBad gap mapping                                                                                 #
#########################################################################################################################
    ### ~ [3] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    def synBadMapping(self,cdb1,cdb2):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            self.printLog('#NOTE','Genome1 ({0}) is always Qry and Genome2 ({1}) is always Hit'.format(base1,base2))
            db = self.db()
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            wanted = ['qry.tdt','hit.tdt','qry.gap.tdt','hit.gap.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt, hit.tdt and *.gap.tdt files found: skipping SynBad mapping (force=F)')
                return True
            if not self.force() and not self.dev():
                self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                self.setBool({'Force':True})

            ### ~ [1] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
            # sorted and `QryGap` and `HitGap` fields added. Any local alignments with flanking hits are then flagged in these
            # new fields with the number of flanking 5', 3' gaps indicated by < and > characters.
            #
            # The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
            # corresponding local hits on the Qry and Hit genomes. If there is also an inversion, `SynSpan` will be negative.
            # If the local hits are against two different sequences from the other genome, the two sequence names will be
            # entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
            # small gaps), `SynSpan` or `GapSpan` will have a value of zero.
            ## ~ [1a] Work through query and hit tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug(self.db().tableNames())
            for qh in ('Qry','Hit'):
                qry = qh
                hit = {'Qry':'Hit','Hit':'Qry'}[qry]
                hitchr = {'Qry':chr1,'Hit':chr2}[hit]
                # gapdb = self.db('{0}gap'.format(qh.lower()))
                gapdb = self.dbTable(qh,'gap')
                spanfield = 'Span{0:d}'.format(self.getInt('SpannedFlank'))
                spancheck = spanfield in gapdb.fields()
                if self.getInt('SynReadSpan') > 0 and not spancheck:
                    self.warnLog('synreadspan={0} but "{1}" not found in fields. Check spannedflank=INT.'.format(self.getInt('SynReadSpan'),spanfield))
                spancheck = spancheck and self.getInt('SynReadSpan') > 0
                if spancheck: gapdb.dataFormat({spanfield:'int'})
                #self.debug(gapdb)
                altdb = {cdb1:cdb2,cdb2:cdb1}[gapdb]
                # locdb = self.db(qh.lower())
                locdb = self.dbTable(qh)
                locdb.index(qry)
                locdb.index(hit)
                for field in ['Hit3','Hit5','Qry3','Qry5']:
                    if field in locdb.fields(): locdb.dropField(field)
                    locdb.addField(field,after='HitGap',evalue='')
                # First, deal with QryGaps
                for seqname in gapdb.indexKeys('SeqName'):
                    if seqname not in locdb.index(qry):
                        self.printLog('#ORPHAN','No alignments found for {0}'.format(seqname))
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = gapdb.index('SeqName')[seqname][0:]
                    loc = [(seqname,0,0)] + locdb.index(qry)[seqname] + [(seqname,-1,-1)]
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        entry5 = None
                        entry3 = None
                        if loc[0][1] > 0:
                            entry5 = locdb.data(loc[0])
                            entry5['{0}Gap'.format(qry)] += '>'
                            if not entry5['{0}3'.format(qry)]: entry5['{0}3'.format(qry)] = []
                            entry5['{0}3'.format(qry)].append(gap[0])
                        if loc[1][1] > 0:
                            entry3 = locdb.data(loc[1])
                            entry3['{0}Gap'.format(qry)] += '<'
                            if not entry3['{0}5'.format(qry)]: entry3['{0}5'.format(qry)] = []
                            entry3['{0}5'.format(qry)].append(gap[0])
                        thisgap = gap.pop(0)
                        gentry = gapdb.data(thisgap)
                        if entry3 and thisgap[2] >= loc[1][1]: # Gap overlapping local alignment
                            gentry['GapSpan'] = 0
                            gentry['SynSpan'] = '{0}:0'.format(entry3[hit])
                            gentry['SynBad'] = 'Aligned'
                            continue
                        # At this point, loc[0][2] < thisgap[1] & thisgap[2] < loc[1][1]
                        if loc[1][1] == -1: gentry['GapSpan'] = gentry['SeqLen'] - loc[0][2]
                        else: gentry['GapSpan'] = loc[1][1] - loc[0][2]
                        if not entry5 or not entry3:
                            gentry['SynBad'] = 'Terminal'
                            if spancheck and gentry[spanfield] >= self.getInt('SynReadSpan'):
                                gentry['SynBad'] = 'Spanned'
                            continue
                        # Assess synteny and classify
                        hitstart = '{0}Start'.format(hit)
                        hitend = '{0}End'.format(hit)
                        gap5 = {'+':entry5[hitend],'-':entry5[hitstart]}[entry5['Strand']]
                        gap3 = {'-':entry3[hitend],'+':entry3[hitstart]}[entry3['Strand']]
                        synspan = max(gap5,gap3) - min(gap5,gap3)
                        #i# Look for picking up synteny by skipping fragments from other sequences
                        #i# Where a local alignment switches the hit sequence, SynBad will look downstream for another
                        #i# aligment to the same sequence. A maximum of `maxsynskip=INT` alignments will be skipped, up to
                        #i# a maximum distance of `maxsynspan=INT` bp.
                        maxskip = self.getInt('MaxSynSkip')     # Maximum number of local alignments to skip for SynTrans classification
                        entryskip = None
                        egap = -1
                        if entry3[hit] != entry5[hit]:
                            for i in range(maxskip):
                                if len(loc) > (2+i) and loc[2+i][1] > 0:
                                    entryskip = locdb.data(loc[2+i])
                                    if not entryskip:
                                        self.warnLog('Problem finding local hit entry for {0}'.format(str(loc[i+2])))
                                        continue
                                    goodskip = True
                                    for field in [hit,'Strand']:
                                        if entryskip[field] != entry5[field]: goodskip = False
                                    if entry5['Strand'] == '+' and entryskip[hitstart] < entry5[hitend]: goodskip = False
                                    elif entry5['Strand'] == '-' and entryskip[hitend] > entry5[hitstart]: goodskip = False
                                    elif entry5['Strand'] == '+': egap = entryskip[hitstart] - entry5[hitend]
                                    elif entry5['Strand'] == '-': egap = entryskip[hitend] - entry5[hitstart]
                                    if egap > self.getInt('MaxSynSpan'): goodskip = False
                                    if loc[2+i][1] - loc[0][2] > self.getInt('MaxSynSpan'): goodskip = False
                                    if goodskip: break
                                    entryskip = None
                        # Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:
                        #
                        # * `Aligned` = Gap is found in the middle of a local alignment to the Hit
                        # * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
                        # * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
                        # * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
                        # * `Duplication` = Overlapping flanking hits on the same strand.
                        # * `Inversion` = Flanking hits are on alternative strands.
                        # * `Translocation` = `SynSpan` indicates matches are on different scaffolds.
                        # * `Terminal` = Gap is between a local alignment and the end of the query sequence.
                        # * `Null` = No mapping between genomes for that gap.
                        if entryskip:
                            gentry['SynSpan'] = '{0}:{1}:{2}'.format(entry5[hit],egap,entryskip[hit])
                            gentry['SynBad'] = 'Insertion'
                        elif entry3[hit] != entry5[hit]:
                            gentry['SynSpan'] = '{0}::{1}'.format(entry5[hit],entry3[hit])
                            if hitchr and (not entry5[hit].startswith(hitchr) or not entry3[hit].startswith(hitchr)):
                                gentry['SynBad'] = 'Fragmentation'
                            else:
                                gentry['SynBad'] = 'Translocation'
                        elif entry3['Strand'] != entry5['Strand']:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Inversion'
                        elif (entry3['Strand'] == '+' and gap3 < gap5) or (entry3['Strand'] == '-' and gap3 > gap5):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],-synspan)
                            gentry['SynBad'] = 'Duplication'
                        elif (synspan - gentry['GapSpan']) > self.getInt('MaxSynSpan'):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Breakpoint'
                        else:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Syntenic'
                        if spancheck and gentry[spanfield] >= self.getInt('SynReadSpan') and gentry['SynBad'] not in ['Aligned','Syntenic']:
                            gentry['SynBad'] = 'Spanned'

                # Shorten classes
                devswap = {'Aligned':'Aln','Syntenic':'Syn','Insertion':'Ins','Breakpoint':'Brk','Duplication':'Dup',
                           'Inversion':'Inv','Translocation':'Tran','Fragmentation':'Frag','Terminal':'Term','Spanned':'Span','Null':'Null'}
                #!# Make this an option? if self.dev():
                for entry in gapdb.entries():
                    if entry['SynBad'] in devswap: entry['SynBad'] = devswap[entry['SynBad']]

                # Then, deal with Hit Gaps
                for seqname in altdb.indexKeys('SeqName'):
                    if seqname not in locdb.index(hit):
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = altdb.index('SeqName')[seqname][0:]
                    loc = [(seqname,0,0)]
                    loc2entry = {}
                    for entry in locdb.indexEntries(hit,seqname):
                        loc.append((entry[hit],entry['{0}Start'.format(hit)],entry['{0}End'.format(hit)]))
                        loc2entry[(entry[hit], entry['{0}Start'.format(hit)], entry['{0}End'.format(hit)])] = entry
                    loc.sort()
                    loc.append((seqname,-1,-1))
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        thisgap = gap.pop(0)
                        if loc[0][1] > 0:
                            entry = loc2entry[loc[0]]
                            entry['{0}Gap'.format(hit)] += '>'
                            if not entry['{0}3'.format(hit)]: entry['{0}3'.format(hit)] = []
                            entry['{0}3'.format(hit)].append(thisgap)
                        if loc[1][1] > 0:
                            entry = loc2entry[loc[1]]
                            entry['{0}Gap'.format(hit)] += '<'
                            if not entry['{0}5'.format(hit)]: entry['{0}5'.format(hit)] = []
                            entry['{0}5'.format(hit)].append(thisgap)
                self.progLog('\r#SYNBAD','{0} SynBad mapping complete.               '.format(qh))
                self.printLog('\r#SYNBAD','{0} SynBad mapping complete.'.format(qh))
                gapdb.indexReport('SynBad')
                #self.deBug('Continue?')

            ### ~ [2] Update tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qrygapdb = self.dbTable('qry','gap')
            hitgapdb = self.dbTable('hit','gap')
            qrydb = self.db('qry')
            hitdb = self.db('hit')
            ex = 0.0; etot = qrydb.entryNum() + hitdb.entryNum()
            for entry in list(qrydb.entries()) + list(hitdb.entries()):
                self.progLog('\r#UPDATE','Updating flanking gap SynBad ratings: %.2f%%' % (ex/etot)); ex += 100
                for field in ('Qry5','Qry3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(qrygapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                for field in ('Hit5','Hit3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(hitgapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                if not entry['QryGap']: entry['QryGap'] = '.'
                if not entry['HitGap']: entry['HitGap'] = '.'
            self.printLog('\r#UPDATE','Updating flanking gap SynBad ratings complete.')

            return True
        except: self.errorLog('%s.synBadMapping error' % self.prog()); return False
#########################################################################################################################
    ### <8> ### Assembly map methods                                                                                    #
#########################################################################################################################
    def dataCheck(self,qh='both'):    ### Performs a full check/summary of data tables and assembly maps
        '''
        Performs a full check/summary of data tables and assembly maps.
        :return: bool [True/False]
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if qh == 'both': return self.dataCheck('qry') and self.dataCheck('hit')
            db = self.db()
            #self.deBug(', '.join(db.tableNames()))
            amap = self.list[qh]
            fdb = self.dbTable(qh,'flanks',expect=False)
            self.headLog('Checking {0} flank data integrity'.format(qh),line='~')
            if not amap and not fdb:
                self.printLog('#CHECK','No assembly map or flanks table: no flank check possible!')
                return True
            if amap: self.printLog('#CHECK','Checking against assembly map.')
            else: self.printLog('#CHECK','No assembly map for checking.')
            aflanks = []
            if fdb:
                self.printLog('#CHECK','Checking against flanks table.')
                aflanks = fdb.dataList(fdb.entries(),'Name',sortunique=True,empties=False)
            else: self.printLog('#CHECK','No flanks table for checking.')
            ### ~ [1] Check flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            problems = 0
            allbad = []
            fcheck = {'blocks':['Flank5','Flank3'],'contigs':['Flank5','Flank3'],'corrections':['Flank1','Flank2'],
                      'flanks':['Name','Pair'],'gap':['GapFlank5','GapFlank3','BestFlank5','BestFlank3'],
                      'hicbest':['Flank1','Flank2'],'hicpairs':['Flank1','Flank2']}
            for ttype in rje.sortKeys(fcheck):
                table = self.dbTable(qh,ttype,expect=False)
                if not table:
                    self.vLog('#CHECK','{0}.{1}: False'.format(qh,ttype))
                    continue
                self.vLog('#CHECK','{0}.{1}: {2} entries'.format(qh,ttype,table.entryNum()))
                for field in fcheck[ttype]:
                    if field not in table.fields(): continue
                    flanks = table.dataList(table.entries(),field,sortunique=True,empties=False)
                    if '?' in flanks: flanks.remove('?')
                    if amap:
                        badflanks = rje.listDifference(flanks,amap)
                        msg = '{0} of {1} {2} {3} flanks missing from {4} assembly map.'.format(rje.iLen(badflanks),rje.iLen(flanks),ttype,field,qh)
                        if badflanks:
                            self.warnLog(msg)
                            problems += 1
                            allbad += badflanks
                        else: self.devLog('#CHECK',msg,debug=False)
                    if aflanks:
                        badflanks = rje.listDifference(flanks,aflanks)
                        msg = '{0} of {1} {2} {3} flanks missing from {4} flanks.'.format(rje.iLen(badflanks),rje.iLen(flanks),ttype,field,qh)
                        if badflanks:
                            self.warnLog(msg)
                            problems += 1
                            allbad += badflanks
                        else: self.devLog('#CHECK',msg,debug=False)
            if problems:
                self.printLog('#CHECK','{0} {1} fields have missing flanks'.format(problems,qh))
                allbad = rje.sortUnique(allbad)
                if len(allbad) > 8:
                    raise ValueError('{0} flanks missing: {1} ... {2}'.format(len(allbad),','.join(allbad[:4]),','.join(allbad[-4:])))
                else:
                    raise ValueError('{0} flanks missing: {1}'.format(len(allbad),','.join(allbad)))
                #if self.i() < 0 or rje.yesNo('Quit SynBad?'): return False
            else: self.printLog('#CHECK','No {0} flank data integrity issues detected.'.format(qh))
            return True
        except: self.errorLog('%s.dataCheck error' % self.prog()); return False
#########################################################################################################################
    def contigMapping(self):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            self.seqObjSetup()
            ## ~ [0a] Check and load outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            wanted = []     #X#'qry.tdt','hit.tdt']
            for qh in ['qry','hit']:
                wanted.append('{0}.flanks.bed'.format(qh))
                wanted.append('{0}.flanks.fasta'.format(qh))
                for tname in ['flanks','contigs','gap']:
                    wanted.append('{0}.{1}.tdt'.format(qh,tname))
            makecontigs = False
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad Flank and Contig files found: loading (force=F)')
                problems = False
                for qh in ['qry','hit']:
                    gdb = db.addTable(name='{0}.gap'.format(qh),mainkeys=['SeqName','Start','End'],expect=True,ignore=[],replace=True)
                    gdb.dataFormat(dbformats)
                    for field in ['GapFlank5','GapFlank3']:
                        if field not in gdb.fields(): problems = True
                    bed = db.addTable(name='{0}.flanks'.format(qh),mainkeys=['SeqName','Start','End'],expect=True,ignore=[],replace=True)
                    bed.dataFormat(dbformats)
                    cdb = db.addTable(name='{0}.contigs'.format(qh),mainkeys=['SeqName','Start','End'],expect=True,ignore=[],replace=True)
                    cdb.dataFormat(dbformats)
                if problems:
                    self.warnLog('Problems found with loaded data: regenerating contigs and flanks')
                    makecontigs = True
            else: makecontigs = True
            ## ~ [0a] Check and load gap tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ['qry','hit']:
                if not self.dbTable(qh,'gap',expect=False):
                    gdb = db.addTable(name='{0}.gap'.format(qh),mainkeys=['SeqName','Start','End'],expect=True,ignore=[],replace=True)
                    gdb.dataFormat(dbformats)

            ### ~ [1] Extract gap flanks to BED and fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if makecontigs:
                if not self.force() and not self.dev():
                    self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                    self.setBool({'Force':True})
                self.headLog('Generating Contigs and Flanks',line='-')
                for qh in ('qry','hit'):
                    self.makeContigsAndFlanks(qh)
                ## ~ [1a]  Update flank pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('Updating Flanks table with contig pairs',line='~')
                if not self.flankPairs(): raise ValueError('Failed to generate flank pairs from contigs')
                for qh in ('qry','hit'):
                    flankpairs = self.dict['FlankPairs'][qh]
                    bed = self.dbTable(qh,'flanks')
                    for entry in bed.entries():
                        entry['Pair'] = flankpairs[entry['Name']]
                    bed.saveToFile()

            ### ~ [2] Generate assembly maps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Maps',line='-')
            for qh in ['qry','hit']:
                wanted = ['map.fasta','map.txt']
                qhbase = '{0}.{1}.'.format(self.baseFile(),qh)
                ### ~ [1] Check/Load Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if not self.force() and rje.checkForFiles(wanted,basename=qhbase,log=self.log,cutshort=False,ioerror=False,missingtext='Not found.'):
                    self.printLog('#MAP','{0} assembly map found: loading (force=F)'.format(qh))
                    if self.loadAssemblyMap(qh): continue
                    maptxt = qhbase + 'map.txt'
                    self.list[qh] = []
                    try:
                        for line in open(maptxt,'r').readlines():
                            #self.bugPrint(line)
                            data = rje.split(line)
                            if not data: continue
                            if data[1] == '=': self.list[qh] += data[2:]
                        continue
                    except:
                        self.errorLog('Problem processing {0}map.txt: will regenerate'.format(qhbase))
                        self.list[qh] = []
                if not self.makeAssemplyMap(qh): return False
                if not self.saveAssemblyMaps(qh,mapname='map'): return False
                if not self.saveTelociraptorMaps(qh,mapname='map'): return False

            ### ~ [3] Extra flanks analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Flank copy number ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog('Flank copy number',line='-')
            for qh in ['qry','hit']:
                cdb = self.flankCNV(qh)
                if not cdb: continue
                fdb = self.dbTable(qh,'flanks')
                for field in ['MeanX','CN']:
                    if field not in fdb.fields(): fdb.addField(field,evalue=0)
                gdb = self.dbTable(qh,'gap')
                for field in ['CN5','CN3']:
                    if field not in gdb.fields(): gdb.addField(field,evalue=0)
                cx = 0
                for ekey in list(fdb.dataKeys()):
                    if cdb.data(ekey):
                        entry = fdb.data(ekey)
                        centry = cdb.data(ekey)
                        for field in ['MeanX','CN']:
                            entry[field] = centry[field]
                        cx += 1
                        for gentry in gdb.indexEntries('GapFlank5',entry['Name']):
                            gentry['CN5'] = centry['CN']
                        for gentry in gdb.indexEntries('GapFlank3',entry['Name']):
                            gentry['CN3'] = centry['CN']
                self.printLog('#CNV','Updated {0} flanks table with {1} Mean X depth and CN estimates'.format(qh,cx))

            ## ~ [3b] KAT kmer analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog('KAT kmer analysis',line='-')
            kat = os.popen('kat --version 2>&1').read()
            if kat and "not found" not in kat:
                self.printLog('#KAT','KAT version detected: {0}'.format(rje.chomp(kat)))
                for qh in ['qry','hit']:
                    kdb = self.kmerFrequencies(qh)
                    kdb.dataFormat({'median':'int','mean':'num'})
                    if not kdb: continue
                    fdb = self.dbTable(qh,'flanks')
                    for field in ['katMed','katMean']:
                        if field not in fdb.fields(): fdb.addField(field,evalue=0)
                    gdb = self.dbTable(qh,'gap')
                    for field in ['kat5','kat3']:
                        if field not in gdb.fields(): gdb.addField(field,evalue=0)
                    kx = 0
                    for entry in fdb.entries():
                        if kdb.data(entry['Name']):
                            kentry = kdb.data(entry['Name'])
                            entry['katMed'] = kentry['median']
                            entry['katMean'] = kentry['mean']
                            kx += 1
                            for gentry in gdb.indexEntries('GapFlank5',entry['Name']):
                                gentry['kat5'] = kentry['median']
                            for gentry in gdb.indexEntries('GapFlank3',entry['Name']):
                                gentry['kat3'] = kentry['median']
                    self.printLog('#KAT','Updated {0} flanks table with {1} kat kmer Median and Mean values'.format(qh,kx))
            else:
                self.warnLog('Cannot run "kat --version": check installation or pre-generation of files')

            return True
        except:
            self.errorLog('%s.contigMapping error' % self.prog())
            return False
#########################################################################################################################
    def loadAssemblyMap(self,qh,mapname='map',check=True): ### Loads map and populates self.list['qry] and self.list['hit']
        '''
        Loads map and populates self.list['qry] and self.list['hit'].
        :param qh:str = qry/hit
        :param mapname: str [map] = map name element of file
        :param check:bool [True] = whether to check integrity against gap, contigs and flanks tables.
        :return:
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Loading {0} "{1}" assembly map'.format(qh,mapname),line='~')
            qh = qh.lower()
            maptxt = '{0}.{1}.{2}.txt'.format(self.baseFile(),qh,mapname)
            self.list[qh] = []

            for line in open(maptxt,'r').readlines():
                #self.bugPrint(line)
                data = rje.split(line)
                if not data: continue
                if '=' in data: self.list[qh] += ['|'] + data[data.index('=')+1:] + ['|']

            seqterm = self.list[qh].count('|')
            self.printLog('#MAP','Loaded assembly map from {0}: {1} sequences'.format(maptxt,seqterm/2))

            if check:
                return self.dataCheck(qh)

            return True
        except:
            self.errorLog('%s.loadAssemblyMap error' % self.prog())
            self.list[qh] = []
            return False
#########################################################################################################################
    def makeContigsAndFlanks(self,qh):  ### Generates contig and flanks data. Updates gaps table.
        '''
        Generates contig and flanks data. Updates gaps table.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            flanklen = self.getInt('GapFlanks')
            db = self.db()
            gdb = self.dbTable(qh,'gap')
            for field in ['GapFlank5','GapFlank3']:
                if field not in gdb.fields(): gdb.addField(field,evalue='')
            bedfile = '{0}.{1}.flanks.bed'.format(db.baseFile(),qh)
            bed = db.addEmptyTable('{0}.flanks'.format(qh),['SeqName','Start','End','Name','Score','Strand','Pair'],['SeqName','Start','End'],log=True)
            cdb = db.addEmptyTable('{0}.contigs'.format(qh),['SeqName','Start','End','CtgLen','Name','Flank5','Flank3','SynBad'],['SeqName','Start','End'],log=True)
            ### ~ [1] Extract Flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            entry = None
            entries = gdb.entries(sorted=True)
            while entries:
                prev = entry
                entry = entries.pop(0)
                #self.bugPrint('\nGap :: %s' % (gdb.entrySummary(entry,collapse=True)))
                # 5' flank
                i = entry['Start'] - flanklen
                if prev and prev['SeqName'] == entry['SeqName'] and self.getBool('PureFlanks'): i = max(prev['End'] + 1,i)
                else: i = max(1,i)
                j = entry['Start'] - 1
                score = 0
                #if 'Span0' in gdb.fields(): score = entry['Span0']
                if i <= j:
                    bentry = bed.addEntry({'SeqName':entry['SeqName'],'Start':i,'End':j,'Name':'{0}.{1}-{2}'.format(entry['SeqName'],i,j),'Score':score,'Strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    entry['GapFlank5'] = '{0}.{1}-{2}'.format(entry['SeqName'],i,j)
                # 3' flank
                i = entry['End'] + 1
                j = min(entry['SeqLen'],entry['End'] + flanklen)
                if entries and entries[0]['SeqName'] == entry['SeqName']:
                    j = min(entries[0]['Start']-1,j)
                if i <= j:
                    bentry = bed.addEntry({'SeqName':entry['SeqName'],'Start':i,'End':j,'Name':'{0}.{1}-{2}'.format(entry['SeqName'],i,j),'Score':score,'Strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    entry['GapFlank3'] = '{0}.{1}-{2}'.format(entry['SeqName'],i,j)
                #i# Extra end of sequence contig
                if prev and entry['SeqName'] != prev['SeqName']:
                    i = max(prev['End']+1,prev['SeqLen']-flanklen+1)
                    j = prev['SeqLen']
                    if i <= j:
                        bed.addEntry({'SeqName':prev['SeqName'],'Start':i,'End':j,
                                      'Name':'{0}.{1}-{2}'.format(prev['SeqName'],i,j),'Score':0,'Strand':'+'},warn=False)
                        centry = {'SeqName':prev['SeqName'],'Start':prev['End']+1,'End':prev['SeqLen'],
                                  'Flank5':prev['GapFlank3'],'Flank3':'{0}.{1}-{2}'.format(prev['SeqName'],i,prev['SeqLen']),
                                  'SynBad':'{0}-End'.format(prev['SynBad'])}
                        centry['Name'] = '{0}.{1}-{2}'.format(centry['SeqName'],centry['Start'],centry['End'])
                        centry['CtgLen'] = centry['End'] - centry['Start'] + 1
                        cdb.addEntry(centry)
                        #self.bugPrint('-> Extra end contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                #i# Start of sequence contig
                if not prev or prev['SeqName'] != entry['SeqName']:
                    i = 1
                    j = min(entry['Start'] - 1,flanklen)
                    if i <= j:
                        bentry = bed.addEntry({'SeqName':entry['SeqName'],'Start':i,'End':j,
                                      'Name':'{0}.{1}-{2}'.format(entry['SeqName'],i,j),'Score':0,'Strand':'+'},warn=False)
                        #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                        centry = {'SeqName':entry['SeqName'],'Start':1,'End':entry['Start'] - 1,
                                  'Flank5':'{0}.{1}-{2}'.format(entry['SeqName'],i,j),'Flank3':entry['GapFlank5'],
                                  'SynBad':'End-{0}'.format(entry['SynBad'])}
                        centry['Name'] = '{0}.{1}-{2}'.format(centry['SeqName'],centry['Start'],centry['End'])
                        centry['CtgLen'] = centry['End'] - centry['Start'] + 1
                        cdb.addEntry(centry)
                        #self.bugPrint('-> Start contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                #i# Middle of sequence contig
                else:
                    centry = {'SeqName':entry['SeqName'],'Start':prev['End']+1,'End':entry['Start'] - 1,
                              'Flank5':prev['GapFlank3'],'Flank3':entry['GapFlank5'],
                              'SynBad':'{0}-{1}'.format(prev['SynBad'],entry['SynBad'])}
                    centry['Name'] = '{0}.{1}-{2}'.format(centry['SeqName'],centry['Start'],centry['End'])
                    centry['CtgLen'] = centry['End'] - centry['Start'] + 1
                    cdb.addEntry(centry)
                    #self.bugPrint('-> Mid contig: %s' % (cdb.entrySummary(centry,collapse=True)))
            prev = entry
            if prev:
                i = max(prev['End']+1,prev['SeqLen']-flanklen+1)
                j = prev['SeqLen']
                if i <= j:
                    bentry = bed.addEntry({'SeqName':prev['SeqName'],'Start':i,'End':j,
                              'Name':'{0}.{1}-{2}'.format(prev['SeqName'],i,j),'Score':0,'Strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    centry = {'SeqName':prev['SeqName'],'Start':prev['End']+1,'End':prev['SeqLen'],
                              'Flank5':prev['GapFlank3'],'Flank3':'{0}.{1}-{2}'.format(prev['SeqName'],i,prev['SeqLen']),
                              'SynBad':'{0}-End'.format(prev['SynBad'])}
                    centry['Name'] = '{0}.{1}-{2}'.format(centry['SeqName'],centry['Start'],centry['End'])
                    centry['CtgLen'] = centry['End'] - centry['Start'] + 1
                    cdb.addEntry(centry)
                    #self.bugPrint('-> Final end contig: %s' % (cdb.entrySummary(centry,collapse=True)))

            #i# Add sequences without gaps to flanks and contigs
            seqdict = seqobj[qh].seqNameDic()
            sx = 0
            for seq in seqobj[qh].seqs():
                seqname = seqobj[qh].shortName(seq)
                if seqname not in gdb.index('SeqName'):
                    seqlen = seqobj[qh].seqLen(seq)
                    bentry = {'SeqName':seqname,'Start':1,'End':min(flanklen,seqlen),'Score':0,'Strand':'+'}
                    bentry['Name'] = '{0}.{1}-{2}'.format(bentry['SeqName'],bentry['Start'],bentry['End'])
                    bed.addEntry(bentry,warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    centry = {'SeqName':seqname,'Start':1,'End':seqlen,'Flank5':bentry['Name'],'SynBad':'End-End'}
                    bentry = {'SeqName':seqname,'Start':max(1,seqlen-flanklen+1),'End':seqlen,'Score':0,'Strand':'+'}
                    bentry['Name'] = '{0}.{1}-{2}'.format(bentry['SeqName'],bentry['Start'],bentry['End'])
                    bed.addEntry(bentry,warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    centry['Flank3'] = bentry['Name']
                    centry['Name'] = '{0}.{1}-{2}'.format(centry['SeqName'],centry['Start'],centry['End'])
                    centry['CtgLen'] = centry['End'] - centry['Start'] + 1
                    cdb.addEntry(centry); sx += 1
                    #self.bugPrint('-> Full-length contig: %s' % (cdb.entrySummary(centry,collapse=True)))
            self.printLog('#CONTIG','Added flanks for {0} sequences without gaps'.format(rje.iStr(sx)))
            cdb.saveToFile()

            #i# Save contigs to fasta file
            fasout = '{0}.{1}.contigs.fasta'.format(db.baseFile(),qh)
            if rje.exists(fasout) and not self.force():
                self.printLog('#CTG','{0} found (force=F)'.format(fasout))
            else:
                FASOUT = open(fasout,'w'); sx = 0
                seqlist = seqobj[qh]
                cx = 0.0; ctot = cdb.entryNum()
                for seq in seqlist.seqs():
                    seqname = seqlist.shortName(seq)
                    if seqname in cdb.index('SeqName'):
                        sequence = seqlist.getSeq(seq)[1]
                        for centry in cdb.indexEntries('SeqName',seqname):
                            self.progLog('\r#CTG','Extracting contig sequences: {0:.1f}%'.format(cx/ctot)); cx += 100.0
                            flankseq = sequence[centry['Start']-1:centry['End']]
                            if len(flankseq) != centry['End'] - centry['Start'] + 1:
                                self.debug(cdb.entrySummary(centry,collapse=True))
                                raise ValueError('Contig sequence length mismatch for {0}'.format(centry['Name']))
                            FASOUT.write('>{0}\n{1}\n'.format(centry['Name'],flankseq)); sx += 1
                    else:
                        self.warnLog('Sequence {0} missing from contigs table'.format(seqname))
                FASOUT.close()
                self.printLog('\r#CTG','Extracted {0} sequences for {1} contigs'.format(rje.iStr(sx),rje.iStr(ctot)))
                if sx != ctot: raise ValueError('Failure to generate full contigs fasta file!')

            #i# Save flanks to fasta file
            fasout = '{0}.{1}.flanks.fasta'.format(db.baseFile(),qh)
            if rje.exists(fasout) and not self.force():
                self.printLog('#FLANK','{0} found (force=F)'.format(fasout))
            else:
                FASOUT = open(fasout,'w'); sx = 0
                for seqname in bed.index('SeqName'):
                    sequence = seqobj[qh].getSeq(seqdict[seqname])[1]
                    #self.bugPrint('{0}: {1} bp'.format(seqname,len(sequence)))
                    for entry in bed.indexEntries('SeqName',seqname):
                        self.progLog('\r#FLANK','{0} {1} gap flanks output to {2}'.format(rje.iStr(sx),qh,fasout)); sx += 1
                        flankseq = sequence[entry['Start']-1:entry['End']]
                        if len(flankseq) != entry['End'] - entry['Start'] + 1:
                            self.debug(bed.entrySummary(entry,collapse=True))
                            raise ValueError('Flanking sequence length mismatch for {0}'.format(entry['Name']))
                        FASOUT.write('>{0}\n{1}\n'.format(entry['Name'],flankseq))
                FASOUT.close()
                self.printLog('\r#FLANK','{0} {1} gap flanks output to {2}'.format(rje.iStr(bed.entryNum()),qh,fasout))
            bed.saveToFile(bedfile,delimit='\t',headers=False,savefields=['SeqName','Start','End','Name'])

            return True
        except:
            self.errorLog('%s.makeContigsAndFlanks error' % self.prog())
            return False
#########################################################################################################################
    def flankPairs(self):   ### Checks/creates FlankPairs dictionaries
        '''
        Checks/creates FlankPairs dictionaries.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if 'FlankPairs' not in self.dict or not self.dict['FlankPairs']: self.dict['FlankPairs'] = {'qry':{},'hit':{}}
            ### ~ [1] Map flank pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                flankpairs = self.dict['FlankPairs'][qh]     # Dictionary of contig flanks pairs for internal paired read removal
                if flankpairs: continue
                cdb = self.dbTable(qh,'contigs',expect=False)
                if not cdb:
                    cdb = db.addTable(name='{0}.contigs'.format(qh),mainkeys=['SeqName','Start','End'],expect=True)
                for centry in cdb.entries():
                    flankpairs[centry['Flank5']] = centry['Flank3']
                    flankpairs[centry['Flank3']] = centry['Flank5']
            return True
        except:
            self.errorLog('%s.flankPairs error' % self.prog())
            return False
#########################################################################################################################
    def mapNewFlanks(self):   ### Sets up old-to-new Contig End flank regions for HiC mapping.
        '''
        Sets up old-to-new Contig End flank regions for HiC mapping.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Set up Flank Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['FlankMap'] = {'qry':{}, 'hit':{}}
            if not self.getBool('PureFlanks'):
                self.setBool({'PureFlanks':True})
                self.warnLog('pureflanks=F not currently supported: setting pureflanks=T')
            ### ~ [1] Map old flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                mdb = self.mapFlanksToGenome(qh)
                if mdb:
                    for entry in mdb.entries():
                        self.dict['FlankMap'][qh][entry['GapFlank']] = entry['OldFlank']
                self.printLog('#FLANKS','{0} {1} flank mappings.'.format(rje.iLen(self.dict['FlankMap'][qh]),qh))
            return True
        except:
            self.errorLog('%s.mapNewFlanks error' % self.prog())
            return False
#########################################################################################################################
    def mapFlanksToGenome(self,qh):   ### Runs GABLAM of flanks against assembly to map old flanks onto new assembly.
        '''
        Runs GABLAM of flanks against assembly to map old flanks onto new assembly.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('{0} Flank GABLAM mapping'.format(qh),line='-')
            db = self.db()
            basefile = self.baseFile()
            gabbase = '{0}.{1}.flanks'.format(basefile,qh)
            wanted = ['.local.tdt','.map.tdt']
            maptable = '{0}.flanks.map'.format(qh)
            oldflanks = {'qry':'MapFlanks1','hit':'MapFlanks2'}[qh]
            if not self.getStrLC(oldflanks):
                self.printLog('#FLANKS','No {0}=FILE flank mapping provided. Will use current flank locations for HiC analysis.'.format(oldflanks.lower()))
                return None
            elif not rje.exists(self.getStr(oldflanks)):
                raise IOError('{0}={1} not found!'.format(oldflanks.lower(), self.getStr(oldflanks)))

            ### ~ [1] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.force() and rje.checkForFiles(wanted,basename=gabbase,log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Flanks local hits mapping found: skipping GABLAM (force=F)')
            #!# Need to rationalise use of force=T
            #if rje.checkForFiles(wanted,basename=gabbase,log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
            #    self.printLog('#SKIP','Flanks local hits mapping found: skipping GABLAM.')
                mdb = db.addTable(name=maptable,mainkeys=['GapFlank'],ignore=[])
                return mdb
            ## ~ [1a] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            loctable = '{0}.local.tdt'.format(gabbase)
            flankfas = '{0}.{1}.flanks.fasta'.format(basefile,qh)
            #!# Add more intelligent minloclen filter based on contig size?
            if self.force() or not rje.exists(loctable):
                gabcmd = ['seqin={0}'.format(flankfas),'searchdb={0}'.format(self.getStr(oldflanks)),'mapper=minimap','minlocid=99','minloclen=0','basefile={0}'.format(gabbase),'uniqueout=F']
                gabobj = gablam.GABLAM(self.log,self.cmd_list+gabcmd+['debug=F'])
                gabobj.gablam()
            ## ~ [2b] Load tables and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mdb = db.addEmptyTable(maptable,['GapFlank','OldFlank'],['GapFlank'])
            ldb = db.addTable(loctable,mainkeys=['Qry','Hit','AlnNum'],name=loctable,ignore=[],expect=True)
            ldb.dataFormat({'AlnNum':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            #i# Qry is the new flank. Hit is the old flank. Should not need to consider the Strand.
            for entry in ldb.sortedEntries('Length',reverse=True):
                if entry['Qry'] in mdb.data(): continue
                mdb.addEntry({'GapFlank':entry['Qry'],'OldFlank':entry['Hit']})
            mdb.saveToFile()
            return mdb
        except: self.errorLog('%s.mapFlanksToGenome error' % self.prog()); return None
#########################################################################################################################
    def flankCNV(self,qh):    ### Runs Diploidocus RegCNV analysis on the flanks and returns regcnv table.
        '''
        Runs Diploidocus RegCNV analysis on the flanks and returns regcnv table. Generates *.checkcnv.tdt
        # SeqName Start   End     Name    Score   Strand  Pair   SeqBP   ReadBP  MeanX   ModeX   CN      CIsyst  CIrand
        :return: (cdb1,cdb2)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qh = qh.lower()
            G = {'qry':1,'hit':2}[qh]   # Genome number
            db = self.db()
            basefile = '{0}.{1}.flanks'.format(self.baseFile(),qh)
            ## ~ [0a] Check for existing results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb = self.addTable(qh,'flanks.checkcnv',make=False)
            #!# Replace some of the force=T/F with date checks
            if cdb:
                self.printLog('#SKIP','{0} flank CNV file found: skipping Diploidocus run.'.format(qh))
                return cdb

            #!# Replace with direct DepthKopy run and maybe even change the table read in?

            ### ~ [1] Run Diploidocus on flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.bugPrint(self.int)
            gensize = self.getInt('GenomeSize{0}'.format(G))
            self.printLog('#GSIZE','{0} genome size: {1}'.format(qh,rje_seqlist.dnaLen(gensize)))
            self.bugPrint(self.num)
            scdepth = self.getNum('SCDepth{0}'.format(G))
            self.printLog('#SCDEP','{0} single copy depth: {1:.2f}X'.format(qh,scdepth))
            dcmd = ['bam=', 'paf=', 'reads=', 'readtype=ont','scdepth=0','genomesize=0']
            #?#dcmd += ['busco=','kmerself=F','homfile=F','seqstats=F']
            gbase = rje.baseFile(self.getStr('Genome{0}'.format(G)), strip_path=True)
            base = basefile
            cmd1 = ['seqin={0}'.format(self.getStr('Genome{0}'.format(G))),
                    'runmode=regcnv'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base),'regcheck={0}.tdt'.format(basefile),
                    'checkfields=SeqName,Start,End','checkflanks={0}'.format(','.join([str(i) for i in self.list['CheckFlanks']]))]
            if self.getStrLC('BAM{0}'.format(G)): cmd1.append('bam={0}'.format(self.getStr('BAM{0}'.format(G))))
            else: cmd1.append('bam={0}.bam'.format(gbase))
            if self.getStrLC('PAF{0}'.format(G)): cmd1.append('paf={0}'.format(self.getStr('PAF{0}'.format(G))))
            else: cmd1.append('paf={0}.paf'.format(gbase))
            if self.getStrLC('BUSCO{0}'.format(G)): cmd1.append('busco={0}'.format(self.getStr('BUSCO{0}'.format(G))))
            if self.getNum('SCDepth{0}'.format(G)): cmd1.append('scdepth={0}'.format(self.getNum('SCDepth{0}'.format(G))))
            if gensize: cmd1.append('genomesize={0}'.format(gensize))
            if self.list['Reads{0}'.format(G)]: cmd1.append('reads={0}'.format(','.join(self.list['Reads{0}'.format(G)])))
            if self.list['ReadType{0}'.format(G)]: cmd1.append('readtype={0}'.format(','.join(self.list['ReadType{0}'.format(G)])))
            self.debug(self.cmd_list + dcmd + cmd1)
            dip = diploidocus.Diploidocus(self.log, self.cmd_list + dcmd + cmd1)
            #i# Check inputs and abort if not provided
            if not rje.exists(dip.getStr('BAM')) and not dip.list['Reads']:
                self.printLog('#READS','No reads=FILELIST files or bam=FILE given: no flank depth analysis.')
                return None
            #i# Run Diploidocus (to run DepthKopy!)
            if dip.run():
                cnvfile = '{0}.regcnv.tsv'.format(basefile)
                newcnvfile = '{0}.{0}.regcnv.tsv'.format(basefile)
                if rje.exists(newcnvfile) and not rje.exists(cnvfile):
                    cnvfile =  newcnvfile
                if rje.exists(cnvfile):
                    os.rename(cnvfile,'{0}.checkcnv.tdt'.format(basefile))
                    self.printLog('#MOVE','{0} -> {1}.checkcnv.tdt'.format(cnvfile,basefile))
                cdb = self.addTable(qh, 'flanks.checkcnv', make=False)
                # !# Replace some of the force=T/F with date checks
                if cdb: return cdb
            self.warnLog('Failed to generate {0} file'.format(basefile))
            return None
        except: self.errorLog('%s.flankCNV error' % self.prog()); raise
#########################################################################################################################
    def kmerFrequencies(self,qh):   ### Generate KAT assembly kmer summary
        '''
        Generate KAT assembly kmer summary.
        :param qh: qry/hit
        :return: KAT kmer frequency table for flanks kmer stats.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qh = qh.lower()
            genome = {'qry':self.getStr('Genome1'),'hit':self.getStr('Genome2')}[qh]
            basefile = '{0}.{1}.flanks'.format(self.baseFile(),qh)
            flanks = '{0}.fasta'.format(basefile)
            katfile = '{0}.kat-stats.tsv'.format(basefile)
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            # NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1 28      111164.29792    0.41661 122971501       122971475       54368   0.04421 122619711       99.71395        99.75805
            katcvg =  '{0}.kat-counts.cvg'.format(basefile)
            # ${PREFIX}-counts.cvg
            # >NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1
            # 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 11 528 585 553 470 111 ...
            rje.checkForFiles(filelist=[katfile,katcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            katcall = 'kat sect -t {0} -o {1}.kat {2} {3}'.format(self.threads(),basefile,flanks,genome)
            #!# Need to re-rationalise the use of force. Have an extract 'recovery' mode to use instead of setting force=T
            if not (rje.exists(katfile) and rje.exists(katcvg)):    #!# or self.force()
                self.printLog('#SYS',katcall)
                #i# Catching completion in case KAT hangs after running
                KAT = os.popen(katcall)
                while not KAT.readline().startswith('Total runtime'): continue
                KAT.close()

            kdb = self.db().addTable(katfile,mainkeys=['seq_name'],name='{0}.kat-stats'.format(qh))
            kdb.dataFormat({'median':'int','mean':'num','gc%':'num','seq_length':'int','kmers_in_seq':'int','invalid_kmers':'int','%_invalid':'num','non_zero_kmers':'int','%_non_zero':'num','%_non_zero_corrected':'num'})
            #!# Add processing! Rename fields and add to flanks table.
            return kdb

        except: self.errorLog('%s.kmerFrequencies error' % self.prog()); raise
#########################################################################################################################
    #!# Will want to move these functions out to rje_assemblies.py to make available for other tools.
    #i# Assembly maps are stored in self.list['qry'] and self.list['hit']
    #i# Form: ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
    #i# Saved to a *.txt file with one line per sequence and a "seqname = " prefix
    #i# Also save the contig sequences themselves to a matching *.contigs.fasta file.
    #i# Have methods for:
    # - (1) Generating the *.txt file from the map, renaming where required.
    # - (2) Taking a *.txt file and *.contigs.fasta file and making a new fasta file.
    #i# Edit this assembly map with the corrections rather than the actual assembly -> then make new fasta
    #i# Can re-map flanks and thus define any new assembly in this context too.
#########################################################################################################################
    def makeAssemplyMap(self,qh):   ### Populates self.list[qh] with the assembly map.
        '''
        Populates self.list[qh] with the assembly map. Assembly map will be in the form:
        ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            qhbase = '{0}.{1}'.format(basefile,qh)
            self.list[qh] = []  # Assembly map
            self.headLog('{0} Assembly Map generation'.format(qh),line='~')
            ### ~ [1] Generate Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.dbTable(qh,'contigs')
            gdb = self.dbTable(qh,'gap')
            gdb.index('GapFlank5')
            cx = 0.0; ctot = cdb.entryNum()
            seqname = None
            seqmap = []
            for entry in cdb.entries(sorted=True):
                self.progLog('\r#MAP','Generating assembly map: {0:.1f}%'.format(cx/ctot)); cx += 100.0
                if seqmap and seqname != entry['SeqName']:
                    #self.bugPrint('\n{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                    self.list[qh] += seqmap
                    seqmap = []
                # Generate seqmap
                seqname = entry['SeqName']
                if not seqmap: seqmap = ['|']
                seqmap += [entry['Flank5'],'>',entry['Name'],'>',entry['Flank3']]
                gtype = entry['SynBad'].split('-')[1]
                if gtype == 'End':
                    seqmap += ['|']
                else:
                    seqmap += [':{0}:{1}:'.format(gtype,gdb.indexEntries('GapFlank5',entry['Flank3'])[0]['GapLen'])]
            if seqmap:
                #self.bugPrint('\n{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                self.list[qh] += seqmap
            self.printLog('\r#MAP','Generated {0} assembly map from {1} contigs.'.format(qh,rje.iStr(ctot)))
            return True
        except:
            self.errorLog('%s.makeAssemplyMap error' % self.prog())
            return False
#########################################################################################################################
    def getContigSeqObj(self,qh,reload=False):   ### Returns contigs.fasta seqlist object
        '''
        Returns contigs.fasta seqlist object.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqobj = '{0}.contigs'.format(qh)
            if seqobj in self.obj and self.obj[seqobj] and self.obj[seqobj].seqNum() and not reload: return self.obj[seqobj]
            ### ~ [1] Load contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fasin = '{0}.{1}.fasta'.format(self.baseFile(),seqobj)
            seqcmd = self.cmd_list + ['seqin={0}'.format(fasin),'summarise=T']
            self.obj[seqobj] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            if not self.obj[seqobj].seqNum(): raise IOError('Failed to load sequences from {0}'.format(fasin))
            return self.obj[seqobj]
        except:
            self.errorLog('%s.getContigSeqObj error' % self.prog()); raise
#########################################################################################################################
    def saveAssemblyMaps(self,qh,mapname='map',hidegaps=['Hide']):    ### Saves map text and fasta files for qry or hit data
        '''
        Saves map text and fasta files for qry or hit data.
        >> hidegaps:list [] = Hides gaps with given SynBad gap types.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newacc = self.getStrLC({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            if newacc:
                newacc = self.getStr({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            mapout = '{0}.{1}.{2}.txt'.format(self.baseFile(),qh,mapname)
            if not self.outputAssemblyMap(self.list[qh],mapout,newacc): return False
            fasout = rje.baseFile(mapout) + '.fasta'
            contigs = self.getContigSeqObj(qh)
            if not self.fastaFromAssemblyMap(contigs,mapout,fasout,hidegaps): return False
            return True
        except:
            self.errorLog('%s.saveAssemblyMaps error' % self.prog()); return False
#########################################################################################################################
    def outputAssemblyMap(self,maplist,mapout,newacc=None):   ### Generate assembly text file using assembly map.
        '''
        Generate assembly text file using assembly map. Will regenerate names and change gap lengths if required.
        This file will then be used with a contigs fasta file to generate an updated assembly fasta output.
        >> maplist:list = ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
        >> mapout:str = output file name for assembly
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Map output',line='~')
            #i# Sort scaffolds by length.
            #i# If not newacc, will use the first seqname from that assembly map sequence.
            #i# Where there is a clash, will add a .X counter.
            # gaplen=INT : Set new standardised gap length for SynBad assembly output [500]
            newgaplen = max(0,self.getInt('GapLen'))

            ### ~ [1] Generate scaffolds from map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds = []  # List of (scafflen, tempname, mapstr) for each scaffold.
            maptxt = ' '.join(maplist)
            #self.bugPrint(maptxt)
            scaffmaps = maptxt.split('|')
            namecount = {}
            sx = 0.0; stot = len(scaffmaps)
            #self.bugPrint('Assembly map split into {0} scaffold chunks'.format(stot))
            for scaff in scaffmaps:
                smap = scaff.split()
                if not scaff or not smap:
                    stot -= 1
                    continue
                self.progLog('\r#SCAFF','Parsing scaffolds from assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                #self.bugPrint(scaff)
                #self.bugPrint(smap[0])
                scaffname = '.'.join(smap[0].split('.')[:-1])
                #self.bugPrint(scaffname)
                scafflen = 0
                #i# Calculate scaffold length, fix gap lengths and check formatting
                #i# Each cycle should be Flank1,'>',CtgName,'>',Flank2 then either ':SynBad:GapLen:' or end
                i = 0
                while i < len(smap):
                    mapel = ' '.join(smap[i:i+6])
                    try:
                        #!# Move check for coherent map elements into a function -> check edits too
                        #i# Check Dirn
                        if smap[i+1] != smap[i+3] or smap[i+1] not in '<>': raise ValueError('Assembly map strand formatting error!')
                        contig = smap[i+2]
                        cspan = list(map(int,contig.split('.')[-1].split('-')))
                        scafflen += (cspan[1] - cspan[0] + 1)
                        #i# Gap?
                        i += 5
                        if i >= len(smap): break
                        part = smap[i]
                        if part[:1] != ':': raise ValueError('Assembly map gap formatting error!')
                        gap = part.split(':')
                        if len(gap) != 4: raise ValueError('Assembly map gap split formatting error!')
                        gaplen = int(gap[2])
                        if newgaplen and newgaplen != gaplen:
                            gap[2] = str(gaplen)
                            smap[i] = ':'.join(gap)
                            gaplen = newgaplen
                        scafflen += gaplen
                        i += 1
                    except:
                        self.errorLog('Assembly map generation error')
                        raise ValueError('Problem with map element: {0}'.format(mapel))
                scaffstr = ' '.join(smap)
                scaffolds.append((scafflen,scaffname,scaffstr))
                if scaffname not in namecount: namecount[scaffname] = 0
                namecount[scaffname] += 1
            self.printLog('\r#SCAFF','Parsing {0} scaffolds from assembly map complete'.format(rje.iStr(stot)))
            if len(scaffolds) != stot: raise ValueError('Scaffold count mismatch!')
            allnames = rje.sortKeys(namecount)

            ## ~ [1b] Sort scaffolds by length and rename if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maxx = max(namecount.values())
            for name in list(namecount.keys()):
                if namecount[name] == 1: namecount.pop(name)
            scaffolds.sort()
            sorted = scaffolds
            scaffolds = []
            while sorted:
                if newacc:
                    nx = len(sorted)
                    newname = '{0}{1}'.format(newacc,rje.preZero(nx,len(scaffolds)))
                    while newname in allnames:
                        nx += 1
                        newname = '{0}{1}'.format(newacc, rje.preZero(nx, len(scaffolds)))
                    allnames.append(newname)
                else:
                    newname = sorted[0][1]
                    nx = 0
                    while sorted[0][1] in namecount and newname in allnames:
                        accx = rje.preZero(namecount[sorted[0][1]]+nx,maxx); nx += 1
                        newname = '{0}{1}'.format(sorted[0][1],accx)
                    if nx:
                        namecount[sorted[0][1]] -= 1
                        allnames.append(newname)
                scaff = sorted.pop(0)
                scaffolds.append((newname,scaff[0],scaff[2]))
            self.printLog('#SCAFF','{0} scaffolds sorted and renamed.'.format(rje.iLen(scaffolds)))

            ## ~ [1c] Sort scaffolds by name and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds.sort()
            sx = 0.0
            rje.backup(self,mapout)
            OUT = open(mapout,'w')
            for scaff in scaffolds:
                self.progLog('\r#MAP','Outputting assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                #self.bugPrint('{0} ({1}) = {2}\n'.format(scaff[0],scaff[1],scaff[2]))
                OUT.write('{0} ({1}) = {2}\n'.format(scaff[0],scaff[1],scaff[2]))
            OUT.close()
            self.printLog('\r#MAP','Outputted assembly map -> {0}'.format(mapout))

            return True
        except:
            self.errorLog('%s.outputAssemblyMap error' % self.prog())
            return False
#########################################################################################################################
    def saveTelociraptorMaps(self,qh,mapname='map'):    ### Saves map text and fasta files for qry or hit data
        '''
        Saves map text and fasta files for qry or hit data.
        >> hidegaps:list [] = Hides gaps with given SynBad gap types.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newacc = self.getStrLC({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            if newacc:
                newacc = self.getStr({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            mapout = '{0}.{1}.{2}.telociraptor.txt'.format(self.baseFile(),qh,mapname)
            if not self.outputTelociraptorMap(self.list[qh],mapout,newacc): return False
            return True
        except:
            self.errorLog('%s.saveTelociraptorMaps error' % self.prog()); return False
#########################################################################################################################
    def outputTelociraptorMap(self,maplist,mapout,newacc=None):   ### Generate assembly text file using assembly map.
        '''
        Generate assembly text file using assembly map. Will regenerate names and change gap lengths if required.
        This file will then be used with Telociraptor to generate an updated assembly fasta output.
        >> maplist:list = ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
        >> mapout:str = output file name for assembly
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Map output',line='~')
            #i# Sort scaffolds by length.
            #i# If not newacc, will use the first seqname from that assembly map sequence.
            #i# Where there is a clash, will add a .X counter.
            # gaplen=INT : Set new standardised gap length for SynBad assembly output [500]
            newgaplen = max(0,self.getInt('GapLen'))

            ### ~ [1] Generate scaffolds from map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds = []  # List of (scafflen, tempname, mapstr) for each scaffold.
            maptxt = ' '.join(maplist)
            scaffmaps = maptxt.split('|')
            namecount = {}
            sx = 0.0; stot = len(scaffmaps)
            for scaff in scaffmaps:
                smap = scaff.split()
                tmap = []
                if not scaff or not smap:
                    stot -= 1
                    continue
                self.progLog('\r#SCAFF','Parsing scaffolds from assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                scaffname = '.'.join(smap[0].split('.')[:-1])
                scafflen = 0
                #i# Calculate scaffold length, fix gap lengths and check formatting
                #i# Each cycle should be Flank1,'>',CtgName,'>',Flank2 then either ':SynBad:GapLen:' or end
                i = 0
                while i < len(smap):
                    mapel = ' '.join(smap[i:i+6])
                    try:
                        #i# Check Dirn
                        if smap[i+1] != smap[i+3] or smap[i+1] not in '<>': raise ValueError('Assembly map strand formatting error!')
                        contig = smap[i+2]
                        cspan = list(map(int,contig.split('.')[-1].split('-')))
                        scafflen += (cspan[1] - cspan[0] + 1)
                        if smap[i+1] == ">":
                          tmap.append('{}:{}-{}:+'.format(scaffname,cspan[0],cspan[1]))
                        else:
                          tmap.append('{}:{}-{}:-'.format(scaffname,cspan[0],cspan[1]))
                        #i# Gap?
                        i += 5
                        if i >= len(smap): break
                        part = smap[i]
                        if part[:1] != ':': raise ValueError('Assembly map gap formatting error!')
                        gap = part.split(':')
                        if len(gap) != 4: raise ValueError('Assembly map gap split formatting error!')
                        gaplen = int(gap[2])
                        if newgaplen and newgaplen != gaplen:
                            gap[2] = str(gaplen)
                            smap[i] = ':'.join(gap)
                            gaplen = newgaplen
                        tmap.append('~{}:{}~'.format(gap[1],gaplen))
                        scafflen += gaplen
                        i += 1
                    except:
                        self.errorLog('Assembly map generation error')
                        raise ValueError('Problem with map element: {0}'.format(mapel))
                scaffstr = '|'.join(tmap)
                scaffolds.append((scafflen,scaffname,scaffstr))
                if scaffname not in namecount: namecount[scaffname] = 0
                namecount[scaffname] += 1
            self.printLog('\r#SCAFF','Parsing {0} scaffolds from assembly map complete'.format(rje.iStr(stot)))
            if len(scaffolds) != stot: raise ValueError('Scaffold count mismatch!')
            allnames = rje.sortKeys(namecount)

            ## ~ [1b] Sort scaffolds by length and rename if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maxx = max(namecount.values())
            for name in list(namecount.keys()):
                if namecount[name] == 1: namecount.pop(name)
            scaffolds.sort()
            sorted = scaffolds
            scaffolds = []
            while sorted:
                if newacc:
                    nx = len(sorted)
                    newname = '{0}{1}'.format(newacc,rje.preZero(nx,len(scaffolds)))
                    while newname in allnames:
                        nx += 1
                        newname = '{0}{1}'.format(newacc, rje.preZero(nx, len(scaffolds)))
                    allnames.append(newname)
                else:
                    newname = sorted[0][1]
                    nx = 0
                    while sorted[0][1] in namecount and newname in allnames:
                        accx = rje.preZero(namecount[sorted[0][1]]+nx,maxx); nx += 1
                        newname = '{0}{1}'.format(sorted[0][1],accx)
                    if nx:
                        namecount[sorted[0][1]] -= 1
                        allnames.append(newname)
                scaff = sorted.pop(0)
                scaffolds.append((newname,scaff[0],scaff[2]))
            self.printLog('#SCAFF','{0} scaffolds sorted and renamed.'.format(rje.iLen(scaffolds)))

            ## ~ [1c] Sort scaffolds by name and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds.sort()
            sx = 0.0
            rje.backup(self,mapout)
            OUT = open(mapout,'w')
            for scaff in scaffolds:
                self.progLog('\r#MAP','Outputting assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                #self.bugPrint('{0} ({1}) = {2}\n'.format(scaff[0],scaff[1],scaff[2]))
                OUT.write('||{0} ({1} bp)>>{2}<<\n'.format(scaff[0],scaff[1],scaff[2]))
            OUT.close()
            self.printLog('\r#MAP','Outputted Telociraptor assembly map -> {0}'.format(mapout))

            return True
        except:
            self.errorLog('%s.outputTelociraptorMap error' % self.prog())
            return False
#########################################################################################################################
    def fastaFromAssemblyMap(self,contigs,maptxt,fasout,hidegaps=['Hide']):   ### Generate assembly fasta from assembly map and contig sequences.
        '''
        Generate assembly fasta from assembly map and contig sequences.
        >> contigs:SeqList = contig sequences.
        >> maptxt:str = Assembly map text file.
        >> fasout:str = output file name.
        >> hidegaps:list [] = Hides gaps with given SynBad gap types.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Map fasta generation',line='~')
            ctgdict = {}    # {name:sequence}
            for seq in contigs.seqs():
                ctgdict[contigs.shortName(seq)] = contigs.getSeq(seq)[1]
            self.printLog('#CTG','Built sequence dictionary for {0} contigs.'.format(rje.iStr(contigs.seqNum())))
            rje.backup(self,fasout)

            ### ~ [1] Load Map and output fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            FASOUT = open(fasout,'w'); sx = 0
            for mapline in open(maptxt,'r').readlines():
                self.progLog('\r#MAPFAS','Converting assembly map to fasta: {0} scaffolds'.format(rje.iStr(sx)))
                mapdata = mapline.split(' = ')
                if len(mapdata) < 2: continue
                elif len(mapdata) > 2: raise ValueError('Input mapping format error!')
                scaffname = mapdata[0].split()[0]
                scafflen = rje.matchExp('\((\d+)\)',mapdata[0].split()[-1])
                if scafflen: scafflen = int(scafflen[0])
                scaffseq = ''
                smap = mapdata[1].split()
                smap = rje.join(smap)
                smap = rje.split(smap)
                #self.deBug('\n"{0}" "{1}" "{2}" ... "{3}" "{4}"'.format(smap[0],smap[1],smap[2],smap[-2],smap[-1]))
                i = 0; cx = 0
                while i < len(smap):
                    if not smap[i] or not rje.matchExp('(\S)',smap[i]): i += 1; continue
                    try:
                        #i# Check Dirn
                        if smap[i+1] != smap[i+3] or smap[i+1] not in '<>': raise ValueError('Assembly map formatting error!')
                        ctgseq = ctgdict[smap[i+2]]
                        if smap[i+1] == '<': ctgseq = rje_sequence.reverseComplement(ctgseq)
                        scaffseq += ctgseq; cx += 1
                        #i# Gap?
                        i += 5
                        if i >= len(smap): break
                        if not smap[i] or not rje.matchExp('(\S)',smap[i]): i += 1; continue
                        if not smap[i].startswith(':'):
                            self.warnLog('Map assembly element should be gap but is not: "{0}"'.format(smap[i]))
                            i += 1
                            continue
                        gap = smap[i].split(':')
                        if len(gap) != 4: raise ValueError('Assembly map formatting error!')
                        gaplen = int(gap[2])
                        if gap[1] in hidegaps:
                            repn = gaplen // 9
                            extran = gaplen - (repn * 9)
                            leftn = extran // 2
                            rightn = extran - leftn
                            gapseq = 'N' * leftn + 'NNNNGNNNN' * repn + 'N' * rightn
                            if len(gapseq) != gaplen: raise ValueError()
                            scaffseq += gapseq
                        else:
                            scaffseq += 'N' * gaplen
                        i += 1
                        #self.bugPrint('{0}: "{1}"'.format(i,' '.join(smap[i-6:i])))
                    except:
                        #self.bugPrint('\n{0}: "{1}"'.format(i,' '.join(smap[i-6:i])))
                        self.errorLog('Assembly map processing error')
                        raise ValueError('Problem with map element: "{0}"'.format(' '.join(smap[i:i+6])))
                #i# Check scaffold length
                if scafflen and scafflen != len(scaffseq):
                    self.warnLog('Scaffold {0} length mismatch: map says {1} bp but sequence parsed is {2} bp'.format(scaffname,rje.iStr(scafflen),rje.iLen(scaffseq)))
                #!# Add checking of length and filtering if too small
                if self.getInt('MinScaffLen') > len(scaffseq):
                    self.printLog('#MINLEN','Scaffold {0} rejected: {1} < minscafflen={2}'.format(scaffname,rje_seqlist.dnaLen(len(scaffseq)),rje_seqlist.dnaLen(self.getInt('MinScaffLen'))))
                    continue
                FASOUT.write('>{0} {1} bp ({2} contigs)\n{3}\n'.format(scaffname,len(scaffseq),cx,scaffseq)); sx += 1
            self.printLog('\r#MAPFAS','Converted assembly map to fasta: {0} scaffolds -> {1}'.format(rje.iStr(sx),fasout))

            ### ~ [2] Summarise output fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list + ['seqin={0}'.format(fasout),'summarise=T']
            rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])

            return True
        except:
            self.errorLog('%s.fastaFromAssemblyMap error' % self.prog())
            return False
#########################################################################################################################
    ### <9> ### HiC mapping methods                                                                                     #
#########################################################################################################################
    def synBadHiCMap(self):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.flankPairs(): raise ValueError('Failed to generate flank pairs from contigs')
            base = {'qry':rje.baseFile(self.getStr('Genome1'),strip_path=True),'hit':rje.baseFile(self.getStr('Genome2'),strip_path=True)}
            if not self.getBool('PureFlanks'):
                self.setBool({'PureFlanks':True})
                self.warnLog('pureflanks=F not currently supported: setting pureflanks=T')
                #!# Will just need to extract full BED file for this mode
            db = self.db()
            self.seqObjSetup()
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            hicbam = {}
            hicbam['qry'] = self.getStr('HiCBAM1')
            hicbam['hit'] = self.getStr('HiCBAM2')
            bamdir = {}
            bamdir['qry'] = self.getStr('HiCDir1')
            bamdir['hit'] = self.getStr('HiCDir2')
            base = {'qry':rje.baseFile(self.getStr('Genome1'),strip_path=True),'hit':rje.baseFile(self.getStr('Genome2'),strip_path=True)}
            bamout = True
            if os.popen('samtools --version').read():
                self.printLog('#SYS',' '.join(os.popen('samtools --version').read().split()))
            else:
                self.printLog('Cannot open samtools: no HiCBAM analysis')
                bamout = False

            ## ~ [0a] Check and load outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            wanted = []  #X#'qry.tdt','hit.tdt']
            for qh in ['qry','hit']:
                for tname in ['hicpairs','hicbest']:
                    wanted.append('{0}.{1}.tdt'.format(qh,tname))
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main HiC Map files found: loading (force=F)')
                for qh in ['qry','hit']:
                    pdb = db.addTable(name='{0}.hicpairs'.format(qh),mainkeys=['Flank1','Flank2'],expect=True,ignore=[],replace=True)
                    pdb.dataFormat(dbformats)
                    pdb = db.addTable(name='{0}.hicpairs'.format(qh),mainkeys=['Flank1','Flank2'],expect=True,ignore=[],replace=True)
                    pdb.dataFormat(dbformats)
                return True
            if not self.force() and not self.dev():
                self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                self.setBool({'Force':True})
            pcount = {}         # HiCPairs counter used for reloading data


            ### ~ [2] Extract regions to BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                bed = self.dbTable(qh,'flanks')
                gdb = self.dbTable(qh,'gap')
                flankpairs = self.dict['FlankPairs'][qh]     # Dictionary of contig flanks pairs for internal paired read removal
                #!# If not flankpairs will need to reload contigs table
                #!# => Make this a method


                extract = True  # Whether to extract reads from BAM file
                hicdir = True   # Whether HiC directory given to find existing files
                if not bamout: extract = False
                if hicbam[qh].lower() in ['','none']: hicbam[qh] = '{0}.HiC.bam'.format(base[qh])
                if not rje.exists(hicbam[qh]):
                    self.printLog('#HICBAM','HiC BAM file for {0} not found: {1}'.format(qh,hicbam[qh]))
                    extract = False
                if bamdir[qh].lower() in ['','none']: bamdir[qh] = '{0}.{1}flanks/'.format(db.baseFile(),qh)
                if rje.exists(bamdir[qh]):
                    self.printLog('#HICDIR','HiC extraction directory for {0} found: {1}'.format(qh,bamdir[qh]))
                else:
                    self.printLog('#HICDIR','HiC extraction directory for {0} not found: {1}'.format(qh,bamdir[qh]))
                    hicdir = False
                    if extract: rje.mkDir(self,bamdir[qh])
                if extract:
                    self.printLog('#HICBAM','HiC BAM file for {0} found: {1}'.format(qh,hicbam[qh]))
                elif hicdir and not self.force():
                    self.printLog('#HICDIR','HiC BAM file for {0} not found but HiC directory exists for reading previous results: {1}'.format(qh,bamdir[qh]))
                elif not hicdir:
                    self.printLog('#HIC','No HiC BAM file or directory for {0}: no HiC best pair analysis'.format(qh,bamdir[qh]))
                    continue
                self.setBool({'FullHiC':False})
                pcount[qh] = 0
                pdb = None
                if not self.force():
                    pdb = self.addTable(qh,'hicpairs')
                    #pdb = self.dbTable(qh,'hicpairs',expect=False)
                    if not pdb:
                        pdb = db.addTable(name='{0}.hicpairs'.format(qh),mainkeys=['Flank1','Flank2'],expect=False,ignore=[],replace=True)
                        if pdb: pdb.dataFormat(dbformats)
                if pdb and pdb.entryNum() > 0:
                    pcount[qh] = pdb.entryNum()
                elif self.dict['FlankMap'][qh]:
                    self.printLog('#FLANK','{0} flank mapping: skipping read extraction'.format(qh))
                    if not pdb:
                        pdb = db.addEmptyTable('{0}.hicpairs'.format(qh),['Flank1','Flank2','ID1','ID2','Pairs','Type','Score','WTScore'],['Flank1','Flank2'])
                else:
                    if not pdb:
                        pdb = db.addEmptyTable('{0}.hicpairs'.format(qh),['Flank1','Flank2','ID1','ID2','Pairs','Type','Score','WTScore'],['Flank1','Flank2'])
                    uniqfile = '{0}{1}.uniq.id'.format(bamdir[qh],qh)
                    if rje.exists(uniqfile) and not self.force():
                        uniqid = rje.split(open(uniqfile,'r').read())
                        self.printLog('\r#HICBAM','{0} {1} BAM flank read IDs read from {2} (force=F).'.format(rje.iLen(uniqid),qh,uniqfile))
                    elif extract:
                        bx = 0; fx = 0; tx = 0; ex = 0.0; etot = bed.entryNum()
                        readids = []
                        for entry in bed.entries():
                            self.progLog('\r#HICBAM','Extracting {0} BAM flank regions and reads: {1:.1f}%'.format(qh,ex/etot)); ex += 100.0
                            region = '{0}:{1}-{2}'.format(entry['SeqName'],entry['Start'],entry['End'])
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['Name'])
                            sampipe = "samtools view %s %s | awk '{print $1;}' | sort | uniq -c | awk '$1 == 1' | awk '{print $2;}' | tee %s" % (hicbam[qh],region,regtmp)
                            regids = os.popen(sampipe).readlines()
                            readids += regids
                            #X#entry['Score'] = min(1000,len(regids))
                            if rje.exists(regtmp): tx += 1
                            else:
                                self.warnLog('Read ID tmp file {0} not created'.format(regtmp))
                                self.debug(sampipe)
                            if self.dev():  #?# Add as options?
                                regbam = '%s%s.bam' % (bamdir[qh],entry['Name'])
                                fasout = '%s%s.reads.fasta' % (bamdir[qh],entry['Name'])
                                samcmd = 'samtools view -h -F 0x104 {0} {1} | samtools view -h -b -o {2} - 2>&1'.format(hicbam[qh],region,regbam)
                                os.popen(samcmd).read()
                                if rje.exists(regbam): bx += 1
                                else: self.warnLog('BAM file {0} not created'.format(regbam))
                                bam2fas = 'samtools fasta {0} > {1} 2>&1'.format(regbam,fasout)
                                os.popen(bam2fas).read()
                                #!# Parse processed 3452 reads
                                if rje.exists(fasout): fx += 1
                                else: self.warnLog('Read fasta file {0} not created'.format(fasout))
                        self.printLog('\r#HICBAM','Extracted {0} read IDs from {1} BAM flank regions: {2} tmp files'.format(rje.iLen(readids),qh,rje.iStr(tx)))
                        if self.dev():
                            self.printLog('\r#HICBAM','Extracted {0} BAM flank regions and reads: {1} BAM and {2} fasta files'.format(qh,rje.iStr(bx),rje.iStr(fx)))
                        ## ~ Identify region-identical paired IDs ~ #
                        #!# NOTE: This would not work with pureflanks=F
                        tmpfile = '{0}{1}.tmp'.format(bamdir[qh],qh)
                        open(tmpfile,'w').writelines(readids)
                        uniqid = rje.split(os.popen("sort %s | uniq -c | awk '$1 == 1' | awk '{print $2;}' | tee %s" % (tmpfile,uniqfile)).read())
                        self.printLog('\r#HICBAM','{0} of {1} {2} BAM flank read IDs found in one flank only.'.format(rje.iLen(uniqid),rje.iLen(readids),qh))
                    ## ~ Reduce to non-unique read IDs per region ~ ##
                    #!# This is very slow. Consider speeding up with diff and/or forking #!#
                    nx = 0; ix = 0; ex = 0.0; etot = bed.entryNum()
                    for entry in bed.entries():
                        self.progLog('\r#HICBAM','Reducing {0} BAM flank IDs to partial pairs: {1:.2f}%'.format(qh,ex/etot)); ex += 100.0
                        regfile = '%s%s.id' % (bamdir[qh],entry['Name'])
                        if entry['Name'] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir[qh],self.dict['FlankMap'][qh][entry['Name']])
                        if rje.exists(regfile) and not self.force():
                            regids = rje.split(open(regfile,'r').read())
                        elif extract:
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['Name'])
                            regids = rje.listDifference(rje.split(open(regtmp,'r').read()),uniqid)
                            if flankpairs[entry['Name']] != entry['Name']:
                                regtmp2 = '%s%s.tmp' % (bamdir[qh],flankpairs[entry['Name']])
                                regids = rje.listDifference(regids,rje.split(open(regtmp2,'r').read()))
                            #regids = rje.split(os.popen("diff %s %s | grep '^<' | awk '{print $2;}' | tee %s" % (regtmp,uniqfile,regfile)).read())
                            open(regfile,'w').write('\n'.join(regids))
                        else:
                            self.warnLog('No ID file for {0} and no BAM for extraction'.format(entry['Name']))
                            regids = []
                        if not regids: nx += 1
                        ix += len(regids)
                        entry['Score'] = min(1000,len(regids))
                    self.printLog('\r#HICBAM','Reduced {0} BAM flank IDs to partial pairs -> {1} read IDs; {2} regions without reads'.format(qh,rje.iStr(ix),rje.iStr(nx)))
                    #i# Cleanup
                    if not self.dev() and not self.debugging():
                        for entry in bed.entries():
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['Name'])
                            if rje.exists(regtmp): os.unlink(regtmp)

            ### ~ [3] Calculate pair overlaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                #X#self.dict['HiCOffProb'][qh] = self.calculateHiCOffProb(qh)
                for field in ['HiCScore','BestFlank5','BestFlank3']: gdb.addField(field)
                hicscore = self.getStr('HiCScore')
                bscore = self.getStr('HiCScore') #'pairs'
                gx = 0.0; gtot = gdb.entryNum()
                for entry in gdb.entries():
                    self.progLog('\r#HICBAM','Calculating {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,gx/gtot)); gx += 100.0
                    regions = [entry['GapFlank5'],entry['GapFlank3']]
                    pentry = self.hicPair(regions,qh,pdb,entry['SynBad'])
                    entry['HiCScore'] = pentry[hicscore]
                    #?# Should this be moved to be part of the correction method?
                    if entry['HiCScore'] == 'NA':
                        entry['BestFlank5'] = '?'
                        entry['BestFlank3'] = '?'
                        continue
                    # elif entry['HiCScore'] >= 0.5:
                    #     entry['BestFlank3'] = entry['GapFlank3']
                    #     entry['BestFlank5'] = entry['GapFlank5']
                    #     flankcheck = []
                    elif entry['SynBad'] in ['Syn', 'Aln', 'Long', 'InvFix', 'Dup', 'InvDupFix'] and self.getStrLC('HiCMode') != 'full':
                        #i# For the same-sequence gaps, just compare to others in same sequence
                        #i# Note that the reciprocal search should still identify a better match for out of place contigs
                        flankcheck = bed.indexEntries('SeqName',entry['SeqName'])
                    else:
                        flankcheck = bed.entries()
                    #i# 5' Join
                    bestscore = pentry[bscore] #X# entry['HiCScore']
                    bestjoin = entry['GapFlank5']
                    for jentry in flankcheck:
                        regions = [entry['GapFlank3'],jentry['Name']]
                        centry = self.hicPair(regions,qh,pdb,'Check')
                        if centry and centry[hicscore] != 'NA' and centry[bscore] > bestscore:
                            bestscore = centry[bscore]
                            bestjoin = jentry['Name']
                    entry['BestFlank5'] = bestjoin
                    #i# 3' Join
                    bestscore = pentry[bscore] #X#entry['HiCScore']
                    bestjoin = entry['GapFlank3']
                    for jentry in flankcheck:
                        regions = [entry['GapFlank5'],jentry['Name']]
                        centry = self.hicPair(regions,qh,pdb,'Check')
                        if centry and centry[hicscore] != 'NA' and centry[bscore] > bestscore:
                            bestscore = centry[bscore]
                            bestjoin = jentry['Name']
                    entry['BestFlank3'] = bestjoin
                self.printLog('\r#HICBAM','Calculating {1} {0} BAM flank ID pair overlaps complete'.format(qh,rje.iStr(pdb.entryNum())))

                #i# Update flanks Score from hicpairs table to have ID count
                hdb = self.dbTable(qh,'hicpairs')
                idict = {}
                for entry in hdb.entries():
                    idict[entry['Flank1']] = entry['ID1']
                    idict[entry['Flank2']] = entry['ID2']
                fdb = self.dbTable(qh,'flanks')
                for entry in fdb.entries():
                    if entry['Name'] in idict: entry['Score'] = idict[entry['Name']]

            ### ~ [4] Add random/full data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                #i# NOTE: Random pairs are very unlikely to share IDs. Need a good way to score match and then compare.
                #i# Need to take into account the number of IDs in each flank. Naively assume half match in each direction.
                #?# Should we have some positive control regions of the adjacent regions within a chunk
                #i# NOTE: Duplication gaps are going to cause issues due to repeated sequences?
                #!# Check the setting for HiCMode
                if self.getStrLC('HiCMode') == 'random':
                    regnames = list(bed.index('Name').keys())
                    #!# Add some kind of safety check, e.g. if total combinations < 2x replicates -> do them all!
                    replicates = 10000
                    randx = 0
                    reglist1 = rje.randomList(regnames)
                    reglist2 = rje.randomList(regnames)
                    while randx < replicates:
                        self.progLog('\r#RANDOM','Generating random {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,100.0*randx/replicates))
                        if not reglist1:
                            reglist1 = rje.randomList(regnames)
                            reglist2 = rje.randomList(regnames)
                        regions = [reglist1.pop(0),reglist2.pop(0)]
                        self.hicPair(regions,qh,pdb,'Rand')
                        randx += 1
                    self.printLog('\r#RANDOM','Generated {1} random {0} BAM flank ID pair overlaps.'.format(qh,rje.iStr(replicates)))
                elif self.getStrLC('HiCMode') == 'full':
                    regnames = list(bed.index('Name').keys())
                    replicates = len(regnames)
                    randx = 0
                    for region1 in regnames:
                        self.progLog('\r#RANDOM','Generating full {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,100.0*randx/replicates))
                        for region2 in regnames:
                            regions = [region1,region2]
                            self.hicPair(regions,qh,pdb,'Rand')
                        randx += 1
                    self.printLog('\r#RANDOM','Generated {1} full BAM flank ID pair overlaps.'.format(qh,rje.iStr(replicates)))
                    self.setBool({'FullHiC':True})

            ### ~ [5] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Without HiC data, this will only have partial data
            for qh in ['qry','hit']:
                for tname in ['hicpairs']:
                    table = self.dbTable(qh,tname,expect=False)
                    if not table: continue
                    if tname == 'hicpairs' and table and pcount[qh] == table.entryNum():
                        self.printLog('#SAVE','Skipping saving of {0} table - no new pairs'.format(tname))
                        continue
                    elif tname == 'hicpairs': #!# Add hiczero=T/F
                        self.progLog('\r#DROP','Dropping zero-overlaps from {0}...'.format(tname))
                        table.dropEntriesDirect('Pairs',[0])
                    table.saveToFile()
                self.bestHiCPairTable(qh)
            return True
        except:
            self.errorLog('%s.synBadHiCMap error' % self.prog())
            return False
#########################################################################################################################
    def calculateHiCOffProb(self,qh,pdb=None):  ### Returns the calculated probability of a HiC read not being in best pair
        '''
        Returns the hicpairs entry for pair of regions.
        >> qh:str qry or hit
        >> pdb:Table hicpairs table (will grab if None)
        :return: prob
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Not currently being used.
            if not pdb: pdb = self.dbTable(qh,'hicpairs')
            gdb = self.dbTable(qh,'gap')
            gx = 0.0; gtot = gdb.entryNum()
            idsum = 0; pairsum = 0
            for entry in gdb.entries():
                self.progLog('\r#HICBAM','Calculating Syn/Aln {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,gx/gtot)); gx += 100.0
                if entry['SynBad'] not in ['Syn', 'Aln']: continue
                regions = [entry['GapFlank5'],entry['GapFlank3']]
                pentry = self.hicPair(regions,qh,pdb,entry['SynBad'],score=False)
                idsum += pentry['ID1']
                idsum += pentry['ID2']
                pairsum += 2 * pentry['Pairs']
            self.printLog('\r#HICBAM','Calculation of Syn/Aln {0} BAM flank ID pair overlaps complete.')
            return float(idsum - pairsum) / float(idsum)
        except:
            self.errorLog('%s.calculateHiCOffProb error' % self.prog())
            raise
#########################################################################################################################
    def hicPair(self,regions,qh,pdb=None,ptype='',score=True):  ### Returns the hicpairs entry for pair of regions.
        '''
        Returns the hicpairs entry for pair of regions.
        >> regions:list [flank1,flank2]
        >> qh:str qry or hit
        >> pdb:Table hicpairs table (will grab if None)
        >> ptype:str [''] = Pair type for entry
        >> score:bool [True] = Whether to calculate scores and add to pdb table
        :return: entry
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if regions[0] == regions[1]: return None
            if not pdb: pdb = self.dbTable(qh,'hicpairs')
            bamdir = {'qry':self.getStr('HiCDir1'),'hit':self.getStr('HiCDir2')}[qh]
            if bamdir.lower() in ['','none']: bamdir = '{0}.{1}flanks/'.format(pdb.baseFile(),qh)
            regions.sort()
            pentry = pdb.data((regions[0],regions[1]))
            if pentry: return pentry
            if self.getBool('FullHiC'): return None
            regfile = '%s%s.id' % (bamdir,regions[0])
            if regions[0] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir,self.dict['FlankMap'][qh][regions[0]])
            idlist1 = rje.split(open(regfile,'r').read())
            regfile = '%s%s.id' % (bamdir,regions[1])
            if regions[1] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir,self.dict['FlankMap'][qh][regions[1]])
            idlist2 = rje.split(open(regfile,'r').read())
            pairs = rje.listIntersect(idlist1,idlist2)
            pentry = {'Flank1':regions[0],'Flank2':regions[1],'ID1':len(idlist1),'ID2':len(idlist2),'Pairs':len(pairs),'Type':ptype,'WTScore':0.0,'Score':0.0}
            if not score: return pentry
            #X#p = self.dict['HiCOffProb'][qh]
            if pentry['Pairs']:
                pentry['Score'] = 2.0 * pentry['Pairs'] / (pentry['ID1']+pentry['ID2'])
                pentry['WTScore'] = (0.5 * pentry['Pairs'] / pentry['ID1']) + (0.5 * pentry['Pairs'] / pentry['ID2'])
                #?# pentry['pscore'] = -> Can we use rje.poisson(observed,expected,exact=False,callobj=None,uselog=True) ? How?
            elif min(pentry['ID1'],pentry['ID2']) < 1:
                pentry['Score'] = 'NA'
                pentry['WTScore'] = 'NA'
            elif self.getBool('FullHiC'): return None
            pdb.addEntry(pentry)
            return pentry
        except:
            self.errorLog('%s.hicPair error' % self.prog())
            return None
#########################################################################################################################
    def doubleFlanks(self,qh):  ### Returns a list of flank regions that are both ends of a contig (have two "best" pairs)
        '''
        Returns a list of flank regions that are both ends of a contig (have two "best" pairs).
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qh = qh.lower()
            if not 'DoubleFlanks' in self.dict: self.dict['DoubleFlanks'] = {'qry':[],'hit':[]}
            if self.dict['DoubleFlanks'][qh]: return self.dict['DoubleFlanks'][qh]
            cdb = self.dbTable(qh,'contigs')

            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Rename flanks as "Ends" for contigs = only flanks for Gaps
            #!# "Pairs" are ends that match another end over a gap (i.e. flank the same gap)
            for entry in cdb.entries():
                if entry['Flank5'] == entry['Flank3']: self.dict['DoubleFlanks'][qh].append(entry['Flank5'])
            return self.dict['DoubleFlanks'][qh]
        except:
            self.errorLog('%s.doubleFlanks error' % self.prog())
            raise
#########################################################################################################################
    def bestHiCPairTable(self,qh):  ### Reduces the HiC pairs table to the best pairs.
        '''
        Reduces the HiC pairs table to the best pairs.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# At some point will want to add up/downstream checks for consistency
            hictable = '{0}.hicbest'.format(qh.lower())
            if not self.dbTable(qh,'hicpairs',expect=False): return None
            hdb = self.db().copyTable(self.dbTable(qh,'hicpairs'),hictable,replace=True,add=True)
            score = self.getStr('HiCScore')
            doubles = self.doubleFlanks(qh)
            #self.debug('{0}: {1}'.format(qh,' '.join(doubles)))

            ### ~ [1] Rate pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdb.addField('Best',evalue='Null')
            allpair = {}   # Dictionary of flank:[(score,id,pair)] -> Sort then assess
            ## ~ [1a] Build pair lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in hdb.entries():
                for flank in ['Flank1','Flank2']:
                    pair = {'Flank1':'Flank2','Flank2':'Flank1'}[flank]
                    if entry[flank] not in allpair: allpair[entry[flank]] = []
                    allpair[entry[flank]].append((entry[score],entry['Pairs'],entry[pair]))
                    #self.bugPrint('{0}: {1}'.format(entry[flank],allpair[entry[flank]]))
            #self.debug('1a||')
            ## ~ [1b] Reduce to best pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bestpair = {}   # Dictionary of flank:[bestpairs]
            weakpair = {}   # Dictionary of flank:[bestpairs that fail to meet HiCMin]
            for flank in allpair:
                allpair[flank].sort(reverse=True)
                bestpair[flank] = []
                weakpair[flank] = []
                #i# Strong first
                for pair in allpair[flank]:
                    if pair[1] >= self.getInt('HiCMin'): bestpair[flank].append(pair[2])
                    else: weakpair[flank].append(pair[2])
                bestn = {True:2,False:1}[flank in doubles]
                bestpair[flank] = bestpair[flank][:bestn]
                weakpair[flank] = weakpair[flank][:(bestn-len(bestpair[flank]))]
                #self.bugPrint('{0}: {1} / {2}'.format(flank,bestpair[flank],weakpair[flank]))
            #self.debug('1b||')
            ## ~ [1c] Update table with best pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in hdb.entries():
                if entry['Type'] not in ['Check','Rand']: entry['Best'] = 'None'
                for flank in ['Flank1','Flank2']:
                    pair = {'Flank1':'Flank2','Flank2':'Flank1'}[flank]
                    if entry[pair] in bestpair[entry[flank]] and entry[flank] in bestpair[entry[pair]]: entry['Best'] = 'Both'
                    elif entry[pair] in weakpair[entry[flank]] and entry[flank] in bestpair[entry[pair]]: entry['Best'] = 'Both'
                    elif entry[pair] in bestpair[entry[flank]] and entry[flank] in weakpair[entry[pair]]: entry['Best'] = 'Both'
                    elif entry[pair] in bestpair[entry[flank]]: entry['Best'] = flank
                    elif entry[flank] in bestpair[entry[pair]]: entry['Best'] = pair
                    elif entry[pair] in weakpair[entry[flank]] and entry[flank] in weakpair[entry[pair]]: entry['Best'] = 'Weak'
                    elif entry[pair] in weakpair[entry[flank]]: entry['Best'] = 'Weak' + flank[-1]
                    elif entry[flank] in weakpair[entry[pair]]: entry['Best'] = 'Weak' + pair[-1]
                #self.bugPrint('{0}'.format(hdb.entrySummary(entry,collapse=True)))
            #self.debug('1c||')

            ### ~ [2] Save table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdb.dropEntriesDirect('Best',['Null'])
            hdb.saveToFile()
            return hdb

        except:
            self.errorLog('%s.bestHiCPairTable error' % self.prog())
            return None
#########################################################################################################################
    def bestHiCPairDict(self,qh):   ### Returns dictionary of {flank:bestpair}
        '''
        Returns dictionary of {flank:bestpair}.
        '''
        besthic = {}
        hicdb = self.dbTable(qh,'hicbest')
        for entry in hicdb.indexEntries('Best','Both'):
            besthic[entry['Flank1']] = entry['Flank2']
            besthic[entry['Flank2']] = entry['Flank1']
        for entry in hicdb.indexEntries('Best','Flank1'):
            if entry['Flank1'] not in besthic:
                besthic[entry['Flank1']] = entry['Flank2']
        for entry in hicdb.indexEntries('Best','Flank2'):
            if entry['Flank2'] not in besthic:
                besthic[entry['Flank2']] = entry['Flank1']
        return besthic
#########################################################################################################################
    ### <10> ### SynBad Compression and Tidy Methods                                                                    #
#########################################################################################################################
    def synBadCompressTable(self,qh,locdb,maxskip=0,quick=False):   ### Main SynBad gap mapping compression for a local hits table.
        '''
        Main SynBad gap mapping compression for a local hits table.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            if 'Alt' not in locdb.fields():
                locdb.addField('Alt',evalue=0,after='Non')
                locdb.addField('Trans',evalue=0,after='Alt')
                locdb.addField('Inv',evalue=0,after='Trans')
                locdb.addField('Dup',evalue=0,after='Inv')
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: quick = True
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = []

            ### ~ [1] Cycle through entries, identifying and merging flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            merged = True
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            while entries:
                self.progLog('\r#CMPRSS','Compress %s table by %s: %.2f%%' % (tname,qh,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if not merged and not lookback:
                        lookback = True
                    elif not merged:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                elif not merged:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break
                elif eskip > 0:
                    # if merging was performed for higher eskip -> will need to reset eskip and cycle again
                    eskip = 0
                    lookback = False
                merged = False

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.': continue #i# Don't merge gaps
                #i# - (2) Lack of internal qry gaps
                if '>' in merge1[qgap] or '<' in merge2[qgap]:
                    if debug: self.bugPrint('QryGap...')
                    continue
                skipgap = False
                for entry in toskip:
                    if entry[qgap] != '.': skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps
                #?# Currently eliminating any internal entries with gaps. May want to refine this.
                for entry in [merge1,merge2] + toskip:
                    if entry[hgap] != '.': skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                if qdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                if hdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                    continue

                ## ~ [1c] Classify and merge intermediates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                if debug: self.bugPrint('-> merge ->')
                merged = True
                mstart = merge1[hend] - self.getInt('MaxOverlap')
                mend = merge2[hstart] + self.getInt('MaxOverlap')
                if bwd:
                    mstart = merge2[hend] - self.getInt('MaxOverlap')
                    mend = merge1[hstart] + self.getInt('MaxOverlap')
                prev = merge1
                merge1['Dup'] += max(0,-hdist)
                for entry in toskip + [merge2]:
                    ex += 100.0
                    if entry == merge2: etype = 'Merge'
                    #i# (i) Different sequence (AltChrom)
                    elif entry[alt] != merge1[alt]: etype = 'Alt'
                    #i# (i) Same seq, collinear
                    elif entry[hstart] > mstart and entry[hend] < mend:
                        if entry['Strand'] == merge1['Strand']: etype = 'Col'
                    #i# (ii) Same seq, inside, inverted (Inv)
                        else: etype = 'Inv'
                    #i# (iii) Same seq, outside (Trans)
                    else: etype = 'Trans'
                    #i# Merge into merge1
                    for field in ['Length','Identity','Alt','Trans','Inv','Dup']:
                        merge1[field] += entry[field]
                    inslen = entry[hend] - entry[hstart] + 1
                    if etype in ['Alt','Trans','Inv']: merge1[etype] += inslen
                    merge1['Non'] += max(0,entry[qstart]-prev[qend]-0)
                    entries.remove(entry)
                    if entry in sorted: sorted.remove(entry)
                    prev = entry

                ## ~ [1d] Update and remove the merged entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                merge1[qend] = merge2[qend]
                if bwd: merge1[hstart] = merge2[hstart]
                else: merge1[hend] = merge2[hend]
                merge1['{0}3'.format(qh)] = merge2['{0}3'.format(qh)]
                if merge1[qgap] != '.' and merge2[qgap] != '.':
                    merge1[qgap] = merge1[qgap] + merge2[qgap]
                elif merge2[qgap] != '.':
                    merge1[qgap] = merge2[qgap]
                locdb.dropEntryList(toskip + [merge2],log=False)
                if debug: self.debug('%s' % (locdb.entrySummary(merge1,collapse=True)))
            self.printLog('\r#CMPRSS','Compressed %s table by %s: %s alignments -> %s' % (tname,qh,rje.iStr(etot),rje.iStr(locdb.entryNum())))
            locdb.remakeKeys()

            return (etot - locdb.entryNum())
        except: self.errorLog('%s.synBadCompressTable error' % self.prog())
#########################################################################################################################
    def synBadCompress(self):   ### Main SynBad gap mapping compression
        '''
        Main SynBad gap mapping compression.
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            self.printLog('#CMPRSS','Maximum distance between local hits for compression: %s' % rje_seqlist.dnaLen(self.getInt('MaxSynSpan')))
            self.printLog('#CMPRSS','Maximum overlap between local hits for compression: %s' % rje_seqlist.dnaLen(self.getInt('MaxOverlap')))
            self.printLog('#CMPRSS','Maximum number of local hits to bypass for compression: %s' % self.getInt('MaxSynSkip'))
            ## ~ [0a] Check for tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #X# wanted = ['qryfull.tdt','qry.tdt','qrypairs.tdt','qrygap.tdt','hitfull.tdt','hit.tdt','hitpairs.tdt','hitgap.tdt']
            #i# Removing pairs from the check: not used for anything!
            wanted = ['qry.full.tdt','qry.tdt','qry.gap.tdt','hit.full.tdt','hit.tdt','hit.gap.tdt']
            if rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                if self.force() and not self.dev():
                    self.printLog('#FORCE','SynBad output files found, but regenerating (force=T).')
                else:
                    for tdtfile in wanted:
                        tname = tdtfile[:-4]
                        if 'full' in tname: continue
                        if 'gap' in tname: continue
                        if self.getBool('Update'):
                            if self.db(tname): db.deleteTable(tname)
                        if tname.startswith('qry'):
                            tkeys = ['Qry','QryStart','QryEnd']
                            if 'pairs' in tname: tkeys.insert(1,'Hit')
                        else:
                            tkeys = ['Hit','HitStart','HitEnd']
                            if 'pairs' in tname: tkeys.insert(1,'Qry')
                        table = db.addTable(mainkeys=tkeys,name=tname,ignore=[],expect=True,replace=True)
                        table.dataFormat(dbformats)
                    for qh in ('Qry','Hit'):
                        self.updateHitGHaps(qh)
                    if self.getBool('Update'):
                        self.printLog('#UPDATE','SynBad output files found and loaded, but updating (force=F update=T).')
                    else:
                        self.printLog('#SKIP','SynBad output files found and loaded (force=F update=F).')
                        return True
            else:
                self.printLog('#CHECK','Not all SynBad output files found: generating')
                if not self.force() and not self.dev():
                    self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                    self.setBool({'Force':True})
            ## ~ [0b] Save full tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                tname = qh.lower()
                locdb = self.db(tname)
                if not locdb:
                    tkeys = ['Hit','HitStart','HitEnd']
                    if tname.startswith('qry'): tkeys = ['Qry','QryStart','QryEnd']
                    locdb = db.addTable(mainkeys=tkeys,name=tname,ignore=[],expect=True,replace=True)
                    locdb.dataFormat(dbformats)
                    self.updateHitGHaps(qh)
                #!# Don't output if qry and qryfull both found?
                locdb.info['Name'] = '{0}.full'.format(qh.lower())
                locdb.saveToFile()
                locdb.info['Name'] = qh.lower()
            ## ~ [0c] Add gaps to tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                gapdb = self.dbTable(qh.lower(),'gap')
                gapdb.dataFormat({'Start':'int','End':'int','GapLen':'int'})
                for gentry in gapdb.entries():
                    if gentry['SynBad'] in ['Aln','Aligned']: continue
                    #  seqname start end seqlen gaplen gap MaxFlank5  MaxFlank3  Span0  Span100  Span1000  Span5000  GapSpan  SynSpan SynBad
                    # Qry           QryStart   QryEnd     Hit           HitStart   HitEnd     Strand  Length   Identity  Non  Alt    Trans   Inv    QryGap  HitGap  Qry5       Qry3      Hit5       Hit3
                    addentry = {qh:gentry['SeqName'],'{0}Start'.format(qh):gentry['Start'],'{0}End'.format(qh):gentry['End'],
                                alt:gentry['SynSpan'],'{0}Start'.format(alt):1,'{0}End'.format(alt):gentry['GapLen'],
                                '{0}Gap'.format(qh):gentry['SynBad'],'{0}Gap'.format(alt):'.','Strand':'.',
                                'Non':gentry['GapLen'],'Alt':0,'Trans':0,'Inv':0,'Dup':0,
                                'Length':gentry['GapLen'],'Identity':0,'Qry5':'-','Qry3':'-','Hit5':'-','Hit3':'-'}
                    locdb.addEntry(addentry,overwrite=False)
                self.printLog('#GAPS','Added {0} {1} gaps to {2} table -> {3} entries'.format(rje.iStr(gapdb.entryNum()),qh,qh.lower(),rje.iStr(locdb.entryNum())))

            ### ~ [1] Quick merge of adjacent flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ('Qry','Hit'):
                locdb = self.dbTable(qh)
                skipx = 0
                self.headLog('SynBad {0} quick compression: merge adjacent local hits'.format(qh),line='~')
                self.synBadCompressTable(qh,locdb,skipx,quick=True)

            ### ~ [2] Remove internal odd entries and merge flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('MaxSynSkip') > 0:
                for qh in ('Qry','Hit'):
                    alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                    locdb = self.dbTable(qh)
                    cyc = 0
                    prex = 0
                    while locdb.entryNum() != prex:
                        prex = locdb.entryNum(); cyc += 1
                        # for skipx in range(self.getInt('MaxSynSkip'),-1,-1):
                        #     self.headLog('SynBad {0} compression cycle {1} (SynSkip={2})'.format(qh,cyc,skipx),line='~')
                        #     self.synBadCompressTable(qh,locdb,skipx)
                        self.headLog('SynBad {0} compression cycle {1} (MaxSynSkip={2})'.format(qh,cyc,self.getInt('MaxSynSkip')),line='~')
                        self.synBadCompressTable(qh,locdb,self.getInt('MaxSynSkip'))
                    continue
            else: self.printLog('#MERGE','Adjacent merging only (maxsynskip=0).')

            return True
        except: self.errorLog('%s.synBadCompress error' % self.prog()); return False
#########################################################################################################################
    def tidyTables(self):   ### SynBad tidying of gaps with clear fixes
        '''
        This method works through a table and identifies regions between pairs of gaps that might possibly be able to be
        fixed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Notes for future updates ...
            #?# Add a gaphide=LIST mode to "hide" gaps with a certain final classification.
            #?# Add a gapshow=FASFILE mode to resinstate hidden gaps in a given file.
            #?# Add a gapignore=LIST option to ignore certain gaps during compression? (Or just encourage hide -> re-run -> show?)
            #i# NOTE: Don't actually want to hide these gaps, I think. It is possible that maxsynskip settings could
            #i# result in a Syn gap that could still contain another mis-placed contig.
            #!# Update the gap tables with a Rating field that gets additional elements like Diploidocus
            #!# SynBad|Span|HiC|FlankDepth|Tidy
            db = self.db()
            for qh in ('qry','hit'):
                self.addTable(qh,'corrections',make=False)
            wanted = ['qry.blocks.tdt','hit.blocks.tdt','qry.corrections.tdt','hit.corrections.tdt']
            if not self.update() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Tidied SynBad *.blocks.tdt and *.corrections.tdt files found: reloading and skipping generation (force=F update=F)')
                for qh in ('qry','hit'):
                    self.addTable(qh,'blocks',expect=True)
                    self.addTable(qh,'corrections',expect=True)
                return True
            if not self.force() and not self.dev():
                self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                self.setBool({'Force':True})

            ### ~ [1] Real translocations -> Divergent/Long Synteny ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Any Tran/Brk gaps that link regions of two scaffolds where scaffold 1 is unambiguously linked up/downstream
            #i# to scaffold 2, should be reclassified as "Div" (Divergent). Make sure the Hit does not have any gaps here too?
            #i# Probably also want to have some kind of proximity filter. Could this be the same as the merging code,
            #i# except re-labelling rather than merging due to the presence of the gap(s)?
            #i# Any Brk gaps that are flanked by collinear hits > maxsynpan bp apart will be re-classified as 'Long'
            for qh in ('Qry','Hit'):
                locdb = self.dbTable(qh)
                gapdb = self.dbTable(qh,'gap')
                self.addTable(qh,'corrections')
                # cordb = self.dbTable(qh,'corrections')
                skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} translocation/breakpoint assessment'.format(qh),line='~')
                self.synBadDivergence(qh,locdb,gapdb,skipx,quick=True)

            ### ~ [2] HiC support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Updates based on best HiC support (BestFlank matching GapFlank)
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                gapdb = self.dbTable(qh,'gap')
                hicdb = self.dbTable(qh,'hicpairs',expect=False)
                if not hicdb or not hicdb.entryNum(): continue
                hx = 0; dx = 0; px = 0
                self.headLog('SynBad {0} HiC pair assessment'.format(qh),line='~')
                #!# Should there be a toggle to allow reclassification of partial HiC support? Would this indicate a duplication problem?!
                for gentry in gapdb.entries():
                    if gentry['BestFlank3'] == gentry['GapFlank3'] and gentry['BestFlank5'] == gentry['GapFlank5']:
                        if gentry['SynBad'] in hicgaps:
                            gentry['SynBad'] = 'HiC'; hx += 1
                            locdb.data((gentry['SeqName'],gentry['Start'],gentry['End']))['QryGap'] = 'HiC'
                        elif gentry['SynBad'] in ['Dup']:
                            gentry['SynBad'] = 'DupHiC'; dx += 1
                            locdb.data((gentry['SeqName'],gentry['Start'],gentry['End']))['QryGap'] = 'DupHiC'
                    elif gentry['SynBad'] not in hicgaps + ['Dup']:
                        px += 1
                    elif gentry['BestFlank3'] == gentry['GapFlank3'] or gentry['BestFlank5'] == gentry['GapFlank5']:
                        self.printLog('#HIC','{0} Gap {1}::{2} has partial best-HiC support'.format(gentry['SynBad'],gentry['GapFlank5'],gentry['GapFlank3']))
                self.printLog('#HIC','{0} {1} gaps with best HiC-scoring flank pairs reclassified to "HiC"'.format(hx,qh))
                self.printLog('#HIC','{0} {1} gaps with best HiC-scoring flank pairs reclassified to "DupHiC"'.format(dx,qh))
                self.printLog('#HIC','{0} {1} gaps with one-directional best HiC-scoring flank pair'.format(px,qh))

            ### ~ [3] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Inversions consist of a region flanked by two Inv QryGaps. In the simplest case, there are no intervening
            #i# gaps, but the following gaps can also be skipped: Syn, Span, Dup. (Aln should not appear!)
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                gapdb = self.dbTable(qh,'gap')
                #skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} inversion assessment'.format(qh),line='~')
                self.synBadInversions(qh,locdb,gapdb,0,quick=True)
                self.synBadInvertedDuplications(qh,locdb,gapdb,0,quick=True)

            ### ~ [4] Update hitgaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ('Qry','Hit'):
                self.updateHitGHaps(qh)

            ### ~ [5] Generate Synteny blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #># Generate a new fragments table that skips through Syn, Aln, Span, InvFix, InvDupFix, Long and Div gaps. Sort by size?
            # and look for out of place neighbours. Build a collinear backbone in size order? Base this on the ends.
            # Identity all that don't fit. Sub classify as cis and trans. Look at whether cis can be moved into any gaps.
            # Need to identify cis and trans gaps too - include the Div gaps as cis, ignoring the non-cis region.
            # This latter bit is basically looking at the pairs tables but including the gaps in every pair?
            # Need to rename the Div gaps so they have the cis scaffolds in their names.
            # Can then add cis gaps to the pairs table and work off that.
            for qh in ('Qry','Hit'):
                self.headLog('SynBad {0} synteny blocks and cis/trans assignment'.format(qh),line='~')
                self.synBadSyntenyBlocks(qh)

            ### ~ [5] HiC "best" translocations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Updates based on best HiC support (BestFlank matching GapFlank)
            for qh in ('Qry','Hit'):
                self.synBadHiCBestTranslocations(qh)



            ### ~ [!] Duplications ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# These need the addition of read depth data. Here, the duplicated regions should be identified and then
            #i# analysed for the relative probability of being <=0.5N versus >=1N. If the region between two gaps,
            #i# including a Dup gap, is <=0.5N then it should be removed: will need to be careful not to remove two!
            #i# If a region is spanned by two Dup gaps then it is the best candidate: focus on these first!
            #?# NOTE: I think this needs to be a different tool. SynBad should output plots that can highlight issues.

            ### ~ [!] Extract debris ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# If the BestFlank3 is downstream and the GapFlank3 has no HiC support, and the BestFlank5 gap also has
            #i# no HiC support, then remove the entire section between GapFlank5 and BestFlank5 and move to "debris", i.e.
            #i# extract into its own sequence.

            ### ~ [!] Relocate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Look for a chunk between flanks where the ends (1) do not have BestFlank support, and (2) are the
            #i# two BestFlank entries for the same gap -> lift out that chunk and insert into the gap. (May need to invert.)

            ### ~ [!] Translocate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Look for a pair of gaps without BestFlank support
            # reciprocal swap where the BestFlank5 of Gap 1 is the GapFlank5 of Gap2 and vice versa.
            #i# Then check that the BestFlank3 for both is also not the current GapFlank3. (Would this indicate a duplication?)
            #i# If the BestFlank is not GapFlank5, identify the BestFlank5 for the current GapFlank3 of the gap that
            #i# has the focal BestFlank5 as its GapFlank5. Then scan downstream for


            #i# ~ [!] Breakpoints ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# When a Brk gap is encountered, look for something that would fit in the gap. If something is found,
            #i# and has Tran/Brk flanks, replace them with 'Mis' ("Misassembly").


            ### ~ [!] Translocations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Find regions that are flanked by two Trans gaps with the same Hit. Then find another Trans gap that is
            #i# adjacent and close enough to the flanked region to move.

            #!# Will want to differentiate edits between transpositions and single break re-join swaps.

            #!# Identify a gap that has a different BestJoin to the current GapFlank
            # -> Identify the current contig that is on that side
            # -> Identify the better contig to be on that side => Check the score for reciprocal improvement
            # -> Q. Should we check the other end of each contig? I think not - let's hope the single rearrangements solve it!


            #i# ~ [!] Dodgy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Identify dodgy regions that have non-syntenic gaps and duplications etc. without a clear fix.
            #!# After other fixes are applies, identify the Dup regions and output to a file. (Maybe with the gaps
            #!# themselves?) Then run rje_lrbridge.regionSpan(self,regfile). If no edits made, can just use the BAM
            #!# file already made of mapped reads. Otherwise, will need to re-map read data to corrected output.
            #!# -> Use the assembly map to identify the regions?

            return True
        except: self.errorLog('%s.tidyTables error' % self.prog()); return False
#########################################################################################################################
    def synBadHiCBestTranslocations(self,qh): ### Updates based on best HiC support (BestFlank matching GapFlank)
        '''
        Updates based on best HiC support (BestFlank matching GapFlank)
        # 1. Generate a list of non-pure gaps as (flank1,flank2) (sort)
        # 2. Generate a list of all flanks of bad gaps.
        # 3. Generate a list of pairs of flanks at the end of blocks -> assess via list [2] for misplaced blocks.
        # 4. Generate a list of "best" pairs for each misplaced block = (up,down) (sort)
        # 5. Generate a list of actual pairs for each misplaced block = (up,down) (sort)
        # 6. Cross-reference list [5] with hicbest pairs -> "extract" and close gap if found
        # 7. Cross-reference list [4] with list [1] -> "relocate" if a block goes somewhere else (NOTE: Could be an extracted block)
        # 8. Add actual pairs of relocated block to list [1]. (Should have already been caught by [6] if a best pair.)
        # 9. Cycle until no more extract/relocate corrections made.
        #10. Identify and report pairs of bad gaps that share hicbest flanks? (Could add as swap corrections and have toggle to execute?)
        >> qh:str = Qry or Hit
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not self.dev(): return True
            self.headLog('Assessing SynBad {0} HiC translocations'.format(qh),line='~')
            locdb = self.db(qh.lower())
            gapdb = self.dbTable(qh,'gap')
            gapdb.index('GapFlank5')
            gapdb.index('GapFlank3')
            bdb = self.dbTable(qh,'blocks')
            cdb = self.dbTable(qh,'corrections')
            qseqname = qh
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            hicdb = self.addTable(qh,'hicbest',expect=False,make=False)
            if not hicdb or not hicdb.entryNum():
                self.printLog('#HIC','No HiC best pairs data: no HiC translocations')
                return False
            #i# Behaviour will depend on gap types
            # goodgaps # Avoid translocations
            # puregaps # Avoid any disruption
            ## ~ [0a] Set up list of best pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bestset = ['Both']   # Could consider: ['Flank1','Flank2','Weak','Weak1','Weak2']
            bestpairs = []
            for entry in hicdb.entries():
                if entry['Type'] == 'Check' and entry['Best'] in bestset: bestpairs.append([entry['Flank1'],entry['Flank2']])

            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Updates based on best HiC support (BestFlank matching GapFlank)
            # 1. Generate a list of non-pure gaps as (flank1,flank2) (sort)
            # 2. Generate a list of all flanks of bad gaps.
            badgaps = []
            badflanks = []
            termflanks = [] # List of all terminal block flanks
            for entry in gapdb.entries():
                if entry['SynBad'] not in puregaps:
                    flanks = [entry['GapFlank5'],entry['GapFlank3']]
                    flanks.sort()
                    badgaps.append(flanks)
                    badflanks += flanks
            ctgdb = self.dbTable(qh,'contigs')
            prev = None
            for contig in ctgdb.entries(sorted=True):
                seqname = contig['SeqName']
                if not prev or prev['SeqName'] != seqname:
                    if prev: badflanks.append(prev['Flank3'])
                    badflanks.append(contig['Flank5'])
                prev = contig
            if prev: badflanks.append(prev['Flank3'])
            # 3. Generate a list of pairs of flanks at the end of blocks -> assess via list [2] for misplaced blocks.
            badblocks = []
            for entry in bdb.entries():
                if entry['SynType'] == 'Gap': continue  #i# Do not want to assess Gaps!
                if entry['Flank5'] in badflanks and entry['Flank3'] in badflanks:
                    badblocks.append(entry)
            # 4. Generate a list of "best" pairs for each misplaced block = (up,down) (sort)
            besthic = self.bestHiCPairDict(qh)
            self.printLog('#BEST','Best HiC pairing for {0} {1} flanks'.format(rje.iLen(besthic),qh))
            blockbest = []
            for block in badblocks:
                if cdb.data((block[qseqname],block[qstart],block[qend],'invert')):
                    #i# This block is already being relocated by inversion
                    continue
                flanks = []
                if block['Flank5'] in besthic: flanks.append(besthic[block['Flank5']])
                else: flanks.append('NA')
                if block['Flank3'] in besthic: flanks.append(besthic[block['Flank3']])
                else: flanks.append('NA')
                flanks.sort()
                blockbest.append(flanks)
            # 5. Cross-reference list [4] with list [1] -> "relocate" if a block goes somewhere else (NOTE: Could be an extracted block)
                if flanks in badgaps:   # Insert this block into this gap!
                    if flanks[1] == besthic[block['Flank5']]: flanks.reverse()
                    #i# NOTE: The first flank in details matches the 5' of the block. The second matches the 3'.
                    centry = {'SeqName':block[qseqname],'Start':block[qstart],'End':block[qend],
                              'Flank1':block['Flank5'],'Flank2':block['Flank3'],'Edit':'relocate',
                              'Details':'{0}:^:{1}'.format(flanks[0],flanks[1])}
                    cdb.addEntry(centry)
                #!# Look for end of sequence
                elif flanks[0] in termflanks or flanks[1] in termflanks:
                    #!# Add code here to append to end of sequence.
                    pass
            # 6. Generate a list of actual pairs for each misplaced block = (up,down) (sort)
            blockflanks = []
            for block in badblocks:
                if cdb.data((block[qseqname],block[qstart],block[qend],'invert')):
                    #i# This block is already being extracted for inversion
                    continue
                flanks = []
                region1 = block['Flank5']
                region2 = block['Flank3']
                if region1 in gapdb.index('GapFlank3'):
                    flanks.append(gapdb.indexEntries('GapFlank3',region1)[0]['GapFlank5'])
                else: flanks.append('NA')
                if region2 in gapdb.index('GapFlank5'):
                    flanks.append(gapdb.indexEntries('GapFlank5',region2)[0]['GapFlank3'])
                else: flanks.append('NA')
                if 'NA' in flanks: continue
                flanks.sort()
                blockflanks.append(flanks)
            # 7. Cross-reference list [6] with hicbest pairs -> "extract" and close gap if found
                if flanks in bestpairs:
                    centry = {'SeqName':block[qseqname],'Start':block[qstart],'End':block[qend],
                              'Flank1':block['Flank5'],'Flank2':block['Flank3'],'Edit':'extract',
                              'Details':'Flanking gaps are mutual HiC best pairing'}
                    cdb.addEntry(centry)
            #!# Future extension...
            # 8. Add actual pairs of relocated blocks to list [1]. (Should have already been caught by [6] if a best pair.)
            # 9. Cycle until no more extract/relocate corrections made.

            #!# Add capacity to move a block to the end of a scaffold

            #10. Identify and report pairs of bad gaps that share hicbest flanks? (Could add as swap corrections and have toggle to execute?)
            #i# Swap edits are skipped if either flank is part of an inversion, extraction or relocation
            cdb.index('Edit',force=True)
            fixflanks = []
            for etype in ['invert','relocate','extract']:
                fixflanks += cdb.indexDataList('Edit',etype,'Flank1')
                fixflanks += cdb.indexDataList('Edit',etype,'Flank2')
            for entry in gapdb.entries():
                if entry['SynBad'] in puregaps: continue
                if 'Fix' in entry['SynBad']: continue
                if entry['BestFlank5'] == entry['GapFlank5']: continue
                if entry['BestFlank3'] == entry['GapFlank3']: continue
                if entry['BestFlank5'] not in besthic: continue
                if entry['BestFlank3'] not in besthic: continue
                if entry['BestFlank5'] in fixflanks or entry['BestFlank3'] in fixflanks: continue
                if entry['GapFlank5'] in fixflanks or entry['GapFlank3'] in fixflanks: continue
                #i# Going to break this gap if either BestFlank is itself in badflanks;
                #i# Will also need to break that gap to re-join, but this should be done (gap not in puregaps, therefore flanks in badflanks
                breakgap = False
                if entry['BestFlank5'] in badflanks and besthic[entry['BestFlank5']] == entry['GapFlank3']:
                    centry = {'SeqName':entry['SeqName'],'Start':entry['Start'],'End':entry['End'],
                              'Flank1':entry['BestFlank5'],'Flank2':entry['GapFlank3'],'Edit':'join',
                              'Details':'Join two gaps flanks that are mutual HiC best pairs'}
                    cdb.addEntry(centry,overwrite=False)
                    breakgap = True
                if entry['BestFlank3'] in badflanks and besthic[entry['BestFlank3']] == entry['GapFlank5']:
                    centry = {'SeqName':entry['SeqName'],'Start':entry['Start'],'End':entry['End'],
                              'Flank1':entry['GapFlank5'],'Flank2':entry['BestFlank3'],'Edit':'join',
                              'Details':'Join two gaps flanks that are mutual HiC best pairs'}
                    cdb.addEntry(centry,overwrite=False)
                    breakgap = True
                if breakgap:
                    centry = {'SeqName':entry['SeqName'],'Start':entry['Start'],'End':entry['End'],
                              'Flank1':entry['GapFlank5'],'Flank2':entry['GapFlank3'],'Edit':'break',
                              'Details':'Two gaps have flanks that are mutual HiC best pairing'}
                    cdb.addEntry(centry,overwrite=False)

            return True
        except: self.errorLog('%s.synBadHiCBestTranslocations error' % self.prog()); return False
#########################################################################################################################
    def updateHitGHaps(self,qh): ### Uses alt gapdb gap ratings to update qh locdb.
        '''
        Uses alt gapdb gap ratings to update qh locdb.
        >> qh:str = Qry or Hit
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Updating SynBad {0} hit gaps'.format(qh),line='~')
            locdb = self.db(qh.lower())
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            gapdb = self.dbTable(qh,'gap')
            altgapdb = self.dbTable(alt,'gap')
            for field in ['Hit3','Hit5','Qry3','Qry5']:
                if field not in locdb.fields(): locdb.addField(field,after='HitGap',evalue='.')
            q5 = '{0}5'.format(qh)
            q3 = '{0}3'.format(qh)
            qgap = '{0}Gap'.format(qh)
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            ### ~ [1] Sort QryGaps and build hitsort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitdict = {}    # {(hitseq,hitpos,start/end/gap}:entry}
            #!# This is replacing hitsort = (hitseq,hitpos,start/end/gap,entry)
            #!# Replace with {(tuple):entry} and [(tuple)] list to sort - can no longer sort dictionaries.
            #!# Problem here is the mixture of entry and gentry['SynBad'] else could just replace entry with ekey
            #x# ekey = locdb.makeKey(entry)
            entries = locdb.entries(sorted=True)
            etot = len(entries) * 4 + gapdb.entryNum() * 3
            ex = 100.0
            entry = entries.pop(0)
            hitdict[(entry[alt],entry[hstart],'Start')] = entry
            hitdict[(entry[alt],entry[hend],'End')] = entry
            while entries:
                self.progLog('\r#SYNGAP','Updating {0} {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                prev = entry
                entry = entries.pop(0)
                entry[hgap] = ''
                entry[h5] = entry[h3] = '-'
                if self.isGap(entry[qgap]): continue
                entry[qgap] = ''
                if self.isGap(prev[qgap]):
                    entry[qgap] += '<'; entry[q5] = prev[qgap]
                if entries and self.isGap(entries[0][qgap]):
                    entry[qgap] += '>'; entry[q3] = entries[0][qgap]
                if not entry[qgap]: entry[qgap] = '.'
                hitdict[(entry[alt],entry[hstart],'Start')] = entry
                hitdict[(entry[alt],entry[hend],'End')] = entry
            ## ~ [1a] Add hit gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gentry in altgapdb.entries():
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                hitdict[(gentry['SeqName'],gentry['Start'],'gapstart')] = gentry['SynBad']
                hitdict[(gentry['SeqName'],gentry['End'],'gapend')] = gentry['SynBad']
            hitsort = list(hitdict.keys())
            hitsort.sort()
            self.bugPrint(hitsort[:30])
            self.debug(hitsort[-30:])
            ### ~ [2] Update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(hitsort)):
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                #i# Find gap and look each side
                if hitsort[i][2] == 'gapstart':
                    if i == 0: continue
                    j = i-1
                    if hitsort[j][2].startswith('gap'): continue
                    jend = hitsort[j][1]
                    while j >= 0 and hitsort[j][1] == jend:
                        if hitsort[j][2] == 'End':
                            hitdict[hitsort[j]][hgap] += '>'
                            hitdict[hitsort[j]][h3] = hitdict[hitsort[i]]
                        j -= 1
                if hitsort[i][2] == 'gapend':
                    j = i+1
                    if j >= len(hitsort): continue
                    if hitsort[j][2].startswith('gap'): continue
                    jstart = hitsort[j][1]
                    while j < len(hitsort) and hitsort[j][1] == jstart:
                        if hitsort[j][2] == 'Start':
                            hitdict[hitsort[j]][hgap] += '<'
                            hitdict[hitsort[j]][h5] = hitdict[hitsort[i]]
                        j += 1
            for entry in locdb.entries():
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                if not entry[hgap]: entry[hgap] = '.'
                if len(entry[hgap]) > 1 and entry[hgap][:1] == '.': entry[hgap] = entry[hgap][1:]
                if len(entry[qgap]) > 1 and entry[qgap][:1] == '.': entry[qgap] = entry[qgap][1:]
                if entry[hgap] == '><': entry[hgap] = '<>'
            self.printLog('\r#UPDATE','Updating {0} flanking {1}, {2} and {3} complete.'.format(qh,hgap,h5,h3))

            return True
        except: self.errorLog('%s.updateHitGHaps error' % self.prog()); return False
#########################################################################################################################
    def isGap(self,gapstr): ### Returns whether gapstr is a gap entry.
        if type(gapstr) == int: return False
        if ('>' in gapstr or '<' in gapstr or gapstr == '.'): return False
        return True
#########################################################################################################################
    def synBadSyntenyBlocks(self,qh,quick=False,skiptypes=None):   ### SynBad collinear chunks and Cis/Trans assessment
        '''
        SynBad collinear chunks and Cis/Trans assessment.
        #># Generate a new fragments table that skips through Syn, Aln, Span, Long and Div gaps. Sort by size?
        # and look for out of place neighbours. Build a collinear backbone in size order? Base this on the ends.
        # Identity all that don't fit. Sub classify as cis and trans. Look at whether cis can be moved into any gaps.
        # Need to identify cis and trans gaps too - include the Div gaps as cis, ignoring the non-cis region.
        # This latter bit is basically looking at the pairs tables but including the gaps in every pair?
        # Need to rename the Div gaps so they have the cis scaffolds in their names.
        # Can then add cis gaps to the pairs table and work off that.

        # This method also performs the minctglen and minbadctg sequence extraction corrections.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        >> skiptypes:list [] = List of gap types to skip when making chunks. (None=Default)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if skiptypes == None:
                skiptypes = skipgaps
            locdb = self.db(qh.lower())
            gapdb = self.dbTable(qh,'gap')
            cdb = self.addTable(qh,'corrections',make=True,replace=False)
            #if 'SynType' not in gapdb.fields(): gapdb.addField('SynType')
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            ## ~ [0a] Generate synteny blocks table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            syntable = '{0}.blocks'.format(qh.lower())
            syndb = self.db().copyTable(locdb,syntable,replace=True,add=True)
            #syndb.addField('GapNum',evalue=0)
            syndb.addField('SynType',evalue="")
            syndb.dropFields(['Alt','Trans','Inv','Dup',hgap,'Qry5','Qry3','Hit5','Hit3'])
            entries = syndb.entries(sorted=True)

            ### ~ [1] Generate synteny blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while entries:
                #i# Start new block
                block = entries.pop(0)
                #i# Check for intervening gap to skip
                if self.isGap(block[qgap]):
                    block['SynType'] = 'Gap'
                    continue
                block[qgap] = 0
                #i# Find end of block
                blockend = block
                prev = block
                while entries and (entries[0][qh] == block[qh]) and ('>' in entries[0][qgap] or '<' in entries[0][qgap] or entries[0][qgap] == '.' or entries[0][qgap] in skiptypes):
                    if entries[0][qgap] in skiptypes: block[qgap] += 1
                    blockend = entries.pop(0)
                    #i# Will probably want to fix or drop these at some point
                    for field in ['Length','Identity','Non']:
                        block[field] += blockend[field]
                    block['Non'] += max(0,blockend[qstart]-prev[qend]-0)
                    prev = blockend
                    syndb.dropEntry(blockend)
                #i# Update entry
                if blockend:
                    block[qend] = blockend[qend]
                    bstart = block[hstart]
                    bend = blockend[hend]
                    if block['Strand'] == '-': bstart = block[hend]
                    if blockend['Strand'] == '-': bend = blockend[hstart]
                    block[hstart] = bstart
                    block[hend] = bend
                    block['Strand'] = block['Strand'] + blockend['Strand']
                    if block[alt] != blockend[alt]:
                        block['SynType'] = 'Trans'
                        block[alt] = '{0}::{1}'.format(block[alt],blockend[alt])
                    else: block['SynType'] = 'Cis'
            #syndb.dropField(qgap)
            syndb.indexReport('SynType')

            ### ~ [2] Add flanks to Blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            syndb.addField('Flank5')
            syndb.addField('Flank3')
            ctgdb = self.dbTable(qh,'contigs')
            termflanks = {}
            for contig in ctgdb.entries(sorted=True):
                seqname = contig['SeqName']
                if seqname not in termflanks:
                    termflanks[seqname] = [contig['Flank5'],contig['Flank3']]
                else:
                    termflanks[seqname][1] = contig['Flank3']
            entries = syndb.entries(sorted=True)
            prev = None
            while entries:
                entry = entries.pop(0)
                #i# New sequence
                try:
                    if not prev or entry[qh] != prev[qh]:
                        #!# Updated code for missing flanks, possibly because of terminal Ns
                        #!# This is a temporary fix -> Need to get to bottom of it
                        if prev:
                            try: prev['Flank3'] = termflanks[prev[qh]][1]
                            except: prev['Flank3'] = ''
                        try: entry['Flank5'] = termflanks[entry[qh]][0]
                        except: entry['Flank5'] = ''
                    elif entry['SynType'] == 'Gap':
                        gap = gapdb.data(syndb.makeKey(entry))
                        prev['Flank3'] = gap['GapFlank5']
                        entry['Flank5'] = gap['GapFlank5']
                        entry['Flank3'] = gap['GapFlank3']
                    else:
                        entry['Flank5'] = prev['Flank3']
                except:
                    self.errorLog('Problem with blocks entry: {0}'.format(syndb.entrySummary(entry,collapse=True)))
                prev = entry
            if prev: prev['Flank3'] = termflanks[prev[qh]][1]

            ### ~ [3] Check for short contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flanks = {}     # 5' : 3' contig flanks
            for contig in ctgdb.entries():
                flanks[contig['Flank5']] = contig['Flank3']
                if contig['CtgLen'] < self.getInt('MinCtgLen'):
                    centry = {'SeqName':contig['SeqName'],'Start':contig['Start'],'End':contig['End'],
                              'Flank1':contig['Flank5'],'Flank2':contig['Flank3'],'Edit':'remove',
                              'Details':'Contig failed to meet MinCtgLen.'}
                    cdb.addEntry(centry,overwrite=False)
            for block in syndb.entries():
                if block['SynType'] == 'Gap': continue
                if block['Flank5'] not in flanks:
                    if block['Flank5']:
                        self.warnLog('Cannot find block 5\' region ({0}) in contig flanks'.format(block['Flank5']))
                    continue
                if block['Flank3'] == flanks[block['Flank5']]:  # contig!
                    if (block[qend] - block[qstart] + 1) < self.getInt('MinBadCtg'):
                        centry = {'SeqName':block[qh],'Start':block[qstart],'End':block[qend],
                              'Flank1':block['Flank5'],'Flank2':block['Flank3'],'Edit':'extract',
                              'Details':'Single-contig block failed to meet MinBadCtg.'}
                        cdb.addEntry(centry,overwrite=False)


            syndb.saveToFile()

            return True
        except: self.errorLog('%s.synBadChunkTable error' % self.prog()); return False
#########################################################################################################################
    def synBadDivergence(self,qh,locdb,gapdb,maxskip=0,quick=False,gaptypes=('Brk','Tran')):   ### Main SynBad gap mapping compression for a local hits table.
        '''
        Main SynBad gap mapping compression for a local hits table.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: quick = True
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = [('NSCUCHR1.01',95640,477499),('NSCUCHR1.01',1312891,1314913),
                       ('NSCUCHR1.01',1320580,1323528),
                       ('NSCUCHR1.01',1316193,1319918)]
            debugme = []

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            while entries:
                self.progLog('\r#DIV','Assessing %s table divergence: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = False # self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.': continue #i# Don't merge gaps
                #i# - (2) Currently only want one internal gap of chosen type
                skipgap = False
                gapnum = 0
                trangap = None
                for entry in toskip:
                    if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if entry[qgap] not in gaptypes: skipgap = True
                    else: gapnum += 1; trangap = entry
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                elif gapnum != 1:
                    if debug: self.bugPrint('NoFocalGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                maxspan = False
                if qdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                    maxspan = True
                if hdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                    maxspan = True
                #i# maxspan=True when maximum length breached: used for converting Brk to Long
                if maxspan:
                    if trangap[qgap] != 'Brk' or eskip != 1: continue

                ## ~ [1c] Re-classify gap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                try:
                    gkey = locdb.makeKey(trangap)
                    newgap = 'Div'
                    if maxspan: newgap = 'Long'
                    self.printLog('\r#DIV','Reclassifying {0} {1}-{2} {3} {4} -> {5}'.format(gkey[0],gkey[1],gkey[2],trangap[alt],trangap[qgap],newgap))
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = newgap
                    trangap[qgap] = newgap
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(trangap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1
                except:
                    self.errorLog('%s.synBadDivergence error' % self.prog())
                    self.warnLog('Failed to update gap class')
            self.printLog('\r#DIV','Assessed %s table divergence: %s gaps reclassified' % (tname,rje.iStr(dx)))

            return dx
        except: self.errorLog('%s.synBadDivergence error' % self.prog()); raise
#########################################################################################################################
    def synBadInversions(self,qh,locdb,gapdb,maxskip=0,quick=False):   ### SynBad inversion correction.
        '''
        SynBad inversion correction.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            tdb = self.dbTable(qh,'corrections')
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: maxskip = etot
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = []#('NSCUCHR1.01',11501373,11553304),('NSCUCHR1.01',58321100,58755424)]

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            focus = None
            while entries:
                if focus and locdb.makeKey(focus) in debugme: self.deBug('?')
                self.progLog('\r#DIV','Assessing simple %s table inversions: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\n\nSkip={0}; Lookback={1}; Focus={2}'.format(eskip,lookback,locdb.makeKey(focus)))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                # if debug:
                #     for entry in [merge1] + toskip:
                #         self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                #     self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (0) Inversions
                if len(toskip) < 3: continue
                if toskip[0][qgap] != 'Inv':
                    eskip = maxskip
                    continue
                if toskip[-1][qgap] != 'Inv': continue
                if toskip[1][alt] != merge1[alt]:
                    eskip = maxskip
                    continue
                if toskip[-2][alt] != merge1[alt]:
                    continue
                eskip = maxskip
                #if not debug: self.bugPrint('\n\nSkip={0}; Lookback={1}; Focus={2}'.format(eskip,lookback,locdb.makeKey(focus)))
                #debug = self.debugging()
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.bugPrint('%s' % (locdb.entrySummary(merge2,collapse=True)))
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.':
                    continue #i# Don't merge gaps
                #i# - (2) Will allow skipping of "good" gaps, including Dup gaps.
                skipgap = False
                trangap = None
                for entry in toskip[1:-1]:
                    #X#if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if self.isGap(entry[qgap]) and entry[qgap] not in skipgaps:
                        skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                #if qdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                #    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                #if hdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                #    continue
                #i# ~ (5) Positions of inverted region
                (flankstart, flankend) = (merge1[hend],merge2[hstart])
                (invstart, invend) = (toskip[-2][hstart],toskip[1][hend])
                if bwd:
                    (flankstart, flankend) = (merge2[hstart],merge1[hend])
                    (invstart, invend) = (toskip[1][hstart],toskip[-2][hend])
                #self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))
                if invstart - flankstart + self.getInt('MaxOverlap') < 0:
                    if debug: self.bugPrint('ToSkip start too small')
                    continue
                if flankend - invend + self.getInt('MaxOverlap') < 0:
                    if debug: self.bugPrint('ToSkip end too big')
                    continue

                ## ~ [1c] Invert sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                for invgap in (toskip[0],toskip[-1]):
                    gkey = locdb.makeKey(invgap)
                    self.printLog('\r#INV','Reclassifying {0} {1}-{2} {3} {4} -> InvFix'.format(gkey[0],gkey[1],gkey[2],invgap[alt],invgap[qgap]))
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = 'InvFix'
                    invgap[qgap] = 'InvFix'
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(invgap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1
                #i# Fix intervening sequences: invert and reverse order
                # i1 = merge1 -> +1 = Inv1
                # i2 = merge2 -> -1 = Inv2
                #!# No longer doing this -> might cause issues. Re-run on the edited fasta files instead.
                invstart = toskip[0][qend] + 1
                invend = toskip[-1][qstart] - 1
                #i# Make sure that this inversion is documented -> will need to check/execute later (or reverse)
                self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                # tdb.addEntry({'seqname':merge1[qh],'Start':invstart,'End':invend,'edit':'invert','details':'Inversion in-place.'})
                flank1 = gapdb.data((toskip[0][qh],toskip[0][qstart],toskip[0][qend]))['GapFlank3']
                flank2 = gapdb.data((toskip[-1][qh],toskip[-1][qstart],toskip[-1][qend]))['GapFlank5']
                tdb.addEntry({'SeqName':merge1[qh],'Start':invstart,'End':invend,'Flank1':flank1,'Flank2':flank2,'Edit':'invert','Details':'Inversion in-place.'})
            locdb.remakeKeys()
            self.printLog('\r#INV','Assessed simple %s table inversions: %s gaps reclassified' % (tname,rje.iStr(dx)))
            if dx: self.warnLog('NOTE: Fixed inversions will not carry over into %s table. (Fix assemblies and re-run SynBad.)' % alt)
            return dx
        except: self.errorLog('%s.synBadInversions error' % self.prog()); raise
#########################################################################################################################
    def synBadInvertedDuplications(self,qh,locdb,gapdb,maxskip=0,quick=False):   ### Extended SynBad inversion correction.
        '''
        Extended SynBad inversion correction.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            tdb = self.dbTable(qh,'corrections')
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: maxskip = etot
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = []#('NSCUCHR1.01',11501373,11553304),('NSCUCHR1.01',58321100,58755424)]

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            while entries:
                self.progLog('\r#DIV','Assessing complex %s table inversions: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                # if debug:
                #     for entry in [merge1] + toskip:
                #         self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                #     self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (0) Inversions
                if len(toskip) < 3: continue
                if toskip[0][qgap] != 'Inv':
                    eskip = maxskip
                    continue
                if toskip[-1][qgap] != 'Inv': continue
                if toskip[1][alt] != merge1[alt]:
                    eskip = maxskip
                    continue
                if toskip[-2][alt] != merge1[alt]:
                    continue
                eskip = maxskip
                debug = self.debugging()
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.bugPrint('%s' % (locdb.entrySummary(merge2,collapse=True)))
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.':
                    continue #i# Don't merge gaps
                #i# - (2) Will allow skipping of "good" gaps, including Dup gaps.
                skipgap = False
                trangap = None
                for entry in toskip[1:-1]:
                    #X#if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if self.isGap(entry[qgap]) and entry[qgap] not in skipgaps:
                        skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                #if qdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                #    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                #if hdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                #    continue
                #i# ~ (5) Positions of inverted region
                (flankstart, flankend) = (merge1[hend],merge2[hstart])
                (invstart, invend) = (toskip[-2][hstart],toskip[1][hend])
                if bwd:
                    (flankstart, flankend) = (merge2[hstart],merge1[hend])
                    (invstart, invend) = (toskip[1][hstart],toskip[-2][hend])
                #self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))

                ## ~ [1c] Invert sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                invert = False
                #i# Next, work through, classify and merge each intermediate entry and merge2
                if invstart - flankstart + self.getInt('MaxOverlap') >= 0:
                    toskip[0][qgap] = 'InvFix'; invert = True
                elif invstart < flankstart < invend:
                    toskip[0][qgap] = 'InvDupFix'; invert = True
                else:
                    toskip[0][qgap] = 'InvBrk'
                if flankend - invend + self.getInt('MaxOverlap') >= 0:
                    toskip[-1][qgap] = 'InvFix'; invert = True
                elif invstart < flankend < invend:
                    toskip[-1][qgap] = 'InvDupFix'; invert = True
                else:
                    toskip[-1][qgap] = 'InvBrk'

                for invgap in (toskip[0],toskip[-1]):
                    gkey = locdb.makeKey(invgap)
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = invgap[qgap]
                    self.printLog('\r#INV','Reclassifying {0} {1}-{2} {3} Inv -> {4}'.format(gkey[0],gkey[1],gkey[2],invgap[alt],invgap[qgap]))
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(invgap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1


                #i# Fix intervening sequences: invert and reverse order
                # i1 = merge1 -> +1 = Inv1
                # i2 = merge2 -> -1 = Inv2
                if invert:
                    invstart = toskip[0][qend] + 1
                    invend = toskip[-1][qstart] - 1
                    #i# Make sure that this inversion is documented -> will need to check/execute later (or reverse)
                    self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                    # tdb.addEntry({'seqname':merge1[qh],'Start':invstart,'End':invend,'edit':'invert','details':'Inversion in-place.'})
                    flank1 = gapdb.data((toskip[0][qh],toskip[0][qstart],toskip[0][qend]))['GapFlank3']
                    flank2 = gapdb.data((toskip[-1][qh],toskip[-1][qstart],toskip[-1][qend]))['GapFlank5']
                    tdb.addEntry({'SeqName':merge1[qh],'Start':invstart,'End':invend,'Flank1':flank1,'Flank2':flank2,'Edit':'invert','Details':'Inversion in-place.'})
                    # tdb.addEntry({'seqname':merge1[qh],'flank1':invstart,'flank2':invend,'edit':'invert','details':'Inversion in-place.'})
            locdb.remakeKeys()
            self.printLog('\r#INV','Assessed complex %s table inversions: %s gaps reclassified' % (tname,rje.iStr(dx)))
            if dx: self.warnLog('NOTE: Fixed inversions will not carry over into %s table. (Fix assemblies and re-run SynBad.)' % alt)

            return dx
        except: self.errorLog('%s.synBadInvertedDuplications error' % self.prog()); raise
#########################################################################################################################
    ### <11> ### SynBad Correction Methods                                                                              #
#########################################################################################################################
    def couldConflict(self,qh,entry): # Whether a correction could conflict with another
        '''
        Whether a correction could conflict with another.
        >> qh:str = qry/hit
        >> correction: dict = Corrections table entry #['seqname','Start','End','Flank1','Flank2','edit']
        :return: True/False whether this edit could cause conflicts
        '''
        #i# Return False if edit is a single contig
        if entry['Flank1'] == entry['Flank2']: return False
        if entry['Edit'] in ['break','join','swap']: return False
        flank1 = entry['Flank1']
        flank2 = entry['Flank2']
        amap = self.list[qh]
        try:
            mapi = [amap.index(flank1), amap.index(flank2)]
        except:
            self.errorLog('Problem mapping flank during conflict check')
            return True
        mapi.sort()
        editchunk = amap[mapi[0]:mapi[1]+1]
        if len(editchunk) == 5: return False
        #i# ELse, return True
        return True
#########################################################################################################################
    def checkSubMap(self,submap,gap5='',gap3=''):   ### Checks whether a submap list is a legitimate chunk
        '''
        Checks whether a submap list is a legitimate chunk.
        >> submap:list = list of assembly map elements
        >> gap5:str = optional upstream gap position to check
        >> gap3:str = optional downstream gap position to check
        '''
        if gap5 and gap5[:1] not in ['|',':']:
            self.bugPrint(' '.join(submap))
            self.debug(gap5)
            raise ValueError('Upstream element not an end or gap')
        if gap3 and gap3[:1] not in ['|',':']:
            self.bugPrint(' '.join(submap))
            self.debug(gap3)
            raise ValueError('Downstream element not an end or gap')
        if len(submap) % 6 != 5:
            self.debug(' '.join(submap))
            raise ValueError('Length of assembly submap inconsistent with whole contigs')
        if '|' in submap:
            self.debug(' '.join(submap))
            raise ValueError('Assembly submap spans multiple scaffolds.')
        return True
#########################################################################################################################
    def mapCorrections(self):  ### Loads corrections table, edits sequences and outputs corrected assembly
        '''
        Loads corrections table, edits sequences and outputs corrected assembly.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Check for existing output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            wanted = ['qry.synbad.txt','hit.synbad.txt','qry.synbad.fasta','hit.synbad.fasta']
            if not self.update() and rje.checkForFiles(wanted,basename=self.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Corrected SynBad assembly files found: skipping error-correction (force=F update=F)')
                for qh in ['qry','hit']:
                    self.loadAssemblyMap(qh,'synbad')
                return True
            if not self.force() and not self.dev():
                self.printLog('#FORCE','Setting force=T: downstream steps will not be skipped.')
                self.setBool({'Force':True})

            self.seqObjSetup()
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            goodedits = self.list['Correct'] + ['remove'] #['invert','extract','relocate']
            for qh in ['qry','hit']:
                self.addTable(qh,'corrections')
                # if cdb: db.deleteTable(cdb)
                # cdb = db.addTable(mainkeys=['SeqName','Start','End','Edit'],name='{0}.corrections'.format(qh),replace=True,ignore=[],expect=True)

            ### ~ [1] Check edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Check for conflicting edits, i.e. those sharing a flank? Might want to have an explicit ordering added.
            #!# Add start and end positions to table for correct sorting. Use this for overlaps? (See older code.)
            #!# Make sure it is reported if edits are skipped. (Filter edit types first and add commandline option.)
            for qh in ['qry','hit']:
                fixed = []  # List of entries already fixed in previous run
                cdb = self.dbTable(qh,'corrections')
                if 't' in self.list['Correct'] or 'true' in self.list['Correct']:
                    goodedits += cdb.indexKeys('Edit')
                prev = None
                for entry in cdb.entries(sorted=True):
                    if entry['SynBad'] == 'Fixed':
                        fixed.append(entry)
                        continue
                    if entry['Edit'] not in goodedits:
                        entry['SynBad'] = 'Not implemented'
                        continue
                    if not self.couldConflict(qh,entry): continue
                    if prev and prev['End'] >= entry['Start'] and prev['SeqName'] == entry['SeqName']:
                        self.warnLog('Dropped correction {0} {1}-{2} due to overlapping edit region.'.format(entry['SeqName'],entry['Start'],entry['End']))
                        entry['SynBad'] = 'Blocked'
                        continue
                    prev = entry
                cdb.indexReport('SeqName')
                cdb.indexReport('Edit')
                for etype in cdb.index('Edit'):
                    if etype not in goodedits:
                        self.warnLog('Edit type "{0}" not in correct=LIST.'.format(etype))

            ### ~ [2] Update assembly map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                etype = 'invert'
                if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                    self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                    continue
                cdb = self.dbTable(qh,'corrections')
                totx = 0    # Total number of edits
                ## ~ [2a] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} inversions.'.format(qh),line='~')
                self.infoLog('Inversions invert the orientation between two flanks.')
                self.progLog('\r#INV','Processing {0} inversions'.format(qh))
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','invert'):
                    if entry in fixed: continue
                    if entry['SynBad'] in ['Blocked']: continue
                    if self.mapInversion(qh,entry['Flank1'],entry['Flank2']):
                        entry['SynBad'] = 'Fixed'
                        fixed.append(entry)
                        editx += 1
                    else:
                        entry['SynBad'] = 'Failed'
                        failx += 1
                totx += editx
                self.printLog('\r#INV','Processed {0} inversions: {1} edits; {2} failed'.format(qh,editx,failx))

                ## ~ [2b] Extractions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} extractions.'.format(qh),line='~')
                self.infoLog('Extractions pull a block out of a pair of gaps and join the gaps back together.')
                    # centry = {'SeqName':block[qseqname],'Start':block[qstart],'End':block[qend],
                    #           'Flank1':block['Flank5'],'Flank2':block['Flank3'],'Edit':'extract',
                    #           'Details':'Flanking gaps are mutual HiC best pairing'}
                self.progLog('\r#EXT','Processing {0} extractions'.format(qh))
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','extract'):
                    etype = 'extract'
                    if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                        self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                        break
                    if entry in fixed: continue
                    if entry['SynBad'] in ['Blocked']: continue
                    if self.mapExtraction(qh,entry['Flank1'],entry['Flank2']):
                        entry['SynBad'] = 'Fixed'
                        fixed.append(entry)
                        editx += 1
                    else:
                        entry['SynBad'] = 'Failed'
                        failx += 1
                totx += editx
                self.printLog('\r#EXT','Processed {0} extractions: {1} edits; {2} failed'.format(qh,editx,failx))

                ## ~ [2c] Relocations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} relocations.'.format(qh),line='~')
                self.infoLog('Relocations will extract a block if required and then insert it into a gap.')
                #i# The first flank in details is where the 5' of the block goes
                #centry = {'SeqName':block[qseqname],'Start':block[qstart],'End':block[qend],
                #              'Flank1':block['Flank5'],'Flank2':block['Flank3'],'Edit':'relocate',
                #              'Details':'{0}:^:{1}'.format(flanks[0],flanks[1])}
                self.progLog('\r#INS','Processing {0} relocations'.format(qh))
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','relocate'):
                    etype = 'relocate'
                    if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                        self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                        break
                    if entry in fixed: continue
                    if entry['SynBad'] in ['Blocked']: continue
                    if self.mapInsertion(qh,entry['Flank1'],entry['Flank2'],entry['Details']):
                        entry['SynBad'] = 'Fixed'
                        fixed.append(entry)
                        editx += 1
                    else:
                        entry['SynBad'] = 'Failed'
                        failx += 1
                totx += editx
                self.printLog('\r#INS','Processed {0} relocations: {1} edits; {2} failed'.format(qh,editx,failx))

                ## ~ [2d] Breaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} breaks and exchanges.'.format(qh),line='~')
                self.infoLog('Breaks prime exchanges of two block flanks that are mutual best hits.')
                self.infoLog('The next phase will look to join terminal flanks that are mutual best hits.')
                self.infoLog('The "break" edit types will have two flanks of a gap.')
                self.infoLog('The "join" edit types will have the two flanks that form a new gap.')
                self.progLog('\r#BREAK','Processing {0} breaks'.format(qh))
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','break'):
                    etype = 'break'
                    if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                        self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                        break
                    if entry in fixed: continue
                    if entry['SynBad'] in ['Blocked']: continue
                    if self.mapBreak(qh,entry['Flank1'],entry['Flank2']): entry['SynBad'] = 'Broken'
                    else: entry['SynBad'] = 'Failed'
                    editx += 1
                self.printLog('\r#BREAK','Processed {0} breaks: {1} edits'.format(qh,editx))
                self.progLog('\r#JOIN','Processing {0} joins'.format(qh))
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','join'):
                    etype = 'join'
                    if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                        self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                        break
                    if entry in fixed: continue
                    if entry['SynBad'] in ['Blocked']: continue
                    joined = self.mapJoin(qh,entry['Flank1'],entry['Flank2'])
                    if joined:
                        entry['SynBad'] = 'Fixed'
                        fixed.append(entry)
                        editx += 1
                    elif joined == False:
                        entry['SynBad'] = 'Failed'
                        failx += 1
                totx += editx
                self.printLog('\r#JOIN','Processed {0} joins: {1} edits; {2} failed'.format(qh,editx,failx))

                ## ~ [2e] Removals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} removals.'.format(qh),line='~')
                self.infoLog('Contigs failing to meet minctglen={0} will be extracted.'.format(self.getInt('MinCtgLen')))
                self.progLog('\r#REM','Processing {0} removals'.format(qh))
                remflanks = []
                failx = editx = 0
                for entry in cdb.indexEntries('Edit','remove'):
                    etype = 'remove'
                    if not 't' in self.list['Correct'] and not 'true' in self.list['Correct'] and not etype in self.list['Correct']:
                        self.printLog('#SKIP','Edit type "{0}" not in correct=LIST.'.format(etype))
                        break
                    if entry in fixed: continue
                    remflanks += [entry['Flank1'],entry['Flank2']]
                    if entry['SynBad'] in ['Blocked']: continue
                    if self.mapExtraction(qh,entry['Flank1'],entry['Flank2']):
                        entry['SynBad'] = 'Fixed'
                        fixed.append(entry)
                        editx += 1
                    else:
                        entry['SynBad'] = 'Failed'
                        failx += 1
                totx += editx
                self.printLog('\r#REM','Processed {0} removals: {1} edits; {2} failed'.format(qh,editx,failx))

                ## ~ [2f] HiC Best End Joins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                bdb = self.dbTable(qh,'hicbest',expect=False)
                self.headLog('{0} end joins.'.format(qh),line='~')
                self.infoLog('Join sequence termini that are mutual HiC best partners.')
                if bdb:
                    self.progLog('\r#HIC','Processing {0} HiC end joins'.format(qh))
                    editx = 0
                    for entry in bdb.indexEntries('Best','Both'):
                        if entry['Flank1'] in remflanks: continue
                        if entry['Flank2'] in remflanks: continue
                        #?# Should we try to work out orientation?
                        joined = self.mapJoin(qh,entry['Flank1'],entry['Flank2'],expect=False)
                        if joined:
                            (seqname,i,j) = rje.matchExp('^(\S+)\.(\d+)-(\d+)',entry['Flank1'])
                            centry = {'SeqName':seqname,'Start':int(i),'End':int(j),
                                      'Flank1':entry['Flank1'],'Flank2':entry['Flank2'],'Edit':'join',
                                      'SynBad':'Fixed',
                                      'Details':'Join two gaps flanks that are mutual HiC best pairs'}
                            cdb.addEntry(centry,overwrite=False)
                            fixed.append(entry)
                            editx += 1
                        #else: failx += 1
                    self.printLog('\r#HIC','Processed {0} HiC end joins: {1} edits'.format(qh,editx))
                    totx += editx
                else: self.printLog('\r#HIC','No {0} HiC best pairs to process.'.format(qh))

                ## ~ [2g] Rejoins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.headLog('{0} rejoins.'.format(qh),line='~')
                self.infoLog('The final step is to rejoin termini to recreate original gaps that were not part of a removal.')
                if self.getBool('Rejoin'):
                    gdb = self.dbTable(qh,'gap')
                    self.infoLog('Re-run on output with fragment=T (or rejoin=F) to fragment these gaps.')
                    self.progLog('\r#JOIN','Processing {0} rejoins'.format(qh))
                    editx = 0
                    for entry in gdb.entries():
                        if entry['GapFlank5'] in remflanks: continue
                        if entry['GapFlank3'] in remflanks: continue
                        joined = self.mapJoin(qh,entry['GapFlank5'],entry['GapFlank3'],expect=False)
                        if joined:
                            centry = {'SeqName':entry['SeqName'],'Start':entry['Start'],'End':entry['End'],
                                      'Flank1':entry['GapFlank5'],'Flank2':entry['GapFlank3'],'Edit':'rejoin',
                                      'SynBad':'Rejoined',
                                      'Details':'Rejoined two gaps flanks ended as sequence termini.'}
                            cdb.addEntry(centry,overwrite=False)
                            editx += 1
                    self.printLog('\r#JOIN','Processed {0} rejoins: {1} edits'.format(qh,editx))
                else: self.printLog('\r#JOIN','No {0} rejoin processing (rejoin=F).'.format(qh))

                ## ~ [2h] Fragment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.printLog('#FIX','{0} {1} edits; total {2} "Fixed" entries'.format(qh,totx,len(fixed)))
                if totx and self.getBool('Fragment'):
                    self.setBool({'Fragment':False})
                    self.printLog('#FRAG','Assembly edits incompatible with fragmentation: setting fragment=F')
                    self.infoLog('Re-run with: genome1={0}.qry.synbad.fasta genome2={0}.hit.synbad.fasta mapflanks1={0}.qry.flanks.fasta mapflanks2={0}.hit.flanks.fasta'.format(self.baseFile()))

            ### ~ [3] Output updated map and fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if not self.saveAssemblyMaps(qh,mapname='synbad',hidegaps=self.list['HideGaps']): return False
                if not self.saveTelociraptorMaps(qh,mapname='synbad'): return False
                cdb.saveToFile(backup=False)

            return True
        except:
            self.errorLog('%s.corrections error' % self.prog())
            return False
#########################################################################################################################
    def reverseAssemblyChunk(self,invchunk):
        invchunk.reverse()
        for i in range(len(invchunk)):
            if invchunk[i] == '>': invchunk[i] = '<'
            elif invchunk[i] == '<': invchunk[i] = '>'
        return invchunk
#########################################################################################################################
    def mapInversion(self,qh,flank1,flank2):    ### Inverts assembly map between flank1 and flank2
        '''
        Inverts assembly map between flank1 and flank2.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh]
            checklen = len(amap)
            ### ~ [1] Invert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapi = [amap.index(flank1), amap.index(flank2)]
            if flank1 == flank2:
                mapi[1] = mapi[1] + 4
                if amap[mapi[1]] != flank2:
                    raise ValueError('Problem finding both ends of contig with identical flanks')
            mapi.sort()
            #!# Add check for correct orientation prior to inversion #!#
            invchunk = amap[mapi[0]:mapi[1]+1]
            try: self.checkSubMap(invchunk,amap[mapi[0]-1],amap[mapi[1]+1])
            except:
                self.printLog('#BADMAP',' '.join(amap[mapi[0]-1:mapi[1]+2]))
                self.errorLog('Inversion failed',quitchoice=False)
                return False
            invchunk.reverse()
            for i in range(len(invchunk)):
                if invchunk[i] == '>': invchunk[i] = '<'
                elif invchunk[i] == '<': invchunk[i] = '>'
            self.list[qh] = amap[:mapi[0]] + invchunk + amap[mapi[1]+1:]
            if len(self.list[qh]) != checklen: raise ValueError('Length of assembly map has changed during inversion')
            self.printLog('#EDIT','Inverted {0}...{1}'.format(flank1,flank2))
            return self.list[qh]
        except:
            self.errorLog('%s.mapInversion error' % self.prog())
            return False
#########################################################################################################################
    def mapExtraction(self,qh,flank1,flank2,keep=True):    ### Deletes assembly map between flank1 and flank2
        '''
        Deletes assembly map between flank1 and flank2.
        >> qh:qry/hit
        >> flank1:str = One end of block to remove
        >> flank2:str = Other end of block to remove
        >> keep:bool = Whether to keep extracted block as a new sequence. ([...,'|',block,'|'])
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh]
            ### ~ [1] Delete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapi = [amap.index(flank1), amap.index(flank2)]
            mapi.sort()
            if flank1 == flank2:
                mapi[1] = mapi[1] + 4
                if amap[mapi[1]] != flank2:
                    raise ValueError('Problem finding both ends of contig with identical flanks')
            try: self.checkSubMap(amap[mapi[0]:mapi[1]+1],amap[mapi[0]-1],amap[mapi[1]+1])
            except:
                self.errorLog('Extraction failed',quitchoice=False)
                self.debug(' '.join(amap[mapi[0]-1:mapi[1]+2]))
                return False
            block = amap[mapi[0]:mapi[1]+1]
            ### ~ [2] Assess and deal with adjacent ends/gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Whole sequence extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if amap[mapi[0]-1] == '|' and amap[mapi[1]+1] == '|':
                if keep:
                    self.warnLog('Trying to extract a complete sequence |{0}...{1}|!'.format(flank1,flank2))
                    return True
                else:
                    #i# Delete block plus both sides
                    self.list[qh] = amap[:mapi[0]-1] + amap[mapi[1]+2:]
            ## ~ [2b] 5' terminal sequence extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif amap[mapi[0]-1] == '|':
                #i# Delete block plus 3' gap
                self.list[qh] = amap[:mapi[0]] + amap[mapi[1]+2:]
            ## ~ [2c] 3' terminal sequence extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif amap[mapi[1]+1] == '|':
                #i# Delete block plus 5' gap
                self.list[qh] = amap[:mapi[0]-1] + amap[mapi[1]+1:]
            ## ~ [2d] Central sequence extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                #i# Merge gaps. Note that a new gap will be created if the sequence is inserted somewhere
                #i# Total assembly length may therefore change
                #i# Delete block plus 3' gap
                self.list[qh] = amap[:mapi[0]] + amap[mapi[1]+2:]
            ### ~ [3] Keep fragment as extra sequence? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if keep:
                self.list[qh] += ['|'] + block + ['|']
                self.printLog('#EDIT','Extracted |{0}...{1}|'.format(flank1,flank2))
            else:
                self.printLog('#EDIT','Removed {0}...{1}'.format(flank1,flank2))
            return self.list[qh]
        except:
            self.errorLog('%s.mapExtraction error' % self.prog())
            return False
#########################################################################################################################
    def mapInsertion(self,qh,flank1,flank2,site,gapsize=0):    ### Relocate an assembly map block
        '''
        Relocate an assembly map block. If the block to be relocated is not a whole sequence (|block|) then it will
        first be extracted using mapExtraction().
        >> qh:qry/hit
        >> flank1:str = One end of block to move
        >> flank2:str = Other end of block to move
        >> site:str = Site of the insertion in the form "5'Flank:^:3'Flank"
        >> gapsize:int = Size of new gap to add, if not using loaded GapSize
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh]
            if not gapsize: gapsize = self.getInt('GapSize')
            ## ~ [0a] Check insertion site ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                [flank5,flank3] = site.split(':^:')
            except:
                raise ValueError('mapInsertion() site format not recognised: {0}'.format(site))
            mapi = [amap.index(flank5), amap.index(flank3)]
            mapi.sort()
            if mapi[1] != mapi[0]+2 or not amap[mapi[0]+1].startswith(':'):
                #i# Map insertion site no longer flanks a gap
                return False
            ## ~ [0b] Check insertion block ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapi = [amap.index(flank1), amap.index(flank2)]
            mapi.sort()
            if flank1 == flank2:
                mapi[1] = mapi[1] + 4
                if amap[mapi[1]] != flank2:
                    raise ValueError('Problem finding both ends of contig with identical flanks')
            try: self.checkSubMap(amap[mapi[0]:mapi[1]+1],amap[mapi[0]-1],amap[mapi[1]+1])
            except:
                self.errorLog('Insertion failed',quitchoice=False)
                return False
            block = amap[mapi[0]:mapi[1]+1]
            if flank5 in block or flank3 in block:
                self.warnLog('Trying to relocate assembly block {0}...{1} into itself! Possible adjacent inversions.'.format(flank1,flank2))
                return False
            ## ~ [0c] Extract if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if amap[mapi[0]-1] != '|' or amap[mapi[1]+1] != '|':
                amap = self.mapExtraction(qh,flank1,flank2,keep=True)
                if not amap:
                    return False
                mapi = [amap.index(flank1), amap.index(flank2)]
                mapi.sort()
                if block != amap[mapi[0]:mapi[1]+1]:
                    self.debug(' '.join(block))
                    raise ValueError('Extracted assembly block has changed during extraction!')
                if amap[mapi[0]-1] != '|' or amap[mapi[1]+1] != '|' or not self.checkSubMap(block):
                    self.warnLog('Extraction of {0}...{1} for re-insertion failed!'.format(flank1,flank2))
                    return False

            ### ~ [1] Delete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = amap[:mapi[0]-1] + amap[mapi[1]+2:]

            ### ~ [2] Insert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapi = [amap.index(flank5), amap.index(flank3)]
            mapi.sort()
            ## ~ [2a] Bwd insertion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if amap[mapi[0]] == flank3:
                block.reverse()
                for i in range(len(block)):
                    if block[i] == '>': block[i] = '<'
                    elif block[i] == '<': block[i] = '>'
            elif amap[mapi[0]] != flank5:
                raise ValueError('Block 5\' flank mismatch!')
            ## ~ [2b] Make insertion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            amap = amap[:mapi[1]] + block + [':Fix:{0}:'.format(gapsize)] + amap[mapi[1]:]
            self.list[qh] = amap
            self.printLog('#EDIT','Inserted {0}...{1} into {2}'.format(flank1,flank2,site))
            return self.list[qh]
        except:
            self.errorLog('%s.mapInsertion error' % self.prog())
            return False
#########################################################################################################################
    def mapBreak(self,qh,flank1,flank2):    ### Breaks assembly map between flank1 and flank2. Must be adjacent.
        '''
        Breaks assembly map between flank1 and flank2. Must be adjacent and either side of a gap.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh.lower()]
            if flank1 == flank2:
                raise ValueError('Cannot break between identical flanks!')
            ### ~ [1] Break ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapi = [amap.index(flank1), amap.index(flank2)]
            mapi.sort()
            if amap[mapi[0]+1][:1] != ':' or mapi[1] != (mapi[0] + 2):
                self.warnLog('Tried to break assembly between non-adjacent flanks {0} and {1}'.format(flank1,flank2))
                return False
            self.list[qh] = amap[:mapi[0]+1] + ['|','|'] + amap[mapi[1]:]
            self.printLog('#EDIT','Split {0} || {1}'.format(amap[mapi[0]],amap[mapi[1]]))
            return self.list[qh]
        except:
            self.errorLog('%s.mapBreak error' % self.prog())
            self.debug(amap[:20])
            return False
#########################################################################################################################
    def mapJoin(self,qh,flank1,flank2,expect=True): ### Joins flank1 3' to 5' of flank2, as long as they are at the end of sequences
        '''
        Joins flank1 3' to 5' of flank2, as long as they are at the end of sequences. Will invert a block if the flank is
        at the wrong end, unless both are in which case they will be reversed.
        :return: map/None/False - None if already joined!
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh.lower()]
            if flank1 == flank2:
                raise ValueError('Cannot join identical flanks!')
            if not flank1 or not flank2: return False
            gapsize = self.getInt('GapSize')
            newgap = [':Join:{0}:'.format(gapsize)]
            ### ~ [1] Check ends and establish blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Block1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapi = amap.index(flank1)
            i1 = j1 = mapi
            while amap[i1] != '|': i1 -= 1
            while amap[j1] != '|': j1 += 1
            if flank1 not in [amap[i1+1],amap[j1-1]]:
                mapi = [amap.index(flank1), amap.index(flank2)]
                mapi.sort()
                if amap[mapi[0]+1][:1] != ':' or mapi[1] != (mapi[0] + 2):
                    if expect:
                        self.warnLog('Tried to join assembly at internal flank {0}.'.format(flank1))
                    return False
                else:
                    #i# These flanks are already joined!
                    return None
            if flank2 in [amap[i1+1],amap[j1-1]]:
                if expect:
                    self.warnLog('Trying to make circular join between ends of same block: {0} and {1}'.format(flank1,flank2))
                return False
            block1 = amap[i1:j1+1]
            self.bugPrint(' '.join(block1))
            fwd1 = flank1 == block1[-2]
            ## ~ [1b] Block2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapi = amap.index(flank2)
            i2 = j2 = mapi
            while amap[i2] != '|': i2 -= 1
            while amap[j2] != '|': j2 += 1
            if flank2 not in [amap[i2+1],amap[j2-1]]:
                if expect:
                    self.warnLog('Tried to join assembly at internal flank {0}'.format(flank2))
                    self.debug(' '.join(amap[i2:j2+1]))
                return False
            block2 = amap[i2:j2+1]
            self.bugPrint(' '.join(block2))
            fwd2 = flank2 == block2[1]
            ## ~ [1c] Assembly chunks outside these blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if i1 < i2:
                chunk1 = amap[:i1]
                chunk2 = amap[j1+1:i2] + amap[j2+1:]
            else:
                chunk1 = amap[:i2]
                chunk2 = amap[j2+1:i1] + amap[j1+1:]
            self.bugPrint(' '.join(chunk1[-20:]))
            self.deBug(' '.join(chunk2[:20]))
            if chunk1 and chunk1[-1] != '|':
                raise ValueError('Problem with adjacent sequence not ending with terminus')
            if chunk2 and chunk2[0] != '|':
                raise ValueError('Problem with adjacent sequence not ending with terminus')
            ### ~ [2] Make joins at flank 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fwd1 and fwd2:
                amap = chunk1 + block1[:-1] + newgap + block2[1:] + chunk2
                self.printLog('#EDIT','Joined >{0}>::>{1}>'.format(flank1,flank2))
            elif fwd1:
                amap = chunk1 + block1[:-1] + newgap + self.reverseAssemblyChunk(block2)[1:] + chunk2
                self.printLog('#EDIT','Joined >{0}>::<{1}<'.format(flank1,flank2))
            elif fwd2:
                amap = chunk1 + self.reverseAssemblyChunk(block1)[:-1] + newgap + block2[1:] + chunk2
                self.printLog('#EDIT','Joined <{0}<::>{1}>'.format(flank1,flank2))
            else:
                amap = chunk1 + block2[:-1] + newgap + block1[1:] + chunk2
                self.printLog('#EDIT','Joined <{1}<::<{0}<'.format(flank1,flank2))
            self.list[qh] = amap
            return self.list[qh]
        except:
            self.errorLog('%s.mapJoin error' % self.prog())
            self.debug(amap[:20])
            return False
#########################################################################################################################
    ### <12> ### SynBad Output Methods                                                                                  #
#########################################################################################################################
    def saveTables(self,backup=True):   ### Saves the output tables, once all processing is complete.
        '''
        Saves the output tables, once all processing is complete.
        '''
        try:### ~ [1] Save gap and local tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.headLog('SAVE SYNBAD TABLES',line='=')
            for qh in ('Qry','Hit'):
                gdb = self.dbTable(qh,'gap')
                for field in ['MaxFlank5','MaxFlank3']:
                    if field in gdb.fields(): gdb.dropField(field)
                gdb.saveToFile(backup=backup)
                self.dbTable(qh,'corrections').saveToFile(backup=backup)
                locdb = self.db(qh.lower())
                locdb.dropField('AlnNum')
                outfields = locdb.fields()[0:]
                try:
                    outfields.remove('{0}5'.format(qh))
                    outfields.remove('{0}3'.format(qh))
                except:
                    pass    # Will need to update pairs tables!
                locdb.saveToFile(savefields=outfields,backup=backup)
                #X# This is now done earlier: check this OK - self.bestHiCPairTable(qh)

            ### ~ [2] Paired scaffold output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [2e] Top-matched pairs only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bestpair = self.getBool('BestPair')
            if not bestpair: return True
            pairs = []
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'top{0}'.format(qh.lower()),replace=True,add=True)
                topdb.keepFields(['Qry','Hit','Length'] + list(topdb.keys()))
                topdb.compress(['Qry','Hit'],default='sum')
                topdb.keepFields(['Qry','Hit','Length'])
                topdb.rankFieldByIndex(qh,'Length',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False,warn=True,highest=False)
                topdb.dropEntriesDirect('Rank',[1],inverse=True,log=True,force=False)
                pairs += list(topdb.dataKeys())
            ## ~ [2b] Output pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'{0}.pairs'.format(qh.lower()),replace=True,add=True)
                if bestpair:
                    ex = 0.0; etot = topdb.entryNum()
                    for ekey in list(topdb.datakeys())[0:]:
                        self.progLog('\r#PAIRS','Reducing %s table to top-aligned pairs: %.2f%%' % (qh,ex/etot)); ex += 100
                        entry = topdb.data(ekey)
                        if (entry['Qry'],entry['Hit']) not in pairs: topdb.dict['Data'].pop(ekey)
                    self.printLog('\r#PAIRS','Reduced %s table to %s alignments between top sequence pairs.' % (qh,rje.iStr(topdb.entryNum())))
                else:
                    if qh == 'Qry': topdb.newKey(['Qry','Hit','QryStart','QryEnd'])
                    else: topdb.newKey(['Hit','Qry','HitStart','HitEnd'])
                topdb.saveToFile(backup=backup)

            return True
        except: self.errorLog('%s.saveTables error' % self.prog()); return False
#########################################################################################################################
    def synBadSummarise(self):  ### Summarises ratings etc.
        '''
        Summarises ratings etc.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db().addEmptyTable('summary',['File','SeqNum','CtgNum','Good','Fix','Dupl','Bad','HiCBest','HiCPart','HiCWeak','HiCPoor','HiCNone','AltBest','AltPart','AltWeak','HiCScore']+rje.sortKeys(gapdesc),['File'],log=True)
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            #i# Add a summary of each class with descriptions and [Good/Fix/Dup/Bad] followed by percentage classes
            base = {'qry':rje.baseFile(self.getStr('Genome1'), strip_path=True),
                    'hit':rje.baseFile(self.getStr('Genome2'), strip_path=True)}
            ### ~ [1] Summary reports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                gapx = {'Total':0,'Good':0,'Fix':0,'Dupl':0,'Bad':0}
                #!# Add stats to a summary table and output: File SeqNum CtgNum Good Fix Dup Bad <SynBad Types>
                sentry = {'File':seqobj[qh].name(),'SeqNum':seqobj[qh].seqNum(),'CtgNum':self.dbTable(qh,'contigs').entryNum()}
                for field in sdb.fields()[3:]: sentry[field] = 0
                self.headLog('SYNBAD {0} {1}'.format(qh,base[qh]),line='-')
                gapdb = self.dbTable(qh,'gap')
                gapdb.index('SynBad',force=True)
                for gtype in gapdb.indexKeys('SynBad'):
                    gx = len(gapdb.index('SynBad')[gtype])
                    if gtype not in gapdesc: gapdesc[gtype] = 'See docs for details'
                    else: sentry[gtype] += gx
                    gclass = 'Bad'
                    if gtype in goodgaps: gclass = 'Good'
                    elif 'Fix' in gtype: gclass = 'Fix'
                    elif 'Dup' in gtype: gclass = 'Dupl'
                    self.printLog('#SYNBAD','{0} {1} ({2}): {3} [{4}]'.format(gx,gtype,qh,gapdesc[gtype],gclass))
                    gapx['Total'] += gx
                    gapx[gclass] += gx
                    sentry[gclass] += gx
                for gclass in ['Good','Fix','Dupl','Bad']:
                    self.printLog('#SYNBAD','{0:.1f}% {1} gaps in {2} SynBad class'.format(100.0*gapx[gclass]/gapx['Total'],base[qh],gclass))

                #i# Report numbers of reciprocal and partial best HiC pairs and sum HiCScore
                #!# Add Weak and None?

                #i# HiCBest = Gaps with Best HiC pairs
                #i# HiCPart = Gaps with partial best HiC pairs (one flank)
                #i# HiCWeak = Gaps with any weak best HiC pairs
                #i# HiCPoor = Gaps with some HiC support but alt best pairs
                #i# HiCNone = Gaps with no HiC support
                #i# AltBest = Number of mutual best alternative HiC pairings
                #i# AltPart = Number of partial best alternative HiC pairings
                #i# AltWeak = Number of weak best alternative HiC pairings

                #i# 'HiCBest','HiCPart','HiCWeak','HiCPoor','HiCNone','AltBest','AltPart','AltWeak'


                hdb = self.addTable(qh,'hicbest',expect=False,make=False)
                if hdb:
                    gapdb.dataFormat({'HiCScore':'num'})
                    hicscores = gapdb.orderedDataList('HiCScore',empties=False)
                    while '' in hicscores: hicscores.remove('')
                    while 'NA' in hicscores: hicscores.remove('NA')
                    while '?' in hicscores: hicscores.remove('?')
                    hicscores = list(map(float,hicscores))
                    sentry['HiCScore'] = sum(hicscores)
                    for gentry in gapdb.entries():
                        flanks = [gentry['GapFlank5'],gentry['GapFlank3']]
                        flanks.sort()
                        hentry = hdb.data((flanks[0],flanks[1]))
                        htype = 'HiCNone'
                        try: 
                            if gentry['HiCScore'] > 0: htype = 'HiCPoor'
                        except:
                            self.warnLog('Non-numerical HiCScore for {0}: {1}'.format(gentry['GapName'],gentry['HiCScore']))
                        if hentry:
                            if hentry['Best'] == 'Both': htype = 'HiCBest'
                            elif hentry['Best'].startswith('Flank'): htype = 'HiCPart'
                            elif hentry['Best'].startswith('Weak'): htype = 'HiCWeak'
                        sentry[htype] += 1
                    for hentry in hdb.indexEntries('Type','Check'):
                        htype = None
                        if hentry['Best'] == 'Both': htype = 'AltBest'
                        elif hentry['Best'].startswith('Flank'): htype = 'AltPart'
                        elif hentry['Best'].startswith('Weak'): htype = 'AltWeak'
                        if htype: sentry[htype] += 1
                    #!# Add overall HiC summary
                    for hclass in ['HiCBest','HiCPart','HiCWeak','HiCPoor','HiCNone']:
                        self.printLog('#SYNBAD','{0:.1f}% ({1}) {2} gaps in {3} HiC class'.format(100.0*sentry[hclass]/gapx['Total'],sentry[hclass],base[qh],hclass))

                # #!# This needs to ignore the Check entries #!# (Report separately?)
                # if hdb:
                #     hdb.index('Best')
                #     if 'Flank1' in hdb.index('Best'): sentry['HiCPart'] += len(hdb.index('Best')['Flank1'])
                #     if 'Flank2' in hdb.index('Best'): sentry['HiCPart'] += len(hdb.index('Best')['Flank2'])
                #     if 'Both' in hdb.index('Best'): sentry['HiCBest'] += len(hdb.index('Best')['Both'])

                sdb.addEntry(sentry)

            sdb.saveToFile()
            return True
        except:
            self.errorLog('%s.synBadSummarise error' % self.prog())
            return False
#########################################################################################################################
    def mapFragment(self):  ### SynBad fragmentation from assembly map
        '''
        SynBad fragmentation from assembly map.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            badgaptypes = self.list['FragTypes']
            ## ~ [0a] Check for existing output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            wanted = ['qry.frag.txt','hit.frag.txt','qry.frag.fasta','hit.frag.fasta','qry.frag.tdt','hit.frag.tdt']
            if not self.update() and rje.checkForFiles(wanted,basename=self.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Fragmented SynBad assembly files found: skipping fragmentation (force=F update=F)')
                return True
            self.printLog('#FRAG','Fragmenting on gap types: {0}'.format(', '.join(badgaptypes)))
            self.printLog('#FRAG','Min read span depth to prevent fragmentation: {0}'.format(self.getInt('MinReadSpan')))

            ### ~ [1] Update corrections table as fragment, update assembly map and output ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                self.headLog('Identifying SynBad {0} gaps to fragment'.format(qh),line='-')
                gapdb = self.addTable(qh,'gap',expect=True,make=False)
                fixdb = db.copyTable(self.addTable(qh,'corrections',expect=True,make=True),'{0}.frag'.format(qh),replace=True,add=True)
                fragx = 0; failx = 0
                for entry in gapdb.entries():
                    self.bugPrint(gapdb.entrySummary(entry,collapse=True))
                    if entry['Span0'] < self.getInt('MinReadSpan') and entry['SynBad'] in badgaptypes:
                        frag = {'SeqName':entry['SeqName'],'Start':entry['Start'],'End':entry['End'],
                                'Flank1':entry['GapFlank5'],'Flank2':entry['GapFlank3'],'Edit':'Fragment',
                                'Details':'{0} ({1})'.format(entry['GapName'],entry['SynBad']),'SynBad':'Fragment'}
                        if not self.mapBreak(qh,entry['GapFlank5'],entry['GapFlank3']): frag['SynBad'] = 'Failed'; failx += 1
                        fixdb.addEntry(frag); fragx += 1
                self.printLog('#FRAG','{0} {1} gaps identified for fragmentation; {2} failed.'.format(fragx,qh,failx))
                self.saveAssemblyMaps(qh,mapname='frag')
                self.saveTelociraptorMaps(qh,mapname='frag')

            return True
        except:
            self.errorLog('%s.mapFragment error' % self.prog())
            return False
#########################################################################################################################
    def fragment(self):   ### SynBad fragmentation
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqObjSetup()
            seqobj = {}
            seqobj['qryfrag'] = self.obj['SeqList1']
            seqobj['hitfrag'] = self.obj['SeqList2']
            db = self.db()
            basefile = self.baseFile()
            qrygap = self.addTable('qry','gap',expect=True,make=False)
            qryfrag = db.copyTable(qrygap,'qry.frag',replace=True,add=True)
            hitgap = self.addTable('hit','gap',expect=True,make=False)
            hitfrag = db.copyTable(hitgap,'hit.frag',replace=True,add=True)
            ### ~ [1] Filter gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
            # * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
            # * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
            # * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
            # * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
            # * `Inv` = `Inversion` = Flanking hits are on alternative strands.
            # * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
            # * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
            # * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.
            # * `Span` = `Spanned` = Any gaps with Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
            # * `Null` = No mapping between genomes for that gap.
            badgaptypes = self.list['FragTypes']
            self.headLog('SYNBAD GAP FILTER',line='=')
            for table in (qryfrag,hitfrag):
                table.addField('FragID',evalue='')
                for entry in table.entries():
                    if entry['Span0'] < self.getInt('MinReadSpan') and entry['SynBad'] in badgaptypes:
                        entry['FragID'] = 'Filter'
                table.dropEntriesDirect('FragID',['Filter'],inverse=True)
                table.newKey(['SeqName','Start','End','FragID'])
                table.addField('Syn5',evalue='')
                table.addField('Syn3',evalue='')
            ### ~ [2] Fragment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('FRAGMENT ON GAPS',line='=')
            for frag in ('qry.frag','hit.frag'):
                table = self.db(frag)
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                seqlist = seqobj[frag]
                prev = {'SeqName':None}
                fragname = ''; fx = 1
                for ekey in list(table.dataKeys()):
                    entry = table.data(ekey)
                    #i# Add extra 5' fragment for each sequence
                    if entry['SeqName'] != prev['SeqName']:
                        prev = rje.combineDict({},entry)
                        prev['Start'] = 1
                        prev['End'] = entry['Start'] - 1
                        fragname = entry['SeqName']; fx = 1
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        prev['FragID'] = '{0}.{1}'.format(fragname,fx); fx += 1
                        prev['Syn5'] = 'Terminal'
                        prev['Syn3'] = 'Terminal'    #i# May be over-written
                        prev = table.addEntry(prev,remake=False)
                    #i# Update details of this entry
                    prev['Syn3'] = entry['SynBad']
                    prev['End'] = entry['Start'] - 1
                    entry['Start'] = entry['End'] + 1
                    entry['FragID'] = '{0}.{1}'.format(fragname,fx); fx += 1
                    entry['Syn5'] = entry['SynBad']
                    entry['Syn3'] = 'Terminal'                          #i# May be over-written
                    entry['End'] = entry['SeqLen']                      #i# May be overwritten
                    prev = entry
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    #if entry['SeqName'] == 'NSCUSCAFF155': self.debug('?')
                #i# Add sequences without any gaps
                table.remakeKeys()
                seqnames = rje.sortKeys(table.index('SeqName'))
                self.printLog('#FRAG','{0} fragments from {1} sequences with gaps'.format(rje.iStr(table.entryNum()),rje.iLen(seqnames)))
                sx = 0
                for seq in seqlist.seqs():
                    seqname = seqlist.shortName(seq)
                    if seqname not in seqnames:
                        fragname = seqname
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        fragid = '{0}.1'.format(fragname)
                        table.addEntry({'FragID':fragid,'SeqName':seqname,'Start':1,'SeqLen':seqlist.seqLen(seq),'End':seqlist.seqLen(seq),'Syn5':'Terminal','Syn3':'Terminal'})
                        sx += 1
                self.printLog('#SEQ','Added {0} full-length sequences without gaps'.format(rje.iStr(sx)))
                #i# Tidy up fragment table
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                #table.newKey(['FragID'])
                table.setFields(['FragID','SeqName','Start','End','SeqLen','Syn5','Syn3'])
                table.saveToFile()
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
            ### ~ [3] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                seqdict = seqlist.seqNameDic()
                fragfas = '{0}.{1}.fasta'.format(basefile,frag)
                rje.backup(self,fragfas)
                FRAGFAS = open(fragfas,'w'); ex = 0
                for entry in table.entrySort():
                    self.progLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(ex),fragfas)); ex += 1
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    seq = seqdict[entry['SeqName']]
                    (seqname,sequence) = seqlist.getSeq(seq)
                    FRAGFAS.write('>{0}\n'.format(entry['FragID']))
                    FRAGFAS.write(sequence[entry['Start']-1:entry['End']]+'\n')
                FRAGFAS.close()
                self.printLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(table.entryNum()),fragfas))
            return True
        except:
            self.errorLog('%s.fragment error' % self.prog())
            return False
#########################################################################################################################
### End of SECTION II: SynBad Class                                                                                     #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: SynBad(mainlog,['basefile=synbad']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
