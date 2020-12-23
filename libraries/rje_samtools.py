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
Version:      1.20.3
Last Edit:    06/06/19
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

    If `depthplot=T` then read depth will be calculated across each locus and additional depth statistics (MinX, MaxX,
    MedianX and Coverage) will be added to the coverage file. If `rgraphics=T` (default) and R is installed, plots per
    locus will be generated. If `readlen=T` then additional outputs will be generated of the maximum read length at each
    point in the genome.

    If `dirlen=X` > 0 (default=500), the maximum 5' and 3' read linkage will be calculated every `dirlen=X` nucleotides
    (plus the end position). The examines the reads spanning that position and reports the maximum distance any read
    extends in the 5' (Len5) or 3' (Len3) direction from that position. This can identify misassemblies resulting in
    breaks in read coverage. Note that reads in this context are contiguous SAM/pileup hits and so sequencing reads with
    multiple or split mapping will be treated as multiple reads.

Commandline:
    ### ~ MPileup Parsing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    batch=FILELIST  : List of MPileup/SAM files to be parsed and filtered (e.g. *.pileup) []
    qcut=X          : Min. quality score for a call to include [30]
    minqn=X         : Min. number of reads meeting qcut (QN) for output [10]
    rid=T/F         : Whether to include Read ID (number) lists for each allele [True]
    snponly=T/F     : Whether to restrict parsing output to SNP positions (will use mincut settings below) [False]
    indels=T/F      : Whether to include indels in "SNP" parsing [True]
    skiploci=LIST   : List of loci to exclude from pileup parsing (e.g. mitochondria) []
    snptableout=T/F : Output filtered alleles to SNP Table [False]
    ### ~ SNP Frequency Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    control=FILE    : MPileup or processed TDT file to be used for control SNP frequencies []
    treatment=FILE  : MPileup or processed TDT file to be used for treatment SNP frequencies []
    labels=X,Y      : Optional labels for Control and Treatment fields in output (other file basename used) []
    mincut=X        : Minimum read count for minor allele (proportion if <1) [1]
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
    rgraphics=T/F   : Whether to generate snpfreq multichromosome plots [True]
    ### ~ Double Genome Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    altcontrol=FILE     : MPileup or processed TDT file to be used for control SNP frequencies in Alt genome []
    alttreatment=FILE   : MPileup or processed TDT file to be used for treatment SNP frequencies in Alt genome []
    ### ~ Read Coverage Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    readcheck=FILE      : BAM/SAM/Pileup/RID file with read mappings [None]
    seqin=FASFILE       : Sequence file for loci in MPileup/SAM files (e.g. matching all relevant RID files) [None]
    depthplot=T/F       : Whether to generate Xdepth plot data for the readcheck FILE. (May be slow!) [False]
    fullcut=X           : Proportion of read to be mapped to count as full-length [0.9]
    readlen=T/F         : Include read length data for the readcheck file (if depthplot=T) [True]
    dirnlen=X           : Include directional read length data at X bp intervals (readlen=T; 0=OFF) [500]
    minsoftclip=INT     : If set, will generate *.clipped.* output for reads with terminal soft clipping above INT [0]
    maxsoftclip=INT     : If set, will generate *.noclip.* output for reads with less terminal soft clipping than INT [0]
    minreadlen=INT      : Minimum read length to include in depth plots (all reads will be in RID file) [0]
    depthsmooth=X       : Smooth out any read plateaus < X nucleotides in length [200]
    peaksmooth=X        : Smooth out Xcoverage peaks < X depth difference to flanks (<1 = %Median) [0.05]
    rgraphics=T/F       : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
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
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
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
    # 1.9.0 - Added depthplot data generation. (Will need to add R function for plot itself.)
    # 1.9.1 - Changed mincut default to 0.1.
    # 1.10.0 - Added readlen output, which is like the depth plot but uses max read length (kb) instead of depth.
    # 1.11.0 - Added dirnlen=X : Include directional read length data at X bp intervals (depthplot=T; 0=OFF) [500]
    # 1.11.1 - Minor tweaks to try and speed up pileup parsing.
    # 1.12.0 - Updated the snpfreq run code to make clearer and check for parsing issues. Set mincut=1 default.
    # 1.13.0 - Added skiploci=LIST - need to screen out mitochondrion from Illumina Pileup parsing!
    # 1.14.0 - Added forking of pileup parsing for SNPFreq analysis.
    # 1.14.1 - Fixed SNPFreq rerunning bug.
    # 1.15.0 - Added rgraphics=T/F   : Whether to generate snpfreq multichromosome plots [True]
    # 1.16.0 - Add coverage calculation per locus to depth plot table output (depthplot=T).
    # 1.16.1 - Added reporting of existing files for parsing Pileup.
    # 1.17.0 - Added parsing of lengths from SAM files to RID file.
    # 1.18.0 - Updated processing of Treatment and Control without Alt to still limit to SNPTable. Fixed SNPFreq filters.
    # 1.19.0 - snptableout=T/F    : Output filtered alleles to SNP Table [False]
    # 1.19.1 - Fixed AltLocus SNP table bug.
    # 1.19.2 - Updated forker parsing to hopefully fix bug.
    # 1.20.0 - Added parsing of BAM file - needs samtools on system. Added minsoftclip=X, maxsoftclip=X and minreadlen=X.
    # 1.20.1 - Fixed mlen bug. Added catching of unmapped reads in SAM file. Fixed RLen bug. Changed softclip defaults.
    # 1.20.2 - Fixed readlen coverage bug and acut bug.
    # 1.20.3 - Fixed RLen bug.
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
    # [ ] : Add checking and resorting of RID file to match sorted BAM file versus raw SAM file.
    # [ ] : Not sure if above strictly possible, so may need to stick with checking RID versus pileup file.
    # [ ] : Should there also be DirnDepth plots of different lengths?
    # [ ] : Add skiploci=LIST - need to screen out mitochondrion from Illumina Pileup parsing!
    # [ ] : snpsubset=LIST  : List of SNPs to restrict analysis to, based on subsetfield=X (comma separated  []
    # [ ] : subsetfield=X   : SNPTable field to be used for matching snpsubset list []
    # [ ] : subsample=X : Option to subsample every X bases for read depth summaries etc.
    # [ ] : Find out why SAMPhaser parsing skip failed.
    # [Y] : Add length proportion filters for read depths from SAM files.
    # [ ] : Add minlength cutoffs for read depth analysis (minrlen and minmlen). [Q. Put in output filename?]
    # [ ] : Tidy, clarify and document re-use of files and use of strip_path. Set OutDir to overcome?
    # [ ] : Swap the order of X versus Y - should be Treatment versus Control.
    # [Y] : Add MinSoftClip=INT = minimum length of terminal soft-clipping allowed for SAM/BAM parsing (Also MaxSoftClip)
    # [ ] : Add scdepth=INT [0] to over-ride use of modal read depth as the 1n line on depth plots.
    # [ ] : Add WIG format:
    variableStep  chrom=chrN
    [span=windowSize]
    chromStartA  dataValueA
    chromStartB  dataValueB
    ... etc ...  ... etc ...
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('rje_samtools', '1.20.3', 'June 2019', '2013')
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
    - DepthPlot=T/F       : Whether to generate Xdepth plot data for the readcheck FILE. (May be slow!) [False]
    - IgnoreN = Whether to exclude "N" calls for major/minor alleles [True]
    - IgnoreRef=T/F : If False will always keep Reference allele and call fixed change as SNP [False]
    - Indels=T/F      : Whether to include indels in "SNP" parsing [True]
    - MajDif = Whether to restrict output and stats to positions with Major Allele differences [True]
    - MajFocus=T/F    : Whether the focus is on Major Alleles (True) or Mutant/Reference Alleles (False) [True]
    - MajMut = Whether to restrict output and stats to positions with non-reference Major Allele [True]
    - MajRef = Whether to restrict output and stats to positions with non-reference Major Allele [False]
    - ReadLen=T/F         : Include read length data for the readcheck file (if depthplot=T) [True]
    - RGraphics=T/F       : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
    - RID=T/F         : Whether to include Read ID (number) lists for each allele [False]
    - SNPFreq=T/F       : Whether this is a snpfreq run: set snponly=F
    - SNPOnly=T/F     : Whether to restrict parsing output to SNP positions (will use mincut settings below) [False]
    - SNPTableOut=T/F    : Output filtered alleles to SNP Table [False]

    Int:integer
    - AbsMinCut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    - DepthSmooth=X   : Smooth out any read plateaus < X nucleotides in length [20]
    - DirnLen=X       : Include directional read length data at X bp intervals (depthplot=T; 0=OFF) [500]
    - FDRCut=X        : Additional FDR cutoff for enriched treatment SNPs [1.0]
    - MajCut=X        : Frequency cutoff for Treatment major allele [0.0]
    - MinCut=X        : Minimum read count for minor allele (proportion if <1) [1]
    - MinReadLen=INT      : Minimum read length to include in depth plots (all reads will be in RID file) [0]
    - MinSoftClip=INT : If set, will generate *.clipped.* output for reads with terminal soft clipping above INT [0]
    - MaxSoftClip=INT : If set, will generate *.noclip.* output for reads with less terminal soft clipping than INT [0]
    - MinQN = Min. number of reads meeting qcut (QN) for output [10]
    - PeakSmooth=X        : Smooth out Xcoverage peaks < X depth difference to flanks [5]
    - QCut = Min. quality score for a call to include [30]
    - SigCut=X        : Significance cutoff for enriched treatment SNPs [0.05]

    Num:float
    - FullCut=X           : Proportion of read to be mapped to count as full-length [0.9]
    - MinFreq = Minor allele(s) frequency correction for zero counts (e.g. Sequencing error) [0.01]
    
    List:list
    - Batch=FILELIST  : List of MPileup files to be parsed and filtered []
    - CheckFields=LIST    : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    - CheckFlanks=LIST       : List of lengths flanking check regions that must also be spanned by reads [0,100,500,1000]
    - Labels=X,Y      : Optional labels for Control and Treatment fields in output (other file basename used) []
    - SkipLoci=LIST   : List of loci to exclude from pileup parsing (e.g. mitochondria) []
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
        self.boollist = ['Biallelic','Child','DepthPlot','ReadLen','RGraphics','IgnoreN','IgnoreRef','Indels','MajDif',
                         'MajFocus','MajMut','MajRef','RID','SNPOnly','SNPFreq','SNPTableOut']
        self.intlist = ['AbsMinCut','CheckFlanks','DepthSmooth','DirnLen','MaxSoftClip','MinSoftClip','MinQN','PeakSmooth','QCut']
        self.numlist = ['FDRCut','FullCut','MinCut','MinFreq','SigCut']
        self.listlist = ['Batch','CheckFields','CheckFlanks','Labels','SkipLoci','SNPTabKeys','SNPTabMap']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'RefSeq':''})
        self.setBool({'Biallelic':False,'Child':False,'DepthPlot':False,'ReadLen':True,'RGraphics':True,'IgnoreN':True,
                      'IgnoreRef':False,'Indels':True,'MajDif':False,'MajFocus':True,'MajMut':True,'MajRef':False,
                      'RID':True,'SNPTableOut':False})
        self.setInt({'DepthSmooth':200,'DirnLen':500,'PeakSmooth':0.05,'QCut':30,'MinQN':10,'AbsMinCut':2,
                     'MaxSoftClip':0,'MinSoftClip':0,'MinReadLen':0})
        self.setNum({'FullCut':0.9,'MajCut':0.0,'MinFreq':0.001,'MinCut':1,'SigCut':0.05,'FDRCut':1.0})
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
                self._cmdReadList(cmd,'bool',['Biallelic','DepthPlot','IgnoreN','IgnoreRef','Indels','MajDif','MajFocus',
                                              'MajMut','MajRef','ReadLen','RGraphics','RID','SNPOnly','SNPTableOut'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['AbsMinCut','DepthSmooth','DirnLen','MaxSoftClip','MinSoftClip','MinReadLen','PeakSmooth','QCut','MinQN'])   # Integers
                self._cmdReadList(cmd,'float',['FullCut','MinCut','FDRCut','MajCut','MinFreq','SigCut']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CheckFields','Labels','SkipLoci','SNPTabKeys','SNPTabMap'])  # List of strings (split on commas or file lines)
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
        #i# Update SoftClip settings
        if self.getInt('MinSoftClip') < 0:
            self.setInt({'MinSoftClip':0})
            self.warnLog('Cannot set minsoftclip<0: reset minsoftclip=0')
        if self.getInt('MaxSoftClip') < 0:
            self.setInt({'MaxSoftClip':0})
            self.warnLog('Cannot set maxsoftclip<0: reset maxsoftclip=0')
        # Forking
        if self.getBool('Win32') or self.getBool('NoForks'): self.setInt({'Forks':1})
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
                else:
                    snpfile = self.pileUpStats(snpdb=self.db('SNP'))
                fdrdb = self.rateSNPs(snpfile)
                self.combineSNPs(fdrdb)
                if self.getBool('RGraphics'):
                    self.rScriptGraphics()
                    self.rScriptGraphics(rargs='fdrplot=T')
            if self.getStrLC('ReadCheck'):
                ridfile = self.getStr('ReadCheck')
                ridbase = rje.baseFile(ridfile,strip_path=True)
                if not self.baseFile(return_none=None): self.baseFile(ridbase)
                if not ridfile.endswith('.rid.tdt'):
                    self.setBool({'RID':True})
                    self.parsePileup(ridfile)
                    ridfile = '%s.rid.tdt' % (ridbase)  #?# Should this strip path? Probably not!
                    if not rje.exists(ridfile):
                        ridbase = rje.baseFile(self.getStr('ReadCheck'),strip_path=False)
                        ridfile = '%s.rid.tdt' % (ridbase)  #?# Should this strip path? Probably not!
                    if not rje.exists(ridfile):
                        self.warnLog('Cannot find "%s" (or path-stripped version)' % ridfile)
                        ridfile = None
                if self.db('rid'): ridfile = None   # Use existing data
                self.coverageFromRID(ridfile,depthplot=self.getBool('DepthPlot'))
                if self.db('rid'): ridfile = None   # Use existing data
                if self.getBool('DepthPlot') and self.getBool('ReadLen'): self.coverageFromRID(ridfile,depthplot=True,readlen=True)
                #i# Soft-clipped read depth
                if self.getInt('MinSoftClip') > 0:
                    if self.db('rid'): ridfile = None   # Use existing data
                    self.coverageFromRID(ridfile,depthplot=self.getBool('DepthPlot'),clip=True)
                    if self.db('rid'): ridfile = None   # Use existing data
                    if self.getBool('DepthPlot') and self.getBool('ReadLen'): self.coverageFromRID(ridfile,depthplot=True,readlen=True,clip=True)
                    self.db().setBasefile(self.basefile())
                #i# No-clipped read depth
                if self.getInt('MaxSoftClip') > self.getInt('MinSoftClip') > 0:
                    self.printLog('#CLIP','maxsoftclip=INT > minsoftclip=INT: will use for soft-clipping range; no *.noclip.* output.')
                elif self.getInt('MaxSoftClip') > 0:
                    if self.db('rid'): ridfile = None   # Use existing data
                    self.coverageFromRID(ridfile,depthplot=self.getBool('DepthPlot'),noclip=True)
                    if self.db('rid'): ridfile = None   # Use existing data
                    if self.getBool('DepthPlot') and self.getBool('ReadLen'): self.coverageFromRID(ridfile,depthplot=True,readlen=True,noclip=True)
                    self.db().setBasefile(self.basefile())
            #!# Add different rGraphics modes/options - currently in coverageFromRID()
            #!# Add HTML output linking to graphics
        except SystemExit:
            if self.getBool('Child'): sys.exit(0)
            else: raise ValueError
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
                    self.setBool({'SNPFreq':True})
                    snpdb = self.loadSNPTable()
                    #i# The SNPOnly/SNPFreq mismatch is now handled by file naming and handling.
                    #if self.getBool('SNPOnly'):
                    #    self.warnLog('snponly=T is incompatible with SNPFreq analysis. Setting snponly=F. Could be issues if pileup files previously parsed with snponly=T')
                    #    self.setBool({'SNPOnly':True})
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
            for locus in self.list['SkipLoci']: self.printLog('#SKIP','Skipping locus: %s' % locus)
            if self.getInt('Forks') > 1:
                outfile = self.forkerParsePileUp()
            else:
                outfile = {}
                for pfile in self.list['Batch']:
                    outfile[pfile] = self.parsePileup(pfile)
            ### ~ [3] SNPTables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add forking of this too #!#
            if self.getBool('SNPTableOut'):
                for pfile in self.list['Batch']:
                    if not outfile or not outfile[pfile]: continue
                    self.snpTable(pfile,outfile[pfile])
            return True
        except SystemExit:
            if self.getBool('Child'): sys.exit(0)
            else: raise ValueError
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def snpTable(self,pfile,snpfile):   ### Generates SNP Table from parsed output file
        '''
        Generates SNP Table from parsed output file
        :param pfile:
        :param snpfile:
        :return:
        '''
        try:## ~ [2b] Filter based on allele frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            snpdb = db.openTable(snpfile,mainkeys=['Locus','Pos'],name='%s.snp' % pfile,expect=True)
            snptabdb = db.addEmptyTable('snptable',['Locus','Pos','REF','ALT','N','Freq'],['Locus','Pos','REF','ALT'],log=True)
            snpx = 0  # Number of SNPs read from SNP table
            snpend = rje.endPos(snpdb.obj['File'])
            while snpdb:
                entry = snpdb.readEntry(add=False,close=True)
                if not entry: break
                snpx += 1
                alleles = {}
                for aseq in string.split(entry['Seq'],'|'):
                    [a,n] = string.split(aseq,':')
                    n = int(n)
                    #!# Add filtering criteria
                    #if self.getNum('PhaseCut') < 1.0:
                    #    if n / float(entry['QN']) < self.getNum('PhaseCut'): continue
                    #    if n < self.getInt('AbsPhaseCut'): continue
                    #elif n < self.getInt('PhaseCut'): continue
                    #if (a[0] == '-' or len(a) > 1) and not self.getBool('PhaseIndels'): continue
                    #!#
                    alleles[a] = n
                #?# Add biallelic toggle?
                #if len(alleles) != 2: continue   # Not biallelic: Need to filter better!
                entry['Seq'] = []
                for a in rje.sortKeys(alleles):
                    if entry['Ref'] != a:
                        allentry = {'Locus':entry['Locus'],'Pos':int(entry['Pos']),'REF':entry['Ref'],'ALT':a,'N':alleles[a],
                                    'Freq':float(alleles[a])/sum(alleles.values())}
                        snptabdb.addEntry(allentry)
                self.progLog('#SNP','Parsing SNPs: %.2f%%; %s Pos -> %s SNP alleles.' % (100.0 * snpdb.obj['File'].tell() / snpend,rje.iStr(snpx),rje.iStr(snptabdb.entryNum())))
            tabfile = '%s.snp.tdt' % rje.baseFile(snpfile,strip_path=True)
            #!# Add more info on filters?
            self.printLog('#SNP','%s SNP alleles parsed from %s SNP into %s.' % (rje.iStr(snptabdb.entryNum()),rje.iStr(snpx),tabfile))
            snptabdb.saveToFile(tabfile)

            db.deleteTable(snpdb)
            db.deleteTable(snptabdb)
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
            if rje.isYounger(filename,ridfile) == filename:     ### Returns younger file or None if either does not exist
                self.printLog('#OLDRID','%s found but older than %s!' % (ridfile,filename))
                rje.backup(self,ridfile)
            if rje.exists(ridfile) and not self.force():
                self.printLog('#SKIP','%s found! (force=F)' % ridfile)
                return True
            if not rje.exists(filename):
                self.warnLog('SAM file "%s" missing. (May not matter if already parsed.)' % filename)
                return False
            self.printLog('#~~#','## ~~~~~ Parsing SAM File: %s ~~~~~ ##' % filename)
            #!# Possible parsing of lengths from read files #!#
            RIDOUT = open(ridfile,'w')
            rje.writeDelimit(RIDOUT,outlist=['RID','Locus','Start','End','RLen','MLen','Clip5','Clip3'],delimit='\t')
            #if filename.endswith('.bam'): SAM = os.popen('samtools view -h %s -o sam -@ 16' % filename)
            if filename.endswith('.bam'): SAM = os.popen('samtools view %s' % filename)
            else: SAM = open(filename,'r')
            rid = 0          # Read counter (ID counter)
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            oldblasrwarn = ''; xwarnx = 0; unmappedx = 0
            for line in SAM:
                ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.startswith('@'): continue
                samdata = string.split(line)
                if len(samdata) < 11: continue
                self.progLog('\r#PARSE','Parsing %s: %s reads (+ %s unmapped)...' % (filename,rje.iStr(rid),rje.iStr(unmappedx)),rand=0.1)
                cigstr = samdata[5]
                if cigstr == '*': unmappedx += 1; continue
                rid += 1
                rname = samdata[0]
                zdata = string.split(string.replace(rname,'/',' '))
                try:
                    [smrt,zmw,pos] = zdata[:3]
                    [beg,end] = string.split(pos,'_')
                    rlen = int(end) - int(beg) + 1
                except: rlen = 0
                locus = samdata[2]
                rpos = int(samdata[3])
                clip5 = 0
                if rje.matchExp('^(\d+)S',cigstr): clip5 = int(rje.matchExp('^(\d+)S',cigstr)[0])
                clip3 = 0
                if rje.matchExp('(\d+)S$',cigstr): clip5 = int(rje.matchExp('(\d+)S$',cigstr)[0])
                mlen = len(samdata[9]) - clip5 - clip3
                if len(samdata) > 22:
                    oldblasr = samdata[13].startswith('XS:i:')
                    oldblasr = oldblasr and samdata[14].startswith('XE:i:')
                    oldblasr = oldblasr and samdata[18].startswith('XL:i:')
                    oldblasr = oldblasr and samdata[22].startswith('XQ:i:')
                    if oldblasr:
                        if not oldblasrwarn: oldblasrwarn = 'Old BLASR SAM format detected: extracted clipping from XS,XE,XL and XQ info.'
                        XS = int(samdata[13][5:])
                        XE = int(samdata[14][5:])
                        XL = int(samdata[18][5:])
                        XQ = int(samdata[22][5:])
                        #if (XS + XL) != XE: xwarnx += 1
                        if rlen:
                            if XQ != rlen: xwarnx += 1
                        else: rlen = XQ
                        clip5 = XS-1
                        clip3 = XQ-XE+1
                        mlen = len(samdata[9])
                    elif not rlen:
                        rlen = len(samdata[9])
                        mlen = len(samdata[9]) - clip5 - clip3
                #!# Add filter for terminal nnnnS softclipping > self.getInt('MinSoftClip')
                try:
                    cigdata = parseCigar(cigstr)
                except:
                    self.errorLog('Problem with cigar string (RID:%d): "%s"' % (rid,cigstr))


                rje.writeDelimit(RIDOUT,outlist=[rid,locus,rpos,rpos+cigarAlnLen(cigdata)-1,rlen,mlen,clip5,clip3],delimit='\t')
            self.printLog('#RID','Parsed %s: %s read start/end positions output to %s' % (filename,rje.iStr(rid),ridfile))
            if unmappedx: self.printLog('#SAM','Parsed and rejected %s unmapped reads.' % (rje.iStr(unmappedx)))
            if oldblasrwarn: self.warnLog(oldblasrwarn)
            if xwarnx: self.warnLog('%s reads did not have compatible XS+XL=XE data' % rje.iStr(xwarnx))
            RIDOUT.close()
            SAM.close()
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
            alt = False     # Whether this is an Alt pileup file (else Ref).
            if filename.endswith('.sam'): return self.parseSAM(filename)
            if filename.endswith('.bam'): return self.parseSAM(filename)
            qtxt = 'Q%d.%d' % (self.getInt('QCut'),self.getInt('MinQN'))
            outbase = '%s.%s' % (rje.baseFile(filename,strip_path=True),qtxt)
            if self.getBool('SNPFreq'):
                outfile = '%s.snpfreq.tdt' % (outbase)
                snpdb = self.loadSNPTable()
                if self.getBool('RID'):
                    self.setBool({'RID':False})
                    self.printLog('#RID','RID allele output switched off for SNPFreq analysis.')
            elif self.getBool('SNPOnly'):
                suffix = 'b%si%s' % (str(self.getBool('Biallelic'))[:1],str(self.getBool('Indels'))[:1])
                outfile = '%s.snponly.%s.tdt' % (outbase,suffix)
            else: outfile = '%s.tdt' % (outbase)
            ridfile = '%s.rid.tdt' % (rje.baseFile(filename,strip_path=True))
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStr(fstr) == filename:
                    self.setStr({fstr:outfile})
                    if fstr.startswith('Alt'): alt = True
            #i# Assess whether the parsing needs to be performed.
            #i# Independently log whether to remake RID file: should be the same for all runs
            skiprun = not self.force()
            ridout = self.getBool('RID') and not rje.exists(ridfile)  # Always output RID file for safety if RID mapping used
            ridout = ridout and not self.getBool('SNPFreq')
            if self.getBool('RID'): self.printLog('#RID','%s: %s' % (ridfile,rje.exists(ridfile)))
            if self.getBool('RID') and rje.exists(ridfile) and not ridout:
                self.warnLog('Using existing RID file. If program crashes try deleting *.rid.tdt or run force=T.')
            skiprun = skiprun and not ridout    #(rje.exists(ridfile) or self.getBool('SNPFreq'))
            skiprun = skiprun and rje.exists(outfile)
            self.printLog('#OUT','%s: %s' % (outfile,rje.exists(outfile)))
            if skiprun:
                #!# Should we add something to check integrity/completeness? #!#
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return True
            self.printLog('#~~#','## ~~~~~ Parsing PileUp File: %s ~~~~~ ##' % filename)
            ## ~ [1a] Calculate mean read depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Currently this is based on unfiltered 'N' values.
            PILEUP = open(filename,'r'); px = 0; rx = 0
            for line in PILEUP:
                data = string.split(rje.chomp(line))
                if not data: break
                if data[0] in self.list['SkipLoci']:
                    self.progLog('\r#DEPTH','Skipping %s depth: %s pos; %s read bases' % (filename,rje.iStr(px),rje.iStr(rx)),rand=0.001)
                    continue
                self.progLog('\r#DEPTH','Parsing %s depth: %s pos; %s read bases.' % (filename,rje.iStr(px),rje.iStr(rx)),rand=0.001)
                rx += int(data[3]); px += 1
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
            if ridout:
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
            slocus = None   # Current SNPTable locus under consideration for SNPFreq skipping
            sposlist = []   # List of positions for current SNPTable locus for SNPFreq skipping
            slx = 0
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in PILEUP:
                prevrid = ridlist[0:]
                ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #self.bugPrint(rje.sortKeys(rstart))
                #self.bugPrint(ridlist)
                #self.bugPrint(line)
                entry = self.parsePileupLine(line,ri,ridlist,rdel,qc)
                #self.bugPrint(entry)
                #self.bugPrint(rje.sortKeys(rstart))
                #self.debug(ridlist)
                if not entry:
                    if slocus:
                        self.printLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),log=False)
                        if prevrid:
                            self.warnLog('Reached end of locus %s but %s unfinished reads!' % (slocus,rje.iLen(prevrid)))
                        if slocus in self.list['SkipLoci']: self.printLog('\r#PARSE','Skipped %s: %s %s SNPs.      ' % (slocus,rje.iStr(slx),qtxt))
                        else: self.printLog('\r#PARSE','Parsed %s: %s %s SNPs.      ' % (slocus,rje.iStr(slx),qtxt))
                    break
                ## ~ [2b] Filter non-SNP positions for SNPFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                freqskip = entry['Locus'] in self.list['SkipLoci']    # Whether to skip this position for QXX output.
                if self.getBool('SNPFreq'):
                    #i# Match entry to SNPTable if snpfreq analysis
                    if entry['Locus'] != slocus:
                        if slocus:
                            self.printLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),log=False)
                            if prevrid:
                                self.warnLog('Reached end of locus %s but %s unfinished reads!' % (slocus,rje.iLen(prevrid)))
                            if slocus in self.list['SkipLoci']: self.printLog('\r#PARSE','Skipped %s: %s %s SNPs.      ' % (slocus,rje.iStr(slx),qtxt))
                            else: self.printLog('\r#PARSE','Parsed %s: %s %s SNPs.      ' % (slocus,rje.iStr(slx),qtxt))
                        slocus = entry['Locus']; slx = 0
                        if alt: sposlist = snpdb.indexDataList('AltLocus',slocus,'AltPos',sortunique=False)
                        else: sposlist = snpdb.indexDataList('Locus',slocus,'Pos',sortunique=False)
                    freqskip = entry['Pos'] not in sposlist
                ## ~ [2c] Update RIDList and output if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ridlist: ri = max(ri,max(ridlist))
                else: self.debug(entry)
                self.progLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),rand=0.01); px += 1
                entry['Dep'] = rje.sf(entry['N']/meandepth)
                for rid in entry.pop('Start'): rstart[rid] = entry['Pos']
                for rid in entry.pop('End'):
                    try:
                        startpos = rstart.pop(rid)
                        if ridout:
                            rje.writeDelimit(RIDOUT,outlist=[rid,entry['Locus'],startpos,entry['Pos']],delimit='\t')
                    except:
                        self.errorLog('RID Ending problem! Started: %s; Ending: %s' % (rje.sortKeys(rstart),rid))
                        self.debug(entry)
                ## ~ [2d] Filter non-SNP positions for SNPFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if freqskip: continue
                ## ~ [2e] Filter low quality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Remove (from back) any reads than do not meet QV cutoff
                entry['QN'] = len(entry['Reads'])
                if entry['QN'] < self.getInt('MinQN'): continue     # Not enough data for output
                ## ~ [2f] Alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.pileupAlleles(entry): continue
                #table.addEntry(entry)
                outlist = []
                for field in outfields: outlist.append(entry[field])
                rje.writeDelimit(PILEOUT,outlist,delimit='\t'); ex += 1; slx += 1
            self.printLog('\r#PARSE','Parsed %s: %s entries from %s lines: %s reads.' % (filename,rje.iStr(ex),rje.iStr(px),rje.iStr(ri)))
            PILEOUT.close()
            PILEUP.close()
            if ridout:
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
            return outfile
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
            #i# SNPOnly=T filtering should NOT be applied if performing snpfreq analysis
            snponly = self.getBool('SNPOnly') and not self.getBool('SNPFreq')
            self.bugPrint('%s: %s' % (entry,snponly))
            acut = 0
            if snponly:             # Parsing allele mincut values for comparison; quick bypass of non-SNP sites.
                aset = set(qread.values())
                self.bugPrint('%s' % aset)
                if len(aset) < 1: return False  # No reads?!
                elif len(aset) < 2:
                    if self.getBool('IgnoreRef') or entry['Ref'] in aset:
                        return False   # All reads identical
                elif len(aset) == 2:
                    if not self.getBool('Indels') and '-' in aset: return False   # All non-gap reads identical
                    elif self.getBool('IgnoreN') and 'N' in aset: return False   # All non-N reads identical
                elif len(aset) == 3:
                    if not self.getBool('Indels') and '-' in aset and self.getBool('IgnoreN') and 'N' in aset: return False   # All non-gap/N reads identical
                acut = self.getNum('MinCut')    # Minimum allele count allowed
                if acut < 1: acut = max(self.getInt('AbsMinCut'),len(qread)*acut)
            elif self.getBool('SNPFreq'):
                for allele in [entry['Ref']] + entry['ALT']:
                    alleles[allele] = 0
                    allrid[allele] = []
            ### ~ [1] Parse Alleles (and corresponding Read IDs) from entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for rid in rje.sortKeys(qread):
                read = qread[rid]
                if read in alleles: alleles[read] += 1; allrid[read].append(str(rid))
                else: alleles[read] = 1; allrid[read] = [str(rid)]
            asort = []  # List of tuples for sorting
            if self.getBool('SNPFreq'):
                for allele in entry['ALT']: asort.append((alleles[allele],allele))
                asort.sort()
                asort.reverse()     # This should now be in order of allele frequency
                allele = entry['Ref']
                asort.append((alleles[allele],allele))
            else:
                for allele in alleles: asort.append((alleles[allele],allele))
                asort.sort()
                asort.reverse()     # This should now be in order of allele frequency
            acheckx = 0         # Count check
            aseq = []           # List of (filtered?) alleles and counts
            akept = []          # List of kept alleles
            self.bugPrint(asort)
            for i in range(len(asort)):
                acheckx += asort[i][0]
                if asort[i][1] == 'N' and snponly and self.getBool('IgnoreN'): continue
                # asort = [(allele count, allele)]
                # asort[i][0] = count
                # asort[i][1] = allele: '-' = deletion; len > 1 = insertion
                if (asort[i][1] == '-' or len(asort[i][1]) > 1) and self.getBool('SNPOnly') and not self.getBool('SNPFreq') and not self.getBool('Indels'): continue
                if asort[i][0] >= acut:
                    aseq.append('%s:%d' % (asort[i][1],asort[i][0]))
                    akept.append(asort[i][1])
            if not aseq: return False   # No alleles kept
            if snponly:
                if len(aseq) == 1:
                    if self.getBool('IgnoreRef') or aseq[0][0] == entry['Ref']: return False         # Not a SNP
                elif self.getBool('Biallelic') and len(aseq) > 2: return False     # Not a biallelic SNP
            entry['Seq'] = string.join(aseq,'|')
            entry['RID'] = []
            akept.sort()
            for allele in akept:    #rje.sortKeys(allrid):
                entry['RID'].append('%s:%s' % (allele,string.join(allrid[allele],',')))
            entry['RID'] = string.join(entry['RID'],'|')
            if acheckx != entry['QN']:  #!# Convert these to test statements?
                self.errorLog('Allele versus Quality count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
            self.deBug(entry)
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
            #i# Skipping loci in skiploci list: no need to parse reads as should be no reads spanning loci
            if data[0] in self.list['SkipLoci']: return entry
            #i# NOTE: entry['Pos'] can be used with self.db('SNP') to bypass irrelevant positions for SNPFreq
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
                    elif rseq[:1] == '*':   # Check for existing deletions - this should be the missing bases
                        if rid not in rdel:
                            self.warnLog('Deletion sequence confusion @ %s:%s (RID:%s not in rdel)' % (data[0],data[1],rid),warntype='del_confusion',quitchoice=True,suppress=True,dev=True)
                        elif ri > -1:
                            if rdel[rid][:1] != entry['Ref']: self.warnLog('Deletion sequence mismatch @ %s:%s (%s vs %s)' % (data[0],data[1],rdel[rid][:1],entry['Ref']),warntype='del_mismatch',quitchoice=True,suppress=True,dev=True)
                            rdel[rid] = rdel[rid][1:]
                            if not rdel[rid]: rdel.pop(rid) # End of deletion
                        reads.append('-'); rseq = rseq[1:]
                    elif rseq[:1] == '^':
                        rseq = rseq[2:]   # Indicates a new read; skip the quality score
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
            #self.bugPrint('%s:%s : %s'% (entry['Locus'],entry['Pos'],ridlist))
            ## ~ [1a] Check data integrity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ri > -1 and len(reads) != len(ridlist):
                raise ValueError('Read allele versus Read ID count mismatch for %s Pos %s: %s vs %s' % (entry['Locus'],entry['Pos'],reads,ridlist))
            if len(reads) != entry['N']:    # Formerly erroenous: (entry['N'] + delx):
                self.deBug('%s = %d' % (data[4],entry['N']))
                self.deBug('%s = %d' % (reads,len(reads)))
                self.errorLog('Read versus Read Count mismatch for %s Pos %s: %s vs %s' % (entry['Locus'],entry['Pos'],reads,entry['N']),printerror=False)
                raise ValueError
            ### ~ [2] Assess Quality Scores and generate alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qread = {}    # Dictionary of {rid:allele} passing qscore
            qual = data[5]
            if len(reads) != len(qual):
                self.deBug('%s = %d' % (reads,len(reads)))
                self.deBug('%s = %d' % (qual,len(qual)))
                self.deBug(data)
                self.errorLog('Read versus Quality length mismatch for %s Pos %s: %s vs %s' % (entry['Locus'],entry['Pos'],reads,qual),printerror=False)
                raise ValueError
            for r in range(len(qual)):
                qrid = ridlist[r]
                q = ord(qual[r]) - 33
                if qc:
                    qc += [0] * (q - len(qc)); qc[q-1] += 1
                if q >= self.getInt('QCut'): qread[qrid] = reads[r]
            while ridend: ridlist.remove(ridend.pop(0))     # Remove ended reads
            entry['Reads'] = qread
            return entry
        except:
            self.errorLog('%s.parsePileupLine() error' % (self.prog()));
            raise
#########################################################################################################################
    def forkerParsePileUp(self,depth=True):  ### Extracts, filters and processes PileUp data
        '''
        Extracts, filters and processes PileUp data.
        >> depth:bool [True] = Whether to generate read depth calculation.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Forking Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setBool({'Child':False})   # Whether this is a child process
            forkx = self.getInt('Forks')    # Number of forks to have running at one time
            forks = []                      # List of active fork PIDs
            tofork = []                     # List of pileup forking dictionaries
            #i# Each locus for each file will be forked off separately. The first pass that calculates read depths will
            #i# also generate these forking dictionaries:
            #i# - File: pileup file name
            #i# - Alt: Whether this is a Ref (False) or Alt (True) file
            #i# - Locus: Locus to be parsed in fork
            #i# - FPos: Position in file that starts locus
            #i# - Dep: Mean depth of locus
            #i# - RIDStart: RID for first locus read
            #i# - PID: PID of fork once started
            forked = {}                             # Dictionary of PID:fork dictionary
            killforks = self.getInt('KillForks')    # Time in seconds to wait after main thread has apparently finished
            killtime = time.time()                  # Baseline time for killtime assessment
            ## ~ [1b] General settings update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outbase = {}    # Dictionary of {filename:outbase}
            outfile = {}    # Dictionary of {filename:outfile}
            ridfile = {}    # Dictionary of {filename:ridfile} - None if ridout=False
            qtxt = 'Q%d.%d' % (self.getInt('QCut'),self.getInt('MinQN'))
            suffix = 'b%si%s' % (str(self.getBool('Biallelic'))[:1],str(self.getBool('Indels'))[:1])
            if self.getBool('SNPFreq'):
                snpdb = self.loadSNPTable()
                if self.getBool('RID'):
                    self.setBool({'RID':False})
                    self.printLog('#RID','RID allele output switched off for SNPFreq analysis.')
            ## ~ [1c] Output file update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfields = ['Locus','Pos','Ref','N','QN','Seq']    # Fields for output file
            if depth: outfields.append('Dep')
            if self.getBool('RID'):
                outfields.append('RID')
            for pfile in self.list['Batch']:
                outbase[pfile] = '%s.%s' % (rje.baseFile(pfile,strip_path=True),qtxt)
                if self.getBool('SNPFreq'): outfile[pfile] = '%s.snpfreq.tdt' % (outbase[pfile])
                elif self.getBool('SNPOnly'): outfile[pfile] = '%s.snponly.%s.tdt' % (outbase[pfile],suffix)
                else: outfile[pfile] = '%s.tdt' % (outbase[pfile])

            ### ~ [2] Setup the forking dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ Setup Forked Parsing of PileUp Files ##')
            for filename in self.list['Batch']:
                ridfile[filename] = '%s.rid.tdt' % (rje.baseFile(filename,strip_path=True))
                alt = False     # Whether this is an Alt pileup file (else Ref)
                for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                    if self.getStr(fstr) == filename:
                        self.setStr({fstr:outfile[filename]})     # Also update the filename for later parsing
                        if fstr.startswith('Alt'): alt = True
                ## ~ [2a] Assess whether the parsing needs to be performed ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Independently log whether to remake RID file: should be the same for all runs
                skiprun = not self.force()
                ridout = self.getBool('RID') or not rje.exists(ridfile[filename])  # Always output RID file for safety if RID mapping used
                ridout = ridout and not self.getBool('SNPFreq')
                if not ridout: ridfile[filename] = None
                skiprun = skiprun and not ridout    #(rje.exists(ridfile) or self.getBool('SNPFreq'))
                skiprun = skiprun and rje.exists(outfile[filename])
                if skiprun:
                    #!# Should we add something to check integrity/completeness? #!#
                    self.printLog('#SKIP','%s found! (force=F)' % outfile[filename])
                    continue
                ## ~ [2b] Parse depth and generate forking dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Currently this is based on unfiltered 'N' values.
                PILEUP = open(filename,'r'); px = 0; rx = 0; bx = 0
                prex = 0        # Number of reads preceding current locus
                prep = 0        # Number of positions preceding current locus
                locus = None    # Previous locus
                fpos = 0        # File position at start of locus
                linepos = 0     # This will be the position of the start of the current line
                line = 'Parsing'# Contents of last line
                while line:
                    line = PILEUP.readline()
                    data = string.split(rje.chomp(line))
                    if locus and (not data or data[0] != locus):    # Forking entry for previous locus
                        forkdict = {'File':filename,'Locus':locus,'Alt':alt,'FPos':fpos,'RIDStart':prex+1,'Pos':px-prep,
                                    'Out':'%s.%s.tdt' % (outbase[filename],locus),
                                    'RID':'%s.%s.rid.tdt' % (outbase[filename],locus),
                                    'Size':linepos-fpos,'Reads':rx-prex}
                        if not ridfile[filename]: forkdict['RID'] = None
                        #!# Could use Size or Reads to sort from big to small
                        if locus not in self.list['SkipLoci']:
                            tofork.append(forkdict)
                            self.printLog('\r#PARSE','Queued %s %s for forking: %s pos; %s reads.' % (filename,locus,rje.iStr(forkdict['Pos']),rje.iStr(forkdict['Reads'])))
                        else:
                            self.printLog('\r#DEPTH','Skipped %s %s depth (%s pos; %s reads parsed).' % (filename,locus,rje.iStr(forkdict['Pos']),rje.iStr(forkdict['Reads'])))
                            if self.dev(): self.printLog('#DEV','%s' % forkdict)
                        prex = rx; prep = px
                        fpos = linepos
                    if not data: break
                    locus = data[0]
                    linepos = PILEUP.tell()
                    if data[0] in self.list['SkipLoci']:
                        self.progLog('\r#DEPTH','Skipping %s depth: %s pos; %s reads; %s bases' % (filename,rje.iStr(px),rje.iStr(rx),rje.iStr(bx)),rand=0.001)
                        continue
                    self.progLog('\r#DEPTH','Parsing %s depth: %s pos; %s reads; %s bases.' % (filename,rje.iStr(px),rje.iStr(rx),rje.iStr(bx)),rand=0.001)
                    bx += int(data[3]); px += 1
                    rx += data[4].count('^')
                meandepth = float(bx) / px
                self.printLog('\r#DEPTH','Parsed %s depth: %s pos; %s reads; %s bases -> mean depth = %.1f' % (filename,rje.iStr(px),rje.iStr(rx),rje.iStr(bx),meandepth))
                PILEUP.close()
                ## ~ [2c] Fill in mean depth values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for forkdict in tofork:
                    if 'Dep' not in forkdict: forkdict['Dep'] = meandepth
                ## ~ [2d] Setup files for main processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rje.backup(self,outfile[filename])
                PILEOUT = open(outfile[filename],'w')
                rje.writeDelimit(PILEOUT,outlist=outfields,delimit='\t')
                PILEOUT.close()
                if ridout:
                    RIDOUT = open(ridfile[filename],'w')
                    rje.writeDelimit(RIDOUT,outlist=['RID','Locus','Start','End'],delimit='\t')
                    RIDOUT.close()
                else: ridfile[filename] = None
            ## ~ [2d] Sort by size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sortfork = []
            for fork in tofork: sortfork.append((fork['Size'],fork))
            sortfork.sort()
            sortfork.reverse()  #i# Sort the list of forks from big to small: run largest first!
            tofork = []
            for (size,fork) in sortfork: tofork.append(fork)
            self.printLog('#FORK','%d loci to fork for pileup processing; sorted by size (big to small).' % len(tofork))

            ### ~ [3] Fork through pileup parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            preforktime = time.time()   # Time at start of forking: used to estimate max time remaining.
            estend = preforktime        # Estimated maximum end time
            complete = 0                # Number of complete forks
            while tofork or forks:
                if complete:
                    #secs = float((len(tofork) + len(forks))) * (time.time() - preforktime) / float(complete)
                    secs = estend - time.time()
                    self.progLog('\r#FORK','%d of %d loci to fork; %d running (est %s)...' % (len(tofork),len(sortfork),len(forks),rje.hhmmss(secs)))
                else: self.progLog('\r#FORK','%d of %d loci to fork; %d running...' % (len(tofork),len(sortfork),len(forks)))
                ## ~ [3a] Set next fork going ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while tofork and len(forks) < forkx:
                    killtime = time.time()      # Reset killtime - still doing stuff
                    nextfork = tofork.pop(0)     # Forking dictionary to process
                    newpid = os.fork()
                    if newpid == 0: # child
                        self.setBool({'Child':True})
                        self.parsePileUpFork(nextfork,outfields)
                        sys.exit()    # Exit process
                    elif newpid == -1: self.errorLog('Problem forking %s:%s.' % (rje.baseFile(nextfork['File']),nextfork['Locus']),printerror=False)  # Error!
                    else:
                        self.printLog('\r#FORK','Forking %s:%s -> PID:%s' % (rje.baseFile(nextfork['File']),nextfork['Locus'],newpid))
                        nextfork['PID'] = newpid
                        forked[newpid] = nextfork
                        forks.append(newpid)    # Add fork to list

                ## ~ [3b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(1)       # Sleep for 1s
                forklist = self._activeForks(forks)
                if len(forklist) != len(forks):
                    self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    forks = forklist[0:]        # Reduce the list of active fork PIDs
                    for pid in forked.keys():   # Go through current forks
                        if pid not in forks:
                            fork = forked.pop(pid)
                            self.printLog('\r#FORK','Fork %s:%s finished.' % (rje.baseFile(fork['File']),fork['Locus']))
                            filename = fork['File']
                            if ridfile[filename]:
                                rje.fileTransfer(fromfile=fork['RID'],tofile=ridfile[filename],deletefrom=not self.dev(),append=True)
                                self.printLog('#RID','%s -> %s transfer complete.' % (fork['RID'],ridfile[filename]))
                            rje.fileTransfer(fromfile=fork['Out'],tofile=outfile[filename],deletefrom=not self.dev(),append=True)
                            self.printLog('#OUT','%s -> %s transfer complete.' % (fork['Out'],outfile[filename]))
                            killtime = time.time()  # Reset killtime - still doing stuff
                    forking = len(tofork) + len(forks)
                    complete = len(sortfork) - forking
                    estend = time.time() + (float(forking) * (time.time() - preforktime) / float(complete))

                ## ~ [3c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if (time.time() - killtime) > killforks:
                    self.warnLog('%d seconds of main thread inactivity. %d forks still active!' % (killforks,len(forks)))
                    for fork in forks:
                        self.warnLog('Fork %s, PID %d still Active!' % (forked[fork],fork))
                    if self.i() < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                    elif rje.yesNo('Kill hanging forks?'):
                        for fork in forks:
                            self.printLog('\r#KILL','Killing Fork %s, PID %d.' % (forked[fork],fork))
                            os.system('kill %d' % fork)
                    else: killtime = time.time()    # Rest killtime

            ### ~ [4] End of forker tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#FORK','Forked pileup parsing complete.')
            self.printLog('#QC','NOTE: No QC summary output for forked parsing.')
            #!# In future, combine QC scores #!#
            return outfile

        except SystemExit:
            if self.getBool('Child'): sys.exit(0)
            else:
                self.errorLog('Error in forkerParsePileUp()')
                raise ValueError('Parent forker raised SystemExit!')
        except:
            if self.getBool('Child'): sys.exit(1)
            self.errorLog('Error in forkerParsePileUp()'); raise
#########################################################################################################################
    def parsePileUpFork(self,fork,outfields): ### Extracts, filters and processes PileUp data for forked locus
        '''
        Extracts, filters and processes PileUp data for forked locus.
        >> fork:dict = Forked dictionary.
        >> outfields:list = Main output fields.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Each locus for each file will be forked off separately. The first pass that calculates read depths will
            #i# also generate these forking dictionaries:
            #i# - File: pileup file name
            #i# - Alt: Whether this is a Ref (False) or Alt (True) file
            #i# - Locus: Locus to be parsed in fork
            #i# - FPos: Position in file that starts locus
            #i# - Dep: Mean depth of locus
            #i# - RIDStart: RID for first locus read
            #i# - PID: PID of fork once started
            #i# - Out: output file name (temporary)
            #i# - RID: RID output file name (temporary) - None if no RID output
            snpdb = self.db('SNP')
            ## ~ [1a] Setup files for main processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            RIDOUT = None
            if fork['RID']: RIDOUT = open(fork['RID'],'w')
            PILEOUT = open(fork['Out'],'w')
            PILEUP = open(fork['File'],'r')
            PILEUP.seek(fork['FPos'])
            qtxt = 'Q%d.%d' % (self.getInt('QCut'),self.getInt('MinQN'))
            ## ~ [1b] Setup counters and lists for processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            px = 0          # Position counter for parsed file
            ex = 0          # Output entry counter for parsed file
            qc = [0]        # List of counts at different quality scores
            ri = fork['RIDStart']   # Read counter (ID counter)
            ridlist = []    # List of read IDs
            rdel = {}       # Dictionary of {rid:current deletion}
            rstart = {}     # Dictionary of {rid:start point}
            slocus = fork['Locus']  # Current SNPTable locus under consideration for SNPFreq skipping
            sposlist = []           # List of positions for current SNPTable locus for SNPFreq skipping
            if snpdb:
                if fork['Alt']: sposlist = snpdb.indexDataList('AltLocus',slocus,'AltPos',sortunique=False)
                else: sposlist = snpdb.indexDataList('Locus',slocus,'Pos',sortunique=False)
                if not sposlist: self.printLog('#SNP','No loaded SNPs for %s:%s!' % (fork['File'],fork['Locus']))
                altalleles = {}
                for entry in snpdb.indexEntries('Locus',slocus):
                    ekey = (entry['Locus'],entry['Pos'])
                    if ekey not in altalleles: altalleles[ekey] = []
                    altalleles[ekey].append(entry['ALT'])
            else: sposlist = '*'
            slx = 0         # SNP counter. Should match ex for forked locus!
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while sposlist:
                ## ~ [2a] Check next locus and end if move on ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                line = PILEUP.readline()
                if line: nextlocus = string.split(line)[0]
                else: nextlocus = None
                if not nextlocus or nextlocus != slocus:   # Reached end of locus
                    if ridlist: self.warnLog('Reached end of locus %s but %s unfinished reads!' % (slocus,rje.iLen(ridlist)))
                    #self.printLog('\r#PARSE','Parsed %s:%s: %s %s positions.      ' % (fork['File'],slocus,rje.iStr(slx),qtxt))
                    break
                ## ~ [2b] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.bugPrint(rje.sortKeys(rstart))
                self.bugPrint(ridlist)
                self.bugPrint(line)
                entry = self.parsePileupLine(line,ri,ridlist,rdel,qc)
                self.bugPrint(entry)
                self.bugPrint(rje.sortKeys(rstart))
                #self.debug(ridlist)
                ## ~ [2c] Filter non-SNP positions for SNPFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                freqskip = False    # Whether to skip this position for QXX output.
                if self.getBool('SNPFreq'):
                    #i# Match entry to SNPTable if snpfreq analysis
                    freqskip = entry['Pos'] not in sposlist
                ## ~ [2d] Update RIDList and output if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ridlist: ri = max(ri,max(ridlist))
                else: self.debug(entry)
                entry['Dep'] = rje.sf(entry['N']/fork['Dep'])
                for rid in entry.pop('Start'): rstart[rid] = entry['Pos']
                for rid in entry.pop('End'):
                    try:
                        startpos = rstart.pop(rid)
                        if RIDOUT:
                            rje.writeDelimit(RIDOUT,outlist=[rid,entry['Locus'],startpos,entry['Pos']],delimit='\t')
                    except:
                        self.errorLog('RID Ending problem! Started: %s; Ending: %s' % (rje.sortKeys(rstart),rid))
                        self.debug(entry)
                ## ~ [2e] Filter non-SNP positions for SNPFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if freqskip: continue
                ## ~ [2f] Filter low quality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Remove (from back) any reads than do not meet QV cutoff
                entry['QN'] = len(entry['Reads'])
                if entry['QN'] < self.getInt('MinQN'): continue     # Not enough data for output
                ## ~ [2g] Alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry['ALT'] = []
                if snpdb:
                    ekey = (entry['Locus'],entry['Pos'])
                    entry['ALT'] = altalleles[ekey]
                if not self.pileupAlleles(entry): continue
                #table.addEntry(entry)
                outlist = []
                for field in outfields: outlist.append(entry[field])
                rje.writeDelimit(PILEOUT,outlist,delimit='\t'); ex += 1; slx += 1
            self.printLog('\r#PARSE','Parsed %s:%s: %s %s entries from %s lines.' % (fork['File'],slocus,rje.iStr(ex),qtxt,rje.iStr(px)))
            PILEOUT.close()
            PILEUP.close()
            if RIDOUT:
                self.printLog('#RID','%s read start/end positions output to %s' % (rje.iStr(ri),fork['RID']))
                RIDOUT.close()
            return True
        except: self.errorLog('Error in parsePileUpFork()'); raise
#########################################################################################################################
    ### <4> ### Pileup SNP stats methods                                                                                #
#########################################################################################################################
    def pileUpStats(self,snpdb=None,alt=False,locfmt='full'):  ### Calculates statistics of genetic differences from parsed PileUp Tables
        '''
        Calculates statistics of genetic differences from parsed PileUp Tables.
        >> snpdb:Table [None] = Loaded SNP Table. If found, will restrict analysis to SNPTable positions and alleles
        >> alt:Bool [False] = Whether to use AltLocus and AltPos as keys.
        >> locfmt:str ['full'] = Locus format for mapping Pileup-parsed loci onto SNPdb.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Print details to log ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if alt:
                self.printLog('#~~#','## ~~~~~ Parsing & Combining Alt PileUp Tables (%s format) ~~~~~ ##' % locfmt)
                self.printLog('#ALT','AltControl and AltTreatment')
            else:
                self.printLog('#~~#','## ~~~~~ Parsing & Combining PileUp Tables (%s format) ~~~~~ ##' % locfmt)
                self.printLog('#BASE','Output basefile: %s' % self.baseFile())
            ## ~ [0b] Set output file name and return/backup if found  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfile = '%s.snp.tdt' % self.baseFile()
            snpindex = '#Locus#|#Pos#'
            if snpdb:
                if alt:
                    outfile = '%s.altcomb.tdt' % self.baseFile()
                    snpindex = '#AltLocus#|#AltPos#'
                else:
                    outfile = '%s.snpcomb.tdt' % self.baseFile()
                    snpindex = '#Locus#|#Pos#'
                snpdb.index(snpindex,make=True)
            if not self.force() and os.path.exists(outfile):
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return outfile
            rje.backup(self,outfile)
            ## ~ [0c] Set output fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# These are relative to the reference used for the fasta, e.g. Ref for ALt analysis will be the Alt Genome
            outfields = ['Locus','Pos','Ref']
            for fstr in self.list['Labels']:
                for field in ['N','Dep','QN','AN','Seq']:
                    outfields.append('%s|%s' % (field,fstr))
            outfields += ['MajFreq','MajDiff','MajProb']
            ## ~ [0d] Open file handles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
            #i# NOTE: if snpdb then both files should have parsed (and contain) the same positions
            #i# If a position is found in only one file, this should indicate a lack of reads mapping in one of the
            #i# pileup files.
            locus = None
            cdata = self.readDelimit('Control'); cx = 1
            tdata = None; tx = 0
            bx = 0      # Positions covered by both Control and Treatment
            sx = 0      # SNPs covered by both
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
                # Check against SNPDB if given. Should be present.
                if snpdb:
                    acc = self.mapLocus(locus,locfmt)
                    pos = int(cdata[1])
                    ikey = '%s|%d' % (acc,pos)
                    #self.debug(ikey)
                    #?# Add warning for multiple alleles in Pileup data #?#
                    if alt: snpfound = ikey in snpdb.index(snpindex)
                    else: snpfound = ikey in snpdb.index(snpindex)
                    if not snpfound:
                        self.warnLog('%s found in Control and Treatment but not SNPTable. Changed SNPTable?' % ikey,warntype='ChangedSNPTable',quitchoice=True,suppress=True)
                        cdata = self.readDelimit('Control'); cx += 1; nosnpx += 1; continue
                # Convert to dictionaries
                # infields = ['Locus','Pos','Ref','N','QN','Seq','Dep']
                #self.bugPrint('C: %s' % cdata)
                #self.bugPrint('T: %s' % tdata)
                cdict = rje.list2dict(cdata,infields)
                tdict = rje.list2dict(tdata,infields)
                #self.bugPrint('C: %s' % cdict)
                #self.bugPrint('T: %s' % tdict)
                ## ~ [1b] Combine and assess alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                totx = 0.0      # Total allele count
                alleles = {}    # Dictionary of {allele:count}
                refseq = cdict['Ref']   # This is the "Reference" allele for this
                #i# If snpdb is given, alleles will be restricted to the REF and ALT alleles for this position.
                #i# These can be extracted from snpdb.index(snpindex).
                #i# Due to one-to-many mapping, there might be multiple alleles. Not all may have counts.
                if snpdb:
                    alleles[refseq] = 0
                    for allele in snpdb.indexDataList(snpindex,ikey,'REF') + snpdb.indexDataList(snpindex,ikey,'ALT'):
                        if not self.getBool('Indels') and (allele == '-' or len(allele) > 1): indx += 1; continue
                        elif self.getBool('IgnoreN') and allele == 'N': continue
                        alleles[allele] = 0
                    #i# Biallelic=T now restricts analysis to biallelic SNPs in SNPTable
                    #!# Added additional filter for multiple loci matching, also filtered by "biallelic=T" filter.
                    if self.getBool('Biallelic') and (len(alleles) != 2 or len(snpdb.index(snpindex)[ikey]) > 1):
                        notbix += 1
                        cdata = self.readDelimit('Control')
                        cx += 1
                        continue
                #i# If no SNPDB, build completely from data
                elif not self.getBool('IgnoreRef'): alleles[refseq] = 0
                for allele in string.split('%s|%s' % (cdict['Seq'],tdict['Seq']),'|'):
                    adata = string.split(allele,':')
                    if self.getBool('IgnoreN') and adata[0] == 'N': continue    # Skip N alleles
                    ax = int(adata[1])
                    if adata[0] not in alleles:
                        if snpdb: continue   #i# Only use snpdb alleles
                        alleles[adata[0]] = 0
                    totx += ax
                    alleles[adata[0]] += ax
                #self.debug(alleles)
                # Check min allele frequency
                #?# Will we want to handle this differently for snpdb #?#
                acut = self.getNum('MinCut')
                if acut < 1: acut = max(self.getInt('AbsMinCut'),totx*acut)
                for aseq in rje.sortKeys(alleles):
                    if aseq == refseq and not self.getBool('IgnoreRef'): continue   # Keep Reference allele
                    if alleles[aseq] < acut:    # Too low frequency: dump
                        cutx += 1
                        alleles.pop(aseq)
                if not alleles: raise ValueError('No alleles survived mincut!')
                # Check for SNP
                if snpdb:
                    if not totx:
                        self.warnLog('No reads with SNPTable or Ref Alleles map to %s' % ikey)
                        cdata = self.readDelimit('Control'); cx += 1; continue
                    if len(alleles) == 2 and refseq in alleles and not alleles[refseq]: difx += 1
                else:
                    if len(alleles) == 1:
                        if refseq not in alleles: difx += 1 # Fixed difference from reference
                        else: refx += 1
                        cdata = self.readDelimit('Control'); cx += 1
                        continue  # Not a SNP!
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
                if not ctot:
                    if snpdb: self.warnLog('No SNPTable allele read counts for Control %s' % ikey,warntype='Allele Zero Count',quitchoice=False,suppress=True)
                    cdata = self.readDelimit('Control'); cx += 1; continue
                cfreq /= ctot
                majx = int(tfreq)
                if not ttot:
                    if snpdb: self.warnLog('No SNPTable allele read counts for Treatment %s' % ikey,warntype='Allele Zero Count',quitchoice=False,suppress=True)
                    cdata = self.readDelimit('Control'); cx += 1; continue
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
                #self.bugPrint(entry)
                #self.debug('%s vs %s' % (entry['MajProb'] * fdrx,(totx - entry['FDR'] + 1)))
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
    def loadSNPTable(self,name='SNP',indels=True):     ### Loads and returns a SNP table for the genome(s) being analysed
        '''
        Loads and returns a SNP table for the genome(s) being analysed.
        >> indels:bool [False] = Whether to keep indels (True) or use self.getBool('Indels') (False)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.db(name): return self.db(name)
            if not self.getStrLC('SNPTable'): self.printLog('\r#SNP','No SNP table to add.'); return None
            snptable = self.getStr('SNPTable')
            snpkeys = self.list['SNPTabKeys']
            for field in ['Pos','Locus']:
                if field in snpkeys: snpkeys.remove(field)
                snpkeys.insert(0,field)
            self.printLog('#KEYS','Loading %s with keys: %s' % (snptable,string.join(snpkeys,'; ')))
            ### ~ [1] Load SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db().addTable(snptable,name=name,expect=True,mainkeys=snpkeys)
            ## ~ [1a] Reformat data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb.dataFormat({'Pos':'int','AltPos':'int'})
            snpdb.indexReport('Locus')
            if 'AltLocus' in snpdb.fields(): snpdb.indexReport('AltLocus')
            ## ~ [1b] Sort out fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'REF' not in snpdb.fields():
                if 'Ref' in snpdb.fields(): snpdb.renameField('Ref','REF')
                else: raise ValueError('SNP Table needs Ref or REF field.')
            if 'ALT' not in snpdb.fields():
                if 'Alt' in snpdb.fields(): snpdb.renameField('Alt','ALT')
                else: raise ValueError('SNP Table needs Alt or ALT field.')
            ## ~ [1c] Reduce data if no indels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if indels and not self.getBool('Indels'):
                indx = snpdb.entryNum()
                snpdb.dropEntriesDirect('REF',['-'])
                snpdb.dropEntriesDirect('ALT',['-'])
                indx -= snpdb.entryNum()
                self.printLog('#INDEL','%s Indels dropped (indels=F).' % rje.iStr(indx))
            return snpdb
        except ValueError: self.errorLog('%s.loadSNPTable() error' % (self)); raise
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
            #snpdb = self.db().addTable(snptable,name='SNP',expect=True,mainkeys=snpkeys)
            snpdb = self.loadSNPTable()
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
                fdb.makeField('#Locus#|#Pos#',warn=False)
                snpjoinfield = '#%s#|#%s#' % (snpjoin[0],snpjoin[1])
                snpdb.makeField(snpjoinfield,warn=False)
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
        except: self.errorLog('%s.combineSNPs() error' % (self)); return None
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
            snpdb = self.loadSNPTable('combined')   #X# ,indels=False) - Use existing SNP Table if loaded
            if not snpdb: raise IOError('Cannot perform AltControl/AltTreatment analysis without SNPTable.')
            #self.debug(snpdb.keys())
            mapfields = ['Locus','Pos','AltLocus','AltPos','#Locus#|#Pos#','#AltLocus#|#AltPos#']
            for field in mapfields[0:4]:
                if field not in snpdb.fields(): raise ValueError('Cannot perform AltControl/AltTreatment analysis without SNPTable field "%s".' % field)
            ## ~ [0b] Drop Indels: cannot match positions via short reads ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Default for SNPFreq analysis should be indels=F Not sure how/if it will work with indels=T!
            if not self.getBool('Indels'):
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
            #self.debug(snpdb.datakeys()[:10])
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
            elif not self.getBool('MajFocus') and self.getBool('MajRef'): self.printLog('#REFSNP','Generating output for all Ref ("Reference") alleles.')
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
            #self.debug(snpdb.keys())
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
    def coverageFromRID(self,ridfile=None,depthplot=False,readlen=False,clip=False,noclip=False):  ### Extracts read data from RID file and summarises read coverage
        '''
        Extracts read data from SAM file.
        >> ridfile:str [None] = RID file name to use. (Will use RID table if None)
        >> depthplot:bool [False] = Whether to calculate and output full read depth table.
        >> readlen:bool [False] = Whether to generate read length rather than read depth data.
        >> clip:bool[False] = Whether to restrict analysis to clipped reads (minsoftclip)
        >> noclip:bool[False] = Whether to restrict analysis to non-clipped reads (maxsoftclip) if clip=False
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            readdesc = 'read'
            if readlen: readdesc += ' length'
            else: readdesc += ' depth'
            if clip: readdesc = 'clipped ' + readdesc
            elif noclip: readdesc = 'noclip ' + readdesc
            self.printLog('#~~#','## ~~~~~ %s coverageFromRID ~~~~~ ##' % readdesc)
            db = self.db()
            #!# Add lengths to basefile and update prolog messages etc.
            if clip: db.setBasefile('%s.clipped' % self.baseFile())
            elif noclip: db.setBasefile('%s.noclip' % self.baseFile())
            else: db.setBasefile(self.baseFile())
            dirfile = None; dirstep = 0
            if readlen:
                covfile = '%s.readlen.tdt' % db.baseFile()
                depfile = '%s.readlenplot.tdt' % db.baseFile()
                if self.getInt('DirnLen') > 0:
                    dirstep = self.getInt('DirnLen')
                    dirfile = '%s.dirnlenplot.tdt' % db.baseFile()
            else:
                covfile = '%s.coverage.tdt' % db.baseFile()
                depfile = '%s.depthplot.tdt' % db.baseFile()
            makedep = self.force() or not rje.exists(covfile)
            if depthplot:
                makedep = makedep or not rje.exists(depfile)
                if not makedep:
                    if self.getStrLC('CheckPos'): self.warnLog('No new checkpos output: depthplot data found and force=F')
                    return self.printLog('#DEPTH','%s and %s found (force=F)' % (covfile,depfile))
            elif not makedep:
                if self.getStrLC('CheckPos'): self.warnLog('No new checkpos output: coverage data found and force=F')
                return self.printLog('#COV','%s found (force=F)' % (covfile))
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
            rdb.dataFormat({'Start':'int','End':'int','RLen':'int','MLen':'int','Clip5':'int','Clip3':'int'})
            if (clip or noclip) and 'Clip5' not in rdb.fields():
                self.printLog('#CLIP','Clip fields missing from RID table: cannot perform clip or noclip analysis')
                return False
            if dirfile: rdb.newKey(['Locus','Start','RID'])
            self.printLog('#~~#','## ~~~~~ Calculating %s coverage ~~~~~ ##' % readdesc)
            ## ~ [1b] Setup Check Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # This is self.getStr('CheckPos')
            if self.getStrLC('CheckPos') and not readlen:
                cdb = db.addTable(self.getStr('CheckPos'),mainkeys=self.list['CheckFields'],name='check',expect=True)
                cdb.dataFormat({startfield:'int',endfield:'int'})
                cdb.setStr({'Delimit':'\t'})
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
            if self.dev():
                rdb.indexReport('Locus')
                #self.debug('>>>')

            ### ~ [2] Calculate read coverage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fullcalc = depthplot and 'RLen' in rdb.fields() and 'MLen' in rdb.fields()
            #self.debug(self.getNum('FullCut'))
            dirdb = None    # DirnlenPlot data
            if readlen: covdb = db.addEmptyTable('readlen',['Locus','Length','MeanX'],['Locus'])
            else: covdb = db.addEmptyTable('coverage',['Locus','Length','MeanX'],['Locus'])
            rx = 0.0; rtot = rdb.entryNum(); lx = 0; ltot = len(rdb.index('Locus'))
            depth = {}    # Dictionary of locus:[depth per site]; Only used if depthplot=True
            depsmooth = self.getInt('DepthSmooth')  # Smooth out any read plateaus < depsmooth nucleotides in length
            maxsmooth = self.getNum('PeakSmooth')   # Max X coverage difference to smooth out for Xcoverage peaks
            depfields = ['Locus','Pos','X']
            if fullcalc:
                depfields += ['FullX','PartX']
            DEPFILE = None
            if depthplot:
                covdb.addFields(['MinX','MaxX','MedianX','Coverage'])
                rje.backup(self,depfile)
                DEPFILE = open(depfile,'w')
                DEPFILE.write('%s\n' % string.join(depfields,'\t'))
                if readlen and dirfile: dirdb = db.addEmptyTable('dirnlenplot',['Locus','Pos','Len5','Len3'],['Locus','Pos'])
            minfiltx = 0
            for locus in rdb.indexKeys('Locus'):
                lx += 1
                self.progLog('\r#COV','Calculating %s coverage: %.2f%% (Locus %d of %d)' % (readdesc,rx/rtot,lx,ltot),rand=0.01)
                centry = {'Locus':locus,'Length':0,'MeanX':0.0}
                if locus in seqdict: centry['Length'] = seqlist.seqLen(seqdict[locus])
                else: centry['Length'] = max(rdb.indexDataList('Locus',locus,'End',sortunique=False))
                #self.debug('%s: %s (%s)' % (locus,centry['Length'],type(centry['Length'])))
                if depthplot: depth[locus] = [0] * centry['Length']
                dirpos = [] # DirnLen assay points
                if dirdb:
                    di = 0; dirpos = [1]
                    # Each dirdb entry stores the longest 5' and 3' distances to ends of reads.
                    dirdb.addEntry({'Locus':locus,'Pos':1,'Len5':0,'Len3':0})
                    while di < centry['Length']:
                        di += dirstep; dirpos.append(di)
                        dirdb.addEntry({'Locus':locus,'Pos':di,'Len5':0,'Len3':0})
                    dirdb.addEntry({'Locus':locus,'Pos':centry['Length'],'Len5':0,'Len3':0})
                    dirpos.append(centry['Length'])
                    #self.debug('%s: %s' % (locus,dirpos))
                fullpos = [[],[]]   # [Start,End] of RID with MLen/RLen meeting the 'FullCut' threshold (1-L)
                partpos = [[],[]]   # [Start,End] of RID with MLen/RLen failing the 'FullCut' threshold (1-L)
                for rentry in rdb.indexEntries('Locus',locus):
                    self.progLog('\r#COV','Calculating %s coverage: %.2f%% (Locus %d of %d)' % (readdesc,rx/rtot,lx,ltot),rand=0.01); rx += 100
                    if 'RLen' in rdb.fields() and rentry['RLen'] < self.getInt('MinReadLen'): minfiltx += 1; continue
                    #!# Filter based on clip or noclip and 'Clip5' not in rdb.fields(): clip = noclip = False
                    if clip:
                        if max(rentry['Clip5'],rentry['Clip3']) < self.getInt('MinSoftClip'): continue
                        if self.getInt('MinSoftClip') < self.getInt('MaxSoftClip') and max(rentry['Clip5'],rentry['Clip3']) > self.getInt('MaxSoftClip'): continue
                    elif noclip:
                        if max(rentry['Clip5'],rentry['Clip3']) > self.getInt('MaxSoftClip'): continue
                    rlen = rentry['End'] - rentry['Start'] + 1
                    centry['MeanX'] += rlen
                    if fullcalc:
                        if rentry['RLen'] and float(rentry['MLen']) / rentry['RLen'] >= self.getNum('FullCut'):
                            fullpos[0].append(rentry['Start'])
                            fullpos[1].append(rentry['End'])
                        else:
                            partpos[0].append(rentry['Start'])
                            partpos[1].append(rentry['End'])
                    # Directional read length data
                    if dirdb:
                        #self.debug(rentry)
                        #x#if dirpos[0] > rentry['Start']: dirdb.addEntry({'Locus':locus,'Pos':rentry['Start'],'Len5':0,'Len3':rlen-1})
                        #x#if dirpos[0] > rentry['End']: dirdb.addEntry({'Locus':locus,'Pos':rentry['End'],'Len5':rlen-1,'Len3':0})
                        for pos in dirpos[0:]:
                            if pos < rentry['Start']: dirpos.remove(pos)
                            elif pos <= rentry['End']:
                                dentry = dirdb.data((locus,pos))
                                dentry['Len5'] = max(dentry['Len5'],pos-rentry['Start'])
                                dentry['Len3'] = max(dentry['Len3'],rentry['End']-pos)
                                #self.debug(dentry)
                            else: break
                        #self.debug('%s: %s' % (locus,dirpos))
                    # Position check data
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
                    #!# NOTE: This is not very efficient. Could be done much better by making depdata (pos,depth) list
                    #!# directly by walking through Start and End positions in order. Would make median calculation more
                    #!# of a challenge but not impossible. (Could regenerate depth list from depdata!)
                    #!# >> convert depdata into a list of [(depth,length)], sort, and walk through till sum(length) >= half centry['Length'].
                    #!# Job for another day!
                    #!# NOTE: Updating like this would make the readlen=TRUE calculations harder!
                    rlen = rje.dp(rlen/1000.0,1)    # Convert read length to kb (1 d.p.)
                    if depthplot:
                        for i in range(rentry['Start']-1,rentry['End']):
                            if readlen: depth[locus][i] = max(depth[locus][i],rlen)
                            else: depth[locus][i] += 1
                if readlen:
                    centry['MeanX'] = sum(depth[locus]) / centry['Length']
                else:
                    centry['MeanX'] /= centry['Length']
                if depthplot:
                    # Calculate Median
                    centry['MaxX'] = max(depth[locus])
                    centry['MinX'] = min(depth[locus])
                    centry['Coverage'] = len(depth[locus]) - depth[locus].count(0)
                    median = depth[locus][0:]
                    median.sort()
                    medlen = len(median)
                    #self.debug('%s >> %d' % (locus,medlen))
                    if medlen:
                        if rje.isOdd(medlen): centry['MedianX'] = median[medlen/2]
                        else: centry['MedianX'] = (median[medlen/2] + median[(medlen-1)/2]) / 2.0
                    else: centry['MedianX'] = 0.0
                    peakx = maxsmooth
                    if peakx < 1: peakx *= centry['MedianX']
                    # Compress depth[locus] and output boundaries
                    depdata = []    # (Pos,X) tuples
                    for i in range(len(depth[locus])):
                        if i in [0,len(depth[locus])-1]:
                            depdata.append((i+1,depth[locus][i]))
                        elif depth[locus][i] != depth[locus][i+1] or depth[locus][i] != depth[locus][i-1]:
                            depdata.append((i+1,depth[locus][i]))
                    # Reduce depdata to "peaks"
                    i = 1   # This is the position being considered for removal
                    while i < (len(depdata) - 1):
                        # Need to account for -/--/ vs -/--\ vs -/--------/ patterns
                        if depdata[i][1] == depdata[i+1][1] and depdata[i+1][0]-depdata[i][0] >= depsmooth:
                            i += 2; continue    # long enough plateau to keep both ends
                        elif depdata[i][1] == depdata[i+1][1]:
                            if i >= (len(depdata) - 2): break   # Reached end, so keep
                            # Only keep if a peak, otherwise remove both
                            up5 = depdata[i][1] > depdata[i-1][1]
                            up3 = depdata[i+1][1] > depdata[i+2][1]
                            if up5 == up3:  # Peak!
                                i += 2; continue
                            else: depdata.pop(i); depdata.pop(i)
                        else:   # Only keep if a peak, otherwise remove i
                            up5 = depdata[i][1] > depdata[i-1][1]
                            up3 = depdata[i][1] > depdata[i+1][1]
                            if up5 == up3:  # Peak!
                                i += 1; continue
                            else: depdata.pop(i)
                    # Next remove small deviations
                    for x in range(1,depsmooth):
                        # Remove any peaks of size x within peakx X coverage of flanking values
                        i = 1
                        while i < (len(depdata) - 1):
                            if i >= (len(depdata) - 2): break   # Reached end, so keep
                            if depdata[i][1] != depdata[i+1][1]: i += 1; continue       # Not a peak
                            if depdata[i+1][0] - depdata[i][0] != x: i += 1; continue   # Wrong size
                            if max(rje.modulus(depdata[i][1] - depdata[i-1][1]), rje.modulus(depdata[i+1][1] - depdata[i+2][1])) > peakx:
                                i += 1; continue
                            depdata.pop(i); depdata.pop(i)
                            while i < (len(depdata) - 1) and depdata[i-1][1] == depdata[i][1] == depdata[i+1][1]:
                                depdata.pop(i); depdata.pop(i-1)
                    # FullCut
                    if fullcalc:
                        fullpos[0].sort()
                        fullpos[1].sort()
                        fullx = 0
                        partpos[0].sort()
                        partpos[1].sort()
                        partx = 0
                        #self.debug(fullpos[0])
                        #self.debug(fullpos[1])
                        #self.debug(partpos[0])
                        #self.debug(partpos[1])
                    # Output
                    for (pos,x) in depdata:
                        posdata = [locus,str(pos),str(x)]
                        if fullcalc:
                            #self.bugPrint('%d <= %d: %s %s' % (fullpos[0],pos,fullpos[0] <= pos,fullpos[0] <= int(pos)))
                            while fullpos[0] and fullpos[0][0] <= pos: fullx += 1; fullpos[0].pop(0)
                            while fullpos[1] and fullpos[1][0] <= pos: fullx -= 1; fullpos[1].pop(0)
                            while partpos[0] and partpos[0][0] <= pos: partx += 1; partpos[0].pop(0)
                            while partpos[1] and partpos[1][0] <= pos: partx -= 1; partpos[1].pop(0)
                            posdata += [str(fullx),str(partx)]
                            #self.debug('%s: %s -> %s | %s -> %s' % (pos,fullpos[0][:10],fullpos[1][:10],partpos[0][:10],partpos[1][:10]))
                        DEPFILE.write('%s\n' % string.join(posdata,'\t'))
                covdb.addEntry(centry)
            self.printLog('\r#COV','Calculating %s coverage: complete (%d of %d loci)' % (readdesc,lx,ltot))
            if self.getInt('MinReadLen') > 0: self.printLog('#MINLEN','%s reads filtered (< %s bp)' % (rje.iStr(minfiltx),rje.iStr(self.getInt('MinReadLen'))))
            if depthplot:
                DEPFILE.close()
                self.printLog('\r#DEPTH','XCoverage %s plot data (depthsmooth=%d; peaksmooth=%.2f) output to %s.' % (readdesc,depsmooth,maxsmooth,depfile))
            if dirdb: dirdb.saveToFile()
            if cdb: cdb.saveToFile(sfdict={'MeanX':3})
            covdb.saveToFile(sfdict={'MeanX':3})
            if depthplot and self.getBool('RGraphics'):
                if readlen: self.rGraphics('samreadlen',basefile=db.baseFile())
                else: self.rGraphics('samdepth',basefile=db.baseFile())
            return True
        except: self.errorLog('%s.coverageFromRID() error' % (self.prog())); return False
#########################################################################################################################
    def rGraphics(self,rtype='samtools',rargs='',basefile=None):    ### Generates SAMTools R Graphics
        '''Generates SAMTools R Graphics.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# See pagsat.PAGSAT.rGraphics() for additional steps, such as file checking and adding options.
            #?# Do we want to generate a summary plot of coverage (e.g. histogram?)
            self.printLog('#~~#','## ~~~~~ SAMTools R Graphics ~~~~~ ##')
            if not basefile: basefile = self.baseFile()
            rcmd = '%s --no-restore --no-save --args "%s" "%s"' % (self.getStr('RPath'),rtype,basefile)
            rdir = '%slibraries/r/' % slimsuitepath
            rtmp = '%s.r.tmp.txt' % basefile
            if rargs: rcmd += ' %s' % rargs
            rcmd += ' "rdir=%s" < "%srje.r" > "%s"' % (rdir,rdir,rtmp)
            self.printLog('#RPNG',rcmd)
            problems = os.popen(rcmd).read()
            if problems:
                for ptxt in problems: self.warnLog(ptxt)
            elif rje.exists(rtmp) and not self.dev() and not self.debugging(): os.unlink(rtmp)
            if not problems: self.printLog('#RPNG','R Graphics generated.')
            return not problems
        except: self.errorLog('%s.rGraphics error' % self.prog()); return False
#########################################################################################################################
    def rScriptGraphics(self,rtype='snpfreq',rargs=''):    ### Generates SAMTools R Graphics
        '''Generates SAMTools R Graphics.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ SAMTools R Graphics ~~~~~ ##')
            basename = self.baseFile(strip_path=True)
            rdir = '%slibraries/r/' % slimsuitepath
            rcmd = 'Rscript "%srje.r" "%s" "%s"' % (rdir,rtype,basename)
            rtmp = '%s.r.tmp.txt' % self.baseFile()
            if rtype == 'snpfreq':
                rcmd += ' freqpath=./ freqcomp=%s multiplot=T' % basename
                snapbase = self.getStr('SNPTable') # *.snpmap.tdt
                snapbase = string.join(string.split(snapbase,'.')[:-2],'.')
                rcmd += ' snapbase=%s' % snapbase
            if rargs: rcmd += ' %s' % rargs
            rcmd += ' "rdir=%s"' % rdir
            #rcmd += ' "rdir=%s" > "%s"' % (rdir,rtmp)
            self.printLog('#RPNG',rcmd)
            #!# Improve handling of R output and progress! #!#
            problems = False
            if self.dev():
                problems = os.popen(rcmd).read()
                if problems:
                    for ptxt in problems: self.warnLog(ptxt)
                elif rje.exists(rtmp) and not self.dev() and not self.debugging(): os.unlink(rtmp)
            else:
                RCMD = os.popen(rcmd)
                rline = True
                while rline:
                    rline = RCMD.readline()
                    if rline: self.vPrint(rje.chomp(rline))
                RCMD.close()
            if not problems: self.printLog('#RPNG','R Graphics generated.')
            return not problems
        except: self.errorLog('%s.rGraphics error' % self.prog()); return False
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
