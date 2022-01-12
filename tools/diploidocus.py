#!/usr/bin/python

# See below for name and description
# Copyright (C) 2017 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       Diploidocus
Description:  Diploid genome assembly analysis toolkit
Version:      1.1.1
Last Edit:    06/01/22
Nala Citation:  Edwards RJ et al. (2021), BMC Genomics [PMID: 33726677]
DipNR Citation: Stuart KC, Edwards RJ et al. (preprint), bioRxiv 2021.04.07.438753; [doi: 10.1101/2021.04.07.438753]
Tidy Citation:  Chen SH et al. & Edwards RJ (preprint), bioRxiv 2021.06.02.444084; [doi: 10.1101/2021.06.02.444084]
GitHub:       https://github.com/slimsuite/diploidocus
Copyright (C) 2020  Richard J. Edwards - See source code for GNU License Notice

Function:
    Diploidocus is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
    The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
    BUSCO gene prediction and contaminant screening for a number of assembly tasks including contamination identification,
    haplotig identification/removal and low quality contig/scaffold trimming/filtering.

    In addition, Diploidocus will use mapped long reads and BUSCO single copy read depths for genome size prediction
    (`gensize`), and coverage (`regcheck`) or copy number estimation (`regcnv`) for user-defined regions. Diploidocus
    also has functions for removing redundancy (`sortnr`), generating a non-redundant pseudo-diploid assembly with primary
    and secondary scaffolds from 10x pseudohap output (`diphap`), and creating an in silico diploid set of long reads from two
    haploid parents (for testing phasing etc.) (`insilico`).

    Please note that Diploidocus is still in development and documentation is currently a bit sparse. Like its (almost)
    namesake, Diploidocus has developed into a large, lumbering beast and is in the process of being split up into a
    number of specialised tools. Core mechanics will still be performed by Diploidocus but some of the documentation will
    be moved to other programs and sites.

    The different run modes are set using `runmode=X`:

    * `diploidocus` default run mode will run `gensize`, `telomeres`, `vecscreen`, `deptrim` and `purgehap` analysis
    * `dipcycle` runs iterative cycles of `diploidocus` mode until convergence (no more filtering) is reached
    * `gensize` uses BUSCO results, a BAM file and read file(s) to predict the genome size of the organism
    * `purgehap` filters scaffolds based on post-processing of purge_haplotigs
    * `telomeres` performs a regex telomere search based on method of https://github.com/JanaSperschneider/FindTelomeres
    * `vecscreen` searches for contaminants and flags/masks/removes identified scaffolds
    * `deptrim` trims sequence termini of at least `mintrim=INT` bp with less than `deptrim=INT` read depth
    * `regcheck` checks reads spanning given regions and also calculates mean depth and estimated copy number (if regcnv=T)
    * `regcnv` calculates mean depth and estimated copy number for regcheck regions from BAM file (no read spanning analysis)
    * `gapspan` is a specialised `regcheck` analysis that checks for reads spanning genome assembly gaps
    * `gapass` extends the `gapspan` mode to extract the reads spanning a gap and re-assemble them using `flye`
    * `gapfill` extends the `gapass` mode to map re-assembled gaps back onto the assembly and replace gaps spanned by the new contigs
    * `sortnr` performs an all-by-all mapping with minimap2 and then removes redundancy
    * `diphap` splits a pseudodiploid assembly into primary and alternative scaffolds
    * `diphapnr` runs `sortnr` followed by `diphap`
    * `insilico` generates balanced diploid combined reads from two sequenced haploid parents

    See <https://slimsuite.github.io/diploidocus/> for details of each mode. General SLiMSuite run documentation can be
    found at <https://github.com/slimsuite/SLiMSuite>.

    Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
    <https://github.com/slimsuite/diploidocus>.

Run Modes:
    ### ~ Main Diploidocus filtering [runmode=diploidocus] ~ ###
    Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
    BUSCO results to reclassify scaffolds and partition them into:

    * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
    * `*.core.fasta` = the same set of scaffolds, minus repeats
    * `*.repeats.fasta` = the repeat scaffolds excluded from core
    * `*.junk.fasta` = low coverage scaffolds, removed as junk
    * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

    **NOTE:** PurgeHaplotigs was not using zero coverage bases in its percentages. This is now fixed by Diploidocus.

    _See main docs for details._

    ---
    ### ~ Cycled Diploidocus filtering [runmode=dipcycle] ~ ###

    Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
    no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
    using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
    scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

    Output for each cycle will be initially generated in the run directory but then moved to a `dipcycle_$BASEFILE`
    directory upon completion.

    Final outputs from the final cycle will then be compiled under the original `$BASEFILE` prefix:

    * `$BASEFILE.diploidocus.tdt` = Final ratings for the input scaffolds. This is like the single cycle output with an additional `Cycle` field containing the last cycle this scaffold was processed in.
    * `$BASEFILE.ratings.tdt` = Reduced final ratings output for the input scaffolds (`SeqName`,`SeqLen`,`ScreenPerc`,`Class`,`Rating`,`Cycle`).
    * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
    * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
    * `$BASEFILE.repeats.fasta` = the repeat scaffolds excluded from core
    * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
    * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

    **NOTE:** Contents for these four `*.fasta` files are summarised in the main log. Individual purge cycles have their own
    log files in the `dipcycle_$BASEFILE` directory.

    See `runmode=diploidocus` documentation for more details.

    ---
    ### ~ Genome size prediction [runmode=gensize] ~ ###

    This works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid read
    depth) regions along with some ooor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
    are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
    region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
    read depth.

    Diploidocus uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
    and extracts the modal read depth for each BUSCO gene along with the overall modal read depth for all gene
    regions. Genome size is then estimated based on a crude calculation using the total combined sequencing length.
    This will be caculated from `reads=FILELIST` unless provided with `readbp=INT`.

    BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
    `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
    ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
    and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.

    **NOTE:** The current genome size prediction appears to be an over-estimate. There is currently no adjustment for
    contamination. The `mapadjust` option attemtps to correct for read mapping and imbalanced insertion:deletion ratios
    etc. but has not been extensively tested.

    ---
    ### ~ Running Purge_haplotigs using BUSCO-guided cutoffs [runmode=purgehap] ~ ###
    _See main docs_

    ---
    ### ~ Telomere finding [runmode=telomere] ~ ###
    _Details coming soon!_

    ---
    ### ~ Vector/contamination screening [runmode=vecscreen] ~ ###
    _See main docs_

    ---
    ### ~ Depth trimming [runmode=deptrim] ~ ###

    Depth trimming (`deptrim`) mode trims sequence termini of at least `mintrim=INT` bp with less than `deptrim=INT`
    read depth. First, samtools `mpileup` or `depth` (`depmethod=X`) is used to find the first and last positions that
    exceed `deptrim=INT`. If no positions meet this criterio, the entire sequence will removed. Otherwise, either
    terminus that exceeds `mintrim=INT` base pairs of insufficent read depth are trimmed.

    ---
    ### ~ Region checking [runmode=regcheck] ~ ###

    Region checking, whether for read spanning analysis (`runmode=regcheck`) or copy number analysis
    (`runmode=regcnv` or `runmode=regcheck regcnv=T`), analyses regions extracted from a delimited file given by:
    `regcheck=FILE`. This can be a GFF file, in which case `gfftype=LIST` will restrict analysis to a specific feature
    type, or regions can be defined by `checkfields=LIST`, which defines the locus, start and end positions. These
    fields default to `SeqName`, `Start` and `End` fields. If these fields cannot be found, the first three fields
    of the `regcheck=FILE` file will be used.

    Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are then mapped onto the assembly
    (`seqin=FILE`) using minimap2 to generate a PAF file. This is then parsed and reads spanning each feature based
    on their positions and the target start and end positions in the PAF file. In addition to absolute spanning of
    regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
    sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
    using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
    due to sequence length constraints.

    If `regcnv=T` then the region copy number analysis (below) will also be performed.

    **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
    parallel calculations is set by `forks=INT`.

    ---
    ### ~ Region copy number analysis [runmode=regcnv] ~ ###

    Copy number analysis uses the same single copy depth profile analysis as the `runmode=gensize` genome size
    prediction. In short, the modal read depth of BUSCO single copy `Complete` genes is calculated using samtools
    `mpileup` (or samtools `depth` if `quickdepth=T`) and used to defined "single copy read depth". BUSCO single-copy
    genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
    `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
    ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with
    `Status` = `Complete` and the `Contig`, `Start` and `End` fields are used to define the regions that should be
    predominantly single copy.

    Single-copy read depth can also be set using `scdepth=NUM` to re-used or over-ride previously calculated values.

    Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are then mapped onto the assembly
    (`seqin=FILE`) using minimap2. This can be skipped by providing a BAM file with `bam=FILE`. For each regcheck
    feature, the same samtools depth calculation is performed as for the BUSCO data. The mean depth across the region
    is then divided by the single copy depth to estimate total copy number across the region. Note that unlike the
    single copy depth estimation itself, this might be biased by repeat sequences etc. that have a different depth
    profile to the main region. One way to control for this might be to restrict analysis to a subset of reads that
    meet a certain minimum length cutoff, e.g. 10kb.

    **Query-based CNV analysis.** If the `regcheck=FILE` file has additional `Qry`, `QryLen`, `QryStart` and `QryEnd`
    fields, the copy number analysi will have an additional query-focus. In this case, each region mapping onto a
    specific query is summed up, adjusted for the proportion of the query covered by that region. For example, 3.5X
    mean depth of a 100% length copy and 3.0X coverage of a 50% length copy would sum to (3.5x1.0 + 3x0.5 = 5 total
    copies). If these fields are not present, each region will be analysed independently.

    **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
    parallel calculations is set by `forks=INT`.

    ---

    ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

    This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
    existing) and then runs the read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
    with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
    This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
    In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
    sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
    using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
    due to sequence length constraints.

    Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
    `seqname.start-end`.

    This mode is now primarily documented and updated through GapSpanner.

    ---

    ### ~ Assembly gap re-assembly [runmode=gapass] ~ ###

    In addition to the `gapspan` analysis, reads identified as spanning each gap are extracted and assembled using `flye`
    in a `$BASEFILE__gapassemble/` output directory.

    This mode is now primarily documented and updated through GapSpanner.

    ---

    ### ~ Re-assembled gap-filling [runmode=gapfill] ~ ###

    In addition to the `gapspan` and `gapass` outputs, re-assembled gap regions are compiled into a single file and then
    mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
    are reduced to unique coverage of the assembly sequences. Gaps are filled if one of the two conditions are met:

    1. A single local alignment spans an entire gap.
    2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

    In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
    re-assembly contig region. For a pair of hits, the region between the two hits is replaced.

    This mode is now primarily documented and updated through GapSpanner.

    ---

    ### ~ Sorted non-redundant assembly cleanup [runmode=sortnr] ~ ###

    The sorted non-redundant assembly cleanup mode (`runmode=sortnr`) screens out any sequences that are 100% gap,
    then removes any sequences that are 100% redundant with other sequences in the input. This includes full and
    partial matches, i.e. if sequence X is wholly contained within sequence Y then X will be removed.

    First, sequences are loaded from the file given with `seqin=FILE` and any [rje_seqlist](http://rest.slimsuite.unsw.edu.au/seqlist)
    filters and sequence sorting are applied to the input. Sequences that are 100% Ns are removed and any gaps
    exceeding 10 nt are reduced to 10 `N`s (`NNNNNNNNNN`) to prevent minimap2 from splitting sequences on long gaps.
    These gap-reduced sequences are output to `$BASEFILE.tmp.fasta` and used for an all-by-all minimap2 search.

    By default, minimap2 is run with the options to generate a `$BASEFILE.tmp.paf` file:

        --cs -p 0.0001 -t 4 -x asm20 -N 250

    To modify minimap2 search settings, please see the [rje_paf](http://rest.slimsuite.unsw.edu.au/rje_paf)
    documentation.

    **NOTE:** These run options can probably be made more stringent to speed up minimap2 without loss of function.
    Future releases may alter defaults accordingly.

    Minimap2 output is parsed to identify scaffold-scaffold matches. Self-hits are ignored.
    The minimum (gap-reduced) sequence length is used as a rapid parsing filter: any minimap2 matches that are less
    than 95% of the query sequence (`Length`+`nn` fields) or less that 100% identity (`Identity`+`nn`)/(`Length`+`nn`)
    are filtered during parsing.

    **NOTE:** Future releases may feature an option to reduce the global percentage identity cut-off. Please contact
    the author if you wish to see this implemented.

    Minimap2 hits are then processed reverse-sorted by reference sequence size (e.g. scaffold length). Any hits
    where either sequence has already been filtered are skipped. Otherwise, if the match (as determined by the
    length of `:` regions in the CS string) matches the query length, the Query sequence will be flagged for
    remove as "identical to" or "contained within" the Hit. (Mutually partial overlapping exact matches are NOT
    filtered.) Filtered IDs and their matches are output to `$BASEFILE.redundant.txt`.

    Once all sequences have been filtered, the remaining sequences are output to: `$BASEFILE.nr.fasta`.

    **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
    add the `sortseq=invsize` command to the Diploidocus run command.

    Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to:
    `$BASEFILE.summarise.tdt`.

    Temporary gap-reduced and minimap2 PAF files are deleted unless running in `debug` or
    `dev` modes.

    ---
    ### ~ Pseudodiploid to primary and alternative haploptigs [runmode=diphap(nr)] ~ ###

    This protocol is based on 10x assemblies made for multiple organisms with supernova v2.0.0 and supernova v2.1.1.
    In each case, some redundancy in output was discovered (a) within pseudohap output, and (b) in terms of fully
    homozygous (identical) scaffolds shared by both haplotigs. It was also not entirely clear on what basis a
    particular haplotype was assigned to pseudohap1 or pseudohap2.

    The general workflow therefore sought to remove redundancy, generate a set of primary scaffolds based on scaffold
    length, and generate a non-redundant set of alternative scaffolds where heterozygosity exists. If `diphapnr` mode
    is used, the full workflow is implemented by first running the `sortnr` workflow described above. In the reduced
    `diphap` mode, redundancy is not removed first.

    Sequences are loaded and matching haplotigs identified based on their names. Sequence names MUST end `HAP(\d+)`,
    where `(\d+)` indicates an integer that matches up haplotigs (as produced by supernova pseudohap2 output, for
    example). This is **not** a pipeline to identify haplotig pairs, it is purely for splitting identified
    haplotigs into primary and alternative assemblies.

    Processing itself is quite simple. Haplotig pairs are identified based on matching `HAP(\d+)` numbers. Where a
    single haplotig is found, it is assigned as `diploid`, under the assumption that the two haplotigs were identical
    and one was removed. (It is possible that only one parent had this scaffold, e.g. sex chromosomes, so some post-
    processing of descriptions may be required.) If two haplotigs with the same number are identified, the longest
    is assigned to `haploidA` and the shorter `haploidB`.

    The **Primary Assemmbly** is then compiled from all `haploidA` and `diploid` sequences. These are given `pri`
    prefixes and output to `$BASEFILE.pri.fasta`. The **Alternative** comprised of all `haploidB` sequences is output
    to `$BASEFILE.alt.fasta`. If redundancy has been removed, this will likely be a subset of the full assembly. The
    combined set of all primary and alternative sequences is output to `$BASEFILE.dipnr.fasta`.

    **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
    add the `sortseq=invsize` command to the Diploidocus run command.

    Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to
    `$BASEFILE.summarise.tdt`:

    * `$BASEFILE.dipnr.fasta` = Combined pseudodiploid with `haploidA`, `haploidB` and `diploid` annotation.
    * `$BASEFILE.pri.fasta` = Primary assembly with `haploidA` and `diploid` sequences.
    * `$BASEFILE.alt.fasta` = Alternative assembly with `haploidB` sequences.


    ---
    ### ~ In silico diploid generator [runmode=insilico] ~ ###

    This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
    parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
    parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
    identifier table.)

    A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
    unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
    selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
    two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
    This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
    no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
    subreads will be added from the other parent if they reduce the difference in cumulative output for each parent, or
    until `lenfilter=X` is reached.

    Final output will be a `*.LXXXRQXX.fasta` file in which each parent has a similar total sequence content and for
    which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
    assemblies, where one parent has higher quality data than the other.

    NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
    higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
    minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
    relaxed. Similarly, only sequences above `lenfilter=X` in length will be output. These are the figures given in the
    `LXXXRQXX` part of the output file, e.g. defaults of RQ>=0.84 and Len>=500 generates `*.L500RQ84.fas`.


Commandline:
    ### ~ Main Diploidocus run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    runmode=X       : Diploidocus run mode (diploidocus/dipcycle/gensize/purgehap/telomeres/vecscreen/regcheck/regcnv/deptrim/sortnr/diphap/diphapnr/insilico) [diploidocus]
    basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    genomesize=INT  : Haploid genome size (bp) [0]
    scdepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
    ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth [True]
    mapadjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [False]
    ### ~ DiploidocusHocusPocus and Purge haplotigs options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    kmerreads=FILELIST : File of high quality reads for KAT kmer analysis []
    10xtrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
    minmedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
    minlen=INT      : Minimum scaffold length to avoid low quality filter [500]
    phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
    phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
    phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
    zeroadjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
    includegaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
    mingap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
    purgemode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
    diploidify=T/F  : Whether to generate alternative diploid output with duplicated diploid contigs and no hpurge [False]
    pretrim=T/F     : Run vectrim/vecmask and deptrim trimming prior to diploidocus run [False]
    maxcycle=INT    : Restrict run to maximum of INT cycles (0=No limit) [0]
    purgecyc=INT    : Minimum number of purged sequences to trigger next round of dipcycle [2]
    ### ~ Depth Trim options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    deptrim=INT     : Trim termini with <X depth [1]
    mintrim=INT     : Min length of terminal depth trimming [1000]
    ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
    telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
    telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
    teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
    ### ~ VecScreen options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    screendb=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
    screenmode=X    : Action to take following vecscreen searching (report/purge) [report]
    minvechit=INT   : Minimum length for a non-identical screendb match [50]
    minidhit=INT    : Minimum length for an identical screendb match (based on NCBI filtering) [27]
    efdr=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
    vecpurge=PERC   : Remove any scaffolds with >= PERC % vector coverage [50.0]
    vectrim=INT     : Trim any vector hits (after any vecpurge) within INT bp of the nearest end of a scaffold [1000]
    vecmask=INT     : Mask any vector hits of INT bp or greater (after vecpurge and vectrim) [900]
    maskmode=X      : Whether to perform full (all bases) or partial (every other base) masking of hits [partial]
    keepnames=T/F   : Whether to keep names unchanged for edited sequences or append 'X' [False]
    veccheck=T/F    : Check coverage of filtered contaminant hits using reads=FILELIST data [False]
    ### ~ Region checking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    regcheck=FILE   : File of SeqName, Start, End positions for read coverage checking [None]
    checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
    checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
    spanid=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
    regcnv=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
    gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    subforks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
    mingapspan=INT  : Minimum number of reads spanning a gap in order to re-assemble [2]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [500]
    ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
    ### ~ In silico diploid input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rqfilter=X      : Minimum RQ for output subreads [0]
    lenfilter=X     : Min read length for filtered subreads [500]
    parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
    ### ~ Advanced/Developmental options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    memperthread=INT: Number of Gb per thread to allocate to samtools etc. [6]
    useqsub=T/F     : Whether to use qsub to queue up system calls [False]
    qsubvmem=INT    : Memory setting (Gb) when queuing with qsub [126]
    qsubwall=INT    : Walltime setting (hours) when queuing with qsub [12]
    modules=LIST    : List of modules that needs to be loaded for running with qsub []
    legacy=T/F      : Run Legacy modes of updated methods that have been farmed out to other programs [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, math, os, re, string, subprocess, sys, time, shutil
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_forker, rje_obj, rje_rmd, rje_seqlist, rje_sequence, rje_paf #, rje_genomics
import rje_kat, rje_readcore, depthkopy, depthsizer
import rje_blast_V2 as rje_blast
import smrtscape
import slimfarmer
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Fixed bugs with parent basefile, genome size default and Sequel data parsing.
    # 0.2.0 - Added sortnr minimap2 run mode, and diphap primary/alternative assembly splitting.
    # 0.3.0 - Added summarising of sequences post-run. Improved documents. Moved to tools/.
    # 0.4.0 - Added VecScreen mode.
    # 0.5.0 - Added Telomere finding mode based on https://github.com/JanaSperschneider/FindTelomeres.
    # 0.6.0 - Added GenomeSize option.
    # 0.7.0 - Added PurgeHaplotigs, minimum vescreen hit length, eFDR and vecscreen coverage. [GitHub]
    # 0.7.1 - Code tidying and documentation improvements. Added dipcycle subdirectory and minlen=INT.
    #       - Fixed purge_haplotigs temp directory issues by running in subdirectory.
    #       - Added diploidify=T/F : Whether to generate alternative diploid output with duplicated diploid contigs and no hpurge [False]
    # 0.8.0 - Added deptrim=INT mode for trimming the ends of contigs based on low depth.
    # 0.8.1 - Fixed failure to screen contaminants in nala, simple and crude purgemode. Add dev minimap2 fix.
    # 0.8.2 - Shifted LOWQUAL to quarantine rather than junk. Fixed eFDR calculation error (inverse ranking).
    #       - Added veccheck=T/F : Check coverage of filtered contaminant hits using reads=FILELIST data [False]
    # 0.9.0 - Added pretrim, vecpurge, vecmask and vectrim options. Reduced screenmode to report/purge.
    # 0.9.1 - Fixed partial vecpurge when setting parameters to zero.
    # 0.9.2 - Added keepnames=T/F : Whether to keep names unchanged for edited sequences or append 'X' [False]
    # 0.9.3 - Fixed termination of program when BUSCO gene mpileup goes wrong (unless debug=T).
    # 0.9.4 - Added capacity to pick up an aborted dipcycle run. Added maxcycle=INT setting to limit run times.
    # 0.9.5 - Added use of qsub for system calls and memperthread=INT to control max memory.
    # 0.9.6 - Dropped default vecmask to 900bp as 1kb length matches picked up by NCBI are sometimes shorter on scaffold.
    # 0.9.7 - Check for | characters in the input sequences. (Will break the program.)
    # 0.9.8 - Added additional Mode of Modes genomes size prediction to log.
    # 0.10.0 - Added regcheck=TDTFILE and regcnv=T/F functions (runmode=regcheck)
    # 0.10.1 - Added gfftype=X optional feature type to use if performing regcheck on GFF file (e.g. gene) []
    # 0.10.2 - Fixed GFF regcheck bug. Add some documentation. Made gfftype=LIST not X.
    # 0.10.3 - Added BUSCO Complete Copy Number calculation and CI estimate to GenomeSize (and RegCNV). Made default gfftype=gene
    # 0.10.4 - Added improved SCdepth calculation option using BUSCO depth density profile.
    # 0.10.5 - Added 95% CI estimates to CNV calculations if BUSCO depth profiles are present.
    # 0.10.6 - Fixed GenomeSize bug with zero coverage BUSCO genes (e.g. using different long reads from assembly)
    # 0.10.7 - Updated the depmode.R script to limit analysis to one depmethod.
    # 0.11.0 - Updated defaults to depdensity=T for genome size prediction.
    # 0.12.0 - Added gapspan function (regcheck but first makes the gaps table, then loads this in, and outputs reads per gap.)
    # 0.13.0 - Added full gap-filling to gapfill function.
    # 0.14.0 - Added HiFi read type.
    # 0.14.1 - Added GapSpanner sequence summary to end.
    # 0.14.2 - Fixed scdepth re-use printing bug.
    # 0.15.0 - Added mapadjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [False]
    # 0.16.0 - Updated vecscreen mode to have two-tier length filter, and partial masking to avoid splitting contigs.
    # 0.16.1 - Added *.repeats.fasta output of the non-score diploidocus sequences.
    # 0.16.2 - Minor bug fix catching Class key error. Added logo to docs.
    # 0.16.3 - Updated some of the documentation.
    # 0.16.4 - Fixed a bug where pretrim of vecscreen results will cause BUSCO genes to be missed during classification.
    # 0.17.0 - Added purgecyc=INT : Minimum number of purged sequences to trigger next round of dipcycle [2]
    # 0.17.1 - Minor tweaks to log output.
    # 0.17.2 - Stopped CSI indexing from crashing Diploidocus but not compatible with PurgeHaplotigs.
    # 0.18.0 - Implementation of density-based CNV estimation (dev=T).
    # 1.0.0 - Updated to use rje_kat and rje_readcore and made v1.0.0 in line with Stuart et al. publication.
    # 1.0.1 - Bug fixes for GFF checkpos.
    # 1.1.0 - Fixed TeloRev sequence for finding 3' telomeres. Add reporting of telomere lengths.
    # 1.1.1 - Fixed DepthSizer object bug for DipCycle.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add download of NCBI Vector Database if vecdb=ncbi.
    # [Y] : Need to add a eukaryote mode and/or minlen for VecScreen - too many expected hits with NCBI rules.
    # [Y] : Add eFDR calculation and filtering to VecScreen.
    # [Y] : Implement long read BUSCO read depth and genome size prediction.
    # [Y] : Implement running of Purge haplotigs using BUSCO read depths to set cutoffs.
    # [Y] : Fix generation of warnings and error when Description is not available!
    # [Y] : Change purgehap to run purge_haplotigs only
    # [Y] : Add saving and reloading of total read count and depth?
    # [ ] : Add renaming of sequences to the dipnr analysis mode.
    # [Y] : Move dipcycle output into a subdirectory.
    # [Y] : Add minlength filter.
    # [Y] : Add diploidify=T/F : Duplicate "pure" diploid scaffolds (DIP, NON) and keep haplotigs [False]
    # [Y] : Change minlength=INT to minlen=INT.
    # [ ] : Add option for simplified sequence descriptions.
    # [ ] : Add auto-recognition of ont and pb subreads for BAM generation.
    # [ ] : Add katana-specific spawning of BUSCO analysis?
    # [Y] : Add vecpurge=PERC cutoff for ScreenPerc contamination filtering
    # [Y] : Add vecmask=INT cutoff to mask contamination above a certain length
    # [Y] : Add vectrim=INT cutoff to trim contamination within INT bp of sequence ends
    # [N] : Replace screenmode=X with vecpurge, vectrim, vecmask hierarchy.
    # [Y] : Add optional and pretrim=T/F for dipcycle to run vecscreen and/or depthtrim
    # [ ] : Add rDNA finding with barrnap. -> Feed into regcheck
    # [ ] : Add (optional?) re-use of the first vecscreen search if purgemode=dipcycle and no additional trimming.
    # [Y] : Add maxcycle=INT to terminate dipcycle after INT cycles.
    # [ ] : Document regcheck process.
    # [Y] : Test regcheck GFF mode.
    # [Y] : Add gapspan mode = same as regcheck but first makes the gaps table, then loads this in, and outputs reads per gap.
    # [Y] : Add gapass mode = gap reassembly, trying to assemble the reads spanning each gap
    # [Y] : Add gapfill mode = runs GABLAM of assembly versus original gap chunk and then tries to fill gaps?
    # [ ] : Add assembler=X = option to set for using different assemblers for gapass and gapfill modes.
    # [Y] : Add paf=FILE option to be used a bit like the bam=FILE input.
    # [Y] : Add minimum number of reads to assemble with gapass mode.
    # [Y] : Add option to switch off additional gapfill read span check.
    # [Y] : Add updated gapfill regcheck output to include gap edges. (For long replacements.)
    # [ ] : Add BUSCOMP generation of new BUSCO ratings using buscofas.
    # [ ] : Fix regcheck bug with provided PAF file.
    # [ ] : Separate out the spanning code from the read depth code and explicitly use PAF and BAM.
    # [ ] : Split out code to DepthSizer etc. so they don't all need all Diploidocus dependencies.
    # [ ] : Update docs to point to individual programs.
    # [Y] : Add purgecyc=INT : Minimum number of purged sequences to trigger next round of dipcycle [2]
    # [ ] : Add final output of input and output to *.tdt (rather than having to run summarise again).
    # [ ] : Add optional additional BAM files to collate depth stats for (e.g. short reads)
    # [ ] : Replace PurgeHaplotigs so that CSI indexing is OK.
    # [ ] : Add DensK and DensDep statistics for each sequence, using new DepthCopy code.
    # [Y] : Update to Version 1.0.
    # [ ] : Added auto-detection of regcheck=FILE and setting correct run mode if needed.
    # [Y] : Check 3' recognition of Telomeres and telomere output positions.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Diploidocus', '1.1.1', 'January 2022', '2017')
    description = 'Diploid genome assembly analysis toolkit.'
    author = 'Dr Richard J. Edwards.'
    comments = ['Citation: Chen SH et al. & Edwards RJ (preprint): bioRxiv 2021.06.02.444084 (doi: 10.1101/2021.06.02.444084)',
                'Please raise bugs or questions at https://github.com/slimsuite/diploidocus.',
                'NOTE: telomere finding rules are based on https://github.com/JanaSperschneider/FindTelomeres',
                rje_obj.zen()]
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
            if rje.yesNo('Show Minimap2 run (rje_paf) commandline options?',default='N'): out.verbose(-1,4,text=rje_paf.__doc__)
            if rje.yesNo('Show SeqList commandline options?',default='N'): out.verbose(-1,4,text=rje_seqlist.__doc__)
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
paf_defaults = {'N':'250','p':'0.0001','x':'asm20'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Diploidocus Class                                                                                       #
#########################################################################################################################
class Diploidocus(rje_readcore.ReadCore,rje_kat.KAT):
    '''
    Diploidocus Class. Author: Rich Edwards (2019).

    Str:str
    - BAM=FILE        : BAM file of reads mapped onto assembly [$BASEFILE.bam]
    - BUSCO=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    - MaskMode=X      : Whether to perform full (all bases) or partial (every other base) masking of hits [partial]
    - PAF=FILE        : PAF file of reads mapped onto assembly [$BASEFILE.paf]
    - Parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    - Parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    - PurgeMode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
    - RegCheck=TDTFILE: File of SeqName, Start, End positions for read coverage checking [None]
    - RunMode=X       : Diploidocus run mode [insilico/sortnr/diphap/vecscreen]
    - ScreenDB=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
    - ScreenMode=X    : Action to take following vecscreen searching (report/purge) [report]
    - SeqIn=FILE      : Input sequence assembly (sortnr/diphap modes) []
    - SeqOut=FILE     : Output sequence assembly [$BASEFILE.fasta]
    - SpanID=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
    - TeloFwd=X      : Basic telomere sequence for search [C{2,4}T{1,2}A{1,3}]
    - TeloRev=X      : Basic telomere sequence for search [T{1,3}A{1,2}G{2,4}]
    - TmpDir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]

    Bool:boolean
    - DepDensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth [True]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Diploidify=T/F  : Whether to generate alternative diploid output with duplicated diploid contigs and no hpurge [False]
    - Diploidocus=T/F : Whether to code is being run from a direct Diploidocus commandline call [False]
    - IncludeGaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
    - KeepNames=T/F   : Whether to keep names unchanged for edited sequences or append 'X' [False]
    - Legacy=T/F      : Run Legacy modes of updated methods that have been farmed out to other programs [False]
    - MapAdjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [False]
    - PreTrim=T/F     : Run vectrim/vecmask and deptrim trimming prior to diploidocus run [False]
    - PurgeCyc=INT    : Minimum number of purged sequences to trigger next round of dipcycle [2]
    - QuickDepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    - RegCNV=T/F      ; Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
    - Summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    - VecCheck=T/F    : Check coverage of filtered contaminant hits using reads=FILELIST data [False]
    - UseQSub=T/F     : Whether to use qsub to queue up system calls (dev only) [False]
    - ZeroAdjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
    - 10xTrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer reads data [False]

    Int:integer
    - DepTrim=INT     : Trim termini with <X depth [1]
    - GenomeSize=INT  : Haploid genome size (bp) [0]
    - LenFilter=X     : Min read length for filtered subreads [500]
    - MaxCycle=INT    : Restrict run to maximum of INT cycles (0=No limit) [0]
    - MemPerThread=INT: Number of Gb per thread to allocate to samtools etc. [6]
    - MinGap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
    - MinGapSpan=INT  : Minimum number of reads spanning a gap in order to re-assemble [2]
    - MinIDHit=INT    : Minimum length for an identical screendb match (based on NCBI filtering) [27]
    - MinLength=INT   : Minimum scaffold lenght to avoid low quality filter [500]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [0]
    - MinMedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
    - MinTrim=INT     : Min length of terminal depth trimming [1000]
    - MinVecHit=INT   : Minimum length for a screendb match [50]
    - PHLow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
    - PHMid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
    - PHHigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
    - QSubPPN=INT     : Number of CPU to use when queuing with qsub [16]
    - QSubVMem=INT    : Memory setting (Gb) when queuing with qsub [126]
    - QSubWall=INT    : Walltime setting (hours) when queuing with qsub [12]
    - ReadBP=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    - RegCNV=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
    - SubForks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
    - TeloSize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
    - VecMask=INT     : Mask any vectore hits of INT bp or greater (after vecpurge and vecmerge) [900]
    - VecTrim=INT     : Trim any vector hits (after any vecpurge) within INT bp of the nearest end of a scaffold [1000]

    Num:float
    - CheckCov=PERC   : Percentage coverage for double-checking partial exact matches [95]
    - eFDR=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    - RQFilter=X      : Minimum RQ for output subreads [0]
    - SCDepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    - TeloPerc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
    - VecPurge=PERC   : Remove any scaffolds with >= PERC % vector coverage [50.0]

    File:file handles with matching str filenames
    
    List:list
    - CheckFields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
    - CheckFlanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
    - GFFType=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    - KmerReads=FILELIST   : File of reads for KAT kmer analysis []
    - Reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database Object
    - Forker = Forker Object
    - SeqIn = SeqList Object for input sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM','BUSCO','GenomeSize','MaskMode','PAF','Parent1','Parent2','PurgeMode','RegCheck','RunMode','ScreenDB','ScreenMode','SeqIn','SeqOut','DebugStr','SpanID','TeloFwd','TeloRev','TmpDir']
        self.boollist = ['DepDensity','Diploidify','Diploidocus','DocHTML','IncludeGaps','KeepNames','Legacy','MapAdjust','PreTrim','QuickDepth','RegCNV','Summarise','UseQSub','VecCheck','ZeroAdjust','10xTrim']
        self.intlist = ['DepTrim','GenomeSize','LenFilter','MaxCycle','MinGap','MinGapSpan','MinIDHit','MinLen','MinMedian','MinTrim','MinVecHit',
                        'QSubPPN','QSubVMem','QSubWall','MemPerThread','MinLocLen','PurgeCyc',
                        'PHLow','PHMid','PHHigh','ReadBP','SubForks','TeloSize','VecMask','VecTrim']
        self.numlist = ['CheckCov','eFDR','MinLocID','RQFilter','SCDepth','TeloPerc','VecPurge']
        self.filelist = []
        self.listlist = ['CheckFields','CheckFlanks','GFFType','KmerReads','Reads','ReadType']
        self.dictlist = []
        self.objlist = ['Forker','SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self._setReadCoreAttributes()   # See rje_readcore
        self._setKatAttributes()        # See rje_kat
        self.setStr({'MaskMode':'partial','PurgeMode':'complex','RunMode':'diploidocus','ScreenMode':'report','TeloFwd':'C{2,4}T{1,2}A{1,3}','TeloRev':'T{1,3}A{1,2}G{2,4}','TmpDir':'./tmpdir/'})
        self.setBool({'DepDensity':True,'Diploidify':False,'DocHTML':False,'IncludeGaps':False,'KeepNames':False,
                      'Legacy':False,
                      'PreTrim':False,'QuickDepth':False,'RegCNV':True,'Summarise':True,'UseQSub':False,'ZeroAdjust':True,'10xTrim':False})
        self.setInt({'DepTrim':0,'LenFilter':500,'MaxCycle':0,'MinMedian':3,'MinIDHit':27,'MinLen':500,'MinTrim':1000,'MinVecHit':50,
                     'QSubPPN':16,'QSubVMem':126,'QSubWall':12,'MemPerThread':6,'MinGapSpan':2,'MinLocLen':500,'PurgeCyc':2,
                     'GenomeSize':0,'ReadBP':0,'SubForks':1,'TeloSize':50,'MinGap':10,'VecMask':900,'VecTrim':1000})
        self.setNum({'CheckCov':95.0,'eFDR':1.0,'RQFilter':0,'SCDepth':0,'TeloPerc':50.0,'VecPurge':50.0,'MinLocID':0.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['CheckFields'] = ['Hit','SbjStart','SbjEnd']
        self.list['CheckFlanks'] = [0,100,1000,5000]
        self.list['ReadType'] = ['ont']
        self.list['GFFType'] = ['gene']
        self.obj['SeqIn'] = None
        self.obj['Forker']  = rje_forker.Forker(self.log,['logfork=F','killmain=F']+self.cmd_list)
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
                self._readCoreCmd(cmd)  # Will set all the core commands recognised.
                self._katCmd(cmd)       # Set kat commands recognised.
                ### Class Options (No need for arg if arg = att.lower()) ###
                self._cmdRead(cmd,type='file',att='ScreenDB',arg='vecscreen')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['GenomeSize','DebugStr','MaskMode','PurgeMode','RunMode','ScreenMode','SpanID','TeloFwd','TeloRev'])   # Normal strings
                self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['BAM','PAF','Parent1','Parent2','RegCheck','ScreenDB','SeqIn','SeqOut','BUSCO'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DepDensity','Diploidify','Diploidocus','DocHTML','IncludeGaps','KeepNames','Legacy','MapAdjust','PreTrim','QuickDepth','RegCNV','Summarise','UseQSub','VecCheck','ZeroAdjust','10xTrim'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['DepTrim','LenFilter','MaxCycle','MemPerThread','MinGap','MinGapSpan','MinIDHit','MinLocLen','MinLen','MinMedian','MinTrim','MinVecHit','PurgeCyc','QSubPPN','QSubVMem','QSubWall','PHLow','PHMid','PHHigh','ReadBP','SubForks','TeloSize','VecMask','VecTrim'])   # Integers
                self._cmdReadList(cmd,'float',['eFDR','RQFilter','SCDepth']) # Floats
                self._cmdReadList(cmd,'perc',['CheckCov','MinLocID','TeloPerc','VecPurge']) # Percentage
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CheckFields','GFFType','ReadType'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['CheckFlanks'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['KmerReads','Reads']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStrLC('GenomeSize'):
            try: self.setInt({'GenomeSize':rje_seqlist.bpFromStr(self.getStrLC('GenomeSize'))})
            except:
                self.errorLog('Problem with GenomeSize. Setting genomesize=0')
                self.setInt({'GenomeSize':0})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main Diploidocus run method
        '''
        # Diploidocus: Diploid genome assembly analysis tools

        <img alt="Diploidocus logo" align="right" style="margin-left:10px;margin-bottom" width="150" src="https://raw.githubusercontent.com/slimsuite/diploidocus/master/docs/diploidocus_logo.png">
        Diploidocus is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
        The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
        BUSCO gene prediction and contaminant screening for a number of assembly tasks including contamination identification,
        haplotig identification/removal and low quality contig/scaffold trimming/filtering.

        In addition, Diploidocus will use mapped long reads and BUSCO single copy read depths for genome size prediction
        (`gensize`), and coverage (`regcheck`) or copy number estimation (`regcnv`) for user-defined regions. Diploidocus
        also has functions for removing redundancy (`sortnr`), generating a non-redundant pseudo-diploid assembly with primary
        and secondary scaffolds from 10x pseudohap output (`diphap`), and creating an in silico diploid set of long reads from two
        haploid parents (for testing phasing etc.) (`insilico`).

        Please note that Diploidocus is still in development and documentation is currently a bit sparse.

        The different run modes are set using `runmode=X`:

        * `diploidocus` default run mode will run `gensize`, `telomeres`, `vecscreen`, `deptrim` and `purgehap` analysis
        * `dipcycle` runs iterative cycles of `diploidocus` mode until convergence (no more filtering) is reached
        * `gensize` uses BUSCO results, a BAM file and read file(s) to predict the genome size of the organism
        * `purgehap` filters scaffolds based on post-processing of purge_haplotigs
        * `telomeres` performs a regex telomere search based on method of https://github.com/JanaSperschneider/FindTelomeres
        * `vecscreen` searches for contaminants and flags/masks/removes identified scaffolds
        * `deptrim` trims sequence termini of at least `mintrim=INT` bp with less than `deptrim=INT` read depth
        * `regcheck` checks reads spanning given regions and also calculates mean depth and estimated copy number (if regcnv=T)
        * `regcnv` calculates mean depth and estimated copy number for regcheck regions from BAM file (no read spanning analysis)
        * `gapspan` is a specialised `regcheck` analysis that checks for reads spanning genome assembly gaps
        * `gapass` extends the `gapspan` mode to extract the reads spanning a gap and re-assemble them using `flye`
        * `gapfill` extends the `gapass` mode to map re-assembled gaps back onto the assembly and replace gaps spanned by the new contigs        * `sortnr` performs an all-by-all mapping with minimap2 and then removes redundancy
        * `diphap` splits a pseudodiploid assembly into primary and alternative scaffolds
        * `diphapnr` runs `sortnr` followed by `diphap`
        * `insilico` generates balanced diploid combined reads from two sequenced haploid parents

        See <https://slimsuite.github.io/diploidocus/> for details of each mode. General SLiMSuite run documentation can be
        found at <https://github.com/slimsuite/SLiMSuite>.

        Diploidocus is available as part of SLiMSuite, or via a standalone GitHub repo at
        <https://github.com/slimsuite/diploidocus>.

        ## Citing Diploidocus

        If using Diploidocus in a publication, please cite: Edwards RJ et al. (2021), BMC Genomics [PMID: 33726677]. Not
        all of the Diploidocus functions were described in this paper. Future versions of the documentation will include
        a more detailed breakdown of appropriate citations. If in doubt, please contact the author.

        ---

        # Running Diploidocus

        Diploidocus is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/diploidocus.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [Diploidocus git repo](https://github.com/slimsuite/diploidocus), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [Diploidocus git repo](https://github.com/slimsuite/diploidocus)
        for running on example data.

        ## Dependencies

        For full functionality, Diploidocus needs a number of additional programs installed on the system:

        * Python 2.x and Python 3.x
        * KAT
        * BedTools
        * R
        * purge_haplotigs
        * bbmap
        * BLAST+
        * [Minimap2](https://github.com/lh3/minimap2) (added to environment `$PATH` or given with the `minimap2=PROG` setting)
        * [Samtools](http://www.htslib.org/)

        To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For Diploidocus documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main Diploidocus run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        runmode=X       : Diploidocus run mode (diploidocus/dipcycle/gensize/purgehap/telomeres/vecscreen/regcheck/regcnv/deptrim/sortnr/diphap/diphapnr/insilico) [diploidocus]
        basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        genomesize=INT  : Haploid genome size (bp) [0]
        scdepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
        bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
        readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
        quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
        depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth [True]
        mapadjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [False]
        ### ~ DiploidocusHocusPocus and Purge haplotigs options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        kmerreads=FILELIST : File of high quality reads for KAT kmer analysis []
        10xtrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
        minmedian=INT   : Minimum median depth coverage to avoid low coverage filter [3]
        minlen=INT      : Minimum scaffold length to avoid low quality filter [500]
        phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
        phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
        phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
        zeroadjust=T/F  : Add zero coverage bases to purge_haplotigs LowPerc and adjust total [True]
        includegaps=T/F : Whether to include gaps in the zero coverage bases for adjustment (see docs) [False]
        mingap=INT      : Minimum length of a stretch of N bases to count as a gap for exclusion [10]
        purgemode=X     : Rules used for purgehap analysis (simple/complex/nala) [complex]
        diploidify=T/F  : Whether to generate alternative diploid output with duplicated diploid contigs and no hpurge [False]
        pretrim=T/F     : Run vectrim/vecmask and deptrim trimming prior to diploidocus run [False]
        maxcycle=INT    : Restrict run to maximum of INT cycles (0=No limit) [0]
        ### ~ Depth Trim options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        deptrim=INT     : Trim termini with <X depth [1]
        mintrim=INT     : Min length of terminal depth trimming [1000]
        ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
        telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
        telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
        teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
        ### ~ VecScreen options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        screendb=FILE   : File of vectors/contaminants to screen out using blastn and VecScreen rules []
        screenmode=X    : Action to take following vecscreen searching (report/purge) [report]
        minvechit=INT   : Minimum length for a screendb match [50]
        minidhit=INT    : Minimum length for an identical screendb match (based on NCBI filtering) [27]
        efdr=NUM        : Expected FDR threshold for VecScreen queries (0 is no filter) [1.0]
        vecpurge=PERC   : Remove any scaffolds with >= PERC % vector coverage [50.0]
        vectrim=INT     : Trim any vector hits (after any vecpurge) within INT bp of the nearest end of a scaffold [1000]
        vecmask=INT     : Mask any vector hits of INT bp or greater (after vecpurge and vectrim) [900]
        maskmode=X      : Whether to perform full (all bases) or partial (every other base) masking of hits [partial]
        keepnames=T/F   : Whether to keep names unchanged for edited sequences or append 'X' [False]
        veccheck=T/F    : Check coverage of filtered contaminant hits using reads=FILELIST data [False]
        ### ~ Region checking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        regcheck=FILE   : File of SeqName, Start, End positions for read coverage checking [None]
        checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking (unless GFF) [Hit,SbjStart,SbjEnd]
        checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
        spanid=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
        regcnv=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
        gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
        mingapspan=INT  : Minimum number of reads spanning a gap in order to re-assemble [2]
        minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [500]
        ### ~ SortNR filtering/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkcov=PERC   : Percentage coverage for double-checking partial exact matches [95]
        seqout=FILE     : Output sequence assembly [$BASEFILE.nr.fasta]
        ### ~ In silico diploid input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        rqfilter=X      : Minimum RQ for output subreads [0]
        lenfilter=X     : Min read length for filtered subreads [500]
        parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
        parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
        See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
        ### ~ Advanced/Developmental options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        memperthread=INT: Number of Gb per thread to allocate to samtools etc. [6]
        useqsub=T/F     : Whether to use qsub to queue up system calls [False]
        qsubvmem=INT    : Memory setting (Gb) when queuing with qsub [126]
        qsubwall=INT    : Walltime setting (hours) when queuing with qsub [12]
        modules=LIST    : List of modules that needs to be loaded for running with qsub []
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        # Diploidocus run modes

        Details for the main Diploidocus run modes are given below.

        **NOTE:** Diploidocus is under development and documentation might be a bit sparse. Please contact the author or
        post an issue on GitHub if you have any questions.

        ---

        ## Cycled Diploidocus filtering [runmode=dipcycle]

        Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
        no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
        using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
        scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

        Output for each cycle will be initially generated in the run directory but then moved to a `dipcycle_$BASEFILE`
        directory upon completion.

        Final outputs from the final cycle will then be compiled under the original `$BASEFILE` prefix:

        * `$BASEFILE.diploidocus.tdt` = Final ratings for the input scaffolds. This is like the single cycle output with an additional `Cycle` field containing the last cycle this scaffold was processed in.
        * `$BASEFILE.ratings.tdt` = Reduced final ratings output for the input scaffolds (`SeqName`,`SeqLen`,`ScreenPerc`,`Class`,`Rating`,`Set`,`Cycle`).
        * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
        * `$BASEFILE.repeats.fasta` = the repeat scaffolds excluded from core
        * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
        * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
        * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

        **NOTE:** Contents for these four `*.fasta` files are summarised in the main log. Individual purge cycles have their own
        log files in the `dipcycle_$BASEFILE` directory.

        See `runmode=diploidocus` documentation for more details.

        ### Re-starting aborted runs

        If Diploidocus is killed before completion, files will not have been moved to the `dipcycle_$BASEFILE` directory.
        Any cycles that have generated a `*.diploidocus.fasta` file will be skipped until an incomplete cycle is found.
        In the rare case that Diploidocus was killed whilst generating fasta output, and `*.diploidocus.fasta` is thus
        incomplete, further cycles will be messed up. In this case, deleted the relevant `*.diploidocus.fasta` file
        before re-running.

        If Diploidocus is run with `runmode=dipcycle` and the `dipcycle_$BASEFILE` directory is found, it will abort.
        Rename or delete the directory to re-run Diploidocus with the same `$BASEFILE` setting, or change `basefile=X`.


        ---

        ## Main Diploidocus filtering [runmode=diploidocus]

        Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
        BUSCO results to reclassify scaffolds and partition them into:

        * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
        * `*.core.fasta` = the same set of scaffolds, minus repeats
        * `*.repeats.fasta` = the repeat scaffolds excluded from core
        * `*.junk.fasta` = low coverage scaffolds, removed as junk
        * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

        Unless `summarise=F`, summary statistics for sequence data before and after processing will be included in the run log.

        In addition, two key output tables are generated:

        * `*.diploidocus.tdt` = full compiled data and classification. (See field descriptions, below.)
        * `*.rating.tdt` = reduced data and classification for easier extraction of sequence names with `grep` etc.

        All output files are named with the `basefile=X` prefix. This is the prefix of the assembly given with `seqin=FILE` if no `basefile=X` setting is provided.

        ### Dependencies

        Diploidocus needs the following programs installed for full functionality:

        * `bbmap`
        * `blast+`
        * `kat`
        * `minimap2`
        * `purge_haplotigs`
        * `samtools`

        ### Input and assembly processing

        The main inputs for Diploidocus rating and filtering are:

        * `seqin=FILE` : Input sequence assembly to tidy [Required].
        * `screendb=FILE` : File of vectors/contaminants to screen out using blastn and VecScreen rules [Optional].
        * `reads=FILELIST`: List of fasta/fastq files containing long reads. Wildcard allowed. Can be gzipped. For a single run (not cycling), a BAM file can be supplied instead with `bam=FILE`. (This will be preferentially used if found, and defaults to `$BASEFILE.bam`.) Read types (pb/ont) for each file are set with `readtype=LIST`, which will be cycled if shorter (default=`ont`). Optionally, the pre-calculated total read length can be provided with `readbp=INT` and/or the pre-calculated (haploid) genome size can be provided with `genomesize=INT`.
        * `busco=TSVFILE` : BUSCO full table [`full_table_$BASEFILE.busco.tsv`] used for calculating single copy ("diploid") read depth. This can be over-ridden by setting `scdepth=NUM`.
        * `kmerreads=FILELIST` : File of high quality (_i.e._ short or error-corrected) reads for KAT kmer analysis [Optional]

        If a BAM file is not provided/found, Diploidocus will use minimap2 to generate a BAM file of `reads=FILELIST` data mapped onto the `seqin=FILE` assembly. Each read file is mapped separately (`--secondary=no -L -ax map-ont` or `--secondary=no -L -ax map-pb`) and converted into a sorted BAM files, before merging the BAM files with `samtools` and indexing the combined file.

        Diploidocus will re-use files where they already exist, providing the downstream files are newer than the upstream files. (If files have been copied and lost their datestamp information, switching `ignoredate=T` will re-use files regardless.) Setting `force=T` should force regeneration of files even if they exist.


        #### Purge Haplotigs

        Diploidocus uses output from purge_haplotigs as on of the core components for its classification system. Running purge_haplotigs is automated, using the single copy (diploid) read depth to set the purge_haplotigs thresholds. If `scdepth=NUM` is not set, this will be calculated using BUSCO single copy genes (see **Genome size prediction [runmode=gensize]** for details). The `purge_haplotigs` mid cutoff (`-m X`) will be set at the halfway point between the diploid depth (SCDepth) and the haploid depth (SCDepth/2). The low depth (`-l X`) and high depth (`-h X`) cutoffs will then be set to SCDepth/4 and SCDepth*2. This can be over-ridden with:

        * `phlow=INT` : Low depth cutoff for purge_haplotigs (`-l X`).
        * `phmid=INT` : Middle depth for purge_haplotigs (`-m X`).
        * `phhigh=INT` : High depth cutoff for purge_haplotigs (`-h X`).

        Diploidocus adjusts the `purge_haplotigs` binning proportions, which do not use the zero coverage bases, adding
        the zero coverage bases to the low coverage bin. This adjustment can be with `zeroadjust=F`. By default, gaps are
        excluded from this adjustment, defined as a run of 10+ `N` bases. Setting `includegaps=T` will include these
        regions in the adjustment, whilst changing the `mingap=INT` setting will redefine gaps.

        #### KAT analysis

        Next, the `kat sect` function of KAT is used to calculate kmer frequencies for (a) raw data using the `kmerreads=FILELIST` file of high quality (_i.e._ short or error-corrected) reads, and (b) the assembly itself. This generates four outputs:

        * `*.kat-stats.tsv` = summary `kmerreads` kmer coverage per sequence
        * `*.kat-counts.cvg` = `kmerreads` kmer coverage per base for each input sequence
        * `*.selfkat-stats.tsv` = summary assembly kmer coverage per sequence
        * `*.selfkat-counts.cvg` = assembly kmer coverage per base for each input sequence

        If the `kmerreads=FILELIST` files are 10x chromium reads, the first 16bp of read 1 (the 10x barcode) can be trimmed by setting `10xtrim=T`. (This sets `--5ptrim 16` in kat for the first read file.)

        #### BBMap depth statistics

        Read depth statistics for the BAM file are calculated per input sequence using the `pileup.sh` program of `bbmap`. This generates two files:

        * `*.depth.tdt` = read depth stats per sequence
        * `*.depstat` = overall read depth summary, generated by `pileup.sh` stdout and stderr

        #### Vector/contaminant coverage

        The `runmode=vecscreen` contaminant screening using vector/contaminants from `screendb=FILE` is run on the assembly - see VecScreen mode documentation for details.

        #### Telomeres

        The `runmode=telomeres` telomere identifcation screen is run on the assembly - see Telomeres mode documentation for details.

        #### BUSCO gene prediction

        The full BUSCO results table (`busco=FILE`) is used to generate counts of `Complete`, `Duplicate` and `Fragmented` BUSCO genes for each input sequence. These are converted into a single BUSCO rating of the dominant (i.e. most abundant) BUSCO class during classification.

        ### Data compilation

        Once all the different elements have been run, Diploidocus compiles results per input sequence, extracting a subset of fields of interest and renaming fields where appropriate. The following results files are compiled, joined by input sequence name (`SeqName`):

        * `*.depth.tdt`: bbmap pileup.sh read mapping statistics
            - `SeqName` = Scaffold sequence name
            - `SeqLen` = Scaffold sequence length (bp)
            - `Median_fold` = Median long read sequence depth for sequence
            - `Avg_fold` = Mean long read sequence depth for sequence
            - `Covered_percent` = Percentage of sequence covered by long reads
            - `Covered_bases` = Number of bp covered by long reads
            - `Plus_reads` = Number of reads mapped to positive strand
            - `Minus_reads` = Number of reads mapped to negative strand
            - `Read_GC` = GC content of reads mapping to sequence
        * `*.purge.coverage_stats.csv`: purge_haplotigs coverage stats
            - `LowPerc` = Proportion of bases with read depth below low coverage cutoff
            - `HapPerc` = Proportion of bases with read depth between low coverage and mid coverage cutoffs
            - `DipPerc` = Proportion of bases with read depth between mid coverage and high coverage cutoffs
            - `HighPerc` = Proportion of bases with read depth above high coverage
        * `*.purge.reassignments.tsv`: purge_haplotigs best sequence hits and reassignments
            - `TopHit` = scaffold with biggest sequence coverage match
            - `SecHit` = scaffold with second biggest sequence coverage match
            - `TopHitCov` = percentage sequence coverage by `TopHit`
            - `MaxHitCov` = combined percentage sequence coverage by all sequence hits
            - `PurgeHap` = purge_haplotigs reassignment rating
            - `TopNum` = number of sequences for which this sequence is the `TopHit`
            - `SecNum` = number of sequences for which this sequence is the `SecHit`
        * `*.selfkat-stats.tsv`: assembly kmer coverage stats
            - `SelfMedK` = Median assembly kmer frequency for sequence
            - `SelfAvgK` = Mean assembly kmer frequency for sequence
        * `*.kat-stats.tsv`: short read kmer coverage stats
            - `MedK` = Median short read kmer frequency for sequence
            - `AvgK` = Mean short read kmer frequency for sequence
            - `SeqGC` = GC content for sequence
            - `KPerc` = Percentage of non-gap kmers in sequence that are found in short read data
            - **NOTE:** In the absence of short read data, these values will be set to -1
        * `*.screencov.tdt`: Contaminant coverage from vecscreen mode
            - `ScreenPerc` = Total combined sequence coverage by contaminants
        * `*.telomeres.tdt`: Telomere prediction results
            - `Tel5` = Whether a 5' telomere is predicted
            - `Tel3` = Whether a 5' telomere is predicted
            - `Tel5Len` = Length in window chunks of 5' telomere
            - `Tel3Len` = Length in window chunks of 3' telomere
            - `TelPerc` = Percentage of sequence predicted to telomeres. (Crude calculation.)
        * `full_table_*.busco.tsv`
            - `Complete` = Number of BUSCO Complete genes in sequence
            - `Duplicated` = Number of BUSCO Duplicated genes in sequence
            - `Fragmented` = Number of BUSCO Fragmented genes in sequence
        * `$SEQIN`: input assembly
            - `N_bases` = total count of all unresolved (`N`) bases
            - `Gap_bases` = total count of all gap (`N`) bases as determined by `mingap=INT`
            - `SeqDesc` = sequence description


        ### Diploidocus classification

        Diploidocus classification adds the following fields to the compiled table:

        * `Class` = Six-part (`PURITY|DEPTH|HOM|TOP|MEDK|BUSCO`) Diploidocus classifcation (see below)
        * `TopClass` = Class of `TopHit` (if any)
        * `SecClass` = Class of `SecHit` (if any)

        Classification is based on a combination of commandline parameters and hard-coded cutoffs:

        * `minmedian=INT` : Minimum median depth coverage to avoid low coverage filter [3]
        * `pureperc=80` : Depth class percentage to classify as `PURE`
        * `rephitcov=250` : Min `MaxHitCov` value to classifty as `REPEAT`

        Each sequence is given classification in six different criteria (`PURITY|DEPTH|HOM|TOP|MEDK|BUSCO`) with two possible suffixes:

        1. `PURITY`: Purity of dominant read depth class
            - `PURE` = At least 80% of sequence in that depth bin
            - `GOOD` = At least 50% of sequence in that depth bin
            - `WEAK` = Under 50% of sequence in that depth bin
        2. `DEPTH`: Dominant read depth class
            - `LOWX` = Median read depth (`Median_fold`) fails to meet `minmedian=INT` criterion (default=3)
            - `LOW` = `LowPerc` read depth bin has highest percentage coverage (ties assigned to other class)
            - `HAP` = `HapPerc` read depth bin has highest percentage coverage (non-`DIP` ties addigned to `HAP`)
            - `DIP` = `DipPerc` read depth bin has highest percentage coverage (ties assigned to `DIP`)
            - `EXS` = `HighPerc` read depth bin has highest percentage coverage
        3. `HOM`: Homology/repeat status based on purge_haplotigs hits
            - `UNIQ` = No `TopHit`
            - `PART` = Partial (<50%) coverage of `TopHit`
            - `HAPL` = 50%+ `TopHit` coverage but no `SecHit`
            - `HOMO` = `TopHit` and `SecHit` but combined `MaxHitCov` < 250%
            - `REPT` = `TopHit` and `SecHit` and 250%+ combined `MaxHitCov`
        4. `TOP`: Top/Sec Hit status for sequence
            - `TOP` = `TopHit` for at least one other sequence
            - `SEC` = Not a `TopHit` but `SecHit` at least one other sequence
            - `NON` = Neither a `TopHit` nor `SecHit` for any other sequence
        5. `MEDK`: Assembly redundancy based on KAT assembly kmers
            - `PRI` = Over 50% unique kmers (`SelfMedK` = 1)
            - `ALT` = Median kmer frequency of two  (`SelfMedK` = 2)
            - `REP` = Median kmer frequency exceeds two  (`SelfMedK` > 2)
        6. `BUSCO`: Dominant BUSCO class (can help decide about risk)
            - `COMP` = 1+ Complete BUSCO genes and more Complete than Duplicated
            - `DUPL` = 1+ Duplicated BUSCO genes and more Duplicated than Complete
            - `FRAG` = 1+ Fragmented BUSCO genes and no Complete or Duplicated
            - `NONE` = No Complete, Duplicated or Fragmented BUSCO genes
        7. `+TEL`: If any telomeres are detected, `+TEL` is added
        8. `+VEC`: If any contamination is detected, `+VEC` is added



        ### Purge modes [purgemode=X]

        Diploidocus rating based on `purgemode=X` adds the following fields to the compiled table:

        * `Rating` = Main Diploidocus rating of sequence (see below)
        * `TopRating` = Rating of `TopHit` (if any)
        * `SecRating` = Rating of `SecHit` (if any)
        * `Set` = Diploidocus output set (`keep`/`repeat`/`quarantine`/`junk`)

        There are three basic `purgemode` implementations:

        * `purgemode=complex` [Default] uses all of the Diploidocus classification data to partition sequences into classes and then partition those classes into the main output sets.
        * `purgemode=simple` uses a simplified set of Diploidocus classification rules data to partition sequences into classes before partitioning those classes into the main output sets. This is conceptually easier than the `complex` mode but unlikely to produce as good results.
        * `purgemode=crude` uses a crude set of Diploidocus classification rules data to partition sequences into classes before partitioning those classes into the main output sets. This is more complicated than `simple` but still conceptually easier than the `complex` mode but unlikely to produce as good results.
        * `purgemode=nala` uses a simplified and conservative set of filtering rules used for the Nala German Shepherd Dog assembly.

        Details for these modes can be found below.

        Ratings are based on a combination of commandline parameters and hard-coded cutoffs:

        * `minlen=INT` : Minimum scaffold length to avoid low quality filter [500]
        * `mindipcov=20` : Minimim DipPerc coverage to assign a "KEEP" rating (`purgemode=crude`)
        * `artefactcov=80` : Min High/Low coverage to assign as artefacts (`purgemode=crude`, `purgemode=nala`)
        * `hpurgecov=80` : Min TopHitCov to assign `HPURGE` rating.
        * `covwarning=95` : Read coverage threshold for warnings

        #### Complex rating

        Diploidocus performs a hierarchical rating of scaffolds, based on their classifications and compiled data:

        1. Initial low quality filter:
            * `CONTAMINATION` = 50%+ identified contamination (`ScreenPerc`)
            * `LOWCOV` = Poor median read coverage (`Median_fold < minmedian=INT`)
            * `LOWQUAL` = Scaffolds below the sequence length set by `minlen=INT`
            * `LOWCOV` = Low read depth and 50%+ sequence without short read kmers (`PURE|LOW|...` or `GOOD|LOW|...`, and `MedK=0`)
            * `LOWCOV` = Low read depth and high frequency assembly kmers (`...|LOW|...|REP|...`)
        2. High quality scaffolds to keep:
            * `QUALITY` = Highest quality scaffolds: pure duploid, complete BUSCOs, no Duplicated BUSCOs (`PURE|DIP|UNIQ|...|PRI|COMP` & `Duplicated` = 0)
            * `FINAL` = As Quality but `PRI|FRAG` and `PRI|NONE` also allowed (`PURE|DIP|UNIQ|...|PRI|...` & `Duplicated` = 0)
            * `CORE` = Predominantly Diploid with <50% covered by `TopHit` and `SelfMedK`=1 (`[PURE/GOOD]|DIP|[UNIQ/PART]|...|PRI|...`)
            * `COREHAP` = Predominantly haploid but <50% covered by `TopHit` and 1+ Complete BUSCOs (`...|HAP|[UNIQ/PART]|...` & `Complete` > 0)
            * `PRIMARY` = Putative primary scaffold but with possible alternative scaffolds still in assembly and/or low quality regions (`...|DIP|[UNIQ/PART/HAPL/HOMO]|...|[PRI/ALT]|...`)
        3. Repetitive scaffolds to keep:
            * `PRIRPT` = Putative primary scaffold but >50% repeated (`...|DIP|[UNIQ/PART]|...|REP|...` or `[PURE/GOOD]|DIP|...|REP|...`)
            * `COLLAPSED` = High coverage scaffolds representing putative collapsed repeats. (Check for other contamination.) (`...|EXS|[UNIQ/PART]|...` or `[GOOD/PURE]|EXS|...|TOP|...`, or `WEAK|EXS|...` and less of the sequence in the Diploid coverage bin than below Diploid coverage.
            * `REPEAT` = Predominantly Diploid scaffolds that have major signs of redundancy, probably due to presence of alternative contigs (`...|DIP|REPT|...` or `WEAK|DIP|[HAPL/HOMO]|...|REP|...`)
        4. Haploid coverage scaffolds to keep (**NOTE**: any `HAPLOTIG` rating will be converted to `HAPRPT` if `...|REP|...`)
            * `HAPLOID` = Predominantly haploid coverage but enough unique sequence to keep for now? (Option to quarantine?) Might be very heterozygous Alternative haplotigs (`...|HAP|[UNIQ/PART]|PRI|...`)
            * `HAPLOTIG` = Predominantly haploid coverage but enough unique sequence to keep for now - possible Alternative haplotig (`...|HAP|[UNIQ/PART]|ALT|...` or `...|HAP|...|TOP|[PRI/ALT]|...`)
            * `HAPRPT` = Low quality scaffold that is probably most repeats, but not bad enough to dump outright (`...|HAP|[UNIQ/PART]|...|REP|...` or `...|HAP|...|TOP|...|REP|...`)
            * `HAPLOTIG` = Low/haploid coverage scaffold with insufficient coverage of a Diploid scaffold to purge (`...|[LOW/HAP]|[HAPL/HOMO/REPT]|...` and `|DIP|` in `TopHit` rating and `TopHitCov` < `hpurgecov` (80%))
            * `HAPLOTIG` = Low/haploid coverage scaffold with insufficient coverage of another scaffold to purge (`[PURE/GOOD]|HAP|HOMO|...' and `TopHitCov` < `hpurgecov` (80%))
            * `HAPRPT` = Other haploid coverage scaffold with insufficient coverage of another scaffold to purge (`...|HAP|...` and `TopHitCov` < `hpurgecov` (80%))
            * `HAPRPT` = Low coverage and collapsed repeats with insufficient coverage of another scaffold to purge (`...|EXS|...` and less of the sequence in the Diploid coverage bin than below Diploid coverage and `TopHitCov` < `hpurgecov` (80%))
            * `REPEAT` = High coverage scaffold with insufficient coverage of another scaffold to purge (`[GOOD/PURE]|EXS|[HAPL/HOMO/REPT]|[SEC/NON]|...` and `TopHitCov` < `hpurgecov` (80%))
        5. Scaffolds to quarantine for removal
            * `RPURGE` = Messy scaffolds that are largely repeats and are sufficiently redundant/low quality to purge (`WEAK|EXS` and less of the sequence in the Diploid coverage bin than below Diploid coverage, or `[GOOD/PURE]|EXS|[HAPL/HOMO/REPT]|[SEC/NON]`, and `TopHitCov` >= `hpurgecov` (80%))
        6. Additional low coverage scaffolds
            * `LOWCOV` = Low depth of coverage and strong matches to other scaffolds, and not TopHit (`...|LOW|[HAPL/HOMO/REPT]|[SEC/NON]`)
        7. Further scaffolds to quarantine for removal
            * `HPURGE` = Clear candidate haplotig to purge (`...|LOW|[HAPL/HOMO/REPT]|...` and `|DIP|` in `TopHit` rating, or `...|HAP|...`, and `TopHitCov` >= `hpurgecov` (80%))
        8. Final low quality filter
            * `LOWQUAL` = Unconvincing scaffolds that do not fall into a clear class but are not bad enough to dump outright as junk (`...|LOW|...`)

        Any scaffolds that escape the above rules will be `UNRATED`. (There should not be any!)


        #### Simple rating

        The simplified rating system is based much more on the dominant read depth class. A subset core scaffolds are identified and the rest are classified on the basis of dominant read depth and `REP` or `REPT` status.

        * `LOWCOV` = Poor median read coverage (`Median_fold < minmedian=INT`)
        * `CONTAMINATION` = 50%+ identified contamination (`ScreenPerc`)
        * `LOWQUAL` = Scaffolds below the sequence length set by `minlen=INT`
        * `LOWQUAL` = Low quality artefacts with less than 50% of the scaffold covered by short-read kmers (`MedK=0`)
        * Keep Scaffolds with Diploid+ depth (`[DIP/EXS]`), or limited homology to other scaffolds (`[UNIQ/PART]`) or Top hits for other scaffolds (`TOP`) or 1+ `Complete` BUSCO genes. These are broken down into:
            - `PRIMARY` = Diploid depth or no `TopHit` or dominant BUSCO rating is `Complete`, and not repetitive
            - `PRIRPT` = As `PRIMARY` but `REPT` or `REP` classifications
            - `COLLAPSED` = Dominant high depth (`EXS`)
            - `REPEAT` = Remaining Scaffolds to keep that have `REPT` or `REP` classifications
            - `HAPLOID` = Remaining Scaffolds to keep without `REPT` or `REP` classifications
        * `HPURGE` = Dominant haploid read depth without `REPT` or `REP` classifications
        * `RPURGE` = Dominant haploid read depth, plus `REPT` or `REP` classifications
        * `LOWCOV` = Any remaining scaffolds (dominant `LOW` read depth)

        #### Crude rating

        The `crude` purge mode is a less nuanced version of the main Diploidocus rating.

        * `LOWCOV` = Poor median read coverage (`Median_fold < minmedian=INT`)
        * `CONTAMINATION` = 50%+ identified contamination (`ScreenPerc`)
        * `LOWQUAL` = Scaffolds below the sequence length set by `minlen=INT`
        * `HPURGE` = Low+Haploid coverage meet `artefactcov` cutoff (>=80%) and `TopHitCov` meets `hpurgecov` cutoff (>=80%)
        * `HAPLOTIG` = Haploid coverage exceeds 50%, `TopHitCov` meets `hpurgecov` cutoff (>=80%), and `SelfMedK=2`
        * Scaffolds with a Diploid depth coverage otherwise meeting the `mindipcov` cutoff (>=20%) are kept:
            - `PRIMARY` = Median assembly kmer frequency of 1
            - `REPEAT` = `MaxHitCov` meets the `rephitcov` cutoff (>=250%)
            - `KEEP` = Remaining diploid scaffolds.
        * `LOWCOV` = Any other scaffolds with 80%+ low coverage bases are filtered as Low Coverage.
        * `COLLAPSED` = Any other scaffolds with 80%+ high coverage bases are flagged as COLLAPSED.
        * `HAPLOTIG` = Scaffolds that are 50%+ Haploid Depth if mostly present 2+ times (median assembly kmer frequency = 2)
        * `HAPLOID` =  Scaffolds that are 50%+ Haploid Depth if mostly present once (median assembly kmer frequency = 1)
        * `ARTEFACT` = Scaffolds <50% Haploid Depth but 80%+ Low or Haploif depth
        * `HAPLOTIG` = Scaffolds with 80%+ hitting other scaffolds but does not meet the `rephitcov` cutoff (<250%)
        * `HAPRPT` = Scaffolds with 80%+ hitting other scaffolds and meets the `rephitcov` cutoff (>=250%)
        * `ARTEFACT` = Remaining scaffolds with 80%+ bases in the low/haploid coverage bins


        #### Nala rating

        The `purgemode=nala` rating scheme was used for the Nala German Shepherd Dog genome assembly, and features a simplified set of ratings:

        * `CONTAMINATION` = 50%+ identified contamination (`ScreenPerc`)
        * `LOWCOV` = Poor median read coverage (`Median_fold < minmedian=INT`)
        * `LOWQUAL` = Scaffolds below the sequence length set by `minlen=INT`
        * `HPURGE` = Any scaffold with 80%+ bases in the low/haploid coverage bins (haplotigs or assembly artefacts).
        * `PRIMARY` = Scaffolds with 20%+ diploid coverage are marked as retention as probable diploids.
        * `COLLAPSED` = Scaffolds with <20% diploid coverage and 50%+ high coverage are marked as probable collapsed repeats.
        * `JUNK`/`HAPLOTIG`/`KEEP` = Remaining Scaffolds are given the PurgeHaplotigs rating (over 80% low/high coverage will be filtered as a probable artefact)


        #### Final sequence sets for output

        The `purgemode` ratings are then converted into subsets of scaffolds for output:

        * `keep` = Scaffolds that form part of the `*.core.fasta` and `*.diploidocus.fasta` outputs:
            - ['QUALITY','FINAL','CORE','COREHAP','PRIMARY','PRIRPT','HAPLOID','HAPLOTIG','KEEP']
        * `repeat` = Scaffolds excluded from the `*.core.fasta` but part of the `*.diploidocus.fasta` output:
            - ['COLLAPSED', 'REPEAT', 'HAPRPT']
        * `quarantine` = Scaffolds removed as probably Haplotigs and/or excess copies of repeats but saved to `*.quarantine.fasta` in case additional checks are required:
            - ['HPURGE', 'RPURGE']
        * `junk` = Scaffolds removed due to low coverage, poor quality, or contamination. These are partitioned into `*.junk.fasta`:
            - ['LOWCOV', 'LOWQUAL', 'CONTAMINATION','JUNK']

        The final ratings for each `TopHit` and `SecHit` are also added to the main Diploidocus table as `TopHitRating` and `SecHitRating` fields. This might flag up additional scaffolds (e.g. with `HAPLOTIG` ratings) that the user decides to remove despite failing to meet the automated criteria.

        **NOTE:** A cutdown output of Scaffold classifications and ratings is saved as `*.ratings.tdt`. This can be used to easily extract a subset of sequence names that can be filtered from the final output using `SeqSuite`, _e.g._

        ```
        python $SLIMSUITE/tools/seqsuite.py seqin=$PREFIX badseq=$(grep HAPLOTIG | awk '{print $1;}' ORS=",") seqout=$PREFIX.filtered.fasta
        ```

        #### Rating warnings

        Once rating is complete, the following warnings are generated in the log file:

        * All `CONTAMINATION` scaffolds
        * Scaffolds not classed as `junk` but:
            - <50% covered by short read kmers
            - any vector contamination flagged
            - long read coverage < 95%
        * Any `COREHAP/HAPLOID/HAPLOTIG` scaffolds with `Duplicated` BUSCO genes (more evidence for redundancy)
        * Any purged (`junk` or `quarantine`) scaffolds that are the PurgeHaplotigs `TopHit` for 1+ other scaffolds


        ### Full Diploidocus outputs

        Temporary files generated during forking will be generated in a subdirectory set by `tmpdir=PATH` (Default: `./tmpdir/`).

        _Details of full outputs will be added. Please see other documentation sections for specific outputs, or contact the author if anything is unclear._

        ### Special output mode: diploidify

        Diploidocus features a special `diploidify=T` mode, which aims to generate fully diploid output. This is achieved by adding `HPURGE` scaffolds to the `keep` set rather than the `quarantine` set. In addition, diploid scaffolds are duplicated. For this purpose, diploid scaffolds are defined as those in the `keep` or `repeat` sets that meet are:

        * Classified as `...|DUP|...`
        * or Classified as `WEAK|...|NON|...`, _i.e._ are neither a `TopHit` nor `SecHit` for any other sequence

        Diploid scaffolds are saved with an `X2` suffix on their sequence name.

        ---

        ## DepthSizer Genome size prediction [runmode=gensize]

        **NOTE:** This mode is now primarily documented and updated through [DepthSizer](https://github.com/slimsuite/depthsizer).

        The main inputs for Diploidocus genome size prediction are:

        * `seqin=FILE` : Input sequence assembly to tidy [Required].
        * `reads=FILELIST`: List of fasta/fastq files containing long reads. Wildcard allowed. Can be gzipped. For a single run (not cycling), a BAM file can be supplied instead with `bam=FILE`. (This will be preferentially used if found, and defaults to `$BASEFILE.bam`.) Read types (pb/ont/hifi) for each file are set with `readtype=LIST`, which will be cycled if shorter (default=`ont`). Optionally, the pre-calculated total read length can be provided with `readbp=INT` and/or the pre-calculated (haploid) genome size can be provided with `genomesize=INT`.
        * `busco=TSVFILE` : BUSCO full table [`full_table_$BASEFILE.busco.tsv`] used for calculating single copy ("diploid") read depth. This can be over-ridden by setting `scdepth=NUM`.
        * `quickdepth=T/F`  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [Default: `False`]

        This works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid read
        depth) regions along with some ooor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
        are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
        region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
        read depth.

        Diploidocus uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
        and extracts the modal read depth for each BUSCO gene along with the overall modal read depth for all gene
        regions. Genome size is then estimated based on a crude calculation using the total combined sequencing length.
        This will be caculated from `reads=FILELIST` unless provided with `readbp=INT`.

        BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
        and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.

        **NOTE:** There is currently no adjustment for contamination, read mapping or imbalanced insertion:deletion
        ratios etc. As a consequence, the current genome size prediction appears to be an over-estimate.

        ---

        ## Running Purge_haplotigs using BUSCO-guided cutoffs [runmode=purgehap]

        Diploidocus can automate the running of purge_haplotigs is automated, using the single copy (diploid) read depth
        to set the purge_haplotigs thresholds. If `scdepth=NUM` is not set, this will be calculated using BUSCO single
        copy genes (see **Genome size prediction [runmode=gensize]** for details).

        The `purge_haplotigs` mid cutoff (`-m X`) will be set at the halfway point between the diploid depth (SCDepth)
        and the haploid depth (SCDepth/2). The low depth (`-l X`) and high depth (`-h X`) cutoffs will then be set to
        SCDepth/4 and SCDepth*2. This can be over-ridden with:

        * `phlow=INT` : Low depth cutoff for purge_haplotigs (`-l X`).
        * `phmid=INT` : Middle depth for purge_haplotigs (`-m X`).
        * `phhigh=INT` : High depth cutoff for purge_haplotigs (`-h X`).

        Output from this run is:

        * `*.purge.coverage_stats.csv`: purge_haplotigs coverage stats
        * `*.purge.reassignments.tsv`: purge_haplotigs best sequence hits and reassignments

        Additional purge_haplotigs files and output will be saved in a subdirectory:

        * `purge_*/`

        This is to stop multiple purge_haplotigs runs from interfering with each other.

        ---

        ## Telomere finding [runmode=telomere]

        Diploidocus has implemented a regex-based search for Telomeres, based on the code at
        https://github.com/JanaSperschneider/FindTelomeres. This looks for a canonical telomere motif of TTAGGG/CCCTAA,
        allowing for some variation. For each sequence, Diploidocus trims off any trailing Ns and then searches for
        telomere-like sequences at sequence ends. For each sequence, the presence/absence and length of trimming are
        reported for the 5' end (`tel5` and `trim5`) and 3' end (`tel3` and `trim3`), along with the total percentage
        telomeric sequence (`TelPerc`).

        By default, Diploidocus searches for a forward telomere regex sequence of `C{2,4}T{1,2}A{1,3}` at the 5' end, and
        a reverse sequence at the 3' end of `T{1,3}A{1,2}G{2,4}`. These can be set with `telofwd=X` and `telorev=X`.
        Telomeres are marked if at least 50% (`teloperc=PERC`) of the terminal 50 bp (`telosize=INT`) matches the
        appropriate regex. If either end contains a telomere, the total percentage of the sequence matching either
        regex is calculated as `TelPerc`. Note that this number neither restricts matches to the termini, not includes
        sequences within predicted telomeres that do not match the regex.

        ---

        ## Vector/contamination screening [runmode=vecscreen]

        This mode screens scaffolds for possible contaminants, given by `screendb=FILE`.

        ### VecScreen overview

        First, a `blastn` search of the `screendb=FILE` sequences is performed, using the NCBI VecScreen search and
        match strategy (see below). These parameters can be modified using `blastopt=X` and/or `optionfile=FILE`,
        which will be appended to the run command. Alternatively, the `$BASEFILE.vecscreen.blast` file produced can be
        pre-made with different parameters (if `force=F`).

        Results are then parsed into a local hits table and rated according to the strength of a match. This is performed
        iteratively, re-assigning internal matches as "proximal" if they are withing 25 bases of another match (of any
        strength) and re-rating the strength of the match, until no ratings changes occur. Two additional fields are
        added to the local hits table during this process: `MatchPos` and `MatchStr`. Once all assignments have been
        made, segments of the assembly between two matches and/or sequence ends are added to the table as `Suspect`.

        `MatchPos` will have a value of:

        * `Terminal` = within 25 bp of either end of the sequence.
        * `Proximal` = within 25 bp of a vecsreen match (`Weak`, `Moderate` or `Strong`).
        * `Internal` = over 25 bp from a sequence end or vecsreen match.
        * `Suspect` = Segments added as `Suspect`.

        `MatchStr` will have a value of:

        * `Strong` = `Terminal`/`Proximal` match with Score >= 24, or `Internal` match with Score >= 30.
        * `Moderate` = `Terminal`/`Proximal` match with Score 19 to 23, or `Internal` match with Score 25 to 29.
        * `Weak` = `Terminal`/`Proximal` match with Score 16 to 18, or `Internal` match with Score 23 to 24.
        * `Suspect` = Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        These parameters seem to be designed for small (Illumina?) assemblies, and give excessive false positives for
        larger genomes. Diploidocus therefore features two additional contaminant identification filters:

        * `eFDR` = The "Expected False Discovery Rate` is calculated for each contaminant sequence. This is simply the
        `Expect` value for that hit, divided by the total number of hits with that `Expect` value (or lower). If `efdr=X`
        is set, any hits with an `eFDR` value exceeding the threshold will be removed. By default, this is set to `1.0`,
        i.e. any contaminants with fewer assembly hits than their `Expect` value will be dropped.
        * `Minimum hit length` = Short hits are inevitable in large assemblies and unlikely to be real for long read
        assemblies (e.g. without cloning etc.). An additional `minvechit=X` setting will remove any really short hits.
        This is set by default to `minvechit=50`, meaning that a hit of at least 50 bp is required.

        Finally, the percentage coverage per scaffold is calculated from the filtered hits. This is performed first for
        each contaminant individually, before being collapsed into total contamination coverage per query.

        Results are output into two main delimited results files:

        * `*.vecscreen.tdt` = Contaminant local hit table with VecScreen ratings and eFDR calculation. (Unfiltered)
        * `*.screencov.tdt` = Query coverage per contaminant
            - `Query` = Contaminant name. For total coverage, this will be `TOTAL`.
            - `Hit` = Scaffold name.
            - `HitLen` = Length of scaffold (bp).
            - `Coverage` = Covered length of scaffold (bp).
            - `CovPC` = Percentage coverage of scaffold (0-100).

        ### VecScreen BLAST+ parameters (from NCBI website)

        The VecScreen BLAST+ parameters are pre-set using blastn options: -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000

        VecScreen Match Categories
        Vector contamination usually occurs at the beginning or end of a sequence; therefore, different criteria are
        applied for terminal and internal matches. VecScreen considers a match to be terminal if it starts within 25
        bases of the beginning of the query sequence or stops within 25 bases of the end of the sequence. Matches that
        start or stop within 25 bases of another match are also treated like terminal matches. Matches are categorized
        according to the expected frequency of an alignment with the same score occurring between random sequences.

        Strong Match to Vector
        (Expect 1 random match in 1,000,000 queries of length 350 kb.)
        Terminal match with Score >= 24.
        Internal match with Score >= 30.

        Moderate Match to Vector
        (Expect 1 random match in 1,000 queries of length 350 kb.)
        Terminal match with Score 19 to 23.
        Internal match with Score 25 to 29.

        Weak Match to Vector
        (Expect 1 random match in 40 queries of length 350 kb.)
        Terminal match with Score 16 to 18.
        Internal match with Score 23 to 24.

        Segment of Suspect Origin
        Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        ### Vector coverage checking

        If `veccheck=T` then long reads given by `reads=FILELIST readtype=LIST` will be mapped onto the assembly using
        minimap2 and complete coverage of each hit reported. This will report the number of reads that span the entrie
        contaminant hit, plus flanking distances each side set by `checkflanks=LIST` (default: 0,100,1000,5000). If the
        flanking distance (bp) or the end of the sequence is reached before the end of the read in both directions, that
        read will be added to the coverage count. Because this can inflate apparent flanking coverage for hits near the
        end of a sequence, the distance from each end of the sequence is also returned as `MaxFlank5` and `MaxFlank3`.

        ---

        ## Depth trimming [runmode=deptrim]

        Depth trimming (`deptrim`) mode trims sequence termini of at least `mintrim=INT` bp with less than `deptrim=INT`
        read depth. First, samtools `mpileup` or `depth` (`depmethod=X`) is used to find the first and last positions that
        exceed `deptrim=INT`. If no positions meet this criterio, the entire sequence will removed. Otherwise, either
        terminus that exceeds `mintrim=INT` base pairs of insufficent read depth are trimmed.

        ---

        ## Assembly region read-spanning and copy number analysis [runmode=regcheck/regcnv]

        Region checking, whether for read spanning analysis (`runmode=regcheck`) or copy number analysis
        (`runmode=regcnv` or `runmode=regcheck regcnv=T`), analyses regions extracted from a delimited file given by:
        `regcheck=FILE`. This can be a GFF file, in which case `gfftype=LIST` will restrict analysis to specific feature
        types, or regions can be defined by `checkfields=LIST`, which defines the locus, start and end positions. These
        fields default to `SeqName`, `Start` and `End` fields. If these fields cannot be found, the first three fields
        of the `regcheck=FILE` file will be used.

        ### Region read-spanning analysis [runmode=regcheck]

        Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly
        (`seqin=FILE`) using minimap2 to generate a PAF file. This is then parsed and reads spanning each feature based
        on their positions and the target start and end positions in the PAF file. In addition to absolute spanning of
        regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
        sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
        using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
        due to sequence length constraints.

        If `regcnv=T` then the region copy number analysis (below) will also be performed.

        **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
        parallel calculations is set by `forks=INT`.

        ---

        ### Region copy number analysis [runmode=regcnv]

        Copy number analysis uses the same single copy depth profile analysis as the [DepthSizer](https://github.com/slimsuite/depthsizer)
        (Diploidocus `runmode=gensize`) genome size prediction. In short, the modal read depth of BUSCO single copy
        `Complete` genes is calculated using samtools
        `mpileup` (or samtools `depth` if `quickdepth=T`) and used to defined "single copy read depth". BUSCO single-copy
        genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with
        `Status` = `Complete` and the `Contig`, `Start` and `End` fields are used to define the regions that should be
        predominantly single copy.

        Single-copy read depth can also be set using `scdepth=NUM` to re-used or over-ride previously calculated values.

        Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are then mapped onto the assembly
        (`seqin=FILE`) using minimap2. This can be skipped by providing a BAM file with `bam=FILE`. For each regcheck
        feature, the same samtools depth calculation is performed as for the BUSCO data. The mean depth across the region
        is then divided by the single copy depth to estimate total copy number across the region. Note that unlike the
        single copy depth estimation itself, this might be biased by repeat sequences etc. that have a different depth
        profile to the main region. One way to control for this might be to restrict analysis to a subset of reads that
        meet a certain minimum length cutoff, e.g. 10kb.

        **Query-based CNV analysis.** If the `regcheck=FILE` file has additional `Qry`, `QryLen`, `QryStart` and `QryEnd`
        fields, the copy number analysis will have an additional query-focus. In this case, each region mapping onto a
        specific query is summed up, adjusted for the proportion of the query covered by that region. For example, 3.5X
        mean depth of a 100% length copy and 3.0X coverage of a 50% length copy would sum to (3.5x1.0 + 3x0.5 = 5 total
        copies). If these fields are not present, each region will be analysed independently.

        **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
        parallel calculations is set by `forks=INT`.

        ---

        ## GapSpanner functions [runmode=gapspan/gapass/gapfill]

        **NOTE:** These modes are now primarily documented and updated through [GapSpanner](https://github.com/slimsuite/gapspanner).


        ### Assembly gap read-spanning analysis [runmode=gapspan]

        This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
        existing) and then runs the read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
        with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
        This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
        In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
        sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
        using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
        due to sequence length constraints.

        Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
        `seqname.start-end`.

        ---

        ### Assembly gap re-assembly [runmode=gapass]

        In addition to the `gapspan` analysis, reads identified as spanning each gap are extracted and assembled using `flye`
        in a `$BASEFILE__gapassemble/` output directory. Only gaps with at least `mingapspan=INT` (default 2) reads are
        re-assembled.

        ---

        ### Re-assembled gap-filling [runmode=gapfill]

        In addition to the `gapspan` and `gapass` outputs, re-assembled gap regions are compiled into a single file and then
        mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
        are reduced to unique coverage of the assembly sequences. Gaps are filled if one of the two conditions are met:

        1. A single local alignment spans an entire gap.
        2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

        In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
        re-assembly contig region. For a pair of hits, the region between the two hits is replaced.

        ---

        ## Sorted non-redundant assembly cleanup [runmode=sortnr]

        The sorted non-redundant assembly cleanup mode (`runmode=sortnr`) screens out any sequences that are 100% gap,
        then removes any sequences that are 100% redundant with other sequences in the input. This includes full and
        partial matches, i.e. if sequence X is wholly contained within sequence Y then X will be removed.

        First, sequences are loaded from the file given with `seqin=FILE` and any [rje_seqlist](http://rest.slimsuite.unsw.edu.au/seqlist)
        filters and sequence sorting are applied to the input. Sequences that are 100% Ns are removed and any gaps
        exceeding 10 nt are reduced to 10 `N`s (`NNNNNNNNNN`) to prevent minimap2 from splitting sequences on long gaps.
        These gap-reduced sequences are output to `$BASEFILE.tmp.fasta` and used for an all-by-all minimap2 search.

        By default, minimap2 is run with the options to generate a `$BASEFILE.tmp.paf` file:

            --cs -p 0.0001 -t 4 -x asm20 -N 250

        To modify minimap2 search settings, please see the [rje_paf](http://rest.slimsuite.unsw.edu.au/rje_paf)
        documentation.

        **NOTE:** These run options can probably be made more stringent to speed up minimap2 without loss of function.
        Future releases may alter defaults accordingly.

        Minimap2 output is parsed to identify scaffold-scaffold matches. Self-hits are ignored.
        The minimum (gap-reduced) sequence length is used as a rapid parsing filter: any minimap2 matches that are less
        than 95% of the query sequence (`Length`+`nn` fields) or less that 100% identity (`Identity`+`nn`)/(`Length`+`nn`)
        are filtered during parsing.

        **NOTE:** Future releases may feature an option to reduce the global percentage identity cut-off. Please contact
        the author if you wish to see this implemented.

        Minimap2 hits are then processed reverse-sorted by reference sequence size (e.g. scaffold length). Any hits
        where either sequence has already been filtered are skipped. Otherwise, if the match (as determined by the
        length of `:` regions in the CS string) matches the query length, the Query sequence will be flagged for
        remove as "identical to" or "contained within" the Hit. (Mutually partial overlapping exact matches are NOT
        filtered.) Filtered IDs and their matches are output to `$BASEFILE.redundant.txt`.

        Once all sequences have been filtered, the remaining sequences are output to: `$BASEFILE.nr.fasta`.

        **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
        add the `sortseq=invsize` command to the Diploidocus run command.

        Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to:
        `$BASEFILE.summarise.tdt`.

        Temporary gap-reduced and minimap2 PAF files are deleted unless running in `debug` or
        `dev` modes.

        ---

        ## Pseudodiploid to primary and alternative haploptigs [runmode=diphap(nr)]

        This protocol is based on 10x assemblies made for multiple organisms with supernova v2.0.0 and supernova v2.1.1.
        In each case, some redundancy in output was discovered (a) within pseudohap output, and (b) in terms of fully
        homozygous (identical) scaffolds shared by both haplotigs. It was also not entirely clear on what basis a
        particular haplotype was assigned to pseudohap1 or pseudohap2.

        The general workflow therefore sought to remove redundancy, generate a set of primary scaffolds based on scaffold
        length, and generate a non-redundant set of alternative scaffolds where heterozygosity exists. If `diphapnr` mode
        is used, the full workflow is implemented by first running the `sortnr` workflow described above. In the reduced
        `diphap` mode, redundancy is not removed first.

        Sequences are loaded and matching haplotigs identified based on their names. Sequence names MUST end `HAP(\d+)`,
        where `(\d+)` indicates an integer that matches up haplotigs (as produced by supernova pseudohap2 output, for
        example). This is **not** a pipeline to identify haplotig pairs, it is purely for splitting identified
        haplotigs into primary and alternative assemblies.

        Processing itself is quite simple. Haplotig pairs are identified based on matching `HAP(\d+)` numbers. Where a
        single haplotig is found, it is assigned as `diploid`, under the assumption that the two haplotigs were identical
        and one was removed. (It is possible that only one parent had this scaffold, e.g. sex chromosomes, so some post-
        processing of descriptions may be required.) If two haplotigs with the same number are identified, the longest
        is assigned to `haploidA` and the shorter `haploidB`.

        The **Primary Assemmbly** is then compiled from all `haploidA` and `diploid` sequences. These are given `pri`
        prefixes and output to `$BASEFILE.pri.fasta`. The **Alternative** comprised of all `haploidB` sequences is output
        to `$BASEFILE.alt.fasta`. If redundancy has been removed, this will likely be a subset of the full assembly. The
        combined set of all primary and alternative sequences is output to `$BASEFILE.dipnr.fasta`.

        **NOTE:** By default, sequences are output in the same order as in the input. To output in reverse size order,
        add the `sortseq=invsize` command to the Diploidocus run command.

        Finally, the input and output files are summarised (unless `summarise=F`) and statistics output to
        `$BASEFILE.summarise.tdt`:

        * `$BASEFILE.dipnr.fasta` = Combined pseudodiploid with `haploidA`, `haploidB` and `diploid` annotation.
        * `$BASEFILE.pri.fasta` = Primary assembly with `haploidA` and `diploid` sequences.
        * `$BASEFILE.alt.fasta` = Alternative assembly with `haploidB` sequences.


        ---

        ## In silico diploid generator [runmode=insilico]

        This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
        parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
        parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
        identifier table.)

        A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
        unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
        selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
        two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
        This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
        no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
        subreads will be added from the other parent if they reduce the difference in cumulative output for each parent, or
        until `lenfilter=X` is reached.

        Final output will be a `*.LXXXRQXX.fasta` file in which each parent has a similar total sequence content and for
        which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
        assemblies, where one parent has higher quality data than the other.

        NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
        higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
        minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
        relaxed. Similarly, only sequences above `lenfilter=X` in length will be output. These are the figures given in the
        `LXXXRQXX` part of the output file, e.g. defaults of RQ>=0.84 and Len>=500 generates `*.L500RQ84.fas`.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.getBool('DocHTML'): return self.docHTML()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#MODE',self.getStrLC('RunMode'))
            if self.getStrLC('RunMode') in ['sortnr','seqnr','nrseq','nr']: return self.sortNR()
            elif self.getStrLC('RunMode') in ['diphap','diphapnr']: return self.dipHap()
            elif self.getStrLC('RunMode') == 'insilico':
                return self.inSilicoHybrid()
            elif self.getStrLC('RunMode') == 'vecscreen':
                if self.getStr('ScreenMode') == 'purge': return self.vecPurge()
                else: return self.vecScreen()
            elif self.getStrLC('RunMode').startswith('telomere'): return self.findTelomeres()
            elif self.getStrLC('RunMode') == 'diploidocus': return self.diploidocusHocusPocus()
            elif self.getStrLC('RunMode').startswith('purgehap'): return self.diploidocusHocusPocus()
            elif self.getStrLC('RunMode') in ['gensize','genomesize']:
                if not self.getBool('Legacy'):
                    self.infoLog('Running DepthSizer for gensize mode (legacy=F)')
                    return self.genomeSize()
                self.infoLog('Running legacy gensize mode (legacy=T)')
                return self.legacyGenomeSize(makebam=True)
            elif self.getStrLC('RunMode') in ['dipcycle','purgecycle']: return self.purgeCycle()
            elif self.getStrLC('RunMode') in ['deptrim']: return self.depthTrim()
            elif self.getStrLC('RunMode') in ['gapspan','gapass','gapfill']: return self.gapSpan()
            elif self.getStrLC('RunMode') in ['regcheck']: return self.regCheck()
            elif self.getStrLC('RunMode') in ['regcnv']:
                if self.getBool('Legacy'):
                    self.infoLog('Running legacy regcnv mode (legacy=T)')
                    return self.regCheck()
                self.infoLog('Running DepthKopy for regcnv mode (legacy=F)')
                depcmd = ['winsize=0'] + self.cmd_list + ['regfile={0}'.format(self.getStr('RegCheck'))]
                depcmd += ['checkfields={0}'.format(','.join(self.list['CheckFields']))]
                if not depthkopy.DepthKopy(self.log,depcmd).run(): raise ValueError('DepthKopy failed')
                return True
            else: raise ValueError('RunMode="%s" not recognised!' % self.getStrLC('RunMode'))
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def purgeCycle(self):  ### Repeat Diploidocus purge cycles to convergence (none removed).                    # v0.7.0
        '''
        ### ~ Cycled Diploidocus filtering [runmode=dipcycle] ~ ###
        Diploidocus can be automated to cycle through repeated rounds of the main purge_haplotigs/filtering method until
        no further scaffolds are removed. Each cycle will run Diploidocus with a numerical suffix, e.g. `$BASEFILE.1.*`,
        using the `*.diploidocus.fasta` output from the previous cycle as input. Cycling will continue until no further
        scaffolds are filtered into either `*.quarantine.fasta` or `*.junk.fasta`.

        Output for each cycle will be initially generated in the run directory but then moved to a `dipcycle_$BASEFILE`
        directory upon completion.

        Final outputs from the final cycle will then be compiled under the original `$BASEFILE` prefix:

        * `$BASEFILE.diploidocus.tdt` = Final ratings for the input scaffolds. This is like the single cycle output with an additional `Cycle` field containing the last cycle this scaffold was processed in.
        * `$BASEFILE.ratings.tdt` = Reduced final ratings output for the input scaffolds (`SeqName`,`SeqLen`,`ScreenPerc`,`Class`,`Rating`,`Cycle`).
        * `$BASEFILE.diploidocus.fasta` = the scaffolds kept from the final Diploidocus cycle
        * `$BASEFILE.core.fasta` = the same set of scaffolds, minus repeats
        * `$BASEFILE.quarantine.fasta` = concatenated purged scaffolds from all Diploidocus cycles.
        * `$BASEFILE.junk.fasta` = concatenated low coverage and low quality scaffolds, removed as junk, from all cycles.

        **NOTE:** Contents for these four `*.fasta` files are summarised in the main log. Individual purge cycles have their own
        log files in the `dipcycle_$BASEFILE` directory.

        See `runmode=diploidocus` documentation for more details.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = self.baseFile(strip_path=True)
            newbase = cycbase = basefile #'dipcycle_{}/{}'.format(basefile,basefile)
            cycdir = 'dipcycle_{}/'.format(basefile)
            if rje.exists(cycdir):
                self.printLog('#CYCLE','Completed dipcycle directory found: {}'.format(cycdir))
                self.printLog('#ABORT','Aborting run - delete or rename existing dipcycle directory to re-run.')
                return False
            cycle = 0
            prevseqx = 0
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            db = self.db()
            dipdb = None
            ## ~ [1a] ~ Special pretim trimming ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# When cycling, this is done once before the cycling
            if self.getBool('PreTrim'):
                seqin = self.preTrim()  ### Performs vecscreen and deptrim trimming, updates self.seqinObj() and returns trimmed fasta file
            else:
                self.warnLog('Diploidocus run with pretrim=F (default): check results for signs of vector contamination.')

            ### ~ [2] ~ Cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maxstop = False; purgestop = False
            while seqlist.seqNum() != prevseqx:
                # Check maxcyle
                if self.getInt('MaxCycle') > 0 and cycle >= self.getInt('MaxCycle'):
                    self.printLog('#CYCLE','Max cycle reached (maxcycle={}). Finishing run.'.format(self.getInt('MaxCycle')))
                    maxstop = True
                    break
                # Check purgecycle
                if self.getInt('PurgeCyc') > 0 and prevseqx:
                    purgex = (seqlist.seqNum() - prevseqx)
                    if purgex < self.getInt('PurgeCyc'):
                        self.printLog('#CYCLE','Min sequence purging not exceeded (purgecycle={}). Finishing run.'.format(self.getInt('PurgeCyc')))
                        purgestop = True
                        break
                # Increment cycle and make new basefile
                prevseqx = seqlist.seqNum()
                cycle += 1
                oldbase = newbase
                newbase = '{}.{}'.format(cycbase,cycle)
                seqout = '{}.diploidocus.fasta'.format(newbase)
                if rje.exists(seqout) and self.force():
                    self.printLog('#ABORT','{} found and force=T: delete/move aborted previous run files and try again.'.format(seqout))
                    return False
                if rje.exists(seqout):
                    self.printLog('#CYCLE','Cycle {} skipped: {} found (force=F)'.format(cycle,seqout))
                    if dipdb:
                        cycdb = db.addTable('{}.diploidocus.tdt'.format(newbase),mainkeys=['SeqName'],expect=True,name='diploidocus.{}'.format(cycle))
                    else:
                        cycdb = db.addTable('{}.diploidocus.tdt'.format(newbase),mainkeys=['SeqName'],expect=True,name='diploidocus')
                else:
                    self.printLog('#SEQIN','Using {} as input for Diploidocus cycle {}.'.format(seqin,cycle))
                    self.printLog('#CYCLE','Cycle {}: running as {}.*'.format(cycle,newbase))
                    # Run Diploidocus
                    info = makeInfo()
                    cyccmd = ['i=-1']+self.cmd_list+['basefile={}'.format(newbase),'runmode=diploidocus','seqin=%s' % seqin]
                    if self.debugging(): cyccmd.append('i=1')
                    #i# Do not perform veccheck if already been performed! (Potentially, some
                    if self.getBool('PreTrim') or cycle > 1:
                        if self.getBool('VecCheck'):
                            cyccmd.append('veccheck=F')
                            self.printLog('#CHECK','NOTE: setting veccheck=F for cycle {} (previously performed)'.format(cycle))
                    cmd_list = rje.getCmdList(cyccmd,info=info)   # Reads arguments and load defaults from program.ini
                    out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
                    out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
                    out.printIntro(info)                                # Prints intro text using details from Info object
                    log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
                    if self.debugging():
                        cycdb = Diploidocus(log,['dna=T']+cmd_list+['i=1']).run()
                    else:
                        cycdb = Diploidocus(log,['dna=T','i=-1']+cmd_list).run()
                if not cycdb:
                    raise IOError('Cycle %s failed! Check %s.log. Aborting run.' % (cycle,newbase))
                cycdb.addField('Cycle',evalue=cycle)
                # Update main table
                if not dipdb: dipdb = cycdb
                else:
                    for cycseq in cycdb.dataKeys(): dipdb.data()[cycseq] = cycdb.data()[cycseq]
                # Check and process output
                if not rje.exists(seqout):
                    raise IOError('Expected %s output for Cycle %s not found! Check %s.log. Aborting run.' % (seqout,cycle,newbase))
                self.printLog('#CYCLE','Cycle {} complete. See {}.log for details.'.format(cycle,newbase))
                seqin = seqout
                seqlist = self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file','seqin=%s' % seqin,'autofilter=F'])

            ### ~ [3] Tidy up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if maxstop:
                if self.i() >= 0 and not rje.yesNo('Tidy up run as if convergence reached?'):
                    self.warnLog('Cycle data not tidied: re-run with higher maxcycle=INT to resume')
                    return True
                self.warnLog('#Cycling terminated: re-run on {}.diploidocus.fasta output to resume tidying'.format(basefile))
            elif purgestop:
                if self.i() >= 0 and not rje.yesNo('Tidy up run as if convergence reached?'):
                    self.warnLog('Cycle data not tidied: re-run with lower purgecyc=INT to resume')
                    return True
                self.warnLog('#Cycling terminated: re-run on {}.diploidocus.fasta output to resume tidying'.format(basefile))
            else:
                self.printLog('#CYCLE','Convergence achieved after cycle {}'.format(cycle))
            rje.mkDir(self,cycdir,log=True)
            ## ~ [3a] Transfer main data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # The main diploidocus ratings table, and reduced *.ratings.tdt table, are compiled during cycling. These
            # are saved as primary $BASFILE.* output.
            dipdb.baseFile(basefile)
            dipdb.saveToFile(sfdict={'LowPerc':4, 'HapPerc':4, 'DipPerc':4, 'HighPerc':4})
            dipdb.saveToFile(filename='%s.ratings.tdt' % basefile, savefields=['SeqName','SeqLen','ScreenPerc','Class','Rating','Cycle'])
            for ext in ['diploidocus.fasta','core.fasta','repeats.fasta']:
                if rje.exists('{}.{}'.format(newbase,ext)):
                    rje.backup(self,'{}.{}'.format(basefile,ext))
                    shutil.copy('{}.{}'.format(newbase,ext),'{}.{}'.format(basefile,ext))
                    self.printLog('#COPY','Output copied: {}.{} -> {}.{}'.format(newbase,ext,basefile,ext))
            for ext in ['quarantine.fasta','junk.fasta']:
                cycfiles = []
                for i in range(1,cycle+1):
                    ibase = '{}.{}'.format(cycbase,i)
                    if rje.exists('{}.{}'.format(ibase,ext)): cycfiles.append('{}.{}'.format(ibase,ext))
                if cycfiles:
                    rje.backup(self,'{}.{}'.format(basefile,ext))
                    os.system('cat {} > {}.{}'.format(' '.join(cycfiles),basefile,ext))
                    self.printLog('#COPY','Output copied: {} *.{} files -> {}.{}'.format(len(cycfiles),ext,basefile,ext))
            ## ~ [3b] Cleanup cycle data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in range(1,cycle+1):
                ibase = '{}.{}'.format(cycbase,i)
                ifiles = glob.glob('{}.*'.format(ibase))
                if rje.exists('purge_{}'.format(ibase)): ifiles.append('purge_{}'.format(ibase))
                for ifile in ifiles: os.rename(ifile,'{}{}'.format(cycdir,ifile))
                self.printLog('#MOVE','{} {}.* files moved to {}.'.format(len(ifiles),ibase,cycdir))

            ### ~ [4] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist.getBool('Summarise'):
                for seqset in ['diploidocus','core','repeats','quarantine','junk']:
                    setfas = '%s.%s.fasta' % (basefile,seqset)
                    if not rje.baseFile(setfas): continue
                    seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % setfas,'autofilter=F']
                    rje_seqlist.SeqList(self.log,seqcmd)


        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def preTrim(self):  ### Performs vecscreen and deptrim trimming, updates self.seqinObj() and returns trimmed fasta file
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS PRE-TIDY TRIMMING',line='=')
            seqin = self.getStr('SeqIn')
            basefile = self.baseFile(strip_path=True)

            #i# Run VecScreen if appropriate
            vecout = '{}.vecscreen.fasta'.format(basefile)
            trimdb = None
            if not self.force() and rje.exists(vecout):
                seqin = vecout
                self.printLog('#PURGE','{} found (force=F): skipping VecPurge'.format(vecout))
            else:
                trimdb = self.vecPurge()
            if trimdb:
                seqin = vecout
                self.setStr({'SeqIn':seqin})
                purged = trimdb.indexDataList('Edit','mask','SeqName')
                if purged:
                    self.warnLog('%s purged sequences will not be in main Diploidocus output: %s' % (rje.iLen(purged),', '.join(purged)))
            #i# Run DepTrim if appropriate
            depdb = None
            depout = '{}.trim.fasta'.format(basefile)
            if not self.force() and rje.exists(depout):
                seqin = depout
                self.printLog('#TRIM','{} found (force=F): skipping DepTrim'.format(depout))
            elif self.getInt('DepTrim'): depdb = self.depthTrim()
            if depdb:
                seqin = depout
                self.setStr({'SeqIn':seqin})
                purged = depdb.indexDataList('DepTrim','dump','SeqName')
                if purged:
                    self.warnLog('%s dumped sequences will not be in main Diploidocus output: %s' % (rje.iLen(purged),', '.join(purged)))
            return seqin
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
                else: self.baseFile('diploidocus')
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
            if self.getNum('SCDepth'):
                self.printLog('#SCDEP','Single copy read depth (scdepth=NUM) = {0:.2f}X'.format(self.getNum('SCDepth')))
            if self.getInt('GenomeSize'):
                self.printLog('#GSIZE','Genome size (genomesize=INT) = {0}'.format(rje_seqlist.dnaLen(self.getInt('GenomeSize'))))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    #i# This method is being added to ReadCore.
    def LEGACYseqinObj(self,summarise=True): ### Returns the a SeqList object for the SeqIn file
        '''
        Returns the a SeqList object for the SeqIn file.
        :return: self.obj['SeqIn']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['SeqIn']:
                seqcmd = self.cmd_list
                if summarise: seqcmd = ['summarise=T']+self.cmd_list
                self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
                sx = 0.0; stot = self.obj['SeqIn'].seqNum()
                for seq in self.obj['SeqIn'].seqs():
                    self.progLog('\r#CHECK','Checking sequences names: %.1f%%' % (sx/stot)); sx += 100.0
                    if '|' in self.obj['SeqIn'].shortName(seq):
                        raise ValueError('Pipe "|" characters found in seqin=FILE names: will break program. Please rename and try again.')
                self.printLog('\r#CHECK','Checking sequences names complete.')
        except ValueError:
            self.printLog('\r#CHECK','Checking sequences names aborted.')
            self.errorLog('Diploidocus input sequence error'); raise
        except:
            self.errorLog('Diploidocus.seqinObj() error')
        return self.obj['SeqIn']
#########################################################################################################################
    #!# Replace with rje_rmd.docHTML(self)
    def docHTML(self):  ### Generate the Diploidocus Rmd and HTML documents.                                        # v0.1.0
        '''Generate the Diploidocus Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2019 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown Diploidocus documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    #!# Try replacing with rje_obj version
    def LEGACYloggedSysCall(self,cmd,syslog=None,stderr=True,append=True,verbosity=1,nologline='WARNING: No run log output!',threaded=True):    ### Makes a system call, catching output in log file
        '''
        Makes a system call, catching output in log file.
        :param cmd:str = System call command to catch
        :param syslog:str [None] = Filename for system log in which to capture output. Will use $BASEFILE.sys.log if None.
        :param stderr:bool [True] = Whether to also capture the STDERR
        :param append:bool [False] = Whether to append the log if it exists
        :param verbosity:int [1] = Verbosity level at which output also goes to screen (tee, not redirect)
        :param nologline:str [None] = Default logline returned if nothing is in the log.
        :param threaded:bool [True] = Whether to use multiple PPN when using dev qsub mode.
        :return: last line of syslog output
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mydir = os.path.abspath(os.curdir)
            #i# System log filename
            if not syslog: syslog = '{}.sys.log'.format(self.baseFile())
            #i# Setup log with command and identify position to capture extra content
            headline = '[{}] {}\n'.format(rje.dateTime(),cmd)
            if append: open(syslog,'a').write(headline)
            else:
                rje.backup(self,syslog,appendable=False)
                open(syslog,'w').write(headline)
            #i# Identify position at end of file
            fend = rje.endPos(filename=syslog)
            #i# Generate system command
            if stderr:
                if ' > ' in cmd: self.warnLog('Cannot capture stderr or stdout for command: {}'.format(cmd))
                else: cmd = '{} 2>&1'.format(cmd)
            if append:
                if self.v() >= verbosity: cmd = '{} | tee -a {}'.format(cmd,syslog)
                else: cmd = '{} >> {}'.format(cmd,syslog)
            else:
                if self.v() >= verbosity: cmd = '{} | tee {}'.format(cmd,syslog)
                else: cmd = '{} > {}'.format(cmd,syslog)
            ### ~ [2] ~ Process System Call ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('UseQSub'):
                if not rje.exists('tmp_qsub'): rje.mkDir(self,'tmp_qsub/',log=True)
                qbase = rje.baseFile(syslog)
                ppn = self.threads()
                vmem = self.getInt('QSubVMem')
                if not threaded: ppn = 1; vmem = 16
                farmcmd = ['farm={}'.format(cmd),'jobwait=T','qsub=T','basefile={}'.format(qbase),'slimsuite=F',
                           'qpath={}'.format(mydir),'monitor=F',
                           'ppn={}'.format(ppn),'vmem={}'.format(vmem),'walltime={}'.format(self.getInt('QSubWall'))]
                self.printLog('#DEV','Using SLiMFarmer to farm system call in tmp_qsub/')
                farmer = slimfarmer.SLiMFarmer(self.log,self.cmd_list+farmcmd)
                self.printLog('#DEV','SLiMFarmer commands: {}'.format(' '.join(farmer.cmd_list)))
                if not farmer.list['Modules']:
                    for mod in os.popen('module list 2>&1').read().split():
                        if '/' in mod: farmer.list['Modules'].append(mod)
                    if farmer.list['Modules']:
                        modlist = ','.join(farmer.list['Modules'])
                        self.printLog('#MOD','Read modules for qsub from environment: {}'.format(modlist))
                        self.cmd_list.append('modules={}'.format(modlist))
                        farmer = slimfarmer.SLiMFarmer(self.log,self.cmd_list+farmcmd)
                        self.printLog('#DEV','SLiMFarmer commands: {}'.format(' '.join(farmer.cmd_list)))
                    else:
                        raise ValueError('Trying to run with useqsub=T but modules=LIST not set!')
                os.chdir('tmp_qsub/')
                excode = farmer.run()
                os.chdir(mydir)
            elif self.dev():
                self.printLog('#CALL',cmd)
                excode = subprocess.call(cmd)
            else:
                self.printLog('#SYS',cmd)
                excode = os.system(cmd)
            if excode > 0: raise ValueError('Non-zero exit status for: {}'.format(cmd))
            logline = nologline
            SYSLOG = open(syslog,'r')
            SYSLOG.seek(fend)
            for readline in SYSLOG.readlines():
                readline = rje.chomp(readline)
                if readline: logline = readline
            self.printLog('#SYSEND',logline)
            return logline
        except:
            self.errorLog('Diploidocus.loggedSysCall() error')
            os.chdir(mydir)
            return None
#########################################################################################################################
    ### <3> ### In silico hybrid run method                                                                             #
#########################################################################################################################
    def inSilicoHybrid(self):  ### Filter and combine subreads from parent and output to fasta file.
        '''
        Filter and combine subreads from parent and output to fasta file.

        This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
        parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
        parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
        identifier table.)

        A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
        unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
        selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
        two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
        This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
        no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
        subreads will be added from the other parent if they reduce the difference in cumulative output for each parent.

        Final output will be a `*.subreads.fasta` file in which each parent has a similar total sequence content and for
        which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
        assemblies, where one parent has higher quality data than the other.

        NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
        higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
        minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
        relaxed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Parent 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 1 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent1: %s' % self.getStr('Parent1'))
            base1 = rje.baseFile(self.getStr('Parent1'))
            parent1 = smrtscape.SMRTSCAPE(self.log,self.cmd_list+['batch=%s' % self.getStr('Parent1'),'basefile=%s' % base1])
            parent1.setup()
            udb1 = parent1.udb()
            cdb = parent1.db('smrt',add=True,mainkeys=['Name'])
            cdb.dataFormat({'SMRT':'int'})
            cx = cdb.entryNum()
            ## ~ [0a] Parent 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 2 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent2: %s' % self.getStr('Parent2'))
            base2 = rje.baseFile(self.getStr('Parent2'))
            parent2 = smrtscape.SMRTSCAPE(self.log,self.cmd_list+['batch=%s' % self.getStr('Parent2'),'basefile=%s' % base2])
            parent2.setup()
            udb2 = parent2.udb()
            cdb2 = parent2.db('smrt',add=True,mainkeys=['Name'])
            cdb2.dataFormat({'SMRT':'int'})
            # Shift all of the Parent2 SMRT IDs to avoid conflict with Parent1
            for entry in cdb2.entries() + udb2.entries(): entry['SMRT'] = entry['SMRT'] + cx
            cdb = parent1.db().mergeTables(cdb,cdb2)
            ## ~ [0c] Output Sequence File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ DIPLOIDOCUS SUBREADS ~~~~~~~~~~~~~~~~~~~~ #')
            minlen = self.getInt('LenFilter')
            minrq = self.getNum('RQFilter')
            rqstr = '%s' % minrq
            filtfile = '%s.L%sRQ%s.fasta' % (self.baseFile(),minlen,rqstr[2:])
            ## ~ [0d] Input Sequence Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqbatch = []   # List of SeqList objects
            self.printLog('#BATCH','%s sequence files to process.' % rje.iLen(parent1.list['Batch']+parent2.list['Batch']))
            for seqfile in parent1.list['Batch']+parent2.list['Batch']:
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=F','seqin=%s' % seqfile,'autofilter=F']
                seqbatch.append(rje_seqlist.SeqList(self.log,seqcmd))
            self.printLog('#BATCH','%s sequence files to summarise.' % rje.iLen(seqbatch))
            if not seqbatch: raise IOError('No batch input fasta files found! Make sure parentN=FILE settings given *.fofn.')
            ## ~ [0e] Setup subread lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elists = [udb1.sortedEntries('Len',reverse=True),udb2.sortedEntries('Len',reverse=True)]
            plen = [0,0]    # Summed lengths for each parent
            pseq = [0,0]    # Total sequence number for each parent
            prq = [0,0]     # Total sequence RQ for each parent (convert to mean)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            lastlen = max(elists[0][0]['Len'],elists[1][0]['Len'])    # Length of last selected read
            for elist in elists:
                while elist and elist[0]['RQ'] < minrq: elist.pop(0)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            nextp = 0       # Index of next parent to use
            if elists[0][0]['Len'] < elists[1][0]['Len']: nextp = 1

            ### ~ [1] Filter and Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Filter Unique Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            zmwlist = []    # List of (smrt,zmw) meeting filtering criteria
            ux = 0.0; utot = len(elists[0])+len(elists[1])
            while lastlen:
                self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq: elist.pop(0); ux += 100.0
                if elist and elist[0]['Len'] < minlen: ux += 100.0 * len(elist); elist = []
                if not elist: nextp = 1 - nextp; break  # Finish
                entry = elist.pop(0); ux += 100.0
                zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                plen[nextp] += entry['Len']
                prq[nextp] += entry['RQ']
                pseq[nextp] += 1
                if plen[1-nextp] <= plen[nextp]: nextp = 1 - nextp
                lastlen = entry['Len']
            ## ~ [1b] Final processing of last reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while elists[nextp]:
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    elist.pop(0); ux += 100.0
                while elist and elist[0]['Len'] >= minlen:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    entry = elist.pop(0); ux += 100.0
                    pdiff = rje.modulus(plen[0]-plen[1])
                    ediff = rje.modulus(plen[nextp]+entry['Len']-plen[1-nextp])
                    if ediff >= pdiff: elists[nextp] = []; break    #Finish!
                    zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                    plen[nextp] += entry['Len']
                    prq[nextp] += entry['RQ']
                    pseq[nextp] += 1
            self.printLog('\r#DIP','Diploidising subreads complete: %s subreads to output.' % rje.iLen(zmwlist))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent1'),rje.iStr(pseq[0]),rje.iStr(plen[0]),1.0*plen[0]/self.getInt('GenomeSize'),prq[0]/pseq[0]))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent2'),rje.iStr(pseq[1]),rje.iStr(plen[1]),1.0*plen[1]/self.getInt('GenomeSize'),prq[1]/pseq[1]))
            ## ~ [1b] Extract Filtered Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,filtfile)
            SEQOUT = open(filtfile,'w')
            sx = 0.0; stot = 0; sn = len(seqbatch); fx = 0
            for seqlist in seqbatch:
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/9/0_3967 RQ=0.784
                si = 100.0/seqlist.seqNum(); stot += seqlist.seqNum()
                for seq in seqlist.seqs():
                    self.progLog('\r#OUT','Extracting subreads: %.2f%%' % (sx/sn)); sx += si
                    (name,sequence) = seqlist.getSeq(seq)
                    try: [smrt,zmw,pos,rq] = string.split(string.replace(name,'/',' '))
                    except:
                        [smrt,zmw,pos] = string.split(string.replace(name,'/',' '))
                        rq = minrq
                    if (cdb.data(smrt)['SMRT'],int(zmw),pos) not in zmwlist: continue
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fx += 1
            self.printLog('\r#OUT','Saved %s filtered subreads to %s.' % (rje.iStr(fx),filtfile))

            ### ~ [2] Summarise Filtered File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % filtfile,'autofilter=F']
            rje_seqlist.SeqList(self.log,seqcmd)

            return True
        except: self.errorLog('%s.run error' % self.prog()); return False
#########################################################################################################################
    ### <4> ### SortNR & DipHap run methods                                                                             #
#########################################################################################################################
    def sortNR(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            debugstr = 'NOTDEBUGGING'
            if self.getStrLC('DebugStr'): debugstr = self.getStr('DebugStr')
            ## ~ [1a] Input Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F','autofilter=F'])
            tmpseq = '%s.tmp.fasta' % self.baseFile()
            TMPSEQ = open(tmpseq,'w')
            #i# seqdict will store all the sequences, mapped by shortname onto PAF results
            #i# As sequences are filtered, they will be removed from seqdict and added to badseq. Remaining sequences will be saved.
            seqdict = seqin.seqNameDic()
            seqin.dict['Filter']['NullSeq'] = 0
            seqin.nextSeq()
            minlen = seqin.seqNonX()
            sx = 0.0; stot = seqin.seqNum()
            while seqin.currSeq():
                self.progLog('\r#LEN','Scanning sequences: %.1f%%' % (sx/stot)); sx += 100
                if seqin.seqNonX():
                    (name,sequence) = seqin.currSeq()
                    sequence = string.join( re.split('[Nn]{10}[Nn]+',sequence), 'NNNNNNNNNN')
                    TMPSEQ.write('>%s\n%s\n' % (name,sequence))
                    minlen = min(minlen,len(sequence))
                else:
                    seqin.dict['Filter']['NullSeq'] += 1
                    seqdict.pop(seqin.shortName())
                if not seqin.nextSeq(): break
            TMPSEQ.close()
            self.printLog('\r#NULL','%s 100%% N sequences filtered.' % (rje.iStr(seqin.dict['Filter']['NullSeq'])))
            self.printLog('#TMP','%s gap-reduced sequences output to temporary file %s for minimap2 search' % (rje.iLen(seqdict),tmpseq))
            self.debug(minlen)

            #!# There is an issue of long gaps causing problems, so need to replace all gaps with 10 Ns prior to PAF
            #!# generation, and make sure that files are subsequently cleaned up (unless debug=T or dev=T).

            ## ~ [1b] MiniMap PAF generation and parsing object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pafin = '%s.tmp.paf' % self.baseFile()
            self.obj['PAF'] = paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'basefile=%s' % self.baseFile(),'seqin=%s' % tmpseq,'reference=%s' % tmpseq,'sortseq=None'])
            paf.obj['DB'] = self.obj['DB']
            paf.setNum({'MinLocID':100})
            paf.setInt({'MinLocLen':int(minlen*self.getPerc('CheckCov'))})
            paf.dict['MapOpt'] = rje.combineDict(paf_defaults,paf.dict['MapOpt'],overwrite=True)
            paf.setup()
            if not rje.exists(pafin) or self.force():
                rje.backup(self,pafin)
                #!# Add use of qsub #!#
                paf.minimap2()
            # Parse PAF with initial minloclen and minlocid filtering
            pafdb = paf.parsePAF(debugstr=debugstr)
            pafdb.newKey(['Hit','Qry','#'])

            ### ~ [2] NR Sequence filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nrfile = '%s.redundant.txt' % self.baseFile()
            rje.backup(self,nrfile)
            NRFILE = open(nrfile,'w')
            #i# pafhead = Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality .. cs
            seqin.dict['Filter']['NR'] = 0
            px = 0.0; ptot = pafdb.entryNum()
            #i# Work through PAF table starting with longest hits and investigate status of queries
            for pentry in pafdb.sortedEntries('SbjLen',reverse=True):
                self.progLog('\r#NR','Filtering redundant sequences from PAF CS strings: %.2f%%' % (px/ptot)); px += 100
                qry = pentry['Qry']
                hit = pentry['Hit']
                if qry.endswith(debugstr): self.bugPrint(pentry)
                #i# Ignore self-hits
                if qry == hit: continue
                #i# Skip any hits already filtered
                if qry not in seqdict: continue
                if hit not in seqdict: continue
                cstats = paf.statsFromCS(pentry['cs'])
                if qry.endswith(debugstr): self.bugPrint(cstats)
                #!# Might need to use an actual value rather than a percentage!
                if self.getPerc('CheckCov') * pentry['QryLen'] <= cstats[':'] < pentry['QryLen']:  # Hit might be redundant with query
                    qryseq = seqin.getSeq(seqdict[qry])[1]
                    hitseq = seqin.getSeq(seqdict[hit])[1]
                    if pentry['Strand'] == '-': qryseq = rje_sequence.reverseComplement(qryseq)
                    if qryseq in hitseq: cstats[':'] = pentry['QryLen']
                    if qry.endswith(debugstr): self.bugPrint('\n%s\n vs \n%s\n = %s' % (qryseq,hitseq,qryseq in hitseq))
                if cstats[':'] == pentry['QryLen']:  # Hit is redundant with query
                    rtxt = qry
                    if pentry['Strand'] == '-': rtxt = rtxt + ' (RevComp)'
                    if cstats[':'] == pentry['SbjLen']: rtxt += ' = identical to '
                    else: rtxt += ' = contained within '
                    rtxt += hit
                    NRFILE.write('%s\n' % rtxt)
                    seqin.dict['Filter']['NR'] += 1
                    seqdict.pop(qry)
                if qry.endswith(debugstr): self.debug('%s in seqdict?: %s' % (qry,qry in seqdict))
            NRFILE.close()
            self.printLog('\r#NR','%s redundant sequences filtered from PAF CS strings: see %s.' % (rje.iStr(seqin.dict['Filter']['NR']),nrfile))

            ### ~ [3] Save filtered sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodseq = seqdict.values()
            goodseq.sort()
            seqin.list['Seq'] = goodseq
            if not seqin.getStrLC('SeqOut'):
                seqin.setStr({'SeqOut':'%s.nr.fasta' % self.baseFile()})
            seqfiles = [seqin.getStr('SeqIn'),seqin.getStr('SeqOut')]
            seqin.saveSeq()
            seqin.setStr({'SeqIn':seqin.getStr('SeqOut')})
            ## ~ [3a] Option to delete PAF files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for tmp in ['paf','paf.cmd','fasta','fasta.index']:
                tmpfile = '%s.tmp.%s' % (self.baseFile(),tmp)
                if rje.exists(tmpfile) and (not (self.dev() or self.debugging()) or (self.i() > 0 and rje.yesNo('Delete %s?' % tmpfile))):
                    os.unlink(tmpfile)
                    self.printLog('#CLEAN','%s deleted' % tmpfile)
            ## ~ [3b] Summarise sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Summarise'): rje_seqlist.batchSummarise(self,seqfiles,save=True,overwrite=self.force())
            return seqin     # SortNR successful
        except: self.errorLog('Problem during %s.sortNR.' % self.prog()); return False
#########################################################################################################################
    def dipHap(self):   ### Update sequence names with haploid[AB] and diploid, then output along with pri and alt files.
        '''Update sequence names.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Setup input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            overwrite = self.force()
            seqfiles = []
            if self.getStrLC('RunMode') == 'diphapnr':
                seqlist = self.sortNR()
                overwrite = False
            else:
                seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
                seqfiles = [seqlist.getStr('SeqIn')]
            ## ~ [0b] Setup output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            SEQOUT = open(self.baseFile() + '.dipnr.fasta','w')
            PRIOUT = open(self.baseFile() + '.pri.fasta','w')
            ALTOUT = open(self.baseFile() + '.alt.fasta','w')
            seqfiles += [self.baseFile() + '.dipnr.fasta',self.baseFile() + '.pri.fasta',self.baseFile() + '.alt.fasta']
            ### ~ [1] Allocate primary and alternative scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            diplist = []
            haplist = []
            dx = px = ax = sx = 0.0; stot = seqlist.seqNum() * 2
            for seq in seqlist.seqs():
                self.progLog('\r#DIPHAP','Pseudodiploid haplotig assignment: %.1f%%' % (sx/stot)); sx += 100
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in diplist: haplist.append(hap)
                else: diplist.append(hap)
            for seq in seqlist.seqs():
                self.progLog('\r#DIPHAP','Pseudodiploid haplotig assignment: %.1f%%' % (sx/stot)); sx += 100
                sname = string.split(seqlist.shortName(seq),'_')
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in haplist:
                    if hap in diplist: haptxt = 'haploidA'; diplist.remove(hap); sname[0] = 'pri%s' % hap; px += 1
                    else: haptxt = 'haploidB'; sname[0] = 'alt%s' % hap; ax += 1
                else: haptxt = 'diploid'; sname[0] = 'pri%s' % hap; dx += 1; px += 1
                sname = string.join(sname,'_')
                SEQOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                if sname[:3] == 'pri':
                    PRIOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                else:
                    ALTOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                #self.printLog('#DIP','%s: %s\n' % (sname,haptxt))
            self.printLog('\r#DIPHAP','Pseudodiploid haplotig assignment complete!')
            SEQOUT.close()
            PRIOUT.close()
            ALTOUT.close()
            self.printLog('#OUT','%s Primary assembly haplotigs (inc. %s Diploid) output to %s.pri.fasta' % (rje.iStr(px),rje.iStr(dx),self.baseFile()))
            self.printLog('#OUT','%s Alternative assembly haplotigs output to %s.alt.fasta' % (rje.iStr(ax),self.baseFile()))
            self.printLog('#OUT','%s Combined Primary+Alternative assembly haplotigs output to %s.dipnr.fasta' % (rje.iStr(ax+px),self.baseFile()))
            ### ~ [2] Summarise sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Summarise'): rje_seqlist.batchSummarise(self,seqfiles,save=True,overwrite=overwrite)
            return True
        except: self.errorLog('Problem during %s.dipHap.' % self.prog()); return False
#########################################################################################################################
    ### <5> ### VecScreen methods                                                                                       #
#########################################################################################################################
    def vecScreen(self,makedb=False):    ### Screen Scaffolds for possible contaminants
        '''
        ## Vector contamination screening [runmode=vecscreen]

        This mode screens scaffolds for possible contaminants, given by `screendb=FILE`.

        ### VecScreen overview

        First, a `blastn` search of the `screendb=FILE` sequences is performed, using the NCBI VecScreen search and
        match strategy (see below). These parameters can be modified using `blastopt=X` and/or `optionfile=FILE`,
        which will be appended to the run command. Alternatively, the `$BASEFILE.vecscreen.blast` file produced can be
        pre-made with different parameters (if `force=F`).

        Results are then parsed into a local hits table and rated according to the strength of a match. This is performed
        iteratively, re-assigning internal matches as "proximal" if they are withing 25 bases of another match (of any
        strength) and re-rating the strength of the match, until no ratings changes occur. Two additional fields are
        added to the local hits table during this process: `MatchPos` and `MatchStr`. Once all assignments have been
        made, segments of the assembly between two matches and/or sequence ends are added to the table as `Suspect`.

        `MatchPos` will have a value of:

        * `Terminal` = within 25 bp of either end of the sequence.
        * `Proximal` = within 25 bp of a vecsreen match (`Weak`, `Moderate` or `Strong`).
        * `Internal` = over 25 bp of a sequence end or vecsreen match.
        * `Suspect` = Segments added as `Suspect`.

        `MatchStr` will have a value of:

        * `Strong` = `Terminal`/`Proximal` match with Score >= 24, or `Internal` match with Score >= 30.
        * `Moderate` = `Terminal`/`Proximal` match with Score 19 to 23, or `Internal` match with Score 25 to 29.
        * `Weak` = `Terminal`/`Proximal` match with Score 16 to 18, or `Internal` match with Score 23 to 24.
        * `Suspect` = Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        These parameters seem to be designed for small (Illumina?) assemblies, and give excessive false positives for
        larger genomes. Diploidocus therefore features two additional contaminant identification filters:

        * `eFDR` = The "Expected False Discovery Rate` is calculated for each contaminant sequence. This is simply the
        `Expect` value for that hit, divided by the total number of hits with that `Expect` value (or lower). If `efdr=X`
        is set, any hits with an `eFDR` value exceeding the threshold will be removed. By default, this is set to `1.0`,
        i.e. any contaminants with fewer assembly hits than their `Expect` value will be dropped.
        * `Minimum hit length` = Short hits are inevitable in large assemblies and unlikely to be real for long read
        assemblies (e.g. without cloning etc.). An additional `minvechit=X` setting will remove any really short hits.
        This is set by default to `minvechit=50`, meaning that a hit of at least 50 bp is required.

        Finally, the percentage coverage per scaffold is calculated from the filtered hits. This is performed first for
        each contaminant individually, before being collapsed into total contamination coverage per query.

        Results are output into two main delimited results files:

        * `*.vecscreen.tdt` = Contaminant local hit table with VecScreen ratings and eFDR calculation. (Unfiltered)
        * `*.screencov.tdt` = Query coverage per contaminant
            - `Query` = Contaminant name. For total coverage, this will be `TOTAL`.
            - `Hit` = Scaffold name.
            - `HitLen` = Length of scaffold (bp).
            - `Coverage` = Covered length of scaffold (bp).
            - `CovPC` = Percentage coverage of scaffold (0-100).

        ### VecScreen BLAST+ parameters (from NCBI website)

        The VecScreen BLAST+ parameters are pre-set using blastn options: -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000

        VecScreen Match Categories
        Vector contamination usually occurs at the beginning or end of a sequence; therefore, different criteria are
        applied for terminal and internal matches. VecScreen considers a match to be terminal if it starts within 25
        bases of the beginning of the query sequence or stops within 25 bases of the end of the sequence. Matches that
        start or stop within 25 bases of another match are also treated like terminal matches. Matches are categorized
        according to the expected frequency of an alignment with the same score occurring between random sequences.

        Strong Match to Vector
        (Expect 1 random match in 1,000,000 queries of length 350 kb.)
        Terminal match with Score >= 24.
        Internal match with Score >= 30.

        Moderate Match to Vector
        (Expect 1 random match in 1,000 queries of length 350 kb.)
        Terminal match with Score 19 to 23.
        Internal match with Score 25 to 29.

        Weak Match to Vector
        (Expect 1 random match in 40 queries of length 350 kb.)
        Terminal match with Score 16 to 18.
        Internal match with Score 23 to 24.

        Segment of Suspect Origin
        Any segment of fewer than 50 bases between two vector matches or between a match and an end.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS VECTOR SCREENING',line='-')
            db = self.db()
            #i#covdb = self.db().addEmptyTable('screencov',['Query','Hit','HitLen','Coverage','CovPC'],['Query','Hit'],log=self.debugging())
            screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and not makedb and rje.checkForFiles(filelist=[screencov],basename='',log=self.log):
                return screencov
            forks = self.getInt('Forks')
            bfile = '%s.vecscreen.blast' % self.baseFile()
            locfile = '%s.vecscreen.Local.tdt' % self.baseFile()
            if not rje.exists(self.getStr('SeqIn')):
                raise IOError('Diploidocus VecScreen mode needs input assembly (seqin=FILE)')
            if not rje.exists(bfile) and not rje.exists(self.getStr('ScreenDB')):
                raise IOError('Diploidocus VecScreen mode needs screendb=FILE set or exisiting %s file' % bfile)
            blast = rje_blast.BLASTRun(log=self.log,cmd_list=['blaste=700','keepblast=T','blastgz=T']+self.cmd_list+['bitscore=F'])
            self.obj['BLAST'] = blast
            mytables = self.db().tables()
            blast.obj['DB'] = self.db()
            blast.setStr({'Type':'blastn','InFile':self.getStr('ScreenDB'),'BlastRes':bfile,'Name':bfile,'DBase':self.getStr('SeqIn')})
            blast.setup(load=False)
            ## ~ [1a] ~ Check existing results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            runblast = self.force() or not rje.exists(locfile)
            if rje.exists(locfile):
                if self.force(): self.printLog('#FORCE','%s found but will be regenerated. (force=T)' % locfile)
                else: self.printLog('#LOCAL','%s found - BLAST will be skipped. (force=F)' % locfile)
            if runblast and not self.force() and blast.checkBLAST():                ### BLAST Results exist
                if blast.getBool('IgnoreDate'): runblast = False       ### Don't check age!
                elif rje.isYounger(blast.getStr('DBase'),blast.getStr('Name')) == blast.getStr('Name') and rje.isYounger(blast.getStr('InFile'),blast.getStr('Name')) == blast.getStr('Name'): runblast = False
            ## ~ [1b] ~ Setup search database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if runblast:
                blast.formatDB(protein=False,force=self.force(),details=self.debugging())
            ## ~ [1c] Input Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F'])
            #i# seqdict will store all the sequences, mapped by shortname
            seqdict = seqin.seqNameDic()

            ### ~ [2] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First, a `blastn` search of the `screendb=FILE` sequences is performed, using the NCBI VecScreen search and
            # match strategy (see below). These parameters can be modified using `blastopt=X` and/or `optionfile=FILE`,
            # which will be appended to the run command. Alternatively, the `$BASEFILE.vecscreen.blast` file produced can be
            # pre-made with different parameters (if `force=F`).
            if runblast:
                veccmd = ' -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -searchsp 1750000000000'
                command = blast.blastPath() + blast.getStr('Type') + veccmd
                command += ' -evalue %e' % blast.getNum('E-Value')
                command += ' -query %s' % blast.getStr('InFile')
                command += ' -db %s' % blast.getStr('DBase')
                command += ' -out %s' % blast.getStr('Name')
                if forks > 1: command += ' -num_threads %d' % forks
                if blast.getStrLC('BLASTOpt'): command = '%s %s' % (command,blast.getStr('BLASTOpt'))
                if rje.exists(blast.getStr('OptionFile')):
                    for line in open(blast.getStr('OptionFile'),'r').readlines(): command = '%s %s' % (command,rje.chomp(line))
                blast.str['BLASTCmd'] = command
                #self.printLog('\r#SYS',command)
                self.loggedSysCall(command,syslog='{}.blastn.log'.format(self.baseFile()))
                #os.system(command)
                if not blast.checkBLAST(): raise IOError('Problem with VecSreen BLAST results file "%s"' % bfile)
            ## ~ [2a] Read Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vecdb = None
            if rje.exists(locfile) and not self.force():
                blast.db().list['Tables'] = mytables
                vecdb = self.db().addTable(locfile,['Query','Hit','AlnID'],name='vecscreen')
                blast.formatTables()
                vecdb.dataFormat({'QRank':'int'})
            else:
                keepaln = False     #!# Might want to add alignment output later
                rmblast = not blast.getBool('KeepBLAST')
                #self.debug('Delete BLAST: %s' % rmblast)
                blast.readBLAST(resfile=blast.getStr('BlastRes'),clear=True,gablam=False,local=True,keepaln=keepaln,unlink=rmblast)
                blast.formatTables()
            ## ~ [2b] Add eFDR calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                vecdb = self.db('Local')    #?# Should this be a copy?
                #i# Keys: ['Query','Hit','AlnID']
                #i# Fields: ['Query','Hit','AlnID','Score','Expect','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd']
                vecdb.keepFields(['Query','Hit','AlnID','Score','Expect','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd'])
                self.progLog('#EFDR','Calculating expected FDR (eFDR) based on Expect scores and observed hits...')
                vecdb.rankFieldByIndex('Query','Expect',newfield='QRank',rev=False,absolute=True,lowest=False,warn=True,highest=True)
                vecdb.makeField('Expect/QRank','eFDR')
                self.printLog('\r#EFDR','Calculation of expected FDR (eFDR) based on Expect scores and observed hits complete.')
                vecdb.saveToFile(locfile)
                vecdb.setStr({'Name':'vecscreen'})
            ## ~ [2c] Filter Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vecdb.dataFormat({'eFDR':'float'})
            if self.getNum('eFDR') > 0:
                #self.printLog('#EFDR','Filtering entries with eFDR<{}'.format(self.getNum('eFDR')))
                #vecdb.dropEntries('eFDR<{}'.format(self.getNum('eFDR')))
                ex = 0.0; etot = vecdb.entryNum()
                for ekey, data in vecdb.data().items():
                    self.progLog('\r#EFDR','Filtering entries with eFDR>{}: {:.2f}%'.format(self.getNum('eFDR'),ex/etot)); ex += 100.0
                    if data['eFDR'] > self.getNum('eFDR'): vecdb.data().pop(ekey)
                    #else: self.debug(data)
                self.printLog('\r#EFDR','Filtered entries with eFDR<{}: {} -> {} entries'.format(self.getNum('eFDR'),rje.iStr(etot),rje.iStr(vecdb.entryNum())))
            if self.getNum('MinVecHit') > 0 or self.getNum('MinIDHit') > 0:
                self.printLog('#HITLEN','Filtering identical hits with Length<{0} bp and non-identical hits <{1} bp'.format(self.getInt('MinIDHit'),self.getInt('MinVecHit')))
                idx = 0
                vecx = 0
                prex = vecdb.entryNum()
                for vkey, ventry in vecdb.data().items():
                    if ventry['Length'] < self.getInt('MinIDHit'):
                        vecdb.dict['Data'].pop(vkey); idx += 1
                    elif ventry['Identity'] != ventry['Length'] and ventry['Length'] < self.getInt('MinVecHit'):
                        vecdb.dict['Data'].pop(vkey); vecx += 1
                if idx or vecx: vecdb.dict['Index'] = {}
                self.printLog('#DROP','Length filtering: %s vecscreen entries reduced to %s entries' % (rje.integerString(prex),rje.integerString(vecdb.entryNum())))
                #x#vecdb.dropEntries('Length<{}'.format(self.getInt('MinVecHit')))

            ### ~ [3] Process Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            proxbp = 25
            suspectbp = 50
            # Results are then parsed into a local hits table and rated according to the strength of a match. This is performed
            # iteratively, re-assigning internal matches as "proximal" if they are withing 25 bases of another match (of any
            # strength) and re-rating the strength of the match, until no ratings changes occur. Two additional fields are
            # added to the local hits table during this process: `MatchPos` and `MatchStr`. Once all assignments have been
            # made, segments of the assembly between two matches and/or sequence ends are added to the table as `Suspect`.
            #
            # `MatchPos` will have a value of:
            #
            # * `Terminal` = within 25 bp of either end of the sequence.
            # * `Proximal` = within 25 bp of a vecsreen match (`Weak`, `Moderate` or `Strong`).
            # * `Internal` = over 25 bp of a sequence end or vecsreen match.
            # * `Suspect` = Segments added as `Suspect`.
            #
            # `MatchStr` will have a value of:
            #
            # * `Strong` = `Terminal`/`Proximal` match with Score >= 24, or `Internal` match with Score >= 30.
            # * `Moderate` = `Terminal`/`Proximal` match with Score 19 to 23, or `Internal` match with Score 25 to 29.
            # * `Weak` = `Terminal`/`Proximal` match with Score 16 to 18, or `Internal` match with Score 23 to 24.
            # * `Suspect` = Any segment of fewer than 50 bases between two vector matches or between a match and an end.
            ## ~ [3a] Setup table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Add Strand and make start/end consistent
            vecdb.addField('Strand',evalue='+')
            for ventry in vecdb.entries():
                if ventry['SbjStart'] > ventry['SbjEnd']: (ventry['SbjStart'],ventry['SbjEnd'],ventry['Strand']) = (ventry['SbjEnd'],ventry['SbjStart'],'-')
            ## Add MatchPos and MatchType
            vecdb.addField('MatchPos',evalue='Internal')
            vecdb.addField('MatchStr',evalue='None')
            #i# MatchStart and MatchEnd have the Subject Start/End if a Match
            vecdb.addField('MatchStart',evalue=0)
            vecdb.addField('MatchEnd',evalue=0)
            #i# Terminal and Internal give the match ratings depending on position
            vecdb.addField('Terminal',evalue='None')
            vecdb.addField('Internal',evalue='None')
            ## ~ [3b] Filter < Weak matches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ventry in vecdb.entries():
                if ventry['Score'] >= 30: ventry['Internal'] = 'Strong'
                elif ventry['Score'] >= 25: ventry['Internal'] = 'Moderate'
                elif ventry['Score'] >= 23: ventry['Internal'] = 'Weak'
                if ventry['Score'] >= 24: ventry['Terminal'] = 'Strong'
                elif ventry['Score'] >= 19: ventry['Terminal'] = 'Moderate'
                elif ventry['Score'] >= 16: ventry['Terminal'] = 'Weak'
            vecdb.dropEntriesDirect('Terminal',['None'])
            ## ~ [3c] Re-order based on Hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# New key is for determining proximity to other matches.
            vecdb.newKey(['Hit','SbjStart','SbjEnd','Query','AlnID'])
            ## ~ [3d] Add Terminal ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x# vecdb.indexReport('Hit')
            for sname in vecdb.index('Hit'):
                hitlen = seqin.seqLen(seqdict[sname])
                for ventry in vecdb.indexEntries('Hit',sname):
                    if ventry['SbjStart'] <= proxbp: ventry['MatchPos'] = 'Terminal'
                    elif hitlen - ventry['SbjEnd'] < proxbp: ventry['MatchPos'] = 'Terminal'
                    if ventry['MatchPos'] == 'Terminal': ventry['MatchStr'] = ventry['Terminal']
                    else: ventry['MatchStr'] = ventry['Internal']
                    if ventry['MatchStr'] != 'None': (ventry['MatchStart'],ventry['MatchEnd']) = (ventry['SbjStart'],ventry['SbjEnd'])
            ## ~ [3e] Cycle and rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cycling = True
            while cycling:
                cycling = False     #i# Set to True if ratings updated.
                sx = 0.0; stot = len(vecdb.index('Hit'))
                for sname in vecdb.index('Hit'):
                    vhits = vecdb.indexEntries('Hit',sname)[0:]
                    sstep = 100.0 / len(vhits)
                    for ventry in vhits[0:]:
                        self.progLog('\r#MATCH','Updating proximal matches: %.1f%%' % (sx/stot)); sx += sstep
                        if ventry['MatchPos'] in ['Terminal','Proximal']: continue
                        #self.bugPrint(vecdb.entrySummary(ventry))
                        #i# Identify a match within the required distance #i#
                        proximal = False
                        proxstart = max(1,ventry['SbjStart'] - proxbp + 1)
                        proxend = ventry['SbjEnd'] + proxbp - 1
                        for vhit in vhits[0:]:
                            if vhit['MatchEnd'] < proxstart: vhits.remove(vhit)
                            if vhit['MatchStart'] <= proxend:
                                #self.debug(vecdb.entrySummary(vhit,collapse=True))
                                proximal = True
                                break
                        if proximal:
                            ventry['MatchPos'] = 'Proximal'
                            ventry['MatchStr'] = ventry['Terminal']
                            (ventry['MatchStart'],ventry['MatchEnd']) = (ventry['SbjStart'],ventry['SbjEnd'])
                            cycling = True
                if cycling: self.printLog('\r#MATCH','Updating proximal matches: still cycling.')
                else: self.printLog('\r#MATCH','Updating proximal matches complete.')
            vecdb.dropEntriesDirect('MatchStr',['None'])
            ## ~ [3g] Add suspect regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0.0; stot = len(vecdb.index('Hit')); suspectx = 0
            for sname in vecdb.index('Hit'):
                self.progLog('\r#REGION','Identifying suspect regions: %.1f%%' % (sx/stot)); sx += 100.0
                hitlen = seqin.seqLen(seqdict[sname])
                matches = []    # (start,end) tuple list
                vhits = vecdb.indexEntries('Hit',sname)[0:]
                for ventry in vhits[0:]:
                    matches.append((ventry['MatchStart'],ventry['MatchEnd']))
                matches.sort()
                matches = rje.collapseTupleList(matches)
                suspects = rje.invertTupleList(matches,minx=1,maxx=hitlen)
                for (starti,endi) in suspects:
                    if (endi - starti) < suspectbp:
                        vecdb.addEntry({'Query':'Suspect','Hit':sname,'AlnID':0,'SbjStart':starti,'SbjEnd':endi,
                                        'Score':0,'Identity':0,'Length':endi - starti + 1,'Strand':'+',
                                        'MatchPos':'Suspect','MatchStr':'Suspect'})
                        suspectx += 1
            self.printLog('\r#REGION','Identifying suspect regions complete: %s regions' % rje.iStr(suspectx))

            ### ~ [4] Save Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vecdb.dropFields(['MatchStart','MatchEnd','Terminal','Internal'])
            vecdb.saveToFile()
            ## ~ [4a] Optional VecCheck output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('VecCheck'):
                vfile = db.dbFileName('vecscreen')
                pfile = self.baseFile() + '.paf'
                pafcmd = self.cmd_list + ['checkpos={}'.format(vfile),'checkfields=Hit,SbjStart,SbjEnd','pafin={}'.format(pfile)]
                paf = rje_paf.PAF(self.log, pafcmd)
                vecdb = paf.checkPos(save=False)
                vecdb.setStr({'Name':'vecscreen'})
                vecdb.saveToFile(backup=False)

            ### ~ [5] Final Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Need to look at numbers per contaminant to help detect false positives.
            vecdb.indexReport('Query')
            ## ~ [5a] Calculate percentage coverage per query/hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            covdb = self.db().addEmptyTable('screencov',['Query','Hit','HitLen','Coverage','CovPC'],['Query','Hit'],log=self.debugging())
            qhcov = {}  # Dictionary of { (Qry,Hit) : [(hitstart,hitend) , (...)]
            vx = 0.0; vtot = vecdb.entryNum()
            for ventry in vecdb.entries():
                self.progLog('\r#COV','Calculating percentage coverage per query/hit pair: {:.2f}%'.format(vx/vtot)); vx += 50.0
                qry = ventry['Query']
                hit = ventry['Hit']
                qh = (qry,hit)
                if qh not in qhcov: qhcov[qh] = []
                qhcov[qh].append((ventry['SbjStart'],ventry['SbjEnd']))
            qx = 0.0; qtot = len(qhcov)
            for qh, covlist in qhcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per query/hit pair: {:.2f}%'.format(50.0+(qx/qtot))); qx += 50.0
                covsum = 0
                for (ci,cj) in rje.collapseTupleList(covlist,joindistance=1,overlaps=False):
                    covsum += cj - ci + 1
                hlen = seqin.seqLen(seqdict[qh[1]])
                covdb.addEntry({'Query':qh[0], 'Hit':qh[1], 'HitLen':hlen, 'Coverage':covsum, 'CovPC':rje.dp(100.0 * covsum /  hlen,2)})
            self.printLog('\r#COV','Calculation of percentage coverage per query/hit pair complete')
            ## ~ [5b] Calculate total contamination coverage for each Hit ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hcov = {}
            qx = 0.0; qtot = len(qhcov)
            for qh, covlist in qhcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per scaffold: {:.2f}%'.format((qx/qtot))); qx += 50.0
                if qh[1] not in hcov: hcov[qh[1]] = []
                hcov[qh[1]] += covlist
            qx = 0.0; qtot = len(hcov)
            for hit, covlist in hcov.items():
                self.progLog('\r#COV','Calculating percentage coverage per scaffold: {:.2f}%'.format(50.0+(qx/qtot))); qx += 50.0
                covsum = 0
                for (ci,cj) in rje.collapseTupleList(covlist,joindistance=1,overlaps=False):
                    covsum += cj - ci + 1
                hlen = seqin.seqLen(seqdict[hit])
                covdb.addEntry({'Query':'TOTAL', 'Hit':hit, 'HitLen':hlen, 'Coverage':covsum, 'CovPC':rje.dp(100.0 * covsum /  seqin.seqLen(seqdict[hit]),2)})
            self.printLog('\r#COV','Calculation of percentage coverage per scaffold complete')
            covdb.saveToFile()
            return screencov

        except IOError:
            self.errorLog('Diploidocus.vecSreen() error'); raise
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def vecPurge(self): ###
        '''
        ### VecScreen purge mode

        When `screenmode=purge`, sequences will be removed, trimmed and/or masked based on identified contamination
        sequences. This is performed after `minvechit=INT` and `efdr=NUM` filtering. If Diploidocus is re-run (without
        `force=T`), existing `BLAST+` results and the local hits table will be re-used, but the `*.vescreend.tdt` and
        `*.screencov.tdt` output will be regenerated prior to purging. This enables purging to be re-run with different
        parameters.

        The three main purge parameter are:

            vecpurge=PERC   : Remove any scaffolds with >= PERC % vector coverage [50.0]
            vectrim=INT     : Trim any vector hits (after any vecpurge) within INT bp of the nearest end of a scaffold [1000]
            vecmask=INT     : Mask any vectore hits of INT bp or greater (after vecpurge and vecmerge) [1000]

        These are applied in order:

        1. First, any scaffolds meeting the vecpurge=PERC total coverage threshold are removed completely.
        2. Next, any hits within `vectrim=INT` bp of either sequence end are trimmed off. This will continue until no hits
        (meeting `minvechit=INT efdr=NUM`) are within `vectrim=INT` bp. If both ends are within `vectrim=INT`, the whole
        sequence will be trimmed and the sequence will be removed.
        3. Finally, any remaining hits meeting the `vecmask=INT` length cutoff are masked by replacing with Ns.

        Output is then saved to:

        * `*.vecscreen.fasta` = The purged, trimmed and masked sequences. If no sequences were identified for editing,
        this will be the same as the input sequences.
        * `*.vecpurge.tdt` = The table of vecsreen edits
            - `SeqName`, `SeqLen`, `Start`, `End`, `Edit`

        :return: trimdb or None if no purging/trimming/masking
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.seqNameDic()
            basefile = self.baseFile(strip_path=True)
            if self.getStrLC('MaskMode') not in ['full','partial']:
                self.warnLog('#MASK','MaskMode "{0}" not recognised: defaulting to "partial".'.format(self.getStrLC('MaskMode')))
                self.setStr({'MaskMode':'partial'})
            maskmode = self.getStrLC('MaskMode')
            self.printLog('#MASK','Mask mode: {0}'.format(maskmode))

            ## ~ [1a] ~ Setup database tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            vecdb = self.db('vecscreen')
            covdb = self.db('screencov')
            ## ~ [1b] ~ Generate database tables if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Re-run vecScreen() rather than load data in case filtering settings have changed.
            #i# The BLAST results or Local hits table will be reloaded unless force=T.
            if not (vecdb and covdb):
                self.vecScreen(makedb=True)
            vecdb = self.db('vecscreen')
            covdb = self.db('screencov')    # ['Query','Hit','HitLen','Coverage','CovPC']
            if not (vecdb and covdb): raise IOError('Problem generating vecscreen and screencov tables')
            covdb.dataFormat({'CovPC':'num'})

            ### ~ [2] ~ Perform purging, trmming and masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            trimdb = db.addEmptyTable('vecpurge',['SeqName','SeqLen','Start','End','Edit'],['SeqName','Start','End'],log=self.debugging())
            ## ~ [2a] ~ Purge scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# 1. First, any scaffolds meeting the vecpurge=PERC total coverage threshold are removed completely.
            purgelist = []
            for entry in covdb.indexEntries('Query','TOTAL'):
                if self.getNum('VecPurge') > 0 and entry['CovPC'] >= self.getNum('VecPurge'):
                    seqname = entry['Hit']
                    seqlen = entry['HitLen']
                    trimdb.addEntry({'SeqName':seqname, 'SeqLen':seqlen, 'Start':1, 'End':seqlen, 'Edit':'purge'})
                    purgelist.append(seqname)
            ## ~ [2b] ~ Generate collapsed [(start,end)] tuple lists for remaining sequences ~~~~~~ ##
            for seqname in vecdb.index('Hit'):
                if seqname in purgelist: continue
                seqlen = seqlist.seqLen(seqdict[seqname])
                vecpos = []
                for entry in vecdb.indexEntries('Hit',seqname): vecpos.append((entry['SbjStart'],entry['SbjEnd']))
                vecpos = rje.collapseTupleList(vecpos)
                trimpos = rje.collapseTupleList(vecpos,joindistance=self.getInt('VecTrim'))
            ## ~ [2c] ~ Identify regions to trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# 2. Next, any hits within `vectrim=INT` bp of either sequence end are trimmed off. This will continue until no hits
            #i# (meeting `minvechit=INT efdr=NUM`) are within `vectrim=INT` bp.
                trim5 = 1; trim3 = seqlen
                if self.getInt('VecTrim') > 0 and trimpos[0][0] <= self.getInt('VecTrim'): trim5 = trimpos[0][1]   # 5' trim
                if (seqlen-trimpos[-1][1]+1) <= self.getInt('VecTrim'): trim3 = trimpos[-1][0]   # 3' trim
                if trim5 > trim3:
                    trimdb.addEntry({'SeqName':seqname, 'SeqLen':seqlen, 'Start':1, 'End':seqlen, 'Edit':'fulltrim'})
                    purgelist.append(seqname)
                    continue
                if trim5 > 1:
                    trimdb.addEntry({'SeqName':seqname, 'SeqLen':seqlen, 'Start':1, 'End':trim5, 'Edit':'trim'})
                if trim3 < seqlen:
                    trimdb.addEntry({'SeqName':seqname, 'SeqLen':seqlen, 'Start':trim3, 'End':seqlen, 'Edit':'trim'})
            ## ~ [2d] ~ Identify regions to mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# 3. Finally, any remaining hits meeting the `vecmask=INT` length cutoff are masked by replacing with Ns.
                if self.getInt('VecMask') <= 0: continue
                for (mask5,mask3) in vecpos:
                    masklen = mask3-mask5+1
                    self.bugPrint('{}: ({},{})={} >= {}? (Trim <{}|{}<)'.format(seqname,mask5,mask3,masklen,self.getInt('VecMask'),trim5,trim3))
                    if masklen < self.getInt('VecMask'): continue
                    if mask3 <= trim5 or mask5 >= trim3: continue
                    trimdb.addEntry({'SeqName':seqname, 'SeqLen':seqlen, 'Start':mask5, 'End':mask3, 'Edit':'mask'})
                self.bugPrint('{0}: {1}'.format(seqname,vecdb.index('Hit')[seqname]))
                self.deBug('{0}: {1}'.format(seqname,vecpos))

            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # If any sequences are purged/trimmed/masked, output is then saved to:
            #
            # * `*.vecscreen.fasta` = The purged, trimmed and masked sequences. If no sequences were identified for editing,
            # this will not be generated.
            # * `*.vecpurge.tdt` = The table of vecsreen edits.
            #     - `SeqName`, `SeqLen`, `Start`, `End`, `Edit`
            trimdb.saveToFile()
            if not trimdb.entryNum(): return None
            trimdb.indexReport('Edit')
            ## ~ [3a] ~ Save sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist.obj['Current'] = None
            fasout = '{}.vecscreen.fasta'.format(basefile)
            rje.backup(self,fasout)
            FASOUT = open(fasout,'w')
            sx = 0.0; stot = seqlist.seqNum(); trimx = 0; dumpx = 0; maskx = 0; seqx = 0
            while seqlist.nextSeq():
                self.progLog('\r#VEC','VecPurge trimming/masking of : {:.2f}%'.format(sx/stot)); sx += 100.0
                seqname = seqlist.shortName()
                sequence = seqlist.seqSequence()
                seqdesc = seqlist.seqDesc()
                seqlen = len(sequence)
                if seqname in purgelist: dumpx += 1; continue
                #i# Mask
                trimmed = masked = False
                for entry in trimdb.indexEntries('SeqName',seqname):
                    if entry['Edit'] == 'mask' and maskmode == 'full':
                        sequence = sequence[:entry['Start']-1] + 'N' * (entry['End']-entry['Start']+1) + sequence[entry['End']:]
                        if len(sequence) != seqlen: raise ValueError('Masking problem!')
                        masked = True
                    elif entry['Edit'] == 'mask':
                        maskseq = sequence[:entry['Start']-1]
                        maskbase = True
                        for i in range(entry['Start']-1,entry['End']):
                            if maskbase: maskseq += 'N'
                            else: maskseq += sequence[i]
                            maskbase = not maskbase
                        maskseq += sequence[entry['End']:]
                        sequence = maskseq
                        if len(sequence) != seqlen: raise ValueError('Masking problem!')
                        masked = True
                    elif entry['Edit'] == 'trim':
                        sequence = sequence[:entry['Start']-1] + '-' * (entry['End']-entry['Start']+1) + sequence[entry['End']:]
                        if len(sequence) != seqlen: raise ValueError('Trimming problem!')
                        trimmed = True
                #i# Trim
                sequence = sequence.replace('-','')
                #i# Desc
                if trimmed:
                    seqdesc += ' (Vecscreen:trimmed)'
                    trimx += 1
                if masked:
                    seqdesc += ' (Vecscreen:masked)'
                    maskx += 1
                if (trimmed or masked) and not self.getBool('KeepNames'): seqname = seqname+'X'
                #i# Save
                FASOUT.write('>%s %s\n%s\n' % (seqname,seqdesc,sequence)); seqx += 1
            FASOUT.close()
            self.printLog('\r#SAVE','{} sequences processed; {} output to {}: {} trimmed; {} masked; {} dumped.'.format(rje.iStr(stot),rje.iStr(stot-dumpx),fasout,rje.iStr(trimx),rje.iStr(maskx),rje.iStr(dumpx)))
            ## ~ [3c] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if seqlist.getBool('Summarise'):
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % fasout,'autofilter=F']
                rje_seqlist.SeqList(self.log,seqcmd)
            return trimdb
        except:
            self.errorLog('Diploidocus.vecPurge() error'); raise
        return None
#########################################################################################################################
    #!# regcnv has been replaced with DepthKopy
    #!# should drop regcnv=T from regcheck mode until tidier
    def regCheck(self): ### Performs read check and/or CNV analysis of supplied region
        '''
        Performs read check and/or CNV analysis of supplied region. Based on VecCheck and SCDepth methods.

        ### ~ Region checking [runmode=regcheck] ~ ###

        Region checking, whether for read spanning analysis (`runmode=regcheck`) or copy number analysis
        (`runmode=regcnv` or `runmode=regcheck regcnv=T`), analyses regions extracted from a delimited file given by:
        `regcheck=FILE`. This can be a GFF file, in which case `gfftype=LIST` will restrict analysis to specific feature
        types, or regions can be defined by `checkfields=LIST`, which defines the locus, start and end positions. These
        fields default to `SeqName`, `Start` and `End` fields. If these fields cannot be found, the first three fields
        of the `regcheck=FILE` file will be used.

        Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are then mapped onto the assembly
        (`seqin=FILE`) using minimap2 to generate a PAF file. This is then parsed and reads spanning each feature based
        on their positions and the target start and end positions in the PAF file. In addition to absolute spanning of
        regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
        sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
        using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
        due to sequence length constraints.

        If `regcnv=T` then the region copy number analysis (below) will also be performed.

        **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
        parallel calculations is set by `forks=INT`.

        ---

        ### ~ Region copy number analysis [runmode=regcnv] ~ ###

        Copy number analysis uses the same single copy depth profile analysis as the `runmode=gensize` genome size
        prediction. In short, the modal read depth of BUSCO single copy `Complete` genes is calculated using samtools
        `mpileup` (or samtools `depth` if `quickdepth=T`) and used to defined "single copy read depth". BUSCO single-copy
        genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with
        `Status` = `Complete` and the `Contig`, `Start` and `End` fields are used to define the regions that should be
        predominantly single copy.

        Single-copy read depth can also be set using `scdepth=NUM` to re-used or over-ride previously calculated values.

        Long read data, given with the `reads=FILELIST` and `readtype=LIST` options, are then mapped onto the assembly
        (`seqin=FILE`) using minimap2. This can be skipped by providing a BAM file with `bam=FILE`. For each regcheck
        feature, the same samtools depth calculation is performed as for the BUSCO data. The mean depth across the region
        is then divided by the single copy depth to estimate total copy number across the region. Note that unlike the
        single copy depth estimation itself, this might be biased by repeat sequences etc. that have a different depth
        profile to the main region. One way to control for this might be to restrict analysis to a subset of reads that
        meet a certain minimum length cutoff, e.g. 10kb.

        **Query-based CNV analysis.** If the `regcheck=FILE` file has additional `Qry`, `QryLen`, `QryStart` and `QryEnd`
        fields, the copy number analysi will have an additional query-focus. In this case, each region mapping onto a
        specific query is summed up, adjusted for the proportion of the query covered by that region. For example, 3.5X
        mean depth of a 100% length copy and 3.0X coverage of a 50% length copy would sum to (3.5x1.0 + 3x0.5 = 5 total
        copies). If these fields are not present, each region will be analysed independently.

        **NOTE:** Depth calculations are performed in parallel in the directory set with `tmpdir=PATH`. The number of
        parallel calculations is set by `forks=INT`.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            rdir = '%slibraries/r/' % slimsuitepath
            if self.getBool('QuickDepth'): depmethod = 'depth'
            ## ~ [1a] ~ Check input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add feature to recognise and change first field if not given #!#
            if not len(self.list['CheckFields']) == 3:
                raise ValueError('checkfields=LIST must have exactly 3 elements: SeqName, Start, End. %d found!' % len(self.list['CheckFields']))
            [locusfield,startfield,endfield] = self.list['CheckFields']
            cdb = None
            #i# Check for GFF file and parse differently
            if self.getStrLC('RegCheck').endswith('.gff') or self.getStrLC('RegCheck').endswith('.gff3'):
                self.printLog('#GFF','GFF file recognised for RegCheck')
                if self.list['GFFType']:
                    self.printLog('#GFF','Parsing "{0}" features only'.format(','.join(self.list['GFFType'])))
                else: self.printLog('#GFF','Parsing all features types (no gfftype=LIST set)')
                gffhead = string.split('seqid source type start end score strand phase attributes')
                cdb = db.addTable(self.getStr('RegCheck'),mainkeys='auto',datakeys='All',delimit='\t',headers=gffhead,ignore=['#'],lists=False,name='check',expect=True)
                cdb.addField('Query')
                if self.list['GFFType']:
                    #?# Change this to a list of types? #?#
                    cdb.dropEntriesDirect('type',self.list['GFFType'],inverse=True)
                for entry in cdb.entries():
                    if rje.matchExp('Parent=([^;$]+)[;$]',entry['attributes']): entry['Query'] = rje.matchExp('Parent=([^;$]+)[;$]',entry['attributes'])[0]
                    elif rje.matchExp('ID=([^;$]+)[;$]',entry['attributes']): entry['Query'] = rje.matchExp('ID=([^;$]+)[;$]',entry['attributes'])[0]
                    elif rje.matchExp('Name=([^;$]+)[;$]',entry['attributes']): entry['Query'] = rje.matchExp('Name=([^;$]+)[;$]',entry['attributes'])[0]
                    entry['Query'] = '%s.%s' % (entry['Query'],entry['type'])
                    #i# Might want to add an ability to parse these from the feature description?
                    entry['QryLen'] = 100
                    entry['QryStart'] = 1
                    entry['QryEnd'] = 100
                cdb.renameField('seqid',locusfield)
                cdb.renameField('start',startfield)
                cdb.renameField('end',endfield)
            #i# Otherwise, load table as normal
            else:
                cdb = db.addTable(self.getStr('RegCheck'),mainkeys='auto',name='check',expect=True)
            if not cdb: raise IOError('Cannot find checkpos file "%s"' % self.getStr('CheckPos'))
            if locusfield not in cdb.fields():
                self.warnLog('Field "%s" not found in checkpos file: will use "%s for sequence name' % (locusfield,cdb.fields()[0]))
                locusfield = cdb.fields()
            if startfield not in cdb.fields():
                newstart = cdb.fields()[2]
                if 'Start' in cdb.fields(): newstart = 'Start'
                self.warnLog('Field "%s" not found in checkpos file: will use "%s for start position' % (startfield,newstart))
                startfield = newstart
            if endfield not in cdb.fields():
                newfield = cdb.fields()[2]
                if 'End' in cdb.fields(): newfield = 'End'
                self.warnLog('Field "%s" not found in checkpos file: will use "%s for end position' % (endfield,newfield))
                endfield = newfield
            self.list['CheckFields'] = [locusfield,startfield,endfield]
            checkcmd = ['checkfields=%s,%s,%s' % (locusfield,startfield,endfield)]
            cdb.dataFormat({startfield:'int',endfield:'int'})
            ## ~ [1b] ~ Check for query data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qrycnv = self.getStrLC('RunMode') in ['regcnv'] or self.getBool('RegCNV')
            for qfield in ['Qry','QryLen','QryStart','QryEnd']:
                if qfield not in cdb.fields(): qrycnv = False
            ## ~ [1c] ~ Update checkfields ~
            compfields = self.list['CheckFields'][0:]
            if qrycnv:
                cdb.dataFormat({'QryLen':'int','QryStart':'int','QryEnd':'int'})
                compfields = ['Qry','QryStart','QryEnd'] + compfields
            cdb.compress(compfields,rules={self.getStr('SpanID'):'list'})

            ### ~ [2] Simple PAF-based method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# This is the same as used by veccheck
            if self.getStrLC('RunMode') in ['regcheck']:
                vfile = self.getStr('RegCheck')
                #pfile = self.baseFile() + 'checkpos.paf'
                pfile = self.getPAFFile()   #baseFile() + 'checkpos.paf'
                pafcmd = self.cmd_list + ['checkpos={}'.format(vfile),'pafin={}'.format(pfile)] + checkcmd
                paf = rje_paf.PAF(self.log, pafcmd)
                paf.list['CheckFields'] = self.list['CheckFields']
                paf.obj['DB'] = self.db()
                #!# Add provision of existing cdb
                cdb = paf.checkPos(save=False,cdb=cdb)
                if cdb not in self.db().list['Tables']:
                    self.db().list['Tables'].append(cdb)
                cdb.setStr({'Name':'checkpos'})
                if not self.getBool('RegCNV'):
                    cdb.saveToFile(backup=False)
                    return True
                else: self.warnLog('RegCNV mode has been improved with DepthKopy - it is recommended to run that instead.')

            ### ~ [3] Complex BAM-based method for CNV calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('RunMode') in ['regcnv']: cdb.setStr({'Name':'checkcnv'})
            cdb.addFields(['SeqBP','ReadBP','MeanX','ModeX','DensX','CN'])
            ## ~ [3a] Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bamfile = self.getBamFile()
            if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            self.printLog('#BAM',bamfile)
            ## ~ [3b] ~ Establish SC read depth using samtools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
            scdepth = self.getNum('SCDepth')
            if self.getNum('SCDepth'):
                self.printLog('#SCDEP','Using loaded single copy read depth = {0:.2f}X'.format(scdepth))
            else:
                self.printLog('#SCDEP','No single copy read depth (scdepth=NUM): will calculate from BAM and BUSCO')
                scdepth = self.genomeSize(scdepth=True)
                self.printLog('#SCDEP','Using BUSCO-derived single copy read depth = {0}X'.format(scdepth))
                if not scdepth: raise ValueError('Failed to establish SC read depth')
                self.setNum({'SCDepth':scdepth})
            busdep = 'busco{}'.format(depmethod)
            bdb = self.db(busdep)
            if not bdb:
                busco = '{0}.{1}.tdt'.format(self.baseFile(),busdep)
                bdb = self.db().addTable(busco,mainkeys=['#'],expect=False,name=busco)
                if bdb: bdb.dataFormat({'MeanX':'num','DensX':'num'})
            if bdb:
                (buscX,buscSD) = rje.meansd(bdb.dataList(bdb.entries(),'MeanX',sortunique=False,empties=False))
                #i# CI of RegCNV estimate given RegX:
                regX = scdepth
                #i# 1. Systematic sequencing bias
                regCIN1 = 1.96 * buscSD * math.sqrt(regX / (buscX ** 3))
                regCIN2 = 1.96 * buscSD * regX / (buscX ** 2)
                cdb.addFields(['CIsyst','CIrand'])
            ## ~ [3c] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)
            forker = self.obj['Forker']
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for bentry in cdb.entries():
                if bentry[endfield] < bentry[startfield]: (bentry[endfield],bentry[startfield]) = (bentry[startfield],bentry[endfield])
                if qrycnv:
                    tmpfile = '{}{}.{}-{}-{}-{}.{}.tmp'.format(tmpdir,basefile,bentry['Qry'],bentry[locusfield],bentry[startfield],bentry[endfield],depmethod)
                else:
                    tmpfile = '{}{}.{}-{}-{}.{}.tmp'.format(tmpdir,basefile,bentry[locusfield],bentry[startfield],bentry[endfield],depmethod)
                if rje.exists(tmpfile):
                    if not self.force() and len(open(tmpfile,'r').readline().split()) == 2: skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                if depmethod == 'mpileup':
                    forker.list['ToFork'].append("samtools view -b -h -F 0x100 %s %s:%s-%s | samtools mpileup -BQ0 - 2> /dev/null | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $4;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry[locusfield],bentry[startfield],bentry[endfield],bentry[startfield],bentry[endfield],tmpfile))
                else:
                    forker.list['ToFork'].append("samtools view -h -F 0x100 %s %s:%s-%s | samtools depth - | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $3;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry[locusfield],bentry[startfield],bentry[endfield],bentry[startfield],bentry[endfield],tmpfile))
            self.printLog('#REGION','{} check regions queued for forking ({} existing files deleted); {} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [3d] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of region depth analysis completed.')
                else:
                    try:
                        self.errorLog('Samtools forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Samtools forking did not complete')
            ## ~ [3e] Load in data and calculate CNV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for bentry in cdb.entries():
                if qrycnv:
                    tmpfile = '{}{}.{}-{}-{}-{}.{}.tmp'.format(tmpdir,basefile,bentry['Qry'],bentry[locusfield],bentry[startfield],bentry[endfield],depmethod)
                else:
                    tmpfile = '{}{}.{}-{}-{}.{}.tmp'.format(tmpdir,basefile,bentry[locusfield],bentry[startfield],bentry[endfield],depmethod)
                try:
                    if not rje.exists(tmpfile):
                        raise IOError('Cannot find {}'.format(tmpfile))
                    try:
                        bentry['ModeX'] = int(open(tmpfile,'r').readline().split()[1])
                    except IndexError:
                        self.warnLog('Possible lack of primary read mapping to %s %s-%s (%s -> "%s")' % (bentry[locusfield],bentry[startfield],bentry[endfield],depmethod,tmpfile))
                        bentry['SeqBP'] = 0
                        bentry['ReadBP'] = 0
                        bentry['MeanX'] = 0.0
                        bentry['DensX'] = 0.0
                        bentry['CN'] = 0.0
                        self.printLog('#REGCNV','%s %s..%s = %.1fX -> %.2fN (1N=%dX)' % (bentry[locusfield],bentry[startfield],bentry[endfield],bentry['MeanX'],bentry['CN'],self.getNum('SCDepth')))
                        continue
                    except:
                        self.errorLog('Something has gone wrong with "%s"' % tmpfile)
                        raise
                    seqbp = 0
                    readbp = 0
                    for dline in open(tmpfile,'r').readlines():
                        data = rje.chomp(dline).split()
                        n = int(data[0])
                        X = int(data[1])
                        seqbp += n
                        readbp += (n * X)
                    bentry['SeqBP'] = seqbp
                    bentry['ReadBP'] = readbp
                    bentry['MeanX'] = (1.0 * readbp) / (bentry[endfield] - bentry[startfield] + 1)
                    try:
                        rcmd = 'Rscript {0}depmode.R {1} pure'.format(rdir, tmpfile)
                        bentry['DensX'] = float(rje.chomp(os.popen(rcmd).readlines()[0]))
                    except:
                        self.warnLog('Problem calling Rscript depmode.R for "%s"' % tmpfile,suppress=True)
                    #!# Add option for choice? Dev option for now.
                    cnfield = 'MeanX'
                    if self.dev(): cnfield = 'DensX'
                    bentry['CN'] = bentry[cnfield]/self.getNum('SCDepth')
                    if bdb:
                        regX = bentry['DensX']
                        bentry['CIsyst'] = 1.96 * buscSD * math.sqrt(regX / (buscX ** 3))
                        bentry['CIrand'] = 1.96 * buscSD * regX / (buscX ** 2)

                        self.printLog('#REGCNV','%s %s..%s = %.1fX -> %.2fN +/- %.3fN (95%% CI) (1N=%.2fX)' % (bentry[locusfield],bentry[startfield],bentry[endfield],bentry[cnfield],bentry['CN'],bentry['CIrand'],self.getNum('SCDepth')))
                    else:
                        self.printLog('#REGCNV','%s %s..%s = %.1fX -> %.2fN (1N=%.2fX)' % (bentry[locusfield],bentry[startfield],bentry[endfield],bentry[cnfield],bentry['CN'],self.getNum('SCDepth')))
                except:
                    self.errorLog('Samtools depth result processing error',quitchoice=self.debugging())
                    continue
            cdb.saveToFile(backup=False,sfdict={'CN':4,'MeanX':4,'DensX':4})

            ### ~ [4] Add QryCNV calculations and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# QryCNV does not currently have densemode option
            if qrycnv:
                cdb.addField('QryFrac')
                for bentry in cdb.entries():
                    bentry['QryFrac'] = (bentry['QryEnd'] - bentry['QryStart'] + 1.0) / bentry['QryLen']
                    bentry['CN'] = bentry['CN'] * bentry['QryFrac']
                cdb.setStr({'Name':'qrycnv'})
                cdb.compress(['Qry'],rules={'QryLen':'max','CN':'sum','QryFrac':'sum','SeqBP':'sum','ReadBP':'sum'})
                cdb.keepFields(['Qry','QryLen','QryFrac','SeqBP','ReadBP','MeanX','CN'])
                for bentry in cdb.entries():
                    bentry['MeanX'] = float(bentry['ReadBP']) / bentry['SeqBP']
                    if bdb:
                        regX = bentry['MeanX']
                        bentry['CIsyst'] = 1.96 * buscSD * math.sqrt(regX / (buscX ** 3))
                        bentry['CIrand'] = 1.96 * buscSD * regX / (buscX ** 2)
                        self.printLog('#QRYCNV','%s = %.1f%% @ %.1fX (%.2fN) -> %.2fN +/- %.3fN (95%% CI) (1N=%.2fX)' % (bentry['Qry'],100.0*bentry['QryFrac'],bentry['MeanX'],bentry['MeanX']/self.getNum('SCDepth'),bentry['CN'],bentry['CIrand'],self.getNum('SCDepth')))
                    else:
                        self.printLog('#QRYCNV','%s = %.1f%% @ %.1fX (%.2fN) -> %.2fN (1N=%.2fX)' % (bentry['Qry'],100.0*bentry['QryFrac'],bentry['MeanX'],bentry['MeanX']/self.getNum('SCDepth'),bentry['CN'],self.getNum('SCDepth')))
                    bentry['QryFrac'] = rje.dp(bentry['QryFrac'],2)
                cdb.saveToFile(backup=False,sfdict={'CN':4,'MeanX':4})

        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def gapSpan(self): ### Performs assembly gap read check and spanning analysis.
        '''
        ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

        This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
        existing) and then runs the read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
        with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
        This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
        In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
        sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
        using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
        due to sequence length constraints.

        Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
        `seqname.start-end`.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            mainbase = self.baseFile()
            basefile = self.baseFile(strip_path=True)
            seqbase =  rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            ## ~ [1a] ~ Check or create gapstats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=self.log):
                self.cmd_list.append('gapstats')    # Try to run automatically if possible
                seqin = self.seqinObj()
                if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=self.log):
                    seqin.setBool({'Raw':False,'GapStats':True,'DNA':True})
                    seqin.str['SeqType'] = 'dna'
                    seqin.summarise()
            if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=None):
                raise ValueError('Problem finding/generating %s.gaps.tdt' % seqbase)
            gapfile = '%s.gaps.tdt' % seqbase
            checkfile = '%s.checkpos.tdt' % basefile
            ## ~ [1b] ~ Setup PAF CheckFlanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.list['CheckFlanks']:
                self.list['CheckFlanks'] = [0]
                self.cmd_list.append('checkflanks=0')
            minflank = self.list['CheckFlanks'][0]
            minspan = 'Span%d' % minflank
            maxflank = self.list['CheckFlanks'][-1]

            ### ~ [2] ~ Load Gaps and run checkpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = None
            if rje.checkForFiles(filelist=['.checkpos.tdt'],basename=basefile,log=self.log) and not self.force():
                cdb = db.addTable(checkfile,mainkeys=['seqname','start','end'],name='checkpos',ignore=[],expect=True)
                cdb.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
                if minspan not in cdb.fields() and 'Span0' in cdb.fields():
                    self.printLog('#SPAN','Cannot find "{0}" field in {1}: will use Span0'.format(minspan,checkfile))
                    minspan = 'Span0'
                elif minspan not in cdb.fields():
                    self.printLog('#SPAN','Cannot find "{0}" field in {1}: will regenerate'.format(minspan,checkfile))
                    cdb = None
            if not cdb:
                cdb = db.addTable(gapfile,mainkeys=['seqname','start','end'],name='gap_pos',ignore=[],expect=True)
                cdb.dataFormat({'seqlen':'int','start':'int','end':'int'})
                cdb.makeField(formula='#seqname#.#start#-#end#',fieldname='gap')
                cdb.saveToFile()
                self.setStr({'RegCheck':cdb.saveToFileName()})
                self.list['CheckFields'] = ['seqname','start','end']
                checkcmd = ['checkfields=seqname,start,end','spanid=gap']
                ## ~ [2a] ~ Run Checkpos and SpanID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                vfile = self.getStr('RegCheck')
                pfile = self.getPAFFile()   #baseFile() + 'checkpos.paf'
                pafcmd = self.cmd_list + ['checkpos={}'.format(vfile),'pafin={}'.format(pfile)] + checkcmd
                paf = rje_paf.PAF(self.log, pafcmd)
                cdb = paf.checkPos(save=False)
                cdb.dataFormat({'MaxFlank3':'int','seqlen':'int','end':'int'})
                for centry in cdb.entries(): centry['MaxFlank3'] = centry['seqlen'] - centry['end']
                cdb.rename('checkpos')
                cdb.saveToFile(backup=False)
                #!# Clean up gap_pos file #!#

            ### ~ [3] ~ Reassemble gaps (runmode=gapass) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
            reassemble = self.getStrLC('RunMode') in ['gapass','gapfill']
            assembler = 'flye' # self.getStr('Assembler')
            filldir = '%s_gapfill/' % basefile
            assfile = '{0}{1}.assembledgaps.fasta'.format(filldir,basefile)
            if self.getStrLC('RunMode') in ['gapfill'] and rje.exists(assfile) and not self.force():
                self.printLog('#GAPASS','{0} found (force=F): skipping assembly'.format(assfile))
                reassemble = False
            if reassemble:
                try:
                    vcheck = string.split(os.popen('{0} --version'.format(assembler)).read())[0]
                    self.printLog('#PROG','{0} version: {1}'.format(assembler,vcheck))
                except:
                    raise ValueError('Assembler check error - failed to run: {0} --version'.format(assembler))
            if reassemble:
                rtype = self.list['ReadType'][0]
                if rtype in ['pacbio','pac']: rtype = 'pb'
                if rtype in ['hifi','ccs']: rtype = 'hifi'
                if rtype not in ['ont','pb','hifi']:
                    self.warnLog('Read Type "%s" not recognised (pb/ont/hifi): check readtype=LIST. Will use "ont".' % rtype)
                    rtype = 'ont'
                if assembler == 'flye':
                    if rtype == 'ont': rtype = 'nano-raw'
                    elif rtype == 'hifi': rtype = 'pacbio-hifi'
                    else: rtype = 'pacbio-raw'
                bamfile = self.getBamFile()
                iddir = '%s_spanid/' % basefile
                assdir = '%s_gapassemble/' % basefile
                rje.mkDir(self,assdir)
                nonx = 0; assx = 0; failx = 0; skipx = 0
                for centry in cdb.entries():
                    spanner = centry['gap']
                    target = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                    if rje.exists(target) and not self.force():
                        self.printLog('#GAPASS','{0} found (force=F): skipped.'.format(target))
                        skipx = 0
                        continue
                    spout = '%s%s.%s.span.id' % (iddir,basefile,spanner)
                    ## ~ [3a] ~ Check for read ID files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if not rje.exists(spout):
                        if centry[minspan] > 0:
                            self.warnLog('No {0} file but {1}={2}. May need to delete {3} and regenerate.'.format(spout,minspan,centry[minspan],checkfile))
                        elif self.v() > 1: self.printLog('No %s file' % spout)
                        nonx += 1
                        continue
                    ## ~ [3b] ~ Extract the reads from the BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #!# Convert this to a forked process too?
                    alog = '{0}{1}.{2}.assemble.log'.format(assdir,basefile,spanner)
                    self.printLog('#BAM','Extracting gap-spanning reads from BAM file.')
                    # Exclude secondary and umapped reads
                    region = '{0}:{1}-{2}'.format(centry['seqname'],max(centry['start']-maxflank,1),centry['end']+maxflank)
                    idbam = '%s%s.%s.bam' % (assdir,basefile,spanner)
                    fasout = '%s%s.%s.reads.fasta' % (assdir,basefile,spanner)
                    if self.force() or not rje.exists(fasout):
                        gap2bam = 'samtools view -h -F 0x104 {0} {1} | grep -e @ -f {2} | samtools view -h -b -o {3} -'.format(bamfile,region,spout,idbam)
                        logline = self.loggedSysCall(gap2bam,alog,append=True,stderr=False,nologline='No stdout from samtools view',threaded=False)
                        bam2fas = 'samtools fasta {0} > {1}'.format(idbam,fasout)
                        logline = self.loggedSysCall(bam2fas,alog,append=True,stderr=False,nologline='No stdout from samtools fasta',threaded=False)
                    gensize = maxflank*2
                    if rje.exists(fasout) and open(fasout,'r').read():
                        fascmd = ['dna=T','raw=T','gapstat=F','summarise=F','autoload=T','seqmode=list','seqin={0}'.format(fasout)]
                        fasta = rje_seqlist.SeqList(self.log,self.cmd_list+fascmd)
                        if fasta.seqNum() < 1:
                            self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                            continue
                        fasdat = fasta.summarise(sumdb=False,save=False)
                        gensize = int(fasdat['MeanLength'])
                    else:
                        self.warnLog('{0} has no reads despite {1} spanning read ID output'.format(fasout,spout))
                        continue
                    ## ~ [3c] Assemble reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #i# Adding forking
                    #i# asslist = []    # List of assembly commands to fork out and outputs (acmd, output, assembly)
                    if assembler == 'flye':
                        acmd = 'flye --{0} {1} --out-dir {2}{3}.{4}_flye --genome-size {5} --threads {6}'.format(rtype,fasout,assdir,basefile,spanner,gensize,self.getInt('SubForks'))
                        assembly = '{0}{1}.{2}_flye/assembly.fasta'.format(assdir,basefile,spanner)
                    if self.v() > 0:
                        acmd = '{0} 2>&1 | tee -a {1}'.format(acmd,alog)
                    else:
                        acmd = '{0} 2>&1 >> {1}'.format(acmd,alog)
                    asslist.append((acmd,assembly,'{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)))
                    continue
                    #!# Old code:
                    assembly = 'Assembly'
                    self.printLog('#GAPASS','Assembling gap-spanning reads for {0}.'.format(spanner))
                    if assembler == 'flye':
                        acmd = 'flye --{0} {1} --out-dir {2}{3}.{4}_flye --genome-size {5} --threads {6}'.format(rtype,fasout,assdir,basefile,spanner,gensize,self.getInt('SubForks'))
                        logline = self.loggedSysCall(acmd,alog,append=True,nologline='No stdout from flye',threaded=False)
                        assembly = '{0}{1}.{2}_flye/assembly.fasta'.format(assdir,basefile,spanner)
                    if rje.exists(assembly):
                        self.printLog('#GAPASS','{0} generated -> {1}'.format(assembly,target))
                        rje.fileTransfer(assembly,target,deletefrom=False,append=False)
                        assx += 1
                    else:
                        self.printLog('#GAPASS','No {0} generated.'.format(assembly))
                        failx += 1
                #self.printLog('#GAPASS','{0} Gap assemblies generated; {1} failed; {2} gaps without read IDs'.format(rje.iStr(assx),rje.iStr(failx),rje.iStr(nonx)))
                self.printLog('#GAPASS','{0} gaps without read IDs'.format(rje.iStr(nonx)))
                if skipx:
                    self.printLog('#GAPASS','{0} existing assemblies skipped (force=F).'.format(rje.iStr(skipx)))

                ### ~ [4] ~ Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.printLog('#FORK','{0} gap assemblies to fork out'.format(rje.iLen(asslist)))
                ## ~ [4a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                forker = self.obj['Forker']
                forker.list['ToFork'] = []
                for afork in asslist:
                    forker.list['ToFork'].append(afork[0])
                self.debug('\n'.join(forker.list['ToFork']))
                ## ~ [4b] Fork out assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if forker.list['ToFork']:
                    if self.getNum('Forks') < 1:
                        #i# Warn lack of forking
                        self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                        for acmd in forker.list['ToFork']:
                            alog = afork[2][:-6] + '.log'
                            logline = self.loggedSysCall(acmd,alog,append=True,nologline='No stdout from flye',threaded=False)
                            #self.printLog('#SYS',forkcmd)
                            #os.system(forkcmd)
                    else:
                        self.printLog('#FORK','Forking assemblies using {0} x {1} threads'.format(self.getNum('Forks'),self.getNum('SubForks')))
                        if forker.run():
                            self.printLog('#FORK','Forking of assemblies completed.')
                        else:
                            try:
                                self.errorLog('Assembly forking did not complete',printerror=False,quitchoice=True)
                            except:
                                raise RuntimeError('Assembly forking did not complete')
                ## ~ [4c] Process assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                assx = 0; failx = 0
                for afork in asslist:
                    assembly = afork[1]
                    target = afork[2]
                    if rje.exists(assembly):
                        self.printLog('#GAPASS','{0} generated -> {1}'.format(assembly,target))
                        rje.fileTransfer(assembly,target,deletefrom=False,append=False)
                        assx += 1
                    else:
                        self.printLog('#FAIL','No {0} generated.'.format(assembly))
                        failx += 1
                self.printLog('#GAPASS','{0} Gap assemblies generated; {1} failed'.format(rje.iStr(assx),rje.iStr(failx)))

            ### ~ [5] ~ Fill gaps (runmode=gapfill) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# The idea here is to use GABLAM-style PAF mapping to map the assembled regions back onto the assembly
            if self.getStrLC('RunMode') in ['gapfill']:
                ## ~ [4a] ~ Compile the assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                assdir = '%s_gapassemble/' % basefile
                filldir = '%s_gapfill/' % basefile
                rje.mkDir(self,filldir)
                assfile = '{0}{1}.assembledgaps.fasta'.format(filldir,basefile)
                if self.force() or not rje.exists(assfile):
                    rje.backup(self,assfile)
                    ASSFILE = open(assfile,'w')
                    nonx = 0; assx = 0; nullx = 0; ctgx = 0
                    for centry in cdb.entries():
                        spanner = centry['gap']
                        assembly = '{0}{1}.{2}.assembly.fasta'.format(assdir,basefile,spanner)
                        if not rje.exists(assembly): nonx += 1; continue
                        aseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=list','autoload=F','seqin={0}'.format(assembly)])
                        aseq.loadSeq()
                        if aseq.list['Seq']:
                            sx = 1; assx += 1
                            for (name, sequence) in aseq.list['Seq']:
                                ASSFILE.write('>{0}.ctg{1} {2}\n{3}\n'.format(spanner,sx,name,sequence))
                                sx += 1; ctgx += 1
                        else: nullx += 1
                    ASSFILE.close()
                    self.printLog('#GAPCTG','Collated {0} assembled contigs for {1} gaps; {2} assemblies w/o contigs; {3} w/o assemblies'.format(rje.iStr(ctgx),rje.iStr(assx),rje.iStr(nullx),rje.iStr(nonx)))
                else:
                    self.printLog('#GAPCTG','Collated assembled contigs found (force=F): {0}'.format(assfile))
                ## ~ [4b] ~ Map gap assemblies to genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fillbase = '%s%s.gapctg' % (filldir,self.baseFile())
                pafin = '%s.paf' % (fillbase)
                paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'%s%s.gapctg' % (filldir,self.baseFile()),'seqin=%s' % assfile,'reference=%s' % self.getStr('SeqIn')])
                defaults = paf_defaults
                paf.dict['MapOpt'] = rje.combineDict(defaults,paf.dict['MapOpt'],overwrite=True)
                paf.setup()
                if not rje.exists(pafin) or self.force():
                    rje.backup(self,pafin)
                    paf.minimap2()
                    #!# Add localaln=T to PAF #!#
                #?# Should paf.db() basename be changed to pafin basefile for the initial output?
                ## ~ [4c] ~ GapFiller local hits table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gapfiller = '%s.gapfiller.tdt' % (fillbase)
                paf.basefile(fillbase)
                if rje.exists(gapfiller) and not self.force():
                    udb = paf.db().addTable(gapfiller,mainkeys=['Hit','SbjStart','SbjEnd'],name='gapfiller',ignore=[],expect=True)
                    udb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','QryLen':'int','SbjLen':'int'})
                else:
                    paf.run()
                    paf.db('hitunique').rename('gapfiller')     # ['Qry','QryStart','QryEnd']
                    paf.db('gapfiller').dropFields(['nn','cs'])
                    for entry in paf.db('gapfiller').entries():
                        if entry['SbjStart'] > entry['SbjEnd']:
                            (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                    paf.db('gapfiller').newKey(['Hit','SbjStart','SbjEnd'])
                    paf.db('gapfiller').saveToFile()
                    #  Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd  QryLen  SbjLen  Strand
                self.printLog('#GAPMAP','%s%s.gapctg.* output generated.' % (filldir,self.baseFile()))
                ## ~ [4d] ~ Identify gap fills using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Work through by qry position, looking for a single hit or pair of hits that spans the gap
                #i# -> For a single hit, replace the whole local hit region
                #i# -> For a pair of hits, replace the gap between the Sbj regions
                #i# => Build a list of (seqname, start, end, gapctg, start, end, strand)
                #i# .. Reverse sort of replace sequences in turn
                #i# => Output to *.gapfill.fasta
                #i# -> Report number of gaps filled and new stats
                paf.basefile(mainbase)
                udb = paf.db('gapfiller')
                preventry = None
                gapfills = []   # Tuples of (seqname, start, end, gapctg, start, end, strand)
                filled = []
                #i# Filter on local alignment length (and identity?)
                minloclen = self.getInt('MinLocLen')
                minlocid = self.getPerc('MinLocID')
                if minloclen > 1: udb.dropEntries('Length<%d' % minloclen)
                if minlocid > 0:
                    udb.makeField('100.0*Identity/Length','LocID')
                    udb.dropEntries('LocID<%s' % minlocid)
                    udb.dropField('LocID')
                #i# Identify gaps to fill
                ex = 0.0; etot = udb.entryNum()
                for entry in udb.entries(sorted=True):
                    self.progLog('\r#FILL','Identifying gap-filling re-assemblies: %.2f%%' % (ex/etot)); ex += 100.0
                    gapid = '.'.join(entry['Qry'].split('.')[:-1])  # PTEXCHR1.01.10697073-10697172.ctg1
                    if gapid in filled: continue
                    #i# Check for matching sequences
                    if not entry['Qry'].startswith('{0}.'.format(entry['Hit'])): continue
                    #i# Check preventry
                    if preventry and preventry['Qry'] != entry['Qry']: preventry = None
                    #i# Check for single hit spanning gap
                    #if entry['Strand'] == '-': (entry['SbjStart'],entry['SbjEnd']) = (entry['SbjEnd'],entry['SbjStart'])
                    [gapx,gapy] = entry['Qry'].split('.')[-2].split('-')
                    gapx = int(gapx)
                    gapy = int(gapy)
                    if entry['SbjStart'] < gapx and entry['SbjEnd'] > gapy:
                        filled.append(gapid)
                        gapfills.append((entry['Hit'],entry['SbjStart'],entry['SbjEnd'],entry['Qry'],entry['QryStart'],entry['QryEnd'],entry['Strand'],'span'))
                        preventry = entry
                        continue
                    if not preventry:
                        preventry = entry
                        continue
                    #i# Check for dual hit spanning gap
                    if preventry['SbjEnd'] <= gapx and entry['SbjStart'] >= gapy and preventry['Strand'] == entry['Strand']:
                        if (entry['Strand'] == '+' and entry['SbjStart'] > preventry['SbjEnd']):
                            filled.append(gapid)
                            gapfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],preventry['QryEnd'],entry['QryStart'],entry['Strand'],'flank'))
                        elif (entry['Strand'] == '-' and preventry['SbjStart'] > entry['SbjEnd']):
                            filled.append(gapid)
                            gapfills.append((entry['Hit'],preventry['SbjEnd'],entry['SbjStart'],entry['Qry'],entry['QryEnd'],preventry['QryStart'],entry['Strand'],'flank'))
                    preventry = entry
                self.progLog('\r#FILL','Identified %s gap-filling re-assemblies from %s.' % (rje.iLen(gapfills),gapfiller))
                ## ~ [4e] ~ GapFill using local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Add a check for overlapping replacements? #!#
                gapfills.sort()
                gapshift = {}
                gdb = db.addEmptyTable('gapfill',['seqname','start','end','gapctg','ctgstart','ctgend','strand','gaplen','type','newname','newstart','newend','newchunklen'],['seqname','start','end'],log=False)
                for filldat in gapfills:
                    fentry = {'seqname':filldat[0],'start':filldat[1],'end':filldat[2],'gapctg':filldat[3],'ctgstart':filldat[4],'ctgend':filldat[5],'strand':filldat[6],'type':filldat[7]}
                    if fentry['seqname'] not in gapshift: gapshift[fentry['seqname']] = 0
                    fentry['newname'] = fentry['seqname']
                    fentry['newstart'] = fentry['start'] + gapshift[fentry['seqname']]
                    gapshift[fentry['seqname']] += (fentry['ctgend'] - fentry['ctgstart']) - (fentry['end'] - fentry['start'])
                    fentry['newend'] = fentry['end'] + gapshift[fentry['seqname']]
                    fentry['newchunklen'] = fentry['newend'] - fentry['newstart'] + 1
                    [gapx,gapy] = fentry['gapctg'].split('.')[-2].split('-')
                    fentry['gaplen'] = int(gapy) - int(gapx) + 1
                    gdb.dict['Data'][(filldat[0],filldat[1],filldat[2])] = fentry
                gdb.index('seqname')
                seqin = self.seqinObj(summarise=False)
                ctgin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(assfile),'gapstats=F'])
                ctgdict = ctgin.seqNameDic()
                fillfas = '{0}.fillcheck.fasta'.format(basefile)
                FILLFAS = open(fillfas,'w')
                fx = 0  # Gap fills
                nx = 0  # Unchanged sequences
                ex = 0  # Edited sequences
                for seq in seqin.seqs():
                    self.progLog('\r#FILL','Gap-filling sequences -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                    (seqname, sequence) = seqin.getSeq(seq)
                    sname = seqin.shortName(seq)
                    if sname not in gdb.index('seqname'):
                        nx += 1
                        FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                        continue
                    seqlen = len(sequence)
                    changes = gdb.indexEntries('seqname',sname)
                    changes.reverse()
                    for centry in changes:
                        ctgseq = ctgin.seqSequence(ctgdict[centry['gapctg']])
                        fillseq = ctgseq[centry['ctgstart']-1:centry['ctgend']]
                        if centry['strand'] == '-': fillseq = rje_sequence.reverseComplement(fillseq)
                        sequence = sequence[:centry['start']] + fillseq + sequence[centry['end']:]
                        fx += 1
                    ex += 1
                    newname = '{0}-{1}fill'.format(sname,len(changes))
                    for centry in changes: centry['newname'] = newname
                    self.printLog('\r#FILL','Gap-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                    FILLFAS.write('>{0} {1} (Diploidocus gap-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
                FILLFAS.close()
                self.printLog('\r#FILL','Gap-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} gaps filled.'.format(nx,ex,fillfas,fx,cdb.entryNum()))
                gdb.saveToFile()
                if fx:
                    rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','seqin={0}'.format(fillfas),'summarise=T','dna=T'])

                ## ~ [4f] Check the newly-filled gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # regcheck=FILE   : File of SeqName, Start, End positions for read coverage checking [None]
                # checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
                # checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
                # spanid=X        : Generate sets of read IDs that span veccheck/regcheck regions, grouped by values of field X []
                # regcnv=T/F      : Whether to calculate mean depth and predicted CNV of regcheck regions based on SCdepth [True]
                # gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
                # subforks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
                #!# Replace this with a new *.fillcheck.tdt file that has full-length filled gaps and just the edges.
                cdb = db.copyTable(gdb,'fillcheck',replace=True,add=True)
                cdb.newKey(['newname','newstart','newend'])
                for centry in cdb.entries():
                    newentry = {'newstart':centry['newstart']-1,'newend':centry['newstart'],'type':'start'}
                    newentry = rje.combineDict(newentry,centry,overwrite=False)
                    cdb.addEntry(newentry)
                    newentry = {'newstart':centry['newend'],'newend':centry['newend']+1,'type':'end'}
                    newentry = rje.combineDict(newentry,centry,overwrite=False)
                    cdb.addEntry(newentry)
                cdb.saveToFile()
                checkfile = '{0}.fillcheck.tdt'.format(basefile)
                checkcmd = self.cmd_list+['seqin={0}'.format(fillfas),'regcheck={0}'.format(checkfile),'checkfields=newname,newstart,newend','runmode=regcheck','paf=None','basefile={0}.fillcheck'.format(basefile)]
                if not self.getNum('SCDepth'):
                    checkcmd.append('regcnv=F')
                    self.printLog('#CHECK','Single-copy read depth (scdepth=NUM) not given: setting regcnv=F')
                if self.getInt('SubForks') > self.getInt('Forks'): checkcmd.append('forks={0}'.format(self.getInt('SubForks')))
                spanchecker = Diploidocus(self.log,checkcmd)
                spanchecker.run()

                ## ~ [4g] Fix any gaps without start and end support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Add reading in of *.fillcheck.checkpos.tdt and reversing gap-filling of gaps without edge support
                checkfile = '{0}.fillcheck.checkpos.tdt'.format(basefile)
                fdb = db.addTable(checkfile,mainkeys=['newname','newstart','newend','type'],name='fillchecked',ignore=[],expect=True)
                fformats = {}
                for field in fdb.fields():
                    if field not in ['seqname','gapctg','strand','type','newname']: fformats[field] = 'int'
                fdb.dataFormat(fformats)
                fdb.dropEntriesDirect('type',['start','end'],inverse=True)
                fdb.dropEntriesDirect('Span0',[0],inverse=True)
                fdb.compress(['newname','newstart','newend'],default='max',rules={'type':'list'},joinchar='&')
                self.printLog('#REVGAP','{0} filled gaps identified for reversion based on inserted chunk ends without spanning reads'.format(fdb.entryNum()))
                if fdb.entryNum():     #?# Add toggle to switch this off?
                    #i# Load sequences for reversion
                    origin = seqin
                    origdict = origin.seqNameDic()
                    seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','dna=T','seqin={0}'.format(fillfas),'gapstats=F'])
                    #i# Revert
                    fillfas = '{0}.gapfill.fasta'.format(basefile)
                    FILLFAS = open(fillfas,'w')
                    rx = 0  # Reverted gap fills
                    nx = 0  # Unchanged sequences
                    ex = 0  # Edited sequences
                    for seq in seqin.seqs():
                        self.progLog('\r#FILL','Reversing unsupported gap-filling -> {2}: {0} unchanged; {1} edited.'.format(nx,ex,fillfas))
                        (seqname, sequence) = seqin.getSeq(seq)
                        sname = seqin.shortName(seq)
                        if sname not in fdb.index('newname'):
                            nx += 1
                            FILLFAS.write('>{0}\n{1}\n'.format(seqname,sequence))
                            continue
                        seqlen = len(sequence)
                        changes = fdb.indexEntries('newname',sname)
                        changes.reverse()
                        for centry in changes:
                            ctgseq = origin.seqSequence(origdict[centry['seqname']])
                            fillseq = ctgseq[centry['start']-1:centry['end']]
                            sequence = sequence[:centry['newstart']] + fillseq + sequence[centry['newend']:]
                            rx += 1
                        ex += 1
                        newname = '{0}-{1}rev'.format(sname,len(changes))
                        #for centry in changes: centry['newname'] = newname
                        if rx: self.printLog('\r#FILL','Corrected gap-filled {0} -> {1}: {2} -> {3}.  '.format(sname,newname,rje_seqlist.dnaLen(seqlen),rje_seqlist.dnaLen(len(sequence))))
                        FILLFAS.write('>{0} {1} (Diploidocus gap-filled)\n{2}\n'.format(newname,seqin.seqDesc(seq),sequence))
                    FILLFAS.close()
                    self.printLog('\r#FASOUT','Corrected Gap-filled sequences output to {2}: {0} unchanged; {1} edited; {3} of {4} filled gaps reinstated.'.format(nx,ex,fillfas,rx,fx))
                else:
                    checkfas = fillfas
                    fillfas = '{0}.gapfill.fasta'.format(basefile)
                    os.rename(checkfas,fillfas)
                    self.printLog('#FASOUT','{0} renamed {1}'.format(checkfas,fillfas))
                rje_seqlist.SeqList(self.log, self.cmd_list + ['seqmode=file', 'autoload=T', 'seqin={0}'.format(fillfas), 'summarise=T', 'dna=T'])

            return True
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    ### <6> ### Genome Size and Purge Haplotigs methods                                                                 #
#########################################################################################################################
# echo EOG MPileup Depth > busco.depth.txt
# for REGION in $( grep Complete $BUSCO | cut -f 1,3-5 | sed 's/\t/:/g'); do
#   EOG=$(echo $REGION | awk -F ':' '{print $1;}')
#   LOC=$(echo $REGION | awk -F ':' '{print $2;}')
#   START=$(echo $REGION | awk -F ':' '{print $3;}')
#   END=$(echo $REGION | awk -F ':' '{print $4;}')
#   echo $REGION
#   DEP=$(samtools view -h -F 0x100 $BAM ${LOC}:${START}-${END} | samtools depth -  | awk -v pos="$START" '$2 >= pos' | awk -v pos="$END" '$2 <= pos' | awk '{print $3;}' | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')
#   PILE=$(samtools view -b -h -F 0x100 $BAM ${LOC}:${START}-${END} | samtools mpileup -BQ0 - | awk -v pos="$START" '$2 >= pos' | awk -v pos="$END" '$2 <= pos' |  awk '{print $4;}' | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')
#   echo $EOG $PILE $DEP | tee -a busco.depth.txt
# done
# echo 'Mpileup:'
# awk '{print $2;}' busco.depth.txt | sort | uniq -c | sort -nr | head
# echo 'Depth:'
# awk '{print $3;}' busco.depth.txt | sort | uniq -c | sort -nr | head
#########################################################################################################################
    def genomeSize(self,scdepth=False):   ### Uses read depth from BUSCO single copy genes to predict genome size
        '''
        Uses read depth from BUSCO single copy genes to predict genome size.
        >> scdepth:bool [False] = Whether to return single copy read depth only (w/o Genome Size prediction)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'DepthSizer' not in self.obj or not self.obj['DepthSizer']:
                self.obj['DepthSizer'] = depthsizer.DepthSizer(self.log, ['basefile=depthsizer'] + self.cmd_list)
                self.obj['DepthSizer'].setup()
            if scdepth: return self.obj['DepthSizer'].getSCDepth()
            estgensize = self.obj['DepthSizer'].depthSizer()
            self.setInt({'EstGenomeSize': estgensize})
            return estgensize
        except:
            self.errorLog('Diploidocus.genomeSize() error')
            return False
#########################################################################################################################
    def genomeSizeFromModeFile(self,dephist,scdepth=False):   ### Uses read depth from BUSCO single copy genes to predict genome size
        '''
        Uses read depth from BUSCO single copy genes to predict genome size.
        >> dephist:str = *.dephist.tdt file previously generated from BUSCO analysis
        >> scdepth:bool [False] = Whether to return single copy read depth only (w/o Genome Size prediction)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            #i# NOTE: readbp must have been pre-calculated for this method. Use genomeSize() if not.
            readbp = self.getInt('ReadBP')
            if scdepth and not self.getInt('GenomeSize') and readbp: scdepth = False
            if self.getBool('MapAdjust'): readbp *= self.getNum('MapAdjust')
            ## ~ [1a] ~ Load table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# headers=['Method','X','n','Modes'],
            mdb = db.addTable(dephist,mainkeys=['Method','X'],expect=True,name='dephist')
            mdb.dataFormat({'X':'int','n':'int','Modes':'int'})
            if depmethod not in mdb.index('Method'):
                self.printLog('#METHOD','Previous mode calculation did not use {}'.format(depmethod))
                return False
            ### ~ [2] Calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Generate mode and mode of modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            modelist = []
            for entry in mdb.indexEntries('Method',depmethod):
                n = entry['Modes']
                X = entry['X']
                modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'ModeOfModes':modelist[0][1]})
            self.printLog('#SCDEP','BUSCO SC mode of samtools {} modal read depth: {}X'.format(depmethod,self.getInt('ModeOfModes')))
            for entry in mdb.indexEntries('Method',depmethod):
                n = entry['n']
                X = entry['X']
                modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'BUSCOMode':modelist[0][1]})
            self.printLog('#SCDEP','BUSCO SC combined samtools {} modal read depth: {}X'.format(depmethod,self.getInt('BUSCOMode')))
            ## ~ [2b] Generate density modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DepDensity'):
                #!# Add check of Rscript #!#
                depfile = '{0}.dephist.tdt'.format(self.baseFile())
                rdir = '%slibraries/r/' % slimsuitepath
                try:
                    rcmd = 'Rscript {0}depmode.R {1} {2}'.format(rdir, depfile, depmethod)
                    self.printLog('#RCMD',rcmd)
                    depmode = float(rje.chomp(os.popen(rcmd).readlines()[0]))
                    self.setNum({'DensityMode':depmode})
                    self.printLog('#SCDEP','BUSCO SC combined samtools {0} modal density read depth: {1:.2f}X'.format(depmethod,depmode))
                    if scdepth: return self.getNum('DensityMode')
                except:
                    raise ValueError('Problem calling Rscript depmode.R')
            elif scdepth: return self.getInt('BUSCOMode')
            ### ~ [3] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.calculateGenomeSize(readbp)
            # estgensize = int(0.5+(float(readbp) / self.getInt('BUSCOMode')))
            # self.setInt({'EstGenomeSize':estgensize})
            # if not self.getInt('GenomeSize'): self.setInt({'GenomeSize':estgensize})
            # self.printLog('#GSIZE','Estimated genome size ({} at {}X): {}'.format(rje_seqlist.dnaLen(readbp,dp=0,sf=3),self.getInt('BUSCOMode'),rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
            # estgensize2 = int(0.5+(float(readbp) / self.getInt('ModeOfModes')))
            # self.printLog('#GSIZE','Mode of Modes Estimated genome size ({} at {}X): {}'.format(rje_seqlist.dnaLen(readbp,dp=0,sf=3),self.getInt('ModeOfModes'),rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
            # return self.getInt('BUSCOMode')
        except:
            self.errorLog('Diploidocus.genomeSizeFromModeFile() error')
            return False
#########################################################################################################################
    def legacyGenomeSize(self,scdepth=False,makebam=False):   ### Uses read depth from BUSCO single copy genes to predict genome size
        '''
        Uses read depth from BUSCO single copy genes to predict genome size.
        >> scdepth:bool [False] = Whether to return single copy read depth only (w/o Genome Size prediction)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            forker = self.obj['Forker']
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            readbp = self.getInt('ReadBP')
            if scdepth and not self.getInt('GenomeSize') and (readbp or self.list['Reads']): scdepth = False
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj()
            seqdict = seqlist.seqNameDic()
            ## ~ [1a] Check programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.popen('samtools --version').read():
                self.printLog('#SYS',' '.join(os.popen('samtools --version').read().split()))
            else:
                raise IOError('Cannot open samtools: check installation and/or module loading')
            ## ~ [1b] Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if makebam: bamfile = self.getBamFile()
            else:
                bamfile = self.getStr('BAM')
                if not self.getStrLC('BAM'): bamfile = self.baseFile(strip_path=True) + '.bam'
            if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            self.printLog('#BAM',bamfile)
            busco = self.getStr('BUSCO')
            if not self.getStrLC('BUSCO'): busco = 'full_table_{}.busco.tsv'.format(self.baseFile(strip_path=True))
            if not rje.exists(busco): raise IOError('Cannot find BUSCO full table "{}" (busco=FILE)'.format(busco))
            self.printLog('#BUSCO',busco)
            if not readbp and not self.list['Reads'] and not scdepth: raise IOError('Cannot find any read files (reads=FILELIST) and no readbp=X set.')
            ## ~ [1c] Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not readbp and self.list['Reads'] and not scdepth:
                for rfile in self.list['Reads']:
                    countfile = '{}.basecount.txt'.format(rfile)
                    if rje.exists(countfile) and not self.force():
                        filebp = int(rje.chomp(open(countfile,'r').readline()))
                        self.printLog('#READBP','%s: %s bases in %s' % (countfile,rje.iStr(filebp),rfile))
                        readbp += filebp
                        continue
                    self.progLog('\r#READBP','Counting bases in {}...'.format(rfile))
                    gzip = rfile.endswith('.gz')
                    fastq = rfile.endswith('q.gz') or rfile.endswith('q')
                    if gzip and fastq:
                        filebp = int(os.popen("zcat %s | grep '^+$' -B1 | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    elif fastq:
                        filebp = int(os.popen("grep '^+$' -B1 %s | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    elif gzip:
                        filebp = int(os.popen("zcat %s | grep -v '^>' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    else:
                        filebp = int(os.popen("grep -v '^>' %s | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                    self.printLog('#READBP','%s bases counted in %s' % (rje.iStr(filebp),rfile))
                    open(countfile,'w').write('%d\n' % filebp)
                    readbp += filebp
                self.setInt({'ReadBP':readbp})
            ## ~ [1d] MapAdjust Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapadjust = 1.0
            if not scdepth:
                self.printLog('#READBP', 'Total base count (unadjusted): %s' % (rje.iStr(readbp)))
                if self.getBool('MapAdjust'):
                    ratiofile = bamfile + '.mapratio.txt'
                    if self.force() or not rje.exists(ratiofile) or not open(ratiofile,'r').readline().split()[0]:
                        self.printLog('\r#ADJUST', 'Calculating Mapadjust ratio (samtools coverage)...')
                        covbases = float(os.popen("samtools coverage {0} | grep -v coverage | awk  '{{sum += ($7 * $5)}} END {{print sum}}'".format(bamfile)).read().split()[0])
                        self.printLog('\r#ADJUST', 'Calculating Mapadjust ratio (samtools fasta)...')
                        mapbases = float(os.popen("samtools fasta {0} | grep -v '^>' | wc | awk '{{ $4 = $3 - $2 }} 1' | awk '{{print $4}}'".format(bamfile)).read().split()[0])
                        open(ratiofile,'w').write('{0} {1}\n'.format(covbases,mapbases))
                    else:
                        self.printLog('\r#ADJUST','Use coverage:mapped ratio from {0} (force=F)'.format(ratiofile))
                        ratio = open(ratiofile, 'r').read().split()
                        covbases = float(ratio[0])
                        mapbases = float(ratio[1])
                    self.printLog('\r#ADJUST','Reference base read coverage: {0}'.format(rje.iStr(int(covbases))))
                    self.printLog('#ADJUST','Mapped read bases: {0}'.format(rje.iStr(int(mapbases))))
                    if not covbases or not mapbases: raise ValueError('Cannot have zero values for mapadjust. Something has gone wrong!')
                    mapadjust = covbases / mapbases
                    self.printLog('#ADJUST', 'Mapadjust ratio: {0:.3f}'.format(mapadjust))
                    readbp = readbp * mapadjust
                    self.printLog('#READBP', 'Total base count (adjusted): %s' % (rje.iStr(readbp)))
            self.setNum({'MapAdjust':mapadjust})
            ## ~ [1e] Load BUSCO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ext = rje.delimitExt(db.getStr('Delimit'))
            scfiles = ['.busco{}.{}'.format(depmethod,ext), '.buscodep.{}'.format(ext), '.dephist.{}'.format(ext)]
            if rje.checkForFiles(filelist=scfiles,basename=db.baseFile(),log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.'):
                runok = self.genomeSizeFromModeFile('{}.dephist.{}'.format(db.baseFile(),ext),scdepth)
                if runok: return runok
                else: self.printLog('#SCDEP','Unable to establish modal SC depth from {}'.format(scfiles[-1]))
            fullhead = ['BuscoID','Status','Contig','Start','End','Score','Length']
            bdb = db.addTable(busco,mainkeys='auto',headers=fullhead,expect=True,name='busco{}'.format(depmethod))
            bdb.dropEntriesDirect('Status',['Complete'],inverse=True)
            dephead = ['Method','BuscoID','X','n']
            depdb = db.addEmptyTable('buscodep',dephead,['Method','BuscoID','X'],log=self.debugging())
            ## ~ [1f] Temp directory for forked depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)

            ### ~ [2] Cycle through and fork out the depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []
            for seqname in bdb.index('Contig'):
                if seqname not in seqdict:
                    bx = rje.iLen(bdb.index('Contig')[seqname])
                    self.warnLog('"%s" (%s BUSCO-complete sequences) not found in "%s". May have been removed in previous purge cycle.' % (seqname,bx,seqin))
                    missing.append(seqname)
            if missing:
                bdb.dropEntriesDirect('Contig',missing,log=False)
            ## ~ [2a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for bentry in bdb.entries():
                # if bentry['Contig'] not in seqdict:
                #     self.warnLog('BUSCO-complete sequence "%s" not found in "%s". May have been removed in previous purge cycle.' % (bentry['Contig'],seqin))
                #     bdb.dropEntry(bentry)
                #     continue
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,bentry['BuscoID'],depmethod)
                if rje.exists(tmpfile):
                    #?# Add checking for completeness #?#
                    if not self.force() and len(open(tmpfile,'r').readline().split()) == 2: skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                if depmethod == 'mpileup':
                    forker.list['ToFork'].append("samtools view -b -h -F 0x100 %s %s:%s-%s | samtools mpileup -BQ0 - 2> /dev/null | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $4;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry['Contig'],bentry['Start'],bentry['End'],bentry['Start'],bentry['End'],tmpfile))
                else:
                    forker.list['ToFork'].append("samtools view -h -F 0x100 %s %s:%s-%s | samtools depth - | awk -v pos=\"%s\" '$2 >= pos' | awk -v pos=\"%s\" '$2 <= pos' | awk '{print $3;}' | sort | uniq -c | sort -nr > %s" % (bamfile,bentry['Contig'],bentry['Start'],bentry['End'],bentry['Start'],bentry['End'],tmpfile))
            self.printLog('#BUSCO','{} Complete BUSCO gene regions queued for forking ({} existing files deleted); {} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [2b] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of BUSCO depth analysis completed.')
                else:
                    try:
                        self.errorLog('Samtools forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Samtools forking did not complete')

            ### ~ [3] Read in and process depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            depmodes = {0:0}   # Dictionary of depth modes
            depcounts = {}  # Dictionary of depth counts
            ## ~ [3a] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb.addField('Mode')
            bdb.addField('MeanX')
            bdb.addField('DensX')
            for bentry in bdb.entries():
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,bentry['BuscoID'],depmethod)
                try:
                    if not rje.exists(tmpfile):
                        raise IOError('Cannot find {}'.format(tmpfile))
                    #?# Add full processing and calculation of mean coverage?
                    try:
                        bentry['Mode'] = int(open(tmpfile,'r').readline().split()[1])
                    except IndexError:
                        self.warnLog('Possible lack of primary read mapping to %s (%s -> "%s")' % (bentry['BuscoID'],depmethod,tmpfile))
                        bentry['MeanX'] = 0.0
                        bentry['DensX'] = 0.0
                        bentry['Mode'] = 0
                        depdb.addEntry({'Method':depmethod,'BuscoID':bentry['BuscoID'],'n':0,'X':0})
                        continue
                    except:
                        self.errorLog('Something has gone wrong with "%s"' % tmpfile)
                        raise
                    if bentry['Mode'] not in depmodes: depmodes[bentry['Mode']] = 0
                    depmodes[bentry['Mode']] += 1
                    seqbp = 0
                    bambp = 0
                    for dline in open(tmpfile,'r').readlines():
                        data = rje.chomp(dline).split()
                        n = int(data[0])
                        X = int(data[1])
                        seqbp += n
                        bambp += (n * X)
                        depdb.addEntry({'Method':depmethod,'BuscoID':bentry['BuscoID'],'n':n,'X':X})
                        if X not in depcounts: depcounts[X] = 0
                        depcounts[X] += n
                    bentry['MeanX'] = (1.0 * bambp) / seqbp
                    try:
                        rdir = '%slibraries/r/' % slimsuitepath
                        rcmd = 'Rscript {0}depmode.R {1} pure'.format(rdir, tmpfile)
                        bentry['DensX'] = float(rje.chomp(os.popen(rcmd).readlines()[0]))
                    except:
                        self.warnLog('Problem calling Rscript depmode.R for "%s"' % tmpfile,suppress=True)
                except:
                    self.errorLog('Samtools depth result processing error',quitchoice=self.debugging())
                    continue
            #!# Use the depdb table to make a box plot of the normalised depth counts?
            depdb.saveToFile()
            ## ~ [3b] Table of total and mode counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            depdb.compress(['Method','X'],rules={'n':'sum'})
            depdb.dropField('BuscoID')
            depdb.addField('Modes',evalue=0)
            for X, n in depmodes.items():
                dentry = depdb.data((depmethod,X))
                if dentry:
                    dentry['Modes'] = n
                elif n:
                    self.warnLog('Cannot find {0} {1}X in dephist table.'.format(depmethod,X))
            depdb.rename('dephist')
            depdb.saveToFile()
            ## ~ [3c] Generate mode and mode of modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            modelist = []
            for X, n in depmodes.items():
                if X:
                    modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'ModeOfModes':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC mode of samtools {} modal read depth: {}X'.format(depmethod,self.getInt('ModeOfModes')))
            modelist = []
            for X, n in depcounts.items(): modelist.append((n,X))
            modelist.sort(reverse=True)
            self.setInt({'BUSCOMode':modelist[0][1]})
            self.printLog('#MODE','BUSCO SC combined samtools {} modal read depth: {}X'.format(depmethod,self.getInt('BUSCOMode')))
            ## ~ [3d] Generate density modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DepDensity'):
                #!# Add check of Rscript #!#
                depfile = '{0}.dephist.tdt'.format(self.baseFile())
                rdir = '%slibraries/r/' % slimsuitepath
                try:
                    rcmd = 'Rscript {0}depmode.R {1} {2}'.format(rdir, depfile, depmethod)
                    self.printLog('#RCMD',rcmd)
                    depmode = float(rje.chomp(os.popen(rcmd).readlines()[0]))
                    self.setNum({'DensityMode':depmode})
                    self.printLog('#MODE','BUSCO SC combined samtools {} modal density read depth: {}X'.format(depmethod,depmode))
                except:
                    raise ValueError('Problem calling Rscript depmode.R')
            ## ~ [3e] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # !# Add option for choice? Dev option for now.
            cnfield = 'MeanX'
            if self.dev(): cnfield = 'DensX'

            bdb.addField('CN')
            bdb.addField('CN-MoM')
            for bentry in bdb.entries():
                bentry['CN'] = bentry[cnfield]/self.getInt('BUSCOMode')
                bentry['CN-MoM'] = bentry[cnfield]/self.getInt('ModeOfModes')
            (mean,se) = rje.meanse(bdb.dataList(bdb.entries(),'CN-MoM',sortunique=False,empties=False))
            median = rje.median(bdb.dataList(bdb.entries(),'CN-MoM',sortunique=False,empties=False))
            self.printLog('#CNV','BUSCO Complete Mode of Modes CNV depth check: %.2fN +/- %.3f (95%% CI); Median = %.2fN' % (mean,1.96 * se,median))
            (mean,se) = rje.meanse(bdb.dataList(bdb.entries(),'CN',sortunique=False,empties=False))
            median = rje.median(bdb.dataList(bdb.entries(),'CN',sortunique=False,empties=False))
            self.printLog('#CNV','BUSCO Complete CNV depth check: %.2fN +/- %.3f (95%% CI); Median = %.2fN' % (mean,1.96 * se,median))
            if self.getBool('DepDensity'):
                bdb.addField('CN-Density')
                for bentry in bdb.entries():
                    bentry['CN-Density'] = bentry[cnfield]/self.getNum('DensityMode')
                (mean,se) = rje.meanse(bdb.dataList(bdb.entries(),'CN-Density',sortunique=False,empties=False))
                median = rje.median(bdb.dataList(bdb.entries(),'CN-Density',sortunique=False,empties=False))
                self.printLog('#CNV','BUSCO Complete Density Mode CNV depth check: %.2fN +/- %.3f (95%% CI); Median = %.2fN' % (mean,1.96 * se,median))
            bdb.saveToFile()
            ## ~ [3f] Generate density modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DepDensity') and scdepth: return self.getNum('DensityMode')
            elif scdepth: return self.getInt('BUSCOMode')

            ### ~ [4] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #estgensize = int(0.5+(float(readbp) / self.getInt('BUSCOMode')))
            #self.setInt({'EstGenomeSize':estgensize})
            #if not self.getInt('GenomeSize'): self.setInt({'GenomeSize':estgensize})
            #self.printLog('#GSIZE','Estimated genome size ({} at {}X): {}'.format(rje_seqlist.dnaLen(readbp,dp=0,sf=3),self.getInt('BUSCOMode'),rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
            # return self.getInt('BUSCOMode')
            n = len(bdb.dataList(bdb.entries(),'CN-MoM',sortunique=False,empties=False))
            sd = se * math.sqrt(n)
            if (mean + 1.96 * se) < 1 and (median + sd) < 1:
                self.warnLog('BUSCO CNV mean < 1 and median + stddev < 1: genome size over-estimated?')
            if (mean - 1.96 * se) > 1 and (median - sd) > 1:
                self.warnLog('BUSCO CNV mean > 1 and median + stddev > 1: genome size under-estimated?')
            return self.calculateGenomeSize(readbp)
        except SystemExit: raise    # Child
        except:
            self.errorLog('Diploidocus.legacyGenomeSize() error')
            return False
#########################################################################################################################
    def calculateGenomeSize(self,readbp):  ### Calculates genome size from stored values and reports
        '''
        Calculates genome size from stored values and reports
        :param readbp: Total number of bp in read data
        :return: scdepth used for Genome Size estimate
        '''
        try:### ~ [0] Set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db('gensize')
            if not gdb: gdb = self.db().addEmptyTable('gensize',['SeqFile','DepMethod','ModeMethod','ReadBP','MapAdjust','SCDepth','EstGenomeSize'],['SeqFile','DepMethod','ModeMethod'])
            sizemethods = [('ModeOfModes','Mode of Modes'),('BUSCOMode','BUSCO Modes'),('DensityMode','BUSCO density mode')]
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            if self.getBool('MapAdjust'): depmethod += '-adjusted'
            ### ~ [1] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for (nkey, desc) in sizemethods:
                if self.getNum(nkey):
                    scdepth = self.getNum(nkey)
                    estgensize = int(0.5+(float(readbp) / scdepth))
                    scprint = scdepth
                    if nkey in ['DensityMode']: scprint = rje.dp(scdepth,3)
                    self.printLog('#GSIZE','{0} ({1}) Estimated genome size ({2} at {3}X): {4}'.format(desc,depmethod,rje_seqlist.dnaLen(readbp,dp=0,sf=3),scprint,rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
                    gdb.addEntry({'SeqFile':os.path.basename(self.getStr('SeqIn')),'DepMethod':depmethod,
                                  'ModeMethod':nkey,'ReadBP':self.getInt('ReadBP'),'MapAdjust':self.getNum('MapAdjust'),
                                  'SCDepth':scdepth,'EstGenomeSize':estgensize})
            self.setInt({'EstGenomeSize':estgensize})
            if not self.getInt('GenomeSize'): self.setInt({'GenomeSize':estgensize})
            return scdepth
        except:
            self.errorLog('Diploidocus.calculateGenomeSize() error')
            raise
#########################################################################################################################
    def assemblyMinimap(self):  ### Performs assembly versus self minimap2 search
        '''
        Performs assembly versus self minimap2 search
        :return: True/False
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

            #!# Modify the code below #!#

            #!# NOTE: This isn't currently used for anything. #!# Purge_haplotigs hits used #!#

            #mapopt=CDICT    : Dictionary of updated minimap2 options [N:250,p:0.0001,x:asm20]
            paf_defaults = {'N':'250','p':'0.0001','x':'asm20'}
            if self.getStr('Mapper') in ['minimap','minimap2']:
                pafin = '%s.paf' % self.baseFile()
                paf = rje_paf.PAF(self.log, self.cmd_list+['pafin=%s' % pafin,'basefile=%s' % self.baseFile(),'seqin=%s' % self.getStr('QueryDB'),'reference=%s' % self.getStr('SearchDB')])
                paf.setInt({'MinLocLen':self.getInt('LocalMin')})
                paf.setNum({'MinLocID':self.getNum('LocalIDMin')})
                defaults = paf_defaults
                paf.dict['MapOpt'] = rje.combineDict(defaults,paf.dict['MapOpt'],overwrite=True)
                paf.setup()
                if not rje.exists(pafin) or self.force():
                    rje.backup(self,pafin)
                    paf.minimap2()
                    #!# Add localaln=T to PAF #!#
                paf.run()
                paf.db('hitunique').rename('unique')
                #!# Change code to handle Qry and AlnNum #!#
                paf.db('unique').renameField('Qry','Query')
                paf.db('unique').renameField('AlnNum','AlnID')
                paf.db('unique').saveToFile()
                self.obj['BLAST'] = rje_blast.BLASTRun(log=self.log,cmd_list=['blastf=F']+self.cmd_list+['checkblast=F'])
                self.obj['BLAST'].obj['DB'] = self.obj['DB'] = paf.db()
                self.debug('PAF to GABLAM conversion complete')
                if self.getStr('NRStat') == 'OrderedAlnID':
                    self.printLog('#NRSTAT','NRStat updated to AlnID for minimap2 run.')
                    self.setStr({'NRStat':'AlnID'})
                return True
        except:
            self.errorLog('Diploidocus.assemblyMinimap() error')
            return False
#########################################################################################################################
    #!# Now ReadCore.longreadMinimap(paf=False)
    def LEGACYlongreadMinimap(self):  ### Performs long read versus assembly minimap2 and converts to BAM file
        '''
        Performs long read versus assembly minimap2 and converts to BAM file
        :return: bamfile/None
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['PAF'] = paf = rje_paf.PAF(self.log, self.cmd_list)
            bamfile = self.baseFile() + '.bam'
            mmv = os.popen('%s --version' % paf.getStr('Minimap2')).read()
            if not mmv: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % paf.getStr('Minimap2'))
            rje.backup(self,bamfile)

            ### ~ [2] ~ Generate individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Check these BAM files have headers! #!#
            bamlist = []; rx = 0
            if not self.list['ReadType']:
                self.warnLog('Read Type not given (pb/ont): check readtype=LIST. Will use "ont".')
            elif len(self.list['ReadType']) == 1 and len(self.list['Reads']) != 1:
                self.printLog('#READS','Using "%s" as read type for all long reads' % self.list['ReadType'][0])
            elif len(self.list['ReadType']) != len(self.list['Reads']):
                self.warnLog('reads=FILELIST vs readtype=LIST length mismatch: check readtype=LIST. Will cycle if needed.')
            for readfile in self.list['Reads']:
                if not rje.exists(readfile): raise IOError('Read file "{}" not found!'.format(readfile))
                if self.list['ReadType']:
                    try: rtype = self.list['ReadType'][rx]; rx +=1
                    except: rtype = self.list['ReadType'][0]; rx = 1
                    if rtype in ['pacbio','pac']: rtype = 'pb'
                    if rtype in ['hifi','ccs']: rtype = 'hifi'
                    if rtype not in ['ont','pb','hifi']:
                        self.warnLog('Read Type "%s" not recognised (pb/ont/hifi): check readtype=LIST. Will use "ont".' % rtype)
                        rtype = 'ont'
                else: rtype = 'ont'
                prefix = '{}.{}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(readfile,strip_path=True))
                maplog = '{}.log'.format(prefix)
                ## ~ [2a] Make SAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                imax = 4
                seqinsize = os.path.getsize(self.getStr('SeqIn')) / 1e9
                while seqinsize > imax: imax += 1
                maprun = '{} -t {} -I {}G --secondary=no -o {}.sam -L -ax map-{} {} {}'.format(paf.getStr('Minimap2'),self.threads(),imax,prefix,rtype,self.getStr('SeqIn'),readfile)
                if rtype in ['hifi']:
                    maprun = '{} -t {} -I {}G --secondary=no -o {}.sam -L -ax asm20 {} {}'.format(paf.getStr('Minimap2'),self.threads(),imax,prefix,self.getStr('SeqIn'),readfile)
                # self.printLog('#SYS',maprun)
                # if self.v() < 1:
                #     mapsys = subprocess.Popen(maprun, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #     stdout, stderr = mapsys.communicate()
                #     if stderr:
                #         logline = rje.chomp(stderr).split('\n')[-1]
                #         self.printLog('#SYSEND',logline)
                #     else: self.printLog('#SYSEND','No STDERR output produced during minimap2 mapping!')
                # else:
                #     #i# Generally useful to see minimap2 progress
                #     os.system(maprun)
                logline = self.loggedSysCall(maprun,maplog,append=False)
                #!# Add check that run has finished #!#
                if not rje.exists('{}.sam'.format(prefix)):
                    if self.i() > -1 and rje.yesNo('{}.sam missing! Pause and make manually?'.format(prefix)) and not rje.yesNo('{}.sam ready? Yes to continue; No to terminate.'.format(prefix)):
                        raise KeyboardInterrupt()
                if not rje.exists('{}.sam'.format(prefix)):
                    raise IOError('{}.sam missing!'.format(prefix))
                ## ~ [2b] Converting SAM to BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #sam2bam = 'samtools view -bo {}.tmp.bam -@ {} -S {}.sam'.format(prefix,self.threads()-1,prefix)
                self.printLog('#BAM','Converting SAM to BAM. Using a single thread due to past issues of missing data.')
                sam2bam = 'samtools view -bo {}.tmp.bam -S {}.sam'.format(prefix,prefix)
                logline = self.loggedSysCall(sam2bam,maplog,append=True,nologline='No stdout from sam2bam',threaded=False)
                #!# Add check that run has finished #!#
                ## ~ [2c] Sorting BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.printLog('#BAM','Sorting BAM file.')
                sortbam = '{}.bam'.format(prefix)
                mgb = self.getInt('MemPerThread')
                if self.getBool('UseQSub'):
                    ppn = self.threads()
                    vmem = self.getInt('QSubVMem')
                    mgb = min(mgb,int(vmem/float(ppn)))
                bamsort = 'samtools sort -@ {} -o {}.bam -m {}G {}.tmp.bam'.format(self.threads()-1,prefix,mgb,prefix)
                logline = self.loggedSysCall(bamsort,maplog,append=True)
                #!# Add check that run has finished #!#
                if not rje.exists(sortbam): raise IOError('Sorted BAM file "%s" not generated' % sortbam)
                os.unlink('{}.sam'.format(prefix))
                os.unlink('{}.tmp.bam'.format(prefix))
                bamlist.append(sortbam)

            ### ~ [3] ~ Merge individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(bamlist) > 1:
                # samtools merge - merges multiple sorted input files into a single output.
                bammerge = 'samtools merge -@ {} {} {}'.format(self.threads()-1,bamfile,' '.join(bamlist))
                logline = self.loggedSysCall(bammerge,append=True)
                if not rje.exists(bamfile): raise IOError('Merged BAM file "%s" not generated' % bamfile)
                for sortbam in bamlist: os.unlink(sortbam)
            else: os.rename(bamlist[0],bamfile)
            ## ~ [3a] Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            return bamfile
        except IOError: raise
        except KeyboardInterrupt: raise
        except:
            self.errorLog('Diploidocus.longreadMinimap() error')
            return None
#########################################################################################################################
    #i# Now part of ReadCore
    def LEGACYgetBamFile(self):  ### Checks/Creates indexed BAM file and returns filename as string
        '''
        Checks/Creates indexed BAM file and returns filename as string.
        :return: bamfile [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for program in ['minimap2','samtools']:
                if not os.popen('{} --version 2>&1'.format(program)).read():
                    self.warnLog('Cannot run "{} --version": check installation or pre-generation of files'.format(program))
            ## ~ [1a] IO names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = self.getStr('SeqIn')
            bamfile = self.getStr('BAM')

            ### ~ [2] BAM file check/generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('BAM') and not rje.exists(bamfile) and self.i() >= 0 and not rje.yesNo('Cannot find BAM file "{}" (bam=FILE): generate?'.format(bamfile)):
                raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            elif not self.getStrLC('BAM'):
                bamfile = self.baseFile(strip_path=True) + '.bam'
            # Generate if missing, or regenerate if inventing name and too old
            if not rje.exists(bamfile):
                self.printLog('#BAM','Cannot find BAM file "{}": generating'.format(bamfile))
                bamfile = self.longreadMinimap()
            elif not self.getStrLC('BAM') and self.needToRemake(bamfile,seqin):
                self.printLog('#BAM','SeqIn younger than BAM file "{}": regenerating'.format(bamfile))
                bamfile = self.longreadMinimap()
            # Check for existence now
            if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            self.setStr({'BAM':bamfile})
            baifile = '{}.bai'.format(bamfile)
            csifile = '{}.csi'.format(bamfile)
            #i# NOTE: If Diploidocus keeps remaking files, switch ignoredate=T
            ## ~ [2a] Index BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.checkForFiles(filelist=[bamfile,baifile,csifile],basename='',log=self.log,cutshort=False,ioerror=False)
            if self.needToRemake(baifile,bamfile) and self.needToRemake(csifile,bamfile):
                makebai = 'samtools index -b {} {}.bai'.format(bamfile,bamfile)
                logline = self.loggedSysCall(makebai,append=True,threaded=False,nologline='No stdout from samtools index')
            return bamfile
        except:
            self.errorLog('Diploidocus.getBamFile() error'); raise
#########################################################################################################################
    #i# Now part of ReadCore
    def LEGACYgetPAFFile(self):  ### Checks for PAF file and returns filename as string
        '''
        Checks for PAF file and returns filename as string.
        :return: paffile [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] IO names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = self.getStr('SeqIn')
            paffile = self.getStr('PAF')

            ### ~ [2] BAM file check/generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('PAF'):
                paffile = self.baseFile(strip_path=True) + '.paf'
                seqpaf = rje.baseFile(self.getStr('SeqIn'),strip_path=True) + '.paf'
                if not rje.exists(paffile) and rje.exists(seqpaf):
                    self.printLog('#PAF','Cannot find PAF file "{0}": using "{1}"'.format(paffile,seqpaf))
                    paffile = seqpaf
                elif not rje.exists(paffile):
                    self.printLog('#PAF','Cannot find PAF file "{0}": will generate "{1}"'.format(paffile,seqpaf))
                    paffile = seqpaf
            self.setStr({'PAF':paffile})
            return paffile
        except:
            self.errorLog('Diploidocus.getPAFFile() error'); raise
#########################################################################################################################
    def diploidocusHocusPocus(self):   ### Combines purge_haplotigs data with other stats to perform assembly filtering
        '''
        Diploidocus builds on the PurgeHaplotigs classifications to use the depth bins, KAT assembly kmer frequencies and
        BUSCO results to reclassify scaffolds and partition them into:

        * `*.diploidocus.fasta` = the scaffolds kept for the next round of PurgeHap
        * `*.core.fasta` = the same set of scaffolds, minus repeats
        * `*.repeats.fasta` = the repeat scaffolds excluded from core
        * `*.junk.fasta` = low coverage scaffolds, removed as junk
        * `*.quarantine.fasta` = putative haplotigs, removed from the assembly but retained for reference.

        **NOTE:** PurgeHaplotigs was not using zero coverage bases in its percentages. This is now fixed by Diploidocus.

        The first step in the classification decision process is to load the data and reduce it to the core data frame
        for analysis. `SeqName` will be used for rownames to help xref the full data.

        NOTE: The current method should be divided into submethods.

        NOTE: Purge_Haplotigs only calculates coverage based on non-zero coverage bases!
        Add includegaps=T/F and mingap=INT settings to include gaps in percentages (including lowcov) or not

        NOTE: Gaps are still not considered for the low coverage (Median depth) filter - set this to 0 if you are
        worried about retention of very gap-dense regions.

        NOTE: After adjustment, this additional filter may not be required!


        :return: dipdb
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [1a] IO names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            basefile = self.getStr('Basefile')
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.seqNameDic()
            bamfile = self.getBamFile()
            purgemode = self.getStrLC('PurgeMode')
            ## ~ [1b] Programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Make more nuanced. Add a partial=T/F mode, which can run with some of the programs? Or useX settings?
            if not os.popen('purge_haplotigs hist 2>&1').read():
                self.warnLog('Cannot run "purge_haplotigs hist": check installation or pre-generation of files')
            for program in ['kat','samtools','pileup.sh']:
                if not os.popen('{} --version 2>&1'.format(program)).read():
                    self.warnLog('Cannot run "{} --version": check installation or pre-generation of files'.format(program))
            ## ~ [1c] BAM file check/generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # if not self.getStrLC('BAM'):
            #     bamfile = self.baseFile(strip_path=True) + '.bam'
            #     if not rje.exists(bamfile):
            #         self.printLog('#BAM','Cannot find BAM file "{}": generating'.format(bamfile))
            #         bamfile = self.longreadMinimap()
            # if not rje.exists(bamfile): raise IOError('Cannot find BAM file "{}" (bam=FILE)'.format(bamfile))
            # self.setStr({'BAM':bamfile})
            # baifile = '{}.bai'.format(bamfile)
            # #i# NOTE: If Diploidocus keeps remaking files, switch ignoredate=T
            # rje.checkForFiles(filelist=[bamfile,baifile],basename='',log=self.log,cutshort=False,ioerror=False)
            # if self.needToRemake(baifile,bamfile):
            #     makebai = 'samtools index -b {} {}.bai'.format(bamfile,bamfile)
            #     logline = self.loggedSysCall(makebai,append=True)
            #     #os.system('samtools index -b {} {}.bai'.format(bamfile,bamfile))

            ### ~ [2] ~ Run PurgeHapolotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# This establishes SC read depth using samtools prior to running purge_haplotigs
            self.purgeHaplotigs()
            #scdepth = self.getNum('SCDepth')
            bamstrip = os.path.basename(bamfile)
            gencov = '{}.gencov'.format(bamstrip)
            covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            purge = '{}.purge.reassignments.tsv'.format(basefile)
            if not rje.checkForFiles(filelist=[gencov,covstats,purge],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: failed!'):
                phdir = 'purge_{}/'.format(basefile)
                raise IOError('Cannot find purge_haplotigs output. Check {}'.format(phdir))

            # ## ~ [2a] ~ Establish SC read depth using samtools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # # scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
            # scdepth = self.getNum('SCDepth')
            # if self.getNum('SCDepth'):
            #     self.printLog('#SCDEP','Using loaded single copy read depth = {}X'.format(scdepth))
            # else:
            #     scdepth = self.genomeSize(scdepth=True)
            #     self.printLog('#SCDEP','Using BUSCO-derived single copy read depth = {}X'.format(scdepth))
            #     if not scdepth: raise ValueError('Failed to establish SC read depth')
            # ## ~ [2b] ~ Setup purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # # phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
            # if self.getInt('PHLow') <= 0: self.setInt({'PHLow': int(float(scdepth)/4.0) })
            # phlow = self.getInt('PHLow')
            # # phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
            # if self.getInt('PHMid') <= 0:
            #     dupdepth = scdepth/2.0
            #     self.setInt({'PHMid': int(1.5 * dupdepth) })
            # phmid = self.getInt('PHMid')
            # # phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
            # if self.getInt('PHHigh') <= 0:
            #     self.setInt({'PHHigh': scdepth * 2 })
            # phhigh = self.getInt('PHHigh')
            # #?# Add checks and warnings of cutoff conflicts
            # ## ~ [2c] ~ Run purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # gencov = '{}.gencov'.format(bamfile)
            # covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            # purge = '{}.purge.reassignments.tsv'.format(basefile)
            # rje.checkForFiles(filelist=[gencov,covstats,purge],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            # #i# The -depth setting will be increased from 200 to 2xphhigh if >100
            # phcmd1 = 'purge_haplotigs hist -b {} -g {} -t {} -d {}'.format(bamfile,seqin,self.threads(),max(200,2*phhigh))
            # if self.needToRemake(gencov,bamfile):
            #     logline = self.loggedSysCall(phcmd1,append=True)
            # #!# Option to update the automatically set cutoffs
            # self.printLog('#PHDEP','Low=%dX; Mid=%dX; High=%dX. (SC=%dX)' % (phlow,phmid,phhigh,scdepth))
            # phcmd2 = 'purge_haplotigs cov -i {}.gencov -l {} -m {} -h {} -o {}.purge.coverage_stats.csv -j 80 -s 80'.format(bamfile,phlow,phmid,phhigh,basefile)
            # if self.needToRemake(covstats,gencov):
            #     logline = self.loggedSysCall(phcmd2,append=True)
            # else: self.printLog('#NOTE','Reusing existing %s on assumption that cutoffs have not changed' % covstats)
            # phcmd3 = 'purge_haplotigs purge -g {} -c {}.purge.coverage_stats.csv -t {} -o {}.purge -a 95'.format(seqin,basefile,self.threads(),basefile)
            # if self.needToRemake(purge,covstats):
            #     logline = self.loggedSysCall(phcmd3,append=True)
            if self.getStrLC('RunMode').startswith('purgehap'):
                self.printLog('#PURGE','purge_haplotigs run complete. Use runmode=diploidocus for additional filtering')
                return True

            ### ~ [3] ~ Run KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Make a simple python wrapper for kat to enable qsubbing? #!#
            if self.getBool('UseQSub'): self.printLog('#KAT','kat currently cannot be qsubbed due to problem with hanging runs.')
            ## ~ [3a] KmerReads data versus Assembly KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Add --5ptrim 16,0 setting, or additional kat options on commandline?
            katfile = '{}.kat-stats.tsv'.format(basefile)
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            # NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1 28      111164.29792    0.41661 122971501       122971475       54368   0.04421 122619711       99.71395        99.75805
            katcvg =  '{}.kat-counts.cvg'.format(basefile)
            rje.checkForFiles(filelist=[katfile,katcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            # Establish SC kmer count using BUSCO too
            # ${PREFIX}-counts.cvg
            # >NALACHR1.01 RevComp HiC_scaffold_38_pilon Len=123Mb; PBdep=39.8X; Nala German Shepherd Dog Chromosome 1
            # 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 11 528 585 553 470 111 ...
            katcall = 'kat sect -t {} -o {}.kat {} {}'.format(self.threads(),basefile,seqin,' '.join(self.list['KmerReads']))
            if self.getBool('10xTrim'):
                trim5 = ['16'] + ['0'] * (len(self.list['KmerReads']) - 1)
                trim5 = ','.join(trim5)
                katcall = 'kat sect -t {} --5ptrim {} -o {}.kat {} {}'.format(self.threads(),trim5,basefile,seqin,' '.join(self.list['KmerReads']))
            if not self.list['KmerReads']:
                self.printLog('#KAT','Cannot use KAT kmer analysis without KmerReads data')
            if self.list['KmerReads'] and (self.force() or not (rje.exists(katfile) and rje.exists(katcvg))):
                self.printLog('#SYS',katcall)
                #i# Catching completion in case KAT hangs after running
                KAT = os.popen(katcall)
                while not KAT.readline().startswith('Total runtime'): continue
                KAT.close()
            ## ~ [3b] Assembly versus self KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            selfkat = '{}.selfkat-stats.tsv'.format(basefile)
            selfcvg =  '{}.selfkat-counts.cvg'.format(basefile)
            rje.checkForFiles(filelist=[selfkat,selfcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            katcall = 'kat sect -t {} -o {}.selfkat {} {}'.format(self.threads(),basefile,seqin,seqin)
            if (self.force() or not (rje.exists(selfkat) and rje.exists(selfcvg))):
                self.printLog('#SYS',katcall)
                #i# Catching completion in case KAT hangs after running
                KAT = os.popen(katcall)
                while not KAT.readline().startswith('Total runtime'): continue
                KAT.close()

            ### ~ [4] ~ Run BBMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # module add java/8u201-jre bbmap/38.51
            depfile = '{}.depth.tdt'.format(basefile)
            depstat = '{}.depth.stats'.format(basefile)
            rje.checkForFiles(filelist=[depfile,depstat],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.')
            bbmap = 'samtools view -h {} | pileup.sh in=stdin out={}.depth.tdt 2>&1 | tee {}.depth.stats'.format(bamfile,basefile,basefile)
            if self.force() or not (rje.exists(depfile) and rje.exists(depstat)):
                logline = self.loggedSysCall(bbmap,append=True,threaded=False)
            #bblines = os.popen('samtools view -h {} | pileup.sh in=stdin out={}.depth.tdt 2>&1 | tee {}.depth.stats'.format(bamfile,basefile,basefile)).readlines()

            ### ~ [5] ~ Add VecScreen, Telomere and BUSCO Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [5a] Vecscreen coverage file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if self.getStrLC('ScreenDB'):
                self.printLog('#VEC','VecScreen screendb=FILE given -> running VecScreen')
                screencov = self.vecScreen()
            elif not rje.checkForFiles(filelist=[screencov],basename='',log=self.log):
                screencov = None
            vecdb = self.db('screencov')
            if screencov:
                if not vecdb: vecdb = db.addTable(screencov,mainkeys=['Query','Hit'],expect=True,name='screencov')
                vecdb.dropEntriesDirect('Query',['TOTAL'],inverse=True)
                vecdb.newKey(['Hit'])
                vecdb.setFields(['Hit','Coverage','CovPC'])
                vecdb.renameField('Coverage','ScreenCov')
                vecdb.renameField('CovPC','ScreenPerc')
            else:
                screencov = '{}.screencov.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            ## ~ [5b] Find Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            teldb = self.findTelomeres()
            ## ~ [5c] BUSCO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            busdb = None
            busco = self.getStr('BUSCO')
            if not self.getStrLC('BUSCO'): busco = 'full_table_{}.busco.tsv'.format(self.baseFile(strip_path=True))
            if rje.checkForFiles(filelist=[busco],basename='',log=self.log,missingtext='Not found: no BUSCO counts incorporated'):
                self.printLog('#BUSCO',busco)
                fullhead = ['BuscoID','Status','Contig','Start','End','Score','Length']
                busdb = db.addTable(busco,mainkeys='auto',headers=fullhead,expect=True,name='busco')
                for rating in ['Complete','Duplicated','Fragmented','Missing']: busdb.addField(rating,evalue=0)
                for bentry in busdb.entries(): bentry[bentry['Status']] = 1
                busdb.keepFields(['#','Contig','Complete','Duplicated','Fragmented'])
                busdb.compress(['Contig'],default='sum')
                busdb.dropField('#')
                missing = 0
                trimmed = 0
                for bentry in busdb.entries():
                    if bentry['Contig'] not in seqdict:
                        contigx = bentry['Contig'] + 'X'
                        if self.getBool('PreTrim') and contigx in seqdict:
                            bentry['Contig'] = contigx; trimmed += 1
                        else:
                            busdb.dropEntry(bentry); missing += 1
                if trimmed:
                    self.printLog('#BUSCO','{} BUSCO contigs mapped onto *X vecscreen trim names (pretrim=T)'.format(rje.iStr(trimmed)))
                self.printLog('#BUSCO','{} BUSCO contigs dropped: not found in {}'.format(rje.iStr(missing),seqin))

            ### ~ [6] ~ Compile Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            joinlist = []
            # >> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
            #     formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
            #     new table. If Field does not exist, it will be added. (Field may be a Formula.)
            # # ~ [6a] Summarise Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # depfile = '{}.depth.tdt'.format(basefile)
            # covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            # purge = '{}.purge.reassignments.tsv'.format(basefile)
            # katfile = '{}.kat-stats.tsv'.format(basefile)
            # selfkat = '{}.selfkat-stats.tsv'.format(basefile)
            rje.checkForFiles(filelist=[depfile,covstats,purge,katfile,selfkat,busco,screencov],basename='',log=self.log,missingtext='Not found: excluded from classification')

            ## ~ [6b] BBMap Depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.depth.tdt <==
            # #ID     Avg_fold        Length  Ref_GC  Covered_percent Covered_bases   Plus_reads      Minus_reads     Read_GC Median_fold     Std_Dev
            # Scaff10x_0  29.6599 14742825        0.0000  99.9948 14742058        28191   27814   0.4111  29      10.10
            depdb = db.addTable(depfile,mainkeys=['#ID'],expect=True,name='depth')
            depdb.renameField('#ID','SeqName')
            depdb.renameField('Length','SeqLen')
            depdb.setFields(['SeqName','SeqLen','Median_fold','Avg_fold','Covered_percent','Covered_bases','Plus_reads','Minus_reads','Read_GC'])
            joinlist.append((depdb,'SeqName'))

            ## ~ [6c] PurgeHaplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# NOTE: Had a problem with a sequence missing from *.coverage_stats.csv
            #  - will need some default values and warnings
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.purge.coverage_stats.csv <==
            # #contig,contig_reassign,bases_hap_dip,bases_low_high,bases_all,perc_low_coverage,perc_hap_coverage,perc_dip_coverage,perc_high_coverage
            # Scaff10x_0,,14719523,37459,14756982,0.004,9.823,89.923,0.250
            phcovdb = db.addTable(covstats,mainkeys=['#contig'],expect=True,name='phcov')
            phcovdb.renameField('#contig','SeqName')
            for cov in ['Low','Hap','Dip','High']:
                phcovdb.renameField('perc_{}_coverage'.format(cov.lower()), '{}Perc'.format(cov))
            phcovdb.setFields(['SeqName','LowPerc','HapPerc','DipPerc','HighPerc'])
            joinlist.append((phcovdb,'SeqName'))

            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.purge.reassignments.tsv <==
            # #reassigned_contig      top_hit_contig  second_hit_contig       best_match_coverage     max_match_coverage      reassignment
            # Scaff10x_1000   Scaff10x_562    Scaff10x_676    98.37   918.99  REPEAT
            purgedb = db.addTable(purge,mainkeys=['#reassigned_contig'],expect=True,name='purge')
            purgedb.renameField('#reassigned_contig','SeqName')
            purgedb.renameField('top_hit_contig','TopHit')
            purgedb.renameField('second_hit_contig','SecHit')
            purgedb.renameField('best_match_coverage','TopHitCov')
            purgedb.renameField('max_match_coverage','MaxHitCov')
            purgedb.renameField('reassignment','PurgeHap')
            purgedb.index('TopHit')
            purgedb.index('SecHit')
            purgedb.addFields(['TopNum','SecNum'],evalue=0)
            for entry in purgedb.entries():
                if entry['SeqName'] in purgedb.index('TopHit'): entry['TopNum'] = len(purgedb.index('TopHit')[entry['SeqName']])
                if entry['SeqName'] in purgedb.index('SecHit'): entry['SecNum'] = len(purgedb.index('SecHit')[entry['SeqName']])
            joinlist.append((purgedb,'SeqName'))

            ## ~ [6d] KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # ==> tiger.wtdbg2v1.racon2.10x.pilon2.scaffolds.kat-stats.tsv <==
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            # Scaff10x_0      1       1854.52173      0.41236 14757069        14757043        0       0.00000 14757043        100.00000       100.00000
            selfdb = db.addTable(selfkat,mainkeys=['seq_name'],expect=True,name='selfkat')
            selfdb.renameField('seq_name','SeqName')
            selfdb.renameField('median','SelfMedK')
            selfdb.renameField('mean','SelfAvgK')
            selfdb.setFields(['SeqName','SelfMedK','SelfAvgK'])
            joinlist.append((selfdb,'SeqName'))
            for entry in selfdb.entries(): entry['SeqName'] = entry['SeqName'].split()[0]
            selfdb.remakeKeys()

            #!# Add capacity for running without KAT -> MedK=-1
            katdb = db.addTable(katfile,mainkeys=['seq_name'],expect=False,name='kat')
            if katdb:
                katdb.renameField('seq_name','SeqName')
                katdb.renameField('median','MedK')
                katdb.renameField('mean','AvgK')
                katdb.renameField('gc%','SeqGC')
                katdb.renameField('%_non_zero_corrected','KPerc')
                katdb.setFields(['SeqName','MedK','AvgK','SeqGC','KPerc'])
                for entry in katdb.entries(): entry['SeqName'] = entry['SeqName'].split()[0]
                katdb.remakeKeys()
            else:
                katdb = db.addEmptyTable('kat',['SeqName','MedK','AvgK','SeqGC','KPerc'],['SeqName'])
                for entry in selfdb.entries(): katdb.addEntry({'SeqName':entry['SeqName'],'MedK':-1,'AvgK':-1,'SeqGC':-1,'KPerc':-1})
            joinlist.append((katdb,'SeqName'))

            ## ~ [6e] Other ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if vecdb:
                vecdb.renameField('Hit','SeqName')
                joinlist.append((vecdb,'SeqName'))
            if teldb:
                teldb.renameField('Name','SeqName')
                teldb.dropField('SeqLen')
                joinlist.append((teldb,'SeqName'))
            if busdb:
                busdb.renameField('Contig','SeqName')
                joinlist.append((busdb,'SeqName'))

            ## ~ [6x] Join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dipdb = db.joinTables(name='diploidocus',join=joinlist,newkey=['SeqName'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True,warnings=True)
            dipdb.remakeKeys()
            mingap = max(1,self.getInt('MinGap'))
            self.printLog('#GAPS','Identifying all runs of %d+ Ns as gaps' % mingap)
            if self.getBool('ZeroAdjust'):
                if self.getBool('IncludeGaps'):
                    self.printLog('#GAPS','Gap regions will be included in zero-coverage purge_haplotigs percentage adjustments (includegaps=T)')
                else:
                    self.printLog('#GAPS','Gap regions will be excluded from zero-coverage purge_haplotigs percentage adjustments (includegaps=F)')
            dipdb.addFields(['N_bases','Gap_bases'],evalue=-1)
            dipdb.addField('SeqDesc')
            ex = 0.0; etot = dipdb.entryNum(); gapwarn = 0; gaptot = 0
            for entry in dipdb.entries():
                self.progLog('\r#GAPS','Adding sequence gap counts and descriptions: %.1f%%' % (ex/etot)); ex += 100.0
                gapn = 0
                try:
                    seqname, seq = seqlist.getSeq(seqdict[entry['SeqName']])
                    seq = seq.upper()
                    entry['N_bases'] = seq.count('N')
                    seqdata = string.split(seqname,maxsplit=1)
                    if len(seqdata) > 1: entry['SeqDesc'] = seqdata[1]
                    gapn = len(''.join(re.findall('N{%d,}' % mingap,seq)))
                    gaptot += gapn
                    if gapn > (len(seq)/2.0): gapwarn += 1
                except:
                    self.warnLog('Could not extract description for %s!' % entry['SeqName'])
                    entry['SeqDesc'] = 'ERROR!'
                    open('seqdict.tmp','w').write('{}\n'.format(seqdict).replace(',','\n'))
                entry['Gap_bases'] = gapn
            self.printLog('\r#GAPS','Added %s bp sequence gaps (%d+ Ns) and descriptions.' % (rje.iStr(gaptot),mingap))
            if gapwarn: self.warnLog('{} sequences have >50% gaps. Check use of minmedian=X'.format(rje.iStr(gapwarn)))
            self.printLog('\r#FIELDS',', '.join(dipdb.fields()))
            #i# Tidy up join
            dipdb.fillBlanks(blank='False',fields=['Tel5','Tel3','Tel5Len','Tel3Len'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=-1,fields=['Trim5','Trim3'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=0,fields=['ScreenCov','Complete','Duplicated','Fragmented'],fillempty=True,prog=True,log=True)
            dipdb.fillBlanks(blank=0.0,fields=['TelPerc','ScreenPerc'],fillempty=True,prog=True,log=True)
            for entry in dipdb.entries():
                if entry['TopHitCov'] == '-': entry['TopHitCov'] = 0.0
                if entry['MaxHitCov'] == '-': entry['MaxHitCov'] = 0.0
            #!# Reorder dipdb fields
            fields = {'str':['SeqName', 'TopHit', 'SecHit', 'PurgeHap'],
                    'int': ['SeqLen', 'Median_fold', 'Covered_bases', 'Plus_reads', 'Minus_reads','TopNum','SecNum', 'SelfMedK', 'MedK', 'ScreenCov', 'Tel5Len', 'Tel3Len', 'Trim5', 'Trim3', 'Complete', 'Duplicated', 'Fragmented','Gap_bases','N_bases'],
                    'num': ['Avg_fold', 'Covered_percent', 'Read_GC', 'LowPerc', 'HapPerc', 'DipPerc', 'HighPerc', 'TopHitCov', 'MaxHitCov', 'SelfAvgK', 'AvgK', 'SeqGC', 'KPerc', 'ScreenPerc', 'TelPerc'],
                    'bool': ['Tel5', 'Tel3']}
            reformat = {}
            for ftype, tfields  in fields.items():
                for field in tfields: reformat[field] = ftype
            dipdb.dataFormat(reformat)

            ## ~ [6z] ~ Zero Adjustment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('ZeroAdjust'):
                ex = 0.0; etot = dipdb.entryNum(); adjx = 0
                for d in dipdb.entries():
                    self.progLog('\r#ZERO','Adjusting purge_hapotigs percentages for zero coverage: %.1f%%' % (ex/etot)); ex += 100.0
                    zerobp = d['SeqLen'] - d['Covered_bases']
                    if not self.getBool('IncludeGaps'): zerobp = max(0, zerobp - d['Gap_bases'])
                    if zerobp < 0: self.warnLog('Covered_bases < SeqLen for {}!'.format(d['SeqName'])); continue
                    if not zerobp: continue
                    if d['Covered_bases'] < 1:
                        d['LowPerc'] = 100.0
                        d['HapPerc'] = d['DipPerc'] = d['HighPerc'] = 0.0
                        adjx += 1
                        continue
                    norm = d['Covered_bases'] / 100.0
                    bpbins = [d['LowPerc'], d['HapPerc'], d['DipPerc'], d['HighPerc']]
                    try:
                        bpbins = [zerobp + (norm * d['LowPerc']), norm * d['HapPerc'], norm * d['DipPerc'], norm * d['HighPerc']]
                    except:
                        self.debug(d)
                        self.errorLog(self.wisdom())
                    newsum = sum(bpbins)
                    norm = 100.0 / newsum
                    d['LowPerc'] = norm * bpbins[0]
                    d['HapPerc'] = norm * bpbins[1]
                    d['DipPerc'] = norm * bpbins[2]
                    d['HighPerc'] = norm * bpbins[3]
                    adjx += 1
                self.printLog('\r#ZERO','Adjusted purge_hapotigs percentages for zero coverage in %s of %s sequences' % (rje.iStr(adjx),rje.iStr(etot)))

            ### ~ [7] ~ Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS SCAFFOLD RATING',line='=')
            dipdb.addField('Class', evalue='-')
            dipdb.addField('TopClass', evalue='-')
            dipdb.addField('SecClass', evalue='-')
            dipdb.addField('Rating', evalue='UNRATED')
            dipdb.addField('TopRating', evalue='-')
            dipdb.addField('SecRating', evalue='-')
            ## ~ [7a] ~ Cutoffs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            minmedian = self.getInt('MinMedian')
            self.printLog('#MINDEP','Low coverage min. median depth filter: {:d}X'.format(minmedian))
            mindipcov = 20  # Minimim DipPerc coverage to assign a "KEEP" rating
            self.printLog('#DIPCOV','Min. DipPerc coverage to assign as PRIMARY: {:d}%'.format(mindipcov))
            artefactcov = 80   # Min High/Low coverage to assign as artefacts
            self.printLog('#LOWCOV','Min. High/Low coverage to assign as artefacts: {:d}%'.format(artefactcov))
            hpurgecov = 80  # Min TopHitCov to assign HPURGE
            self.printLog('#HPURGE','Min. TopHitCov for HPURGE rating: {:d}%'.format(hpurgecov))
            rephitcov = 250  # Min MaxHitCov for REPEAT
            self.printLog('#REPEAT','Min. MaxHitCov for REPEAT rating: {:d}%'.format(rephitcov))
            pureperc = 80  # Depth class percentage to classify as pure
            self.printLog('#PURE','Depth class percentage to classify as pure: {:d}%'.format(pureperc))
            covwarning = 95
            self.printLog('#COV','Read coverage threshold for warnings: {:d}%'.format(covwarning))
            purge = 'normal'
            ## ~ [7b] ~ Pre-Rating Classifcation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                #i# Set up classification list to join
                eclass = ['WEAK','Dip']
                #i# Establish primary depth class
                for dep in ['Dip','Hap','High','Low']:
                    if entry['{}Perc'.format(dep)] >= pureperc: eclass[0] = 'PURE'
                    elif entry['{}Perc'.format(dep)] >= 50: eclass[0] = 'GOOD'
                if entry['Median_fold'] < minmedian: eclass[0] = 'LOWX'
                for dep in ['Hap','High','Low']:
                    if entry['{}Perc'.format(dep)] > entry['{}Perc'.format(eclass[1])]: eclass[1] = dep
                if eclass[1] == 'High': eclass[1] = 'EXS'
                else: eclass[1] = eclass[1].upper()
                #i# TopHit, SecHit and HitCov - Homology/Repeat status
                if not entry['TopHit'] or entry['TopHit'] == '-': eclass.append('UNIQ')
                elif entry['TopHitCov'] < 50: eclass.append('PART')
                elif not entry['SecHit'] or entry['SecHit'] == '-': eclass.append('HAPL')
                elif entry['MaxHitCov'] < rephitcov: eclass.append('HOMO')
                else: eclass.append('REPT')
                #i# TopHitNum and SecHitNum (v bad sequences will not be top hits)
                if entry['TopNum'] > 0: eclass.append('TOP')
                elif entry['SecNum'] > 0: eclass.append('SEC')
                else: eclass.append('NON')
                #i# KAT based self-assessment
                if entry['SelfMedK'] == 1: eclass.append('PRI')
                elif entry['SelfMedK'] == 2: eclass.append('ALT')
                else: eclass.append('REP')
                #i# BUSCO Genes (can help decide about risk)
                if entry['Complete'] > 0 and entry['Complete'] > entry['Duplicated']: eclass.append('COMP')
                elif entry['Duplicated'] > 0: eclass.append('DUPL')
                elif entry['Fragmented'] > 0: eclass.append('FRAG')
                else: eclass.append('NONE')
                #i# Extras
                if entry['TelPerc'] > 0: eclass.append('+TEL')
                if entry['ScreenCov'] > 0: eclass.append('+VEC')
                #i# Join
                entry['Class'] = '|'.join(eclass)
            dipdb.indexReport('Class')

            ## ~ [7b] ~ Default Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                tentry = dipdb.data(entry['TopHit'],expect=False)
                sentry = dipdb.data(entry['SecHit'],expect=False)
                if tentry:
                    try: entry['TopClass'] = tentry['Class']
                    except:
                        self.warnLog('Could not extract Class for TopHit "%s"' % entry['TopHit'])
                        tentry = None
                if sentry:
                    try: entry['SecClass'] = sentry['Class']
                    except:
                        self.warnLog('Could not extract Class for SecHit "%s"' % entry['SecHit'])
                        sentry = None
                # PURITY | DEPTH | HOM | TOP | MEDK | BUSCO
                cdata = entry['Class'].split('|')
                rep = cdata[4] == "REP"

                ### --- INITIAL LOW QUALITY FILTER --- ###
                # * `CONTAMINATION` = Scaffolds with 50%+ identified contamination
                if entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'
                #  - `Median_fold` < 3
                elif entry['Class'].startswith('LOWX'): entry['Rating'] = 'LOWCOV'
                elif entry['SeqLen'] < self.getInt('MinLen'): entry['Rating'] = 'LOWQUAL'
                #  - `[GOOD/PURE]|LOW` and `MedK`<1
                elif entry['Class'].startswith('PURE|LOW') and entry['MedK'] == 0: entry['Rating'] = 'LOWCOV'
                elif entry['Class'].startswith('GOOD|LOW') and entry['MedK'] == 0: entry['Rating'] = 'LOWCOV'
                # - `LOW` & `REP`
                elif cdata[1] == 'LOW' and rep: entry['Rating'] = 'LOWCOV'

                ### --- KEEP --- ###
                # * `QUALITY` = Highest quality scaffolds: pure duploid, complete BUSCOs, no Duplicated BUSCOs
                #  - `PURE|DIP|UNIQ` & `PRI|COMP` & `Duplicated` = 0
                elif entry['Class'].startswith('PURE|DIP|UNIQ') and "PRI|COMP" in entry['Class'] and entry['Duplicated'] == 0: entry['Rating'] = 'QUALITY'
                #* `FINAL` = As Quality but `PRI|FRAG` and `PRI|NONE` also allowed
                elif entry['Class'].startswith('PURE|DIP|UNIQ') and cdata[4] == 'PRI' and entry['Duplicated'] == 0: entry['Rating'] = 'FINAL'
                #* `CORE` = Predominantly Diploid with <50% covered by TopHit and SelfMedK=1
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "PRI": entry['Rating'] = 'CORE'
                #* `COREHAP` = Predominantly haploid but <50% covered by TopHit and 1+ Complete BUSCOs
                #  - `[PURE/GOOD]|HAP` & `UNIQ/PART` & `Complete` > 0
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and entry['Complete'] > 0: entry['Rating'] = 'COREHAP'
                #* `PRIMARY` = Putative primary scaffold but with possible alternative scaffolds still in assembly and/or low quality regions
                elif cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART","HAPL","HOMO") and not rep: entry['Rating'] = 'PRIMARY'

                #* `PRIRPT` = Putative primary scaffold but >50% repeated
                #  - DIP with UNIQ/PART and REP
                elif cdata[1] == "DIP" and cdata[2] in ("UNIQ","PART") and rep: entry['Rating'] = 'PRIRPT'
                #   - `[PURE/GOOD]|DIP|[HAPL/HOMO]` & `REP`
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "DIP" and rep: entry['Rating'] = 'PRIRPT'

                #* `COLLAPSED` = High coverage scaffolds representing putative collapsed repeats. Check for other contamination.
                #  - `EXS` & `UNIQ`
                elif cdata[1] == "EXS" and cdata[2] in ("UNIQ","PART"): entry['Rating'] = 'COLLAPSED'
                # - `[GOOD/PURE]|EXS` & `TOP`
                elif cdata[0] in ("PURE","GOOD") and cdata[1] == "EXS" and cdata[3] == "TOP": entry['Rating'] = 'COLLAPSED'
                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)<=`DipPerc`
                elif cdata[1]  == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) <= entry['DipPerc']: entry['Rating'] = 'COLLAPSED'

                # * `REPEAT` = Predominantly Diploid scaffolds that have major signs of redundancy, probably due to presence of alternative contigs
                #   - `DIP|REPT`
                elif "|DIP|REPT|" in entry['Class']: entry['Rating'] = 'REPEAT'
                #- `WEAK|DIP|[HAPL/HOMO]` & `REP`
                elif cdata[0] == 'WEAK' and cdata[1] == 'DIP' and cdata[2] in ("HAPL","HOMO") and rep: entry['Rating'] = 'REPEAT'

                #* `HAPLOID` = Predominantly haploid coverage but enough unique sequence to keep for now? (Option to quarantine?) Might be very heterozygous Alternatice haplotigs
                #  - `HAP` `UNIQ/PART` & `PRI`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "PRI": entry['Rating'] = 'HAPLOID'
                #* `HAPLOTIG` = Predominantly haploid coverage but enough unique sequence to keep for now and/or a TOP hit? Probable Alternative haplotig
                #  - `HAP` `UNIQ/PART` & `ALT`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and cdata[4] == "ALT": entry['Rating'] = 'HAPLOTIG'
                elif cdata[1] == "HAP" and cdata[3] == "TOP" and not rep: entry['Rating'] = 'HAPLOTIG'
                #* `HAPRPT` = Low quality scaffold that is probably most repeats, but not bad enough to dump outright
                #  - `HAP` & `UNIQ/PART` & `REP`
                elif cdata[1] == "HAP" and cdata[2] in ("UNIQ","PART") and rep: entry['Rating'] = 'HAPRPT'
                elif cdata[1] == "HAP" and cdata[3] == "TOP" and rep: entry['Rating'] = 'HAPRPT'
                ## Other HAPLOTIGs
                #  - `HAP/LOW` & `HAPL/HOMO/REPT` & `TopClass`-`DIP`
                elif cdata[1] in ("LOW","HAP") and cdata[2] not in ("UNIQ","PART") and '|DIP|' in entry['TopClass'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'
                elif cdata[0] in ('PURE','GOOD') and '|HAP|HOMO|' in entry['Class'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'
                # Assess remaining HAP scaffolds
                # - `HAP` and dipdb$TopHitCov < hpurgecov = HAPREPEAT
                elif cdata[1] == "HAP" and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPLOTIG'

                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)>`DipPerc`
                elif cdata[1] == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) > entry['DipPerc'] and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'HAPRPT'
                # - `[GOOD/PURE]|EXS[HAPL/HOMO/REPT]` & not `TOP`
                elif cdata[1] == "EXS" and cdata[0] != "WEAK" and cdata[2] not in ("UNIQ","PART") and cdata[3] != "TOP" and entry['TopHitCov'] < hpurgecov: entry['Rating'] = 'REPEAT'

                ### --- DUMP --- ###
                #* `RPURGE` = Messy scaffolds that are largely repeats and are sufficiently redundant/low quality to purge
                #  - `WEAK|EXS` & (`HapPerc`+`LowPerc`)>`DipPerc`
                elif cdata[1] == "EXS" and cdata[0] == "WEAK" and (entry['HapPerc'] + entry['LowPerc']) > entry['DipPerc'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'RPURGE'
                # - `[GOOD/PURE]|EXS[HAPL/HOMO/REPT]` & not `TOP`
                elif cdata[1] == "EXS" and cdata[0] != "WEAK" and cdata[2] not in ("UNIQ","PART") and cdata[3] != "TOP" and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'RPURGE'

                # - `LOW` & `HAPL/HOMO/REPT` & NOT `TOP`
                elif cdata[1] == 'LOW' and cdata[2] not in ['UNIQ','PART'] and cdata[3] != 'TOP': entry['Rating'] = 'LOWCOV'

                #* `HPURGE` = Clear candidate haplotig to purge
                #  - `HAP/LOW` & `HAPL/HOMO/REPT` & `TopClass`-`DIP`
                elif cdata[1] in ("LOW","HAP") and cdata[2] not in ("UNIQ","PART") and '|DIP|' in entry['TopClass'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                elif cdata[0] in ('PURE','GOOD') and '|HAP|HOMO|' in entry['Class'] and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                # - `HAP` and dipdb$TopHitCov >= hpurgecov = HPURGE
                elif cdata[1] == "HAP" and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'

                ### ---- LEFTOVERS ---- ###
                # * `LOWQUAL` = Unconvincing scaffolds that do not fall into a clear class but are not bad enough to dump outright
                #   - `LOW` & otherwise escaping `LOWCOV` filtering
                elif cdata[1] == "LOW": entry['Rating'] = 'LOWQUAL'

                ### Change HAPLOTIG to HAPRPT
                if entry['Rating'] == "HAPLOTIG" and rep: entry['Rating'] = 'HAPRPT'

                ### Special output
                #if self.getBool('Diploidify') and entry['Rating'] == 'HPURGE': entry['Rating'] = 'HAPLOTIG'

            ## ~ [7b] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if purgemode == 'simple':
                for entry in dipdb.entries():
                    # PURITY | DEPTH | HOM | TOP | MEDK | BUSCO
                    cdata = entry['Class'].split('|')
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    elif entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'; continue
                    elif entry['SeqLen'] < self.getInt('MinLen'): entry['Rating'] = 'LOWQUAL'; continue
                    #i# Low quality artefacts (e.g. less than 50% of the scaffold covered by short-read kmers)
                    elif entry['MedK'] == 0: entry['Rating'] = 'LOWQUAL'; continue
                    #i# KEEP set are those with DIP, EXS, UNIQ, PART, TOP or COMP = REPT/REP = REPEAT
                    if cdata[1] in ['DIP','EXS'] or cdata[2] in ['UNIQ','PART'] or cdata[3] == 'TOP' or entry['Complete'] > 0:
                        # PRIMARY
                        if cdata[1] in ['DIP'] or cdata[2] in ['UNIQ'] or cdata[5] == 'COMP':  #?# entry['Complete'] > 0:
                            if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'PRIRPT'
                            else: entry['Rating'] = 'PRIMARY'
                        # COLLAPSED
                        elif cdata[1] == 'EXS': entry['Rating'] = 'COLLAPSED'
                        else:
                            if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'REPEAT'
                            else: entry['Rating'] = 'HAPLOID'
                    #i# QUARANTINE any HAP
                    elif cdata[1] == 'HAP':
                        if cdata[2] == 'REPT' or cdata[4] == 'REP': entry['Rating'] = 'RPURGE'
                        else: entry['Rating'] = 'HPURGE'
                    #i# JUNK anything else LOW
                    elif cdata[1] == 'LOW': entry['Rating'] = 'LOWCOV'
                    else: raise ValueError('No rating: %s' % entry)
                    if entry['SeqName'] == 'cam006_MBG479__MBG479C13.006':
                        self.deBug('Rating: %s' % entry)

            ## ~ [7b] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'dev':
                for entry in dipdb.entries():
                    tentry = dipdb.data(entry['TopHit'],expect=False)
                    sentry = dipdb.data(entry['SecHit'],expect=False)
                    if tentry: entry['TopClass'] = tentry['Class']
                    if sentry: entry['SecClass'] = sentry['Class']
                    #i# INITIAL FILTERS #i#
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    elif entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'; continue
                    elif entry['SeqLen'] < self.getInt('MinLen'): entry['Rating'] = 'LOWQUAL'; continue
                    #i# Any scaffolds with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment.
                    # Scaffolds with 80%+ bases in the low/haploid coverage bins and 95%+ of their length mapped by
                    # PurgeHaplotigs onto another scaffold were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                    elif (entry['HapPerc']) >= 50 and entry['TopHitCov'] >= hpurgecov and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    # Any other scaffolds exceeding the mindipcov DipPerc coverage assigned to KEEP for further assignment
                    elif (entry['DipPerc']) >= mindipcov and entry['SelfMedK'] == 1: entry['Rating'] = 'PRIMARY'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'REPEAT'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov: entry['Rating'] = 'KEEP'  # 'PRIMARY' ?
                    # Any other scaffolds with 80%+ low coverage bases were filtered as Low Coverage.
                    elif (entry['LowPerc']) >= artefactcov: entry['Rating'] = 'LOWCOV'
                    # Any other scaffolds with 80%+ high coverage bases are flagged as COLLAPSED.
                    elif (entry['HighPerc']) >= artefactcov: entry['Rating'] = 'COLLAPSED'
                    # Scaffolds that are mostly Haploid Depth are marked as HAPLOTIGS if mostly present 2+ times, or HAPLOID if present once
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 1: entry['Rating'] = 'HAPLOID'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif entry['HapPerc'] < 50 and (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'
                    # Scaffolds with at least 80% hitting other scaffolds are HAPLOTIG or REPEAT
                    elif entry['TopHitCov'] >= hpurgecov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'HAPREPEAT'
                    elif entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HAPLOTIG'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'

            ## ~ [7x] ~ Basic Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'crude':
                for entry in dipdb.entries():
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    elif entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'; continue
                    elif entry['SeqLen'] < self.getInt('MinLen'): entry['Rating'] = 'LOWQUAL'; continue
                    #i# Any scaffolds with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment.
                    # Scaffolds with 80%+ bases in the low/haploid coverage bins and 95%+ of their length mapped by
                    # PurgeHaplotigs onto another scaffold were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov and entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HPURGE'
                    elif (entry['HapPerc']) >= 50 and entry['TopHitCov'] >= hpurgecov and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    # Any other scaffolds exceeding the mindipcov DipPerc coverage assigned to KEEP for further assignment
                    elif (entry['DipPerc']) >= mindipcov and entry['SelfMedK'] == 1: entry['Rating'] = 'PRIMARY'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'REPEAT'  # 'PRIMARY' ?
                    elif (entry['DipPerc']) >= mindipcov: entry['Rating'] = 'KEEP'  # 'PRIMARY' ?
                    # Any other scaffolds with 80%+ low coverage bases were filtered as Low Coverage.
                    elif (entry['LowPerc']) >= artefactcov: entry['Rating'] = 'LOWCOV'
                    # Any other scaffolds with 80%+ high coverage bases are flagged as COLLAPSED.
                    elif (entry['HighPerc']) >= artefactcov: entry['Rating'] = 'COLLAPSED'
                    # Scaffolds that are mostly Haploid Depth are marked as HAPLOTIGS if mostly present 2+ times, or HAPLOID if present once
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 2: entry['Rating'] = 'HAPLOTIG'
                    elif (entry['HapPerc']) >= 50 and entry['SelfMedK'] == 1: entry['Rating'] = 'HAPLOID'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif entry['HapPerc'] < 50 and (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'
                    # Scaffolds with at least 80% hitting other scaffolds are HAPLOTIG or REPEAT
                    elif entry['TopHitCov'] >= hpurgecov and entry['MaxHitCov'] >= rephitcov: entry['Rating'] = 'HAPRPT'
                    elif entry['TopHitCov'] >= hpurgecov: entry['Rating'] = 'HAPLOTIG'
                    # Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'ARTEFACT'

            ## ~ [7x] ~ Add Nala Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode == 'nala1':
                # Purge Haplotigs analysis - round 1
                #
                # Subreads were re-mapped on to the remaining 837 scaffolds and processed with PurgeHaplotigs v20190612 [REF]
                # (implementing Perl v5.28.0, BEDTools v2.27.1 [REF], R v3.5.3, and SAMTools v1.9 [REF]). Based on the
                # PurgeHaplotigs depth histogram, low-, mid- and high-depth thresholds were set to 5X, 30X and 80X. Any scaffolds
                # with <80% at diploid read depth were identified by PurgeHaplotigs for reassignment. Scaffolds with 80%+ bases in
                # the low/haploid coverage bins and 95%+ of their length mapped by PurgeHaplotigs onto another scaffold were
                # filtered as haplotigs or assembly artefacts. Any other scaffolds with 80%+ low coverage bases were filtered as
                # Low Coverage.
                #
                self.errorLog('purgemode=nala1 not yet implemented!')
            elif purgemode == 'nala':
                for entry in dipdb.entries():
                    #i# Low coverage artefacts (e.g. less than 50% of the scaffold covered by at least X reads)
                    if entry['Median_fold'] < minmedian: entry['Rating'] = 'LOWCOV'; continue
                    elif entry['ScreenPerc'] >= 50: entry['Rating'] = 'CONTAMINATION'; continue
                    elif entry['SeqLen'] < self.getInt('MinLen'): entry['Rating'] = 'LOWQUAL'
                    #i# Any scaffold with 80%+ bases in the low/haploid coverage bins were filtered as haplotigs or
                    #i# assembly artefacts.
                    elif (entry['LowPerc']+entry['HapPerc']) >= artefactcov: entry['Rating'] = 'HPURGE'
                    #i# Scaffolds with 20%+ diploid coverage were marked as retention as probable diploids.
                    elif (entry['DipPerc']) >= 20: entry['Rating'] = 'PRIMARY'
                    #i# Scaffolds with <20% diploid coverage and 50%+ high coverage were marked as probable collapsed repeats.
                    elif (entry['HighPerc']) >= 50: entry['Rating'] = 'COLLAPSED'
                    #i# A single remaining Scaffold marked as JUNK by PurgeHaplotigs (over 80% low/high coverage) was also filtered as a probable artefact.
                    else: entry['Rating'] = entry['PurgeHap']

            ## ~ [7x] ~ Add Nala Rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif purgemode not in ['','none','default','diploidocus','complex','nala']:
                self.errorLog('purgemode={} not recognised! Using purgemode=complex'.format(purgemode))

            ## ~ [7c] ~ Add Top and SecHit Ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in dipdb.entries():
                tentry = dipdb.data(entry['TopHit'],expect=False)
                sentry = dipdb.data(entry['SecHit'],expect=False)
                if tentry: entry['TopRating'] = tentry['Rating']
                if sentry: entry['SecRating'] = sentry['Rating']

            ## ~ [7d] ~ Add Set for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dipdb.addField('Set',evalue='-')
            keep = ['QUALITY','FINAL','CORE','COREHAP','PRIMARY','PRIRPT','HAPLOID','HAPLOTIG','KEEP']
            repeat = ['COLLAPSED', 'REPEAT', 'HAPRPT']
            quarantine = ['HPURGE', 'RPURGE', 'LOWQUAL']
            if self.getBool('Diploidify'):
                keep.append('HPURGE')
                quarantine = ['RPURGE', 'LOWQUAL']
            junk = ['LOWCOV', 'CONTAMINATION','JUNK','ARTEFACT']

            for entry in dipdb.entries():
                if entry['Rating'] in keep: entry['Set'] = 'keep'
                elif entry['Rating'] in repeat: entry['Set'] = 'repeat'
                elif entry['Rating'] in quarantine: entry['Set'] = 'quarantine'
                elif entry['Rating'] in junk: entry['Set'] = 'junk'
                else: self.warnLog('Unknown Rating for %s: %s' % (entry['SeqName'],entry['Rating']))

            ## ~ [7d] ~ Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Will need to rationalise these warnings: so many!
            dipdb.index('TopHit')
            # * All `CONTAMINATION` scaffolds
            if 'CONTAMINATION' in dipdb.index('Rating'):
                self.warnLog('{} scaffolds had 50%+ contamination: {}'.format(rje.iLen(dipdb.index('Rating')['CONTAMINATION']),', '.join(dipdb.index('Rating')['CONTAMINATION'])))
            for entry in dipdb.entries():
                if entry['Set'] != 'junk':
                    # * `MedK` < 1 = Scaffolds where most of the scaffold is not in Illumina data.
                    if entry['MedK'] == 0:
                        self.warnLog('MedK warning: less than half of {} {} has kmerread support.'.format(entry['Rating'],entry['SeqName']))
                    # * `ScreenCov > 0` = Scaffolds with possible vector contamination
                    if entry['ScreenPerc'] > 0:
                        self.warnLog('{} {} has {}% contamination detected'.format(entry['Rating'],entry['SeqName'],entry['ScreenPerc']))
                    # * Covered_percent < 95
                    if entry['Covered_percent'] < covwarning:
                        self.warnLog('{} {} only has {}% read coverage'.format(entry['Rating'],entry['SeqName'],entry['Covered_percent']))
                # * Any `HAPLOID` scaffolds with `Duplicated`
                if entry['Rating'] in ['COREHAP','HAPLOID','HAPLOTIG'] and entry['Duplicated'] > 0:
                    self.warnLog('{} {} has {} Duplicated BUSCO: evidence for redundant haplotig?'.format(entry['Rating'],entry['SeqName'],entry['Duplicated']))
                if entry['Set'] in ['junk','quarantine']:
                    # * Any purged scaffolds with `TOP`
                    if entry['TopNum'] > 0:
                        self.warnLog('{} ({}) {} is the Best PurgeHap hit for {} Scaffolds: {}'.format(entry['Rating'],entry['Set'],entry['SeqName'],entry['TopNum'],', '.join(dipdb.index('TopHit')[entry['SeqName']])))
                #?# Think about adding warning for:
                # * Haploid scaffolds that could be Alternative haplotigs, or they could be Haploid alternative regions,
                # e.g. Sex chromosomes. QUESTION: Should be extract one sex chromosome as an alternative contig?!

            ## ~ [7e] ~ Save Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog('Rating Summary',line='-')
            dipdb.indexReport('Rating')
            self.headLog('Sequence Sets',line='-')
            dipdb.indexReport('Set')
            dipdb.saveToFile(sfdict={'LowPerc':4, 'HapPerc':4, 'DipPerc':4, 'HighPerc':4})
            dipdb.saveToFile(filename='%s.ratings.tdt' % basefile, savefields=['SeqName','SeqLen','ScreenPerc','Class','Rating','Set'])

            for seqset in ['diploidocus','core','repeats','quarantine','junk']:
                rje.backup(self,'%s.%s.fasta' % (basefile,seqset),appendable=False)
                open('%s.%s.fasta' % (basefile,seqset),'w')
            seqx = {'diploidocus':0,'core':0,'repeats':0,'quarantine':0,'junk':0,'diploidify':0}
            for rating in keep:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.diploidocus.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    open('{}.core.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['diploidocus'] += 1; seqx['core'] += 1
            for rating in repeat:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.diploidocus.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    open('{}.repeats.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['diploidocus'] += 1
                    seqx['repeats'] += 1
            for rating in quarantine:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.quarantine.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['quarantine'] += 1
            for rating in junk:
                for entry in dipdb.indexEntries('Rating',rating):
                    (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                    seqname = '{} Diploidocus:{}'.format(seqname,rating)
                    while rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname):
                        oldrate = rje.matchExp('Diploidocus:(\S+) Diploidocus:',seqname)[0]
                        seqname = seqname.replace('Diploidocus:{} Diploidocus:'.format(oldrate),'Diploidocus:{}|'.format(oldrate))
                    open('{}.junk.fasta'.format(basefile),'a').write('>{}\n{}\n'.format(seqname, sequence))
                    seqx['junk'] += 1
            for seqset in ['diploidocus','core','repeats','quarantine','junk']:
                self.printLog('#SEQ','%s sequences saved to %s.%s.fasta' % (rje.iStr(seqx[seqset]),basefile,seqset))
            ### Special diploidify output
            if self.getBool('Diploidify'):
                seqset = 'diploidify'
                rje.backup(self,'%s.%s.fasta' % (basefile,seqset),appendable=False)
                open('%s.%s.fasta' % (basefile,seqset),'w')
                for entry in dipdb.entries():
                    cdata = entry['Class'].split('|')
                    if cdata[1] == 'DIP' or (cdata[0] == 'WEAK' and cdata[3] == 'NON') and entry['Rating'] not in junk + quarantine:
                        (seqname, sequence) = seqlist.getSeq(seqdict[entry['SeqName']])
                        if rje.matchExp('^(\S+)_(\S+)__(\S+)',seqname):
                            sdat = rje.matchExp('^(\S+)_(\S+)__(\S+)',seqname)
                            seqname = '%sX2_%s__%sX2 Diploidify: %s' % (sdat[0],sdat[1],sdat[2],seqname)
                        else:
                            sdat = string.split(seqname)
                            seqname = '%sX2 Diploidify: %s' % (sdat[0],seqname)
                        open('{}.{}.fasta'.format(basefile,seqset),'a').write('>{}\n{}\n'.format(seqname, sequence))
                        seqx[seqset] += 1

            if seqlist.getBool('Summarise'):
                for seqset in ['diploidocus','core','repeats','quarantine','junk','diploidify']:
                    if not seqx[seqset]: continue
                    setfas = '%s.%s.fasta' % (basefile,seqset)
                    seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % setfas,'autofilter=F']
                    rje_seqlist.SeqList(self.log,seqcmd)
            return dipdb
        except:
            self.errorLog('Diploidocus.diploidocusHocusPocus() error')
            return None
#########################################################################################################################
    def purgeHaplotigs(self):   ### Runs purge_haplotigs in a subdirectory, using SC read depths to set parameters.
        '''
        Runs purge_haplotigs in a subdirectory, using SC read depths to set parameters.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = self.baseFile(strip_path=True)
            bamfile = os.path.abspath(self.getStr('BAM'))
            bamstrip = os.path.basename(bamfile)
            seqin = os.path.abspath(self.getStr('SeqIn'))
            #i# Need to run PH in subdirectory to avoid conflicts between runs/cycles
            mydir = os.path.abspath(os.curdir)
            ## ~ [1a] ~ Check need to run purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gencov = '{}.gencov'.format(bamstrip)
            covstats = '{}.purge.coverage_stats.csv'.format(basefile)
            purge = '{}.purge.reassignments.tsv'.format(basefile)
            if rje.checkForFiles(filelist=[gencov,covstats,purge],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.') and not self.force():
                return True
            ## ~ [1b] ~ Establish SC read depth using samtools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # scdepth=INT     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
            scdepth = self.getNum('SCDepth')
            if self.getNum('SCDepth'):
                self.printLog('#SCDEP','Using loaded single copy read depth = {0:.2f}X'.format(scdepth))
            else:
                scdepth = self.genomeSize(scdepth=True)
                self.printLog('#SCDEP','Using BUSCO-derived single copy read depth = {0}X'.format(scdepth))
                if not scdepth: raise ValueError('Failed to establish SC read depth')
                self.setNum({'SCDepth':scdepth})
            ## ~ [1c] ~ Setup purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # phlow=INT       : Low depth cutoff for purge_haplotigs (-l X). Will use SCDepth/4 if zero. [0]
            if self.getInt('PHLow') <= 0: self.setInt({'PHLow': int(float(scdepth)/4.0) })
            phlow = self.getInt('PHLow')
            # phmid=INT       : Middle depth for purge_haplotigs (-m X). Will derive from SCDepth if zero. [0]
            if self.getInt('PHMid') <= 0:
                dupdepth = scdepth/2.0
                self.setInt({'PHMid': int(1.5 * dupdepth) })
            phmid = self.getInt('PHMid')
            # phhigh=INT      : High depth cutoff for purge_haplotigs (-h X). Will use SCDepth x 2 if zero. [0]
            if self.getInt('PHHigh') <= 0:
                self.setInt({'PHHigh': int(scdepth * 2+0.5) })
            phhigh = self.getInt('PHHigh')
            #?# Add checks and warnings of cutoff conflicts
            ## ~ [1d] ~ Make and enter new directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            phdir = 'purge_{}/'.format(basefile)
            rje.mkDir(self,phdir,log=True)
            self.printLog('#PURGE','Running purge_haplotigs in: {}'.format(phdir))
            os.chdir(phdir)

            ### ~ [2] ~ Run purge haplotigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# The -depth setting will be increased from 200 to 2xphhigh if >100
            phcmd1 = 'purge_haplotigs hist -b {} -g {} -t {} -d {}'.format(bamfile,seqin,self.threads(),max(200,2*phhigh))
            if self.needToRemake(gencov,bamfile):
                logline = self.loggedSysCall(phcmd1,append=True)
            #!# Option to update the automatically set cutoffs
            self.printLog('#PHDEP','Low=%dX; Mid=%dX; High=%dX. (SC=%dX)' % (phlow,phmid,phhigh,scdepth))
            phcmd2 = 'purge_haplotigs cov -i {} -l {} -m {} -h {} -o {}.purge.coverage_stats.csv -j 80 -s 80'.format(gencov,phlow,phmid,phhigh,basefile)
            if self.needToRemake(covstats,gencov):
                logline = self.loggedSysCall(phcmd2,append=True,threaded=False)
            else: self.printLog('#NOTE','Reusing existing %s on assumption that cutoffs have not changed' % covstats)
            phcmd3 = 'purge_haplotigs purge -g {} -c {}.purge.coverage_stats.csv -t {} -o {}.purge -a 95'.format(seqin,basefile,self.threads(),basefile)
            if self.needToRemake(purge,covstats):
                logline = self.loggedSysCall(phcmd3,append=True,threaded=True)
            os.chdir(mydir)
            ## ~ [2a] ~ Link output files back to main directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ofile in (gencov,covstats,purge):
                phfile = '{}{}'.format(phdir,ofile)
                os.system('ln {} .'.format(phfile))
            # if not rje.checkForFiles(filelist=[gencov,covstats,purge],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: failed!'):
            #     raise IOError('Cannot find purge_haplotigs output. Check {}'.format(phdir))
            return True
        except:
            self.errorLog('Diploidocus.purgeHaplotigs() error')
            os.chdir(mydir)
            raise
#########################################################################################################################
    def depthTrim(self):    ### Trims ends off contigs based on minimum depth
        '''
        Trims ends off contigs based on minimum depth.
        :return:
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            dephead = ['SeqName','SeqLen','Trim5','Trim3','DepTrim']
            depdb = db.addEmptyTable('deptrim',dephead,['SeqName'],log=self.debugging())
            forker = self.obj['Forker']
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            deptrim = self.getInt('DepTrim')
            mintrim = self.getInt('MinTrim')
            if not deptrim: return False
            self.printLog('#TRIM','Trimming %s+ <%dX depth' % (rje_seqlist.dnaLen(mintrim),deptrim))
            ## ~ [1a] IO names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqin = self.getStr('SeqIn')
            seqlist = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.seqNameDic()
            bamfile = self.getBamFile()
            ## ~ [1b] Temp directory for forked depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)
            ## ~ [1c] Check programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.popen('samtools --version').read():
                self.printLog('#SYS',' '.join(os.popen('samtools --version').read().split()))
            else:
                raise IOError('Cannot open samtools: check installation and/or module loading')

            ### ~ [2] Cycle through and fork out the depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for seqname in seqdict:
                tmpfile = '{}{}.{}.{}.deptrim.tmp'.format(tmpdir,basefile,seqname,depmethod)
                if rje.exists(tmpfile):
                    #?# Add checking for completeness #?#
                    if not self.force() and len(open(tmpfile,'r').readline().split()) > 1: skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                if depmethod == 'mpileup':
                    forker.list['ToFork'].append("samtools view -b -h -F 0x100 %s %s | samtools mpileup -BQ0 - 2> /dev/null | awk '$4 >= %d' | (head -n1 && tail -n1) | (awk '{print $2;}' ORS=\" \" && echo) > %s" % (bamfile,seqname,deptrim,tmpfile))
                else:
                    forker.list['ToFork'].append("samtools view -h -F 0x100 %s %s | samtools depth - | awk '$3 >= %d' | (head -n1 && tail -n1) | (awk '{print $2;}' ORS=\" \" && echo) > %s" % (bamfile,seqname,deptrim,tmpfile))
            self.printLog('#DEPTH','{} sequences queued for forking ({} existing files deleted); {} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [2b] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of depth trim analysis completed.')
                else:
                    try:
                        self.errorLog('Samtools forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Samtools forking did not complete')

            ### ~ [3] Read in and process depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for seqname in seqdict:
                tmpfile = '{}{}.{}.{}.deptrim.tmp'.format(tmpdir,basefile,seqname,depmethod)
                seqlen = seqlist.seqLen(seqdict[seqname])
                try:
                    if not rje.exists(tmpfile):
                        raise IOError('Cannot find {}'.format(tmpfile))
                    #?# Add full processing and calculation of mean coverage?
                    try:
                        trimdata = open(tmpfile,'r').readline().split()
                        trim5 = int(trimdata[0])
                        trim3 = int(trimdata[1])
                        if not trim5 or not trim3:  # Trim entire sequence!
                            trim5 = trim3 = seqlen
                    except:
                        self.errorLog('Something has gone wrong with "%s"' % tmpfile)
                        raise
                except:
                    self.errorLog('Samtools depth result processing error',quitchoice=True)
                    continue
                depdb.addEntry({'SeqName':seqname,'SeqLen':seqlen,'Trim5':trim5,'Trim3':trim3,'DepTrim':'none'})
            ## ~ [3b] Trim Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist.obj['Current'] = None
            fasout = '{}.trim.fasta'.format(basefile)
            rje.backup(self,fasout)
            FASOUT = open(fasout,'w')
            sx = 0.0; stot = seqlist.seqNum(); trimx = 0; dumpx = 0; seqx = 0   #?# Add shortx = 0 warning for a shorter length?
            while seqlist.nextSeq():
                self.progLog('\r#TRIM','Trimming {} sequences for low depth ends: {:.2f}%'.format(rje.iStr(stot),sx/stot)); sx += 100.0
                seqname = seqlist.shortName()
                sequence = seqlist.seqSequence()
                seqlen = len(sequence)
                entry = depdb.data(seqname)
                if seqlen != entry['SeqLen']: raise ValueError('%s length mismatch!' % seqname)
                if entry['Trim5'] == entry['SeqLen']:
                    entry['DepTrim'] = 'dump'
                    entry['Trim5'] = entry['Trim3'] = entry['SeqLen']
                    dumpx += 1; continue
                start = 1
                if entry['Trim5'] >= mintrim: start = entry['Trim5']
                end = seqlen
                if seqlen - entry['Trim3'] >= mintrim: end = entry['Trim3']
                if start == 1 and end == seqlen:
                    FASOUT.write('>%s\n%s\n' % (seqlist.seqName(),sequence)); seqx += 1
                    entry['DepTrim'] = 'keep'
                else:
                    trimx += 1; seqx += 1
                    FASOUT.write('>%sX (DepTrim:%s-%s)\n%s\n' % (seqlist.seqName(),rje.iStr(start),rje.iStr(end),sequence[start-1:end]))
                    entry['DepTrim'] = 'trim'
                entry['Trim5'] = entry['Trim5'] - 1
                entry['Trim3'] = seqlen - entry['Trim3']
            FASOUT.close()
            self.printLog('\r#SAVE','{} sequences output to {}: {} trimmed; {} dumped.'.format(rje.iStr(stot),fasout,rje.iStr(trimx),rje.iStr(dumpx)))
            ## ~ [3c] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if seqlist.getBool('Summarise'):
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % fasout,'autofilter=F']
                rje_seqlist.SeqList(self.log,seqcmd)
            ## ~ [3d] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            depdb.saveToFile()
            if rje.iStr(trimx) or rje.iStr(dumpx): return depdb
            else: return False

        except:
            self.errorLog('Diploidocus.depthTrim() error')
            return None
#########################################################################################################################
    # def joinTables(self,name='',join=[],newkey=[],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True,warnings=True):   ### Makes a new table using join of [(Table,Field[,Fieldlist])]
    #     '''
    #     Makes a new table by joining existing tables, using the given fields, in priority order listed. If a Field from a
    #     latter table already exists, it will be added as "Table_Field" instead. If multiple combinations arise for the
    #     same new key, only the first will be kept.
    #     >> name:str [''] = Name for new table. If not given will become "TableX"
    #     >> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
    #         formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
    #         new table. If Field does not exist, it will be added. (Field may be a Formula.)
    #     >> newkey:list [] = If None, will make a new "AutoID" key field. Else, will use the given Fields.
    #     >> cleanup:bool [True] = If True will delete any Fields generated just to make the join
    #     >> delimit:str ['\t'] = Delimiter to be used to join the key fields
    #     >> empties:bool [True] = Whether to keep entries that do not link to 1+ tables with empty values or delete them.
    #     >> check:bool [False] = Whether to check for entries that are not being joined.
    #     >> keeptable:bool [True] = Whether to add new table to self.list['Tables']
    #     '''
#########################################################################################################################
    ### <7> ### Telomere methods                                                                                        #
#########################################################################################################################
### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#     telomere=X      : Basic telomere sequence for search [TTAGGG]
#     telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
#     teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
# parser.add_argument("-w", "--window", type=int, help="This defines the number of first and last nucleotides that will get scanned for telomeric repeats (default: 50).")
# parser.add_argument("-c", "--cutoff", type=float, help='''A telomere is detected if >= c%% of the first (last) nucleotides are telomeric repeats (default: 50%%).''')
#########################################################################################################################
    def findTelomere(self,sequence):    ### Looks for telomeres in nucleotide sequence using telomere regex
        '''
        Looks for telomeres in nucleotide sequence using telomere regex search of sequence ends. Returns a dictionary of
        whether the ends have telomeres and how much was trimmed off as being Ns (5' and 3').
        Based on https://github.com/JanaSperschneider/FindTelomeres.
        >> sequence:str = DNA sequence to search
        << returns dictionary of {'tel5':T/F,'tel3':T/F,'trim5':INT,'trim3':INT,'tel5len':INT,'tel3len':INT}
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tel_forward, tel_reverse = self.getStrUC('TeloFwd'), self.getStrUC('TeloRev')
            sequence = sequence.upper()
            WINDOW = self.getInt('TeloSize')
            REPEAT_CUTOFF = self.getNum('TeloPerc')
            ## ~ [1a] Terminal N-trimming ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            trim5 = 0
            for index, position in enumerate(sequence):
                if position != 'N':
                    trim5 = index
                    break
            start_of_sequence_withoutNs = trim5
            trim3 = 0
            for index, position in enumerate(reversed(sequence)):
                if position != 'N':
                    trim3 = index
                    break
            end_of_sequence_withoutNs = len(sequence) - trim3

            ### ~ [2] Look for Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Look for telomeric repeats at the start of the sequence ~~~~~~~~~~~~~~~~~~~~~ ##
            tel5 = 0    # Keep cycling through bigger windows until it breaks down
            while start_of_sequence_withoutNs < end_of_sequence_withoutNs:
                telomeric_repeats = re.findall(tel_forward, sequence[start_of_sequence_withoutNs:start_of_sequence_withoutNs+WINDOW])
                # Calculate the % of nucleotides that are part of telomeric repeats
                percent_telomeric_repeats_start = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)
                # If more than half of nucleotides at the start/end are telomeric repeats
                if percent_telomeric_repeats_start >= REPEAT_CUTOFF:
                    tel5 += 1
                    start_of_sequence_withoutNs += WINDOW
                else:
                    break
            telomere_at_start = tel5 > 0
            ## ~ [2b] Look for telomeric repeats at the end of the sequence ~~~~~~~~~~~~~~~~~~~~~ ##
            tel3 = 0
            while end_of_sequence_withoutNs > trim5:
                telomeric_repeats = re.findall(tel_reverse, sequence[(end_of_sequence_withoutNs-WINDOW):end_of_sequence_withoutNs])
                # Calculate the % of nucleotides that are part of telomeric repeats
                percent_telomeric_repeats_end = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)
                if percent_telomeric_repeats_end >= REPEAT_CUTOFF:
                    tel3 += 1
                    end_of_sequence_withoutNs -= WINDOW
                else:
                    break
            telomere_at_end = tel3 > 0

            ## ~ [2c] Calculate total percentage telomeres (does not enforce terminal sequences) ~~ ##
            telperc = 0.0
            if telomere_at_start or telomere_at_end:
                telomeric_repeats = re.findall(tel_forward, sequence) + re.findall(tel_reverse, sequence)
                telperc = 100.0 * sum([len(repeat) for repeat in telomeric_repeats]) / float(len(sequence) - sequence.count('N'))

            #!# Update to be more sophisticated and mark end position
            return {'Tel5':telomere_at_start, 'Tel3':telomere_at_end,
                    'Tel5Len':WINDOW*tel5, 'Tel3Len':WINDOW*tel3,
                    'Trim5':trim5, 'Trim3':trim3, 'TelPerc':telperc}
        except:
            self.errorLog('Diploidocus.findTelomere() error'); raise
#########################################################################################################################
    def findTelomeres(self,save=True,keepnull=False):
        '''
        ## Telomere finding [runmode=telomere]

        Canonical motif is TTAGGG/CCCTAA, but one might see variation, so this searches for regex based on the code at
        https://github.com/JanaSperschneider/FindTelomeres.


        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            #i#teldb = self.db().addEmptyTable('telomeres',['Name','SeqLen','Tel5','Tel3','Tel5Len','Tel3Len','Trim5','Trim3','TelPerc'],['Name'],log=self.debugging())
            telfile = '{}.telomeres.{}'.format(db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and rje.checkForFiles(filelist=[telfile],basename='',log=self.log):
                teldb = db.addTable(telfile,name='telomeres',mainkeys=['Name'])
                teldb.dataFormat({'SeqLen':'int','Trim5':'int','Trim3':'int','TelPerc':'num'})
                return teldb
            forks = self.getInt('Forks')
            ## ~ [1a] ~ Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.exists(self.getStr('SeqIn')):
                raise IOError('Diploidocus Telomere mode needs input assembly (seqin=FILE)')
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','summarise=F','autofilter=F'])
            ## ~ [1b] ~ Results table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            teldb = self.db().addEmptyTable('telomeres',['Name','SeqLen','Tel5','Tel3','Tel5Len','Tel3Len','Trim5','Trim3','TelPerc'],['Name'],log=self.debugging())
            telomeres = []  # List of sequences with telomeres
            tel5 = tel3 = telboth = 0

            ### ~ [2] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tel_forward, tel_reverse = self.getStrUC('TeloFwd'), self.getStrUC('TeloRev')
            self.printLog('#TEL','Forward (5\') telomere sequence: {0}'.format(tel_forward))
            self.printLog('#TEL','Reverse (3\') telomere sequence: {0}'.format(tel_reverse))
            sx = 0.0; stot = seqin.seqNum()
            while seqin.nextSeq():
                self.progLog('\r#TELO','Analysing {} sequences for telomeric repeats: {:.2f}%'.format(rje.iStr(stot),sx/stot)); sx += 100.0
                sname = seqin.shortName()
                sequence = seqin.seqSequence()
                tentry = teldb.addEntry(rje.combineDict({'Name':sname,'SeqLen':len(sequence)},self.findTelomere(sequence)))
                # Add reporting in verbose mode?
                if tentry['Tel5'] or tentry['Tel3']:
                    telomeres.append(sname)
                    if tentry['Tel5']: tel5 += 1
                    if tentry['Tel3']: tel3 += 1
                    if tentry['Tel5'] and tentry['Tel3']: telboth += 1
                    if self.v() > 0 and tentry['Tel5']:
                        self.printLog('\r#TELO','5\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim5']))
                    if self.v() > 0 and tentry['Tel3']:
                        self.printLog('\r#TELO','3\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim3']))
            if telomeres: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences: {} at 5\', {} at 3\' ({} both)'.format(rje.iLen(telomeres),rje.iStr(stot),tel5,tel3,telboth))
            else: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences.'.format(rje.iLen(telomeres),rje.iStr(stot)))

            ### ~ [3] Save Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not keepnull:
                for ekey in list(teldb.dict['Data'].keys()):
                    if ekey not in telomeres: teldb.dict['Data'].pop(ekey)
            if save: teldb.saveToFile(sfdict={'TelPerc':4})

            return teldb
        except:
            self.errorLog('Diploidocus.findTelomeres() error')
            return False
#########################################################################################################################
### End of SECTION II: Diploidocus Class                                                                                #
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
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Diploidocus(mainlog,['dna=T','diploidocus=T']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
