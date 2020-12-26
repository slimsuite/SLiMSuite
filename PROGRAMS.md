# Master SLiMSuite and SeqSuite programs

SLiMSuite and SeqSuite are the master control programs that can run all the main SLiMSuite programs, plus some accessory functions. Nominally, SLiMSuite is focused on short linear motif (SLiM) discovery and anlaysis tools, whilst SeqSuite wraps the genomics and general sequence analysis tools. All the SeqSuite tools can also be run by SLiMSuite. In each case the first commandline argument is parsed as the tool to be run, or `prog=X` can set it explicitly.

Help for the selected tool can be accessed using the `help=T` option. Note that `-h`, `-help` or `help` alone will
trigger the SeqSuite help (this!). As `-help` or `help` will also set `help=T`, these commands will trigger both the
SeqSuite help and the selected program help (unless over-ruled by `help=F`). An explicit `help=T` command will only
trigger the selected program help.


## SLiMSuite: Short Linear Motif analysis Suite

* Citation: Edwards et al. (2020), Methods Mol Biol. 2141:37-72. [PMID: 32696352]
* Old citation: Edwards RJ & Palopoli N (2015): Methods Mol Biol. 1268:89-141. [PMID: 25555723]

SLiMSuite is designed to be a front end for the SLiMSuite set of sequence analysis tools. The relevant tool is given
by the first system command, or selected using `prog=X` (or `program=X`). As much as possible, SLiMSuite will emulate
running that tool from the commandline, adding any matching `X.ini` file to the default commandline options read in
(*before* settings read from slimsuite.ini itself). By default, the SLiMCore tool will be called
(libraries/rje_slimcore.py) and read in commands from slimcore.ini.

Help for the selected tool can be accessed using the `help=T` option. Note that `-h`, `-help` or `help` alone will
trigger the SLiMSuite help (this!). As `-help` or `help` will also set `help=T`, these commands will trigger both the
SLiMSuite help and the selected program help (unless over-ruled by `help=F`). An explicit `help=T` command will only
trigger the selected program help.

Running SLiMSuite should also try importing all the main SLiMSuite modules, testing for download errors etc.

SLiMSuite tools recognised by `prog=X`:

- SLiMCore = rje_slimcore.SLiMCore. SLiMSuite core module with MegaSLiM and UPC functions.
- SLiMFarmer = slimfarmer.SLiMFarmer. SLiMSuite job forking/HPC controller.
- QSLiMFinder = qslimfinder.QSLiMFinder. Query-based Short Linear Motif Finder - de novo SLiM prediction.
- SLiMFinder = slimfinder.SLiMFinder. Short Linear Motif Finder - de novo SLiM prediction.
- SLiMList = rje_slimlist.SLiMList. Short Linear Motif manipulation/filtering module.
- SLiMProb = slimprob.SLiMProb. Short Linear Motif Probability - known SLiM prediction.
- SLiMBench = slimbench.SLiMBench. SLiM discovery benchmarking module.
- SLiMMaker = slimmaker.SLiMMaker. Simple SLiM generation from aligned peptide sequences.
- SLiMMutant = slimmutant.SLiMMutant. Short Linear Motif Mutation Analysis Tool
- SLiMParser - slimparser.SLiMParser. SLiMSuite REST output parsing tool.
- PeptCluster = peptcluster.PeptCluster. Peptide alignment, pairwise distance and clustering tool.

Please also see the SeqSuite documentation for additional utilities, which can be run from SLiMSuite or SeqSuite.


## SeqSuite:  Miscellaneous biological sequence analysis toolkit

SeqSuite is a wrapper for the SeqList tool (`libraries/rje_seqlist.py`), which includes a number of sequence parsing and formatting functions. The other tools and functions wrapped by SeqSuite are:

- Apollo = rje_apollo.Apollo GABLAM wrapper for searching against apollo genomes.
- BLAST = rje_blast_V2.BLASTRun. BLAST+ Control Module.
- BUSCOMP = buscomp.BUSCOMP. BUSCO Compiler and Comparison tool.
- DBase = rje_dbase.DatabaseController. Database downloading and processing.
- Diploidocus = diploidocus.Diploidocus. Diploid genome assembly analysis toolkit.
- Ensembl = rje_ensembl.EnsEMBL. EnsEMBL Processing/Manipulation.
- ExTATIC = extatic.ExTATIC. !!! Development only. Not available in main download. !!!
- FIESTA = fiesta.FIESTA. Fasta Input EST Analysis. Transcriptome annotation/querying.
- GABLAM = gablam.GABLAM. Global Analysis of BLAST Local AlignMents
- GASP = gasp.GASP. Gapped Ancestral Sequence Prediction.
- Genbank = rje_genbank.GenBank. Genbank fetching/parsing module.
- GFF = rje_gff.GFF. GFF File Parser and Manipulator.
- GOPHER = gopher.GOPHER. Generation of Orthologous Proteins from Homology-based Estimation of Relationships.
- HAQESAC = haqesac.HAQESAC. Homologue Alignment Quality, Establishment of Subfamilies and Ancestor Construction.
- MITAB = rje_mitab.MITAB. MITAB PPI parser.
- MultiHAQ = multihaq.MultiHAQ. Multi-Query HAQESAC controller.
- PAF = rje_paf.PAF. Minimap2 PAF to GABLAM format converter.
- PAFScaff = pafsaff.PAFScaff. Pairwise mApping Format reference-based scaffold anchoring and super-scaffolding.
- PAGSAT = pagsat.PAGSAT. Pairwise Assembled Genome Sequence Analysis Tool. (Development only)
- PINGU = pingu_V4.PINGU. Protein Interaction Network & GO Utility.
- PPI = rje_ppi.PPI. Protein-Protein Interaction Module.
- PyDocs = rje_pydocs.PyDoc. Python Module Documentation & Distribution.
- RJE_Seq = rje_seq.SeqList. Fasta file sequence manipulation/reformatting.
- SAAGA = saaga.SAAGA. Summarise, Annotate & Assess Genome Annotations
- SAMPhaser = samphaser.SAMPhaser. Diploid chromosome phasing from SAMTools Pileup format.
- SAMTools = rje_samtools.SAMTools. SAMTools mpileup analysis tool. (Development only)
- SeqList = rje_seqlist.SeqList. Fasta file sequence manipulation/reformatting.
- SeqMapper = seqmapper.SeqMapper. Sequence Mapping Program.
- SMRTSCAPE (/PacBio) = smrtscape.SMRTSCAPE. SMRT Subread Coverage & Assembly Parameter Estimator. (Development only)
- Snapper = snp_mapper.SNPMap. SNV to feature annotations mapping and rating tool. (Development only)
- Taxonomy = rje_taxonomy.Taxonomy. Taxonomy download/conversion tool.
- Tree = rje_tree.Tree. Phylogenetic Tree Module.
- Uniprot = rje_uniprot.Uniprot. Uniprot download and parsing module.
- XRef = rje_xref.XRef. Identifier cross-referencing module.
- Zen - rje_zen.Zen. Random Zen wisdom generator and test code.


# Main SLiMSuite and SeqSuite programs

The following tools are under active use or development. Please report any bugs or requests for improved documentation and/or functionality. Testers for some of the newer tools and functions are most welcome! Get in touch through the [SLiMSuite GitHub community][1]. More information on these tools can also be found through the main [SLiMSuite blog][3] and in the module docstrings, available through the `--help` command or [Edwards Lab Software page][6].

## BADASP: Burst After Duplication with Ancestral Sequence Prediction

* Citation: Edwards & Shields (2005), Bioinformatics 21(22):4190-1. [PMID: [16159912][5]]

[BADASP][8] implements the previously published Burst After Duplication (BAD) algorithm, plus two variants that have been used
successfully in identifying functionally interesting sites in platelet signalling proteins and can identify Type I and
Type II divergence. In addition, several other measures of functional specificity and conservation are calculated and
output in plain text format for easy import into other applications. 

See the [BADSASP Manual][13] for further details.

## BUSCOMP: BUSCO Compiler and Comparison tool

* Citation: Edwards RJ (2019). F1000Research 8:995 (slides) (doi: [10.7490/f1000research.1116972.1](https://doi.org/10.7490/f1000research.1116972.1))
* GitHub: <https://github.com/slimsuite/buscomp>

[BUSCOMP][11] is designed to overcome some of the non-deterministic limitations of BUSCO to:

1. compile a non-redundant maximal set of complete BUSCOs from a set of assemblies, and
2. use this set to provide a "true" comparison of completeness between different assemblies of the same genome
with predictable behaviour.

For each BUSCO gene, BUSCOMP will extract the best "Single Complete" sequence from those available, using the
`full_table_*.tsv` results table and `single_copy_busco_sequences/` directory of hit sequences. BUSCOMP ranks all
the hits across all assemblies by Score and keeps the top-ranking hits. Ties are then resolved by Length, keeping
the longest sequence. Ties for Score and Length will keep an arbitrary entry as the winner. Single complete hits
are given preference over Duplicate hits, even if they have a lower score, because only Single hit have their
sequences saved by BUSCO in the `single_copy_busco_sequences/` directory. This set of predicted gene sequences
forms the "BUSCOMPSeq" gene set.

BUSCOMP uses minimap2 to map BUSCOSeq predicted CDS sequences onto genome/transcriptome assemblies, including
those not included in the original BUSCO compilation. This way, the compiled set of species-specific BUSCO
sequences can also be used to generate a quick-and-dirty assessment of completeness for a new genome assembly.
Hits are converted into percentage coverage stats, which are then used to reclassify the BUSCO gene on the basis
of coverage and identity. BUSCOMP ratings are designed to mimic the original BUSCO ratings but have different
definitions. In addition, two extra classes of low quality hit have been added: "Partial" and "Ghost".

For more details, see the main [BUSCOMP documentation][12].

## CompariMotif: Motif vs Motif Comparison Software

* Citation: Edwards, Davey & Shields (2008), Bioinformatics 24(10):1307-9. [PMID: [18375965][16]]
* Webserver: <http://www.slimsuite.unsw.edu.au/servers/comparimotif.php>

[CompariMotif][15] is a piece of software with a single objective: to take two lists of regular expression protein motifs
(typically SLiMs) and compare them to each other, identifying which motifs have some degree of overlap, and
identifying the relationships between those motifs. It can be used to compare a list of motifs with themselves, their
reversed selves, or a list of previously published motifs, for example (e.g. ELM (http://elm.eu.org/)). CompariMotif
outputs a table of all pairs of matching motifs, along with their degree of similarity (information content) and
their relationship to each other.

See the [CompariMotif Manual][14] for more details.

## Diploidocus: Diploid genome assembly analysis toolkit

* Citation:     Edwards RJ et al. (2020), bioRxiv (doi: [10.1101/2020.11.11.379073](https://doi.org/10.1101/2020.11.11.379073))
* GitHub:       <https://github.com/slimsuite/diploidocus>

[Diploidocus][18] is a sequence analysis toolkit for a number of different analyses related to diploid genome assembly.
The main suite of analyses combines long read depth profiles, short read kmer analysis, assembly kmer analysis,
BUSCO gene prediction and contaminant screening for a number of assembly tasks including contamination identification,
haplotig identification/removal and low quality contig/scaffold trimming/filtering.

In addition, Diploidocus will use mapped long reads and BUSCO single copy read depths for genome size prediction
(`gensize`), and coverage (`regcheck`) or copy number estimation (`regcnv`) for user-defined regions. Diploidocus
also has functions for removing redundancy (`sortnr`), generating a non-redundant pseudo-diploid assembly with primary
and secondary scaffolds from 10x pseudohap output (`diphap`), and creating an in silico diploid set of long reads from two
haploid parents (for testing phasing etc.) (`insilico`).

For more information, please see the [main Diploidocus documentation][17].

## GABLAM: Global Analysis of BLAST Local AlignMents

* Citation: Davey, Shields & Edwards (2006), Nucleic Acids Res. 34(12):3546-54. [PMID: 16855291]
* Manual: <http://bit.ly/GABLAMManual>

[GABLAM][19] takes one or two sequence datasets, peforming an intensive All by All BLAST or Minimap2 search and then tabulating
the results as a series of local hits, pairwise comparisons, and hit summaries per query. GABLAM can also generate a SNP table based on local hits, restricted by default to unique "best" hits be query region. All hits can be filtered by length and percentage identity at the local or global level. 

Hit sequences can be output as fasta files, either as full-length sequences or restricted to the local hit fragments. This output can be per query sequence, or combined for the whole query dataset. SAM and GFF3 outputs for local hits are also available.

See the [GABLAM manual](http://bit.ly/GABLAMManual) and [GABLAM docs][19] for details.

## GASP: Gapped Ancestral Sequence Prediction

* Citation: Edwards & Shields (2004), BMC Bioinformatics 5(1):123. [PMID: 15350199]

GASP runs the GASP gapped ancestral sequence prediction algorithm that forms part of HAQESAC as a standalone program. 

See the [GASP Manual](https://github.com/slimsuite/SLiMSuite/blob/master/docs/manuals/GASP%20Manual.pdf) for more details.

## GOPHER: Generation of Orthologous Proteins from Homology-based Estimation of Relationships

* Citation:     Davey, Edwards & Shields (2007), Nucleic Acids Res. 35(Web Server issue):W455-9. [PMID: 17576682]
* Manual:       <http://bit.ly/GOPHERManual>

GOPHER is designed to take in two sequences files and generate datasets of orthologous sequence alignments.
The first `seqin` sequence set is the 'queries' around which orthologous datasets are to be assembled. This is now
optimised for a dataset consisting of one protein per protein-coding gene, although splice variants should be dealt
with OK and treated as paralogues. This will only cause problems if the postdup=T option is used, which restricts
orthologues returned to be within the last post-duplication clade for the sequence.

The second `orthdb` is the list of proteins from which the orthologues will be extracted. The `seqin` sequences are
then BLASTed against the orthdb and processed (see below) to retain putative orthologues using an estimation of the
phylogenetic relationships based on pairwise sequences similarities.

Please cleanup the input data into a desired non-redundant dataset before running GOPHER. Duplicate accession numbers will not be tolerated by GOPHER and (arbitrary) duplicates will be deleted if the sequences are the same, or renamed otherwise. Renaming may cause problems later. It is highly desirable not to have two proteins with the same
accession number but different amino acid sequences. Note that unknown species are also not permitted. Version 3.0 has improved directory organisation for multi-species and multi-orthdb GOPHER runs on the same system.

See the [GOPHER Manual](http://bit.ly/GOPHERManual) for more information.

## HAQESAC: Homologue Alignment Quality, Establishment of Subfamilies and Ancestor Construction

* Citation: Edwards et al. (2007), Nature Chem. Biol. 3(2):108-112. [PMID: 17220901]

HAQESAC is a tool designed for processing a dataset of potential homologues into a trustworthy gobal alignment and well
bootstrap-supported phylogeny. By default, an intial dataset consisting of homologues (detected by BLAST, for example)
is processed into a dataset consisting of the query protein and its orthologues in other species plus paralogous
subfamilies in the same gene family.

Individual sequences are therefore screened fulfil two criteria:
   1. They must be homologous to (and alignable with) the query sequence of interest.
   2. They must be a member of a subfamily within the gene family to which the query sequence belongs. 

HAQ. The first stage of data cleanup is therefore to remove rogue sequences that either do not fit in the gene family at
all or are too distantly related to the query protein for a decent alignment that can be used for useful further
analysis. This is achieved firstly by a simple identity cut-off, determined by pairwise alignments of sequences, and
then by a more complex procedure of removing sequences for whom the overall alignability is poor. During this procedure,
sequences that have too many gaps are also removed as too many gapped residues can cause problems for downstream
evolutionary analyses. Further screening is achieved based on phylogenetic information.

ES. Once the dataset has been 'cleaned up' (and, indeed, during processing), HAQESAC can be used to assign sequences to
subgroups or subfamilies, if such information is needed for downstream analyses.

AC. The final step that HAQESAC is able to perform is ancestral sequence prediction using the GASP (Gapped Ancestral
Sequence Prediction) algorithm (Edwards & Shields 2005)

By default, HAQESAC will perform all these operations. However, it is possible to turn one or more off and only, for
example, reject individually badly aligned sequences, if desired. Details can be found in the [HAQESAC Manual](https://github.com/slimsuite/SLiMSuite/blob/master/docs/manuals/HAQESAC%20Manual.pdf).

## MultiHAQ: Multi-Query HAQESAC controller

* Citation: Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: [20924652][9]]

MultiHAQ is a wrapper for multiple HAQESAC runs where different query proteins are to be BLASTed against the same
search database(s) and run through HAQESAC with the same settings. The default expectation is that some queries will
be returned by the HAQESAC runs of other queries and may therefore be skipped as a result, although this can be
switched off using screenqry=F. For large runs, the first phase of MulitHAQ will take a long time to run. In these
cases, it may be desirable to set the second, interactive, phase running before it has finished. This is achieved
using the "chaser" option, which will set the second phase in motion, "chasing" the progress of the first. To avoid
jumbled log output, this should be given a different log file using log=FILE.

Note: that all options will be output into a haqesac.ini file in the haqdir path, for both HAQESAC runs within the
framework of MultiHAQ itself and also for later runs using the batch file produced. Any generic HAQESAC options
should therefore be placed into a multihaq.ini file, not a haqesac.ini file and multiple runs with different settings
using the same haqdir should be avoided.

Note: Because HAQESAC makes use of `RJE_SEQ` filtering options, they will NOT be applied to the MultiHAQ query input
file prior to analysis. To filter this input, run it through `rje_seq.py` separately in advance of running multihaq.

## PAFScaff: Pairwise mApping Format reference-based scaffold anchoring and super-scaffolding.

* Citation: Field et al. (2020), GigaScience 9(4):giaa027. [PMID: [32236524][22]]
* GitHub: <https://github.com/slimsuite/pafscaff>

[PAFScaff][20] is designed for mapping genome assembly scaffolds to a closely-related chromosome-level reference genome
assembly. It uses (or runs) [Minimap2](https://github.com/lh3/minimap2) to perform an efficient (if rough) all-
against-all mapping, then parses the output to assign assembly scaffolds to reference chromosomes.

Mapping is based on minimap2-aligned assembly scaffold ("Query") coverage against the reference chromosomes.
Scaffolds are "placed" on the reference scaffold with most coverage. Any scaffolds failing to map onto any
chromosome are rated as "Unplaced". For each reference chromosome, PAFScaff then "anchors" placed assembly
scaffolds starting with the longest assembly scaffold. Each placed scaffold is then assessed in order of
decreasing scaffold length. Any scaffolds that do not overlap with already anchored scaffolds in terms of the
Reference chromosome positions they map onto are also considered "Anchored". if `newprefix=X` is set, scaffolds
are renamed with the Reference chromosome they match onto. The original scaffold name and mapping details are
included in the description. Unplaced scaffolds are not renamed.

Finally, Anchored scaffolds are super-scaffolded by inserting gaps of `NnNnNnNnNn` sequence between anchored
scaffolds. The lengths of these gaps are determined by the space between the reference positions, modified by
overhanging query scaffold regions (min. length 10). The alternating case of these gaps makes them easy to
identify later.

For details, see the [PAFScaff documentation][21].

## PAGSAT: Pairwise Assembled Genome Sequence Analysis Tool

* Citation: Pérez-Bercoff et al. (2015), [F1000Research 4:1022 (poster)][23].

[PAGSAT][24] is for the assessment of an assembled genome versus a suitable reference. For optimal results, the
reference genome will be close to identical to that which should be assembled. However, comparative analyses should
still be useful when different assemblies are run against a related genome - although there will not be the same
expectation for 100% coverage and accuracy, inaccuracies would still be expected to make an assembly less similar
to the reference.

For more details, see the [PAGSAT documentation][24] and [Yeast SMRT sequencing poster][23].

## PeptCluster: Peptide Clustering Module

* Webserver: <http://bioware.soton.ac.uk/peptcluster.html>

PeptCluster is for simple sequence-based clustering of short (aligned) peptide sequences. First, a pairwise distance
matrix is generated from the peptides. This distance matrix is then used to generate a tree using a distance method
such as Neighbour-Joining or UPGMA.

Default distances are amino acid property differences loaded from an amino acid property matrix file.

Version 1.5.0 incorporates a new peptide alignment mode to deal with unaligned peptides. This is controlled by the
`peptalign=T/F/X` option, which is set to True by default. If given a regular expression, this will be used to guide
the alignment. Otherwise, the longest peptides will be used as a guide and the minimum number of gaps added to
shorter peptides. Pairwise peptide distance measures are used to assess different variants, starting with amino acid
properties, then simple sequence identity (if ties) and finally PAM distances. One of the latter can be set as
the priority using `peptdis=X`. Peptide alignment assumes that peptides have termini (^ & $) or flanking wildcards
added. If not, set `termini=F`.

## QSLiMFinder: Query Short Linear Motif Finder

* Citation: Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]
* SLiMFinder: Edwards, Davey & Shields (2007), PLoS ONE 2(10): e967. [PMID: 17912346]
* Webserver: <http://www.slimsuite.unsw.edu.au/servers/qslimfinder.php>
* Manual: <http://bit.ly/SFManual>

QSLiMFinder is a modification of the basic SLiMFinder tool to specifically look for SLiMs shared by a query sequence
and one or more additional sequences. To do this, SLiMBuild first identifies all motifs that are present in the query
sequences before removing it (and its UPC) from the dataset. The rest of the search and stats takes place using the
remainder of the dataset but only using motifs found in the query. The final correction for multiple testing is made
using a motif space defined by the original query sequence, rather than the full potential motif space used by the
original SLiMFinder. This is offset against the increased probability of the observed motif support values due to the
reduction of support that results from removing the query sequence but could potentially still identify SLiMs will
increased significance.

Note that minocc and ambocc values *include* the query sequence, e.g. minocc=2 specifies the query and ONE other UPC. 

See the [SLiMFinder Manual](http://bit.ly/SFManual) for more details.

## SAAGA: Summarise, Annotate & Assess Genome Annotations

* Citation: Edwards RJ et al. (2020), bioRxiv <https://doi.org/10.1101/2020.11.11.379073>
* GitHub: <http://github.com/slimsuite/saaga>

[SAAGA][25] is a tool for summarising, annotating and assessing genome annotations, with a particular focus on annotation
generated by GeMoMa. The core of SAAGA is reciprocal MMeqs searches of the annotation and reference proteomes. These
are used to identify the best hits for protein product identification and to assess annotations based on query and
hit coverage. SAAGA will also generate annotation summary statistics, and extract the longest protein from each gene
for a representative non-redundant proteome (e.g. for BUSCO analysis).

See the [SAAGA documentation][26] for more information.

## SAMPhaser: Diploid chromosome phasing from SAMTools Pileup format.

* Citation: Song, Thomas & Edwards (2019), [Marine Genomics 48:100687](https://doi.org/10.1016/j.margen.2019.05.002).
* Documentation: <https://github.com/slimsuite/SLiMSuite/wiki/SAMPhaser>

[SAMPhaser][27] is a tool designed to take an input of long read (e.g. PacBio) data mapped onto a genome assembly and
phase the data into haplotype blocks before "unzipping" the assembly into phased "haplotigs". Unphased regions are
also output as single "collapsed" haplotigs. This is designed for phasing PacBio assemblies of diploid organisms.
By default, only SNPs are used for phasing, with indel polymorphisms being ignored. This is because indels are more
likely to be errors. In particular, mononucleotide repeats could have indels that look like false well-supported
polymorphisms.

See the [SAMPhaser Documentation][27] for more details.

## SeqMapper: Sequence Mapping Program

SeqMapper is for mapping one set of protein sequences onto a different sequence database, using Accession Numbers
etc where possible and then using GABLAM when no direct match is possible. The program gives the following outputs:

- *.*.mapped.fas = Fasta file of successfully mapped sequences
- *.*.missing.fas = Fasta file of sequences that could not be mapped
- *.*.mapping.tdt = Delimited file giving details of mapping (Seq, MapSeq, Method)

See the SeqMapper documentation for more details.

## SLiMBench: Short Linear Motif prediction Benchmarking

* Citation: Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]

SLiMBench has two primary functions:

1. Generating SLiM prediction benchmarking datasets from ELM (or other data in a similar format). This includes
options for generating random and/or simulated datasets for ROC analysis etc.

2. Assessing the results of SLiM predictions against a Benchmark. This program is designed to work with SLiMFinder
and QSLiMFinder output, so some prior results parsing may be needed for other methods.

If `generate=F benchmark=F`, SLiMBench will check and optionally download the input files but perform no additional
processing or analysis.

Please see the SLiMBench manual for more details.

## SLiMFarmer: SLiMSuite HPC and multi-CPU job farming control program

* Citation: Edwards et al. (2020), Methods Mol Biol. 2141:37-72. [PMID: [32696352][28]]

SLiMFarmer is designed to control and execute parallel processing jobs on an HPC cluster using PBS and QSUB. If
qsub=T it will generate a job file and use qsub to place that job in the queue using the appropriate parameter
settings. If `slimsuite=T` and `farm=X` gives a recognised program (below) or `hpcmode` is not `fork` then the qsub
job will call SLiMFarmer with the same commandline options, plus `qsub=F i=-1 v=-1`. If `seqbyseq=T`, this will be
run in a special way. (See SeqBySeq mode.) Otherwise `slimsuite=T` indicates that `farm=X` is a SLiMSuite program,
for which the python call and `pypath` will be added. If this program uses forking then it should parallelise over a
single multi-processor node. If `farm=X` contains a `/` path separator, this will be added to `pypath`, otherwise it
will be assumed that `farm` is in `tools/`.

If `slimsuite=F` then farm should be a program call to be queued in the PBS job file instead. In this case, the
contents of jobini=FILE will be added to the end of the farm call. This enables commands with double quotes to be
included in the farm command, for example.

Currently recognised SLiMSuite programs for farming: SLiMFinder, QSLiMFinder, SLiMProb, SLiMCore.

Currently recognised SLiMSuite programs for rsh mode only qsub farming: GOPHER, SLiMSearch, UniFake.

NOTE: Any commandline options that need bracketing quotes will need to be placed into an ini file. This can either
be the ini file used by SLiMFarmer, or a `jobini=FILE` that will only be used by the farmed programs. Note that
commands in `slimfarmer.ini` will not be passed on to other SLiMSuite programs unless `ini=slimfarmer.ini` is given
as a commandline argument.

The runid=X setting is important for SLiMSuite job farming as this is what separates different parameter setting
combinations run on the same data and is also used for identifying which datasets have already been run. Running
several jobs on the same data using the same SLiMSuite program but with different parameter settings will therefore
cause problems. If runid is not set, it will default to the job=X setting.

The hpcmode=X setting determines the method used for farming out jobs across the nodes. hpcmode=rsh uses rsh to spawn
the additional processes out to other nodes, based on a script written for the IRIDIS HPC by Ivan Wolton.
hpcmode=fork will restrict analysis to a single node and use Python forking to distribute jobs. This can be used even
on a single multi-processor machine to fork out SLiMSuite jobs. basefile=X will set the log, RunID, ResFile, ResDir
and Job: RunID and Job will have path stripped; ResFile will have .csv appended.

Initially, it will call other programs but, in time, it is envisaged that other programs will make use of SLiMFarmer
and have parallelisation built-in.

## SLiMFinder: Short Linear Motif Finder

* Citation: Edwards RJ, Davey NE & Shields DC (2007), PLoS ONE 2(10): e967. [PMID: 17912346]
* ConsMask Citation: Davey NE, Shields DC & Edwards RJ (2009), Bioinformatics 25(4): 443-50. [PMID: 19136552]
* SigV/SigPrime Citation: Davey NE, Edwards RJ & Shields DC (2010), BMC Bioinformatics 11: 14. [PMID: 20055997]
* SLiMScape/REST Citation: Olorin E, O'Brien KT, Palopoli N, Perez-Bercoff A & Shields DC, Edwards RJ (2015), F1000Research 4:477.
* SLiMMaker Citation: Palopoli N, Lythgow KT & Edwards RJ (2015), Bioinformatics 31(14): 2284-2293. [PMID: 25792551]
* Webserver: <http://www.slimsuite.unsw.edu.au/servers/slimfinder.php>
* Manual: <http://bit.ly/SFManual>

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

Please see [SLiMFinder Manual][29] for more details.

## SLiMMaker: SLiM generator from aligned peptide sequences

* Citation: Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]
* Webserver: http://www.slimsuite.unsw.edu.au/servers/slimmaker.php

SLiMMaker has a fairly simple function of reading in a set of sequences and generating a regular expression motif
from them. It is designed with protein sequences in mind but should work for DNA sequences too. Input sequences can
be in fasta format or just plain text (with no sequence headers) and should be aligned already. If varlength=F then
gapped positions will be ignored (treated as Xs) and variable length wildcards are not returned. If varlength=T, any
gapped positions will be assessed based on the ungapped peptides at that position and a variable length inserted.
This variable-length position may be a wildcard or it may be a defined position if there is sufficient signal in the
peptides with amino acids at that position.

See the SLiMMaker docs for more details.

## SLiMParser: SLiMSuite REST output parsing tool.

* Citation: Edwards et al. (2020), Methods Mol Biol. 2141:37-72. [PMID: [32696352][28]]

This module is for parsing the full REST output of a program into a couple of dictionaries, with options to output
the data to files or convert to/from json format.

If `restin=FILE` is a URL, this will be interpreted as a REST command for API access. Use with `rest=X` and
`pureapi=T` to print the output to STDOUT once the run is complete. Use in conjunction with `v=-1` to avoid
additional STDOUT output and `silent=T` to avoid log generation.

REST URLs can include files to be uploaded. These must be prefixed with `file:`, e.g. `&seqin=file:input.fas`. If the
specified file exists then the content will replace the file name in the REST call.

## SLiMProb: Short Linear Motif Probability tool

* Citation: Davey, Haslam, Shields & Edwards (2010), Lecture Notes in Bioinformatics 6282: 50-61.
* Webserver: <http://www.slimsuite.unsw.edu.au/servers/slimprob.php>
* Manual: <http://bit.ly/SProbManual>

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

## SLiMSearch: Short Linear Motif Search tool

* Citation: Davey, Haslam, Shields & Edwards (2010), Lecture Notes in Bioinformatics 6282: 50-61. 

SLiMSearch is a tool for finding pre-defined SLiMs (Short Linear Motifs) in a protein sequence database. SLiMSearch
can make use of corrections for evolutionary relationships and a variation of the SLiMChance alogrithm from
SLiMFinder to assess motifs for statistical over- and under-representation. SLiMSearch is a replacement for PRESTO
and uses many of the same underlying modules.

Benefits of SLiMSearch that make it more useful than a lot of existing tools include:

* searching with mismatches rather than restricting hits to perfect matches.
* optional equivalency files for searching with specific allowed mismatched (e.g. charge conservation)
* generation or reading of alignment files from which to calculate conservation statistics for motif occurrences.
* additional statistics, including protein disorder, surface accessibility and hydrophobicity predictions
* recognition of "n of m" motif elements in the form `<X:n:m>`, where X is one or more amino acids that must occur n+
times across which m positions. E.g. `<IL:3:5>` must have 3+ Is and/or Ls in a 5aa stretch.

Main output for SLiMSearch is a delimited file of motif/peptide occurrences but the motifaln=T and proteinaln=T also
allow output of alignments of motifs and their occurrences. The primary outputs are named *.csv for the occurrence
data and *.summary.csv for the summary data for each motif/dataset pair. 

NOTE: SLiMSearch has now been largely superseded by SLiMProb for motif statistics.

## Snapper: Genome-wide SNP Mapper

Snapper is designed to generate a table of SNPs from a BLAST comparison of two genomes, map those SNPs onto genome
features, predict effects and generate a series of output tables to aid exploration of genomic differences.

For more details, see the Snapper docs.

## SynBad: Synteny-based scaffolding adjustment

* GitHub: <https://github.com/slimsuite/synbad>

[SynBad][30] is a tool for comparing two related genome assemblies and identify putative translocations and inversions
between the two that correspond to gap positions. These positions could indicate misplaced scaffolding. SynBad can also 
fragment assemblies on gaps that are not syntenic, unless more than `minreadspan=INT` reads span the gap.

For more details, see the [SynBad documentation][30].

---

# Unsupported SLiMSuite/SeqSuite programs

The following tools have not been functionally updated and/or tested for several years. They should hopefully still be functional but are not currently under ongoing use or development. As a consequence, they have limited support. That said, if something has broken, feel free to [raise an issue][2] or get in touch through the [SLiMSuite GitHub community][1].

## APHID:  Automated Processing of High-resolution Intensity Data

* Citation: Raab et al. (2010), Proteomics 10: 2790-2800. [PMID: [20486118][4]]

[APHID][7] takes for input the partially processed results of MS analysis, with intensity data, filters based on
scores thresholds, removes redundancy (using PINGU) and calculates relative intensity scores. PINGU is then used to
generate outputs for use with Cytoscape and other visualisation tools.

Input takes the form of a delimited text file with the following column headers: Expt, Subpop, Identifier,
logInt, Score. In addition, other columns may be present (and may be used to filter data). The "unique" column headers allow
individual identifications to be isolated, which is important for data filtering and intensity mapping.

Intensities are converted into relative intensities with two options: (1) the redundancy level determines which
results are combined. This may be simply at the identified peptide level, the gene level, or even at the protein
family level (as determined through BLAST homology). If "fam" is used then GABLAM will be used to generate families
for a given level of identity.


## BUDAPEST: Bioinformatics Utility for Data Analysis of Proteomics on ESTs

* Citation: Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: [20924652][9]]

Proteomic analysis of EST data presents a bioinformatics challenge that is absent from standard protein-sequence
based identification. EST sequences are translated in all six Reading Frames (RF), most of which will not be
biologically relevant. In addition to increasing the search space for the MS search engines, there is also the added
challenge of removing redundancy from results (due to the inherent redundancy of the EST database), removing spurious
identifications (due to the translation of incorrect reading frames), and identifying the true protein hits through
homology to known proteins.

BUDAPEST (Bioinformatics Utility for Data Analysis of Proteomics on ESTs) aims to overcome some of these problems by
post-processing results to remove redundancy and assign putative homology-based identifications to translated RFs
that have been "hit" during a MASCOT search of MS data against an EST database. Peptides assigned to "incorrect" RFs
are eliminated and EST translations combined in consensus sequences using FIESTA (Fasta Input EST Analysis). These
consensus hits are optionally filtered on the number of MASCOT peptides they contain before being re-annotated using
BLAST searches against a reference database. Finally, HAQESAC can be used for automated or semi-automated phylogenetic
analysis for improved sequence annotation.

## FIESTA: Fasta Input EST Analysis

* Citation: Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: [20924652][9]]

FIESTA has three primary functions: <1>	Discovery, assembly and evolutionary analysis of candidate genes in an EST
library; <2> Assembly of an EST library for proteomics analysis; <3> Translation/Annotation of an EST library for
proteomics analysis. 

See FIESTA docs for more details.

## GFESSA: Genome-Free EST SuperSAGE Analysis

* Citation: Johansson et al. (2020), [Algal Research 48:101917](https://www.sciencedirect.com/science/article/pii/S2211926420300369).

GFESSA is for the automated processing, mapping and identification-by-homology for SuperSAGE tag data for
organisms without genome sequences, relying predominantly on EST libraries etc. Although designed for genome-free
analysis, there is no reason why transcriptome data from genome projects cannot be used in the pipeline. 

GFESSA aims to take care of the following main issues:
1. Removal of unreliable tag identification/quantification based on limited count numbers.
2. Converting raw count values into enrichment in one condition versus another.
3. Calculating mean quantification for genes based on all the tags mapping to the same sequence.
4. The redundancy of EST libraries, by mapping tags to multiple sequences where necessary and clustering sequences
on shared tags.

The final output is a list of the sequences identified by the SAGE experiment along with enrichment data and
clustering based on shared tags.

## HAPPI: Hyperlinked Analysis of Protein-Protein Interactions

* Citation: Edwards et al. (2011), Molecular Biosystems DOI: 10.1039/C1MB05212H. [PMID: 21879107]

HAPPI makes a set of linked webpages for the quick and dirty analysis of PPI networks around one or more
sets of genes. A table of genes and their corresponding classes is used to make the initial front page tables.
Additional data tables can also be loaded and linked via the "Gene" field for individual Gene Pages. These must
be named in the format X.N.tdt, where N will be used as the tab name.

## PICSI: Proteomics Identification from Cross-Species Inference

* Citation: Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: [20924652][9]]

This module is for cross-species protein identifications using searches against NCBInr, for example. MASCOT
processing uses BUDAPEST. Hits are then converted into peptides for redundancy removal. Hits from a given (known)
query species are preferentially kept and any peptides belonging to those hits are purged from hits returned by
other species. All hits are then classified:
- UNIQUE    : Contains 2+ peptides, including 1+ unique peptides
- NR        : Contains 2+ peptides; None unique but 1+ peptides only found in other "NR" proteins
- REDUNDANT : Contains 2+ peptides but all found in proteins also containing UNIQUE peptides
- REJECT    : Identified by <2 peptides once query-species peptides subtracted

## PINGU: Protein Interaction Network & GO Utility

PINGU (Protein Interaction Network & GO Utility) is designed to be a general utility for Protein Protein Interaction
(PPI) and Gene Ontology (GO) analysis. Earlier versions of PINGU contained a lot of the code for processing PPI and
GO data, which have subsequently been moved to `rje_ppi.py` and `rje_go.py` libraries.

## PRESTO: Protein Regular Expression Search Tool

Note: This program has been superceded in most functions by SLiMSearch.

PRESTO is what the acronym suggests: a search tool for searching proteins with peptide sequences or motifs using an
algorithm based on Regular Expressions. The simple input and output formats and ease of use on local databases make
PRESTO a useful alternative to web resources for high throughput studies.

The additional benefits of PRESTO that make it more useful than a lot of existing tools include:

* PRESTO can be given alignment files from which to calculate conservation statistics for motif occurrences.
* searching with mismatches rather than restricting hits to perfect matches.
* additional statistics, inlcuding protein disorder, surface accessibility and hydrophobicity predictions
* production of separate fasta files containing the proteins hit by each motif.
* production of both UniProt format results and delimited text results for easy incorporation into other applications.
* inbuilt tandem Mass Spec ambiguities. 

PRESTO recognises "n of m" motif elements in the form `<X:n:m>`, where X is one or more amino acids that must occur n+
times across which m positions. E.g. `<IL:3:5>` must have 3+ Is and/or Ls in a 5aa stretch.

Main output for PRESTO is a delimited file of motif/peptide occurrences but the motifaln=T and proteinaln=T also allow
output of alignments of motifs and their occurrences. PRESTO has an additional motinfo=FILE output, which produces a
summary table of the input motifs, inlcuding Expected values if searchdb given and information content if motifIC=T.
Hit proteins can also be output in fasta format (fasout=T) or UniProt format with occurrences as features (uniprot=T).
    
## SLiMMutant: Short Linear Motif Mutation Analysis Tool

SLiMMutant is a Short Linear Motif Mutation Analysis Tool, designed to identify and assess mutations that create
and/or destroy SLiMs. There are three main run modes:

- Mode 1. Generating mutant datasets for analysis [`generate=T`]

The main input is: (1) a file of protein sequence mutations [`mutfile=FILE`] in a delimited text format with aa
substitution [`mutfield=X`] and protein accession number [`protfield=X`] data; (2) a corresponding sequence file
[`seqin=FILE`]; (3) a file of SLiMs to analyse [`motifs=FILE`]. This will process the data and generate two sequence
files: `*.wildtype.fas` and `*.mutant.fas`. These files will be named after the input `mutfile` unless `basefile=X`
is used.

- Mode 2. Run SLiMProb on datasets [`slimprob=T`]

This will run SLiMProb on the two datasets, once per `*.ini` file given by `slimini=LIST`. These runs should have
distinct `runid=X` settings. If no `*.ini` files are given, as single run will be made using commandline settings.

- Mode 3. Compile results of SLiMProb runs. [`analyse=T`]

This will compare the `*.wildtype.fas` and `*.mutant.fas` results from the `*.occ.csv` file produced by SLiMProb.
All mutations analysed will be identified from `*.mutant.fas`. SLiM occurrences are then matched up between wildtype
and mutant versions of the same sequence. If none of the mutations have effected the SLiM prediction, then the
wildtype and all mutant sequences will return the motif. If, on the other hand, mutations have created/destroyed
motifs, occurrences will be missing from the wildtype and/or 1+ mutant sequences. All unaffected SLiM instances are
first removed and altered SLiM instances  output to `*.MutOcc.csv`. Differences between mutants and wildtypes are
calculated for each `RunID`-`Motif` combination and summary results output to `*.Mut_vs_WT.csv`. If `motlist=LIST` is
given, analysis is restricted to a subset of motifs.

Unless `basefile=FILE` is given, output files will be named after `mutfile=FILE` but output into the current run
directory. If running in batch mode, basefile cannot be used.

NOTE: SLiMMutant is still in development and has not been thoroughly tested or benchmarked.

## SMRTSCAPE: SMRT Subread Coverage & Assembly Parameter Estimator

`SMRTSCAPE` (SMRT Subread Coverage & Assembly Parameter Estimator) is a tool for analysis of PacBio raw sequencing
data to assist the design and execution of PacBio genomics projects. It has a number of functions concerned with
predicting and/or assessing the quantity and quality of useable data produced:

1. **Genome Coverage (`coverage=T`).** This method tries to predict genome coverage and accuracy for different depths of
PacBio sequencing. This is useful for estimating genome coverage and/or required numbers of SMRT cells from predicted
read outputs or emprical (average) SMRT cell data if the `BASEFILE.unique.tdt` output (generated by `summary=T`, below) is
found. NOTE: Default settings for SMRT cell output are not reliable and you should speak to your sequencing provider
for their up-to-date figures. By default, output for this mode is incremented by `XCoverage` but this can be switched
to numbers of SMRT cells with `bysmrt=T`. 

2. **Summarise subreads (`summarise=T`).** This function summarises subread data from a given `seqin=FILE` fasta
file, or a set of subread fasta files given with `batch=FILELIST` (or listed in `*.fofn`). This produces sequence
summary data (read lengths, N50 etc.) for each sequence file, SMRT cell and the combined dataset (`*.summary.tdt`).
In addition, tables are generated that summarise each read individually, which can then be used for further read
filtering or calculations. Summaries are produced for _all_ data (`*.zmw.tdt`) and just the **Unique** subread data,
which is the _longest_ read from each ZMW (`*.unique.tdt`). FALCON only uses unique read data, and so it is these
data that are used for the rest of SMRTSCAPE functions. A summary of **Read Quality (RQ)** data is also output
(`*.rq.tdt`).

3. **Calculate length cutoffs (`calculate=T`).** Calculates length cutoffs for different XCoverage and RQ
combinations from subread summary data. Generates `*.cutoffs.tdt`.

4. **Optimise (`optimise=T`).** This function attempts to generate predicted optimum assembly settings from the
`summary=T` and `calculate=T` table data. NOTE: In `V1.x`, this option was `parameters=T`.

5. **Filter subreads (`readfilter=T`).** This filters *unique* subreads into a new fasta file (`*.LXXXRQXXX.fasta`)
based on min. read quality (`rqfilter=X`) and min. read length (`lenfilter=X`). NOTE: These filters are not
available in FALCON, so sequence input must be pre-filtered in this way.

6. **Preassembly fragmentation analysis (`preassembly=FILE`).** Processes a Preassembly fasta file to assess/correct
over-fragmentation. Corrected preassembly reads are output to `*.smrtscape.fas` and summary data output to
`*.fragment.tdt`.

7. **Feature coverage (`ftxcov=INTLIST`).** Calculates full-length coverage for a list of feature lengths, as well as
their probability of detection (in raw subread data) if present at different population frequencies
(`ftxfreq=NUMLIST`).  If seqin=FILE is given, this will be used directly unless summarise=T, calculate=T or
optimise=T. Otherwise, unique reads from (`*.unique.tdt`) will be used.

## UniFake: Fake UniProt DAT File Generator

This program runs a number of in silico predication programs and converts protein sequences into a fake UniProt DAT
flat file. Additional features may be given as one or more tables, using the features=LIST option. Please see the
UniFake Manual for more details. 



---

# Extra SeqSuite tools and functions

_Details of additional tools and functions will be added here. If you have any questions about any of the tools provided in the `extras/` directory, please get in touch through the [SLiMSuite GitHub community][1].


[1]: https://github.com/slimsuite/SLiMSuite/discussions
[2]: https://github.com/slimsuite/SLiMSuite/issues
[3]: http://slimsuite.blogspot.com/
[4]: https://pubmed.ncbi.nlm.nih.gov/20486118/
[5]: https://pubmed.ncbi.nlm.nih.gov/16159912/
[6]: http://www.slimsuite.unsw.edu.au/software.php
[7]: http://rest.slimsuite.unsw.edu.au/aphid
[8]: http://rest.slimsuite.unsw.edu.au/badasp
[9]: https://pubmed.ncbi.nlm.nih.gov/20924652/
[10]: http://rest.slimsuite.unsw.edu.au/budapest
[11]: https://github.com/slimsuite/buscomp
[12]: https://slimsuite.github.io/buscomp/
[13]: https://github.com/slimsuite/SLiMSuite/blob/master/docs/manuals/BADASP%20Manual.pdf
[14]: http://bit.ly/CompariMotifManual
[15]: http://rest.slimsuite.unsw.edu.au/comparimotif
[16]: https://pubmed.ncbi.nlm.nih.gov/18375965/
[17]: https://slimsuite.github.io/diploidocus/
[18]: https://github.com/slimsuite/diploidocus
[19]: http://rest.slimsuite.unsw.edu.au/gablam
[20]: https://github.com/slimsuite/pafscaff
[21]: https://slimsuite.github.io/pafscaff/
[22]: https://www.ncbi.nlm.nih.gov/pubmed/32236524/
[23]: http://f1000research.com/posters/4-1022
[24]: http://rest.slimsuite.unsw.edu.au/pagsat
[25]: https://github.com/slimsuite/saaga
[26]: https://slimsuite.github.io/saaga/
[27]: https://github.com/slimsuite/SLiMSuite/wiki/SAMPhaser
[28]: https://pubmed.ncbi.nlm.nih.gov/32696352/
[29]: http://bit.ly/SFManual
[30]: https://slimsuite.github.io/synbad/
