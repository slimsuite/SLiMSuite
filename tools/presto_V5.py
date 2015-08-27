#!/usr/local/bin/python

# PRESTO - Peptide Regular Expression Search Tool
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      PRESTO
Description:  Protein Regular Expression Search Tool
Version:      5.0
Last Edit:    22/01/07
Note:         This program has been superceded in most functions by SLiMSearch.
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
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

    PRESTO recognises "n of m" motif elements in the form <X:n:m>, where X is one or more amino acids that must occur n+
    times across which m positions. E.g. <IL:3:5> must have 3+ Is and/or Ls in a 5aa stretch.

    Main output for PRESTO is a delimited file of motif/peptide occurrences but the motifaln=T and proteinaln=T also allow
    output of alignments of motifs and their occurrences. PRESTO has an additional motinfo=FILE output, which produces a
    summary table of the input motifs, inlcuding Expected values if searchdb given and information content if motifIC=T.
    Hit proteins can also be output in fasta format (fasout=T) or UniProt format with occurrences as features (uniprot=T).
    
Release Notes:
    Expectation scores have now been modified since PRESTO Version 1.x. In addition to the expectation score for the no.
    of occurrences of a given motif (given the number of mismatches) in the entire dataset ("EXPECT"), there is now an
    estimation of the probability of the observed number of occurrences, derived from a Poisson distribution, which is
    output in the log file ("#PROB"). Further more, these values are now also calculated per sequence individually
    ("SEQ_EXP" and "SEQ_PROB").

    Note on MS-MS mode: The old Perl version of Presto had a handy MS-MS mode for searching peptides sequenced from tandem
    mass-spec data. (In this mode [msms=T], amino acids of equal mass (Leu-Ile [LI], Gln-Lys [QK], MetO-Phe [MF]) are
    automatically placed as possible variants and additional output columns give information of predicted tryptic fragment
    masses etc.) Implementation of MS-MS mode has been started in this version but discontinued due to lack of demand. As a
    result, extra tryptic fragment data is not produced. If you would like to use it, contact me at richard.edwards@ucd.ie
    and I will finish implementing it.

    Note for compare=T mode: This is still fully functional but main documentation has been moved to comparimotif.py.

    !!!NEW!!! for version 3.7, PRESTO has an additional domfilter=FILE option. This is quite crude and will read in domains
    to be filtered from the FILE given. This file MUST be tab-delimited and must have at least three columns, with headers
    'Name','Start' and 'Stop', where Name matches the short name of the Hit and 'Start' and 'End' are the positions of the
    domain 1-N. This will output two additional columns, plus a further two if iupred=T:
    * DOM_MASK = Gives the motif a score of the length of the domain if it would be masked out by masking domains or 0 if not
    * DOM_PROP = Gives the proportion of motif positions in a domain
    * DOM_DIS = Gives the motif the mean disorder score for the *domain* if in the domain, else 1.0 if not
    * DOM_COMB = Gives positions in the domain the mean disorder score for the domain, else they keep their own scores

    !!!NEW!!! for version 4.0, PRESTO has a Peptide design mode (peptides=T), using winsize=X to set size of peptides around
    occurrences. This will output peptide sequences into a fasta file and additional columns to the main PRESTO output file:
    * PEP_SEQ = Sequence of peptide
    * PEP_DESIGN = Peptide design comments. "OK" if all looking good, else warnings bad AA combos (DP, DC, DG, NG, NS or PP)

Development Notes: (To be assimilated with release notes etc. when version is fully functional.)
    Main output is now determined by outfile=X and/or basefile=X, which will set the self.info['Basefile'] attribute,
    using standard rje module commands. If it is not set (i.e. is '' or 'None'), it will be generated using the motif and
    searchdb files as with the old PRESTO. Main search output will use this file leader and add the appropriate extension based
    on the output type and delimiter:
    * resfile = Main PRESTO search = *.presto.tdt
    * motifaln = Produce fasta files of local motif alignments *.motifaln.fas
    * uniprot = Output of hits as a uniprot format file = *.uniprot.presto 
    * motinfo = Motif summary table = *.motinfo.tdt
    * ftout = Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    * peptides = Peptides designed around motifs = *.peptides.fas

    Other special output will generate their names using protein and/or motif names using the root PATH of basefile
    (e.g. the PATH will be stripped and ProteinAln/ or HitFas/ directories made for output):
    * proteinaln=T/F  : Search for alignments of proteins containing motifs and produce new file containing motifs in [False]
    * fasout=T/F      : Whether to output hit sequences as a fasta format file motif.fas [False]
    
    Reformatting and ouputting motifs require a file name to be given:
    * motifout=FILE   : Filename for output of reformatted (and filtered?) motifs in PRESTO format [None]
    
    

PRESTO Commands:
    ## Basic Input Parameters ##
    motifs=FILE     : File of input motifs/peptides [None]
                      Single line per motif format = 'Name Sequence #Comments' (Comments are optional and ignored)
                      Alternative formats include fasta, SLiMDisc output and raw motif lists.
    minpep=X        : Min length of motif/peptide X aa [2]
    minfix=X        : Min number of fixed positions for a motif to contain [0]
    minic=X         : Min information content for a motif (1 fixed position = 1.0) [2.0]
    trimx=T/F       : Trims Xs from the ends of a motif [False]
    nrmotif=T/F     : Whether to remove redundancy in input motifs [False]
    searchdb=FILE   : Protein Fasta file to search (or second motif file to compare) [None]
    xpad=X          : Adds X additional Xs to the flanks of the motif (after trimx if trimx=T) [0]
    xpaddb=X        : Adds X additional Xs to the flanks of the search database sequences (will mess up alignments) [0]
    minimotif=T/F   : Input file is in minimotif format and will be reformatted (PRESTO File format only) [False]
    goodmotif=LIST  : List of text to match in Motif names to keep (can have wildcards) []

    ## Basic Output Parameters ##
    outfile=X       : Base name of results files, e.g. X.presto.tdt. [motifsFILE-searchdbFILE.presto.tdt]
    expect=T/F      : Whether to give crude expect values based on AA frequencies [True]
    nohits=T/F      : Save list of sequence IDs without motif hits to *.nohits.txt. [False]
    useres=T/F      : Whether to append existing results to *.presto.txt and *.nohits.txt (continuing afer last sequence)
                      and/or use existing results in to search for conservation in alignments if usealn=T. [False]
    mysql=T/F       : Output results in mySQL format - lower case headers and no spaces [False]
    hitname=X       : Format for Hit Name: full/short/accnum [short]
    fasout=T/F      : Whether to output hit sequences as a fasta format file motif.fas [False]
    datout=T/F      : Whether to output hits as a uniprot format file *.uniprot.presto [False]
    motinfo=T/F     : Whether to output motif summary table *.motinfo.tdt [None]
    motifout=FILE   : Filename for output of reformatted (and filtered?) motifs in PRESTO format [None]

    ## Advanced Output Options ##
    winsa=X         : Number of aa to extend Surface Accessibility calculation either side of motif [0]
    winhyd=X        : Number of aa to extend Eisenberg Hydrophobicity calculation either side of motif [0]
    windis=X        : Extend disorder statistic X aa either side of motif (use flanks *only* if negative) [0]
    winchg=X        : Extend charge calculations (if any) to X aa either side of motif [0]
    winsize=X       : Sets all of the above window sizes (use flanks *only* if negative) [0]
    slimchg=T/F     : Calculate Asolute, Net and Balance charge statistics (above) for occurrences [False]
    iupred=T/F      : Run IUPred disorder prediction [False]
    foldindex=T/F   : Run FoldIndex disorder prediction [False]
    iucut=X         : Cut-off for IUPred results (0.0 will report mean IUPred score) [0.0]
    iumethod=X      : IUPred method to use (long/short) [short]
    iupath=PATH     : The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    domfilter=FILE  : Use the DomFilter options, reading domains from FILE [None]
    runid=X         : Adds an additional Run_ID column identifying the run (for multiple appended runs [None]
    restrict=LIST   : List of files containing instances (hit,start,end) to output (only) []
    exclude=LIST    : List of files containing instances (hit,start,end) to exclude []
    peptides=T/F    : Peptide design mode, using winsize=X to set size of peptides around motif [False]
    newscore=LIST   : Lists of X:Y, create a new statistic X, where Y is the formula of the score. []

    ## Basic Search Parameters ##
    mismatch=X,Y    : Peptide must be >= Y aa for X mismatches
    ambcut=X        : Cut-off for max number of choices in ambiguous position to be shown as variant [10]
    expcut=X        : The maximum number of expected occurrences allowed to still search with motif [0] (if -ve, per seq)
    alphabet=X,Y,.. : List of letters in alphabet of interest [AAs]
    reverse=T/F     : Reverse the motifs - good for generating a test comparison data set [False]
    *** No longer outputs *.rev.txt - use motifout=X instead! ***

    msms=T/F        : Whether searching Tandem Mass Spec peptides [False]
    ranking=T/F     : Whether to rank hits by their rating in MSMS mode [False]
    memsaver=T/F    : Whether to store all results in Objects (False) or clear as search proceeds (True) [True]    
    startfrom=X     : Accession Number / ID to start from. (Enables restart after crash.) [None]

    ## Conservation Parameters ##
    usealn=T/F      : Whether to search for and use alignemnts where present. [False]
    gopher=T/F      : Use GOPHER to generate missing orthologue alignments in alndir - see gopher.py options [False]
    fullforce=T/F   : Force GOPHER to re-run even if alignment exists [False]
    alndir=PATH     : Path to alignments of proteins containing motifs [./] * Use forward slashes (/)
    alnext=X        : File extension of alignment files, accnum.X [aln.fas]
    alngap=T/F      : Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore
                      as putative sequence fragments [False]  (NB. All X regions are ignored as sequence errors.)
    conspec=LIST    : List of species codes for conservation analysis. Can be name of file containing list. [None]
    conscore=X      : Type of conservation score used:  [pos]
                        - abs = absolute conservation of motif using RegExp over matched region
                        - pos = positional conservation: each position treated independently 
                        - prop = conservation of amino acid properties
                        - all = all three methods for comparison purposes
    consamb=T/F     : Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    consinfo=T/F    : Weight positions by information content (does nothing for conscore=abs) [True]
    consweight=X    : Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
                        - 0 gives equal weighting to all. Negative values will upweight distant sequences.
    posmatrix=FILE  : Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) [None]
    aaprop=FILE     : Amino Acid property matrix file. [aaprop.txt]
    consout=T/F     : Outputs an additional result field containing information on the conservation score used [False]

    ## Additional Output for Extracted Motifs ##
    motific=T/F     : Output Information Content for motifs [False]
    motifaln=T/F    : Produce fasta files of local motif alignments [False]  
    proteinaln=T/F  : Search for alignments of proteins containing motifs and produce new file containing motifs [False]
    protalndir=PATH : Directory name for output of protein aligments [ProteinAln/]
    flanksize=X     : Size of sequence flanks for motifs [30]
    xdivide=X       : Size of dividing Xs between motifs [10]
    ftout=T/F       : Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    unipaths=LIST   : List of additional paths containing uniprot.index files from which to look for and extract features ['']
    statfilter=LIST : List of stats to filter (*discard* occurrences) on, consisting of X*Y where:
                      - X is an output stat (the column header),
                      - * is an operator in the list >, >=, !=, =, >= ,<    !!! Remember to enclose in "quotes" for <> !!!
                      - Y is a value that X must have, assessed using *.
                      This filtering is crude and may behave strangely if X is not a numerical stat!

    ## Motif Comparison Parameters ##
    compare=T/F     : Compare the motifs from the motifs FILE with the searchdb FILE (or self if None) [False]
    minshare=X      : Min. number of non-wildcard positions for motifs to share [2]
    matchfix=X      : If >0 must exactly match *all* fixed positions in the motifs from:  [0]
                        - 1: input (motifs=FILE) motifs
                        - 2: searchdb motifs
                        - 3: *both* input and searchdb motifs
    matchic=T/F     : Use (and output) information content of matched regions to asses motif matches [True]
    motdesc=X       : Sets which motifs have description outputs (0-3 as matchfix option) [3]
    outstyle=X      : Sets the output style for the resfile [normal]
                        - normal = all standard stats are output
                        - multi = designed for multiple appended runs. File names are also output
                        - single = designed for searches of a single motif vs a database. Only motif2 stats are output
                        - normalsplit/multisplit = as normal/multi but stats are grouped by motif rather than by type

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_aaprop, rje_disorder, rje_motif_V3, rje_motif_cons, rje_scoring, rje_seq, rje_sequence,
    rje_blast, rje_pam, rje_uniprot
Other modules needed: rje_dismatrix, 
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_motif_V3, rje_motif_stats, rje_scoring, rje_seq, rje_sequence
import rje_motiflist 
#########################################################################################################################
### History
# 4.2 - See presto_V4.2 for history upto and including version 4.2.
# 5.0 - Major reworking with new module & class structure to unify better with SLiMPickings 3.0 (& CompariMotif 2.0?)
#########################################################################################################################
### Major Functionality to Add
# [ ] : Rework structure in line with SLiMPickings 3.0
# [ ] : Possible reintroduction of MSMS modes statistics (Add SNT calculations?)
# [ ] : Add proper use of existing results, i.e. read into MotifOcc objects (Currently just makes append=T (+nohits))
# [ ] : Expand use of restrict and exclude to be more flexible w.r.t. input file format.
# [ ] : Change Expectation Filter to use statfilter methods? (Or just fix.)
# [ ] : Go through and identify necessary variables and reduce accordingly
# [ ] : Can UniProt output be improved to be more like SLiMPickings, i.e. embedded in real UniProt entries? (MotifList?)
# [ ] : Check <KR:3:5> vs <[KR]:3:5>
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit) = ('PRESTO', '5.0', 'January 2007')  
    description = 'Peptide Regular Expression Search Tool'
    author = 'Dr Richard J. Edwards.'
    return rje.Info(program,version,last_edit,description,author,time.time())
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program,info.version,time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,-1,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
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
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### Constants                                                                                                           #
#########################################################################################################################
outfile = {'Search':'.presto.tdt','DatOut':'.uniprot.presto','MotInfo':'.motinfo.tdt','FTOut':'.features.tdt',
           'MotifAln':'.motifaln.fas','Peptides':'.peptides.fas'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Presto Class:                                                                                           #
#########################################################################################################################
class Presto(rje.RJE_Object):     
    '''
    PRESTO Class. Author: Rich Edwards (2005).

    Handles main thread for PRESTO motif searches.    

    Info:str
    - Name = Name of input motif/peptide file
    - SearchDB = Name of protein fasta file to search for motifs
    - AlnDir = Path to alignment files
    - AlnExt = File extensions of alignments: AccNum.X
    - ResFile = Base for results files *.presto.txt and *.nohits.txt
    - StartFrom = Accession Number / ID to start from. (Enables restart after crash.) 
    - UniPaths = List of paths containing uniprot.index files from which to look for and extract features ['']
    - ConScore = Type of conservation score used:  [abs]
        - abs = absolute conservation of motif: reports percentage of homologues in which conserved
        - prop = conservation of amino acid properties
        - pam = calculated a weighted conservation based on PAM distance between sequence
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix)
    - OutStyle = Sets the output style for the resfile [normal]
        - normal = all standard stats are output
        - multi = designed for multiple appended runs. File names are also output
        - single = designed for searches of a single motif vs a database. Only motif2 stats are output
        - normalsplit/multisplit = as normal/multi but stats are grouped by motif rather than by type
    - ProtAlnDir = Directory name for output of protein aligments [ProteinAln/]
    - HitName = Format for Hit Name: full/short/accnum [short]
    - DomFilter = Use the DomFilter options, reading domains from FILE [None]
    - RunID = Adds an additional Run_ID column identifying the run (for multiple appended runs [None]
    - MotifOut = Filename for output of reformatted (and filtered?) motifs in PRESTO format [None]
    
    Opt:boolean
    - MotInfo = Output motif summary table [False]
    - NRMotif = Whether to remove redundancy in input motifs [False]
    - Expect = Whether to calculate crude 'expected' values based on AA composition.
    - MSMS = Whether to run in MSMS mode
    - Ranking = Whether to rank hits in MSMS mode. [Currently not implemented.]
    - UseAln = Whether to look for conservation in alignments
    - UseRes = Whether to use existing results files
    - NoHits = Whether to save list of sequences that lack hits
    - FasOut = whether to output hit sequences to fasta file
    - DatOut = whether to output hits to uniprot format file
    - Reverse = Reverse the motifs - good for generating a test comparison data set [False]
    - FTOut = Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [True]
    - MotifAln = Produce fasta files of local motif alignments [True]
    - ProteinAln = Search for alignments of proteins containing motifs and produce new file containing motifs [False]
    - Compare = Compare the motifs from the motifs FILE with the searchdb FILE (or self if None) [False]
    - MatchIC = Use (and output) information content of matched regions to asses motif matches [True]
    - MotifIC = Output Information Content for motifs [False]
    - IUPred = Run IUPred disorder prediction [False]
    - FoldIndex = Run FoldIndex disorder prediction [False]
    - ConsInfo = Weight positions by information content [True]
    - AlnGap = Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore as putative sequence fragments [True]
    - Gopher = Use GOPHER to generate missing orthologue alignments in outdir/Gopher - see gopher.py options [False]
    - FullForce = Force GOPHER to re-run even if alignment exists [False]
    - ConsAmb = Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    - ConsOut = Outputs an additional result field containing information on the conservation score used [False]
    - TrimX = Trims Xs from the ends of a motif
    - Searched = Whether PRESTO has been used to search a database with a list of motifs
    - SlimChg = Calculate Asolute, Net and Balance charge statistics (above) for occurrences [False]
    - Peptides = Peptide design mode, using winsize=X to set size of peptides around motif [False]
    - MiniMotif = Input file is in minimotif format and will be reformatted [False]
    - Search = Whether to search database with motifs [True]
    
    Stat:numeric
    - AmbCut = Cut-off for max number of choices in ambiguous position to be shown as variant [10]
        For mismatches, this is the max number of choices for an ambiguity to be replaced with a mismatch wildcard
    - ExpCut = The maximum number of expected occurrences allowed to still search with motif [0]
    - MinPep = Minimum length of motif/peptide (non-X characters)
    - MinFix = Min number of fixed positions for a motif to contain [0]
    - MinIC = Min information content for a motif (1 fixed position = 1.0) [2.0]
    - FlankSize = Size of sequence flanks for motifs in MotifAln [30]
    - XDivide = Size of dividing Xs between motifs [10]
    - MinShare = Min. number of non-wildcards for motifs to share in Compare [2]
    - MatchFix = If >0 must match *all* fixed positions in the motifs from:  [0]
        - 1: input (motifs=FILE) motifs
        - 2: searchdb motifs
        - 3: *both* input and searchdb motifs
    - ConsWeight = Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
    - MotDesc = Sets which motifs have description outputs (0-3 as matchfix option) [3]
    - XPad = Adds X additional Xs to the flanks of the motif (after trimx if trimx=T) [0]
    - XPadDB = Adds X additional Xs to the flanks of the search database sequences [0]

    List:list
    - Alphabet = List of letters in alphabet of interest
    - Headers = Column headers for data output. Stored for ease of results output   #!# Eliminate when dict functioning #!#
    - Motifs = List of rje_motif_V3.Motif objects
    - SeqHits = List of PrestoSeqHit objects (if memsaver=F)
    ## Advanced Filtering/Ranking Options ##
    - StatFilter = List of stats to filter on, consisting of X*Y where:
          - X is an output stat (the column header),
          - * is an operator in the list >, >=, =, <= ,<, !=
          - Y is a value that X must have, assessed using *.
          This filtering is crude and may behave strangely if X is not a numerical stat!
    - Restrict = List of files containing instances (hit,start,end) to output (only) []
    - Exclude = List of files containing instances (hit,start,end) to exclude []
    - GoodMotif = List of text to match in Motif names to keep (can have wildcards) []
    - NewScore = self.dict['NewScore'] keys() in order they were read in

    Dict:dictionary
    - Expect = Dictionary of {Motif:{mm:exp}}
    - MisMatch = Dictionary of mismatches X:Y
    - ConsSpecLists = Dictionary of {BaseName:List} lists of species codes for special conservation analyses
    - OccCount = Dictionary of {Motif:Number of occurrences in dataset}
    - MotifOcc = Dictionary of {Sequence:{Motif:{Pos:Match}}}   (Pos is from 0 to L-1)
    - MotOccVar = Dictionary of {Sequence:{Motif:{Pos:Variant}}}   (Pos is from 0 to L-1)
    - ElementIC = Dictionary of {Position Element:IC}
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) {}
    - StatFilter = Dictionary of stat filters made from self.list['StatFilter']
    - DomFilter = Dictionary of {HitName:list of domains [(start,end)] arranged in length order}
    - Restrict = Dictionary of instances {motif:[(hit,start,end)]} to output (only) []
    - Exclude = Dictionary of instances {motif:[(hit,start,end)]} to exclude []
    - Headers = Dictionary of {Outputfiletype:[Headers]}    # See module output dictionary for types and extensions #
    - NewScore = dictionary of {X:Y} for new statistic X, where Y is the formula of the score. []

    Obj:RJE_Objects
    - AAPropMatrix = rje_aaprop.AAPropMatrix object
    - PAM = rje_pam.PamCtrl object for PAM conservation
    '''
    ### Attributes
    def motifNum(self): return len(self.motifs())
    def motifs(self): return self.obj['MotifList'].motifs()
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object'''
        ### <a> ### Basics 
        self.infolist = ['Name','SearchDB','AlnDir','AlnExt','ResFile','StartFrom','UniPaths','ConScore',
                         'PosMatrix','OutStyle','ProtAlnDir','HitName','DomFilter','RunID','MotifOut']
        self.statlist = ['AmbCut','ExpCut','MinPep','FlankSize','XDivide','XPad','MinShare','MinFix','MatchFix',
                         'ConsWeight','MotDesc','XPad','XPadDB','MinIC']
        self.optlist = ['Expect','Ranking','MSMS','UseAln','UseRes','NoHits','FasOut','DatOut','Reverse','AlnGap',
                        'FTOut','MotifAln','ProteinAln','NRMotif','Compare','MatchIC','MotifIC','IUPred','FoldIndex',
                        'ConsInfo','ConsAmb','Gopher','FullForce','ConsOut','TrimX','Searched','SlimChg','Peptides',
                        'MiniMotif','MotInfo','Search']
        self.listlist = ['Headers','Motifs','Alphabet','SeqHits','StatFilter','Restrict','Exclude','GoodMotif',
                         'NewScore']
        self.dictlist = ['ConsSpecLists','Expect','MisMatch','OccCount','MotifOcc','MotOccVar','ElementIC','PosMatrix',
                         'StatFilter','DomFilter','Restrict','Exclude','NewScore']
        self.objlist = ['PAM','AAPropMatrix']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'MinFix':0,'MinPep':2,'FlankSize':30,'XDivide':10,'MinShare':2,'AmbCut':10,'MotDesc':3,
                      'XPad':0,'XPadDB':0,'MinIC':2.0})
        self.setOpt({'Expect':True,'MatchIC':True,'MemSaver':True,'ConsAmb':True,'ConsInfo':True,'Search':True})
        self.setInfo({'AlnDir':'Gopher/ALN/','AlnExt':'orthaln.fas','ConScore':'pos','UniPaths':'','OutStyle':'normal',
                      'ProtAlnDir':rje.makePath('ProteinAln/'),'HitName':'short'})
        self.list['Alphabet'] = string.split('A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y',',')
        self.obj['MotifList'] = rje_motiflist.MotifList(self.log,self.cmd_list)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)

                ### Basic Input/Search Options ###
                self._cmdRead(cmd,type='info',att='Name',arg='motifs')  

                self._cmdReadList(cmd,'file',['SearchDB'])
                self._cmdReadList(cmd,'int',['XPad','XPadDB','AmbCut','MinPep','MinFix'])
                self._cmdReadList(cmd,'list',['GoodMotif','Alphabet'])
                self._cmdReadList(cmd,'stat',['ExpCut','MinIC'])
                self._cmdReadList(cmd,'opt',['NRMotif','Expect','MiniMotif','Reverse','AlnGap','MSMS','Search'])
                mismatch = rje.matchExp('mismatch=(\d+),(\d+)',cmd)
                if mismatch:
                    self.dict['MisMatch'][mismatch[0]] = string.atoi(mismatch[1])
                self._cmdRead(cmd,type='info',att='StartFrom')

                ## Output ##
                self._cmdReadList(cmd,'info',['ResFile','HitName','RunID'])
                self._cmdReadList(cmd,'file',['DomFilter','MotifOut'])
                self._cmdReadList(cmd,'opt',['NoHits','UseRes','DatOut','Ranking','FasOut','SlimChg','Peptides','MotInfo'])
                self._cmdReadList(cmd,'glist',['Restrict','Exclude'])
                ## Special NewScore ##
                self._cmdRead(cmd,'cdictlist','NewScore')
                
                ## Conservation Stats ##
                self._cmdRead(cmd,type='opt',att='UseAln')
                self._cmdRead(cmd,type='info',att='AlnDir')
                self._cmdRead(cmd,type='info',att='AlnExt')
                self._cmdRead(cmd,type='info',att='ConScore')
                self._cmdRead(cmd,type='opt',att='ConsInfo')
                self._cmdRead(cmd,type='opt',att='ConsAmb')
                self._cmdReadList(cmd,'opt',['Gopher','FullForce'])
                self._cmdRead(cmd,type='info',att='PosMatrix')
                self._cmdRead(cmd,type='opt',att='ConsOut')
                self._cmdRead(cmd,type='stat',att='ConsWeight')

                ## Additional Output ##
                self._cmdRead(cmd,type='opt',att='MotifAln')
                self._cmdRead(cmd,type='int',att='FlankSize')
                self._cmdRead(cmd,type='int',att='XDivide')
                self._cmdRead(cmd,type='opt',att='FTOut')
                self._cmdRead(cmd,type='clist',att='UniPaths')
                self._cmdRead(cmd,type='opt',att='ProteinAln')
                self._cmdRead(cmd,type='path',att='ProtAlnDir')
                self._cmdRead(cmd,type='opt',att='IUPred')
                self._cmdRead(cmd,type='opt',att='FoldIndex')
                self._cmdRead(cmd,type='list',att='StatFilter')

                ## Motif Compare ##
                self._cmdReadList(cmd,'info',['OutStyle'])
                self._cmdReadList(cmd,'int',['MinShare','MatchFix','MotDesc'])
                self._cmdReadList(cmd,'opt',['Compare','MatchIC','MotifIC','TrimX'])

                ### <c> ### Special
                if rje.matchExp('^conspec=(.+)',cmd):
                    conspec = rje.matchExp('^conspec=(.+)',cmd)[0]
                    consfiles = glob.glob(conspec)
                    if len(consfiles) > 0:     # File
                        for cfile in consfiles:
                            self.dict['ConsSpecLists'][rje.baseFile(cfile,True)] = rje.listFromCommand(cfile)
                    else:
                        self.dict['ConsSpecLists']['SPEC'] = rje.listFromCommand(conspec,checkfile=False)
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.opt['UseRes']:  #!# Needs improvement #!#
            self.opt['Append'] = True
        self._cmdMismatch()
#########################################################################################################################
    def _cmdMismatch(self):      ### Sets up the search Mismatchdic from cmd_list
        '''Sets up the search Mismatchdic from cmd_list.'''
        try:
            mmkeys = self.dict['MisMatch'].keys()
            ## Catch dodgy input ##
            for mk in mmkeys[0:]:
                try:
                    mm = string.atoi(mk)
                    if mm < 1:
                        mmkeys.remove(mk)
                except:
                    mmkeys.remove(mk)
            needval = []
            mm = 1
            while len(mmkeys) > 0:
                mmkey = '%d' % mm
                if mmkey in mmkeys:
                    for key in needval:
                        self.dict['MisMatch'][key] = self.dict['MisMatch'][mmkey]
                    needval = []
                    mmkeys.remove(mmkey)
                else:
                    needval.append(mmkey)
                mm += 1
        except:
            self.log.errorLog('Problem with Presto._cmdMisMatch: %s' % self.dict['MisMatch'])
        self.obj['MotifList'].dict['MisMatch'] = self.dict['MisMatch']
#########################################################################################################################
    def run(self):      ### General Run method for Presto
        '''
        General Run method for Presto:
        0. Check relevant files
        1. Load motifs, Reformat and check for redundancy, Make (mismatch) variants
        2. Search Database and output results (or output MotInfo table)
        '''
        try:
            ### Setup PRESTO Run ###
            if not self.setupRun():
                return

            ### Menu ###
            if self.stat['Interactive'] > 0:
                self.runMenu()   # Fully interactive run including individual options #
                #!# Add return when this exists #!#

            ### Load Motifs ### Loads, checks for redundancy, reverses (if desired) and makes variants #
            self.obj['MotifList'].loadMotifs(self.info['Name'])
            if self.stat['ExpCut'] > 0:
                self.opt['Expect'] = True
            if self.opt['Expect']:
                self.obj['MotifList'].seqListExp(filename=self.info['SearchDB'],cutoff=self.stat['ExpCut'])
                self.log.printLog('\r#EXP','Expect value calculations complete.')
            self.setupResults()

            ### Optional Output ###
            if self.info['MotifOut'] and self.info['MotifOut'].lower() != 'none':
                rje.backup(self,self.info['MotifOut'])
                OUT = open(self.info['MotifOut'],'a')
                for motif in self.motifs():
                    OUT.write('%s  %s  # %s\n' % (motif.info['Name'],motif.info['Sequence'],motif.info['Description']))
                OUT.close()   
            
            ### Motif Search ###
            if self.opt['Search']:
                self.searchDB()

            ### Non-Memsaver output ###
            if self.opt['MotifAln'] and not self.opt['MemSaver']:    
                resfile = self.info['Basefile'] + 'motifaln.fas'
                self.obj['MotifList'].motifAlignments(resfile)

            ### Motif Information File. Now contains Hit Information Too ###
            if self.opt['MotInfo']:
                self.motifInfo(expfile=self.info['SearchDB'])

        except SystemExit:
            raise
        except:
            self.log.errorLog('Error in Presto.run()',quitchoice=True) 
#########################################################################################################################
    def runMenu(self):  #!# Set this up using rje_menu and all options #!#
        self.verbose(0,1,'\n\nPRESTO Menu will be added shortly. Please use help option for command options.\n\n',1)
        return
#########################################################################################################################
    def setupRun(self):      ### Sets up input and output files etc.
        '''Sets up input and output files etc.'''
        try:
            ### Setup Input Files: I - Motif File ###
            while not os.path.exists(self.info['Name']):
                self.log.errorLog('Motif file (%s) missing.' % self.info['Name'],printerror=False,quitchoice=True)
                choice = ''
                if self.stat['Interactive'] == 0:
                    self.stat['Interactive'] = 1
                if self.stat['Interactive'] > 0:    #!# When menu is added, change this to only exit if i=0 #!#
                    self.info['Name'] = rje.choice('Enter motif file name (blank to exit): ')
                if self.info['Name'] == '':
                    sys.exit()
            
            ### Setup Input Files: II - SearchDB File ###
            self.opt['Searched'] = False
            while not os.path.exists(self.info['SearchDB']) and self.info['MotifOut'].lower() in ['','none']:
                self.log.errorLog('Search database (%s) missing.' % self.info['SearchDB'],printerror=False,quitchoice=True)
                choice = ''
                if self.stat['Interactive'] >= 0:
                    choice = rje.choice('Enter search database file name (blank to exit): ')
                    if choice:
                        self.info['SearchDB'] == choice
                if choice == '':
                    sys.exit()

            ### Setup Default Result Filenames ###
            if self.info['Basefile'].lower() in ['','none']:
                mname = rje.baseFile(self.info['Name'],strip_path=True,extlist=['txt','motif','motifs','fas'])
                sname = rje.baseFile(self.info['SearchDB'],strip_path=True,extlist=['txt','motif','motifs','fas'])
                self.info['Basefile'] = '%s-%s' % (mname,sname)

            return True
        except:
            self.log.errorLog('Error in run Setup')
            sys.exit()
#########################################################################################################################
    def setupResults(self):  ### Sets up headers for relative output files into self.dict['Headers']
        '''Sets up headers for relative output files into self.dict['Headers'].'''
        try:
            ### Setup Headers ###  (NB. Feature table should be handled by UniProt object in MotifList)
            self.dict['Headers'] = {'Search':self.prestoHeaders(),'MotInfo':self.motInfoHeaders()}
            
            ### Backups and Header Outputs ###
            delimit = rje.getDelimit(self.cmd_list)
            for out in outfile.keys():
                if self.opt[out]:
                    resfile = '%s%s' % (self.info['Basefile'],outfile[out])
                    #X#print out, self.info['Basefile'], outfile[out], resfile
                    rje.backup(self,resfile)
                    if self.dict['Headers'].has_key(out):
                        rje.delimitedFileOutput(self,resfile,self.dict['Headers'][out],delimit=delimit)

            ### Additional Output Files ###
            basedir = os.path.split(self.info['Basefile'])[0]
            resbase = rje.baseFile(self.info['Basefile'],True)
            if basedir:
                basedir = basedir + '/'
            ## ProteinAln ##
            if self.opt['ProteinAln'] and not os.path.exists(rje.makePath('%sProteinAln/' % basedir)):
                os.mkdir(rje.makePath('%sProteinAln/' % basedir))
            self.info['ProtAlnDir'] = rje.makePath('%sProteinAln/' % basedir)
            ## MotifAln ##
            if self.opt['MotifAln'] and not os.path.exists(rje.makePath('%sMotifAln/' % basedir)):
                os.mkdir(rje.makePath('%sMotifAln/' % basedir))
            self.info['MotifAlnDir'] = rje.makePath('%sMotifAln/' % basedir) + resbase + '.'
            for Motif in self.motifs():
                resfile = self.info['MotifAlnDir'] + Motif.info['Name'] + '.motifaln.fas'
                rje.backup(self,resfile)
            ## FasOut ##
            if self.opt['FasOut'] and not os.path.exists(rje.makePath('%sHitFas/' % basedir)):
                os.mkdir(rje.makePath('%sHitFas/' % basedir))
            for motif in self.motifs():
                if self.opt['FasOut']:
                    motif.info['FasOut'] = rje.makePath('%sHitFas/%s.%s.hits.fas' % (basedir,resbase,motif.info['Name']),wholepath=True)
                    rje.backup(self,motif.info['FasOut'])

        except:
            self.log.errorLog('Error in setupResults()')
#########################################################################################################################
    def prestoHeaders(self):    ### Sets up main results file headers from Options
        '''Sets up main results file headers from Options.'''
        ### Basics ###
        header = ['MOTIF','LEN','VARIANT','MATCHSEQ','MATCH_ID','MOTIF_CONS','HOM_NUM','GLOB_ID','LOC_ID','SA','HYD']
        ### MotifList Options ###
        mlist = self.obj['MotifList']
        if mlist.opt['SlimChg']:
            header += ['CHG_ABS','CHG_NET','CHG_BAL']
        for o in ['IUPred','FoldIndex']:
            if mlist.opt[o]:
                header.append(o.upper())
        if mlist.info['DomFilter'] not in ['','None']:
            mlist.setupDomFilter()
            header += ['DOM_MASK','DOM_PROP']
            if mlist.opt['IUPred']:
                header += ['DOM_DIS','DOM_COMB']
        if self.opt['MySQL']:
            header.append('ACC_NUM')
        else:
            header.append('HIT')
        ### Presto (self) options) ###
        header += ['START_POS','END_POS','PROT_LEN']
        if self.opt['MotifIC']:
            header += ['MOTIF_IC']
        if self.opt['Expect']:
            header += ['EXPECT','SEQ_EXP','SEQ_PROB']
        ### Conservation (MotifList) options ###
        for key in rje.sortKeys(mlist.dict['ConsSpecLists']):
            for h in ['CONS','HOM','GLOB_ID','LOC_ID']:
                header.append('%s_%s' % (key,h))
        if mlist.info['ConScore'] == 'all':
            header += ['CONS_ABS','CONS_POS','CONS_PROP']   # MOTIF_CONS takes mean value!
            for key in rje.sortKeys(mlist.dict['ConsSpecLists']):
                for h in ['CONS_ABS','CONS_POS','CONS_PROP']:
                    header.append('%s_%s' % (key,h))
        if self.opt['ConsOut']:
            header.append('CONS_METHOD')
        ### More Presto (self) options ###
        if self.info['RunID'] not in ['','None']:
            header.append('RUN_ID')
        if self.opt['MSMS']:
            header += ['SNT','RATING','ORFMwt','M-ORFMWt','FRAGMWt','N[KR]','C[KR]']
            if self.opt['Ranking']:
                header += ['RANK']
        if self.opt['Peptides']:
            header += ['PEP_SEQ','PEP_DESIGN']
        ## Custom scores ##
        (header,self.list['NewScore'],self.dict['NewScore']) = rje_scoring.setupCustomScores(self,header[0:],self.list['NewScore'],self.dict['NewScore'])
                    
        ### Finish ###
        return header
#########################################################################################################################
    def motInfoHeaders(self):   ### Sets up Motif Info Table headers
        '''Sets up Motif Info Table headers.'''
        headers = ['Motif','Pattern','Description','MaxLength','MinLength','FixLength','FullLength']
        if self.opt['Expect'] and os.path.exists(self.info['SearchDB']):
            headers.append('Expect')
        if self.opt['MotifIC']:
            headers.append('IC')
        if os.path.exists(self.info['SearchDB']) and self.opt['Search']:
            headers += ['OccNum','OccSeq']
        return headers
#########################################################################################################################
    ### <2> ### Database Searching                                                                                      #
#########################################################################################################################
    def searchDB(self):      ### Main method for searching database with motifs. 
        '''Main method for searching database with motifs.'''
        try:
            ### Setup Filtering ###
            self.dict['StatFilter'] = rje_scoring.setupStatFilter(self,self.dict['Headers']['Search'],self.list['StatFilter'])
            self.setupRestrictExclude() #!# Needs improvement #!#

            ### Sequences without hits ###
            nohitsfile = '%s.nohits.txt' % self.info['Basefile']
            if self.opt['NoHits']:
                rje.backup(self,nohitsfile)
            nohits = []     ### List of sequences without any hits
            if self.opt['UseRes']:
                nhlines = self.loadFromFile(nohitsfile)
                nohits = string.split(rje.chomp(string.join(nhlines)))

            ### Search database ### 
            startfrom = self.info['StartFrom']
            seqnum = rje_seq.SeqCount(self,self.info['SearchDB'])
            searchdb = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=F','seqnr=F','autofilter=F'])
            SEARCHDB = open(self.info['SearchDB'],'r')
            seqx = 0.0    # Counter for Log output
            (nextseq, lastline) = searchdb.nextFasSeq(fileobject=SEARCHDB,lastline='')
            while nextseq:  # Sequence returned: search!
                seqx += 100.0 
                stxt = 'Searching %s sequences: %.1f%%.' % (rje.integerString(seqnum),(seqx/seqnum))
                self.log.printLog('\r#PRESTO',stxt,log=False,newline=False)
                searchdb.seq = [nextseq]
                ## Skip if appropriate ##
                if startfrom != 'None':
                    if nextseq.info['Name'].find('%s ' % self.info['StartFrom']) >= 0:
                        startfrom = 'None'
                    else:
                        (nextseq, lastline) = searchdb.nextFasSeq(fileobject=SEARCHDB,lastline=lastline)
                        continue
                if nextseq.shortName() not in nohits:   # Previously looked at and no results
                    ## Search Sequence ##
                    seqhit = self.searchSeq(seq=nextseq,logtext=stxt)
                (nextseq, lastline) = searchdb.nextFasSeq(fileobject=SEARCHDB,lastline=lastline)

            ### Finish ###
            SEARCHDB.close()
            self.opt['Searched'] = True
            self.log.printLog('\n#PRESTO','%s sequences searched for %s motifs.' % (rje.integerString(seqnum),rje.integerString(self.motifNum())))

        except SystemExit:
            sys.exit()
        except:
            self.log.errorLog('Error in Presto.searchDB()')
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def searchSeq(self,seq,logtext=''):   ### Searches a single sequence with motifs.
        '''
        Searches a single sequence with motifs. Based on old newSeqHit method and PrestoSeqHit object. First, a
        list of MotifOcc objects is built, then this is given to rje_motif_stats.
        >> seq:Sequence object to search
        >> logtext:str = Leading text of progress log output
        '''
        try:
            ### Setup ###
            seq.deGap()
            ### XPad ###
            xpaddb = self.obj['MotifList'].stat['XPadDB']
            sequence = 'X' * xpaddb + seq.info['Sequence'] + 'X' * xpaddb
            mx = 0.0    # Motif progress counter
            ## AAFreq ##
            aafreq = rje_sequence.aaFreq(sequence,aafreq={},newkeys=True)
            aafreq = rje.dictFreq(aafreq)
            
            ### Search ###
            occlist = []
            for Motif in self.motifs():
                stxt = '%s %.1f%%' % (logtext,(mx/self.motifNum()))
                mx += 100.0
                hitlist = Motif.searchSequence(sequence,stxt)    # {Pos,Variant,Match,ID,MisMatch}
                seqexp = Motif.expectation(aafreq,aanum=len(sequence),seqnum=1)
                for hitdict in hitlist:
                    pos = hitdict['Pos'] - xpaddb
                    if pos < 1:
                        pos -= 1
                    occdata = {'Stat':{'Pos':pos,'ID':hitdict['ID'],'MisMatch':hitdict['MisMatch']},
                               'Info':{'Variant':hitdict['Variant'],'Match':hitdict['Match'],'SearchVar':hitdict['SearchVar']},
                               'Data':{'SEQ_EXP':rje_motif_V3.expectString(seqexp[hitdict['MisMatch']]),
                                       'SEQ_PROB':rje_motif_V3.expectString(rje_motif_V3.occProb(observed=len(hitlist),expected=seqexp[hitdict['MisMatch']]))}}
                    occlist.append(self.obj['MotifList'].addOcc(seq,Motif,occdata))
            stxt = '%s 100%% [n mismatch]:' % logtext

            ### Calculate MotifOcc Stats ###
            rje_motif_stats.occStats(self.obj['MotifList'],occlist,logtext=('\r#PRESTO',stxt))    
            occlist = self.processSeqOccs(occlist)   ### Filtering occurs here ###
        
            ### Finish ###        
            if not occlist and self.opt['NoHits']:
                open('%s.nohits.txt' % self.info['Basefile'],'a').write('%s\n' % seq.shortName())

            #!# What's this? #!#
            #!# self.uniProtHits()

        except:
            self.log.errorLog('Fatal error searching sequence "%s".' % seq.shortName())
            raise          
#########################################################################################################################
    def processSeqOccs(self,occlist):    ### Processes Occurrences after search - output to results file(s)
        '''
        Processes Occurrences after search - output to results file(s).
        >> occlist:list of MotifOcc objects - must have same obj['Seq']
        << returns filtered occlist
        '''
        try:
            ### Setup ###
            if not occlist:
                return []
            Seq = occlist[0].obj['Seq']
            seqname = Seq.shortName()
            for key in Seq.info.keys():
                if key.lower() == self.info['HitName'].lower():
                    seqname = Seq.info[key]
            delimit = rje.getDelimit(self.cmd_list)

            ### Occurrence results ###
            mymotifs = []
            results = {}
            for hit in occlist[0:]:
                motif = hit.obj['Motif']
                if motif not in mymotifs:
                    mymotifs.append(motif)
                start_pos = hit.stat['Pos'] 
                end_pos = hit.stat['Pos'] + len(hit.info['Match']) - 1
                #X#self.deBug(hit.stat)
                ## Results ##
                results[hit] = {'MOTIF':motif.info['Name'],
                                'LEN':'%d' % len(hit.info['Variant']),
                                'VARIANT':hit.info['Variant'],
                                'MATCHSEQ':hit.info['Match'],
                                'MATCH_ID':'%.2f' % hit.stat['ID'],
                                'MOTIF_CONS':'%.3f' % hit.getStat('ALL_CONS'),
                                'HOM_NUM':'%d' % hit.getStat('ALL_HOM'),
                                'GLOB_ID':'%.3f' % hit.getStat('ALL_GLOB_ID'),
                                'LOC_ID':'%.3f' % hit.getStat('ALL_LOC_ID'),
                                'SA':'%.3f' % hit.stat['SA'],
                                'HYD':'%.3f' % hit.stat['Hyd'],
                                'ACC_NUM':seqname, 'HIT':seqname,
                                'START_POS':'%d' % start_pos,
                                'END_POS':'%d' % end_pos,
                                'PROT_LEN':'%d' % Seq.aaLen(),
                                'EXPECT':motif.dict['Expect'][self.info['SearchDB']],
                                'SEQ_EXP':hit.getData('SEQ_EXP'),
                                'SEQ_PROB':hit.getData('SEQ_PROB'),
                                'MOTIF_IC':'%.3f' % motif.stat['IC'],
                                'RUN_ID':self.info['RunID']}
                if self.opt['SlimChg']:
                    results[hit]['CHG_ABS'] = '%d' % hit.stat['CHG_ABS']
                    results[hit]['CHG_NET'] = '%d' % hit.stat['CHG_NET']
                    results[hit]['CHG_BAL'] = '%d' % hit.stat['CHG_BAL']
                for o in ['IUPred','FoldIndex']:
                    if self.obj['MotifList'].opt[o]:
                        results[hit][o.upper()] = '%.2f' % hit.stat[o]
                for dom in ['MASK','PROP','DIS','COMB']:
                    results[hit]['DOM_%s' % dom] = '%.2f' % hit.stat[dom]
                if self.obj['MotifList'].opt['UseAln']:
                    for key in rje.sortKeys(self.obj['MotifList'].dict['ConsSpecLists']):
                        for h in ['CONS','HOM','GLOB_ID','LOC_ID']:
                            results[hit]['%s_%s' % (key,h)] = '%.3f' % hit.stat['%s_%s' % (key,h)]
                    if self.obj['MotifList'].info['ConScore'] == 'all':
                        try:
                            for method in ['ABS','POS','PROP']:
                                results[hit]['CONS_%s' % method] = hit.stat['ALL_CONS_%s' % method]
                                for key in self.obj['MotifList'].dict['ConsSpecLists'].keys():
                                    results[hit]['%s_CONS_%s' % (key,method)] = hit.stat['%s_CONS_%s' % (key,method)]
                        except:
                            self.verbose(-1,-1,hit.stat,1)
                    if self.opt['ConsOut']:
                        a = '%s' % self.obj['MotifList'].opt['ConsAmb']
                        i = '%s' % self.obj['MotifList'].opt['ConsInfo']
                        results[hit]['CONS_METHOD'] = '%s_w%d_a%s_i%s' % (self.obj['MotifList'].info['ConScore'],self.obj['MotifList'].stat['ConsWeight'],a[0],i[0])
                if self.opt['Peptides']:
                    results[hit]['PEP_SEQ'] = hit.info['PepSeq']
                    results[hit]['PEP_DESIGN'] = rje_sequence.peptideDetails(hit.info['PepSeq'],self)
                ## Custom scores ##
                for new in self.list['NewScore']:
                    results[hit][new] = rje.formula(self,formula=self.dict['NewScore'][new],data=results[hit])

                #!#if self.opt['MSMS']:
                #!#   header += ['SNT','RATING','ORFMwt','M-ORFMWt','FRAGMWt','N[KR]','C[KR]']
                #!#    if self.opt['Ranking']:
                #!#        header += ['RANK']

            ### Filtering & Output ###
            results = rje_scoring.statFilter(self,results,self.dict['StatFilter'])  #!# Add restrict/exclude #!#
            for Occ in occlist[0:]:
                if results.has_key(Occ):
                    rje.delimitedFileOutput(self,'%s%s' % (self.info['Basefile'],outfile['Search']),self.dict['Headers']['Search'],delimit=delimit,datadict=results[Occ])
                    ## Peptide fasta output ##  here #!# What abiut UniProtHits?
                    if self.opt['Peptides']:
                        pepname = '%s-%s-%s %daa (%s)' % (results[Occ]['MOTIF'],results[Occ]['HIT'],results[Occ]['START_POS'],len(results[Occ]['PEP_SEQ']),results[Occ]['PEP_DESIGN'])
                        open('%s%s' % (self.info['Basefile'],outfile['Peptides']),'a').write('>%s\n%s\n' % (pepname,results[Occ]['PEP_SEQ']))
                else:
                    occlist.remove(Occ)
                    if Occ in self.obj['MotifList'].list['MotifOcc']:
                        self.obj['MotifList'].list['MotifOcc'].remove(Occ)

            ## Fasta Files ##
            if self.opt['ProteinAln'] and self.info['HitName'] == 'short':  #!# Gopher should have been run before, if appropriate #!#
                self.obj['MotifList'].singleProteinAlignment(Seq,occlist,alndir=self.info['ProtAlnDir'],hitname='short',gopher=False)
            elif self.opt['ProteinAln']:
                self.obj['MotifList'].singleProteinAlignment(Seq,occlist,alndir='',hitname='AccNum',gopher=False)
            
            ## Fasta per motif ##
            mymotifs = []
            for hit in occlist[0:]:
                motif = hit.obj['Motif']
                if motif not in mymotifs:
                    mymotifs.append(motif)
            if self.opt['FasOut']:
                for motif in mymotifs:
                    open(motif.info['FasOut'],'a').write('>%s\n%s\n' % (Seq.info['Name'],Seq.info['Sequence']))
            if self.opt['MotifAln'] and self.opt['MemSaver']:    
                for motif in mymotifs:
                    resfile = self.info['MotifAlnDir'] + motif.info['Name'] + '.motifaln.fas'
                    if os.path.exists(resfile):
                        self.obj['MotifList'].motifAlnLong(motif,{Seq:occlist},append=True,memsaver=False,resfile=resfile)
                    else:
                        self.obj['MotifList'].motifAlnLong(motif,{Seq:occlist},append=False,memsaver=False,resfile=resfile)


            ## MemSaver ##
            if self.opt['MemSaver']:
                for Occ in occlist:
                    Occ.memSaver()

            return occlist                    

        except:
            self.log.errorLog('Error in processSeqOccs()')
            return []
#########################################################################################################################


    def setupRestrictExclude(self):    ### Sets up Dictionaries of Annotated Motif occurrences {Motif:[Hit,Start,Stop]}
        '''
        Sets up Dictionary of Annotated Motif occurrences {Motif:[Hit,Start,Stop]}. These files must have a single
        header row and four columns containing: Motif, Hit, Start, End. Positions are given 1-N.
        '''
        try:
            ### Setup ###
            for type in ['Restrict','Exclude']:
                self.dict[type] = {}
                if not self.list[type]:
                    continue
                self.log.printLog('#%s' % type.upper(),'Setting up "%s" list of Annotated Motif Instances' % type,log=False)
                ### Read ###
                rx = 0
                for file in self.list[type]:
                    filetype = None
                    for line in self.loadFromFile(file,chomplines=True):
                        if not filetype:    # Check file type
                            if string.split(line,',')[:4] == ['Motif','Len','Hit','Start_Pos']:     # SlimPickings OccRes
                                filetype = 'OccResCSV'
                            else:
                                filetype = 'Simple'
                            continue
                        if filetype == 'OccResCSV':
                            [hit,start,end,matchseq] = string.split(line,',')[2:6]
                            motif = '=%s' % matchseq
                        else:
                            [motif,hit,start,end] = string.split(line)[:4]
                        if not self.dict[type].has_key(motif):
                            self.dict[type][motif] = []
                        self.dict[type][motif].append([hit,start,end])
                        rx += 1
                ### Summary ###
                self.log.printLog('#%s' % type.upper(),'"%s" list of %s Annotated Motif Instances read for %d motifs' % (type,rje.integerString(rx),len(self.dict[type])))
        except:
            self.log.errorLog('setupRestrictExclude fell over',quitchoice=True)
#########################################################################################################################
    ### <4> ### Additional Outputs                                                                                      #
#########################################################################################################################
## These outputs were originally designed for slim_pickings but have been adapted here for PRESTO results:
## - MotifAln produces a single file for the dataset of alignments of each motif occurrence [Memsaver=F only]
## - ProteinAln produces a single file per sequence of all the motifs positioned on the sequence alignmnet
## - FTOut produces a single file for the dataset of extracted UniProt features with the motif positions added
## - MotifInfo produces a summary table of motif information
#########################################################################################################################
### End of SECTION II: Presto                                                                                           #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
    def uniProtHits(self):  ### Processes PrestoSeqHit after search - output to uniprot-format results file
        '''
        Processes PrestoSeqHit after search - output to uniprot-format results file.
        #!# This has not been altered since Presto Version 1.8. #!#
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            RESFILE = open('%s.uniprot.presto' % self.info['ResFile'],'a')
            seq = self.obj['Sequence']

            ### <1> ### SeqDetails
            _stage = '<1> SeqDetails'
            RESFILE.write('ID   %s     STANDARD;     PRT;   %d AA.\n' % (seq.shortName(),seq.aaLen()))
            RESFILE.write('AC   %s;\n' % (seq.info['AccNum']))
            RESFILE.write('DE   %s\n' % (seq.info['Description']))

            ### <2> ### NoHits
            _stage = '<2> NoHits'
            _nohits = []
            for motif in self.list['Motifs']:
                if len(self.dict['Hits'][motif]) == 0:
                    _nohits.append(motif.info['Name'])
            RESFILE.write('CC   -!- PRESTO: No hits for %s.\n' % string.join(_nohits,sep=', '))

            ### <3> ### Hit Features
            _stage = '<3> Hit Features'
            mhit = []
            for motif in self.list['Motifs']:
                for hit in self.dict['Hits'][motif]:
                    (p1,p2) = (hit.stat['Pos'],(hit.stat['Pos'] + len(hit.info['Variant']) - 1))
                    hittxt = 'FT   MOTIF'
                    hittxt += ' ' * (14 - len('%d' % p1)) + '%d' % p1
                    hittxt += ' ' * (7 - len('%d' % p2)) + '%d' % p2
                    hittxt += ' %s (%s)\n' % (motif.info['Name'],motif.info['Sequence'])
                    mhit.append(hittxt)
            mhit.sort()
            for hit in mhit:
                RESFILE.write(hit)
            self.verbose(1,2,' => %d hit features saved to %s.' % (len(mhit),'%s.uniprot.presto' % self.info['ResFile']),1)

            ### <4> ### Sequence/End
            _stage = '<4> Sequence/End'
            RESFILE.write('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE05;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen(),rje_seq.MWt(seq.info['Sequence'])))
            uniseq = seq.info['Sequence'][0:]
            while len(uniseq) > 0:
                RESFILE.write('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                uniseq = uniseq[60:]
            RESFILE.write('//\n')
            RESFILE.close()
        except:
            self.log.errorLog('Error in uniProtHits(%s):' % _stage)
#########################################################################################################################
## End of SECTION III: PrestoSeqHit                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:        
        presto = Presto(log=mainlog,cmd_list=cmd_list)
        presto.run()

    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
