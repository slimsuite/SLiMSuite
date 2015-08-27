#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       slimpickings
Description:  SLiMDisc results compilation and extraction
Version:      3.0
Last Edit:    18/01/07
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This is a basic results compiler for multiple SLiMDisc motif discovery datasets. There are currently the following
    functional elements to the module:

    1. Basic compilation of results from multiple datasets into a single file. This will search through the current
        directory and any subdirectories (unless subdir=F) and pull out results into a single comma-separated file
        (slimdisc_results.csv or outfile=FILE). With the basic run, the following statistics are output:
            ['Dataset','SeqNum','TotalAA','Rank','Score','Pattern','Occ','IC','Norm','Sim']
        This file can then be imported into other applications for analysis. (E.g. rje_mysql.py can be run on the
        file to construct a BUILD statement for MySQL, or StatTranfer can convert the file for STATA analysis etc.)
        !!! NB: If multiple datasets (e.g. in subdirectories) have the same name, slim_pickings will become confused and
        may generate erroneous data later. Please ensure that all datasets are uniquely named. !!!

    2. Additional optional stats based on the motifs sequences themselves to help rank and filter interesting results.
        These are:
            - AbsChg = Number of charged positions [KRDE]
            - NetChg = Net charge of motif [KR] - [DE]
            - BalChg = Balance of charge = Net charge in first half motif - Net charge in second half
            - AILMV = Whether all positions in the motif are A,I,L,M or V.
            - Aromatic = Count of F+Y+W
            - Phos = Potential phosphorylated residues X (none) or [S][T][Y], whichever are present

    3. Calculation of additional statistics from the input sequences, using PRESTO. These are:
            - Mean IUPred/FoldIndex Protein Disorder around the motif occurrences (including extended window either side)
            - Mean Surface Accessibility around the motif occurrences (including an extended window either side)
            - Mean Eisenberg Hydrophobicity around the motif occurrences (including an extended window either side)
            - SLiM conservation across orthologous proteins. (This calculation needs improving.)
        The mean for all occurrences of a motif will be output. In addition, percentile steps can be used to assess
        motifs according to selected threshold criteria (in another package). This will return the threshold at a
        given percentile, e.g. SA_pc75=2.0 would mean that 75% of occurrences have a mean Surface Accessibility value
        of 2.0 or greater. (For hydrophobicity, Hyd_pc50=0.3 would be 50% of occurrences have mean Hydrophobicity of
        0.3 or *less*. This is because low hydrophobicity is good for a (non-structural) functional motif.)

    4. Collation and extraction of key data for specific results. These may be by any combination of protein, motif
        or dataset. If a list of datasets is not given, then all datasets will be considered. (Likewise proteins and
        motifs.) To be very specific, all three lists may be specified (slimlist=LIST protlist=LIST datalist=LIST).
        Information is pulled out in a two-step process:
            (1) The slimdisc.*.index files are consulted for the appropriate list of datasets. If missing, these will
                be regenerated. (slimdisc.motif.index and slimdisc.protein.index both point to dataset names.
                slimdisc.dataset.index points these names to the full path of the results.) Only datasets returned by
                all appropriate lists will be analysed for data extraction.
            (2) The appropriate data on the motifs will be extracted into a directory as determined by outdir=PATH.
                Depending on the options selected, the following (by default all) data is returned:
                - *.motifaln.fas = customised fasta file with motifs aligned in different sequences, ready for dotplots
                    and manual inspection for homology not detected by BLAST.
                - *.dat = UniProt DAT file for as many parents as possible.
        These files will be saved in the directory set by outdir=PATH.

    5. Re-ranking of results. rerank=X will now re-rank the results for each dataset according to the statistic set by
    rankstat=X, and output the top X results only. By default, this is the "R-score" = ic * norm * occ / exp.  The output
    "Rank" will be replaced with the new rank and a new column "OldRank" added to the ouput. zscore=T/F turns on and off
    a simple Z-score calculation based on the slimranks read in. Version 2.5 added a new option for a crude length
    correction of the RScore, dividing by 20 to the power of the motif IC (as calculated by SLiM Pickings on a scale of
    1.0 per fixed position). This is controlled by the lencorrect=T/F option. By default this is False (for backwards
    compatibility) but with future versions this may become the default as it is assumed (by me) that it will improve
    performance. However, there is currently no justification for this, so use with caution!

    6. Filtering of results using the statfilter=LIST option, allowing results to be filtered according to a set of
    rules: LIST should be (a file containing) a comma-separated list of stats to filter on, consisting of X*Y where X is
    an output stat (the column header); * is an operator in the list >, >=, =, =< ,< ; and Y is a value that X must have,
    assessed using *. This filtering is crude and may behave strangely if X is not a numerical stat (although Python does
    seem to assess these alphabetically, so it may be OK)! This filtering is performed before the reranking of the motifs
    if rerank=X is used. This can make run times quite long as many more  motifs need stats calculations. (If rerank=X is
    used without statfilter, re-ranking is done earlier to save time.) See the manual for details.

    7. !!!NEW!!! with version 3.0, customised scores can be created using the newscore=LIST option, where LIST is in the
    form X:Y,X:Y, where in turn X is the name for the new score (a column with this name will be produced) and Y is the
    formula of the score. This formula may contain any output column names, numbers and the operators +-*/^ (^ is "to
    the power of"), using brackets to set the order of calculation. Without brackets, a strict left to right hierarchy is
    observed. e.g. newscore=Eg:3+2*6 will generate a column called "Eg" containing the value 30.0. Custom scores can
    feature previously defined custom scores in the command options, so a second newscore call could be
    newscore=Eg:3+2*6,Eg2:Eg^2 (= Eg squared = 900.0). This can be used in conjunction with statfilter,
    e.g. newscore=UDif:UHS-UP statfilter=UDif>1.

Commandline:
    ## Basic compilation options ##
    outfile=FILE    : Name of output file. [slimdisc_results.csv]
    dirlist=LIST    : List of directories from which to extract files (wildcards OK) [./]
    compile=T/F     : Compile motifs from SliMDisc rank files into output file. (False=index only) [True]
    append=T/F      : Append file rather than over-writing [False]
    slimranks=X     : Maximum number of SlimDisc ranks to exract from any given dataset [5000]
    rerank=X        : Re-ranks according to RScore (if expect=T) and only outputs top X new ranks (if > 0) [5000]
    rankstat=X      : Stat to use to re-rank data [RScore]
    motific=T/F     : Recalulate IC using PRESTO. Used for re-ranking. OldIC also output. [False]
    lencorrect=T/F  : Implements crude length correction in RScore [False]
    delimit=X       : Change delmiter to X [,]

    ## Advanced compilation options ##
    subdir=T/F      : Whether to search subdirectories for rank files [True]
    webid=LIST      : List of SLiMDisc webserver IDs to compile. (Works only on bioware!) []
    slimversion=X   : SLiMDisc results version for compiled output [1.4]

    ## Additonal statistics ##
    abschg=T/F      : Whether to output number of charged positions (KRDE) [True]
    netchg=T/F      : Whether to output net charge of motif (KR) - (DE) [True]
    balchg=T/F      : Whether to output the *balance* of charge (netNT - netCT) [True]
    ailmv=T/F       : Whether to output if all positions in the motif are A,I,L,M or V. [True]
    aromatic=T/F    : Whether to output count of F+Y+W [True]
    phos=T/F        : Whether to output potential phosphorylated residues X (none) or [S][T][Y], if present [True]
    expect=T/F      : Calculate min. expected occurrence of motif in search dataset [True]
    zscore=T/F      : Calculate z-scores for each motif using the entire dataset (<=slimranks) [True]
    newscore=LIST   : Lists of X:Y, create a new statistic X, where Y is the formula of the score. []
    custom=LIST     : Calulate Custom score as a produce of stats in LIST []

    ## Additional calculations to make: see PRESTO for additional relevant commandline options ##
    slimsa=T/F      : Calculate SA information for SLiMDisc Results [True]
    winsa=X         : Number of aa to extend Surface Accessibility calculation either side of motif [0]
    slimhyd=T/F     : Calculate Eisenbeg Hydophobicity for SLiMDisc Results [True]
    winhyd=X        : Number of aa to extend Eisenberg Hydrophobicity calculation either side of motif [0]
    slimcons=T/F    : Calculate Conservation stats for SLiMDisc results [False]
                    - See PRESTO conservation options. (NB. consamb does nothing.)
    slimchg=T/F     : Calculate selected charge statistics (above) for occurrences in addition to pattern [False]
    slimfold=T/F    : Calculate disorder using FoldIndex over the internet [False]
    slimiup=T/F     : Calculate disorder using local IUPred [True]
    windis=X        : Number of aa to extend disorder prediction each side of occurrence [0]
    iucut=X         : Cut-off for IUPred results [0.2]
    iumethod=X      : IUPred method to use (long/short) [short]
    iupath=PATH     : The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    percentile=X    : Percentile steps to return in addition to mean [0]

    ## Collation and Extraction of specific results ##
    index=T/F       : Whether to create index files (slimpicks.*.index) for proteins, motifs and datasets [True]
    bigindex=T/F    : Whether to use the special makeBigIndexFiles() method [False]
    fullpath=T/F    : Whether to use full path (else relative) for dataset index [True]
    slimpath=PATH   : Path to place (or find) index files. *Cannot be used for extraction if fullpath=F* [./]
    slimlist=LIST   : List (A,B,C) or FILE containing list of SLiMs (motifs) to extract []
    protlist=LIST   : List (A,B,C) or FILE containing list of proteins for which to extract results  []
    datalist=LIST   : List (A,B,C) or FILE containing list of datasets for which to extract results []
    strict=T/F      : Only extract protein/occurrence details for those proteins in protlist [False]
                      (False = extract details for all proteins in datalist datasets containing slimlist motifs)
    outdir=PATH     : Directory into which extracted data will be placed. [./]
    picksid=X       : Outputs an extra 'PicksID' column containg the identifier X []
    inputext=LIST   : List of file extensions for original input files. (Should be in same dir as *.rank, or one dir above)
                      [dat,fas,fasta,faa]
    indexre=LIST    : List of alternative regular expression patterns to try for index retrieval []
                      - ipi     : 'ipi_HUMAN__(\S+)-*\d*=(\S.+)',           # IPI Human sequence
                      - ipi_sv  : '^ipi_HUMAN__([A-Za-z0-9]+)-*\d*=(\S.+)', # IPI Human UniProt splice variant
                      - ft      : '^(\S+)_HUMAN=(\S+)',                     # SLiMDisc FullText (UniProt format) retrieval
                      - ft_sv   : '^([A-Za-z0-9]+)-*\d*_HUMAN=(\S+)'        # SLiMDisc FullText (UniProt format) splice variant

    ## Additional Output for Extracted Motifs ##
    occres=FILE     : Output individual occurrence data in FILE [None]
    extract=T/F     : Extract additional data for motifs [True if datasets/SLiMs/accnums given, else False]
    motifaln=T/F    : Produce fasta files of local motif alignments [True]
    flanksize=X     : Size of sequence flanks for motifs [30]
    xdivide=X       : Size of dividing Xs between motifs [10]
    datout=FILE     : Extract UniProt entries from parent proteins where possible into FILE [uniprot_extract.dat]
    unitab=T/F      : Make tables of UniProt data using rje_uniprot.py [True]
    ftout=FILE      : Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [None]
    unipaths=LIST   : List of additional paths containing uniprot.index files from which to look for and extract features ['']
    peptides=T/F    : Peptide design around discovered motifs [False]

    ## Additional Output for Proteins ##
    proteinaln=T/F  : Search for alignments of proteins containing motifs and produce new file containing motifs [True]
    gopher=T/F      : Use GOPHER to generate missing orthologue alignments in alndir - see gopher.py options [False]
    alndir=PATH     : Path to alignments of proteins containing motifs [./] * Use forward slashes (/) [Gopher/ALN/]
    alnext=X        : File extension of alignment files, accnum.X [orthaln.fas]

    ## Advanced Filtering Options ##
    statfilter=LIST : List of stats to filter (remove matching motifs) on, consisting of X*Y where:
                      - X is an output stat (the column header),
                      - * is an operator in the list >, >=, !=, =, >= ,<    !!! Remember to enclose in "quotes" for <> !!!
                      - Y is a value that X must have, assessed using *.
                      This filtering is crude and may behave strangely if X is not a numerical stat!
    zfilter=T/F     : Calculate the Z-score on the filtered dataset (True) or the whole dataset (False) [False]
    rankfilter=T/F  : Re-ranks the filtered dataset (True) rather than the whole (pre-filtered) dataset (False) [True]
                      - NB. If zfilter=T then rankfilter=T.

    ## Old/obselete options ##
    advprob=T/F     : Calculate advanced probability based on actual sequences containing motifs [False] #!# Not right yet!! #!#
    advmax=X        : Max number of sequences to use computationally intensive advanced probability [35]

    *** See RJE_UNIPROT options for UniProt settings ***

Uses general modules: copy, math, os, re, string, sys, time
Uses RJE modules: rje, gopher, presto, rje_disorder, rje_motif, rje_scoring, rje_seq, rje_sequence, rje_uniprot
Other modules needed: rje_blast, rje_dismatrix, rje_pam 
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, math, os, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import gopher, presto, slimpicks_index
import rje, rje_disorder, rje_motif, rje_scoring, rje_seq, rje_sequence, rje_uniprot
#########################################################################################################################
### History
# 3.0 - Stable version for use with SLiMDisc webserver. Tidied code and added newscore. See V2.9 for previous history.
#########################################################################################################################
### Major Functionality to Add
# [ ] : Modify Expect calculation to account for similarities and fact that found at least once?
# [ ] : In interactive mode (i=1), ask user to confirm each higher level directory to go through (e.g. DIR/*/*.rank).
# [ ] : Check for multiple datasets with same name.
# [ ] : Add output of features in given sequence that overlap with/near motif (as in PIRATE locft)?
# [ ] : Add option to count Histidine as charged?
# [ ] : Add output of motifs integrated into UniProt extracted data (rather than PRESTO UniProt output)
# [~] : Tidy implementation of Gopher: Handle uniprot input files and redirect output correctly.
# [ ] : Replace global attributes from main run method with object attribute - safer!
# [ ] : Move index methods into slimpicks_index.py
# [ ] : Check whether more presto and rje_motif methods can be used.
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit) = ('SLiMPickings', '3.0', 'January 2007')  
    description = 'SLiMDisc results compilation and extraction'
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
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show PRESTO commandline options?'):
                out.verbose(-1,4,text=presto.__doc__)
            if rje.yesNo('Show RJE_UNIPROT commandline options?'):
                out.verbose(-1,4,text=rje_uniprot.__doc__)
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
alt_index_re = {'ipi':'^ipi_HUMAN__(\S+)=(\S.+)',   # IPI Human sequence
                'ipi_sv':'^ipi_HUMAN__([A-Za-z0-9]+)-*\d*=(\S.+)',  # IPI Human UniProt splice variant
                'ft':'^(\S+)_HUMAN=(\S+)',       # SLiMDisc FullText (UniProt format) retrieval
                'ft_sv':'^([A-Za-z0-9]+)-*\d*_HUMAN=(\S+)'       # SLiMDisc FullText (UniProt format) splice variant
                }
pstats = ['Dataset', 'SeqNum', 'FullMST', 'TotalAA', 'Rank', 'Score', 'Pattern', 'Occ', 'IC', 'Norm', 'Sim', 'MST',
          'UHS', 'UP', 'Expect', 'Prob', 'RScore', 'OldRank', 'ZScore', 'OldIC', 'AbsChg', 'NetChg', 'BalChg', 'AILMV',
          'PicksID','Aromatic','Phos']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SlimPicker Class:                                                                                       #
#########################################################################################################################
class SlimPicker(rje.RJE_Object):     
    '''
    SlimPicker Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of output file [slimdisc_results.csv]
    - OutDir = Directory into which extracted data will be placed [./]
    - DatOut = Extract UniProt entries from parent proteins where possible into FILE [uniprot_extract.dat]
    - OccRes = Output individual occurrence data in FILE [None]
    - SlimPath = Path to place (or find) index files. *Cannot be used for extraction if fullpath=F* [./]
    - FTOut = Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [None]
    - UniPaths = List of paths containing uniprot.index files from which to look for and extract features [/home/richard/Databases/UniProt/]
    - AlnDir = Path to alignments of proteins containing motifs [./] * Use forward slashes (/) [Gopher/ALN/]
    - AlnExt = File extension of alignment files, accnum.X [aln.fas]
    - SlimVersion = SLiMDisc results version for compiled output [1.4]
    - Extract = Stores command for extraction of additional data for motifs, for setting default value in self.pickSetup()
    - RankStat = Stat to use to re-rank data [RScore]
    - PicksID = Outputs an extra 'PicksID' column containg the identifier X []
    
    Opt:boolean
    * Basic Stats Output *
    - Compile = Compile motifs from SliMDisc rank files into output file. (False=index only) [True]
    - Extract = Extract additional data for motifs [True if datasets/SLiMs/accnums given, else False]
    - AbsChg = Whether to output number of charged positions (KRDE) [True]
    - NetChg = Whether to output net charge of motif (KR) - (DE) [True]
    - BalChg = Whether to output the *balance* of charge (netNT - netCT) [True]
    - AILMV = Whether to output if all positions in the motif are A,I,L,M or V. [True]
    - Aromatic = Whether to output count of F+Y+W [True]
    - Phos = Whether to output potential phosphorylated residues X (none) or [S][T][Y], whichever are present [True]
    - Expect = Calculate min. expected occurrence of motif in search dataset [True]
    - AdvProb = Calculate advanced probability based on actual sequences containing motifs [True]
    - SubDir = Whether to search subdirectories for rank files [True]
    - MotifIC = Recalulate IC using PRESTO. Used for re-ranking. OldIC also output. [False]
    * PRESTO Stats Output *
    - SlimSA = Calculate SA information for SLiMDisc Results [True]
    - SlimHyd = Calculate Eisenbeg Hydophobicity for SLiMDisc Results [True]
    - SlimCons = Calculate Conservation stats for SLiMDisc results [False]
    - SlimChg = Calculate selected charge statistics (above) for occurrences in addition to pattern [False]
    - SlimFold = Calculate disorder using FoldIndex over the internet [False]
    - SlimIUP = Calculate disorder using local IUPred [True]
    * Specific Results extraction *
    - Index = Whether to create index files (slimpicks.*.index) for proteins, motifs and datasets [True]
    - BigIndex = Whether to use the special makeBigIndexFiles() method
    - FullPath = Whether to use full path (else relative) for dataset index [True]
    - MotifAln = Produce fasta files of local motif alignments [True]
    - UniTab = Make tables of UniProt data using rje_uniprot.py [True]
    - Strict = Only extract protein details for those proteins in protlist [False]
              (False = extract details for all proteins in datalist datasets containing slimlist motifs)
    - ProteinAln = Search for alignments of proteins containing motifs and produce new file containing motifs [True]
    - Gopher = Use GOPHER to generate missing orthologue alignments in outdir/Gopher - see gopher.py options [False]
    * Special *
    - ZScore = Calculate z-scores for each motif using the entire dataset (<=slimranks) [True]
    - ZFilter = Calculate the Z-score on the filtered dataset (True) or the whole dataset (False) [False]
    - RankFilter = Re-ranks the filtered dataset (True) rather than the whole (pre-filtered) dataset (False) [True]
    - LenCorrect = Implements crude length correction in RScore [False]

    Stat:numeric
    - SlimRanks = Maximum number of SlimDisc ranks to exract from any given dataset (0=all) [0]
    - Percentile = Percentile steps to return in addition to mean [25]
    - FlankSize = Size of sequence flanks for motifs [30]
    - XDivide = Size of dividing Xs between motifs [10]
    - AdvMax = Max number of sequences to use computationally intensive advanced probability [35]
    - ReRank = Re-ranks according to RScore (if expect=T) and only outputs top X new ranks (if > 0) [0]

    List:list
    - DirList = List of directories from which to extract files (wildcards OK) [./]
    * Specific Results Extraction *
    - SlimList = List (A,B,C) or FILE containing list of SLiMs (motifs) to extract []
    - ProtList = List (A,B,C) or FILE containing list of proteins for which to extract results []
    - DataList = List (A,B,C) or FILE containing list of datasets for which to extract results []
    - InputExt = List of file extensions for original input files. (Should be in same dir as *.rank, or one dir above) [dat,fas]
    - FileList = List of actual rank files to check out (paths to these files) #~# Setup in run() #~#
    - Specials = List of special (SlimX) statistics to calculate. #~# Setup in run() #~#
    - IndexRE = List of alternative regular expression patterns to try for index retrieval []
    ## Advanced Filtering/Ranking Options ##
    - StatFilter = List of stats to filter on, consisting of X*Y where:
          - X is an output stat (the column header),
          - * is an operator in the list >, >=, =, <= ,<, !=
          - Y is a value that X must have, assessed using *.
          This filtering is crude and may behave strangely if X is not a numerical stat!
    - Custom = Calulate Custom score as a produce of stats in LIST []
    - WebID = List of SLiMDisc webserver IDs to compile. (Works only on bioware!) []
    - NewScore = self.dict['NewScore'] keys() in order they were read in

    Dict:dictionary
    - IndexFiles = Dictionary of index file names and types {'Slim':PATH, etc.} #~# Setup in run() #~#
    - StatFilter = Dictionary of stat filters made from self.list['StatFilter']
    - NewScore = dictionary of {X:Y} for new statistic X, where Y is the formula of the score. []

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object'''
        ### Basics ###
        self.infolist = ['Name','DatOut','OccRes','OutDir','SlimPath','FTOut','UniPaths','AlnDir','AlnExt','SlimVersion',
                         'Extract','RankStat','PicksID']
        self.optlist = ['Compile','AbsChg','NetChg','BalChg','AILMV','Expect','Index','FullPath','Expect',
                        'BigIndex','SlimSA','SlimHyd','SlimCons','SlimChg','SubDir','SlimIUP','SlimFold','MotifAln',
                        'UniTab','Strict','AdvProb','ProteinAln','Gopher','ZScore','Extract','MotifIC','ZFilter',
                        'RankFilter','Aromatic','Phos','LenCorrect']
        self.statlist = ['Percentile','SlimRanks','FlankSize','XDivide','AdvMax','ReRank']
        self.listlist = ['DirList','SlimList','ProtList','DataList','InputExt','IndexRE','StatFilter','Custom','WebID',
                         'NewScore']
        self.dictlist = ['StatFilter','NewScore']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Name':'slimdisc_results.csv','Path':rje.makePath('./',return_blank=False),'RankStat':'RScore',
                      'OutDir':rje.makePath('./',return_blank=False),'SlimPath':rje.makePath('./',return_blank=False),
                      'AlnDir':'Gopher/ALN/','AlnExt':'orthaln.fas','FTOut':'None','UniPaths':'','PicksID':''})
        self.setStat({'SlimRanks':5000,'ReRank':5000,'Percentile':0,'FlankSize':30,'XDivide':10,'AdvMax':35})
        for o in ['BigIndex','Strict','AdvProb','SlimCons','Gopher','SlimFold','MotifIC','LenCorrect']:
            self.opt[o] = False
        self.list['InputExt'] = ['dat','fas','fasta','faa']
        #X#self.list['DirList'] = [rje.makePath('./',return_blank=False)]
        self.list['PStats'] = pstats[0:]    # List of P-Statistics (pattern-only, no occurrence data needed)
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
                ### Compilation Options ### 
                self._cmdRead(cmd,type='info',att='Name',arg='outfile')
                self._cmdReadList(cmd,type='info',attlist=['PicksID','RankStat','SlimVersion','Extract'])
                self._cmdReadList(cmd,type='list',attlist=['DirList','WebID','Custom'])
                self._cmdReadList(cmd,'int',['SlimRanks','ReRank','AdvMax'])
                self._cmdReadList(cmd,'opt',['Compile','SubDir','AbsChg','NetChg','BalChg','AILMV','Expect','ZScore',
                                             'AdvProb','ZFilter','RankFilter','Aromatic','Phos','LenCorrect','Extract',
                                             'MotifIC'])
                ### PRESTO Stats Options ### 
                self._cmdRead(cmd,type='info',att='OccRes')
                self._cmdReadList(cmd,'opt',['SlimSA','SlimHyd','SlimCons','SlimChg','SlimIUP','SlimFold'])
                self._cmdRead(cmd,type='int',att='Percentile')
                self._cmdRead(cmd,type='list',att='StatFilter')
                ## Collation and Extraction of specific results ##
                self._cmdReadList(cmd,'path',['OutDir','SlimPath','AlnDir'])
                self._cmdReadList(cmd,'info',['DatOut','FTOut','AlnExt'])
                self._cmdReadList(cmd,'opt',['Index','BigIndex','FullPath','Strict','MotifAln','UniTab','ProteinAln',
                                             'Gopher'])
                self._cmdReadList(cmd,'list',['SlimList','ProtList','DataList','IndexRE','InputExt'])
                self._cmdReadList(cmd,'int',['FlankSize','XDivide'])
                self._cmdRead(cmd,type='clist',att='UniPaths')
                ## Old/obselete methods ##
                self._cmdRead(cmd,type='list',att='RankMeans')
                ## Special NewScore ##
                self._cmdRead(cmd,'cdictlist','NewScore')
                continue
                if cmd.lower().find('newscore=') == 0 and cmd.lower().find(':') > 0:
                    new = cmd[cmd.find('=')+1:cmd.find(':')]
                    self.dict['NewScore'][new] = cmd[cmd.find(':')+1:]
                    if new in self.list['NewScore']:
                        self.list['NewScore'].remove(new)
                    self.list['NewScore'].append(new)
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        #X#print self.dict['NewScore'], self.list['NewScore']
        try:
            for r in range(len(self.list['RankMeans'])):
                self.list['RankMeans'][r] = string.atoi(self.list['RankMeans'][r])
        except:
            return               
#########################################################################################################################
    ### <2> ### Revised controlling run method                                                                          #
#########################################################################################################################
    def run(self):      ### Main controlling run method for new Slim Pickings 1.0
        '''Main controlling run method for new Slim Pickings 1.0.'''
        try:
            ### Setup Index Files etc. ###
            _stage = 'Setup Index Files'
            ## Check need for index files & create name dictionary ##
            need_index = False      # Whether index files are needed for extraction
            index_exists = False    # Whether any index files exist
            use_index = True        # Whether index files should be used for extraction (if needed)
            self.dict['IndexFiles'] = {}
            for type in ['Slim','Prot','Data']:
                if self.list['%sList' % type]:
                    #X#print type, self.list['%sList' % type]
                    need_index = True
                self.dict['IndexFiles'][type] = self.info['SlimPath'] + 'slimpicks.%s.index' % type.lower()
                self.dict['IndexFiles']['tmp%s' % type] = self.info['SlimPath'] + '%s.%s.tmp' % (type,rje.randomString(10))
                if os.path.exists(self.dict['IndexFiles'][type]):
                    index_exists = True
                elif need_index:
                    self.log.printLog('#INDEX','Index file %s missing.' % self.dict['IndexFiles'][type])
                    use_index = False
            ## Create index files? ##
            if not self.opt['Index'] and need_index and not use_index:    # Only partial extract #
                if self.stat['Interactive'] < 0 or rje.yesNo('(Re)generate index files for partial results extraction?'):
                    self.opt['Index'] = True
                else:
                    self.log.printLog('#INDEX','Index files missing. Will scan all rank files.')
            if self.opt['Index'] and index_exists:  # Index will be over-written #
                if need_index and use_index and (self.stat['Interactive'] < 0 or rje.yesNo('Use existing Index files exist for extraction?')):
                    self.log.printLog('#INDEX','Index files exist. Cancelled index generation.')
                    self.opt['Index'] = False
                elif self.stat['Interactive'] < 0 or rje.yesNo('One or more index files exist. Overwrite?'):
                    self.log.printLog('#INDEX','One or more index files exist. Will be overwritten.')
                else:
                    self.log.printLog('#INDEX','One or more index files exist. Cancelled index generation.')
                    self.opt['Index'] = False

            ### Build Index ###
            self.list['FileList'] = []
            if self.opt['Index'] and self.makeIndexFiles():
                use_index = True

            ### Rank file list ###
            _stage = 'Setup Rank file list'
            if need_index and use_index:    ### Get list of files from index
                self._getRankFilesFromIndexFiles()
            elif not need_index and not self.opt['Compile']:
                return
            else:
                if not self.list['FileList']:
                    self._getAllRankFiles()
                self.log.printLog('\n#FILES','%s rank files identified for data extraction.' % rje.integerString(len(self.list['FileList'])))

            ### Get SLiM Data from Rank Files ###
            self.picks()
            
        except:
            self.log.errorLog('Fatal error during main SlimPicker.run(%s)' % _stage,quitchoice=True)
            raise
#########################################################################################################################
    ### <3> ### SLiM Pickings Index Methods                                                                             #
#########################################################################################################################
# Currently, three index files will be used:
# slimpicks.data.index: Dataset=PATH
# slimpicks.slim.index: Pattern=Dataset List
# slimpicks.prot.index: ProteinID=Dataset List
#########################################################################################################################
    def makeIndexFiles(self):    ### Generates index files from rank files
        return slimpicks_index.makeIndexFiles(self)
#########################################################################################################################
    def makeBigIndexFiles(self):    ### Generates index files from rank files
        return slimpicks_index.makeBigIndexFiles(self)
#########################################################################################################################
    def _getListFromIndex(self,type,searchlist):    ### Returns a list of returned index values for searchlist and dict
        '''
        Returns a list of returned index values for searchlist.
        >> type:str = Type of slim index search (Data,Prot or Slim)
        >> searchlist:list of target strings
        << (returnlist,replacedict) = list of returned targets and a dictionary of {search:[found]} for alternative REs.
        '''
        return slimpicks_index._getListFromIndex(self,type,searchlist)
#########################################################################################################################
    def _getRankFilesFromIndexFiles(self):   ### Populates self.list['FileList'] using index files and self.lists
        '''
        Populates self.list['FileList'] using index files and self.lists. This method is only called if one of the
        Data/Prot/Slim List lists is populated. This method works by the following:
        1. Make a list of OK datasets from ProtList and SlimList - only keep common datasets
        2. Extract full path from those datasets that also occur in DataList
        - if any list is empty, then it is ignored as a potential filter!
        '''
        return slimpicks_index._getRankFilesFromIndexFiles(self)
#########################################################################################################################
    def _getAllRankFiles(self):  ### Populates self.list['FileList'] with all rank files
        '''Populates self.list['FileList'] with all rank files.'''
        self.list['FileList'] = []
        for id in self.list['WebID']:
            tmproot = '/home/slimdisc/public_html/SLiMDisc_dev4/slimdiscBrowser/datasets/SLiMDisc_temp_'
            self.list['DirList'].append(tmproot+id)
        if not self.list['DirList']:
            self.list['DirList'] = ['./']
        for dir in self.list['DirList']:
            sfiles = rje.getFileList(callobj=self,folder=rje.makePath(dir,return_blank=False),filelist=['*.rank'],subfolders=self.opt['SubDir'])
            for slimdisc in sfiles:
                if slimdisc[-9:] != '.dat.rank':
                    self.list['FileList'].append(slimdisc)
            self.log.printLog('\r#FILES','Rank files extracted from %s: %s total.' % (dir,rje.integerString(len(self.list['FileList']))))
#########################################################################################################################
    ### <4> ### Main SLiMDisc results complation Methods                                                                #
#########################################################################################################################
    def picks(self):      ### SLiMDisc Results compiler
        '''SLiMDisc Results compiler.'''
        ### Setup ###
        try:
            if not self.pickSetup():    ### Sets up variables and output etc. for picks() SLiMDisc Results compiler
                return False
            delimit = self.info['Delimit']                  # text delimiter for output
            presto_calc = self.opt['PRESTO Calculations']   # whether to perform occurrence-specific PRESTO calculations
            extract = self.opt['Extract']                   # whether to extract occurrence data for results 
            headers = self.list['Headers']                  # list of headers for compiled output
            special_stats = self.list['Special Stats']      # list of "special statistics" to calculate
            slim_presto = self.obj['SlimPRESTO']            # PRESTO object for storing motifs and making calculations
            slim_seq = self.obj['SlimSeq']                  # SeqList object for handling search dataset sequences
            presto_stats = self.list['PRESTO Stats']        # PRESTO occurrence stats headers
        except:
            self.log.errorLog('Error in SlimPicker.picks() Setup',printerror=True,quitchoice=True)
            return False

        ### Main Processing ###
        try:
            ### Setup ###
            _stage = 'Setup'
            fx = 0  # Number of files
            rx = 0  # Number of ranked motifs analysed
            px = 0  # Number of proteins for which details have been extracted
            uni_extract = []    # List of accession numbers to extract from UniProt
            slim_ft = {}        # Dictionary of {AccNum:[Extra features (SLiMs) to add to table}
           
            ### Work Through Files ###
            for slimdisc in self.list['FileList'][0:]:
                fx += 1
                dataset = rje.baseFile(slimdisc,True)

                ## Proteins and Motifs ##
                _stage = 'Proteins and Motifs'
                datlines = self.loadFromFile(filename=string.replace(slimdisc,'rank','dat.rank'),v=1)
                motiflist = []      # List of motif strings to consider (filtered by motif and protein) *in rank order*
                slim_presto.list['Motifs'] = []
                slim_prot = {}      # Dictionary of ID number (str): shortName()
                slim_motifs = {}    # Dictionary of motif: list of occurrances (ID:Pos)
                extract_prots = []  # List of proteins to extract details for *in input data order*
                readprot = False
                readmot = False
                ## Process File ##
                for line in datlines:
                    if rje.matchExp('^(#Proteins)',line):
                        readprot = True
                    elif rje.matchExp('^(#Motif)',line):
                        readmot = True
                        readprot = False
                    elif readprot and rje.matchExp('^(\d+)\s+(\S+)',line):
                        dmatch = rje.matchExp('^(\d+)\s+(\S+)',line)
                        slim_prot[dmatch[0]] = dmatch[1]
                    elif readmot and rje.matchExp('^(\d+\s+\S+\s+\S.+)$',line):
                        dmatch = string.split(rje.matchExp('^(\d+\s+\S+\s+\S.+)$',line)[0])
                        if len(dmatch) > 2 and (string.atoi(dmatch[0]) <= self.stat['SlimRanks'] or self.stat['SlimRanks'] <= 0):
                            ## Check motif and proteins ##
                            if self.list['SlimList'] and dmatch[1] not in self.list['SlimList']:    # Filtered by motif
                                continue
                            protfound = False
                            for occ in dmatch[2:]:
                                myprot = slim_prot[rje.matchExp('(\d+):(\d+)',occ)[0]]
                                if myprot not in extract_prots:
                                    px += 1
                                    extract_prots.append(myprot)
                                if not self.list['ProtList'] or myprot in self.list['ProtList']:     
                                    protfound = True
                            if not extract_prots or not protfound:   # Filtered by protein
                                continue
                            ## Add ##
                            rx += 1
                            motiflist.append(dmatch[1])
                            slim_motifs[dmatch[1]] = dmatch[2:]
                            slim_presto._addMotif(name=dmatch[1],seq=dmatch[1])
                            self.log.printLog('\r#DAT','Dataset %s: Processing (%d patterns, %d proteins)' % (dataset,len(motiflist),len(extract_prots)),log=False,newline=False)
                ## Check filtered Motifs remain ##
                #X#print slimdisc,motiflist,extract_prots,slim_motifs
                if motiflist:
                    self.log.printLog('\r#DAT','Dataset %s: Processing (%d patterns, %d proteins)' % (dataset,len(motiflist),len(extract_prots)))
                else:
                    self.log.printLog('\r#DAT','Dataset %s: No patterns to process after protein/motif filtering.' % dataset)
                    continue

                ### P-Stat Calculations ### >-------------------------------------------------------------<### P-Stat ###
                ## Basic Compilation Stats ##
                _stage = 'Basic Stats'
                ranklines = self.loadFromFile(filename=slimdisc,v=1)
                self.dict['BasicStats'] = {}  # Dictionary of {motif:{stat:value}}
                fullset = {'Dataset':dataset}
                rankheads = []
                for line in ranklines[1:]:
                    ## General dataset stats ##
                    if line[0] == '#':  ## Comment line
                        if rje.matchExp('No. of proteins:\s+(\d+)',line):
                            fullset['SeqNum'] = rje.matchExp('No. of proteins:\s+(\d+)',line)[0]
                        if rje.matchExp('No. of residues:\s+(\d+)',line):
                            fullset['TotalAA'] = rje.matchExp('No. of residues:\s+(\d+)',line)[0]
                        if rje.matchExp('Overall MST:\s+(\d\S+)',line):
                            fullset['FullMST'] = rje.matchExp('Overall MST:\s+(\d\S+)',line)[0]
                        continue
                    ## SLiMDisc rank file stats headers ###
                    if line.find('Rank') == 0:
                        rankheads = string.split(rje.chomp(line))
                    ## SLiMDisc Motif stats ##
                    mydata = rje.matchExp('^\((\d+)\)\s+(\S+)\s+(\S+)\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)',rje.chomp(line))
                    if not mydata or mydata[2] not in motiflist:    # Not interested
                        continue
                    if 'MST' in rankheads:
                        mydata = [mydata[0]] + string.split(line)[1:]
                    pattern = mydata[2]
                    #X#self.dict['BasicStats'][pattern] = {'Pattern':pattern}
                    self.dict['BasicStats'][pattern] = {}
                    for key in fullset.keys():
                        self.dict['BasicStats'][pattern][key] = fullset[key]
                    for key in ['Rank','Score','Pattern','Occ','IC','Norm','Sim','MST','UHS','UP','Expect','Prob']:
                        if key in rankheads:
                            try:
                                self.dict['BasicStats'][pattern][key] = mydata[rankheads.index(key)]
                            except:
                                self.log.errorLog('Bugger')
                                print key
                                print rankheads
                                print mydata
                                print line
                                raise
                    ### PicksID ###
                    if self.info['PicksID']:
                        self.dict['BasicStats'][pattern]['PicksID'] = self.info['PicksID']

                ### OldStats
                for pattern in motiflist:
                    self.dict['BasicStats'][pattern]['OldRank'] = self.dict['BasicStats'][pattern]['Rank']
                    self.dict['BasicStats'][pattern]['OldIC'] = self.dict['BasicStats'][pattern]['IC']

                ### Input Sequence Files ###
                _stage = 'Input Sequence Files'
                inputseq = None     # Original input sequence file
                if 'fasta' in self.list['InputExt']:
                    while 'fasta' in self.list['InputExt']:
                        self.list['InputExt'].remove('fasta')
                    self.list['InputExt'].append('fasta')
                for ext in self.list['InputExt']:   # Look at each extension in directory above and same directory
                    for tryme in ['%s.%s' % (string.join(string.split(slimdisc,os.path.sep)[:-1],os.path.sep),ext),string.replace(slimdisc,'rank',ext)]:
                        #X#self.deBug('%s:%s' % (tryme,os.path.exists(tryme)))
                        if not inputseq and os.path.exists(tryme):
                            inputseq = tryme
                        #X#self.log.printLog('#TRY','Try: %s; Input: %s; Exists: %s' % (tryme,inputseq,os.path.exists(tryme)))
                ## Make Fasta file for Gopher? ##
                if self.opt['Gopher'] and inputseq[-4:] != '.fas':
                    self.log.printLog('#FAS','Making new fasta file for GOPHER')   #!# Change to test if Fasta #!#
                    fcmd = ['seqin=%s' % inputseq,'seqout=%s' % string.replace(slimdisc,'rank','fas'),'reformat=fasta','memsaver=T']
                    rje_seq.SeqList(self.log,fcmd)
                    inputseq = string.replace(slimdisc,'rank','fas')

                ## Expected Calculation AA Frequency file ##
                fasta = string.replace(slimdisc,'rank','fasta')     # TEIRESIAS input #
                if not os.path.exists(fasta):
                    fasta = inputseq

                ### Generate motif_seq dictionary from slim_prot and slim_motifs ###
                _stage = 'Motif:Sequence dictionary'
                self.dict['FullSet'] = fullset
                self.dict['SlimMotif'] = slim_motifs        # Dictionary of motif: list of occurrances (ID:Pos)
                self.dict['SlimProt'] = slim_prot           # Dictionary of ID number (str): shortName()
                self.motifOccDict(inputseq,extract_prots)   # Creates self.dict['Motif PNames'] and self.dict['Motif Occ']
                motif_pnames = self.dict['Motif PNames']
                                
                ### Expected Occurrence ### 
                _stage = 'Expected Occurrence'  # Also calculates RScores #
                prob_occ = {}
                if (self.opt['Compile'] or self.stat['ReRank'] > 0) and self.opt['Expect']:
                    prob_occ = self.expectation(slim_presto,motif_pnames,fasta)  # For OccRes

                ### Additional motif stats ###
                if self.opt['Compile'] or self.stat['ReRank'] > 0:
                    _stage = 'Additional Stats'
                    for pattern in motiflist:
                        self.log.printLog('\r#STAT','SLiM Stats: %.f%%' % (100.0*motiflist.index(pattern)/len(motiflist)),log=False,newline=False)
                        patstats = self.patternStats(pattern) # Returns a dictionary of stats (AbsChg,NetChg,BalChg,AILMV)
                        for key in patstats.keys():    
                            if self.opt[key]:
                                self.dict['BasicStats'][pattern][key] = patstats[key]
                    self.log.printLog('\r#STAT','SLiM Stats: Complete.')

                ### Custom Scores ###
                for new in self.list['NewScore'][0:]:   # self.dict['NewScore'] keys() in order they were read in
                    if not self.dict['CustomOcc'][new]:
                        self.customScore(motiflist,new)     
                
                ### Post P-Stat Filtering, Ranking and Z-Scores ### >--------------------------------<### Post P-Stat ###
                if self.info['RankStat'] in self.list['PStats'] and not self.opt['FilterOcc'] and not self.opt['ZFilter'] and not self.opt['RankFilter']:
                    ## Z-Scores ##
                    self.zScores(motiflist)
                    ## Re-ranking and reduction by slimranks ##
                    motiflist = self.reRank(motiflist)
                    ## P-Stat Filtering ##
                    motiflist = self.statFilter(motiflist,type='P')
                elif self.info['RankStat'] in self.list['PStats'] and not self.opt['ZFilter'] and self.opt['RankFilter']:
                    ## Z-Scores ##
                    self.zScores(motiflist)
                    ## P-Stat Filtering ##
                    motiflist = self.statFilter(motiflist,type='P')
                elif self.info['RankStat'] in self.list['PStats'] and not self.opt['FilterOcc'] and self.opt['ZFilter']:
                    ## P-Stat Filtering ##
                    motiflist = self.statFilter(motiflist,type='P')
                    ## Z-Scores ##
                    self.zScores(motiflist)
                    ## Re-ranking and reduction by slimranks ##
                    motiflist = self.reRank(motiflist)
                elif self.info['RankStat'] or self.opt['ZFilter']:
                    ## P-Stat Filtering ##
                    motiflist = self.statFilter(motiflist,type='P')

                ### O-Stat Calculations ### >-------------------------------------------------------------<### O-Stat ###

                ### GOPHER Alignment Generation and Alteration of Attributes for prestoStats() and proteinAln() ###
                self.deBug('%s => Gopher=%s' % (inputseq,self.opt['Gopher']))
                if self.opt['Gopher'] and inputseq:     #!# Add to PRESTO #!#
                    #X#self.gopher(inputseq) #!# Add this method to copy existing alignments and perform search otherwise
                    slimpickdir = os.getcwd()
                    gcmd = self.cmd_list + ['gopher=%s' % os.path.abspath(inputseq), 'gnspacc=T', 'orthalign']
                    self.log.printLog('#GCMD',string.join(gcmd))
                    os.chdir(self.info['GopherDir'])
                    try:
                        mygopher = gopher.Gopher(log=self.log,cmd_list=gcmd)
                        mygopher.run()
                    except:
                        self.log.errorLog('Problem with Gopher run!')
                    os.chdir(slimpickdir)

                ### Read and calculate Occurrence data ###
                _stage = 'Occurrence Data'
                if presto_calc and inputseq:
                    presto_seqhit = self.prestoStats(slim_presto,self.dict['Motif Occ'],prob_occ)  ### Populates and returns presto_seqhit list of objects storing data
                    
                ## Calculate statistics for all motif occurrences ##
                _stage = 'Combined Occurrence Data'
                if presto_calc and inputseq and self.opt['Compile']:
                    for motif in slim_presto.list['Motifs']:
                        ## Setup list of hits for this motif ##
                        self.log.printLog('\r#COMB','Combined Occurence Stats: %.f%%' % (100.0*slim_presto.list['Motifs'].index(motif)/slim_presto.motifNum()),log=False,newline=False)
                        pattern = motif.info['Name']
                        if pattern not in motiflist:
                            continue
                        motif_hits = []
                        for seqhit in presto_seqhit:
                            if seqhit.dict['Hits'].has_key(motif):
                                motif_hits += seqhit.dict['Hits'][motif]
                            else:
                                seqhit.dict['Hits'][motif] = []

                        ## Special calculations ##
                        for special in special_stats:   ## Setup during "Main Compilation Setup": should only contain selected stats
                            speclist = []   # List of given stat for each occurrence
                            for hit in motif_hits:
                                hit.stat['Cons'] = hit.stat['ALL_CONS'] #!# Needs better unifying? #!#
                                #X#self.deBug(hit.stat)
                                if special != 'Cons' or hit.stat['ALL_HOM'] > 0:    # If Cons, must have Homologues
                                    speclist.append(string.atof(hit.stat[special]))
                            if speclist:    # Only add to dictionary if actual value: output will insert '' for missing values
                                self.dict['BasicStats'][pattern]['%s_mean' % special] = '%.2f' % (sum(speclist)/len(speclist)) # Mean
                            if self.stat['Percentile'] > 0:
                                if special == 'Hyd':
                                    speclist.sort(reverse=True)
                                else:
                                    speclist.sort()
                                pc = 100.0
                                pclist = speclist[0:]
                                while pc > 0:
                                    while len(pclist) > 1 and len(pclist) > (pc * len(speclist) / 100.0):
                                        pclist.pop(0)
                                    if pclist:
                                        self.dict['BasicStats'][pattern]['%s_pc%d' % (special,pc)]  = '%.2f' % pclist[0]
                                    pc -= self.stat['Percentile']
                    self.log.printLog('\r#COMB','Combined Occurence Stats: Complete.')

                ### Custom Scores ###
                for new in self.list['NewScore'][0:]:   # self.dict['NewScore'] keys() in order they were read in
                    if self.dict['CustomOcc'][new]:
                        self.customScore(motiflist,new)     

                ### Post O-Stat Filtering, Ranking and Z-Scores ### >--------------------------------<### Post O-Stat ###
                rankocc = not self.info['RankStat'] in self.list['PStats']
                if (rankocc or self.opt['FilterOcc']) and not self.opt['ZFilter'] and not self.opt['RankFilter']:
                    ## Z-Scores ##
                    self.zScores(motiflist)
                    ## Re-ranking and reduction by slimranks ##
                    motiflist = self.reRank(motiflist)
                elif rankocc and not self.opt['ZFilter'] and self.opt['RankFilter']:
                    ## Z-Scores ##
                    self.zScores(motiflist)
                ## O-Stat Filtering ##
                motiflist = self.statFilter(motiflist,type='O')
                if not self.opt['ZFilter'] and self.opt['RankFilter']:
                    ## Re-ranking and reduction by slimranks ##
                    motiflist = self.reRank(motiflist)
                elif self.opt['ZFilter'] and self.opt['RankFilter']:
                    ## Z-Scores ##
                    self.zScores(motiflist)
                    ## Re-ranking and reduction by slimranks ##
                    motiflist = self.reRank(motiflist)

                ### Reduce motif_occ to match reduced motiflist ###
                self.reduceMotifOcc(motiflist)      # Use self.dict['Motif Occ']

                ### UniProt Extraction ### >-------------------------------------------------------------<### UniProt ###
                self.dict['SlimFT'] = slim_ft
                self.list['UniExtract'] = uni_extract
                self.slimUniProt(dataset,motiflist)
                slim_ft = self.dict['SlimFT']
                uni_extract = self.list['UniExtract']

                ### OutPut ### >--------------------------------------------------------------------------<### Output ###
                ## Compiled Output ##
                if self.opt['Compile']:
                    _stage = 'Output'
                    cx = 0.0
                    for pattern in motiflist:
                        self.log.printLog('\r#OUT','Compilation of results: %.1f%%' % (cx/len(motiflist)),log=False,newline=False)
                        cx += 100.0
                        self.compileOut(self.info['Name'],self.list['Headers'],delimit,datadict=self.dict['BasicStats'][pattern])
                    self.log.printLog('\r#OUT','Compilation of results: Complete.')
                    if self.info['OccRes'].lower() not in ['','none'] and presto_calc and inputseq:
                        for motif in slim_presto.list['Motifs']:
                            for seqhit in presto_seqhit:
                                seqocc = False
                                for seq in self.dict['Motif Occ'].keys():
                                    if seqhit.info['Name'] in [seq.shortName(),seq.info['AccNum']] and self.dict['Motif Occ'][seq].has_key(motif):
                                        seqocc = True
                                        break
                                if not seqocc:
                                    continue
                                for hit in seqhit.dict['Hits'][motif]:
                                    hitdict = hit.stat
                                    hitdict['Motif'] = motif.info['Name']
                                    hitdict['Len'] = '%d' % len(hit.info['Match'])
                                    hitdict['Hit'] = seqhit.info['Name']
                                    hitdict['Start_Pos'] = '%d' % hit.stat['Pos']
                                    hitdict['End_Pos'] = '%d' % (hit.stat['Pos'] + len(hit.info['Match']) - 1)
                                    hitdict['MatchSeq'] = hit.info['Match']
                                    hitdict['MotifCons'] = hit.stat['ALL_CONS']
                                    hitdict['HomNum'] = hit.stat['ALL_HOM']
                                    hit.stat['Cons'] = hit.stat['ALL_CONS']
                                    if self.obj['SlimPRESTO'].opt['Peptides']:
                                        hitdict['PepSeq'] = hit.info['PepSeq']
                                        hitdict['PepDesign'] = rje_sequence.peptideDetails(hit.info['PepSeq'],self)
                                        pepfas = rje.baseFile(self.info['Name'])+'.peptides.fas'
                                        pepname = '%s-%s-%s-%s %daa (%s)' % (dataset,hitdict['Motif'],hitdict['Hit'],hitdict['Start_Pos'],len(hitdict['PepSeq']),hitdict['PepDesign'])
                                        open(pepfas,'a').write('>%s\n%s\n' % (pepname,hitdict['PepSeq']))
                                    self.compileOut(self.info['OccRes'],presto_stats,delimit,datadict=hitdict)

                        #self.deBug(mydata)

                ### Special Fasta file for Conservation ###
                _stage = 'Protein Fasta Aln'
                if extract and inputseq and self.opt['ProteinAln']:
                    self.proteinAlignments(slim_presto,self.dict['Motif Occ'])

                ### Special Fasta file for dotplots etc. ###
                _stage = 'Motif Fasta Aln'
                if extract:
                    extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
                    extract_seq.opt['ReplaceChar'] = False
                if extract and inputseq and self.opt['MotifAln']:
                    self.motifAlignments(dataset,slim_presto,self.dict['Motif Occ'],motiflist)

                ### End of this dataset ###
                self.log.printLog('#SLIM','Processed: %s files. %s SLiMs. %s proteins.' % (rje.integerString(fx),rje.integerString(rx),rje.integerString(px)),log=False)
                self.log.printLog('>>>','',log=False)
                    
        except:
            self.log.errorLog('Error in SlimPicker.picks Main Processing(%s)' % _stage,printerror=True,quitchoice=True)

        ### UniProt Processing ### >-----------------------------------------------------------------------<### UniProt ###
        try:
            ### Extract uni_extract from UniProt ###
            if (self.info['DatOut'] or self.info['FTOut']) and extract:
                _stage = 'Extract UniProt'
                ## Setup UniProt File ##
                my_entries = []
                uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
                if self.info['DatOut']:
                    uniprot.info['DATOut'] = self.info['DatOut']
                    if self.opt['UniTab']:
                         uniprot.setInfo({'TabOut':'%s.uniprot.%s' % (rje.baseFile(self.info['DatOut']),rje.delimitExt(delimit)),
                                          'LinkOut':'%s.links.%s' %  (rje.baseFile(self.info['DatOut']),rje.delimitExt(delimit))})
                uniprot.list['Extract'] = rje.sortUnique(uni_extract)    
                if self.info['FTOut']:
                    uniprot.opt['MemSaver'] = False
                ## Dat File and Table Output ##
                unipaths = []
                if self.info['UniPaths'] not in ['','None']:
                    unipaths = string.split(self.info['UniPaths'],',')
                #X#self.deBug(unipaths)
                #X#self.deBug(uniprot.info['UniPath'])
                uniprot.run()
                my_entries += uniprot.list['Entry'][0:]
                for path in unipaths:
                    uniprot.opt['Append'] = True
                    uniprot.info['UniPath'] = rje.makePath(path,return_blank=False)
                    uniprot.run()   #!# Add DomTable = Makes a table of domains from uniprot file [False] #!#
                    my_entries += uniprot.list['Entry'][0:]     #!# Add features to existing entries at some point #!#
                ## Features Table incorporating SLiMs ##
                if self.info['FTOut']:
                    # Extract list is AccNums from slim_seq seqlist.
                    _stage = 'UniProt Features'
                    uniprot.opt['Append'] = self.opt['Append']
                    accout = []
                    for entry in my_entries:
                        acc = entry.obj['Sequence'].info['AccNum']
                        if slim_ft.has_key(acc):
                            entry.list['Feature'] += slim_ft[acc]
                        accout.append(acc)
                    for acc in slim_ft.keys():
                        svacc = None
                        if rje.matchExp('^(\S+)\-(\d+)',acc):
                            svacc = rje.matchExp('^(\S+)\-(\d+)',acc)[0]
                        if acc not in accout and (not svacc or svacc not in accout):
                            _entry = rje_uniprot.UniProtEntry(log=self.log,cmd_list=self.cmd_list)
                            _entry.obj['Sequence'].info['AccNum'] = acc
                            _entry.list['Feature'] = slim_ft[acc]
                            my_entries.append(_entry)
                    uniprot.list['Entry'] = my_entries[0:]
                    uniprot.ftTable(self.info['FTOut'])

            ### Finish ###
            self.log.printLog('#SLIM','Processed: %s files. %s SLiMs. %s proteins.' % (rje.integerString(fx),rje.integerString(rx),rje.integerString(px)))
                
        except:
            self.log.errorLog('Error in SlimPicker.picks UniProt Processing(%s)' % _stage,printerror=True,quitchoice=True)
#########################################################################################################################
    def pickSetup(self):      ### Sets up variables and output etc. for picks() SLiMDisc Results compiler
        '''
        Sets up variables and output etc. for picks() SLiMDisc Results compiler.
        Generates additonal class attributes (former method-specific attributes in brackets):
        - self.info['Delimit'] = text delimiter for output (delimit)
        - self.opt['PRESTO Calculations'] = whether to perform occurrence-specific PRESTO calculations (presto_calc)
        - self.list['PRESTO Stats'] = PRESTO occurrence stats headers (presto_calc)
        - self.opt['Extract'] = whether to extract occurrence data for results (alignments, UniProt etc.)
        - self.list['Headers'] = list of headers for compiled output (headers)
        - self.list['Special Stats'] = list of "special statistics" to calculate (special_stats)
        - self.obj['SlimPRESTO'] = slim_presto    # PRESTO object for storing motifs and making calculations
        - self.obj['SlimSeq'] = slim_seq          # SeqList object for handling search dataset sequences
        '''
        ### Setup ###
        try:
            ### General Setup ##
            _stage = 'General'
            self.opt['PRESTO Calculations'] = False # whether to perform occurrence-specific PRESTO calculations (presto_calc)
            self.list['PRESTO Stats'] = []  # PRESTO occurrence stats headers (presto_calc)
            self.list['Headers'] = []       # list of headers for compiled output (headers)
            self.list['Special Stats'] = [] # list of "special statistics" to calculate (special_stats)
            self.obj['SlimPRESTO'] = None   # PRESTO object for storing motifs and making calculations
            self.obj['SlimSeq'] = None      # SeqList object for handling search dataset sequences
            ## Check Files for extraction ##
            if not self.list['FileList']:
                self.log.errorLog('No files for SLiM extraction.',printerror=False)
                return False
            ## Assess extraction option ##
            extract = False
            for type in ['Slim','Prot','Data']:
                if self.list['%sList' % type]:
                    extract = True
            if extract and self.info['Extract'].lower() in ['f','false']:
                self.opt['Extract'] = False
            elif not extract and self.info['Extract'].lower() in ['t','true']:
                self.opt['Extract'] = True
            else:
                self.opt['Extract'] = extract
            if not self.opt['Extract'] and not self.opt['Compile']:
                return False
            
            ### Folder for Extracted Data ###
            _stage = 'OutDir'
            if not os.path.exists(self.info['OutDir']):
                os.mkdir(self.info['OutDir'])
                self.log.printLog('#DIR','%s made for extracted output.' % self.info['OutDir'])
            for file in ['Name','OccRes','DatOut','FTOut']:
                if self.info[file] not in ['','None']:
                    self.info[file] = self.info['OutDir'] + self.info[file]
                else:
                    self.info[file] = ''    # For ease of checking later

            ### Main Compilation Setup ###
            _stage = 'Compilation'
            if self.info['Name'] in ['','None']:
                self.info['Delimit'] = rje.getDelimit(self.cmd_list,default=',')
            else:
                self.info['Delimit'] = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['Name']))
            if self.opt['Compile']:
                rje.backup(self,self.info['Name'],unlink=True)
            self.opt['PRESTO Calculations'] = False     # whether to perform occurrence-specific PRESTO calculations
            special_stats = []

            ## Basic Compilation Stats ##
            headers = ['Dataset','SeqNum','TotalAA','Rank','Score','Pattern','Occ','IC','Norm','Sim']
            if self.info['SlimVersion'] != '1.3':
                headers = ['Dataset','SeqNum','FullMST','TotalAA','Rank','Score','Pattern','Occ','IC','Norm','Sim','MST','UHS','UP','Expect','Prob']
            ## Expect ##
            if self.opt['Expect']:
                for h in ['Expect','RScore','Prob']:
                    if h not in headers:
                        headers.append(h)
                if self.opt['AdvProb']:
                    headers.append('AdvProb')
                ## NewRanks ##
                if self.stat['ReRank'] > 0:
                    headers.append('OldRank')
            ## Z-Score ##
            if self.opt['ZScore']:
                headers.append('ZScore')
            ## PrestoIC ##
            if self.opt['MotifIC']:
                headers.append('OldIC')
            ## Additional motif stats ##
            special_stats = ['SA','Hyd','Cons','Fold','IUP']
            for mstat in ['AbsChg','NetChg','BalChg','AILMV','Aromatic','Phos']:
                if self.opt[mstat]:
                    headers.append(mstat)
                    if self.opt['SlimChg'] and mstat[-3:] == 'Chg':
                        special_stats.append(mstat)
                        self.opt['Slim%s' % mstat] = True
            ## Special calculations ##
            for special in special_stats[0:]:
                if self.opt['Slim%s' % special]:
                    self.opt['PRESTO Calculations'] = True
                    headers.append('%s_mean' % special)
                    if self.stat['Percentile'] > 0:
                        pc = 100
                        while pc > 0:
                            headers.append('%s_pc%d' % (special,pc))
                            pc -= self.stat['Percentile']
                else:
                    special_stats.remove(special)
            if self.info['PicksID']:
                headers.append('PicksID')

            ### Setup Custom Scores ###
            ## Backwards compatibility ##
            if self.list['Custom'] and not self.dict['NewScore'].has_key('Custom'):
                self.dict['NewScore']['Custom'] = string.join(self.list['Custom'],'*')
                self.list['NewScore'].append('Custom')
            ## General Setup ##
            (self.list['Headers'],self.list['NewScore'],self.dict['NewScore']) = rje_scoring.setupCustomScores(self,headers[0:],self.list['NewScore'],self.dict['NewScore'])
            ## CustomOcc dictionary ##
            self.dict['CustomOcc'] = {}         # Whether the Custom Score needs occurrence stats
            for new in self.list['NewScore']:   # self.dict['NewScore'] keys() in order they were read in
                self.dict['CustomOcc'][new] = True  # Need occurrence statistics
                if rje.formula(self,self.dict['NewScore'][new],varlist=self.list['PStats'],check=False,calculate=False):
                    self.dict['CustomOcc'][new] = False
                if not self.dict['CustomOcc'][new]:
                    self.list['PStats'].append(new)

            ### Output Compilation Headers ###
            self.compileOut(self.info['Name'],self.list['Headers'],self.info['Delimit'])

            ### Update global attributes ###                
            self.list['Special Stats'] = special_stats
            #X#self.deBug(self.list['Headers'])
            
            ### Setup RankStat and StatFilter based on headers ###
            if self.opt['ZFilter']:
                self.opt['RankFilter'] = True
            ## RankStat ##
            if self.info['RankStat'] not in self.list['Headers']:
                for h in self.list['Headers']:
                    if h.lower() == self.info['RankStat'].lower():
                        self.info['RankStat'] = h
                        break
            if self.info['RankStat'] not in self.list['Headers']:
                self.log.errorLog('RankStat "%s" not recognised. Will use (R)Score.' %  self.info['RankStat'] == 'RScore',printerror=False)
                self.info['RankStat'] == 'RScore'
            if self.info['RankStat'] == 'RScore' and not self.opt['Expect']:
                self.info['RankStat'] == 'Score'
            ## Make StatFilter dictionary from list ##
            self.dict['StatFilter'] = rje_scoring.setupStatFilter(self,self.list['Headers'],self.list['StatFilter'])
            self.opt['FilterOcc'] = False
            for stat in self.dict['StatFilter'].keys():
                if stat not in self.list['PStats']:
                    self.opt['FilterOcc'] = True

            ### GOPHER Alignment Generation and Alteration of Attributes for prestoStats() and proteinAln() ###
            if self.opt['Gopher']:
                if string.split(self.info['AlnDir'],os.sep)[-2] == 'ALN':
                    self.info['GopherDir'] = string.join(string.split(self.info['AlnDir'],os.sep)[:-2],os.sep)
                    #X#self.deBug('GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                else:
                    self.info['GopherDir'] = '%s' % self.info['AlnDir']
                    self.info['AlnDir'] = rje.makePath(self.info['AlnDir'] + 'ALN')
                    #X#self.deBug('GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                #X#self.deBug('GopherDir = %s' % self.info['GopherDir'])
                if not os.path.exists(self.info['GopherDir']):
                    os.mkdir(self.info['GopherDir'])
                if not os.path.exists(self.info['AlnDir']):
                    os.mkdir(self.info['AlnDir'])
                self.info['AlnExt'] = 'orthaln.fas'

            ### PRESTO and SeqList Objects ###
            _stage = 'Objects'
            if self.opt['SlimCons'] or self.opt['ProteinAln']:
                cons_cmd = ['usealn=T']
                if self.opt['Gopher']:
                    cons_cmd.append('gopher=T')
            else:
                cons_cmd = ['usealn=F']
            cons_cmd += ['alndir=%s' % self.info['AlnDir'],'alnext=%s' % self.info['AlnExt']]
            cons_cmd += ['resfile=%s' % self.info['OccRes'],'newlog=F']
            if self.opt['LenCorrect']:
                cons_cmd += ['motific=T']
            #X#self.deBug(cons_cmd)
            self.log.printLog('#PASS',string.join(self.cmd_list+cons_cmd))
            slim_presto = presto.Presto(self.log,self.cmd_list+cons_cmd)
            self.log.printLog('#PRESTO',string.join(slim_presto.cmd_list))
            slim_presto.posMatrix()    
            slim_presto.opt['ConsAmb'] = True   #!# Cannot use this option #!#
            slim_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F','extract=None'])
            if self.opt['SlimCons'] and not slim_presto.opt['UseAln']:
                self.log.printLog('#ERR','Conservation score slimcons=True but PRESTO option usealn=False. Cancelled SlimCons.')
                self.log.printLog('#CMD','Use PRESTO options "usealn=T alndir=PATH alnext=X alngap=T/F" for conservation calculations.')
                self.opt['SlimCons'] = False
            self.obj['SlimPRESTO'] = slim_presto    # PRESTO object for storing motifs and making calculations
            self.obj['SlimSeq'] = slim_seq          # SeqList object for handling search dataset sequences

            ### Setup File of Data per occurrence ###
            _stage = 'OccRes'
            self.list['PRESTO Stats'] = []
            if slim_presto.opt['Peptides'] and self.info['OccRes'].lower() in ['','none']:
                self.info['OccRes'] = rje.baseFile(self.info['Name']) + '.occ.' + rje.delimitExt(self.info['Delimit'])
            if self.info['OccRes'].lower() not in ['','none']:      
                rje.backup(self,self.info['OccRes'],unlink=True)
                self.opt['PRESTO Calculations'] = True
                self.list['PRESTO Stats'] = ['Motif','Len','Hit','Start_Pos','End_Pos','MatchSeq','MotifCons','HomNum','GlobID','LocID','SA','Hyd']
                if self.opt['SlimFold']:
                    self.list['PRESTO Stats'].append('Fold')
                if self.opt['SlimIUP']:
                    self.list['PRESTO Stats'].append('IUP')
                if self.opt['SlimChg']:
                    self.list['PRESTO Stats'] += ['AbsChg','NetChg','BalChg']
                for key in rje.sortKeys(slim_presto.dict['ConsSpecLists']):
                    for h in ['Cons','Hom','GlobID','LocID']:
                        self.list['PRESTO Stats'].append('%s_%s' % (key,h))
                if slim_presto.opt['Peptides']:
                    self.list['PRESTO Stats'] += ['PepSeq','PepDesign']
                    rje.backup(self,rje.baseFile(self.info['Name'])+'.peptides.fas')
                self.compileOut(self.info['OccRes'],self.list['PRESTO Stats'],self.info['Delimit'])               
                
            return True
        except:
            self.log.errorLog('Error in SlimPicker.picks Setup(%s)' % _stage,printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def compileOut(self,filename,headers,delimit,datadict={}):   ### Outputs a single line of compiled data to file
        '''
        Outputs a single line of compiled data to file.
        >> filename:str = Name of output file
        >> headers:list of field headers
        >> delimit:str = text delimiter
        >> datadict:dictionary of {header:str} = data to output. If none, will output headers themselves.
        '''
        rje.delimitedFileOutput(self,filename,headers,delimit,datadict)
        return
#########################################################################################################################
    def motifOccDict(self,inputseq=None,extract_prots=[]):  ### Creates self.dict['Motif PNames'] and self.dict['Motif Occ']
        '''
        Creates self.dict['Motif PNames'] and self.dict['Motif Occ'].
        >> inputseq:str = Original input sequence file [None]
        >> extract_prots = List of proteins to extract details for *in input data order*
        '''
        try:
            ### Setup ###
            self.dict['Motif PNames'] = {}
            self.dict['Motif Occ'] = {}
            if not inputseq:
                return
            motif_occ = {}      # Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
            motif_pnames = {}   # Dictionary of {Motif:[pnames]}
            slim_seq = self.obj['SlimSeq']
            slim_presto = self.obj['SlimPRESTO']
            fullset = self.dict['FullSet']
            slim_prot = self.dict['SlimProt']
            
            ### Generate seq_dict of name:Sequence object mappings ###
            slim_seq.seq = []
            if self.list['ProtList'] and self.opt['Strict']:
                slim_seq.list['GoodSeq'] = self.list['ProtList']
                extract_prots = self.list['ProtList'][0:]
                self.log.printLog('#FILT','Strict filtering of proteins not in ProtList.')
            elif len(slim_prot.keys()) > len(extract_prots):
                slim_seq.list['GoodSeq'] = extract_prots[0:]    # Filter by extract_prots
                self.log.printLog('#FILT','Filtering proteins not containing 1+ motifs of interest')
            else:
                slim_seq.list['GoodSeq'] = []
            slim_seq.list['GoodAcc'] = slim_seq.list['GoodSeq']    # Filter by extract_prots
            #X#self.deBug('Loading from: %s\nFiltering to keep %d seqs: %s' % (inputseq,len(slim_seq.list['GoodSeq']),slim_seq.list['GoodSeq']))
            slim_seq.loadSeqs(seqfile=inputseq)
            if slim_seq.list['GoodSeq']:
                self.log.printLog('#FILT','%s sequences retained after filter.' % rje.integerString(slim_seq.seqNum()))
            seq_dict = slim_seq.seqNameDic('Max')    # Dictionary of name:Sequence object

            #!# Added a fudge for Pushker's re-coded SLiMDisc #!#
            seq_dict_ncbi = slim_seq.seqNameDic(key='NCBI')     # Dictionary of name:Sequence object
            #!# Add a fudge for annoying reformatted UniProt entries #!#
            seq_dict_ft = slim_seq.seqNameDic(key='UniProt')    # Dictionary of name:Sequence object
            for pname in slim_prot.values():
                #X#if pname not in extract_prots:  # No need to extract this sequence
                #X#    self.deBug('%s filtered out as unnecessary!' % pname)
                if pname in extract_prots and not seq_dict.has_key(pname):
                    if seq_dict_ft.has_key(pname):
                        seq_dict[pname] = seq_dict_ft[pname]
                    elif seq_dict_ncbi.has_key(pname):                    
                        seq_dict[pname] = seq_dict_ncbi[pname]
                    else:
                        self.log.errorLog('Could not map protein "%s" onto sequences from "%s"!' % (pname,inputseq),printerror=False)
                        self.deBug(rje.sortKeys(seq_dict))
                        self.deBug(rje.sortKeys(seq_dict_ft))
                        self.deBug(rje.sortKeys(seq_dict_ncbi))
                        #X#self.deBug(slim_seq.list['GoodSeq'])
                        seq_dict[pname] = None
            #!# End of fudge #!#

            ### Populate motif_occ with occurrences and motif_pnames with links ###
            for motif in slim_presto.list['Motifs']:
                pattern = motif.info['Name']
                if not self.dict['BasicStats'].has_key(pattern):
                    self.log.errorLog('Pattern %s missing from basicstats keys!' % pattern, printerror=False)
                    self.dict['BasicStats'][pattern] = {'Pattern':pattern}
                    for key in fullset.keys():
                        self.dict['BasicStats'][pattern][key] = fullset[key]
                    self.dict['BasicStats'][pattern]['Rank'] = '-1'
                    self.dict['BasicStats'][pattern]['Score'] = '-1'
                    self.dict['BasicStats'][pattern]['Occ'] = '1'
                    self.dict['BasicStats'][pattern]['IC'] = '0'
                    self.dict['BasicStats'][pattern]['Norm'] = '1'
                    self.dict['BasicStats'][pattern]['Sim'] = '0'
                motif_pnames[motif] = []
                for occ in self.dict['SlimMotif'][pattern]:
                    (id,pos) = rje.matchExp('^(\d+):(\d+)$',occ)
                    pname = self.dict['SlimProt'][id]
                    if pname not in seq_dict.keys():    # PName filtered (probably by strict=T)
                        continue
                    motif_pnames[motif].append(pname)
                    seq = seq_dict[pname]   # ID -> shortName() -> Sequence Object
                    if not seq: #!# Get to bottom of this #!#
                        continue
                    if not motif_occ.has_key(seq):
                        motif_occ[seq] = {}
                    if not motif_occ[seq].has_key(motif):
                        motif_occ[seq][motif] = []
                    motif_occ[seq][motif].append(string.atoi(pos)-1)

            ### Finish ###
            self.dict['Motif PNames'] = motif_pnames
            self.dict['Motif Occ'] = motif_occ
            self.dict['SeqDict'] = seq_dict

        except:
            self.log.errorLog('Error in SlimPicker.motifOccDict')
            raise
#########################################################################################################################
    def reduceMotifOcc(self,motiflist):     ### Reduce self.dict['Motif Occ'] to match reduced motiflist
        '''
        Reduce self.dict['Motif Occ'] to match reduced motiflist.
        >> motiflist:list of retained motifs following filters etc.
        '''
        try:
            ### Remove filtered motifs and seqs without motifs left ###
            self.deBug(self.dict['Motif Occ'])
            self.deBug(motiflist)
            for seq in self.dict['Motif Occ'].keys()[0:]:
                for motif in self.dict['Motif Occ'][seq].keys()[0:]:
                    if motif.info['Name'] not in motiflist:
                        self.dict['Motif Occ'][seq].pop(motif)
                if not self.dict['Motif Occ'][seq]:
                    self.dict['Motif Occ'].pop(seq)
            self.deBug(self.dict['Motif Occ'])
        except:
            self.log.errorLog('Error in SlimPicker.reduceMotifOcc()')
#########################################################################################################################
    def slimUniProt(self,dataset,motiflist):    ### Creates self.dict['SlimFT'] and self.list['UniExtract']
        '''
        Creates self.dict['SlimFT'] and self.list['UniExtract'].
        >> dataset:str = current dataset being processed.
        >> motiflist:list of patterns to output
        '''
        try:
            ### Setup ###
            slim_ft = self.dict['SlimFT']
            uni_extract = self.list['UniExtract']
            slim_presto = self.obj['SlimPRESTO']
            slim_seq = self.obj['SlimSeq']
            
            ### Slim FT ###
            if self.info['FTOut'] and self.opt['Extract']:
                for motif in slim_presto.list['Motifs']:
                    pattern = motif.info['Name']
                    if pattern not in motiflist:
                        continue
                    for occ in self.dict['SlimMotif'][pattern]:
                        (id,pos) = rje.matchExp('^(\d+):(\d+)$',occ)
                        pname = self.dict['SlimProt'][id]
                        if self.dict['SeqDict'].has_key(pname):
                            seq = self.dict['SeqDict'][pname]
                        else:
                            continue    # PName filtered (probably by strict=T)
                    
                        ## slim_ft ##
                        acc = seq.info['AccNum']
                        try:
                            my_ft = {'Type':'SLIM','Start':string.atoi(pos),'End':(string.atoi(pos)-1+motif.stat['Length']+string.count(pattern,'.')),
                                     'Desc':'%s from %s SLiMDisc (%s).' % (pattern,dataset,self.dict['BasicStats'][pattern]['Rank'])}
                        except:
                            my_ft = {'Type':'SLIM','Start':string.atoi(pos),'End':(string.atoi(pos)-1+motif.stat['Length']+string.count(pattern,'.')),
                                     'Desc':'%s from %s SLiMDisc (???).' % (pattern,dataset)}
                        if not slim_ft.has_key(acc):
                            slim_ft[acc] = []   # Dictionary of {AccNum:[Extra features (SLiMs) to add to table}
                        slim_ft[acc].append(my_ft)
                        if rje.matchExp('^(\S+)\-(\d+)',acc):   # Splice variant
                            svacc = rje.matchExp('^(\S+)\-(\d+)',acc)[0]
                            try:
                                sv_ft = {'Type':'SLIM','Start':string.atoi(pos),'End':(string.atoi(pos)-1+motif.stat['Length']+string.count(pattern,'.')),
                                         'Desc':'%s from %s SLiMDisc (%s). (Splicevar %s)' % (pattern,dataset,self.dict['BasicStats'][pattern]['Rank'],acc)}
                            except:
                                sv_ft = {'Type':'SLIM','Start':string.atoi(pos),'End':(string.atoi(pos)-1+motif.stat['Length']+string.count(pattern,'.')),
                                         'Desc':'%s from %s SLiMDisc (???). (Splicevar %s)' % (pattern,dataset,acc)}
                            if not slim_ft.has_key(svacc):
                                slim_ft[svacc] = []   # Dictionary of {AccNum:[Extra features (SLiMs) to add to table}
                            slim_ft[svacc].append(sv_ft)

            ### Update uni_extract with list of accession numbers to extract from UniProt ###
            if (self.info['DatOut'] or self.info['FTOut']) and self.opt['Extract']:
                for seq in slim_seq.seq:
                    if seq.info['AccNum'] not in uni_extract:
                        uni_extract.append(seq.info['AccNum'])

            ### Finish ###
            self.dict['SlimFT'] = slim_ft
            self.list['UniExtract'] = uni_extract

        except:
            self.log.errorLog('Error in SlimPicker.slimUniProt()')
            raise
#########################################################################################################################
    def expectation(self,slim_presto,motif_pnames,fasta):    ### Calculates expected support and probability
        '''
        Calculates expected support and probability. This part of the program has three functions:
            1. Calculate the expected *support* of a motif from the input dataset
            2. Use this with the observed support to calculate the probability of seeing at least that much support
            3. Create a dictionary of probabilities that each motif will occur in each protein (for occres)
        >> slim_presto:PRESTO object containing Motif objects (pattern = motif.info['Name'])
        >> motif_pnames:dictionary of {Motif:[pnames]}
        >> fasta:str = Name of TEIRESIAS input file. Link sequences via shortName()
        << prob_occ:Dictionary of {pname:{Motif:prob}}
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            ## Missing fasta ##
            if not fasta:
                for motif in slim_presto.list['Motifs']:
                    pattern = motif.info['Name']
                    self.dict['BasicStats'][pattern]['Expect'] = -1
                    self.dict['BasicStats'][pattern]['RScore'] = 0
                    slimdisc = self.dict['BasicStats'][pattern]['Dataset']
                self.log.errorLog('No input files found for %s: Cannot calculate expectations.' % slimdisc,printerror=False)
                return {}
            ## Input sequences ##
            #X#seq_cmd = self.cmd_list + ['seqin=%s' % fasta,'accnr=F','seqnr=F','autofilter=F','autoload=T','memsaver=F','i=-1']
            #X#self.deBug(seq_cmd)
            inseq = rje_seq.SeqList(self.log,cmd_list=self.obj['SlimSeq'].cmd_list)
            inseq.loadSeqs(seqfile=fasta)
            if inseq.seqNum() < 1:
                self.log.errorLog('No sequences loaded from %s! No expectation calculations possible.' % fasta,printerror=False)
                return {}
            ## Dictionaries ##
            expect = {}         # Dictionary of {Motif:expected support}
            poisson_prob = {}   # Dictionary of {Motif:probability of observed support}
            prob_occ = {}       # Dictionary of {pname:{Motif:prob}}
            seq_frag = {}       # Dictionary of {pname:[Inseq Sequence objects]}
            prot_aafreq = {}    # Dictionary of {pname:{aa:freq}}
            sx = 0
            for seq in inseq.seq:
                sx += 1
                self.log.printLog('\r#AA','Amino acid frequencies: %s fragments' % rje.integerString(sx),log=False,newline=False)
                # Setup
                pname = seq.shortName()
                if not seq_frag.has_key(pname):
                    seq_frag[pname] = []
                    prob_occ[pname] = {}
                    prot_aafreq[pname] = {'Total':0.0}
                    for aa in rje_seq.alph_protx:
                        prot_aafreq[pname][aa] = 0.0
                # Append details
                seq_frag[pname].append(seq)
                for aa in rje_seq.alph_protx:
                    x = string.count(seq.info['Sequence'].upper(),aa)
                    prot_aafreq[pname][aa] += x
                    prot_aafreq[pname]['Total'] += x
                
            ### AA Frequencies ###
            _stage = 'AA Freq'
            px = 0
            totx = 0
            for pname in prot_aafreq.keys():
                px += 1
                self.log.printLog('\r#AA','Amino acid frequencies: %s fragments; %s proteins.' % (rje.integerString(sx),rje.integerString(px)),log=False,newline=False)
                if prot_aafreq[pname]['Total']:
                    totx += prot_aafreq[pname]['Total']
                    for aa in rje_seq.alph_protx:
                        prot_aafreq[pname][aa] /= prot_aafreq[pname]['Total']
                else:
                    self.log.printLog('\n#ERR','No amino acids read in for %s!' % pname)
                    prob_occ.pop(pname) # Do not calculate if no aa!
            self.log.printLog('\r#AA','Amino acid frequencies: %s fragments; %s proteins; %s AA.' % (rje.integerString(sx),rje.integerString(px),rje.integerString(totx)))
           
            ### Expectations ###
            _stage = 'Expectation'
            # motif.info['Sequence'] = Presto format motifs are strings of elements separated by '-', where each element is:
            # > a single AA letter
            # > a choice of letters in the form [ABC]
            # > a choice of combinations in the form (AB|CD)    (*Should have none of these!*)
            # prob_occ = {}       # Dictionary of {pname:{Motif:prob}}
            calcx = len(prob_occ) * slim_presto.motifNum()
            cx = 0
            for motif in slim_presto.list['Motifs']:
                pattern = motif.info['Name']
                self.dict['BasicStats'][pattern]['Expect'] = 0
                for pname in prob_occ.keys():
                    ## Setup ##
                    cx += 100.0
                    self.log.printLog('\r#EXP','Expected support calculation: %.1f%%' % (cx/calcx),log=False,newline=False)
                    elements = motif.list['PRESTO']
                    motiflen = len(elements)    # Length of motif
                    faa = prot_aafreq[pname]    # AA frequencies for this protein
                    faa['X'] = 1.0  # Wildcard does not affect probability
                    positions = faa['Total'] - (len(seq_frag[pname]) * (motiflen - 1))  # No. of possible motif starting positions
                    ## Calculate probability per position ##
                    p = 1.0
                    for pos in elements:     # Position in motif
                        pos_p = 0.0     # Probability of this position in the motif
                        for aa in pos:  # Take each aa in this position
                            if faa.has_key(aa):     # Ignores [] and unconventional aa symbols
                                pos_p += faa[aa]    # Sum probability over possible motifs
                        p *= pos_p     # Overall probability is product of different positions
                    ## Calculate expected frequency of motif in sequence ##
                    expect = p * positions
                    ## Calculate probability of occurrence in this sequence ##
                    prob_occ[pname][motif] = 1 - math.exp(-expect)
                    ## Update expected support for motif across all sequences ##
                    self.dict['BasicStats'][pattern]['Expect'] += prob_occ[pname][motif]
            self.log.printLog('\r#EXP','Expected support calculation: Complete.')

            ### New Information Content ###
            _stage = 'New IC'   
            if self.opt['MotifIC']:     #!# Move this - will not be called if expect=F #!#
                for motif in slim_presto.list['Motifs']:
                    self.dict['BasicStats'][motif.info['Name']]['IC'] = '%.3f' % motif.stat['IC']
            if self.opt['LenCorrect']:
                for motif in slim_presto.list['Motifs']:
                    self.dict['BasicStats'][motif.info['Name']]['LenIC'] = motif.stat['IC']

            ### Calculate probability of Observed support ###
            _stage = 'Observed Prob'
            mx = 0
            for motif in slim_presto.list['Motifs']:
                mx += 100.0
                self.log.printLog('\r#OBS','Probability of observed support calculation: %.1f%%' % (mx/slim_presto.motifNum()),log=False,newline=False)
                pattern = motif.info['Name']
                observed = string.atoi(self.dict['BasicStats'][pattern]['Occ'])
                expected = self.dict['BasicStats'][pattern]['Expect']
                self.dict['BasicStats'][pattern]['Expect'] = rje_motif.expectString(expected)
                self.dict['BasicStats'][pattern]['RScore'] = string.atof(self.dict['BasicStats'][pattern]['IC']) * string.atof(self.dict['BasicStats'][pattern]['Norm']) * observed / expected
                if self.opt['LenCorrect']:
                    self.dict['BasicStats'][pattern]['RScore'] /= (20 ** self.dict['BasicStats'][pattern]['LenIC'])
                prob = 0
                for x in range(0,observed):
                    try:    #!# Fudge for OverflowError: long int too large to convert to float
                        prob += (math.exp(-expected) * pow(expected,x) / rje.factorial(x))
                    except:
                        break
                self.dict['BasicStats'][pattern]['Prob'] = rje_motif.expectString(1 - prob)
            self.log.printLog('\r#OBS','Probability of observed support calculation: Complete.')

            ### Advanced Probability Calculation ###
            _stage = 'Advanced Prob'        
            if not self.opt['AdvProb']:
                return prob_occ
            mx = 0
            if len(prob_occ) > self.stat['AdvMax']:
                for motif in slim_presto.list['Motifs']:
                    mx += 100.0
                    self.log.printLog('\r#OBS','Advanced Probability of observed support calculation: %.1f%%' % (mx/slim_presto.motifNum()),log=False,newline=False)
                    pattern = motif.info['Name']
                    plist = motif_pnames[motif]
                    observed = len(plist)
                    ## Calculate probability of first occurrence for each prot in plist ##
                    first_occ = {}
                    total_occ = 0.0
                    for pname in plist:
                        total_occ += prob_occ[pname][motif]
                    ## Calculate probabilities of reduced observed support in other proteins ##
                    advprob = 0.0
                    for firstname in plist:
                        first_occ = prob_occ[firstname][motif] / total_occ  # This is the probability of firstname being the first occurrence
                        redsup = 0.0
                        for pname in prob_occ.keys():
                            if pname != firstname:
                                redsup += prob_occ[pname][motif]
                        prob = 0    # The probability of the other occurrences if firstname is the first occurrence
                        for x in range(0,observed-1):
                            prob += (math.exp(-redsup) * pow(redsup,x) / rje.factorial(x))
                        advprob += ((1 - prob) * first_occ)
                    self.dict['BasicStats'][pattern]['AdvProb'] = rje_motif.expectString(advprob)
                self.log.printLog('\r#OBS','Advanced Probability of observed support calculation: Complete.')
                return prob_occ

            ### Advanced Probability Calculation 2 ###
            _stage = 'Advanced Prob 2'        
            mx = 0
            for motif in slim_presto.list['Motifs']:
                mx += 100.0
                self.log.printLog('\r#OBS','Advanced Probability of observed support calculation: %.1f%%' % (mx/slim_presto.motifNum()),log=False,newline=False)
                pattern = motif.info['Name']
                plist = motif_pnames[motif]
                observed = len(plist)
                ## Calculate probability of first occurrence for each prot in plist ##
                first_occ = {}
                total_occ = 0.0
                for pname in plist:
                    total_occ += prob_occ[pname][motif]
                ## Calculate probabilities of reduced observed support in other proteins ##
                advprob = 0.0
                for firstname in plist:
                    first_occ = prob_occ[firstname][motif] / total_occ  # This is the probability of firstname being the first occurrence
                    prob = 0.0
                    binlist = [0] * (len(prob_occ) - 1)
                    while sum(binlist) < (len(prob_occ) - 1):
                        binlist = rje.binaryCount(binlist)
                        if sum(binlist) >= (observed - 1):      # This counts!
                            this_prob = 1.0
                            elist = rje.sortKeys(prob_occ)
                            elist.remove(firstname)
                            for i in range(len(binlist)):
                                if binlist[i] > 0:
                                     this_prob *= prob_occ[elist[i]][motif]
                            prob += this_prob
                    advprob += (prob * first_occ)
                self.dict['BasicStats'][pattern]['AdvProb'] = rje_motif.expectString(advprob)
            self.log.printLog('\r#OBS','Advanced Probability of observed support calculation: Complete.')

            ### Finish ###
            return prob_occ
        
        except:
            self.log.errorLog('Error in SlimPicker.expectation(%s)' % _stage,printerror=True,quitchoice=True)
            return prob_occ
#########################################################################################################################
    def customScore(self,motiflist,new):  ### Calculates custom score for all motifs in dataset
        '''
        Calculates custom score for all motifs in dataset.
        >> motiflist:list of patterns (in rank order)
        >> new:str = Key for self.dict['NewScore']
        '''
        try:
            ### Setup ###
            if not motiflist or not self.dict['NewScore'].has_key(new):
                return False

            ### Score ###
            for pattern in motiflist:
                #X#self.deBug(self.dict['BasicStats'][pattern])
                value = rje.formula(self,formula=self.dict['NewScore'][new],data=self.dict['BasicStats'][pattern])
                self.dict['BasicStats'][pattern][new] = value

            ### Finish ###
            return True                    
        except:
            self.log.errorLog('Error in SlimPicker.zScores()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def zScores(self,motiflist):  ### Calculates Z-Scores for all motifs in dataset
        '''
        Calculates Z-Scores for all motifs in dataset.
        >> motiflist:list of patterns (in rank order)
        '''
        try:
            ### Setup ###
            if not self.opt['Compile'] or not self.opt['ZScore'] or not motiflist:
                return False

            ### Mean ###
            meanscore = 0.0
            for pattern in motiflist:
                try:
                    meanscore += self.dict['BasicStats'][pattern][self.info['RankStat']]
                except:
                    self.dict['BasicStats'][pattern][self.info['RankStat']] = string.atof(self.dict['BasicStats'][pattern][self.info['RankStat']])
                    meanscore += self.dict['BasicStats'][pattern][self.info['RankStat']]
            meanscore /= len(motiflist)

            ### SD ###
            sdscore = 0.0
            for pattern in motiflist:
                dif = (self.dict['BasicStats'][pattern][self.info['RankStat']] - meanscore)
                sdscore += (dif * dif)
            sdscore = math.sqrt(sdscore/len(motiflist))

            ### Z ###                    
            for pattern in motiflist:
                if not sdscore:
                    self.dict['BasicStats'][pattern]['ZScore'] = 0.0
                else:
                    self.dict['BasicStats'][pattern]['ZScore'] = (self.dict['BasicStats'][pattern][self.info['RankStat']] - meanscore) / sdscore

            ### Finish ###
            return True                    
        except:
            self.log.errorLog('Error in SlimPicker.zScores()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def reRank(self,motiflist):     ### ReRanks according to the new RScore and imposes second rank cut-off
        '''
        ReRanks according to the new RScore and imposes second rank cut-off.
        >> motiflist:list of patterns (in rank order)
        << motiflist:list of patterns, re-ranked and filtered.
        '''
        try:
            ### Setup ###
            if not self.opt['Expect'] or self.stat['ReRank'] < 1:
                return motiflist

            ### NewRanks ###
            scores = []
            for pattern in motiflist:
                scores.append(self.dict['BasicStats'][pattern][self.info['RankStat']])
            newranks = rje.rankList(scorelist=scores,rev=True,absolute=True,lowest=True)
            for pattern in motiflist:
                self.dict['BasicStats'][pattern]['Rank'] = newranks.pop(0)
            if newranks:
                self.log.errorLog('New Ranks remain after all patterns processed!',printerror=False)
                raise ValueError

            ### Rank cut-off ###
            rankdict = {}
            for pattern in motiflist:
                if self.dict['BasicStats'][pattern]['Rank'] <= self.stat['ReRank']:   # Proceed!
                    r = self.dict['BasicStats'][pattern]['Rank']
                    if rankdict.has_key(r):
                        rankdict[r].append(pattern)
                    else:
                        rankdict[r] = [pattern]
                    self.dict['BasicStats'][pattern]['Rank'] = '%d' % self.dict['BasicStats'][pattern]['Rank']
            newlist = []
            for r in rje.sortKeys(rankdict):
                newlist += rankdict[r]

            ### Finish ###
            self.log.printLog('#RANK','%s motifs re-ranked; %s retained.' % (rje.integerString(len(motiflist)),rje.integerString(len(newlist))))
            return newlist
        except:
            self.log.errorLog('Error in SlimPicker.reRank()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def statFilter(self,motiflist,type='All'): ### Filters motifs according to self.list['StatFilter']. 
        '''
        Filters motifs according to self.list['StatFilter']. 
        >> motiflist:list of patterns 
        >> type:str = P/O Stats or All for both (default)
        << motiflist:list of filtered patterns.
        '''
        try:
            ### New Filtering Using rje_scoring ###
            statfilter = {}
            for stat in self.dict['StatFilter'].keys():
                if type not in ['O','P'] or (type == 'P' and stat in self.list['PStats']) or (type == 'O' and stat not in self.list['PStats']):
                    statfilter[stat] = self.dict['StatFilter'][stat]
            #X#print type, statfilter
            if not statfilter:
                return motiflist
            mx = len(motiflist)
            ### Filter ###
            self.dict['BasicStats'] = rje_scoring.statFilter(self,self.dict['BasicStats'],statfilter)
            for pattern in motiflist[0:]:
                if not self.dict['BasicStats'].has_key(pattern):
                    motiflist.remove(pattern)
            ### Finish ###
            self.log.printLog('#FILTER','%s motifs filtered by %d %s stats; %s retained.' % (rje.integerString(mx),len(statfilter),type,rje.integerString(len(motiflist))))
            return motiflist

            
            ### Setup ###
            statfilter = self.dict['StatFilter'].keys()
            if type in ['O','P']:
                for stat in statfilter[0:]:
                    if (type == 'P' and stat not in self.list['PStats']) or (type == 'O' and stat in self.list['PStats']):
                        statfilter.remove(stat)
            #X#print type, self.list['PStats'], statfilter
            if not statfilter:
                return motiflist
            mx = len(motiflist)

            ### Filter patterns ###
            for stat in statfilter:     # {stat:(op,cutoff,numcut)}
                (op,strcut,numcut) = self.dict['StatFilter'][stat]
                for pattern in motiflist[0:]:
                    ## Check for stat ##
                    if not self.dict['BasicStats'][pattern].has_key(stat):
                        self.log.errorLog('Pattern "%s" missing stat "%s"!' % (pattern,stat),printerror=False)
                        continue
                    value = self.dict['BasicStats'][pattern][stat]
                    ## Numeric? ##
                    numeric = None
                    if numcut:
                        try:
                            numeric = float(self.dict['BasicStats'][pattern][stat])
                        except:
                            numeric = None
                    ## Evaluate and Filter ##
                    if numcut and numeric:
                        (value,cutoff) = (numeric,numcut)
                    else:
                        cutoff = strcut
                    #X#print stat, cutoff, op, pattern, value, numeric
                    try:
                        if op == '==' and value == cutoff:
                            motiflist.remove(pattern)
                        elif op == '!=' and value != cutoff:
                            motiflist.remove(pattern)
                        elif op == '>=' and value >= cutoff:
                            motiflist.remove(pattern)
                        elif op == '>' and value > cutoff:
                            motiflist.remove(pattern)
                        elif op == '<=' and value <= cutoff:
                            motiflist.remove(pattern)
                        elif op == '<' and value < cutoff:
                            motiflist.remove(pattern)
                    except:
                        self.log.errorLog('Problem filtering %s by %s %s' % (pattern,stat,self.dict['StatFilter'][stat]),printerror=False)
                        break

            ### Finish ###
            self.log.printLog('#FILTER','%s motifs filtered by %d %s stats; %s retained.' % (rje.integerString(mx),len(statfilter),type,rje.integerString(len(motiflist))))
            return motiflist
        except:
            self.log.errorLog('Error in SlimPicker.reRank()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def patternStats(self,pattern):  ### Returns a dictionary of stats (AbsChg,NetChg,BalChg,AILMV)
        '''
        Returns a dictionary of stats (AbsChg,NetChg,BalChg,AILMV,Aromatic,Phos).
        >> pattern:str = Pattern - could be a motif or an occurrence of a motif.
        '''
        ### Setup ###
        charge = []
        ailmv = True
        aromatic = 0
        phos = []
        ### Calculations ###
        motseq = string.split(pattern.upper(),']')
        for part in motseq:
            ## Aromatic & Phos ##
            if part.find('F') >= 0 or part.find('W') >= 0 or part.find('Y') >= 0:
                aromatic += 1
            for p in 'STY':
                if part.find(p) >= 0 and p not in phos:
                    phos.append(p)
            ## Charge and AILMV ##
            if part[:1] == '[':  # Ambiguity
                if part in ['[KR','[RK']:   # +ve
                    charge.append(1)
                    ailmv = False
                elif part in ['[DE','[ED']:
                    charge.append(-1)
                    ailmv = False
                elif ailmv:
                    for aa in 'CDEFGHKNPQRSTWY':
                        if part.find(aa) > 0:
                            ailmv = False
                            break
            else:   # Fixed positions
                for a in part: # Motif
                    if a in ['K','R']:
                        charge.append(1)
                        ailmv = False
                    elif a in ['D','E']:
                        charge.append(-1)
                        ailmv = False
                    elif a in ['C','F','G','H','N','P','Q','S','T','W','Y']:
                        charge.append(0)
                        ailmv = False
                    else:
                        charge.append(0)
        ### Return Data ###
        if not phos:
            phos = ['X']
        phos.sort()
        return {'AbsChg':'%d' % (charge.count(1) + charge.count(-1)), 'NetChg':'%d' % sum(charge),
                'BalChg':'%d' % (sum(charge[:int(len(charge)/2)]) - sum(charge[-int(len(charge)/2):])),
                'AILMV':'%s' % ailmv,'Aromatic':'%d' % aromatic,'Phos':string.join(phos,'')}
#########################################################################################################################
    def prestoStats(self,slim_presto,motif_occ,prob_occ):  ### Populates and returns presto_seqhit list of objects storing data
        '''
        Populates and returns presto_seqhit list of objects storing data.
        >> slim_presto:Presto object containing motifs
        >> motif_occ:Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
        >> prob_occ:Dictionary of {pname:{Motif:prob}}
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            seqhit = []

            ### Generate PrestoSeqHit objects ###
            _stage = 'PrestoSeqHit'
            sx = 0
            for nextseq in motif_occ.keys():
                self.log.printLog('\r#PRESTO','PRESTO Occurence Stats: %.1f%%' % (100.0*sx/len(motif_occ)),log=False,newline=False)
                sx += 1
                ## Setup SeqHit ##
                _stage = 'Setup PrestoSeqHit'
                pname = nextseq.shortName()
                if not prob_occ.has_key(pname):
                    prob_occ[pname] = {}
                #X#newseqhit = slim_presto.newSeqHit(nextseq,varnum,expect=prob_occ[pname],log=self.log,cmd_list=self.cmd_list+['newlog=F'])
                newseqhit = slim_presto.newSeqHit(nextseq)
                newseqhit.dict['Expect'] = {}
                newseqhit.dict['SeqExp'] = {}
                for motif in slim_presto.list['Motifs']:
                    if motif in prob_occ[pname].keys():
                        newseqhit.dict['SeqExp'][motif] = {0:prob_occ[pname][motif]}
                        newseqhit.dict['Expect'][motif] = {0:prob_occ[pname][motif]}
                    else:
                        newseqhit.dict['SeqExp'][motif] = {0:-1}
                        newseqhit.dict['Expect'][motif] = {0:-1}
                seqhit.append(newseqhit)
                nextseq.deGap()
                sequence = nextseq.info['Sequence']
                seq_sa = rje_seq.surfaceAccessibility(sequence,returnlist=True)
                seq_eis = rje_seq.eisenbergHydropathy(sequence,returnlist=True)
                seq_dis = {}
                dismethod = {'IUP':'iupred','Fold':'foldindex'}
                for o in ['IUP','Fold']:
                    if self.opt['Slim%s' % o]:
                        dcmd = self.cmd_list+['disorder=%s' % dismethod[o],'sequence=%s' % sequence]
                        self.log.stat['Verbose'] -= 1
                        disorder = rje_disorder.Disorder(log=self.log,cmd_list=dcmd)
                        disorder.flatten()   # Converts to 1/0 #
                        self.log.stat['Verbose'] += 1
                        seq_dis[o] = disorder.list['ResidueDisorder'][0:]

                ## Occurrences ##
                _stage = 'PrestoSeqHit Occurrences'
                newseqhit.hit = {}
                for motif in motif_occ[nextseq].keys():
                    for r in motif_occ[nextseq][motif]:     # r is the start position of the match (0 -> L-1)
                        match = rje.matchExp('(%s)\S*$' % motif.info['Name'],sequence[r:])[0]
                        hitdict = {'Pos':(r+1), 'Variant':motif.info['Name'], 'SearchVar':motif.info['Sequence'], 'Match':match, 'ID':1.0, 'MisMatch':0}
                        newhit = newseqhit._addHit(motif,hitdict)
                        newhit.info['Name'] = motif.info['Name']
                        newhit.info['Variant'] = motif.info['Name']
                        newhit.info['Match'] = match
                        w = r - newseqhit.stat['WinSA']
                        if w < 0:
                            w = 0
                        newhit.stat['SA'] = sum(seq_sa[w:r+len(match)+newseqhit.stat['WinSA']]) / len(seq_sa[w:r+len(match)+newseqhit.stat['WinSA']])
                        w = r - newseqhit.stat['WinHyd']
                        if w < 0:
                            w = 0
                        newhit.stat['Hyd'] = sum(seq_eis[w:r+len(match)+newseqhit.stat['WinHyd']]) / len(seq_eis[w:r+len(match)+newseqhit.stat['WinHyd']])
                        ## Disorder = proportion of motif region identified as disordered ##
                        w = r - newseqhit.stat['WinDis']
                        if w < 0:
                            w = 0
                        for o in ['IUP','Fold']:
                            if self.opt['Slim%s' % o]:
                                newhit.stat[o] = sum(seq_dis[o][w:r+len(match)+newseqhit.stat['WinDis']]) / len(seq_dis[o][w:r+len(match)+newseqhit.stat['WinDis']])
                        ## Additional stats for SLiM Pickings ##
                        if self.opt['SlimChg']:
                            patstats = self.patternStats(match) # Returns a dictionary of stats (AbsChg,NetChg,BalChg,AILMV)
                            for key in patstats.keys():    
                                if self.opt[key]:
                                    newhit.stat[key] = patstats[key]
                        ## Peptide Design ##
                        pepwin = rje.modulus(newseqhit.stat['WinSize'])
                        w = r - pepwin
                        if w < 0:
                            w = 0
                        newhit.info['PepSeq'] = sequence[w:r+len(match)+newseqhit.stat['WinSize']]     # Region for calculation

                                    
                        #!# Attempted fix for stupid fucking Python #!#
                        addhit = copy.deepcopy(newhit)
                        addhit.info = copy.deepcopy(addhit.info)
                        addhit.stat = copy.deepcopy(addhit.stat)
                        newseqhit.dict['Hits'][motif][-1] = addhit
                        #!# End of stupid bug fix #!#

                ## AlnCon ##
                if self.opt['SlimCons']:        #!# Needs improving for SLiMDisc results #!#
                    newseqhit.hitAlnCon()

                ## End ##
                #X#newseqhit.verbose(0,2,'%d hits from %d motifs.' % newseqhit.hitNum(),1)
                #X#newseqhit.processSeqHit(mysql=newseqhit.opt['mySQL'])
                #X#if newseqhit.opt['UniProt']:
                #X#    newseqhit.uniProtHits()

            ### End ###
            self.log.printLog('\r#PRESTO','PRESTO Occurence Stats: Complete.')
            return seqhit
        except:            
            self.log.errorLog('Major problem with SlimPicker.prestoStats(%s)' % _stage,quitchoice=True)
            return seqhit
#########################################################################################################################
    def proteinAlignments(self,slim_presto,motif_occ):    ### Generates copies of alignments including SLIMs
        '''
        Generates copies of alignments including SLIMs.
        >> slim_presto:Presto object containing motifs  #!# Not used! #!#
        >> motif_occ:Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            alnout = rje.makePath(self.info['OutDir'] + 'ProteinAln/')
            if not os.path.exists(alnout):
                os.mkdir(alnout)

             ### Generate Alignments ###
            _stage = 'Generate Alignments'
            for seq in motif_occ.keys():
                ## Try to find alignment ##
                _stage = 'Find Alignment'
                alndir = rje.makePath(self.info['AlnDir'])
                if self.info['AlnExt'][0] != '.':
                    self.info['AlnExt'] = '.%s' % self.info['AlnExt']
                alnstart = [seq.info['AccNum'],seq.info['ID'],seq.shortName(),None]
                for file in alnstart:
                    if file:
                        file = '%s%s%s' % (alndir,file,self.info['AlnExt'])
                        if rje.checkForFile(file):   # File found
                            break
                ## Load of Make SeqList ##
                _stage = 'SeqList'
                if file:    
                    alncmd = ['seqin=None','query=%s' % seq.shortName(),'accnr=F','seqnr=F','autofilter=F','newlog=F']
                    aln = rje_seq.SeqList(self.log,cmd_list=self.cmd_list+alncmd)
                    aln.loadSeqs(seqfile=file,seqtype='Protein',aln=True,nodup=None)
                else:
                    alncmd = ['seqin=None','newlog=F']
                    aln = rje_seq.SeqList(self.log,cmd_list=self.cmd_list+alncmd)
                    aln.seq = [seq]
                    aln.obj['QuerySeq'] = seq
                    self.log.printLog('#ALN','No *%s alignment file found for %s in %s.' % (self.info['AlnExt'],seq.shortName(),alndir))
                    #X#print seq
                ## Check Query and File ##
                _stage = 'Check SeqList'
                if not aln.obj['QuerySeq'] and not aln.querySeq(query=seq.info['AccNum']):
                    self.log.printLog('#ERR','Problem finding %s in %s.' % (self.info['Name'],file))
                    continue
                qry = aln.obj['QuerySeq']
                ## New Motif Seq ##
                _stage = 'MotifSeq'
                motifseq = ['-'] * qry.seqLen()
                sequence = seq.info['Sequence']
                for motif in motif_occ[seq].keys():
                    for m in motif_occ[seq][motif]:     # r is the start position of the match (0 -> L-1)
                        match = rje.matchExp('^(%s)' % motif.info['Name'],sequence[m:])[0]
                        ## Map onto motifseq ##
                        (r,a,i) = (0,-1,0)
                        while r < qry.seqLen():
                            aa = qry.info['Sequence'][r]
                            if aa in rje_seq.alph_protx:        # Amino acid
                                a += 1
                            if a == m:  # Found start of motif 
                                while i < len(match):   # Found start of motif but not reached end
                                    if aa in rje_seq.alph_protx:
                                        if motif.list['PRESTO'][i] == 'X':
                                            if motifseq[r] == '-':
                                                motifseq[r] = 'X'
                                        else:
                                            motifseq[r] = match[i]
                                        i += 1
                                        if i == len(match):
                                            break
                                    r += 1
                                    if r == qry.seqLen():
                                        self.log.errorLog('Problem finding match "%s" for pattern "%s": overshot end of sequence' % (match,motif.info['Name']),printerror=False)
                                        break
                                    aa = qry.info['Sequence'][r]
                                break   # Go on to next occurrence
                            r += 1
                newseq = aln._addSeq('Motifs',string.join(motifseq,''))
                #X#print newseq, newseq.info
                aln.seq = [newseq] + aln.seq[:-1]
                aln.saveFasta(seqfile='%s%s.proteinaln.fas' % (alnout,qry.info['AccNum']))

        except:            
            self.log.errorLog('Major problem with SlimPicker.proteinAlignments(%s)' % _stage,quitchoice=True)
#########################################################################################################################
    def motifAlignments(self,dataset,slim_presto,motif_occ,motiflist):  ### Makes alignments of the occurrences of each motif
        '''
        Makes alignments of the occurrences of each motif.
        >> dataset:str = Dataset name
        >> slim_presto:Presto object containing motifs
        >> motif_occ:Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
        >> motiflist:List of motif strings to consider (filtered by motif and protein) *in rank order*
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            ltxt = 'Constructing motif alignments: %s motifs & %s seqs.' % (rje.integerString(slim_presto.motifNum()),
                                                                            rje.integerString(len(motif_occ)))
            self.log.printLog('#ALN',ltxt,log=False,newline=False)
            extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
            extract_seq.opt['ReplaceChar'] = False
            alndir = rje.makePath(self.info['OutDir'] + 'MotifAln/')
            if not os.path.exists(alndir):
                os.mkdir(alndir)

            ### Setup positions in sequences ###
            ## Work off motif_occ = Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
            _stage = 'SeqPos'
            mymotifs = {}   # Dictionary of (pattern:Motif}
            max_occ = {}    # Dictionary of Sequence:Max occurrences for any motif
            self.deBug(motif_occ)
            max_pos = 0      # Max motif position
            occ_seq = {'Motif':{}}    # Dictionary of {Sequence:{Motif:[sequences]}}
            for motif in slim_presto.list['Motifs']:
                mymotifs[motif.info['Name']] = motif
                occ_seq['Motif'][motif] = []
            for seq in motif_occ.keys():
                occ_seq[seq] = {}
                max_occ[seq] = 0
                for motif in slim_presto.list['Motifs']:
                    occ_seq[seq][motif] = []
                    if not motif_occ[seq].has_key(motif):
                        motif_occ[seq][motif] = []
                    if len(motif_occ[seq][motif]) > max_occ[seq]:
                        max_occ[seq] = len(motif_occ[seq][motif])
                    for pos in motif_occ[seq][motif]:
                        if pos >= max_pos:
                            max_pos = pos + 1
            self.log.printLog('\r#ALN','%s: setup complete.' % ltxt,log=False)

            ### Make sequences ###
            _stage = 'Sequences'
            for pattern in motiflist:
                self.log.printLog('\r#ALN','%s: %.f%%.' % (ltxt,100.0*motiflist.index(pattern)/len(motiflist)),log=False,newline=False)
                motif = mymotifs[pattern]
                ## Pattern info ##
                patseq = '-%s-' % rje.preZero(self.dict['BasicStats'][pattern]['Rank'],max_pos)      # Rank of motif
                overlap = len(pattern) - motif.stat['FullLength']                       # Extra length of pattern vs. longest occurrence
                #x#print pattern, motif.stat['FullLength'], overlap, int(overlap/2), int((overlap+1)/2)
                patseq += '-' * self.stat['FlankSize']
                patseq += pattern       
                patseq += '-' * (self.stat['FlankSize'] - overlap)
                occ_seq['Motif'][motif].append(patseq[0:])
                ## Occurrences ##
                for seq in motif_occ.keys():
                    for x in range(max_occ[seq]):  # Need an entry for each potential occurrence
                        if x < len(motif_occ[seq][motif]):  # Actual entry
                            r = motif_occ[seq][motif][x]
                            patseq = '-%s-' % rje.preZero(r+1,max_pos)
                            (left,right) = (r-self.stat['FlankSize'],r+motif.stat['FullLength']-1+self.stat['FlankSize'])
                            if left < 0:
                                patseq += '-' * -left + seq.info['Sequence'][:right]
                            else:
                                patseq += seq.info['Sequence'][left:right]
                            patseq += '-' * (len(occ_seq['Motif'][motif][0]) - len(patseq))
                        else:   # Add blank one
                            patseq = '-' * len(occ_seq['Motif'][motif][0])
                        occ_seq[seq][motif].append(patseq[0:])
            self.log.printLog('\r#ALN','%s: 100.0%%.' % (ltxt))

            ### Output sequence file ###
            _stage = 'Output'
            ## Motif Data ##
            name = 'SlimPicked %s Motifs' % dataset
            outlist = []
            for pattern in motiflist:
                motif = mymotifs[pattern]
                outlist.append(occ_seq['Motif'][motif][0])
            extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ## Occurrences ##
            for seq in motif_occ.keys():
                name = seq.info['Name']
                for x in range(max_occ[seq]):
                    outlist = []
                    for pattern in motiflist:
                        motif = mymotifs[pattern]
                        outlist.append(occ_seq[seq][motif][x])
                    extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ### Output ###
            extract_seq.info['Name'] = alndir + '%s.slimpicks.fas' % dataset
            extract_seq.saveFasta()

        except:            
            self.log.errorLog('Major problem with SlimPicker.motifAlignments(%s)' % _stage,quitchoice=True)
#########################################################################################################################
## End of SECTION II: SlimPicker Class                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
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
        SlimPicker(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
