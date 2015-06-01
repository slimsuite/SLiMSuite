#!/usr/local/bin/python

# GABLAM - Global Analysis of BLAST Local AlignMents 
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
Program:      GABLAM
Description:  Global Analysis of BLAST Local AlignMents
Version:      2.19.2
Last Edit:    15/03/15
Citation:     Davey, Shields & Edwards (2006), Nucleic Acids Res. 34(12):3546-54. [PMID: 16855291]
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for taking one or two sequence datasets, peforming an intensive All by All BLAST and then tabulating
    the results as a series of pairwise comparisons (`*.gablam.*`):

    * Qry       = Query Short Name (or AccNum)
    * Hit       = Hit Short Name
    * Rank      = Rank of that Hit vs Query (based on Score)
    * Score     = BLAST Score (one-line)
    * E-Value   = BLAST E-value
    * QryLen    = Length of Query Sequence
    * HitLen    = Length of Hit Sequence
    * Qry_AlnLen        = Total length of local BLAST alignment fragments in Query (Unordered)
    * Qry_AlnID         = Number of Identical residues of Query aligned against Hit in local BLAST alignments (Unordered)
    * Qry_AlnSim        = Number of Similar residues of Query aligned against Hit in local BLAST alignments (Unordered)
    * Qry_OrderedAlnLen = Total length of local BLAST alignment fragments in Query (Ordered) 
    * Qry_OrderedAlnID  = Number of Identical residues of Query aligned against Hit in local BLAST alignments (Ordered)
    * Qry_OrderedAlnSim = Number of Similar residues of Query aligned against Hit in local BLAST alignments (Ordered)
    * Hit_AlnLen        = Total length of local BLAST alignment fragments in Hit (Unordered)
    * Hit_AlnID         = Number of Identical residues of Hit aligned against Query in local BLAST alignments (Unordered)
    * Hit_AlnSim        = Number of Similar residues of Hit aligned against Query in local BLAST alignments (Unordered)
    * Hit_OrderedAlnLen = Total length of local BLAST alignment fragments in Hit (Ordered)
    * Hit_OrderedAlnID  = Number of Identical residues of Hit aligned against Query in local BLAST alignments (Ordered)
    * Hit_OrderedAlnSim = Number of Similar residues of Hit aligned against Query in local BLAST alignments (Unordered)
    * ALIGN_ID          = Number of Identical residues as determined by pairwise ALIGN
    * ALIGN_Len         = Length of pairwise ALIGN

    By default, all BLAST hits will return alignments. (blastv=N blastb=N, where `N` is the size of `searchdb`.) This can
    be over-ridden by the blastv=X and blastb=X options to limit results to the top `X` hits.

    GABLAM will also produce a single table of summary statistics for all non-self hits (`*.hitsum.*`) (self hits included
    if selfhit=T selfsum=T):

    * Qry       = Query Short Name (or AccNum)
    * Hits      = Number of Hits
    * MaxScore  = Max non-self BLAST Score (one-line)
    * E-Value   = BLAST E-value for max score

    Version 2.8 onwards features explicit extra functionality for all-by-all searches, where the QueryDB (seqin=FILE) and
    SearchDB (searchdb=FILE) are the same. (Failing to give a searchdb will run in this mode.)

    Version 2.16 introduces a new "fullblast" mode, which performs a full BLAST search (using forks=X to set the number
    of processors for the BLAST search) followed by the blastres=FILE multiGABLAM processing. This should be faster for
    large datasets but precludes any appending of results files. This is incompatible with the missing=LIST advanced
    update option. (missing=LIST should only be required for aborted fullblast=F runs.)

Commandline:
    ### ~ Input/Search Options ~ ###
    seqin=FILE      : Query dataset file [infile.fas]
    searchdb=FILE   : Database to search. [By default, same as seqin]
    blastres=FILE   : BLAST results file for input (over-rides seqin and searchdb) [None]
    fullblast=T/F   : Whether to perform full BLAST followed by blastres analysis [False]
    blastp=X        : Type of BLAST search to perform (blastx for DNA vs prot; tblastn for Prot vs DNA) [blastp]
    gablamcut=X     : Min. percentage value for a GABLAM stat to report hit [0.0]  (GABLAM from FASTA only)
    cutstat=X       : Stat for gablamcut (eg. AlnLen or OrderedAlnSim. See above for full list) [OrderedAlnID]
    cutfocus=X      : Focus for gablamcut. Can be Query/Hit/Either/Both. [Either]
    localcut=X      : Cut-off length for local alignments contributing to global GABLAM stats) [0]

    ### ~ General Output Options ~ ###
    append=T/F      : Whether to append to output file or not. (Not available for blastres=FILE or fullblast=F) [False]
    fullres=T/F     : Whether to output full GABLAM results table [True]
    hitsum=T/F      : Whether to output the BLAST Hit Summary table [True]
    local=T/F       : Whether to output local alignment summary stats table [True]
    qassemble=T/F   : Whether to fully assemble query stats from all hits in HitSum [False]
    localmin=X      : Minimum length of local alignment to output to local stats table [0]
    selfhit=T/F     : Whether to include self hits in the fullres output [True] * See also selfsum=T/F *
    selfsum=T/F     : Whether to also include self hits in hitsum output [False] * selfhit must also be T *
    qryacc=T/F      : Whether to use the Accession Number rather than the short name for the Query [True]
    keepblast=T/F   : Whether to keep the blast results files rather than delete them [False]
    blastdir=PATH   : Path for blast results file (best used with keepblast=T) [./]
    percres=T/F     : Whether output is a percentage figures (2d.p.) or absolute numbers [True]
                      - Note that enough data is output to convert one into the other in other packages
    reduced=LIST    : List of terms that must be included in reduced output headers (e.g. Hit or Qry_Ordered) []

    ### ~ All-by-all Output Options ~ ###
    maxall=X        : Maximum number of sequences for all-by-all outputs [100]
    dismat=T/F      : Whether to output compiled distance matrix [True]
    diskey=X        : GABLAM Output Key to be used for distance matrix ['Qry_AlnID']
    distrees=T/F    : Whether to generate UPGMA tree summaries of all-by-all distances [True]
    treeformats=LIST: List of output formats for generated trees (see rje_tree.py) [nwk,text,png]
    disgraph=T/F    : Whether to output a graph representation of the distance matrix (edges = homology) [True]
    graphtypes=LIST : Formats for graph outputs (svg/xgmml/png/html) [xgmml,png]
    clusters=T/F    : Whether to output a list of clusters based on shared BLAST homology [True]
    bycluster=X     : Generate separate trees and distance matrix for clusters of X+ sequences [0]
    clustersplit=X  : Threshold at which clusters will be split (e.g. must be < distance to cluster) [1.0]
    singletons=T/F  : Whether to include singleton in main tree and distance matrix [False]
    saveupc=T/F     : Whether to output a UPC file for SLiMSuite compatibility [False]

    ### ~ Sequence output options ~ ###
    localalnfas=T/F : Whether to output local alignments to *.local.fas fasta file (if local=T) [False]
    fasout=T/F      : Output a fasta file per input sequence "ACCNUM.DBASE.fas" [False]   (GABLAM from FASTA only)
    fasdir=PATH     : Directory in which to save fasta files [BLASTFAS/]
    fragfas=T/F     : Whether to output fragmented Hits based on local alignments [False]
    gablamfrag=X    : Length of gaps between mapped residues for fragmenting local hits [100]
    addflanks=X     : Add flanking regions of length X to fragmented hits [0]
    combinedfas=T/F : Whether to generate a combined fasta file [False]

    ### ~ Advanced/Obselete Search/Output Options ~ ###
    dotplots=T/F    : Whether to use gablam.r to output dotplots. (Needs R installed an setup) [False]
    dotlocalmin=X   : Minimum length of local alignment to output to local hit dot plots [1000]
    mysql=T/F       : Whether to output column headers for mysql table build [False]
    missing=LIST    : This will go through and add missing results for AccNums in FILE (or list of AccNums X,Y,..) [None]
    startfrom=X     : Accession number to start from [None]
    alnstats=T/F    : Whether to output GABLAM stats or limit to one-line stats (blastb=0) [True]
    posinfo=T/F     : Output the Start/End limits of the BLAST Hits [True]
    outstats=X      : Whether to output just GABLAM, GABLAMO or All [All]

    ### ~ GABLAM Non-redundancy options. NOTE: These are different to rje_seq NR options. ~ ###
    nrseq=T/F       : Make sequences Non-Redundant following all-by-all. [False]
    nrcut=X         : Cut-off for non-redundancy filter, uses nrstat=X for either query or hit [100.0]
    nrstat=X        : Stat for nrcut (eg. AlnLen or OrderedAlnSim. See above for full list) [OrderedAlnID]
    nrchoice=LIST   : Order of decisions for choosing NR sequence to keep. Otherwise keeps first sequence. (swiss/nonx/length/spec/name/acc/manual) [swiss,nonx,length]
    nrsamespec=T/F  : Non-Redundancy within same species only. [False]
    nrspec=LIST     : List of species codes in order of preference (good to bad) []

    ### ~ BLAST Options ~ ###
    blastpath+=PATH : path for blast+ files [c:/bioware/blast+/] *Use fwd slashes
    blastpath=PATH  : path for blast files [c:/bioware/blast/] *Use fwd slashes
    blaste=X        : E-Value cut-off for BLAST searches (BLAST -e X) [1e-4]
    blastv=X        : Number of one-line hits per query (BLAST -v X)
    blastb=X        : Number of hit alignments per query (BLAST -b X)
    blastf=T/F      : Complexity Filter (BLAST -F X) [False]
    checktype=T/F   : Whether to check sequence types and BLAST program selection [True]

    ### ~ Additional ALIGN Global Identity ~ ###
    globid=T/F  : Whether to output Global %ID using ALIGN [False]
    rankaln=X   : Perform ALIGN pairwise global alignment for top X hits [0]
    evalaln=X   : Perform ALIGN pairwise global alignment for all hits with e <= X [1000]
    alncut=X    : Perform ALIGN pairwise global alignment until < X %ID reached [0]
    
    ### ~ Forking Options ~ ###
    noforks=T/F     : Whether to avoid forks [False]
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, re, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
#########################################################################################################################
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_db, rje_ppi, rje_seq, rje_seqlist, rje_tree
import rje_blast_V2 as rje_blast
import rje_dismatrix_V2 as rje_dismatrix
#########################################################################################################################
### History
#########################################################################################################################
### Major Functionality to Add
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - First working version based on BAM 1.4
    # 1.1 - Added blastres=FILE option
    # 1.2 - Added percres=T/F option
    # 1.3 - Added option to use a GABLAM stat as a cut-off
    # 1.4 - Added more output options for webserver
    # 2.0 - Major tidy up of module.
    # 2.1 - Added DNA GABLAM
    # 2.2 - Added reduced=LIST : List of terms that must be included in reduced output headers (e.g. Hit or Qry_Ordered)
    # 2.3 - Added local alignment stat output.
    # 2.4 - Added distance matrix output and visualisation.
    # 2.5 - Miscellaneous cleanup and bug fixes. Updated output file names to use basefile.
    # 2.6 - Added full name output for NSF tree.
    # 2.7 - Added cluster output.
    # 2.8 - Added graph output and replaced dispng and disnsf with distrees. Replaced dismatrix PNG with just tree.
    # 2.9 - Added LocalMin and LocalCut for controlling how local alignments are output and/or contribute to GABLAM.
    # 2.10- Added dot plot PNG output using gablam.r. Added sequence checking with rje_seqlist -> correct BLAST type.
    # 2.11- Altered to use BLAST+ and rje_blast_V2.
    # 2.12- Consolidated use of BLAST V2.
    # 2.13- Fixed Protein vs DNA GABLAM. Modified sequence extraction to handle larger sequences. Add blastdir=PATH/.
    # 2.14- Added checktype=T/F option to check sequence/BLAST type.
    # 2.15.0 - Added seqnr function. Add run() method.
    # 2.16.0 - Added fullblast=T/F : Whether to perform full BLAST followed by blastres analysis [False]
    # 2.16.1 - Fixed a bug where the fullblast option was failing to return scores and evalues.
    # 2.17.0 - Added localalnfas=T/F : Whether to output local alignments to *.local.fas fasta file (if local=T) [False]
    # 2.17.1 - Fixed bug where query and hit lengths were not being output for fullblast.
    # 2.18.0 - Added blaste filtering to be applied to existing BLAST results.
    # 2.19.0 - Added maxall=X limits to all-by-all analyses. Added qassemble=T.
    # 2.19.1 - Fixed handling of basefile and results generation for blastres=FILE.
    # 2.19.2 - Modified output to be in rank order.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Implement a DNA version
    # [Y] : blastres=FILE   : Generate table from file of BLAST results rather than input sequences (for webserver)
    # [Y] : Figure(s) for pairwise comparisons (for webserver)
    # [ ] : Reduced Log Output?
    # [Y] : Add posinfo=T/F for start/end of BLAST alignments
    # [ ] : Add output=Ordered/Unordered/All to restrict columns
    # [Y] : Add orientation for DNA BLAST
    # [Y] : Add extraction (and assembly?) of local alignment sequences into fasta file.
    # [Y] : Tidy up naming of files etc. using basefile.
    # [ ] : Convert to new object style (V3)
    # [ ] : Replace old rje_seq SeqList (and memsaver) with new rje_seqlist
    # [X] : Dismatrix PNG output needs some work!
    # [Y] : Consider converting cluster output to UPC format!
    # [ ] : Check sequences and auto-detect right BLAST method.
    # [ ] : Add basic assembly function (or program) to stick contigs together in right order.
    # [ ] : Add UPC clustering and check XGMML output. Replace UPC method in rje_slimcore.
    # [ ] : Add citation.
    # [ ] : Look to generate searchGABLAM output directly from rje_blast_V2?
    # [ ] : Massively improve naming of gablam fragments and carry through of original protein names.
    # [ ] : Add additional NR assessment based on annotation and/or gene type.
    # [?] : Add     nrkeepann=T/F   : Append annotation of redundant sequences onto NR sequences [False] ?
    # [ ] : Add an flankorf=T/F function for extending fragfas protein hits to up/downstream Met/STOP codons. (5' M?)
    # [Y] : Add a sequence number limit for dismat output to avoid excessively long runtimes.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('GABLAM', '2.19.2', 'March 2015', '2006')
    description = 'Global Analysis of BLAST Local AlignMents'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please cite: Davey, Shields & Edwards (2006), Nucleic Acids Res. 34(12):3546-54. [PMID: 16855291]']
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc   __ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GABLAM Class:                                                                                           #
#########################################################################################################################
class GABLAM(rje.RJE_Object):     
    '''
    Global Analysis of BLAST Local AlignMents (GABLAM) Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of comparison querydb.vs.searchdb
    - QueryDB = Query dataset file [infile.fas]
    - SearchDB = Database to search. [By default, same as seqin]
    - BlastDir = Path for blast results file (best used with keepblast=T) [./]
    - BlastRes = BLASTP results file for input [None]
    - StartFrom = Accession number to start from [None]
    - FasDir = Directory in which to save fasta files [BLASTFAS/]
    - CutStat = Stat for gablamcut (eg. AlnLen or OrderedAlnSim. See above for full list) [OrderedAlnID]
    - CutFocus = Focus for gablamcut. Can be Qry/Hit/Either/Both. [Either]
    - GablamOut = Name of GABLAM Output file (by default seqin.vs.searchdb.gablam.tdt) [None]
    - HitSumOut = Name of GABLAM Output file (by default seqin.vs.searchdb.gablam.tdt) [None]   
    - LocalOut = Name of Local alignment stat output file (by default seqin.vs.searchdb.local.tdt) [None]
    - DisMatOut = Name of Local alignment stat output file (by default seqin.vs.searchdb.dis.tdt) [None]
    - DisKey = GABLAM Output Key to be used for distance matrix ['Qry_AlnID']
    - OutStats = Whether to output just GABLAM, GABLAMO or All [All]
    - NRStat=X        : Stat for nrcut (eg. AlnLen or OrderedAlnSim. See above for full list) [OrderedAlnID]

    Opt:boolean
    - AlnStats = Whether to output alignment stats or limit to one-line stats (blastb=0) [True]
    - ByCluster = Whether trees and distance matrix should be by cluster [False]
    - CheckType=T/F   : Whether to check sequence types and BLAST program selection [True]
    - Clusters = Whether to output a list of clusters based on shared BLAST homology [False]
    - CombinedFas = Whether to generate a combined fasta file [False]
    - DisMat = Whether to output compiled distance matrix [True]
    X- DisNSF = Whether to output NSF text tree based on distance matrix [True]
    X- DisPNG = Whether to generate PNG visualisation of summary distance matrix [True]
    - DisTrees = Whether to generate UPGMA tree summaries of all-by-all distances [True]
    - DisGraph = Whether to output a graph representation of the distance matrix (edges = homology) [True]
    - DotPlots = Whether to use gablam.r to output dotplots. (Needs R installed an setup) [False]
    - FullBlast = Whether to perform full BLAST followed by blastres analysis [False]
    - FragFas = Whether to output fragmented Hits based on local alignments [False]
    - GlobID = Whether to output Global %ID using ALIGN [False]
    - LocalAlnFas = Whether to output local alignments to *.local.fas fasta file (if local=T) [False]
    - SaveUPC = Whether to output a UPC file for SLiMSuite compatibility [False]
    - SelfHit = Whether to include self hits in the output [True]
    - SelfSum = Whether to also include self hits in hitsum output [False] * selfhit must also be T *
    - Singletons = Whether to include singleton in main tree and distance matrix [False]
    - QAssemble = Whether to fully assemble query stats from all hits [False]
    - QryAcc = Whether to use the Accession Number rather than the short name for the Query [True]
    - MySQL = Whether to output column headers for mysql table build [False]
    - FasOut = Output a fasta file per input sequence "AccNum.blast.fas" [False]
    - FullRes = Whether to output full GABLAM results table [True]
    - HitSum = Whether to output the BLAST Hit Summary table [True]
    - Local = Whether to output local alignment summary stats table [True]
    - PercRes = Whether output is a percentage figures (3d.p.) or absolute numbers [True]
    - KeepBlast = Whether to keep the blast results files rather than delete them [False]
    - PosInfo = Output the Start/End limits of the BLAST Hits [True]
    - NRSeq=T/F       : Make sequences Non-Redundant following all-by-all. [False]
    - NRSameSpec=T/F      : Non-Redundancy within same species only. [False]
    - NRKeepAnn=T/F   : Append annotation of redundant sequences onto NR sequences [False]

    Stat:numeric
    - RankAln = Perform ALIGN pairwise global alignment for top X hits [0]
    - EvalAln = Perform ALIGN pairwise gloabl alignment for all hits with e <= X [-1]
    - AddFlanks = Add flanking regions of length X to fragmented hits [0]
    - AlnCut =  Perform ALIGN pairwise global alignment until < X %ID reached [0]
    - BlastE=X            # BLAST evalue cut-off. Could/should be used for filtering even if re-using GABLAM [1e-4]
    - ClusterSplit = Threshold at which clusters will be split (e.g. must be < distance to cluster) [1.0]
    - DBSize =  Size of BLAST Database
    - DotLocalMin = Minimum length of local alignment to output to local hit dot plots [1000]
    - MaxAll = Maximum number of sequences for all-by-all outputs [100]
    - NRCut = Cut-off for non-redundancy filter, uses cutstat=X for either query or hit [100.0]
    - GABLAMCut = Min. percentage value for a GABLAM stat to report hit [0.0]
    - LocalCut = Cut-off length for local alignments contributing to global GABLAM stats) [0]
    - LocalMin = Minimum length of local alignment to output to local stats table [0]

    List:list
    - FullResHeaders = List of headers for FullRes output
    - GraphTypes = Formats for graph outputs (tdt/r/svg/xgmml/png) [xgmml,png]
    - HitSumHeaders = List of headers for FullRes output
    - Missing = This will go through and add missing results for AccNums in FILE (or list of AccNums X,Y,..) [None]
    - Reduced = List of terms that must be included in reduced output headers (e.g. Hit or Qry_Ordered) []
    - NRChoice=LIST   : Order of decisions for choosing NR sequence to keep. Otherwise keeps first sequence. (db/nonx/length/spec/name/acc/manual) [spec,db,nonx,length,manual]
    - NRSource=LIST     : List of databases in order of preference (good to bad) [sprot,ipi,uniprot,trembl,ens_known,ens_novel,ens_scan]
    - NRSpec=LIST   : List of species codes in order of preference (good to bad) []

    Dict:dictionary
    - NameMap = Mapping dictionary of {AccNum/Short:Full name}

    Obj:RJE_Objects
    - SeqList = SeqList object containing GABLAM query sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','QueryDB','SearchDB','BlastRes','StartFrom','FasDir','CutStat','CutFocus','BlastDir',
                         'GablamOut','HitSumOut','OutStats','BLASTP','LocalOut','DisMatOut','DisKey','NRStat']
        self.optlist = ['AlnStats','Singletons','GlobID','SelfHit','SelfSum','QryAcc','MySQL','FasOut','FullRes',
                        'PercRes','HitSum','KeepBlast','PosInfo','Local','LocalAlnFas','DisMat','FragFas','CheckType',
                        'CombinedFas','Clusters','SaveUPC','DisTrees','DisGraph','DotPlots','FullBlast',
                        'NRSeq','NRSameSpec','NRKeepAnn','QAssemble']
        self.statlist = ['AddFlanks','BlastE','RankAln','EvalAln','AlnCut','DBSize','GABLAMCut','ByCluster','LocalCut',
                         'LocalMin','DotLocalMin','NRCut','MaxAll']
        self.listlist = ['Missing','FullResHeaders','HitSumHeaders','Reduced','GraphTypes','NRChoice','NRSpec','NRSource']
        self.dictlist = ['NameMap']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'QueryDB':'infile.fas','FasDir':rje.makePath('BLASTFAS'),'CutStat':'OrderedAlnID','BlastDir':'',
                      'CutFocus':'Either','OutStats':'All','BLASTP':'blastp','DisKey':'Qry_AlnID','NRStat':'OrderedAlnID'})
        self.setOpt({'AlnStats':True,'MemSaver':True,'SelfHit':True,'QryAcc':True,'FullRes':True,'HitSum':True,'CheckType':True,
                     'PercRes':True,'PosInfo':True,'Local':True,'DisMat':True,'DisTrees':True,'DisGraph':True,
                     'NRSeq':False,'NRSameSpec':False,'NRKeepAnn':False})
        self.setStat({'AddFlanks':0,'EvalAln':1000,'GABLAMCut':0.0,'ByCluster':0,'ClusterSplit':1.0,'LocalCut':0,
                      'LocalMin':0,'DotLocalMin':1000,'NRCut':100.0,'BlastE':1e-4,'MaxAll':100})
        self.list['GraphTypes'] = ['xgmml','png']
        self.list['NRChoice'] = ['db','nonx','length','manual']
        self.list['NRSource'] = ['sprot','ipi','uniprot','trembl','ens_known','ens_novel','ens_scan']
        ### Other Attributes ###
        self._setForkAttributes()
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### Read Commands ###
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                self._forkCmd(cmd)  
                ### Class Options ###
                self._cmdRead(cmd,type='file',att='QueryDB',arg='seqin')  
                self._cmdReadList(cmd,'path',['FasDir','BlastDir'])
                self._cmdReadList(cmd,'file',['QueryDB','SearchDB','BlastRes'])
                self._cmdReadList(cmd,'info',['StartFrom','CutStat','CutFocus','OutStats','BLASTP','DisKey','NRStat'])
                self._cmdReadList(cmd,'opt',['AlnStats','Singletons','GlobID','SelfHit','SelfSum','QryAcc','MySQL','DotPlots',
                                             'FasOut','FullRes','FragFas','SaveUPC','HitSum','PercRes','KeepBlast','CheckType',
                                             'PosInfo','Local','LocalAlnFas','DisMat','DisTrees','DisGraph','CombinedFas','Clusters',
                                             'NRSeq','NRSameSpec','NRKeepAnn','FullBlast','QAssemble'])
                self._cmdReadList(cmd,'stat',['RankAln','EvalAln','AlnCut','GABLAMCut','ClusterSplit','NRCut','BlastE'])
                self._cmdReadList(cmd,'int',['AddFlanks','ByCluster','LocalCut','LocalMin','DotLocalMin','MaxAll'])
                self._cmdReadList(cmd,'list',['Missing','Reduced','GraphTypes','NRChoice'])
                self._cmdReadList(cmd,'lclist',['NRChoice','NRDB'])
                self._cmdReadList(cmd,'uclist',['NRSpec'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        ### Additional Processing ###
        if self.info['SearchDB'].lower() in ['','none']: self.info['SearchDB'] = self.info['QueryDB']
        if self.getBool('FragFas'): self.opt['FasOut'] = True; self.cmd_list = ['gablamfrag=100'] + self.cmd_list
        if self.opt['FasOut'] and not os.path.exists(self.info['FasDir']): rje.mkDir(self,self.info['FasDir'])
        if self.stat['GABLAMCut'] < 1.0: self.stat['GABLAMCut'] *= 100.0
        ### Update options for FullBlast=True ###
        if self.getBool('Local') and self.getBool('LocalAlnFas') and not self.getBool('FullBlast'):
            if self.i() < 0 or rje.yesNo('localnfas=T incompatible with fullblast=F. Switch fullblast=T?'):
                self.printLog('#CMD','Switching fullblast=T. Incompatible with localnfas=T.')
                self.setBool({'FullBlast':True})
            else:
                self.printLog('#CMD','Switching localnfas=F. Incompatible with fullblast=F.')
                self.setBool({'LocalAlnFas':False})
        if self.getBool('FullBlast') and self.list['Missing']:
            if self.i() >= 0 and rje.yesNo('Missing=LIST incompatible with fullblast=T. Switch fullblast=F?',default='N'):
                self.printLog('#CMD','Switching fullblast=F. Incompatible with missing=LIST.')
                self.setBool({'FullBlast':False,'LocalAlnFas':False})
            else:
                self.printLog('#CMD','Switching off missing=LIST. Incompatible with fullblast=T.')
                self.list['Missing'] = []
        if self.getBool('FullBlast') and self.getStrLC('StartFrom'):
            if self.i() < 0 or rje.yesNo('startfrom=X incompatible with fullblast=T. Switch fullblast=F?'):
                self.printLog('#CMD','Switching fullblast=F. Incompatible with startfrom=X.')
                self.setBool({'FullBlast':False,'LocalAlnFas':False})
            else:
                self.printLog('#CMD','Switching off startfrom=X. Incompatible with fullblast=T.')
                self.setStr({'StartFrom':'None'})
        if self.getBool('FullBlast') and self.getBool('Append'):
            self.printLog('#CMD','Switching append=F. Incompatible with fullblast=T.')
            self.setBool({'Append':False})
        #self.deBug(self.info)
        #self.debug(self.opt)
#########################################################################################################################
    ### <2> ### Run Setup Section                                                                                       #
#########################################################################################################################
    def setupSearchFiles(self): ### Sets up Results files, reformatting searchdb to fasta if necessary
        '''Sets up Results files, reformatting searchdb to fasta if necessary.'''
        try:### ~ [0] ~ QueryDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while not rje.checkForFile(self.getStr('QueryDB')):
                self.errorLog('IOError: Cannot find %s!' % self.getStr('QueryDB'),False,False)
                if self.i() < 0: raise ValueError
                else: self.setStr({'QueryDB':rje.choice('Input File name for query sequences?:',default=self.getStr('SearchDB'))})
            ### ~ [1] ~ SearchDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while not rje.checkForFile(self.getStr('SearchDB')):
                self.errorLog('IOError: Cannot find %s!' % self.getStr('SearchDB'),False,False)
                if self.i() < 0: raise ValueError
                else: self.setStr({'SearchDB':rje.choice('Input File name for search database?:',default=self.getStr('QueryDB'))})
            ## ~ [1a] ~ Reformat to Fasta if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if open(self.getStr('SearchDB'), 'r').readline()[:1] != '>':     # Not fasta format!
                self.printLog('#DB','Reformatting search database to fasta.')
                rlist = ['seqin=%s' % self.getStr('SearchDB'), 'seqout=%s.fas' % self.getStr('SearchDB'), 'reformat=fasta',
                         'seqmode=file', 'autoload=T', 'autofilter=T']
                sdbseq = rje_seqlist.SeqList(log=self.log,cmd_list=rlist)
                self.setStr({'SearchDB':'%s.fas' % self.getStr('SearchDB')})
            else: 
                rlist = ['seqin=%s' % self.getStr('SearchDB'), 'seqmode=file', 'autoload=T', 'autofilter=F']
                sdbseq = rje_seqlist.SeqList(log=self.log,cmd_list=rlist)
            ## ~ [1b] ~ BLAST Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qlist = ['seqin=%s' % self.getStr('QueryDB'), 'seqmode=file', 'autoload=T', 'autofilter=F']
            qseq = rje_seqlist.SeqList(log=self.log,cmd_list=qlist)
            #if self.dev(): blast = rje_blast.blastObj(self.log,['blastf=F','gablamfrag=100']+self.cmd_list,type='New')   #BLASTRun(log=self.log,cmd_list=['blastf=F','gablamfrag=100']+self.cmd_list)
            #else: blast = rje_blast.blastObj(self.log,['blastf=F','gablamfrag=100']+self.cmd_list,type='Old')
            blast = rje_blast.blastObj(self.log,['blastf=F','gablamfrag=100']+self.cmd_list,type='dev')
            btype = blast.getStr('Type')
            if self.getBool('CheckType'): blast.checkProg(qtype=qseq.readSeqType(log=False),stype=sdbseq.readSeqType(log=False))
            if btype != blast.getStr('Type'): self.cmd_list += ['blastp=%s' % blast.getStr('Type')]
            if blast.getStr('Type') in ['blastn','tblastn','tblastx']:  blast.formatDB(fasfile=self.getStr('SearchDB'),protein=False,force=False)
            else: blast.formatDB(fasfile=self.getStr('SearchDB'),protein=True,force=False)
            self.setInt({'DBSize':rje_seq.SeqCount(self,self.getStr('SearchDB'))})
        except: self.errorLog('Error in setupSearchFiles()',printerror=True); raise
#########################################################################################################################
    def setupSeqList(self): ### Sets up SeqList, including any changes due to self.list['Missing']
        '''Sets up SeqList, including any changes due to self.list['Missing']. Returns SeqList.'''
        try:### ~ [0] ~ Basic SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Needs tidying and V3 conversion #!#
            self.obj['SeqList'] = None
            if self.getBool('FullBlast'): seqlist = rje_seqlist.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=T'])
            elif not self.opt['MemSaver']:
                seqlist = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list)  #+['autofilter=F'])
            else:
                seqlist = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None'])    #,'autofilter=F'])
            seqlist.setStr({'Name':self.info['QueryDB']})
            self.obj['SeqList'] = seqlist

            ### Missing AccNums ###     #!# Does this still happen? #!#
            if self.list['Missing']:
                # Make sequence object #
                SEQFILE = open(self.info['QueryDB'], 'r')
                lastline = 'Start'
                (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                while nextseq:
                    if nextseq.shortName() in self.list['Missing']: # OK
                        self.list['Missing'].remove(nextseq.shortName())
                    elif nextseq.info['AccNum'] in self.list['Missing']: # OK
                        self.list['Missing'].remove(nextseq.info['AccNum'])
                    else:
                        seqlist.seq.remove(nextseq) # Not in list
                    (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                self.log.printLog('#SEQ','%d missing sequences read from %s.' % (seqlist.seqNum(),seqlist.info['Name']))
                self.log.printLog('#SEQ','%d missing sequences not found in %s.' % (len(self.list['Missing']),seqlist.info['Name']))
                for missing in self.list['Missing']:
                    self.log.printLog('#ERR','%s not found in %s!' % (missing,seqlist.info['Name']))
                if seqlist.seqNum() == 0:
                    self.log.printLog('No sequences to process in GABLAM!',False,False)
                    return None
                self.opt['MemSaver'] = False
                if not self.opt['Append'] and self.stat['Interactive'] >= 0 and rje.yesNo('Adding missing sequences. Append results?'):
                    self.opt['Append'] = True
                self.obj['SeqList'] = seqlist

        except:
            self.log.errorLog('Error in setupSeqList()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def setupResultsFiles(self,saveheads=True):    ### Sets up Results Files and self.list headers
        '''Sets up Results Files and self.list headers.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.basefile().lower() in ['','none']:
                if self.getStr('QueryDB') == self.getStr('SearchDB'): self.setBasefile(rje.baseFile(self.getStr('QueryDB'),strip_path=True))
                else: self.setBasefile('%s.vs.%s' % (rje.baseFile(self.getStr('QueryDB'),strip_path=True),rje.baseFile(self.getStr('SearchDB'),strip_path=True)),cascade=False)
            delimit_ext = rje.delimitExt(rje.getDelimit(self.cmd_list))
            if self.getBool('DisMat') and self.getStr('QueryDB') != self.getStr('SearchDB'):
                self.printLog('#DIS','Cannot generate Distance Matrix outputs if QueryDB and SearchDB differ.')
                self.setBool({'DisMat':False,'DisTrees':False,'DisGraph':False})
            elif self.getStr('QueryDB') != self.getStr('SearchDB'): pass
            elif not self.getBool('DisMat') and (self.getBool('DisTrees') or self.getBool('DisGraph')):
                self.printLog('#DIS','Cannot generate all-by-all trees or graph without DisMat=T')
                self.setBool({'DisMat':self.i() >=0 and rje.yesNo('Switch on Distance Matrix output?','N')})
                if not self.getBool('DisMat'): self.setBool({'DisMat':False,'DisTrees':False,'DisGraph':False})
                #self.setBool({'DisMat':self.getBool('DisMat') or self.getBool('DisTrees') or self.getBool('DisGraph')})
            self.setBool({'QryAcc':self.getBool('QryAcc') and not self.getBool('DisMat')})
            self.setBool({'Gablam':self.getBool('FullRes') or self.getBool('DisMat')})
            ### ~ [1] ~ Setup Results Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for res in ['Gablam','HitSum','Local','DisMat']:
                res_out = '%sOut' % res
                if not self.getBool(res): continue
                if self.getStr(res_out).lower() in ['false','f']: self.setBool({res:False}); continue
                ## ~ [1a] ~ Filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                file_ext = '%s.%s' % (res.lower(),delimit_ext)
                if self.getStr(res_out).lower() in ['','none','t','true']: self.setStr({res_out:'%s.%s' % (self.basefile(),file_ext)})
                ## ~ [1b] ~ Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if res == 'Gablam':
                    headlist = ['Qry','Hit','Rank','Score','EVal','QryLen','HitLen']
                    outs = []
                    if self.getStr('OutStats').upper() != 'GABLAMO': outs.append('')
                    if self.getStr('OutStats').upper() != 'GABLAM': outs.append('Ordered')
                    hds = ['AlnLen','AlnID','AlnSim']
                    if self.info['BLASTP'].lower() != 'blastp': hds += ['Dirn']     # Context-specific use of Dirn
                    if self.opt['PosInfo']: hds += ['Start','End']
                    if self.opt['AlnStats']:
                        for qh in ['Qry','Hit']:
                            for o in outs:
                                for hd in hds:
                                    newhead = '%s_%s%s' % (qh,o,hd)
                                    if not self.list['Reduced']: headlist.append(newhead)
                                    else:
                                        for red in self.list['Reduced']:
                                            if newhead.find(red) >= 0:
                                                headlist.append(newhead)
                                                break
                    if self.opt['GlobID']: headlist += ['ALIGN_ID','ALIGN_Len']
                    self.list['FullResHeaders'] = headlist[0:]
                elif res == 'HitSum':
                    headlist = ['Qry','HitNum','MaxScore','EVal']
                    if self.getBool('QAssemble'): headlist += ['Length','Coverage','Identity','Positives']
                    self.list['HitSumHeaders'] = headlist[0:]
                    #self.debug(self.list['HitSumHeaders'])
                elif res == 'Local':
                    headlist = ['Qry','Hit','AlnNum','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd']
                    self.list['LocalHeaders'] = headlist[0:]
                elif res == 'DisMat':
                    self.obj['DisMat'] = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
                    self.obj['DisMat'].info['Name'] = '%s vs %s GABLAM' % (rje.baseFile(self.getStr('QueryDB'),strip_path=True),rje.baseFile(self.getStr('SearchDB')))
                    if not self.getBool('PercRes'):
                        self.setBool({'PercRes':True})
                        self.printLog('#PERC','Switched percres=T for distance matrix output')
                if saveheads and res != 'DisMat' and not self.getBool('FullBlast'): rje.delimitedFileOutput(self,self.getStr(res_out),headlist,rje_backup=True)
        except: self.errorLog('Error in setupResultsFiles()',printerror=True); raise
#########################################################################################################################
    ### <3> ### Main Run Methods Section                                                                                #
#########################################################################################################################
    def run(self):  ### Main Run Method. Selects NRSeq versus GABLAM run modes. Handles data return.
        '''Main Run Method. Selects NRSeq versus GABLAM run modes. Handles data return.'''
        try:### ~ [1] ~ NRSeq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('NRSeq'): return self.nrSeq()
            ### ~ [2] ~ GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else: return self.gablam()
        except SystemExit: sys.exit()
        except: self.errorLog('Error in GABLAM.run()',printerror=True,quitchoice=True)
#########################################################################################################################
    def nrSeq(self):  ### Main Run Method. Selects NRSeq versus GABLAM run modes. Handles data return.
        '''Main Run Method. Selects NRSeq versus GABLAM run modes. Handles data return.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('NRSeq') and self.getStr('QueryDB') != self.getStr('SearchDB'):
                self.printLog('#NRSEQ','Cannot perform NR filter if QueryDB and SearchDB differ.')
                return False
            if not self.basefile(): self.setBasefile(rje.baseFile(self.getStr('QueryDB'),strip_path=True))
            delimit_ext = rje.delimitExt(rje.getDelimit(self.cmd_list))
            gfile = '%s.gablam.%s' % (self.basefile(),delimit_ext)
            gqfile = '%s.gablam.%s' % (rje.baseFile(self.getStr('QueryDB'),strip_path=True),delimit_ext)
            gablam = not rje.exists(gfile) and not rje.exists(gqfile)
            db = rje_db.Database(self.log,self.cmd_list)
            gablam = self.force() or gablam
            qcut = 'Qry_%s' % self.getStr('NRStat')
            hcut = 'Hit_%s' % self.getStr('NRStat')
            nrcut = self.getNum('NRCut')
            if nrcut < 1.0: nrcut = 100.0 * nrcut
            self.printLog('#NR','GABLAM NR filter: %s >= %.1f%%' % (self.getStr('NRStat'),nrcut))
            ### ~ [2] ~ Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gablam: self.gablam()
            else:
                if not rje.exists(gfile): gfile = gqfile
                self.printLog('#GABLAM','Results file %s found (force=F)' % gfile)
            ### ~ [3] ~ Load/Process Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = db.addTable(gfile,['Qry','Hit'],name='gablam')
            if qcut not in gdb.fields(): raise ValueError('Cannot find nrstat field "%s" in %s' % (qcut,gfile))
            if hcut not in gdb.fields(): raise ValueError('Cannot find nrstat field "%s" in %s' % (hcut,gfile))
            gdbformat = {}
            gscores = {True:'float',False:'int'}[self.getBool('PercRes')]
            gfields = []    # Fields containing GABLAM Scores (could be used for filtering)
            for field in gdb.fields():
                if field in ['Rank','QryLen','HitLen'] or field[-3:] in ['End','art']: gdbformat[field] = 'int'
                elif field not in ['Qry','Hit']:
                    gdbformat[field] = gscores
                    if field[:3] in ['Qry','Hit']: gfields.append(field)
                else: gdbformat[field] = 'str'
            gdb.dataFormat(gdbformat)
            if not self.getBool('PercRes'):
                ex = 0.0; etot = gdb.entryNum()
                for entry in gdb.entries():
                    self.progLog('\r#PERC','Converting GABLAM results to percentage figures: %.1f%%' % (ex/etot)); ex += 100.0
                    for field in gfields:
                        qh = field[:3]
                        entry[field] = 100.0 * entry[field] / entry['%sLen' % qh]
                self.printLog('\r#PERC','Converted GABLAM results to percentage figures using Qry/Hit lengths.')
            ex = 0.0; etot = gdb.entryNum(); qx = 0; hx = 0
            for entry in gdb.entries():
                self.progLog('\r#ACC','Checking query/hit accnum: %.1f%%' % (ex/etot)); ex += 100.0
                qry = string.split(entry['Qry'],'__')[-1]
                if qry != entry['Qry']: entry['Qry'] = qry; qx += 1
                hit = string.split(entry['Hit'],'__')[-1]
                if hit != entry['Hit']: entry['Hit'] = hit; hx += 1
            self.printLog('\r#ACC','%s entry queries converted from shortname to accnum.' % rje.iStr(qx))
            self.printLog('\r#ACC','%s entry hits converted from shortname to accnum.' % rje.iStr(hx))
            if qx or hx: gdb.remakeKeys()
            ### ~ [4] ~ SeqNR process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('QueryDB'),'autoload=F'])
            seqin.loadSeq()
            seqdict = seqin.makeSeqNameDic('acc')
            #self.debug('%s ... %s' % (rje.sortKeys(seqdict)[:10],rje.sortKeys(seqdict)[-10:]))
            nrbad = []
            gdb.index('Qry'); gdb.index('Hit')
            ex = 0.0; etot = gdb.entryNum()
            for ekey in gdb.dataKeys():
                entry = gdb.data(ekey)
                self.progLog('\r#NR','GABLAM NR filtering: %.1f%%' % (ex/etot)); ex += 100.0
                if entry['Qry'] in nrbad or entry['Hit'] in nrbad: continue
                if entry['Qry'] == entry['Hit']: continue
                if entry[qcut] < nrcut or entry[hcut] < nrcut: continue
                qseq = seqdict[entry['Qry']]
                hseq = seqdict[entry['Hit']]
                if self.getBool('NRSameSpec') and seqin.seqSpec(qseq) != seqin.seqSpec(hseq): continue
                # Establish which sequence to keep!
                qname = seqin.shortName(qseq); hname = seqin.shortName(hseq)
                nrtext = '%s (%.1f%%) vs %s (%.1f%%)' % (qname,entry[qcut],hname,entry[hcut])
                #self.debug(nrtext)
                nrdump = (None,None)   # Tuple of (Reject,Reason)
                for nrchoice in self.list['NRChoice']:
                    # swiss/nonx/length/spec/name/acc/manual
                    if nrchoice == 'swiss': # Favour SwissProt over uniprot
                        if seqin.isSwiss(qseq) == seqin.isSwiss(hseq): continue
                        if seqin.isSwiss(qseq): nrdump = ('Hit','Swissprot vs non Swissprot')
                        else: nrdump = ('Qry','Swissprot vs non Swissprot')
                    if nrchoice == 'nonx':
                        qlen = seqin.seqNonX(qseq); hlen = seqin.seqNonX(hseq)
                        if qlen == hlen: continue
                        elif qlen > hlen: nrdump = ('Hit','%s vs %s non-X positions' % (rje.iStr(hlen),rje.iStr(qlen)))
                        else: nrdump = ('Qry','%s vs %s non-X positions' % (rje.iStr(qlen),rje.iStr(hlen)))
                    if nrchoice == 'len':
                        qlen = seqin.seqLen(qseq); hlen = seqin.seqLen(hseq)
                        if qlen == hlen: continue
                        elif qlen > hlen: nrdump = ('Hit','%s vs %s positions' % (rje.iStr(hlen),rje.iStr(qlen)))
                        else: nrdump = ('Qry','%s vs %s positions' % (rje.iStr(qlen),rje.iStr(hlen)))
                    if nrchoice == 'spec':
                        qspec = seqin.seqSpec(qseq); hspec = seqin.seqSpec(hseq)
                        if qspec == hpsec: continue
                        if qspec not in self.list['NRSpec'] and hspec not in self.list['NRSpec']: continue
                        if qspec not in self.list['NRSpec'] or self.list['NRSpec'].index(hspec) < self.list['NRSpec'].index(qspec): nrdump = ('Qry','Favoured species %s vs %s' % (hspec,qspec))
                        else: nrdump = ('Hit','Favoured species %s vs %s' % (qspec,hspec))
                    if nrchoice == 'name':
                        qname = seqin.shortName(qseq); hname = seqin.shortName(hseq);
                        nrnames = [qname,hname]; nrnames.sort()
                        if nrnames[0] == qname: nrdump = ('Hit','Decision based on name sorting')
                        else: nrdump = ('Qry','Decision based on name sorting')
                    if nrchoice in ['acc','accnum']:
                        qname = seqin.seqAcc(qseq); hname = seqin.seqAcc(hseq);
                        nrnames = [qname,hname]; nrnames.sort()
                        if nrnames[0] == qname: nrdump = ('Hit','Decision based on accnum sorting')
                        else: nrdump = ('Qry','Decision based on accnum sorting')
                    if nrchoice == 'manual':
                        while not nrdump[0]:
                            choice = rje.choice('%s: Keep <Q>ry %s or <H>it %s?' % (nrtext,qname,hname)).upper()[:1]
                            if choice not in 'QH': continue
                            nrdump = ({'Q':'Hit','H':'Qry'}[choice],'Manual selection')
                            if not rje.yesNo('Reject %s %s?' % (nrdump[0],entry[nrdump[0]])): nrdump = (None,None)
                    if nrdump[0]: break
                if not nrdump[0]:   # Keep first sequence
                    if seqin.seqs().index(qseq) < seqin.seqs().index(hseq): nrdump = ('Hit','Kept first sequence')
                    else: nrdump = ('Qry','Kept first sequence')
                # Add to rejects and report
                nrbad.append(entry[nrdump[0]])
                self.printLog('\r#REJ','%s: %s rejected - %s' % (nrtext,entry[nrdump[0]],nrdump[1]))
            self.printLog('\r#NR','GABLAM NR filter, %s >= %.1f%%: %s of %s sequences rejected.' % (self.getStr('NRStat'),nrcut,rje.iLen(nrbad),rje.iStr(seqin.seqNum())))
            ### ~ [5] ~ Save NR Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin.list['BadAcc'] = nrbad
            seqin.filterSeqs()
            seqin.saveSeq(seqfile='%s.nr%d.fas' % (self.basefile(),int(self.getNum('NRCut'))))
            return True
        except SystemExit: sys.exit()
        except: self.errorLog('Error in GABLAM.nrSeq()',printerror=True,quitchoice=True); return False
#########################################################################################################################
    def gablam(self):      ### Main GABLAM Method
        '''
        Main GABLAM Method. This makes use of two forking sets: 1 for running BLAST and 1 for processing BLAST results.
        Returns True if successful, else False.
        '''
        try:### ~ [0] ~ Special Results From existing BLAST Results File run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('BlastRes').lower():
                while not rje.checkForFile(self.getStr('BlastRes')):
                    self.errorLog('IOError: Cannot find BLAST results file %s!' % self.getStr('BlastRes'),False,False)
                    if self.i() < 0: sys.exit()
                    else: self.setStr({'BlastRes':rje.choice('Input File name for BLAST results file:')})
                return self.GABLAMFromBLASTRes()
            ### ~ [1] ~ Setup Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Input Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupSearchFiles()
            ## ~ [1b] ~ SeqList Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupSeqList()
            if not self.obj['SeqList']: return False
            ## ~ [1c] ~ StartFrom ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            startfrom = None
            if self.getStrLC('StartFrom'):
                startfrom = self.getStr('StartFrom')
                if not self.getBool('Append') and self.i() >= 0: self.setBool({'Append':rje.yesNo('Starting from %s. Append results?' % startfrom)})
            self.setStr({'StartFrom': startfrom})
            ## ~ [1d] ~ Results Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupResultsFiles()
            ### ~ [2] ~ Run GABLAM and Output results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if (self.getBool('FullBlast') and self.fullBlastGABLAM()) or (not self.getBool('FullBlast') and self.bamForker()):
                if self.getBool('DotPlots'): self.pngDotPlots()
                if self.getBool('DisMat'): self.disMatrixOut()
                if self.getBool('CombinedFas') and self.getBool('FasOut'): self.combinedFas()
                return True
        except SystemExit: sys.exit()
        except: self.errorLog('Error in GABLAM()',printerror=True,quitchoice=True)
        return False
#########################################################################################################################
    def bamForker(self):    ### Main GABLAM Forker Method = Controls GABLAMMing of sequences
        '''Main GABLAM Forker Method = Controls GABLAMMing of sequences. Returns True/False.'''
        try:
            ## Forks ##
            forkx = int(self.stat['Forks'])   # Number of forks to have running at one time
            if self.opt['Win32'] or forkx < 1: self.opt['NoForks'] = True
            killforks = int(self.stat['KillForks']) # Time in seconds to wait after main thread has apparently finished

            ### Setup ###
            seqlist = self.obj['SeqList']
            startfrom = self.info['StartFrom']
            qseqx = 0   # Counter for sequences
            if self.opt['MemSaver']:
                SEQFILE = open(self.info['QueryDB'], 'r')
                lastline = 'Start'
            readingseq = True   # Still reading in sequences from file
            activedic = {}  # Dictionary of active forks
            blastforks = []  # List of active BLAST fork PIDs
            readforks = []  # List of active results reading fork PIDs

            ### Master Loop ###
            while readingseq or len(blastforks) or len(readforks) or activedic:
                prevtime = time.time()  # Used to check for background sleeping

                ## BLAST forks ## (Or whole lot if noforks) ##
                while readingseq and (len(blastforks) < forkx or self.opt['NoForks']):     
                    killtime = time.time()  # Reset killtime - still doing stuff
                    nextseq = None
                    if self.opt['MemSaver']:
                        (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                        seqlist.seq = [nextseq]
                    elif qseqx < seqlist.seqNum():
                        nextseq = seqlist.seq[qseqx]
                    if not nextseq: # Finished Reading!
                        readingseq = False
                        break
                    qseqx +=1 

                    if startfrom:
                        if nextseq.info['Name'].find('%s ' % startfrom) >= 0: startfrom = None
                        else: continue

                    new_fork_id = '%s.%s' % (rje.fileSafeString(nextseq.info['AccNum']),string.split(self.info['SearchDB'],os.sep)[-1])
                    new_fork_id = string.replace(new_fork_id,'.fas','')
                    self.dict['NameMap'][nextseq.info['AccNum']] = self.dict['NameMap'][nextseq.shortName()] = nextseq.info['Name']
                    if self.opt['NoForks']:     ## BLAST *and* read results ##
                        self.verbose(2,3,'No forks: Straight processing of %s' % new_fork_id)
                        self.bamBLAST(seqlist,nextseq,new_fork_id)
                        self.bamRead(new_fork_id)
                    else:   ## BLAST Forks ##
                        newpid = os.fork() 
                        if newpid == 0: # child
                            self.stat['Interactive'] = -1
                            self.log.stat['Interactive'] = -1
                            self.verbose(2,3,'=> BLAST Forked %s' % new_fork_id)
                            self.log.info['LogFile'] = '%s.log' % new_fork_id
                            self.bamBLAST(seqlist,nextseq,new_fork_id)
                            sys.exit()    # Exit process 
                        elif newpid == -1: # error
                            self.log.errorLog('Problem forking %s.' % new_fork_id)
                        else: 
                            activedic[new_fork_id] = newpid     # This controls BLAST read forks
                            blastforks.append(newpid)           # Add fork to list 

                ## Monitor and remove finished BLAST forks ##
                forklist = self._activeForks(blastforks)
                for key in activedic.keys():
                    if activedic[key] in blastforks and activedic[key] not in forklist:   # Fork finished!
                        activedic[key] = 0
                if len(forklist) != len(blastforks):
                    self.verbose(2,2,' => %d of %d BLAST forks finished!' % (len(blastforks) - len(forklist),len(blastforks)),1)
                blastforks = forklist[0:]

                ## Read Forks ##
                while 0 in activedic.values() and len(readforks) < forkx:   
                    for key in activedic.keys():
                        if activedic[key] == 0:   # New readfork
                            new_fork_id = key
                            break
                    newpid = os.fork() 
                    if newpid == 0: # child
                        self.stat['Interactive'] = -1
                        self.log.stat['Interactive'] = -1
                        self.verbose(2,3,'=> Read Forked %s' % new_fork_id)
                        self.log.info['LogFile'] = '%s.log' % new_fork_id
                        self.bamRead(new_fork_id)
                        sys.exit()    # Exit process 
                    elif newpid == -1: # error
                        self.log.errorLog('Problem forking %s.' % new_fork_id)
                    else: 
                        activedic[new_fork_id] = newpid
                        readforks.append(newpid)    # Add fork to list 
                    
                ## Monitor and remove finished Read forks ##
                forklist = self._activeForks(readforks)
                for key in activedic.keys():
                    if activedic[key] in readforks and activedic[key] not in forklist:   # Fork finished!
                        finito = activedic.pop(key)
                        if self.opt['FullRes']: rje.fileTransfer(fromfile='%s.gablam' % key,tofile=self.info['GablamOut'])
                        if self.opt['HitSum']: rje.fileTransfer(fromfile='%s.hitsum' % key,tofile=self.info['HitSumOut'])
                        if self.opt['Local']: rje.fileTransfer(fromfile='%s.local' % key,tofile=self.info['LocalOut'])
                        rje.fileTransfer(fromfile='%s.log' % key,tofile=self.log.info['LogFile'])
                        self.verbose(2,2,'%s.gablam and %s.log copied and deleted.' % (key,key))
                if len(forklist) != len(readforks):
                    self.verbose(2,2,' => %d of %d Read forks finished!' % (len(readforks) - len(forklist),len(readforks)),1)
                readforks = forklist[0:]

                ## Look for eternal hanging of forks ##
                if time.time() - killtime > killforks:
                    self.verbose(0,1,'\n%d seconds of main program inactivity. %d forks still active!' % (killforks,len(activedic)),1)
                    for fork in activedic.keys(): self.verbose(0,2,' => Fork PID %d still Active!' % (activedic[fork]),1)
                    if (time.time() - prevtime) > killforks:
                        if self.i() >= 0 and rje.yesNo('Possible background sleep detected. Kill anyway?',default='N'): break
                    elif self.i() < 0 or rje.yesNo('Kill?'): break   #!# killing options
                    killtime = time.time()

            ### Finish ###
            try: SEQFILE.close()
            except: self.log.printLog('#LOG','No SEQFILE to close in bamForker(). Everything should be fine, though. I hope! :o)')
            if len(activedic) > 0:
                #self.errorLog('%d Forks still active after %d seconds of main program inactivity' % (len(activedic),killforks))
                raise ValueError('%d Forks still active after %d seconds of main program inactivity' % (len(activedic),killforks))
            else: self.verbose(1,1,'Forks have finished.',2)
            return True

        except SystemExit: sys.exit()
        except:
            self.errorLog('Error in bamForker()',printerror=True,quitchoice=True)
            return False   
#########################################################################################################################
    def bamBLAST(self,seqlist,nextseq,new_fork_id):      ### Performs BLAST of single sequence
        '''
        Performs BLAST of single sequence.
        >> seqlist:rje_seq.SeqList object containing one sequence that needs to be BLASTed.
        >> nextseq:Sequence Object of sequence to search
        >> new_fork_id:str = Base of output name for file *.blast
        '''
        try:
            if self.opt['QryAcc']: seqlist.saveFasta(name='AccNum',seqs=[nextseq],seqfile='%s.qry' % new_fork_id,screen=self.v()>1)
            else: seqlist.saveFasta(name='short',seqs=[nextseq],seqfile='%s.qry' % new_fork_id,screen=self.v()>1)
            blastcmd = ['blastf=F','blastv=%d' % self.stat['DBSize'],'blastb=%d' % self.stat['DBSize'],'ignoredate=T'] + self.cmd_list + ['i=-1']
            if not self.opt['AlnStats']: blastcmd += ['blastb=0']
            if self.getBool('FragFas'): blastcmd = ['gablamfrag=100'] + blastcmd
            blast = rje_blast.blastObj(self.log,blastcmd,type='dev')
            #else: blast = rje_blast.blastObj(self.log,blastcmd,type='Old')
            #self.deBug(blast.cmd_list)
            blast.setInfo({'Name':'%s%s.blast' % (self.getStr('BlastDir'),new_fork_id),
                           'DBase':self.info['SearchDB'],'InFile':'%s.qry' % new_fork_id})
            if self.opt['DeBug']: blast.setInt({'Interactive':1})
            #if seqlist.info['Type'] in ['RNA','DNA']: blast.info['Type'] = 'blastx'
            #else: blast.info['Type'] = 'blastp'
            blast.blast(log=self.opt['DeBug'],use_existing=not self.force())
        except: self.log.errorLog('Error in bamBLAST()',printerror=True,quitchoice=False)
#########################################################################################################################
    def bamRead(self,new_fork_id):      ### Reads BLAST of single sequence and outputs results
        '''
        Performs BLAST of single sequence.
        >> new_fork_id:str = Base of output name for file *.blast
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blastres = '%s%s.blast' % (self.getStr('BlastDir'),new_fork_id)
            ## ~ [0a] ~ BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#BLAST','Reading %s ...' % blastres,screen=self.v()>0)
            if not rje.checkForFile(blastres): return self.errorLog('%s.blast Missing!' % (new_fork_id),printerror=False)
            if rje.checkForFile('%s.qry' % new_fork_id): os.unlink('%s.qry' % new_fork_id)
            #blast = rje_blast.BLASTRun(log=self.log,cmd_list=self.cmd_list+['i=-1'])
            blastcmd = self.cmd_list+['i=-1']
            blast = rje_blast.blastObj(self.log,blastcmd,type='dev')
            if self.opt['DeBug']: blast.setInt({'Interactive':1})
            blast.setStr({'Name':blastres,'DBase':self.info['SearchDB'],'InFile':'%s.qry' % new_fork_id})
            ## ~ [0b] ~ Output Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('NoForks'): self.setStr({'BamFile':self.info['GablamOut'],'SumFile':self.info['HitSumOut'],'LocFile':self.info['LocalOut']})
            else: self.setStr({'BamFile':'%s.gablam' % (new_fork_id),'SumFile':'%s.hitsum' % (new_fork_id),'LocFile':'%s.local' % (new_fork_id)})
            ### ~ [1] ~ Read BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.readBLAST(resfile=blastres,clear=True,gablam=True,local=self.opt['Local'],unlink=not self.getBool('KeepBlast'))
            if 'rje_blast_V1' in '%s' % blast:
                search = blast.search[0]
                self.searchGABLAMLegacy(search,blast)
            else: self.searchGABLAM(blast)
        except: self.errorLog('Error in bamRead(%s.blast)' % new_fork_id,printerror=True,quitchoice=False)
#########################################################################################################################
    def fullBlastGABLAM(self):  ### Performs BLAST of all sequences and then GABLAMFromBLASTRes
        '''Performs BLAST of all sequences and then GABLAMFromBLASTRes.'''
        try:### [0] Setup BLAST search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList']
            if self.getBool('QryAcc'):
                seqfile = '%s.acc.fas' % rje.baseFile(self.getStr('QueryDB'),strip_path=False)
                #!# Consider adding self.force() with additional bool switch to cover programs with parallel runs (e.g. REVERT)
                if not rje.exists(seqfile): seqlist.saveSeq(reformat='accfas',seqfile=seqfile)
            else: seqfile = self.getStr('QueryDB')
            blastfile = '%s%s.blast' % (self.getStr('BlastDir'),os.path.basename(self.basefile()))
            blastcmd = ['blastf=F','blastv=%d' % self.stat['DBSize'],'blastb=%d' % self.stat['DBSize'],'ignoredate=T'] + self.cmd_list + ['i=-1']
            if not self.opt['AlnStats']: blastcmd += ['blastb=0']
            if self.getBool('FragFas'): blastcmd = ['gablamfrag=100'] + blastcmd
            if self.getInt('Forks') > 0: blastcmd += ['blasta=%s' % self.getInt('Forks')]
            blast = rje_blast.blastObj(self.log,blastcmd,type='dev')
            blast.setInfo({'Name':blastfile,'DBase':self.getStr('SearchDB'),'InFile':seqfile})
            if self.debugging(): blast.setInt({'Interactive':1})
            ### [1] Run BLAST Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.blast(log=True,use_existing=not self.force())
            self.setStr({'BlastRes':blastfile})
            ### [2] Perform GABLAM Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.GABLAMFromBLASTRes()
        except: self.log.errorLog('Error in bamBLAST()',printerror=True,quitchoice=False)
#########################################################################################################################
    def GABLAMFromBLASTRes(self):   ### Produces output from BLAST results file alone
        '''Produces output from BLAST results file alone.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RES','Generating GABLAM output from full BLAST results: %s' % self.getStr('BlastRes'))
            if self.basefile().lower() in ['','none']: self.basefile(self.getStr('BlastRes'))
            blast = rje_blast.BLASTRun(log=self.log,cmd_list=['blastf=F']+self.cmd_list)
            self.setupResultsFiles(saveheads=False)
            self.setStr({'BamFile':self.getStr('GablamOut'),'SumFile':self.getStr('HitSumOut'),
                         'LocFile':self.getStr('LocalOut'),'DisFile':self.getStr('DisMatOut')})
            ### ~ [1] Read Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Local'): self.setBool({'LocalAlnFas':False})
            blast.readBLAST(resfile=self.getStr('BlastRes'),clear=True,gablam=True,local=self.getBool('Local'),keepaln=self.getBool('LocalAlnFas'))
            if self.dev():  # Not sure whether these can be used
                for table in blast.db().list['Tables']: table.saveToFile()
            blast.formatTables()
            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'rje_blast_V1' in '%s' % blast:
                for search in blast.search: self.searchGABLAMLegacy(search,blast)
            else: return self.multiGABLAM(blast)
            return True
        except:
            self.errorLog('Error in GABLAMFromBLASTRes()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def multiGABLAM(self,blast):  ### Performs GABLAM manipulation and output for a single search
        '''
        Performs GABLAM manipulation and output for a mutliple searches. Will ultimately replace single search GABLAM
        with a switch to perform the full BLAST search first and then use multiGABLAM reading of the the output. GABLAM
        Output file names should be in self.info 'FullFile' and 'SumFile' as set by bamRead() and GABLAMFromBLASTRes().
        Align is performed if desired with searchAlign().
        >> blast:BLASTRun object needed for Fasta output
        '''
        try:### [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ores = {'GABLAM':'Aln','GABLAMO':'OrderedAln'}
            if self.getBool('GlobID'): self.warnLog('ALIGN Gobal identity not yet available with BLAST+','align',suppress=True)
            ## blast.db() now contains:
            ## Table('Run',['Run','Type','E-Value','DBase','InFile','BLASTCmd','Complexity Filter','Composition Statistics','SoftMask','GappedBLAST','OneLine','HitAln','DBLen','DBNum'],['Run'])
            ## Table('Search',['Query','Length','Hits']+['Coverage','Identity','Positives'],['Query'])
            sdb = blast.db('Search')
            ## Table('Hit',['Query','Rank','Hit','Description','BitScore','E-Value','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Query','Hit'])
            hdb = blast.db('Hit')
            ## Table('Local',['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],['Query','Hit','AlnID'])
            ldb = blast.db('Local')
            ## Table('GABLAM',['Query','Hit','QryHit','GABLAM Start','GABLAM End','GABLAM Dirn','GABLAM Len','GABLAM Sim','GABLAM ID','GABLAM Frag','GABLAMO Start','GABLAMO End','GABLAMO Dirn','GABLAMO Len','GABLAMO Sim','GABLAMO ID','GABLAMO Frag'],['Run','Query','Hit','QryHit'])
            gdb = blast.db('GABLAM')
            resdb = blast.db().addEmptyTable('Full',self.list['FullResHeaders'],['Qry','Hit'])
            ## [0a] GABLAMCut settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qfields = ['Coverage','Identity','Positives']
            gstat = gtype = None
            gcut = self.getNum('GABLAMCut')
            if gcut > 0:
                gstat = self.getStr('CutStat')
                if self.getBool('QAssemble') and gstat in qfields + ['Local']: gtype = 'qassemble'
                else:
                    gstat = string.replace(gstat,'Ordered','O')
                    gstat = string.replace(gstat,'Aln',' ')
                    gstat = 'GABLAM' + gstat
                    if gstat not in gdb.fields():
                        self.errorLog('GABLAM Cut-off stat "%s" not recognised. Should be generic (use with cutfocus for Query/Hit).' % gstat,False,False)
                        self.debug(gstat)
                        if self.i() < 0 or rje.yesNo('Proceed without GABLAMCut?'):
                            self.setStat({'GABLAMCut':0})
                        else: return False
                    gtype = self.getStrLC('CutFocus')
                    if gtype not in ['query','qry','hit','either','both']:
                        self.errorLog('GABLAM Cut-off type "%s" not recognised. Will use "Either"' % self.getStr('CutType'),False,False)
                        self.setStr({'CutFocus':'Either'})
                        gtype = 'either'

            ### [1] Filter Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## [1a] Remove Self Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getBool('SelfHit'):
                for table in [hdb,ldb]:
                    hx = 0
                    for entry in table.entries()[0:]:
                        (Qry,Hit) = (entry['Query'],entry['Hit'])
                        self_hit = Qry == Hit
                        if self.getBool('QryAcc'): self_hit = self_hit or string.split(Hit,'__')[-1] == Qry
                        if self_hit: table.dropEntry(entry); hx += 1    #? Also cleanup Local and GABLAM datasets?
                    self.printLog('\r#SELF','%s %s self-hits removed: %s remain' % (rje.iStr(hx),table.name(),rje.iStr(hdb.entryNum())))
            ## [1b] Apply GABLAM Cutoffs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if gcut > 0 and gtype == 'qassemble':
                qdump = []
                for sentry in sdb.entries()[0:]:
                    if sentry['Hits'] == 0: continue
                    if gstat == 'Local': cut = float(sentry['Coverage']) <=0 or (100.0 * float(sentry['Identity']) / float(sentry['Coverage'])) < gcut
                    else: cut = 100.0 * int(sentry[gstat]) / sentry['Length'] < gcut
                    if cut:
                        for field in qfields + ['Hits']: sentry[field] = 0
                        qdump.append(sentry['Query'])
                if qdump:
                    hdb.dropEntriesDirect('Query',qdump)
                    ldb.dropEntriesDirect('Query',qdump)
                    gdb.dropEntriesDirect('Query',qdump)
                self.printLog('#CUT','Results for %s of %s queries rejected (< %.1f%% QAssemble.%s)' % (rje.iLen(qdump),rje.iStr(sdb.entryNum()),gcut,gstat))
            if self.getBool('AlnStats'):
                ex = 0.0; etot = hdb.entryNum(); fx = 0
                for entry in hdb.entries()[0:]:
                    self.progLog('\r#GABLAM','Processing GABLAM Data: %.1f%%' % (ex/etot)); ex += 100.0
                    (Qry,Hit) = (entry['Query'],entry['Hit'])
                    plen = {'Hit':entry['Length'],'Query':sdb.data(Qry)['Length']}
                    gdict = blast.gablamData(Hit,Qry)
                    cut = False
                    if self.getNum('GABLAMCut') > 0 and gtype != 'qassemble':
                        cutq = 100.0 * float(gdict['Query'][gstat]) / plen['Query']
                        cuth = 100.0 * float(gdict['Hit'][gstat]) / plen['Hit']
                        if gtype in ['query','qry','both'] and cutq < gcut: cut = True
                        elif gtype in ['hit','both'] and cuth < gcut: cut = True
                        elif gtype in ['either'] and cuth < gcut and cutq < gcut: cut = True
                        self.verbose(2,2,'%.2f%% & %.2f%% vs %.2f%% %s (%s) => Cut=%s' % (cutq,cuth,gcut,gstat,gtype,cut),1)
                    #self.debug(entry)
                    cut = cut or float(entry['E-Value']) > self.getNum('BlastE')
                    if cut: hdb.dropEntry(entry); fx += 1    #? Also cleanup Local and GABLAM datasets?
                    else:   # Add to full results table
                        resdict = {}
                        for qh in ['Query','Hit']:
                            rqh = {'Query':'Qry','Hit':'Hit'}[qh]
                            for order in ['GABLAM','GABLAMO']:
                                if self.getStrUC('OutStats') not in [order,'ALL']: continue
                                for stat in ['Len', 'ID', 'Sim']:
                                    dstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order],stat)
                                    if self.getBool('PercRes'):
                                        if plen[qh] > 0: resdict[rstat] = '%.2f' % (100.0 * float(gdict[qh][dstat]) / plen[qh])
                                        else: resdict[rstat] = '0.0'
                                    else: resdict[rstat] = gdict[qh][dstat]
                                for stat in ['Dirn','Start','End']:
                                    dstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order][:-3],stat)
                                    resdict[rstat] = gdict[qh][dstat]
                        fentry = rje.combineDict(resdict,entry)
                        fentry['Qry'] = fentry.pop('Query')
                        #self.debug(fentry)
                        fentry['Score'] = fentry.pop('BitScore')
                        fentry['EVal'] = fentry.pop('E-Value')
                        fentry['QryLen'] = plen['Query']
                        fentry['HitLen'] = plen['Hit']
                        resdb.addEntry(fentry)
                self.printLog('\r#GABLAM','Processing GABLAM Data complete.',log=False)
                if gcut > 0 and gtype != 'qassemble':
                    self.printLog('\r#FILT','%s of %s hits filtered by %s cutoff (%s): %s remain' % (rje.iStr(fx),rje.iStr(etot),gstat,gtype,rje.iStr(hdb.entryNum())))

            ### [2] Output results tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## [2a] Full GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # outstats=X      : Whether to output just GABLAM, GABLAMO or All [All]
            #Qry	Hit	Rank	Score	EVal	QryLen	HitLen	Qry_AlnLen	Qry_AlnID	Qry_AlnSim	Qry_Dirn	Qry_Start	Qry_End	Qry_OrderedAlnLen	Qry_OrderedAlnID	Qry_OrderedAlnSim	Qry_OrderedDirn	Qry_OrderedStart	Qry_OrderedEnd	Hit_AlnLen	Hit_AlnID	Hit_AlnSim	Hit_Dirn	Hit_Start	Hit_End	Hit_OrderedAlnLen	Hit_OrderedAlnID	Hit_OrderedAlnSim	Hit_OrderedDirn	Hit_OrderedStart	Hit_OrderedEnd
            if self.getBool('FullRes'):
                ex = resdb.entryNum()
                keydict = {}
                for ekey in resdb.dataKeys():
                    entry = resdb.data(ekey)
                    keysort = '%s.%s' % (entry['Qry'],rje.preZero(entry['Rank'],ex))
                    if keysort not in keydict: keydict[keysort] = []
                    keydict[keysort].append(ekey)
                savekeys = []
                for ksort in rje.sortKeys(keydict): savekeys += keydict[ksort]
                resdb.saveToFile(self.getStr('BamFile'),savekeys=savekeys)
            ## [2b] Local Stat Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Local'):
                ldb.dropEntries(['Length<%d' % self.getInt('LocalMin')],logtxt='LocalMin',log=True)
                if self.getBool('LocalAlnFas') and ldb.entries():
                    lalnfile = self.dict['Output']['localalnfas'] = '%s.local.aln.fas' % self.basefile()
                    rje.backup(self,lalnfile)
                    self.progLog('\r#FAS','Local alignment output to %s...' % lalnfile)
                    LALN = open(lalnfile,'w')
                    for lkey in rje.sortKeys(ldb.data()):
                        lentry = ldb.data(lkey)
                        try: qid = '%.1f%%' % (100.0 * int(lentry['Identity']) / (int(lentry['QryEnd']) - int(lentry['QryStart']) + 1))
                        except: qid = 'ERROR'
                        try: qpos = '%.1f%%' % (100.0 * int(lentry['Positives']) / (int(lentry['QryEnd']) - int(lentry['QryStart']) + 1))
                        except: qpos = 'ERROR'
                        LALN.write('>%s|%s|%s|Qry:%s-%s Identity:%s Positives:%s\n' % (lentry['Query'],lentry['Hit'],lentry['AlnID'],lentry['QryStart'],lentry['QryEnd'],qid,qpos))
                        LALN.write('%s\n' % lentry['QrySeq'])
                        LALN.write('>%s|%s|%s|Aln:%s BitScore:%s Expect:%s Identity:%s Positives:%s\n' % (lentry['Query'],lentry['Hit'],lentry['AlnID'],lentry['Length'],lentry['BitScore'],lentry['Expect'],lentry['Identity'],lentry['Positives']))
                        LALN.write('%s\n' % lentry['AlnSeq'])
                        try: hid = '%.1f%%' % (100.0 * int(lentry['Identity']) / (int(lentry['SbjEnd']) - int(lentry['SbjStart']) + 1))
                        except: hid = 'ERROR'
                        try: hpos = '%.1f%%' % (100.0 * int(lentry['Positives']) / (int(lentry['SbjEnd']) - int(lentry['SbjStart']) + 1))
                        except: hpos = 'ERROR'
                        LALN.write('>%s|%s|%s|Hit:%s-%s Identity:%s Positives:%s\n' % (lentry['Query'],lentry['Hit'],lentry['AlnID'],lentry['SbjStart'],lentry['SbjEnd'],hid,hpos))
                        LALN.write('%s\n' % lentry['SbjSeq'])
                    LALN.close()
                    self.printLog('\r#FAS','Local alignment output to %s complete.' % lalnfile)
                ldb.dropFields(['QrySeq','SbjSeq','AlnSeq'])
                ldb.renameField('Query','Qry')
                ldb.renameField('AlnID','AlnNum')
                ldb.saveToFile(self.getStr('LocFile'))
            ## [2c] Hit Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('HitSum'):
                resdb.setStr({'Name':'HitSum'})
                headlist = ['Qry','Hit','Score','EVal']
                if self.getBool('QAssemble'): headlist += ['Length','Coverage','Identity','Positives']
                resdb.keepFields(headlist,log=False)
                if self.getBool('SelfHit') and not self.getBool('SelfSum'):
                    for entry in resdb.entries()[0:]:
                        if entry['Qry'] == entry['Hit']: resdb.dropEntry(entry)
                        elif self.getBool('QryAcc') and string.split(entry['Hit'],'__')[-1] == entry['Qry']: resdb.dropEntry(entry)
                resdb.addField('HitNum',evalue=1,after='Qry')
                resdb.compress(['Qry'],{'HitNum':'sum','Score':'max','EVal':'min'})
                for field in headlist:
                    if field not in resdb.fields(): resdb.addField(field)
                resdb.dropField('Hit')
                resdb.renameField('Score','MaxScore')
                nullx = 0
                for entry in sdb.entries():
                    if resdb.data(entry['Query']):
                        qentry = resdb.data(entry['Query'])
                        if self.getBool('QAssemble'):
                            for field in ['Length','Coverage','Identity','Positives']: qentry[field] = entry[field]
                    else:
                        qentry = {'Qry':entry['Query'],'HitNum':entry['Hits'],'MaxScore':entry['MaxScore'],'EVal':entry['TopE']}
                        if self.getBool('QAssemble'):
                            for field in ['Length','Coverage','Identity','Positives']: qentry[field] = entry[field]
                        resdb.addEntry(qentry); nullx += 1
                self.printLog('#NULL','%s BLAST searches w/o hits added to HitSum.' % rje.iStr(nullx))
                resdb.saveToFile(self.info['SumFile'])

            ### ~ [3] Output Additional Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Note that hdb has the full list of Query,Hit pairs to process
            ## ~ [3a] Fasta Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug('FasOut')
            addflanks = self.getInt('AddFlanks')
            gabfrag = max(blast.getInt('GablamFrag') - addflanks,0)     # Want to join fragments within gabfrag distance
            if self.getBool('FasOut') and hdb.entryNum():               # Data to output
                outfas = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('SearchDB'),'autoload=T','autofilter=F','mode=file'])
                qx = 0; qtot = len(hdb.index('Query')); fasx = 0
                for Qry in hdb.index('Query'):
                    self.progLog('\r#FAS','Generating fasta output: %.1f%%' % (qx/qtot)); qx += 100.0
                    bfas = os.path.join(rje.makePath(self.getStr('FasDir')),'%s.fas' % os.path.basename(Qry))
                    hitlist = hdb.indexDataList('Query',Qry,'Hit')
                    if not hitlist: continue
                    if self.getBool('FragFas'):
                        fragfas = []    # Tuple of (name,sequence) for fragments
                        blast.setOpt({'DeBug':self.getBool('DeBug')})
                        hitseq = blast.hitToSeq(outfas,hitlist)
                        hx = 0.0; htot = len(hitlist); hfx = 0
                        #self.deBug(hitseq)
                        for hit in hitlist:
                            #self.bugPrint(hit.dict['GABLAM']['Hit'])
                            if hit not in hitseq: continue
                            gdata = blast.gablamData(hit,Qry)
                            ftot = len(gdata['Hit']['GABLAM Frag'])
                            if self.v()>1: self.progLog('\r#FRAG','Generating FragFas output: %.2f%%' % (hx/htot)); hx += (100.0 / (ftot+1))
                            #if addflanks > 0:
                            gfrag = []
                            for (fragstart,fragend) in gdata['Hit']['GABLAM Frag']:
                                fragstart = max(0,fragstart-addflanks)
                                fragend += addflanks
                                if not gfrag or fragstart > (gfrag[-1][1] + 1 + gabfrag): gfrag.append((fragstart,fragend))
                                else: gfrag[-1] = (gfrag[-1][0],fragend)
                            if len(gfrag) != ftot:
                                if addflanks > 0: self.printLog('\r#FLANK','%s nt flanks added. %d %s fragments merged to %d.' % (rje.iStr(addflanks),ftot,hit,len(gfrag)))
                                else: self.printLog('\r#FLANK','%d %s fragments merged to %d.' % (ftot,hit,len(gfrag)))
                                ftot = len(gfrag)
                            #else: gfrag = gdata['Hit']['GABLAM Frag']
                            seq = hitseq[hit]
                            (fullseqname,fullseq) = outfas.getSeq(seq,format='tuple')
                            seqlen = len(fullseq)
                            for fx in range(ftot):
                                if self.v()>1: self.progLog('\r#FRAG','Generating FragFas output: %.2f%%' % (hx/htot)); hx += (100.0 / (ftot+1))
                                #self.bugPrint(hit.dict['GABLAM']['Hit']['GABLAM Frag'][fx])
                                (fragstart,fragend) = gfrag[fx]
                                fragstart = max(1,fragstart+1)
                                fragend = min(fragend+1,seqlen)
                                sequence = fullseq[fragstart-1:fragend]
                                #seq = hitseq[hit]
                                #(fullseqname,sequence) = outfas.getSeqFrag(seq,fragstart+1,fragend+1)
                                seqname = string.split(fullseqname)
                                seqname[0] = '%s-%s.%s' % (seqname[0],rje.preZero(fragstart,seqlen),rje.preZero(fragend,seqlen))
                                if len(seqname) == 1: seqname.append('No description')
                                seqname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                                seqname.insert(1,'%s GABLAM Hit:' % Qry)
                                seqname = string.join(seqname)
                                fragfas.append((seqname,sequence)); hfx += 1
                                #self.debug(fullseqname)
                                #self.debug(Qry)
                                #self.debug(seqname)
                        self.printLog('\r#FRAG','Generated %s fragments for %s %s GABLAM hits => %s.' % (rje.iStr(hfx),rje.iStr(htot),Qry,bfas),screen=self.v()>1)
                        if fragfas:
                            outfas.setStr({'SeqMode':'list'})
                            outfas.saveSeq(fragfas,bfas,reformat='fas',append=None,log=False,screen=self.dev() or self.v() > 1,backup=True); fasx += 1
                            outfas.setStr({'SeqMode':'file'})
                        else: self.printLog('#ERR','0 sequences from %d hits output to %s!' % (len(hitlist),bfas))
                        #self.deBug(fragfas)
                    else: blast.hitToSeq(outfas,hitlist,filename=bfas,appendfile=False); fasx += 1
                self.progLog('\r#FAS','Generated %s fasta files for %s queries.' % (rje.iStr(fasx),rje.iStr(qtot)))

            return True
        except:
            self.errorLog('GABLAM.multiGABLAM() is banjaxed',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def searchGABLAM(self,blast):  ### Performs GABLAM manipulation and output for a single search
        '''
        Performs GABLAM manipulation and output for a single search. Output file names should be in self.info 'FullFile'
        and 'SumFile' as set by bamRead() and GABLAMFromBLASTRes(). Align is performed if desired with searchAlign().
        >> blast:BLASTRun object needed for Fasta output
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## blast.db() now contains:
            ## Table('Run',['Run','Type','E-Value','DBase','InFile','BLASTCmd','Complexity Filter','Composition Statistics','SoftMask','GappedBLAST','OneLine','HitAln','DBLen','DBNum'],['Run'])
            ## Table('Search',['Run','Query','Length','Hits'],['Run','Query'])
            ## Table('Hit',['Run','Query','Rank','Hit','Description','BitScore','E-Value','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Run','Query','Hit'])
            ## Table('Local',['Run','Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],['Run','Query','Hit','AlnID'])
            ## Table('GABLAM',['Run','Query','Hit','QryHit','GABLAM Start','GABLAM End','GABLAM Dirn','GABLAM Len','GABLAM Sim','GABLAM ID','GABLAM Frag','GABLAMO Start','GABLAMO End','GABLAMO Dirn','GABLAMO Len','GABLAMO Sim','GABLAMO ID','GABLAMO Frag'],['Run','Query','Hit','QryHit'])
            delimit = rje.getDelimit(self.cmd_list)
            resdict = {}
            ores = {'GABLAM':'Aln','GABLAMO':'OrderedAln'}
            sentry = blast.db('Search').entries()[0]
            Qry = sentry['Query']
            QryLen = sentry['Length']
            Rank = 1
            last_score = 0
            max_score = 0
            max_eval = 1000
            hit_num = 0
            hitlist = []
            ## ~ [0a] Qassemble filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qfields = ['Coverage','Identity','Positives']
            if self.getBool('Qassemble') and self.getStr('CutStat') in qfields + ['Local']:
                if self.getStr('CutStat') == 'Local': cut = float(sentry['Coverage']) <=0 or (100.0 * float(sentry['Identity']) / float(sentry['Coverage'])) < self.getNum('GABLAMCut')
                else: cut = 100.0 * int(sentry[self.getStr('CutStat')]) / sentry['Length'] < self.getNum('GABLAMCut')
                if cut:
                    blast.db('Hit').clear()
                    for field in qfields: sentry[field] = 0
            ### ~ [1] Process Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.db('Hit').index('BitScore',log=False)
            blasthits = len(blast.db('Hit').sortedEntries('BitScore',reverse=True)) > 0
            for hentry in blast.db('Hit').sortedEntries('BitScore',reverse=True):
                Hit = hentry['Hit']
                try:## ~ [1a] Basics Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #self.bugPrint(Qry)
                    #self.bugPrint(Hit)
                    self_hit = Qry == Hit
                    if self.getBool('QryAcc'): self_hit = self_hit or string.split(Hit,'__')[-1] == Qry
                    #self_hit = Qry.find(Hit) >= 0 or (self.opt['QryAcc'] and Hit.find(Qry) >= 0)
                    #self.bugPrint('%s vs %s self hit: %s' % (Qry,Hit,self_hit))
                    if self_hit and not self.opt['SelfHit']:
                        blast.db('Hit').dropEntry(hentry)
                        continue
                    Score = hentry['BitScore']
                    if Score < last_score: Rank += 1
                    last_score = Score
                    EValue = hentry['E-Value']
                    HitLen = hentry['Length']
                    resdict[Hit] = {'Qry':Qry,'Hit':Hit,'Rank':Rank,'Score':Score,'EVal':EValue,'QryLen':QryLen,'HitLen':HitLen}
                    hitlist.append(Hit)
                    ## ~ [1b] Complex Align Stuff (GABLAM) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    plen = {'Query':QryLen,'Hit':HitLen}                    
                    if self.opt['AlnStats']:
                        gdict = blast.gablamData(Hit,Qry)
                        #self.deBug(gdict)
                        for qh in ['Query','Hit']:
                            rqh = qh
                            if qh == 'Query': rqh = 'Qry'
                            if plen[qh] <= 0:
                                self.log.errorLog('Error in bamRead(%s, Hit %s): %s len < 1!' % (Qry,Hit,qh),printerror=False)
                                if qh == 'Hit': self.log.errorLog('Check for multiple occurrences of %s in search database!' % (Hit),printerror=False)
                                self.bugPrint(gdict)
                            for order in ['GABLAM','GABLAMO']:
                                for stat in ['Len', 'ID', 'Sim']:
                                    gstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order],stat)
                                    if self.opt['PercRes']:
                                        if plen[qh] > 0: resdict[Hit][rstat] = '%.2f' % (100.0 * float(gdict[qh][gstat]) / plen[qh])
                                        else: resdict[Hit][rstat] = '0.0'
                                    else: resdict[Hit][rstat] = gdict[qh][gstat]
                                for stat in ['Dirn','Start','End']:
                                    gstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order][:-3],stat)
                                    resdict[Hit][rstat] = gdict[qh][gstat]

                    ## Check GABLAM Cut ##
                    if self.stat['GABLAMCut'] > 0:
                        gcut = self.stat['GABLAMCut'] 
                        gstat = self.info['CutStat']
                        gstat = string.replace(gstat,'Ordered','O')
                        gstat = string.replace(gstat,'Aln',' ')
                        gstat = 'GABLAM' + gstat
                        if gstat not in gdict['Query'].keys():
                            print gstat, 'not cool. :o('
                            self.errorLog('GABLAM Cut-off stat "%s" not recognised. Should be generic (use with cutfocus for Query/Hit).' % gstat,False,False)
                            self.stat['GABLAMCut'] = 0
                            self.info['CutStat'] = 'AlnLen'                       
                        gtype = self.info['CutFocus'].lower()
                        if gtype not in ['query','qry','hit','either','both']:
                            self.log.errorLog('GABLAM Cut-off type "%s" not recognised. Will use "Either"' % self.info['CutType'],False,False)
                            self.info['CutFocus'] = 'Either'
                            gtype = 'either'
                        cutq = 100.0 * float(gdict['Query'][gstat]) / plen['Query']
                        cuth = 100.0 * float(gdict['Hit'][gstat]) / plen['Hit']
                        cut = False
                        if gtype in ['query','qry','both'] and cutq < gcut: cut = True
                        elif gtype in ['hit','both'] and cuth < gcut: cut = True
                        elif gtype in ['either'] and cuth < gcut and cutq < gcut: cut = True
                        self.verbose(2,2,'%.2f%% & %.2f%% vs %.2f%% %s (%s) => Cut=%s' % (cutq,cuth,gcut,gstat,gtype,cut),1)
                        if cut:
                            hitlist.remove(Hit)
                            reslist = []
                            continue

                    ## Summary - now happens after GABLAM Cut-off ##
                    if self.opt['SelfSum'] or not self_hit:
                        hit_num += 1
                        if Score > max_score :
                            max_score = Score
                            max_eval = EValue

                except:
                    self.errorLog('Major problem in searchGABLAM(%s, Hit=%s).' % (Qry,Hit))
                    hitlist.remove(Hit)
            if blasthits: self.printLog('#HIT','%s: %s hits following GABLAM filtering.' % (Qry,rje.iLen(hitlist)))

            ### ~ [2] ALIGN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['GlobID']: self.warnLog('ALIGN Gobal identity not yet available with BLAST+','align',suppress=True)
            #resdict = self.searchALIGN(search,resdict)
                    
            ### ~ [3] Output Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Full GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug('FullRes')
            if self.opt['FullRes']:
                for Hit in hitlist:
                    rje.delimitedFileOutput(self,self.info['BamFile'],self.list['FullResHeaders'],delimit,resdict[Hit])
            ## ~ [3b] Local Stat Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug('Local: %s' % self.getBool('Local'))
            if self.getBool('Local'):
                for Hit in hitlist:
                    #self.debug("%s,%s: %s" % (Hit,Qry,blast.localData(Hit,Qry)))
                    for lentry in blast.localData(Hit,Qry):
                        #self.deBug(lentry)
                        if lentry['Length'] < self.stat['LocalMin']: continue
                        rje.delimitedFileOutput(self,self.info['LocFile'],self.list['LocalHeaders'],delimit,rje.combineDict({'Qry':Qry,'Hit':Hit,'AlnNum':lentry['AlnID']},lentry))
            ## ~ [3c] Fasta Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug('FasOut')
            addflanks = self.getInt('AddFlanks')
            gabfrag = max(blast.getInt('GablamFrag') - addflanks,0)     # Want to join fragments within gabfrag distance
            if self.opt['FasOut'] and blast and hitlist:
                bfas = os.path.join(rje.makePath(self.info['FasDir']),'%s.fas' % os.path.basename(Qry))
                #outfas = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                outfas = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('SearchDB'),'autoload=T','autofilter=F','mode=file'])
                newseqlist = True
                if self.getBool('FragFas') and newseqlist:
                    fragfas = []    # Tuple of (name,sequence) for fragments
                    blast.setOpt({'DeBug':self.getBool('DeBug')})
                    hitseq = blast.hitToSeq(outfas,hitlist)
                    hx = 0.0; htot = len(hitlist); hfx = 0
                    #self.deBug(hitseq)
                    for hit in hitlist:
                        #self.bugPrint(hit.dict['GABLAM']['Hit'])
                        if hit not in hitseq: continue
                        gdata = blast.gablamData(hit,Qry)
                        ftot = len(gdata['Hit']['GABLAM Frag'])
                        self.progLog('\r#FRAG','Generating FragFas output: %.2f%%' % (hx/htot)); hx += (100.0 / (ftot+1))
                        #if addflanks > 0:
                        gfrag = []
                        for (fragstart,fragend) in gdata['Hit']['GABLAM Frag']:
                            fragstart = max(0,fragstart-addflanks)
                            fragend += addflanks
                            if not gfrag or fragstart > (gfrag[-1][1] + 1 + gabfrag): gfrag.append((fragstart,fragend))
                            else: gfrag[-1] = (gfrag[-1][0],fragend)
                        if len(gfrag) != ftot:
                            if addflanks > 0: self.printLog('\r#FLANK','%s nt flanks added. %d %s fragments merged to %d.' % (rje.iStr(addflanks),ftot,hit,len(gfrag)))
                            else: self.printLog('\r#FLANK','%d %s fragments merged to %d.' % (ftot,hit,len(gfrag)))
                            ftot = len(gfrag)
                        #else: gfrag = gdata['Hit']['GABLAM Frag']
                        seq = hitseq[hit]
                        (fullseqname,fullseq) = outfas.getSeq(seq,format='tuple')
                        seqlen = len(fullseq)
                        for fx in range(ftot):
                            self.progLog('\r#FRAG','Generating FragFas output: %.2f%%' % (hx/htot)); hx += (100.0 / (ftot+1))
                            #self.bugPrint(hit.dict['GABLAM']['Hit']['GABLAM Frag'][fx])
                            (fragstart,fragend) = gfrag[fx]
                            fragstart = max(1,fragstart+1)
                            fragend = min(fragend+1,seqlen)
                            sequence = fullseq[fragstart-1:fragend]
                            #seq = hitseq[hit]
                            #(fullseqname,sequence) = outfas.getSeqFrag(seq,fragstart+1,fragend+1)
                            seqname = string.split(fullseqname)
                            seqname[0] = '%s-%s.%s' % (seqname[0],rje.preZero(fragstart,seqlen),rje.preZero(fragend,seqlen))
                            if len(seqname) == 1: seqname.append('No description')
                            seqname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                            seqname.insert(1,'%s GABLAM Hit:' % Qry)
                            seqname = string.join(seqname)
                            fragfas.append((seqname,sequence)); hfx += 1
                            #self.debug(fullseqname)
                            #self.debug(Qry)
                            #self.debug(seqname)
                    self.printLog('\r#FRAG','Generated %s fragments for %s GABLAM hits.' % (rje.iStr(hfx),rje.iStr(htot)))
                    if fragfas:
                        outfas.setStr({'SeqMode':'list'})
                        outfas.saveSeq(fragfas,bfas,reformat='fas',append=None,log=True,screen=self.dev() or self.v() > 0,backup=True)
                        outfas.setStr({'SeqMode':'file'})
                    elif hitlist: self.printLog('#ERR','0 sequences from %d hits output to %s!' % (len(hitlist),bfas))
                    #self.deBug(fragfas)
                elif self.getBool('FragFas'):
                    blast.setOpt({'DeBug':self.getBool('DeBug')})
                    hitseq = blast.hitToSeq(outfas,hitlist)
                    #self.deBug(hitseq)
                    outfas.seq = []
                    for hit in hitlist:
                        #self.bugPrint(hit.dict['GABLAM']['Hit'])
                        if hit not in hitseq: continue
                        gdata = blast.gablamData(hit,Qry)
                        ftot = len(gdata['Hit']['GABLAM Frag'])
                        for fx in range(len(gdata['Hit']['GABLAM Frag'])):
                            #self.bugPrint(hit.dict['GABLAM']['Hit']['GABLAM Frag'][fx])
                            (fragstart,fragend) = gdata['Hit']['GABLAM Frag'][fx]
                            seq = hitseq[hit]
                            fullseq = seq.getSequence()
                            seqname = string.split(seq.info['Name'])
                            seqname[0] = '%s-%s.%s' % (seqname[0],rje.preZero(fragstart+1,seq.seqLen()),rje.preZero(fragend+1,seq.seqLen()))
                            seqname.insert(1,'GABLAM Fragment %d of %d (%s - %s)' % (fx+1,ftot,rje.iStr(fragstart+1),rje.iStr(fragend+1)))
                            seqname = string.join(seqname)
                            sequence = fullseq[fragstart:fragend+1]
                            outfas._addSeq(seqname,sequence)
                    if outfas.seq: outfas.saveFasta(seqfile=bfas,append=False); #self.deBug(outfas.seq[0].info)
                    elif hitlist: self.printLog('#ERR','0 sequences from %d hits output to %s!' % (len(hitlist),bfas))
                    #self.deBug(outfas.seq)
                else: blast.hitToSeq(outfas,hitlist,filename=bfas,appendfile=False)
            #elif self.opt['FasOut'] and blast:
            ## ~ [3d] Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['HitSum']:
                sumres = {'Qry':Qry,'HitNum':hit_num,'MaxScore':max_score,'EVal':max_eval}
                if self.getBool('QAssemble'):
                    for field in ['Length','Coverage','Identity','Positives']: sumres[field] = sentry[field]
                #self.bugPrint(sumres)
                rje.delimitedFileOutput(self,self.info['SumFile'],self.list['HitSumHeaders'],delimit,sumres)
                self.printLog('#HIT','%s => %d Database hits.' % (Qry,hit_num),screen=self.v()>1)
        except:
            self.errorLog('GABLAM.searchGABLAM() is banjaxed',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def searchGABLAMLegacy(self,search,blast=None):  ### Performs GABLAM manipulation and output for a single search
        '''
        Performs GABLAM manipulation and output for a single search. Output file names should be in self.info 'FullFile'
        and 'SumFile' as set by bamRead() and GABLAMFromBLASTRes(). Align is performed if desired with searchAlign().
        >> search:BLAST Search object with gablam already performed on hits
        >> blast:BLASTRun object needed for Fasta output
        '''
        try:
            ### Setup ###
            delimit = rje.getDelimit(self.cmd_list)
            resdict = {}
            ores = {'GABLAM':'Aln','GABLAMO':'OrderedAln'}
            Qry = search.info['Name']
            #print 'Search name =', Qry
            QryLen = search.stat['Length']
            Rank = 1
            last_score = 0
            max_score = 0
            max_eval = 1000
            hit_num = 0

            ### Process ###
            for hit in search.hit[0:]:
                ## Basics ##
                Hit = hit.info['Name']
                try:
                    self_hit = Qry.find(Hit) >= 0 or (self.opt['QryAcc'] and Hit.find(Qry) >= 0)

                    #if Qry.find(Hit)  >= 0 or (self.opt['QryAcc'] and Hit.find(Qry) >= 0):    # Self
                    #    self_hit = True
                    #    if self.opt['SelfSum']: hit_num += 1
                    #else:
                    #    self_hit = False
                    #    hit_num += 1
                    if self_hit and not self.opt['SelfHit']:
                        search.hit.remove(hit)
                        continue
                    Score = hit.stat['BitScore']
                    if Score < last_score: Rank = search.hit.index(hit) + 1
                    last_score = Score
                    EValue = hit.stat['E-Value']
                    HitLen = hit.stat['Length']
                    resdict[hit] = {'Qry':Qry,'Hit':Hit,'Rank':Rank,'Score':Score,'EVal':EValue,'QryLen':QryLen,'HitLen':HitLen}
                    #self.deBug(resdict[hit])
                                
                    ## Complex Align Stuff (GABLAM) ##
                    plen = {'Query':QryLen,'Hit':HitLen}                    
                    if self.opt['AlnStats']:
                        gdict = hit.dict['GABLAM']
                        #self.deBug(gdict)
                        for qh in ['Query','Hit']:
                            rqh = qh
                            if qh == 'Query': rqh = 'Qry'
                            if plen[qh] <= 0:
                                self.log.errorLog('Error in bamRead(%s, Hit %s): %s len < 1!' % (Qry,Hit,qh),printerror=False)
                                if qh == 'Hit': self.log.errorLog('Check for multiple occurrences of %s in search database!' % (Hit),printerror=False)
                                if self.opt['DeBug']:
                                    #X#os.system('cp %s.blast err.%s.blast' % (new_fork_id,new_fork_id))
                                    print hit.stat
                                    print gdict
                            for order in ['GABLAM','GABLAMO']:
                                for stat in ['Len', 'ID', 'Sim']:
                                    gstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order],stat)
                                    if self.opt['PercRes']:
                                        if plen[qh] > 0: resdict[hit][rstat] = '%.2f' % (100.0 * float(gdict[qh][gstat]) / plen[qh])
                                        else: resdict[hit][rstat] = '0.0'
                                    else: resdict[hit][rstat] = gdict[qh][gstat]
                                for stat in ['Dirn','Start','End']:
                                    gstat = '%s %s' % (order,stat)
                                    rstat = '%s_%s%s' % (rqh,ores[order][:-3],stat)
                                    resdict[hit][rstat] = gdict[qh][gstat]

                    ## Check GABLAM Cut ##
                    if self.stat['GABLAMCut'] > 0:
                        gcut = self.stat['GABLAMCut'] 
                        gstat = self.info['CutStat']
                        gstat = string.replace(gstat,'Ordered','O')
                        gstat = string.replace(gstat,'Aln',' ')
                        gstat = 'GABLAM' + gstat
                        if gstat not in gdict['Query'].keys():
                            print gstat, 'not cool. :o('
                            self.log.errorLog('GABLAM Cut-off stat "%s" not recognised. Should be generic (use with cutfocus for Query/Hit).' % gstat,False,False)
                            self.stat['GABLAMCut'] = 0
                            self.info['CutStat'] = 'AlnLen'                       
                        gtype = self.info['CutFocus'].lower()
                        if gtype not in ['query','qry','hit','either','both']:
                            self.log.errorLog('GABLAM Cut-off type "%s" not recognised. Will use "Either"' % self.info['CutType'],False,False)
                            self.info['CutFocus'] = 'Either'
                            gtype = 'either'
                        cutq = 100.0 * float(gdict['Query'][gstat]) / plen['Query']
                        cuth = 100.0 * float(gdict['Hit'][gstat]) / plen['Hit']
                        cut = False
                        if gtype in ['query','qry','both'] and cutq < gcut: cut = True
                        elif gtype in ['hit','both'] and cuth < gcut: cut = True
                        elif gtype in ['either'] and cuth < gcut and cutq < gcut: cut = True
                        self.verbose(1,2,'%.2f%% & %.2f%% vs %.2f%% %s (%s) => Cut=%s' % (cutq,cuth,gcut,gstat,gtype,cut),1)
                        if cut:
                            search.hit.remove(hit)
                            reslist = []
                            continue

                    ## Summary - now happens after GABLAM Cut-off ##
                    if not self_hit or self.opt['SelfSum']:
                        hit_num += 1
                        if Score > max_score :
                            max_score = Score
                            max_eval = EValue

                except:
                    self.log.errorLog('Major problem in searchGABLAM(%s, Hit=%s).' % (Qry,Hit))
                    search.hit.remove(hit)

            ### ALIGN ###
            resdict = self.searchALIGN(search,resdict,blast)
                    
            ### Output Results ###
            if self.opt['FullRes']:
                for hit in search.hit[0:]:
                    #X#print self.list['FullResHeaders'], resdict[hit]
                    rje.delimitedFileOutput(self,self.info['BamFile'],self.list['FullResHeaders'],delimit,resdict[hit])

            ### Local Stat Output ###
            if self.opt['Local']:
                for hit in search.hit[0:]:
                    for i in rje.sortKeys(hit.dict['Local']):
                        if hit.dict['Local'][i]['Length'] < self.stat['LocalMin']: continue
                        rje.delimitedFileOutput(self,self.info['LocFile'],self.list['LocalHeaders'],delimit,rje.combineDict({'Qry':Qry,'Hit':hit.info['Name'],'AlnNum':i},hit.dict['Local'][i]))

            ### Fasta Output ###
            if self.opt['FasOut'] and blast:
                bfas = os.path.join(rje.makePath(self.info['FasDir']),'%s.fas' % os.path.basename(Qry))
                outfas = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                if self.getBool('FragFas'):
                    blast.setOpt({'DeBug':self.getBool('DeBug')})
                    hitseq = blast.hitToSeq(seqlist=outfas,searchlist=[search],filename=None,appendfile=False)
                    #self.deBug(hitseq)
                    outfas.seq = []
                    for hit in search.hit[0:]:
                        #self.deBug(hit.cmd_list)
                        #self.bugPrint(hit.dict['GABLAM']['Hit'])
                        if hit not in hitseq: continue
                        ftot = len(hit.dict['GABLAM']['Hit']['GABLAM Frag'])
                        for fx in range(len(hit.dict['GABLAM']['Hit']['GABLAM Frag'])):
                            #self.bugPrint(hit.dict['GABLAM']['Hit']['GABLAM Frag'][fx])
                            (fragstart,fragend) = hit.dict['GABLAM']['Hit']['GABLAM Frag'][fx]
                            seq = hitseq[hit]
                            fullseq = seq.getSequence()
                            seqname = string.split(seq.info['Name'])
                            seqname[0] = '%s-%s.%s' % (seqname[0],rje.preZero(fragstart+1,seq.seqLen()),rje.preZero(fragend+1,seq.seqLen()))
                            seqname.insert(1,'GABLAM Fragment %d of %d (%s - %s)' % (fx+1,ftot,rje.iStr(fragstart+1),rje.iStr(fragend+1)))
                            seqname = string.join(seqname)
                            sequence = fullseq[fragstart:fragend+1]
                            outfas._addSeq(seqname,sequence)
                    if outfas.seq: outfas.saveFasta(seqfile=bfas,append=False); #self.deBug(outfas.seq[0].info)
                    elif search.hit: self.printLog('#ERR','0 sequences from %d hits output to %s!' % (len(search.hit),bfas))
                    #self.deBug(outfas.seq)
                else: blast.hitToSeq(seqlist=outfas,searchlist=[search],filename=bfas,appendfile=False)

            ### Summary ###
            if self.opt['HitSum']:
                sumres = {'Qry':Qry,'HitNum':hit_num,'MaxScore':max_score,'EVal':max_eval}
                self.debug('Legacy: %s' % sumres)
                rje.delimitedFileOutput(self,self.info['SumFile'],self.list['HitSumHeaders'],delimit,sumres)
                self.log.printLog('#HIT','%s => %d Database hits.' % (search.info['Name'],search.hitNum()))

        except:
            self.log.errorLog('GABLAM.searchGABLAMLegacy() is banjaxed',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def searchALIGN(self,search,resdict,blast):  ### Performs ALIGN for a single search (if desired) and updates resdict
        '''
        Performs ALIGN for a single search (if desired) and updates resdict.
        >> search:BLAST Search object with gablam already performed on hits
        >> resdict:dictionary of {Hit:{Header:string}}
        << resdict:dictionary with Align stats added
        '''
        try:### Setup ###
            if not self.opt['GlobID']:
                return resdict
            Qry = search.info['Name']
            QryLen = search.stat['Length']

            ### GlobID and QrySeq Object ###
            queryseq = None
            bseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
            blast.hitToSeq(seqlist=bseq,searchlist=[search],filename=None,appendfile=False)
            for h in range(len(search.hit)):
                hit = search.hit[h]
                Hit = hit.info['Name']
                if Qry.find(Hit)  >= 0 or (self.opt['QryAcc'] and Hit.find(Qry) >= 0):    # Self
                    queryseq = bseq.seq[h]
                    break

            ## Process ##
            for hit in search.hit[0:]:      # Already reduced by search GABLAM
                ## Basics ##
                Hit = hit.info['Name']
                try:
                    Score = hit.stat['BitScore']
                    Rank = string.atoi(resdict[hit]['Rank'])
                    EValue = hit.stat['E-Value']
                    HitLen = hit.stat['Length']
                    ## Global %ID from ALIGN ##
                    if (self.stat['RankAln'] > 0 and Rank < self.stat['RankAln']) or EValue > self.stat['EvalAln']:
                        continue
                    h = search.hit.index(hit)
                    hitseq = bseq.seq[h]
                    pwaln = bseq.pwAln(queryseq,hitseq)
                    if self.stat['AlnCut'] > 0:
                        qvh = float(pwaln.stat['Identity'] * pwaln.stat['Length']) / pwaln.stat['QryEnd']
                        hvq = float(pwaln.stat['Identity'] * pwaln.stat['Length']) / pwaln.stat['SbjEnd']
                        if self.stat['AlnCut'] > qvh or self.stat['AlnCut'] > hvq: continue
                    resdict[hit]['ALIGN_ID'] = '%d' % (pwaln.stat['Identity'] * pwaln.stat['Length'] / 100.0)
                    resdict[hit]['ALIGN_Len'] = '%d' % pwaln.stat['Length']
                except:
                    self.log.errorLog('Major problem in searchALIGN(%s, Hit=%s).' % (Qry,Hit))

            return resdict
        except:
            self.log.errorLog('GABLAM.searchALIGN() is banjaxed',printerror=True,quitchoice=True)
            return resdict
#########################################################################################################################
    def disMatrixOut(self): ### Output BLAST GABLAM distance matrix, based on rje_slimcore code.
        '''Output BLAST GABLAM distance matrix, based on rje_slimcore code.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['SeqList'].seqNum() > self.getInt('MaxAll') > 0:
                self.warnLog('Cancelling all-by-all output as %s sequences > maxall=%d' % (self.obj['SeqList'].seqNum(),self.getInt('MaxAll')))
                return
            gablam = self.obj['DisMat']
            gablam.basefile(self.basefile())
            diskey = self.getStr('DisKey')
            if diskey not in self.list['FullResHeaders']:
                if string.replace(diskey,'Ordered','') in self.list['FullResHeaders']:
                    diskey = string.replace(diskey,'Ordered','')
                    self.printLog('#DIS','Replaced distance matrix key (%s) with %s for output compatibility' % (self.getStr('DisKey'),diskey))
                elif string.replace(diskey,'_','_Ordered') in self.list['FullResHeaders']:
                    diskey = string.replace(diskey,'_','_Ordered')
                    self.printLog('#DIS','Replaced distance matrix key (%s) with %s for output compatibility' % (self.getStr('DisKey'),diskey))
                else: self.errorLog('#ERR','Cannot use "%s" as distance matrix key: not in output headers' % diskey); raise ValueError
            ### ~ [1] ~ Load GABLAM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#DIS','Reading Qry-Hit pairs from %s' % self.getStr('GablamOut'))
            gablam.loadFromDataTable(self.getStr('GablamOut'),key1='Qry',key2='Hit',distance=diskey,normalise=100.0,inverse=True,checksym=False)
            gablam.forceSymmetry(method='min',missing=1.0)
            ### ~ [2] ~ Output Matrix and associated all-by-all outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gablam.saveMatrix(filename=self.getStr('DisMatOut'),delimit=rje.delimitFromExt(filename=self.getStr('DisMatOut')))
            ## ~ [2a] ~ Generate tree outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gablam.rename(self.dict['NameMap'])    #!# Make tree-happy names #!#
            upgma = rje_tree.Tree(self.log,['savetype=none','treeformats=nwk,text,png']+self.cmd_list+['autoload=F'])
            if self.getBool('Clusters') or self.getBool('SaveUPC'):
                gclusters = gablam.cluster(maxdis=self.getNum('ClusterSplit'),singletons=True)
                if self.getBool('Clusters'): gablam.saveClusters(gclusters)
                if self.getBool('SaveUPC'): gablam.saveClusters(gclusters,upc=True)
            elif self.getInt('ByCluster'): gclusters = gablam.cluster(maxdis=1.0,singletons=True)
            else: gclusters = []
            self.printLog('#MAT','%s distance matrix objects. (%s sequences.)' % (rje.iStr(gablam.objNum()),rje.iStr(self.obj['SeqList'].seqNum())))
            if gablam.objNum() < 3: self.printLog('#TREE','Too few sequences for GABLAM tree output.')
            else: gablam.savePNG(upgma,savensf=self.getBool('DisTrees'),bycluster=self.getInt('ByCluster'),singletons=self.getBool('Singletons'),byclusters=gclusters)
            ## ~ [2b] ~ Generate graph outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DisGraph'):
                self.list['GraphTypes'] = string.split(string.join(self.list['GraphTypes']).lower())
                ppi = rje_ppi.PPI(self.log,self.cmd_list)
                basefile = ppi.info['Basefile'] = '%s.graph' % self.basefile()
                G = ppi.dict['PPI'] = gablam.makeGraph(cutoff=self.getNum('ClusterSplit'),singletons=False)    #self.getBool('Singletons')) Why not working?
                ppi.purgeOrphans()
                #!# Add output per byclusters sub-cluster! #!#
                npos = ppi.rjeSpringLayout(G)
                ppi.opt['ColByDeg'] = True
                if 'svg' in self.list['GraphTypes'] or 'html' in self.list['GraphTypes']:
                    svghtm = ppi.saveSVG(npos,basefile=basefile,G=G,font=0,width=1600,ntype='ellipse',backups=True)
                if 'png' in self.list['GraphTypes'] or 'html' in self.list['GraphTypes']:
                    ppi.saveR(npos,basefile=basefile,G=G,cleantdt=False,backups=True)
                if 'xgmml' in self.list['GraphTypes']:
                    xgmml = ppi.ppiToXGMML(G,basefile)
                    xgmml.saveXGMML('%s.xgmml' % basefile)
                if 'html' in self.list['GraphTypes']:
                    import rje_html        
                    html = rje_html.htmlHead(self.info['Name'],stylesheets=[],tabber=False) + svghtm
                    html += '\n<hr>\n<img src=%s.png width=1600>\n' % basefile
                    html += rje_html.htmlTail(copyright='RJ Edwards 2012',tabber=False)
                    open('%s.htm' % basefile,'w').write(html)
                self.printLog('#OUT','Graph output complete')
        except: self.errorLog('Major problem with %s.disMatrixOut' % self)
#########################################################################################################################
    def combinedFas(self):  ### Generates a combined sequence set for output
        '''Generates a combined sequence set for output.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            combinedfas = []    # List of (seqname,sequence) tuples
            bugger = '0'
            hitseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=list'])
            seqlist = rje_seqlist.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=T'])   #self.obj['SeqList']
            qryx = qseqx = 0   # Counter for sequences
            qtot = seqlist.seqNum()
            ### ~ [1] ~ Load BLAST mapped sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while seqlist.nextSeq():
                self.progLog('\r#COMB','Combining Hits: %.1f%%' % (100.0*qryx/qtot),rand=0.1); qryx +=1
                if self.opt['QryAcc']:
                    Qry = seqlist.seqAcc()
                    # Temp fix:
                    if rje.matchExp('^(sp|tr|sw|uniprot)\|(\S+)\|(\S+)_(\S+)',Qry):
                        Qry = rje.matchExp('^(sp|tr|sw|uniprot)\|(\S+)\|(\S+)_(\S+)',Qry)[1]
                else: Qry = seqlist.shortName()
                bfas = rje.makePath(self.getStr('FasDir')) + '%s.fas' % Qry
                self.debug('%s: %s' % (bfas,rje.exists(bfas)))
                if os.path.exists(bfas): hitseq.loadSeq(seqfile=bfas,nodup=True,clearseq=False,screen=self.v()>1); qseqx += 1
            self.printLog('\r#COMB','%s Combined Hits loaded for %s of %s Queries' % (rje.iLen(hitseq.seqs()),rje.iStr(qseqx),rje.iStr(qryx)))
            if not hitseq.seqNum(): return False    # Nothing to combine!
            if not hitseq.seqNum(): return False    # Nothing to combine!
            ### ~ [2] ~ Sort out GABLAM Fragments, if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('FragFas'): hitseq.saveSeq(seqfile='%s.fas' % self.basefile(),append=False); return
            #!# Debug with AAGE02013358.1 ?
            fragdict = {}; fx = 0
            for seq in hitseq.seqs(True):
                self.progLog('\r#FRAG','Compiling fragments: %s redundant fragments removed' % rje.iStr(fx))
                bugger = '1'
                try:
                    (name,fragstart,fragend) = rje.matchExp('^(\S+)-(\d+)\.(\d+)$',string.split(seq[0])[0])
                    maxlen = int('9' * len(fragend))
                    if name not in fragdict: fragdict[name] = []
                    addfrag = True
                    fragstart = int(fragstart); fragend = int(fragend)
                except:
                    self.errorLog('Shite')
                    self.printLog('#DEBUG',seq[0])
                    sys.exit()
                for (prevstart,prevend,prevseq) in fragdict[name]:
                    if fragstart >= prevstart and fragend <= prevend: addfrag = False
                if addfrag:
                    for (prevstart,prevend,prevseq) in fragdict[name]:
                        if fragstart <= prevstart and fragend >= prevend:
                            fragdict[name].remove((prevstart,prevend,prevseq)); fx += 1
                            bugger = '2'
                            prevqry = string.split(string.split(prevseq[0])[1],'|')
                            bugger = '3'
                            for qry in string.split(string.split(seq[0])[1],'|'):
                                if qry not in prevqry: prevqry.insert(1,qry)
                            prevqry.sort()
                            bugger = '4'
                            sname = string.split(seq[0])
                            sname[1] = string.join(prevqry,'|')
                            seq = (string.join(sname),seq[1])
                        # Combine overlaps! #!# Possibly combine within gabfrag distance too! #!#
                        elif prevstart <= fragstart <= prevend or prevstart <= fragend <= prevend: # Overlap
                            #i seq/prevseq are (name,seq) tuples
                            fragdict[name].remove((prevstart,prevend,prevseq)); fx += 1
                            # Extend fragstart/end and remake sequence and name #
                            startext = max(0,fragstart - prevstart)
                            endext = max(0,prevend - fragend)
                            if startext: seq = (seq[0],prevseq[1][:startext] + seq[1])
                            if endext: seq = (seq[0],seq[1] + prevseq[1][-endext:])
                            fragstart = min(fragstart,prevstart)
                            fragend = max(fragend,prevend)
                            if len(seq[1]) != fragend - fragstart + 1: raise ValueError()
                            bugger = '5'
                            prevqry = string.split(string.split(prevseq[0])[1],'|')
                            bugger = '6'
                            for qry in string.split(string.split(seq[0])[1],'|'):
                                if qry not in prevqry: prevqry.insert(1,qry)
                            prevqry.sort()
                            bugger = '7'
                            sname = string.split(seq[0])
                            sname[1] = string.join(prevqry,'|')
                            bugger = '8'
                            sname[-1] = string.split(sname[-1],'|')[0]
                            sname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                            seq = ('%s-%s.%s %s' % (name,rje.preZero(fragstart,maxlen),rje.preZero(fragend,maxlen),string.join(sname[1:])),seq[1])
                    fragdict[name].append((fragstart,fragend,seq))
                else: hitseq.list['Seq'].remove(seq); fx += 1
            self.printLog('\r#FRAG','Compiled fragments: %s redundant fragments removed.' % rje.iStr(fx))
            ## ~ [2a] ~ Rename and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newseq = []
            for name in rje.sortKeys(fragdict):
                fragdict[name].sort()
                for i in range(len(fragdict[name])):
                    seq = fragdict[name][i][-1]
                    #seqname.insert(1,'GABLAM Fragment %d of %d (%s - %s)' % (fx+1,ftot,rje.iStr(fragstart+1),rje.iStr(fragend+1)))
                    #seqname = string.split(seq[0])
                    #seqname[3] = '%d' % (i + 1)
                    #seqname[5] = '%d' % (len(fragdict[name]) + 1)
                    #newseq.append((string.join(seqname),seq[1]))
                    newseq.append(seq)
            hitseq.list['Seq'] = newseq
            hitseq.saveSeq(seqfile='%s.fas' % self.basefile(),append=False)
        except: self.errorLog('Major problem with %s.combinedFas(point %s)' % (self,bugger))
#########################################################################################################################
    def OLDcombinedFas(self):  ### Generates a combined sequence set for output
        '''Generates a combined sequence set for output.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            combinedfas = []    # List of (seqname,sequence) tuples
            bugger = '0'
            hitseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=list'])
            seqlist = rje_seqlist.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=T'])   #self.obj['SeqList']
            qseqx = 0   # Counter for sequences
            SEQFILE = open(self.info['QueryDB'], 'r')
            nextseq = lastline = 'Start'
            ### ~ [1] ~ Load BLAST mapped sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while nextseq:
                (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                if not nextseq: break
                seqlist.seq = [nextseq]
                qseqx +=1
                #self.progLog('\r#COMB','Combining Hits for %s Queries' % rje.iStr(qseqx))
                if self.opt['QryAcc']: Qry = nextseq.info['AccNum']
                else: Qry = nextseq.shortName()
                bfas = rje.makePath(self.info['FasDir']) + '%s.fas' % Qry
                if os.path.exists(bfas): hitseq.loadSeq(seqfile=bfas,nodup=True,clearseq=False)
            self.printLog('\r#COMB','%s Combined Hits loaded for %s Queries' % (rje.iLen(hitseq.seqs()),rje.iStr(qseqx)))
            if not hitseq.seqNum(): return False    # Nothing to combine!
            if not hitseq.seqNum(): return False    # Nothing to combine!
            ### ~ [2] ~ Sort out GABLAM Fragments, if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('FragFas'): hitseq.saveSeq(seqfile='%s.fas' % self.basefile(),append=False); return
            #!# Debug with AAGE02013358.1 ?
            fragdict = {}; fx = 0
            for seq in hitseq.seqs(True):
                self.progLog('\r#FRAG','Compiling fragments: %s redundant fragments removed' % rje.iStr(fx))
                bugger = '1'
                try:
                    (name,fragstart,fragend) = rje.matchExp('^(\S+)-(\d+)\.(\d+)$',string.split(seq[0])[0])
                    maxlen = int('9' * len(fragend))
                    if name not in fragdict: fragdict[name] = []
                    addfrag = True
                    fragstart = int(fragstart); fragend = int(fragend)
                except:
                    self.errorLog('Shite')
                    self.printLog('#DEBUG',seq[0])
                    sys.exit()
                for (prevstart,prevend,prevseq) in fragdict[name]:
                    if fragstart >= prevstart and fragend <= prevend: addfrag = False
                if addfrag:
                    for (prevstart,prevend,prevseq) in fragdict[name]:
                        if fragstart <= prevstart and fragend >= prevend:
                            fragdict[name].remove((prevstart,prevend,prevseq)); fx += 1
                            bugger = '2'
                            prevqry = string.split(string.split(prevseq[0])[1],'|')
                            bugger = '3'
                            for qry in string.split(string.split(seq[0])[1],'|'):
                                if qry not in prevqry: prevqry.insert(1,qry)
                            prevqry.sort()
                            bugger = '4'
                            sname = string.split(seq[0])
                            sname[1] = string.join(prevqry,'|')
                            seq = (string.join(sname),seq[1])
                        # Combine overlaps! #!# Possibly combine within gabfrag distance too! #!#
                        elif prevstart <= fragstart <= prevend or prevstart <= fragend <= prevend: # Overlap
                            #i seq/prevseq are (name,seq) tuples
                            fragdict[name].remove((prevstart,prevend,prevseq)); fx += 1
                            # Extend fragstart/end and remake sequence and name #
                            startext = max(0,fragstart - prevstart)
                            endext = max(0,prevend - fragend)
                            if startext: seq = (seq[0],prevseq[1][:startext] + seq[1])
                            if endext: seq = (seq[0],seq[1] + prevseq[1][-endext:])
                            fragstart = min(fragstart,prevstart)
                            fragend = max(fragend,prevend)
                            if len(seq[1]) != fragend - fragstart + 1: raise ValueError()
                            bugger = '5'
                            prevqry = string.split(string.split(prevseq[0])[1],'|')
                            bugger = '6'
                            for qry in string.split(string.split(seq[0])[1],'|'):
                                if qry not in prevqry: prevqry.insert(1,qry)
                            prevqry.sort()
                            bugger = '7'
                            sname = string.split(seq[0])
                            sname[1] = string.join(prevqry,'|')
                            bugger = '8'
                            sname[-1] = string.split(sname[-1],'|')[0]
                            sname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                            seq = ('%s-%s.%s %s' % (name,rje.preZero(fragstart,maxlen),rje.preZero(fragend,maxlen),string.join(sname[1:])),seq[1])
                    fragdict[name].append((fragstart,fragend,seq))
                else: hitseq.list['Seq'].remove(seq); fx += 1
            self.printLog('\r#FRAG','Compiled fragments: %s redundant fragments removed.' % rje.iStr(fx))
            ## ~ [2a] ~ Rename and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newseq = []
            for name in rje.sortKeys(fragdict):
                fragdict[name].sort()
                for i in range(len(fragdict[name])):
                    seq = fragdict[name][i][-1]
                    #seqname.insert(1,'GABLAM Fragment %d of %d (%s - %s)' % (fx+1,ftot,rje.iStr(fragstart+1),rje.iStr(fragend+1)))
                    #seqname = string.split(seq[0])
                    #seqname[3] = '%d' % (i + 1)
                    #seqname[5] = '%d' % (len(fragdict[name]) + 1)
                    #newseq.append((string.join(seqname),seq[1]))
                    newseq.append(seq)
            hitseq.list['Seq'] = newseq
            hitseq.saveSeq(seqfile='%s.fas' % self.basefile(),append=False)
        except: self.errorLog('Major problem with %s.combinedFas(point %s)' % (self,bugger))
#########################################################################################################################
    def pngDotPlots(self):
        '''Generate dotplot outputs using R.'''
        try:### ~ [1] ~ Try to run R to generate PNG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RPATH',self.info['RPath'],log=False,screen=self.v()>0)
            self.info['PathR'] = rje.makePath(os.path.abspath(string.join(string.split(sys.argv[0],os.sep)[:-2]+['libraries','r',''],os.sep)))
            rscript = rje.makePath('%sgablam.r' % (self.info['PathR']),wholepath=True)
            rcmd = '%s --no-restore --no-save --args "%s" "%s"' % (self.info['RPath'],self.baseFile(),self.getInt('DotLocalMin'))
            rcmd += ' < "%s" > "%s.dot.r.run"' % (rscript,self.baseFile())
            dotpath = '%s.DotPlots/' % self.baseFile()
            #self.deBug('%s:%s' % (dotpath, os.path.exists(dotpath)))
            rje.mkDir(self,dotpath,True)
            self.printLog('#RPLOT',rcmd)
            problems = os.popen(rcmd).read()
            if not self.getBool('DeBug') and not self.getBool('Test') and os.path.exists('%s.dot.r.run' % self.baseFile()): os.unlink('%s.dot.r.run' % self.baseFile())
            if problems: self.errorLog(problems,printerror=False)
            return rcmd
        except: self.log.errorLog('Major Problem with Tree.pngTree()')
#########################################################################################################################
### End of SECTION II: GABLAM Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: BAM Class                                                                                              #
#########################################################################################################################
class BAM(GABLAM):     
    '''For backwards compatibility only'''
    def bam(self): return self.gablam()
#########################################################################################################################
### END OF SECTION III: BAM Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: GABLAM(log=mainlog,cmd_list=cmd_list).run()
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
