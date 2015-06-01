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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      PINGU
Description:  Protein Interaction Network & GO Utility
Version:      3.9
Last Edit:    16/07/13
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This utility was originally created for handling proteomics data with EnsEMBL peptide IDs. The data needed to be
    mapped onto Genes, overlaps and redundancies identified, gene lists output for GO analysis with FatiGO, and PPI data
    from HPRD and BioGRID to identify potential complexes.

    See rje_ensembl documentation for details of what to download for EnsGO files and how to make EnsLoci files etc. The
    ens_SPECIES.GO.tdt file used for go mapping should be suitable for the ensmap=FILE file.

    Giving a ppioutdir=PATH will produce combined PPI sequence files for all genes. A pingu.combinedppi.tdt summary file
    will be placed in resdir.

    QSLiMFinder=FILE will perform an analysis for shared motifs between primary interactors of those genes identified
    in a given sample and the original sequence used for the pulldown. FILE should be a fasta file where names of the
    sequences match Sample names. Datasets will then be formed that contain that sequence plus the primary PPI of each
    gene in that sample as a dataset named SAMPLE_GENE.fas in a directory RESDIR/SLiMFinder.

Commandline:
    ### Main Input Options ###
    data=LIST       : List of files of results containing "Sample" and "Identifier" columns
    ensmap=FILE     : Mappings from EnsEMBL - peptides, genes and HGNC IDs (from BioMart)
    ipilinks=FILE   : IPI Links file with 'IPI', 'Symbol' and 'EnsG' fields []
    ensloci=FILE    : File of EnsEMBL genome EnsLoci treatment []
    baits=LIST      : List of genes of interest for overlap analysis []
    addbaits=T/F    : Whether to add primary interactors of baits as additional samples [False]
    combaits=X      : Whether to combine bait PPIs into single sample (X) (if addbaits=T) []
    controls=LIST   : List of sample names that correspond to controls [Control]
    experiments=LIST: List of sample names that correspond to key samples of interest []
    exponly=T/F     : Limited analysis to samples listed as experiments (before baits added etc.) [False]
    addalias=FILE   : Extra (manual?) aliases to add to GeneMap object following loading of pickles etc. [None]

    ### Processing Options ###
    hgnconly=T/F    : Whether to restrict PPI data to only those proteins with Gene Symbol links [False]
    pickle=T/F      : Whether to save/load pickle of parsed/combined data rather than regenerating each time [True]
    pingupickle=FILE: Full path to Pingu pickle file to look for/use/save [pingu.pickle]
    nocontrols=T/F  : Whether to remove genes found in designated controls from designated experiments [False]
    gablam=T/F      : Whether to run all-by-all GABLAM on EnsLoci and add homology to networks [True]
    ppitype=LIST    : List of acceptable interaction types to parse out []
    badtype=LIST    : List of bad interaction types, to exclude [indirect_complex,neighbouring_reaction]
    makefam=X       : GABLAM Percentage identity threshold for grouping sequences into families [0.0]
    gofilter=LIST   : List of GO IDs to filter out of gene lists []
    goexcept=LIST   : List of GO ID exceptions to filtering []
    remsticky=X     : Remove "sticky" hubs as defined by >X known PPI [0]
    stickyhubs=T/F  : Only remove "sticky" spokes but keep sticky hubs [False]
    stickyppi=T/F   : Only remove "sticky" hubs from samples, not from total PPI [False]
    addlinks=T/F    : Add linking proteins (linking two Sample proteins) [False]
    
    ### Main Output Options ###
    resdir=PATH     : Redirect output files to specified directory [./]
    basefile=X      : Results file prefix if no data file given with data=FILE [pingu]
    fulloutput=T/F  : Generate all possible outputs from one input [False]
    genelists=T/F   : Generate lists of genes for each sample (e.g. for FatiGO upload) [False]
    gosummary=T/F   : Make a GO summary table [False]
    summaryhgnc=T/F : Generate a summary table of genes in dataset, including peptide lists for each sample [False]
    mapout=T/F      : Generate a summary table of full peptide mapping [False]
    dbcomp=T/F      : Comparison of PPI databases [False]
    dbsizes=T/F     : Outputs a file of PPI dataset sizes (histogram) [False]
    allbyall=X      : Generates an all-by-all table of PPI links upto X degrees of separation (sample only) [0]
    pathfinder=X    : Perform (lengthy) PathFinder analysis to link genes upto X degree separation (-1 = no limit) [0]
    pathqry=LIST    : Limit PathFinder analysis to start with given queries []
    overlap=T/F     : Produce a table of the overlap (mapped through HGNC) between samples (and bait 1y PPI) [False]
    cytoscape=T/F   : Produce old cytoscape input files from allbyall table (reads back in) [False]
    xgmml=T/F       : Produce an XGMML file with all Cytoscape data and more [False]
    xgformat=T/F    : Whether to add colour/shape formatting to XGMML output [False]
    xgexpand=X      : Expand XGMML network with additional levels of interactors [0]
    xgcomplex=T/F   : Restrict XGMML output (and expansion) to protein complex edges [False]
    compresspp=T/F  : Whether to compress multiple samples of interest into ShareX for cytoscape [False]
    seqfiles=T/F    : Whether to generate protein sequence fasta files using EnsLoci [False]
    goseqdir=PATH   : Path to output full GO fasta files (No output if blank/none) []
    ppioutdir=PATH  : Path to output combined PPI files (No output if blank/none) []
    acconly=T/F     : Whether to output lists of Accession numbers only, rather than full fasta files [False]
    ensdat=PATH     : Path to EnsDAT files to use for making combined PPI datasets [None]
    qslimfinder=FILE: File containing sequences matching Sample names for Query SLiMFinder runs [None]
    screenddi=FILE  : Whether to screen out probably domain-domain interactions from file [None]
    domppidir=PATH  : Produce domain-based PPI files and output into PATH (No output if blank/none) []
    nocomplex=T/F   : Perform crude screening of complexes (PPI triplets w/o homodimers) [False]
    fasid=X         : Text ID for PPI fasta files ['ppi']
    association=T/F : Perform experiment association analysis [False]
    asscombo=T/F    : Whether to subdivide genes further based on combinations of experiments containing them [False]
    noshare=T/F     : Whether to exclude those genes that are shared between samples when comparing those samples [True]
    selfonly=T/F    : Whether to only look at associations within experiments, not between [False]
    randseed=X      : Seed for randomiser [0]
    randnum=X       : Number of randomisations [1000]
    
    ### Database/Path options ###
    enspath=PATH    : Path to EnsEMBL downloads
    ensgopath=PATH  : Path to EnsGO files   (!!! Restricted to Humans Currently !!!)
    unipath=PATH    : Path to UniProt files [UniProt/]
    hprd=PATH       : Path to HPRD flat files [None]
    biogrid=FILE    : BioGRID flat file [None]
    intact=FILE     : IntAct flat file [None]
    mint=FILE       : MINT flat file [None]
    reactome=FILE   : Reactome interactions flat file [None]
    dip=FILE        : DIP interactions flat file [None]
    domino=FILE     : Domino interactions flat file [None]
    pairwise=FILE   : Load interaction data from existing Pingu Pairwise file [None]
    addppi=FILE     : Add additional PPI from a simple delimited file IDA,IDB,Evidence [None]
    genepickle=FILE : Pickled GeneMap object. Alternatively, use below commands to make GeneMap object [None]
    - hgncdata/sourcedata/pickledata/aliases : See rje_genemap docstring.
    pfamdata=FILE   : Delimited files containing domain organisation of sequences [None]
    evidence=FILE   : Mapping file for evidence terms [None]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_biogrid, rje_genemap, rje_go, rje_hprd, rje_ppi, rje_seq, rje_uniprot, rje_xgmml, rje_zen
import rje_dismatrix_V2 as rje_dismatrix
import gablam, qslimfinder, rje_slimcore, rje_tree
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.1 - Initial Compilation. Basic GO mapping for EnsEMBL data.
    # 0.2 - Mapping of EnsEMBL genes onto Gene Symbols and summary data table output.
    # 0.3 - Reading and collation of PPI data.
    # 0.4 - All-by-all PPI and sample overlap analyses.
    # 0.5 - Cytoscape output from all-by-all analysis.
    # 0.6 - Option to add interactors of baits as additional samples. Added resdir=PATH option. Added pickling.
    # 0.7 - Added generation of EnsLoci datasets
    # 0.8 - Combined PPI Dataset output
    # 0.9 - Added DAT output option.
    # 1.0 - Full working version with XGMML output including GABLAM relationships.
    # 1.1 - Added Reactome and DIP.
    # 1.2 - Added an output of shared PPIs for clustering. (*.cluster.tdt)
    # 2.0 - Replace GeneCards with GeneMap. Improved compatibility with APHID and functioning of new options.
    # 2.1 - Added rje_go.GO Object to store GO mappings etc.
    # 2.2 - Altered the Pingu data input to be a list of files, not just one file.
    # 2.3 - Added acclist and GO dataset outputs.
    # 2.4 - Added tracking of source databases and evidence codes.
    # 2.5 - Added PNG visualisations.
    # 2.6 - Added loading interactions from pingu.pairwise.tdt.
    # 2.7 - Expanded XGMML options and added pathfinder output.
    # 3.0 - Major tidying of gene/peptide mapping. Added extra bait, experiment and protein family options.
    # 3.1 - Removal of sticky spokes/hubs
    # 3.2 - Added Domain-based PPI output.
    # 3.3 - Added crude Complex filtering.
    # 3.4 - Updated GO stuff.
    # 3.5 - Added Experiment association output.
    # 3.6 - Added addlinks=T/F option.
    # 3.7 - Improved XGMML output.
    # 3.8 - Hopefully fixed issue of Fasta file generation log output writing to wrong log file.
    # 3.9 - Tidied imports.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Replace use of self.dict['GeneMap'] with self.obj['GeneMap']
    # [ ] : Add different species - C. elegans, Drosophila, Yeast.
    # [ ] : Alter overlap output to have PPI in rows only and to have Baits in columns to (for bait vs bait).
    # [ ] : Add GO / GO_SLiM to network output.
    # [ ] : Add GO-HGNC mapping from somewhere.
    # [ ] : Add SLiMFinder runs (without query)
    # [ ] : Improve/add a specific run for pulling out PPI for baits alone and running SLiMFinder etc.
    # [ ] : Identify source of plain number identifiers. (HPRD?)
    # [ ] : Add an output of shared PPIs for clustering. (*.shared.tdt)
    # [ ] : Add a trace function for tracing certain genes through whole program for debugging?!
    # [ ] : Add filtering of GO categories to PPI genes too
    # [ ] : Add additional directory for Pingu Construction output (*.missing.txt etc.)?
    # [ ] : Modify mapout and hgncmapping to cope with family grouping.
    # [ ] : Adding reading of pseudogenes from EnsEMBL map file and identify rather than report no mapping.
    # [ ] : Found out where "CHEBI and P16638? come from and fix.
    # [ ] : Modify Domino reading to use separate binary and complex files.
    # [ ] : Needs a massive overhaul of functions and import classes that now do functions.
    # [ ] : Move some of the bigger functions to slimsuite and/or seqsuite modules.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('PINGU', '3.9', 'July 2013', '2007')
    description = 'Protein Interaction Network & GO Utility'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program needs an R installation for PNG visualisations.',
                'WARNING: This program remains in test phase only',rje_zen.Zen(None,[]).wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show GeneMap commandline options?'): out.verbose(-1,4,text=rje_genemap.__doc__)
            if rje.yesNo('Show BLAST commandline options (for GABLAM)?'): out.verbose(-1,4,text=gablam.__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PINGU Class                                                                                             #
#########################################################################################################################
class PINGU(rje.RJE_Object):     
    '''
    Main Protein Interaction Network & GO Utility Class. Author: Rich Edwards (2007).

    Info:str
    - AddAlias = Extra (manual?) aliases to add to GeneMap object following loading of pickles etc. [None]
    - AddPPI = Add additional PPI from a simple delimited file IDA,IDB,Evidence [None]
    - EnsGOPath = Path to EnsGO files 
    - EnsMap = File of mappings from EnsEMBL - peptides, genes and HGNC IDs
    - EnsLoci = File of EnsEMBL genome EnsLoci treatment []
    - Evidence = Mapping file for evidence terms [None]
    - Name = File of proteomics results containing "Sample" and "Identifier" columns
    - BioGrid = BioGRID flat file [None]
    - Domino = Domino interactions flat file [None]
    - DIP = DIP interactions flat file [None]
    - HPRD = Path to HPRD flat files [None]
    - IntAct = IntAct flat file [None]
    - MINT = MINT flat file [None]
    - Reactome = Reactome interactions flat file [None]
    - IPILinks = IPI Links file with 'IPI', 'Symbol' and 'EnsG' fields []
    - Pairwise = Load interaction data from existing Pingu Pairwise file [None]
    - PinguPickle = Full path to Pingu pickle file to look for/use/save [pingu.pickle]
    - PPIOutDir = Path to output combined PPI files (No output if blank/none) []
    - GOSeqDir = Path to output full GO fasta files (No output if blank/none) []    
    - ResDir = Redirect output files to specified directory [./]
    - EnsDat = Path to EnsDAT files to use for making combined PPI datasets [None]
    - QSLiMFinder = File containing sequences matching Sample names for Query SLiMFinder runs [None]
    - ScreenDDI = Whether to screen out probably domain-domain interactions from file [None]
    - ComBaits = Whether to combine bait PPIs into single sample (X) (if addbaits=T) []
    - DomPPIDir = Produce domain-based PPI files and output into PATH (No output if blank/none) []
    - FasID = Text ID for PPI fasta files ['ppi']
        
    Opt:boolean
    - AccOnly = Whether to output lists of Accession numbers only, rather than full fasta files [True]
    - AddLinks = Add linking proteins (linking two Sample proteins) [False]
    - Association = Perform experiment association analysis [False]
    - AssCombo = Whether to subdivide genes further based on combinations of experiments containing them [False]
    - SelfOnly = Whether to only look at associations within experiments, not between [False]
    - FullOutput = Generate all possible outputs from one input [False]
    - GeneLists = Generate lists of genes for each sample (e.g. for FatiGO upload) [False]
    - GOSummary = Make a GO summary table [False]
    - SummaryHGNC = Generate a summary table of genes in dataset, including peptide lists for each sample [False]
    - MapOut = Generate a summary table of full peptide mapping [False]
    - DBComp = Comparison of PPI databases [False]
    - HGNCOnly = Whether to restrict PPI data to only those proteins with Gene Symbol links [False]
    - DBSizes = Outputs a file of PPI dataset sizes (histogram) [False]
    - Overlap = Produce a table of the overlap (mapped through HGNC) between samples (and 1y PPI) [False]
    - Cytoscape = Produce cytoscape input files from allbyall table (reads back in) [False]
    - AddBaits = Whether to add primary interactors of baits as additional samples [False]
    - CompressPP = Whether to compress multiple samples of interest into pXp for cytoscape [False]
    - Pickle = Whether to save/load pickle of parsed/combined data rather than regenerating each time [True]
    - SeqFiles = Whether to generate protein sequence fasta files using EnsLoci [False]
    - NoControls = Whether to remove genes found in designated controls from designated experiments [False]
    - NoShare = Whether to exclude those genes that are shared between samples when comparing those samples [True]
    - DatFiles = Whether to generate protein DAT files using EnsEMBL gene AccNum and UniPath [False]
    - XGMML = Output XGMML-format network file [False]
    - XGFormat = Whether to add colour/shape formatting to XGMML output [False]
    - XGComplex = Restrict XGMML output (and expansion) to protein complex edges [False]
    - GABLAM = Generate GABLAM relationships and add to networks [False]
    - GenericData = Whether genericData() method has been run [False]
    - ExpOnly = Limited analysis to samples listed as experiments (before baits added etc.) [False]
    - StickyHubs = Only remove "sticky" spokes but keep sticky hubs [False]
    - StickyPPI = Only remove "sticky" hubs from samples, not from total PPI [False]
    - NoComplex = Perform crude screening of complexes (PPI triplets w/o homodimers) [False]

    Stat:numeric
    - AllByAll = Generates an all-by-all table of PPI links (sample only) [0]
    - MakeFam = GABLAM Percentage identity threshold for grouping sequences into families [0.0]
    - RandSeed = Seed for randomiser [0]
    - RandNum = Number of randomisations [1000]
    - RemSticky = Remove "sticky" hubs as defined by >X known PPI [0]
    - XGExpand = Expand XGMML network with additional levels of interactors [0]

    List:list
    - Baits = List of genes of interest for overlap analysis []
    - Controls = List of sample names that correspond to controls [Control]
    - Data = List of input files []
    - EnsGenes = Full list of EnsEMBL genes with GO data
    - Experiments = List of sample names that correspond to key samples of interest []
    - Genes = Full list of final NR genes
    - AmbigAlias = List of Ambiguous aliases from HGNC - use to assess mapping problems []
    - PathFinder = Perform (lengthy) PathFinder analysis to link genes upto X degree separation (-1 = no limit) [0]
    - PPIType = List of acceptable interaction types to parse out []
    - BadType = List of bad interaction types, to exclude [indirect_complex,neighbouring_reaction]
    - GOFilter = List of GO IDs to filter out of gene lists []
    - GOExcept = List of GO ID exceptions to filtering []

    Dict:dictionary
    - Samples = Dictionary of {sample:[peptide IDs]}
    - Datasets = Dictionary of {sample:[NR mapped genes]}
    - Evidence = Evidences for interactions {ID:{ID:[(DBase,Evidence types)]}}
    - BackMap = Dictionary of {NR gene:[prot IDs]} backwards mappings to original data
    - GO = Dictionary of Gene-GO Mappings {Gene:{[GO]}}
    - PathQry = Limit PathFinder analysis to start with given queries []
    - PPIDB = Dictionary of PPIDB objects {Type:Obj}, where type matches self.info and Obj is rje_biogrid or rje_hprd
    - PPI = Dictionary of compiled PPI data {Gene:[Interactors]}
    - EnsSeqNameDict = Dictionary of short names mapped to sequences for EnsLoci
    - DDI = Dictionary of {Domain:[Interacting domains]}
    - Families = Dictionary of {FamName:[Genes]}
    - FamMap = Dictionary of {Gene:FamName}
    - FamPPI = Dictionary of {Fam:[FamPPIs]}

    Obj:RJE_Objects
    - DisMatrix = rje_dismatrix_V2.DisMatrix object
    - EnsEMBL = rje_ensembl.EnsEMBL object
    - EnsLoci = SeqList object containing EnsLoci sequences
    - GABLAM = GABLAM object used for homology
    - GeneMap = rje_genemap.GeneMap object
    - GO = rje_go.GO Object
    - Homology = rje_dismatrix_V2.DisMatrix containing Homology relationships
    '''
    def go(self,id=None): return self.obj['GO'].go(id)
    def mapper(self): return self.obj['GeneMap']
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['EnsGOPath','EnsMap','BioGrid','HPRD','IntAct','MINT','IPILinks','ResDir','EnsLoci','Domino',
                         'PPIOutDir','EnsDat','Reactome','DIP','QSLiMFinder','AddAlias','GOSeqDir','ScreenDDI',
                         'Evidence','Pairwise','AddPPI','ComBaits','PinguPickle','DomPPIDir','PFamData','FasID']
        self.optlist = ['FullOutput','GOSummary','SummaryHGNC','MapOut','DBComp','HGNCOnly','DBSizes','Overlap','XGMML',
                        'Cytoscape','AddBaits','CompressPP','Pickle','SeqFiles','NoControls','DatFiles','GABLAM',
                        'AccOnly','GenericData','XGFormat','GeneLists','ExpOnly','StickyHubs','NoComplex','Association',
                        'AssCombo','SelfOnly','NoShare','AddLinks','StickyPPI','XGComplex']
        self.statlist = ['AllByAll','XGExpand','PathFinder','MakeFam','RemSticky','RandSeed','RandNum']
        self.listlist = ['Baits','Controls','Experiments','Genes','AmbigAlias','Data','PPIType','BadType','PathQry',
                         'GOFilter','GOExcept','EnsGenes']
        self.dictlist = ['Samples','Datasets','BackMap','GO','PPIDB','DDI','PPI',
                         'Evidence','Families','FamMap','FamPPI']
        self.objlist = ['EnsEMBL','GABLAM','GeneMap','EnsLoci','DisMatrix','Homology']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
        self.setInfo({'ResDir':rje.makePath('./'),'Basefile':'pingu','PinguPickle':'pingu.pickle','FasID':'ppi'})
        self.setStat({'AllByAll':0,'XGExpand':0,'PathFinder':0,'RemSticky':0,'RandSeed':0,'RandNum':1000})
        self.setOpt({'Pickle':True,'GABLAM':False,'AccOnly':False,'NoShare':True,'StickyPPI':False})
        self.list['Controls'] = ['Control']
        self.list['BadType'] = ['indirect_complex','neighbouring_reaction']
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                #x#self._cmdRead(cmd,type='info',att='Name',arg='data')
                self._cmdReadList(cmd,'int',['AllByAll','XGExpand','PathFinder','RemSticky','RandSeed','RandNum'])
                self._cmdReadList(cmd,'stat',['MakeFam'])
                self._cmdReadList(cmd,'info',['ComBaits','FasID'])
                self._cmdReadList(cmd,'file',['EnsMap','BioGrid','HPRD','IntAct','MINT','IPILinks','EnsLoci','Domino',
                                              'Reactome','DIP','QSLiMFinder','AddAlias','ScreenDDI','Evidence',
                                              'Pairwise','AddPPI','PinguPickle','PFamData'])
                self._cmdReadList(cmd,'path',['EnsGOPath','HPRD','ResDir','PPIOutDir','EnsDat','GOSeqDir','DomPPIDir'])
                self._cmdReadList(cmd,'opt',['FullOutput','GOSummary','SummaryHGNC','MapOut','DBComp','HGNCOnly','XGMML',
                                             'DBSizes','Overlap','Cytoscape','AddBaits','CompressPP','Pickle','SeqFiles',
                                             'NoControls','DatFiles','GABLAM','AccOnly','XGFormat','GeneLists','ExpOnly',
                                             'StickyHubs','NoComplex','Association','AssCombo','SelfOnly','NoShare',
                                             'AddLinks','StickyPPI','XGComplex'])
                self._cmdReadList(cmd,'list',['Baits','Controls','Experiments','PPIType','BadType','PathQry','GOFilter',
                                              'GOExcept'])
                self._cmdReadList(cmd,'glist',['Data'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.list['PPIType']: self.list['PPIType'] = string.split(string.join(self.list['PPIType'],'|').lower(),'|')
        if self.list['BadType']: self.list['BadType'] = string.split(string.join(self.list['BadType'],'|').lower(),'|')
#########################################################################################################################
    ### <2> ### Main Class Run Methods                                                                                  #
#########################################################################################################################
    def run(self):  ### Main run method. Loads generic data and then processes specific samples.
        '''Main run method. Loads generic data and then processes specific samples.'''
        try:### ~ Run Full Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.genericData()
            self.process()
        except: self.errorLog('Pingu.run() Error')
#########################################################################################################################
    def processPickle(self,newme):  ### Changes attributes according to new Pickle
        '''Changes attributes accordingly.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict = newme.dict
            self.obj = newme.obj
            newme.setLog(self.log,cascade=True)
            return newme
        except: self.errorLog('Problem during %s.processPickle()' % self); return None
#########################################################################################################################
    def genericData(self):  ### Loads generic data: ID mappings, GO terms and PPI data.
        '''Loads generic data: ID mappings, GO terms and PPI data..'''
        try:
            ### ~ [1] Look for Pickle of Mappings and PPI Data etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'GenericData' in self.opt and self.opt['GenericData']: return self.printLog('#SKIP','Skipping Pingu.genericData() - already run.')
            mypickle = None
            try: pfile = self.info['PinguPickle']
            except: pfile = 'pingu'
            if pfile[-3:] == '.gz': pfile = pfile[:-3]
            if pfile[-7:] == '.pickle': pfile = pfile[:-7]
            if self.opt['Pickle']: mypickle = self.unpickleMe(pfile)
            ### ~ [2] ~ Load fresh data if Pickle not found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not mypickle:
                ## ~ [2a] Look for GeneMap data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.obj['GeneMap'] = rje_genemap.GeneMap(self.log,self.cmd_list)
                self.obj['GeneMap'].run()
                ## ~ [2b] Parse basic data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                       
                self.mapEnsEMBL()   # Parses mappings from EnsEMBL peptides to genes to HGNC
                self.mapIPI()       # Parses mappings from IPI links table to Ens/HGNC
                self.parseGO()      # Parses GO terms, types and descriptions etc. using rje_go
                self.mapEnsGO()     # Maps genes onto GO terms using BioMart download
                ## ~ [2c] Read in data from PPI Databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if rje.checkForFile(self.info['Pairwise']): self.parsePairwisePPI()
                else:
                    self.parsePPI()     # Makes self.dict['PPIDB'] of PPI objects and parses using rje_hprd & rje_biogrid
                    self.convertPPI()   # Converts PPI objects into HGNC symbol-based interactions rather than prot IDs
                    self.evidence()     # If appropriate, will convert and output a summary of the interaction types read 
                    self.combinePPI()   # Combines converted PPI objects into self.dict['PPI']
            ### ~ [3] ~ Add additional aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.addAlias()
            self.addPPI()
            ### ~ [4] ~ Save Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Pickle'] and not mypickle: self.pickleMe(pfile)
            self.opt['GenericData'] = True
        except: self.errorLog('Pingu.genericData() Error')
#########################################################################################################################
    def process(self):  ### Processes sample data after reading in generic PPI etc. data
        '''Processes sample data after reading in generic PPI etc. data.'''
        try:
            ### ~ [1] Process Gene data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['ResDir'])
            for dfile in self.list['Data']: self.readData(dfile)    # Reads datasets of samples & EnsEMBL peptides / IPI / HGNC Genes
            if len(self.list['Data']) == 1: self.info['Name'] = self.list['Data'][0]
            else: self.info['Name'] = self.info['Basefile']
            self.info['Name'] = rje.baseFile(self.info['Name'],strip_path=True)
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = self.info['Name']
            self.info['Basefile'] = self.info['ResDir'] + self.info['Basefile']
            if self.opt['AddBaits']: self.addBaits()

            ### ~ [2] Convert EnsEMBL to NR/HGNC gene symbols ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.makeNRGenes()      # Converts peptide datasets into NR gene datasets
            if self.stat['RemSticky'] > 0: self.remSticky(purgeppi=not self.opt['StickyPPI'])     # Remove "sticky" hub proteins/families
            if self.opt['AddLinks']: self.addLinks()
            #x# self.convertToHGNC()    # Convert NR gene datasets to symbols and update links etc. #!# REDUNDANT? #!#
            self.filterGOSamples()  # Filter genes using GO data
            self.filterGOPPI()      # Filter genes using GO data
            fullout = self.opt['FullOutput']
            ## ~ [2a] Sequence file output and possible family clustering ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.stat['MakeFam']: self.opt['GABLAM'] = True
            if fullout or self.opt['SeqFiles'] or self.opt['GABLAM']: self.sampleSeqFiles()
            if self.opt['GABLAM']: self.gablamMatrix()
            if self.stat['MakeFam'] > 0: self.makeFam()         # Generate families based on GABLAM data

            ### ~ [3] Basic PINGU Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fullout or self.opt['SummaryHGNC']: self.summaryHGNC()    # Summary table of genes in dataset, including peptide lists for each sample
            if fullout or self.opt['MapOut']: self.mapOut()              # Summary table of full peptide mapping
            if fullout or self.opt['GOSummary']: self.sampleGO()         # Perform sample GO Analysis
            self.goSeqFiles()

            ### ~ [4] PPI Dataset outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['NoControls']: self.purgeControls()
            if self.opt['DBSizes'] or fullout: self.sizeDB()    # Dataset size output (for histogram)
            ## ~ [4a] Generic PPI Output - sample independent ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['DBComp'] or fullout: self.compDB()     # Generates database comparison output
            self.screenPPI()                                    # Screen based on evidence and predicted complexes
            self.screenDDI()
            self.combinedPPISeqFiles()
            self.domPPI()
                
            ### ~ [5] Generate All-by-All outputs and analyses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fullout or self.opt['Overlap']: self.sampleOverlap()
            if self.stat['AllByAll']: self.allByAll()
            if fullout or self.opt['Association']: self.association()
            if self.stat['PathFinder']: self.pathFinder()
            else: self.printLog('#PATH','No PathFinder analysis')
            if self.opt['Cytoscape']: self.sampleCytoscape()
            if fullout or self.opt['XGMML']:
                xgmml = self.sampleXGMML()
                ppiobj = rje_ppi.PPI(self.log,self.cmd_list)
                ppiobj.xgmmlToPPI(xgmml)
                ppiobj.addCol()
                ppiobj.main()

            ### ~ [6] Special SLiMFinder analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.qSLiMFinder()
        except: self.log.errorLog('Pingu.process() Error')
#########################################################################################################################
    def debugGeneCheck(self):
        gene = 'C3orf62'.upper()
        if not self.obj['GeneMap'].getGeneData(gene): self.deBug('%s >> Fuck' % gene)
        else: self.deBug('Yo!: %s' % self.obj['GeneMap'].getGeneData(gene))
#########################################################################################################################
    ### <3> ### Data reading and parsing methods                                                                        #
#########################################################################################################################
    def readData(self,filename):     ### Reads experimental data into dictionary
        '''Reads experimental data into dictionary.'''
        try:### ~ [1] ~ Check whether sample file has been given ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename.lower() in ['','none']: return self.log.printLog('#DATA','No samples to read/process')
            elif not os.path.exists(filename): return self.log.printLog('#DATA','Input file "%s" missing!' % filename)
            base = rje.baseFile(filename,strip_path=True)
            ### ~ [2] ~ Read Data from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ix = 0
            indata = rje.dataDict(self,filename,['Sample'],['Identifier'],lists=True)
            for sample in indata:
                idlist = indata[sample]['Identifier']
                idlist.sort(); self.printLog('#READ','Read %d "%s" identifiers' % (len(idlist),sample))
                if self.opt['ExpOnly'] and self.list['Experiments'] and sample not in self.list['Experiments']: continue
                if sample not in self.dict['Samples']: self.dict['Samples'][sample] = []
                self.dict['Samples'][sample] += idlist
                ix += len(idlist)
            ### ~ [3] Summarise Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#DATA','%s: %s protein IDs read for %d samples' % (base,rje.integerString(ix),len(self.dict['Samples'])))
            self.printLog('#SAMPLE','%d %s samples: %s.' % (len(self.dict['Samples']),base,string.join(rje.sortKeys(self.dict['Samples']),'; ')))
        except: self.errorLog('Problem loading data from "%s"' % filename); raise
#########################################################################################################################
    def addBaits(self): ### Adds PPI of baits as new samples (named after bait genes)
        '''Adds PPI of baits as new samples (named after bait genes).'''
        try:
            ### ~ Pull out combined PPI for each bait and add as sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for bait in self.list['Baits']:
                if bait not in self.dict['PPI']:
                    self.printLog('#BAIT','No PPI parsed for bait "%s" - no sample added' % bait)
                    continue
                if self.info['ComBaits'].lower() not in ['','none']: sample = self.info['ComBaits']
                else: sample = '%s-PPI' % bait
                if sample not in self.dict['Samples']: self.dict['Samples'][sample] = []
                for ppi in self.dict['PPI'][bait]:
                    gene = self.geneMap(ppi)
                    if gene not in self.dict['Samples'][sample]: self.dict['Samples'][sample].append(gene)
                self.printLog('#BAIT','%d PPI parsed for bait "%s" - new sample added' % (len(self.dict['PPI'][bait]),bait))
        except: self.errorLog(rje_zen.Zen().wisdom())  
#########################################################################################################################
    def mapIPI(self):   ### Parses mappings from IPI links table to Ens/HGNC                                        |3.0|
        '''Parses mappings from IPI links table to Ens/HGNC.'''
        ### ~ [0] ~ Setup paths and files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.info['IPILinks'].lower() in ['','none']: return
        if not os.path.exists(self.info['IPILinks']):
            self.errorLog('IPI map file "%s" missing' % self.info['IPILinks'],printerror=False)
            return False
        genemap = self.obj['GeneMap']
        ### ~ [1] ~ Parse Sequence IDs Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ipidata = rje.dataDict(self,self.info['IPILinks'])  # This is now a dictionary of IPI:Data
        ix = 0.0
        for ipi in ipidata:
            ## ~ [1a] ~ Parse line of IPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#IPI','Mapping IPI: %.1f%%' % (ix/len(ipidata))); ix += 100.0
            genelist = string.split(ipidata[ipi]['EnsG'],',')
            while '' in genelist: genelist.remove('')
            hgnclist = string.split(ipidata[ipi]['Symbol'].upper(),',')     #!# Upper case GeneList #!#
            while '' in hgnclist: hgnclist.remove('')
            ## ~ [1b] ~ Update GeneMap object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gene in hgnclist + genelist: genemap.addAlias(gene,ipi)
        self.printLog('\r#IPI','Mapping IPI complete: %s proteins.' % rje.integerString(len(ipidata)))
#########################################################################################################################
    def mapEnsEMBL(self):   ### Extracts EnsEMBL mapping data from a BioMart download                               |3.0|
        '''Extracts EnsEMBL mapping data from a BioMart download.'''
        ### ~ [0] ~ Setup paths and files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(self.info['EnsMap']): return self.errorLog('EnsEMBL map file "%s" missing' % self.info['EnsMap'],printerror=False)
        ### ~ [1] ~ Parse Sequence IDs Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [1a] ~ Read in data and adjust for version of download (headers) ~~~~~~~~~~~~~~~~~~~~~ ##
        if open(self.info['EnsMap'],'r').readline().find('HGNC Symbol') > 0:        
            mainkeys = ['Ensembl Gene ID','Ensembl Peptide ID','HGNC Symbol']
        else: mainkeys = ['Ensembl Gene ID','Ensembl Peptide ID','HGNC symbol']     
        if open(self.info['EnsMap'],'r').readline().find('Ensembl Protein ID') > 0: 
            mainkeys[1] = 'Ensembl Protein ID'
        ensdata = rje.dataDict(self,self.info['EnsMap'],mainkeys,mainkeys)
        ## ~ [1b] ~ Extract useful data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        epx = 0
        for map in ensdata:
            self.progLog('\r#ENS','%s EnsEMBL mappings' % rje.integerString(epx)); epx += 1
            gene = ensdata[map]['Ensembl Gene ID']
            try: pept = ensdata[map]['Ensembl Peptide ID']
            except: pept = ensdata[map]['Ensembl Protein ID']
            try: hgnc = ensdata[map]['HGNC Symbol']
            except: hgnc = ensdata[map]['HGNC symbol']
            hgnc = hgnc.upper()     # Make upper case
            self.mapper().addAlias(gene,pept)
            self.mapper().addAlias(hgnc,gene)
        self.printLog('\r#ENS','%s EnsEMBL mappings' % (rje.integerString(epx)))
#########################################################################################################################
    def addAlias(self): ### Adds manual aliases to GeneMap object from non-conventional format file
        '''
        Adds manual aliases to GeneMap object from non-conventional format file. Multiple lines of:
        alias = ID1 ID2 ...
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(self.info['AddAlias']): return
            if not self.obj['GeneMap']: return
            self.obj['GeneMap'].dict['BestMap'] = {}
            ### ~ [2] ~ Add Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in self.loadFromFile(self.info['AddAlias'],chomplines=True):
                ids = string.split(line,' = ')
                if len(ids) < 2: continue
                alias = ids[0]
                for id in string.split(ids[1]): self.obj['GeneMap'].addAlias(id,alias)
            self.printLog('#ALIAS','Added aliases from %s' % self.info['AddAlias'])
        except: self.errorLog('Pingu.addAlias error')
#########################################################################################################################
    def identifiers(self):  ### Returns list of identifiers used for Datasets                                       |3.0|
        '''Returns list of identifiers used for Datasets.'''
        idlist = []
        for sample in self.dict['Datasets']: idlist += self.dict['Datasets'][sample]
        return rje.sortUnique(idlist,xreplace=False)
#########################################################################################################################
    def geneList(self,idlist=[]):     ### Returns full list of genes, expanding families if necessary               |3.0|
        '''
        Returns full list of genes, expanding families if necessary.
        >> idlist:list [] = list of identifiers - could be genes or families
        '''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not idlist: return self.list['Genes'][0:]
        if not self.dict['Families']: return idlist[0:]
        ### ~ [2] ~ Expand families and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        genes = idlist[0:]
        for gene in idlist:
            if gene in self.dict['Families']:
                i = genes.index(gene)
                genes = genes[:i] + self.dict['Families'][gene] + genes[i:]
        genes.sort()
        return genes
#########################################################################################################################
    def fullGeneList(self): return self.geneList(self.identifiers())    ### Returns full gene list from identifiers |3.0|
#########################################################################################################################
    def samples(self): return rje.sortKeys(self.dict['Datasets'])   ### Returns keys of dataset dictionary          |3.0|
#########################################################################################################################
    ### <4> ### NR/HGNC Data processing methods                                                                         #
#########################################################################################################################
    def makeNRGenes(self):  ### Converts (peptide) datasets into NR gene datasets (EnsEMBL)                         |3.0|
        '''Converts peptide datasets into NR gene datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Genes'] = []     # Final NR GeneList
            self.dict['Datasets'] = {}  # Final NR by dataset
            self.dict['BackMap'] = {}   # Mapping back from final NR genes to original proteins
            ### ~ [1] ~ Protein/Gene conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sample in self.dict['Samples']:
                self.dict['Datasets'][sample] = []
                for prot in self.dict['Samples'][sample]:
                    ## Map onto gene ##
                    gene = self.geneMap(prot)
                    ## Check for duplicity ##
                    if gene not in self.dict['BackMap']: self.dict['BackMap'][gene] = [prot]
                    elif prot not in self.dict['BackMap'][gene]:
                        self.printLog('\r#DUP','%s %s = Gene %s - already mapped to %s.' % (sample,prot,gene,string.join(self.dict['BackMap'][gene],', ')))
                        self.dict['BackMap'][gene].append(prot)
                    ## Add to dataset ##
                    if gene not in self.list['Genes']: self.list['Genes'].append(gene)
                    if gene not in self.dict['Datasets'][sample]: self.dict['Datasets'][sample].append(gene)
                ## Finish sample ##
                self.dict['Datasets'][sample].sort()
                self.printLog('#SAMPLE','%s proteins from "%s" mapped onto %s genes.' % (len(self.dict['Samples'][sample]),sample,len(self.dict['Datasets'][sample])))
                if self.opt['GeneLists'] or self.opt['FullOutput']:
                    gfile = '%s.%s.nrgenes.txt' % (self.info['Basefile'],sample)
                    open(gfile,'w').write(string.join(self.dict['Datasets'][sample],'\n'))
                    self.printLog('#GFILE','%d %s genes output to %s' % (len(self.dict['Datasets'][sample]),sample,gfile))
            self.list['Genes'].sort()
        except: self.errorLog('Pingu.makeNRGenes() Error')
#########################################################################################################################
    def convertToHGNC(self):  ### Converts EnsEMBL datasets into NR HGNC gene datasets                              |3.0|
        '''Converts EnsEMBL datasets into NR HGNC gene datasets.'''
        try:### ~ [1] ~ Generate NR EnsEMBL-HGNC conversions (just use one HGNC symbol) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hgncgenes = []
            for gene in self.list['Genes']:
                ## Map Gene onto a single NR HGNC symbol ##
                symbol = self.geneMap(gene)
                ## Update dictionaries ##
                if symbol not in hgncgenes: hgncgenes.append(symbol)
                if symbol not in self.dict['BackMap']: self.dict['BackMap'][symbol] = []
                for mapped in self.dict['BackMap'][gene]:
                    if mapped not in self.dict['BackMap'][symbol]: self.dict['BackMap'][symbol].append(mapped)
                if gene[:4] == 'ENSG':      #!# Add EnsEMBL gene data to GeneMap data. Is this necessary? #!#
                    rje.combineDict(self.mapper().getGeneData(symbol),self.mapper().getGeneData(gene),overwrite=False,replaceblanks=True)
                    rje.combineDict(self.mapper().getGeneData(symbol),{'EnsEMBL':gene},overwrite=False,replaceblanks=True)
                if gene in self.dict['GO']:                                     # Combine EnsEMBL GO mappings 
                    for go in self.dict['GO'][gene]: self.addGeneGO(symbol,go)  # Add HGNC to GO dictionary
            self.printLog('#GENE','%s EnsEMBL genes mapped to %s NR HGNC symbols' % (len(self.list['Genes']),len(hgncgenes)))
            hgncgenes.sort()
            self.list['Genes'] = hgncgenes[0:]                  # Final NR GeneList
            ### ~ [2] ~ Convert Datasets to HGNC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sample in rje.sortKeys(self.dict['Datasets']):
                oldgenes = self.dict['Datasets'].pop(sample)
                self.dict['Datasets'][sample] = []  # Final NR by dataset
                for gene in oldgenes:
                    hgnc = self.geneMap(gene)
                    if hgnc not in self.dict['Datasets'][sample]: self.dict['Datasets'][sample].append(hgnc)
                self.dict['Datasets'][sample].sort()
                self.log.printLog('#HGNC','%s: %d EnsEMBL -> %d HGNC' % (sample,len(oldgenes),len(self.dict['Datasets'][sample])))
        except: self.log.errorLog('Pingu HGNC conversion error.')
#########################################################################################################################
    def geneMap(self,symbol):   ### Gets best gene symbol for gene and updates GeneMap dictionary
        '''Gets best gene symbol for gene and updates GeneMap dictionary.'''
        try: return self.mapper().bestMap(symbol)
        except: self.log.errorLog('Sunnofabitch. Pingu has misbehaved.')
#########################################################################################################################
    def summaryHGNC(self):  ### Generate a summary table of genes in dataset, inc. peptide lists for each sample    |3.0|
        '''
        Generate a summary table of genes in dataset, including peptide lists for each sample. This method does not
        currently incorporate additional GeneMap data but could/should be expanded to do so.
        '''
        try:### ~ [0] ~ Setup output columns and file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s.hgnc.tdt' % self.info['Basefile']
            headers = ['Gene','Aliases','HGNC','Entrez']
            for sample in rje.sortKeys(self.dict['Samples']): headers.append(sample)
            if self.dict['Families']: headers.append('Family')
            headers += ['PPI','UniProt','HPRD','OMIM','EnsEMBL','EnsLoci','Desc']
            rje.delimitedFileOutput(self,hfile,headers,rje_backup=True)
            ### ~ [1] ~ Output data per gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            genes = self.list['Genes'][0:]
            for gene in self.list['Genes']:
                if gene in self.dict['Families']:
                    i = genes.index(gene)
                    genes = genes[:i] + self.dict['Families'][gene] + genes[i:]
            hx = 0.0
            for hgnc in genes:
                try:
                    self.progLog('\r#HGNC','Summary HGNC table: %.1f%%' % (hx/len(self.list['Genes']))); hx += 100.0
                    ## ~ [1a] ~ Basic Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    data = {'Gene':hgnc}    #,'GO_C':0,'GO_F':0,'GO_P':0,'GO_X':0}
                    rje.combineDict(data,self.obj['GeneMap'].getGeneData(hgnc),overwrite=False,replaceblanks=True)
                    if rje.dictValues(data,'EnsDesc',valtype='str'): data['Desc'] = data['EnsDesc']
                    data['Aliases'] = string.join(rje.sortUnique(self.mapper().redundancy(hgnc)),',')
                    if self.dict['Families']:
                        data['Family'] = self.dict['FamMap'][hgnc]
                        if data['Family'] in self.dict['PPI']: data['PPI'] = len(self.dict['PPI'][data['Family']])
                    if hgnc in self.dict['PPI']: data['PPI'] = len(self.dict['PPI'][hgnc])
                    ## ~ [1b] ~ Sample Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for sample in self.dict['Samples']:
                        peplist = []
                        for pep in self.dict['BackMap'][hgnc]:
                            if pep in self.dict['Samples'][sample]: peplist.append(pep)
                        data[sample] = string.join(rje.sortUnique(peplist),',')
                    ## ~ [1c] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    rje.delimitedFileOutput(self,hfile,headers,datadict=data)
                except: self.errorLog('SummaryHGNC Problem with gene "%s" (= "%s")' % (hgnc,self.geneMap(hgnc)))
            self.printLog('\r#HGNC','Summary HGNC table for %d genes: %s' % (len(self.list['Genes']),hfile))                
        except: self.errorLog('Error during Pingu.summaryHGNC()')            
#########################################################################################################################
    def mapOut(self):   ### Generates a summary table of full peptide mapping
        '''Generates a summary table of full peptide mapping.'''
        try:### ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '%s.hgnc_map.tdt' % self.info['Basefile']
            headers = ['Protein','Gene','Redundancy','BackMap']
            if self.dict['FamMap']: headers.append('Family')
            peplist = []
            for sample in rje.sortKeys(self.dict['Samples']):
                headers.append(sample)
                peplist += self.dict['Samples'][sample]
            peplist = rje.sortUnique(peplist)
            rje.delimitedFileOutput(self,mfile,headers,rje_backup=True)

            ### ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            px = 0.0
            for pept in peplist:
                self.log.printLog('\r#MAP','Peptide mapping output: %.1f%%' % (px/len(peplist)),newline=False,log=False)
                px += 100.0
                data = {'Protein':pept,'Gene':self.geneMap(pept),'Family':''}
                data['Gene'] = self.geneMap(pept)
                try:
                    data['BackMap'] = string.join(self.dict['BackMap'][data['Gene']],',')
                    data['Redundancy'] = len(self.dict['BackMap'][data['Gene']])
                except: data['Redundancy'] = 'N/A'
                if data['Gene'] in self.dict['FamMap']: data['Family'] = self.dict['FamMap'][data['Gene']]
                for sample in rje.sortKeys(self.dict['Samples']):
                    data[sample] = ''
                    if pept in self.dict['Samples'][sample]: data[sample] = 'P'
                    else:
                        for alt in self.dict['BackMap'][data['Gene']]:
                            if alt in self.dict['Samples'][sample]: data[sample] = 'G'
                    if not data[sample] and data['Family'] in self.dict['Datasets'][sample]: data[sample] = 'F'
                rje.delimitedFileOutput(self,mfile,headers,datadict=data)
            self.printLog('\r#MAP','Peptide mapping output to %s' % mfile)
        except: self.log.errorLog('Wah, wah, wah! Cry like a baby! Pingu is disturbed')            
#########################################################################################################################
    def makeFam(self):  ### Groups sequences into families and combines PPI                                         |3.0|
        '''Groups sequences into families and combines PPI.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            clusdis = self.stat['MakeFam']
            if clusdis > 1.0: clusdis /= 100.0
            clusdis = 1.0 - clusdis
            if clusdis <= 0.0: return self.printLog('#FAM','MakeFam=%s -> maxdis < 0!' % self.stat['MakeFam'])
            self.dict['Families'] = {}
            self.dict['FamMap'] = {}
            self.dict['FamPPI'] = {}
            ### ~ [2] ~ Combine and name into families ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['Homology'].cluster(clusdis)
            for cluster in self.obj['Homology'].list['Clusters']:
                if len(cluster) == 1: fam = cluster[0]
                else:
                    fgene = cluster[0]
                    while cluster[0][:4] == 'ENSG':
                        cluster = cluster[1:] + cluster[:1]
                        if cluster[0] == fgene: break
                    fam = '%s.fam%d' % (cluster[0],len(cluster))
                    #x#self.deBug('%s: %s' % (fam,cluster))
                self.dict['Families'][fam] = cluster[0:]
                for gene in cluster:
                    self.dict['FamMap'][gene] = fam
                    if fam not in self.dict['BackMap']: self.dict['BackMap'][fam] = []
                    if gene not in self.dict['BackMap'][fam]: self.dict['BackMap'][fam].append(gene)
            #x#self.deBug(self.dict['FamMap'])
            for gene in self.list['Genes']:
                if gene not in self.obj['Homology'].dict['Matrix']:
                    self.dict['Families'][gene] = [gene]
                    self.dict['FamMap'][gene] = gene
            self.printLog('#FAM','%s family clusters from %s genes' % (len(self.dict['Families']),len(self.list['Genes'])))
            ### ~ [3] ~ Combine PPI into Family PPI dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (fx,ftot) = (0.0,len(self.dict['PPI']))
            for gene in rje.sortKeys(self.dict['PPI']):
                self.progLog('\r#PPI','Combining PPI for families: %.1f%%' % (fx/ftot)); fx += 100.0
                try: fam = self.dict['FamMap'][gene]
                except: fam = gene
                if fam not in self.dict['FamPPI']: self.dict['FamPPI'][fam] = []
                for ppi in self.dict['PPI'][gene]:
                    try: ppifam = self.dict['FamMap'][ppi]
                    except: ppifam = ppi
                    if ppifam not in self.dict['FamPPI'][fam]: self.dict['FamPPI'][fam].append(ppifam)
                self.dict['FamPPI'][fam].sort()
            self.printLog('\r#PPI','Combining of PPI for families complete')
            self.dict['PPI'] = self.dict['FamPPI']  #!# Do we want this? Includes all PPI! #!#
            ## ~ [3a] ~ Update Samples and Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fmhead = ['Family','Genes','PPI']
            fmdata = {}
            for fam in rje.sortKeys(self.dict['Families']):
                fmdata[fam] = {'Family':fam, 'Genes':string.join(self.dict['Families'][fam])}
                if fam in self.dict['FamPPI']: fmdata[fam]['PPI'] = len(self.dict['FamPPI'][fam])
                else: fmdata[fam]['PPI'] = 0
                for sample in rje.sortKeys(self.dict['Samples']): fmdata[fam][sample] = []
            for sample in rje.sortKeys(self.dict['Samples']):
                fmhead.append(sample)
                sgenes = self.dict['Datasets'][sample][0:]
                self.dict['Datasets'][sample] = []
                for gene in sgenes:
                    try: ppifam = self.dict['FamMap'][gene]
                    except:
                        print self.dict['FamMap']
                        self.deBug(gene)
                        self.deBug(rje.sortKeys(self.dict['FamMap']))
                        ppifam = gene
                    if ppifam not in self.dict['Datasets'][sample]: self.dict['Datasets'][sample].append(ppifam)
                    if gene not in fmdata[ppifam][sample]: fmdata[ppifam][sample].append(gene)
                self.dict['Datasets'][sample].sort()
            ### ~ [4] ~ Output of family data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fmfile = '%s.FamMap.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,fmfile,fmhead,rje_backup=True)
            for fam in rje.sortKeys(self.dict['Families']):
                for sample in rje.sortKeys(self.dict['Samples']):
                    fmdata[fam][sample].sort()
                    fmdata[fam][sample] = string.join(fmdata[fam][sample])
                rje.delimitedFileOutput(self,fmfile,fmhead,datadict=fmdata[fam])
            self.printLog('#FAM','Family-Gene mapping output to %s' % fmfile)
            fmfile = '%s.Families.tdt' % self.info['Basefile']
            fmhead = ['Sample','Family']
            rje.delimitedFileOutput(self,fmfile,fmhead,rje_backup=True)
            for sample in rje.sortKeys(self.dict['Samples']):
                for fam in self.dict['Datasets'][sample][0:]:
                    rje.delimitedFileOutput(self,fmfile,fmhead,datadict={'Family':fam,'Sample':sample})
            self.printLog('#FAM','Sample families output to %s' % fmfile)
            #x#self.list['Genes'] = rje.sortKeys(self.dict['Families'])
        except: self.log.errorLog('Pingu.process() Error')
#########################################################################################################################
    def remSticky(self,purgeppi=True):  ### Removed "sticky" hub proteins from datasets (and PPI dictionary)        |3.0|
        '''Removed "sticky" hub proteins from datasets (and PPI dictionary if purgeppi=True).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            autorem = self.i() < 0 or not rje.yesNo('Manually verify sticky sample hub removal?')
            sticky = []     # List of sticky genes from PPI
            if self.opt['StickyPPI']: stickcheck = self.list['Genes'][0:]
            else: stickcheck = rje.sortKeys(self.dict['PPI'])
            (px,ptot) = (0.0,len(stickcheck))
            for hub in stickcheck:
                self.progLog('\r#STICKY','Generating "sticky" hub list: %.2f%%' % (px/ptot)); px += 100.0
                if hub not in self.dict['PPI']: continue
                if len(self.dict['PPI'][hub]) > self.stat['RemSticky']:
                    if hub not in self.list['Genes'] or autorem or rje.yesNo('Remove Sticky Hub %s (%d ppi)?' % (hub,len(self.dict['PPI'][hub]))): sticky.append(hub)
            self.printLog('\r#STICKY','Generated "sticky" hub list: %s hubs >%d PPI' % (rje.integerString(len(sticky)),self.stat['RemSticky']))
            ### ~ [1] ~ Purge Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot,rx) = (0.0,len(sticky),0)
            for hub in sticky:
                self.progLog('\r#STICKY','Removing "sticky" hubs from samples: %.2f%%' % (sx/stot)); sx += 100.0
                samples = []
                for sample in rje.sortKeys(self.dict['Samples']):
                    if hub in self.dict['Datasets'][sample]:
                        self.dict['Datasets'][sample].remove(hub)
                        samples.append(sample)
                if samples:
                    self.printLog('\r#STICKY','Removed "sticky" hub gene "%s" from: %s' % (hub,string.join(samples,'; '))); rx += 1
                    if hub in self.list['Genes']: self.list['Genes'].remove(hub)
                    elif hub in self.dict['Families']:
                        for gene in self.dict['Families'].pop(hub): self.list['Genes'].remove(gene)
            self.printLog('\r#STICKY','Removed %d "sticky" hubs from samples (%d remain).' % (rx,len(self.identifiers())))
            ### ~ [2] ~ Purge PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not purgeppi: return self.printLog('\r#STICKY','Not purging PPI of "sticky" hubs.')
            (px,ptot,rx) = (0.0,len(self.dict['PPI']),0)
            for hub in rje.sortKeys(self.dict['PPI']):
                self.progLog('\r#STICKY','Purging PPI of "sticky" hubs: %.2f%%' % (px/ptot)); px += 100.0
                if hub in sticky and not self.opt['StickyHubs']: rx += len(self.dict['PPI'].pop(hub))
                else:
                    for spoke in self.dict['PPI'][hub][0:]:
                        if spoke in sticky: self.dict['PPI'][hub].remove(spoke); rx += 1
            self.printLog('\r#STICKY','Purged %s "sticky" hub PPI (%s hubs >%d PPI).' % (rje.integerString(rx/2),rje.integerString(len(sticky)),self.stat['RemSticky']))            
        except: self.errorLog('Things got sticky during Pingu.remSticky()')
#########################################################################################################################
    ### <5> ### Gene Ontology Analysis Methods                                                                          #
#########################################################################################################################
    def mapEnsGO(self):   ### Extracts EnsEMBL GO mapping data from a BioMart download                              |3.0|
        '''Extracts EnsEMBL GO mapping data from a BioMart download.'''
        self.obj['GO'].mapEnsGO(spec='HUMAN',gokey='EnsGO',fixhead=True)
        self.dict['GO'] = self.obj['GO'].dict['EnsGO']
        self.list['EnsGenes'] = rje.sortKeys(self.dict['GO'])
#########################################################################################################################
    def parseGO(self):  ### Sets up self.obj['GO'] using rje_go
        '''Sets up self.obj['GO'] using rje_go.'''
        try:### ~ Sets up self.obj['GO'] using rje_go ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            self.obj['GO'].readGO()
        except: self.errorLog('Problem setting up rje_go.GO object')
#########################################################################################################################
    def getGO(self,gene):   ### Gets the full list of GO ID for a given gene (or protein) using EnsGO and GeneMap
        '''Gets the full list of GO ID for a given gene (or protein) using EnsGO and GeneMap.'''
        try:### ~ [1] ~ Look for gene, then try symbol, then try EnsEMBL mapped to symbol ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gene in self.dict['GO']: return self.dict['GO'][gene]        # Gene found in dictionary
            gene = self.geneMap(gene)                                       # Try mapping to NR gene/HGNC
            if gene in self.dict['GO']: return self.dict['GO'][gene]        # Mapped gene found in dictionary
            if not self.obj['GeneMap'].getGeneData(gene): return []         # No mapping to EnsEMBL via GeneMap
            ens = rje.dictValues(self.obj['GeneMap'].getGeneData(gene),'EnsEMBL',valtype='str')     # Try to map EnsEMBL
            if ens and ens in self.dict['GO']: return self.dict['GO'][ens]  # Found EnsEMBL mapping in GO dictionary
        except: self.errorLog('Pingu.getGO(%s) error' % gene)
        return []                                                           # Return empty list if failed
#########################################################################################################################
    def filterGOSamples(self):  ### Filters genes from samples based on GO annotation                               |3.0|
        '''Filters genes from samples based on GO annotation.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            filtered = self.filterGO(self.list['Genes'],'Sample')
            if not filtered: return
            gofilter = self.list['GOFilter']
            goexcept = self.list['GOExcept']
            ### ~ [2] ~ Output list of filtered genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            GFILE = open('%s.gofilter.txt' % self.info['Basefile'],'w')
            GFILE.write('### GO Filter ###\n%s\n\n' % string.join(gofilter,'\n'))
            if goexcept: GFILE.write('### GO Exceptions ###\n%s\n\n' % string.join(goexcept,'\n'))
            GFILE.write('### Filtered Genes ###\n')
            for gene in rje.sortKeys(filtered):
                try: desc = self.obj['GeneMap'].getGeneData(gene)['Desc']
                except: desc = '???'
                GFILE.write('%s\t-\t%s\t-\t%s\n' % (gene,filtered[gene],desc))
            GFILE.close()
            ### ~ [3] ~ Filter out relevant genes from samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gene in filtered:
                self.list['Genes'].remove(gene)
                for sample in rje.sortKeys(self.dict['Datasets']):
                    if gene in self.dict['Datasets'][sample]: self.dict['Datasets'][sample].remove(gene)
        except: self.errorLog('Pingu.filterGOSamples() error'); raise
#########################################################################################################################
    def filterGOPPI(self):      ### Filters genes from PPI dictionary based on GO annotation                        |3.0|
        '''Filters genes from samples based on GO annotation.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            filtered = self.filterGO(rje.sortKeys(self.dict['PPI']),'PPI')
            if not filtered: return 
            ### ~ [2] ~ Filter out relevant PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (px,ptot,fx) = (0.0,len(self.dict['PPI']),0)
            for hub in rje.sortKeys(self.dict['PPI']):
                self.progLog('\r#PPI','Filtering PPI based on GO annotation: %.2f%%' % (px/ptot)); px += 100.0
                if hub in filtered: fx += len(self.dict['PPI'].pop(hub)); continue
                for spoke in self.dict['PPI'][hub][0:]:
                    if spoke in filtered:
                        self.dict['PPI'][hub].remove(spoke)
                        fx += 1
            self.printLog('\r#PPI','%s PPI filtered based on GO annotation of interactors.' % (rje.integerString(fx/2)))
        except: self.errorLog('Pingu.filterGOPPI() error'); raise
#########################################################################################################################
    def filterGO(self,genelist,gtext):     ### Filters genes based on GO annotation                                 |3.0|
        '''
        Filters genes based on GO annotation.
        >> genelist:list of genes to filter (e.g. self.list['Genes'] or self.dict['PPI'] keys)
        >> gtext:str = Description of genelist for filtering
        << filtered:dict of {filtered gene:Reason for filtering}
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not genelist: return {}
            gofilter = self.list['GOFilter']
            goexcept = self.list['GOExcept']
            if not gofilter: self.printLog('#FILTER','No %s GO filtering.' % gtext); return {}
            safe = []       # List of "safe" genes (from GOExcept)
            filtered = {}   # Dictionary of {gene:Reason}
            ### ~ [2] ~ Except ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if goexcept:
                (gx, gtot) = (0.0, len(self.list['Genes']))
                for gene in self.list['Genes']:
                    self.progLog('\r#GOFILT','Checking genes for GO filtering exceptions: %.1f%%' % (gx/gtot)); gx += 100
                    for goid in self.getGO(gene):
                        go = self.go(goid)['name']
                        for skipgo in goexcept:
                            if skipgo == goid or go.lower().find(skipgo.lower()) >= 0: safe.append(gene); break
                        if gene in safe: break
                self.printLog('\r#GOFILT','Checked genes for GO filtering exceptions: %d "safe"' % (len(safe)))
            ### ~ [3] ~ Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (gx, gtot) = (0.0, len(self.list['Genes']))
            for gene in self.list['Genes']:
                self.progLog('\r#GOFILT','Filtering genes on GO terms: %.1f%%' % (gx/gtot)); gx += 100
                if gene in safe: continue
                for goid in self.getGO(gene):
                    go = self.go(goid)['name']
                    for badgo in gofilter:
                        if badgo == goid or go.lower().find(badgo.lower()) >= 0: filtered[gene] = '%s %s' % (goid,go); break
                    if gene in filtered: break
            self.printLog('\r#GOFILT','Filtered %d genes using GO terms.' % (len(filtered)))
        except: self.errorLog('Pingu.filterGO() error')
#########################################################################################################################
    def countGO(self,subset=[]):    ### Returns GO count dictionary for subset (or all if subset missing)
        '''
        Returns GO count dictionary for subset (or all if subset missing).
        >> subset:list of GO IDs. Will use self.go() if empty list
        << dictionary of {ID+['n','bp','cc','mf']:Count}
        '''
        try:### ~ [1] ~ Setup dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pgo = {}                # Dictionary of {GO ID:{Sample:Count}}
            if not subset: subset = self.go().keys()
            for id in subset + ['n','bp','cc','mf']:
                pgo[id] = {}
                for sample in self.samples() + ['Pingu','EnsEMBL']: pgo[id][sample] = 0
            ### ~ [2] ~ Add samples/pingu to dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            allgenes = self.fullGeneList()  #x#self.list['Genes'][0:]
            self.deBug(allgenes)
            (sx,stot) = (0.0,len(self.list['EnsGenes']))
            for sample in self.dict['Samples']: stot += len(self.geneList(self.dict['Datasets'][sample]))
            for sample in self.dict['Samples']:
                for gene in self.geneList(self.dict['Datasets'][sample]):
                    self.progLog('#GOX','Counting GO term occurrences for %s IDs: %.1f%%' % (rje.iLen(subset),sx/stot))
                    sx += 100.0
                    types = {'n':1,'bp':0,'cc':0,'mf':0}
                    for id in self.getGO(gene):         # Take each GO ID in turn
                        if id not in subset: continue   # Check part of subset
                        go = self.go(id)                # Get GO data
                        types[go['type']] = 1           # This gene has 1+ terms of this type
                        pgo[id][sample] += 1            # Add 1 to sample count
                        if gene in allgenes: pgo[id]['Pingu'] += 1  # Count for whole pingu set too
                    for type in types:                              # Update each total
                        pgo[type][sample] += types[type]
                        if gene in allgenes: pgo[type]['Pingu'] += types[type]
                    #self.deBug('%s:%s' % (gene, gene in allgenes))
                    if gene in allgenes: allgenes.remove(gene)      # Only count once
            ### ~ [3] ~ Add EnsEMBL data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('\r#GOX','Counting GO term occurrences for %s IDs complete. (Bypassed EnsEMBL.)' % (rje.iLen(subset)),log=False)
            #self.printLog('#ZEN',rje_zen.Zen().wisdom())
            #return pgo
            sample = 'EnsEMBL'
            for gene in self.list['EnsGenes']:     # Dictionary made during self.mapEnsEMBL()
                self.progLog('#GOX','Counting GO term occurrences for %s IDs: %.1f%%' % (rje.iLen(subset),sx/stot))
                sx += 100.0
                if gene[:4] != 'ENSG': continue                     # Not an EnsEMBL gene
                types = {'n':1,'bp':0,'cc':0,'mf':0}
                for id in self.getGO(gene):                         # Take each GO ID in turn
                    if id not in subset: continue                   # Check part of subset
                    go = self.go(id)                                # Get GO data
                    types[go['type']] = 1                           # This gene has 1+ terms of this type
                    pgo[id][sample] += 1                            # Add 1 to sample count
                for type in types: pgo[type][sample] += types[type] # Update each total
            ### ~ [4] ~ Finish & Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#GOX','Counting GO term occurrences for %s IDs complete.' % (rje.iLen(subset)),log=False)
            self.deBug(pgo['n'])
            return pgo
        except: self.errorLog('Pingu.countGO Error'); raise
#########################################################################################################################
    def sampleGO(self):     ### Performs new sample GO analysis using rje_go mappings etc.                          |3.0|
        '''Performs new sample GO analysis using rje_go mappings etc.'''
        try:### ~ [1] ~ Add Samples and Counts to GO object dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['EnsGenes']: self.mapEnsGO()
            self.deBug(len(self.list['EnsGenes']))
            pgo = {'go':self.countGO()}                                 # Dictionary of {Subset:{GO ID:{Sample:Count}}}
            for id in self.go(): self.go(id)['pingu'] = pgo['go'][id]   # Add counts to GO dictionary

            ### ~ [2] ~ Generate Additional GO slim datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            subset = self.obj['GO'].powerGO(pgo['go']['n'],sig=0.001,samples='all',total='Pingu',countkey='pingu',ignore=['EnsEMBL'])
            if subset: self.obj['GO'].dict['Subset']['goslim_power'] = {'name':'Power-derived GO slim','terms':subset}
            subset = self.obj['GO'].topTerms(slimx=15,total='Pingu',countkey='pingu')
            if subset: self.obj['GO'].dict['Subset']['goslim_top15'] = {'name':'Top15 GO slim','terms':subset}
            for subset in ['goslim_generic','goslim_user','goslim_power','goslim_top10']:
                if subset in self.obj['GO'].dict['Subset']: pgo[subset] = self.countGO(self.obj['GO'].dict['Subset'][subset]['terms'])
            self.obj['GO'].dict['Subset']['go'] = rje.sortKeys(self.go())   # This "subset" contains All GO terms
            
            ### ~ [3] ~ Output full and GO slim files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for subset in rje.sortKeys(pgo):
                self.deBug('%s::%s' % (subset,pgo[subset]['n']))
                if not pgo[subset]['n']['EnsEMBL']: self.printLog('#ENSGO','No EnsEMBL GO total (n) for "%s"' % subset)
                (dx,dtot) = (0.0,len(pgo[subset]))
                ## ~ [3a] ~ Setup file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gfile = '%s.%s.tdt' % (self.info['Basefile'],subset)
                headers = ['GO','Desc','Type','EnsEMBL','Pingu'] + rje.sortKeys(pgo) + rje.sortKeys(self.dict['Samples'])
                headers.remove('go')
                for sample in rje.sortKeys(self.dict['Samples']): headers.append('e_%s' % sample)
                for sample in rje.sortKeys(self.dict['Samples']): headers.append('p_%s' % sample)
                for sample in rje.sortKeys(self.dict['Samples']): headers.append('ense_%s' % sample)
                for sample in rje.sortKeys(self.dict['Samples']): headers.append('ensp_%s' % sample)
                rje.delimitedFileOutput(self,gfile,headers,rje_backup=True)
                ## ~ [3b] ~ Output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gotypes = {'bp':'biological_process','cc':'cellular_component','mf':'molecular_function','n':'genes'}
                N = pgo[subset]['n']['Pingu']
                for id in rje.sortKeys(pgo[subset]):
                    try: 
                        self.progLog('#GOTDT','%s output: %.1f%%' % (subset,dx/dtot))
                        dx += 100.0
                        if id in ['n','bp','mf','cc']: datadict = {'GO':id,'Desc':gotypes[id],'Type':id}
                        else:
                            go = self.go(id)
                            datadict = {'GO':'GO:%s' % id,'Desc':go['name'],'Type':go['type']}
                        for subgo in rje.sortKeys(pgo):
                            if subgo != 'go' and id in pgo[subgo]: datadict[subgo] = len(pgo[subgo]) - 4
                            else: datadict[subgo] = 0
                        for sample in pgo[subset][id]: datadict[sample] = pgo[subset][id][sample]
                        for sample in self.samples():
                            n = pgo[subset]['n'][sample]
                            k = pgo[subset][id][sample]
                            p = float(pgo[subset][id]['Pingu']) / N
                            datadict['e_%s' % sample] = p * n
                            datadict['p_%s' % sample] = min(rje.binomial(k,n,p,callobj=self),(1 - rje.binomial(k+1,n,p,callobj=self)))
                            if pgo[subset]['n']['EnsEMBL']:
                                ep = float(pgo['go'][id]['EnsEMBL']) / pgo[subset]['n']['EnsEMBL']
                                datadict['ense_%s' % sample] = ep * n
                                datadict['ensp_%s' % sample] = min(rje.binomial(k,n,ep,callobj=self),(1 - rje.binomial(k+1,n,ep,callobj=self)))
                        rje.delimitedFileOutput(self,gfile,headers,datadict=datadict)
                    except: self.errorLog('GO Problem %s::%s' % (id,go))
                self.printLog('\r#GOTDT','GO output for %s terms to %s' % (rje.iStr(dtot),gfile))
        except: self.errorLog('Pingu.sampleGO has had a little trouble')            
#########################################################################################################################
    def goSeqFiles(self):   ### Outputs sequence files for GO categories
        '''Outputs sequence files for GO categories.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['GOSeqDir'].lower() in ['none','','none/','none\\']: return
            if not self.obj['GeneMap']: return self.errorLog('Cannot map EnsLoci without GeneMap.', printerror=False)
            rje.mkDir(self,self.info['GOSeqDir'])
            #x#if self.info['EnsDat'].lower() not in ['none','','none/','none\\']: return self.goDatFiles()
            ## ~ [1a] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapper = self.obj['GeneMap']
            try: ensloci = mapper.obj['EnsLoci']
            except:
                try:
                    ensloci = mapper.info['EnsLoci']
                    if ensloci.seqNum() < 1: raise ValueError
                except:
                    mapper.loadEnsLoci(self.info['EnsLoci'])
                    ensloci = mapper.obj['EnsLoci']
            if not ensloci: return self.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            ensloci.log = self.log
            ### ~ [2] Generate Sequence Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goseq = {}
            (gx,gtot) = (0.0,len(mapper.dict['Data']))
            for gene in rje.sortKeys(mapper.dict['Data']):
                self.progLog('\r#GOSEQ','Mapping GO to EnsLoci: %.2f%%' % (gx/gtot))
                gx += 100
                ## ~ [2a] ~ Check for GO terms and EnsLoci sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if gene not in self.dict['GO']: continue
                try: seq = mapper.dict['EnsLoci'][mapper.dict['Data'][gene]['EnsLoci']]
                except: continue
                ## ~ [2b] ~ Add sequence to Go lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for go in self.dict['GO'][gene]:
                    if go not in goseq: goseq[go] = []
                    if seq not in goseq[go]: goseq[go].append(seq)
            gtot = len(goseq)
            self.printLog('\r#GOSEQ','Mapping GO to EnsLoci: %s terms with 1+ sequences' % rje.integerString(gtot))
            ### ~ [3] Output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for go in rje.sortKeys(goseq):
                sfile = self.info['GOSeqDir'] + 'GO_%s.fas' % go
                try:
                    if self.opt['AccOnly']: ensloci.saveAcc(seqs=goseq[go],accfile='%s.acc' % sfile[:-4])
                    else: ensloci.saveFasta(seqs=goseq[go],seqfile=sfile)
                except: self.errorLog(rje_zen.Zen().wisdom())
            self.printLog('\r#GOSEQ','Mapping GO to EnsLoci: %s files.' % rje.integerString(gtot))
        except: self.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    ### <6> ### PPI Database Methods                                                                                    #
#########################################################################################################################
    def addPPI(self): ### Parses interaction data from simple pairwise file
        '''Parses interaction data from simple pairwise file.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(self.info['AddPPI']): return
            db = rje.baseFile(self.info['AddPPI'],strip_path=True)
            ### ~ [1] Read PPI from AddPPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('\r#PPI','Processing %s: 0.00%%' % (self.info['AddPPI']))
            pdata = rje.dataDict(self,self.info['AddPPI'],['IDA','IDB','Evidence'],datakeys=['IDA','IDB','Evidence'])
            (px,ptot) = (0.0,len(pdata))
            for pkey in pdata.keys()[0:]:
                self.progLog('\r#PPI','Processing %s: %.2f%%' % (self.info['AddPPI'],(px/ptot))); px += 100.0
                idata = pdata.pop(pkey)
                p1 = self.geneMap(idata.pop('IDA'))
                p2 = self.geneMap(idata.pop('IDB'))
                if p1 not in self.dict['PPI']: self.dict['PPI'][p1] = []
                if p2 not in self.dict['PPI'][p1]: self.dict['PPI'][p1].append(p2)
                if p2 not in self.dict['PPI']: self.dict['PPI'][p2] = []
                if p1 not in self.dict['PPI'][p2]: self.dict['PPI'][p2].append(p1)
                if p1 not in self.dict['Evidence']: self.dict['Evidence'][p1] = {}
                if p2 not in self.dict['Evidence']: self.dict['Evidence'][p2] = {}
                if p2 not in self.dict['Evidence'][p1]: self.dict['Evidence'][p1][p2] = []
                if p1 not in self.dict['Evidence'][p2]:self.dict['Evidence'][p2][p1] = []
                for type in string.split(idata['Evidence'],'|'):
                    if (db,type) not in self.dict['Evidence'][p1][p2]: self.dict['Evidence'][p1][p2].append((db,type))
                    if (db,type) not in self.dict['Evidence'][p2][p1]: self.dict['Evidence'][p2][p1].append((db,type))
            self.printLog('\r#PPI','Processing %s: %s hubs; %s PPI read.' % (self.info['AddPPI'],rje.integerString(len(self.dict['PPI'])),rje.integerString(ptot/2)))
        except: self.log.errorLog('Pingu.addPPI() my arse!')
#########################################################################################################################        
    def parsePairwisePPI(self): ### Parses interaction data from pingu pairwise file
        '''Parses interaction data from pingu pairwise file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['PPI'] = {}
            self.dict['Evidence'] = {}
            map = self.obj['GeneMap']
            ### ~ [2] Read PPI from Pairwise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('\r#PPI','Processing %s: 0.00%%' % (self.info['Pairwise']))
            if open(self.info['Pairwise'],'r').readline().find('SpokePPI'):
                pdata = rje.dataDict(self,self.info['Pairwise'],['Hub','Spoke'],datakeys=['Hub','HubUni','Spoke','SpokePPI','SpokeUni','SpokeSeq','Evidence'],headers=['Hub','HubUni','Spoke','SpokePPI','SpokeUni','SpokeSeq','Evidence'])
            else: pdata = rje.dataDict(self,self.info['Pairwise'],['Hub','Spoke'],datakeys=['Hub','HubUni','Spoke','SpokeUni','SpokeSeq','Evidence'],headers=['Hub','HubUni','Spoke','SpokeUni','SpokeSeq','Evidence'])
            (px,ptot) = (0.0,len(pdata))
            for pkey in pdata.keys()[0:]:
                self.progLog('\r#PPI','Processing %s: %.2f%%' % (self.info['Pairwise'],(px/ptot))); px += 100.0
                idata = pdata.pop(pkey)
                p1 = idata.pop('Hub')
                p2 = idata.pop('Spoke')
                if p1 not in self.dict['PPI']: self.dict['PPI'][p1] = []
                if p2 not in self.dict['PPI'][p1]: self.dict['PPI'][p1].append(p2)
                if p2 not in self.dict['PPI']: self.dict['PPI'][p2] = []
                if p1 not in self.dict['PPI'][p2]: self.dict['PPI'][p2].append(p2)
                if p1 not in self.dict['Evidence']: self.dict['Evidence'][p1] = {}
                if p2 not in self.dict['Evidence']: self.dict['Evidence'][p2] = {}
                if p2 not in self.dict['Evidence'][p1]: self.dict['Evidence'][p1][p2] = []
                if p1 not in self.dict['Evidence'][p2]:self.dict['Evidence'][p2][p1] = []
                for evidence in string.split(idata['Evidence'],'|'):
                    try:
                        [db,type] = string.split(evidence,':')
                        if (db,type) not in self.dict['Evidence'][p1][p2]: self.dict['Evidence'][p1][p2].append((db,type))
                        if (db,type) not in self.dict['Evidence'][p2][p1]: self.dict['Evidence'][p2][p1].append((db,type))
                    except: self.errorLog('Problem with %s<->%s evidence "%s"' % (p1,p2,evidence))
                if idata['HubUni'] not in ['','-']:
                    map.addData(p1,{'UniProt':idata['HubUni']})
                    map.addAlias(p1,idata['HubUni'])            
                if idata['SpokeUni'] not in ['','-']:
                    map.addData(p2,{'UniProt':idata['SpokeUni']})
                    map.addAlias(p2,idata['SpokeUni'])            
                if idata['SpokeSeq'] not in ['','-']:
                    map.addData(p2,{'EnsLoci':idata['SpokeSeq']})
                    map.addAlias(p2,idata['SpokeSeq'])            
            self.printLog('\r#PPI','Processing %s: %s hubs; %s PPI read.' % (self.info['Pairwise'],rje.integerString(len(self.dict['PPI'])),rje.integerString(ptot/2)))
        except: self.log.errorLog('Pingu.parsePairwisePPI() my arse!')
#########################################################################################################################        
    def parsePPI(self):     ### Parses PPI Database files into objects
        '''Parses PPI Database files into objects.'''
        try:
            ### ~ [1] Setup HPRD-type databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['HPRD'].lower() not in ['','none']:
                hcmd = ['ppitype=','badtype=','hprdfas=F','alliso=F','domainfas=F','complexfas=F']    # Defaults
                hcmd = hcmd + self.cmd_list + ['genecards=F','hprdpath=%s' % self.info['HPRD']]
                self.dict['PPIDB']['HPRD'] = rje_hprd.HPRD(self.log,hcmd)
            ### ~ [2] Setup BioGRID-type databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ensloci = self.getEnsLoci()
            ensseq = ensloci.seq[0:]
            for db in ['BioGrid','MINT','IntAct','Reactome','DIP','Domino']:
                if self.info[db].lower() in ['','none']: continue
                elif not rje.checkForFile(self.info[db]):
                    self.errorLog('Cannot find %s: %s!' % (db,self.info[db]),printerror=False,quitchoice=True)
                dcmd = ['ppitype=','badtype=indirect_complex,neighbouring_reaction','symmetry=T','species=human','ppifas=F','ppitab=F','alltypes=F']  # Defaults
                dcmd = dcmd + self.cmd_list + ['ppifile=%s' % self.info[db],'seqin=None','dbsource=%s' % db.lower()]
                #x#dcmd.append('seqin=%s' % self.info['EnsLoci']) #?#
                self.dict['PPIDB'][db] = rje_biogrid.BioGRID(self.log,dcmd)
                self.dict['PPIDB'][db].obj['SeqList'] = ensloci
            ### ~ [3] Run objects to parse in data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in rje.sortKeys(self.dict['PPIDB']):
                #X# self.dict['PPIDB'][db].run()   # Parse data + any extra execution
                self.dict['PPIDB'][db].parse()
                #X#if db != 'HPRD': self.dict['PPIDB'][db].updateGenes(genecards=self.obj['GeneMap'],species='HUMAN')
            ensloci.seq = ensseq[0:]    # Remove extra sequences added during PPI parsing!
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def convertPPI(self,species=['HUMAN']):   ### Converts PPI objects into HGNC symbol-based interactions rather than protein IDs
        '''
        Converts PPI objects into HGNC symbol-based interactions rather than protein IDs.
        BioGRID dictionaries used:
        - Protein = Data for each protein {ID:{Seq object,Gene,Aliases}}
        - PPI = Interaction data {ID:[IDs]}
        - Mapping = Dictionary of backwards mapping {Link:[IDs]}
        HPRD Dictionaries used:
        - HPRD = Data for each HPRD protein {HPRD_ID:{Seq object,GB,UniProt,Gene,Desc,Entrez,OMIM}}
        - PPI = Interaction data {HPRD_ID:[HPRD_IDs]}
        - Mapping = Dictionary of backwards mapping {Gene:HPRD_ID}  # Could try going through GB,UniProt,Gene or Entrez? #
        >> species:list ['HUMAN'] = species codes to check for during UniProt extraction
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map = self.obj['GeneMap']
            try: uniprot = self.obj['UniProt']
            except: uniprot = self.obj['UniProt'] = rje_uniprot.UniProt(self.log,self.cmd_list)
            missing = []

            ### ~ [1] ~ Map all possible protein IDs and make a missing list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in rje.sortKeys(self.dict['PPIDB']): missing += self.geneMapPPI(db)

            ### ~ [2] ~ Look for missing in UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mx, mnum) = (0.0,len(missing))
            rejected = []   # List of entries rejected as wrong species: remove from ppi dictionaries!
            if missing:     
                ## ~ [2a] Extract missing sequences from UniProt, where possible ~~~~~~~~~~~~~~~~~~ ##
                self.log.printLog('#ACC','Looking for UniProt AccNum...',newline=False,log=False)
                accdict = uniprot.accDict(acc_list=missing)
                for acc in missing[0:]:
                    mx += 100.0
                    if acc in accdict:
                        entry = accdict[acc]
                        missing.remove(acc)
                        if species and not entry.isSpecies(speclist=species):
                            #x#self.printLog('#SPEC','Rejected %s - not a %s sequence!' % (acc,species),screen=False)
                            rejected.append(acc)
                            continue
                        symbol = None
                        for db in ['Symbol','EnsG','HGNC','Entrez']:
                            if db not in entry.dict['DB']: continue
                            for link in entry.dict['DB'][db]:
                                if db in ['HGNC','Entrez']: g1 = '%s%s' % (db.upper(),link)
                                else: g1 = link
                                if db == 'EnsG': map.addData(g1,{'UniProt':acc,'EnsEMBL':g1})
                                else: map.addData(g1,{'UniProt':acc,db:g1})
                                map.addAlias(g1,acc)
                    self.progLog('\r#SEQ','Mapping missing proteins to UniProt: %.1f%%; %s missing' % ((mx/mnum),rje.integerString(len(missing))))
                self.printLog('\r#SEQ','Mapping missing proteins to UniProt complete: %s missing' % (rje.integerString(len(missing))))
            if rejected:
                rejfile = '%s.spec-rejected.txt' % self.info['Basefile']
                open(rejfile,'w').write('%s\n' % string.join(rejected,'\n'))
                self.printLog('#SPEC','%s sequences rejected - not %s: see %s' % (rje.integerString(len(rejected)),species,rejfile))

            ### ~ [3] ~ Update TempData given new links from databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map.checkTempAliases()
            
            ### ~ [4] ~ Convert PPI using mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in rje.sortKeys(self.dict['PPIDB']): self.geneMapConvertPPI(db,missing,rejected)

        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def getEvidence(self,g1,g2):    ### Gets evidence for interaction as string
        '''Gets evidence for interaction as string.'''
        try:### ~ [1] Get Evidence, join and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elist = []
            if g1 in self.dict['Evidence'] and g2 in self.dict['Evidence'][g1]:
                try: elist = string.split(self.dict['Evidence'][g1][g2],'|')
                except: elist = self.dict['Evidence'][g1][g2]
                for ev in elist:
                    if ev not in elist: elist.append(ev)
            if g2 in self.dict['Evidence'] and g1 in self.dict['Evidence'][g2]:
                try: elist = string.split(self.dict['Evidence'][g2][g1],'|')
                except: elist = self.dict['Evidence'][g2][g1]
                for ev in elist:
                    if ev not in elist: elist.append(ev)
            ## ~ [1a] Convert if not already converted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in range(len(elist)):
                if len(elist[i]) == 2: elist[i] = string.join(elist[i],':')
            if elist: return string.join(elist,'|')
            else: return 'Unknown'
        except KeyboardInterrupt: raise
        except SystemExit: raise
        except: return 'Unknown'            
#########################################################################################################################
    def evidence(self):  ### If appropriate, will output a summary of the interaction types read in
        '''If appropriate, will output a summary of the interaction types read in.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Evidence']: return
            headers = ['Interaction_Type'] + rje.sortKeys(self.dict['PPIDB'])
            efile = '%s%s.evidence.tdt' % (self.info['ResDir'],self.info['Basefile'])
            edata = {}
            if self.list['BadType']: self.list['BadType'] = string.split(string.join(self.list['BadType'],';').lower(),';')
            if self.list['PPIType']: self.list['PPIType'] = string.split(string.join(self.list['PPIType'],';').lower(),';')
            self.deBug('BadType: %s' % self.list['BadType'])
            self.deBug('PPIType: %s' % self.list['PPIType'])
            ## ~ [1a] ~ Load conversion data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            conversion = {}
            if self.info['Evidence'].lower() not in ['','none']:
                convert = rje.dataDict(self,self.info['Evidence'],['Interaction_Type'],['Grouping']+headers)
                for itype in convert:
                    group = convert[itype].pop('Grouping')
                    if group.lower() in ['?','*minor*']: group = itype
                    for db in convert[itype]:
                        if not convert[itype][db]: continue
                        try: string.atoi(convert[itype][db])
                        except: conversion[convert[itype][db]] = group
                self.printLog('#ITYPE','Converted %s evidence codes to groups' % len(conversion))
                headers.append('Grouping')
            self.deBug('Conversion: %s' % conversion)
            ### ~ [2] ~ Generate data dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (px,ptot) = (0.0,len(self.dict['Evidence']))
            for p1 in self.dict['Evidence']:
                self.progLog('\r#ITYPE','Compiling interaction evidence codes: %.2f%%' % (px/ptot)); px += 100
                for p2 in self.dict['Evidence'][p1]:
                    converted = []
                    ## ~ [2a] ~ Update data for summary output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for (db,type) in self.dict['Evidence'][p1][p2]:
                        if type not in edata:
                            edata[type] = {'Interaction_Type':type}
                            if type in conversion: edata[type]['Grouping'] = conversion[type]
                        if db in edata[type]: edata[type][db] += 1
                        else: edata[type][db] = 1
                    ## ~ [2b] ~ Reformat and filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if type in conversion: self.deBug('%s >> %s' % (type,conversion[type])); type = conversion[type]
                        if type.lower() in self.list['BadType'] or (self.list['PPIType'] and type.lower() not in self.list['PPIType']): self.deBug('Not %s %s' % (db,type)); continue
                        type = '%s:%s' % (db,type)
                        if type not in converted: converted.append(type)
                    ## ~ [2c] ~ Update evidence and PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.dict['Evidence'][p1][p2] = converted[0:]
                    self.deBug('%s:%s = %s' % (p1,p2,self.dict['Evidence'][p1][p2]))
            self.progLog('\r#ITYPE','Compiling interaction evidence codes: %.2f%%' % (px/ptot))
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.delimitedFileOutput(self,efile,headers,rje_backup=True)
            for type in rje.sortKeys(edata): rje.delimitedFileOutput(self,efile,headers,datadict=edata[type])
            self.printLog('\r#ITYPE','%d interaction evidence codes output to %s' % (len(edata),efile))
            
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def geneMapPPI(self,db):  ### Uses GeneMap object and UniProt to extract gene and alias information
        '''
        Uses GeneMap object and UniProt to extract gene and alias information.
        >> db:str = Key of self.dict['PPIDB']
        << returns list of "missing" (unmappable) protein IDs
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Summarise and clarify existing data storage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppidb = self.dict['PPIDB'][db]  ##> rje_biogrid or rje_hprd object
            ppi = ppidb.dict['PPI']         ##> PPI data is currently stored as a dictionary of PROTEIN ids 
            map = 'Protein'
            if db in ['HPRD']: map = 'HPRD'
            dbdata = ppidb.dict[map]        ##> Each PROTEIN has additional 'Gene' (str) and 'Alias' (list) data stored
            ppx = rje.integerString(len(ppi))   # No. of hubs in original dictionary
            ## ~ [0b] ~ Tools for mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            map = self.obj['GeneMap']
            ## ~ [0c] ~ New data storage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                
            missing = []                    ##> Want to keep track of missing protein IDs to work out how to improve
            ## ~ [0d] ~ Debugging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#self.deBug('%s : %s' % (rje.sortKeys(ppi),db))
            rlist = [] #x# rje.randomList(ppi)[:10]
            if db.lower() in ['biogrid','mint']: rlist = rlist[:3]

            ### ~ [1] Convert Proteins to Unified Gene Set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for p1 in rje.sortKeys(ppi):    # Take each hub in turn, alphabetically
                self.log.printLog('\r#MAP','Mapping %s %s IDs: %s missing' % (ppx,db,rje.integerString(len(missing))),newline=False,log=False)
                ## ~ [1a] ~ Adjust HPRD to match other DB data structures ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if db in ['HPRD']:
                    try: dbdata[p1]['Alias'] = []
                    except: dbdata[p1] = {'Alias':[]}
                    try: dbdata[p1]['Alias'].append(dbdata[p1]['sp'])
                    except: pass  # There may be trouble ahead.... self.errorLog('No SP entry for HPRD %s' % p1)
                    try: dbdata[p1]['Alias'].append(dbdata[p1]['gb'].upper())
                    except: pass  # There may be trouble ahead.... self.errorLog('No GB entry for HPRD %s' % p1)
                    try: dbdata[p1]['Alias'].append('ENTREZ%s' % dbdata[p1]['entrez'])
                    except: pass  # There may be trouble ahead.... self.errorLog('No ENTREZ entry for HPRD %s' % p1)
                    dbdata[p1]['Alias'].append('HPRD%s' % p1)
                if p1 in rlist:
                    self.deBug(dbdata[p1])
                    self.deBug(ppi[p1])
                ## ~ [1b] ~ Get gene for p1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                g1 = rje.getFromDict(dbdata[p1],'Gene',returnkey=False,case=False,default='')
                if db in ['HPRD']: alias = ['HPRD%s' % p1]
                else: alias = [p1]
                if g1: alias.append(g1)
                alias += rje.getFromDict(dbdata[p1],'Alias',returnkey=False,case=False,default=[])
                for id in alias:
                    id = string.join(string.split(id,':'),'')
                    kids = map.getKeyID(id)     # List of GeneMap IDs matching given alias
                    if kids: break              # Successful mappings of ID
                ## ~ [1c] ~ Store linked IDs for later replacement, or add to missing ~~~~~~~~~~~~~ ##
                if p1 in rlist: self.deBug('%s >> %s >> %s' % (p1,alias,kids))
                if kids:
                    for k in kids:
                        for a in alias: map.addAlias(k,a)
                else: missing.append(p1)
            self.printLog('\r#MAP','Mapping %s %s IDs: %s missing' % (ppx,db,rje.integerString(len(missing))),log=False)
        except: self.log.errorLog('Crap sandwich. Pingu.geneMapPPI fuck-up!')
        return missing
#########################################################################################################################
    def geneMapConvertPPI(self,db,unmapped,rejected):  ### Uses GeneMap object and UniProt to extract gene and alias information
        '''
        Uses GeneMap object and UniProt to extract gene and alias information.
        >> db:str = Key of self.dict['PPIDB']
        >> unmapped:list = full list of unmapped (missing) proteins from all DB
        >> rejected:list = list of proteins rejected for being wrong species
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Summarise and clarify existing data storage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppidb = self.dict['PPIDB'][db]  ##> rje_biogrid or rje_hprd object
            ppi = ppidb.dict['PPI']         ##> PPI data is currently stored as a dictionary of PROTEIN ids 
            map = 'Protein'
            if db in ['HPRD']: map = 'HPRD'
            dbdata = ppidb.dict[map]        ##> Each PROTEIN has additional 'Gene' (str) and 'Alias' (list) data stored
            ppx = len(ppi)                  # No. of hubs in original dictionary
            ## ~ [0b] ~ Tools for mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            map = self.obj['GeneMap']
            ## ~ [0c] ~ New data storage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                
            new_ppi = {}                    ##> Will store new PPI data using GeneMap IDs where possible
            missing = []                    ##> Want to keep track of missing protein IDs to work out how to improve
            ## ~ [0d] ~ Debugging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#self.deBug('%s : %s' % (rje.sortKeys(ppi),db))
            rlist = [] #x# rje.randomList(ppi)[:3]
            #if db.lower() in ['biogrid','mint']: rlist = rlist[:3]

            ### ~ [1] Convert Proteins to Unified Gene Set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            px = 0.0
            for p1 in rje.sortKeys(ppi):    # Take each hub in turn, alphabetically
                self.log.printLog('\r#MAP','Converting %s %s PPI data: %.1f%%' % (rje.integerString(ppx),db,px/ppx),newline=False,log=False)
                px += 100.0
                ## ~ [1a] ~ Checking missing/rejected ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if p1 in rejected: continue
                if p1 in unmapped and p1 not in missing: missing.append(p1)
                ## ~ [1b] ~ Get gene for p1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if db in ['HPRD']: i1 = 'HPRD%s' % p1
                else: i1 = string.join(string.split(p1,':'),'')
                ga = map.getKeyID(i1,restricted=False)               # ga is now a list of gene IDs
                if not ga: ga = [self.geneMap(i1)]
                ## ~ [1c] ~ Convert spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if p1 in rlist: self.deBug('%s >> %s : %s' % (p1,ga,ppi[p1]))
                for p2 in ppi[p1]:
                    if p2 in rejected: continue
                    if p2 in unmapped and p2 not in missing: missing.append(p2)
                    if db in ['HPRD']: i2 = 'HPRD%s' % p2
                    else: i2 = string.join(string.split(p2,':'),'')
                    gb = map.getKeyID(i2,restricted=False)           # gb is now a list of gene IDs
                    if not gb: gb = [self.geneMap(i2)]
                    for g1 in ga:
                        if g1 not in new_ppi: new_ppi[g1] = []
                        for g2 in gb:
                            if g2 not in new_ppi: new_ppi[g2] = []
                            if g2 not in new_ppi[g1]: new_ppi[g1].append(g2)
                            if g1 not in new_ppi[g2]: new_ppi[g2].append(g1)
                            ## ~ [1d] ~ Add to Evidence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            if g1 not in self.dict['Evidence']: self.dict['Evidence'][g1] = {}
                            if g2 not in self.dict['Evidence']: self.dict['Evidence'][g2] = {}
                            if g2 not in self.dict['Evidence'][g1]: self.dict['Evidence'][g1][g2] = []
                            if g1 not in self.dict['Evidence'][g2]: self.dict['Evidence'][g2][g1] = []
                            try:
                                for type in ppidb.dict['Evidence'][p1][p2]:
                                    if (db,type) not in self.dict['Evidence'][g1][g2]: self.dict['Evidence'][g1][g2].append((db,type))
                                    if (db,type) not in self.dict['Evidence'][g2][g1]: self.dict['Evidence'][g2][g1].append((db,type))
                            except: pass
                            try:
                                for type in ppidb.dict['Evidence'][p2][p1]:
                                    if (db,type) not in self.dict['Evidence'][g1][g2]: self.dict['Evidence'][g1][g2].append((db,type))
                                    if (db,type) not in self.dict['Evidence'][g2][g1]: self.dict['Evidence'][g2][g1].append((db,type))
                            except: pass
                    if p1 in rlist:
                        try: self.deBug('%s >>> %s (Symbol:%s)' % (p2,g2,self.obj['GeneMap'].dict['Data'][g2]['Symbol']))
                        except: self.deBug('%s >>> %s (Symbol:-)' % (p2,g2))
                        try: self.deBug('Evidence: %s' % self.dict['Evidence'][g2])
                        except: self.deBug('No evidence for %s' % g2)
            ppidb.dict['PPI'] = new_ppi
            ppidb.dict.pop('Evidence')
            self.log.printLog('\r#PPI','Converted %s PPI data: %s prots -> %s genes.' % (db,rje.integerString(ppx),rje.integerString(len(new_ppi))))                

            ### ~ Output Missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('#MAP','%s IDs: %s missing' % (db,rje.integerString(len(missing))))
            if missing: open('%s%s.%s.missing.txt' % (self.info['ResDir'],self.info['Basefile'],db),'w').write(string.join(missing,'\n'))

        except:
            self.log.errorLog('Pingu.geneMapConvertPPI is stuffed!')
            raise
#########################################################################################################################
    def combinePPI(self):   ### Combines converted PPI objects into self.dict['PPI']
        '''Combines converted PPI objects into self.dict['PPI']'''
        try:
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            new_ppi = {}    # Convert to self.dict['PPI'] at end
            new_evi = {}    # Convert to self.dict['Evidence'] at end
            ### ~ [2] Combine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ix = 0
            for db in self.dict['PPIDB']:
                ## Add ##
                for ga in self.dict['PPIDB'][db].dict['PPI']:
                    g1 = self.geneMap(ga)
                    #x#self.deBug('%s >> %s :: %s' % (ga,g1,self.obj['GeneMap'].getKeyID(ga,restricted=False)))
                    for gb in self.dict['PPIDB'][db].dict['PPI'][ga]:
                        g2 = self.geneMap(gb)
                        ## ~ [2a] ~ Update New PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if g1 not in new_ppi: new_ppi[g1] = []
                        if g2 not in new_ppi: new_ppi[g2] = []
                        if g2 not in new_ppi[g1]:
                            new_ppi[g1].append(g2)
                            ix += 1
                        if g1 not in new_ppi[g2]:
                            new_ppi[g2].append(g1)
                            ix += 1
                        ## ~ [2b] ~ Update New Evidence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if g1 not in new_evi: new_evi[g1] = {}
                        if g2 not in new_evi: new_evi[g2] = {}
                        if g2 not in new_evi[g1]: new_evi[g1][g2] = []
                        if g1 not in new_evi[g2]: new_evi[g2][g1] = []
                        for ev in self.dict['Evidence'][ga][gb]:
                            if ev not in new_evi[g1][g2]: new_evi[g1][g2].append(ev)
                            if ev not in new_evi[g2][g1]: new_evi[g2][g1].append(ev)
                        self.deBug('%s >> %s :: %s %s' % (ga,g1,g2,new_evi[g2][g1]))
                    self.log.printLog('\r#PPI','Combining PPI: %s genes -> %s PPI' % (rje.integerString(len(new_ppi)),rje.integerString(ix)),newline=False,log=False)
            ## ~ [2c] ~ Clean and tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for g1 in new_evi:
                for g2 in new_evi[g1]:
                    self.deBug('%s-%s="%s"' % (g1,g2,new_evi[g1][g2]))
                    try:
                        new_evi[g1][g2].sort()
                        if self.dict['Evidence'] and not new_evi[g1][g2]:
                            new_ppi[g1].remove(g2); ix -= 1
                            self.progLog('\r#PPI','Filtering PPI: %s genes -> %s PPI' % (rje.integerString(len(new_ppi)),rje.integerString(ix)))
                        new_evi[g1][g2] = string.join(new_evi[g1][g2],'|')
                    except: self.errorLog('Shiiiite: %s-%s="%s"' % (g1,g2,new_evi[g1][g2]))
            for g1 in new_ppi: new_ppi[g1].sort()
            self.dict['PPI'] = new_ppi
            self.dict['Evidence'] = new_evi
            self.log.printLog('\r#PPI','PPI Combined: %s genes -> %s PPI.' % (rje.integerString(len(new_ppi)),rje.integerString(ix)))
            ### ~ [3] ~ Output interactions into files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppinum = {}     # Dictionary of Genes with given numbers of PPI 
            #!# Add additional sequence sources? #!# UniProt #!#
            ## ~ [3a] ~ Summary file for each gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chead = ['Gene','EnsLoci','PPINum','PPI']
            cfile = '%s%s.ppi.tdt' % (self.info['ResDir'],self.info['Basefile'])
            self.log.printLog('\r#PPI','Output %s PPI to %s ...' % (rje.integerString(ix),cfile),newline=False,log=False)
            rje.delimitedFileOutput(self,cfile,chead,rje_backup=True)
            for g1 in rje.sortKeys(new_ppi):
                if len(new_ppi[g1]) not in ppinum: ppinum[len(new_ppi[g1])] = []
                ppinum[len(new_ppi[g1])].append(g1)
                cdata = {'Gene':g1,'PPI':string.join(new_ppi[g1],','),'PPINum':len(new_ppi[g1])}
                try: cdata['EnsLoci'] = self.obj['GeneMap'].getGeneData(g1)['EnsLoci']
                except: cdata['EnsLoci'] = '-'
                rje.delimitedFileOutput(self,cfile,chead,datadict=cdata)
            self.printLog('\r#PPI','Output of %s PPI to %s complete.' % (rje.integerString(ix),cfile))
            nfile = '%s%s.ppinum.tdt' % (self.info['ResDir'],self.info['Basefile'])
            nhead = ['PPINum','Count']
            rje.delimitedFileOutput(self,nfile,nhead,rje_backup=True)
            for pn in rje.sortKeys(ppinum): rje.delimitedFileOutput(self,nfile,nhead,datadict={'PPINum':pn,'Count':len(ppinum[pn])})
            ## ~ [3b] ~ Pairwise output with evidence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chead = ['Hub','HubUni','Spoke','SpokePPI','SpokeUni','SpokeSeq','Evidence']
            cfile = '%s%s.pairwise.tdt' % (self.info['ResDir'],self.info['Basefile'])
            rje.delimitedFileOutput(self,cfile,chead,rje_backup=True)
            px = 0
            for g1 in rje.sortKeys(new_ppi):
                self.progLog('\r#PPI','Output of %s pairwise PPI' % rje.integerString(px))
                for g2 in new_ppi[g1]:
                    px += 1
                    cdata = {'Hub':g1,'Spoke':g2,'SpokePPI':len(new_ppi[g2])}
                    try: cdata['HubUni'] = self.obj['GeneMap'].getGeneData(g1)['UniProt']
                    except: cdata['HubUni'] = '-'
                    try: cdata['SpokeUni'] = self.obj['GeneMap'].getGeneData(g2)['UniProt']
                    except: cdata['SpokeUni'] = '-'
                    try: cdata['SpokeSeq'] = self.obj['GeneMap'].getGeneData(g2)['EnsLoci']
                    except: cdata['SpokeSeq'] = '-'
                    try: cdata['Evidence'] = self.dict['Evidence'][g1][g2]
                    except: cdata['Evidence'] = '?'
                    rje.delimitedFileOutput(self,cfile,chead,datadict=cdata)
            self.printLog('\r#PPI','Output of %s pairwise PPI to %s complete.' % (rje.integerString(px),cfile))
        except: self.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def ppiSource(self,p1,p2):  ### Returns list of source databases for PPI p1:p2
        '''Returns list of source databases for PPI p1:p2.'''
        ppisource = []
        for db in self.dict['PPIDB']:
            if p1 in self.dict['PPIDB'][db].dict['PPI'] and p2 in self.dict['PPIDB'][db].dict['PPI'][p1]: ppisource.append(db) 
        return ppisource
#########################################################################################################################
    def compDB(self):   ### Outputs a comparison of interaction datasets and their overlaps
        '''Outputs a comparison of interaction datasets and their overlaps.'''
        try:
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            geneshare = {}  # Dictionary of shared genes
            ppishare = {}   # Dictionary of shared PPI
            self.log.printLog('\r#DB','Comparing databases...',newline=False,log=False)

            ### ~ [2] Compare ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in self.dict['PPIDB']:
                ## Setup ##
                geneshare[db] = {'Data':'Genes'}
                for db2 in self.dict['PPIDB']: geneshare[db][db2] = 0
                for dx in range(len(self.dict['PPIDB'])+1): geneshare[db]['%d' % dx] = 0
                ppishare[db] = {'Data':'PPI'}   # Dictionary of shared PPI
                for db2 in self.dict['PPIDB']: ppishare[db][db2] = 0
                for dx in range(len(self.dict['PPIDB'])+1): ppishare[db]['%d' % dx] = 0
                ## DB compare ##
                for g1 in self.dict['PPIDB'][db].dict['PPI']:
                    if self.list['Genes'] and g1 not in self.list['Genes']: continue    # Only care about samples
                    ## PPI ##
                    for g2 in self.dict['PPIDB'][db].dict['PPI'][g1]:
                        dx = 0
                        for db2 in self.dict['PPIDB']:
                            if g1 in self.dict['PPIDB'][db2].dict['PPI'] and g2 in self.dict['PPIDB'][db2].dict['PPI'][g1]:
                                ppishare[db][db2] += 1
                                dx += 1
                        ppishare[db]['%d' % dx] += 1
                    ## Genes ##
                    dx = 0
                    for db2 in self.dict['PPIDB']:
                        if g1 in self.dict['PPIDB'][db2].dict['PPI']:
                            geneshare[db][db2] += 1
                            dx += 1
                    geneshare[db]['%d' % dx] += 1
            self.log.printLog('\r#DB','Comparing databases complete.')

            ### ~ [3] DB Comp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Setup ##
            cfile = '%s.dbcomp.tdt' % self.info['Basefile']
            headers = ['Database','Data','Total'] + rje.sortKeys(self.dict['PPIDB'])
            for dx in range(len(self.dict['PPIDB'])): headers.append('%d' % (dx+1))
            rje.delimitedFileOutput(self,cfile,headers,rje_backup=True)
            ## Output ##
            for db in rje.sortKeys(self.dict['PPIDB']):
                geneshare[db]['Database'] = db
                geneshare[db]['Total'] = geneshare[db][db]
                rje.delimitedFileOutput(self,cfile,headers,datadict=geneshare[db])
                ppishare[db]['Database'] = db
                ppishare[db]['Total'] = ppishare[db][db]
                rje.delimitedFileOutput(self,cfile,headers,datadict=ppishare[db])
                
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def sizeDB(self):   ### Outputs sizes of interaction datasets for histogram 
        '''Outputs sizes of interaction datasets for hisogram.'''
        try:
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppisize = {}    # Dictionary of PPI dataset sizes
            ppidb = {'Combined':self}
            for db in self.dict['PPIDB']: ppidb[db] = self.dict['PPIDB'][db]
            self.log.printLog('\r#DB','Databases size histograms...',newline=False,log=False)

            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in ppidb:
                ## Setup ##
                ppisize[db] = {}
                ## DB compare ##
                for g1 in ppidb[db].dict['PPI']:
                    dx = len(ppidb[db].dict['PPI'][g1])
                    if dx in ppisize[db]: ppisize[db][dx] += 1
                    else: ppisize[db][dx] = 1
            self.log.printLog('\r#DB','Databases size histograms complete.')

            ### ~ [3] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Setup ##
            cfile = '%s.dbsize.tdt' % self.info['Basefile']
            headers = ['DBSize','Combined'] + rje.sortKeys(self.dict['PPIDB'])
            rje.delimitedFileOutput(self,cfile,headers,rje_backup=True)
            ## Output ##
            maxdx = rje.sortKeys(ppisize['Combined'])[-1]
            dx = 0
            while dx < maxdx:
                dx += 1
                datadict = {}
                for db in headers[1:]:
                    if dx in ppisize[db]: datadict[db] = ppisize[db][dx]
                if datadict:
                    datadict['DBSize'] = dx
                    rje.delimitedFileOutput(self,cfile,headers,datadict=datadict)
                
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def loadDDI(self):  ### Reads Domain-domain interactions into self.dict['DDI']
        '''Reads Domain-domain interactions into self.dict['DDI'].'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['GeneMap'].loadPFamData()
            if not self.obj['GeneMap'].dict['PFam']: return False
            data = rje.dataDict(self,self.info['ScreenDDI'],mainkeys=['Name1'],datakeys=['Name2'],
                                headers=['Pfam1','Pfam2','Name1','Name2','Acc1','Acc2','Code1','Code2'],lists=True)
            ### ~ [2] ~ Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (dx,dtot) = (0.0,len(data))
            self.deBug(data)
            try: rje.sortKeys(data)
            except: self.errorLog('Fuck',quitchoice=True)
            for p1 in rje.sortKeys(data):
                self.progLog('\r#DDI','Parsing DDI from iPFam: %.1f%%' % (dx/dtot))
                if p1 not in self.dict['DDI']: self.dict['DDI'][p1] = []
                for p2 in data[p1]['Name2']:
                    if p2 not in self.dict['DDI']: self.dict['DDI'][p2] = []
                    if p2 not in self.dict['DDI'][p1]: self.dict['DDI'][p1].append(p2)
                    if p1 not in self.dict['DDI'][p2]: self.dict['DDI'][p2].append(p1)
            self.printLog('\r#DDI','Parsing DDI from iPFam: %s domains' % (rje.integerString(dtot)))
            ### ~ [3] ~ Convert to genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ddi = {}    # Will replace DDI dictionary
            for hub in self.dict['PPI']:
                for hseq in self.obj['GeneMap'].getEnsLoci(hub):
                    hens = hseq.shortName()
                    if hens not in self.obj['GeneMap'].dict['PFam']: continue
                    for spoke in self.dict['PPI'][hub]:
                        for sseq in self.obj['GeneMap'].getEnsLoci(spoke):
                            sens = sseq.shortName()
                            if sens not in self.obj['GeneMap'].dict['PFam']: continue
                            for hdom in self.obj['GeneMap'].dict['PFam'][hens]:
                                if hdom not in self.dict['DDI']: continue
                                for sdom in self.obj['GeneMap'].dict['PFam'][sens]:
                                    if sdom in self.dict['DDI'][hdom]:
                                        if hub not in ddi: ddi[hub] = []
                                        if spoke not in ddi[hub]: ddi[hub].append(spoke)
            self.dict['DDI'] = ddi
            self.printLog('#DDI','Converted DDI: %s genes with DDI' % rje.integerString(len(ddi)))
        except: self.errorLog('Problem with Pingu.loadDDI()',quitchoice=True)
#########################################################################################################################
    def screenDDI(self):    ### Screens out self.dict['DDI'] from self.dict['PPI']
        '''Screens out self.dict['DDI'] from self.dict['PPI'].'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(self.info['ScreenDDI']): self.dict['DDI'] = {}; return
            self.loadDDI()
            dx = 0
            ### ~ [2] ~ Screen ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in self.dict['DDI']:
                for spoke in self.dict['DDI'][hub]:
                    if spoke in self.dict['PPI'][hub]:
                        self.dict['PPI'][hub].remove(spoke)
                        dx += 1
            self.printLog('#DDI','%s domain-domain interactions screened from PPI' % rje.integerString(dx))
        except: self.errorLog('Problem with Pingu.screenDDI()')
#########################################################################################################################
    def screenPPI(self):    ### Screens out PPI based on evidence and/or complex predictions
        '''Screens out PPI based on evidence and/or complex predictions.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['NoComplex'] and not self.list['PPIType'] and not self.list['BadType']: return
            fullppi = self.dict['PPI']  # Replacement PPI dictionary
            self.dict['PPI'] = {}       # Replacement PPI dictionary
            ### ~ [2] Filter complexes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            px = 0.0; ptot = len(fullppi); ppix = 0
            for hub in fullppi:
                self.progLog('\r#PPI','Filtering PPI: %.2f%% (%s hubs; %s ppi)' % (px/ptot,rje.integerString(len(self.dict['PPI'])),rje.integerString(ppix/2))); px +=100.0
                self.dict['PPI'][hub] = []
                for spoke in fullppi[hub]:
                    try: evidence = self.dict['Evidence'][hub][spoke]
                    except: evidence = 'Unknown'
                    self.deBug('%s <-> %s = %s' % (hub,spoke,evidence))
                    ## ~ [2a] Good evidence first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    goodspoke = False
                    for ptype in self.list['PPIType']:
                        if rje.matchExp(':(%s)($|\|)' % ptype, evidence): goodspoke = True; break
                    self.deBug('%s:%s - %s = %s' % (hub,spoke,evidence,goodspoke))
                    if goodspoke: self.dict['PPI'][hub].append(spoke); continue
                    ## ~ [2b] Bad evidence filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    goodspoke = False
                    if self.list['BadType']:
                        evlist = string.split(evidence,'|')
                        for ev in evlist:
                            if string.split(ev,':')[1] not in self.list['BadType']: goodspoke = True; break
                    if not goodspoke: continue
                    ## ~ [2c] Complex filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if not self.opt['NoComplex']: goodspoke = True     # GoodPPI only
                    elif hub == spoke: continue     # Do not add self-interactions without PPIType support
                    else:
                        for spoke2 in fullppi[hub]:
                            if spoke2 in [spoke,hub]: continue      
                            if spoke2 in fullppi[spoke]: goodspoke = False; break
                    if goodspoke: self.dict['PPI'][hub].append(spoke)
                ppix += len(self.dict['PPI'][hub])
                self.deBug('%s :: %s ' % (hub,self.dict['PPI'][hub]))
                if not self.dict['PPI'][hub]: self.dict['PPI'].pop(hub)
            self.printLog('\r#PPI','Filtered PPI: (%s -> %s hubs; %s ppi)' % (rje.integerString(len(fullppi)),rje.integerString(len(self.dict['PPI'])),rje.integerString(ppix/2)))
        except: self.errorLog('Error with pingu.screenPingu()')
#########################################################################################################################
    def domPPI(self):   ### Constructs domain-based PPI data
        '''Constructs domain-based PPI data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['DomHubs'] = {}   # Dictionary of {Domain:[Hubs genes]}
            self.dict['DomPPI'] = {}    # Dictionary of {Domain:[Spoke genes]}
            mindomseq = 2               # Min number of sequences that a domain must have PPI for to be included
            data = rje.dataDict(self,self.info['PFamData'],datakeys=['Name'],lists=True)
            ### ~ [2] ~ Make Domain Hub Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(data)
            for pfam in rje.sortKeys(data):
                self.progLog('\r#DOM','Making domain-hub links: %.2f%%' % (sx/stot)); sx += 100.0
                for ens in data[pfam]['Name']:
                    hub = self.geneMap(ens)
                    if hub not in self.dict['PPI']: continue
                    if pfam not in self.dict['DomHubs']: self.dict['DomHubs'][pfam] = []
                    if hub not in self.dict['DomHubs'][pfam]: self.dict['DomHubs'][pfam].append(hub)
            self.progLog('\r#DOM','Reducing domain-hub links (%d+ hubs) ...' % (mindomseq))
            for pfam in rje.sortKeys(self.dict['DomHubs']):
                if len(self.dict['DomHubs'][pfam]) < mindomseq: self.dict['DomHubs'].pop(pfam)
            self.printLog('\r#DOM','%s domains have %d+ hubs with PPI data.' % (rje.integerString(len(self.dict['DomHubs'])),mindomseq))
            ### ~ [3] ~ Make Domain PPI Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(self.dict['DomHubs'])
            for dom in self.dict['DomHubs']:
                self.progLog('\r#DOM','Generating domain-based PPI Datasets: %.2f%%' % (sx/stot)); sx += 100.0
                if dom not in self.dict['DomPPI']: self.dict['DomPPI'][dom] = []
                for hub in self.dict['DomHubs'][dom]:
                    for spoke in self.dict['PPI'][hub]:
                        if spoke not in self.dict['DomPPI'][dom]: self.dict['DomPPI'][dom].append(spoke)
                self.deBug('%s: %s' % (dom,self.dict['DomPPI'][dom]))
            self.printLog('\r#DOM','Generated %s domain-based PPI Datasets for %s domains.' % (rje.integerString(len(self.dict['DomPPI'])),rje.integerString(len(data))))
            ### ~ [3] ~ Output Domain PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.domainPPISeqFiles()
        except: self.errorLog('Problem with Pingu.domPPI()',quitchoice=True)
#########################################################################################################################
    def OLDdomPPI(self):   ### Constructs domain-based PPI data
        '''Constructs domain-based PPI data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['DomPPI'] = {}    # Dictionary of {Domain:[Spoke genes]}
            if not self.obj['GeneMap'].dict['PFam']:
                if not rje.exists(self.obj['GeneMap'].info['PFamData']) and rje.exists(self.info['PFamData']):
                    self.obj['GeneMap'].info['PFamData'] = self.info['PFamData']
                self.obj['GeneMap'].loadPFamData()
            self.printLog('#PFAM','PFam data loaded for %s sequences' % rje.integerString(len(self.obj['GeneMap'].dict['PFam'])))
            if not self.obj['GeneMap'].dict['PFam']: return False
            ### ~ [2] ~ Make Domain PPI Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(self.obj['GeneMap'].dict['PFam'])
            for seq in self.obj['GeneMap'].dict['PFam']:
                self.progLog('\r#DOM','Generating domain-based PPI Datasets: %.2f%%' % (sx/stot)); sx += 100.0
                hub = self.geneMap(seq)
                self.deBug('%s: %s' % (hub,hub in self.dict['PPI']))
                if hub not in self.dict['PPI']: continue
                for dom in self.obj['GeneMap'].dict['PFam'][seq]:
                    if dom not in self.dict['DomPPI']: self.dict['DomPPI'][dom] = []
                    for spoke in self.dict['PPI'][hub]:
                        if spoke not in self.dict['DomPPI'][dom]: self.dict['DomPPI'][dom].append(spoke)
                self.deBug('%s: %s' % (dom,self.dict['DomPPI'][dom]))
            self.printLog('\r#DOM','Generated %s domain-based PPI Datasets for %s domains.' % (rje.integerString(len(self.dict['DomPPI'])),rje.integerString(stot)))
            ### ~ [3] ~ Output Domain PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.domainPPISeqFiles()
        except: self.errorLog('Problem with Pingu.domPPI()',quitchoice=True)
#########################################################################################################################
    def addLinks(self,iterate=False): ### Add proteins linking Sample proteins through PPI
        '''Add proteins linking Sample proteins through PPI.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Links' not in self.dict['Datasets']:
                self.dict['Datasets']['Links'] = []
            allgenes = []
            for sample in self.dict['Datasets']:
                if sample == 'Links':
                    if not iterate: continue
                elif self.list['Experiments'] and sample not in self.list['Experiments']: continue
                for gene in self.dict['Datasets'][sample]:
                    if gene not in allgenes: allgenes.append(gene)
            ## ~ [0a] ~ Setup Orphan list that needs to be joined ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.opt['LinkOrphans'] = True
            if self.opt['LinkOrphans']:
                orphans = []
                for gene in allgenes:
                    if gene not in self.dict['PPI']: continue
                    orphans.append(gene)
                    for ppi in self.dict['PPI'][gene]:
                        if ppi in allgenes: orphans.remove(gene); break
            else: orphans = allgenes[0:]
            ### ~ [1] ~ Add Links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            added = True; cycx = 0
            while added:
                added = False; cycx +=1; ax = 0.0; atot = len(allgenes)
                for gene in allgenes[0:]:
                    self.progLog('\r#ADD','Add PPI Links (Cycle %d): %.2f%%' % (cycx,ax/atot)); ax += 100.0
                    if gene not in self.dict['PPI']: continue
                    for ppi in self.dict['PPI'][gene]:
                        if ppi in allgenes + self.dict['Datasets']['Links']: continue
                        if self.stat['RemSticky'] > 0 and len(self.dict['PPI'][ppi]) > self.stat['RemSticky']: continue
                        for ppi2 in self.dict['PPI'][ppi]:
                            if ppi2 == gene: continue
                            if ppi2 in allgenes:
                                self.dict['Datasets']['Links'].append(ppi)
                                self.deBug('Adding %s-> %s <-%s' % (gene,ppi,ppi2))
                                if iterate: allgenes.append(ppi)
                                if iterate and not self.opt['LinkOrphans']: orphans.append(ppi)
                                added = True
                                break
                self.printLog('\r#ADD','Added PPI Links (Cycle %d): %d links' % (cycx,len(self.dict['Datasets']['Links'])))
            self.dict['Samples']['Links'] = self.dict['Datasets']['Links']
            #self.deBug(self.dict['Datasets'])
        except: self.errorLog('Problem with Pingu.addLinks()')
#########################################################################################################################
    ### <7> ### All-by-All analysis Methods                                                                             #
#########################################################################################################################
    def gablamMatrix(self):   ### Uses GABLAM to generate relationships based on sequence identity.                 |3.0|
        '''Uses GABLAM to generate relationships based on sequence identity.'''
        try:### ~ [1] Setup objects etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = '%s.Combined.fas' % self.info['Basefile']
            if not os.path.exists(seqfile): return False
            self.obj['GABLAM'] = famgablam = gablam.GABLAM(self.log,self.cmd_list)
            famgablam.opt['QryAcc'] = False
            gfile = famgablam.info['GablamOut'] = '%s.gablam.tdt' % self.info['Basefile']
            famgablam.info['HitSumOut'] = '%s.hitsum.tdt' % self.info['Basefile']
            famgablam.info['LocalOut'] = '%s.local.tdt' % self.info['Basefile']
            famgablam.info['QueryDB'] = famgablam.info['SearchDB'] = seqfile
            #gablamcut=X

            ### ~ [2] Use GABLAM to generate relationships ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(gfile) or (self.interactive() >= 2 and not rje.yesNo('Use existing %s?' % gfile)):
                famgablam.gablam()
            ensloci = self.getEnsLoci()
            seqdict = self.dict['EnsSeqNameDict']
            ensname = {}    # Mapping of sequence short names to gene names
            for seq in seqdict.values():
                try: gene = rje.matchExp(' gene:(\S+)\]',seq.info['Name'])[0]
                except:
                    self.log.errorLog('Problem extracting gene from "%s"' % seq.info['Name']) 
                    gene = seq.shortName()
                symbol = self.geneMap(gene)
                ensname[seq.shortName()] = symbol
            matrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            matrix.loadFromDataTable(gfile,normalise=100.0,inverse=True,checksym=False,distance='Qry_%s' % famgablam.info['CutStat'])  ### Loads matrix from GABLAM
            matrix.forceSymmetry()          ### Use min distance for each pair of values
            matrix.rename(ensname)          ### Goes through matrix and renames objects using given dictionary
            self.obj['Homology'] = matrix   ### Should now be mapped 
            matrix.saveCytoscape(self.info['Basefile'],type='gablam',cutoff=1.0,inverse=False)  ### Output for Cytoscape
            return True
        except:
            self.log.errorLog('PINGU.gablamMatrix() failure')
            raise
        return False
#########################################################################################################################        
    def purgeControls(self):    ### Removes genes found in Controls from Experiments
        '''Removes genes found in Controls from Experiments.'''
        try:### ~ [1] Setup list of Control genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            controls = []
            for c in self.list['Controls']:
                if c in self.dict['Datasets']: controls += self.dict['Datasets'][c]
            self.log.printLog('#CTRL','%d control genes to remove from %d experiments' % (len(controls),len(self.list['Experiments'])))
            if not controls or not self.list['Experiments']: return 

            ### ~ [2] Cut down Experiments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            experiments = self.list['Experiments']
            if string.join(experiments).lower() == 'all' and experiments[0] not in self.dict['Samples']:
                experiments = rje.sortKeys(self.dict['Samples'])
            for sample in self.list['Experiments']:
                if sample not in self.dict['Datasets']:
                    self.log.printLog('#EXP','No experimental sample "%s"' % sample)
                    continue
                ok = []
                for gene in self.dict['Datasets'][sample]:
                    if gene not in controls: ok.append(gene)
                self.log.printLog('#CTRL','%d %s genes -> %d not in controls' % (len(self.dict['Datasets'][sample]),sample,len(ok)))
                self.dict['Datasets'][sample] = ok[0:]
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def allByAll(self): ### Generates All-by-all outputs
        '''Generates All-by-all outputs.'''
        try:
            ### ~ [1] Setup data structures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            allbyall = {}
            gx = 0.0

            ### ~ [2] Generate data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for g1 in self.list['Genes']:
                self.log.printLog('\r#ALL','All-by-all PPI matrix: %.1f%%' % (gx/len(self.list['Genes'])),newline=False,log=False)
                gx += 100.0
                ## Setup ##
                allbyall[g1] = {'Gene':g1}
                prev = []
                this = []
                next = []
                dx = 1      # Degrees of separation
                if g1 in self.dict['PPI']: next = self.dict['PPI'][g1][0:]
                ## Expand ##
                while dx <= self.stat['AllByAll'] and next:
                    ## This degree of separation ##
                    this = next[0:]
                    for g2 in self.list['Genes']:
                        if g2 not in allbyall[g1] and g2 in this: allbyall[g1][g2] = dx
                    ## Next degree of separation ##
                    if dx == self.stat['AllByAll']: break
                    dx += 1
                    next = []
                    for g2 in this:
                        if g2 not in prev: prev.append(g2)
                        for g3 in rje.dictValues(self.dict['PPI'],g2):
                            if g3 not in prev + this + next: next.append(g3)
            self.log.printLog('\r#ALL','All-by-all PPI matrix output...',newline=False,log=False)

            ### ~ [3] Output All-by-All data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                            
            ## Setup ##
            cfile = '%s.allbyall.tdt' % self.info['Basefile']
            headers = ['Gene'] + self.list['Genes']
            rje.delimitedFileOutput(self,cfile,headers,rje_backup=True)
            ## Output ##
            for gene in headers[1:]: rje.delimitedFileOutput(self,cfile,headers,datadict=allbyall[gene])
            self.log.printLog('\r#ALL','All-by-all PPI matrix output to %s' % cfile)
            
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def pathFinder(self):   ### Establishes and outputs links between proteins in Samples
        '''Establishes and outputs links between proteins in Samples.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pathfinder = []     # List of pathfinder genes
            linked = {}         # Dictionary of linked proteins {gene:{linked:[path]}}
            targets = {}        # Dictionary of current targets for each protein to "find"
            deadends = {}       # Dictionary of {protein:dead end links}
            ## ~ [1b] ~ Output file headers and any read in data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pfile = '%s.pathfinder.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,pfile,['Gene1','Gene2','Links','Path'],rje_backup=True)
            ## ~ [1b] ~ Setup dictionaries for queries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.list['PathQry']: self.list['PathQry'] = self.list['Genes'][0:]
            pathqry = self.list['PathQry'][0:]
            alltargets = rje.sortUnique(self.list['PathQry'] + self.list['Genes'],False)
            for gene in self.list['PathQry']:
                if gene not in self.list['Genes']: pathfinder.append(gene)
                linked[gene] = {}
                targets[gene] = alltargets[0:]
                targets[gene].remove(gene)
                deadends[gene] = [gene]       
                if gene in self.dict['PPI'] and gene in self.dict['PPI'][gene]:
                    open(pfile,'a').write('%s\n' % string.join([gene,gene,'1',self.getEvidence(gene,gene)],'\t'))
                elif gene in self.dict['PPI']:
                    open(pfile,'a').write('%s\n' % string.join([gene,gene,'0','Self'],'\t'))
            ## ~ [1c] ~ Account for no-PPI genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ax = len(alltargets)
            for gene in alltargets[0:]:
                if gene not in self.dict['PPI']:
                    open(pfile,'a').write('%s\n' % string.join([gene,'N/A','N/A','NoPPI'],'\t'))
                    alltargets.remove(gene)
                    if gene in pathqry: pathqry.remove(gene)
                    if gene in targets: targets.pop(gene)
            self.printLog('#NOPPI','%d genes removed from PathFinder: No PPI. %d remain' % (ax-len(alltargets),len(alltargets)))
            ppgenes = alltargets[0:]    #self.list['Genes'][0:]
            for gene in alltargets:      #self.list['Genes']:
                if gene not in self.dict['PPI']:
                    ppgenes.remove(gene)
                    if gene in pathqry:
                        pathqry.remove(gene)
                        for target in targets.pop(gene):
                            if target not in self.dict['PPI']: ev = 'NoHubPPI|NoSpokePPI'
                            else: ev = 'NoHubPPI'
                            open(pfile,'a').write('%s\n' % string.join([gene,target,'-1',ev],'\t'))
                    for qry in pathqry:
                        try: targets[qry].remove(gene)
                        except: self.errorLog('%s from %s - %s' % (gene,qry,targets))
                        if qry not in self.dict['PPI']: ev = 'NoHubPPI|NoSpokePPI'
                        else: ev = 'NoSpokePPI'
                        open(pfile,'a').write('%s\n' % string.join([qry,gene,'-1',ev],'\t'))
            ## ~ [1d] ~ Create first-order links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gene in pathqry:   #ppgenes:
                for ppi in self.dict['PPI'][gene]:
                    linked[gene][ppi] = [ppi]
                    if gene in targets and ppi in targets[gene]:
                        path = '%s::%s' % (ppi,self.getEvidence(gene,ppi))
                        open(pfile,'a').write('%s\n' % string.join([gene,ppi,'1',path],'\t'))
                        targets[gene].remove(ppi)
                        if not targets[gene]: targets.pop(gene); linked.pop(gene)

            ### ~ [2] ~ Establish paths and output as found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdebug = self.opt['DeBug']
            px = 1              # Length of current paths under investigation
            while targets:      # Still need to find targets
                self.opt['DeBug'] = pdebug 
                px += 1
                if px > self.stat['PathFinder'] and self.stat['PathFinder'] > 0: break
                (tx, ttot) = (-100.0, len(targets))
                self.progLog('\r#PATH','Pathfinder %d gene links: (%d)           ' % (px, ttot)); 
                for gene in rje.sortKeys(targets):
                    self.deBug('%s: %d links; %d links' % (gene,px,len(linked[gene])))
                    tx += 100.0; self.progLog('\r#PATH','Pathfinder %d gene links: %.2f%%         ' % (px, tx/ttot))
                    (lx,lj) = (0.0, (100.0/(len(linked[gene])*ttot)))
                    for linker in rje.sortKeys(linked[gene]):
                        if gene not in linked: break
                        if self.opt['DeBug'] or self.stat['Verbose'] > 0:
                            self.progLog('\r#PATH','Pathfinder %d gene links: %.3f%% (%s|%d|%d)   ' % (px, (tx+lx)/ttot,gene,len(linked[gene]),len(deadends[gene]))); lx += lj
                        linkpath = linked[gene].pop(linker)
                        deadends[gene].append(linker)
                        ## ~ [2a] ~ Expand into PPIs of linker ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                        
                        for ppi in self.dict['PPI'][linker]:
                            if ppi not in linked[gene].keys() + deadends[gene]:  # Not looked at before
                                linked[gene][ppi] = linkpath + [ppi]
                                if ppi in targets[gene]:
                                    path = [gene] + linked[gene][ppi]
                                    for i in range(len(path),0,-1)[1:]:
                                        if path[i] not in pathfinder + self.list['Genes']: pathfinder.append(path[i])
                                        path[i] = '%s::%s' % (path[i],self.getEvidence(path[i-1],path[i]))
                                    path = string.join(path[1:],'|>>|')
                                    open(pfile,'a').write('%s\n' % string.join([gene,ppi,'%d' % px,path],'\t'))
                                    targets[gene].remove(ppi)
                                    if not targets[gene]: targets.pop(gene); linked.pop(gene); break
                    if gene in linked and not linked[gene]:    # Run out of places to go
                        for target in targets[gene]: open(pfile,'a').write('%s\n' % string.join([gene,target,'-1','NoLinkage[%d]' % px],'\t'))
                        targets.pop(gene); linked.pop(gene); deadends.pop(gene)
            self.printLog('\r#PATH','Pathfinder gene links complete: <=%d linking genes.' % (px))
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
        try:### ~ [3] ~ Add data as it is, even if killed? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pathfinder.sort()
            self.dict['Datasets']['PathFinder'] = self.dict['Samples']['PathFinder'] = pathfinder
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def sampleOverlap(self):    ### Produces a table of overlap (mapped through HGNC) between samples (and 1y PPI)  |3.0|
        '''
        Produces a table of the overlap (mapped through HGNC) between samples (and 1y PPI?). Each sample and PPI database
        is considered for the genes it contains and compared to each sample and each the 1y interactors of each given
        bait for self.list['Baits']. In addition, 'Combined' sample and 'PPI' sets are also compared and an 'AllBaits'
        comparison is made.
        '''
        try:
            ### ~ [1] Setup file headers and data structures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['Sample','Total'] + rje.sortKeys(self.dict['Samples']) + ['Combined','PPI'] + rje.sortKeys(self.dict['PPIDB'])
            cfile = '%s.overlap.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,cfile,headers,rje_backup=True)
            ### ~ [2] Setup data Bait 1ry PPI GeneLists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            addbaits = self.list['Baits']
            if self.opt['AddBaits']: addbaits = []
            baitlist = {'AllBaits':[]}
            for bait in addbaits:
                baitlist[bait] = []
                if bait not in self.dict['PPI']: continue
                for gene in self.dict['PPI'][bait]:
                    if gene not in baitlist[bait]: baitlist[bait].append(gene)
                    if gene not in baitlist['AllBaits']: baitlist['AllBaits'].append(gene)
            if addbaits: addbaits.append('AllBaits')           
            ### ~ [3] Sample comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (vx,vxx) = (0.0,(len(addbaits) + len(self.dict['Samples']) + 1) * len(headers[2:]))
            for sample in rje.sortKeys(self.dict['Samples']) + ['Combined'] + addbaits:
                ## Get genelist 1 ##
                if sample in self.dict['Samples']: genes = self.dict['Datasets'][sample]
                elif sample == 'Combined': genes = self.identifiers()
                else: genes = baitlist[sample]
                overlap = {'Sample':sample,'Total':len(genes)}
                ## Compare ##
                for comp in headers[2:]:
                    self.progLog('\r#VS','Sample vs Sample/Bait comparisons: %.1f%%' % (vx/vxx)); vx += 100.0
                    overlap[comp] = 0
                    ## Get genelist 2 ##
                    if comp in self.dict['Samples']: cgenes = self.dict['Datasets'][comp]
                    elif comp == 'Combined': cgenes = self.identifiers()
                    elif comp == 'PPI': cgenes = self.dict['PPI'].keys()
                    else: cgenes = self.dict['PPIDB'][comp].dict['PPI'].keys()
                    ## Calculate ##
                    for g1 in genes:
                        if g1 in cgenes: overlap[comp] += 1
                ## Output ##
                rje.delimitedFileOutput(self,cfile,headers,datadict=overlap)
            self.printLog('\r#VS','Sample vs Sample/Bait comparisons output to %s' % cfile)
            self.clusterPPI()
            ### ~ [4] ~ Add Venn Diagram Figure output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vfile = '%s.venn.tdt' % self.info['Basefile']
            vhead = ['Grouping','Count'] + rje.sortKeys(self.dict['Samples'])
            rje.delimitedFileOutput(self,vfile,vhead,rje_backup=True)
            vdata = {}
            for gene in self.identifiers():
                grp = []
                for sample in self.samples():
                    if gene in self.dict['Datasets'][sample]: grp.append(sample)
                group = string.join(grp,'::')
                if group not in vdata:
                    vdata[group] = {'Grouping':group,'Count':0}
                    for sample in grp: vdata[group][sample] = 0
                if not grp: self.deBug('No Group for %s' % gene)
                vdata[group]['Count'] +=1 
                for sample in grp: vdata[group][sample] += 1
            for group in rje.sortKeys(vdata): rje.delimitedFileOutput(self,vfile,vhead,datadict=vdata[group])
            self.printLog('#VENN','Data for Venn diagrams output to %s' % vfile)
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def clusterPPI(self):    ### Produces a table of the 1ry PPI distance for each sample/bait and gene.
        '''Produces a table of the 1ry PPI distance for each sample/bait and gene.'''
        try:
            ### ~ [1] Setup file headers and data structures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] List of row/column headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cgenes = []     # List of genes for comparison
            for sample in rje.sortKeys(self.dict['Samples']): cgenes.append('%s*' % sample)
            cgenes += self.list['Baits']
            for gene in self.list['Genes']:
                if gene not in cgenes: cgenes.append(gene)
            ## ~ [1b] Output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = ['Identifier'] + cgenes
            cfile = '%s.cluster.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,cfile,headers,rje_backup=True)
            ## ~ [1c] Distance dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdist = {}
            for g in cgenes: cdist[g] = {}

            ### ~ [2] Calculate distances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (vx,vxx) = (0.0,len(cgenes))
            for g1 in cgenes:
                self.progLog('#DIS','Sample vs Sample/Bait/Gene PPI distances: %.1f%%' % (vx/vxx))
                vx += 100.0
                if g1[-1] == '*': pp1 = self.dict['Datasets'][g1[:-1]]
                elif g1 in self.dict['PPI']: pp1 = self.dict['PPI'][g1]
                else: pp1 = []
                for g2 in cgenes:
                    if g2[-1] == '*': pp2 = self.dict['Datasets'][g2[:-1]]
                    elif g2 in self.dict['PPI']: pp2 = self.dict['PPI'][g2]
                    else: pp2 = []
                    if not pp1 or not pp2:      # No interactions -> None shared!
                        cdist[g1][g2] = 1.0
                        continue
                    ppx = 0.0
                    for pp in pp2:
                        if pp in pp1: ppx += 1.0
                    cdist[g1][g2] = ((1.0 - (ppx / len(pp1))) + (1.0 - (ppx / len(pp2)))) / 2.0     # Mean shared proportion
                if len(cdist[g1]) != len(cgenes): self.errorLog('Problem with length of cluster dict for %s' % g1,printerror=False)
                cdist[g1]['Identifier'] = g1
                rje.delimitedFileOutput(self,cfile,headers,datadict=cdist[g1])
            self.log.printLog('\r#DIS','Sample vs Sample/Bait/Gene PPI distances output to %s' % cfile)
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def sampleXGMML(self):      ### Outputs XGMML format file for Cytoscape
        '''Outputs XGMML format file for Cytoscape.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['XGFormat']: return self.simpleSampleXGMML()
            self.getEnsLoci()
            seqdict = self.dict['EnsSeqNameDict']
            xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            experiments = self.list['Experiments']
            if string.join(experiments).lower() == 'all' and experiments[0] not in self.dict['Samples']:
                experiments = rje.sortKeys(self.dict['Samples'])
            ## ~ [1a] Node attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodeatt = xgmml.dict['NodeAtt'] = {'SAMPLE':'string','CONTROL':'string','HUBSIZE':'real',
                                               'DESCRIPTION':'string','ENSLOCI':'string','ENSGENE':'string'}
            for cyt in ['shape','fillColor','size','fontSize','borderColor','font']: xgmml.dict['NodeAtt']['node.%s' % cyt] = 'string' 
            edgeatt = xgmml.dict['EdgeAtt'] = {'GABLAM':'real'}
            for cyt in ['lineStyle','color']: xgmml.dict['EdgeAtt']['edge.%s' % cyt] = 'string' 
            ## ~ [1b] Setup list of control & experiment genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cgenes = []
            for sample in self.list['Controls']:
                if sample in self.dict['Datasets']:
                    cgenes += self.dict['Datasets'][sample]
            egenes = []
            for sample in experiments:
                if sample in self.dict['Datasets']:
                    egenes += self.dict['Datasets'][sample]
                
            ### ~ [2] Output XGMML files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Datasets']['Combined'] = rje.sortUnique(egenes)
            if 'PathFinder' in self.dict['Datasets']: self.dict['Datasets']['Combined'] += self.dict['Datasets']['PathFinder']
            if 'Links' in self.dict['Datasets']: self.dict['Datasets']['Combined'] = rje.sortUnique(self.dict['Datasets']['Combined']+self.dict['Datasets']['Links'])
            for sample in rje.sortKeys(self.dict['Samples']) + ['Combined']:
                ## ~ [2a] Setup sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if sample in self.dict['Samples']: sgenes = self.dict['Datasets'][sample][0:]
                elif sample == 'Combined': sgenes = self.list['Genes'][0:]
                sgenes.sort()
                xgmml.info['Name'] = '%s:%s' % (rje.baseFile(self.info['Name']),sample)
                xfile = '%s.%s.xgmml' % (self.info['Basefile'],sample)
                ## ~ [2b] Update node data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nodes = xgmml.dict['Node'] = {}
                for gene in sgenes:
                    desc = hgnc = gene
                    data = {}
                    rje.combineDict(data,self.obj['GeneMap'].getGeneData(hgnc),overwrite=False,replaceblanks=True)
                    #x#else: self.log.printLog('#CARD','No GeneCard data for %s' % hgnc)
                    #x#for ens in self.dict['EnsGenes'][hgnc]: rje.combineDict(data,self.obj['GeneMap'].getGeneData(ens),overwrite=False,replaceblanks=True)
                    if rje.dictValues(data,'EnsDesc',valtype='str'): data['Desc'] = data['EnsDesc']
                    if rje.dictValues(data,'Desc',valtype='str'): desc = data['Desc']
                    try: hubx = len(self.dict['PPI'][gene])
                    except: hubx = 0
                    nodes[gene] = {'HUBSIZE':hubx,'DESCRIPTION':desc}
                    try: nodes[gene]['ENSLOCI'] = data['EnsLoci']
                    except: self.log.printLog('#ENS','No EnsLoci data for %s' % gene)
                    try: nodes[gene]['ENSGENE'] = data['EnsEMBL']
                    except: self.log.printLog('#ENS','No EnsEMBL gene data for %s' % gene)
                    if gene in cgenes:
                        nodes[gene]['CONTROL'] = 'Control'
                        if gene in egenes: nodes[gene]['node.shape'] = 'roundrect'
                        else: nodes[gene]['node.shape'] = 'rectangle'
                    elif gene not in egenes:
                        nodes[gene]['CONTROL'] = 'Extra'
                        nodes[gene]['node.shape'] = 'diamond'
                    else:
                        nodes[gene]['CONTROL'] = 'Experiment'
                        nodes[gene]['node.shape'] = 'ellipse'
                    nodes[gene]['node.font'] = 'Impact,plain,10'
                    nodes[gene]['node.fontSize'] = '10.0'
                    nsize = 30
                    if hubx >= 100: nsize = 60
                    elif hubx >= 50: nsize = 50
                    elif hubx >= 10: nsize = 40
                    elif hubx >= 0: nsize = 35
                    nodes[gene]['node.size'] = '%d.0' % nsize
                ## ~ [2b-i] Sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for gene in sgenes[:-1]:
                    samples = []
                    for sample in rje.sortKeys(self.dict['Samples']):
                        if gene in self.dict['Datasets'][sample]: samples.append(sample)
                    ## Simple case of only in one sample ##
                    if len(samples) == 1:
                        nodes[gene]['SAMPLE'] = samples[0]
                        if samples[0] in self.list['Controls']: nodes[gene]['node.fillColor'] = '153,153,153'
                        elif samples[0] in experiments: nodes[gene]['node.fillColor'] = '153,153,255'
                        else: nodes[gene]['node.fillColor'] = '153,255,153'
                        continue
                    ## Remove controls ##
                    for c in self.list['Controls']:
                        if c in samples: samples.remove(c)
                    if len(samples) == 1:
                        nodes[gene]['SAMPLE'] = samples[0]
                        nodes[gene]['node.fillColor'] = '255,153,255'
                        continue
                    ## Remove non-experiments, else return "extras"                    
                    for s in samples[0:]:
                        if s not in experiments: samples.remove(s)
                    nodes[gene]['node.fillColor'] = '255,153,153'
                    if not samples:
                        nodes[gene]['SAMPLE'] = 'extras'
                        nodes[gene]['node.fillColor'] = '153,255,153'
                    elif not self.opt['CompressPP'] or len(samples) == 1: nodes[gene]['SAMPLE'] = string.join(samples,'-')
                    else: nodes[gene]['SAMPLE'] = 'Share%d' % len(samples)
                ## ~ [2b-ii] Additional Cytoscape details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Colours: Pale yellow = 255,255,153; pale blue = 153,153,255; green = 0,255,0; [RGB]
                #<att type="string" name="node.fillColor" label="node.fillColor" value="153,153,255"/>
                #<att type="string" name="node.size" label="node.size" value="50.0"/>
                #<att type="string" name="node.borderColor" label="node.borderColor" value="0,0,255"/>
                    
                ## ~ [2c] Update edge data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                edges = xgmml.dict['Edge'] = {'pp':{},'hom':{}}
                for g1 in sgenes[:-1]:
                    ## ~ [2c-i] PPI edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for g2 in sgenes[sgenes.index(g1)+1:]:
                        if g1 in self.dict['PPI'] and g2 in self.dict['PPI'][g1]:
                            edges['pp'][(g1,g2)] = {}       #!# Add attributes with time? (PPIDB etc?)
                            edges['pp'][(g1,g2)]['edge.lineStyle'] = 'SOLID'
                            edges['pp'][(g1,g2)]['edge.color'] = '0,0,0'
                    ## ~ [2c-ii] Homology ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if self.opt['GABLAM']:
                            dis = self.obj['Homology'].getDis(g1,g2)
                            if dis:
                                dis = 100.0 * (1.0 - dis)
                                if dis >= self.obj['GABLAM'].stat['GABLAMCut']:
                                    #edges['hom'][(g2,g1)] = {'GABLAM':dis,'edge.lineStyle':'SOLID','edge.color':'255,0,0'}
                                    edges['hom'][(g2,g1)] = {'GABLAM':dis,'edge.lineStyle':'LONG_DASH','edge.color':'0,0,255'}
                #<att type="string" name="edge.targetArrowShape" label="edge.targetArrowShape" value="DELTA"/>
                #<att type="string" name="edge.targetArrowColor" label="edge.targetArrowColor" value="0,0,255"/>

            
                ## ~ [2d] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                xgmml.saveXGMML(xfile); return xgmml
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def xgGeneSample(self,gene,default='extras'):    ### Returns sample(s) for gene for 'SAMPLE' attribute
        '''Returns sample(s) for gene for 'SAMPLE' attribute.'''
        samples = []
        for sample in rje.sortKeys(self.dict['Samples']):
            if gene in self.dict['Datasets'][sample]: samples.append(sample)
        ## Simple case of only in one sample ##
        if len(samples) == 1: return samples[0]
        ## Remove controls ##
        for c in self.list['Controls']:
            if c in samples: samples.remove(c)
        #if len(samples) == 1:
        #    nodes[gene]['SAMPLE'] = samples[0]
        #    continue
        ## Remove non-experiments, else return "extras"
        if not self.opt['AddLinks']:
            for s in samples[0:]:
                if s not in experiments: samples.remove(s)
        if not samples: return default
        elif not self.opt['CompressPP'] or len(samples) == 1: return string.join(samples,'-')
        else: return 'Share%d' % len(samples)
#########################################################################################################################
    def simpleSampleXGMML(self):    ### Outputs XGMML format file for Cytoscape with formatting (use Styles)
        '''Outputs XGMML format file for Cytoscape with formatting (use Styles).'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.getEnsLoci()
            seqdict = self.dict['EnsSeqNameDict']
            xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            experiments = self.list['Experiments']
            if string.join(experiments).lower() == 'all' and experiments[0] not in self.dict['Samples']:
                experiments = rje.sortKeys(self.dict['Samples'])
            if 'Links' in self.dict['Datasets'] and 'Links' not in experiments: experiments.append('Links') 
            ## ~ [1a] Node attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodeatt = xgmml.dict['NodeAtt'] = {'SAMPLE':'string','CONTROL':'string','HUBSIZE':'real',
                                               'DESCRIPTION':'string','ENSLOCI':'string','ENSGENE':'string'}
            edgeatt = xgmml.dict['EdgeAtt'] = {'GABLAM':'real','EVIDENCE':'string'}
            ## ~ [1b] Setup list of control & experiment genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cgenes = []
            for sample in self.list['Controls']:
                if sample in self.dict['Datasets']:
                    cgenes += self.dict['Datasets'][sample]
            egenes = []
            for sample in experiments:
                if sample in self.dict['Datasets']:
                    egenes += self.dict['Datasets'][sample]
                
            ### ~ [2] Output XGMML files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Datasets']['Combined'] = rje.sortUnique(egenes)
            if 'PathFinder' in self.dict['Datasets']: self.dict['Datasets']['Combined'] += self.dict['Datasets']['PathFinder']
            if 'Links' in self.dict['Datasets']: self.dict['Datasets']['Combined'] = rje.sortUnique(self.dict['Datasets']['Combined']+self.dict['Datasets']['Links'])
            for sample in rje.sortKeys(self.dict['Samples']) + ['Combined']:
                ## ~ [2a] Setup sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #if sample in self.dict['Samples']: sgenes = self.dict['Datasets'][sample][0:]
                #elif sample == 'Combined': sgenes = self.list['Genes'][0:]
                sgenes = self.dict['Datasets'][sample][0:]
                sgenes.sort()
                xgmml.info['Name'] = '%s:%s' % (rje.baseFile(self.info['Name']),sample)
                xfile = '%s.%s.xgmml' % (self.info['Basefile'],sample)
                ## ~ [2b] Update node data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nodes = xgmml.dict['Node'] = {}
                for gene in sgenes:
                    desc = hgnc = gene
                    data = {}
                    rje.combineDict(data,self.obj['GeneMap'].getGeneData(hgnc),overwrite=False,replaceblanks=True)
                    #x#else: self.log.printLog('#CARD','No GeneCard data for %s' % hgnc)
                    #x#if hgnc in self.dict['EnsGenes']:
                    #x# for ens in self.dict['EnsGenes'][hgnc]: rje.combineDict(data,self.obj['GeneMap'].getGeneData(ens),overwrite=False,replaceblanks=True)
                    if rje.dictValues(data,'EnsDesc',valtype='str'): data['Desc'] = data['EnsDesc']
                    if rje.dictValues(data,'Desc',valtype='str'): desc = data['Desc']
                    try: hubx = len(self.dict['PPI'][gene])
                    except: hubx = 0
                    nodes[gene] = {'HUBSIZE':hubx,'DESCRIPTION':desc}
                    try: nodes[gene]['ENSLOCI'] = data['EnsLoci']
                    except: self.log.printLog('#ENS','No EnsLoci data for %s' % gene)
                    try: nodes[gene]['ENSGENE'] = data['EnsEMBL']
                    except: self.log.printLog('#ENS','No EnsEMBL gene data for %s' % gene)
                    if gene in cgenes: nodes[gene]['CONTROL'] = 'Control'
                    elif gene not in egenes: nodes[gene]['CONTROL'] = 'Extra'
                    else: nodes[gene]['CONTROL'] = 'Experiment'
                ## ~ [2b-i] Sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for gene in sgenes:
                    nodes[gene]['SAMPLE'] = self.xgGeneSample(gene); continue
                    samples = []
                    for sample in rje.sortKeys(self.dict['Samples']):
                        if gene in self.dict['Datasets'][sample]: samples.append(sample)
                    ## Simple case of only in one sample ##
                    if len(samples) == 1:
                        nodes[gene]['SAMPLE'] = samples[0]
                        continue
                    ## Remove controls ##
                    for c in self.list['Controls']:
                        if c in samples: samples.remove(c)
                    #if len(samples) == 1:
                    #    nodes[gene]['SAMPLE'] = samples[0]
                    #    continue
                    ## Remove non-experiments, else return "extras"
                    if not self.opt['AddLinks']:
                        for s in samples[0:]:
                            if s not in experiments: samples.remove(s)
                    if not samples:
                        nodes[gene]['SAMPLE'] = 'extras'
                    elif not self.opt['CompressPP'] or len(samples) == 1: nodes[gene]['SAMPLE'] = string.join(samples,'-')
                    else: nodes[gene]['SAMPLE'] = 'Share%d' % len(samples)
                ## ~ [2b-ii] Add Seconday/Tertiary etc. interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                next = xnodes = sgenes[0:]
                for ex in range(self.stat['XGExpand']):
                    prev = next[0:]
                    next = []
                    gx = 0.0; gtot = len(prev)
                    for g1 in prev[0:]:
                        self.progLog('\r#PPI','Expanding PPI network (Lvl %d): %.2f%%' % (ex+1,gx/gtot)); gx += 100.0
                        if g1 not in self.dict['PPI']: continue
                        for g2 in self.dict['PPI'][g1]:
                            if self.opt['XGComplex'] and self.getEvidence(g1,g2).lower().find('complex') < 1: continue
                            if g2 not in xnodes:     # Add this gene
                                next.append(g2)
                                xnodes.append(g2)
                                desc = hgnc = gene = g2
                                data = {}
                                rje.combineDict(data,self.obj['GeneMap'].getGeneData(hgnc),overwrite=False,replaceblanks=True)
                                if rje.dictValues(data,'EnsDesc',valtype='str'): data['Desc'] = data['EnsDesc']
                                if rje.dictValues(data,'Desc',valtype='str'): desc = data['Desc']
                                try: hubx = len(self.dict['PPI'][gene])
                                except: hubx = 0
                                nodes[gene] = {'HUBSIZE':hubx,'DESCRIPTION':desc,'CONTROL':'Extra'}
                                sample = 'X%d' % (ex + 1)
                                while sample in rje.sortKeys(self.dict['Samples']): sample = 'X%s' % sample
                                nodes[gene]['SAMPLE'] = self.xgGeneSample(gene,sample)
                                try: nodes[gene]['ENSLOCI'] = data['EnsLoci']
                                except: self.log.printLog('#ENS','No EnsLoci data for %s' % gene)
                                try: nodes[gene]['ENSGENE'] = data['EnsEMBL']
                                except: self.log.printLog('#ENS','No EnsEMBL gene data for %s' % gene)
                    self.printLog('\r#PPI','Expanded PPI network (Lvl %d): %s nodes added' % (ex+1,len(next)))
                    
                ## ~ [2c] Update edge data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                edges = xgmml.dict['Edge'] = {'pp':{},'hom':{}}
                for g1 in xnodes[:-1]:
                    ## ~ [2c-i] PPI edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for g2 in xnodes[xnodes.index(g1)+1:]:
                        if g1 in self.dict['PPI'] and g2 in self.dict['PPI'][g1]:
                            if self.opt['XGComplex'] and self.getEvidence(g1,g2).lower().find('complex') < 1: continue
                            edges['pp'][(g1,g2)] = {'EVIDENCE':self.getEvidence(g1,g2)}       
                    ## ~ [2c-ii] Homology ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if self.opt['GABLAM'] and not self.stat['XGExpand']:
                            dis = self.obj['Homology'].getDis(g1,g2)
                            if dis:
                                dis = 100.0 * (1.0 - dis)
                                if dis >= self.obj['GABLAM'].stat['GABLAMCut']:
                                    edges['hom'][(g2,g1)] = {'GABLAM':dis,'EVIDENCE':'GABLAM'}
            
                ## ~ [2d] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                xgmml.saveXGMML(xfile)
            return xgmml
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def sampleCytoscape(self):  ### Outputs sample-specific SIF files (and combined)
        '''Outputs sample-specific SIF files (and combined).'''
        return self.log.printLog('#SIF','Old Cytoscape SIF output discontinued in V2.0')
#########################################################################################################################
    def association(self):  ### Performs Association analyses
        '''Performs Association analyses.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#ASSOC','Setting up Association analyses.')
            experiments = self.list['Experiments']
            if not experiments: experiments = rje.sortKeys(self.dict['Samples'])
            ## ~ [0a] ~ Expt/Gene dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            exptgene = {}   # Make dictionary of {expt:gene}
            geneexpt = {}   # Make dictionary of {gene:expt}
            for sample in experiments:
                exptgene[sample] = self.dict['Datasets'][sample][0:]
                for gene in exptgene[sample]:
                    if gene not in geneexpt: geneexpt[gene] = []
                    geneexpt[gene].append(sample)
            if not geneexpt: return self.printLog('#ASSOC','No Experiment Genes for association')
            if self.opt['AssCombo']: (exptgene,geneexpt) = self.assCombo(exptgene,geneexpt)
            allbyall = {}; afile = '%s.allbyall.tdt' % self.info['Basefile']
            if os.path.exists(afile): allbyall = rje.dataDict(self,afile,['Gene'])
            self.printLog('#AVA','Read AllByAll data for %s genes' % rje.integerString(len(allbyall)))
            go2go = {}; gx = 0.0; gtot = len(geneexpt)
            for g1 in geneexpt:
                self.progLog('\r#GO','Making GO sharing dictionary: %.2f%%' % (gx/gtot)); gx += 100.0
                g1go = self.getGO(g1)   #!# Replace this with precalc dictionary #!#
                go2go[g1] = {}
                for g2 in geneexpt:
                    go2go[g1][g2] = 0
                    for go in self.getGO(g2):
                        if go in g1go: go2go[g1][g2] += 1     
            self.printLog('\r#GO','Making GO sharing dictionary: %.2f%%' % (gx/gtot),log=False)
            ## ~ [0b] ~ Reduced PPI dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppi = {} # restricted PPI dictionary of just genes being considered
            ppi2 = {} # Genes with common interaction partner (not including one of them) - does not need to be in expt genes
            gx = 0.0; gtot = len(geneexpt)
            for gene in rje.sortKeys(geneexpt):
                self.progLog('\r#ASSOC','Reducing PPI dictionaries: %.2f%%' % (gx/gtot)); gx += 100.0
                ppi[gene] = []
                ppi2[gene] = []
                if gene not in self.dict['PPI']: continue
                for spoke in self.dict['PPI'][gene]:
                    if gene == spoke: continue
                    if spoke in geneexpt: ppi[gene].append(spoke)
                    for common in self.dict['PPI'][spoke]:
                        if common in [gene,spoke]: continue
                        if common in geneexpt: ppi2[gene].append(common); break
            self.printLog('\r#ASSOC','Reducing PPI dictionaries: %.2f%%' % (gx/gtot),log=False)
            ## ~ [0c] ~ Setup randomiser ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            random.seed(self.stat['RandSeed'])  # Random seed
            shareprofile = {}   # {No. expts: No. genes}
            for i in range(len(exptgene)): shareprofile[i+1] = 0
            for gene in geneexpt:
                shareprofile[len(geneexpt[gene])] += 1
            ## ~ [0d] ~ Setup output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfile = '%s.association.tdt' % self.info['Basefile']
            outhead = ['Expt1','Expt2','N','Overlap','Direct','Indirect','Complex','PPI','GO']
            for head in ['Overlap','Direct','Indirect','Complex','PPI','GO']:
                outhead = outhead[:outhead.index(head)+1] + ['r%s' % head, 'se%s' % head, 'p%s' % head] + outhead[outhead.index(head)+1:]
            ### ~ [1] ~ Calculate Gene and PPI associations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            realdata = self.associationStats(exptgene,geneexpt,ppi,ppi2,allbyall,go2go,log=True)
            outdata = {}
            for e1 in rje.sortKeys(exptgene):
                for e2 in rje.sortKeys(exptgene):
                    if e1 != e2 and self.opt['SelfOnly']: continue
                    outdata[(e1,e2)] = realdata[e1][e2]
                    outdata[(e1,e2)]['Expt1'] = e1
                    outdata[(e1,e2)]['Expt2'] = e2
                    for head in ['Overlap','Direct','Indirect','Complex','PPI','GO']: outdata[(e1,e2)]['r%s' % head] = []
                    for head in outhead:
                        if head not in outdata[(e1,e2)]: outdata[(e1,e2)][head] = 0
            ### ~ [2] ~ Randomise genes within experiments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            r = 0
            while r < self.stat['RandNum']:
                self.progLog('\r#RAND','Randomising experiments: %.2f%%' % (100.0*r/self.stat['RandNum'])); r += 1
                retry = 5
                while retry:
                    try:
                        elist = []
                        for e in exptgene: elist += [e] * len(exptgene[e])
                        self.deBug(elist)
                        elist = rje.randomList(elist)
                        glist = rje.randomList(geneexpt.keys())
                        self.deBug(elist); self.deBug(glist); self.deBug(shareprofile)
                        newgeneexpt = {}
                        newexptgene = {}
                        for e in exptgene: newexptgene[e] = []
                        for i in rje.sortKeys(shareprofile,revsort=True):
                            self.deBug(i)
                            for x in range(shareprofile[i]):
                                newg = glist.pop(0)
                                newgeneexpt[newg] = []
                                for n in range(i):
                                    err = len(elist)
                                    while elist[0] in newgeneexpt[newg]: 
                                        elist.append(elist.pop(0))
                                        err -= 1
                                        if err < 1:
                                            self.deBug(elist); self.deBug(newgeneexpt[newg]); self.deBug(glist)
                                            raise ValueError  # Cycled through whole list!
                                    e = elist.pop(0)
                                    newgeneexpt[newg].append(e)
                                    newexptgene[e].append(newg)
                        break
                    except ValueError: retry -= 1
                    except: raise
                if not retry: raise ValueError
                # Calculate Random Expectation # 
                rdata = self.associationStats(newexptgene,newgeneexpt,ppi,ppi2,allbyall,go2go)
                for e1 in rje.sortKeys(exptgene):
                    for e2 in rje.sortKeys(exptgene):
                        if e1 != e2 and self.opt['SelfOnly']: continue
                        for stat in ['Overlap','Direct','Indirect','Complex','PPI','GO']:
                            outdata[(e1,e2)]['r%s' % stat].append(rdata[e1][e2][stat])
                            if outdata[(e1,e2)][stat] <= rdata[e1][e2][stat]: outdata[(e1,e2)]['p%s' % stat] += 1
            self.printLog('\r#RAND','Randomising experiments complete.')
            ### ~ [3] ~ Convert to p-Values and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.delimitedFileOutput(self,outfile,outhead,rje_backup=True)
            for e1 in rje.sortKeys(exptgene):
                for e2 in rje.sortKeys(exptgene):
                    if e1 != e2 and self.opt['SelfOnly']: continue
                    for stat in ['Overlap','Direct','Indirect','Complex','PPI','GO']:
                        (m,s) = rje.meanse(outdata[(e1,e2)]['r%s' % stat])
                        outdata[(e1,e2)]['r%s' % stat] = m
                        outdata[(e1,e2)]['se%s' % stat] = s
                        outdata[(e1,e2)]['p%s' % stat] /= float(self.stat['RandNum'])
                    rje.delimitedFileOutput(self,outfile,outhead,datadict=outdata[(e1,e2)])
            self.printLog('\r#ASSOC','Association output to %s.' % outfile)
        except: self.errorLog('Error during Pingu Association analysis')            
#########################################################################################################################
    def assCombo(self,exptgene,geneexpt):   ### Combines genes sharing samples into new samples
        '''Combines genes sharing samples into new samples.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newexptgene = {}
            outhead = ['Sample','Identifier']
            outfile = '%s.asscombo.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,outfile,outhead,rje_backup=True)
            ### ~ [1] ~ Redefine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gene in rje.sortKeys(geneexpt):
                geneexpt[gene].sort(); expt = string.join(geneexpt[gene],'|')
                geneexpt[gene] = [expt]
                if expt not in newexptgene: newexptgene[expt] = []
                newexptgene[expt].append(gene)
                rje.delimitedFileOutput(self,outfile,outhead,datadict={'Sample':expt,'Identifier':gene})
            self.opt['SelfOnly'] = True
            return (newexptgene,geneexpt)
        except: self.errorLog('Pingu.assCombo error'); raise
#########################################################################################################################
    def associationStats(self,exptgene,geneexpt,ppi,ppi2,allbyall,go2go,log=False):  ### Calculates association stats
        '''Calculates association stats: Overlap, Direct, Indirect'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = {}
            for expt in exptgene:
                data[expt] = {}
                for expt2 in exptgene:
                    if expt != expt2 and self.opt['SelfOnly']: continue
                    data[expt][expt2] = {'Overlap':len(rje.listIntersect(exptgene[expt],exptgene[expt2])),
                                         'Direct':0,'Indirect':0,'Complex':0,'GO':0,'PPI':0.0,
                                         'N':max(len(exptgene[expt]),len(exptgene[expt2]))}
            ### ~ [1] ~ PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                    gx = 0.0; gtot = len(exptgene[expt])
                    for g1 in exptgene[expt]:
                        if log: self.progLog('\r#ASSOC','Association stats %s-%s: %.2f%%' % (expt,expt2,gx/gtot)); gx += 100.0
                        if g1 in exptgene[expt2]: continue  # Only genes not shared
                        for g2 in exptgene[expt2]:
                            if self.opt['NoShare'] and expt != expt2 and g2 in exptgene[expt]: continue   # Only genes not shared
                            if g2 in ppi[g1]: data[expt][expt2]['Direct'] += 1
                            if g2 in ppi2[g1]: data[expt][expt2]['Indirect'] += 1
                            if g2 in ppi[g1] or g2 in ppi2[g1]: data[expt][expt2]['Complex'] += 1
                            if g1 in allbyall and g2 in allbyall[g1]:
                                try: data[expt][expt2]['PPI'] += 1 / string.atof(allbyall[g1][g2])
                                except: self.deBug(allbyall); pass
            ### ~ [2] ~ GO Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                            if g1 in go2go and g2 in go2go[g1]: data[expt][expt2]['GO'] += go2go[g1][g2]
                    if log: self.printLog('\r#ASSOC','Association stats %s-%s: %.2f%%' % (expt,expt2,gx/gtot)); gx += 100.0
            return data
        except: self.errorLog('Error during Pingu AssociationStats'); raise
#########################################################################################################################
    ### <8> ### Sequence data generation                                                                                #
#########################################################################################################################
    def getEnsLoci(self):   ### Returns EnsLoci SeqList object, creating if necessary
        '''Returns EnsLoci SeqList object, creating if necessary.'''
        scmd = self.cmd_list+['seqin=%s' % self.info['EnsLoci'],'autoload=T','seqnr=F','accnr=F','align=False']
        if not 'EnsLoci' in self.obj or not self.obj['EnsLoci']: self.obj['EnsLoci'] = rje_seq.SeqList(self.log,scmd)
        if not 'EnsSeqNameDict' in self.obj or not self.obj['EnsSeqNameDict']: self.dict['EnsSeqNameDict'] = self.obj['EnsLoci'].seqNameDic()
        return self.obj['EnsLoci']
#########################################################################################################################
    def sampleSeqFiles(self):   ### Generates a fasta file per sample *.sample.fas
        '''Generates a fasta file per sample *.sample.fas.'''
        try:### ~ [1] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['GeneMap']: return self.log.errorLog('Cannot map EnsLoci without GeneMap.', printerror=False)
            #ensloci = self.getEnsLoci()
            ensloci = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F'])
            #seqdict = self.dict['EnsSeqNameDict']
            #if not seqdict: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)

            ### ~ [2] Generate Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sample in self.dict['Datasets']:
                sfile = '%s.%s.fas' % (self.info['Basefile'],sample)
                sseqs = []
                for gene in self.dict['Datasets'][sample]:
                    try: sseqs += self.obj['GeneMap'].getEnsLoci(gene)
                    except:
                        self.errorLog('Oh, man!')
                        if self.obj['GeneMap'].getGeneData(gene): self.log.printLog('#SEQ','No EnsLoci sequence for "%s"' % gene)
                        else: self.log.printLog('#MAP','No GeneMAP for "%s"' % gene)
                        self.deBug('%s: %s' % (gene,self.obj['GeneMap'].getGeneData(gene)))
                        self.deBug(self.obj['GeneMap'].dict['EnsLoci'])
                self.log.printLog('#SEQ','%d %s genes mapped to %d EnsLoci proteins' % (len(self.dict['Datasets'][sample]),sample,len(sseqs)))
                if sseqs: ensloci.saveFasta(seqs=sseqs,seqfile=sfile)

            ### ~ [3] Combined Fasta File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfile = '%s.Combined.fas' % self.info['Basefile']
            sseqs = []
            for gene in self.list['Genes']:
                try: sseqs += self.obj['GeneMap'].getEnsLoci(gene)   ## Update with GeneMap EnsLoci dictionary? 
                except: pass
            self.log.printLog('#SEQ','%d Combined genes mapped to %d EnsLoci proteins' % (len(self.list['Genes']),len(sseqs)))
            if sseqs: ensloci.saveFasta(seqs=sseqs,seqfile=sfile)
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def domainPPISeqFiles(self):    ### Generates a fasta file per domain PPI dataset *.ppidom.fas
        '''Generates a fasta file per domain PPI dataset *.ppidom.fas.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['DomPPIDir'].lower() in ['none','','none/','none\\']: return
            if not self.obj['GeneMap']: return self.log.errorLog('Cannot map EnsLoci without GeneMap.', printerror=False)
            rje.mkDir(self,self.info['DomPPIDir'])
            #!# if self.info['EnsDat'].lower() not in ['none','','none/','none\\']: return self.domainPPIDatFiles() #!#
            ## ~ [1a] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapper = self.obj['GeneMap']
            try: ensloci = mapper.obj['EnsLoci']
            except:
                try:
                    ensloci = mapper.info['EnsLoci']
                    if ensloci.seqNum() < 1: raise ValueError
                except:
                    mapper.loadEnsLoci(self.info['EnsLoci'])
                    ensloci = mapper.obj['EnsLoci']
            if not ensloci: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            ensloci.log = self.log
            ## ~ [1b] Setup Summary file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sumfile = '%s%s.domainppi.tdt' % (self.info['ResDir'],self.info['Basefile'])
            shead = ['Domain','Hubs','PPI','SeqNum','Missing']
            rje.delimitedFileOutput(self,sumfile,shead,rje_backup=True)
            missing = []

            ### ~ [2] Generate Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in rje.sortKeys(self.dict['DomPPI']):
                sfile = self.info['DomPPIDir'] + '%s.%sdom.fas' % (string.replace(hub,'/','-'),self.info['FasID'])
                sseqs = []
                data = {'Domain':hub,'Hubs':len(self.dict['DomHubs'][hub]),'PPI':len(self.dict['DomPPI'][hub]),'Missing':[]}
                for gene in self.dict['DomPPI'][hub]:
                    try:
                        eseq = mapper.getEnsLoci(gene)[0]
                        if eseq not in sseqs: sseqs.append(eseq)
                    except: data['Missing'].append(gene)
                data['SeqNum'] = len(sseqs)
                for gene in data['Missing']:
                    if gene not in missing: missing.append(gene)
                data['Missing'].sort()
                data['Missing'] = string.join(data['Missing'],',')
                rje.delimitedFileOutput(self,sumfile,shead,datadict=data)
                self.log.printLog('#SEQ','%d %s genes mapped to %d EnsLoci proteins (%d hubs)' % (len(self.dict['DomPPI'][hub]),hub,len(sseqs),len(self.dict['DomHubs'][hub])))
                if sseqs:
                    try:
                        if self.opt['AccOnly']: ensloci.saveAcc(seqs=sseqs,accfile='%s.acc' % sfile[:-4])
                        else: ensloci.saveFasta(seqs=sseqs,seqfile=sfile)
                    except: self.errorLog(rje_zen.Zen().wisdom())
            if missing:
                missing.sort()
                self.errorLog('Problem mapping %s DomPPI genes onto EnsLoci' % rje.integerString(len(missing)),printerror=False)
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def combinedPPISeqFiles(self):   ### Generates a fasta file per combined PPI dataset *.ppi.fas
        '''Generates a fasta file per combined PPI dataset *.ppi.fas.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['PPIOutDir'].lower() in ['none','','none/','none\\']: return
            if not self.obj['GeneMap']: return self.log.errorLog('Cannot map EnsLoci without GeneMap.', printerror=False)
            rje.mkDir(self,self.info['PPIOutDir'])
            if self.info['EnsDat'].lower() not in ['none','','none/','none\\']: return self.combinedPPIDatFiles()
            ## ~ [1a] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapper = self.obj['GeneMap']
            try: ensloci = mapper.obj['EnsLoci']
            except:
                try:
                    ensloci = mapper.info['EnsLoci']
                    if ensloci.seqNum() < 1: raise ValueError
                except:
                    mapper.loadEnsLoci(self.info['EnsLoci'])
                    ensloci = mapper.obj['EnsLoci']
            if not ensloci: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            ensloci.log = self.log
            ## ~ [1b] Setup Summary file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sumfile = '%s%s.combinedppi.tdt' % (self.info['ResDir'],self.info['Basefile'])
            shead = ['Gene','PPI','SeqNum','Missing']
            rje.delimitedFileOutput(self,sumfile,shead,rje_backup=True)
            missing = []

            ### ~ [2] Generate Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in rje.sortKeys(self.dict['PPI']):
                sfile = self.info['PPIOutDir'] + '%s.%s.fas' % (string.replace(hub,'/','-'),self.info['FasID'])
                sseqs = []
                data = {'Gene':hub,'PPI':len(self.dict['PPI'][hub]),'Missing':[]}
                for gene in self.dict['PPI'][hub]:
                    try:
                        eseq = mapper.getEnsLoci(gene)[0]
                        if eseq not in sseqs: sseqs.append(eseq)
                    except: data['Missing'].append(gene)
                data['SeqNum'] = len(sseqs)
                for gene in data['Missing']:
                    if gene not in missing: missing.append(gene)
                data['Missing'].sort()
                data['Missing'] = string.join(data['Missing'],',')
                rje.delimitedFileOutput(self,sumfile,shead,datadict=data)
                self.log.printLog('#SEQ','%d %s genes mapped to %d EnsLoci proteins' % (len(self.dict['PPI'][hub]),hub,len(sseqs)))
                if sseqs:
                    try:
                        if self.opt['AccOnly']: ensloci.saveAcc(seqs=sseqs,accfile='%s.acc' % sfile[:-4])
                        else: ensloci.saveFasta(seqs=sseqs,seqfile=sfile)
                    except: self.log.errorLog(rje_zen.Zen().wisdom())
            if missing:
                missing.sort()
                self.log.errorLog('Problem mapping %s PPI genes onto EnsLoci' % rje.integerString(len(missing)),printerror=False)
                open('%s%s.no_ensloci.txt' % (self.info['ResDir'],self.info['Basefile']),'w').write(string.join(missing,'\n'))
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def combinedPPIDatFiles(self):   ### Generates a DAT file per combined PPI dataset *.ppi.dat
        '''Generates a DAT file per combined PPI dataset *.ppi.dat.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(self.info['EnsDat']): return self.log.errorLog('EnsDat path not found!', printerror=False)
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
            uniprot.info['UniPath'] = self.info['EnsDat']

            ### ~ [2] Generate DAT Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in self.dict['PPI']:
                dfile = self.info['PPIOutDir'] + '%s.%s.dat' % (hub,self.info['FasID'])
                sacc = []
                for gene in self.dict['PPI'][hub]:
                    try: sacc.append(self.obj['GeneMap'].getGeneData(gene)['EnsEMBL'])
                    except: pass
                self.printLog('#SEQ','%d %s genes mapped to %d EnsEMBL entries' % (len(self.dict['PPI'][hub]),hub,len(sacc)))
                if sacc:
                    try:
                        uniprot.readUniProt(clear=True,acclist=sacc)
                        uniprot.saveUniProt(dfile)
                    except: self.log.errorLog(rje_zen.Zen().wisdom())            
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def geneSeqFile(self,filename=None):    ### Generates a single fasta file with Gene names as protein IDs
        '''Generates a single fasta file with Gene names as protein IDs.'''
        try:### ~ [1] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not filename: filename = '%s.genes.fas' % self.info['Basefile']
            if not self.obj['GeneMap']: return self.errorLog('Cannot map EnsLoci without GeneMap.', printerror=False)
            mapper = self.obj['GeneMap']
            try: ensloci = mapper.obj['EnsLoci']
            except:
                try:
                    ensloci = mapper.info['EnsLoci']
                    if ensloci.seqNum() < 1: raise ValueError
                except:
                    mapper.loadEnsLoci(self.info['EnsLoci'])
                    ensloci = mapper.obj['EnsLoci']
            if not ensloci: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            ensloci.log = self.log
            geneseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
            
            ### ~ [2] Generate Fasta File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gene in self.list['Genes']:
                try:
                    enseq = mapper.getEnsLoci(gene)[0]
                    geneseq._addSeq('%s %s' % (gene,enseq.info['Description']),enseq.info['Sequence'])
                except: self.log.printLog('#MISSING!','Cannot map sequence for %s' % gene)
            try: geneseq.saveFasta(seqfile=filename)
            except: self.log.errorLog(rje_zen.Zen().wisdom())            
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    ### <9> ### Special Analysis Methods                                                                                #
#########################################################################################################################
    def qSLiMFinder(self):  ### Performs an analysis for shared motifs between primary interactors & samples
        '''
        Performs an analysis for shared motifs between primary interactors of those genes identified in a given sample
        and the original sequence used for the pulldown. FILE should be a fasta file where names of the sequences match
        Sample names. Datasets will then be formed that contain that sequence plus the primary PPI of each gene in that
        sample as a dataset named SAMPLE_GENE.fas in a directory RESDIR/SLiMFinder.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = self.cmd_list+['seqin=%s' % self.info['QSLiMFinder'],'autoload=T','autofilter=F','accnr=F','seqnr=F']
            qseq = rje_seq.SeqList(self.log,scmd)
            if not qseq.seq: return self.log.printLog('#QSF','Not performing QSLiMFinder analysis')
            qdict = qseq.seqNameDic(key='Max')
            sfdir = rje.makePath(self.info['ResDir'] + 'QSLiMFinder/')
            rje.mkDir(self,sfdir)
            ## ~ [1a] Setup EnsLoci file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mapper = self.obj['GeneMap']
            try: ensloci = mapper.obj['EnsLoci']
            except:
                try:
                    ensloci = mapper.info['EnsLoci']
                    if ensloci.seqNum() < 1: raise ValueError
                except:
                    mapper.loadEnsLoci(self.info['EnsLoci'])
                    ensloci = mapper.obj['EnsLoci']
            if not ensloci: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            
            ### ~ [2] Make Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sample in rje.sortKeys(self.dict['Samples']):
                sfcmd = ['resdir=%s' % sfdir,'batch=%s*.%s.fas' % (sfdir,sample),'query=%s' % sample,'pickle=T']
                open('%s_qslimfinder.ini' % sample,'w').write(string.join(sfcmd,'\n'))
                if sample not in qdict:
                    if sample in ['Control','PathFinder']: continue
                    self.log.errorLog('Cannot find %s in %s' % (sample,self.info['QSLiMFinder']),printerror=False)
                    continue
                for gene in self.dict['Datasets'][sample]:
                    base = '%s.%s' % (gene,sample)
                    sfile = sfdir + base + '.fas'
                    if not gene in self.dict['PPI']:
                        self.log.printLog('#QSF','No QSLiMFinder analysis for %s - no %s PPI' % (base,gene))
                        continue
                    pseqs = []
                    for ppi in self.dict['PPI'][gene]:
                        try:
                            eseq = mapper.getEnsLoci(ppi)[0]
                            if eseq not in pseqs: pseqs.append(eseq)
                        except:
                            if self.obj['GeneMap'].getGeneData(ppi): self.log.printLog('#SEQ','No EnsLoci sequence for %s' % ppi)
                            else: self.log.printLog('#MAP','No GeneMap for %s' % ppi)
                            self.deBug('!')
                    if pseqs:
                        pseqs = [qdict[sample]] + pseqs
                        ensloci.saveFasta(seqs=pseqs,seqfile=sfile)

            ### ~ [3] Run SLiMFinder? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['Interactive'] < 0 or rje.yesNo('Run SLiMFinder?'):
                first = True
                for sample in rje.sortKeys(self.dict['Samples']):
                    sfcmd = self.cmd_list+['ini=%s_qslimfinder.ini' % sample]
                    if first: first = False
                    else: sfcmd += ['append=T']
                    cmd_list = rje.getCmdList(sfcmd,info=qslimfinder.makeInfo())
                    try: qslimfinder.QSLiMFinder(self.log,cmd_list).run()
                    except: self.errorLog('QSLiMFinder problem')
                    #x#sf2cmd = cmd_list + ['batch=%s*%s*fas' % (self.info['ResDir'],sample),'query=None','append=T']
                    #x#try: slimfinder.SLiMFinder(self.log,sf2cmd).run()
                    #x#except: self.errorLog('SLiMFinder problem')
        except: self.errorLog('Ahh, shucks. Pingu.QSLiMFinder() busted.')            
#########################################################################################################################
    ### <10> ### R Visualisation Methods                                                                                #
#########################################################################################################################
    def interactomePNG(self,sample):   ### Generates interactome visualisation for given sample                     |3.0|
        '''Generates interactome visualisation for given sample.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # This visualisation needs a UPC file and a dis.tdt file
            ## ~ [1a] ~ Generate UPC groupings and DIS file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            basefile = '%s.%s' % (self.info['Basefile'],sample)
            sfile = '%s.fas' % basefile
            dataset = rje.baseFile(sfile,strip_path=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return self.printLog('#PNG','%s.png exists. (Force=F)' % basefile)
            if not sfile: self.sampleSeqFiles()
            scmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F','seqin=%s' % sfile]
            slimcore = rje_slimcore.SLiMCore(self.log,scmd)
            slimcore.info['ResDir'] = self.info['ResDir']
            slimcore.obj['SeqList'] = rje_seq.SeqList(self.log,scmd)
            slimcore.setupBasefile()
            slimcore.makeUPC()
            disfile = '%s.dis.tdt' % basefile       # Make tree from this
            ## ~ [1b] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s.upc' % basefile           # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = string.split(uline)[3:]
                upc.sort()
                uplist.append(string.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = string.split(upc)
                reorder += uprots     #!# Then reorder sequences! #!#
                if len(uprots) > 1:
                    for u in uprots: upid[self.geneMap(u)] = upx
                    upx += 1
                else: upid[self.geneMap(uprots[0])] = 0

            ### ~ [2] ~ Load distance matrix and make Tree/Heatmap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disdat = rje.dataDict(self,disfile,['SEQ'])
            dismat = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            for obj1 in reorder:
                for obj2 in reorder:
                    dis = (string.atof(disdat[obj1][obj2]) + string.atof(disdat[obj2][obj1])) / 2.0
                    dismat.addDis(self.geneMap(obj1),self.geneMap(obj2),dis)
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            nsftree = dismat.upgma()
            open('%s.nsf' % basefile,'w').write(nsftree)
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            upgma.rTree('%s.tree.csv' % basefile,seqname='short')
            reduced = upgma._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            dismat.info['Name'] = '%s interactome GABLAM' % sample
            dismat.saveMatrix(reduced,basefile+'.heatmap.tdt',delimit='\t')

            ### ~ [3] ~ Add PPI links between spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            names = []
            for prot in reorder: names.append(self.geneMap(prot))
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for p1 in names:
                datadict = {}
                for p2 in names:
                    try: ppi = p2 in self.dict['PPI'][p1]  #False
                    except: ppi = False
                    if ppi: datadict[p2] = 1
                    else: datadict[p2] = 0
                rje.delimitedFileOutput(self,pfile,names,datadict=datadict)
            datadict = {}
            for p1 in names: datadict[p1] = upid[p1]
            rje.delimitedFileOutput(self,pfile,names,datadict=datadict)

            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "interactome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%srje_call.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt" 2>&1' % (rslimjim,basefile)
            self.printLog('#RCALL',rcmd)
            problems = os.popen(rcmd).read()
            if problems: self.errorLog(problems,printerror=False)
            if not os.path.exists('%s.png' % basefile): self.errorLog('%s.png not created' % basefile,printerror=False)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt','.nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

        except: self.errorLog('Pingu.interactomePNG() stuffed up')            
#########################################################################################################################
### End of SECTION II: PINGU Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
class GeneMap(rje_genemap.GeneMap):
    '''Just for pickling.'''
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
    try: PINGU(mainlog,cmd_list).run()
        
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
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
