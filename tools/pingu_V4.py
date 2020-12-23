#!/usr/bin/python

# See below for name and description
# Copyright (C) 2013 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       PINGU
Description:  Protein Interaction Network & GO Utility
Version:      4.10.0
Last Edit:    21/05/19
Copyright (C) 2013  Richard J. Edwards - See source code for GNU License Notice

Function:
    PINGU (Protein Interaction Network & GO Utility) is designed to be a general utility for Protein Protein Interaction
    (PPI) and Gene Ontology (GO) analysis. Earlier versions of PINGU contained a lot of the code for processing PPI and
    GO data, which have subsequently been moved to `rje_ppi.py` and `rje_go.py` libraries.

    PINGU 3.x was dominated by code to compile PPI data from multiple databases and map onto sequences from different
    sources, in combination with `rje_dbase.py` database downloads and processing. There was a substantial amount of
    code for mapping data with different IDs, including peptide-based MS Ensembl identifications on to HGNC gene
    identifiers. Some of this code is now handled by `rje_genemap.py` whilst some of it has been depracated due to newer
    (better) datasets and/or a shift in focus of the Edwards lab.

    PINGU 4.x is designed to work in a more streamlined fashion with a more controlled subset of data, making
    documentation and re-use a bit simpler and clearer. Many of the older functions are therefore run using
    `pingu_V3.py`, in which case a `#PINGU` log statement will be generated. Some older functions will only be possible
    by running `pingu_V3.py` directly.

    PINGU 4.0 is designed to work with HINT interaction data and UniProt sequences. The initial PPI download and
    compilation is based on `SLiMBench`. This has been updated in 4.9.0 following some changed to HINT downloads
    (http://hint.yulab.org/download/) - there may be some additional unexpected/unwelcome consequences of these changes.

    PINGU 4.1 updated the PPI compilation methods of PINGU 3.x, which can be triggered using `ppicompile=T`. This will
    need a database cross-reference file (`xrefdata=LIST`).

    PINGU 4.2 will download and use HGNC as a database xref file if xrefdata=HGNC. Clearly, this will only work for
    human data. Note that HINT is mapped to genes via Uniprot entries and does not use the xrefdata table.

    PINGU 4.3 add domain-based domppi dataset generation (`domppi=T`). This uses Pfam domain composition from Uniprot to
    generate datasets of proteins that interact with hubs sharing a domain.

    PINGU 4.4.x replaced `ppicompile=T` with `ppicompile=CDICT`. `HPRD`, `HINT` and `Reactome` will be recognised and parsed
    using custom methods: `hprd=PATH` and `reactome=FILE` must be set; HINT data will be read from `sourcepath=PATH/`.
    Otherwise, entries will be treated as files (wildcard lists allowed) and either parsed as a pairwise PPI file (`Hub`
    and `Spoke` fields found) else a MITAB file (see rje_mitab for advanced field settings). The compiled PPI data
    will be output to `BASEFILE.pairwise.tdt` and used as `ppisource=X` for additional processing/output.

    Default mapping fields for XRef mapping are: `Secondary,Ensembl,Aliases,Accessions,RefSeq,Previous Symbols,Synonyms`.
    `Secondary` will be added from Uniprot data if missing from the XRef table. `unifield=X` will also be added to the
    map fields if not included.

    PINGU 4.6.x fixed/updated the PPI Fasta output methods (`ppifas=T`). These will output to a directory named after the
    `ppisource` file and `ppispec`. Each hub gene will produce a fasta file, `gene.fasid.fas` where `fasid` is set by the
    `ppisource` unless changed with `fasid=X`. If `combineppi=T`, a single `spec.fasid.fas` file will be created. The
    `xhubppi=T` setting will generate a set of files containing spoke proteins that have x+ Hub interactors. Note that
    `*.1hub.fas` is essentially the same as the `combineppi=T` `spec.fasid.fas` file.

    PINGU 4.7 added `ppidbreport=T/F` output for PPI compilation, summarising the evidence codes and PPITypes read from
    different sources.

Commandline:
    ### ~ SOURCE DATA OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sourcepath=PATH/    : Will look in this directory for input files if not found ['SourceData/']
    sourcedate=DATE     : Source file date (YYYY-MM-DD) to preferentially use [None]
    ppisource=X         : Source of PPI data. (HINT/FILE) FILE needs 'Hub', 'Spoke' and 'SpokeUni' fields. ['HINT']
    ppispec=LIST        : List of PPI files/species/databases to generate PPI datasets from [HUMAN]
    download=T/F        : Whether to download files directly from websites where possible if missing [True]
    integrity=T/F       : Whether to quit by default if source data integrity is breached [True]
    xrefdata=LIST       : List of files with delimited data of identifier cross-referencing (see rje_xref) []
    unifield=X          : Uniprot accession number field identifier for xrefdata ['Uniprot']
    mapfields=LIST      : List of XRef fields to use for identifier mapping (plus unifield) [see docs]

    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    resdir=PATH         : Redirect output files/directories to specified directory [./]
    basefile=X          : Results file prefix [pingu]
    acconly=T/F         : Whether to output lists of Accession numbers only, rather than full fasta files [False]
    fasid=X             : Text ID for fasta files (*.X.fas) [default named after ppisource(+'-dom')]

    ### ~ PPI OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ppiout=FILE         : Save pairwise PPI file following processing (if rest=None) [None]
    ppifas=T/F          : Whether to output PPI fasta files [False]
    domppi=T/F          : Whether to generate Pfam Domain-based PPI files instead of protein-based PPI files [False]
    minppi=X            : Minimum number of PPI for file output [0]
    combineppi=T/F      : Whether to combine all spokes into a single fasta file [False]
    xhubppi=T/F         : Whether to generate PPI files of spokes interacting with X+ hubs [False]
    queryppi=FILE       : Load a file of 'Query','Hub' PPI and generate expanded PPI Datasets in PPI.*/ [None]
    queryseq=FILE       : Fasta file containing the Query protein sequences corresponding to QuerySeq [*.fas]
    allquery=T/F        : Whether to include all the new Queries from QueryPPI in all files for a given hub [True]

    ### ~ ADVANCED/DEV OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sourceurl=CDICT     : Dictionary of Source URL mapping (see code)

    ### ~ PPI COMPILATION/FILTERING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hublist=LIST        : List of hub genes to restrict pairwise PPI to []
    hubonly=T/F         : Whether to restrict pairwise PPI to those with both hub and spoke in hublist [False]
    hubfield=X          : Hub field to use for hublist=LIST [Hub]
    spokefield=X        : Spoke field to use for hublist=LIST hubonly=T [Spoke]
    ppicompile=CDICT    : List of db:file PPI Sources to compile and generate *.pairwise.tdt []
    ppidbreport=T/F     : Summary output for PPI compilation of evidence/PPIType/DB overlaps [True]
    symmetry=T/F        : Whether to enforce Hub-Spoke symmetry during PPI compilation [True]
    hprd=PATH           : Path to HPRD flat files [None]
    taxid=LIST          : List of NCBI Taxa IDs to use [9606]
    badppi=LIST         : PPI Types to be removed. Will only remove PPI if no support remains []
    goodppi=LIST        : Reduce PPI to those supported by listed types []
    baddb=LIST          : PPI Types to be removed. Will only remove PPI if no support remains []
    gooddb=LIST         : Reduce PPI to those supported by listed types []

    ### ~ OBSOLETE OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    biogrid=FILE        : BioGRID flat file [None]
    intact=FILE         : IntAct flat file [None]
    mint=FILE           : MINT flat file [None]
    reactome=FILE       : Reactome interactions flat file [None]
    dip=FILE            : DIP interactions flat file [None]
    domino=FILE         : Domino interactions flat file [None]
    evidence=FILE       : Mapping file for evidence terms [None]    #!# Not currently implemented! #!#
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_go, rje_obj, rje_ppi, rje_seqlist, rje_uniprot, rje_xgmml, rje_xref
import pingu_V3, rje_hprd, rje_biogrid, rje_mitab
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 4.0 - Initial Compilation based on code from SLiMBench and PINGU 3.9 (inherited as pingu_V3).
    # 4.1 - Adding compilation of PPI databases using new rje_xref V1.1 and older objects from PINGU V3.
    # 4.2 - Bug fixes for use of PPISource to create PPI databases. Add HGNC to sourcedata (xrefdata=HGNC)
    # 4.3 - Modified to use Pfam as hub field for DomPPI generation. Modified naming of PPI output after ppisource.
    # 4.4.0 - Converted ppicompile=T to ppicompile=LIST.
    # 4.5.0 - Added hublist=LIST : List of hub genes to restrict pairwise PPI to, and pairwise parsing.
    # 4.5.1 - Debugging missing identifiers and indexing speed. Added good and bad DB.
    # 4.5.2 - Fixed SIF output and changed names to sif-* for opening in browser.
    # 4.5.3 - Updated REST output.
    # 4.6.0 - Added hubonly=T/F : Whether to restrict pairwise PPI to those with both hub and spoke in hublist [False]
    # 4.6.1 - Fixed some ppifas=T/F bugs and added combineppi=T/F : Whether to combine all spokes into a single fasta file [False]
    # 4.6.2 - Added check/filter for multiple SpokeUni pointing to same sequence. (Compilation redundancy mapping failure!)
    # 4.6.3 - Fixed issue with 1:many SpokeUni:Spoke mappings messing up XHub.
    # 4.7.0 - Added ppidbreport=T/F : Summary output for PPI compilation of evidence/PPIType/DB overlaps [True]
    # 4.8.0 - Fixed report duplication issue and added additional summary output.
    # 4.9.0 - Updated HINT download and parsing details.
    # 4.9.1 - Fixed Pairwise parsing and filtering for more flexibility of input. Fixed fasid=X bug and ppiseqfile names.
    # 4.10.0 - Added hubfield and spokefield options for parsing hublist.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Create initial working version of program.
    # [ ] : screenddi=FILE  : Whether to screen out probably domain-domain interactions from file [None]
    # [ ] : nocomplex=T/F   : Perform crude screening of complexes (PPI triplets w/o homodimers) [False]
    # [Y] : Add addition of queries from a second PPI file.
    # [X] : Can the individual PPI data be replaced with PPISource/ downloads etc? (Some perhaps.)
    # [X] : Add automatic download of IntAct, BioGRID, MINT?, Reactome, Domino and 3DID. (May need unzipping.)
    # [ ] : Add parsing of 3DID.
    # [Y] : Add catching and summarising of missing accnum when ppisource is used rather than assuming all acc match spec
    # [Y] : Test running with unipath=PATH setting (rather than URL). -> Dropped!
    # [ ] : Need to get the CAEEL gene mapping working OK!
    # [ ] : Implement/Update the evidence=FILE evidence code mapping.
    # [Y] : Add goodppi and badppi for evidence filtering.
    # [Y] : Replace ppicompile with list of database types/files.
    # [ ] : Add making a rje_ppi.PPI object from PINGU.db('pairwise') as PPI.db('Edge') for outputs.
    # [Y] : Test with REST output.
    # [ ] : Properly sort out the setup and use of MapFields, incorporating XRef KeyID etc.
    # [ ] : Add expandppi=X to increase hublist. (Might be easier to read whole in and then remove)
    # [ ] : Rationalise use of ResDir.
    # [ ] : Tidy up obselete V3 code.
    # [ ] : Add hostpathogen=LIST PPI mode to extract one taxon match and one pathogen taxon match in the PPI pairs.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('PINGU', '4.10.0', 'May 2019', '2013')
    description = 'Protein Interaction Network & GO Utility'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',
                'See PINGU 3.x (pingu_V3.py) for additional functions',
                rje_obj.zen()]
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
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
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
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class PINGU(rje_obj.RJE_Object):
    '''
    PINGU Class. Author: Rich Edwards (2013).

    Str:str
    - Evidence = Mapping file for evidence terms [None]
    - FasID = Text ID for PPI fasta files (*.X.fas) ['ppi/domppi']
    - HubField=X          : Hub field to use for hublist=LIST [Hub]
    - SpokeField=X        : Spoke field to use for hublist=LIST hubonly=T [Spoke]
    - PPISource = Source of PPI data. (See documentation for details.) (HINT/FILE) ['HINT']
    - QueryPPI = Load a file of 'Query','Hub' PPI and generate expanded PPI Datasets in PPI.*/ [None]
    - QuerySeq = Fasta file containing the Query protein sequences corresponding to QuerySeq [*.fas]
    - ResDir = Redirect output files/directories to specified directory [./]
    - SourceDate = Source file date (YYYY-MM-DD) to preferentially use [None]
    - SourcePath = Will look in this directory for input files if not found ['SourceData/']
    - UniField = Uniprot accession number field identifier for xrefdata ['UniProt']
    - DIP = DIP interactions flat file [None]
    - Domino = Domino interactions flat file [None]
    - HPRD = Path to HPRD flat files [None]
    - BioGRID = BioGRID flat file [None]
    - IntAct = IntAct flat file [None]
    - MINT = MINT flat file [None]
    - Reactome = Reactome interactions flat file [None]
    - XRefData = Use to recognise XRefData=HGNC over-ride of XRefData=LIST

    Bool:boolean
    - AccOnly = Whether to output lists of Accession numbers only, rather than full fasta files [False]
    - AllQuery = Whether to include all the new Queries from QueryPPI in all files [True]
    - CombinePPI=T/F      : Whether to combine all spokes into a single fasta file [False]
    - DomPPI = Whether to generate Domain-based PPI files [False]
    - Download = Whether to download files directly from websites where possible if missing [True]
    - HubOnly = Whether to restrict pairwise PPI to those with both hub and spoke in hublist [False]
    - Integrity = Whether to quit by default if input integrity is breached [True]
    - PPICompile = Whether to compile PPI datasets from individual sources and generate *.pairwise.tdt [False]
    - PPIFas = Whether to output PPI files [False]
    - PPIDBReport = Summary output for PPI compilation of evidence/PPIType/DB overlaps [True]
    - Symmetry = Whether to enforce Hub-Spoke symmetry [True]
    - XHubPPI=T/F         : Whether to generate PPI files of spokes interacting with X+ hubs [False]

    Int:integer
    - MinPPI = Minimum number of PPI for file output [0]

    Num:float

    List:list
    - BadDB=LIST         : DB Types to be removed. Will only remove PPI if no support remains []
    - BadPPI=LIST         : PPI Types to be removed. Will only remove PPI if no support remains []
    - GoodDB=LIST        : Reduce PPI to those supported by listed DB []
    - GoodPPI=LIST        : Reduce PPI to those supported by listed types []
    - HubList=LIST        : List of hub genes to restrict pairwise PPI to []
    - MapFields=LIST      : List of XRef fields to use for identifier mapping (plus unifield) [see docs]
    - PPISpec = List of PPI files/species/databases to generate PPI datasets from ['HUMAN','YEAST']
    - TaxID=LIST          : List of NCBI Taxa IDs to use [9606]

    Dict:dictionary
    - PPI = dictionary of {Species:PPI object}
    - PPICompile = Whether to compile PPI datasets from individual sources and generate *.pairwise.tdt [False]
    - SourceURL = dictionary of {Data Source: download URL}
    - UniSpec = dictionary of {Species:UniProt Object}

    Obj:RJE_Objects
    - DB = rje_db.DBase object
    - UniProt = rje_uniprot.UniProt object
    - XRef = rje_xref.XRef object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FasID','HubField','SpokeField','QueryPPI','QuerySeq','ResDir','PPISource','SourceDate','SourcePath','XRefData',
                        'Evidence','DIP','Domino','HPRD','BioGRID','IntAct','MINT','Reactome']
        self.boollist = ['AccOnly','AllQuery','CombinePPI','DomPPI','Download','HubOnly','Integrity','PPICompile','PPIDBReport','PPIFas','Symmetry','XHubPPI']
        self.intlist = ['MinPPI']
        self.numlist = []
        self.listlist = ['BadDB','BadPPI','GoodDB','GoodPPI','HubList','MapFields','PPISpec','RestDB','TaxID']
        self.dictlist = ['PPI','PPICompile','SourceURL','UniSpec']
        self.objlist = ['DB','UniProt','XRef']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'PPISource':'HINT','SourcePath':'SourceData/','UniField':'Uniprot','ResDir':rje.makePath('./'),'HubField':'Hub','SpokeField':'Spoke'})
        self.setBool({'Download':True,'HubOnly':False,'Integrity':True,'PPIDBReport':True,'Symmetry':True})
        self.setInt({})
        self.setNum({})
        self.list['MapFields'] = string.split('Gene,Uniprot,UniprotID,Secondary,Ensembl,Aliases,Accessions,RefSeq,GenPept,Previous Symbols,Synonyms',',')
        self.list['PPISpec'] = ['HUMAN']
        self.list['TaxID'] = ['9606']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.dict['SourceURL'] = {'ELMClass':'http://elm.eu.org/elms/browse_elms.html?q=&submit=tsv',
             'ELMInstance':'http://elm.eu.org/elms/browse_instances.tsv?q=*&taxon=&instance_logic=',
             'ELMInteractors':'http://elm.eu.org/infos/interactions_as_tsv',
             'ELMDomains':'http://www.elm.eu.org/infos/browse_elm_interactiondomains.tsv',
             'DMIFile':'http://www.elm.eu.org/infos/browse_elm_interactiondomains.tsv',
             'HINT.HUMAN':'http://hint.yulab.org/download/HomoSapiens/binary/hq/',
             'HINT.YEAST':'http://hint.yulab.org/download/SaccharomycesCerevisiaeS288C/binary/hq/',
             'HINT.MOUSE':'http://hint.yulab.org/download/MusMusculus/binary/hq/',
             'HINT.DROME':'http://hint.yulab.org/download/DrosophilaMelanogaster/binary/hq/',
             'HINT.CAEEL':'http://hint.yulab.org/download/CaenorhabditisElegans/binary/hq/',
             'HGNC':'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_pub_refseq_ids&col=md_mim_id&col=md_prot_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit',
             'Uniprot.HUMAN':'http://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes&format=txt',
             'Uniprot.YEAST':'http://www.uniprot.org/uniprot/?query=organism:559292+AND+reviewed:yes&format=txt',
             'Uniprot.MOUSE':'http://www.uniprot.org/uniprot/?query=organism:10090+AND+reviewed:yes&format=txt',
             'Uniprot.DROME':'http://www.uniprot.org/uniprot/?query=organism:7227+AND+reviewed:no&format=txt',
             'Uniprot.CAEEL':'http://www.uniprot.org/uniprot/?query=organism:6239+AND+reviewed:no&format=txt'
             }
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['UniProt'] = rje_uniprot.UniProt(self.log,['unipath=url']+self.cmd_list)
        self.baseFile('pingu')
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
                self._cmdReadList(cmd,'str',['FasID','HubField','SpokeField','UniField','XRefData'])   # Normal strings
                self._cmdReadList(cmd,'date',['SourceDate'])   # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'path',['HPRD','SourcePath','ResDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['PPISource','QueryPPI','QuerySeq',
                                              'Evidence','DIP','Domino','BioGRID','IntAct','MINT','Reactome'])
                self._cmdReadList(cmd,'bool',['AccOnly','AllQuery','CombinePPI','DomPPI','Download','HubOnly',
                                              'Integrity','PPICompile','PPIDBReport','PPIFas','Symmetry','XHubPPI'])
                self._cmdReadList(cmd,'int',['MinPPI'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['HubList','MapFields','PPISpec','TaxID'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'lclist',['BadDB','BadPPI','GoodDB','GoodPPI'])  # List of LC strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['PPICompile','SourceURL']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStr('UniField') not in self.list['MapFields']: self.list['MapFields'].insert(0,self.getStr('UniField'))
        self.baseFile(self.baseFile())  # Make sure that Database etc. get right basefile
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False    #i# Includes PPICompile code
            self.debug('Run...')
            self.filterPPI()
            if self.getBool('XHubPPI'): self.xHubPPI()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ PPI Fasta output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('PPIFas'): self.ppiFas()
            ## ~ [2b] ~ Query PPI output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('QueryPPI'): self.queryPPI()

            self.generateRestOutput()
            ### ~ [3] ~ PINGU V3.x ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #else:
            #    self.printLog('#PINGU','No PINGU 4.x operations recognised: calling PINGU 3.x!')
            #    pingu_V3.PINGU(self.log,self.cmdlist).run()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def xrefDB(self):   ### Returns (and/or creates if required) XRef database table.
        '''Returns (and/or creates if required) XRef database table.'''
        if not self.obj['XRef'] or not self.obj['XRef'].db('xref'): self.setupXRef()
        return self.obj['XRef'].db('xref')
#########################################################################################################################
    def setupXRef(self):    ### Main class setup method.                                                      #V2.0
        '''Main class setup method.'''
        try:### ~ [0] Setup Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP XREF ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            xref = self.obj['XRef'] = rje_xref.XRef(self.log,self.cmd_list)
            ## ~ [0a] ~ Special HGNC processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('XRefData') == 'hgnc':
                self.setStr({'HGNC':'hgnc_download.tdt'})
                hgncfile = self.sourceDataFile('HGNC',expect=True,ask=False)
                if rje.exists(hgncfile):
                    xref.setStr({'KeyID':'Approved Symbol'})
                    xref.list['XRefData'] = [hgncfile]
                    xref.list['NewHeaders'] = string.split('HGNC,Gene,Name,Status,Previous Symbols,Synonyms,Chromosome,Accessions,Entrez,Ensembl,RefSeq,OMIM,Uniprot',',')
                    xref.list['MapFields'] = string.split('HGNC,Gene,Uniprot,Ensembl,RefSeq,Previous Symbols,Synonyms,Entrez',',')
            ## ~ [0b] ~ Update XRef entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sdb.addEntry({'Name':'XRef','Status':xref.setup()})
            sdb.data('XRef')['File'] = xref.basefile()
            xdb = xref.db('xref')
            if xdb: sdb.data('XRef')['Entries'] = xref.db('xref').entryNum()
            ufield = self.getStr('UniField')
            ## ~ [0c] ~ Add mapping of Uniprot secondary accnum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newxfile = None
            if xdb and ('Secondary' not in xdb.fields() or 'UniprotID' not in xdb.fields()):
                newxfile = '%s.xref.tdt' % self.basefile()
                if not self.force() and rje.exists(newxfile):
                    self.printLog('#XREF','Field "Secondary" or "UniprotID" not in XRef: loading %s (force=F).' % newxfile)
                    xref.list['XRefData'] = [newxfile]
                    xref.db().deleteTable(xdb)
                    xref.setup()
                    sdb.data('XRef')['File'] = newxfile
                    xdb = xref.db('xref')
            if xdb and ('Secondary' not in xdb.fields() or 'UniprotID' not in xdb.fields()):
                self.printLog('#XREF','Field "Secondary" or "UniprotID" not in XRef: adding Uniprot data.')
                xdb.addField('Secondary',after=ufield)
                xdb.addField('UniprotID',after=ufield)
                self.dict['UniSpec'] = {}
                accmap = {}     # Mapping of UniProt accession numbers to Primary AccNum
                for spec in self.list['PPISpec']:   # Download UniProt Data
                    self.setStr({'Uniprot.%s' % spec: 'uniprot.%s.dat' % spec})
                    # Download/Check UniProt DAT and fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    ppidat = self.sourceDataFile('Uniprot.%s' % spec,expect=True,ask=False)  # If this exists, will use
                    if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                    else:
                        ppuni = rje_uniprot.UniProt(self.log,self.cmd_list)
                        ppuni.setStr({'Name':ppidat})
                        ppuni.readUniProt()
                        if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                    accmap = rje.combineDict(accmap,ppuni.accDictFromEntries())     # Map all accession numbers onto species UniProt entries
                for entry in xdb.entries():
                    try: uni = accmap[entry[ufield]]
                    except: uni = None
                    if not uni:
                        if entry[ufield]: self.warnLog('Cannot map "%s" to %s Uniprot entry.' % (entry[ufield],string.join(self.list['PPISpec'],'/')))
                        continue
                    entry[ufield] = string.split(uni.info['Name'],'__')[-1]
                    entry['UniprotID'] = string.split(uni.info['Name'],'__')[0]
                    entry['Secondary'] = string.join(uni.obj['Sequence'].list['Secondary ID'],'|')
                xdb.index(ufield,force=True); xdb.index('Secondary',force=True); xdb.index('UniprotID',force=True)
                xdb.saveToFile(newxfile)
            if xdb:
                for field in self.list['MapFields'][0:]:
                    if field in xdb.fields():
                        xdb.index(field)
                        if field not in xref.list['MapFields']: xref.list['MapFields'].append(field)
                    else: self.warnLog('MapField "%s" not in XRef data.' % field)
            return xdb
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupSourceData(self,setup_ppi=True,setup_dmi=True):    ### Main class setup method.                                                      #V2.0
        '''Main class setup method.'''
        try:### ~ [0] Setup Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            sdb = db.addEmptyTable('Source',['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            rje.mkDir(self,self.getStr('SourcePath'),True)
            ### ~ [1] Compile PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dict['PPICompile']: return self.compilePPI()

            ### ~ [2] Download and Process Pairwise PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP SOURCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('SETUP SOURCE',line='=')
            supported_db = ['HINT']
            setup_ppi = setup_ppi or setup_dmi
            if self.getStr('PPISource') not in supported_db:
                ppifile = self.sourceDataFile('PPISource',expect=setup_ppi and self.getStrLC('PPISource'),ask=True)    # If this exists, use for everything!
                ppdb = self.parsePairwise(ppifile)

            ### ~ [3] Download and Process PPI Data from special DB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Restrict PPISource to be a pairwise file ONLY with Hub,Spoke,SpokeUni - maybe other keys
            #!# >> Ideally Hub,Spoke,HubTaxID,SpokeTaxID
            #!# >> This is HINT: move to special HINT method? #!# (Already made)
            #!# >> This section is just for old domain-based PPI processing until a new method is in place.
            if self.getStr('PPISource') in supported_db:   # Download sources are all found in sourceDataFile() for SOURCE.SPEC
                for spec in self.list['PPISpec']:
                    self.printLog('#SPEC','Looking for %s files.' % spec)
                    ## ~ [4a] Setup new species-specific attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.setStr({'%s.%s' % (self.getStr('PPISource'),spec): '%s.%s.ppi.tdt' % (self.getStr('PPISource'),spec),
                                 'PPI.%s' % spec: '%s.%s.pairwise.tdt' % (self.getStr('PPISource'),spec),
                                 'DomPPI.%s' % spec: '%s.%s.domppi.tdt' % (self.getStr('PPISource'),spec),
                                 'Uniprot.%s' % spec: 'uniprot.%s.dat' % spec})
                    ## ~ [4b] Download/Check UniProt DAT and fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    ppidat = self.sourceDataFile('Uniprot.%s' % spec,expect=setup_ppi,ask=False)  # If this exists, will use
                    if not ppidat and setup_ppi: self.warnLog('Something went wrong making/finding %s Uniprot file.' % spec, quitchoice=self.getBool('Integrity'))
                ## ~ [4c] Download/Check/Generate PPI Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    ppifile = self.sourceDataFile('PPI.%s' % spec,expect=False,ask=False)  # If this exists, will use
                    domfile = self.sourceDataFile('DomPPI.%s' % spec,expect=False,ask=False)  # If this exists, will use
                    needsource = (setup_ppi and not ppifile) or (setup_dmi and not domfile)
                    ppisource = self.sourceDataFile('%s.%s' % (self.getStr('PPISource'),spec),expect=needsource,ask=False)  # If this exists, will use
                ## ~ [4d] Download/Check/Generate Pairwise PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    ppuni = None
                    if ppifile:
                        if not ppisource: self.warnLog('%s exists but source %s not found.' % (self.getStr('PPI.%s' % spec),self.getStr('PPISource')),quitchoice=self.getBool('Integrity'))
                        ppdb = self.db().addTable(ppifile,mainkeys=['Hub','Spoke'],name='PPI.%s' % spec)
                    elif not ppifile and setup_ppi:  ## Generate Pairwise PPI file
                        ppifile = string.replace(ppisource,'ppi','pairwise')
                        ppi = rje_ppi.PPI(self.log,self.cmd_list)
                        ppi.loadPairwisePPI(ppisource)
                        ppdb = ppi.db('Edge')
                        if not ppdb: raise ValueError
                        ppdb.setStr({'Name':'PPI.%s' % spec})

                        # HINT using Gene_A and Gene_B fields with Gene Symbols.
                        for field in ['SpokeUni','HubUni']: ppdb.addField(field,after='Spoke',evalue='')
                        if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                        else:
                            ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=flybase,pfam,wormbase'])   #!# Make sure to add more as needed #!#
                            ppuni.setStr({'Name':ppidat})
                            ppuni.readUniProt()
                            sdb.data('Uniprot.%s' % spec)['Entries'] = ppuni.entryNum()
                            if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                        for entry in ppuni.entries():
                            gene = entry.gene().upper()
                            acc = entry.acc()
                            if 'FlyBase' in entry.dict['DB']:
                                for fbg in entry.dict['DB']['FlyBase']:
                                    for field in ['Id_A','Gene_A']:
                                        if field in ppdb.fields():
                                            for pentry in ppdb.indexEntries(field,fbg): pentry['HubUni'] = acc
                                    for field in ['Id_B','Gene_B']:
                                        if field in ppdb.fields():
                                            for pentry in ppdb.indexEntries(field,fbg): pentry['SpokeUni'] = acc
                            #if spec == 'CAEEL': self.deBug(entry.dict['DB'])
                            if 'WormBase' in entry.dict['DB']:
                                #self.deBug(entry.dict['DB']['WormBase'])
                                for fbg in entry.dict['DB']['WormBase']:
                                    for field in ['Id_A','Gene_A']:
                                        if field in ppdb.fields():
                                            for pentry in ppdb.indexEntries(field,fbg): pentry['HubUni'] = acc
                                    for field in ['Id_B','Gene_B']:
                                        if field in ppdb.fields():
                                            for pentry in ppdb.indexEntries(field,fbg): pentry['SpokeUni'] = acc
                            #if gene in ppdb.index('Hub'):
                            for pentry in ppdb.indexEntries('Hub',gene): pentry['HubUni'] = acc
                            for pentry in ppdb.indexEntries('Spoke',gene): pentry['SpokeUni'] = acc
                        #ppdb.dropEntriesDirect('Hub',[''])
                        #ppdb.dropEntriesDirect('Spoke',[''])
                        ppdb.saveToFile(ppifile); sdb.data('PPI.%s' % spec)['Status'] = 'Generated'
                        self.db().list['Tables'].append(ppdb)
                    else: self.printLog('#PPI','No %s PPI found.' % spec); continue
                    sdb.data('PPI.%s' % spec)['Entries'] = ppdb.entryNum()
                ## ~ [4d] Download/Check/Generate Domain PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if domfile:
                        domdb = self.db().addTable(domfile,mainkeys=['Pfam','Spoke'],name='DomPPI.%s' % spec)
                    elif not domfile and setup_dmi:  ## Generate Pairwise PPI file
                        domfile = string.replace(ppisource,'ppi','domppi')
                        domdb = self.db().addEmptyTable('DomPPI.%s' % spec,['Hub','Pfam','Spoke','SpokeUni'],keys=['Pfam','Spoke'])
                        if not ppuni:
                            if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                            else:
                                ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=flybase,pfam'])   #!# Make sure to add more as needed #!#
                                ppuni.setStr({'Name':ppidat})
                                ppuni.readUniProt()
                                sdb.data('Uniprot.%s' % spec)['Entries'] = ppuni.entryNum()
                                if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                        pfx = 0; epx = 0; ux = 0; utot = ppuni.entryNum()
                        for entry in ppuni.entries():
                            self.progLog('#PFAM','Adding DomPPI for %s UniProt entries: %.2f%%' % (rje.iStr(utot),ux/utot)); ux += 100.0
                            acc = entry.acc()
                            if 'Pfam' not in entry.dict['DB']: continue
                            epx += 1
                            for pfam in entry.dict['DB']['Pfam']:
                                pfx += 1
                                hub = string.split(entry.dict['DB']['Pfam'][pfam],';')[0]
                                for pentry in ppdb.indexEntries('HubUni',acc): domdb.addEntry({'Hub':hub,'Pfam':pfam,'Spoke':pentry['Spoke'],'SpokeUni':pentry['SpokeUni']},warn=False)
                        self.printLog('\r#PFAM','Added %s DomPPI for %s Pfam domains in %s of %s UniProt entries.' % (rje.iStr(domdb.entryNum()),rje.iStr(pfx),rje.iStr(epx),rje.iStr(utot)))
                        domdb.saveToFile(domfile); sdb.data('DomPPI.%s' % spec)['Status'] = 'Generated'
                    else: domdb = None
                    if domdb: sdb.data('DomPPI.%s' % spec)['Entries'] = domdb.entryNum()

            ### ~ [2] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SOURCE SUMMARY ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('SOURCE SUMMARY')
            for entry in sdb.entries():
                if entry['Entries'] == '': self.printLog('#SOURCE','%s (%s): %s' % (entry['Name'],entry['File'],entry['Status']))
                else: self.printLog('#SOURCE','%s %s (%s): %s' % (rje.iStr(entry['Entries']),entry['Name'],entry['File'],entry['Status']))
            #!# Add integrity/date check using sb entries #!#
            if not self.force() and 'Download' in self.db('Source').index('Status'):
                self.setBool({'Force':rje.yesNo('New source download detected: switch force=T?')})
            return True     # Setup successful
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem during PINGU.setupSourceData().'); return False  # Setup failed
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            setup_ppi = self.getBool('PPIFas')
            setup_dmi = self.getBool('DomPPI')
            if not self.setupSourceData(setup_ppi,setup_dmi): return False

            ### ~ [2] Setup Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('FasID'):
                fasid = string.split(rje.stripPath(self.getStrLC('PPISource')),'.')[0]
                self.debug(fasid)
                if setup_dmi: fasid += '-dom'
                self.setStr({'FasID':fasid})
            rje.mkDir(self,self.getStr('ResDir'))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### PPI Compilation/Parsing Methods                                                                         #
#########################################################################################################################
    def compilePPI(self):   ### Compilation of different PPI sources based on PINGU V3 code                         #V4.4
        '''Compilation of different PPI sources based on PINGU V3 code.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not self.dev(): return self.OLDcompilePPI()
            dbconvert = {'cpdb':'consensuspathdb'}
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ COMPILE PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            pdb = self.db().addEmptyTable('pairwise',['#','Hub','Spoke','HubUni','SpokeUni','HubTaxID','SpokeTaxID','Evidence','IType'],['#'])
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            xdb = self.xrefDB()
            if self.getBool('PPIDBReport'):
                self.db().addEmptyTable('report',['DB','Evidence','IType','PPI','Hub','Spoke'],['DB','Evidence','IType'])
                #self.db().addEmptyTable('report-evidence',['Evidence','PPI','Hubs','Spokes'],['Evidence'])
                #self.db().addEmptyTable('report-ppidb',['DB','PPI','Hubs','Spokes'],['DB'])
                #self.db().addEmptyTable('report-ppitype',['PPIType','PPI','Hubs','Spokes'],['DB'])
            repdb = self.db('report',add=False) #['DB','Evidence','IType','PPI','Hub','Spoke'],['DB','Evidence','IType']
            #if self.getBool('PPIDBReport'):
                #repdb = {}
                #repdb['METHOD'] = repdbev = self.db('report-evidence') #['Evidence','PPI','Hubs','Spokes'],['Evidence']
                #repdb['DB'] = repdbdb = self.db('report-ppidb')    #['DB','PPI','Hubs','Spokes'],['DB']
                #repdb['ITYPE'] = repdbtype = self.db('report-ppitype')#['PPIType','PPI','Hubs','Spokes'],['DB']
            if not xdb:
                self.printLog('#XREF','Cannot compile PPI data w/o XRef data. Check xrefdata=FILE setting.')
                return False
            xfields = self.list['MapFields'][0:]
            ## ~ [0a] Check Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for db in rje.sortKeys(self.dict['PPICompile']):
                dbfile = self.dict['PPICompile'][db]
                if not rje.exists(dbfile) and dbfile.upper() not in ['HPRD','HINT']:
                    raise IOError('%s not found!' % dbfile)

            ### ~ [1] ~ Parse PPI database files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in rje.sortKeys(self.dict['PPICompile']):
                dbfile = self.dict['PPICompile'][db]
                self.printLog('#PARSE','%s => Parse %s' % (db,dbfile))
                ## ~ [1a] Special database parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if dbfile.upper() in ['HPRD']: ddb = self.parseHPRD()
                elif dbfile.upper() in ['HINT']: ddb = self.parseHINT()
                ## ~ [1b] Pairwise PPI table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:
                    headers = rje.readDelimit(open(dbfile,'r').readline())
                    if 'Hub' in headers and 'Spoke' in headers:  # Pairwise PPI
                        ddb = self.db().addTable(dbfile,'#',name=db)
                        sdb.addEntry({'Name':db,'File':dbfile,'Status':'Found','Entries':ddb.entryNum()})
                        if 'Evidence' not in ddb.fields(): ddb.addField('Evidence',evalue='%s:unknown' % db)
                        if 'IType' not in ddb.fields(): ddb.addField('Evidence',evalue=db)
                        if 'HubUni' not in ddb.fields():
                            ddb.addField('HubUni')
                            for entry in ddb.entries():
                                #entry['HubUni'] = xref.xref(entry['Hub'],self.getStr('UniField'),mapfields=xfields,unique=True,usedict=True)
                                entry['HubUni'] = self.getUniXRef(entry['Hub'],xfields)
                        if 'SpokeUni' not in ddb.fields():
                            ddb.addField('SpokeUni')
                            for entry in ddb.entries():
                                #entry['SpokeUni'] = xref.xref(entry['Spoke'],self.getStr('UniField'),mapfields=xfields,unique=True,usedict=True)
                                entry['SpokeUni'] = self.getUniXRef(entry['Spoke'],xfields)
                ## ~ [1c] Parse Special files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    elif db.lower() == 'reactome': ddb = self.parseReactome(dbfile)
                    elif db.lower() in ['cpdb','consensuspathdb']: ddb = self.parseConsensusPathDB(dbfile)
                    #!# Remove reactome=FILE
                ## ~ [1d] Parse MITAB file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    else:
                        mcmd = ['mitab=%s' % dbfile,'dbsource=%s' % db]
                        mitab = rje_mitab.MITAB(self.log,self.cmd_list+mcmd)
                        mitab.obj['XRef'] = self.obj['XRef']
                        mitab.run(save=False)
                        ddb = mitab.db('pairwise')
                        sdb.addEntry({'Name':db,'File':dbfile,'Status':'Parsed','Entries':ddb.entryNum()})
                ## ~ [1e] PPIDBReport ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getBool('PPIDBReport'):
                    try:
                        drdb = self.db().copyTable(ddb,'%s-report' % db,replace=True,add=True)
                        drdb.addField('PPI',evalue=1)
                        drdb.keepFields(['#','Evidence','IType','PPI','Hub','Spoke'])
                        for dentry in drdb.entries():
                            #if db == 'hint': self.debug(dentry)
                            if db.lower() in ['cpdb','consensuspathdb']: dentry['Evidence'] = 'cpdb'
                            evidence = string.split(dentry['Evidence'],'|')
                            dentry['Evidence'] = []
                            for istr in evidence: dentry['Evidence'].append(string.split(istr,':',maxsplit=1)[-1])  # Strip db
                            dentry['Evidence'] = string.join(dentry['Evidence'],'|')
                            itypes = string.split(dentry['IType'],'|')
                            dentry['IType'] = []
                            for istr in itypes: dentry['IType'].append(string.split(istr,':',maxsplit=1)[-1])  # Strip db
                            dentry['IType'] = string.join(dentry['IType'],'|')
                        drdb.compress(['Evidence','IType'],rules={'Hub':'list','Spoke':'list','PPI':'sum'},joinchar='|')
                        drdb.dropField('#')
                        drdb.addField('DB',evalue=db)
                        for dentry in drdb.entries():
                            #if db == 'hint': self.bugPrint(dentry)
                            dentry['Hub'] = len(string.split(dentry['Hub'],'|'))
                            dentry['Spoke'] = len(string.split(dentry['Spoke'],'|'))
                        #['DB','Evidence','IType','PPI','Hub','Spoke'],['DB','Evidence','IType']
                        drdb.newKey(['DB','Evidence','IType'])
                        #self.db().mergeTables(repdb,drdb,overwrite=False,matchfields=False)
                        drx = drdb.entryNum(); repx = repdb.entryNum()
                        rje.combineDict(repdb.dict['Data'],drdb.dict['Data'])   # Should be no key conflicts
                        self.printLog('#MERGE','%s + %s %s -> %s report entries.' % (rje.iStr(repx),rje.iStr(drx),db,rje.iStr(repdb.entryNum())))
                    except: self.errorLog('PPIDBReport error')
                ## ~ [1f] Combine Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ddb:
                    #self.db().mergeTables(pdb,ddb,overwrite=False,matchfields=False)
                    ddx = ddb.entryNum(); ppx = pdb.entryNum();
                    if ppx:
                        pmax = max(pdb.datakeys())
                        for dentry in ddb.entries(): dentry['#'] += pmax
                    ddb.remakeKeys()
                    rje.combineDict(pdb.dict['Data'],ddb.dict['Data'])   # Should be no key conflicts now
                    self.printLog('#MERGE','%s + %s %s -> %s pairwise entries.' % (rje.iStr(ppx),rje.iStr(ddx),db,rje.iStr(pdb.entryNum())))
                else: raise ValueError('Parsing %s failed!' % db)

            ### ~ [2] ~ Compress and save Pairwise Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb.compress(['Hub','Spoke','HubTaxID','SpokeTaxID'],rules={'Evidence':'list','IType':'list'},joinchar='|')
            pdb.dropField('#')
            self.printLog('#PPI','%s unique pairwise PPI (Symmetry=%s)' % (rje.iStr(pdb.entryNum()),self.getBool('Symmetry')))
            pdb.saveToFile('%s.pairwise.tdt' % self.baseFile())
            ## ~ [2a] Log Report & Summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idb = None
            if self.v() > 0 or self.getBool('PPIDBReport'):
                idb = self.db().addEmptyTable('ireport',['Data','Value','PPI'],['Data','Value'])
                report = {'DB':{},'METHOD':{},'ITYPE':{}}
                rx = 0.0; rtot = pdb.entryNum()
                for entry in pdb.entries():
                    self.progLog('\r#REPORT','Parsing DB/Evidence report: %.2f%%' % (rx/rtot),rand=0.01); rx += 100.0
                    for dbmethod in string.split(entry['Evidence'],'|'):
                        (db,ev) = string.split(dbmethod,':',1)
                        if db in report['DB']: report['DB'][db] += 1
                        else: report['DB'][db] = 1
                        if ev in report['METHOD']: report['METHOD'][ev] += 1
                        else: report['METHOD'][ev] = 1
                    for itype in string.split(entry['IType'],'|'):
                        if itype in report['ITYPE']: report['ITYPE'][itype] += 1
                        else: report['ITYPE'][itype] = 1
                self.printLog('\r#REPORT','Parsed DB/Evidence report for %s PPI.' % rje.iStr(rtot))
                for ri in ['DB','METHOD','ITYPE']:
                    for rk in rje.sortKeys(report[ri]):
                        self.printLog('#%s' % ri,'%s: %s' % (rk,rje.iStr(report[ri][rk])))
                        idb.addEntry({'Data':ri,'Value':rk,'PPI':report[ri][rk]})
            ## ~ [2b] PPIDBReport ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('PPIDBReport'):
                idb.saveToFile('%s.ireport.tdt' % self.baseFile())
                repdb.saveToFile('%s.report.tdt' % self.baseFile())
                dbdb = self.db().addEmptyTable('dbreport',['DB'],['DB'])
                for db in rje.sortKeys(self.dict['PPICompile']): dbdb.addEntry({'DB':db})
                for db in rje.sortKeys(self.dict['PPICompile']): dbdb.addField(db,evalue=0)
                for db in ['Unique','Total']: dbdb.addField(db,evalue=0)
                for entry in pdb.entries():
                    pdblist = []
                    for dbmethod in string.split(entry['Evidence'],'|'):
                        evdb = string.split(dbmethod,':')[0]
                        if evdb not in self.dict['PPICompile'] and evdb in dbconvert: evdb = dbconvert[evdb]
                        if evdb not in self.dict['PPICompile']:
                            self.warnLog('Database %s not recognised' % evdb,warntype='db_missing',quitchoice=True,suppress=True,dev=False,screen=True)
                            continue
                        pdblist.append(evdb)
                    for db1 in pdblist:
                        for db2 in pdblist:
                            dbdb.data(db1)[db2] += 1
                        dbdb.data(db1)['Total'] += len(pdblist)
                        if len(pdblist) == 1: dbdb.data(db1)['Unique'] += 1
                dbdb.saveToFile('%s.dbreport.tdt' % self.baseFile())

            return pdb
        except: self.errorLog('%s.compilePPI error' % self); return False
#########################################################################################################################
    def getUniXRef(self,seqid,mapfields):   ### Returns first Uniprot ID for given sequence ID and mapfields
        '''Returns first Uniprot ID for given sequence ID and mapfields.'''
        xdb = self.xrefDB()
        xref = self.obj['XRef']
        unixref = xref.xref(seqid,self.getStr('UniField'),mapfields=mapfields,unique=False,usedict=True)
        if unixref: return unixref[0]
        else: return ''
#########################################################################################################################
    def filterPPI(self):   ### Compilation of different PPI sources based on PINGU V3 code                         #V4.4
        '''Compilation of different PPI sources based on PINGU V3 code.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb = self.db('pairwise')
            if not pdb: return None
            if 'HubTaxID' not in pdb.fields() and 'Evidence' not in pdb.fields() and 'IType' not in pdb.fields(): return pdb
            filt = {'TaxID':0,'BadDB':0,'GoodDB':0,'BadPPI':0,'GoodPPI':0}; fx = 0.0
            needtofilter = self.list['TaxID'] and not self.dict['PPICompile']
            for fkey in rje.sortKeys(filt)[:-1]:
                if self.list[fkey]: needtofilter = True
            if not needtofilter: return pdb
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ PPI FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            prex = pdb.entryNum()
            ### ~ [1] PPI Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pkey in pdb.dataKeys():
                self.progLog('\r#FILT','Filtering PPI: %.2f%%' % (fx/prex)); fx += 100.0
                entry = pdb.data(pkey)
                ## ~ [1a] Filter by TaxID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if 'HubTaxID' in pdb.fields():
                    if self.list['TaxID'] and not self.dict['PPICompile']:
                        if entry['HubTaxID'] not in self.list['TaxID'] and entry['SpokeTaxID'] not in self.list['TaxID']:
                            filt['TaxID'] += 1; pdb.dict['Data'].pop(pkey); continue
                ## ~ [1b] Filter by PPIType ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if 'IType' in pdb.fields():
                    itypes = string.split(entry['IType'],'|')
                    if self.list['GoodPPI'] and not rje.listIntersect(self.list['GoodPPI'],itypes):
                        filt['GoodPPI'] += 1; pdb.dict['Data'].pop(pkey); continue
                    if self.list['BadPPI']:
                        itypes = rje.listDifference(itypes,self.list['BadPPI'])
                        if itypes: entry['IType'] = string.join(itypes,'|')
                        else: filt['BadPPI'] += 1; pdb.dict['Data'].pop(pkey); continue
                ## ~ [1c] Filter by DB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if 'Evidence' in pdb.fields():
                    evidence = string.split(entry['Evidence'],'|')
                    goodev = []
                    for dbmethod in evidence[0:]:
                        db = string.split(dbmethod,':')[0]
                        if self.list['GoodDB'] and db not in self.list['GoodDB']: continue
                        if self.list['BadDB'] and db in self.list['BadDB']: continue
                        goodev.append(dbmethod)
                    if goodev: entry['Evidence'] = string.join(goodev,'|')
                    elif self.list['GoodDB']: filt['GoodDB'] += 1; pdb.dict['Data'].pop(pkey); continue
                    elif self.list['BadDB']: filt['BadDB'] += 1; pdb.dict['Data'].pop(pkey); continue
                    else: self.warnLog('Fiddlesticks!')
            self.printLog('\r#FILT','Filtered PPI: %s of %s PPI retained.' % (rje.iStr(pdb.entryNum()),rje.iStr(prex)),log=False)
            for fkey in rje.sortKeys(filt):
                if filt[fkey]: self.printLog('#FILT','%s PPI filtered on %s' % (rje.iStr(filt[fkey]),fkey))
            if prex != pdb.entryNum():
                self.printLog('#FILT','%s -> %s PPI after filtering.' % (rje.iStr(prex),rje.iStr(pdb.entryNum())))
                pdb.dict['Index'] = {}
                pdb.saveToFile('%s.filteredppi.tdt' % self.baseFile())
            return pdb
        except: self.errorLog('%s.filterPPI error' % self); return False
#########################################################################################################################
    def xHubPPI(self):  ### Converts pairwise table into numbers of interacting hubs for each spoke
        '''Converts pairwise table into numbers of interacting hubs for each spoke.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb = self.db('pairwise')
            if not pdb: return None
            ### ~ [1] ~ Calculate Numbers of Hubs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First compress to Hub:SpokeUni pairs. NOTE: SpokeUni can map to multiple spokes, so this is needed!
            pdb.compress(['Hub','SpokeUni'],rules={},default='sum',best=[],joinchar='|')
            # Now count the unique hubs per spoke
            for entry in pdb.entries(): entry['Hub'] = 1
            pdb.compress(['SpokeUni'],rules={},default='sum',best=[],joinchar='|')
            pdb.newKey(['Hub','SpokeUni'])  # Just for compatibility with later code
            #!# This assumes that different SpokeUni are not ultimately mapping on to the same Uniprot entry #!#
            ## ~ [1a] ~ Convert to cumulative counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hubx = pdb.indexKeys('Hub')
            hubx.reverse()
            for i in range(len(hubx)-1):
                x = hubx[i]
                ix = len(pdb.indexEntries('Hub',x))
                prex = len(pdb.indexDataList('Hub',hubx[i+1],'SpokeUni'))
                for entry in pdb.indexEntries('Hub',x): pdb.addEntry(rje.combineDict({'Hub':hubx[i+1]},entry,overwrite=False))
                postx = len(pdb.indexDataList('Hub',hubx[i+1],'SpokeUni'))
                self.printLog('#XHUB','Added %s x %dhub to %dhub: %s -> %s.' % (rje.iStr(ix),x,hubx[i+1],rje.iStr(prex),rje.iStr(postx)))
            return pdb
        except: self.errorLog('%s.xHubPPI error' % self); return False
#########################################################################################################################
    def parsePPI(self):     ### Parses PPI Database files into objects
        '''Parses PPI Database files into objects.'''
        try:### ~ [1] Setup HPRD-type databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            if self.getStrLC('HPRD'):
                hcmd = ['ppitype=','badtype=','hprdfas=F','alliso=F','domainfas=F','complexfas=F']    # Defaults
                hcmd = hcmd + self.cmd_list + ['genecards=F','hprdpath=%s' % self.getStr('HPRD')]
                self.dict['PPIDB']['HPRD'] = rje_hprd.HPRD(self.log,hcmd)
                sdb.addEntry({'Name':'HPRD','File':self.getStr('HPRD'),'Status':'Found'})
            ### ~ [2] Setup BioGRID-type databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in ['BioGrid','MINT','IntAct','Reactome','DIP','Domino']:
                if not self.getStrLC(db): continue
                elif not rje.checkForFile(self.getStr(db)):
                    self.errorLog('Cannot find %s: %s!' % (db,self.info[db]),printerror=False,quitchoice=True)
                dcmd = ['ppitype=','badtype=indirect_complex,neighbouring_reaction','symmetry=T','species=human','ppifas=F','ppitab=F','alltypes=F']  # Defaults
                dcmd = dcmd + self.cmd_list + ['ppifile=%s' % self.getStr(db),'seqin=None','dbsource=%s' % db.lower(),'debug=F']
                self.dict['PPIDB'][db] = rje_biogrid.BioGRID(self.log,dcmd)
                sdb.addEntry({'Name':db,'File':self.getStr(db),'Status':'Found'})
            ### ~ [3] Run objects to parse in data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for db in rje.sortKeys(self.dict['PPIDB']): self.dict['PPIDB'][db].parse()
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def parseReactome(self,dbfile): ### Parses Reactome Database file into objects
        '''Parses Reactome Database file into objects.'''
        try:### ~ [1] Setup database parser ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            db = 'Reactome'
            if not rje.checkForFile(dbfile): self.errorLog('Cannot find %s: %s!' % (db,self.info[db]),printerror=False,quitchoice=True); return False
            dcmd = ['ppitype=','badtype=indirect_complex,neighbouring_reaction','symmetry=T','species=human','ppifas=F','ppitab=F','alltypes=F']  # Defaults
            dcmd = dcmd + self.cmd_list + ['ppifile=%s' % dbfile,'seqin=None','dbsource=%s' % db.lower(),'debug=F']
            biogrid = rje_biogrid.BioGRID(self.log,dcmd)
            biogrid.parse()
            sdb.addEntry({'Name':db,'File':self.getStr(db),'Status':'Parsed'})
            return self.parsedDBtoPairwise(db,biogrid)
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def parseConsensusPathDB(self,dbfile): ### Parses ConsensusPathDB Database file into objects
        '''Parses ConsensusPathDB Database file into objects.'''
        try:### ~ [1] Setup database parser ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            db = 'ConsensusPathDB'
            if not rje.checkForFile(dbfile): self.errorLog('Cannot find %s: %s!' % (db,self.info[db]),printerror=False,quitchoice=True); return False
            pdb = self.db().addEmptyTable('pairwise',['#','Hub','Spoke','HubUni','SpokeUni','HubTaxID','SpokeTaxID','Evidence','IType'],['#'])
            xref = self.obj['XRef']
            ufield = self.getStr('UniField')
            try: taxid = {'human':'9606','mouse':'MOUSE','yeast':'YEAST'}[string.split(dbfile,'_')[1]]  # Add more nos.
            except: self.errorLog('Cannot extract TaxID for %s: %s!' % (db,self.info[db]),printerror=False,quitchoice=True); return False
            if db.lower() == 'hprd': xfields = ['HPRD']
            else: xfields = self.list['MapFields'][0:]
            ### ~ [2] Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #  ConsensusPathDB (version 30) list of human protein interactions
            #  source_databases     interaction_publications        interaction_participants        interaction_confidence
            #Biogrid,PhosphoPOINT,DIP,HPRD   22307056,7687743,17081983,12840032,862676,93330262,7867743,1322499,17192257,16464493,8626767,15466476,23602568,11971971,19847302        MK03_HUMAN      NA
            #BIND,HPRD       12834348,8267636,10779411,9733480,8355279,15952226,7045378,12437104,8663221,8939944,10966741,12039587,11342132,10346927,12593649        FA10_HUMAN      NA
            #BIND,HPRD       9003757 PROC_HUMAN      NA
            self.progLog('\r#PARSE','Parsing %s PPI' % db)
            dlines = open(dbfile,'r').readlines()
            dx = 0.0; dtot = len(dlines); ex = 0; fx = 0
            failed_id = []
            for dline in dlines:
                self.progLog('\r#PARSE','Parsing %s PPI: %.2f%%; %s PPI; %s failed.' % (db,dx/dtot,rje.iStr(ex),rje.iStr(fx))); dx += 100.0
                if dline[:1] == '#': continue
                try: (source,pub,prots,conf) = string.split(dline,'\t')
                except: continue
                evidence = string.join(['']+string.split(source,','),'|cpdb:')[1:]
                prots = string.split(string.replace(prots,'.',','),',')
                pairs = []
                if len(prots) == 1: itype = 'homodimer'; pairs = [prots * 2]
                elif len(prots) == 2: itype = 'binary'
                else: itype = 'cpdb complex'
                if not pairs: pairs = rje.listCombos([prots,prots])
                for pair in pairs:
                    #self.debug(pair)
                    [hub,spoke] = pair
                    hgene = xref.xref(hub,xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    #huni = xref.xref(hub,xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    huni = self.getUniXRef(hub,xfields)
                    if huni and not hgene: hgene = huni
                    if not hgene and hub not in failed_id: failed_id.append(hub)
                    sgene = xref.xref(spoke,xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    #suni = xref.xref(spoke,xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    suni = self.getUniXRef(spoke,xfields)
                    if suni and not sgene: sgene = suni
                    if not sgene and spoke not in failed_id: failed_id.append(spoke)
                    if not hgene or not sgene: fx += 1; continue
                    entry = {'#':ex,'Hub':hgene,'Spoke':sgene,'HubUni':huni,'SpokeUni':suni,
                             'HubTaxID':taxid,'SpokeTaxID':taxid,'Evidence':evidence,'IType':itype}
                    #if 'PCNA' in [hub,spoke,hgene,sgene]: self.bugPrint('%s >> %s' % (pair,entry))
                    #elif 'BRCA1' in [hub,spoke,hgene,sgene]: self.bugPrint('%s >> %s' % (pair,entry))
                    #elif self.dev(): continue
                    pdb.addEntry(entry); ex += 1
            self.printLog('\r#PARSE','Parsing %s PPI complete: %s PPI; %s failed.' % (db,rje.iStr(ex),rje.iStr(fx)))
            failed_id.sort()
            open('%s.%s.failed.id' % (self.basefile(),db),'w').write(string.join(failed_id,'\n'))
            self.printLog('#FAIL','%s failed IDs outout to %s.%s.failed.id.' % (rje.iLen(failed_id),self.basefile(),db))
            return pdb
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def parseHPRD(self):    ### Parses HPRD Database files into objects
        '''Parses HPRD Database files into objects.'''
        try:### ~ [1] Setup HPRD-type databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source') #,['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            if not self.getStrLC('HPRD'): return False
            hcmd = ['ppitype=','badtype=','hprdfas=F','alliso=F','domainfas=F','complexfas=F']    # Defaults
            hcmd = hcmd + self.cmd_list + ['genecards=F','hprdpath=%s' % self.getStr('HPRD')]
            hprd = rje_hprd.HPRD(self.log,hcmd)
            hprd.parse()
            sdb.addEntry({'Name':'HPRD','File':self.getStr('HPRD'),'Status':'Parsed'})
            return self.parsedDBtoPairwise('HPRD',hprd)
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def parsedDBtoPairwise(self,db,dbobj): ### Converts parsed HPRD/BioGRID object to Pairwise PPI table.
        '''Converts parsed HPRD/BioGRID object to Pairwise PPI table.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source')
            ppi = dbobj.dict['Evidence']
            pdb = self.db().addEmptyTable('pairwise',['#','Hub','Spoke','HubUni','SpokeUni','HubTaxID','SpokeTaxID','Evidence','IType'],['#'])
            xdb = self.xrefDB()
            xref = self.obj['XRef']
            ufield = self.getStr('UniField')
            if db.lower() == 'hprd': xfields = ['HPRD']
            else: xfields = self.list['MapFields'][0:]
            ### ~ [2] Parse/convert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0; ex = 0; fax = 0; fx = 0
            for hub in rje.sortKeys(ppi):
                hgene = xref.xref(hub,xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                #huni = xref.xref(hub,xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                huni = self.getUniXRef(hub,xfields)
                if huni and not hgene: hgene = huni
                #if 'PCNA' in [hub,hgene]: self.bugPrint('%s >> %s' % (hub,hgene))
                #elif 'BRCA1' in [hub,hgene]: self.bugPrint('%s >> %s' % (hub,hgene))
                if hgene == False: mx += len(ppi[hub]); fax += len(ppi[hub]); continue
                elif not hgene: mx += len(ppi[hub]); fx += len(ppi[hub]); continue
                for spoke in ppi[hub]:
                    self.progLog('\r#PPI','Converting %s: %s pairs; %s ppi; %s ambiguous; %s failed.' % (db,rje.iStr(mx),rje.iStr(ex),rje.iStr(fax),rje.iStr(fx)))
                    mx += 1
                    sgene = xref.xref(spoke,xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    #suni = xref.xref(spoke,xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True)
                    suni = self.getUniXRef(spoke,xfields)
                    if suni and not sgene: sgene = suni
                    if sgene == False: fax += 1; continue
                    elif not sgene: fx += 1; continue
                    entry = {'#':pdb.entryNum(),'Hub':hgene,'Spoke':sgene,'HubUni':huni,'SpokeUni':suni,
                                  'HubTaxID':'9606','SpokeTaxID':'9606',
                                  'Evidence':string.join(['']+ppi[hub][spoke],'|%s:' % db.lower())[1:],
                                  'IType':string.join(['']+ppi[hub][spoke],'|%s:' % db.lower())[1:]}
                    #if 'PCNA' in [spoke,sgene,hub,hgene]: self.bugPrint('%s,%s >> %s' % (hub,spoke,entry))
                    #elif 'BRCA1' in [spoke,sgene,hub,hgene]: self.bugPrint('%s,%s >> %s' % (hub,spoke,entry))
                    #elif self.dev(): continue
                    pdb.addEntry(entry); ex += 1
                    if self.getBool('Symmetry'):
                        pdb.addEntry({'#':pdb.entryNum(),'Hub':entry['Spoke'],'Spoke':entry['Hub'],
                                      'HubUni':entry['SpokeUni'],'SpokeUni':entry['HubUni'],
                                      'HubTaxID':entry['SpokeTaxID'],'SpokeTaxID':entry['HubTaxID'],
                                      'Evidence':entry['Evidence'],'IType':entry['IType']})
            self.printLog('\r#PPI','Converting %s complete: %s pairs; %s ppi; %s ambiguous; %s failed.' % (db,rje.iStr(mx),rje.iStr(ex),rje.iStr(fax),rje.iStr(fx)))
            sdb.data(db)['Entries'] = pdb.entryNum()
            return pdb
        except: self.log.errorLog('Pingu.parse() my arse!')
#########################################################################################################################
    def parseHINT(self):    ### Parsing HINT PPI file(s) into pairwise
        '''Parsing HINT PPI file(s) into pairwise.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('Source')
            taxid = {'HUMAN':'9606','YEAST':'559292','DROME':'7227','CAEEL':'6239'}
            for spec in rje.sortKeys(taxid):
                if taxid[spec] in self.list['TaxID'] and spec not in self.list['PPISpec']: self.list['PPISpec'].append(spec)
            self.printLog('#HINT','Parsing HINT for %s. (Use pairwise file input if already parsed.)' % string.join(self.list['PPISpec'],', '))
            pdb = self.db().addEmptyTable('hint',['#','Hub','Spoke','HubUni','SpokeUni','HubTaxID','SpokeTaxID','Evidence','IType'],['#'])
            xdb = self.xrefDB()
            xref = self.obj['XRef']
            xfields = self.list['MapFields'][0:]
            self.bugPrint(xfields)
            ### ~ [1] Parse HINT Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['PPISpec']:
                self.printLog('#SPEC','Looking for %s files.' % spec)
                ## ~ [1a] Setup new species-specific attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.setStr({'HINT.%s' % spec: '%s.%s.ppi.tdt' % (self.getStr('PPISource'),spec),
                             'Uniprot.%s' % spec: 'uniprot.%s.dat' % spec})
                ## ~ [2b] Download/Check UniProt DAT file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ppidat = self.sourceDataFile('Uniprot.%s' % spec,expect=True,ask=False)
                if not ppidat: raise IOError('Something went wrong making/finding %s UniProt file.' % spec)
                ppisource = self.sourceDataFile('HINT.%s' % spec,expect=True,ask=False)
                if not ppisource: raise IOError('Something went wrong making/finding %s HINT file.' % spec)
                ## ~ [2c] Parse HINT PPI Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ppuni = None
                ppi = rje_ppi.PPI(self.log,self.cmd_list)
                ppi.loadPairwisePPI(ppisource)
                ppdb = ppi.db('Edge')
                if not ppdb: raise ValueError
                ppdb.setStr({'Name':'PPI.%s' % spec})
                # HINT using Gene_A and Gene_B fields with Gene Symbols.
                for field in ['SpokeUni','HubUni']: ppdb.addField(field,after='Spoke',evalue='')
                if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                else:
                    ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=flybase,pfam,wormbase'])   #!# Make sure to add more as needed #!#
                    ppuni.setStr({'Name':ppidat})
                    ppuni.readUniProt()
                    sdb.data('Uniprot.%s' % spec)['Entries'] = ppuni.entryNum()
                    if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                for entry in ppuni.entries():
                    gene = entry.gene().upper()
                    acc = entry.acc()
                    if 'FlyBase' in entry.dict['DB']:
                        for fbg in entry.dict['DB']['FlyBase']:
                            for field in ['Id_A','Gene_A']:
                                if field in ppdb.fields():
                                    for pentry in ppdb.indexEntries(field,fbg): pentry['HubUni'] = acc
                            for field in ['Id_B','Gene_B']:
                                if field in ppdb.fields():
                                    for pentry in ppdb.indexEntries(field,fbg): pentry['SpokeUni'] = acc
                    #if spec == 'CAEEL': self.deBug(entry.dict['DB'])
                    if 'WormBase' in entry.dict['DB']:
                        #self.deBug(entry.dict['DB']['WormBase'])
                        for fbg in entry.dict['DB']['WormBase']:
                            for field in ['Id_A','Gene_A']:
                                if field in ppdb.fields():
                                    for pentry in ppdb.indexEntries(field,fbg): pentry['HubUni'] = acc
                            for field in ['Id_B','Gene_B']:
                                if field in ppdb.fields():
                                    for pentry in ppdb.indexEntries(field,fbg): pentry['SpokeUni'] = acc
                    #if gene in ppdb.index('Hub'):
                    for pentry in ppdb.indexEntries('Hub',gene): pentry['HubUni'] = acc
                    for pentry in ppdb.indexEntries('Spoke',gene): pentry['SpokeUni'] = acc
                ## ~ [2d] Update HINT pairwise table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sdb.data('HINT.%s' % spec)['Entries'] = ppdb.entryNum()
                for pentry in ppdb.entries():
                    entry = {'#':pdb.entryNum(),
                             'Hub':xref.xref(pentry['Hub'],xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True),
                             'Spoke':xref.xref(pentry['Spoke'],xfield=None,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True),
                             #'HubUni':xref.xref(pentry['HubUni'],xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True),
                             'HubUni':self.getUniXRef(pentry['HubUni'],xfields),
                             #'SpokeUni':xref.xref(pentry['SpokeUni'],xfield=ufield,mapfields=xfields,altfields=[],fullmap=False,unique=True,usedict=True),
                             'SpokeUni':self.getUniXRef(pentry['SpokeUni'],xfields),
                             'HubTaxID':taxid[spec],'SpokeTaxID':taxid[spec],'IType':'hint','Evidence':'hint:hint'}
                    try:
                        entry['Evidence'] = []
                        for hintev in string.split(pentry['pmid:method:quality'],'|'):
                            entry['Evidence'].append('hint:%s' % string.split(hintev,':')[1])
                            entry['Evidence'].append('hint:%s' % string.split(hintev,':')[2])
                        entry['Evidence'] = string.join(rje.sortUnique(entry['Evidence']),'|')
                    except: entry['Evidence'] = 'hint:unparsed'
                    #if 'PCNA' in [pentry['Hub'],pentry['Spoke'],entry['Hub'],entry['Spoke']]: self.bugPrint('%s >> %s' % (pentry,entry))
                    #elif 'BRCA1' in [pentry['Hub'],pentry['Spoke'],entry['Hub'],entry['Spoke']]: self.bugPrint('%s >> %s' % (pentry,entry))
                    #elif self.dev(): continue
                    if 'None' in [entry['Hub'],entry['Spoke']]: continue
                    if not entry['Hub'] or not entry['Spoke']: continue
                    pdb.addEntry(entry)
                    if self.getBool('Symmetry'):
                        pdb.addEntry({'#':pdb.entryNum(),'Hub':entry['Spoke'],'Spoke':entry['Hub'],
                                      'HubUni':entry['SpokeUni'],'SpokeUni':entry['HubUni'],
                                      'HubTaxID':entry['SpokeTaxID'],'SpokeTaxID':entry['HubTaxID'],
                                      'Evidence':entry['Evidence'],'IType':entry['IType']})
            return pdb
        except: self.log.errorLog('Pingu.parse(HINT) my arse!'); return None
#########################################################################################################################
    def parsePairwise(self,dbfile): ### Parse PINGU pairwise file into 'pairwise' Database table.
        '''
        Parse PINGU pairwise file into 'pairwise' Database table.
        >> dbfile:str = File from which to parse PPI data.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ PARSE PAIRWISE PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            sdb = self.db('Source')
            self.list['HubList'].sort()
            delimit = rje.delimitFromExt(filename=dbfile)
            headers = rje.readDelimit(open(dbfile,'r').readline(),delimit)
            hubfield = self.getStr('HubField')
            if hubfield not in headers and hubfield.lower() in ['hubuni','uniprot','accnum']:
                self.printLog('#FIELD','Hub field "%s" -> "HubUni"' % hubfield)
                hubfield = 'HubUni'
            if hubfield not in headers: self.warnLog('HubField "%s" not found in headers' % hubfield)
            spokefield = self.getStr('SpokeField')
            if spokefield not in headers and spokefield.lower() in ['spokeuni','uniprot','accnum']:
                self.printLog('#FIELD','Spoke field "%s" -> "SpokeUni"' % spokefield)
                spokefield = 'SpokeUni'
            if spokefield not in headers: self.warnLog('SpokeField "%s" not found in headers' % spokefield)
            ### ~ [1] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add expandppi=X looping to increase hublist. (Might be easier to read whole in and then remove)
            if self.list['HubList'] and headers[0] == hubfield:
                ex = 0
                #self.file['pairwise'] = open(dbfile,'r')
                self.setStr({'pairwise':dbfile})
                ppdb = self.db().addEmptyTable('pairwise',headers,['Hub','Spoke'])
                self.open('pairwise')
                dline = self.findlines('pairwise',self.list['HubList'],asdict=False,wrap=False,chomp=True,warnings=True,nextonly=True)
                while dline:
                    self.progLog('\r#PPI','Reading Pairwise PPI for %s hubs (hubonly=%s): %s PPI' % (rje.iLen(self.list['HubList']),self.getBool('HubOnly'),rje.iStr(ex)))
                    pentry = rje.listToDict(rje.readDelimit(dline,delimit),headers)
                    if not self.getBool('HubOnly') or pentry[spokefield] in self.list['HubList']: ppdb.addEntry(pentry); ex += 1
                    dline = self.findlines('pairwise',self.list['HubList'],asdict=False,wrap=False,chomp=True,warnings=True,nextonly=True)
                self.close('pairwise')
                self.printLog('\r#PPI','Read in %s pairwise PPI for %s hubs (hubonly=%s) from %s.' % (rje.iStr(ppdb.entryNum()),rje.iLen(self.list['HubList']),self.getBool('HubOnly'),dbfile))
                ppdb.saveToFile('%s.pairwise.tdt' % self.baseFile())
            else:
                #?# Symmetry cause clashes and don't actually want redundancy?
                if 'InteractionID' in headers:  #APID data
                    ppdb = self.db().addTable(dbfile,mainkeys=['Hub','InteractionID'],name='pairwise')
                else:
                    ppdb = self.db().addTable(dbfile,mainkeys=['Hub','Spoke'],name='pairwise')
                self.printLog('\r#PPI','Reading %s pairwise PPI from %s.' % (rje.iStr(ppdb.entryNum()),dbfile))
                if self.list['HubList']:
                    ppdb.dropEntriesDirect(hubfield,self.list['HubList'],inverse=True)
                    if self.getBool('HubOnly'):
                        ppdb.dropEntriesDirect(spokefield,self.list['HubList'],inverse=True)
            ### ~ [2] Entry data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.data('PPISource')['Entries'] = ppdb.entryNum()
            return ppdb
        except: self.errorLog('Pingu.parse() my arse!'); return None
#########################################################################################################################
    ### <4> ### Class Output Methods                                                                                    #
#########################################################################################################################
    def ppiFas(self,pdb=None):   ### Outputs PPI data to Fasta Files
        '''
        Outputs PPI data to Fasta Files.
        >> pdb:Table with "Hub" and "SpokeUni" fields.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # PPI data should be stored in PPI.SPEC or DomPPI.SPEC tables.
            force = self.force()
            db = self.obj['DB']
            if not pdb: pdb = self.db('pairwise',add=False)
            self.debug(pdb)
            ppitype = 'PPI'; hubtype = 'Hub'
            if self.getBool('DomPPI'): ppitype = 'DomPPI'; hubtype = 'Pfam'
            ### ~ [1] ~ Generate PPI Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['PPISpec']:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ %s PPI DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % spec)
                if not pdb: pdb = self.db('%s.%s' % (ppitype,spec))
                ppisource = rje.baseFile(self.getStr('PPISource'),strip_path=True)
                datadir = rje.makePath('%s%s.%s.%s/' % (self.getStr('ResDir'),ppisource,ppitype,spec))
                self.debug(datadir)

                #if self.db('PPISource',add=False):
                #    pdb = self.db('PPISource')
                #    datadir = rje.makePath('%s%s.%s/' % (self.getStr('ResDir'),rje.baseFile(pdb.getStr('Source'),True),spec))
                #elif self.getBool('DomPPI'):
                #    pdb = self.db('DomPPI.%s' % spec)
                #    datadir = rje.makePath('%sDomPPI.%s/' % (self.getStr('ResDir'),spec))
                #else:
                #    pdb = self.db('PPI.%s' % spec)
                #    datadir = rje.makePath('%sPPI.%s/' % (self.getStr('ResDir'),spec))
                if os.path.exists(datadir) and force: rje.deleteDir(self,datadir)
                rje.mkDir(self,datadir,True)
                self.setStr({'Uniprot.%s' % spec: 'uniprot.%s.dat' % spec})
                ppidat = self.sourceDataFile('Uniprot.%s' % spec,expect=True,ask=False)
                if not ppidat: raise IOError('Something went wrong making/finding %s UniProt file.' % spec)
                #ppidat = self.getStr('Uniprot.%s' % spec)
                self.debug(ppidat)
                if spec in self.dict['UniSpec']:
                    ppuni = self.dict['UniSpec'][spec]
                    unimap = ppuni.accDictFromEntries()     # Map all accession numbers onto species UniProt entries
                else: unimap = {}
                ## ~ [1a] Take each hub in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.debug(self.getBool('CombinePPI'))
                if self.getBool('CombinePPI'):
                    #basename = rje.baseFile(self.baseFile(),strip_path=True)
                    ppiseqfile = '%s%s.%s.fas' % (datadir,spec,self.getStr('FasID'))
                    ppiseqfile = rje.fileSafeString(ppiseqfile,'-')
                    if not self.getBool('CombinePPI') and rje.exists(ppiseqfile) and not force:
                        self.printLog('#SKIP','%s file %s found (force=F).' % (spec,ppiseqfile))
                        continue
                seqfilex = 0; unihubx = 0; minx = 0; sfail = []; hfail = []
                for hub in rje.sortKeys(pdb.index(hubtype)):  # This is now a hub protein: find in ppi
                    self.progLog('\r#PPI','%s PPI files made from %s hubs. (%s < %d MinPPI)' % (rje.iStr(seqfilex),rje.iStr(unihubx),rje.iStr(minx),self.getInt('MinPPI')))
                    unihubx += 1
                    if not self.getBool('CombinePPI'):
                        if self.getBool('XHubPPI'): ppiseqfile = '%s%s.%s.%shub.fas' % (datadir,spec,self.getStr('FasID'),hub)
                        else: ppiseqfile = '%s%s.%s.fas' % (datadir,hub,self.getStr('FasID'))
                        ppiseqfile = rje.fileSafeString(ppiseqfile,'-')
                        if rje.exists(ppiseqfile) and not force: seqfilex += 1; continue
                    if not unimap:    # This is done here to save time if all files are made
                        #ppuni = rje_uniprot.UniProt(self.log,self.cmd_list)
                        ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=%s' % self.getStr('SourceData')])
                        self.debug(self.getStr('SourceData'))
                        ppuni.setStr({'Name':ppidat})
                        self.debug(ppidat)
                        ppuni.readUniProt()
                        if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                        unimap = ppuni.accDictFromEntries()     # Map all accession numbers onto species UniProt entries
                        self.debug(len(unimap))
                    ppiseqs = []
                    if self.getBool('CombinePPI'):
                        spokelist = pdb.indexKeys('SpokeUni')
                        hub = 'combined'
                    else: spokelist = pdb.indexDataList(hubtype,hub,'SpokeUni')
                    while '' in spokelist: spokelist.remove('')
                    for spoke in spokelist:
                        spacc = spoke.split('-')[0]     # PPI can have isoforms?
                        try: seq = (unimap[spacc].obj['Sequence'],spoke)   # SeqAcc can be isoform
                        except:
                            if self.db('pairwise',add=False):
                                if spoke not in sfail: sfail.append(spoke)
                            else: self.printLog('#ERR','Failed to get sequence for %s.' % (spoke))
                            continue
                        if seq not in ppiseqs: ppiseqs.append(seq)
                    if ppiseqs:
                        if len(ppiseqs) < self.getInt('MinPPI'): minx += 1; continue
                        ppifas = []; redseq = []
                        for seq in ppiseqs:
                            seqfas = seq[0].fasta(seq[1])
                            if seqfas not in ppifas: ppifas.append(seqfas)
                            else:
                                if seq[0] not in redseq: redseq.append(seq[0])
                                self.warnLog('%s already detected in output. Secondary AccNum in pairwise PPI?' % seq[0].shortName())
                        if len(ppifas) < self.getInt('MinPPI'): minx += 1; continue
                        ppifas.sort()
                        PPISEQ = open(ppiseqfile,'w')
                        PPISEQ.write(string.join(ppifas,''))
                        PPISEQ.close()
                        self.printLog('\r#SEQ','%s sequences output to %s' % (len(ppifas),ppiseqfile),screen=self.v()>1)
                        for rseq in redseq:
                            redspoke = []
                            for seq in ppiseqs:
                                if seq[0] == rseq: redspoke.append(seq[1])
                            redspoke.sort()
                            self.warnLog('Check %s for multiple mapping to %s: %s.' % (ppisource,rseq.shortName(),string.join(redspoke,'; ')))
                        seqfilex += 1
                    else:
                        if self.db('pairwise',add=False):
                            if hub not in hfail: hfail.append(hub)
                        else: self.printLog('\r#PPI','No sequences to output for %s PPI. (No mapped spokes.)' % (hub))
                    if self.getBool('CombinePPI'): unihubx = len(pdb.index('Hub')); break
                self.printLog('\r#PPI','%s PPI files made from %s hubs. (%s < %d MinPPI)' % (rje.iStr(seqfilex),rje.iStr(unihubx),rje.iStr(minx),self.getInt('MinPPI')))
                if self.db('pairwise'):
                    self.printLog('#PPI','%s hubs and %s spokes missing %s PPI data in %s' % (rje.iLen(hfail),rje.iLen(sfail),spec,self.db('pairwise').getStr('Source')))
            return
        except: self.errorLog('%s.ppiFas error' % self)
#########################################################################################################################
    def queryPPI(self):   ### Outputs Query PPI data to Fasta Files
        '''Outputs Query PPI data to Fasta Files'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # PPI data should be in PPI.SPEC/ or DomPPI.SPEC/ directories.
            force = self.force()
            db = self.obj['DB']
            pdb = self.db('pairwise',add=False)
            spec = self.list['PPISpec'][0]
            self.printLog('#QSPEC','## ~~~~~~~~~~~~~~~~~~~~~ %s QueryPPI ~~~~~~~~~~~~~~~~~~~ ##' % spec,log=False)
            if len(self.list['PPISpec']) > 1: self.printLog('#QPPI','Can only process one PPISpec at a time for QueryPPI: will use %s.' % spec)
            ppitype = 'PPI'; hubtype = 'Hub'
            if self.getBool('DomPPI'): ppitype = 'DomPPI'; hubtype = 'Pfam'
            ppisource = rje.baseFile(self.getStr('PPISource'))
            if not pdb: pdb = self.db('%s.%s' % (ppitype,spec))
            if not pdb: self.debug('%s.%s' % (ppitype,spec));
            #datadir = rje.makePath('%s%s.%s/' % (self.getStr('ResDir'),ppitype,spec))
            ## ~ [1a] Load Query PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qdb = db.addTable(self.getStr('QueryPPI'),mainkeys=['Query','Hub'],name='QueryPPI',expect=True)
            useqseq = self.getStrLC('QuerySeq')
            if not self.getStrLC('QuerySeq'): self.setStr({'QuerySeq':rje.baseFile(self.getStr('QueryPPI'))+'.fas'})
            useqseq = useqseq or rje.exists(self.getStr('QuerySeq'))
            if useqseq:
                qseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('QuerySeq'),'autoload=T'])
                qseqdict = qseqlist.seqNameDic()
            ## ~ [1b] Check for sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for qry in qdb.index('Query'):
                    if qry not in qseqdict: self.warnLog('Cannot find Query "%s" in %s' % (qry,self.getStr('QuerySeq')),warntype='Missing Query',suppress=self.dev())
            else: self.printLog('#QSEQ','No query sequence file: will look for queries within PPI datasets.')
            ## ~ [1c] Make output directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppisource = string.split(self.getStr('PPISource'),'.')[0]
            datadir = rje.makePath('%s%s.%s.%s.%s/' % (self.getStr('ResDir'),ppisource,ppitype,spec,rje.baseFile(self.getStr('QueryPPI'),strip_path=True)))
            if os.path.exists(datadir) and force: rje.deleteDir(self,datadir)
            rje.mkDir(self,datadir,True)

            ### ~ [2] ~ Generate PPI Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hubx = 0; qppix = 0; seqfilex = 0; minx = 0
            #self.debug(qdb.data())
            for qhub in rje.sortKeys(qdb.index('Hub'))[0:]:
                #self.progLog('\r#QPPI','%s Query PPI Datasets for %s hubs, %s Queries. (%s < %d MinPPI)' % (rje.iStr(seqfilex),rje.iStr(hubx),rje.iStr(qppix),rje.iStr(minx),self.getInt('MinPPI')))
                if qhub in pdb.index('Hub'): qhub = hub
                else: hub = self.obj['XRef'].xref(qhub,xfield=None,fullmap=False,unique=True)
                if not hub: self.warnLog('Failed to map Hub %s' % qhub); continue
                ## ~ [2a] Find Hub PPI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                hubppifile = None
                #for spec in self.list['PPISpec']:
                ppiseqfile = '%s%s.%s.%s/%s.%s.fas' % (self.getStr('ResDir'),ppisource,ppitype,spec,hub.upper(),self.getStr('FasID'))
                if os.path.exists(ppiseqfile): hubppifile = ppiseqfile#; break
                if not hubppifile:
                    #self.printLog('#SEQ','Unable to find %s%s.*/%s.%s.fas' % (self.getStr('ResDir'),ppitype,hub,self.getStr('FasID')))
                    self.printLog('#SEQ','Unable to find %s!' % (hubppifile))
                    if self.dev(): self.deBug('...')
                    continue
                hseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqindex=F'])
                hseqlist.loadSeq(hubppifile)
                hseqdict = hseqlist.makeSeqNameDic(keytype='max')
                hubx += 1
                ## ~ [2b] Generate a file per query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                hubqueries = qdb.indexDataList('Hub',qhub,'Query')
                for qry in hubqueries:
                    qentry = qdb.data(keylist=[qry,qhub])
                    #self.debug(qentry)
                    qid = None
                    if useqseq and qry not in qseqdict: continue
                    elif not useqseq and qry not in hseqdict:
                        qid = self.obj['XRef'].xref(qry,xfield='Uniprot',fullmap=False,unique=True)
                        if qid not in hseqdict:
                            self.warnLog('Cannot find Query "%s" (%s) in %s' % (qry,qid,hubppifile),warntype='Missing Query',suppress=self.dev())
                            continue
                    qppix += 1
                    if 'QueryID' in qentry and qentry['QueryID']: ppiseqfile = '%s%s.%s.%s.fas' % (datadir,hub.upper(),rje.fileSafeString(qentry['QueryID']),self.getStr('FasID'))
                    else: ppiseqfile = '%s%s.%s.%s.fas' % (datadir,hub.upper(),rje.fileSafeString(qry),self.getStr('FasID'))
                    if useqseq:
                        ppiseqs = [qseqlist.getSeq(qseqdict[qry])]
                        if self.getBool('AllQuery'):
                            for allqry in hubqueries:
                                if allqry not in [qry] + hseqdict.keys(): ppiseqs.append(qseqlist.getSeq(qseqdict[allqry]))
                    elif qid: ppiseqs = [hseqlist.getSeq(hseqdict[qid])]
                    else: ppiseqs = [hseqlist.getSeq(hseqdict[qry])]
                    pepmatch = True
                    if 'QRegion' in qdb.fields():
                        qregion = qentry['QRegion']
                        qsequence = ppiseqs[-1][1].lower()
                        try: pepseq = string.split(qentry['Peptide'],'|')
                        except: pepseq = []
                        for qregion in string.split(qregion,'|')[0:]:
                            if not qregion: continue
                            qregion = string.split(qregion,',')
                            qregion[0] = int(qregion[0]) - 1
                            qregion[1] = int(qregion[1])
                            #self.debug(qsequence)
                            regpeptide = qsequence[qregion[0]:qregion[1]].upper()  #!# Modify to split on | and match by index to QRegion
                            if pepseq:
                                peptide = pepseq.pop(0).upper()
                                if regpeptide != peptide:
                                    pepmatch = False
                                    self.warnLog('Peptide sequence "%s" does not match QRegion %s..%s = %s.' % (peptide,qregion[0],qregion[1],regpeptide))
                            qsequence = qsequence[:qregion[0]] + regpeptide + qsequence[qregion[1]:]
                            ppiseqs[-1] = (ppiseqs[-1][0],qsequence)
                            #self.debug(ppiseqs[-1][1])
                    for seq in hseqlist.seqs()[0:]:
                        (name,sequence) = hseqlist.getSeq(seq)
                        if ppiseqs[0][0] != name: ppiseqs.append((name,sequence))
                    if ppiseqs:
                        if len(ppiseqs) < self.getInt('MinPPI'): minx += 1; continue
                        PPISEQ = open(ppiseqfile,'w')
                        for seq in ppiseqs: PPISEQ.write('>%s\n%s\n' % (seq[0],seq[1]))
                        PPISEQ.close()
                        self.printLog('\r#SEQ','%s sequences output to %s' % (len(ppiseqs),ppiseqfile),screen=self.v()>0)
                        seqfilex += 1
            self.printLog('\r#QPPI','%s Query PPI Datasets for %s hubs, %s Queries. (%s < %d MinPPI)' % (rje.iStr(seqfilex),rje.iStr(hubx),rje.iStr(qppix),rje.iStr(minx),self.getInt('MinPPI')))
        except: self.errorLog('%s.queryPPI error' % self)
#########################################################################################################################
    ### <5> ### REST server output methods. (Should be replaced in mature Classes)                                      #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT:

        pairwise = main table of identified PPI for given `hublist=LIST` proteins. [tdt]
        spokes = non-redundant list of "spoke" genes that interact with hubs [list]
        uniprot = non-redundant list of uniprot accession numbers for proteins that interact with hubs [list]
        sif-gene = simple interaction file (SIF) format of gene identifiers for PPI [sif]
        sif-uni = simple interaction file (SIF) format of uniprot identifiers for PPI [sif]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in ['spokes','uniprot','sif-gene','sif-uni']: self.dict['Output'][outfmt] = 'No output generated.'
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):  ### Defines REST Output order
        '''Defines REST Output order.'''
        rfmt = ['pairwise','spokes','uniprot','sif-gene','sif-uni'] + self.list['RestDB']
        for fmt in rje.sortKeys(self.dict['Output']):
            if fmt not in rfmt: rfmt.append(fmt)
        return rfmt
#########################################################################################################################
    def generateRestOutput(self):   ### Generate REST output dictionary
        '''Generate REST output dictionary.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Rest'): return
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ REST OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            if not self.list['HubList']: self.dict['Output']['pairwise'] = 'HubList needed for REST output.'; return
            ### ~ [1] ~ Generate Rest Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # pairwise output links through self.db('pairwise')
            pdb = self.db('pairwise')
            if not pdb or not pdb.entryNum(): self.dict['Output']['pairwise'] = 'No PPI parsed/retained.'
            else:
                self.dict['Output']['spokes'] = string.join(pdb.indexKeys('Spoke'),'\n')
                uacc = pdb.indexKeys('SpokeUni')
                if '' in uacc: uacc.remove('')
                #afile = '%s.uniprot' % self.baseFile()
                #self.dict['Output']['uniprot'] = afile
                #open(afile,'w').write(string.join(uacc,'\n'))
                self.dict['Output']['uniprot'] = string.join(uacc,'\n')
            ## ~ [1a] ~ SIF Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gfile = '%s.gene.sif' % self.baseFile(); GFILE = open(gfile,'w')
            self.dict['Output']['sif-gene'] = gfile
            ufile = '%s.uniprot.sif' % self.baseFile(); UFILE = open(ufile,'w')
            self.dict['Output']['sif-uni'] = ufile
            for ekey in pdb.dataKeys():
                entry = pdb.data(ekey)
                if entry['Hub'] and entry['Spoke']: GFILE.write('%s pp %s\n' % (entry['Hub'],entry['Spoke']))
                if entry['HubUni'] and entry['SpokeUni']: UFILE.write('%s pp %s\n' % (entry['HubUni'],entry['SpokeUni']))
            GFILE.close()
            UFILE.close()
            ## ~ [1b] ~ Database PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['RestDB'] = []
            if self.dev():
                phead = pdb.fields()
                pdb.index('Evidence',splichar='|')
                for evidence in pdb.indexKeys('Evidence'):
                    (db,etype) = string.split(evidence,':',1)
                    if not self.db(db): dbdb = self.db().addEmptyTable(db,phead,['Hub','Spoke']); self.list['RestDB'].append(db)
                    else: dbdb = self.db(db)
                    for entry in pdb.indexEntries('Evidence',evidence): dbdb.addEntry(entry)
                self.list['RestDB'].sort()
        except: self.errorLog('PINGU.generateRestOutput() error')
#########################################################################################################################
### End of SECTION II: PINGU Class                                                                                      #
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
    try: PINGU(mainlog,cmd_list).run()

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
