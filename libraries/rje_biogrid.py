#!/usr/local/bin/python

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
Module:       rje_biogrid
Description:  BioGRID Database processing module
Version:      1.6
Last Edit:    07/05/10
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed primarily for parsing the plain text ORGANISM downloads from the BioGRID database. These
    have names in the form: BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.27.tab.txt.

    BioGRID tables contain useful information that can be used for cross-referencing to other sources, namely the protein
    names and gene symbols/aliases. The latter will be added to the dict['Mapping'] links dictionary of the BioGRID
    object, linking each symbol to the primary protein ID. These protein IDs will be used for storing the PPI data (in
    dict['PPI']) and extracting gene data from external sequence databases. These sequence databases need to be provided
    separately. This will be read in and added to the dict['Protein'] which will also store gene symbol data etc.

    The selection of sequence files might turn out to be quite tricky, as different species have very different protein
    identifiers used. I will add a list of recommended sequence sources as I find them:
    * Yeast = EnsLoci treatment of the EnsEMBL yeast genome

    BioGRID contains data for a number of experimental types. Those of interest can be specified with the ppitype=LIST
    option. Choices include: Affinity Capture-MS; Affinity Capture-Western; Biochemical Activity; Co-crystal Structure;
    Co-fractionation; Co-purification; Dosage Lethality; Dosage Rescue; Far Western; FRET; Phenotypic Enhancement;
    Phenotypic Suppression; Protein-peptide; Reconstituted Complex; Synthetic Growth Defect; Synthetic Lethality;
    Synthetic Rescue; Two-hybrid;

    IntAct has the following:  anti bait coip | pull down | two hybrid pooling | two hybrid | tap | x-ray diffraction |
    anti tag coip | fluorescence imaging | cosedimentation | elisa | protein kinase assay | coip | biochemical |
    antibody array | confocal microscopy | beta galactosidase | imaging techniques | two hybrid array |
    molecular sieving | ion exchange chrom | affinity chrom | protein array | enzymatic study | inferred by curator |
    far western blotting | spr | fps | phosphatase assay | fret | dhfr reconstruction | bn-page | peptide array | nmr |
    facs | affinity techniques | crosslink | itc | one hybrid | fluorescence | solution sedimentati | ch-ip | emsa |
    complementation | density sedimentatio | comig non denat gel | filter binding | chromatography | ub reconstruction |
    reverse phase chrom | emsa supershift | electron microscopy | protein crosslink | competition binding | mappit |
    gallex | gtpase assay | in gel kinase assay | spa | biophysical | radiolabeled methyl | experimental interac |
    fluorescence spectr | cd | bret | protein tri hybrid | transcription compl | deacetylase assay | footprinting |
    yeast display | saturation binding | protease assay | lambda phage | light scattering | htrf | fcs | toxcat |
    phage display | t7 phage | kinase htrf | methyltransferase as

    MINT txt downloads can also be parsed. Experiment types for MINT include:  affinity chromatography technologies |
    affinity technologies | anti bait coimmunoprecipitation | anti tag coimmunoprecipitation |
    beta galactosidase complementation | beta lactamase complementation | biochemical |
    bioluminescence resonance energy transfer | biophysical | chromatography technologies | circular dichroism |
    classical fluorescence spectroscopy | coimmunoprecipitation | colocalization by fluorescent probes cloning |
    colocalization by immunostaining | colocalization/visualisation technologies | competition binding | copurification |
    cosedimentation | cosedimentation in solution | cosedimentation through density gradients | cross-linking studies |
    electron microscopy | enzymatic studies | enzyme linked immunosorbent assay | experimental interaction detection |
    far western blotting | filter binding | fluorescence-activated cell sorting | fluorescence microscopy |
    fluorescence polarization spectroscopy | fluorescence technologies | fluorescent resonance energy transfer |
    gst pull down | his pull down | imaging techniques | isothermal titration calorimetry | lambda phage display |
    mass spectrometry studies of complexes | molecular sieving | nuclear magnetic resonance | peptide array |
    phage display | protease assay | protein array | protein complementation assay | protein kinase assay | pull down |
    saturation binding | surface plasmon resonance | t7 phage display | two hybrid | two hybrid array |
    two hybrid fragment pooling approach | two hybrid pooling approach | ubiquitin reconstruction | unknown |
    x-ray crystallography

    Reactome interactions are restricted to those of the "reaction" type. There are also "neighbouring_reaction" and
    "direct_complex" and "indirect_complex"

    DIP interactions are restricted to those with two uniprotkb IDs. DIP has similar annotation to MINT, with MI nos.

    Domino  interactions are restricted to those with two uniprotkb IDs. Has similar annotation to MINT, with MI nos.

Commandline:
    ### ~ BioGRID parsing and PPI Dataset Generation Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ppifile=FILE    : PPI database flat file [None]
    seqin=FILE      : Sequence file containing protein sequences with appropriate Accession Numbers/IDs [None]
    genecards=FILE  : File of links between IDs. For human, should have HGNC and EnsLoci columns. [None]
    ppitype=LIST    : List of acceptable interaction types to parse out []
    badtype=LIST    : List of bad interaction types, to exclude [indirect_complex,neighbouring_reaction]
    symmetry=T/F    : Enforce symmetry in interaction datasets [True]
    dbsource=X      : Source database (biogrid/dip/intact/mint/reactome) [biogrid]
    mitab=T/F       : Whether source file is in MITAB flat file format [True]
    species=X       : Name of species to use data for (will be read from file if BioGRID) [human]
    taxid=LIST      : List of NCBI Taxa IDs to use (for DIP and Domino) [9606]
    unipath=PATH    : Path to UniProt files [UniProt/]

    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ppifas=T/F      : Whether to output PPI datasets as fasta files into Species/BIOGRID_Datasets/ [True]
    minseq=X        : Minimum number of PPI sequences in order to output fasta file [3]
    ppitab=T/F      : Whether to output PPI table with aliases etc. [True]
    alltypes=T/F    : Output a full list of PPITypes. (Will populate the PPIType list) [False]

    ### ~ Special Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hostvirus=T/F   : Whether to pull out host-virus interactions only (MINT/IntAct only) [False]
    vcodes=LIST     : List/File of viral species codes for IntAct hostvirus=T []
    
Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_seq, rje_dismatrix
Other modules needed: None
"""
# TaxaID: C elegans = 6239
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_uniprot, rje_zen
#import rje_ppi, rje_genecards
#import rje_dismatrix_V2 as rje_dismatrix        #!# Check V2 is OK! (Might need V1)
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation with BioGRID Flat-File parsing information and sequence extraction for yeast.
    # 1.0 - Added cross-referencing via GeneCards output to generate Human Datasests.
    # 1.1 - Added IntAct and MINT parsing.
    # 1.2 - Add option to pull out host-virus interactions.
    # 1.3 - Added Reactome & DIP parsing.
    # 1.4 - Added rje_genemap object functionality.
    # 1.5 - Added Domino parsing and tracking of evidence codes.
    # 1.6 - Updated BioGRID parsing to use mitab format.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Get basic parser working.
    # [ ] : Ensure all organisms can be parsed and mapped to sequences. Recognised different taxa names.
    # [ ] : Handle interactions where proteins are from different organisms.
    # [ ] : Improve parser to spot and handle errors in file.
    # [X] : Incorporate HPRD for human interaction data? - See RJE_HPRD
    # [ ] : Cytoscape output.
    # [Y] : Incorporate MINT data.
    # [Y] : Change opt['IntAct'] and opt['MINT'] to dbsource=X and try to auto-detect.
    # [ ] : Update to allow parsing of the full data table.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_BioGRID', '1.6', 'May 2010', '2007')
    description = 'BioGRID Database processing module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
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
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: BioGRID Class                                                                                           #
#########################################################################################################################
class BioGRID(rje.RJE_Object):     
    '''
    BioGRID parsing and data management class. Author: Rich Edwards (2007).

    Info:str
    - PPIFile = PPI database flat file [None]
    - DBSource = Source database (biogrid/intact/mint) [biogrid]    
    - GeneCards = File of links between IDs. For human, should have HGNC and EnsLoci columns. [None]    
    - Species = BioGRID organism, read from file name [Unknown]
    
    Opt:boolean
    - IntAct = Whether the file is actually an IntAct download [False]
    - MITAB = Whether source file is in MITAB flat file format [True]
    - PPIFas = Whether to output PPI datasets as fasta files into Species/BIOGRID_Datasets/ [True]
    - PPITab = Whether to output PPI table with aliases etc. [True]
    - Symmetry = Enforce symmetry in interaction datasets [True]
    - AllTypes = Output a full list of PPITypes. (Will populate the PPIType list) [False]
    - HostVirus = Whether to pull out host-virus interactions only (MINT only) [False]
    
    Stat:numeric
    - MinSeq = Minimum number of PPI sequences in order to output fasta file [3]

    List:list
    - PPIType = List of acceptable interaction types to parse out []
    - BadType = List of bad interaction types, to exclude [indirect_complex,neighbouring_reaction]
    - TaxID = List of NCBI Taxa IDs to use (for DIP) [9606]
    - VCodes = List/File of viral species codes for IntAct hostvirus=T []

    Dict:dictionary
    - Evidence = Evidences for interactions {ID:{ID:[Evidence types]}}
    - Protein = Data for each protein {ID:{Seq object,Gene,Aliases}}
    - PPI = Interaction data {ID:[IDs]}
    - Mapping = Dictionary of backwards mapping {Link:[IDs]}  
    
    Obj:RJE_Object
    - GeneMap = GeneMap object for mapping onto HGNC etc. genes where possible
    - SeqList = SeqList object containing actual sequences, loaded from seqin=FILE.
    - UniProt = rje_uniprot.UniProt object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object'''
        ### Basics ###
        self.infolist = ['PPFile','GeneCards','Species','DBSource']
        self.optlist = ['PPIFas','PPITab','Symmetry','AllTypes','HostVirus','MITAB']
        self.statlist = ['MinSeq']
        self.listlist = ['PPIType','BadType','VCodes','TaxID']
        self.dictlist = ['Protein','PPI','Mapping','Evidence']
        self.objlist = ['GeneMap','SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Species':'human','DBSource':'biogrid'})
        self.setStat({'MinSeq':3})
        self.setOpt({'HostVirus':False,'AllTypes':False})
        self.list['PPIType'] = []   #['Two-hybrid','Protein-peptide','Co-crystal Structure']
        self.list['BadType'] = ['indirect_complex','neighbouring_reaction']
        self.list['TaxID'] = ['9606']
        self.obj['SeqList'] = rje_seq.SeqList(self.log,['seqnr=F','accnr=F']+self.cmd_list+['autoload=T','type=protein','align=False'])
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
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['PPIFile','GeneCards'])
                self._cmdReadList(cmd,'info',['Species','DBSource'])
                self._cmdReadList(cmd,'opt',['PPIFas','PPITab','Symmetry','AllTypes','HostVirus','MITAB'])
                self._cmdReadList(cmd,'int',['MinSeq'])
                self._cmdReadList(cmd,'list',['PPIType','BadType','VCodes','TaxID'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        self.info['DBSource'] = self.info['DBSource'].lower()
        if self.opt['AllTypes']: self.list['PPIType'] = []
        if self.list['PPIType']: self.list['PPIType'] = rje.split(rje.join(self.list['PPIType'],'|').lower(),'|')
        if self.list['BadType']: self.list['BadType'] = rje.split(rje.join(self.list['BadType'],'|').lower(),'|')
        if self.opt['HostVirus']: self.info['Species'] = '%s-virus' % self.info['Species']
#########################################################################################################################
    ### <2> ### Main Class Run Methods                                                                                  #
#########################################################################################################################
    def run(self):  ### Main Run Method
        '''Main Run Method.'''
        try:
            ### ~ Parse BioGRID Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.parse()    # Generates a dictionary of PPI using Proteins as keys/lists with some mapping in self.dict
            self.mapSeq()   # Maps sequences from SeqList onto proteins in self.dict  
            #Q# Do we want to add aliases to dict['Mapping'] if not already used.
            
            ### ~ Standard Run Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['PPITab']: self.ppiTabOut()
            if self.opt['PPIFas'] and self.obj['SeqList'].seqNum(): self.saveFasta()
            if self.opt['AllTypes']: self.allTypesOut()
        except: self.log.errorLog('Major problem with rje_biogrid')
#########################################################################################################################
    def parse(self):    ### Main Parsing method. Generates Mappings and PPI data dictionary
        '''BioGRID Parsing method. Generates Mappings between IDs and PPI dictionary.'''
        try:
            ### ~ Open file and read one line at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(self.info['PPIFile']):
                self.log.errorLog('PPI database file "%s" missing!' % self.info['PPIFile'],printerror=False)
                return False
            if self.info['DBSource'] == 'biogrid':
                try: self.info['Species'] = rje.split(self.info['PPIFile'],'-')[2]
                except: pass
            if self.info['DBSource'] == 'reactome':
                try: self.info['Species'] = rje.split(rje.baseFile(self.info['PPIFile'],strip_path=True),'.')[0]
                except: pass
            ## Open file and parse through headers ##
            GRID = open(self.info['PPIFile'],'r')
            line = GRID.readline()
            if self.info['DBSource'] == 'biogrid' and not self.opt['MITAB']:
                while line and line.find('INTERACTOR_A') != 0: line = GRID.readline()
            ## Read in actual data ##
            px = 0
            self.log.printLog('\r#PPI','Parsing %s interactions: %s prot; %s ppi.' % (self.info['DBSource'],rje.integerString(len(self.dict['PPI'])),rje.integerString(px)),newline=False,log=False)
            while line:
                ## Parse line and check ##
                data = rje.readDelimit(line,delimit='\t')
                line = GRID.readline()
                if self.info['DBSource'] == 'reactome':
                    if data[0][:1] == '#': continue
                    if len(data) < 8: continue
                elif len(data) < 11: continue
                if self.info['DBSource'] == 'biogrid' and data[0] == 'INTERACTOR_A': continue
                if self.info['DBSource'] in ['dip','intact','domino','biogrid'] and data[0] == 'ID interactor A': continue
                if self.info['DBSource'] in ['dip','intact','domino','biogrid'] and data[0] == '#ID Interactor A': continue
                if self.info['DBSource'] == 'mint' and data[0].find('ID interactor') == 0: continue
                if self.info['DBSource'] == 'reactome': [p1,g1,a1,p2,g2,a2,type] = data[:7]
                else: [p1,p2,g1,g2,a1,a2,type] = data[:7]
                type = type.lower()
                if self.info['DBSource'] == 'intact':
                    if self.opt['HostVirus']:   # Must have one virus and one human
                        try:
                            s1 = rje.matchExp('TAXID:\d+\((\S+)\)',data[9].upper())[0]
                            s2 = rje.matchExp('TAXID:\d+\((\S+)\)',data[10].upper())[0]
                            if s1 not in self.list['VCodes'] and s2 not in self.list['VCodes']: continue
                            if s1 != self.getStr('Species').upper() and s2 != self.getStr('Species').upper(): continue    
                        except: continue
                    elif data[9].lower().find(self.info['Species'].lower()) < 0 or data[10].lower().find(self.info['Species'].lower()) < 0: continue
                    (p1,p2,g1,g2,a1,a2,type) = self.convertMITAB(p1,p2,g1,g2,a1,a2,type)
                elif self.info['DBSource'] == 'biogrid' and self.opt['MITAB']:
                    (p1,p2,g1,g2,a1,a2,type) = self.convertMITAB(p1,p2,g1,g2,a1,a2,type)
                elif self.info['DBSource'] == 'mint' and self.opt['MITAB']:     # For data parsed from the TXT file
                    if p1 == '-' or p2 == '-' or not p1 or not p2: continue
                    p1 = rje.replace(p1,';','|')
                    if p1[:1] == '|': p1 = p1[1:]
                    p2 = rje.replace(p2,';','|')
                    if p2[:1] == '|': p2 = p2[1:]
                    g1 = rje.replace(g1,';','|')
                    if g1[:1] == '|': g1 = g1[1:]
                    g2 = rje.replace(g2,';','|')
                    if g2[:1] == '|': g2 = g2[1:]
                    a1 = rje.replace(a1,';','|')
                    if a1[:1] == '|': a1 = a1[1:]
                    a2 = rje.replace(a2,';','|')
                    if a2[:1] == '|': a2 = a2[1:]
                    (p1,p2,g1,g2,a1,a2,type) = self.convertMITAB(p1,p2,g1,g2,a1,a2,type)
                elif self.info['DBSource'] == 'mint':      # For data parsed from the TXT file
                    ## Check type ##
                    try: type = rje.matchExp('MI:\d+\((.+)\)',data[8])[0].lower()
                    except: type = 'unknown'
                    #if not self.opt['AllTypes'] and self.list['PPIType'] and type not in self.list['PPIType']: continue
                    #elif type in self.list['BadType']: continue
                    ## Check Taxa ##
                    # a1 and a2 are taxa IDs   taxid:9606(Homo sapiens)        taxid:9606(Homo sapiens)
                    try:
                        if self.opt['HostVirus']:   # Must have one virus and one human
                            if a1.lower().find('virus') < 0 and a2.lower().find('virus') < 0: continue
                            elif a1 != 'taxid:9606(Homo sapiens)' and a2 != 'taxid:9606(Homo sapiens)': continue
                        elif rje.matchExp('taxid:(\d+)',a1)[0] not in self.list['TaxaID'] or rje.matchExp('taxid:(\d+)',a2)[0] not in self.list['TaxaID']: continue       #!# Update to treat species properly with species dictionary #!#
                    except: self.deBug(a1); continue
                    a1 = a2 = ''
                    ## Convert Protein AccNum ##
                    p1 = rje.split(rje.split(p1,':')[1],'-')[0]   # Splice isoforms are removed
                    p1 = rje.split(p1,';')[0]
                    if rje.matchExp('^[pqo](\d\S+)',p1): p1 = p1.upper()
                    p2 = rje.split(rje.split(p2,':')[1],'-')[0]   # Splice isoforms are removed
                    p2 = rje.split(p2,';')[0]
                    if rje.matchExp('^[pqo](\d\S+)',p2): p2 = p2.upper()
                    if g1.find('/') >= 0: g1 = ''
                    if g2.find('/') >= 0: g2 = ''
                    if g1.find('_human') > 0: g1 = rje.split(rje.replace(g1,'_human',''),'-')[0].upper()  # Convert to UniProt, no splice isoforms
                    if g2.find('_human') > 0: g2 = rje.split(rje.replace(g2,'_human',''),'-')[0].upper()  # Convert to UniProt, no splice isoforms
                    if g1 == '-' or g2 == '-': self.deBug([p1,p2,g1,g2,a1,a2])
                elif self.info['DBSource'] in ['dip','domino']:
                    self.deBug(data)
                    ## Check type ##
                    try:
                        types = []
                        for mi in rje.split(data[6],'|'):
                            try: types.append(rje.matchExp('MI:\d+\((.+)\)',mi)[0].lower())
                            except: pass
                    except: types = ['unknown']
                    #for type in types[0:]:
                    #    if not self.opt['AllTypes'] and self.list['PPIType'] and type not in self.list['PPIType']: types.remove(type)
                    #    elif type in self.list['BadType']: types.remove(type)
                    #if not types: continue
                    type = rje.join(types,'|')
                    ## Check Taxa ##
                    badtaxa = False
                    for taxa in data[9:11]:
                        self.deBug('%s = %s' % (taxa,rje.matchExp('taxid:(\S+)\(',taxa)))
                        try:
                            if self.info['DBSource'] in ['dip'] and rje.matchExp('taxid:(\S+)\(',taxa)[0] not in self.list['TaxID']: badtaxa = True
                            if self.info['DBSource'] in ['domino'] and rje.matchExp('taxid:(\S+)\(',taxa)[0] not in self.list['TaxID']: badtaxa = True
                        except: badtaxa = True
                    if badtaxa: continue
                    if self.info['DBSource'] == 'domino' and self.opt['MITAB']: (p1,p2,g1,g2,a1,a2,type) = self.convertMITAB(p1,p2,g1,g2,a1,a2,type)
                    else: #!# UniProt only #!#
                        self.deBug('%s v %s' % (p1,p2))
                        g1 = g2 = a1 = a2 = ''
                        for p in rje.split(p1,'|'):
                            if rje.matchExp('uniprotkb:(\S+)',p): g1 = rje.matchExp('uniprotkb:(\S+)',p)[0]
                        p1 = g1
                        for p in rje.split(p2,'|'):
                            if rje.matchExp('uniprotkb:(\S+)',p): g2 = rje.matchExp('uniprotkb:(\S+)',p)[0]
                        p2 = g2
                        self.deBug('-> %s v %s' % (p1,p2))
                        if not p1 or not p2: continue
                elif self.info['DBSource'] == 'reactome':
                    ## Check Type = "reaction" ##
                    #x#if type != 'reaction': continue
                    ## Split and process each gene ##
                    if not g1: g1 = a1
                    if not g1: g1 = p1
                    if not g2: g2 = a2
                    if not g2: g2 = p2
                    p1 = rje.split(p1,':')[1]
                    p2 = rje.split(p2,':')[1]
                    a1 = a2 = ''
                    glist1 = rje.split(g1,'|')
                    glist2 = rje.split(g2,'|')
                    for gs1 in glist1:
                        g1 = rje.split(gs1,':')[1]
                        for gs2 in glist2: g2 = rje.split(gs2,':')[1]
                        #x#print p1,p2,g1,g2,a1,a2
                        self.updatePPI(p1,p2,g1,g2,a1,a2,type)
                    px += 1
                    self.log.printLog('\r#PPI','Parsing %s interactions: %s prot; %s ppi.' % (self.info['DBSource'],rje.integerString(len(self.dict['PPI'])),rje.integerString(px)),newline=False,log=False)
                    continue
                #else:
                #    if not self.opt['AllTypes'] and self.list['PPIType'] and type not in self.list['PPIType']: continue
                #    elif type in self.list['BadType']: continue
                #x#self.deBug([p1,p2,g1,g2,a1,a2])
                self.updatePPI(p1,p2,g1,g2,a1,a2,type)
                px += 1
                self.log.printLog('\r#PPI','Parsing %s interactions: %s prot; %s ppi.' % (self.info['DBSource'],rje.integerString(len(self.dict['PPI'])),rje.integerString(px)),newline=False,log=False)
            self.log.printLog('\r#PPI','Parsing %s interactions: %s prot; %s ppi.' % (self.info['DBSource'],rje.integerString(len(self.dict['PPI'])),rje.integerString(px)))
            GRID.close()
            if self.opt['AllTypes']: self.list['PPIType'].sort(); self.log.printLog('#TYPE',rje.join(self.list['PPIType'],' | '))
            return True
        except:
            self.log.errorLog('Error in BioGRID.parse()',printerror=True,quitchoice=False)
            return False
#########################################################################################################################
    def convertIntAct(self,p1,p2,g1,g2,a1,a2):  ### Converts Intact formatting into standard BioGrid
        '''Converts Intact formatting into standard BioGrid.'''
        try:
            p1 = rje.split(rje.split(p1,'|')[0],':')[1]
            p2 = rje.split(rje.split(p2,'|')[0],':')[1]
            if g1.find(':') > 0: g1 = rje.split(rje.split(g1,':')[1],'(')[0]
            if g2.find(':') > 0: g2 = rje.split(rje.split(g2,':')[1],'(')[0]
            if a1.find(':') > 0:
                alist = []
                for a in rje.split(a1,'|'): alist.append(rje.split(a,':')[1])
                a1 = rje.join(alist,'|')
            if a2.find(':') > 0:
                alist = []
                for a in rje.split(a2,'|'): alist.append(rje.split(a,':')[1])
                a2 = rje.join(alist,'|')
        except: self.log.errorLog('Problem during convertIntAct()')
        return (p1,p2,g1,g2,a1,a2)
#########################################################################################################################
    def convertMITAB(self,p1,p2,g1,g2,a1,a2,type):  ### Converts Intact formatting into standard BioGrid
        '''Converts Intact formatting into standard BioGrid.'''
        try:
            p1 = rje.split(rje.split(p1,'|')[0],':')[1]
            p2 = rje.split(rje.split(p2,'|')[0],':')[1]
            if self.info['Species'] == 'Caenorhabditis_elegans':
                if g1.find(':') > 0: g1 = rje.split(rje.split(g1,'|')[0],':')[1]
                if g2.find(':') > 0: g2 = rje.split(rje.split(g2,'|')[0],':')[1]
            else:
                if g1.find(':') > 0: g1 = rje.split(rje.split(g1,':')[1],'(')[0]
                if g2.find(':') > 0: g2 = rje.split(rje.split(g2,':')[1],'(')[0]
            if a1.find(':') > 0:
                alist = []
                for a in rje.split(a1,'|'): alist.append(rje.split(rje.split(a,':')[1],'(')[0])
                a1 = rje.join(alist,'|')
            if a2.find(':') > 0:
                alist = []
                for a in rje.split(a2,'|'): alist.append(rje.split(rje.split(a,':')[1],'(')[0])
                a2 = rje.join(alist,'|')
            types = []
            for mitype in rje.split(type,'|'):
                try: types.append(rje.matchExp('\((.+)\)',mitype)[0])
                except: continue
            type = rje.join(types,'|')
        except:
            self.errorLog('Problem during convertMITAB(%s)' % rje.join([p1,p2,g1,g2,a1,a2,type],' || '))
        return (p1,p2,g1,g2,a1,a2,type)
#########################################################################################################################
    def updatePPI(self,p1,p2,g1,g2,a1,a2,types=None):   ### Adds an interaction to the dictionaries
        '''Adds an interaction to the dictionaries.'''
        ### ~ [1] ~ Check type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not types: types = 'unknown'
        types = rje.split(types,'|')
        for type in types[0:]:
            if self.opt['AllTypes']:
                if type not in self.list['PPIType']: self.list['PPIType'].append(type)
            else:
                if type in self.list['BadType'] or (self.list['PPIType'] and type not in self.list['PPIType']):  types.remove(type)
        if not types: return    # No acceptable experiment types remain
        ### ~ [2] ~ Make upper case and strip splice variants of proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try: p1 = rje.matchExp('^([OPQ]\S+)\-\d+',p1.upper())[0]
        except: p1 = p1.upper()
        try: p2 = rje.matchExp('^([OPQ]\S+)\-\d+',p2.upper())[0]
        except: p2 = p2.upper()
        g1 = g1.upper()
        if g1 in ['','-','N/A']: g1 = p1
        g2 = g2.upper()
        if g2 in ['','-','N/A']: g2 = p2
        ### ~ [3] ~ Update dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if p1 not in self.dict['PPI']: self.dict['PPI'][p1] = []
        if p2 not in self.dict['PPI']: self.dict['PPI'][p2] = []
        if p2 not in self.dict['PPI'][p1]: self.dict['PPI'][p1].append(p2)
        if self.opt['Symmetry'] and p1 not in self.dict['PPI'][p2]: self.dict['PPI'][p2].append(p1)
        ## ~ [3a] ~ Evidence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if p1 not in self.dict['Evidence']: self.dict['Evidence'][p1] = {}
        if p2 not in self.dict['Evidence'][p1]: self.dict['Evidence'][p1][p2] = []
        for type in types:
            if type not in self.dict['Evidence'][p1][p2]: self.dict['Evidence'][p1][p2].append(type)
        ## ~ [3b] ~ Links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if p1 not in self.dict['Protein']: self.dict['Protein'][p1] = {'Gene':g1,'Alias':rje.split(a1,'|')}
        if g1 not in self.dict['Mapping']: self.dict['Mapping'][g1] = []
        if p1 not in self.dict['Mapping'][g1]: self.dict['Mapping'][g1].append(p1)
        if p2 not in self.dict['Protein']: self.dict['Protein'][p2] = {'Gene':g2,'Alias':rje.split(a2,'|')}
        if g2 not in self.dict['Mapping']: self.dict['Mapping'][g2] = []
        if p2 not in self.dict['Mapping'][g2]: self.dict['Mapping'][g2].append(p2)
#########################################################################################################################
    def updateGenes(self,genecards=None,species='HUMAN'):  ### Uses Genecards and UniProt to extract gene and alias information
        '''
        Uses Genecards and UniProt to extract gene and alias information.
        >> genecards:Genecards object [None] = Should have dict['CardMap'] with UniProt etc. mappings
        >> species:str [HUMAN] = species code to check (UniProt)
        '''
        try:### ~ [1] Look in Genecards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []    # List of missing protein IDs
            if genecards:
                try: cards = genecards.dict['GeneCard']     # Old GeneCards object
                except:
                    cards = genecards.dict['Data']      # New GeneMap object
                    genecards.makeGeneMap()             # For compatibility with old GeneCards
                cardmap = genecards.dict['CardMap']
                (px,pnum) = (0.0,len(self.dict['Protein']))
                for p1 in rje.sortKeys(self.dict['Protein']):       # Contains all proteins
                    symbol = None
                    for g1 in [self.dict['Protein'][p1]['Gene']] + self.dict['Protein'][p1]['Alias']:
                        try: symbol = cards[g1]['Symbol']
                        except:
                            try: symbol = cardmap[g1]
                            except: symbol = None
                        #if g1 in cards: symbol = cards[g1]['Symbol']          #x#continue    # No need to try and map - already in GeneCards
                        #elif g1 in cardmap: symbol = cardmap[g1]
                        #if symbol == '!FAILED!': symbol = g1
                        if symbol: break                # Found a GeneCard entry, so stop
                    if symbol:
                        if symbol == self.dict['Protein'][p1]['Gene']: continue
                        else:
                            if self.dict['Protein'][p1]['Gene'] not in self.dict['Protein'][p1]['Alias']: self.dict['Protein'][p1]['Alias'].append(self.dict['Protein'][p1]['Gene'])
                            self.dict['Protein'][p1]['Gene'] = symbol
                    else: missing.append(p1)                    # Will look for UniProt entry
                    px += 100.0
                    self.log.printLog('\r#CARDS','Mapping %s proteins to GeneCards: %.1f%% - %s missing' % (rje.integerString(pnum),px/pnum,rje.integerString(len(missing))),newline=False,log=False)
            else: missing = rje.sortKeys(self.dict['Protein'])  # Look for all proteins!
            self.log.printLog('\r#CARDS','Mapping %s proteins to GeneCards complete: %s missing' % (rje.integerString(pnum),rje.integerString(len(missing))))
            ### ~ [2] Look in UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mx, mnum) = (0.0,len(missing))
            rejected = []
            if missing:     #!# Add mapping onto Gene Symbols through UniProt??? Return Entry and use DBLinks! #!#
                ## ~ [2a] Extract missing sequences from UniProt, where possible ~~~~~~~~~~~~~~~~~~ ##
                self.obj['UniProt'] = rje_uniprot.UniProt(self.log,self.cmd_list)
                self.log.printLog('#ACC','Looking for UniProt AccNum...',newline=False,log=False)
                accdict = self.obj['UniProt'].accDict(acc_list=missing)
                for acc in missing[0:]:
                    mx += 100.0
                    if acc in accdict:
                        entry = accdict[acc]
                        if species and not entry.isSpecies(species):
                            self.log.printLog('#SPEC','Rejected %s - not a %s sequence!' % (acc,species),screen=False)
                            missing.remove(acc)
                            self.dict['PPI'].pop(acc)
                            rejected.append(acc)
                            for p1 in self.dict['PPI']:
                                if acc in self.dict['PPI'][p1]: self.dict['PPI'][p1].remove(acc)
                            continue
                        symbol = None
                        for db in ['Symbol','EnsG','HGNC','Entrez']:
                            if db not in entry.dict['DB']: continue
                            for link in entry.dict['DB'][db]:
                                if db == 'HGNC': g1 = 'HGNC%s' % link   #!# Colon removed for GeneMap
                                else: g1 = link
                                try: symbol = cards[g1]['Symbol']
                                except:
                                    try: symbol = cardmap[g1]
                                    except: symbol = None
                                #if genecards: and g1 in cards: symbol = cards[g1]['Symbol']
                                #elif genecards and g1 in cardmap: symbol = cardmap[g1]
                                #if symbol == '!FAILED!': symbol = g1
                                if symbol: break                # Found a GeneCard entry, so stop
                                elif g1 not in self.dict['Protein'][acc]['Alias']: self.dict['Protein'][acc]['Alias'].append(g1)
                            if symbol: break                # Found a GeneCard entry, so stop
                        if symbol:
                            if symbol == self.dict['Protein'][acc]['Gene']: continue
                            else:
                                if self.dict['Protein'][acc]['Gene'] not in self.dict['Protein'][acc]['Alias']: self.dict['Protein'][acc]['Alias'].append(self.dict['Protein'][acc]['Gene'])
                                self.dict['Protein'][acc]['Gene'] = symbol
                                missing.remove(acc)
                    self.log.printLog('\r#SEQ','Mapping missing proteins to UniProt: %.1f%%; %s missing' % ((mx/mnum),rje.integerString(len(missing))),newline=False,log=False)
                self.log.printLog('\r#SEQ','Mapping missing proteins to UniProt complete: %s missing' % (rje.integerString(len(missing))))
            ### ~ Output Missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if missing: open('%s.%s.missing.txt' % (self.info['DBSource'],self.info['Species']),'w').write(rje.join(missing,'\n'))
            if rejected: self.log.printLog('#SPEC','%s sequences rejected - not %s.' % (rje.integerString(len(rejected)),species))
        except: self.log.errorLog('Major problem with BioGRID.updateGenes()')
#########################################################################################################################
    def mapSeq(self):   ### Maps sequences from SeqList onto proteins in self.dict
        '''Maps sequences from SeqList onto proteins in self.dict.'''
        try:
            ### ~ Setup: check sequences loaded and make name dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#if self.obj['SeqList'].seqNum() < 1: return self.log.printLog('#SEQ','No sequences loaded')  # UniProt
            namedict = self.obj['SeqList'].seqNameDic('Max')
            ### ~ Optional use of GeneCards to flesh out links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Species'].lower() in ['homo_sapiens','human'] and os.path.exists(self.info['GeneCards']):
                update = rje.dataDict(self,self.info['GeneCards'],getheaders=True)
                if not 'EnsLoci' in update.pop('Headers'):
                    self.log.errorLog('EnsLoci column missing from GeneCards data',printerror=False)
                    update = {}
                gx = 0
                for gene in update:
                    ensloci = update[gene]['EnsLoci']
                    if not ensloci in namedict: continue    # Cannot map
                    if update[gene]['Symbol']: namedict[update[gene]['Symbol']] = namedict[ensloci]
                    if update[gene]['HGNC']: namedict['HGNC%s' % update[gene]['HGNC']] = namedict[ensloci]
                    if update[gene]['UniProt']: namedict[update[gene]['UniProt']] = namedict[ensloci] 
                    gx += 1
                    self.log.printLog('\r#LINK','Linking protein names using GeneCard file: %s links' % rje.integerString(gx),newline=False,log=False)
                self.log.printLog('\r#LINK','Linked protein names using GeneCard file: %s links made.' % rje.integerString(gx))                    
            elif self.info['Species'].lower() in ['homo_sapiens','human'] and not os.path.exists(self.info['GeneCards']):
                self.log.printLog('#LINK','No linking with GeneCard file: "%s" not found' % self.info['GeneCards'])
            ### ~ Make sequence links, where possible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []    # List of missing protein IDs
            (px,pnum) = (0.0,len(self.dict['Protein']))
            for p1 in self.dict['Protein']:
                if p1 in namedict: self.dict['Protein'][p1]['Seq'] = namedict[p1]
                else:
                    try: self.dict['Protein'][p1]['Seq'] = namedict[self.dict['Protein'][p1]['Gene']]
                    except: missing.append(p1)
                px += 100.0
                self.log.printLog('\r#SEQ','Mapping %s proteins to sequences: %.1f%% - %s missing' % (rje.integerString(pnum),px/pnum,rje.integerString(len(missing))),newline=False,log=False)
            self.log.printLog('\r#SEQ','Mapping %s proteins to sequences complete: %s missing' % (rje.integerString(pnum),rje.integerString(len(missing))))
            ### ~ Extract missing sequences from UniProt, where possible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if missing:     #!# Add mapping onto Gene Symbols through UniProt??? Return Entry and use DBLinks! #!#
                self.obj['UniProt'] = rje_uniprot.UniProt(self.log,self.cmd_list)
                self.log.printLog('#ACC','Looking for UniProt AccNum...',newline=False,log=False)
                (accnamedict,accseqdict) = self.obj['UniProt'].accNameSeq(missing,justsequence=False)
                for acc in missing[0:]:
                    if acc in accseqdict:
                        seq = accseqdict[acc]
                        if self.dict['Protein'][acc]['Gene'] == acc: self.dict['Protein'][acc]['Gene'] = seq.info['Gene'].upper()
                        else: self.dict['Protein'][acc]['Alias'].append(seq.info['Gene'].upper())   # Add to Alias for Pingu mapping
                        self.dict['Protein'][acc]['Alias'].append(seq.info['AccNum'].upper())       # Add to Alias for Pingu mapping
                        if seq.shortName() in namedict: self.dict['Protein'][acc]['Seq'] = namedict[seq.shortName()]
                        else:
                            for i in ['Gene','AccNum']:
                                if seq.info[i].upper() in namedict: self.dict['Protein'][acc]['Seq'] = namedict[seq.info[i].upper()]
                    if 'Seq' in self.dict['Protein'][acc]:
                        missing.remove(acc)
                        continue
                    if acc in accnamedict:
                        self.dict['Protein'][acc]['Seq'] = self.obj['SeqList']._addSeq(accnamedict[acc],accseqdict[acc].info['Sequence'])
                        missing.remove(acc)
                    self.log.printLog('\r#SEQ','Mapping missing proteins to UniProt: %s missing' % (rje.integerString(len(missing))),newline=False,log=False)
                self.log.printLog('\r#SEQ','Mapping missing proteins to UniProt complete: %s missing' % (rje.integerString(len(missing))))
            ### ~ Output Missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if missing: open('%s.%s.missing.txt' % (self.info['DBSource'],self.info['Species']),'w').write(rje.join(missing,'\n'))
        except: self.log.errorLog('Error during BioGRID.mapSeq()')
#########################################################################################################################
    def ppiTabOut(self):    ### Outputs table of protein links etc.
        '''Outputs table of protein links etc.'''
        try:### ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tfile = '%s.%s.tdt' % (self.info['DBSource'],self.info['Species'])
            headers = ['Protein','Symbol','Aliases','PPI','Prot','Desc']
            rje.delimitedFileOutput(self,tfile,headers,'\t',rje_backup=True)
            ### ~ Output Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (px,pnum) = (0.0,len(self.dict['Protein']))
            for p in rje.sortKeys(self.dict['Protein']):
                datadict = {'Protein':p,'Symbol':self.dict['Protein'][p]['Gene'],'PPI':len(self.dict['PPI'][p]),'Prot':'!FAILED!'}
                if self.dict['Protein'][p]['Alias'] not in [[''],['N/A'],[],['-']]: datadict['Aliases'] = rje.join(self.dict['Protein'][p]['Alias'],'; ')+';'
                if 'Seq' in self.dict['Protein'][p]:
                    datadict['Prot'] = self.dict['Protein'][p]['Seq'].shortName()
                    datadict['Desc'] = self.dict['Protein'][p]['Seq'].info['Description']
                rje.delimitedFileOutput(self,tfile,headers,'\t',datadict)
                px += 100.0
                self.log.printLog('\r#TAB','Tabulating data for %s proteins to %s: %.1f%%' % (rje.integerString(pnum),tfile,px/pnum),newline=False,log=False)
            self.log.printLog('\r#TAB','Tabulating data for %s proteins to %s complete.' % (rje.integerString(pnum),tfile))
        except: self.log.errorLog('Problem during BioGRID.ppiTabOut()')            
#########################################################################################################################
    def saveFasta(self):    ### Outputs parsed PPI datasets in Fasta format
        '''Outputs parsed PPI datasets in Fasta format.'''
        try:
            ### ~ Output PPI Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datpath = rje.makePath('%s/%s_Datasets/' % (self.info['Species'],{'biogrid':'BioGRID','intact':'IntAct','mint':'MINT','dip':'DIP'}[self.info['DBSource']]))
            rje.mkDir(self,datpath)
            fx = 0
            for p1 in rje.sortKeys(self.dict['PPI']):
                mylist = []
                for p2 in self.dict['PPI'][p1]:
                    try: mylist.append(self.dict['Protein'][p2]['Seq'])
                    except: self.log.printLog('#MISS','Protein sequence for %s missing (%s PPI)' % (p2,p1))
                if self.dict['Protein'][p1]['Gene'] == 'N/A': sfile = '%s%s_biogrid.fas' % (datpath,p1)
                else: sfile = '%s%s_biogrid.fas' % (datpath,self.dict['Protein'][p1]['Gene'])
                sfile = rje.replace(sfile,'biogrid',self.info['DBSource'])
                if len(mylist) >= self.stat['MinSeq']:
                    self.obj['SeqList'].saveFasta(seqs=mylist,seqfile=sfile)
                    fx += 1
                else: self.log.printLog('#MIN','Not enough PPIs (%d) for %s output.' % (len(mylist),p1))
            self.log.printLog('#FAS','%s PPI Fasta files output to %s' % (rje.integerString(fx),datpath))
        except:
            self.log.errorLog('Error in BioGRID.saveFasta()',printerror=True,quitchoice=False)
            raise
#########################################################################################################################
    def allTypesOut(self):  ### Outputs a list of possible PPI Types
        '''Outputs a list of possible PPI Types.'''
        tfile = '%s.ppi_types.txt' % self.info['DBSource']
        open(tfile,'w').write(rje.join(rje.sortUnique(self.list['PPIType']),'\n'))
        self.log.printLog('#TYPE','%s PPI Types output to %s' % (len(self.list['PPIType']),tfile))
#########################################################################################################################
### End of SECTION III: BioGRID Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MAIN PROGRAM                                                                                             #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: BioGRID(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################
