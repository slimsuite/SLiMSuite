#!/usr/bin/python

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
Module:       rje_ensembl
Description:  EnsEMBL Processing/Manipulation Module 
Version:      2.15.2
Last Edit:    20/04/15
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for processing EnsEMBL data for the rje_dbase module. The main class is an EnsEMBL class, which stores
    information on EnsEMBL proteins in terms of their gene IDs, loci and descriptions. This generates the "EnsLoci"
    dataset for each genome, consisting of the "best" peptide for a given locus. For known genes, UniProt accession
    numbers will be used in place of the EnsEMBL accession number. If the EnsEMBL sequence maps to a SwissProt sequence
    but is of really low quality (20+ consecutive Xs with less non-X residues than the SwissProt sequence) then the
    SwissProt sequence itself will replace the EnsEMBL sequence. This is the only time that the relationship between
    EnsEMBL peptide ID and sequence will break down.

    Version 1.7 introduced a new "EnsGO" function for making GO datasets for the species codes listed. This mode will
    need, for each SPECIES, the EnsLoci file enspath/ens_SPECIES.loci.fas, the GO mapping enspath/ens_SPECIES.GO.tdt
    and the GO ID file [GO.terms_ids_obs]. GO mapping files can be created for the relevant species using EnsEMBL's
    BioMart tool (http://www.ensembl.org/biomart/martview/), while the ID file can be downloaded from GO
    (http://www.geneontology.org/doc/GO.terms_ids_obs). From BioMart, the following columns should be downloaded:
    "Ensembl Gene ID","Ensembl Transcript ID","Ensembl Peptide ID","GO ID","GO description","GO evidence code",
    "EntrezGene ID","HGNC Symbol". Other fields can also be downloaded if desired. This function has been further
    updated in version 1.8 & 1.9. From Version 2.8, the columns should be: "Ensembl Gene ID", "Ensembl Transcript ID",
    "Ensembl Protein ID", "GO Term Accession", "GO Term Evidence Code", "EntrezGene ID", "HGNC symbol"

    Version 2.0 introduced a new "EnsDat" function for generating fake UniProt format entries for EnsLoci data using
    PFam HMM domain prediction, TMHMM transmembrane topology prediction, SIGNALP signal peptide prediction and IUPRED
    disorder prediction. Assumes that the EnsLoci files have been created. (Use download=T ensloci=T if not!) Sequences
    should be extracted from the file created by this method using Accession Numbers only.

    Version 2.11 is the start of a major reworking in preparation for V3.0. Species codes are now read in automatically
    and Ensembl species alone downloaded from Uniprot for EnsLoci processing. (This can be quite slow depending on
    connection etc.) This avoids the need for pre-processing Uniprot in order to make EnsLoci sequences. Modified Uniprot
    downloads and data extraction is used for db xref mapping in place of manual biomart tables. Species data is now
    split into subsets according to Ensembl sets (main, metazoa, protists etc.) and EnsLoci files are similarly split
    within an `ensloci/` subdirectory of `enspath/`.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Primary Module Functions ###
    download=T/F    : Download EnsEMBL databases [False]
    makeuniprot=T/F : Whether to generate an Ensembl.dat file of UniProt entries for species [False]
    ensloci=T/F     : Create EnsEMBL datasets "reduced by loci" [False]
    enspep=T/F      : Create full gnspacc EnsEMBL peptide datasets [False]
    hgncmap=FILE    : File to be used for HGNC ID mapping []
    resume=X        : Species or species code to pickup run from [None]
    sections=LIST   : List of Ensembl sections to use for run (else All) []
    speclist=LIST   : List of species to use for run (else All) []
    chromspec=LIST  : List of species codes to download chromosomes for [HUMAN,DROME,CAEEL,YEAST,MOUSE,DANRE,CHICK,XENTR]
    speedskip=T/F   : Whether to assume download is fine if pep.all/cdna.all/dna.toplevel file found [True]
    ### Advanced UniProt Mapping Options ###
    mapstat=X       : GABLAM Stat to use for mapping assessment (ID/Sim/Len) [ID]
    automap=X       : Minimum value of mapstat for mapping to occur [80.0]
    unispec=FILE    : Alternative UniProt species file [None]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### EnsGO Options ###    
    ensgo=LIST      : List of species codes to make EnsGO Datasets for []
    mingo=X         : Minumum number of genes to output GO category [0]
    obsgo=T/F       : Whether to include obselete terms [False]
    splicego=T/F    : Whether to include all splice variants (EnsEMBL peptides) in GO datasets [False]
    goids=FILE      : File containing GO IDs [GO.terms_ids_obs]
    goevidence=LIST : List of acceptable GO evidence codes. (Will use all if blank.) []
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### EnsDat Options ###
    ensdat=LIST     : Perform EnsDat construction of predicted UniProt data for the species listed []
    tmhmm=FILE      : Path to TMHMM program [/home/richard/Bioware/TMHMM2.0c/bin/tmhmm]
    signalp=FILE    : Path to SIGNALP program [/home/richard/Bioware/signalp-3.0/signalp]
    hmmerpath=PATH  : Path for hmmer files [/home/richard/Bioware/hmmer-2.3.2/src/]
    pfam=FILE       : Path to PFam LS file [/home/richard/Databases/PFam/Pfam_ls]
    datpickup=FILE  : Text file containing names of proteins already processed (skip and append) [ensdat.txt]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### System Parameters ###
    enspath=PATH    : Path to EnsEMBL file [EnsEMBL/]
    unipath=PATH    : Path to UniProt files [enspath=PATH/uniprot/]
    specsleep=X     : Sleep for X seconds between species downloads [60]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Uses general modules: copy, glob, urllib, os, string, sys, time
Uses RJE modules: rje, rje_uniprot
Other modules needed: rje_sequence
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, urllib, os, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_db, rje_forker, rje_seq, rje_seqlist, rje_taxonomy, rje_tm, rje_uniprot, seqmapper, rje_zen
import rje_hmm_V1 as rje_hmm
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial compilation.
    # 1.0 - Initial working version with download and EnsLoci functions.
    # 1.1 - Added reformatDB() method to be called by rje_dbase
    # 1.2 - Fixed known-ccds bug
    # 1.3 - SwissProt bug still remained. Fixed and added speclist=LIST
    # 1.4 - Improved EnsLoci mapping to ignore Xs and use SeqMapper for better match to SwissProt
    # 1.5 - Fixed bug that results in multiple occurrences of some sequence names
    # 1.6 - Added crap-sequence catching
    # 1.7 - Added GO dataset partitioning
    # 1.8 - Modification to the GO dataset generation, and GO adaptation for use with PINGU.
    # 1.9 - Added the possibility to restrict GO data to certain evidence codes
    # 2.0 - Added new EnsDat functionality. (Now part of "UniFake" module.)
    # 2.1 - Made changes in line with new EnsEMBL setup.
    # 2.2 - Made changes in line with new EnsEMBL setup. Again. Grrrr. Stop changing things, EnsEMBL!
    # 2.3 - Added new species, metlist and option to download chromosome sequences.
    # 2.4 - Modified to allow HGNC evidence for Human EnsEMBL.
    # 2.5 - Modified to allow additional species-specific evidence in *.map.tdt file. (Mouse, Yeast, Zebrafish)
    # 2.6 - Added additional EnsEMBL sites to metazoa: fungi, plants, protists
    # 2.7 - Updated for new EnsEMBL format with ID mapping in gene.txt file.
    # 2.8 - Bug fixes for updated EnsEMBL release.
    # 2.9 - Reduced DNA chromosome downloads. Updated some species data. Added "known_by_projection" handling.
    # 2.10- Miscellaneous fixes.
    # 2.11- Added rje_taxonomy and makeuniprot=T/F. Removed metlist. Moved release and species data extraction.
    # 2.12- Changed chromspec to enable downloads of all species but also download toplevel files, not chromosomes.
    # 2.13- Added speedskip=T/F [True] that will skip when pep.all, cdna.all and dna.toplevel are found.
    # 2.14- Add enspep=T/F      : Create full gnspacc EnsEMBL peptide datasets [False]
    # 2.15.0 - Added capacity to download/process a section of Ensembl with speclist=LIST.
    # 2.15.1 - Improved error handling for too many FTP connections: still need to fix problem!
    # 2.15.2 - Trying to improve speed of Uniprot parsing for EnsLoci.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add checking date of files vs *.loci.fas and force=T/F option.
    # [Y] : Add specsleep option?
    # [Y] : Add resume=X option - starts from given species or species code: rje.sortKeys(enspec)[X:]
    # [Y] : Make it work with known only
    # [Y] : Add storage and output of stats on types of genes etc. (append=T/F) & .clear()
    # [ ] : Add incorporation of actual UniProt data into EnsDat function.
    # [ ] : Add PUTATIVE gene type in addition to known and novel?
    # [ ] : Add EnsEMBL Bacteria! (Like other esite but much more and with _collection layer.)
    # [Y] : Split of subsections into subdirectories
    # [Y] : Add Species CODES for new species as they are read in and output to file for easy additonal to module.
    # [ ] : Update species lists and common names. Check some of the previously missing species codes.
    # [ ] : Makes a reduced set of EnsEMBL species for making alignments etc. Or make Taxa-specific EnsEMBL data?
    # [Y] : Method to check species codes against UniProt.
    # [ ] : Get working with new SeqMapper.
    # [ ] : Add REST access: http://beta.rest.ensembl.org/sequence/id/ENST00000299130?
    # [Y] : makeuniprot=T/F : Whether to generate an Ensembl.dat file of UniProt entries for species [False]
    # [Y] : Incorporate rje_taxonomy for spcode mapping and TaxID output for UniProt download.
    # [ ] : Upgrade to V3.0 with new RJE_Object and fully replace species dictionaries with rje_taxonomy.
    # [ ] : Use rje_seqlist for speed.
    # [Y] : Direct downloads into subdirectory for data management.
    # [ ] : Tidy common names from species and then use to over-ride rje_taxonomy extraction.
    # [ ] : Contemplate adding to EnsTaxID rather than replacing it from the commandline.
    # [Y] : Add extraction of DB Xref table from Uniprot following download for EnsLoci construction.
    # [Y] : Switch self.dict['Stats'] to be a db table. For EnsLoci, EnsEMBL objects can use same DB object.
    # [?] : Split Uniprot download into sections?
    # [Y] : Get rid of stable_gene_id: not in current release and not required.
    # [Y] : Use forking to generate EnsLoci in parallel.
    # [ ] : Add switch for downloading selection of toplevel, chromosomes, masked, assembly etc.
    #       |-- See: ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/README
    # [Y] : Add speedskip=T/F that will skip when pep.all, cdna.all and dna.toplevel are found.
    # [ ] : Add an rsync-based method, e.g. rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_fasta/*/dna/*.dna.toplevel.fa.gz.
    #       |-- Might only work for main species? ensemblgenomes address does not work?
    # [ ] : Add option to tidy up genome sequence name.
    # [ ] : Add option to download but not unzip genome files.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('RJE_EnsEMBL', '2.15.2', 'April 2015', '2007')
    description = 'EnsEMBL Processing/Manipulation Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is under development and may contain bugs!',
                'If things go badly wrong, check whether EnsEMBL has made changes.',
                rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CONSTANTS                                                                                               #
#########################################################################################################################
### Dictionary of species linked by underscore as in EnsEMBL file name and corresponding UniProt species code ###
enspec = {'Acyrthosiphon_pisum':'ACYPI',
          'Aedes_aegypti':'AEDAE',
          'Ailuropoda_melanoleuca':'AILME',
          'Anas_platyrhynchos':'ANAPP',
          'Anolis_carolinensis':'ANOCA',
          'Anopheles_gambiae':'ANOGA',
          'Apis_mellifera':'APIME',
          'Arabidopsis_lyrata':'ARALY',
          'Arabidopsis_thaliana':'ARATH',
          'Aspergillus_clavatus':'ASPCL',
          'Aspergillus_flavus':'ASPFL',
          'Aspergillus_fumigatus':'ASPFU',	# Not in UniProt
          'Aspergillus_nidulans':'ASPNI',	# Not in UniProt
          'Aspergillus_niger':'ASPNC',	# Not ASPNI
          'Aspergillus_oryzae':'ASPOR',
          'Aspergillus_terreus':'ASPTE',
          'Bos_taurus':'BOVIN',
          'Brachypodium_distachyon':'BRADI',
          'Caenorhabditis_elegans':'CAEEL',
          'Caenorhabditis_brenneri':'CAEBE',
          'Caenorhabditis_briggsae':'CAEBR',
          'Caenorhabditis_japonica':'CAEJA',
          'Caenorhabditis_remanei':'CAERE', 
          'Callithrix_jacchus':'CALJA',
          'Canis_familiaris':'CANFA',
          'Cavia_porcellus':'CAVPO',
          'Choloepus_hoffmanni':'CHOHO',
          'Chlamydomonas_reinhardtii':'CHLRE',
          'Ciona_intestinalis':'CIOIN',
          'Ciona_savignyi':'CIOSA',
          'Culex_quinquefasciatus':'CULQU', 
          'Danio_rerio':'DANRE',
          'Daphnia_pulex':'DAPPU',
          'Dasypus_novemcinctus':'DASNO',
          'Dictyostelium_discoideum':'DICDI',
          'Dipodomys_ordii':'DIPOR',
          'Drosophila_melanogaster':'DROME',
          'Drosophila_ananassae':'DROAN',
          'Drosophila_erecta':'DROER',
          'Drosophila_grimshawi':'DROGR',
          'Drosophila_mojavensis':'DROMO',
          'Drosophila_persimilis':'DROPE',
          'Drosophila_pseudoobscura':'DROPS',
          'Drosophila_sechellia':'DROSE',
          'Drosophila_simulans':'DROSI',
          'Drosophila_virilis':'DROVI',
          'Drosophila_willistoni':'DROWI',
          'Drosophila_yakuba':'DROYA',
          'Echinops_telfairi':'ECHTE',
          'Equus_caballus':'EQUPR',
          'Erinaceus_europaeus':'ERIEU',
          'Felis_catus':'FELCA',
          'Fusarium_oxysporum':'FUSOX',
          'Gallus_gallus':'CHICK',
          'Gasterosteus_aculeatus':'GASAC',
          'Gibberella_moniliformis':'GIBMO',
          'Gibberella_zeae':'GIBZA',	# Not GIBZE
          'Glycine_max':'SOYBN',	# Not GLYMA
          'Gorilla_gorilla':'GORGO',
          'Homo_sapiens':'HUMAN',
          'Ixodes_scapularis':'IXOSC',
          'Leishmania_major':'LEIMA',
          'Loxodonta_africana':'LOXAF',
          'Macaca_mulatta':'MACMU',
          'Macropus_eugenii':'MACEU',
          'Meleagris_gallopavo':'MELGA',
          'Microcebus_murinus':'MICMU',
          'Monodelphis_domestica':'MONDO',
          'Mus_musculus':'MOUSE',
          'Myotis_lucifugus':'MYOLU',
          'Mycosphaerella_graminicola':'MYCGR',
          'Nectria_haematococca':'NECH7',	# Not NECHA
          'Neosartorya_fischeri':'NEOFI',
          'Neurospora_crassa':'NEUCR',
          'Nematostella_vectensis':'NEMVE',
          'Nomascus_leucogenys':'NOMLE',
          'Ochotona_princeps':'OCHPR',
          'Ornithorhynchus_anatinus':'ORNAN',
          'Oryctolagus_cuniculus':'RABIT',
          'Oryza_glaberrima':'ORYGL',
          'Oryza_indica':'ORYIN',	# Not in UniProt
          'Oryza_sativa':'ORYSA',
          'Oryzias_latipes':'ORYLA',
          'Otolemur_garnettii':'OTOGA',
          'Ovis_aries':'SHEEP',
          'Pan_troglodytes':'PANTR',
          'Pediculus_humanus':'PEDHC',
          'Petromyzon_marinus':'PETMA',
          'Phaeodactylum_tricornutum':'PHATC',	# Not PHATR
          'Phaeosphaeria_nodorum':'PHAND',	# Not PHANO
          'Physcomitrella_patens':'PHYPA',	# Not in UniProt
          'Phytophthora_infestans':'PHYIN',
          'Phytophthora_ramorum':'PHYRM',	# Not PHYRA
          'Phytophthora_sojae':'PHYSO',	# Not in UniProt
          'Plasmodium_berghei':'PLABA',	# Not PLABE
          'Plasmodium_chabaudi':'PLACH',
          'Plasmodium_falciparum':'PLAF1',	# Not PLAFA
          'Plasmodium_knowlesi':'PLAKH',	# Not PLAKN
          'Plasmodium_vivax':'PLAVB',	# Not PLAVI
          'Pongo_abelii':'PONAB',
          'Pongo_pygmaeus':'PONPY', # Renamed! #
          'Populus_trichocarpa':'POPTR',
          'Pristionchus_pacificus':'PRIPA',
          'Procavia_capensis':'PROCA',
          'Pteropus_vampyrus':'PTEVA',
          'Puccinia_graministritici':'PUCGR',	# Not in UniProt
          'Puccinia_triticina':'PUCTR',	# Not in UniProt
          'Pythium_ultimum':'PYTUL',	# Not in UniProt
          'Rattus_norvegicus':'RAT',
          'Saccharomyces_cerevisiae':'YEAST',
          'Sarcophilus_harrisii':'SARHA',
          'Schizosaccharomyces_pombe':'SCHPM',	# Not SCHPO
          'Schistosoma_mansoni':'SCHMA',
          'Selaginella_moellendorffii':'SELML',	# Not SELMO
          'Sorghum_bicolor':'SORBI',
          'Sorex_araneus':'SORAR',
          'Spermophilus_tridecemlineatus':'SPETR',
          'Strongylocentrotus_purpuratus':'STRPU',
          'Sus_scrofa':'PIG',
          'Taeniopygia_guttata':'TAEGU',
          'Takifugu_rubripes':'TAKRU',
          'Tarsius_syrichta':'TARSY',
          'Tetraodon_nigroviridis':'TETNG',
          'Thalassiosira_pseudonana':'THAPS',
          'Trypanosoma_brucei':'TRYBR',	# Not in UniProt
          'Tuber_melanosporum':'TUBMM',	# Not TUBME
          'Trichoplax_adhaerens':'TRIAD',
          'Tupaia_belangeri':'TUPGB',
          'Tursiops_truncatus':'TURTR',
          'Ustilago_maydis':'USTMA',
          'Vitis_vinifera':'VITVI',
          'Vicugna_pacos':'LAMPA',
          'Xenopus_tropicalis':'XENTR',
          'Zea_mays':'MAIZE',	# Not ZEAMA
          ### Check and update below ###
          "Amphimedon_queenslandica":"AMPQUI",    # Does not actually exist yet
          "Ashbya_gossypii":"ASHGO",              # read for "Ashbya_gossypii".
          "Atta_cephalotes":"ATTCE",
          "Bombyx_mori":"BOMMO",    # read for "Bombyx_mori".
          "Cyanidioschyzon_merolae":"CYAME",   # read for "Cyanidioschyzon_merolae".
          "Gadus_morhua":"GADMO",              # read for "Gadus_morhua".
          "Latimeria_chalumnae":"LATCH",       # read for "Latimeria_chalumnae".
          "Tribolium_castaneum":"TRICA",       # read for "Tribolium_castaneum".
          'Albugo_laibachii':'ALBLA',	# Not in UniProt
          'Botryotinia_fuckeliana':'BOTF4',	# Not BOTFU
          'Brassica_rapa':'BRARA',
          'Danaus_plexippus':'DANPL',
          'Entamoeba_histolytica':'ENTHI',
          'Gaeumannomyces_graminis':'GAEGR',	# UniProt
          'Heliconius_melpomene':'HELME',	# Not in UniProt
          'Hyaloperonospora_arabidopsidis':'HYAARA',	# Not in UniProt, Not HYARA
          'Magnaporthe_oryzae':'MAGO7',	# Not MAGOR
          'Magnaporthe_poae':'MAGPO',
          'Oreochromis_niloticus':'ORENI',
          'Oryza_brachyantha':'ORYBR',	# Not in UniProt
          'Puccinia_graminis':'PUCGR',
          'Sclerotinia_sclerotiorum':'SCLS1',	# Not SCLSC
          'Setaria_italica':'SETIT',
          'Solanum_lycopersicum':'SOLLC',	# Not SOLLY
          'Tetrahymena_thermophila':'TETTH',
          'Toxoplasma_gondii':'TOXGO',
          'Trichinella_spiralis':'TRISP'
          }
enscom = {'CAEBE':'C. brenneri',
          'ACYPI':'Pea aphid',
          'AEDAE':'Mosquito',
          'AILME':'Panda',
          'ANAPP':'Duck',
          'ANOCA':'Anole lizard',
          'ANOGA':'Anopheles',
          'APIME':'Honeybee',
          'BOVIN':'Cow',
          'CAEBR':'C. briggsae',
          'CAEEL':'Nematode',
          'CAEJA':'C. japonica',
          'CAERE':'C. remanei',
          'CALJA':'Marmoset',
          'CANFA':'Dog',
          'CAVPO':'Guinea Pig',
          'CHICK':'Chicken',
          'CHOHO':'Sloth',
          'CIOIN':'Ciona Int.',
          'CIOSA':'Ciona Sav.',
          'CULQU':'Southern house mosquito',
          'DANRE':'Zebrafish','BRARE':'Zebrafish',
          'DASNO':'Armadillo',
          'DIPOR':'Kangaroo rat',
          'DROAN':'Drosophila ananassae',
          'DROER':'Drosophila erecta',
          'DROGR':'Drosophila grimshawi',
          'DROME':'Fruitfly',
          'DROMO':'Drosophila mojavensis',
          'DROPE':'Drosophila persimilis',
          'DROPS':'Drosophila pseudoobscura',
          'DROSE':'Drosophila sechellia',
          'DROSI':'Drosophila simulans',
          'DROVI':'Drosophila virilis',
          'DROWI':'Drosophila willistoni',
          'DROYA':'Drosophila yakuba',
          'ECHTE':'Tenrec',
          'EQUPR':'Horse',
          'ERIEU':'Hedgehog',
          'FELCA':'Cat',
          'FUGRU':'Fugu',
          'GASAC':'Stickleback',
          'GORGO':'Gorilla',
          'HUMAN':'Human',
          'IXOSC':'Deer tick',
          'LAMPA':'Alpaca',
          'LOXAF':'Elephant',
          'MACEU':'Wallaby',
          'MACMU':'Macaque',
          'MELGA':'Turkey',
          'MICMU':'Mouse lemur',
          'MONDO':'Opossum',
          'MOUSE':'Mouse',
          'MYOLU':'Bat',
          'NEMVE':'Sea anemone',
          'NOMLE':'Gibbon', # Northern white-cheeked gibbon
          'OCHPR':'Pika',
          'ORYLA':'Ricefish',
          'ORNAN':'Platypus',
          'OTOGA':'Bushbaby',
          'PANTR':'Chimp',
          'PEDHC':'Human body louse',
          'PETMA':'Lamprey',
          'PIG':'Pig',
          'PONAB':'Orangutan',
          'PONPY':'Orangutan',
          'PROCA':'Hyrax',
          'PRIPA':'Pristionchus pacificus',
          'PTEVA':'Megabat',
          'RAT':'Rat',
          'RABIT':'Rabbit',
          'SCHMA':'Blood fluke',
          'SHEEP':'Sheep',
          'SORAR':'Shrew',
          'SPETR':'Squirrel',
          'STRPU':'Sea urchin',
          'TAEGU':'Zebra finch',
          'TARSY':'Tarsier',
          'TETNG':'Pufferfish',
          'TRIAD':'Trichoplax adhaerens',
          'TUPGB':'Treeshrew',
          'TURTR':'Dolphin',
          'XENTR':'Xenopus',
          'YEAST':'Yeast',
          "AMPQUI":"Amphimedon queenslandica",    # Code Does not actually exist yet
          "ASHGO":'Ashbya gossypii (Yeast)',              # read for "Ashbya_gossypii".
          "ATTCE":'Leafcutter ant',
          "BOMMO":'Silk moth',    # read for "Bombyx_mori".
          "CYAME":'Red alga',   # read for "Cyanidioschyzon_merolae".
          "GADMO":'Atlantic cod',              # read for "Gadus_morhua".
          "LATCH":'West Indian ocean coelacanth',       # read for "Latimeria_chalumnae".
          "TRICA":'Red flour beetle',       # read for "Tribolium_castaneum".
          'ARALY':'A. lyrata',	# Arabidopsis_lyrata'
          'ARATH':'A. thaliana',	# Arabidopsis_thaliana'
          'ASPCL':'A. clavatus',	# Aspergillus_clavatus'
          'ASPFL':'A. flavus',	# Aspergillus_flavus'
          'ASPFU':'A. fumigatus',	# Aspergillus_fumigatus'	# Not in UniProt
          'ASPNI':'A. nidulans',	# Aspergillus_nidulans'	# Not in UniProt
          'ASPNC':'A. niger',	# Aspergillus_niger'	# Not ASPNI
          'ASPOR':'A. oryzae',	# Aspergillus_oryzae'
          'ASPTE':'A. terreus',	# Aspergillus_terreus'
          'BRADI':'B. distachyon',	# Brachypodium_distachyon'
          'CHLRE':'C. reinhardtii',	# Chlamydomonas_reinhardtii'
          'DAPPU':'D. pulex',	# Daphnia_pulex'
          'DICDI':'D. discoideum',	# Dictyostelium_discoideum'
          'FUSOX':'F. oxysporum',	# Fusarium_oxysporum'
          'GIBMO':'G. moniliformis',	# Gibberella_moniliformis'
          'GIBZA':'G. zeae',	# Gibberella_zeae'	# Not GIBZE
          'SOYBN':'Soybean',	# Glycine_max'	# Not GLYMA
          'LEIMA':'L. major',	# Leishmania_major'
          'MYCGR':'M. graminicola',	# Mycosphaerella_graminicola'
          'NECH7':'N. haematococca',	# Nectria_haematococca'	# Not NECHA
          'NEOFI':'N. fischeri',	# Neosartorya_fischeri'
          'NEUCR':'N. crassa',	# Neurospora_crassa'
          'ORYGL':'O. glaberrima',	# Oryza_glaberrima'
          'ORYIN':'O. indica',	# Oryza_indica'	# Not in UniProt
          'ORYSA':'O. sativa',	# Oryza_sativa'
          'PHATC':'P. tricornutum',	# Phaeodactylum_tricornutum'	# Not PHATR
          'PHAND':'P. nodorum',	# Phaeosphaeria_nodorum'	# Not PHANO
          'PHYPA':'P. patens',	# Physcomitrella_patens'	# Not in UniProt
          'PHYIN':'P. infestans',	# Phytophthora_infestans'
          'PHYRM':'P. ramorum',	# Phytophthora_ramorum'	# Not PHYRA
          'PHYSO':'P. sojae',	# Phytophthora_sojae'	# Not in UniProt
          'PLABA':'P. berghei',	# Plasmodium_berghei'	# Not PLABE
          'PLACH':'P. chabaudi',	# Plasmodium_chabaudi'
          'PLAF1':'P. falciparum',	# Plasmodium_falciparum'	# Not PLAFA
          'PLAKH':'P. knowlesi',	# Plasmodium_knowlesi'	# Not PLAKN
          'PLAVB':'P. vivax',	# Plasmodium_vivax'	# Not PLAVI
          'POPTR':'P. trichocarpa',	# Populus_trichocarpa'
          'PUCTR':'P. triticina',	# Puccinia_triticina'	# Not in UniProt
          'PYTUL':'P. ultimum',	# Pythium_ultimum'	# Not in UniProt
          'SARHA':'S. harrisii',	# Sarcophilus_harrisii'
          'SCHPM':'S. pombe',	# Schizosaccharomyces_pombe'	# Not SCHPO
          'SELML':'S. moellendorffii',	# Selaginella_moellendorffii'	# Not SELMO
          'SORBI':'S. bicolor',	# Sorghum_bicolor'
          'THAPS':'T. pseudonana',	# Thalassiosira_pseudonana'
          'TRYBR':'T. brucei',	# Trypanosoma_brucei'	# Not in UniProt
          'TUBMM':'T. melanosporum',	# Tuber_melanosporum'	# Not TUBME
          'USTMA':'U. maydis',	# Ustilago_maydis'
          'VITVI':'Grape',	# Vitis_vinifera'
          'MAIZE':'Maize',	# Zea_mays'	# Not ZEAMA
          'ALBLA':'A. laibachii',	# Albugo_laibachii'	# Not in UniProt
          'BOTF4':'B. fuckeliana',	# Botryotinia_fuckeliana'	# Not BOTFU
          'BRARA':'B. rapa',	# Brassica_rapa'
          'DANPL':'D. plexippus',	# Danaus_plexippus'
          'ENTHI':'E. histolytica',	# Entamoeba_histolytica'
          'GAEGR':'G. graminis',	# Gaeumannomyces_graminis'	# Not in UniProt
          'HELME':'H. melpomene',	# Heliconius_melpomene'	# Not in UniProt
          'HYAAR':'H. arabidopsidis',	# Hyaloperonospora_arabidopsidis'	# Not in UniProt
          'MAGO7':'M. oryzae',	# Magnaporthe_oryzae'	# Not MAGOR
          'MAGPO':'M. poae',	# Magnaporthe_poae'
          'ORENI':'O. niloticus',	# Oreochromis_niloticus'
          'ORYBR':'O. brachyantha',	# Oryza_brachyantha'	# Not in UniProt
          'PUCGR':'P. graminis',	# Puccinia_graminis'
          'SCLS1':'S. sclerotiorum',	# Sclerotinia_sclerotiorum'	# Not SCLSC
          'SETIT':'S. italica',	# Setaria_italica'
          'SOLLC':'S. lycopersicum',	# Solanum_lycopersicum'	# Not SOLLY
          'TETTH':'T. thermophila',	# Tetrahymena_thermophila'
          'TOXGO':'Toxoplasma',	# Toxoplasma_gondii'
          'TRISP':'T. spiralis'	# Trichinella_spiralis'
          }
add_common = '''
'Aspergillus_fumigatus':'ASPFU',	# Not in UniProt
'Aspergillus_nidulans':'ASPNI',	# Not in UniProt
'Aspergillus_niger':'ASPNC',	# Not ASPNI
'Gibberella_zeae':'GIBZA',	# Not GIBZE
'Glycine_max':'SOYBN',	# Not GLYMA
'Nectria_haematococca':'NECH7',	# Not NECHA
'Oryza_indica':'ORYIN',	# Not in UniProt
'Phaeodactylum_tricornutum':'PHATC',	# Not PHATR
'Phaeosphaeria_nodorum':'PHAND',	# Not PHANO
'Physcomitrella_patens':'PHYPA',	# Not in UniProt
'Phytophthora_ramorum':'PHYRM',	# Not PHYRA
'Phytophthora_sojae':'PHYSO',	# Not in UniProt
'Plasmodium_berghei':'PLABA',	# Not PLABE
'Plasmodium_falciparum':'PLAF1',	# Not PLAFA
'Plasmodium_knowlesi':'PLAKH',	# Not PLAKN
'Plasmodium_vivax':'PLAVB',	# Not PLAVI
'Puccinia_graministritici':'PUCGR',	# Not in UniProt
'Puccinia_triticina':'PUCTR',	# Not in UniProt
'Pythium_ultimum':'PYTUL',	# Not in UniProt
'Schizosaccharomyces_pombe':'SCHPM',	# Not SCHPO
'Selaginella_moellendorffii':'SELML',	# Not SELMO
'Trypanosoma_brucei':'TRYBR',	# Not in UniProt
'Tuber_melanosporum':'TUBMM',	# Not TUBME
'Zea_mays':'MAIZE',	# Not ZEAMA
'Albugo_laibachii':'ALBLA',	# Not in UniProt
'Botryotinia_fuckeliana':'BOTF4',	# Not BOTFU
'Gaeumannomyces_graminis':'GAEGR',	# Not in UniProt
'Heliconius_melpomene':'HELME',	# Not in UniProt
'Hyaloperonospora_arabidopsidis':'HYAAR',	# Not in UniProt
'Magnaporthe_oryzae':'MAGO7',	# Not MAGOR
'Oryza_brachyantha':'ORYBR',	# Not in UniProt
'Sclerotinia_sclerotiorum':'SCLS1',	# Not SCLSC
'Solanum_lycopersicum':'SOLLC',	# Not SOLLY
'''
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: EnsEMBL Class                                                                                          #
#########################################################################################################################
class EnsEMBL(rje.RJE_Object):     
    '''
    EnsEMBL Processing Class. Author: Rich Edwards (2006).

    Info:str
    - Name = EnsEMBL Species 
    - EnsPath = Path to EnsEMBL file [EnsEMBL/]
    - GOIDs = File containing GO IDs [GO.terms_ids_obs]
    - HGNCMap = File to be used for HGNC ID mapping []
    - UniPath = Path to UniProt files [UniProt/]
    - UniSpec = Alternative UniProt species file
    - Resume = Species or species code to pickup run from [None]
    - HMMerPath = Path for hmmer files [/home/richard/Bioware/hmmer-2.3.2/src/]    
    - TMHMM = Path to TMHMM program [/home/richard/Bioware/TMHMM2.0c/bin/tmhmm]
    - SignalP = Path to SIGNALP program [/home/richard/Bioware/signalp-3.0/signalp]
    - PFam = Path to PFam LS file [/home/richard/Databases/PFam/Pfam_ls]
    - DatPickup = Text file containing names of proteins already processed (skip and append) [ensdat.txt]
    
    Opt:boolean
    - DownLoad = Download EnsEMBL databases [False]
    - EnsLoci = Create EnsEMBL datasets "reduced by loci" [False]
    - EnsPep = Create full gnspacc EnsEMBL peptide datasets [False]
    - MakeUniprot = Whether to generate an Ensembl.dat file of UniProt entries for species [False]
    - ObsGO = Whether to include obselete terms [False]
    - SpeedSkip = Whether to assume download is fine if pep.all/cdna.all/dna.toplevel file found [True]
    - SpliceGO = Whether to include all splice variants (EnsEMBL peptides) in GO datasets [False]

    Stat:numeric
    - SpecSleep = Sleep for X seconds between species downloads [60]
    - MinGO = Minumum number of genes to output GO category [3]
    
    List:list
    - ChromSpec = List of species codes to download chromosomes for [HUMAN,DROME,CAEEL,YEAST]
    - EnsDat = Perform EnsDat construction of predicted UniProt data for the species listed []
    - EnsGO = List of species codes to make EnsGO Datasets for []
    - GOEvidence = List of acceptable GO evidence codes. (Will use all if blank.) []
    - Sections = List of Ensembl sections to use for run (else All) []
    - SpecList = List of species to use for run (else All) []

    Dict:dictionary
    - AccName = Dictionary of {AccNum:New short name} from UniProt
    - NewAccName = Dictionary of {ProtID:New short name} from UniProt (avoids double naming)
    - AccSeq = Dictionary of {AccNum:Sequence} from UniProt
    - EnsLoci = Dictionary of {ProtID:Full Name} for reduced dataset
    - EnsSpec = Dictionary like global enspec {Genus_species:code} but can be updated
    - EnsTaxID = Dictionary of {Genus_species:taxid} to over-ride value read from Taxonomy files. {'Saccharomyces_cerevisiae':'559292'}
    - GeneAcc = Dictionary of {GeneNum:AccNum}
    - GeneDesc = Dictionary of {GeneNum:Desc}
    - GeneMap = Dictionary of {GeneID:GeneNum}
    - GO = Dictionary of GO data read in during EnsGO
    - Loci = Dictionary of {Type:{GeneID:ProtIDList}}
    - ProtLen = Dictionary of {Type:{ProtID:Length}}
    - ProtSeq = Dictionary of {ProtID:Sequence} for those ProtIDs with AccNums in AccSeq
    - Stats = Stores a number of stats for output into EnsEMBL/ens_loci.tdt

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - Taxonomy = rje_taxonomy.Taxonomy object
    - UniProt = rje_uniprot.UniProt object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','EnsPath','UniPath','Resume','GOIDs','TMHMM','SignalP','PFam','DatPickup','HMMerPath',
                         'HGNCMap','UniSpec']
        self.optlist = ['Download','EnsLoci','EnsPep','ObsGO','SpliceGO','MakeUniprot','SpeedSkip']
        self.statlist = ['SpecSleep','MinGO']
        self.listlist = ['SpecList','EnsGO','GOEvidence','EnsDat','ChromSpec','Sections']
        self.dictlist = ['AccName','NewAccName','AccSeq','EnsLoci','GeneAcc','GeneDesc','GeneMap','Loci','ProtLen',
                         'ProtSeq','Stats','EnsSpec','GO']
        self.objlist = ['UniProt']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0,obj=None,setlist=True,setdict=True)
        self.setInfo({'EnsPath':rje.makePath('EnsEMBL/'),'UniPath':'','GOIDs':'GO.terms_ids_obs',
                      'HMMerPath':'/home/richard/Bioware/hmmer-2.3.2/src/',
                      'TMHMM':'/home/richard/Bioware/TMHMM2.0c/bin/tmhmm',
                      'SignalP':'/home/richard/Bioware/signalp-3.0/signalp',
                      'PFam':'/home/richard/Databases/PFam/Pfam_ls','DatPickup':'ensdat.txt'})
        self.setStat({'MinGO':0,'SpecSleep':60})
        self.setOpt({'SpeedSkip':True})
        self.list['ChromSpec'] = rje.split('HUMAN,DROME,CAEEL,YEAST,MOUSE,DANRE,CHICK,XENTR',',')
        self.list['SpecList'] = []  #rje.sortKeys(enspec)
        self.dict['EnsTaxID'] = {'Saccharomyces_cerevisiae':'559292','Oryza_indica':'39946','Gorilla_gorilla':'9595',
                                 'Pyrenophora_triticirepentis':'426418'}
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['Taxonomy'] = rje_taxonomy.Taxonomy(self.log,self.cmd_list)
        self.log.no_suppression += ['Missing Taxon',"Invented SpCode"]
        self.log.warnings += ['SpCode Missing TaxID']
        self._setForkAttributes()   # Used by Uniprot
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
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options ###
                self._cmdReadList(cmd,'info',['Resume'])
                self._cmdReadList(cmd,'path',['EnsPath','UniPath','HMMerPath'])
                self._cmdReadList(cmd,'file',['TMHMM','SignalP','PFam','DatPickup','HGNCMap','UniSpec'])
                self._cmdReadList(cmd,'int',['SpecSleep','MinGO'])
                self._cmdReadList(cmd,'opt',['Download','EnsLoci','EnsPep','ObsGO','SpliceGO','MakeUniprot','SpeedSkip'])
                self._cmdReadList(cmd,'list',['SpecList','EnsGO','GOEvidence','EnsDat','ChromSpec','Sections'])
                self._cmdReadList(cmd,'cdict',['EnsTaxID'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        if not self.getStr('UniPath'): self.setStr({'UniPath':rje.makePath('%suniprot/' % self.getStr('EnsPath'))})
        self.list['ChromSpec'] = rje.split(rje.join(self.list['ChromSpec']).upper())
#########################################################################################################################
    ### <2> ### Main Run Methods                                                                                        #
#########################################################################################################################
    def run(self):  ### Downloads and reformats if desired
        '''Downloads and reformats if desired.'''
        ### ~ General Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.obj['Taxonomy'].setupSourceData(): self.obj['Taxonomy'] = None
        if not self.getStrLC('Resume'): self.setStr({'Resume':''})
        if self.getBool('EnsLoci'): self.summaryStats(setup=True)
        if not self.getStrLC('UniSpec'): self.info['UniSpec'] = self.info['UniPath'] + 'uniprot.spec.tdt'
        else: self.warnLog('unispec=FILE discontinued in favour of rje_taxonomy.')
        ensforker = EnsForker(self.log,self.cmd_list)
        #self.debug(self.getInt('Forks'))

        ### ~ EnsEMBL Downloads and EnsLoci processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('Download'):    #!# Remove EnsLoci from download
            if not self.downloadEnsEMBL(self.list['SpecList']): return False
        if self.getBool('EnsPep'): self.makeEnsPepFasta()
        if self.getBool('EnsLoci') or self.getBool('MakeUniprot'):
            rdb = self.makeEnsemblUniprot()
        if self.getBool('EnsLoci') and not self.getBool('Download'):
            if not self.list['SpecList']: self.list['SpecList'] = rdb.dataKeys()
            if self.getInt('Forks') > 1:
                ensforker.list['ToFork'] = []
                ensforker.list['ResFile'] = ['ens_loci.tdt']
                ensforker.dict['SiteSpec'] = self.dict['SiteSpec']
                ensforker.obj['DB'] = self.obj['DB']
            for spec in self.list['SpecList']:
                species = self.getSpecies(spec)
                if self.getStr('Resume'):
                    if self.getStr('Resume') in [species,self.speciesCode(species)]: self.setStr({'Resume':''})
                    else: continue
                if self.getInt('Forks') > 1:
                    if not self.force() and rje.exists(self.ensLociFas(species)):
                        self.printLog('#LOCI','%s found! (force=F)' % self.ensLociFas(species))
                        continue
                    ensforker.list['ToFork'].append(species)
                else: self.ensLoci(species)
            if self.getInt('Forks') > 1:
                self.debug(ensforker.list['ToFork'])
                ensforker.run()
                #!# Replace rdb with reloaded data? #!#

        ### ~ Additional functions and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.ensGO()
        self.ensDat()
        return True
#########################################################################################################################
    def speciesCode(self,getspecies=None):  ### Returns species code from Species = self.info['Name']
        '''Returns species code from Species = self.info['Name'].'''
        ### ~ [1] Try to pull from dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not getspecies: getspecies = self.info['Name']
        if getspecies in self.dict['EnsSpec']: return self.dict['EnsSpec'][getspecies]
        elif self.db('release') and self.db('release').data(getspecies) and self.db('release').data(getspecies)['spcode']:
            return self.db('release').data(getspecies)['spcode']
        elif self.obj['Taxonomy']:
            spcode = self.obj['Taxonomy'].getSpCode(getspecies,invent=False)
            self.dict['EnsSpec'][getspecies] = spcode
            if not spcode:
                spcode = self.obj['Taxonomy'].getSpCode(getspecies,invent=True)
                self.dict['EnsSpec'][getspecies] = spcode
                if self.i() < 1: return spcode
                else: self.dict['EnsSpec'][getspecies] = code = rje.choice('Species code for "%s"?' % getspecies,spcode,True)
                self.printLog('\n#CODE','Species code "%s" created for "%s".' % (code,getspecies))
                open('rje_ensembl.newspec','a').write("'%s':'%s',\t# Not in UniProt\n" % (getspecies,code))
            return self.dict['EnsSpec'][getspecies]
        elif getspecies in enspec:
            self.dict['EnsSpec'][getspecies] = enspec[getspecies]
            return enspec[getspecies]
        else:
            for species in enspec.keys():
                if enspec[species] == getspecies:
                    if self.info['Name'] == getspecies: self.info['Name'] = species
                    getspecies = species
                    self.dict['EnsSpec'][getspecies] = enspec[getspecies]
                    return enspec[getspecies]
#########################################################################################################################
    def getSpecies(self,species):   ### Returns the EnsEMBL species for a given species or code (from enspec dict)
        '''Returns the EnsEMBL species for a given species or code (from enspec dict).'''
        species = species.replace(' ','_')
        if species in self.dict['EnsSpec']: return species
        elif self.db('release') and self.db('release').data(species): return species
        elif self.db('release') and species in self.db('release').index('spcode'):
            return self.db('release').indexEntries('spcode',species)[0]['species']
        elif self.obj['Taxonomy']:
            species = self.obj['Taxonomy'].getSpecies(species)
            species = rje.join(species.split()[:2],'_')
            return species
        elif species in enspec: return species          # Species is OK as it is.
        else:
            for spec in enspec.keys():
                if enspec[spec] == species: return spec     # Species code converted
        self.errorLog('Species "%s" not found!' % species,printerror=False) 
        return ''                                           # No species found
#########################################################################################################################
    def ensType(self,spec=None): ### Returns Ensembl subtype for species
        '''Returns Ensembl subtype for species.'''
        if not spec: spec = self.getStr('Name')
        try: return self.db('release').data(spec)['section']
        except:
            if self.dev(): self.warnLog("Still using self.dict['SiteSpec'].",dev=True)
            for esite in ['main'] + rje.sortKeys(self.dict['SiteSpec']):
                if spec in self.dict['SiteSpec'][esite]: return esite
#########################################################################################################################
    def clear(self):    ### Clears all the dictionaries except Stats
        '''Clears all the dictionaries except Stats.'''
        for key in self.dict.keys():
            if key != 'Stats': self.dict['Stats'] = {}
#########################################################################################################################
    def summaryStats(self,setup=False): ### Outputs summary stats
        '''Outputs summary stats.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['Species','Code','Common','Release','KnownProt','KnownLoci','NovelProt','NovelLoci','PutativeProt','PutativeLoci',
                       'SProt','TrEMBL','Known','Novel','Redundant','Replaced']
            sumfile = self.getStr('EnsPath') + 'ens_loci.tdt'
            sdb = self.db('ensloci')
            if not sdb: sdb = self.db().addTable(sumfile,mainkeys=['Species'],name='ensloci',expect=False)
            if not sdb: sdb = self.db().addEmptyTable('ensloci',headers,['Species'])
            else: self.log.no_suppression.append('entry_overwrite')
            #if setup and not self.info['Resume'] and not self.opt['Append']:
            #    rje.backup(self,sumfile)
            #    rje.delimitedFileOutput(self,sumfile,headers,'\t',datadict={})
            ### ~ [2] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not setup:
                datadict = self.dict['Stats']
                datadict['Species'] = self.info['Name']
                datadict['Code'] = self.speciesCode()
                try: datadict['Common'] = self.db('release').data(datadict['Species'])['common']
                except:
                    if datadict['Code'] in enscom: datadict['Common'] = enscom[datadict['Code']]
                #rje.delimitedFileOutput(self,sumfile,headers,'\t',datadict)
                sdb.addEntry(datadict)
                sdb.saveToFile(sumfile,backup=False)
        except: self.errorLog('Error in EnsEMBL.summaryStats')
#########################################################################################################################
    ### <3> ### EnsEMBL download subroutines                                                                            #
#########################################################################################################################
    def ensemblRelease(self,save=True):   ### Extracts Ensembl release data, or reads from existing table.
        '''Extracts Ensembl release data, or reads from exisiting table.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            rje.mkDir(self,'%srelease/' % self.getStr('EnsPath'))
            rfields = ['species','section','release']     # taxid, spcode and common name added in self.addTaxonomy()
            rfile = '%srelease/ensembl.release.tdt' % self.getStr('EnsPath')
            #i# Previous current list is now self.db('release').dataKeys() = List of current Ensembl species
            #i# Previous sitespec dictionary is now self.db('release').index('section') = Dictionary of species lists for each section
            ### ~ [1] Check for existing file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setBool({'RDBMod':False})
            if self.force() or not rje.exists(rfile): rdb = None
            else: rdb = db.addTable(rfile,mainkeys=['species'],name='release')
            if rdb: self.dict['SiteSpec'] = rdb.index('section'); return rdb
            self.setBool({'RDBMod':True})
            ### ~ [2] Extract data from FTP site ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            release = {}    # Dictionary of section:release
            rdb = db.addEmptyTable('release',rfields,['species'])
            ## ~ [2a] Extract species list from current_fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for line in urllib.urlopen('ftp://ftp.ensembl.org/pub/current_fasta/').readlines():
                line = '%s ' % rje.chomp(line)
                if rje.matchExp(' ([A-Za-z]+_[A-Za-z]+) ',line):
                    newspec = rje.matchExp(' ([A-Za-z]+_[A-Za-z]+) ',line)[0]
                    newspec = newspec[:1].upper() + newspec[1:]
                    rdb.addEntry({'species':newspec,'section':'main','release':'current'})
            ## ~ [2b] Get species from subsections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for esite in ['metazoa','fungi','plants','protists']:
                release[esite] = 10
                for line in urllib.urlopen('ftp://ftp.ensemblgenomes.org/pub/%s/' % esite).readlines():
                    if rje.matchExp('release-(\d+)',line):
                        r = rje.atoi(rje.matchExp('release-(\d+)',line)[0])
                        release[esite] = max(r,release[esite])
                self.printLog('#REL','Ensembl %s Release %d' % (esite,release[esite]))
                for line in urllib.urlopen('ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/fasta/' % (esite,release[esite])).readlines():
                    line = '%s ' % rje.chomp(line)
                    if rje.matchExp(' ([A-Za-z]+_[A-Za-z]+) ',line):
                        newspec = rje.matchExp(' ([A-Za-z]+_[A-Za-z]+) ',line)[0]
                        newspec = newspec[:1].upper() + newspec[1:]
                        if newspec not in rdb.dataKeys(): rdb.addEntry({'species':newspec,'section':esite,'release':release[esite]})
                        else: rdb.data(newspec)['section'] = esite; rdb.data(newspec)['release'] = release[esite]
            rdb.dropEntriesDirect('species',['Mart','Multi_species','Ancestral_alleles'])
            self.printLog('#ENSPEC','%s species identified from latest Ensembl release.' % rje.iStr(rdb.entryNum()))
            self.dict['SiteSpec'] = rdb.index('section')
            if save: rdb.saveToFile(rfile)  #!# Should only save if modified
            return rdb
        except: self.errorLog('Fatal error with ensemblRelease().',True); raise
#########################################################################################################################
    def ensemblReleaseTaxonomy(self,save=True):  ### Extracts Ensembl release data, or reads from existing table.
        '''Extracts Ensembl release data, or reads from existing table.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tobj = self.obj['Taxonomy']
            rdb = self.ensemblRelease(save=False)
            rdbmod = self.getBool('RDBMod')     # Whether RDB modified and needs saving
            rfile = '%srelease/ensembl.release.tdt' % self.getStr('EnsPath')
            if not rdb: raise ValueError('Ensembl release extraction failed')
            for field in ['taxid','spcode','name','common']:
                if field in rdb.fields() and self.force(): rdb.dropField(field)
                if field not in rdb.fields(): rdb.addField(field); rdbmod = True
            ### ~ [1] Map species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in rdb.entries():
                #!# Add some counters? #!#
                ## ~ [1a] TaxID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if entry['taxid']: taxid = entry['taxid']
                elif entry['species'] in self.dict['EnsTaxID']:
                    taxid = self.dict['EnsTaxID'][entry['species']]
                    if taxid: entry['taxid'] = taxid; rdbmod = True
                else:
                    try:
                        taxid = tobj.mapToTaxID(entry['species'],nodeonly=True)
                        if len(taxid) > 1:
                            for tax in taxid[0:]:
                                if tax in tobj.dict['TaxDict']: self.debug(tobj.dict['TaxDict'][tax])
                                tspec = tobj.taxDict(tax,skipuni=True)['name']
                                if tspec and tspec.startswith(rje.replace(entry['species'],'_',' ')):
                                    taxid.remove(tax)
                                    taxid.insert(0,tax)
                                    self.printLog('#TAXID','TaxID %s = "%s" in NCBI.' % (tax,tspec))
                                    break
                            self.printLog('#TAXID','Using %s (not %s) for "%s" TaxID.' % (taxid[0],rje.join(taxid[1:],'|'),entry['species']))
                        if taxid: taxid = taxid[0]
                        else: taxid = ''
                    except:
                        if self.dev() or self.test(): self.errorLog('?!')
                        taxid = ''
                    if taxid: entry['taxid'] = taxid; rdbmod = True
                ## ~ [1b] SpecDat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                specdat = True
                for field in ['spcode','name','common']: specdat = specdat and entry[field]
                if not specdat:
                    if tobj: taxdict = tobj.taxDict(taxid,store=not self.getBool('MemSaver'))
                    else: taxdict = {'spcode':self.speciesCode(entry['species']),
                                     'name':entry['species'].replace('_',' '),
                                     'common':'%s. %s' % (entry['species'].upper()[:1],entry['species'].split()[1])}
                    for field in ['spcode','name','common']:
                        if not entry[field]:
                            try: entry[field] = taxdict[field]; rdbmod = rdbmod or taxdict[field]
                            except: entry[field] = ''
                if entry['spcode'] == entry['taxid']: entry['spcode'] = ''
                if not entry['spcode']:
                    if entry['taxid']: entry['spcode'] = tobj.getSpCode(entry['taxid'])
                    else: entry['spcode'] = self.speciesCode(entry['species'])
                if not entry['taxid']:
                    self.warnLog('Failed to extract TaxID for %s: check for hyphens etc. in real name.' % entry['species'])
            if save and rdbmod: rdb.saveToFile(rfile)  #!# Modify to only save if modified.
            ### ~ [2] Reduce species by section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sections = []
            for spec in self.list['Sections'][0:]:
                if spec in rdb.index('section'):
                    sections.append(spec)
            if sections:    # Reduce species to these sections
                self.printLog('#SPEC','Reduce processing to %s species.' % rje.join(sections,'/'))
                rdb.dropEntriesDirect('section',sections,inverse=True)
                self.dict['SiteSpec'] = rdb.index('section')
            return rdb
        except: self.errorLog('Fatal error with ensemblReleaseTaxonomy().',True); raise
#########################################################################################################################
    def makeEnsemblUniprot(self):  ### Downloads uniprot for identified Ensembl species.
        '''Downloads uniprot for identified Ensembl species.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            makeindex = not self.list['SpecList']   # Whether index is safe to make
            rdb = self.ensemblReleaseTaxonomy()     # Downloads release data and extracts species information
            if not rdb: raise ValueError('Ensembl release/species extraction failed')
            if self.list['SpecList']:
                if not self.obj['Taxonomy'].getBool('Setup'): self.obj['Taxonomy'].setup()
                taxidsubset = self.obj['Taxonomy'].mapToTaxID(self.list['SpecList'],nodeonly=True)
            ### ~ [1] Download Uniprot entries for species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ifile = '%suniprot/uniprot.index' % self.getStr('EnsPath')
            umod = not rje.exists(ifile)    # Whether index needs to be made
            for section in rdb.index('section'):
                ufile = '%suniprot/%s.dat' % (self.getStr('EnsPath'),section)
                if rje.exists(ufile): usize = os.path.getsize(ufile)
                else: usize = 0
                pfile = '%s.pickup' % ufile
                sfile = '%s.partial' % ufile
                if self.getBool('MakeUniprot') and (self.force() or not rje.exists(ufile) or rje.exists(pfile) or rje.exists(sfile)):
                    if rje.exists(sfile) and not rje.exists(pfile): rje.fileTransfer(sfile,pfile,False,False)
                    rje.mkDir(self,'%suniprot/' % self.getStr('EnsPath'))
                    if self.obj['Taxonomy']:
                        uniprot = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=url'])
                        uniprot.setStr({'DatOut':ufile,'UniFormat':'txt'})
                        #uniprot.list['Taxonomy'] = rdb.dataList(rdb.entries(),'taxid',empties=False)
                        uniprot.list['Taxonomy'] = rdb.indexDataList('section',section,'taxid')
                        #self.debug(self.list['SpecList'])
                        if self.list['SpecList']:
                            uniprot.list['Taxonomy'] = taxidsubset
                            if rje.exists(sfile) and not rje.exists(pfile): rje.fileTransfer(fromfile=sfile,tofile=pfile,deletefrom=False,append=False)
                        #self.debug(uniprot.list['Taxonomy'])
                        uniprot._extractProteomesFromURL(uniprot.list['Taxonomy'],cleardata=True,logft=False,reformat=False,taxonomy=True,pickup=not self.force())
                        if not rje.exists(pfile):
                            if self.list['SpecList']:
                                if rje.exists(sfile): partial = rje.sortUnique(rje.listFromCommand(sfile)+taxidsubset)
                                else: partial = taxidsubset
                                open(sfile,'w').write(rje.join(partial+[''],'\n'))
                            elif rje.exists(sfile): os.unlink(sfile)
                    else: self.warnLog('Cannot make UniProt file without rje_taxonomy setup')
                if rje.exists(ufile):
                    umod = umod or os.path.getsize(ufile) != usize
                else: self.errorLog('%s not created!' % ufile,printerror=False); return None
            ### ~ [2] Process downloaded Uniprot files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setStr({'UniPath':rje.makePath('%suniprot/' % self.getStr('EnsPath'))})
            if makeindex:
                if rje.exists(ifile): os.unlink(ifile); self.printLog('#INDEX','Remaking Uniprot index file.')
                else: self.printLog('#INDEX','Making Uniprot index file.')
                rje_uniprot.processUniProt(self,makeindex=True,makespec=False,makefas=False)
            return rdb
        except: self.errorLog('Fatal error with makeEnsemblUniprot().',True); return None
#########################################################################################################################
    def downloadEnsEMBL(self,speclist=[]):   ### Downloads the relevant files for the desired species (and runs EnsLoci)
        '''
        Downloads the relevant files for the desired species. If no species is given, all current species will be
        downloaded. Each file is tried 5 times. If all 5 attempts fail, errors are reported. If ensloci=T, this method
        will run EnsLoci following a successful download, to help space out the downloads. Note that the UniProt DAT
        files should therefore be downloaded and indexed first.
        >> speclist:list of species names, matching keys of enspec
        << returns list of successful species processed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = self.ensemblReleaseTaxonomy()     # Downloads release data and extracts species information
            if not rdb: raise ValueError('Ensembl release/species extraction failed')
            goodspec = []   # List of species downloaded and processed successfully
            sitespec = self.dict['SiteSpec']
            release = {}
            for esite in sitespec:
                if esite != 'main': release[esite] = int(rdb.indexDataList('section',esite,'release')[0])
                else: release[esite] = 0
            release['main'] = max(release.values())
            ### ~ [1] Get SpecList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            current = rdb.dataKeys()    # List of current Ensembl species
            if not speclist: speclist = current[0:]  # All Ensembl species recognised
            else:
                for spec in speclist[0:]:
                    espec = self.getSpecies(spec)
                    if espec not in current:
                        self.warnLog('Species "%s" missing from current species (case-sensitive): %s removed' % (espec,spec))
                        speclist.remove(spec)
                    elif spec != espec: speclist.remove(spec); speclist.append(espec)
            speclist = rje.sortUnique(speclist)
            #i# Removed determination of old and new species using enspec.
            ### ~ [2] Create directories and Download Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Resume'): self.info['Resume'] = ''
            mydir = os.path.abspath(os.curdir)
            problemfiles = []; downloaded = False
            rje.mkDir(self,self.getStr('EnsPath'))
            sx = 0.0; stot = len(speclist)
            for spec in speclist:
                esite = rdb.data(spec)['section']
                self.bugPrint(spec)
                ## ~ [2a] Get species code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                code = self.speciesCode(spec)
                if self.info['Resume']:
                    if self.info['Resume'] in [spec,code]: self.info['Resume'] = ''
                    else: stot -= 1; continue
                goodspec.append(spec); sx += 1
                ## ~ [2b] Create directory & Sleep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                spdir = rje.makePath('%s%s/%s/' % (self.info['EnsPath'],esite,code))
                rje.mkDir(self,spdir)
                #self.deBug(spdir)
                if self.stat['SpecSleep'] and speclist.index(spec) > 0 and downloaded:
                    try:
                        for sec in range(self.stat['SpecSleep'],0,-1):
                            self.progLog('\r#WAIT','SpecSleep for %d seconds!' % sec)
                            time.sleep(1)
                        self.progLog('\r#WAIT','SpecSleep for 0 seconds!   ')
                    except: self.printLog('\r#WAIT','SpecSleep cancelled!       ',log=False)
                downloaded = False
                self.printLog('\r#~~~#','#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#',timeout=False,log=False)
                self.printLog('\r#NEXT','Next species (%d of %d): %s [%s] (Ensembl %s)' % (sx,stot,spec,code,self.ensType(spec)))
                ## ~ [2c] Download protein fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fasdir = 'Cannot establish FTP download site for %s pep!' % spec
                try:
                    self.deBug('%s*.pep.all.fa' % spdir)
                    self.deBug(glob.glob('%s*.pep.all.fa' % spdir))
                    if not self.force() and self.getBool('SpeedSkip') and glob.glob('%s*.pep.all.fa' % spdir):
                        self.printLog('#DBASE','Skipping downloaded %s peptide sequences (speedskip=T).' % spec)
                    else:
                        #!# Should make this a method that can be used for all data types: will be easier to fix errors etc.
                        if esite == 'main':
                            fasdir = 'ftp://ftp.ensembl.org/pub/current_fasta/%s/pep/' % spec.lower()
                        else:
                            fasdir = 'ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/fasta/%s/pep/' % (esite,release[esite],spec.lower())
                            #try: urllib.urlopen(fasdir).readlines()
                            #except: fasdir = 'ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/fasta/%s/' % (esite,release[esite],spec.lower())
                        fasgzlist = []
                        for attempt in range(5):
                            try: fasgzlist = urllib.urlopen(fasdir).readlines()
                            except:
                                self.errorLog(fasdir)
                                try:
                                    for sec in range(300 * (attempt+1),self.stat['SpecSleep']-1,-10):
                                        self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                        time.sleep(10)
                                    continue
                                except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                except: raise
                        if not fasgzlist: raise ValueError('Unable to read %s' % fasdir)
                        for fasgz in fasgzlist:
                            for ftype in ['all','known','known-ccds','novel','abinitio']:
                                if rje.matchExp('(%s\..+%s.fa.gz)' % (spec,ftype),fasgz):
                                    balls = False
                                    for attempt in range(5):
                                        balls = False
                                        fa = rje.matchExp('(%s\.\S+%s.fa.gz)' % (spec,ftype),fasgz)[0]
                                        wget = '%s%s' % (fasdir,fa)
                                        checkget = rje.makePath('%s/%s' % (spdir,fa[:-3]),wholepath=True)
                                        #self.deBug('%s: %s' % (checkget,os.path.exists(checkget)))
                                        if not self.opt['Force'] and os.path.exists(checkget):  # rje.isYounger(fa[:-3],ensloci) == fa[:-3]:   ### Check date of files
                                            self.printLog('#DBASE','Skipping downloaded %s' % (fa))
                                            break
                                        self.log.printLog('\r#DBASE','Downloading %s' % fa); downloaded = True
                                        os.chdir(spdir)
                                        if self.getBool('OSX'): os.system('curl -O %s' % wget)
                                        else: os.system('wget %s -t 5' % wget)
                                        if not os.path.exists(fa):
                                            problemfiles.append(wget)
                                            balls = True
                                        elif not self.opt['Win32']:
                                            if os.system('gunzip %s -f' % fa):
                                                problemfiles.append(wget)
                                                balls = True
                                        os.chdir(mydir)
                                        if balls:
                                            self.errorLog('Something wrong with %s - check download!' % fa,printerror=False)
                                            try:
                                                for sec in range(120 * (attempt+1),self.stat['SpecSleep']-1,-1):
                                                    self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                                    time.sleep(1)
                                                continue
                                            except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                            except: raise
                                        else:
                                            self.log.printLog('\r#GZ','%s downloaded and unzipped.' % (fa))
                                            break
                                    if balls and spec in goodspec: goodspec.remove(spec)
                except:
                    self.errorLog(fasdir)
                    if spec in goodspec: goodspec.remove(spec)
                ## ~ [2c] Download cDNA fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #fasdir = 'ftp://ftp.ensembl.org/pub/current_%s/data/fasta/pep/' % spec.lower()
                fasdir = 'Cannot establish FTP download site for %s cdna!' % spec
                if esite == 'main':
                    fasdir = 'ftp://ftp.ensembl.org/pub/current_fasta/%s/cdna/' % spec.lower()
                else:
                    fasdir = 'ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/fasta/%s/cdna/' % (esite,release[esite],spec.lower())
                try:
                    if not self.force() and self.getBool('SpeedSkip') and glob.glob('%s*.cdna.all.fa' % spdir):
                        self.printLog('#DBASE','Skipping downloaded %s cDNA sequences (speedskip=T).' % spec)
                    else:
                        fasgzlist = []
                        for attempt in range(5):
                            try: fasgzlist = urllib.urlopen(fasdir).readlines()
                            except:
                                self.errorLog(fasdir)
                                try:
                                    for sec in range(300 * (attempt+1),self.stat['SpecSleep']-1,-10):
                                        self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                        time.sleep(10)
                                    continue
                                except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                except: raise
                        if not fasgzlist: raise ValueError('Unable to read %s' % fasdir)
                        for fasgz in fasgzlist:
                            for ftype in ['all','known','known-ccds','novel','abinitio']:
                                if rje.matchExp('(%s\..+%s.fa.gz)' % (spec,ftype),fasgz):
                                    balls = False
                                    for attempt in range(5):
                                        balls = False
                                        fa = rje.matchExp('(%s\.\S+%s.fa.gz)' % (spec,ftype),fasgz)[0]
                                        wget = '%s%s' % (fasdir,fa)
                                        checkget = rje.makePath('%s/%s' % (spdir,fa[:-3]),wholepath=True)
                                        if not self.opt['Force'] and os.path.exists(checkget):  # rje.isYounger(fa[:-3],ensloci) == fa[:-3]:   ### Check date of files
                                            self.printLog('#DBASE','Skipping downloaded %s' % (fa))
                                            break
                                        self.log.printLog('\r#DBASE','Downloading %s' % fa); downloaded = True
                                        os.chdir(spdir)
                                        if self.getBool('OSX'): os.system('curl -O %s' % wget)
                                        else: os.system('wget %s -t 5' % wget)
                                        if not os.path.exists(fa):
                                            problemfiles.append(wget)
                                            balls = True
                                        elif not self.opt['Win32']:
                                            if os.system('gunzip %s -f' % fa):
                                                problemfiles.append(wget)
                                                balls = True
                                        os.chdir(mydir)
                                        if balls:
                                            self.errorLog('Something wrong with %s - check download!' % fa,printerror=False)
                                            try:
                                                for sec in range(120 * (attempt+1),self.stat['SpecSleep']-1,-1):
                                                    self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                                    time.sleep(1)
                                                continue
                                            except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                            except: raise
                                        else:
                                            self.log.printLog('\r#GZ','%s downloaded and unzipped.' % (fa))
                                            break
                                    if balls and spec in goodspec: goodspec.remove(spec)
                except:
                    self.errorLog(fasdir)
                ## ~ [2d] Download chromosome fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fasdir = 'Chromosome download error for %s' % spec
                if not self.force() and self.getBool('SpeedSkip') and glob.glob('%s*.dna.toplevel.fa' % spdir):
                    self.printLog('#DBASE','Skipping downloaded %s TopLevel DNA sequences (speedskip=T).' % spec)
                elif code in self.list['ChromSpec'] or 'ALL' in self.list['ChromSpec'] or '*' in self.list['ChromSpec'] or esite.upper() in self.list['ChromSpec']:
                    self.printLog('#CHROM','Downloading chromosomes for %s' % spec)
                    try:
                        if esite == 'main': fasdir = 'ftp://ftp.ensembl.org/pub/current_fasta/%s/dna/' % spec.lower()
                        else: fasdir = 'ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/fasta/%s/dna/' % (esite,release[esite],spec.lower())
                        fasgzlist = []
                        for attempt in range(5):
                            try: fasgzlist = urllib.urlopen(fasdir).readlines()
                            except:
                                self.errorLog(fasdir)
                                try:
                                    for sec in range(300 * (attempt+1),self.stat['SpecSleep']-1,-10):
                                        self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                        time.sleep(10)
                                    continue
                                except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                except: raise
                        if not fasgzlist: raise ValueError('Unable to read %s' % fasdir)
                        for fasgz in fasgzlist:
                            #fa = rje.matchExp('(%s\..+\.dna\.chromosome\.\S{1,3}\.fa\.gz)' % (spec),fasgz)
                            fa = rje.matchExp('(%s\..+\.dna\.toplevel\.fa\.gz)' % (spec),fasgz)
                            if fa:
                                fa = fa[0]
                                balls = False
                                for attempt in range(5):
                                    balls = False
                                    wget = '%s%s' % (fasdir,fa)
                                    checkget = rje.makePath('%s/%s' % (spdir,fa[:-3]),wholepath=True)
                                    if not self.opt['Force'] and os.path.exists(checkget):  #rje.isYounger(fa[:-3],ensloci) == fa[:-3]:   ### Check date of files
                                        self.printLog('#DBASE','Skipping downloaded %s' % (fa))
                                        break
                                    self.printLog('\r#DBASE','Downloading %s' % fa); downloaded = True
                                    os.chdir(spdir)
                                    if self.getBool('OSX'): os.system('curl -O %s' % wget)
                                    else: os.system('wget %s -t 5' % wget)
                                    if not os.path.exists(fa):
                                        problemfiles.append(wget)
                                        balls = True
                                    elif not self.opt['Win32']:
                                        if os.system('gunzip %s -f' % fa):
                                            problemfiles.append(wget)
                                            balls = True
                                    os.chdir(mydir)
                                    if balls:
                                        self.errorLog('Something wrong with %s!' % fa,printerror=False)
                                        if os.path.exists('%s%s' % (spdir,fa)):
                                            os.unlink('%s%s' % (spdir,fa))
                                            self.printLog('#DEL','Deleted incomplete/corrupted %s.' % fa)
                                        try:
                                            for sec in range(120 * (attempt+1),self.stat['SpecSleep']-1,-1):
                                                self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                                time.sleep(1)
                                            continue
                                        except KeyboardInterrupt: self.printLog('\r#WAIT','Sleeping cancelled!       ',log=False)
                                        except: raise
                                    else:
                                        self.log.printLog('\r#GZ','%s downloaded and unzipped.' % (fa))
                                        break
                    except: self.errorLog(fasdir)
                ## ~ [2e] Download MySQL Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #sqldir = 'ftp://ftp.ensembl.org/pub/current_%s/data/mysql/' % spec.lower()
                sqldir = 'Trouble finding %s SQL files' % spec
                if esite == 'main': sqldir = 'ftp://ftp.ensembl.org/pub/current_mysql/'
                else: sqldir = 'ftp://ftp.ensemblgenomes.org/pub/%s/release-%d/mysql/' % (esite,release[esite])
                try:
                    if not self.force() and self.getBool('SpeedSkip') and glob.glob('%sgene.txt' % spdir):
                        self.printLog('#DBASE','Skipping downloaded %s gene table (speedskip=T).' % spec)
                    else:
                        coredir = ''
                        for line in urllib.urlopen(sqldir).readlines():
                            if rje.matchExp('(%s_core_\S+)' % spec.lower(),line):
                                coredir = sqldir + rje.matchExp('(%s_core_\S+)' % spec.lower(),line)[0]
                                break
                        if coredir:
                            for table in ['gene']:  #,'xref','gene_stable_id']:
                                balls = False
                                for attempt in range(5):
                                    balls = False
                                    wget = coredir + '/%s.txt.gz' % table   #!# table removed from file names! (Dec '07)
                                    checkget = rje.makePath('%s/%s.txt' % (spdir,table),wholepath=True)
                                    if not self.opt['Force'] and os.path.exists(checkget):  # and rje.isYounger('%s.txt' % table,ensloci) != ensloci:   ### Check date of files
                                        self.printLog('#DBASE','Skipping downloaded %s %s' % (spec,table))
                                        break
                                    self.log.printLog('#DBASE','Downloading %s %s' % (spec,table))
                                    os.chdir(spdir)
                                    if self.getBool('OSX'): os.system('curl -O %s' % wget)
                                    else: os.system('wget %s -t 5' % wget)
                                    if os.path.exists('%s.txt.table.gz' % table): os.rename('%s.txt.table.gz' % table,'%s.txt.gz' % table)
                                    if not os.path.exists('%s.txt.gz' % table): balls = True
                                    elif not self.opt['Win32']:
                                        if os.system('gunzip %s.txt.gz -f' % table): balls = True
                                    os.chdir(mydir)
                                    if balls:
                                        if table == 'gene_stable_id': self.printLog('#GENE','Assuming new gene table format ("gene" file enough).'); break
                                        downloaded = True
                                        problemfiles.append(wget)
                                        self.errorLog('Something wrong with %s - check download!' % table,printerror=False)
                                        try:
                                            for sec in range(120 * (attempt+1),self.stat['SpecSleep']-1,-1):
                                                self.progLog('\r#WAIT','Attempt %d failed; will try again in %d seconds!' % (attempt+1,sec))
                                                time.sleep(1)
                                        except KeyboardInterrupt: self.printLog('\r#GO!!','Cancelled wait with %d seconds remaining!' % sec,log=False)
                                        except: raise
                                        continue
                                    else:
                                        self.log.printLog('\r#GZ','%s downloaded and unzipped.' % (table))
                                        downloaded = True
                                        break
                                if balls and spec in goodspec: goodspec.remove(spec)
                        else: self.log.errorLog('No core directory found in %s' % sqldir,printerror=False)
                except:
                    self.errorLog(sqldir)
                    if spec in goodspec: goodspec.remove(spec)
                ## ~ [2f] EnsLoci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.devPrint('GoodSpec: %s' % goodspec)
                if spec in goodspec and self.opt['EnsLoci']: self.ensLoci(spec)
            ### ~ [3] Tidy up and report problems ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Windows Unzip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['Win32']:
                self.printLog('#WIN','Cannot unzip in Windows, sorry.')
                return []
            ## ~ [3b] Report Problem Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if problemfiles:
                problemfiles = rje.sortUnique(problemfiles)
                for pfile in problemfiles: self.errorLog('Failed: %s' % pfile,printerror=False)
                self.printLog('#WGET','%d file(s) experienced problems during downloading.' % len(problemfiles))
                sugtxt = '\nSuggested action:\n1. Halt program with CTRL+Z\n2. Manually download and unzip files.\n3. Resume with "fg" and hit "Y"\n\n'
                if self.i() < 0 or not rje.yesNo('%sProceed regardless?' % sugtxt,default="N"): return goodspec
            else: self.printLog('#WGET','All files apparently downloaded successfully.')
            return speclist
        except: self.errorLog('Problem with rje_ensembl.downloadEnsEMBL()',quitchoice=True); return []
#########################################################################################################################
    ### <4> ### Loading of Data into data structure                                                                     #
#########################################################################################################################
    def uniprotMapDB(self): ### Generates a mapping table by extracting the data direct from Uniprot
        '''Generates a mapping table by extracting the data direct from Uniprot.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            spcode = self.speciesCode()
            parsedb = ['ensembl']   #,'gene','id','HGNC','MGI','FlyBase','WormBase','VectorBase','SGD','ZFIN']
            ucmd = ['specdat=%s' % spcode,'dbparse=%s' % rje.join(parsedb,','),'memsaver=F','uparse=OS,OC,OX,DR,SY']
            ucmd += ['reviewed=T']
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list+ucmd)
            ### ~ [1] ~ Read from Uniprot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniprot.extractSpecies()
            uniprot.setupDB(xref=True,ft=False)
            uniprot.readUniProt()
            xdb = uniprot.db('xref')    # ['accnum','db','xref'],['accnum','db']
            ## ~ [1a] ~ Convert Uniprot xref table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            udb = self.db().addEmptyTable('unixref',['Uniprot','UniprotID','ENSG','ENSP','Symbol','Gene'],['Uniprot'])
            x = 0.0; xtot = len(xdb.index('accnum'))
            for accnum in xdb.index('accnum'):
                self.progLog('\r#XREF','Converting Uniprot XRef: %.2f%%' % (x/xtot)); x += 100.0
                uentry = {'Uniprot':accnum,'UniprotID':'','ENSG':'','ENSP':'','Symbol':[],'Gene':''}
                for xentry in xdb.indexEntries('accnum',accnum):
                    if xentry['db'] in ['ENSG','ENSP','Gene']: uentry[xentry['db']] = xentry['xref']
                    elif xentry['db'] == 'ID': uentry['UniprotID'] = xentry['xref']
                    else: uentry['Symbol'].append(xentry['xref'])
                uentry['Symbol'] = rje.join(uentry['Symbol'],'|')
                udb.addEntry(uentry)
            self.printLog('\r#XREF','Converted XRef for %s of %s Uniprot entries.' % (rje.iStr(udb.entryNum()),rje.iStr(xtot)))
            ## ~ [1b] ~ Generate and store dictionary of Uniprot acc to UniprotEntry object ~~~~~~~ ##
            self.dict['UniprotObj'] = uniprot.accDictFromEntries()
            return True

            #!# >>> DELETE BELOW HERE ONCE HAPPY #!#
            mdb = self.db().addEmptyTable('Map',['ENSG','UniProt','Symbol'],['ENSG'])
            spcode = self.speciesCode()
            ucmd = ['specdat=%s' % spcode,'dbparse=ensembl,HGNC,FlyBase,MGI,ZFIN,gene,id','memsaver=F']
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list+ucmd)
            ### ~ [1] ~ Read from Uniprot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniprot.extractSpecies()
            uniprot.setupDB(xref=True,ft=False)
            uniprot.readUniProt()
            acc2uni = uniprot.accDictFromEntries()
            ### ~ [2] ~ Generate MapDB from xref table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = uniprot.db('xref')    # ['accnum','db','xref'],['accnum','db']
            self.printLog('#XREF','%s AccNum-DB links parsed for %s.' % (rje.iStr(xdb.entryNum()),spcode))
            for accnum in xdb.indexDataList('db','ENSG','accnum'):
                entry = {'UniProt':accnum}
                ensg = rje.split(xdb.data('%s\tENSG' % accnum)['xref'],'|')
                for db in ['Symbol','Gene','FlyBase']:
                    ## ~ [2a] ~ Extract IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    symb = xdb.data('%s\t%s' % (accnum,db))
                    if not symb: continue
                    symb = rje.split(symb['xref'],'|')
                    if len(symb) > 1:
                        if len(ensg) > 1:
                            self.warnLog('%s Uniprot %s maps to multiple genes/symbols. Cannot match to ENSG via Uniprot. (%s vs %s)' % (spcode,accnum,rje.join(ensg,'|'),rje.join(symb,'|')))
                        else: self.warnLog('%s Uniprot %s maps to %d Gene symbols. Will use %s. (Not %s)' % (spcode,accnum,len(symb),symb[0],rje.join(symb[1:],'|')))
                    entry['Symbol'] = symb[0]
                    ## ~ [2b] ~ Check for ENSG:many mapping and pick "best" Uniprot ~~~~~~~~~~~~~~~ ##
                    for egene in ensg:
                        if egene in mdb.data():
                            olduni = acc2uni[mdb.data(egene)['UniProt']]
                            newuni = acc2uni[accnum]
                            if olduni.reviewed() and not newuni.reviewed(): continue
                            elif olduni.reviewed() == newuni.reviewed() and olduni.length() >= newuni.length(): continue
                            self.printLog('#RED','%s: %s (%s; %saa) replacing %s (%s; %saa).' % (egene,newuni.id(),newuni.reviewed(),newuni.length(),olduni.id(),olduni.reviewed(),olduni.length()))
                        mdb.addEntry(rje.combineDict({'ENSG':egene},entry))
            self.printLog('#MAP',"Mapped Uniprot data for %s ENSG." % rje.iStr(mdb.entryNum()))
            return True
        except: self.errorLog('%s.uniprotMapDB(%s) problem.' % (self,self.speciesCode())); return False
#########################################################################################################################
    def mapDB(self):    ### Look for and load *.map.tdt file, if appropriate
        '''Look for and load *.map.tdt file, if appropriate.'''
        ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        db = self.db(); hdb = None
        mapfile = '%s%s.map.tdt' % (self.info['EnsPath'],self.speciesCode())
        ### ~ [1] ~ Load special HGNC Mapping file if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if rje.exists(self.getStr('HGNCMap')) and self.speciesCode() == 'HUMAN':
            hdb = db.addTable(self.getStr('HGNCMap'),['HGNC ID'],['Approved Symbol','UniProt ID (mapped data supplied by UniProt)'],name='hgnc')
            if hdb:
                hdb.renameField('HGNC ID','HGNC')
                hdb.renameField('Approved Symbol','Symbol')
                hdb.renameField('UniProt ID (mapped data supplied by UniProt)','Uniprot')
            else:
                hdb = db.addTable(self.getStr('HGNCMap'),['Gene'],name='hgnc')
                if hdb: hdb.renameField('Ensembl','ENSG')
        ### ~ [2] ~ Generic mapping file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        mdb = db.addTable(mapfile,['Ensembl Gene ID'],name='map',expect=False)
        if not mdb: return hdb
        mdb.renameField('Ensembl Gene ID','ENSG')
        mdb.renameField('UniProt/SwissProt Accession','Uniprot')
        for spfield in ['MGI Curated Symbol','MGI symbol','ZFIN symbol','FlyBaseName gene','UniProt Gene Name']:
            if spfield in mdb.fields(): mdb.renameField(spfield,'Symbol')
        return True
#########################################################################################################################
    def uniFromMap(self,ensg):   ### Returns mapped UniProt AccNum from gene, if mapping loaded
        '''Returns mapped UniProt AccNum from mapping, if found.'''
        try: return self.db('Map').data()[ensg]['UniProt']
        except: return ''
#########################################################################################################################
    def symbolFromMap(self,ensg):   ### Returns mapped UniProt AccNum from gene, if mapping loaded
        '''Returns mapped UniProt AccNum from mapping, if found.'''
        try:
            gene = self.db('Map').data()[ensg]['Symbol']
            if gene.find(':') > 0: return ''
            if self.speciesCode() == 'DROME' and gene[:2] == 'CG': return ''
            symbol = ''
            for x in gene:
                if x not in string.punctuation: symbol += x
            return symbol
        except: return ''
#########################################################################################################################
    def uniFromHGNC(self,acc):  ### Returns mapped UniProt AccNum from HGNCID, if found
        '''Returns mapped UniProt AccNum from HGNCID, if found.'''
        try:
            if self.speciesCode() != 'HUMAN': return ''
            try: return self.db('hgnc').indexDataList('HGNC',acc,'Uniprot')[0]
            except: return self.db('hgnc').indexDataList('HGNC','HGNC:%s' % acc,'Uniprot')[0]
        except: return ''
#########################################################################################################################
    def symbolFromHGNC(self,acc):  ### Returns mapped UniProt AccNum from HGNCID, if found
        '''Returns mapped Gene Symbol HGNCID, if found.'''
        try:
            try: return self.db('hgnc').indexDataList('HGNC',acc,'Symbol')[0]
            except: return self.db('hgnc').indexDataList('HGNC','HGNC:%s' % acc,'Symbol')[0]
        except: return ''
#########################################################################################################################
    def parseEnsEMBL(self,species='',makensloci=True):  ### Parse data from EnsEMBL files
        '''
        Parse data from EnsEMBL files.
        >> species:str = Species to parse data for
        >> makensloci:bool = Whether data is being parsed for EnsLoci treatment (need UniProt data) [True]
        '''
        try:### ~ [0] ~ Setup parsing. Make UniProt object and get file names etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not makensloci: self.warnLog('Defunct makensloci parameter being set',dev=True)
            self.setBool({'Parsed':False})
            if not self.obj['UniProt']: self.obj['UniProt'] = rje_uniprot.UniProt(self.log,self.cmd_list)
            self.setStr({'Name':species})
            ## ~ [0a] ~ Check for species path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            specpath = rje.makePath('%s%s/%s' % (self.getStr('EnsPath'),self.ensType(species),self.speciesCode()))
            if not os.path.exists(specpath):
                try: common = enscom[self.speciesCode()]
                except: common = '???'
                self.printLog('#ENS','Directory %s not found for %s (%s). Check Ensembl download/availability.' % (specpath,species,common))
                return False
            ## ~ [0b] ~ Identify fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Stats']['Release'] ='?'
            ens = {}
            for stype in ['pep','cdna']:
                ens[stype] = glob.glob('%s%s.*.%s.all.fa' % (specpath,species,stype))
                if not ens[stype]: raise IOError('Cannot find %s%s.*.%s.all.fa' % (specpath,species,stype))
                if len(ens[stype]) > 1: self.warnLog('Multiple %s%s.*.%s.all.fa files found!' % (specpath,species,stype))
                ens[stype] = ens[stype][0]
                self.printLog('#%s' % stype.upper(),ens[stype])
                try: self.dict['Stats']['Release'] = rje.matchExp('%s.(\S+).%s.all.fa' % (species,stype),ens[stype])[0]
                except: pass
            ## ~ [0c] ~ SQL tables with gene data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gene_table = specpath + 'gene.txt'      #x#.table'
            if not os.path.exists(gene_table) and os.path.exists('%s.table' % gene_table): gene_table += '.table'
            stable_id = specpath + 'gene_stable_id.txt' #x#.table'
            if not os.path.exists(stable_id) and os.path.exists('%s.table' % stable_id): stable_id += '.table'
            for path in [gene_table]:   # V2.7 update! ,stable_id]:
                if not os.path.exists(path):
                    self.errorLog('%s missing!' % path,printerror=False)
                    return None
            ## ~ [0d] ~ Check Dates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ensloci = self.ensLociFas()
            needtoparse = self.force() or not makensloci
            for path in [specpath,gene_table] + ens.values():
                if not path: continue
                #self.bugPrint('%s = %s' % (path,rje.isYounger(ensloci,path) != ensloci))
                needtoparse = needtoparse or (path and rje.isYounger(ensloci,path) != ensloci)  # Check date of files
            self.bugPrint('%s NeedtoParse = %s' % (species,needtoparse))
            if not needtoparse:
                if not os.path.exists(ensloci): self.errorLog('No good files found for %s' % species); return None
                self.setBool({'Parsed':True})
                self.printLog('#LOCI','%s already exists' % ensloci)
                return self                

            ### ~ [1] ~ Parse Gene Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Generate db tables with some of: 'ENSG','ENSP','Uniprot','Symbol','UniprotID' and 'Gene' fields.
            # 'map' = 'ENSG' to 'Uniprot','Symbol'
            # 'hgnc' = 'HGNC ID' or 'Gene' to 'ENSG','Uniprot' and 'Symbol' or 'Gene'
            self.mapDB()                # Uses mapping file to populate self.db('map') and/or self.db('hgnc')
            # 'unixref',['Uniprot','UniprotID','ENSG','ENSP','Symbol','Gene'],['Uniprot'])
            self.uniprotMapDB()     # Generates self.db('unixref') and self.dict['UniprotObj'] = {acc:UniProtObject}
            # 'gene' table = 'ENSG' to other details, including 'gene_id','biotype','status','protein' (Y/N),'sourcedb' and 'sourceacc'
            self.parseGeneDetails(gene_table)
            #i# No longer used self.parseStableIDs(stable_id)

            ### ~ [2] ~ Parse protein sequence details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','mode=index'])
            seqlist.loadSeq(ens['pep'],filetype=None,seqtype=None,nodup=True,clearseq=True,mode=None)
            ## ~ [2a] ~ Make peptide table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = self.db().addEmptyTable('peptide',['ENSP','ENSG','Desc','Seq','NonX'],['ENSP'])
            crapx = 0
            #while seqlist.nextSeq() != None:
            #    seq = seqlist.currSeq()
            for seq in seqlist.seqs():
                (name,sequence) = seqlist.getSeq(seq,'tuple')
                ensp = rje.split(name)[0]
                desc = rje.join(rje.split(name)[1:])
                try: ensg = rje.matchExp('gene:(\S+)\s',name)[0]
                except: crapx += 1; ensg = ''
                nonx = len(sequence) - rje.count(sequence.upper(),'X')
                pdb.addEntry({'ENSP':ensp,'ENSG':ensg,'Desc':desc,'Seq':sequence.upper(),'NonX':nonx})
            self.printLog('#ENSP','%s Ensembl peptide sequences added to peptide table.' % (pdb.iNum()))
            if crapx: self.warnLog('%s %s peptide sequences without recognisable gene annotation.' % (rje.iStr(crapx),species))
            pdb.index('ENSG')

            ### ~ [3] ~ Update parsing stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for btype in ['Known','Novel','Putative']:
                bpep = 0; bloci = 0
                for ensg in self.db('gene').indexDataList('biotype',btype.upper(),'ENSG'):
                    bloci += 1
                    if ensg in pdb.index('ENSG'): bpep += 1
                self.dict['Stats']['%sProt' % btype] = bpep
                self.dict['Stats']['%sLoci' % btype] = bloci

            ### ~ [4] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self
        except: self.errorLog('Error in EnsEMBL.parseEnsEMBL(%s)' % self.getStr('Name'),printerror=True,quitchoice=False)
        return None
#########################################################################################################################
    def parseGeneDetails(self,gene_table):  ### Loads EnsEMBL Gene Table into db table
        '''Parses out gene and protein IDs from EnsEMBL Gene Table.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ### ~ [1] ~ Add gene table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# http://asia.ensembl.org/info/docs/api/core/core_schema.html#gene
            gfields = ['gene_id', # Primary key, internal identifier.	primary key
                       'biotype', # Biotype, e.g. protein_coding.
                       'analysis_id',       #	SMALLINT		Foreign key references to the analysis table.	key: analysis_idx
                       'seq_region_id',     #	INT(10)		Foreign key references to the seq_region table.	key: seq_region_idx
                       'seq_region_start',  #	INT(10)		Sequence start position.	key: seq_region_idx
                       'seq_region_end',    #	INT(10)		Sequence end position.
                       'seq_region_strand', #	TINYINT(2)		Sequence region strand: 1 - forward; -1 - reverse.
                       'display_xref_id',   #	INT(10)		External reference for EnsEMBL web site. Foreign key references to the xref table.	key: xref_id_index
                       'source',    #	VARCHAR(20)		e.g ensembl, havana etc.
                       'status',    #	ENUM('KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN', 'ANNOTATED')		Status, e.g.'KNOWN', 'NOVEL', 'PUTATIVE', 'PREDICTED', 'KNOWN_BY_PROJECTION', 'UNKNOWN'.
                       'description',   #	TEXT		Gene description
                       'is_current',    #	BOOLEAN	1	1 - gene is current. Always set to 1 in ensembl dbs, but needed for otterlace dbs
                       'canonical_transcript_id',   #	INT(10)		Foreign key references to the transcript table.	key: canonical_transcript_id_idx
                       'stable_id',     #	VARCHAR(128)	NULL	Release-independent stable identifier.	key: stable_id_idx
                       'version	SMALLINT',  #	1	Stable identifier version number.	key: stable_id_idx
                       'created_date',  #	DATETIME	'0000-00-00	Date created.
                       'modified_date'] #	DATETIME	'0000-00-00	Date modified.
            gdb = self.db().addTable(gene_table,['stable_id'],delimit='\t',headers=gfields,name='gene')
            if not gdb: raise ValueError('Failed to load gene table')
            gdb.renameField('stable_id','ENSG')
            # Source:UniProtKB/TrEMBL;Acc:Q7YW39
            # Source:VB External Description;Acc:AGAP004725
            ## ~ [1a] ~ Identiy protein-coding genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb.addField('protein',evalue='N')
            for gtype in rje.sortKeys(gdb.index('biotype')):
                if gtype == 'protein_coding':   # or gtype[:3] in ['IG_','TR_']:    #!# Don't think IG_ and TR_ genes are wanted.
                    for entry in gdb.indexEntries('biotype',gtype): entry['protein'] = 'Y'
                    self.printLog('#TYPE','%s %s genes (Protein: True).' % (rje.iLen(gdb.index('biotype')[gtype]),gtype))
                else: self.printLog('#TYPE','%s %s genes (Protein: False).' % (rje.iLen(gdb.index('biotype')[gtype]),gtype))
            ## ~ [1b] ~ Update gene details from Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb.addField('sourcedb'); gdb.addField('sourceacc')
            for entry in gdb.entries():
                if entry['description'] == '\\N': entry['description'] = 'No description'
                else:
                    match = rje.matchExp('\[Source:(\S.+);Acc:(\S+)\]',entry['description'])
                    if match:
                        entry['sourcedb'] = match[0]
                        if match[0] == 'HGNC Symbol':  entry['sourceacc'] = 'HGNC:%s' % match[1]; entry['sourcedb'] = 'HGNC'
                        else: entry['sourceacc'] = match[1]
                        if entry['sourcedb'] == 'UniProtKB/TrEMBL': entry['sourcedb'] = 'Uniprot'
            #!# Two columns before the KNOWN/NOVEL gives an xref ID that matches the xref.txt table
            #!# This in turn has the xref primary accession number as column 3. (Column 1 is xrefid)
            #!# Uniprot external DB references are 2000,2001,2010,2200,2201,2202 and 2250
            #!# Would it be easier to parse EnsG-Uniprot links (and gene names?) directly from Uniprot?
            #!# This data should now be available through uniprotMapDB(), which gets uniprot and possible Symbol links.
            return True
        except:
            self.errorLog('parseGeneDetails failed. Will try old method.')
            return self.OLDparseGeneDetails(gene_table)
#########################################################################################################################
    def OLDparseGeneDetails(self,gene_table):  ### Parses out gene IDs from EnsEMBL Gene Table
        '''Parses out gene and protein IDs from EnsEMBL Gene Table.'''
        genes = self.loadFromFile(gene_table)
        gtypes = {}; no_match = []
        for line in genes:
            try:
                gtype = rje.matchExp('^(\d+)\s+(\S+)',line)[1]
                if gtype not in gtypes: gtypes[gtype] = 0
                gtypes[gtype] += 1
            except: continue
            for replace_type in ['C','V','J','D']:
                line = rje.replace(line,'IG_%s_gene' % replace_type,'protein_coding')
                line = rje.replace(line,'TR_%s_gene' % replace_type,'protein_coding')
            match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s*(\S.*)\[Source:(\S.+);Acc:(\S+)\]',line)

            #!# Two columns before the KNOWN/NOVEL gives an xref ID that matches the xref.txt table
            #!# This in turn has the xref primary accession number as column 3. (Column 1 is xrefid)
            #!# Uniprot external DB references are 2000,2001,2010,2200,2201,2202 and 2250
            #!# Would it be easier to parse EnsG-Uniprot links (and gene names?) directly from Uniprot?
            #!# This data should now be available through uniprotMapDB(), which gets uniprot and possible Symbol links.

            if match:
                if match[-2] == 'HGNC Symbol': self.dict['GeneAcc'][match[0]] = 'HGNC:%s' % match[-1]
                else: self.dict['GeneAcc'][match[0]] = match[-1]
            else:
                match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s*(\S.*)\[.+Acc:(\S+)\]',line)
                if match: self.dict['GeneAcc'][match[0]] = match[-1]
            #self.deBug(rje.split(line,'\t'))
            if len(rje.split(line,'\t')) > 14: ens_id = rje.split(line,'\t')[14]; line = rje.join(rje.split(line,'\t')[:14],'\t')
            else: ens_id = None
            #self.deBug(ens_id)
            if not match: match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s*(\S.+)\s+\d+\s+\d+\s+transcript with longest CDS',line)
            if not match: match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s+\d+\s+\d+\s+transcript with longest CDS',line)
            if not match: match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s*(\S.+)\s+\d+$',line)
            if not match: match = rje.matchExp('^(\d+)\s+protein_coding\s.+(KNOWN|NOVEL|PUTATIVE)\s*(\S.*)$',line)
            if not match:
                match = rje.matchExp('^(\d+)\s+protein_coding\s.+(\S.*)$',line)
                if match: match = match + ('Unknown',)
            if not match and line.find('protein_coding') > 0:
                self.log.errorLog('Cannot parse: %s' % line,printerror=False)       #!# Add PUTATIVE #!#
                self.deBug(line)
            if match:
                try: self.dict['GeneDesc'][match[0]] = match[2]
                except: self.dict['GeneDesc'][match[0]] = 'No description'
                self.printLog('\r#GENE','Gene details parsed for %s genes.' % rje.integerString(len(self.dict['GeneDesc'])),log=False,newline=False)
                if ens_id: self.dict['GeneMap'][ens_id] = match[0]
            elif ens_id and line.find('protein_coding') > 0: no_match.append(ens_id)
        #self.deBug(no_match)
        self.printLog('\r#GENE','Gene details parsed for %s genes.' % rje.iLen(self.dict['GeneDesc']))
        self.printLog('\r#GENE','Gene IDs mapped for %s protein coding genes.' % rje.iLen(self.dict['GeneMap']))
        self.printLog('\r#GENE','%d Gene Types (only protein_coding, IG & TR genes used).' % len(gtypes))
        self.printLog('#ENSID','%s Ensembl IDs without successful matching gene ID' % rje.iLen(no_match))
        for gtype in rje.sortKeys(gtypes):
            if rje.matchExp('TR_(\S)_gene',gtype) or rje.matchExp('IG_(\S)_gene',gtype) or 'protein_coding' in gtype: self.printLog('#GTYPE','%s %s: %s (True)' % (self.speciesCode(),gtype,rje.iStr(gtypes[gtype])))
            else: self.printLog('#GTYPE','%s %s: %s (False)' % (self.speciesCode(),gtype,rje.iStr(gtypes[gtype])))
#########################################################################################################################
    ### <5> ### EnsLoci processing methods                                                                              #
#########################################################################################################################
    def ensLociFas(self,species=None,spcode=None,read=False):   ### Returns path to EnsLoci fasta file.
        '''
        Returns path to EnsLoci fasta file.
        >> species:str [None] = Ensembl species. Will use self.getStr('Name') if None.
        >> read:bool [False] = whether file is being located for reading. (Will look in several places.)
        << ensloci:str = Returns filename.
        '''
        if not species:
            if spcode: species = self.getSpecies(spcode)
            else: species = self.getStr('Name')
        if not spcode: spcode = self.speciesCode(species)
        efile = 'ens_%s.loci.fas' % spcode
        ensloci = rje.makePath('%sensloci/%s/' % (self.getStr('EnsPath'),self.ensType(species))) + efile
        if read:
            if not rje.exists(ensloci):
                for epath in [rje.makePath('%sensloci/' % (self.getStr('EnsPath'))),
                              rje.makePath('%s%s/' % (self.getStr('EnsPath'),self.ensType(species))),
                              self.getStr('EnsPath')]:
                     if rje.exists(epath + efile): return epath + efile
        else: rje.mkDir(self,ensloci)
        return ensloci
#########################################################################################################################
    def ensLoci(self,species=''):   ### Performs the EnsLoci cleanup for given species
        '''Performs the EnsLoci cleanup for given species.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.force() and rje.exists(self.ensLociFas(species)):
                self.printLog('#LOCI','%s found! (force=F)' % self.ensLociFas(species))
                return False
            self.bugPrint(species)
            if not species: return False
            if not self.db('release'): self.ensemblRelease(save=False)
            ens = EnsEMBL(self.log,self.cmd_list)
            # Want parent tables but do not want to add gene mapping tables etc. to parent when created.
            ens.obj['DB'].list['Tables'] = self.db().tables()[0:]
            ens.setStr({'Name':species})
            self.dict['Stats'] = {}
            ens = ens.parseEnsEMBL(species)
            if ens == None: self.errorLog('Something disasterous happened during %s.parseEnsEMBL()' % species,printerror=False)
            if not ens: return False
            if ens.getBool('Parsed'): return True # Done already: no need
            ### ~ [1] Process Loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ens.ensLociSeq(): ens.summaryStats(); return True
        except: self.errorLog('Problem with rje_ensembl.ensLoci(%s)' % species)
        return False
#########################################################################################################################
    def ensLociProtein(self,entry): ### Determine the best protein for a given gene and update entry
        '''Determine the best protein for a given gene and update entry.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db('gene')       # 'ENSG' to other details, including 'gene_id','biotype','status','protein' (Y/N),'sourcedb' and 'sourceacc'
            pdb = self.db('peptide')    # ['ENSP','ENSG','Seq','NonX'],['ENSP']
            mdb = self.db('map')        # 'ENSG', 'Uniprot' and 'Symbol' and/or 'Gene' fields.
            hdb = self.db('hgnc')       # As map
            udb = self.db('unixref')    # ['Uniprot','UniprotID','ENSG','ENSP','Symbol','Gene'],['Uniprot'])
            udict = self.dict['UniprotObj'] # = {acc:UniProtObject}
            ensg = entry['ENSG']
            gentry = gdb.data(ensg)
            ## ~ [0a] ~ Generate protein lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ensp = pdb.index('ENSG')[ensg]      # List of ensembl peptide IDs
            if not ensp:
                self.warnLog('No Peptides parsed for %s (%s)' % (ensg,gentry['biotype']))
                entry['ensloci'] += '|NoPepSeq'
                self.bugPrint(gentry)
                self.debug()
            unip = []
            for db in [mdb,hdb,udb]:
                if not db: continue
                if 'ENSG' in db.fields() and ensg in db.index('ENSG',splitchar='|'):
                    for uacc in db.index('ENSG')[ensg]: unip.append(uacc)
                if 'ENSP' in db.fields():
                    for p in ensp:
                        if p not in db.index('ENSP',splitchar='|'): continue
                        for uacc in db.index('ENSP')[p]:  unip.append(uacc)
                if gentry['sourceacc'] and gentry['sourcedb'] in db.fields() and gentry['sourceacc'] in db.index(gentry['sourcedb'],splitchar='|'):
                    for uacc in db.index(gentry['sourcedb'])[gentry['sourceacc']]: unip.append(uacc)
            unip = rje.sortUnique(unip)
            ## ~ [0b] ~ Dump Trembl sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uobj = []
            if unip: entry['ensloci'] += '|Uniprot(%d)' % len(unip)
            for uacc in unip:
                if udict[uacc].reviewed() and udict[uacc] not in uobj: uobj.append(udict[uacc])
                elif self.dev() and udict[uacc].isSwissprot() and udict[uacc] not in uobj:
                    self.warnLog('Uniprot entry %s not reviewed but looks like SwissProt!' % udict[uacc].id())
                    uobj.append(udict[uacc])
            if uobj: entry['ensloci'] += '|Reviewed(%d)|SPROT' % len(uobj)
            elif unip: entry['ensloci'] += '|Reviewed(%d)|TREMBL' % len(uobj)

            ### ~ [1] ~ Allocate best protein sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Reviewed SwissProt is best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for u in uobj:
                if not u.isSpecies(self.speciesCode()):
                    self.warnLog('Uniprot entry %s is not %s!' % (u.id(),self.speciesCode()))
                    continue
                swiss_seq = u.seqi('Sequence').upper()
                for prot in ensp:
                    ## ~ Exactly matching sequence ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if pdb.data(prot)['Seq'] == swiss_seq:
                        entry['accnum'] = prot
                        entry['newname'] = '%s__%s' % (u.id(),u.accNum())
                        entry['newdesc'] = u.seqi('Description')
                        entry['ensloci'] += '|Match'
                        entry['seq'] = swiss_seq
                        return True
            if len(uobj) == 1:  # SeqMapper mapping of closest match to unambiguous Reviewed Uniprot
                if len(ensp) == 1: ensmap = ensp[0]     # Skip all the bother if only one protein!
                else:
                    ensmap = None
                    swiss_seq = uobj[0].seqi('Sequence').upper()
                    ## ~ Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    tmp = 'tmp_%s' % rje.randomString(8)
                    swissqry = '%s.qry.fas' % tmp
                    open(swissqry,'w').write('>%s\n%s\n' % (uobj[0].acc(),swiss_seq))
                    ENS = open('%s.fas' % tmp,'w')
                    for prot in ensp:
                        try: ENS.write('>%s\n%s\n' % (prot,pdb.data(prot)['Seq']))
                        except: self.errorLog('Problem saving sequence for %s' % prot); self.debug(pdb.data(prot))
                    ENS.close()
                    ## ~ Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    mapcmd = ['mapstat=id','minmap=90','automap=90'] + self.cmd_list + ['debug=F']
                    mapcmd += ['i=-1','v=-1','seqin=%s' % swissqry,'mapdb=%s.fas' % tmp,'startfrom=','resfile=%s' % tmp,'combine=T',
                               'gablamout=F','append=T','mapspec=None','mapping=Sequence,GABLAM','ordered=T','mapfocus=query',
                               'blastv=%d' % len(ensp)]
                    try:
                        self.log.opt['Silent'] = True
                        mapdict = seqmapper.SeqMapper(self.log,mapcmd).run(imenu=False,outputmap=False,returndict=True)
                        ensmap = mapdict[uobj[0].acc()]
                    except:
                        self.log.opt['Silent'] = False
                        self.errorLog('Problem with EnsEMBL SeqMapper(%s)' % uobj[0].acc())
                    self.log.opt['Silent'] = False
                    ## ~ Cleanup and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for tmpfile in glob.glob('%s*' % tmp): os.unlink(tmpfile)
                    if os.path.exists('None'): os.unlink('None')
                if ensmap:
                    entry['accnum'] = ensmap
                    entry['newname'] = '%s__%s' % (u.id(),ensmap)
                    entry['newdesc'] = u.seqi('Description')
                    entry['ensloci'] += '|SeqMap'
                    entry['seq'] = pdb.data(ensmap)['Seq']
                    return True
            ## ~ [1b] ~ Identify longest protein ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            badseq = 'X' * 20
            protid = None; badx = False
            for epept in ensp:
                if badseq in pdb.data(epept)['Seq']:
                    if uobj: continue
                    elif not protid: protid = epept; badx = True; continue
                    elif (protid and not badx): continue
                elif badx or not protid: protid = epept; badx = False; continue
                if pdb.data(epept)['NonX'] > pdb.data(protid): protid = epept
            ## ~ [1c] ~ No decent sequences without run of Xs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if uobj and badx or not protid:
                for u in uobj:
                    if u.seqobj().nonX() > (len(entry['seq']) - rje.count(entry['seq'],'X')):
                        entry['accnum'] = u.accNum()
                        entry['newname'] = '%s__%s' % (u.id(),u.accNum())
                        entry['newdesc'] = u.seqi('Description')
                        entry['ensloci'] += '|REPLACED'
                        entry['seq'] = u.seqi('Sequence').upper()
                if entry['seq']: return
            ## ~ [1d] ~ Use longest peptide ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            entry['accnum'] = protid
            entry['ensloci'] += '|Longest'
            entry['seq'] = pdb.data(protid)['Seq']
            ## ~ [1e] ~ Use SwissProt sequence information if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if uobj:
                entry['newname'] = '%s__%s' % (uobj[0].id(),protid)
                entry['newdesc'] = uobj[0].seqi('Description')
                return
            ## ~ [1f] ~ Construct new sequence name if no SwissProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if len(gentry['description']) >= 3 and gentry['description'] != 'No description': entry['newdesc'] = gentry['description']
            elif pdb.data(protid)['Desc']: entry['newdesc'] = pdb.data(protid)['Desc']
            else: entry['newdesc'] = 'Ensembl %s gene %s.' % (gentry['status'],ensg)
            gene = ''
            for db in [mdb,hdb,udb]:
                if not db: continue
                for field in ['Symbol','Gene']:
                    if gene and len(rje.split(gene)) == 1: break
                    try:
                        genes = db.indexDataList('ENSG',ensg,field)
                        for g1 in genes[0:]:
                            for g2 in genes[0:]:
                                if g1 == g2: continue
                                if g1 in g2 and g2 in genes: genes.remove(g2)
                                elif g2 in g1 and g1 in genes: genes.remove(g1)
                        gene = rje.join(genes)
                    except: pass
            if len(rje.split(gene)) > 1: self.warnLog('Bad gene: "%s"' % gene); gene = ''
            if gene.startswith(ensg) or not gene:
                try: gene = {'KNOWN':'ens','NOVEL':'nvl','PUTATIVE':'put','PREDICTED':'pred','KNOWN_BY_PROJECTION':'proj','UNKNOWN':'nvl', 'ANNOTATED':'ens'}[gentry['status']]
                except: gene = 'pep'
            entry['newname'] = '%s_%s__%s' % (gene.lower(),self.speciesCode(),protid)
        except: self.errorLog('Error in EnsEMBL.ensLociProtein(%s)' % self.getStr('Name'),quitchoice=False); entry['ensloci'] += '|ERROR'
#########################################################################################################################
    def ensLociSeq(self):  ### Produces a reduced EnsEMBL dataset of one protein (the longest known) per locus
        '''Produces a reduced EnsEMBL dataset of one protein (the longest known) per locus.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edb = self.db().addEmptyTable('ensloci',['ENSG','accnum','newname','newdesc','seq','ensloci'],['ENSG'])
            gdb = self.db('gene')       # 'ENSG' to other details, including 'gene_id','biotype','status','protein' (Y/N),'sourcedb' and 'sourceacc'
            pdb = self.db('peptide')    # ['ENSP','ENSG','Seq','NonX'],['ENSP']

            ### ~ [1] ~ Work through loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = gdb.entryNum()
            for ensg in gdb.dataKeys():
                self.progLog('\r#LOCI','Building EnsLoci sequences: %.2f%%' % (gx/gtot)); gx += 100.0
                entry = {'ENSG':ensg,'ensloci':gdb.data(ensg)['biotype']}
                ## ~ [1a] ~ Check for protein-coding ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not gdb.data(ensg)['protein'] == 'Y':
                    edb.addEntry(entry)
                    continue
                ## ~ [1b] ~ Get best protein sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.ensLociProtein(entry)
                if entry['newname'] in edb.index('newname'): entry['ensloci'] += 'REDUNDANT'
                edb.addEntry(entry)
                if self.dev() and gx <= 2000: self.bugPrint('\r>%s %s' % (entry['newname'],entry['newdesc']))
            self.printLog('\r#LOCI','Building %s EnsLoci sequences complete.' % edb.iNum())

            ### ~ [2] ~ Check all peptide genes have been processed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ensg in rje.sortKeys(pdb.index('ENSG')):
                if ensg and not edb.data(ensg): raise ValueError('Peptide gene %s not found in EnsLoci table!' % ensg)

            ### ~ [3] ~ Output EnsLoci data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Update stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for skey in ['SProt','TrEMBL','Known','Novel','Putative','Redundant']:
                try: self.dict['Stats'][skey] = len(edb.index('ensloci',splitchar='|')[skey.upper()])
                except: self.dict['Stats'][skey] = 0
            ## ~ [3b] ~ Remove redundant loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'REDUNDANT' in edb.index('ensloci',splitchar='|'):
                rx = 0
                for ensg in edb.index('ensloci')['REDUNDANT']: edb.data().pop(ensg); rx += 1
                self.printLog('#RED','%s REDUNDANT entries dropped: %s EnsLoci remain.' % (rx,edb.iNum()))
            ## ~ [3c] ~ Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            efas = self.ensLociFas()
            ENSLOCI = open(efas,'w'); ex = 0.0; etot = edb.entryNum()
            for ensg in edb.dataKeys():
                self.progLog('\r#FAS','Output EnsLoci sequences to %s: %.1f%%' % (efas,ex/etot)); ex += 100.0
                entry = edb.data(ensg)
                ENSLOCI.write('>%s %s\n%s\n' % (entry['newname'],entry['newdesc'],entry['seq']))
            ENSLOCI.close()
            self.printLog('\r#FAS','%s %s EnsLoci sequences output to %s' % (edb.iNum(),self.getStr('Name'),efas))
            ## ~ [3d] ~ Output ensloci stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            efile = rje.makePath('%s%s/%s' % (self.getStr('EnsPath'),self.ensType(),self.speciesCode())) + 'ensloci.tdt'
            edb.dropField('seq')
            edb.saveToFile(efile)
            return True
        except: self.errorLog('Error in EnsEMBL.ensLociSeq(%s)' % self.getStr('Name'),quitchoice=True)
        return False
#########################################################################################################################
    def bestProt(self,protlist,type,locus): ### Returns the best ProtID for the locus - matching sequence or longest
        '''
        Returns the best ProtID for the locus - matching sequence or longest.
        >> protlist:list of protein IDs
        >> type:str = known/novel/putative
        >> locus:GeneID
        '''
        try:
            ## Match by sequence to SwissProt? ##
            try: acc = self.dict['GeneAcc'][self.dict['GeneMap'][locus]]
            except: acc = '-'
            h2u = self.uniFromHGNC(acc); mapacc = self.uniFromMap(locus)
            if h2u: acc = h2u
            elif mapacc: acc = mapacc
            if acc in self.dict['AccSeq'].keys():
                swiss_seq = self.dict['AccSeq'][acc]
                #self.deBug('%s => %s = %s' % (locus,acc,swiss_seq))
                if swiss_seq:
                    try: (gene,species,newacc,desc) = rje.matchExp('^(\S+)_(\S+)__(\S+)\s*(\S*.*)$',self.dict['AccName'][acc])
                    except:
                        self.log.errorLog('Problem with RegExp for %s' % self.dict['AccName'][acc])
                        protid = self.longestProt(protlist,type)
                        self.dict['NewAccName'][protid] = '%s_%s__%s %s' % (type,self.speciesCode(),protid,self.dict['GeneDesc'][self.dict['GeneMap'][locus]])
                        return protid
                    #X#self.verbose(0,4,'%s => %s = %s' % (locus,acc,self.dict['AccName'][acc]),1)
                    if gene == newacc:  # TrEMBL
                        protid = self.longestProt(protlist,type)
                        if self.dict['ProtSeq'][protid] == swiss_seq:
                            self.dict['NewAccName'][protid] = self.dict['AccName'][acc]
                            return protid
                        self.dict['NewAccName'][protid] = '%s_%s__%s %s' % (gene,self.speciesCode(),protid,self.dict['GeneDesc'][self.dict['GeneMap'][locus]])
                        return protid
                    elif gene == gene.upper():   # SwissProt
                        for prot in protlist:
                            #self.deBug('%s\n vs \n%s\n' % (swiss_seq,self.dict['ProtSeq'][prot]))
                            if self.dict['ProtSeq'][prot] == swiss_seq:
                                #self.deBug('Best=%s' % prot)
                                self.dict['NewAccName'][prot] = self.dict['AccName'][acc]
                                return prot
                            #else: self.deBug('Not %s!' % prot)
                        # No exact sequence match: try SeqMapper, else Longest Protein #
                        protid = self.seqMapper(swiss_seq,acc,protlist,type)
                        #self.deBug('Mapped %s = %s' % (protlist,protid))
                        if not protid:
                            protid = self.longestProt(protlist,type)
                            #self.deBug('Longest %s = %s' % (protlist,protid))
                        # No exact match - change accnum to protid! #
                        self.dict['NewAccName'][protid] = '%s_%s__%s %s' % (gene,self.speciesCode(),protid,desc)
                        #self.deBug('=> %s' % self.dict['NewAccName'][protid])
                        #!# Fudge for dodgy EnsEMBL sequences. :o( #!#
                        if self.dict['ProtSeq'][protid].find('X' * 20) >= 0:    # Possible shite sequence, even though best #
                            p = len(self.dict['ProtSeq'][protid]) - rje.count(self.dict['ProtSeq'][protid],'X')
                            s = len(swiss_seq) - rje.count(swiss_seq,'X')
                            if s >= p:  # Replace crap EnsEMBL sequence with SwissProt
                                self.dict['ProtSeq'][protid] = swiss_seq
                                self.log.printLog('#REPLACE','Replaced rubbish %s sequence with SwissProt %s sequence!' % (protid,acc))
                                self.dict['NewAccName'][protid] = '%s_%s__%s !REPLACED! %s' % (gene,self.speciesCode(),acc,desc)
                        return protid
            return self.longestProt(protlist,type)
        except:
            self.log.errorLog('Error in EnsEMBL.bestProt(%s)' % locus,printerror=True,quitchoice=False)
            return None
#########################################################################################################################
    def seqMapper(self,swiss_seq,acc,protlist,type):    ### Returns the ProtID for the best match to swiss_seq in protlist
        '''Returns the ProtID for the best match to swiss_seq in protlist.'''
        try:
            ### Skip all the bother if only one protein! ###
            if len(protlist) == 1: return protlist[0]
            
            ### Setup Files ###
            tmp = 'tmp_%s' % rje.randomString(8)
            swissqry = '%s.qry.fas' % tmp
            open(swissqry,'w').write('>%s\n%s\n' % (acc,swiss_seq))
            ENS = open('%s.fas' % tmp,'w')
            for prot in protlist: ENS.write('>%s\n%s\n' % (prot,self.dict['ProtSeq'][prot]))
            ENS.close()

            ### Mapping ###
            mapcmd = ['mapstat=id','minmap=80','automap=80'] + self.cmd_list
            mapcmd += ['i=-1','v=-1','seqin=%s' % swissqry,'mapdb=%s.fas' % tmp,'startfrom=','resfile=%s' % tmp,'combine=T',
                       'gablamout=F','append=T','mapspec=None','mapping=Sequence,GABLAM','ordered=T','mapfocus=query',
                       'blastv=%d' % len(protlist)]
            try:
                self.log.opt['Silent'] = True
                mapdict = seqmapper.SeqMapper(self.log,mapcmd).run(imenu=False,outputmap=False,returndict=True)
                #self.deBug(mapdict)
            except:
                self.log.opt['Silent'] = False
                self.log.errorLog('Problem with EnsEMBL SeqMapper(%s)' % acc)
                return None
            self.log.opt['Silent'] = False

            ### Cleanup and finish ###
            for tmpfile in glob.glob('%s*' % tmp): os.unlink(tmpfile)
            #self.deBug(mapdict)
            return mapdict[acc]
        except:
            self.log.errorLog('Error in EnsEMBL.seqMapper(%s)' % protlist,printerror=True,quitchoice=False)
            return None
#########################################################################################################################
    def longestProt(self,protlist,type):    ### Returns the ProtID for the longest prot in protlist
        '''Returns the ProtID for the longest prot in protlist.'''
        try:
            longest = (protlist[0],self.dict['ProtLen'][type][protlist[0]])
            for prot in protlist[1:]:
                if self.dict['ProtLen'][type][prot] > longest[1]: longest = (prot,self.dict['ProtLen'][type][prot])
            return longest[0]
        except:
            self.log.errorLog('Error in EnsEMBL.longestProt(%s)' % protlist,printerror=True,quitchoice=False)
            return None
#########################################################################################################################
    def newName(self,geneid,protid,type):    ### Returns new name for protein based on GeneID and ProtID
        '''Returns new name for protein based on GeneID and ProtID.'''
        try:
            ### AccNum ###
            try:
                genenum = self.dict['GeneMap'][geneid]      # Clears to save memory. Change if want to keep data
                accnum = self.dict['GeneAcc'][genenum]
                geneid = rje.split(geneid,':')[-1]; protid = rje.split(protid,':')[-1]; accnum = rje.split(accnum,':')[-1]
                symbol = self.symbolFromHGNC(accnum)
                descacc = self.uniFromMap(geneid)
                if symbol: accnum = self.uniFromHGNC(accnum); descacc = accnum
                if not symbol: symbol = self.symbolFromMap(geneid)
                if len(accnum) < 3: accnum = protid
                if not descacc: descacc = accnum
                #self.deBug('%s = ID %s; Acc %s; Symbol %s' % (geneid,genenum,accnum,symbol))
            except:
                self.errorLog('Problem getting AccNum/GeneNum from GeneID %s' % geneid)
                accnum = protid; genenum = None
            ## Mapped AccNum but not kept sequence ##
            if protid in self.dict['NewAccName'] and not rje.matchExp('^\S+_\S+__(%s)' % accnum,self.dict['NewAccName'][protid]):   # Primary != secondary
                newacc = rje.matchExp('^(\S+)_(\S+)__(\S+)',self.dict['NewAccName'][protid])[2]
            ## Mapped Primary AccNum ##
            elif protid in self.dict['NewAccName']: newacc = accnum
            ## Could not map AccNum ##
            else:
                newacc = protid
                if accnum != protid: type = 'refseq'
            if newacc == 'pep': newacc = protid = geneid; type = 'pep'  #!# Needs better handling! #!#
            #self.deBug('%s: %s = %s => %s =>> %s' % (type,geneid,protid,accnum,newacc))
            newacc = rje.replace(newacc,'__','-')
            ## Check for double-mapping ##
            if self.dict['GeneAcc'].values().count(newacc) > 1:
                self.log.printLog('#MAP','AccNum %s mapped to %d Genes, including %s (%s)' % (newacc,self.dict['GeneAcc'].values().count(newacc),geneid,protid),screen=False)
                newacc = rje.replace(protid,'__','-')
            ## Special ##
            if self.speciesCode() == 'YEAST': newacc = protid

            ### Description ###
            try: desc = self.dict['GeneDesc'][genenum]
            except: self.errorLog('Problem getting Desc for %s' % accnum); desc = ''
            if len(desc) < 3 and protid in self.dict['NewAccName']:
                (gene,species,descacc,desc) = rje.matchExp('^(\S+)_(\S+)__(\S+)\s*(\S*.*)$',self.dict['NewAccName'][protid])
            if len(desc) < 3: desc = 'EnsEMBL %s' % type
            if self.speciesCode() == 'CAEEL' and not symbol:   # and geneid == protid:
                symbol = rje.split(desc)[0]
                if symbol.find('-') < 0: symbol = ''
            if symbol: desc = rje.join(rje.split('%s [acc:%s pep:%s symbol:%s gene:%s]' % (desc,descacc,protid,symbol,geneid)))
            else: desc = rje.join(rje.split('%s [acc:%s pep:%s gene:%s]' % (desc,descacc,protid,geneid)))
                
            ### Default ID ###
            if symbol: id = '%s_%s__%s' % (symbol.lower(),self.speciesCode(),newacc)
            elif type == 'refseq': id = 'ref_%s__%s' % (self.speciesCode(),newacc)
            elif type == 'known': id = 'ens_%s__%s' % (self.speciesCode(),newacc)
            elif type == 'putative': id = 'put_%s__%s' % (self.speciesCode(),newacc)
            elif type == 'pep': id = 'pep_%s__%s' % (self.speciesCode(),newacc)
            else: id = 'nvl_%s__%s' % (self.speciesCode(),newacc)

            ### Name ###
            if protid in self.dict['NewAccName']:
                (gene,species,descacc,newdesc) = rje.matchExp('^(\S+)_(\S+)__(\S+)\s*(\S*.*)$',self.dict['NewAccName'][protid])
                if species != self.speciesCode():
                    desc = rje.replace(desc,'acc:%s' % accnum,'acc:%s !%s!' % (accnum,species))
                    id = rje.replace(id,newacc,protid)
                    return '%s %s' % (id,desc)
                if newacc == gene:  # TrEMBL
                    return '%s %s' % (id,desc)
                return '%s_%s__%s %s' % (gene,species,newacc,desc)
            return '%s %s' % (id,desc)
        except:
            self.log.errorLog('Error in EnsEMBL.newName(%s,%s)' % (geneid,protid))
            return protid
#########################################################################################################################
    ### <6> ### EnsEMBL reformatting from rje_dbase                                                                     #
#########################################################################################################################
    def makeEnsPepFasta(self):  ### Generates reformatted fasta files of Ensembl peptides
        '''Generates reformatted fasta files of Ensembl peptides.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            rdb = self.db('release')
            if not rdb: rdb = self.ensemblReleaseTaxonomy() # Downloads release data and extracts species information
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ GENERATE PEPTIDE FASTA ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            pepdir = rje.makePath('%senspep/' % (self.getStr('EnsPath')))
            rje.mkDir(self,pepdir)
            ### ~ [2] Generate new fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for section in rdb.indexKeys('section'):
                pepfas = '%sens_%s.pep.fas' % (pepdir,section)
                if rje.exists(pepfas) and not self.force():
                    self.printLog('#FAS','%s found!' % pepfas)
                    continue
                PEPFAS = open(pepfas,'w')
                for entry in rdb.indexEntries('section',section):
                    species = entry['species']
                    spcode = entry['spcode']
                    enspep = glob.glob('%s%s/%s/%s.*.pep.all.fa' % (self.getStr('EnsPath'),section,spcode,species))
                    if not enspep: self.errorLog('Cannot find %s%s/%s/%s.*.pep.all.fa' % (self.getStr('EnsPath'),section,spcode,species)); continue
                    enspep = enspep[0]
                    ## ~ [2a] Load and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    scmd = ['seqin=%s' % enspep,'seqmode=file','seqindex=F','autoload=F']
                    enseq = rje_seqlist.SeqList(self.log,self.cmd_list+scmd)
                    enseq.loadSeq()
                    ex = nx = px = sx = 0
                    while enseq.nextSeq():
                        (name,sequence) = enseq.currSeq()
                        if 'pep:known' in name: gene = 'ens'; ex += 1
                        elif 'pep:novel' in name: gene = 'nvl'; nx += 1
                        else: gene = 'pep'; px += 1
                        PEPFAS.write('>%s_%s__%s\n%s\n' % (gene,spcode,name,sequence)); sx += 1
                    self.printLog('#FAS','%s %s peptides output to %s: %s known; %s novel; %s other.' % (rje.iStr(sx),species,pepfas,rje.iStr(ex),rje.iStr(nx),rje.iStr(px)))
                PEPFAS.close()
        except: self.errorLog('Problem with makeEnsPepFasta')
#########################################################################################################################
    def reformatDB(self,speclist=[],usetypes=['known','known-ccds','novel','abinitio']):   ### This method is for reformatting the EnsEMBL fasta files into a single file
        '''
        This method is for reformatting the EnsEMBL fasta files into a single file.
        >> speclist: list of species to use (keys of self.dict['EnsSpec'])
        >> usetypes: list of types to use ['known','known-ccds','novel','abinitio']
        '''
        try:
            ### Setup ###
            if not speclist: speclist = rje.sortKeys(self.dict['EnsSpec'])
            enspath = rje.makePath(self.info['EnsPath'] + self.dict['EnsSpec'][species])

            ### Make Combined File ###
            for species in speclist:
                ## File names ##
                newens = enspath + 'ens_%s.new.fas' % self.dict['EnsSpec'][species]
                outens = enspath + 'ens_%s.all.fas' % self.dict['EnsSpec'][species]
                ensfiles = []  
                for utype in usetypes:
                    for ufile in glob.glob('%s*.fa' % enspath):
                        if ufile.find(utype) > 0:
                            ensfiles.append(ufile)
                            break
                #X#self.deBug(ensfiles)
                ## Check Dates ##
                skip = not self.opt['Force']
                for ens in ensfiles:
                    if rje.isYounger(newens,ens) != newens: skip = False
                if skip:
                    self.log.printLog('#ENS','%s already exists and younger than components (Force=F)' % (newens))
                    continue    # Next species
                ## Process ##
                NEW = open(newens,'w')
                for ens in ensfiles:
                    INFILE = open(ens,'r')
                    sx = 0
                    name = ''
                    sequence = ''
                    nextname = ''
                    while 1:
                        line = INFILE.readline()
                        if not line:
                            if name:
                                NEW.write('>%s\n%s\n' % (name,sequence))
                                sx += 1
                            break
                        line = rje.chomp(line)
                        if line[:1] == '>':     # Description line
                            if name:
                                NEW.write('>%s\n%s\n' % (name,sequence))
                                sx += 1
                                self.log.printLog('\r#SEQ','%s sequences reformatted from %s.' % (rje.integerString(sx),ens),log=False,newline=False)
                            # Reformat
                            if ens.find('known') > 0: name = 'ens_%s__%s' % (species,line[1:])
                            elif ens.find('novel') > 0: name = 'nvl_%s__%s' % (species,line[1:])
                            else: name = 'scan_%s__%s' % (species,line[1:])
                            sequence = ''
                        else: sequence = '%s%s' % (sequence,line)
                    INFILE.close()
                    self.log.printLog('\r#SEQ','%s sequences reformatted from %s.' % (rje.integerString(sx),ens))
                NEW.close()
                seqcmd = ['reformat=fas','seqout=%s' % outens,'seqin=%s' % newens,'memsaver=T','autoload=T','autofilter=T']
                rje_seq.SeqList(self.log,self.cmd_list+seqcmd)
                os.unlink(newens)
        except: self.log.errorLog('Error in EnsEMBL.reformatDB()',printerror=True)
#########################################################################################################################
    ### <7> ### EnsGO Gene Ontology Datasets
#########################################################################################################################
    def ensGO(self):    ### Parses GO data and generates EnsEMBL GO Datasets.
        '''Parses GO data and generates EnsEMBL GO Datasets from EnsEMBL Loci data and BioMart download.'''
        try:
            ### ~ [1] Setup files and GO dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['EnsGO']: return False     # No species to parse
            self.dict['GO'] = {'GTypes':{'C':'CC','P':'BP','F':'MF','X':'X'},
                               'ID':{'C':[],'P':[],'F':[],'X':[]},         # Type:[IDs]
                               'Desc':{},       # ID:Desc
                               'Obselete':[]}   # List of obselete GO Terms

            ### ~ [2] Parse GO data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.parseGOIDs()
           
            ### ~ [3] Read in EnsEMBL GO downloads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['EnsGO']:
                ## Setup ##
                sfile = 'ens_%s.GO_summary.tdt' % spec
                headers = ['GO ID','GO Desc','GO Type','GeneNum','ProtNum']
                rje.delimitedFileOutput(self,sfile,headers,rje_backup=True)
                self.dict['GO']['EnsGene'] = {}         # Dictionary of GO:[GeneIDs]
                #X#self.dict['GO']['AltGene'] = {}     # Alternative Gene ID for each EnsEMBL gene
                #X#self.dict['GO']['Entrez'] = {}      # Entrez Gene ID for EnsEMBL gene (if downloaded from BioMart)
                ## Parse details ##
                self.parseEnsGO(spec)
                self.parseEnsEMBL(self.getSpecies(spec),False)
                ix = 0
                ## Create files ##
                for type in self.dict['GO']['ID']:
                    if not self.dict['GO']['ID'][type]: continue
                    godir = rje.makePath('%s_%s/' % (spec,self.dict['GO']['GTypes'][type]))
                    rje.mkDir(self,godir)
                    for id in self.dict['GO']['ID'][type]:
                        #self.deBug(id)
                        if not id in self.dict['GO']['EnsGene']: continue
                        if self.stat['MinGO'] > len(self.dict['GO']['EnsGene'][id]): continue
                        GOFAS = open('%s%s.fas' % (godir,id),'w')
                        (gx,px) = (0,0)
                        for gene in self.dict['GO']['EnsGene'][id]:
                            for ens in self.dict['Loci']:
                                if gene in self.dict['Loci'][ens]: break
                            e = {'known':'ens','novel':'nvl'}[ens]
                            gx += 1
                            for prot in self.dict['Loci'][ens][gene]:
                                try: GOFAS.write('>%s_%s__%s %s\n%s\n' % (e,spec,prot,self.dict['GeneDesc'][self.dict['GeneMap'][gene]],self.dict['ProtSeq'][prot]))
                                except: GOFAS.write('>%s_%s__%s %s\n%s\n' % (e,spec,prot,gene,self.dict['ProtSeq'][prot]))
                                px += 1
                        GOFAS.close()
                        ## Output summary ##
                        gdict = {'GO ID':rje.replace(id,'_',':'),'GO Desc':self.dict['GO']['Desc'][id],'GO Type':type,'GeneNum':gx,'ProtNum':px}
                        rje.delimitedFileOutput(self,sfile,headers,datadict=gdict)
                        ix += 1
                        self.log.printLog('\r#GO','Making %s %s GO datasets.' % (rje.integerString(ix),spec),log=False,newline=False)
                self.log.printLog('\r#GO','Making %s %s GO datasets complete.' % (rje.integerString(ix),spec))

        except: self.log.errorLog('Wah wah wah! Cry like a baby - ensLociGo() is busted.')
#########################################################################################################################
    def parseGOIDs(self):   ### Parses general GO ID data in self.dict['GO']
        '''Parses general GO ID data in self.dict['GO'].'''
        #?# os.system('wget http://www.geneontology.org/doc/GO.terms_ids_obs') #?#
        if os.path.exists(self.info['GOIDs']):  go_lines = self.loadFromFile(self.info['GOIDs'])
        else: go_lines = urllib.urlopen('http://www.geneontology.org/doc/GO.terms_ids_obs').readlines()
        if not go_lines: self.log.errorLog('Cannot find or download GO.terms_ids_obs file')
        else:
            lx = 0.0
            for line in go_lines:
                lx += 100.0
                if line[:1] == '!': continue
                go = rje.split(rje.chomp(line),'\t')
                if len(go) >= 3:
                    ## Get GO Information #
                    desc = go[1]
                    id = rje.replace(go[0],':','_')
                    type = go[2]
                    if go[-1] == 'obs':
                        if id not in self.dict['GO']['Obselete']: self.dict['GO']['Obselete'].append(id)
                        if not self.opt['ObsGO']: continue
                    if type not in ['P','C','F']:
                        self.log.errorLog('GO Type Error. What is "%s"?' % line,printerror=False)
                        type = 'X'
                    ## Update dictionary ##
                    self.dict['GO']['ID'][type].append(id)
                    self.dict['GO']['Desc'][id] = desc
                self.log.printLog('\r#GO','Reading GO IDs %.1f%%' % (lx/len(go_lines)),log=False,newline=False)
            gtxt = ''
            for type in ['C','P','F']: gtxt += '%s %s; ' % (rje.integerString(len(self.dict['GO']['ID'][type])),self.dict['GO']['GTypes'][type])
            self.log.printLog('\r#GO',gtxt)
            if self.opt['ObsGO']: self.log.printLog('#OBS','%s Obselete terms included' % rje.integerString(len(self.dict['GO']['Obselete'])))
            else: self.log.printLog('#OBS','%s Obselete terms excluded' % rje.integerString(len(self.dict['GO']['Obselete'])))
#########################################################################################################################
    def parseEnsGO(self,spec):  ### Parses data from EnsEMBL Mart download
        '''Parses data from EnsEMBL Mart download.'''
        ### ~ Setup filename for EnsEMBL BioMart download of GO mappings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        efile = 'ens_%s.GO.tdt' % spec
        if not os.path.exists(efile): efile = 'ens_%s.GO.csv' % spec
        if not os.path.exists(efile):
            self.log.errorLog('Cannot find %s (or tdt)' % efile,printerror=False)
            return

        ### ~ Read data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        godata = rje.dataDict(self,efile,['Ensembl Gene ID','GO ID','GO evidence code'])
        delimit = rje.delimitFromExt(filename=efile)

        ### ~ Process data to generate GO dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        gx = 0.0
        for go in godata:
            self.log.printLog('\r#GO','Processing %s: %.1f%%' % (efile,gx/len(godata)),log=False,newline=False)
            gx += 100.0
            ## ~ Get Gene, GO ID and GO evidence code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try: (gene,id,evidence) = rje.split(go,delimit)
            except: continue    # no GO!
            ## ~ Check evidence code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['GOEvidence'] and evidence not in self.list['GOEvidence']: continue
            ## ~ Reformat and check ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            id = rje.replace(id,':','_')
            if not id or (id in self.dict['GO']['Desc'] and self.dict['GO']['Desc'][id] in ['Gene_Ontology','molecular_function','cellular_component','biological_process']): continue    # Too broad
            ## ~ Update dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if id not in self.dict['GO']['EnsGene']: self.dict['GO']['EnsGene'][id] = []
            self.dict['GO']['EnsGene'][id].append(gene)
        self.log.printLog('\r#GO','Processing complete: %s %s genes with GO annotation.' % (rje.integerString(len(self.dict['GO']['EnsGene'])),spec))
#########################################################################################################################
    def OLDensGO(self):    ### Parses Gene Ontology Datasets
        '''Parses Gene Ontology Datasets from EnsEMBL Loci data and BioMart download.'''
        try:
            ### Setup ###
            if not self.list['EnsGO']:      # No species to parse
                return False

            ### Read in GO dictionary ###
            gtypes = {'C':'CC','P':'BP','F':'MF'}
            go_desc = {}    # Type: {ID:Desc}
            go_type = {}    # ID: Type
            go_obs = []     # List of obselete GO Terms
            go_lines = self.loadFromFile(rje.makePath('EnsGO/GO.terms_ids_obs',True),chomplines=True)
            gx = {'P':0,'C':0,'F':0}            
            lx = 0.0
            for line in go_lines:
                lx += 100.0
                if line[:1] == '!':
                    continue
                go = rje.split(line,'\t')
                if len(go) >= 3:
                    desc = go[1]
                    id = rje.replace(go[0],':','_')
                    type = go[2]
                    if go[-1] == 'obs':
                        if id not in go_obs:
                            go_obs.append(id)
                        if not self.opt['ObsGO']:
                            continue
                    if type not in ['P','C','F']:
                        self.log.errorLog('GO Type Error. What is "%s"?' % line,printerror=False)
                        continue
                    go_desc[id] = desc
                    go_type[id] = gtypes[type]
                    gx[type] += 1
                self.log.printLog('\r#GO','%s CC; %s MF; %s BP; %.1f%%' % (rje.integerString(gx['C']),rje.integerString(gx['F']),rje.integerString(gx['P']),lx/len(go_lines)),log=False,newline=False)
            self.log.printLog('\r#GO','%s CC; %s MF; %s BP; %.1f%%' % (rje.integerString(gx['C']),rje.integerString(gx['F']),rje.integerString(gx['P']),lx/len(go_lines)))
            if self.opt['ObsGO']:
                self.log.printLog('#OBS','%s Obselete terms included' % rje.integerString(len(go_obs)))
            else:
                self.log.printLog('#OBS','%s Obselete terms excluded' % rje.integerString(len(go_obs)))

            ### Make GO Desc DO File for STATA ###
            DO = open(rje.makePath('EnsGO/go_map.do',True),'w')
            gx = 0
            for go in rje.sortKeys(go_desc):
                gx += 1
                DO.write('replace desc = "%s" if dataset == "%s" /* %s */\n' % (go_desc[go],go,go_type[go]))
                self.log.printLog('\r#GO','%s GO terms.' % rje.integerString(gx),newline=False,log=False)
            self.log.printLog('\r#GO','%s GO terms.' % rje.integerString(gx))
            DO.close()

            ### Read in EnsLoci GO dictionary ###
            for spec in self.list['EnsGO']:
                gene_go = {}     # Dictionary of gene:[gocodelist]
                golines = self.loadFromFile(rje.makePath('EnsGO/ens_%s.GO.csv' % spec,True))[1:]
                lx = 0.0
                for line in golines:
                    lx += 100.0
                    go = rje.split(line,',')
                    if len(go) >= 5:
                        gene = go[0]
                        id = rje.replace(go[2],':','_')
                        if not id or (id in go_desc and go_desc[id] in ['Gene_Ontology','molecular_function','cellular_component','biological_process']):
                            continue
                        if gene not in gene_go:
                            gene_go[gene] = []
                        if id not in gene_go[gene]:
                            gene_go[gene].append(id)
                    self.log.printLog('\r#GO','%s Genes; %.1f%%' % (rje.integerString(len(gene_go)),lx/len(golines)),log=False,newline=False)
                self.log.printLog('\r#GO','%s Genes.' % rje.integerString(len(gene_go)))

                ### Make GO Files ###
                locifas = self.ensLociFas(spcode=spec,read=True)
                self.devPrint(locifas)
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                seqx = rje_seq.SeqCount(self,locifas)
                SEQFILE = open(locifas, 'r')
                lastline = ''
                sx = 0.0
                fx = {}     # Dictionary of {file:seqx}
                missing = []      # Missing GO_IDs
                ox = 0      # Obselete GO_IDs
                gx = 0      # Counter
                while 1:
                    (seq, lastline) = seqlist.nextFasSeq(fileobject=SEQFILE,lastline=lastline)
                    if not seq:
                        break
                    seqlist.seq = [seq]
                    gene = rje.matchExp('gene:(\S+)\]',seq.info['Name'])[0]
                    #self.deBug(gene)
                    if gene in gene_go:
                        for go in gene_go.pop(gene):
                            #self.deBug(go)
                            if go not in go_type:
                                if go in go_obs:
                                    ox += 1
                                elif go not in missing:
                                    missing.append(go)
                                continue
                            godir = rje.makePath('EnsGO/%s_%s/' % (spec,go_type[go]))
                            if not os.path.exists(godir):
                                os.mkdir(godir)
                            gofile = rje.makePath('EnsGO/%s_%s/%s.fas' % (spec,go_type[go],go),True)
                            if gofile in fx:
                                fx[gofile] += 1
                            else:
                                fx[gofile] = 1
                                if os.path.exists(gofile):
                                    os.unlink(gofile)
                            open(gofile,'a').write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
                            gx += 1
                            self.log.printLog('\r#GO','Making GO datasets: %s; %.1f%%' % (rje.integerString(gx),(sx/seqx)),log=False,newline=False)
                    sx += 100.0
                    self.log.printLog('\r#GO','Making GO datasets: %s; %.1f%%' % (rje.integerString(gx),(sx/seqx)),log=False,newline=False)
                self.log.printLog('\r#GO','Making GO datasets complete: %s GO-gene pairs' % rje.integerString(gx))
                SEQFILE.close()

                if self.opt['ObsGO']:
                    self.log.printLog('#OBS','%s obselete GO-Gene pairs;' % (rje.integerString(ox)))
                self.log.printLog('#MISS','%s GO-Gene pairs missing.' % (rje.integerString(len(missing))))
                if missing:
                    open(rje.makePath('EnsGO/%s.missing_GO.fas' % spec,True),'w').write(rje.join(missing+[''],'\n'))

                ### MinGO ###
                if self.stat['MinGO'] > 0:
                    sx = 0.0
                    rx = 0
                    for file in fx.keys():
                        sx += 100.0
                        if fx[file] < self.stat['MinGO'] and os.path.exists(file):
                            os.unlink(file)
                            rx += 1
                        self.log.printLog('\r#MIN','Screening %s Datasets <%d seqs; %s deleted; %.1f%%' % (spec,self.stat['MinGO'],rje.integerString(rx),sx/len(fx)),newline=False,log=False)
                    self.log.printLog('\r#MIN','Screening %s Datasets <%d seqs: %s deleted.' % (spec,self.stat['MinGO'],rje.integerString(rx)))

        except: self.log.errorLog('Wah wah wah! Cry like a baby - ensLociGo() is busted.',quitchoice=True)
#########################################################################################################################
    ### <8> ### EnsDAT custom UniProt file generator - see also unifake.py
#########################################################################################################################
    def ensDat(self): ### Runs accessory prediction programs to make fake UniProt DAT file.
        '''Runs accessory prediction programs to make fake UniProt DAT file.'''
        try:### ~ [1] Setup - check file paths etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['EnsDat']: return      # No species for EnsDat processing
            for path in ['HMMerPath','TMHMM','SignalP','PFam']:
                if not os.path.exists(self.info[path]):
                    self.log.errorLog('%s path "%s" does not exist. Will not use.' % (path,self.info[path]),printerror=False)
                    self.info[path]= 'none'

            ### ~ [2] Perform EnsDat processing and generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['EnsDat']:
                self.info['Name'] = spec
                self.makeEnsDat(self.speciesCode())
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def makeEnsDat(self,spec):  ### Generates an EnsDat file for the given species code
        '''Generates an EnsDat file for the given species code.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup list/file of previously run sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datpickup = []
            if os.path.exists(self.info['DatPickup']):
                if not self.opt['Force']: datpickup = rje.listFromCommand(self.info['DatPickup'])
                else: os.unlink(self.info['DatPickup'])
            ## ~ [1b] Setup SeqList for ens_loci sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            enseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % self.ensLociFas(spcode=spec,read=True),'autoload=T'])
            seqdict = enseq.seqNameDic()
            (sx,seqnum) = (0,enseq.seqNum())
            ## ~ [1c] Setup UniProt object and output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)   # UniProt object for saving data
            datfile = '%sens_%s.ensdat.dat' % (self.info['EnsPath'],spec)
            if not datpickup and os.path.exists(datfile): rje.backup(self,datfile)
            ## ~ [1d] Setup RJE_HMM object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hmm = rje_hmm.HMMRun(self.log,self.cmd_list)
            hmmfile = '%sens_%s.pfam.tdt' % (self.info['EnsPath'],spec)
            if not datpickup and os.path.exists(hmmfile): rje.backup(self,hmmfile)
            hmm.list['HMM'] = [self.info['PFam']]
            hmm.opt['HMMPFam'] = True
            hmm.setInfo({'SearchDB':'tmp.fas','HMMOut':'tmp.hmm.out'})      # This will be made for each sequence
            ## ~ [1e] Setup RJE_TM object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tm = rje_tm.TM(self.log,self.cmd_list)

            ### ~ [2] Perform EnsDat processing and generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for name in rje.sortKeys(seqdict):      ## Process and output in alphabetical order! ##
                sx += 1
                seq = seqdict[name]
                if seq.shortName() in datpickup:
                    self.log.printLog('#ENSDAT','Skipping %s: already present in %s' % (name,self.info['DatPickup']),log=False)
                    continue   # Already done!
                self.log.printLog('#SEQ','Processing %s (%s aa) %s...' % (name,rje.integerString(seq.aaLen()),seq.info['Description'][:50]))
                try:
                    open('tmp.fas','w').write('>%s\n%s\n' % (seq.shortName(),seq.info['Sequence']))
                    udata = {'OS':[self.getSpecies(spec)],'CC':['-!- Features generated using rje_ensembl.py']}
                    #!# Consider adding more accession numbers at this point #!#
                    ft = []     # List of features for sequence

                    ## ~ [2a] IUPRED disorder prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try:
                        seq.disorder()
                        for disorder in seq.obj['Disorder'].list['RegionDisorder']:
                            ft.append({'Type':'DISORDER','Desc':'IUPred predicted disorder region','Start':disorder[0],'End':disorder[1]})
                        for fold in seq.obj['Disorder'].list['RegionFold']:
                            ft.append({'Type':'ORDER','Desc':'IUPred predicted ordered region','Start':fold[0],'End':fold[1]})
                    except: self.log.errorLog('EnsDat disorder problem for %s.' % name)

                    ## ~ [2b] PFam HMM domain prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try:
                        hmm.search = []
                        #x#hmm.readHMMPFamSearch(resfile=hmm.hmmSearch(self.info['PFam'],outfile='tmp.hmm.out'),readaln=False)
                        hmm.list['HMMRes'] = [hmm.hmmSearch(self.info['PFam'],outfile='tmp.hmm.out')]   # Used in hmmTable
                        #hmm.readHMMSearch(resfile=hmm.hmmSearch(self.info['PFam'],outfile='tmp.hmm.out'),readaln=False)
                        hmm.hmmTable(outfile=hmmfile,append=True)
                        disorder = seq.obj['Disorder'].list['ResidueDisorder'] # individual (IUPRed) residue results
                        iucut = seq.obj['Disorder'].stat['IUCut']
                        for search in hmm.search:
                            for hit in search.hit:
                                for aln in hit.aln:
                                    region = disorder[aln.stat['SbjStart']-1:aln.stat['SbjEnd']]
                                    hmmdisorder = float(sum(region)) / len(region)
                                    ft.append({'Start':aln.stat['SbjStart'],'End':aln.stat['SbjEnd'],
                                               'Desc':'%s PFam HMM Eval: %.2e; Score: %.1f; IUPRed: %.2f' % (search.info['Name'],aln.stat['Expect'],aln.stat['BitScore'],hmmdisorder)})
                                    if hmmdisorder >= iucut: ft[-1]['Type'] = 'PFAM'
                                    else: ft[-1]['Type'] = 'DOMAIN'
                    except: self.log.errorLog('EnsDat PFam HMM problem for %s.' % name)
                    
                    ## ~ [2c] TMHMM transmembrane topology prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try:
                        tmdat = os.popen('%s tmp.fas -short' % self.info['TMHMM']).readlines()
                        domlist = rje_tm.domainList(rje_tm.parseTMHMM(tmdat[0]))
                        for tmdom in domlist:
                            ft.append(tmdom)
                            ft[-1]['Desc'] = 'TMHMM topology prediction'
                            ft[-1]['Start'] = rje.atoi(ft[-1]['Start'])
                            ft[-1]['End'] = rje.atoi(ft[-1]['End'])
                        if len(domlist) > 1: udata['CC'].append('TMHMM: %d TM domains; N-Term %s' % ((len(domlist)-1)/2,domlist[0]['Type']))
                        else: udata['CC'].append('TMHMM: 0 TM domains')
                    except: self.log.errorLog('EnsDat TMHMM problem for %s.' % name)
                    
                    ## ~ [2d] SIGNALP signal peptide prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try:
                        os.system('%s -f short -t euk tmp.fas > tmp.signalp' % self.info['SignalP'])
                        tm.signalp = {}
                        tm.parseSignalP('tmp.signalp')
                        sigp = tm.signalp.pop(seq.shortName())
                        cpos = 0
                        if sigp['nn_ymax?'] == 'Y':
                            cpos = rje.atoi(sigp['nn_ymaxpos'])
                            desc = 'SignalP NN prediction'
                        if sigp['hmm_cmax?'] == 'Y':
                            hmm_c = rje.atoi(sigp['hmm_cmaxpos'])
                            if cpos == 0:
                                cpos = hmm_c
                                desc = 'SignalP HMM prediction'
                            else:
                                if hmm_c < cpos:
                                    cpos = hmm_c
                                    desc = 'SignalP HMM prediction (NN also Y)'
                                else: desc += ' (HMM also Y)'
                        if cpos > 0: ft.append({'Type':'SIGNALP','Desc':desc,'Start':1,'End':cpos})
                    except: self.log.errorLog('EnsDat SignalP problem for %s.' % name)

                    ## ~ [2e] Convert to UniProt and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #self.deBug(ft)
                    uniprot.list['Entry'] = []
                    if uniprot.addFromSeq(seq,data=udata,ft=ft):    ### Converts into UniProtEntry object 
                        uniprot.saveUniProt(datfile,append=True)
                        open(self.info['DatPickup'],'a').write('%s\n' % seq.shortName())

                ## ~ [2f] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                except: self.log.errorLog('Problem during EnsDat(%s)' % name)
                for tmp in glob.glob('TMHMM*'):
                    if os.path.isdir(tmp): os.rmdir(tmp)
                for tmp in glob.glob('tmp*'): os.unlink(tmp)
                self.log.printLog('#ENSDAT','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(sx),rje.integerString(seqnum-sx)),log=False)

        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION III: EnsEMBL Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: EnsForker Class                                                                                         #
#########################################################################################################################
class EnsForker(rje_forker.Forker):
    '''
    EnsForker Class. Author: Rich Edwards (2014).

    Str:str
    - ForkDir = Alternative directory for forking output. [`./`]

    Bool:boolean
    - LogFork = Whether to log forking (if quick, might simply be a counter in main script)
    - PIDCheck = Whether to generate a `*.pid` file to check PID status.
    - RjePy = Whether the jobs being run are rje_python jobs and therefore need log files.

    Int:integer
    - IOLimit = Limit of number of IOErrors before termination [50]

    Num:float
    - ForkSleep = Sleep time (seconds) between cycles of forking out more process [0]
    - KillTime = Time to monitor for killing of hanging forks [3600]
    - MemFree = Minimum % memory that should be free to start new fork. (Not yet implemented.)

    List:list
    - Forked = List of dictionaries containing data regarding forked processes.
    - ResFile = List of results files (BASEFILE.X) that will need transferring []
    - ToFork = List of data to be given to runFork() method. (runFork() usually replaced by inheriting class)

    Dict:dictionary

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    #i# See rje_forker.Forker
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def startFork(self,fdict):  ### Sets a new fork going using the data in fdict.
        '''Sets a new fork going using the data in fdict.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdict['species'] = self.list['ToFork'].pop(0)
            fdict['ID'] = 'Fork %d' % self.list['Forked'].index(fdict)
            fdict['FID'] = 'f_%s' % rje.randomString(6)
            fdict['Log'] = '%s%s.log' % (self.getStr('RunPath'),fdict['FID'])
            fdict['cmd'] = 'basefile=%s' % fdict['FID']
            fdict['ResFile'] = self.list['ResFile'][0:]
            try: open(fdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting fork.'); return self.endJob(fdict)
            ### ~ [2] ~ Add Fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setNum({'KillTime':time.time()})
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                fdict['PID'] = cpid
                if self.getBool('Win32'): self.printLog('#FORK','Forking %s as %s: %d remain; No memory check (win32=T).' % (fdict['species'],cpid,len(self.list['ToFork'])))
                elif self.getBool('OSX'): self.printLog('#FORK','Forking %s as %s: %d remain; No memory check (osx=T).' % (fdict['species'],cpid,len(self.list['ToFork'])))
                else: self.printLog('#FORK','Forking %s as %s: %d remain; %.1f%% mem free' % (fdict['species'],cpid,len(self.list['ToFork']),fdict['Mem']))
            else:                   # child process
                newlog = rje.setLog(info=makeInfo(),out=self,cmd_list=[fdict['cmd']])
                newens = EnsEMBL(newlog,self.cmd_list)
                newens.dict['SiteSpec'] = self.dict['SiteSpec']
                newens.obj['DB'] = self.obj['DB']
                newens.ensLoci(fdict['species'])
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('EnsForker.startFork error')
#########################################################################################################################
### End of SECTION IV: EnsForker Class                                                                                  #
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
    try: EnsEMBL(mainlog,cmd_list).run()
        
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

