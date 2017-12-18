#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       REVERT
Description:  Retrovirus and Endogenous Viral Element Reconstruction Tool
Version:      0.9.1
Last Edit:    08/09/16
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    REVERT (the Retrovirus and Endogenous Viral Element Reconstruction Tool) is an automated utility for discovery
    candidate viruses and endogenous viral element (EVE) sequences in host genomes. (Transcriptome data, such as EST
    libraries, could also be used as input and REVERT should be able to retrieve expressed contemporary viruses too.)

    NOTE: Version 0.7.0 has reworked the pipeline slightly. Updated docs coming soon. (REST service output might also be
    temporarily affected.)

REVERT Analysis Pipeline:

    1. A set of viral genomes is (downloaded to and) loaded from a genbank file, as determined by `virusgb=FILE/LIST`.
       This is then used to create a `*.full.fas` file of full-length genomes (for later mapping) and `*.prot.fas` file
       of proteins for initial searching against a fasta-format host genome, set by `genome=FILE`. The output basefile
       defaults to `virusgb.genome` but can be over-ruled by `basefile=FILE`. NB. A set of files can be set up for a
       batch run using `vbatch=FILELIST` and `gbatch=FILELIST` (wildcards allowed).

    2. GABLAM tblastn search of viral proteins against host genome DNA. This will produce `*.gablam.tdt`, `*.hitsum.tdt`
       and `*.local.tdt` summary files, and a fasta file per protein in a directory set by `fasdir=PATH`. The fixed
       settings for GABLAM are `fasout=T fragfas=T combined=T`. Other GABLAM settings can be altered on the commandline.
       By default, GABLAM will use the BLAST complexity filter and composition-based statistics (`blastf=T blastcf=T`).
       This GABLAM search uses the fullblast=T switch by default. (Use fullblast=F to switch off.) Only Ordered GABLAM
       output is returned (outstats=GABLAMO) with a default cut-off of 50% viral protein (gablamcut=0.5
       cutstat=OrderedAlnLen cutfocus=Query).

    3. The combined hits (if any) from the GABLAM search are mapped back onto the viral proteins and genome using a
       second BLAST search without the complexity filter: (1) a blastn search, aligning genome hits against the viral
       genome; (2) a blastx search aligning the viral proteins against the genome hits; (3) a tblastn search, aligning
       the translated genome hits against the original viral proteins. The latter is converted into a set of fasta
       alignments. The assemblies use the BLAST+ flat query-anchored output and can include overlapping hits.
       NOTE: The second BLAST may identify additional protein hits between the viral proteins and the genome regions that
       were missed by the original BLAST! Such examples might have HitNum of 0 but FragNum > 0.

    4. The protein alignments are summarised for each virus-genome comparison and converted into % coverage and %
       identity (of the viral proteins), which are summed up across all proteins in the viral proteome. If running in
       batch mode, these will be compiled into a single protein table and single virus summary table:
       - *.revert.details.tdt = each viral protein against each genome
       - *.revert.tdt = each virus against each genome

    5. Compiled results for the batch file run will then be generated and used to make additional outputs. For the virus
       and genome tables, an extra Alias field will be added to use for visualisation outputs. If `aliasfile=FILE` is
       given then aliases will be read from this file. Where missing, commandline prompts for aliases will be given. If
       no alias is given (or `aliases=F`) then the full virus or genome name will be used. Outputs are:
       - *.revert.virus.tdt = a single line for each virus, summarising the total number of genomes hit.
       - *.revert.protein.tdt = a single line for each viral protein, summarising the total number of genomes hit.
       - *.revert.genome.tdt = a single line for each genome, summarising the total number of viruses hit.

    6. TO BE ADDED:
       A final compilation mode will GABLAM all (non-redundant) genome fragments against hit proteins, establish the best
       hits and apply stricter cut-offs for both protein and virus hits.
       mincov is also applied to the initial GABLAM.


Commandline:
    ### ~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    virusgb=FILE/LIST   # Either a genbank download (must end *.gb) or a list of GenBank viral nuccore UIDs [virus.gb]
    vbatch=FILELIST     # Run REVERT on a list of virusgb files (over-rides virusgb=FILE). Wildcards allowed. []
    genome=FILE         # Fasta file containing contigs or chromosomes of host genome [genome.fas]
    gbatch=FILELIST     # Run REVERT on a list of genome files (over-rides genome=FILE). Wildcards allowed. []
    searchdb=FILE       # Fasta file of viral (& host?) proteins for assembly annotation. Use viral proteins if None. []
    taxdir=PATH/        # Will look in this directory for taxonomy input files if not found ['SourceData/']
    virusdir=PATH/      # Will look in this directory for single viral gb files ['VirusGB/']
    sourcedate=DATE     # Source file date (YYYY-MM-DD) for taxonomy files to preferentially use [None]

    ### ~ Basic Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=FILE       # Root of output files. Defaults to `virusgb.genome` if blank/None [None]
    revertdir=PATH/     # Directory to output run files into. (Prefixes basefile for single runs.) [REVERT/]
    fasdir=PATH         # Directory in which to save fasta files [RevertDir/GABLAMFAS/genome/]
    blastdir=PATH       # BLAST directory for GABLAM BLAST files if keepblast=T [RevertDir/BLAST/]
    keepblast=T/F       # Whether to keep GABLAM BLAST files for reference/reuse [True]
    blaste=X            # BLAST evalue cut-off [1e-10]

    ### ~ Results Compilation/Filtering Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    revertnr=T/F        # Compile batch run results with non-redundancy and quality filtering [True]
    gablamfrag=X        # Length of gaps between mapped residues for fragmenting local hits [100]
    minpcov=X           # Min. %coverage for viral proteins during NR reduction [40.0]
    minplocid=X         # Min. local %identity for viral proteins during NR reduction [30.0]
    minvcov=X           # Min. %coverage for viral genomes during NR reduction [40.0]
    minvlocid=X         # Min. local %identity for viral genomes during NR reduction [30.0]

    ### ~ Advanced Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    aliases=T/F         # Whether to use aliases in additional outputs [True]
    aliasfile=FILE      # Delimited file of 'Name','Alias' to use in place of batch summary files [*.alias.tdt]
    vspcode=T/F         # Whether to use viral species codes (even if invented) for VirusID aliases [True]
    gspcode=T/F         # Whether to use genome species codes for Genome field if able to parse name [False]
    vgablam=T/F         # Whether to compile viral genomes/proteomes and conduct all-by-all GABLAM for graphs [True]
    graphformats=LIST   # Formats for virus-genome graph outputs (svg/xgmml/png/html) (dev=T only) [xgmml,png]
    treeformats=LIST    # List of output formats for generated trees (see rje_tree.py) (dev=T only) [nwk,png]

    ### ~ Adavanced Run Options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    farmgablam=X        # Whether to run a pre-REVERT farming of batch BLAST searches using X forks each [0]

    ### ~ Special Run Options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    vgbparse=T/F        # Whether to parse virus IDs from vbatch tables and output BASEFILE.HOST.acc files [False]
    vhost=LIST          # List of viral hosts to output new files for. (Blank for all) []

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_blast_V2, rje_db, rje_forker, rje_genbank, rje_menu, rje_obj, rje_ppi, rje_seqlist, rje_taxonomy, rje_tree
import rje_dismatrix_V2 as rje_dismatrix
import gablam, fiesta
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Improved pickup of existing results.
    # 0.2.0 - Additional output for visualisation of results.
    # 0.3.0 - Added extra mode for generating viral accnum lists from http://www.ncbi.nlm.nih.gov/genome/viruses/.
    # 0.4.0 - Altered *.protein.tdt output to *.details.tdt output and added new *.protein.tdt file.
    # 0.4.1 - Modified to use the new FullBlast GABLAM mode for speed.
    # 0.4.2 - Modified to farm out GABLAM BLAST searches as forks.
    # 0.4.3 - Temp fix for QSub runs with large wall-times being used for PPI Spring Layout walltime.
    # 0.4.4 - Removed limitation of farmgablam being < 1/2 forks.
    # 0.5.0 - Altered defaults and added some extra default GABLAM filtering.
    # 0.6.0 - Addition of revertnr=T/F method to remove redundancy across searches. Made default.
    # 0.7.0 - Reworking and tidying to make use of virus directory the default.
    # 0.7.1 - Minor bug fixing for REVERT REST Server.
    # 0.7.2 - Minor bug fix for updated RJE_DB module.
    # 0.8.0 - Added Local Table output to REST server and fixed virus parse bug.
    # 0.8.1 - Fixed virus genbank from ID file bug. (Need to re-check REST server.)
    # 0.9.0 - Major reworking of Alias generation/use. Deletion of some OLD methods.
    # 0.9.1 - Fixed up a few bugs and outputs. Ready for REST testing.
    # 0.9.2 - Added gablamfrag=100 default to GABLAM object call to hopefully counter GABLAM default switch to 1.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Add addition of addflanks=X to GABLAM fragment hits.
    # [Y] : Improve file organisation. Add revertdir=PATH/ for run directory.
    # [Y] : Rename the final protein assembly files. Do not need such a long filename!
    # [?] : Add some filters including min length of match and extra complexity filter for repetitive sequences.
    # [ ] : Allow SeqIn over-ride of VirusGB for analysis without viral genome mapping. (Or different proteins.)
    # [ ] : Improve REVERT forking: currently forking within GABLAM.
    # [ ] : Fix odd bug when possible incomplete BLAST causes a GABLAM forking runtime error when forks=1.
    # [Y] : Add additional visualisation output.
    # [ ] : Document additional visualisation output.
    # [Y] : Add EdgeTypes for XGMML output.
    # [ ] : Add optional pre-screening of vector sequence database.
    # [ ] : Add viral all-by-all protein and genome GABLAM searches and incorpotate in XGMML output.
    # [ ] : Add option to remove certain text from genomes (".dna.toplevel")
    # [ ] : Add organism to *.revert.tdt output? (Not just genome name.)
    # [ ] : Add zipping of blast results. (And other results?)
    # [ ] : Add compilation=T/F method.
    # [ ] : Add actual consensus sequence to consensus.tdt - put method in rje_seqlist.
    # [ ] : Improved handling (warning not error) for viruses without protein sequence data.
    # [ ] : Add ORF prediction for unannotated sequences (and SeqIn override)? May have overlaps issue. (May already?!)
    # [ ] : Look for genome file and abort earlier if missing.
    # [ ] : Add combined local hit tables to REST/main output.
    # [Y] : Why does the *.alias.tdt file constantly want backing up? After "MAKING ALIASES".
    # [ ] : Add warning if mutliple genomes get the same alias.
    # [ ] : Add genome fragment fasta file to main/REST output.
    # [ ] : Should we consider making a DNA-based mode, rather than protein-based? (Will people criticise the protein aspect)
    # [ ] : Reorganise and tidy.
    # [ ] : Update Docs.
    # [ ] : Check gablamfrag, fragmerge and addflanks in the light of GABLAM updates.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('REVERT', '0.9.2', 'February 2017', '2014')
    description = 'Retrovirus and Endogenous Viral Element Reconstruction Tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
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
### SECTION II: REVERT Class                                                                                            #
#########################################################################################################################
class REVERT(rje_obj.RJE_Object):
    '''
    REVERT Class. Author: Rich Edwards (2014).

    Str:str
    - AliasFile=FILE      # Delimited file of 'Name','Alias' to use in place of batch summary files [None]
    - FasDir=PATH       # Directory in which to save fasta files [None]
    - Genome=FILE       # Fasta file containing contigs or chromosomes of host genome [genome.fas]
    - RevertBlastDir=PATH/   # Directory to output initial GABLAM BLAST files into. [REVERT/BLAST]
    - RevertDir=PATH/   # Directory to output run files into. (Prefixes basefile for single runs.) [REVERT/]
    - SearchDB=FILE     # Fasta file of viral (& host?) proteins for assembly annotation. Use viral proteins if None. []
    - SeqIn=FILE        # Viral protein sequence set [VirusGB.prot.fas]
    - VirusDir=PATH/      # Will look in this directory for single viral gb files ['VirusGB/']
    - VirusGB=FILE/LIST # Either a genbank download (must end *.gb) or a list of GenBank viral nuccore UIDs [virus.gb]

    Bool:boolean
    - Aliases=T/F         # Whether to use aliases in additional outputs [True]
    - AliasBackup=T/F     # Whether to ask about backing up alias file [True]
    - GSpCode=T/F         # Whether to use genome species codes for aliases if able to parse name [True]
    - RevertNR=T/F        # Compile batch run results with non-redundancy and quality filtering [True]
    - VGABLAM=T/F         # Whether to compile viral genomes/proteomes and conduct all-by-all GABLAM for graphs [True]
    - VGbParse=T/F        # Whether to parse virus IDs from vbatch tables and output new *.acc files [False]
    - VSpCode=T/F         # Whether to use viral species codes (even if invented) for aliases [True]

    Int:integer
    - FarmGABLAM=X         # Whether to run a pre-REVERT farming of batch BLAST searches using X forks each [0]
    - GablamFrag=X        # Length of gaps between mapped residues for fragmenting local hits [100]

    Num:float
    - MinPCov=X            # Min. %coverage for viral proteins during NR reduction [40.0]
    - MinPLocID=X          # Min. local %identity for viral proteins during NR reduction [30.0]
    - MinVCov=X            # Min. %coverage for viral genomes during NR reduction [40.0]
    - MinVLocID=X          # Min. local %identity for viral genomes during NR reduction [30.0]

    List:list
    - GBatch=FILELIST     # Run REVERT on a list of genome files (over-rides genome=FILE). Wildcards allowed. []
    - GraphFormats=LIST   # Formats for virus-genome graph outputs (svg/xgmml/png/html) [xgmml,png]
    - TreeFormats=LIST    # List of output formats for generated trees (see rje_tree.py) [nwk,png]
    - VBatch=FILELIST     # Run REVERT on a list of virusgb files (over-rides virusgb=FILE). Wildcards allowed. []
    - VHost=LIST          # List of viral hosts to output new files for. (Blank for all) []

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database Object
    - FIESTA = fiesta.FIESTA Object
    - GABLAM = gablam.GABLAM Object
    - Taxonomy = rje_taxonomy.Taxonomy Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['AliasFile','FasDir','Genome','GBatch','VBatch','RevertBlastDir','RevertDir','SearchDB','SeqIn','VirusDir','VirusGB']
        self.boollist = ['Aliases','AliasBackup','GSpCode','RevertNR','VSpCode','VGABLAM']
        self.intlist = ['FarmGABLAM','GablamFrag']
        self.numlist = ['MinPCov','MinPLocID','MinVCov','MinVLocID']
        self.listlist = ['GBatch','VBatch','VHost','GraphFormats','TreeFormats']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'FasDir':'','RevertDir':rje.makePath('REVERT/'),'VirusGB':'virus.gb','Genome':'genome.fas',
                     'RevertBlastDir':rje.makePath('REVERT/BLAST/'),'VirusDir':rje.makePath('VirusGB/')})
        self.setBool({'Aliases':True,'AliasBackup':True,'GSpCode':False,'RevertNR':True,'VSpCode':True,'VGbParse':False,'VGABLAM':True})
        self.setInt({'GablamFrag':100})
        self.setNum({'MinPCov':40.0,'MinPLocID':30.0,'MinVCov':40.0,'MinVLocID':30.0})
        self.list['GraphFormats'] = ['xgmml','png']
        self.list['TreeFormats'] = ['nwk','png']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.cmd_list = ['fullblast=T','tuplekeys=T'] + self.cmd_list
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['Taxonomy'] = rje_taxonomy.Taxonomy(self.log,self.cmd_list)
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
                self._cmdRead(cmd,type='str',att='RevertBlastDir',arg='blastdir')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['GBatch','VBatch'])   # Normal strings
                self._cmdReadList(cmd,'path',['FasDir','RevertBlastDir','RevertDir','VirusDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['AliasFile','Genome','SearchDB','VirusGB'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['Aliases','AliasBackup','GSpCode','RevertNR','VGbParse','VSpCode','VGABLAM'])
                self._cmdReadList(cmd,'int',['FarmGABLAM','GablamFrag'])   # Integers
                self._cmdReadList(cmd,'perc',['MinPCov','MinPLocID','MinVCov','MinVLocID']) # Percentages
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['GraphFormats','TreeFormats','VHost'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['GBatch','VBatch']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        for lkey in ['GraphFormats','TreeFormats']:
            self.list[lkey] = string.split(string.join(self.list[lkey]).lower())
            if len(self.list[lkey]) == 1 and self.list[lkey] == 'none': self.list[lkey] = []
        #self.setNum({'MinCov':self.getNum('MinVCov'),'MinLocID':self.getNum('MinVLocID')})
        if not self.getStr('RevertDir')[:1] == '/':
            self.setStr({'RevertDir':rje.makePath('%s%s' % (self.getStr('RunPath'),self.getStr('RevertDir')))})
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        `REVERT` is a high throughput pipeline for finding evidence of Endogenous viral elements (EVEs) within a genome.
        Input is one or more genbank IDs corresponding to annotated viral genomes (e.g. NC_001542), and a host genome.
        Available genomes can be accessed at http://rest.slimsuite.unsw.edu.au/alias&genome. Please contact us to add
        more genomes. (Alternatively, the standalone version of REVERT can be run on multiple bespoke genomes.)

        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT:

        nr = Main REVERT summary table of non-redundant (NR) hits between viruses and the search genome. [tdt]
        nr.details = Breakdown by protein of main REVERT results. [tdt]
        consensus = reconstructed consensus sequences from genome fragments for each NR protein [fas]
        revert = Full REVERT summary table before NR/QC filtering [tdt]
        details = Breakdown of full summary table by protein [tdt]
        local = Details of local alignments positions in genome and gfrag sequences [tdt]
        gfrag = Extracted genomic regions with putative EVEs [fas]
        virusgb = Genbank format download of input viruses [gbk]
        virus.locus = Summary table of viruses extracted from genbank [tdt]
        virus.feature = Table of annotated viral features. (Protein/CDS needed for REVERT analysis) [tdt]
        virus.prot = Fasta file of viral proteins used in analysis [fas]
        virus.full = Full genome sequences of viruses used in analysis [fas]
        genome = Genome sequence used for analysis. (May be too big to download) [fas]

        More details of outputs can be found in the [REVERT Help](http://rest.slimsuite.unsw.edu.au/revert) and/or
        [SLiMSuite blog](http://slimsuite.blogspot.com.au/). Please get in touch if anything is not clear.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            self.dict['Output']['virusgb'] = 'VirusGB'
            self.dict['Output']['genome'] = 'Genome'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        return ['nr','nr.details','consensus','revert','details','local','gfrag','virusgb',
                'virus.locus','virus.feature','virus.prot','virus.full','genome']
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,batch=True):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('VGbParse'): return self.vgbParse()
            #if self.list['GBatch'] or self.list['VBatch']: return self.batchRun()
            if batch: return self.batchRun()    # Always use this for file consistency.
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.gablamSearch(): return False
            #if not self.virusAssembly('host'): return False
            #if not self.virusAssembly(): return False
            ## ~ [2a] ~ Convert protein assemblies into alignments and statistics ~~~~~~~~~~~~~~~~~ ##
            #if not self.virusAssembly('protein'): return False
            if not self.revertSummary(): return False
            return True
        except SystemExit: raise
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    def batchRun(self): ### Performs a batch run against multiple input files.
        '''Performs a batch run against multiple input files.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['VBatch']:
                if self.getStrLC('VBatch'): raise IOError('Cannot find vbatch files!')
                if not self.getStrLC('VirusGB'): raise IOError('Cannot find virusgb file!')
                self.list['VBatch'] = [self.getStr('VirusGB')]
            if not self.list['GBatch']:
                if self.getStrLC('GBatch'): raise IOError('Cannot find gbatch files!')
                if not self.getStrLC('Genome'): raise IOError('Cannot find genome file!')
                self.list['GBatch'] = [self.getStr('Genome')]
            self.restSetup()
            self.dict['Output']['genome'] = 'Genome file: %s' % self.getStr('Genome')
            ## ~ [1a] ~ Setup DB object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.log.no_suppression += ['entry_overwrite','Invented SpCode']
            self.log.warnings += ['SpCode Missing TaxID']
            db = self.obj['DB']
            if not self.getStrLC('Basefile'): self.setBaseFile('revert')
            self.setBaseFile(self.basefile(runpath=True))
            #!# Note that the VirusGB field has been removed. This might be useful to track for output.
            #!# Might be better to just deal with this mapping entirely separately.
            #!# The real question is how it impacts the NR reduction of hits.
            # revert.details table, which contains the individual protein vs genome data
            rhead = ['VirusAcc','Virus','Genome','GFile','Protein','ProtAcc','Product','HitNum','MaxScore','EVal','Length','FragNum','Coverage','Identity','Local']
            rkeys = ['VirusAcc','GFile','Protein']
            rdb = db.addTable(name='revert.details',mainkeys=rkeys,expect=False)    # Load if present
            if not rdb: rdb = db.addEmptyTable('revert.details',rhead,rkeys)        # Create new table if absent
            # revert table, which contains the compiled virus vs genome data
            vhead = ['VirusAcc','Virus','Genome','GFile','HitNum','Length','ProtNum','ProtHits','ProtAcc','FragNum','Coverage','Identity','Local']
            vkeys = ['VirusAcc','GFile']
            vdb = db.addTable(name='revert',mainkeys=vkeys,expect=False)            # Load if present
            if not vdb: vdb = db.addEmptyTable('revert',vhead,vkeys)                # Create new table if absent
            ## ~ [1b] ~ Setup viral genomes/proteomes for searching ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupViruses()
            self.list['VBatch'] = self.db('data.virus').indexKeys('locus')  #!# Should just be able to use keys?
            #i# The virus batch list is now individual viruses.
            ## ~ [1c] ~ Setup viral/genome aliases? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# This is currently performed AFTER the batch run. (Not required for main processing.)

            ### ~ [2] ~ Run REVERT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] FarmGABLAM pre-BLAST mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getInt('FarmGABLAM') > 0: self.farmGABLAM()
            ## ~ [2b] Main REVERT processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Now run individual virus vs genome comparisons and compile into main tables
            rcmd = self.cmd_list + ['vbatch=','gbatch=']
            rx = 0; rtot = len(self.list['VBatch']) * len(self.list['GBatch']); fx = 0; sx = 0
            self.printLog('#VBATCH','%d VBatch viral genomes.' % len(self.list['VBatch']))
            self.printLog('#GBATCH','%d GBatch host genomes.' % len(self.list['GBatch']))
            self.printLog('#BATCH','Total batch runs: %s (%d x %d).' % (rje.iStr(rtot),len(self.list['VBatch']),len(self.list['GBatch'])))
            for virus in self.list['VBatch'][0:]:
                run_genomes = vdb.indexDataList('VirusAcc',virus,'Genome')
                for gfile in self.list['GBatch']:
                    gfa = rje.baseFile(gfile,strip_path=True) #+ '.fa'
                    if gfa in run_genomes and not self.force():
                        self.printLog('#SKIP','Skipping batch run %d of %d: virusgb="%s"; genome="%s".' % (rx,rtot,virus,gfile))
                        rx += 1; sx += 1
                        continue
                    if rcmd[-1].startswith('basefile='): rcmd = rcmd[:-1]
                    self.printLog('#~~#','### ==================================================================== ###',timeout=False); rx += 1
                    self.printLog('#BATCH','Batch running %d of %d: virusgb="%s"; genome="%s".' % (rx,rtot,virus,gfile))
                    rcmd.append('basefile=%s.%s' % (virus,rje.baseFile(gfile,strip_path=True)))
                    revert = REVERT(self.log,rcmd)
                    revert.db().list['Tables'] = [self.db('data.virus'),self.db('VirusFT')]
                    revert.setStr({'VirusGB':virus,'Genome':gfile,'BatchBase':self.basefile()})
                    revert.list['VBatch'] = []
                    revert.list['GBatch'] = []
                    if revert.run(batch=False) and revert.db('revert'):
                        db.mergeTables(vdb,revert.db('revert'))
                        db.mergeTables(rdb,revert.db('revert.details'))
                    else: fx += 1
                    self.printLog('#BATCH','Batch run %d of %d complete; %s skipped; %d failed.' % (rx,rtot,sx,fx))
            self.printLog('#~~#','### ==================================================================== ###',timeout=False)
            self.printLog('#BATCH','%s batch runs complete; %s skipped; %d failed.' % (rx,sx,fx))
            if rx == fx: self.errorLog('No successful runs!', printerror=False); return False
            ## ~ [2c] Main REVERT output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rdb.saveToFile(sfdict={'Coverage':3,'Identity':3,'Local':3})    # Detailed table of protein vs genome hits (*.revert.details.tdt)
            vdb.saveToFile(sfdict={'Coverage':3,'Identity':3,'Local':3})    # Pairwise virus vs genome table (*.revert.tdt)
            self.bugPrint('NC_006577: %s' % vdb.indexEntries('VirusAcc','NC_006577'))

            ### ~ [3] ~ RevertNR QC and NR filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('RevertNR'): return True
            fformat = {'HitNum':'int','MaxScore':'num','EVal':'num','Length':'int','FragNum':'int',
                       'Coverage':'num','Identity':'num','Local':'num','ProtNum':'int','ProtHits':'int'}
            rdb.dataFormat(fformat)
            vdb.dataFormat(fformat)
            self.bugPrint('NC_006577: %s' % vdb.indexEntries('VirusAcc','NC_006577'))
            ## ~ [3a] - Sort out aliases prior to NR etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add better processing of aliases here prior to NR etc. Then only need to deal with it once! #!#
            #!# Need to document when/where the aliases are actually used and what stays as unique keys
            #i# VirusAcc for viruses and GFile for genomes are unique
            #i# Virus and Genome are user-friendly aliases.
            if not self.makeAliases(): return False

            ## ~ [3b] - Split results tables by genome for processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            grdict = db.splitTable(rdb,'GFile',asdict=True,keepfield=True)
            db.deleteTable(rdb)
            gvdict = db.splitTable(vdb,'GFile',asdict=True,keepfield=True)
            #db.deleteTable(vdb) - used for XGMML output
            ## ~ [3c] - Coverage/Identity QC Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            append = False
            for genome in rje.sortKeys(grdict):
                grdb = grdict[genome]
                gvdb = gvdict[genome]
                self.bugPrint('NC_006577: %s' % gvdb.indexEntries('VirusAcc','NC_006577'))
                #!# Should move this QC to a method! #!#
                # QC Filtering of viruses and proteins based on coverage and local identity. This is now done early to
                # make subsequent steps more efficient. Filtering is done based on genome-specifc viral hits - by genome
                # First, Viral hits failing QC and their proteins are removed.
                self.bugPrint('NC_006577: %s' % self.num)
                badvirus = []   # List of viruses failing QC
                for ventry in gvdb.entries():
                    #if ventry['VirusAcc'] == 'NC_006577': self.debug(ventry)
                    if ventry['Coverage'] < self.getNum('MinVCov') or ventry['Local'] < self.getNum('MinVLocID'):
                        badvirus.append(ventry['VirusAcc'])
                gvdb.dropEntriesDirect('VirusAcc',badvirus)
                grdb.dropEntriesDirect('VirusAcc',badvirus)
                self.printLog('#QC','%s %s viruses removed (low coverage and/or identity): %s remain' % (rje.iLen(badvirus),genome,rje.iStr(gvdb.entryNum())))
                if not gvdb.entryNum(): continue
                # Next, proteins failing the protein QC threshold are removed. Their coverage is reduced to 0.0 and any
                # viruses now failing the viral QC are removed.
                badvirus = []   # List of viruses with proteins failing QC - will need full stats to be recalculated
                for pentry in grdb.entries():
                    #if pentry['HitNum']: self.debug(pentry)
                    if pentry['Coverage'] < self.getNum('MinPCov') or pentry['Local'] < self.getNum('MinPLocID'):
                        pentry['Coverage'] = pentry['Identity'] = pentry['Local'] = 0.0
                        pentry['FragNum'] = pentry['HitNum'] = 0
                        if pentry['VirusAcc'] not in badvirus: badvirus.append(pentry['VirusAcc'])
                self.printLog('#QC','%s %s viruses had low coverage and/or identity proteins removed' % (rje.iLen(badvirus),genome))
                badvirus.sort()
                for virus in badvirus[0:]:
                    #self.debug(virus)
                    vlen = vcov = vloc = vid = 0
                    for pentry in grdb.indexEntries('VirusAcc',virus):
                        vlen += pentry['Length']
                        if min(pentry['Length'],pentry['Coverage']) <= 0: continue
                        vcov += (pentry['Coverage'] * pentry['Length'] / 100.0)
                        #vloc += (pentry['Local'] * vcov / 100.0)
                        vid += (pentry['Identity'] * pentry['Length'] / 100.0)
                    #self.debug('%s: Len=%s; Cov=%s; Loc=%s; ID=%s.' % (virus,vlen,vcov,vloc,vid))
                    if not vcov or not vlen: continue
                    if 100.0 * vcov / vlen >= self.getNum('MinVCov') or 100.0 * vid / vcov >= self.getNum('MinVLocID'):
                        badvirus.remove(virus)
                gvdb.dropEntriesDirect('VirusAcc',badvirus)
                grdb.dropEntriesDirect('VirusAcc',badvirus)
                self.printLog('#QC','%s %s viruses removed (low coverage and/or identity proteins): %s remain' % (rje.iLen(badvirus),genome,rje.iStr(gvdb.entryNum())))
                if not grdb.entryNum(): continue
                #i# revertNR need only deal with redundancy removal
                #i# It should use the remaining viruses and proteins (with coverage) in rdb to filter local database
                #i# Only count decent hits towards coverage.
                if not self.revertNR(genome,grdb,append): return False
                append = True

            ## ~ [3d] ~ Compile local and gfrag files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            locfile = '%s.revert.local.tdt' % self.baseFile()
            self.dict['Output']['local'] = locfile
            rje.backup(self,locfile,appendable=False)
            locdb = None
            gfragfile = '%s.revertnr.genomefrag.fas' % self.baseFile()
            rje.backup(self,gfragfile,appendable=False)
            lx = fx = 0
            for gfile in rje.sortKeys(grdict):
                gbase = gfile
                basefile = os.path.basename(self.baseFile())
                nrdir = rje.makePath('%sNR/%s/' % (self.getStr('RevertDir'),basefile))
                nrbase = '%s%s' % (nrdir,gbase) #os.path.basename(self.baseFile()))
                nrloc = '%s.local.tdt' % nrbase
                if rje.exists(nrloc):
                    nrlocdb = db.addTable(nrloc,name='nrlocdb',mainkeys=['Qry','Hit','AlnNum'],expect=True)
                    if locdb: db.mergeTables(locdb,nrlocdb)
                    else: locdb = nrlocdb
                    #if rje.exists(locfile): open(locfile,'a').writelines(open(nrloc,'r').readlines()[1:])
                    #else: open(locfile,'w').write(open(nrloc,'r').read())
                    lx += 1
                nrfrag = '%s.gfrag.fas' % nrbase
                if rje.exists(nrfrag):
                    open(gfragfile,'w').write(open(nrfrag,'r').read())
                    fx += 1
            # Process
            if locdb:
                locdb.makeField(formula='#Qry#',fieldname='GFrag',after='Qry')
                locdb.renameField('Qry','GenSeq')
                locdb.renameField('Hit','VirProt')
                locdb.renameField('QryStart','FragStart')
                locdb.renameField('QryEnd','FragEnd')
                locdb.renameField('SbjStart','ProtStart')
                locdb.renameField('SbjEnd','ProtEnd')
                locdb.addFields(['GenStart','GenEnd'])
                for entry in locdb.entries():
                    ldata = string.split(entry['GenSeq'],'.')
                    entry['GenSeq'] = string.join(ldata,'.')[:-1]
                    ldata = string.split(ldata[-1],'-')
                    entry['FragStart'] = int(entry['FragStart'])
                    entry['FragEnd'] = int(entry['FragEnd'])
                    if entry['FragStart'] > entry['FragEnd']:
                        entry['GenStart'] = int(ldata[0]) + entry['FragEnd'] - 1
                        entry['GenEnd'] = int(ldata[0]) + entry['FragStart'] - 1
                    else:
                        entry['GenStart'] = int(ldata[0]) + entry['FragStart'] - 1
                        entry['GenEnd'] = int(ldata[0]) + entry['FragEnd'] - 1
                locdb.saveToFile(locfile)
            # Report
            self.printLog('#LOCAL','Local NR BLAST hits for %s of %s genomes output to %s' % (lx,len(grdict),locfile))
            self.printLog('#GFRAG','Putative EVE Genome Fragments for %s of %s genomes output to %s' % (fx,len(grdict),gfragfile))
            if fx: self.dict['Output']['gfrag'] = gfragfile

            for table in ['nr','nr.details','consensus','revert','details']:
                if os.path.exists('%s.%s.tdt' % (self.basefile(),table)): self.dict['Output'][table] = '%s.%s.tdt' % (self.basefile(),table)
                elif os.path.exists('%s.revert%s.tdt' % (self.basefile(),table)): self.dict['Output'][table] = '%s.revert%s.tdt' % (self.basefile(),table)
                elif os.path.exists('%s.revert.%s.tdt' % (self.basefile(),table)): self.dict['Output'][table] = '%s.revert.%s.tdt' % (self.basefile(),table)
                elif table in self.dict['Output']: self.dict['Output'].pop(table)

            ### ~ [4] ~ Additional Batch Summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #?# Can this be made into a method? #?#
            #i# Review of all aliases to control compression is handled above in makeAliases()

            #if len(self.list['GBatch']) > 1:
            #    self.printLog('#NR','RevertNR currently unavailable for multiple genomes.')
            #    return False
            #if not self.revertNR(): return False

            #?# Add VBatch files (VirusGB) to summary files from self.db('data.virus')
            vlocdb = self.db('data.virus')
            rdb = self.db('revertnr.details',add=True,mainkeys=['VirusAcc','GFile','Protein'])
            vdb = self.db('revertnr',add=True,mainkeys=['VirusAcc','GFile'])
            self.debug(vdb.data())
            for table in [rdb,vdb]:
                table.addField('VirusGB','VirusID')
                for entry in table.entries():
                    vacc = entry['VirusAcc']
                    entry['VirusGB'] = vlocdb.data(vacc)['vfiles']
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~ ADDITIONAL BATCH SUMMARIES ~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            ## ~ [3b] ~ Virus Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vdb.dataFormat({'HitNum':'int','Length':'int','FragNum':'int','Coverage':'num','Identity':'num','ProtNum':'int','ProtHits':'int'})
            vsumdb = db.copyTable(vdb,'summary.virus')     # ['VirusGB','VirusAcc','Virus','Genome','HitNum','Length','ProtNum','ProtHits','Protein','FragNum','Coverage','Identity','Local']
            try: vsumdb.renameField('ProtAcc','Protein')
            except: pass
            for entry in vsumdb.entries(): entry['GFile'] = {True:1,False:0}[entry['FragNum']>0]
            vsumdb.compress(['VirusID'],rules={'Length':'max','Genome':'list','VirusGB':'list','VirusAcc':'list','Protein':'list','Identity':'max','Coverage':'max','Local':'max'},default='sum',joinchar='|')
            for ekey in vsumdb.dataKeys():
                entry = vsumdb.data(ekey)
                entry['Protein'] = rje.sortUnique(string.split(entry['Protein'],'|'))
                if '' in entry['Protein']: entry['Protein'].remove('')
                entry['ProtHits'] = len(entry['Protein'])
                entry['Protein'] = string.join(entry['Protein'],'|')
            vsumdb.dropField('Protein')
            ## ~ [3c] ~ Genome Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gsumdb = db.copyTable(vdb,'summary.genome')     # ['VirusAcc','Virus','Genome','HitNum','Length','ProtNum','ProtHits','Protein','FragNum','Coverage','Identity','Local']
            self.debug(gsumdb.fields())
            try: gsumdb.renameField('ProtAcc','Protein')
            except: pass
            gsumdb.index('VirusGB',splitchar='|')
            vgbfields = gsumdb.indexKeys('VirusGB')
            for vfile in vgbfields: gsumdb.addField(vfile,evalue=0)
            for entry in gsumdb.entries(): entry['Virus'] = {True:1,False:0}[entry['FragNum']>0]
            for vfile in gsumdb.indexKeys('VirusGB'):
                for entry in gsumdb.indexEntries('VirusGB',vfile): entry[vfile] = entry['VirusAcc']
            gsumdb.compress(['Genome','VirusAcc'],rules={'Length':'max','Protein':'list','VirusID':'list','Virus':'max'},default='sum')
            gsumdb.compress(['Genome'],rules={'Length':'max','Protein':'list','VirusID':'list','Virus':'sum'},default='max')
            gsumdb.setFields(['Genome','GFile','Virus','VirusID','HitNum','MaxScore','EVal','Protein','ProtHits','FragNum','BestProt','Best']+vgbfields)
            for ekey in gsumdb.dataKeys():
                entry = gsumdb.data(ekey)
                entry['Protein'] = rje.sortUnique(string.split(entry['Protein'],'|'))
                if '' in entry['Protein']: entry['Protein'].remove('')
                entry['ProtHits'] = len(entry['Protein'])
                entry['Protein'] = string.join(entry['Protein'],'|')
            gsumdb.dropField('Protein')
            ## ~ [3d] ~ Protein Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            psumdb = db.copyTable(rdb,'summary.protein')     # ['VirusAcc','Virus','Genome','Protein','ProtAcc','Product','HitNum','MaxScore','EVal','Length','FragNum','Coverage','Identity','Local']
            # Output fields: ['VirusAcc','Virus','Genome','ProtAcc','Product','HitNum','Length','FragNum','Coverage','Identity','Local']
            for entry in psumdb.entries(): entry['GFile'] = {True:1,False:0}[entry['FragNum']>0]
            psumdb.compress(['VirusAcc','ProtAcc'],rules={'Length':'max','VirusAcc':'list','VirusID':'list','Genome':'list','Identity':'max','Coverage':'max','Local':'max'},default='sum',joinchar='|')
            psumdb.setFields(['VirusAcc','Virus','VirusID','ProtAcc','Product','GFile','Genome','HitNum','Length','FragNum','Coverage','Identity','Local'])

            ## ~ [3e] ~ Save tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vsumdb.saveToFile()
            gsumdb.saveToFile()
            psumdb.saveToFile()

            ### ~ [4] ~ Additional Batch Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ### >>> DEV ONLY AT PRESENT >>> ###
            #x#if not self.dev(): return True
            self.batchOutput()
            return True
        except SystemExit: raise
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    def setupViruses(self): ### Sets up viral genomes and proteomes for searching.
        '''Main class setup method.'''
        try:### ~ [0] General Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            virusdir = rje.makePath('%s/' % os.path.abspath(self.getStr('VirusDir')))    # Directory for viral data
            rje.mkDir(self,self.getStr('VirusDir'),log=True)
            db = self.db()
            ldb = None
            fdb = None
            #!# Should we have a virusGB table of file:accnum?
            if not self.setupAlias(): return False
            ### ~ [1] Check Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~~ VIRUS GENBANK SETUP ~~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            savelocus = True    # Whether to save a combined table of viral loci
            for vfile in self.list['VBatch'][0:]:   # Each file might be a file or IUD list.
                ## ~ [1a] Deal with accnum list if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                virusgb = None
                vfetch = vfile
                if not os.path.exists(vfile):  # Interpret as IUD list
                    savelocus = len(self.list['VBatch']) > 1    # Cancel save locus for a single UID list
                    if ',' in vfile:   # List given: need a new name
                        virusgb = vfile
                        vfile = '%svirus.input' % self.getStr('RunPath')   # This is only really for the REST server
                    else: vfile = '%s%s' % (virusdir,vfile)     # This is only really for the REST server
                if not virusgb: virusgb = rje.baseFile(vfile,strip_path=True)
                ## ~ [1b] Create file if required and parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# This will create a file for the full set of viruses. Separate files to be created later.
                #i# This file will be in the source directory of the viral gb files
                #gbcmd = self.cmd_list + ['fasout=full,gene,prot','tabout=T','basefile=%s%s' % (virusdir,virusgb)]
                gbcmd = self.cmd_list + ['fasout=full,gene,prot','tabout=T','basefile=%s' % rje.baseFile(vfile)]
                if vfile.lower().endswith('.gb') or vfile.lower().endswith('.gbk'): gbcmd += ['seqin=%s' % vfile]
                else: gbcmd += ['fetchuid=%s' % vfetch]
                gbcmd += ['locusout=T','locusdir=%s' % self.getStr('VirusDir'),'addtags=T']
                #self.debug(gbcmd)
                gb = rje_genbank.GenBank(self.log,gbcmd)
                gb.obj['Taxonomy'] = self.obj['Taxonomy']
                #!# Should check for need to run #!#
                rungb = self.force() or self.getStrLC('Rest')
                if not rungb:
                    # Look for gb.db('Locus')
                    locdb = gb.db('Locus',add=True,mainkeys=['locus'])
                    # Look to compile gb.db('Feature') using Locus db and VirusDir
                    if locdb:
                        ftdb = None
                        for locus in locdb.dataKeys():
                            ftfile = '%s%s.Feature.tdt' % (self.getStr('VirusDir'),locus)
                            if not rje.exists(ftfile): rungb = True; break
                            if ftdb:
                                ftdb.obj['File'] = open(ftfile,'r')
                                ftdb.obj['File'].readline()
                                ftdb.readSet(['locus'],clear=False)
                                if ftdb.obj['File']: ftdb.obj['File'].close()
                            else:
                                ftdb = gb.db().addTable(ftfile,mainkeys=['locus','feature','position'],name='Feature',expect=True)
                    else: rungb = True
                if rungb: gb.run()
                if vfile.lower().endswith('.gb') or vfile.lower().endswith('.gbk'): self.dict['Output']['virusgb'] = vfile
                else: self.dict['Output']['virusgb'] = '%s.gb' % rje.baseFile(vfile)
                #i# This should now have generated virus-specific output in virusdir
                ## ~ [1c] Sort out REST output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.dict['Output']['virus.locus'] = '%s.Locus.tdt' % rje.baseFile(vfile)
                self.dict['Output']['virus.feature'] = '%s.Feature.tdt' % rje.baseFile(vfile)
                self.dict['Output']['virus.prot'] = '%s.prot.fas' % rje.baseFile(vfile)
                self.dict['Output']['virus.full'] = '%s.full.fas' % rje.baseFile(vfile)
                if not rje.exists(self.dict['Output']['virus.feature']) and self.getStrLC('Rest'):
                    gb.db('Feature').saveToFile()
                    gb.saveFasta(logskip=self.v()>0,bylocus=False)
                ## ~ [1d] Grab parsed genbank feature data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ldb:
                    #X#db.mergeTables(ldb,gb.db('Locus'),overwrite=False)
                    vbase = rje.baseFile(vfile,strip_path=True)
                    for entry in gb.db('Locus').entries():
                        ekey = ldb.makeKey(entry)
                        if ekey in ldb.data(): ldb.data(ekey)['vfiles'] += '|%s' % vbase
                        else: entry['vfiles'] = vbase; ldb.addEntry(entry)
                else:
                    ldb = gb.db('Locus'); db.list['Tables'].append(ldb); ldb.basefile(db.basefile())
                    ldb.addField('vfiles',evalue=rje.baseFile(vfile,strip_path=True))
                if fdb:
                    if not db.mergeTables(fdb,gb.db('Feature'),overwrite=False): raise ValueError
                    self.printLog('#MERGE','%s features merged -> %s ViralFT' % (rje.iStr(gb.db('Feature').entryNum()),rje.iStr(fdb.entryNum())))
                    db.debug(fdb.data())
                else: fdb = gb.db('Feature'); db.list['Tables'].append(fdb); fdb.basefile(db.basefile())

            #!# NB. Might now need to compile prot.fas data for REST server.

            ### ~ [2] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #virus table:
            #locus   length  type    definition      accession       version gi      organism        spcode vfiles
            #NC_001434       7176    ss-RNA  Hepatitis E virus, complete genome.     NC_001434       NC_001434.1     9626440 Hepatitis E virus       HEV
            for entry in ldb.entries():
                if entry['vfiles'].startswith('|'): entry['vfiles'] = entry['vfiles'][1:]
            ldb.setStr({'Name':'data.virus'}); ldb.index('locus')
            if savelocus: ldb.saveToFile()
            #locus	feature	position	start	end	locus_tag	protein_id	details
            #NC_001434	CDS	4..5085	4	5085	HEVgp01	NP_056779.1	/note="ORF 1" /codon_start="1" /product="polyprotein" /db_xref="GI:9626448" /db_xref="GeneID:1494415"
            fdb.setStr({'Name':'VirusFT'})
            fdb.dropEntriesDirect('feature',['CDS','Protein'],inverse=True,log=True,force=True)    ### Drops certain entries from Table
            #self.debug(fdb.index('feature'))
            fdb.newKey(['locus','protein_id'])
            fdb.index('protein_id'); fdb.index('locus')

            ### ~ [3] Check/create sequence input files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vtot = len(ldb.index('locus')); vx = 0
            for virus in ldb.indexKeys('locus'):
                seqfile = '%s%s.prot.acc.fas' % (virusdir,virus)
                if self.force() or not rje.exists(seqfile):
                    vprotfile = '%s%s.prot.fas' % (virusdir,virus)
                    if rje.exists(vprotfile):
                        seqlist = rje_seqlist.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=%s' % vprotfile,'autoload=T'])
                        seqlist.saveSeq(reformat='accfas',seqfile=seqfile)
                    else: self.warnLog('Protein file %s not found. (No CDS features?)' % vprotfile)
                if rje.exists(seqfile): vx += 1
                else:
                    self.warnLog('Failed to create %s: removing virus "%s" from analysis' % (seqfile,virus))
                    ldb.dropEntriesDirect('locus',[virus])
                    fdb.dropEntriesDirect('locus',[virus])
            self.printLog('#SEQ','Successful file creation/discovery for %s of %s viruses.' % (rje.iStr(vx),rje.iStr(vtot)))
            return True     # Setup successful
        except: self.errorLog('REVERT setupViruses failed.'); return False  # Setup failed
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Check Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            genome = rje.baseFile(self.getStr('Genome'),strip_path=True)
            virusgb = self.getStr('VirusGB')
            self.basefile('%s.%s' % (virusgb,genome))
            rje.mkDir(self,self.getStr('RevertDir'))     # Make RevertDir - has runPath included
            if not os.path.exists(self.getStr('RevertDir')): raise IOError('Failed to find/make RevertDir: %s' % self.getStr('RevertDir'))
            self.basefile('%s%s' % (self.getStr('RevertDir'),self.basefile()))  # Add RevertDir path
            db.basefile('%s%s' % (self.getStr('RevertDir'),self.basefile()))  # Add RevertDir path
            self.cmd_list += ['basefile=%s' % self.basefile()]  # Ensure that subsequent objects have same basefile
            if not self.getStrLC('FasDir'): self.setStr({'FasDir':rje.makePath('%sGABLAMFAS/%s/' % (self.getStr('RevertDir'),genome))})
            ## ~ [1a] Check for previous run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rpfile = '%s.revert.details.tdt' % self.basefile()
            rsfile = '%s.revert.tdt' % self.basefile()
            if rje.exists(rpfile) and rje.exists(rsfile) and not self.force():
                self.printLog('#REVERT','%s REVERT summary tables found.' % self.basefile())
                rkeys = ['VirusAcc','GFile','Protein']
                db.addTable(name='revert.details',mainkeys=rkeys,expect=True)
                db.addTable(name='revert',mainkeys=rkeys[:2],expect=True)
                self.printLog('#SKIP','Skipping %s (force=F).' % self.basefile())
                return False    # Should be no need to rerun anything as should also be in summary files
            return True     # Setup successful
        except: self.errorLog('REVERT setup failed.'); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Main REVERT Search Methods                                                                              #
#########################################################################################################################
    def farmGABLAM(self): # Farms out virus x genome GABLAM searches for vbatch x gbatch                          #V0.7.0
        '''Farms out virus x genome GABLAM searches for vbatch x gbatch.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vbatch = self.list['VBatch']; gbatch = self.list['GBatch']
            #i# vbatch is now a list of viruses with seqin files: '%s%s.prot.acc.fas' % (virusdir,virus)
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~ FARM GABLAM SEARCH ~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            batchx = len(vbatch) * len(gbatch)
            farmforks = self.getInt('Forks') / self.getInt('FarmGABLAM')
            while farmforks > batchx and self.getInt('Forks') > self.getInt('FarmGABLAM'):
                self.setInt({'FarmGABLAM':self.getInt('FarmGABLAM')+1})
                farmforks = self.getInt('Forks') / self.getInt('FarmGABLAM')
            if farmforks == 1: self.setInt({'FarmGABLAM':self.getInt('Forks')})
            self.printLog('#FARM','Farming off %d x %d fork GABLAM runs (%d total forks)' % (farmforks,self.getInt('FarmGABLAM'),self.getInt('Forks')))
            forker = rje_forker.Forker(self.log,self.cmd_list+['forks=%d' % farmforks,'rjepy=F'])
            ### ~ [1] Make the Fork job list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = self.cmd_list + ['vbatch=','gbatch=']
            for virus in vbatch:
                for gfile in gbatch:
                    revert = REVERT(self.log,rcmd)
                    revert.setStr({'VirusGB':virus,'Genome':gfile})
                    gabcmd = revert.gablamSearch(runcmd=True)
                    if not gabcmd: continue
                    gabcmd = 'python %stools/gablam.py %s forks=%d fullblast=T' % (rje.slimsuitepath,string.join(gabcmd),self.getInt('FarmGABLAM'))
                    forker.list['ToFork'].append(gabcmd)
            ### ~ [2] Fork out! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forker.run()
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    def gablamSearch(self,runcmd=False): ### Step 1: Search viral proteins against genome using tblastn.
        '''Step 1: Search viral proteins against genome using tblastn.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gablamdir = rje.makePath('%sGABLAM/' % self.getStr('RevertDir'))
            rje.mkDir(self,gablamdir,log=True)
            virus = self.getStr('VirusGB'); vfile = '%s%s.prot.fas' % (self.getStr('VirusDir'),virus)
            gfile = self.getStr('Genome'); genome = rje.baseFile(gfile,strip_path=True)
            gbase = '%s%s.%s' % (gablamdir,virus,genome)
            if not runcmd: self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~ GABLAM SEARCH ~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            ## ~ [1a] Check for existing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rungablam = self.force()
            for gabfile in ['gablam','hitsum','local']: rungablam = rungablam or not rje.exists('%s.%s.tdt' % (gbase,gabfile))
            if not rungablam:
                self.printLog('#SKIP','%s GABLAM run files found (force=F).' % gbase)
                return not runcmd
            ## ~ [1b] Setup GABLAM run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # GablamCut is not being used. Should be taken care of by the QC filter later
            bcmd = ['blaste=1e-10','keepblast=T','fullblast=T','outstats=GABLAMO','gablamfrag=100']
            gcmd = ['fasout=F','fragfas=F','combinedfas=F',     # Obsolete with qassemble=T - Time and Space guzzler.
                    'seqin=%s' % vfile,'searchdb=%s' % gfile,'blastp=tblastn',
                    'local=T','dismat=F','blastf=T','blastcf=T','qryacc=T','blastdir=%s' % self.getStr('RevertBlastDir'),
                    'qassemble=T','basefile=%s' % gbase]
            gcmd = rje.tidyArgs(bcmd + self.cmd_list + gcmd)
            ### ~ [2] ~ Return command or run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if runcmd: return gcmd  #bcmd + self.cmd_list + gcmd
            self.obj['GABLAM'] = gablam.GABLAM(self.log,gcmd)   #bcmd+self.cmd_list + gcmd)
            return self.obj['GABLAM'].gablam()
            #?# Could be made faster by:
            #?# Modify to not perform GABLAM calculation! (Maybe bypass GABLAM and just use BLAST?) #!#
            #?# Add sequence extraction to rje_seqlist (don't read sequence into memory until extraction point reached)
            #?# Use blast local hit table and rje_seqlist for sequence extraction directly
            #?# Modify GABLAM to use this code instead of its own.
            #?# Modify GABLAM calculation to truncate lists and max possible start and ends for pair, given local hit table
            # ... but the slow point seems to be reading the fragmented sequences using rje_seqlist.
        except: self.errorLog('REVERT.gablamSearch() error.'); return False
#########################################################################################################################
    def virusAssembly(self,seqtype='genome',basefile=None,vbase=None): ### Step 2: Search viral genome against GABLAM hits and assemble.
        '''Step 2: Search viral genome against GABLAM hits and assemble.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug('>>> ASSEMBLY?')
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~ VIRUS ASSEMBLY ~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            if not basefile: basefile = self.basefile()
            if not vbase: vbase = rje.baseFile(self.getStr('VirusGB'))
            gfile = '%s.fas' % basefile      # Output sequence fragments
            if seqtype == 'genome':
                vfile = '%s%s.full.fas' % (self.getStr('VirusDir'),vbase)
                if not rje.exists(vfile): raise IOError('Viral full genome (%s) missing!' % vfile)
                blastp = 'blastn'
            elif seqtype == 'protein':
                vfile = '%s.prot.fas' % vbase
                if not rje.exists(vfile): raise IOError('Viral protein file (%s) missing!' % vfile)
                blastp = 'tblastn'
            elif seqtype == 'host':
                gfile = '%s.prot.fas' % vbase
                if not rje.exists(gfile): raise IOError('Viral protein file (%s) missing!' % gfile)
                vfile = '%s.fas' % basefile
                if not rje.exists(vfile): self.printLog('#SKIP','No GABLAM hits found (%s): no host assembly.' % vfile); return True
                blastp = 'blastx'
            else: raise ValueError('Sequence type "%s" not recognised.' % seqtype)
            #bfile = rje.makePath('%s.ASSEMBLY/' % self.basefile()) + '%s.%s.assembly' % (os.path.basename(self.basefile()),seqtype)
            bfile = rje.makePath('%sASSEMBLY/' % self.getStr('RevertDir')) + '%s.%s.assembly' % (os.path.basename(basefile),seqtype)
            if rje.exists(bfile) and not self.force():
                self.printLog('#BLAST','BLAST assembly file "%s" found (force=F).' % bfile)
                return True
            if not rje.exists(gfile): self.printLog('#SKIP','No GABLAM hits found (%s): no %s assembly.' % (gfile,seqtype)); return True
            self.printLog('#ALIGN','%s assembly of %s vs %s' % (blastp,vfile,gfile))
            bcmd = ['blastopt=-outfmt "4"','blastp=%s' % blastp,'formatdb=T','blastf=F','blastcf=T',
                    'blasti=%s' % vfile,     # Input file (BLAST -i FILE) [None]
                    'blastd=%s' % gfile,     # BLAST database (BLAST -d FILE) [None]
                    'blasto=%s' % bfile]   # Output file (BLAST -o FILE) [*.blast]
            # Read from GABLAM?    blastv=X        : Number of one-line hits per query (BLAST -v X) [500]
            # Read from GABLAM?    blastb=X        : Number of hit alignments per query (BLAST -b X) [250]
            blast = rje_blast_V2.blastObj(self.log,['blaste=1e-10']+self.cmd_list+bcmd)
            blast.log.warnings.append("comp_based_statswarn")
            ### ~ [2] ~ Run BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.blast(wait=True,cleandb=True,use_existing=True,log=True)
            ### ~ [3] ~ Convert BLAST Assembly into fasta alignment and statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'%sASSEMBLYFAS/' % self.getStr('RevertDir'))
            if seqtype in ['protein']:  #,'host']:    #?# Add host too? (Instead?) Needs to deal with RF: will break at present.
                blast.assemblyToAlignments(bfile,fullquery=True)
                fx = 0
                for pafile in glob.glob('%s.*' % bfile):
                    afasfile = '%sASSEMBLYFAS/%s' % (self.getStr('RevertDir'),os.path.basename(string.replace(pafile,'.%s.assembly.' % seqtype,'.')))
                    os.rename(pafile,afasfile); fx += 1
                if fx: self.printLog('#MOVE','%s fasta files moved from ASSEMBLY/ to ASSEMBLYFAS/' % rje.iStr(fx))
            #!# Impose a minimum % identity filter: Should GABLAM be re-run for the viral protein versus hit fragments
            #.. and/or their translations?
            #!# Need more than a simple GABLAM of initial hits, as will want to combine by species for summary stats
            #!# Probably want to add an over-arching summary results file to batch mode that single mode can append
            return rje.exists(bfile)
        except: self.errorLog('REVERT.virusAssembly() error.'); return False
#########################################################################################################################
    def revertSummary(self):    ### Combine the different results into a summary of REVERT success.
        '''Combine the different results into a summary of REVERT success.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug('>>> SUMMARY?')
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~ REVERT SUMMARY ~~~~~~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            db = self.db()
            rje.mkDir(self,rje.makePath('%sPAIRWISE/' % self.getStr('RevertDir')))
            rpfile = rje.makePath('%sPAIRWISE/' % self.getStr('RevertDir')) + '%s.revert.details.tdt' % self.basefile(strip_path=True)
            rsfile = rje.makePath('%sPAIRWISE/' % self.getStr('RevertDir')) + '%s.revert.tdt' % self.basefile(strip_path=True)
            if rje.exists(rpfile) and rje.exists(rsfile) and not self.force():
                self.printLog('#REVERT','%s REVERT summary tables found.' % self.basefile())
                rkeys = ['VirusAcc','GFile','Protein']
                db.addTable(rpfile,name='revert.details',mainkeys=rkeys,expect=True)
                db.addTable(rsfile,name='revert',mainkeys=rkeys[:2],expect=True)
                return True
            virus = self.getStr('VirusGB')
            virusdir = self.getStr('VirusDir')
            ## ~ [0a] ~ Virus tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #locus	feature	position	start	end	locus_tag	protein_id	details
            #NC_001434	CDS	4..5085	4	5085	HEVgp01	NP_056779.1	/note="ORF 1" /codon_start="1" /product="polyprotein" /db_xref="GI:9626448" /db_xref="GeneID:1494415"
            vftdb = self.db('VirusFT')
            #locus   length  type    definition      accession       version gi      organism        spcode
            #NC_001434       7176    ss-RNA  Hepatitis E virus, complete genome.     NC_001434       NC_001434.1     9626440 Hepatitis E virus       HEV
            vlocdb = self.db('data.virus')
            ## ~ [0b] ~ GABLAM HitSum table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #Qry     HitNum  MaxScore        EVal
            #NP_056779.1     4       3455.0  0.0
            self.printLog('#BASE', self.baseFile())
            db.basefile(self.basefile())
            gbase = rje.makePath('%sGABLAM/' % self.getStr('RevertDir')) + rje.stripPath(self.basefile())
            ghfile = '%s.hitsum.tdt' % gbase
            ghdb = db.addTable(ghfile,['Qry'],name='HitSum',expect=True)
            glfile = '%s.local.tdt' % gbase
            gldb = db.addTable(glfile,['Qry','Hit','AlnNum'],name='local',expect=True)
            ## ~ [0c] ~ Protein Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            vprotfile = '%s%s.prot.fas' % (virusdir,virus)
            vseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % vprotfile,'seqmode=tuple','autoload=T','autofilter=F'])
            vseqdict = vseqlist.makeSeqNameDic(keytype='accnum',clear=True,warnings=True)
            ## ~ [0d] ~ Output table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rhead = ['VirusAcc','Virus','Genome','GFile','Protein','ProtAcc','Product','HitNum','MaxScore','EVal','Length','FragNum','Coverage','Identity','Local']
            rkeys = ['VirusAcc','GFile','Protein']
            rdb = db.addEmptyTable('revert.details',rhead,rkeys)    #!# Load and return if present and not self.force()
            #adb = db.addEmptyTable('revert.consensus',['Protein','Product','Genome','Consensus'],['Protein','Product','Genome'])    #!# Load and return if present and not self.force()
            ### ~ [1] ~ Generate REVERT Protein table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ventry = vlocdb.data(virus)
            if not vftdb.indexDataList('locus',virus,'protein_id'):
                self.warnLog('No proteins found in Feature Table for %s' % virus)
            for prot in vftdb.indexDataList('locus',virus,'protein_id'):
                if not prot: self.warnLog('Locus "%s" is returning a blank protein_id' % virus,'blank_protein_id'); continue
                if prot not in vseqdict: self.warnLog('Locus "%s" is missing from Viral sequence dictionary' % virus,'vseq_missing'); continue
                #fentry = vftdb.data('%s\t%s' % (virus,prot))
                fentry = vftdb.data(vftdb.makeKey({'locus':virus,'protein_id':prot}))
                fdata = string.split(' %s' % fentry['details'],' /')
                while fdata and not fdata[0].startswith('product'): fdata.pop(0)
                if fdata:
                    try: product = rje.matchExp('product="(.+)"',fdata[0])[0]
                    except: self.debug(fdata)
                else: product = 'unknown'
                (pname,pseq) = vseqlist.getSeq(vseqdict[prot])
                pname = string.split(pname)[0]
                ## ~ [1a] Global Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                plen = len(pseq)
                fragnum = len(gldb.indexEntries('Qry',prot))
                coverage = 0.0
                identity = 0.0
                gentry = ghdb.data(prot)
                if gentry:
                    #coverage = 100.0 * float(gentry['Coverage']) / float(gentry['Length'])
                    #identity = 100.0 * float(gentry['Identity']) / float(gentry['Length'])
                    # Now convert to % later so that genome summing is weighted for length
                    coverage = float(gentry['Coverage'])
                    identity = float(gentry['Identity'])
                else:
                    self.warnLog('Protein %s missing from %s!' % (prot,ghfile),quitchoice=True)
                    gentry = {'HitNum':-1,'MaxScore':-1,'EVal':-1}
                paln = rje.makePath('%sASSEMBLYFAS/%s.%s.%s.fas' % (self.getStr('RevertDir'),virus,rje.baseFile(self.getStr('Genome'),strip_path=True),pname),wholepath=True)
                if os.path.exists(paln):
                    aseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % paln,'seqmode=tuple','autoload=T','autofilter=F'])
                    qseq = aseqlist.seqs()[0][1]
                    consensus = []  #'.'] * len(qseq)
                    for i in range(len(qseq)):
                        consensus.append('.')
                        if qseq[i] == '-': consensus[-1] = '-'; continue
                        for hit in aseqlist.seqs()[1:]:
                            hseq = hit[1]
                            if hseq[i] == '-': continue
                            if hseq[i] == qseq[i]:
                                consensus[-1] = qseq[i]
                                #identity += 100.0
                                break
                            consensus[-1] = 'x'     #?# Should an actual consensus be produced too? Harder!
                        #if consensus[-1] not in '.-': coverage += 100.0
                    consensus = string.join(consensus,'')
                    self.verbose(0,1,'%s %s\n%s' % (pname,product,consensus))
                    #adb.addEntry({'Genome':rje.baseFile(self.getStr('Genome'),strip_path=True),'Protein':pname,'Product':product,'Consensus':consensus})
                rentry = {'VirusAcc':virus,'Virus':ventry['organism'],'Genome':rje.baseFile(self.getStr('Genome'),strip_path=True),
                          'Protein':pname,'ProtAcc':prot,'Product':product,'GFile':rje.baseFile(self.getStr('Genome'),strip_path=True),
                          'HitNum':int(gentry['HitNum']),'MaxScore':gentry['MaxScore'],'EVal':gentry['EVal'],
                          'FragNum':fragnum,'Length':plen,'Coverage':coverage,'Identity':identity,'Local':0.0}
                #if coverage > 0: rentry['Local'] = 100.0 * identity / coverage
                #if rentry['HitNum'] > 0: self.debug(rentry)
                rdb.addEntry(rentry)
            # Copy this table for compressing into the summary table with coverage and identity as raw values (not %)
            vsumdb = db.copyTable(rdb,'revert')
            #for rentry in rdb.entries():
            #        rentry['Coverage'] = 100.0 * rentry['Coverage'] / rentry['Length']
            #        rentry['Identity'] = 100.0 * rentry['Identity'] / rentry['Length']

            ### ~ [2] ~ Genome Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vsumdb.addField('ProtHits',evalue=0)
            vsumdb.addField('ProtNum',evalue=1)
            for entry in vsumdb.entries():
                if entry['HitNum'] > 0: entry['ProtHits'] = 1
                else: entry['Protein'] = ''
            vsumdb.compress(['VirusAcc','GFile'],rules={'ProtAcc':'list','Protein':'list','GFile':'str'},joinchar='|',default='sum')    #
            for entry in vsumdb.entries():
                if len(string.split(entry['ProtAcc'],'|')) != entry['ProtNum']: self.debug('ProtAcc/number mismatch! %s vs %s' % (entry['ProtNum'],entry['ProtAcc']))
                if (entry['Protein'] or entry['ProtHits']) and len(string.split(entry['Protein'],'|')) != entry['ProtHits']: self.debug('Protein/hits mismatch! %s vs %s' % (entry['ProtHits'],entry['Protein']))
            vsumdb.keepFields(['VirusAcc','Virus','Genome','GFile','HitNum','Length','ProtNum','ProtHits','ProtAcc','FragNum','Coverage','Identity','Local'])
            vsumdb.list['Fields'] = ['VirusAcc','Virus','Genome','GFile','HitNum','Length','ProtNum','ProtHits','ProtAcc','FragNum','Coverage','Identity','Local']

            ### ~ [3] Save tables to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Convert to percentage values
            for entry in vsumdb.entries() + rdb.entries():
                if entry['Coverage']: entry['Local'] = 100.0 * entry['Identity']/entry['Coverage']
                if entry['Length']:
                    entry['Identity'] = 100.0 * entry['Identity']/entry['Length']
                    entry['Coverage'] = 100.0 * entry['Coverage']/entry['Length']
            rdb.saveToFile(rpfile)
            vsumdb.saveToFile(rsfile)
            #if adb.entries(): adb.saveToFile()
            return True
        except: self.errorLog('Problem with REVERT.revertSummary()'); return False
#########################################################################################################################
    def batchOutput(self):  ### Additional all-by-all comparisons and graphics
        '''Additional all-by-all comparisons and graphics.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~ ALL-BY-ALL BATCH OUTPUT ~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            db = self.db()
            vdb = self.db('revertnr',add=False)   # Table of virus-genome summary data
            rdb = self.db('revert',add=False)
            vlocdb = self.db('data.virus',add=False)   # Table of viruses
            #psumdb = self.db('summary.protein',add=False)
            vsumdb = self.db('summary.virus',add=False)   # Table of viruses
            gsumdb = self.db('summary.genome',add=False)  # Table of genomes
            vgabdb = None               # Table of viral genome all-by-all gablam (if run).

            ## ~ [1] ~ All-by-all GABLAMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('VGABLAM'):
                ## Viral Genomes ##
                vgabdir = rje.makePath('%sVGABLAM/' % self.getStr('RevertDir'))
                rje.mkDir(self,vgabdir)
                vbase = '%s%s.vgenome' % (vgabdir,self.basefile(strip_path=True))
                vfullfile = '%s.fas' % vbase
                if self.force() or not rje.exists(vfullfile):
                    vfullseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=list'])
                    for vfile in self.list['VBatch']:   # Was replaced with virus locus
                        vfull = '%s%s.full.fas' % (self.getStr('VirusDir'),vfile)
                        vfullseq.loadSeq(vfull,nodup=True,clearseq=False,mode='tuple')
                    vfullseq.saveSeq(seqfile=vfullfile,reformat=None,append=False,log=True,screen=None,backup=True)
                rungablam = self.force()
                for gabfile in ['gablam','hitsum','local']: rungablam = rungablam or not rje.exists('%s.%s.tdt' % (vbase,gabfile))
                if rungablam:
                    bcmd = ['blaste=1e-10','keepblast=T','fullblast=T']
                    gcmd = ['local=F','seqin=%s' % vfullfile,'searchdb=%s' % vfullfile,'blastp=blastn','saveupc=T',
                            'qryacc=F','basefile=%s' % vbase,'blastdir=%sVBLAST/' % self.getStr('RevertDir')]
                    gablam.GABLAM(self.log,bcmd+self.cmd_list + gcmd).run()
                vgabdb = self.db().addTable('%s.gablam.tdt' % vbase,mainkeys=['Qry','Hit'])
                ## Viral Proteomes ##
                vbase = '%s%s.vproteome' % (vgabdir,self.basefile(strip_path=True))
                vprotfile = '%s.fas' % vbase
                if self.force() or not rje.exists(vprotfile):
                    vfullseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=list'])
                    for vfile in self.list['VBatch']:
                        vfull = '%s%s.prot.fas' % (self.getStr('VirusDir'),vfile)
                        vfullseq.loadSeq(vfull,nodup=True,clearseq=False,mode='tuple')
                    vfullseq.saveSeq(seqfile=vprotfile,reformat=None,append=False,log=True,screen=None,backup=True)
                rungablam = self.force()
                for gabfile in ['gablam','hitsum']: rungablam = rungablam or not rje.exists('%s.%s.tdt' % (vbase,gabfile))
                if rungablam:
                    bcmd = ['blaste=1e-10','keepblast=T','fullblast=T']
                    gcmd = ['local=F','seqin=%s' % vprotfile,'searchdb=%s' % vprotfile,'blastp=blastp','saveupc=T',
                            'qryacc=F','basefile=%s' % vbase,'blastdir=%sVBLAST/' % self.getStr('RevertDir')]
                    gablam.GABLAM(self.log,bcmd+self.cmd_list + gcmd).run()

            ### ~ [2] ~ XGMML network ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Need to be wary of multiple viruses with same alias. Currently these will map to same nodes. #!#
            #!# Also true of genomes with the same alias! This is a "feature" that allows one to combine results! :op
            #!# Maybe add a warning and option to rejig aliases? Want to do this at the initial compression, which is
            #!# currently in the batchRun() method (following dev call) but needs to be moved
            ppi = self.obj['PPI'] = rje_ppi.PPI(self.log,self.cmd_list)
            ppi.setNum({'Walltime':min(ppi.getNum('Walltime'),0.2)})
            ppi.basefile(self.basefile())
            ppi.obj['DB'] = self.obj['DB']
            ndb = db.addEmptyTable('Node',['Node','Name','Type','Size'],['Node'])
            #Add VirusGB from VBatch (virus Table)
            #for vgb in vsumdb.indexKeys('VirusGB'): ndb.addEntry({'Node':vgb,'Name':vgb,'Type':'VirusGB','Size':len(vsumdb.index('VirusGB')[vgb])})

            #i# VirusID and Genome are already aliases
            #!# Need to add VirusGB info back in for network generation

            for entry in vsumdb.entries(): ndb.addEntry({'Node':entry['VirusID'],'Name':entry['VirusID'],'Type':'Virus','Size':entry['GFile']})
            for vacc in rdb.index('VirusAcc'):
                virus = self.alias(vacc)
                if virus not in ndb.data(): ndb.addEntry({'Node':virus,'Name':virus,'Type':'Virus','Size':0})
            for entry in gsumdb.entries(): ndb.addEntry({'Node':entry['Genome'],'Name':entry['Genome'],'Type':'Genome','Size':entry['Virus']})
            for gfile in rdb.index('GFile'):
                genome = self.alias(gfile)
                if genome not in ndb.data(): ndb.addEntry({'Node':genome,'Name':genome,'Type':'Genome','Size':0})
            edb = db.copyTable(vdb,'Edge')
            edb.renameField('VirusID','Hub'); edb.renameField('Genome','Spoke')
            edb.compress(['Hub','Spoke'],rules={'VirusAcc':'list','Virus':'list','VirusGB':'list'},default='max')
            edb.addField('Type',evalue='REVERT'); edb.addField('Evidence',evalue='REVERT');
            edb.dropEntriesDirect('FragNum',[0])
            for vfile in vlocdb.index('vfiles',splitchar='|'):
                self.deBug(vfile)
                ndb.addEntry({'Node':vfile,'Name':vfile,'Type':'VirusGB','Size':len(vlocdb.index('vfiles')[vfile])})
                for vacc in vlocdb.index('vfiles')[vfile]:
                    self.debug('%s = %s' % (vacc,self.alias(vacc)))
                    edb.addEntry({'Hub':vfile,'Spoke':self.alias(vacc),'Evidence':'GenBank','Type':'VirusGB'})
            self.debug(self.db('Node').dataKeys())
            for edge in edb.entries():
                self.debug(edge)
                #edge['Hub'] = self.alias(edge['Hub'])
                if edge['Hub'] not in ndb.data(): raise ValueError('Hub %s missing from Nodes!' % edge['Hub'])
                #edge['Spoke'] = self.alias(edge['Spoke'])
                if edge['Spoke'] not in ndb.data(): raise ValueError('Spoke %s missing from Nodes!' % edge['Spoke'])
            edb.remakeKeys()

            #!# Add GABLAM links #!#
            if self.getBool('VGABLAM'):
                vgabdb.renameField('Qry','Hub'); vgabdb.renameField('Hit','Spoke')
                for entry in vgabdb.entries():
                    entry['Hub'] = self.alias(string.split(entry['Hub'],'__')[-1])
                    entry['Spoke'] = self.alias(string.split(entry['Spoke'],'__')[-1])
                    entry['Type'] = 'GABLAM'; entry['Evidence'] = 'GABLAM'
                    entry['Identity'] = (float(entry['Qry_AlnID'])+float(entry['Hit_AlnID']))/2.0
                    edb.addEntry(entry)

            ## VirusGB-Virus links
            #for vgb in vsumdb.indexKeys('VirusGB'):
            #    for virus in vsumdb.indexDataList('VirusGB',vgb,'Alias'):
            #        edb.addEntry({'Hub':vgb,'Spoke':virus,'Type':'VirusGB','Identity':100.0,'Evidence':'Genbank'})
            #        #if not self.dev(): ppi.addPPI(vgb,virus,evidence=1,asdict=None,sym=True)
            if self.dev(): ndb.saveToFile()
            if self.dev(): edb.saveToFile()
            self.debug(self.db('Node').dataKeys())
            ## Virus-Genome links ##
            ppi.ppiFromEdges(sym=False,clear=True)
            vglist = rje.sortUnique(edb.indexKeys('Hub')+edb.indexKeys('Spoke'))
            #vglist = vsumdb.indexKeys('VirusGB')
            #for entry in vdb.entries():
            #    if not entry['Identity']: continue
            #    (hub,spoke) = (adb.data(entry['Virus'])['Alias'],adb.data(entry['Genome'])['Alias'])
            #    if hub not in vglist: vglist.append(hub)
            #    if spoke not in vglist: vglist.append(spoke)
            #    #if not self.dev(): ppi.addPPI(hub,spoke,evidence=entry['Identity'],asdict=None,sym=True)
            fullG = ppi.dict['PPI'] # Main PPI dictionary: {hub:{spoke:evidence}}
            if self.dev():
                #!# Not really sure what the use of this is!
                xgmml = ppi.ppiToXGMML(fullG,'full.%s' % self.basefile())
                xgmml.saveXGMML('%s.dev.full.xgmml' % self.basefile())
            G = rje_ppi.subGraph(ppi.dict['PPI'],vglist)
            G = ppi.symmetry(G,inplace=False)
            G = rje_ppi.subGraph(G,rje_ppi.kCore(G,k=1))     # ppi.purgeOrphans() on subgraph
            if not G and self.dev(): self.printLog('#GRAPH','No virus-genome hits for graph/tree output.')
            if self.list['GraphFormats'] and G:
                npos = ppi.rjeSpringLayout(G)
                ppi.addCol(coldict={'VirusGB':6,'Virus':5,'Genome':3},ckey='Type')
                if 'svg' in self.list['GraphFormats'] or 'html' in self.list['GraphFormats']:
                    svghtm = ppi.saveSVG(npos,basefile='%s.revertnr' % self.basefile(),G=G,font=0,width=1600,ntype='ellipse',backups=True)
                if 'png' in self.list['GraphFormats'] or 'html' in self.list['GraphTypes']:
                    ppi.saveR(npos,basefile='%s.revertnr' % self.basefile(),G=G,cleantdt=not self.dev(),backups=True)
                #!# Add EdgeTypes for XGMLL output.
                if 'xgmml' in self.list['GraphFormats']:
                    xgmml = ppi.ppiToXGMML(G,self.basefile())
                    xgmml.saveXGMML('%s.revertnr.xgmml' % self.basefile())
                if 'html' in self.list['GraphFormats']:
                    import rje_html
                    html = rje_html.htmlHead(self.info['Name'],stylesheets=[],tabber=False) + svghtm
                    html += '\n<hr>\n<img src=%s.png width=1600>\n' % self.basefile()
                    html += rje_html.htmlTail(copyright='RJ Edwards 2014',tabber=False)
                    open('%s.htm' % self.basefile(),'w').write(html)
                self.printLog('#OUT','REVERT Virus-Genome graph output complete')
            ## ~ [4b] ~ Virus Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['TreeFormats'] and G and vsumdb.entryNum() > 2 and self.dev():
                vtree = rje_tree.Tree(self.log,['treeformats=nsf,png']+self.cmd_list+['autoload=F'])
                vdm = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
                vdis = {}
                for entry in vdb.entries():
                    if entry['VirusID'] not in vdis: vdis[entry['VirusID']] = {}
                    if not entry['FragNum']: continue
                    dis1 = vdis[entry['VirusID']]
                    if entry['Genome'] in dis1: dis1[entry['Genome']] = max(dis1[entry['Genome']],entry['Identity'])    # Might have several versions of same virus
                    else: dis1[entry['Genome']] = entry['Identity']
                for virus1 in vsumdb.dataKeys():
                    if self.alias(virus1) not in G: continue    # Exclude viruses w/o hits from tree
                    id1 = virus1    # adb.data(virus1)['Alias']
                    dis1 = vdis[virus1]
                    for virus2 in vsumdb.dataKeys():
                        if self.alias(virus2) not in G: continue    # Exclude viruses w/o hits from tree
                        id2 = virus2    #adb.data(virus2)['Alias']
                        dis2 = vdis[virus2]
                        pdis = []
                        for genome in rje.listUnion(dis1.keys(),dis2.keys()):
                            if genome in dis1 and genome in dis2 and max(dis1[genome],dis2[genome]):
                                pdis.append((max(dis1[genome],dis2[genome])-min(dis1[genome],dis2[genome]))/max(dis1[genome],dis2[genome]))
                            elif genome in dis1 and genome in dis2: continue    # Both zero: ignore
                            else: pdis.append(1.0)
                        if pdis: vdm.addDis(id1,id2,sum(pdis)/len(pdis))
                        else: vdm.addDis(id1,id2,1.0)
                if vdm.objNum() > 1: vdm.savePNG(vtree,'%s.viruses' % self.basefile(),singletons=False)
            ## ~ [4c] ~ Genome Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['TreeFormats'] and G and gsumdb.entryNum() > 2 and self.dev():
                gtree = rje_tree.Tree(self.log,['treeformats=nsf,png']+self.cmd_list+['autoload=F'])
                gdm = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
                gdis = {}
                for entry in vdb.entries():
                    if entry['Genome'] not in gdis: gdis[entry['Genome']] = {}
                    if not entry['FragNum']: continue
                    dis1 = gdis[entry['Genome']]
                    if entry['VirusID'] in dis1: dis1[entry['VirusID']] = max(dis1[entry['VirusID']],entry['Identity'])    # Might have several versions of same virus
                    else: dis1[entry['VirusID']] = entry['Identity']
                for genome1 in gsumdb.dataKeys():
                    id1 = genome1   #adb.data(genome1)['Alias']
                    dis1 = gdis[genome1]
                    for genome2 in gsumdb.dataKeys():
                        id2 = genome2   #adb.data(genome2)['Alias']
                        dis2 = gdis[genome2]
                        pdis = []
                        for virus in rje.listUnion(dis1.keys(),dis2.keys()):
                            if virus in dis1 and virus in dis2 and max(dis1[virus],dis2[virus]):
                                pdis.append((max(dis1[virus],dis2[virus])-min(dis1[virus],dis2[virus]))/max(dis1[virus],dis2[virus]))
                            elif virus in dis1 and virus in dis2: continue    # Both zero: ignore
                            else: pdis.append(1.0)
                        if pdis: gdm.addDis(id1,id2,sum(pdis)/len(pdis))
                        else: gdm.addDis(id1,id2,1.0)
                if gdm.objNum() > 1: gdm.savePNG(gtree,'%s.genomes' % self.basefile(),singletons=False)
            ## ~ [4d] ~ Summarise failures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Should this be independent of G?
            #self.debug(G)
            if G:
                rnull = {}  # Dictionary of {type:null list}
                for rkey in ['VirusID','Genome']:
                    rnull[rkey] = []
                    for rname in rje.sortKeys(vdb.index(rkey)):
                        if rname in G or self.alias(rname) in G: continue
                        rnull[rkey].append(rname)
                        self.printLog('#NULL','No REVERT reconstructions for %s "%s".' % (rkey,rname))
            else: self.printLog('#NULL','No REVERT reconstructions!')
            ## ~ [4e] ~ Add combined HTML output. Move to new directory? ~~~~~~~~~ ##
            return True
        except: self.errorLog('Problem with REVERT.batchOutput()'); return False
#########################################################################################################################
    ### <4> ### Special Viral Genome ID parsing                                                                         #
#########################################################################################################################
    def vgbParse(self):     # Parses viral IDs and hosts from vbatch files and outputs *.acc files
        '''Parses viral IDs and hosts from vbatch files and outputs *.acc files.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['VBatch']:
                if not self.getStrLC('VirusGB'): raise IOError('Cannot find vbatch or virusgb files!')
                self.list['VBatch'] = [self.getStr('VirusGB')]
            self.log.no_suppression += ['entry_overwrite','Invented SpCode']
            self.log.warnings += ['SpCode Missing TaxID']
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            if not self.getStrLC('Basefile'): self.setBaseFile('revert')
            ### ~ [1] ~ Parse Genbank viral genome info downloads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gbhead = ["Representative","Neighbor","Host","Lineage","Taxonomy","Segment"]
            gbkeys = ["Representative",'Neighbor']
            gbdb = db.addTable(filename=self.list['VBatch'][0],name='genomes',mainkeys=gbkeys,headers=gbhead,expect=True)
            for gbfile in self.list['VBatch'][1:]:
                db.mergeTables(gbdb,db.addTable(filename=gbfile,name='tmp',mainkeys=gbkeys,headers=gbhead,expect=True),overwrite=False,matchfields=True)
            gbdb.index('Representative')
            vhosts = rje.sortKeys(gbdb.index('Host',splitchar=','))
            gbdb.indexReport('Host')
            ### ~ [2] Generate output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['VHost']: self.list['VHost'] = vhosts
            for host in ['all'] + self.list['VHost']:
                if host.lower() == 'none': continue
                if host.lower() == 'all':
                    viruses = rje.sortKeys(gbdb.index('Representative'))
                else:
                    if host not in vhosts: self.warnLog('No viruses found for host "%s".' % host); continue
                    viruses = gbdb.indexDataList('Host',host,'Representative')
                hfile = '%s.%s.acc' % (self.basefile(),host)
                rje.backup(self,hfile)
                open(hfile,'w').write(string.join(viruses+[''],'\n'))
                self.printLog('#ACC','%s representative Genbank IDs output to %s.' % (rje.iLen(viruses),hfile))
            diffx = rje.listDifference(self.list['VHost'],vhosts)
            if diffx: self.printLog('#VHOST','%s vhost=LIST hosts not found in vbatch files.' % rje.iStr(diffx))
            diffx = rje.listDifference(vhosts,self.list['VHost'])
            if diffx: self.printLog('#VHOST','%s vbatch hosts excluded by vhost=LIST.' % rje.iStr(diffx))
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    ### <5> ### Species Alias Methods                                                                                   #
#########################################################################################################################
    def alias(self,name):   ### Returns output alias for name (if there is one) else name itself
        '''Returns output alias for name (if there is one) else name itself.'''
        adb = self.db('alias')
        if not adb: return name
        if name in adb.data(): return adb.data(name)['Alias']
        return name
#########################################################################################################################
    def setupAlias(self,vdb=None):  ### Sets up alias file and, if required, taxonomy data
        '''
        Sets up alias file and, if required, taxonomy data.
        >> vdb:Table [None] = Optional revert table from which to extract 'Virus' and 'Genome' entries.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~~~~~ ALIAS SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            db = self.db()
            ### ~ [1] ~ Check/Load ALias Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('AliasFile'): self.setStr({'AliasFile':'%s.alias.tdt' % self.basefile()})
            if self.getBool('Aliases'): adb = db.addTable(self.getStr('AliasFile'),['Name'],name='alias',expect=False)
            else: adb = None
            if not adb:
                afile = '%s.alias.tdt' % self.basefile()    # Look for basefile.alias.tdt if given file not found
                if self.getBool('Aliases') and rje.exists(afile) and self.i() >= 0 and rje.yesNo('Use %s for aliases?' % afile,default='Y'):
                    adb = db.addTable(afile,['Name'],name='alias',expect=True)
                else: adb = db.addEmptyTable('alias',['Name','Alias'],['Name'])
            return adb
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    def makeAliases(self):  ### Generates/Checks Alias data for all viruses and genomes. Sets up alias file if required.
        '''
        Generates/Checks Alias data for all viruses and genomes. Sets up alias file if required. This is now run once for
        all results.
        @return: True/False depending on success/failure.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~ MAKING ALIASES ~~~~~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            ## ~ [0a] ~ Setup tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rdb = self.db('revert.details',add=False)
            vdb = self.db('revert',add=False)
            if not rdb and vdb: return False
            adb = self.db('alias')
            if not adb: adb = self.setupAlias()
            ## ~ [0b] ~ Check/Use existing alias data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# VirusAcc for viruses and GFile for genomes are unique & used for NR compression
            #i# Virus, VirusID and Genome are user-friendly aliases.
            #i# VirusID and Genome are used for all-by-all output
            #i# The "Alias" field for viral data is no longer used: use VirusID instead.
            needtax = False     # Whether taxonomy data will be required as new aliases are needed.
            for table in [rdb,vdb]:
                table.makeField(formula='#Virus#',fieldname='VirusID',after='Virus')
                for entry in table.entries():
                    vacc = entry['VirusAcc']
                    if vacc in adb.data(): entry['VirusID'] = self.alias(vacc)
                    else: needtax = True
                    genome = entry['GFile']
                    if genome in adb.data(): entry['Genome'] = self.alias(genome)
                    else: needtax = True
                table.remakeKeys()
            review = self.i() >= 0 and rje.yesNo('Review Virus/Genome aliases used for batch summary and all-by-all output?',default='N')
            ## ~ [0c] ~ Check/Load taxonomy data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if needtax or review:
                tax = self.obj['Taxonomy']  #? Why is this not working?
                if not tax.getBool('Setup'): tax.setup()
            else: return True

            #!# Add checks for duplicating existing aliases when adding new ones!

            ### ~ [1] ~ Process Virus Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Now making sure that each VirusAcc has an alias as VirusID field.
            for virus in vdb.indexKeys('Virus'):
                ventries = vdb.indexDataList('Virus',virus,'VirusAcc')
                if len(ventries) > 1: self.warnLog('Virus "%s" mapped to multiple VirusAcc' % virus)
                for vacc in ventries:
                    if not review and vacc in adb.data(): continue
                    if not review and virus in adb.data():
                        adb.addEntry({'Name':vacc,'Alias':self.alias(virus)})
                        continue
                    spcode = tax.getSpCode(virus)
                    usespcode = self.i() < 0 or (self.getBool('VSpCode') and self.i() < 1 and not review)
                    if usespcode: alias = spcode
                    else:
                        menulist = [('S','Species Code: %s' % spcode,'return',spcode),
                                    ('A','Accession Number: %s' % vacc,'return',vacc),
                                    ('N','Virus Name: %s' % virus,'return',virus),
                                    ('E','Enter Alias','return','ENTER'),
                                    ('X','No Alias','return','')]
                        mdefault = 'S'
                        if virus in adb.data():
                            menulist.append(('V','Use virus alias: %s' % adb.data(virus),'return',adb.data(virus)))
                            mdefault = 'V'
                        if review and vacc in adb.data():
                            menulist.append(('K','Keep existing alias: %s' % adb.data(vacc),'return',adb.data(vacc)))
                            mdefault = 'K'

                        alias = rje_menu.menu(self,'Alias options for %s (%s):' % (virus,vacc),menulist,choicetext='Please select:',changecase=True,default=mdefault)
                        #x#alias = rje_menu.menu(self,'Alias options for %s:' % (virus),menulist,choicetext='Please select:',changecase=True,default='S')
                        if alias == 'ENTER': alias = rje.choice('Enter new Alias:',confirm=True)
                    if alias and alias != adb.data(vacc):
                        #?# Consider not warning about overwrites and adding bespoke message.
                        self.printLog('#ALIAS','%s (%s) alias updated from "%s" to "%s"' % (vacc,virus,adb.data(vacc),alias))
                        adb.addEntry({'Name':vacc,'Alias':alias},warn=False)

            ### ~ [2] ~ Process Genome Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Now making sure that each GFile has an alias as Genome field.
            for genome in vdb.indexKeys('GFile'):
                if genome not in adb.data():
                    if not review and genome in adb.data(): continue
                    # Setup possible Genome aliases
                    taxon = string.split(genome,'.')[0]
                    if taxon.count('_') == 1:
                        spcode = tax.getSpCode(taxon,invent=False)
                        spabbrev = string.split(taxon,'_')
                        try: spabbrev = '%s.%s' % (spabbrev[0].upper()[:1],spabbrev[1].lower())
                        except: spabbrev = spabbrev[0]
                    else:
                        spcode = tax.getSpCode(string.replace(string.split(taxon,'_')[0],'-',' '),invent=False)
                        spabbrev = string.split(string.split(taxon,'_')[0],'-')[:2]
                        try: spabbrev = '%s.%s' % (spabbrev[0].upper()[:1],spabbrev[1].lower())
                        except: spabbrev = spabbrev[0]
                    usespcode = self.getBool('GSpCode') and (self.i() < 0 or (self.i() < 1 and not review))
                    usespcode = usespcode and spcode
                    if usespcode: alias = spcode
                    elif self.i() < 0 or (self.i() < 1 and not review): alias = spabbrev
                    else:
                        menulist = [('S','Species Code: %s' % spcode,'return',spcode),
                                    ('G','Genome file: %s' % genome,'return',genome),
                                    ('N','Parsed Organism Name: %s' % spabbrev,'return',spabbrev),
                                    ('E','Enter Alias','return','ENTER'),
                                    ('X','No Alias','return','')]
                        mdefault = 'N'
                        if spcode:
                            if self.getBool('GSpCode'): mdefault = 'S'
                        else: menulist.pop(0)
                        if review and genome in adb.data():
                            menulist.append(('K','Keep existing alias: %s' % adb.data(genome)))
                            mdefault = 'K'

                        alias = rje_menu.menu(self,'Alias options for %s:' % (genome),menulist,choicetext='Please select:',changecase=True,default=mdefault)
                        if alias == 'ENTER': alias = rje.choice('Enter new Alias:',confirm=True)
                    if alias and alias != adb.data(genome):
                        self.printLog('#ALIAS','%s alias updated from "%s" to "%s"' % (genome,adb.data(genome),alias))
                        adb.addEntry({'Name':genome,'Alias':alias},warn=False)

            ### ~ [3] ~ Update Aliases and Save updated table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb.saveToFile(backup=self.getBool('AliasBackup'))
            self.setBool({'AliasBackup':False})
            #!# Make sure that VirusID replaces Alias in code #!#
            for table in [rdb,vdb]:
                for entry in table.entries():
                    vacc = entry['VirusAcc']
                    entry['VirusID'] = self.alias(vacc)
                    genome = entry['GFile']
                    entry['Genome'] = self.alias(genome)
                table.remakeKeys()
            return True
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    ### <6> ### RevertNR combining and filtering                                                                        #
#########################################################################################################################
    def revertNR(self,gfile,grdb,append=False):   ### Performs non-redundancy and quality filtering for combined data of a single genome
        '''Performs non-redundancy and quality filtering for combined batch data.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~ REVERT NR/QC FILTERING ~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            self.printLog('#GENQC','%s: %s' % (self.alias(gfile),gfile))
            gbase = gfile; gfile = None
            #!# Tidy this up a bit? Should gfile come directly as file name?
            for gseqfile in self.list['GBatch']:
                if gbase == rje.baseFile(gseqfile,strip_path=True): gfile = gseqfile
            if not gfile: self.warnLog('Cannot find %s in GBatch!' % gbase); return False
            rungablam = self.force()
            db = self.db()
            #rhead = ['VirusAcc','Virus','Genome','Protein','ProtAcc','Product','HitNum','MaxScore','EVal','Length','FragNum','Coverage','Identity','Local']
            #rkeys = ['VirusAcc','Genome','Protein']
            basefile = os.path.basename(self.baseFile())
            #gbase = rje.baseFile(gfile,strip_path=True)
            rdb = grdb #self.db('revert.details')
            rdb.index('ProtAcc')
            #rdb = db.copyTable(rdb,gbase)
            #rdb.dropEntriesDirect('GFile',gfile,inverse=True,log=self.debugging())
            nrdir = rje.makePath('%sNR/%s/' % (self.getStr('RevertDir'),basefile))
            rje.mkDir(self,nrdir)
            nrbase = '%s%s' % (nrdir,gbase) #os.path.basename(self.baseFile()))

            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #1. Use the *.local.tdt files to assemble fasta files of (a) combined non-overlapping genome fragments and (b) combined viral proteins.
            vprotseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=list'])
            virusdir = self.getStr('VirusDir')
            # Load local tables into database for initial protein culling based on mincov and minlocid using rdb
            ldb = None; gdb = None; qcx = 0
            for virus in self.list['VBatch']:
                vfas = '%s%s.prot.fas' % (virusdir,virus)
                vfasseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % vfas,'seqmode=file'])
                vseqdict = vfasseq.makeSeqNameDic('accnum')
                blocal = '%sGABLAM/%s.%s.local.tdt' % (self.getStr('RevertDir'),virus,gbase)
                if not ldb: mdb = db.addTable(blocal,mainkeys=['Qry','Hit','QryStart','SbjStart','SbjEnd'],name='local')
                else: mdb = db.addTable(blocal,mainkeys=['Qry','Hit','QryStart','SbjStart','SbjEnd'],name='merge')
                #x#if not self.getStrLC('Rest'): mdb.keepFields(['Qry','Hit','QryStart','SbjStart','SbjEnd'])
                mdb.dataFormat({'QryStart':'int','SbjStart':'int','SbjEnd':'int'})
                mdb.addField('GenHit')
                for mentry in mdb.entries(): mentry['GenHit'] = '%s|%s' % (gfile,mentry['Hit'])

                for vprot in mdb.indexKeys('Qry'):
                    if vprot in rdb.index('ProtAcc'):   # Note that GABLAM must be using qryacc=T
                        (name,sequence) = vfasseq.getSeq(vseqdict[vprot])
                        vprotseq._addSeq(name,sequence,nodup=True)
                    else:
                        for ekey in mdb.index('Qry')[vprot]: mdb.dict['Data'].pop(ekey); qcx += 1

                if ldb: db.mergeTables(ldb,mdb,overwrite=False,matchfields=True)
                else: ldb = mdb

                #bgablam = '%s%s.%s.gablam.tdt' % (self.getStr('RevertDir'),rje.baseFile(vfile,strip_path=True),rje.baseFile(gfile,strip_path=True))
                #if not gdb: mdb = db.addTable(bgablam,mainkeys=['Qry','Hit'],name='gablam')
                #else: mdb = db.addTable(bgablam,mainkeys=['Qry','Hit'],name='merge')
                #mdb.dataFormat({'Qry_OrderedAlnID':'num','Qry_OrderedAlnLen':'num'})
                #for gkey in mdb.dataKeys():
                #    entry = mdb.data(gkey)
                #    if entry['Qry_OrderedAlnLen'] < self.getNum('MinCov'): mdb.dict['Data'].pop(gkey); qcx += 1
                #    elif 100.0 * entry['Qry_OrderedAlnID'] / entry['Qry_OrderedAlnLen'] < self.getNum('MinLocID'): mdb.dict['Data'].pop(gkey); qcx += 1

                #!# Filtered these in the batchRun() method and this rdb is pre-filtered #!#
                #!# Just need to get the Queries: can come from ldb
                #!# Remember to delete tables
                #bgablam = '%s%s.%s.hitsum.tdt' % (self.getStr('RevertDir'),rje.baseFile(vfile,strip_path=True),rje.baseFile(gfile,strip_path=True))
                #if not gdb: mdb = db.addTable(bgablam,mainkeys=['Qry'],name='gablam')
                #else: mdb = db.addTable(bgablam,mainkeys=['Qry'],name='merge')
                #mdb.dataFormat({'Coverage':'int','Identity':'int','Length':'int'})
                #for gkey in mdb.dataKeys():
                #    entry = mdb.data(gkey)
                #    if not entry['Coverage']: mdb.dict['Data'].pop(gkey); qcx += 1
                #    elif 100.0 * entry['Coverage'] / entry['Length'] < self.getNum('MinCov'): mdb.dict['Data'].pop(gkey); qcx += 1
                #    elif 100.0 * entry['Identity'] / entry['Coverage'] < self.getNum('MinLocID'): mdb.dict['Data'].pop(gkey); qcx += 1

                #for vprot in mdb.indexKeys('Qry'):
                #    (name,sequence) = vfasseq.getSeq(vseqdict[vprot])
                #    vprotseq._addSeq(name,sequence,nodup=True)
                #if gdb: db.mergeTables(gdb,mdb,overwrite=False,matchfields=True)
                #else: gdb = mdb

            #self.printLog('#QC','%s Local hits removed (low protein/virus coverage and/or identity): %s remain' % (rje.iStr(qcx),rje.iStr(ldb.entryNum())))
            if not ldb.entryNum(): return False
            #prex = ldb.entryNum()


            #for lkey in ldb.dataKeys():
            #    entry = ldb.data(lkey)
            #    gkey = gdb.makeKey(entry)
            #    if not gdb.data(gkey): ldb.dict['Data'].pop(lkey)
            #self.printLog('#QC','%s local hits -> %s' % (rje.iStr(prex),rje.iStr(ldb.entryNum())))
            #if not ldb.entryNum(): return False
            # Compress to non-overlapping genome fragments & viral proteins
            redx = 0; nrx = 0
            gfrags = []                      # List of (contig,start,stop)
            for ghit in ldb.indexKeys('GenHit'):
                #self.debug(ghit)
                (genome,contig) = string.split(ghit,'|',1)
                frags = []
                for lentry in ldb.indexEntries('GenHit',ghit):
                    frags.append((min(lentry['SbjStart'],lentry['SbjEnd']),max(lentry['SbjStart'],lentry['SbjEnd'])))
                    redx += 1
                frags.sort()
                #self.debug(frags)
                i = 1
                while (i < len(frags)):
                    #self.bugPrint('%d : %s' % (i,frags[i-1:i+1]))
                    if frags[i][0] <= frags[i-1][1] + max(1,self.getInt('GablamFrag')+1): frags[i-1] = (frags[i-1][0],max(frags[i-1][1],frags.pop(i)[1]))    # Merge
                    else: i += 1
                #self.debug(frags)
                for frag in frags: gfrags.append((contig,frag[0],frag[1])); nrx += 1
            self.printLog('#FRAG','%s fragments of %s reduced to %s NR fragments.' % (rje.iStr(redx),gbase,rje.iStr(nrx)))
            ldb.dropField('GenHit')
            ldb.dataFormat({'QryEnd':'int','AlnNum':'int'})
            ldb.newKey(['Qry','Hit','AlnNum'])
            locfile = '%s.local.tdt' % nrbase
            ldb.saveToFile(locfile)
            db.deleteTable(ldb)
            # Extract sequences into fasta files
            vprotfas = '%s.prot.fas' % nrbase
            if self.force() or not rje.exists(vprotfas) or rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % vprotfas,'seqmode=file','seqindex=T']).seqNum() != vprotseq.seqNum():
                rungablam = True
                vprotseq.saveSeq(seqfile=vprotfas)

            gfragfas = '%s.gfrag.fas' % nrbase
            #!# Need to merge these gfrag.fas files
            gseq = None
            if os.path.exists(gfragfas) and not self.force() and open(gfragfas,'r').readline()[:1] == '>':
                gseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % gfragfas,'seqmode=file'])
                self.printLog('\r#FAS','%s fragment sequences read from %s.' % (rje.iStr(gseq.seqNum()),gfragfas))
                if gseq.seqNum() != nrx:
                    self.printLog('#FAS','%s exists but contains wrong number of sequences. Will remake.' % gfragfas); gseq = None
                else: self.printLog('#FAS','%s exists and contains right number of sequences (force=F).' % gfragfas)
            if not gseq:
                self.progLog('#FAS','Fragment sequences output to %s...' % (gfragfas))
                GFRAG = open(gfragfas,'w'); sx = 0
                genome = rje.baseFile(gfile,strip_path=True)
                gseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % gfile,'seqmode=file'])
                gdict = gseq.makeSeqNameDic('short')
                (prev,name,sequence) = (None,None,None)
                for (contig,startx,stopx) in gfrags:
                     if contig != prev: (name,sequence) = gseq.getSeq(gdict[contig])
                     prev = string.split(name)[0]; desc = string.join(string.split(name)[1:])
                     GFRAG.write('>%s.%s-%s %s|(Pos: %s..%s)\n' % (prev,rje.preZero(startx,len(sequence)),rje.preZero(stopx,len(sequence)),desc,rje.iStr(startx),rje.iStr(stopx)))
                     GFRAG.write('%s\n' % sequence[startx-1:stopx]); sx += 1
                GFRAG.close()
                self.printLog('\r#FAS','%s fragment sequences output to %s.' % (rje.iStr(sx),gfragfas))
                gseq = None

            # All-by-all vprot GABLAM -> assess which multiple hits to judge against each other
            vbase = rje.baseFile(vprotfas)
            if rungablam and rje.exists('%s.blast' % vbase): os.unlink('%s.blast' % vbase) # Remade sequences
            if rungablam: gseq = None
            for gabfile in ['gablam','hitsum']: rungablam = rungablam or not rje.exists('%s.%s.tdt' % (vbase,gabfile))
            if rungablam:
                bcmd = ['blaste=1e-10','keepblast=T']
                gcmd = ['local=F','seqin=%s' % vprotfas,'searchdb=%s' % vprotfas,'blastp=blastp','fullblast=T',
                        'qryacc=F','basefile=%s' % vbase,'blastdir=%s' % nrdir,'selfhit=T',
                        'dismat=F','distrees=F','disgraph=F','clusters=F','saveupc=F']
                gablam.GABLAM(self.log,bcmd+self.cmd_list + gcmd).run()
            vgab = db.addTable('%s.gablam.tdt' % (vbase),mainkeys=['Qry','Hit'],datakeys=['Qry','Hit'],name='vgablam')
            vgab.index('Qry')

            #2. Perform a genomefrag vs vprot blastx GABLAM search
            if not gseq and rje.exists('%s.blast' % nrbase): os.unlink('%s.blast' % nrbase) # Remade sequences
            bcmd = ['blaste=1e-10','fullblast=T','keepblast=T']
            gcmd = ['local=T','seqin=%s' % gfragfas,'searchdb=%s' % vprotfas,'blastp=blastx','saveupc=F','dismat=F',
                    'qryacc=F','basefile=%s' % nrbase,'blastdir=%s' % nrdir,'outstats=GABLAMO','percres=F']
            gablam.GABLAM(self.log,bcmd+self.cmd_list + gcmd).run()

            #3. Load GABLAM results for NR and further QC filtering
            gdb = db.addTable('%s.gablam.tdt' % nrbase,mainkeys=['Qry','Hit'],name='revertnr.best')
            gdb.dataFormat({'Qry_OrderedAlnID':'int','Hit_OrderedAlnID':'int','HitLen':'int','Hit_OrderedAlnLen':'int'})
            gdb.addField('Best',evalue=0)
            if not gdb.entryNum(): return False

            #4. Filter poorer hits using logPoisson(obs+1,best) >= 0.95? (i.e. < 5% chance of that few identities if actually the same)
            redx = 0#; redprot = []
            for gfrag in gdb.indexKeys('Qry'):
                # This is complicated by possibility of having unrelated proteins hitting same fragment side-by-side
                idhit = {}      # Dictionary of {IDPos:[Hits]}
                hitid = {}      # Dictionary of {Hit:ID}
                for entry in gdb.indexEntries('Qry',gfrag):
                    if entry['Qry_OrderedAlnID'] not in idhit: idhit[entry['Qry_OrderedAlnID']] = []
                    idhit[entry['Qry_OrderedAlnID']].append(entry['Hit'])
                    hitid[entry['Hit']] = entry['Qry_OrderedAlnID']
                # Process and reduce idhit and hitid dictionaries until all hits processed
                while idhit:
                    maxid = rje.sortKeys(idhit,revsort=True)[0]
                    goodid = maxid  # Min number of identities reached that are good
                    badid = 0       # Max number of identities reached that are bad
                    best = idhit[maxid].pop(0)
                    if not idhit[maxid]: idhit.pop(maxid)
                    if best not in hitid: continue  # Already processed
                    hom = vgab.indexDataList('Qry',best,'Hit')
                    for entry in gdb.indexEntries('Qry',gfrag):
                        if entry['Hit'] not in hitid: continue
                        if entry['Hit'] not in hom: continue
                        hitid.pop(entry['Hit'])
                        eid = entry['Qry_OrderedAlnID']
                        if eid == maxid: entry['Best'] = 1
                        elif eid <= badid: entry['Best'] = -1; redx += 1#; redprot.append(entry['Hit'])
                        elif eid >= goodid: entry['Best'] = 0
                        else:
                            if rje.poisson(eid+1,maxid,callobj=self) >= 0.95:   # Bad. Too few IDs to be chance missing out
                                badid = max(badid,eid)
                                entry['Best'] = -1#; redprot.append(entry['Hit'])
                                redx += 1
                            else: goodid = min(goodid,eid)
            self.printLog('#RED','%s hits rejected as redundant (better protein hit)' % rje.iStr(redx))
            gdb.index('Best',force=True)
            gdb.indexReport('Best')
            gdb.saveToFile('%s.best.tdt' % nrbase,append=append)
            gdb.dropEntriesDirect('Best',[-1])
            goodfrag = gdb.indexKeys('Qry',force=True)
            if not gdb.entryNum(): return False

            #5. Count the number of best genome hits per protein too - include in output
            gdb.renameField('Hit','ProtAcc')
            gdb.addField('Genome')
            for entry in gdb.entries():
                entry['ProtAcc'] = string.split(entry['ProtAcc'],'__')[-1]
                entry['Genome'] = self.alias(gbase)
            gdb.remakeKeys()
            gdb.keepFields(['Genome','ProtAcc','Qry','Best'],log=False)
            gdb.compress(['Genome','ProtAcc'],rules={'Best':'max'})
            gdb.dropField('Qry')
            #?# Is gdb output anywhere? Not here?
            # This actually marks the best score per protein per genome - previous is per genome region?


            #6. Regenerate the revert output with the reduced protein hits.
            rdb.addField('BestProt',evalue=0)   # Sum up number of best proteins
            rdb.addField('Best',evalue=-1)      # Max Best score for whole virus
            bestx = 0
            for entry in rdb.entries():
                #if float(entry['Local']) < self.getNum('MinLocID'): poorx += 1; continue
                #elif float(entry['Coverage']) < self.getNum('MinCov'): poorx += 1; continue
                gkey = gdb.makeKey(entry)
                gentry = gdb.data(gkey)
                if gentry: entry['BestProt'] = entry['Best'] = gentry['Best']; bestx += 1
            self.printLog('#BEST','Best Protein data added for %s of %s revert.details entries.' % (rje.iStr(bestx),rje.iStr(rdb.entryNum())))
            db.deleteTable(gdb)
            db.deleteTable(vgab)
            #self.printLog('#QC','%s protein hits rejected due to low coverage/identity' % rje.iStr(poorx))
            rdb.setStr({'Name':'revertnr.details'})
            badvirus = []   # List of viruses to filter
            for virus in rdb.indexKeys('VirusAcc'):     # Should this be virus?!
                if max(rdb.indexDataList('VirusAcc',virus,'Best',sortunique=False)) < 0: badvirus.append(virus)
            self.printLog('#REVNR','%s of %s viruses rejected: no quality proteins hits' % (rje.iLen(badvirus),rje.iLen(rdb.indexKeys('VirusAcc'))))
            rdb.dropEntriesDirect('VirusAcc',badvirus)
            if not rdb.entryNum(): return False

            #7. Filter according to viral coverage (40%) and local identity (30%)
            ddb = rdb
            rdb = db.copyTable(rdb,'revertnr')
            rdb.addField('ProtNum',evalue=1,after='Length')
            rdb.addField('ProtHits',evalue=1,after='ProtNum')
            for entry in rdb.entries():
                entry['Length'] = int(entry['Length'])
                entry['Coverage'] = entry['Length'] * float(entry['Coverage']) / 100.0
                entry['Identity'] = entry['Length'] * float(entry['Identity']) / 100.0
                #try: entry['ProtNum'] = vdb.indexDataList('VirusAcc',entry['VirusAcc'],'ProtNum')[0]
                #except: entry['ProtNum'] = 0
                if entry['Best'] < 0: entry['ProtAcc'] = ''; entry['ProtHits'] = 0
            rdb.compress(['VirusAcc','Genome'],rules={'ProtAcc':'list','Best':'max'},joinchar='|',default='sum')    #
            rdb.dropFields(['Protein','Product'],log=False)
            badvirx = 0   # Count of viruses to filter
            for entry in rdb.entries():
                if len(string.split(entry['ProtAcc'],'|')) != entry['ProtHits']: self.debug('Protein ID/number mismatch! %s vs %s' % (entry['ProtHits'],entry['ProtAcc']))
                entry['Local'] = 100.0 * entry['Identity'] / entry['Coverage']
                if entry['Local'] < self.getNum('MinVLocID'): entry['Best'] = -1; badvirx += 1; continue
                entry['Coverage'] = 100.0 * entry['Coverage'] / entry['Length']
                if entry['Coverage'] < self.getNum('MinVCov'): entry['Best'] = -1; badvirx += 1; continue
                entry['Identity'] = 100.0 * entry['Identity'] / entry['Length']
                for field in ['Identity','Coverage','Local']: entry[field] = '%.2f' % entry[field]
            self.printLog('#QC','%s viruses rejected due to low coverage/identity' % rje.iStr(badvirx))
            rdb.dropEntriesDirect('Best',[-1])
            ddb.dropEntriesDirect('VirusAcc',rdb.indexKeys('VirusAcc'),inverse=True)
            if not rdb.entryNum(): return False
            for table in [rdb,ddb]:
                table.saveToFile(append=append)
                db.deleteTable(table)


            ### [2] Generate Viral Assemblies and Consensus Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = ['autoload=T','autofilter=T','seqin=%s' % gfragfas,'goodseq=%s' % (string.join(goodfrag,','))]
            fragseq = rje_seqlist.SeqList(self.log,self.cmd_list+scmd)
            fragseq.saveSeq(seqfile='%s.fas' % nrbase)
            ## ~ [2a] Assembly BLASTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            assbase = '%s.%s' % (rje.baseFile(self.baseFile(),strip_path=True),gfile)
            self.virusAssembly('protein',basefile=assbase,vbase=nrbase)
            self.virusAssembly('host',basefile=assbase,vbase=nrbase)
            for virus in self.list['VBatch']:
                self.virusAssembly('genome',basefile=assbase,vbase=virus)
            ## ~ [2b] Protein consensus output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Why were we looking for this file and then overwriting the entries #?#
            cdb = self.db('revert.consensus',add=False,mainkeys=['Genome','Protein'])
            if not cdb:
                cdb = db.addEmptyTable('revert.consensus',['Genome','Protein','Desc','Consensus'],['Genome','Protein'])
                append = False
            vprotdict = vprotseq.makeSeqNameDic()
            for vprot in rje.sortKeys(vprotdict):
                paln = '%sASSEMBLYFAS/%s.%s.fas' % (self.getStr('RevertDir'),gbase,vprot)
                if os.path.exists(paln):
                    aseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % paln,'seqmode=tuple','autoload=T','autofilter=F'])
                    (qname,qseq) = aseqlist.seqs()[0]
                    consensus = []  #'.'] * len(qseq)
                    for i in range(len(qseq)):
                        consensus.append('.')
                        if qseq[i] == '-': consensus[-1] = '-'; continue
                        for hit in aseqlist.seqs()[1:]:
                            hseq = hit[1]
                            if hseq[i] == '-': continue
                            if hseq[i] == qseq[i]:
                                consensus[-1] = qseq[i]
                                break
                            consensus[-1] = 'x'     #?# Should an actual consensus be produced too? Harder!
                    consensus = string.join(consensus,'')
                    (pname,pdesc) = string.split(qname,maxsplit=1)
                    self.verbose(0,1,'%s %s\n%s' % (pname,pdesc,consensus))
                    cdb.addEntry({'Genome':gbase,'Protein':pname,'Desc':pdesc,'Consensus':consensus})
            if cdb.entryNum(): cdb.saveToFile(append=append)

            return True
        except: self.errorLog('REVERT.revertNR() error.'); return False
#########################################################################################################################
### End of SECTION II: REVERT Class                                                                                     #
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
    try: REVERT(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
    #mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
