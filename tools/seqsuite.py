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
Module:       SeqSuite
Description:  Genomics and biological sequence analysis toolkit
Version:      1.28.0
Last Edit:    19/09/24
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    SeqSuite is designed to be a front end for the expanding range of miscellaneous sequence utilities that are found
    within the SLiMSuite libraries/ modules. The relevant tool is given by the first system command, or selected using
    `prog=X` (or `program=X`). As much as possible, SeqSuite will emulate running that tool from the commandline, adding
    any matching `X.ini` file to the default commandline options read in (*before* settings read from seqsuite.ini
    itself). By default, the SeqList tool will be called (libraries/rje_seqlist.py) and read in commands from seqlist.ini.

    Help for the selected tool can be accessed using the `help=T` option. Note that `-h`, `-help` or `help` alone will
    trigger the SeqSuite help (this!). As `-help` or `help` will also set `help=T`, these commands will trigger both the
    SeqSuite help and the selected program help (unless over-ruled by `help=F`). An explicit `help=T` command will only
    trigger the selected program help.

SeqSuite tools:
    The list of tools recognised by `prog=X` will be added here as the relevant code is added:
    - Apollo = rje_apollo.Apollo GABLAM wrapper for searching against apollo genomes.
    - BLAST = rje_blast_V2.BLASTRun. BLAST+ Control Module.
    - BUSCOMP = buscomp.BUSCOMP. BUSCO Compiler and Comparison tool.
    - DBase = rje_dbase.DatabaseController. Database downloading and processing.
    - DepthKopy = depthkopy.DepthKopy. Single-copy read-depth based copy number analysis.
    - DepthSizer = depthsizer.DepthSizer. Read-depth based genome size prediction.
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

Special run modes:
    - summarise = run seqlist summarise on each file in `batchrun=FILELIST` [*.fasta]

Commandline:
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    prog=X      # Identifies the tool to be used. Will load defaults from X.ini (before seqsuite.ini) [seqlist]
    help=T/F    # Show the help documentation for program X. [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    batchrun=FILELIST   # Batch run the program on selected files (wildcards allowed) []
    batcharg=X          # Commandline argument to use for batchrun files ['seqin']
    batchlog=X          # Generate separate basefile.X log files for each batch run file (None for single log) [log]
    batchbase=T/F       # Whether to give each batch run a separate basefile=X command in place of log=X [True]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
### NOTE: INSTRUCTIONS TO ADD MODULE...
##  1. Add text, class and description to docstring above.
##  2. Add import statement for module.
##  3. Update mod dictionary
##  4. Add object creation to self.setup().
##  5. Check that the added class has a run() method.
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time, traceback
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
sys.path.append(os.path.join(slimsuitepath,'extras/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_zen
import rje_paf, rje_seq, rje_seqlist, rje_tree
import rje_dbase, rje_ensembl, rje_genbank, rje_gff, rje_mitab, rje_ppi, rje_pydocs, rje_taxonomy, rje_uniprot, rje_xref
# try:
#     import rje_dbase
# except:
#     rje.printf('Import error (python3?): {0}'.format(sys.exc_info()[0]))
#     rje.printf(traceback.extract_tb(sys.exc_info()[2]))
#     rje.printf('WARNING: some functions may not work')
import rje_blast_V2 as rje_blast
import buscomp, fiesta, gablam, gasp, gopher, haqesac, multihaq, seqmapper
import pingu_V4 as pingu
import pagsat, smrtscape, snapper, snp_mapper, rje_samtools, rje_apollo
import saaga, diploidocus, synbad, samphaser
import depthkopy, depthsizer
#########################################################################################################################
# Dev modules (not in main download):
try:
    extatic=None; revert=None; rje_genomics=None
    sys.path.append(os.path.join(slimsuitepath,'dev/'))
    import extatic, revert, rje_genomics
except:
    if len(sys.argv) != 2 or sys.argv[1] not in ['version', '-version', '--version', 'details', '-details', '--details', 'description', '-description', '--description']:
        print('Note: Failed to import dev modules etc. Some experimental features unavailable.\n|--------------------------------------|')
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added rje_seq and FIESTA. Added Uniprot.
    # 1.0 - Moved to tools/ for general release. Added HAQESAC and MultiHAQ. Moved mod to enable easy external access.
    # 1.1 - Added XRef = rje_xref.XRef. Identifier cross-referencing module.
    # 1.2 - Added taxonomy.
    # 1.3.0 - Added rje_zen.Zen. Modified code to work with REST services.
    # 1.4.0 - Added rje_tree.Tree, GABLAM and GOPHER.
    # 1.5.0 - Added extatic.ExTATIC and revert.REVERT. NOTE: Dev only.
    # 1.5.1 - Added 'seq' as alias for 'rje_seq' - want to avoid rje_ prefix requirements.
    # 1.6.0 - Added mitab and rje_mitab for MITAB parsing.
    # 1.6.1 - Added extra error messages.
    # 1.7.0 - Added pingu_V4.PINGU.
    # 1.8.0 - Added rje_pacbio.PacBio.
    # 1.9.0 - Added PAGSAT and SMRTSCAPE.
    # 1.9.1 - Fixed HAQESAC setobjects=True error.
    # 1.10.0 - Added batchrun=FILELIST batcharg=X batch running mode.
    # 1.11.0 - Added SAMTools and Snapper/SNP_Mapper.
    # 1.11.1 - Redirected PacBio to call SMRTSCAPE.
    # 1.11.2 - Fixed batchrun batchlog=False log error.
    # 1.12.0 - Added Snapper.
    # 1.12.1 - Fixed possible download run bug. (No dev/ modules?)
    # 1.13.0 - Changed GABLAM defaults to fullblast=T keepblast=T.
    # 1.14.0 - Added BLAST run object, particularly for qassemblefas running.
    # 1.14.1 - Added zentest for testing the REST servers.
    # 1.15.0 - Added GASP to REST servers.
    # 1.16.0 - Add rje_gff.GFF to REST servers.
    # 1.17.0 - Added batch summarise mode.
    # 1.18.0 - Added rje_apollo.Apollo to REST servers.
    # 1.19.0 - Tweaked the output of batch summarise, adding Gap% and reducing dp for some fields.
    # 1.19.1 - Fixed GapPC summarise output to be a percentage, not a fraction.
    # 1.20.0 - Added rje_paf.PAF.
    # 1.21.0 - Added NG50 and LG50 to batch summarise.
    # 1.22.0 - Added BUSCOMP to programs.
    # 1.23.0 - Added rje_ppi.PPI to programs.
    # 1.23.1 - Fixed GCPC bug in summarise().
    # 1.23.2 - Dropped pacbio as synonym for smrtscape. (Causing REST server issues.)
    # 1.24.0 - Added SAAGA, Diploidocus, SynBad and rje_genomics to programs.
    # 1.25.0 - Added SeqMapper and SAMPhaser.
    # 1.26.0 - Added DepthSizer. Added gzip support to batch summarise.
    # 1.27.0 - Added DepthKopy.
    # 1.27.1 - Py3 fixes & fudges to get some programs working.
    # 1.27.2 - Made batchrun=*.fasta the default for summarise.
    # 1.28.0 - Multiple updates to genomics tools. See release notes and indiviudal GitHubs.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program with SeqList & DBase.
    # [ ] : Add testing options.
    # [ ] : Add menu-driven operation.
    # [ ] : Make sure that programs run use the appropriate log if read from default ini file.
    # [ ] : Add seqbatch=FILELIST and batchcmd="commands with basefile replacement"
    # [ ] : Find all the dict.keys() usages that need list() for python3
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SeqSuite', '1.28.0', 'June 2024', '2014')
    description = 'Genomics and biological sequence analysis toolkit'
    author = 'Dr Richard J. Edwards.'
    comments = ['NOTE: Some tools are still in development and have not been published.',
                'Please see individual programs for citation details.',rje_obj.zen()]
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
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
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
mod = {'seqlist':rje_seqlist,'rje_seqlist':rje_seqlist,'rje_seq':rje_seq,'seq':rje_seq,'xref':rje_xref,'rje_xref':rje_xref,
       'rje_blast':rje_blast,'blast':rje_blast,'blast+':rje_blast,'buscomp':buscomp,
       'rje_zen':rje_zen,'zen':rje_zen,'zentest':rje_zen,'rje_mitab':rje_mitab,'mitab':rje_mitab,'pingu':pingu,
       'dbase':rje_dbase,'database':rje_dbase,'uniprot':rje_uniprot,'rje_uniprot':rje_uniprot,
       'rje_taxonomy':rje_taxonomy,'taxonomy':rje_taxonomy,'rje_tree':rje_tree,'tree':rje_tree,
       'paf':rje_paf,'rje_paf':rje_paf,'ppi':rje_ppi,'rje_ppi':rje_ppi,
       'pydocs':rje_pydocs,'genbank':rje_genbank,'gopher':gopher,'gablam':gablam,'pagsat':pagsat,'smrtscape':smrtscape,
       'gasp':gasp,'gff':rje_gff,'summarise':None,'apollo':rje_apollo,'rje_apollo':rje_apollo,'webapollo':rje_apollo,
       'samtools':rje_samtools,'snapper':snapper,'snp_mapper':snp_mapper,'seqmapper':seqmapper,'samphaser':samphaser,
       'saaga':saaga, 'diploidocus':diploidocus, 'synbad':synbad,'rje_genomics':rje_genomics,'genomics':rje_genomics,
       'ensembl':rje_ensembl,'rje_ensembl':rje_ensembl,'extatic':extatic,
       'depthsizer':depthsizer,'depthkopy':depthkopy,
       'fiesta':fiesta,'haqesac':haqesac,'multihaq':multihaq,'revert':revert}
#########################################################################################################################
# List of queries to remove from command list in REST server output
purgelist = ['mafft','clustalo','clustalw','maxforks','serverload','muscle','prog','debug','v','i','webserver','newlog',
             'ini']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SeqSuite Class                                                                                          #
#########################################################################################################################
class SeqSuite(rje_obj.RJE_Object):
    '''
    SeqSuite Class. Author: Rich Edwards (2014).

    Str:str
    - BatchArg=X          # Commandline argument to use for batchrun files ['seqin']
    - BatchLog=X          # Generate separate basefile.X log files for each batch run file (None for single log) [log]
    - Name = Identifies the tool to be used. Will load defaults from X.ini (after seqsuite.ini) [seqlist]

    Bool:boolean
    - BatchBase=T/F       # Whether to give each batch run a separate basefile=X command in place of log=X [True]
    - Help = Show the help documentation for program X. (Note that  [False]

    Int:integer

    Num:float
    
    List:list
    - BatchRun=FILELIST   # Batch run the program on selected files (wildcards allowed) []

    Dict:dictionary    

    Obj:RJE_Objects
    - Prog = Main program Object
    - ProgInfo = Info object for Program being run.
    '''
#########################################################################################################################
    def prog(self): ### Returns program name
        if self.obj['ProgInfo']: return self.obj['ProgInfo'].program
        else: return self.getStr('Name')
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BatchArg','BatchLog','Name']
        self.boollist = ['BatchBase','Help']
        self.intlist = []
        self.numlist = []
        self.listlist = ['BatchRun']
        self.dictlist = []
        self.objlist = ['Prog','ProgInfo']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'BatchArg':'seqin','BatchLog':'log','Name':'seqlist'})
        self.setBool({'BatchBase':True})
        if sys.argv[1:] and sys.argv[1] in mod: self.setStr({'Name':sys.argv[1]})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                ### Class Options (No need for arg if arg = att.lower()) ### 
                self._cmdRead(cmd,type='str',att='Name',arg='prog')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='Name',arg='program')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='Name',arg='tool')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['BatchArg','BatchLog'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BatchBase','Help'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['BatchRun']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Name') == 'summarise': return self.batchSummarise()
            if self.list['BatchRun']: return self.batchRun()
            if not self.setup(): return False
            slimobj = self.obj['Prog']
            self.obj['Info'] = self.log.obj['Info']
            slimobj.log.obj['Info'] = info = self.obj['ProgInfo']
            slimobj.log.info['Name'] = info.program
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Name') == 'haqesac': slimobj.run(setobjects=True)
            else: slimobj.run()
            self.printLog('#RUN','%s V%s run finished.' % (info.program,info.version))
            slimobj.log.obj['Info'] = self.obj['Info']
            slimobj.log.info['Name'] = self.obj['Info'].program
            return slimobj
        except SystemExit: raise
        except: self.errorLog('Problem during %s run.' % self.prog()); return False
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] ~ Setup Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['Prog'] = None
            prog = self.getStrLC('Name')
            if prog in mod:
                i = self.obj['ProgInfo'] = mod[prog].makeInfo()
                self.printLog('#PROG','%s V%s: %s' % (i.program,i.version,i.description))
                progcmd = rje.getCmdList([],info=i) + self.cmd_list + ['newlog=F']
                out = rje.Out(cmd_list=progcmd)
                out.printIntro(i)
                #self.debug(prog); self.debug(progcmd)
                if self.getBool('Help'): progcmd = mod[prog].cmdHelp(i,out,['help']+progcmd)
                self.printLog('#CMD','Full %s CmdList: %s' % (i.program,rje.argString(rje.tidyArgs(progcmd,nopath=self.getStrLC('Rest') and not self.dev(),purgelist=purgelist))),screen=False)
                #self.debug(prog); self.debug(progcmd)
            ## ~ [1a] ~ Make self.obj['Prog'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if prog in ['seqlist','rje_seqlist']: self.obj['Prog'] = rje_seqlist.SeqList(self.log,['reformat=rest'] + progcmd)
                elif prog in ['uniprot','rje_uniprot']: self.obj['Prog'] = rje_uniprot.UniProt(self.log,progcmd)
                elif prog in ['taxonomy','rje_taxonomy']: self.obj['Prog'] = rje_taxonomy.Taxonomy(self.log,progcmd)
                elif prog in ['tree','rje_tree']: self.obj['Prog'] = rje_tree.Tree(self.log,progcmd)
                elif prog in ['xref','rje_xref']: self.obj['Prog'] = rje_xref.XRef(self.log,progcmd)
                elif prog in ['seq','rje_seq']: self.obj['Prog'] = rje_seq.SeqList(self.log,progcmd)
                elif prog in ['mitab','rje_mitab']: self.obj['Prog'] = rje_mitab.MITAB(self.log,progcmd)
                elif prog in ['dbase','database']: self.obj['Prog'] = rje_dbase.DatabaseController(self.log,progcmd)
                elif prog in ['pydocs']: self.obj['Prog'] = rje_pydocs.PyDoc(self.log,progcmd)
                elif prog in ['ensembl','rje_ensembl']: self.obj['Prog'] = rje_ensembl.EnsEMBL(self.log,progcmd)
                elif prog in ['genbank','rje_genbank']: self.obj['Prog'] = rje_genbank.GenBank(self.log,progcmd)
                elif prog in ['extatic']: self.obj['Prog'] = extatic.ExTATIC(self.log,progcmd)
                elif prog in ['revert']: self.obj['Prog'] = revert.REVERT(self.log,progcmd)
                elif prog in ['fiesta']: self.obj['Prog'] = fiesta.FIESTA(self.log,progcmd)
                elif prog in ['gablam']: self.obj['Prog'] = gablam.GABLAM(self.log,['fullblast=T','keepblast=T']+progcmd)
                elif prog in ['gasp']: self.obj['Prog'] = gasp.GASP(self.log,progcmd)
                elif prog in ['gff']: self.obj['Prog'] = rje_gff.GFF(self.log,progcmd)
                elif prog in ['gopher']: self.obj['Prog'] = gopher.Gopher(self.log,progcmd)
                elif prog in ['haqesac']: self.obj['Prog'] = haqesac.HAQESAC(self.log,progcmd)
                elif prog in ['multihaq']: self.obj['Prog'] = multihaq.MultiHAQ(self.log,progcmd)
                elif prog in ['paf','rje_paf']: self.obj['Prog'] = rje_paf.PAF(self.log,progcmd)
                elif prog in ['ppi','rje_ppi']: self.obj['Prog'] = rje_ppi.PPI(self.log,progcmd)
                elif prog in ['pingu']: self.obj['Prog'] = pingu.PINGU(self.log,progcmd)
                #X# Disabled: elif prog in ['pacbio']: self.obj['Prog'] = rje_pacbio.PacBio(self.log,progcmd)
                elif prog in ['pagsat']: self.obj['Prog'] = pagsat.PAGSAT(self.log,progcmd)
                elif prog in ['smrtscape','pacbio']: self.obj['Prog'] = smrtscape.SMRTSCAPE(self.log,progcmd)
                elif prog in ['samtools']: self.obj['Prog'] = rje_samtools.SAMtools(self.log,progcmd)
                elif prog in ['samphaser']: self.obj['Prog'] = samphaser.SAMPhaser(self.log,progcmd)
                elif prog in ['seqmapper']: self.obj['Prog'] = seqmapper.SeqMapper(self.log,progcmd)
                elif prog in ['snapper']: self.obj['Prog'] = snapper.Snapper(self.log,progcmd)
                elif prog in ['snp_mapper']: self.obj['Prog'] = snp_mapper.SNPMap(self.log,progcmd)
                elif prog in ['rje_zen','zen','zentest']: self.obj['Prog'] = rje_zen.Zen(self.log,progcmd)
                elif prog in ['rje_blast','blast','blast+']: self.obj['Prog'] = rje_blast.BLASTRun(self.log,progcmd)
                elif prog in ['apollo','rje_apollo','webapollo']: self.obj['Prog'] = rje_apollo.Apollo(self.log,progcmd)
                elif prog in ['buscomp']: self.obj['Prog'] = buscomp.BUSCOMP(self.log,progcmd)
                elif prog in ['diploidocus']: self.obj['Prog'] = diploidocus.Diploidocus(self.log,progcmd)
                elif prog in ['saaga']: self.obj['Prog'] = saaga.SAAGA(self.log,progcmd)
                elif prog in ['synbad']: self.obj['Prog'] = synbad.SynBad(self.log,progcmd)
                elif prog in ['rje_genomics','genomics']: self.obj['Prog'] = rje_genomics.Genomics(self.log,progcmd)
                elif prog in ['depthsizer','gensize']: self.obj['Prog'] = depthsizer.DepthSizer(self.log,progcmd)
                elif prog in ['depthkopy','depthcopy']: self.obj['Prog'] = depthkopy.DepthKopy(self.log,progcmd)

            ### ~ [2] ~ Failure to recognise program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['Prog']:
                self.printLog('#ERR','Program "%s" not recognised.' % self.getStr('Name'))
                if self.i() < 0: return False
                if rje.yesNo('Show SeqSuite help with program options?'):
                    extracmd = cmdHelp(cmd_list=['help'])[1:]
                    if extracmd:
                        self.cmd_list += extracmd
                        self._cmdList()
                        if prog != self.getStrLC('Name'): return self.setup()
                self.setStr({'Name':rje.choice('Give program name (Blank or CTRL+C to quit)')})
                if self.getStrLC('Name'): return self.setup()
                else: return False
            return self.obj['Prog']     # Setup successful
        except KeyboardInterrupt: return False
        except SystemExit: raise
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### BatchRun Methods                                                                                        #
#########################################################################################################################
    def batchRun(self,returnobj=False):     ### Execute batch mode runs
        '''Execute batch mode runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            barg = self.getStrLC('BatchArg')
            if not barg: raise ValueError('Cannot use batchrun=FILELIST if batcharg=None.')
            batchfiles = self.list['BatchRun'][0:]
            self.list['BatchRun'] = []  # Avoid recursive running!
            blog = self.getStr('BatchLog')
            if not blog.startswith('.'): blog = '.%s' % blog
            if not blog.endswith('.log'): blog = '%s.log' % blog
            rawcmd = self.cmd_list[0:]
            rawlog = self.log
            batchobj = []
            ### ~ [1] Batch Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bx = 0
            for bfile in batchfiles:
                bx += 1
                self.printLog('#BATCH','Batch running %s of %s: %s=%s' % (rje.iStr(bx),rje.iLen(batchfiles),barg,bfile))
                ## Setup parameters
                bbase = rje.baseFile(bfile,strip_path=True)
                bcmd = ['%s=%s' % (barg,bfile)]
                if self.getBool('BatchBase'):
                    if blog == '.log': bcmd += ['basefile=%s' % bbase]
                    else: bcmd += ['basefile=%s%s' % (bbase,rje.baseFile(blog))]
                elif self.getStrLC('BatchLog'): bcmd += ['log=%s%s' % (bbase,blog)]
                else: bcmd += ['newlog=F']
                #self.debug(bcmd)
                ## Setup Seqsuite object
                self.cmd_list = rawcmd + bcmd
                self.log = rje.setLog(self.log.obj['Info'],self,self.cmd_list)                 # Sets up Log object for controlling log file output
                ## Run
                batchobj.append(self.run())
                ## Finish and Tidy
                self.log = rawlog
                runobj =  batchobj[-1]
                if runobj:
                    if not returnobj: batchobj[-1] = True
                    info = runobj.log.obj['Info']
                    self.printLog('#RUN','%s V%s run finished.' % (info.program,info.version))
                else: self.warnLog('Batch run failed (%s=%s).' % (barg,bfile))
            ### ~ [2] Finish and Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            failx = batchobj.count(False)
            self.printLog('#BATCH','%s batch runs complete: %s failed.' % (rje.iLen(batchfiles),rje.iStr(failx)))
            self.list['BatchRun'] = batchfiles
            return batchobj
        except: self.errorLog('%s.batchRun error' % self); return False
#########################################################################################################################
    def batchSummarise(self):   ### Batch run seqlist summarise on batchrun=LIST files and output table of results
        '''
        Batch run seqlist summarise on batchrun=LIST files and output table of results
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['BatchRun']: self._cmdReadList('batchrun=*.fasta','glist',['BatchRun'])
            if not self.list['BatchRun']: raise ValueError('Need to provide batchrun=LIST files for summarise mode.')
            db = rje_db.Database(self.log,self.cmd_list)
            if not db.baseFile(): db.baseFile('seqsuite')
            self.printLog('#BASE',db.baseFile())
            sdb = None
            if not self.force():
                sdb = db.addTable(mainkeys=['File'],name='summarise',expect=False)
            if not sdb: sdb = db.addEmptyTable('summarise',['File'],['File'])
            ### ~ [2] Run Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#BATCH','Batch summarising %s input files' % rje.iLen(self.list['BatchRun']))
            for file in self.list['BatchRun']:
                seqin = file
                if file.endswith('.gz'):
                    self.printLog('#UNZIP','unpigz -v {0}'.format(file))
                    os.popen('unpigz -v {0}'.format(file)).read()
                    seqin = file[:-3]
                seqdata = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % seqin,'autoload=T','summarise=F']).summarise()
                if file.endswith('.gz'):
                    self.printLog('#GZIP','pigz -v {0}'.format(seqin))
                    os.popen('pigz -v {0}'.format(seqin)).read()
                #!# For gapstats and fracstats, add retrieval and integration of SeqList.db() 'fracstats' and 'gaps' tables
                if seqdata:
                    if 'GC' in seqdata and 'GCPC' in seqdata:
                        seqdata.pop('GC')
                        seqdata['GCPC'] = '%.2f' % seqdata['GCPC']
                    if 'GapLength' in seqdata: seqdata['GapPC'] = '%.2f' % (100.0*seqdata['GapLength']/seqdata['TotLength'])
                    seqdata['MeanLength'] = '%.1f' % seqdata['MeanLength']
                    for field in rje.split('SeqNum, TotLength, MinLength, MaxLength, MeanLength, MedLength, N50Length, L50Count, NG50Length, LG50Count, N50Ctg, L50Ctg, GapLength, GapPC, GCPC',', '):
                        if field in seqdata and field not in sdb.fields(): sdb.addField(field)
                    for field in seqdata.keys():
                        if field not in sdb.fields(): sdb.addField(field)
                    sdb.addEntry(seqdata)
                else: self.errorLog('Summarise failed for %s' % file,printerror=False)
            ### ~ [3] Output Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.saveToFile()
            return True
        except: self.errorLog('%s.batchSummarise error' % self); return False
#########################################################################################################################
    ### <4> ### Testing Methods                                                                                         #
#########################################################################################################################
    def test(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except: self.errorLog('%s.test error' % self)
#########################################################################################################################
### End of SECTION II: SeqSuite Class                                                                                   #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: SeqSuite(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error:', sys.exc_info()[0])
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
