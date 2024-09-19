#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      BUDAPEST
Description:  Bioinformatics Utility for Data Analysis of Proteomics on ESTs
Version:      2.3
Last Edit:    31/07/13
Citation:     Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: 20924652]
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
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

Input:
    BUDAPEST takes three main files as input:
    * A MASCOT results file, specified by mascot=FILENAME.
    * The EST sequences used (or, at least, hit by) the MASCOT search, in fasta format, specified by seqin=FILENAME.
    * A protein database for BLAST-based annotation in fasta format, specified by searchdb=FILENAME.

Output:
    BUDAPEST produces the following main output files, where X is set by basefile=X:
    * X.budapest.tdt = main output table of results
    * X.budapest.fas = BLAST-annotated clustered consensus EST translations using FIESTA
    * X.summary.txt = summary of results from BUDAPEST pipeline
    * X.details.txt = full details of processing for each original MASCOT hit.

    Additional information can also be obtained from the additional sequence files:        
    * X.est.fas = subset of EST sequences from EST database that have 1+ hits in MASCOT results. 
    * X.translations.fas = fasta format of translated RF Hits that are retained after BUDAPEST cleanup.
    * X.fiesta.fas = BLAST-annotated consensus EST translations using FIESTA (pre min. peptide filtering)
    * X_HAQESAC/X.* = HAQESAC results files for annotating translated ESTs (haqesac=T only)
    * X_seqfiles/X.cluster*.fas = fasta files of translations and BLAST hits in NR clusters (clusterfas=T only)

    Lastly, reformatted MASCOT files are produced, named after the original input file (Y):
    * Y.mascot.txt = header information from the MASCOT file.
    * Y.mascot.csv = the delimited data portion of the MASCOT file.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    mascot=FILE     : Name of MASCOT csv file [None]
    seqin=FILE      : Name of EST fasta file used for search [None]
    searchdb=FILE   : Fasta file for GABLAM search of EST translations [None]
    partial=T/F     : Whether partial EST data is acceptable (True) or all MASCOT hits must be found (False) [True]
    itraq=T/F       : Whether data is from an iTRAQ experiment [False]
    empai=T/F       : Whether emPAI data is present in MASCOT file [True]
    samples=LIST    : List of X:Y, where X is an iTRAQ isotag and Y is a sample []
    ### ~ PROCESSING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minpolyat=X     : Min length of poly-A/T to be counted. (-1 = ignore all) [10]
    fwdonly=T/F     : Whether to treat EST/cDNA sequences as coding strands (False = search all 6RF) [False]
    minorf=X        : Min length of ORFs to be considered [10]
    topblast=X      : Report the top X BLAST results [10]
    minaln=X        : Min length of shared region for FIESTA consensus assembly [20]
    minid=X         : Min identity of shared region for FIESTA consensus assembly [95.0]
    minpep=X        : Minimum number of different peptides mapped to final translation/consensus [2]
    ### ~ SEQUENCE FORMATTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    newacc=X        : New base for sequence accession numbers ['BUD']
    gnspacc=T/F     : Convert sequences into gene_SPECIES__AccNum format wherever possible. [True]
    spcode=X        : Species code for EST sequences [None]
    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : "Base" name for all results files, e.g. X.budapest.tdt [MASCOT file basename]
    hitdata=LIST    : List of hit data to add to main budapest table [prot_mass,prot_pi]
    seqcluster=T/F  : Perform additional sequence (BLAST/GABLAM) clustering [True]
    clusterfas=T/F  : Generate fasta files of translations and BLAST hits in NR clusters [False]
    clustertree=LIST: List of formats for cluster tree output (3+ seqs only) [text,nsf,png]
    fiestacons=T/F  : Use FIESTA to auto-construct consensi from BUDAPEST RF translations [True]
    haqesac=T/F     : HAQESAC analysis of identified EST translations [True]
    blastcut=X      : Reduced the number of sequences in HAQESAC runs to X (0 = no reduction) [50]
    multihaq=T/F    : Whether to run HAQESAC in two-phases with second, manual phase [False]
    cleanhaq=T/F    : Delete excessive HAQESAC results files [True]
    haqdb=FILE      : Optional additional search database for MultiHAQ analysis [None]

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import fiesta, haqesac, rje, rje_db, rje_mascot, rje_menu, rje_seq, rje_seqlist, rje_sequence, rje_tree, rje_zen
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Reworked the pipeline in the light of discoveries made from version 0.0 runs.
    # 1.0 - Working version for basic analysis.
    # 1.1 - Modified to work with new MASCOT column headers.
    # 1.2 - Added tracking of MASCOT data, results tables and division of EST-RFs.
    # 1.3 - Split clustering into two levels: peptide and sequence clustering
    # 1.4 - Added FIESTA auto-construction of consensi from BUDAPEST RF translations [True]
    # 1.5 - Added MinPep filtering.
    # 1.6 - Improved tracking of peptides to final consensus sequences and output details.
    # 1.7 - Added menu and extra control of interactivity. Removed rfhits=F option.
    # 1.8 - Added preliminary iTRAQ handling.
    # 1.9 - Bug fixed for new MASCOT output.
    # 2.0 - Revised version using rje_mascot object for loading.
    # 2.1 - Improved handling of iTRAQ data using rje_mascot V1.2.
    # 2.2 - Removed unrequired rje_dismatrix import.
    # 2.3 - Updated to use rje_blast_V2. Needs further updates for BLAST+. Deleted obsolete OLDreadMascot() method.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Perform translation and BLAST *before* redundancy. Only keep one RF where possible as in FIESTA.
    # [Y] : Where peptides occur on multiple RFs, keep them all unless reduced by BLAST as above.
    # [?] : Add PFam HMM searches to convertESTtoRF()
    # [?] : Alter BLAST for EST to RFs to (a) use whole EST list at once, and (b) retain to re-use.
    # [Y] : Add output of topX BLAST hits to details file. (Make X an option.)
    # [ ] : Fix tree-drawing thing to make tree and then try outputs.
    # [ ] : Fix FastTree reading of trees.
    # [Y] : Map peptides onto final consensi and add minpep filter. (Number of different RF peptides assigned to consensus.)
    # [Y] : Cluster Consensi and rename? BUDXXXa BUDXXXb etc.
    # [ ] : Update to use BLAST+ and move/standardise the EST2RF and GABLAM clustering methods.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('BUDAPEST', '2.3', 'July 2013', '2008')
    description = 'Bioinformatics Utility for Data Analysis of Proteomics on ESTs'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show FIESTA commandline options?',default='N'): out.verbose(-1,4,text=fiesta.__doc__)
            if rje.yesNo('Show HAQESAC commandline options?',default='N'): out.verbose(-1,4,text=haqesac.__doc__)
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Budapest Class                                                                                               #
#########################################################################################################################
class Budapest(rje.RJE_Object):     
    '''
    Budapest Class. Author: Rich Edwards (2008).

    Info:str
    - Name = Name of input mascot file
    - HAQDB = Optional additional search database for MultiHAQ analysis [None]
    - SearchDB = Name of protein sequence file against which to search ESTs 
    
    Opt:boolean
    - CleanHAQ = Delete excessive HAQESAC results files [True]
    - ClusterFas = Generate fasta files of translations and BLAST hits in NR clusters [True]
    - emPAI = Whether emPAI data is present in MASCOT file [True]
    - FiestaCons = Use FIESTA to auto-construct consensi from BUDAPEST RF translations [True]
    - FwdOnly = Whether to treat EST/cDNA sequences as coding strands (False = search all 6RF) [False]
    - HAQESAC = HAQESAC analysis of identified EST translations [True]
    - iTRAQ = Whether data is from an iTRAQ experiment [False]
    - MultiHAQ = Whether to run HAQESAC in two-phases [True]
    - Partial = Whether partial EST data is acceptable (True) or all MASCOT hits must be found (False) [True]
    - RFHits = Whether to convert Hits to EST-RFs following "Bad" peptide removal [True]
    - SeqCluster = Peform additional sequence (BLAST/GABLAM) clustering [True]

    Stat:numeric
    - MinORF = Min length of ORFs to be considered [10]
    - MinPep = Minimum number of different peptides mapped to final translation/consensus [2]
    - MinPolyAT = Min length of ploy-A/T to be counted. (-1 = ignore all) [10]
    - TopBlast = Report the top X BLAST results [10]
    
    List:list
    - ClusterTree = List of formats for cluster tree output (3+ seqs only) [text,nsf,png]
    - HitData = List of hit data to add to main budapest table [prot_mass,prot_pi]

    Dict:dictionary
    - Details = For each Hit, a list of text details to be output to *.details.txt
    - Hits = Dictionary of protein hits from MASCOT file. #!# May replace with own Class at some point #!#
    - Peptides = Dictionary of {Translation shortName:[Peptide list]}
    - PepSeq = Dictionary of {peptide:[Sequence objects]}
    - PepTypes = Dictionary of {'Common':[peplist],'Cluster':[peplist],'Unique':[peplist]}
    - RF = Dictionary of {EST ID:[Acceptable RF translations]}
    - Support = Dictionary of {no. peptides:[List of IDs with this many peptides]}

    Obj:RJE_Objects
    - GABLAM = GABLAM object used for processing EST translations
    - SeqList = SeqList object containing input ESTs
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','SeqIn','SearchDB','HAQDB']
        self.optlist = ['OldPipe','RFHits','ClusterFas','SeqCluster','FiestaCons','HAQESAC','MultiHAQ','Partial',
                        'CleanHAQ','iTRAQ','emPAI','FwdOnly']
        self.statlist = ['TopBlast','MinPolyAT','MinORF','MinPep']
        self.listlist = ['HitData']
        self.dictlist = ['Details','Hits','RF','Support','Peptides']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setOpt({'RFHits':True,'ClusterFas':False,'SeqCluster':True,'FiestaCons':True,'HAQESAC':True,'MultiHAQ':False,'Partial':True,'CleanHAQ':True,'emPAI':True})
        self.setStat({'MinPolyAT':10,'MinORF':10,'TopBlast':10,'MinPep':2})
        self.list['HitData'] = ['prot_mass','prot_pi']
        self.list['ClusterTree'] = ['text','nsf','png']
        ### CmdList defaults ###
        for cmd in ['gnspacc=T']:
            if cmd not in self.cmd_list: self.cmd_list.insert(0,cmd)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdRead(cmd,type='info',att='Name',arg='mascot')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'opt',['OldPipe','RFHits','ClusterFas','SeqCluster','FiestaCons','HAQESAC',
                                             'MultiHAQ','Partial','CleanHAQ','iTRAQ','emPAI','FwdOnly'])
                self._cmdReadList(cmd,'int',['TopBlast','MinORF','MinPolyAT','MinPep'])
                self._cmdReadList(cmd,'file',['SeqIn','SearchDB','HAQDB'])
                self._cmdReadList(cmd.lower(),'list',['HitData','ClusterTree'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Run Methods                                                                                  #
#########################################################################################################################
    def run(self):  ### Main run method. Included menu for interactive setting of attributes
        '''Main run method. Included menu for interactive setting of attributes.'''
        ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['OldPipe']: return self.errorLog('OldPipe no longer supported - run V1.1 or earlier',printerror=False)
        if not self.opt['RFHits']: self.errorLog('RFHits=F no longer supported - run V1.6 or earlier')
        ### ~ [1] Main menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.i() < 0: return self.budapest()
        self.fileCheck()
        ## ~ [1a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        mainhead = 'Welcome to BUDAPEST %s' % self.log.obj['Info'].version
        #menu(callobj,headtext='',menulist=[],choicetext='Please select:',changecase=True,default='')
        #(edit code,description,optiontype,optionkey)
        mainmenu = [('','### ~ HELP ~ ###','',''),
                    ('B','BUDAPEST Help & Commands','showtext',__doc__),
                    ('F','FIESTA Help & Commands','showtext',fiesta.__doc__),
                    ('H','HAQESAC Help & Commands','showtext',haqesac.__doc__),
                    ('G','General Help & Commands','showtext',rje.__doc__),
                    ('','\n### ~ PARAMETERS ~ ###','',''),
                    #('E','Edit/Review Parameter settings','return','E'),
                    ('V','View commandline options','showtext','%s' % self.cmd_list),
                    ('A','Add commandline options','addcmd',''),
                    ('C','Check input files','return','C'),
                    ('','\n### ~ RUN ~ ###','',''),
                    ('Q','Quit BUDAPEST','return','Q'),
                    ('I','Interactive run','return','I'),
                    ('R','Run BUDAPEST (full auto)','return','R')]
        maintext = 'Select choice (or <ENTER> to continue)'
        ## ~ [1b] Menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        while True:
            choice = rje_menu.menu(self,mainhead,mainmenu,maintext,changecase=True,default='R')
            if choice == 'Q' and rje.yesNo('Quit BUDAPEST without running? Are you sure?'): return
            if choice == 'R': self.stat['Interactive'] = -1; self.cmd_list += ['i=-1']
            if choice in ['R','I'] and self.fileCheck(): return self.budapest()
            if choice == 'C': self.fileCheck()
#########################################################################################################################
    def fileCheck(self):    ### Checks input files etc.
        '''Checks input files etc.'''
        if os.path.exists(self.info['Name']): self.vPrint('MASCOT file "%s" found.' % self.info['Name'],v=0)
        else: self.vPrint('WARNING: MASCOT file "%s" does not exist!' % self.info['Name'],v=0)
        if os.path.exists(self.info['SeqIn']): self.vPrint('EST file "%s" found.' % self.info['SeqIn'],v=0)
        else: self.vPrint('WARNING: EST file "%s" does not exist!' % self.info['SeqIn'],v=0)
        if os.path.exists(self.info['SearchDB']): self.vPrint('SearchDB "%s" found.' % self.info['SearchDB'],v=0)
        else: self.vPrint('WARNING: SearchDB "%s" does not exist!' % self.info['SearchDB'],v=0)
        if os.path.exists(self.info['HAQDB']): self.vPrint('HAQDB "%s" found.' % self.info['HAQDB'],v=0)
        elif self.getStr('HAQDB').lower() not in ['','none']: self.vPrint('WARNING: HAQDB "%s" does not exist!' % self.info['HAQDB'],v=0)
        check = os.path.exists(self.info['Name']) and os.path.exists(self.info['SeqIn']) and os.path.exists(self.info['SearchDB'])
        if check: self.vPrint('\n*** All input files found. ***',v=0)
        return check
#########################################################################################################################
    def budapest(self): ### Main BUDAPEST pipeline
        '''
        Main run method.
        1. Read in MASCOT data.
        2. Extract all ESTs and translate into 6RFs. BLAST against search database and retain only RFs with BLAST hits
        (or all) and only peptides in good RFs.      
        3. Compile into protein hit entries and (child) peptide hit entries.
            3a. Perform shared-peptide comparisons to generate similarity matrix between hits
        4. Extract hits from EST database, assign frame and translate into sequence object
        5. Save hit sequences and perform all-by-all GABLAM redundancy check
        6. Perform GABLAM against metazoan sequence for possible annotation. Also PFam domain search etc.?
        '''
        try:### ~ [1] Setup & Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(self.info['SearchDB']): return self.log.errorLog('SearchDB "%s" does not exist!' % self.info['SearchDB'],printerror=False)
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = rje.baseFile(self.info['Name'])
            if not self.readMascot(): return

            ### ~ [2] Extract ESTs, translate, BLAST and identify good RFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            accfile = '%s.acc' % self.basefile()
            open(accfile,'w').write(rje.join(rje.sortKeys(self.dict['Hits']),'\n'))
            rje_seqlist.SeqList(self.log,self.cmd_list+['goodacc=%s.acc' % self.basefile(),'seqout=%s.hit_est.fas' % self.basefile(),'autoload','autofilter','seqmode=file'])
            self.info['SeqIn'] = '%s.hit_est.fas' % self.basefile()
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['accnr=F','seqnr=F','gnspacc=F','seqin=%s' % self.info['SeqIn']])
            if seqlist.seqNum() < len(self.dict['Hits']):
                return self.errorLog('Problem with EST file (seqin=FILE): not enough ESTs.',printerror=False)
            ## ~ [2a] ~ Match hits to ESTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['RF'] = {}
            estlist = []
            self.deBug(rje.sortKeys(self.dict['Hits']))
            for seq in seqlist.seq:
                if seq.info['Format'] == 'NCBI':
                    acc = rje.join(rje.split(seq.shortName(),'|')[:2],'|')
                    if acc not in self.dict['Hits']: acc = seq.info['AccNum']   # Alternative sequencing identifier
                else: acc = seq.info['AccNum']    # Can get gi number from seq.info['NCBI'] if seq.info['Format'] == 'NCBI'
                self.deBug(seq.info)
                self.deBug(acc in self.dict['Hits'])
                if acc in self.dict['Hits']: self.dict['Hits'][acc]['EST'] = seq
                elif seq.shortName() in self.dict['Hits']:
                    self.dict['Hits'][seq.shortName()]['EST'] = seq
                    self.dict['Hits'][acc] = self.dict['Hits'].pop(seq.shortName())
                else: continue
                estlist.append(seq)
                ## ~ [2b] ~ Translate and extract frames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.printLog('>>>','%d of %d <<<' % (len(estlist),len(self.dict['Hits'])),log=False)
                self.printLog('#EST','%s: %dnt (%d Peptides)' % (acc,seq.aaLen(),len(self.dict['Hits'][acc]['Peptides'])))
                self.dict['Details'][acc] = []
                try: self.convertESTtoRF(acc,seq)
                except KeyboardInterrupt:
                    if self.opt['Test']: break
                    else: raise
                except: raise
            RED = open('%s.summary.txt' % self.info['Basefile'],'w')
            RED.write('%d Hits read from %s => %s ESTs matched from %s' % (len(self.dict['Hits']),self.info['Name'],len(estlist),seqlist.info['Name']))
            RED.write('\n________________________________________\n\n')
            if estlist: seqlist.saveFasta(seqs=estlist,seqfile='%s.est.fas' % self.info['Basefile'])
            else: return self.errorLog('No EST sequences found in Hits!',printerror=False)
            ## ~ [2b] ~ Reduce hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if len(estlist) < len(self.dict['Hits']):
                if not self.opt['Partial']: return self.errorLog('Only %d EST sequences found for %d Hits! (partial=F)' % (len(estlist),len(self.dict['Hits'])),printerror=False)
                else:
                    missing = []
                    for acc in rje.sortKeys(self.dict['Hits']):
                        if 'PeptRF' not in self.dict['Hits'][acc]: self.dict['Hits'].pop(acc); missing.append(acc)
                    for acc in missing: self.printLog('#PARTIAL','EST %s missing from %s' % (acc,seqlist.info['Name']))
                    self.printLog('#HITS','%d EST hits not found in %s. (partial=T)' % (len(missing),seqlist.info['Name']))

            ### ~ [3] ~ Remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Regenerate Support Dictionary from reduced peptide list ~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Support'] = {}
            badlist = []
            for prot in rje.sortKeys(self.dict['Hits']):
                px = len(self.dict['Hits'][prot]['Peptides'])
                if not px: self.log.printLog('#REJ','Rejected %s: no peptides in best RFs' % prot)
                if px not in self.dict['Support']: self.dict['Support'][px] = []
                self.dict['Support'][px].append(prot)
                if self.dict['Hits'][prot]['Bad Peptides']: badlist.append(prot)
                #!# Add output of removed peptides #!#
            ## ~ [3b] Output removed peptides and proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            RED.write('Bad Peptides, removed as failing to match "good" RFs (%s Proteins):\n\n' % rje.integerString(len(badlist)))
            badlist.sort()
            for prot in badlist:
                RED.write('|--- %s: %s\n' % (prot,self.pepList(self.dict['Hits'][prot]['Bad Peptides'])))
            RED.write('\n________________________________________\n\n')
            if 0 in self.dict['Support']:
                RED.write('Bad Hits, removed as no peptides match "good" RFs (%s Proteins):\n\n' % rje.integerString(len(self.dict['Support'][0])))
                self.dict['Support'][0].sort()
                for prot in self.dict['Support'].pop(0):
                    RED.write('|--- %s: %s\n' % (prot,self.pepList(self.dict['Hits'][prot]['Bad Peptides'])))
                    self.dict['Details'][prot].append('REJECTED: No peptides match "good" RFs')
            else:
                #self.deBug(self.dict['Support'])
                RED.write('No Bad Hits removed with no peptides matching "good" RFs!\n\n')
            RED.write('\n________________________________________\n\n')            
            ## ~ [3c] Remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            px = len(self.dict['Hits'])
            self.convertHitsToRF()
            RED.write('%d EST Hits converted into %d RF Hits\n\n' % (px,len(self.dict['Hits'])))
            RED.write('\n________________________________________\n\n')            
            pxlist = rje.sortKeys(self.dict['Support'])
            pxlist.reverse()
            self.list['NR'] = []     # Non-redundant protein IDs
            for px in pxlist:
                pxnr = []           # Non-redundant proteins with px peptides
                pxp = self.dict['Support'][px][0:]
                while pxp:      # More proteins to process
                    nextprot = pxp.pop(0)
                    peps = self.dict['Hits'][nextprot]['Peptides']
                    ## ~ First check whether a longer match is a parent containing all peptides ~ ##
                    parent = False
                    for nr in self.list['NR']:      # Proteins with more proteins
                        parent = True
                        for pep in peps:
                            if pep not in self.dict['Hits'][nr]['Peptides']:    # Not a subset of this protein!
                                parent = False
                                break
                        if parent:
                            self.dict['Hits'][nr]['NRGrp'].append(nextprot)
                            break
                    if parent: continue
                    ## ~ Then find biggest protein with px peptides that are all the same ~ ##
                    (best,hitnum,eval,nrgrp) = (nextprot,self.dict['Hits'][nextprot]['HitNum'],self.dict['Hits'][nextprot]['E-Value'],[nextprot])
                    for pxp2 in pxp[0:]:    # Second protein to compare
                        if self.dict['Hits'][pxp2]['Peptides'] == peps:
                            pxp.remove(pxp2)
                            nrgrp.append(pxp2)
                            if self.dict['Hits'][nextprot]['E-Value'] < eval:   # Keep larger protein in preference
                                (best,hitnum,eval) = (nextprot,self.dict['Hits'][nextprot]['HitNum'],self.dict['Hits'][nextprot]['E-Value'])
                            elif self.dict['Hits'][nextprot]['HitNum'] > hitnum and self.dict['Hits'][nextprot]['E-Value'] == eval:
                                (best,hitnum,eval) = (nextprot,self.dict['Hits'][nextprot]['HitNum'],self.dict['Hits'][nextprot]['E-Value'])
                    ## ~ Update pxnr and continue ~ ##
                    pxnr.append(best)
                    self.dict['Hits'][best]['NRGrp'] = nrgrp[0:]
                self.list['NR'] += pxnr[0:]
            self.printLog('#NR','%s NR protein identifications from %s original hits' % (rje.integerString(len(self.list['NR'])),rje.integerString(len(self.dict['Hits']))))
            myrfs = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','autofilter=F'])
            myrfs.info['Name'] = '%s.translations.fas' % self.info['Basefile']
            myrfs.opt['ReplaceChar'] = False
            self.dict['NameToAcc'] = {}
            for nr in rje.sortKeys(self.dict['Hits']):  #!# No longer - self.list['NR']:
                for rf in rje.sortKeys(self.dict['RF'][nr]):
                    rfseq = self.dict['RF'][nr][rf]
                    myrfs.seq.append(rfseq); self.dict['NameToAcc'][rfseq.shortName()] = nr
                    ## ~ Add peptides as upper case ~ ##
                    oldseq = newseq = rfseq.info['Sequence'].lower()
                    for pep in self.dict['Hits'][nr]['Peptides']:
                        pi = oldseq.find(pep.lower())
                        if pi < 0: continue
                        newseq = newseq[:pi] + pep.upper() + newseq[pi+len(pep):]
                    rfseq.addSequence(newseq,case=True)
            myrfs.saveFasta(case=True)

            ## ~ [3d] Output redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Unique'] = []        # Whether the NR group has 100% unique (to group) peptides
            RED.write('Peptide Redundancy List (%s Groups):\n\n' % rje.integerString(len(self.list['NR'])))
            for nr in self.list['NR']:
                unique = True
                for pep in self.dict['Hits'][nr]['Peptides']:
                    for nr2 in self.list['NR']:
                        if nr != nr2 and pep in self.dict['Hits'][nr2]['Peptides']:
                            unique = False
                            break
                    if not unique: break
                if unique: self.list['Unique'].append(nr); RED.write('>%s [UNIQUE]\n' % nr)
                else: RED.write('>%s\n' % nr)
                for prot in self.dict['Hits'][nr]['NRGrp']:
                    RED.write('|--- %s: %s\n' % (prot,self.pepList(self.dict['Hits'][prot]['Peptides'])))
                    if prot != nr: self.dict['Details'][prot].append('REDUNDANT: All peptides found in %s' % nr)
                RED.write('\n')
            RED.write('________________________________________\n\n')
            self.printLog('#RED','%s Unique proteins (no shared peptides)' % rje.integerString(len(self.list['Unique'])))
            ## ~ [3f] Cluster Remaining Proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.clusterHits(myrfs)
            clusters = self.list['PepClusters'] 
            self.printLog('#RED','%s Shared Peptide Clusters' % rje.integerString(len(clusters)))
            RED.write('PepClusters:\n\n')
            for i in range(len(clusters)):
                RED.write('%s: %s\n' % (rje.preZero(i+1,len(clusters)),rje.join(clusters[i],', ')))
                for prot in clusters[i]:
                    RED.write('|--- %s: %s\n' % (prot,self.pepList(self.dict['Hits'][prot]['Peptides'])))
                RED.write('\n')
            RED.write('________________________________________\n\n')
            if self.opt['SeqCluster']: 
                clusters = self.list['SeqClusters'] 
                self.printLog('#RED','%s Sequence-based Clusters' % rje.integerString(len(clusters)))
                RED.write('SeqClusters:\n\n')
                for i in range(len(clusters)):
                    RED.write('%s: %s\n' % (rje.preZero(i+1,len(clusters)),rje.join(clusters[i],', ')))
                    for prot in clusters[i]:
                        RED.write('|--- %s: %s\n' % (prot,self.dict['Hits'][prot]['TopBlast']))
                    RED.write('\n')
                RED.write('________________________________________\n\n')

            ### ~ [4] ~ Output final matching RFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            RED.write('Translations:\n')
            for nr in self.list['NR']:
                for rf in rje.sortKeys(self.dict['RF'][nr]):
                    RED.write('|--- %s: %s\n' % (nr,self.dict['RF'][nr][rf].info['Name']))
            RED.write('________________________________________\n\n')

            ### ~ [5] ~ FIESTA Consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['FiestaCons']:
                ## ~ [5a] ~ Generate FIESTA consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.fiestaConsensi()
                consmap = rje.dataDict(self,'%s.bud.tdt' % self.info['Basefile'],['BUD'],'All',lists=True)
                conseq = self.obj['FIESTA'].obj['SeqList']; conseq.opt['ReplaceChar'] = False
                conseq.loadSeqs('%s.fiesta.fas' % self.info['Basefile'],seqtype='Protein')
                seqdict = conseq.seqNameDic('Max',proglog=True)
                RED.write('FIESTA Consensi:\n\n')
                for bud in rje.sortKeys(consmap):
                    RED.write('>%s %s\n' % (bud,seqdict[bud].info['Description']))
                    self.dict['Peptides'][bud] = []
                    for trans in consmap[bud]['Trans']:
                        RED.write('|--- %s [%s]\n' % (trans,self.pepList(self.dict['Peptides'][trans])))
                        for pep in self.dict['Peptides'].pop(trans):
                            if pep not in self.dict['Peptides'][bud]: self.dict['Peptides'][bud].append(pep)
                    RED.write('|=== %d different peptides: %s\n\n' % (len(self.dict['Peptides'][bud]),self.pepList(self.dict['Peptides'][bud])))
                RED.write('________________________________________\n\n')
                ## ~ [5b] ~ Filter by MinPep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pepfiltered = []
                if self.stat['MinPep'] > 1:
                    for rfname in rje.sortKeys(self.dict['Peptides']):
                        if len(self.dict['Peptides'][rfname]) < self.stat['MinPep']:
                            pepfiltered.append(rfname)
                            conseq.removeSeq(text='Excluded by MinPep filter (%d peptides)' % len(self.dict['Peptides'][rfname]),seq=seqdict[rfname])
                    RED.write('MinPep Filter: %d sequences < %d peptides removed\n' % (len(pepfiltered),self.stat['MinPep']))
                    for rfname in pepfiltered:
                        RED.write('|--- %s (%d) [%s]\n' % (rfname,len(self.dict['Peptides'][rfname]),self.pepList(self.dict['Peptides'][rfname])))
                    RED.write('________________________________________\n\n')
                else: RED.write('~ %s ~ \n\n' % self.printLog('#MINPEP','No MinPep filter'))
                ## ~ [5c] ~ Cluster Consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.clusterConsensi(RED)
                conseq.info['Name'] = '%s.budapest.fas' % self.info['Basefile']
                conseq.saveFasta(case=True)
                for delfile in ['bud.fas','bud.tdt','fiesta.tdt']:
                    if os.path.exists('%s.%s' % (self.info['Basefile'],delfile)): os.unlink('%s.%s' % (self.info['Basefile'],delfile))
                ## ~ [5d] Update details dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                myfiesta = self.obj['FIESTA']
                for trans in myfiesta.dict['Mapping']:
                    rfacc = '%s_RF%s' % (myfiesta.dict['Mapping'][trans]['EST'],myfiesta.dict['Mapping'][trans]['RF'])
                    self.dict['Details'][rfacc].append('FIESTA Consensus Sequence: %s' % myfiesta.dict['Mapping'][trans]['BUD'])
                    if myfiesta.dict['Mapping'][trans]['BUD'] in pepfiltered:
                        myfiesta.dict['Mapping'][trans]['Consensus'] = 'Removed by MinPep filter'
                        self.dict['Details'][rfacc].append('Removed by MinPep filter')
                    else:
                        try: self.dict['Details'][rfacc].append('Clustered FIESTA Consensus Sequence: %s %s' % (myfiesta.dict['Mapping'][trans]['Consensus'],myfiesta.dict['Mapping'][trans]['Description']))
                        except: self.dict['Details'][rfacc].append('Clustered FIESTA Consensus Failure!')
                    if rfacc not in self.dict['Hits']: self.errorLog('FIESTA RFAcc "%s" missing from Hit dictionary!' % rfacc,printerror=False); continue
                    for mkey in rje.sortKeys(myfiesta.dict['Mapping'][trans]): self.dict['Hits'][rfacc][mkey.lower()] = myfiesta.dict['Mapping'][trans][mkey]
            RED.close()

            ### ~ [6] ~ Final Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            self.outputDetails()    #!# Check FIESTA stuff #!#
            self.outputBUDAPEST()

            ### ~ [7] ~ Run Multiple rounds of HAQESAC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['FiestaCons']: self.multiHAQ(conseq,batchonly=False)
        except: self.errorLog('Error during main BUPAPEST pipeline.'); raise
#########################################################################################################################
    def pepList(self,peplist): return rje.join(peplist,', ')
#########################################################################################################################
    def convertESTtoRF(self,acc,seq,pepscreen=True,blastonly=False):   ### Converts a given EST sequence to RFs and BLASTs to find best sequence.
        '''
        Converts a given EST sequence to RFs and BLASTs to find best sequence. First translate, then BLAST against given
        search database. Read in alignments and map onto sequence. Identify ORFs with some coverage and reduce sequence
        to not include any ORFs without BLAST coverage. Give relevant sequence a new description using annotation and
        add sequences to self.dict['RF'].
        >> acc:str = Hit identifier corresonding to self.dict['Hit']
        >> seq:Sequence object = EST sequences as rje_sequence.Sequence object
        >> pepscreen:bool [True] = whether to screen according to peptide coverage
        >> blastonly:bool [False] = whether to only keep those RF with 1+ BLAST hit
        '''
        try:### ~ [1] ~ First translate, then BLAST against given search database. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: est = rje.split(acc,'|')[1]
            except: est = acc
            try: est = self.dict['Hits'][acc]['EST'].info['AccNum']
            except: pass
            #self.deBug('%s (%s) -> %s' % (seq.shortName(),acc,est))
            ## ~ [1a] ~ Generate sequence files of appropriate translated sequences ~~~~~~~~~~~~~~~ ##
            myrfs = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','autofilter=F','gnspacc=F'])
            myrfs.info['Name'] = '%s.%s.tmp.fas' % (est,rje.randomString(6))
            myrfs.opt['ReplaceChar'] = False
            #self.deBug(self.opt)
            rftrans = rje_sequence.estTranslation(seq.info['Sequence'],self.stat['MinPolyAT'],fwdonly=self.opt['FwdOnly'])
            #self.deBug(seq.info); self.deBug(rftrans)
            self.dict['Details'][acc].append('EST Translation: %d RF %s' % (len(rftrans),rje.sortKeys(rftrans)))
            self.dict['RF'][acc] = {}
            for rf in rje.sortKeys(rftrans):
                trans = rftrans[rf]
                if len(rje_sequence.bestORF(trans)) < self.stat['MinORF']:
                    self.dict['Details'][acc].append('RF%d removed: no ORF >= %daa' % (rf,self.stat['MinORF']))
                    self.log.printLog('#REM','Removed %s RF%s as no ORF >= %daa' % (acc,rf,self.stat['MinORF']))
                    continue     # Too crap!
                self.dict['RF'][acc][rf] = myrfs._addSeq('%s_RF%s' % (est,rf),trans)
            ## ~ [1b] ~ BLAST sequence against search database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            myrfs.saveFasta()
            #self.deBug(open(myrfs.info['Name'],'r').read())
            bcmd = ['blasti=%s' % myrfs.info['Name'],'blastd=%s' % self.info['SearchDB'],'blasto=%s.blast' % rje.baseFile(myrfs.info['Name'])]  #,'blastf=T']
            bcmd += ['blastb=%d' % self.stat['TopBlast'],'blastv=%d' % self.stat['TopBlast']]
            rfblast = rje_blast.blastObj(self.log,self.cmd_list+bcmd,type='Old')    #!# Need to convert to new BLAST
            rfblast.formatDB(force=False)
            rfblast.blast()
            os.unlink(myrfs.info['Name'])
            searchdbase = os.path.basename(rfblast.info['DBase'])
            ## ~ [1c] ~ Also search with PFam HMMs? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            
            ### ~ [2] ~ Read in alignments and map onto sequence. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rfblast.readBLAST(unlink=True)
            self.dict['Hits'][acc]['HitNum'] = hitx = rfblast.hitNum()         # Total number of BLAST hits
            self.printLog('#BLAST','%s: %d BLAST hits over %d RF' % (acc,hitx,len(self.dict['RF'][acc])))
            self.dict['Details'][acc].append('%d BLAST hits over %d RF' % (hitx,len(self.dict['RF'][acc])))
            self.dict['Hits'][acc]['E-Value'] = 1000
            self.dict['Hits'][acc]['TopBlast'] = {}
            self.dict['Hits'][acc]['BlastSeq'] = {}
            self.dict['Hits'][acc]['RFEvalue'] = {}
            searchdict = rfblast.searchSeq(myrfs,proglog=False,inverse=True)
            gindex = ['-','X','|','+']
            for rf in rje.sortKeys(self.dict['RF'][acc]):
                rfseq = self.dict['RF'][acc][rf]
                self.dict['Hits'][acc]['RFEvalue'][rf] = 1000
                ## ~ [2a] ~ Map BLAST GABLAM Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try: search = searchdict[rfseq]
                except:
                    self.errorLog('BLAST search for %s missing!' % rfseq.shortName(),printerror=False)
                    self.dict['Details'][acc].append('!!!RF%d BLAST error!!!' % rf)
                    #print rfseq.info['Name']
                    tester = rfblast.searchSeq(myrfs,proglog=False)
                    #for search in tester: print search.info['Name'], '>>>', tester[search].info['Name']
                    if self.opt['Test'] or self.opt['DeBug'] and rje.yesNo('Stop?'): break
                    continue
                if search.hitNum():
                    self.log.printLog('#RFHIT','%s RF%s: %d hits' % (acc,rf,search.hitNum()))
                    self.dict['Details'][acc].append('RF%s: %d hits' % (rf,search.hitNum()))
                elif hitx or blastonly:
                    myrfs.removeSeq(text='RF lacking BLAST hits.',seq=rfseq)
                    self.dict['Details'][acc].append('RF%d removed: lacking BLAST hits.' % rf)
                    self.dict['RF'][acc].pop(rf)
                    continue
                search.gablam(keepaln=True)  # GABLAM Alignments in hit.dict['GABLAM']['Aln']['Qry','Hit','QryO','HitO']
                galn = ['-'] * rfseq.aaLen()    # Replace with '|', '+' or 'X'
                if rfseq.info['Sequence'][-1] == '*': galn.pop(0)
                for hit in search.hit:
                    self.dict['Hits'][acc]['E-Value'] = min(hit.stat['E-Value'],self.dict['Hits'][acc]['E-Value'])
                    self.dict['Hits'][acc]['RFEvalue'][rf] = min(hit.stat['E-Value'],self.dict['Hits'][acc]['RFEvalue'][rf])
                    for i in range(len(galn)):
                        try:
                            if gindex.index(hit.dict['GABLAM']['Aln']['Qry'][i]) > gindex.index(galn[i]): galn[i] = hit.dict['GABLAM']['Aln']['Qry'][i]
                        except:
                            self.log.errorLog('Problem with %s vs %s GABLAM Aln position %d' % (search.info['Name'],hit.info['Name'],i+1))
                            #print rfseq.info['Sequence'], len(rfseq.info['Sequence'])
                            #print rje.join(galn,''), len(galn)
                            #print hit.dict['GABLAM']['Aln']['Qry'], len(hit.dict['GABLAM']['Aln']['Qry'])
                            self.deBug('')
                ## ~ [2b] ~ Also search with PFam HMMs? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rfseq.info['GAln'] = rje.join(galn,'')
            
            ### ~ [3] ~ Identify ORFs with some coverage and reduce sequences to covered ORFs ~~~~~~~~~~~~~~~~~~~~~~~ ###
                orfs = rje.split(rfseq.info['Sequence'],'*')
                (start,end) = (0,0)
                (firstcov,lastcov) = (-1,-1)
                for i in range(len(orfs)):
                    orf = orfs[i]
                    end += len(orf)
                    orfgaln = rfseq.info['GAln'][start:end]
                    if rje.count(orfgaln,'-') != len(orfgaln):   # At least partial coverage by 1+ hits
                        if firstcov < 0: firstcov = i
                        lastcov = i
                    start = end = (end + 1)
                lastcov += 1
                    
            ### ~ [4] ~ Give relevant sequence a new description using annotation and add to self.dict['RF'] ~~~~~~~~ ###
                ## ~ [4a] ~ Generate new reduced sequence based on ORF coverage ~~~~~~~~~~~~~~~~~~~ ##
                if hitx and firstcov < 0:   # Reject this RF as lacking BLAST hits
                    myrfs.removeSeq(text='RF lacking BLAST hits.',seq=rfseq)
                    self.dict['Details'][acc].append('RF%d removed: lacking BLAST hits.' % rf)
                    self.dict['RF'][acc].pop(rf)
                    continue
                elif hitx:
                    rfseq.info['Sequence'] = rje.join(orfs[firstcov:lastcov],'*')
                    self.dict['Details'][acc].append('RF%d BLAST hits cover ORFs %d-%d of %d' % (rf,firstcov+1,lastcov,len(orfs)))
                    if firstcov > 0 or lastcov < len(orfs):
                        self.dict['Details'][acc].append('RF%d sequence reduced to covered ORFs' % (rf))
                        self.log.printLog('#ORF','%s RF%d sequence reduced to BLAST-hit ORFs %d-%d of %d' % (acc,rf,firstcov+1,lastcov,len(orfs)))
                ## ~ [4b] ~ Add description to sequence based on BLAST hits (if any) ~~~~~~~~~~~~~~ ##
                if hitx:
                    bestseq = myrfs.seqFromFastaCmd(search.hit[0].info['Name'],rfblast.info['DBase'])
                    if bestseq: rfseq.info['Description'] = 'Similar (e=%s) to %s' % (rje.expectString(search.hit[0].stat['E-Value']),bestseq.info['Name'])
                    else: rfseq.info['Description'] = 'No BLAST hit (e<%s) to %s' % (rje.expectString(rfblast.stat['E-Value']),searchdbase)
                else: rfseq.info['Description'] = 'No BLAST hit (e<%s) to %s' % (rje.expectString(rfblast.stat['E-Value']),searchdbase)
                rfseq.info['Name'] = '%s %s' % (rfseq.info['Name'],rfseq.info['Description'])
                rfseq.info['Sequence'] = rfseq.info['Sequence'].lower()     # Peptides mapped onto upper case
                self.log.printLog('#RFSEQ',rfseq.info['Name'])
                self.dict['Details'][acc].append(rfseq.info['Name'])
                self.dict['Hits'][acc]['TopBlast'][rf] = []
                self.dict['Hits'][acc]['BlastSeq'][rf] = []
                for hit in search.hit[:self.stat['TopBlast']]:
                    x = search.hit.index(hit) + 1
                    hseq = myrfs.seqFromFastaCmd(hit.info['Name'],rfblast.info['DBase'])
                    self.dict['Hits'][acc]['BlastSeq'][rf].append(hseq)
                    self.dict['Hits'][acc]['TopBlast'][rf].append('RF%s %s [%s] %s' % (rf,rje.preZero(x,self.stat['TopBlast']),rje.expectString(hit.stat['E-Value']),hseq.info['Name']))
                
            ### ~ [5] ~ Map peptides and reduce accordingly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not pepscreen: return
            self.dict['Hits'][acc]['Bad Peptides'] = []
            self.dict['Hits'][acc]['Non-RF Peptides'] = []
            self.dict['Hits'][acc]['PeptRF'] = {}
            for pep in self.dict['Hits'][acc]['Peptides'][0:]:
                self.dict['Hits'][acc]['PeptRF'][pep] = []
                for rf in rje.sortKeys(self.dict['RF'][acc]):
                    rfseq = self.dict['RF'][acc][rf]
                    if rfseq.info['Sequence'].find(pep):
                        i = rfseq.info['Sequence'].lower().find(pep.lower())
                        if i < 0: continue
                        rfseq.info['Sequence'] = rfseq.info['Sequence'][:i] + pep + rfseq.info['Sequence'][i+len(pep):]
                        self.dict['Hits'][acc]['PeptRF'][pep].append(rf)
                        self.printLog('#PEPT','%s peptide %s => RF%d' % (acc,pep,rf))
                        self.dict['Details'][acc].append('Peptide %s found in RF%d' % (pep,rf))
                if not self.dict['Hits'][acc]['PeptRF'][pep]:
                    self.dict['Hits'][acc]['Peptides'].remove(pep)
                    self.dict['Hits'][acc]['Bad Peptides'].append(pep)
                    self.printLog('#PEPT','%s peptide %s => No RF' % (acc,pep))
                    self.dict['Details'][acc].append('Peptide %s found in No RFs' % (pep))

            ### ~ [6] ~ Get rid of RFs without any peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for rf in rje.sortKeys(self.dict['RF'][acc]):
                if self.dict['RF'][acc][rf].info['Sequence'].lower() == self.dict['RF'][acc][rf].info['Sequence']:
                    self.dict['RF'][acc].pop(rf)    # No peptides match this RF, so not interesting
                    self.printLog('#REM','Removed %s RF%s as no peptides match sequence' % (acc,rf))
                    self.dict['Details'][acc].append('Removed RF%s as no peptides match sequence.' % (rf))
        except KeyboardInterrupt:
            if rje.yesNo('Quit BUDAPEST?'): raise
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def convertHitsToRF(self):  ### Converts all EST hits to RF hits
        '''Converts all EST hits to RF hits.'''
        try:### ~ [1] ~ Convert Hit Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for acc in rje.sortKeys(self.dict['Hits']):
                if not self.dict['Hits'][acc]['Peptides']: continue
                convtxt = 'EST %s converted into %d RF:' % (acc,len(self.dict['RF'][acc]))
                for rf in rje.sortKeys(self.dict['RF'][acc]): convtxt += ' %s;' % rf
                self.dict['Details'][acc].append(convtxt)
                for rf in self.dict['RF'][acc]:
                    newacc = '%s_RF%s' % (acc,rf)
                    self.dict['Hits'][newacc] = rje.combineDict({'EST':acc},self.dict['Hits'][acc])
                    self.dict['Details'][newacc] = self.dict['Details'][acc][0:]
                    self.dict['Hits'][newacc]['E-Value'] = self.dict['Hits'][acc]['RFEvalue'][rf]
                    self.dict['Hits'][newacc]['BlastSeq'] = self.dict['Hits'][acc]['BlastSeq'][rf][0:]
                    self.dict['Hits'][newacc]['ESTTopBlast'] = rje.combineDict({},self.dict['Hits'][acc]['TopBlast'])
                    if self.dict['Hits'][acc]['TopBlast'][rf]:
                        self.dict['Hits'][newacc]['TopBlast'] = self.dict['Hits'][acc]['TopBlast'][rf][0]
                    else: self.dict['Hits'][newacc]['TopBlast'] = 'None'
                    for listkey in ['Peptides','Non-RF Peptides','Bad Peptides']:
                        self.dict['Hits'][newacc][listkey] = self.dict['Hits'][acc][listkey][0:]
                    for pep in self.dict['Hits'][newacc]['Peptides'][0:]:
                        if rf not in self.dict['Hits'][newacc]['PeptRF'][pep]:
                            self.dict['Hits'][newacc]['Non-RF Peptides'].append(pep)
                            self.dict['Hits'][newacc]['Peptides'].remove(pep)
                            self.dict['Details'][newacc].append('Peptide %s found in different RF(s)' % (pep))
                    self.dict['RF'][newacc] = {rf:self.dict['RF'][acc][rf]}
                    self.dict['Details'][newacc].append('%d of %d peptides found in RF%s' % (len(self.dict['Hits'][newacc]['Peptides']),len(self.dict['Hits'][acc]['Peptides']),rf))
                self.printLog('#CONV',convtxt)
                self.dict['Hits'].pop(acc)

            ### ~ [2] ~ Remake Support Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###   
            self.dict['Support'] = {}
            badlist = []
            for prot in rje.sortKeys(self.dict['Hits']):
                px = len(self.dict['Hits'][prot]['Peptides'])
                if not px: continue
                if px not in self.dict['Support']: self.dict['Support'][px] = []
                self.dict['Support'][px].append(prot)
            self.printLog('#RFHIT','EST hits converted to RF hits')                
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def clusterHits(self,rfseq=None):     ### Performs peptide and optional sequence clustering on Hits
        '''
        Performs peptide and optional sequence clustering on Hits.
        >> rfseq:SeqList object containing non-redundant translations 
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tocluster = self.list['NR'][0:]                 # List of proteins to cluster
            for prot in tocluster[0:]:
                self.dict['Hits'][prot]['PepCluster'] = {}  # Dictionary of {prot2:[no. shared peptides]}
                self.dict['Hits'][prot]['SeqCluster'] = {}  # Dictionary of {prot2:[shared BLAST Hits/GABLAM]}
                if prot in self.list['Unique']: tocluster.remove(prot)
            ### ~ [1] ~ Cluster using Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Identify shared Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for prot in tocluster[:-1]:
                for prot2 in tocluster[tocluster.index(prot)+1:]:
                    shared = []
                    for pep in self.dict['Hits'][prot]['Peptides']:     #!# Replace with Python sets #!#
                        if pep in self.dict['Hits'][prot2]['Peptides']: shared.append(pep)
                    if shared:
                        self.dict['Hits'][prot]['PepCluster'][prot2] = self.dict['Hits'][prot2]['PepCluster'][prot] = shared[0:]
                        self.dict['Details'][prot].append('Shares %d common peptides with %s: %s' % (len(shared),prot2,self.pepList(shared)))
                        self.dict['Details'][prot2].append('Shares %d common peptides with %s: %s' % (len(shared),prot,self.pepList(shared)))
            ## ~ [1b] ~ Perform clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            clusters = []
            while tocluster:
                prot = tocluster.pop(0)
                pcluster =  self.dict['Hits'][prot]['PepCluster'].keys()
                cloop = True
                while cloop:
                    cloop = False
                    for prot2 in pcluster:
                        if prot2 in tocluster:
                            tocluster.remove(prot2)
                            cloop = True
                            for prot3 in self.dict['Hits'][prot2]['PepCluster'].keys():
                                if prot3 not in pcluster: pcluster.append(prot3)
                if prot not in pcluster: pcluster.append(prot)
                pcluster.sort()
                clusters.append(pcluster)
            self.list['PepClusters'] = clusters
            ### ~ [2] ~ Cluster using Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['SeqCluster']: return
            ## ~ [2a] ~ Clustering on Shared BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tocluster = self.list['NR'][0:]                 # List of proteins to cluster
            for prot in tocluster[:-1]:
                prot_hits = []
                for seq in self.dict['Hits'][prot]['BlastSeq']: prot_hits.append(seq.shortName())
                for prot2 in tocluster[tocluster.index(prot)+1:]:
                    shared = []
                    for seq in self.dict['Hits'][prot2]['BlastSeq']:
                        if seq.shortName() in prot_hits: shared.append(seq.shortName())
                    if shared:
                        shared.sort()
                        self.dict['Hits'][prot]['SeqCluster'][prot2] = shared[0:]
                        self.dict['Hits'][prot2]['SeqCluster'][prot] = shared[0:]
                        self.dict['Details'][prot].append('Shares %d common BLAST hits with %s: %s' % (len(shared),prot2,self.pepList(shared)))
                        self.dict['Details'][prot2].append('Shares %d common BLAST hits with %s: %s' % (len(shared),prot,self.pepList(shared)))
            ## ~ [2b] ~ GABLAM of translations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                #!# Replace with actual GABLAM clustering and reading in (or conversion) of results. UPC clustering? #!#
                seqfile = rfseq.info['Name']        # File containing NR translations
                seqdict = rfseq.seqNameDic()
                ### Setup BLAST ###
                blast = rje_blast.blastObj(self.log,['blastf=T','blaste=1e-4']+self.cmd_list,type='Old')
                blast.setInfo({'InFile':seqfile,'DBase':seqfile,'Name':'%s.self.blast' % self.info['Basefile']})
                blast.setStat({'OneLine':rfseq.seqNum(),'HitAln':rfseq.seqNum()})
                blast.formatDB(fasfile=seqfile,force=self.opt['Force'],protein=True)
                if not rje_blast.checkForDB(dbfile=seqfile,checkage=False,log=self.log,protein=True,oldblast=blast.getBool('OldBLAST')):
                    self.log.errorLog('FormatDB failed for unknown reasons.',printerror=False)
                    raise IOError
                ### Perform and read BLAST ###
                blast.blast(cleandb=True,use_existing=False)
                if not blast.readBLAST(gablam=True,unlink=True):
                    self.log.errorLog('Major problem with BLAST for unknown reasons.')
                    raise IOError
                ### Add to SeqClusters ##
                for search in blast.search:
                    prot = self.dict['NameToAcc'][search.info['Name']]
                    if 'SeqCluster' not in self.dict['Hits'][prot]:
                        if prot in self.list['NR']: self.printLog('#ERR','SeqCluster missing from %s hits dictionary!' % prot,screen=False)
                        self.dict['Hits'][prot]['SeqCluster']= {}
                    for hit in search.hit:
                        prot2 = self.dict['NameToAcc'][hit.info['Name']]
                        if 'SeqCluster' not in self.dict['Hits'][prot2]:
                            if prot2 in self.list['NR']: self.printLog('#ERR','SeqCluster missing from %s hits dictionary!' % prot2,screen=False)
                            self.dict['Hits'][prot2]['SeqCluster']= {}
                        if prot != prot2:
                            if prot2 not in self.dict['Hits'][prot]['SeqCluster']: self.dict['Hits'][prot]['SeqCluster'][prot2] = []
                            if prot not in self.dict['Hits'][prot2]['SeqCluster']: self.dict['Hits'][prot2]['SeqCluster'][prot] = []
                            gdis = 100.0 * max( (hit.dict['GABLAM']['Query']['GABLAMO ID'] / float(search.stat['Length'])), (hit.dict['GABLAM']['Hit']['GABLAMO ID'] / float(hit.stat['Length'])) )
                            self.dict['Hits'][prot]['SeqCluster'][prot2].append('GABLAM: %.2f%% ID' % gdis)
                            self.dict['Hits'][prot2]['SeqCluster'][prot].append('GABLAM: %.2f%% ID' % gdis)
            except: self.errorLog('GABLAM clustering error')
            ## ~ [2c] ~ Perform clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            clusters = []
            while tocluster:
                prot = tocluster.pop(0)
                pcluster =  self.dict['Hits'][prot]['SeqCluster'].keys()
                cloop = True
                while cloop:
                    cloop = False
                    for prot2 in pcluster:
                        if prot2 in tocluster:
                            tocluster.remove(prot2)
                            cloop = True
                            for prot3 in self.dict['Hits'][prot2]['SeqCluster'].keys():
                                if prot3 not in pcluster: pcluster.append(prot3)
                if prot not in pcluster: pcluster.append(prot)
                pcluster.sort()
                clusters.append(pcluster)
            self.list['SeqClusters'] = clusters
            ## ~ [2d] ~ Make new SeqCluster entry for proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for acc in rje.sortKeys(self.dict['Hits']):
                cx = 0
                if acc in self.list['NR']: nr = acc
                elif not self.dict['Hits'][acc]['Peptides']: nr = None; cx = -1
                else: nr = None
                if cx >= 0 and not nr:
                    for nracc in self.list['NR']:
                        if acc in self.dict['Hits'][nracc]['NRGrp']: nr = nracc
                if nr:
                    for i in range(len(clusters)):
                        if nr in clusters[i]: cx = i + 1; break
                if cx >= 0: self.dict['Hits'][acc]['SeqCluster'] = rje.preZero(cx,len(clusters))
                else: self.dict['Hits'][acc]['SeqCluster'] = '-1'
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def outputDetails(self):    ### Outputs contents of self.dict['Details'] into file
        '''Outputs contents of self.dict['Details'] into file.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            detfile = '%s.details.txt' % self.info['Basefile']
            rje.backup(self,detfile)
            DET = open(detfile,'w')
            DET.write('Analysis details for BUDAPEST analysis of %s @ %s\n\n' % (self.info['Name'],time.asctime(time.localtime(time.time()))))
            ### ~ [1] Output Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for acc in rje.sortKeys(self.dict['Hits']):
                try:
                    data = self.dict['Hits'][acc]
                    (px,bx) = (len(data['Peptides']), len(data['Bad Peptides']))
                    if acc in self.list['Unique']: type = 'UNIQUE'
                    elif acc in self.list['NR']: type = 'NON-REDUNDANT'
                    elif not px: type = 'REJECTED'
                    else: type = 'REDUNDANT'
                    DET.write('>%s: %s; %d peptides (%d Good; %d Bad).\n' % (acc,type,(px+bx),px,bx))
                    if px: DET.write('| > Good: %s\n' % self.pepList(data['Peptides']))
                    if bx: DET.write('| > Bad: %s\n' % self.pepList(data['Bad Peptides']))
                    DET.write('| > Top BLAST hits:\n')
                    try:
                        for rf in self.dict['RF'][acc]:
                            for line in self.dict['Hits'][acc]['ESTTopBlast'][rf]: DET.write('|-- %s\n' % line)
                    except: self.errorLog(rje_zen.Zen().wisdom())
                    DET.write('| > BUDAPEST Details:\n')
                    for line in self.dict['Details'][acc]: DET.write('|-- %s\n' % line)
                    if self.opt['FiestaCons']:
                        if acc in self.obj['FIESTA'].dict['Mapping']: DET.write('| >> FIESTA consensus %s\n' % self.obj['FIESTA'].dict['Mapping'][acc]['BUD'])
                    DET.write('\n\n')
                except: DET.write('| !!!ERROR!!!\n\n')
                
            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            DET.write('End of Analysis. Please e-mail r.edwards@southampton.ac.uk with queries.')
            DET.close()
            self.log.printLog('#DETAILS','Details output to %s' % detfile)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def outputBUDAPEST(self):   ### Outputs main BUDAPEST summary table
        '''Outputs main BUDAPEST summary table.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.budapest.tdt' % self.info['Basefile']
            headers = ['hit','type','pepcluster'] + self.list['HitData'] + ['pepseq'] + ['bad_pep','pep_num']
            headers += ['nonrf_pep','rf_pep','pepscore']
            fhead = ['est','rf','trans','bud','consensus','description','conspep','uniqpep','cluspep','commpep','exact','coverage']
            if self.opt['SeqCluster']: headers.insert(3,'seqcluster')
            #else: return self.printLog('#OUT','BUDAPEST Table output not supported if rfhits=F')# headers += ['rf']
            if self.opt['FiestaCons']: headers = headers[:1] + fhead[:3] + headers[1:] + fhead[3:]
            else: headers += ['topblast','evalue','species','id','accnum']
            rje.delimitedFileOutput(self,outfile,headers,rje_backup=True)
            clusters = {}       # Store clusters for ClusterFas output
            ### ~ [1] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            for acc in rje.sortKeys(self.dict['Hits']):
                try:
                    ## ~ [1a] ~ Basic data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    data = {'hit':acc}
                    for h in self.list['HitData']:
                        try: data[h] = self.dict['Hits'][acc][h]
                        except: self.printLog('#ERR','Data "%s" missing for %s' % (h,acc),screen=False)
                    ## ~ [1b] ~ Hit Type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if acc in self.list['Unique']: data['type'] = 'UNIQUE'
                    elif acc in self.list['NR']: data['type'] = 'NON-REDUNDANT'
                    elif not self.dict['Hits'][acc]['Peptides']: data['type'] = 'REJECTED'
                    else: data['type'] = 'REDUNDANT'
                    ## ~ [1b] ~ Cluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if self.opt['SeqCluster']: data['seqcluster'] = self.dict['Hits'][acc]['SeqCluster']
                    data['pepcluster'] = 0
                    if acc in self.list['Unique']: data['pepcluster'] = self.list['Unique'].index(acc) + 1
                    elif acc in self.list['NR']: nr = acc
                    elif not self.dict['Hits'][acc]['Peptides']: data['pepcluster'] = -1
                    else: nr = None
                    if not data['pepcluster'] and not nr:
                        for nracc in self.list['NR']:
                            if acc in self.dict['Hits'][nracc]['NRGrp']: nr = nracc
                    if not data['pepcluster'] and nr:
                        if nr in self.list['Unique']: data['pepcluster'] = self.list['Unique'].index(nr) + 1
                        else:
                            for i in range(len(self.list['PepClusters'])):
                                if nr in self.list['PepClusters'][i]:
                                    data['pepcluster'] = len(self.list['Unique']) + i + 1
                                    break
                    if data['pepcluster'] not in clusters: clusters[data['pepcluster']] = []
                    clusters[data['pepcluster']].append(acc)
                    data['pepcluster'] = rje.preZero(data['pepcluster'],len(self.list['Unique'])+len(self.list['NR']))
                    ## ~ [1c] ~ Peptides/RF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    data['bad_pep'] = len(self.dict['Hits'][acc]['Bad Peptides'])
                    data['nonrf_pep'] = len(self.dict['Hits'][acc]['Non-RF Peptides'])
                    data['rf_pep'] = len(self.dict['Hits'][acc]['Peptides'])
                    data['pepseq'] = rje.join(self.dict['Hits'][acc]['Peptides'],',')
                    data['pep_num'] = sum([data['nonrf_pep'],data['rf_pep']])
                    if data['pep_num']: data['pepscore'] = '%.2f' % (data['rf_pep'] * data['rf_pep'] / float(data['pep_num']))
                    else: data['pepscore'] = 0.0
                    if self.opt['FiestaCons']:
                        for h in fhead:
                            try: data[h] = self.dict['Hits'][acc][h]
                            except: pass
                    ## ~ [1d] ~ BLAST Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try: data['topblast'] = rje.join(rje.split(self.dict['Hits'][acc]['TopBlast'])[4:])
                    except: data['topblast'] = 'None'
                    if not data['topblast']: data['topblast'] = self.dict['Hits'][acc]['TopBlast']
                    data['evalue'] = self.dict['Hits'][acc]['E-Value']
                    if data['topblast'] == 'None': data['evalue'] = '-'
                    try:
                        data['species'] = rje.matchExp('Tax=(\S.+)\s+RepID',data['topblast'])[0]
                        data['topblast'] = rje.matchExp('^(.+)\s+Tax=',data['topblast'])[0]
                    except: data['species'] = '-'
                    try: data['id'] = rje.split(self.dict['Hits'][acc]['TopBlast'])[3]
                    except: data['id'] = '-'
                    data['accnum'] = rje.split(data['id'],'__')[-1]
                    data['id'] = rje.split(data['id'],'__')[0]
                    if data['species'] == '-' and len(rje.split(data['id'],'_')) > 1: data['species'] = rje.split(data['id'],'_')[1]
                    rje.delimitedFileOutput(self,outfile,headers,datadict=data)
                except: self.errorLog('Problem with BUDAPEST table output for %s' % acc)
            self.printLog('#OUT','Main summary table output to %s' % outfile)
            if self.opt['ClusterFas'] and self.opt['SeqCluster']: self.clusterFas()    #clusters)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def clusterFas(self):  ### Aligns and Outputs EST translations and BLAST hits to fasta files
        '''
        Aligns and Outputs EST translations and BLAST hits to fasta files.
        >> clusters:dict of {cluster ID:[acclist]}
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            clusters = {}
            for acc in rje.sortKeys(self.dict['Hits']):
                cx = self.dict['Hits'][acc]['SeqCluster']
                if cx not in clusters: clusters[cx] = []
                clusters[cx].append(acc)
            if not clusters: return self.printLog('#FAS','No clusters for ClusterFas output')
            cdir = rje.makePath('%s_seqfiles' % self.info['Basefile'])
            rje.mkDir(self,cdir)
            cbase = rje.makePath(cdir + self.info['Basefile'],wholepath=True)
            cseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','usecase=T','gnspacc=F'])
            coretreecmd = ['autoload=T','outnames=long','root=mid','maketree=clustalw','bootstraps=100']
            ### ~ [1] ~ Rejects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if '-1' in clusters:
                cseq.seq = []
                cseq.info['Name'] = '%s.rejects.fas' % cbase
                rejects = clusters.pop('-1')
                rejects.sort()
                for acc in rejects:
                    for rf in rje.sortKeys(self.dict['RF'][acc]): cseq.seq.append(self.dict['RF'][acc][rf])
                cseq.saveFasta()
            ### ~ [2] ~ Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cx in rje.sortKeys(clusters):
                cseq.seq = []
                cseq.info['Name'] = '%s.%s.fas' % (cbase,cx)
                clusters[cx].sort()
                for acc in clusters[cx]:
                    for rf in rje.sortKeys(self.dict['RF'][acc]): cseq.seq.append(self.dict['RF'][acc][rf])
                blastseqs = []
                for acc in clusters[cx]:
                    for seq in self.dict['Hits'][acc]['BlastSeq']:
                        if seq.shortName() not in blastseqs: cseq.seq.append(seq); blastseqs.append(seq.shortName())
                #self.deBug(cseq.info)
                if cseq.seqNum() > 1: cseq.align(outfile=None,mapseq=True)
                cseq.saveFasta()
                if cseq.seqNum() > 2:
                    try:#!# Needs improving #!#
                        treecmd = ['seqin=%s' % cseq.info['Name'],'treeformats=%s' % rje.join(self.list['ClusterTree'],','),'savetree=%s.%s.nsf' % (cbase,cx)]
                        tree = rje_tree.Tree(self.log,coretreecmd+self.cmd_list+treecmd)
                        #self.deBug(tree.info)
                    except: self.errorLog('%s tree error' % type)
                    if os.path.exists('%s.phb' % rje.baseFile(cseq.info['Name'],strip_path=True)):
                        os.unlink('%s.phb' % rje.baseFile(cseq.info['Name'],strip_path=True))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def fiestaConsensi(self,peptides=True):   ### Generates and annotates consensus EST translations using FIESTA
        '''
        Aligns and Outputs EST translations and BLAST hits to fasta files.
        >> clusters:dict of {cluster ID:[acclist]}
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fcmd = ['gnspacc=T','newacc=BUD','pickup=F']+self.cmd_list+['basefile=%s' % self.info['Basefile'],'est2rf=F','blastann=F','dna=F','replacechar=F']
            self.deBug(fcmd)
            for cmd in fcmd:
                try: x = len(rje.split(cmd,'='))
                except: self.deBug(cmd)                             
            myfiesta = fiesta.FIESTA(self.log,fcmd[0:])
            myfiesta.obj['BLAST'] = rje_blast.blastObj(self.log,fcmd[0:],type='Old')
            myfiesta.obj['BLAST'].stat['Verbose'] = self.stat['Verbose'] - 1
            myfiesta.list['Headers'] = ['EST','RF','Trans']
            myfiesta.stat['MinORF'] = self.stat['MinORF']
            cseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','usecase=T'])
            cseq.seq = []
            for acc in rje.sortKeys(self.dict['Hits']):  #!# No longer - self.list['NR']:
                for rf in rje.sortKeys(self.dict['RF'][acc]):
                    seq = self.dict['RF'][acc][rf]
                    rfname = seq.shortName()
                    cseq.seq.append(seq)
                    if rje.matchExp('(\S+)_RF',acc): est = rje.matchExp('(\S+)_RF',acc)[0]
                    else: est = acc
                    myfiesta.dict['Mapping'][rfname] = {'EST':est,'Trans':rfname,'RF':rf}
                    if peptides:
                        if rfname not in self.dict['Peptides']: self.dict['Peptides'][rfname] = []
                        for pep in self.dict['Hits'][acc]['Peptides']:
                            if pep not in self.dict['Peptides'][rfname]: self.dict['Peptides'][rfname].append(pep)
                    #!# Consider converting Hits dictionary to this? #!#
            #myfiesta.(myfiesta.dict['Mapping'])
            myfiesta.setupSeqNames()
            ### ~ [1] ~ Generate consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug(cseq.seqNum())
            myfiesta.makeConsensi(cseq,log='BUDAPEST',type='BUD')
            #self.deBug(cseq.seqNum())
            cseq.info['Name'] = '%s.bud.fas' % self.info['Basefile']
            cseq.saveFasta(append=False,log=True)
            myfiesta.saveMapping('BUD')
            myfiesta.annotateTrans(cseq,type='fiesta')
            myfiesta.obj['SeqList'] = cseq
            self.obj['FIESTA'] = myfiesta
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def clusterConsensi(self,RED):  ### Clusters and renames sequences in FIESTA seqlist
        '''Clusters and renames sequences in FIESTA seqlist.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            myfiesta = self.obj['FIESTA']
            seqlist = myfiesta.obj['SeqList']
            if not seqlist.seqNum():
                RED.write('No sequences remain for consensus clustering')
                RED.write('\n________________________________________\n\n')
                return self.printLog('#NULL','No sequences remain for consensus clustering')
            seqfile = '%s.tmpdb' % self.info['Basefile']
            seqlist.saveFasta(seqfile=seqfile)
            seqdict = seqlist.seqNameDic()
            ## ~ [1a] Setup BLAST etc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = rje_blast.blastObj(self.log,['blastf=T','blaste=1e-4']+self.cmd_list+['dna=F'],type='Old')
            blast.setInfo({'InFile':seqfile,'DBase':seqfile,'Name':'%s.tmp.blast' % self.info['Basefile'],'Type':'blastp'})
            blast.setStat({'OneLine':seqlist.seqNum(),'HitAln':0})
            blast.formatDB(fasfile=seqfile,force=True,protein=True)

            ### ~ [2] Perform BLAST and generate hit matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.blast(cleandb=False,use_existing=False,log=True)
            blast.readBLAST(gablam=False,unlink=True,log=True)
            rje_blast.cleanupDB(self,seqfile,deletesource=True)
            ## ~ [2a] Cluster by BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cluster = {}    # Dictionary of {seq:hit seqs} for clustering
            for search in blast.search:
                seq = seqdict[search.info['Name']]
                cluster[seq] = []
                for hit in search.hit: cluster[seq].append(seqdict[hit.info['Name']])
            #self.deBug(cluster)
            ## ~ [2b] Combine clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            clusters = []   # List of [seqs] in clusters
            seqclusters = {}    # Dictionary of {seq:[cluster]}
            for seq in seqlist.seqs():
                if seq not in cluster: continue
                newcluster = [seq]; seqclusters[seq] = newcluster
                hits = cluster.pop(seq)
                while hits:
                    hit = hits.pop(0)
                    if hit not in newcluster: newcluster.append(hit)
                    if hit in seqclusters and seqclusters[hit] != newcluster:
                        prevcluster = seqclusters[hit]
                        for newseq in newcluster:
                            if newseq not in prevcluster: prevcluster.append(newseq)
                        newcluster = prevcluster
                        for seq in newcluster: seqclusters[seq] = newcluster
                    else: seqclusters[hit] = newcluster
                    if hit in cluster: hits += cluster.pop(hit)
                if newcluster not in clusters: clusters.append(newcluster)
            #self.deBug(clusters)

            ### ~ [3] Assign peptides to consensi as "Common", "Cluster" or "Unique" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Match peptides to sequence lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pepcons = {}
            for seq in seqlist.seqs():
                rfname = seq.info['AccNum']
                for pep in self.dict['Peptides'][rfname]:
                    if pep not in pepcons: pepcons[pep] = []
                    pepcons[pep].append(seq)
            self.dict['PepSeq'] = pepcons
            ## ~ [3b] Classify peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['PepTypes'] = {'Common':[],'Cluster':[],'Unique':[]}
            for pep in pepcons:
                if len(pepcons[pep]) == 1: self.dict['PepTypes']['Unique'].append(pep); continue
                pepclus = []
                for seq in pepcons[pep]:
                    for cluster in clusters:
                        if seq in cluster and cluster not in pepclus: pepclus.append(cluster)
                if len(pepclus) == 1: self.dict['PepTypes']['Cluster'].append(pep)
                else: self.dict['PepTypes']['Common'].append(pep)
            ## ~ [3c] Summarise Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            RED.write('%d different Peptide sequences:\n' % len(pepcons))
            for ptype in ['Common','Cluster','Unique']: self.dict['PepTypes'][ptype].sort()
            RED.write('|--- %d Unique to one consensus: %s\n' % (len(self.dict['PepTypes']['Unique']),self.pepList(self.dict['PepTypes']['Unique'])))
            RED.write('|--- %d Resticted to one cluster: %s\n' % (len(self.dict['PepTypes']['Cluster']),self.pepList(self.dict['PepTypes']['Cluster'])))
            RED.write('|--- %d Common to multiple clusters: %s\n' % (len(self.dict['PepTypes']['Common']),self.pepList(self.dict['PepTypes']['Common'])))
            RED.write('\n________________________________________\n\n')

            ### ~ [4] Rename sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            calph = string.ascii_lowercase; cx = 0
            RED.write('%d Consensi Clusters:\n\n' % len(clusters))
            for type in ['Consensus','Description','ConsPep','UniqPep','ClusPep','CommPep']:
                if type not in myfiesta.list['Headers']: myfiesta.list['Headers'].append(type)
            newseq = []
            for cluster in clusters:
                cx += 1
                base = '%s%s' % (myfiesta.info['NewAcc'],rje.preZero(cx,len(clusters)))
                RED.write('>%s [%s]\n' % (base,rje.join(seqlist.accList(cluster),', ')))
                for seq in cluster:
                    newseq.append(seq)
                    oldacc = seq.info['AccNum']
                    try: seq.info['AccNum'] = '%s%s' % (base,calph[cluster.index(seq)])
                    except: seq.info['AccNum'] = '%sz%s' % (base,cluster.index(seq)-len(calph))
                    seq.info['ID'] = 'bud_%s' % myfiesta.info['SpCode']
                    seq.info['Name'] = '%s %s' % (seq.info['AccNum'],seq.info['Description'])
                    if myfiesta.opt['GnSpAcc']: seq.info['Name'] = '%s__%s' % (seq.info['ID'],seq.info['Name'])
                    RED.write('|--- %s\n' % seq.info['Name'])
                    for trans in myfiesta.dict['Mapping']:
                        if myfiesta.dict['Mapping'][trans]['BUD'] == oldacc:
                            myfiesta.dict['Mapping'][trans]['Consensus'] = seq.info['AccNum']
                            myfiesta.dict['Mapping'][trans]['Description'] = seq.info['Description']
                            myfiesta.dict['Mapping'][trans]['ConsPep'] = len(self.dict['Peptides'][oldacc])
                            myfiesta.dict['Mapping'][trans]['Exact'] = 0
                            myfiesta.dict['Mapping'][trans]['Coverage'] = 0
                            covseq = oldseq = caseseq = seq.info['Sequence'].lower()
                            for ptype in ['Common','Cluster','Unique']:
                                myfiesta.dict['Mapping'][trans][ptype] = 0
                                for pep in self.dict['PepTypes'][ptype]:
                                    if pep in self.dict['Peptides'][oldacc]:
                                        myfiesta.dict['Mapping'][trans][ptype] += 1
                                        pi = oldseq.find(pep.lower())
                                        if pi < 0: continue
                                        caseseq = caseseq[:pi] + pep.upper() + caseseq[pi+len(pep):]
                                        covseq = covseq[:pi] + '+' * len(pep) + covseq[pi+len(pep):]
                                        myfiesta.dict['Mapping'][trans]['Exact'] += 1
                            myfiesta.dict['Mapping'][trans]['Coverage'] = float(covseq.count('+')) / seq.aaLen()
                            myfiesta.dict['Mapping'][trans]['UniqPep'] = myfiesta.dict['Mapping'][trans].pop('Unique')
                            myfiesta.dict['Mapping'][trans]['ClusPep'] = myfiesta.dict['Mapping'][trans].pop('Cluster')
                            myfiesta.dict['Mapping'][trans]['CommPep'] = myfiesta.dict['Mapping'][trans].pop('Common')
                            seq.addSequence(caseseq,case=True)                            
                RED.write('\n')
            RED.write('________________________________________\n\n')
            seqlist.seq = newseq
            #x#myfiesta.saveMapping('Consensus')

            ### ~ [5] Peptide Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            peptdt = '%s.peptides.tdt' % self.info['Basefile']
            pephead = ['Peptide','Classification','Consensi']
            rje.delimitedFileOutput(self,peptdt,pephead,rje_backup=True)
            for ptype in ['Common','Cluster','Unique']:
                for pep in self.dict['PepTypes'][ptype]:
                    data = {'Peptide':pep,'Classification':ptype,'Consensi':seqlist.accList(self.dict['PepSeq'][pep])}
                    data['Consensi'].sort()
                    data['Consensi'] = rje.join(data['Consensi'],'|')
                    rje.delimitedFileOutput(self,peptdt,pephead,datadict=data)
            self.printLog('#PEP','Peptide details output to %s' % peptdt)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def multiHAQ(self,seqlist,batchonly=False): ### Perform multiple rounds of HAQESAC analysis on final consensi
        '''
        Perform multiple rounds of HAQESAC analysis on final consensi.
        >> seqlist:SeqList object containing sequences for HAQESAC analysis
        '''
        try:
            ### ~ [1] ~ Setup files for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqlist.seqNum(): return self.printLog('#HAQ','No sequences for HAQESAC analysis.')
            seqlist.list['Blast2Fas'] = [seqlist.info['Name'],self.info['SearchDB']]
            if self.info['HAQDB'].lower() not in ['','none']: seqlist.list['Blast2Fas'].append(self.info['HAQDB'])
            self.list['HitSeq'] = seqlist.seqs()
            batchfile = '%s_HAQESAC/haqesac.bat' % self.info['Basefile']
            rje.backup(self,batchfile)
            if rje.checkForFile(self.info['HAQDB']):
                seqlist.opt['Append'] = True
                rje_seq.Blast2Fas(seqlist)
                rje.mkDir(self,'%s_HAQESAC/' % self.info['Basefile'])
                for seq in seqlist.seqs():
                    try: os.rename('%s.blast.fas' % seq.info['AccNum'],'%s%s.fas' % (rje.makePath('%s_HAQESAC/' % self.info['Basefile']),seq.info['AccNum']))
                    except: self.errorLog('Problem moving %s.blast.fas' % seq.info['AccNum'],quitchoice=False)
                rje_blast.cleanupDB(seqlist.info['Name'],deletesource=False)
            else:
                if self.opt['HAQESAC']:
                    self.printLog('#HAQ','No HAQESAC analysis without SearchDB')
                    self.opt['HAQESAC'] = False
            ### ~ [2] ~ Run HAQESAC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            if not self.opt['HAQESAC']: batchonly = True    # return self.printLog('#HAQ','No HAQESAC analysis')
            if self.opt['MultiHAQ'] and not batchonly: hcyc = [0,1]
            else: hcyc = [0]
            haqcovered = []                 # List of Hits included in HAQESAC results
            ## ~ [2a] ~ Generate Batch file and INI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            haqcmd = ['ac=F','blastcut=50','maketree=fasttree','bootstraps=1000','alnprog=mafft','qcover=0','qryid=0','allowvar=F'] + self.cmd_list + ['newlog=F']
            if not self.opt['MultiHAQ']: haqcmd.insert(0,'i=-1') 
            for cmd in haqcmd[0:]:
                if cmd[:4].lower() == 'ini=': haqcmd.remove(cmd)
            open('%s_HAQESAC/haqesac.ini' % self.info['Basefile'],'w').write(rje.join(haqcmd,'\n'))
            for seq in self.list['HitSeq']:
                acc = seq.info['AccNum']
                haqcmd = ['root=mid','group=dup','qryvar=T','seqin=%s.fas' % acc, 'query=%s' % acc, 'basefile=%s' % acc]
                if self.opt['MultiHAQ']: haqcmd += ['multihaq=T']
                open(batchfile,'a').write('python %shaqesac.py %s\n' % (self.info['Path'],rje.join(haqcmd)))
            self.printLog('#HAQBAT','HAQESAC Batch file output to %s' % batchfile)
            if batchonly: return self.printLog('#HAQ','No HAQESAC analysis - run using %s' % batchfile)
            ## ~ [2b] ~ Run Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            for cyc in hcyc:
                self.printLog('#HAQ','Running HAQESAC - will have own log files etc.')
                for seq in self.list['HitSeq']:
                    if cyc == hcyc[-1] and seq in haqcovered and (self.opt['MultiHAQ'] or self.i() >= 0):
                        if self.i() < 0 or rje.yesNo('%s already covered by previous Hit HAQESAC. Skip?' % seq.shortName()):
                            self.printLog('#SKIP','%s already covered by previous Hit HAQESAC: Skipped' % seq.shortName()); continue
                    acc = seq.info['AccNum']
                    haqcmd = ['ac=F','blastcut=50','maketree=fasttree','bootstraps=1000','alnprog=maaft','qcover=0','qryid=0','allowvar=F','treeformats=bud,text,nsf']
                    haqcmd += self.cmd_list + ['root=mid','group=dup','qryvar=T','seqin=%s.fas' % acc, 'query=%s' % acc, 'basefile=%s' % acc, 'newlog=F']
                    if self.opt['MultiHAQ']: haqcmd += ['multihaq=T']
                    os.chdir('%s_HAQESAC' % self.info['Basefile'])
                    info = haqesac.makeInfo()
                    out = rje.Out(cmd_list=haqcmd)                    # Sets up Out object for controlling output to screen
                    out.printIntro(info)                                # Prints intro text using details from Info object
                    haqlog = rje.setLog(info,out,haqcmd)                 # Sets up Log object for controlling log file output
                    try:
                        haqesac.HAQESAC(log=haqlog, cmd_list=haqcmd).run(setobjects=True)
                        os.chdir(self.info['RunPath'])
                        self.printLog('#HAQ','HAQESAC run for %s' % seq.shortName())
                    except: 
                        os.chdir(self.info['RunPath'])
                        self.errorLog('Problem running HAQESAC for %s' % seq.shortName())
                    if cyc != hcyc[-1]: continue
                    if self.opt['CleanHAQ']:
                        keeplist = ['tree.txt','fas','png','pickle.gz','log','nsf']
                        for file in glob.glob('%s_HAQESAC/%s.*' % (self.info['Basefile'],acc)):
                            #x#if rje.replace(file,'%s_HAQESAC/%s.' % (self.info['Basefile'],acc),'') in keeplist: continue
                            if rje.split(file,'.')[-1] in keeplist: continue
                            if rje.join(rje.split(file,'.')[-2:],'.') in keeplist: continue
                            os.unlink(file)
                    seqlist.loadSeqs(rje.makePath('%s_HAQESAC/%s.fas' % (self.info['Basefile'],acc),True))
                    acclist = seqlist.accList()
                    for hseq in self.list['HitSeq']:
                        if hseq.info['AccNum'] in acclist: haqcovered.append(hseq)
            self.printLog('#HAQ','HAQESAC run - each dataset will have own log files etc.')
        except: os.chdir(self.info['RunPath']); self.errorLog('Problem with BUDAPEST.multiHAQ()')
#########################################################################################################################
    ### <3> ### Data setup Methods                                                                                      #
#########################################################################################################################
    def readMascot(self):  ### Reads the MASCOT file into self.dict['Hit'] dictionary
        '''
        Reads the MASCOT file into self.dict['Hit'] dictionary.
        1. Split MASCOT file into a header text file and delimited data file. (Check for and add emPAI headers if needed.)
        2. Read delimited data file into dictionary (lists for each header).
        '''
        try:### ~ [1] Split MASCOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log, self.cmd_list)
            mascot = rje_mascot.MASCOT(self.log,self.cmd_list+['basefile=None'])
            mascot.run()
            dfile = '%s.mascot.csv' % rje.baseFile(self.info['Name'])
            #dlines = open(dfile,'r').readlines()
            #!# Improve usinng rje_db.Database table #!#
            ### ~ [2] Read delimited data file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb = db.addTable(dfile,name='mascot',mainkeys=['prot_hit_num','pep_query'],datakeys='All')
            mdb.dataFormat({'prot_hit_num':'int'})
            #headers = rje.readDelimit(dlines.pop(0),',')    # Mascot delimited headers
            headers = mdb.fields()
            self.deBug(headers)
            self.dict['Hits'] = {}
            for protid in rje.sortKeys(mdb.index('prot_hit_num')):
                protacc = mdb.indexDataList('prot_hit_num',protid,'prot_acc')
                while '' in protacc: protacc.remove('')
                protacc = protacc[0]
                self.dict['Hits'][protacc] = {'Peptides':{}}
                for entry in mdb.indexEntries('prot_hit_num',protid):
                    pep = entry['pep_query']
                    self.dict['Hits'][protacc]['Peptides'][pep] = {}
                    for field in mdb.fields():
                        if field[:4] == 'prot' and entry[field]: self.dict['Hits'][protacc][field] = entry[field]
                        elif field != 'pep_query': self.dict['Hits'][protacc]['Peptides'][pep][field] = entry[field]
            self.printLog('#MASCOT','Mascot data read in for %s hits' % rje.integerString(len(self.dict['Hits'])))
            if not self.dict['Hits']:
                self.printLog('#MASCOT','Cannot proceed with No hits! Check MASCOT formatting and itraq=T/F setting.')
                return False
            ### ~ [3] Compress peptides to unique sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Support'] = {}       # Dictionary of {numbers of peptides:[proteins]}
            for prot in rje.sortKeys(self.dict['Hits']):
                pepdata = self.dict['Hits'][prot].pop('Peptides')
                self.dict['Hits'][prot]['Peptides'] = []
                for num in rje.sortKeys(pepdata):
                    #x#pepseq = '%s%s%s' % (pepdata[num]['pep_res_before'],pepdata[num]['pep_seq'],pepdata[num]['pep_res_after'])
                    #x#self.deBug(pepdata)
                    pepseq = pepdata[num]['pep_seq']
                    if pepseq not in self.dict['Hits'][prot]['Peptides']: self.dict['Hits'][prot]['Peptides'].append(pepseq)
                self.dict['Hits'][prot]['Peptides'].sort()
                px = len(self.dict['Hits'][prot]['Peptides'])
                if px not in self.dict['Support']: self.dict['Support'][px] = []
                self.dict['Support'][px].append(prot)
            self.printLog('#MASCOT','MASCOT hits compressed to unique peptide sequences')
            if self.i() < 0 or rje.yesNo('Continue?'): return True
            else: return False
        except: self.errorLog('BUDAPEST error reading MASCOT file'); return False
#########################################################################################################################
### End of SECTION II: Budapest Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
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
    try: Budapest(mainlog,cmd_list).run()
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
