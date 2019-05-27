#!/usr/local/bin/python

# HAQESAC - Homologue Alignment Quality, Establishment of Subfamilies and Ancestor Construction
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
Program:      HAQESAC
Description:  Homologue Alignment Quality, Establishment of Subfamilies and Ancestor Construction
Version:      1.13.0
Last Edit:    24/05/19
Citation:     Edwards et al. (2007), Nature Chem. Biol. 3(2):108-112. [PMID: 17220901]
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    HAQESAC is a tool designed for processing a dataset of potential homologues into a trustworthy gobal alignment and well
    bootstrap-supported phylogeny. By default, an intial dataset consisting of homologues (detected by BLAST, for example)
    is processed into a dataset consisting of the query protein and its orthologues in other species plus paralogous
    subfamilies in the same gene family.

    Individual sequences are therefore screened fulfil two criteria:
       1. They must be homologous to (and alignable with) the query sequence of interest.
       2. They must be a member of a subfamily within the gene family to which the query sequence belongs. 

    HAQ. The first stage of data cleanup is therefore to remove rogue sequences that either do not fit in the gene family at
    all or are too distantly related to the query protein for a decent alignment that can be used for useful further
    analysis. This is achieved firstly by a simple identity cut-off, determined by pairwise alignments of sequences, and
    then by a more complex procedure of removing sequences for whom the overall alignability is poor. During this procedure,
    sequences that have too many gaps are also removed as too many gapped residues can cause problems for downstream
    evolutionary analyses. Further screening is achieved based on phylogenetic information.

    ES. Once the dataset has been 'cleaned up' (and, indeed, during processing), HAQESAC can be used to assign sequences to
    subgroups or subfamilies, if such information is needed for downstream analyses.

    AC. The final step that HAQESAC is able to perform is ancestral sequence prediction using the GASP (Gapped Ancestral
    Sequence Prediction) algorithm (Edwards & Shields 2005)

    By default, HAQESAC will perform all these operations. However, it is possible to turn one or more off and only, for
    example, reject individually badly aligned sequences, if desired. Details can be found in the Manual:
    https://github.com/slimsuite/SLiMSuite/blob/master/docs/manuals/HAQESAC%20Manual.pdf.

    HAQESAC outputs can be controlled by d=X:
    0 - *.fas, *.grp, *.nwk, *.tree.txt, *.png, *.anc.*
    1 - *.bak (seqin), *.clean.fas, *.haq.fas, *.degap.fas
    2 - *.saq.*.fas, *.paq.*.fas, *.scanseq.fas
    3 - *.saqx.*.fas

Commandline Options:
    # General Dataset Input/Output Options # 
    seqin=FILE  : Loads sequences from FILE (fasta, phylip, aln, uniprot format, or list of fastacmd names for fasdb option)
    query=X     : Selects query sequence by name (or part of name, e.g. Accession Number)
    acclist=X   : Extract only AccNums in list. X can be FILE or list of AccNums X,Y,.. [None]
    fasdb=FILE  : Fasta format database to extract sequences from [None]
    basefile=X  : Basic 'root' for all files X.* [By default will use 'root' of seqin=FILE if given or haq_AccNum if qblast]
    v=X         : Sets verbosity (-1 for silent) [0]
    i=X         : Sets interactivity (-1 for full auto) [0]
    d=X         : Data output level (0-3, see docstring) [1]
    resdir=PATH : Output directory for d>0 outputs [./HAQESAC/]
    log=FILE    : Redirect log to FILE [Default = calling_program.log or basefile.log]
    newlog=T/F  : Create new log file. [Default = False: append log file]
    multihaq=T/F: If pickle present, will load and continue, else will part run then pickle and stop [False]
    
    # Pre-HAQESAC Data selection #
    qseqfile=FILE   : Sequence file from which to extract query (query=X) and peform BLAST [None]
    qblast=FILE     : BLAST database against which to BLAST query (see rje_blast.py options) [None] ['blaste=1e-7','blastv=1000,blastb=0']
    qfastacmd=X     : Extract query X from qblast database using fastacmd (may also need query=X) [None]

    # Pre-HAQESAC Sequence Filtering #
    gnspacc=T/F     : Convert sequences into gene_SPECIES__AccNum format wherever possible. [True] 
    backup=T/F      : Whether to backup initial fasta file. Will overwrite existing *.fas.bak. [True]
    accnr=T/F       : Check for redundant Accession Numbers/Names on loading sequences. [False]
    #!# filterseq=FILE  : Filters out sequences in given file (Log, Fasta or list of names)
    #!# filterspec=FILE : Filters out sequences according to species (codes) listed in given file
    unkspec=T/F     : Whether sequences of unknown species are allowed [True]
    9spec=T/F       : Whether to treat 9XXXX species codes as actual species (generally higher taxa) [False]
    dblist=X,Y,..,Z : List of databases in order of preference (good to bad)
    dbonly=T/F      : Whether to only allow sequences from listed databases [False]
    #!# keepdesc=FILE   : Only keeps sequences with 1+ of text listed in given file in sequence description
    minlen=X        : Minimum length of sequences [0]
    maxlen=X        : Maximum length of sequences (<=0 = No maximum) [0]
    goodX=LIST      : Filters where only sequences meeting the requirement of LIST are kept.
                      - LIST may be a list X,Y,..,Z or a FILE which contains a list [None]
                        - goodacc  = list of accession numbers
                        - goodseq  = list of sequence names
                        - goodspec = list of species codes
                        - gooddb   = list of source databases
                        - gooddesc = list of terms that, at least one of which must be in description line
    badX=LIST       : As goodX but excludes rather than retains filtered sequences

    # Pre-HAQ Data Cleanup #
    cleanup=T/F : Initial data cleanup [True]
    seqnr=T/F   : Make sequence Non-Redundant [True]
    specnr=T/F  : Non-Redundancy within same species only [True]
    nrid=X      : %Identity cut-off for Non-Redundancy [99.0]
    blastcut=X  : Maximum number of sequences to have in dataset (BLAST query against NR dataset.) [0]
    blaste=X    : E-Value cut-off for BLAST searches (BLAST -e X) [1e-4]
    blastv=X    : Number of one-line hits per query (BLAST -v X) [500]
    blastf=T/F  : Complexity Filter (BLAST -F X) [True]
    qcover=X    : Min. % (ordered) BLAST coverage of Query vs Hit or Hit vs Query [60.0]
    gablam=T/F  : Whether to use GABLAMO Global Alignment from BLAST Local Alignment Matrix (Ordered) rather than ALIGN [True]
    gabsim=T/F  : Whether to use %Similarity for GABLAMO comparisons, rather than %Identity [True]
    pairid=X    : Fasta Alignment ID cut-off for any pair of sequences [0.0]
    qryid=X     : Fasta Alignment ID with Query cut-off [40.0]
    maxgap=X    : Maximum proportion of sequence that may be gaps compared to nearest neighbour (<=0 = No maximum) [0.5]

    # General HAQ Options #
    qregion=X,Y     : Concentrate on the region of the query from (and including) residue X to residue Y [1,-1]
    haq=T/F         : Homologue Alignment Quality [True]
    noquery=T/F     : No Query for SAQ, Random Query for PAQ (else query=X or first sequence) [False]
    keep=T/F        : Keep all sequences (saqkl=0, paqkl=0) [False]
    usealn=T/F      : Use current alignment (do not realign, degap=F) [False]
    cwcut=X         : Total number of residues above which to use ClustalW for alignment in place of alnprog=X [0]
    haqmatrix=FILE  : File of AA vs AA scores used in SAQ and PAQ. [None]

    # Single Sequence Alignment Quality (SAQ) Options #
    saq=T/F     : Single Sequence AQ [True]
    saqc=X      : Min score for a residue in SAQ (Default matrix: no. seqs to share residue). [2]
    saqb=X      : SAQ Block length. [10]
    saqm=X      : No. residues to match in SAQ Block. [7]
    saqks=X     : Relative Weighting of keeping Sequences in SAQ. [3]
    saqkl=X     : Relative Weighting of keeping Length in SAQ. [1]
    mansaq=T/F  : Manual over-ride of sequence rejection decisions in SAQ [False]

    # Pairwise Alignment Quality (PAQ) Options #
    prepaq=T/F  : PrePAQ tree grouping [True]
    paq=T/F     : Pairwise AQ [True]
    paqb=X      : PAQ Block length. [7]
    paqm=X      : Min score in PAQ Block (Default matrix: No. residues to match). [3]
    paqks=X     : Relative Weighting of keeping Sequences in PAQ. [3]
    paqkl=X     : Relative Weighting of keeping Length in PAQ. [1]
    manpaq=T/F  : Manual over-ride of sequence rejection decisions in PAQ [False]

    # Establishment of Subfamilies (and PrePAQ) Options #
    es=T/F      : Establishment of Subfamilies [True]
    root=X      : Rooting of tree (rje_tree.py), where X is:
        - mid = midpoint root tree. [Default]
        - ran = random branch.
        - ranwt = random branch, weighted by branch lengths.
        - man = always ask for rooting options (unless i<0).
        - FILE = with seqs in FILE as outgroup. (Any option other than above)
    bootcut=X   : cut-off percentage of tree bootstraps for grouping.
    mfs=X       : minimum family size [3]
    fam=X       : minimum number of families (If 0, no subfam grouping) [0]
    orphan=T/F  : Whether orphans sequences (not in subfam) allowed. [True]
    allowvar=T/F: Allow variants of same species within a group. [False]
    qryvar=T/F  : Keep variants of query species within a group (over-rides allowvar=F). [False]
    groupspec=X : Species for duplication grouping [None]
    qspec=T/F   : Whether to highlight query species in PNG tree files [True]
    specdup=X   : Minimum number of different species in clade to be identified as a duplication [2]
    group=X     : Grouping of tree
        - man = manual grouping (unless i<0).
        - dup = duplication (all species unless groupspec specified).
        - qry = duplication with species of Query sequence (or Sequence 1) of treeseq
        - one = all sequences in one group
        - None = no group (case sensitive)
        - FILE = load groups from file
    phyoptions=FILE : File containing extra Phylip tree-making options ('batch running') to use [None]
    protdist=FILE   : File containing extra Phylip PROTDIST ('batch running') to use [None]
    maketree=X      : Program for making tree [None]
        - None = Do not make tree from sequences 
        - clustalw = ClustalW NJ method
        - neighbor = PHYLIP NJ method (NB. Bootstraps not yet supported)
        - upgma    = PHYLIP UPGMA (neighbor) method (NB. Bootstraps not yet supported)
        - fitch    = PHYLIP Fitch method (NB. Bootstraps not yet supported)
        - kitsch   = PHYLIP Kitsch (clock) method (NB. Bootstraps not yet supported)
        - protpars = PHYLIP MP method (NB. Bootstraps not yet supported)
        - proml    = PHYLIP ML method (NB. Bootstraps not yet supported)
        - PATH     = Alternatively, a path to a different tree program/script can be given. This should accept ClustalW parameters.

    # Ancestor Construction (GASP) Options
    ac=T/F      : Ancestor Construction (GASP) [True]
    pamfile=FILE: Sets PAM1 input file [jones.pam]
    pammax=X    : Initial maximum PAM matrix to generate [100]
    pamcut=X    : Absolute maximum PAM matrix [1000]
    fixpam=X    : PAM distance fixed to X [0].
    rarecut=X   : Rare aa cut-off [0.05].
    fixup=T/F   : Fix AAs on way up (keep probabilities) [True].
    fixdown=T/F : Fix AAs on initial pass down tree [False].
    ordered=T/F : Order ancestral sequence output by node number [False].
    pamtree=T/F : Calculate and output ancestral tree with PAM distances [True].
    desconly=T/F: Limits ancestral AAs to those found in descendants [True].
    xpass=X     : How many extra passes to make down & up tree after initial GASP [1].

    # System Info Options #
    See rje_seq and general commandline options for external program system path settings.

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_ancseq, rje_haq, rje_seq, rje_tree
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added clustalw for early stages
    # 1.0 - Fully working version. (Improvements still to be made!) Added BLAST to front-end.
    # 1.1 - Added backup.
    # 1.2 - Added basename and interactive menu to start if nothing given
    # 1.3 - Added specification of query subportion. Changed scoring to use sub matrix. Added GASP.
    # 1.4 - Added option to use BLAST coverage rather than ALIGN
    # 1.5 - Added PHYLIP for tree drawing (rje_tree.py 1.6)
    # 1.6 - Added multiHAQ for more rapid treatment of multiple consecutive datasets (and "pickup")
    # 1.7 - Added qryvar=T/F  : Allow variants of query species within a group (over-rides allowvar=F). [False]
    # |---- Fixed minor typo.
    # 1.8 - Added qspec tree option
    # 1.9 - Added rje_blast_V2 implementation and BLAST+. Use oldblast=T for old BLAST.
    # 1.10- Added exceptions for BLAST failure.
    # 1.10.1 - Tweaked QryVar interactivity.
    # 1.10.2 - Corrected typos and disabled buggy post-HAQESAC data reduction.
    # 1.10.3 - Added catching of bad query when i=-1.
    # 1.11.0 - Added resdir=PATH [./HAQESAC/] for d>0 outputs.
    # 1.12.0 - 9spec=T/F   : Whether to treat 9XXXX species codes as actual species (generally higher taxa) [False]
    # 1.13.0 - Modified qregion=X,Y to be 1-L numbering.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Fix bug when returning to group menu from group review.
    # [ ] : General tidy up.
    # [ ] : Think about using entropy calculations somehow in place of crude SAQ/PAQ calculations.
    # [ ] : Make sure that HAQESAC is checking sequence numbers and quitting if appropriate. (Too few for alignment/tree)
    # [ ] : Fix the bug that gives dodgy output after saving "full" copies of files. (Remove this function.)
    # [ ] : Add memory of selection choices for group sequence selection.
    # [ ] : Fix SeqNR bug.
    # [ ] : Fix i=1 bug.
    # [ ] : Add REST output.
    # [ ] : Convert to new object type.
    # [Y] : Add option to screen out 9XXXX species codes for duplication identification.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('HAQESAC', '1.13.0', 'May 2019', '2005')
    description = 'Homologue Alignment Quality, Establishment of Subfamilies and Ancestor Construction'
    author = 'Dr Richard J. Edwards.'
    comments = ['Cite: Edwards RJ et al. (2007) Nature Chem. Biol. 3(2):108-112.']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
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
            if rje.yesNo('Show rje_seq commandline options?',default='N'): out.verbose(-1,4,text=rje_seq.__doc__)
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
### SECTION II: HAQESAC Class                                                                                           #
#########################################################################################################################
class HAQESAC(rje.RJE_Object):     
    '''
    HAQESAC Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of sequence file (also temp seqin)
    - Basefile = root for all files. root.*
    - QrySeqFile = Sequence file from which to extract query (query=X) and peform BLAST [None]
    - QryBlast = BLAST database against which to BLAST query (see rje_blast.py options) [None]
    - QryFastaCmd = Extract query X from qblast database using fastacmd (may also need query=X) [None]
    - QryRegion = X,Y : Concentrate on the region of the query from (and including) residue X to residue Y
    - ResDir=PATH : Output directory for d>0 outputs [./HAQESAC/]

    Opt:boolean
    - Backup = Whether to store backup of original input file (*.fas.bak)
    - MultiHAQ = If pickle present, will load and continue, else will part run then pickle and stop [False]
    - KeepAll = Keep all sequences (saqkl=0, paqkl=0) [False]    
    - UseAln = Use current alignment (do not realign, degap=F) [False]
    - GABLAM = Whether to use GABLAMO Global Alignment from BLAST Local Alignment Matrix (Ordered) rather than ALIGN [True]
    - GABSim = Whether to use %Similarity for GABLAMO comparisons, rather than %Identity [True]
    - CleanUp = Whether to perform pre-HAQ cleanup [True]
    - HAQ = Homologue Alignment Quality [True]
    - SAQ = Single Sequence AQ [True]
    - PrePAQ = PrePAQ tree grouping [True]
    - PAQ = Pairwise AQ [True]
    - QSpec = Whether to highlight query species in PNG tree files [True]
    - ES = Establishment of Subfamilies [True]
    - AC = Ancestor Construction (GASP) [True]    

    Stat:numeric
    - BlastCut = Maximum number of sequences to have in dataset (BLAST query against NR dataset.)
    - QCover = Min. % (ordered) BLAST coverage of Query vs Hit or Hit vs Query [60.0]
    - PairID = Fasta Alignment ID cut-off for any pair of sequences [0.0]
    - QryID = Fasta Alignment ID with Query cut-off [40.0]
    - CWCut = Total number of residues above which to use ClustalW for alignment, not MUSCLE [0]
    - DataLevel = Data Output level from 0-3

    Obj:RJE_Objects
    - SeqList:rje_seq.SeqList Object
    - Tree:rje_tree.Tree Object
    - HAQ:rje_haq.HAQ Object
    - GASP:rje_ancseq.Gasp Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','Basefile','QrySeqFile','QryBlast','QryFastaCmd','QryRegion','ResDir']
        self.statlist = ['BlastCut','PairID','QryID','CWCut','DataLevel','QCover']
        self.optlist = ['Backup','KeepAll','UseAln','GABLAM','GABSim','CleanUp','HAQ','SAQ','PrePAQ','PAQ','ES','AC','MultiHAQ','QSpec']
        self.objlist = ['SeqList','Tree','HAQ','GASP']
        self.listlist = []
        self.dictlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None)   #,setlist=True,setdict=True)
        self.setInfo({'QryRegion':'1,-1','ResDir':rje.makePath('./HAQESAC/')})
        self.setStat({'QCover':60.0,'QryID':40.0,'CWCut':0,'DataLevel':1})
        self.setOpt({'KeepAll':False,'UseAln':False,'MultiHAQ':False})
        ### Other Attributes ###
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
                self._cmdRead(type='info',att='Name',arg='seqin',cmd=cmd)
                self._cmdRead(type='info',att='QrySeqFile',arg='qseqfile',cmd=cmd)
                self._cmdRead(type='info',att='QryBlast',arg='qblast',cmd=cmd)
                self._cmdRead(type='info',att='QryFastaCmd',arg='qfastacmd',cmd=cmd)
                self._cmdRead(type='info',att='Basefile',cmd=cmd)
                self._cmdRead(type='int',att='BlastCut',arg='blastcut',cmd=cmd)
                self._cmdReadList(cmd,'perc',['QryID','PairID','QCover'])
                self._cmdReadList(cmd,'path',['ResDir'])
                self._cmdRead(type='stat',att='CWCut',arg='cwcut',cmd=cmd)
                self._cmdRead(type='int',att='DataLevel',arg='d',cmd=cmd)
                self._cmdReadList(cmd,'opt',['Backup','MultiHAQ','QSpec'])
                self._cmdRead(type='opt',att='KeepAll',arg='keep',cmd=cmd)
                self._cmdRead(type='opt',att='UseAln',arg='usealn',cmd=cmd)
                self._cmdRead(cmd,type='opt',att='GABLAM')
                self._cmdRead(cmd,type='opt',att='GABLAM',arg='gablamo')
                self._cmdRead(cmd,type='opt',att='GABSim')
                self._cmdRead(type='opt',att='CleanUp',arg='cleanup',cmd=cmd)
                self._cmdRead(type='info',att='QryRegion',arg='qregion',cmd=cmd)
                self._cmdRead(type='opt',att='HAQ',arg='haq',cmd=cmd)
                self._cmdRead(type='opt',att='SAQ',arg='saq',cmd=cmd)
                self._cmdRead(type='opt',att='PrePAQ',arg='prepaq',cmd=cmd)
                self._cmdRead(type='opt',att='PAQ',arg='paq',cmd=cmd)
                self._cmdRead(type='opt',att='ES',arg='es',cmd=cmd)
                self._cmdRead(type='opt',att='AC',arg='ac',cmd=cmd)                
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def setObjects(self):   ### Sets HAQESAC objects (SeqList, Tree, HAQ)
        '''
        Sets HAQESAC objects (SeqList, Tree, HAQ).
        '''
        try:
            ### <0> ### Setup: pre HAQESAC BLAST 
            _stage = '<0> Setup'
            self.preBLAST()
            
            ### <1> ### SeqList
            _stage = '<1> SeqList'
            ## <1a> ## Setup
            seqcmd = ['maxgap=0.5','gapfilter=F','accnr=T','nrid=99.0','specnr=T','degap=F','gnspacc=T'] + self.cmd_list + ['seqin=None'] 
            if self.opt['UseAln']:  #!# MaxGap problem - temp change?
                seqcmd.append('tidygap=F')
            self.obj['SeqList'] = rje_seq.SeqList(log=self.log, cmd_list=seqcmd)
            #print seqcmd, self.obj['SeqList'].opt
            try:
                #print self.info['Name']
                #self.obj['SeqList'].info['Name'] = self.info['Name']
                self.obj['SeqList'].loadSeqs(seqfile=self.info['Name'])
            except:
                self.obj['SeqList'].seq = []
            ## <1b> ## Interactive sequence name provision
            while self.obj['SeqList'].seqNum() == 0:
                if self.obj['SeqList'].info['Name'] == 'None':
                    self.verbose(0,4,'\nNo sequence file given.',1)
                else:
                    self.log.errorLog('Sequence file (%s) contains no sequences.' % self.obj['SeqList'].info['Name'],printerror=False)
                if self.stat['Interactive'] < 0:
                    sys.exit()
                else:
                    self.info['Name'] = rje.choice(text='\nInput filename? (Blank to exit.)')
                    if self.info['Name'] == '':
                        sys.exit()
                    try:
                        self.cmd_list = self.cmd_list + rje.inputCmds(self,self.cmd_list)
                        self.obj['SeqList'].loadSeqs(seqfile=self.info['Name'])
                    except:
                        self.obj['SeqList'].seq = []
                        
            ## <1c> ## Setup rest of Seqlist
            self.obj['SeqList']._checkAln(aln=self.opt['UseAln'])
            self.obj['SeqList'].addMatrix('PWAln ID')
            self.obj['SeqList'].addMatrix('MSA ID')
            self.obj['SeqList'].addMatrix('MSA Gaps')
            if self.info['Basefile'] == 'None':
                self.obj['SeqList'].makeBaseFile()
                self.info['Basefile'] = self.obj['SeqList'].info['Basefile']
                #if re.search('^(\S+)\.(\S+)', self.obj['SeqList'].info['Name']):
                #    self.info['Basefile'] = rje.matchExp('^(\S+)\.(\S+)', self.obj['SeqList'].info['Name'])[0]
                #else:
                #    self.info['Basefile'] = rje.matchExp('^(\S+)', self.obj['SeqList'].info['Name'])[0]
            else:
                self.obj['SeqList'].info['Name'] = '%s.fas' % self.info['Basefile']
                self.obj['SeqList'].info['Basefile'] = self.info['Basefile']
            self.info['Name'] = self.obj['SeqList'].info['Name']
            ## <1d> ## Backup
            _stage = '<1d> Backup'
            if self.opt['Backup'] and rje.checkForFile(self.obj['SeqList'].info['Name']) and self.stat['DataLevel'] > 0:
                try:
                    bakfile = '%s%s.bak' % (self.getStr('ResDir'),self.obj['SeqList'].info['Name'])
                    BAK = open(bakfile, 'w')
                    FAS = open(self.obj['SeqList'].info['Name'], 'r')
                    BAK.writelines(FAS.readlines())
                    BAK.close()
                    FAS.close()
                except:
                    self.log.errorLog('Problem making backup file, %s.bak' % self.obj['SeqList'].info['Name'], printerror=False, quitchoice=True)
            
            ### <2> ### HAQ
            _stage = '<3> HAQ'
            if self.opt['KeepAll']:
                seqcmd += ['saqkl=0','paqkl=0']
            haq = self.obj['HAQ'] = rje_haq.HAQ(log=self.log, cmd_list=seqcmd)

            ### <3> ### Tree    
            _stage = '<2> Tree'
            treecmd = ['maketree=clustalw','root=mid','bootstraps=1000','specdup=2']
            if self.opt['QSpec'] and not haq.opt['NoQuery']: treecmd.append('treeformats=nwk,text,qspec')
            else: treecmd.append('treeformats=nwk,text,png')
            self.obj['Tree'] = rje_tree.Tree(log=self.log, cmd_list=treecmd+seqcmd+['disin=None','autoload=F'])            
            self.obj['Tree'].obj['SeqList'] = self.obj['SeqList']                               

            ### <4> ### GASP
            _stage = '<4> GASP'
            if self.opt['AC']:
                self.obj['GASP'] = rje_ancseq.Gasp(tree=self.obj['Tree'],ancfile='%s%s' % (self.getStr('ResDir'),self.info['Basefile']),cmd_list=seqcmd+['autofilter=F'],log=self.log)
        except SystemExit:
            sys.exit()
        except:
            self.log.errorLog('Major Problem with HAQESAC setObject (%s).' % _stage, quitchoice=True)            
#########################################################################################################################
    def preBLAST(self):     ### Makes a sequence file from a BLAST search                                           #V1.9
        '''
        Makes a sequence file from a BLAST search, using self.cmd_list.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('QryBlast').lower() in ['','none']: return False
            if self.getStr('QrySeqFile').lower() in ['','none'] and self.getStr('QryFastaCmd').lower() in ['','none']: return False
            blast_cmd = ['blaste=1e-7','blastv=1000','blastb=0'] + self.cmd_list + ['blastd=%s' % self.info['QryBlast']]
            blast = rje_blast.blastObj(log=self.log, cmd_list=blast_cmd,type='Dev')
            blast.formatDB(fasfile=self.info['QryBlast'],protein=True,force=False)
            
            ### ~ [1] Query Sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## This either pulls out the given query sequence form the QrySeqFile if given (using rje_seq)
            ## -OR- uses fastacmd to pull out qfastacmd from the qblast database
            if self.getStr('QrySeqFile').lower() in ['','none']:
                qseq = rje_seq.SeqList(log=self.log, cmd_list=['seqin=%s' % self.info['QrySeqFile'],'autofilter=F'])
                qseq.cmd_list = self.cmd_list
                if qseq.querySeq(): qseq.seq = [qseq.obj['QuerySeq']]   # Query found
                else: qseq.seq = qseq.seq[0:1]                          # No Query found - use first sequence
            else:   # FastaCmd
                qseq = rje_seq.SeqList(log=self.log, cmd_list=self.cmd_list+['seqin=None'])
                if blast.getBool('OldBLAST'): qseq.seq = [qseq.seqFromFastaCmd(self.info['QryFastaCmd'],self.info['QryBlast'])]
                else: qseq.seq = [qseq.seqFromBlastDBCmd(self.info['QryFastaCmd'],self.info['QryBlast'])]
            self.printLog('#QRY','Query sequence %s read for BLAST search against %s.' % (qseq.seq[0].shortName(),self.info['QryBlast']))

            ### ~ [2] BLAST Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qseq.saveFasta(seqfile='%s.qry' % qseq.seq[0].info['AccNum'])   # Save sequence
            blast.setStr({'InFile':'%s.qry' % qseq.seq[0].info['AccNum'],'Name':'%s.blast' % qseq.seq[0].info['AccNum']})
            blast.blast(); os.unlink(blast.getStr('InFile'))
            if not blast.readBLAST(unlink=True): raise ValueError

            ### ~ [3] Make new sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Name'] = 'None'
            for cmd in self.cmd_list: self._cmdRead(type='info',att='Name',arg='seqin',cmd=cmd) # Look for instructions
            if self.info['Basefile'] != 'None':         # If basefile given, will use this
                self.info['Name'] = '%s.fas' % self.info['Basefile']
                self.cmd_list.append('seqin=%s' % self.info['Name'])
            if self.info['Name'].lower() in ['','none']:    # Make new filename = ID or accnum
                if qseq.seq[0].info['DBase'] == 'sprot':  # If SwissProt will use the gene ID for the file name
                    self.info['Name'] = '%s.fas' % qseq.seq[0].info['ID'].lower()
                else:                    
                    self.info['Name'] = 'haq_%s.fas' % qseq.seq[0].info['AccNum'].lower()
                self.cmd_list.append('seqin=%s' % self.info['Name'])
            self.cmd_list = ['query=%s' % qseq.seq[0].info['AccNum']] + self.cmd_list
            blastseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None'])
            blastseq.info['Name'] = self.info['Name']
            if blast.getBool('OldBLAST'): 
                for hit in blast.search[0].hit:
                    blastseq.seqFromFastaCmd(hit.info['Name'],self.info['QryBlast'])
            else: blastseq.seqFromBlastDBCmd(blast.queryHits(self,qseq.seq[0].shortName()),self.info['QryBlast'])
            blastseq.saveFasta()    # Does this have the Query Sequence in it?
            #!# Should this check that the query sequence is in blastseq and warn if not?
        except SystemExit: sys.exit()
        except: self.errorLog('Major Problem with HAQESAC.preBLAST()', quitchoice=True)
#########################################################################################################################
    ### <2> ### Main HAQESAC Run
#########################################################################################################################
    def _aln(self,seqlist=None):    ### Aligns sequences
        '''
        Aligns sequences.
        >> seqlist:SeqList object
        << aligned SeqList object
        '''
        try:
            if seqlist.info['AlnProg'].lower() == 'muscle' and self.stat['CWCut'] > 0 and seqlist.aaCount() > self.stat['CWCut']:
                self.printLog('#CWCUT','%s aa => Using ClustalW' % rje.integerString(seqlist.aaCount()))
                return seqlist.align(alnprog='clustalw',outfile=seqlist.info['Name'],clustalbackup=False)
            else: return seqlist.align()
        except: self.errorLog('Major Problem with _aln().', quitchoice=True); raise
#########################################################################################################################
    def run(self,setobjects=False):  ### Main Run Backbone
        '''
        Main HAQESAC Run Backbone.
        <1> Cleanup
        <2> Homologue Alignment Quality
         <a> SAQ
         <b> PrePAQ
         <c> PAQ
        <3> Establishment of Subfamilies
        <4> Ancestor Construction
        <5> Final Adjustments
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            #if self.stat['Interactive'] >= 1:
            #    self.runMenu() #!' Include Manual sequence review!
            rje.mkDir(self,self.getStr('ResDir'))
            haqpickle = self.unpickleMe()
            if haqpickle and not self.opt['MultiHAQ'] and self.i() >= 0 and not rje.yesNo('HAQESAC pickle found and extracted. Would you like to use this and continue the previous analysis?'): haqpickle = None
            self.opt['MultiKill'] = self.opt['MultiHAQ'] and not haqpickle  # Conditions to pickle then die!
            if setobjects and not haqpickle: self.setObjects()
            elif haqpickle:
                try: haqpickle.obj['Tree']._cutCmdList(self.cmd_list)
                except: self.errorLog('Could not update HAQESAC Tree Object commandline options')
                for opt in ['PrePAQ','PAQ','ES','AC']: haqpickle.opt[opt] = self.opt[opt]
                self = haqpickle; self.opt['MultiKill'] = False
            seqlist = self.obj['SeqList']
            haq = self.obj['HAQ']
            tree = self.obj['Tree']
            #!#
            if 'SpecDup' not in tree.stat:
                tree.stat['SpecDup'] = 2
                for node in tree.node: node.opt['SpecDup'] = False
            #!#
            if 'TreeFormats' not in tree.list: tree.list['TreeFormats'] = ['te','nsf','text','png']
            if not haqpickle and not haq.opt['NoQuery']:
                if seqlist.obj['QuerySeq'] == None:
                    if self.i() < 0:
                        raise IOError('Query "%s" not found and noquery=F. Quitting HAQESAC.' % seqlist.getStr('Query'))
                    if rje.yesNo('No Query given but noquery=F. Use %s as query sequence?' % seqlist.seq[0].shortName()):
                        seqlist.obj['QuerySeq'] = seqlist.seq[0]
                        self.log.printLog('#QRY','Query Sequence = %s.' % seqlist.obj['QuerySeq'].shortName())
                    else:
                        newquery = rje.choice('Enter ID or AccNum of query. (Blank to exit.)')
                        seqlist.cmd_list.append('query=%s' % newquery)
                        if seqlist.querySeq() == False:
                            self.log.errorLog('Query %s not found. Quitting HAQESAC.' % newquery,printerror=False)
                            sys.exit()

            ### <1> ### Pre-HAQESAC Cleanup
            _stage = '<1> Cleanup'
            if self.opt['CleanUp'] and not haqpickle:
                seqlist = self.cleanUp()
                if self.stat['DataLevel'] > 0:
                    seqlist.saveFasta(seqfile='%s%s.clean.fas' % (self.getStr('ResDir'),self.info['Basefile']))
                
            self.stat['QryID'] = 0  # No Need to repeat this!

            ### <2> ### Homologue Alignment Quality
            _stage = '<2> HAQ'

            ## <a> ## SAQ: Single Sequence Alignment Quality => to remove true rogues. Can work with or without master sequence.
            if self.opt['HAQ'] and self.opt['SAQ'] and not haqpickle:
                _stage = '<2a> SAQ'
                #seqlist._checkAln(aln=True,realign=True)
                seqnum = [0,seqlist.seqNum()]
                haq.stat['SAQCyc'] = 0
                for seq in seqlist.seq: seq.info['SAQ'] = seq.info['Sequence'][0:]      
                while seqnum[-2] != seqnum[-1] and seqlist.seqNum() > 1:
                    if (haq.stat['SAQCyc'] > 0 and self.opt['UseAln'] == False) or seqlist.opt['Aligned'] == False:
                        self._IDGapFilter(seqlist)  # Handles realignment, %ID and GapFilter
                    if seqnum[-1] != seqlist.seqNum():
                        seqnum.append(seqlist.seqNum())
                        continue
                    haq.stat['SAQCyc'] += 1
                    self.log.printLog('#INF','Single Sequence Alignment Quality (SAQ) %d for %s (%d Seqs).' % (haq.stat['SAQCyc'], seqlist.info['Name'],seqlist.seqNum()))
                    haq.singleSeqAQ(seqlist,focus=self._queryFocus())
                    seqnum.append(seqlist.seqNum())
                    self.log.printLog('#HAQ','SAQ%d: %d -> %d sequences.' % (haq.stat['SAQCyc'],seqnum[-2],seqnum[-1]))
                    if self.stat['DataLevel'] > 2:
                        if self.opt['UseAln'] == False:
                            seqlist.tidyGaps(key='SAQX')
                        haq.saveHAQ(seqlist,filename='%s%s.saqx.%d.fas' % (self.getStr('ResDir'),self.info['Basefile'],haq.stat['SAQCyc']),key='SAQX')
                    if self.stat['DataLevel'] > 1:
                        if self.opt['UseAln'] == False:
                            seqlist.tidyGaps(key='SAQ')
                        haq.saveHAQ(seqlist,filename='%s%s.saq.%d.fas' % (self.getStr('ResDir'),self.info['Basefile'],haq.stat['SAQCyc']),key='SAQ')
                if self.opt['UseAln'] == False:
                    seqlist.tidyGaps()
                if self.stat['DataLevel'] > 0:
                    if self.opt['UseAln'] == False:
                        seqlist.tidyGaps(key='SAQ')
                    haq.saveHAQ(seqlist,filename='%s%s.haq.fas' % (self.getStr('ResDir'),self.info['Basefile']),key='SAQ')
                self.verbose(0,1,'\nFinished Single Sequence Alignment Quality (%d Cycles):\n => %d of %d sequences retained.' % (haq.stat['SAQCyc'],seqnum[-1],seqnum[1]),1)
                seqlist.saveFasta()

            ### <b> ### PAQ: Pairwise Alignment Quality
            _stage = '<2b> PAQ'
            if self.opt['MultiKill']: self.pickleMe(); sys.exit(0)                
            if self.opt['HAQ'] and (self.opt['PAQ'] or self.opt['PrePAQ']):
                _stage = 'PAQ'
                seqnum = [0,seqlist.seqNum()]
                if self.opt['PAQ']:
                    for seq in seqlist.seq: seq.info['PAQ'] = seq.info['Sequence'][0:]      
                while seqnum[-2] != seqnum[-1] and seqlist.seqNum() > 1:
                    ## <i> ## PrePAQ
                    if self.opt['PrePAQ']: self.prePAQ()   #!# Make part of rje_haq? Treat special the first time (famtrim)
                    self.pickleMe()
                    seqnum.append(seqlist.seqNum())
                    ## <ii> ## SAQ - *if* prePAQ removes sequences
                    if seqnum[-2] != seqnum[-1] and self.opt['KeepAll']:
                        self.printLog('#PREQ','PrePAQ Tree Edit: %d -> %d sequences.' % (seqnum[-2],seqnum[-1]))
                        self._IDGapFilter(seqlist)  # Handles realignment, %ID and GapFilter
                        if seqnum[-1] != seqlist.seqNum(): seqnum.append(seqlist.seqNum())
                        continue
                    if seqnum[-2] != seqnum[-1] and self.opt['SAQ']:
                        self._IDGapFilter(seqlist)  # Handles realignment, %ID and GapFilter
                        if seqnum[-1] != seqlist.seqNum(): seqnum.append(seqlist.seqNum()); continue
                        haq.stat['SAQCyc'] += 1
                        self.log.printLog('#INF','Single Sequence Alignment Quality (SAQ) %d for %s (%d Seqs).' % (haq.stat['SAQCyc'], seqlist.info['Name'],seqlist.seqNum()))
                        haq.singleSeqAQ(seqlist,focus=self._queryFocus())
                        seqnum.append(seqlist.seqNum())
                        self.log.printLog('#HAQ','SAQ%d: %d -> %d sequences.' % (haq.stat['SAQCyc'],seqnum[-2],seqnum[-1]))
                        if self.stat['DataLevel'] > 2:
                            if self.opt['UseAln'] == False: seqlist.tidyGaps(key='SAQX')
                            haq.saveHAQ(seqlist,filename='%s%s.saqx.%d.fas' % (self.getStr('ResDir'),self.info['Basefile'],haq.stat['SAQCyc']),key='SAQX')
                        if self.stat['DataLevel'] > 1:
                            if self.opt['UseAln'] == False: seqlist.tidyGaps(key='SAQ')
                            haq.saveHAQ(seqlist,filename='%s%s.saq.%d.fas' % (self.getStr('ResDir'),self.info['Basefile'],haq.stat['SAQCyc']),key='SAQ')
                        seqnum[-1] = seqlist.seqNum()
                        if seqnum[-2] != seqnum[-1]: continue   # Dumped seqs at SAQ - continue                            
                    ## <iii> ## PAQ
                    if self.opt['PAQ']:
                        haq.stat['PAQCyc'] += 1
                        self.log.printLog('#INF','Pairwise Alignment Quality (PAQ) %d for %s (%d Seqs).' % (haq.stat['PAQCyc'],seqlist.info['Name'],seqlist.seqNum()))
                        haq.pairwiseAQ(seqlist,focus=self._queryFocus())
                        seqnum[-1] = seqlist.seqNum()
                        self.log.printLog('#HAQ','PAQ%d: %d -> %d sequences.' % (haq.stat['PAQCyc'],seqnum[-2],seqnum[-1]))
                        if self.stat['DataLevel'] > 1:
                            if self.opt['UseAln'] == False: seqlist.tidyGaps(key='PAQ')
                            haq.saveHAQ(seqlist,filename='%s%s.paq.%d.fas' % (self.getStr('ResDir'),self.info['Basefile'],haq.stat['PAQCyc']),key='SAQ')
                    ## <iv> ## PostPAQ
                    while self.opt['PrePAQ'] and seqnum[-2] == seqnum[-1]:
                        if self.i() < 0 or rje.yesNo('*** No Sequences lost. Accept %d seqs (%d Groups)?' % (seqlist.seqNum(),tree.groupNum())): break
                        elif rje.yesNo('\nSave Sequences in %s.allfam.fas?' % self.info['Basefile']):
                            seqlist.saveFasta(seqfile='%s.allfam.fas' % self.info['Basefile'])  # &SaveDataSeq("allfam.", "Post PAQ before manual cutdown")
                            self.prePAQ(review=True)
                            seqnum.append(seqlist.seqNum())
                if self.opt['UseAln'] == False:
                    seqlist.tidyGaps()
                if self.stat['DataLevel'] > 0:
                    if self.opt['UseAln'] == False and self.opt['PAQ']: seqlist.tidyGaps(key='PAQ')
                    haq.saveHAQ(seqlist,filename='%s%s.haq.fas' % (self.getStr('ResDir'),self.info['Basefile']),key='PAQ')
                self.verbose(0,1,'\nFinished Pairwise Alignment Quality (%d Cycles):\n => %d of %d sequences retained.' % (haq.stat['PAQCyc'],seqnum[-1],seqnum[1]),1)
                seqlist.saveFasta()
                    
            ### <3> ### ES: Establishment of Subfamilies (Based on AFD)
            _stage = '<3> ES'
            if (self.opt['ES'] or (self.opt['HAQ'] and self.opt['AC'])) and seqlist.seqNum() > 1:
                self.log.printLog('#INF','Establishment of Subfamilies for %s (%d Seqs, Group=%s).' % (seqlist.info['Name'],seqlist.seqNum(),tree.info['Grouping']))
                if seqlist.opt['Aligned'] == False:
                    seqlist = self._aln(seqlist)
                seqnum = 0
                while seqlist.seqNum() != seqnum:
                    self._makeTree(boot=1000,kimura=True) #!# Add more tree options (different drawing methods)
                    tree.treeGroup(callmenu=True); self.setStat({'Interactive':tree.stat['Interactive']})
                    tree._reviewGroups(interactive=self.stat['Interactive'])  
                    if tree._orphanCount() > 0 and tree.opt['Orphans'] == False:
                        if self.stat['Interactive'] < 0 or rje.yesNo('Remove Sequences not in a subgroup?'):
                            tree._purgeOrphans()
                    tree._saveGroups(filename='%s.grp' % self.info['Basefile'],groupnames=True)
                    tree._reorderSeqToGroups()
                    seqnum = seqlist.seqNum()
                seqlist.saveFasta()
                tree.info['Basefile'] = self.info['Basefile'] 
                try: tree.saveTrees()
                except:
                    tree.saveTree(filename='%s.nsf' % self.info['Basefile'],type='nsf',seqnum=False,seqname='long',blen='Length',bootstraps='boot',multiline=True)
                    tree.textTree(filename='%s.tree.txt' % self.info['Basefile'],seqnum=True,seqname='short',nodename='short',showboot=True,showlen='branch',scale=4,spacer=1,compress=False)
                
            ### <4> ### AC: Ancestor Construction (GASP)
            _stage = '<4> AC (GASP)'
            if self.opt['AC'] and self.obj['SeqList'].seqNum() > 2:
                self.verbose(0,2,'\n%s' % self.obj['GASP'].details(),1)
                if self.stat['Interactive'] > 0 and not rje.yesNo('Use these parameters?'):
                    self.obj['GASP'].edit()
                if not self.opt['ES'] and not self.opt['HAQ']:
                    if rje.exists('%s.nsf' % self.info['Basefile']):
                        tree.loadTree(file='%s.nsf' % self.info['Basefile'],seqlist=seqlist)
                    elif rje.exists('%s.nwk' % self.info['Basefile']):
                        tree.loadTree(file='%s.nwk' % self.info['Basefile'],seqlist=seqlist)
                    else:
                        self.errorLog('Cannot finf *.nwk or *.nsf: no GASP analysis',printerror=False)
                        self.setBool({'AC':False})
                if self.getBool('AC') and self.obj['SeqList'].seqNum() > 2: self.obj['GASP'].gasp()

            ### <5> ### Final Data Adjustments
            _stage = '<5> Final Data Adjustments'
            #!# Add option to save full and cut-down #!#
            #!# NOTE: this is not really working
            if self.dev() and self.opt['ES'] and self.stat['Interactive'] > -1 and self.obj['SeqList'].seqNum() > 2 and rje.yesNo('Save *.full.* backup and reduce groups for further analysis (BADASP etc.)?',default='N'):
                #!# Add output of this to self.getStr('ResDir')
                seqlist.saveFasta(seqfile='%s.full.fas' % self.info['Basefile'])
                tree.info['Basefile'] = '%s.full' % self.info['Basefile'] 
                try: tree.saveTrees()
                except:
                    tree.saveTree(filename='%s.full.nsf' % self.info['Basefile'],type='nsf',seqnum=False,seqname='long',blen='Length',bootstraps='boot',multiline=True)
                    tree.textTree(filename='%s.full.tree.txt' % self.info['Basefile'],seqnum=True,seqname='short',nodename='short',showboot=True,showlen='branch',scale=4,spacer=1,compress=False)
                tree._saveGroups(filename='%s.full.grp' % self.info['Basefile'],groupnames=True)
                if self.opt['AC']:
                    tree.ancSeqOut(file='%s.full.anc.fas' % self.info['Basefile'],ordered=self.obj['GASP'].opt['Ordered'])
                tree.treeGroup(callmenu=True); self.setStat({'Interactive':tree.stat['Interactive']})
                tree._reviewGroups(interactive=self.stat['Interactive'])  
                tree._saveGroups(filename='%s.grp' % self.info['Basefile'],groupnames=True)
                tree._reorderSeqToGroups()
                seqnum = seqlist.seqNum()
                seqlist.saveFasta()
                tree.info['Basefile'] = self.info['Basefile'] 
                try: tree.saveTrees()
                except:
                    tree.saveTree(filename='%s.nsf' % self.info['Basefile'],type='nsf',seqnum=False,seqname='long',blen='Length',bootstraps='boot',multiline=True)
                    tree.textTree(filename='%s.tree.txt' % self.info['Basefile'],seqnum=True,seqname='short',nodename='short',showboot=True,showlen='branch',scale=4,spacer=1,compress=False)
                if self.opt['AC']:
                    tree.ancSeqOut(file='%s.anc.fas' % self.info['Basefile'],ordered=self.obj['GASP'].opt['Ordered'])
            ### Adjust other data output ###
            if self.stat['DataLevel'] > 1:
                seqlist.saveScanSeq(seqfile='%s%s.scanseq.fas' % (self.getStr('ResDir'),self.baseFile()))
                seqlist.degapSeq(log=False)
                seqlist.saveFasta(seqfile='%s%s.degap.fas' % (self.getStr('ResDir'),seqlist.info['Basefile']),name='AccNum')
            for ext in ['psq','psi','psd','pin','phr']:
                if rje.checkForFile('%s.%s' % (self.info['Basefile'],ext)):
                    os.unlink('%s.%s' % (self.info['Basefile'],ext))
            if rje.checkForFile('%s_cw.phb' % self.info['Basefile'].lower()):
                os.unlink('%s_cw.phb' % self.info['Basefile'].lower())
            if rje.checkForFile('%s.phb' % self.info['Basefile']):
                os.rename('%s.phb' % self.info['Basefile'],'%s%s.phb' % (self.getStr('ResDir'),self.info['Basefile']))
            if self.stat['DataLevel'] > 0 and self.stat['Interactive'] > 0 and not rje.yesNo('Keep all output files?'):
                filelist = glob.glob('%s.*' % self.info['Basefile']) + glob.glob('%s%s.*' % (self.getStr('ResDir'),self.baseFile()))
                for file in filelist:
                    if not rje.yesNo('Keep %s?' % file):
                        os.unlink(file)

            ### <6> ### End
            _stage = '<6> End'
            email = 'richard.edwards@unsw.edu.au'
            self.verbose(0,0,'\n\nThank you for using HAQESAC. As always, it is recommended that you check your alignments manually for errors. Please send any comments or suggestions to %s.' % email,2)

        except SystemExit: self.printLog('#EXIT','Main HAQESAC run terminated.')
        except: self.errorLog('Fatal problem with %s during main HAQESAC Run().' % _stage,printerror=True, quitchoice=True)
#########################################################################################################################
    def _queryFocus(self):  ### Returns start and end position in alignment of focus for query
        '''
        Returns start and end position in alignment of focus for query [X,Y].
        '''
        try:
            _focus = [0,-1]
            query = self.obj['SeqList'].obj['QuerySeq']
            qfocus = rje.matchExp('(\d+),(\d+)',self.info['QryRegion'])
            if query and qfocus:
                (q,r) = (0,0)
                if int(qfocus[0]) > int(qfocus[1]):
                    self.log.errorLog('Beginning of Query focus (%s) exceeds end (%s)! Will use all Query.' % (qfocus[0],qfocus[1]),printerror=False,quitchoice=True)
                    return _focus
                elif int(qfocus[0]) > query.aaLen():
                    self.log.errorLog('Beginning of Query focus (%s) exceeds length of Query %s (%d aa)! Will use all Query.' % (qfocus[0],query.shortName(),query.aaLen()),printerror=False,quitchoice=True)
                    return _focus
                elif int(qfocus[1]) > query.aaLen():
                    self.log.errorLog('End of Query focus (%s) exceeds length of Query %s (%d aa)! Will use end.' % (qfocus[1],query.shortName(),query.aaLen()),printerror=False,quitchoice=True)
                    return _focus
                while q < int(qfocus[0]) - 1 and r < query.seqLen():
                    #!# NOTE: This is now 1-L numbering.
                    if query.info['Sequence'][r] != '-':
                        q += 1
                        _focus[0] = r
                    r += 1
                while q < int(qfocus[1]) and r < query.seqLen():
                    if query.info['Sequence'][r] != '-':
                        q += 1
                        _focus[1] = r 
                    r += 1
            return _focus
        except:
            self.log.errorLog('Major Problem with _queryFocus()')
            return [0,-1]
#########################################################################################################################
    def cleanUp(self,seqlist=None):  ### Cleans up data and does prelim cleanup of low %ID sequences                #V1.9
        '''
        Cleans up data and does prelim cleanup of low %ID sequences.
        >> seqlist:SeqList object (else use self.obj['SeqList'])
        << returns cleaned up seqlist
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqlist: seqlist = self.obj['SeqList']
            seqnum = [seqlist.seqNum()]
            seqlist.printLog('#INF','Cleaning up data for %s (%d Sequences).' % (seqlist.info['Name'],seqlist.seqNum()))

            ### ~ [2] Aln Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['UseAln']: seqlist.makeNR(text='Cleanup (Aligned).',blast=0,pw_aln=False)
            else: seqlist.degapSeq()
            self.cmd_list += ['accnr=F','seqnr=F']
            seqlist.opt['AccNR'] = False
            seqlist.opt['SeqNR'] = False            
                
            ### ~ [3] BLAST Reduction & GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            query = seqlist.obj['QuerySeq']
            if self.stat['BlastCut'] == 0:
                self.stat['BlastCut'] = seqlist.seqNum()
            if query and (self.stat['BlastCut'] > 0 or self.opt['GABLAM'] or self.stat['QCover'] > 0):  # BLAST reduction
                ## ~ [3a] BlastCut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                prenum = seqlist.seqNum()
                # Setup BLAST
                seqlist.saveFasta()
                seqlist.saveFasta(seqs=[query],seqfile='%s.qry' % query.info['AccNum'])
                blastcmd = self.cmd_list + ['blastp=blastp','blastd=%s' % seqlist.info['Name'],'formatdb=T','blastv=%d' % self.stat['BlastCut'],'i=0','blasti=%s.qry' % query.info['AccNum']]
                if self.opt['GABLAM'] or self.stat['QCover'] > 0:   ### Need alignments
                    blastcmd += ['blastb=%d' % self.stat['BlastCut'],'blastf=F']
                else:
                    blastcmd += ['blastb=0']
                # Perform BLAST
                qblast = rje_blast.blastObj(log=self.log,cmd_list=blastcmd,type='Dev')
                qblast.setStr({'Name':'%s.blast' % query.info['AccNum']})
                qblast.blast(cleandb=True)
                os.unlink(qblast.getStr('InFile'))
                if not qblast.readBLAST(unlink=True,clear=True): raise ValueError
                # Read in hit sequences (in Rank order)
                if qblast.oldBLAST():
                    search = qblast.search[0]
                    hitdic = search.hitSeq(seqlist)
                    qbseq = []
                    for hit in search.hit:
                        qbseq.append(hitdic[hit])
                else:
                    qbseq = qblast.hitToSeq(seqlist,asdict=False)
                    hitdic = seqlist.seqNameDic(proglog=False)
                # Check for missing sequences, remove and reorder within seqlist
                for seq in seqlist.seq[0:]:
                    if seq not in qbseq:
                        if len(qbseq) < self.stat['BlastCut']:
                            seqlist.removeSeq(text='Filtered by BLAST: Hit < %e vs %s.' % (qblast.getNum('E-Value'),query.shortName()),seq=seq,checkAln=False)
                        else:
                            seqlist.removeSeq(text='Filtered by BLAST: not in top %d hits vs %s.' % (self.stat['BlastCut'],query.shortName()),seq=seq,checkAln=False)
                seqlist.seq = qbseq[0:]
                if prenum > seqlist.seqNum():
                    self.log.printLog('#SEQ','BLAST Cut-off Filter vs %s: %d sequences remain.' % (query.shortName(), seqlist.seqNum()))
                    prenum = seqlist.seqNum()
                ## ~ [3b] Perform GABLAM vs Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if qblast.oldBLAST():
                    qlen = search.stat['Length']
                    for hit in search.hit[0:]:
                        if hitdic[hit] == query: continue   # Self
                        if not hitdic[hit]:
                            self.errorLog('BLAST hit "%s" not found in original sequences for some reason!' % hit.info['Name'],printerror=False)
                            continue        # Already removed somehow?
                        gdict = hit.globalFromLocal(qlen)
                        # Coverage #
                        if (float(100 * gdict['Query']['GABLAMO Len']) / qlen) < self.stat['QCover'] and (float(100 * gdict['Hit']['GABLAMO Len']) / hit.stat['Length']) < self.stat['QCover']:
                            seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% coverage vs Query.' % self.stat['QCover'],seq=hitdic[hit],checkAln=False)    
                        # GABLAMO #
                        elif self.opt['GABLAM'] and self.stat['QryID'] > 0 and self.opt['GABSim']:
                            if (float(100 * gdict['Query']['GABLAMO Sim']) / qlen) < self.stat['QryID'] and (float(100 * gdict['Hit']['GABLAMO Sim']) / hit.stat['Length']) < self.stat['QryID']:
                                seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% similarity vs Query.' % self.stat['QryID'],seq=hitdic[hit],checkAln=False)    
                        elif self.opt['GABLAM'] and self.stat['QryID'] > 0:
                            if (float(100 * gdict['Query']['GABLAMO ID']) / qlen) < self.stat['QryID'] and (float(100 * gdict['Hit']['GABLAMO ID']) / hit.stat['Length']) < self.stat['QryID']:
                                seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% identity vs Query.' % self.stat['QryID'],seq=hitdic[hit],checkAln=False)
                else:
                    qlen = query.aaLen()
                    for hseq in qbseq:
                        if hseq == query: continue  # Self
                        hit = hseq.shortName()
                        hlen = hseq.aaLen()
                        gdict = qblast.gablamData(hit,query.shortName())
                        # Coverage #
                        if (float(100 * gdict['Query']['GABLAMO Len']) / qlen) < self.stat['QCover'] and (float(100 * gdict['Hit']['GABLAMO Len']) / hlen) < self.stat['QCover']:
                            seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% coverage vs Query.' % self.stat['QCover'],seq=hseq,checkAln=False)    
                        # GABLAMO #
                        elif self.opt['GABLAM'] and self.stat['QryID'] > 0 and self.opt['GABSim']:
                            if (float(100 * gdict['Query']['GABLAMO Sim']) / qlen) < self.stat['QryID'] and (float(100 * gdict['Hit']['GABLAMO Sim']) / hlen) < self.stat['QryID']:
                                seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% similarity vs Query.' % self.stat['QryID'],seq=hseq,checkAln=False)    
                        elif self.opt['GABLAM'] and self.stat['QryID'] > 0:
                            if (float(100 * gdict['Query']['GABLAMO ID']) / qlen) < self.stat['QryID'] and (float(100 * gdict['Hit']['GABLAMO ID']) / hlen) < self.stat['QryID']:
                                seqlist.removeSeq(text='Ordered BLAST Local hits < %.1f%% identity vs Query.' % self.stat['QryID'],seq=hseq,checkAln=False)
                if prenum > seqlist.seqNum():
                    self.printLog('#SEQ','GABLAMO Filter vs %s: %d sequences remain.' % (query.shortName(), seqlist.seqNum()))

            ### ~ [4] Paired %ID - PairID and QryID - and GapFilter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self._IDGapFilter(seqlist)
            
            ### ~ [5] Finish and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist.printLog('#INF','Cleanup of %s complete: %d of %d Sequences remain.' % (seqlist.info['Name'],seqlist.seqNum(),seqnum[0]))
            seqlist.saveFasta()
            self.verbose(0,1,'',0)  # Pause if i=1
            return seqlist            
        except: self.errorLog('Major Problem with Data Cleanup.',quitchoice=True)
#########################################################################################################################
    def _IDGapFilter(self,seqlist):     ### Performs loops of pairID and GapFilter (if appropriate) until no seqs lost
        '''
        Performs loops of pairID and GapFilter (if appropriate) until no seqs lost.
        >> seqlist:SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prevx = seqlist.seqNum()    # Store current sequence number to assess change
            if prevx <2: return         # No need/point in comparing if only one sequence!
            seqnum = [0,prevx]          # Set up sequence count for loops
            qry = seqlist.obj['QuerySeq']   # Set up query sequence for focus of comparison
            ### ~ [1] ~ Loop through comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while seqnum[-2] != seqnum[-1]:
                self._pairedIDFilter(seqlist)
                if self.opt['UseAln'] == False:
                    seqlist.info['Name'] = '%s.fas' % self.info['Basefile']
                    seqlist = self._aln(seqlist)
                seqlist.gapSeqFilter(relative='neighbour')  # Relative to neighbour (Aligned)
                if qry and qry not in seqlist.seqs():
                    seqlist.seq = [qry] + seqlist.seq[0:]
                    self.printLog('#GAP','Reinstated very gappy Query %s! (Consider changing maxgap=X.)' % qry.shortName())
                seqnum.append(seqlist.seqNum())
            if prevx != seqlist.seqNum(): self.printLog('#SEQ','%%ID/Gap Filter: %d -> %d sequences.' % (prevx, seqlist.seqNum()))
        except: self.errorLog('Major Problem with _IDGapFilter')
#########################################################################################################################
    def _pairedIDFilter(self,seqlist):  ### Populated DisMatrix Object and removes bad sequences
        '''
        Populated DisMatrix Object and removes bad sequences.
        >> seqlist=SeqList Object
        '''
        try:
            ### <0> ### Setup
            _stage = '0'
            if self.stat['QryID'] <= 0 and self.stat['PairID'] <= 0:
                return
            query = seqlist.obj['QuerySeq']
            idkey = 'PWAln ID'
            blastcmd = self.cmd_list + ['blastp=blastp','blastd=%s' % seqlist.info['Name'],'formatdb=T','i=-1','v=-1']
            if self.opt['GABLAM']:   ### Need alignments
                blastcmd += ['blastb=%d' % self.stat['BlastCut'],'blastf=F']
            else:
                blastcmd += ['blastb=0']
            if self.opt['UseAln']:
                idkey = 'MSA ID'
                blastseq = copy.deepcopy(seqlist)
                blastseq.info['Name'] = '%s.degap.fas' % self.info['Basefile']
                blastseq.degapSeq(log=False)
                blastseq.saveFasta()
                blastcmd.append('blastd=%s' % blastseq.info['Name'])
            else:
                seqlist.degapSeq(log=False)
                seqlist.saveFasta()
            seqlist.addMatrix(idkey,sym=False)
            _stage = '0b'
            for seq in seqlist.seq:
                #print seqlist.seq.index(seq), seqlist.seqNum(), seq
                if 'PairID' not in seq.obj.keys():
                    seq.obj['PairID'] = None
            
            ### <1> ### Query
            _stage = '1'
            if self.stat['QryID'] > 0 and not self.opt['GABLAM'] and not self.obj['HAQ'].opt['NoQuery']:
                _stage = '1a'
                query.obj['PairID'] = query
                for seq in seqlist.seq[0:]:
                    _stage = '1b'
                    if seqlist.getDis(seq,query,key=idkey) < self.stat['QryID']:    # Bad!
                        _stage = '1c'
                        seqlist.removeSeq(text='< %.1f%%ID with Query (%s)' % (self.stat['QryID'],query.shortName()),seq=seq,checkAln=False)
                    elif seqlist.getDis(seq,query,key=idkey) >= self.stat['PairID'] and seq.info['SpecCode'] != query.info['SpecCode']: # Good!
                        _stage = '1d'
                        self.verbose(0,2,'%d of %d: %s >= %.1f%%ID with %s.' % (seqlist.seq.index(seq),seqlist.seqNum(),seq.shortName(),self.stat['PairID'],query.shortName()),1)
                        seq.obj['PairID'] = query
                _stage = '1c'
                seqlist.saveFasta()

            ### <2> ### Others
            _stage = '2'
            if self.stat['PairID'] <= 0:
                return
            blast = rje_blast.blastObj(log=self.log,cmd_list=blastcmd,type='Dev')
            cycle = True
            while cycle:
                _stage = '2a'
                cycle = False
                for seq in seqlist.seq[0:]:
                    _stage = '2b'
                    self.verbose(2,2,'%s: PairID=%s' % (seq.shortName(),seq.obj['PairID']),1)
                    if seq.obj['PairID'] and seq.obj['PairID'] in seqlist.seq:    # OK for now
                        continue
                    seq.obj['PairID'] = None
                    seqlist.saveFasta(seqs=[seq],seqfile='%s.fas' % seq.info['AccNum'],log=False)
                    blast.setStr({'InFile':'%s.fas' % seq.info['AccNum'],'Name':'%s.blast' % seq.info['AccNum']})
                    blast.blast(cleandb=True); os.unlink(blast.getStr('InFile'))
                    if not blast.readBLAST(unlink=True,clear=True): raise ValueError

                    if blast.oldBLAST():
                        search = blast.search[-1]
                        hitdic = search.hitSeq(seqlist)
                        for hit in search.hit:
                            # Check whether comparison needed #
                            _stage = '2c'
                            otherseq = hitdic[hit]
                            if otherseq == None or otherseq == query or otherseq not in seqlist.seq or seq.info['SpecCode'] == otherseq.info['SpecCode']:
                                continue
                            if seq.obj['PairID'] in seqlist.seq and otherseq.obj['PairID'] in seqlist.seq:   # Next hit
                                continue
                            # Calculate pairwise comparisons #
                            if self.opt['GABLAM']:
                                qlen = seq.aaLen()
                                gdict = hit.globalFromLocal(qlen)
                                if self.opt['GABSim']:
                                    pairid = (float(100 * gdict['Hit']['GABLAMO Sim']) / hit.stat['Length'],float(100 * gdict['Query']['GABLAMO Sim']) / qlen)
                                else:
                                    pairid = (float(100 * gdict['Hit']['GABLAMO ID']) / hit.stat['Length'],float(100 * gdict['Query']['GABLAMO ID']) / qlen)
                            else:
                                pairid = (seqlist.getDis(otherseq,seq,key=idkey),seqlist.getDis(seq,otherseq,key=idkey))
                            # Process and keep if OK #
                            self.verbose(3,2,'Hit %d = %s: %.1f%%ID <> %.1fID' % (search.hit.index(hit),otherseq.shortName(),pairid[0],pairid[1]),1)
                            if seq.obj['PairID'] and pairid[0] < self.stat['PairID']:   # Next search
                                break
                            if otherseq.obj['PairID'] == None and pairid[0] >= self.stat['PairID']:
                                cycle = True
                                otherseq.obj['PairID'] = seq
                                self.verbose(0,2,'%d of %d: %s >= %.1f%%ID with %s.' % (seqlist.seq.index(otherseq),seqlist.seqNum(),otherseq.shortName(),self.stat['PairID'],seq.shortName()),1)
                            if seq.obj['PairID'] == None and pairid[1] >= self.stat['PairID']:
                                cycle = True
                                seq.obj['PairID'] = otherseq
                                self.verbose(0,2,'%d of %d: %s >= %.1f%%ID with %s.' % (seqlist.seq.index(seq),seqlist.seqNum(),seq.shortName(),self.stat['PairID'],otherseq.shortName()),1)
                                #break
                    else:   #!# Needs a bit of tidying, I think! #!#
                        qlen = seq.aaLen()
                        hitlist = blast.queryHits(query.shortName())
                        hitdic = blast.hitToSeq(seqlist,hitlist)
                        for hit in hitlist:
                            # Check whether comparison needed #
                            _stage = '2c'
                            otherseq = hitdic[hit]
                            if otherseq == None or otherseq == query or otherseq not in seqlist.seq or seq.info['SpecCode'] == otherseq.info['SpecCode']:
                                continue
                            if seq.obj['PairID'] in seqlist.seq and otherseq.obj['PairID'] in seqlist.seq:   # Next hit
                                continue
                            # Calculate pairwise comparisons #
                            if self.opt['GABLAM']:
                                hlen = otherseq.aaLen()
                                gdict = blast.gablamData(hit,query.shortName())
                                if self.opt['GABSim']:
                                    pairid = (float(100 * gdict['Hit']['GABLAMO Sim']) / hlen,float(100 * gdict['Query']['GABLAMO Sim']) / qlen)
                                else:
                                    pairid = (float(100 * gdict['Hit']['GABLAMO ID']) / hlen,float(100 * gdict['Query']['GABLAMO ID']) / qlen)
                            else:
                                pairid = (seqlist.getDis(otherseq,seq,key=idkey),seqlist.getDis(seq,otherseq,key=idkey))
                            # Process and keep if OK #
                            self.verbose(3,2,'Hit %d = %s: %.1f%%ID <> %.1fID' % (hitlist.index(hit),otherseq.shortName(),pairid[0],pairid[1]),1)
                            if seq.obj['PairID'] and pairid[0] < self.stat['PairID']:   # Next search
                                break
                            if otherseq.obj['PairID'] == None and pairid[0] >= self.stat['PairID']:
                                cycle = True
                                otherseq.obj['PairID'] = seq
                                self.verbose(0,2,'%d of %d: %s >= %.1f%%ID with %s.' % (seqlist.seq.index(otherseq),seqlist.seqNum(),otherseq.shortName(),self.stat['PairID'],seq.shortName()),1)
                            if seq.obj['PairID'] == None and pairid[1] >= self.stat['PairID']:
                                cycle = True
                                seq.obj['PairID'] = otherseq
                                self.verbose(0,2,'%d of %d: %s >= %.1f%%ID with %s.' % (seqlist.seq.index(seq),seqlist.seqNum(),seq.shortName(),self.stat['PairID'],otherseq.shortName()),1)

                    _stage = '2d'
                    #print seq.shortName(), seq.obj, self.info
                    if seq.obj['PairID'] == None and seq.info['SpecCode'] != 'UNK':
                        seqlist.removeSeq(text='< %.1f%%ID with any other (non %s) sequence.' % (self.stat['PairID'],seq.info['SpecCode']),seq=seq,checkAln=False)
        except:
           self.log.errorLog('Major Problem with _pairedIDFilter <%s>.' % _stage, printerror=True, quitchoice=True)
#########################################################################################################################
    def _makeTree(self,boot=1000,kimura=True):    ### Makes a tree and populates self.obj['Tree']
        '''
        Makes a tree and populates self.obj['Tree'].
        >> type:str = type of tree to make. #!# Add more tree options (different drawing methods)
        >> boot:int = number of bootstraps 
        >> kimura:boolean = whether to correct for multiple hits
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList']
            tree = self.obj['Tree']
            tree.setInt({'Bootstraps':boot})
            tree.setOpt({'Kimura':kimura})
            ### ~ [2] ~ Make tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist.seqNum() > 2: tree.makeTreeMenu(interactiveformenu=1,force=True,make_seq=seqlist)
            else:   # Assume there are two sequences, otherwise raise and Error and fix code calling this!
                idkey = 'MSA ID'
                seqlist.addMatrix(idkey,sym=False)
                seq = seqlist.seqs()
                branchlen = seqlist.getDis(seq[0],seq[1],key='MSA ID') / 200.0
                self.printLog('#TREE','Made tree of two sequences from MSA %ID')
                nsftree = '(%s:%s,%s:%s);' % (seq[0].shortName(),branchlen,seq[1].shortName(),branchlen)
                tree.buildTree(nsftree,seqlist,type='nsf',postprocess=False)
        except: self.errorLog('Major Problem with makeTree().', printerror=True, quitchoice=True)
#########################################################################################################################
    def prePAQ(self,review=False):    ### Prelim to Pairwise Alignment Assessment for single sequences relative to Query
        '''
        Prelim to Pairwise Alignment Assessment for single sequences relative to Query.
        >> review:boolean = whether to force review. [False]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tree = self.obj['Tree']
            bootnum = tree.stat['Bootstraps']
            if bootnum > 100: bootnum = 100
            ### ~ [2] ~ Make tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self._makeTree(boot=bootnum,kimura=False)
            ### ~ [3] ~ Manual pre-PAQ gubbins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['HAQ'].stat['PAQCyc'] == 0 or tree.info['Grouping'] == 'man':
                tree.treeGroup(callmenu=True); self.setStat({'Interactive':tree.stat['Interactive']})
                review = True
            else: tree.treeGroup()
            if review or tree._checkGroups() == False:     #!# or manrev - add manrev opt that can be turned off!
                tree._checkGroupNames()
                tree._reviewGroups(interactive=self.stat['Interactive'])  #!# Look at interactive thing!
            tree._reorderSeqToGroups()        
            if not tree.opt['Orphans'] or (self.stat['Interactive'] > 0 and not rje.yesNo('Keep Sequences not in a subgroup?')): tree._purgeOrphans
        except SystemExit: raise
        except: self.errorLog('Major Problem with prePAQ().',quitchoice=True)
#########################################################################################################################
### End of SECTION II: HAQESAC Class                                                                                    #
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
    try:HAQESAC(log=mainlog, cmd_list=cmd_list).run(setobjects=True)

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: pass     # Killed by MultiHAQ (I presume!)
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
