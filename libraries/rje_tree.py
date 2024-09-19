#!/usr/bin/python

# rje_tree.py - Phylogenetic Tree module
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 99 Wilton Road, Southampton SO15 5JH, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_tree
Description:  Phylogenetic Tree Module
Version:      2.17.1
Last Edit:    11/07/19
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    Reads in, edits and outputs phylogenetic trees. Executes duplication and subfamily determination. More details available
    in documentation for HAQESAC, GASP and BADASP at http://www.bioinformatics.rcsi.ie/~redwards/

General Commands:
    nsfin=FILE      : load NSF tree from FILE
    phbin=FILE      : load ClustalW Format *.phb NSF tree from FILE
    seqin=FILE      : load sequence list from FILE (not compatible with useanc)
    disin=FILE      : load distance matrix from FILE (Phylip format for use with distance matrix methods) [None]
    useanc=FILE     : load sequences from ancestral sequence FILE (not compatible with seqin)
    deflen=X        : Default length for branches when no lengths given [0.1] (or 0.1 x longest branch)
    *Note that in the case of conflicts (e.g. seqin=FILE1 useanc=FILE2), the latter will be used.*
    autoload=T/F    : Whether to automatically load sequences upon initiating object [True]

Rooting Commands:
    root=X  : Rooting of tree (rje_tree.py):
        - mid = midpoint root tree.
        - ran = random branch.
        - ranwt = random branch, weighted by branch lengths.
        - man = always ask for rooting options (unless i<0).
        - none = unrooted tree
        - FILE = with seqs in FILE as outgroup. (Any option other than above)
    rootbuffer=X    : Min. distance from node for root placement (percentage of branch length)[0.1]

Grouping/Subfamily Commands:
    bootcut=X   : cut-off percentage of tree bootstraps for grouping.
    mfs=X       : minimum family size [3]
    fam=X       : minimum number of families (If 0, no subfam grouping) [0]
    orphan=T/F  : Whether orphans sequences (not in subfam) allowed. [True]
    allowvar=T/F: Allow variants of same species within a group. [False]
    qryvar=T/F  : Keep variants of query species within a group (over-rides allowvar=F). [False]
    keepvar=LIST: Automatically keep variants for listed species (over-rides allowvar=F). []
    groupspec=X : Species for duplication grouping [None]
    specdup=X   : Minimum number of different species in clade to be identified as a duplication [1]
    9spec=T/F   : Whether to treat 9XXXX species codes as actual species (generally higher taxa) [False]
    group=X     : Grouping of tree
        - man = manual grouping (unless i<0).
        - dup = duplication (all species unless groupspec specified).
        - qry = duplication with species of Query sequence (or Sequence 1) of treeseq
        - one = all sequences in one group
        - None = no group (case sensitive)
        - FILE = load groups from file

Tree Making Commands:
    cwtree=FILE     : Make a ClustalW NJ Tree from FILE (will save *.ph or *.phb) [None]
    kimura=T/F      : Whether to use Kimura correction for multiple hits [True]
    bootstraps=X    : Number of bootstraps [0]
    clustalw=CMD    : Path to CLUSTALW (and including) program [''] * Use forward slashes (/)
    fasttree=PATH   : Path to FastTree (and including) program ['']
    iqtree=FULLPATH : Path IQTree program including program ['iqtree']
    phylip=PATH     : Path to PHYLIP programs [''] * Use forward slashes (/)
    phyoptions=FILE : File containing extra Phylip tree-making options ('batch running') to use [None]
    protdist=FILE   : File containing extra Phylip PROTDIST options ('batch running') to use [None]
    maketree=X      : Program for making tree [None]
        - None = Do not make tree from sequences 
        - clustalw = ClustalW NJ method
        - neighbor = PHYLIP NJ method 
        - upgma    = PHYLIP UPGMA (neighbor) method 
        - fitch    = PHYLIP Fitch method 
        - kitsch   = PHYLIP Kitsch (clock) method 
        - protpars = PHYLIP MP method 
        - proml    = PHYLIP ML method
        - fasttree = Use FastTree
        - iqtree   = Use IQTree
        - PATH     = Alternatively, a path to a different tree program/script can be given. This should accept ClustalW parameters.

Tree Display/Saving Commands
    savetree=FILE   : Save a generated tree as FILE [seqin.maketree.nsf]
    savetype=X      : Format for generated tree file (nsf/nwk/text/r/png/bud/qspec/cairo/te/svg/html) [nwk]
    treeformats=LIST: List of output formats for generated trees [nwk]
    outnames=X      : 'short'/'long' names in output file [short]
    truncnames=X    : Truncate names to X characters (0 for no truncation) [123]
    branchlen=T/F   : Whether to use branch lengths in output tree [True]
    deflen=X        : Default branch length (when none given, also for tree scaling) [0.1]
    textscale=X     : Default scale for text trees (no. of characters per deflen distance) [4]
    seqnum=T/F      : Output sequence numbers (if making tree from sequences) [True]
    
Classes:
    Tree(rje.RJE_Object):     
        - Phylogenetic Tree class.
    Node(rje.RJE_Object): 
        - Individual nodes (internal and leaves) for Tree object.
    Branch(rje.RJE_Object):
        - Individual branches for Tree object.
        
Uses general modules: copy, os, random, re, string, sys, time
Uses RJE modules: rje, rje_ancseq, rje_seq, rje_tree_group
Other module needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, random, re, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_ancseq, rje_html, rje_seq, rje_svg, rje_tree_group
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Changed CladeList to lists inside tuple
    # 0.2 - Completely reworked with a much greater OO focus
    # 0.3 - No Out Object in Objects
    # 1.0 - Better
    # 1.1 - Added tree-drawing with ClustalW
    # 1.2 - Bug fixing.
    # 1.3 - Added separate maketree=PATH command to enable replacement of ClustalW tree drawing
    # 1.4 - Added toggling of full sequence description in TextTree/EditTree
    # 1.5 - Added ability to read in integer branch lengths
    # 1.6 - Added PHYLIP tree making
    # 1.7 - Modified text/log output. Added commandline scale option for text trees.
    # 1.8 - Updating reporting of redundant duplication finding.
    # 1.9 - Modified mapSeq() method to be more robust to different formats
    # 1.10- Fixed some bugs and had a minor tidy.
    # 2.0 - Updated code to be more inline with newer RJE modules. Fixed some bugs.
    # 2.1 - Added tree savetype
    # 2.2 - Added treeformats=LIST option to eventually replace savetype. Added te (TreeExplorer) format.
    # 2.3 - Added specdup=X   : Minimum number of different species in clade to be identified as a duplication [1]
    # 2.4 - Added fasttree generation of trees.
    # 2.5 - Added qryvar=T/F  : Allow variants of query species within a group (over-rides allowvar=F). [False]
    # 2.6 - Added PNG R variants.
    # 2.7 - Added Improved NSF Tree reading.
    # 2.8 - Added SVG and HTML Tree output.
    # 2.9 - Added NWK output (=NSF output with different extension for MEGA!)
    # 2.10- Added cleanup of *.r.csv file following R-based PNG generation.
    # 2.11.0 - Modified for standalone running as part of SeqSuite.
    # 2.11.1 - Tweaked QryVar interactivity.
    # 2.11.2 - Updated tree paths.
    # 2.12.0 - Added treeLen() method.
    # 2.13.0 - Updated PNG saving with R to use newer code.
    # 2.14.0 - Added cladeSpec().
    # 2.14.1 - Fixed clustalw2 makeTree issue.
    # 2.15.0 - Added IQTree.
    # 2.16.0 - 9spec=T/F   : Whether to treat 9XXXX species codes as actual species (generally higher taxa) [False]
    # 2.16.1 - Modified NSF reading to cope with extra information beyond the ";".
    # 2.17.0 - keepvar=LIST: Automatically keep variants for listed species (over-rides allowvar=F). []
    # 2.17.1 - Added forking to multithread IQ-TREE runs. Tweaked BUDAPEST R call.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [y] Interactive menus for module and major activities
    # [y] read in from files*
    # [y] re-root/unroot as desired* (midpoint/random/outgroup/manual)
    # [y] save to files
    # [y] output text tree
    # [y] prune (remove/compress leaves/clades)
    # [ ] Add option for combining sequences (into a consensus?) - manipulates treeseqs object?
    # [ ] .. Also replace clade with ancestral sequence (GASP)
    # [y] identify duplications
    # [y] establish subfamilies/groups
    # [y] .. group=X commandline option
    # [ ] .. Modify Load Groups to allow both pure and fuzzy mapping (allows for some homoplasy?)
    # [Y] make trees from sequences
    # [Y] - CW NJ tree (as HAQESAC)
    # [ ] - Add different formats of distance matrix for disin=FILE
    # [ ] - distance matrix methods in SeqList class
    # [ ] .. Add UPGMA generation from aligned sequences?
    # [y] Tidy and __doc__ all methods
    # [y] Load ancestral sequences
    # [Y] Add a scale commandline option for textTrees
    # [ ] Add amino acid probabilities to nodes (for anc seq reconstruction)
    # [ ] Add better Exception trapping and check use of interactive and verbose settings.
    # [ ] Add a treein=X argument - will determine format from extension
    # [ ] Improve handling of duplication calculation - store status in self.opt and change only if pruned or rooting changed
    # [ ] Upgrade AncMapSeq to match modified MapSeq
    # [ ] Give module complete facelift in line with other modules! (Use list and dict attribute dictionaries.)
    # [ ] Sort out problem of multiple dup finding
    # [ ] Add MrBayes for making trees
    # [ ] Phase out savetype in favour of treeformats
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_TREE', '2.17.1', 'July 2019', '2007')
    description = 'RJE Phylogenetic Tree Module'
    author = 'Dr Richard J. Edwards.'
    comments = []
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
        print('Problem during initial setup.')
        raise
#########################################################################################################################
treeformats = {'nsf':'Newick Standard Format (NSF)','text':'Plain text','r':'table for R PNG maker',
               'png':'Portable Network Graphic (PNG)','te':'TreeExplorer-compatible NSF',
               'nwk':'Newick Standard Format (MEGA compatible)',
               'te.nsf':'TreeExplorer-compatible NSF','qspec':'PNG with query species highlighted',
               'bud':'PNG graphic for BUDAPEST runs','Cairo':'PNG using R with Cairo library',
               'svg':'Support Vector Graphic (SVG)','html':'Support Vector Graphic (SVG) embedded in HTML'}
formatext = {'nsf':'nsf','text':'tree.txt','r':'r','png':'png','te':'te.nsf','bud':'png','cairo':'png','svg':'svg',
             'html':'htm','te.nsf':'te.nsf','qspec':'png','nwk':'nwk'}
#########################################################################################################################
def treeName(name,nospace=False):     ### Reformats name to be OK in NSF file
    '''Reformats name to be OK in NSF file.'''
    rep = name[0:]
    rep = re.sub(',',' -',rep)
    if nospace: rep = re.sub('\s','_',rep)
    rep = re.sub('[\(\[]','{',rep)
    rep = re.sub('[\)\]]','}',rep)
    rep = re.sub('[:;]','-',rep)
    return rep
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Tree Class                                                                                              #
#########################################################################################################################
class Tree(rje.RJE_Object):     ### Class for handling phylogenetic trees.
    '''
    Phylogenetic Tree class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of tree (usually filename)
    - Type = Method of construction (if known)
    - Rooting = Rooting strategy [man]
    - RootMethod = Method used to determine root
    - GroupSpecies = Species to be used for Grouping
    - Grouping = Method used to determine groups
    - ClustalW = Path to ClustalW
    - Phylip = Path to PHYLIP programs ['c:/bioware/phylip3.65/exe/'] * Use forward slashes (/)
    - PhyOptions = File containing extra Phylip options ('batch running') to use [None]
    - ProtDist = File containing extra Phylip PROTDIST ('batch running') to use [None]
    - MakeTree = Tree drawing program (ClustalW options)
    - CWTree = Whether to use ClustalW NJ to make tree in self.makeTree()
    - FastTree = Path to FastTree (and including) program [./FastTree]
    - DisIn = load distance matrix from FILE (for use with distance matrix methods) [None]
    - SaveTree = Save a generated tree as FILE [None]
    - SaveType = Format for generated tree file (nsf/text/r/png/bud/cairo) [nsf]
    - OutNames = 'short'/'long' names in output file [short]
    - IQTree = Path IQTree program including program ['iqtree']

               
    Opt:boolean
    - Rooted = Whether tree is rooted
    - ReRooted = Whether tree has had its root altered
    - Branchlengths = Whether tree has branchlengths
    - Bootstrapped = Whether tree has bootstraps
    - QueryGroup = Group using specified Query Species
    - Orphans = Whether 'orphan' sequences allowed
    - AllowVar = Allow variants of same species within a group. [False]
    - QryVar = Allow variants of query species within a group (over-rides allowvar=F). [False]
    - Kimura = Whether to use Kimura multiple hit correction
    - OutputBranchLen = Whether to use branch lengths in output tree [True]
    - AutoLoad = Whether to automatically load sequences upon initiating object [True]
    - SeqNum = Output sequence numbers for general make tree thing [True]
    - 9SPEC=T/F   : Whether to treat 9XXXX species codes as actual species (generally higher taxa) [False]

    Stat:numeric
    - DefLen = Default length for branches when no lengths given [0.1]
    - TextScale = Default scale for text trees (no. of characters per deflen distance) [4]
    - RootBuffer = Min. distance from node for root placement (percentage of branch length)[0.1]
    - SeqNum = Number of seqs (termini)
    - Bootstraps = Number of bootstraps (if known)
    - BootCut = cut-off percentage of tree bootstraps for grouping.
    - MinFamSize = minimum family size [2]
    - MinFamNum : minfamnum = 0   # minimum number of families (If 0, no subfam grouping)
    - SpecDup = Minimum number of different species in clade to be identified as a duplication [1]
    - TruncNames = Truncate names to X characters (0 for no truncation) [123]

    List:list
    - KeepVar=LIST: Automatically keep variants for listed species (over-rides allowvar=F). []
    - TreeFormats = List of output formats for generated trees [nsf]

    Dict:dictionary    

    Obj:RJE_Objects
    - SeqList = rje_seq.SeqList object
    - PAM = rje_pam.PamCtrl object

    Other:
    - node = List of Node Objects
    - branch = List of Branch Objects
    - subfam = List of Node Objects that specifiy subgroup clades

    Additional Commands:
    - nsfin=FILE = load NSF tree from FILE
    '''
    ### Additional Attributes    
    node = []   # List of Node Objects
    branch = [] # List of Branch Objects
    subfam = [] # List of Node Objects that form subfamily clades
#########################################################################################################################
    ### <0> ### Basic Attribute methods                                                                                 #
#########################################################################################################################
    def groupNum(self): return len(self.subfams())
    def nodeNum(self): return len(self.nodes())
    def branchNum(self): return len(self.branches())
#########################################################################################################################
    def nodes(self): return self.node       # At some point, make a self.list['Node']
    def branches(self): return self.branch  # At some point, make a self.list['Branch']
    def subfams(self): return self.subfam   # At some point, make a self.list['SubFam']
#########################################################################################################################
    def seqNum(self):
        if self.obj['SeqList']: return self.obj['SeqList'].seqNum()
        elif self.opt['Rooted']: return int((self.nodeNum() + 1.0) / 2)
        else: return int((self.nodeNum() + 2.0) / 2)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','Rooting','RootMethod','GroupSpecies','Grouping','ClustalW','Phylip','PhyOptions',
                         'ProtDist','CWTree','MakeTree','DisIn','SaveTree','OutNames','SaveType','FastTree','IQTree']
        self.statlist = ['DefLen','RootBuffer','SeqNum','Bootstraps','BootCut','MinFamSize','MinFamNum','TruncNames',
                         'TextScale','SpecDup']
        self.optlist = ['Rooted','ReRooted','Branchlengths','Bootstrapped','QueryGroup','Orphans','AllowVar','Kimura',
                        'OutputBranchLen','AutoLoad','SeqNum','QryVar','9SPEC']
        self.listlist = ['KeepVar','TreeFormats']
        self.dictlist = []
        self.objlist = ['SeqList','PAM']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Type':'Unknown','Rooting':'man','ClustalW':'',
                      'Phylip':'','OutNames':'short','SaveType':'nwk',
                      'FastTree':'','IQTree':'iqtree'})
        self.setStat({'DefLen':0.1,'TextScale':4,'RootBuffer':0.1,'BootCut':0.7,'MinFamSize':3,'TruncNames':123,'SpecDup':1})
        self.setOpt({'Orphans':True,'Kimura':True,'OutputBranchLen':True,'AutoLoad':True,'SeqNum':False,'9SPEC':False})
        self.list['TreeFormats'] = ['nwk']
        self.obj['SeqList'] = None
        ### Other Attributes ###
        self.node = []   # List of Node Objects
        self.branch = [] # List of Branch Objects
        self.subfam = [] # List of Node Objects specifiying subfamilies
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [0] ~ Setup Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqfile = None
        ancseq = None
        ### ~ [1] ~ Read commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:## ~ [1a] ~ General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self._generalCmd(cmd)
                self._forkCmd(cmd)  # Delete if no forking
                ## ~ [1b] Normal Class Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self._cmdRead(cmd,type='info',att='Rooting',arg='root')
                self._cmdReadList(cmd,'stat',['DefLen','TextScale','RootBuffer','BootCut'])
                self._cmdReadList(cmd,'int',['MinFamSize','MinFamNum','Bootstraps','TruncNames','SpecDup'])
                self._cmdRead(cmd,type='int',att='MinFamSize',arg='mfs')
                self._cmdRead(cmd,type='int',att='MinFamNum',arg='fam')
                self._cmdRead(cmd,type='opt',att='Orphans',arg='orphan')
                self._cmdReadList(cmd,'opt',['Orphans','AllowVar','Kimura','AutoLoad','SeqNum','QryVar','9SPEC'])
                self._cmdRead(cmd,type='info',att='Grouping',arg='group')
                self._cmdRead(cmd,type='info',att='GroupSpecies',arg='groupspec')
                self._cmdReadList(cmd,'file',['CWTree','ClustalW','PhyOptions','ProtDist','MakeTree','DisIn','SaveTree','FastTree','IQTree'])
                self._cmdRead(cmd,type='path',att='Phylip')
                self._cmdRead(cmd,type='int',att='Bootstraps',arg='bootstrap')
                self._cmdReadList(cmd,'info',['OutNames','SaveType'])
                self._cmdRead(cmd,type='opt',att='OutputBranchLen',arg='branchlen')
                self._cmdReadList(cmd,'list',['KeepVar','TreeFormats'])
                ## ~ [1c] ~ Special ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if cmd.find('seqin=') == 0:
                    seqfile = cmd[len('seqin='):]
                    ancseq = None
                if cmd.find('useanc=') == 0:
                    ancseq = cmd[len('useanc='):]
                    seqfile = None
                self._cmdRead(cmd,type='file',att='Name',arg='nsfin')
                self._cmdRead(cmd,type='file',att='Name',arg='phbin')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [2] ~ Adjustments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.info['SaveType'] = self.info['SaveType'].lower()
        self.list['TreeFormats'] = [self.info['SaveType']] + rje.split(rje.join(self.list['TreeFormats']).lower())
        for format in self.list['TreeFormats'][0:]:
            if format not in treeformats or self.list['TreeFormats'].count(format) > 1: self.list['TreeFormats'].remove(format)
        if self.info['SaveTree'] in ['','none'] and self.list['TreeFormats'] and self.info['Basefile'].lower() not in ['','none']:
            self.info['SaveTree'] = '%s.%s' % (self.info['Basefile'],self.list['TreeFormats'][0])
        if self.info['CWTree'] != 'None':
            seqfile = self.info['CWTree']
            self.info['MakeTree'] = 'clustalw'
        ### ~ [3] ~ Actions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.opt['AutoLoad']: return
        if seqfile and seqfile.lower() != 'none':
            self.obj['SeqList'] = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=%s' % seqfile])
        if (self.obj['SeqList'] or self.info['DisIn'].lower() != 'none') and self.info['MakeTree'] != 'None':
            self.makeTree(make_seq=self.obj['SeqList'])
            #ltype = self.makeTree()
            #self.loadTree(file=self.info['Name'],seqlist=treeseq,type=ltype)
            if self.opt['OutputBranchLen']: withbranchlengths = 'Length'
            else: withbranchlengths = 'none'
            outnames = self.info['OutNames']
            maxnamelen = self.stat['TruncNames'] 
            self.saveTrees(seqname=outnames,blen=withbranchlengths)
            #if self.obj['SeqList']: self.saveTree(seqnum=self.opt['SeqNum'],type=self.info['SaveType'],seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            #else: self.saveTree(seqnum=False,type=self.info['SaveType'],seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            #x#if self.opt['SoapLab']: self.rTree('output.png',seqname=outnames,blen=withbranchlengths,compress=False)
        elif self.info['Name'] != 'None':
            if 'phbin=%s' % self.info['Name'] in self.cmd_list:
                self.loadTree(file=self.info['Name'],seqlist=self.obj['SeqList'],type='phb')
            else: self.loadTree(file=self.info['Name'],seqlist=self.obj['SeqList'])
            if self.info['SaveTree'].lower() not in ['','none'] or self.list['TreeFormats']: 
                if self.opt['OutputBranchLen']: withbranchlengths = 'Length'
                else: withbranchlengths = 'none'
                outnames = self.info['OutNames']
                maxnamelen = self.stat['TruncNames']
                if self.getStrLC('SaveTree'): self.info['Basefile'] = rje.baseFile(self.info['SaveTree'])
                else: self.info['Basefile'] = rje.baseFile(self.info['Name'])
                self.saveTrees(seqname=outnames,blen=withbranchlengths)
                #if self.obj['SeqList']:
                #    self.saveTree(seqnum=self.opt['SeqNum'],type=self.info['SaveType'],seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
                #else: self.saveTree(seqnum=False,type=self.info['SaveType'],seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
        if ancseq and ancseq.lower() != 'none': self.mapAncSeq(ancseq)
        if self.info['Grouping'] != 'None' and self.obj['SeqList']: self._autoGroups(self.info['Grouping'])
#########################################################################################################################
    def _cutCmdList(self,cmds=None):     ### Sets reduced Attribute set from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [0] ~ Setup Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if cmds == None: cmds = self.cmd_list
        ### ~ [1] ~ Read commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in cmds:
            try:## ~ [1a] ~ General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self._cmdRead(cmd,type='info',att='Rooting',arg='root')
                self._cmdReadList(cmd,'stat',['DefLen','TextScale','RootBuffer','BootCut'])
                self._cmdReadList(cmd,'int',['MinFamSize','MinFamNum','Bootstraps','TruncNames','SpecDup'])
                self._cmdRead(cmd,type='int',att='MinFamSize',arg='mfs')
                self._cmdRead(cmd,type='int',att='MinFamNum',arg='fam')
                self._cmdRead(cmd,type='opt',att='Orphans',arg='orphan')
                self._cmdReadList(cmd,'opt',['Orphans','AllowVar','QryVar'])
                self._cmdRead(cmd,type='info',att='Grouping',arg='group')
                self._cmdRead(cmd,type='info',att='GroupSpecies',arg='groupspec')
                self._cmdRead(cmd,type='int',att='Bootstraps',arg='bootstrap')
                self._cmdReadList(cmd,'info',['OutNames','SaveType'])
                self._cmdRead(cmd,type='opt',att='OutputBranchLen',arg='branchlen')
                self._cmdReadList(cmd,'list',['TreeFormats'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Tree Reading/Building                                                                                   #
    ### => Read from various formats and convert to NSF string with branch lengths (these may be arbitrary)             #
    ### => NSF string is then processed to tree                                                                         #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### [0] Basic run options handled in self._cmdList(). ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.i() >= 0: treeMenu(self,self.log,self.cmd_list,self)
        except: self.errorLog('Error during Tree.run()')
#########################################################################################################################
    def loadTree(self,file='tree',type='nsf',boot=0,seqlist=None,postprocess=True):   ### Loads trees from saved files
        '''
        Calls appropriate method to load a tree from a certain format of file.
        >> file:str ['tree'] = filename
        >> type:str ['nsf'] = Type of file (e.g. nsf)
        - 'nsf' = Newick Standard Format
        - 'phb' = Newick Standard Format from ClustalW
        - 'sim' = build from *.sim.anc.fas file of GASP test simulation data
        >> boot:int [0] = number of bootstraps, if important (0 = ignore bootstraps)
        >> seqlist:SeqList object to map onto tree
        >> postprocess:boolean = Whether to re-root tree and identify duplications [True]
        '''
        self.stat['Bootstraps'] = boot
        self.info['Name'] = file
        self.node = []
        self.branch = []
        if type == 'sim':
            self.verbose(0,0,'treeFromSim() currently not implemented. Sorry.',1)
            raise ValueError    #X# self.treeFromSim(file)
        nsftree = self.treeFromNSF(file)
        self.log.printLog('#TREE','Loaded %s tree data from %s.' % (type,file))
        self.verbose(1,3,'\nTree: %s' % nsftree,2)
        self.buildTree(nsftree,seqlist,type=type,postprocess=postprocess)
#########################################################################################################################
    def treeFromNSF(self,nsf_file):     ### Reads tree from NSF file and returns nsftree
        '''
        Reads from Newick Standard Format file and returns NSF tree.
        * NOTE: whitespace is removed - do not have whitespace in sequence names
        >> nsf_file:str = filename
        << nsftree:str
        '''
        ### ~ [1] ~ Open file & Read Lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        file_lines = self.loadFromFile(filename=nsf_file,v=0,checkpath=False,chomplines=True)
        if not file_lines: raise ValueError
        ### ~ [2] ~ Remove whitespace and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        multi_seq_per_line = False
        for line in file_lines:
            if line.count(':') > 1: multi_seq_per_line = True; break
        if multi_seq_per_line: return rje.join(rje.split(rje.join(file_lines)),'')
        nsftree = ''
        for line in file_lines:
            if line.find(':') > 0:
                start = rje.split(rje.replace(line,':',' '))[0]
                end = rje.split(line,':')[-1]
                nsftree += '%s:%s' % (start,end)
                #self.deBug(nsftree)
            else: nsftree += line
        nsftree = nsftree[:nsftree.find(';')+1]
        #self.deBug(nsftree)
        return nsftree        
#########################################################################################################################
    def buildTree(self,nsftree,seqlist=None,type='nsf',postprocess=True):        ### Builds Tree Object from nsf string
        '''
        Builds tree from NSF string.
        >> nsftree:str = Newick Standard Format tree
        >> seqlist:SeqList object to map onto tree
        >> type:str = Type (Format) of tree
        >> postprocess:boolean = Whether to re-root tree and identify duplications [True]
        '''
        try:
        ### <a> ### Determine attributes of tree
            _stage = '<a> Determination of tree attributes'
            ## <i> ##  Sequence names
            seqnames = False    # Whether sequence names (as opposed to numbers) are used in tree
            re_txtname = re.compile('[A-Za-z]')
            if re_txtname.search(nsftree):
                seqnames = True
            ## <ii> ## Branchlengths
            nsftree = self.reformatIntegerBranchLengths(nsftree)
            re_branchlength = re.compile(':(\-?\d+\.\d+)[,\)\[]') 
            if re_branchlength.search(nsftree):
                self.opt['Branchlengths'] = True
            ## <iii> ## Bootstraps
            myboot = self.stat['Bootstraps']  # Empirical number of bootstraps
            if (type != 'phb' and re.search('\)\d+',nsftree)) or (type == 'phb' and re.search('\[\d+\][,\)]',nsftree)):
                self.opt['Bootstrapped'] = True
        except:
            self.verbose(0,3,nsftree,1)
            self.log.errorLog('Problem with BuildTree(%s).' % _stage)
            raise

        ### <b> ### Read in terminal branches and nodes to leave numbers only
        _stage = '<b> Read termini'
        self.node = []
        self.branch = []
        try:
            ## <i> ## Make RE pattern
            if self.opt['Branchlengths']:
                re_tipre = '([\(,])([=<>\w\.\'\#\-\/\[\]\{\}]+):(\-?\d+\.\d+)([,\)])'
                re_tipre = '([\(,])([^\(\),:]+):(\-?\d+\.\d+)([,\)])'
            else:
                re_tipre = '([\(,])([=<>\w\.\'\#\-\/\[\]\{\}]+)([,\)])'
                re_tipre = '([\(,])([^\(\),:]+)([,\)])'
            re_tip = re.compile(re_tipre)
            #self.deBug(nsftree)
            #self.deBug(re_tipre)
            #self.deBug(rje.matchExp(re_tipre,nsftree))
            ## <ii> ## Replacements
            tmptree = '%s' % nsftree
            self.stat['SeqNum'] = 0
            while re_tip.search(tmptree):
                m_tip = re_tip.search(tmptree)
                n_tip = m_tip.groups()
                name = n_tip[1]
                blen = self.stat['DefLen']
                if self.opt['Branchlengths']:
                    blen = float(n_tip[2])
                deltxt = m_tip.group()[:-1] + m_tip.group()[-1]
                tmptree = rje.strReplace(tmptree,deltxt,n_tip[0] + n_tip[-1])[0]
                nsftree = rje.strReplace(nsftree,deltxt,n_tip[0] + '%d' % self.stat['SeqNum'] + n_tip[-1])[0]
                #self.deBug(nsftree)
                ## <iii> ## New node and branch                
                newnode = Node(log=self.log,cmd_list=self.cmd_list)
                newnode.info['Name'] = name
                newnode.stat['ID'] = self.stat['SeqNum']
                newbranch = Branch(log=self.log,cmd_list=self.cmd_list)
                newbranch.stat['Length'] = blen
                newbranch.node = [newnode]
                self.node.append(newnode)
                newnode.branch = [newbranch]
                self.branch.append(newbranch)
                self.stat['SeqNum'] += 1
            ## <iv> ## Check treeseq matches seqnum
            self.log.printLog('#SEQ','%d sequences in Tree.' % self.stat['SeqNum'])
            #self.debugTree()
            self.verbose(2,3,nsftree,2)
        except:
            self.verbose(0,3,nsftree,1)
            self.verbose(0,3,tmptree,1)
            self.log.errorLog('Problem with reading terminal nodes.')
            raise
        ### Now have a tree with termini as numbers only

        ### <c> ### Build internal nodes and branches
        try:
            ## <i> ## Make RE pattern
            if self.opt['Branchlengths']:
                if type == 'phb':
                    re_clade = re.compile('\((\d+),(\d+)\):(\-?\d+\.\d+)\[*(\d*)\]*[,\)]')
                else:
                    re_clade = re.compile('\((\d+),(\d+)\)([\.\d]*):(\-?\d+\.\d+)[,\)]')
            else:
                if type == 'phb':
                    re_clade = re.compile('\((\d+),(\d+)\)\[*(\d*)\]*[,\)]')
                else:
                    re_clade = re.compile('\((\d+),(\d+)\)([\.\d]*)[,\)]')
            ## <ii> ## Search and replace
            nodex = self.stat['SeqNum']
            while re_clade.search(nsftree):
                m_clade = re_clade.search(nsftree)
                n_clade = m_clade.groups()
                if type == 'phb': n_clade = [n_clade[0],n_clade[1],n_clade[-1],n_clade[2]]
                n0 = int(n_clade[0])
                n1 = int(n_clade[1])
                if self.opt['Bootstrapped']:
                    try:
                        if '.' in n_clade[2]: boot = int(self.stat['Bootstraps']*rje.atof(n_clade[2])+0.5)
                        else: boot = rje.atoi(n_clade[2])
                    except: boot = 0
                    myboot = self._maxBoot(boot,myboot)
                else: boot = -1
                blen = self.stat['DefLen']
                if self.opt['Branchlengths']:
                    blen = float(n_clade[3])
                deltxt = '\\' + m_clade.group()[:-1]
                deltxt = re.sub('\)','\\)',deltxt)
                deltxt = re.sub('\]','\\]',deltxt)
                deltxt = re.sub('\[','\\[',deltxt)
                nsftree = re.sub(deltxt, '%d' % nodex, nsftree)
                ## <iii> ## New nodes and branches
                newnode = Node(log=self.log,cmd_list=self.cmd_list)
                newnode.info['Name'] = 'Node %d (%d,%d)' % (nodex+1,n0+1,n1+1) 
                newnode.stat['ID'] = nodex
                newbranch = Branch(log=self.log,cmd_list=self.cmd_list)
                newbranch.stat['Length'] = blen
                if self.opt['Bootstrapped']:
                    newbranch.stat['Bootstrap'] = boot
                newbranch.node = [newnode]
                self.node.append(newnode)
                self.branch.append(newbranch)
                newnode.branch = [self.branch[n0],self.branch[n1],newbranch]
                ## <iv> ## Existing nodes and branches                
                self.branch[n0].node.append(newnode)
                self.branch[n1].node.append(newnode)
                if int(self.branch[n0].node[0].stat['ID']) != n0 or int(self.branch[n1].node[0].stat['ID']) != n1:
                    self.log.errorLog('List index problems during Tree Build!',printerror=False)
                    raise ValueError
                nodex += 1
                #self.deBug(nsftree)
            self.verbose(1,3,"%d nodes in Tree. " % (nodex+1),0)
            #x#self.debugTree()
            self.stat['Bootstraps'] = myboot
            for s in range(self.stat['SeqNum']): self.branch[s].stat['Bootstrap'] = myboot            
        except:
            self.verbose(0,3,nsftree,1)
            self.log.errorLog('Problem with reading internal nodes.')
            raise
        
        ### <d> ### Add root or trichotomy
        try:
            ## <i> ## Rooting Status
            re_last = re.compile('\((\d+),(\d+),(\d+)\)')
            re_root = re.compile('\((\d+),(\d+)\)')
            if re_last.search(nsftree):     # Trichotomy
                self.opt['Rooted'] = False
            elif re_root.search(nsftree):   # Rooted
                self.info['RootMethod'] = 'Loaded'
                self.opt['Rooted'] = True
                re_last = re_root
            else:
                self.log.errorLog("Problem with tree formatting. %s remains." % nsftree,printerror=False)
                raise ValueError
            #self.deBug(nsftree)
            ## <ii> ## Add node
            m_last = re_last.search(nsftree)
            n_last = m_last.groups()
            n0 = int(n_last[0])
            n1 = int(n_last[1])
            n2 = int(n_last[-1])
            newnode = Node(log=self.log,cmd_list=self.cmd_list)
            newnode.stat['ID'] = nodex
            ## <iii> ## Existing nodes and branches                
            self.branch[n0].node.append(newnode)
            self.branch[n1].node.append(newnode)
            if self.opt['Rooted']:
                newnode.branch = [self.branch[n0],self.branch[n1]]
                self.verbose(0,3,'(Rooted)',1)
                newnode.info['Name'] = 'Root Node %d (%d,%d)' % (nodex+1,n0+1,n1+1)
            else:
                newnode.branch = [self.branch[n0],self.branch[n1],self.branch[n2]]
                self.branch[n2].node.append(newnode)
                self.verbose(0,3,'(Unrooted)',1)
                newnode.info['Name'] = 'Node %d (%d,%d,%d)' % (nodex+1,n0+1,n1+1,n2+1)
            self.node.append(newnode)
            if int(self.branch[n0].node[0].stat['ID']) != n0 or int(self.branch[n1].node[0].stat['ID']) != n1 or int(self.branch[n2].node[0].stat['ID']) != n2:
                self.log.errorLog('List index problems during Tree Build!', False, False)
                raise ValueError
            nodex+=1
            #X#self.debugTree()
            if len(self.node) != nodex or len(self.branch) != (nodex-1):
                self.log.errorLog('List size problems during Tree Build!',printerror=False)
                raise ValueError
        except:
            self.verbose(0,3,nsftree,1)
            self.log.errorLog('Problem with final node.')
            raise

        ### <e> ### Raise IDs to be 1-N not 0-N and Map Sequences if Given
        try:
            _stage = '<e> Raising IDs'
            ## Check Node numbers in names ##
            nodenameid = True
            for node in self.node:
                if seqnames and rje.matchExp('^(\d+)_(.+)$',node.info['Name']):
                    mname = rje.matchExp('^(\d+)_(.+)$',node.info['Name'])
                    if rje.atoi(mname[0]) > self.nodeNum(): nodenameid = False
                else: nodenameid = False
                if not nodenameid: break
            for node in self.node:
                node.stat['ID'] += 1
                #!# Check that this next change does not cause problems!
                if seqnames and nodenameid:
                    mname = rje.matchExp('^(\d+)_(.+)$',node.info['Name'])
                    node.info['Name'] = mname[1]
                    node.stat['ID'] = rje.atoi(mname[0])
                    #self.deBug('%s:%s' % (node.info, node.stat))
            #!# Check that this next change does not cause problems!
            if seqnames : self._reorderNodes()
            self._renumberNodes()

        ### <f> ### Rescale default drawing
            _stage = '<f> Rescaling'
            max_branch_len = 0.0
            for branch in self.branch:
                if branch.stat['Length'] > max_branch_len: max_branch_len = branch.stat['Length'] 
            if max_branch_len > 2 or max_branch_len < 0.1: self.stat['DefLen'] = max_branch_len / 10.0
            for cmd in self.cmd_list: self._cmdRead(cmd,type='stat',att='DefLen')

        ### <g> ### Option to re-root
            _stage = '<g> Mapping sequences'
            if seqlist or self.obj['SeqList']: self.mapSeq(seqlist)
            if self.stat['Verbose'] >= 0: self.textTree()
            if postprocess: self.treeRoot()

        except:
            self.verbose(0,3,nsftree,1)
            self.log.errorLog('Problem with BuildTree(%s).' % _stage)
            raise
#########################################################################################################################
    def reformatIntegerBranchLengths(self,nsftree):     ### Reformats any integer branch lengths to floats
        '''
        Reformats any integer branch lengths to floats.
        >> nsftree:str = NSF tree
        << newtree:str = Returned tree with reformatted branchlengths
        '''
        try:
            newtree = nsftree
            mtree = rje.matchExp('^(.+:)(\d+)([,\)\[].+)$',newtree)
            while mtree:
                newtree = '%s%f%s' % (mtree[0],float(mtree[1]),mtree[2])
                mtree = rje.matchExp('^(.+:)(\d+)([,\)\[].+)$',newtree)
            return newtree
        except:
            self.log.errorLog('Problem with reformatIntegerBranchLengths(). Keeping input tree format.')
            return nsftree
#########################################################################################################################
    def debugTree(self): # ! # Tmp
        b = 1
        print()
        for branch in self.branch:
            print(b, branch.node)
            b += 1
        raw_input()
#########################################################################################################################
    def _checkTree(self):       ### Checks integrity of tree - that nodes and branches match
        '''
        Checks integrity of tree
        - that nodes and branches match
        - nodes have correct number of branches (1-3)
        - branches have correct number of nodes (2)
        << True if OK. False if Bad.
        '''
        try:
            dandy = True
            ### <a> ### Check Numbers
            if len(self.node) != (len(self.branch)+1):
                self.verbose(0,2,'Wrong number of nodes (%d) to branches (%d)!' % (len(self.node),len(self.branch)),1)
                dandy = False 
            ### <b> ### Check Nodes
            for node in self.node:
                if node.stat['ID'] <= self.stat['SeqNum'] and len(node.branch) != 1:
                    self.verbose(0,2,'Wrong number of branches (%d) for terminal node %d (%s)!' % (len(node.branch),node.stat['ID'],node.info['Name']),1)
                    dandy = False
                if len(node.branch) > 3 or len(node.branch) < 1:
                    self.verbose(0,2,'Wrong number of branches (%d) for node %d!' % (len(node.branch),node.stat['ID']),1)
                    dandy = False
                for branch in node.branch:
                    if node in branch.node:
                        continue    # OK
                    else:
                        self.verbose(0,2,'Node %d missing from one of its branches (%s)!' % (node.stat['ID'],branch.show()),1)
                        dandy = False
            ### <c> ### Check Branches
            for branch in self.branch:
                if len(branch.node) != 2:
                    self.verbose(0,2,'Wrong number of nodes (%d) for branch %s!' % (len(branch.node),branch.show()),1)
                    dandy = False
                for node in branch.node:
                    if branch in node.branch:
                        continue    # OK
                    else:
                        self.verbose(0,2,'Branch %s missing from one of its nodes (%d)!' % (branch.show(),node.stat['ID']),1)
                        dandy = False
            return dandy
        except:
            self.log.errorLog('Major problem during _checkTree()')
            return False
#########################################################################################################################
    def _renumberNodes(self):       ### Renumbers nodes from root
        '''Renumbers nodes from root.'''
        try:
            ### Setup ###
            ptext = 'Renumbering Internal Nodes (Root:%s)' % self.opt['Rooted']
            self.log.printLog('\n#TREE',ptext,log=False,newline=False)
            intorder = []
            inttarget = (self.nodeNum() - self.stat['SeqNum'])
            next = []
            if self.opt['Rooted']:
                next = [self._getRootNode()]
            else:
                maxid = 0
                for node in self.node:
                    if node.stat['ID'] > maxid:
                        maxid = node.stat['ID']
                        next = [node]

            ### Renumber internal nodes ###
            while len(next) > 0:
                self.log.printLog('\r#TREE','%s ... %.1f%%' % (ptext,(100.0 * len(intorder))/inttarget),log=False,newline=False)
                node = next.pop(0)
                intorder.append(node)
                for branch in node.branch:
                    nextnode = branch.link(node)
                    if nextnode in intorder or nextnode.stat['ID'] <= self.stat['SeqNum']: continue
                    else:
                        #self.deBug('%s %s' % (nextnode.info,nextnode.stat))
                        next.append(nextnode)
            self.log.printLog('\r#TREE','%s ... %.1f%%' % (ptext,(100.0 * len(intorder))/inttarget),log=False)
            if len(intorder) != inttarget:
                a = len(intorder)
                b = (self.nodeNum() - self.stat['SeqNum'])
                self.log.errorLog('Wrong number of internal nodes! (%d not %d-%d = %d)' % (a,self.nodeNum(),self.stat['SeqNum'],b),printerror=False)
                print(self.nodeList(intorder))
                print(self.nodeList(self.node))
                raise ValueError
            else:
                n = self.nodeNum()
                for node in intorder:
                    node.stat['ID'] = n
                    n -= 1
            if n != self.stat['SeqNum']:
                self.log.errorLog('Next (terminal) node should be numbered %d but %d sequences!' % (n,self.stat['SeqNum']),printerror=False)
                raise ValueError
            self._reorderNodes()

            ### Rename ###
            for node in self.node:
                node.rename(rooting=self.info['RootMethod'])
            for branch in self.branch:
                if branch.node[0].stat['ID'] > branch.node[1].stat['ID']:
                    branch.info['Name'] = '%s -> %s' % (branch.node[0].shortName(),branch.node[1].shortName())
                else:
                    branch.info['Name'] = '%s -> %s' % (branch.node[1].shortName(),branch.node[0].shortName())
        except:
            self.log.errorLog('Problem in renumbering internal nodes: check there is <=1 trichotomy.')
            raise
#########################################################################################################################
    def _reorderNodes(self):       ### Reorders nodes in ID order (for output clarity)
        '''Renumbers nodes from root.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            n = 1
            ptext = 'Reordering Internal Nodes (Root:%s) ...' % self.opt['Rooted']
            orderednodes = []
            self.log.printLog('\r#TREE','%s %.1f%%' % (ptext,(100.0 * n)/self.nodeNum()),log=False,newline=False)
            ### ~ [2] Reorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while n <= self.nodeNum():
                pren = n
                for node in self.node:
                    #X#self.log.printLog('\r#TREE','%s ... %d' % (ptext,n),log=False,newline=False)
                    if node.stat['ID'] == n:
                        orderednodes.append(node)
                        n += 1
                        self.log.printLog('\r#TREE','%s %.1f%%' % (ptext,(100.0 * n)/self.nodeNum()),log=False,newline=False)
                if n == pren:
                    self.log.printLog('\r#TREE','%s %.1f%%' % (ptext,(100.0 * n)/self.nodeNum()),log=False)
                    self.log.errorLog('Problem reordering nodes! Stuck on %d' % n,printerror=False)
                    raise ValueError
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('\r#TREE','%s %.1f%%' % (ptext,100.0),log=False)
            self.node = orderednodes
        except:
            self.log.errorLog('Problem reordering nodes (%d)' % n)
            raise
#########################################################################################################################
    def _maxBoot(self,b,maxb):   ### Returns predicted no. Bootstraps
        '''
        Returns predicted number of bootstraps assuming a power of ten.
        >> b:int = branch bootstrap
        >> max:int = current max
        '''
        boot = max(10,maxb)
        while b > boot: boot *= 10
        return boot
#########################################################################################################################
    def _nodesNumbered(self,check1toN=False):   ### Whether sequence names are purely numbers
        '''
        Whether sequence names are purely numbers.
        >> check1toN:boolean = whether to also check that the numbers are 1 to N [False]
        << Returns True or False
        '''
        try:
            numlist = []
            for node in self.node:
                if node.stat['ID'] < self.stat['SeqNum'] and re.search('\D',node.info['Name']): return False
                elif node.stat['ID'] <= self.stat['SeqNum'] and check1toN: numlist.append(rje.atoi(node.info['Name']))
            numlist.sort()
            if check1toN and numlist != range(1,len(numlist)+1): return False
            return True
        except:
            self.log.errorLog('Error during _nodesNumbered()')
            raise
#########################################################################################################################
    def mapAncSeq(self,filename=None):  ### Maps Ancestral Sequences to tree
        '''
        Maps Ancestral Sequences to tree.
        >> filename:str = name of ancestral sequence file. Should be in GASP output format.
        '''
        try:
        ### <a> ### Load sequences in SeqList object
            _stage ='<a> Load Sequences'
            if self.stat['SeqNum'] < 1:
                self.log.errorLog('Cannot Map Ancestral Sequences from %s - no tree to map to!' % filename,printerror=False)
                raise ValueError
            anclist = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=%s' % filename,'accnr=F'])

        ### <b> ### Map Terminal sequences
            _stage ='<b> Map Termini'
            ## <i> ##  Setup
            _stage ='<b-i> Map Termini - Setup'
            #X#self.verbose(0,3,'Matching Ancestral Sequence Names to Tree Names...',1)
            seqdic = anclist.seqNameDic(key='NumName')
            #self.deBug(seqdic)
            ## <ii> ## Work through each leaf in turn
            _stage ='<b-ii> Map Termini - Leaves'
            mapseq = []
            mapn = 0
            for node in self.node[:self.stat['SeqNum']]:
                self.log.printLog('\r#ANC','Mapping Ancestral Sequences from %s: %.1f%%' % (filename,100.0 * mapn / self.nodeNum()),newline=False,log=False)
                mapped = None
                node.obj['Sequence'] = None  # Clear sequence
                ## <iii> ## Text Match
                _stage ='<b-iii> Map Termini - Text Match'
                if node.info['Name'] in seqdic:   # Exact match
                    mapped = seqdic[node.info['Name']]
                else:   # Partial Match
                    for seq in anclist.seq:
                        if seq.info['Name'].find(node.info['Name']) > 1:   # Matches somewhere in Name (but not at start - that's the node number!)
                            mapped = seq
                            break
                        elif rje.matchExp('^\d+\s+(\S+)',seq.info['Name']) and node.info['Name'].find(rje.matchExp('^\d+\s+(\S+)',seq.info['Name'])[0]) == 0:  ## Shortname at start of nodename
                            mapped = seq
                            break
                            
                ## <iv> ## Make match
                _stage ='<b-iv> Map Termini - Make match'
                #print mapped
                if mapped:
                    #self.deBug('%s = %s' % (seq.info['Name'],node.info['Name']))
                    if mapped in mapseq:
                        self.log.errorLog('WARNING!: %s maps to multiple nodes!' % (mapped.info['Name']),printerror=False)
                        raise ValueError
                    ## <v> ## Remake sequence name
                    _stage ='<b-v> Map Termini - SeqName'
                    if re.search('^(\d+)\s(\S.+)',mapped.info['Name']):
                        details = rje.matchExp('^(\d+)\s(\S.+)',mapped.info['Name'])
                        #print details
                        s = rje.atoi(details[0])
                        mapped.info['Name'] = details[1]
                        mapped.extractDetails()
                    else:
                        self.log.printLog('AncSeqMatch Error - Node %d (%s) mapped to Sequence %s without ID.' % (node.stat['ID'],node.info['Name'],mapped.info['ID']))
                        raise ValueError
                    self.verbose(1,4,'\rNode %d (%s) mapped to Sequence %d (%s).' % (node.stat['ID'],node.info['Name'],s,mapped.info['ID']),1)
                    #print mapped,anclist.seq.index(mapped)
                    ## <vi> ## Map Seq
                    _stage ='<b-vi> Map Termini - MapSeq'
                    node.mapSeq(seq=mapped,id=s)
                    mapseq.append(mapped)
                    mapn += 1
                    #print mapn
                else:
                    self.log.errorLog('Unable to map Node %d (%s). Check tree names within sequence names (or numbered).' % (node.stat['ID'],node.info['Name']),printerror=False)
                    raise ValueError
            ## <vii> ## Finish off
            _stage ='<b-vii> Map Termini - Finish'
            #X#self.verbose(0,2,'%d of %d leaf nodes mapped to %d sequences.' % (mapn,self.stat['SeqNum'],len(mapseq)),1)
            self.log.printLog('\r#ANC','Mapping Ancestral Sequences from %s: %d Termini mapped.' % (filename,mapn),log=False)
            if mapn != len(mapseq):
                self.log.errorLog('Mapping Error: Number of mapped nodes (%d) != mapped sequences (%d)!' % (mapn,len(mapseq)),printerror=False,quitchoice=True)

        ### <c> ### Map ancestral sequences based on IDs            
            _stage ='<c> Map Ancestrals'
            ## <i> ## Work through each internal node in turn
            for node in self.node[self.stat['SeqNum']:]:
                self.log.printLog('\r#ANC','Mapping Ancestral Sequences from %s: %.1f%%' % (filename,100.0 * mapn / self.nodeNum()),newline=False,log=False)
                mapped = None
                node.obj['Sequence'] = None  # Clear sequence
                ## <iii> ## Match according to links
                for seq in anclist.seq:
                    if re.search('^\d.+\S*Node \d+ \((\S+)\)',seq.info['Name']):
                        links = rje.split(rje.matchExp('^\d.+\S*Node \d+ \((\S+)\)',seq.info['Name'])[0],',') # List of link IDs
                        mapped = seq
                        for desc in node.neighbours(ignore=[node.ancNode()]):
                            d = '%d' % desc.stat['ID']
                            #print d, links
                            if d not in links:
                                mapped = None
                        if mapped:
                            break
                    else:
                        pass
                        #self.deBug('%s not internal?' % seq.info['Name'])
                ## <iv> ## Make match
                #print mapped
                if mapped:
                    if mapped in mapseq:
                        self.log.errorLog('WARNING!: %s maps to multiple nodes!' % (mapped.info['Name']),printerror=False)
                        raise ValueError
                    ## <v> ## Remake sequence name
                    if re.search('^(\d+)\s(\S.+)',mapped.info['Name']):
                        details = rje.matchExp('^(\d+)\s(\S.+)',mapped.info['Name'])
                        #print details
                        s = rje.atoi(details[0])
                        mapped.info['Name'] = details[1]
                        mapped.extractDetails()
                    else:
                        self.log.printLog('AncSeqMatch Error - Node %d (%s) mapped to Sequence %s without ID.' % (node.stat['ID'],node.info['Name'],mapped.info['ID']),printerror=False)
                        raise ValueError
                    self.verbose(1,4,'\rNode %d (%s) mapped to Sequence %d: %s.' % (node.stat['ID'],node.info['Name'],s,self.nodeList(node.neighbours(ignore=[node.ancNode()]))),1)
                    #print mapped,anclist.seq.index(mapped)
                    ## <vi> ## Map Seq
                    node.mapSeq(seq=mapped,id=s)
                    mapseq.append(mapped)
                    mapn += 1
                    #print mapn
                else:
                    self.log.errorLog('Unable to map Node %d (%s). Check tree names within sequence names (or numbered).' % (node.stat['ID'],node.info['Name']),printerror=False)
                    raise ValueError
            ## <vi> ## Finish off
            self.log.printLog('\r#ANC','Mapping Ancestral Sequences from %s: %d nodes mapped.' % (filename,mapn))
            #X#self.verbose(0,2,'%d of %d nodes mapped to %d sequences.' % (mapn,len(self.node),len(mapseq)),1)
            if mapn != len(mapseq):
                self.log.errorLog('Mapping Error: Number of mapped nodes (%d) != mapped sequences (%d)!' % (mapn,len(mapseq)),printerror=False,quitchoice=True)

        ### <d> ### Reorder self.node (for output clarity only)
            _stage ='<d> Reorder nodes'
            self.obj['SeqList'] = copy.deepcopy(anclist)
            self.obj['SeqList'].seq = mapseq
            self.obj['SeqList'].querySeq()
            self._renumberNodes()
            if self.opt['Rooted']: self.findDuplications(duptext='New sequence info')
            #self.textTree()
           
        except:
            self.obj['SeqList'] = None
            for node in self.node:
                node.obj['Sequence'] = None
            self.log.errorLog('Major Problem in mapAncSeq(%s):' % _stage)
#########################################################################################################################
    def clearSeq(self):     ### Clears current sequence information
        '''Clears current sequence information.'''
        self.obj['SeqList'] = None
        for node in self.node: node.obj['Sequence'] = None
#########################################################################################################################
    def mapSeq(self,seqlist=None,clearseq=True):      ### Maps SeqList object onto tree
        '''
        Maps SeqList object onto Tree.
        >> seqlist:rje_seq.SeqList Object.
        >> clearseq:boolean = whether to clear any current node sequence before matches
            - if matching from multiple SeqLists then set to False for 2nd and subsequent matches
            *** It is best to combine into a single SeqList object, which can be linked to the tree object. ***

        Names should either include tree sequence names (at start of name) or be ordered according to tree sequences (1-N).
        Will take sequences from a larger seqlist than the tree and try to find each node, first with an exact match and then partial.
        '''
        try:
            ### Setup SeqList ###
            if seqlist:
                self.obj['SeqList'] = seqlist
            else:
                seqlist = self.obj['SeqList']
            if not seqlist:
                self.log.errorLog('Major error with Tree.mapSeq() - no SeqList given!',printerror=False)
                self.clearSeq()
                return False
            if not seqlist.seq:
                self.log.errorLog('Major error with Tree.mapSeq() - SeqList has no sequences!',printerror=False)
                self.clearSeq()
                return False

            ### Setup Naming Formats ###
            seqnumbers = seqlist.numbersForNames(check1toN=True)    # Sequence names are pure numbers 1 to N - use EXACT match
            treenumbers = self._nodesNumbered(check1toN=True)       # Leaves are pure numbers 1 to N - match numbers unless seqnumbers
            numbernames = True      # Sequence names have a node number followed by the sequence names: Match number or Name
            numlist = []
            for seq in seqlist.seq:
                if rje.matchExp('^(\d+)\s+\S',seq.info['Name']):
                    numlist.append(rje.atoi(rje.matchExp('^(\d+)\s+\S',seq.info['Name'])[0]))
                else:
                    numbernames = False
            numlist.sort()
            if numbernames and numlist != range(1,len(numlist)+1):
                numbernames = False
            matchnumbers = False    # Whether to match using terminal node numbers [True] or sequence name [False]
            #!#print treenumbers, numbernames, seqnumbers
            if treenumbers and not numbernames and not seqnumbers:
                matchnumbers = True
                self.verbose(0,3,'Matching Sequence Names to Terminal Node Numbers...',1)
            else:
                self.verbose(0,3,'Matching Sequence Names to Tree Names...',1)
                #X#seqdic = seqlist.seqNameDic(key='Name')

            ### Map each leaf in turn ###
            mapseq = []     # List of Mapped sequence objects
            mapn = 0        # Number of mapped leaf nodes
            for node in self.node[:self.stat['SeqNum']]:    # Each leaf node in turn
                mapped = {'EXACT':[],'NUMBERS':[],'UNDERSCORE':[]}  # Dictionary of Method:Mapped sequence(s)
                if clearseq:
                    node.obj['Sequence'] = None  # Clear sequence
                ## Numbered match ##
                if matchnumbers:    ### Get appropriate sequence
                    num = int(node.info['Name'])
                    if seqlist.seqNum() < num:  # Not enough sequences!
                        self.log.errorLog('Not enough sequences to map leaf numbered %d' % num,printerror=False)
                        raise ValueError
                    else:   # Map to number
                        mapped['NUMBERS'].append(seqlist.seq[num-1])
                ## Name Match ##
                else:
                    name = rje.split(node.info['Name'])[0]
                    if treenumbers:
                        num = rje.atoi(name)
                    for seq in seqlist.seq:
                        seqname = seq.shortName()
                        if numbernames and len(rje.split(seq.info['Name'])) > 1:
                            seqname = rje.split(seq.info['Name'])[1]
                        ## EXACT ##
                        if seqname == name:
                            mapped['EXACT'].append(seq)
                        ## NUMBERS ##
                        elif numbernames and treenumbers and num == rje.atoi(seq.shortName()):
                            mapped['NUMBERS'].append(seq)
                        ## UNDERSCORE ##
                        else:
                            if name.find(seq.shortName() + '_') == 0:
                                mapped['UNDERSCORE'].append(seq)
                            elif numbernames and name.find(seq.shortName() + '_' + seqname + '_') == 0:
                                mapped['UNDERSCORE'].append(seq)
                ## Extract best match ##
                nmapped = None
                for method in ['EXACT','NUMBERS','UNDERSCORE']:
                    #X#print method, mapped[method]
                    if mapped[method]:
                        if len(mapped[method]) > 1:
                            self.log.errorLog('Node "%s" has %d %s mappings!' % (name,len(mapped['EXACT']),method),printerror=False)
                            raise ValueError
                        else:
                            nmapped = mapped[method][0]
                            break
                ## Make match ##
                if nmapped:
                    s = seqlist.seq.index(nmapped) + 1
                    self.verbose(1,4,'Node %d (%s) mapped to Sequence %d (%s).' % (node.stat['ID'],node.info['Name'],s,nmapped.shortName()),1)
                    #print mapped,seqlist.seq.index(mapped)
                    node.mapSeq(seq=nmapped,id=s)
                    #print mapseq
                    if nmapped in mapseq:
                        self.log.errorLog('Sequence %d (%s) maps to multiple nodes!' % (s,nmapped.shortName()),printerror=False)
                        raise ValueError
                    else:
                        mapseq.append(nmapped)
                    mapn += 1
                    #print mapn
                else:
                    self.log.errorLog('Unable to map Node %d (%s). Check tree names match sequence names (or are numbered).' % (node.stat['ID'],node.info['Name']))
                    raise ValueError
            ## <vi> ## Finish off
            self.verbose(0,2,'%d of %d leaf nodes mapped to %d sequences.' % (mapn,self.stat['SeqNum'],len(mapseq)),1)
            if mapn != len(mapseq):
                self.log.errorLog('Mapping Error: Number of mapped nodes (%d) != mapped sequences (%d)!' % (mapn,len(mapseq)),printerror=False)
                raise ValueError
            self.obj['SeqList'] = seqlist
        ### <c> ### Reorder self.node (for output clarity only)
            self._renumberNodes()
            if self.opt['Rooted']: self.findDuplications(duptext='New sequence info')
            return True
        except:
            self.clearSeq()
            self.log.errorLog('Major problem mapping SeqList object onto tree.',quitchoice=True)
            return False
#########################################################################################################################
    def _prune(self,branch=None,remtext='Manual Tree Pruning'):       ### Removes branch and all descendants
        '''
        Removes branch and all descendants.
        >> branch:Branch object
        >> remtext:str = text description for reason of sequence removal
        '''
        try:
            ### <a> ### Lists of stuff to remove
            delnode = []        # Nodes to remove
            delbranch = []      # Branches to remove
            for node in self.node:
                if branch in self.rootPath(node):   # If deleted branch is ancestral
                    delnode.append(node)
                    delbranch.append(node.ancBranch())
            ### <b> ### Combine other branches
            fami = -1
            anc = branch.ancNode()
            if anc in self.subfam:
                fami = self.subfam.index(anc)
                self.subfam[fami] = None
            delnode.append(anc)
            if anc == self.node[-1] and self.opt['Rooted']:    # Root
                delbranch.remove(branch)
                delbranch += anc.branch
            else:   # Internal
                mergeb = copy.copy(anc.branch)
                mergeb.remove(branch)
                ancb = mergeb[0]
                descb = mergeb[1]
                if ancb.commonNode(descb) != anc:
                    print('Problem with ancestral node not linking correct branches!')
                    raise ValueError
                if anc == self.node[-1]:  # Trichotomy
                    ancn = ancb.link(anc)
                elif anc.ancBranch() == descb:
                    [ancb,descb] = [descb,ancb]
                    ancn = ancb.ancNode()
                else:
                    ancn = ancb.ancNode()
                self.verbose(0,3,"Remove node %d and combine branches %s and %s." % (anc.stat['ID'],ancb.show(),descb.show()),1)
                ancb.combine(anc,descb)
                ancb.link(ancn).branch.append(ancb)
                delbranch.append(descb)
                if fami >= 0:   # Move subfam to new node
                    self.subfam[fami] = ancb.descNode()
                    fami = -1
            ### <c> ### Remove nodes, seqs and branches
            ## <i> ## Tidy = reduce delnode and delbranch to single copies of each node/branch
            for node in delnode:
                while delnode.count(node) > 1:
                    delnode.remove(node)
            for branch in delbranch:
                while delbranch.count(branch) > 1:
                    delbranch.remove(branch)
            self.verbose(1,3,'Removing %d nodes and %d branches...' % (len(delnode),len(delbranch)),1)
            ## <ii> ## Remove nodes
            for node in delnode:
                if self.obj['SeqList'] != None and node.obj['Sequence'] != None:
                    self.obj['SeqList'].removeSeq(remtext,node.obj['Sequence'])    
                #if node.stat['ID'] <= self.seqNum():
                if node in self.subfam:
                    self.subfam.remove(node)
                for othernode in self.node:
                    if node.stat['ID'] < othernode.stat['ID']:
                        othernode.stat['ID'] -= 1
                self.node.remove(node)
                self.stat['SeqNum'] = self.seqNum()
            ## <iii> ## Remove branches
            for branch in delbranch:
                for node in self.node:
                    if branch in node.branch:
                        node.branch.remove(branch)
                self.branch.remove(branch)
            ### <d> ### Adjust
            if self._checkTree() == False:
                raise ValueError
            self._renumberNodes()
            if fami >= 0:
                self.subfam[fami] = self.node[-1]
        except:
            self.log.errorLog('Major Problem in _prune().')
            raise
#########################################################################################################################
    ### <3> ### Tree Making                                                                                             #
#########################################################################################################################
    def makeTreeMenu(self,interactiveformenu=0,force=False,make_seq=None): ### Menu for making tree
        '''
        Menu for making tree.
        >> interactiveformenu:int = Interactive level at which to give menu options [0]
        >> force:boolean = Whether to force making of tree (not allow exit without) [False]
        >> make_seq:SeqList object from which to make tree
        '''
        try:
            ### Setup ###
            if self.stat['Interactive'] < interactiveformenu:
                self.makeTree(make_seq)
                return
            ### Prepare Menu Parameters ##
            seqin = 'None'
            if make_seq:
                seqin = make_seq.info['Name']
            ### Show Settings ###
            print('\n\n### Tree Making Menu ###')
            print('Tree Generation Method: %s' % self.info['MakeTree'])
            print('Bootstraps: %d' % self.stat['Bootstraps'])
            print('Sequence File: %s' % seqin)
            if self.info['MakeTree'] in ['neighbor','upgma','fitch','kitch']:    ### PHYLIP
                print('Distance Matrix File: %s' % self.info['DisIn'])
            if self.info['MakeTree'] in ['neighbor','upgma','protpars','proml','fitch','kitch']:    ### PHYLIP
                print('Phylip Path: %s' % self.info['Phylip'])
                print('Phylip Options File: %s' % self.info['PhyOptions'])
                print('Protdist Options File: %s' % self.info['ProtDist'])
            if self.info['MakeTree'] in ['clustalw']:
                print('ClustalW Path: %s' % self.info['ClustalW'])
                print('Use Kimura multiple hit correction: %s' % self.opt['Kimura'])
            if self.info['MakeTree'] in ['fasttree']:
                print('FastTree Path: %s' % self.info['FastTree'])
            if self.info['MakeTree'] in ['iqtree']:
                print('IQTree Path: %s' % self.info['IQTree'])
            ### Options ###
            choicetext = '\n<M>ake Tree, Change <O>ptions'
            choices = ['M','O']
            if not force:
                choicetext += ', <Q>uit MakeTree'
                choices.append('Q')
            choice = ''
            while choice not in choices:
                choice = rje.choice('%s\n\nChoice?:' % choicetext,default='M').upper()[:1]
            if choice == 'M':
                self.makeTree(make_seq)
            elif choice == 'Q' and not force:
                return
            elif choice == 'O': ### Change Options ###
                self.info['MakeTree'] = rje.choice('Tree Generation Method:',self.info['MakeTree'],True)
                self.stat['Bootstraps'] = rje.getInt('Bootstraps:',default=self.stat['Bootstraps'],confirm=True)
                newseqin = rje.choice('Sequence File:',seqin,True)
                if newseqin != seqin:
                    make_seq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=%s' % newseqin,'autoload=T'])
                if self.info['MakeTree'] in ['neighbor','upgma','fitch','kitch']:    ### PHYLIP
                    self.info['DisIn'] = rje.choice('Distance Matrix File:',self.info['DisIn'],True)
                if self.info['MakeTree'] in ['neighbor','upgma','protpars','proml','fitch','kitch']:    ### PHYLIP
                    self.info['Phylip'] = rje.choice('Phylip Path:',self.info['Phylip'],True)
                    self.info['PhyOptions'] = rje.choice('Phylip Options File:',self.info['PhyOptions'],True)
                    self.info['ProtDist'] = rje.choice('Protdist Options File:',self.info['ProtDist'],True)
                if self.info['MakeTree'] in ['clustalw']:
                    self.info['ClustalW'] = rje.choice('ClustalW Path:', self.info['ClustalW'],True)
                    self.opt['Kimura'] = rje.yesNo('Use Kimura multiple hit correction?')
                if self.info['MakeTree'] in ['fasttree']:
                    self.info['FastTree'] = rje.choice('FastTree Path:', self.info['FastTree'],True)
                if self.info['MakeTree'] in ['iqtree']:
                    self.info['FastTree'] = rje.choice('IQTree Path:', self.info['IQTree'],True)
            self.makeTreeMenu(interactiveformenu,force,make_seq)
        except:
            self.log.errorLog('Major Problem in makeTreeMenu().',True,True)
#########################################################################################################################
    def makeTree(self,make_seq=None,keepfile=True):     ### Uses attributes to call program and make tree
        '''
        Uses attributes to call program and make tree.
        >> make_seq:SeqList object from which to make tree
        >> keepfile:bool = whether to keep the tree file produced by makeTree or delete
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not make_seq: make_seq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=T'])
            if make_seq.info['Name'] == 'None' and self.info['CWTree'].lower() != 'none':
                make_seq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=%s' % self.info['CWTree'],'autoload=T'])
            if make_seq.seqNum() < 2: return self.errorLog('Cannot make tree with %d sequences!' % make_seq.seqNum())
            ## ~ [0a] Setup tree parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.opt['ReRooted'] = False
            while self.info['MakeTree'].lower() in ['','none']:
                if self.stat['Interactive'] >= 0: self.info['MakeTree'] = rje.choice('Enter method or path for makeTree:',default='clustalw',confirm=True)
                else: self.info['MakeTree'] = 'clustalw'
            if not os.path.exists(self.info['MakeTree']): self.info['MakeTree'] = self.info['MakeTree'].lower()
            if self.info['DisIn'].lower() != 'none' and not os.path.exists(self.info['DisIn']):
                if make_seq:
                    self.errorLog('Distance Matrix file %s missing. Will make tree from %s.' % (self.info['DisIn'],make_seq.info['Name']),printerror=False)
                    self.info['DisIn'] = 'None'
                else:
                    self.errorLog('Distance Matrix file %s missing and no sequence file given!' % self.info['DisIn'],printerror=False)
                    raise IOError
            ### ~ [1] Generate tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] PHYLIP Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['MakeTree'] in ['neighbor','upgma','protpars','proml','fitch','kitch']:    ### PHYLIP
                ## Make Directory ##
                make_id = 'phylip_tmp_%s' % rje.randomString(8)
                os.mkdir(make_id)
                os.chdir(make_id)
                ## Run ##
                self.phylipTree(make_id,make_seq,self.info['DisIn'])
                ## Process results ##
                os.chdir('..')
                rje.fileTransfer(fromfile='%s%s%s' % (make_id,os.sep,self.log.info['LogFile']),tofile=self.log.info['LogFile'],deletefrom=True)
                if make_seq.info['Name'].lower() == 'none':
                    self.info['Name'] = '%s.%s.ph' % (rje.baseFile(self.info['DisIn'],True),self.info['MakeTree'])
                else:
                    self.info['Name'] = '%s.%s.ph' % (rje.baseFile(make_seq.info['Name'],True),self.info['MakeTree'])
                if os.path.exists(self.info['Name']):
                    self.verbose(2,1,'%s exists and will be overwritten.' % self.info['Name'],1)
                    os.unlink(self.info['Name'])
                self.log.printLog('#TREE','Unbootstrapped %s tree saved as %s.' % (self.info['MakeTree'].upper(),self.info['Name']))
                os.rename('%s%souttree' % (make_id,os.sep),self.info['Name'])
                ## Bootstraps ##
                if self.stat['Bootstraps'] > 0 and make_seq:    ### Make boostraps
                    tempstat = {'Interactive':self.stat['Interactive'],'Bootstraps':self.stat['Bootstraps']}
                    self.stat['Interactive'] = -1
                    self.opt['Bootstrapped'] = False
                    self.loadTree(file=self.info['Name'],seqlist=None,type='phb',postprocess=False)
                    self.setStat(tempstat)
                    os.chdir(make_id)
                    #self.verbose(0,2,'Generating %d bootstraps...' % self.stat['Bootstraps'],0)
                    for b in range(self.stat['Bootstraps']):
                        # Make Sequences #
                        bootseq = rje_seq.SeqList(log=self.log,cmd_list=['i=-1','v=%d' % (self.stat['Verbose']-1)])
                        bootseq.seq = []
                        boot_order = []
                        self.verbose(0,3,'Bootstrap %d of %d...' % ((b+1),self.stat['Bootstraps']),0)
                        for r in range(make_seq.seq[0].seqLen()):
                            boot_order.append(random.randint(0,make_seq.seq[0].seqLen()-1))
                        for seq in make_seq.seq:
                            randseq = ''
                            for i in boot_order:
                                randseq += seq.info['Sequence'][i]
                            bootseq._addSeq(seq.info['Name'],randseq)
                        self.verbose(0,3,'Sequences made. Making Phylip tree...',2)
                        # Run phylip #
                        self.phylipTree(make_id,bootseq,'none')
                        # Process Tree #
                        os.rename('outtree','boottree_%d' % b)
                        #self.verbose(0,2,rje.progressPrint(self,b+1,10,100),0)
                    #self.verbose(0,1,'Done!',1)
                    # Compile Bootstraps #
                    self.log.printLog('#TREE','Mapping %d bootstraps onto tree.' % self.stat['Bootstraps'])
                    self.verbose(0,2,'Compiling %d bootstraps...' % self.stat['Bootstraps'],0)
                    branch_clades = {}
                    for branch in self.branch:
                        branch.stat['Bootstrap'] = 0
                        real_clades = self.branchClades(branch)   # Node objects
                        branch_clades[branch] = ([],[])
                        for i in [0,1]:
                            for node in real_clades[i]:
                                branch_clades[branch][i].append(node.info['Name'])
                            branch_clades[branch][i].sort()
                    for b in range(self.stat['Bootstraps']):
                        boottree = Tree(log=self.log,cmd_list=['i=-1','v=-1'])
                        boottree.loadTree(file='boottree_%d' % b,seqlist=None,type='phb',postprocess=False)
                        for bbranch in boottree.branch:
                            bclades = boottree.branchClades(bbranch)   # Node objects
                            comp = ([],[])
                            for i in [0,1]:
                                for node in bclades[i]:
                                    comp[i].append(node.info['Name'])
                                comp[i].sort()
                            for branch in self.branch:
                                if comp[0][0:] in branch_clades[branch]:   # Boot +
                                    if comp[1][0:] in branch_clades[branch]: # OK
                                        branch.stat['Bootstrap'] += 1
                                        break
                                    else:   # Buggered
                                        self.log.errorLog('Knackers. Branch shares one clade but not the other!',printerror=False)
                        os.unlink('boottree_%d' % b)
                        rje.progressPrint(self,(b+1),1,10)
                    self.opt['Bootstrapped'] = True
                    self.verbose(0,1,'Done!',1)
                    # Tidy up #
                    os.chdir('..')
                    rje.fileTransfer(fromfile='%s%s%s' % (make_id,os.sep,self.log.info['LogFile']),tofile=self.log.info['LogFile'],deletefrom=True)
                    self.info['Name'] = '%s.nsf' % rje.baseFile(self.info['Name'])
                    self.opt['Bootstrapped'] = True
                    self.textTree()
                    self.saveTree(filename=self.info['Name'],type=self.info['SaveType'],seqname='num',maxnamelen=125,multiline=False)
                    self.mapSeq(make_seq)
                    self.textTree()
                    self.treeRoot()
                    #rtype = 'nsf'
                else:
                    self.loadTree(file=self.info['Name'],seqlist=make_seq,type='phb')
                ## Finish! ##                
                os.rmdir(make_id)
            ## ~ [1b] FastTree Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif self.info['MakeTree'] in ['fasttree','iqtree']:
                infile = '%s_tree.fas' % rje.baseFile(make_seq.info['Name'],True)
                make_seq.saveFasta(seqfile=infile,name='Short')
                self.info['Name'] = '%s.nsf' % rje.baseFile(make_seq.info['Name'],True)
                if self.info['MakeTree'] in ['fasttree']:
                    cpath = rje.makePath(self.info['FastTree'],wholepath=True)
                    seed = random.randint(1,999)
                    command = cpath + ' -boot %d -seed %d %s' % (self.stat['Bootstraps'],seed,infile)
                    if self.v() < 0: command = '%s -quiet' % command
                    self.printLog('#TREE',command)
                    try: self.buildTree(os.popen(command).read(),make_seq,type='nsf')
                    except:
                        self.errorLog('FastTree build failure. Will try ClustalW instead')
                        self.info['MakeTree'] = 'clustalw'
                        self.makeTree(make_seq,keepfile)
                        self.info['MakeTree'] = 'fasttree'
                else:
                    cpath = rje.makePath(self.info['IQTree'],wholepath=True)
                    seed = random.randint(1,999)
                    command = cpath + ' -s %s ' % infile
                    if self.stat['Bootstraps'] >= 100: command = '%s-b %d ' % (command,self.stat['Bootstraps'])
                    elif self.stat['Bootstraps']:
                        self.printLog('#BOOT','Boostraps increased from %s to 100 for IQTree.' % self.stat['Bootstraps'])
                        self.setStat({'Bootstraps':100})
                        command = '%s-b %d ' % (command,self.stat['Bootstraps'])
                    if not self.getBool('NoForks'):
                        if self.getInt('Forks') > 0:
                            command += ' -nt %d' % self.getInt('Forks')
                        else:
                            command += ' -nt AUTO'
                    if self.v() < 0: command = '%s -quiet' % command
                    self.printLog('#TREE',command)
                    if self.v() < 0: os.popen(command).read()
                    else: os.system(command)
                    try:
                        self.debug('%s.treefile' % infile)
                        self.debug(rje.exists('%s.treefile' % infile))
                        os.rename('%s.treefile' % infile,self.info['Name'])
                        self.buildTree(open(self.info['Name']).read(),make_seq,type='nsf')
                    except:
                        #!# Need to add file checks and/or cleanup and "redo" option.
                        self.errorLog('IQTree build failure. Will try ClustalW instead')
                        self.info['MakeTree'] = 'clustalw'
                        self.makeTree(make_seq,keepfile)
                        self.info['MakeTree'] = 'iqtree'
                if infile != make_seq.info['Name']: os.unlink(infile)
                if not keepfile and os.path.exists(self.info['Name']): os.unlink(self.info['Name'])
            ## ~ [1c] ClustalW Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                if self.info['MakeTree'] in ['clustalw','clustalw2']:
                    if make_seq.seqNum() < 4:
                        self.printLog('#TREE','%s has too few sequences (<4) for ClustalW Tree' % make_seq.info['Name'])
                        return self.upgmaIDTree(make_seq)
                    path = self.info['ClustalW']
                    infile = '%s_cw.fas' % rje.baseFile(make_seq.info['Name'],True)
                    make_seq.saveFasta(seqfile=infile,name='Number')
                else:
                    infile = make_seq.info['Name']
                    path = self.info['MakeTree']
                cpath = rje.makePath(path,wholepath=True)
                seed = random.randint(1,999)
                command = cpath + ' '
                if self.opt['Win32'] == False:
                    command += '-infile='
                command += '%s -bootstrap=%d -seed=%d' % (infile,self.stat['Bootstraps'],seed)
                if self.opt['Kimura']:
                    command += ' -kimura'
                self.log.printLog('#TREE',command)
                if self.v() < 0: os.popen(command).read()
                else: os.system(command)
                self.info['Name'] = '%s.phb' % rje.baseFile(make_seq.info['Name'],True)
                if os.path.exists(self.info['Name']):
                    os.unlink(self.info['Name'])
                os.rename('%s_cw.phb' % rje.baseFile(make_seq.info['Name'],True),self.info['Name'])
                #!# Warning, this replaces self Name with *.phb
                self.loadTree(file=self.info['Name'],seqlist=make_seq,type='phb')
                if infile != make_seq.info['Name']:
                    os.unlink(infile)
                if not keepfile:
                    os.unlink(self.info['Name'])
            #return rtype
        except: self.errorLog('Problem with makeTree().'); raise
#########################################################################################################################
    def upgmaIDTree(self,seqlist): ### Makes UPGMA tree based on MSA %ID
        '''Makes UPGMA tree based on MSA %ID.'''
        try:### ~ [1] Setup distance matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist.obj['MSA ID'] = None
            seqlist.addMatrix('MSA ID')
            seqlist._checkAln(aln=True,realign=True)
            for seq1 in seqlist.seqs():
                for seq2 in seqlist.seqs(): seqlist.getDis(seq1,seq2,'MSA ID')
            ### ~ [2] Make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.buildTree(seqlist.obj['MSA ID'].upgma(nosim=100.0),seqlist,type='nsf')
        except: self.errorLog('Problem during Tree.upgmaIDTree()'); raise
#########################################################################################################################
    def phylipTree(self,make_id,make_seq,make_dis):   ### Makes a Phylip tree from make_seq and make_dis
        '''
        >> make_id:random string for temp phylip directory
        >> make_seq:SeqList object from which to make tree
        >> make_dis:Distance Matrix file to build tree from
        '''
        ## Basic Options ##
        phy_opt = []
        if self.info['PhyOptions'] != 'None':
            if os.path.exists(self.info['PhyOptions']):
                PHYOPT = open(self.info['PhyOptions'],'r')
                phy_opt += PHYOPT.readlines()
                PHYOPT.close()
            else:
                self.log.errorLog('Phylip tree options file "%s" missing! Proceeding with defaults.' % self.info['PhyOptions'],True,False)
        #if make_seq.stat['Verbose'] < 0:
        #    phy_opt += ['2\n']
        phy_opt += ['3\ny\n']
        ## ProtDist Options ##
        pdis_opt = []
        if self.info['ProtDist'] != 'None' and make_dis.lower() == 'none':
            if os.path.exists(self.info['ProtDist']):
                PHYOPT = open(self.info['ProtDist'],'r')
                pdis_opt += PHYOPT.readlines()
                PHYOPT.close()
            else:
                self.log.errorLog('Phylip ProtDist options file "%s" missing! Proceeding with defaults.' % self.info['ProtDist'],True,False)
        #if make_seq.stat['Verbose'] < 0:
        #    pdis_opt += ['2\n']
        pdis_opt += ['y\n','\n','\n']
        ## Make Input files ##
        if make_dis.lower() == 'none' or self.info['MakeTree'] in ['protpars','proml']:
            make_seq.savePhylip(seqfile='infile')
        if self.info['MakeTree'] in ['protpars','proml','fitch','kitsch']:
            phy_opt = ['1\n'] + phy_opt
        phy_opt = ['j\n','%d\n' % (4 * random.randint(1,8000) + 1)] + phy_opt
        ## Make dismatrix ##
        if self.info['MakeTree'] in ['neighbor','upgma','fitch','kitsch']:
            if make_dis.lower() != 'none':
                rje.fileTransfer('..%s%s' % (os.sep,make_dis),'infile',False)
                self.log.printLog('#TREE','Using %s for %s tree distance matrix.' % (make_dis,self.info['MakeTree'].upper()))
            else:
                OPTFILE = open('%s.txt' % make_id,'w')
                OPTFILE.writelines(pdis_opt)
                OPTFILE.close()
                command = '%sprotdist < %s.txt' % (rje.makePath(self.info['Phylip']),make_id)
                self.log.printLog('#CMD',command)
                os.system(command)
                os.unlink('infile')
                os.rename('outfile','infile')
        ## Run tree program ##
        if self.info['MakeTree'] == 'upgma':
            phy_opt = ['n\n'] + pdis_opt
            command = '%sneighbor < %s.txt' % (rje.makePath(self.info['Phylip']),make_id)
        else:
            command = '%s%s < %s.txt' % (rje.makePath(self.info['Phylip']),self.info['MakeTree'],make_id)
        OPTFILE = open('%s.txt' % make_id,'w')
        OPTFILE.writelines(phy_opt)
        OPTFILE.close()
        self.log.printLog('#CMD',command)
        os.system(command)
        os.unlink('infile')
        os.unlink('outfile')
        os.unlink('%s.txt' % make_id)                        
#########################################################################################################################
    ### <4> ### Tree Output                                                                                             #
#########################################################################################################################
    def fullDetails(self,nodes=True,branches=True):  ### Displays details of Tree, Nodes and Branches
        '''Displays details of SeqList and all Sequences.'''
        self.verbose(0,1,self.details(),0)
        if nodes:
            for node in self.node: self.verbose(0,1,node.details(),0)
        if branches:
            for branch in self.branch: self.verbose(0,1,branch.details(),0)
#########################################################################################################################
    def sumTree(self):      ### Print summary of tree
        '''Prints Summary of Tree.'''
        try:
            if self.opt['Rooted']: sumtxt = '\nRooted Tree '
            else: sumtxt = '\nUnrooted Tree '
            if self.opt['Bootstrapped']: sumtxt += '(%d bootstraps).' % self.stat['Bootstraps']
            else: sumtxt += '(No bootstraps).'
            if self.opt['Branchlengths']: sumtxt += ' Branch Lengths given (total=%s).' % (rje.sf(self.treeLen(),3))
            else: sumtxt += ' No Branch Lengths.'
            sumtxt += ' %d nodes.' % self.nodeNum()
            self.verbose(0,1,sumtxt,2)
        except:
            self.log.errorLog('Problem with sumTree(). Execution Continued.')
#########################################################################################################################
    def treeLen(self): return self.pathLen(self.branch)
#########################################################################################################################
    def nodeList(self,nodelist,id=False,space=True):    ### Returns a text summary of listed nodes '(X,Y,Z)' etc.
        '''
        Returns a text summary of listed nodes '(X,Y,Z)' etc.
        >> nodelist:list of Node Objects
        >> id:Boolean [False] = whether to use node IDs (True) or names (False)
        >> space:Boolean [True] = whether to have a space after commas for clarity
        << sumtxt:str shortnames of listed nodes (X,Y,Z...)
        '''
        try:
            sumtxt = []
            for node in nodelist:
                if id: sumtxt.append('%d' % node.stat['ID'])
                else: sumtxt.append(node.shortName())
            if space: return '(%s)' % rje.join(sumtxt,', ')
            else: return '(%s)' % rje.join(sumtxt,',')
        except:
            self.log.errorLog('Problem with nodeList()')
            raise
#########################################################################################################################
    def _getNode(self,id):  ### Returns Node Object with given ID
        '''
        Returns Node Object with given ID.
        >> id:int = node.stat['ID']
        << node:Node Object
        '''
        if self.node[id-1].stat['ID'] == id: return self.node[id-1]
        else:
            for node in self.node:
                if node.stat['ID'] == id: return node
        return None
#########################################################################################################################
    def _makeNSFTree(self,seqnum=False,seqname='short',maxnamelen=123,blen='Length',bootstraps='boot',multiline=True,te=True):
        """
        Generates an NSF Tree from tree data.
        >> seqnum:boolean [False] = whether to print sequence numbers
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names;  
        >> maxnamelen:int [123] = truncate names to this length for compatability with other programs.
        >> blen:str ['Length'] = stat to use for branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> bootstraps:str ['boot'] = what to print for bootstraps
        - 'none' = nothing; 'boot' = boostraps if given (else none); 'node' = node numbers
        >> multiline:boolean [True] = whether to spread file over multiple lines (1 per node)
        >> te:boolean [True] = make TreeExplorer compatible
        - whitespaces are replaced with underscores, brackets with '' and colons/semicolons with -
        << outtree:str = NSF tree
        """
        try:### ~ [0] ~ Seed Tree with root or trichotomy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nodes = self.node[-1].neighbours()
            outtree = '%s;' % self.nodeList(nodes,id=True,space=False)
            nodeout = [self.node[-1]] 
            ### ~ [1] ~ Expand internal nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            n = 2
            while self.node[-n].stat['ID'] > self.stat['SeqNum']:    # Still nodes to be expanded
                node = self.node[-n]
                re_node = re.compile('([\(,\n])%d([\),])' % node.stat['ID'])
                if re_node.search(outtree):
                    ## ~ [1a] ~ New Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    nodeout.append(node)
                    nodes = node.neighbours(ignore=nodeout)
                    rep = self.nodeList(nodes,id=True,space=False)
                    branch = node.ancBranch()
                    if (bootstraps == 'boot') and self.opt['Bootstrapped']: rep += '%d' % branch.stat['Bootstrap']
                    elif bootstraps == 'node': rep += '%d' % node.stat['ID']
                    if blen == 'fix': rep += ':%f' % self.stat['DefLen']
                    elif blen.lower() == 'pam': rep += ':%f' % (float(branch.stat['PAM']) / 100)
                    elif blen.lower() != 'none': rep += ':%f' % branch.stat[blen]
                    ## ~ [1b] ~ Replacement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    m_node = re_node.search(outtree)
                    g_node = m_node.groups()
                    rep = '%s%s%s' % (g_node[0],rep,g_node[1])
                    if multiline == 1: rep += '\n'
                    outtree = re_node.sub(rep,outtree)
                    n += 1
                else:
                    self.log.errorLog("Node %d missing from tree: %s" % (node.stat['ID'],outtree))
                    raise
            ### ~ [2] ~ Termini ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in self.node:
                if node.stat['ID'] > self.stat['SeqNum']: continue    # Internal
                else:
                    re_node = re.compile('([\(,\n])%d([\),])' % node.stat['ID'])
                    if re_node.search(outtree):
                        ## ~ [2a] ~ Name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        rep = ''
                        if seqnum:
                            rep = rje.preZero(node.stat['ID'],self.nodeNum())
                            if seqname != 'num': rep += ' '
                        if seqname == 'short': rep += node.shortName()
                        elif seqname == 'long': rep += node.info['Name']
                        if te:
                            rep = re.sub('\s','_',rep)
                            rep = re.sub('[\(\)\[\]]','"',rep)
                            rep = re.sub('[:;]','-',rep)
                            rep = re.sub(',','_-',rep)
                        else:
                            rep = re.sub('\(','{',rep)
                            rep = re.sub('\)','}',rep)
                            rep = re.sub('[:;]','-',rep)
                            rep = re.sub(',',' -',rep)
                        if len(rep) > maxnamelen and maxnamelen > 3:
                            mx = maxnamelen - 3
                            rep = '%s...' % (rep[:mx])
                        ## <ii> ##  Branch    
                        branch = node.branch[0]
                        if blen == 'fix': rep += ':%f' % self.stat['DefLen']
                        elif blen.lower() == 'pam': rep += ':%f' % (float(branch.stat['PAM']) / 100)
                        elif blen.lower() != 'none': rep += ':%f' % branch.stat[blen]
                        m_node = re_node.search(outtree)
                        g_node = m_node.groups()
                        rep = '%s%s%s' % (g_node[0],rep,g_node[1])
                        if multiline == 1: rep += '\n'
                        outtree = re_node.sub(rep,outtree)
                    else:
                        self.log.errorLog("Node %d missing from tree: %s" % (node.stat['ID'],outtree))
                        raise
            return outtree
        except:
            self.log.errorLog("Problems making tree in _makeNSFTree()!")
            raise
#########################################################################################################################
    def saveTrees(self,seqname='long',blen='Length',bootstraps='boot'):    ### Generates all tree file formats selected in normal format
        '''
        Generates all tree file formats selected in normal format.
        >> seqname:str ['long'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names; 
        - whitespaces are replaced with underscores, brackets with '' and colons/semicolons with -
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> bootstraps:str ['boot'] = what to print for bootstraps
        - 'none' = nothing; 'boot' = bootstraps if given (else none); 'node' = node numbers
        '''
        try:### ~ [1] ~ Try each tree format in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Replace type!
            for type in self.list['TreeFormats']:
                treefile = '%s.%s' % (self.info['Basefile'],formatext[type])
                if type == 'text':
                    self.saveTree(filename=treefile,type=type,seqnum=True,seqname=seqname,blen=blen,bootstraps=bootstraps,multiline=True,maxnamelen=0)
                elif type == 'te': self.saveTree(filename=treefile,type=type,seqnum=False,seqname=seqname,blen=blen,bootstraps=bootstraps,multiline=True,maxnamelen=120)
                elif self.obj['SeqList']: self.saveTree(filename=treefile,type=type,seqnum=self.opt['SeqNum'],seqname=seqname,blen=blen,bootstraps=bootstraps,multiline=True,maxnamelen=self.stat['TruncNames'])
                else: self.saveTree(filename=treefile,type=type,seqnum=False,seqname=seqname,blen=blen,bootstraps=bootstraps,multiline=True,maxnamelen=self.stat['TruncNames'])
        except: self.errorLog('Major problem with saveTrees.')
#########################################################################################################################
    def saveTree(self,filename='None',type='nsf',seqnum=False,seqname='short',maxnamelen=123,blen='Length',bootstraps='boot',multiline=True):
        '''
        Saves tree to a file.
        >> filename:str ['None']
        >> type:str ['nsf']
        - 'nsf' = Newick Standard Format; 'text' = text; 'r' = for r graphic, 'png' = r png, 'te' = TreeExplorer NSF
        >> seqnum:boolean [False] = whether to print sequence numbers (from zero if bootstraps = 'num')
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names; 
        - whitespaces are replaced with underscores, brackets with '' and colons/semicolons with -
        >> maxnamelen:int [123] = truncate names to this length for compatability with other programs.
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> bootstraps:str ['boot'] = what to print for bootstraps
        - 'none' = nothing; 'boot' = boostraps if given (else none); 'node' = node numbers
        >> multiline:boolean = whether to spread file over multiple lines (1 per node)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename.lower() in ['none','']: filename = self.info['SaveTree']
            if filename.lower() in ['none','']: filename = '%s.%s' % (rje.baseFile(self.info['Name'],True),type)
            if filename.lower() in ['none','']: filename = 'out.nsf'
            if seqname == 'num': seqnum = True
            try: self.printLog('\r#TREE','Saving tree to %s as %s.' % (filename,treeformats[type]))
            except: self.printLog('\r#TREE','Saving tree to %s in %s format.' % (filename,type))
            self.log.printLog("#OUT","filename='%s', type='%s', seqnum=%s, seqname='%s', maxnamelen=%d, blen='%s', bootstraps='%s'\n" % (filename,type,seqnum,seqname,maxnamelen,blen,bootstraps))
            ### ~ [2] ~ Generate Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ NSF Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
            if type in ['nsf','te','nwk']: outtree = self._makeNSFTree(seqnum,seqname,maxnamelen,blen,bootstraps,multiline,te=type=='te')
            ## ~ [2b] ~ Plain text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif type == 'text': return self.textTree(filename,seqnum,seqname,maxnamelen,blen=blen,showboot=bootstraps=='boot',compress=False)
            ## ~ [2c] ~ R Graphic  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif type == 'r': return self.rTree(filename,seqname,blen=blen,compress=False)
            elif type in ['png','bud','qspec','cairo']: return self.pngTree(filename,seqname,blen=blen,compress=False,type=type)
            ## ~ [2d] ~ SVG Graphic  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif type == 'svg': return self.svgTree(filename,seqname,blen=blen,compress=False,save=True)
            elif type == 'html': return self.htmlTree(filename,seqname,blen=blen,compress=False)
            ## ~ [2x] ~ Unsupported Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                self.log.errorLog('Sorry. Output type %s not supported.' % type,printerror=False)
                raise ValueError
            ### ~ [3] ~ Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if multiline: outtree = rje.replace(outtree,'(','(\n')
            open(filename, 'w').write(outtree)
        except: self.errorLog('Major problem with saveTree(%s).' % type)
#########################################################################################################################
    def editTree(self,groupmode=False,fromnode=None,reroot=True):     ### Options to alter tree details and compress nodes etc.
        '''
        Options to alter tree details and compress nodes etc.
        - Expand/Collapse nodes
        - Flip Clade
        - Zoom to Clade (show subset of tree)
        - Prune Tree (permanently delete branches and nodes)
        - Edit Branch
        - Edit Terminal Node
        - ReRoot/Unroot
        - Display Options
        - Save Text Tree
        >> groupmode:Boolean [False] = whether node compression defines groups
        >> fromnode:Node Object [None] = start by drawing tree from given Node only (Zoom)
        >> reroot:Boolean [True] = whether to allow rerooting of tree
        '''
        try:
            ### <a> ### Set Options 
            seqname = 'short'
            maxnamelen = 0
            showboot = True
            showlen = 'branch'
            showdesc = False    # Show descriptions
            blen = 'Length'
            scale = self.stat['TextScale']
            spacer = 1
            pause = 50
            ### <b> ### Choices 
            while self.stat['Interactive'] >= 0:
                self.textTree(filename=None,seqnum=True,seqname=seqname,maxnamelen=maxnamelen,nodename='num',showboot=showboot,showlen=showlen,blen=blen,scale=scale,spacer=spacer,pause=pause,fromnode=fromnode)
                if groupmode:
                    rje.printf('\nDisplay: <C>ollapse/Group, <E>xpand/Ungroup, ',newline=False)
                else:
                    rje.printf('\nDisplay: <C>ollapse, <E>xpand, ',newline=False)
                rje.printf('<F>lip, <N>ame Clade, <Z>oom, <O>ptions, <D>escriptions [%s], <S>ave to file, <Q>uit' % (showdesc))
                rje.printf('Edit: <B>ranch, <T>erminal node, <P>rune',newline=False)
                if reroot:
                    print(', <R>oot on Branch, <U>nroot, Root <M>enu')
                else:
                    print()
                choice = rje.choice('Choice?: ',default='Q').upper()
                while not reroot and choice in ['R','U','M']:
                    choice = rje.choice('Choice?: ',default='Q').upper()
                neednode = ['C','E','F','N','T']
                needbranch = ['B','P','R']
                branch = None
                ## Add automatic collapse/expand if number given
                try:
                    n = rje.atoi(choice) # Will raise ValueError if choice not an integer string
                    node = self._getNode(n)
                    if node.opt['Compress']:
                        choice = 'E'
                    else:
                        choice = 'C'
                except: 
                    node = None
                #print choice, node
                ## Choose node/branch if appropriate
                while choice[0] in neednode and node == None:   # Choose Node
                    n = rje.getInt('Node: ')
                    node = self._getNode(n)
                    if node == None:
                        print('Cannot find Node %d!' % n)
                while choice[0] in needbranch and branch == None:   # Choose Branch
                    n = rje.getInt('Descendant Node: ')
                    node = self._getNode(n)
                    if node == None:
                        print('Cannot find Node %d!' % n)
                    else:
                        branch = node.ancBranch()
                ## <i> ## Collapse/Name node
                if choice.find('C') == 0 or choice.find('N') == 0:
                    if choice.find('C') == 0:
                        node.opt['Compress'] = True
                        if groupmode:
                            self._addGroup(node)
                    if node.info['CladeName'] == 'None':
                        clades = self._descClades(node)
                        clades = clades[0] + clades[1]
                        node.info['CladeName'] = '%d Seqs %s' % (len(clades),self.nodeList(clades))
                    #print 'Clade Name (Blank to Keep)',node.info['CladeName']
                    node.info['CladeName'] = self._editChoice(text='Clade Name (Blank to Keep)',value=node.info['CladeName'])
                ## <ii> ## Expand node
                if choice.find('E') == 0: 
                    node.opt['Compress'] = False
                    if groupmode and node in self.subfam:
                        self.subfam.remove(node)
                ## <iii> ## Flip node
                if choice.find('F') == 0: 
                    node.flipDesc()
                ## <iv> ## Options
                if choice.find('D') == 0:
                    showdesc = not showdesc
                    if showdesc:
                        seqname = 'long'
                    else:
                        seqname = 'short'
                if choice.find('O') == 0:
                    print(self.textTree.__doc__)
                    seqname = self._editChoice('Sequence name display (num/short/full)',seqname)
                    showdesc = False
                    if seqname == 'long':
                        showdesc = True
                    maxnamelen = int(self._editChoice('Maximum Name Length',maxnamelen,numeric=True))
                    showboot = self._editChoice('Show Bootstraps',showboot,boolean=True)
                    showlen = self._editChoice('Show Lengths',showlen)
                    blen = self._editChoice('Branch Lengths',blen)
                    scale = int(self._editChoice('Scale',scale,numeric=True))
                    spacer = int(self._editChoice('Spacer',spacer,numeric=True))
                    pause = int(self._editChoice('Pause',pause,numeric=True))
                ## <v> ## Save to file
                if choice.find('S') == 0: 
                    filename = rje.choice('Filename: ',default=self.info['Name'])
                    if rje.yesNo('Save text tree as %s?' % filename):
                        self.textTree(filename=filename,seqnum=True,seqname=seqname,maxnamelen=maxnamelen,nodename='num',showboot=showboot,showlen=showlen,blen=blen,scale=scale,spacer=spacer,pause=pause,fromnode=fromnode)
                ## <vi> ## Edit Branch
                if choice.find('B') == 0: 
                    branch.edit()
                ## <vii> ## Save to file
                if choice.find('T') == 0: 
                    print('Currently Unavailable!')
                ## <viii> ## Prune Tree
                if choice.find('P') == 0:
                    if rje.yesNo('Remove %d and all descendant nodes and branches?' % n):
                        self._prune(branch)
                ## <ix> ## Root on Branch
                if choice.find('R') == 0:
                    self.unRoot()
                    self.placeRoot('Manual',branch)
                ## <x> ## Unroot
                if choice.find('U') == 0: 
                    self.unRoot()
                ## <xi> ## Root Menu
                if choice.find('M') == 0:
                    self.info['Rooting'] = 'man'
                    self.treeRoot()
                ## <xii> ## Quit
                if choice.find('Q') == 0 and rje.yesNo('Quit Tree Edit?'): # Quit
                    return                    
                ## <xiii> ## Zoom
                if choice.find('Z') == 0:
                    n = rje.getInt('Node (0 to Zoom Out): ')
                    if n == 0:
                        fromnode = None
                    else:
                        fromnode = self._getNode(n)
                        if fromnode == None:
                            print('Cannot find Node %d!' % n)
                        elif fromnode == self.node[-1]:
                            fromnode = None
                    
#        except KeyboardInterrupt:
#            if rje.yesNo('Quit Tree Edit?'): # Quit
#                return
#            else:
#                self.editTree()
        except:
            self.log.errorLog('Major Problem during editTree().',True)
#########################################################################################################################
    def textTree(self,filename=None,seqnum=True,seqname='short',maxnamelen=50,nodename='num',showboot=True,showlen='branch',blen='Length',scale=-1,spacer=1,pause=50,fromnode=None,compress=True):
        '''
        Outputs tree as ASCII text.
        >> filename:str [None]
        >> seqnum:boolean [True] = whether to print sequence numbers (from zero if nodename <> 'none')
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names
        >> maxnamelen:int [50] = truncate names to this length 
        >> nodename:str ['num'] = how to label nodes
        - 'none' = no label, 'num' = numbers only, 'short' = short names, 'long' = long names
        >> showboot:boolean [True] = whether to show boostraps if given
        >> showlen:str ['branch'] = whether to show lengths of branches leading to nodes
        - 'none' = no lengths, 'branch' = branch lengths, 'sum' = summed length from root
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> scale:int [-1] = no. of characters per self._deflen distance
        >> spacer:int [1] = 'blank' rows between nodes
        >> pause:int [50] = Number of lines printed before pausing for user <ENTER>
        >> fromnode:Node Object [None] = draw subtree from Node
        >> compress:boolean = whether to compress compressed nodes
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Setup Scale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.stat['DefLen'] <= 0.0:
                self.printLog('#DEF','Default branch length cannot be <= 0.0 - reset to 0.1')
                self.stat['DefLen'] = 0.1
            xpos = {}       # Dictionary of node positions on x-axis
            if scale <= 0: scale = self.stat['TextScale']
            if scale <= 0:
                self.log.errorLog('Scale must be > 0! Changed to 1.',printerror=False)
                scale = 1
                self.stat['TextScale'] = 1
            ## ~ [0b] ~ Generate PAM BranchLengths if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (blen == 'pam') and (self.branch[0].stat['PAM'] < 0): self.branchPam()

            ### ~ [1] ~ Establish Base of Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ignore = []     # List of nodes to ignore when expanding tree (not in tree or already expanded)
            tric = False    # Whether base of tree is a trichotomy
            if not fromnode or fromnode == self.node[-1]:
                basenode = self.node[-1]
                if len(basenode.branch) == 3: tric = True   # Trichotomy
                xpos[basenode] = scale + 1
            else:
                basenode = fromnode
                anc = basenode.ancNode()
                ignore.append(anc)
                xpos[anc] = 1
                branch = basenode.ancBranch()
                if blen == 'fix': mylen = self.stat['DefLen']
                elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                else: mylen = branch.stat[blen]
                xpos[basenode] = int(mylen * scale / self.stat['DefLen']) + 2 + xpos[anc]

            ### ~ [2] ~ Make a list of order of node output and calculate X-pos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ynode = [basenode]      # List of ynodes in order
            replacements = True     # Controls continuing of loop until no node replaced with clade
            while replacements:     # Loop while replacements being made
                replacements = False
                newy = []           # New ynode list to replace ynode
                for node in ynode:
                    ## ~ [2a] ~ Check if node already replaced or terminus/outside subtree ~~~~~~~~ ##
                    if node in ignore or node.stat['ID'] <= self.stat['SeqNum']:
                        newy.append(node)
                        continue
                    ## ~ [2b] ~ Process node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    replacements = True     # Keep Looping
                    ignore.append(node)     # Ignore this node next time
                    if node.opt['Compress'] and compress: newy.append(node)   # Compress clade to node (do not expand)
                    elif node.stat['ID'] > self.stat['SeqNum']:
                        ## ~ [2c] ~ Calculate ypos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        nodes = node.neighbours(ignore=ignore)
                        newy += [nodes[0],node] + nodes[1:]
                        ignore.append(node)
                        ## ~ [2d] ~ Calculate xpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        for i in range(len(nodes)):
                            branch = nodes[i].link(node)
                            if blen == 'fix': mylen = self.stat['DefLen']
                            elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                            else: mylen = branch.stat[blen]
                            xpos[nodes[i]] = int(mylen * scale / self.stat['DefLen']) + 2 + xpos[node]
                    else: raise ValueError
                ynode = newy        # Update ynode list

            ### ~ [3] ~ Generate Text Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maxpos = 0      # Establish maximum x position of node
            vline = {}      # Dictionary used for drawing vertical lines
            for node in ynode:
                vline[node] = False    # List of on/off values for vertical line
                maxpos = max(maxpos,xpos[node])
            treelines=[]    # List of text lines to be output
            ## ~ [3a] ~ Build lines, node by node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for node in ynode:
                branch = node.ancBranch()
                ## ~ [3b] ~ Branch leading to node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not branch:      # Root/Trichotomy
                    line = '+'      # Start of branch
                    while len(line) < (xpos[node] - 1): line += '-'
                else:
                    anc = branch.link(node)
                    line = ' ' * (xpos[anc] - 1)    # Space before Branch
                    if anc.opt['Duplication']:      # Duplication-specific branch characters
                        line += '#'     
                        while len(line) < (xpos[node] - 1): line += '='
                    else:
                        if anc.opt['SpecDup']: line += '*'
                        else: line += '+'
                        while len(line) < (xpos[node] - 1): line += '-'
                if node.opt['Duplication']: line += '#'     # End of branch
                elif node.opt['SpecDup']: line += '*'
                else: line += '+'
                ## ~ [3c] ~ Node Name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                txtname = ''
                if node.stat['ID'] <= self.stat['SeqNum']:      # Terminal node
                    if seqnum or (seqname == 'None') or (seqname == 'num'):
                        txtname = '%d' % node.stat['ID']
                        if seqname != 'num' and seqname != 'None': txtname += ': '
                    if seqname == 'short': txtname += node.shortName()
                    elif seqname == 'long': txtname += node.info['Name']
                    if len(txtname) > maxnamelen and maxnamelen > 0: txtname = txtname[:maxnamelen] + '...'                
                else:                                           # Internal node
                    if nodename == 'none': txtname = ''
                    elif node.opt['Compress'] and compress:
                        if nodename == 'num': txtname = '*%d*: ' % node.stat['ID']
                        txtname += node.info['CladeName']
                    elif nodename == 'num': txtname = '%d' % node.stat['ID']
                    elif nodename == 'short': txtname += node.shortName()
                    elif nodename == 'long': txtname += node.info['Name']
                ## ~ [3d] ~ Bootstrap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if showboot and self.opt['Bootstrapped'] and branch and node.stat['ID'] > self.stat['SeqNum']:
                    txtname += ' [%d]' % branch.stat['Bootstrap']
                ## ~ [3e] ~ Branch lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if branch:
                    if showlen == 'branch':
                        if blen == 'fix': mylen = self.stat['DefLen']
                        elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                        else: mylen = branch.stat[blen]
                        txtname += ' {%f}' % mylen
                    elif showlen == 'sum': txtname += ' {%f}' % self.pathLen(self.rootPath(node))
                line += ' ' + txtname
                ## ~ [3f] ~ Add vertical lines preceeding clade ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for vnode in ynode:
                    if vline[vnode] and vnode != anc and vnode != node:
                        vx = xpos[vnode] - 1
                        line = line[:vx] + '|' + line[(vx+1):]
                if branch != None and node != fromnode:
                    vline[anc] = not vline[anc]
                    if anc == basenode and tric and vline[anc] == False:
                        vline[anc] = True
                        tric = False
                treelines.append('%s\n' % line)
                ## ~ [3g] ~ Add spacer lines between node lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for y in range(spacer):
                    line = ' ' * xpos[node]
                    for vnode in ynode:
                        if vline[vnode]:
                            vx = xpos[vnode] - 1
                            line = line[:vx] + '|' + line[(vx+1):]
                treelines.append('%s\n' % line)

            ### ~ [4] ~ Save and/or Print Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'd'
            if filename and filename.lower() != 'none':
                try:
                    open(filename, 'w').write(rje.join(treelines,''))
                    self.log.printLog("#OUT","Text tree saved to %s" % filename,1)
                except:
                    self.log.errorLog("Problem saving text tree to %s" % filename)
                    raise
            else:
                p = 1
                self.verbose(0,5,'\n',1)
                for line in treelines:
                    if p < pause:
                        self.verbose(0,5,line,0)
                        p += 1
                    else:
                        self.verbose(0,1,line,0)
                        p = 0
            return rje.join(treelines,'')
        except: self.errorLog('Major Problem with Tree.textTree()')
#########################################################################################################################
    def rTree(self,filename,seqname='short',nodename='long',blen='Length',fromnode=None,compress=True,title=None,qry='qry'):
        '''
        Outputs details of tree to a file to be read in and processed by R.
        >> filename:str [None]
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names
        >> nodename:str ['num'] = how to label nodes
        - 'none' = no label, 'num' = numbers only, 'short' = short names, 'long' = long names
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> fromnode:Node Object [None] = draw subtree from Node
        >> compress:boolean = whether to compress compressed nodes
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xpos = {}       # Dictionary of node positions on x-axis
            if qry == 'qspec': qseq = self.obj['SeqList'].obj['QuerySeq']
            else: qseq = None
            ## ~ [0a] ~ Generate PAM BranchLengths if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (blen == 'pam') and (self.branch[0].stat['PAM'] < 0): self.branchPam()
            ## ~ [0b] ~ Setup file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = ['nodenum','ypos','xpos','anc','ancy','ancx','name','boot','family']
            rje.delimitedFileOutput(self,filename,headers,',',rje_backup=True)
            ## ~ [0c] ~ Setup families ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            subfams = {}
            for subfam in self.subfams():
                famname = subfam.info['CladeName']
                if rje.matchExp('(\d+) Seqs \(',famname): famname = 'SubFam%d' % self.subfams().index(subfam) 
                subfams[subfam] = famname
                for node in self._nodeClade(subfam,internal=True): subfams[node] = famname
            for node in self.nodes():
                if rje.split(node.info['Name'],'_')[0] == qry: subfams[node] = qry
                elif node.obj['Sequence'] and qseq and node.obj['Sequence'].sameSpec(qseq): subfams[node] = 'qry'
            #self.debug(subfams)
                    
            ### ~ [1] ~ Establish Base of Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ignore = []     # List of nodes to ignore when expanding tree (not in tree or already expanded)
            tric = False    # Whether base of tree is a trichotomy
            if not fromnode or fromnode == self.node[-1]:
                basenode = self.node[-1]
                if len(basenode.branch) == 3: tric = True   # Trichotomy
                xpos[basenode] = 0
            else:
                basenode = fromnode
                anc = basenode.ancNode()
                ignore.append(anc)
                xpos[anc] = 0
                branch = basenode.ancBranch()
                if blen == 'fix': mylen = self.stat['DefLen']
                elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                else: mylen = branch.stat[blen]
                xpos[basenode] = max(0,mylen) + xpos[anc]

            ### ~ [2] ~ Make a list of order of node output and calculate X-pos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ynode = [basenode]      # List of ynodes in order
            replacements = True     # Controls continuing of loop until no node replaced with clade
            while replacements:     # Loop while replacements being made
                replacements = False
                newy = []           # New ynode list to replace ynode
                for node in ynode:
                    ## ~ [2a] ~ Check if node already replaced or terminus/outside subtree ~~~~~~~~ ##
                    if node in ignore or node.stat['ID'] <= self.stat['SeqNum']:
                        newy.append(node)
                        continue
                    ## ~ [2b] ~ Process node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    replacements = True     # Keep Looping
                    ignore.append(node)     # Ignore this node next time
                    if node.opt['Compress'] and compress: newy.append(node)   # Compress clade to node (do not expand)
                    elif node.stat['ID'] > self.stat['SeqNum']:
                        ## ~ [2c] ~ Calculate ypos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        nodes = node.neighbours(ignore=ignore)
                        newy += [nodes[0],node] + nodes[1:]
                        ignore.append(node)
                        ## ~ [2d] ~ Calculate xpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        for i in range(len(nodes)):
                            branch = nodes[i].link(node)
                            if blen == 'fix': mylen = self.stat['DefLen']
                            elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                            else: mylen = max(0,branch.stat[blen])
                            xpos[nodes[i]] = mylen + xpos[node]
                    else: raise ValueError
                ynode = newy        # Update ynode list

            ### ~ [3] ~ Output Tree Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in ynode:
                data = {'nodenum':node.stat['ID'],'ypos':ynode.index(node),'xpos':xpos[node],'name':node.info['Name']}
                if node == self.node[-1]:
                    data['ancx'] = -0.1
                    data['ancy'] = data['ypos']
                else:
                    anc = node.ancNode()
                    data['anc'] = anc.stat['ID']
                    data['ancx'] = xpos[anc]
                    data['ancy'] = ynode.index(anc)
                if node in subfams: data['family'] = subfams[node]
                #self.debug(data)
                data['boot'] = ''
                if self.opt['Bootstrapped'] and node.stat['ID'] > self.stat['SeqNum']:
                    branch = node.ancBranch()
                    if branch: data['boot'] = branch.stat['Bootstrap']
                if node.stat['ID'] <= self.stat['SeqNum']:      # Terminal node
                    if seqname == 'num': data['name'] = node.stat['ID']
                    if seqname == 'short': data['name'] = node.shortName()
                    elif seqname == 'long': data['name'] = node.info['Name']
                    if not qseq and node.obj['Sequence'] and node.obj['Sequence'].info['AccNum'][:1] == 'E' and node.obj['Sequence'].info['SpecCode'] == 'EMIHU':
                        #if node.obj['Sequence'].info['Description'] == 'unnamed protein product':
                        #    data['name'] = '*** %s Consensus EST translation %s *** ' % (node.shortName(),node.obj['Sequence'].info['AccNum'])
                        #else: data['name'] = '*** %s *** ' % data['name'] 
                        data['family'] = 'EHUX'
                else:                                           # Internal node
                    #data['family'] = ''    # Why???
                    if nodename == 'none': txtname = ''
                    elif node.opt['Compress'] and compress: data['name'] = node.info['CladeName']
                    elif nodename == 'num': data['name'] = node.stat['ID']
                    elif nodename == 'short': data['name'] = node.shortName()
                    elif nodename == 'long': data['name'] = node.info['Name']
                data['name'] = rje.replace(data['name'],'"',"'") 
                rje.delimitedFileOutput(self,filename,headers,',',data)
            self.printLog('#RTREE','Data for %d nodes output to %s' % (len(ynode),filename))

        except: self.log.errorLog('Major Problem with Tree.rTree()')
#########################################################################################################################
    def pngTree(self,filename,seqname='short',nodename='long',blen='Length',fromnode=None,compress=True,title=None,type='png'):
        '''
        Outputs details of tree to a file to be read in and processed by R.
        >> filename:str [None]
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names
        >> nodename:str ['num'] = how to label nodes
        - 'none' = no label, 'num' = numbers only, 'short' = short names, 'long' = long names
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> fromnode:Node Object [None] = draw subtree from Node
        >> compress:boolean = whether to compress compressed nodes
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rtreefile = '%s.r.csv' % rje.baseFile(filename)
            self.rTree(rtreefile,seqname,nodename,blen,fromnode,compress,title,type)
            if not os.path.exists(rtreefile):
                self.errorLog('PNG Tree error: CSV file not made!',printerror=False)
                raise IOError
            ### ~ [1] ~ Try to run R to generate PNG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RPATH',self.info['RPath'],log=False,screen=self.v()>0)
            #self.info['RPath'] = 'c:\\"Program Files"\\R\\R-2.6.2\\bin\\R.exe'
            #print self.info['RPath']
            rfile = {'png':'rje','bud':'budapest','qspec':'budapest','cairo':'rje_tree.cairo'}[type]
            self.info['PathR'] = rje.makePath(os.path.abspath(rje.join(rje.split(sys.argv[0],os.sep)[:-1]+[''],os.sep)))
            rtree = rje.makePath('%s%s.r' % (self.info['PathR'],rfile),wholepath=True)
            if not os.path.exists(rtree):
                self.info['PathR'] = rje.makePath(os.path.abspath(rje.join(rje.split(sys.argv[0],os.sep)[:-2]+['libraries','r',''],os.sep)))
                rtree = rje.makePath('%s%s.r' % (self.info['PathR'],rfile),wholepath=True)
            if not os.path.exists(rtree):
                self.info['PathR'] = rje.makePath(os.path.abspath(rje.join(rje.split(sys.argv[0],os.sep)[:-1]+['libraries','r',''],os.sep)))
                rtree = rje.makePath('%s%s.r' % (self.info['PathR'],rfile),wholepath=True)
            self.printLog('#%s' % type.upper(),rtree)
            if rfile == 'rje':
                rcmd = '%s --no-restore --no-save --args "tree" "%s" "rdir=%s"' % (self.info['RPath'],rje.baseFile(filename),self.info['PathR'])
                if title: rcmd += ' "treetitle=%s"' % title
                else: rcmd += ' "treetitle=%s"' % rje.baseFile(filename,strip_path=True)
            else:
                rcmd = '%s --no-restore --no-save --args "%s" "%s" "%s"' % (self.info['RPath'],rtreefile,filename,self.info['PathR'])
                if title: rcmd += ' "%s"' % title
            #rcmd += ' < "%s" > "%s.r.tmp.txt" 2>&1' % (rtree,rje.baseFile(filename))
            rcmd += ' < "%s" > "%s.r.run"' % (rtree,rje.baseFile(filename))
            #x#rcmd = rje.replace(rcmd,'\\','\\\\')
            self.printLog('#RTREE',rcmd)
            problems = os.popen(rcmd).read()
            if problems:
                self.errorLog(problems,printerror=False)
                for ptxt in problems: self.warnLog(ptxt)
            elif rje.exists(filename) and rje.exists(rtreefile) and not (self.getBool('DeBug') or self.getBool('Test')): os.unlink(rtreefile)
            if not self.getBool('DeBug') and not self.getBool('Test') and os.path.exists('%s.r.run' % rje.baseFile(filename)): os.unlink('%s.r.run' % rje.baseFile(filename))
            return rcmd
        except: self.log.errorLog('Major Problem with Tree.pngTree()')
#########################################################################################################################
    def dictTree(self,seqname='short',nodename='long',blen='Length',fromnode=None,compress=True,title=None,qry='qry'):
        '''
        Returns tree as data dictionary for SVG (and R) output.
        >> filename:str [None]
        >> seqname:str ['short'] = name to use for sequences
        - 'num' = numbers only; 'short' = short sequence names; 'long' = long names
        >> nodename:str ['num'] = how to label nodes
        - 'none' = no label, 'num' = numbers only, 'short' = short names, 'long' = long names
        >> blen:str ['Length'] = branch lengths 
        - 'none' = do not use branch lengths; 'fix' = fix at _deflen;
        - 'pam' = replace with PAM distances if calculated (else none)
        - other = use node.stat[blen]
        >> fromnode:Node Object [None] = draw subtree from Node
        >> compress:boolean = whether to compress compressed nodes
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xpos = {}       # Dictionary of node positions on x-axis
            ## ~ [0a] ~ Generate PAM BranchLengths if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (blen == 'pam') and (self.branch[0].stat['PAM'] < 0): self.branchPam()
            ## ~ [0b] ~ Setup file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datadict = {}
            headers = ['nodenum','ypos','xpos','anc','ancy','ancx','name','boot','family']
            ## ~ [0c] ~ Setup families ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            subfams = {}
            for subfam in self.subfams():
                famname = subfam.info['CladeName']
                if rje.matchExp('(\d+) Seqs \(',famname): famname = 'SubFam%d' % self.subfams().index(subfam) 
                subfams[subfam] = famname
                for node in self._nodeClade(subfam,internal=True):
                    subfams[node] = famname
                    if rje.split(node.info['Name'],'_')[0] == type: subfams[node] = type
                    
            ### ~ [1] ~ Establish Base of Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ignore = []     # List of nodes to ignore when expanding tree (not in tree or already expanded)
            tric = False    # Whether base of tree is a trichotomy
            if not fromnode or fromnode == self.node[-1]:
                basenode = self.node[-1]
                if len(basenode.branch) == 3: tric = True   # Trichotomy
                xpos[basenode] = 0
            else:
                basenode = fromnode
                anc = basenode.ancNode()
                ignore.append(anc)
                xpos[anc] = 0
                branch = basenode.ancBranch()
                if blen == 'fix': mylen = self.stat['DefLen']
                elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                else: mylen = branch.stat[blen]
                xpos[basenode] = max(0,mylen) + xpos[anc]

            ### ~ [2] ~ Make a list of order of node output and calculate X-pos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ynode = [basenode]      # List of ynodes in order
            replacements = True     # Controls continuing of loop until no node replaced with clade
            while replacements:     # Loop while replacements being made
                replacements = False
                newy = []           # New ynode list to replace ynode
                for node in ynode:
                    ## ~ [2a] ~ Check if node already replaced or terminus/outside subtree ~~~~~~~~ ##
                    if node in ignore or node.stat['ID'] <= self.stat['SeqNum']:
                        newy.append(node)
                        continue
                    ## ~ [2b] ~ Process node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    replacements = True     # Keep Looping
                    ignore.append(node)     # Ignore this node next time
                    if node.opt['Compress'] and compress: newy.append(node)   # Compress clade to node (do not expand)
                    elif node.stat['ID'] > self.stat['SeqNum']:
                        ## ~ [2c] ~ Calculate ypos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        nodes = node.neighbours(ignore=ignore)
                        newy += [nodes[0],node] + nodes[1:]
                        ignore.append(node)
                        ## ~ [2d] ~ Calculate xpos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        for i in range(len(nodes)):
                            branch = nodes[i].link(node)
                            if blen == 'fix': mylen = self.stat['DefLen']
                            elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                            else: mylen = max(0,branch.stat[blen])
                            xpos[nodes[i]] = mylen + xpos[node]
                    else: raise ValueError
                ynode = newy        # Update ynode list

            ### ~ [3] ~ Make Tree Data Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in ynode:
                data = {'nodenum':node.stat['ID'],'ypos':ynode.index(node),'xpos':xpos[node],'name':node.info['Name']}
                if node == self.node[-1]:
                    data['ancx'] = -0.1
                    data['ancy'] = data['ypos']
                else:
                    anc = node.ancNode()
                    data['anc'] = anc.stat['ID']
                    data['ancx'] = xpos[anc]
                    data['ancy'] = ynode.index(anc)
                if node in subfams: data['family'] = subfams[node]
                data['boot'] = ''
                if self.opt['Bootstrapped'] and node.stat['ID'] > self.stat['SeqNum']:
                    branch = node.ancBranch()
                    if branch: data['boot'] = branch.stat['Bootstrap']
                if node.stat['ID'] <= self.stat['SeqNum']:      # Terminal node
                    if seqname == 'num': data['name'] = node.stat['ID']
                    if seqname == 'short': data['name'] = node.shortName()
                    elif seqname == 'long': data['name'] = node.info['Name']
                    if node.obj['Sequence'] and node.obj['Sequence'].info['AccNum'][:1] == 'E' and node.obj['Sequence'].info['SpecCode'] == 'EMIHU':
                        if node.obj['Sequence'].info['Description'] == 'unnamed protein product':
                            data['name'] = '*** %s Consensus EST translation %s *** ' % (node.shortName(),node.obj['Sequence'].info['AccNum'])
                        else: data['name'] = '*** %s *** ' % data['name'] 
                        data['family'] = 'EHUX'
                else:                                           # Internal node
                    if nodename == 'none': txtname = ''
                    elif node.opt['Compress'] and compress: data['name'] = node.info['CladeName']
                    elif nodename == 'num': data['name'] = node.stat['ID']
                    elif nodename == 'short': data['name'] = node.shortName()
                    elif nodename == 'long': data['name'] = node.info['Name']
                data['name'] = rje.replace(data['name'],'"',"'")
                for h in headers:
                    if h not in data: data[h] = ''
                datadict[node.stat['ID']] = data
            return datadict
        except: self.log.errorLog('Major Problem with Tree.dictTree()')
#########################################################################################################################
    def svgTree(self,filename,seqname='short',nodename='long',blen='Length',fromnode=None,compress=True,title=None,qry='qry',save=True,treesplit=0.5,font=12,maxfont=20,width=1600,height=0,xoffset=0,yoffset=0):
        '''
        Generates a data dictionary similar to that output to a file to be read in and processed by R for SVG generation.
        See rje_tree.dictTee() and rje_svg.svgTree() for parameters.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datadict = self.dictTree(seqname,nodename,blen,fromnode,compress,title,qry)
            ### ~ [1] ~ Generate SVG Tree code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svg = rje_svg.SVG(self.log,self.cmd_list).svgTree(rje.baseFile(filename,extlist=['.svg']),datadict,treesplit,font,maxfont,width,height,save,xoffset,yoffset)
            return svg
        except: self.log.errorLog('Major Problem with Tree.svgTree()')
#########################################################################################################################
    def htmlTree(self,filename,seqname='short',nodename='long',blen='Length',fromnode=None,compress=True,title=None,qry='qry',treesplit=0.5,font=12,maxfont=20,width=1600,height=0,xoffset=0,yoffset=0):
        '''
        Generates a data dictionary similar to that output to a file to be read in and processed by R for SVG generation.
        See rje_tree.dictTee() and rje_svg.svgTree() for parameters.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            svgfile = rje.baseFile(filename,strip_path=False,extlist=['.svg']) + '.svg'
            svglink = rje.baseFile(filename,strip_path=True,extlist=['.svg']) + '.svg'
            svg = self.svgTree(filename,seqname,nodename,blen,fromnode,compress,title,qry,True,treesplit,font,maxfont,width,height,xoffset,yoffset)
            if not title: title = 'Tree "%s"' % self.info['Name']
            html = rje_html.htmlHead(title,stylesheets=[],tabber=False,frontpage=False,nobots=True)
            html += rje_svg.SVG(self.log,self.cmd_list).svgHTML(svgfile,title,svgfile)
            html += rje_html.htmlTail(tabber=False)
            open(filename,'w').write(html)
            self.printLog('#HTML','Tree %s output to %s and %s' % (self.info['Name'],svgfile,filename))
        except: self.log.errorLog('Major Problem with Tree.svgTree()')
#########################################################################################################################
    def _vertOrder(self,fromnode=None,compress=False,internal=True,namelist=False):     ### Returns vertical (tree) ordering of nodes
        '''
        Returns vertical (tree) ordering of nodes.
        >> fromnode:Node Object = root of tree (Actual root if None)
        >> compress:bool [False] = whether to compress 'compressed' nodes
        >> internal:bool [True] = whether to return internal nodes (or just leaves)
        >> namelist:bool [False] = whether to return list of names rather than node objects
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ignore = []     # List of nodes to ignore from compression etc.
            tric = False    # Whether base is trichotomy
            if not fromnode or fromnode == self.node[-1]:
                basenode = self.node[-1]                    # Start from root
                if len(basenode.branch) == 3: tric = True   # Trichotomy
            else:
                basenode = fromnode
                anc = basenode.ancNode()
                ignore.append(anc)
            ynode = [basenode]
            ### ~ [2] ~ Make a list of order of node output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            replacements = True
            while replacements:
                replacements = False
                newy = []
                for node in ynode:          # ! # Add in future 'clade' option
                    if node in ignore or node.stat['ID'] <= self.stat['SeqNum']:  
                        newy.append(node)    # Already replaced/terminus - already handled 
                    elif compress and node.opt['Compress']:  # Compress clade to node
                        replacements = True     # Keep Looping
                        newy.append(node)       # Do not expand
                        ignore.append(node)     
                    elif node.stat['ID'] > self.stat['SeqNum']:
                        replacements = True     # Keep Looping
                        nodes = node.neighbours(ignore=ignore)
                        newy += [nodes[0],node] + nodes[1:]
                        ignore.append(node)
                    else: raise ValueError
                ynode = newy
            ### ~ [3] ~ Tidy up if desired ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in ynode[0:]:
                if not internal and node.stat['ID'] > self.stat['SeqNum']: ynode.remove(node)
                elif namelist: ynode[ynode.index(node)] = node.info['Name']
            return ynode
        except: self.errorLog('Fatal Error with _vertOrder().')
#########################################################################################################################
    def _regenerateSeqList(self,seqlist=None,nodelist=[],id=False):   ### Reorders/resizes SeqList to match given nodelist
        '''
        Reorders/resizes SeqList to match given nodelist.
        >> seqlist:rje_seq.SeqList Object
        >> nodelist:list of Node Objects, which may or may not have Sequence objects associated with them.
        >> id:boolean [False] = Whether to append Node ID to front of name (for AncSeqOut)
        '''
        try:
            if seqlist == None:
                self.log.errorLog('Called _regenerateSeqList without SeqList!',printerror=False)
                raise
            seqlist.seq = []
            
            for node in nodelist:
                seq = node.obj['Sequence']
                if seq != None:
                    seqlist.seq.append(seq)
                    if id:
                        seq.info['Name'] = '%s %s' % (rje.preZero(node.stat['ID'],self.nodeNum()),seq.info['Name'])
        except:
            self.log.errorLog('Major Problem with _regenerateSeqList().')
            raise
#########################################################################################################################
    def ancSeqOut(self,file=None,type='fas',ordered=False):    ### Output of ancestral sequences.
        '''
        Output of ancestral sequences.
        Makes a new SeqList object and outputs it.
        >> file:str [None] = filename
        >> type:str ['fas'] = type of file
        >> ordered:int [0] = whether to order nodes after seqs. (Default intersperses.)
        '''
        try:
            if file == None:
                self.verbose(0,1,"No file output.",2)
                return 0
            seqlist = self.obj['SeqList']
            ### <a> ### Determine output order
            self.verbose(0,1,"Saving Ancestral Sequences in %s..." % file,1)
            if ordered:     # Use order given in self.node = existing SeqList
                self._regenerateSeqList(seqlist,self.node,id=True)
            else:
                self._regenerateSeqList(seqlist,self._vertOrder(),id=True)
            ### <b> ### Output
            seqlist.saveFasta(seqfile=file) #!# Change this to saveSeqList when this exists!!
            self._regenerateSeqList(seqlist,self.node)
            for seq in seqlist.seq:
                if rje.matchExp('^(\d+\s+)\S',seq.info['Name']):
                    try: seq.info['Name'] = rje.replace(seq.info['Name'],rje.matchExp('^(\d+\s+)\S',seq.info['Name'])[0],'',1)
                    except: seq.info['Name'] = seq.info['Name'].replace(rje.matchExp('^(\d+\s+)\S',seq.info['Name'])[0],'',1)
        except:
            self.log.errorLog('Major Problem with ancSeqOut().')
#########################################################################################################################
    ### <5> ### Tree Rooting/Editing                                                                                    #
#########################################################################################################################
    def _getRootNode(self):  ### Returns root node or None if not rooted
        '''Returns root node or None if not rooted.'''
        try:
            root = None
            for node in self.node:
                node.setType()
                if node.info['Type'].find('Root') == 0:
                    root = node
            return root
        except:
            self.log.errorLog('Problem with _getRootNode()')
            raise
#########################################################################################################################
    def treeRoot(self):     ### Tree Rooting loop.
        '''Tree Rooting loop.'''
        while self._changeRoot():
            self.opt['ReRooted'] = True     # Distinguishes loaded root from method rooting
            self.verbose(1,2," => Root Changed!",1)
            self.textTree()
#########################################################################################################################
    def _changeRoot(self):  ### Gives options to Re-root Tree
        '''
        Gives options to Re-root (or unroot) Tree.
        Returns True if root changes or False if no change.
        Should call treeRoot() from other modules - will loop this method as appropriate.
        '''
        try:
            self.sumTree()
            rooting = self.info['Rooting']
            self.verbose(1,2,'(Root=%s. Rooted:%s. ReRooted:%s.)' % (rooting,self.opt['Rooted'],self.opt['ReRooted']),1)
            ## <a> ## Want unrooted tree
            if rooting.lower() in ['no','none','unrooted']:
                if self.opt['Rooted']: return self.unRoot()
                else: return False
            ## <b> ## First time rooting => automatic options
            elif self.opt['ReRooted'] == False:
                if rooting.lower() == 'mid':       # Midpoint root
                    return self.midRoot()    
                elif rooting.lower() == 'ran':     # Random root, no branch weighting
                    return self.randomRoot()    
                elif rooting.lower() == 'ranwt':   # Random root, branch weighting
                    return self.randomWtRoot()    
                elif rooting.lower() == 'yes':    
                    if self.opt['Rooted'] == False:    # Must root but give choice
                        if self.stat['Interactive'] >= 0:
                            return self.manRoot()  # Give manual choice
                        else:
                            return self.midRoot()  # If no interaction, must root somehow!
                    elif self.stat['Interactive'] > 0:
                        return self.manRoot()      # Give manual choice    # ! # Change to >1 and have a second menu for i=1   
                    else:   # Rooted & no need to ask
                        return False
                elif rooting.lower() == 'man':    # Ask whether to root
                    if self.stat['Interactive'] >= 0:      # Option to re-root tree
                        return self.manRoot()    
                    else: return False
                elif os.path.exists(rooting):   # Root from outgroup file
                    return self.fileRoot(rooting)
                else:
                    self.printLog('#ROOT','Rooting option "%s" unrecognised and not found as file. No re-rooting.' % rooting)
                    return False
            ## <c> ## Re-rooting options
            else:
                if self.stat['Interactive'] >= 1:      # Option to re-root tree
                    return self.manRoot()
                elif (rooting == 'man') & (self.stat['Interactive'] >= 0):      # Option to re-root tree
                    return self.manRoot()    
                else: return False
        except:
            self.log.errorLog('Major Problem with _changeRoot().')
            raise
#########################################################################################################################
    def unRoot(self):   ### UnRoots the tree.
        ''''Unroots' the tree => Changes node number and branches.'''
        try:
            ### <a> ### Setup
            root = self._getRootNode()
            if root == None:
                self.log.errorLog('Trying to unroot tree but no root found!',printerror=False)
                return False
            self.log.printLog('#ROOT','Unrooting Tree.')
            tric = root.branch[0].link(root)
            self.verbose(1,3,"Remove node %d and combine branches %s and %s." % (root.stat['ID'],root.branch[0].show(),root.branch[1].show()),1)
            ### <b> ### Combine branches
            root.branch[1].combine(root,root.branch[0])
            tric.branch.remove(root.branch[0])
            tric.branch.append(root.branch[1])
            self.branch.remove(root.branch[0])
            self.node.remove(root)
            self.opt['Rooted'] = False
            self.info['RootMethod'] = 'None'
            self._renumberNodes()
            if self._checkTree() == False:
                raise ValueError
            self.findDuplications(duptext='Tree unrooted')
            return True
        except:
            self.log.errorLog("Big Problem unrooting Tree!")
            raise
#########################################################################################################################
    def midRoot(self,fixlen=0.1):       ### Midpoint roots tree
        '''
        Midpoint roots tree.
        Returns True if successful and False if fails.
        '''
        ### <a> ### Prepare for Rooting
        try:
            if self.opt['Rooted']: self.unRoot()                    
            #X#self.verbose(0,3,'\nMidpoint rooting tree...',1)
            maxdis = 0
            extremetip = (None,None)
            rootpath = []
            max_prog = (self.nodeNum() - 1) * self.nodeNum()
            ### <b> ### Find extreme tips (furthest distance)
            for n0 in range(self.nodeNum()-1):
                for n1 in range(n0 + 1,self.nodeNum()):
                    node = (self.node[n0],self.node[n1])
                    if node[0].stat['ID'] <= self.stat['SeqNum'] and node[1].stat['ID'] <= self.stat['SeqNum']:
                        path = self.pathLink(node[0],node[1])
                        dis = self.pathLen(path,'Length')
                        if dis > maxdis:
                            maxdis = dis
                            extremetip = (node[0],node[1])
                            rootpath = path[0:]
                            #X#self.verbose(0,3,'%d..%d = %f' % (node[0].stat['ID'],node[1].stat['ID'],maxdis),0)
                    progtxt = 'Midpoint root - Establishing Extreme Tips %.1f%%' % (100.0*(n0*self.nodeNum()+n1)/max_prog)
                    if extremetip[0] and extremetip[1]:
                        self.log.printLog('\r#ROOT','%s: MaxDis = %f, %3d..%3d' % (progtxt,maxdis,extremetip[0].stat['ID'],extremetip[1].stat['ID']),log=False,newline=False)
                    else:
                        self.log.printLog('\r#ROOT','%s: MaxDis = %f, None' % (progtxt,maxdis),log=False,newline=False)
            progtxt = 'Midpoint root - Established Extreme Tips 100.0%%'
            self.log.printLog('\r#ROOT','%s: MaxDis = %f, %3d..%3d' % (progtxt,maxdis,extremetip[0].stat['ID'],extremetip[1].stat['ID']))
        except:
            self.log.errorLog('Problem establishing Extreme Tip Distance.')
            return False
        
        ### <b> ### Root
        try:
            ##  <i>  ## identify branch
            self.verbose(0,3,'=> %s' % self.pathSummary(rootpath),1)
            rootperc = 0
            rootdis = 0
            b = -1  # Counter for branch in rootpath
            while rootdis < (maxdis/2): # Find halfway point of extremes = root
                b += 1
                root = rootpath[b]  # Current branch
                rootdis += root.stat['Length'] # Add length until half maxdis exceeded
            ##  <ii>  ## Pick spot
            rootdis -= (maxdis/2)   # Overhang in current branch
            #rootdis = root.stat['Length'] - rootdis
            try:
                rootperc = rootdis / root.stat['Length']   # Distance along branch anc -> desc
            except(ZeroDivisionError):
                self.log.errorLog('Branch %s has zero length?! (%d)' % (root.show(),root.stat['Length']))
            # Work out which end of branch to root from
            firstnode = None
            lastnode = None
            if b > 0:  # Occurred on first branch
                firstnode = root.commonNode(rootpath[b-1])
            if b < (len(rootpath)-1):
                lastnode = root.commonNode(rootpath[b+1])
            if firstnode == root.descNode() or lastnode == root.ancNode():   # Rooted 'backwards'
                rootperc = 1 - rootperc
            #self.log.printLog('#ROOT','Midpoint root found! (%.1f of %s)' % (rootperc,root.show()),2)
            self.placeRoot('Midpoint',root,fraction=rootperc)
            return True
        except:
            self.log.errorLog('Midpoint Rooting Problem!')
            raise
#########################################################################################################################
    def randomRoot(self):   ### Places root on random branch, ignoring branch lengths.
        """Places root on random branch, ignoring branch lengths."""
        try:
            if self.opt['Rooted']:
                self.unRoot()
            rb = random.randint(0, len(self.branch))
            self.verbose(1,3,'Placing random root on branch: %s' % self.branch[rb].show(),1)
            self.placeRoot('Random',self.branch[rb])
            return True
        except:
            self.log.errorLog('Major Problem during random rooting.')
            return True     # Still want to repeat even if failure
#########################################################################################################################
    def randomWtRoot(self):     ### Places root on random branch, weighted by branch lengths.
        """Places root on random branch, weighted by branch lengths."""
        try:
            if self.opt['Rooted']:
                self.unRoot()
            tblen = 0
            for branch in self.branch:
                tblen += branch.stat['Length']
            if tblen == 0:
                return self.randomRoot()   # No branch lengths to use!
            else:
                rlen = random.random() * tblen
                for branch in self.branch:
                    if rlen <= branch.stat['Length']:  # Root
                        try:
                            rf = rlen/branch.stat['Length']
                        except(ZeroDivisionError):
                            self.log.errorLog("Branch %s has zero length?! (%f)" % (branch.show(),branch.stat['Length']))
                            rf = 0.0
                        self.verbose(1,3,'Placing random root at %f of branch: %s' % (rf,branch.show()),1)
                        self.placeRoot('Random (Weighted)',branch,fraction=rf)
                        return True
                    else:
                        rlen -= branch.stat['Length']
            self.log.errorLog("Problem with randomWtRoot - unable to find root!",printerror=False)
            return True     # Still want to repeat even if failure
        except:
            self.log.errorLog('Major Problem during weighted random rooting.')
            return True     # Still want to repeat even if failure
#########################################################################################################################
    def manRoot(self): ### Gives manual options for rooting.
        '''
        Gives manual choices for rooting.
        << False if unchanged, True if changed.
        '''
        ### <a> ### Show Choices
        try:
            print('\n#*# Rooting Options [%s] #*#\n' % self.opt['Rooted'])
            current = self._getRootNode()
            if current != None:
                current = current.info['Name']
            print('Current rooting: %s\n' % current)
            keepok = True
            if self.info['Rooting'] == 'yes': 
                keepok = self.opt['Rooted']
            if keepok:
                print(' <0> Keep current rooting.')
                print(' --- ')
            print(' <1> Midpoint Root.')
            print(' <2> Root from Outgroup File.')
            print(' <3> Make Outgroup File.')
            print(' <4> Root on Branch.')
            print(' <5> Manual rooting using text Tree (with editing options).')
            print(' <6> Random Root.')
            print(' <7> Random Root (branch-length weighted).')
            if self.opt['Rooted']:
                print(' --- \n <8> Unroot.')
        except:
            self.log.errorLog('Problem with manRoot() choice display.')
        ### <b> ### Make Choice
        try:
            choice = rje.getInt(text='\nChoice for Rooting?',blank0=True)
            if (choice == 0) and keepok: # Keep current rooting
                return False
            elif choice == 1: # Midpoint Root
                return self.midRoot()
            elif choice == 2: # Root from file
                try:
                    rootfile = raw_input('Name of Outgroup file: ')
                    TEST = open(rootfile, 'r')
                    TEST.close()
                    return self.fileRoot(rootfile)
                except(IOError):    # If file does not exist then give fresh choice
                    print('%s not found!' % rootfile)
            elif choice == 3: # Make outgroup file. Give option afterwards to root from it!
                return self.makeRootFile()
            elif choice == 4: # Root on Branch.
                return self.branchRoot()  # Options for terminal or internal. <R>oot, <N>ext, <C>ancel
            elif choice == 5: # editTree() rooting
                self.editTree()
                return False
            elif choice == 6: # Random root
                return self.randomRoot()     # Root on random branch
            elif choice == 7: # Random root (weighted)
                return self.randomWtRoot()     # Option for weighting by branch lengths
            elif (choice == 8) and self.opt['Rooted']: # Unroot
                return self.unRoot()
            else:   # Failed to choose
                return self.manRoot()
        except:
            self.log.errorLog('Problem with manRoot() choice.')
            return False
#########################################################################################################################
    def branchRoot(self,nextb=0,inroot=None):    ### Manually places root on branch. Option to save outgroup.
        '''
        Manually places root on branch. Option to save outgroup.
        '''
        try:
            ### <a> ### Setup & Unroot if rooted
            b = nextb
            if inroot == None:
                inroot = self.opt['Rooted']
                if self.opt['Rooted']:
                    self.unRoot()
            branch = self.branch[b]
            clades = self.branchClades(branch)
            ### <b> ### Choices
            print('\nBranch %s:\n%s\n vs\n%s\n' % (branch.show(),self.nodeList(clades[0]),self.nodeList(clades[1])))
            rje.printf('<R>oot, <N>ext, <P>revious, <I>nternal branches, <Q>uit',newline=False)
            choice = rje.choice(': ',default='N').upper()
            print(choice)
            if choice.find('R') == 0:
                self.placeRoot('Manual',branch)
                return True
            elif choice.find('P') == 0:
                b -= 1
            elif choice.find('I') == 0:
                b = self.stat['SeqNum']
            elif choice.find('Q') == 0:
                if inroot:
                    return True
                else:
                    return False
            elif choice[0] == 'N':
                b += 1
            ### <c> ### 'Wraparound' branches and iterate
            if b < 0:
                b = len(self.branch) - 1
            if b >= len(self.branch):
                b = 0
            return self.branchRoot(nextb=b,inroot=inroot)
        except:
            self.log.errorLog('Major Problem with branchRoot()')
            if inroot:
                return True
            else:
                return False
#########################################################################################################################
    def fileRoot(self,rootfile,inroot=False):    ### Places root on unrooted tree based on outgroup sequences
        '''
        Places root on unrooted tree based on outgroup sequences.
        >> rootfile:str = filename
        >> inroot:boolean = whether tree already rooted (for Exception return)
        << returns True if rooting changed, False if not
        '''
        try:
            ### Read sequence names for 'outgroup' ###
            self.verbose(0,3,'Rooting tree using %s...' % rootfile,1)
            rootseqs = self.loadFromFile(rootfile,chomplines=True)
            ### Unroot if rooted ###
            if self.opt['Rooted']:
                inroot = True
                self.unRoot()

            ### Make list of outgroup nodes from sequence names/numbers in file ###
            # - Numbers must match order of input sequences
            # - Names must be included in names of nodes
            outgroup = []
            ## <i> ## Numbered only
            named = False
            for name in rootseqs:
                if re.search('\D.*',name):
                    named = True
                elif name and rje.matchExp('(\d+)',name):
                    num = int(rje.matchExp('(\d+)',name)[0])
                    if num <= self.stat['SeqNum']:
                        outgroup.append(self.node[num-1])
                    else:
                        self.log.errorLog('Sequence numbers from %s exceed number of sequences (%d)! Trying as names...' % (rootfile,self.stat['SeqNum']),printerror=False)
                        named = True
                        break
            ## <ii> ## Named
            if named:
                outgroup = []
                for name in rootseqs:
                    if rje.matchExp('^(\S.+)\s*$',name):
                        match = rje.matchExp('^(\S.+)\s*$',name)[0]
                        for n in range(self.stat['SeqNum']):
                            if self.node[n].info['Name'].find(match) == 0:
                                outgroup.append(self.node[n])

            ### Place Root using outgroup ###
            # - include all outgroup seqs in clade but minimise others
            root = None
            if len(outgroup) == 0:
                self.log.errorLog("No %s Root Sequences found in tree!" % rootfile,printerror=False)
                raise
            elif len(outgroup) == 1:  # Single outgroup sequence
                root = outgroup[0].branch[0]
            else:
                minseq = len(self.node)
                for branch in self.branch:
                    clades = self.branchClades(branch)
                    outbranch = [True,True]
                    for node in outgroup:
                        if node in clades[0]:
                            outbranch[1] = False
                        elif node in clades[1]:
                            outbranch[0] = False
                        else:
                            self.log.errorLog('Major Problem in fileRoot(). Outgroup node missing from Branch Clades.',printerror=False)
                            raise
                    if outbranch[0] and (len(clades[0]) < minseq):
                        minseq = len(clades[0])
                        root = branch
                    elif outbranch[1] and (len(clades[1]) < minseq):
                        minseq = len(clades[1])
                        root = branch

            ### Place Root ###
            if root == None:
                self.log.errorLog("Cannot place root in tree! Check outgroups in single clade." % rootfile,printerror=False)
                raise
            else:
                self.verbose(1,3,"Rooting on %d of %d input seqs %s" % (len(outgroup),len(rootseqs),self.nodeList(outgroup)),1)
                self.placeRoot(method='File (%s)' % rootfile, root=root)
        except:
            self.log.errorLog('Major Problem with fileRoot()')
            if inroot:
                return True
            else:
                return False
#########################################################################################################################
    def makeRootFile(self):     ### Makes outgroup list and places in file. Gives option to root from it.
        '''
        Makes outgroup list and places in file. Gives option to root from it.
        << False if no rooting, True if rooted.
        '''
        try:
            ### <a> ### Make file
            while 1:
                rje.printf('Name of Outgroup file: ',newline=False)
                outfile = raw_input('')
                print(outfile)
                if outfile != '' and rje.yesNo("Save as '%s'" % outfile):
                    break
            OUT = open(outfile,'w')
            ### <b> ### Choose sequence
            outgroup = []
            snum = self.stat['SeqNum']
            s = 0
            while len(outgroup) < (snum-1):  # Must have at least one sequence not in outgroup!
                if self.node[s] in outgroup:
                    s += 1
                else:
                    print('%s: %s' % (rje.preZero((s+1),snum),self.node[s].info['Name']))
                    print(' <O>utgroup, <N>ext, <P>revious, <Q>uit')
                    action = raw_input(': ').lower()
                    if action.find('o') == 0:
                        OUT.write('%s\n' % self.node[s].info['Name'])
                        outgroup.append(self.node[s])
                    elif action.find('p') == 0:
                        s -=1
                    elif action.find('q') == 0:
                        if rje.yesNo('Finish file with %d outgroup seqs %s?' % (len(outgroup),self.nodeList(outgroup))):
                            break
                    else:
                        s +=1
                    if s < 0:
                        s = snum - 1
                    elif s > snum:
                        s = 0
            ### <c> ### Finish
            OUT.close()
            self.log.printLog('#OUT','%d Outgroup seqs saved in %s %s\n' % (len(outgroup),outfile,self.nodeList(outgroup)),1)
            if rje.yesNo('Root using %s' % outfile):
                return self.fileRoot(outfile)
            else:
                return False
        except:
            self.log.errorLog('Major Problem with makeRootFile()')
            return False
#########################################################################################################################
    def placeRoot(self,method,root,fraction=None):     ### Places root on unrooted tree.
        '''
        Places root on unrooted tree.
        >> method:str = Method of rooting
        >> root:Branch = branch on which to place root
        >> fraction:float = fraction along branch to place root anc -> desc
        '''
        try:
            ### <a> ### Setup
            calcfrac = False
            if fraction == None:
                calcfrac = True
            elif (fraction < 0) | (fraction > 1):
                self.log.errorLog("Root was to be placed outside boundaries of root branch!",printerror=False)
                calcfrac = True
            rootbuffer = self.stat['RootBuffer']
            rootlen = root.stat['Length']
            self.verbose(2,3,'\nPlacing root on %s' % root.show(),1)

            ### <b> ### Position root along branch if placement not given
            if calcfrac:
                ## <i> ## Calculate mean distance from tips to middle of branch
                self.verbose(1,3,'Calculating position on branch to place root...',1)
                rdis = [0.0,0.0]    # Mean distance to middle of branch
                rx = [0,0]          # Number of seqs from middle of branch
                rn = root.node
                for s in range(self.stat['SeqNum']):
                    s_path = self.pathLink(self.node[s],rn[0])
                    if root in s_path:  # Path always includes root to give at least one branch
                        n = 1
                    else:
                        n = 0
                        s_path.append(root)
                    rdis[n] = rdis[n] + float(self.pathLen(s_path))
                    rx[n] += 1
                rdis[0] = (rdis[0]/rx[0]) - (rootlen/2)
                rdis[1] = (rdis[1]/rx[1]) - (rootlen/2)
                self.verbose(2,3,'%d Tips through %d = %f\n%d Tips through %d = %f' % (rx[0],rn[0].stat['ID'],rdis[0],rx[1],rn[1].stat['ID'],rdis[1]),2)
                ## <ii> ## Place root nearest the deeper node
                if rdis[0] >= rdis[1]:   # Place Root nearer rn[0]
                    dif = rdis[0] - rdis[1]
                    dis = (rootlen - dif) / 2
                    if dis < rootbuffer:
                        dis = rootlen * rootbuffer
                else:   # Place Root nearer rn[1]
                    dif = rdis[1] - rdis[0]
                    dis = (rootlen + dif) / 2
                    if (rootlen - dis) < rootbuffer:
                        dis = rootlen * (1 - rootbuffer)
                try:
                    fraction = dis / rootlen
                except(ZeroDivisionError):
                    self.log.errorLog("Root Branch has zero length?! (%d)" % rootlen)
                    fraction = 0.5  # Will not matter if zero length!

            ### <c> ### Place root
            self.log.printLog('\r#ROOT','%s Root placed. (%.1f%% of %s, len %f)' % (method,100.0*fraction,root.show(),rootlen),1)
            ## <i> ## New node and branch
            newnode = Node(log=self.log,cmd_list=self.cmd_list)
            newnode.info['Name'] = '%s Root' % method
            newnode.stat['ID'] = self.nodeNum() + 1
            newbranch = Branch(log=self.log,cmd_list=self.cmd_list)
            newbranch.stat['Length'] = rootlen * (1 - fraction)
            newbranch.stat['Bootstrap'] = root.stat['Bootstrap']
            newbranch.node = [root.node[1],newnode]
            newnode.branch = [root,newbranch]
            ## <ii> ## Adjust existing nodes and branches
            ## Node 'beyond' root
            root.node[1].branch.remove(root)
            root.node[1].branch.append(newbranch)
            ## Rooted branch
            root.stat['Length'] = rootlen * fraction
            root.node.remove(root.node[1])
            root.node.append(newnode)
            self.node.append(newnode)
            self.branch.append(newbranch)
                
            self.opt['Rooted'] = True
            self.info['RootMethod'] = method
            self._renumberNodes()
            if self._checkTree() == False:
                raise ValueError
            self.findDuplications(duptext='Root Placed')
        except:
            self.log.errorLog('Major problem with placeRoot()')
            raise
#########################################################################################################################
    ### <6> ### Path Links etc.                                                                                         #
#########################################################################################################################
    def pathLink(self,node1,node2):     ### Returns list of branches linking nodes
        '''
        Returns list of branches linking nodes.
        >> node1:Node Object
        >> node2:Node Object
        >> retry:boolean  = whether to retry if fails
        << list of branches in order node1 -> node2,
        '''
        try:
            ### <a> ### Find Common branch
            path1 = self.rootPath(node1)
            path2 = self.rootPath(node2)
            path2.reverse()
            for branch in path1:
                if branch in path2:
                    path = path1[:path1.index(branch)] + path2[path2.index(branch)+1:]
                    return path
            ### <b> ### Retry if failure
            path = path1 + path2
            return path
        except:
            self.log.errorLog('pathLink() Problem!')
            raise
#########################################################################################################################
    def rootPath(self,node):     ### Returns path to root or 'trichotomy'
        '''
        Returns path to root or 'trichotomy'.
        << list of branch objects
        '''
        try:
            ### <a> ### First Branch
            path = []
            nextnode = node
            while nextnode.stat['ID'] < self.nodeNum():
                prevnode = nextnode
                for branch in nextnode.branch:
                    if branch.link(nextnode).stat['ID'] > nextnode.stat['ID']:
                        path.append(branch)
                        nextnode = branch.link(nextnode)
                        break
                if prevnode == nextnode:    ## No change!
                    self.log.errorLog('Major problem with node numbering.',printerror=False)
                    raise ValueError('Major problem with node numbering for rootPath()')
            return path
        except:
            self.log.errorLog('Major problem with rootPath().')
            raise
#########################################################################################################################
    def pathLen(self,path,stat='Length'):     ### Returns length of path (list of branches)
        '''
        Returns length of path as defined by given stat.
        >> path:list of Branch objects
        >> stat:key for branch.stat dictionary
        << dis:float = total length as float
        '''
        try:
            dis = 0.0
            for branch in path:
                dis += float(branch.stat[stat])
            return dis
        except:
            self.log.errorLog('pathLen Problem!')
            raise            
#########################################################################################################################
    def pathSummary(self,pathin):     ### Returns summary of path as Node numbers X -> Y etc.
        '''
        Returns summary of path as Node numbers X -> Y etc.
        >> path:list of branch objects
        << summary:str = text summary of path
        '''
        try:
            path = pathin[0:]
            if len(path) == 1:
                return path[0].show()
            else:
                firstlink = path[0].commonNode(path[1])
                summary = '%d' % path[0].link(firstlink).stat['ID']
                while len(path) > 1:
                    nextlink = path[0].commonNode(path[1])
                    summary += ' -> %d' % nextlink.stat['ID']
                    path.pop(0)
                summary += ' -> %d' % path[0].link(nextlink).stat['ID']
                return summary
        except:
            self.log.errorLog('pathSummary() Problem!')
            raise            
#########################################################################################################################
    def branchClades(self,branch,internal=False):  ### Return tuple of lists of sequences either side of branch ('anc','desc')
        '''
        Returns lists of sequences either side of branch, (anc,desc)
        >> branch:Branch
        >> internal:bool [False] = whether to return internal nodes as well as termini
        << tuple:lists of Node Objects
        '''
        try:
            ntuple = ([],[])
            for node in self.node:
                if internal or len(node.branch) == 1:   ### Terminal node
                    if node.stat['ID'] > self.seqNum() and not internal:
                        self.log.errorLog('Terminal node %s not really terminal? ID:%d, %d termini?' % (node.info['Name'],node.stat['ID'],self.seqNum()),True,printerror=False)
                        print(self.obj['SeqList'])
                        print(self.node)
                    link = self.pathLink(branch.ancNode(),node)
                    if branch in link: ntuple[1].append(node) ### Node is in 'descendant' clade
                    else: ntuple[0].append(node)
            return ntuple
        except:
            self.log.errorLog('Major Problem with branchClades(%s)' % branch)
            print(branch.stat, branch.node)
            raise
#########################################################################################################################
    def _descClades(self,node,internal=False):     ### Returns tuple of lists of both descendant clades
        '''
        Returns tuple of lists of both descendant clades.
        >> node:Node Object = ancestral node
        >> internal:bool [False] = whether to return internal nodes as well as termini
        << tuple:lists of node objects in descendant clades
        '''
        try:
            if len(node.branch) == 1: return ([node],[]) # Terminal Sequence
            if len(node.branch) == 2: # Root
                desc = node.neighbours()
            else:    
                ancb = node.ancBranch()
                ancn = ancb.link(node)
                desc = node.neighbours(ignore=[ancn])
            return (self.branchClades(node.link(desc[0]),internal=internal)[1],self.branchClades(node.link(desc[1]),internal=internal)[1])
        except:
            self.log.errorLog('Fatal Error with _descClades().')
            raise
#########################################################################################################################
    def _nodeClade(self,node,internal=False):  ### Returns list of descendant clade
        '''Returns list of descendant clade.'''
        if node == self.node[-1] or node.ancBranch() == None:
            if internal: return self.nodes()
            else: return self.node[:self.stat['SeqNum']]
        else: return self.branchClades(node.ancBranch(),internal=internal)[1]
#########################################################################################################################
    def _nodeSeqs(self,nodelist):  ### Returns list of seqs from list of nodes
        '''Returns list of descendant clade.'''
        if self.obj['SeqList']:
            seqs = []
            for node in nodelist: seqs.append(node.obj['Sequence'])
            return seqs
        else: return []
#########################################################################################################################
    def _bestCladeNode(self,nodelist):  ### Returns the best ancestral node for the nodelist clade given
        '''
        Returns the best ancestral node for the nodelist clade given.
        - include all seqs in clade but minimise others.
        >> nodelist:list of node objects
        << cladenode:Node object
        '''
        try:
            cladenode = None
            if len(nodelist) == 0:  # No sequences
                cladnode = None
            if len(nodelist) == 1:  # Single sequence
                cladnode = nodelist[0]
            else:
                minseq = len(self.node)
                for cnode in self.node:
                    clade = self._nodeClade(cnode)
                    ok = True
                    extra = 0
                    for node in self.node[:self.stat['SeqNum']]:
                        if node in nodelist and node not in clade:
                            ok = False
                        elif node in clade and node not in nodelist:
                            extra += 1
                    if ok and extra < minseq:
                        minseq = extra
                        cladenode = cnode
            return cladenode
        except:
            self.log.errorLog('Problem with _bestCladeNode().')
            raise
#########################################################################################################################
    def outGroupNode(self,node):    ### Returns 'outgroup' to node = other descendant of ancestral node
        '''
        Returns 'outgroup' to node = other descendant of ancestral node.
        >> node:Node Object
        '''
        anc = node.ancNode()
        out = anc.neighbours(ignore=[node,anc.ancNode()])
        if len(out) != 1:
            self.log.errorLog('Problem with outGroupNode(). Not returning a single node. (%s)' % out,printerror=False)
            return None
        return out[0]
#########################################################################################################################
    def cladeSpec(self,node):   ### Returns dictionary of {species code:count} for descendant clade (assumes gn_sp__acc)
        '''Returns list of species codes for descendant clade (assumes gn_sp__acc).'''
        specdict = {}
        for node in self._nodeClade(node,internal=False):
            try: spec = rje.split(node.shortName(),'_')[1]
            except: spec = 'Unknown'
            if spec not in specdict: specdict[spec] = 0
            specdict[spec] += 1
        return specdict
#########################################################################################################################
    ### <7> ### Duplications                                                                                            #
#########################################################################################################################
    def _isAnc(self,an,dn):     ### Returns boolean whether an is ancestral to dn
        '''
        Returns boolean whether an is descendant of dn.
        >> an:Node = putative ancestral node 
        >> dn:Node = putative descendant node
        << anc:boolean = whether an is ancestral of dn
        '''
        try:
            if an == dn:
                return False
            if an.stat['ID'] < dn.stat['ID']:
                return False
            #path = self.pathLink(an,dn)
            #for branch in self.node[-1].branch:
            path = self.pathLink(dn,self.node[-1])  # Path to root
            for branch in an.branch:
                if branch in path:  
                    return True
            return False
        except:
            self.log.errorLog('Problem with _isAnc().')
            raise
#########################################################################################################################
    def findDuplications(self,species=None,duptext=''):    ### Finds Duplication Nodes
        '''
        Finds Duplication nodes in tree.
        >> species:str = mark duplication for this species only [None]
        '''
        try:
            ### <a> ### Make sure sequences are present and have species
            if self.obj['SeqList'] == None:
                self.log.printLog('#DUP','Cannot find Duplications without sequence information.')
                for node in self.node: node.opt['Duplication'] = False
                return False
            elif self.opt['Rooted'] == False:
                self.log.printLog('#DUP','Unrooted tree: clearing Duplications.')
                for node in self.node: node.opt['Duplication'] = False
                return
            else:
                for seq in self.obj['SeqList'].seq:
                    if seq.info['Species'] == 'Unknown':
                        seq.info['Species'] = seq.info['SpecCode']
            ### <b> ### Identify Duplications using species data
            if duptext: duptext = ' (%s)' % duptext
            for node in self.node:
                self.log.printLog('\r#DUP','Finding Duplications %s: %.1f%%' % (duptext,100.0*self.node.index(node)/len(self.node)),newline=False,log=False)
                node.opt['Duplication'] = False
                if len(node.branch) > 1:    # Internal
                    clade = self._descClades(node)
                    cladespec = []
                    for s in clade[0] + clade[1]:
                        spcode = s.obj['Sequence'].info['SpecCode']
                        if spcode.startswith('9') and not self.getBool('9SPEC'): spcode = 'UNK'
                        if spcode not in cladespec: cladespec.append(spcode)
                    for s in clade[0]:
                        if node.opt['Duplication']: break
                        elif (species == None) or (species == s.obj['Sequence'].info['Species']) or (species == s.obj['Sequence'].info['SpecCode']):
                            myspec = s.obj['Sequence'].info['SpecCode']
                            if myspec.startswith('9') and not self.getBool('9SPEC'): myspec = 'UNK'
                            for t in clade[1]:
                                tspec = t.obj['Sequence'].info['SpecCode']
                                if tspec.startswith('9') and not self.getBool('9SPEC'): tspec = 'UNK'
                                if node.opt['Duplication']: break
                                elif (tspec != 'UNK') and (tspec == myspec): node.opt['Duplication'] = True
                if node.opt['Duplication'] and len(cladespec) < self.stat['SpecDup'] and 'UNK' not in cladespec:
                    node.opt['SpecDup'] = True; node.opt['Duplication'] = False
            self.log.printLog('\r#DUP','Finding Duplications %s: 100.0%%' % duptext)
        except: self.log.errorLog('Problem in findDuplications().',printerror=True,quitchoice=False)
#########################################################################################################################
    def findDuplicationsNoSeq(self,species=None,duptext=''):    ### Finds Duplication Nodes
        '''
        Finds Duplication nodes in tree.
        >> species:str = mark duplication for this species only [None]
        '''
        try:
            ### <a> ### Make sure sequences are present and have species
            if self.opt['Rooted'] == False:
                self.log.printLog('#DUP','Unrooted tree: clearing Duplications.')
                for node in self.node: node.opt['Duplication'] = False
                return
            ### <b> ### Identify Duplications using species data
            if duptext: duptext = ' (%s)' % duptext
            for node in self.node:
                self.log.printLog('\r#DUP','Finding Duplications %s: %.1f%%' % (duptext,100.0*self.node.index(node)/len(self.node)),newline=False,log=False)
                node.opt['Duplication'] = False
                if len(node.branch) > 1:    # Internal
                    clade = self._descClades(node)
                    cladespec = []
                    for s in clade[0] + clade[1]:
                        spcode = rje.split(s.info['Name'],'_')[1]
                        if spcode.startswith('9') and not self.getBool('9SPEC'): spcode = 'UNK'
                        if spcode not in cladespec: cladespec.append(spcode)
                    for s in clade[0]:
                        spcode = rje.split(s.info['Name'],'_')[1]
                        if spcode.startswith('9') and not self.getBool('9SPEC'): spcode = 'UNK'
                        if node.opt['Duplication']: break
                        elif (species == None) or (species == spcode):
                            myspec = spcode
                            for t in clade[1]:
                                tspec = rje.split(t.info['Name'],'_')[1]
                                if tspec.startswith('9') and not self.getBool('9SPEC'): tspec = 'UNK'
                                if node.opt['Duplication']: break
                                elif (tspec != 'UNK') and (tspec == myspec): node.opt['Duplication'] = True
                if node.opt['Duplication'] and len(cladespec) < self.stat['SpecDup'] and 'UNK' not in cladespec:
                    node.opt['SpecDup'] = True; node.opt['Duplication'] = False
            self.log.printLog('\r#DUP','Finding Duplications %s: 100.0%%' % duptext)
        except: self.log.errorLog('Problem in findDuplicationsNoSeq().',printerror=True,quitchoice=False)
#########################################################################################################################
    ### <8> ### Tree Subfamilies and Groupings                                                                          #
    ### Note that the contents of these methods have been moved to rje_tree_group.py                                    #
#########################################################################################################################
    ## <8A> ## Master Subfamily/Grouping Methods                                                                        #
#########################################################################################################################
    def _checkGroupNames(self):    ### Automatically names groups after given gene unless already renamed
        return rje_tree_group._checkGroupNames(self)
#########################################################################################################################
    def treeGroup(self,callmenu=False):     ### Master Tree Grouping loop.
        '''See rje_tree_group.treeGroup().'''
        rje_tree_group.treeGroup(self,callmenu)
#########################################################################################################################
    def _autoGroups(self,method='man'):     #### Automatically goes into grouping routines
        '''See rje_tree_group._autoGroups().'''
        rje_tree_group._autoGroups(self,method)
#########################################################################################################################
    def _checkGroups(self):     ### Checks that Group selection does not break 'rules'
        '''See rje_tree_group._checkGroups().'''
        return rje_tree_group._checkGroups(self)
#########################################################################################################################
    ## <8B> ## Grouping Gubbins                                                                                         #
#########################################################################################################################
    def _orphanCount(self):     ### Returns number of orphan sequences
        '''See rje_tree_group._orphanCount().'''
        return rje_tree_group._orphanCount(self)
#########################################################################################################################
    def _purgeOrphans(self): ### Removes orphan nodes
        '''See rje_tree_group._purgeOrphans().'''
        return rje_tree_group._purgeOrphans(self)
#########################################################################################################################
    def _resetGroups(self):     ### Clears compression of non-group nodes
        '''See rje_tree_group._resetGroups().'''
        rje_tree_group._resetGroups(self)
#########################################################################################################################
    def _addGroup(self,node):   ### Adds group based at node, removing existing descendant groups
        '''See rje_tree_group._addGroup().'''
        rje_tree_group._addGroup(self,node)
#########################################################################################################################
    def _grpVarDel(self,delseq=[],kept=None,fam=None):     ### Deletes all variants in list
        '''See rje_tree_group._grpVarDel().'''
        rje_tree_group._grpVarDel(self,delseq,kept,fam)
#########################################################################################################################
    def _grpSeqSort(self,seqs=[],compseq=None):  ### Reorders seqs according to %ID, Gaps and Extra   
        '''See rje_tree_group._grpSeqSort().'''
        rje_tree_group._grpSeqSort(self,seqs,compseq)
#########################################################################################################################
    def _reorderSeqToGroups(self):  ### Reorders sequences according to groups #!# Add query? And _grpSeqSort()?
        '''See rje_tree_group._reorderSeqToGroups().'''
        rje_tree_group._reorderSeqToGroups(self)
#########################################################################################################################
    ## <8C> ## Interactive Grouping                                                                                     #
#########################################################################################################################
    def _groupChoice(self):     ### Gives manual options for grouping and updates self.info['Grouping'].
        '''See rje_tree_group._groupChoice().'''
        return rje_tree_group._groupChoice(self)
#########################################################################################################################
    def _sumGroups(self):   ### Prints summary of Groups
        '''See rje_tree_group._sumGroups().'''
        return rje_tree_group._sumGroups(self)
#########################################################################################################################
    def _groupRules(self):     ### Options to change grouping options
        '''See rje_tree_group._groupRules().'''
        rje_tree_group._groupRules(self)
#########################################################################################################################
    def _saveGroups(self,filename='rje_tree.grp',groupnames=True):   # Saves sequence names in Groups
        '''See rje_tree_group._saveGroups().'''
        rje_tree_group._saveGroups(self,filename,groupnames)
#########################################################################################################################
    ## <8D> ## Grouping Methods                                                                                         #
#########################################################################################################################
    def _clearGroups(self):     ### Clears current group selection
        '''See rje_tree_group._clearGroups().'''
        return rje_tree_group._clearGroups(self)
#########################################################################################################################
    def _loadGroups(self,filename='rje_tree.grp'):   # Saves sequence names in Groups
        '''See rje_tree_group._loadGroups().'''
        rje_tree_group._loadGroups(self,filename)
#########################################################################################################################
    def _dupGroup(self):     ### Duplication grouping
        '''See rje_tree_group._dupGroup().'''
        rje_tree_group._dupGroup(self)
#########################################################################################################################
    ## <8E> ## Review Grouping                                                                                          #
#########################################################################################################################
    def _reviewGroups(self,interactive=1):    ### Summarise, scan for variants (same species), edit group
        '''See rje_tree_group._reviewGroups().'''
        rje_tree_group._reviewGroups(self,interactive)
#########################################################################################################################
    def _groupDisSum(self,seq1,seq2,text=''):    ### Prints ID, Gaps and Extra Summary of seq1 vs seq2
        '''See rje_tree_group._groupDisSum().'''
        rje_tree_group._groupDisSum(self,seq1,seq2,text)
#########################################################################################################################
    ### <9> ### Miscellaneous Tree Methods                                                                              #
#########################################################################################################################
    def branchPam(self,pam=None):   ### Calculate PAM distances for branches
        '''
        Calculates PAM distances for branches.
        >> pam:rje_pam.PamCtrl Object [None]
        '''
        try:### ~ [1] ~ Set up PAM Control Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if pam: self.obj['PAM'] = pam
            if not self.obj['PAM']: self.obj['PAM'] = rje_pam.PamCtrl(self,self.log,self.cmd_list)
            ### ~ [2] ~ Calculate PAM distances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (bx,btot) = (0.0,self.branchNum())
            for branch in self.branch:
                self.printLog('\r#PAM','Calculating branch ML PAM: %.1f%%' % (bx/btot),newline=False,log=False)
                bx += 100.0
                n = [branch.node[0].stat['ID'],branch.node[1].stat['ID']]
                if branch.node[0].obj['Sequence'] == None or branch.node[1].obj['Sequence'] == None:
                    self.log.errorLog('Attempting to calculate PAM distance with missing sequence (%d or %d)' % (n[0],n[1]),printerror=False)
                else:
                    seq = [branch.node[0].obj['Sequence'].info['Sequence'],branch.node[1].obj['Sequence'].info['Sequence']]
                    if n[0] > n[1]: branch.stat['PAM'] = self.obj['PAM'].pamML('Branch %d->%d' % (n[0],n[1]),seq[0],seq[1])
                    else: branch.stat['PAM'] = self.obj['PAM'].pamML('Branch %d->%d' % (n[1],n[0]),seq[1],seq[0])
                self.verbose(2,4,'Branch %s = PAM%d' % (branch.show(),branch.stat['PAM']),1)
            self.printLog('\r#PAM','Calculation of branch ML PAM complete.',log=False)
        except:
            self.log.errorLog("Major problem calculating branch PAM distances")
            raise
#########################################################################################################################
    def indelTree(self,filename='indel.txt'):  ### Produces a text tree with indels marked
        '''
        Produces a text tree with indels marked.
        > filename:str = output file name
        '''
        try:
            ### <a> ### Make branch indel dictionaries
            indel = {}
            indeltxt = {}
            keys = ['ins_num','ins_aa','del_num','del_aa']
            self.verbose(0,3,'Calculating indel stats...',0)
            for branch in self.branch:
            ## <i> ## Setup
                indel[branch] = {}
                for key in keys:
                    indel[branch][key] = 0
                status = 'None'
                anc = branch.ancNode().obj['Sequence']
                desc = branch.descNode().obj['Sequence']
                if anc == None or desc == None:
                    self.log.errorLog('Ancestral or Descendant sequences missing!',printerror=False)
                    raise ValueError
                seqlen = anc.seqLen()
                if desc.seqLen() != seqlen:
                    self.log.errorLog('Ancestral (%s) and Descendant (%s) Sequences have different lengths!' % (anc.info['ID'],desc.info['ID']),printerror=False)
                    raise ValueError
            ## <ii> ## Calculate
                for r in range(seqlen):
                    if anc.info['Sequence'][r] == '-' and desc.info['Sequence'][r] == '-':  # Gap
                        status = 'None'
                    elif anc.info['Sequence'][r] != '-' and desc.info['Sequence'][r] != '-': # Seq
                        status = 'None'
                    elif anc.info['Sequence'][r] == '-' and desc.info['Sequence'][r] != '-': # Insertion
                        indel[branch]['ins_aa'] += 1
                        if status != 'Ins':    # New insertion
                            indel[branch]['ins_num'] += 1
                        status = 'Ins'
                    elif anc.info['Sequence'][r] != '-' and desc.info['Sequence'][r] == '-': # Deletion
                        indel[branch]['del_aa'] += 1
                        if status != 'Del':    # New insertion
                            indel[branch]['del_num'] += 1
                        status = 'Del'
            ## <iii> ## Make indeltxt
                indeltxt[branch] = '{+%d(%d),-%d(%d)}' % (indel[branch]['ins_aa'],indel[branch]['ins_num'],indel[branch]['del_aa'],indel[branch]['del_num'])
                        
            self.verbose(0,3,'Done!',1)
        except:
            self.log.errorLog('Major problem with indelTree() Calculations.')
            raise
            
            ### <b> ### Tree Output #!# Copies textTree for now!
        try:
            ## <i> ## Setup
            seqnum = True
            seqname = 'short'
            maxnamelen = 50
            nodename = 'num'
            showboot = True
            showlen = 'branch'
            blen = 'Length'
            scale = 4
            spacer = 1
            pause = 50
            fromnode = None

            ## <ii> ## textTree() #!# Confusing numbering of comments ahead! :o)
            ### <a> ### Setup
            ignore = []
            ## <i> ## Scale
            xpos = {}
            if scale <= 0:
                self.log.errorLog("Scale must be > 0! Changed to 1.",printerror=False)
                scale = 1
            ## <ii> ## BranchLengths
            if (blen == 'pam') & (self.branch[0].stat['PAM'] < 0):
                self.branchPam()
            ## <iii> ### Base of Tree
            tric = False
            if fromnode == None or fromnode == self.node[-1]:
                basenode = self.node[-1]
                if len(basenode.branch) == 3:  # Trichotomy
                    print('Trichotomy')
                    tric = True
                xpos[basenode] = scale + 1
            else:
                basenode = fromnode
                anc = basenode.ancNode()
                ignore.append(anc)
                xpos[anc] = 1
                branch = basenode.ancBranch()
                if blen == 'fix':
                    mylen = self.stat['DefLen']
                elif blen.lower() == 'pam':
                    mylen = float(branch.stat['PAM']) / 100
                else:
                    mylen = branch.stat[blen]
                xpos[basenode] = int(mylen * scale / self.stat['DefLen']) + 2 + xpos[anc]
            ynode = [basenode]
            ### <b> ### Make a list of order of node output and calculate X-pos
            replacements = True
            while replacements:
                replacements = False
                newy = []
                for node in ynode:          # ! # Add in future 'clade' option
                    if node in ignore or node.stat['ID'] <= self.stat['SeqNum']:  
                        newy.append(node)    # Already replaced/terminus - already handled 
                    elif node.opt['Compress']:  # Compress clade to node
                        replacements = True     # Keep Looping
                        newy.append(node)       # Do not expand
                        ignore.append(node)     
                    elif node.stat['ID'] > self.stat['SeqNum']:
                        replacements = True     # Keep Looping
                        ## <i> ## ypos
                        nodes = node.neighbours(ignore=ignore)
                        newy += [nodes[0],node] + nodes[1:]
                        ignore.append(node)
                        ## <ii> ## xpos
                        for i in range(len(nodes)):
                            branch = nodes[i].link(node)
                            if blen == 'fix':
                                mylen = self.stat['DefLen']
                            elif blen.lower() == 'pam':
                                mylen = float(branch.stat['PAM']) / 100
                            else:
                                mylen = branch.stat[blen]
#!#
                            mylen = int(mylen * scale / self.stat['DefLen'])
                            mylen = (2 * mylen) + len(indeltxt[branch])  #!# Sandwich indel data!                           
                            xpos[nodes[i]] = mylen + 2 + xpos[node]

                    else:   
                        raise ValueError
                ynode = newy
            #print self.nodeList(ynode)
            ### <c> ### Draw Tree            
            maxpos = 0
            vline = {}
            for node in ynode:
                #print node.stat['ID'], xpos[node]
                vline[node] = False    # List of on/off values for vertical line
                maxpos = max(xpos[node],maxpos)
            self.verbose(0,5,'\n',0)
            treelines=[]
            for node in ynode:
                #print node.stat['ID']
                branch = node.ancBranch()
                ## <i> ## Branch
                if branch == None:  # Root/Trichotomy
                    line = '+'                       # Space before Branch
                    while len(line) < (xpos[node] - 1): line += '-'
                else:
                    anc = branch.link(node)
                    line = ' ' * (xpos[anc] - 1)    # Space before Branch

#!#                    # Start of branch
                    line += '+'
                    if blen == 'fix': mylen = self.stat['DefLen']
                    elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                    else: mylen = branch.stat[blen]
                    mylen = int(mylen * scale / self.stat['DefLen'])
                    line = line + ('-' * mylen) + indeltxt[branch] + ('-' * mylen)
                line += '-+'

                ## <ii> ## Name
                txtname = ''
                if node.stat['ID'] <= self.stat['SeqNum']:      # Terminal node
                    if seqnum or (seqname == 'None') or (seqname == 'num'):
                        txtname = '%d' % node.stat['ID']
                        if seqname != 'num' and seqname != 'None': txtname += ': '
                    if seqname == 'short': txtname += node.shortName()
                    elif seqname == 'long': txtname += node.info['Name']
                    if len(txtname) > maxnamelen: txtname = txtname[:maxnamelen] + '...'                
                else:   # Internal node
                    if nodename == 'none': txtname = ''
                    elif node.opt['Compress']:
                        if nodename == 'num': txtname = '*%d*: ' % node.stat['ID']
                        txtname += node.info['CladeName']
                    elif nodename == 'num': txtname = '%d' % node.stat['ID']
                    elif nodename == 'short': txtname += node.shortName()
                    elif nodename == 'long': txtname += node.info['Name']
                if showboot and self.opt['Bootstrapped'] and branch != None and node.stat['ID'] > self.stat['SeqNum']:
                    txtname += ' [%d]' % branch.stat['Bootstrap']
                # display branchlengths?
                if branch != None:
                    if showlen == 'branch':
                        if blen == 'fix': mylen = self.stat['DefLen']
                        elif blen.lower() == 'pam': mylen = float(branch.stat['PAM']) / 100
                        else: mylen = branch.stat[blen]
                        txtname += ' {%f}' % mylen
                    elif showlen == 'sum': txtname += ' {%f}' % self.pathLen(self.rootPath(node))
                    txtname += ' (%d aa)' % branch.descNode().obj['Sequence'].aaLen()
                else: txtname += ' (%d aa)' % self.node[-1].obj['Sequence'].aaLen()
                line += ' ' + txtname 
                ## <iii> ## vlines anc of xn
                for vnode in ynode:
                    if vline[vnode] and vnode != anc and vnode != node:
                        vx = xpos[vnode] - 1
                        line = line[:vx] + '|' + line[(vx+1):]
                if branch != None and node != fromnode:
                    vline[anc] = not vline[anc]
                    if anc == basenode and tric and vline[anc] == False:
                        vline[anc] = True
                        tric = False
                #self.verbose(0,5,line,1)
                treelines.append('%s\n' % line)
                ## <iv> ## spacer
                for y in range(spacer):
                    line = ' ' * xpos[node]
                    for vnode in ynode:
                        if vline[vnode]:
                            vx = xpos[vnode] - 1
                            line = line[:vx] + '|' + line[(vx+1):]
                #self.verbose(0,5,line,1)
                treelines.append('%s\n' % line)
            ### <c> ### Save/Print Tree
            if filename != None:
                try:
                    open(filename, 'w').write(rje.join(treelines,'\n'))
                    self.log.printLog("#OUT","Text tree saved to %s" % filename,1)
                except:
                    self.log.errorLog("Problem saving text tree to %s" % filename)
                    raise
            else:
                p = 1
                self.verbose(0,5,'\n',0)
                for line in treelines:
                    if p < pause:
                        self.verbose(0,5,line,0)
                        p += 1
                    else:
                        self.verbose(0,1,line,0)
                        p = 0
        except: self.log.errorLog('Major Problem with textTree().')
#########################################################################################################################
## End of Tree Class                                                                                                    #
#########################################################################################################################
    

                                                    ### ~ ### ~ ###

#########################################################################################################################
##  Node Class: Individual Tree Nodes                                                                                   #
#########################################################################################################################
class Node(rje.RJE_Object): 
    '''
    Individual nodes (internal and leaves) for Tree object. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of Node
    - CladeName = Name to be used if describing descendant clade
    - Type = Type of node = Internal/Terminal/Root/(Duplication)
    
    Opt:boolean
    - Duplication = Whether the node is a duplication node
    - Compress = Whether to compress clade in Tree.textree()
    - SpecDup = Wheter node is a species-specific(ish) duplication

    Stat:numeric
    - ID = node number

    Obj:RJE_Objects
    - Sequence = rje_seq.Sequence object

    Other:
    - branch = list of branches (1 for terminal, 2 for root, 3 for internal)
    '''
    ### Old Node Attributes
    ## Info:
    # Name : name=None           # Short name of node
    # - longname=None       # Full name of node
    # CladeName : cladename=None      # Name of node to be used if describing descendant clade
    # Type: type=None           # Type of node = Internal/Terminal/Root/Trichotomy/Duplication
    ## Stat:
    # ID : id=0                # Node ID number
    ## Obj:
    # Sequence : sequence=None       # Node sequence
    ## None:
    # - nodelink=[-1,-1,-1] # Nodes linked to self: [0]&[1]=desc, [2]=anc
    # - nodepath=None       # Path from Node to root

    ### Other attributes
    branch = []   # List of branches linking node to other nodes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','CladeName']
        self.statlist = ['ID']
        self.optlist = ['Duplication','Compress','SpecDup']
        self.objlist = ['Sequence']
        self.listlist = []
        self.dictlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.info['Name'] = 'Node'
        self.branch = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try: self._generalCmd(cmd)
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        return
#########################################################################################################################
    ### <2> ### General Class Methods
#########################################################################################################################
    def setType(self):     ### Sets the node's type as a string
        '''
        Sets the node's type based on link and duplication.
        '''
        try:
            try: type = ['Unformed','Terminal','Root','Internal'][len(self.branch)]
            except: type = 'Erroneous'
            if self.opt['Duplication']: type += ' (Duplication)'
            elif self.opt['SpecDup']: type += ' (Lineage-specific duplication)'
            self.info['Type'] = type
            return type
        except:
            self.log.errorLog('Major problem establishing %s (%d) type' % (self.info['Name'],self.stat['ID']))
            raise
#########################################################################################################################
    def mapSeq(self,seq=None,id=0):  ### Maps a Sequence object onto the node
        '''
        Maps a rje_seq.Sequence object onto Node.
        >> seq:rje_seq.Sequence object
        >> id:int = order of sequence in SeqList (for output clarity only)
        '''
        try:
            if seq == None: return
            self.obj['Sequence'] = seq
            self.stat['ID'] = id
            self.info['CladeName'] = self.info['Name'] = seq.info['Name']
        except:
            self.log.errorLog('Major problem mapping Sequence %d to Node' % id)
            raise
#########################################################################################################################
    def shortName(self):    ### Returns short name.
        '''Returns short name = first word of name.'''
        try:
            word = self.info['Name'].find(' ')
            if self.info['Type'] == 'Terminal' and word > 0: return self.info['Name'][:word]
            else: return self.info['Name']
        except:
            self.log.errorLog('Major problem with shortName(%s)' % self.info['Name'])
            raise
#########################################################################################################################
    def spCode(self):   ### Return species code
        try: return self.obj['Sequence'].spCode()
        except: return 'UNK'
#########################################################################################################################
    def rename(self,rooting=None):   ### Gives node a good name
        '''
        Gives internal nodes a good name based on numbering.
        >> rooting:str = method of rooting if to be added.
        '''
        try:
            if len(self.branch) == 1: return  # Terminus
            links = []
            for branch in self.branch: links.append('%d' % branch.link(self).stat['ID'])
            newname = 'Node %d (%s)' % (self.stat['ID'],rje.join(links,','))
            if self.opt['Duplication']: newname = 'Duplication ' + newname
            elif self.opt['SpecDup']: newname = 'Lineage-specific duplication' + newname
            if len(self.branch) == 2: # Root
                newname = 'Root ' + newname
                if rooting != None: newname = '%s ' % rooting + newname
            self.info['Name'] = newname
        except:
            self.log.errorLog('Major problem with rename(%s)' % self.info['Name'])
            raise
#########################################################################################################################
    def link(self,othernode):   ### Returns branch that links node with self or None if none
        '''
        Returns branch that links node with self or None if none.
        >> othernode:Node Object to link
        << link:Branch Object
        '''
        try:
            for branch in self.branch:
                if othernode in branch.node: return branch
            return None
        except:
            self.log.errorLog('Major Problem with Node link().')
            raise
#########################################################################################################################
    def ancBranch(self):   ### Returns branch that links node with 'ancestor' or None if none
        '''
        Returns branch that links node with 'ancestor' or None if none.
        << branch:Branch Object
        '''
        try:
            for branch in self.branch:
                if branch.link(self).stat['ID'] > self.stat['ID']: return branch
            #self.verbose(1,4,'No ancbranch for node %d - check number of nodes!' % (self.stat['ID']),1)
            return None
        except:
            self.log.errorLog('Major Problem with Node.ancBranch().')
            print(self.branch)
            print(self.stat)
            raise
#########################################################################################################################
    def ancNode(self):   ### Returns linked node that is 'ancestral' or None if none
        '''
        Returns linked node that is 'ancestral' or None if none
        << node:Node Object
        '''
        try:
            for branch in self.branch:
                node = branch.link(self)
                if node.stat['ID'] > self.stat['ID']: return node
            return None
        except:
            self.log.errorLog('Major Problem with Node.ancNode().')
            raise
#########################################################################################################################
    def neighbours(self,ignore=[]):   ### Returns list of Node objects linked by a single branch
        '''
        Returns list of Node objects linked by a single branch.
        >> ignore:list of Node objects
        << neighbours:list of Node Objects
        '''       
        try:
            neighbours = []
            for branch in self.branch:
                node = branch.link(self)
                if node not in ignore:
                    neighbours.append(node)
            return neighbours
        except:
            self.log.errorLog('Major Problem with Node link().')
            raise
#########################################################################################################################
    def flipDesc(self):     ### Flips 'Descendant' Nodes
        '''Flips 'Descendant' Nodes by swapping Branches in self.branch.'''
        try:
            if len(self.branch) < 2: return
            elif len(self.branch) > 2:
                anc = self.ancBranch()
                newb = [anc]
                self.branch.remove(anc)
                newb = newb + [self.branch[1],self.branch[0]]
                self.branch = newb
            else: self.branch = [self.branch[1],self.branch[0]]
        except:
            self.log.errorLog('Major Problem with Node flipDesc().')
            raise            
#########################################################################################################################
## End of Node Class                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
##  Branch Class: Individual Tree Branches                                                                              #
#########################################################################################################################
class Branch(rje.RJE_Object):
    '''
    Individual branches for Tree object. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of Branch
    
    Stat:numeric
    - Bootstrap = Bootstrap Support (-1 = none)
    - Length = Branch Length (-1 = none)
    - PAM = Branch PAM (-1 = none)
    
    Obj:RJE_Objects
    - Sequence = rje_seq.Sequence object

    Other:
    - node = list of two nodes connected to branch
    '''
    ### Old Branch Attributes
    ## Stat:
    # Bootstrap : bootstrap = 0 # Bootstrap Support
    # Length : length = 0    # Branch Length
    # PAM : pam = -1        # Branch PAM (-1 = none)
    ## None: 
    # ancnode = -1  # Ancestral node (parent tree object)
    # descnode = -1 # Descendant node (parent tree object)
    ### Other attributes
    node = []   
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type']
        self.statlist = ['Bootstrap','Length','PAM']
        self.optlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=-1.0,obj=None,setlist=True,setdict=True)
        self.info['Type'] = 'Branch'
        ### Other Attributes ###
        self.node = []
#########################################################################################################################
    ### <2> ### General Class Methods                                                                                   #
#########################################################################################################################
    def link(self,node):    ### Returns other end of branch (Node Object)
        '''
        Returns other end of branch.
        >> node:Node object
        << link:Node object
        '''
        try:
            if node in self.node:   # OK
                if len(self.node) < 2 or len(self.node) > 2:
                    self.log.errorLog('Branch has wrong number of nodes!',printerror=False)
                    raise
                else:
                    for link in self.node:
                        if link != node: return link
            else:
                self.log.errorLog('link() called for wrong branch for given node (%s vs %s)!\n' % (node,self.node),printerror=False)
                raise
        except:
            self.log.errorLog('Major problem with link().')
            raise
#########################################################################################################################
    def combine(self,node,branch):    ### Combines data from another branch (during unrooting)
        '''
        Combines data from another branch, usually during unrooting.
        Will warn if bootstraps are not compatible and die if node is incorrect.
        >> node:Node object = common node between branches
        >> branch:Branch object = other branch
        '''
        try:### ~ [1] ~ Check Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if node in self.node and node in branch.node:   # OK
                if self.stat['Bootstrap'] != branch.stat['Bootstrap']:
                    self.log.printLog('#TREE','%s:[%s] vs %s:[%s]' % (self.info['Name'],self.stat['Bootstrap'],branch.info['Name'],branch.stat['Bootstrap']))
                    self.log.printLog('#TREE','Bootstrap disagreement during branch combining. Will retain higher.')
                    if branch.stat['Bootstrap'] == None: self.log.errorLog('Bootstrap Missing! Retaining non-missing value.',printerror=False)
                    elif self.stat['Bootstrap'] < branch.stat['Bootstrap']: self.stat['Bootstrap'] = branch.stat['Bootstrap']
            ### ~ [2] ~ Combine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.info['Name'] = 'Combined %s and %s' % (self.info['Name'],branch.info['Name'])
                self.stat['Length'] += branch.stat['Length']
                if self.stat['PAM'] >= 0 and branch.stat['PAM'] >= 0: self.stat['PAM'] += branch.stat['PAM']
                else: self.stat['PAM'] = -1
                self.node.remove(node)
                branch.node.remove(node)
                self.node += branch.node
            ### ~ [3] ~ Big Problem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                self.log.errorLog('Incorrect common node given for combining branches!',printerror=False)
                raise
        except:                
            self.log.errorLog('Major problem in combine().')
            raise
#########################################################################################################################
    def show(self):     ### Shows nodes X -> Y
        '''Returns Node numbers X -> Y.'''
        if self.node[0].stat['ID'] > self.node[1].stat['ID']: return '%d -> %d' % (self.node[0].stat['ID'],self.node[1].stat['ID'])
        else: return '%d -> %d' % (self.node[1].stat['ID'],self.node[0].stat['ID'])
#########################################################################################################################
    def ancNode(self):  ### Returns 'ancestral' node
        '''Returns 'ancestral' node.'''
        try:
            if self.node[0].stat['ID'] > self.node[1].stat['ID']: return self.node[0]
            else: return self.node[1]
        except:
            self.log.errorLog('Problem with ancNode().')
            raise
#########################################################################################################################
    def descNode(self):  ### Returns 'descendant' node
        '''Returns 'descendant' node.'''
        try:
            if self.node[0].stat['ID'] < self.node[1].stat['ID']: return self.node[0]
            else: return self.node[1]
        except:
            self.log.errorLog('Problem with descNode().')
            raise
#########################################################################################################################
    def commonNode(self,branch):    ### Common node with another branch
        '''
        Returns common node with other branch.
        >> branch:Branch Object
        << node:Node Object or None if not common.
        '''
        try:
            if self.node[0] in branch.node: return self.node[0]
            elif self.node[1] in branch.node: return self.node[1]
            else: return None
        except:
            self.log.errorLog('Major problem finding Common Node (%s vs %s)' % (self.node,branch.node))
            return None
#########################################################################################################################
## End of Branch Class                                                                                                  #
#########################################################################################################################

#########################################################################################################################
## End of SECTION II: Module Classes                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def treeMenu(out,mainlog,cmd_list,tree=None):        ### Menu for handling Tree Functions in 'standalone' running
    '''
    Menu for handling Tree Functions in 'standalone' running.
    - load sequences
    - load tree
    - save tree
    - edit tree
    - define subfams
    '''
    try:
        ### <a> ### Setup parameters etc.
        if tree == None:
            print('No Tree!')
            return
            #tree = Tree(log=mainlog,cmd_list=cmd_list)
        ### <b> ### Menu
        while out.stat['Interactive'] >= 0:
            ## <i> ## Options
            print('\n\n *** Tree Menu *** \n')
            print('<L>oad Tree')
            print('<M>ake Tree')
            if len(tree.node) > 0:
                print('<S>ave Tree')
                print('<I>mport Sequence Data')
                print(' --- \n<R>oot Options')
                print('<E>dit Tree')
                print('<G>rouping Options')
            if len(tree.node) > 0 and tree.obj['SeqList']:
                print(' --- \n<A>ncestral Sequence Prediction (GASP)')
                print('E<x>port Sequence Data')
            print(' --- \n<Q>uit')
            ## <ii> ## Choice
            choice = rje.choice('\nChoice: ',default='Q').upper()
            if choice.find('L') == 0:   # Load Tree
                filename = rje.choice('Filename: ')
                useseq = False
                if tree.obj['SeqList']: useseq = rje.yesNo('Use currently loaded sequence list?')
                try:
                    if useseq: tree.loadTree(file=filename,seqlist=tree.obj['SeqList'])
                    else: tree.loadTree(file=filename)
                except: continue
            elif choice.find('M') == 0:   # Load Tree
                tree.makeTreeMenu(interactiveformenu=0,force=False,make_seq=tree.obj['SeqList'])
            elif choice.find('S') == 0 and len(tree.node) > 0: # Save
                filename = rje.choice('Filename: ',default=tree.info['SaveTree'])
                if tree.opt['OutputBranchLen']: withbranchlengths = 'Length'
                else: withbranchlengths = 'none'
                outnames = tree.info['OutNames']
                maxnamelen = tree.stat['TruncNames']
                useseqnum = tree.opt['SeqNum'] and tree.obj['SeqList']
                if not rje.yesNo('Save with branch lengths?'): withbranchlengths = 'none'
                if rje.yesNo('Use full sequence names?'):
                    outnames = 'long'
                    if not rje.yesNo('Truncate long names for program compatability?'):
                        maxnamelen = 0
                if tree.obj['SeqList']: useseqnum = rje.yesNo('Save with sequence numbers?',default='N')
                #!# Add ...  ,type=self.info['SaveType'] #!#
                tree.saveTree(filename=filename,seqnum=useseqnum,seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            elif choice.find('I') == 0 and len(tree.node) > 0: # Import Seqs
                filename = rje.choice('Filename: ')
                seqs = rje_seq.SeqList(log=mainlog,cmd_list=cmd_list+['accnr=F'])
                seqs.loadSeqs(seqfile=filename,nodup=False)
                tree.mapSeq(seqlist=seqs)
                tree.textTree()
            elif choice.find('X') == 0 and tree.obj['SeqList']: # Import Seqs
                filename = rje.choice('Filename: ')
                tree.obj['SeqList'].saveFasta(seqfile=filename)
            elif choice.find('R') == 0 and len(tree.node) > 0: # Root Options
                tree.info['Rooting'] = 'man'
                tree.treeRoot()
            elif choice.find('G') == 0 and len(tree.node) > 0: # Grouping Options
                if tree.stat['MinFamNum'] < 1:
                    tree.stat['MinFamNum'] = 1
                tree.treeGroup(callmenu=True)
            elif choice.find('E') == 0 and len(tree.node) > 0: # Edit Options
                tree.editTree(reroot=True)
            elif choice.find('A') == 0 and len(tree.node) > 0 and tree.obj['SeqList'] != None: # GASP
                filename = rje.choice('Root filename (FILE.anc.fas, FILE.anc.nsf, FILE.txt): ',default=tree.obj['SeqList'].info['Basefile'])
                if rje.yesNo('Use %s as root filename?' % filename):
                    mygasp = rje_ancseq.Gasp(tree=tree,ancfile=filename,cmd_list=cmd_list,log=mainlog)
                    out.verbose(0,2,'%s' % mygasp.details(),1)
                    if not rje.yesNo('Use these parameters?'): mygasp.edit()
                    mygasp.gasp()
            elif choice.find('Q') == 0 and rje.yesNo('Quit Tree Menu?'): # Quit
                return tree
        return tree
    except IOError:
        out.verbose(0,1,'Problem with given File!',1)
        treeMenu(out,mainlog,cmd_list,tree)
    except KeyboardInterrupt:
        if rje.yesNo('Quit Tree Menu?'): raise 
        else: treeMenu(out,mainlog,cmd_list,tree)
    except:
        mainlog.errorLog('Fatal Error in main TreeMenu()')
        raise    
#########################################################################################################################
### END OF SECTION III                                                                                                  #
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
        print('Unexpected error during program setup:', sys.exc_info()[0])
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        Tree(log=mainlog,cmd_list=cmd_list).run()
        #if 'reroot' in cmd_list: tree.saveTree(tree.info['Name'],tree=self.info['SaveType'],seqname='long',maxnamelen=1000)
        #else: treeMenu(out,mainlog,cmd_list,tree)

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error:', sys.exc_info()[0])
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
