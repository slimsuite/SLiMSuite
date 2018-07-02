#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_supertree
Description:  RJE SuperTree results tidier
Version:      0.2.0
Last Edit:    22/06/18
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The initial function of this module is to take the consensus tree generated from a supertree algorithm and overlap
    information from the source trees:

    1. Rate each branch for its support in the original trees (absolute or percentage).
    2. Calculate the mean branch length for each branch from when it is observed in the original trees.
    3. Map new node labels onto the tree.

Commandline:

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    supertree=NWK   : Supertree in Newick format []
    sourcetrees=NWK : Source trees in a single file (Newick format) []
    translate=FILE  : Text file containing translations of supertree tips into OTU names []
    cladecounts=FILE: Text file containing lists of tip combinations to counts (clades versus complete sets in trees) []
    annotate=FILE   : Text file of new names for supertree. Will map and count taxa only []
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje
import rje_obj
import rje_tree
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version, mapping source trees onto supertrees and calculating coverage.
    # 0.2.0 - annotate=FILE   : Text file of new names for supertree. Will map and count taxa only []
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add clade counts
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_SUPERTREE', '0.2.0', 'June 2018', '2016')
    description = 'RJE SuperTree results tidier'
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
### SECTION II: SuperTree Class                                                                                         #
#########################################################################################################################
class SuperTree(rje_obj.RJE_Object):
    '''
    SuperTree Class. Author: Rich Edwards (2018).

    Str:str
    - Annotate=FILE   : Text file of new names for supertree. Will map and count taxa only []
    - CladeCounts=FILE: Text file containing lists of tip combinations to counts (clades versus complete sets in trees) []
    - SuperTree=NWK   : Supertree in Newick format []
    - SourceTrees=NWK : Source trees in a single file (Newick format) []
    - Translate=FILE  : Text file containing translations of supertree tips into OTU names []

    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Annotate','CladeCounts','SuperTree','SourceTrees','Translate']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
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
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Annotate','CladeCounts','SuperTree','SourceTrees','Translate'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Annotate'): self.annotate()
            else: self.mapTrees()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        There is currently no specific help available on REST output for this program.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def mapTrees(self):      ### Generic method
        '''
        mapTrees method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #1# Load the master tree and translation table
            tree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            tree.loadTree(self.getStr('SuperTree'))
            #tree.unRoot()
            self.deBug(tree._checkTree())
            #tree.debugTree()
            self.deBug(tree.sumTree())
            #2# Convert names in tree
            if self.getStrLC('Translate'):
                trans = {}
                for line in open(self.getStr('Translate'),'r').readlines():
                    tdata = string.split(line)
                    try: trans[tdata[0]] = tdata[1]
                    except: pass
                self.printLog('#TRANS','%d node translations loaded' % len(trans))
                # Map
                for node in tree.node[:tree.stat['SeqNum']]:
                    try: node.info['Name'] = trans[node.info['Name']]
                    except: self.warnLog('Failed to find "%s" in translations from %s' % (node.info['Name'],self.getStr('Translate')))
            else: self.printLog('#TRANS','No node translations loaded')
            # Clade counts
            cladecount = {}     # Dictionary of {(clade):[clade count, full OTU in tree count]}
            if self.getStrLC('CladeCounts'):
                for line in open(self.getStr('CladeCounts'),'r').readlines():
                    if line[:1] == '#': continue
                    clade = string.split(string.join(string.split(line,',')))
                    clade.sort()
                    if not clade: continue
                    self.printLog('#CLADE','Counting: %s' % string.join(clade,'|'))
                    cladecount[tuple(clade)] = [0,0]
                self.printLog('#CLADE','%d clades loaded to count' % len(cladecount))
            # Full OTU list (in case supertree is an edited version)
            superotu = []
            for node in tree.node[:tree.stat['SeqNum']]:
                superotu.append(node.shortName())


            #3# Add for each branch:
            branchtuples = {}
            for branch in tree.branch:
                # a counter for support and consistency
                branch.stat['Support'] = 0      # Number of trees with a complete clade match on one side
                branch.stat['Consistency'] = 0  # Number of trees with consistent OTUs on both sides of branch
                branch.stat['Coverage'] = 0     # Number of trees with 1+ OTU from each side of branch
                branch.list['Lengths'] = []              # branchlengthlist = []  # List of lengths read for that branch
                # cladetuple = sorted (clade1 OTU,clade2 OTU)
                clades = tree.branchClades(branch,internal=False)  ### Return tuple of lists of sequences either side of branch ('anc','desc')
                clade0 = []; clade1 = []
                for node in clades[0]: clade0.append(node.shortName())
                for node in clades[1]: clade1.append(node.shortName())
                clade0.sort(); clade1.sort()
                bclades = [clade0,clade1]
                bclades.sort()
                btuple = (tuple(bclades[0]),tuple(bclades[1]))
                self.bugPrint('%s: "%s"' % (branch.name(),btuple))
                branchtuples[btuple] = branch   # Make a dictionary of {(clade1,clade2):branch object} from branches
                branch.obj['Tuple'] = btuple

            #4# Read each source tree and
            treex = 0
            for line in open(self.getStr('SourceTrees'),'r').readlines():
                nsftree = line[:line.find(';')+1]
                if not nsftree: continue
                treex += 1
                self.printLog('#TREE','Source tree %d: %s' % (treex,nsftree))
                itree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
                itree.buildTree(nsftree,seqlist=None,type='nsf',postprocess=False)
                itree.unRoot()
                otu = []
                #a# Generate a full list of OTUs.
                for node in itree.node[:itree.stat['SeqNum']]:
                    otu.append(node.shortName())
                otu.sort()
                self.printLog('#OTU','Source tree %d: %d OTU = %s' % (treex,len(otu),string.join(otu,'; ')))
                # Assess Coverage
                ibranchtuples = {}  # Dict of branches with coverage: (reduced clades)
                for (clade0,clade1) in branchtuples:
                    otu0 = rje.listIntersect(otu,clade0)
                    otu0.sort()
                    otu1 = rje.listIntersect(otu,clade1)
                    otu1.sort()
                    # "Coverage" = having OTUs on each side of branch
                    if (otu0 and otu1): # or tuple(otu0) in (clade0,clade1) or tuple(otu1) in (clade0,clade1):
                        sbranch = branchtuples[(clade0,clade1)]
                        sbranch.stat['Coverage'] += 1
                        ibranchtuples[(clade0,clade1)] = (tuple(otu0),tuple(otu1))

                #b# Cycle through each branch and generate sorted (clade1,clade2)
                consistent = [] # List of consistent branches
                supported = []  # List of supported branches
                counted = []    # List of counted clades
                for clade in cladecount:
                    if len(rje.listIntersect(otu,clade)) == len(clade): cladecount[clade][1] += 1
                for branch in itree.branch:
                    # cladetuple = sorted (clade1 OTU,clade2 OTU)
                    clades = itree.branchClades(branch,internal=False)  ### Return tuple of lists of sequences either side of branch ('anc','desc')
                    clade0 = []; clade1 = []
                    for node in clades[0]: clade0.append(node.shortName())
                    for node in clades[1]: clade1.append(node.shortName())
                    [clade0,clade1] = [rje.listIntersect(superotu,clade0),rje.listIntersect(superotu,clade1)]
                    clade0.sort(); clade1.sort()
                    bclades = [clade0,clade1]
                    bclades.sort()
                    btuple = (tuple(bclades[0]),tuple(bclades[1]))
                    for clade in cladecount:
                        if clade in btuple and clade not in counted: counted.append(clade)

                    #c# if (clade1,clade2) in cladetuple:
                    #self.bugPrint('%d::%s: %s' % (treex,branch.name(),btuple))
                    if btuple in branchtuples:
                        # Add length to length lists
                        sbranch = branchtuples[btuple]
                        # Add 1 to branch support and consistency
                        if sbranch not in supported: supported.append(sbranch)
                        else: self.warnLog('Branch supported twice!')
                        if sbranch not in consistent: consistent.append(sbranch)
                        sbranch.list['Lengths'].append(branch.stat['Length'])
                        if len(btuple[0]) > 1 and len(btuple[1]) > 1:
                            self.printLog('#SUPPORT','%s (%s),(%s) = %s' % (branch.name(),string.join(btuple[0],'|'),string.join(btuple[1],'|'),sbranch.name()))
                            self.debug('...')
                    # - else
                    else:
                        # cycle through cladetuple and generate reduced OTU lists based on full tree OTU
                        bx = 0
                        for (clade0,clade1) in ibranchtuples:
                            #self.bugPrint('0: %s' % str(clade0))
                            if 'oldcode' == 'run':
                                otu0 = rje.listIntersect(otu,clade0)
                                if not otu0: continue       # Outside source tree. (No coverage)
                                otu0.sort()
                                #self.bugPrint(otu0)
                                #self.bugPrint('1: %s' % str(clade1))
                                otu1 = rje.listIntersect(otu,clade1)
                                if not otu1: continue       # Outside source tree. (No coverage)
                                otu1.sort()
                                #self.bugPrint(otu1)
                                # assess whether this branch is now supported - +1 consistency if so

                                #!# Check this is right and check the Coverage calculation above - currently mismatch #!#

                                if tuple(otu0) in btuple and tuple(otu1) in btuple:
                                    sbranch = branchtuples[(clade0,clade1)]
                                    if sbranch not in consistent: consistent.append(sbranch)
                                    if tuple(otu0) == clade0 or tuple(otu1) == clade1:  # One complete clade in tree
                                        if sbranch not in supported: supported.append(sbranch)
                                        else:
                                            self.warnLog('Branch supported twice!')
                                        sbranch.list['Lengths'].append(branch.stat['Length'])
                                    #self.deBug('Yes!');
                                    bx += 1
                            #self.bugPrint('%s -> %d'  % (branch.name(),bx))

                            (otu0,otu1) = ibranchtuples[(clade0,clade1)]
                            # assess whether this branch is now supported - +1 consistency if so
                            if otu0 in btuple and otu1 in btuple:
                                sbranch = branchtuples[(clade0,clade1)]
                                if sbranch not in consistent: consistent.append(sbranch)
                                if otu0 in (clade0,clade1) or otu1 in (clade0,clade1):  # One complete clade in tree
                                    if sbranch not in supported: supported.append(sbranch)
                                    else:
                                        self.warnLog('Branch supported twice!')
                                    sbranch.list['Lengths'].append(branch.stat['Length'])
                                    if len(otu1) > 1 and len(otu0) > 1:
                                        self.printLog('#SUPPORT','%s (%s),(%s) =~ (%s),(%s) %s => (%s),(%s)' % (branch.name(),string.join(btuple[0],'|'),string.join(btuple[1],'|'),string.join(clade0,'|'),string.join(clade1,'|'),sbranch.name(),string.join(otu0,'|'),string.join(otu1,'|')))
                                        self.debug('...')
                                #self.deBug('Yes!');
                                bx += 1
                        #self.bugPrint('%s -> %d'  % (branch.name(),bx))
                for sbranch in supported:
                    sbranch.stat['Support'] += 1
                for sbranch in consistent:
                    sbranch.stat['Consistency'] += 1
                for clade in counted: cladecount[clade][0] += 1

            #tree._checkTree()
            #tree.debugTree()
            self.deBug(tree.sumTree())

            #5# Convert lists into mean branch lengths and update branch objects. Root.
            for branch in tree.branch:
                if branch.list['Lengths']: branch.stat['Length'] = sum(branch.list['Lengths'])/float(len(branch.list['Lengths']))
                else: branch.stat['Length'] = 0
                branch.info['Name'] = '%s [%d/%d/%d]' % (branch.info['Name'],branch.stat['Support'],branch.stat['Consistency'],branch.stat['Coverage'])
                self.printLog('#BRANCH',branch.info['Name'])
                self.bugPrint('%s\n%s\n%s' % (branch.name(),branch.list,branch.stat))
            self.deBug(tree.sumTree())
            tree.midRoot()
            self.deBug(tree.sumTree())

            #6# Output tree with (a) support, and (b) consistency counts in place of bootstrap values.
            tree.stat['BootStraps'] = treex
            tree.opt['Bootstrapped'] = True
            for out in ['Support','Consistency','Coverage']:
                tree.info['Basefile'] = '%s.%s' % (self.baseFile(),out.lower())
                for branch in tree.branch:
                    try: branch.stat['Bootstrap'] = branch.stat[out]
                    except: branch.stat['Bootstrap'] = 0
                tree.saveTrees()
            #6b# Update Bootstrap to be percentage of coverage consistent
            #i# NOTE: This will currently round down to integer in output.
            tree.info['Basefile'] = '%s.perc' % (self.baseFile())
            tree.stat['BootStraps'] = 100
            for branch in tree.branch:
                try:
                    if branch.stat['Coverage'] > 0:
                        branch.stat['Bootstrap'] = (100.0 * branch.stat['Consistency'] / branch.stat['Coverage'])
                    else: branch.stat['Bootstrap'] = -1
                except: branch.stat['Bootstrap'] = 0
            tree.saveTrees()


            ### Clade counts
            for clade in rje.sortKeys(cladecount):
                self.printLog('#CLADE','(%s) = %d / %d trees with OTUs.' % (string.join(clade,'|'),cladecount[clade][0],cladecount[clade][1]))

            tree.info['Basefile'] = '%s.edit' % (self.baseFile())
            tree.editTree()
            tree.saveTrees()

            return
        except: self.errorLog('%s.mapTrees error' % self.prog())
#########################################################################################################################
    def annotate(self):      ### Generic method
        '''
        mapTrees method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #1# Load the master tree and translation table
            tree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            tree.loadTree(self.getStr('SuperTree'))
            #tree.unRoot()
            self.deBug(tree._checkTree())
            #tree.debugTree()
            self.deBug(tree.sumTree())
            #2# Convert names in tree
            if self.getStrLC('Annotate'):
                trans = {}
                for line in open(self.getStr('Annotate'),'r').readlines():
                    tdata = string.split(line)
                    try: trans[tdata[0]] = string.join(tdata[1:])
                    except: pass
                self.printLog('#TRANS','%d node translations loaded' % len(trans))
                # Map
                for node in tree.node[:tree.stat['SeqNum']]:
                    try: node.info['Name'] = trans[node.info['Name']]
                    except: self.warnLog('Failed to find "%s" in translations from %s' % (node.info['Name'],self.getStr('Translate')))
            else: self.printLog('#TRANS','No node translations loaded')

            nodecount = {}

            #4# Read each source tree and
            treex = 0
            for line in open(self.getStr('SourceTrees'),'r').readlines():
                nsftree = line[:line.find(';')+1]
                if not nsftree: continue
                treex += 1
                self.printLog('#TREE','Source tree %d: %s' % (treex,nsftree))
                itree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
                itree.buildTree(nsftree,seqlist=None,type='nsf',postprocess=False)
                for node in itree.node:
                    otu = node.name()
                    if otu in trans: otu = trans[otu]
                    if otu not in nodecount: nodecount[otu] = 0
                    nodecount[otu] += 1
            #tree._checkTree()
            #tree.debugTree()
            self.deBug(tree.sumTree())

            for node in tree.node:
                otu = node.name()
                if otu in nodecount: node.info['Name'] = '%s [%d / %.1f%%]' % (otu,nodecount[otu],100.0*nodecount[otu]/treex)

            tree.info['Basefile'] = '%s.annotated' % (self.baseFile())
            tree.editTree()
            tree.saveTrees()

            return
        except: self.errorLog('%s.annotate error' % self.prog())
#########################################################################################################################
### End of SECTION II: SuperTree Class                                                                                  #
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
    try: SuperTree(mainlog,['treeformats=nwk,text,png']+cmd_list).run()

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
