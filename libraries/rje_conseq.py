#!/usr/local/bin/python

# rje_conseq - Sequence Conservation Methods
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_conseq
Description:  Sequence Conservation Methods
Version:      0.1
Last Edit:    18/05/06
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Sequence Conservation Methods.

Commandline:
    rank=T/F        : Whether to calculate ranks
    trimtrunc=T/F   : Whether to trim the leading and trailing gaps (within groups) -> change to X [False]
    winsize=X       : Window size for window scores

Uses general modules: math, os, string, sys, time
Uses RJE modules: rje
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import math
import os
import string
import sys
import time
#############################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#############################################################################################################################
### History
# 0.0 - Initial Compilation.
# 0.1 - Modified general methods
#############################################################################################################################
### Major Functionality to Add
# [ ] More conservation measures?
#############################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_CONSEQ'
        version = "0.0"
        last_edit = "July 05"  
        description = "Sequence Conservation Module"
        author = "Dr Richard J. Edwards."
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#############################################################################################################################
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
#############################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    ### <0> ### Objects setup
    try:
        ## <a> ## Initial Command Setup & Info
        cmd_list = sys.argv[1:]
        info = makeInfo()
        cmd_list = rje.getCmdList(cmd_list,info=info)      ### Load defaults from program.ini
        ## <b> ## Out object
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#############################################################################################################################
methodlist = ['Info','PCon_Abs','PCon_Mean','QPCon_Mean','QPCon_Abs'] #,'QPCon_Mean_All']
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#############################################################################################################################

#############################################################################################################################
### SeqStat Class: 
#############################################################################################################################
class SeqStat(rje.RJE_Object):     
    '''
    SeqStat (Seq Conservation) Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of Method
    
    Opt:boolean
    - TrimTrunc = whether to trim the leading and trailing gaps (within groups) -> change to X [False]
        >> trimtrunc=T/F
    - Rank = whether to also calculate ranks for stats

    Stat:numeric
    - WinSize = Window size for window scores
        >> winsize=X   

    Obj:RJE_Objects
    - Tree = Tree Object with sequences mapped and groupings made
    - AAProp = AAProp Object
    '''
    ### Attributes
    score = {}      # dictionary of scores (keys = Methods, values = lists of scores for residues)
    alnwin = {}     # window scores - (a) the whole alignment,
    qrywin = {}     # window scores - (b) the Query protein of interest (if given) and
    rank = {}       # ranks of scores (keys = Methods, values = lists of scores for residues)
    alnrankwin = {}     # window scores - (a) the whole alignment,
    qryrankwin = {}     # window scores - (b) the Query protein of interest (if given) and
    nodeseq = {}    # Dictionary of sequences (keys = Node Objects (or Sequence Objects if no tree), values = sequences)
#############################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#############################################################################################################################
    def __init__(self,log=None,cmd_list=[],tree=None,aaprop=None,rank=False):  
        '''
        RJE_Object:
        > log:Log = rje.Log object
        > cmd_list:List = List of commandline variables

        On intiation, this object:
        - sets the Log object
        - sets verbosity and interactive attributes
        - calls the _setAttributes() method to setup class attributes
        - calls the _cmdList() method to process relevant Commandline Parameters       
        '''
        ### <0> ### Log Object
        self.cmd_list = cmd_list
        if log == None:
            self.log = rje.Log(cmd_list=cmd_list)
        else:
            self.log = log
            
        ### <1> ### Attributes
        self._setGeneralAttributes()
        self._setAttributes()
        self.obj['Tree'] = tree
        self.obj['AAProp'] = aaprop
        self.opt['Rank'] = rank
        self._setNodeSeq()
        
        ### <2> ### Commandline Options
        self._cmdList()
        if self.opt['TrimTrunc']:
            self._trimTrunc()
#############################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Stats:float ['WinSize']
        - Opt:boolean ['TrimTrunc','Rank']
        - Obj:RJE_Object ['Tree','AAProp']
        '''
        ### <a> ### Basics 
        self.infolist = ['Name']
        for info in self.infolist:
            self.info[info] = 'None'
        self.statlist = ['WinSize']
        for stat in self.statlist:
            self.stat[stat] = 0.0
        self.optlist = ['TrimTrunc','Rank']
        for opt in self.optlist:
            self.opt[opt] = False
        self.objlist = ['Tree','AAProp']
        for obj in self.objlist:
            self.obj[obj] = None
        ### <b> ### Defaults
        ### <c> ### Other Attributes
        self.score = {}
        self.rank = {} 
        self.nodeseq = {} 
#############################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='opt',att='TrimTrunc')
                self._cmdRead(cmd,type='opt',att='Rank')
                self._cmdRead(cmd,type='int',att='WinSize')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#############################################################################################################################
    def _setNodeSeq(self):  ### Sets up self.nodeseq
        '''Sets up self.nodeseq.'''
        if self.obj['Tree']:
            for node in self.obj['Tree'].node:
                if node.obj['Sequence']:
                    self.nodeseq[node] = node.obj['Sequence'].info['Sequence']
#############################################################################################################################
    def _trimTrunc(self):   ### Trims the leading and trailing gaps (within groups) -> change to X 
        '''Trims the leading and trailing gaps (within groups) -> change to X.'''
        try:
            if self.obj['Tree'] and self.obj['Tree'].obj['SeqList']:
                for fam in self.obj['Tree'].subfam:
                    clade = self.obj['Tree']._nodeClade(fam)
                    ## <i> ## Leading spaces
                    r = 0
                    while r < self.obj['Tree'].obj['SeqList'].seq[0].seqLen():
                        gap = False
                        res = False
                        for node in clade:
                            if node in self.nodeseq.keys():
                                if self.nodeseq[node][r] == '-':
                                    gap = True
                                else:
                                    res = True
                        if gap == False:    # End of leading spaces
                            break
                        for node in clade:
                            if res and node in self.nodeseq.keys() and self.nodeseq[node][r] == '-':
                                self.nodeseq[node][r] = 'X'
                        r += 1
                    ## <ii> ## Trailing spaces
                    r = self.obj['Tree'].obj['SeqList'].seq[0].seqLen()
                    while r > 0:
                        r -= 1
                        gap = False
                        res = False
                        for node in clade:
                            if node in self.nodeseq.keys():
                                if self.nodeseq[node][r] == '-':
                                    gap = True
                                else:
                                    res = True
                        if gap == False:    # End of leading spaces
                            break
                        for node in clade:
                            if res and node in self.nodeseq.keys() and self.nodeseq[node][r] == '-':
                                self.nodeseq[node][r] = 'X'
        except:
            self.log.errorLog('Major problem during _trimTrunc().')
#############################################################################################################################
    ### <2> ### Class Methods
#############################################################################################################################
    def nodeScore(self,node=None,query=None,methods=[]):   ### Sends appropriate seqlist to self.calcScore()
        '''
        Sends appropriate seqlist to self.calcScore().
        >> node:Node Object - if None, sends all sequences.
        >> query:key from self.nodeseq that corresponds to Query of group (if appropriate - for QPCon)
        >> methods:list of str = methods of conservation calculation
        '''
        if method:
            self.info['Name'] = method
        else:
            method = self.info['Name']
        if query in self.nodeseq.keys():
            qseq = self.nodeseq[query]
        if self.obj['Tree'] == None or node == None:
            self.calcScore(self.nodeseq.values(),qseq,methods)
        else:
            seqlist = []
            clade = self.obj['Tree']._nodeClade(node)
            for cnode in clade:
                seqlist.append(self.nodeseq[cnode])
            self.calcScore(seqlist,qseq,methods)
#############################################################################################################################
    def calcScore(self,seqlist=[],query=None,methods=[],nonself=[]):     ### Calculates Specificity Scores determined by methods 
        '''
        Calculates Specificity Scores determined by methods and populates dictionaries.
        >> seqlist:list of Sequence Objects for Conservation Calculation
        >> query:Sequence Object (if appropriate for method, e.g. QPCon)
        >> method:str = method of conservation calculation
        >> nonself:list of Sequence Objects that are to be compared to seqlist (if they are not it seqlist).
            - PConMean only
        '''
        if methods == [] or 'all' in methods:
            methods = methodlist
        try:
            ### <1> ### Set up
            ## <a> ## Integrity checks
            _stage ='<1a> Setup - integrity checks'
            if seqlist == []:
                if self.obj['Tree'] == None or self.obj['Tree'].obj['SeqList'] == None:
                    self.log.errorLog('No Sequences Given and Tree or Tree Sequences missing!',printerror=False)
                    raise ValueError
                else:
                    seqlist = self.obj['Tree'].obj['SeqList'].seq
            if seqlist == []:
                self.log.errorLog('No Sequences Given!',printerror=False)
                raise ValueError
            for method in methods[0:]:
                if method not in methodlist:
                    self.log.errorLog('Method \'%s\' not recognised!' % method,printerror=False)
                    methods.remove(method)
                if method == 'QPCon_Mean' and query == None:
                    self.log.errorLog('Method QPCon called but no Query Sequence!',printerror=False)
                    methods.remove(method)
                if method.find('PCon') >= 0 and self.obj['AAProp'] == None:
                    self.log.errorLog('Method PCon called but no AA Property Matrix!',printerror=False)
                    methods.remove(method)
            ## <b> ## Setup variables
            _stage ='<1b> Setup variables'
            tree = self.obj['Tree']
            #print 'Tree', tree
            aaprop = self.obj['AAProp']
            #print 'AAProp', aaprop
            propx = len(aaprop.prop.keys())
            #print propx
            for method in methods:
                self.score[method] = []      # List of scores to return
            #print seqlist
            seqlen = seqlist[0].seqLen()
            ## <c> ## seqlist and nonself
            for seq in seqlist:
                if seq in nonself:
                    nonself.remove(seq)
            if nonself == []:
                nonself = seqlist[0:]
            
            self.verbose(0,3,'Calculating %s scores...' % methods,0)

            ### <2> ### Conservation methods
            ### <a> ### Information Content Calculations
            _stage ='<2a> Conservation - information content'
            if string.join(methods).find('Info') >= 0:
                ## <i> ## Setup
                self.verbose(1,3,'...Information Content (%s)' % seqlist[0].info['Type'],0)
                aax = 21    #!# Note that Gaps are treated as a separate AA and Xs are ignored!
                if seqlist[0].info['Type'] == 'DNA':
                    aax = 5 #!# Note that Gaps are treated as a separate nt and Xs are ignored!
                ## <ii> ## Each residue
                for r in range(seqlen):
                    aacount = {}
                    info = 1.0
                    seqx = 0
                    for seq in seqlist:
                        a = seq.info['Sequence'][r]
                        if a != 'X':
                            seqx += 1
                            if aacount.has_key(a):
                                aacount[a] += 1
                            else:
                                aacount[a] = 1.0
                    ilog = math.log(aax)
                    if len(aacount.keys()) > aax:
                        ilog = math.log(len(aacount.keys()))
                    for aa in aacount.keys():
                        if seqx > 0:
                            info += (aacount[aa] / seqx) * math.log(aacount[aa] / seqx) / ilog
                    self.score['Info'].append(info)

            ### <b> ### Property Conservation                    
            _stage ='<2b> Conservation - aa properties'
            if string.join(methods).find('PCon') >= 0:
                ## <i> ## Setup
                self.verbose(1,3,'...Property Conservation',0)
                propx = len(self.obj['AAProp'].prop.keys())
                ## <ii> ## Work through residues
                for r in range(0,seqlen):
                    if string.join(methods).find('PCon_Abs') >= 0:
                        pcon = propx    # Abs PCon
                        qpcon= propx
                        for property in self.obj['AAProp'].prop.keys():
                            p = -1
                            qp = -1
                            if query:
                                qp = self.obj['AAProp'].prop[property][query.info['Sequence'][r]]
                            for seq in seqlist:
                                a = seq.info['Sequence'][r]
                                if a != 'X' and a != '-':
                                    if p < 0:
                                        p = self.obj['AAProp'].prop[property][a]
                                    elif self.obj['AAProp'].prop[property][a] != p:
                                        pcon -= 1
                                        break
                            for seq in seqlist:
                                a = seq.info['Sequence'][r]
                                if a != 'X' and a != '-' and self.obj['AAProp'].prop[property][a] != qp:
                                    qpcon -= 1
                                    break
                        if 'PCon_Abs' in methods:
                            self.score['PCon_Abs'].append(pcon)
                        if string.join(methods).find('QPCon_Abs') >= 0:
                            self.score['QPCon_Abs'].append(qpcon)
                    if string.join(methods).find('PCon_Mean') >= 0:
                        mpcon = 0.0     # Mean PCon
                        mpx = 0
                        qpcon = 0.0     # Query Mean PCon
                        qpx = 0
                        for seq in seqlist:
                            a1 = seq.info['Sequence'][r]
                            for otherseq in nonself:
                                if seq == otherseq:
                                    continue
                                a2 = otherseq.info['Sequence'][r]
                                if a1 != 'X' and a2 != 'X':         #!# Xs are ignored
                                    mpx += 1
                                    mpcon += self.obj['AAProp'].pdif['%s%s' % (a1,a2)]
                                    if query and seq == query:
                                        qpx += 1
                                        qpcon += self.obj['AAProp'].pdif['%s%s' % (a1,a2)]
                        if mpx > 0:
                            mpcon /= mpx
                            mpcon = len(self.obj['AAProp'].prop.keys()) - mpcon
                        if 'PCon_Mean' in methods:
                            self.score['PCon_Mean'].append(mpcon)                            
                        if qpx > 0:
                            qpcon /= mpx
                            qpcon = len(self.obj['AAProp'].prop.keys()) - qpcon
                        if string.join(methods).find('QPCon_Mean') >= 0:
                            self.score['QPCon_Mean'].append(qpcon)                            

            ### <3> ### Ranks
            _stage ='<3> Ranks'
            if self.opt['Rank']:
                self.verbose(1,1,'...ranks',0)
            for method in methods[0:]:
                if len(self.score[method]) > 0:
                    self.rank[method] = self.rankScore(self.score[method])
                else:   # Remove methods that have no scores
                    methods.remove(method)
                    self.score.pop(method)
            
            ### <4> ### Windows
            _stage ='<4> Windows'
            self.verbose(1,1,'...win(%d)' % self.stat['WinSize'],0)
            for method in methods:
                ## <i> ## Alignment
                self.alnwin[method] = self.winScore(self.score[method])
                self.alnrankwin[method] = self.winScore(self.rank[method])
                ## <ii> ## Query
                qnode = None
                if query:
                    for node in tree.node:
                        if node.obj['Sequence'] == query:
                            qnode = node
                if qnode:
                    self.qrywin[method] = self.winScore(self.score[method],node=qnode)
                    self.qryrankwin[method] = self.winScore(self.rank[method],node=qnode)

            self.verbose(0,1,'...Done!',1)
            
        except:
            self.log.errorLog('Major problem during calcScore(%s): %s.' % (_stage,methods))
#############################################################################################################################
    def winScore(self,scorelist,node=None):  ### Returns window scores as list for given node
        '''
        Returns rank of scores as list.
        >> scorelist = list of scores
        >> node = Node Object (whole alignment if None) (Can be sequence object if no tree)
        << winlist:list of window mean scores (None if no residue at given position)
        '''
        try:
            ### <a> ### Setup
            self.stat['WinSize'] = int(self.stat['WinSize'])
            #print self.stat['WinSize'] 
            ## <i> ## Sequence
            if self.stat['WinSize'] < 1:    # No windows!
                return [None] * len(scorelist)
            elif self.stat['WinSize'] == 1:    # No need for windows!
                return scorelist
            if node:
                sequence = self.nodeseq[node]
            else:
                sequence = 'X' * len(scorelist)
            if len(sequence) != len(scorelist):
                self.log.errorLog('Sequence length (%d) does not match scorelist length (%d)!' % (len(sequence),len(scorelist)),printerror=False)
                raise ValueError
                           
            ### <b> ### Windows
            half = int(self.stat['WinSize']/2)
            start = 0
            end = 0     # End of window
            winlist = []
            i = 0
            winpeaked = False
            while i < len(scorelist):   # Start of window
                ## <i> ## Next residue needing score
                if sequence[i] == '-':  # Gap - no window score!
                    winlist.append(None)
                    i += 1
                    continue
                ## <ii> ## End of window
                x = 0   # Count of residues in window
                end = i + 1 # Count of positions in window
                while end < len(scorelist) and x < half:
                    if sequence[end] != '-':
                        x += 1
                    end += 1
                ## <iii> ## Window score
                x = 0
                win = 0
                for w in range(start,end):
                    if sequence[w] != '-':
                        x += 1
                        win += scorelist[w]
                if x > 0:
                    win /= x
                #if i < 30 or (len(scorelist)-i) < 30:
                self.verbose(3,4,'Res %d, Start %d, End %d, Win %d, Mean = %f' % (i,start,end,x,win),1)
                winlist.append(win)
                ## <iv> ## New start
                if x == self.stat['WinSize']:
                    winpeaked = True
                while winpeaked and (start+1) < len(scorelist):   # Start of window
                    start += 1
                    if sequence[start] != '-':  
                        break
                ## <v> ## Next residue
                i += 1
            #self.verbose(1,4,sequence[-30:],2)
            #print node, len(scorelist), len(winlist)
            if len(scorelist) != len(winlist):
                self.log.errorLog('len(scorelist) (%d) != len(winlist) (%d).' % (len(scorelist), len(winlist)),printerror=False)
                raise ValueError
            return winlist    
        except:
            self.log.errorLog('Major problem during winScore().')
            return [None] * len(scorelist)
#############################################################################################################################
    def rankScore(self,scorelist):  ### Returns rank of scores as list
        '''
        Returns rank of scores as list.
        >> scorelist = list of scores
        << ranklist:list of ranks (0 = Lowest, 1 = Highest)
        '''
        if self.opt['Rank']:
            return rje.rankList(scorelist)
        else:
            return [0] * len(scorelist)
#############################################################################################################################
### End of SeqStat
#############################################################################################################################
    
#############################################################################################################################
## End of SECTION II
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                           #
#############################################################################################################################
def methodDic(log=None,cmd_list=[],tree=None,aaprop=None):     ### Returns a dictionary of SeqStat objects
    '''
    Returns a dictionary of full SeqStat objects.
    >> log:Log Object
    >> cmd_list:list of commandline arguments
    >> tree:Tree Object
    >> aaprop:AA Prop Object
    << methoddic:dictionary of SeqStat objects with methods as keys.
    '''
    methoddic = {}
    seqstat = SeqStat(log=log,cmd_list=cmd_list,tree=tree,aaprop=aaprop,rank=True)
    for method in methodlist:
        methoddic[method] = copy.deepcopy(seqstat)
        methoddic[method].info['Name'] = method
        methoddic[method].nodeScore(node=None)
    return methoddic
#############################################################################################################################
### END OF SECTION III
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                               #
#############################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:        
        print 'Not for standalone running.'
        print 'grep rje_conseq *.py for other modules that will call conseq.py'

    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV
#############################################################################################################################
