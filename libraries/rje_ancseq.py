#!/usr/local/bin/python

# rje_ancseq - Ancestral Sequence Prediction Module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_ancseq
Description:  Ancestral Sequence Prediction Module
Version:      1.3
Last Edit:    24/07/13
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the objects and methods for ancestral sequence prediction. Currently, only GASP (Edwards & Shields
    2004) is implemented. Other methods may be incorporated in the future.

GASP Commandline:
    fixpam=X\t: PAM distance fixed to X [0].
    rarecut=X\t: Rare aa cut-off [0.05].
    fixup=T/F\t: Fix AAs on way up (keep probabilities) [True].
    fixdown=T/F\t: Fix AAs on initial pass down tree [False].
    ordered=T/F\t: Order ancestral sequence output by node number [False].
    pamtree=T/F\t: Calculate and output ancestral tree with PAM distances [True].
    desconly=T/F\t: Limits ancestral AAs to those found in descendants [True].
    xpass=X\t: How many extra passes to make down & up tree after initial GASP [1].

Classes:
    Gasp(log=None,cmd_list=[],tree=None,ancfile='gasp'):
        - Handles main GASP algorithm.
        >> log:Log = rje.Log object
        >> cmd_list:List = List of commandline variables
        >> tree:Tree = rje_tree.Tree Object
        >> ancfile:str = output filename (basefile))        
    GaspNode(realnode,alphabet,log):
        - Used by Gasp Class to handle specific node data during GASP.
        >> realnode:Node Object (rje_tree.py)
        >> alphabet:list of amino acids for use in GASP
        >> log:Log Object

Uses general modules: copy, sys, time
Uses RJE modules: rje, rje_pam
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import string
import os, sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_pam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - No Out Object in Objects
    # 1.0 - Better documentation for GASP V:1.2
    # 1.1 - Neatened up
    # 1.2 - Added 'RST Text' to gaspnode.rst list containing RST file-style data
    # 1.3 - Changed "biproblem" error handling in gaspProbs()
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [?] : GASP Menu method (select options etc.)
    # [ ] : CODEML and others?
    # [ ] : Finish neatening output (gasp() and gaspProbs())
    # [ ] : Upgrade to new module type.
    # [ ] : Integrate FastML predictions.
    # [ ] : Update tree output styles.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_ANCSEQ'
        version = '1.3'
        last_edit = 'July 13'
        description = 'Ancestral Sequence Prediction Module'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print('Problem making Info object.')
        raise
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
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
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
        print('Major Problem with cmdHelp()')
#########################################################################################################################
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
        print('Problem during initial setup.')
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Gasp Class - Ancestral Prediction Object                                                                #
#########################################################################################################################
class Gasp(rje.RJE_Object):     ### Class for performing GASP.
    '''
    GASP Class. Author: Rich Edwards (2005).
    GASP: Gapped Ancestral Sequence Prediction.
    
    Info:str
    - Name = Name (basefile)
    
    Opt:boolean
    - FixUp [True] = Fix AAs on way up (keep probabilities)
    - FixDown [False] = Fix AAs on initial pass down tree
    - Ordered [False] = Order ancestral sequence output by node number
    - PamTree [True] = Whether to output *.nsf & *.txt trees
    - DescOnly [False] = Limits ancestral AAs to those found in descendants
    - RST [False] = Produces RST-style AA probability text

    Stat:numeric
    - FixPam [0] = PAM distance fixed to X
    - RareCut [0.5] = Rare aa cut-off
    - XPass [1] =  How many extra passes to make down & up tree after initial GASP.  
    
    Obj:RJE_Objects
    - Tree:rje_tree.Tree Object (Uses SeqList and PAM objects of Tree)
    '''
    ## Other
    aafreq = []
    gaspnode = {}   # Dictionary of GaspNode objects relating to Tree Nodes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def __init__(self,log=None,cmd_list=[],tree=None,ancfile='gasp'):  
        '''
        RJE_Object:
        > log:Log = rje.Log object
        > cmd_list:List = List of commandline variables
        > tree:Tree = rje_tree.Tree Object
        > ancfile:str = output filename (basefile)

        On intiation, this object:
        - sets the Log object
        - calls the _setAttributes() method to setup class attributes
        - calls the _cmdList() method to process relevant Commandline Parameters       
        '''
        ### <0> ### Log Object
        self.cmd_list = cmd_list
        if log: self.log = log
        else: self.log = rje.Log(cmd_list=cmd_list)
            
        ### <1> ### Attributes
        self._setGeneralAttributes()
        self._setAttributes()
        self.obj['Tree'] = tree
        self.info['Name'] = ancfile        
        
        ### <2> ### Commandline Options
        self._cmdList()
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Stats:float ['FixPam','RareCut','XPass']
        - Opt:boolean ['FixUp','FixDown','Ordered','PamTree','DescOnly','RST']
        - Obj:RJE_Object ['Tree']
        '''
        ### <a> ### Basics
        try:
            self.infolist = ['Name']
            self.statlist = ['FixPam','RareCut','XPass']
            self.optlist = ['FixUp','FixDown','Ordered','PamTree','DescOnly','RST']
            self.objlist = ['Tree']
            ### <b> ### Defaults
            self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=False,setdict=False)
            self.info['Name'] = 'gasp'
            self.stat['FixPam'] = 0
            self.stat['RareCut'] = 0.05
            self.stat['XPass'] = 1
            self.opt['FixUp'] = True
            self.opt['PamTree'] = True
            ### <c> ### Other Attributes
            self.aafreq = []
            self.gaspnode = {}  # Dictionary using tree.nodes as key
        except:
            raise
#########################################################################################################################
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
                self._cmdRead(cmd,type='int',att='FixPam')
                self._cmdRead(cmd,type='stat',att='RareCut')
                self._cmdReadList(cmd,'opt',['FixUp','FixDown','Ordered','PamTree','DescOnly','RST'])
                self._cmdRead(cmd,type='int',att='XPass')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        return
#########################################################################################################################
    def details(self):  ### Prints Details to screen
        '''Returns object details as text.'''
        try:
            ### <a> ### Summary
            details = 'GASP (%s)\n' % (self.info['Name'])
            ### <b> ### Info
            for info in self.infolist:
                if info != 'Name' and info != 'Type' and self.info[info] != '':
                    details += '%s: %s\n' % (info,self.info[info])
            ### <c> ### Stats
            if len(self.statlist) > 0:
                details += 'Stats: [ '
                for stat in self.statlist:
                    details += '%s:%s; ' % (stat,self.stat[stat])
                details += ']\n'
            ### <d> ### Options
            if len(self.optlist) > 0:
                details += 'Options: [ '
                for option in self.optlist:
                    details += '%s:%s; ' % (option,self.opt[option])
                details += ']\n'
            details += '\n'
            ### <f> ### Return
            return details
        except:
            self.log.errorLog('Major problem with Gasp.details().')
            return "Major problem with details(): %s\n" % sys.exc_info()[0]
#########################################################################################################################
    def gasp(self):         ### Performs GASP: Gapped Ancestral Sequence Prediction
        '''Performs GASP: Gapped Ancestral Sequence Prediction.'''
        try:
            ### <a> ### Preparation
            self.obj['Tree'].cmd_list.append('unkspec=T')
            self.obj['Tree'].obj['SeqList'].opt['UnkSpec'] = True
            ## <i> ## Screen Output
            self.verbose(0,3,"\nMaking Ancestral Sequences",0)
            if self.stat['FixPam'] > 0:
                self.verbose(0,3,"- Fixed PAM%d" % self.stat['FixPam'],1)
            else:
                self.verbose(0,3,"- Variable PAM Weighting",1)
            ## <ii> ## PAM Matrix Setup
            try:
                if self.obj['Tree'].obj['PAM'] == None:
                    self.obj['Tree'].obj['PAM'] = rje_pam.PamCtrl(log=self.log,cmd_list=self.cmd_list)
                if self.stat['FixPam'] <= 0:
                    maxblen = 0
                    for b in self.obj['Tree'].branch:
                        if b.stat['Length'] > maxblen:
                            maxblen = b.stat['Length']
                    self.verbose(1,3,'Max Branch Length = %f: ' % maxblen,0)
                    maxblen = int(maxblen*100)+1
                else:
                    maxblen = self.stat['FixPam']
                self.verbose(1,3,'Max PAM = %d' % maxblen,1)
                #print tree.pam.getPamMax(), maxblen
                if self.obj['Tree'].obj['PAM'].stat['PamMax'] < maxblen:
                    #print 'Upping PAM!'
                    self.obj['Tree'].obj['PAM'].stat['PamMax'] = maxblen
                    self.obj['Tree'].obj['PAM'].pamUp()
            except:
                self.log.errorLog("Fatal run Exception during PAM Matrix Setup\n")
                raise
            ##<iii> ## AA Freqs
            aalist = self.obj['Tree'].obj['PAM'].alphabet
            self.verbose(1,3,aalist,1)
            if aalist.count('-') == 0:
                aalist.append('-')
            if aalist.count('X') == 0:
                aalist.append('X')
            self.aafreq = self.obj['Tree'].obj['SeqList'].aaFreq(alphabet=aalist)
            self.aafreq['-'] = 0.0
            self.aafreq['X'] = 0.0
            #tree.deBug(aafreq)
            
            ### <b> ### Terminal sequences - probabilities etc. are known (sequences are known!)
            self.gaspnode = {}   # Array of GaspNode objects
            for node in self.obj['Tree'].node:
                ## <i> ## Check Sequence Exists
                if node.stat['ID'] > self.obj['Tree'].stat['SeqNum']:
                    if node.obj['Sequence'] == None:
                        self.obj['Tree'].obj['SeqList']._addSeq(node.info['Name'],'X' * self.obj['Tree'].obj['SeqList'].seq[0].seqLen())
                        node.obj['Sequence'] = self.obj['Tree'].obj['SeqList'].seq[-1]
                ## <ii> ## Create GaspNode object
                self.gaspnode[node] = GaspNode(node,aalist,self.log)
                ## <iii> ## Termini
                if node.stat['ID'] <= self.obj['Tree'].stat['SeqNum']:
                    self.gaspnode[node].probFromSeq()
                    #print s, len(gaspnode[s].sequence), gaspnode[s].ancfix
                    self.gaspnode[node].ancfix = [True] * len(node.obj['Sequence'].info['Sequence'])    
            
            ### <c> ### GASP 1: Gap Status
            self._gapStatus()

            ##  <d>  ## From tips to root
            #X#self.verbose(0,4,"GASP",0)
            aalist.remove('-')
            if aalist.count('X') > 0:
                aalist.remove('X')
            self._gaspProbs(aalist=aalist,useanc=False,dir='down',aaprobs=True,aasub=self.opt['FixDown'],aafix=self.opt['FixDown'])
            if self.opt['FixDown']:
                self.obj['Tree'].ancSeqOut(file='%s.anc.fas' % self.info['Name'],ordered=self.opt['Ordered'])
                return  
            # Should now have matrix of aa probabilities right back to root...
            ##  <b>  ## Fix Root
            self._gaspProbs(aalist=aalist,useanc=False,dir='root',aaprobs=False,aasub=True,aafix=self.opt['FixUp'])
            ##  <c>  ## Back up tree using all 3 branches
            self._gaspProbs(aalist=aalist,useanc=True,dir='up',aaprobs=True,aasub=True,aafix=self.opt['FixUp'])

            ##  <d>  ## Back down tree with all 3 branches to soften 'outgroup sweep' near root
            for x in range(self.stat['XPass']):
                #X#self.verbose(0,4,":%d:" % (x+1),0)
                self._gaspProbs(aalist=aalist,useanc=True,dir='down',aaprobs=True,aasub=False,aafix=False,gpass=(x+1))
                self._gaspProbs(aalist=aalist,useanc=True,dir='down',aaprobs=True,aasub=True,aafix=True,gpass=(x+1))
                        
            ### <4> ### Finished => Save
            for node in self.obj['Tree'].node:
                node.obj['Sequence'].info['Sequence'] = self.gaspnode[node].sequence
            #X#self.verbose(0,2,"Done!",1)
            self.log.printLog('\r#GASP','Gapped Ancestral Sequence Prediction Complete.')
            self.obj['Tree'].ancSeqOut(file='%s.anc.fas' % self.info['Name'],ordered=self.opt['Ordered'])

            ### <5> ### PAM Distances & PAM Tree
            if self.opt['PamTree']:
                try:
                    self.obj['Tree'].branchPam()
                    self.obj['Tree'].saveTree(filename='%s.anc.nsf' % self.info['Name'],type='nsf',seqnum=1,seqname='short',maxnamelen=127,blen='pam',bootstraps='node',multiline=1)
                    self.obj['Tree'].textTree(seqnum=1,seqname='short',maxnamelen=30,nodename='short',showboot=1,showlen='branch',blen='pam',scale=4,spacer=1,compress=False)
                    self.obj['Tree'].textTree(filename='%s.anc.txt' % self.info['Name'],seqnum=1,seqname='short',maxnamelen=30,nodename='short',showboot=1,showlen='branch',blen='pam',scale=4,spacer=1,compress=False)
                except:
                    self.log.errorLog("Major Problem with PAM Tree.")
                    raise

            ### <6> ### RST Output
            if self.opt['RST']:
                rstfile ='%s.rst' % self.info['Name']
                rje.backup(self,rstfile)
                RST = open(rstfile,'a')
                RST.write('Supplemental results for GASP - main output %s.anc.fas\n\n' % self.info['Name'])
                for node in self.obj['Tree'].node[self.obj['Tree'].stat['SeqNum']:]:
                    gn = self.gaspnode[node]
                    RST.write('%s\n\n' % string.join(gn.rst,'\n'))
                RST.close()
                self.log.printLog('RST output %s.rst complete.' % self.info['Name'])
        except:
            self.log.errorLog('Fatal Error during GASP.')
            raise
#########################################################################################################################
    def _gaspProbs(self,aalist,useanc=False,dir='down',aaprobs=True,aasub=False,aafix=False,gpass=0):  ### Work through tree calculating AA probs
        '''
        Work through tree calculating AA probs.
        >> aalist:list = alphabet
        >> useanc:boolean [False] = whether to consider predicted ancestor
        >> dir:str ['down'] = direction to move through tree ('up','down','rootonly')
        >> aaprobs:boolean [True] = whether to calculate AA probabilities
        >> aasub:boolean [False] = whether to change sequence to most likely AA
        >> aafix:boolean [False] = whether to fix most likely AA as probability 1.0
        >> gpass:int [0] = Extra pass (for log output only)
        '''
        try:
            ### Setup ###
            tree = self.obj['Tree']
            seqlen = tree.node[0].obj['Sequence'].seqLen()
            forloop = tree.node[tree.stat['SeqNum']:]
            if dir == 'up':
                forloop.reverse()
                forloop = forloop[1:]
            elif dir == 'root':
                forloop = [forloop[-1]]
            logtxt = '(%s)'% dir
            if aasub:
                logtxt = 'substitutions (%s)'% dir
            elif aaprobs:
                logtxt = 'probabilities (%s)'% dir
            if gpass > 0:
                logtxt = '[Extra Pass %d] %s' % (gpass,logtxt)
            rstseq = [''] * seqlen  # RST sequences for each column of alignment
            for node in tree.node[:tree.stat['SeqNum']]:
                for r in range(seqlen):
                    rstseq[r] += node.obj['Sequence'].info['Sequence'][r]
                            
            ### <b> ### Each gaspnode in turn
            for node in forloop:
                self.log.printLog('\r#GASP','Calculating GASP %s: %.1f%%' % (logtxt,100.0 * float(forloop.index(node))/len(forloop)),log=False,newline=False)
                an = self.gaspnode[node]
                if node == tree.node[-1]:   # Root
                    desc_seq = tree.node
                else:                    
                    clades = tree.branchClades(node.ancBranch())
                    desc_seq = clades[1]
                ## <i> ## AA probabilities
                if aaprobs:
                    # Establish relevant nodes and PAM distances
                    anc = node.ancNode()
                    desc = node.neighbours(ignore=[anc])
                    treenodes = desc
                    n = [self.gaspnode[desc[0]],self.gaspnode[desc[1]]] # Gasp nodes
                    p = {}
                    if anc != None:
                        treenodes.append(anc)   # Tree Nodes
                        n.append(self.gaspnode[anc])
                    for tn in treenodes:
                        p[tn] = self.stat['FixPam']  # pam distances
                        if self.stat['FixPam'] == 0:
                            p[tn] = int(node.link(tn).stat['Length'] * 100) + 1
                    # Calculate Probabilities for each residue
                    bigproblem = False
                    for r in range(seqlen):
                        if (an.ancgap[r] == 0) and (an.ancfix[r] == False):   # No gap & not fixed
                            if n[0].ancfix[r] and n[1].ancfix[r] and (n[0].sequence[r] == n[1].sequence[r]):    # Fixed in descendants
                                a = n[0].sequence[r]
                                an.sequence = rje.strSub(an.sequence,r,r,a)
                                an.ancaap[r][a] = 1.0
                                an.ancfix[r] = True
                            else:   # Variation
                                for a in aalist:
                                    aap = [0.0] * 3
                                    #print a, aap, 'PAM', p[desc[0]], p[desc[1]]
                                    for d in aalist:
                                        #print a, '->', d, 'PAM', p[desc[0]], tree.obj['PAM'].matrix[p[desc[0]]].pamp[a+d], '*' ,n[0].ancaap[r][d]
                                        try: aap[0] += tree.obj['PAM'].matrix[p[desc[0]]].pamp[a+d] * n[0].ancaap[r][d]
                                        except:
                                            self.errorLog('Big Problem with GASP Prob calculation',quitchoice=True); raise
                                            if bigproblem: self.errorLog('Big Problem with GASP Prob calculation')
                                            else: bigproblem = True
                                        #print a, '->', d, 'PAM', p[desc[1]], tree.obj['PAM'].matrix[p[desc[1]]].pamp[a+d], '*' ,n[1].ancaap[r][d]
                                        try: aap[1] += tree.obj['PAM'].matrix[p[desc[1]]].pamp[a+d] * n[1].ancaap[r][d]
                                        except:
                                            self.errorLog('Big Problem with GASP Prob calculation',quitchoice=True); raise
                                            if bigproblem: self.errorLog('Big Problem with GASP Prob calculation')
                                            else: bigproblem = True
                                        try:
                                            if useanc and anc != None: aap[2] += tree.obj['PAM'].matrix[p[anc]].pamp[d+a] * n[2].ancaap[r][d]
                                        except:
                                            self.errorLog('Big Problem with GASP Prob calculation',quitchoice=True); raise
                                            if bigproblem: self.errorLog('Big Problem with GASP Prob calculation')
                                            else: bigproblem = True
                                    if useanc and anc != None:
                                        an.ancaap[r][a] = (aap[0] + aap[1] + aap[2]) / 3
                                    else:
                                        an.ancaap[r][a] = (aap[0] + aap[1]) / 2
                    # Adjust Probabilities to total 1, using RareCut and DescOnly
                    if dir != 'up':
                        desc_seq = tree.node
                    an.adjustAncAAP(self.stat['RareCut'],desconly=self.opt['DescOnly'],desc=desc_seq)
                    if self.opt['RST']:
                        an.makeRST(rstseq)
                if aasub:
                    for r in range(seqlen):
                        if (an.sequence[r] != '-') & (an.ancfix[r] == False):
                            maxp = 0.0
                            for a in aalist:
                                if an.ancaap[r][a] > maxp:
                                    an.sequence = rje.strSub(an.sequence,r,r,a)
                                    maxp = an.ancaap[r][a]
                                elif (an.ancaap[r][a] == maxp) & (self.aafreq[a] > self.aafreq[an.sequence[r]]):
                                    an.sequence = rje.strSub(an.sequence,r,r,a)
                            if maxp == 0:   # No AA => X
                                an.sequence = rje.strSub(an.sequence,r,r,'X')
                            if aafix:
                                for a in aalist:
                                    if a == an.sequence[r]:
                                        an.ancaap[r][a] = 1.0
                                    else:
                                        an.ancaap[r][a] = 0.0

            ### Finish ###
            self.log.printLog('\r#GASP','Calculated GASP %s: 100.0%%.' % logtxt)
        except:
            self.log.errorLog('Fatal Error during _gaspProbs().')
            raise
#########################################################################################################################
    def _gapStatus(self):   ### Predicts ancestral gap status
        '''Predicts ancestral gap status.'''
        try:
            ### Setup ###
            tree = self.obj['Tree']
            gaspnode = self.gaspnode
            progress = 0.01 * 2 * len(tree.node[tree.stat['SeqNum']:]) + len(gaspnode[tree.node[-1]].sequence)
            p = 0
            
            ###  From tips to root ###
            self.log.printLog('\r#GASP','Calculating Gap Status (Down) ... 0.0%',log=False,newline=False)
            for node in tree.node[tree.stat['SeqNum']:]:
                an = gaspnode[node]
                anc = node.ancNode()
                desc = node.neighbours(ignore=[anc])
                n = [gaspnode[desc[0]],gaspnode[desc[1]]]
                for r in range(len(an.sequence)):
                    an.ancgap[r] = (n[0].ancgap[r] + n[1].ancgap[r])/2
                an.nodegap = 1
                p += 100.0
                self.log.printLog('\r#GASP','Calculating Gap Status (Down) ... %.1f%%' % (p/progress),log=False,newline=False)

            ### Fix Root ###
            an = gaspnode[tree.node[-1]]
            for r in range(len(an.sequence)):
                if an.ancgap[r] <= 0.5:
                    an.sequence = rje.strSub(an.sequence,r,r,'X')
                    an.ancgap[r] = 0.0
                else:
                    an.sequence = rje.strSub(an.sequence,r,r,'-')
                    an.ancgap[r] = 1.0
                p += 100.0
                self.log.printLog('\r#GASP','Calculating Gap Status (Root) ... %.1f%%' % (p/progress),log=False,newline=False)

            ### Back up tree using all 3 branches ###
            uptree = tree.node[tree.stat['SeqNum']:-1]
            uptree.reverse()
            for node in uptree:
                an = gaspnode[node]
                anc = node.ancNode()
                desc = node.neighbours(ignore=[anc])
                n = [gaspnode[desc[0]],gaspnode[desc[1]],gaspnode[anc]]
                for r in range(len(an.sequence)):
                    an.ancgap[r] = (n[0].ancgap[r] + n[1].ancgap[r] + n[2].ancgap[r])/3
                    if an.ancgap[r] <= 0.5:
                        an.sequence = rje.strSub(an.sequence,r,r,'X')
                        an.ancgap[r] = 0.0
                    else:
                        an.sequence = rje.strSub(an.sequence,r,r,'-')
                        an.ancgap[r] = 1.0
                p += 100.0
                self.log.printLog('\r#GASP','Calculating Gap Status ( Up ) ... %.1f%%' % (p/progress),log=False,newline=False)

            self.log.printLog('\r#GASP','Calculation of Gap Status complete.      ')
        except:
            self.log.errorLog('Fatal Error during _gapStatus().')
            raise
#########################################################################################################################
## End of SECTION II: Gasp Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: GaspNode Class - Individual Tree Nodes                                                                 #
#########################################################################################################################
class GaspNode(object):     ### Class for use by GASP only
    '''
    Class for use by GASP only.
    >> node:rje_tree.Node   = 'real' node associated with GaspNode
    >> alphabet:list =[]    = list of letters in sequence - matches PAM matrix
    >> log:rje.Log  = Log object for error messages
    ancaap:list     = list of dictionaries [r]{a}, probability of amino acid a @ residue r
    ancgap:list     = [r], probablility of gap @ residue r
    ancfix:list     = [r], boolean whether residue fixed in descendants
    rst:list        = [r], string containing RST-style probability information
    '''    
    node = None     # rje_tree.Node object
    alphabet=[]     # List of letters in sequence - matches PAM matrix
    log = None      # Log for error messages
    sequence = ''   # Sequence
    ancaap = []     # [r]{a}  = probability of amino acid a @ residue r
    ancgap = []     # [r] = probablility of gap @ residue r
    ancfix = []
    rst = []
#########################################################################################################################
    ### <1> ### Class Initiation: sets attributes                                                                       #
#########################################################################################################################
    def __init__(self,realnode,alphabet,log):
        '''
        Used by Gasp Class to handle specific node data during GASP.
        >> realnode:Node Object (rje_tree.py)
        >> alphabet:list of amino acids for use in GASP
        >> log:Log Object
        '''
        self.node = realnode
        self.alphabet = alphabet
        self.sequence = realnode.obj['Sequence'].info['Sequence']
        self.name = realnode.obj['Sequence'].shortName()
        self.ancaap = []
        self.ancgap = []
        self.ancfix = []
        self._buildAncLists()
        self.log = log
#########################################################################################################################
    def _errorLog(self,errortxt):   ### Tries to report error if log exists, else prints
        '''Tries to report error if log exists, else prints.'''
        if self.log:
            try: self.log.errorLog(errortxt)
            except: self.log.errorLog(errortxt,printerror=False) 
        else: print('!!! Error but no log !!! %s !!!' % errorttxt)
#########################################################################################################################
    def _buildAncLists(self):   ### Builds initial ancaap and ancgap lists.
        '''Builds initial ancaap and ancgap lists.'''
        try:
            for r in range(len(self.sequence)):
                rdic = {}
                for letter in self.alphabet: rdic[letter] = 0.0
                self.ancaap.append(rdic)
                self.ancgap.append(0.0)
                self.ancfix.append(False)
        except:
            self._errorLog('Fatal Error with GaspNode._buildAncLists(%s).' % self.name)
            raise
#########################################################################################################################
    def probFromSeq(self):      ### Makes probability matrices based on own sequence (1 or 0)
        '''Makes probability matrices based on own sequence (1 or 0).'''
        try:
            self._buildAncLists()
            for r in range(len(self.sequence)):
                aa = self.sequence[r]
                self.ancaap[r][aa] = 1.0
                if aa == '-': self.ancgap[r] = 1.0
        except:
            self._errorLog('Fatal Error with GaspNode.probFromSeq().')
            raise
#########################################################################################################################
    def adjustAncAAP(self,rarecut,desconly=False,desc=[]):    ### Adjusts AA frequencies according to rarecut
        '''
        Adjusts AA frequencies according to rarecut.
        >> rarecut:float = Frequency cut-off
        >> desconly:boolean [False] = whether to limit probabilities to AAs found in descendants.
        >> desc:list of descendant terminal nodes
        '''
        for r in range(len(self.sequence)):
            try:
                if self.ancgap[r] == 0:
                    probtot = 0.0
                    desc_alpha = []
                    if desconly:
                        for node in desc:
                            desc_alpha.append(node.obj['Sequence'].info['Sequence'][r])
                    for letter in self.alphabet:
                        if desconly == False or letter in desc_alpha:
                            probtot += self.ancaap[r][letter]
                        else:
                            self.ancaap[r][letter] = 0.0
                    if probtot > 0:
                        for letter in self.alphabet:
                            self.ancaap[r][letter] /= probtot
                            if self.ancaap[r][letter] < rarecut:    # Remove low frequency AAs
                                self.ancaap[r][letter] = 0.0
                    else: self._errorLog('Total AncAAP for node %d, res %d = 0!' % (self.node.stat['ID'],r))
            except: self._errorLog('Major Problem with adjustAncAAP for node %d, res %d = 0!.' % (self.node.stat['ID'],r))
        if rarecut > 0: self.adjustAncAAP(0.0)    # Run first portion again to normalise frequencies
#########################################################################################################################
    def makeRST(self,rstseq):   ### Makes RST text using sequence list and own probablities
        '''
        Makes RST text using sequence list and own probablities and populates self.rst.
        >> rstseq:list = Terminal node sequences
        '''
        try:
            ### Setup ###
            self.rst = [''] * len(self.ancaap)

            ### Make RST Text ###
            for r in range(len(self.ancaap)):
                self.rst[r] = '%7d%7d   %s:' % (r+1,len(rstseq[r]),rstseq[r])
                for letter in self.alphabet:
                    self.rst[r] += ' %s(%.3f)' % (letter,self.ancaap[r][letter])

            ### Finish up ###
            self.rst = ['Prob distribution at node %d, by site' % (self.node.stat['ID']),'','   site  Freq   Data'] + self.rst

        except:
            self._errorLog('Major Problem with makeRST for node %d!' % (self.node.stat['ID']))
            self.rst = ['Major Problem with makeRST for node %d!' % (self.node.stat['ID'])]
#########################################################################################################################
### END OF SECTION III: GaspNode Class                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print('Unexpected error during program setup:', sys.exc_info()[0])
        return
        
    ### Rest of Functionality... ###
    try:        
        print('\n\n *** No standalone functionality! Must be given a tree etc. Run rje_tree.py or gasp.py. *** \n\n')
        
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
        print('Cataclysmic run error:', sys.exc_info()[0])
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
