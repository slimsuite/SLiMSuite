#!/usr/local/bin/python

# rje_specificity - Functional Specificity Methods
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_specificity
Description:  Functional Specificity Methods
Version:      0.1
Last Edit:    18/05/06
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function: Functional Specificity methods for BADASP.

Commandline:
    
Uses general modules: os, string, sys, time
Uses RJE modules: rje, rje_conseq
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import os
import string
import sys
import time
#############################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_conseq
#############################################################################################################################
### History
# 0.0 - Initial Working Compilation.
# 0.1 - Working version for BADASP
#############################################################################################################################
### Major Functionality to Add
# [ ] Mutual Information
# [ ] Add an extra attribute = dictionary of funcspec scores rather than just a single one
#############################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_SPECIFICITY'
        version = "0.1"
        last_edit = 'May 06'
        description = "Functional Specificity Prediction Module"
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
methodlist = ['BAD','BADN','BADX','SSC','PDAD','ETA','ETAQ']
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#############################################################################################################################

#############################################################################################################################
### FuncSpec Class: 
#############################################################################################################################
class FuncSpec(rje_conseq.SeqStat):     
    '''
    FuncSpec (Functional Specificity Class. Author: Rich Edwards (2005).
    Based on SeqStat (Seq Conservation) Class or rje_conseq. 

    Info:str
    - Name = Name 
    
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
    famwin = {}     # window scores  - (c) the ancestral sequence of each subfamily;
    famrankwin = {}     # window scores  - (c) the ancestral sequence of each subfamily;
#############################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
    ### .. All in rje_conseq.SeqStat Class
#############################################################################################################################
    ### <2> ### Class Methods
#############################################################################################################################
    def calcScore(self,query=None,methods=[]):     ### Calculates Specificity Scores determined by methods and populates dictionaries
        '''
        Calculates Specificity Scores determined by methods and populates dictionaries.
        >> query:Sequence Object 
        >> methods:list of str = methods of functional specificity to calculate
        '''
        if methods == [] or 'all' in methods:
            methods = methodlist
        try:
            ### <1> ### Set up
            _stage = '<1> Setup'
            ## <a> ## Integrity checks
            _stage = '<1a> Setup - Integrity Checks'
            if self.obj['Tree'] == None or self.obj['Tree'].obj['SeqList'] == None:
                self.log.errorLog('Tree or Tree Sequences missing!',printerror=False)
                raise ValueError
            for method in methods[0:]:
                if method not in methodlist:
                    self.log.errorLog('Method \'%s\' not recognised!' % self.info['Name'],printerror=False)
                    methods.remove(method)
                if self.obj['AAProp'] == None and method in ['BAD','BADN','BADX','SSC','PDAD']:
                    self.log.errorLog('Method %s called but no AA Property Matrix!' % method,printerror=False)
                    methods.remove(method)
                if method != 'BAD' and self.obj['Tree'].groupNum() < 2:     ### Original BAD needs 2 groups, others 2+
                    self.log.errorLog('Method %s needs 2+ subfam defined!' % method,printerror=False)
                    methods.remove(method)
                elif method == 'BAD' and self.obj['Tree'].groupNum() != 2:  ### Original BAD needs 2 groups
                    self.log.errorLog('Method %s needs exactly 2 subfam defined!' % method,printerror=False)
                    methods.remove(method)
                if method in ['BADX','BADN'] and query == None:
                    self.log.errorLog('Method %s needs Query defined!' % method,printerror=False)
                    methods.remove(method)

            ## <b> ## Setup variables
            _stage = '<1b> Setup - variables'
            tree = self.obj['Tree']
            aaprop = self.obj['AAProp']
            propx = len(aaprop.prop.keys())
            for method in methods:
                self.score[method] = []      # List of scores to return
            #print self.stat
            seqlen = len(self.nodeseq.values()[0])
            self.verbose(0,3,'Calculating %s scores... (%d residues)' % (methods,seqlen),0)

            ## <c> ## Subfams etc.
            _stage = '<1c> Setup - Subfams etc.'
            subfam = tree.subfam    # List of nodes designating groups
            ancnode = tree.node[-1] # Ancestral node to Query Group
            outnode = None          # Outgroup to ancestor of Query Node
            clade = {}              # Dictionary of lists of nodes, key = group-defining nodes
            seqs = {}               # Dictionary of lists of sequences, key = group-defining nodes
            queryfam = None         # Query subfamily
            for fam in subfam:
                #self.verbose(1,3,'...%s' % fam.info['Name'],0)
                clade[fam] = tree._nodeClade(fam)
                seqs[fam] = []
                for node in clade[fam]:
                    #self.verbose(1,3,':%s' % node.info['Name'],0)
                    seqs[fam].append(self.nodeseq[node])
                    if query and query == node.obj['Sequence'] and fam != tree.node[-1]:
                        ancnode = fam.ancNode()
                        outnode = tree.outGroupNode(fam)
                        queryfam = fam
            #print 'BADX between %s and %s?' % (queryfam.info['Name'],outnode.info['Name'])
            self.verbose(1,3,'...fams',0)
            clade[ancnode] = tree._nodeClade(ancnode)
            seqs[ancnode] = []
            for node in clade[ancnode]:
                seqs[ancnode].append(self.nodeseq[node])
                                
            ### <2> ### Specificity methods
            _stage = '<2> Specificity Scores'
            ### <a> ### Burst After Duplication (BAD)
            _stage = '<2a> Specificity Scores - BAD'
            if string.join(methods).find('BAD') >= 0:
                ## <i> ## Setup
                self.verbose(1,3,'...BAD',0)
                ## <ii> ## Each residue, common to all BAD
                for r in range(0,seqlen):
                    bad = 0.0
                    badx = 0.0
                    badn = 0.0
                    anc = self.nodeseq[ancnode][r]  # Ancestral AA
                    gx = 0  # Subfam count
                    if anc == 'X':
                        continue
                    for fam in subfam:
                        dup = self.nodeseq[fam][r]
                        sx = 0     # Sequence Count (i.e. all bar Xs)
                        gbad = 0.0
                        if dup == 'X':
                            continue
                        ac = propx - aaprop.pdif['%s%s' % (anc,dup)]    # Anc Conservation of Properties
                        dc = 0.0    # Desc Conservation of Properties
                        for node in clade[fam]:
                            desc = self.nodeseq[node][r]
                            #if desc == 'X':
                            #    continue
                            dc += propx - aaprop.pdif['%s%s' % (dup,desc)]
                            sx += 1
                        if sx > 0:
                            gbad = (dc / sx)
                        if fam == queryfam:
                            #badx = gbad - (propx - aaprop.pdif['%s%s' % (anc,self.nodeseq[outnode][r])])
                            badx = gbad - (propx - aaprop.pdif['%s%s' % (dup,self.nodeseq[outnode][r])])
                        gbad -= ac
                        badn += gbad
                        bad += gbad
                        gx += 1
                    if gx > 0:
                        badn /= (gx - 1)
                    ## <iii> ## Original BAD
                    #self.verbose(1,0,'%d: %f, X %f, N %f' % (r,bad,badx,badn),1)
                    if 'BAD' in methods:
                        self.score['BAD'].append(bad)
                    if 'BADX' in methods:
                        self.score['BADX'].append(badx)
                    if 'BADN' in methods:
                        self.score['BADN'].append(badn)
                        
            ### <b> ### SSC (Livingstone and Barton) & PDAD
            _stage = '<2b> Specificity Scores - SSC & PDAD'
            # .. SS = Total number of conserved properties in subfam clade - ancestral clade
            # .. SSC = Mean across all subfams
            # .. PDAD = Mean PCon in subfams vs mean PCon in all
            if 'SSC' in methods:
                ## <i> ## Setup
                self.verbose(1,3,'...SSC',0)
                pcon = 'PCon_Abs'
                aseqs = []
                fseqs = {}
                for fam in subfam:
                    fseqs[fam] = []
                    for node in clade[fam]:
                        aseqs.append(node.obj['Sequence'])
                        fseqs[fam].append(node.obj['Sequence'])
                seqstat = rje_conseq.SeqStat(log=self.log,aaprop=self.obj['AAProp'],cmd_list=self.cmd_list+['v=-1'])
                seqstat.calcScore(seqlist=aseqs,methods=[pcon])
                ac = seqstat.score[pcon][0:]
                seqstat.score[pcon] = []
                sc = {}
                for fam in subfam:
                    seqstat.calcScore(seqlist=fseqs[fam],methods=[pcon])
                    sc[fam] = seqstat.score[pcon][0:]
                    seqstat.score[pcon] = []
                ## <ii> ## Each residue
                for r in range(seqlen):
                    ssc = 0.0
                    for fam in subfam:
                        ssc += sc[fam][r]
                    ssc = (ssc / len(subfam)) - ac[r]
                    self.score['SSC'].append(ssc)
            if 'PDAD' in methods:
                ## <i> ## Setup
                self.verbose(1,3,'...PDAD',0)
                pcon = 'PCon_Mean'
                seqstat = rje_conseq.SeqStat(log=self.log,aaprop=self.obj['AAProp'],cmd_list=self.cmd_list+['v=-1'])
                sc = {}
                aseqs = []
                fseqs = {}
                # SC
                for fam in subfam:
                    fseqs[fam] = []
                    for node in clade[fam]:
                        aseqs.append(node.obj['Sequence'])
                        fseqs[fam].append(node.obj['Sequence'])
                    seqstat.calcScore(seqlist=fseqs[fam],methods=[pcon])
                    sc[fam] = seqstat.score[pcon][0:]
                    seqstat.score[pcon] = []
                # AC    
                ac = {}
                for fam in subfam:
                    seqstat.calcScore(seqlist=fseqs[fam],methods=[pcon],nonself=aseqs[0:])
                    ac[fam] = seqstat.score[pcon][0:]
                    seqstat.score[pcon] = []
                    
                #seqstat.calcScore(seqlist=aseqs,methods=[pcon])
                #ac = seqstat.score[pcon][0:]
                #seqstat.score[pcon] = []
                ## <ii> ## Each residue
                for r in range(0,seqlen):
                    pdad = 0.0
                    for fam in subfam:
                        pdad += sc[fam][r] - ac[fam][r]
                    pdad /= len(subfam)
                    self.score['PDAD'].append(pdad)
                        
            ### <c> ### ETA (Evolutionary Trace Analysis
            _stage = '<2c> Specificity Scores - ETA'
            # .. (Pure - 1 for all groups, zero for none)
            # .. ETAQ (Quantitative 0-1 stepped by groups)
            if string.join(methods).find('ETA') >= 0:
                self.verbose(1,3,'...ETA',0)
                for r in range(seqlen):
                    eta = 1
                    etaq = 0.0
                    fx = 0
                    raa = []
                    for fam in subfam:
                        et = True
                        aa = 'X'
                        for seq in seqs[fam]:
                            if seq[r] != aa and seq[r] != 'X':
                                if aa != 'X':
                                    et = False
                                    break
                                aa = seq[r]
                        if aa != 'X':
                            fx += 1
                            if et and aa not in raa:
                                raa.append(aa)
                                etaq += 1
                    if fx > 0:
                        etaq /= fx
                    if etaq < 1.0:
                        eta = 0
                    if 'ETA' in methods:
                        self.score['ETA'].append(eta)
                    if 'ETAQ' in methods:
                        self.score['ETAQ'].append(etaq)
            #print self.score.keys()
                    
            ### <3> ### Ranks
            _stage = '<3> Ranks'
            if self.opt['Rank']:
                self.verbose(1,2,'...ranks',0)
            for method in methods:
                self.rank[method] = self.rankScore(self.score[method])

            ### <4> ### Windows
            _stage = '<4> Windows'
            self.verbose(1,2,'...win(%d)' % self.stat['WinSize'],0)
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
                ## <iii> ## Subfams
                self.famwin[method] = {}
                self.famrankwin[method] = {}
                #print subfam
                for fam in subfam:
                    #print method, fam
                    self.famwin[method][fam] = self.winScore(self.score[method],node=fam)
                    self.famrankwin[method][fam] = self.winScore(self.rank[method],node=fam)                   

            self.verbose(0,1,'...Done!',1)
            
        except:
            self.log.errorLog('Major problem during calcScore(%s) %s:' % (_stage,methods))            
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
    funcspec = FuncSpec(log=log,cmd_list=cmd_list,tree=tree,aaprop=aaprop,rank=True)
    for method in methodlist:
        methoddic[method] = copy.deepcopy(funcspec)
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
        print 'grep rje_specificity *.py for other modules that will call conseq.py'

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
