#!/usr/local/bin/python

# rje_haq - Homologue Alignment Quality module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_haq
Description:  Homologue Alignment Quality module
Version:      1.3
Last Edit:    21/02/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    HAQ of HAQESAC: SAQ and PAQ methods from haqesac.

    NB. The classes in this module are designed to take sequence objects and perform analyses, not create sequence objects
    themselves if none are given. The module haqesac.py will do this.

Commandline:
    noquery=T/F : No Query for SAQ, Random Query for PAQ
    saqc=X      : Min no. seqs to share residue in SAQ. [2]
    saqb=X      : SAQ Block length. [10]
    saqm=X      : No. residues to match in SAQ Block. [7]
    saqks=X     : Relative Weighting of keeping Sequences in SAQ. [3]
    saqkl=X     : Relative Weighting of keeping Length in SAQ. [1]
    mansaq=T/F  : Manual over-ride of sequence rejection decisions in SAQ [False]
    paqb=X      : PAQ Block length. [7]
    paqm=X      : No. residues to match in PAQ Block. [3]
    paqks=X     : Relative Weighting of keeping Sequences in PAQ. [3]
    paqkl=X     : Relative Weighting of keeping Length in PAQ. [1]
    manpaq=T/F  : Manual over-ride of sequence rejection decisions in PAQ [False]
    anchors=T/F	: Whether to use conserved 'anchors' to extend well-aligned regions in PAQ	[True]

Uses general modules: copy, os, random, sys, time
Uses RJE modules: rje
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import copy
import os
import random
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
# 0.1 - Added query subsequence focus and optional scoring matrices
# 1.0 - Fully working version for functional HAQESAC
# 1.1 - Added mansaq and manpaq options
# 1.2 - Altered the filling in of anchors to ignore gaps at the end of sequences (e.g. sequence fragments)
# 1.3 - Tidied up some of the code. Increased output during SAQ and PAQ.
#############################################################################################################################
### Major Functionality to Add
# [?] Add altering of info['SAQX'] following realignemnt in PAQ
#############################################################################################################################
default_aa_score_matrix = { 'AA':1.0, 'CC':1.0, 'DD':1.0, 'EE':1.0, 'FF':1.0,
                            'GG':1.0, 'HH':1.0, 'LL':1.0, 'RR':1.0, 'VV':1.0,
                            'MM':1.0, 'II':1.0, 'PP':1.0, 'SS':1.0, 'WW':1.0,
                            'NN':1.0, 'KK':1.0, 'QQ':1.0, 'TT':1.0, 'YY':1.0 }
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#############################################################################################################################

#############################################################################################################################
### HAQ Class: 
#############################################################################################################################
class HAQ(rje.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name
    - HAQMatrix
    
    Opt:boolean
    - NoQuery = No Query for SAQ, Random Query for PAQ
    - Anchors = Whether to use conserved 'anchors' to extend well-aligned regions in PAQ	
    - ManSAQ = Manual over-ride of sequence rejection decisions in SAQ [False]
    - ManPAQ = Manual over-ride of sequence rejection decisions in PAQ [False]
    
    Stat:numeric
    - SAQCyc = SAQ cycle number - used for fulldata output
    - PAQCyc = PAQ cycle number - used for fulldata output
    - SAQCon = Min no. seqs to share residue in SAQ. [2]
    - SAQBlock = SAQ Block length. [10]
    - SAQMatch = No. residues to match in SAQ Block. [7]
    - SAQKeepSeq = Relative Weighting of keeping Sequences in SAQ. [3]
    - SAQKeepLen = Relative Weighting of keeping Length in SAQ. [1]
    - PAQBlock = PAQ Block length. [7]
    - PAQMatch = No. residues to match in PAQ Block. [3]
    - PAQKeepSeq = Relative Weighting of keeping Sequences in PAQ. [3]
    - PAQKeepLen = Relative Weighting of keeping Length in PAQ. [1]

    List:list

    Dict:dictionary
    - AA Score Matrix = Dictionary of scores for amino acid combos {combo:score}

    Obj:RJE_Objects
    '''
    ### Attributes
#############################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                      #
#############################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics 
        self.infolist = ['Name','HAQMatrix']
        self.statlist = ['SAQCyc','PAQCyc','SAQCon','SAQBlock','SAQMatch','SAQKeepSeq','SAQKeepLen','PAQBlock','PAQMatch','PAQKeepSeq','PAQKeepLen']
        self.optlist = ['NoQuery','Anchors','ManSAQ','ManPAQ']
        self.listlist = []
        self.dictlist = ['AA Score Matrix']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'SAQCon':2,'SAQBlock':10,'SAQMatch':7,'SAQKeepSeq':3.0,'SAQKeepLen':1.0,
                      'PAQBlock':7,'PAQMatch':3,'PAQKeepSeq':3.0,'PAQKeepLen':1.0,'PAQCyc':0})
        self.setOpt({'Anchors':True})
        self.dict['AA Score Matrix'] = default_aa_score_matrix
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
                self._cmdRead(cmd,type='opt',att='NoQuery')
                self._cmdRead(cmd,type='int',att='SAQCon',arg='saqc')
                self._cmdRead(cmd,type='int',att='SAQBlock',arg='saqb')
                self._cmdRead(cmd,type='int',att='SAQMatch',arg='saqm')
                self._cmdRead(cmd,type='stat',att='SAQKeepSeq',arg='saqks')
                self._cmdRead(cmd,type='stat',att='SAQKeepLen',arg='saqkl')
                self._cmdRead(cmd,type='int',att='PAQBlock',arg='paqb')
                self._cmdRead(cmd,type='int',att='PAQMatch',arg='paqm')
                self._cmdRead(cmd,type='stat',att='PAQKeepSeq',arg='paqks')
                self._cmdRead(cmd,type='stat',att='PAQKeepLen',arg='paqkl')
                self._cmdRead(cmd,type='info',att='HAQMatrix')
                self._cmdRead(cmd,type='opt',att='Anchors')
                self._cmdRead(cmd,type='opt',att='ManSAQ')
                self._cmdRead(cmd,type='opt',att='ManPAQ')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.info['HAQMatrix'] != 'None':
            self.loadMatrix()
#############################################################################################################################
    def loadMatrix(self,filename=None):     ### Builds HAQ Matrix from file
        '''
        Builds HAQ Matrix from file
        >> filename:str = Name of file (if None, use self.info['HAQMatrix'])
        '''
        try:
            _stage = '<0> Setup'
            if not filename:
                filename = self.info['HAQMatrix']
            else:
                self.info['HAQMatrix'] = filename
            self.verbose(0,3,'Loading HAQ matrix...',0)

            _stage = '<1> Open file & Read Lines'
            filematrix = self.loadFromFile(filename)
            if len(filematrix) < 1:
                self.log.errorLog('Loading HAQ Matrix (%s) failed.' % filename,printerror=False,quitchoice=True)
                return
            
            _stage = '<2a> Read in alphabet'
            self.verbose(1,3,filematrix[0],1)
            _alphabet = filematrix[0].split()
            _stage = '<2b> Read in Matrix'
            self.dict['AA Score Matrix'] = {}
            for a in _alphabet:
                line = filematrix[_alphabet.index(a)+1].split()
                if len(line) != (len(self.alphabet)+1):
                    self.log.errorLog('%s has wrong format! Does not match %s' % (line, _alphabet))
                    self.dict['AA Score Matrix'] = default_aa_score_matrix
                    raise ValueError
                for i in range(1,len(line)):
                    haqscore = float(line[i])
                    if haqscore:
                        self.dict['AA Score Matrix']['%s%s' % (a,_alphabet[i-1])] = haqscore
            self.verbose(0,2,'Done!',2)
        except:
            self.log.errorLog('Oh dear. Error in loadMatrix(). Will use default matrix.',quitchoice=True)
#############################################################################################################################
    ### <2> ### HAQ Methods                                                                                                 #
#############################################################################################################################
    def _saqCon(self, qres, ores):    ### Returns score of qres vs ores 
        '''
        Returns score of qres vs ores 
        >> qres:str = amino acid (or nt?)
        >> ores:str = amino acid (or nt?)
        << Float
        '''
        if qres == 'X':
            return 0.0
        for _comp in ['%s%s' % (qres,ores),'%s%s' % (ores,qres)]:
            if _comp in self.dict['AA Score Matrix'].keys():
                return self.dict['AA Score Matrix'][_comp]
        return 0.0
#############################################################################################################################
    def saveHAQ(self,seqlist,filename,key='SAQ'):     ### Saves appropriate Key of seq.info in Fasta File
        '''
        Saves using seq.info[key] rather than seq.info['Sequence'].
        >> seqlist:rje_seq.SeqList Object
        >> key:str = key for seq.info 
            - SAQX = SAQ sequences with individual Xs
            - SAQ = SAQ sequences with aligment Xs
            - PAQ = PAQ sequences with aligment Xs
        '''
        try:
            haqlist = seqlist   # SeqList Object to store individually Xd sequences
            for seq in haqlist.seq:
                (seq.info['Sequence'],seq.info[key]) = (seq.info[key],seq.info['Sequence'])
            haqlist.saveFasta(seqfile=filename)
            for seq in haqlist.seq:
                (seq.info[key],seq.info['Sequence']) = (seq.info['Sequence'],seq.info[key])
        except:
            self.log.errorLog('Problem with saveHAQ(%s)' % key)
#############################################################################################################################
    def singleSeqAQ(self,seqlist,focus=[0,-1]):     ### Performs SAQ on seqlist, adding seq.info['SAQ']
        '''
        Performs SAQ on seqlist, adding seq.info['SAQ'].
        >> seqlist:rje_seq.SeqList Object
        - NB. This object will itself have sequences removed from it, so beware!
        - A new info key will be added: SAQX = SAQ sequences with individual Xs
        - A new info key will be added: SAQ = SAQ sequences with aligment Xs
        >> focus:list of range positions [X:Y] to look at. If Y=0 then [X:].
        '''
        ### <SAQ1> ### Setup
        try:
            _stage = '<1> Setup'
            haqlist = seqlist   # SeqList Object to store individually Xd sequences
            query = haqlist.obj['QuerySeq']
            if self.opt['NoQuery']:
                query = None
            badres = [-1,0]     # List of how many bad residues in total dataset
            block_align = {}    # Dictionary of whether residue in block of sequence that is well-aligned or not
            res_align = {}      # Dictionary of whether residue of sequence is well-aligned or not
            res_gap = {}        # Dictionary of whether residue of sequence is a gap or not
            gap_align = {}      # Dictionary of whether residue of sequence is a gap in a well-aligned block or not
            for seq in haqlist.seq:
                seq.info['SAQ'] = seq.info['Sequence'][0:]      # Note! Sequence is modified and SAQ not, then they are swapped at end!
                block_align[seq] = [False] * seq.seqLen()
                res_align[seq] = [False] * seq.seqLen()
                res_gap[seq] = [False] * seq.seqLen()
                gap_align[seq] = [False] * seq.seqLen()

        ### <SAQ2> ### Repeated cycles of defining well- and badly-aligned blocks
            #X#self.deBug(self.stat)
            _stage = '<2> BlockID'
            while badres[-1] != badres[-2]:     # Change in number of bad residues
                total_res = 0
                badres.append(0)    # badres[-1] is the current number of bad residues
                infotxt = 'SAQ%d-%d: Calculating "bad" residues ...' % (self.stat['SAQCyc'],len(badres)-2)
                for seq in haqlist.seq:
                    myinfo = '%s %.1f%%' % (infotxt,(100.0 * haqlist.seq.index(seq) / haqlist.seqNum()))
                    self.log.printLog('\r#SAQ',myinfo,log=False,newline=False)
                    #self.verbose(0,3,'\r%45s' % myinfo,0)

                    ## <SAQ2a> ## For each sequence, mark residues as aligned or gapped
                    _stage = '<2a> Mark Residues'
                    for r in range(seq.seqLen()):
                        gap_align[seq][r] = False
                        res_align[seq][r] = False
                        if block_align[seq][r] or len(badres) == 3:     # After first cycle, look only at well-aligned blocks (well-aligned for sequence not whole alignment)
                            a = seq.info['Sequence'][r]
                            res_gap[seq][r] = False
                            if a == '-':
                                res_gap[seq][r] = True
                                gap_align[seq][r] = True
                            else:   # 'X' handled by self._saqCon
                                conx = 0  # Matches with good regions of otherseqs (*including self*)
                                for otherseq in haqlist.seq[0:]:
                                    if otherseq == seq:     # > so self not counted!
                                        continue
                                    if len(otherseq.info['Sequence']) != len(seq.info['Sequence']):
                                        self.log.errorLog('Sequence lengths do not match - should be aligned!',printerror=False)
                                        raise ValueError
                                    if (block_align[otherseq][r] or len(badres) == 3):
                                        conx += self._saqCon(a, otherseq.info['Sequence'][r])
                                #if seq == query and r > 590:
                                #    print seq.shortName(),r,conx,'vs',self.stat['SAQCon'],
                                if conx >= self.stat['SAQCon']:    
                                    res_align[seq][r] = True
                        #if seq == query and r > 590:
                        #    print r, res_align[seq][r]

                    ## <SAQ2b> ## Marked regions of well-aligned residues for each sequence
                    _stage = '<2b> Mark Regions'
                    ## <i> ## Clear first
                    _stage = '<2b-i> Mark Regions'
                    for r in range(seq.seqLen()):
                        block_align[seq][r] = False
                    ## <ii> ## Recalculate
                    _stage = '<2b-ii> Mark Regions'
                    for r in range(seq.seqLen()):
                        _stage = '<2b-ii> Blocks'
                        if res_align[seq][r]:   # Start of potential block
                            blen = 0    # Block length (SAQBlock) = AAs
                            win = 0     # Window length = all sequence
                            matchx = 1  # Good residues in window (first residue must be good!) (SAQMatch)
                            while blen < self.stat['SAQBlock'] and matchx < self.stat['SAQMatch']:
                                win += 1
                                if (r + win) >= seq.seqLen() or seq.info['Sequence'][r+win] == 'X':     # Hit Bad Region: Abort
                                    break
                                else:   # Better region
                                    if gap_align[seq][r+win]:   # Decent gap
                                        continue
                                    else:
                                        blen += 1   # Increase Block
                                        if res_align[seq][r+win]:   # Good residue
                                            matchx += 1
                            #if seq == query and r > 590:
                            #    print seq.shortName(),r,matchx,'vs',self.stat['SAQMatch'],
                            if matchx >= self.stat['SAQMatch']:
                                for w in range((win+1)):
                                    block_align[seq][r+w] = True
                        #if seq == query and r > 590:
                        #    print r, block_align[seq][r]
                    ## <iii> ## Update bad residue count
                    for r in range(seq.seqLen()):
                        _stage = '<2b-iii> Mark Regions'
                        #print seq.shortName(), r, seq.seqLen(), block_align[seq][r], res_gap[seq][r], badres[-1]   # Bad residue
                        if not block_align[seq][r] and not res_gap[seq][r]:   # Bad residue
                            badres[-1] += 1
                        if not res_gap[seq][r]:
                            total_res += 1
                myinfo = '%s 100.0%%' % infotxt
                myinfo += ' => %s bad of %s total residues' % (rje.integerString(badres[-1]),rje.integerString(total_res))
                self.log.printLog('\r#SAQ',myinfo)
                #self.verbose(0,3,'\r%45s' % myinfo,0)
                if badres[-1] == total_res:
                    self.log.errorLog('All residues marked as bad in SAQ!',printerror=False,quitchoice=True)
                # Now have all residues in all sequences marked as good (block_align=True) or bad (block_align=False)

        ### <SAQ3> ### X out badly-aligned blocks
            _stage = '<3> X-Out'
            self.log.printLog('#SAQ','SAQ%d-%d: Masking "bad" residues ...' % (self.stat['SAQCyc'],len(badres)-2),log=False,newline=False)
            #self.verbose(0,3,'SAQ%d-%d: Masking "bad" residues ...' % (self.stat['SAQCyc'],len(badres)-2),0)
            for seq in haqlist.seq:
                newseq = ''
                for r in range(seq.seqLen()):
                    if block_align[seq][r] or seq.info['Sequence'][r] == '-':   #!# Was backwards? res_gap[seq][r] == False:
                        newseq += seq.info['Sequence'][r]
                    else: # Bad residue
                        newseq += 'X'
                seq.info['Sequence'] = newseq[0:]
                seq.info['SAQX'] = newseq[0:]       # Stores Xd sequences for individuals for use in PAQ
            #!# Add saving of data in 'datafull' option

        ### <SAQ4> ### Remove sequences and/or badly-aligned regions
            _stage = '<4> Removal'
            self.log.printLog('\r#SAQ','SAQ%d-%d: Removing bad sequences and/or dodgy regions...' % (self.stat['SAQCyc'],len(badres)-2),log=False,newline=False)
            #self.verbose(0,3,'\rSAQ%d-%d: Removing bad sequences and/or dodgy regions...' % (self.stat['SAQCyc'],len(badres)-2),0)
            ## <SAQ4a> ## Process Query first - only interested in good regions within query
            _stage = '<4a> Query Removal'
            if self.opt['NoQuery'] or query == None:  # No preprocessing of Query
                self.verbose(0,4,'no Master Query processing...',0)
            else:
                haqlist.mapX(query, qtrim=True, focus=focus) # Replaces other sequence ends and query X columns with Xs
                self.verbose(0,4,'Query (%s) processed...' % query.shortName(),0)
            self.verbose(0,3,'',1)
            if self.opt['ManSAQ']:
                haqlist.saveFasta(seqfile='%s.mansaq.fas' % haqlist.info['Basefile'])

            ## <SAQ4b> ## Cycle through other sequences (worst first) until no more good residues or sequences are lost
            _stage = '<4b> Seq Removal'
            goodres = [0, self._getGood(haqlist.seq)]   # List of number of 'good' residues
            goodseq = [0, haqlist.seqNum()]
            while goodres[-1] != goodres[-2] or goodseq[-1] != goodseq[-2]:
                colgood = [0] * haqlist.seq[0].seqLen()    # Good residues per column
                for r in range(haqlist.seq[0].seqLen()):
                    for seq in haqlist.seq:
                        if seq.info['Sequence'][r] != '-' and seq.info['Sequence'][r] != 'X':
                            colgood[r] += 1
                ## <i> ## Compare relative loss of masking and losing each sequence
                keepx = {}  # Dictionary of seq:number of lost residues if seq kept
                losex = {}  # Dictionary of seq:number of lost residues if seq lost
                badkx = -1  # Biggest loss if kept
                badlx = -1  # Biggest loss if lost
                bads = None # Worst sequence
                for seq in haqlist.seq:
                    if seq == query and self.opt['NoQuery'] == False:
                        continue    # Next sequence
                    # Calculate keepx and losex
                    keepx[seq] = 0
                    for r in range(seq.seqLen()):
                        if seq.info['Sequence'][r] == 'X':
                            keepx[seq] += colgood[r]
                    losex[seq] = self._getGood([seq])
                    # Update bads if worse
                    if keepx[seq] > badkx:
                        badkx = keepx[seq]
                        badlx = losex[seq]
                        bads = seq
                    elif keepx[seq] == badkx and losex[seq] < badlx:
                        badlx = losex[seq]
                        bads = seq
                ## <ii> ## Remove bad sequences and/or regions
                if badkx > 0:
                    if self.opt['ManSAQ']:
                        default = 'N'
                        if badkx * self.stat['SAQKeepLen'] > badlx * self.stat['SAQKeepSeq']:   # Lose sequence!
                            default = 'Y'
                        if rje.yesNo('%s worst: -%s aa if kept vs -%s aa if lost. Remove?' % (bads.shortName(),rje.integerString(badkx),rje.integerString(badlx)),default):
                            seqlist.removeSeq(text='SAQ%d: -%s aa if kept vs -%s aa if lost. (Manual decision.)' % (self.stat['SAQCyc'],rje.integerString(badkx),rje.integerString(badlx)),seq=bads)
                        else:   # X out
                            haqlist.mapX(bads)
                    else:
                        self.verbose(1,3,'%s worst: -%s aa if kept vs -%s aa if lost.' % (bads.shortName(),rje.integerString(badkx),rje.integerString(badlx)),1)
                        #!# Add option for upweighting certain sequence type? (e.g. vs fragment or hypothetical?)
                        if badkx * self.stat['SAQKeepLen'] > badlx * self.stat['SAQKeepSeq']:   # Lose sequence!
                            haqlist.removeSeq(text='SAQ%d: -%s aa if kept vs -%s aa if lost.' % (self.stat['SAQCyc'],rje.integerString(badkx),rje.integerString(badlx)),seq=bads)
                        else:   # X out
                            haqlist.mapX(bads)
                ### <iii> ### Recalculate goodres
                goodres.append(self._getGood(haqlist.seq))
                goodseq.append(haqlist.seqNum())
                #X#self.verbose(1,3,'%d -> %d "good" aa' % (goodres[-2],goodres[-1]),1)

            ### <SAQ5> ### Reinstate UnX'd sequence:
            _stage = '<4b> Seq Removal'
            for seq in haqlist.seq:
                #print seq.info
                [seq.info['SAQ'],seq.info['Sequence']] = [seq.info['Sequence'],seq.info['SAQ']]
            if self.opt['ManSAQ'] and rje.checkForFile('%s.mansaq.fas' % haqlist.info['Basefile']):
                os.unlink('%s.mansaq.fas' % haqlist.info['Basefile'])

        except:
            self.log.errorLog('Problem with singleSeqAQ() %s.' % _stage, quitchoice=True)
#############################################################################################################################
    def pairwiseAQ(self,seqlist=None,query=None,focus=[0,0]):     ### Performs PAQ on seqlist, adding seq.info['PAQ']
        '''
        Performs PAQ on seqlist, adding seq.info['PAQ']
        >> seqlist:rje_seq.SeqList Object
        - NB. This object will itself have sequences removed from it, so beware!
        - A new info key will be added: PAQ = PAQ sequences with alignment Xs
        >> focus:list of range positions [X:Y] to look at. If Y=0 then [X:]. 
        '''
        ### <PAQ0> ### Setup
        try:
            _stage = '<0> Setup'
            haqlist = seqlist   # SeqList Object to store individually Xd sequences
            if not query:
                query = haqlist.obj['QuerySeq']
            if self.opt['NoQuery'] or not query:
                query = haqlist.seq[random.randint(0,haqlist.seqNum()-1)]
                self.log.printLog('#QRY','Temp (random) query %s assigned for PAQ' % query.shortName())
            #!# paqx = [False] * seqlist.seq[0].seqLen()    # List of whether a column of the alignment is bad (has an X) [True] or not [False]
            #!# - make this a method?!

            pwaq = {}    # Dictionary of lists of pairwise alignements
            block_align = {}    # Dictionary of whether residue in block of sequence that is well-aligned or not
            for seq in haqlist.seq:
                block_align[seq] = [False] * seq.seqLen()
                seq.info['PAQ'] = seq.info['Sequence'][0:]
                if seq.info.has_key('SAQX') and len(seq.info['SAQX']) == seq.seqLen():   #!# Should no longer be issues due to length changes following realignment
                    seq.info['Sequence'] = seq.info['SAQX'][0:]
                elif seq.info.has_key('SAQX'):
                    self.log.errorLog('Cannot use SAQX for %s in PAQ as wrong length.' % seq.shortName(),printerror=False)
                for otherseq in haqlist.seq:
                    pwaq[(seq,otherseq)] = [False] * seq.seqLen()

        ### <PAQ1> ### Directional Pairwise Comparisons of sequences
            _stage = '<1> Pairwise Comparisons'
            infotxt = 'PAQ%d: Pairwise Comparisons ...' % self.stat['PAQCyc']
            #print self.stat
            for seq in haqlist.seq:
                for otherseq in haqlist.seq:
                    myinfo = '%s %.1f%% %.1f%%   ' % (infotxt,(100.0 * haqlist.seq.index(seq) / haqlist.seqNum()),(100.0 * haqlist.seq.index(otherseq) / haqlist.seqNum()))
                    self.log.printLog('\r#PAQ',myinfo,log=False,newline=False)
                    for r in range(seq.seqLen()):
                        ar = seq.info['Sequence'][r]
                        ## <i> ## Look for PW aligned block
                        _stage = '<1-i> Pairwise Comparisons'
                        if ar not in ['-','X']: # Start of test block
                            blen = 0    # Block length (PAQBlock) = AAs
                            win = 0     # Window length = all sequence
                            matchx = 0  # Score for residues in window 
                            while blen < self.stat['PAQBlock'] and (r+win) < seq.seqLen():     # This time we allow overshoots in both directions
                                ar = seq.info['Sequence'][r+win]
                                at = otherseq.info['Sequence'][r+win]
                                if 'X' in [ar,at]:     # Hit Bad Region: Abort
                                    break
                                else:   # Better region
                                    if ar != '-':   
                                        blen += 1   # Increase Block
                                        matchx += self._saqCon(ar,at)
                                win += 1
                        ## <ii> ## Update pwaq if block good
                            _stage = '<1-ii> Pairwise Comparisons'
                            if matchx >= self.stat['PAQMatch']:
                                for w in range(win):
                                    if seq.info['Sequence'][r+w] in ['-','X']:
                                        pwaq[(seq,otherseq)][r+w] = False
                                    else:
                                        pwaq[(seq,otherseq)][r+w] = True           
            self.log.printLog('\r#PAQ','%s 100.0% 100.0%.   ' % infotxt,log=False)
                
        ### <PAQ2> ### Link back to Query
            _stage = '<2> Linking to Query'
            ### <PAQ2a> ### Network of Pairwise Quality alignments
            _stage = '<2a> Linking to Query'
            #self.verbose(1,3,'PAQ%d: Linking Residues to Query (%s)' % (self.stat['PAQCyc'],query.shortName()),0)
            infotxt = 'PAQ%d: Linking Residues to Query (%s) ...' % (self.stat['PAQCyc'],query.shortName())
            for r in range(query.seqLen()):
                _stage = '<2a> Linking to Query'
                self.log.printLog('\r#PAQ','%s %.1f%%' % (infotxt,(100.0 * r / query.seqLen())),log=False,newline=False)
                qok = {}    # Dictionary of whether residue in seq OK, i.e. linked to query
                for seq in haqlist.seq:
                    qok[seq] = False
                qok[query] = True
                sok = [0,1] # List of OK sequence for residue
                while sok[-2] != sok[-1]:
                    ## <i> ## Match pairs, starting with query
                    _stage = '<2a-i> Linking to Query'
                    for seq in haqlist.seq:
                        if qok[seq]:
                            for otherseq in haqlist.seq:
                                if pwaq[(seq,otherseq)][r] or pwaq[(otherseq,seq)][r]:
                                    qok[otherseq] = True
                    ## <ii> ## Update sok
                    _stage = '<2a-ii> Linking to Query'
                    sok.append(0)
                    for seq in haqlist.seq:
                        if qok[seq]:
                            sok[-1] += 1
                            block_align[seq][r] = True
                _stage = '<2a-iii> Linking to Query'
                if sok[-1] == 1:    # Only query OK!
                    block_align[query][r] = False
            self.log.printLog('\r#PAQ','%s 100.0%%' % infotxt,log=False)
            
            ### <PAQ2b> ### Allow for divergence (Conserved Anchors)
            _stage = '<2b> Anchors'
            if self.opt['Anchors']:
                infotxt = 'PAQ%d: Accounting for divergence within aligned regions ...' % self.stat['PAQCyc']
                ## <i> ## Setup gapped list
                gapped = [False] * query.seqLen()   # Whether column of alignment is gapped
                for seq in haqlist.seq:
                    self.log.printLog('\r#PAQ','%s %.1f%%  ' % (infotxt,(50.0 * haqlist.seq.index(seq) / haqlist.seqNum())),log=False,newline=False)
                    (start,end) = (0,seq.seqLen())
                    while seq.info['Sequence'][start] == '-':
                        start += 1
                    while seq.info['Sequence'][end-1] == '-':
                        end -=1
                    for r in range(start,end):
                        if seq.info['Sequence'][r] == '-':
                            gapped[r] = True
                ## <ii> ## Correction
                for seq in haqlist.seq:
                    self.log.printLog('\r#PAQ','%s %.1f%%  ' % (infotxt,(50 + (50.0 * haqlist.seq.index(seq) / haqlist.seqNum()))),log=False,newline=False)
                    for r in range(seq.seqLen()):
                        if block_align[seq][r] or gapped[r]:    # No need for correction
                            continue
                        # Move in both directions: if good residues (or sequence end) reached before gaps then reinstate
                        winf = 0
                        fwd = True
                        fok = False
                        winb = 0
                        bwd = True
                        bok = False
                        while fwd or bwd:
                            # End of seqs
                            if (r + winf) >= seq.seqLen():
                                fwd = False
                            if (r - winb) < 0:
                                bwd = False
                            # Gaps/OK
                            if fwd:
                                if gapped[r+winf]:
                                    fok = False
                                    fwd = False
                                elif block_align[seq][r+winf]:
                                    fwd = False
                                else:
                                    winf += 1
                            if bwd:
                                if gapped[r-winb]:
                                    bok = False
                                    bwd = False
                                elif block_align[seq][r-winb]:
                                    bwd = False
                                else:
                                    winb += 1
                        if fok and bok: # Reinstate
                            for w in range(r-winb,r+winf+1):
                                block_align[seq][w] = True
                self.log.printLog('\r#PAQ','%s 100.0%%  ' % infotxt,log=False)

        ### <PAQ3> ### X out badly-aligned blocks
            _stage = '<3> Making bad sequence blocks'
            for seq in haqlist.seq:
                newseq = ''
                for r in range(seq.seqLen()):
                    if block_align[seq][r] or seq.info['Sequence'][r] == '-':
                        newseq += seq.info['Sequence'][r]
                    else: # Bad residue
                        newseq += 'X'
                seq.info['Sequence'] = newseq[0:]
            #!# Add saving of data in 'datafull' option

        ### <PAQ4> ### Remove sequences and/or badly-aligned regions
            _stage = '<4> Removing sequences/regions'
            self.verbose(0,4,'PAQ%d: Removing bad sequences and/or dodgy regions...' % self.stat['PAQCyc'],0)
            ## <PAQ4a> ## Process Query first - only interested in good regions within query
            if self.opt['NoQuery']:  # No preprocessing of Query
                self.verbose(0,4,'no Master Query processing...',0)
            else:
                haqlist.mapX(query, qtrim=True, focus=focus) # Replaces other sequence ends and query X columns with Xs
                self.verbose(0,4,'Query (%s) processed...' % query.shortName(),0)
            self.verbose(0,3,'',1)
            if self.opt['ManPAQ']:
                haqlist.saveFasta(seqfile='%s.manpaq.fas' % haqlist.info['Basefile'])

            ## <PAQ4b> ## Cycle through other sequences (worst first) until no more good residues are lost
            goodres = [0, self._getGood(haqlist.seq)]   # List of number of 'good' residues
            goodseq = [0, haqlist.seqNum()]
            while goodres[-1] != goodres[-2] or goodseq[-1] != goodseq[-2]:
                colgood = [0] * haqlist.seq[0].seqLen()    # Good residues per column
                for r in range(haqlist.seq[0].seqLen()):
                    for seq in haqlist.seq:
                        if seq.info['Sequence'][r] != '-' and seq.info['Sequence'][r] != 'X':
                            colgood[r] += 1
                ## <i> ## Compare relative loss of masking and losing each sequence
                keepx = {}  # Dictionary of seq:number of lost residues if seq kept
                losex = {}  # Dictionary of seq:number of lost residues if seq lost
                badkx = -1  # Biggest loss if kept
                badlx = -1  # Biggest loss if lost
                bads = None # Worst sequence
                for seq in haqlist.seq:
                    if seq == query and self.opt['NoQuery'] == False:
                        continue    # Next sequence
                    # Calculate keepx and losex
                    keepx[seq] = 0
                    for r in range(seq.seqLen()):
                        if seq.info['Sequence'][r] == 'X':
                            keepx[seq] += colgood[r]
                        #?# In Perl HAQESAC there was an option to ignore Orphans in this calculation. Reinstate?
                    losex[seq] = self._getGood([seq])
                    # Update bads if worse
                    if keepx[seq] > badkx:
                        badkx = keepx[seq]
                        badlx = losex[seq]
                        bads = seq
                    elif keepx[seq] == badkx and losex[seq] < badlx:
                        badlx = losex[seq]
                        bads = seq
                ## <ii> ## Remove bad sequences and/or regions
                if badkx > 0:
                    if self.opt['ManPAQ']:
                        default = 'N'
                        if badkx * self.stat['PAQKeepLen'] > badlx * self.stat['PAQKeepSeq']:   # Lose sequence!
                            default = 'Y'
                        if rje.yesNo('%s worst: -%s aa if kept vs -%s aa if lost. Remove?' % (bads.shortName(),rje.integerString(badkx),rje.integerString(badlx)),default):
                            seqlist.removeSeq(text='PAQ%d: -%s aa if kept vs -%s aa if lost. (Manual decision.)' % (self.stat['PAQCyc'],rje.integerString(badkx),rje.integerString(badlx)),seq=bads)
                        else:   # X out
                            haqlist.mapX(bads)
                    else:
                        self.verbose(1,3,'%s worst: -%s aa if kept vs -%s aa if lost.' % (bads.shortName(),rje.integerString(badkx),rje.integerString(badlx)),1)
                        #!# Add option for upweighting certain sequence type? (e.g. vs fragment or hypothetical?)
                        if badkx * self.stat['PAQKeepLen'] > badlx * self.stat['PAQKeepSeq']:   # Lose sequence!
                            seqlist.removeSeq(text='PAQ%d: -%s aa if kept vs -%s aa if lost.' % (self.stat['PAQCyc'],rje.integerString(badkx),rje.integerString(badlx)),seq=bads)
                        else:   # X out
                            haqlist.mapX(bads)
                ### <iii> ### Recalculate goodres
                goodres.append(self._getGood(haqlist.seq))
                goodseq.append(haqlist.seqNum())
                self.verbose(1,3,'%d -> %d "good" aa' % (goodres[-2],goodres[-1]),1)
                        
        ### <PAQ5> ### Reinstate UnX'd sequence:
            _stage = '<5> Replacing sequences'
            for seq in haqlist.seq:
                [seq.info['PAQ'],seq.info['Sequence']] = [seq.info['Sequence'],seq.info['PAQ']]
            if self.opt['ManPAQ'] and rje.checkForFile('%s.manpaq.fas' % haqlist.info['Basefile']):
                os.unlink('%s.manpaq.fas' % haqlist.info['Basefile'])

        except:
            self.log.errorLog('rje_haq.py ~ Problem with pairwiseAQ %s.' % _stage, True)
#############################################################################################################################
    def _getGood(self,seqlist):     ### Returns number of good (non-X) characters in sequence list 
        '''
        Returns number of good (non-X) characters in sequence list (relative to query).
        >> seqlist:list of Sequence Objects
        '''
        nonx = 0
        for seq in seqlist:
            nonx += seq.aaLen() - seq.info['Sequence'].count('X')
        return nonx
#############################################################################################################################
### End of HAQ Class                                                                                                        #
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                               #
#############################################################################################################################
if __name__ == "__main__":      
    try:
        print 'Not for standalone running.'
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#############################################################################################################################
### END OF SECTION III                                                                                                      #
#############################################################################################################################
