#!/usr/bin/python

# See below for name and description
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_slimcalc
Description:  SLiM Attribute Calculation Module
Version:      0.10.0
Last Edit:    20/03/18
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is based on the old rje_motifstats module. It is primarily for calculating empirical attributes of SLiMs
    and their occurrences, such as Conservation, Hydropathy and Disorder. 

Commandline:
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Motif Occurrence Attribute Options ###
    slimcalc=LIST   : List of additional attributes to calculate for occurrences - Cons,SA,Hyd,Fold,IUP,Chg,Comp []
    winsize=X       : Used to define flanking regions for calculations. If negative, will use flanks *only* [0]
    relconwin=X     : Window size for relative conservation scoring [30]
    iupath=PATH     : The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    iucut=X         : Cut-off for IUPred results (0.0 will report mean IUPred score) [0.0]
    iumethod=X      : IUPred method to use (long/short) [short]
    percentile=X    : Percentile steps to return in addition to mean [0]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Alignment Settings ###
    usealn=T/F      : Whether to search for and use alignments where present. [False]
    alndir=PATH     : Path to pre-made alignment files [./]
    alnext=X        : File extension of alignment files, AccNum.X (checked before Gopher used) [orthaln.fas]
    gopherdir=PATH  : Path from which to look for GOPHER alignments (if not found in alndir) and/or run GOPHER [./] 
    usegopher=T/F   : Use GOPHER to generate orthologue alignments missing from alndir - see gopher.py options [False]
    fullforce=T/F   : Whether to force regeneration of alignments using GOPHER [False]
    orthdb=FILE     : File to use as source of orthologues for GOPHER []
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Conservation Parameters ###   
    conspec=LIST    : List of species codes for conservation analysis. Can be name of file containing list. [None]
    conscore=X      : Type of conservation score used:  [rlc]
                        - abs = absolute conservation of motif using RegExp over matched region
                        - pos = positional conservation: each position treated independently
                        - prob = conservation based on probability from background distribution
                        - prop = conservation of amino acid properties
                        - rlc = relative local conservation
                        - all = all methods for comparison purposes
    consamb=T/F     : Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    consinfo=T/F    : Weight positions by information content (does nothing for conscore=abs) [True]
    consweight=X    : Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
                        - 0 gives equal weighting to all. Negative values will upweight distant sequences.
    minhom=X        : Minimum number of homologues for making conservation score [1]
    homfilter=T/F   : Whether to filter homologues using seqfilter options [False]
    alngap=T/F      : Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore
                      as putative sequence fragments [False]  (NB. All X regions are ignored as sequence errors.)
    posmatrix=FILE  : Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) [None]
    aaprop=FILE     : Amino Acid property matrix file. [aaprop.txt]
    masking=T/F     : Whether to use seq.info['MaskSeq'] for Prob cons, if present (else 'Sequence') [True]
    vnematrix=FILE  : BLOSUM matrix file to use for VNE relative conservation []
    relgappen=T/F   : Whether to invoke the "Gap Penalty" during relative conservation calculations [True]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SLiM/Occ Filtering Options ###
    slimfilter=LIST : List of stats to filter (remove matching) SLiMs on, consisting of X*Y  []
                      - X is an output stat (the column header),
                      - * is an operator in the list >, >=, !=, =, >= ,<    
                      - Y is a value that X must have, assessed using *.
                      This filtering is crude and may behave strangely if X is not a numerical stat!
                      !!! Remember to enclose in "quotes" for <> filtering !!!
    occfilter=LIST  : Same as slimfilter but for individual occurrences []
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, math, pickle, os, re, string, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
import rje_aaprop, rje_disorder, rje_seq, rje_sequence, rje_slim, rje_scoring
import gopher
import ned_eigenvalues
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on rje_motifstats methods.
    # 0.1 - Added new probability-based method, inspired by (but different to) Dinkel & Sticht (2007)
    # 0.2 - Mended OccPos finding for wildcards. Added new relative conservation score.
    # 0.3 - Added von Neumann entropy code.
    # 0.4 - Added Webserver pickling of RLC lists.
    # 0.5 - Altered to use GOPHER V3 and handle nested alignment directories.
    # 0.6 - Minor tweak to avoid unwanted GOPHER directory generation.
    # 0.7 - Added RLC to "All" conscore running.
    # 0.8 - Made RLC the default.
    # 0.9 - Improvements to use of GOPHER.
    # 0.9.1 - Modified combining of motif stats to handle expectString format for individual values.
    # 0.9.2 - Changed default conscore in docstring to RLC.
    # 0.9.3 - Changed fudge error to warning.
    # 0.10.0 - Added extra disorder methods to slimcalc.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add and fully integrate the new/custom score options.
    # [ ] : Add new probabilistic Conservation score.
    # [ ] : Expand probabilistic Conservation score to other slimcalc scores.
    # [ ] : Split into rje_seqcons, which should then be inherited by this module but usable by rje_seq alone.
    # [ ] : Add tree-weighted conservation. Generally tidy code.
    # [ ] : Add ANCHOR to slimcalc score output.
    # [Y] : Add general disorder output "Dis" to scores to calculate and additional disorder scores
    '''
#########################################################################################################################
### CONSTANTS ###                                                                                                     
occstats = {'con':'Cons','cons':'Cons','sa':'SA','hyd':'Hyd','fold':'Fold','iup':'IUP','dis':'Dis','comp':'Comp',
            'chg':'Chg','charge':'Chg'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMCalc Class                                                                                          #
#########################################################################################################################
class SLiMCalc(rje.RJE_Object):     
    '''
    SLiMCalc Class. Author: Rich Edwards (2007).

    Info:str
    - AlnDir = Path to alignment files [./]
    - AlnExt = File extensions of alignments: AccNum.X [aln.fas]
    - ConScore = Type of conservation score used:  [abs]
        - abs = absolute conservation of motif: reports percentage of homologues in which conserved
        - prop = conservation of amino acid properties
        - prob = conservation based on probability from background distribution
        - pos = positional conservation: each position treated independently 
        - all = all three methods for comparison purposes
    - GopherDir = Directory from which GOPHER will be run []
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix)
    - VNEMatrix = BLOSUM matrix file to use for VNE relative conservation []
    
    Opt:boolean
    - AlnGap = Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore as putative sequence fragments [True]
    - ConsAmb = Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    - ConsInfo = Weight positions by information content [True]
    - FullForce = Whether to force regeneration of alignments using GOPHER
    - HomFilter = Whether to filter homologues using seqfilter options [False]
    - Masking = Whether to use seq.info['MaskSeq'] for Prob cons, if present (else 'Sequence') [True]
    - RelGapPen = Whether to invoke the "Gap Penalty" during relative conservation calculations [True]
    - UseAln = Whether to look for conservation in alignments
    - UseGopher = Use GOPHER to generate missing orthologue alignments in outdir/Gopher - see gopher.py options [False]
    
    Stat:numeric
    - ConsWeight = Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
    - MinHom = Minimum number of homologues for making conservation score [1]
    - Percentile = Percentile steps to return in addition to mean [25]
    - RelConWin = Determines window size for RelCons calculation [30]
    - WinSize = Used to define flanking regions for stats. If negative, will use flanks *only* [0]

    List:list
    - Alphabet = List of letters in alphabet of interest
    - GopherRun = List of sequence objects that have already been put through GOPHER. (Do not keep retrying!)
    - Headers - List of Headers for combined occurrence stats
    - OccHeaders - List of Headers for individual occurrences
    - Percentile - List of Percentiles to return for combined stats (set up with Headers)
    - SLiMCalc - List of occurrence statistics to calculate []
    - SlimFilter - List of stats to filter (remove matching) SLiMs on, consisting of X*Y  []
                      - X is an output stat (the column header),
                      - * is an operator in the list >, >=, !=, =, >= ,<    
                      - Y is a value that X must have, assessed using *.
                      This filtering is crude and may behave strangely if X is not a numerical stat!
                      !!! Remember to enclose in "quotes" for <> filtering !!!
    - OccFilter = Same as slimfilter but for individual occurrences []

    Dict:dictionary
    - ConsSpecLists = Dictionary of {BaseName:List} lists of species codes for special conservation analyses
    - ElementIC = Element IC values
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) {}
    - BLOSUM = BLOSUM matrix loaded from VNE calculation {}

    Obj:RJE_Objects
    - AAPropMatrix = rje_aaprop.AAPropMatrix object
    - Gopher = Gopher Fork object for alignment generation
    - IUPred = Disorder object for running IUPred disorder
    - FoldIndex = Disorder object for running FoldIndex disorder
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AlnDir','AlnExt','ConScore','GopherDir','PosMatrix','VNEMatrix']
        self.optlist = ['AlnGap','ConsAmb','ConsInfo','FullForce','UseAln','UseGopher','Masking','HomFilter','RelGapPen']
        self.statlist = ['ConsWeight','Percentile','WinSize','RelConWin','MinHom']
        self.listlist = ['Alphabet','Headers','OccHeaders','Percentile','SLiMCalc','SlimFilter','OccFilter','GopherRun']
        self.dictlist = ['ConsSpecLists','PosMatrix','ElementIC','BLOSUM']
        self.objlist = ['AAPropMatrix','Gopher','IUPred','FoldIndex']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0,obj=None,setlist=True,setdict=True)
        self.setInfo({'AlnDir':rje.makePath('',return_blank=True),'GopherDir':rje.makePath('./',return_blank=False),
                      'AlnExt':'aln.fas','ConScore':'rel'})
        self.setOpt({'AlnGap':True,'ConsAmb':True,'ConsInfo':True,'RelGapPen':True})
        self.setStat({'RelConWin':30})
        self.obj['IUPred'] = rje_disorder.Disorder(self.log,self.cmd_list+['disorder=iupred'])
        self.obj['FoldIndex'] = rje_disorder.Disorder(self.log,self.cmd_list+['disorder=foldindex'])
        self.list['Alphabet'] = rje_slim.default_aas
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] Read in commands from cmd_list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._cmdReadList(cmd.lower(),'info',['ConScore'])  
                self._cmdReadList(cmd,'info',['AlnExt'])  
                self._cmdReadList(cmd,'file',['ConScore','PosMatrix','VNEMatrix'])  
                self._cmdReadList(cmd,'path',['AlnDir','GopherDir'])
                self._cmdReadList(cmd,'opt',['AlnGap','ConsAmb','ConsInfo','FullForce','Masking','UseAln','UseGopher',
                                             'HomFilter','RelGapPen'])  
                self._cmdReadList(cmd,'int',['ConsWeight','WinSize','RelConWin','MinHom'])
                self._cmdReadList(cmd,'stat',['Percentile'])
                self._cmdReadList(cmd,'list',['SlimFilter','OccFilter'])
                self._cmdReadList(cmd.lower(),'list',['SLiMCalc'])
                self._cmdReadList(cmd.upper(),'list',['Alphabet'])
                ### Special ###
                if rje.matchExp('^conspec=(.+)',cmd.lower()):
                    conspec = rje.matchExp('^conspec=(.+)',cmd.lower())[0]
                    consfiles = glob.glob(conspec)
                    if len(consfiles) > 0:     # File
                        for cfile in consfiles:
                            self.dict['ConsSpecLists'][rje.baseFile(cfile,True)] = rje.listFromCommand(cfile)
                    else: self.dict['ConsSpecLists']['SPEC'] = rje.listFromCommand(conspec,checkfile=False)
            except: self.log.errorLog('Problem with cmd:%s' % cmd)

        ### ~ [2] Adjust Conservation Attributes etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if 'cons' in self.list['SLiMCalc']:
            self.opt['UseAln'] = True
            if self.info['ConScore'].lower() in ['all','prop']: self.obj['AAPropMatrix'] = rje_aaprop.AAPropMatrix(self.log,self.cmd_list)
            if self.info['ConScore'].lower() in ['all','prop','pos']: self.posMatrix()
            if self.info['VNEMatrix'].lower() not in ['','none']: self.loadBLOSUM()

        ### ~ [3] Setup GopherDir and AlnDir for alignments and conservation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.setupGopher()
#########################################################################################################################
    def setupGopher(self):  ### Sets up GOPHER directory etc.                                                       #V0.5
        '''Sets up GOPHER directory etc.'''
        try:### ~ [1] Create/Modify Directory info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['Gopher'] = None
            if not os.path.exists(self.info['AlnDir']): self.info['AlnDir'] = ''
            if self.getStr('GopherDir').lower() not in ['','none']:
                gcmd = ['orthalign'] + self.cmd_list + ['gnspacc=T','i=-1','autofilter=F']
                self.obj['Gopher'] = gopher.Gopher(self.log,gcmd)
                if self.opt['UseGopher']:
                    self.obj['Gopher'].obj['BLAST'] = self.obj['Gopher'].setupBlast()  
                    self.obj['Gopher'].obj['BLAST'].log = self.log
                    if self.obj['Gopher'].opt['FullForce']: self.log.printLog('#ALN','FullForce=T. Will call Gopher regardless of existing files')
            try:
                if not self.info['AlnExt'][0] == '.': self.info['AlnExt'] = '.%s' % self.info['AlnExt']
            except: pass    # Presumably no alnext!
        except:
            self.errorLog('Problem with SLiMCalc.setupGopher(). UseGopher cancelled')
            self.opt['UseGopher'] = False
#########################################################################################################################
    def posMatrix(self):    ### Loads and builds PosMatrix for Conservation Scoring
        '''Loads and builds PosMatrix for Conservation Scoring.'''
        try:
            ### ~ [1] Setup PosMatrix dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['PosMatrix'] = {}
            ## ~ [1a] Generate from property matrix if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['ConScore'].lower() in ['all','prop']:
                propx = len(self.obj['AAPropMatrix'].prop)
                for a1 in self.list['Alphabet']:
                    for a2 in self.list['Alphabet']:
                        self.dict['PosMatrix']['%s%s' % (a1,a2)] = float(propx - self.obj['AAPropMatrix'].pdif['%s%s' % (a1,a2)]) / float(propx)
                self.dict['PropPosMatrix'] = self.dict.pop('PosMatrix')
                self.dict['PosMatrix'] = {}
                if self.info['ConScore'].lower() in ['prop']: return
            ## ~ [1b] Make default PosMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for aa in self.list['Alphabet']: self.dict['PosMatrix']['%s%s' % (aa,aa)] = 1.0

            ### ~ [2] Load PosMatrix file if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Look for PosMatrix File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            posfile = self.info['PosMatrix']
            if posfile.lower() in ['','none']: return
            if not os.path.exists(posfile): return self.log.errorLog('PosMatrix file "%s" not found!' % posfile,printerror=False)
            ## ~ [2b] Load Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            plines = self.loadFromFile(posfile)
            _alphabet = plines[0].split()
            if len(_alphabet) > 1:  # Treat as matrix
                for a in _alphabet:
                    line = plines[_alphabet.index(a)+1].split()
                    if len(line) != (len(self.alphabet)+1):
                        self.log.errorLog('%s has wrong format! Does not match %s' % (line, _alphabet))
                        self.info['PosMatrix'] = 'None'
                        self.posMatrix()
                        raise ValueError
                    for i in range(1,len(line)):
                        score = float(line[i])
                        if score: self.dict['PosMatrix']['%s%s' % (a,_alphabet[i-1])] = score
                self.log.printLog('#CONS','Position Scoring Matrix set from %s as matrix.' % posfile)
            else:   # Treat as equivalence files
                for line in plines:
                    match = rje.matchExp('^(\S+)',line)
                    if match:
                        for a1 in match[0]:
                            for a2 in match[0]: self.dict['PosMatrix']['%s%s' % (a1,a2)] = 1.0
                self.log.printLog('#CONS','Position Scoring Matrix set from %s as equivalences.' % posfile)

        except: self.log.errorLog('Major problem in Presto.posMatrix(%s)' % self.info['PosMatrix'],quitchoice=True)
#########################################################################################################################
    ### <2> ### Setup Headers etc.                                                                                      #
#########################################################################################################################
    def setupHeaders(self): ### Sets up Headers and OccHeaders lists based on attribute settings (and filters)
        '''Sets up Headers and OccHeaders lists based on attribute settings.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Headers'] = []
            self.list['OccHeaders'] = []
            for i in range(len(self.list['SLiMCalc'])):
                att = self.list['SLiMCalc'][i]
                try: self.list['SLiMCalc'][i] = occstats[att].lower()
                except:
                    self.list['SLiMCalc'][i] = self.list['SLiMCalc'][i].lower()
                    self.printLog('#CALC','Treating "%s" as a disorder score method' % self.list['SLiMCalc'][i])
                    #self.list['SLiMCalc'][i] = 'err'
                    #self.log.errorLog('SLiMCalc stat "%s" not recognised!' % att)
            self.list['SLiMCalc'] = rje.sortUnique(self.list['SLiMCalc'],False,False)
            if 'err' in self.list['SLiMCalc']: self.list['SLiMCalc'].remove('err')
            ## ~ [1a] Percentiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Percentile'] = []
            if self.stat['Percentile'] > 0:
                pc = 100.0
                while pc > 0:
                    self.list['Percentile'].append(pc)
                    pc -= self.stat['Percentile']
                self.list['Percentile'].append(0)
            
            ### ~ [2] Occurrence Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for att in self.list['SLiMCalc']:
                if att in ['chg']: self.list['OccHeaders'] += ['AbsChg','NetChg']
                elif att in ['cons']:
                    for taxa in ['ALL'] + rje.sortKeys(self.dict['ConsSpecLists']):
                        if taxa == 'ALL': self.list['OccHeaders'] += ['Cons','HomNum','GlobID','LocID']
                        else: self.list['OccHeaders'] += ['%s_Cons' % taxa,'%s_HomNum' % taxa,'%s_GlobID' % taxa,'%s_LocID' % taxa]
                        if self.info['ConScore'] == 'all':
                            for method in ['Abs','Pos','Prob','Prop','Rel']: self.list['OccHeaders'].append('%s_Cons_%s' % (taxa,method))
                elif att in occstats: self.list['OccHeaders'] += [occstats[att]]
                else: self.list['OccHeaders'].append(att)

            ### ~ [3] Combined attribute Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for head in self.list['OccHeaders']:
                self.list['Headers'].append('%s_mean' % head)
                for pc in  self.list['Percentile']: self.list['Headers'].append('%s_pc%d' % (head,int(pc+0.5)))

        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <3> ### Occurrence Attribute Calculation Methods I: General                                                     #
#########################################################################################################################
    def occStats(self,occlist,xpad=0,progress=False,silent=False):   ### Calculates general occurrence stats for occlist
        '''
        Calculates general occurrence stats for occlist - should all have same Seq object. 
        >> occlist:list of MotifOcc objects to calculate stats for (must all have same Seq)
        >> xpad:int [0] = Xs to be added to either side of sequence
        >> progress:bool [False] = whether to print progress to screen or not
        >> silent:bool [False] = whether to make verbosity -1 for duration of occStats 
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not occlist: return
            v = self.log.stat['Verbose']
            if silent: self.log.stat['Verbose'] = -1
            calc = self.list['SLiMCalc']
            ## ~ [1a] Setup sequence object and attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            Seq = occlist[0]['Seq']
            sequence = 'X' * xpad + Seq.getSequence(gaps=False) + 'X' * xpad
            calcdict = {'chg':sequence,'comp':sequence}
            if 'sa' in calc: calcdict['sa'] = rje_seq.surfaceAccessibility(sequence,returnlist=True)
            if 'hyd' in calc: calcdict['hyd'] = rje_seq.eisenbergHydropathy(sequence,returnlist=True)
            ## ~ [1b] Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'iup' in calc:
                if Seq.list.has_key('IUPred'): calcdict['iup'] = Seq.list['IUPred'][0:]
                else:
                    self.obj['IUPred'].disorder(sequence=sequence,name=Seq.shortName())
                    calcdict['iup'] = self.obj['IUPred'].list['ResidueDisorder'][0:]
                    Seq.list['IUPred'] = self.obj['IUPred'].list['ResidueDisorder'][0:]
            if 'fold' in calc:
                if Seq.list.has_key('FoldIndex'): calcdict['fold'] = Seq.list['FoldIndex'][0:]
                else:
                    self.obj['FoldIndex'].disorder(sequence=sequence,name=Seq.shortName())
                    calcdict['fold'] = self.obj['FoldIndex'].list['ResidueDisorder'][0:]
                    Seq.list['FoldIndex'] = self.obj['FoldIndex'].list['ResidueDisorder'][0:]
            #i# Additional disorder methods
            for att in calc:
                if att == 'cons': continue  # Conservation is special and handled separately below.
                if att not in calcdict:
                    if att == 'dis': self.obj[att] = rje_disorder.Disorder(self.log,self.cmd_list)
                    else: self.obj[att] = rje_disorder.Disorder(self.log,self.cmd_list+['disorder=%s' % att])
                    self.debug(att)
                    self.obj[att].disorder(sequence=sequence,name=Seq.shortName())
                    calcdict[att] = self.obj[att].list['ResidueDisorder'][0:]
                    Seq.list[att] = self.obj[att].list['ResidueDisorder'][0:]

            #!# seq_dom = seqDom(self,Seq,seq_dis) #!# Not yet implementd. See rje_motif_stats. #!#
            
            ### ~ [2] Calculation of Occurrence Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ox = 0.0
            for Occ in occlist:
                ## ~ [2a] Setup occurrence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if progress: self.log.printLog('\r#CALC','Calculating %s occurrence stats %.1f%%' % (Seq.shortName(),ox/len(occlist)),newline=False,log=False)
                ox += 100.0
                r = int(Occ['Pos'] - 1) + xpad
                if Occ['Pos'] < 1: r += 1  # XPadDB
                match = Occ['Match']
                ## ~ [2b] Calculate SA, Hyd, Charge, Complexity and disorder ~~~~~~~~~~~~~~~~~~~~~~ ##
                for att in calc:
                    if att == 'cons': continue  # Conservation is special and handled separately below.
                    if att in occstats:
                        stat = occstats[att]
                    else:
                        stat = att
                    ## Set window ##
                    win = self.stat['WinSize']    # Will return 0 if missing
                    if win >= 0: w = r - win
                    else: w = r + win
                    if w < 0: w = 0
                    ## Get appropriate Data region ##
                    if win >= 0: winreg = calcdict[att][w:r+len(match)+win]     # Region includes occurrence
                    else: winreg = calcdict[att][w:r] + calcdict[att][r+len(match):r+len(match)-win]    # Flanks only
                    if att == 'chg':
                        chgdict = rje_sequence.chargeDict(winreg)       #!# Check this #!#
                        for chg in chgdict: Occ[chg] = chgdict[chg]
                    elif att == 'comp':     # Complexity of region in terms of maximum complexity
                        Occ[stat] = float(len(rje.sortUnique(rje.strList(winreg)))) / min(20,len(winreg))
                    else:                        
                        Occ[stat] = 0   # Note that the stat is now Hyd not Hydropathy!
                        if winreg: Occ[stat] = sum(winreg) / len(winreg)
                ## ~ [2c] Additional calculations for possible use ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Domain Filtering has not been implemented. See rje_motif_stats to implement. #!#
                ## Peptide Design ##
                pepwin = rje.modulus(self.stat['WinSize'])
                w = r - pepwin
                if w < 0: w = 0
                Occ['PepSeq'] = sequence[w:r+len(match)+self.stat['WinSize']]     # Region for calculation
            if progress: self.log.printLog('\r#CALC','Calculating %s occurrence stats complete.' % Seq.shortName())

            ### ~ [3] Special calculations of SLiM conservation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'cons' in calc: self.hitAlnCon(occlist,progress)

        except: self.log.errorLog(rje_zen.Zen().wisdom())
        self.log.stat['Verbose'] = v
#########################################################################################################################
    ### <4> ### Occurrence Attribute Calculation Methods II: Conservation                                               #
#########################################################################################################################
    def hitAlnCon(self,occlist,progress=True):    ### Looks for alignment and, if appropriate, calculate conservation stats.
        '''
        Looks for alignment and, if appropriate, calculate conservation stats.
        
        Any homologues with masked (X) residues that coincide to non-wildcard positions of the motif occurrence will be
        ignored from conservation calculations. Gaps, however, shall be treated as divergence. The exception is that when
        the alngap=F option is used, 100% gapped regions of homologues are also ignored.

        This method deals with all the occurrences of all motifs for a single sequence and its alignment. Global
        alignment statistics are calculated first, then each occurrence for each motif is processed. Since version 1.1,
        subtaxa are treated the same as all taxa to reduce the coding: the default all taxa is now effectively an
        additional subtaxa set.
        
        >> occlist:list of MotifOcc objects to calculate stats for (must all have same Seq)
        >> progress:bool [True] = whether to print progress to screen or not
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not occlist: return
            OccSeq = occlist[0]['Seq']
            ## ~ [1a] Default Values (0.0) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for Occ in occlist:
                for h in self.list['OccHeaders']:
                    if h.find('Cons') >= 0: 
                        if self.info['ConScore'] in ['rel','rlc']: Occ[h] = 1.0   #!# Default 0.0 if no homologues #!#
                        else: Occ[h] = 1.0   #!# Default 1.0 if no homologues #!#
                    if h.find('HomNum') >= 0: Occ[h] = 0
                    if h.find('GlobID') >= 0 or h.find('LocID') >= 0: Occ[h] = 0.0
            ## ~ [1b] Identify/Create and Load Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aln = self.loadOrthAln(OccSeq)  # Try to load alignment, using GOPHER if appropriate.
            if not aln:                     # Alignment rejected during loadOrthAln. No conservation statistics here.
                self.printLog('#ALN','No alignment file found for %s.' % OccSeq.shortName(),screen=progress)
                return      
            ## ~ [1c] Setup alignment attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = aln.obj['QuerySeq']
            seqs = aln.seq[0:]      # This is a list of Sequence objects, minus the query, which is stored in qry
            seqs.remove(qry)
            self.printLog('#ALN','Comparing %d homologues for %s.' % (len(seqs),qry.shortName()),screen=progress)
            if len(seqs) < self.stat['MinHom']:
                self.printLog('#HOM','Too few (< %d) homologues for %s.' % (self.stat['MinHom'],qry.shortName()),screen=progress)
                return

            ### ~ [2] Calculate Global Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Sequence-specific statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqid = {}      # Dictionary of global %ID vs query for each sequence (consweight<>0)
            seqwt = {}      # Dictionary of global %ID weighting for each sequence (consweight<>0)
            for seq in seqs:
                igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'],seq.info['Sequence'],nomatch=['X'])
                seqid[seq] = float(igedic['ID'][0]) / igedic['Len'][0]
                seqwt[seq] = seqid[seq] ** self.stat['ConsWeight']
            ## ~ [2b] Taxa-specific Sequence lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            conseqs = {'ALL':seqs[0:]}    # Dictionary of sequence lists {speclist:seqs}
            for taxa in rje.sortKeys(self.dict['ConsSpecLists']):
                conseqs[taxa] = []
                codelist = self.dict['ConsSpecLists'][taxa]
                for seq in seqs:
                    if seq.info['SpecCode'] in codelist or seq.info['Species'] in codelist: conseqs[taxa].append(seq)
            ## ~ [2c] Taxa-specific Global ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            globid = {}     # Dictionary of mean global %ID for each taxonomic grouping
            for taxa in conseqs.keys():
                globid[taxa] = 0.0
                for seq in conseqs[taxa]: globid[taxa] += seqid[seq]
                ## Average Global identity ##
                if len(conseqs[taxa]) > 0: globid[taxa] /= len(conseqs[taxa])
            ## ~ [2d] Probability-based measure setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['ConScore'] in ['all','prob']:
                subdict = aln.setupSubDict(masking=self.opt['Masking'],alphabet=self.list['Alphabet'])
                #x#print seqConsList(OccSeq,self)
            ## ~ [2e] Setup Relative conservation measure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['ConScore'] in ['all','rel','rlc']:
                taxarel = {'ALL':self.relConList(qry,seqs,seqwt,window=self.stat['RelConWin'])}   # List of relative conservation scores for different taxa groups
                for taxa in rje.sortKeys(self.dict['ConsSpecLists']): taxarel[taxa] = self.relConList(qry,conseqs[taxa],seqwt,window=self.stat['RelConWin']) 

            ### ~ [3] Process Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
            for Occ in occlist:
                Motif = Occ['Motif']

                ## ~ [3a] Find position in aligment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (start,end) = self.findOccPos(Occ,qry)      
                if (start,end) == (-1,-1): continue # Could not find!
                    
                ## ~ [3b] Homologues, Sequence Fragments and Local Identity ~~~~~~~~~~~~~~~~~~~~~~~ ##
                hithom = {}     # Dictionary of {taxa:homologues} for this occurrence
                locid = {}      # Dictionary of {seq:local %ID} and {taxa:local %ID}
                alnfrag = {qry:qry.info['Sequence'][start:end]} # Partial alignments across occurrence
                seqfrag = {}    # Dictionary of {seq:degapped sequence frag across occurrence}
                for taxa in conseqs.keys():
                    hithom[taxa] = []
                    locid[taxa] = 0.0
                for seq in seqs:
                    igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'][start:end],seq.info['Sequence'][start:end],nomatch=['X'])
                    locid[seq] = float(igedic['ID'][0]) / igedic['Len'][0]
                    alnfrag[seq] = seq.info['Sequence'][start:end]
                    seqfrag[seq] = rje_seq.deGap(seq.info['Sequence'][start:end])
                    if (len(seqfrag[seq]) >= 1 or self.opt['AlnGap']) and len(seqfrag[seq]) > string.count(seqfrag[seq].upper(),'X'):
                        for taxa in conseqs.keys():
                            if seq in conseqs[taxa]:
                                hithom[taxa].append(seq)    # Seq only in hithom if right taxa and not all gap/X
                                locid[taxa] += locid[seq]
                for taxa in conseqs.keys():
                    if hithom[taxa]:  locid[taxa] /= len(hithom[taxa])

                ## ~ [3c] Reduced alignment for Pos and Prop methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.info['ConScore'] in ['all','pos','prop']:
                    red_aln = {}    # Dictionary of {seq:sequence to consider}
                    for seq in [qry]+seqs:
                        red_aln[seq] = ''
                        m = 0   # Number of positions checked
                        for r in range(len(alnfrag[qry])):
                            a = alnfrag[qry][r]
                            if a == '-': continue   # Skip
                            if alnfrag[seq][r] == '-': red_aln[seq] += '-'
                            elif Occ['Variant'][m] in ['X','*']: red_aln[seq] += 'X'   # Not important for match
                            else: red_aln[seq] += alnfrag[seq][r]
                            m += 1
                        if m != len(Occ['Match']):
                            self.log.errorLog('Something wrong with Reduced alignment: %s' % red_aln.values(),printerror=False)
                            raise ValueError
                    #self.deBug(red_aln)
                    # >> red_aln sequence should now have an X at each wildcard position, otherwise the relevant residue from the alignment << #

                ## ~ [3d] Conservation Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for taxa in conseqs.keys():
                    hitcons = {'ABS':0.0,'POS':0.0,'PROP':0.0,'PROB':0.0,'REL':0.0}
                    if self.info['ConScore'] in ['all','abs']:
                        hitcons['ABS'] = self.absCons(Occ,hithom[taxa],seqfrag,seqwt)
                    if self.info['ConScore'] in ['all','pos']:                        
                        hitcons['POS'] = self.posCons(Occ,hithom[taxa],red_aln,seqwt,self.dict['PosMatrix'])
                    if self.info['ConScore'] in ['all','prop']:
                        hitcons['PROP'] = self.posCons(Occ,hithom[taxa],red_aln,seqwt,self.dict['PropPosMatrix'])
                    if self.info['ConScore'] in ['all','prob']:
                        hitcons['PROB'] = self.probCons(Occ,hithom[taxa],red_aln,seqwt,subdict)
                    if self.info['ConScore'] in ['all','rel','rlc']:
                        hitcons['REL'] = self.relCons(Occ,taxarel[taxa])
                    #x#print taxa, hitcons

                ## ~ [3e] Update Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if taxa == 'ALL': 
                        Occ['HomNum'] = len(hithom[taxa])
                        Occ['GlobID'] = rje_slim.expectString(globid[taxa])
                        Occ['LocID'] = rje_slim.expectString(locid[taxa])
                    else:
                        Occ['%s_HomNum' % taxa] = len(hithom[taxa])
                        Occ['%s_GlobID' % taxa] = rje_slim.expectString(globid[taxa])
                        Occ['%s_LocID' % taxa] = rje_slim.expectString(locid[taxa])
                    if self.info['ConScore'] == 'all':
                        for method in hitcons.keys():
                            outmethod = method.upper()[:1] + method.lower()[1:]
                            Occ['%s_Cons_%s' % (taxa,outmethod)] = hitcons[method]
                        prob = hitcons.pop('PROB')  #!# Don't know how to add this to other methods #!#
                        rel = hitcons.pop('REL')    #!# Don't know how to add this to other methods #!#
                        if taxa == 'ALL': Occ['Cons'] = sum(hitcons.values()) / len(hitcons)
                        else: Occ['%s_Cons' % taxa] = sum(hitcons.values()) / len(hitcons)
                    else:
                        if taxa == 'ALL': Occ['Cons'] = hitcons[self.info['ConScore'].upper()]
                        else: Occ['%s_Cons' % taxa] = hitcons[self.info['ConScore'].upper()]

        except: self.errorLog('Error in rje_slimcalc.hitAlnCon()',quitchoice=True); raise
#########################################################################################################################
    def absCons(self,Occ,hithom,seqfrag,seqwt):    ### Absolute conservation score.
        '''Absolute conservation score.'''
        try:### ~ [1] Absolute matching of motif in corresponding homologous region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = Occ['Motif']
            hitcon = {}     # Dictionary of {seq:conservation}
            for seq in hithom:
                hitcon[seq] = 0.0
                ## ~ [1a] Select variants allowed to match ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.opt['ConsAmb']: vlist = Motif.dict['Search'][0]     # Search degenerate motif
                else: vlist = [Occ['Variant']]                      # Search with matched variant
                ## ~ [1b] Look for matching variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for variant in vlist:
                    searchvar = '(%s)' % string.replace(variant,'X','[A-Z]')
                    if rje.matchExp(searchvar,seqfrag[seq]):
                        hitcon[seq] = 1.0
                        break
            ### ~ [2] Weight by distance? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.consWeight(hitcon,seqwt)
        except: return self.log.errorLog('Error in rje_slimcalc.absCons()',quitchoice=True) 
#########################################################################################################################
    def consWeight(self,hitcon,seqwt):    ### Weights conservation and returns final conservation score
        '''
        Weights conservation and returns final conservation score.
        >> hitcon: raw dictionary of {seq:conservation}
        >> seqwt: weighting dictionary of {seq:weighting}
        << cons: conservation score
        '''
        try:### ~ [1] Weight conservation of sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cons = 0.0
            totweight = 0.0
            for seq in hitcon.keys():
                totweight += seqwt[seq]
                cons += hitcon[seq] * seqwt[seq]
                if hitcon[seq] > 1.0: self.log.errorLog('HitCon = %.3f > 1 (%s)!' % (hitcon[seq],seq.shortName()),printerror=False)
            if totweight: return cons / totweight
            return 0.0
        except:
            self.log.errorLog('Error in rje_slimcalc.absCons()',quitchoice=True) 
            return -1.0
#########################################################################################################################
    def relCons(self,Occ,relcon):    ### Positional conservation score.
        '''
        Relative conservation score.
        >> Occ:Motif occurrence dictionary.
        >> relcon:list = Relative conservation score for each residue
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infowt = []     # IC at each position for weighting
            ## ~ [1a] Setup Information Content Weighting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['ConsInfo']: searchvar = Occ['SearchVar'][0:]
            else: searchvar = Occ['Variant'][0:]
            tvar = searchvar[0:]
            while searchvar:
                if searchvar[0] in ['$','^']: searchvar = searchvar[1:]; continue   #!# Deal with better later? #!#
                elif searchvar.find('[A-Z]') == 0:
                    el = 'X'
                    searchvar = searchvar[5:]
                elif searchvar[0] == '[':
                    el = searchvar[1:searchvar.find(']')]
                    searchvar = searchvar[searchvar.find(']')+1:]
                else:
                    el = searchvar[:1]
                    searchvar = searchvar[1:]
                if el not in ['(',')']:
                    if not self.dict['ElementIC'].has_key(el): self.dict['ElementIC'][el] = rje_slim.elementIC(el)
                    infowt.append(self.dict['ElementIC'][el])
            infosum = sum(infowt)
            if len(infowt) != len(Occ['Match']):
                self.log.errorLog('Search variant (%s) does not match length of match (%s)' % (tvar,Occ['Match']),printerror=False)
                return -1.0
            ### ~ [2] Relative conservation score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cons = 0.0
            #print relcon
            #print infowt
            #print Occ
            for i in range(len(infowt)): cons += relcon[Occ['Pos']-1+i] * infowt[i]   ### Each residue in turn
            return cons / infosum       # Return weighted mean
        except:
            self.errorLog('Error in rje_slimcalc.posCons() %s' % Occ,quitchoice=False) 
            return -1.0
#########################################################################################################################
    def posCons(self,Occ,hithom,red_aln,seqwt,posmatrix):    ### Positional conservation score.
        '''Positional conservation score.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infowt = []     # IC at each position for weighting
            ## ~ [1a] Setup Information Content Weighting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['ConsInfo']: searchvar = Occ['SearchVar'][0:]
            else: searchvar = Occ['Variant'][0:]
            tvar = searchvar[0:]
            while searchvar:
                if searchvar[0] in ['$','^']: searchvar = searchvar[1:]; continue   #!# Deal with better later? #!#
                elif searchvar.find('[A-Z]') == 0:
                    el = 'X'
                    searchvar = searchvar[5:]
                elif searchvar[0] == '[':
                    el = searchvar[1:searchvar.find(']')]
                    searchvar = searchvar[searchvar.find(']')+1:]
                else:
                    el = searchvar[:1]
                    searchvar = searchvar[1:]
                if el not in ['(',')']:
                    if not self.dict['ElementIC'].has_key(el): self.dict['ElementIC'][el] = rje_slim.elementIC(el)
                    infowt.append(self.dict['ElementIC'][el])
            infosum = sum(infowt)
            ## ~ [1b] Setup VarList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            varlist = []    # List of AAs at each position for checking (just use hit.info['Variant'] if not self.opt['ConsAmb']
            if self.opt['ConsAmb']:       
                searchvar = Occ['SearchVar'][0:]
                while searchvar:
                    if searchvar.find('[A-Z]') == 0:
                        el = 'X'
                        searchvar = searchvar[5:]
                    elif searchvar[0] == '[':
                        el = searchvar[1:searchvar.find(']')]
                        searchvar = searchvar[searchvar.find(']')+1:]
                    else:
                        el = searchvar[:1]
                        searchvar = searchvar[1:]
                    if el not in ['(',')']: varlist.append(el)
            else:       # Search matched variant only
                varlist = Occ['Motif'].variant(Occ)   # Can use "for X in Y" for both lists and strings, so this is OK!
            #X#self.deBug('%s varlist = %s => %s' % (hit.info['Name'],varlist,infowt))
            ### ~ [1c] Check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if len(infowt) != len(Occ['Match']):
                self.log.errorLog('Search variant (%s) does not match length of match (%s)' % (tvar,Occ['Match']),printerror=False)
                return -1.0

            ### ~ [2] Positional conservation score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitcon = {}     # Dictionary of {seq:conservation}
            for seq in hithom:
                hitcon[seq] = 0.0
                posvarlist = varlist[0:]
                if '^' in posvarlist: posvarlist.remove('^')
                if '$' in posvarlist: posvarlist.remove('$')
                try:
                    for p in range(len(posvarlist)):
                        a = red_aln[seq][p]
                        hitcon[seq] += self.bestScore(a,varlist[p],posmatrix) * infowt[p]
                    hitcon[seq] /= infosum
                except:
                    self.deBug(varlist)
                    self.deBug(posvarlist)
                    self.deBug('Occ: %s' % Occ)
                    self.deBug('hithom: %s' % hithom)
                    self.deBug('red_aln: %s' % red_aln)
                    self.deBug('seqwt: %s' % seqwt)
                    self.deBug('posmatrix: %s' % posmatrix)
                    self.log.errorLog('Problem with %s posCons(%s) aln "%s" pos %d' % (Occ['Match'],seq.shortName(),red_aln[seq],p),printerror=True,quitchoice=True)
                    hitcon[seq] = 0.0

            ### ~ [3] Weight by distance? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.consWeight(hitcon,seqwt)
        except:
            self.errorLog('Error in rje_slimcalc.posCons()',quitchoice=True) 
            return -1.0
#########################################################################################################################
    def bestScore(self,aa,aalist,posmatrix):  ### Best score for compared AAs.
        '''
        >> aa:str = Amino acid to be compared
        >> aalist:str = List of aas to compare to
        >> posmatrix:dict of {'a1a2':score} to get score from
        '''
        if aa == 'X' or aalist == 'A-Z': return 0.0
        if aa in aalist: return 1.0
        best = 0.0
        for a in aalist:
            if posmatrix.has_key('%s%s' % (aa,a)) and posmatrix['%s%s' % (aa,a)] > best: best = posmatrix['%s%s' % (aa,a)]
            elif posmatrix.has_key('%s%s' % (a,aa)) and posmatrix['%s%s' % (a,aa)] > best: best = posmatrix['%s%s' % (a,aa)]
        return best
#########################################################################################################################
    def probCons(self,Occ,hithom,red_aln,seqwt,subdict):    ### Probabilistic conservation score.
        '''Positional conservation score.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup VarList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qryvar = Occ['Motif'].variant(Occ)      # Actual variant matching query
            varlist = []    # List of AAs at each position for checking (just use hit.info['Variant'] if not self.opt['ConsAmb']
            if self.opt['ConsAmb']: varlist = rje_slim.prestoFromCode(rje_slim.slimFromPattern(Occ['Variant']))
            else:       # Search matched variant only
                varlist = qryvar   # Can use "for X in Y" for both lists and strings, so this is OK!
            #X#print '%s (%s) => %s' % (Occ['Variant'],self.opt['ConsAmb'],varlist)

            ### ~ [2] Probabilistic conservation score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitcon = {}     # Dictionary of {seq:conservation}
            for seq in hithom:
                ## ~ [2a] Generate probability of observed conservation ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                probcon = 1.0
                posvarlist = varlist[0:]
                if '^' in posvarlist: posvarlist.remove('^')
                if '$' in posvarlist: posvarlist.remove('$')
                try:
                    for i in range(len(posvarlist)):
                        aa = qryvar[i]  # Actual query AA
                        if aa in ['X','.']: continue  # Not interested in wildcards
                        sa = red_aln[seq][i]    # Aligned sequence
                        va = posvarlist[i]         # Variant options
                        if sa in va:    # Add to probability if actually conserved
                            pv = 0.0    # Probability of substitution with observed variant
                            for qa in va: pv += rje.getFromDict(subdict[seq][aa],qa,False,False,default=0.0)
                            probcon *= min(1.0,pv)
                except:
                    self.log.errorLog('Problem with %s posCons(%s) aln "%s" pos %d' % (Occ['Match'],seq.shortName(),red_aln[seq],i),printerror=True,quitchoice=True)
                    probcon = 1.0
                ## ~ [2b] Generate probability of observed divergence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                probdiv = 1.0
                try:
                    for i in range(len(posvarlist)):
                        aa = qryvar[i]  # Actual query AA
                        if aa in ['X','.']: continue  # Not interested in wildcards
                        sa = red_aln[seq][i]    # Aligned sequence
                        va = posvarlist[i]         # Variant options
                        if sa not in va:    # Add to probability if not conserved
                            pv = 1.0        # Probability of substitution with no variant
                            for qa in va: pv -= rje.getFromDict(subdict[seq][aa],qa,False,False,default=0.0)
                            probdiv *= max(0.0,pv)
                except:
                    self.log.errorLog('Problem with %s posCons(%s) aln "%s" pos %d' % (Occ['Match'],seq.shortName(),red_aln[seq],i),printerror=True,quitchoice=True)
                    probdiv = 1.0
                ## ~ [2c] Calculate ratio for final score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if probcon: hitcon[seq] = probdiv / probcon
                else:
                    self.log.errorLog('Problem with %s posCons(%s): ProbCon = 0.0' % (Occ['Match'],seq.shortName()),printerror=True,quitchoice=True)
                    hitcon[seq] = probdiv * probdiv

            ### ~ [3] Weight by distance? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return rje.geoMean(hitcon.values())
        except:
            self.log.errorLog('Error in rje_slimcalc.probCons()',quitchoice=True) 
            return -1.0
#########################################################################################################################
    def loadOrthAln(self,seq,usegopher=True,log=True,screen=True,alnfile=False):    ### Identifies file, loads and checks alignment.                      #V0.5
        '''
        Identifies file, loads and checks alignment. If the identified file is not actually aligned, then RJE_SEQ will try to
        align the proteins using MUSCLE or ClustalW.
        >> seq:Sequence being analysed.
        >> usegopher:bool [True] = whether to try and use GOPHER if appropriate (set False if already tried)
        >> alnfile:bool [False] = Whether to return the alignment file rather than the SeqList object.
        << aln = SeqList object containing alignment with queryseq
        '''
        try:### ~ [1] Setup Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gobj = self.obj['Gopher']
            filelist = []
            if self.opt['UseAln']:
                for alnstart in [seq.info['AccNum'],seq.info['ID'],seq.shortName()]:
                    filelist.append('%s%s%s' % (self.info['AlnDir'],alnstart,self.info['AlnExt']))
            if gobj:
                gopherdir = gobj.gopherDir(seq,mkdir=False)
                gfork = gopher.GopherFork(self.log,self.cmd_list)
                gfork.setInfo({'Name':seq.shortName(),'OutPath':gopherdir})
                gfork.obj['Sequence'] = seq
                if not gobj.getBool('FullForce'): filelist.append(gfork.gFile('ALN/','orthaln.fas'))
            filelist.append(None)   # Indicates that previous files not found!
            ### ~ [2] Work through files and load if found, else use Gopher if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gobj.getBool('FullForce'): self.deBug(filelist)
            for afile in filelist:
                ## ~ [2a] Check for existing file and break loop if found ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if afile and os.path.exists(afile): break    # Existing alignment File found
                elif afile: continue
                ## ~ [2b] Call Gopher to make alignment if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if usegopher and seq not in self.list['GopherRun'] and self.opt['UseGopher']:  
                    #self.deBug('Run GOPHER in %s' % self.info['GopherDir'])
                    self.printLog('\n#GOPHER','Running GOPHER on %s' % seq.shortName(),log=log,screen=screen)
                    try:    #!# Add log.silent() method? #!#
                        gobj.gopherDir(seq)
                        self.list['GopherRun'].append(seq)
                        gfork.obj['BLAST'] = gobj.obj['BLAST']
                        gfork.info['Name'] = seq.shortName()
                        gfork.obj['Sequence'] = seq
                        gfork.obj['Parent'] = gobj
                        gfork.run('orthalign')    
                    except:
                        self.errorLog('Problem with Gopher run!')
                        return None
                    afile = gfork.gFile('ALN/','orthaln.fas')
                    if not os.path.exists(afile): afile = None
                ## ~ [2c] Abandon due to failure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not afile: return None
            if alnfile: return afile
            
            ### ~ [3] Load Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            v =  self.log.stat['Verbose']
            self.log.stat['Verbose'] = v - 1
            try:
                alncmd = ['seqin=None','query=%s' % seq.shortName(),'accnr=F','seqnr=F','autofilter=F','align=T','gnspacc=F'] 
                aln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
                aln.loadSeqs(seqfile=afile,seqtype='Protein',aln=True,nodup=None)
            except:
                self.log.stat['Verbose'] = v
                raise
            self.log.stat['Verbose'] = v 
            ## ~ [3a] Check Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = aln.obj['QuerySeq']
            if not qry:
                if aln.querySeq(query=seq.info['AccNum']): qry = aln.obj['QuerySeq']
                else:
                    self.printLog('#ALN','Problem finding %s in %s.' % (seq.shortName(),afile),screen=False,log=log)
                    return None
            if qry.aaLen() != seq.aaLen() and not gobj.getBool('FullForce') and usegopher:  # Problem with query sequence lengths. (Maybe an updated query.) Replace GOPHER alignment
                self.printLog('#ERR','Query %s (%d aa) is a different length (%d aa) in %s. Will try replacing with GOPHER.' % (seq.shortName(),seq.aaLen(),qry.aaLen(),afile),screen=True,log=log)
                gcmd = self.cmd_list[0:]; usealn = self.getBool('UseAln')
                gobj.setBool({'FullForce':True}); self.cmd_list.append('fullforce=T')
                self.setBool({'UseAln':False})  # Ignore existing alignment and replace!
                newaln = self.loadOrthAln(seq,usegopher,log,screen)
                gobj.setBool({'FullForce':False}); self.cmd_list = gcmd
                self.setBool({'UseAln':usealn}) # Return to use of existing alignments
                return newaln
            elif qry.aaLen() != seq.aaLen():
                self.printLog('#ERR','Query %s (%d aa) is a different length (%d aa) in %s.' % (seq.shortName(),seq.aaLen(),qry.aaLen(),afile),screen=True,log=log)
                return None
            ## ~ [3b] Check sequences for gaps and Xs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['HomFilter']:
                aln.autoFilter(['gapfilter=F']+self.cmd_list+['accnr=F','seqnr=F'])
                aln.gapSeqFilter(relative='query')
            ## ~ [3c] Check Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if aln.seqNum() < 2:
                self.printLog('#ALN','Not enough sequences for %s in %s.' % (seq.shortName(),afile),screen=False,log=log)
                return None
            if aln._checkAln(aln=True,realign=True): return aln
            else: self.log.printLog('#ERR','%s not aligned!!!' % (afile),screen=screen)
            return None       
        except:
            self.errorLog('Something bad has happened in rje_motif_stats.loadOrthAln()')
            raise
#########################################################################################################################
    def findOccPos(self,Occ,qry,fudge=0):     ### Finds Motif Occurence in alignment
        '''
        Finds Motif Occurence in alignment.
        >> Occ = MotifOcc object
        >> qry = query Sequence object from alignment file
        >> fudge = amount to try shifting match to find occurrence is non-matching sequence
        << (start,end) = start and end position in aligment to allow sequence[start:end]
        '''
        try:
            ### ~ [1] Find Hit in Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (start,end) = (-1,qry.seqLen())             # Start and end positins of match *in alignment*
            qpos = Occ['Pos'] + fudge                   # Starting position of hit (from 1->L)
            qmatch = Occ['Match']
            qend = qpos + len(qmatch) - 1               # Ending position of hit (1->L)
            (r,a) = (0,0)                               # Counters for aln residues (r) and amino acid positions (a)
            while r < qry.seqLen():     # Keep looking
                if qry.info['Sequence'][r] != '-': a += 1   # Not a gap: increment a by 1
                if a == qpos and start < 0: start = r       # Start of match (not yet r+=1 because pos is 1->L)
                r += 1                                      # Move on through aligned sequences
                if a == qend:                               # End of match
                    end = r
                    break

            ### ~ [2] Assess whether hit is right! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Can this be replaced with the rje_sequence method? #!#
            amatch = string.replace(qry.info['Sequence'][start:end],'-','')
            #x#if amatch == qmatch: return (start,end)     # Everything is OK!
            if re.search('^%s$' % string.replace(qmatch,'X','\S'),amatch): return (start,end) 
            if fudge != 0: raise ValueError             # Problem: already fudging!

            ### ~ [3] Something is wrong! Try to find real match! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            etxt = 'Alignment sequence (%s) does not match occurrence (%s)' % (amatch,qmatch)
            ## Try to find match by moving start (using fudge) ##
            if self.stat['Interactive'] < 1 or rje.yesNo('%s. Try to find closest correct match?' % etxt):
                fudge = self.findFudge(string.replace(qry.info['Sequence'],'-',''),qmatch,Occ['Pos']-1)
                if fudge:
                    if self.stat['Verbose'] > 0: self.warnLog('%s in alignment differs from input: Fudged %s by %d aa!' % (qry.shortName(),qmatch,fudge),'fudge')
                    else: self.warnLog('%s in alignment differs from input: Fudged %s by %d aa!' % (qry.shortName(),qmatch,fudge),'fudge',screen=False)
                    return self.findOccPos(Occ,qry,fudge)
                self.log.errorLog('%s in alignment differs from input: Cannot find %s anywhere' % (qry.shortName(),qmatch),printerror=False)
                return (-1,-1)
            self.log.errorLog(etxt,printerror=False)
            return (-1,-1)
        except:
            self.log.errorLog('Something bad has happened in rje_motif_stats.findOccPos()')
            return (-1,-1)
#########################################################################################################################
    def findFudge(self,qryseq,match,pos):    ### Returns fudge to closest match
        '''
        >> qryseq: degapped query sequence
        >> match:match sequence 
        >> pos:position match is meant to be from 0 to L
        << fudge:int = amount to move pos to find match in qryseq. 0 = not there!
        '''
        ### No match! ###
        rmatch = string.replace(match,'X','\S')
        #x#if qryseq.find(match) < 0: return 0
        if not re.search(rmatch,qryseq): return 0
        f = 1
        while (pos - f) >= 0  or (pos + f) < len(qryseq):
            #if (pos - f) >= 0 and qryseq[(pos - f):].find(match) == 0: return -f
            #elif (pos + f) < len(qryseq) and qryseq[(pos + f):].find(match) == 0: return f
            if (pos - f) >= 0 and re.search('^%s' % rmatch,qryseq[(pos - f):]): return -f
            elif (pos + f) < len(qryseq) and re.search('^%s' % rmatch,qryseq[(pos + f):]): return f
            f += 1
        raise ValueError
#########################################################################################################################
    ### <5> ### Combine Occurrence Calculations for a motif                                                             #
#########################################################################################################################
    def combMotifOccStats(self,occlist,revlist=['Hyd']):    ### Combines mean and percentile stats for the Occurrences of a single Motif
        '''
        Combines mean and percentile stats for the Occurrences of a Motif. OccList should all be for same SLiM.
        >> occlist:list of MotifOcc objects to calculate stats for (must all have same Motif)
        >> revlist:list of stats that should be ordered from low(best) to high(worst) rather than the other way round
        '''
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not occlist: return          # No occurrences to combine
        SLiM = occlist[0]['Motif']     # Should all have this SLiM object

        ### ~ [2] Combine into Motif.stat dictionary matching self.list['Headers'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for stat in self.list['OccHeaders']:
            ## ~ [2a] Generate list of values for each occurrence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            occval = []     # List of given stat for each occurrence
            for occ in occlist:
                try: occval.append(float(occ[stat]))
                except: self.deBug('%s occurrence missing %s value' % (SLiM.info['Name'],stat))
            if not occval: continue     # Nothing to combine
            ## ~ [2b] Basic Mean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try: SLiM.stat['%s_mean' % stat] = float(sum(occval)) / len(occval)
            except: self.log.errorLog('Problem combining %s values for %s' % (stat,SLiM.info['Name']))
            ## ~ [2c] Percentiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.stat['Percentile'] <= 0: continue
            if stat in revlist: occval.sort(reverse=True)
            else: occval.sort()
            pc = 100.0
            ovx = len(occval)
            for i in range(ovx):
                ipc = 100.0 * float((ovx-1) - i) / (ovx-1)  # This position
                jpc = 100.0 * float(ovx - i) / (ovx-1)      # Next position
                if ipc == pc: SLiM.stat['%s_pc%d' % (stat,int(pc+0.5))] = occval[i]    # Exact percentile
                elif ipc > pc and jpc < pc and (i+1) < ovx: SLiM.stat['%s_pc%d' % (stat,int(pc+0.5))] = float(occval[i] + occval[i+1]) / 2
                elif (i+1) == ovx: SLiM.stat['%s_pc0' % stat] = occval[i]
                else: continue
                pc -= self.stat['Percentile']
#########################################################################################################################
    ### <6> ### Occurrence and Motif Filtering                                                                          #
#########################################################################################################################
    def setupFilters(self,slimheaders=[],occheaders=[]):    ### Sets up SLiM/Occ Filters
        '''
        Sets up SLiM/Occ Filters.
        >> slimheaders:list [] = List of SLiM Attribute headers that can be used as filters.
        >> occheaders:list [] = List of SLiM Occurrence Attribute headers that can be used as filters.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not slimheaders: slimheaders = self.list['Headers']
            if not occheaders: occheaders = self.list['OccHeaders']
            ### ~ [2] Setup Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['SlimFilter'] = rje_scoring.setupStatFilter(self,slimheaders,self.list['SlimFilter'])
            self.dict['OccFilter'] = rje_scoring.setupStatFilter(self,occheaders,self.list['OccFilter'])
        except: self.log.errorLog('Problem during rje_slimcalc.setupFilters()')
#########################################################################################################################
    ### <7> ### Protein alignments                                                                                      #
#########################################################################################################################
    def singleProteinAlignment(self,seq,occlist,alndir='',hitname='AccNum',usegopher=True,savefasta=True,wintuple=0):    ### Generates copies of protein alignments, with motif hits marked.
        '''
        Generates copies of protein alignment, with motif hits marked. Return SeqList object.
        >> seq:Sequence object for alignment
        >> occlist:list of MotifOcc objects for sequence
        >> alndir:str = Alignment directory for output
        >> hitname:str = Format of hitnames, used for naming files
        >> usegopher:boolean = whether to look to use Gopher if settings correct (set this False if done already)
        >> savefasta:boolean [True] = whether to save fasta file or simply return SeqList object alone.
        >> wintuple:int [0] = Add a "WinTuple" list to aln Object containing [WinStart,WinEnd,VariantAln] (A list!)
        '''
        try:### ~ [1] Try to find alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['UseAln']: aln = self.loadOrthAln(seq,usegopher=usegopher)
            if not self.opt['UseAln'] or not aln:
                alncmd = ['seqin=None','query=%s' % seq.shortName(),'accnr=F','seqnr=F','autofilter=F','align=T','gnspacc=F']
                aln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
                aln.seq = [seq]
                aln.obj['QuerySeq'] = seq

            ### ~ [2] Check Query and File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not aln.obj['QuerySeq'] and not aln.querySeq(query=seq.info['AccNum']):
                self.log.printLog('#ERR','Problem finding %s in %s.' % (self.info['Name'],file))
                return None
            qry = aln.obj['QuerySeq']

            ### ~ [3] Create new sequence containing motifs, mapped onto Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if wintuple < 0: wintuple = 0
            aln.list['WinTuple'] = []
            motifseq = ['-'] * qry.seqLen()
            for Occ in occlist:
                (start,end) = self.findOccPos(Occ,qry)
                r = i = 0
                v = Occ['Motif'].variant(Occ)
                if v[:1] == '^': v = v[1:]
                w = ''
                for a in qry.info['Sequence'][start:end]:
                    r += 1
                    if a == '-':
                        w += '-'
                        continue
                    if v[i] not in ['.','X']: motifseq[start+r-1] = v[i]    #x# Occ['Variant'][i]
                    w += v[i]
                    i += 1
                aln.list['WinTuple'].append([max(0,start-wintuple),min(qry.seqLen(),end+wintuple),w])
                aln.list['WinTuple'][-1].append(start-aln.list['WinTuple'][-1][0])
            aln._addSeq('Motifs',string.join(motifseq,''))
            aln.seq = aln.seq[-1:] + aln.seq[:-1]
            if savefasta:
                if qry.info.has_key(hitname): aln.saveFasta(seqfile='%s%s.proteinaln.fas' % (alndir,qry.info[hitname]),log=False)
                else: aln.saveFasta(seqfile='%s%s.proteinaln.fas' % (alndir,qry.shortName()),log=False)
            return aln

        except:            
            self.log.errorLog('Major problem with rje_slimcalc.proteinAlignments()')
            return None
#########################################################################################################################
    ### <8> ### Whole sequence stats for background                                                                     #
#########################################################################################################################
    def seqAlnConList(self,seq):    ### Looks for alignment and, if appropriate, calculate conservation stats.
        '''
        Looks for alignment and, if appropriate, calculate conservation stats. Returns a list of the pos/prop
        conservation score for each residue, incorporating sequence weighting.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Identify/Create and Load Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aln = self.loadOrthAln(seq)
            if not aln: # Alignment rejected during loadOrthAln. No conservation statistics here.
                self.log.printLog('#ALN','No alignment file found for %s.' % OccSeq.shortName())  #!#,screen=False)
                return [1.0] * seq.aaLen()     
            ## ~ [1b] Setup alignment attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = aln.obj['QuerySeq']
            seqs = aln.seq[0:]      # This is a list of Sequence objects, minus the query, which is stored in qry
            seqs.remove(qry)
            self.log.printLog('#ALN','Comparing %d homologues for %s.' % (len(seqs),qry.shortName()))  #!#,screen=False)   #!# Log? #!#
            if len(seqs) < self.stat['MinHom']:
                self.printLog('#HOM','Too few (< %d) homologues for %s.' % (self.stat['MinHom'],qry.shortName()))  #!#,screen=False)   #!# Log? #!#
                return [1.0] * seq.aaLen()     

            ### ~ [2] Calculate Global Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqid = {}      # Dictionary of global %ID vs query for each sequence (consweight<>0)
            seqwt = {}      # Dictionary of global %ID weighting for each sequence (consweight<>0)
            for seq in seqs:
                igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'],seq.info['Sequence'],nomatch=['X'])
                seqid[seq] = float(igedic['ID'][0]) / igedic['Len'][0]
                seqwt[seq] = seqid[seq] ** self.stat['ConsWeight']
            
            ### ~ [3] Calculate Conservation List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            posmatrix = self.dict['PosMatrix']
            conslist = []
            for i in range(qry.seqLen()):   # Each residue in turn
                qa = qry.info['Sequence'][i]
                if qa == '-': continue
                elif qa == 'X':
                    conslist.append(-1.0)   # Why -1?! Surely conservation should be 0.5 or something? Check later calc.
                    continue
                ## ~ [3a] Positional conservation score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                hitcon = {}     # Dictionary of {seq:conservation}
                for seq in seqs:
                    aa = seq.info['Sequence'][i]
                    if aa != 'X': hitcon[seq] = self.bestScore(aa,[qa],posmatrix)
                ## ~ [3b] Weight by distance? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if hitcon: conslist.append(self.consWeight(hitcon,seqwt))
                else: conslist.append(1.0)
            return conslist
        except: self.errorLog('Error in rje_slimcalc.seqAlnConList()',quitchoice=True); raise
        return [-1.0] * seq.aaLen()     
#########################################################################################################################
    def relConListFromSeq(self,seq,window=30,store=False,log=True,screen=True):   ### Looks for alignment and returns list of relative conservation per residue.
        '''
        Looks for alignment and returns list of relative conservation per residue, incorporating sequence weighting.
        >> seq:Sequence object to analyse
        >> window:int [30] = Size of +/- window for calculation
        >> store:bool [False] = Whether to store results in seq.list['Cons'] and seq.list['RelCons']
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rlcp = '%s%s.rlc.pickle' % (self.info['AlnDir'],seq.info['AccNum'])
            if self.opt['Webserver'] and os.path.exists(rlcp): return pickle.load(open(rlcp,'r'))
            ## ~ [1a] Identify/Create and Load Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aln = self.loadOrthAln(seq,log=log,screen=screen)
            if not aln:                     # Alignment rejected during loadOrthAln. No conservation statistics here.
                self.log.printLog('#ALN','No alignment file found for %s.' % seq.shortName(),log=log,screen=screen)
                return [0.0] * seq.aaLen()  # Relative alignments are centred around a mean of zero   
            ## ~ [1b] Setup alignment attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = aln.obj['QuerySeq']   # This is the query Sequence object
            qry.obj['Disorder'] = seq.obj['Disorder']   # For MegaSLiM compatibility
            seqs = aln.seq[0:]          # This is a list of Sequence objects, minus the query, which is stored in qry
            seqs.remove(qry)
            self.log.printLog('#ALN','Comparing %d homologues for %s.' % (len(seqs),qry.shortName()),log=log,screen=screen)
            if len(seqs) < self.stat['MinHom']:
                self.printLog('#HOM','Too few (< %d) homologues for %s.' % (self.stat['MinHom'],qry.shortName()),log=log,screen=screen)
                return [0.0] * seq.aaLen()     
            ### ~ [2] Calculate Global Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqid = {}      # Dictionary of global %ID vs query for each sequence (consweight<>0)
            seqwt = {}      # Dictionary of global %ID weighting for each sequence (consweight<>0)
            for aseq in seqs:
                if self.stat['ConsWeight'] != 0:
                    igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'],aseq.info['Sequence'],nomatch=['X'])
                    seqid[aseq] = float(igedic['ID'][0]) / igedic['Len'][0]
                    seqwt[aseq] = seqid[aseq] ** self.stat['ConsWeight']
                else: seqwt[aseq] = 1
            relcon = self.relConList(qry,seqs,seqwt,window,store)
            if store:
                seq.list['RelCons'] = relcon[0:]
                seq.list['Cons'] = qry.list.pop('Cons')
            return relcon                
        except: self.errorLog('Error in rje_slimcalc.relConListFromSeq()',quitchoice=True); raise
        return [0.0] * seq.aaLen()              # Return flat nothing if failure.
#########################################################################################################################
    def relConList(self,qry,seqs,seqwt,window=30,store=False):   ### Returns list of relative conservation per residue using qry and seqs.
        '''
        Returns list of relative conservation per residue, incorporating sequence weighting, using qry and seqs.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rlcp = '%s%s.rlc.pickle' % (self.info['AlnDir'],qry.info['AccNum'])
            if self.opt['Webserver'] and os.path.exists(rlcp): return pickle.load(open(rlcp,'r'))
            ### ~ [1] Calculate Conservation List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xse = -1                            # Entropy for use with 100% X columns
            conslist = []
            for i in range(qry.seqLen()):       # Each residue in turn
                qa = qry.info['Sequence'][i]
                if qa == '-': continue          # Not returning any score for gapped residues in QUERY
                ## ~ [1a] Weighted Shannon Entropy score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cfreq = {qa:1.0}                # Dictionary of {aa:column frequency}
                for seq in seqs:
                    aa = seq.info['Sequence'][i]
                    if aa not in cfreq: cfreq[aa] = 0.0
                    cfreq[aa] += seqwt[seq] # Weight presence of AA by sequence weighting
                ## ~ [1b] ~ Gap Penalty and Masked residues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gap_pen = 1
                if '-' in cfreq:
                    gap_pen = 1 - (cfreq['-'] / sum(cfreq.values()))    # Proportion of ungapped sequences
                    cfreq.pop('-')
                if 'X' in cfreq: cfreq.pop('X') # X's count towards ungapped positions but contribute nothing to entropy
                if qry.dna() and 'N' in cfreq: cfreq.pop('N')
                if not self.opt['RelGapPen']: gap_pen = 1
                ## ~ [1c] Calculate entropy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #x#self.opt['DeBug'] = True
                base = {True:4,False:20}[qry.dna()]
                cfreq = self.vneEigen(rje.dictFreq(cfreq,total=False))
                se = 1.0    
                for aa in cfreq:
                    if cfreq[aa] > 1e-10: se += cfreq[aa] * math.log(cfreq[aa],base)     # Shannon entropy is -ve
                #x#self.deBug('%d: %s = %.3f' % (base,cfreq,se))
                if not cfreq:
                    if xse < 0: xse = self.aaFreqEntropy(qry,seqs)
                    se = xse
                if se < 0.0: raise ValueError
                conslist.append(gap_pen * se)   # Final entropy, weighted by non-gap proportion
            ### ~ [2] Convert into relative conservation across window ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if store: qry.list['Cons'] = conslist[0:]
            relconlist = self.relWinCon(conslist,window,qry)
            if self.opt['Webserver']: pickle.dump(relconlist,open(rlcp,'w'))
            return relconlist 
        except: self.log.errorLog('Error in rje_slimcalc.relConList()',quitchoice=True)
        return [0.0] * seq.aaLen()              # Return flat nothing if failure.
#########################################################################################################################
    def relWinCon(self,conslist,window=30,qry=None):    ### Returns relative conservation (x-mean)/sd for windows
        '''Returns relative conservation (x-mean)/sd for +/- windows.'''
        try:### ~ [1] ~ Calculate relative conservation (x-mean)/sd for +/- windows ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            relcon = []
            for i in range(len(conslist)):      # Each residue in turn
                if qry:
                    win = []
                    for w in range(max(0,i-window),min(len(conslist),i+window+1)):
                        if qry.isDisordered(w) == qry.isDisordered(i): win.append(conslist[w])
                else: win = conslist[max(0,i-window):i+window+1]
                (wmean,wsd) = rje.meansd(win)
                if wsd >= 1e-6: relcon.append((conslist[i] - wmean) / wsd)
                else: relcon.append(0.0)
            return relcon
        except: self.log.errorLog('Error in rje_slimcalc.relCon()',quitchoice=True)
        return [0.0] * len(conslist)            # Return flat nothing if failure.
#########################################################################################################################
    def aaFreqEntropy(self,qry,seqs):   ### Returns entropy based on whole alignment AAFreqs
        '''Returns entropy based on whole alignment AAFreqs.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if qry.dna(): alph = rje_slim.default_nts[0:]
            else: alph = rje_slim.default_aas[0:]
            freq = {}
            for x in alph: freq[x] = 0.0
            ### ~ [2] Add Freqs from Seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in [qry] + seqs: freq = seq.aaFreq(freq,newkeys=False)  # Adds to aafreq dictionary
            freq = rje.dictFreq(freq,total=False)                           # Convert to frequencies
            ### ~ [3] Calculate entropy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            base = len(alph)
            se = 1.0    
            for aa in freq:
                if freq[aa] > 0.0: se += freq[aa] * math.log(freq[aa],base)        # Shannon entropy is -ve
            if se < 0.0: raise ValueError
            return se
        except: self.errorLog('Problem with slimcalc.aaFreqEntropy()')
        return 0.0
#########################################################################################################################
    ### <9> ### VNE Methods                                                                                             #
#########################################################################################################################
    def loadBLOSUM(self):   ### Loads BLOSUM matrix for VNE calculations
        '''Loads BLOSUM matrix for VNE calculations.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['BLOSUM'] = {}
            if self.info['VNEMatrix'].lower() in ['','none']: return
            blosum = self.loadFromFile(self.info['VNEMatrix'],chomplines=True)
            while blosum[0][0] == '#': blosum.pop(0)
            ### ~ [2] ~ Read in BLOSUM matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aalist = blosum.pop(0)
            while aalist[0] == ' ': aalist = aalist[1:]
            aalist = string.split(aalist)
            #x#self.deBug(aalist)
            for aa in aalist:
                self.dict['BLOSUM'][aa] = {}
                val = string.split(blosum.pop(0))
                if val[0] == '': val = val[1:]
                #x#self.deBug('%s::%s' % (aa,val))
                for j in range(len(aalist)): self.dict['BLOSUM'][aa][aalist[j]] = string.atof(val[j])
            #self.deBug(self.dict['BLOSUM'])
        except:
            self.errorLog('Problem loading BLOSUM matrix "%s". Will not use VNE.' % self.info['VNEMatrix'])
            self.dict['BLOSUM'] = {}
#########################################################################################################################
    def vneEigen(self,aafreq):  ### Converts amino acid frequencies into VNE eigen values
        '''Converts amino acid frequencies into VNE eigen values.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['BLOSUM'] or not aafreq: return aafreq
            dm = []
            for ai in rje.sortKeys(aafreq):
                dm.append([])
                for aj in rje.sortKeys(aafreq):
                    dm[-1].append(aafreq[ai]*self.dict['BLOSUM'][ai][aj])
                    #self.deBug('%s [%s] x %s [%s] = %s' % (ai,aafreq[ai],aj,self.dict['BLOSUM'][ai][aj],dm[-1][-1]))
            #self.deBug(dm)
            ### ~ [2] ~ Convert to Eigen values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ev = ned_eigenvalues.Eigenvalues().get_eigenvalues(dm)
            #self.deBug(ev)
            eigfreq = {}
            for ai in rje.sortKeys(aafreq): eigfreq[ai] = ev.pop(0)
            eigfreq = rje.dictFreq(eigfreq,total=False)
            #self.deBug('%s =>\n%s' % (aafreq,eigfreq))
            return eigfreq 
        except:
            self.errorLog('vneEigen error!')
            return aafreq
#########################################################################################################################
### End of SECTION II: SLiMCalc Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def seqConsList(seq,callobj=None):   ### Creates SLiMCalc object and calls seqAlnConList(seq) method
    '''
    Creates SLiMCalc object and calls seqAlnConList(seq) method.
    >> seq:Sequence object
    >> callobj:Object containing relevant log and cmd_list objects. If None, will use seq
    '''
    ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if not callobj: callobj = seq
    return SLiMCalc(callobj.log,callobj.cmd_list+['usealn=T','conscore=pos']).seqAlnConList(seq)
#########################################################################################################################
def runMain():
    ### Rest of Functionality... ###
    try: print('\n\n *** No standalone functionality! *** \n\n')
    ### End ###
    except SystemExit: return  # Fork exit etc.
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
