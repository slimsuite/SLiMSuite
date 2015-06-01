#!/usr/bin/python

# RJE SLiM Class and Methods Module
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       RJE_SLiM
Description:  Short Linear Motif class module
Version:      1.12.0
Last Edit:    30/05/15
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the new SLiM class, which replaces the old Motif class, for use with both SLiMFinder and
    SLiMSearch. In addition, this module encodes some general motif methods. Note that the new methods are not
    designed with Mass Spec data in mind and so some of the more complicated regexp designations for unknown amino acid
    order etc. have been dropped. Because the SLiM class explicitly deals with *short* linear motifs, wildcard gaps are
    capped at a max length of 9.

    The basic SLiM class stores its pattern in several forms:
    - info['Sequence'] stores the original pattern given to the Motif object
    - info['Slim'] stores the pattern as a SLiMFinder-style string of defined elements and wildcard spacers
    - dict['MM'] stores lists of Slim strings for each number of mismatches with flexible lengths enumerated. This is
      used for actual searches in SLiMSearch.
    - dict['Search'] stores the actual regular expression variants used for searching, which has a separate entry for
      each length variant - otherwise Python RegExp gets confused! Keys for this dictionary relate to the number of
      mismatches allowed in each variant and match dict['MM'].

    The following were previously used by the Motif class and may be revived for the new SLiM class if needed:    
    - list['Variants'] stores simple strings of all the basic variants - length and ambiguity - for indentifying the
      "best" variant for any given match

    The SLiM class is designed for use with the SLiMList class. When a SLiM is added to a SLiMList object, the
    SLiM.format() command is called, which generates the 'Slim' string. After this - assuming it is to be kept -
    SLiM.makeVariants() makes the 'Variants' list. If creating a motif object in another module, these method should be
    called before any sequence searching is performed. If mismatches are being used, the SLiM.misMatches() method must
    also be called.

    SLiM occurrences are stored in the dict['Occ'] attribute. The keys for this are Sequence objects and values are
    either a simple list of positions (1 to L) or a dictionary of attributes with positions as keys.

Commandline:
    These options should be listed in the docstring of the module using the motif class:
    - alphabet=LIST     : List of letters in alphabet of interest [AAs]
    - ambcut=X          : Cut-off for max number of choices in ambiguous position to be shown as variant (0=All) [10]
    - trimx=T/F         : Trims Xs from the ends of a motif [False]
    - dna=T/F           : Whether motifs should be considered as DNA motifs [False]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import math, os, string, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_scoring
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial working version.
    # 1.1 - Added DNA option.
    # 1.2 - Added "N of M or B" format options.
    # 1.3 - Fixed the terminal variant bug.
    # 1.4 - Added makeSLiM method for converting a list of instances into a regexp
    # 1.5 - Added method to report whether motif splitting is necessary.
    # 1.6 - Fixed splitting bug introduced by lower case motifs.
    # 1.7 - Fixed import slimFix(slim) error that was reporting slimProb().
    # 1.8 - Modified use of aa/dna defaults to (hopefully) not break when using extended alphabets.
    # 1.9 - Reinstated ambcut for slimToPattern()
    # 1.10.0 - Added varlength option to makeSlim() method.
    # 1.10.1 - Fixed varlength and terminal position compatibility.
    # 1.10.2 - Fixed issue of [] returns.
    # 1.10.3 - Fixed makeSlim bug with variable length wildcards at start of sequence.
    # 1.11.0 - Added splitMotif() function.
    # 1.12.0 - Added equiv to makeSlim() function.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add new %x either/or regular expression. E.g. %1K..P..P..%1K = K..P..P or P..P..K
    # [ ] : Add weighted IC
    # [ ] : Add capability to handle wildcard spacers > 9.
    '''
#########################################################################################################################
default_aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
default_nts = ['A', 'C', 'G', 'T']
dna_ambig = {'U':'T','N':'.','B':'GTC','D':'GAT','H':'ACT','V':'GCA',   # official IUB/IUPAC abbreviations:
             'Y':'CT',          # pyrimidine
             'R':'AG',          # purine
             'K':'GT',          # keto
             'M':'AC',          # amino
             'S':'GC',          # strong
             'W':'AT'}          # weak
even_aafreq = {'A': 0.05, 'C': 0.05, 'E': 0.05, 'D': 0.05, 'G': 0.05, 'F': 0.05, 'I': 0.05, 'H': 0.05, 'K': 0.05,
               'M': 0.05, 'L': 0.05, 'N': 0.05, 'Q': 0.05, 'P': 0.05, 'S': 0.05, 'R': 0.05, 'T': 0.05, 'W': 0.05,
               'V': 0.05, 'Y': 0.05}
even_ntfreq = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiM Class                                                                                              #
#########################################################################################################################
class SLiM(rje.RJE_Object):     
    '''
    Short Linear Motif class. Author: Rich Edwards (2007). Based on old rje_motif_V3.Motif object.

    Info:str
    - Name = Name of motif
    - Description = Description of motif
    - Sequence = *Original* pattern given to motif
    - Slim = Reformatted sequence
    
    Opt:boolean
    - DNA = Whether motifs should be considered as DNA motifs [False]
    - TrimX = Trims Xs from the ends of a motif

    Stat:numeric
    - AmbCut = Cut-off for max number of choices in ambiguous position to be shown as variant [10]
        For mismatches, this is the max number of choices for an ambiguity to be replaced with a mismatch wildcard
    - IC = Information Content of motif

    List:list
    - Alphabet = List of letters in alphabet of interest [AAs]

    Dict:dictionary    
    - MM = stores lists of Slim strings for each number of mismatches. This is used for actual searches in SLiMSearch.
    - Occ = Dictionary of occurrences {Sequence:[Pos/{Hit stats}]} (Cannot use pos as key due to variable wildcards
    - Search = stores the actual regular expression variants used for searching, which has a separate entry for
      each length variant - otherwise Python RegExp gets confused! Keys for this dictionary relate to the number of
      mismatches allowed in each variant and match dict['MM'].
      
    Obj:RJE_Objects
    - SLiMList = Parent SLiMList object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Sequence','Description','Slim']
        self.optlist = ['DNA','TrimX']
        self.statlist = ['AmbCut','IC']
        self.listlist = ['Alphabet']
        self.dictlist = ['MM','Search','Occ']
        self.objlist = ['SLiMList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'AmbCut':10})
        self.list['Alphabet'] = default_aas
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._cmdReadList(cmd,'opt',['DNA','TrimX'])
                self._cmdReadList(cmd,'int',['AmbCut'])
                self._cmdReadList(cmd,'list',['Alphabet'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Basic Attribute Methods                                                                                 #
#########################################################################################################################
    def slimPos(self): return slimPosFromCode(self.info['Slim'])  ### Length of motif in terms of defined positions
    def slimLen(self): return slimLenFromCode(self.info['Slim'])  ### Maximum length of the motif (all positions)
    def slimMinLen(self): return slimMinLenFromCode(self.info['Slim'])  ### Min. length of motif (all positions)
    def slimFix(self): return slimFixFromCode(self.info['Slim'])  ### Length of motif in terms of fixed positions
    def seqs(self): return self.dict['Occ'].keys()  ### Sequences containing motif
    def seqNum(self): return len(self.dict['Occ'])  ### No. sequences containing motif
    def slim(self,wings=False): ### SLiMFinder format i-X-j. Wings adds extra '-' bracketting for subsequence matching
        if wings: return '-%s-' % self.info['Slim'] 
        return self.info['Slim']        
    def pattern(self): return patternFromCode(self.info['Slim'])
    def ambCutPattern(self): return patternFromCode(self.info['Slim'],self.stat['AmbCut'])
    def varLength(self): return '{' in self.pattern()
#########################################################################################################################
    def weightedIC(self,aafreq):   ### Replaces IC with weighted IC given aafreq, and returns.
        '''
        Replaces IC with weighted IC given aafreq, and returns.
        >> aafreq:dict = AA frequency dictionary
        '''
        self.stat['IC'] = 0.0
        for el in string.split(self.info['Slim'],'-'): self.stat['IC'] += weightedElementIC(el,aafreq)
        return self.stat['IC']
#########################################################################################################################
    def occNum(self):   ### Returns total number of occurrences
        '''Returns total number of occurrences.'''
        ox = 0
        for seq in self.dict['Occ']: ox += len(self.dict['Occ'][seq])
        return ox
#########################################################################################################################
    def occList(self):  ### Returns simple list of occurrences
        '''Returns simple list of occurrences.'''
        occlist = []
        for seq in self.dict['Occ']: occlist += self.dict['Occ'][seq]
        return occlist
#########################################################################################################################
    def variant(self,Occ):  ### Returns actual variant used for match in Occ
        '''Returns actual variant used for match in Occ.'''
        try:
            pattern = ''
            slist = string.split(slimFromPattern(Occ['Variant']),'-')
            for i in range(0,len(slist),2):
                ## ~ Deal with defined position first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                a = string.replace(slist[i].upper(),'X','.')
                if len(a) == 1: pattern += a
                elif pattern[:1] == '^': pattern += Occ['Match'][len(pattern)-1]
                else: pattern += Occ['Match'][len(pattern)]
                ## ~ Deal with wildcard spacer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                x = '0'
                if (i+1) < len(slist): x = slist[i+1]
                pattern += '.' * int(x)
            return pattern
        except:
            self.errorLog('SLiM error returning Variant pattern for %s' % Occ['Variant'])
            #if self.opt['Test']: self.opt['DeBug'] = True
            #self.deBug(Occ)
            return Occ['Variant']
#########################################################################################################################
    ### <3> ### Motif Formatting Methods                                                                                #
#########################################################################################################################
    def format(self,reverse=False):     ### Formats motif to generate self.info['Slim']
        '''
        Formats motif to generate self.info['Slim'].
        >> reverse:boolean [False] = whether to reverse sequence.
        << returns True if successful, or False if reformatting fails.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['MM'] = {}
            self.dict['Search'] = {}
            ### ~ [2] Reformat sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Slim'] = slimFromPattern(self.info['Sequence'][0:],reverse,trimx=self.opt['TrimX'],motif=self,dna=self.opt['DNA'])
            if self.info['Slim'] == 'WildLenErr':
                self.progLog('\r#REM',' ' * 100)
                self.printLog('\r#REJECT','Rejected %s: Max wildcard length (9) exceeded.' % self.info['Name'])
                return False
            if reverse: self.info['Sequence'] = patternFromCode(self.info['Slim'])
            ### ~ [3] Information Content & End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.stat['IC'] = 0.0
            for el in string.split(self.info['Slim'],'-'): self.stat['IC'] += elementIC(el)
            return True
        except:
            self.log.errorLog('Error in %s format(%s)' % (self.info['Name'],self.info['Sequence']))     
            return False
#########################################################################################################################
    def wildScram(self):    ### Performs wildcard scrambling
        '''Performs wildcard scrambling.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            wild = False
            slim = []
            ### ~ [2] Reformat sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for el in string.split(self.info['Slim'],'-'):
                if wild and el == '0': el = '1'
                elif wild: el = '0'
                slim.append(el)
                wild = not wild
            self.info['Slim'] = string.join(slim,'-')
            self.info['Sequence'] = patternFromCode(self.info['Slim'])
        except: self.log.errorLog('Error in %s wildscram(%s)' % (self.info['Name'],self.info['Sequence']))     
#########################################################################################################################
    def misMatch(self,mismatch={}):     ### Makes self.dict['MM'] using mismatch dictionary
        '''
        Makes self.dict['MM'] using mismatch dictionary.
        >> mismatch:dictionary of {mm 'X':Y aa}
        '''
        try:### ~ [1] No mismatches - expand %XORs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slim = string.split(self.info['Slim'],'-')
            mm = [slim] # List of lists used to make basic SLiMs
            ## ~ [1a] Generate XOR dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xor = {}    # Dictionary of {ID number:list of indices}
            for el in slim:
                if el[0] == '%':
                    if el[1] not in xor: xor[el[1]] = []
                    xor[el[1]].append(slim.index(el))
            ## ~ [1b] Enumerate XOR dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for x in xor:
                xvar = []
                for i in xor[x]:
                    for base in mm:
                        xnew = base[0:]
                        xnew[i] = xnew[i][2:]   # This is the variant this time
                        for j in xor[x]:
                            if i != j: xnew[j] = 'X'    # This variant is not wanted this time
                        xvar.append(xnew)
                mm = xvar[0:]
            ## ~ [1c] Convert mm list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['MM'] = {0:[]}
            for newvar in mm:
                newslim = slimFromPattern(patternFromCode(string.join(newvar,'-')))
                if newslim not in self.dict['MM'][0]: self.dict['MM'][0].append(newslim)
                
            ### ~ [2] Add mismatches according to mismatch dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 1      # No. of mismatches
            while mx in mismatch and self.slimPos() >= mismatch[mx]:    # Potentially have that many mismatches
                self.dict['MM'][mx] = []
                mvar = []   # List of slim lists
                for prevar in self.dict['MM'][mx-1]:    # Take each previous variant and add one more mismatch
                    slim = string.split(prevar,'-')
                    for i in range(0,len(slim),2):      # Take each defined position in turn
                        if slim[i] == 'X': continue     # Mismatched in previous round
                        newvar = slim[0:]
                        nofm = rje.matchExp('<(\D+):(\d+):(\d+)>',slim[i])
                        nofmorb = rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',slim[i])
                        if nofm and nofm[1] == '0': continue     # Mismatched in previous round
                        if nofmorb and nofmorb[1] == '0': continue     # Mismatched in previous round
                        if nofm: newvar[i] = '<%s:%d:%s>' % (nofm[0],max(0,int(nofm[1])- 1),nofm[2])
                        elif nofmorb: newvar[i] = '<%s:%d:%s:%s>' % (nofmorb[0],max(0,int(nofmorb[1])- 1),nofmorb[2],nofmorb[3])
                        else: newvar[i] = 'X'
                        newslim = slimFromPattern(patternFromCode(string.join(newvar,'-')))
                        if newslim not in self.dict['MM'][mx]: self.dict['MM'][mx].append(newslim)
                mx += 1
        except: self.log.errorLog('Major problem during %s mismatch (%s)' % (self.info['Name'],mismatch))
#########################################################################################################################
    def searchDict(self):   ### Makes 'Search' dictionary from mismatch dictionary
        '''Makes 'Search' dictionary from mismatch dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['MM']: self.misMatch()
            self.dict['Search'] = {}
            ### ~ [2] Make Search dictionary, expanding nofm and length variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mx in rje.sortKeys(self.dict['MM']):
                self.dict['Search'][mx] = []
                for mslim in self.dict['MM'][mx]:
                    slim = string.split(mslim,'-')
                    bases = [[]]
                    for i in range(0,len(slim),2):
                        ## ~ [2a] Defined position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        newvar = []
                        nofm = rje.matchExp('<(\D+):(\d+):(\d+)>',slim[i])
                        nofmorb = rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',slim[i])
                        if nofm:
                            vlist = [nofm[0]] * int(nofm[1]) + ['.'] * (int(nofm[2])-int(nofm[1]))
                            for nvar in rje.listRearrange(vlist):
                                for base in bases: newvar.append(base+nvar)
                        elif nofmorb:
                            vlist = [nofmorb[0]] * int(nofmorb[1]) + [nofmorb[3]] * (int(nofmorb[2])-int(nofmorb[1]))
                            for nvar in rje.listRearrange(vlist):
                                for base in bases: newvar.append(base+nvar)                            
                        elif len(slim[i]) > 1:
                            for base in bases: newvar.append(base+['[%s]' % string.replace(slim[i],'X','.')])
                        else:
                            for base in bases: newvar.append(base+[string.replace(slim[i],'X','.')])
                        bases = newvar[0:]
                        ## ~ [2b] Wildcards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if (i+1) < len(slim):
                            newvar = []
                            for x in slim[i+1]:
                                for base in bases: newvar.append(base + ['.'] * int(x))
                            bases = newvar[0:]
                    ## ~ [2c] Convert to proper RegExp and add ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for var in bases:
                        newre = string.join(var,'')
                        if newre not in self.dict['Search'][mx]: self.dict['Search'][mx].append(newre)
                        #self.deBug('%s >> %s' % (self.info['Name'],self.dict['Search'][mx]))
        except: self.log.errorLog('Major problem during %s searchDict' % (self.info['Name']))
#########################################################################################################################
    def varList(self,aafreq={},var_ic={}):  ### Makes variant list for CompariMotif searches etc.
        '''Makes variant list for CompariMotif searches etc.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            varlist = []
            self.weightedIC(aafreq)
            self.searchDict()
            for searchvar in self.dict['Search'][0][0:]:
                addvar = []
                while searchvar:
                    if searchvar[0] == '[':
                        el = rje.strSort(searchvar[1:searchvar.find(']')])
                        searchvar = searchvar[searchvar.find(']')+1:]
                    else:
                        el = searchvar[:1]
                        searchvar = searchvar[1:]
                    if not var_ic.has_key(el):
                        var_ic[el] = weightedElementIC(el,aafreq)    ### Can add aafreq and wild_pen etc. here if desired.
                    addvar.append(el)
                varlist.append(addvar[0:])
            return varlist
        except: self.errorLog('Error in SLiM.varDict()'); return []
#########################################################################################################################
    ### <4> ### Motif Searching                                                                                         #
#########################################################################################################################
    def searchSequence(self,seq=None,sequence='',logtext=''):  ### Searches the given sequence for occurrences of self and hits 
        '''
        Searches the given sequence for occurrences of self and returns a list of hit dictionaries: Pos,Variant,Match
        >> seq:Sequence object [None] = add hits to occurrence dictionary with this as key
        >> sequence:str = sequence to be searched if no sequence object given
        >> logtext:str [''] = text to precede progress printing. If '', no progress printing!
        << hitlist:list of dictionaries with hit information: Pos,Variant,Match,ID,MisMatch
        '''
        try:### ~ [1] Setup search objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Search']: self.searchDict()
            hitlist = []    # List of hit dictionaries to return
            rawhits = []    # List of ['pos:match'] to remove redundancy when mismatches used
            if logtext:
                varx = 0
                for mm in rje.sortKeys(self.dict['Search']): varx += len(self.dict['Search'][mm])
                sx = 0.0
                self.log.printLog('\r#SEARCH','%s [0 mismatch]: %.1f%%' % (logtext,(sx/varx)),log=False,newline=False)
        
            ### ~ [2] Search with self.dict['Search'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq: sequence = seq.info['Sequence'][0:]
            for mm in rje.sortKeys(self.dict['Search']):
                # Each motvar is a fixed length variant. The matchExp should therefore uniquely match the first...
                # ... occurrence of motvar in the sequence given each time, allowing one to "walk along" the sequence
                for motvar in self.dict['Search'][mm]:
                    ## ~ [2a] Setup search regexp string for this search variant ~~~~~~~~~~~~~~~~~~ ##
                    searchvar = string.replace(motvar,'X','.')
                    if searchvar[0] == '^': searchvar = '^(' + searchvar[1:]
                    else: searchvar = '(' + searchvar[0:]
                    if searchvar[-1] == '$': searchvar = searchvar[:-1] + ')$' 
                    else: searchvar += ')'
                    ## ~ [2b] Search sequence, walking along as matches are found ~~~~~~~~~~~~~~~~~ ##
                    r = 0   # Starting point of search
                    while rje.matchExp(searchvar,sequence[r:]):           # Found a match
                        match = rje.matchExp(searchvar,sequence[r:])[0]
                        if searchvar[-1] == '$': r = len(sequence) - len(match) + 1     # Must be at end
                        elif searchvar[0] == '^': r = 0                                 # Must be at beginning
                        else: r = r + sequence[r:].find(match) + 1        # r is now start of next search
                        if '%d:%s' % (r, match) not in rawhits:     # New hit
                            rawhits.append('%d:%s' % (r, match))
                            hitlist.append({'Match':match,'SearchVar':searchvar,'Variant':patternFromCode(slimFromPattern(motvar),ambcut=self.stat['AmbCut']),
                                            'MisMatch':mm,'Pos':r})      # Pos is aa number 1->L
                        if searchvar[0] == '^' or searchvar[-1] == '$': break    # Only match once!
                        #self.deBug(match)
                        #self.deBug(hitlist[-1])
                    ## ~ [2c] Progress ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if logtext:
                        sx += 100.0
                        self.log.printLog('\r#SEARCH','%s [%d mismatch]: %.1f%%' % (logtext,mm,(sx/varx)),log=False,newline=False)

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq and hitlist:
                if seq not in self.dict['Occ']: self.dict['Occ'][seq] = []
                self.dict['Occ'][seq] += hitlist[0:]
            return hitlist
        except:
            self.log.errorLog('Problem during SLiM.searchSequence(%s)' % self.info['Name'])
            return hitlist
#########################################################################################################################
    def occFilter(self,occfilter):  ### Filters occurrences according to given OccFilter dictionary
        '''Filters occurrences according to given OccFilter dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            occdict = {}
            occlist = self.occList()
            #x#print occlist
            for occ in occlist: occdict[occlist.index(occ)] = occ
            ### ~ [2] Perform Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje_scoring.statFilter(self,occdict,occfilter)
            occlist = occdict.values()
            for seq in self.dict['Occ'].keys()[0:]:
                for occ in self.dict['Occ'][seq][0:]:
                    if occ not in occlist: self.dict['Occ'][seq].remove(occ)
                if not self.dict['Occ'][seq]: self.dict['Occ'].pop(seq)
        except:
            self.log.errorLog('Error during SLiM.occFilter()')
#########################################################################################################################
### End of SECTION II: SLiM Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: General SLiM Methods                                                                                   #
#########################################################################################################################
def expectString(_expect):  ### Returns formatted string for _expect value
    '''Returns formatted string for _expect value.'''
    try:
        if _expect >= 10: return '%.1f' % _expect
        elif _expect >= 0.1: return '%.2f' % _expect
        elif _expect >= 0.001: return '%.3f' % _expect
        else: return '%.2e' % _expect
    except:
        print _expect
        raise
#########################################################################################################################
def checkSlimFormat(sequence):  ### Checks whether given sequence is in correct "Slim" format
    '''Checks whether given sequence is in correct "Slim" format.'''
    if len(string.split(sequence,'-')) < 3: return False
    wild = False
    for el in string.split(sequence,'-'):
        if wild and rje.matchExp('(\D)',el): return False       # Should be numbers only
        elif not wild and not rje.matchExp('(\D)',el): return False       # Should have non-numbers
        wild = not wild
    return wild     # Should end wanting another wildcard position
#########################################################################################################################
def patternIC(pattern,aafreq={}): return slimIC(slimFromPattern(pattern),aafreq) ### Returns (weighted) IC given aafreq.
def slimIC(slim,aafreq={}): ### Returns (weighted) IC given aafreq.
    '''
    Returns (weighted) IC given aafreq.
    >> aafreq:dict = AA frequency dictionary
    '''
    if not checkSlimFormat(slim): slim = slimFromPattern(slim)
    ic = 0.0
    for el in string.split(slim,'-'):
        if aafreq: ic += weightedElementIC(el,aafreq)
        else: ic += elementIC(el)
    return ic
#########################################################################################################################
def elementIC(element=''):     ### Calculates simplified IC for a given pattern element
    '''
    Calculates the IC for a given pattern element. This is much simplified from the rje_Motif method as (a) amino acid
    frequencies are never used, and (b) only wildcards can have flexible lengths, which does not affect this score.
    >> element:str = part of motif pattern
    << returns calculated IC for element
    '''
    ### ~ [1] Wildcards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if element in ['.','X','x'] or not rje.matchExp('(\D)',element): return 0.0
    ### ~ [2] NofM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if rje.matchExp('^<(\D+):(\d+):(\d+)>',element):    # "n of m" format
        nofm = rje.matchExp('^<(\D+):(\d+):(\d+)>',element)
        m = string.atoi(nofm[1])
        n = string.atoi(nofm[2])
        return elementIC(nofm[0]) * m
    if rje.matchExp('^<(\D+):(\d+):(\d+):(\D+)>',element):    # "n of m" format
        nofm = rje.matchExp('^<(\D+):(\d+):(\d+):(\D+)>',element)
        m = string.atoi(nofm[1])
        n = string.atoi(nofm[2])
        return elementIC(nofm[0]) * m + elementIC(nofm[3]) * (n-m)
    ### ~ [3] Defined positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    return ((math.log(0.05,2) - math.log(1.0/len(element),2)) / math.log(0.05,2))
#########################################################################################################################
def weightedElementIC(element,aafreq):     ### Calculates simplified IC for a given pattern element
    '''
    Calculates the IC for a given pattern element, weighted by amino acid frequencies are never used. Only wildcards can
    have flexible lengths, which does not affect this score.
    >> element:str = part of motif pattern
    >> aafreq:dict = dictionary of amino acid (or DNA) frequencies: length of dictionary defines no. of possibilities
    << returns calculated IC for element
    '''
    ### ~ [1] Wildcards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if element in ['.','X','x'] or not rje.matchExp('(\D)',element): return 0.0
    ### ~ [2] NofM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if rje.matchExp('^<(\D+):(\d+):(\d+)>',element):    # "n of m" format
        nofm = rje.matchExp('^<(\D+):(\d+):(\d+)>',element)
        m = string.atoi(nofm[1])
        n = string.atoi(nofm[2])
        return weightedElementIC(nofm[0],aafreq) * m
    ## ~ [2a] NofMorB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if rje.matchExp('^<(\D+):(\d+):(\d+):(\D+)>',element):    # "n of m or b" format
        nofm = rje.matchExp('^<(\D+):(\d+):(\d+):(\D+)>',element)
        m = string.atoi(nofm[1])
        n = string.atoi(nofm[2])
        return weightedElementIC(nofm[0],aafreq) * m + weightedElementIC(nofm[3],aafreq) * (n-m)
    ### ~ [3] Defined positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    evenfreq = 1.0 / len(aafreq)
    elfreq = 0.0
    try:
        for a in element: elfreq += aafreq[a]
    except: elfreq = evenfreq   # Allows terminal characters etc to get IC of 1.0
    return math.log(elfreq,2) / math.log(evenfreq,2)
#########################################################################################################################
def patternFromCode(slim,ambcut=20,dna=False):  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
    '''Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG).'''
    try:
        defaults = {False:default_aas,True:default_nts}[dna]
        pattern = ''
        if not slim: return pattern
        slist = string.split(slim,'-')
        for i in range(0,len(slist),2):
            ## ~ Deal with defined position first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            a = string.replace(slist[i].upper(),'X','.')
            if a[0] == '%':
                pattern += a[:2]
                a = a[2:]
            nofm = rje.matchExp('<(\D+):(\d+):(\d+)>',a)
            nofmorb = rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',a)
            if len(a) == 1: pattern += a
            elif rje.matchExp('<\D+:0:(\d+)>',a): pattern += '.' * int(rje.matchExp('<\D+:0:(\d+)>',a)[0])
            elif rje.matchExp('<\D+:0:(\d+):(\D+)>',a): pattern += rje.matchExp('<\D+:0:(\d+):(\D+)>',a)[1] * int(rje.matchExp('<\D+:0:(\d+):(\D+)>',a)[0])
            elif nofm: pattern += '<%s:%s:%s>' % (ambigPrint(nofm[0],defaults),nofm[1],nofm[2])
            elif nofmorb: pattern += '<%s:%s:%s:%s>' % (ambigPrint(nofmorb[0],defaults),nofmorb[1],nofmorb[2],ambigPrint(nofmorb[3],defaults))
            elif a[0] == '<': pattern += a
            elif a[0] == '(': pattern += a
            else: pattern += ambigPrint(a,defaults,ambcut)
            #x#elif ambcut and ambcut < len(a): pattern += 'X'     # Keep X for ambcut-altered position
            #elif len(a) > (len(defaults)+1) / 2:                # Inverse
            #    newaa = defaults[0:]
            #    for aa in a: newaa.remove(aa)
            #    pattern += '[^%s]' % string.join(newaa,'')
            #else: pattern += '[%s]' % a
            ## ~ Deal with wildcard spacer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            x = '0'
            if (i+1) < len(slist): x = slist[i+1]
            if len(x) == 1: pattern += '.' * int(x)
            else: pattern += '.{%d,%d}' % (int(x[0]),int(x[-1]))
        return pattern
    except:
        print slim, pattern
        raise
#########################################################################################################################
def prestoFromCode(slim):   ### Returns old PRESTO list from new SLiM code
    '''Returns old PRESTO list from new SLiM code. (Cannot have variable wildcard - will use longer.'''
    try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        presto = []
        slist = string.split(slim,'-')
        ### ~ [2] Extend PRESTO list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while slist:
            presto.append(slist.pop(0))
            if slist: presto += ['X'] * string.atoi(slist.pop(0)[-1])
        return presto
    except:
        print slim, presto
        raise
#########################################################################################################################
def convertDNA(slimcode):   ### Splits, converts to DNA and recombines
    '''Splits, converts to DNA and recombines.'''
    slim = []
    for s in string.split(slimcode,'-'):
        if s in dna_ambig: s = dna_ambig[s]
        slim.append(s)
    return string.join(slim,'-')
#########################################################################################################################
def needToSplitPattern(inseq):  ### Checks whether splitting needed to reformat motif to slim code
    '''
    Checks whether splitting needed to reformat motif to slim code.
    >> inseq:str = Sequence to reformat. Can be slimcode or pattern (inc. PROSITE)
    << returns True if split looks necessary.
    '''
    ### ~ [1] Setup sequence for reformatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    inseq = string.replace(inseq.upper(),'X','.')   # For clarity, all wildcards are . (and all upper case)
    inseq = string.replace(inseq,'[A-Z]','.')       # Keep all wildcards as .s
    inseq = string.join(string.split(inseq,'-'),'') # Join prosite-format
    if inseq.find('|') > 0: return True
    inseq = stripMotifBrackets(inseq)
    if rje.matchExp('([\[A-Za-z\^\]]+)\{(\d+),(\d+)\}',inseq): return True
    return False
#########################################################################################################################
def slimFromPattern(inseq,reverse=False,trimx=False,motif=None,dna=False):     ### Formats motif to generate slim code
    '''
    Formats motif to generate self.info['Slim'].
    >> inseq:str = Sequence to reformat. Can be slimcode or pattern (inc. PROSITE)
    >> reverse:boolean [False] = whether to reverse sequence.
    >> trimx:boolean [False] = whether to remove leading and trailing wildcards
    >> dna:boolean [False] = whether motif is a DNA motif
    >> motif:Motif object
    << returns slim code.
    '''
    ### ~ [1] Setup sequence for reformatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    defaults = {False:default_aas,True:default_nts}[dna]
    ## ~ [1a] Check whether the format is correct already ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if checkSlimFormat(inseq):  # Already formatted! Check aa ordering only.
        slim = []
        for el in string.split(inseq,'-'): slim.append(rje.strSort(el))
        if dna: return convertDNA(string.join(slim,'-'))
        else: return string.join(slim,'-')

    ## ~ [1b] Setup sequence string for reformatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    inseq = string.replace(inseq.upper(),'X','.')   # For clarity, all wildcards are . (and all upper case)
    if dna:
        for n in dna_ambig:     # Keep all wildcards as .s, Convert RNA to DNA, Convert DNA ambiguities
            a = dna_ambig[n]
            if len(a) == 1: inseq = string.replace(inseq,n,a)   # Convert RNA & wildcards
            else: inseq = string.replace(inseq,n,'[%s]' % a)   # Convert DNA ambiguities 
    inseq = string.replace(inseq,'[A-Z]','.')       # Keep all wildcards as .s
    #x#for unwelcome in ')(': inseq = string.replace(inseq,unwelcome,'')   # Get rid of these
    inseq = string.join(string.split(inseq,'-'),'') # Join prosite-format
    inseq = stripMotifBrackets(inseq)
    
    ### ~ [2] Reformat sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    slim = []       # Stores elements that will ultimately be joined to make new info['Slim']
    xor = False     # Whether previous element was a new and/or format
    w = (0,0)       # Min,Max wildcard
    while inseq:
        ## ~ [2z] Special convertion of variable non-wildcard to wildcard ~~~~~~~~~~~~~~~~~ ##
        if rje.matchExp('^((\[[A-Z]+\])\{(\d+),(\d+)\})',inseq[0:]) or rje.matchExp('^(([A-Z])\{(\d+),(\d+)\})',inseq[0:]):
            if rje.matchExp('^((\[[A-Z]+\])\{(\d+),(\d+)\})',inseq[0:]):
                (full,aa,mn,mx) = rje.matchExp('^((\[[A-Z]+\])\{(\d+),(\d+)\})',inseq[0:])
            else: (full,aa,mn,mx) = rje.matchExp('^(([A-Z])\{(\d+),(\d+)\})',inseq[0:])
            mn = string.atoi(mn)
            mx = string.atoi(mx)
            if mn == mx: newfull = aa * mx
            elif mn < 1: newfull = '.{%d,%d}' % (mn,mx)
            else: newfull = '%s.{%d,%d}' % (aa,mn-1,mx-1)
            inseq = newfull + inseq[len(full):]
            if motif:
                try: motif.log.printLog('\r#FORMAT','"%s" replaced with "%s" for Motif %s            ' % (full,newfull,motif.info['Name']))
                except: pass
            continue                              
        ## ~ [2a] Wildcard ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if inseq[0] == '.':     # Wildcard
            if rje.matchExp('^(\.\{(\d+),(\d+)\})',inseq[0:]):      # .{m,n}
                varwild = rje.matchExp('^(\.\{(\d+),(\d+)\})',inseq[0:])
                w = (w[0] + int(varwild[1]), w[1] + int(varwild[2]))
                inseq = inseq[len(varwild[0]):]
            elif rje.matchExp('^(\.\{(\d+)\})',inseq[0:]):      # .{m}
                varwild = rje.matchExp('^(\.\{(\d+)\})',inseq[0:])
                w = (w[0] + int(varwild[1]), w[1] + int(varwild[1]))
                inseq = inseq[len(varwild[0]):]
            else:
                w = (w[0]+1,w[1]+1)
                inseq = inseq[1:]
            continue
        if not slim and trimx:  w = (0,0)   # Lose N-term wildcards
        elif w[1] > 9:                      # Wildcard spacer too long
            try: motif.printLog('\r#FORMAT','Motif %s contains wildcard spacer > 9 residues: REJECTED.' % (motif.info['Name']))
            except: pass
            return 'WildLenErr'
        elif w[1] > 0 and not slim:         # Add wildcard spacer
            slim.append('X')    
            w = (max(0,w[0]-1),w[1]-1)
        if (w[1] > 0 or slim) and not xor:
            slim.append('')    
            for i in range(w[0],w[1]+1): slim[-1] += str(i)
        w = (0,0)
        ## ~ [2b] New and/or format %x ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if inseq[0] == '%':
            slim.append(inseq[:2])
            inseq = inseq[2:]
            xor = True
        ## ~ [2c] Fixed amino acid position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif inseq[0] in '^$' + string.uppercase:
            if xor: slim[-1] = slim[-1] + inseq[0]
            else: slim.append(inseq[0])
            inseq = inseq[1:]
            xor = False
        ## ~ [2d] Ambiguous amino acid position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif inseq[0] == '[':
            aa = inseq[1:inseq.find(']')]
            if aa[0] == '^':    # Inverse
                newaa = defaults[0:]
                for a in aa[1:]: newaa.remove(a)
                aa = string.join(newaa,'')
            if xor: slim[-1] = slim[-1] + aa
            else: slim.append(aa)
            inseq = inseq[inseq.find(']')+1:]
            xor = False
        ## ~ [2e] Special n of m position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif inseq[0] == '<':
            nm = inseq[:inseq.find('>')+1]      
            nofm = rje.matchExp('<(\D+):(\d+):(\d+)>',nm)
            if not nofm: nofm = rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',nm)
            if not nofm: print 'NofM ("<") Problem: %s' % nm; raise ValueError
            #!# Reformat n of m ambiguities #!#
            m1 = ambigPos(nofm[0],defaults)
            if len(nofm) == 4:
                b1 = ambigPos(nofm[3],defaults)
                nm = '<%s:%s:%s:%s>' % (m1,nofm[1],nofm[2],b1)
            else: nm = '<%s:%s:%s>' % (m1,nofm[1],nofm[2])
            if xor: slim[-1] = slim[-1] + nm
            else: slim.append(nm)
            inseq = inseq[inseq.find('>')+1:]
            xor = False
        ## ~ [2f] Either/or options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif inseq[0] == '(':
            eo = inseq[:inseq.find(')')+1]
            if eo.count('(') > 1 or eo.count(')') > 1: raise ValueError
            #!!#
        ## ~ [2g] Multiple characters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif rje.matchExp('^(\{(\d+)\})',inseq):
            mult = rje.matchExp('^(\{(\d+)\})',inseq)
            for m in range(int(mult[1]) - 1): slim += slim[-2:]
            slim = slim[:-1]
            inseq = inseq[len(mult[0]):]
        else:
            print 'Problem: %s' % inseq
            raise ValueError
    ## ~ [2h] Final Wildcards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if w[1] > 0 and not trimx:
        w = (max(0,w[0]-1),w[1]-1)
        slim.append('')    
        for i in range(w[0],w[1]+1): slim[-1] += str(i)
        slim.append('X')    
    ## ~ [2i] Stick elements together ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if reverse:
        slim.reverse()
        if slim[-1] == '^': slim[-1] = '$'
        if slim[0] == '$': slim[0] = '^'
    #x#raw_input(string.join(slim,'-'))
    if dna: return convertDNA(string.join(slim,'-'))
    else: return string.join(slim,'-')
#########################################################################################################################
def stripMotifBrackets(pattern):    ### Strips unneccessary brackets from motifs
    '''Strips unneccessary brackets from motifs.'''
    ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    stripped = ''
    base = pattern
    ### ~ [1] Cycle and split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    while base and base.find('|') > 0:
        i = base.find('|')
        ## ~ [1a] Find selection spot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        pre = post = 0  # Number of pre/post parentheses
        x = y = i       # x & y = ends of selection to vary
        while pre < 1:
            x -= 1
            if x < 0: raise ValueError
            if base[x] == '(': pre += 1
            if base[x] == ')': pre -= 1
            #self.bugPrint('%s:%s:%s - %s>>%s' % (x,i,y,base[x],base[y]))
        while post < 1:
            y += 1
            if y >= len(base): raise ValueError
            if base[y] == '(': post -= 1
            if base[y] == ')': post += 1
            #self.bugPrint('%s:%s:%s - %s>>%s' % (x,i,y,base[x],base[y]))
        ## ~ [1b] Replace swap region with variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        #self.bugPrint('%s >> %s & %s' % (base,base[x+1:i],base[i+1:y]))
        varlist = []
        for var in string.split(base[x+1:y],'|'):
            if var[0] == '(' and var[-1] == ')': var = var[1:-1]
            varlist.append(var)
        stripped += '%s(%s)' % (rje.stringStrip(base[:x+1],'()'),string.join(varlist,'|'))
        base = base[y+1:]
        #var1 = base[x+1:i]; var2 = base[i+1:y]
        #if var1[0] == '(' and var1[-1] == ')': var1 = var1[1:-1]
        #if var2[0] == '(' and var2[-1] == ')': var2 = var2[1:-1]
        #stripped += '(%s|%s)' % (var1,var2)
        #self.deBug('%s >> %s' % (bases,newmotifs))
    return stripped + rje.stringStrip(base,'()')
#########################################################################################################################
def ambigPos(ambig,defaults,ambcut=20):   ### Returns formatted ambiguities.
    '''Returns formatted ambiguities.'''
    if len(ambig) == 1: return ambig
    if ambig[0] == '[' and ambig[-1] == ']': ambig = ambig[1:-1]
    if ambig[0] == '^':
        alist = defaults[0:]
        for x in ambig[1:]:
            if x in alist: alist.remove(x)  # This will ignore odd alphabets, unfortunately
    else:
        alist = []
        for x in ambig[0:]: alist.append(x)
    alist.sort()
    if len(alist) > ambcut: return 'X'
    return '[%s]' % string.join(alist,'')
#########################################################################################################################
def ambigPrint(ambig,defaults,ambcut=20):   ### Returns formatted ambiguities.
    '''Returns formatted ambiguities.'''
    a = ambigPos(ambig,defaults,ambcut)    # Convert to all positive
    if len(a) == 1: return a
    a = a[1:-1]
    try:
        if len(a) > (len(defaults)+1) / 2:                # Inverse
            newaa = defaults[0:]
            for aa in a: newaa.remove(aa)
            return '[^%s]' % string.join(newaa,'')
    except: pass    #!# Odd alphabet mucking things up. Might need more robust fix!
    return '[%s]' % a
#########################################################################################################################
def slimPosFromCode(slim):  ### Returns the number of defined positions in a slim
    '''Returns the number of positions in a slim.'''
    px = 0
    for el in string.split(slim,'-'):
        if not rje.matchExp('(\D)',el): continue    # Wildcard
        elif rje.matchExp('<\D+:(\d+):\d+>',el): px += int(rje.matchExp('<\D+:(\d+):\d+>',el)[0])   # N of M
        elif rje.matchExp('<\D+:(\d+):(\d+):\D+>',el): px += int(rje.matchExp('<\D+:\d+:(\d+):\D+>',el)[0])   # N of M or B
        elif el != 'X': px += 1     # Non-wildcard position
    return px
#########################################################################################################################
def slimPos(slim):  ### Returns the number of defined positions in a SLiM of any format
    return slimPosFromCode(slimFromPattern(slim))
#########################################################################################################################
def slimLenFromCode(slim): ### Returns maximum length of SLiM in slimcode format
    '''Returns maximum length of SLiM in slimcode format.'''
    px = 0
    for el in string.split(slim,'-'):
        if not rje.matchExp('(\D)',el): px += int(el[-1])    # Wildcard
        elif rje.matchExp('<\D+:(\d+):\d+>',el): px += int(rje.matchExp('<\D+:\d+:(\d+)>',el)[0])   # N of M
        elif rje.matchExp('<\D+:(\d+):\d+:\D+>',el): px += int(rje.matchExp('<\D+:\d+:(\d+):\D+>',el)[0])   # N of M or B
        else: px += 1     # Single position
    return px
#########################################################################################################################
def slimMinLenFromCode(slim): ### Returns maximum length of SLiM in slimcode format
    '''Returns maximum length of SLiM in slimcode format.'''
    px = 0
    for el in string.split(slim,'-'):
        if not rje.matchExp('(\D)',el): px += int(el[0])    # Wildcard
        elif rje.matchExp('<\D+:(\d+):\d+>',el): px += int(rje.matchExp('<\D+:\d+:(\d+)>',el)[0])   # N of M
        elif rje.matchExp('<\D+:(\d+):\d+:\D+>',el): px += int(rje.matchExp('<\D+:\d+:(\d+):\D+>',el)[0])   # N of M or B
        else: px += 1     # Single position
    return px
#########################################################################################################################
def slimLen(slim):  ### Returns length of SLiM of any format
    if slim in ['','-']: return 0
    return slimLenFromCode(slimFromPattern(slim))
#########################################################################################################################
def slimFixFromCode(slim):  ### Returns the number of fixed positions in a slim
    '''Returns the number of positions in a slim.'''
    px = 0
    for el in string.split(slim,'-'):
        if not rje.matchExp('(\D)',el): continue    # Wildcard
        elif rje.matchExp('<(\D+):(\d+):\d+>',el):  # N of M
            if len(rje.matchExp('<(\D+):(\d+):\d+>',el)[0]) == 1: px += int(rje.matchExp('<\D+:(\d+):\d+>',el)[0])
        elif rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',el):  # N of M or B
            if len(rje.matchExp('<(\D+):(\d+):\d+:\D+>',el)[0]) == 1: px += int(rje.matchExp('<\D+:(\d+):\d+\D+>',el)[0])
            if len(rje.matchExp('<\D+:\d+:\d+:(\D+)>',el)[0]) == 1:
                px += int(rje.matchExp('<\D+:\d+:(\d+)\D+>',el)[0]) - int(rje.matchExp('<\D+:(\d+):\d+\D+>',el)[0])
        elif el != 'X' and len(el) == 1: px += 1     # Non-wildcard position
    return px
#########################################################################################################################
def slimFix(slim):  ### Returns the number of fixed positions in a SLiM of any format
    return slimFixFromCode(slimFromPattern(slim))
#########################################################################################################################
def makeSlim(peptides,minseq=3,minfreq=0.75,maxaa=5,callobj=None,ignore='X-',varlength=False,equiv=[]):   ### Generates a regexp SLiM from a peptide list.
    '''
    Generates a regexp SLiM from a peptide list (no variable wildcards unless varlength=True).
    >> peptides:list of str = List of peptide sequences (aligned, no gaps)
    >> minseq:int [3] = Min. no. of sequences for an aa to be in
    >> minfreq:num [0.75] = Min. combined freq of accepted aa to avoid wildcard
    >> maxaa:int [5] = Max. no. different amino acids for one position
    >> ignore:str ['X-'] = Amino acid(s) to ignore. (If nucleotide, would be N)
    >> varlength:bool [False] = Whether to recognise gaps and generate variable-length motifs.
    >> equiv:list [] = List of AA equivalencies to use to extend positions.
    '''
    try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        minseq = max(1,minseq)
        if len(peptides) < minseq: return ''
        mincount = minfreq * len(peptides)
        poslist = []
        ### ~ [1] Generate lists of each position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for frag in peptides:
            for i in range(len(frag)):
                if i >= len(poslist): poslist.append([])
                if frag[i] not in ignore: poslist[i].append(frag[i])
                elif varlength and frag[i] == '-': poslist[i].append('-')
        ### ~ [2] Build motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        motif = ''
        for pos in poslist:
            if varlength: mincount = minfreq * (len(pos) - pos.count('-'))
            gapped = pos.count('-') > 0
            aadict = {}
            for a in pos: 
                aadict[a] = pos.count(a)
            aasum = 0
            for a in rje.sortKeys(aadict):
                if aadict[a] < minseq or a == '-': aadict.pop(a)
                else: aasum += aadict[a]
            #if len(aadict) > 1 and '$' in aadict:   # Messy end of peptides: stop and truncate!
            #    break
            ## Incorporate equiv ambiguity
            eqsum = {}; best = ''
            for eqset in equiv:
                if rje.listDifference(aadict.keys(),eqset):    ### Returns the elements of list1 that are not found in list 2
                    continue    # Position has aas not in eqset
                elif not rje.listDifference(eqset,aadict.keys()): continue  # No aas to add
                else:
                    eqsum[eqset] = aasum
                    for a in rje.listDifference(eqset,aadict.keys()): eqsum[eqset] += pos.count(a)
                    if eqsum[eqset] > aasum:
                        if not best or eqsum[eqset] > eqsum[best] or (eqsum[eqset] >= eqsum[best] and len(eqset) < len(best)): best = eqset
            if best:    # equiv extended
                for a in rje.listDifference(best,aadict.keys()):
                    if a in pos: aadict[a] = pos.count(a); aasum += pos.count(a)
            ## Define position
            if aasum < mincount: motif += '.'
            elif len(aadict) > maxaa: motif += '.'
            elif len(aadict) == 1: motif += rje.sortKeys(aadict)[0]
            elif aadict: motif += '[%s]' % string.join(rje.sortKeys(aadict),'')
            if varlength and gapped:
                if '^' in aadict:
                    if callobj: callobj.warnLog('Gapped positions (-) aligned with sequence start (^). Replace gap with X or re-align!')
                elif '$' in aadict:
                    if callobj: callobj.warnLog('Gapped positions (-) aligned with sequence end ($). Replace gap with X or re-align!')
                elif pos.count('-') or motif.endswith('.') >= minseq: motif += '{0,1}'
                else: callobj.warnLog('Gapped position below minseq threshold: ignored.')
        ### ~ [3] Compress motif to varlength ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if varlength:
            mpos = []
            while motif:
                if motif.startswith('['): mpos.append(motif[:motif.find(']')+1]); motif = motif[motif.find(']')+1:]
                elif motif.startswith('{'): mpos.append(motif[:motif.find('}')+1]); motif = motif[motif.find('}')+1:]
                else: mpos.append(motif[:1]); motif = motif[1:]
            i = 1
            callobj.debug(mpos)
            while i < (len(mpos)-1):
                if mpos[i].startswith('{') and mpos[i-1] == mpos[i+1]:     # Try to combine adjacent positions
                    (m1,n1) = rje.matchExp('{(\d+),(\d+)}',mpos[i])
                    m1 = int(m1); n1 = int(n1)
                    if i < (len(mpos)-2) and mpos[i+2].startswith('{'):     # Combine two varlength
                        (m2,n2) = rje.matchExp('{(\d+),(\d+)}',mpos[i+2])
                        m2 = int(m2); n2 = int(n2)
                        mpos[i+2] = '{%d,%d}' % (m1+m2,n1+n2)
                        mpos = mpos[:i] + mpos[i+2:]
                    else:   # Combine varlength with neighbour
                        mpos[i] = '{%d,%d}' % (m1+1,n1+1)
                        mpos = mpos[:i+1] + mpos[i+2:]
                elif mpos[i] == mpos[i-1] == '.' and mpos[i+1].startswith('{'):
                    (m1,n1) = rje.matchExp('{(\d+),(\d+)}',mpos[i+1])
                    m1 = int(m1); n1 = int(n1)
                    mpos[i] = '{%d,%d}' % (m1+1,n1+1)
                    mpos = mpos[:i+1] + mpos[i+2:]
                else: i += 1
            if len(mpos) > 1 and mpos[1].startswith('{') and mpos[0] == '.': mpos = mpos[2:]
            if mpos and mpos[-1].startswith('{') and mpos[-2] == '.': mpos = mpos[:-2]
            motif = string.join(mpos,'')
        ### ~ [4] Trim end of motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while motif[:1] == '.': motif = motif[1:]
        if motif[:1] == '{': motif = motif[motif.find('}')+1:]
        while motif[-1:] == '.': motif = motif[:-1]
        if len(motif) == motif.count('^') + motif.count('.') + motif.count('$'): return ''
        return motif
    except:
        if callobj: callobj.errorLog('Problem with rje_slim.makeSlim()')
        raise
#########################################################################################################################
def slimDefPos(slim,pos,seqcheck=None,callobj=None):     ### Returns the list of positions of defined points
    '''
    Returns the list of positions of defined points. Will return empty list if seqcheck given and does not match.
    Flexible-length wildcards will split into variants and combine positions.
    >> slim:str = SLiM in SLiMSuite code format.
    >> pos:int = Starting position (0 < L) of occurrence.
    >> seqcheck:str [None] = Full-length sequence (optional) to check matches against.
    >> callobj:obj [None] = Optional object for debugging etc.
    '''
    ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    slimpos = []
    slimvar = []
    slimel = string.split(slim,'-')
    slimvar.append(slimel.pop(0))
    while slimel:
        wild = slimel.pop(0)
        varcount = len(slimvar); vx = 0
        slimvar = slimvar * len(wild)
        for w in wild:
            for v in range(varcount):
                slimvar[vx] += '-%s' % w
                vx += 1
        aa = slimel.pop(0)
        for vx in range(len(slimvar)): slimvar[vx] += '-%s' % aa
    if callobj: callobj.debug(slim)
    ### ~ [1] ~ Generate position lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    for slimcode in slimvar:
        ## ~ [1a] ~ Check sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if seqcheck:
            if callobj: callobj.debug(seqcheck[pos:])
            if callobj: callobj.debug('^(%s)' % patternFromCode(slimcode))
            slimmatch = rje.matchExp('^(%s)' % patternFromCode(slimcode),seqcheck[pos:])
            if callobj: callobj.debug(slimmatch)
            if not slimmatch: continue
        ## ~ [1b] ~ Identify positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        px = 0
        slimel = string.split(slim,'-')[1:]
        slimpos.append(pos)
        while slimel:
            px += int(slimel.pop(0)) + 1
            slimpos.append(pos+px)
            slimel.pop(0)
    if callobj: callobj.debug(slimpos)
    return rje.sortUnique(slimpos)
#########################################################################################################################
def varWild(slim):  ### Whether the SLiM has variable-length wildcards
    '''Whether the SLiM has variable-length wildcards.'''
    slimel = string.split(slim,'-')[1:]
    while slimel:
        if len(slimel.pop(0)) > 1: return True
        slimel.pop(0)
    return False
#########################################################################################################################
def splitPattern(regex):     ### Splits complex motifs on "|" and returns list
    '''Splits complex motifs on "|" and adds each separately.'''
    try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        bases = []          # This is a list of basic strings to be added to with each variation
        newmotifs = [regex]   # New motif variants
        splitting = True
        ### ~ [1] Cycle and split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while splitting:
            bases = newmotifs[0:]
            newmotifs = []
            splitting = False
            for base in bases:
                i = base.find('|')
                if i < 0:
                    vardef = rje.matchExp('([\[A-Za-z\^\]]+)\{(\d+),(\d+)\}',base)
                    if not vardef: newmotifs.append(base); continue
                    splitting = True
                    (nonwild,m,n) = vardef
                    m = int(m)
                    n = int(n)
                    if nonwild[-1] == ']': nonwild = '[%s' % string.split(nonwild,'[')[-1]
                    else: nonwild = nonwild[-1]
                    if m == n: newmotifs.append(string.replace(base,'%s{%d,%d}' % (nonwild,m,n),nonwild*m))
                    else:
                        for x in range(m,n+1): newmotifs.append(string.replace(base,'%s{%d,%d}' % (nonwild,m,n),nonwild*x))
                    continue
                splitting = True
                pre = 0 # Number of pre parentheses - check X|Y split has surrounding parentheses (or add)
                x = i   # x = start of selection to vary
                while pre < 1:
                    x -= 1
                    if x < 0:   # The whole pattern is an X|Y split with no outside brackets.
                        base = '(%s)' % base
                        i += 1
                        break
                    if base[x] == '(': pre += 1
                    if base[x] == ')': pre -= 1
                ## ~ [1a] Find selection spot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pre = post = 0  # Number of pre/post parentheses
                x = y = i       # x & y = ends of selection to vary
                while pre < 1:
                    x -= 1
                    if x < 0:   # The whole pattern is an X|Y split with no outside brackets.
                        raise ValueError
                    if base[x] == '(': pre += 1
                    if base[x] == ')': pre -= 1
                while post < 1:
                    y += 1
                    if y >= len(base): raise ValueError
                    if base[y] == '(': post -= 1
                    if base[y] == ')': post += 1
                ## ~ [1b] Replace swap region with variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for var in (base[x+1:i],base[i+1:y])[0:]:
                    if var[0] == '(' and var[-1] == ')': var = var[1:-1]
                    newmotifs.append(base[:x] + var + base[y+1:])

        ### ~ [2] Clean up motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return rje.sortUnique(newmotifs)
    except: raise
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def test(motif,cmd_list=[],mm={}):    ### Temp test method
    '''Temp test method.'''
    slim = SLiM(cmd_list=cmd_list)
    slim.info['Name'] = slim.info['Sequence'] = motif
    slim.format()
    print slim.info
    slim.misMatch(mm)
    print slim.dict['MM']
    slim.searchDict()
    print slim.dict['Search']
    return slim
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: print 'This module is not for standalone running.'
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
