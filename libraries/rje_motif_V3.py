#!/usr/local/bin/python

# RJE Motif Class and Methods Module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_motif
Description:  Motif Class and Methods Module
Version:      3.1
Last Edit:    04/06/14
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the Motif class for use with both Slim Pickings and PRESTO, and associated methods. This basic
    Motif class stores its pattern in several forms:
    - info['Sequence'] stores the original pattern given to the Motif object
    - list['PRESTO'] stores the pattern in a list of PRESTO format elements, where each element is a discrete part of
      the motif pattern
    - list['Variants'] stores simple strings of all the basic variants - length and ambiguity - for indentifying the "best"
      variant for any given match
    - dict['Search'] stores the actual regular expression variants used for searching, which has a separate entry for
      each length variant - otherwise Python RegExp gets confused! Keys for this dictionary relate to the number of
      mismatches allowed in each variant.

    The Motif Class is designed for use with the MotifList class. When a motif is added to a MotifList object, the
    Motif.format() command is called, which generates the 'PRESTO' list. After this - assuming it is to be kept -
    Motif.makeVariants() makes the 'Variants' list. If creating a motif object in another module, these method should be
    called before any sequence searching is performed. If mismatches are being used, the Motif.misMatches() method must
    also be called.

Commandline:
    These options should be listed in the docstring of the module using the motif class:
    - alphabet=LIST     : List of letters in alphabet of interest [AAs]
    - ambcut=X          : Cut-off for max number of choices in ambiguous position to be shown as variant (0=All) [10]
    - trimx=T/F         : Trims Xs from the ends of a motif [False]

Uses general modules: copy, math, os, re, string, sys
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import math
import os
import re
import string
import sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 2.0 - Revised version based on RJE_MOTIF 1.1 and bits of PRESTO 1.8.
    # 2.1 - Added "n of m" style elements in the form <X:n:m>
    # 2.2 - AmbCut now limits the number of ambiguous choices when compare=T
    # 2.3 - Added Information Content Methods
    # 2.4 - Added defineMotif method
    # 2.5 - Added reformatMiniMotif method
    # 3.0 - Reworked for use with MotifList and MotifOcc objects. Add Expect dictionary.
    # 3.1 - Fixed minor code bugs.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] General tidy with updated SLiMSuite.
    '''
#########################################################################################################################
basic_aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
basic_aafreq = {'A': 0.05, 'C': 0.05, 'E': 0.05, 'D': 0.05, 'G': 0.05, 'F': 0.05, 'I': 0.05, 'H': 0.05, 'K': 0.05,
                'M': 0.05, 'L': 0.05, 'N': 0.05, 'Q': 0.05, 'P': 0.05, 'S': 0.05, 'R': 0.05, 'T': 0.05, 'W': 0.05,
                'V': 0.05, 'Y': 0.05}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Motif Class                                                                                             #
#########################################################################################################################
class Motif(rje.RJE_Object):     
    '''
    Motif Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of motif
    - Description = Description of motif
    - Sequence = *Original* pattern given to motif
    
    Opt:boolean
    - TrimX = Trims Xs from the ends of a motif
    - Compare = Compare the motifs from the motifs FILE with the searchdb FILE (or self if None) [False]
    - MatchIC = Use (and output) information content of matched regions to asses motif matches [True]
    - MotifIC = Output Information Content for motifs [False]

    Stat:numeric
    - AmbCut = Cut-off for max number of choices in ambiguous position to be shown as variant [10]
        For mismatches, this is the max number of choices for an ambiguity to be replaced with a mismatch wildcard
    - Length = Maximum length of the motif in terms of non-wildcard positions
    - MinLength = Minimum length of the motif in terms of non-wildcard positions
    - FixLength = Maximum Length of motif in terms of fixed positions
    - FullLength = Maximum length of the motif, including wildcard positions
    - IC = Information Content of motif
    - OccNum = Number of occurrences in search database
    - OccSeq = Number of different sequences it occurs in in search database

    List:list
    - Alphabet = List of letters in alphabet of interest
    - PRESTO = Presto format motifs are strings of elements separated by '-', where each element is:
        > a single AA letter
        > a wildcard 'X'
        > a choice of letters in the form [ABC] *** NB. an "except" [^ABC] is converted to inclusive ambiguity ***
        > a choice of combinations in the form (AB|CD)
        > a start ^ or end $ of sequence marker
        > may be combined with variable numbers of positions {m,n}
    - Variants = List of string variants lists, incorporating length variation and different combos.
        => This is primarily used to determine the best "variant" match for an actual match but also as the base for mismatches.

    Dict:dictionary
    - Expect = dictionary of {key:expected number of occurrences}, where key could be a filename or Sequence object etc.
    - ExpectMM = same as Expect but for each number of mismatches {key:{mm:expect}} - PRESTO only.
    - Search = dictionary of {no. mismatches:list of variant regexps to search}

    Obj:RJE_Objects
    - MotifList = "Parent" MotifList object - contains objects of use to Motif without need to duplicate
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Sequence','Description']
        self.optlist = ['TrimX','Compare','MatchIC','MotifIC']
        self.statlist = ['AmbCut','Length','MinLength','FullLength','FixLength','IC','OccNum','OccSeq']
        self.listlist = ['Alphabet','PRESTO','Variants']
        self.dictlist = ['Expect','ExpectMM','Search']
        self.objlist = ['MotifList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Pattern':''})
        self.setStat({'AmbCut':10,'OccNum':0,'OccSeq':0})
        self.list['Alphabet'] = basic_aas
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ### 
                self._cmdRead(cmd,type='list',att='Alphabet')
                self._cmdRead(cmd,type='int',att='AmbCut')
                self._cmdReadList(cmd,'opt',['TrimX','Compare','MatchIC','MotifIC'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Motif Formatting                                                                                        #
#########################################################################################################################
    def format(self,msmode=False,reverse=False):  ### Generates list.['PRESTO'] and list.['Variants'] from info['Sequence']
        '''
        Generates list.['PRESTO'] and list.['Variants'] from info['Sequence']. See docstring for details of PRESTO format.
        >> msmode:boolean [False] = whether to interpret motif as MSMS peptide sequencing.
        >> reverse:boolean [False] = whether to reverse sequence
        '''
        try:
            ### Setup Reformatting Tools ###
            _stage = 'Setup Tools'
            ## These RegExps are used to identify elements of the pattern from its sequence. ##
            ## They all match the start of the expression, as this is then chopped off to make the next element ##
            regexps = ['^([A-Z])',          # Single AA : A
                       '^(\[[A-Z]+\])',     # AA choice : [ABC]
                       '^(\[\^[A-Z]+\])',   # ELM AA aversion : [^ABC]
                       '^(\(([A-Z])\))',    # ELM single focus : (A)
                       '^(\((\[[A-Z]+\])\))',   # ELM ambiguous focus : ([AB])
                       '^(\([A-Z\|]+\))',       # Special combo : (AB|CD) or (AB)
                       '^(\^)',             # Start of sequence : ^
                       '^(\$)',             # End of sequence   : $
                       '^(\{(\d+)\})',      # Single number : {m}
                       '^(\{(\d+),(\d+)\})',    # Ambiguous numbers : {m,n}
                       '^(<(\D+):(\d+):(\d+)>)',    # "n of m" format
                       '^([\-\s])']       # Spacer : - or space
            ## This dictionary is used to put a text name to the element matched for easier processing ##
            regtype = {'^([A-Z])':'single',
                       '^(\[[A-Z]+\])':'choice',
                       '^(\[\^[A-Z]+\])':'not',
                       '^(\(([A-Z])\))':'focus',
                       '^(\((\[[A-Z]+\])\))':'amb_focus',
                       '^(\([A-Z\|]+\))':'combo',
                       '^(\^)':'nterm',
                       '^(\$)':'cterm',
                       '^(\{(\d+)\})': 'multiple',
                       '^(\{(\d+),(\d+)\})':'numbers',
                       '^(<(\D+):(\d+):(\d+)>)':'n_of_m',
                       '^([\-\s])':'spacer'}

            ### Setup Attributes ###
            _stage = 'Setup Attributes'
            self.list['PRESTO'] = []
            
            ### Setup Sequence ###
            _stage = 'Setup Sequence'
            inseq = self.info['Sequence'][0:]       # This is the string that will be matched and reduced during formatting
            inseq = string.replace(inseq.upper(),'.','X')   # For clarity, all wildcards are X (and all upper case)
            inseq = string.replace(inseq,'[A-Z]','X')       # Keep all wildcards as Xs
            inseq = string.replace(inseq,')',') ')  # Not sure why this is here but there must be a reason! (Old match?)
            inseq = string.replace(inseq,'/','|')   # This allows '/' to be used as "OR"   

            ### TrimX ###
            while self.opt['TrimX']:
                if inseq[0] == 'X':
                    inseq = inseq[1:]
                elif inseq[-1] == 'X':
                    inseq = inseq[:-1]
                elif inseq[0] == '{':
                    inseq = inseq[inseq.find('}')+1:]
                else:
                    break

            ### XPad ###
            inseq = 'X' * self.obj['MotifList'].stat['XPad'] + inseq + 'X' * self.obj['MotifList'].stat['XPad']

            ### Format Sequence ###
            _stage = 'Formatting'
            presto = []     # Temporary store of formatted motif
            while len(inseq) > 0:
                ## Identify appropriate regular expression match ##
                match = None
                for reg in regexps:
                    if rje.matchExp(reg,inseq):    
                        match = rje.matchExp(reg,inseq)     # match has the actual element matched
                        break                               # reg now has the regular expression that matched
                ## Modify inseq if match OK, else Error
                if match:
                    inseq = inseq[len(match[0]):]       # Reduce inseq ready for next match
                else:
                    self.log.errorLog('Unrecognised character " %s ". Check allowed formatting in manual.' % inseq[0],printerror=False)
                    return False
                ## Add element, reformatting as necessary ##
                rtype = regtype[reg]    # Type of element matched by regexp
                element = match[0]      # This is the element that will be added to PRESTO
                if rtype == 'spacer':   # Ignore spacers and continue
                    continue
                if rtype == 'n_of_m' and rje.matchExp('^<\[(\D+)\]:(\d+):(\d+)>',element):
                    element = '<%s:%s:%s>' % rje.matchExp('^<\[(\D+)\]:(\d+):(\d+)>',element)
                ## Reformat ELM-type RegExp ##
                if rtype == 'not':      # Reverse this into an inclusive ("standard") ambiguity
                    newchoice = ''
                    for aa in self.list['Alphabet']:
                        if element.find(aa) < 0:        # Not excluded
                            newchoice += aa
                    if len(newchoice) > 1:
                        element = '[%s]' % newchoice    # Replace element with new ambiguity
                        rtype = 'choice'
                    elif len(newchoice) == 1:
                        element = newchoice    # Replace element with new fixed position
                        rtype = 'single'
                    else:   # Removed all possibles!
                        self.log.errorLog('Exclusive pattern "%s" removed entire alphabet!' % element,printerror=False)
                        raise ValueError
                    match = (element)
                if rtype in ['focus','amb_focus']:      # Remove additional brackets identifying focal AA
                    #X#match = match[1:2]
                    element = match[1]  #X#[0]
                    if rtype == 'focus':
                        rtype = 'single'   # Single
                    else:
                        rtype = 'choice'   # Choice
                ## Reformat AmbCut to Wildcard if self.opt['Compare']
                #X#self.deBug('')
                if rtype == 'choice' and self.opt['Compare'] and self.stat['AmbCut'] > 0 and len(element) > (self.stat['AmbCut'] + 2):
                    self.log.printLog('\r#AMB','%s Ambiguity "%s" > ambcut (%d) => Changed to wildcard "X"' % (self.info['Name'],element,self.stat['AmbCut']))
                    #X#self.deBug('%s -> X' % element)
                    element = 'X'
                    rtype = 'single'
                ## Special Reformatting for MSMS mode ##
                if msmode:      ### Replace user-defined known ambiguities to avoid problems later
                    x = {}
                    for a in 'QKIL':    # All Is become Ls and Qs become Ks
                        x[a] = string.count(element,a)
                    if rtype in ['single','combo']:
                        if x['I']:
                            element = string.replace(element,'I','L')
                        if x['Q']:
                            element = string.replace(element,'Q','K')
                    elif rtype == 'choice':
                        if x['I'] and x['L']:
                            element = string.replace(element,'I','')
                        elif x['I']:
                            element = string.replace(element,'I','L')
                        if x['K'] and x['Q']:
                            element = string.replace(element,'Q','')
                        elif x['Q']:
                            element = string.replace(element,'Q','K')
                        if len(element) == 3:
                            element = element[1]
                            rtype = 'single'
                ## Reformat Combos with special MS Mode option ##
                if rtype == 'combo':
                    if element.find('|') < 0:   # No explicit option, so change (AB) to (AB|BA)
                        aas = rje.matchExp('\(([A-Z]+)\)',element)[0]
                        element = '(%s|%s)' % (aas,rje.strReverse(aas))
                    elif msmode:        # Complicated! Want all options (AB|CD|BA|DC)
                        aalist = []     # List of different AA options [AB,CD,BA,DC]
                        heart = element[1:-1]       # AAs without brackets
                        while heart.find('|') > 0:  # Options remain
                            hx = heart.find('|')    # Find first set of options
                            aas = heart[:hx]        # Pull out aas
                            aalist += rje.strRearrange(aas) # Add all possible orders of AAs to aalist
                            heart = heart[hx+1:]    # Redefine remainder of pattern to search
                        aalist += rje.strRearrange(heart)   # Add all possible combos of remaining AAs to aalist
                        element = '(%s)' % string.join(aalist,sep='|')  # Make new combo element with all options!
                ## Multiple occurrences of elements ##
                if rtype == 'multiple':     ## Convert {m} to {m,n} (where m=n)
                    rtype = 'numbers'
                    match = ('{%s,%s}' % (match[1],match[1]),match[1],match[1])
                    element = match[0]
                ## Add to PRESTO - special if numbers ##
                if rtype == 'numbers':   ### {m,n} format
                    if match[1] == match[2]:    # One number! => add multiple occurrences of element to PRESTO
                        presto = presto[:-1] + (presto[-1:] * string.atoi(match[1]))
                    else:
                        presto[-1] = presto[-1] + element
                else:   # Singles, NTerms and CTerms, n_of_m.
                    presto.append(element)
                
            ### Finish Off ###
            _stage = 'Finish'
            self.list['PRESTO'] = presto[0:]
            if reverse:
                self.list['PRESTO'].reverse()
                if presto[0] == '^':
                    self.list['PRESTO'][-1] = '$'
                if presto[-1] == '$':
                    self.list['PRESTO'][0] = '^'
                self.info['Sequence'] = string.replace(string.join(self.list['PRESTO'],''),'X','.')
            #X#print self.info['Name'], self.info['Sequence'], self.list['PRESTO']
            self._calculateLength()             # Calculates Length Statistics

            ### Information Content ###
            _stage = 'IC'
            presto = self.obj['MotifList']
            self.stat['IC'] = 0.0
            for el in self.list['PRESTO']:
                #X#self.deBug(el)
                if el not in presto.dict['ElementIC'].keys():
                    presto.dict['ElementIC'][el] = elementIC(el,callobj=self)
                #X#self.deBug(presto.dict['ElementIC'])
                self.stat['IC'] += presto.dict['ElementIC'][el]    # Can add aafreq later if desired
            return True

        except:
            self.log.errorLog('Error in %s format(%s) %s: ' % (self.info['Name'],self.info['Sequence'],_stage),quitchoice=True)     
            return False
#########################################################################################################################
    def _msVar(self,varlist): ### Returns MS-altered variant list
        '''
        Returns MS-altered variant list.
        >> varlist:list of variant regexp sequence strings
        '''
        try:
            newlist = []
            for var in varlist[0:]:
                ## Replace special MSMS ambiguities ##
                var = string.replace(var,'L','[IL]')
                var = string.replace(var,'K','[KQ]')
                var = string.replace(var,'B','[KR]')
                var = string.replace(var,'F','[MF]')
                ## Remove nested [] caused by replacement: [A[BC]] becomes [ABC] ##
                while rje.matchExp('\[\S*(\[(\S+)\])',string.replace(var,']','] ')):
                    nesting = rje.matchExp('(\[\S*)(\[(\S+)\])',string.replace(var,']','] '))
                    var = string.replace(var,nesting[0]+nesting[1],nesting[0]+nesting[2])
                newlist.append(var)
            return newlist
        except:
            self.log.errorLog('Error in _msVar:')
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def nofm(self,element): ### Returns list of variants for "n of m" format
        '''Returns list of variants for "n of m" format.'''
        ### Setup ###
        if not rje.matchExp('^(<(\D+):(\d+):(\d+)>)',element):    # "n of m" format
            return []
        match = rje.matchExp('^<(\D+):(\d+):(\d+)>',element)
        aas = match[0]
        if len(aas) > 1:
            aas = '[%s]' % aas
        n = string.atoi(match[1])
        m = string.atoi(match[2])
        if n > m:   # Assume reversed format
            (m,n) = (n,m)
        elif n == m:
            return [aas * m]
        ### Make variants ###
        varlist = []
        binlist = [0] * m
        while sum(binlist) < len(binlist):  # Count in binary!
            if sum(binlist) == n:   # n of m variant
                newvar = ''
                for i in binlist:
                    if i:
                        newvar += aas
                    else:
                        newvar += 'X'
                varlist.append(newvar)
            binlist = rje.binaryCount(binlist)
            #X#self.deBug('%s: %s' % (binlist,varlist))
        return varlist
#########################################################################################################################
    def makeVariants(self,msmode=False,ambvar=True):  ### Makes self.list['Variants'] of length variants (Still PRESTO format)
        '''
        Makes self.list['Variants'] of variants and basic self.dict['Search'] with no mismatches.
        - self.list['Variants'] = with non-regexp variants of different ambiguities etc.
        - self.dict['Search'][0] = regular expressions (length variants) for searching.
        >> msmode:boolean [False] = whether to interpret motifs as MSMS peptides.
        >> ambvar:boolean [True] = whether to make full ambiguity variants for PRESTO search or just dict['Search'][0]
            - this should be set to False for CompariMotif, especially when ambcut is high
        '''
        try:
            ### Length Variations ###
            bases = ['']    # This is a list of basic strings to be added to with each variation
            #X#print self.list['PRESTO']
            for element in self.list['PRESTO']:
                nofm = self.nofm(element) # Returns list of variants for "n of m" format
                nums = rje.matchExp('(\{(\d+),(\d+)\})',element)
                variants = []   # List of variants, comprising of each base and its options!
                for base in bases:
                    if nums:
                        el = element[:-len(nums[0])]
                        for i in range(string.atoi(nums[1]),string.atoi(nums[2])+1):
                            #X#print base, element, i
                            variants.append(base + (el * i))
                    elif nofm:
                        for nmv in nofm:
                            variants.append(base + nmv)
                    else:
                        variants.append(base + element)
                bases = variants[0:]
                #X#self.deBug(bases)
            ## MS Conversion if necessary ##
            variants = bases[0:]
            if msmode:
                variants = self._msVar(variants)
            ## Update Attributes ##
            self.dict['Search'] = {0:variants[0:]}
            
            ### String Variants ###
            changes = ambvar      # Looping boolean
            while changes:
                changes = False
                bases = variants[0:]
                variants = []
                for base in bases:
                    ## Remove termini ##
                    if base[:1] == '^':
                        base = base[1:]
                    if base[-1:] == '$':
                        base = base[:-1]
                    ## Identify ambiguities ##    
                    if base.find('(') >= 0:
                        changes = True
                        choice = base[base.find('('):base.find(')')+1]
                        cvar = string.split(choice[1:-1],'|')
                        for var in cvar:    # Add a variant with each variant!
                            variants.append(base.replace(choice,var,1))
                    elif base.find('[') >= 0:
                        changes = True
                        choice = base[base.find('['):base.find(']')+1]
                        cvar = choice[1:-1]
                        if self.stat['AmbCut'] > 0 and len(cvar) > self.stat['AmbCut']:
                            variants.append(base.replace(choice,'X',1))
                        else:
                            for var in cvar:    # Add a variant with each variant!
                                variants.append(base.replace(choice,var,1))
                    else:
                        variants.append(base)   # Nothing new!
            self.list['Variants'] = variants[0:]
            #x#self.deBug(self.list['Variants'])
        except:
            self.log.errorLog('Error in Motif.makeVariants()',quitchoice=True)
            raise   
#########################################################################################################################       
    def misMatches(self,mismatch={},msmode=False,trimx=False,basevar=False): ### Populates attributes with variants and regular expressions 
        '''
        Populates attributes with variants and regular expressions.
        - self.dict['Search'][X] = regular expressions (length variants) for searching with X mismatches.
        >> mismatch:dictionary of {mm 'X':Y aa}
        >> msmode:boolean = whether to interpret motifs as MSMS peptides.
        >> trimx:boolean = whether to trim Xs from ends of motif (for compare)
        >> basevar:boolean = whether to use self.list['Variants'] as bases rather than self.dict['Search'][0]
        '''
        try:
            ### Setup ### 
            mm = 1
            if basevar:
                bases = self.list['Variants'][0:]
            else:
                bases = self.dict['Search'][0][0:]

            ### Mismatches ###
            while '%d' % mm in mismatch.keys() and self.stat['Length'] >= mismatch['%d' % mm]:  # Potentially have that many mismatches
                variants = []
                for base in bases:
                    temp = base[0:]
                    while rje.matchExp('([A-WY])',temp):
                        aa = rje.matchExp('([A-WY])',temp)[0]
                        i = temp.find(aa)
                        temp = string.replace(temp,aa,'X',1)
                        newvar = string.replace(base[:i] + 'X' + base[i+1:],']','] ')
                        while rje.matchExp('(\[\S*X\S*\])',newvar):
                            newvar = string.replace(newvar,rje.matchExp('(\[\S*X\S*\])',newvar)[0],'X')
                        newvar = string.replace(newvar,' ','')
                        while trimx and newvar[:1] == 'X':
                            newvar = newvar[1:]
                        while trimx and newvar[-1:] == 'X':
                            newvar = newvar[:-1]
                        if newvar not in variants:
                            variants.append(newvar)
                bases = variants[0:]
                if msmode:
                    variants = self._msVar(variants)
                self.dict['Search'][mm] = variants[0:]
                mm += 1
            #x#self.deBug(self.dict['Search'][0])
        except:
            self.log.errorLog('Error in Motif.misMatches()',quitchoice=True)    
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def slimCode(self): ### Makes a SLiMFinder slimcode and stores in self.info['Slim']. Returns True/False.
        '''Makes a SLiMFinder slimcode and stores in self.info['Slim']. Returns code or empty string.'''
        try:
            ### Setup ###
            if 'Slim' in self.info and self.info['Slim']: return self.info['Slim']
            aa = basic_aas + ['^','$']
            checklist = []
            slim = []
            wild = False    # Whether expecting wildcard
            ### Process each element ###
            for el in self.list['PRESTO']:
                if el[0] == '(' or (el.find('{') >= 0 and el[0] != 'X'):
                    self.log.errorLog('Cannot convert SLiM "%s"' % self.info['Sequence'],printerror=False)
                    return False
                elif el[0] in aa:
                    if wild:
                        slim += ['0',el]
                    else:
                        slim.append(el)
                        wild = True
                elif rje.matchExp('^\[(\S+)\]',el):
                    amb = rje.strSort(el[1:-1])
                    if wild:
                        slim += ['0',amb]
                    else:
                        slim.append(amb)
                        wild = True
                elif el[0] == 'X':
                    (m,n) = (1,1)
                    if rje.matchExp('X\{(\d+),(\d+)\}',el):
                        (m,n) = (int(rje.matchExp('X\{(\d+),(\d+)\}',el)[0]), int(rje.matchExp('X\{(\d+),(\d+)\}',el)[1]))
                    if wild and m == n:
                        slim.append('%s' % m)
                    elif wild:
                        w = ''
                        for i in range(m,n+1):
                            w += '%s' % i
                        slim.append(w)
                    else:
                        old = rje.strList(slim[-1])
                        new = []
                        if m == n:
                            for i in range(len(old)):
                                new.append(int(old[i])+m)
                        else:
                            for i in range(len(old)):
                                for x in range(m,n+1):
                                    new.append(int(old[i])+x)
                        new = rje.sortUnique(new,num=True)
                        w = ''
                        for i in new:
                            w += '%s' % i
                        slim[-1] = w
                    wild = False
                else:
                    self.log.errorLog('Cannot convert SLiM "%s"' % self.info['Sequence'],printerror=False)
                    return False
            ### Finish ###
            self.info['Slim'] = string.join(slim,'-')
            return self.info['Slim']
        except:
            self.log.errorLog('Problem with Motif.slimCode(%s)' % self.info['Name'])
            return ''
#########################################################################################################################
    ### <3> ### Motif Statistics                                                                                        #
#########################################################################################################################
    def _calculateLength(self):     ### Calculates Length Statistics
        '''Calculates Length Statistics.'''
        try:
            ### Setup ###
            self.stat['Length'] = 0
            self.stat['MinLength'] = 0
            self.stat['FullLength'] = 0
            self.stat['FixLength'] = 0

            ### Calculate ###
            for element in self.list['PRESTO']:
                ## Ignore ends ##
                if element in ['^','$']:    # Not part of length
                    continue
                ## Work out number of positions ##
                len_add = 0     # Max no. Non-X positions
                min_add = -1    # Min no. Non-X positions
                full_add = 0    # Max no. of all positions
                fix_add = 0     # Max no. of fixed positions
                if element[0] == '(':   # Choice
                    aalist = element[0:]
                    if element.find('{') > 0:
                        aalist = element[:element.find('{')]    # Exclude {m,n} from aalist
                    aalist = string.split(aalist[1:-1],sep='|')
                    for aas in aalist:
                        aax = len(aas) - string.count(aas,'X')
                        if aax > len_add:
                            len_add = aax
                        if min_add < 0 or aax < min_add:
                            min_add = aax
                        if len(aas) > full_add:
                            full_add = len(aas)
                            fix_add = len(aas)  # These are all fixed
                elif rje.matchExp('^(<(\D+):(\d+):(\d+)>)',element):    # "n of m" format
                    match = rje.matchExp('^<(\D+):(\d+):(\d+)>',element)
                    len_add = string.atoi(match[1])
                    min_add = len_add
                    full_add = string.atoi(match[2])
                    if len(match[0]) == 1 and match[0] != 'X':  # Fixed
                        fix_add = string.atoi(match[1])
                else:   # Single or ambiguity
                    len_add = 1 - string.count(element,'X')
                    min_add = len_add
                    full_add = 1
                    if len(element) == 1 and element != 'X':
                        fix_add = 1
                ## Multiple occurrences ##
                nums = rje.matchExp('(\{(\d+),(\d+)\})',element)
                if nums:
                    len_add *= string.atoi(nums[2])
                    min_add *= string.atoi(nums[1])
                    full_add *= string.atoi(nums[2])
                    fix_add *= string.atoi(nums[2])
                self.stat['Length'] += len_add
                self.stat['MinLength'] += min_add
                self.stat['FullLength'] += full_add
                self.stat['FixLength'] += fix_add 
                debugtxt = '%s (%s): Length=%d; MinLength=%d; FullLength=%d; FixLength=%d' % (self.info['Name'],self.info['Sequence'],self.stat['Length'],self.stat['MinLength'],self.stat['FullLength'],self.stat['FixLength'])
            #X#self.deBug(debugtxt)

        except:
            self.log.errorLog('Error in Motif._calculateLength(%s)' % string.join(self.list['PRESTO'],'-'),quitchoice=True)
#########################################################################################################################
    def patternStats(self):      ### Performs calculations based on basic pattern (info['Sequence'])
        '''Performs calculations based on basic pattern (info['Sequence']), adding to self.stat/info/opt.'''
        ### Setup ###
        charge = []
        ailmv = True
        aromatic = 0
        phos = []
        ### Calculations ###
        motseq = string.split(self.info['Sequence'].upper(),']')
        for part in motseq:
            ## Aromatic & Phos ##
            if part.find('F') >= 0 or part.find('W') >= 0 or part.find('Y') >= 0:
                aromatic += 1
            for p in 'STY':
                if part.find(p) >= 0 and p not in phos:
                    phos.append(p)
            ## Charge and AILMV ##
            if part[:1] == '[':  # Ambiguity
                if part in ['[KR','[RK']:   # +ve
                    charge.append(1)
                    ailmv = False
                elif part in ['[DE','[ED']:
                    charge.append(-1)
                    ailmv = False
                elif ailmv:
                    for aa in 'CDEFGHKNPQRSTWY':
                        if part.find(aa) > 0:
                            ailmv = False
                            break
            else:   # Fixed positions
                for a in part: # Motif
                    if a in ['K','R']:
                        charge.append(1)
                        ailmv = False
                    elif a in ['D','E']:
                        charge.append(-1)
                        ailmv = False
                    elif a in ['C','F','G','H','N','P','Q','S','T','W','Y']:
                        charge.append(0)
                        ailmv = False
                    else:
                        charge.append(0)
        ### Update Data ###
        if not phos:
            phos = ['X']
        phos.sort()
        self.setStat({'AbsChg':(charge.count(1) + charge.count(-1)), 'NetChg':sum(charge),
                      'BalChg':(sum(charge[:int(len(charge)/2)]) - sum(charge[-int(len(charge)/2):])),
                      'Aromatic':aromatic})
        self.setOpt({'AILMV':ailmv})
        self.setInfo({'Phos':string.join(phos,'')})
#########################################################################################################################
    ### <4> ### Motif Searching                                                                                         #
#########################################################################################################################
    def searchSequence(self,sequence,logtext=''):  ### Searches the given sequence for occurrences of self and returns hit list
        '''
        Searches the given sequence for occurrences of self and returns a list of hit dictionaries: Pos,Variant,Match
        >> sequence:str = sequence to be searched
        >> logtext:str [''] = text to precede progress printing. If '', no progress printing!
        << hitlist:list of dictionaries with hit information: Pos,Variant,Match,ID,MisMatch
        '''
        try:
            ### Setup ###
            hitlist = []    # List of hit dictionaries to return
            rawhits = []    # List of ['pos:match'] to remove redundancy when mismatches used
            if logtext:
                varx = 0
                for mm in rje.sortKeys(self.dict['Search']):
                    varx += len(self.dict['Search'][mm])
                sx = 0.0
                self.log.printLog('\r#PRESTO','%s [0 mismatch]: %.1f%%' % (logtext,(sx/varx)),log=False,newline=False)
        
            ### Search with self.dict['Search'] ###
            for mm in rje.sortKeys(self.dict['Search']):
                ## Search ##                
                for motvar in self.dict['Search'][mm]:
                    # Each motvar is a fixed length variant. The matchExp should therefore uniquely match the first
                    # occurrence of motvar in the sequence given each time, allowing one to "walk along" the sequence
                    searchvar = string.replace(motvar,'X','[A-Z]')
                    if searchvar[0] == '^':
                        searchvar = '^(' + searchvar[1:]
                    else:
                        searchvar = '(' + searchvar[0:]
                    if searchvar[-1] == '$':
                        searchvar = searchvar[:-1] + ')$' 
                    else:
                        searchvar += ')'
                    #X#print motvar,searchvar
                    r = 0   # Starting point of search
                    while rje.matchExp(searchvar,sequence[r:]):           # Found a match
                        match = rje.matchExp(searchvar,sequence[r:])[0]
                        if searchvar[-1] == '$':
                            r = len(sequence) - len(match) + 1
                        else:
                            r = r + sequence[r:].find(match) + 1        # r is now start of next search
                        if '%d:%s' % (r, match) not in rawhits:     # New hit
                            rawhits.append('%d:%s' % (r, match))
                            hitlist.append(self._hitStats(match,searchvar))
                            hitlist[-1]['Pos'] = r      # Pos is aa number 1->L
                            hitlist[-1]['MisMatch'] = mm
                        if searchvar[0] == '^' or searchvar[-1] == '$':     # Only match once!
                            break
                    ## Progress ##
                    if logtext:
                        sx += 100.0
                        self.log.printLog('\r#PRESTO','%s [%d mismatch]: %.1f%%' % (logtext,mm,(sx/varx)),log=False,newline=False)

            ### Finish ###
            #X# print hitlist #!# This is fine: problem must be in PRESTO conversion #!#
            return hitlist
        except:
            self.log.errorLog('Problem during Motif.searchSequence(%s)' % self.info['Name'])
            return hitlist
#########################################################################################################################
    def _hitStats(self,match,searchvar):   ### Calculates best Variant & ID for Match
        '''
        Calculates best Variant, ID for Match.
        >> match:str = Matched part of sequence
        >> searchvar:str = regular expression variant used in original match
        << histats:dictionary of stats
        '''
        try:
            ### Setup ###
            hitstats = {'Match':match,'Variant':'#ERR','SearchVar':searchvar}
            hitstats['ID'] = 0.0
            while self.stat['AmbCut'] > 0 and rje.matchExp('(\[[A-Z]{%d}[A-Z]*\])' % self.stat['AmbCut'],searchvar):    # AmbCut!
                reptext = rje.matchExp('(\[[A-Z]{%d}[A-Z]*\])' % self.stat['AmbCut'],searchvar)[0]
                searchvar = string.replace(searchvar,reptext,'[A-Z]')
            ### Try each variant and calculate identity to match ###
            for variant in self.list['Variants']:
                variant = string.strip(variant,'$^')
                ## Check possible match ##
                if len(variant) != len(match):
                    continue
                if not rje.matchExp('(%s)' % searchvar,variant):    # Variant is not a possible (mis)match for searchvar
                    continue
                ## Calculate ID and position mismatches (*) ##
                id = 0.0
                r = 0
                while r < len(variant):
                    if variant[r] != 'X' and variant[r] == match[r]:
                        id += 1
                    elif variant[r] != 'X':     # Mismatch! #
                        variant = rje.strSub(variant,r,r,'*')
                    r += 1
                ## Assess ##
                vlen = float(len(variant) - string.count(variant,'X'))
                if (id/vlen) > hitstats['ID']:
                    hitstats['ID'] = (id/vlen)
                    hitstats['Variant'] = variant
                if hitstats['ID'] == 1.0:   # Cannot improve without adding SNT back
                    break
            ### Finish ###
            if hitstats['Variant'] == '#ERR':
                self.log.errorLog('Could not find variant for match "%s" using "%s"' % (match,searchvar),printerror=False)
            return hitstats
        except:
            self.log.errorLog('Problem during Motif.hitStats(%s:%s)' % (self.info['Name'],match))
            return hitstats
#########################################################################################################################
    ### <5> ### Expectation methods                                                                                     #
#########################################################################################################################
    def expectation(self,aafreq={},aanum=0,seqnum=0):   ### Returns a dictionary of {mismatch:expectation}
        '''
        Returns a dictionary of {mismatch:expectation}.
        >> aafreq:dictionary of AA frequencies {aa:freq}
        >> aanum:int = sum total of positions in dataset
        >> seqnum:int = number of different sequence fragments searched
        << expdict:dictionary of {mismatch:expectation}
        '''
        try:
            ### Setup ###
            expdict = {}    # Expectation dictionary
            prevexp = 0.0   # Sum of previous expectation for less mismatches
            ### Expect Scores ###
            for mm in rje.sortKeys(self.dict['Search']):
                mmexp = 0.0 # Expectation for this number of mismatches
                for pattern in self.dict['Search'][mm]:
                    mmexp += expect(pattern,aafreq,aanum,seqnum)
                expdict[mm] = mmexp - prevexp
                prevexp += mmexp
            ### Finish ###
            return expdict
        except:
            self.log.errorLog('Problem during Motif.expectation (%s; %d seq; %d aa)' % (self.info['Sequence'],seqnum,aanum),quitchoice=True)
            return expdict
#########################################################################################################################
### END OF SECTION II: Motif Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: EXPECTATION/PROBABILITY METHODS                                                                         #
#########################################################################################################################
def expect(pattern,aafreq,aanum,seqnum,binomial=False,adjustlen=True):    ### Returns the expected number of occurrences for that pattern
    '''
    Returns the expected number of occurrences for given pattern. Xs and .s both count as wildcards.
    >> aafreq:dictionary of AA frequencies {aa:freq}
    >> aanum:int = sum total of positions in dataset
    >> seqnum:int = number of different sequence fragments searched
    >> binomial:bool [False] = Whether to return n & p data for binomial rather than expectation for poisson
    >> adjustlen:bool [True] = Whether to adjust no. of sites by length of motif
    << expected number of occurrences *or* (prob_per_site,num_sites) if binomial=True
    '''
    ### Setup ###
    terminal_constraint = False
    expvar = string.replace(pattern,')',') ')   # String to analyses
    expvar = string.replace(expvar,'[A-Z]','X')
    expvar = string.replace(expvar,'.','X')
    expvar = string.replace(expvar,']','] ')
    if expvar[0] == '^' or expvar[-1] == '$':
        terminal_constraint = True
    patlen = 0          # Length of pattern for calculating possible number of positions
    aafreq['X'] = 1.0   # Wildcards do not affect expectation (prob=1)
    prob_per_site = 1.0 # Probability of the pattern occurring at any given site.

    ### Calculate prob_per_site ###    
    while expvar:       # Still some pattern to look at
        ## Deal with spacers inserted for ease of pattern matching ##
        if expvar[0] in [' ','^','$']:
            expvar = expvar[1:]
            continue
        ## Wildcard ##
        if expvar[:1] == 'X':   
            expvar = expvar[1:]
            patlen += 1
            continue
        ## Update probability per site ##
        if expvar[0] == '[' and expvar.find(']') > 0:   # Choices
            csum = 0.0      # Summed frequency over choices
            for c in expvar[1:expvar.find(']')]:    # Region between []
                csum += rje.getFromDict(aafreq,c,returnkey=False,case=True,default=0.0)    # Prob = sum of all choices
            if csum < 1.0:
                prob_per_site *= csum
            expvar = expvar[expvar.find(']')+1:]
            patlen += 1
        elif expvar[0] == '(' and expvar.find(')') > 0: # Complex choice!
            csum = 0.0      # Sum prob of whole choice
            cvar = string.split(expvar[1:expvar.find(')')],'|')     # List of different options
            for cv in cvar:
                cvexp = 1.0 # Probability of just one portion of choice
                while rje.matchExp('(\[(\S+)\])',cv):
                    msum = 0.0  # sum for choice within portion!
                    cvm = rje.matchExp('(\[(\S+)\])',cv)
                    cv = string.replace(cv,cvm[0],'',1)
                    for m in cvm[1]:
                        msum += rje.getFromDict(aafreq,m,returnkey=False,case=True,default=0.0)   # Prob = sum of all choices
                    cvexp *= msum
                cv = string.replace(cv,' ','')
                for c in cv:
                    cvexp *= aafreq[c]
                csum += cvexp
            if csum < 1.0:
                prob_per_site *= csum
            expvar = expvar[expvar.find(')')+1:]
            patlen += float(len(string.join(cvar,'')))/len(cvar)    # Add mean length of options
        else:   # Simple
            prob_per_site *= rje.getFromDict(aafreq,expvar[0],returnkey=False,case=True,default=0.0) 
            expvar = expvar[1:]

    ### Convert to Expectation ###
    if terminal_constraint: num_sites = seqnum
    elif adjustlen: num_sites = aanum - (seqnum * (patlen - 1))
    else: num_sites = aanum 
    if binomial: return (prob_per_site,num_sites)
    return prob_per_site * num_sites
    #!# Could use these to calculate binomial probability of 1+ occurrences, i.e. seqprob? #!#
#########################################################################################################################
def occProb(observed,expected):     ### Returns the poisson probability of observed+ occurrences, given expected
    '''Returns the poisson probability of observed+ occurrences, given expected.'''
    prob = 0
    for x in range(0,observed):
        try:        #!# Fudge for OverflowError: long int too large to convert to float
            prob += (math.exp(-expected) * pow(expected,x) / rje.factorial(x))
        except:
            break
    return 1 - prob
#########################################################################################################################
def expectString(_expect):  ### Returns formatted string for _expect value
    '''Returns formatted string for _expect value.'''
    try:
        if _expect > 10:
            return '%.1f' % _expect
        elif _expect > 0.1:
            return '%.2f' % _expect
        elif _expect > 0.001:
            return '%.3f' % _expect
        else:
            return '%.2e' % _expect
    except:
        print expect
        raise
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: IC METHODS                                                                                              #
#########################################################################################################################
def maxInfo(aafreq):    ### Calculates the maximum information content score given aafreq
    '''
    Calculates the maximum information content score given aafreq. Note that this is slightly misleading as this is the
    largest *negative* value, which is then subtracted from the IC measure to give the actual IC value. For fixed
    positions, the IC becomes 0.0 - max_info = - max_info
    '''
    max_info = 0.0
    for a in aafreq:
        max_info += aafreq[a] * math.log(aafreq[a],2)
    return max_info
#########################################################################################################################
def elementIC(element='',aafreq={},wild_pen=0.0,max_info=0.0,callobj=None):      ### Calculates the IC for a given pattern element
    '''
    Calculates the IC for a given pattern element. See Motif.__doc__ for description of elements:
    >> pattern:str = motif pattern
    >> aafreq:dict = aa frequencies
    >> wild_pen:float = wildcard penalty
    '''
    try:
        ### Variable numbers {m,n} - always multiply by m (lowest IC) ###
        multiplier = 1
        if rje.matchExp('\{(\d+),\d+\}',element):
            minmult = string.atoi(rje.matchExp('\{(\d+),\d+\}',element)[0])
            if minmult > 1:
                multiplier = minmult
        ### Wildcard ###
        if element in ['.','X','x']:
            return wild_pen * multiplier
        ### End of Sequence ###
        if element in ['^','$']: return 1.0
        ### Setup AA Freqs & Max Info ###
        if not aafreq:
            aafreq = basic_aafreq
            max_info = math.log(0.05,2)
        if not max_info:
            max_info = maxInfo(aafreq)
        ### Special combos ###
        if element.find('|') > 0:
            ic = 0.0
            combo = string.split(element[1:-1],'|')
            for el in combo:
                for aa in el:
                    ic += elementIC(aa,aafreq,wild_pen,max_info)
            if combo:
                ic /= len(combo)    # Returns mean IC of combos
        ### Special n of m format ###
        elif rje.matchExp('^<(\D+):(\d+):(\d+)>',element):    # "n of m" format
            nofm = rje.matchExp('^<(\D+):(\d+):(\d+)>',element)
            m = string.atoi(nofm[1])
            n = string.atoi(nofm[2])
            ic = elementIC(nofm[0],aafreq,wild_pen,max_info,callobj) * m + elementIC('X',aafreq,wild_pen,max_info,callobj) * (n-m)
            return ic * multiplier
        ### Information Content ###
        else:
            re_aa = re.compile('([A-Z])')
            aas = re_aa.findall(element)
            ic = 0.0   
            denom = 0.0     # Denominator for IC calculation
            for a in aas:
                if aafreq.has_key(a):
                    denom += aafreq[a]
            for a in aas:
                if aafreq.has_key(a):
                    ic += (aafreq[a]/denom) * math.log((aafreq[a]/denom),2)
        return multiplier * (max_info - ic) / max_info
    except:
        if callobj:
            callobj.log.errorLog('Major error during elementIC(%s)' % element,quitchoice=True)
            return 0.0
        raise
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MOTIF DEFINITION METHODS                                                                                 #
#########################################################################################################################
def defineMotif(callobj=None,occlist=[],profile=False,minfreq=0.2,minocc=2,ambcut=19): ### Takes occurrences and makes motif(s) from them
    '''
    Takes occurrences and makes motif(s) from them.
    >> callobj:Object to handle errors etc.
    >> occlist:list of instances. Can be variants (with wildcards) or not
    >> profile:boolean  = whether to return profile-esque patterns with numbers [False]
    >> minfreq:float = min freq of any aa for position to be non-wildcard [0.2]
    >> minocc:int = min number of any aa for position to be non-wildcard (in addition to minfreq) [2]
    >> ambcut:int = number of ambiguities allowed before position marked as wildcard [19]
    << redefined:str = redefined motif (or csv motifs if lengths cannot be compressed using wildcards)
    '''
    try:
        ### Setup ###
        if not occlist:
            return ['']
        elif len(occlist) == 1:
            return occlist
        for i in range(len(occlist)):
            occlist[i] = string.replace(occlist[i],'.','X')
            occlist[i] = string.replace(occlist[i],'x','X')       
        
        ### String Lists ###
        strlist = []
        for occ in occlist:
            newlist = []
            for a in occ.upper():
                if not newlist:
                    newlist = [a]
                elif a == 'X' and newlist[-1][:1] == 'X':
                    newlist[-1] = newlist[-1] + 'X'
                else:
                    newlist.append(a)
            strlist.append(newlist[0:])
        #X#print strlist

        ### Add spacers for {0,1} variable lengths ###        
        #for occ in strlist:
        #    for compocc in strlist:
        #        if len(compocc) > occ:  # Look to add ['']
        #            newlist = []
        #            (o,c) = (0,0)
        #            while o < len(occ) and c < len(compocc):
        #                if occ[o][:1] != 'X' and newocc[o][:1] != 'X'

        ### Compile by length, incorporating variable wildcards ###
        ## Construct ##
        lendict = {}
        for occ in strlist:
            olen = len(occ)
            if not lendict.has_key(olen):
                lendict[olen] = [[] for o in range(olen)]
                #for o in range(olen):
                #    lendict[olen].append([])
                #X#print lendict
            for p in range(len(occ)):
                lendict[olen][p].append(occ[p])
                #X#print occ, p, occ[p], lendict[olen][p]
        #X#print lendict
        ## Compress ##
        cdict = {}  # Compressed sequences
        for olen in rje.sortKeys(lendict):
            cdict[olen] = ''
            for alist in lendict[olen]:
                if alist[0][:1] == 'X' or alist[0] == '':    # Wildcard
                    minx = len(alist[0])
                    maxx = len(alist[0])
                    for aa in alist:
                        if len(aa) < minx:
                            minx = len(aa)
                        elif len(aa) > maxx:
                            maxx = len(aa)
                    if minx == maxx:
                        cdict[olen] += '.' * minx
                    else:
                        cdict[olen] += '.{%d,%d}' % (minx,maxx)
                else:
                    adict = {}
                    maxaa = 1
                    for aa in alist:
                        if adict.has_key(aa):
                            adict[aa] += 1
                            if adict[aa] > maxaa:
                                maxaa = adict[aa]
                        else:
                            adict[aa] = 1
                    if len(adict) == 1:    # Fixed
                        cdict[olen] += rje.sortKeys(adict)[0]
                    elif len(adict) > ambcut or float(maxaa) / len(alist) < minfreq or maxaa < minocc:
                        cdict[olen] += '.'
                    else:   # Ambiguity?! #
                        pos = '['
                        for aa in rje.sortKeys(adict):
                            pos += aa
                            if profile:
                                pos += '%d' % adict[aa]
                        pos += ']'
                        cdict[olen] += pos
        #X#print cdict
                        
        ### Return ###
        rlist = []
        for olen in rje.sortKeys(cdict):
            rlist.append(cdict[olen])
        return rlist

    except:
        if callobj:
            callobj.log.errorLog('Major problem with rje_motif.defineMotif()',quitchoice=True)
            return ''
        raise
#########################################################################################################################
def reformatMiniMotif(callobj=None,pattern=''): ### Reformats minimotif into standard motif format
    '''
    Reformats minimotif into standard motif format.
    >> callobj:Object to handle errors etc.
    >> pattern:str = MiniMotif pattern
    << redefined:str = redefined motif 
    '''    
    try:
        ### Setup ###
        old = pattern
        if pattern[0] == '<':
            pattern = string.replace(pattern,'<','^')
        else:
            pattern = '-%s' % pattern
        pattern = string.replace(pattern,'>',' $')
        pattern = string.replace(pattern,'?','.')
        pattern = string.replace(pattern,'-',' -')
        while rje.matchExp('(\-(\S+) )',pattern):
            match = rje.matchExp('(\-(\S+) )',pattern)
            if string.count(match[1],'/') > 0:
                amb = string.replace(match[1],'/','')
                pattern = string.replace(pattern,match[0],'[%s]' % amb)
            else:
                pattern = string.replace(pattern,match[0],match[1])
        while rje.matchExp('(\-(\S+))$',pattern):
            match = rje.matchExp('(\-(\S+))$',pattern)
            if string.count(match[1],'/') > 0:
                amb = string.replace(match[1],'/','')
                pattern = string.replace(pattern,match[0],'[%s]' % amb)
            else:
                pattern = string.replace(pattern,match[0],match[1])
        pattern = string.replace(pattern,' ','')
        if callobj:
            callobj.log.printLog('\r#MOT','%s => %s' % (old,pattern))
        return pattern
    except:
        if callobj:
            callobj.log.errorLog('Major problem with rje_motif.reformatMiniMotif(%s)' % pattern,quitchoice=True)
            return ''
        raise
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION VI: 'MAIN' PROGRAM                                                                                          #
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        print 'This module is not for standalone running.'
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################
