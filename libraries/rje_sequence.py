#!/usr/bin/python

# rje_sequence - DNA/Protein sequence object
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_sequence
Description:  DNA/Protein sequence object
Version:      2.6.0
Last Edit:    15/11/17
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the Sequence Object used to store sequence data for all PEAT applications that used DNA or
    protein sequences. It has no standalone functionality.

    This modules contains all the methods for parsing out sequence information, including species and source database,
    based on the format of the input sequences. If using a consistent but custom format for fasta description lines,
    please contact me and I can add it to the list of formats currently recognised.

Uses general modules: copy, os, random, re, sre_constants, string, sys, time
Uses RJE modules: rje, rje_disorder
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_disorder
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Separated Sequence object from rje_seq.py
    # 1.1 - Rudimentary opt['GeneSpAcc'] added
    # 1.2 - Modified RegExp for sequence detail extraction
    # 1.3 - Added list of secondary accession numbers and hasID() method to check ID and all AccNum (and gnspacc combos)
    # 1.4 - Added Peptide Design methods
    # 1.5 - Added storing of case in a dictionary self.dict['Case'] = {'Upper':[(start,stop)],'Lower':[(start,stop)]}
    # 1.6 - Added disorder and case masking
    # 1.7 - Added FudgeFactor and AA codes
    # 1.8 - Added position-specific AA masking
    # 1.9 - Added EST translation functions. Fixed fudging. Added dna() method.
    # 1.10- Fixed sequence name bug
    # 1.11- Added recognition of UniRef
    # 1.12- Added AA masking
    # 1.13- Added Taxonomy list and UniProt dictionary for UniProt sourced sequences (primarily).
    # 1.14- Added maskRegion()
    # 1.15- Added disorder proportion calculations.
    # 1.16- Added additional Genbank and EnsEMBL BioMart sequence header recognition.
    # 1.17- Added nematode sequence conversion.
    # 2.0 - Replaced RJE_Object with RJE_ObjectLite.
    # 2.1 - Added re_unirefprot = re.compile('^([A-Za-z0-9\-]+)\s+([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
    # 2.2 - Added more yeast species.
    # 2.3 - Added alternative self.info keys for sequence (for UniProt splice variants). Added SpliceVar dict.
    # 2.4 - Added recognition of modified IPI format. Added standalone low complexity masking.
    # 2.4.1 - Moved the gnspacc fragment recognition to reduce issues. Should perhaps remove completely?
    # 2.5.0 - Added yeast genome renaming.
    # 2.5.1 - Modified reverse complement code.
    # 2.5.2 - Tried to speed up dna2prot code.
    # 2.5.3 - Fixed genetic code warning error.
    # 2.6.0 - Added mutation dictionary to Ks calculation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add descriptions of recognised formats to documentation.
    # [ ] : Make more efficient sequence data extraction for larger sequence datasets.
    # [ ] : What uses gnspacc fragment format? (Can't imagine anything using it!)
    # [ ] : Add alternative initiator codons to genetic codes.
    '''
#########################################################################################################################
### Regular Expressions for Sequence Details Extraction
re_acconly = re.compile('^(\S+)\s*$')
re_plain = re.compile('^(\S+)\s+(\S.*)$')
re_gnspec = re.compile('^(\S+)_(\S+)')
re_gnspacc = re.compile('^(\S+)_([A-Za-z0-9]+)\/(\S+)')
re_gn_sp__acc = re.compile('^(\S+)_([A-Za-z0-9]+)__(\S+)')
re_uniprot = re.compile('^(sp|tr|sw|uniprot)\|(\S+)\|(\S+)_(\S+)')
re_uniref = re.compile('^(UniRef\S+)_(\S+)\s.*\Tax=(\S.+)\sRepID=(\S+)')
re_uniprot2 = re.compile('^(\S+)_(\S+)\s+\((\S+)\)')
re_uniprot3 = re.compile('^([A-Za-z0-9\-]+)\|([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
re_unirefprot = re.compile('^([A-Za-z0-9\-]+)\s+([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
re_disprot = re.compile('^DisProt\|(DP\d+)\|')
re_genbank = re.compile('^(\S+)\|(\d+)\|')
re_ncbi = re.compile('^gi\|(\d+)\|(\S+)\|(\S+)\|')
re_ncbi2 = re.compile('^gi\|(\d+)\|(\S+)\|(\S+)')
re_unigene = re.compile('^\S+\|UG\|(\S+)')
re_ensemblpep = re.compile('^(\S+)\s+pep:(\S+)')
re_ipi = re.compile('^>*IPI:(IPI\d+\.\d+)\|.+Tax_Id=(\d+)')
re_ipimod= re.compile('^IPI:(IPI\d+\.\d+)\s+Gene_Symbol=(\S+)')
re_hprd = re.compile('^ID_HPRD_(\d+)')
re_hprd2 = re.compile('^(\d+)\|(\d+_\d+)\|(\S+)\|')
re_pipe = re.compile('^(\S+)_(\S+)\|(\S+)')
re_db_pipe = re.compile('^(\S+)\|(\S+)')
re_tigr = re.compile('^(\S+)\s.+\{(\S.+)\}')
re_flybase = re.compile('^(\S+) type=.+name=(\S+);.+species=(\S+);')
re_jgi = re.compile('^jgi\|(\S+)\|(\d+)\|(\S.+)')
re_enst = re.compile('^(ENS\S+)\|ENS(\S+)\|')
ipi_taxa = {'3702':'ARATH', '7955':'BRARE', '9031':'CHICK', '9606':'HUMAN', '10090':'MOUSE', '10116':'RAT'}
#########################################################################################################################
genetic_code = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'S','UCA':'S','UCG':'S','UAU':'Y','UAC':'Y','UAA':'*','UAG':'*','UGU':'C','UGC':'C','UGA':'*','UGG':'W',
                'CUU':'L','CUC':'L','CUA':'L','CUG':'L','CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R','CGC':'R','CGA':'R','CGG':'R',
                'AUU':'I','AUC':'I','AUA':'I','AUG':'M','ACU':'T','ACC':'T','ACA':'T','ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R',
                'GUU':'V','GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G','GGA':'G','GGG':'G'}
aa_code_3 = {'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His','I':'Ile','L':'Leu','M':'Met',
             'P':'Pro','Q':'Gln','R':'Arg','S':'Ser','T':'Thr','V':'Val','X':'Unk','Y':'Tyr','*':'STOP'}
aa_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'MET': 'M', 'THR': 'T', 'PRO': 'P', 'STOP': '*', 'HIS': 'H',
           'PHE': 'F', 'ALA': 'A', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'UNK': 'X', 'VAL': 'V', 'GLU': 'E',
           'TYR': 'Y'}
transl_table = {'1':genetic_code,   # Standard
                '2':rje.combineDict({'AGA':'*','AGG':'*','AUA':'M','UGA':'W'},genetic_code,overwrite=False),# Vert. mito.
                '3':rje.combineDict({'AUA':'M','CUU':'T','CUC':'T','CUA':'T','CUG':'T','UGA':'W','CGA':'!','CGC':'!'},genetic_code,overwrite=False), # Yeast mito.
                '11':rje.combineDict({},genetic_code,overwrite=False)} # Bacterial, Archaeal and Plant Plastid Code (transl_table=11). Multiple initiator codons.
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Sequence Class for Individual sequences                                                                 #
#########################################################################################################################
class Sequence(rje.RJE_ObjectLite):
    '''
    Individual Sequence Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of sequence
    - Type = Type of sequences (RNA/DNA/Protein)
    - Sequence = Actual sequence
    - Description = Description of sequence list (if desired)
    - ID = Sequence identifier (unique 'word')
    - AccNum = Accession Number
    - DBase = Source database (for AccNum)
    - Gene = Gene Symbol
    - Species = species name
    - SpecCode = SwissProt Species Code
    - Format = Name Format
    - NCBI = NCBI GenBank gi number
    
    Stat:numeric

    Opt:boolean
    - RevComp = Whether sequence has been reverse complemented or not [False]

    List:list
    - IsDisordered = List of True/False values for each residue, whether disordered or not
    - Secondary ID = List of secondary IDs and Accession numbers
    - Taxonomy = List of taxonomic levels for source species
    
    Dict:dictionary
    - Case = Stores case of original input sequence as tuples {'Upper':[(start,stop)],'Lower':[(start,stop)]}
    - SpliceVar = Dictionary of {Splice AccNum:Isoform Name} (key also used in self.info to store splice variant sequence
    - UniDAT = UniProt data dictionary (if sequence read from UniProt entry)

    Obj:RJE_Objects
    - Disorder = rje_disorder.Disorder object containing disorder prediction results. (Must run self.disorder() first!)
    '''
    ### Attributes
    def seqLen(self): return len(self.getInfo('Sequence'))
    def aaLen(self): return len(self.getInfo('Sequence')) - string.count(self.getInfo('Sequence'),'-')
    def aaNum(self): return self.aaLen()
    def nonX(self):
        if self.dna(): return self.nonN()
        return self.aaLen() - string.count(self.getInfo('Sequence').upper(),'X')
    def nonN(self): return self.aaLen() - string.count(self.getInfo('Sequence').upper(),'N')
    def dna(self): return self.seqType() in ['RNA','DNA']   # Whether a DNA (or RNA) sequence
    def unit(self): {True:'nt',False:'aa'}[self.dna()]      # Return appropriate unit
    def MWt(self): return MWt(self.getSequence(gaps=False))
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','Description','Sequence','ID','AccNum','DBase','Gene','Species','SpecCode','Format',
                         'PreMask','MaskSeq','NCBI']
        self.statlist = ['Verbose','Interactive']
        self.optlist = ['RevComp','Special']
        self.listlist = ['IsDisordered','Secondary ID','Taxonomy']
        self.dictlist = ['Case','SpliceVar','UniDAT']
        self.objlist = ['Disorder']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'SpecCode':'UNK','NCBI':''})
        self.dict['Case'] = {'Upper':[],'Lower':[]}
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'opt',['Yeast','Nematode','Special'])
                self._cmdRead(cmd,'info','SpecCode','spcode')
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.setInfo({'SpecCode':self.info['SpecCode'].upper()})
        if not self.getInfo('SpecCode'): self.setInfo({'SpecCode':'UNK'})
        return
#########################################################################################################################
    ### <2> ### Sequence Data (from name) extraction and updates                                                        #
#########################################################################################################################
    def addSequence(self,sequence='',case=True,caselist=[],stripnum=False): ### Adds sequence to object, extracting case information
        '''
        Adds sequence to object, extracting case information.
        >> sequence:str [''] = sequence to add to self.info['Sequence']
        >> case:bool [True] = whether to store case information in self.dict['Case']
        >> caselist:list [] = list of case positions to over-ride read in case (self.list['Case'] from rje_seq)
        >> stripnum:boolean [False] = whether to strip numbers from sequences
        '''
        try:### ~ [0] ~ Strip Space and Numbers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sequence = string.join(string.split(sequence),'')
            if stripnum:
                for badness in ['0','1','2','3','4','5','6','7','8','9']: sequence = string.replace(sequence,badness,'')
            ### ~ [1] ~ Update Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setInfo({'Sequence':sequence.upper()})
            ## ~ [1a] ~ Case Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Case'] = {'Upper':[],'Lower':[]}
            if caselist:
                switch = {'Upper':'Lower','Lower':'Upper'}
                seq = {'Upper':sequence.upper(),'Lower':sequence.lower()}
                clist = [0]
                for i in caselist:
                    try: clist.append(string.atoi(i))
                    except: self.errorLog('Cannot use "%s" from caselist!' % i)
                c = 'Upper'
                while clist:
                    r = clist.pop(0)
                    sequence = sequence[:r] + seq[c][r:]
                    c = switch[c]
            if case: self.dict['Case'] = caseDict(sequence)
        except: self.errorLog('Problem with rje_sequence.addSequence()')
#########################################################################################################################
    def fasta(self,isoform=None,case=False,gaps=True,strict=False): ### Returns fasta formatted text
        '''
        Returns fasta formatted text using self.name() and optional isoform information.
        >> isoform:str [None] = isoform to output. All if 'All'. Nothing returned if missing,
        >> case:bool [False] = whether to use self.dict['Case'] information.
        >> gaps:bool [True] = whether to leave gaps in the sequence
        >> strict:bool [False] = whether to raise error if isoform is missing.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not isoform or isoform in [self.getStr('AccNum'),'%s-1' % self.getStr('AccNum')]: ikey = 'Sequence'
            elif not strict and isoform not in self.info: ikey = 'Sequence'
            elif strict and isoform not in self.info: raise ValueError('%s is not a recognised isoform of %s! (Secondary AccNum?)' % (isoform,self.getStr('AccNum')))
            else: ikey = isoform
            ### ~ [1] All sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if isoform.lower() == 'all':
                fastatxt = ''
                if '%s-1' % self.getStr('AccNum') not in self.dict['SpliceVar']:
                    self.dict['SpliceVar']['%s-1' % self.getStr('AccNum')] = '1'
                for isoform in rje.sortKeys(self.dict['SpliceVar']): fastatxt += self.fasta(isoform,case,gaps)
                return fastatxt
            ### ~ [2] One isoform ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
            return '>%s\n%s\n' % (self.name(isoform),self.getSequence(case,gaps,ikey))
        except ValueError: raise
        except: self.errorLog('Problem with rje_sequence.fasta()')
#########################################################################################################################
    def getSequence(self,case=False,gaps=True,ikey='Sequence'):   ### Returns sequence
        '''
        Returns sequence.
        >> case:bool [False] = whether to use self.dict['Case'] information.
        >> gaps:bool [True] = whether to leave gaps in the sequence
        >> ikey:str [Sequence] = self.info key pointing to sequence
        '''
        if ikey not in self.info: 
            if ikey and ikey in self.info.values(): return self.getSequence(case,gaps)
            return ''
        if ikey == '%s-1' % self.getStr('AccNum')  and ikey not in self.info: ikey = 'Sequence'
        sequence = self.info[ikey].upper()
        if case:
            for (start,stop) in self.dict['Case']['Lower']:
                sequence = sequence[:start] + sequence[start:(stop+1)].lower() + sequence[(stop+1):]
        if gaps: return sequence
        else: return string.replace(self.info[ikey],'-','')
#########################################################################################################################
    def reverseComplement(self,rna=False):  ### Converts to reverse complement of sequence (upper case)
        '''Converts to reverse complement of sequence (upper case).'''
        self.info['Sequence'] = reverseComplement(self.info['Sequence'],rna=rna)
        self.opt['RevComp'] = not self.opt['RevComp']
#########################################################################################################################
    def trimPolyA(self):    ### Removes 3' As
        '''Removes 3' As.'''
        while self.info['Sequence'][-1:] == 'A': self.info['Sequence'] = self.info['Sequence'][:-1]
#########################################################################################################################
    def extractDetails(self,gnspacc=False):   ### Extracts ID, Acc# etc. from sequence name
        '''Extracts ID, Acc# etc. from sequence name.'''
        ### ~ [1] ~ General Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.info['FullName'] = self.info['Name']
        try:
            self.setInfo(extractNameDetails(self.info['Name'],self))
            if self.info['Format'] == 'gn_sp__acc':
                self.info['Name'] = string.join(['%s__%s' % (self.info['ID'],self.info['AccNum'])] + string.split(self.info['Name'])[1:])
        except: self.errorLog('Error extracting basic sequence name information.'); return 0
        ### ~ [2] ~ Special details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getAtt('opt','Special',default=False) or self.getAtt('opt','Yeast',default=False): self.specialDetails()
        #print self.info, self.opt
        #raw_input('...')
        ### ~ [3] ~ Gene_Species__Acc formatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            for g in self.gene()[0:]:
                if not rje.matchExp('([A-Za-z0-9_-])',g) and g not in ['.','#']: self.gene(string.replace(self.gene(),g,''))
            if gnspacc and self.info['Format'] != 'gn_sp__acc':
                if self.gene().lower() in ['','none']:
                    self.setInfo({'ID':'%s_%s' % (self.getInfo('AccNum'),self.getInfo('SpecCode'))})
                else: self.setInfo({'ID':'%s_%s' % (self.gene(),self.getInfo('SpecCode'))})
                self.info['Name'] = '%s__%s' % (self.info['ID'],self.info['AccNum'])
                if self.info['Description'].lower() not in ['','none']: self.info['Name'] += ' %s' % self.info['Description']
                self.info['Format'] == 'gn_sp__acc'
        except: self.errorLog('Gene_Species__Acc formatting error (%s).' % self.info['Name'])
        return 1
#########################################################################################################################
    def specialDetails(self):   ### Access special information from Names
        '''Access special information from Names.'''
        try:### ~ [1] ~ Special Kate Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'].find('TP_') >= 0 and self.info['Name'].find('{Treponema pallidum Nichols}') >= 0:
                self.setInfo({'Gene':'p','SpecCode':'TREPA','ID':'%s_%s' % (self.getInfo('Gene'),self.getInfo('SpecCode'))})
            ### ~ [2] ~ Yeast ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getAtt('opt','Yeast',default=False):
                if rje.matchExp('gi\|\d+\|gb\|(\S+)\| Saccharomyces cerevisiae (\S+) (\S.+\S)$',self.info['Name']):
                    #gi|803443638|gb|CP005241.1| Saccharomyces cerevisiae YJM1304 chromosome VII sequence
                    #gi|785784932|gb|CP004567.1| Saccharomyces cerevisiae YJM1615 plasmid 2 micron, complete sequence
                    #gi|768457878|gb|CP006507.1| Saccharomyces cerevisiae YJM996 mitochondrion, complete genome
                    #gi|763972742|gb|CP006294.1| Saccharomyces cerevisiae YJM1342 chromosome III
                    ydata = rje.matchExp('gi\|\d+\|gb\|(\S+)\| Saccharomyces cerevisiae (\S+) (\S.+\S)$',self.info['Name'])
                    gene = string.split(string.split(ydata[2],',')[0])
                    if len(gene) > 1: gene = gene[0][:3] + gene[1]
                    else: gene = gene[0][:3]
                    self.setInfo({'Gene':gene,'SpecCode':ydata[1],'AccNum':ydata[0],'Description':self.info['Name']})
                    self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                    #self.obj['Parent'].deBug('%s -> %s %s' % (self.info['Name'],self.info['ID'],self.info['AccNum']))
                    return
                self.obj['Parent'].deBug(self.info['Name'])
                #!# Add yeast genome recognition here!
                if self.info['Description'].lower() in ['','none']: self.info['Description'] = '[YGOB: %s]' % self.shortName()
                else: self.info['Description'] += ' [YGOB: %s]' % self.shortName()
                if self.shortName()[4:5] == '_': self.info['AccNum'] = self.shortName()[:4].upper() + self.shortName()[5:]
                if self.shortName()[:4] == 'Kpol': self.setInfo({'SpecCode':'VANPO','Gene':'Kpol'})
                elif self.shortName()[:4] == 'TPHA': self.setInfo({'SpecCode':'TETPH','Gene':'Tpha'})
                elif self.shortName()[:4] == 'TBLA': self.setInfo({'SpecCode':'TETBL','Gene':'Tbla'})
                elif self.shortName()[:4] == 'NDAI': self.setInfo({'SpecCode':'NAUDC','Gene':'Ndai'})
                elif self.shortName()[:4] == 'NCAS': self.setInfo({'SpecCode':'NAUCA','Gene':'Ncas'})
                elif self.shortName()[:4] == 'KNAG': self.setInfo({'SpecCode':'KAZNA','Gene':'Knag'})
                elif self.shortName()[:4] == 'KAFR': self.setInfo({'SpecCode':'KAZAF','Gene':'Kafr'})
                elif self.shortName()[:4] == 'CAGL': self.setInfo({'SpecCode':'CANGA','Gene':'Cgla','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Cgla': self.setInfo({'SpecCode':'CANGA','Gene':'Cgla'})
                elif self.shortName()[:4] == 'Sbay': self.setInfo({'SpecCode':'SACBA','Gene':'Sbay'})
                elif self.shortName()[:4] == 'Suva': self.setInfo({'SpecCode':'SACBA','Gene':'Sbay'})
                elif self.shortName()[:4] == 'Skud': self.setInfo({'SpecCode':'SACKU','Gene':'Skud'})
                elif self.shortName()[:4] == 'Smik': self.setInfo({'SpecCode':'SACMI','Gene':'Smik'})
                elif self.shortName()[:4] == 'Scer': self.setInfo({'SpecCode':'YEAST','Gene':'Scer'})
                elif self.shortName()[:4] == 'Anc_': self.setInfo({'SpecCode':'YANC','Gene':'Anc','AccNum':'Anc%s' % self.shortName()[4:]})
                elif self.shortName()[:3] == 'ZYR': self.setInfo({'SpecCode':'ZYGRO','Gene':'Zrou','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Zrou': self.setInfo({'SpecCode':'ZYGRO','Gene':'Zrou'})
                elif self.shortName()[:4] == 'TDEL': self.setInfo({'SpecCode':'LACDE','Gene':'Ldel'})
                elif self.shortName()[:4] == 'Klac': self.setInfo({'SpecCode':'KLULA','Gene':'Klac'})
                elif self.shortName()[:4] == 'KLLA': self.setInfo({'SpecCode':'KLULA','Gene':'Klac','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Agos': self.setInfo({'SpecCode':'ASHGO','Gene':'Egos'})
                elif self.shortName()[:4] == 'Ecym': self.setInfo({'SpecCode':'ERECY','Gene':'Ecym'})
                elif self.shortName()[:4] == 'Sklu': self.setInfo({'SpecCode':'SACKL','Gene':'Sklu'})
                elif self.shortName()[:4] == 'SAKL': self.setInfo({'SpecCode':'SACKL','Gene':'Sklu','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'KLTH': self.setInfo({'SpecCode':'LACTH','Gene':'Kthe','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Kthe': self.setInfo({'SpecCode':'LACTH','Gene':'Kthe'})
                elif self.shortName()[:4] == 'Kwal': self.setInfo({'SpecCode':'KLUWA','Gene':'Kwal'})
                elif self.shortName()[:4] == 'Scas': self.setInfo({'SpecCode':'SACCA','Gene':'Scas'})
                elif rje.matchExp('(Y\S+\d+[WC])',self.shortName()): self.setInfo({'SpecCode':'YEAST','Gene':'Scer','AccNum':self.shortName()})
                elif rje.matchExp('(A\S+\d+[WC])',self.shortName()): self.setInfo({'SpecCode':'ASHGO','Gene':'Egos','AccNum':self.shortName()})
                self.info['AccNum'] = string.replace(self.info['AccNum'],'YGOB_','')
                self.info['AccNum'] = string.replace(self.info['AccNum'],'Anc_','Anc-')
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.obj['Parent'].deBug('%s -> %s %s' % (self.shortName(),self.info['ID'],self.info['AccNum']))
                return
            ## ~ [2a] ~ NEMBASE4 Nematode sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getAtt('opt','Nematode',default=False) and (rje.matchExp('^(\D\D\D\S+)_\d+',self.info['Name']) or rje.matchExp('^(\D\D\S\d+)\s*',self.info['Name'])):
                if rje.matchExp('^(\D\D\D\S+)_\d+',self.info['Name']):
                    acc = rje.matchExp('^(\D\D\D\S+_\d+)',self.info['Name'])[0]
                    spcode = species = sp = acc[:3]
                else: 
                    acc = rje.matchExp('^(\D\D\S\d+)\s*',self.info['Name'])[0]
                    spcode = species = sp = acc[:2]
                if sp == 'AAP': species = 'Angiostrongylus cantonensis'
                if sp == 'ABP': species = 'Ancylostoma braziliense'
                if sp == 'ACP': species = 'Ancylostoma caninum'
                if sp == 'AIP': species = 'Anisakis simplex'
                if sp == 'ALP': species = 'Ascaris lumbricoides'
                if sp == 'ASP': species = 'Ascaris suum'
                if sp == 'AYP': species = 'Ancylostoma ceylanicum'
                if sp in ['BMP','BM']: species = 'Brugia malayi'
                if sp == 'BPP': species = 'Brugia pahangi'
                if sp == 'BUP': species = 'Bursaphelenchus mucronatus'
                if sp == 'BXP': species = 'Bursaphelenchus xylophilus'
                if sp == 'CBP': species = 'Caenorhabditis brenneri'; spcode = 'CAEBE'
                if sp == 'CEP': species = 'Caenorhabditis elegans'
                if sp == 'CGP': species = 'Caenorhabditis briggsae'
                if sp == 'CJP': species = 'Caenorhabditis japonica'
                if sp == 'CRP': species = 'Caenorhabditis remanei'
                if sp == 'CSP': species = 'Caenorhabditis sp.'
                if sp == 'DAP': species = 'Ditylenchus africanus'
                if sp == 'DIP': species = 'Dirofilaria immitis'
                if sp == 'DVP': species = 'Dictyocaulus viviparus'
                if sp == 'GMP': species = 'Globodera mexicana'; spcode = 'GLOMX'
                if sp == 'GPP': species = 'Globodera pallida'
                if sp == 'GRP': species = 'Globodera rostochiensis'
                if sp == 'HBP': species = 'Heterorhabditis bacteriophora'
                if sp == 'HCP': species = 'Haemonchus contortus'
                if sp == 'HGP': species = 'Heterodera glycines'
                if sp == 'HSP': species = 'Heterodera schachtii'
                if sp == 'LLP': species = 'Loa loa'
                if sp == 'LSP': species = 'Litomosoides sigmodontis'
                if sp == 'MAP': species = 'Meloidogyne arenaria'
                if sp == 'MCP': species = 'Meloidogyne chitwoodi'
                if sp == 'MHP': species = 'Meloidogyne hapla'
                if sp == 'MIP': species = 'Meloidogyne incognita'
                if sp == 'MJP': species = 'Meloidogyne javanica'
                if sp == 'MPP': species = 'Meloidogyne paranaensis'
                if sp == 'NAP': species = 'Necator americanus'
                if sp == 'NBP': species = 'Nippostrongylus brasiliensis'
                if sp == 'OCP': species = 'Onchocerca ochengi'
                if sp == 'ODP': species = 'Oesophagostomum dentatum'
                if sp == 'OFP': species = 'Onchocerca flexuosa'
                if sp == 'OOP': species = 'Ostertagia ostertagi'
                if sp == 'OVP': species = 'Onchocerca volvulus'
                if sp == 'PAP': species = 'Parelaphostrongylus tenuis'; spcode = 'PARTN'
                if sp == 'PEP': species = 'Pratylenchus penetrans'
                if sp in ['PPP','PP']: species = 'Pristionchus pacificus'
                if sp == 'PSP': species = 'Panagrolaimus superbus'
                if sp == 'PTP': species = 'Parastrongyloides trichosuri'
                if sp == 'PVP': species = 'Pratylenchus vulnus'
                if sp == 'RSP': species = 'Radopholus similis'
                if sp == 'SCP': species = 'Steinernema carpocapsae'; spcode = 'STECR'
                if sp == 'SFP': species = 'Steinernema feltiae'
                if sp == 'SRP': species = 'Strongyloides ratti'; spcode = 'STORA'
                if sp == 'SSP': species = 'Strongyloides stercoralis'; spcode = 'STOST'
                if sp == 'TCP': species = 'Toxocara canis'
                if sp == 'TDP': species = 'Teladorsagia circumcincta'
                if sp == 'TIP': species = 'Trichostrongylus vitrinus'; spcode = 'TRIVN'
                if sp == 'TLP': species = 'Toxascaris leonina'; spcode = 'TOXLN'
                if sp == 'TMP': species = 'Trichuris muris'; spcode = 'TRIMR'
                if sp == 'TSP': species = 'Trichinella spiralis'
                if sp == 'TVP': species = 'Trichuris vulpis'; spcode = 'TRIVL'
                if sp == 'WBP': species = 'Wuchereria bancrofti'
                if sp == 'XIP': species = 'Xiphinema index'
                if sp == 'ZPP': species = 'Zeldia punctata'
                if sp == species: species = 'Nematode sp.'
                if spcode == sp: spcode = species[:3].upper() + string.split(species)[-1][:2].upper()
                self.info['Gene'] = 'nem'   #sp.lower()
                self.info['Species'] = species
                self.info['SpecCode'] = spcode
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = acc
                self.info['DBase'] = 'NEMBASE4'
                self.info['Description'] = string.replace(self.info['Name'],'->',' ')
                return
            ### ~ [3] ~ Special E hux EST translation consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if (self.shortName()[:3] == 'EHC' or self.shortName()[:4] == 'EHUX') and self.shortName() == self.info['Name']:
                self.info['Gene'] == 'est'
                self.info['SpecCode'] = 'EMIHU'
                self.info['Species'] = 'Emiliania huxleyi'
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            ### ~ [4] ~ eORF analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            eorf_re = '^(\S+.+)\s(\S+\|.+)gene:(\S+)'
            if rje.matchExp(eorf_re,self.info['Name']):
                eorf = rje.matchExp(eorf_re,self.info['Name'])
                acc = string.split(eorf[0])[0]
                if len(string.split(acc,'_')) > 2: return
                gene = string.split(eorf[1],'|')[0]
                if '_CAEEL' in gene: (gene,spcode) = string.split(gene,'_')
                elif 'WS200' in eorf[1]: spcode = 'CAEEL'
                elif acc[:7] == 'ENSGALT': spcode = 'CHICK'
                elif acc[:7] == 'ENSBTAT': spcode = 'BOVIN'
                elif acc[:4] == 'FBtr': spcode = 'DROME'
                elif acc[:4] == 'ENST': spcode = 'HUMAN'
                elif acc[:7] == 'ENSMUST': spcode = 'MOUSE'
                elif acc[:7] == 'ENSSSCT': spcode = 'PIG'
                elif acc[:7] == 'ENSXETT': spcode = 'XENTR'
                elif acc[:7] == 'ENSDART': spcode = 'DANRE'
                elif 'SGD' in eorf[1]: spcode = 'YEAST'
                if len(gene) < 2 or len(gene) > 6: gene = 'eorf'
                self.info['Gene'] = gene.lower()
                self.info['AccNum'] = acc
                self.info['SpecCode'] = spcode
            elif self.info['SpecCode'] == 'T00': (self.info['SpecCode'],self.info['Species']) = ('HUMAN','Homo sapiens')
            ### ~ [5] ~ Chlam Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'chlam' in self.cmd_list:
                chlam_re = '^lcl\|(\S+)_cdsid_\S+ \[gene=(\S+)\] \[protein=(\S[^\]]+)\].+'
                chlam_re2 = '^lcl\|(\S+)_cdsid_(\S+) \[protein=(\S[^\]]+)\].+'
                chlam_re3 = '^lcl\|(\S+)_cdsid_(\S+) \[gene=\S[^\]]+\] \[protein=(\S[^\]]+)\].+'
                cdata = rje.matchExp(chlam_re,self.info['Name'])
                if not cdata: cdata = rje.matchExp(chlam_re2,self.info['Name'])
                if not cdata: cdata = rje.matchExp(chlam_re3,self.info['Name'])
                if not cdata: print self.info['Name'], cdata; raw_input('???'); return
                self.info['Gene'] = cdata[1]
                if not rje.matchExp(chlam_re,self.info['Name']) and string.split(cdata[2])[-1][:2] == 'CT': self.info['Gene'] = string.split(cdata[2])[-1]; self.info['DBase'] = 'ct'
                elif rje.matchExp(chlam_re,self.info['Name']): self.info['DBase'] = 'gene'
                else: self.info['DBase'] = 'gb'
                self.info['AccNum'] = self.info['Gene'] # cdata[0]
                self.info['SpecCode'] = 'CHLT2'
                self.info['Species'] = 'Chlamydia trachomatis'
                
        except: self.errorLog('Special Details problem for %s' % self.info['Name'])
#########################################################################################################################
    def newGene(self,gene='p',keepsp=False,gnspacc=True):     ### Gives sequence new gene
        '''
        Gives sequence new gene.
        >> gene:str = new gene
        >> keepsp:str = whether to keep an UPPER CASE gene identifier
        '''
        if keepsp:
            if self.info['Gene'].find(gene.upper()) >= 0:
                if self.info['Name'].find('acc:%s' % self.info['Gene']) < 0: return
            elif self.info['DBase'] == 'sprot': return  # SPROT
        self.info['Gene'] = gene
        self.info['ID'] = '%s_%s' % (gene,self.info['SpecCode'])
        if self.info['Format'] == 'gnspacc':
            self.info['Name'] = '%s/%s %s' % (self.info['ID'],self.info['AccNum'],self.info['Description'])
        elif self.info['Format'] in ['gn_sp__acc','uniprot2'] or gnspacc:
            self.info['Name'] = '%s__%s %s' % (self.info['ID'],self.info['AccNum'],self.info['Description'])
        else:
            self.verbose(0,3,'Gene but not Name changed for %s, (%s format).' % (self.shortName(),self.info['Format']),1)
#########################################################################################################################
    def deGap(self):    ### Degaps sequence
        '''Degaps sequence.'''
        self.info['Sequence'] = string.join(string.split(self.info['Sequence'],'-'),'')
#########################################################################################################################
    ### <3> ### Sequence information methods                                                                            #
#########################################################################################################################
    def gene(self,newgene=None):    ### Returns gene or sets new gene if newgene given
        if newgene: self.setInfo({'Gene':newgene})
        return self.getInfo('Gene')
#########################################################################################################################
    def name(self,isoform=None):  ### Returns name with optional splice isoform accnum replacement 
        '''Returns name with optional splice isoform accnum replacement.'''
        if isoform in self.dict['SpliceVar']:
            if self.info['Format'] == 'gnspacc':
                name = '%s/%s %s' % (self.info['ID'],isoform,self.info['Description'])
            elif self.info['Format'] in ['gn_sp__acc','uniprot2'] or re_gn_sp__acc.match(self.name()):
                name = '%s__%s %s' % (self.info['ID'],isoform,self.info['Description'])
            else: name = '%s %s' % (string.replace(self.shortName(),self.info['AccNum'],isoform),self.info['Description'])
            name += ' SpliceVar: isoform %s' % self.dict['SpliceVar'][isoform]
            return name
        return self.getInfo('Name')
#########################################################################################################################
    def shortName(self):    ### Returns short name.
        '''Returns short name = first word of name.'''
        try: return string.split(self.info['Name'])[0]
        except:
            self.log.errorLog('Major problem with shortName(%s)' % self.info['Name'])
            raise
#########################################################################################################################
    def seqType(self):  ### Returns (and possible guesses) Sequence Type - Protein/DNA/RNA
        '''
        Returns (and possible guesses) Sequence Type
        - Protein if non-ATGCUN
        - DNA if ATGCN only
        - RNA if AUGCN only
        '''
        if self.info['Type'] == 'None': # Work it out
            if re.search('[DEFHIKLMPQRSVWY]',self.info['Sequence'].upper()): self.info['Type'] = 'Protein'
            elif re.search('U',self.info['Sequence'].upper()): self.info['Type'] = 'RNA'
            else: self.info['Type'] = 'DNA'
        return self.info['Type']            
#########################################################################################################################
    def sameSpec(self,otherseq,unk=False):    ### Returns true if same spec (or SpecCode) or False if not
        '''
        Returns true if same spec (or SpecCode) or False if not.
        >> otherseq:Sequence Object.
        >> unk:bool [False] = value to be returned if either sequence is of unknown species.
        '''
        try:
            if not otherseq: return False
            for unktest in ['Unknown','UNK','UNKSP']:
                if unktest in [self.info['SpecCode'],otherseq.info['SpecCode']]: return unk
            if self.info['SpecCode'] == otherseq.info['SpecCode']: return True
            return False
        except:
            self.log.errorLog('Major problem with sameSpec()')
            return False
#########################################################################################################################
    def hasID(self,id=None,gnspacc=True,uniprot=True):  ### Checks ID, AccNum and secondary ID for match
        '''
        Checks ID, AccNum and secondary ID for match.
        >> id:str [None] = desired ID.
        >> gnspacc:boolean [True] = whether to try generating GnSpAcc formats for testing.
        >> uniprot:boolean [True] = whether to try making "ID (Acc)" combos for testing.
        << Returns True if ID found or False if not.
        '''
        if not id: return False
        if id in [self.info['ID'],self.info['AccNum'],self.list['Secondary ID'],self.shortName()]: return True
        fullacc = [self.info['AccNum']] + self.list['Secondary ID']
        if gnspacc:
            gnsp = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            for acc in fullacc:
                if id == '%s__%s' % (gnsp,acc): return True
        if uniprot:
            for acc in fullacc:
                if id in ['%s (%s)' % (self.info['ID'],acc),'%s_%s' % (acc,self.info['SpecCode'])]: return True
        return False
#########################################################################################################################
    def fudgeFactor(self,posdict={},case=False,gaps=True,wildcards=True):   ### Returns the necessary fudge factor to match posdict
        '''
        Returns the necessary fudge factor to match posdict. Will raise error if it cannot do it!
        >> posdict:dict = {pos (1 to L):sequence}
        >> case:bool [False] = whether to use self.dict['Case'] information.
        >> gaps:bool [True] = whether to leave gaps in the sequence
        >> wildcards:bool [True] = whether to account for masked/wildcard Xs in posdict sequence
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not posdict: return 0
            maxf = max(rje.sortKeys(posdict)[0],self.aaLen() - rje.sortKeys(posdict)[-1])
            f = 0   # Fudge factor
            sequence = self.getSequence(case,gaps)
            ### ~ [2] Check fudging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while f < maxf:
                ## ~ [2a] Try forwards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fudged = True
                for pos in posdict:
                    if wildcards: pos_re = string.replace(posdict[pos],'X','\S')
                    else: pos_re = posdict[pos]
                    if not re.search('^%s' % pos_re,sequence[pos-1+f:]):
                        fudged = False
                        break
                if fudged: return f
                ## ~ [2b] Try backwards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fudged = True
                for pos in posdict:
                    if wildcards: pos_re = string.replace(posdict[pos],'X','\S')
                    else: pos_re = posdict[pos]
                    if not re.search('^%s' % pos_re,sequence[pos-1-f:]):
                        fudged = False
                        break
                if fudged: return -f
                f += 1
        except: self.log.errorLog('Major problem with %s fudgeFactor()' % self.shortName())
        raise ValueError
#########################################################################################################################
    ### <4> ### Sequence calculation methods                                                                            #
#########################################################################################################################
    def disorder(self,returnlist=False,reset=False,ikey='Sequence'): ### Adds disorder object and predicts disorder. See rje_disorder.py for details.
        '''Adds disorder object and predicts disorder. See rje_disorder.py for details.'''
        ### ~ [1] ~ Check for existing disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        done = self.obj['Disorder'] and self.obj['Disorder'].list['ResidueDisorder'] and not reset
        ### ~ [2] ~ Add Object and calculate disorder, if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not done:
            self.obj['Disorder'] = rje_disorder.Disorder(self.log,self.cmd_list)
            done = self.obj['Disorder'].disorder(sequence=self.info[ikey],name=self.info['Name'])
        ### ~ [3] ~ Return list of disordered residues, or simply whether the disorder prediction worked ~~~~~~~~~~~~ ###
        if returnlist:
            if done: return self.obj['Disorder'].list['ResidueDisorder'][0:]    # List of scores?
            else: return []
        else: return done
#########################################################################################################################
    def isDisordered(self,pos=-1,reset=False): ### Returns whether a particular residue is disordered
        '''
        Returns whether a particular residue is disordered.
        >> pos:int [-1] = position of residue (0->(L-1)). If < 0 will return True/False for status only.
        >> reset:bool [False] = Reset disorder list before assessing.
        '''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if reset or 'IsDisordered' not in self.list or not self.list['IsDisordered']:
            if not self.disorder(reset=reset):
                self.errorLog('Disorder prediction failed for %s' % (self.shortName()),printerror=False)
                raise ValueError
            disorder = [True] * self.aaLen()
            self.list['IsDisordered'] = [False] * self.aaLen()
            for region in self.obj['Disorder'].list['RegionDisorder']:
                self.list['IsDisordered'] = self.list['IsDisordered'][:region[0]-1] + disorder[region[0]-1:region[1]] + self.list['IsDisordered'][region[1]:]
        ### ~ [2] ~ Return relevant information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if pos < 0: return len(self.list['IsDisordered']) > 0
        return self.list['IsDisordered'][pos]
#########################################################################################################################
    def gappedDisorder(self,gap='prev'):   ### Returns gapped disorder list
        '''Returns gapped disorder list.'''
        ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        disorder = self.disorder(returnlist=True)
        gapdis = []
        ### ~ [1] ~ Add gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for r in self.getSequence():
            if r == '-':
                if gapdis and gap == 'prev': gapdis.append(gapdis[-1])
                elif gap == 'prev': gapdis.append(disorder[0])
                else: gapdis.append(gap)
            else: gapdis.append(disorder.pop(0))
        return gapdis
#########################################################################################################################
    def globProportion(self,absolute=False):     ### Returns proportion defined as globular
        '''Returns proportion defined as globular.'''
        self.disorder()
        return self.obj['Disorder'].globProportion(absolute)
#########################################################################################################################
    def aaFreq(self,aafreq={},newkeys=True):    ### Adds to aafreq dictionary (if given) and returns
        '''
        Adds to aafreq dictionary (if given) and returns.
        >> aafreq:dictionary of {aa:freq(count)}
        >> newkeys:boolean [True] = whether to add new AA keys if missing from aafreq
        << aafreq:new dictionary of values
        '''
        try: return aaFreq(self.info['Sequence'],aafreq,newkeys)
        except:
            self.log.errorLog('Problem with rje_sequence.aaFreq(%s)' % self.shortName(),quitchoice=True)
            return aafreq
#########################################################################################################################
    def sixFrameTranslation(self):  ### Translates DNA in all six reading frames into 'Translation' dictionary
        '''Translates DNA in all six reading frames into 'Translation' dictionary.'''
        try:self.dict['Translation'] = sixFrameTranslation(self.info['Sequence'][0:])
        except: self.log.errorLog('Problem with Sequence.sixFrameTranslation()')
#########################################################################################################################
    ### <5> ### Sequence masking methods                                                                                #
#########################################################################################################################
    def maskAA(self,maskaa,mask='X',log=True):    ### Masks given residues by type
        '''
        Adds disorder object with prediction and masks disorder.
        >> maskaa:list of AAs to be masked
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        << returns number of AAs masked
        '''
        try:### ~ [1] ~ Straight masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0
            for aa in maskaa:
                mx += self.info['Sequence'].count(aa)
                self.info['Sequence'] = self.info['Sequence'].replace(aa,mask)
            if log: self.printLog('#MASK','AA Mask (%s): %s %s added to %s' % (string.join(maskaa,';'),rje.integerString(mx),mask,self.shortName()))
            return mx
        except: self.errorLog('Problem masking AAs in %s' % (self.shortName())); return 0
#########################################################################################################################
    def maskDisorder(self,inverse=False,mask='X',log=True,ikey='Sequence'):   ### Adds disorder object with prediction and masks disorder
        '''
        Adds disorder object with prediction and masks disorder.
        >> inverse:bool [False] = Masks out predicted ordered regions (i.e. keeps disorder only)
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        >> ikey:str ['Sequence'] = String key to be used for sequence data. (e.g. PreMask)
        '''
        try:### ~ [1] ~ Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deGap()    #!# Update to ignore gaps at some point? #!#
            if not self.disorder(ikey): return self.errorLog('%s Disorder prediction failed for %s' % (ikey,self.shortName()),printerror=False)
            oldseq = self.info['Sequence'][0:]
            newseq = mask * len(oldseq)
            prex = oldseq.count(mask)
            if inverse: (newseq,oldseq) = (oldseq,newseq)
            ### ~ [2] ~ Mask relevant regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for region in self.obj['Disorder'].list['RegionDisorder']: oldseq = oldseq[:region[0]-1] + newseq[region[0]-1:region[1]] + oldseq[region[1]:]
            #x#self.deBug('%s::%s\n%s' % (self.obj['Disorder'].list['RegionDisorder'],self.obj['Disorder'].list['RegionFold'],oldseq))               
            ### ~ [3] ~ Update sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0:
                if inverse: mtxt = '%s masked %d ' % (self.shortName(),len(self.obj['Disorder'].list['RegionFold']))
                else: mtxt = '%s masked %d dis' % (self.shortName(),len(self.obj['Disorder'].list['RegionDisorder']))
                self.printLog('#MASK','%sordered regions. (%d %s added.)' % (mtxt,maskx,mask),screen=log)
        except: self.errorLog('Problem masking disorder from %s' % (self.shortName()))
#########################################################################################################################
    def maskLowComplexity(self,lowfreq=5,winsize=10,mask='X',log=True):     ### Masks low complexity regions of sequence
        '''
        Masks low complexity regions of sequence.
        >> lowfreq:int = Number of same aas in window size to mask
        >> winsize:int = Size of window to consider
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        '''
        try:### Setup ###
            self.deGap()
            if lowfreq > winsize or lowfreq < 3: return
            oldseq = self.info['Sequence'][0:]
            prex = oldseq.count(mask)

            ### Mask ###
            for r in range(len(oldseq)):
                a = self.info['Sequence'][r]
                if a == mask: continue
                x = 1
                for i in range(1,winsize):
                    if (r+i) >= len(oldseq): break
                    if self.info['Sequence'][r+i] == a: x += 1
                    if x >= lowfreq:
                        oldseq = oldseq[:r+1] + string.replace(oldseq[r+1:r+i],a,mask) + oldseq[r+i:]
                        #X#self.deBug('%s => %s' % (self.info['Sequence'][r+1:r+i],string.replace(oldseq[r+1:r+i],a,mask)))
                        break
                
            ### Update ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s "low complexity" regions. (%d %s added.)' % (self.shortName(),maskx,mask),screen=log)
                    
        except: self.log.errorLog('Problem masking low complexity for %s' % self.shortName())
#########################################################################################################################
    def maskCase(self,case='lower',mask='X',log=True):   ### Masks sequence of given case
        '''
        Returns sequence.
        >> case:str ['lower'] = whether mask lower/upper case information.
        >> mask:str ['X'] = character to replace sequence with
        >> log:bool [True] = whether to log the amount of masking
        '''
        try:
            ### Setup ###
            if case.lower()[:1] == 'l': case = 'Lower'
            elif case.lower()[:1] == 'u': case = 'Upper'
            else: return False

            ### Mask ###
            prex = self.info['Sequence'].upper().count(mask.upper())
            sequence = self.info['Sequence'].upper()
            for (start,stop) in self.dict['Case'][case]: sequence = sequence[:start] + mask * (stop + 1 - start) + sequence[(stop+1):]
            self.info['Sequence'] = sequence

            ### Finish ###
            maskx = self.info['Sequence'].upper().count(mask.upper()) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s %s case regions. (%d %s added.)' % (self.shortName(),case,maskx,mask),screen=log)
        except:
            self.log.errorLog('Problem masking case "%s" for %s' % (case,self.shortName()))
#########################################################################################################################
    def maskRegion(self,maskregion=[],inverse=False,mask='X',log=True): ### Masks region(s) of protein
        '''
        Masks region of protein.
        >> maskregion:list [] = Pairs of positions to (inclusively) mask
        >> inverse:bool [False] = Masks out predicted ordered regions (i.e. keeps disorder only)
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        '''
        try:### ~ [1] ~ Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not maskregion: return 0
            maskregion = maskregion[0:]
            self.deGap()    #!# Update to ignore gaps at some point? #!#
            oldseq = self.info['Sequence'][0:]
            newseq = mask * len(oldseq)
            prex = oldseq.count(mask)
            if inverse: (newseq,oldseq) = (oldseq,newseq)
            ### ~ [2] ~ Mask relevant regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0
            while maskregion:
                start = maskregion.pop(0)
                try: end = maskregion.pop(0)
                except: end = -1
                if start < 0: start = self.aaLen() - start
                if end < 0: end = self.aaLen() - end
                end += 1
                oldseq = oldseq[:start] + newseq[start:end] + oldseq[end:]
                mx += 1
            ### ~ [3] ~ Update sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0: self.printLog('#MASK','%d region(s). (%d %s added.)' % (mx,maskx,mask),screen=log)
            return maskx
        except: self.errorLog('Problem with Region masking for %s' % (self.shortName()))
#########################################################################################################################
    def maskPosAA(self,maskdict={},mask='X',log=True):  ### Masks position-specific amino acids
        '''
        Masks position-specific amino acids. Returns and updates sequence.
        >> maskdict:dictionary of {pos:AAs} where pos is 1->L and AAs is a string of the AAs to mask
        >> mask:str ['X'] = character to replace sequence with
        >> log:bool [True] = whether to log the amount of masking
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = self.info['Sequence'].upper().count(mask.upper())
            sequence = self.info['Sequence'].upper()[0:]

            ### ~ [2] Mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pos in rje.sortKeys(maskdict):
                aas = maskdict[pos].upper()
                i = int(pos)
                if i > 0: i -= 1    # Can have backwards counting too and that's the same!
                try:
                    if aas in ['*','X','.'] or sequence[i] in rje.strList(aas): sequence = rje.strSub(sequence,i,i,mask)
                except: pass
            self.info['Sequence'] = sequence

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maskx = self.info['Sequence'].upper().count(mask.upper()) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s position-specific AAs. (%d %s added.)' % (self.shortName(),maskx,mask),screen=log)
        except: self.log.errorLog('Problem masking positions for %s' % (self.shortName()))
        return self.info['Sequence']
#########################################################################################################################
## End of SECTION II: Sequence Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def specCodeFromName(name):  ### Returns the species code from a sequence name
    '''
    Returns the species code from a sequence name. (Stripped down version of Sequence.extractDetails())
    >> name:str = Sequence name
    << spcode:str = Species code
    '''
    seqname = string.replace(name,'||','|')  #!# Dodgy!
    accnum = seqname
    spcode = 'UNK'
    if rje.regExp(re_gnspacc,seqname) and '|' not in string.split(seqname,'_')[0]:
        (gene,spcode,accnum) = rje.regExp(re_gnspacc,seqname)
    elif rje.regExp(re_gn_sp__acc,seqname) and '|' not in string.split(seqname,'_')[0]:
        (gene,spcode,accnum) = rje.regExp(re_gn_sp__acc,seqname)
    elif rje.regExp(re_uniprot,seqname):
        (db,accnum,gene,spcode) = rje.regExp(re_uniprot,seqname)
    elif rje.regExp(re_uniprot2,seqname):
        (gene,spcode,accnum) = rje.regExp(re_uniprot2,seqname)
    elif rje.regExp(re_ipi,seqname):
        (accnum,_taxa) = rje.regExp(re_ipi,seqname)
        if _taxa in ipi_taxa.keys():
            spcode = ipi_taxa[_taxa]
        else:
            spcode = 'UNK'
    elif rje.regExp(re_ncbi,seqname):
        (id,db,accnum) = rje.regExp(re_ncbi,seqname)
        if rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',seqname):
            [accnum,gene,spcode] = rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',seqname)
        else:
            if rje.matchExp('\[(.+)\]\s*$',seqname):
                species = rje.matchExp('\[(.+)\]\s*$',seqname)[0]
                while re.search('\[(.+)$',species):
                    species = rje.matchExp('\[(.+)$',species)[0]
                spcode = getSpecCode(species)
    elif rje.regExp(re_ncbi2,seqname):
        (id,db,accnum) = rje.regExp(re_ncbi2,seqname)
        if rje.matchExp('\[(.+)\]\s*$',seqname):
            species = rje.matchExp('\[(.+)\]\s*$',seqname)[0]
            while re.search('\[(.+)$',species): species = rje.matchExp('\[(.+)$',species)[0]
            spcode = getSpecCode(species)
    elif re_genbank.match(seqname): (db,accnum) = rje.regExp(re_genbank,seqname)
    elif re_unigene.match(seqname): accnum = rje.regExp(re_unigene,seqname)[0]
    elif re_flybase.match(seqname):
        (acc,gene,spec) = rje.regExp(re_flybase,seqname)
        spcode = '%sRO%s' % (spec[0].upper(),spec[1:3].upper())
    elif re_ensemblpep.match(seqname):
        accnum = rje.regExp(re_ensemblpep,seqname)[0]
        spcode = 'UNK'
        if accnum.find('ENSAPMP') == 0: spcode = 'APIME'
        elif accnum.find('ENSBTAP') == 0: spcode = 'BOVIN'
        elif accnum.find('ENSGALP') == 0: spcode = 'CHICK'
        elif accnum.find('ENSPTRP') == 0: spcode = 'PANTR'
        elif accnum.find('ENSCINP') == 0: spcode = 'CIOIN'
        elif accnum.find('ENSCAV') == 0: spcode = 'CIOSA'
        elif accnum.find('ENSCAFP') == 0: spcode = 'CANFA'
        elif accnum.find('CG') == 0: spcode = 'DROME'
        elif accnum.find('ENSDNOP') == 0: spcode = 'DASNO'
        elif accnum.find('ENSETEP') == 0: spcode = 'ECHTE'
        elif accnum.find('SINFRUP') == 0: spcode = 'FUGRU'
        elif accnum.find('ENSGACP') == 0: spcode = 'GASAC'
        elif accnum.find('ENSDARP') == 0: spcode = 'BRARE'
        elif accnum.find('ENSP0') == 0: spcode = 'HUMAN'
        elif accnum.find('ENSLAFP') == 0: spcode = 'LOXAF' # Elephant
        elif accnum.find('ENSANGP') == 0: spcode = 'ANOGA'
        elif accnum.find('ENSMUSP') == 0: spcode = 'MOUSE'
        elif accnum.find('ENSOCUP') == 0: spcode = 'RABIT'  # Rabbit
        elif accnum.find('ENSRNOP') == 0: spcode = 'RAT'
        elif accnum.find('GSTENP') == 0: spcode = 'TETNG'
        elif accnum.find('ENSXETP') == 0: spcode = 'XENLA'
        elif accnum.find('AAEL') == 0: spcode = 'AEDAE'
        elif seqname.find('chromosome:SGD') > 0: spcode = 'YEAST'
        elif seqname.find('chromosome:CEL') > 0: spcode = 'CAEEL'
        elif seqname.find('ENSMMUP') == 0: spcode = 'MACMU'
        elif seqname.find('scaffold:MMUL') > 0: spcode = 'MACMU'
        elif seqname.find('ENSMODP') == 0: spcode = 'MONDO'
        elif seqname.find('scaffold:BROAD') > 0: spcode = 'MONDO'
    elif re_pipe.match(seqname): spcode = rje.regExp(re_pipe,seqname)[1]

    ## Special Databases ##
    if seqname.find('TP_') >= 0 and seqname.find('{Treponema pallidum Nichols}') >= 0: spcode = 'TREPA'
    if spcode == 'UNK' and (seqname.find('Emiliania huxleyi cDNA') >= 0 or accnum[:4] == 'EHUX'): spcode = 'EMIHU'
    if seqname.find('jgi|Thaps3|') == 0: spcode = 'THAPS'

    return spcode
#########################################################################################################################
def extractNameDetails(name,callobj=None):   ### Extracts details from sequence name and returns as dictionary
    '''Extracts details from sequence name and returns as dictionary.'''
    ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    name = string.replace(name,'||','|')
    data = {'Name':name,'SpecCode':'UNK'}
    ### ~ [1] ~ General Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if re_plain.match(name):
        match = re_plain.match(name)
        data['Format'] = "plain"
        grp = match.groups()
        data['ID'] = grp[0]
        data['AccNum'] = grp[0]
        data['Description'] = grp[1]
        if rje.matchExp('gi:(\d+)',data['Description']): data['NCBI'] = rje.matchExp('gi:(\d+)',data['Description'])[0]
    elif re_acconly.match(name):
        match = re_acconly.match(name)
        data['Format'] = "acconly"
        grp = match.groups()
        data['ID'] = grp[0]
        data['AccNum'] = grp[0]
        data['Description'] = ''
    else:
        if callobj: callobj.printLog('#ERR','Major problem recognising sequence name "%s"!' % name)
        raise ValueError
    ### ~ [2] ~ Special Formats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        if re_gn_sp__acc.match(name) and '|' not in string.split(name,'_')[0]:
            data['Format'] = 'gn_sp__acc'
            [data['Gene'],data['SpecCode']] = string.split(name,'_')[:2]
            data['AccNum'] = string.join(string.split(string.split(name)[0],'__')[1:],'-')
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            if data['Gene'].upper() == data['Gene']:
                if name.find('acc:%s' % data['Gene']) > 0: data['DBase'] = 'trembl'
                elif data['Gene'] != data['AccNum']: data['DBase'] = 'sprot'
                else: data['DBase'] = 'trembl'
        elif re_gnspacc.match(name) and '|' not in string.split(name,'_')[0]:
            data['Format'] = 'gnspacc'
            match = re_gnspacc.match(name)
            grp = match.groups()
            data['Gene'] = grp[0]
            data['SpecCode'] = grp[1]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = grp[2]
            if data['Gene'].upper() == data['Gene']:
                if name.find('acc:%s' % data['Gene']) > 0: data['DBase'] = 'trembl'
                elif data['Gene'] != data['AccNum']: data['DBase'] = 'sprot'
                else: data['DBase'] = 'trembl'
        elif re_enst.match(name):
            details = string.split(name,'|')   # re_jgi = re.compile('^jgi\|Thaps3\|(\d+)\|(\S+)')
            data['DBase'] = data['Format'] = 'ENST'
            data['Gene'] = details[-1].lower()
            data['AccNum'] = details[0]
            data['Description'] = details[2]
            data['EnsG'] = details[1]
            data['SpecCode'] = details[0][3:6]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
        elif re_jgi.match(name):
            details = string.split(name,'|')   # re_jgi = re.compile('^jgi\|Thaps3\|(\d+)\|(\S+)')
            data['Gene'] = 'jgi'
            if details[1] == 'Thaps3':
                data['SpecCode'] = 'THAPS'
                data['Species'] = 'Thalassiosira pseudonana'
                if details[1][-2:] == 'bd': data['AccNum'] = 'TPSBD%s' % details[2]
                else: data['AccNum'] = 'TPSJGI%s' % details[2]
            elif details[1] == 'Emihu1':
                data['SpecCode'] = 'EMIHU'
                data['Species'] = 'Emiliania huxleyi'
                data['AccNum'] = 'EHUXJGI%s' % details[2]
            elif details[1] == 'Nemve1':
                data['SpecCode'] = 'NEMVE'
                data['Species'] = 'Nematostella vectensis'
                data['AccNum'] = 'NEMVEJGI%s' % details[2]
            else:
                data['SpecCode'] = details[1].upper()
                data['AccNum'] = 'JGI%s' % details[2]
            data['Description'] = details[3]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
        elif re_uniprot.match(name):
            data['Format'] = 'uniprot'
            detail = rje.regExp(re_uniprot,name)
            data['Gene'] = detail[2]
            data['SpecCode'] = detail[3]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[1]
            if data['Gene'] == data['AccNum']:
                data['DBase'] = 'trembl'
                data['Gene'] = 'tr'
            else:
                data['DBase'] = 'sprot'
        elif re_disprot.match(name):
            data['Format'] = 'DisProt'
            data['DBase'] = 'disprot'
            data['AccNum'] = re_disprot.match(name)[0]
            data['Gene'] = 'dp'
        elif re_hprd.match(name):
            data['DBase'] = data['Format'] = 'HPRD'
            data['ID'] = name
            data['AccNum'] = 'HPRD%s' % rje.regExp(re_hprd,name)[0]
            data['Gene'] = 'hprd'
            data['SpecCode'] = 'HUMAN'
        elif re_hprd2.match(name):
            data = rje.combineDict(data,{'DBase':'HPRD','Format':'HPRD_Raw','Gene':'hprd','SpecCode':'HUMAN'})
            (data['ID'],sv,data['AccNum']) = rje.regExp(re_hprd2,name)
        elif re_uniref.match(name):    # re.compile('^(UniRef\S+)_(\S+)\s.*\Tax=(\S.+)\sRepID=(\S+)')
            data['Format'] = 'uniref'
            detail = rje.regExp(re_uniref,name)
            if detail[3].find('_') > 0: [data['Gene'],data['SpecCode']] = string.split(detail[3],'_')
            else:
                data['Gene'] = 'uref90'
                data['SpecCode'] = getSpecCode(detail[2])
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[1]
            data['DBase'] = detail[0]
        elif re_uniprot2.match(name):
            data['Format'] = 'uniprot2'
            detail = rje.regExp(re_uniprot2,name)
            data['Gene'] = detail[0]
            data['SpecCode'] = detail[1]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[2]
            if data['Gene'] == data['AccNum']:
                data['DBase'] = 'trembl'
                data['Gene'] = 'tr'
            else:
                data['DBase'] = 'sprot'
        elif re_uniprot3.match(name):
            data['Format'] = 'uniprot3'
            detail = rje.regExp(re_uniprot3,name)
            data['Gene'] = detail[1]
            data['SpecCode'] = detail[2]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[0]
            if data['Gene'] == data['AccNum']:
                data['DBase'] = 'trembl'
                data['Gene'] = 'tr'
            else:
                data['DBase'] = 'sprot'
        elif re_unirefprot.match(name): # re.compile('^([A-Za-z0-9\-]+)\s+([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
            data['Format'] = 're_unirefprot'
            detail = rje.regExp(re_unirefprot,name)
            data['Gene'] = detail[1]
            data['SpecCode'] = detail[2]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[0]
            if data['Gene'] == data['AccNum']:
                data['DBase'] = 'trembl'
                data['Gene'] = 'tr'
            else:
                data['DBase'] = 'sprot'
            data['Description'] = string.join(string.split(data['Description'])[1:])
        elif re_ipi.match(name):
            data['DBase'] = data['Format'] = 'IPI'
            data['ID'] = rje.regExp(re_ipi,name)[0]
            data['AccNum'] = data['ID']
            _taxa = rje.regExp(re_ipi,name)[1]
            if _taxa in ipi_taxa.keys():
                data['SpecCode'] = ipi_taxa[_taxa]
            else:
                data['SpecCode'] = 'UNK'
        elif re_ipimod.match(name): # re.compile('^>IPI:(IPI\d+\.\d+)\s+Gene_Symbol=(\S+)')
            data['DBase'] = 'IPI'; data['Format'] = 'IPIMod'
            (data['ID'],gene) = rje.regExp(re_ipimod,name)
            data['AccNum'] = data['ID']
            data['SpecCode'] = 'HUMAN'
            if gene == '-': data['Gene'] = 'ipi'
            else: data['Gene'] = string.split(gene,';')[0]
        elif re_gnspacc.match(data['Description']): #!# This could cause problems - try other databases first?
            data['Format'] = 'gnspacc fragment'
            match = re_gnspacc.match(data['Description'])
            grp = match.groups()
            data['Gene'] = grp[0]
            data['SpecCode'] = grp[1]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = grp[2]
            if data['Gene'].upper() == data['Gene']:
                if name.find('acc:%s' % data['Gene']) > 0: data['DBase'] = 'trembl'
                elif data['Gene'] != data['AccNum']: data['DBase'] = 'sprot'
                else: data['DBase'] = 'trembl'
        elif re_ncbi.match(name):
            data['Format'] = 'NCBI'
            [data['NCBI'],data['DBase'],data['AccNum']] = rje.regExp(re_ncbi,name)
            data['Gene'] = data['DBase']
            if rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',data['ID']):
                [data['AccNum'],data['Gene'],data['SpecCode']] = rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',data['ID'])
                data['DBase'] = 'sprot'
            else:
                data['ID'] = '%s_%s' % (data['DBase'],data['AccNum'])
                if rje.matchExp('\[(.+)\]\s*$',name):
                    data['Species'] = rje.matchExp('\[(.+)\]\s*$',name)[0]
                    while re.search('\[(.+)$',data['Species']):
                        data['Species'] = rje.matchExp('\[(.+)$',data['Species'])[0]
                    data['SpecCode'] = getSpecCode(data['Species'])
            data['Description'] = string.split(name)
            data['Description'][0] = string.split(data['Description'][0],'|')[-1]
            data['Description'] = string.join(data['Description'])
            data['Description'] = 'gi:%s %s' % (data['NCBI'],data['Description'])
        elif re_ncbi2.match(name):
            data['Format'] = 'NCBI'
            [data['NCBI'],data['DBase'],data['AccNum']] = rje.regExp(re_ncbi2,name)
            data['Gene'] = data['DBase']
            data['ID'] = '%s_%s' % (data['DBase'],data['AccNum'])
            if rje.matchExp('\[(.+)\]\s*$',name):
                data['Species'] = rje.matchExp('\[(.+)\]\s*$',name)[0]
                while re.search('\[(.+)$',data['Species']):
                    data['Species'] = rje.matchExp('\[(.+)$',data['Species'])[0]
                data['SpecCode'] = getSpecCode(data['Species'])
            data['Description'] = 'gi:%s %s' % (data['NCBI'],data['Description'])
        elif re_genbank.match(name):
            data['Format'] = 'GenBank'
            [data['DBase'],data['ID']] = rje.regExp(re_genbank,name)
            data['AccNum'] = data['ID']
        elif re_unigene.match(name):
            data['DBase'] = data['Format'] = 'UniGene'
            data['ID'] = rje.regExp(re_unigene,name)[0]
            data['AccNum'] = data['ID']
        elif re_flybase.match(name):
            (acc,gene,spec) = rje.regExp(re_flybase,name)
            data['DBase'] = 'FlyBase'
            data['ID'] = data['AccNum'] = acc
            if gene[-3:] in ['-RA','-PA']: gene = gene[:-3]
            data['Gene'] = string.split(gene,'\\')[-1]
            data['SpecCode'] = '%sRO%s' % (spec[0].upper(),spec[1:3].upper())
        elif re_ensemblpep.match(name):
            data['Format'] = 'ensemblpep'
            detail = rje.regExp(re_ensemblpep,name)
            data['AccNum'] = detail[0]
            data['Gene'] = 'p'
            if detail[1].lower() == 'known': data['Gene'] = 'ens'
            elif detail[1].lower() == 'novel': data['Gene'] = 'nvl'
            elif detail[1].lower() == 'genscan': data['Gene'] = 'scan'
            data['SpecCode'] = 'UNK'
            if data['AccNum'].find('ENSAPMP') == 0: data['SpecCode'] = 'APIME'
            elif data['AccNum'].find('ENSBTAP') == 0: data['SpecCode'] = 'BOVIN'
            elif data['AccNum'].find('ENSGALP') == 0: data['SpecCode'] = 'CHICK'
            elif data['AccNum'].find('ENSPTRP') == 0: data['SpecCode'] = 'PANTR'
            elif data['AccNum'].find('ENSCINP') == 0: data['SpecCode'] = 'CIOIN'
            elif data['AccNum'].find('ENSCAFP') == 0: data['SpecCode'] = 'CANFA'
            elif data['AccNum'].find('CG') == 0: data['SpecCode'] = 'DROME'
            elif data['AccNum'].find('FBpp') == 0: data['SpecCode'] = 'DROME'
            elif data['AccNum'].find('SINFRUP') == 0: data['SpecCode'] = 'FUGRU'
            elif data['AccNum'].find('ENSDARP') == 0: data['SpecCode'] = 'DANRE'
            elif data['AccNum'].find('ENSP0') == 0: data['SpecCode'] = 'HUMAN'
            elif data['AccNum'].find('ENSANGP') == 0: data['SpecCode'] = 'ANOGA'
            elif data['AccNum'].find('AGAP') == 0: data['SpecCode'] = 'ANOGA'
            elif name.find('group:AMEL') > 0: data['SpecCode'] = 'YEAST'
            elif data['AccNum'].find('ENSMUSP') == 0: data['SpecCode'] = 'MOUSE'
            elif data['AccNum'].find('ENSRNOP') == 0: data['SpecCode'] = 'RAT'
            elif data['AccNum'].find('GSTENP') == 0: data['SpecCode'] = 'TETNG'
            elif data['AccNum'].find('ENSXETP') == 0: data['SpecCode'] = 'XENLA'
            elif name.find('chromosome:SGD') > 0: data['SpecCode'] = 'YEAST'
            elif name.find('chromosome:CEL') > 0: data['SpecCode'] = 'CAEEL'
            elif name.find('chromosome:WS') > 0: data['SpecCode'] = 'CAEEL'
            elif name.find('ENSMMUP') == 0: data['SpecCode'] = 'MACMU'
            elif name.find('scaffold:MMUL') > 0: data['SpecCode'] = 'MACMU'
            elif name.find('ENSMODP') == 0: data['SpecCode'] = 'MONDO'
            elif name.find('scaffold:BROAD') > 0: data['SpecCode'] = 'MONDO'
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
        elif re_pipe.match(name):
            data['Format'] = 'pipe'
            detail = rje.regExp(re_pipe,name)
            data['Gene'] = detail[0]
            data['SpecCode'] = detail[1]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = detail[2]
            if data['Gene'] == data['AccNum']:
                data['DBase'] = 'trembl'
                data['Gene'] = 'tr'
            else: data['DBase'] = 'sprot'
        elif re_db_pipe.match(name):
            data['Format'] = 'db_pipe'
            (data['DBase'],data['AccNum']) = rje.regExp(re_db_pipe,name)
            data['Gene'] = data['DBase'].lower()
            data['ID'] = data['AccNum']
        elif re_gnspec.match(name):
            data['Format'] = 'gnspec'
            detail = rje.regExp(re_gnspec,name)
            data['Gene'] = detail[0]
            data['SpecCode'] = detail[1]
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
            data['AccNum'] = data['ID']
            data['DBase'] = 'custom'
        elif re_tigr.match(name):
            data['DBase'] = data['Format'] = 'TIGR'
            (data['AccNum'],data['Species']) = rje.matchExp(re_tigr,name)
            while rje.matchExp('\S+\|(\S+)',data['AccNum']):
                data['AccNum'] = rje.matchExp('\S+\|(\S+)',data['AccNum'])[0]
            data['Gene'] = 'p'
            data['SpecCode'] = getSpecCode(data['Species'])
            data['ID'] = '%s_%s' % (data['Gene'],data['SpecCode'])
        ## ~ [2a] ~ Additional Ensembl processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if 'Gene' not in data: data['Gene'] = data['AccNum']
        if data['Gene'] == 'ens': data['DBase'] = 'ens_known'
        elif data['Gene'] == 'nvl': data['DBase'] = 'ens_novel'
        elif data['Gene'] == 'scan': data['DBase'] = 'ens_scan'
        elif data['Gene'] == 'ipi': data['DBase'] = 'IPI'
        elif data['Gene'] == 'ref': data['DBase'] = 'ens_known'
        #!# Add more Pattern matching and the like for different formats!
    except:
        if callobj: callobj.errorLog('Advanced name ("%s") detail extraction error' % name)
    return data
#########################################################################################################################
def getSpecCode(species):   ### Returns spec_code for given species
    '''Returns spec_code for given species. This should be moved later and expanded to allow to read from speclist.'''
    spec_dic = {'Mus musculus':'MOUSE',
                'Homo sapiens':'HUMAN',
                'Gallus gallus':'CHICK',
                'Rattus norvegicus':'RAT',
                'Bos taurus':'BOVIN',
                'Danio rerio':'BRARE',
                'Saccharomyces cerevisiae':'YEAST',
                'Human immunodeficiency virus 1':'HIV1','Human immunodeficiency virus 2':'HIV2'
                }
    if species in spec_dic.keys(): return spec_dic[species]
    tax = string.split(species.upper()); t = 1
    while t < len(tax):
        if tax[t][:1] == '(': tax.pop(t)
        else: t += 1
    if len(tax) > 1: return tax[0][:3] + tax[1][:2]
    else: return string.replace(species,' ','').upper()[:10]
#########################################################################################################################
def eisenbergHydropathy(sequence,returnlist=False):  ### Returns the Eisenberg Hydropathy for the sequence
    '''
    Returns the Eisenberg Hydropathy for the sequence.
    >> sequence:str = AA sequence
    >> returnlist:boolean [False] = Returns list of hydropathies rather than total
    << hyd:float = Eisenberg Hydropathy for the sequence
    '''
    ### Setup ###
    eishyd = {'A':0.62, 'R':-2.53, 'N':-0.78, 'D':-0.9, 'C':0.29, 'Q':-0.85, 'E':-0.74, 'G':0.48, 'H':-0.4, 'I':1.38,
              'L':1.06, 'K':-1.5, 'M':0.64, 'F':1.19, 'P':0.12, 'S':-0.18, 'T':-0.05, 'W':0.81, 'Y':0.26, 'V':1.08,
              '-':0, 'X':0}
    eislist = []
    ### Calculate ###
    for aa in sequence: eislist.append(eishyd[aa])
    if returnlist: return eislist
    return sum(eislist)
#########################################################################################################################
def aaFreq(sequence,aafreq={},newkeys=True):    ### Adds to aafreq dictionary (if given) and returns
    '''
    Adds to aafreq dictionary (if given) and returns.
    >> sequence:str = Sequence for AAFreq calculation
    >> aafreq:dictionary of {aa:freq(count)}
    >> newkeys:boolean [True] = whether to add new AA keys if missing from aafreq
    << aafreq:new dictionary of values
    '''
    for aa in sequence:
        if aafreq.has_key(aa): aafreq[aa] += 1
        elif newkeys: aafreq[aa] = 1.0
    return aafreq
#########################################################################################################################
def codons(sequence,codonfreq={},newkeys=True,rna=True,code_only=True): ### Adds to codonfreq dictionary (if given) and returns
    '''
    Adds to codonfreq dictionary (if given) and returns.
    >> sequence:str = Sequence for AAFreq calculation
    >> codonfreq:dictionary of {codon:freq(count)}
    >> newkeys:boolean [True] = whether to add new keys if missing from codonfreq
    >> rna:bool [True] = whether to convert to RNA (True) or DNA (False)
    >> code_only:bool [True] = whether to only allow codons in (RNA!) Genetic Code (True)
    << codonfreq:new dictionary of values
    '''
    if rna: sequence = string.replace(sequence.upper(),'T','U')
    else: sequence = string.replace(sequence.upper(),'U','T')
    while sequence:
        codon = sequence[:3]; sequence = sequence[3:]
        if len(codon) < 3: break
        if rna and code_only and codon not in genetic_code: continue
        if codonfreq.has_key(codon): codonfreq[codon] += 1
        elif newkeys: codonfreq[codon] = 1
    return codonfreq
#########################################################################################################################
def trypDigest(sequence=''):    ### Returns trypsin digestion of sequence
    '''Returns trypsin digestion of sequence.'''
    digest = string.join(string.split(sequence.upper(),'K'),'K ')
    digest = string.join(string.split(digest,'R'),'R ')
    digest = string.split(digest)
    return digest
#########################################################################################################################
def MWt(sequence='',type='raw'):   ### Returns Molecular Weight of Sequence
    '''Returns Molecular Weight of Sequence.'''
    mwt_raw = {'A':89, 'V':117, 'L':131, 'I':131, 'P':115,
              'F':165, 'W':204, 'M':149, 'G':75, 'S':105,
              'T':119, 'C':121, 'Y':181, 'N':132, 'Q':146,
              'D':133, 'E':147, 'K':146, 'R':174, 'H':155,
              'X':136.75}
    mwt_mono = {'A':71.03711,'R':156.10111,'N':114.04293,'D':115.02694,'C':103.00919,'E':129.04259,'Q':128.05858,
                'G':57.02146,'H':137.05891,'I':113.08406,'L':113.08406,'K':128.09496,'M':131.04049,'F':147.06841,
                'P':97.05276,'S':87.03203,'T':101.04768,'W':186.07931,'Y':163.06333,'V':99.06841,
                'U':150.953636,'O':237.147727,'X':136.75}
    mwt_ave = {'A':71.0788,'R':156.1875,'N':114.1038,'D':115.0886,'C':103.1388,'E':129.1155,'Q':128.1307,
               'G':57.0519,'H':137.1411,'I':113.1594,'L':113.1594,'K':128.1741,'M':131.1926,'F':147.1766,'P':97.1167,
               'S':87.0782,'T':101.1051,'W':186.2132,'Y':163.1760,'V':99.1326,
               'U':150.0388,'O':237.3018,'X':136.75}
    if type == 'mono': mwtdic = mwt_mono
    elif type == 'ave': mwtdic = mwt_ave
    else: mwtdic = mwt_raw
    _mwt = 0.0
    for aa in sequence.upper():
        if aa in mwtdic.keys():
            _mwt += mwtdic[aa]
            if type not in ['mono','ave']: _mwt -= 18  # H2O
    if _mwt: _mwt += 18
    return _mwt
#########################################################################################################################
def chargeDict(sequence,callobj=None):   ### Performs absolute, net and charge balance calculations 
    '''
    Performs absolute, net and charge balance calculations.
    >> sequence:str = sequence to calculate stats on.
    >> callobj:Object = calling object for error messages etc.
    << dictionary of {Type:Value)
    '''
    try:
        charge = []
        for a in sequence: # Motif
            if a in ['K','R']:
                charge.append(1)
            elif a in ['D','E']:
                charge.append(-1)
            else:
                charge.append(0)
        return {'AbsChg':(charge.count(1) + charge.count(-1)),
                'NetChg':sum(charge),
                'BalChg':(sum(charge[:int(len(charge)/2)]) - sum(charge[-int(len(charge)/2):]))}
    except:
        if callobj:
            callobj.log.errorLog('Error in rje_sequence.chargeDict()')
            return {}
        raise
#########################################################################################################################
def peptideDetails(sequence,callobj=None):  ### Returns OK if sequence alright, or warning if bad aa combos
    '''
    Returns OK if sequence alright, or warning if bad aa combos.
    >> sequence:str = sequence to calculate stats on.
    >> callobj:Object = calling object for error messages etc.
    << string of peptide assessment
    '''
    try:
        bad_combo = []
        for combo in ['DP','DC','DG','NG','NS','PPP']:
            if string.strip(sequence.upper(),'-').find(combo) >= 0: bad_combo.append(combo)
        if bad_combo: return 'Warning: %s' % string.join(bad_combo,', ')
        else: return 'OK'
    except:
        if callobj:
            callobj.log.errorLog('Error in rje_sequence.peptideDetails()')
            return {}
        raise
#########################################################################################################################
def mapGaps(inseq,gapseq,callobj=None):  ### Returns inseq with gaps inserted as found in gapseq
    '''Returns inseq with gaps inserted as found in gapseq.'''
    try:
        newseq = ''
        while len(newseq) < len(gapseq):
            if gapseq[len(newseq)] == '-': newseq += '-'
            elif inseq:
                newseq += inseq[:1]
                inseq = inseq[1:]
            else:
                if callobj: callobj.errorLog('Sequence length mismatch during mapGaps()',printerror=False)
                while len(newseq) < len(gapseq): newseq += '-'
        return newseq
    except:
        print newseq, inseq
        print gapseq
        return newseq
#########################################################################################################################
def geneticCode(transl='1',warnobj=None):   ### Returns a specific genetic code.
    '''
    Returns a specific genetic code.
    >> transl:str ['1'] = NCBI translation table. Default (1) = "Universal" genetic code
    '''
    transl = str(transl)
    if transl in transl_table: return transl_table[transl]
    if warnobj: warnobj.warnLog('Translation table %s not (yet) implemented. Contact the author. Standard code used.' % transl)
    return genetic_code
#########################################################################################################################
def dna2prot(dnaseq,case=False,transl='1',warnobj=None):   ### Returns a protein sequence for a given DNA sequence
    '''Returns a protein sequence for a given DNA sequence.'''
    gencode = geneticCode(transl,warnobj)
    prot = ['']
    i = 0
    while i < len(dnaseq):
        codon = dnaseq[i:i+3]
        i += 3
        lowc = case and codon == codon.lower()
        if lowc:
            try: prot.append(gencode[string.replace(codon.upper(),'T','U')].lower())
            except: prot.append('x')
        else:
            try: prot.append(gencode[string.replace(codon.upper(),'T','U')])
            except: prot.append('X')
    return string.join(prot,'')
#########################################################################################################################
def codonKs(codon,mutdict={}): ### Returns the proportion of substitutions that would be synonymous
    '''
    Returns the proportion of substitutions that would be synonymous.
    '''
    codon = string.replace(codon.upper(),'T','U')
    try: aa = genetic_code[codon]
    except: raise ValueError('%s not an acceptable codon' % codon)
    Ks = 0.0
    for i in range(3):
        for n in 'ACGU':
            if codon[i] == n: continue
            if genetic_code[codon[:i]+n+codon[i+1:]] == aa:
                if mutdict:
                    n1 = codon[i]
                    if n1 == 'U': n1 = 'T'
                    n2 = n
                    if n2 == 'U': n2 = 'T'
                    Ks += mutdict[n1][n2]
                else:
                    Ks += 1.0
    return Ks / 9.0
#########################################################################################################################
def mutationDict(mutations=[],gaps=False,nosub=False,callobj=None):    ### Generates dictionary of relative mutation frequencies
    '''
    Generates dictionary of relative mutation frequencies. These are 1.0 if all mutations are equally likely.
    >> mutations:list [] = List of (From,To[,Count) mutation tuples
    >> gaps:bool [False] = Whether to include indels (mutations to/from gaps)
    >> nosub:bool [False] = Whether to include "mutations" that do not result in a substitution
    >> callobj:Object [None] = Calling object to control log output
    '''
    try:### ~ [0] Setup the dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        mutdict = {}
        for n1 in 'AGCT-':
            if n1 == '-' and not gaps: continue
            mutdict[n1] = {}
            for n2 in 'AGCT-':
                if n2 == '-' and not gaps: continue
                if n1 == n2 and not nosub: continue
                if mutations or n1 == n2: mutdict[n1][n2] = 0.0
                else: mutdict[n1][n2] = 1.0
        ### ~ [1] Normalise dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for mut in mutations:
            n1 = mut[0].upper()
            n2 = mut[1].upper()
            if n1 == n2 and not nosub: continue
            if n1 not in mutdict: continue
            if n2 not in mutdict: continue
            nn = 1
            if len(mut) > 2: nn = mut[2]
            mutdict[n1][n2] += nn
        ### ~ [2] Normalise dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for n1 in 'AGCT-':
            if n1 == '-' and not gaps: continue
            dsum = float(sum(mutdict[n1].values())) / len(mutdict[n1])
            for n2 in 'AGCT-':
                if n2 == '-' and not gaps: continue
                if n1 == n2 and not nosub: continue
                mutdict[n1][n2] = mutdict[n1][n2] / dsum
        return mutdict
    except:
        if callobj: callobj.errorLog('mutationFreqDict() error')
        raise
#########################################################################################################################
def sequenceKs(sequence,callobj=None,ksdict={},mutdict={}):    ### Returns the proportion of possible synonymous substitutions
    '''Returns the proportion of possible synonymous substitutions.'''
    if not ksdict and callobj and 'Ks' in callobj.dict: ksdict = callobj.dict['Ks']
    if not ksdict: ksdict = kSDict(callobj,mutdict)
    ks = 0.0; ki = 0.0
    sequence = string.replace(sequence.upper(),'T','U')
    while sequence:
        codon = sequence[:3]
        sequence = sequence[3:]
        if len(codon) < 3: break
        elif codon not in genetic_code: continue
        ks += ksdict[codon]
        ki += 1
    return ks/ki
#########################################################################################################################
def kSDict(callobj=None,mutdict={}):   ### Returns the Ks frequency for each codon as a dictionary. Adds to callobj if given.
    '''Returns the Ks frequency for each codon as a dictionary. Adds to callobj if given.'''
    ksdict = {}
    for codon in genetic_code: ksdict[codon] = codonKs(codon,mutdict)
    if callobj: callobj.dict['Ks'] = ksdict
    return ksdict
#########################################################################################################################
def complement(dnaseq,rna=False):  ### Returns the complement of the DNA sequence given (mixed case)
    '''Returns the complement of the DNA sequence given.'''
    revcomp = dnaseq
    pairs = [('C','G'),('A','T')]
    if rna: revcomp = string.replace(revcomp,'U','T')
    if rna: revcomp = string.replace(revcomp,'u','t')
    for (n1,n2) in pairs:
        revcomp = string.replace(revcomp,n1,'!')
        revcomp = string.replace(revcomp,n2,n1)
        revcomp = string.replace(revcomp,'!',n2)
        n1 = n1.lower(); n2 = n2.lower()
        revcomp = string.replace(revcomp,n1,'!')
        revcomp = string.replace(revcomp,n2,n1)
        revcomp = string.replace(revcomp,'!',n2)
    if rna: revcomp = string.replace(revcomp,'T','U')
    if rna: revcomp = string.replace(revcomp,'t','u')
    return revcomp
#########################################################################################################################
def reverseComplement(dnaseq,rna=False):  ### Returns the reverse complement of the DNA sequence given (mixed case)
    '''Returns the reverse complement of the DNA sequence given.'''
    revcomp = rje.strReverse(dnaseq)
    pairs = [('C','G'),('A','T')]
    if rna: revcomp = string.replace(revcomp,'U','T')
    if rna: revcomp = string.replace(revcomp,'u','t')
    for (n1,n2) in pairs:
        revcomp = string.replace(revcomp,n1,'!')
        revcomp = string.replace(revcomp,n2,n1)
        revcomp = string.replace(revcomp,'!',n2)
        n1 = n1.lower(); n2 = n2.lower()
        revcomp = string.replace(revcomp,n1,'!')
        revcomp = string.replace(revcomp,n2,n1)
        revcomp = string.replace(revcomp,'!',n2)
    if rna: revcomp = string.replace(revcomp,'T','U')
    if rna: revcomp = string.replace(revcomp,'t','u')
    return revcomp
#########################################################################################################################
def OLDreverseComplement(dnaseq,rna=False):  ### Returns the reverse complement of the DNA sequence given (upper case)
    '''Returns the reverse complement of the DNA sequence given.'''
    revcomp = ''
    #!# Use rje.strReverse() followed by string.replace!
    repdict = {'G':'C','C':'G','T':'A','A':'T'}
    if rna: repdict['A'] = 'U'
    for aa in string.replace(dnaseq.upper(),'U','T'):   # Convert RNA to DNA
        try: revcomp = '%s%s' % (repdict[aa],revcomp)    
        except: revcomp = '%s%s' % ('N',revcomp)    
    return revcomp 
#########################################################################################################################
def sixFrameTranslation(dnaseq):  ### Translates DNA in all six reading frames into dictionary
    '''Translates DNA in all six reading frames into dictionary.'''
    sixrf = {}
    for rf in range(3): sixrf[rf+1] = dna2prot(dnaseq[rf:])
    dnaseq = reverseComplement(dnaseq)
    for rf in range(3): sixrf[-(rf+1)] = dna2prot(dnaseq[rf:])
    return sixrf
#########################################################################################################################
def threeFrameTranslation(dnaseq,minpoly=0):  ### Translates DNA in all three reading frames into dictionary
    '''Translates DNA in all six reading frames into dictionary.'''
    a = 1
    while dnaseq and minpoly >0 and dnaseq[-a].upper() == 'A': a += 1
    a -= 1
    if a >= minpoly and minpoly > 0: dnaseq = dnaseq[:-a]
    rf3 = {}
    for rf in range(3): rf3[rf+1] = dna2prot(dnaseq[rf:])
    return rf3
#########################################################################################################################
def bestORF(protseq,startm=False,nonx=True):  ### Returns longest ORF (or first of equal length)
    '''
    Returns longest ORF (or first of equal length).
    >> sequence:str = Protein sequence, i.e. translation of DNA/RNA.
    >> startm:bool [False] = Whether the ORF should start with a methionine
    >> nonx:bool [True] = Whether ORFs should only be assessed in terms of their non-X content.
    '''
    (bestorf,bestlen) = ('',0)
    for orf in string.split(protseq,'*'):
        if startm and orf.count('M') < 1: continue
        if startm and orf.find('M') >= 0: orf = orf[orf.find('M'):]
        orflen = len(orf)
        if nonx: orflen -= string.count(orf.upper(),'X')
        if orflen > bestlen: (bestorf,bestlen) = (orf,orflen)
    return bestorf
#########################################################################################################################
def estTranslation(dnaseq,minpoly=10,fwdonly=False):  ### Returns translations of EST into protein reading frames as dictionary
    '''
    Returns translation of EST into protein reading frames.
    >> dnaseq:str = DNA sequence to translate
    >> minpoly:str = Min. length of poly-A or poly-T to be recognised and removed.
    '''
    dnaseq = dnaseq.upper()
    t = 0       # Look for run of 5' Ts (i.e. rev complement of poly-A tail
    while dnaseq and dnaseq[t] == 'T' and not fwdonly: t += 1
    a = 1
    while dnaseq and dnaseq[-a] == 'A': a += 1
    a -= 1
    poly = ''
    if a >= minpoly and minpoly > 0:
        dnaseq = dnaseq[:-a]
        poly = 'A'
    elif t >= minpoly and minpoly > 0:
        dnaseq = dnaseq[t:]
        poly = 'T'
    sixrf = sixFrameTranslation(dnaseq)
    if fwdonly or poly in ['A','T']:
        for rf in rje.sortKeys(sixrf):
            if rf > 0  and poly == 'T': sixrf.pop(rf)
            elif rf < 0  and (fwdonly or poly == 'A'): sixrf.pop(rf)
            elif poly: sixrf[rf] = '%s*' % string.join(string.split(sixrf[rf],'*')[:-1],'*')
    return sixrf
#########################################################################################################################
def estTrunc(dnaseq,minpoly=10,fwdonly=False):  ### Returns truncation of EST along with poly-AT as tuple (5',seq,3')
    '''
    Returns translation of EST into protein reading frames.
    >> dnaseq:str = DNA sequence to translate
    >> minpoly:str = Min. length of poly-A or poly-T to be recognised and removed.
    '''
    dnaseq = dnaseq.upper()
    t = 0       # Look for run of 5' Ts (i.e. rev complement of poly-A tail
    while dnaseq[t] == 'T' and not fwdonly: t += 1
    a = 1
    while dnaseq[-a] == 'A': a += 1
    a -= 1
    poly = ''
    if a >= minpoly and minpoly >= 0: return ('',dnaseq[:-a],dnaseq[-a:])
    elif t >= minpoly and minpoly >= 0: return (dnaseq[:t],dnaseq[t:],'')
    return ('',dnaseq,'')
#########################################################################################################################
def caseDict(sequence): ### Return dictionary of {'Upper':[(start,end)],'Lower':[(start,end)]}
    '''Return dictionary of {'Upper':[(start,end)],'Lower':[(start,end)]}.'''
    case = {'Upper':[],'Lower':[]}      # Dictionary
    r = 0                               # Residue position
    for block in re.findall(re.compile('[A-Z\-\*]+|[a-z\-\*]+'),sequence):
        ckey = 'Lower'
        if block == block.upper(): ckey = 'Upper'
        case[ckey].append((r,r+len(block)-1))
        r += len(block)
    return case
#########################################################################################################################
def maskLowComplexity(sequence,lowfreq=5,winsize=10,mask='X'):     ### Masks low complexity regions of sequence
    '''
    Masks low complexity regions of sequence, keeping "outside" residues.
    >> sequence:str = Sequence to be masked. (Should not have gaps.)
    >> lowfreq:int [5] = Number of same aas in window size to mask
    >> winsize:int [10] = Size of window to consider
    >> mask:str ['X'] = character to use for masking
    << maskseq:str = returns masked sequence in upper case.
    '''
    ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if lowfreq > winsize or lowfreq < 3: return sequence    # No masking
    if len(mask) != 1: raise ValueError('Must have single mask=X character for maskLowComplexity()!')
    sequence = sequence.upper()
    mask = mask.upper()
    maskseq = sequence[0:]
    ### ~ [2] Mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    for r in range(len(sequence)):
        a = sequence[r]
        if a == mask: continue
        x = 1
        for i in range(1,winsize):
            if (r+i) >= len(sequence): break
            if sequence[r+i] == a: x += 1
            if x >= lowfreq:
                maskseq = maskseq[:r+1] + string.replace(maskseq[r+1:r+i],a,mask) + maskseq[r+i:]
                break
    if len(maskseq) != len(sequence): raise ValueError('Masked sequence length != input sequence length!')
    return maskseq
#########################################################################################################################
def maskPosAA(sequence,maskdict={},mask='X'):  ### Masks position-specific amino acids
    '''
    Masks position-specific amino acids. Returns masked sequence.
    >> sequence:str = Sequence to be masked. (Should not have gaps.)
    >> maskdict:dictionary of {pos:AAs} where pos is 1->L and AAs is a string of the AAs to mask
    >> mask:str ['X'] = character to replace sequence with
    << maskseq:str = returns masked sequence in upper case.
    '''
    ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sequence = sequence.upper()[0:]
    if len(mask) != 1: raise ValueError('Must have single mask=X character for maskPosAA()!')
    ### ~ [2] Mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    for pos in rje.sortKeys(maskdict):
        aas = maskdict[pos].upper()
        i = int(pos)
        if i > 0: i -= 1    # Can have backwards counting too and that's the same!
        try:
            if aas in ['*','X','.'] or sequence[i] in aas: sequence = rje.strSub(sequence,i,i,mask)
        except: pass
    ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    return sequence
#########################################################################################################################
def maskAA(sequence,maskaa,mask='X'):    ### Masks given residues by type
    '''
    Adds disorder object with prediction and masks disorder.
    >> maskaa:list of AAs to be masked
    >> mask:str ['X'] = character to use for masking
    << returns masked sequence (UC)
    '''
    ### ~ [1] ~ Straight masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sequence = sequence.upper()     #!# Could make this a maskSetup method that raises exception?
    if len(mask) != 1: raise ValueError('Must have single mask=X character for maskPosAA()!')
    for aa in maskaa: sequence = sequence.replace(aa,mask)
    return sequence
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: print 'Not for standalone running!'
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV
#########################################################################################################################
