#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       ExTATIC
Description:  Extensions and Truncations from Alternative Translation Initiation Codons
Version:      0.2.0
Last Edit:    28/01/15
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    ExTATIC predicts alternative initiation sites in cDNA (mRNA) sequences based on a set of rules regarding start codon
    efficiency. Input is in the form of two sequence files:

    1. A file of cDNA sequences (cdna=FASFILE).

    2. A file of CDS sequences (cds=FASFILE).

    For full functionality, Biomart Ensembl downloads are recommended. Sequence names should be in the format:

    >TranscriptID|GeneID[|GeneSymbol][|Description]

    The GeneID will be used to match different transcripts from the same gene. GeneSymbol is used purely for output and
    could be any external database xref of choice. Likewise, Description is purely for interpretation of results and is
    not used by the program. If fewer than four fields are found by splitting on '|', the third will be assumed to be the
    description and GeneID will be duplicated for GeneSymbol. If fewer than three, GeneID will also be duplicated for
    Description. If fewer that two, the first word of each sequence name will be used for both TranscriptID and GeneID
    (and GeneSymbol) with the remaining name forming the Description. If the Description contains "Gene:X" and no GeneID
    or GeneSymbol is given, X will be used for these values. Additional naming formats can be added on request.

    For stringency, 5' flanking regions will NOT be used to expand UTRs where missing/short: this should be performed by
    an upstream program/analysis if necessary.

    The context=DICT provides IUPAC contexts for AUG start codons and their strengths. Each context must contain XXX that
    indicates the start codon position. DNA or RNA notation can be used but contexts must not clash, i.e. each codon must
    uniquely match one context. If no contexts are provided, default values will be used:
    * XXXG = Strong
    * RnnXXXH = Mid ([AG]nnXXX[ACU])
    * CnnXXXU = Mid
    * UnnXXXC = Mid
    * CnnXXXM = Weak (CnnXXX[AC])
    * UnnXXXW = Weak (UnnXXX[AU])
    * ^XXXH = Unknown (XXX[ACU] at start of sequence)
    * ^nXXXH = Unknown (nXXX[ACU] at start of sequence)
    * ^nnXXXH = Unknown (nnXXX[ACU] at start of sequence)

    Alternative codons, given by altcodons=LIST, can match any contexts given by altcontext=LIST. By default, only
    "Strong" contexts are considered, for CTG and GTG codons. Annotated start sites that do not meet any allowed context
    (i.e. are non-canonical codons not in altcodons=LIST) will be given a context that is just the start codon.

    ### ~ ORFs and ORFTypes ~ ###
    The main output from ExTATIC is a series of tables identifying possible start codons and the corresponding Open
    Reading Frames (ORFs). A fasta file of predicted ORFs is also output. In this file, protein sequences are in lower
    case upto the first "Strong" ATG start codon. Upstream predicted start codon positions are also given in upper case.

    AIC are classified according to their position and reading frame relative to the annotated start codon (if CDS are
    given). The orftypes=LIST option sets which set of AIC will be returned. ORFs are classified according to their most
    5' AIC.

    * uORF = Upstream ORF that terminates upstream of the annotated start.
    * eORF = Extended annotated ORF. (Upstream in-frame start site with no stop codon before the annotated start.)
    * oORF = Upstream ORF that overlaps the annotated start.
    * aORF = Annotated ORF start site.
    * tORF = Truncated annotated ORF. (Downstream in-frame start site before stop codon or "Strong" ATG.
    * dORF = Downstream ORF that starts before a stop codon or "Strong" ATG in annotated ORF.

    ### ~ Single Sequence Analysis ~ ###
    If singleseq=X is given then ExTATIC will analyse a single sequence only. SeqIn and CDSIn should be sequence strings
    rather than fasta files. The sequence name will be taken from singleseq=X. Use append=T to add results to previous
    singleseq analyses and set basefile=X. (Default basefile for singleseq analysis is the first word from singleseq.)

    ### ~ REST Output ~ ###
    For details of available REST output for ExTATIC, see: http://rest.slimsuite.unsw.edu.au/extatic&rest=outfmt

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : cDNA input seqence file (*.cdna.fas). See docs for naming advice. [ExTATIC.cdna.fas]
    cdna=FASFILE    : Alternative cDNA input sequence file option.
    cds=FASFILE     : CDS input sequence file. (or cdsin=FASFILE) [*.cds.fas based on cdna or singleseq]
    singleseq=X     : Analyse a single sequence only (named X). cDNA and CDS should be sequence strings. [None]

    ### ~ PROCESSING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    orftypes=LIST   : List of alternative ORF types to generate (eORF/tORF/uORF/oORF/dORF) [e,t,u,o]
    cdsonly=T/F     : Whether to restrict analysis to cDNA with matching CDS [True]
    context=DICT    : Dictionary file of codon context *XXX* and AUG strength rating (Context:Strength) [see docs]
    altcodons=LIST  : List of acceptable non-AUG start codons (RNA or DNA) [CTG,GTG]
    altcontext=LIST : List of contexts to be considered for altcodons [XXXG]
    nrflanks=X      : Flanking sequence length (added 5' & 3') for analysing for redundancy within a gene [10]
    fullnr=T/F      : Perfrom NR Flank analysis across all genes [False]
    minorf=X        : Minimum ORF lengths to be considered [0]
    minutr=X        : Minimum 5' UTR lengths to be included in analysis [1]
    mincds=X        : Minimum CDS lengths to be included in analysis [0]

    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : Root name for output files. Path will be stripped and resdir used. [* based on cdna or singleseq]
    resdir=PATH     : Path for results output [./]

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_sequence, rje_slim
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.1.0 - Initial Compilation based on PATIS V0.3.
    # 0.2.0 - Added tabular and fasta output. Basic REST service output but not tested.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Add functionality for single sequences.
    # [ ] : Add processing of aic table into full (filtered) ExTATIC table.
    # [Y] : Add eORF, uORF and dORF sequence output. Mark internal AIC with case until Strong then all UC.
    # [Y] : Sort out REST output for single sequences only. Add context and altcodons etc. to REST outputs.
    # [Y] : Replace Context File with dictionary.
    # [?] : Do we want a minimum extension filter? Post-filtering option? (Add AltLen field.)
    # [?] : Replace &rest=X format for seqin/cdna/cds to return fasta generated by seqobj
    # [Y] : Add human genes to an alias table for REST servers.
    # [Y] : Add nrflanks=X for analysing for redundancy within a gene. (fullnr=T/F switch to count all as same gene?)
    # [Y] : Add QC step before main ExTATIC analysis: dump transcripts with dodgy CDS.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('ExTATIC', '0.2.0', 'January 2015', '2014')
    description = 'Extensions and Truncations from Alternative Translation Initiation Codons'
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
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
default_contexts = {'XXXG':'Strong','RnnXXXH':'Mid','CnnXXXU':'Mid','UnnXXXC':'Mid','CnnXXXM':'Weak','UnnXXXW':'Weak',
                    '^XXXH':'Unknown','^nXXXH':'Unknown','^nnXXXH':'Unknown'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class ExTATIC(rje_obj.RJE_Object):
    '''
    ExTATIC Class. Author: Rich Edwards (2014).

    Str:str
    - CDSIn=FILE      : CDS input sequence file (*.cds.fas). See docs for naming advice. [ExTATIC.cds.fas]
    - ResDir=PATH     : Path for results output [./]
    - SeqIn=FILE      : cDNA input seqence file (*.cdna.fas). See docs for naming advice. [ExTATIC.cdna.fas]
    - SingleSeq=X     : Analyse a single sequence only. SeqIn and CDSIn should be sequence strings. [None]

    Bool:boolean
    - CDSOnly=T/F     : Whether to restrict analysis to cDNA with matching CDS [True]
    - FullNR=T/F      : Perfrom NR Flank analysis across all genes [False]

    Int:integer
    - MinCDS=X        : Minimum CDS lengths to be included in analysis [0]
    - MinORF=X        : Minimum ORF lengths to be considered [0]
    - MinUTR=X        : Minimum 5' UTR lengths to be included in analysis [1]
    - NRFlanks=X      : Flanking sequence length (added 5' & 3') for analysing for redundancy within a gene [10]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - AltCodons=LIST  : List of acceptable non-AUG start codons (RNA or DNA) [CTG,GTG]
    - AltContext=LIST : List of context to be considered for altcodons [XXXG]
    - ORFTypes=LIST   : List of alternative ORF types to generate (eORF/tORF/uORF/oORF/dORF) [e,t,u,o]

    Dict:dictionary    
    - Context=DICT    : Dictionary file of codon context *XXX* and AUG strength rating (Context,Strength) [see docs]

    Obj:RJE_Objects
    - cDNA = SeqList object of cDNA sequences.
    - CDS = SeqList object of CDS sequences.
    - DB = Database object of key data. (Populated during analysis.)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['CDSIn','ResDir','SeqIn','SingleSeq']
        self.boollist = ['CDSOnly','FullNR']
        self.intlist = ['MinCDS','MinORF','MinUTR','NRFlanks']
        self.numlist = []
        self.filelist = []
        self.listlist = ['AltCodons','AltContext','ORFTypes']
        self.dictlist = ['Context']
        self.objlist = ['DB','CDS','cDNA']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'SeqIn':'ExTATIC.cdna.fas','ResDir':rje.makePath('./')})
        self.setBool({'CDSOnly':True})
        self.setInt({'NRFlanks':10,'MinUTR':1})
        self.setNum({})
        self.list['AltCodons'] = ['CTG','GTG']
        self.list['AltContext'] = ['XXXG']
        self.list['ORFTypes'] = ['e','t','u','o']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
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
                self._cmdRead(cmd,type='file',att='CDSIn',arg='cds')
                self._cmdRead(cmd,type='file',att='SeqIn',arg='cdna')
                self._cmdReadList(cmd,'str',['SingleSeq'])   # Normal strings
                self._cmdReadList(cmd,'path',['ResDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['SeqIn','CDSIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['CDSOnly','FullNR'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinCDS','MinORF','MinUTR','NRFlanks'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['AltCodons','AltContext','ORFTypes'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['Context']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [1] Process Commands where required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        orftypes = []
        for orf in self.list['ORFTypes']:
            if orf.lower()[:1] in 'etuod': orftypes.append(orf.lower()[:1])
            else: self.warnLog('ORFType "%s" not recognised!' % orf)
        self.list['ORFTypes'] = orftypes
        for orf in ['Extended','Truncated','Upstream','Overlapping','Downstream']:
            self.printLog('#%sORF' % orf[0],'Analyse %s ORFs (%sORFs): %s' % (orf,orf[0].lower(),orf[0].lower() in orftypes))
        if not orftypes: self.warnLog('No alternative initiation ORF types selected for analysis!')
        if self.getInt('MinUTR') < 1: self.warnLog('Switching off 5\' UTR requirement (minutr=0) not recommended.')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.extatic(): return False
            self.orfFasta()
            #!# Sort, document and check REST output & add table fields to docstring.
            #!# Add identification of redundant and internal AIC
            #!# Add SignalP, TargetP and Domain Prediction
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [0] Setup Method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Rest'): self.restSetup()
            ### ~ [1] Setup File Names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = self.setupFileNames()            # Will return basefile from sequence files, or None if problem.
            ### ~ [2] Setup Database Object/Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupCodons(): return False     # Will return True/False
            ### ~ [3] Load Sequence Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupSequences(): return False  # Will return number of transcripts in db Table
            ### ~ [4] Summarise and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SETUP','ExTATIC %s setup complete.' % basefile)
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def setupFileNames(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup File Names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            need_cds = self.getBool('CDSOnly')
            if self.getStrLC('SingleSeq'):
                basefile = rje.split(self.getStr('SingleSeq'))[0]
                if not self.getStrLC('SeqIn'): raise ValueError('No SeqIn cDNA sequence given!')
                if need_cds and not self.getStrLC('CDSIn'): raise ValueError('No CDSIn CDS sequence given (cdsonly=T)!')
            else:
                if not rje.exists(self.getStr('SeqIn')): raise IOError('SeqIn "%s" not found!' % self.getStr('SeqIn'))
                self.printLog('#CDNA','cDNA File "%s" found.' % self.getStr('SeqIn'))
                basefile = rje.baseFile(self.getStr('SeqIn'))
                if basefile.endswith('.cds'): raise ValueError('CDS (*.cds.fas) given for SeqIn rather than cDNA!')
                #if not basefile.endswith('.cdna'): basefile += '.cdna'
                if basefile.endswith('.cdna'): basefile = rje.baseFile(basefile)
                if not rje.exists(self.getStr('CDSIn')):
                    #cdsfile = '%s.cds.fas' % rje.baseFile(basefile)
                    cdsfile = '%s.cds.fas' % basefile
                    if rje.exists(cdsfile): self.setStr({'CDSIn':cdsfile})
                    elif need_cds: raise IOError('CDSIn "%s" or "%s" not found!' % (self.getStr('CDSIn'),cdsfile))
                self.printLog('#CDS','CDS File "%s" found.' % self.getStr('CDSIn'))
            if self.getStrLC('Basefile'): basefile = self.basefile(runpath=True)
            self.basefile('%s%s' % (self.getStr('ResDir'),os.path.basename(basefile)))
            self.printLog('#BASE','Output basefile set: %s' % self.basefile(runpath=True))
            return basefile
        except: self.errorLog('Problem during %s.setupFileNames()' % self.prog()); return None  # Setup failed
#########################################################################################################################
    def setupCodons(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup Database Object/Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            self.db().addEmptyTable('aic',['ENST','Pos','RF','Codon','Context','Type','Strength','Stop'],['ENST','Pos'])
            ## ~ [1a] Setup Context Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb = db.addEmptyTable('context',['Context','Strength'],['Context','Strength'])
            if self.dict['Context']:
                for context in self.dict['Context']: cdb.addEntry({'Context':context,'Strength':self.dict['Context'][context]})
            else:
                for context in default_contexts:  cdb.addEntry({'Context':context,'Strength':default_contexts[context]})
            for context in self.list['AltContext']: cdb.addEntry({'Context':context,'Strength':'Alt'})
            #?# Extend contexts where necessary to make same length +/- AUG position. (Why?)
            #!# Add code to check integrity of contexts: cannot have conflicts
            #!# Add contexts to rest output. This could be useful w/o sequences for checking input
            conlen = [0,0]    # Length of [5',3'] context
            for ckey in cdb.dataKeys():
                entry = cdb.data(ckey)
                consplit = rje.split(entry['Context'].upper(),'XXX')
                if len(consplit) != 2: raise ValueError('Context "%s" format error!' %  entry['Context'])
                consplit[0] = rje_slim.patternFromCode(rje_slim.slimFromPattern(consplit[0],dna=True),dna=True)
                conlen[0] = max(conlen[0],rje_slim.slimLen(consplit[0]))
                consplit[1] = rje_slim.patternFromCode(rje_slim.slimFromPattern(consplit[1],dna=True),dna=True)
                conlen[1] = max(conlen[1],rje_slim.slimLen(consplit[1]))
                self.printLog('#START','%s = %s.' % (entry['Context'],entry['Strength']))
            self.printLog('#START','Start site context ranges from -%d to +%d' % (conlen[0],conlen[1]))
            #self.printLog('#AIC','Alternative codons %s allowed in %s context.' % (rje.join(self.list['AltCodons'],'/'),rje.join(self.list['AltContext'],'/')))
            self.dict['Output']['altcodons'] = rje.join(self.list['AltCodons'],'\n')
            #?# Add check of integrity of cdb, i.e. non-overlapping contexts? Do during analysis?
            return True
        except: self.errorLog('Problem during %s.setupCodons().' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def setupSequences(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Load Sequence Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if self.getStrLC('SingleSeq'):
                self.dict['Output']['name'] = self.getStr('SingleSeq')
                edb = None
                scmd = ['autoload=F','seqmode=tuple']
                self.obj['cDNA'] = rje_seqlist.SeqList(self.log,self.cmd_list+scmd)
                self.obj['cDNA'].list['Seq'] = [(self.getStr('SingleSeq'),self.getStr('SeqIn'))]
                self.dict['Output']['cdna'] = '>%s\n%s\n' % (self.getStr('SingleSeq'),self.getStr('SeqIn'))
                if self.getStrLC('CDSIn'):
                    self.obj['CDS'] = rje_seqlist.SeqList(self.log,self.cmd_list+scmd)
                    self.obj['CDS'].list['Seq'] = [(self.getStr('SingleSeq'),self.getStr('CDSIn'))]
                self.dict['Output']['cds'] = '>%s\n%s\n' % (self.getStr('SingleSeq'),self.getStr('CDSIn'))
            else:
                scmd = ['autoload=T','seqmode=file']
                self.obj['cDNA'] = rje_seqlist.SeqList(self.log,self.cmd_list+scmd+['seqin=%s' % self.getStr('SeqIn')])
                #self.dict['Output']['cdna'] = '%s sequences loaded from %s' % (self.obj['cDNA'].seqNum(),self.getStr('SeqIn'))
                self.dict['Output']['cdna'] = self.getStr('SeqIn')
                if rje.exists(self.getStr('CDSIn')):
                    self.obj['CDS'] = rje_seqlist.SeqList(self.log,self.cmd_list+scmd+['seqin=%s' % self.getStr('CDSIn')])
                    #self.dict['Output']['cds'] = '%s sequences loaded from %s' % (self.obj['CDS'].seqNum(),self.getStr('CDSIn'))
                    self.dict['Output']['cds'] = self.getStr('CDSIn')
                edb = self.db('transcripts')
            for otype in 'etod':
                if not self.obj['CDS'] or not self.obj['CDS'].seqNum() and otype in self.list['ORFTypes']:
                    self.warnLog('Cannot identify %sORF without CDS data.' % otype)
            ### ~ [2] Generate Sequence Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if edb: edb.dataFormat({'cDNA':'int','CDS':'int','TSS':'int','Stop':'int'})
            else:
                edb = db.addEmptyTable('transcripts',['ENST','ENSG','Gene','Description','cDNA','CDS','Length'],['ENST'])
                for stype in ['cDNA','CDS']:
                    sobj = self.obj[stype]
                    if not sobj: continue
                    sx = 0; px = 0.0; ptot = sobj.seqNum()
                    while sobj.nextSeq():
                        self.progLog('\r#SEQ','Processing transcript %s: %.1f%%' % (stype,px/ptot)); px += 100.0
                        (name,sequence) = sobj.currSeq()
                        if 'Sequence unavailable' in sequence: continue
                        namedata = rje.split(name,'|',3)
                        if len(namedata) < 3: enst = ensg = gene = namedata[0]
                        else: [enst,ensg,gene] = namedata[:3]
                        if len(namedata) > 3: desc = namedata[3]
                        else: desc = gene
                        #if stype == 'cDNA' and enst in edb.dict['Data']:
                        #    self.bugPrint(name); self.bugPrint(sequence)
                        #    self.debug(edb.data(enst))
                        sx += 1
                        if stype == 'cDNA': edb.addEntry({'cDNA':sobj.obj['Current'],'ENST':enst,'ENSG':ensg,'Gene':gene,'Description':desc,'CDS':-1,'Length':len(sequence)})
                        else: edb.data(enst)['CDS'] = sobj.obj['Current']
                        if self.getStrLC('SingleSeq'): edb.data(enst)[stype] = edb.data(enst)[stype][1]
                    self.printLog('\r#SEQ','%s of %s transcripts with %s sequence data.' % (rje.iStr(sx),rje.iStr(sobj.seqNum()),stype))
                #if not self.getStrLC('SingleSeq'): edb.saveToFile()
            if self.getBool('CDSOnly'): edb.dropEntriesDirect('CDS',[-1])
            return edb.entryNum()     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return 0  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        context = Table of codon contexts and strength ratings. [tdt]
        altcodons = List of alternative (non-AUG) codons hopefully.
        orfs = output table of predicted open reading frames. [tdt]
        aic = output table of predicted alternative initiation sites. [tdt]
        nr = output of non-redundant AIC, grouped by codons plus nrflanks=X flanking sequence. (By gene unless fullnr=T) [tdt]
        transcripts = output table summarising ORFs and AIC for each transcript. [tdt]
        fas = ORF fasta sequences with upper case AIC and lower case sequence until the first Strong ATG codon. [fas]
        name = Sequence name (&singleseq=X).
        cdna = Input cDNA sequence(s). [fas]
        cds = [Optional] Input coding sequence(s). [fas]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Output formats that can directly fetch data from self.str ~~~~~~~~~~~~~~~~~ ##
            self.dict['Output'] = {'cds':'CDSIn','cdna':'SeqIn','name':'SingleSeq'}
            ## ~ [0b] ~ These outputs need to be set during data processing ~~~~~~~~~~~~~~~~~~~~~~~ ##
            for outfmt in ['fas','altcodons']: self.dict['Output'][outfmt] = 'No output generated.'
            ## ~ [0c] ~ Output formats that can directly fetch data from self.str ~~~~~~~~~~~~~~~~~ ##
            # aic, orfs, transcripts, nr, context

            # The following database tables do not need to be in self.dict['Output'] and should get processed anyway:
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            #?# Only have singleseq=T output #?#
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['context','altcodons','orfs','aic','nr','transcripts','fas','name','cdna','cds']
#########################################################################################################################
    def single(self): return self.getStrLC('SingleSeq')
#########################################################################################################################
    ### <3> ### AIC Identification Methods                                                                              #
#########################################################################################################################
    def extatic(self):  ### Main ExTATIC processing method. Send each sequence to extaticSeq.
        '''Main ExTATIC processing method. Send each sequence to extaticSeq.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tdb = self.db('transcripts')
            if not self.force():
                adb = self.db().addTable(name='aic',mainkeys=['ENST','Pos'],expect=False)
                if adb: adb.dataFormat({'Pos':'int','RF':'int','ORFPos':'int','ORFLen':'int','NRX':'int','Stop':'int'})
                odb = self.db().addTable(name='orfs',mainkeys=['ENST','ORF'],expect=False)
                if odb: odb.dataFormat({'Start':'int','RF':'int','Stop':'int','aORF':'int','eORF':'int','oORF':'int',
                                        'tORF':'int','uORF':'int','dORF':'int','ORFLen':'int'})
                ndb = self.db().addTable(name='nr',mainkeys=['NRFlanks'],expect=False)
                if tdb and adb and odb and ndb: return True
            cds = self.obj['CDS']
            cdna = self.obj['cDNA']
            needtss = 'TSS' not in tdb.fields()
            if needtss: tdb.addField('TSS',evalue=-1); tdb.addField('Stop',evalue=-1)
            if 'QC' not in tdb.fields(): tdb.addField('QC',evalue='OK')
            #!# Should not need this but do! Somewhere, output is getting populated
            #for table in ['aic','transcripts','orfs']: self.dict['Output'][table] = '%s.%s.tdt' % (self.basefile(runpath=True),table)
            ### ~ [1] Process Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rejx = tdb.entryNum()
            badcdsx = 0
            nocdstopx = 0
            frameshiftx = 0
            nocdsx = 0
            mincdsx = 0
            minutrx = 0
            self.printLog('#AIC','Generating AIC for %s transcripts...' % rje.iStr(rejx),log=False)
            for entry in tdb.entries():
                cdsseq = ''
                if self.single(): (name,cdnaseq) = (entry['ENST'],entry['cDNA'])
                else: (name,cdnaseq) = cdna.getSeq(entry['cDNA'])
                entry['Stop'] = len(cdnaseq)
                if entry['CDS'] >= 0:
                    if self.single(): (name,cdsseq) = (entry['ENST'],entry['CDS'])
                    else: (name,cdsseq) = cds.getSeq(entry['CDS'])
                    if cdsseq not in cdnaseq: entry['QC'] = 'Rejected: CDS not found in cDNA.'; badcdsx +=1 ; continue
                    if cdsseq[-3:] not in ['TGA','TAA','TAG']: entry['QC'] = 'Rejected: CDS does not end with STOP codon (%s).' % cdsseq[-3:]; nocdstopx +=1 ; continue
                    if len(cdsseq) % 3: entry['QC'] = 'Rejected: CDS sequence %s...%s has frameshift (Len=%dnt)' % (cdsseq[:3],cdsseq[-3:],len(cdsseq)); frameshiftx +=1 ; continue
                    if needtss: entry['TSS'] = cdnaseq.find(cdsseq) + 1
                elif self.getBool('CDSOnly'): entry['QC'] = 'Rejected: No CDS.'; nocdsx +=1 ; continue
                if len(cdsseq) < self.getInt('MinCDS'):
                    entry['QC'] =  'Rejected: CDS (%d nt) < MinCDS=%d' % (len(cdsseq),self.getInt('MinCDS'))
                    if self.getStrLC('SingleSeq'): self.dict['Output']['transcripts'] = 'Rejected: CDS (%d nt) < MinCDS=%d' % (len(cdsseq),self.getInt('MinCDS'))
                    mincdsx +=1 ; continue
                if entry['TSS'] <= self.getInt('MinUTR') > 0:
                    entry['QC'] = 'Rejected: 5\' UTR (%d nt) < MinUTR=%d' % (entry['TSS']-1,self.getInt('MinUTR'))
                    if self.getStrLC('SingleSeq'): self.dict['Output']['transcripts'] = 'Rejected: 5\' UTR (%d nt) < MinUTR=%d' % (entry['TSS']-1,self.getInt('MinUTR'))
                    minutrx +=1 ;continue
                self.extaticSeq(cdnaseq,cdsseq,entry['ENST']); rejx -= 1
            if badcdsx: self.printLog('#REJ','%s transcripts rejected for bad CDS (not found in cDNA).' % (rje.iStr(badcdsx)))
            if nocdstopx: self.printLog('#REJ','%s transcripts rejected for bad CDS (no STOP codon).' % (rje.iStr(nocdstopx)))
            if frameshiftx: self.printLog('#REJ','%s transcripts rejected for bad CDS (frameshift).' % (rje.iStr(frameshiftx)))
            if nocdsx: self.printLog('#REJ','%s transcripts rejected for no CDS (cdsonly=T).' % (rje.iStr(nocdsx)))
            if mincdsx: self.printLog('#REJ','%s transcripts rejected for short CDS (mincds=%d).' % (rje.iStr(mincdsx),self.getInt('MinCDS')))
            if minutrx: self.printLog('#REJ','%s transcripts rejected for short 5\' UTR (minutr=%d).' % (rje.iStr(minutrx),self.getInt('MinUTR')))
            self.printLog('#REJ','%s of %s transcripts rejected for QC reasons.' % (rje.iStr(rejx),rje.iStr(tdb.entryNum())))
            ### ~ [2] Compile ORF Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb = self.db('aic')
            odb = self.db().copyTable(adb,'orfs')
            odb.dropFields(['Codon','Context','Strength'])
            ## ~ [2a] Use STOP positions to identify AIC in same ORF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            odb.compress(['ENST','Stop','Type'],rules={'RF':'min','Pos':'list'},joinchar='|')
            odb.reshapeWide('Type',reshape=['Pos'],evalue='')   # Keys now ENST and Stop
            for otype in ['a']+self.list['ORFTypes']:
                if 'Pos|%sORF' % otype in odb.fields(): odb.renameField('Pos|%sORF' % otype,'%sORF' % otype)
                else: odb.addField('%sORF' % otype,evalue='')
            ## ~ [2b] Identify START and convert AIC positions to ORF positions ~~~~~~~~~~~~~~~~~~~ ##
            odb.addField('Pos')     # List of all AIC AA positions in cDNA
            odb.addField('ORFPos')  # List of all AIC AA positions in ORF
            odb.addField('Start')   # Start position in cDNA of 5' codon
            odb.addField('Type')    # Type of 5' codon
            odb.addField('ORFLen')  # Length of ORF
            for entry in odb.entries():
                entry['Pos'] = []
                for otype in 'tdaoeu':  # Will go in order and over-write to most 5'
                    if otype not in ['a']+self.list['ORFTypes']: continue
                    ox = 0
                    for pos in rje.split('%s' % entry['%sORF' % otype],'|'):
                        try: entry['Pos'].append(int(pos)); entry['Type'] = '%sORF' % otype; ox += 1
                        except: pass
                    entry['%sORF' % otype] = ox
                entry['Pos'].sort()
                entry['ORFPos'] = entry['Pos'][0:]
                entry['Start'] = entry['Pos'][0]
                # Convert Pos to AA pos
                entry['ORFLen'] = (entry['Stop'] - entry['Start']) // 3
                for i in range(len(entry['ORFPos'])): entry['ORFPos'][i] = 1 + (entry['Pos'][i] - entry['Start']) // 3
                entry['Pos'] = rje.replace('%s' % entry['Pos'],', ','|')[1:-1]
                entry['ORFPos'] = rje.replace('%s' % entry['ORFPos'],', ','|')[1:-1]
            ## ~ [2c] Convert to use Start as Key and create ORF identifiers ~~~~~~~~~~~~~~~~~~~~~~ ##
            odb.newKey(['ENST','Type','Start'],startfields=True)    # Include Type for sorting
            odb.addField('ORF')     # ORF Identifier (numbered within Type in position order)
            enst = None; otype = None; orfx = 0
            for okey in odb.dataKeys():
                entry = odb.data(okey)
                if entry['ENST'] != enst or entry['Type'] != otype: orfx = 0
                enst = entry['ENST']; otype = entry['Type']; orfx += 1
                entry['ORF'] = '%s%d' % (otype,orfx)
            odb.newKey(['ENST','ORF'],startfields=True)

            ### ~ [3] Update AIC table with ORF Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb.addField('ORF'); adb.addField('ORFPos'); adb.addField('ORFLen')
            for entry in odb.entries():
                orf = entry['ORF']
                orfpos = rje.split(entry['ORFPos'],'|')
                dnapos = rje.split(entry['Pos'],'|')
                for i in range(len(dnapos)):
                    akey = adb.makeKey({'ENST':entry['ENST'],'Pos':dnapos[i]})
                    try:
                        adb.data(akey)['ORF'] = orf
                        adb.data(akey)['ORFPos'] = orfpos[i]
                        adb.data(akey)['ORFLen'] = entry['ORFLen'] - int(orfpos[i]) + 1
                    except: self.errorLog('ORF/AIC Mismatch')

            ### ~ [4] Generate NR AIC by gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db().copyTable(adb,'nr')
            gdb.index('ENST')
            gdb.makeField('#ENST#:#ORF#:#Type#:#ORFPos#','AIC')
            gdb.makeField('#AIC#','NRFlanks')         # Use AIC for redundancy unless NRFlanks given
            gdb.dropFields(['Type','Stop','ORFLen','RF'])
            # Add genes
            gdb.addField('Gene',after='ENST')
            for entry in gdb.entries(): entry['Gene'] = tdb.data(enst)['Gene']
            # Add flanking sequences
            if self.getInt('NRFlanks') > -1:
                sobj = self.obj['cDNA']
                sobj.obj['Current'] = None
                px = 0.0; ptot = sobj.seqNum(); ox = 0
                while sobj.nextSeq():
                    self.progLog('\r#NR','Generating nrflanks: %.1f%%' % (px/ptot)); px += 100.0
                    (name,sequence) = sobj.currSeq()
                    if 'Sequence unavailable' in sequence: continue
                    enst = rje.split(name,'|',3)[0]
                    for entry in gdb.indexEntries('ENST',enst):
                        entry['NRFlanks'] = sequence[max(0,entry['Pos']-self.getInt('NRFlanks')-1):entry['Pos']+self.getInt('NRFlanks')+3]
                self.printLog('\r#NR','Generated nrflanks +/- %d nt' % self.getInt('NRFlanks'))
            #gdb.saveToFile(); self.deBug('...')
            # Compress to unique AIC
            gdb.newKey(['AIC'])
            gdb.dropFields(['ENST','Pos','ORF','ORFPos'])
            # >> ENST Gene Codon Context AIC NRFlanks
            if self.getBool('FullNR'): gdb.compress(['NRFlanks'],default='list',joinchar='|')
            else: gdb.compress(['Gene','NRFlanks'],default='list',joinchar='|')
            gdb.list['Fields'] = ['NRFlanks','Codon','Context','Strength','Gene','AIC']
            ## ~ [4a] Update AIC table with NR data? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb.index('AIC',splitchar='|')
            adb.addField('NR'); adb.addField('NRX')
            for entry in adb.entries():
                nkey = rje.join([entry['ENST'],entry['ORF'],entry['Type'],entry['ORFPos']],':')
                try: nentry = gdb.indexEntries('AIC',nkey)[0]
                except: self.bugPrint(nkey); continue
                aiclist = rje.split(nentry['AIC'],'|')
                entry['NRX'] = len(aiclist)     # This might be enough!
                # Rate as NR, Redundant, Internal ...?
                if len(aiclist) == 1: entry['NR'] = 'NR'            # NR = Unique
                elif entry['ORFPos'] > 1: entry['NR'] = 'Redundant' # Redundant = 2+ Start or Internal
                else:
                    entry['NR'] = 'Redundant'                       # ORF internal that is same as a start is also "Redundant"
                    for aic in aiclist:                             # Internal = ORF start that is inside another ORF
                        if int(rje.split(aic,':')[-1]) > 1: entry['NR'] = 'Internal'
                #X#self.debug(entry)

            #?# Add AIC number to AIC output, e.g. 1 to N in 5' to 3' order: get from odb. key = 'ENST','ORF'

            ### ~ [5] Save Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Rest'): tdb.dropFields(['cDNA','CDS'])
            if needtss: tdb.saveToFile(append=self.getBool('Append'))
            else: tdb.saveToFile(append=False)
            adb.saveToFile(append=self.getBool('Append'))
            odb.saveToFile(append=self.getBool('Append'))
            gdb.saveToFile(append=self.getBool('Append'))
            return True
        except: self.errorLog('%s.extatic() error' % self.prog()); return False
#########################################################################################################################
    def extaticSeq(self,cdnaseq,cdsseq='',enst=''):  ### Main ExTATIC processing method for one sequence.
        '''
        Main ExTATIC processing method for one sequence. Positions added to table are 1 to L.
        >> cdnaseq:str = cDNA sequence.
        >> cdsseq:str [''] = Optional CDS sequence. (uORF until first strong only if not given.)
        >> enst:str [''] = Optional sequence identifier.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NOTE: Will need to add additional strategy options regarding which AIC to return and when to stop looking.
            # Initial rule = first Strong downstream in-frame ATG in CDS if given, or first strong ATG if no CDS
            # Might want to add additional criteria for overlapping uORFs etc.
            ex = 0  # Number of entries being added
            odb = self.db('aic')    # ['ENST','Pos','RF','Codon','Context','Type','Strength','Stop'], keys:['ENST','Pos']
            cdb = self.db('context')   # ['Context','Strength']
            edb = self.db('transcripts')
            cdnaseq = cdnaseq.upper()
            cdsseq = cdsseq.upper()
            cdstop = len(cdnaseq)
            ## ~ [0a] Find stop codons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            stops = self.codonPositions('(TGA|TAA|TAG)',cdnaseq)    #?# Not sure if this is required
            stops = rje.sortKeys(stops)
            ## ~ [ba] TSS from CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tss = -1
            if cdsseq:
                tss = cdnaseq.find(cdsseq)
                if tss < 0: self.warnLog('%s CDS not found in cDNA! Rejected.' % enst); return 0
            if cdsseq and stops:
                cdstop = stops[0:]
                while cdstop and (cdstop[0] < tss or (cdstop[0] - tss) % 3): cdstop.pop(0)
                if cdstop: cdstop = cdstop[0]
                else:
                    self.bugPrint('%s:%s' % (tss,cdnaseq[tss:][:3]))
                    self.bugPrint(stops)
                    self.deBug('%d:%d' % (len(cdsseq),len(cdsseq)%3))
                    self.warnLog('%s CDS does not have in-frame (TGA|TAA|TAG) STOP codon! Rejected.' % enst); return 0
            elif cdsseq:
                self.warnLog('%s CDS does not have (TGA|TAA|TAG) STOP codon! Rejected.' % enst)
                return 0
            edb.data(enst)['Stop'] = cdstop + 1     # Will be > cDNA Length if no stop codon!
            ### ~ [1] Process Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            codons = {}     # Dictionary of {Strength:{Codon position:match}}
            ## ~ [1a] Establish "end" point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Identify the "end" point = First strong ATG > TSS, in-frame if CDS.
            codon = 'ATG'
            for centry in cdb.entries():
                context = centry['Context']
                strength = centry['Strength']
                if strength == 'Alt': continue
                if strength not in codons: codons[strength] = {}
                rje.combineDict(codons[strength],self.codonPositions(codon,cdnaseq,context),overwrite=False)
            try: end = rje.sortKeys(codons['Strong'])[0]     # If not TSS, Endpoint will be first strong ATG
            except: end = len(cdnaseq)
            if tss > -1:    # End point is first in-frame Strong ATG codon 3' of TSS (including TSS itself)
                end = cdstop
                for cpos in rje.sortKeys(codons['Strong']):
                    if cpos > end: break            # Overshot stop codon
                    if (tss - cpos) % 3: continue   # Wrong frame
                    if cpos >= tss: end = cpos; break   # Found an in-frame 3' strong ATG
            ## ~ [1c] Add alternative codon positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for codon in self.list['AltCodons']:    # Could combine with: '(%s)' % rje.join(self.list['AltCodons'])
                if codon not in codons: codons[codon] = {} # If combined, extract actual codon using position and remake dictionary
                for context in self.list['AltContext']:
                    rje.combineDict(codons[codon],self.codonPositions(codon,cdnaseq[:end],context),overwrite=False)
            ### ~ [2] Update table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cstrengths = rje.sortKeys(codons)
            cstrengths.remove('Strong')
            cstrengths.insert(0,'Strong')
            stopstrong = {}     # Dictionary of 5' Strong codon positions for each ORF (stop position)
            tssfound = tss < 0
            #self.debug(codons)
            for strength in cstrengths:
                #self.bugPrint('%s=%s' % (strength,rje.sortKeys(codons[strength])))
                for cpos in rje.sortKeys(codons[strength]):
                    if cpos > end: break
                    ## ~ [2a] Establish position of first in-frame stop codon ~~~~~~~~~~~~~~~~~~~~~ ##
                    # Can be used to identify codons in same ORF as well as establish ORF length
                    stop = stops[0:]
                    while stop and (stop[0] < cpos or (stop[0] - cpos) % 3): stop.pop(0)
                    if stop: stop = stop[0]
                    else:
                        stop = len(cdnaseq)
                        while (stop-cpos) % 3: stop += 1
                    if stop == cdstop and cpos < tss: pass      # Do not stop annotated ORF until TSS reached.
                    elif stop not in stopstrong and strength == 'Strong': stopstrong[stop] = cpos
                    elif stop in stopstrong and cpos > stopstrong[stop]: continue   # 3' of strong codon in same ORF
                    ## ~ [2b] Establish reading frame relative to CDS. Will be -1 if no CDS ~~~~~~~ ##
                    if tss < 0: rf = -1
                    else: rf = (tss-cpos) % 3
                    ## ~ [2c] Establish type of codon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if tss < 0 or cpos < tss: ctype = 'uORF'                            # Upstream ORF
                    elif cpos == tss: ctype = 'aORF'; tssfound = True                   # Annotated ORF
                    elif rf: ctype = 'dORF'                                             # Downstream ORF (out of frame)
                    elif stop == cdstop: ctype = 'tORF'                                 # Truncated ORF
                    else:
                        self.bugPrint(cdnaseq)
                        self.bugPrint(cdsseq)
                        self.warnLog('Impossible %s ORF detected: %d (%s) RF%s; TSS=%s; Stop=%d; CDStop=%d. Something has gone wrong!' % (enst,cpos,codons[strength][cpos],rf,tss,stop,cdstop))
                        self.debug('%s...%s' % (cdsseq[:3],cdsseq[-3:]))
                        self.bugPrint(cdnaseq[:tss])
                        self.debug('%s%s' % (cdnaseq[tss-3:tss].lower(),cdnaseq[tss:tss+4].upper()))
                        continue
                    if ctype == 'uORF' and stop > 0 and stop == cdstop: ctype = 'eORF'  # Extended ORF
                    elif ctype == 'uORF' and stop > tss and cpos < tss: ctype = 'oORF'  # Overlapping ORF
                    if ctype[:1] not in self.list['ORFTypes']+['a']: continue
                    if (stop - cpos) // 3 < self.getInt('MinORF'): continue
                    ## ~ [2d] Update table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    entry = {'ENST':enst,'Pos':cpos+1,'RF':rf,'Codon':cdnaseq[cpos:][:3],'Type':ctype,'Strength':strength,'Stop':stop+1,'Context':codons[strength][cpos]}
                    #ekey = '%s\t%s' % (enst,cpos)
                    #self.bugPrint(entry)
                    #if odb.data(ekey): self.deBug(odb.data(ekey));
                    odb.addEntry(entry); ex += 1
            if not tssfound:
                entry = {'ENST':enst,'Pos':tss+1,'RF':0,'Codon':cdnaseq[tss:][:3],'Type':'aORF','Strength':cdnaseq[tss:][:3],'Stop':cdstop+1,'Context':cdnaseq[tss:][:3]}
                self.warnLog('Transcript %s has unusual annotated start codon (pos %d): %s%s%s' % (enst,tss+1,cdnaseq[max(0,tss-3):tss].lower(),entry['Codon'],cdnaseq[tss+3:tss+4].lower()))
                #self.bugPrint(cdnaseq)
                #self.bugPrint(cdsseq)
                #self.debug('%s...%s' % (cdsseq[:3],cdsseq[-3:]))
                #self.bugPrint(cdnaseq[:tss])
                #self.debug('%s%s%s' % (cdnaseq[max(0,tss-3):tss].lower(),cdnaseq[tss:tss+3].upper(),cdnaseq[tss+3:tss+4].lower()))
                odb.addEntry(entry); ex += 1
                try: self.db('transcripts').data(enst)['QC'] = 'Unusual annotated start codon (%s%s%s)' % (cdnaseq[max(0,tss-3):tss].lower(),entry['Codon'],cdnaseq[tss+3:tss+4].lower())
                except: pass
        except: self.errorLog('%s.extaticSeq error' % self.prog())
        return ex
#########################################################################################################################
    def codonPositions(self,codon,sequence,context='XXX',frames=[],startpos=0):    ### Returns a list of codon positions (first base of codon)
        '''
        Returns a list of codon positions (first base of codon) in sequence. To put a cap on the position, simply
        truncate the sequence. Variable-length contexts may not work.
        >> codon:str = 3-letter codon to look for, or (AAA|BBB|CCC) options.
        >> sequence:str = Sequence in which to search for codons.
        >> context:str ['XXX'] = Codon context. XXX will be replaced with codon. IUPAC codes will be used to match.
        >> frames:list [] = Limit to matching certain reading frames (0/1/2)
        >> startpos:int [0] = Position to start looking (0 to L-1).
        << codons:dict = Dictionary of {position:context}, positions numbers 0 to L-1.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            i = startpos                # Position in sequence, will increment as matches found
            sequence = sequence.upper() # Sequence to search through.
            codons = {}                 # Codon dictionary to return
            ## ~ [0a] Convert context to regular expression for matching ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            context_regex = []          # Context with IUPAC conversion and codon replacement
            for n in context.upper():
                try:
                    base = rje_slim.dna_ambig[n]
                    if len(base) > 1: context_regex.append('[%s]' % base)
                    else: context_regex.append(base)
                except: context_regex.append(n)
            if 'X' in context_regex: offset = context_regex.index('X')      # Position of codon within match
            context_regex = '(%s)' % rje.join(context_regex,'')
            if 'XXX' not in context_regex: raise ValueError('Cannot find XXX for codon in context %s -> %s' % (context,context_regex))
            context_regex = context_regex.replace('XXX',codon)
            if 'X' in context_regex: raise ValueError('Too many X in codon context %s -> %s' % (context,context_regex))
            ### [1] Perform codon search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while i < len(sequence):    # Continue until end of sequence reached. (Will actually break first!)
                match = rje.matchExp(context_regex,sequence[i:])
                if not match: break     # No more hits
                if '^' in context_regex: i = -1               # Offset will be inflated by ^
                else: i = sequence[i:].find(match[0]) + i     # Adjust position to actual start of sequence
                if frames and i % 3 not in frames: continue   # Invalid reading frame
                codons[i+offset] = match[0]
                if '^' in context_regex: break
                i += 1
            ### [2] Return codon matches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return codons
        except: self.errorLog('%s.codonPositions() error' % self.prog())
#########################################################################################################################
    ### <4> ### ORF Sequence Methods                                                                                    #
#########################################################################################################################
    def orfFasta(self): ### Generates ORF Fasta files
        '''Generates ORF Fasta files.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            orffas = '%s.orf.fas' % self.basefile(runpath=True)
            rje.backup(self,orffas)
            self.dict['Output']['fas'] = orffas
            odb = self.db('orfs'); odb.index('ENST')
            sobj = self.obj['cDNA']
            sobj.obj['Current'] = None
            px = 0.0; ptot = sobj.seqNum(); ox = 0
            FAS = open(orffas,'w')
            ### ~ [1] Generate ORF Sequences to fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while sobj.nextSeq():
                self.progLog('\r#FAS','Generating ORFs: %.1f%%' % (px/ptot)); px += 100.0
                (name,sequence) = sobj.currSeq()
                if 'Sequence unavailable' in sequence: continue
                enst = rje.split(name,'|',3)[0]
                if enst not in odb.index('ENST'): continue
                for okey in odb.index('ENST')[enst]:
                    entry = odb.data(okey)
                    name = '%s.%s AIC:%s ORFLen=%s' % (enst,entry['ORF'],entry['ORFPos'],entry['ORFLen'])
                    try:
                        orfseq = rje_sequence.dna2prot(sequence[entry['Start']-1:entry['Stop']+3]).lower()
                        i = -1
                        for pos in rje.split(entry['ORFPos'],'|'):
                            i = int(pos)-1
                            try: orfseq = orfseq[:i] + orfseq[i].upper() + orfseq[i+1:]
                            except:
                                self.bugPrint(orfseq)
                                self.bugPrint('SeqLen=%d' % len(orfseq))
                                self.debug(entry)
                                raise
                        if i >=0: orfseq = orfseq[:i] + orfseq[i:].upper()
                        FAS.write('>%s\n%s\n' % (name,orfseq)); ox += 1
                    except: self.warnLog('ORFFas output failed for: %s' % name,'orffas_fail',suppress=True)
            FAS.close()
            self.printLog('\r#FAS','Generated ORFs output for %s/%s ORFs' % (rje.iStr(ox),rje.iStr(odb.entryNum())))
        except: self.errorLog('%s.codonPositions() error' % self.prog())
#########################################################################################################################
### End of SECTION II: ExTATIC Class                                                                                    #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: ExTATIC(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
