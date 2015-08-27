#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <software@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      FIESTA
Description:  Fasta Input EST Analysis
Version:      1.9.0
Last Edit:    26/11/14
Citation:     Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: 20924652]
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    FIESTA has three primary functions: <1>	Discovery, assembly and evolutionary analysis of candidate genes in an EST
    library; <2> Assembly of an EST library for proteomics analysis; <3> Translation/Annotation of an EST library for
    proteomics analysis. These functions are outlined below.

Candidate Gene Discovery:
    After optional pre-assembly of the EST library using the DNA FIESTA pipeline (see below), the candidate protein
    dataset (the QueryDB) is used for translation and annotation (see below) of those ESTs with BLAST homology to 1+
    candidates. These translations are then assembled into consensus sequences where appropriate and alignments and
    trees made for the hits to each candidate protein. Finally, if desired, HAQESAC is run for each candidate protein.

EST Assembly:
    All EST assembly experiences two trade-offs: one between speed and accuracy, and a second between redundancy and
    accuracy. In particular, distinguishing sequencing errors from sequence variants (alleles) from different gene family
    members is not trivial.
    
    FIESTA is designed to provide straightforward assembly and BLAST-based annotation of EST sequences in Fasta format.
    The rationale behind its design was to try and optimise quality & redundancy versus comprehensive coverage for
    Proteomics identifications from Mass Spec data. FIESTA is designed to function in a relatively standalone capacity,
    with BLAST being the only other tool necessary. Due to this simplicity, FIESTA has some limitations; the main one
    being its inability to identify and deal with frameshift (indel) sequencing errors.
    
    FIESTA has two assembly and annotation pipelines: a protein pipeline based loosely on BUDAPEST and a DNA pipeline for
    "true" EST assembly. Details can be found in the Manual.

EST library assembly/annotation:
    In addition to the main functions, parts of the main FIESTA assembly/annotation pipeline can be run as standalone
    functions.

    ESTs can be converted to Reading Frames (RF) with est2rf=T:
    1. Identify orientation using 5' poly-T or 3' poly-A.
    - 1a. Where poly-AT tail exists, remove, translate in 3 forward RF and truncate at terminal stop codon.
    - 1b. Where no poly-AT tail exists, translate in all six RF.
    2. BLAST translations vs. search database with complexity filter on.
    - 2a. If EST has BLAST hits, retain RFs with desired e-value or better. 
    - 2b. If no BLAST hits, retain all RFs.

    Alternatively, translated RFs or other unannotated protein sequences can be given crude BLAST-based annotations using
    searchdb=FILE sequences with blastann=T. Note that these are simply the top BLAST hit and better annotation would be
    achieved using HAQESAC (or MultiHAQ for many sequences).
    
Commandline:
    ### ~ GENERAL INPUT ~ ###
    seqin=FILE      : EST file to be processed [None]
    fwdonly=T/F     : Whether to treat EST/cDNA sequences as coding strands (False = search all 6RF) [False]
    minpolyat=X     : Min length of poly-AT to be considered a poly AT [10]
    minorf=X        : Min length of ORFs to be considered [20]   
    blastopt=FILE   : File containing additional BLAST options for run, e.g. -B F [None]
    ntrim=X         : Trims of regions >= X proportion N bases [0.5]

    ### ~ SEQUENCE FORMATTING ~ ###
    gnspacc=T/F     : Convert sequences into gene_SPECIES__AccNum format wherever possible. [False]
    spcode=X        : Species code for EST sequences [None]
    species=X       : Species for EST sequences [None]
    newacc=X        : New base for sequence accession numbers ['' or spcode]

    ### ~ EST ASSEMBLY ~ ###
    minaln=X        : Min length of shared region for consensus assembly [40]
    minid=X         : Min identity of shared region for consensus assembly [95.0]
    bestorf=T/F     : Whether to use the "Best" ORF only for ESTs without BLAST Hits [True]
    pickup=T/F      : Whether to read in partial results and skip those sequences [True]
    annotate=T/F    : Annotate consensus sequences using BLAST-based approach [False]
    dna=T/F         : Implement DNA-based GABLAM assembly [True]
    resave=X        : Number of ESTs to remove before each resave of GABLAM searchdb [200]
    gapblast=T/F    : Whether to allow gaps during BLAST identification of GABLAM homologues [False]
    assmode=X       : Mode to use for EST assembly (nogab,gablam,oneqry) [oneqry]
    gabrev=T/F      : Whether to use GABLAM-based reverse complementation [True]

    ### ~ ANNOTATION ~ ###
    est2rf=T/F      : Execute BLAST-based EST to RF translation/annotation only, on seqin [False]
    est2haq=T/F     : Execute BLAST-based EST to RF translation/annotation on seqin followed by HAQESAC analysis [False]
    blastann=T/F    : Execute BLAST-based annotation of conensus translations only, on seqin [False]
    truncnt=T/F     : Whether to truncate N-terminal to Met in final BLAST annotation (if hit) [False]
    searchdb=FILE   : Fasta file for GABLAM search of EST translations [None]

    ### ~ QUERY SEARCH ~ ###
    batch=LIST      : List of EST libraries to search (will use seqin if none given) []
    querydb=FILE    : File of query sequences to search for in EST library [None]
    qtype=X         : Sequence "Type" to be used with NewAcc for annotation of translations [hit]
    assembly=T/F    : Assemble EST sequences prior to search [False]
    consensi=T/F    : Assemble hit ORF into consensus sequences [False]

    ### ~ HAQESAC OPTIONS ~ ###
    haqesac=T/F     : HAQESAC analysis of identified EST translations [True]
    multihaq=T/F    : Whether to run HAQESAC in two-phases [True]
    blastcut=X      : Reduced the number of sequences in HAQESAC runs to X (0 = no reduction) [50]
    cleanhaq=T/F    : Delete excessive HAQESAC results files [True]
    haqdb=FILELIST  : Optional extra databases to search for HAQESAC analysis []
    haqbatch=T/F    : Whether to only generate HAQESAC batch file (True) or perform whole run (False) [False]
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import multihaq, rje, rje_seq, rje_sequence, rje_tree, rje_zen
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added FIESTA pipeline for protein-based clustering in addition to TIGR based (partial) method
    # 0.2 - Removed TIGR pipeline and replaced with DNA version of FIESTA.
    # 0.3 - Added annotation method and extra mapping.
    # 0.4 - Added oneqry method for GABLAM consensus generation
    # 1.0 - Added querydb search option
    # 1.1 - Added assmode option = Mode to use for EST assembly (nogab,gablam,oneqry) [oneqry]
    # 1.2 - Added FwdOnly option for EST annotation.
    # 1.3 - Add HAQESAC run following annotation.
    # 1.4 - Tidied and modified QueryESTs analysis.
    # 1.5 - Bug removal and additional tidying for MultiHAQ and annotateEST methods.
    # 1.6 - Removed HAQESAC import (uses MultiHAQ).
    # 1.7 - Updated to use rje_blast_V2. Needs work to make function with BLAST+.
    # 1.8 - Minor crash fixes. Updated more functions to work with BLAST+.
    # 1.8.1 - Replaced type with stype throughout to try and avoid TypeError crashes.
    # 1.9.0 - Altered HAQDB to be a list of files rather than just one.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Make first functioning version
    # [?] : ntruncblast=T/F : Whether to truncate ORF N-termini at the extremes of BLAST local matches [True]
    # [Y] : Scale up to run on IRIDIS
    # [Y] : Add startfrom/pickup
    # [?] : Search with PFam HMMs?
    # [Y] : Better tracing of ESTs through consensus/assembly process
    # [Y] : Implement resave option.
    # [Y] : Implement DNA option.
    # [Y] : Implement Mapping
    # [Y] : Add annotate=T/F and options and split annotate from fiesta() as individual method.
    # [Y] : Make a more efficient DNA GABLAM consensi that assembles in chunks and re-searches with refseq?
    # [Y] : Check correct usage of reverse complementation in GABLAM One-Query
    # [ ] : Add options for GABLAM-free, GABLAM & GABLAM One-Query and benchmark
    # |-- [ ] : Add to non-DNA assembly version.
    # [ ] : Add proper PickUp option for DNA assembly
    # [Y] : Add trimming on N-regions from the ends.
    # [Y] : Add proper multihaq run following annotation as performed by BUDAPEST in full auto mode.
    # [ ] : Update MultiHAQ method to run using multihaq module, as for queryEST multiHAQ.
    # [ ] : Add forking
    # [ ] : Add protein sequence assembly and annotation
    # [ ] : Upgrade to use BLAST+
    # [ ] : Improve documentation.
    # [ ] : Add an option to split hits if spanning an ORF w/o any local hits. (newblast only?)
    # [ ] : Make sure that seqin=ACCLIST is acceptable input.
    # [ ] : Add option for a BLAST database to be purely a taxid: download from Uniprot.
    # [ ] : Modify to handle ESTs from multiple species.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('FIESTA', '1.9.0', 'November 2014', '2009')
    description = 'Fasta Input EST Analysis'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: FIESTA Class                                                                                            #
#########################################################################################################################
class FIESTA(rje.RJE_Object):     
    '''
    FIESTA Class. Author: Rich Edwards (2009).

    Info:str
    - AssMode = Mode to use for EST assembly (nogab,gablam,oneqry) [oneqry]
    - NewAcc = New base for sequence accession numbers ['']
    - NewGene = New gene to be used for EST sequences that have accession numbers only ['est']
    - QType = Sequence "Type" to be used with NewAcc for annotation of translations [HIT]
    - QueryDB = File of query sequences to search for in EST library [None]
    - SearchDB = Fasta file for GABLAM search of EST translations [None]
    - SeqIn = Input EST file [None]
    - SpCode = Species code for EST sequences [None]
    - Species = Species for EST sequences [None]
    
    Opt:boolean
    - Annotate = Annotate consensus sequences using BLAST-based approach [False]
    - Assembly = Assemble EST sequences prior to search [False]
    - BestORF = Whether to use the "Best" ORF only for ESTs without BLAST Hits [True]
    - BLASTAnn = Execute BLAST-based annotation of conensus translations only, on seqin [False]
    - Consensi = Assemble hit ORF into consensus sequences [False]
    - DNA = Implement DNA-based GABLAM assembly [True]
    - EST2HAQ = Execute BLAST-based EST to RF translation/annotation on seqin followed by HAQESAC analysis [False]
    - EST2RF = Execute BLAST-based EST to RF translation/annotation only on seqin [False]
    - FwdOnly = Whether to treat EST/cDNA sequences as coding strands (False = search all 6RF) [False]
    - GabRev = Whether to use GABLAM-based reverse complementation [True]
    - GapBLAST = Whether to allow gaps during BLAST identification of GABLAM homologues [False]
    - GnSpAcc = Convert sequences into gene_SPECIES__AccNum format wherever possible. [False]
    - HAQESAC = HAQESAC analysis of identified EST translations [True]
    - HAQBatch = Whether to only generate HAQESAC batch file (True) or perform whole run (False) [False]
    - MultiHAQ = Whether to run HAQESAC in two-phases [True]
    - CleanHAQ = Delete excessive HAQESAC results files [True]
    - PickUp = Whether to read in partial results and skip those sequences [True]
    - TruncNT = Whether to truncate N-terminal to Met in final BLAST annotation [False]

    Stat:numeric
    - BlastCut = Reduced the number of sequences in HAQESAC runs to X (0 = no reduction) [50]
    - ESTNum = Number of input EST sequences (used for naming preZero)
    - MinPolyAT = Min length of poly-AT to be considered a poly AT [10]
    - MinORF = Min length of ORFs to be considered [20]   
    - MinID = Min identity of shared BLAST region [95.0]
    - MinAln = Min length of shared BLAST region [40]
    - NTrim = Trims of regions >= X proportion N bases [0.5]
    - ReSave = Number of ESTs to remove before each resave [200]

    List:list
    - Batch = List of EST libraries to search (will use seqin if none given) []
    - HAQDB = Optional extra database to search for HAQESAC analysis [None]

    Dict:dictionary
    - Counter = Dictionary counting consensi of different types
    - Mapping = Dictionary used to map between IDs {ID:{Type:ID}}
    - PickUp = Partial results (ESTs) to skip {EST:[Trans]}

    Obj:RJE_Objects
    - HAQESAC = HAQESAC object for analysis of identified EST translations 
    - SeqList = rje_seq.SeqList object containing the EST sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['NewAcc','SearchDB','QueryDB','SeqIn','SpCode','Species','QType','AssMode']
        self.optlist = ['Annotate','DNA','EST2RF','BestORF','PickUp','BLASTAnn','TruncNT','Assembly','HAQESAC',
                        'MultiHAQ','GnSpAcc','GapBLAST','GabRev','FwdOnly','EST2HAQ','CleanHAQ','HAQBatch','Consensi']
        self.statlist = ['MinPolyAT','MinID','MinAln','Unmatched','ReSave','MinORF','ESTNum','NTrim','BlastCut']
        self.listlist = ['Headers','Batch','HAQDB']
        self.dictlist = ['Counter','Mapping','PickUp']
        self.objlist = ['SeqList','HAQESAC']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'QType':'hit','AssMode':'oneqry'})
        self.setOpt({'BestORF':True,'PickUp':True,'TruncNT':False,'HAQESAC':True,'MultiHAQ':True,'GabRev':True,
                     'DNA':True,'CleanHAQ':True})
        self.cmd_list = ['dna=T','ntrim=0.5'] + self.cmd_list     # Add this default for other objects.
        self.setStat({'MinPolyAT':10,'MinID':95.0,'MinAln':40,'Unmatched':25,'ReSave':200,'MinORF':20,'NTrim':0.5,
                      'BlastCut':50})
        ### Other Attributes ###
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ###
                self._cmdReadList(cmd,'info',['NewAcc','Species','SpCode','QType','AssMode'])
                self._cmdReadList(cmd,'file',['SearchDB','QueryDB','SeqIn'])
                self._cmdReadList(cmd,'int',['MinPolyAT','MinAln','Unmatched','ReSave','BlastCut'])
                self._cmdReadList(cmd,'stat',['MinID','NTrim'])
                self._cmdReadList(cmd,'opt',['Annotate','DNA','BestORF','PickUp','EST2RF','BLASTAnn','TruncNT','EST2HAQ',
                                             'HAQESAC','Assembly','MultiHAQ','GnSpAcc','GapBLAST','GabRev','FwdOnly',
                                             'CleanHAQ','HAQBatch','Consensi'])
                self._cmdReadList(cmd,'glist',['Batch','HAQDB'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Pipeline Methods                                                                                   #
#########################################################################################################################
    def run(self):  ### Main run method.
        '''Main run method. Selects TIGR or FIEST pipeline.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MinID'] > 1:  self.stat['MinID'] /= 100.0
            self.setupSeqNames()
            ### ~ [2] ~ Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('QueryDB'): self.queryESTs()
            elif self.opt['DNA'] and not (self.opt['EST2RF'] or self.opt['BLASTAnn'] or self.getBool('EST2HAQ')): self.fiestaDNA() 
            else: self.fiesta()
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupSeqNames(self):    ### Setup species codes and NewAcc base
        '''Setup species codes and NewAcc base.'''
        try:### ~ [1] ~ Species codes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Species'].lower() == 'none': self.info['Species'] = ''
            if self.info['NewAcc'].lower() == 'none': self.info['NewAcc'] = ''
            if self.info['SpCode'].lower() in ['none','']: self.info['SpCode'] = rje_sequence.getSpecCode(self.info['Species'])
            if not self.info['SpCode'] and self.info['NewAcc']: self.info['SpCode'] = self.info['NewAcc']
            if not self.info['Species'] and self.info['SpCode']: self.info['Species'] = self.info['SpCode']
            ### ~ [2] ~ NewAcc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.info['NewAcc']: 
                if len(string.split(self.info['Species'])) > 1:
                    [genus,species] = string.split(self.info['Species'].upper())[:2]
                    self.info['NewAcc'] = '%s%s' % (genus[:1],species[:3])
                elif self.info['SpCode']: self.info['NewAcc'] = self.info['SpCode']
                else: self.info['Species'] = 'Unknown'; self.info['SpCode'] = 'UNK'
            ### ~ [3] ~ QType gene naming ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['QType'].lower() in ['','none'] or self.info['QueryDB'].lower() in ['','none']: self.info['QType'] = self.info['NewAcc'].lower()
        except: self.errorLog('Problem setting up species/NewAcc')
#########################################################################################################################
    def queryESTs(self):    ### Analysis pipeline to query ESTs for specific proteins
        '''Analysis pipeline to query ESTs for specific proteins.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RUN','Executing FIESTA candidate EST discovery pipeline')
            if not rje.checkForFile(self.getStr('QueryDB')): return self.errorLog('QueryDB "%s" not found!' % self.getStr('QueryDB'),printerror=False)
            if not self.getStrLC('Basefile'): self.basefile(rje.baseFile(self.getStr('QueryDB'),strip_path=True))
            basefile = self.basefile()
            self.setBool({'Annotate':False,'DNA':True})
            ## ~ [0a] ~ Sort out files to be used for various searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('SearchDB'):
                searchdb = self.getStr('SearchDB')
                if self.getBool('HAQESAC'): self.printLog('#HAQDB','SearchDB %s will be used for HAQESAC.' % searchdb)
            else: searchdb = None
            if self.list['HAQDB'] and self.getBool('HAQESAC'): self.printLog('#HAQDB','%d HAQDB will be used for HAQESAC.' % len(self.list['HAQDB']))
            #if self.info['HAQDB'].lower() not in ['','none'] and self.info['HAQDB'] != self.info['SearchDB']:
            #    self.printLog('#HAQDB','HAQDB replaced with SearchDB for Query EST analysis')
            self.setStr({'SearchDB':self.getStr('QueryDB')})
            self.printLog('#SEARCH','QueryDB %s will be used for EST annotation SearchDB' % self.getStr('QueryDB'))
            if not self.list['Batch']: self.list['Batch'] = [self.getStr('SeqIn')]

            #!# Add looking for previous steps and skipping if force=F (can ask if i>=1) #!#

            ### ~ [1] ~ Optional EST Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Assembly']:
                batchin = self.list['Batch'][0:]; self.list['Batch'] = []
                for seqin in batchin:
                    self.basefile('%s.%s' % (basefile,rje.baseFile(seqin)))
                    consfas = '%s.cons.fas' % self.basefile()
                    if self.force() or not rje.exists(consfas) or (self.i() > 0 and rje.yesNo('%s found. Repeat assembly? (force=F)' % consfas,default='N')):
                        self.fiestaDNA()
                    self.list['Batch'].append(consfas)

            ### ~ [2] ~ Search with query sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batchres = []   # List of Batch results files
            for seqin in self.list['Batch']:
                self.basefile('%s.%s' % (basefile,rje.baseFile(seqin,strip_path=True)))
                self.setStr({'SeqIn':seqin})
                self.printLog('#SEQIN',seqin)
                ## ~ [2a] ~ Generate *.est_hits.fas from BLAST search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                esthits = '%s.est_hits.fas' % self.basefile()
                seqlist = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','dna=T','autofilter=F','seqnr=F','logrem=F','replacechar=T'])
                seqlist.setStr({'Name':esthits})
                self.obj['BLAST'] = blast = rje_blast.blastObj(log=self.log,cmd_list=self.cmd_list+['blastd=%s' % self.getStr('SeqIn'),'blasti=%s' % self.getStr('QueryDB')])
                if self.force() or not rje.exists(esthits) or (self.i() > 0 and rje.yesNo('%s found. Repeat BLAST extraction? (force=F)' % esthits,default='N')):
                    blast.setStr({'Name':'%s.blast' % self.basefile(),'Type':'tblastn'})
                    blast.setInt({'Verbose': self.v() - 1,'HitAln':0})
                    blast.formatDB(protein=False,force=False)
                    blast.blast(use_existing=not self.force())
                    blast.readBLAST(clear=True,unlink=True,gablam=False,local=False)
                    blast.hitToSeq(seqlist)
                    seqlist._checkForDup()
                    seqlist.saveFasta()
                else: seqlist.loadSeqs(esthits,seqtype='DNA')
                ## ~ [2b] ~ Annotate EST translations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## >> Sort out seqnames ... ###
                if seqlist.seqs()[0].getStr('SpecCode') != 'UNK':
                    self.setStr({'SpCode':seqlist.seqs()[0].getStr('SpecCode'),'Species':seqlist.seqs()[0].getStr('Species')})
                if seqlist.seqs()[0].getStr('AccNum')[4:8] == 'CONS': self.setStr({'NewAcc':seqlist.seqs()[0].getStr('AccNum')[:4]})
                self.setupSeqNames()
                ## >> Annotate ... ###
                self.annotateEST(seqlist)   # This should pickup incomplete runs and finish only.
                self.loadMapping('Trans')   # Load the mapping from the *.tdt file
                #!# Add pre- post- and both- assembly options #!#
                ## ~ [2c] ~ Optional post-translation consensus assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.setBool({'DNA':False})
                seqlist.setBool({'DNA':False})
                seqlist.loadSeqs('%s.trans.fas' % self.basefile(),seqtype='Protein')
                qtype = self.getStrLC('QType')
                if self.getBool('Consensi'):
                    self.warnLog('Post-translation assembly not yet tested!')
                    self.dict['Counter'][qtype.upper()] = 0
                    self.makeConsensi(seqlist,log=self.getStr('QType'),stype=qtype.upper())
                    seqlist.saveFasta(seqfile='%s.%s.fas' % (self.basefile(),qtype),append=False,log=True)
                    self.saveMapping(qtype.upper())
                    seqlist.loadSeqs('%s.%s.fas' % (self.basefile(),qtype),seqtype='Protein')
                self.annotateTrans(seqlist,stype='final')    # This should pickup incomplete runs and finish only.
                if qtype != 'final':
                    rje.fileTransfer('%s.final.fas' % self.basefile(),'%s.%s.fas' % (self.basefile(),qtype),deletefrom=True,append=False)
                    self.printLog('#%s' % qtype.upper(),'%s sequences copied to from %s to %s.' % (seqlist.seqNum(),'%s.final.fas' % self.basefile(),'%s.%s.fas' % (self.basefile(),qtype)))
                batchres.append('%s.%s.fas' % (self.basefile(),qtype))
                #!# Add option to delete files? #!#
                #for bfile in ['trans.tdt','final.tdt'.'finals.fas']:
                #    if os.path.exists('%s.%s' % (self.basefile(),bfile)): os.unlink('%s.%s' % (self.basefile(),bfile))

            ### ~ [3] ~ Generate Alignments and Trees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#TREES','## ~~~~~~~ ALIGNMENT & TREE GENERATION ~~~~~~~ ##')
            self.basefile(basefile)
            seqlist.setStr({'HAQBAT':''})
            seqlist.seq = []
            seqlist.loadSeqs(self.info['QueryDB'],seqtype='Protein')
            #self.deBug(seqlist.seq[0].info['Sequence'])
            seqlist.list['Blast2Fas'] = [self.info['QueryDB']] + batchres
            seqlist.setBool({'DNA':False})
            seqlist.cmd_list += ['blastp=blastp','dna=F']
            #!# Add capacity to skip if already existing! Make optional? #!#
            rje_seq.Blast2Fas(seqlist)
            rje.mkDir(self,'%s_ALN/' % self.info['Basefile'])
            rje.mkDir(self,'%s_TREE/' % self.info['Basefile'])
            for seq in seqlist.seqs():
                ## ~ [2a] ~ Align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                infile = '%s.blast.fas' % seq.info['AccNum']
                hitseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % infile,'autoload=T','autofilter=F','query=%s' % seq.info['AccNum'],'dna=F'])
                hitseq.seq = hitseq.seq
                afile = '%s%s.aln.fas' % (rje.makePath('%s_ALN/' % self.info['Basefile']),seq.info['AccNum'])
                hitseq.align(outfile=afile)
                os.unlink(infile)
                ## ~ [2b] ~ Make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if hitseq.seqNum() < 4: self.printLog('#TREE','%s: too few hits for tree construction.' % seq.shortName())
                else:
                    tfile = '%s%s.nsf' % (rje.makePath('%s_TREE/' % self.info['Basefile']),seq.info['AccNum'])
                    tree = rje_tree.Tree(self.log,cmd_list=['root=mid']+self.cmd_list+['autoload=F','savetree=%s' % tfile])
                    tree.makeTree(hitseq,keepfile=False)
                    if tree.opt['OutputBranchLen']: withbranchlengths = 'Length'
                    else: withbranchlengths = 'none'
                    outnames = tree.info['OutNames']
                    maxnamelen = tree.stat['TruncNames']
                    tree.info['Basefile'] = rje.baseFile(tfile)
                    try: tree.saveTrees(seqname=outnames,blen=withbranchlengths)
                    except: tree.saveTree(seqnum=True,seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            self.printLog('#HIT','Alignments and trees made of Queries and translated EST hits')

            ### ~ [4] ~ Perform MultiHAQ with Queries and SearchDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#TREES','## ~~~~~~~ MULTIHAQ CANDIDATE ANALYSIS ~~~~~~~ ##')
            haqcmd = ['ac=F','blastcut=50','maketree=fasttree','bootstraps=1000','alnprog=maaft','qcover=0','qryid=0','allowvar=F','treeformats=qspec,text,nwk','group=dup']
            haqcmd += self.cmd_list + ['root=mid','qryvar=T','dna=F']
            #mhaq = multihaq.MultiHAQ(self.log,self.cmd_list+['dna=F','blastcut=%d' % self.getInt('BlastCut')])
            mhaq = multihaq.MultiHAQ(self.log,haqcmd)
            mhaq.obj['SeqList'] = seqlist
            mhaq.info['HaqDir'] = rje.makePath('%s_HAQESAC/' % self.info['Basefile'])
            rje.mkDir(self,mhaq.info['HaqDir'])
            seqlist.seq = []
            for resfile in batchres: seqlist.loadSeqs(resfile,seqtype='Protein',nodup=False,clearseq=False)
            seqlist.info['Name'] = '%s.final.fas' % self.info['Basefile']
            seqlist._checkForDup(True)
            seqlist.saveFasta()
            mini = mhaq.cmd_list + ['seqin=%s' % seqlist.name(),'blast2fas=','haqdir=%s' % mhaq.info['HaqDir']]
            open('%s.multihaq.ini' % self.baseFile(),'w').write(string.join(mini,'\n'))
            self.printLog('#INI','MultiHAQ INI File %s.multihaq.ini generated (for manual re-run w/o blast2fas).' % self.baseFile())
            if rje.checkForFile(searchdb): batchres.append(searchdb)
            batchres += self.list['HAQDB']
            seqlist.list['Blast2Fas'] = [self.info['QueryDB']] + batchres
            mhaq.blast2fas()
            mhaq.haqBatch()
            mhaq.multiHAQ()
        except: os.chdir(self.info['RunPath']); self.errorLog('Error during queryESTs()')
#########################################################################################################################
    ### <3> ### Annotation Pipeline                                                                                     #
#########################################################################################################################
    def fiesta(self):   ### Main FIESTA analysis pipeline
        '''
        FIESTA analysis pipeline based loosely on BUDPAEST:
        1. Identify orientation using 5' poly-T or 3' poly-A.
        1a. Where poly-AT tail exists, remove, translate in 3 forward RF and truncate at terminal stop codon
        1b. Where no poly-AT tail exists, translate in all six RF
        2. BLAST translations vs. EMBL metazoa with complexity filter on.
        3. If EST has BLAST hits, retain RFs with desired e-value or better. If no BLAST hits, retain all RFs.
        4. All-by-all GABLAM to identify relationships between translations.
        5. Assemble in order of best e-vlaue using BLAST alignment and %ID cutoff.
        6. Align assembled translations. (Or use BLAST alignments?)  #!# Could be tricky #!#
        7. Generate a consensus for each assembly.
        8. Generate XML file with all assembly data plus a FASTA file for searching.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['EST2HAQ']:
                self.printLog('#RUN','Executing FIESTA EST translation/annotation and HAQESAC only')
                #!# Check species code etc. #!#
            elif self.opt['EST2RF']: self.printLog('#RUN','Executing FIESTA EST translation/annotation only'); self.info['QType'] = ''
            elif self.opt['BLASTAnn']: self.printLog('#RUN','Executing FIESTA BLAST-based translation annotation only')
            else: self.printLog('#RUN','Executing FIESTA assembly/annotation pipeline')
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F','accnr=F','seqnr=F','replacechar=F'])
            if not rje.exists(seqlist.info['Name']): return self.errorLog('Input file "%s" not found' % seqlist.info['Name'],printerror=False)
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = '%s.fiesta' % rje.baseFile(seqlist.info['Name'])
            blast = self.obj['BLAST'] = rje_blast.blastObj(self.log,self.cmd_list,type='Dev')
            if not blast.getBool('OldBLAST'): self.warnLog('FIESTA does not really work with oldblast=F. Use old BLAST and oldblast=T.',quitchoice=True)
            blast.setStat({'Verbose':self.stat['Verbose'] - 1})
            self.list['Headers'] = ['Trans','EST','TopHit']

            ### ~ [2] Convert ESTs to RF(s) and make consensus translations based on TopHit clusters ~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Convert to RFs with BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['BLASTAnn']: self.annotateTrans(seqlist); return
            if self.opt['EST2HAQ'] or self.opt['EST2RF']: self.annotateEST(seqlist,stype=self.getStr('NewAcc'))
            else: self.annotateEST(seqlist,stype='trans')
            if self.opt['EST2HAQ']: return self.multiHAQ(seqlist,batchonly=self.getBool('HAQBatch'))
            if self.opt['EST2RF']: return
            ## ~ [2b] ~ Cluster based on shared BLAST hits and make consensi ~~~~~~~~~~~~~~~~~~~~~~ ##
            self.clusterTrans()     # Read in fiestatdt to get shared hits and BLAST names for FASTACMD etc.
            (prev,next) = ('trans','top')
            fasin = '%s.%s.fas' % (self.info['Basefile'],prev)
            fasout = '%s.%s.fas' % (self.info['Basefile'],next); tdtout = '%s.%s.tdt' % (self.info['Basefile'],next)
            if self.opt['Force']: rje.backup(self,fasout); rje.backup(self,tdtout)
            if not os.path.exists(fasout) or self.opt['Force']:     
                (cx,ctot,prex,sx) = (0.0,len(self.dict['Clusters']),0,0)
                blast.formatDB(fasin,protein=True,log=self.log)  
                for cluster in self.dict['Clusters']:
                    self.progLog('\r#CONS','Assembling consensi via shared top BLAST hits: %.2f%% (%s -> %s)' % (cx/ctot,rje.integerString(prex),rje.integerString(sx))); cx += 100.0
                    seqlist.seq = []
                    for trans in self.dict['Clusters'][cluster]['Trans']: seqlist.seqFromFastaCmd(trans,fasin)
                    prex += seqlist.seqNum(); self.makeConsensi(seqlist,stype='TOP'); sx += seqlist.seqNum()
                    seqlist.saveFasta(seqfile=fasout,append=True,log=False)
                self.printLog('\r#CONS','%s consensi assembled from %s translations via shared top BLAST hits.' % (rje.integerString(sx),rje.integerString(prex)))
                self.saveMapping('TOP')
            else: self.loadMapping('TOP')

            ### ~ [5] ~ GABLAM and assemble Hits, NoHits, then all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [5a] ~ GABLAM Consensi with BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (prev,next) = ('top','hit')
            fasin = '%s.%s.fas' % (self.info['Basefile'],prev)
            fasout = '%s.%s.fas' % (self.info['Basefile'],next); tdtout = '%s.%s.tdt' % (self.info['Basefile'],next)
            if self.opt['Force']: rje.backup(self,fasout); rje.backup(self,tdtout)
            if not os.path.exists(fasout) or self.opt['Force']:     
                seqlist.loadSeqs(fasin,filetype='fas',seqtype='protein',aln=False,nodup=False)
                self.gablamConsensiOneQry(seqlist,fasin,fasout,stype='HIT')
                self.saveMapping('HIT')
            else: self.loadMapping('HIT')
            ## ~ [5b] ~ GABLAM Consensi without BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            (prev,next) = ('nohit','non')
            fasin = '%s.%s.fas' % (self.info['Basefile'],prev)
            fasout = '%s.%s.fas' % (self.info['Basefile'],next); tdtout = '%s.%s.tdt' % (self.info['Basefile'],next)
            if self.opt['Force']: rje.backup(self,fasout); rje.backup(self,tdtout)
            if not os.path.exists(fasout) or self.opt['Force']:     
                seqlist.seq = []
                (tx, ttot) = (0.0, len(self.list['NoHits']))
                #for trans in self.list['NoHits']:
                    #self.progLog('\r#NOHIT','Extracting RFs w/o BLAST Hits: %.2f%% (%s seq)' % (tx/ttot,rje.integerString(seqlist.seqNum()))); tx += 100.0
                    #try: seqlist.seqFromFastaCmd(trans,'%s.trans.fas' % (self.info['Basefile']))
                    #except: self.errorLog('Problem with NoHits fastacmd for "%s"' % trans)
                seqlist.loadSeqs('%s.trans.fas' % (self.info['Basefile']),filetype='fas',seqtype='protein',aln=False,nodup=False)
                seqlist.opt['LogRem'] = False
                seqlist.list['GoodAcc'] = self.list['NoHits']
                self.progLog('\r#NOHIT','Extracting RFs w/o BLAST Hits ...')
                seqlist._filterSeqs(); seqlist.list['GoodAcc'] = []
                self.printLog('\r#NOHIT','Extracted %s RFs w/o BLAST Hits from %s.' % (rje.integerString(seqlist.seqNum()),'%s.trans.fas' % (self.info['Basefile'])))
                seqlist.saveFasta(seqfile=fasin)
                self.gablamConsensiOneQry(seqlist,fasin,fasout,stype='NON')
                self.saveMapping('NON')
            else: self.loadMapping('NON')
            ## ~ [5c] ~ GABLAM Consensi for all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            next = 'con'
            fasout = '%s.%s.fas' % (self.info['Basefile'],next); tdtout = '%s.%s.tdt' % (self.info['Basefile'],next)
            if self.opt['Force']: rje.backup(self,fasout); rje.backup(self,tdtout)
            if not os.path.exists(fasout) or self.opt['Force']:     
                combcons = '%s.combcons.fas' % self.info['Basefile']
                open(combcons,'w').write(open('%s.hit.fas' % self.info['Basefile'],'r').read())
                open(combcons,'a').write(open('%s.non.fas' % self.info['Basefile'],'r').read())
                seqlist.loadSeqs(combcons,filetype='fas',seqtype='protein',aln=False,nodup=False)
                self.gablamConsensiOneQry(seqlist,combcons,fasout)    #!# Consider complexity filter at this stage? #!#
                self.saveMapping('CON')
            else: self.loadMapping('CON')

            ### ~ [6] ~ Final BLAST-based identification ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Annotate']: return
            seqlist.info['Name'] = fasout
            self.annotateTrans(seqlist,stype='final')
            #!# Final BLAST search should truncate N-terminal if BLAST hit starts with M and included in alignment? #!#            
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def convertESTtoRF(self,seq,stype='trans'):       ### Converts a given EST sequence to RFs and BLASTs to find best sequence.
        '''
        Converts a given EST sequence to RFs and BLASTs to find best sequence. First translate, then BLAST against given
        search database. Read in alignments and map onto sequence. Identify ORFs with some coverage and reduce sequence
        to not include any ORFs without BLAST coverage. Give relevant sequence a new description using annotation and
        add sequences to self.dict['RF'].
        >> acc:str = Hit identifier corresonding to self.dict['Hit']
        >> seq:Sequence object = EST sequences as rje_sequence.Sequence object
        '''
        try:### ~ [1] ~ First translate, then BLAST against given search database. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            est = seq.info['AccNum']
            if self.info['QType']: seq.info['Name'] = '%s_%s__%s' % (self.info['QType'],self.info['SpCode'],est)
            rfblast = self.obj['BLAST']
            rfblast.setBool({'GappedBLAST':True})
            fiestafas = '%s.%s.fas' % (self.info['Basefile'],stype.lower())
            fiestatdt = '%s.%s.tdt' % (self.info['Basefile'],stype.lower())
            ## ~ [1a] ~ Generate sequence files of appropriate translated sequences ~~~~~~~~~~~~~~~ ##
            myrfs = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','autofilter=F','dna=F'])
            myrfs.info['Name'] = '%s.%s.tmp.fas' % (est,rje.randomString(6))
            myrfs.opt['ReplaceChar'] = False
            if len(seq.getSequence()) < (3 * self.stat['MinORF']):
                open(fiestatdt,'a').write('%s\n' % string.join([est,'0','-','-','-'],'\t'))
                return 0
            rftrans = rje_sequence.estTranslation(seq.getSequence(),self.stat['MinPolyAT'],self.opt['FwdOnly'])
            for rf in rje.sortKeys(rftrans):
                trans = rftrans[rf]
                if len(rje_sequence.bestORF(trans)) < self.stat['MinORF']: continue     # Too crap!
                myrfs._addSeq('%s_RF%s' % (seq.shortName(),rf),trans)
            ## ~ [1b] ~ BLAST sequence against search database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not myrfs.seq:
                open(fiestatdt,'a').write('%s\n' % string.join([est,'%d' % myrfs.seqNum(),'-','-','-'],'\t'))
                return 0
            myrfs.saveFasta(log=False)
            rfblast.setInfo({'InFile':myrfs.info['Name'],'Name':'%s.blast' % rje.baseFile(myrfs.info['Name']),'Type':'blastp'})
            rfblast.blast(log=False); os.unlink(myrfs.info['Name'])
            try: rfblast.readBLAST(clear=True,screen=False,unlink=True,local=True,gablam=False)
            except: self.errorLog('%s problem.' % rfblast.getStr('Name'))

            ### ~ [2] ~ Read in alignments and map onto sequence. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rfblast.getBool('OldBLAST'):
                hitx = rfblast.hitNum()                                             # Total number of BLAST hits
                bestrf = None; beste = 1000
                for rfseq in myrfs.seqs():
                    ### ~ [2a] ~ Deal with lack of BLAST Hits first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if not hitx:        # No BLAST Hits at all
                        if self.opt['BestORF']:     # Reduce to best ORF and pick best of the best!
                            rfseq.info['Sequence'] = rje_sequence.bestORF(rfseq.getSequence(),startm=False,nonx=True)
                            if not bestrf or rfseq.nonX() > bestrf.nonX(): bestrf = rfseq
                        rfseq.info['Name'] = '%s No BLAST hit (e<%s) to %s' % (rfseq.shortName(),rje.expectString(rfblast.getNum('E-Value')),os.path.basename(rfblast.getStr('DBase')))
                        continue
                    ## ~ [2b] ~ Deal with BLAST Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    query = rfseq.shortName()
                    if rfblast.queryHits(query): hit = rfblast.queryHits(query)[0]    # Should be best hit
                    else: continue   # Not best RF
                    hdict = rfblast.hitData(hit,query)
                    rfseq.info['Name'] = '%s Similar (e=%s) to %s' % (rfseq.shortName(),rje.expectString(hdict['E-Value']),hit)
                    if not bestrf or hdict['E-Value'] < beste: (bestrf,beste) = (rfseq,hdict['E-Value'])
                    ## ~ [2c] ~ Identify ORFs with some coverage and reduce sequences to covered ORFs ~ ##
                    orfseq =  rfseq.getSequence()
                    qstart = len(orfseq)+1; qend = -1
                    for ldict in rfblast.localData(hit,query):  # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq']
                        #self.debug(ldict)
                        qstart = min(qstart,ldict['QryStart']-1)     # 0-L
                        qend = max(qend,ldict['QryEnd'])
                        #!# Add code (an option) to catch an unmatched ORF in the middle: 2+ hits in same RF?!
                    if qend < qstart: raise ValueError('Failed to map local hit data for %s vs %s!' % (query,hit))
                    prevorf = string.split(orfseq[:qstart],'*')[-1]
                    postorf = string.split(orfseq[qend:],'*')[0]
                    #self.debug(qstart); self.debug(qend)
                    rfseq.info['Sequence'] = prevorf + orfseq[qstart:qend] + postorf
                    #self.deBug(orfseq)
                    #self.deBug(rfseq.info['Sequence'])

                ### ~ [3] ~ Select/Output the best (if there is one) or all if no best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if bestrf: myrfs.seq = [bestrf]
                for rfseq in myrfs.seqs(): open(fiestafas,'a').write('>%s\n%s\n' % (rfseq.info['Name'],rfseq.getSequence()))
                if bestrf and hitx:     # Best RF has BLAST hit
                    query = bestrf.shortName()
                    hdict = rfblast.hitData(rfblast.queryHits(query)[0],query)
                    open(fiestatdt,'a').write('%s\n' % string.join([est,'%d' % myrfs.seqNum(),query,hdict['Hit'],rje.expectString(hdict['E-Value'])],'\t'))
                else:
                    for rfseq in myrfs.seqs():
                        open(fiestatdt,'a').write('%s\n' % string.join([est,'%d' % myrfs.seqNum(),rfseq.shortName(),'-','-'],'\t'))

            else:
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','autofilter=F','dna=F'])
                searchdict = rfblast.searchSeq(myrfs,proglog=False,inverse=True)    # Return dictionary of sequence:search
                hitx = rfblast.hitNum()                                             # Total number of BLAST hits
                gindex = ['-','X','|','+']
                bestrf = None; beste = 1000
                for rfseq in myrfs.seqs():
                    try: search = searchdict[rfseq]
                    except:
                        print searchdict
                        rfseq.info['Name'] = '%s (BLAST error)' % rfseq.shortName()
                        self.errorLog('BLAST search for %s missing!' % rfseq.shortName(),printerror=False)
                        search = None
                    ### ~ [2a] ~ Deal with lack of BLAST Hits first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if not hitx:        # No BLAST Hits at all
                        if self.opt['BestORF']:     # Reduce to best ORF and pick best of the best!
                            rfseq.info['Sequence'] = rje_sequence.bestORF(rfseq.getSequence(),startm=False,nonx=True)
                            if not bestrf or rfseq.nonX() > bestrf.nonX(): bestrf = rfseq
                        rfseq.info['Name'] = '%s No BLAST hit (e<%s) to %s' % (rfseq.shortName(),rje.expectString(rfblast.getStat('E-Value')),os.path.basename(rfblast.info['DBase']))
                        continue
                    ## ~ [2b] ~ Deal with BLAST Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if hitx and not search or not search.hitNum(): continue   # Not best RF
                    search.gablam(keepaln=True)  # GABLAM Alignments in hit.dict['GABLAM']['Aln']['Qry','Hit','QryO','HitO']
                    hit = search.hit[0]
                    #self.deBug(hit.aln)
                    #for gkey in rje.sortKeys(hit.dict['GABLAM']): self.deBug('%s: %s' % (gkey,hit.dict['GABLAM'][gkey]))
                    try: hitname = seqlist.seqFromFastaCmd(hit.info['Name'],rfblast.info['DBase']).info['Name']
                    except: hitname = hit.info['Name']
                    rfseq.info['Name'] = '%s Similar (e=%s) to %s' % (rfseq.shortName(),rje.expectString(hit.stat['E-Value']),hitname)
                    if not bestrf or hit.stat['E-Value'] < beste: (bestrf,beste) = (rfseq,hit.stat['E-Value'])
                    galn = hit.dict['GABLAM']['Aln']['Qry']    # Replace with '|', '+' or 'X'
                    rfseq.info['GAln'] = string.join(galn,'')
                    #self.deBug(rfseq.info['GAln'])
                    ## ~ [2c] ~ Identify ORFs with some coverage and reduce sequences to covered ORFs ~ ##
                    orfs = string.split(rfseq.getSequence(),'*')
                    (start,end) = (0,0)
                    (firstcov,lastcov) = (-1,-1)
                    for i in range(len(orfs)):
                        orf = orfs[i]
                        end += len(orf)
                        orfgaln = rfseq.info['GAln'][start:end]
                        #self.deBug('%d: %s' % (i,orfgaln))
                        if string.count(orfgaln,'-') != len(orfgaln):   # At least partial coverage by 1+ hits
                            if firstcov < 0: firstcov = i
                            lastcov = i
                        #!# else: # Add code (an option) to catch an unmatched ORF in the middle: 2+ hits in same RF?!
                        start = end = (end + 1)
                    lastcov += 1
                    rfseq.info['Sequence'] = string.join(orfs[firstcov:lastcov],'*')
                    #self.deBug('%s' % orfs)
                    #self.deBug(rfseq.info['Sequence'])
                    
                ### ~ [3] ~ Select/Output the best (if there is one) or all if no best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if bestrf: myrfs.seq = [bestrf]
                for rfseq in myrfs.seqs(): open(fiestafas,'a').write('>%s\n%s\n' % (rfseq.info['Name'],rfseq.getSequence()))
                if bestrf and hitx:     # Best RF has BLAST hit
                    hit = searchdict[bestrf].hit[0]
                    open(fiestatdt,'a').write('%s\n' % string.join([est,'%d' % myrfs.seqNum(),bestrf.shortName(),hit.info['Name'],rje.expectString(hit.stat['E-Value'])],'\t'))
                else:
                    for rfseq in myrfs.seqs():
                        open(fiestatdt,'a').write('%s\n' % string.join([est,'%d' % myrfs.seqNum(),rfseq.shortName(),'-','-'],'\t'))

            return myrfs.seqNum()

        except KeyboardInterrupt:
            if rje.yesNo('Quit FIESTA?'): raise
            else: return 0
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadMapping(self,stype): ### Loads previous mappings from file
        '''Loads previous mappings from file.'''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '%s.%s.tdt' % (self.info['Basefile'],stype.lower())
            self.dict['Mapping'] = rje.dataDict(self,mfile,['Trans'],'All',getheaders=True,lists=False)
            if '-' in self.dict['Mapping']: self.dict['Mapping'].pop('-')
            self.list['Headers'] = self.dict['Mapping'].pop('Headers')
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def clusterTrans(self):     ### Read in fiestatdt to get shared hits and BLAST names for FASTACMD etc.
        '''Read in fiestatdt to get shared hits and BLAST names for FASTACMD etc.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fiestatdt = '%s.trans.tdt' % self.info['Basefile']
            searchdb = self.obj['BLAST'].info['DBase']
            ### ~ [2] ~ Read in Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#HITS','Compiling BLAST hit data...')
            self.loadMapping('Trans')
            clusters = rje.dataDict(self,fiestatdt,['TopHit'],['EST','Trans'],lists=True)
            if clusters.has_key('-'): nohits = clusters.pop('-')
            else: nohits = {'EST':[],'Trans':[]}
            nohits['EST'] = rje.sortUnique(nohits['EST'],xreplace=False)
            self.printLog('\r#HITS','%s different BLAST hits. %s ESTs (%s Trans) w/o hits.' % (rje.integerString(len(clusters)),rje.integerString(len(nohits['EST'])),rje.integerString(len(nohits['Trans']))))
            self.list['NoHits'] = nohits['Trans']
            while '-' in self.list['NoHits']: self.list['NoHits'].remove('-')
            self.dict['Clusters'] = clusters
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def annotateEST(self,seqlist,stype='trans'):    ### Uses BLAST etc. to translate/annotate EST sequences
        '''
        Uses BLAST etc. to translate/annotate EST sequences.
        >> seqlist:SeqList object => seqlist.info['Name'] defines input file.
        >> stype:str ['trans'] = defines output files (basefile.stype.fas and basefile.stype.tdt)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#EST','## ~~~~~~~ EST TRANSLATION/ANNOTATION ~~~~~~~ ##')
            if not stype: stype = 'ann'
            ## ~ [1a] ~ Output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fiestafas = '%s.%s.fas' % (self.info['Basefile'],stype.lower())
            fiestatdt = '%s.%s.tdt' % (self.info['Basefile'],stype.lower())
            append = self.opt['Append']
            ## ~ [1b] ~ Check for data from previous (aborted) run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # self.opt['Force'] will force re-running from scratch
            # self.opt['PickUp'] will look for existing data and then look to fill it in
            # if neither of the above and pickup data is found, a complete run will be assumed
            self.dict['PickUp'] = {}
            if rje.exists(fiestatdt) and not self.opt['Force']: 
                self.dict['PickUp'] = rje.dataDict(self,fiestatdt,['EST'],['Trans'],lists=True)
                self.printLog('#EST','%s EST sequences already processed (pickup=%s)' % (rje.integerString(len(self.dict['PickUp'])),str(self.opt['PickUp'])[0]))
                self.stat['ESTNum'] = len(self.dict['PickUp'])
            else:
                rje.backup(self,fiestafas)      
                rje.delimitedFileOutput(self,fiestatdt,['EST','RFX','Trans','TopHit','E'],'\t',rje_backup=True)
            ### ~ [2] Convert ESTs to RF(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['BLAST'].setStr({'DBase':self.info['SearchDB']})
            self.obj['BLAST'].setInt({'OneLine':1,'HitAln':1})          # Top Hit only (include more in future?)
            self.obj['BLAST'].formatDB(force=False)
            SEQFILE = open(seqlist.info['Name'], 'r')
            (nextseq,lastline) = (None,'Start')
            self.progLog('\r#EST','Counting EST sequences...')
            (sx,seqx,rfx) = (0,rje_seq.SeqCount(self,seqlist.info['Name']),0)
            self.stat['ESTNum'] = seqx
            self.printLog('\r#EST','%s EST sequences to process' % rje.integerString(seqx))
            while lastline:
                self.progLog('\r#EST','Processing ESTs: %.2f%%' % ((100.0 *sx)/seqx)); sx += 1
                ## ~ [2a] ~ Read next sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                seqlist.seq = [nextseq]     # Keep object small
                if nextseq.info['AccNum'] in self.dict['PickUp']: rfx += len(self.dict['PickUp'][nextseq.info['AccNum']]['Trans'])
                else: rfx += self.convertESTtoRF(nextseq,stype=stype)
                #!# Add option to truncate N-terminal if BLAST hit starts with M and included in alignment #!#
                #self.deBug(nextseq.info)
            SEQFILE.close()
            self.printLog('\r#EST','Processed ESTs: %s RF translations >> %s' % (rje.integerString(rfx),fiestafas))
            seqlist.info['Name'] = fiestafas
        except: self.errorLog('Problem during fiesta.annotate(%s)' % stype); return False
#########################################################################################################################
    def annotateTrans(self,seqlist,stype='final'):   ### Uses BLAST etc. to translate/annotate EST translations
        '''
        Uses BLAST etc. to annotate EST translations.
        >> seqlist:SeqList object => seqlist.info['Name'] defines input file.
        >> stype:str ['trans'] = defines output files (basefile.stype.fas and basefile.stype.tdt)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#TRANS','## ~~~~~~~ TRANSLATION ANNOTATION ~~~~~~~ ##')
            ## ~ [1a] ~ Output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fiestafas = '%s.%s.fas' % (self.info['Basefile'],stype)
            fiestatdt = '%s.%s.tdt' % (self.info['Basefile'],stype)
            rfblast = self.obj['BLAST']
            rfblast.setStr({'DBase':self.getStr('SearchDB'),'Type':'blastp'})
            rfblast.setInt({'OneLine':1,'HitAln':1})    # Top Hit only
            rfblast.setBool({'GappedBLAST':True})
            ## ~ [1b] ~ Check for data from previous (aborted) run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # self.opt['Force'] will force re-running from scratch
            # self.opt['PickUp'] will look for existing data and then look to fill it in
            # if neither of the above and pickup data is found, a complete run will be assumed
            self.dict['PickUp'] = {}
            if rje.exists(fiestatdt) and not self.opt['Force']: 
                self.dict['PickUp'] = rje.dataDict(self,fiestatdt,['SeqID'],['TopHit'],lists=True)
                self.printLog('#EST','%s EST translations already processed (pickup=%s)' % (rje.integerString(len(self.dict['PickUp'])),str(self.opt['PickUp'])[0]))
                self.stat['ESTNum'] = len(self.dict['PickUp'])
            else:
                rje.backup(self,fiestafas)      
                rje.delimitedFileOutput(self,fiestatdt,['SeqID','TopHit','E'],'\t',rje_backup=True)
            ### ~ [2] BLAST Annotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['PickUp'] or not self.dict['PickUp'] or self.opt['Force']:
                rfblast.formatDB(force=False)
                SEQFILE = open(seqlist.info['Name'], 'r')
                (nextseq,lastline) = (None,'Start')
                self.progLog('\r#SEQ','Counting sequences...')
                (sx,seqx,bx) = (0,rje_seq.SeqCount(self,seqlist.info['Name']),0)
                self.printLog('\r#SEQ','%s sequences to process' % rje.integerString(seqx))
                while lastline:
                    self.progLog('\r#SEQ','Processing sequences: %.2f%%' % (100.0*sx/seqx)); sx += 1
                    ## ~ [2a] ~ Read next sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                    if not nextseq: self.debug('Done?'); continue
                    seqlist.seq = [nextseq]     # Keep object small
                    if nextseq.shortName() in self.dict['PickUp']:
                        bx += len(self.dict['PickUp'][nextseq.shortName()]['TopHit'])
                        continue
                    ## ~ [2b]~ BLAST sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    qfile = '%s.%s.tmp.fas' % (nextseq.shortName(),rje.randomString(6))
                    seqlist.saveFasta(seqfile=qfile,log=False)
                    rfblast.setStr({'InFile':qfile,'Name':'%s.blast' % rje.baseFile(qfile)})
                    rfblast.blast(log=False); os.unlink(qfile)
                    #try: rfblast.readBLAST(clear=True,screen=self.debugging(),unlink=not self.debugging(),local=True,gablam=False)
                    try: rfblast.readBLAST(clear=True,screen=False,unlink=True,local=True,gablam=False)
                    except: self.errorLog('%s problem.' % rfblast.getStr('Name'))
                    ## ~ [2c] ~ Read in alignment and map onto sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if not rfblast.getBool('OldBLAST'):
                        query = nextseq.shortName()
                        try:
                            if rfblast.queryHits(query):
                                hit = rfblast.queryHits(query)[0]
                                hdict = rfblast.hitData(hit=hit,query=query)
                                nextseq.info['Name'] = '%s Similar (e=%s) to %s %s' % (nextseq.shortName(),rje_blast.expectString(hdict['E-Value']),hit,hdict['Description'])
                                #!# Add option to truncate N-terminal if BLAST hit starts with M and included in alignment #!#
                                #!# This needs improvement. What if there is more than one M in the query and/or STOPs? #!#
                                if self.opt['TruncNT']:
                                    try:
                                        hit.globalFromLocal(search.stat['Length'],keepaln=True)
                                        qstart = hit.dict['GABLAM']['Aln']['QryStartO'] - hit.dict['GABLAM']['Aln']['SbjStartO']
                                        #if qstart > 0:
                                            #for i in range(qstart,0,-1):
                                            #    if nextseq.info['Sequence'][i-1:i] == 'M': nextseq.info['Sequence'][i-1:i] = nextseq.info['Sequence'][i-1:]; break
                                        mstart = nextseq.getSequence().find('M')
                                        if mstart <= qstart and mstart > 0: nextseq.info['Sequence'] = nextseq.getSequence()[mstart:]
                                    except: pass
                            else:
                                hit = None
                                nextseq.info['Name'] = '%s No BLAST hit (e<%s) to %s' % (nextseq.shortName(),rje_blast.expectString(rfblast.getNum('E-Value')),os.path.basename(rfblast.getStr('DBase')))
                        except:
                            if self.dev() or self.debugging(): self.errorLog('Arse'); raise
                            nextseq.info['Name'] = '%s (BLAST error)' % nextseq.shortName()
                        ## ~ [2d] ~ Select/Output the best (if there is one) or all if no best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                        open(fiestafas,'a').write('>%s\n%s\n' % (nextseq.info['Name'],nextseq.getSequence()))
                        if hit:
                            hdict = rfblast.hitData(hit=hit,query=query)
                            open(fiestatdt,'a').write('%s\n' % string.join([nextseq.shortName(),hit,rje_blast.expectString(hdict['E-Value'])],'\t')); bx += 1
                        else: open(fiestatdt,'a').write('%s\n' % string.join([nextseq.shortName(),'-','-'],'\t'))
                    else:
                        hit = None
                        try:
                            search = rfblast.search[0]
                            if search.hit:
                                hit = search.hit[0]
                                try:
                                    hitseq = seqlist.seqFromFastaCmd(hit.info['Name'],rfblast.info['DBase'])
                                    nextseq.info['Name'] = '%s Similar (e=%s) to %s' % (nextseq.shortName(),rje_blast.expectString(hit.stat['E-Value']),hitseq.info['Name'])
                                except: nextseq.info['Name'] = '%s Similar (e=%s) to %s' % (nextseq.shortName(),rje_blast.expectString(hit.stat['E-Value']),hit.info['Name'])
                                #!# Add option to truncate N-terminal if BLAST hit starts with M and included in alignment #!#
                                #!# This needs improvement. What if there is more than one M in the query and/or STOPs? #!#
                                if self.opt['TruncNT']:
                                    try:
                                        hit.globalFromLocal(search.stat['Length'],keepaln=True)
                                        qstart = hit.dict['GABLAM']['Aln']['QryStartO'] - hit.dict['GABLAM']['Aln']['SbjStartO']
                                        #if qstart > 0:
                                            #for i in range(qstart,0,-1):
                                            #    if nextseq.info['Sequence'][i-1:i] == 'M': nextseq.info['Sequence'][i-1:i] = nextseq.info['Sequence'][i-1:]; break
                                        mstart = nextseq.getSequence().find('M')
                                        if mstart <= qstart and mstart > 0: nextseq.info['Sequence'] = nextseq.getSequence()[mstart:]
                                    except: pass
                            else: nextseq.info['Name'] = '%s No BLAST hit (e<%s) to %s' % (nextseq.shortName(),rje_blast.expectString(rfblast.stat['E-Value']),os.path.basename(rfblast.info['DBase']))
                        except: nextseq.info['Name'] = '%s (BLAST error)' % nextseq.shortName()
                        ## ~ [2d] ~ Select/Output the best (if there is one) or all if no best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                        open(fiestafas,'a').write('>%s\n%s\n' % (nextseq.info['Name'],nextseq.getSequence()))
                        if hit: open(fiestatdt,'a').write('%s\n' % string.join([nextseq.shortName(),hit.info['Name'],rje_blast.expectString(hit.stat['E-Value'])],'\t')); bx += 1
                        else: open(fiestatdt,'a').write('%s\n' % string.join([nextseq.shortName(),'-','-'],'\t'))
                SEQFILE.close()
                self.printLog('\r#SEQ','Processed sequences: %s with BLAST hits >> %s' % (rje.integerString(bx),fiestafas))
        except: self.errorLog('Problem during fiesta.annotateTrans(%s)' % stype); return False
#########################################################################################################################
    ### <4> ### Consensus Methods                                                                                       #
#########################################################################################################################
    def makeConsensi(self,seqlist,log='',stype='CON',revcomp=False,oneqry=None): ### Combines sequences in seqlist into consensus sequences where possible
        '''
        Combines sequences in seqlist into consensus sequences where possible.
        >> seqlist:rje_seq.SeqList object - modified in situ
        >> returns True if changes made, else False if not
        >> log:bool [False] = whether to log progress
        >> stype:str ['CON'] = consensus stype for renaming and mapping
        >> revcomp:bool [False] = whether to also try reverse complemented sequences for DNA consensi
        >> oneqry:Sequence [None] = selected query for single (partial) assembly 
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            revcomp = revcomp and self.opt['DNA']
            seqlist.sortByLen(proglog=log and not oneqry)                  # Now ordered with longest sequence first
            inseq = seqlist.seq[0:]         # List of sequence objects to manipulate
            newseq = []                     # List of sequence objects converted to consensi
            if self.stat['MinID'] > 1:  self.stat['MinID'] /= 100.0
            minid = max(0,1.0 - self.stat['MinID'])                 # Min non-identity of shared region for mismatches
            if self.opt['DNA']: minaln = self.stat['MinAln']        # Min length of shared region 
            else: minaln = min(self.stat['MinORF'],self.stat['MinAln']) 
            maxm = 0                                                # Max number of mismatches - will scale with overlap
            if stype not in self.dict['Counter']: self.dict['Counter'][stype] = 0
            ### ~ [2] ~ Compile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = len(inseq)
            ctxt = '%s Consensi' % log
            while inseq:    # If inseq contains sequence, more consensus clusters must be made
                if log: self.progLog('\r#CONS','%s: %.2f%% (%s->%s) |>>>>>|<|     ' % (ctxt,100.0*(prex-len(inseq))/prex,rje.integerString(prex-len(inseq)),rje.integerString(len(newseq))))
                if oneqry:
                    inseq.remove(oneqry); qry = oneqry
                    subseq = qry.list['SubSeq']          # Use this to store a list of sequences to make final consensus
                else:
                    qry = inseq.pop(0)
                    #self.deBug(qry.shortName())
                    #if not self.opt['DNA'] and qry.nonX() < int(self.stat['MinID']*self.stat['MinORF']): inseq = []; break    # Too short - ordered by length so all other too
                    #elif self.opt['DNA'] and qry.nonN() < int(self.stat['MinID']*self.stat['MinAln']): inseq = []; break      # Too short - ordered by length so all other too
                    subseq = [qry.getSequence()] # Use this to store a list of sequences to make final consensus
                if 'Consensus' not in qry.list: qry.list['Consensus'] = [qry.shortName()]
                refseq = qry.getSequence()   # Use this to compile consensus reference sequence
                concycle = True; checkrev = revcomp; cycx = 0
                self.vPrint('\n%s' % refseq,v=2)
                while concycle:                 # More subsequences added - must now check again against earlier seqs
                    checkrev = revcomp and not checkrev     # Check reverse comp every other cycle
                    if not checkrev: cycx = len(subseq); concycle = False       # First (possibly only) loop of this cycle
                    (sx,stot) = (0.0,len(inseq))
                    for seq in inseq[0:]:
                        if log: self.progLog('\r#CONS','%s: %.2f%% (%s->%s) |%.1f%%|%d|   ' % (ctxt,100.0*(prex-len(inseq))/prex,rje.integerString(prex-len(inseq)),rje.integerString(len(newseq)+1),sx/stot,len(subseq))); sx += 100.0
                        qseq = refseq; qlen = len(qseq)
                        sseq = seq.getSequence(); slen = len(sseq)
                        if checkrev: sseq = rje_sequence.reverseComplement(sseq)
                        #self.deBug('%s:%d vs %d' % (seq.shortName(),seq.nonX(),int(self.stat['MinID']*self.stat['MinORF'])))
                        ## ~ [2a] ~ Check for complete redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if sseq in refseq:
                            qry.list['Consensus'].append(seq.shortName())
                            si = refseq.find(sseq)
                            subseq.append('%s%s' % ('-' * si,sseq))
                            inseq.remove(seq)
                            self.vPrint('\n%s + %d' % (qry.list['Consensus'],len(inseq)),v=2)
                            continue
                        if not self.opt['DNA'] and seq.nonX() < int(self.stat['MinID']*self.stat['MinORF']): break   # Too short - ordered by length so all other too
                        elif self.opt['DNA'] and seq.nonN() < int(self.stat['MinID']*self.stat['MinAln']): break
                        ## ~ [2b] ~ Find partial overlap and combine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        alnx = 0; bestoff = 0; bestmm = qlen
                        offset = -(qlen - minaln + 1)
                        while alnx < qlen and offset < (slen - minaln +1):
                            offset += 1     # Move along one
                            ix = 0          # Identical residue count
                            mm = 0          # Mismatch count
                            if offset < 0: qi = -offset; si = 0
                            else: qi = 0; si = offset
                            overlap = min(qlen-qi,slen-si)
                            if overlap < alnx: continue
                            if minid: maxm = int(minid * overlap) + 1     # Max no. mismatches (round up)
                            for i in range(overlap):
                                if qseq[qi+i] == sseq[si+i]: ix += 1
                                else: mm += 1
                                if mm > maxm or (overlap - mm) < alnx: break
                            if mm > maxm or ix < alnx: continue         # Not good enough
                            if ix == alnx and mm >= bestmm: continue    # No better than existing
                            alnx = ix; bestoff = offset; bestmm = mm    # Update best alignment (and good enough)
                            bestmaxm = maxm; bestolap = overlap
                        ## ~ [2c] ~ Combine if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if alnx < 1: continue                           # Do not combine
                        #!# Improve to handle mismatches better #!#
                        if bestoff < 0:
                            self.vPrint('\n%dMM vs %.2f x %d = %d' % (bestmm,minid,min(qlen+bestoff,slen),bestmaxm),v=2)
                            self.vPrint('\n%s' % qseq,v=2)
                            self.vPrint('%s%s' % ('-' * -bestoff,sseq),v=2)
                            refseq = qseq + sseq[qlen+bestoff:]
                            subseq.append('%s%s' % ('-' * -bestoff,sseq))
                        else:
                            self.vPrint('\n%dMM vs %.2f x %d = %d' % (bestmm,minid,min(slen-bestoff,qlen),bestmaxm),v=2)
                            self.vPrint('\n%s%s' % ('-' * bestoff,qseq),v=2)
                            self.vPrint(sseq,v=2)
                            refseq = sseq[:bestoff] + qseq
                            for i in range(len(subseq)): subseq[i] = '%s%s' % ('-' * bestoff,subseq[i])
                            subseq.append(sseq)
                        self.vPrint(refseq,v=2)
                        if bestoff >= 0 and self.v() >= 2: self.deBug('')
                        qry.list['Consensus'].append(seq.shortName())
                        inseq.remove(seq)
                        self.vPrint('%s + %d' % (qry.list['Consensus'],len(inseq)),v=2)
                    ## ~ [2d] ~ Assess next cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    concycle = cycx < len(subseq)          # Checkpoint for whether cycle should continue
                    if revcomp and not checkrev: concycle = True    # Loop to check revcomp
                ## ~ [2e] ~ Make consensus and append sequence to new list ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if oneqry:
                    self.vPrint('%d -> %d + %d' % (seqlist.seqNum(),len(subseq),len(inseq)),v=2)
                    prex = seqlist.seqNum() 
                    seqlist.seq = [oneqry] + inseq[0:]
                    oneqry.info['Sequence'] = refseq
                    oneqry.list['SubSeq'] = subseq
                    return prex != seqlist.seqNum()
                newseq.append(qry)
                qry.info['Sequence'] = self.makeConsensus(refseq,subseq)
            if log and log != 'GABLAM': self.printLog('\r#CONS','%s: %s sequences -> %s consensi.   ' % (ctxt,rje.integerString(prex),rje.integerString(len(newseq))))
            ### ~ [3] ~ Finish and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,len(newseq))
            for seq in newseq:
                if oneqry: raise ValueError
                if log: self.progLog('\r#CONS','%s renaming: %.2f%%' % (ctxt,sx/stot)); sx += 100.0
                self.dict['Counter'][stype] += 1
                seq.info['AccNum'] = '%s%s%s' % (self.info['NewAcc'],stype,rje.preZero(self.dict['Counter'][stype],self.stat['ESTNum']))
                seq.info['ID'] = '%s_%s' % (stype.lower(),self.info['SpCode'])
                seq.info['Name'] = '%s %s' % (seq.info['AccNum'],seq.info['Description'])
                if self.opt['GnSpAcc']: seq.info['Name'] = '%s__%s' % (seq.info['ID'],seq.info['Name'])
                for oldid in seq.list['Consensus']:
                    if oldid not in self.dict['Mapping']: self.dict['Mapping'][oldid] = {}
                    self.dict['Mapping'][oldid][stype] = seq.info['AccNum']
                    self.vPrint('%s -> %s' % (oldid,self.dict['Mapping'][oldid]),v=2)
            if log and log != 'GABLAM': self.printLog('\r#CONS','%s renaming complete.' % ctxt,log=False)
            seqlist.seq = newseq[0:]
            if seqlist.seqNum() == stot: return False       # No reduction in sequence number
            return True
        except: self.errorLog(rje_zen.Zen().wisdom()); raise
#########################################################################################################################
    def makeConsensus(self,refseq,subseq):  ### Makes a consensus sequence from reference and subsequences
        '''
        Makes a consensus sequence from reference and subsequences.
        >> refseq:str = Full length Reference sequence
        >> subseq:list of strings making consensus (may be missing Cterm/3')
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(subseq) == 1: return refseq
            conseq = ''
            wild = {True:'N',False:'X'}[self.opt['DNA']]
            ### ~ [1] ~ Compile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(refseq)):
                r = refseq[i]
                ilist = []
                for seq in subseq:
                    if seq[i:i+1] not in ['-','*',wild,'']: ilist.append(seq[i])
                for x in ilist:
                    if ilist.count(x) > ilist.count(r): r = x
                conseq = conseq + r
            return conseq
        except: self.errorLog(rje_zen.Zen().wisdom()); return refseq
#########################################################################################################################
    def gablamConsensi(self,seqlist,infile,outfile,stype='CON'):   ### Iterative GABLAM and consensus generation
        '''Iterative GABLAM and consensus generation.'''
        try:### ~ [1] ~  Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            processed = []      # List of sequences that have been processed
            seqdict = seqlist.seqNameDic(key='short',proglog=True)
            inseq = seqlist.seq[0:]
            gblast = self.obj['BLAST']
            gblast.info['DBase'] = infile
            gblast.stat['OneLine'] = gblast.stat['HitAln'] = seqlist.seqNum()
            gblast.opt['GappedBLAST'] = self.opt['GapBLAST']
            if self.opt['DNA']: gblast.info['stype'] = 'blastn'
            else: gblast.info['Type'] = 'blastp'
            #self.deBug(gblast.info)
            gblast.formatDB(force=False,protein=not self.opt['DNA'])
            rje.backup(self,outfile)
            resavefile = '%s.resave.fas' % rje.baseFile(infile)
            savex = prex = seqlist.seqNum()
            conx = 0
            ### ~ [2] ~ Loop, GABLAM and Make Consensus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while inseq:
                glead = '%s GABLAM consensi (%s' % (infile,rje.integerString(prex-len(inseq)))
                self.progLog('\r#GCON','%s->%s): %.2f%% ||     ' % (glead,rje.integerString(conx),100.0*(prex-len(inseq))/prex))
                seq = inseq.pop(0)
                if seq in processed: continue
                ## ~ [2a] ~ Use GABLAM to make a cluster of sequence to make consensi from ~~~~~~~~ ##
                gx = 0; glist = [seq]   # List of sequences to GABLAM and expand GABLAM cluster
                while gx < len(glist):  # glist growing
                    ## ~ [2a-i] ~ Check for sequence reduction via resave ~~~~~~~~~~~~~~~ ##
                    if savex - len(inseq) >= self.stat['ReSave']:    ### Resave & format search data
                        self.printLog('#RESAVE',resavefile,log=False)
                        seqlist.saveFasta(inseq,seqfile=resavefile,log=False,screen=True)
                        gblast.info['DBase'] = resavefile
                        savex = gblast.stat['OneLine'] = gblast.stat['HitAln'] = len(inseq)
                        gblast.formatDB(force=True,protein=not self.opt['DNA'],log=False)
                    gx = len(glist)
                    self.progLog('\r#GCON','%s->%s): %.2f%% |%d|     ' % (glead,rje.integerString(conx),100.0*(prex-len(inseq))/prex,len(glist)))
                    ## ~ [2a-ii] ~ GABLAM each sequence in cluster unless alread done ~~~~ ##
                    for qry in glist[0:]:
                        if qry in processed: continue   # Done in a previous cycle
                        processed.append(qry)
                        qfile = '%s.%s.tmp.fas' % (qry.info['AccNum'],rje.randomString(6))
                        seqlist.saveFasta([qry],seqfile=qfile,log=False)
                        gblast.info['InFile'] = qfile
                        gblast.info['Name'] = '%s.blast' % rje.baseFile(qfile)
                        self.progLog('\r#GCON','%s~>%s): %.2f%% |%d|     ' % (glead,rje.integerString(conx),100.0*(prex-len(inseq))/prex,len(glist)))
                        gblast.blast(); os.unlink(qfile)
                        self.progLog('\r#GCON','%s=>%s): %.2f%% |%d|     ' % (glead,rje.integerString(conx),100.0*(prex-len(inseq))/prex,len(glist)))
                        try: gblast.readBLAST(clear=True,screen=self.opt['DeBug'],unlink=True,gablam=False)
                        except: self.errorLog('%s problem.' % gblast.info['Name'])
                        try: search = gblast.search[0]
                        except: continue
                        for hit in search.hit:
                            try: hseq = seqdict[hit.info['Name']]
                            except: self.errorLog('Problem with hit "%s" - not found in seqdict' % hit.info['Name']); continue
                            if hseq in processed: continue
                            hit.globalFromLocal(search.stat['Length'])
                            glen = hit.dict['GABLAM']['Query']['GABLAMO Len']
                            gid = hit.dict['GABLAM']['Query']['GABLAMO ID']
                            if glen < self.stat['MinAln'] or (float(gid) / glen) < self.stat['MinID']: continue
                            if hseq not in glist: glist.append(hseq)
                            if hseq in inseq: inseq.remove(hseq)
                ## ~ [2b] ~ Generate consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#GCON','%s#>%s): %.2f%% |%d|     ' % (glead,rje.integerString(conx),100.0*(prex-len(inseq))/prex,len(glist)))
                clog = {True:'GABLAM',False:''}[self.opt['DeBug'] or self.stat['Verbose'] >= 1]
                seqlist.seq = glist
                self.makeConsensi(seqlist,stype=stype,revcomp=self.opt['DNA'],log=clog)
                seqlist.saveFasta(seqfile=outfile,append=True,log=False)
                conx += seqlist.seqNum()
            self.printLog('\r#GCON','%s %s sequences -> %s %s GABLAM consensi.' % (rje.integerString(prex),infile,rje.integerString(conx),stype))
            if os.path.exists(resavefile): rje_blast.cleanupDB(self,resavefile,deletesource=True)
        except:  self.errorLog('Error during FIESTA GABLAM consensi',quitchoice=True)
#########################################################################################################################
    def gablamConsensiOneQry(self,seqlist,infile,outfile,stype='CON'):   ### Iterative GABLAM and consensus generation
        '''Iterative GABLAM and consensus generation.'''
        try:### ~ [1] ~  Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            processed = []      # List of sequences that have been processed
            seqlist.sortByLen(proglog=True)                  # Now ordered with longest sequence first
            seqdict = seqlist.seqNameDic(key='short',proglog=True)
            inseq = seqlist.seq[0:]
            gblast = self.obj['BLAST']
            gblast.info['DBase'] = infile
            gblast.stat['OneLine'] = gblast.stat['HitAln'] = seqlist.seqNum()
            gblast.opt['GappedBLAST'] = self.opt['GapBLAST']
            if self.opt['DNA'] and self.opt['Test']: gblast.info['BLASTOpt'] = '-n T'
            if self.opt['DNA']: gblast.info['Type'] = 'blastn'
            else: gblast.info['Type'] = 'blastp'
            self.bugPrint(gblast.info)
            gblast.formatDB(force=False,protein=not self.opt['DNA'])
            rje.backup(self,outfile)
            resavefile = '%s.resave.fas' % rje.baseFile(infile)
            savex = prex = seqlist.seqNum()
            conx = 0; procx = 0; procseq = []
            if stype not in self.dict['Counter']: self.dict['Counter'][stype] = 0
            ### ~ [2] ~ Loop, GABLAM and Make Consensus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while inseq:    ### Keep cycling to make consensi
                processed = procseq[0:]
                glead = '%s GABLAM consensi (%s|%s' % (infile,rje.integerString(procx),rje.integerString(prex-len(inseq)))
                self.progLog('\r#GCON','%s->%s): %.2f%% ||     ' % (glead,rje.integerString(conx),100.0*(procx)/prex))
                #self.bugPrint('\n%d inseq; %d processed; %d consensi; next: %s' % (len(inseq),len(processed),conx,inseq[0].shortName()))
                ## ~ [2a] ~ Check for sequence reduction via resave ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if savex - len(inseq) >= self.stat['ReSave']:    ### Resave & format search data
                    self.printLog('\n#RESAVE','%s -> %s' % (rje.integerString(len(inseq)),resavefile),log=False)
                    seqlist.saveFasta(inseq,seqfile=resavefile,log=False,screen=False)
                    gblast.info['DBase'] = resavefile
                    savex = gblast.stat['OneLine'] = gblast.stat['HitAln'] = len(inseq)
                    gblast.formatDB(force=True,protein=not self.opt['DNA'],log=False)
                ## ~ [2b] ~ Setup Query for next oneqry consensus-making ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qry = seq = inseq.pop(0)      # Sequence to act as onequery 
                #self.deBug('%s >> %s [%s]' % (processed,seq,seq in processed))
                if seq in processed: continue
                seq.list['SubSeq'] = [seq.getSequence()]
                seq.list['Consensus'] = [seq.shortName()]
                ## ~ [2b] ~ Use GABLAM to make a cluster of sequence to make consensi from ~~~~~~~~ ##
                glist = [seq]
                while qry:  # Perform GABLAM, make consensus and cycles while consensus changed
                    self.bugPrint('\n%d inseq; %d processed; %d consensi; query=%s; %d subseq; %d GList' % (len(inseq),len(processed),conx,seq.shortName(),len(seq.list['SubSeq']),len(glist)))
                    gx = len(glist)
                    self.progLog('\r#GCON','%s->%s): %.2f%% |%d:%d|     ' % (glead,rje.integerString(conx),100.0*(procx+len(seq.list['SubSeq']))/prex,len(glist),len(seq.list['SubSeq'])))
                    ## ~ [2b-i] ~ GABLAM each sequence in cluster unless alread done ~~~~ ##
                    processed.append(qry)
                    qfile = '%s.%s.tmp.fas' % (qry.info['AccNum'],rje.randomString(6))
                    seqlist.saveFasta([qry],seqfile=qfile,log=False)
                    gblast.info['InFile'] = qfile
                    gblast.info['Name'] = '%s.blast' % rje.baseFile(qfile)
                    self.progLog('\r#GCON','%s~>%s): %.2f%% |%d:%d|     ' % (glead,rje.integerString(conx),100.0*(procx+len(seq.list['SubSeq']))/prex,len(glist),len(seq.list['SubSeq'])))
                    gblast.blast(); os.unlink(qfile)
                    self.progLog('\r#GCON','%s=>%s): %.2f%% |%d:%d|     ' % (glead,rje.integerString(conx),100.0*(procx+len(seq.list['SubSeq']))/prex,len(glist),len(seq.list['SubSeq'])))
                    try: gblast.readBLAST(clear=True,screen=self.opt['DeBug'],unlink=True,gablam=False)
                    except: self.errorLog('%s problem.' % gblast.info['Name'])
                    try: search = gblast.search[0]
                    except: seqlist.seq = [qry]; break
                    for hit in search.hit:
                        if hit.info['Name'] in seq.list['Consensus']: continue     # Already added
                        try: hseq = seqdict[hit.info['Name']]
                        except: self.errorLog('Problem with hit "%s" - not found in seqdict' % hit.info['Name']); continue
                        if hseq in processed: continue
                        hit.globalFromLocal(search.stat['Length'])
                        glen = hit.dict['GABLAM']['Query']['GABLAMO Len']
                        gid = hit.dict['GABLAM']['Query']['GABLAMO ID']
                        if glen < self.stat['MinAln'] or (float(gid) / glen) < self.stat['MinID']: continue
                        if hseq not in glist: glist.append(hseq)
                        ## ~ Special DNA manipulation ~ ##
                        if self.opt['GabRev'] and self.opt['DNA']:
                            if hit.dict['GABLAM']['Query']['GABLAMO Dirn'] != hit.dict['GABLAM']['Hit']['GABLAMO Dirn']:   # Whether to reverse complement
                                if hseq.opt['RevComp'] == qry.opt['RevComp']: hseq.reverseComplement()
                            elif hseq.opt['RevComp'] != qry.opt['RevComp']: hseq.reverseComplement()
                    ## ~ [2b-ii] ~ Cycle or assess query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    seqlist.seq = glist[0:]         # Qry + GABLAM hits (squared)
                    self.bugPrint('\n%d inseq; %d processed; %d consensi; query=%s; %d subseq; %d GList; blastqry=%s' % (len(inseq),len(processed),conx,seq.shortName(),len(seq.list['SubSeq']),len(glist),qry.shortName()))
                    if len(glist) == 1: break   # No GABLAM homologues found
                    self.progLog('\r#GCON','%s#>%s): %.2f%% |%d:%d|     ' % (glead,rje.integerString(conx),100.0*(procx+len(seq.list['SubSeq']))/prex,len(glist),len(seq.list['SubSeq'])))
                    self.bugPrint('\nMake consensi from sequences so far...')
                    prelen = seq.seqLen(); clog = {True:'Query',False:''}[self.opt['DeBug'] or self.stat['Verbose'] >= 1]
                    qmod = self.makeConsensi(seqlist,stype=stype,revcomp=self.opt['DNA'] and not self.opt['GabRev'],oneqry=seq,log=clog)    # Query modified enough to re-BLAST
                    self.bugPrint('Query modified = %s' % qmod)
                    self.bugPrint('\n%d inseq; %d processed; %d consensi; query=%s; %d subseq; %d in -> %d out' % (len(inseq),len(processed),conx,seq.shortName(),len(seq.list['SubSeq']),len(glist),seqlist.seqNum()))
                    if qmod:
                        for gseq in glist:
                            if gseq not in seqlist.seq and gseq in inseq: inseq.remove(gseq)    # Assimilated into Consensus
                        for hseq in seqlist.seq:
                            if hseq in processed: processed.remove(hseq)    # Need to process all remaining again, including query
                        glist = [seq]
                        self.bugPrint('=> %d inseq; %d processed; %d consensi; query=%s; %d subseq' % (len(inseq),len(processed),conx,seq.shortName(),len(seq.list['SubSeq'])))
                        if prelen == seq.seqLen(): break    # No need to Re-BLAST (no increase in length)
                    else: break
                ## ~ [2c] ~ Compile finished Query consensus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                glead = '%s GABLAM consensi (%s|%s' % (infile,rje.integerString(procx),rje.integerString(prex-len(inseq)))
                self.progLog('\r#GCON','%s->%s): %.2f%% |%d:%d|     ' % (glead,rje.integerString(conx),100.0*(procx+len(seq.list['SubSeq']))/prex,len(glist),len(seq.list['SubSeq'])))
                self.bugPrint('=> Final consensus ... %d inseq; %d processed; %d consensi; query=%s; %d subseq' % (len(inseq),len(processed),conx,seq.shortName(),len(seq.list['SubSeq'])))
                procx += len(seq.list['SubSeq'])
                seq.info['Sequence'] = self.makeConsensus(seq.getSequence(),seq.list['SubSeq'])
                self.dict['Counter'][stype] += 1
                seq.info['AccNum'] = '%s%s%s' % (self.info['NewAcc'],stype,rje.preZero(self.dict['Counter'][stype],self.stat['ESTNum']))
                seq.info['ID'] = '%s_%s' % (stype.lower(),self.info['SpCode'])
                seq.info['Name'] = '%s %s' % (seq.info['AccNum'],seq.info['Description'])
                if self.opt['GnSpAcc']: seq.info['Name'] = '%s__%s' % (seq.info['ID'],seq.info['Name'])
                self.bugPrint(seq.info['Name'])
                self.bugPrint(seq.list['Consensus'])
                for oldid in seq.list['Consensus']:
                    procseq.append(seqdict[oldid])
                    if oldid not in self.dict['Mapping']: self.dict['Mapping'][oldid] = {}
                    self.dict['Mapping'][oldid][stype] = seq.info['AccNum']
                    self.bugPrint('%s -> %s' % (oldid,self.dict['Mapping'][oldid]))
                #for cseq in seqlist.seq[1:]:
                #    if cseq in inseq: self.deBug('%s already in inseq?!' % cseq.shortName())
                #    elif cseq in procseq: self.deBug('%s already in procseq?!' % cseq.shortName())
                seqlist.saveFasta([seq],seqfile=outfile,append=True,log=False)
                conx += 1
                self.bugPrint('Total = %d' % (len(inseq)+procx))
            self.printLog('\r#GCON','%s %s sequences -> %s %s GABLAM consensi.' % (rje.integerString(prex),infile,rje.integerString(conx),stype))
            if os.path.exists(resavefile): rje_blast.cleanupDB(self,resavefile,deletesource=True)
        except:  self.errorLog('Error during FIESTA GABLAM consensi (OneQry)',quitchoice=True)
#########################################################################################################################
    def saveMapping(self,stype):     ### Loads in previous TDT file(s) and outputs mapping
        '''
        Loads in previous TDT file(s) and outputs mapping.
        >> stype:str = Output consensus stype
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: prevtypes = {'TOP':['Trans'],'HIT':['TOP'],'NON':['Trans'],'CON':['HIT','NON']}[stype]
            except: prevtypes = []
            next = stype.lower()
            if stype not in self.list['Headers']: self.list['Headers'].append(stype)
            ### ~ [1] ~ Add current data and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.%s.tdt' % (self.info['Basefile'],next)
            rje.delimitedFileOutput(self,outfile,self.list['Headers'],rje_backup=True)
            (tx,px) = (0.0,len(self.dict['Mapping']))
            for trans in rje.sortKeys(self.dict['Mapping']):
                self.progLog('\r#MAP','Mapping %s data: %.2f%%' % (stype,tx/px)); tx += 100.0
                if 'Trans' not in self.dict['Mapping'][trans] or self.dict['Mapping'][trans]['Trans'] != trans: continue
                ## ~ [1a] ~ Map to mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for prev in prevtypes:
                    if prev in self.dict['Mapping'][trans]:
                        oldid = self.dict['Mapping'][trans][prev]
                        try: self.dict['Mapping'][trans][stype] = self.dict['Mapping'][oldid][stype]
                        except: pass #X#self.errorLog('%s-%s missing from mapping dictionary' % (oldid, stype),printerror=False)
                ## ~ [1b] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rje.delimitedFileOutput(self,outfile,self.list['Headers'],datadict=self.dict['Mapping'][trans])
            self.printLog('\r#MAP','Mapped %s data output to %s' % (stype,outfile))
        except: self.errorLog('Problem with fiesta.saveMapping(%s)' % stype)
#########################################################################################################################
    ### <5> ### EST Assembly Methods                                                                                    #
#########################################################################################################################
    def fiestaDNA(self):  ### Main run method.
        '''
        Main run method. This pipeline is based on the TIGR rules: 95%+ identity over 40+nt with <25nt unmatched overlap
        at each end. Each EST will be taken in turn and BLASTed. Any hits will then be BLASTed first to build up
        clusters. Each BLAST cluster is then processed to make actual UniGene-style clusters and alignments. The next EST
        is then used for BLAST and the process continues. After a certain number of ESTs are clustered and removed, the
        remains are resaved to speed up subsequent BLASTs.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F','accnr=F','seqnr=F','replacechar=F'])
            if not rje.exists(seqlist.info['Name']): return self.errorLog('Input file "%s" not found' % seqlist.info['Name'],printerror=False)
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = '%s.fiesta' % rje.baseFile(seqlist.info['Name'])
            self.obj['BLAST'] = rje_blast.BLASTRun(self.log,self.cmd_list)
            self.obj['BLAST'].stat['Verbose'] = self.stat['Verbose'] - 1
            self.list['Headers'] = ['EST','CONS']
            ## ~ [1a] ~ Sequence file and EST sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            estfile = '%s.est.fas' % self.info['Basefile']
            if os.path.exists(estfile) and not self.opt['Force']: seqlist.loadSeqs(estfile,filetype='fas',seqtype='dna',aln=False,nodup=False)
            else: seqlist.loadSeqs(filetype='fas',seqtype='dna',aln=False,nodup=False)
            (sx,seqx) = (0.0,seqlist.seqNum())
            self.stat['ESTNum'] = seqx
            self.dict['Mapping'] = {}
            shortx = 0
            for seq in seqlist.seq[0:]:
                self.progLog('\r#MAP','Setting up EST->Consensus mapping: %.1f%%' % (sx/seqx)); sx += 100.0
                if seq.nonN() < int(self.stat['MinID']*self.stat['MinAln']): shortx += 1; seqlist.seq.remove(seq); continue
                seq.info['Name'] = '%s %s' % (seq.info['AccNum'],seq.info['Description'])
                self.dict['Mapping'][seq.shortName()] = {'EST':seq.shortName()}
            self.printLog('\r#MAP','Setting up EST->Consensus mapping complete.',log=False)
            if shortx: self.printLog('#REM','%s EST sequences removed - too few non-N nt' % rje.integerString(shortx))
            seqlist.info['Name'] = estfile
            if not os.path.exists(estfile) or self.opt['Force'] or shortx:
                rje.backup(self,estfile)
                seqlist.saveFasta()
            self.printLog('\r#EST','%s EST sequences to process' % rje.integerString(seqlist.seqNum()))
            ### ~ [2] GABLAM-based DNA consensus generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fasout = '%s.cons.fas' % self.info['Basefile']
            tdtout = '%s.cons.tdt' % self.info['Basefile']
            if self.opt['Force']: rje.backup(self,fasout); rje.backup(self,tdtout)
            if not os.path.exists(fasout) or self.opt['Force']:
                if self.info['AssMode'].lower() in ['query','oneqry','onequery']:
                    self.gablamConsensiOneQry(seqlist,seqlist.info['Name'],fasout,stype='CONS')
                elif self.info['AssMode'].lower() in ['gab','gablam']:
                    self.gablamConsensi(seqlist,seqlist.info['Name'],fasout,stype='CONS')
                else:
                    self.makeConsensi(seqlist,log='FiestaDNA',stype='CONS',revcomp=True)
                    seqlist.saveFasta(seqfile=fasout,log=True)
                ## >>> ------------ <<< ##
                rje.delimitedFileOutput(self,tdtout,self.list['Headers'],rje_backup=True)
                (sx,seqx) = (0.0,len(self.dict['Mapping']))
                for est in rje.sortKeys(self.dict['Mapping']):
                    self.progLog('\r#MAP','EST->Consensus mapping: %.2f%%' % (sx/seqx))
                    rje.delimitedFileOutput(self,tdtout,self.list['Headers'],datadict=self.dict['Mapping'][est])
                self.printLog('\r#MAP','EST->Consensus mapping output to %s' % tdtout)
            else: self.dict['Mapping'] = rje.dataDict(self,tdtout,['EST'],'All')
            ### ~ [3] Annotate consensus ESTs with BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist.info['Name'] = fasout
            seqlist.loadSeqs(filetype='fas',seqtype='dna',aln=False,nodup=False)
            if self.opt['Annotate']: self.annotateEST(seqlist,stype='trans')
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <6> ### HAQESAC Annotation Methods                                                                              #
#########################################################################################################################
    def multiHAQ(self,seqlist,batchonly=False): ### Perform multiple rounds of HAQESAC analysis on final consensi
        '''
        Perform multiple rounds of HAQESAC analysis on final consensi.
        >> seqlist:SeqList object containing sequences for HAQESAC analysis
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist.cmd_list.append('dna=F')
            seqlist.loadSeqs(seqtype='Protein')
            if not seqlist.seqNum(): return self.printLog('#HAQ','No sequences for HAQESAC analysis.')
            ## ~ [0b] ~ Setup MultiHAQ object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            haqcmd = ['ac=F','blastcut=50','maketree=fasttree','bootstraps=1000','alnprog=maaft','qcover=0','qryid=0','allowvar=F','treeformats=qspec,text,nwk','group=dup']
            haqcmd += self.cmd_list + ['root=mid','qryvar=T','dna=F']
            if self.i() < 0 or rje.yesNo('Run HAQESAC in full auto mode?'): haqcmd += ['i=-1']
            mhaq = multihaq.MultiHAQ(self.log,haqcmd)
            mhaq.obj['SeqList'] = seqlist
            mhaq.info['HaqDir'] = rje.makePath('%s_HAQESAC/' % self.info['Basefile'])
            rje.mkDir(self,mhaq.info['HaqDir'])
            ## ~ [0c] ~ Setup BLAST2FAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist.list['Blast2Fas'] = [seqlist.info['Name'],self.info['SearchDB']] + self.list['HAQDB']
            ### ~ [1] ~ Run using MultiHAQ module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mhaq.blast2fas()
            mhaq.haqBatch()
            if not batchonly: mhaq.multiHAQ()
        except: self.errorLog('Problem with FIESTA.multiHAQ()')
#########################################################################################################################
### End of SECTION II: FIESTA Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
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
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: FIESTA(mainlog,['oldblast=F']+cmd_list).run()
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
