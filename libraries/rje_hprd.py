#!/usr/local/bin/python

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
Module:       rje_hprd
Description:  HPRD Database processing module
Version:      1.2.1
Last Edit:    16/05/16
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed for specific PPI Database manipulations:

    1. Parsing HPRD Flat Files. [This is the default and is run if no other option is selected.]
    Upon downloading the HPRD FLAT_FILES, this module can be run to parse out binary protein interactions and make
    sequence-specific interaction datasets. Currently, it is not clear how (if at all) the "isoform" data is used in
    HPRD for distinguishing interactions, so all isoform_1 sequences will be used for Fasta datasets. Sequences will be
    reformatted into:

    >Gene_HUMAN__AccNum Description [Gene:WWWW HPRD:XXXX; gb:YYYY; sp:ZZZZ]

    The Gene will be the HUGO gene name (as parsed from HPRD) where available, else it will be the HPRD ID. This will be
    unique for each protein and will correspond to a dataset of the same name: HPRD_Datasets/Gene_hprd.fas. All proteins
    will also be saved in a file hprd.fas. The AccNum will be UniProt if possible, else GenBank. If the option alliso=T
    is used, then all isoforms will be included and the AccNum will be X-Y where X is the HPRD ID and Y is the isoform.

    2. Converting a table of interactions into a distance matrix.
    This table should be a plain text file in which the first column is the interacting protein name and the subsequent
    columns are for the proteins (hubs) to be clustered. The first row contains their name. The rows for each spoke
    protein should be empty (or value 0) if there is no interaction and have a non-zero value if there is an interaction:
    Gene	Beta	Epsilon	Eta	Gamma	Sigma	Theta	Zeta
    AANAT							1
    
    A distance matrix is then produced (outfile=FILE => FILE.ppi_dis.txt) consisting of the number of unique interactors
    for each pairwise comparison. (The format is set by outmatrix=X : text / mysql / phylip)

    3. Incorporation of data from the GeneCards website (and Human EnsLoci) using rje_genecards. This will create a file
    called HPRD.genecards.tdt by default but this can be over-ridden using cardout=FILE. EnsLoci data will also be looked
    for in /home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas but this can be over-ridden with ensloci=FILE.

Commandline:
    ### HRPD Options ###
    hprdpath=PATH   : Path to HPRD Flat Files [./]
    genecards=T/F   : Make the HRPD.genecards.tdt file using rje_genecards (and its options) [False]
    hprdfas=T/F     : Whether to generate HPRD fasta files [False]
    alliso=T/F      : Whether to include all isoforms in the output [False]
    ppitype=LIST    : List of acceptable interaction types to parse out [in vitro;in vivo;yeast 2-hybrid]
    badtype=LIST    : List of bad interaction types, to exclude []
    domainfas=T/F   : Whether to output Domain fasta files [False]
    complexfas=T/F  : Whether to output Protein Complex fasta files [False]
    outdir=PATH     : The output directory for the files produced [./]

    ### Distance Matrix Options ###    
    ppitab=FILE     : File containing PPI data (see 2 above)
    scaled=T/F      : Whether distance matrix is to be scaled by total number of interactors in pairwise comparison [F]
    
Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_genecards, rje_ppi, rje_seq, rje_dismatrix
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_genecards, rje_seq
import rje_dismatrix_V2 as rje_dismatrix        #!# Check V2 is OK! (Might need V1)
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Working version based on rje_ppi.
    # 1.1 - Added protein complexes and PPI Types.
    # 1.2 - Added tracking of evidence.
    # 1.2.1 - Fixed "PROTEIN_ARCHITECTURE" bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : May add incorporation of POST_TRANSLATIONAL_MODIFICATIONS.txt into UniProt format output?
    # [ ] : Add data from other flat files?
    # [ ] : Read in other (custom) PPI data for focused proteins and combine with HPRD data
    # [ ] : Combine with GeneCards to make hprd_genecards file
    # [ ] : Cytoscape output
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit) = ('RJE_HPRD', '1.2.1', 'May 2016')
    description = 'HPRD Database processing module'
    author = 'Dr Richard J. Edwards.'
    return rje.Info(program,version,last_edit,description,author,time.time())
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info:
            info = makeInfo()
        if not out:
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
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: HPRD Class                                                                                               #
#########################################################################################################################
class HPRD(rje.RJE_Object):     
    '''
    HPRD Class. Author: Rich Edwards (2005).

    Info:str
    - HPRDPath = Path to HPRD Flat Files []
    - PPITab = File containing PPI data (see 2 above)
    - OutDir = The output directory for the files produced [./]

    Opt:boolean
    - ComplexFas = Whether to output Protein Complex fasta files [False]
    - DomainFas = Whether to output Domain fasta files [False]
    - GeneCards = Make the HRPD.genecards.tdt file using rje_genecards (and its options) [False]
    - HPRDFas = Whether to generate HPRD fasta files [False]
    - Scaled = Whether distance matrix is to be scaled by total number of interactors in pairwise comparison [F]
    - AllIso = Whether to include all isoforms in the output [False]

    Stat:numeric

    List:list
    - BadType = List of bad interaction types, to exclude [indirect_complex,neighbouring_reaction]
    - PPIType = List of acceptable interaction types to parse out [in vitro;in vivo;yeast 2-hybrid]

    Dict:dictionary
    - Complex = Dictionary of {Complex Name:[List of HPRD_IDs]}
    - Domains = Dictionary of {Domain Name:[List of HPRD_IDs]}
    - DomainSource = Dictionary of {Domain Name:[Sources of Domain Data]}
    - Evidence = Evidences for interactions {ID:{ID:[Evidence types]}}
    - HPRD = Data for each HPRD protein {HPRD_ID:{Seq object,GB,UniProt,Gene,Desc,Entrez,OMIM}}
    - PPI = Interaction data {HPRD_ID:[HPRD_IDs]}
    - Mapping = Dictionary of backwards mapping {Gene:HPRD_ID}  # Could try going through GB,UniProt,Gene or Entrez? #
    
    Obj:RJE_Object
    - SeqList = SeqList object containing actual sequences from PROTEIN_SEQUENCES.txt
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object'''
        ### Basics ###
        self.infolist = ['HPRDPath','PPITab','OutDir']
        self.optlist = ['ComplexFas','DomainFas','GeneCards','HPRDFas','Scaled','AllIso']
        self.statlist = []
        self.listlist = ['BadType','PPIType']
        self.dictlist = ['Complex','Domains','Evidence','HPRD','PPI','Mapping','DomainSource']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'HPRDPath':'','OutDir':''})
        self.list['PPIType'] = string.split('in vitro;in vivo;yeast 2-hybrid',';')
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
                self._cmdReadList(cmd,'path',['HPRDPath','OutDir'])
                self._cmdReadList(cmd,'file',['PPITab'])
                self._cmdReadList(cmd,'opt',['ComplexFas','DomainFas','GeneCards','HPRDFas','Scaled','AllIso'])
                self._cmdReadList(cmd,'list',['BadType','PPIType'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.list['PPIType']: self.list['PPIType'] = string.split(string.join(self.list['PPIType'],'|').lower(),'|')
        if self.list['BadType']: self.list['BadType'] = string.split(string.join(self.list['BadType'],'|').lower(),'|')
#########################################################################################################################
    ### <2> ### Main Class Run Methods                                                                                  #
#########################################################################################################################
    def run(self):  ### Main Run Method
        '''Main Run Method.'''
        try:
            ### ~ Special Run Methods First ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['PPITab'].lower() not in ['','none']: return self.ppiDisMatrix()
            ### ~ Parse HPRD Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['PPIType'] = string.split(string.join(self.list['PPIType'],';').lower(),';')
            #if self.list['PPIType'] in [[],['']]: self.list['PPIType'] = string.split('in vitro;in vivo;yeast 2-hybrid;complex',';')
            self.list['BadType'] = string.split(string.join(self.list['BadType'],';').lower(),';')
            parsecomplex = self.opt['ComplexFas'] or 'complex' in self.list['PPIType']
            if self.opt['DomainFas'] or self.opt['HPRDFas']: self.parse(parsedom=self.opt['DomainFas'],parsecomplex=parsecomplex)
            else: self.parse(parsedom=False,parseseq=False,parsecomplex='complex' in self.list['PPIType'])
            ### ~ Standard Run Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['HPRDFas']: self.saveFasta()
            if self.opt['DomainFas']: self.domainFasta()
            if self.opt['ComplexFas']: self.complexFasta()
            if self.opt['GeneCards']: self.geneCards()
        except:
            self.log.errorLog('Major problem with rje_ppi')
#########################################################################################################################
    def parse(self,parsedom=True,parseseq=True,parsecomplex=True):  ### HPRD Parsing method. Generates Mappings, HPRD data dictionary, Domain dictionary & Sequences
        '''HPRD Parsing method. Generates Mappings, HPRD data dictionary, Domain dictionary & Sequences.'''
        try:
            ### ~ Parse HPRD Mappings onto other database IDs from HPRD_ID_MAPPINGS.txt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['HPRD'] = {}
            self.dict['Mapping'] = {}
            hprd = self.loadFromFile('%sHPRD_ID_MAPPINGS.txt' % self.info['HPRDPath'],v=1,checkpath=True,chomplines=True)
            hx = float(len(hprd))
            while hprd:
                entry = hprd.pop(0)
                px = 100.0 * (hx - len(hprd)) / hx 
                self.log.printLog('\r#HPRD','Parsing HPRD_ID_MAPPINGS: %.1f%%' % px,newline=False,log=False)
                data = string.split(entry)
                ## Check ##
                if len(data) < 7: continue
                if self.dict['HPRD'].has_key(data[0]):
                    self.log.errorLog('HPRD ID %s duplicated! Aaargh!' % data[0],printerror=False)
                ## Update ##
                self.dict['HPRD'][data[0].upper()] = {'gene':data[1].upper(), 'gb':data[3], 'entrez':data[4], 'omim':data[5],
                                              'sp':data[6].upper(), 'desc':string.join(data[7:])}
                for i in [1,3,6]: self.dict['Mapping'][data[i].upper()] = data[0]
            self.log.printLog('\r#HPRD','Parsing HPRD_ID_MAPPINGS complete!')

            ### ~ Parse HPRD Domain Mappings from PROTEIN_Architecture.txt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Domains'] = {}
            self.dict['DomainSource'] = {}
            if parsedom:
                hprd = self.loadFromFile('%sPROTEIN_ARCHITECTURE.txt' % self.info['HPRDPath'],v=1,checkpath=True,chomplines=True)
                hx = float(len(hprd))
                while hprd:
                    entry = hprd.pop(0)
                    px = 100.0 * (hx - len(hprd)) / hx 
                    self.log.printLog('\r#HPRD','Parsing PROTEIN_ARCHITECTURE: %.1f%%' % px,newline=False,log=False)
                    data = string.split(entry)
                    ## Check ##
                    if len(data) < 9: continue
                    (hid,domain,type,source) = (data[0],data[4],data[5],data[8])
                    if type != 'Domain': continue
                    ## Update ##
                    if domain not in self.dict['Domains']: self.dict['Domains'][domain] = [hid]
                    elif hid not in self.dict['Domains'][domain]: self.dict['Domains'][domain].append(hid)
                    if domain not in self.dict['DomainSource']: self.dict['DomainSource'][domain] = [source]
                    elif source not in self.dict['DomainSource'][domain]: self.dict['DomainSource'][domain].append(source)
                self.log.printLog('\r#HPRD','Parsing PROTEIN_ARCHITECTURE complete!')

            ### ~ Make SeqList from PROTEIN_SEQUENCES.txt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if parseseq:
                scmd = self.cmd_list + ['autoload=T','gnspacc=F','seqin=%sPROTEIN_SEQUENCES.txt' % self.info['HPRDPath'],
                                        'autofilter=F','accnr=F','seqnr=F']
                self.obj['SeqList'] = rje_seq.SeqList(self.log,scmd)
                self.obj['SeqList'].info['Name'] = self.info['OutDir'] + 'hprd.fas'
                sx = 0.0
                for seq in self.obj['SeqList'].seq[0:]:     # seq.info['ID'] should be the HPRD ID #
                    ## Initial processing of sequence. Only keep if AllIso or isoform 1 ##
                    self.log.printLog('\r#SEQ','Processing HPRD Sequences: %.1f%%' % (sx/self.obj['SeqList'].seqNum()),newline=False,log=False)
                    iso = 'X'
                    h = seq.info['ID']
                    try: iso = rje.matchExp('^\d+\|\d+_(\d+)\|',seq.info['Name'])[0]
                    except: self.deBug(seq.info['Name'])
                    try:
                        if h not in self.dict['HPRD']:
                            self.printLog('\r#ERR','Missing from HPRD_ID_MAPPINGS?: %s' % seq.info['Name'])
                            data = string.split(seq.info['Name'],'|')
                            self.dict['HPRD'][h] = {'gene':'-', 'gb':data[2], 'entrez':'', 'omim':'','sp':'',
                                                    'desc':string.join(data[3:],'|')}      
                        if not self.opt['AllIso'] and self.dict['HPRD'][h].has_key('Seq') and iso != '1':
                            self.obj['SeqList'].seq.remove(seq)
                            continue
                        #x#if h == '00001': self.deBug('%s = %s' % (h,iso))
                        sx += 100.0
                        seq.setInfo({'Gene':self.dict['HPRD'][h]['gene'],
                                    'Description':self.dict['HPRD'][h]['desc'] + ' [Gene:%s HPRD:%s; gb:%s; sp:%s]' % (self.dict['HPRD'][h]['gene'],h,self.dict['HPRD'][h]['gb'],self.dict['HPRD'][h]['sp']),
                                    'AccNum':self.dict['HPRD'][h]['sp']})
                        ## AllIso options ##
                        if self.opt['AllIso']:
                            if 'Seq' not in self.dict['HPRD'][h]: self.dict['HPRD'][h]['Seq'] = [seq]
                            else: self.dict['HPRD'][h]['Seq'].append(seq)
                            seq.setInfo({'AccNum':'%s-%s' % (h,iso)})
                        else: self.dict['HPRD'][h]['Seq'] = seq
                        #x#print h, self.dict['HPRD'][h]['Seq']
                        ## Finish formatting ##
                        if seq.info['Gene'] == '-': self.dict['HPRD'][h]['gene'] = seq.info['Gene'] = 'HPRD' + h
                        if seq.info['AccNum'] == '-': seq.info['AccNum'] = self.dict['HPRD'][h]['gb']
                        seq.info['ID'] = '%s_HUMAN' % seq.info['Gene']
                        seq.info['Name'] = '%s__%s %s' % (seq.info['ID'],seq.info['AccNum'],seq.info['Description'])
                    except: self.errorLog('Protein Parse Error (%s)' % seq.info['Name'])
                self.log.printLog('\r#SEQ','Processing HPRD Sequences complete!')
            
            ### ~ Make PPI Data from BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []
            self.dict['PPI'] = {}
            ppi = self.loadFromFile('%sBINARY_PROTEIN_PROTEIN_INTERACTIONS.txt' % self.info['HPRDPath'],v=1,checkpath=True,chomplines=True)
            hx = float(len(ppi))
            ix = 0
            while ppi:
                entry = ppi.pop(0)
                px = 100.0 * (hx - len(ppi)) / hx 
                self.log.printLog('\r#PPI','Parsing BINARY_PROTEIN_PROTEIN_INTERACTIONS: %.1f%%' % px,newline=False,log=False)
                data = string.split(entry,'\t')
                ## Check ##
                if len(data) < 7: continue
                types = string.split(data[6],';')
                if not types: types = ['unknown']
                for type in types[0:]:
                    if type in self.list['BadType'] or (self.list['PPIType'] and type not in self.list['PPIType']): types.remove(type)
                if not types: continue
                ix += 1
                ## Update ##
                (p1,p2) = (data[1].upper(),data[4].upper())
                if p1 not in self.dict['HPRD']:
                    if p1 not in missing:
                        missing.append(p1)
                        self.log.printLog('#ERR','HPRD ID "%s" missing from HPRD_ID_MAPPINGS!' % p1,screen=False)
                    continue
                if p2 not in self.dict['HPRD']:
                    if p2 not in missing:
                        missing.append(p2)
                        self.log.printLog('#ERR','HPRD ID "%s" missing from HPRD_ID_MAPPINGS!' % p1,screen=False)
                    continue
                if not self.dict['PPI'].has_key(p1): self.dict['PPI'][p1] = []
                if p2 not in self.dict['PPI'][p1]: self.dict['PPI'][p1].append(p2)
                if not self.dict['PPI'].has_key(p2): self.dict['PPI'][p2] = []
                if p1 not in self.dict['PPI'][p2]: self.dict['PPI'][p2].append(p1)
                if p1 not in self.dict['Evidence']: self.dict['Evidence'][p1] = {}
                if p2 not in self.dict['Evidence'][p1]: self.dict['Evidence'][p1][p2] = []
                for type in types:
                    if type not in self.dict['Evidence'][p1][p2]: self.dict['Evidence'][p1][p2].append(type)
                #x#if p1 == '12422': self.deBug(self.dict['PPI'][p1])
            self.log.printLog('\r#PPI','Parsing BINARY_PROTEIN_PROTEIN_INTERACTIONS complete!')

            ### ~ Parse protein Complex data from PROTEIN_COMPLEXES.txt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Complex'] = {}
            ppi = self.loadFromFile('%sPROTEIN_COMPLEXES.txt' % self.info['HPRDPath'],v=1,checkpath=True,chomplines=True)
            hx = float(len(ppi))
            while ppi:
                entry = ppi.pop(0)
                px = 100.0 * (hx - len(ppi)) / hx 
                self.log.printLog('\r#PPI','Parsing PROTEIN_COMPLEXES: %.1f%%' % px,newline=False,log=False)
                data = string.split(entry)
                ## Check ##
                if len(data) < 5: continue
                ## Update ##
                (complex,hprd) = (data[0],data[1])
                if hprd == 'None': continue
                if not self.dict['Complex'].has_key(complex): self.dict['Complex'][complex] = []
                if hprd not in self.dict['Complex'][complex]: self.dict['Complex'][complex].append(hprd)
                #x#if p1 == '12422': self.deBug(self.dict['PPI'][p1])
            self.log.printLog('\r#PPI','Parsing PROTEIN_COMPLEXES complete!')

            ### ~ Update PPI from protein Complex data if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            type = 'complex'
            if type not in self.list['BadType'] and (not self.list['PPIType'] or type in self.list['PPIType']):
                cx = 0.0
                for complex in self.dict['Complex']:
                    self.log.printLog('\r#PPI','Adding protein complex data to PPI: %.1f%%' % (cx/len(self.dict['Complex'])),newline=False,log=False)
                    cx += 100.0
                    for p1 in self.dict['Complex'][complex]:
                        for p2 in self.dict['Complex'][complex]:
                            if not self.dict['PPI'].has_key(p1): self.dict['PPI'][p1] = []
                            if p2 not in self.dict['PPI'][p1]: self.dict['PPI'][p1].append(p2)
                            if p1 not in self.dict['Evidence']: self.dict['Evidence'][p1] = {}
                            if p2 not in self.dict['Evidence'][p1]: self.dict['Evidence'][p1][p2] = []
                            if type not in self.dict['Evidence'][p1][p2]: self.dict['Evidence'][p1][p2].append(type)
                self.log.printLog('\r#PPI','Added protein complex data to PPI for %s complexes' % rje.integerString(len(self.dict['Complex'])))
            ptxt = '%s proteins; %s interactions' % (rje.integerString(len(self.dict['PPI'])),rje.integerString(ix))
            self.log.printLog('\r#PPI','Parsing interactions complete: %s.' % ptxt)
            if missing: open('HPRD.missing.txt','w').write(string.join(missing,'\n'))
        except:
            self.log.errorLog('Error in HPRD.parse()',printerror=True,quitchoice=False)
            raise
#########################################################################################################################
    def saveFasta(self):    ### Outputs parsed PPI datasets in Fasta format
        '''Outputs parsed PPI datasets in Fasta format.'''
        try:
            ### Setup ###
            datpath = self.info['OutDir'] + rje.makePath('HPRD_Datasets/')
            rje.mkDir(self,datpath)
            ## Check Seqs ##
            for p1 in rje.sortKeys(self.dict['PPI']):
                if 'Seq' not in self.dict['HPRD'][p1]:      #!# KeyError #!#
                    print p1, self.dict['HPRD'][p1]
                    self.deBug('No Seq for %s' % p1)

            ### All sequences ###
            self.obj['SeqList'].saveFasta()
            ### Output PPI Datasets ###
            for p1 in rje.sortKeys(self.dict['PPI']):
                mylist = []
                for p2 in self.dict['PPI'][p1]:
                    if self.opt['AllIso']: mylist += self.dict['HPRD'][p2]['Seq']
                    else: mylist.append(self.dict['HPRD'][p2]['Seq'])
                sfile = '%s%s_hprd.fas' % (datpath,self.dict['HPRD'][p1]['gene'])
                if mylist: self.obj['SeqList'].saveFasta(seqs=mylist,seqfile=sfile)
            self.log.printLog('#FAS','HPRD PPI fasta output complete.')
        except: self.log.errorLog('Error in HPRD.saveFasta()',printerror=True,quitchoice=False)
#########################################################################################################################
    def domainFasta(self):    ### Outputs parsed domain and domain PPI datasets in Fasta format
        '''Outputs parsed PPI datasets in Fasta format.'''
        try:
            ### ~ Tab delimited domain-HPRD pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['Domain','HPRD','Gene']
            dfile = self.info['OutDir'] + 'HPRD.domains.tdt'
            rje.delimitedFileOutput(self,dfile,headers,'\t')
            sfile = self.info['OutDir'] + 'HPRD.domsource.tdt'
            shead = ['Domain','Source']
            rje.delimitedFileOutput(self,sfile,shead,'\t')
            dx = 0.0
            for domain in rje.sortKeys(self.dict['Domains']):
                self.log.printLog('\r#DOM','HPRD Domain output (%s): %.1f%%' % (dfile,dx/len(self.dict['Domains'])),newline=False,log=False)
                dx += 100.0
                for hid in self.dict['Domains'][domain]:
                    datadict = {'Domain':domain,'HPRD':hid,'Gene':self.dict['HPRD'][hid]['gene']}
                    rje.delimitedFileOutput(self,dfile,headers,'\t',datadict)
                for source in self.dict['DomainSource'][domain]:
                    datadict = {'Domain':domain,'Source':source}
                    rje.delimitedFileOutput(self,sfile,shead,'\t',datadict)
            self.log.printLog('\r#DOM','HPRD Domain output (%s): %s domains.' % (dfile,rje.integerString(len(self.dict['Domains']))))
                       
            ### ~ Domain PPI Dataset Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datpath = self.info['OutDir'] + rje.makePath('HPRD_Domain_Datasets/')
            rje.mkDir(self,datpath)
            for domain in rje.sortKeys(self.dict['Domains']):
                ## Generate a list of all interactors with domain-containing proteins ##
                plist = []
                for p1 in self.dict['Domains'][domain]:
                    if p1 not in self.dict['PPI']: continue
                    for p2 in self.dict['PPI'][p1]:
                        if p2 not in plist: plist.append(p2)
                plist.sort()
                ## Generate Sequence list and output ##
                mylist = []
                for p in plist:
                    if self.opt['AllIso']: mylist += self.dict['HPRD'][p]['Seq']
                    else: mylist.append(self.dict['HPRD'][p]['Seq'])
                sfile = '%s%s_hprd.fas' % (datpath,domain)
                if mylist: self.obj['SeqList'].saveFasta(seqs=mylist,seqfile=sfile)
                else: self.log.printLog('#DOM','No PPI partners for domain "%s"' % domain)
            self.log.printLog('\r#DOM','HPRD Domain fasta output complete.')
        except:
            self.log.errorLog('Error in HPRD.saveFasta()',printerror=True,quitchoice=False)
            raise
#########################################################################################################################
    def complexFasta(self):     ### Outputs parsed complex datasets in Fasta format
        '''Outputs parsed complex datasets in Fasta format.'''
        try:
            ### Setup ###
            datpath = self.info['OutDir'] + rje.makePath('HPRD_Complexes/')
            rje.mkDir(self,datpath)

            ### Output PPI Datasets ###
            for complex in rje.sortKeys(self.dict['Complex']):
                mylist = []
                for p2 in self.dict['Complex'][complex]:
                    if self.opt['AllIso']: mylist += self.dict['HPRD'][p2]['Seq']
                    else: mylist.append(self.dict['HPRD'][p2]['Seq'])
                sfile = '%s%s_hprd.fas' % (datpath,complex)
                if mylist: self.obj['SeqList'].saveFasta(seqs=mylist,seqfile=sfile)
            self.log.printLog('#FAS','HPRD complex fasta output complete.')
        except:
            self.log.errorLog('Error in HPRD.complexFasta()',printerror=True,quitchoice=False)
            raise
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def geneCards(self):    ### Combines parsed HPRD data with Genecards into HPRD.genecards.tdt
        '''Combines parsed HPRD data with Genecards into HPRD.genecards.tdt.'''
        try:
            ### Setup Basic GeneCards ###
            cfile = 'HPRD.genecards.tdt'
            gcmd = ['species=Human','cardout=%s' % cfile,'ensloci=/home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas']
            cards = rje_genecards.GeneCards(self.log,gcmd+self.cmd_list)
            cards.setup()

            ### Reconfigure parsed HPRD for genecards ###
            self.addToGeneCards(cards)

            ## Add GeneCards Data and Output ##
            cards.run(setup=False)

        except:
            self.log.errorLog('Error in HPRD.parse()',printerror=True,quitchoice=False)
            raise
#########################################################################################################################
    def addToGeneCards(self,cards,addcards=True): ### Reconfigures and adds parsed HPRD data to GeneCards
        '''
        Reconfigures and adds parsed HPRD data to GeneCards.
        >> cards:rje_genecards.GeneCards object
        >> addcards:boolean [True] = whether to add genes from HPRD to the GeneCards dictionary
        '''
        ### Add relevant headers for future output ###
        for h in ['HPRD','OMIM','EntrezCheck','Desc']:
            if h not in cards.list['Headers']:
                cards.list['Headers'].append(h)
            for gene in cards.list['Genes']:
                if h not in cards.dict['GeneCard'][gene]: cards.dict['GeneCard'][gene][h] = ''
        ### Add to GeneCards ###
        (hx,htot) = (0.0,len(self.dict['HPRD']))
        for hprd in self.dict['HPRD']:
            self.log.printLog('\r#HPRD','Adding HPRD to GeneCards: %.1f%%' % (hx/htot),newline=False,log=False)
            hx += 100.0
            self.deBug(self.dict['HPRD'][hprd])
            gene = self.dict['HPRD'][hprd]['gene']
            omim = self.dict['HPRD'][hprd]['omim']
            entrez = self.dict['HPRD'][hprd]['entrez']
            if gene in cards.list['Genes']:
                if cards.dict['GeneCard'][gene]['HPRD'] == '': cards.dict['GeneCard'][gene]['HPRD'] = hprd
                elif hprd not in string.split(cards.dict['GeneCard'][gene]['HPRD'],','):
                    cards.dict['GeneCard'][gene]['HPRD'] = string.join(string.split(cards.dict['GeneCard'][gene]['HPRD'],',')+[hprd],',')
                if cards.dict['GeneCard'][gene]['OMIM'] == '': cards.dict['GeneCard'][gene]['OMIM'] = omim
                elif omim not in string.split(cards.dict['GeneCard'][gene]['OMIM'],','):
                    cards.dict['GeneCard'][gene]['OMIM'] = string.join(string.split(cards.dict['GeneCard'][gene]['OMIM'],',')+[omim],',')
                if cards.dict['GeneCard'][gene]['EntrezCheck'] == '': cards.dict['GeneCard'][gene]['EntrezCheck'] = entrez
                elif entrez not in string.split(cards.dict['GeneCard'][gene]['EntrezCheck'],','):
                    cards.dict['GeneCard'][gene]['EntrezCheck'] = string.join(string.split(cards.dict['GeneCard'][gene]['EntrezCheck'],',')+[entrez],',')
            elif addcards:
                if gene == '-': gene = 'HPRD' + hprd
                cards.list['Genes'].append(gene)
                cards.dict['GeneCard'][gene] = {'Symbol':'!FAILED!','HPRD':hprd,'OMIM':omim,'EntrezCheck':entrez,'Desc':self.dict['HPRD'][hprd]['desc']} 
        self.log.printLog('\r#HPRD','Added %s HPRD genes to GeneCards.' % (rje.integerString(htot)))
#########################################################################################################################
    def ppiDisMatrix(self): ### Converts PPI Table into distance matrix
        '''Converts PPI Table into distance matrix.'''
        try:
            ### Check File ###
            if not os.path.exists(self.info['PPITab']):
                self.log.errorLog('PPI Table file "%s" missing!' % self.info['PPITab'],printerror=False)
                return False

            ### Setup ###            
            data = rje.dataDict(self,self.info['PPITab'],getheaders=True)
            headers = data.pop('Headers')
            ppidis = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            ppidis.opt['Symmetric'] = True
            ppidis.setInfo({'Name':'%s.ppi_dis.txt' % self.info['Basefile'],'Type':'PPI'})

            ### Make DisMatrix ###
            for p1 in headers[1:]:
                ppidis.addDis(p1,p1,0)
                for p2 in headers[headers.index(p1)+1:]:
                    ppi = 0
                    unique = 0
                    for i in data.keys():
                        try:
                            v1 = int(data[i][p1])
                        except:
                            v1 = data[i][p1]
                        try:
                            v2 = int(data[i][p2])
                        except:
                            v2 = data[i][p2]
                        if v1 or v2:
                            ppi += 1
                            if not (v1 and v2):
                                unique += 1
                    if self.opt['Scaled']:
                        ppidis.addDis(p1,p2,float(unique)/float(ppi))
                    else:
                        ppidis.addDis(p1,p2,unique)

            ### Output ###
            delimit = rje.getDelimit(self.cmd_list,default=',')
            ppidis.saveMatrix(headers[1:],ppidis.info['Name'],delimit)
            
        except:
            self.log.errorLog('Major problem with rje_ppi.ppiDisMatrix')
            return False
#########################################################################################################################
### End of SECTION II: HPRD Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: PPI Dictionary  METHODS                                                                                #
#########################################################################################################################
# See rje_ppi
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try: HPRD(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
