#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_genbank
Description:  RJE GenBank Module
Version:      1.5.2
Last Edit:    05/09/16
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for parsing information out of GenBank files and converting them to other formats.

Input Options:
    seqin=FILE      : Input Genbank file []
    fetchuid=LIST   : Genbank retrieval to of a list of nucleotide entries to generate seqin=FILE []
    spcode=X        : Overwrite species read from file (if any!) with X [None]
    taxdir=PATH     : Path to taxonomy files for species code extraction. (Will not use if blank or None) [./SourceData/]
    addtags=T/F     : Add locus_tag identifiers if missing - needed for gene/cds/prot fasta output [False]

Output Options:
    basefile=FILE   : Root of output file names (same as input file by default) []
    tabout=T/F      : Delimited table output of features [False]
    features=LIST   : Subset of features to extract from Genbank file (blank for all) []
    details=LIST    : List of feature details to extract into own columns []
    detailskip=LIST : Subset of feature details to exclude from extraction [translation]
    fasout=LIST     : Types of sequences to output into files (full/gene/cds/prot) as *.*.fas []
    geneacc=X       : Feature detail to use for gene sequence accession number (added to details) [locus_tag]
    protacc=X       : Feature detail to use for protein sequence accession number (added to details) [protein_id]
    locusout=T/F    : Whether to generate output by locus (True, locus as basefile) or combined (False) [False]
    locusdir=PATH   : Directory in which to generate output by locus [./]
    
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_sequence, rje_taxonomy, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Modified and Tidied output a little.
    # 0.2 - Added details to skip and option to use different detail for protein accession number.
    # 0.3 - Added reloading of features.
    # 1.0 - Basic functioning version. Added fetchuid=LIST Genbank retrieval to generate seqin=FILE.
    # 1.1 - Added use of rje_taxonomy for getting Species Code from TaxID.
    # 1.2 - Modified to deal with genbank protein entries.
    # 1.2.1 - Fixed feature bug that was breaking parser and removing trailing '*' from protein sequences.
    # 1.2.2 - Fixed more features that were breaking parser.
    # 1.3.0 - Added split viral output.
    # 1.3.1 - Fixed bug in split viral output.
    # 1.3.2 - Fixed bug in reverse complement sequences with introns.
    # 1.4.0 - Added addtags=T/F : Add locus_tag identifiers if missing - needed for gene/cds/prot fasta output [False]
    # 1.4.1 - Fixed genetic code warning.
    # 1.5.0 - Added setupRefGenome() method based on PAGSAT code.
    # 1.5.1 - Fixed logskip append locus sequence file bug.
    # 1.5.2 - Fixed addtag(s) bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add conversion of TaxID to SpCode using rje_taxonomy.
    # [ ] : Add dealing with alternative genetic code.
    # [ ] : Fix appending to protein sequence (by locus) error.
    # [ ] : Add more codon tables - some just have different initiator codons.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_GenBank', '1.5.2', 'September 2016', '2011')
    description = 'RJE GenBank Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
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
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
nonsplit_features = ['proviral','pseudo','gene','focus','ribosomal_slippage','trans_splicing']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GenBank Class                                                                                           #
#########################################################################################################################
class GenBank(rje_obj.RJE_Object):     
    '''
    GenBank Class. Author: Rich Edwards (2011).

    Str:str
    - FetchUID = If FetchUID is a file, this file will be stored here and can be used to name SeqIn []
    - GeneAcc = Feature detail to use for gene sequence accession number [locus_tag]
    - ProtAcc = Feature detail to use for protein sequence accession number [locus_tag]
    - SeqIn = Input Genbank file []
    - SpCode = Overwrite species read from file (if any!) with X [None]
    - TaxDir = Path to taxonomy files for species code extraction. (Will not use if blank or None) [./SourceData/]
    - LocusDir = Directory in which to generate output by locus [./]

    Bool:boolean
    - AddTags=T/F      : Add locus_tag identifiers if missing - needed for gene/cds/prot fasta output [False]
    - TabOut = Delimited table output of features [False]
    - LocusOut = Whether to generate output by locus (True, accnum basefile) or combined (False) [False]

    Int:integer

    Num:float
    
    List:list
    - Details = List of feature details to extract into own columns []
    - DetailSkip = Subset of feature details to exclude from extraction [translation]
    - FasOut = Types of sequences to output into files (full/gene/cds/prot) as *.*.fas
    - Features = Subset of features to extract from Genbank file []
    - FetchUID = Genbank retrieval to of a list of nucleotide entries to generate seqin=FILE []

    Dict:dictionary
    - Sequence = Dictionary of {Locus:Sequence}

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - Taxonomy = rje_taxonomy.Taxonomy object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['GeneAcc','ProtAcc','SeqIn','SpCode','FetchUID','TaxDir','LocusDir']
        self.boollist = ['AddTags','TabOut','LocusOut']
        self.intlist = []
        self.numlist = []
        self.listlist = ['Details','DetailSkip','FasOut','Features','FetchUID']
        self.dictlist = ['Sequence']
        self.objlist = ['DB','Taxonomy']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'GeneAcc':'locus_tag','ProtAcc':'protein_id','TaxDir':'./SourceData/','LocusDir':''})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        #self.list['FasOut'] = string.split('full,gene,cds,prot',',')
        self.list['DetailSkip'] = ['translation']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
        self._setForkAttributes()   # Delete if no forking
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
                ### Class Options ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'file',['SeqIn','FetchUID'])
                self._cmdReadList(cmd,'path',['LocusDir'])
                self._cmdReadList(cmd,'str',['GeneAcc','ProtAcc','SpCode','TaxDir'])
                self._cmdReadList(cmd,'bool',['AddTags','TabOut','LocusOut'])
                self._cmdReadList(cmd,'list',['Details','DetailSkip','FasOut','Features','FetchUID'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStr('SpCode').lower() in ['','none']: self.str['SpCode'] = ''
        else:
            self.printLog('#SPCODE',self.getStr('SpCode'))
            if self.getStrLC('TaxDir') and not rje.exists(self.getStr('TaxDir')) and self.i() >= 0:
                if rje.yesNo('Species code given and TaxDir not found. Set taxdir=None? (Otherwise, will create.)'):
                    self.setStr({'TaxDir':'None'}); self.cmd_list.append('taxdir=None')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.fetchGenbank(): return self.printLog('#FAIL','Failed to fetch uid from Genbank nuccore.')
            if not self.parseGenbank(): return self.printLog('#FAIL','Failed to parse Genbank file.')
            if not self.addTags(): return self.printLog('#FAIL','Failed to deal with locus_tags for fasta output.')
            if self.getBool('TabOut'):
                self.processFeatures()
                if self.getBool('LocusOut'):
                    tables = self.db().tables()[0:]
                    ldb = self.db('Locus')
                    tables.remove(ldb)
                    ldb.saveToFile()
                    self.db().saveDB(tables,splitfield='locus',outdir=self.getStr('LocusDir'),logskip=self.v()>0,replace=self.force(),log=self.v()>0)
                    self.printLog('#OUT','%s tables output to %s for %s loci' % (string.join(self.db().tableNames(),', '),self.getStr('LocusDir'),ldb.entryNum()))
                else: self.db().saveDB(logskip=self.v()>0,replace=self.force())
            if self.list['FasOut']: self.saveFasta(logskip=self.v()>0)
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Basefile').lower() in ['','none']:
                self.setBasefile(rje.baseFile(self.getStr('SeqIn'),strip_path=False,extlist=['gbk','gb']))
                self.cmd_list.append('basefile=%s' % self.getStr('Basefile'))
            self.list['Features'] = string.split(string.join(self.list['Features']).lower())
            locdb = self.db().addEmptyTable('Locus',['locus','length','type','definition','accession','version','gi','organism','spcode'],['locus'])
            #self.list['Details'] = string.split(string.join(self.list['Details']).lower())
            if self.getStr('GeneAcc') and self.getStr('GeneAcc') not in self.list['Details']: self.list['Details'].append(self.getStr('GeneAcc'))
            if self.getStr('ProtAcc') and self.getStr('ProtAcc') not in self.list['Details']: self.list['Details'].append(self.getStr('ProtAcc'))
            ftfields = ['locus','feature','position','start','end'] + self.list['Details'] + ['details']
            ftdb = self.db().addEmptyTable('Feature',ftfields,['locus','feature','position'])
            if not self.getStrLC('SpcCode'):
                tax = self.obj['Taxonomy']
                if not tax: tax = self.obj['Taxonomy'] = rje_taxonomy.Taxonomy(self.log,self.cmd_list)
                #if self.getStrLC('TaxDir'): tax.setup()
                #else:
                #tax.setBool({'Setup':False})
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def fetchGenbank(self,uid=[],force=None):  ### Fetches GenBank entries into SeqIn file.
        '''
        Fetches GenBank entries into SeqIn file.
        >> uid:list [] = List of Nucleotide IDs. Will use self.list['FetchUID'] if empty.
        >> force:bool [None] = Whether to force download and replace SeqIn. Will use self.force() if None.
        << True if ran OK (including returning nothing for empty list), False if issues.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if force == None: force = self.force()
            if not uid: uid = self.list['FetchUID']
            while '' in uid: uid.remove('')
            if not uid: return True     # Nothing to fetch!
            uidunique = rje.sortUnique(uid)
            if len(uid) != len(uidunique):
                self.warnLog('%s duplicate UID identified. Download order may change.' % (len(uid) - len(uidunique)))
                uid = uidunique
            uidstr = string.join(uid,',')
            if len(string.split(uidstr)) > 1:
                self.warnLog('Spaces found in 1+ FetchUID accession numbers and will be removed.')
                uidstr = string.join(string.split(uidstr),'')
            if not self.getStrLC('SeqIn'):
                if not self.getStrLC('Basefile'):
                    if rje.exists(self.getStr('FetchUID')): self.setBasefile(rje.baseFile(self.getStr('FetchUID'),strip_path=True))
                    elif len(uid): self.setBasefile(uid[0])
                    else: self.setBasefile('fetchuid')
                self.setStr({'SeqIn':self.baseFile()+'.gb'})
            if rje.exists(self.getStr('SeqIn')) and not force:
                self.printLog('#SEQIN','%s already found (force=F). Skipping Genbank Fetch' % self.getStr('SeqIn'))
                return True
            rje.backup(self,self.getStr('SeqIn'))
            ## ~ [1a] EFetch URL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#FETCH','%s Genbank IDs to download into %s.' % (rje.iLen(uid),self.getStr('SeqIn')))
            if len(uid) > 200: return self.fetchSplitGenbank(uid)
            baseurl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=text&id='
            fullurl = baseurl + uidstr
            if len(string.split(fullurl)) > 1:
                self.errorLog('Spaces present in UIDs. Cannot fetch.',printerror=False); return False
            ### ~ [2a] Fetch UID into file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('OSX'): os.system('curl "%s" -o "%s"' % (fullurl,self.getStr('SeqIn')))
            else: os.system('wget "%s" -O "%s"' % (fullurl,self.getStr('SeqIn')))
            self.printLog('#FETCH','%s Genbank IDs downloaded into %s.' % (rje.iLen(uid),self.getStr('SeqIn')))
            if rje.exists(self.getStr('SeqIn')) and open(self.getStr('SeqIn'),'r').read().startswith('ID list is empty!'):
                os.unlink(self.getStr('SeqIn'))
                self.warnLog('Nothing downloaded for %s' % uidstr)
                raise ValueError('GenBank reported that ID list is empty!')
            return rje.exists(self.getStr('SeqIn'))
            #curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002031&rettype=gb&retmode=text" -o test.gb
        except: self.errorLog('Problem during %s.fetchGenbank().' % self); return False  # Setup failed
#########################################################################################################################
    def fetchSplitGenbank(self,uid):  ### Fetches GenBank entries into SeqIn file.
        '''
        Fetches GenBank entries into SeqIn file.
        >> uid:list [] = List of Nucleotide IDs.
        << True if ran OK (including returning nothing for empty list), False if issues.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.getStr('SeqIn')
            splitx = 0; uidx = len(uid)
            ### ~ [2] Fetch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while uid:
                if splitx: self.setStr({'SeqIn':'%s.tmp.%d.gb' % (self.basefile(),splitx)})
                self.fetchGenbank(uid[:200],force=True)
                if splitx: rje.fileTransfer(fromfile=self.getStr('SeqIn'),tofile=seqin,deletefrom=True,append=True)
                uid = uid[200:]
                splitx += 1
            self.setStr({'SeqIn':seqin})
            self.printLog('#FETCH','%s Genbank IDs downloaded in %s batches into %s.' % (rje.iStr(uidx),splitx,self.getStr('SeqIn')))
            return rje.exists(self.getStr('SeqIn'))
        except: self.errorLog('Problem during %s.fetchSplitGenbank().' % self); return False  # Setup failed
#########################################################################################################################
    def addTags(self):  ### Adds locus_tags to features if required
        '''Adds locus_tags to features if required.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdb = self.db('Feature')
            ## ~ [1a] Establish presence of locus_tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ftags = {}  # Dictionary of feature type and whether all have locus_tags
            for ftentry in ftdb.entries():
                ftype = ftentry['feature']
                if ftype not in ftags: ftags[ftype] = True
                ftags[ftype] = ftags[ftype] and ftentry['locus_tag']
            ## ~ [1b] Establish need to add tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            addtags = self.getBool('AddTags')
            for (ftype,fasout) in [('CDS','prot'),('CDS','cds'),('gene','gene')]:
                if addtags: break
                if fasout in self.list['FasOut'] and ftype in ftags and not ftags[ftype]:
                    if not self.getBool('AddTags'):
                        self.setBool({'AddTags':rje.yesNo('%s locus_tags required for %s fasta output. Set addtags=T?' % (ftype,fasout))})
                        addtags = self.getBool('AddTags')
                        if not addtags and rje.yesNo('Switch off %s fasta output. (May get strange behaviour otherwise.)' % fasout):
                            self.list['FasOut'].remove(fasout)
                        if not addtags: self.warnLog('%s fasta output but addtags=False.' % fasout)
            if not addtags: return True
            ### ~ [2] Add Tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.processFeatures()
            ftdb.newKey(['locus','feature','start','end'])
            for ftype in ftdb.index('feature'):
                loctot = {}; locx = {}; addx = 0
                floci = ftdb.indexDataList('feature',ftype,'locus',sortunique=False)
                for locus in floci: loctot[locus] = floci.count(locus); locx[locus] = 0
                for ftentry in ftdb.indexEntries('feature',ftype):
                    locus = ftentry['locus']
                    locx[locus] += 1
                    if not ftentry['locus_tag']: ftentry['locus_tag'] = '%s%s%s' % (locus,ftype.upper(),rje.preZero(locx[locus],loctot[locus])); addx += 1
                self.printLog('#TAG','Added %s %s locus_tag.' % (rje.iStr(addx),ftype))
            return True
        except: self.errorLog('Problem during %s.addTags().' % self); return False  # Setup failed
#########################################################################################################################
    def parseGenbank(self,filename=None): ### Parses details from GenBank file into attributes
        '''Parses details from GenBank file into attributes.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locdb = self.db('Locus'); ftdb = self.db('Feature')
            locus = ''; features = False; entry = {}; ex = 0; fx = 0
            if not filename: filename = self.getStr('SeqIn')
            if not rje.checkForFile(filename):
                try: raise IOError
                except: self.errorLog('Genbank file "%s" not found!' % filename); return False
            loadtxt = 'GenBank file %s' % rje.baseFile(filename,strip_path=True,extlist=['gbk'])
            ftdic = {}  # Counts of feature types
            if self.getStrLC('TaxDir'):
                self.obj['Taxonomy'].log.no_suppression.append("Invented SpCode")
                self.obj['Taxonomy'].setup()
            ### ~ [2] Load file contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            FILE = open(filename,'r'); line = FILE.readline(); nx = 0; features = False; readft = False; fentry = {}
            while line:
                self.progLog('#GB','Loading %s: %s entries; %s features; %s nt' % (loadtxt,rje.iStr(ex),rje.iStr(fx),rje.iStr(nx)))
                data = string.split(rje.chomp(line))
                if not data: line = FILE.readline(); continue
                cntd = line[:12] == '            '
                if not cntd: ltype = data[0]
                if line.startswith(ltype): readft = False     # Stop reading FEATURES
                if ltype == 'LOCUS':
                    locus = data[1]
                    if locus in self.dict['Sequence']:
                        self.debug(locus)
                        self.debug(self.dict['Sequence'])
                        raise ValueError
                    entry = {'locus':data[1],'length':string.atoi(data[2]),'type':data[4]}
                    if features or fentry: raise ValueError('Feature processing error')
                elif ltype in ['DEFINITION','ACCESSION']:
                    if cntd: entry[ltype.lower()] = string.join([entry[ltype.lower()]] + data)
                    else: entry[ltype.lower()] = string.join(data[1:])
                elif ltype == 'VERSION' and not cntd:
                    entry['version'] = data[1]
                    try: entry['gi'] = string.split(data[2],'GI:')[1]
                    except: entry['gi'] = '-'
                elif ltype in ['ORGANISM'] and not cntd: entry[ltype.lower()] = string.join(data[1:])
                elif ltype == 'FEATURES': features = True; readft = True; ftentry = {}
                elif ltype == 'ORIGIN': self.dict['Sequence'][locus] = ''
                elif ltype == '//':
                    if self.getStr('SpCode'): entry['spcode'] = self.getStr('SpCode')
                    elif self.obj['Taxonomy'].getBool('Setup'): entry['spcode'] = self.obj['Taxonomy'].getSpCode(entry['organism'])
                    else: entry['spcode'] = rje_sequence.getSpecCode(entry['organism'])
                    if features and ftentry:
                        while string.count(ftentry['position'],'(') != string.count(ftentry['position'],')'):
                            if ftentry['locus'] in ['NC_006146','NC_007044'] and ftentry['details'].find('CeHV15gLMP1')>0:
                                self.debug('Pos continued: "%s"' % ftentry['details'])
                            details = string.split(ftentry['details'],' /')
                            if '..' not in details[0]: break
                            ftentry['position'] += string.replace(details.pop(0),' ','')
                            ftentry['details'] = string.join(['']+details,' /')
                            #if ftentry['locus'] in ['NC_006146','NC_007044']: self.debug(ftentry)
                        if string.count(ftentry['position'],'(') != string.count(ftentry['position'],')'):
                            self.warnLog('Unmatched parentheses in Feature position: %s' % ftentry['position'])
                        if ftentry['locus'] in ['NC_006146','NC_007044'] and ftentry['details'].find('CeHV15gLMP1')>0: self.debug(ftentry)
                        ftdb.addEntry(ftentry); fx += 1
                        if ftentry['feature'] not in ftdic: ftdic[ftentry['feature']] = 0
                        ftdic[ftentry['feature']] += 1
                        ftentry = {}
                    locdb.addEntry(entry); ex += 1; locus = ''; features = False; entry = {}
                elif locus in self.dict['Sequence']:
                    seqline = string.join(data[1:],'')
                    self.dict['Sequence'][locus] += seqline
                    nx += len(seqline)
                elif readft:
                    if cntd and ftentry:
                        ftentry['details'] = string.join([ftentry['details']] + data)
                        #self.debug('"%s"' % ftentry['details'])
                    elif ftentry:
                        while string.count(ftentry['position'],'(') != string.count(ftentry['position'],')'):
                            if ftentry['locus'] in ['NC_006146','NC_007044'] and ftentry['details'].find('CeHV15gLMP1')>0:
                                self.debug('Pos continued: "%s"' % ftentry['details'])
                            details = string.split(ftentry['details'],' /')
                            if '..' not in details[0]: break
                            ftentry['position'] += string.replace(details.pop(0),' ','')
                            ftentry['details'] = string.join(['']+details,' /')
                            #if ftentry['locus'] in ['NC_006146','NC_007044']: self.debug(ftentry)
                        if string.count(ftentry['position'],'(') != string.count(ftentry['position'],')'):
                            self.warnLog('Unmatched parentheses in Feature position: %s' % ftentry['position'])
                        if ftentry['locus'] in ['NC_006146','NC_007044'] and ftentry['details'].find('CeHV15gLMP1')>0: self.debug(ftentry)
                        ftdb.addEntry(ftentry); fx += 1
                        if ftentry['feature'] not in ftdic: ftdic[ftentry['feature']] = 0
                        ftdic[ftentry['feature']] += 1
                        ftentry = {}
                    if not cntd:
                        if self.list['Features'] and ltype.lower() not in self.list['Features']: ftentry = {}
                        else: ftentry = {'locus':locus,'feature':ltype,'position':data[1],'details':''}
                        #self.debug('"%s"' % ftentry['details'])
                line = FILE.readline()
            FILE.close()
            self.printLog('\r#GB','Loaded %s: %s entries; %s features; %s nt in total.' % (loadtxt,rje.iStr(ex),rje.iStr(fx),rje.iStr(nx)))
            for ft in rje.sortKeys(ftdic): self.printLog('#FT','%s: %s parsed' % (ft,rje.iStr(ftdic[ft])))
            return True
        except: self.errorLog('%s.parseGenBank error' % self); return False
#########################################################################################################################
    def processFeatures(self):  ### Extracts relevant information from features for *.Feature.tdt
        '''Extracts relevant information from features for *.Feature.tdt.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('Feature')
            for ftkey in fdb.dataKeys():
                ftentry = fdb.data(ftkey)
                ftdic = self.featureDict(ftentry)
                for field in self.list['Details'] + ['start','end','details']: ftentry[field] = ftdic[field]
                for field in ['start','end']: ftentry[field] = int(ftentry[field])
        except: self.errorLog('%s.processFeatures error' % self)
#########################################################################################################################
    def featureDict(self,ftentry):  ### Converts feature entry to dictionary of data (including sequence)
        '''Converts feature entry to dictionary of data (including sequence).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ftentry['locus'] in ['NC_006146','NC_007044'] and ftentry['details'].find('CeHV15gLMP1')>0: self.debug(ftentry)
            locdb = self.db('Locus')
            locus = ftentry['locus']
            entry = locdb.data(locus)
            ### ~ [1] ~ Make dictionary of data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdic = rje.combineDict(rje.combineDict({},ftentry),entry)
            if ftdic['details'] and not ftdic['details'].startswith(' '): ftdic['details'] = ' %s' % ftdic['details']
            if ftdic['details'] and not ftdic['details'].startswith(' /'):
                self.warnLog('%s %s feature missing leading " /": "%s"' % (locus,ftdic['feature'],ftdic['details']),"bad_feature",suppress=True)
                self.debug(ftdic)
            details = ftdic.pop('details')  # Should now start with " /"
            ftdic['details'] = []
            ftdata = []
            #self.debug('"%s"' % details)
            datasplit = string.split(details,' /')
            #self.debug(datasplit)
            if datasplit[0]: self.warnLog('%s %s feature details "%s" not properly recognised for parsing.' % (locus,ftdic['feature'],details))
            else: datasplit = datasplit[1:]
            for data in datasplit: # This starts with 1: assuming that the first feature starts with /
                if data in nonsplit_features or '=' in data: ftdata.append(data)
                elif not ftdata:
                    self.warnLog('Odd %s %s feature details "%s" not properly recognised for parsing.' % (locus,ftdic['feature'],details))
                    if data and rje.yesNo('Add %s to non-split (no "=") feature list?' % data): nonsplit_features.append(data); ftdata.append(data)
                else: ftdata[-1] = '%s/%s' % (ftdata[-1],data)
            for data in ftdata:
                if data in nonsplit_features: dtype = content = data
                else:
                    try:
                        sdata = string.split(data,'=')
                        (dtype,content) = (sdata[0],string.join(sdata[1:],'='))
                    except: self.errorLog('Problem with detail: %s' % data); continue
                content = string.join(string.split(content))
                if content[:1] == '"': content = content[1:]
                if content[-1:] == '"': content = content[:-1]
                if dtype in self.list['DetailSkip']:
                    if dtype in ftdic: ftdic.pop(dtype)
                    continue
                elif dtype in ftdic and ftdic[dtype]: ftdic[dtype] = '%s; %s' % (ftdic[dtype],content)
                else: ftdic[dtype] = content
                if dtype not in self.list['Details']:
                    if dtype in nonsplit_features: ftdic['details'].append('/%s' % (dtype))
                    else: ftdic['details'].append('/%s="%s"' % (dtype,content))
            ftdic['details'] = string.join(ftdic['details'])
            ## ~ [1a] ~ Special Protein feature details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ftdic['feature'] == 'Protein':
                if 'protein_id' not in ftdic or not ftdic['protein_id']: ftdic['protein_id'] = locus
            ### ~ [2] ~ Add sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Make this a function ftdic['sequence'] = self.featureSequence(ftdict,locus)
            ftdic['sequence'] = self.featureSequence(ftdic,force=True)
            return ftdic
        except: self.errorLog('%s.featureDict error' % self,quitchoice=not self.debugging()); self.deBug(ftentry); self.deBug(ftdic)
#########################################################################################################################
    def featureSequence(self,ftdic,force=False): ### Adds sequence to ftdic using locus, position and self.dict['Sequence']
        '''
        Adds sequence to ftdic using locus, position and self.dict['Sequence'].
        >> ftdic:dict of feature information. This is updated with ftdic['sequence'], which is also returned.
        >> force:bool [False] = Whether to regenerate sequence even if found.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug(self.dict['Sequence'].keys())
            if 'sequence' in ftdic and ftdic['sequence'] and not force: return ftdic['sequence']
            try:
                ftstart = ftdic['start']
                ftend = ftdic['end']
            except: ftstart = ftend = 0
            ftdic['sequence'] = ''
            locus = ftdic['locus']
            if locus not in self.dict['Sequence']: ftlocus = string.split(locus,'_')[-1]
            else: ftlocus = locus
            ### ~ [1] Generate list of positions to stick together ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pos = ftdic['position'][0:]
            if pos.find('complement') == 0:
                complement = True
                pos = pos[len('complement('):-1]
            else: complement = False
            if pos.find('order') == 0: pos = pos[len('order('):-1]
            if pos.find('join') == 0: pos = pos[len('join('):-1]
            pos = string.split(pos,',')
            ### ~ [2] Process list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try:
                if '..' not in pos[0]:
                    if '^' in pos[0]: (ftdic['start'],ftdic['end']) = string.split(pos[0],'^')
                    elif pos[0][:1] in '<>': ftdic['start'] = ftdic['end'] = pos[0][1:]
                    else: ftdic['start'] = ftdic['end'] = '%d' % string.atoi(pos[0])
                else: (ftdic['start'],ftdic['end']) = string.split(pos[0],'..')
            except: self.debug('%s\n-> %s' % (ftdic['position'],pos)); raise
            if ftdic['start'][:1] in '<>': ftdic['start'] = ftdic['start'][1:] # Adjust for exons outside coding region
            if ftdic['end'][:1] in '<>': ftdic['end'] = ftdic['end'][1:]       # Adjust for exons outside coding region
            for frag in pos:
                if frag.find('complement') == 0:
                    fragcomp = True
                    frag = frag[len('complement('):-1]
                    if fragcomp: self.warnLog('Odd complementation for %s %s' % (ftdic['locus'],ftdic['position']))
                else: fragcomp = False #complement
                if '..' not in frag:
                    if '^' in frag: (start,end) = string.split(frag,'^')
                    elif frag[:1] in '<>': start = end = frag[1:]
                    else: start = end = '%d' % string.atoi(frag)
                else: (start,end) = string.split(frag,'..')
                if start[:1] in ['<','>']: start = start[1:]
                if end[:1] in ['>','>']: end = end[1:]
                try: string.atoi(end); ftdic['end'] = end
                except:
                    self.debug('%s\n-> %s' % (ftdic['position'],pos))
                    self.printLog('#ERR','Problem with locus %s position: %s' % (ftdic['locus'],ftdic['position']))
                    ftdic['end'] = rje.matchExp('(\d+)',end)[0]
                    self.printLog('#END','Corrected end position: %s -> %d' % (end,string.atoi(ftdic['end'])))
                if fragcomp: ftdic['sequence'] += rje_sequence.reverseComplement(self.dict['Sequence'][ftlocus][string.atoi(start)-1:string.atoi(ftdic['end'])])
                else: ftdic['sequence'] += self.dict['Sequence'][ftlocus][string.atoi(start)-1:string.atoi(ftdic['end'])]
            if complement: ftdic['sequence'] = rje_sequence.reverseComplement(ftdic['sequence'])
            #x# V1.3.1: Complement was always handled for each fragment - this would give sequences stuck together
            #X# in the wrong order!
            if 'join(' in ftdic['position']:
                #self.debug('%s => %s to %s = %d nt => %d nt sequence' % (ftdic['position'],ftdic['start'],ftdic['end'],int(ftdic['end'])-int(ftdic['start'])+1,len(ftdic['sequence'])))
                if int(ftdic['end'])-int(ftdic['start'])+1 == len(ftdic['sequence']): self.warnLog('Exonic join failure for %s' % ftdic['position'])
            ### ~ [3] Check and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ftstart:
                if int(ftstart) != int(ftdic['start']): self.warnLog('Feature start mismatch (%s %s -> %s vs %s)' % (locus,ftdic['position'],ftstart,ftdic['start']))
                ftdic['start'] = ftstart
            if ftend:
                if int(ftend) != int(ftdic['end']): self.warnLog('Feature start mismatch (%s %s -> %s vs %s)' % (locus,ftdic['position'],ftend,ftdic['end']))
                ftdic['end'] = ftend
            return ftdic['sequence']
        except: self.errorLog('%s.featureSequence error' % self,quitchoice=not self.debugging()); self.deBug(ftdic)
#########################################################################################################################
    def featurePos(self,ftdic,genpos): ### Returns position relative to feature. (Both positions 1-L, not 0<L)
        '''
        Returns position relative to feature. (Both positions 1-L, not 0<L).
        >> ftdic:dict of feature information. This is updated with ftdic['sequence'] if missing.
        >> genpos:int = Genomic position to convert to feature position (1-L).
        << ftpos:int = Position relative to ftdic['sequence'] (1-L). Will return 0 if outside or -1 for intron.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locus = ftdic['locus']
            if locus not in self.dict['Sequence']: ftlocus = string.split(locus,'_')[-1]
            else: ftlocus = locus
            ftseq = '' #self.featureSequence(ftdic)
            if ftdic['start'] < ftdic['end'] and (genpos < ftdic['start'] or genpos > ftdic['end']): return 0
            elif ftdic['start'] > ftdic['end'] and (ftdic['start'] > genpos > ftdic['end']): return 0    # Feature over split circle
            pos = ftdic['position'][0:]
            if pos.find('complement') == 0:
                complement = True
                pos = pos[len('complement('):-1]
            else: complement = False
            if pos.find('order') == 0: pos = pos[len('order('):-1]
            if pos.find('join') == 0: pos = pos[len('join('):-1]
            pos = string.split(pos,',')
            if 'join' in ftdic['position']: self.debug(pos)
            # NB. start and end are sorted by featurePos even for strange positions
            ### ~ [1] ~ Work through fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for frag in pos:
                if frag.find('complement') == 0:
                    fragcomp = True
                    frag = frag[len('complement('):-1]
                    if fragcomp: self.warnLog('Odd complementation for %s %s' % (ftdic['locus'],ftdic['position']))
                else: fragcomp = False #complement
                if '..' not in frag:
                    if '^' in frag: (start,end) = string.split(frag,'^')
                    elif frag[:1] in '<>': start = end = frag[1:]
                    else: start = end = '%d' % string.atoi(frag)
                else: (start,end) = string.split(frag,'..')
                if start[:1] in ['<','>']: start = start[1:]
                if end[:1] in ['>','>']: end = end[1:]
                try: end = string.atoi(end)
                except:
                    self.debug('%s\n-> %s' % (ftdic['position'],pos))
                    self.printLog('#ERR','Problem with locus %s position: %s' % (ftdic['locus'],ftdic['position']))
                    end = string.atoi(rje.matchExp('(\d+)',end)[0])
                ## Map position
                start = int(start); end = int(end)
                #!# might struggle with features spanning split circles. (Check and debug)
                self.debug('%s <= %s <= %s?' % (start,genpos,end))
                if genpos < start: return -1  # Intron
                if genpos > end:
                    if fragcomp: ftseq += rje_sequence.reverseComplement(self.dict['Sequence'][ftlocus][start-1:end])
                    else: ftseq += self.dict['Sequence'][ftlocus][start-1:end]
                    continue
                if fragcomp:
                    ftpos = end - genpos + 1 + len(ftseq)
                else:
                    ftpos = genpos - start + 1 + len(ftseq)
                if complement: ftpos = len(ftdic['sequence']) - ftpos + 1
                return ftpos
            #x# V1.3.1: Complement was always handled for each fragment - this would give sequences stuck together
            #X# in the wrong order!
            raise ValueError('Genomic Position %s failed to map to %s' % (genpos,ftdic['position']))
        except: self.errorLog('%s.featurePos error' % self,quitchoice=not self.debugging()); return 0
#########################################################################################################################
    def loadFeatures(self,expect=False):     ### Load features and loci into database tables
        '''
        Load features and loci into database tables.
        >> Expect:bool [False] = Whether to expect files and raise error (True) or just None (False)
        '''
        try:### ~ [1] Full entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locdb = self.db().addTable(name='Locus',expect=expect)
            ftdb = self.db().addTable(name='Feature',expect=expect,mainkeys=['locus','feature','position'])
            return locdb and ftdb
        except: self.errorLog('%s.featureDict error' % self)
#########################################################################################################################
    def saveFasta(self,logskip=True,bylocus=None):    ### Saves sequences and features in fasta format
        '''Saves sequences and features in fasta format.'''
        try:### ~ [1] Full entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locdb = self.db('Locus'); ftdb = self.db('Feature')
            if bylocus == None: bylocus = self.getBool('LocusOut')
            if 'full' in self.list['FasOut']:
                fasfile = self.getStr('Basefile') + '.full.fas'; rje.backup(self,fasfile,appendable=False)
                if not bylocus: FASOUT = open(fasfile,'w'); fx = 0
                for locus in rje.sortKeys(self.dict['Sequence']):
                    entry = locdb.data(locus)
                    code = entry['spcode']
                    name = 'gb_%s__%s %s|gi:%s %s' % (code,entry['accession'],locus,entry['gi'],entry['definition'])
                    addorg = '[%s]' % entry['organism']
                    if addorg not in name: name = '%s %s' % (name,addorg)
                    if bylocus:
                        fasfile = '%s%s.full.fas' % (self.getStr('LocusDir'),locus)
                        if os.path.exists(fasfile) and not self.force():
                            if logskip: self.printLog('#SKIP','Skipping %s output (force=F)' % fasfile)
                            continue
                        FASOUT = open(fasfile,'w'); fx = 0
                    FASOUT.write('>%s\n%s\n' % (name,self.dict['Sequence'][locus])); fx += 1
                    if bylocus:
                        FASOUT.close()
                        self.printLog('#FAS','%s sequences output to %s' % (rje.iStr(fx),fasfile))
                if not bylocus:
                    FASOUT.close()
                    self.printLog('#FAS','%s sequences output to %s' % (rje.iStr(fx),fasfile))
            ### ~ [2] Save features as sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seqtype in self.list['FasOut']:
                if seqtype == 'full': continue
                if seqtype == 'prot':
                    if 'Protein' in ftdb.index('feature'):
                        #if 'CDS' in ftdb.index('feature'): self.warnLog('"Protein" and "CDS" features found! Mixed sequence types? Will only output "Protein" features to *.prot.fas')
                        ftype = 'Protein'
                    else: ftype = 'CDS'
                    self.printLog('#PROT','%s feature type used for *.prot.fas output.' % ftype)
                else: ftype = seqtype
                fasfile = self.getStr('Basefile') + '.%s.fas' % seqtype; rje.backup(self,fasfile,appendable=False)
                if not ftdb.indexEntries('feature',ftype):
                    self.printLog('#FAS','No sequences to output to %s.' % (fasfile))
                    continue
                fx = 0
                if not bylocus: FASOUT = open(fasfile,'w')
                elif self.force():
                    loc2save = ftdb.indexKeys('locus')
                    for locus in loc2save:
                        fasfile = '%s%s.%s.fas' % (self.getStr('LocusDir'),locus,seqtype)
                        rje.backup(self,fasfile,appendable=False)
                else:
                    loc2save = []
                    for locus in ftdb.indexKeys('locus'):
                        fasfile = '%s%s.%s.fas' % (self.getStr('LocusDir'),locus,seqtype)
                        if os.path.exists(fasfile):
                            if logskip: self.printLog('#SKIP','Skipping %s output (force=F)' % fasfile)
                        else: loc2save.append(locus)
                for ftentry in ftdb.indexEntries('feature',ftype):
                    locus = ftentry['locus']
                    if bylocus:
                        if locus not in loc2save: continue
                        fasfile = '%s%s.%s.fas' % (self.getStr('LocusDir'),locus,seqtype)
                        FASOUT = open(fasfile,'a')
                    ftdic = self.featureDict(ftentry)
                    try:
                        if ftentry['gene'] == 'COX3': self.debug(ftentry); self.debug(ftdic)
                    except: pass
                    #self.bugPrint(ftentry['details'])
                    #self.debug(ftdic['details'])
                    transl = rje.matchExp('/transl_table="?(\d+)"?',ftdic['details'])
                    if transl: transl = transl[0]; self.bugPrint('%s %s => /transl_table="%s"' % (locus,ftdic['position'],transl))
                    else:
                        transl = '1'
                        if 'transl_table' in ftdic['details']: self.debug(ftdic['details'])
                        if 'transl_table' in ftentry['details']: self.debug(ftentry['details'])
                    try:
                        name = '%s_%s' % (string.split(seqtype.lower(),'_')[0],ftdic['spcode'])
                        if seqtype in ['prot','CDS'] and self.getStr('ProtAcc') in ftdic and ftdic[self.getStr('ProtAcc')]:
                            accnum = ftdic[self.getStr('ProtAcc')]
                        elif seqtype == 'gene' and self.getStr('GeneAcc') in ftdic and ftdic[self.getStr('GeneAcc')]:
                            accnum = ftdic[self.getStr('GeneAcc')]
                        elif seqtype == 'mobile_element':
                            accmatch = rje.matchExp('/note="(\S+);',ftdic['details'])
                            if accmatch: accnum = accmatch[0]
                            else: accnum = 'TE%s' % rje.preZero(fx+1,len(ftdb.indexEntries('feature',ftype)))
                        elif seqtype in ['mobile_element','LTR']:
                            accmatch = rje.matchExp('/note="(\S+)\s.+/db_xref="SGD:(\S+)"',ftdic['details'])
                            if accmatch: accnum = '%s.%s' % accmatch
                            else: accnum = '%s%s' % ({'mobile_element':'TE','LTR':'LTR'}[seqtype],rje.preZero(fx+1,len(ftdb.indexEntries('feature',ftype))))
                        else: accnum = ftdic['locus_tag']
                        name = '%s__%s' % (name,accnum)
                    except: self.errorLog('Cannot output %s %s %s' % (ftentry['locus'],ftentry['feature'],ftentry['position'])); continue
                    if transl != '1': self.bugPrint('%s -> transl_table %s' % (name,transl))
                    for data in ['pseudo','protein_id','db_xref','product']:
                        if data in ftdic and ftdic[data] != accnum: name = '%s %s' % (name,ftdic[data])
                    addorg = '[%s]' % ftdic['organism']
                    if addorg not in name: name = '%s %s' % (name,addorg)
                    if seqtype == 'prot' and ftype == 'CDS':
                        try: sequence = string.join(ftdic['translation'],'')
                        except: sequence = rje_sequence.dna2prot(ftdic['sequence'],transl=transl,warnobj=self)
                        if self.test() and transl != '1':
                            translseq = rje_sequence.dna2prot(ftdic['sequence'],transl=transl,warnobj=self)
                            if translseq == rje_sequence.dna2prot(ftdic['sequence']):
                                self.warnLog('/transl_table="%s" does nothing for %s?' % (transl,name))
                            if translseq != sequence: self.warnLog('%s seems to have ignored /transl_table="%s"?' % (name,transl))
                        if '*' in sequence[:-1]: self.warnLog('Internal STOP codon in %s: possible alternative genetic code?' % name,'altcode',suppress=True)
                    else: sequence =  ftdic['sequence']
                    if sequence[-1:] == '*': sequence = sequence[:-1]
                    FASOUT.write('>%s\n%s\n' % (name,sequence)); fx += 1
                    if bylocus: FASOUT.close()
                if not bylocus:
                    FASOUT.close()
                    self.printLog('#FAS','%s sequences output to %s' % (rje.iStr(fx),fasfile))
                else: self.printLog('#FAS','%s sequences output to %s locus files' % (rje.iStr(fx),rje.iLen(loc2save)))
        except: self.errorLog('%s.saveFasta error' % self)
#########################################################################################################################
### End of SECTION II: GenBank Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def setupRefGenome(callobj):    ### This is a generic reference genome setup method for callobj.
    '''
    This is a generic reference genome setup method for callobj.
    >> callobj should have str['RefGenome']. This could be a fasta file, genbank file (*.gb or *.gbk) or basefile.
    << callobj will get str['RefBase'], str['RefGenome'], bool['Features'] and str['FTFile'] set.
    '''
    self = callobj
    try:### ~ [1] Setup Reference Genome with Genbank extraction if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [1a] Reference Basename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        gbfile = genbankFile(self.getStr('RefGenome'))
        if gbfile:
            self.setStr({'RefBase':rje.baseFile(gbfile)})
            self.printLog('#GBK','GenBank file found: %s' % gbfile)
        else: self.setStr({'RefBase':rje.baseFile(self.getStr('RefGenome'))})
        # >>> #
        ## ~ [1b] Determine whether GenBank Parsing required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        checkfiles = ['.full.fas','.gene.fas','.prot.fas','.Feature.tdt']
        rungb = gbfile and not rje.checkForFiles(filelist=checkfiles,basename=self.getStr('RefBase'),log=self.log,cutshort=False,ioerror=False,missingtext='')
        self.debug('Run Genbank: %s' % rungb)
        ## ~ [1c] GenBank Parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if rungb:
            gcmd = ['protacc=locus_tag','details=product,gene_synonym,note,db_xref']   # Defaults
            gcmd += self.cmd_list   # Can over-ride/add. This include spcode=X
            gcmd += ['seqin=%s.gb' % self.getStr('RefBase'),'taxdir=','tabout=T','fasout=full,gene,cds,prot']
            GenBank(self.log,gcmd).run()
            rje.checkForFiles(filelist=checkfiles,basename=self.getStr('RefBase'),log=self.log,cutshort=False,ioerror=True,missingtext=' Genbank parsing failed?')
        ## ~ [1d] Features file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ftfile = '%s.Feature.tdt' % self.getStr('RefBase')
        self.setBool({'Features':os.path.exists(ftfile)})   # Whether features table generated
        if not rje.exists(self.getStr('FTFile')):
            self.setStr({'FTFile':ftfile})
            self.printLog('#FTFILE','FTFile set: %s' % ftfile)
        self.debug('Features: %s' % self.getBool('Features'))
        ## ~ [1e] RefGenome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not rje.exists(self.getStr('RefGenome')) or self.getStr('RefGenome') == gbfile:
            self.setStr({'RefGenome':'%s.full.fas' % self.getStr('RefBase')})
            self.printLog('#SETREF','RefGenome: %s.full.fas' % self.getStr('RefBase'))
        if not rje.exists(self.getStr('RefGenome')): raise IOError('Cannot find RefGenome: %s!' % self.getStr('RefGenome'))
    except: self.errorLog('%s setupRefGenome() error.' % self.prog())
#########################################################################################################################
def genbankFile(filename):  ### Looks for associated *.gb or *.gbk file and returns, else None.
    '''
    Looks for associated *.gb or *.gbk file and returns, else None.
    >> filename:basefile or filename for which extension will be removed.
    << genbank file name or None if not found.
    '''
    files = [filename,rje.baseFile(filename)]
    exts = ['gb','gbk']
    for gbase in files:
        for ext in exts:
            gfile = '%s.%s' % (gbase,ext)
            if rje.exists(gfile): return gfile
    return None
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
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: GenBank(mainlog,cmd_list).run()

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
