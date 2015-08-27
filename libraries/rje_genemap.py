#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_genemap
Description:  RJE Gene & Database ID Mapping Module
Version:      1.5
Last Edit:    16/12/13
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to replace rje_genecards, which has become a bit unwieldy since its conception. Some of the
    original functions of rje_genecards will still be maintained by (hopefully) in a simplified format. Some of the
    additional mapping functions of Pingu will be added to this module for easier implementation across packages.

    The main functions of this module are:
    1. To map, store and retrieve database cross-references for a key dataset of gene IDs, usually HGNC symbols.
    2. To store a number of aliases for the key gene IDs, including old versions of accession numbers etc.
    3. To retrieve sequences from given datasets for stored aliases/genes.

    The main processing pipeline is as follows:
    1. Read in key data and generate data structure OR load pickle.
    2. Repeat 1 until all data/pickles integrated.
    3. Save pickle and output new data flat files, if desired.

    This is the limit of the standalone functionality of the program. However, the GeneMap class with have a number of
    additional methods for data retrieval by other programs that use it.

GeneMap Class:
    The GeneMap Class stores two main data dictionaries:
    1. A dictionary for the key Gene IDs that contains mappings to other databases.
    2. A dictionary that maps aliases onto other Gene IDs.

    In addition, sequence files may be loaded and used to map IDs onto Sequence objects.    

Key Input:
    There are four primary input files that are processed into the mapping:
    1. Designed for human data, an HGNC download file is one of the key input files. Headers will be converted into those
    from the original rje_genecards files, which are now replaced by sourcedata=FILE.
    2. Source Data files are delimited text files containing mapping to various databases. In each case, the first column
    should be unique for each line. This will be treated as an Alias. If a Symbol column is found, this will be treated
    as a key identifier (unless keyid=X has been changed).
    3. Alias files containing simple lists of ID:Alias to populate Alias dictionary. The first column can have any header
    but must be the identifier to map *to*. Another column must have the header "Aliases" and be a comma-separated list
    of aliases.
    4. Pickle data containing a pickled GeneMap object.

    In addition, sequence files may be loaded that have additional links and can be used to map to sequences. These are:
    1. EnsLoci = This is used to add additional mapping of genes to proteins and to EnsLoci protein IDs.

Input Processing:
    As data is loaded, either from a pickle or a text file, its data is integrated. NOTE: These commands can be repeated
    several times and, unlike normal, subsequent commands will not replace earlier ones but simply add to the list. If
    there is danger of additional unwanted commands in the command argument list, then the loadData() method should be
    called with a specified list of commandline options rather than using the default system arguments.

    If a given set of data has a "Symbol" (KeyID) Header then this is added to the main Data dictionary as a key and all
    column headers as stored data. (These column headers are stored in the "Header" list.) The Alias - the original key
    of the dictionary - is added to self.dict['Alias']. Note that each Alias can be involved in many-to-many and circular
    referencing, which will need to be dealt with by the class when mapping. And headers that are in the "XRef" list of
    database cross-references will also be added to the Alias dictionary. If there is no KeyID, the data is stored in a
    "TempData" dictionary, and XRef headers are aliased to the Alias rather than the KeyID.

    If a KeyID already exists in the Data dictionary, then any blank entries will be overwritten but data loaded from a
    previous file will not be. All aliases will be mapped to the KeyID, however, even if they do not end up in the Data
    dictionary itself. If a KeyID is missing from the Data dictionary but present in the TempData dictionary, it will be
    overwritten in the same way and moved to the Data dictionary. If an Alias without a KeyID is already present in the
    TempData dictionary, the same will happen without any transferral.

    Sequence data will be processed according to the specifics of the type of sequence file it is. EnsLoci sequences
    will be converted into an EnsLoci dictionary of {ID:Sequence} but also key gene-protein mappings will be extracted.
    The protein to gene aliases will be added to the Alias dictionary, while the EnsLoci ID will be added as an XRef to
    the appropriate Data or TempData dictionary element.

    Once all data has been read in, each TempData entry will be assessed using the Alias mappings to see if it, or any of
    its XRef entries, is an alias for a KeyID. If so, its data will be combined with that in Data and it will be removed
    from TempData. If not, it will be assessed for being an alias of another TempData entry and will be combined if so.
    After this final stage of processing, any entries still in TempData will be promoted to KeyIDs and added to the main
    dictionary, though they will not appear in the "KeyID" list.

Commandline:
    ### ~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hgncdata=FILE   : Download file containing HGNC data. []
    mgidata=FILE    : Download file containing MGI data (ftp://ftp.informatics.jax.org/pub/reports/MGI_MouseHumanSequence.rpt) []
    sourcedata=FILE : File containing data in order of preference regarding conflicting data. []
    aliases=FILE    : Files containing aliases only. []
    pickledata=FILE : Genemap pickle to import and use. []
    ensloci=FILE    : File of EnsLoci genome to incorporate [None]
    genepickle=FILE : Use pickle of GeneMap data without additional loading/processing etc. [None]
    pfamdata=FILE   : Delimited files containing domain organisation of sequences [None]
    approved=LIST   : Approved HGNC gene symbols to avoid over-zealous alias mapping (will add to from HGNC) []

    ### ~ Processing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    keyid=X         : Key field header to be used in main Data dictionary - aliases map to this [Symbol]
    xref=LIST       : Headers in Data dictionaries that are used for aliases [EnsEMBL,Entrez,HGNC,HPRD,UniProt]
    useweb=T/F      : Whether to try and extract missing data from GeneCards website [False]
    skiplist=LIST   : Skip genes matching LIST when using GeneCards website (e.g. XP_*) ['HPRD*']
    
    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : Root for output files [genemap]
    flatout=T/F     : Whether to output flatfiles (*.data.tdt & *.aliases.tdt) [False]
    pickleout=T/F   : Whether to output pickle (*.pickle.gz) [False]

Uses general modules: copy, glob, pickle, os, string, sys, time, urllib2
Uses RJE modules: rje, rje_seq, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, pickle, os, string, sys, time, urllib2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based loosely on rje_genecards V0.4 (28-Mar-08).
    # 1.0 - Standalone working version with basic functions.
    # 1.1 - Added bestMap() function for better compatibility with PINGU (V3.0)
    # 1.2 - Fixed bug of mapping current Approved Gene Symbols to other Gene Symbols due to redundant Aliases.
    # 1.3 - Add reduction of data to gene list.
    # 1.4 - Modified to read in MOUSE data.
    # 1.5 - Minor tweak of expected HGNC input following change to downloads.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Could add a ForceWeb option, which would parseCard for ALL entries and OVERWRITE data.
    # [ ] : Add tracking of source files?
    # [ ] : Fix bug that incorrectly maps aliases when the alias is a valid gene symbol in its own right. (Remove these!)
    # [ ] : Modify to make appropriate use of rje_db
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_GeneMap', '1.5', 'December 2013', '2008')
    description = 'RJE Gene & Database ID Mapping Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.'] #,rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        cmdhelp = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmdhelp > 0:
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
badalias = ['!FAILED!','','NONE','-','N/A']     # Entries to most definitely ignore!
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GeneMap Class                                                                                           #
#########################################################################################################################
class GeneMap(rje.RJE_Object):     
    '''
    GeneMap Class. Author: Rich Edwards (2008).

    Info:str
    - EnsLoci = File of EnsLoci genome to incorporate [None]
    - PFamData = Delimited files containing domain organisation of sequences [None]
    - GenePickle = Use pickle of GeneMap data without additional loading/processing etc. [None]
    - KeyID = Key field header to be used in main Data dictionary - aliases map to this [Symbol]

    Opt:boolean
    - FlatOut = Whether to output flatfiles (*.data.tdt & *.aliases.tdt) [False]
    - PickleOut = Whether to output pickle (*.pickle.gz) [False]
    - UseWeb = Whether to try and extract missing data from GeneCards website [False]

    Stat:numeric

    List:list
    - Approved = Full list of Approved HGNC gene symbols (from HGNC) to avoid over-zealous alias mapping
    - Headers = Full list of headers in Data dictionaries []
    - KeyID = List of KeyIDs from original data, i.e. not promoted from TempData []
    - SkipList = Skip genes matching LIST when using GeneCards website (e.g. XP_*) [HPRD*]
    - WebCheck = List of IDs already checked against website (do avoid repetition) []
    - XRef = List of headers in Data dictionaries that are used for aliases ['EnsEMBL','Entrez','HGNC','HPRD','UniProt']

    Dict:dictionary
    - Alias = Dictionary of {Alias:["Better" Aliases]}
    - BestMap = Dictionary of {id:self.bestMap(id)} to save time during many bestMap() calls
    - Data = Main data dictionary for KeyIDs
    - EnsLoci = Dictionary of {EnsLoci ID:Sequence object}
    - GeneMap = Dictionary of successful gene mappings to accelerate Alias mapping {Alias:[KeyIDs]}
    - PFam = Dictionary of {EnsLoci ID:[PFam domains]}
    - Redundancy = Dictionary of {BestID:[AltIDs]} generated from bestMap() method
    - TempData = Temporary dictionary during data processing.

    Obj:RJE_Objects
    - DB = Database object
    - EnsLoci = SeqList containing EnsLoci data (just the last one added)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['EnsLoci','KeyID','GenePickle','PFamData','AddAlias']
        self.optlist = ['FlatOut','PickleOut','UseWeb']
        self.statlist = []
        self.listlist = ['Headers','KeyID','SkipList','WebCheck','XRef','Approved']
        self.dictlist = ['Alias','EnsLoci','Data','TempData','PFam','Redundancy','BestMap']
        self.objlist = ['DB']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'KeyID':'Symbol','Basefile':'genemap'})
        self.setOpt({'FlatOut':False,'PickleOut':False})
        self.list['XRef'] = ['EnsEMBL','Entrez','HGNC','HPRD','UniProt']
        self.list['SkipList'] = ['HPRD*']
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
                self._cmdReadList(cmd,'info',['KeyID'])
                self._cmdReadList(cmd,'file',['GenePickle','PFamData'])
                self._cmdReadList(cmd,'opt',['FlatOut','PickleOut','UseWeb'])
                self._cmdReadList(cmd,'list',['XRef','SkipList','Approved'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self.unPickle()
#########################################################################################################################
    def run(self):  ### Main standalone run method
        '''Main standalone run method.'''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['GenePickle']: return self.log.printLog('#RUN','No additional processing following unpickling of %s' % self.info['GenePickle'])
            self.loadData(self.cmd_list)
            self.deBug(self.list['Headers'])
            ### ~ [2] Save Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['FlatOut']: self.flatOut()
            if self.opt['PickleOut']: self.pickleOut()
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def unPickle(self):     ### Attempts to replace self with pickled version
        '''Attempts to replace self with pickled version.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sourcefile = self.info['GenePickle']
            if sourcefile.lower() in ['','none']:
                self.info['GenePickle'] = ''
                return
            if not os.path.exists(sourcefile):
                self.info['GenePickle'] = ''
                return self.log.errorLog('Pickle file "%s" not found' % (sourcefile),printerror=False)
            ### ~ [2] ~ Deal with Zip and load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gunzip = sourcefile[-3:] == '.gz' and not self.opt['Win32']
            if gunzip:
                self.gUnzip(sourcefile)
                sourcefile = sourcefile[:-3]
            try:
                self.log.printLog('#LOAD','Attempting to load GeneMap pickle.',log=False)
                mypickle = loadGeneMapPickle(sourcefile) 
                self.log.printLog('#LOAD','Pickle loaded: %s.' % (sourcefile))
            except: return self.log.errorLog('Problem loading %s' % sourcefile)
            if gunzip: self.gZip(sourcefile)
            ### ~ [3] Process Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mypickle.log = self.log
            for obj in mypickle.obj:
                try: obj.log = self.log
                except: pass
            mypickle.info['GenePickle'] = sourcefile
            self.setInfo(mypickle.info)
            self.setStat(mypickle.stat)
            self.setList(mypickle.list)
            if len(self.list['Headers']) < 7:
                for h in ['Symbol','Desc','Entrez','OMIM','UniProt','EnsEMBL','HGNC','EnsLoci']:
                    if h not in self.list['Headers']: self.list['Headers'].append(h)
            self.setDict(mypickle.dict)
            #x#self.setOpt(mypickle.opt)
            self.setObj(mypickle.obj)
            #x#self = mypickle
            self.log.printLog('#PICKLE','%s: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <2> ### Data Extraction Methods                                                                                 #
#########################################################################################################################
    def getEnsLoci(self,id):    ### Returns list of sequences for given ID
        '''Returns list of sequences for given ID.'''
        kid = self.getKeyID(id)
        ens = []
        for k in kid:
            try: seq = self.dict['EnsLoci'][self.dict['Data'][k]['EnsLoci']]
            except: continue
            if seq not in ens: ens.append(seq)
        return ens
#########################################################################################################################
    def getGeneData(self,id):   ### Returns dictionary for given ID
        '''Returns dictionary for given ID. Maps aliases if necessary. Returns first keyID if multiples.'''
        kid = self.getKeyID(id)
        if not kid: return {}
        for k in kid:
            if k in self.list['KeyID']: return self.dict['Data'][k]
        return self.dict['Data'][kid[0]]
#########################################################################################################################
    def getKeyID(self,id,restricted=False):  ### Maps id onto KeyIDs if possible and returns keys or empty list
        '''
        Maps id onto KeyIDs if possible and returns keys or empty list.
        >> restricted:bool = whether to restrict to "original" KeyIDs, in self.list['KeyID']
        '''
        try:### ~ [1] ~ Check ID straight against KeyIDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            id = id.upper()
            if id in badalias: return []
            if id in self.dict['Data'] and (id in self.list['KeyID'] or not restricted): return [id] 
            ### ~ [2] ~ Look for aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            keyids = []
            for k in self.mapGene(id):
                if k in self.dict['Data'] and (k in self.list['KeyID'] or not restricted): keyids.append(k)
            return keyids
        except: self.errorLog(rje_zen.Zen().wisdom())
        return []
#########################################################################################################################
    ### <3> ### Alias/Mapping Methods                                                                                   #
#########################################################################################################################
    def addAlias(self,id,alias):  ### Add id to self.dict['Alias'][alias]
        '''Add id to self.dict['Alias'][alias].'''
        try:
            if not alias or not id: return
            alias = alias.upper()
            id = id.upper()
            if alias == id or alias in badalias or id in badalias: return
            if alias not in self.dict['Alias']: self.dict['Alias'][alias] = []
            if id not in self.dict['Alias'][alias]: self.dict['Alias'][alias].append(id)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def redundancy(self,id):    ### Returns redundancy identified during bestMap()                                  |1.1|
        '''Returns redundancy identified during bestMap().'''
        try: return self.dict['Redundancy'][id]
        except: return []
#########################################################################################################################
    def bestMap(self,id,addaliases=False,usedict=True):   ### Returns single best ID from mapGene                   |1.1|
        '''
        Returns single best ID from mapGene.
        >> id:str = Gene/Protein ID to map using Aliases.
        >> addaliases:bool [False] = whether to add other legitimate IDs as aliases for BestGene.
        >> usedict:bool [True] = Whether to read/write from self.dict['BestMap'] to save time.
        '''
        ### ~ [1] ~ Select best gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if usedict and id in self.dict['BestMap']: return self.dict['BestMap'][id]
        mid = self.mapGene(id)
        rid = self.getKeyID(id,restricted=True)
        kid = self.getKeyID(id,restricted=False)
        bestid = mid[0]
        for alt in mid[1:]:
            if alt in rid and bestid not in rid: bestid = alt
            elif alt in kid and bestid not in kid: bestid = alt
            elif len(alt) < len(bestid): bestid = alt     #!# Add additional criteria with time (pseudogenes?)
        ### ~ [2] ~ Update dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if bestid not in self.dict['Redundancy']: self.dict['Redundancy'][bestid] = []
        for alt in mid:
            if alt == id: continue
            if addaliases: self.addAlias(bestid,alt)
            if alt not in self.dict['Redundancy'][bestid]: self.dict['Redundancy'][bestid].append(alt)
        if usedict: self.dict['BestMap'][id] = bestid
        return bestid
#########################################################################################################################
    def mapGene(self,id):   ### Maps an ID onto aliases with entries in Data dictionaries (if possible)             |1.1|
        '''Maps an ID onto aliases with entries in Data dictionaries (if possible).'''
        try:
            id = id.upper()
            if id in badalias: return []
            aliases = self.aliases(id)
            mapgene = []
            for a in aliases:
                if a in self.dict['Data'] or a in self.dict['TempData']: mapgene.append(a)
            #x#if mapgene: self.deBug('%s: %s' % (id,mapgene))
            if mapgene: return mapgene
            else: return [id]
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def aliases(self,id,full=True):   ### Returns full list of aliases for a given gene                                       |1.2|
        '''Returns full list of aliases for a given gene.'''
        try:### ~ [1] ~ Cycle through adding all aliases, and aliases of aliases etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            id = id.upper()
            if id in badalias: return []
            aliases = [id]
            process = [id]
            while process:                              # New alias to add aliases of...
                id = process.pop(0)                     # Reduce process list, so each looked at just once
                if id in self.list['Approved']: continue# Do not map extra aliases of Approved Gene Symbols
                if id in self.dict['Alias']:            # May have no aliases
                    for a in self.dict['Alias'][id]:    # Look at each alias in turn
                        if a not in aliases:            # New - add to lists
                            aliases.append(a)
                            process.append(a)
            if not full:
                for a in aliases[0:]:
                    if a in self.dict['Alias'] and rje.listIntersect(aliases,self.dict['Alias'][a]): aliases.remove(a)
            return aliases
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def reverseAliases(self):   ### Returns a dictionary with all the aliases pointing to each ID
        '''Returns a dictionary with all the aliases pointing to each ID.'''
        rev = {}; ax = 0.0; atot = len(self.dict['Alias'])
        for alias in self.dict['Alias']:
            self.progLog('\r#REV','Reversing aliases: %.2f%%' % (ax/atot)); ax += 100.0
            for id in self.aliases(alias)[1:]:  # Not self!
                if id not in rev: rev[id] = [alias]
                else: rev[id].append(alias)
        return rev
#########################################################################################################################
    def makeGeneMap(self):  ### Makes CardMap dictionary For compatibility with old GeneCards
        '''Makes CardMap dictionary For compatibility with old GeneCards. NOTE: Currently only one mapping per alias.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['CardMap'] = {}
            mtxt = 'Generating CardMap dictionary for GeneCards compatibility:'
            (ax,atot) = (0.0,len(self.dict['Alias']))

            ### ~ [2] ~ Generate dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for alias in self.dict['Alias']:
                self.log.printLog('\r#MAP','%s %.1f%%' % (mtxt,ax/atot),newline=False,log=False)
                ax += 100.0
                if alias in self.dict['CardMap']: continue
                kid = self.getKeyID(alias,restricted=True)
                if kid:
                    self.dict['CardMap'][alias] = kid[0]        #!# Consider an ambiguous version #!#
                    continue
                kid = self.getKeyID(alias)
                if kid:
                    self.dict['CardMap'][alias] = kid[0]        #!# Consider an ambiguous version #!#
            self.log.printLog('\r#MAP','Generation of CardMap dictionary for GeneCards compatibility complete.')
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### Input loading/processing Methods                                                                        #
#########################################################################################################################
    def loadData(self,cmd_list):    ### Loads data into self and processes.                                         |1.2|
        '''
        Loads data into self and processes.
        >> cmd_list = List of commandline arguments to use to determine source data.
        hgncdata=FILE   : Download file containing HGNC data. []
        mgidata=FILE    : Download file containing mouse MGI data. []
        sourcedata=FILE : File containing data in order of preference regarding conflicting data. []
        aliases=FILE    : Files containing aliases only. []
        pickledata=FILE : Genemap pickle to import and use. []
        ensloci=FILE    : File of EnsLoci genome to incorporate [None]
        '''
        try:### ~ [1] ~ Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cmd in cmd_list:
                try:
                    if cmd.find('hgncdata=') == 0: self.loadHGNC(cmd[cmd.find('=')+1:])
                    if cmd.find('mgidata=') == 0: self.loadMGI(cmd[cmd.find('=')+1:])
                    if cmd.find('sourcedata=') == 0: self.loadSource(cmd[cmd.find('=')+1:])
                    if cmd.find('aliases=') == 0: self.loadAlias(cmd[cmd.find('=')+1:])
                    if cmd.find('pickledata=') == 0: self.loadPickle(cmd[cmd.find('=')+1:])
                    if cmd.find('ensmap=') == 0: self.loadEnsMap(cmd[cmd.find('=')+1:])
                    if cmd.find('ensloci=') == 0: self.loadEnsLoci(cmd[cmd.find('=')+1:])
                    if cmd.find('addalias=') == 0: self.addAliasFile(cmd[cmd.find('=')+1:])
                except: self.log.errorLog('Problem with cmd:%s' % cmd)

            ### ~ [2] ~ Check TempData for Alias of KeyID or other Alias ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.checkTempAliases()
                    
            ### ~ [3] ~ Look up missing entried at GeneCards website ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['UseWeb']:
                (gx,fx,tx) = (0,0,len(self.dict['TempData']))
                try:
                    for temp in rje.sortKeys(self.dict['TempData']):
                        if self.parseCard(temp): gx += 1
                        else: fx += 1
                        ctxt = 'Parsing GeneCards for %s Temp IDs: %.1f%%;' % (rje.integerString(tx),100.0*(gx+fx)/tx)
                        self.progLog('\r#CARD','%s %d parsed; %d failed.' % (ctxt,gx,fx))
                        self.verbose(1,4,' [%s]        ' % temp,0)
                    self.printLog('\r#CARD','Parsing GeneCards for %s Temp IDs complete: %s parsed; %s failed.' % (rje.integerString(tx),rje.integerString(gx),rje.integerString(fx)))
                except KeyboardInterrupt: self.log.printLog('\r#CARD','Parsing GeneCards for %s genes stopped: %s parsed; %s failed.' % (rje.integerString(tx),rje.integerString(gx),rje.integerString(fx)))
                except: raise
                self.checkTempAliases()     # Check again following update

            ### ~ [4] ~ Finally, upgrade TempData to KeyIDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (tx,tot) = (0.0,len(self.dict['TempData']))
            for temp in rje.sortKeys(self.dict['TempData']):
                self.progLog('\r#DATA','Moving temp IDs to main Data: %.1f%%' % (tx/tot))
                tx += 100.0
                self.dict['Data'][temp] = self.dict['TempData'].pop(temp)
            self.printLog('\r#DATA','Temp IDs moved: %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def checkTempAliases(self):     ### Check TempData for Alias of KeyID or other Alias                            |1.2|
        '''Check TempData for Alias of KeyID or other Alias.'''
        try:### ~ [1] ~ Check TempData for Alias of KeyID or other Alias ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (tx,tot) = (0.0,len(self.dict['TempData']))
            for temp in rje.sortKeys(self.dict['TempData']):
                self.progLog('\r#TEMP','Looking for aliases of temp IDs: %.1f%%' % (tx/tot))
                tx += 100.0
                ## ~ [2a] ~ Check for Alias to KeyID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                kids = self.getKeyID(temp)
                try: tdata = self.dict['TempData'][temp]
                except:
                    self.printLog('\r#TEMP','Gene "%s" already gone from TempData! (Aliases: %s)' % (temp,string.join(self.aliases(temp,full=False))))
                    continue
                for xref in self.list['XRef']:
                    if xref in tdata and tdata[xref]: kids += self.getKeyID(tdata[xref])
                kids = rje.sortUnique(kids,xreplace=False)
                if kids:
                    for k in kids: self.addData(k,tdata)     # Add to aliased KeyIDs                  
                    #self.deBug('%s KeyID: %s' % (temp,kids))
                ## ~ [2b] ~ Check for alias to TempData ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:
                    kids = self.aliases(temp)[1:]    # No self
                    for k in kids:
                        if k in self.dict['TempData']: self.dict['TempData'][k] = rje.combineDict(self.dict['TempData'][k],self.dict['TempData'][temp],overwrite=False)
                    #self.deBug('%s Temp: %s' % (temp,kids))
                if kids and temp in self.dict['TempData']: self.dict['TempData'].pop(temp)
            self.printLog('\r#TEMP','Aliases of temp ID mapped: %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadHGNC(self,sourcefile):  ### Loads HGNC data                                                             |1.2|
        '''
        Loads HGNC data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Read into dictionary with HGNC ID as key ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.log.errorLog('HGNC file "%s" not found' % (sourcefile),printerror=False)
            hgncdata = rje.dataDict(self,sourcefile,['HGNC ID'])
            aliaii = {}     # Dictionary of withdrawn symbols to map
            for h in ['Symbol','Desc','Entrez','OMIM','UniProt','EnsEMBL','HGNC']:
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            ### ~ [2] Parse out information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot) = (0.0,len(hgncdata))
            for hgnc in rje.sortKeys(hgncdata):
                self.progLog('\r#HGNC','Processing HGNC: %.1f%%' % (hx/htot))
                hx += 100.0
                ## ~ [2a] Adjust headers for new vs old HGNC compatibility ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                data = hgncdata[hgnc]
                if 'Synonyms' in data and 'Aliases' not in data: data['Aliases'] = data.pop('Synonyms')
                for hkey in rje.sortKeys(data):
                    if rje.matchExp('^(\S.+\S)\s*\(mapped data supplied by \S+\)',hkey):
                        data['%s (mapped data)' % rje.matchExp('^(\S.+\S)\s*\(mapped data supplied by \S+\)',hkey)[0]] = data.pop(hkey)
                    if rje.matchExp('^(\S.+\S)\s*\(supplied by \S+\)',hkey):
                        data['%s (mapped data)' % rje.matchExp('^(\S.+\S)\s*\(supplied by \S+\)',hkey)[0]] = data.pop(hkey)
                gene = data['Approved Symbol'].upper()
                ## ~ [2b] Bypass obselete symbols ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if gene.find('~WITHDRAWN') > 0: continue
                if gene not in self.list['Approved']: self.list['Approved'].append(gene)
                ## ~ [2c] Add aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for alias in string.split(data['Aliases'],', '): self.addAlias(gene,alias)
                #!#are these reused?#!# for alias in string.split(data['Previous Symbols'],', '): self.addAlias(gene,alias)
                ## ~ [2d] Make dictionary of Genecards data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gdict = {}
                gdict['Symbol'] = gene 
                gdict['Desc'] = data['Approved Name']
                gdict['Entrez'] = data['Entrez Gene ID']
                if not gdict['Entrez']: gdict['Entrez'] = data['Entrez Gene ID (mapped data)']
                gdict['OMIM'] = data['OMIM ID (mapped data)']
                gdict['UniProt'] = data['UniProt ID (mapped data)']
                gdict['EnsEMBL'] = ensgene = data['Ensembl ID (mapped data)']
                gdict['HGNC'] = string.replace(hgnc,'HGNC:','')
                ## ~ [2e] Update self.dict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.addData(gene,gdict)
            self.list['Approved'].sort()
            self.printLog('\r#HGNC','Processed HGNC: %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
            self.printLog('\r#HGNC','%s Approved HGNC symbols' % (rje.integerString(len(self.list['Approved']))))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadMGI(self,sourcefile):   ### Loads MGI mouse data                                                        |1.4|
        '''
        Loads MGI mouse data (download from ftp://ftp.informatics.jax.org/pub/reports/MGI_MouseHumanSequence.rpt)
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Read into dictionary with MGI ID as key ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.errorLog('MGI file "%s" not found' % (sourcefile),printerror=False)
            db = rje_db.Database(self.log,self.cmd_list)
            mdb = db.addTable(sourcefile,mainkeys=['MGI Marker Accession ID'],name='mgi')
            aliaii = {}     # Dictionary of withdrawn symbols to map
            ## ~ [1a] Cleanup table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rename = {'MGI Marker Accession ID':'MGI','Mouse Symbol':'Symbol','Mouse Name':'Desc',
                      'Mouse Entrez Gene ID':'Entrez','Ensembl Gene ID':'EnsEMBL','VEGA Gene ID':'VEGA',
                      'SwissProt IDs':'UniProt','Mouse Synonyms':'Aliases'}
            reject = ['cM Position','NCBI Gene Chromosome','NCBI Gene Start','NCBI Gene End','NCBI Gene Strand',
                      'Ensembl Gene Chromosome','Ensembl Gene Start','Ensembl Gene End','Ensembl Gene Strand',
                      'VEGA Gene Chromosome','VEGA Gene Start','VEGA Gene End','VEGA Gene Strand','Ensembl Transcript IDs',
                      'Ensembl Protein IDs','VEGA Transcript IDs','VEGA Protein IDs','InterPro IDs',
                      'Human Entrez Gene ID','Human Name','Human Chr','Human RefSeq IDs','Human Synonyms']
            extras = ['GenBank IDs','Unigene IDs','RefSeq IDs','Human Symbol']
            mdb.dropFields(reject)
            for field in rename:
                mdb.renameField(field,rename[field])
                if rename[field] not in self.list['Headers']: self.list['Headers'].append(rename[field])
            for field in extras:
                if field not in self.list['Headers']: self.list['Headers'].append(field)
            ### ~ [2] Parse out information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mx,mtot) = (0.0,mdb.entryNum())
            for data in mdb.entries():
                self.progLog('\r#MGI','Processing MGI: %.1f%%' % (mx/mtot)); mx += 100.0
                gene = data['Symbol'].upper()
                if gene not in self.list['Approved']: self.list['Approved'].append(gene)
                ## ~ [2a] Add aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for alias in string.split(data['Aliases'],', '): self.addAlias(gene,alias)
                ## ~ [2b] Make dictionary of Genecards data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gdict = {}
                for field in rename.values() + extras: gdict[field] = data[field]
                gdict['MGI'] = string.replace(data['MGI'],'MGI:','')
                ## ~ [2e] Update self.dict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.addData(gene,gdict)
            self.list['Approved'].sort()
            self.printLog('\r#MGI','Processed MGI: %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
            self.printLog('\r#MGI','%s MGI symbols' % (rje.integerString(len(self.list['Approved']))))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadSource(self,sourcefile):  ### Loads Source data                                                         |1.2|
        '''
        Loads Source data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.log.errorLog('Source file "%s" not found' % (sourcefile),printerror=False)
            data = rje.dataDict(self,sourcefile,getheaders=True)
            for h in data.pop('Headers'):
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            ### ~ [2] Parse out Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot) = (0.0,len(data))
            for gene in rje.sortKeys(data):
                self.progLog('\r#DATA','Processing %s: %.1f%%' % (sourcefile,hx/htot))
                hx += 100.0
                ## ~ [2a] Update self.dict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.addData(gene,data[gene])
            self.printLog('\r#DATA','Processed %s: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))           
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadAlias(self,sourcefile):  ### Loads Alias data
        '''
        Loads Alias data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.log.errorLog('Alias file "%s" not found' % (sourcefile),printerror=False)
            data = rje.dataDict(self,sourcefile,datakeys=['Aliases'],lists=True)
            ### ~ [2] Parse out Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot) = (0.0,len(data))
            for id in data:
                self.log.printLog('\r#ALIAS','Processing %s: %.1f%%' % (sourcefile,hx/htot),newline=False,log=False)
                hx += 100.0
                ## ~ [2a] Update self.dict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for alist in data[id]['Aliases']:
                    for alias in string.split(alist,','): self.addAlias(id,alias)
            self.log.printLog('\r#ALIAS','Processed %s: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))           
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadPickle(self,sourcefile):  ### Loads Pickle data                                                         |1.2|
        '''
        Loads Pickle data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.errorLog('Pickle file "%s" not found' % (sourcefile),printerror=False)
            gunzip = sourcefile[-3:] == '.gz' and not self.opt['Win32']
            if gunzip:
                self.gUnzip(sourcefile)
                sourcefile = sourcefile[:-3]
            try:
                self.printLog('#LOAD','Attempting to load pickle.',log=False)
                mypickle = pickle.load(open(sourcefile,'r'))
                self.printLog('#LOAD','Pickle loaded: %s.' % (sourcefile))
            except: return self.log.errorLog('Problem loading %s' % sourcefile)
            if gunzip: self.gZip(sourcefile)
            ### ~ [2] Process Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for h in mypickle.list['Headers']:
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            (ix,itot) = (0.0,len(mypickle.dict['Data']))
            for k in rje.sortKeys(mypickle.dict['Data']):
                self.printLog('\r#PICKLE','Updating Data from pickle: %.1f%%' % (ix/itot),newline=False,log=False)
                ix += 100.0
                self.addData(k,mypickle.dict['Data'].pop(k))
            self.printLog('\r#PICKLE','Processed %s Data: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))           
            (ix,itot) = (0.0,len(mypickle.dict['TempData']))
            for k in rje.sortKeys(mypickle.dict['TempData']):
                self.printLog('\r#PICKLE','Updating TempData from pickle: %.1f%%' % (ix/itot),newline=False,log=False)
                ix += 100.0
                self.addData(k,mypickle.dict['TempData'].pop(k))
            self.printLog('\r#PICKLE','Processed %s TempData: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))           
            (ix,itot) = (0.0,len(mypickle.dict['Alias']))
            for k in rje.sortKeys(mypickle.dict['Alias']):
                self.printLog('\r#PICKLE','Updating Aliases from pickle: %.1f%%' % (ix/itot),newline=False,log=False)
                ix += 100.0
                for a in mypickle.dict['Alias'].pop(k): self.addAlias(k,a)
            self.printLog('\r#PICKLE','Processed %s Aliases: %s key IDs; %s temp IDs; %s aliases' % (sourcefile,rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
            if mypickle.obj['EnsLoci']: self.obj['EnsLoci'] = mypickle.obj['EnsLoci']
            self.dict['EnsLoci'] = rje.combineDict(self.dict['EnsLoci'],mypickle.dict['EnsLoci'])
            for a in mypickle.list['Approved']:
                if a not in self.list['Approved']: self.list['Approved'].append(a)
            self.list['Approved'].sort()
            self.printLog('\r#HGNC','%s Approved HGNC symbols' % (rje.integerString(len(self.list['Approved']))))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadEnsMap(self,sourcefile):   ### Extracts EnsEMBL mapping data from a BioMart download                               |3.0|
        '''Extracts EnsEMBL mapping data from a BioMart download.'''
        ### ~ [0] ~ Setup paths and files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if sourcefile.lower() in ['','none']: return 
        if not os.path.exists(sourcefile): return self.errorLog('EnsEMBL map file "%s" missing' % sourcefile,printerror=False)
        ### ~ [1] ~ Parse Sequence IDs Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [1a] ~ Read in data and adjust for version of download (headers) ~~~~~~~~~~~~~~~~~~~~~ ##
        if open(sourcefile,'r').readline().find('HGNC Symbol') > 0: mainkeys = ['Ensembl Gene ID','Ensembl Peptide ID','HGNC Symbol']
        elif open(sourcefile,'r').readline().find('MGI symbol') > 0: mainkeys = ['Ensembl Gene ID']
        else: mainkeys = ['Ensembl Gene ID','Ensembl Peptide ID','HGNC symbol']     
        if open(sourcefile,'r').readline().find('Ensembl Protein ID') > 0: mainkeys[1] = 'Ensembl Protein ID'
        ensdata = rje.dataDict(self,sourcefile,mainkeys,'all')
        ## ~ [1b] ~ Extract useful data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        epx = 0
        for map in ensdata:
            self.progLog('\r#ENS','%s EnsEMBL mappings' % rje.integerString(epx)); epx += 1
            gene = ensdata[map]['Ensembl Gene ID']
            if len(mainkeys) > 1:
                try: pept = ensdata[map]['Ensembl Peptide ID']
                except: pept = ensdata[map]['Ensembl Protein ID']
                self.addAlias(gene,pept)
                try: hgnc = ensdata[map]['HGNC Symbol']
                except: hgnc = ensdata[map]['HGNC symbol']
                hgnc = hgnc.upper()     # Make upper case
                self.addAlias(hgnc,gene)
                mgi = ''
            else:
                mgi = ensdata[map]['MGI symbol'].upper()
                mgid = ensdata[map]['MGI ID']
                if mgi: self.addAlias(mgi,gene)
                if mgid: self.addAlias(mgi,mgid)
                hgnc = ''
            if hgnc: self.addData(hgnc,{'EnsEMBL':gene})
            elif mgi: self.addData(mgi,{'EnsEMBL':gene})
            else: self.addData(gene,{'EnsEMBL':gene})
        self.printLog('\r#ENS','%s EnsEMBL mappings; %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(epx),rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))
#########################################################################################################################
    def loadEnsLoci(self,sourcefile):  ### Loads EnsLoci data
        '''
        Loads EnsLoci data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.log.errorLog('EnsLoci file "%s" not found' % (sourcefile),printerror=False)
            ecmd = ['accnr=F','seqnr=F'] + self.cmd_list + ['seqin=%s' % sourcefile,'autoload=T']
            ensloci = self.obj['EnsLoci'] = rje_seq.SeqList(self.log,ecmd)
            for h in ['EnsLoci','EnsEMBL']:
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            ### ~ [2] Parse Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ex,etot) = (0.0,ensloci.seqNum())
            for eseq in ensloci.seqs():
                self.progLog('\r#ENS','Processing EnsLoci: %.1f%%' % (ex/etot)); ex += 100.0
                try:
                    loci = eseq.shortName()
                    edata = rje.matchExp('\[acc:(\S+) pep:(\S+) gene:(\S+)\]',eseq.info['Name'])
                    if not edata:
                        edata = rje.matchExp('\[acc:(\S+) pep:(\S+) symbol:(\S+) gene:(\S+)\]',eseq.info['Name'])
                        self.addAlias(edata[-2],edata[-1])
                    ensg = edata[-1]
                    for alias in [loci,edata[0],edata[1]]: self.addAlias(ensg,alias)
                    if edata[0][:1] in ['P','Q']: self.addData(ensg,{'EnsLoci':loci,'EnsEMBL':ensg,'UniProt':edata[0]})
                    elif 'MOUSE' in sourcefile and edata[0][:1] in '0123456789': self.addData(ensg,{'EnsLoci':loci,'EnsEMBL':ensg,'MGI':edata[0]})
                    else: self.addData(ensg,{'EnsLoci':loci,'EnsEMBL':ensg})
                    if loci not in self.dict['EnsLoci']: self.dict['EnsLoci'][loci] = eseq
                except: self.errorLog('Problem processing %s' % eseq.info['Name'])
            self.printLog('\r#ENS','Processed EnsLoci: %s key IDs; %s temp IDs; %s aliases' % (rje.integerString(len(self.dict['Data'])),rje.integerString(len(self.dict['TempData'])),rje.integerString(len(self.dict['Alias']))))           
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def addAliasFile(self,sourcefile): ### Adds manual aliases to GeneMap object from non-conventional format file
        '''
        Adds manual aliases to GeneMap object from non-conventional format file. Multiple lines of:
        alias = ID1 ID2 ...
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(sourcefile): return
            ### ~ [2] ~ Add Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in self.loadFromFile(sourcefile,chomplines=True):
                ids = string.split(line,' = ')
                if len(ids) < 2: continue
                alias = ids[0]
                for id in string.split(ids[1]): self.addAlias(id,alias)
            self.printLog('#ALIAS','Added aliases >> %s aliases' % rje.integerString(len(self.dict['Alias'])))
        except: self.errorLog('Pingu.addAlias error')
#########################################################################################################################
    def addData(self,id,data,overwrite=False):  ### Adds data to self.dict
        '''
        Adds data to self.dict.
        >> id:str = ID to be used for dictionary key. If missing, will try to extract from data using self.list['XRef']
        >> data:dict = Data read in as dictionary for ID
        >> overwrite:bool [False] = Whether to overwrite existing data.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['KeyID'] in data and data[self.info['KeyID']].upper() in badalias: data.pop(self.info['KeyID'])
            if id and id in self.dict['Alias']:
                klist = self.getKeyID(id)
                #x#self.deBug('%s already mapped to %s' % (id,klist))
            elif self.info['KeyID'] in data: klist = [data[self.info['KeyID']].upper()]
            else: klist = []
            for appkey in ['Symbol','Approved Symbol']:
                if appkey in data and data[appkey] not in self.list['Approved']: self.list['Approved'].append(data[appkey])
            ## ~ [0a] ~ Get ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not id:
                for x in self.list['XRef']:
                    if x in data and data[x]:
                        id = data[x]
                        break
            #!# Add looking for id in known aliases mapping to KeyID and getting k #!#
            if not klist:
                if id in self.dict['Data']: klist = [id]
                elif id not in self.dict['TempData']: klist = self.getKeyID(id)     # Speed things up at his stage
            
            ### ~ [1] ~ Look for KeyID and update Data dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for k in klist:
                ## ~ [1a] ~ Add aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                
                if id and k != id.upper(): self.addAlias(k,id)
                try: 
                    if data['Entrez']: self.addAlias(k,'LOC%s' % (data['Entrez']))
                except: pass
                for x in self.list['XRef']:
                    if x in data and data[x]:
                        if x in ['Entrez','OMIM','HGNC','HPRD','MGI']: self.addAlias(k,'%s%s' % (x,data[x]))
                        elif ',' in data[x]:
                            for alias in string.split(data[x],','): self.addAlias(id,alias)
                        else: self.addAlias(k,data[x])
                ## ~ [1b] ~ Update data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if k not in self.list['KeyID']: self.list['KeyID'].append(k)
                if k not in self.dict['Data']:
                    if k in self.dict['TempData']: self.dict['Data'][k] = self.dict['TempData'].pop(k)
                    else: self.dict['Data'][k] = {}
                self.dict['Data'][k] = rje.combineDict(self.dict['Data'][k],data,overwrite=overwrite)
            if klist: return

            ### ~ [2] ~ Update Temporary Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not id: return self.log.errorLog('Failure to get gene identifier from data dictionary: %s' % (data),printerror=False)
            id = id.upper()
            ## ~ [2a] ~ Add aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                if data['Entrez']: self.addAlias('LOC%s' % (data['Entrez']),id)
            except: pass
            for x in self.list['XRef']:
                if x in data and data[x]:
                    if x in ['Entrez','OMIM','HGNC','HPRD','MGI']: self.addAlias(id,'%s%s' % (x,data[x]))
                    elif ',' in data[x]:
                        for alias in string.split(data[x],','): self.addAlias(id,alias)
                    else: self.addAlias(id,data[x])
            ## ~ [2b] ~ Update data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if id not in self.dict['TempData']: self.dict['TempData'][id] = {}
            self.dict['TempData'][id] = rje.combineDict(self.dict['TempData'][id],data,overwrite=overwrite)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def parseCard(self,alias):  ### Parses relevant GeneCard given Alias                                            |1.2|
        '''
        Parses relevant GeneCard given Alias.
        >> alias:str = Gene symbol or alias
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            url = 'http://www.genecards.org/cgi-bin/carddisp.pl'
            params = 'gene=%s&alias=yes' % alias
            data = {}
            for skipper in self.list['SkipList']:
                if rje.matchExp('(%s)' % string.replace(skipper,'*','\S+'),alias): return False
            gc_regexp = {'HGNC':'href="http:\/\/www.gene.ucl.ac.uk\/cgi-bin\/nomenclature\/get_data.pl\?hgnc_id=(\d+)"',
                         'Entrez':'href="http:\/\/www.ncbi.nlm.nih.gov\/entrez\/query.fcgi.+list_uids=(\S+)"',
                         'UniProt':'href="http:\/\/www\.expasy\.org\/uniprot\/(\S+)"',
                         'EnsEMBL':'href="http\:\/\/www.ensembl.org\/Homo_sapiens\/geneview\?gene=(\S+)"'}

            ### ~ [2] ~ Download GeneCard ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: flines = urllib2.urlopen(url, params).readlines()
            except KeyboardInterrupt: raise
            except: return False

            ### ~ [3] ~ Parse Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for html in flines:
                ## ~ [3a] ~ Primary Gene Symbol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if rje.matchExp('<title>GeneCard for (\S+)<\/title>',html): data['Symbol'] = rje.matchExp('<title>GeneCard for (\S+)<\/title>',html)[0]
                if rje.matchExp('<title>(\S+) GeneCard<\/title>',html): data['Symbol'] = rje.matchExp('<title>(\S+) GeneCard<\/title>',html)[0]
                ## ~ [3b] ~ Additional data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for gtype in gc_regexp:
                    if gtype in data: continue  # Just take first entry
                    if rje.matchExp(gc_regexp[gtype],html): data[gtype] = rje.matchExp(gc_regexp[gtype],html)[0]
                ## ~ [3c] ~ Stop if all read ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if len(data) > len(gc_regexp) or html.find('outside databases for aliases') > 0: break    # Extracted all possible data
            
            ### ~ [4] ~ Finish and update data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Symbol' in data and data['Symbol'] not in self.list['Approved']: self.list['Approved'].append(data['Symbol'])
            if data: 
                self.addData(alias,data)
                return True
            else: return False
        except KeyboardInterrupt:
            self.log.errorLog('Parsing (%s) cancelled.' % alias)
            raise
        except:
            self.log.errorLog('Error in parseCard(%s)' % alias)
            return False
#########################################################################################################################
    def loadPFamData(self): ### Loads PFamData from a file in self.dict['PFam']
        '''Loads PFamData from a file in self.dict['PFam'].'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = rje.dataDict(self,self.info['PFamData'],datakeys=['Name'],lists=True)
            ### ~ [2] ~ Add to dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pfam in rje.sortKeys(data):
                for ens in data[pfam]['Name']:
                    if ens not in self.dict['PFam']: self.dict['PFam'][ens] = []
                    if pfam not in self.dict['PFam'][ens]: self.dict['PFam'][ens].append(pfam)
        except: self.errorLog('GeneMap.loadPFamData error.')            
#########################################################################################################################
    ### <5> ### Data output Methods                                                                                     #
#########################################################################################################################
    def flatOut(self):  ### Saves flat files
        '''Saves flat files.'''
        try: ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            delimit = self.info['Delimit']
            dfile = '%s.data.%s' % (self.info['Basefile'],rje.delimitExt(delimit))
            afile = '%s.alias.%s' % (self.info['Basefile'],rje.delimitExt(delimit))

            ### ~ [1] ~ Save Main Data File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('#OUT','Saving main data to %s ...' % dfile,newline=False,log=False)
            headers = ['Gene'] + self.list['Headers']
            rje.delimitedFileOutput(self,dfile,headers,delimit,rje_backup=True)
            for gene in rje.sortKeys(self.dict['Data']):
                data = rje.combineDict({'Gene':gene},self.dict['Data'][gene],overwrite=False)
                rje.delimitedFileOutput(self,dfile,headers,delimit,data)
            self.log.printLog('\r#OUT','Saved main data for %s genes to %s.' % (rje.integerString(len(self.dict['Data'])),dfile))

            ### ~ [2] ~ Save list of Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            rev = self.reverseAliases()
            self.log.printLog('\r#REV','Reversed: %s IDs with aliases' % rje.integerString(len(rev)))
            self.log.printLog('#OUT','Saving aliases to %s ...' % afile,newline=False,log=False)
            headers = ['Gene','Aliases']
            rje.delimitedFileOutput(self,afile,headers,delimit,rje_backup=True)
            for gene in rje.sortKeys(rev):
                rev[gene].sort()
                data = {'Gene':gene,'Aliases':string.join(rev[gene],',')}
                rje.delimitedFileOutput(self,afile,headers,delimit,data)
            self.log.printLog('\r#OUT','Saved aliases for %s genes to %s.' % (rje.integerString(len(rev)),afile))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def pickleOut(self):  ### Saves self as pickle (and gzips)
        '''Saves self as pickle (and gzips).'''
        pfile = '%s.pickle' % self.info['Basefile']
        try:
            self.log.printLog('#SAVE','Attempting to save pickle.',log=False)
            pickle.dump(self,open(pfile,'w'))
            self.log.printLog('#SAVE','Python Pickle saved as %s.' % (pfile))
        except: self.log.errorLog('Problem pickling to %s' % pfile)
        try:
            if os.path.exists(pfile) and not self.opt['Win32']:
                self.log.printLog('#GZIP','Zipping %s ...' % (pfile),newline=False,log=False)
                if os.path.exists('%s.gz' % pfile): os.unlink('%s.gz' % pfile)
                os.system('gzip %s' % pfile)
                self.log.printLog('\r#GZIP','%s zipped -> %s.gz' % (pfile,pfile))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <6> ### Data output Methods                                                                                     #
#########################################################################################################################
    def reduceGeneData(self,genelist,data=True,aliases=True):   ### Reduces data to genes given in genelist
        '''
        Reduces data to genes given in genelist.
        >> genelist:list = gene identifiers to (map and) retain.
        >> data:bool [True] = whether to reduce self.dict['Data']
        >> aliases:bool [True] = whether to reduce self.dict['Alias']
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            keepgenes = []
            gx = 0.0; gtot = len(genelist)
            for gene in genelist:
                self.progLog('\r#MAP','Mapping genes for reduceGeneData: %.2f%%' % (gx/gtot)); gx += 100.0
                keepgenes.append(self.bestMap(gene))
            self.printLog('\r#MAP','Mapped %d of %d genes for reduceGeneData.' % (len(keepgenes),gtot)); gx += 100.0
            #keepgenes.sort(); self.deBug(keepgenes[-20:])
            #self.deBug(rje.sortKeys(self.dict['Data'])[-20:])
            ### ~ [1] Gene Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = len(self.dict['Data'])
            for gene in rje.sortKeys(self.dict['Data']):
                self.progLog('\r#DATA','Reducing mapped genes: %.2f%%' % (gx/gtot)); gx += 100.0
                if gene not in keepgenes: self.dict['Data'].pop(gene)
            self.printLog('\r#DATA','Reduced GeneMap data to %d genes' % (len(self.dict['Data'])))
            ### ~ [2] Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = len(self.dict['Alias'])
            for alias in rje.sortKeys(self.dict['Alias']):
                self.progLog('\r#ALIAS','Reducing mapped aliases: %.2f%%' % (gx/gtot)); gx += 100.0
                keepme = False
                for gene in self.mapGene(alias):
                    if gene in keepgenes: keepme = True; break
                if not keepme: self.dict['Alias'].pop(alias)
            self.printLog('\r#ALIAS','Reduced Alias data to %d aliases' % (len(self.dict['Alias'])))
        except: self.errorLog('GeneMap.reduceGeneData error')
#########################################################################################################################
### End of SECTION II: GeneMap Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def loadGeneMapPickle(sourcefile): return pickle.load(open(sourcefile,'r'))
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
    try:
        genemap = GeneMap(mainlog,cmd_list)
        genemap.run()
        if genemap.opt['Test']:
            for id in ['CDC42','CDC42P2','ENSG00000070831']:
                print id, ':getEnsLoci:', genemap.getEnsLoci(id)
                print id, ':getKeyID:', genemap.getKeyID(id)
                print id, ':getGeneData:', genemap.getGeneData(id)
                print id, ':mapGene:', genemap.mapGene(id)
                print id, ':aliases:', genemap.aliases(id)
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
