#!/usr/bin/python

# rje_uniprot.py - RJE Module to Handle Uniprot Files
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_uniprot
Description:  RJE Module to Handle Uniprot Files
Version:      3.25.3
Last Edit:    21/04/20
Copyright (C) 2007 Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains methods for handling UniProt files, primarily in other rje modules but also with some
    standalone functionality. To get the most out of the module with big UniProt files (such as the downloads from EBI),
    first index the UniProt data using the rje_dbase module.

    This module can be used to extract a list of UniProt entries from a larger database and/or to produce summary tables
    from UniProt flat files. Version 3.14 introduced direct querying from the UniProt website if unipath=None or
    unipath=URL.

    In addition to method associated with the classes of this module, there are a number of methods that are called from
    the rje_dbase module (primarily) to download and process the UniProt sequence database.

    Version 3.19 has seen an over-haul of the dbxref extraction. `dblist=LIST` and `dbparse=LIST` can now be used
    largely synonymously. Rather than extract all db xref by default, there is now a default list of databases to parse:
    ['UniProtKB/Swiss-Prot','ensembl','REFSEQ','HGNC','Entrez Gene','FlyBase','Pfam','GO','MGI','ZFIN'].

Input/Output:
    ### ~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    unipath=PATH    : Path to UniProt Datafile (looks here for DB Index file; url = use Web downloads) [url]
    dbindex=FILE    : Database index file [uniprot.index]
    uniprot=FILE    : Name of UniProt file [None]
    extract=LIST    : Extract IDs/AccNums in list. LIST can be FILE or list of IDs/AccNums X,Y,.. []
    acclist=LIST    : As extract=LIST.
    uniprotid=LIST  : As extract=LIST.
    acctable=FILE   : Delimited text file from which to retrieve a list of accession numbers [None]
    accfield=X      : Accession number field for acctable=FILE extraction [UniProt]
    specdat=LIST    : Make a UniProt DAT file of the listed species from the index (over-rules extract=LIST) []
    proteome=LIST   : Extract complete proteomes for listed Taxa ID (e.g. 9606 for human) []
    taxonomy=LIST   : Extract all entries for listed Taxa ID (e.g. 4751 for fungi) []
    usebeta=T/F     : Whether to use beta.uniprot.org rather than www.uniprot.org [False]
    splicevar=T/F   : Whether to search for AccNum allowing for splice variants (AccNum-X) [True]
    tmconvert=T/F   : Whether to convert TOPO_DOM features, using first description word as Type [False]
    reviewed=T/F    : Whether to restrict input to reviewed entries only [False]
    complete=T/F    : Whether to restrict proteome downloads to "complete proteome" sets [False]
    uniformat=X     : Desired UniProt format for proteome download. Append gz to compress. [txt]
                        - html | tab | xls | fasta | gff | txt | xml | rdf | list | rss
    onebyone=T/F    : Whether to download one entry at a time. Slower but should maintain order [False]

    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : If set, can use "T" or "True" for other `*out` options. (Will default to `datout` if given) [None]
    datout=T/F/FILE : Name of new (reduced) UniProt DAT file of extracted sequences [None]
    tabout=T/F/FILE : Table of extracted UniProt details [None]
    ftout=T/F/FILE  : Outputs table of features into FILE [None]
    pfamout=T/F/FILE: Whether to output a "long" table of (accnum, pfam, name, num) [False]
    domtable=T/F    : Makes a table of domains from uniprot file [False]
    gotable=T/F     : Makes a table of AccNum:GO mapping [False]
    cc2ft=T/F       : Extra whole-length features added for TISSUE and LOCATION (not in datout) [False]
    xrefout=T/F/FILE: Table of extracted Database xref (Formerly linkout=FILE) [None]
    longlink=T/F    : Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]
    dblist=LIST     : List of databases to extract (extract all if empty or contains 'all') [see above]
    dbsplit=T/F     : Whether to generate a table per dblist database (basefile.dbase.tdt) [False]
    dbdetails=T/F   : Whether to extract full details of DR line rather than parsing DB xref only [False]
    append=T/F      : Append to results files rather than overwrite [False]

    ### ~ UniProt Conversion Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ucft=X          : Feature to add for UpperCase portions of sequence []
    lcft=X          : Feature to add for LowerCase portions of sequence []
    maskft=LIST     : List of Features to mask out []
    invmask=T/F     : Whether to invert the masking and only retain maskft features [False]
    caseft=LIST     : List of Features to make upper case with rest of sequence lower case []

Specialist Options:
    ### ~ Parsing Options (Programming Only) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    fullref=T/F     : Whether to store full Reference information in UniProt Entry objects [False]
    dbparse=LIST    : List of databases to parse from DR lines in UniProtEntry object [see code]
    uparse=LIST     : Restricted lines to parse to accelerate parsing of large datasets for limited information []

    ### ~ General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    memsaver=T/F    : Memsaver option to save memory usage - does not retain entries in UniProt object [False]
    cleardata=T/F   : Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]
    specsleep=X     : Sleep for X seconds between species downloads [60]

    ### ~ UniProt Download Processing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    makeindex=T/F   : Generate UniProt index files [False]
    makespec=T/F    : Generate species table [False]
    makefas=T/F     : Generate fasta files [False]
    grepdat=T/F     : Whether to use GREP in attempt to speed up processing [False]

Uses general modules: glob, os, re, string, sys, time
Uses RJE modules: rje, rje_db, rje_sequence
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, re, string, sys, time
try:
    import urllib2 as urllib2
except:
    import urllib.request as urllib2
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_db, rje_sequence
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial working version for interaction_motifs.py
    # 1.1 - Minor tidying and modification
    # 2.0 - Moved functions to rje_dbase. Added option to extract using index files.
    # 2.1 - Added possibility to extract splice variants
    # 2.2 - Added output of feature table for the entries in memory (not compatible with memsaver mode)
    # 2.3 - Added ID to tabout and also added accShortName() method to extract dictionary of {acc:ID__PrimaryAcc}
    # 2.4 - Add method for converting Sequence object and dictionary into UniProt objects... and saving
    # 2.5 - Added cc2ft Extra whole-length features added for TISSUE and LOCATION [False] and ftout=FILE
    # 2.6 - Added features based on case of sequence. (Uses seq.dict['Case'])
    # 2.7 - Added masking of features - Entry.maskFT(type='EM',inverse=False)
    # 2.8 - Added making of Taxa-specific databases using a list of UniProt Species codes
    # 2.9 - Added extraction of EnsEMBL, HGNC, UniProt and EntrezGene from IPI DAT file.
    # 3.0 - Added some module-level methods for use with rje_dbase.
    # 3.1 - Added extra linking of databases from UniProt entries
    # 3.2 - Added feature masking and TM conversion.
    # 3.3 - Added DBase processing options.
    # 3.4 - Made modifications to allow extended EMBL functionality as part of rje_embl.
    # 3.5 - Added SplitOut to go with rje_embl V0.1
    # 3.6 - Added longlink=T/F  : Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]
    # 3.7 - Added cleardata=T/F : Whether to clear unprocessed Entry data or retain in Entry & Sequence objects [True]
    # 3.8 - Added extraction of NCBI Taxa ID.
    # 3.9 - Added grepdat=T/F     : Whether to use GREP in attempt to speed up processing [False]
    # 3.10- Added forking for speeding up of processing.
    # 3.11- Added storing of Reference information in UniProt entries.
    # 3.12- Added addition accdict extraction method for all entries read in.
    # 3.13- Minor bug fix for link table output.
    # 3.14- Added direct retrieval of UniProt entries from URL, including full proteomes. Updated output file naming.
    # 3.14- Added dblist=LIST and dbsplit=T/F for additional DB link output control. Set unipath default to url.
    # 3.15- Added extraction of taxonomic groups. Add UniFormat to improve pure downloads.
    # 3.16- Added WBGene ID's from WormBase as one of the recognised DB XRef to parse.
    # 3.17- Efficiency tweak to URL-based extraction of acclist.
    # 3.18- Minor modification to database parsing.
    # 3.19- Updated and consolidated dbxref table generation (formerly linkout) using rje_db. Changed acc_num to accnum.
    #     - Added gotable=T generation of GO table. Fixed makeindex to use a single fork if needed.
    # 3.20- Updated dbsplit=T output and checked function with Pfam. Probably needs work for other databases.
    # 3.20.1 - Added uniprotid=LIST as an alias to acclist=LIST and extract=LIST.
    # 3.20.2 - Added extra sequence return methods to UniprotEntry. Added fasta REST output.
    # 3.20.3 - Fixed bug if new uniprot extraction method fails.
    # 3.20.4 - Fixed bug introduced by REST access modifications.
    # 3.20.5 - Improved handling of downloads for uniprot IDs that have been updated (i.e. no direct mapping).
    # 3.20.6 - Improved handling of zero accession numbers for extraction.
    # 3.20.7 - Fixed uniformat default error.
    # 3.21.0 - Added uparse=LIST option to try and accelerate parsing of large datasets for limited information.
    # 3.21.1 - FullText is no longer stored in Uniprot object. Will need special handling if required.
    # 3.21.2 - Fixed single uniprot extraction bug.
    # 3.21.3 - Added REST datout to proteomes extraction.
    # 3.21.4 - Fixed Feature masking. Should this be switched off by default?
    # 3.22.0 - Tweaked REST table output.
    # 3.23.0 - Added accnum map table output. Fixed REST output bug when bad IDs given. Added version and about output.
    # 3.24.0 - Added pfam out and changed map table headers.
    # 3.24.1 - Fixed process Uniprot error when uniprot=FILE given.
    # 3.24.2 - Updated HTTP to HTTPS. Having some download issues with server failures.
    # 3.25.0 - Fixed new Uniprot batch query URL. Added onebyone=T/F    : Whether to download one entry at a time. Slower but should maintain order [False].
    # 3.25.1 - Fixed proteome download bug following Uniprot changes.
    # 3.25.2 - Fixed Uniprot protein extraction issues by using curl. (May not be a robust fix!)
    # 3.25.3 - Fixed some problems with new Uniprot feature format.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Lots of functionality to add! Look also to BioPython.
    # [Y] : Modify the searching for entry in acclist to cope with partial matches (exclude them)
    # [ ] : Improve reading/generation of basefile and modify DomTable to work with Memsaver
    # [ ] : Modify additional outputs to use basefile if set to T/True
    # [ ] : Convert DB Links etc. into rje_db tables.
    # [Y] : Add specific database detail extraction, first for IPI DAT files and later for UniProt
    # [Y] : Add a database mapping method for extracting DB cross-refs.
    # [ ] : Change the non-forking processing method to match forked one, which is generally faster!
    # [Y] : Add extraction of complete/reviewed proteome from TaxaID.
    # [ ] : Add extraction of complete/reviewed proteome from SpecCode.
    # [ ] : Update to new module type, use rje_db (have mode=X) and document!
    # [ ] : Rename dblist to dblinks - need to change existing self.list['DBLinks'] first!
    # [Y] : Add dbparse=LIST to parse db into entry.dict['DB'] = {db:list} instead of parsing for table output..
    # [ ] : Add options to compress UniProt URL downloads. Improve downloads of acclist etc.
    # [ ] : Add extraction of GO categories from UniProt entries.
    # [ ] : Properly fix DB Xref parsing.
    # [ ] : Add option to use PICR REST service to fetch (and then parse XLM for) Uniprot IDs.
    # [Y] : When safe, replace fttable and domtable acc_num with accnum
    # [ ] : Add ftlist=LIST subset of features to extract.
    # [ ] : Code up proper zero fork index generation.
    # [ ] : Need to add retrieval in batches to extractProteinsFromURL method.
    # [ ] : Why does extract from a Uniprot file no longer work?!
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('RJE_UNIPROT', '3.25.3', 'April 2020', '2007')
    description = 'RJE Uniprot Parsing/Extraction Module'
    author = 'Dr Richard J. Edwards.'
    comments = []
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit: sys.exit(0)
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        if len(sys.argv) == 2 and sys.argv[1] in ['about','-about','--about']: sys.exit(0)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit(0)
    except KeyboardInterrupt: sys.exit()
    except:
        print('Problem during initial setup.')
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: MODULE CONSTANTS                                                                                        #
#########################################################################################################################
### UniProt Parsing dictionary: Add crucial information to parse out here. Used by UniProtEntry.process()   ###
uniparse = {
    'ID' : '\s+'.join(['^(\S+)','(\S.+);','(\d+)\s+(AA|BP)']),    # ID, Type, Length
    'AC' : '(\S+);',     # Primary Acc
    'DE' : '\s*(\S.+)',  # Description
    'GN' : 'Name=(\S+);?',   # Gene Name
    'SY' : 'Synonyms=(\S.+);',   # Gene Synonyms
    'OS' : '^(\S.+)\s*$',   # Species
    'OX' : '^OX\s+NCBI_TaxID=(\d+)',    # NCBI Taxa ID
    'RN' : '\[(\d+)\]',     # Reference number
    'RX' : 'PubMed=(\d+);', # PubMed ID
    'RC' : 'TISSUE=(.+);',  # Tissue(s)
    'DR' : '^(\S+);\s+(\S.+)$',  # Database links (Dbase,Details)
    'CC' : '^-!-\s+(\S.+):\s+(\S.+)$',  # Comments (Type, Details)
    'FT' : '\s+'.join(['(\S+)','<*(\d+)','>*(\d+)\.*','(\S.+)\s*$'])   # Feature type, start, stop, desc
    }
#########################################################################################################################
useful_data = ['ID','AC','DE','GN','OS','OC','OX','RX','CC','DR','RC','KW','FT']     # Data to retain following parsing # ?? #
#!# NB. This list is not currently used! #!#
#########################################################################################################################
featurelist = ['LIPID','TRANSMEM','MOD_RES','DOMAIN']   #!# Features for function table. Add more! #!#
#########################################################################################################################
# Set defaults for specialDB() parsing during entry.process()                                                      #V3.14
dbparse = ['UniProtKB/Swiss-Prot','ensembl','REFSEQ','HGNC','Entrez Gene','FlyBase','Pfam','GO','MGI','ZFIN','gene','id']
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: UniProt Class                                                                                          # 
#########################################################################################################################
class UniProt(rje.RJE_Object):     
    '''
    UniProt Download Class. Author: Rich Edwards (2005).

    Info:str
    - AccTable = Delimited text file from which to retrieve a list of accession numbers [None]
    - AccField = Accession number field for acctable=FILE extraction [UniProt]
    - Name = Name of UniProt File 
    - UniPath = Path to UniProt Datafile (will look here for DB Index file) [UniProt/]
    - DBIndex = Database index file [uniprot.index]
    - DATOut = Name of new (reduced) UniProt DAT file of extracted sequences [None]
    - PfamOut = Whether to output a "long" table of (accnum, pfam, desc) [False]
    - TabOut = Name of table of extracted sequence details [None]
    - XRefOut = Table of extracted Database links [None]
    - FTOut = Outputs table of features into FILE [None]
    - SplitOut = If path given, will split output into individual files per entry into PATH []
    - UCFT = Feature to add for UpperCase portions of sequence []
    - LCFT = Feature to add for LowerCase portions of sequence []
    - UniFormat = URL format to download (see UniProt REST services) [txt]
    
    Opt:boolean
    - ClearData = Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]
    - Complete = Whether to restrict proteome downloads to "complete proteome" sets [False]
    - DBDetails = Whether to extract full details of DR line rather than parsing DB xref only [False]
    - DBSplit = Whether to generate a table per database (basefile.dbase.tdt) [False]
    - DomTable = Makes a table of domains from uniprot file [False]
    - FullRef = Whether to store full Reference information in UniProt Entry objects [False]
    - GOTable = Makes a table of AccNum:GO mapping [False]
    - GrepDat = Whether to use GREP in attempt to speed up processing [False]
    - LongLink = Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]    
    - MakeIndex = Generate UniProt index files [False]
    - MakeSpec = Generate species table [False]
    - MakeFas = Generate fasta files [False]
    - OneByOne=T/F    : Whether to download one entry at a time. Slower but should maintain order [False]
    - PfamOut = Whether to output a "long" table of (accnum, pfam, name, num) [None]
    - Reviewed = Whether to restrict input to reviewed entries only [False]
    - SpliceVar = Whether to search for AccNum allowing for splice variants (AccNum-X) [True]
    - UseBeta = Whether to use beta.uniprot.org rather than www.uniprot.org [False]

    Stat:numeric
    - SpecSleep = Sleep for X seconds between species downloads [60]

    List:list
    - DBList = List of databases to extract (extract all if empty) []
    - Entry = list of UniProt Entries
    - Extract = Extract AccNums/IDs in list. LIST can be FILE or list of AccNums X,Y,.. []
    - Proteome = Extract complete proteomes for listed Taxa ID (e.g. 9606 for human) []
    - SpecDat = Make a UniProt DAT file of the listed species from the index []
    - UParse = Restricted lines to parse to accelerate parsing of large datasets for limited information []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
    ### Attributes
    def entryNum(self): return len(self.list['Entry'])
    def entries(self): return self.list['Entry']
#####################4####################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AccTable','AccField','Name','UniPath','DBIndex','DATOut','TabOut','XRefOut','FTOut','PfamOut',
                         'UCFT','LCFT','SplitOut','UniFormat']
        self.optlist = ['Complete','DBSplit','DBDetails','DomTable','LongLink','MakeIndex','MakeSpec','MakeFas',
                        'SpliceVar','ClearData','GrepDat','FullRef','Reviewed','UseBeta','GOTable',
                        'DATOut','TabOut','XRefOut','FTOut','MapOut','PfamOut','OneByOne']
        self.statlist = ['SpecSleep']
        self.listlist = ['DBList','Extract','Entry','Proteome','SpecDat','Taxonomy','UParse']
        self.dictlist = []
        self.objlist = ['DB']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'UniPath':rje.makePath('url'),'DBIndex':'uniprot.index','UCFT':'','LCFT':'','AccField':'UniProt',
                      'UniFormat':'txt'})
        self.setOpt({'SpliceVar':True,'LongLink':True,'ClearData':True,'DBDetails':False,'OneByOne':False})
        self.setInt({'SpecSleep':60})
        self.list['DBList'] = dbparse[0:]
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] ~ Read in commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:### General Options ###
                self._generalCmd(cmd)
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options ###
                self._cmdReadList(cmd,'path',['UniPath'])
                self._cmdReadList(cmd,'file',['AccTable','DBIndex','DATOut','TabOut','XRefOut','FTOut','PfamOut'])
                self._cmdRead(cmd,'file','XRefOut','linkout')
                self._cmdReadList(cmd,'info',['AccField','UCFT','LCFT'])
                self._cmdReadList(cmd.lower(),'info',['UniFormat'])
                self._cmdReadList(cmd,'opt',['Complete','DomTable','LongLink','MakeIndex','SpliceVar','MakeSpec','OneByOne',
                                             'MakeFas','ClearData','GrepDat','FullRef','Reviewed','DBSplit','DBDetails',
                                             'UseBeta','GOTable','DATOut','TabOut','XRefOut','FTOut','MapOut','PfamOut'])
                self._cmdReadList(cmd,'int',['SpecSleep'])
                self._cmdReadList(cmd,'list',['DBList','Extract','Proteome','SpecDat','Taxonomy','UParse'])
                self._cmdRead(cmd,type='list',att='DBList',arg='dbparse')
                self._cmdRead(cmd,type='file',att='Name',arg='uniprot')
                self._cmdRead(cmd,type='list',att='Extract',arg='acclist')
                self._cmdRead(cmd,type='list',att='Extract',arg='uniprotid')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [2] ~ Special processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('GOTable'): self.list['DBList'].append('GO')
        if 'all' in string.join(self.list['DBList']).lower(): self.list['DBList'] = []
        else: self.list['DBList'] = rje.listLower(self.list['DBList'])
        delimit = rje.getDelimit(self.cmd_list)
        for ofile in ['DBIndex','DATOut','TabOut','XRefOut','FTOut','MapOut','PfamOut']:
            if self.getStr(ofile).lower() in ['none','f','false']: self.setStr({ofile:''})
            # Where a boolean has been added too, a straight call of the parameter will switch on output
            if ofile in self.opt and self.getOpt(ofile) and not self.getStr(ofile): self.setStr({ofile:'t'})
        if self.basefile().lower() in ['','none'] and self.getStr('DATOut').lower() not in ['','t','true']: self.basefile(rje.baseFile(self.getStr('DATOut')))
        if self.basefile().lower() in ['','none'] and self.getStr('Name'): self.basefile(rje.baseFile(self.getStr('Name')))
        for ofile in ['DATOut','TabOut','XRefOut','FTOut','MapOut','PfamOut']:
            newfile = '%s.%s' % (self.basefile(),ofile[:-3].lower())
            if not newfile.endswith('dat'): newfile = '%s.%s' % (newfile,rje.delimitExt(delimit))
            if self.getStr(ofile).lower() in ['t','true'] or ofile == 'MapOut':
                self.setStr({ofile:newfile})
                self.printLog('#%s' % ofile.upper()[:-3],'New %s filename set: %s' % (ofile,newfile))
#########################################################################################################################
    def run(self):  ### Main Run Method if called direct from commandline
        '''Main Run Method if called direct from commandline. Returns True if no Errors, else False.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Rest'): self.restSetup()
            ### ~ [1] ~ Make Index file or general processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: self.opt['NoForks'] = self.opt['Win32'] or self.opt['NoForks'] or self.stat['Forks'] < 1
            except: self.opt['NoForks'] = True
            if self.opt['MakeIndex'] or self.opt['MakeSpec'] or self.opt['MakeFas']:
                processUniProt(self,makeindex=self.opt['MakeIndex'],makespec=self.opt['MakeSpec'],makefas=self.opt['MakeFas'])
                newuni = UniProt(self.log,self.cmd_list)
                newuni.setBool({'MakeIndex':False,'MakeSpec':False,'MakeFas':False})
                return newuni.run()
            ### ~ [2] ~ Setup for reading and processing Uniprot entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Entry/Proteome extraction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Get Extract List from species #
            if self.list['SpecDat']: self.extractSpecies()
            # Get Extract List from Table #
            elif rje.exists(self.getStr('AccTable')):
                self.list['Extract'] = rje.sortKeys(rje.dataDict(self,self.getStr('AccTable'),[self.getStr('AccField')],[self.getStr('AccField')],lists=True,ignore=['#']))
            ## ~ [2b] ~ Check for input needs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('Name') and not (self.list['Extract'] + self.list['Proteome'] + self.list['Taxonomy']):   # No AccNum!
                self.errorLog('No input file, acclist/proteome/taxonomy list given. Use "help" option for parameters.',printerror=False)
                return False
            ## ~ [2c] ~ Setup Output DAT File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,self.getStr('DATOut'))  # Will skip if empty value
            ## ~ [2d] ~ Setup special Tabular output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('TabOut') and self.getBool('MemSaver'):
                memtext = 'TabOut will not function with MemSaver mode.'
                if self.i() >= 0 and rje.yesNo('%s Switch Memsaver off?' % memtext):
                    self.setBool({'MemSaver':False})
                    self.cmd_list += ['memsaver=F']
                else: self.printLog('#MEM',memtext); self.setStr({'TabOut':''})
            rje.backup(self,self.getStr('TabOut'))  # Will skip if empty value
            ## ~ [2d] ~ Setup general tabular output file db tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupDB(xref=self.getStr('XRefOut') or self.getBool('DBSplit') or self.getStr('PfamOut'),ft=self.getStr('FTOut') or self.getBool('DomTable'))
            #i# DomTable is just a cut-down FTOut with DOMAIN features only! Store in ft table and then make at end.

            ### >>> Delete when safe <<< ###
#            outfiles = ['DATOut','TabOut','XRefOut','FTOut']
#            if self.getBool('DBSplit'):
#                if not self.list['DBList']:
###                    self.warnLog('Cannot use dbsplit=T without dblist=LIST database listing')
#                    self.setBool({'DBSplit':False})
#                else:
#                    if not self.basefile(): self.warnLog('Need to set basefile for proper function!')
#                    outfiles.remove('XRefOut')
#            for file in outfiles:
#                # if self.info[file].lower() == 'none': self.info[file] = ''    # Now handled in _cmdList()
#                rje.backup(self,self.getStr(file))  # Will skip if empty value

            ### ~ [3] ~ Read UniProt File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.readUniProt()

            ### ~ [4] ~ Tabular output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.tableOutput()
            return True
        except:
            self.log.errorLog('Fundamental error during UniProt.run()')
            return False
#########################################################################################################################
    def setupDB(self,xref=True,ft=True):    ### Set up empty database tables to be populated following _readSingleEntry().
        '''
        Set up empty database tables to be populated following _readSingleEntry().
        >> xref:bool [True] = setup xref table.
        >> ft:bool [True] = setup ft table.
        '''
        self.db().addEmptyTable('map',['uniprotid','accnum'],['uniprotid'])  # Mapping of input accnum to output accnum
        if xref: self.db().addEmptyTable('xref',['accnum','db','xref'],['accnum','db'])
        if ft: self.db().addEmptyTable('ft',['accnum','feature','ft_start','ft_end','description'],['accnum','feature','ft_start','ft_end','description'])
#########################################################################################################################
    def restOutput(self,outfmt=None,maxfilesize=0):   ### Controls what is returned by REST server output
        '''Controls what is returned by REST server output.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not outfmt: outfmt = self.getStrLC('Rest')
            if outfmt == 'default': outfmt = self.restOutputDefault()
            if not outfmt: return self  # If no REST output (rest=None), the object itself is returned.
            #!# This catch should be obselete.
            for restfmt in ['dom','go']:
                if restfmt not in self.dict['Output'] or not self.dict['Output'][restfmt] or not rje.exists(self.dict['Output'][restfmt]):
                    self.dict['Output'][restfmt] = 'No output generated.'
            for restfmt in ['ft','tab','xref','map']:
                if restfmt not in self.dict['Output'] or not self.dict['Output'][restfmt] or not rje.exists(self.getStr(self.dict['Output'][restfmt])):
                    self.dict['Output'][restfmt] = 'No output generated.'
            ### ~ [1] RESTful output to return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfmt == 'log': return open(self.log.info['LogFile'],'r').read()
            if outfmt == 'loghtml': return self.restLogHTMLOutput()
            if outfmt in ['full','text']: return self.restFullOutput()
            if outfmt in ['fullhtml','html']: return self.restHTMLOutput()
            if outfmt == 'version': return rje.Out(self.log,['v=-1']).printIntro(self.log.obj['Info'])
            if outfmt == 'outfmt': return self.restSetup.__doc__
            if outfmt == 'warnings': return string.join(self.log.list['WarnLog'],'\n')  # List of log warning messages.
            if outfmt == 'errors': return string.join(self.log.list['ErrorLog'],'\n')   # List of log error messages.
            if outfmt in ['fas','fasta']:
                fastxt = ''
                for entry in self.entries(): fastxt += entry.fasta()
                return fastxt
            if outfmt in self.dict['Output']:
                outdata = self.dict['Output'][outfmt]
                if outdata in self.info: outdata = self.info[outdata]
                if rje.exists(outdata):
                    nbytes = os.path.getsize(outdata)
                    if nbytes > maxfilesize > 0:
                        otext = '%s is too large to return (> %s)' % (outdata,rje.humanByteSize(nbytes))
                        if outfmt == self.getStrLC('Rest'): return 'ERROR: %s.' % otext
                        return '%s in full output. Try retrieve&rest=%s.' % (otext,outfmt)
                    return open(outdata,'r').read()
                else: return '%s\n' % outdata

            # summary = Reduced output
            # Default: return full text
            return '# SERVER ERROR: "&rest=%s" not a recognised output\n%s' % (outfmt,self.restFullOutput())
        except: return self.errorLog('Uniprot.restOutput(%s) error.' % outfmt)
#########################################################################################################################
    def restOutputOrder(self):
        outlist = []
        for outfmt in ['dat','map','failed','tab','ft','dom','pfam','xref','go']:
            if outfmt in self.dict['Output']: outlist.append(outfmt)
        return outlist
#########################################################################################################################
    def restOutputDefault(self): return 'dat'
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        dat = main Uniprot flat file. [dat]
        map = table of `uniprotid` elements given and the accession numbers they mapped to. [tdt]
        failed = list of identifiers that failed to map to uniprot. (Missing if none.) [list]
        tab = tabular summary of uniprot data (&tabout=T). [tdt]
        ft = table of parsed uniprot features (&ftout=T). [tdt]
        dom = table of parsed `DOMAIN` features and their positions (&domtable=T). [tdt]
        pfam = table of parsed `Pfam` domains and their counts for each protein (&pfamout=T). [tdt]
        xref = table of extracted database cross-references (&xrefout=T). [tdt]
        go = table of GO categories extracted for each protein (&gotable=T). [tdt]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ext = {}
            ### ~ [1] ~ Possible rest service outputs: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # dat = DAT download (or any alternative format using uniformat=X
            self.dict['Output']['dat'] = 'DATOut'; ext['dat'] = 'dat'
            # tab = Table of extracted UniProt details [None]
            self.dict['Output']['tab'] = 'TabOut'; ext['tab'] = 'tdt'
            # map = Table of mapped accession numbers [None]
            self.dict['Output']['map'] = 'MapOut'; ext['map'] = 'map.tdt'
            # ft = Outputs table of features into FILE [None]
            self.dict['Output']['ft'] = 'FTOut'; ext['ft'] = 'features.tdt'
            # xref = Table of extracted Database xref (Formerly linkout=FILE) [None]
            self.dict['Output']['xref'] = 'XRefOut'; ext['xref'] = 'xref.tdt'
            # pfam = Table of PFam domains from uniprot file [False]
            self.dict['Output']['pfam'] = '%s.pfam.tdt' % self.basefile()     #!# This needs modifying/updating
            # dom = Table of domains from uniprot file [False]
            self.dict['Output']['dom'] = '%s.domains.tdt' % self.basefile()     #!# This needs modifying/updating
            # go = Makes a table of AccNum:GO mapping [False]
            self.dict['Output']['go'] = '%s.go.tdt' % self.basefile()           #!# This needs modifying/updating
            # acc = Primary accession numbers
            # fas = Fasta format of entries [?]
            ## ~ [1a] ~ Update output files names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('Rest') in ext:
                skey = self.dict['Output'][self.getStr('Rest')]
                self.setStr({skey:'%s.%s' % (self.basefile(runpath=True),ext[self.getStr('Rest')])})
            elif self.getStr('Rest') in ['full']:
                for restout in ['tab','dat','ft','xref']:
                    skey = self.dict['Output'][restout]
                    self.setStr({skey:'%s.%s' % (self.basefile(runpath=True),ext[restout])})
            elif self.getStr('Rest') in ['default','seqin'] and not self.getStrLC('DATOut'):
                self.setStr({'DATOut':'%s.dat' % (self.basefile(runpath=True))})
            ### ~ [2] ~ Special REST output options for other programs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Rest') in ['seqin']: self.dict['Output'][self.getStr('Rest')] = 'DATOut'
        except: self.errorLog('Problem during Uniprot.restSetup()')
#########################################################################################################################
    def extractSpecies(self):   ### Uses index file to convert species codes into list of Accession Numbers
        '''Uses index file to convert species codes into list of Accession Numbers.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('UniPath').lower() in ['url','url/']:
                self.errorLog('Cannot make Taxa DAT file if unipath=url!',printerror=False)
                sys.exit()
            indexfile = self.info['UniPath'] + self.info['DBIndex']
            if not os.path.exists(indexfile):
                self.log.errorLog('Index file "%s" missing. Cannot make Taxa DAT file!' % indexfile,printerror=False)
                sys.exit()
            self.list['Extract'] = []
            ### ~ [1] ~ Read IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            line = INDEX.readline()
            self.progLog('\r#ACC','Extracting IDs for %d taxa: %s found.' % (len(self.list['SpecDat']),rje.integerString(len(self.list['Extract']))))
            while line:
                if rje.matchExp('(\S+_(%s));' % string.join(self.list['SpecDat'],'|'),line):
                    self.list['Extract'].append(rje.matchExp('(\S+_(%s));' % string.join(self.list['SpecDat'],'|'),line)[0])
                    pc = 100.0 * INDEX.tell() / end_pos
                    self.progLog('\r#ACC','Extracting IDs for %d taxa %.1f%%: %s found.' % (len(self.list['SpecDat']),pc,rje.integerString(len(self.list['Extract']))))
                line = INDEX.readline()
            INDEX.close()                    
            self.printLog('\r#ACC','Extracting IDs for %d taxa complete: %s found.' % (len(self.list['SpecDat']),rje.integerString(len(self.list['Extract']))))
        except: self.errorLog('Problem with %s.extractSpecies()' % self,quitchoice=False); raise
#########################################################################################################################
    ### <2> ### Reading Uniprot Entry
#########################################################################################################################
    def readUniProt(self,filename=None,clear=True,acclist=[],logft=False,use_index=True,cleardata=None,reformat=False):    ### Reads UniProt download into UniProtEntry objects
        '''
        Reads UniProt download into UniProtEntry objects.
        >> filename:str = UniProt download filename [None]
        >> clear:boolean = Whether to clear self.list['Entry'] before reading [True]
        >> acclist:list of str objects = UniProt accnum or id list to read
        >> logft:boolean [False] = whether to write number of features to log file
        >> use_index:boolean [True] = whether to use index file if present
        >> cleardata:boolean [None] = whether to clear processed data to save memory 
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << True if success, False if fail
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.getBool('ClearData')
            ## ~ [0a] ~ Index and DAT Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            indexfile = ''
            if self.info['DBIndex'].lower() not in ['','none'] and use_index: indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## ~ [0b] ~ Alternative File Name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not filename: filename = self.info['Name']
            if rje.exists(filename) and indexfile and not os.path.exists(indexfile):
                #self.log.printLog('#DB','Index file %s not found but DAT file given: ignoring index.' % indexfile)
                indexfile = ''
            if self.list['Proteome'] or self.list['Taxonomy']: filename = ''; indexfile = 'URL'
            elif not rje.exists(filename) and self.getStr('UniPath').lower() in ['','url','none/','url/']: indexfile = 'URL'
            elif not rje.exists(filename) and (not indexfile or not os.path.exists(indexfile)):
                self.printLog('#ERR','UniProt file "%s" and index file "%s" not found!' % (filename,indexfile))
                if rje.yesNo('Try to extract from UniProt website?'): indexfile = 'URL'
                else: return False
            ## ~ [0c] ~ Setup Extract list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not acclist: acclist = self.list['Extract'][0:]  
            acclist.sort()
            if not acclist and rje.exists(filename): indexfile = ''   # Process whole file
            if indexfile == 'URL':
                if self.list['Proteome']: return self._extractProteomesFromURL(self.list['Proteome'],cleardata=cleardata,logft=logft,reformat=reformat)
                elif self.list['Taxonomy']: return self._extractProteomesFromURL(self.list['Taxonomy'],cleardata=cleardata,logft=logft,reformat=reformat,taxonomy=True)
                elif acclist: return self._extractFromURL(acclist,cleardata=cleardata,logft=logft,reformat=reformat,onebyone=self.getBool('OneByOne'))
                self.printLog('#ERR','Cannot extract from URL without list of accession numbers or TaxaIDs')
                return False
            ## ~ [0d] ~ DatFiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datfiles = []
            if indexfile:
                datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
                if not datfiles:
                    self.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                    indexfile = ''
            if clear: self.list['Entry'] = []

            ### ~ [1] ~ Process whole file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not indexfile: return self._processWholeFile(filename,cleardata=cleardata,logft=logft,reformat=reformat)

            ### ~ [2] ~ Process from Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Index Dictionaries Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## ~ [2b] ~ Read and locate DataFiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.printLog('#SOURCE',line,screen=self.v()>0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                line = INDEX.readline()
            ## ~ [2c] Build index dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acclist:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc): splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else: missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else: dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else: acc_dict[acc] = [{key:pos}]   
                self.progLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            self.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            ## ~ [2d] ~ Splice Variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    
                if self.list['Extract']: self.list['Extract'] += splicevar
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else: dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else: acc_dict[acc] = [{key:pos}]   
                    self.progLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
                self.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            ## ~ [2e] ~ Missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for acc in missing: self.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### ~ [3] ~ Extract From DAT Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                if not rje.checkForFile(unifile): self.errorLog('%s entries skipped due to missing file' % unifile,printerror=False); continue     
                UNIFILE = open(unifile,'r')
                ex = 0; rx = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)   #?# [1]
                    # Read
                    _entry = self._readSingleEntry(UNIFILE,logft=logft,cleardata=cleardata,reformat=reformat,write=self.getStrLC('DATOut'))
                    if _entry == True: rx += 1 # Rejected: non-Reviewed entry etc.
                    elif _entry:
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        new_id = self.list['Entry'][-1].obj['Sequence'].info['ID']
                        new_acc = self.list['Entry'][-1].obj['Sequence'].info['AccNum']
                        for acc in wanted_acc:
                            tracc = string.split(acc,'_')[0]
                            if acc not in [new_id,new_acc] and tracc != new_acc:
                                self.printLog('\n#ACC','Secondary AccNum %s mapped to %s (%s).' % (acc,new_id,new_acc),screen=self.v() > 0)
                            try: self.db('map').addEntry({'uniprotid':acc,'accnum':new_acc},warn=False)
                            except: self.debug(self.db('map'))
                        # Update
                        ex += 1
                        self._entryTableData(self.list['Entry'][-1])
                        if self.opt['MemSaver']: self.list['Entry'] = []     # Delete for memsaver
                        #if self.dev() and ex >= 1000: self.warnLog('Dev break still in code!!',dev=True); break
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.progLog('\r#DAT','%s entries extracted from %s; %s rejected.' % (rje.iStr(ex),unifile,rje.iStr(rx)))
                UNIFILE.close()
                self.printLog('\r#DAT','%s entries extracted from %s; %s rejected.' % (rje.iStr(ex),unifile,rje.iStr(rx)))
            return True
        except: self.errorLog('UniProt.readUniProt() Failed. Check format.'); return False
#########################################################################################################################
    def _extractProteomesFromURL(self,taxalist,logft=True,cleardata=None,reformat=False,log=True,taxonomy=False,pickup=False):    ### Extracts acclist from UniProt website
        '''
        Extracts acclist from UniProt website.
        >> taxalist:list = List of TaxaIDs.
        >> logft:boolean = whether to write number of features to log file
        >> reformat:boolean = whether to save to DatOut using cut-down method
        >> log:bool [True] = whether to print results to Log
        >> taxonomy:bool [False] = whether to query taxonomy rather than organism
        >> pickup:bool [False] = whether to read/write to pickup file.
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            rx = 0      # Number of taxa read
            ex = 0      # Number of entries extracted
            fx = 0
            extracted = {}  # Taxa extracted successfully
            ## ~ [0a] Special Download Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uniformat = None; pfile = None
            if self.getStr('UniFormat') not in ['','none']:
                datout = self.getStr('DatOut')
                uniformat = self.getStr('UniFormat')
                if uniformat.endswith('gz'): uniformat = uniformat[:-2]; compress = True
                else: compress = False
                if datout.lower() in ['','none']:
                    if uniformat == 'txt': datout = '%s.dat' % (self.baseFile())
                    else: datout = '%s.%s' % (self.baseFile(),uniformat)
                if compress and datout[-3:] != '.gz': datout += '.gz'
                compress = {True:'yes',False:'no'}[datout[-3:] == '.gz']
                self.dict['Output']['dat'] = datout
                #self.deBug('%s (&compress=%s)' % (datout,compress))
                ## ~ [0a] Special pickup setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pfile = '%s.pickup' % datout
                if pickup:
                    if not rje.exists(pfile): pickup = False
                    else: pickup = rje.listFromCommand(pfile)
                elif rje.exists(pfile): os.unlink(pfile)
                if not pickup: rje.backup(self,datout)
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtext = 'Extracting proteomes from https://www.uniprot.org/uniprot/'
            if self.getBool('UseBeta'): logtext = string.replace(logtext,'www','beta')
            if self.opt['Complete']:
                logtext = string.replace(logtext,'proteomes','complete proteomes')
                self.warnLog('Specifying complete=T may not be working after Uniprot change. Use with caution.')
            if self.opt['Reviewed']: logtext = string.replace(logtext,'proteomes','reviewed proteomes')
            #self.deBug(logtext)
            #if log: self.progLog('\r#URL','%s ...' % logtext)
            ttot = len(taxalist); tx = 0
            self.printLog('#TAXID','%s TaxID for proteome extraction.' % (rje.iStr(ttot)))
            for taxid in taxalist:
                extracted[taxid] = 0; tx += 1
                if taxonomy: uniurl = 'https://www.uniprot.org/uniprot/?query=taxonomy:%s' % taxid
                else: uniurl = 'https://www.uniprot.org/uniprot/?query=taxonomy:%s+AND+proteome:*' % taxid
                try:
                    if self.getBool('UseBeta'): uniurl = string.replace(uniurl,'www','beta')
                    if self.opt['Complete']: uniurl += '+AND+keyword:"Complete"'
                    if self.opt['Reviewed']: uniurl += '+AND+reviewed:yes'
                    ## ~ [1a] Special download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if uniformat:   #!# This won't work for multiple compressed files!
                        if pickup and taxid in pickup: self.printLog('#SKIP','Skipping TaxID %s: found in %s.' % (taxid,pfile)); continue
                        uniurl += '&format=%s&compress=%s' % (uniformat, compress)
                        if log: self.printLog('#URL',uniurl,log=False)
                        if log: self.progLog('\r#URL','Downloading Taxa %s (%s) to %s (%s of %s)' % (taxid,uniformat,datout,rje.iStr(tx),rje.iStr(ttot)))
                        open(datout,'a').write(urllib2.urlopen(uniurl).read())
                        if log: self.printLog('\r#URL','Downloaded Taxa %s (%s) to %s (%s of %s).' % (taxid,uniformat,datout,rje.iStr(tx),rje.iStr(ttot)))
                        open(pfile,'a').write('%s\n' % taxid)
                        if self.getInt('SpecSleep') and taxid != taxalist[-1]:
                            try:
                                for sec in range(self.getInt('SpecSleep'),0,-1):
                                    self.progLog('\r#WAIT','SpecSleep for %d seconds!' % sec)
                                    time.sleep(1)
                                self.progLog('\r#WAIT','SpecSleep for 0 seconds!   ')
                            except: self.printLog('\r#WAIT','SpecSleep cancelled!       ',log=False)
                        continue
                    ## ~ [1b] Regular processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    acclist = string.split(urllib2.urlopen('%s&format=list' % uniurl).read())
                    self.printLog('#ACC','%s accession numbers identified for TaxaID %s proteome.' % (rje.iLen(acclist),taxid))
                    uniurl += '&format=txt'     # NB. Could use list to get accnum only
                    #!# Add option to download only as compressed file #!#
                    #self.deBug(uniurl)
                    UNIPROT = urllib2.urlopen(uniurl)
                    rx += 1
                    while self._readSingleEntry(UNIPROT,logft=logft,cleardata=cleardata,reformat=reformat):
                        self._entryTableData(self.list['Entry'][-1])
                        if self.opt['MemSaver']:
                            if self.list['Entry']: ex += 1; extracted[taxid] += 1
                            self.list['Entry'] = []     # Delete for memsaver
                        elif self.entryNum() != ex: ex += 1; extracted[taxid] += 1
                        if log: self.progLog('\r#URL','%s: %s read, %s extracted; %s failed.' % (logtext,rje.iStr(rx),rje.iStr(ex),rje.iStr(fx)))
                    UNIPROT.close()
                    if log: self.printLog('\r#PROT','%s of %s UniProt entries extracted for TaxaID %s Proteome.' % (rje.iStr(extracted[taxid]),rje.iLen(acclist),taxid))
                    if not extracted[taxid]: fx += 1
                except urllib2.HTTPError: self.errorLog('UniProt proteome "%s" not found!' % taxid); continue
                except ValueError: self.errorLog('Something went wrong parsing %s' % uniurl); fx += 1
                except: raise
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log and not uniformat: self.printLog('\r#URL','%s: %s read, %s extracted; %s failed.' % (logtext,rje.iStr(rx),rje.iStr(ex),rje.iStr(fx)))
            if rje.exists(pfile): os.unlink(pfile)
            return True     
        except: self.errorLog('UniProt._extractProteomesFromURL() Failed. Check format.'); return False
#########################################################################################################################
    def _extractProteinsFromURL(self,acclist,logft=True,cleardata=None,reformat=False,log=True):    ### Extracts acclist from UniProt website
        '''
        Extracts acclist from UniProt website.
        >> acclist:list = List of UniProt IDs and/or accession numbers
        >> logft:boolean = whether to write number of features to log file
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            rx = 0      # Number of entries read
            ## ~ [0a] ~ Strip splice variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            extract = []
            for acc in acclist:
                extract.append(string.split(acc,'-')[0])
            extract = rje.sortUnique(extract)
            if self.warn() and len(extract) != len(acclist): self.warnLog('AccList for _extractProteinsFromURL contains duplicates.')
            if self.warn() and len(extract) > 1e5: self.warnLog('Very long extraction list (%s identifiers) might fail!' % rje.iLen(extract))
            if not extract: self.printLog('#ERR','No accession numbers for Uniprot extraction!'); return False
            #!# Consider downloading in batches using while extract (with diminishing extract list)
            batchn = 100    # Number of entries to extract in one go
            batch = extract[0:]
            singles = []
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtext = 'Extracting %s AccNum from https://www.uniprot.org/uniprot/' % rje.iLen(extract)
            if log: self.progLog('\r#URL','%s...' % (logtext))
            while batch:
                tmpdat = '%s.%s.tmp' % (rje.baseFile(self.getStr('DATOut')),rje.randomString(6))
                #uniurl = "https://www.uniprot.org/uniprot/?query=%s&format=txt" % string.join(batch[:batchn],',')
                #!# Note, this will not maintain the order #!#
                uniurl = "https://www.uniprot.org/uniprot/?query=accession:%s&format=txt" % string.join(batch[:batchn],'+OR+accession:')
                #self.debug(uniurl)
                try: UNIPROT = urllib2.urlopen(uniurl)
                except:
                    self.errorLog('Uniprot URL access failure: %s' % uniurl,printerror=False)
                    try:
                        syscmd = 'curl -k -o %s "%s"' % (tmpdat,uniurl)
                        self.warnLog('Regular URL download failed. Will try to download using curl.')
                        self.printLog('#CURL','%s' % (syscmd))
                        os.system(syscmd)
                        if rje.exists(tmpdat):
                            UNIPROT = open(tmpdat,'r')
                        else:
                            # Try again
                            self.printLog('#WAIT','Uniprot curl failure: will try again in %d seconds' % self.getInt('SpecSleep'))
                            time.sleep(self.getInt('SpecSleep'))
                            os.system(syscmd)
                            if rje.exists(tmpdat):
                                UNIPROT = open(tmpdat,'r')
                            else:
                                raise ValueError(syscmd)
                    except:
                        self.printLog('#FAIL','UniProt._extractProteinsFromURL() failures. Will try one-by-one for failures.')
                        singles = batch[:batchn]
                        UNIPROT = urllib2.urlopen("https://www.uniprot.org/uniprot/%s.txt" % singles.pop(0))
                batch = batch[batchn:]
                try:
                    while self._readSingleEntry(UNIPROT,logft=logft,cleardata=cleardata,reformat=reformat,expect=False,write=True):
                        rx += 1     # Increment read entries
                        _entry = self.list['Entry'][-1]
                        [new_id,new_acc] = [_entry.obj['Sequence'].info['ID'],_entry.obj['Sequence'].info['AccNum']]
                        for acc in [new_id,new_acc]:
                            if acc in extract:
                                extract.remove(acc)
                                try: self.db('map').addEntry({'uniprotid':acc,'accnum':new_acc},warn=False)
                                except: self.debug(self.db('map'))
                            if acc in batch: batch.remove(acc)
                        try: secondaries = _entry.obj['Sequence'].list['Secondary ID']
                        except: secondaries = []
                        for acc in secondaries:
                            old_id = '%s_%s' % (acc, _entry.obj['Sequence'].info['SpecCode'])
                            if acc in batch: batch.remove(acc)
                            if old_id in batch: batch.remove(acc)
                            if acc in extract:
                                extract.remove(acc)
                                self.printLog('\n#ACC','Secondary AccNum %s mapped to %s (%s).' % (acc,new_id,new_acc),screen=self.v() > 0)
                                try: self.db('map').addEntry({'uniprotid':acc,'accnum':new_acc},warn=False)
                                except: self.debug(self.db('map'))
                            if old_id in extract:
                                extract.remove(old_id)
                                self.printLog('\n#ACC','Secondary AccNum %s mapped to %s (%s).' % (old_id,new_id,new_acc),screen=self.v() > 0)
                                try: self.db('map').addEntry({'uniprotid':old_id,'accnum':new_acc},warn=False)
                                except: self.debug(self.db('map'))
                        #self._writeSingleEntry(_entry,reformat=reformat)
                        self._entryTableData(_entry)
                        if self.opt['MemSaver']: self.list['Entry'] = []     # Delete for memsaver
                        if log: self.progLog('\r#URL','%s: %s read & extracted...' % (logtext,rje.iStr(rx)))
                        if singles:
                            UNIPROT.close()
                            UNIPROT = urllib2.urlopen("https://www.uniprot.org/uniprot/%s.txt" % singles.pop(0))
                except:
                    self.errorLog('Something went wrong parsing entries from https://www.uniprot.org/uniprot/?query=ACCLIST&format=txt')
                    #i# Extract list will now be longer
                x = UNIPROT.readlines()
                #self.debug('>>>\n%d\n%s\n???' % (len(x),string.join(x[-5:])))
                UNIPROT.close()
                if rje.exists(tmpdat): os.unlink(tmpdat)
                singles = []
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fx = len(extract)               # Number of AccNum failures
            dx = len(acclist) - (rx + fx)   # Number of duplicates ignored
            if log:
                self.printLog('\r#URL','%s: %s read & extracted; %s failed; %s duplicates.' % (logtext,rje.iStr(rx),rje.iStr(fx),rje.iStr(dx)))
                for acc in extract:
                    self.warnLog('Uniprot extraction for %s failed.' % acc,'extract_fail',suppress=True)
                    try:
                        if not self.db('map').data(acc):
                            self.db('map').addEntry({'uniprotid':acc,'accnum':'!FAILED!'},warn=False)
                    except: self.debug(self.db('map'))
            if extract: #x# and not self.getStrLC('Rest'):
                if self.getStrLC('DATOut'): mfile = rje.baseFile(self.getStr('DATOut')) + '.missing.acc'
                else: mfile = None
                if mfile and (self.i() < 0 or rje.yesNo('Save missing accnum to %s?' % mfile)):
                    open(mfile,'w').write(string.join(extract,'\n'))
                    self.printLog('#FAIL','%s missing accnum output to %s' % (rje.iLen(extract),mfile))
                    self.dict['Output']['failed'] = mfile
                if not self.getStrLC('Rest') and (self.i() < 0 or rje.yesNo('Try one-by-one for the %s failures?' % rje.iStr(fx))):
                    self.printLog('#FAIL','UniProt._extractProteinsFromURL() failures. Will try one-by-one for failures.')
                    return rx > 0 or self._extractFromURL(extract,logft,False,reformat,log,onebyone=True)
            return True
        except:
            self.errorLog('UniProt._extractProteinsFromURL() Failed. Check format. Will try one-by-one for full set.')
            rje.backup(self,self.getStr('DATOut'))  # Will skip if empty value
            return self._extractFromURL(acclist,logft,cleardata,reformat,log,onebyone=True)
#########################################################################################################################
    def _extractFromURL(self,acclist,logft=True,cleardata=None,reformat=False,log=True,onebyone=False):    ### Extracts acclist from UniProt website
        '''
        Extracts acclist from UniProt website.
        >> acclist:list = List of UniProt IDs and/or accession numbers
        >> logft:boolean = whether to write number of features to log file
        >> cleardata:bool [None] =
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(acclist) > 1 and not onebyone: return self._extractProteinsFromURL(acclist,logft,cleardata,reformat,log)
            #!# Can delete rest of this once _extractProteinsFromURL() is fully tested.
            if cleardata == None: cleardata = self.opt['ClearData']
            rx = 0      # Number of entries read
            ex = 0      # Number of entries extracted
            fx = 0      # Number of AccNum failures
            dx = 0      # Number of duplicates ignored
            extracted = {}  # Accession numbers extracted successfully
            primaries = []  # Primary Accession numbers (avoid multiple extraction)
            secondaries = []# Secondary Accession numbers (avoid multiple extraction)
            if self.warn() and len(rje.sortUnique(acclist)) != len(acclist): self.warnLog('AccList for _extractFromURL contains duplicates.')
            if not acclist: self.printLog('#ERR','No accession numbers for Uniprot extraction!'); return False
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtext = 'Extracting %s AccNum from https://www.uniprot.org/uniprot/' % rje.iLen(acclist)
            for acc in acclist[0:]:
                if self.getBool('SpliceVar'): acc = string.split(acc,'-')[0]
                if acc in primaries + secondaries: dx += 1; continue
                extracted[acc] = 0
                try:
                    #!# Update to use https://www.uniprot.org/uniprot/?query=P40688,O43521&format=txt
                    #!# Will need to handle splice variants differently. - Safest to strip prior to download
                    UNIPROT = urllib2.urlopen("https://www.uniprot.org/uniprot/%s.txt" % acc)
                    while self._readSingleEntry(UNIPROT,logft=logft,cleardata=cleardata,reformat=reformat,expect=not extracted[acc],write=False,fulltext=True):
                        # Check for replaced AccNum/ID
                        rx += 1
                        extracted[acc] += 1
                        new_id = self.list['Entry'][-1].obj['Sequence'].info['ID']
                        new_acc = self.list['Entry'][-1].obj['Sequence'].info['AccNum']
                        try: self.db('map').addEntry({'uniprotid':acc,'accnum':new_acc},warn=False)
                        except: self.debug(self.db('map'))
                        if acc not in [new_id,new_acc]: self.printLog('\n#ACC','Secondary AccNum %s mapped to %s (%s).' % (acc,new_id,new_acc),screen=self.v() > 0)
                        if new_acc in primaries: self.list['Entry'].pop(-1); dx += 1 # Duplicate!
                        else:
                            ex += 1; primaries.append(new_acc)
                            _entry = self.list['Entry'][-1]
                            self._writeSingleEntry(_entry,reformat=reformat)
                            try: secondaries += _entry.obj['Sequence'].list['Secondary ID']
                            except: pass
                            self._entryTableData(self.list['Entry'][-1])
                            if self.opt['MemSaver']: self.list['Entry'] = []     # Delete for memsaver
                    UNIPROT.close()
                    if not extracted[acc]: fx += 1; raise ValueError
                    if log: self.progLog('\r#URL','%s: %s read, %s extracted; %s failed; %s duplicates.' % (logtext,rje.iStr(rx),rje.iStr(ex),rje.iStr(fx),rje.iStr(dx)))
                except urllib2.HTTPError: self.errorLog('UniProt entry "%s" not found! Possibly withdrawn or demerged' % acc); fx += 1; continue
                except ValueError: self.errorLog('Something went wrong parsing https://www.uniprot.org/uniprot/%s.txt' % acc); fx += 1
                except: raise
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log: self.printLog('\r#URL','%s: %s read, %s extracted; %s failed; %s duplicates.' % (logtext,rje.iStr(rx),rje.iStr(ex),rje.iStr(fx),rje.iStr(dx)))
            return True
        except: self.errorLog('UniProt._extractFromURL() Failed. Check format.'); return False
#########################################################################################################################
    def _processWholeFile(self,filename,logft=True,cleardata=None,reformat=False,log=True,url=False):    ### Processes whole file into entries
        '''
        Processes whole file into entries.
        >> filename:str = UniProt filename
        >> logft:boolean = whether to write number of features to log file
        >> reformat:boolean = whether to save to DatOut using cut-down method
        >> url:bool [False] = whether filename is actually a URL query to be parsed.
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            rx = 0      # Number of entries read
            ex = 0      # Number of entries extracted
            if url: UNIPROT = urllib2.urlopen(filename)
            else: UNIPROT = open(filename, 'r')
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtext = 'Extracting entries from %s' % string.split(filename,os.sep)[-1]
            while self._readSingleEntry(UNIPROT,logft=logft,cleardata=cleardata,reformat=reformat):
                rx += 1
                self._entryTableData(self.list['Entry'][-1])
                if self.opt['MemSaver'] and self.entryNum() > 0:    # Kept sequence
                    ex += 1
                    self.list['Entry'] = []     # Delete for memsaver
                else: ex = self.entryNum()
                if log: self.progLog('\r#DAT','%s: %s read, %s extracted.' % (logtext,rje.integerString(rx),rje.integerString(ex)))
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log: self.printLog('\r#DAT','%s: %s read, %s extracted.' % (logtext,rje.integerString(rx),rje.integerString(ex)))
            UNIPROT.close()
            return True     
        except: self.errorLog('Cataclysmic error during %s._processWholeFile()!' % self); return False
#########################################################################################################################
    def _newEntryObject(self):
        ecmd = self.cmd_list
        # Add Pfam to DBParse list if required for PfamOut
        if self.getStrLC('PfamOut') and 'pfam' not in self.list['DBList']: ecmd.append('dblist=%s,pfam' % string.join(self.list['DBList'],','))
        return UniProtEntry(log=self.log,cmd_list=ecmd,parent=self)
#########################################################################################################################
    def _writeSingleEntry(self,_entry,reformat=False): ### Handles writing of entry data from _readSingleEntry
        '''Handles writing of entry data from _readSingleEntry.'''
        if self.getStrLC('DATOut'):
            #self.debug(_entry.info['ID'])
            if self.info['SplitOut'].lower() not in ['','none']:
                datapp = False
                datout = rje.makePath(self.info['SplitOut']) + self.info['DATOut']
                rje.mkDir(self,datout)
                datout = string.join([datout,_entry.info['Type'],_entry.info['ID'],'dat'],'.')
            else: datout = self.info['DATOut']; datapp = True
            #self.debug(datout)
            #self.debug(datapp)
            if reformat: self.saveUniProt(datout,[_entry],append=datapp)
            else: open(datout,{True:'a',False:'w'}[datapp]).write(_entry.info['FullText'])
            #!# >>>>>> This is a fudge but it's OK for now <<<<<<<<<<< #!#
            #X#self.tableOutput(_entry)
            for k in ['DBLinks']:
                if not self.list.has_key(k): self.list[k] = []
            for k in ['DBLinks']:
                for i in _entry.dict[k].keys():
                    if i not in self.list[k]: self.list[k].append(i)
            #!# ^^^^^^^^^^ This is a fudge but it's OK for now ^^^^^^^ #!#
            _entry.info['FullText'] = ''
#########################################################################################################################
    def _readSingleEntry(self,UNIPROT,logft=True,cleardata=None,reformat=False,expect=False,write=True,fulltext=False):    ### Reads a single entry from current point in file and processes
        '''
        Processes whole file into entries.
        >> UNIPROT:FileHandle = Open UniProt file for reading *at start of entry*
        >> logft:boolean = whether to write number of features to log file
        >> cleardata:boolean = whether to clear data to save memory after processing
        >> reformat:boolean = whether to save to DatOut using cut-down method
        >> expect:boolean = whether to raise an error if not found
        >> write:boolean [True] = whether to write entry to DATOut if given. If False, _entry['FullText'] is not stored.
        << returns _entry/False dependent on success
        '''
        try:### ~ [0] ~ Check for end of file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fulltext = fulltext or write
            if cleardata == None: cleardata = self.opt['ClearData']
            line = UNIPROT.readline()
            if not line:
                if expect: raise ValueError('No (more) lines read from <UNIPROT>!')
                return False
            revskip = False     # Whether to speedskip an unreviewed entry
            ### ~ [1] ~ Process Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _entry = None       # Current entry
            _reading = False    # Whether currently reading an entry
            while line:
                #self.bugPrint(line)
                ## ~ [1a] ~ Full Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if _entry and fulltext: _entry.info['FullText'] = '%s%s' % (_entry.info['FullText'],line)
                line = rje.chomp(line)
                ltype = line[0:2]
                if revskip and ltype == '//': return True
                elif revskip or ltype in ['','XX']: line = UNIPROT.readline(); continue
                elif self.list['UParse'] and ltype not in self.list['UParse'] + ['ID','AC','DE','GN','SQ','  ','//']:
                    line = UNIPROT.readline(); continue
                rest = line[5:]
                ## ~ [1b] ~ New Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                    
                if ltype == 'ID':
                    _reading = True
                    _entry = self._newEntryObject() 
                    if fulltext: _entry.info['FullText'] = '%s\n' % line
                elif not _entry:
                    self.errorLog('Expected ID entry, got "%s"' % line,printerror=False)
                    raise ValueError
                ## ~ [1c] ~ End of Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ltype == '//':
                    if self._add_entry(_entry,acclist=self.list['Extract']):     # Entry is one of the desired entries
                        self.list['Entry'].append(_entry)   # Add entry to list
                        _entry.process(logft=logft,cleardata=cleardata and not reformat)         # Extracts details from uniprot
                        if self.opt['Reviewed'] and not _entry.info['Type'] == 'Reviewed':
                            self.list['Entry'].remove(_entry)
                            return True
                        if not write: return _entry
                        self._writeSingleEntry(_entry,reformat=reformat)
                        return _entry
                    else:
                        try: self.warnLog('Entry not found in extract list: %s' % _entry.dict['Data']['ID'],warntype='bad_entry',suppress=True)
                        except: pass
                        self.deBug('Bad Entry'); self.deBug(_entry); self.deBug(self.list['Extract']); return None
                ## ~ [1d] ~ Entry Details read into a dictionary within the entry ~~~~~~~~~~~~~~~~~ ##
                if not rest: self.deBug(line)
                if _entry.dict['Data'].has_key(ltype):   # Append list
                    if rest[:1] != ' ': _entry.dict['Data'][ltype].append(rest)   # New entry
                    else:
                        while rest[:1] == ' ': rest = rest[1:]
                        _entry.dict['Data'][ltype][-1] = '%s %s' % (_entry.dict['Data'][ltype][-1], rest)
                elif ltype == '  ':
                    if _entry.dict['Data'].has_key('SEQ'):
                        _entry.dict['Data']['SEQ'][0] = '%s %s' % (_entry.dict['Data']['SEQ'][0],rest)
                    elif _entry.dict['Data'].has_key('SQ'): _entry.dict['Data']['SEQ'] = [rest]
                elif ltype not in ['XX']: _entry.dict['Data'][ltype] = [rest]
                if self.opt['FullRef'] and ltype in ['RN','RP','RC','RX','RA','RT','RL','RG']: _entry.list['References'].append(line)
                ## ~ [1e] ~ Quick skip if reviewed constraint not met ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.opt['Reviewed'] and ltype == 'ID':
                    parse = _entry._uniParse('ID')
                    if parse and parse[1] != 'Reviewed': revskip = True
                ## ~ [1f] ~ Next Line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                line = UNIPROT.readline()
            ### ~ [2] ~ Reached EOF before // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if _reading: self.errorLog('Started UniProt Entry but EOF reached before "//". Truncated input file?',printerror=False)
            return False                
        except: self.errorLog('UniProt._readSingleEntry() parsing error!',quitchoice=False); raise
#########################################################################################################################
    def _add_entry(self,_entry,acclist):    ### Returns True if _entry in acclist or false if not
        '''
        Returns True if _entry in acclist or false if not.
        >> _entry:uniProtEntry object (unprocessed)
        >> acclist:list of accession numbers and/or IDs
        '''
        try:### ~ Check for entry in acclist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not acclist: return True
            elif self.getStr('UniPath').lower() in ['url','url/']: return True
            for acc in acclist[0:]:
                if self.getBool('SpliceVar'): acc = string.split(acc,'-')[0] 
                if _entry.dict['Data'].has_key('AC') and string.join([' '] + _entry.dict['Data']['AC'],' ').find(' %s;' % acc) > 0: return True
                if _entry.dict['Data'].has_key('AC') and string.join([' '] + _entry.dict['Data']['AC'],' ').find(' %s;' % string.split(acc,'_')[0]) > 0: return True
                if _entry.dict['Data'].has_key('ID') and string.join(_entry.dict['Data']['ID'],' ').find('%s ' % acc) == 0: return True
        except: self.errorLog('Cataclysmic error during _add_entry')
        return False
#########################################################################################################################
    def addFromSeq(self,seq=None,sequence='',name='',data={},ft=[]):    ### Converts into UniProtEntry object 
        '''
        Converts into UniProtEntry object and adds to self.
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns entry if successful or None if fails
        '''
        ### ~ [1] ~ Add case features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        addft = ft[0:]
        if self.info['UCFT'] and seq:
            for uc in seq.dict['Case']['Upper']:
                addft.append({'Type':self.info['UCFT'].upper(),'Desc':self.info['UCFT'],'Start':uc[0]+1,'End':uc[1]+1})
        if self.info['LCFT'] and seq:
            for lc in seq.dict['Case']['Lower']:
                addft.append({'Type':self.info['LCFT'].upper(),'Desc':self.info['LCFT'],'Start':lc[0]+1,'End':lc[1]+1})
        ### ~ [2] ~ Make Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
        newentry = self._newEntryObject()
        entry = newentry.uniProtFromSeq(seq,sequence,name,data,addft)
        if entry: self.list['Entry'].append(entry)
        return entry
#########################################################################################################################
    ### <3> ### Uniprot Info Output
#########################################################################################################################
    def _entryTableData(self,entry):    ### Updates db xref and ft db tables from entry data
        '''Updates db xref and ft db tables from entry data.'''
        try:### ~ [0] Setup Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            fdb = self.db('ft')
            ### ~ [1] DB XRef Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xdb:     # ['accnum','db','xref'],['accnum','db']
                #self.debug(entry.dict['DB'])
                for db in entry.dict['DB']:
                    if db in ['Pfam']:
                        idlist = []
                        for dbid in rje.sortKeys(entry.dict['DB'][db]): idlist.append('%s; %s' % (dbid,entry.dict['DB'][db][dbid]))
                        xentry = xdb.addEntry({'accnum':entry.obj['Sequence'].info['AccNum'],'db':db,'xref':string.join(idlist,'|')})
                    else: xentry = xdb.addEntry({'accnum':entry.obj['Sequence'].info['AccNum'],'db':db,'xref':string.join(rje.sortUnique(entry.dict['DB'][db]),'|')})
                    #self.debug(xentry)
            ### ~ [2] Features Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fdb:     # ['accnum','feature','ft_start','ft_end','description'],['accnum','feature','ft_start','ft_end']
                for ft in entry.list['Feature']:
                    fdb.addEntry({'accnum':entry.obj['Sequence'].info['AccNum'],'feature':ft['Type'],'ft_start':ft['Start'],'ft_end':ft['End'],'description':ft['Desc']})
        except: self.errorLog('Problem updating db tables for %s' % entry)
#########################################################################################################################
    def tableOutput(self):   ### Tabulated output of UniProt information
        '''Tabulated output of UniProt information. Divided into TabOut (UniProt summary) and XRefOut (database links)'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if (self.db('xref') or self.db('ft')) and not self.getBool('MemSaver'):
            #    for entry in self.entries(): self._entryTableData(self,entry)

            ### ~ [1] ~ Output xref and feature tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Database Links Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('GOTable'): self.goTable()
            if self.getStrLC('XRefOut') or self.getBool('DBSplit') or self.getStrLC('PfamOut'): self.xrefOutput()
            ## ~ [1b] ~ Feature/Domain output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('FTOut') or self.getBool('DomTable'): self.ftOutput()
            ## ~ [1c] ~ Mapping information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                if self.getStrLC('DATOut') or self.getStrLC('TabOut'):
                    if self.db('map').entryNum(): self.db('map').saveToFile(self.getStr('MapOut'))
            except: self.errorLog('Problem saving AccNum map table.')

            ### ~ [2] ~ Special Table Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['TabOut'].lower() in ['','none']:
                return


            #!# Add secondary accnum! #!#

            ### Setup Groups and Columns ###
            headers = {'#1# Basic Details #':['#','AccNum','ID','Len','Description','Gene','Species'],
                       '#2# Function & Activity #':['Function','GO_MF','GO_BP','Activity','Interactions','Phenotype',
                                                    'Similarity'],
                       '#3# Expression, Location & Structure #':['Tissue','Cell_Loc','PDB','InterPro','Pfam','PROSITE',
                                                                 'Isoforms','Ensembl'],
                       '#4# References & Links #':['GeneCards','PubMed','Keywords','Comments']
                       }

            ### Open File and write header ###
            delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['TabOut']))
            TABOUT = open(self.info['TabOut'],'a')   # Already deleted if append=F
            if not self.getStrLC('Rest'):
                TABOUT.write('# Generated by %s: %s\n' % (self.log.info['Name'],time.asctime(time.localtime(time.time()))))
                if self.info['Name'] != 'None':
                    TABOUT.write('# Source: %s\n' % os.path.abspath(self.info['Name']))
                else:
                    TABOUT.write('#Source: %s\n' % os.path.abspath(self.info['UniPath']+'*.dat'))
                TABOUT.write('# Seqnum: %s\n\n' % rje.integerString(self.entryNum()))
            ## Headers ##
            head1 = []
            head2 = []
            for h in rje.sortKeys(headers):
                head1.append(h)
                head1 += [''] * (len(headers[h])-1)
                head2 += headers[h]
            if not self.getStrLC('Rest'):
                rje.writeDelimit(TABOUT,head1,delimit)  #!# Consider removing these or making them an option? #!#
            rje.writeDelimit(TABOUT,head2,delimit)

            ### Write Data for Entries ###
            ex = 0
            for entry in self.list['Entry']:
                seq = entry.obj['Sequence']
                ex += 1
                data = []
                comments = rje.sortKeys(entry.dict['Comments'])
                for h in head2:     # Column headers:
                    #1# Basic Details #
                    if h == '#':
                        data.append(rje.preZero(ex,self.entryNum()))
                    elif h in ['AccNum','ID','Description']:
                        data.append(seq.info[h])
                    elif h == 'Len':
                        data.append('%d' % seq.aaLen())
                    elif h == 'Gene':
                        data.append(string.join([seq.info['Gene']] + entry.list['Synonyms'],'; '))
                    elif h == 'Species':
                        data.append('%s [%s]' % (seq.info['Species'],seq.info['SpecCode']))
                    #2# Function & Activity #
                    elif h == 'Function':
                        text = ''
                        for cc in ['FUNCTION','PATHWAY','DOMAIN']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'GO_MF':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; F:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'GO_BP':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; P:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'Activity':   #!# Join to Function #!#
                        text = ''
                        for cc in ['CATALYTIC ACTIVITY','COFACTOR']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Interactions':
                        text = ''
                        for cc in ['INTERACTION','ENZYME REGULATION','SUBUNIT','PTM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Phenotype':
                        text = ''
                        for cc in ['DISEASE','POLYMORPHISM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Similarity':
                        text = ''
                        for cc in ['SIMILARITY']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #3# Expression, Location & Structure #
                    elif h == 'Tissue':
                        text = ''
                        if entry.list['Tissues']:
                            text = 'TISSUES: %s' % string.join(entry.list['Tissues']+[''],'; ')
                        for cc in ['TISSUE SPECIFICITY','DEVELOPMENTAL STAGE','INDUCTION']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Cell_Loc':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; C:') < 0:
                                    go.remove(g)
                        cc = 'SUBCELLULAR LOCATION'
                        if entry.dict['Comments'].has_key(cc):
                            comments.remove(cc)
                            go = entry.dict['Comments'][cc] + go
                        data.append(string.join(go,' >> '))
                    elif h in ['PDB','InterPro','Pfam','PROSITE','Ensembl']:
                        if entry.dict['DBLinks'].has_key(h):
                            data.append(string.join(entry.dict['DBLinks'][h],' >> '))
                        else:
                            data.append('')
                    elif h == 'Isoforms':
                        text = ''
                        for cc in ['ALTERNATIVE PRODUCTS']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #4# References & Links
                    elif h == 'GeneCards':
                        if entry.dict['DBLinks'].has_key('HGNC'):
                            data.append(string.join(entry.dict['DBLinks']['HGNC'],' >> '))
                        else:
                            data.append('')
                    elif h in ['PubMed']:
                        if entry.list[h]:
                            data.append('http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=%s' % string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h in ['Keywords']:
                        if entry.list[h]:
                            data.append(string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h == 'Comments':
                        text = ''
                        for cc in comments:
                            if text:
                                text += ' '
                            text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                            if text[-1] != '.':
                                text += '.'
                        data.append(text)

                rje.writeDelimit(TABOUT,data,delimit)
                self.log.printLog('\r#OUT','UniProt summary output: %.1f%%' % (100.0 * ex / self.entryNum()),log=False,newline=False)
            self.log.printLog('\r#OUT','UniProt summary output for %s entries.' % rje.integerString(self.entryNum()))
            TABOUT.close()

        except: self.errorLog('Major error during UniProt.tableOutput()!',quitchoice=True)
#########################################################################################################################
    def goTable(self):   ### Delimited output of UniProt:GO database links
        '''Delimited output of UniProt:GO database links.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            if not xdb.entryNum(): self.printLog('#XREF','No GO xref parsed for output!')
            gdb = self.db().addEmptyTable('go',['accnum','go','type','desc','evidence','source'],['accnum','go'])
            ### ~ [1] ~ Add Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #DR   GO; GO:0021987; P:cerebral cortex development; IMP:MGI.
            #DR   GO; GO:0021766; P:hippocampus development; IMP:MGI.
            #DR   GO; GO:0035308; P:negative regulation of protein dephosphorylation; IDA:MGI.
            for entry in xdb.indexEntries('db','GO'):
                for gox in string.split(entry['xref'],'|'):
                    godata = string.split(gox,'; ')
                    gentry = {'accnum':entry['accnum'],
                              'go':string.split(godata[0],':')[1],
                              'type':string.split(godata[1],':')[0],
                              'desc':string.split(godata[1],':')[1],
                              'evidence':string.split(godata[2],':')[0],
                              'source':string.split(godata[2],':')[1]}
                    gdb.addEntry(gentry)
            gdb.saveToFile()
            ### ~ [2] ~ Memsaver clearup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb.dropEntriesDirect('db',['GO'])
            if self.getBool('MemSaver'): self.db().deleteTable('go')
        except: self.errorLog('Major error during UniProt.xrefOutput()!',quitchoice=True)
#########################################################################################################################
    def xrefOutput(self):   ### Delimited output of UniProt database links
        '''Delimited output of UniProt database links.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            if self.getStrLC('PfamOut'):
                pfamdata = xdb.subset('db','Pfam',copy=False)
                self.debug(pfamdata)
            if 'pfam' not in self.list['DBList'] and 'Pfam' in xdb.index('db'): xdb.dropEntriesDirect('db',['Pfam'])
            if not self.getStrLC('XRefOut'): self.printLog('#XREF','No xref output (xrefout=F/None).')
            elif not xdb.entryNum(): self.printLog('#XREF','No database xref parsed for output!')
            ### ~ [1] ~ Split output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elif self.getBool('DBSplit'): self.linkDBSplitOutput()
            ### ~ [2] ~ Simple output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elif self.getBool('LongLink'): xdb.saveToFile(self.getStr('XRefOut'))
            ### ~ [3] ~ Wide output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                ## ~ [3a] ~ Setup wide output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                xdb.index('accnum')
                dhead = ['Uniprot'] + rje.sortKeys(xdb.index('db'))
                rje.backup(self,self.getStr('XRefOut'))
                delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['XRefOut']))
                rje.delimitedFileOutput(self,self.getStr('XRefOut'),dhead,delimit)
                LINKFILE = open(self.info['XRefOut'],'a')
                ## ~ [3b] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ex = 0.0; etot = len(xdb.index('accnum'))
                for acc in rje.sortKeys(xdb.index('accnum')):
                    datadict = {'Uniprot':acc}
                    for db in xdb.index('db'): datadict[db] = ''
                    for entry in xdb.indexEntries('accnum',acc): datadict[entry['db']] = entry['xref']
                    ex += 1
                    rje.delimitedFileOutput(self,LINKFILE,dhead,delimit,datadict)
                    self.progLog('\r#OUT','Database Links output: %.1f%%' % (ex/etot)); ex += 100.0
                self.printLog('\r#OUT','Database Links output for %s entries.' % rje.integerString(self.entryNum()))
                LINKFILE.close()
            ### ~ [4] ~ Pfam Table Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('PfamOut'):
                pfamdb = self.db().addEmptyTable('pfam',['accnum','pfam','name','num'],['accnum','pfam'])
                acclist = []
                for (acc,db) in rje.sortKeys(pfamdata):
                    pfamdata[acc] = pfamdata.pop((acc,db))
                    acclist.append(acc)
                self.debug(pfamdata)
                for entry in self.entries():
                    acc = entry.accNum()
                    if acc not in acclist: acclist.append(acc)
                if self.db('map'): acclist = rje.listUnion(acclist,self.db('map').indexKeys('accnum'))
                acclist.sort()
                if '!FAILED!' in acclist: acclist.remove('!FAILED!')
                self.debug(acclist)
                for acc in acclist:
                    if acc in pfamdata:
                        for pfamdom in string.split(pfamdata.pop(acc)['xref'],'|'):
                            domdata = string.split(pfamdom,'; ')
                            pentry = {'accnum':acc,'pfam':domdata[0],'name':pfamdom,'num':0}
                            if len(domdata) > 1: pentry['name'] = domdata[1]
                            if len(domdata) > 2: pentry['num'] = int(domdata[2][:-1])
                            pfamdb.addEntry(pentry)
                    else: pfamdb.addEntry({'accnum':acc,'pfam':'','name':'No Pfam domains parsed.','num':0})
                pfamdb.saveToFile(self.getStr('PfamOut'))
            ### ~ [5] ~ Memsaver clearup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MemSaver'): self.db().deleteTable('xref')
        except: self.errorLog('Major error during UniProt.xrefOutput()!',quitchoice=True)
#########################################################################################################################
    def linkOutput(self):   ### Delimited output of UniProt database links                                      #!#OLD#!#
        '''Delimited output of UniProt database links.'''
        try:
            ### Setup ##
            if self.getBool('DBSplit'): return self.linkDBSplitOutput()    # Replace with rje_db tables at some point
            if self.opt['LongLink']: return self.linkOutputLong()  # self.obj['Sequence'].list['Secondary ID']
            delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['XRefOut']))
            dblist = rje.dictValues(self.list,'DBLinks')
            dblist.sort()   #!# Group Databases by type later #!#
            
            ### Open File and write header ###
            LINKFILE = open(self.info['XRefOut'],'a')   # Already deleted if append=F
            LINKFILE.write('# Generated by %s: %s\n' % (self.log.info['Name'],time.asctime(time.localtime(time.time()))))
            if self.info['Name'] != 'None':
                LINKFILE.write('# Source: %s\n' % os.path.abspath(self.info['Name']))
            else:
                LINKFILE.write('# Source: %s\n' % os.path.abspath(self.info['UniPath']+'*.dat'))
            LINKFILE.write('# Seqnum: %s\n\n' % rje.integerString(self.entryNum()))
            rje.writeDelimit(LINKFILE,['#','AccNum']+dblist,delimit)
            
            ### Write Data for Entries ###
            ex = 0
            for entry in self.list['Entry']:
                ex += 1
                data = [rje.preZero(ex,self.entryNum()),entry.obj['Sequence'].info['AccNum']]
                for db in dblist:
                    if entry.dict['DBLinks'].has_key(db):
                        data.append(string.join(entry.dict['DBLinks'][db],' >> '))
                    else:
                        data.append('')
                rje.writeDelimit(LINKFILE,data,delimit)
                self.log.printLog('\r#OUT','Database Links output: %.1f%%' % (100.0 * ex / self.entryNum()),log=False,newline=False)
            self.log.printLog('\r#OUT','Database Links output for %s entries.' % rje.integerString(self.entryNum()))
            LINKFILE.close()
            
        except: self.log.errorLog('Major error during UniProt.linkOutput()!',quitchoice=True)
#########################################################################################################################
    def linkDBSplitOutput(self):   ### Delimited output of UniProt database links
        '''Delimited output of UniProt database links.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            #if not self.list['DBList']: self.list['DBList'] = rje.sortKeys(xdb.index('db'))
            for db in rje.sortKeys(xdb.index('db')):    #self.list['DBList']:
                #if db not in xdb.index('db'):
                #    self.printLog('#DB','No %s xref to output.' % db)
                #    continue
                if db == 'GO' and self.getBool('GOTable'): continue
                if db == 'Pfam': dbhead = ['Uniprot',db,'Domain','N']
                else: dbhead = ['Uniprot',db]
                dfile = '%s.%s.tdt' % (self.basefile(),db)
                rje.backup(self,dfile)
                OUTFILE = open(dfile,'a')
                rje.writeDelimit(OUTFILE,dbhead)
            ### ~ [1] ~ Write Data for Entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                ex = 0
                for entry in xdb.indexEntries('db',db):
                    ex += 1
                    for dbid in string.split(entry['xref'],'|'):
                        data = {'Uniprot':entry['accnum'],db:dbid}
                        if db == 'Pfam':
                            #self.debug(entry)
                            [data[db],data['Domain'],data['N']] = string.split(data[db],'; ')
                            data['N'] = data['N'][:-1]
                            rje.delimitedFileOutput(self,OUTFILE,dbhead,datadict=data,delimit='\t')
                        else: OUTFILE.write('%s\t%s\n' % (entry['accnum'],dbid))
                OUTFILE.close()
                self.printLog('#DB','%s %s xref output to %s.' % (rje.iStr(ex),db,dfile))
            self.printLog('\r#OUT','Database Links output complete (dbsplit=T).')
        except: self.errorLog('Major error during UniProt.linkOutput()!',quitchoice=True)
#########################################################################################################################
    def linkOutputLong(self):   ### Delimited output of UniProt database links                                  #!#OLD#!#
        '''Delimited output of UniProt database links.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            listhead = ['AccNum','DBase','LinkAcc']
            rje.delimitedFileOutput(self,self.info['XRefOut'],listhead,rje_backup=True)
            ### ~ [1] ~ Write Data for Entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0
            for entry in self.list['Entry']:
                ex += 1
                for secid in entry.obj['Sequence'].list['Secondary ID']:
                    data = {'AccNum':entry.obj['Sequence'].info['AccNum'],'DBase':'UniProt','LinkAcc':secid}
                    rje.delimitedFileOutput(self,self.info['XRefOut'],listhead,datadict=data)
                for db in rje.sortKeys(entry.dict['DBLinks']):
                    for dbid in entry.dict['DBLinks'][db]:
                        data = {'AccNum':entry.obj['Sequence'].info['AccNum'],'DBase':db,'LinkAcc':dbid}
                        rje.delimitedFileOutput(self,self.info['XRefOut'],listhead,datadict=data)
                self.progLog('\r#OUT','Database Links output: %.1f%%' % (100.0 * ex / self.entryNum()))
            self.printLog('\r#OUT','Database Links output for %s entries.' % rje.integerString(self.entryNum()))
        except: self.errorLog('Major error during UniProt.linkOutput()!',quitchoice=True)            
#########################################################################################################################
    def ftOutput(self): ### Outputs feature and/or domain tables.
        '''Outputs feature and/or domain tables.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('ft')
            if not fdb.entryNum(): self.printLog('#FT','No features (or domains) parsed for output!')
            ### ~ [1] ~ Feature Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('FTOut'): fdb.saveToFile(self.getStr('FTOut'))
            ### ~ [2] ~ Domain Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DomTable'):
                if not self.getBool('MemSaver'): fdb = self.db().copyTable(fdb,'domains')
                else: fdb.setStr({'Name':'domains'})
                for entry in fdb.entries():
                    if entry['feature'].upper() != 'DOMAIN': fdb.dropEntry(entry)
                fdb.dropField('feature')
                fdb.renameField('description','domain')
                fdb.renameField('ft_start','dom_start')
                fdb.renameField('ft_end','dom_end')
                fdb.list['Fields'] = ['accnum','domain','dom_start','dom_end']
                fdb.saveToFile()
                self.db().deleteTable(fdb)
            elif self.getBool('MemSaver'): self.db().deleteTable(fdb)
        except: self.errorLog('Cataclysmic error during ftOutput!')
#########################################################################################################################
    def domTable(self): return self.ftOutput() ### Outputs domain info into a table                    #!#OLD#!#    #V2.9
#########################################################################################################################
    def ftTable(self,outfile):  ### Outputs features into a table                   #!#OLD#!# Used by motiflist/slimlist!
        '''
        Outputs features into a table.
        >> outfile:str = Name of output file
        '''
        try:
            ### Setup ###
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            if self.opt['Append']:
                FT = open(outfile,'a')
            else:
                FT = open(outfile,'w')
                rje.writeDelimit(FT,['acc_num','feature','ft_start','ft_end','description'],delimit)

            ### Output ###
            (fx,ex) = (0,0.0)
            for entry in self.list['Entry']:
                ex += 100.0
                acc = entry.obj['Sequence'].info['AccNum']
                ## Make dictionary of {start:{end:[features]}}
                ft_dict = {}
                for ft in entry.list['Feature']:
                    ft_start = ft['Start']
                    if not ft_dict.has_key(ft_start):
                        ft_dict[ft_start] = {}
                    ft_end = ft['End']
                    if not ft_dict[ft_start].has_key(ft_end):
                        ft_dict[ft_start][ft_end] = []
                    ft_dict[ft_start][ft_end].append(ft)
                ## Sort and output ##
                for ft_start in rje.sortKeys(ft_dict):
                    for ft_end in rje.sortKeys(ft_dict[ft_start]):
                        for ft in ft_dict[ft_start][ft_end]:
                            outlist = [acc]
                            for fk in ['Type','Start','End','Desc']:
                                outlist.append('%s' % ft[fk])
                            rje.writeDelimit(FT,outlist,delimit)
                        fx += 1
                self.log.printLog('\r#FT','Feature output: %.1f%% (%s features)' % (ex/len(self.list['Entry']),rje.integerString(fx)),log=False,newline=False)

            ### End ###
            FT.close()
            self.log.printLog('\r#FT','Feature output complete: %s features, %s entries.' % (rje.integerString(fx),rje.integerString(len(self.list['Entry']))))
        except: self.log.errorLog('Program error during rje_uniprot.ftTable()',quitchoice=True)
#########################################################################################################################
    def saveUniProt(self,outfile,entries=[],append=False):    ### Saves self as a DAT file
        '''
        Saves self as a DAT file.
        >> outfile:str = Name of output file
        >> entries:list of entries (self.list['Entry'] if none given)
        >> append:boolean = whether to append file
        '''
        try:### ~ [0] ~ Setup Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not append: rje.backup(self,outfile)
            if not entries: entries = self.list['Entry'][0:]
            ### ~ [1] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT = open(outfile,'a')                    
            for entry in entries[0:]:
                if not entry.dict['Data'] and not entry.uniProtFromSeq():
                    entries.remove(entry)
                    self.errorLog('Problem with %s (%s) - cannot output' % (entry,entry.info['Name']),printerror=False)
                    continue
                seq = entry.obj['Sequence']
                ## ~ [1a] ~ Standard info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for key in ['ID','AC','DT','DE','GN','OS']:
                    if entry.dict['Data'].has_key(key):
                        for rest in entry.dict['Data'][key]: OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## ~ [1b] ~ Other data, except Features and sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for key in rje.sortKeys(entry.dict['Data']):
                    if key not in ['ID','AC','DT','DE','GN','OS','FT','SQ','SEQ','//']:
                        for rest in entry.dict['Data'][key]: OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## ~ [1c] ~ Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry.orderFT()
                for ftdict in entry.list['Feature']:
                    (p1,p2) = (ftdict['Start'],ftdict['End'])
                    ftxt = 'FT   %s' % ftdict['Type']
                    while len(ftxt) < 14 or ftxt[-1] != ' ': ftxt += ' '
                    ftxt += '%6s' % ('%d' % p1)
                    while len(ftxt) > 20 and ftxt[-(len('%d' % p1)+2):-len('%d' % p1)] == '  ': ftxt = ftxt[:-(len('%d' % p1)+1)] + ftxt[-len('%d' % p1):]
                    ftxt += '%7s' % ('%d' % p2)
                    while len(ftxt) > 27 and ftxt[-(len('%d' % p2)+2):-len('%d' % p2)] == '  ': ftxt = ftxt[:-(len('%d' % p2)+1)] + ftxt[-len('%d' % p2):]
                    ftxt += ' %s\n' % ftdict['Desc']
                    OUT.write(ftxt)
                ## ~ [1d] ~ Sequence/End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if seq.dna(): OUT.write('SQ   SEQUENCE%s%d BP;  XXX MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen()))
                else: OUT.write('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen(),rje_sequence.MWt(seq.info['Sequence'])))
                uniseq = seq.info['Sequence'][0:]
                while len(uniseq) > 0:
                    OUT.write('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                    uniseq = uniseq[60:]
                OUT.write('//\n')
            OUT.close()
            if not append or len(entries) > 1: self.printLog('#OUT','UniProt format for %d entries saved to %s' % (len(entries),outfile))            
        except: self.errorLog('Major problem with %s.saveUniProt()' % self)
#########################################################################################################################
    ### <4> ### Additional UniProt Tools                                                                                #
#########################################################################################################################
    def accNameSeq(self,acc_list=[],spec=None,justsequence=True):  ### Method to extract dictionaries of {acc:'ID__PrimaryAcc Desc'} & {acc:seq}
        '''
        Method to extract dictionary of {acc:'ID__PrimaryAcc Desc'} & {acc:seq} using index.
        >> acclist:list of accession numbers. Will use self.list['Extract'] if none given
        >> spec:Limit to species code
        >> justsequence:bool [True] = Whether to just return the sequence data (True) or a sequence object (False)
        << tuple of dictionaries of ({acc:ID__PrimaryAcc},{acc:uniprot sequence})
        '''
        try:
            ### Setup ###
            accshortname = {}
            accseq = {}
            ## Index and DAT Files ##
            indexfile = ''
            if self.info['DBIndex'] != 'None':
                indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## Extract list ##
            if acc_list:
                acc_list.sort() #X#
            else:
                acc_list = self.list['Extract'][0:]  #X# 
                acc_list.sort()
            ## DatFiles ##
            datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
            if not datfiles:
                self.log.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                return {}

            ### Process from Index ###
            _stage = 'Index Dictionaries (Setup)'
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## Read DataFiles ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.verbose(1,3,line,0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.log.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                #X#print rje.fileLineFromSeek(INDEX,start_pos,reseek=False,next=False)
                line = INDEX.readline()

            ### Build index dictionary ###
            _stage = 'Build Index Dictionaries'
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acc_list:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc):
                        splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else:
                        missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else:
                    dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else:                    
                    acc_dict[acc] = [{key:pos}]   
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
            self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Splice Variants ###
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    #X#                splicevar = rje.sortUnique(splicevar)
                if self.list['Extract']:
                    #self.deBug(self.list['Extract'])
                    self.list['Extract'] += splicevar
                    #self.deBug(self.list['Extract'])
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else:
                        dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else:                    
                        acc_dict[acc] = [{key:pos}]   
                    self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Missing acc have no entry in returned dictionary ###
            #X#for acc in missing:            
            #X#    self.log.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### Extract From UniProt using dictionaries ###
            _stage = 'Extract using Dictionaries'
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                UNIFILE = open(unifile,'r')
                ex = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)#[1]
                    # Read
                    if self._readSingleEntry(UNIFILE,logft=False,cleardata=True):
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        new_id = self.list['Entry'][-1].obj['Sequence'].info['ID']
                        new_acc = self.list['Entry'][-1].obj['Sequence'].info['AccNum']
                        new_desc = self.list['Entry'][-1].obj['Sequence'].info['Description']
                        if spec and self.list['Entry'][-1].obj['Sequence'].info['SpecCode'] != spec:
                            self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                            continue                            
                        #X#print new_id, new_acc, wanted_acc
                        for acc in wanted_acc:
                            accshortname[acc] = '%s__%s %s' % (new_id,new_acc,new_desc)
                            if justsequence: accseq[acc] = self.list['Entry'][-1].obj['Sequence'].info['Sequence']
                            else: accseq[acc] = self.list['Entry'][-1].obj['Sequence']
                        # Update
                        ex += 1
                        self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.log.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile),log=False,newline=False)
                UNIFILE.close()
                self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile))

            return (accshortname,accseq)
            
        except:
            self.log.errorLog('Program error during rje_uniprot.accShortName()',quitchoice=True)
#########################################################################################################################
    def accDictFromEntries(self,acc_list=[],empties=False,remdup=True):   ### Generates dictionary of {acc:UniProtEntry} from self.list['Entries']
        '''
        Generates dictionary of {acc:UniProtEntry} from self.list['Entries']. Will be all if no acc_list given.
        >> acc_list:list [] = Optional subset of accnum as list
        >> empties:bool [False] = whether to return None for missing acc in acc_list if not found (otherwise missing)
        >> remdup:bool [True] = whether to remove duplicate accession numbers that point to 2+ entries
        << dictionary of {acc:UniProtEntry} for all accnum (including secondary)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            accdict = {}
            dupwarn = []     # Possibly add an option to reject duplicately mapped AccNum
            accmap = {}
            for acc in acc_list: accmap[acc] = acc.split('-')[0]    # Strip splice variants
            ### ~ [1] Extract ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in self.entries():
                e_acc = string.split(entry.info['Name'],'__')
                try: e_acc += entry.obj['Sequence'].list['Secondary ID']
                except: pass
                for acc in e_acc:
                    #self.deBug(acc in acc_list)
                    if acc_list and acc not in accmap.values(): continue
                    if acc in accdict:
                        if acc not in dupwarn: dupwarn.append(acc)
                        self.printLog('\r#DUP','Dup = %s from %s already mapped to %s' % (acc,entry.shortName(),accdict[acc].shortName()),screen=self.v() > 0)
                    else: accdict[acc] = entry
            if dupwarn:
                if remdup:
                    for dupacc in dupwarn: accdict.pop(dupacc)
                    self.printLog('\r#WARN','%s duplicate (one-to-many) accession numbers removed from mapping.' % rje.iLen(dupwarn))
                else: self.printLog('\r#WARN','%s duplicate (one-to-many) accession numbers! (Single entry mapping.)' % rje.iLen(dupwarn))
            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if acc_list and empties:
                for acc in acc_list:
                    if acc not in accdict: accdict[acc] = None
            return accdict
        except: self.errorLog('UniProt.accDictFromEntries() failure.'); raise
#########################################################################################################################
    def accDict(self,acc_list=[],cleardata=None):      ### Method to extract dictionaries of {acc:UniProtEntry}
        '''
        Method to extract dictionaries of {acc:UniProtEntry}.
        >> acclist:list of accession numbers. Will use self.list['Extract'] if none given
        << dictionary of {acc:UniProtEntry}
        '''
        try:
            ### Setup ###
            if cleardata == None: cleardata = self.opt['ClearData']
            accentry = {}
            ## Index and DAT Files ##
            indexfile = ''
            if self.info['DBIndex'] != 'None': indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## Extract list ##
            if acc_list: acc_list.sort() #X#
            else:
                acc_list = self.list['Extract'][0:]  #X# 
                acc_list.sort()
            ## DatFiles ##
            datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
            if not datfiles:
                self.log.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                return {}

            ### Process from Index ###
            _stage = 'Index Dictionaries (Setup)'
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## Read DataFiles ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.verbose(1,3,line,0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.log.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                #X#print rje.fileLineFromSeek(INDEX,start_pos,reseek=False,next=False)
                line = INDEX.readline()

            ### Build index dictionary ###
            _stage = 'Build Index Dictionaries'
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acc_list:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc):
                        splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else:
                        missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else:
                    dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else:                    
                    acc_dict[acc] = [{key:pos}]   
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
            self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Splice Variants ###
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    #X#                splicevar = rje.sortUnique(splicevar)
                if self.list['Extract']:
                    #self.deBug(self.list['Extract'])
                    self.list['Extract'] += splicevar
                    #self.deBug(self.list['Extract'])
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else:
                        dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else:                    
                        acc_dict[acc] = [{key:pos}]   
                    self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Missing acc have no entry in returned dictionary ###
            #X#for acc in missing:            
            #X#    self.log.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### Extract From UniProt using dictionaries ###
            _stage = 'Extract using Dictionaries'
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                UNIFILE = open(unifile,'r')
                ex = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)#[1]
                    # Read
                    if self._readSingleEntry(UNIFILE,logft=False,cleardata=cleardata):
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        for acc in wanted_acc:
                            accentry[acc] = self.list['Entry'][-1]
                        # Update
                        ex += 1
                        self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.log.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile),log=False,newline=False)
                UNIFILE.close()
                self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile))

            return (accentry)
            
        except:
            self.log.errorLog('Program error during rje_uniprot.accEntry()',quitchoice=True)
#########################################################################################################################
## End of SECTION III: UniProt Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: UniProtEntry Class                                                                                      # 
#########################################################################################################################
class UniProtEntry(rje.RJE_Object):     
    '''
    UniProt Entry Class. Author: Rich Edwards (2005).

    Info:str
    - Name = UniProt ID of Entry
    - Type = Preliminary or Standard
    - FullText = Full Text of Entry
    
    Opt:boolean
    - CC2FT = Extra whole-length features added for TISSUE and LOCATION [False]
    - ClearData = Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]
    - InvMask = Whether to invert the masking and only retain maskft features [False]    
    - FullRef = Whether to store full Reference information in UniProt Entry objects [False]
    - TMConvert = Whether to convert TOPO_DOM features, using first description word as Type [False]

    Stat:numeric
    - Length = Length of Sequence as annotated

    List:list
    - CaseFT = List of Features to make upper case with rest of sequence lower case []
    - Feature = List of feature dictionaries: [Type,Start,End,Desc]
    - MaskFT = List of Features to mask out []
    - PubMed = List of PubMed IDs (as strings)
    - Keywords = List of UniProt Keywords
    - Tissues = List of UniProt Tissues
    - References = List of reference dictionaries [RP,RC,RX,RA,RT,RL,RG]
    - Synonyms = List of Gene synonyms
    
    Dict:dictionary
    - Data = Dictionary of lists of UniProt data (Keys are line headers ID/AC/CC etc.)
    - DB = Specific extractions from DR lines for use in other programs. {DB:[AccNum/ID]}
    - Comments = Dictionary of comments: {Type:List of Comments}
    - DBLinks = List of Database Link dictionaries {Dbase,List of Details} for dblinks output

    Obj:RJE_Objects
    - Sequence = rje_sequence.Sequence object
    '''
    ### Attributes
    def seqobj(self): return self.obj['Sequence']
    def gene(self): return self.seqobj().getStr('Gene')
    def id(self): return self.seqobj().getStr('ID')
    def acc(self): return self.seqobj().getStr('AccNum')
    def accNum(self): return self.seqobj().getStr('AccNum')
    def shortName(self): return '%s (%s)' % (self.id(),self.acc())
    def seqi(self,ikey): return self.seqobj().getStr(ikey)
    def reviewed(self): return self.getStr('Type') == 'Reviewed'
    def length(self):
        if self.getInt('Length'): return self.getInt('Length')
        else: return self.seqobj().aaLen()
#########################################################################################################################
    def fasta(self): return '>%s__%s %s\n%s\n' % (self.id(),self.accNum(),self.seqi('Description'),self.seqi('Sequence'))
    def seqname(self): return '%s__%s %s' % (self.id(),self.accNum(),self.seqi('Description'))
    def sequence(self): return self.seqi('Sequence')
#########################################################################################################################
    def isSwissprot(self):  ### Return whether or not entry is perceived to by Swissprot (as opposed to Trembl)
        uid = self.id()
        return '_' in uid and uid.upper() == uid and not uid.startswith(self.acc())
#########################################################################################################################
    def isSpecies(self,spec=None,speclist=[]):  ### Returns True if entry corresponds to listed species
        '''
        Returns True if entry corresponds to listed species (or species code).
        >> spec:str = single species that MUST be be right
        >> speclist:list = can match any species in list
        '''
        if spec: return spec in [self.obj['Sequence'].info['Species'],self.obj['Sequence'].info['SpecCode']]
        for sp in speclist:
            if sp in [self.obj['Sequence'].info['Species'],self.obj['Sequence'].info['SpecCode']]: return True
        return False
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','FullText']
        self.statlist = ['Length']
        self.optlist = ['CC2FT','InvMask','TMConvert','ClearData','FullRef']
        self.listlist = ['Feature','PubMed','Keywords','Tissues','Synonyms','MaskFT','CaseFT','References','DBParse']
        self.dictlist = ['Data','Comments','DBLinks','DB']
        self.objlist = ['Sequence']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.obj['Sequence'] = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list)
        self.info['FullText'] = ''
        self.list['Feature'] = []   # List of features = {'Type':str,'Start':int,'End':int,'Desc':str}
        self.list['MaskFT'] = []
        self.list['DBParse'] = dbparse[0:]
        self.setOpt({'ClearData':True})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'opt',['CC2FT','InvMask','TMConvert','ClearData','FullRef'])
                self._cmdReadList(cmd,'list',['CaseFT','MaskFT','DBParse'])
                self._cmdRead(cmd,'list','DBParse','dblist')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self.list['DBParse'] = rje.listLower(self.list['DBParse'])
#########################################################################################################################
    ### <2> ### Attribute Processing
#########################################################################################################################
    def _uniParse(self,key):    ### Parses list of elements from self.dict['Data'] (and pops) 
        '''
        Parses list of elements from self.dict['Data'] (and pops).
        >> key:str = Key of UniProt entry type
        << List of matched elements or False if failure.
        '''
        try:### ~ [1] Check key and parse data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if key not in self.dict['Data'].keys(): return False
            return rje.matchExp(uniparse[key],self.dict['Data'][key][0])
            #!# This no longer pops elements from dictionary: check backwards compatibility #!#
        except:
            self.log.errorLog('UniProtEntry: Cataclysmic error during _uniParse(%s)!' % key)
            return False
#########################################################################################################################
    def process(self,logft=True,cleardata=None):  ### Extract Details from self.dict['Data'] to Sequence object
        '''
        Extract Details from self.dict['Data'] to Sequence object.
        >> logft:boolean = whether to write number of features to log file
        >> cleardata=Whether to clear self.dict['Data'] after processing (to save memory) [True]
        >> uparse:list = Restricted list of data keys to
        << True if OK, False if not.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            if cleardata == None: cleardata = self.opt['ClearData']
            seq_info = ['Name','Type','Description','Sequence','ID','AccNum','DBase','Gene','Species','SpecCode','Format','TaxaID']
            seqi = self.obj['Sequence'].info
            for i in seq_info:
                if i not in seqi: seqi[i] = ''
            seqi['Type'] = 'Protein'
            #X#print '\nAcc:', self.dict['Data']['AC']
            _stage = 'Compress Certain uniprot lists'
            for key in ['DE','GN','OC','AC','RA','RT']:
                if self.dict['Data'].has_key(key): self.dict['Data'][key] = [string.join(self.dict['Data'][key])]
                                    
            ### Basic Sequence Details ###
            _stage = 'Basic Details (ID)'
            parse = self._uniParse('ID')
            if parse:
                (self.info['Name'],self.info['Type'],self.stat['Length']) = parse[:3]
                self.stat['Length'] = string.atoi(self.stat['Length'])
            else: self.stat['Length'] = -1
            self.info['ID'] = seqi['ID'] = self.info['Name']

            _stage = 'AccNum (AC)'
            full_acc = string.split(string.join(self.dict['Data']['AC']))
            self.obj['Sequence'].list['Secondary ID'] = string.split(string.join(full_acc,''),';')[1:-1]
            parse = self._uniParse('AC')
            if parse: seqi['AccNum'] = parse[0]
            if string.split(self.info['Name'])[0] == seqi['AccNum']:
                seqi['ID'] = seqi['AccNum']
                seqi['DBase'] = 'custom'
            elif self.info['Name'].find(seqi['AccNum']) == 0: seqi['DBase'] = 'trembl'
            else: seqi['DBase'] = 'sprot'
            if seqi['ID'] != seqi['AccNum']: self.info['Name'] = '%s__%s' % (seqi['ID'],seqi['AccNum'])

            _stage = 'Description (DE)'
            parse = self._uniParse('DE')
            if parse: seqi['Description'] = parse[0]
                
            _stage = 'Gene (GN)'
            parse = self._uniParse('GN')
            if parse: seqi['Gene'] = parse[0]
            if seqi['Gene'][-1:] == ';': seqi['Gene'] = seqi['Gene'][:-1]
            if self.dict['Data'].has_key('GN'):
                if rje.matchExp(uniparse['SY'],self.dict['Data']['GN'][0]):
                    syn = string.split(rje.matchExp(uniparse['SY'],self.dict['Data']['GN'][0])[0],';')[0]
                    self.list['Synonyms'] = string.split(syn,', ')
                    
            _stage = 'Species (OS)'
            parse = self._uniParse('OS')
            if parse: seqi['Species'] = parse[0]
            parse = self._uniParse('OX')
            if parse: seqi['TaxaID'] = parse[0]
            seqi['SpecCode'] = seqi['ID'][seqi['ID'].find('_')+1:]
            
            _stage = 'Name'
            if seqi['DBase'] == 'trembl' and seqi['Gene'] != 'None':
                for g in seqi['Gene'][0:]:
                    if not rje.matchExp('([A-Za-z0-9_-])',g) and g not in ['.','#']: seqi['Gene'] = string.replace(seqi['Gene'],g,'')
                seqi['ID'] = '%s_%s' % (seqi['Gene'].lower(),seqi['SpecCode'])
            if seqi['ID'] == seqi['AccNum']:
                seqi['Name'] = '%s %s' % (seqi['ID'],seqi['Description'])
                seqi['Format'] = 'gnspec'
            else:
                seqi['Name'] = '%s__%s %s' % (seqi['ID'],seqi['AccNum'],seqi['Description'])
                seqi['Format'] = 'gn_sp__acc'
            
            _stage = 'Sequence'
            if self.dict['Data'].has_key('SEQ'): seqi['Sequence'] = self.dict['Data'].pop('SEQ')[0]
            seqi['Sequence'] = re.sub('\s+','',seqi['Sequence']).upper()
            if not self.stat['Length']: self.stat['Length'] = len(seqi['Sequence'])

            ### Comments (CC) ###
            _stage = 'Comments'
            if self.dict['Data'].has_key('CC'):
                self.dict['Comments'] = {}  # Dictionary of comments: {Type:List of Comments}
                for cc in self.dict['Data']['CC']:
                    if cc.find('-----') == 0: break
                    csplit = string.split(cc[4:],': ')
                    ctype = csplit[0]
                    cdetail = string.join(csplit[1:],': ')
                    if self.dict['Comments'].has_key(ctype): self.dict['Comments'][ctype].append(cdetail)
                    else: self.dict['Comments'][ctype] = [cdetail]

            ### Features ###
            _stage = 'Features'
            if self.dict['Data'].has_key('FT'):
                for ft in self.dict['Data']['FT']:
                    # Remove '?'
                    while rje.matchExp('(\?\d)',ft) or rje.matchExp('(\d\?)',ft): ft = string.replace(ft,'?','')
                    ft = string.replace(ft,'?','0')
                    parse = rje.matchExp(uniparse['FT'],ft)
                    parse_nodesc = rje.matchExp('(\S+)\s+<*(\d+)\s+>*(\d+)\.*',ft)
                    parse_onepos = rje.matchExp('(\S+)\s+(\d+)\.*\s+(\S+)',ft)
                    parse_newpos = rje.matchExp('(\S+)\s+<*(\d+)\.\.>*(\d+) /note="(.+)"',ft)
                    parse_newid = rje.matchExp('(\S+)\s+<*(\d+)\.\.>*(\d+) .*/id="(.+)"',ft)
                    parse_newev = rje.matchExp('(\S+)\s+<*(\d+)\.\.>*(\d+) /evidence="(.+)"',ft)
                    parse_newnull = rje.matchExp('(\S+)\s+<*(\d+)\.\.>*(\d+)',ft)
                    parse_cntd = rje.matchExp('\s+(\S.*)$',ft)
                    if parse:
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : parse[3]
                            }
                        if rje.matchExp(string.join(['(\S+)','<(\d+)','>*(\d+)\.*','(\S.+)\s*$'], '\s+'),ft) or rje.matchExp(string.join(['(\S+)','<*(\d+)','>(\d+)','(\S.+)\s*$'], '\s+'),ft):
                            ftdic['Desc'] = ftdic['Desc'] + ' (Truncated?)'
                        ftdic['Desc'] = re.sub('\s+',' ',ftdic['Desc'])
                        if self.opt['TMConvert'] and ftdic['Type'] == 'TOPO_DOM':
                            ftdic['Type'] = string.strip(string.split(ftdic['Desc'])[0].upper(),'.')
                            ftdic['Desc'] = 'TOPO_DOM %s' % ftdic['Desc']
                        self.list['Feature'].append(ftdic)
                    elif parse_nodesc:
                        parse = parse_nodesc
                        if parse:
                            ftdic = {
                                'Type' : parse[0],
                                'Start' : string.atoi(parse[1]),
                                'End' : string.atoi(parse[2]),
                                'Desc' : parse[0]
                                }
                            if rje.matchExp('(\S+)\s+<(\d+)\s+>*(\d+)\.*',ft) or rje.matchExp('(\S+)\s+<*(\d+)\s+>(\d+)\.*',ft):
                                ftdic['Desc'] = ftdic['Desc'] + ' (Truncated?)'
                            self.list['Feature'].append(ftdic)
                    elif parse_onepos:
                        parse = parse_onepos
                        if parse:
                            ftdic = {
                                'Type' : parse[0],
                                'Start' : string.atoi(parse[1]),
                                'End' : string.atoi(parse[1]),
                                'Desc' : parse[2]
                                }
                            self.list['Feature'].append(ftdic)
                    elif parse_newpos:
                        parse = parse_newpos
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : parse[3]
                            }
                        self.list['Feature'].append(ftdic)
                    elif parse_newid:
                        parse = parse_newid
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : parse[3]
                            }
                        self.list['Feature'].append(ftdic)
                    elif parse_newev:
                        parse = parse_newev
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : 'evidence='+parse[3]
                            }
                        self.list['Feature'].append(ftdic)
                    elif parse_newnull:
                        parse = parse_newnull
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : ''
                            }
                        self.list['Feature'].append(ftdic)
                    elif parse_cntd:
                        try:
                            ftdic = self.list['Feature'][-1]
                            ftdic['Desc'] = re.sub('\s+',' ','%s %s' % (ftdic['Desc'],parse_cntd[0]))
                        except: self.printLog('#ERR','No feature parsed for addition of "%s".' % ft)
                    else: self.log.printLog('#ERR','Cannot parse feature details from %s.' % ft)

            ### VAR_SEQ (FT) ###
            _stage = 'Splice Variants'
            splicevar = {}      # Dictionary read from comments of {isoform:IsoID}
            if 'ALTERNATIVE PRODUCTS' in self.dict['Comments']:
                #self.bugPrint('\nALTERNATIVE PRODUCTS\n')
                #self.bugPrint(self.dict['Comments']['ALTERNATIVE PRODUCTS'][0])
                for isoform in self.dict['Comments']['ALTERNATIVE PRODUCTS'][0].split('Name=')[1:]:
                    #self.bugPrint(isoform)
                    try:
                        if rje.matchExp('^(\S)\s',isoform):
                            splicevar[rje.matchExp('^(\S)\s',isoform)[0]] = rje.matchExp('IsoId=(\S+)[;,]',isoform)[0]
                        elif rje.matchExp('^(\S+)\s\{',isoform):
                            splicevar[rje.matchExp('^(\S+)\s\{',isoform)[0]] = rje.matchExp('IsoId=(\S+)[;,]',isoform)[0]
                        else: splicevar[isoform.split(';')[0]] = rje.matchExp('IsoId=(\S+)[;,]',isoform)[0]
                    except: self.deBug('\nIsoform recognition failure...\n%s' % isoform)
                #self.bugPrint(splicevar)
            if self.list['Feature']:
                # Generate a list of VAR_SEQ tuples
                varseq = {}
                for ftdic in self.list['Feature']:
                    if ftdic['Type'] != 'VAR_SEQ': continue
                    #self.bugPrint(ftdic['Desc'])
                    for isoform in string.split(ftdic['Desc'],'isoform')[1:]:
                        #self.bugPrint(isoform)
                        isoform = string.replace(isoform,' and','')
                        try: var = splicevar[string.replace(isoform,' ','')]
                        except:
                            #self.bugPrint('%s' % rje.matchExp('^\s*(\d+)[\W]',isoform))
                            try: var = splicevar[rje.matchExp('^\s*(\d+)[\W]',isoform)[0]]
                            except:
                                var = None
                                for isomap in splicevar:
                                    var = rje.matchExp('^\s*(%s)\W' % rje.strEscape(isomap,'()[]+'),isoform)
                                    if not var: var = rje.matchExp('^(%s)\W' % rje.strEscape(isomap,'()[]+'),string.replace(isoform,' ',''))
                                    if var: var = splicevar[isomap]; break
                                if not var:
                                    self.deBug('\nIsoform mapping failure...\n%s' % isoform)
                                    self.deBug(splicevar)
                                    continue
                        if var not in varseq: varseq[var] = []
                        if ftdic['Desc'].startswith('Missing'): varseq[var].append((ftdic['Start']-1,ftdic['End'],''))
                        else:
                            vardat = rje.matchExp('^(\S+)->(\w+)\(',string.replace(ftdic['Desc'],' ',''))
                            if vardat:
                                if seqi['Sequence'][ftdic['Start']-1:ftdic['End']] != vardat[0]:
                                    self.errorLog('%s-%s: VAR_SEQ replacement sequence error' % (seqi['AccNum'],var),printerror=False); varseq[var].append('')
                                else: varseq[var].append((ftdic['Start']-1,ftdic['End'],vardat[1]))
                            else: self.errorLog('Failed to parse seqvar from %s' % ftdic,printerror=False); varseq[var].append('')
                # Generate splice variant sequences
                for var in varseq:
                    seqi[var] = seqi['Sequence']
                    while varseq[var]:
                        try:
                            vardat = varseq[var].pop(-1)    # Work backwards to keep sequence positions OK
                            #self.deBug('...%s %s %s...' % (seqi[var][:vardat[0]][-10:], vardat[2], seqi[var][vardat[1]:][:10]))
                            seqi[var] = seqi[var][:vardat[0]] + vardat[2] + seqi[var][vardat[1]:]                            
                        except: seqi.pop(var); break
                seq = self.obj['Sequence']
                for isoform in splicevar:   # Store the reverse mapping in seq.dict to get isoform name from AccNum
                    if splicevar[isoform] in seqi: seq.dict['SpliceVar'][splicevar[isoform]] = isoform

            ### Tissues (RC) ###
            _stage = 'Tissues'
            if self.dict['Data'].has_key('RC'):
                tissues = []
                self.list['Tissues'] = []
                for rc in self.dict['Data']['RC']:
                    parse = rje.matchExp(uniparse['RC'],rc)
                    if parse: tissues += string.split(parse[0],', ')
                for tissue in tissues:
                    if tissue[:4] == 'and ': self.list['Tissues'].append(tissue[4:])
                    else: self.list['Tissues'].append(tissue)

            ### Keywords (KW) ###
            _stage = 'Keywords'
            if self.dict['Data'].has_key('KW'):
                keywords = string.join(self.dict['Data']['KW'])
                if keywords[-1:] == '.': keywords = keywords[:-1]    # Remove full stop
                self.list['Keywords'] = string.split(keywords,'; ')

            ### References (RX) ###
            _stage = 'PubMed IDs'
            if self.dict['Data'].has_key('RX'):
                self.list['PubMed'] = []
                for rx in self.dict['Data']['RX']:
                    parse = rje.matchExp(uniparse['RX'],rx)
                    if parse: self.list['PubMed'].append(parse[0])
            if self.list['References']:
                newreflist = []
                while self.list['References']:
                    rline = self.list['References'].pop(0)
                    type = string.split(rline)[0]
                    rline = string.join(string.split(rline)[1:])
                    if type == 'RN': newreflist.append({'RN':rje.matchExp(uniparse['RN'],rline)[0]})
                    elif type in newreflist[-1]: newreflist[-1][type] += ' %s' % rline
                    else: newreflist[-1][type] = rline
                self.list['References'] = newreflist
                #self.deBug(self.list['References'])

            ### Database Links ### (DR)
            _stage = 'Database Links'
            if self.dict['Data'].has_key('DR'):
                self.dict['DBLinks'] = {}   # Database Link dictionary {Dbase,List of Details}
                for dr in self.dict['Data']['DR']:
                    if rje.matchExp(uniparse['DR'],dr):
                        (ctype,cdetail) = rje.matchExp(uniparse['DR'],dr)
                    elif rje.matchExp('^(Entrez Gene); (\S+)$',dr):
                        (ctype,cdetail) = rje.matchExp('^(Entrez Gene); (\S+)$',dr)
                    else: continue
                    #if self.dbToParse(ctype) or self.dbLinks(ctype): self.parseDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']
                    if self.dbToParse(ctype): self.parseDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']
                    #if not self.dbLinks(ctype): continue
                    #if self.dict['DBLinks'].has_key(ctype): self.dict['DBLinks'][ctype].append(cdetail)
                    #else: self.dict['DBLinks'][ctype] = [cdetail]
            if self.dbToParse('id'): self.parseDB('ID',self.id())
            if self.dbToParse('gene'): self.parseDB('Gene',self.gene())


            ### CC to FT ###
            _stage = 'End'
            if self.opt['CC2FT']: self.cc2ft()

            ### FT Masking/Case Change ###
            if self.list['MaskFT']: self.maskFT(self.list['MaskFT'],inverse=self.opt['InvMask'])
            if self.list['CaseFT']: self.caseFT(self.list['CaseFT'])
            #X#for k in ['Synonyms','PubMed','Keywords','Tissues']:
            #X#    print k, self.list[k]
            #X#for k in ['Comments','DBLinks']:
            #X#    print k, self.dict[k]         

            ### Cleanup ###
            #X#for key in self.dict['Data'].keys():
            #X#    if key not in useful_data:
            #X#        self.dict['Data'].pop(key)      # Save memory!
            #X#print self.dict['Data']
            #X#print self.info, self.obj['Sequence'].info
            if cleardata: self.dict['Data'] = {}    # Save memory!
            else: self.obj['Sequence'].dict['UniDAT'] = self.dict['Data']
            if logft:
                self.printLog('#FT','%d features for %s.' % (len(self.list['Feature']),self.obj['Sequence'].info['AccNum']))
            return True
        except:
            self.log.errorLog('Cataclysmic error during UniProtEntry: process() %s!' % _stage)
            return False
#########################################################################################################################
    def dbToParse(self,dbase):  ### Whether to parse db or not
        '''Whether to parse db or not.'''
        if dbase.lower().startswith('ensembl'): return 'ensembl' in self.list['DBParse']
        return dbase.lower() in self.list['DBParse']
#########################################################################################################################
    def dbLinks(self,dbase):    ### Whether to parse db links or not
        #!# Sort out case issues! #!# REDUDNANT METHOD #!#
        if not self.obj['Parent'] or not self.obj['Parent'].list['DBList'] or dbase in self.obj['Parent'].list['DBList']: return True
        else: return False
#########################################################################################################################
    def parseDB(self,dbase,details):  ### Extracts specific information to self.dict['DB'] (was specialDB())
        '''
        Extracts specific information to self.dict['DB'].
        >> dbase:str = Database identifier extracted from DR line of DAT file - '^(\S+);\s+(\S.+)$'[0]
        >> details:str = Database links extracted from DR line of DAT file - '^(\S+);\s+(\S.+)$'[1]
        << returns True is parsed DB is added or False if full (or no) details added.
        '''
        try:### ~ [0] Setup and short-cut parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if details.lower() in ['','none']: return False
            if self.parent() and self.parent().getBool('DBDetails'):
                if dbase not in self.dict['DB']: self.dict['DB'][dbase] = []
                self.dict['DB'][dbase].append(details)
                return False

            ### ~ [1] Each Database is dealt with individually. If the database is not here, nothing will happen ~~~~ ###
            ###~UniProt (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   UniProtKB/Swiss-Prot; O95793-1; STAU1_HUMAN; M.
            if dbase == "UniProtKB/Swiss-Prot":
                try: (uni_sv,uni_id) = rje.matchExp('^(\S+); (\S+);',details)
                except: return
                uni_acc = string.split(uni_sv,'-')[0]
                for db in ['UniAccNum','UniID','UniSV']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                self.dict['DB']['UniAccNum'].append(uni_acc)
                self.dict['DB']['UniID'].append(uni_id)
                if uni_acc.find('-') > 0: self.dict['DB']['UniSV'].append(uni_sv)
                return True

            ###~EnsEMBL (From IPI/UniProt DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   ENSEMBL; ENSP00000346163; ENSG00000124214; -.       = IPI
            ## DR   Ensembl; ENSG00000112530; Homo sapiens.             = UniProt
            ## Also EnsemblMetazoa; etc.
            if dbase.lower().startswith('ensembl'):
                for db in ['ENSG','ENSP','ENST']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                if rje.matchExp('(ENS\S*P\S+);',details): self.dict['DB']['ENSP'].append(rje.matchExp('(ENS\S*P\S+);',details)[0])
                if rje.matchExp('(ENS\S*T\S+);',details): self.dict['DB']['ENST'].append(rje.matchExp('(ENS\S*T\S+);',details)[0])
                if rje.matchExp('(ENS\S*G\S+);',details): self.dict['DB']['ENSG'].append(rje.matchExp('(ENS\S*G\S+);',details)[0])
                if rje.matchExp('(ENS\S*G\S+)\.',details): self.dict['DB']['ENSG'].append(rje.matchExp('(ENS\S*G\S+)\.',details)[0])
                if rje.matchExp('\s(\S+-P\S)',details):  # Ensembl subset sequences
                    for ens in string.split(details):
                        ensm = rje.matchExp('(\S+)-(P\S)',ens)
                        if ensm: self.dict['DB']['ENSP'].append('%s-%s' % ensm)
                        else:
                            ensm = rje.matchExp('(\S+)-(R\S)',ens)
                            if ensm: self.dict['DB']['ENST'].append('%s-%s' % ensm)
                        if ensm and ensm[0] not in self.dict['DB']['ENSG']: self.dict['DB']['ENSG'].append(ensm[0])
                # Cleanup empty entries for xref output
                for db in ['ENSG','ENSP','ENST']:
                    if not self.dict['DB'][db]: self.dict['DB'].pop(db)
                return True

            ###~RefSeq/NCBI (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   REFSEQ_REVIEWED; NP_059347; GI:82659087; -.
            ## DR   REFSEQ_PROVISIONAL; NP_004520; GI:4758720; -.
            ## DR   REFSEQ_VALIDATED; NP_001014313; GI:62122851; -.
            ## DR   RefSeq; XP_006499972.1; XM_006499909.1. [Q9CQV8-1]
            ## DR   RefSeq; NP_033562.3; NM_009536.4.
            if dbase.upper().startswith('REFSEQ'):
                if not self.dict['DB'].has_key('RefSeq'): self.dict['DB']['RefSeq'] = []
                self.dict['DB']['RefSeq'].append(rje.matchExp('(\S+_\S+);',details)[0])
                return True

            ###~HGNC (From IPI/UniProt DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   HGNC; 11370; STAU1; -.      = IPI
            ## DR   HGNC; HGNC:19152; PACRG.    = UniProt
            if dbase == 'HGNC':
                for db in ['Symbol','HGNC']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                if rje.matchExp('^HGNC:(\d+); (\S+).',details): (hgnc,symbol) = rje.matchExp('^HGNC:(\d+); (\S+).',details)
                else: (hgnc,symbol) = rje.matchExp('^(\d+); (\S+);',details)
                self.dict['DB']['HGNC'].append(hgnc)
                self.dict['DB']['Symbol'].append(symbol)
                return True

            ###~MGI (From UniProt DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   MGI; MGI:1891917; Ywhab.
            if dbase == 'MGI':
                for db in ['Symbol','MGI']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                if rje.matchExp('^MGI:(\d+); (\S+).',details): (hgnc,symbol) = rje.matchExp('^MGI:(\d+); (\S+).',details)
                else: (hgnc,symbol) = rje.matchExp('^(\d+); (\S+);',details)
                self.dict['DB']['MGI'].append(hgnc)
                self.dict['DB']['Symbol'].append(symbol)
                return True

            ###~Entrez Gene (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   Entrez Gene; 6780; STAU1; -.
            if dbase == 'Entrez Gene':
                if not self.dict['DB'].has_key('Entrez'): self.dict['DB']['Entrez'] = []
                self.dict['DB']['Entrez'].append(rje.matchExp('^(\d+);',details)[0])
                return True
                
            ###~FlyBaseGene (From UniProt)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   FlyBase; FBgn0010339; 128up.
            if dbase == 'FlyBase':
                if not self.dict['DB'].has_key('FlyBase'): self.dict['DB']['FlyBase'] = []
                self.dict['DB']['FlyBase'].append(rje.matchExp('^(FBgn\d+);',details)[0])
                return True

            ###~WormBaseGene (From UniProt)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   WormBase; F15E11.1; CE16999; WBGene00017490; pud-2.1.
            if dbase == 'WormBase':
                if not self.dict['DB'].has_key('WormBase'): self.dict['DB']['WormBase'] = []
                try: self.dict['DB']['WormBase'].append(rje.matchExp('(WBGene\d+);',details)[0])
                except: self.errorLog('Cannot parse WBGene from "%s"!' % details)
                return True

            ###~VectorBaseGene (From UniProt)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   VectorBase; AGAP001357; Anopheles gambiae.
            if dbase in ['VectorBase']:
                if dbase not in self.dict['DB']: self.dict['DB'][dbase] = []
                self.dict['DB'][dbase].append(string.split(details,';')[0])
                return True

            ###~PFam (From UniProt)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   Pfam; PF07714; Pkinase_Tyr; 1.
            if dbase == 'Pfam':     ## Makes a dictionary of {PfamID:'name; number'}
                if 'Pfam' not in self.dict['DB']: self.dict['DB']['Pfam'] = {}
                details = string.split(details,'; ')
                self.dict['DB']['Pfam'][details[0]] = string.join(details[1:],'; ')
                return True

            ### ~ [2] ~ No Special DB recognised: add full details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if dbase not in self.dict['DB']: self.dict['DB'][dbase] = []
            self.dict['DB'][dbase].append(details)
            return False
        except: self.errorLog('Problem with UniProtEntry.specialDB(%s)' % dbase); return False
#########################################################################################################################
    def cc2ft(self):    ### Adds extra full-length features based on LOCATION and TISSUE
        '''Adds extra full-length features based on LOCATION and TISSUE.'''
        try:
            ### Setup ###
            pos = (1,self.obj['Sequence'].aaLen())
            go = []
            if self.dict['DBLinks'].has_key('GO'):
                go = self.dict['DBLinks']['GO'][0:]
                for g in go[0:]:
                    if g.find('; C:') < 0:
                        go.remove(g)
            cc = 'SUBCELLULAR LOCATION'
            if self.dict['Comments'].has_key(cc):
                go = string.split(string.join(self.dict['Comments'][cc],'; '),'; ') + go
            ftlist = {'TISSUE':self.list['Tissues'],'LOCATION':go}
            ### Add FT ##
            for type in ftlist.keys():
                for desc in ftlist[type]:
                    if {'Type':type,'Desc':desc,'Start':pos[0],'End':pos[1]} not in self.list['Feature']:
                        self.list['Feature'].append({'Type':type,'Desc':desc,'Start':pos[0],'End':pos[1]})
            return True                
        except:
            self.log.errorLog('UniProtEntry.cc2ft() has gone wrong!')
            return False
#########################################################################################################################
    def maskFT(self,types=['EM'],inverse=False,mask='X',log=True):   ### Masks given feature types
        '''
        Masks given feature types.
        >> types:list of str [['EM']] = types of feature to mask
        >> inverse:bool [False] = whether to mask all sequence *except* listed types
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to log the affects of masking
        '''
        try:### Setup ###
            types = rje.listUpper(types[0:])
            oldseq = self.obj['Sequence'].info['Sequence'][0:]
            newseq = mask * len(oldseq)
            mx = 0
            prex = oldseq.count(mask)
            if inverse:
                (newseq,oldseq) = (oldseq,newseq)
            ### Mask ###
            for ft in self.list['Feature']:
                if ft['Type'].upper() in types:
                    oldseq = oldseq[:ft['Start']-1] + newseq[ft['Start']-1:ft['End']] + oldseq[ft['End']:]
                    mx += 1
            ### Update ###
            self.obj['Sequence'].info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if log:
                mtxt = self.info['Name']
                if inverse:
                    mtxt += ' inverse'
                self.printLog('#MASK','%s masked %d features. (%d %s added.)' % (mtxt,mx,maskx,mask))
        except: self.errorLog('Problem masking %s features from %s' % (types,self.shortName()))
#########################################################################################################################
    def caseFT(self,types=[]):   ### Masks given feature types
        '''
        Masks given feature types.
        >> types:list of str [] = types of feature to be upper case
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sequence = self.obj['Sequence'].info['Sequence'][0:].lower()
            ### ~ [2] ~ Change Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ft in self.list['Feature']:
                if ft['Type'] in types:
                    sequence = sequence[:ft['Start']-1] + sequence[ft['Start']-1:ft['End']].upper() + sequence[ft['End']:]
            self.obj['Sequence'].dict['Case'] = rje_sequence.caseDict(sequence)
        except: self.errorLog('Problem masking %s features from %s' % (types,self.shortName()))
#########################################################################################################################
    def orderFT(self):  ### Orders features by start, end, type
        '''Orders features by start, end, type.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newft = []
            ### ~ [2] ~ Order ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ft in self.list['Feature']:
                i = 0
                while i < len(newft):
                    if ft['Start'] < newft[i]['Start']: break
                    elif ft['Start'] > newft[i]['Start']: i += 1; continue
                    if ft['End'] < newft[i]['End']: break
                    elif ft['End'] > newft[i]['End']: i += 1; continue
                    if ft['Type'] < newft[i]['Type']: break
                    elif ft['Type'] > newft[i]['Type']: i += 1; continue
                    i += 1
                newft.insert(i,ft)
            self.list['Feature'] = newft
        except: self.errorLog('Problem ordering %s features' % self.info['Name'])            
#########################################################################################################################
    ### <3> ### UniProt Conversion and Saving                                                                           #
#########################################################################################################################
    def uniProtFromSeq(self,seq=None,sequence='',name='',data={},ft=[]): ### Converts into UniProtEntry object (self!)
        '''
        Converts into UniProtEntry object (self!).
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns self if successful or None if fails
        '''
        try:
            ### Setup ###
            if not seq and (sequence and name):
                seq = rje_sequence.Sequence(log=self.log)
                seq.info['Name'] = name
                seq.info['Sequence'] = sequence.upper()  #!# Change at some point to allow mixed case!
                seq.info['Type'] = 'Protein'
                seq.extractDetails()    #!#gnspacc=self.opt['GeneSpAcc'])
            if not seq:
                seq = self.obj['Sequence']
            if not seq:
                raise ValueError('No sequence information given')
            self.obj['Sequence'] = seq

            ### Update self ###
            for key in data:
                self.dict['Data'][key] = data[key]
            self.list['Feature'] += ft
            if seq.info['DBase'] != 'trembl':
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['ID'],seq.info['Type'],seq.aaLen())]
            elif seq.info['SpecCode'] not in ['None','UNK']:
                self.dict['Data']['ID'] = ['%s_%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['SpecCode'],seq.info['Type'],seq.aaLen())]
            else:
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['Type'],seq.aaLen())]
            if self.dict['Data'].has_key('AC'):
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']] + self.dict['Data']['AC']
            else:
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']]
            if seq.info['Description'].lower() != 'none':
                self.dict['Data']['DE'] = [seq.info['Description']]
            else:
                self.dict['Data']['DE'] = ['']
            if seq.info['Species'] not in ['None','Unknown',seq.info['SpecCode'],'']:
                if self.dict['Data'].has_key('OS'):
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']] + self.dict['Data']['OS']
                else:
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']]
            dt = string.split(time.ctime())
            if self.dict['Data'].has_key('DT'):
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])] + self.dict['Data']['DT']
            else:
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])]
            if self.dict['Data'].has_key('CC'):
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()] + self.dict['Data']['CC']
            else:
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()]

            #X#self.deBug(self.dict['Data'])
            if self.process(logft=False,cleardata=False):
                return self
            return None
        except:
            self.log.errorLog('UniProtEntry.uniProtFromSeq() has gone wrong.',quitchoice=True)
            return None
#########################################################################################################################
## End of SECTION III : UniProtEntry Class                                                                              #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: General UniProt Methods                                                                                 #
#########################################################################################################################
def downloadUniProt(callobj):   ### Downloads the UniProt database using the attributes of callobj
    '''Downloads the UniProt database using the attributes of callobj.'''
    try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        outdir = callobj.info['UniPath']
        mydir = os.path.abspath(os.curdir)
        if not os.path.exists(outdir): rje.mkDir(callobj,outdir)
        ## ~ [1a] Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not callobj.force():
            sfile = outdir + 'uniprot_sprot.dat'
            tfile = outdir + 'uniprot_trembl.dat'
            sfound = os.path.exists(sfile) or os.path.exists(sfile+'.gz')
            tfound = os.path.exists(tfile) or os.path.exists(tfile+'.gz')
            if sfound and tfound and (callobj.stat['Interactive'] < 1 or rje.yesNo('Files found. Skip UniProt download?')):
                callobj.printLog('#DBASE','UniProt DAT files found and force=F - no download.')
                return True
        ## ~ [1b] Prepare for download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ftproot = 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/'
        dbtxt = 'UniProt from %s -> %s' % (ftproot,outdir)
        if callobj.stat['Interactive'] >= 1 and not rje.yesNo('Download %s?' % dbtxt): return False
        elements = ['knowledgebase/complete/*dat.gz','README','relnotes.txt']

        ### ~ [2] Download files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        callobj.log.printLog('#DBASE','Downloading %s: %d elements.' % (dbtxt,len(elements)))
        for element in elements:
            callobj.log.printLog('#DBASE','Downloading %s' % (element))
            os.chdir(outdir)
            os.system('wget %s' % rje.makePath(ftproot+element,wholepath=True))
            os.chdir(mydir)
        return True
    except: callobj.errorLog('Major error during UniProt download')
#########################################################################################################################
def processUniProt(callobj,makeindex=True,makespec=True,makefas=True,temp=False):     ### Processes UniProt making index file and spectable as appropriate
    '''Processes UniProt making index file and spectable as appropriate.'''
    try:### ~ [1] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not makeindex and not makespec and not makefas:
            return callobj.log.printLog('#DB','No call to make UniProt index/spectable/fasta files.')
        unipath = callobj.info['UniPath']
        try: indexfile = unipath + callobj.info['DBIndex']
        except: indexfile = unipath + 'uniprot.index'
        if callobj.getStrLC('Name'): datfiles = [callobj.getStr('Name')]
        else:
            datfiles = glob.glob('%s*.dat' % unipath)
            datfiles.sort()
        if not datfiles: return callobj.errorLog('No *.dat files found in %s!' % unipath,printerror=False)
        ## ~ [1a] Check index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makeindex or callobj.stat['Interactive'] < 1 or rje.yesNo('Check/Make UniProt index file?'):
            if not os.path.exists(indexfile) or callobj.opt['Force']: makeindex = True
            else:
                makeindex = False
                for dat in datfiles:
                    if rje.isYounger(indexfile,dat) != indexfile:
                        makeindex = True
                        break
        ## ~ [1b] SpecTable Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        # - can then use "grep 'Metazoa' uniprot.spec.tdt | gawk '{print $1 }' | sort -u" to extract, e.g. Metazoan species codes
        spectable = unipath + 'uniprot.spec.tdt'
        spectemp = unipath + 'uniprot.spec.txt'
        if makespec or (callobj.stat['Interactive'] > 0 and rje.yesNo('Check/Make UniProt species file?')):
            if callobj.opt['Force'] or not os.path.exists(spectable): makespec = True
            else:
                makespec = False
                for dat in datfiles:
                    if rje.isYounger(spectable,dat) != spectable:
                        makespec = True
                        break
        ## ~ [1c] Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makefas or (callobj.stat['Interactive'] > 0 and rje.yesNo('Reformat UniProt to fasta files?')):
            if callobj.opt['Force']: makefas = True
            else:
                makefas = False
                for dat in datfiles:
                    fas = rje.baseFile(dat) + '_all.fas'
                    if not os.path.exists(fas) or rje.isYounger(fas,dat) != fas:
                        makefas = True
                        break
        ## ~ [1d] Abort if no need to process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not makeindex and not makespec and not makefas:
            return callobj.log.printLog('#DB','No need to make UniProt index/spectable/fasta files. (Force=F)')
        callobj.deBug('Index: %s; Spec: %s; Fas: %s' % (makeindex,makespec,makefas))

        ### ~ [2] Setup Output and Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ix = -1 # No. index lines
        tx = 0  # No. taxa lines
        sx = 0  # No. seqs
        dx = 0  # No. dat files
        callobj.deBug('%d forks: %s' % (callobj.getStat('Forks'),not callobj.getBool('NoForks')))
        temp = False
        try:
            if callobj.getStat('Forks') > 1 and not callobj.getBool('NoForks'): forkProcess(callobj,datfiles,makeindex,makespec,makefas)
            else:
                callobj.warnLog('Uniprot Index generation currently does not support NoForks. Will use forks=1.')
                callobj.setStat({'Forks': 1})
                forkProcess(callobj,datfiles,makeindex,makespec,makefas)
        except: pass
        if temp:
            callobj.warnLog('Using old Uniprot index generation. Not sure why you want to do that!',dev=True)
            if makeindex:
                INDEX = open('%s.head' % indexfile,'w')
                for d in range(len(datfiles)):
                    dat = datfiles[d]
                    INDEX.write('#%d=%s\n' % (d+1,os.path.basename(dat)))
                INDEX.close()
                INDEX = open('%s.temp' % indexfile,'w')     # Clear it
            if makespec: SPEC = open(spectemp,'w')          # Clear it
            if makefas:
                (append,callobj.opt['Append']) = (callobj.opt['Append'],False)
                for dat in datfiles: rje.backup(callobj,rje.baseFile(dat) + '_all.fas')
                callobj.opt['Append'] = append

            ### ~ [3] Process Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            wanted = ['ID','AC','//']
            if makespec: wanted += ['OS','OC','OX']
            if makefas: wanted += ['DE','  ']
            for dat in datfiles:        #!# Add Forking to accelerate process #!#
                ## ~ [3a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                grepindex = False
                fasfile = rje.baseFile(dat) + '_all.fas'
                if callobj and 'GrepDat' in callobj.opt and callobj.opt['GrepDat']:
                    grepindex = True
                    wanted.remove('//'); wanted.append('//')
                    wantgrep = string.join(wanted,' |')
                    callobj.printLog('#GREP',"grep '^(%s)' -E %s" % (wantgrep,dat))
                    DAT = os.popen("grep '^(%s)' -E %s" % (wantgrep,dat))
                else: DAT = open(dat, 'r')
                ix = seq_index = -1
                dx += 1
                dbtext = 'Processing %d of %d files (%s):' % (dx,len(datfiles),dat)
                id = acc = code = species = seq = taxid = desc = ''
                taxonomy = []
                file_pos = 0; lx = 0; sx = 0
                ## ~ [3b] Read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while 'Reading File':
                    line = DAT.readline()
                    if not line:
                        DAT.close()
                        break
                    ix += 1
                    if line[:2] not in wanted: continue
                    ## ~ [3c] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if line[0:2] == 'ID':
                        id = string.split(line)[1]
                        try: code = string.split(id,'_')[1]
                        except: pass
                        seq_index = file_pos
                    elif line[0:2] == 'AC':
                        extras = rje.matchExp('^AC\s+(\S.+)\s*$',line)[0]
                        if acc: acc = '%s %s' % (acc,extras)
                        else: acc = extras
                    elif line[:2] == 'DE':
                        if desc: desc += rje.chomp(line[4:])
                        else: desc = rje.chomp(line[5:])
                    elif line[0:2] == 'OS' and not species:
                        try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)\s+\(+',line)[0]
                        except:
                            try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)[\.,]',line)[0]
                            except:
                                try: species = rje.matchExp('^OS\s+(\S.+irus)',line)[0]
                                except: species = rje.matchExp('^OS\s+(\S.+\S)',line)[0]
                        if not code: code = rje_sequence.getSpecCode(species)
                    elif line[0:2] == 'OC':
                        taxonomy = taxonomy + string.split(re.sub('\s+','',rje.matchExp('^OC\s+(\S.+)$',line)[0]),';')
                        if '' in taxonomy: taxonomy.remove('')
                    elif line[0:2] == 'OX':
                        taxid = rje.matchExp('^OX\s+NCBI_TaxID=(\d+)',line)[0]
                        taxonomy.append(species)
                    elif line[:2] == '  ': seq += string.replace(rje.chomp(line[5:]),' ','')
                    elif line[0:2] == '//':
                        sx += 1
                        if makefas and seq: open(fasfile,'a').write('>%s__%s %s\n%s\n' % (id,string.split(acc,';')[0],desc,seq))
                        if code and makespec:
                            tx += 1
                            SPEC.write('%s\t:%s:\t:%s:\n' % (code, string.join(taxonomy,':'), taxid))
                        if makeindex and not grepindex:
                            file_pos = DAT.tell()
                            ilist = string.split(string.join(string.split(acc) + [id],''),';')
                            while '' in ilist: ilist.remove('')
                            for acc in ilist: INDEX.write('%s;%d:%d\n' % (acc,dx,seq_index))
                        id = acc = code = species = seq = taxid = desc = ''
                        taxonomy = []
                        seq_index = -1
                        lx += 1
                        if lx == 1000: callobj.progLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx))); lx = 0
                callobj.printLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx)))
                if grepindex and makeindex:
                    DAT = open(dat, 'r'); sx = 0
                    ix = seq_index = -1
                    dbtext = 'Indexing %d of %d files (%s):' % (dx,len(datfiles),dat)
                    id = acc
                    taxonomy = []
                    file_pos = 0; lx = 0
                    ## ~ [3b] Read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    while 'Reading File':
                        line = DAT.readline()
                        if not line:
                            DAT.close()
                            break
                        ix += 1
                        ## ~ [3c] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        if line[0:2] == 'ID':
                            seq_index = file_pos
                            id = string.split(line)[1]
                            INDEX.write('%s;%d:%d\n' % (id,dx,seq_index))
                        elif line[0:2] == 'AC':
                            for acc in string.split(rje.matchExp('^AC\s+(\S.+)\s*$',line)[0]): INDEX.write('%s;%d:%d\n' % (string.replace(acc,';',''),dx,seq_index))
                        elif line[0:2] == '//':
                            sx += 1; file_pos = DAT.tell(); seq_index = -1; lx += 1
                            if lx == 1000: callobj.progLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx))); lx = 0
                    callobj.printLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx)))                
        
            ### ~ [4] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dbtext = 'Processed %d *.dat files:' % len(datfiles)
            callobj.log.printLog('#DB','%s %s lines; %s entries.' % (dbtext,rje.integerString(ix),rje.integerString(sx)))
        ### ~ [4a] Sort out index files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if makeindex:
                try: INDEX.close()
                except: pass
               
                ### Subsort within subsets ##
                subset = 3  # Cannot use 2 even though I want to because of fucking dumb UNIX sort P_ problem. Grrrr!
                TMP = open('%s.temp' % indexfile,'r')
                if callobj and callobj.opt['GrepDat']:
                    ## Read first sort and split into subsorts ##
                    subindex = {}   # Dictionary of subindex lines
                    reading = True
                    px = 0
                    while reading:
                        line = TMP.readline()
                        if line in ['',None]: break
                        px += 1
                        start = line[:subset]
                        if start not in subindex: subindex[start] = []
                        subindex[start].append(string.replace(line,';;',';'))
                        callobj.progLog('\r#INDEX','Processing index lines: %s' % rje.integerString(px))
                    callobj.printLog('\r#INDEX','Processed %s index lines.' % rje.integerString(px),log=False,newline=False)
                    TMP.close()
                    ### Process SubFiles ###
                    tmplines = callobj.loadFromFile('%s.head' % indexfile)
                    INDEX = open(indexfile,'w')
                    INDEX.write(string.join(tmplines,''))
                    sx = 0.0; stot = len(subindex)
                    for start in rje.sortKeys(subindex):
                        subindex[start].sort()
                        INDEX.write(string.join(subindex[start],''))
                        callobj.progLog('\r#INDEX','Saving index: %.2f%%' % (sx/stot)); sx += 100.0
                    INDEX.close()
                    callobj.printLog('\r#INDEX','Saving index %s complete.' % indexfile)
                else:
                    subindex = []   # List of temp subindex files
                    ## Read first sort and split into subsorts ##
                    tmplines = []
                    reading = True
                    px = 0
                    while reading:
                        line = TMP.readline()
                        if line in ['',None]: break
                        px += 1
                        start = 'index_%s.tmp' % line[:subset]
                        if start not in subindex: subindex.append(start)
                        open(start,'a').write(line)
                        callobj.log.printLog('\r#INDEX','Processing index lines: %s' % rje.integerString(px),log=False,newline=False)
                    callobj.log.printLog('\r#INDEX','Processed %s index lines into sub-files' % rje.integerString(px),log=False,newline=False)
                    TMP.close()

                    ### Process SubFiles ###
                    subindex.sort()
                    tmplines = callobj.loadFromFile('%s.head' % indexfile)
                    INDEX = open(indexfile,'w')
                    INDEX.write(string.join(tmplines,''))
                    for file in subindex:
                        TMP = open(file,'r')
                        tmplines = []
                        reading = True
                        while reading:
                            line = TMP.readline()
                            if line in ['',None]: break
                            tmplines.append(line)
                        indexdict = {}
                        pxx = len(tmplines)
                        px = 0
                        while tmplines:
                            px += 50.0
                            indata = string.split(rje.chomp(tmplines.pop(0)),';')
                            if indexdict.has_key(indata[0]): indexdict[indata[0]] 
                            else: indexdict[indata[0]] = []
                            for idat in indata[1:]:
                                if idat not in indexdict[indata[0]]: indexdict[indata[0]].append(idat)
                            callobj.log.printLog('\r#INDEX','Processing %s index lines: %.2f%%' % (file,(px/pxx)),log=False,newline=False)
                        for key in rje.sortKeys(indexdict):
                            px += 50.0
                            outdata = [key,string.join(indexdict.pop(key),';')]
                            INDEX.write('%s\n' % string.join(outdata,';'))
                            callobj.log.printLog('\r#INDEX','Processing %s index lines: %.2f%%' % (file,(px/pxx)),log=False,newline=False)
                        callobj.log.printLog('\r#INDEX','Processing %s index lines: Complete.' % file,log=False)
                        TMP.close()
                        os.unlink(file)
                    INDEX.close()
                    os.unlink('%s.temp' % indexfile)
                    os.unlink('%s.head' % indexfile)
                callobj.log.printLog('\r#INDEX','Generation of UniProt index complete.')
                            
                ### Old way using UNIX sort, which is quicker but doesn't work!! ###
                if 'you_want_to_do_it_the_quick_but_dangerous_way_with_UNIX_sort' in callobj.cmd_list:               
                    if callobj.opt['Win32']: callobj.log.printLog('#DB','Cannot cleanup %s.temp with UNIX sort! Talk to Rich!' % indexfile)
                    else:
                        callobj.log.printLog('#DB','Cleanup of %s.temp with UNIX sort ...' % indexfile,log=False)
                        if not os.system('sort %s.temp > %s.sort' % (indexfile,indexfile)):
                            os.unlink('%s.temp' % indexfile)
                            callobj.log.printLog('#DB','%s.temp sorted with UNIX sort. Concatenating...' % indexfile,log=False)
                            if not os.system('cat %s.head %s.sort > %s' % (indexfile,indexfile,indexfile)):
                                os.unlink('%s.sort' % indexfile)
                                os.unlink('%s.head' % indexfile)
                                callobj.log.printLog('#DB','Sort and Concatenation complete. Temporary files deleted. %s created.' % indexfile)
                            else: callobj.log.printLog('#ERR','Concatenation failed! Temporary files not deleted. %s not created.' % indexfile)                           
                        else: callobj.log.printLog('#ERR','Sorting of %s.temp with UNIX sort Failed!' % indexfile)

        ### ~ [4a] Sort out spectable files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makespec:
            specgrep = 'sort -u %s | grep -v \'^9\' > %s' % (spectemp,spectable)
            specgrep2 = 'sort -u %s > %suniprot.9spec.tdt' % (spectemp,unipath)
            try: SPEC.close()
            except: pass
            if callobj.opt['Win32']:
                callobj.log.printLog('#DB','Cannot cleanup %s with UNIX sort. Very large file!' % spectemp)
                callobj.log.printLog('#DB','Use: "%s"' % (specgrep))
            else:
                callobj.progLog('#DB','Cleanup with UNIX: %s ...' % specgrep)
                if not os.system(specgrep):
                    callobj.printLog('#DB','Cleanup with UNIX: %s - complete' % specgrep)
                    os.system(specgrep2)
                    if callobj.stat['Interactive'] < 1 or rje.yesNo('Delete %s?' % spectemp): os.unlink(spectemp)
                else:
                    callobj.log.printLog('#DB','Cleanup of %s with UNIX sort failed. Very large file!' % spectemp)
                    callobj.log.printLog('#DB','Tried: "%s"' % specgrep)
        
    except: callobj.log.errorLog('Error in processUniProt()',printerror=True)
#########################################################################################################################
def forkProcess(callobj,alldatfiles,makeindex,makespec,makefas,forkbytes=1e8):   ### Forks out DAT download processing to speed up (hopefully!)
    '''Forks out DAT download processing to speed up (hopefully!).'''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self = callobj
        unipath = callobj.getStr('UniPath')
        forkx = self.getStat('Forks')          # Number of forks to have running at one time
        forks = {}; forktot = 0             # List of active forks {pids:id}
        killforks = self.getStat('KillForks')  # Time in seconds to wait after main thread has apparently finished
        killtime = time.time()
        try: indexfile = unipath + callobj.getStr('DBIndex')
        except: indexfile = unipath + 'uniprot.index'
        spectable = unipath + 'uniprot.spec.tdt'
        spectemp = unipath + 'uniprot.spec.txt'
        datfiles = alldatfiles[0:]             # Work through each of these in turn
        datpos = long(0)                    # File position
        DAT = open(datfiles[0],'r')
        DAT.seek(0,2); datend = DAT.tell(); DAT.seek(0)
        if forkx == 1: forkbytes = datend
        totline = 0; totseq = 0
        ## ~ [0a] ~ Setup general files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makeindex:
            INDEX = open('%s.head' % indexfile,'w')
            for d in range(len(datfiles)):
                dat = datfiles[d]
                INDEX.write('#%d=%s\n' % (d+1,os.path.basename(dat)))
            INDEX.close()
            open('%s.temp' % indexfile,'w')     # Clear it
        if makespec: open(spectemp,'w')          # Clear it
        if makefas:
            (append,callobj.opt['Append']) = (callobj.opt['Append'],False)
            for dat in datfiles: rje.backup(callobj,rje.baseFile(dat) + '_all.fas')
            callobj.opt['Append'] = append
        callobj.deBug('Index: %s; Spec: %s; Fas: %s' % (makeindex,makespec,makefas))
        ### ~ [1] ~ Fork out processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while datfiles or forks:
            ### ~ [1a] ~ Start a new fork if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while datfiles and (len(forks) < forkx):
                dat = datfiles[0]
                file_pos = datpos
                datpos += long(forkbytes)
                if datpos >= datend:
                    self.printLog('\r#DAT','',log=False)
                    self.printLog('#DAT','End of %s reached.' % datfiles[0])
                    endpos = datend; datpos = 0
                    DAT.close()
                    datfiles.pop(0)
                    if datfiles:
                        DAT = open(datfiles[0],'r')
                        DAT.seek(0,2); datend = DAT.tell(); DAT.seek(0)
                        if forkx == 1: forkbytes = datend
                else: endpos = datpos
                ## ~ [1b] ~ Create fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newpid = os.fork()
                if newpid == 0: # child
                    ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORK CONTENTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ###
                    if makeindex: INDEX = open('%s.%d.temp' % (indexfile,forktot),'w')      # Clear it
                    if makespec: SPEC = open('%s.%d.txt' % (spectemp[:-4],forktot),'w')      # Clear it
                    if makefas: FAS = open('%s.%d.fas' % (rje.baseFile(dat),forktot),'w')
                    ### ~ [3] Process Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                    wanted = ['ID','AC','//']
                    if makespec: wanted += ['OS','OC','OX']
                    if makefas: wanted += ['DE','  ']
                    ## ~ [3a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    FDAT = open(dat, 'r'); FDAT.seek(file_pos)
                    while FDAT.readline()[:2] != 'ID': file_pos = FDAT.tell()     # Progress to next entry
                    FDAT.seek(file_pos)
                    id = acc = code = species = seq = taxid = desc = lastspec = thisspec = ''
                    taxonomy = []
                    lx = 0; sx = 0; seq_index = -1; dx = alldatfiles.index(dat) + 1
                    ## ~ [3b] Read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    while 'Reading File':
                        line = FDAT.readline()
                        if not line:
                            FDAT.close()
                            break
                        lx += 1
                        if line[:2] not in wanted: continue
                        ## ~ [3c] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        try:
                            if line[0:2] == 'ID' and line[:3] != 'ID=':
                                id = string.split(line)[1]
                                try: code = string.split(id,'_')[1]
                                except: pass
                                seq_index = file_pos
                            elif line[0:2] == 'AC':
                                extras = rje.matchExp('^AC\s+(\S.+)\s*$',line)[0]
                                if acc: acc = '%s %s' % (acc,extras)
                                else: acc = extras
                            elif line[:2] == 'DE':
                                if desc: desc += rje.chomp(line[4:])
                                else: desc = rje.chomp(line[5:])
                            elif line[0:2] == 'OS' and not species:
                                try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)\s+\(+',line)[0]
                                except:
                                    try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)[\.,]',line)[0]
                                    except:
                                        try: species = rje.matchExp('^OS\s+(\S.+irus)',line)[0]
                                        except: species = rje.matchExp('^OS\s+(\S.+\S)',line)[0]
                                if not code: code = rje_sequence.getSpecCode(species)
                            elif line[0:2] == 'OC':
                                taxonomy = taxonomy + string.split(re.sub('\s+','',rje.matchExp('^OC\s+(\S.+)$',line)[0]),';')
                                if '' in taxonomy: taxonomy.remove('')
                            elif line[0:2] == 'OX':
                                taxid = rje.matchExp('^OX\s+NCBI_TaxID=(\d+)',line)[0]
                                taxonomy.append(species)
                            elif line[:2] == '  ': seq += string.replace(rje.chomp(line[5:]),' ','')
                            elif line[0:2] == '//':
                                sx += 1
                                if makefas and seq and id: FAS.write('>%s__%s %s\n%s\n' % (id,string.split(acc,';')[0],desc,seq))
                                elif makefas and seq and not id: self.errorLog('Protein ID missing!',printerror=False)
                                if code and makespec:
                                    thisspec = '%s\t:%s:\t:%s:\n' % (code, string.join(taxonomy,':'), taxid)
                                    if thisspec != lastspec: SPEC.write(thisspec)
                                else: thisspec = ''
                                file_pos = FDAT.tell()
                                if makeindex and id:
                                    ilist = string.split(string.join(string.split(acc) + [id],''),';')
                                    while '' in ilist: ilist.remove('')
                                    for acc in ilist: INDEX.write('%s;%d:%d\n' % (acc,dx,seq_index))
                                id = acc = code = species = seq = taxid = desc = ''
                                (lastspec,thisspec) = (thisspec,'')
                                taxonomy = []
                                seq_index = -1
                                if file_pos >= endpos: break
                        except: self.errorLog('Fork %s problem with line: "%s"' % (forktot, rje.chomp(line)))
                    open('%s.%d.prog' % (rje.baseFile(dat),forktot),'w').write('%d %d' % (lx,sx))
                    try:
                        if makeindex: INDEX.close()
                        if makespec: SPEC.close()
                        if makefas: FAS.close()
                    except:
                        self.errorLog('Fork %d had trouble closing file handles' % forktot)
                        fbool = [makeindex,makespec,makefas]
                        fcheck = ['%s.%d.temp' % (indexfile,forktot),'%s.%d.txt' % (spectemp[:-4],forktot),'%s.%d.fas' % (rje.baseFile(dat),forktot),'%s.%d.prog' % (rje.baseFile(dat),forktot)]
                        for f in range(3):
                            if fbool[f] and not rje.checkForFile(fcheck[f]): self.errorLog('%s not found!' % fcheck[f],printerror=False,quitchoice=True)
                    ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORK CONTENTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ###
                    os._exit(0)    # Exit process 
                elif newpid == -1: self.errorLog('Problem forking %s.' % id,printerror=False)  
                else:
                    forks[newpid] = (dat,forktot); forktot += 1
                    self.printLog('#PID','Fork %d - PID %s: %s -> %s' % (forktot,newpid,file_pos,endpos))
                #self.deBug('Forked PID %s (%s,%d)' % (newpid,dat,forktot-1))

            ### ~ [2] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('\r#DAT','%d active forks; %d finished; %s lines processed; %s sequences >>>' % (len(forks),forktot - len(forks),rje.iStr(totline),rje.iStr(totseq)))
            time.sleep(10)       # Sleep for 1s 
            forklist = self._activeForks(forks.keys())
            if len(forklist) != len(forks):
                self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                for pid in rje.sortKeys(forks):    # Go through current forks
                    if pid not in forklist:
                        (dat,fork) = forks.pop(pid)
                        if dat == 'sort':
                            rje.fileTransfer(fromfile='%s.%d.sort' % (spectemp[:-4],fork),tofile=spectemp,deletefrom=True,append=True)
                            os.unlink('%s.%d.txt' % (spectemp[:-4],fork))
                        else:                            
                            killtime = time.time()  # Reset killtime - still doing stuff
                            fcheck = ['%s.%d.temp' % (indexfile,fork),'%s.%d.txt' % (spectemp[:-4],fork),'%s.%d.fas' % (rje.baseFile(dat),fork),'%s.%d.prog' % (rje.baseFile(dat),fork)]
                            fbool = [makeindex,makespec,makefas]
                            for f in range(3):
                                if fbool[f] and not rje.checkForFile(fcheck[f]): self.errorLog('%s not found!' % fcheck[f],printerror=False,quitchoice=True)
                            prog = string.split(rje.chomp(open('%s.%d.prog' % (rje.baseFile(dat),fork),'r').readline()))
                            totline += string.atoi(prog[0]); totseq += string.atoi(prog[1])
                            os.unlink('%s.%d.prog' % (rje.baseFile(dat),fork))
                            #print pid, dat, fork, prog,
                            if makeindex: rje.fileTransfer(fromfile='%s.%d.temp' % (indexfile,fork),tofile='%s.temp' % (indexfile),deletefrom=True,append=True)
                            if makefas: rje.fileTransfer(fromfile='%s.%d.fas' % (rje.baseFile(dat),fork),tofile='%s_all.fas' % rje.baseFile(dat),deletefrom=True,append=True)
                            #print
                            self.progLog('\r#DAT','%d active forks; %d finished; %s lines processed; %s sequences >>>' % (len(forks),forktot - len(forks),rje.iStr(totline),rje.iStr(totseq)))
                            if makespec:
                                #x#rje.fileTransfer(fromfile='%s.%d.txt' % (spectemp[:-4],fork),tofile=spectemp,deletefrom=True,append=True)
                                newpid = os.fork() 
                                if newpid == 0: # child
                                    ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORK CONTENTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ###
                                    os.system('sort -u %s.%d.txt > %s.%d.sort' % (spectemp[:-4],fork,spectemp[:-4],fork))
                                    ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FORK CONTENTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ###
                                    os._exit(0)    # Exit process 
                                elif newpid == -1: self.errorLog('Problem forking sort %s.' % fork,printerror=False)  
                                else: forks[newpid] = ('sort',fork)
                                #self.deBug('Forked sort PID %s (%d)' % (newpid,fork))
            ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if time.time() - killtime > killforks:
                self.printLog('#KILL','%d seconds of main thread inactivity. %d forks still active!' % (killforks,len(forks)))
                for fork in forks: self.printLog('#FORK','Fork %s, PID %d still Active!' % (forks[fork],fork))
                if self.stat['Interactive'] < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                elif rje.yesNo('Kill hanging forks?'):
                    for fork in forks:
                        self.printLog('#KILL','Killing Fork %s, PID %d.' % (forks[fork],fork))
                        os.system('kill %d' % fork)
                else: killtime = time.time()
        self.printLog('\r#DAT','%d active forks; %d finished; %s lines processed; %s sequences >>>' % (len(forks),forktot - len(forks),rje.iStr(totline),rje.iStr(totseq)))

        ### ~ [3] ~ Process Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if makeindex:
            ### Subsort within subsets ##
            subset = 3  # Cannot use 2 even though I want to because of fucking dumb UNIX sort P_ problem. Grrrr!
            TMP = open('%s.temp' % indexfile,'r')
            ## Read first sort and split into subsorts ##
            subindex = {}   # Dictionary of subindex lines
            reading = True
            px = 0
            while reading:
                line = TMP.readline()
                if line in ['',None]: break
                px += 1
                start = line[:subset]
                if start not in subindex: subindex[start] = []
                subindex[start].append(string.replace(line,';;',';'))
                callobj.progLog('\r#INDEX','Processing index lines: %s' % rje.integerString(px))
            callobj.printLog('\r#INDEX','Processed %s index lines.' % rje.integerString(px),log=False,newline=False)
            TMP.close()
            ### Process SubFiles ###
            tmplines = callobj.loadFromFile('%s.head' % indexfile)
            INDEX = open(indexfile,'w')
            INDEX.write(string.join(tmplines,''))
            sx = 0.0; stot = len(subindex)
            for start in rje.sortKeys(subindex):
                subindex[start].sort()
                INDEX.write(string.join(subindex[start],''))
                callobj.progLog('\r#INDEX','Saving index: %.2f%%' % (sx/stot)); sx += 100.0
            INDEX.close()
            callobj.printLog('\r#INDEX','Saving index %s complete.' % indexfile)
            callobj.printLog('\r#INDEX','Generation of UniProt index complete.')
    except SystemExit: os._exit(0)
    except: self.errorLog('%s.forkProcess error' % self)
#########################################################################################################################
## End of SECTION IV : Generic Module Methods                                                                           #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MAIN PROGRAM                                                                                             #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print('Unexpected error during program setup:', sys.exc_info()[0])
        return
        
    ### Rest of Functionality... ###
    try:        
        uniprot = UniProt(mainlog,cmd_list)
        uniprot.run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
    #mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################

