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
Module:       rje_dbase
Description:  RJE Module to Handle Database manipulations and generations
Version:      2.3
Last Edit:    09/01/13
Copyright (C) 2007 Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to control generic database manipulations routinely used for me to generate customised
    databases:
    1. Download commonly used databases, primarily UniProt, EnsEMBL, PFam and PPI databases
    2. Reformat and index crucial UniProt data from uniprot.dat and trembl.dat for ease of extraction. (UNIX platform)
    3. Generate taxa-specific databases from input databases
    4. Generate non-redundant species-specific databases using EnsEMBL gene locus informtation

    By default, database paths are relative. To peform an update it is advised that a new directory is created and
    RJE_DBASE run in this directory with the dbdownload=LIST dbprocess=LIST and taxadb=FILE options. Once download and
    formatting is complete, the new files can be copied over the old files.

    Database download is controlled in two ways. UniProt and EnsEMBL are managed by their own respective modules. Other
    databases are currently read from a file, which is in (an attempt of) XML format of the basic form:
    <dbxml>
      <database name="EnsEMBL" ftproot="ftp://ftp.ensembl.org/pub/" outdir="EnsEMBL/Current-release">
        <file path="current_aedes_aegypti/data/fasta/pep/*.gz">Yellow Fever Mosquito</file>
      </database>
    </dbxml>

    The pre-version 1.2 options for making IPI-centred datasets can still be called using the makedb=FILE option along
    with its associated options: screenipi=T screenens=F ensloci=F.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Primary Database download and processing options ###
    dbdownload=LIST : List of EnsEMBL/UniProt/XML files containing details of databases to download []
    dbprocess=LIST  : List of EnsEMBL/UniProt/IPI database types to process []
    datindex=T/F    : Create an index file for the Uniprot DAT files in unipath if UniProt in dbprocess [True]
    spectable=T/F   : Makes a table of species codes, taxonomy and taxon_id from DAT files if dbprocess UniProt [True]
    taxadb=FILE     : File containing details of taxanomic sub-databases to make [None]
    formatdb=T/F    : Whether to BLAST format database after making [True]
    force=T/F       : Whether to force regeneration of existing files [False]
    ignoredate=T/F  : Whether to ignore the relative timestamps of files when assessing whether to regenerate [False]
    ensloci=T/F     : Reduce EnsEMBL to a single protein per locus, mapping UniProt where possible [True]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Database Path Details ###
    unipath=PATH    : Path to UniProt files [UniProt/]
    ipipath=PATH    : Path to IPI files [IPI/]
    enspath=PATH    : Path to EnsEMBL file [EnsEMBL/]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Database and sub-database manufacture ###
    dbformat=LIST   : Reformats UniProt, IPI or EnsEMBL databases using RJE_SEQ []
    makedb=FILE     : Makes a database from combined databases [None]
                        - Note that rje_seq commandline options will be applied to this database with the addition of a
                        goodspec=X filter applied from the taxalist=LIST
    useX=T/F        : Whether to use certain aspects of databases,
                      where X is:   uniprot/sprot/trembl/ensembl/known/novel/abinitio/ipi [All but ipi True]
    taxalist=LIST   : List of taxanomic groups to extract spec_codes for reduced database (else all) [None]
    speconly=T/F    : Will simply output a list of SPECIES codes to the makedb file, rather than making dbase [False]
    inversedb=T/F   : TaxaList is a list of taxanomic groups *NOT* to be in database [False]
    screenipi=T/F   : Species represented by IPI databases will be screened out of UniProt and EnsEMBL. [False]
    screenens=T/F   : Species represented by EnsEMBL will be screened out of UniProt [True]
    seqfilter=T/F   : Use rje_seq to filter sequences (True) or simply filter on Species Codes (False) [False]
    ensfilter=T/F   : Run EnsEMBL genomes through rje_seq to apply filters, rather than just concatenating [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_ensembl, rje_seq, rje_seqlist, rje_uniprot, rje_xml, rje_zen
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Initial working version for use with rje_uniprot V2.0
    # 1.1 - Added automated downloading of databases from file
    # 1.2 - Incorporated RJE_ENSEMBL for handling EnsEMBL genomes and making a one-protein-per-gene proteome.
    # 1.3 - Tidied code a little and improved comments/docstrings.
    # 2.0 - Heavily reorganised and modified module.
    # 2.1 - Added seqfilter=T/F to speed up TaxaDB manufacture.
    # 2.2 - Added use of rje_seqlist for TaxaDB manufacture.
    # 2.3 - Added construction of EnsEMBL TaxaDB sets during TaxaDB construction.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Look to BioPython.
    # [ ] : Use "grep 'Metazoa' uniprot.spec.tdt | gawk '{print $1 }' | head" to extract, e.g. Metazoan species codes
    # [ ] : subdb => Makes a reduced uniprot database of reduced details (not so much need with index?)
    # [ ] : Option to extract AccNum etc from other databases and make super-index. Would be able to use with dbindex=FILE.
    # [Y] : Add way to download EnsEMBL either without XML file or within XML file instructions.
    # [ ] : Add generation of IPI Links table from ipi.HUMAN.xrefs and EnsLoci.
    # [ ] : Need to add option to process UniProt prior to EnsEMBL download.
    # [ ] : Tidy up code and use new rje_seqlist module more where possible.
    # [ ] : Change to rje_blast_V2.
    # [ ] : Replace uniprot spec file with rje_taxonomy.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_DBASE', '2.3', 'January 2013', '2007')
    description = 'RJE Database Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is under development and may contain bugs!',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
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
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)
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
### END OF SECTION I
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: DatabaseController                                                                                      #
#########################################################################################################################
class DatabaseController(rje.RJE_Object):     
    '''
    DatabaseController Class. Author: Rich Edwards (2007).

    Info:str
    - EnsPath = Path to EnsEMBL file [EnsEMBL/]
    - IPIPath = Path to IPI files [IPI/]
    - MakeDB = Makes a database from combined databases [None]
    - UniPath = Path to UniProt files [UniProt/]
    - TaxaDB = File containing details of taxanomic sub-databases to make [None]
        
    Opt:boolean
    - DatIndex = Create an index file (uniprot.index) for the Uniprot DAT files in unipath [True]
    - EnsFilter = Run EnsEMBL genomes through rje_seq to apply filters, rather than just concatenating [False]
    - Force = Whether to force regeneration of existing files [False]
    - FormatDB = Whether to BLAST format database after making [True]
    - IgnoreDate = Whether to ignore the relative timestamps of files when assessing whether to regenerate [False]
    - InverseDB = TaxaList is a list of taxanomic groups *NOT* to be in database [False]
    - ScreenIPI = species represented by IPI databases will be screened out of UniProt and EnsEMBL. [False]
    - ScreenEns = species represented by EnsEMBL databases will be screened out of UniProt. [True]
    - SeqFilter = Use rje_seq to filter sequences (True) or simply filter on Species Codes (False) [False]
    - SpecOnly = Will simply output a list of SPECIES codes to the makedb file, rather than making dbase [False]
    - SpecTable = Makes a table of species codes, taxonomy and taxon_id [True]
        - can then use "grep 'Metazoa' uniprot.spec.tdt | gawk '{print $1 }' | sort -u" to extract, e.g. Metazoan species codes
    - UseX = Whether to use certain aspects of databases,
        - UseUniprot,UseSprot,UseTrembl,UseEnsembl,UseKnown,UseNovel,UseAbinitio,UseIPI [all True]

    Stat:numeric

    List:list
    - DBDownload = List of EnsEMBL/UniProt/XML files containing details of databases to download []
    - DBFormat = Reformats UniProt, IPI or EnsEMBL databases using RJE_SEQ []
    - DBProcess = List of EnsEMBL/UniProt/IPI database types to process []
    - TaxaList = List of taxanomic groups to extract spec_codes for reduced database [None]
    
    Dict:dictionary    

    Obj:RJE_Objects
    - EnsEMBL = rje_ensembl.EnsEMBL object
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['MakeDB','UniPath','IPIPath','EnsPath','TaxaDB']
        self.optlist = ['DatIndex','SpecTable','UseUniprot','UseSprot','UseTrembl','UseEnsembl','UseKnown',
                        'UseNovel','UseAbinitio','UseIPI','SpecOnly','Force','IgnoreDate','InverseDB','ScreenIPI',
                        'ScreenEns','FormatDB','EnsFilter','SeqFilter']
        self.statlist = []
        self.listlist = ['DBDownload','DBFormat','DBProcess','TaxaList']
        self.dictlist = []
        self.objlist = ['EnsEMBL']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'UniPath':rje.makePath('UniProt/'),'IPIPath':rje.makePath('IPI/'),
                      'EnsPath':rje.makePath('EnsEMBL/')})
        for opt in ['UseUniprot','UseSprot','UseTrembl','UseEnsembl','UseKnown','UseNovel','UseAbinitio',
                    'ScreenEns','FormatDB','DatIndex','SpecTable']:
            self.opt[opt] = True
        self.obj['EnsEMBL'] = rje_ensembl.EnsEMBL(self.log,['ensloci=T']+self.cmd_list+['download=F'])
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
                self._cmdReadList(cmd,'info',['MakeDB'])
                self._cmdReadList(cmd,'file',['TaxaDB'])
                self._cmdReadList(cmd,'path',['UniPath','IPIPath','EnsPath'])
                self._cmdReadList(cmd,'opt',['DatIndex','SpecTable','UseUniprot','UseSprot','UseTrembl','SeqFilter',
                                             'UseEnsembl','UseKnown','UseNovel','UseAbinitio','UseIPI','SpecOnly',
                                             'Force','IgnoreDate','InverseDB','ScreenIPI','ScreenEns','EnsFilter'])
                self._cmdReadList(cmd,'list',['DBDownload','DBProcess','DBFormat','TaxaList'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
            self.list['DBFormat'] = string.split(string.join(self.list['DBFormat']).lower())
            if 'none' in self.list['DBFormat']: self.list['DBFormat'].remove('none')
#########################################################################################################################
    ### <2> ### Main Run Method                                                                                         #
#########################################################################################################################
    def run(self): self.mainRun()
    def mainRun(self):  ### Main DatabaseController Run Method
        '''
        Main DatabaseController Run Method. This controls the main methods for the DatabaseController class, including
        user-interaction where appropriate. The roles of rje_dbase are sorted into a hierarchy of precedence:
        1. Download new files.
        2. Process database files, starting with UniProt (needed for EnsLoci etc.)
        3. Reformat databases with RJE_SEQ if desired
        4. Make taxa-specific databases with TaxaDB
        5. Make custom databases with MakeDB
        '''
        try:### ~ [1] Download New files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.downloads()
            ### ~ [2] Process download files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.process()
            #!# >>>> After here may need editing <<<< #!#
            ### ~ [3] Check paths for additional reformatting and database generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.checkSetup()
            ### ~ [4]  Additional database functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Needs updating with new list >>> if self.opt['DBFormat']: self.reformatDB()  <<< #!#
            if self.info['TaxaDB'].lower() not in ['','none']: self.taxaDB()
            if self.info['MakeDB'].lower() not in ['','none']: self.makeDB()
        except: self.log.errorLog('Error in RJE_DBASE mainRun()',printerror=True)
#########################################################################################################################
    def checkSetup(self):   ### Check Files, Stages and Paths
        '''
        Check Files, Stages and Paths. All file and path checking will be done here, so there is no need in the
        individual sections? Attributes will be updated accordingly.
        '''
        try:### ~ [1] Check SpecTable needs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['TaxaList'] and not self.opt['SpecTable']:
                self.log.printLog('#CMD','SpecTable needed for TaxaList (spectable=T).')
                self.opt['SpecTable'] = True
            ### ~ [2] Check Paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            unipath = self.info['UniPath'] 
            ipipath = self.info['IPIPath'] 
            enspath = self.info['EnsPath'] 
            ## ~ [2a] UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (self.opt['SpecTable'] or self.opt['UseUniprot']) and not os.path.exists(unipath): rje.mkDir(self,unipath)
            while (self.opt['SpecTable'] or self.opt['UseUniprot']) and not os.path.exists(unipath):
                self.log.printLog('#ERR','UniProt path "%s" does not exist!' % unipath)
                if self.stat['Interactive'] >= 0:
                    self.info['UniPath'] = rje.choice('Enter new UniProt path (or keep to exit)',default=unipath)
                    if unipath == rje.makePath(self.info['UniPath'],return_blank=False): raise IOError
                    unipath = rje.makePath(self.info['UniPath'],return_blank=False)
                else: raise IOError
            ## ~ [2b] IPIPath ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ('ipi' in self.list['DBFormat'] or (self.info['MakeDB'] != 'None' and self.opt['UseIPI'])) and not os.path.exists(ipipath): rje.mkDir(self,ipipath)
            while ('ipi' in self.list['DBFormat'] or (self.info['MakeDB'] != 'None' and self.opt['UseIPI'])) and not os.path.exists(ipipath):
                self.log.printLog('#ERR','IPI path "%s" does not exist!' % ipipath)
                if self.stat['Interactive'] >= 0:
                    self.info['IPIPath'] = rje.choice('Enter new IPI path (or keep to exit)',default=ipipath)
                    if ipipath == rje.makePath(self.info['IPIPath'],return_blank=False): raise IOError
                    ipipath = rje.makePath(self.info['IPIPath'],return_blank=False)
                else: raise IOError
            ## ~ [2c] EnsPath ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ('ensembl' in self.list['DBFormat'] or (self.info['MakeDB'] != 'None' and self.opt['UseEnsembl'])) and not os.path.exists(enspath): rje.mkDir(self,enspath)
            while ('ensembl' in self.list['DBFormat'] or (self.info['MakeDB'] != 'None' and self.opt['UseEnsembl'])) and not os.path.exists(enspath):
                self.log.printLog('#ERR','EnsEMBL path "%s" does not exist!' % enspath)
                if self.stat['Interactive'] >= 0:
                    self.info['EnsPath'] = rje.choice('Enter new EnsEMBL path (or keep to exit)',default=enspath)
                    if enspath == rje.makePath(self.info['EnsPath'],return_blank=False): raise IOError
                    enspath = rje.makePath(self.info['EnsPath'],return_blank=False)
                else: raise IOError
        except:
            self.log.errorLog('Error in RJE_DBASE checkSetup()',printerror=True)
            raise   
#########################################################################################################################
    ### <3> ### Database Download Methods                                                                               #
#########################################################################################################################
    def downloads(self):    ### Reads list of databases from file, if given, and downloads if able
        '''Reads list of databases from file, if given, and downloads if able.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['DBDownload']: return
            dbdirlist = []          # Make a list of directories for later unzipping of downloaded files
            downloadXML = rje_xml.XML(self.log,self.cmd_list)   # XML parsing object

            ### ~ [2] Process download list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for download in self.list['DBDownload']:
                try:
                    ## ~ [2a] Deal with special cases of EnsEMBL and UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if download.lower() == 'ensembl':
                        (ensloci,self.obj['EnsEMBL'].opt['EnsLoci']) = (self.obj['EnsEMBL'].opt['EnsLoci'],False)
                        self.obj['EnsEMBL'].downloadEnsEMBL()
                        self.obj['EnsEMBL'].opt['EnsLoci'] = ensloci
                    elif download.lower() == 'uniprot':
                        rje_uniprot.downloadUniProt(self)
                        dbdirlist.append(self.info['UniPath'])
                    ## ~ [2b] Deal with databases found in XML file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    else:
                        downloadXML.list['XML'] = []
                        downloadXML.parseXML(download)
                        for xml in downloadXML.list['XML']:     # Each of these should be a database
                            ## Check Database ##
                            if xml.info['Name'] != 'database':
                                self.log.errorLog('"%s" read in where "database" expected!' % xml.info['Name'],printerror=False)
                                continue
                            ## Check OutDir and Enter ##
                            outdir = rje.makePath(xml.dict['Attributes']['outdir'])
                            if not os.path.exists(outdir): rje.mkDir(self,outdir)
                            dbdirlist.append(outdir)
                            dbtxt = '%s from %s -> %s' % (xml.dict['Attributes']['name'],xml.dict['Attributes']['ftproot'],outdir)
                            if self.stat['Interactive'] >= 1 and not rje.yesNo('Download %s?' % dbtxt): continue
                            self.log.printLog('#DBASE','Downloading %s: %d elements.' % (dbtxt,len(xml.list['XML'])))
                            mydir = os.path.abspath(os.curdir)
                            #!# Check for contents and backup #!#
                            ## Download elements ##
                            for element in xml.list['XML']:
                                if element.info['Name'] != 'file':
                                    self.log.errorLog('"%s" read in where "file" expected!' % element.info['Name'],printerror=False)
                                    continue
                                wget = rje.makePath(xml.dict['Attributes']['ftproot']+element.dict['Attributes']['path'],wholepath=True)
                                self.log.printLog('#DBASE','Downloading %s %s' % (xml.dict['Attributes']['name'],element.info['Content']))
                                os.chdir(outdir)
                                if self.getBool('OSX'): os.system('curl -O %s' % wget)
                                else: os.system('wget %s' % wget)
                                if os.path.exists('gene_ontology_edit.obo'): os.rename('gene_ontology_edit.obo','GO_OBO_%s.txt' % rje.dateTime(yymmdd=True))
                                os.chdir(mydir)
                except: self.errorLog('Problem with "%s" database download.' % download)

            ### ~ [3] Unzip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Win32']: self.log.printLog('#WIN','Cannot unzip in Windows, sorry.')
            else:
                for dir in dbdirlist:   ## *.gz ##
                    zipfiles = rje.getFileList(callobj=self,folder=dir,filelist=['*.gz'],subfolders=False,summary=False)
                    self.deBug(zipfiles)
                    gtxt = 'Unzipping %d files in %s' % (len(zipfiles),dir)
                    for gz in zipfiles:
                        self.log.printLog('\r#GZ','%s: %s        ' % (gtxt,os.path.basename(gz)),log=False,newline=False)
                        os.system('gunzip %s' % gz)
                    self.log.printLog('\r#GZ','%s complete.        ' % (gtxt))
            ### ~ [4] Finished! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#DBASE','Downloading of Databases complete.')
        except:
            self.log.errorLog('Error in DatabaseController.downloads()',quitchoice=True)
            return False
#########################################################################################################################
    ### <4> ### Process/Reformat Database Methods                                                                       #
#########################################################################################################################
    def process(self):  ### Processes desired databases from self.list['DBProcess']
        '''Processes desired databases from self.list['DBProcess'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['DBProcess']: return
            dblist = string.split(string.join(self.list['DBProcess']).lower())
            ### ~ [2] Process databases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(dblist)
            ## ~ [2a] UniProt should be processed first! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self._setForkAttributes()   # Delete if no forking
            for cmd in self.cmd_list: self._forkCmd(cmd)
            #x#self.stat['Forks'] = 8
            if 'uniprot' in dblist: rje_uniprot.processUniProt(self,makeindex=self.opt['DatIndex'],makespec=self.opt['SpecTable'],makefas=True)
            ## ~ [2b] Next EnsEMBL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'ensembl' in dblist:
                self.obj['EnsEMBL'].run()
                self.obj['EnsEMBL'].opt['Force'] = False    # May call EnsLoci again and do not want to recreate! #
            ## ~ [2c] Last, IPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'ipi' in dblist: self.reformatIPI()
        except:
            self.log.errorLog('Error in DatabaseController.process()',quitchoice=True)
            return False
#########################################################################################################################
    def reformatIPI(self):  ### Reformats raw IPI files
        '''Reformats raw IPI files.'''
        try: #Clean this up!#
            ipipath = self.info['IPIPath']
            ipifiles = []
            ipispec = []
            if self.opt['UseIPI']:
                ipifiles = glob.glob('%sipi.*.fasta' % ipipath)
            for ipi in ipifiles:
                newipi = rje.baseFile(ipi) + '.new'
                species = string.split(ipi.upper(),'.')[1]
                ipispec.append(species)
                outipi = ipipath + 'ipi_%s.fas' % species
                if not self.opt['Force'] and rje.isYounger(outipi,ipi) == outipi:
                    self.log.printLog('#IPI','%s already exists and younger than %s (Force=F)' % (outipi,ipi))
                    continue
                INFILE = open(ipi,'r')
                OUTFILE = open(newipi,'w')
                sx = 0
                name = ''
                sequence = ''
                nextname = ''
                while 1:
                    line = INFILE.readline()
                    if not line:
                        if name:
                            OUTFILE.write('>%s\n%s\n' % (name,sequence))
                            sx += 1
                        break
                    line = rje.chomp(line)
                    if line[:1] == '>':     # Description line
                        if name:
                            OUTFILE.write('>%s\n%s\n' % (name,sequence))
                            sx += 1
                            self.log.printLog('\r#SEQ','%s sequences reformatted from %s.' % (rje.integerString(sx),ipi),log=False,newline=False)
                        # Reformat
                        newname = string.replace(line,'|',' ')
                        newname = string.replace(newname,';',' ; ') # Enables matching of first AccNum only.
                        accs = {'ipi':rje.matchExp('IPI:(\S+)',newname)[0]}
                        if rje.matchExp('SWISS-PROT:(\S+)',newname):
                            accs['sprot'] = rje.matchExp('SWISS-PROT:(\S+)',newname)[0]
                        if rje.matchExp('TREMBL:(\S+)',newname):
                            accs['trembl'] = rje.matchExp('TREMBL:(\S+)',newname)[0]
                        if rje.matchExp('ENSEMBL:(\S+)',newname):
                            accs['ens'] = rje.matchExp('ENSEMBL:(\S+)',newname)[0]
                        if rje.matchExp('REFSEQ:(\S+)',newname):
                            accs['refseq'] = rje.matchExp('REFSEQ:(\S+)',newname)[0]
                        #!# Look up case-insensitive matching #!# s/^.+tax_id=\d+//i;  # Description only
                        t = newname.find('Tax_Id=') + len('Tax_Id=')
                        newname = rje.matchExp('%s\d+\s+(\S.*)$' % newname[:t],newname)[0]
                        newname = string.replace(newname,' ; ',';')     
                        newacc = 'ipi_%s__' % species
                        if accs.has_key('sprot'):
                            newacc = '%s%s ' % (newacc,accs['sprot'])
                        elif accs.has_key('trembl'):
                            newacc = '%s%s ' % (newacc,accs['trembl'])
                        if accs.has_key('ens'):
                            newacc = '%s%s ' % (newacc,accs['ens'])
                        if accs.has_key('refseq'):
                            newacc = '%s%s ' % (newacc,accs['refseq'])
                        newacc = '%s%s ' % (newacc,accs['ipi'])
                        # Replace #
                        name = '%s%s' % (newacc,newname)
                        sequence = ''
                    else:
                        sequence = '%s%s' % (sequence,line)
                INFILE.close()
                OUTFILE.close()
                self.log.printLog('\r#SEQ','%s sequences reformatted from %s.' % (rje.integerString(sx),ipi))
                seqcmd = ['reformat=fas','seqout=%s' % outipi,'seqin=%s' % newipi,'memsaver=T','autoload=T','autofilter=T']
                rje_seq.SeqList(self.log,self.cmd_list+seqcmd)
                os.unlink(newipi)
            if ipispec:
                SPECLIST = open(ipipath + 'ipi.spec_code','w')
                SPECLIST.write(string.join(ipispec,'\n'))
                SPECLIST.close()
        except: self.log.errorLog('Error in rje_dbase.reformatIPI',printerror=True,quitchoice=False)
#########################################################################################################################
    def reformatDB(self):      ### This method is for reformatting raw IPI, EnsEMBL and UniProt files using rje_seq
        '''
        This method is for reformatting raw IPI, EnsEMBL and UniProt files. For all other files, use rje_seq.py.
        '''
        try:### ~ [1] UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            unipath = self.info['UniPath']
            unifiles = []
            if self.opt['UseUniprot'] and self.opt['UseSprot']: unifiles.append(unipath + 'uniprot_sprot.dat')  #!# Check this
            if self.opt['UseUniprot'] and self.opt['UseTrembl']: unifiles.append(unipath + 'uniprot_trembl.dat')  #!# Check this
            for uni in unifiles:
                uniout = rje.baseFile(uni) + '_all.fas'
                if self.opt['Force'] or rje.isYounger(uniout,uni) != uniout:
                    seqcmd = ['reformat=fasta','seqout=%s' % uniout,'seqin=%s' % uni,'memsaver=T','gnspacc=T']
                    rje_seq.SeqList(self.log,self.cmd_list+seqcmd)
                else: self.log.printLog('#UNI','%s already exists and younger than %s (Force=F)' % (uniout,uni))

            ### ~ [2] IPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.reformatIPI()
                
            ### ~ [3] EnsEMBL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['UseEnsembl']:
                usetypes = []
                for type in ['known','novel','abinitio']:
                    checktype = type[:1].upper() + type[1:].lower()
                    if self.opt['Use%s' % checktype]:
                        usetypes.append(type)
                        if type == 'known': usetypes.append('known-ccds')
                    else: self.log.printLog('#DB','*.%s.fa skipped: Use%s=False.' % (type,checktype))
                self.obj['EnsEMBL'].reformatDB(usetypes=usetypes)
        except:
            self.log.errorLog('Error in rje_dbase.reformatDB',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <5> ### Make Compiled Database Methods                                                                          #
#########################################################################################################################
    def taxaDB(self):    ### Reads list of taxa makedb from file, if given, and downloads if able
        '''Reads list of databases from file, if given, and downloads if able.'''
        try:
            ### Setup ###
            rje.mkDir(self,'TaxaDB/')
            if not os.path.exists(self.info['TaxaDB']):
                self.log.printLog('#TAXA','File "%s" not found. Will use taxadb.xml' % self.info['TaxaDB'])
                self.info['TaxaDB'] = 'taxadb.xml'
            force = self.opt['Force']
            
            ### Read Databases to Make into XML object ###
            taxadbXML = rje_xml.XML(self.log,self.cmd_list)
            taxadbXML.parseXML(self.info['TaxaDB'])
            taxaxml = taxadbXML     #.list['XML'][0]
            if taxaxml.info['Name'] != 'taxadb':
                self.log.errorLog('"%s" read in where "taxadb" expected!' % taxaxml.info['Name'],printerror=False)
                return False

            ### Work through and download Databases ###
            for xml in taxaxml.list['XML']:     # Each of these should be a database
                ## Check Database ##
                if xml.info['Name'] != 'makedb':
                    self.log.errorLog('"%s" read in where "makedb" expected!' % xml.info['Name'],printerror=False)
                    continue
                ## MakeDB ##
                self.printLog('#DBASE','Generating TaxaDB "%s"' % xml.info['Content'])
                makecmd = self.cmd_list[0:]
                for key in xml.dict['Attributes']: makecmd.append('%s=%s' % (key,xml.dict['Attributes'][key]))
                #x#print makecmd
                dbase = DatabaseController(self.log,makecmd)
                dbase.setInfo(self.info)
                dbase.setOpt(self.opt)
                dbase._cmdList()
                dbase.setInfo({'Download':'None','TaxaDB':'None'})
                dbase.setOpt({'DBFormat':False,'Index':False,'SpecTable':True,'Force':self.opt['Force']})
                dbase.obj['EnsEMBL'] = self.obj['EnsEMBL']
                dbase.checkSetup()
                dbase.makeDB()
            self.opt['Force'] = force
        except:
            self.log.errorLog('Error in taxaDB()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def makeDB(self):      ### This compiles different databases into a concatenated database use taxanomic groups etc.
        '''
        This compiles different databases into a concatenated database use taxanomic groups etc. By default, species
        represented by EnsEMBL databases will be screened out of UniProt. This can be switched off with screenens=F.
        Any sequence filtering options (such as minumum lengths etc.) will be applied to these new databases.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            if self.opt['Win32']:
                self.log.printLog('#WIN32','Cannot use makeDB() method in Windows: UNIX only.')
                return False
            unipath = self.info['UniPath']
            enspath = self.info['EnsPath']

            ### Get Species codes ###
            _stage = 'Species codes'
            specfile = self.makeSpecFile()
            if self.opt['SpecOnly']: return True
            SPEC = open(specfile,'r')
            speclist = string.split(SPEC.read(),'\n')
            SPEC.close()
            self.log.printLog('#SPEC','%s Species codes for makeDB(%s).' % (rje.integerString(len(speclist)),rje.baseFile(self.info['MakeDB'])))
            
            ### Make Reduced UniProt ###
            uniscreen = ''
            if self.opt['UseUniprot']:
                ## Get species (if any) to screen ##
                screenspec = []
                uni_ext = 'combined.fas'
                if self.opt['ScreenEns'] and self.opt['UseEnsembl']:
                    screenspec = [] # rje_ensembl.enspec.values()
                    for ensloci in glob.glob('%sens_*.loci.fas' % enspath):
                        try: screenspec.append(rje.matchExp('ens_(\S+).loci.fas',ensloci)[0])
                        except: self.errorLog(ensloci)
                    self.printLog('#ENS','%d EnsEMBL species to screen' % len(screenspec))
                    uni_ext = 'no_ens.fas'
                if self.opt['ScreenIPI'] and self.opt['UseIPI'] and os.path.exists(self.info['IPIPath'] + 'ipi.spec_code'):
                    screenspec = string.split(open(self.info['IPIPath'] + 'ipi.spec_code','r').read(),'\n')
                    uni_ext = 'no_ipi.fas'
                ## Produce screened/filtered combined UniProt ##
                uniscreen = unipath + 'uniprot.%s' % uni_ext
                unifiles = []
                if self.opt['UseUniprot'] and self.opt['UseSprot']:
                    unifiles.append(unipath + 'uniprot_sprot_all.fas')      #!# Check this
                if self.opt['UseUniprot'] and self.opt['UseTrembl']:
                    unifiles.append(unipath + 'uniprot_trembl_all.fas')     #!# Check this
                ## Check Date ##
                makescreen = self.opt['Force']
                for file in unifiles:
                    if rje.isYounger(uniscreen,file) != uniscreen: makescreen = True
                if makescreen:
                    if unifiles: rje.backup(self,uniscreen)
                    else:
                        self.log.errorLog('Cannot find uniprot_sprot_all.fas or uniprot_trembl_all.fas in %s' % unipath,printerror=False)
                        return False
                    if screenspec:
                        self.printLog('#UNIPROT','Making species-screened UniProt file "%s"' % uniscreen)
                        for file in unifiles:
                            seqcmd = ['seqin=%s' % file,'seqout=%s' % uniscreen,'reformat=fasta','gnspacc=T','seqmode=file',
                                      'logrem=F','append=T','badspec=%s' % string.join(screenspec,','),'autoload=T']
                            rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd)
                    else:
                        self.printLog('#UNIPROT','Making combined UniProt file "%s"' % uniscreen)
                        os.system('cat %s > %s' % (string.join(unifiles),uniscreen))
                else: self.printLog('#DB','%s already exists and younger that UniProt files (force=F)' % uniscreen)

            ### Get files to make DB from ###
            _stage = 'MakeDB Files'
            dbfiles = []
            ## IPI ##
            ipipath = self.info['IPIPath']
            if self.opt['UseIPI']:
                for ipi in glob.glob('%sipi_*.fas' % ipipath):
                    species = rje.matchExp('ipi_(\S+)\.fas',ipi)[0].upper()
                    ## Check for species in specfile ##
                    if (not self.opt['InverseDB'] and species in speclist) or (self.opt['InverseDB'] and species not in speclist):
                        dbfiles.append(ipi)
            ## EnsEMBL ##
            if self.opt['UseEnsembl']:  #!# Subsets of EnsEMBL now chosen during reformatDB #!#
                for species in rje.sortKeys(rje_ensembl.enspec):
                    code = rje_ensembl.enspec[species]
                    if code not in screenspec: self.printLog('#SPEC','%s (%s) not found in EnsEMBL/' % (species,code)); continue
                    if (not self.opt['InverseDB'] and code in speclist) or (self.opt['InverseDB'] and code not in speclist):
                        if self.obj['EnsEMBL'].opt['EnsLoci']:
                            dbfiles.append(enspath + 'ens_%s.loci.fas' % code)
                        else:
                            dbfiles.append(rje.makePath(enspath + code) + 'ens_%s.all.fas' % code)
            ## EnsFilter ##
            catcmd = ''
            if not self.opt['EnsFilter'] and not self.opt['Win32'] and dbfiles:
                #catcmd = 'cat %s > %s' % (string.join(dbfiles),self.info['MakeDB'])
                taxacode = string.split(rje.baseFile(self.info['MakeDB'],True),'_')[-1]
                ens_taxadb = rje.makePath(enspath) + 'ens_%s.fas' % taxacode
                catcmd = 'cat %s > %s' % (string.join(dbfiles),ens_taxadb)
            ## UniProt ##
            if uniscreen: dbfiles.append(uniscreen)
            ## Check whether run is needed ##
            dbexists = os.path.exists(self.info['MakeDB']) and not self.force()
            for file in dbfiles:
                dbexists = dbexists and rje.isYounger(self.info['MakeDB'],file) == self.info['MakeDB']
            if dbexists: 
                self.printLog('#MAKE','Output file "%s" exists and younger than source files (force=F)' % self.info['MakeDB'])
                return True
            dbfiles = []
            if uniscreen: dbfiles.append(uniscreen)
            elif self.opt['UseUniprot']: dbfiles = unifiles[0:]
            
            ### Make new database ###
            self.deBug(self.info['MakeDB'])
            _stage = 'MakeDB'
            rje.backup(self,self.info['MakeDB'])
            if catcmd:
                self.printLog('#CAT',catcmd)
                os.system(catcmd)
                open(self.info['MakeDB'],'w').write(open(ens_taxadb,'r').read())
                self.printLog('#COPY','Copied %s to %s' % (ens_taxadb,self.info['MakeDB']))
            seqcmd = ['seqout=%s' % self.info['MakeDB'],'reformat=fasta','gnspacc=T','seqmode=file','logrem=F','append=T','autoload=T']
            if self.opt['InverseDB']:
                seqcmd.append('badspec=%s' % specfile)
            else:
                seqcmd.append('goodspec=%s' % specfile)
            for file in dbfiles:
                rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd+['seqin=%s' % file])

            ### BLAST Format and Finish ###
            if self.opt['FormatDB']:
                myblast = rje_blast.BLASTRun(self.log,self.cmd_list+['blastd=%s' % self.info['MakeDB']])
            return True
            
        except:
            self.log.errorLog('Error in makeDB(%s)' % _stage,printerror=True,quitchoice=True)
            return False
#########################################################################################################################
    def makeSpecFile(self):      ### This extracts species codes for the desired taxanomic groups and returns specfile
        '''
        This extracts species codes for the desired taxanomic groups and returns specfile.
        << specfile:str = File name containing species codes, or None if failed.
        '''
        try:
            ### Setup ###
            unipath = self.info['UniPath']
            spectable = unipath + 'uniprot.spec.tdt'
            if self.opt['SpecOnly']:
                specfile = self.info['MakeDB']
            else:
                specfile = rje.baseFile(self.info['MakeDB']) + '.spec_code'
            if not self.opt['Force'] and rje.isYounger(specfile,spectable) == specfile:
                self.log.printLog('#SPEC','%s exists and younger than %s (force=F)' % (specfile,spectable))
                return specfile
            tempfile = '%s.tmp' % rje.randomString(10)

            ### Extract All lines into Temp File ###            
            TEMP = open(tempfile,'w')
            for taxa in self.list['TaxaList']:
                #X# TAXA = os.popen("grep '%s' %s | gawk '{print $1 }' | sort -u" % (taxa,specfile))
                #X# TEMP.write(string.join(TAXA.readlines(),''))
                #X# TAXA.close()
                #!# Gawk not installed on bioware. Therefore need a fix #!#
                if taxa == taxa.upper(): taxalines = os.popen("grep '%s' %s " % (taxa,spectable)).readlines()
                else: taxalines = os.popen("grep -i '%s' %s " % (taxa,spectable)).readlines()
                for t in taxalines:
                    #X#print t, string.split(t)[0]
                    TEMP.write('%s\n' % string.split(t)[0])
            TEMP.close()
            self.deBug('%s made.' % tempfile)

            ### CleanUp ###
            if os.system('sort -u %s > %s' % (tempfile,specfile)):
                self.printLog('#ERR','Something has gone wrong during grep and/or sort')
                raise ValueError
            os.unlink(tempfile)     # Species list now in specfile
            self.log.printLog('#SPEC','Species codes for %d taxanomic groups output to %s.' % (len(self.list['TaxaList']),specfile))
            return specfile
        except:
            self.log.errorLog('Error in makeSpecFile(%s)' % self.info['MakeDB'],printerror=True,quitchoice=True)
            return None
#########################################################################################################################
## End of SECTION II: DatabaseController                                                                                #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain(): ### ~ Basic Setup of Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
    ### ~ Primary Functionality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:  DatabaseController(mainlog,['forks=8']+cmd_list).mainRun()
    ### ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG','%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV
#########################################################################################################################

