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
Module:       rje_genecards
Description:  RJE Genecards Parsing Module 
Version:      0.4
Last Edit:    28/03/08
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

NOTE:
    This module has now been superceded somewhat by the rje_genemap module but is still used with rje_hprd to compile
    links from HPRD. This module may also still be of use for smaller sets of genes that need to me mapped to HGNC, e.g.
    manually compiled lists from experiments.

Function:
    This is a prototype module, which aims to take in a list of Gene Symbols and/or Aliases, find the relevant GeneCard
    entry, download it and extract the relevant protein gene/protein links into a table.

    The ultimate goal is to generate a table pulling in identifiers from EnsEMBL, GeneCards and HPRD to allow easy
    cross-referencing across datasets and compilation of data from different sources. When using the altsource=LIST
    option, subsequent files will overwrite the data read from files earlier in the list. If update=T and the cardout
    file exists, this will be appended to the altsource list.

    To save time, a full download of HGNC symbols can be downloaded from HGNC (http://www.genenames.org/index.html)
    and imported using the hgncdata=FILE option. This file should be delimited and contain the following fields (others
    are allowed):
    - HGNC ID, Approved Symbol, Approved Name, Previous Symbols, Aliases, Entrez Gene ID, RefSeq IDs,
    Entrez Gene ID (mapped data), OMIM, UniProt ID (mapped data), Ensembl ID (mapped data)

Commandline:
    ### Input Options ###
    genes=LIST      : List of gene symbols/aliases to download []
    update=T/F      : Whether to read in any data from cardout file (if present) and add to it [True]
    skiplist=LIST   : Skip genes matching LIST (e.g. XP_*) []
    useweb=T/F      : Whether to try and extract missing data from GeneCards website [True]
    altsource=LIST  : List of alternative sources of data (Delimited files with appropriate headers) []
    hgncdata=FILE   : HGNC download file containing data []

    ### Output Options ###
    species=X       : Species to output in table [Human]
    cardout=FILE    : File for output of genecard data [genecards.tdt]
    ensloci=FILE    : File of EnsLoci genome to incorporate [/home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas]
    restrict=T/F    : Whether to only output lines for gene in the original gene=LIST [False]
    purify=T/F      : Only output lines where the Alias and the Symbol are the same [False]

    ### Special execution options ###
    fullens=T/F     : Incorporate all EnsLoci EnsEMBL genes into cardout file (long run!) [False]
    fullhgnc=T/F    : Output all HGNC codes and unambiguous aliases into file [False]

Uses general modules: copy, glob, os, string, sys, time, urllib2
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time, urllib2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added EnsLoci processing.
    # 0.2 - Added altsource and generally improved function and commenting.
    # 0.3 - Added more interactivity and options.
    # 0.4 - Added reading of HGNC download.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Link to EnsEMBL Human sequence.
    # [Y] : Use HPRD ID information if present. (Read from a table of HPRD GeneCards to use in favour of website.)
    # [Y] : Get the fullens option working properly.
    # [ ] : Add option to reduce output to one-line per entry (i.e. if Alias elsewhere in fields for symbol)
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_GENECARDS', '0.4', 'March 2008', '2007')
    description = 'RJE GeneCards/HGNC Parsing Module'
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
### Constants ###
gc_headers = ['Symbol','HGNC','Entrez','UniProt','EnsEMBL','HPRD','OMIM']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GeneCards Class                                                                                         #
#########################################################################################################################
class GeneCards(rje.RJE_Object):     
    '''
    GeneCards Class. Author: Rich Edwards (2005).

    Info:str
    - CardOut = File for output of genecard data [genecards.tdt]
    - EnsLoci = File of EnsLoci-processed genome to incorporate [None]
    - HGNCData = HGNC download file containing data []
    - Species = Species to output in table [Human]
    
    Opt:boolean
    - FullEns = Incorporate all EnsLoci EnsEMBL genes into cardout file (long run!) [False]
    - FullHGNC = Output all HGNC codes and unambiguous aliases into file [False]
    - Purify = Only output lines where the Alias and the Symbol are the same [False]
    - Restrict = Whether to only output lines for gene in the original gene=LIST [False]
    - Update = Whether to read in any data from cardout file (if present) and add to it [True]
    - UseWeb = Whether to try and extract missing data from GeneCards website [True]

    Stat:numeric

    List:list
    - AltSource = List of alternative sources of data (Delimited files with appropriate headers) []
    - Genes = Input List of Genes/Aliases to extract
    - SkipList = Skip genes matching LIST  []
    - TestGenes = List of gene IDs for special debugging attention

    Dict:dictionary
    - CardMap = Dictionary of {Alias:Symbol} where alias can be gene symbol, UniProt, EntrezGene etc.
    - EnsLoci = Dictionary of {EnsGene:shortName()}
    - EnsDesc = Dictionary of {EnsGene:Description}
    - GeneCard = Dictioary of input Gene Symbol and relevant extracted data

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['CardOut','EnsLoci','HGNCData','Species']
        self.optlist = ['FullEns','FullHGNC','Update','Purify','Restrict','UseWeb']
        self.statlist = []
        self.listlist = ['AltSource','Genes','SkipList','TestGenes']
        self.dictlist = ['CardMap','EnsDesc','EnsLoci','GeneCard']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'CardOut':'genecards.tdt','Species':'Human',
                      'EnsLoci':rje.makePath('/home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas',True)})
        self.setOpt({'FullEns':False,'Purify':False,'Restrict':False})
        #x#self.list['SkipList'] = string.split('XP_*,NP_*,IPI00*',',')
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
                self._cmdReadList(cmd,'file',['CardOut','EnsLoci','HGNCData'])
                self._cmdReadList(cmd,'info',['Species'])
                self._cmdReadList(cmd,'opt',['FullEns','FullHGNC','Update','Purify','Restrict','UseWeb'])
                self._cmdReadList(cmd,'list',['Genes','SkipList','TestGenes'])
                self._cmdReadList(cmd,'glist',['AltSource'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def setup(self):    ### Sets up headers and reads in existing data if present
        '''Sets up headers and reads in existing data if present.'''
        try:
            ### ~ Setup Basic Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #X#headers = ['Alias','Species','Symbol','HGNC','Entrez','UniProt','EnsEMBL','HPRD','OMIM','EnsLoci','Desc']
            headers = ['Alias','Species'] + gc_headers  # All other headers added from altsource list
            ### ~ Read in data from existing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.readHGNC()
            if self.opt['Update'] and os.path.exists(self.info['CardOut']): self.list['AltSource'].append(self.info['CardOut'])
            for altsource in self.list['AltSource']:
                sourcefile = rje.makePath(altsource,True)
                if not os.path.exists(sourcefile):
                    self.log.errorLog('Alternative source "%s" missing!' % sourcefile,printerror=False,quitchoice=True)
                    continue
                update = rje.dataDict(self,sourcefile,getheaders=True,ignore=['#'])
                for h in update.pop('Headers'):
                    if h not in headers:
                        headers.append(h)
                self.log.printLog('#DATA','Read GeneCards data for %d genes.' % (len(update)))
                for gene in rje.sortKeys(update):     # Each source will overwrite data from the file before
                    ## ~ Convert to Upper Case for consistency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if gene != gene.upper() and gene.upper() in update: continue    # Only use upper case one!
                    elif gene != gene.upper():
                        update[gene.upper()] = update.pop(gene)
                        gene = gene.upper()
                    if gene == '!FAILED!': continue
                    ## ~ Update main dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if self.opt['Update'] and altsource == self.info['CardOut'] and gene not in self.list['Genes']: self.list['Genes'].append(gene)
                    if gene in self.dict['GeneCard']: rje.combineDict(self.dict['GeneCard'][gene],update[gene])
                    else: self.dict['GeneCard'][gene] = update[gene]
                    ## ~ Temp Debugging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if gene in self.list['TestGenes']:
                        print gene
                        print update[gene]
                        self.deBug(self.dict['GeneCard'][gene])
                    ## ~ Check Aliases etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if 'Symbol' in self.dict['GeneCard'][gene]: self.dict['GeneCard'][gene]['Symbol'] = self.dict['GeneCard'][gene]['Symbol'].upper()
                    if 'Symbol' in update[gene] and update[gene]['Symbol'] != '!FAILED!':
                        symbol = update[gene]['Symbol']
                        if symbol in self.dict['GeneCard']: rje.combineDict(self.dict['GeneCard'][symbol],update[gene],overwrite=False,replaceblanks=True)
                        else: self.dict['GeneCard'][symbol] = update[gene]
                    self.log.printLog('\r#CARD','Extracted GeneCards data for %d genes.' % (len(self.dict['GeneCard'])),newline=False,log=False)
                    if len(string.split(gene)) > 1: print '!!!', gene, '!!!'
            ### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('\r#CARD','Extracted GeneCards data for %d genes.' % (len(self.dict['GeneCard'])))
            self.list['Headers'] = headers[0:]
            if self.opt['Update']: self.opt['Append'] = False
            #x#if 'TASP1' in self.dict['GeneCard']: self.deBug(self.dict['GeneCard']['TASP1'])
            #x#else: self.deBug(rje.sortKeys(self.dict['GeneCard']))
        except:
            self.log.errorLog('Problem during GeneCards.setup()')
            raise
#########################################################################################################################
    def run(self,setup=True):  ### Main Run Method
        '''
        Main Run Method
        >> setup:bool [True] = Sets up headers and reads in existing data if present.
        '''
        try:
            ### ~ Setup & Read existing data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if setup: self.setup()
            headers = self.list['Headers']
            delimit = rje.delimitFromExt(filename=self.info['CardOut'])
            if os.path.exists(self.info['EnsLoci']):
                for h in ['EnsLoci','EnsDesc']:
                    if h not in headers: headers.append(h)
            rje.delimitedFileOutput(self,self.info['CardOut'],headers,delimit,rje_backup=True)

            ### ~ Read EnsLoci for incorporation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.ensLoci()
                        
            ### ~ Parse data from GeneCards website and/or previously read aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.processGenes(self.list['Genes'])
            self.interactiveUpdate()
        
            ### ~ Add EnsEMBL EnsLoci data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.addEnsLoci()

            ### ~ Output GeneCards data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.outputCards()
            
        except:
            self.log.errorLog('Apocalyptic error with GeneCards.run()')
            raise
#########################################################################################################################
    ### <3> ### Parsing data methods                                                                                    #
#########################################################################################################################
    def readHGNC(self):     ### Read links from HGNC into data structure
        '''Read links from HGNC into data structure.'''
        try:### ~ [1] Read into dictionary with HGNC ID as key ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['HGNCData'].lower() in ['','none']: return
            if not os.path.exists(self.info['HGNCData']): return self.log.errorLog('HGNC file "%s" not found' % (self.info['HGNCData']),printerror=False)
            hgncdata = rje.dataDict(self,self.info['HGNCData'],['HGNC ID'])
            aliaii = {}     # Dictionary of withdrawn symbols to map

            ### ~ [2] Parse out information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot) = (0.0,len(hgncdata))
            for hgnc in rje.sortKeys(hgncdata):
                self.log.printLog('\r#HGNC','Processing HGNC: %.1f%%' % (hx/htot),newline=False,log=False)
                hx += 100.0
                ## ~ [2a] Adjust headers for new vs old HGNC compatibility ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                data = hgncdata[hgnc]
                for hkey in rje.sortKeys(data):
                    if rje.matchExp('^(\S.+\S)\s*\(mapped data supplied by \S+\)',hkey):
                        data['%s (mapped data)' % rje.matchExp('^(\S.+\S)\s*\(mapped data supplied by \S+\)',hkey)[0]] = data.pop(hkey)
                    if rje.matchExp('^(\S.+\S)\s*\(supplied by \S+\)',hkey):
                        data['%s (mapped data)' % rje.matchExp('^(\S.+\S)\s*\(supplied by \S+\)',hkey)[0]] = data.pop(hkey)
                ## ~ [2b] Make dictionary of Genecards data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gdict = {}
                gdict['Symbol'] = gene = data['Approved Symbol'].upper()
                gdict['Desc'] = data['Approved Name']
                ## ~ [2c] Special treatment of obselete symbol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if gene.find('~withdrawn') > 0:     ### Obselete symbol
                    try:
                        gene = gene[:gene.find('~WITHDRAWN')]
                        alias = rje.matchExp(', see (\S+)',gdict['Desc'])[0]
                        if len(string.split(alias)) > 1: continue   # Ambiguous
                        if gene in aliaii and aliaii[gene] != alias: aliaii[gene] = 'AMBIGUOUS'
                        else: aliaii[gene] = alias
                    except: pass
                    continue
                ## ~ [2d] Add additional aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if 'Synonyms' in data and 'Aliases' not in data: data['Aliases'] = data.pop('Synonyms')
                for alias in string.split(data['Previous Symbols'].upper(),', ') + string.split(data['Aliases'].upper(),', '):
                    #x# if alias.upper() != alias: continue     # Not really a symbol
                    if alias in self.dict['GeneCard']: aliaii[alias] = 'AMBIGUOUS'
                    if alias in aliaii and aliaii[alias] != gene: aliaii[alias] = 'AMBIGUOUS'
                    else: aliaii[alias] = gene
                if gene in aliaii: aliaii[gene] = 'AMBIGUOUS'
                gdict['Entrez'] = data['Entrez Gene ID']
                if not gdict['Entrez']: gdict['Entrez'] = data['Entrez Gene ID (mapped data)']
                gdict['OMIM'] = data['OMIM ID (mapped data)']
                gdict['UniProt'] = data['UniProt ID (mapped data)']
                gdict['EnsEMBL'] = ensgene = data['Ensembl ID (mapped data)']
                gdict['HGNC'] = string.replace(hgnc,'HGNC:','')
                if not gene: gene = ensgene
                if not gene:
                    self.log.errorLog('HGNC has no gene for %s: %s' % (gdict['HGNC'],data),printerror=False)
                    continue
                #x#self.deBug(data)
                self.dict['GeneCard'][gene] = gdict
                if self.opt['FullHGNC'] and gene not in self.list['Genes']: self.list['Genes'].append(gene)
                ## ~ [2b] Deal with EnsGene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.opt['FullEns'] and ensgene:
                    if ensgene not in self.list['Genes']: self.list['Genes'].append(ensgene)
                    if ensgene not in self.dict['GeneCard']: self.dict['GeneCard'][ensgene] = {}
                    rje.combineDict(self.dict['GeneCard'][ensgene],gdict,overwrite=False,replaceblanks=True)
            #x#self.deBug(aliaii)
            self.log.printLog('\r#HGNC','Processed HGNC: %s genes & %s aliases' % (rje.integerString(len(self.dict['GeneCard'])),rje.integerString(len(aliaii))))

            ### ~ [3] Deal with aliaii ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ambig = []
            (hx,htot) = (0.0,len(aliaii))
            for alias in aliaii:
                self.log.printLog('\r#HGNC','Processing HGNC aliases: %.1f%%' % (hx/htot),newline=False,log=False)
                hx += 100.0
                if aliaii[alias] == 'AMBIGUOUS':
                    ambig.append(alias)
                    continue       # Alias mapped to multiple genes 
                while aliaii[alias] not in self.dict['GeneCard'] and aliaii[alias] in aliaii: aliaii[alias] = aliaii[aliaii[alias]]     # Map through several aliases if needed
                if aliaii[alias] not in self.dict['GeneCard']: continue     # Alias is not a valid Gene, so ignore
                if alias not in self.dict['GeneCard']: self.dict['GeneCard'][alias] = self.dict['GeneCard'][aliaii[alias]]
                if self.opt['FullHGNC'] and alias not in self.list['Genes']: self.list['Genes'].append(alias)
            self.log.printLog('\r#HGNC','Processed HGNC: %s genes & aliases' % (rje.integerString(len(self.dict['GeneCard']))))
            if ambig:
                self.log.printLog('#AMB','%s ambiguous aliases were not mapped' % rje.integerString(len(ambig)))
                open('hgnc.ambiguities.txt','w').write(string.join(ambig,'\n'))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def ensLoci(self):  ### Reads from EnsLoci file if it exists and parses into dictionaries.
        '''Reads from EnsLoci file if it exists and parses into dictionaries.'''
        self.dict['EnsLoci'] = {}    # Dictionary of {EnsGene:shortName()}
        self.dict['EnsDesc'] = {}    # Dictionary of {EnsGene:Description}
        self.dict['UniEns'] = {}     # Dictionary of {UniProt?:EnsGene}
        if os.path.exists(self.info['EnsLoci']):
            elines = self.loadFromFile(self.info['EnsLoci'])
            (ex,etot) = (0.0,len(elines))
            while elines:
                ex += 100.0
                line = elines.pop(0)
                if line[:1] != '>': continue
                if rje.matchExp('^>(\S+).+ gene:(\S+)\]',line): (name,gene) = rje.matchExp('^>(\S+).+ gene:(\S+)\]',line)
                else:
                    self.log.errorLog('Problem with EnsLoci line: %s' % line,printerror=False)
                    continue
                try: acc = rje.matchExp('\[acc:(\S+)',line)[0]
                except: acc = ''
                if acc: self.dict['UniEns'][acc] = gene
                self.dict['EnsLoci'][gene] = name
                self.dict['EnsDesc'][gene] = string.join(string.split(string.split(line,' [acc:')[0][1:])[1:])
                if self.opt['FullEns'] and gene not in self.list['Genes']:
                    self.list['Genes'].append(gene)
                if self.opt['FullEns'] and gene not in self.dict['GeneCard']:
                    self.dict['GeneCard'][gene] = {'EnsEMBL':gene,'Symbol':'!FAILED!'}
                self.log.printLog('\r#ENS','Parsing EnsLoci %.1f%%: %s genes' % (ex/etot,rje.integerString(len(self.dict['EnsLoci']))),newline=False,log=False)
            self.log.printLog('\r#ENS','Parsing EnsLoci complete: %s genes' % (rje.integerString(len(self.dict['EnsLoci']))))
#########################################################################################################################
    def addEnsLoci(self):   ### Adds EnsLoci data to Gene Information
        '''Adds EnsLoci data to Gene Information.'''
        if not self.dict['EnsLoci']: return
        ex = 0
        for gene in rje.sortKeys(self.dict['GeneCard']):
            if not self.dict['GeneCard'][gene].has_key('EnsEMBL') or not self.dict['GeneCard'][gene]['EnsEMBL']:
                if self.dict['GeneCard'][gene].has_key('UniProt') and self.dict['GeneCard'][gene]['UniProt'] in self.dict['UniEns']:
                    self.dict['GeneCard'][gene]['EnsEMBL'] = self.dict['UniEns'][self.dict['GeneCard'][gene]['UniProt']]
            if self.dict['GeneCard'][gene].has_key('EnsEMBL') and self.dict['GeneCard'][gene]['EnsEMBL'] in self.dict['EnsLoci']:
                ex += 1
                self.dict['GeneCard'][gene]['EnsLoci'] = self.dict['EnsLoci'][self.dict['GeneCard'][gene]['EnsEMBL']]
                self.dict['GeneCard'][gene]['EnsDesc'] = self.dict['EnsDesc'][self.dict['GeneCard'][gene]['EnsEMBL']]
            # EnsEMBL genes might be missing as they might be pseudogenes etc.
            #x#elif self.dict['GeneCard'][gene].has_key('EnsEMBL'): self.log.errorLog('EnsEMBL Gene "%s" missing from EnsLoci dictionary!' % self.dict['GeneCard'][gene]['EnsEMBL'],printerror=False)
            self.log.printLog('\r#ENS','Adding EnsLoci data: %d of %d genes' % (ex,len(self.dict['GeneCard'])),newline=False,log=False)
        self.log.printLog('\r#ENS','Added EnsLoci data for %d of %d genes' % (ex,len(self.dict['GeneCard'])))
#########################################################################################################################
    def processGenes(self,genelist):  ### Tries to extract data for genes in genelist
        '''Tries to extract data for genes in genelist.'''
        ### ~ [1] Parse data from GeneCards (or existing data) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###    
        self.deBug(self.list['Genes'])
        (gx,fx) = (0,0)
        try:
            for gene in genelist:
                if self.parseCard(gene): gx += 1
                else: fx += 1
                self.log.printLog('\r#CARD','Parsing GeneCards for %d genes: %d parsed; %d failed.' % (len(genelist),gx,fx),newline=False,log=False)
            self.log.printLog('\r#CARD','Parsing GeneCards for %d genes complete: %d parsed; %d failed.' % (len(genelist),gx,fx))
        except KeyboardInterrupt: self.log.printLog('\r#CARD','Parsing GeneCards for %d genes stopped: %d parsed; %d failed.' % (len(genelist),gx,fx))
        except: raise
        ### ~ [2] Tidy for mixed success ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        (gx,gtot,cx) = (0.0,len(self.dict['GeneCard']),0)
        for alias in rje.sortKeys(self.dict['GeneCard']):
            self.log.printLog('\r#CARD','Checking and correcting partial successes: %.1f%%' % (gx/gtot),newline=False,log=False)
            gx += 100.0
            if 'HPRD' in self.dict['GeneCard'][alias] and self.dict['GeneCard'][alias]['HPRD'] == alias:
                newalias = 'HPRD' + alias
                self.dict['GeneCard'][newalias] = self.dict['GeneCard'].pop(alias)
                alias = newalias
            try: symbol = self.dict['GeneCard'][alias]['Symbol']
            except:
                #x#print 'Fuck >> ', alias, self.dict['GeneCard'][alias], '<< Fuck!!'
                self.log.errorLog('Problem with alias "%s"' % alias)
                continue
            if symbol in self.dict['GeneCard'] and self.dict['GeneCard'][symbol]['Symbol'] == '!FAILED!':
                self.dict['GeneCard'][symbol] = self.dict['GeneCard'][alias]
                cx += 1
        self.log.printLog('\r#CARD','Checking and correcting partial successes: %d entries corrected.' % (cx))
#########################################################################################################################
    def parseCard(self,alias):  ### Parses relevant GeneCard given Alias
        '''
        Parses relevant GeneCard given Alias.
        >> alias:str = Gene symbol or alias
        '''
        try:
            ### Setup ###
            if self.dict['GeneCard'].has_key(alias) and self.dict['GeneCard'][alias]['Symbol'] != '!FAILED!':
                done = True
                for h in ['Symbol','HGNC','EnsEMBL']:
                    if h not in self.dict['GeneCard'][alias] or not self.dict['GeneCard'][alias][h]: done = False
                if alias == 'TASP1': self.deBug(done)
                if done: return True
            url = 'http://www.genecards.org/cgi-bin/carddisp.pl'
            params = 'gene=%s&alias=yes' % alias
            symbol = '!FAILED!'
            if alias in self.dict['GeneCard']:
                for d in ['Symbol','HGNC','Entrez','UniProt','EnsEMBL']:
                    if d in self.dict['GeneCard'][alias] and not self.dict['GeneCard'][alias][d]: self.dict['GeneCard'][alias].pop(d)
                if 'Symbol' not in self.dict['GeneCard'][alias]: self.dict['GeneCard'][alias]['Symbol'] = symbol
            else: self.dict['GeneCard'][alias] = {'Symbol':symbol}
            skip = not self.opt['UseWeb']
            for skipper in self.list['SkipList']:
                if rje.matchExp('(%s)' % string.replace(skipper,'*','\S+'),alias): skip = True

            ### Download GeneCard ###
            try:
                if skip: flines = []
                else: flines = urllib2.urlopen(url, params).readlines()
            except KeyboardInterrupt: raise
            except: return False

            ### Parse Data ###
            for html in flines:
                ## Primary Gene Symbol ##
                #x#if html.find(alias) > 0: self.deBug(html)
                if rje.matchExp('<title>GeneCard for (\S+)<\/title>',html):
                    symbol = rje.matchExp('<title>GeneCard for (\S+)<\/title>',html)[0]
                    if self.dict['GeneCard'].has_key(symbol) and self.dict['GeneCard'][symbol]['Symbol'] == symbol:
                        self.dict['GeneCard'][alias] = self.dict['GeneCard'][symbol]
                        return True
                    self.dict['GeneCard'][alias]['Symbol'] = symbol
                if rje.matchExp('<title>(\S+) GeneCard<\/title>',html):
                    symbol = rje.matchExp('<title>(\S+) GeneCard<\/title>',html)[0]
                    if self.dict['GeneCard'].has_key(symbol) and self.dict['GeneCard'][symbol]['Symbol'] == symbol:
                        self.dict['GeneCard'][alias] = self.dict['GeneCard'][symbol]
                        return True
                    self.dict['GeneCard'][alias]['Symbol'] = symbol
                elif rje.matchExp('<title>(\S+) Gene.+<\/title>',html):
                    symbol = rje.matchExp('<title>(\S+) Gene.+<\/title>',html)[0]
                    if self.dict['GeneCard'].has_key(symbol) and self.dict['GeneCard'][symbol]['Symbol'] == symbol:
                        self.dict['GeneCard'][alias] = self.dict['GeneCard'][symbol]
                        return True
                    self.dict['GeneCard'][alias]['Symbol'] = symbol
                ## HGNC Number ##
                if not self.dict['GeneCard'][alias].has_key('HGNC') and rje.matchExp('href="http:\/\/www.gene.ucl.ac.uk\/cgi-bin\/nomenclature\/get_data.pl\?hgnc_id=(\d+)"',html):
                    self.dict['GeneCard'][alias]['HGNC'] = rje.matchExp('href="http:\/\/www.gene.ucl.ac.uk\/cgi-bin\/nomenclature\/get_data.pl\?hgnc_id=(\d+)"',html)[0]
                if not self.dict['GeneCard'][alias].has_key('HGNC') and rje.matchExp('href="http:\/\/www.genenames.org\/data\/hgnc_data.php\?hgnc_id=(\d+)"',html):
                    self.dict['GeneCard'][alias]['HGNC'] = rje.matchExp('href="http:\/\/www.genenames.org\/data\/hgnc_data.php\?hgnc_id=(\d+)"',html)[0]
                ## Entrez Gene ##                    
                if not self.dict['GeneCard'][alias].has_key('Entrez') and rje.matchExp('href="http:\/\/www.ncbi.nlm.nih.gov\/entrez\/query.fcgi.+list_uids=(\S+)"',html):
                    self.dict['GeneCard'][alias]['Entrez'] = rje.matchExp('href="http:\/\/www.ncbi.nlm.nih.gov\/entrez\/query.fcgi.+list_uids=(\S+)"',html)[0]
                ## UniProt ##
                if not self.dict['GeneCard'][alias].has_key('UniProt') and rje.matchExp('href="http:\/\/www\.uniprot\.org\/uniprot\/(\S+)"',html):
                    self.dict['GeneCard'][alias]['UniProt'] = rje.matchExp('href="http\:\/\/www\.uniprot\.org\/uniprot\/(\S+)"',html)[0]
                ## EnsEMBL ##
                if not self.dict['GeneCard'][alias].has_key('EnsEMBL') and rje.matchExp('href="http\:\/\/www.ensembl.org\/Homo_sapiens\/geneview\?gene=(\S+)"',html):
                    self.dict['GeneCard'][alias]['EnsEMBL'] = rje.matchExp('href="http\:\/\/www.ensembl.org\/Homo_sapiens\/geneview\?gene=(\S+)"',html)[0]
                    if self.dict['GeneCard'][alias]['Symbol'] != '!FAILED!' and self.dict['GeneCard'][alias].has_key('UniProt'):
                        break
            
            ### Finishing ###
            if self.dict['GeneCard'][alias]['Symbol'] != '!FAILED!' and self.dict['GeneCard'][alias].has_key('UniProt') and self.dict['GeneCard'][alias].has_key('EnsEMBL'):
                if not self.dict['GeneCard'].has_key(symbol): self.dict['GeneCard'][symbol] = self.dict['GeneCard'][alias]
                self.deBug('%s: %s' % (alias,self.dict['GeneCard'][alias]))
                self.deBug('%s: %s' % (symbol,self.dict['GeneCard'][symbol]))
                return True
            else:
                return False
        except KeyboardInterrupt:
            self.log.errorLog('Parsing (%s) cancelled.' % alias)
            raise
        except:
            self.log.errorLog('Error in parseCard(%s)' % alias)
            return False
#########################################################################################################################
    def interactiveUpdate(self):    ### Interactive method for updating failed genes
        '''Interactive method for updating failed genes.'''
        try:
            ### ~ Setup failed lists and check interactivity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['Interactive'] < 0: return
            failures = []
            for gene in self.list['Genes']:
                if self.dict['GeneCard'][gene]['Symbol'] == '!FAILED!': failures.append(gene)
            if not failures or not rje.yesNo('Try manual mapping of %d failures?' % len(failures)): return

            ### ~ Manually map failures onto new gene list and try extracting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mymapping = {}
            newgenes = []
            try:
                for gene in failures:
                    new = rje.choice('New gene symbol for > %s <?' % gene)
                    if not new: continue
                    mymapping[gene] = new
                    if new not in newgenes: newgenes.append(new)
            except KeyboardInterrupt:
                if rje.yesNo('Quit GeneCards?',default='N'): raise
            except: raise
            self.processGenes(newgenes)
            for gene in mymapping: self.dict['GeneCard'][gene] = self.dict['GeneCard'][mymapping[gene]]
            return self.interactiveUpdate()                                             
        except: self.log.errorLog('Problem during rje_GeneCards.interactiveUpdate()')
#########################################################################################################################
    def cardMap(self):  ### Builds self.dict['CardMap'] from GeneCards read in
        '''Builds self.dict['CardMap'] from GeneCards read in.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['CardMap'] = {}
            (ax,anum) = (0.0,len(self.dict['GeneCard']))
            ambmap = []         # Ambiguous mapping terms
            ### ~ [2] Make dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for alias in self.dict['GeneCard']:
                self.log.printLog('\r#MAP','Mapping alternative IDs to Symbols: %.1f%%' % (ax/anum),newline=False,log=False)
                ax += 100.0
                symbol = self.dict['GeneCard'][alias]['Symbol']
                if symbol == '!FAILED!': symbol = alias
                for d in ['HGNC','HPRD','Entrez','UniProt','EnsEMBL']:
                    ## ~ [2a] Try to get relevant data term ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    try:
                        if d in ['HGNC','HPRD']: map = '%s:%s' % (d,self.dict['GeneCard'][alias][d])
                        else: map = self.dict['GeneCard'][alias][d]
                    except: contine
                    ## ~ [2b] Check for ambiguous mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if map in ambmap: continue
                    if map in self.dict['CardMap'] and self.dict['CardMap'][map] != symbol:
                        ambmap.append(map)
                        self.dict['CardMap'].pop(map)
                        continue
                    ## ~ [2c] Add mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.dict['CardMap'][map] = symbol
            ### ~ [3] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('\r#MAP','Mapping alternative IDs to Symbols complete: %s mappings.' % rje.integerString(len(self.dict['CardMap'])))
            return self.dict['CardMap']
        except: self.log.errorLog('Problem with GeneCards.cardMap()')            
#########################################################################################################################
    ### <4> ### Output methods                                                                                          #
#########################################################################################################################
    def outputCards(self):  ### Outputs cards to delimited file
        '''Outputs cards to delimited file.'''
        ### ~ Setup for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        genelist = self.list['Genes']
        if self.opt['Purify'] and self.opt['Restrict']:
            for gene in genelist[0:]:
                if self.dict['GeneCard'][gene]['Symbol'] not in [gene,'!FAILED!']:  # Replace with symbol
                    genelist.remove(gene)
                    if self.dict['GeneCard'][gene]['Symbol'] not in genelist: genelist.append(self.dict['GeneCard'][gene]['Symbol'])
        delimit = rje.delimitFromExt(filename=self.info['CardOut'])
        CARDOUT = open(self.info['CardOut'],'a')
        ### ~ Generate output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        (noens,noloci,ox) = (0,0,0)
        for gene in rje.sortKeys(self.dict['GeneCard']):
            if self.opt['Restrict'] and gene not in genelist: continue
            elif self.opt['Purify'] and self.dict['GeneCard'][gene]['Symbol'] not in [gene,'!FAILED!']: continue
            self.progLog('\r#OUT','Output for %s parsed genes' % rje.iStr(ox)); ox += 1
            self.dict['GeneCard'][gene]['Alias'] = gene
            self.dict['GeneCard'][gene]['Species'] = self.info['Species']
            rje.delimitedFileOutput(self,CARDOUT,self.list['Headers'],delimit,self.dict['GeneCard'][gene])
            if self.dict['GeneCard'][gene]['Symbol'] == gene:   # Not an alias
                if 'EnsEMBL' not in self.dict['GeneCard'][gene] or not self.dict['GeneCard'][gene]['EnsEMBL']: noens += 1
                if 'EnsLoci' not in self.dict['GeneCard'][gene] or not self.dict['GeneCard'][gene]['EnsLoci']: noloci += 1
        CARDOUT.close()
        self.printLog('\r#OUT','Parsed info for %d genes output to %s' % (len(self.list['Genes']),self.info['CardOut']))
        self.printLog('#ENS','%s without EnsGene; %s without EnsLoci' % (rje.integerString(noens),rje.integerString(noloci)))
#########################################################################################################################
### End of SECTION II: GeneCards Class                                                                                  #
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
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try: GeneCards(mainlog,cmd_list).run()
        
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
