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
Module:       rje_1433
Description:  14-3-3 interaction/motif analysis module
Version:      0.4
Last Edit:    14/08/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a collection of methods specifically for the analysis of 14-3-3 interaction data and motifs. In a
    sense, this is a glorified rje_misc module, containing only methods for the 14-3-3 analysis. This module was
    created for analysis of interaction data in conjunction with proteomics data from James McRedmond.

    # XTvSQ #
    The first analysis is an X-Tandem vs Sequest analysis. The goals here are as follows:
    * Map IPI onto single genes, identifying synonyms.
    * Produce a table with links for all IPI proteins in dataset to the genes used in the analysis.
    * Produce a table of results compressed into gene identifiers with best scores
    * Compare the overlap for each isoform (including combined "1433" isoform) at differenct cut-offs.

    # XTvSQvPPI #
    The analysis is similar to XTvSQ but analyses the overlap of the combined identifications with known interactions.

    # Pingu #
    This calls the Pingu application to make files from MS and Literature samples.

    # SLiMFinder #
    This generates a number of sequence datasets by processing the MS and Literature data in different ways ready for
    SLiMFinder motif discovery.

    # SLiMSearch #
    This finds known 14-3-3 motifs in the data.

Commandline:
    ### ~ Main analysis controlling parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    analyses=LIST   : List of analyses to perform []  (XTvSQ,Pingu,XTvSQvPPI,SLiMFinder,SLiMSearch)
    msfile=FILE     : File containing MS Data [1433_msdata_full.txt]
    litfile=FILE    : File of literature interactions [1433_Lit.tdt]
    xtcap=X         : XTandem score cap ("Best" X-Tandem score) [-10]
    xtcut=X         : XTandem score cutoff. Scores must be less than or equal to this [-2.5]
    sqcut=X         : Sequest score cutoff. Scores must be greater than or equal to this [0.95]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_genecards, rje_uniprot
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, math, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, pingu, slimsearch
import rje_genecards, rje_hprd, rje_ppi, rje_uniprot, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added use of Pingu.
    # 0.2 - Added literature.
    # 0.3 - Added SLiMFinder dataset generation.
    # 0.4 - Added SLiMSearch for known 14-3-3 motifs.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Incorporate HPRD Data
    # [Y] : Output Files for Cytoscape - use a common data structure and rje_ppi
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_1433', '0.4', 'August 2007','2007')
    description = '14-3-3 interaction/motif analysis module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Data Protection issues: Please do not use this without permission.',rje_zen.Zen(None,[]).wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
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
### Constants
iso1433 = {'beta':'YWHAB','eta':'YWHAH','epsilon':'YWHAE','gamma':'YWHAG','sigma':'SFN','theta':'YWHAQ','zeta':'YWHAZ'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: RJE_1433 Class                                                                                          #
#########################################################################################################################
class RJE_1433(rje.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2007).

    Info:str
    - LitFile = File of literature interactions [1433_Lit.tdt]
    - MSFile = File containing MS Data [1433_msdata_full.txt]
    
    Opt:boolean

    Stat:numeric
    - XTCap = XTandem score cap ("Best" X-Tandem score) [-10]
    - XTCut = XTandem score cutoff. Scores must be less than or equal to this [-2.5]
    - SQCut = Sequest score cutoff. Scores must be greater than or equal to this [0.95]

    List:list
    - Analyses = List of analyses to perform [XTvSQ]

    Dict:dictionary
    - IPI = {Gene:[IPIs]}
    - IPI Links = dictionary of {IPI:External links}
    - MSData = Dictionary of {Gene:[Presence]}
    - Synonyms = Dictionary of {IPI:[Synonyms]}
    - MSFull = Dictionary of {Iso:{Gene:{'XT':best score,'SQ':best score}}}

    Obj:RJE_Objects
    - GeneCards = rje_genecards.GeneCards object
    - Pingu = pingu.PINGU Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['LitFile','MSFile']
        self.optlist = []
        self.statlist = ['XTCap','XTCut','SQCut']
        self.listlist = ['Analyses']
        self.dictlist = ['IPI','MSData','IPI Links','Synonyms','MSFull']
        self.objlist = ['GeneCards','Pingu']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'LitFile':'1433_Lit.tdt','MSFile':'1433_msdata_full.txt'})
        self.setStat({'XTCap':-10,'XTCut':-2.5,'SQCut':0.95})
        self.list['Analyses'] = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)       ### General Options ### 
                self._cmdReadList(cmd,'file',['LitFile','MSFile'])
                self._cmdReadList(cmd,'stat',['XTCap','XTCut','SQCut'])  
                self._cmdReadList(cmd,'list',['Analyses'])  
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self.list['Analyses'] = string.split(string.join(self.list['Analyses']).lower())
#########################################################################################################################
    ### <3> ### Main Workflow run methods                                                                               #
#########################################################################################################################
    def run(self):  ### Main Run method for RJE_1433
        '''Main Run method for RJE_1433. For now, most of this is hard-coded. Options may be added later.'''
        try:
            ###~ [1] Initial analysis of X-Tandem vs Sequest data to look for overlap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'XTvSQ'.lower() in self.list['Analyses']:
                self.parseIPI()         # Creates self.dict['IPI Links']
                self.parseFullMSData()  # Parses and compacts MS Data into self.dict['MSFull']
                for iso in self.dict['MSFull']: self.compareXTvSQ(iso)
                
            ###~ [2] Basic Pingu analysis of data. Will need to use cut-offs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'pingu' in self.list['Analyses']: self.pingu()
            elif 'slimfinder' in self.list['Analyses'] or 'slimsearch' in self.list['Analyses']:
                self.cmd_list = ['fulloutput=F','allbyall=0'] + self.cmd_list
                self.pingu()            
            
            ###~ [3] XTandem vs SeQuest vs PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'XTvSQvPPI'.lower() in self.list['Analyses']: self.compareXTvSQvPPI()

            ###~ [4] SLiMFinder dataset generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'slimfinder' in self.list['Analyses']:  self.slimFinder()

            ###~ [4] SLiMFinder dataset generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'slimsearch' in self.list['Analyses']:  self.slimSearch()

            ###~ [X] Old analysis of Proteomics data: mapping HPRD and combining experiments ~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'old' in self.list['Analyses']:
                self.parseIPI()     # Creates self.dict['IPI Links']
                self.parseMSData()  # Generates self.dict['MSData'] and self.dict['IPI']
                self.synonyms()     # Generates self.dict['Synonyms'] and synonyms file
                self.msGeneCards()  # Sets up self.obj['GeneCards'] and outputs file
            
        except:
            self.log.errorLog('Major problem with rje_1433.run()')
#########################################################################################################################
    ### <4> ### General MS Data and links methods                                                                       #
#########################################################################################################################
    def parseIPI(self):     ### Parses IPI links for combining data. Creates self.dict['IPI Links'].
        '''Parses IPI links for combining data. Creates self.dict['IPI Links'].'''
        try:
            ###~[1]~Parsing DB Links from IPI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ipi_dat = 'ipi.HUMAN.v3.16.dat'
            ilinkfile = 'ipi_links.tdt'
            ### Make the file if it does not exist ###
            if not os.path.exists(ilinkfile):
                ## Read Data into UniProtEntry objects ##
                ipi = rje_uniprot.UniProt(self.log,self.cmd_list+['uniprot=%s' % ipi_dat])
                ipi.readUniProt()
                ## Make File ##
                headers = ['IPI','Symbol','HGNC','Entrez','EnsG','EnsP','UniID','UniAcc','UniSV','RefSeq']
                rje.delimitedFileOutput(self,ilinkfile,headers,'\t')
                ox = 0.0
                for entry in ipi.list['Entry']:
                    for db in entry.dict['DB']: entry.dict['DB'][db] = string.join(entry.dict['DB'][db],',')
                    entry.dict['DB']['IPI'] = entry.info['ID']
                    rje.delimitedFileOutput(self,ilinkfile,headers,'\t',entry.dict['DB'])
                    ox += 100.0
                    self.log.printLog('\r#IPI','Making %s %.1f%%' % (ilinkfile,ox/ipi.entryNum()),newline=False,log=False)
                self.log.printLog('\r#IPI','Making %s complete: %s entries.' % (ilinkfile,rje.integerString(ipi.entryNum())))
                ## Cleanup ##
                del ipi
            ### Read file into data dictionary ###
            self.dict['IPI Links'] = rje.dataDict(self,ilinkfile)  # This is now a dictionary of IPI:Data
        except:
            self.log.errorLog('Error in RJE_1433.parseIPI()',printerror=True,quitchoice=False)
#########################################################################################################################
    def parseFullMSData(self):  ### Parses and compacts MS Data into self.dict['MSFull']
        '''Parses and compacts MS Data into self.dict['MSFull'] = {Iso:{Gene:{'XT':best score,'SQ':best score}}}'''
        try:
            ### ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ipi_links = self.dict['IPI Links']
            msfile = '1433_msdata_full.txt'
            msfull = {'1433':{}}

            ### ~ Read in Full MS Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in self.loadFromFile(msfile,chomplines=True):
                ## Get data for given line of MS data ##
                data = string.split(line,'\t')
                if data[0] == 'Experiment': continue
                try: [iso,alg,ipi,desc,score] = data
                except: continue
                if iso in ['Vector','Cell']: iso = 'Control'
                if iso not in msfull: msfull[iso] = {}
                score = float(score)
                if alg == 'XT':
                    if self.stat['XTCap'] < 0 and score < self.stat['XTCap']: score = self.stat['XTCap']
                    score = -score
                ## Convert to gene(s) from IPI ##
                if ipi[:4] == 'IPI:': ipi = ipi[4:]
                elif ipi[:3] != 'IPI':
                    self.log.printLog('#REJ','Rejected non-IPI %s protein "%s"' % (iso,ipi))
                    continue
                genelist = [ipi]
                for db in ['Symbol','EnsG','Entrez','UniAcc','RefSeq']:
                    if ipi_links[ipi][db]:
                        genelist = string.split(ipi_links[ipi][db],',')
                        break
                ## Update dictionaries ##
                for gene in genelist:
                    ## Gene-IPI links ##
                    if self.dict['IPI'].has_key(gene) and ipi not in self.dict['IPI'][gene]: self.dict['IPI'][gene].append(ipi)
                    elif not self.dict['IPI'].has_key(gene): self.dict['IPI'][gene] = [ipi]
                    ## MS Results ##
                    if not gene in msfull[iso]: msfull[iso][gene] = {'XT':0.0,'SQ':0.0}
                    if iso != 'Control' and gene not in msfull['1433']: msfull['1433'][gene] = {'XT':0.0,'SQ':0.0}
                    if score > msfull[iso][gene][alg]: msfull[iso][gene][alg] = score
                    if iso != 'Control' and score > msfull['1433'][gene][alg]: msfull['1433'][gene][alg] = score
            self.log.printLog('#MS','MS Data read and combined')

            ### ~ Finish and output reduced data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## IPI Mapping ##
            MAPIPI = open('1433_ipi_map.tdt','w')
            MAPIPI.write('IPI\tGene\n')
            for gene in self.dict['IPI']:
                for ipi in self.dict['IPI'][gene]:
                    MAPIPI.write('%s\t%s\n' % (gene,ipi))
            MAPIPI.close()
            self.log.printLog('#MAP','IPI mapping for %s genes output to 1433_ipi_map.tdt' % rje.integerString(len(self.dict['IPI'])))
            ## Combined MS Data ##
            newfile = '1433_msdata_combined.tdt'
            headers = ['Isoform','Gene','XT','SQ']
            rje.delimitedFileOutput(self,newfile,headers,'\t',{},True)
            for iso in msfull:
                for gene in msfull[iso]:
                    data = {'Isoform':iso,'Gene':gene,'XT':msfull[iso][gene]['XT'],'SQ':msfull[iso][gene]['SQ']}
                    rje.delimitedFileOutput(self,newfile,headers,'\t',data)
            self.dict['MSFull'] = msfull
            self.log.printLog('#MS','MS combined data output to %s' % newfile)
        except: self.log.errorLog('Bad bummer during rje_1433.parseFullMSData()')
#########################################################################################################################
    def parseMSData(self):  ### Generates self.dict['MSData'] and self.dict['IPI']
        '''Generates self.dict['MSData'] and self.dict['IPI'].'''
        try:
            ###~[1]~Parsing MS Data from Text File~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ipi_links = self.dict['IPI Links']
            experiments = ['ctl1','ctl2','beta','eta','epsilon','gamma','sigma','theta','zeta']
            msfile = '1433_msdata.tdt'
            ### Read from File ###
            ms = rje.dataDict(self,msfile,mainkeys=['IPI','Method'],datakeys=['Present'],delimit='\t')
            ### Make into self.dict['MSData'] and self.dict['IPI'] ###
            for ident in ms:
                (ipi,method) = string.split(ident,'\t')
                score = {'Sequest':1,'Xtandem':2}[method]
                present = []
                for i in range(len(ms[ident]['Present'])):
                    if ms[ident]['Present'][i] == '0': present.append(0)
                    else: present.append(score)
                if len(present) != len(experiments):
                    self.log.errorLog('%s only %d presence calls' % (ident,len(present)),printerror=False)
                    continue
                ## Get list of genes ##
                genelist = [ipi]
                for db in ['Symbol','EnsG','Entrez','UniAcc','RefSeq']:
                    if ipi_links[ipi][db]:
                        genelist = string.split(ipi_links[ipi][db],',')
                        break
                ## Update dictionaries ##
                for gene in genelist:
                    ## Gene-IPI links ##
                    if self.dict['IPI'].has_key(gene) and ipi not in self.dict['IPI'][gene]: self.dict['IPI'][gene].append(ipi)
                    elif not self.dict['IPI'].has_key(gene): self.dict['IPI'][gene] = [ipi]
                    ## MS Results ##
                    if not self.dict['MSData'].has_key(gene):
                        self.dict['MSData'][gene] = [0] * len(experiments)
                    for i in range(len(experiments)):
                        if self.dict['MSData'][gene][i] in [score,3]: continue  # Already found
                        elif present[i]: self.dict['MSData'][gene][i] += score
        except:
            self.log.errorLog('Error in RJE_1433.parseMSData()',printerror=True,quitchoice=False)
#########################################################################################################################
    def synonyms(self):     ### Generates self.dict['Synonyms'] and synonyms file
        '''Generates self.dict['Synonyms'] and synonyms file.'''
        try:
            ###~[3]~Synonyms and Multiplicity File~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ### Output self.dict['IPI'] into synonym file! ###
            mfile = '1433.synonyms.txt'
            MFILE = open(mfile,'w')
            MFILE.write('# Genes/Proteins with multiple IPI IDs in dataset #\n')
            temp_ipi = {}   # Dictionary of IPI mappings to genes
            ### Genes with multiple IPI ###
            for gene in self.dict['IPI']:
                if len(self.dict['IPI'][gene]) > 1: MFILE.write('%s = %s\n' % (gene,string.join(self.dict['IPI'][gene],', ')))
                for ipi in self.dict['IPI'][gene]:
                    if ipi in temp_ipi: temp_ipi[ipi].append(gene)
                    else: temp_ipi[ipi] = [gene]
            MFILE.write('\n# IPI IDs in dataset with multiple Genes/Proteins #\n')
            ### IPI with multiple mapping ###
            for ipi in temp_ipi:
                if len(temp_ipi[ipi]) > 1: MFILE.write('%s = %s\n' % (ipi,string.join(temp_ipi[ipi],', ')))
            ### Make Synonyms dictionary ###
            self.dict['Synonym'] = {}
            ## First pass to fill basics ##
            for ipi in temp_ipi:
                for g1 in temp_ipi[ipi][:-1]:
                    if g1 not in self.dict['Synonym']: self.dict['Synonym'][g1] = []
                    for g2 in temp_ipi[ipi][1:]:
                        if g1 == g2: continue
                        if g2 not in self.dict['Synonym']: self.dict['Synonym'][g2] = []
                        if g2 not in self.dict['Synonym'][g1]: self.dict['Synonym'][g1].append(g2)
                        if g1 not in self.dict['Synonym'][g2]: self.dict['Synonym'][g2].append(g1)
            ## Second pass to expand to make all ##
            for gene in self.dict['Synonym']:
                gx = 0
                while gx != len(self.dict['Synonym'][gene]):
                    gx = len(self.dict['Synonym'][gene])
                    for syn in self.dict['Synonym'][gene][0:]:
                        for g2 in self.dict['Synonym'][syn]:
                            if g2 == gene: continue
                            if g2 not in self.dict['Synonym'][gene]: self.dict['Synonym'][gene].append(g2)
            ## Output synonyms ##
            MFILE.write('\n# Gene/Protein Synonyms from IPI Linkage #\n')
            for gene in self.dict['Synonym']:
                if len(self.dict['Synonym'][gene]) > 0: MFILE.write('%s = %s\n' % (gene,string.join(self.dict['Synonym'][gene],', ')))
            MFILE.close()
            print open(mfile,'r').read()
        except:
            self.log.errorLog('Error in RJE_1433.parseMSData()',printerror=True,quitchoice=False)
#########################################################################################################################
    def msGeneCards(self):  ### Sets up self.obj['GeneCards'] and outputs file
        '''Sets up self.obj['GeneCards'] and outputs file.'''
        try:
            ###~[4]~Output GeneCards Data File~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            gfile = '1433.ms_interactions.tdt'
            altsource = '/home/richard/TestBed/GeneCards/EnsLoci-HPRD-Full-2.tdt'
            gcmd = ['cardout=%s' % gfile,'ensloci=/home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas','altsource=%s' % altsource,'pureout=T']
            cards = rje_genecards.GeneCards(self.log,gcmd+self.cmd_list)
            experiments = ['ctl1','ctl2','beta','eta','epsilon','gamma','sigma','theta','zeta']
            if not os.path.exists(gfile):
                ### Output MS Data ###
                headers = ['Alias','Species','Symbol','HGNC','Entrez','UniProt','EnsEMBL','EnsLoci','IPI']
                for ms in experiments: headers.append('MS_%s%s' % (ms.upper()[:1],ms[1:]))
                rje.delimitedFileOutput(self,gfile,headers,'\t',rje_backup=True)
                for gene in self.dict['MSData']:
                    datadict = {'Alias':gene,'IPI':string.join(self.dict['IPI'][gene],',')}
                    for ms in experiments:
                        datadict['MS_%s%s' % (ms.upper()[:1],ms[1:])] = self.dict['MSData'][gene][experiments.index(ms)]
                    rje.delimitedFileOutput(self,gfile,headers,'\t',datadict)
                ### Run RJE_GENECARDS to fill in extra data! ###
                cards.list['Genes'] = rje.sortKeys(self.dict['MSData'])
                cards.run()
            else:
                cards.setup()
                cards.list['Genes'] = rje.sortKeys(cards.dict['GeneCard'])
            self.obj['GeneCards'] = cards
        except:
            self.log.errorLog('Error in RJE_1433.msGeneCards()',printerror=True,quitchoice=False)
#########################################################################################################################
    ### <5> ### Data comparison/analysis methods                                                                        #
#########################################################################################################################
    def compareXTvSQ(self,iso='1433'): ### Generate tables of XT & SQ cut-offs versus overlap
        '''Generate tables of XT & SQ cut-offs versus overlap.'''
        try:
            ### ~ [1] Setup lists for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            genes = rje.sortKeys(self.dict['MSFull'][iso])
            ## Make a sorted list of XT and SQ scores ##
            xgenes = {}
            sgenes = {}
            for gene in genes:
                x = self.dict['MSFull'][iso][gene]['XT']
                #if x > 0: x = 0.1 + (int((100.0 * math.log(x,10))+0.5) / 100.0)
                #self.dict['MSFull'][iso][gene]['XT'] = x
                xgenes[x] = []
                s = self.dict['MSFull'][iso][gene]['SQ']
                #s = int((100.0 * self.dict['MSFull'][iso][gene]['SQ'])+0.5) / 100.0
                #self.dict['MSFull'][iso][gene]['SQ'] = s
                sgenes[s] = []
            self.log.printLog('#XT','%s XTandem scores for %s' % (rje.integerString(len(xgenes)),iso))
            self.log.printLog('#XT','%s Sequest scores for %s' % (rje.integerString(len(sgenes)),iso))
            ## Make a dictionary of score:genes ##
            for x in rje.sortKeys(xgenes):
                for gene in genes:
                    if self.dict['MSFull'][iso][gene]['XT'] >= x: xgenes[x].append(gene)
            for s in rje.sortKeys(sgenes):
                for gene in genes:
                    if self.dict['MSFull'][iso][gene]['SQ'] >= s: sgenes[s].append(gene)
                    
            ### ~ Perform comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cfile = '1433.%s.XTvSQ.tdt' % iso
            headers = ['XT','SQ','Overlap']
            if os.path.exists(cfile): os.unlink(cfile)
            rje.delimitedFileOutput(self,cfile,headers,'\t')
            xgenes.pop(0.0)
            sgenes.pop(0.0)
            (c,cx) = (len(xgenes)*len(sgenes),0.0)
            for x in rje.sortKeys(xgenes):
                for s in rje.sortKeys(sgenes):
                    self.log.printLog('\r#XvS','Comparing overlap for %s: %.1f%%' % (iso,(cx/c)),newline=False,log=False)
                    cx += 100.0
                    ox = 0
                    gx = len(sgenes[s])
                    for gene in xgenes[x]:
                        if gene in sgenes[s]: ox += 1
                        else: gx += 1
                    px = float(ox) / float(gx)
                    rje.delimitedFileOutput(self,cfile,headers,'\t',{'XT':x,'SQ':s,'Overlap':px})
            self.log.printLog('\r#XvS','Comparing overlap for %s complete (%s points)' % (iso,rje.integerString(c)))                
        except: self.log.errorLog('Nooo! Problem with rje_1433.compareXTvSQ(%s)' % iso)
#########################################################################################################################
    def readMSData(self):   ### Reads data into dictionary, using cut-offs, and returns
        '''Reads data into dictionary, using cut-offs, and returns.'''
        try:
            ### ~ [1] Read in Full MS Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            msfile = self.info['MSFile']
            msdata = {'1433':{}}
            for line in self.loadFromFile(msfile,chomplines=True):
                ## Get data for given line of MS data ##
                data = string.split(line,'\t')
                if data[0] == 'Experiment': continue
                try: [iso,alg,ipi,desc,score] = data
                except: continue
                if ipi[:4] == 'IPI:': ipi = ipi[4:]
                if ipi[:3] != 'IPI':
                    try: ipi = rje.matchExp('^\S+_HUMAN \((\S+)\)',ipi)[0]
                    except: continue
                if iso in ['Vector','Cell']: iso = 'Control'
                if iso not in msdata: msdata[iso] = {}                
                if ipi not in msdata[iso]: msdata[iso][ipi] = {'XT':0.0,'SQ':0.0}
                msdata[iso][ipi][alg] = float(score)
                if iso == 'Control': continue
                if ipi not in msdata['1433']: msdata['1433'][ipi] = {'XT':0.0,'SQ':0.0}
                if alg == 'XT': msdata['1433'][ipi][alg] = min(msdata['1433'][ipi][alg],float(score))
                else: msdata['1433'][ipi][alg] = max(msdata['1433'][ipi][alg],float(score))
            self.log.printLog('#IN','MS data read from %s' % msfile)
            return msdata
        except: self.log.errorLog(rje_zen.Zen().wisdom())                        
#########################################################################################################################
    def pingu(self):    ### Generates input for Pingu and then runs
        '''Generates input for Pingu and then runs.'''
        try:
            ### ~ [1] Generate samples input for Pingu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfile = '1433_jm.tdt'
            headers = ['Sample','Identifier','XT','SQ']
            ## ~ Read in Full MS Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            msdata = self.readMSData()
            ## ~ Output in Pingu-recognisable format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.delimitedFileOutput(self,sfile,headers,rje_backup=True)
            for iso in msdata:
                for ipi in msdata[iso].keys()[0:]:
                    if msdata[iso][ipi]['XT'] > self.stat['XTCut'] and msdata[iso][ipi]['SQ'] < self.stat['SQCut']: continue
                    data = msdata[iso].pop(ipi)
                    data['Sample'] = iso
                    data['Identifier'] = ipi
                    rje.delimitedFileOutput(self,sfile,headers,datadict=data)
            self.log.printLog('#OUT','Sample data output to %s' % sfile)

            ### ~ [2] Add Literature PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['LitFile'] and os.path.exists(self.info['LitFile']):
                lit = rje.dataDict(self,self.info['LitFile'],['Paper','Bait'],['Paper','Bait','Identification'],lists=True)
                for entry in lit.values():
                    sample = string.join(entry['Bait']+entry['Paper'],'-')
                    for gene in entry['Identification']:
                        data['Sample'] = sample
                        data['Identifier'] = gene
                        rje.delimitedFileOutput(self,sfile,headers,datadict=data)
                self.log.printLog('#LIT','%d literature samples added to %s' % (len(lit),sfile))
                    
            ### ~ [3] Run Pingu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pcmd = ['fulloutput=T','baits=YWHAB,YWHAH,YWHAE,YWHAG,SFN,YWHAQ,YWHAZ','allbyall=3','addbaits=T',
                    'controls=Control','experiments=Beta,Eta,Epsilon,Gamma,Sigma,Theta,Zeta']
            pcmd = pcmd + self.cmd_list + ['data=%s' % sfile]
            self.obj['Pingu'] = pingu.PINGU(self.log,pcmd)
            self.obj['Pingu'].run()                   
        except: self.log.errorLog(rje_zen.Zen().wisdom())                        
#########################################################################################################################
    def compareXTvSQvPPI(self): ### Compares XTandem and Sequest for overlap with known 1433 interactors
        '''Compares XTandem and Sequest for overlap with known 1433 interactors.'''
        try:
            ### ~ [1] Generate samples input for Pingu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfile = '1433_all.tdt'
            headers = ['Sample','Identifier','XT','SQ']
            ## ~ Read in Full MS Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            msdata = self.readMSData()
            ## ~ Output in Pingu-recognisable format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.delimitedFileOutput(self,sfile,headers,rje_backup=True)
            for iso in msdata:
                for ipi in msdata[iso].keys()[0:]:
                    data = msdata[iso].pop(ipi)
                    data['Sample'] = iso
                    data['Identifier'] = ipi
                    rje.delimitedFileOutput(self,sfile,headers,datadict=data)
            self.log.printLog('#OUT','Sample data output to %s' % sfile)

            ### ~ [2] Run Pingu to combine data and parse PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pcmd = ['fulloutput=F','allbyall=0'] + self.cmd_list + ['data=%s' % sfile]
            self.obj['PINGU'] = pingu.PINGU(self.log,pcmd)
            self.obj['PINGU'].run()
            
            ### ~ [3] Make XT vs SQ vs PPI comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            msdata = self.readMSData()
            for iso in msdata.keys():
                ## ~ Make a sorted list of XT and SQ scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                xgenes = {}     # XTandem genes and scores
                sgenes = {}     # Sequest genes and scores
                for ipi in msdata[iso]:
                    gene = self.obj['PINGU'].geneMap(ipi)
                    x = msdata[iso][ipi]['XT']
                    if self.stat['XTCap'] < 0: x = max(x,self.stat['XTCap'])
                    xgenes[x] = []
                    s = msdata[iso][ipi]['SQ']
                    sgenes[s] = []
                self.log.printLog('#XT','%s XTandem scores for %s' % (rje.integerString(len(xgenes)),iso))
                self.log.printLog('#SQ','%s Sequest scores for %s' % (rje.integerString(len(sgenes)),iso))
                ## Make a dictionary of score:genes ##
                xgenes.pop(0.0)
                sgenes.pop(0.0)
                for x in rje.sortKeys(xgenes):
                    for ipi in msdata[iso]:
                        gene = self.obj['PINGU'].geneMap(ipi)
                        if msdata[iso][ipi]['XT'] <= x: xgenes[x].append(gene)
                for s in rje.sortKeys(sgenes):
                    for ipi in msdata[iso]:
                        gene = self.obj['PINGU'].geneMap(ipi)
                        if msdata[iso][ipi]['SQ'] >= s: sgenes[s].append(gene)
                        
                ## ~ Setup comparison files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ifile = '1433.%s.XTvSQvIso_PPI.tdt' % iso
                afile = '1433.%s.XTvSQvAll_PPI.tdt' % iso
                headers = ['XT','SQ','Overlap']
                for cfile in [afile,ifile]:
                    if os.path.exists(cfile): os.unlink(cfile)
                    rje.delimitedFileOutput(self,cfile,headers,'\t')

                ## ~ Perform comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if iso.lower() in iso1433: igene = iso1433[iso.lower()]
                else: igene = None
                print igene
                (c,cx) = (len(xgenes)*len(sgenes),0.0)
                for x in rje.sortKeys(xgenes):
                    for s in rje.sortKeys(sgenes):
                        self.log.printLog('\r#XvSvP','Comparing PPI overlap for %s: %.1f%%' % (iso,(cx/c)),newline=False,log=False)
                        cx += 100.0
                        (ix,ax) = (0,0)
                        cgenes = rje.sortUnique(sgenes[s]+xgenes[x])    # Combined gene dataset
                        gx = len(cgenes)
                        for gene in cgenes:
                            if igene and gene in rje.dictValues(self.obj['PINGU'].dict['PPI'],igene): ix += 1
                            for agene in iso1433.values():   # Any
                                if gene in rje.dictValues(self.obj['PINGU'].dict['PPI'],agene):
                                    ax += 1
                                    break
                        ipx = float(ix) / float(gx)
                        rje.delimitedFileOutput(self,ifile,headers,'\t',{'XT':-x,'SQ':s,'Overlap':ipx})
                        apx = float(ax) / float(gx)
                        rje.delimitedFileOutput(self,afile,headers,'\t',{'XT':-x,'SQ':s,'Overlap':apx})
                self.log.printLog('\r#XvSvP','Comparing PPI overlap for %s complete (%s points)' % (iso,rje.integerString(c)))                

        except: self.log.errorLog(rje_zen.Zen().wisdom())                        
#########################################################################################################################
    def slimFinder(self):   ### Generates SLiMFinder datasets using Pingu data
        '''Generates SLiMFinder datasets using Pingu data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ping = self.obj['Pingu']
            if not ping.obj['GeneCards']: return self.log.errorLog('Cannot map EnsLoci without GeneCards.', printerror=False)
            #ensloci = ping.getEnsLoci()
            #seqdict = ensloci.seqNameDic()
            #genecards = ping.obj['GeneCards'].dict['GeneCard']
            #if not seqdict: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)

            ### ~ [2] Generate lists of genes for datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ### This section needs to generate the following dataset types:
            ### 1. James MS samples - full (& combined)     =JM
            ### 2. James MS samples - unique isoform interactors only =JMU
            ### 3. The above without any Control genes  =NoC
            ### 4. 1&2 for literature 14-3-3 PPI =Lit
            ### 5. 1&2&3 for combined MS+Lit PPI =Comb
            ## ~ [2a] Make the JM & Lit and Comb lists for each isoform plus '1433' (All) ~~~~~~~~~ ##
            controls = ping.dict['Datasets'].pop('Control')
            ping.dict['Datasets'].pop('1433')
            ping.dict['Datasets']['1433-MS'] = []
            ping.dict['Datasets']['1433-Lit'] = []
            for sample in ping.dict['Datasets'].keys()[0:]:
                if sample.lower() in iso1433:
                    newsample = '%s-MS' % iso1433[sample.lower()]
                    nocontrol = '%s-MSNoC' % iso1433[sample.lower()]
                    ping.dict['Datasets'][nocontrol] = []
                    for gene in ping.dict['Datasets'][sample]:
                        if gene not in ping.dict['Datasets']['1433-MS']: ping.dict['Datasets']['1433-MS'].append(gene)
                        if gene not in controls: ping.dict['Datasets'][nocontrol].append(gene)
                    ping.dict['Datasets'][newsample] = ping.dict['Datasets'].pop(sample)
                else:
                    for gene in ping.dict['Datasets'][sample]:
                        if gene not in ping.dict['Datasets']['1433-Lit']: ping.dict['Datasets']['1433-Lit'].append(gene)
            print rje.sortKeys(ping.dict['Datasets'])

            ## ~ [2b] Make the Comb lists for each isoform plus '1433' (All) ~~~~~~~~~~~~~~~~~~~~~~ ##
            for iso in iso1433.values() + ['1433']: ping.dict['Datasets']['%s-Comb' % iso] = []
            for sample in ping.dict['Datasets']:
                (iso,type) = string.split(sample,'-')
                if type == 'Comb': continue
                for gene in ping.dict['Datasets'][sample]:
                    if gene not in ping.dict['Datasets']['%s-Comb' % iso]: ping.dict['Datasets']['%s-Comb' % iso].append(gene)
            print rje.sortKeys(ping.dict['Datasets'])

            ## ~ [2c] Generate Unique Isoform Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            redundants = rje.sortKeys(ping.dict['Datasets'])
            for sample in redundants:
                (iso,type) = string.split(sample,'-')
                if iso == '1433': continue
                ping.dict['Datasets']['%s-Unique' % sample] = ping.dict['Datasets'][sample][0:]
                for sample2 in redundants:
                    try: (iso2,type2) = string.split(sample2,'-')
                    except: continue
                    #x#print iso, iso2, type, type2, type2 != type or iso2 == iso
                    if type2 != type or iso2 in [iso,'1433']: continue
                    p1 = len(ping.dict['Datasets']['%s-Unique' % sample])
                    for gene in ping.dict['Datasets']['%s-Unique' % sample][0:]:
                        if gene in ping.dict['Datasets'][sample2]: ping.dict['Datasets']['%s-Unique' % sample].remove(gene)
                    self.log.printLog('#UNIQ','%s %s: Unique %d -> %s %s -> %d' % (iso,type,p1,iso2,type2,len(ping.dict['Datasets']['%s-Unique' % sample])))
            print rje.sortKeys(ping.dict['Datasets'])

            ### ~ [3] Generate Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ping.sampleSeqFiles()
        
        except: self.log.errorLog(rje_zen.Zen().wisdom())                        
#########################################################################################################################
    def slimSearch(self):   ### Search proteins with known motifs and generate cytoscape files
        '''Search proteins with known motifs and generate cytoscape files.'''
        try:### ~ [1] Make Gene File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['Pingu'].geneSeqFile('1433_jm.genes.fas')
            
            ### ~ [2] Perform SLiMSearch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = ['motifs=1433_rje.motifs'] + self.cmd_list + ['seqin=1433_jm.genes.fas','gnspacc=F']
            search = slimsearch.SLiMSearch(self.log,scmd)
            search.run()

            ### ~ [3] Convert to Cytoscape Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for slim in search.slims():
                print slim
                self.deBug(slim.dict)
                genelist = self.obj['Pingu'].list['Genes'][0:]
                NODE = open('1433_jm.%s.noa' % (slim.info['Name']),'w')
                NODE.write('%s\n' % slim.info['Name'])
                for seq in slim.dict['Occ']:
                    NODE.write('%s = %d\n' % (seq.shortName(),len(slim.dict['Occ'][seq])))
                    if seq.shortName() in genelist: genelist.remove(seq.shortName())
                #x#for gene in genelist: NODE.write('%s = -1\n' % (gene))
                NODE.close()
                self.log.printLog('#NOA','Node attributes output to 1433_jm.%s.noa' % (slim.info['Name']))
            
        except: self.log.errorLog(rje_zen.Zen().wisdom())                        
#########################################################################################################################
    ### <X> ### Interaction database/literature links methods                                                           #
#########################################################################################################################
    def combineHPRD(self):  #!# Old HPRD combine and output method
        try:
            ###~[5]~Combine with HPRD Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            cards = self.obj['GeneCards']
            hcmd = ['hprdpath=/home/richard/Databases/HPRD/HPRD_Release_6_01012007/FLAT_FILES/']
            hprd = rje_hprd.HPRD(self.log,hcmd+self.cmd_list)
            hprd.parse()
            hprd.addToGeneCards(cards,addcards=False)
            cards.info['CardOut'] = '1433.ms_hprdcards.tdt'
            delimit = rje.delimitFromExt(filename=cards.info['CardOut'])
            rje.delimitedFileOutput(self,cards.info['CardOut'],cards.list['Headers'],delimit,rje_backup=True)
            cards.outputCards()  ### Outputs cards to delimited file
            ## Account for EnsEMBL mapped onto GeneCards ##
            for ms in rje.sortKeys(self.dict['MSData']):
                gene = cards.dict['GeneCard'][ms]['Symbol']
                if gene not in ['!FAILED!',ms]:
                    if gene in self.dict['MSData']:
                        self.log.errorLog('Unexpected duplicity! %s & %s' % (ms,gene),printerror=False)
                        print self.dict['MSData'][ms], self.dict['MSData'][gene]
                    self.dict['MSData'][gene] = self.dict['MSData'].pop(ms)
                    
            ###~[7]~Output Files for Cytoscape~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            #x#self.OLDCytoscape(hprd,cards)    #!# experiments
            # self.dict['MSData'] contains the best identifier for each MS protein now. #
            # cards.dict['GeneCard'] contains links information for all identifiers
            # self.dict['Synonym'] contains synonymous identifiers
            # hprd.dict['HPRD'] contains links information for lit IDs, mapped using hprd.dict['Mapping'][gene] or cards.dict['GeneCard'][gene]['HPRD']

            ### Generate self.dict['PPI'], self.dict['Node'] and self.dict['Edge'] ###
            # self.dict['PPI'] = ppi dictionary of {gene:{gene2:[types]}}   : Initial asymmetric
            # self.dict['Node'] = {Type:{gene:annotation}}
            # self.dict['Edge'] = {Type:{gene1:{intType:{gene2:annotation}}}}
            self.dict['PPI'] = {}

            ### Generate lists of 1433 interactors ###
            map1433 = {'ctl1':'CTRL','ctl2':'CTRL','beta':'YWHAB','eta':'YWHAH','epsilon':'YWHAE','gamma':'YWHAG','sigma':'SFN','theta':'YWHAQ','zeta':'YWHAZ'}
            for iso in map1433.values(): self.dict['PPI'][iso] = {}

            ## MS Data First ##
            for ms in self.dict['MSData']:
                for exp in map1433:
                    iso = map1433[exp]
                    if self.dict['MSData'][ms][experiments.index(exp)]:
                        if ms not in self.dict['PPI'][iso]: self.dict['PPI'][iso][ms] = ['MS']

            ## Literature data - 1433 ##
            addtoppi = rje.sortKeys(self.dict['PPI'])
            while addtoppi:
                ## Convert to HPRD (Literature) ID ##
                hub = addtoppi.pop(0)
                if hub in hprd.dict['Mapping']: hid = hprd.dict['Mapping'][hub]
                elif hub in cards.dict['GeneCard'] and cards.dict['GeneCard'][hub]['HPRD']: hid = cards.dict['GeneCard'][hub]['HPRD']
                else: continue  # No mapping of interactor onto HPRD
                ## Add literature interactions ##
                for litint in hprd.dict['PPI'][hid]:
                    gene = hprd.dict['HPRD'][litint]['gene']
                    if hub not in self.dict['PPI']: self.dict['PPI'][hub] = {gene:['Lit']}
                    elif gene not in self.dict['PPI'][hub]: self.dict['PPI'][hub][gene] = ['Lit']
                    elif 'Lit' not in self.dict['PPI'][hub][gene]: self.dict['PPI'][hub][gene].append('Lit')
                    if gene not in self.dict['PPI'] and gene not in addtoppi:   # Add if a 14-3-3 or CTRL interactor
                        for iso in map1433.values():
                            if gene in self.dict['PPI'][iso]:
                                addtoppi.append(gene)
                                break
                self.log.printLog('\r#PPI','%s Hub interactors' % rje.integerString(len(self.dict['PPI'])),newline=False,log=False)
            self.log.printLog('\r#PPI','%s Hub interactors' % rje.integerString(len(self.dict['PPI'])))

            ## Remove boring secondaries = only interacting with one primary/1433/control ##
            minlink = 2     # Min. no. of primaries/must link too
            safelist = ['YWHAE', 'YWHAB', 'YWHAH', 'YWHAQ', 'YWHAZ', 'SFN', 'YWHAG']
            mustlink = rje.sortKeys(self.dict['PPI'])
            # Setup boring list = reduced as good links found
            boring = []     
            for hub in rje.sortKeys(self.dict['PPI']):
                for spoke in self.dict['PPI'][hub]:
                    if spoke not in self.dict['PPI'] and spoke not in boring: boring += [spoke] * minlink
            # Reduce boring list as good links found
            for hub in rje.sortKeys(self.dict['PPI']):
                for spoke in self.dict['PPI'][hub]:
                    if spoke in boring: boring.remove(spoke)
                    while hub in safelist and spoke in boring: boring.remove(spoke)     # Cannot happen as would be in PPI?
            # Remove boring list from PPI
            bx = 0
            for hub in rje.sortKeys(self.dict['PPI']):
                for spoke in rje.sortKeys(self.dict['PPI'][hub]):
                    if spoke in boring:
                        self.dict['PPI'][hub].pop(spoke)
                        bx += 1
                if not self.dict['PPI'][hub]: self.dict['PPI'].pop(hub)
            self.log.printLog('#PPI','Removed %s boring secondary interactors' % rje.integerString(bx))

            ### Cleanup Interaction Dictionary ###
            ## Combine MS and Lit where appropriate ##
            self.dict['PPI'] = rje_ppi.combineTypes(self,self.dict['PPI'])
            ## Reciprocate interactions ##
            self.dict['PPI'] = rje_ppi.reciprocate(self,self.dict['PPI'])
            ## Merge synonyms ##
            self.dict['PPI'] = rje_ppi.mergeSynonyms(self,self.dict['PPI'],self.dict['Synonym'])
            ## Remove boring secondaries and boring control primaries ##
            self.dict['PPI'] = rje_ppi.removeSingletons(self,self.dict['PPI'],safelist=['YWHAE', 'YWHAB', 'YWHAH', 'YWHAQ', 'YWHAZ', 'SFN', 'YWHAG'])
            ## Merge by topology ##
            self.dict['PPI'] = rje_ppi.compressTopology(self,self.dict['PPI'],safelist=map1433.values())

            ### Define Node and Edge dictionaries of attributes ###
            self.dict['Node'] = {'Type':{},'Hub':{},'Size':{}}      # Node "Presence" has been abandoned for now
            self.dict['Edge'] = {}                  # Not having edge confidences because of combine topologies!
            for hub in self.dict['PPI']:    # What about spoke-onlies! ?
                if hub == 'CTRL':
                    self.dict['Node']['Type'][hub] = 'Control'
                    self.dict['Node']['Hub'][hub] = 'CTRL'
                    self.dict['Node']['Size'][hub] = 2
                elif hub in map1433.values():
                    self.dict['Node']['Type'][hub] = 'Symbol'
                    self.dict['Node']['Hub'][hub] = '1433'
                    self.dict['Node']['Size'][hub] = 2
                else:
                    self.dict['Node']['Hub'][hub] = 'Secondary'
                    for spoke in self.dict['PPI'][hub]:
                        if spoke in map1433.values(): self.dict['Node']['Hub'][hub] = 'Primary'
                    self.dict['Node']['Size'][hub] = len(string.split(hub,'/'))
                    symbol = False
                    protein = False
                    for gene in string.split(hub,'/'):
                        if gene in cards.dict['GeneCard'] and gene == cards.dict['GeneCard'][gene]['Symbol']: symbol = True
                        else: protein = True
                    if symbol and protein: self.dict['Node']['Type'][hub] = 'Mixed'
                    elif symbol: self.dict['Node']['Type'][hub] = 'Symbol'
                    else: self.dict['Node']['Type'][hub] = 'Protein'
                for spoke in self.dict['PPI'][hub]:
                    if spoke in self.dict['PPI']: continue
                    self.dict['Node']['Hub'][spoke] = 'Secondary'
                    if hub in map1433.values(): self.dict['Node']['Hub'][spoke] = 'Primary'
                    self.dict['Node']['Size'][spoke] = len(string.split(spoke,'/'))
                    symbol = False
                    protein = False
                    for gene in string.split(spoke,'/'):
                        if gene in cards.dict['GeneCard'] and gene == cards.dict['GeneCard'][gene]['Symbol']: symbol = True
                        else: protein = True
                    if symbol and protein: self.dict['Node']['Type'][spoke] = 'Mixed'
                    elif symbol: self.dict['Node']['Type'][spoke] = 'Symbol'
                    else: self.dict['Node']['Type'][spoke] = 'Protein'

            ### Output to files for Cytoscape ###
            rje_ppi.cytoscapeSave(self,basefile='ms_1433',ppi=self.dict['PPI'],node=self.dict['Node'],edge=self.dict['Edge'],bipolar=map1433.values())

        except:
            self.log.errorLog('Error in RJE_1433.combineHPRD()',printerror=True,quitchoice=False)
#########################################################################################################################
### End of SECTION II: RJE_1433 Class                                                                                   #
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
    try: RJE_1433(mainlog,cmd_list).run()
        
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
