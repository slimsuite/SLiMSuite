#!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMFRAP
Description:  SLiMFinder Results Analysis Pipeline
Version:      0.1
Last Edit:    06/10/09
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    Customised analysis pipeline for High Throughput SLiMFinder runs.

Commandline:
    ### ~ INPUT OPTIONS ~ ###
    resfiles=LIST   : List of results files to read in and process []
    resdir=PATH     : Path to extra SLiMFinder results files (*.occ.csv) []
    pairwise=FILE   : Pingu Pairwise PPI file [None]
    elm=FILE        : File of ELM Motifs [None]
    nandr=FILE      : File of Neduva & Russell motifs [None]
    pickledata=FILE : Genemap pickle to import and use. (See other rje_genemap options) [None]
    domtable=FILE   : Table of domain predictions [None]
    direct=LIST     : Evidence codes for direct interactions []
    indirect=LIST   : Evidence codes for indirect interactions []
    usego=T/F       : Whether to use GO datasets [False]
    usenandr=T/F    : Whether to use Neduva & Russell results [False]
    ### ~ PROCESSING OPTIONS ~ ###
    generics=LIST   : List of pattern annotation to screen as generic motifs []
    knowns=FILE     : File containing known SLiMs with columns "Pattern" & "Known"
    tp=FILE         : File containing known SLiMs with columns "Known", "Hub" & "Type"
    annotate=T/F    : Whether to manually annotate/classify SLiMs [False]
    skipknown=T/F   : Whether to skip known SLiMs in manual annotation (True), or ask to check (False) [False]
    skipann=LIST    : List of auto-annotations to skip from manual annotation [rand] 
    screenelm=LIST  : List of ELMs to screen from CompariMotif [LIG_Sin3_3,TRG_NES_1]
    fdranncut=X     : FDR cut-off for manual annotation [0.05]
    slimjim=T/F     : Whether to process data for SLiMJIM output [True]
    ### ~ OUTPUT OPTIONS ~ ###
    frapout=T/F     : Whether to generate SLiMFRAP 0.0 output tables [False]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genemap, rje_go, rje_slim, rje_zen, bob, slimjim
import comparimotif_V3 as comparimotif
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added SLiMJIM output.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Improved generics & classification output.
    # [ ] : ... Add auto-classification of generics?
    # [Y] : Put FDR numbers into main tables.
    # [Y] : Fix problems of comdomppi and domppi not appearing in output (vs y2hdom and bindom)
    # [ ] : Replace screening of specific known motifs.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SLiMFRAP', '0.1', 'October 2009', '2009')
    description = 'SLiMFinder Results Analysis Pipeline'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMFRAP Class                                                                                          #
#########################################################################################################################
class SLiMFRAP(bob.Bob):    #rje.RJE_Object):     
    '''
    SLiMFRAP Class. Author: Rich Edwards (2009).

    Info:str
    - DomTable = Table of domain predictions [None]
    - ELM = File of ELM Motifs []
    - NandR = File of Neduva & Russell motifs []
    - Pairwise = Pingu Pairwise PPI file []
    - ResDir = Path to extra SLiMFinder results files (*.occ.csv) []
    - Knowns = File containing known SLiMs with columns "Pattern" & "Known"
    - TP = File containing known SLiMs with columns "Known" & "Hub"
    
    Opt:boolean
    - Annotate = Whether to manually annotate/classify SLiMs [False]
    - FRAPOut = Whether to generate SLiMFRAP 0.0 output tables [False]
    - SkipAnn = List of auto-annotations to skip from manual annotation [rand] 
    - SkipKnown = Whether to skip known SLiMs in manual annotation (True), or ask to check (False) [False]
    - SLiMJIM = Whether to perform SLiMJIM processing [True]
    - UseGO = Whether to use GO datasets [False]
    - UseNandR = Whether to use Neduva & Russell results [False]

    Stat:numeric
    - FDRAnnCut = FDR cut-off for manual annotation [0.05]

    List:list
    - Direct = Evidence codes for direct interactions []
    - Generics = List of pattern annotation to screen as generic motifs []
    - Indirect = Evidence codes for indirect interactions []
    - ResFiles = List of results files to read in and process []
    - ScreenELM = List of ELMs to screen from CompariMotif [LIG_Sin3_3,TRG_NES_1]

    Dict:dictionary
    - Splits = Main, Clouds and TopRank tables split by Analysis

    Obj:RJE_Objects
    - DB = rje_DB.Database Object
    - GeneMap = rje_genemap.GeneMap object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Pairwise','ResDir','Knowns','TP','ELM','NandR','DomTable']
        self.optlist = ['Annotate','SkipKnown','FRAPOut','SLiMJIM','UseGO','UseNandR']
        self.statlist = ['FDRAnnCut']
        self.listlist = ['Generics','ResFiles','ScreenELM','SkipAnn','Direct','Indirect']
        self.dictlist = ['Splits']
        self.objlist = ['DB','GeneMap']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.info['Basefile'] = 'slimfrap'
        self.list['SkipAnn'] = []
        self.list['ScreenELM'] = ['LIG_Sin3_3','TRG_NES_1']
        self.setOpt({'SLiMJIM':True})
        self.obj['GeneMap'] = rje_genemap.GeneMap(self.log,self.cmd_list)
#########################################################################################################################
    def _cmdList(self): self._cmdListFull()     ### Sets Attributes from commandline
    def _cmdListFull(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['Pairwise','Knowns','TP','ELM','NandR','DomTable'])
                self._cmdReadList(cmd,'path',['ResDir'])
                self._cmdReadList(cmd,'opt',['Annotate','SkipKnown','FRAPOut','SLiMJIM','UseGO','UseNandR'])
                self._cmdReadList(cmd,'stat',['FDRAnnCut'])
                self._cmdReadList(cmd,'list',['Generics','ScreenELM','SkipAnn'])
                self._cmdReadList(cmd,'glist',['ResFiles'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Basic data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['FRAPOut']:
                self.summariseResults()
                self.patternInfo()
                self.annotationInfo()
                self.datasetInfo()
                self.cloudInfo()
            ## ~ [2b] ~ UPNum = SeqNum only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.i() >= 1 and rje.yesNo('Run 1-to-1?',default='N'): self.run1to1()
            ## ~ [2c] ~ Screen generics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Generics']:
                self.info['Basefile'] = string.replace(self.info['Basefile'],'.1to1','') + '.screened'
                self.setup(reformat_main=False)
                #self.reformatMain(screen_generics=True,classify=False)
                if self.opt['FRAPOut']:
                    self.summariseResults()
                    self.patternInfo()
                    self.annotationInfo()
                    self.datasetInfo()
                    self.cloudInfo()
            ### ~ [3] ~ SLiMJIM Processing and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SLiMJIM']: self.slimJIM()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def run1to1(self):  ### Reduce data to UPNum=SeqNum and re-run
        '''Reduce data to UPNum=SeqNum and re-run.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Basefile'] = self.info['Basefile'] + '.1to1'
            self.setup()
            self.dbdata().makeField('SeqNum-UPNum','UPXS')
            self.dbdata().dropEntries(['UPXS>0'])
            self.dbdata().deleteField('UPXS')
            self.reformatMain(classify=False)
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.summariseResults()
            self.patternInfo()
            self.datasetInfo()
            self.cloudInfo()
            self.annotationInfo()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def db(self): return self.obj['DB']
    def dbdata(self): return self.db().getTable('Main')
    def dbtop(self): return self.db().getTable('TopRank')
    def dbclouds(self): return self.db().getTable('Clouds')
    def dbppi(self): return self.db().getTable('PPI')
    def dbmap(self): return self.db().getTable('SeqMap')
    def mapGene(self,id): return self.obj['GeneMap'].bestMap(id)
#########################################################################################################################
    def setup(self,reformat_main=True):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            fieldformat = {}    #'Dataset','RunID','Masking','Build','RunTime','Pattern',
            for f in ['IC','ExpUP','Prob','Sig','Cons_mean','HomNum_mean','GlobID_mean','LocID_mean','Hyd_mean','IUP_mean','SA_mean']: fieldformat[f] = 'num'
            for f in ['SeqNum','UPNum','AANum','MotNum','Rank','Occ','Support','UP','CloudSeq','CloudUP']: fieldformat[f] = 'int'
            self.dict['FieldFormat'] = fieldformat
            ### ~ [2] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Main results data from previous run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            needmain = self.opt['Force'] or self.opt['Annotate']
            needmain = needmain or not os.path.exists('%s.full.csv' % self.info['Basefile'])
            needmain = needmain or not os.path.exists('%s.clouds.csv' % self.info['Basefile'])
            needmain = needmain or not os.path.exists('%s.toprank.csv' % self.info['Basefile'])
            if not needmain:
                self.db().addTable('%s.full.csv' % self.info['Basefile'],mainkeys=['Dataset','RunID','Rank','Pattern'],name='Main').dataFormat(fieldformat)
                self.db().addTable('%s.clouds.csv' % self.info['Basefile'],mainkeys=['KCloud'],name='Clouds').dataFormat(fieldformat)
                self.db().addTable('%s.toprank.csv' % self.info['Basefile'],mainkeys=['Dataset'],name='TopRank').dataFormat(fieldformat)
                reformat_main = False
            ## ~ [2b] Main results data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                mainres = []
                for file in self.list['ResFiles']:
                    mainres.append(self.db().addTable(file,mainkeys=['Dataset','RunID','Rank','Pattern']))
                mainres[0].info['Name'] = 'Main'
                while len(mainres) > 1: self.db().mergeTables(mainres[0],mainres.pop(1),matchfields=False)
            ## ~ [2c] Neduva & Russell data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nedfile = '%s.nandr.csv' % self.info['Basefile']
                if rje.checkForFile(nedfile):
                    ndb = self.db().addTable(nedfile,mainkeys=['Dataset','RunID','Rank','Pattern'],name='NandR')
                    if self.opt['UseNandR']: self.db().mergeTables(mainres[0],ndb,matchfields=False)
                elif rje.checkForFile(self.info['NandR']):
                    ndb = self.db().addEmptyTable('NandR',['Dataset','RunID','Rank','Pattern','SeqNum','UPNum','MotNum','Pattern','IC','Support','UP','Prob','Sig','Cloud'])
                    ndb.newKey(['Dataset','RunID','Rank','Pattern'])
                    bonf = [0,20]
                    while len(bonf) < 10: bonf.append(bonf[-1] * 20 * 4)  ## Approx bonf correction
                    #x#self.obj['GeneMap'] = rje_genemap.GeneMap(self.log,self.cmd_list)
                    nlines = self.loadFromFile(filename=self.info['NandR'],chomplines=True)
                    (nx,ntot) = (0.0,len(nlines))
                    for line in nlines:
                        self.progLog('\r#NED','Processing Neduva results: %.2f%%' % (nx/ntot)); nx += 100.0
                        edata = string.split(line)
                        entry = {'Dataset':edata[0],'Pattern':edata[1],'Support':edata[3],'Prob':edata[5][2:-1]}
                        entry['Pattern'] = string.replace(entry['Pattern'],'x','.')
                        entry['IC'] = len(entry['Pattern']) - string.count(entry['Pattern'],'.')
                        try:
                            if string.atof(entry['Prob']): entry['Sig'] = rje.logBinomial(1,bonf[entry['IC']],string.atof(entry['Prob']),callobj=self)
                            else: entry['Sig'] = 0.0
                        except:
                            try: entry['Sig'] = rje.logPoisson(1,bonf[entry['IC']],string.atof(entry['Prob']),callobj=self)
                            except: entry['Sig'] = 0.0; self.deBug(line); self.deBug('%s: %s' % (newkey,entry))
                        entry['Sig'] = max(entry['Sig'],string.atof(entry['Prob']))
                        entry['IC'] = float(entry['IC'])
                        dset = string.split(entry['Dataset'],'_')
                        entry['RunID'] = dset[0]
                        entry['Cloud'] = entry['Rank'] = dset[-1]
                        if dset[0] == 'NRHsD': entry['Dataset'] = string.join(dset[1:-1],'_') + '.neddom'
                        else: entry['Dataset'] = self.mapGene(string.join(dset[1:3],'')) + '.ned'
                        (entry['Support'],entry['SeqNum']) = rje.matchExp('\((\d+)/(\d+)\)',entry['Support'])
                        entry['UP'] = entry['Support']; entry['UPNum'] = entry['SeqNum']
                        newkey = ndb.makeKey(entry)
                        ndb.dict['Data'][newkey] = entry
                        #self.deBug('%s: %s' % (newkey,entry))
                    self.printLog('\r#NED','Processing Neduva results complete.',log=False)
                    ndb.saveToFile(nedfile)
                    if self.opt['UseNandR']: self.db().mergeTables(mainres[0],ndb,matchfields=False)
            ## ~ [2c] Pairwise PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.checkForFile(self.info['Pairwise']):
                self.db().addTable(self.info['Pairwise'],['Hub','Spoke'],name='PPI')
                self.db().addTable(self.info['Pairwise'],['SpokeSeq'],['SpokeUni','Spoke'],name='SeqMap',lists=True)
            ## ~ [2d] Annotation data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.checkForFile(self.info['Knowns']):
                kdb = self.db().addTable(self.info['Knowns'],mainkeys=['Pattern','Known'],name='Knowns')
            else:
                kdb = self.db().addEmptyTable('Knowns',['Pattern','Known'])
                kdb.newKey(['Pattern','Known'])
            kdb.index('Pattern')
            if rje.checkForFile(self.info['TP']):
                tdb = self.db().addTable(self.info['TP'],mainkeys=['Known','Hub'],name='TP')
            else:
                tdb = self.db().addEmptyTable('TP',['Known','Hub','Type'])
                tdb.newKey(['Known','Hub'])
            tdb.index('Known'); tdb.index('Hub'); tdb.index('Type')
            ## ~ [2e] Domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.checkForFile(self.info['DomTable']):
                domdb = self.db().addTable(self.info['DomTable'],mainkeys=['Type','Name','Start','End'],name='Domains')
            ### ~ [3] CompariMotif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            comfile = '%s.elmcomp.compare.tdt' % self.info['Basefile']
            patfile = '%s.patterns.txt' % self.info['Basefile']
            if rje.checkForFile(self.info['ELM']) and (self.opt['Force'] or not rje.checkForFile(comfile)):
                self.dbdata().index('Pattern')
                patterns = rje.sortKeys(self.dbdata().dict['Index']['Pattern'])
                if '!' in patterns: patterns.remove('!')
                if '-' in patterns: patterns.remove('-')
                open(patfile,'w').write(string.join(patterns+[''],'\n'))
                comcmd = self.cmd_list+['motifs=%s' % patfile,'searchdb=%s' % self.info['ELM'],'resfile=%s.elmcomp' % self.info['Basefile']]
                comcmd += ['matchfix=2','motdesc=2','xgmml=F']
                comparimotif.CompariMotif(self.log,comcmd).run()
            if rje.checkForFile(comfile):
                cdb = self.db().addTable(comfile,mainkeys=['Name1','Name2'],name='CompariMotif')
                for bad in self.list['ScreenELM']:
                    self.printLog('#ELM','Screened out hits from ELM "%s"' % bad)
                    cdb.dropEntries(['Name2==%s' % bad])
                cdb.index('Motif1'); cdb.index('Motif2'); cdb.index('Name2')
                edb = self.db().copyTable('CompariMotif','BestELM')
                edb.rankFieldByIndex('Motif1','Score','Rank',rev=True,absolute=True,lowest=True)
                edb.dropEntries(['Rank<1'])
            ### ~ [4] Reformat main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if reformat_main: self.reformatMain()
            else: self.reformatMain(screen_generics=True,classify=False)
            ### ~ [5] Read in occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Force'] and os.path.exists('%s.occ.csv' % self.info['Basefile']):
                self.db().addTable('%s.occ.csv' % self.info['Basefile'],mainkeys=['Dataset','Rank','Pattern','Seq','Start_Pos','End_Pos'],name='Occ')                
            else:                
                rdb = self.db().getTable('TopRank')
                odb = self.db().addEmptyTable('Occ',['Dataset','Rank','Pattern','Sig','Seq','Start_Pos','End_Pos','Prot_Len','Match','Variant','MisMatch','Desc','Cons','HomNum','GlobID','LocID','Hyd','IUP','SA','PepSeq','PepDesign'])
                orig2new = {}
                (ex,enum,rx) = (0.0,rdb.entryNum(),0)
                for entry in rdb.entries():
                    if entry['Rank'] == 0: continue
                    occfile = rje.makePath(self.info['ResDir'] + entry['Original'] + '.occ.csv',wholepath=True)
                    #x#self.deBug('%s: %s' % (occfile,os.path.exists(occfile)))
                    if not os.path.exists(occfile): continue
                    rx += 1
                    orig2new[entry['Original']] = entry['Dataset']
                    odb.loadDataDict(occfile,mainkeys=['Dataset','Rank','Pattern','Seq','Start_Pos','End_Pos'],add=True,screen=False)
                    self.progLog('\r#OCC','Loading %s occurrences: %.2f%%' % (rje.integerString(odb.entryNum()),ex/enum)); ex += 100.0
                self.printLog('\r#OCC','Loaded %s occurrences from %s datasets.' % (rje.integerString(odb.entryNum()),rje.integerString(rx)))
                odb.list['Fields'].insert(0,'Hub')
                odb.list['Fields'].insert(2,'Original')
                for entry in odb.entries():
                    entry['Original'] = entry['Dataset']
                    entry['Dataset'] = orig2new[entry['Dataset']]
                    entry['Hub'] = string.join(string.split(entry['Dataset'],'.')[:-1],'.')
                odb.saveToFile('%s.occ.csv' % self.info['Basefile'])
        except: self.errorLog('Problem during %s setup.' % self); raise
#########################################################################################################################
    ### <3> ### Data reformatting Methods                                                                               #
#########################################################################################################################
    def reformatMain(self,screen_generics=False,classify=True):     ### Reformats the main data table
        '''Reformats the main data table.'''
        try:### ~ [1] Basic reformatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fieldformat = {}    #'Dataset','RunID','Masking','Build','RunTime','Pattern',
            for f in ['IC','ExpUP','Prob','Sig','Cons_mean','HomNum_mean','GlobID_mean','LocID_mean','Hyd_mean','IUP_mean','SA_mean']: fieldformat[f] = 'num'
            for f in ['SeqNum','UPNum','AANum','MotNum','Rank','Occ','Support','UP','CloudSeq','CloudUP']: fieldformat[f] = 'int'
            self.dbdata().dataFormat(fieldformat)
            self.dict['Knowns'] = {}
            self.list['Annotation'] = []
            ### ~ [2] Full data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Annotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb = self.db().getTable('CompariMotif')
            edb = self.db().getTable('BestELM')
            iwant = ['bin','y2h','ppi','com','ned','bindom','y2hdom','ppidom','comdom','neddom','go']
            self.dbdata().list['Fields'] += ['Original','DType','Hub','PPISet','Annotation','ELM','BestELM','ELMSim','ELMPattern','Score','Match','MatchIC']
            for entry in self.dbdata().entries():
                #if entry['Dataset'][:4] == 'rupc' and entry['Dataset'][:5] != 'rupc_': entry['Dataset'] = 'rupc_' + entry['Dataset'][4:]
                y2h = entry['Dataset'].find('y2h') > 0
                #x#if y2h and entry['Dataset'][:4] == 'rupc': self.deBug(entry)
                entry['Original'] = entry['Dataset']
                if entry['Dataset'][:4] in ['rupc','rseq']:
                    entry['DType'] = entry['Dataset'][:4]
                    etype = string.split(entry['Dataset'][4:],'_')[0]
                    edset = string.join([entry['Dataset'][:4]]+string.split(entry['Dataset'],'_')[1:],'_')
                    if etype in ['dom','']: etype = 'ppi' + etype
                else:
                    entry['DType'] = 'real'
                    etype = string.split(entry['Dataset'],'.')[-1]
                    edset = string.join(string.split(entry['Dataset'],'.')[:-1],'.')
                if etype == 'domppi': etype = 'ppidom'
                if etype == 'comppi': etype = 'com'
                if etype == 'comdomppi': etype = 'comdom'
                if entry['Dataset'][:3] == 'GO_': etype = 'go'; edset = entry['Dataset']
                entry['Dataset'] = '%s.%s' % (edset,etype)
                (entry['Hub'],entry['PPISet']) = (edset,etype)                
                if entry['Pattern'] not in ['!','-']:
                    alist = self.annotationList(entry)
                    for a in alist:
                        if a not in self.list['Annotation']: self.list['Annotation'].append(a)#; self.deBug(self.list['Annotation'])
                else: alist = ['']
                #x#if y2h and entry['Dataset'][:4] == 'rupc': self.deBug(entry)
                if edb and entry['Pattern'] in edb.dict['Index']['Motif1']:
                    entry['ELM'] = len(cdb.dict['Index']['Motif1'][entry['Pattern']])
                    bestelm = edb.data()[edb.dict['Index']['Motif1'][entry['Pattern']][0]]
                    entry['BestELM'] = bestelm['Name2']
                    entry['ELMSim'] = bestelm['Sim2']
                    entry['ELMPattern'] = bestelm['Motif2']
                    for field in ['Score','Match','MatchIC']: entry[field] = bestelm[field]
                else:
                    for field in ['ELM','BestELM','ELMSim','ELMPattern','Score','Match','MatchIC']: entry[field] = ''
                    for field in ['ELM','Score','MatchIC']: entry[field] = 0
                entry['Annotation'] = string.join(alist,',')
                if entry['PPISet'] not in iwant: self.deBug(entry)
            if not self.opt['UseGO']: self.dbdata().dropEntries(['PPISet==go','PPISet==rseqgo','PPISet==rupcgo'])
            if not self.opt['UseNandR']:  self.dbdata().dropEntries(['PPISet==ned','PPISet==neddom'])
            ## ~ [2b] Screen Generics / Classification and extra annotation ~~~~~~~~~~~~~~~~~~~~~~~ ##
            if screen_generics: self.screenGenerics()
            self.dbdata().joinFields('KCloud',['Dataset','Cloud'],join='|')   #x#self.db().info['Delimit'])
            self.dbdata().index('KCloud')
            if classify: self.classifyClouds()
            ## ~ [2c] Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Annotation'].sort()#; self.deBug(self.list['Annotation'])
            self.dbdata().joinFields('DCloud',['Hub','Cloud'],join='|')   #x#self.db().info['Delimit'])
            self.dbdata().index('DCloud')
            #self.dbdata().saveToFile('%s.full.csv' % self.info['Basefile'])    ### Saves data to delimited file
            self.dbdata().joinFields('Analysis',['DType','PPISet'],join='.')
            ### ~ [3] Generate Cloud Data Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.db().copyTable('Main','Clouds')
            cdb.dropEntries(['Rank<1'])
            cdb.index('KCloud')
            newdata = {}
            for cloud in cdb.dict['Index']['KCloud']:
                newdata[cloud] = cdb.data()[cdb.dict['Index']['KCloud'][cloud][0]]
                cloudann = string.split(newdata[cloud]['Annotation'],',')
                for dkey in cdb.dict['Index']['KCloud'][cloud][1:]:
                    entry = cdb.data()[dkey]
                    if entry['Rank'] < newdata[cloud]['Rank']: newdata[cloud] = entry
                    for a in string.split(entry['Annotation'],','):
                        if a not in cloudann: cloudann.append(a)
                cloudann.sort()
                newdata[cloud]['Annotation'] = string.join(cloudann,',')
            cdb.dict['Data'] = newdata
            cdb.list['Keys'] = ['KCloud']
            cdb.list['Fields'].remove('KCloud')
            cdb.list['Fields'].insert(0,'KCloud')
            self.printLog('#CLOUD','%s Cloud table entries from %s Main entries' % (rje.integerString(cdb.entryNum()),rje.integerString(self.dbdata().entryNum())))            
            ### ~ [4] Generate TopRank Data Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = self.db().copyTable('Clouds','TopRank')
            rdb.dropEntries(['Rank>1'])
            rdb.newKey(['Dataset'])
            cdb.list['Fields'].remove('Dataset')
            cdb.list['Fields'].insert(0,'Dataset')
            self.printLog('#TOP','%s TopRank table entries from %s Cloud entries' % (rje.integerString(rdb.entryNum()),rje.integerString(cdb.entryNum())))
            ### ~ [5] FDR calculations and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.splitFDR(force=True)
            self.dbdata().list['Fields'].insert(self.dbdata().list['Fields'].index('Sig')+1,'FDR')
            self.dbdata().index('Analysis')
            cdb.list['Fields'].insert(cdb.list['Fields'].index('Sig')+1,'FDR')
            cdb.index('Analysis')
            rdb.list['Fields'].insert(rdb.list['Fields'].index('Sig')+1,'FDR')
            rdb.index('Analysis')
            splits = self.makeSplits()
            for analysis in splits['Clouds']:
                ctable = splits['Clouds'][analysis]
                ctable.index('KCloud')
                for entry in cdb.entryList(cdb.dict['Index']['Analysis'][analysis]):
                    try:
                        if entry['KCloud'] in ctable.dict['Index']['KCloud']: entry['FDR'] = ctable.data()[ctable.dict['Index']['KCloud'][entry['KCloud']][0]]['FDR']
                    except:
                        self.deBug(entry)
                        self.deBug(ctable.dict['Index']['KCloud'][entry['KCloud']])
                        self.deBug(ctable.data()[ctable.dict['Index']['KCloud'][entry['KCloud']][0]])
                for entry in self.dbdata().entryList(self.dbdata().dict['Index']['Analysis'][analysis]):
                    try:
                        if entry['KCloud'] in ctable.dict['Index']['KCloud']: entry['FDR'] = ctable.data()[ctable.dict['Index']['KCloud'][entry['KCloud']][0]]['FDR']
                    except: pass
                rtable = splits['TopRank'][analysis]
                rtable.index('KCloud')
                for entry in rdb.entryList(rdb.dict['Index']['Analysis'][analysis]):
                    try:
                        if entry['KCloud'] in rtable.dict['Index']['KCloud']: entry['FDR'] = rtable.data()[ctable.dict['Index']['KCloud'][entry['KCloud']][0]]['FDR']
                    except: pass
            self.dbdata().saveToFile('%s.full.csv' % self.info['Basefile'])    ### Saves data to delimited file
            cdb.saveToFile('%s.clouds.csv' % self.info['Basefile'])    ### Saves data to delimited file
            rdb.saveToFile('%s.toprank.csv' % self.info['Basefile'])    ### Saves data to delimited file
            return
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def patternCompariMotifText(self,pattern):  ### Returns summary comparimotif text for pattern
        '''Returns summary comparimotif text for pattern.'''
        cdb = self.db().getTable('CompariMotif'); ptxt = ''
        if cdb and pattern in cdb.dict['Index']['Motif1']:
            for ckey in cdb.dict['Index']['Motif1'][pattern]:
                centry = cdb.data()[ckey]
                ptxt += '|--- %s (%s) %s vs ||%s|| %s %s\n' % (centry['Score'],centry['NormIC'],centry['Match'],centry['Motif2'],centry['Name2'],centry['Desc2'])
        return ptxt
#########################################################################################################################
    def patternHubs(self,pattern,withsig=True):  ### Returns list of Hubs returning pattern
        '''Returns list of Hubs returning pattern.'''
        db = self.dbdata(); hubs = []
        if 'Pattern' not in db.dict['Index']: db.index('Pattern')
        for ekey in db.dict['Index']['Pattern'][pattern]:
            try:
                if withsig: hubs.append('%s (%s)' % (db.data()[ekey]['Hub'],rje_slim.expectString(db.data()[ekey]['Sig'])))
                else: hubs.append(db.data()[ekey]['Hub'])
            except:
                self.errorLog('Error during patternHubs(%s)' % pattern)
                self.deBug(ekey)
                self.deBug(db.data()[ekey])
        return rje.sortUnique(hubs)        
#########################################################################################################################
    def patternScore(self,pattern,score='Sig',smallest=True,return_entry=False):   ### Returns best score for pattern
        '''Returns best score for pattern.'''
        try:
            db = self.dbdata(); first = True
            if 'Pattern' not in db.dict['Index']: db.index('Pattern')
            for ekey in db.dict['Index']['Pattern'][pattern]:
                entry = db.data()[ekey]
                escore = entry[score]
                if first: pscore = escore; first = False; pentry = entry
                elif smallest == (escore < pscore): pscore = escore; pentry = entry
            if return_entry: return pentry
            return pscore        
        except: self.errorLog('Bugger'); return 0.0
#########################################################################################################################
    def classifyClouds(self):   ### Classify clouds and add extra annotation
        '''Classify clouds and add extra annotation.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.db().getTable('CompariMotif')
            ## ~ [0a] Known motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            knowns = {}     # Dictionary of {pattern:known}
            kdb = self.db().getTable('Knowns')
            for pattern in kdb.dict['Index']['Pattern']:
                for pkey in kdb.dict['Index']['Pattern'][pattern]:
                    known = kdb.data()[pkey]['Known']
                    if pattern not in knowns: knowns[pattern] = known
                    elif known != knowns[pattern]:
                        self.errorLog('Known "%s" matched to "%s" pattern %s' % (known,knowns[pattern],pattern),printerror=False)
                        if self.i() >= 0 and rje.yesNo('Replace "%s" for %s with "%s"?' % (knowns[pattern],pattern,known)): knowns[pattern] = known
            self.deBug('Knowns: %s' % knowns)
            self.dbdata().index('Pattern')
            patcheck = rje.sortKeys(self.dbdata().dict['Index']['Pattern'])     # List of patterns to check/ask about
            for bad in ['!','-']:
                if bad in patcheck: patcheck.remove(bad)
            ## ~ [0b] Motif classification ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tp = {}         # Dictionary of {known:[datasets]}
            ot = {}         # Dictionary of {known:[datasets]}
            ### ~ [1] Manual annotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try:
                if '?' not in self.list['SkipAnn'] and self.i() >= 0 and self.opt['Annotate'] and rje.yesNo('Skip patterns annotated as "?"?'): self.list['SkipAnn'].append('?')
            except: patcheck = []
            self.printLog('#SKIP','Autoskip: %s' % string.join(self.list['SkipAnn'],', '),log=False)
            px = len(patcheck)
            pattern = ''
            while patcheck:
                prev = pattern
                pattern = patcheck.pop(0)
                if self.i() < 0 or not self.opt['Annotate']: break
                if pattern in knowns and knowns[pattern] in self.list['SkipAnn']: continue
                generic = False
                autoann = self.annotationList({'Pattern':pattern})
                bestentry = self.patternScore(pattern,return_entry=True)
                if bestentry['DType'] in ['rupc','rseq']: autoann.append('rand')
                for a in autoann:
                    if a in self.list['Generics'] and pattern not in knowns: knowns[pattern] = a
                skipme = False
                for a in autoann: 
                    if a in self.list['SkipAnn']    : skipme = True; continue
                if skipme: continue
                ptxt = '_____\n%s of %s: %s << %s\n' % (rje.integerString(px-len(patcheck)),rje.integerString(px),pattern,string.join(self.patternHubs(pattern),'; '))
                ptxt += '|--- %s\n' % string.join(autoann,' / ')
                ptxt += self.patternCompariMotifText(pattern)
                ## ~ [1a] Skip if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    if self.stat['FDRAnnCut'] > 0.0 and self.stat['FDRAnnCut'] < self.patternScore(pattern,score='Sig'): continue
                    if pattern in knowns:
                        if knowns[pattern] != '?' and self.opt['SkipKnown']: continue
                        if knowns[pattern] != '?' and rje.yesNo('%s%s = "%s"?' % (ptxt,pattern,knowns[pattern])): continue
                        else: knowns.pop(pattern)
                ## ~ [1b] Manually annotate pattern ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    known = rje.choice('%s|\nManual annotation for "%s"?:' % (ptxt,pattern),default='',confirm=True)
                    if known == '^' and prev: patcheck.insert(0,pattern); pattern = ''; continue
                    if not known: continue
                    knowns[pattern] = known
                ## ~ [1c] Map onto rest of cloud ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    cloudpat = [pattern]; clouddone = []
                    while cloudpat:
                        cpat = cloudpat.pop(0)
                        try:
                            for ekey in self.dbdata().dict['Index']['Pattern'][cpat]:
                                dcloud = self.dbdata().data()[ekey]['KCloud']
                                if dcloud in clouddone: continue
                                clouddone.append(dcloud)
                                for dkey in self.dbdata().dict['Index']['KCloud'][dcloud]:
                                    pat2 = self.dbdata().data()[dkey]['Pattern']
                                    if pat2 in knowns and pat2 not in patcheck and (knowns[pat2] != '?' or knowns[pat2] == known): continue
                                    ptxt = '_____\n%s: %s << %s\n' % (dcloud,pat2,string.join(self.patternHubs(pat2),'; '))
                                    ptxt += '|--- %s\n' % string.join(self.annotationList({'Pattern':pat2}),' / ')
                                    ptxt += self.patternCompariMotifText(pat2)
                                    if not rje.yesNo('%s%s: %s = %s?' % (ptxt,dcloud,pat2,known)): continue
                                    knowns[pat2] = known; cloudpat.append(pat2)
                                    if pat2 in patcheck: patcheck.remove(pat2)
                        except KeyboardInterrupt: break
                        except: raise
                except KeyboardInterrupt:
                    if rje.yesNo('Exit manual annotation?'): break
                    elif rje.yesNo('Quit SLiMFrap?'): sys.exit()
                except: self.errorLog('Problem during manual annotation (%s)' % pattern,quitchoice=False); break
            ## ~ [1d] Update Database table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for pattern in knowns:
                known = knowns[pattern]
                newentry = {'Known':known,'Pattern':pattern}
                newkey = kdb.makeKey(newentry)
                if newkey in kdb.data(): continue
                kdb.data()[newkey] = newentry
            kdb.saveToFile('%s.knowns.tdt' % self.info['Basefile'])
            ### ~ [2] Annotate Known SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pattern in knowns:
                if pattern not in self.dbdata().dict['Index']['Pattern']: continue
                bestentry = self.patternScore(pattern,return_entry=True)
                for ekey in self.dbdata().dict['Index']['Pattern'][pattern]:
                    entry = self.dbdata().data()[ekey]
                    alist = string.split(entry['Annotation'],',')
                    if bestentry['DType'] in ['rupc','rseq'] and 'rand' not in alist: alist.append('rand')
                    known = knowns[pattern]
                    if known not in alist: alist.append(known)
                    if known not in self.list['Annotation']: self.list['Annotation'].append(known)#; self.deBug(self.list['Annotation'])
                    entry['Annotation'] = string.join(alist,',')
            ### ~ [3] Classify SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            classification = {}     # Dictionary of {annotation:{hub:class}}
            ## ~ [3a] Annotation / Pattern Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            adb = self.db().addEmptyTable('Annotation',['Annotation','Pattern'])
            adb.newKey(['Annotation','Pattern'])
            (ax,atot) = (0.0,self.dbdata().entryNum())
            for entry in self.dbdata().entries():
                self.progLog('\r#ANN','Cross-linking Annotation & Patterns: %.2f%%' % (ax/atot)); ax += 100.0
                pattern = entry['Pattern']
                for a in string.split(entry['Annotation'],','):
                    if a not in self.list['Annotation']: self.list['Annotation'].append(a)#; self.deBug(self.list['Annotation'])
                    akey = adb.makeKey({'Annotation':a,'Pattern':pattern})
                    adb.data()[akey] = {'Annotation':a,'Pattern':pattern}
            adb.index('Annotation')
            self.printLog('\r#ANN','Cross-linked Annotation & Patterns: %s pairs' % rje.integerString(adb.entryNum()))
            ## ~ [3b] TP/FP/OT classification ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Annotation'].sort()#; self.deBug(self.list['Annotation'])
            self.dbdata().index('Pattern')
            tdb = self.db().getTable('TP')
            for entry in tdb.entries():
                known = entry['Known']
                hub = entry['Hub']
                type = entry['Type']
                if known not in classification: classification[known] = {}
                classification[known][hub] = type
            ## ~ [3c] Add Manual stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for annotation in self.list['Annotation']:
                if self.i() < 0 or not self.opt['Annotate']: break
                try:
                    if not annotation: continue
                    if annotation in self.list['Generics'] and rje.yesNo('Classify all "%s" results as "Generic"?' % annotation):
                        if annotation not in classification: classification[annotation] = {}
                        for pkey in adb.dict['Index']['Annotation'][annotation]:
                            pattern = adb.data()[pkey]['Pattern']
                            for ekey in self.dbdata().dict['Index']['Pattern'][pattern]:
                                entry = self.dbdata().data()[ekey]
                                if entry['Hub'] in classification[annotation]: continue
                                classification[annotation][entry['Hub']] = 'Generic'
                        continue
                    if not rje.yesNo('Classify "%s" results?' % annotation): continue
                    if annotation not in classification: classification[annotation] = {}
                    for pkey in adb.dict['Index']['Annotation'][annotation]:
                        pattern = adb.data()[pkey]['Pattern']
                        if self.stat['FDRAnnCut'] > 0.0 and self.stat['FDRAnnCut'] < self.patternScore(pattern): continue
                        ptxt = '_____\n%s: %s << %s\n' % (annotation,pattern,string.join(self.patternHubs(pattern),'; '))
                        ptxt += '|--- %s\n' % string.join(autoann,' / ')
                        ptxt += self.patternCompariMotifText(pattern)
                        for ekey in self.dbdata().dict['Index']['Pattern'][pattern]:
                            entry = self.dbdata().data()[ekey]
                            if entry['DType'] != 'real': continue
                            if entry['Hub'] in classification[annotation]: continue
                            try:
                                if annotation in self.list['Generics']: known = 'Generic'
                                else: known = ''
                                known = rje.choice('SLiM Class for "%s" in %s?:' % (annotation,entry['Hub']),default=known,confirm=True)
                                if not known: continue
                                classification[annotation][entry['Hub']] = known
                            except KeyboardInterrupt:
                                if rje.yesNo('Exit manual annotation?'): raise
                                elif rje.yesNo('Quit SLiMFrap?'): sys.exit()
                            except: raise
                except KeyboardInterrupt: break
                except: raise
            ## ~ [3d] Update Database table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for known in classification:
                for hub in classification[known]:
                    newentry = {'Known':known,'Hub':hub,'Type':classification[known][hub]}
                    newkey = tdb.makeKey(newentry)
                    if newkey in tdb.data(): continue
                    tdb.data()[newkey] = newentry
            tdb.saveToFile('%s.tp.tdt' % self.info['Basefile'])
            ## ~ [3e] Update Main Database table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pref = ['TP','FP','OT','Generic']
            if 'Class' not in self.dbdata().fields():
                self.dbdata().list['Fields'].insert(self.dbdata().fields().index('Annotation')+1,'Class')
            (ex, etot) = (0.0, self.dbdata().entryNum())
            for entry in self.dbdata().entries():
                self.progLog('#CLASS','Updating entry classification: %.1f%%' % (ex/etot)); ex += 100.0
                if 'Class' not in entry: entry['Class'] = ''
                for a in string.split(entry['Annotation'],','):
                    if a in classification and entry['Hub'] in classification[a]:
                        c = classification[a][entry['Hub']]
                        if entry['Class'] == c: continue
                        if not entry['Class'] or c in pref and (entry['Class'] not in pref or pref.index(c) < pref.index(entry['Class'])): entry['Class'] = c
            self.printLog('\r#CLASS','Updating entry classification complete.')
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def tableFDR(self,table,N):     ### Calculates FDR given table and dataset number, N
        '''Calculates FDR given table and dataset number, N.'''
        table = self.db().getTable(table)
        table.rankField('Sig',lowest=True)
        for entry in table.entries():
            #if entry['Sig']:
            entry['FDR'] = (N * entry['Sig']) / entry['Sig.Rank']
            #else: entry['FDR'] = 0.0
#########################################################################################################################
    def screenGenerics(self):   ### Screens generic motifs and re-clouds
        '''Screens generic motifs and re-clouds.'''
        try:### ~ [1] Remove Generics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dkey in self.dbdata().datakeys():
                entry = self.dbdata().data()[dkey]
                alist = string.split(entry['Annotation'],',')
                for a in alist:
                    if a in self.list['Generics']: self.dbdata().data().pop(dkey); break
            ### ~ [2] Re-rank and recloud motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Rerank motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dbdata().index('Dataset')
            for dset in self.dbdata().dict['Index']['Dataset']:
                if len(self.dbdata().dict['Index']['Dataset'][dset]) == 1 and self.dbdata().data()[self.dbdata().dict['Index']['Dataset'][dset][0]]['Rank'] == 0: continue
                scorelist = []
                ekeys = self.dbdata().dict['Index']['Dataset'][dset]
                for e in ekeys: scorelist.append(self.dbdata().data()[e]['Rank'])
                ranks = rje.rankList(scorelist,rev=False,absolute=True,lowest=True)
                for e in ekeys: self.dbdata().data()[e]['Rank'] = ranks.pop(0)
            ## ~ [2b] Recloud motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            return #!# Rest of code needs occurrence data #!#
            self.dbdata().joinFields('KCloud',['Dataset','Cloud'],join='|')   #x#self.db().info['Delimit'])
            self.dbdata().index('KCloud')
            newx = 0
            for cloud in self.dbdata().dict['Index']['KCloud']:
                entries = []
                ekeys = self.dbdata().dict['Index']['Dataset'][cloud]
                for e in ekeys: entries.append(self.dbdata().data()[e])
                newx += self.reCloud(entries)
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <4> ### Data Summary Methods                                                                                    #
#########################################################################################################################
    def makeSplits(self,force=False):   ### Makes split tables by Analysis (self.dict['Splits'])
        '''Makes split tables by Analysis (self.dict['Splits']).'''
        if self.dict['Splits'] and not force: return self.dict['Splits']
        splits = {}
        for table in ['Main','Clouds','TopRank']: splits[table] = self.db().splitTable(table,'Analysis',True) #self.db().getTable(table).index('Analysis')
        self.deBug(rje.sortKeys(splits['Main']))
        for analysis in splits['Main']:
            if analysis not in splits['Clouds']: splits['Clouds'][analysis] = self.db().addEmptyTable('Clouds_%s' % analysis,self.db().getTable('Clouds').list['Fields'])
            if analysis not in splits['TopRank']: splits['TopRank'][analysis] = self.db().addEmptyTable('TopRank_%s' % analysis,self.db().getTable('TopRank').list['Fields'])
        self.dict['Splits'] = splits
        return splits
#########################################################################################################################
    def splitFDR(self,force=False):     ### Calculates FDR for Cloud & TopRank Splits
        '''Calculates FDR for Cloud & TopRank Splits.'''
        splits = self.makeSplits(force)
        for analysis in splits['Main']:
            table = splits['Main'][analysis]
            table.index('Dataset')
            try:
                ctable = splits['Clouds'][analysis]
                if analysis in ['ned','neddom']: N = 1986
                else: N = len(table.dict['Index']['Dataset'])
                self.tableFDR(ctable,N)
                self.tableFDR(splits['TopRank'][analysis],N)
            except:
                self.errorLog('?')
                ctable = None
#########################################################################################################################
    def summariseResults(self):     ### Summarises results for Main, TopRank and Clouds tables
        '''Summarises results for Main, TopRank and Clouds tables.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            summary = {}    # Summary data {key:datadict}
            sumhead = ['Analysis','DType','PPISet','DNum','NoSig']
            sumfile = '%s.summary.tdt' % self.info['Basefile']
            for t in ['p','c','f','e']:
                for p in [0.05,0.01,0.005,0.001,0.0001,0.00001]: sumhead.append('%s%s' % (t,p))
            rje.delimitedFileOutput(self,sumfile,sumhead,rje_backup=True)
            splits = self.makeSplits()
            ### ~ [2] Simple Summary Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Dataset counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for analysis in rje.sortKeys(splits['Main']):
                table = splits['Main'][analysis]
                table.index('Dataset')
                table.index('Rank')
                entry = table.entries()[0]
                summary[analysis] = datadict = {}
                for dkey in sumhead[:3]: datadict[dkey] = entry[dkey]
                try: ctable = splits['Clouds'][analysis]
                except: ctable = None
                datadict['DNum'] = len(table.dict['Index']['Dataset'])
                try: datadict['NoSig'] = len(table.dict['Index']['Rank'][0])
                except: datadict['NoSig'] = 0
                try: datadict['Sig'] = len(table.dict['Index']['Rank'][1])
                except: datadict['Sig'] = 0
                for p in [0.05,0.01,0.005,0.001,0.0001,0.00001]:
                    datadict['p%s' % p] = 0
                    datadict['c%s' % p] = 0
                    datadict['f%s' % p] = 0
                    datadict['e%s' % p] = datadict['DNum'] * p
                    if 1 not in table.dict['Index']['Rank']: continue
                    for ikey in table.dict['Index']['Rank'][1]:
                        entry = table.data()[ikey]
                        if entry['Sig'] <= p: datadict['p%s' % p] += 1
                    if not ctable: continue
                    for entry in ctable.entries():
                        if entry['Sig'] <= p: datadict['c%s' % p] += 1
                        if entry['FDR'] <= p: datadict['f%s' % p] += 1
                rje.delimitedFileOutput(self,sumfile,sumhead,datadict=datadict)
        except: self.errorLog(rje_zen.Zen().wisdom())#; raise   # Delete this if method error not terrible
#########################################################################################################################
    def patternInfo(self):  ### Outputs information on shared patterns
        '''Outputs information on shared patterns.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.db().getTable('CompariMotif')
            edb = self.db().getTable('BestELM')
            patfile = '%s.pattern.tdt' % self.info['Basefile']
            pathead = ['Pattern','Annotation']
            if edb: pathead += ['ELM','BestELM','ELMSim','ELMPattern','Score','Match','MatchIC']
            splits = self.dict['Splits']
            self.dbdata().index('Pattern')
            for analysis in splits['Main']:
                pathead += [analysis,'%s_Sig' % analysis,'%s_cSig' % analysis,'%s_cFDR' % analysis]
                splits['Main'][analysis].index('Pattern')
                try: splits['Clouds'][analysis].newKey(['DCloud'])
                except: pass
            rje.delimitedFileOutput(self,patfile,pathead,rje_backup=True)
            ### ~ [2] Summarise Pattern Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,len(self.dbdata().dict['Index']['Pattern']))
            for pattern in rje.sortKeys(self.dbdata().dict['Index']['Pattern']):
                self.progLog('\r#SLIM','Summarising %s SLiMs: %.2f%%' % (rje.integerString(stot),sx/stot)); sx += 100.0
                if pattern in ['!','-']: continue
                pdat = {'Pattern':pattern}
                if edb and pattern in edb.dict['Index']['Motif1']:
                    pdat['ELM'] = len(cdb.dict['Index']['Motif1'][pattern])
                    bestelm = edb.data()[edb.dict['Index']['Motif1'][pattern][0]]
                    pdat['BestELM'] = bestelm['Name2']
                    pdat['ELMSim'] = bestelm['Sim2']
                    pdat['ELMPattern'] = bestelm['Motif2']
                    for field in ['Score','Match','MatchIC']: pdat[field] = bestelm[field]
                else:
                    for field in ['ELM','BestELM','ELMSim','ELMPattern','Score','Match','MatchIC']: pdat[field] = ''
                    for field in ['ELM','Score','MatchIC']: pdat[field] = 0
                alist = []
                for analysis in splits['Main']:
                    table = splits['Main'][analysis]
                    pdat[analysis] = 0
                    if not table.dict['Index']['Pattern'].has_key(pattern): continue
                    pdat[analysis] = len(table.dict['Index']['Pattern'][pattern])
                    sig = 1; csig = 1; cfdr = 1
                    for tkey in table.dict['Index']['Pattern'][pattern]:
                        entry = table.data()[tkey]
                        for a in string.split(entry['Annotation'],','):
                            if a and a not in alist: alist.append(a)
                        sig = min(sig,entry['Sig'])
                        cloud = entry['DCloud']
                        if cloud not in splits['Clouds'][analysis].data():
                            self.deBug('%s: %s' % (pdat,cloud))
                            self.deBug('%s: %s' % (analysis,splits['Clouds'][analysis].data().keys()[:30]))
                        else:
                            csig = min(csig,splits['Clouds'][analysis].data()[cloud]['Sig'])
                            cfdr = min(cfdr,splits['Clouds'][analysis].data()[cloud]['FDR'])
                    pdat['%s_Sig' % analysis] = sig
                    pdat['%s_cSig' % analysis] = csig
                    pdat['%s_cFDR' % analysis] = cfdr
                    alist.sort(); pdat['Annotation'] = string.join(alist,',')
                rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            self.printLog('\r#SLIM','Summary for %s SLiMs output to %s' % (rje.integerString(stot),patfile))
        except: self.errorLog(rje_zen.Zen().wisdom())#; raise   # Delete this if method error not terrible
#########################################################################################################################
    def annotationInfo(self):  ### Outputs information on annotation types
        '''Outputs information on shared patterns.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            patfile = '%s.annotation.tdt' % self.info['Basefile']
            pathead = ['Annotation']
            splits = self.dict['Splits']
            for analysis in rje.sortKeys(splits['Main']): pathead += [analysis]
            for analysis in rje.sortKeys(splits['Main']): pathead += ['%s_Clouds' % analysis]
            rje.delimitedFileOutput(self,patfile,pathead,rje_backup=True)
            ### ~ [2] Summarise Pattern Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,len(self.list['Annotation'])+3)
            ## ~ [2a] Total SLiM/Cloud numbers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#ANN','Summarising annotation: %.2f%%' % (sx/stot)); sx += 100.0
            pdat = {'Annotation':'Datasets'}
            for analysis in splits['Main']:
                table = splits['Main'][analysis]
                pdat[analysis] = pdat['%s_Clouds' % analysis] = len(table.dict['Index']['Dataset'])
            rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            self.progLog('\r#ANN','Summarising annotation: %.2f%%' % (sx/stot)); sx += 100.0
            pdat = {'Annotation':'Total'}
            for analysis in splits['Main']:
                table = splits['Main'][analysis]
                try: pdat[analysis] = table.entryNum() - len(table.dict['Index']['Rank'][0])
                except: pdat[analysis] = table.entryNum()
                pdat['%s_Clouds' % analysis] = splits['Clouds'][analysis].entryNum()
            rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            ## ~ [2b] Annotations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for a in self.list['Annotation']+['None']:
                self.progLog('\r#ANN','Summarising annotation: %.2f%%' % (sx/stot)); sx += 100.0
                pdat = {'Annotation':a}
                if a == 'None': a = ''
                for analysis in splits['Main']:
                    pdat[analysis] = 0
                    pdat['%s_Clouds' % analysis] = 0
                    table = splits['Main'][analysis]
                    for entry in table.entries():
                        if entry['Pattern'] in ['-','!']: continue
                        if a in string.split(entry['Annotation'],','): pdat[analysis] += 1
                    table = splits['Clouds'][analysis]
                    for entry in table.entries():
                        if a in string.split(entry['Annotation'],','): pdat['%s_Clouds' % analysis] += 1
                rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            self.printLog('\r#ANN','Summary for annotation output to %s' % (patfile))
        except: self.errorLog(rje_zen.Zen().wisdom())#; raise   # Delete this if method error not terrible
#########################################################################################################################
    def datasetInfo(self):  ### Outputs information on real datasets
        '''Outputs information on real datasets.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = self.db().copyTable('Main','Real')
            rdb.dropEntries(['DType=rupc','DType=rseq','PPISet=go'])
            rdb.index('Hub')
            patfile = '%s.dataset.tdt' % self.info['Basefile']
            pathead = ['Hub','Meta','TP','OT','Generic','Total','MinFDR']
            splits = self.dict['Splits']
            for analysis in splits['Main']:
                if analysis[:4] in ['rupc','rseq']: continue
                if analysis[-2:] == 'go': continue
                splits['Main'][analysis].index('Hub',force=True)
                splits['Clouds'][analysis].index('Hub',force=True)
                splits['TopRank'][analysis].newKey(['Hub'])
                if analysis[-3:] == 'dom' or analysis[-6:-3] == 'dom': continue
                pathead += ['%s_SeqNum' % analysis,'%s_SLiM' % analysis,'%s_Clouds' % analysis,'%s_TopRank' % analysis,'%s_Sig' % analysis,'%s_FDR' % analysis]
            rje.delimitedFileOutput(self,patfile,pathead,rje_backup=True)
            ## ~ [1a] True Positives etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tdb = self.db().getTable('TP')
            tdb.index('Hub')
            ### ~ [2] Summarise Pattern Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,len(rdb.index('Hub')))
            for hub in rje.sortKeys(rdb.index('Hub')):
                self.progLog('\r#HUB','Summarising %s Hubs: %.2f%%' % (rje.integerString(stot),sx/stot)); sx += 100.0
                pdat = {'Hub':hub,'Total':0,'MinFDR':1.0}
                if hub[:4] in ['rupc','rseq']: continue
                if hub in tdb.dict['Index']['Hub']:
                    for entry in tdb.entryList(tdb.dict['Index']['Hub'][hub]):
                        if entry['Type'] not in pathead: continue
                        if entry['Type'] in pdat: pdat[entry['Type']] += '; %s' % entry['Known']
                        else: pdat[entry['Type']] = entry['Known']
                for analysis in splits['Main'].keys()[0:]:
                    try:
                        if analysis[:4] in ['rupc','rseq']: continue
                        if analysis[-2:] == 'go': continue
                        table = splits['Main'][analysis]
                        if not table.dict['Index']['Hub'].has_key(hub): pdat[analysis] = 0; continue
                        split = analysis
                        if analysis[-3:] == 'dom': pdat['Meta'] = 'Domain'; analysis = analysis[:-3]
                        elif analysis[-6:-3] == 'dom': pdat['Meta'] = 'Domain'; analysis = string.replace(analysis,'.dom','.')
                        else: pdat['Meta'] = 'Protein'
                        pdat['%s_SeqNum' % analysis] = table.data()[table.dict['Index']['Hub'][hub][0]]['SeqNum']
                        if len(table.dict['Index']['Hub'][hub]) == 1 and table.data()[table.dict['Index']['Hub'][hub][0]]['Pattern'] in ['-','!']:
                            pdat['%s_SLiM' % analysis] = 0
                            pdat['%s_TopRank' % analysis] = table.data()[table.dict['Index']['Hub'][hub][0]]['Pattern']
                            continue
                        pdat['%s_SLiM' % analysis] = len(table.dict['Index']['Hub'][hub])
                        pdat['Total'] += pdat['%s_SLiM' % analysis]
                        pdat['%s_Clouds' % analysis] = len(splits['Clouds'][split].dict['Index']['Hub'][hub])
                        pdat['%s_TopRank' % analysis] = splits['TopRank'][split].data()[hub]['Pattern']
                        pdat['%s_Sig' % analysis] = splits['TopRank'][split].data()[hub]['Sig']
                        pdat['%s_FDR' % analysis] = splits['TopRank'][split].data()[hub]['FDR']
                        if analysis[:3] != 'ned': pdat['MinFDR'] = min(pdat['MinFDR'],pdat['%s_FDR' % analysis])
                    except: self.errorLog('Problem with datasetInfo(%s - %s)' % (hub,analysis))
                rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            self.printLog('\r#HUB','Summary for %s Hubs output to %s' % (rje.integerString(stot),patfile))
        except: self.errorLog(rje_zen.Zen().wisdom())   #; raise   # Delete this if method error not terrible
#########################################################################################################################
    def cloudInfo(self):  ### Outputs information on real datasets
        '''Outputs information on real clouds.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = self.db().getTable('Real')
            rdb.index('DCloud')
            patfile = '%s.clouds.tdt' % self.info['Basefile']
            pathead = ['Cloud','Meta']
            splits = self.dict['Splits']
            for analysis in splits['Main']:
                if analysis[:4] in ['rupc','rseq']: continue
                if analysis[-2:] == 'go': continue
                splits['Main'][analysis].index('DCloud')
                splits['Clouds'][analysis].newKey('DCloud')
                if analysis[-3:] == 'dom' or analysis[-6:-3] == 'dom': continue
                pathead += ['%s_Pattern' % analysis,'%s_Sig' % analysis,'%s_FDR' % analysis]
            rje.delimitedFileOutput(self,patfile,pathead,rje_backup=True)
            ### ~ [2] Summarise Pattern Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,len(rdb.dict['Index']['DCloud']))
            for cloud in rje.sortKeys(rdb.dict['Index']['DCloud']):
                self.progLog('\r#CLOUD','Summarising %s Clouds: %.2f%%' % (rje.integerString(stot),sx/stot)); sx += 100.0
                pdat = {}   
                for analysis in splits['Main'].keys()[0:]:
                    if analysis[:4] in ['rupc','rseq']: continue
                    if analysis[-2:] == 'go': continue
                    table = splits['Clouds'][analysis]
                    if not table.data().has_key(cloud): continue
                    if analysis[-3:] == 'dom': pdat['Meta'] = 'Domain'; analysis = analysis[:-3]
                    elif analysis[-6:-3] == 'dom': pdat['Meta'] = 'Domain'; analysis = string.replace(analysis,'.dom','.')
                    else: pdat['Meta'] = 'Protein'
                    entry = table.data()[cloud]
                    for x in ['Pattern','Sig','FDR']: pdat['%s_%s' % (analysis,x)] = entry[x]
                if not pdat: continue
                pdat['Cloud'] = cloud
                rje.delimitedFileOutput(self,patfile,pathead,datadict=pdat)
            self.printLog('\r#CLOUD','Summary for %s Clouds output to %s' % (rje.integerString(stot),patfile))
        except: self.errorLog(rje_zen.Zen().wisdom())#; raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <5> ### Data reformatting Methods                                                                               #
#########################################################################################################################
    def slimJIM(self):  ### Make tables for SLiMJIM processing and then generate output using SLiMJIM
        '''Make tables for SLiMJIM processing and then generate output using SLiMJIM.'''
        try:### ~ [0] Setup and pre-processing information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = {}    # Dictionary of {'Hub':{'Spoke':[]}}
            # Table 'PPI' = ['Hub','Spoke'] : Hub	HubUni	Spoke	SpokeUni	SpokeSeq	Evidence
            ppidb = self.db().getTable('PPI')
            ppx = 0
            for entry in ppidb.entries():
                hub = entry['Hub']; spoke = entry['Spoke']; evidence = string.split(entry['Evidence'],'|')
                if hub not in data: data[hub] = {'ppi':{},'type':'gene'}
                data[hub]['ppi'][spoke] = evidence; ppx += 1
                self.progLog('\r#PPI','%s PPI for %s hubs' % (rje.integerString(ppx),rje.integerString(len(data))))
            self.printLog('\r#PPI','%s PPI for %s hubs' % (rje.integerString(ppx),rje.integerString(len(data))))
            # Table 'SeqMap' = ['SpokeSeq'] : ['SpokeUni','Spoke']
            mapdb = self.db().getTable('SeqMap')
            # Table 'Domains' =  ['Type,'Name,'Start,'End'] : Type	Name	Start	End	Eval	Score
            domains = {}    # Dictionary of {'Domain':[hubs]}
            domdb = self.db().getTable('Domains')
            ddx = 0
            for entry in domdb.entries():
                dom = entry['Type']; seq = entry['Name']
                try: genelist = mapdb.data()[seq]['Spoke']
                except: continue    # This sequence is not in the PPI data
                if dom not in domains: domains[dom] = []
                for gene in genelist:
                    if gene not in domains[dom]: domains[dom].append(gene)
                    if gene not in data: data[gene] = {'ppi':{},'type':'gene'}
                    if 'domains' not in data[gene]: data[gene]['domains'] = []
                    if dom not in data[gene]['domains']: data[gene]['domains'].append(dom); ddx += 1
                self.progLog('\r#DOM','%s examples of %s domains' % (rje.integerString(ddx),rje.integerString(len(domains))))
            self.printLog('\r#DOM','%s examples of %s domains' % (rje.integerString(ddx),rje.integerString(len(domains))))
            genelist = rje.sortKeys(data)
            domlist = rje.sortKeys(domains)

            ### ~ [1] Generate tables of datasets {hub:{type:[spokes]}} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Genes as hubs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (gx,gtot) = (0.0,len(genelist)); bx = 0
            for gene in genelist:
                self.progLog('\r#GENE','Genes as hubs %.2f%%' % (gx/gtot)); gx += 100.0
                for type in ['y2h','bin','com']: data[gene][type] = {}
                for spoke in data[gene]['ppi']:
                    p = {'y2h':False,'bin':True,'com':False,'bum':True}
                    for e in data[gene]['ppi'][spoke]:
                        if string.split(e,':')[1] in self.list['Direct']: p['y2h'] = True
                        if string.split(e,':')[1] in self.list['Indirect']: p['com'] = True
                    if not p['y2h']:
                        for i in data[gene]['ppi']:
                            if i != spoke and i in data[spoke]['ppi']: p['bum'] = False
                            if i not in [spoke,gene] and i in data[spoke]['ppi']: p['bin'] = False
                    for type in ['y2h','bin','com']:
                        if p[type]: data[gene][type][spoke] = data[gene]['ppi'][spoke]
                    #if p['bin'] and not p['bum']: self.printLog('#BUM','%s >> %s' % (gene,spoke),screen=False); bx += 1
            self.printLog('\r#GENE','%s genes as hubs' % rje.integerString(gtot))
            #self.printLog('#BUM','%s potential problem binary PPI messed up' % rje.integerString(bx))
            ## ~ [1b] Domains as hubs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (gx,gtot) = (0.0,len(domlist))
            for dom in domlist:
                self.progLog('\r#DOM','Domains as hubs %.1f%%' % (gx/gtot)); gx += 100.0
                data[dom] = {'y2h':{},'bin':{},'com':{},'ppi':{},'type':'domain'}
                for hub in domains[dom]:
                    for type in ['y2h','bin','com','ppi']:
                        for spoke in data[hub][type]:
                            if spoke not in data[dom][type]: data[dom][type][spoke] = ['Domain']
            self.printLog('\r#DOM','%s domains as hubs.' % rje.integerString(gtot))
            if 'Filament' not in data: self.deBug(rje.sortKeys(data))
            # >> Currently have all genes and domains as keys in data dictionary linking to different ppi << #
            # >> {gene:{type:{spokegene:[evidence]}}} or {domain:{type:{spokegenes:[]}}} << #
            #!# Old SLiMJIM MakeTables should be replaced with code here once new website needs are defined #!#

            ## ~ [0a] ~ General Data dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ##!## OLD SLiMJIM Dictionaries ##!##
            #self.dict['SLiMFinder'] = {}    # *All* results for specified RunIDs
            #self.dict['PPI'] = {}           # Dictionary of {HubGene:{SpokeGene:{Data}}}
            #self.dict['SeqMap'] = {}        # Dictionary of {SpokeSeq:[Genes]}
            ## ~ [0b] ~ Results data lists and dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.list['Hubs'] = []          # List of all hub genes with datasets for processing
            #self.dict['HubSlims'] = {}      # Dictionary of {Hub:[SlimList]}
            #self.dict['SpokeSlims'] = {}    # Dictionary of {SpokeSeq:[SlimList]}

            ### ~ [2] Setup SLiMJIM Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sj = self.obj['SLiMJIM'] = slimjim.SLiMJIM(self.log,self.cmd_list)
            sj.obj['GeneMap'] = self.obj['GeneMap']
            self.obj['GeneMap'].run()   #?# Is this necessary? (From SLiMJIM)
            sj.dict['Data'] = data      # Replace some of SLiMJIM dictionaries
            sj.obj['DB'] = self.db()    # Replace some of SLiMJIM dictionaries
            ## ~ [2a] Read GO Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sj.obj['GO'] = self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            self.obj['GO'].readGO()
            self.obj['GO'].mapEnsGO()
            sj.dict['GO'] = self.obj['GO'].dict['EnsGO']
            ## ~ [2b] ~ Results data lists and dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sj.list['RunTypes'] = ['ppi','y2h','bin','com']
            sj.list['Spokes'] = genelist    # List of all spoke genes contained in one or more PPI datasets
            sj.dict['EnsLoci'] = {}         # Dictionary of {SpokeGene:SpokeSeq}
            for seq in mapdb.data():
                for gene in mapdb.data()[seq]['Spoke']:
                    sj.dict['EnsLoci'][gene] = seq
                    if gene in data: data[gene]['seq'] = seq

            ### ~ [3] Update data dictionary with extra data for SLiMJIM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dbdata().dropEntries(['DType!=real','PPISet==ned','PPISet==neddom'],inverse=False,log=True)
            fulldata = rje.combineDict({},self.dbdata().data())
            sj.list['Hubs'] = rje.sortKeys(self.dbdata().index('Hub'))  # Processed datasets
            self.dbdata().dropEntries(['Rank==0'],inverse=False,log=True)
            odb = self.db().getTable('Occ')
            ## ~ [3a] SlimX counts of returned SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for hub in data:
                if hub in self.dbdata().index('Hub'): data[hub]['slimx'] = len(self.dbdata().index('Hub')[hub])
                else: data[hub]['slimx'] = 0
                try: seq = sj.dict['EnsLoci'][hub]
                except: continue
                if seq in odb.index('Seq'): data[hub]['spokeslimx'] = len(odb.index('Seq')[seq])
                else: data[hub]['spokeslimx'] = 0
            self.dbdata().dict['Data'] = fulldata
            self.dbdata().dict['Index'] = {}
            ## ~ [3b] SLiMHubs and SLiMSpokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sj.dict['SlimHubs'] = {}      # Dictionary of {Pattern:{Hub:[SpokeSeq]}}
            sj.dict['SlimSpokes'] = {}    # Dictionary of {Pattern:{SpokeSeq:[Hubs]}}
            for pattern in rje.sortKeys(odb.index('Pattern')):
                sj.dict['SlimHubs'][pattern] = {'ALL':{}}
                sj.dict['SlimSpokes'][pattern] = {}
                for entry in odb.entryList(odb.index('Pattern')[pattern]):
                    ppitype = string.split(entry['Dataset'],'.')[-1][:3]
                    if ppitype not in sj.dict['SlimHubs'][pattern]: sj.dict['SlimHubs'][pattern][ppitype] = {}
                    if entry['Hub'] not in sj.dict['SlimHubs'][pattern][ppitype]: sj.dict['SlimHubs'][pattern][ppitype][entry['Hub']] = rje_slim.expectString(float(entry['Sig']))
                    elif float(entry['Sig']) < sj.dict['SlimHubs'][pattern][ppitype][entry['Hub']]: sj.dict['SlimHubs'][pattern][ppitype][entry['Hub']] = rje_slim.expectString(float(entry['Sig']))
                    if entry['Hub'] not in sj.dict['SlimHubs'][pattern]['ALL']: sj.dict['SlimHubs'][pattern]['ALL'][entry['Hub']] = rje_slim.expectString(float(entry['Sig']))
                    elif float(entry['Sig']) < sj.dict['SlimHubs'][pattern]['ALL'][entry['Hub']]: sj.dict['SlimHubs'][pattern]['ALL'][entry['Hub']] = rje_slim.expectString(float(entry['Sig']))
                    if entry['Seq'] not in sj.dict['SlimSpokes'][pattern]: sj.dict['SlimSpokes'][pattern][entry['Seq']] = []
                    if entry['Hub'] not in sj.dict['SlimSpokes'][pattern][entry['Seq']]: sj.dict['SlimSpokes'][pattern][entry['Seq']].append(entry['Hub']) 
            
            ### ~ [4] Run modified SLiMJIM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sj.opt['MakeJIM']: sj.makeJIM()    #!# Not yet modified #!#
            if sj.opt['MakeHTML']: sj.makeHTML()  #!# Not yet modified #!#

            # Compile tables of required information
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: SLiMFRAP Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
class GeneMap(rje_genemap.GeneMap):
    pass
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
    try:SLiMFRAP(mainlog,cmd_list).run()

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
