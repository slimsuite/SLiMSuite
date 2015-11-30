#!/usr/bin/python

# See below for name and description
# Copyright (C) 2010 Richard J. Edwards <software@cabbagesofdoom.co.uk>
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
# Author contact: <software@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       RJE_HM_HTML
Description:  Module for generating HTML for APHID enrichment study
Version:      0.3
Last Edit:    04/11/10
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module makes a set of linked webpages for the H&M integrin tail analysis. It is based on, and heavily uses,
    the code for rje_slimhtml.

Commandline:
    datapath=PATH       : Path to parent data directory [./]
    htmlpath=PATH       : Path of parent html directory [./aa_ff]
    stylesheets=LIST    : List of CSS files to use ['../example.css','../redwards.css']
    border=X            : Border setting for tables [0]
    dropfields=LIST     : Fields to exclude from summary tables []
    badgenes=LIST       : Genes to be moved to Lysate dataset []
    fakehtml=T/F        : Whether to make UniFake HTML [True]
    unifake=PATH        : Path to UniFake dat file(s) [./unifake/]
    unipath=PATH        : Path to real UniProt dat file(s) [./uniprot/]
    unireal=LIST        : Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]
    makepng=T/F         : Whether to (look for and) make PNG files with R [True]
    makepages=LIST      : Types of pages to make [front,gene,nested,interactome,trees,complex]
    titletext=FILE      : File containing (Page,ID,Title) [titletext.tdt]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_slim, rje_uniprot, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_go, rje_slim, rje_zen, rje_slimhtml, rje_ppi
import rje_dismatrix_V2, rje_tree
import slimfinder
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation. Based on rje_slimhtml.
    # 0.1 - Added Output for additional Pingu Analysis.
    # 0.2 - Added protein complex output.
    # 0.3 - Modified output for HAPPI analysis.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Incorporate new HAPPI code.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_HM_HTML', '0.3', 'November 2010', '2010')
    description = 'RJE H&M APHID HTML Module'
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
### SECTION II: SLiMHTML Class                                                                                          #
#########################################################################################################################
class HMHTML(rje_slimhtml.SLiMHTML):     
    '''
    HMHTML Class. Author: Rich Edwards (2010).

    Info:str
    - DataPath = Path to parent data directory [./]
    - HTMLPath = Path of parent html directory [./aa_ff]
    - TitleText = File containing (Page,ID,Title) [titletext.tdt]
    - UniFake = Path to UniFake dat file(s) [./unifake/]
    - UniPath = Path to real UniProt dat file(s) [./uniprot/]
    
    Opt:boolean
    - FakeHTML = Whether to make UniFake HTML [True]
    - MakePNG = Whether to (look for and) make PNG files with R [True]

    Stat:numeric
    - Border = Border setting for tables [0]

    List:list
    - DropFields = Fields to exclude from summary tables []
    - MakePages = Types of pages to make [front,gene,nested,trees]
    - UniReal = Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['DataPath','HTMLPath','TitleText']
        self.optlist = ['FakeHTML','MakePNG']
        self.statlist = ['Border']
        self.listlist = ['StyleSheets','UniReal','MakePages','BadGenes']
        self.dictlist = ['GeneMap']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'DataPath':rje.makePath('/home/re1u06/researchfiles/SBSBINF/Databases/DBase_100505/'),
                      'HTMLPath':rje.makePath('./aa_ff/'),'Basefile':'hm_html',
                      'UniFake':rje.makePath('./HumSF09_UniFake/'),'UniPath':rje.makePath('./HumSF09_UniProt/'),
                      'TitleText':'titletext.tdt'})
        self.setStat({'Border':2})
        self.setOpt({'FakeHTML':True,'MakePNG':True,'Iridis':False})
        self.list['StyleSheets'] = ['../example.css','../slimhtml.css']
        self.list['DropFields'] = []
        self.list['MakePages'] = string.split('front,gene,nested,trees,interactome,complex',',')
        self.list['UniReal'] = string.split('AC,GN,RC,RX,CC,DR,PE,KW,FT',',')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'file',['TitleText'])
                self._cmdReadList(cmd,'path',['DataPath','HTMLPath'])
                self._cmdReadList(cmd,'int',['Border'])
                self._cmdReadList(cmd,'opt',['FakeHTML','MakePNG'])
                self._cmdReadList(cmd,'list',['StyleSheets','DropFields','UniReal','MakePages','BadGenes'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            hpath = self.info['HTMLPath']
            genelist = self.list['Genes']
            ### ~ [2] ~ Make Front Page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Front page = tab per AA/FF/AA-FF/LYS/ITGA2B/ITGAL/ITGA5
            # - Summary results tab
            # - Tree tab
            # - SLiMFinder tab?
            if 'front' in self.list['MakePages']:
                html = self.frontPage(genelist)
                hfile = hpath + 'index.htm'
                open(hfile,'w').write(html['front'])
                self.printLog('#HTML','Made main %s page.' % hfile)
            ### ~ [3] ~ Make Gene Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Gene page = Details + GO (Add one row table of results)
            # - PPI (Just other genes in sample) = Spoke, Evidence, Description, Class, AA, FF
            # - Nested interactome: ITGA2B, AA, FF (inside), spokes in genelist outside
            if 'gene' in self.list['MakePages']:
                gx = 0.0; gtot = len(genelist)
                for gene in genelist:
                    self.progLog('\r#GENE','Making Gene Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.genePage(gene)
                self.printLog('\r#GENE','%s Gene HTML Pages made.' % rje.integerString(gtot))
            ## ~ [3a] UniFake Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if 'fake' in self.list['MakePages']: self.uniFakeHTML()     #??#
            ## ~ [3b] GO Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if 'go' in self.list['MakePages']:
            #    gdb = self.db().getTable('GO')
            #    gx = 0.0; gtot = len(gdb.index('GO_Desc'))
            #    for go in rje.sortKeys(gdb.index('GO_Desc')):
            #        self.progLog('\r#GO','Making GO Pages: %.2f%%' % (gx/gtot)); gx += 100.0
            #        goid = gdb.dataList(gdb.indexEntries('GO_Desc',go),'GO_ID')[0]
            #        hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
            #        if not self.opt['Force'] and self.checkHTML(hfile): continue
            #        self.goPage(go,goid)
            #    self.printLog('\r#GO','%s GO HTML Pages made.' % rje.integerString(gtot))
            ### ~ [4] ~ Make Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['MakePNG']:
                if 'nested' in self.list['MakePages']:
                    for gene in genelist: self.nestedPNG(gene)
                if 'complex' in self.list['MakePages']:
                    for gene in genelist: self.complexPNG(gene)
                if 'interactome' in self.list['MakePages']:
                    for gene in genelist: self.interactomePNG(gene)
                    for pclass in self.db().getTable('Main').index('Class'):  self.interactomePNG(pclass)
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dpath = rje.makePath(self.info['DataPath'])
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            ### ~ [2] General Database Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] PPI table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = dpath + 'Pingu/pingu.100505.pairwise.tdt'   #!# Update this #!# (Pingu?)
            ppi = self.db().addTable(dfile,['Hub','Spoke'],'All',name='PPI')
            #Hub     HubUni  Spoke   SpokeUni        SpokeSeq        Evidence        ppi     bin     com     y2h
            #1082356 -       GSTK1   Q9Y2Q3  GSTK1_HUMAN__Q9Y2Q3     IntAct:anti bait coip   Y       Y       -       -
            ## ~ [2b] XRef table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = dpath + 'HGNC/genemap_100505.data.tdt'
            dbxref = self.db().addTable(dfile,['Gene'],'All',name='DBXRef')
            #Gene    Entrez  HPRD    OMIM    UniProt EnsEMBL EnsLoci EnsDesc ppi     bin     com     y2h
            #A1BG    1       00726   138670  P04217  ENSG00000121410 A1BG_HUMAN__ENSP00000263100     Alpha-1B-glycoprotein Precursor (Alpha-1-B glycoprotein)        1       1       0       0
            ## ~ [2c] GO table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gofile = 'test.go.tdt'
            if not os.path.exists(gofile): rje_go.GO(self.log,self.cmd_list).test()
            godb = self.db().addTable(gofile,['EnsG','GO_ID'],'All',name='GO')
            self.dict['GO'] = {}
            for gkey in godb.data():
                (ensg,go) = string.split(gkey)
                gtype = godb.data()[gkey]['GO_Type'].upper()
                if ensg not in self.dict['GO']: self.dict['GO'][ensg] = {}
                if gtype not in self.dict['GO'][ensg]: self.dict['GO'][ensg][gtype] = []
                self.dict['GO'][ensg][gtype].append((go,godb.data()[gkey]['GO_Desc']))
            #EnsG    GO_ID   GO_Type GO_Desc
            #ENSG00000000003 0007186 bp      G-protein coupled receptor protein signaling pathway
            ## ~ [2d] PFam table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #dfile = string.replace(dpath,'MainData','UniFake') + 'ens_HUMAN.unifake.pfam.tdt'   #?#
            #pfam = self.db().addTable(dfile,['Type','Name','Start','End'],'All',name='PFam')
            #pfam.index('Type'); pfam.index('Name')
            #Type    Name    Start   End     Eval    Score
            #Dysbindin       DBND1_HUMAN__ENSP00000002501    14      153     2.30e-92        317.8
            #http://pfam.sanger.ac.uk/family/Dysbindin

            ### ~ [3] APHID Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NRID	Description	pi.RA.v.RL	pi.AA.v.AL	pi.RF.v.RL	pi.AF.v.AL	pi.RA.v.XL	pi.RF.v.XL	pj.RA.v.RL	pj.AA.v.AL	pj.RF.v.RL	pj.AF.v.AL	pj.RA.v.XL	pj.RF.v.XL
            # Identifier	HGNC	EnsG	Rep|AA	IntNorm|AA	Rep|AF	IntNorm|AF	Rep|AL	IntNorm|AL	Rep|RF	IntNorm|RF	Rep|RA	IntNorm|RA	Rep|RL	IntNorm|RL	Rep|XL	IntNorm|XL
            # AA.v.AL	AF.v.AL	RF.v.RL	RA.v.RL	RF.v.XL	RA.v.XL
            intfield = string.split('Rep|AA	 Rep|AF	Rep|AL Rep|RF Rep|RA Rep|RL Rep|XL')
            numfield = string.split('pi.RA.v.RL pi.AA.v.AL	pi.RF.v.RL	pi.AF.v.AL	pi.RA.v.XL	pi.RF.v.XL	pj.RA.v.RL	pj.AA.v.AL	pj.RF.v.RL	pj.AF.v.AL	pj.RA.v.XL	pj.RF.v.XL')
            numfield += string.split('IntNorm|AA IntNorm|AF IntNorm|AL IntNorm|RF IntNorm|RA IntNorm|RL IntNorm|XL AA.v.AL	AF.v.AL	RF.v.RL	RA.v.RL	RF.v.XL	RA.v.XL')
            dfile = '%s.bootstrap.tdt' % self.info['Basefile']
            main = self.db().addTable(dfile,['NRID'],'All',name='Main')            
            ## ~ [3a] Classify Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            main.addField('PPI',after='NRID'); main.addField('Class',after='NRID')
            main.addField('enrFF',after='NRID'); main.addField('FF',after='NRID'); main.addField('enrAA',after='NRID'); main.addField('AA',after='NRID');
            main.addField('Entrez',after='HGNC')
            pi = 0.5
            for entry in main.entries():
                try: entry['Entrez'] = dbxref.data()[entry['NRID']]['Entrez']
                except: entry['Entrez'] = '-'
                for f in numfield:
                    try: entry[f] = string.atof(entry[f])
                    except: pass
                for f in intfield:
                    try: entry[f] = string.atoi(entry[f])
                    except: pass
                entry['Class'] = ''
                entry['Identifier'] = string.replace(entry['Identifier'],';','; ')
                # ~ Convert to single AA and FF values ~ #
                if entry['pi.AA.v.AL'] != 'n/a' and entry['pi.RA.v.RL'] != 'n/a': entry['AA'] = max(entry['pi.AA.v.AL'],entry['pi.RA.v.RL'])
                elif entry['pi.AA.v.AL'] != 'n/a': entry['AA'] = entry['pi.AA.v.AL']
                elif entry['pi.RA.v.RL'] != 'n/a': entry['AA'] = entry['pi.RA.v.RL']
                else: entry['AA'] = 0.0
                if entry['AA.v.AL'] and entry['RA.v.RL']: entry['enrAA'] = max(entry['AA.v.AL'],entry['RA.v.RL'])
                elif entry['AA.v.AL']: entry['enrAA'] = entry['AA.v.AL']
                elif entry['RA.v.RL']: entry['enrAA'] = entry['RA.v.RL']
                else: entry['enrAA'] = 0.0
                if entry['pi.AF.v.AL'] != 'n/a' and entry['pi.RF.v.RL'] != 'n/a': entry['FF'] = max(entry['pi.AF.v.AL'],entry['pi.RF.v.RL'])
                elif entry['pi.AF.v.AL'] != 'n/a': entry['FF'] = entry['pi.AF.v.AL']
                elif entry['pi.RF.v.RL'] != 'n/a': entry['FF'] = entry['pi.RF.v.RL']
                else: entry['FF'] = 0.0
                if entry['AF.v.AL'] and entry['RF.v.RL']: entry['enrFF'] = max(entry['AF.v.AL'],entry['RF.v.RL'])
                elif entry['AF.v.AL']: entry['enrFF'] = entry['AF.v.AL']
                elif entry['RF.v.RL']: entry['enrFF'] = entry['RF.v.RL']
                else: entry['enrFF'] = 0.0
                # ~ Classify ~ #
                if entry['NRID'] in self.list['BadGenes']: entry['Class'] = 'LYS'
                elif entry['AA'] > pi and entry['enrAA'] > 1.0:
                    entry['Class'] = 'AA'
                    if entry['FF'] > pi and entry['enrFF'] > 1.0: entry['Class'] = 'AA-FF'
                elif entry['FF'] > pi and entry['enrFF'] > 1.0: entry['Class'] = 'FF'
                else: entry['Class'] = 'LYS'
                #!# Add PPI count if not already in table #!#
                entry['PPI'] = 0
                if entry['NRID'] in ppi.index('Hub'): entry['PPI'] = len(ppi.indexEntries('Hub',entry['NRID']))
            ## ~ [3b] Add ITGA2B interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for addgene in ['ITGA2B','ITGA5','ITGAL']:
                for gene in ppi.dataList(ppi.indexEntries('Hub',addgene),'Spoke'):
                    igene = self.gene(gene)
                    if gene in main.index('NRID') or igene in main.index('NRID'): continue
                    gentry = {'NRID':gene,'Identifier':igene,'HGNC':'','Desc':igene,'EnsG':'-','Entrez':'-'}
                    if igene in dbxref.data():
                        gentry['HGNC'] = dbxref.data()[igene]['Symbol']
                        gentry['Desc'] = dbxref.data()[igene]['EnsDesc']
                        gentry['EnsG'] = dbxref.data()[igene]['EnsEMBL']
                        gentry['Entrez'] = dbxref.data()[igene]['Entrez']
                    for f in numfield + intfield: gentry[f] = 'n/a'
                    gentry['AA'] = gentry['FF'] = gentry['enrAA'] = gentry['enrFF'] = 0.0
                    gentry['Class'] = 'Integrin'    #addgene
                    gentry['PPI'] = len(ppi.indexEntries('Hub',gene))
                    main.dict['Data'][main.makeKey(gentry)] = gentry
                    main.dict['Index']['NRID'][gene] = [main.makeKey(gentry)]
            self.list['Genes'] = rje.sortKeys(main.index('NRID',force=True))
            main.saveToFile('%s.full.tdt' % self.info['Basefile'])
            ## ~ [3c] Add as new PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppi_pi = 1.0    # Add as PPI if ANY enrichment.
            for entry in main.entries():
                if entry['enrAA'] > ppi_pi:
                    try:
                        pentry = {'Hub':'AA','Spoke':entry['NRID'],'HubUni':'','SpokeUni':'','SpokeSeq':'-','Evidence':'Rest:%s|Act:%s' % (str(entry['RA.v.RL'])[:5],str(entry['AA.v.AL'])[:5])}
                        ppi.dict['Data'][ppi.makeKey(pentry)] = pentry
                        pentry = {'Spoke':'AA','Hub':entry['NRID'],'HubUni':'','SpokeUni':'','SpokeSeq':'-','Evidence':'Rest:%s|Act:%s' % (str(entry['RA.v.RL'])[:5],str(entry['AA.v.AL'])[:5])}
                        ppi.dict['Data'][ppi.makeKey(pentry)] = pentry
                    except: self.deBug(entry)
                if entry['enrFF'] > ppi_pi:
                    pentry = {'Hub':'FF','Spoke':entry['NRID'],'HubUni':'','SpokeUni':'','SpokeSeq':'-','Evidence':'Rest:%s|Act:%s' % (str(entry['RF.v.RL'])[:5],str(entry['AF.v.AL'])[:5])}
                    ppi.dict['Data'][ppi.makeKey(pentry)] = pentry
                    pentry = {'Spoke':'FF','Hub':entry['NRID'],'HubUni':'','SpokeUni':'','SpokeSeq':'-','Evidence':'Rest:%s|Act:%s' % (str(entry['RF.v.RL'])[:5],str(entry['AF.v.AL'])[:5])}
                    ppi.dict['Data'][ppi.makeKey(pentry)] = pentry
            ppi.dict['Index'] = {}
            #self.deBug('AA:%s' % len(ppi.indexEntries('Hub','AA')))
            #self.deBug('FF:%s' % len(ppi.indexEntries('Hub','FF')))
            for field in string.split('IntNorm|AA IntNorm|AF IntNorm|AL IntNorm|RF IntNorm|RA IntNorm|RL IntNorm|XL'): main.renameField(field,string.replace(field,'IntNorm','Int'))
            sampdb = self.db().copyTable(main,'Samples')
            sampdb.addField('Sample'); sampdb.addField('Source'); sampdb.addField('AvF',after='enrFF'); sampdb.addField('MeanEnr')
            for entry in sampdb.entries():
                entry['Sample'] = entry['Class']
                entry['Source'] = entry['Identifier']
                entry['Identifier'] = entry['NRID']
                entry['AvF'] = entry['AA'] - entry['FF']
                entry['MeanEnr'] = []
                for field in ['AA.v.AL','AF.v.AL','RF.v.RL','RA.v.RL']:
                    try: entry['MeanEnr'].append(float(entry[field]))
                    except: entry['MeanEnr'].append(1.0)
                entry['MeanEnr'] = rje.geoMean(entry['MeanEnr'])
                for field in sampdb.fields():
                    if field not in ['Sample','Class','NRID','Class','PPI','Identifier','HGNC','Entrez','EnsG','Desc'] and field[:3] != 'Rep':
                        try: entry[field] = rje.expectString(float(entry[field]))
                        except: pass
            sampdb.saveToFile('%s.samples.tdt' % self.info['Basefile']) 
            
            ## ~ [3d] Raw Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = '%s.data_nr_map.tdt' % self.info['Basefile']
            raw = self.db().addTable(dfile,['Expt','Subpop','Slice','Identifier'],'All',name='Raw')
            
            ### ~ [4] GABLAM Distance Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dfile = '%s.gablam.tdt' % self.info['Basefile']
            self.obj['Dis'] = dis = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
            dis.loadFromDataTable(dfile)
            dis.forceSymmetry(method='mean')
            ## ~ [4a] Convert and add ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newmatrix = {}
            for obj1 in rje.sortKeys(dis.dict['Matrix']):
                newmatrix[self.gene(obj1)] = {}
                for obj2 in dis.dict['Matrix'][obj1]: newmatrix[self.gene(obj1)][self.gene(obj2)] = dis.dict['Matrix'][obj1][obj2]
            dis.dict['Matrix'] = newmatrix
            for gene in self.list['Genes']:
                if gene not in dis.dict['Matrix']: dis.dict['Matrix'][gene] = {}
            self.printLog('#DIS','GABLAM matrix loaded and converted to Gene Identifiers')

            ### ~ [5] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['HTMLPath'])
            subdirlist = ['gene','nested_png','interactome_png']
            for subdir in subdirlist: rje.mkDir(self,rje.makePath(self.info['HTMLPath']+subdir))
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Results pages                                                                                           #
#########################################################################################################################
    def seqDetailsHTML(self,gene,dbxref):    #gene,seqid,dbxref,desc,godata):  ### Returns HTML text for seq details table
        '''
        Returns HTML text for seq details table.
        >> gene:str = Gene symbol
        >> seqid:str = Sequence Identifier
        >> dbxref:dict = Dictionary of {db:id} for GeneCards, EBI, EnsEMBL, HPRD, OMIM
        >> desc:str = Sequence description
        >> godata:dict = {CC/BP/MF:[(id,name)] list}
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqid = '-'; desc = '<i>No sequence mapping.</i>'; godata = {}
            if dbxref:      #gene in self.data('DBXRef'):
                #dbxref = self.data('DBXRef')[gene]
                seqid = dbxref['EnsLoci']
                desc = dbxref['EnsDesc']
                ensg = dbxref['EnsEMBL']
                try: godata = self.dict['GO'][ensg]
                except: pass
            ### ~ [2] ~ Gene Name and Links Out ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html = ['','<!-- ~~~~ %s Gene Summary details table ~~~~ -->' % gene,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="30%%"><a href="../gene/%s.html">%s<b>%s</b></FONT></a></td>' % (gene,gfont,gene),
                    '<td width="50%%">%s<b>%s</b></FONT></td>' % (ifont,seqid),
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>']
            #x#if 'Gene' not in dbxref: dbxref['Gene'] = gene
            html += ['<tr><td colspan=2>']
            for db in ['Gene','UniProt','EnsEMBL','HPRD','OMIM']:
                if db not in dbxref: continue
                if db == 'Gene':
                    href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' % dbxref[db]
                    alt = 'GeneCards %s' % dbxref[db]
                    src = '../resources/genecard.png'
                    title = self.titleText('xref','gene')
                if db == 'UniProt':
                    href = 'http://www.uniprot.org/uniprot/%s' % dbxref[db]
                    alt = 'UniProt %s' % dbxref[db]
                    src = '../resources/logo_ebi.png'
                    title = self.titleText('xref','uniprot')
                if db == 'EnsEMBL':
                    href = 'http://www.ensembl.org/Homo_sapiens/geneview?gene=%s' % dbxref[db]
                    alt = 'EnsEMBL %s' % dbxref[db]
                    src = '../resources/e-bang.gif'
                    title = self.titleText('xref','ensembl')
                if db == 'HPRD':
                    #href = 'http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1' % (dbxref[db],dbxref[db])
                    href = 'http://www.hprd.org/summary?hprd_id=%s&isoform_id=%s_1&isoform_name=Isoform_1' % (dbxref[db],dbxref[db])
                    alt = 'HPRD %s' % dbxref[db]
                    src = '../resources/hprd.png'
                    title = self.titleText('xref','hprd')
                if db == 'OMIM':
                    href = 'http://www.ncbi.nlm.nih.gov/omim/%s' % dbxref[db]
                    alt = 'OMIM %s' % dbxref[db]
                    src = '../resources/omim.png'
                    title = self.titleText('xref','omim')
                html += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)] 
            html += ['</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,desc)]
            html += ['</table>','<!-- ~~~~ End %s Gene Summary details table ~~~~ -->' % gene,'']
            ### ~ [3] ~ Summary Data Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NRID	Description	pi.RA.v.RL	pi.AA.v.AL	pi.RF.v.RL	pi.AF.v.AL	pi.RA.v.XL	pi.RF.v.XL	pj.RA.v.RL	pj.AA.v.AL	pj.RF.v.RL	pj.AF.v.AL	pj.RA.v.XL	pj.RF.v.XL
            # Identifier	HGNC	EnsG	Rep|AA	IntNorm|AA	Rep|AF	IntNorm|AF	Rep|AL	IntNorm|AL	Rep|RF	IntNorm|RF	Rep|RA	IntNorm|RA	Rep|RL	IntNorm|RL	Rep|XL	IntNorm|XL
            # AA.v.AL	AF.v.AL	RF.v.RL	RA.v.RL	RF.v.XL	RA.v.XL
            main = self.db().getTable('Main')
            ## ~ [3a] ~ Class and Mapping Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            droplist = []
            for field in main.fields():
                if field not in ['NRID','Identifier','HGNC','EnsG','Desc','AA','enrAA','FF','enrFF','Class','PPI']: droplist.append(field)
            fwidth = {'NRID':120,'Identifier':150,'HGNC':90,'EnsG':150,'Desc':250,'AA':50,'FF':50,'enrAA':50,'enrFF':50,'Class':50,'PPI':50}
            html += [self.resultTableHTML(main.fields(),main.subset('NRID',gene),drop=droplist,width=fwidth)]
            ## ~ [3b] ~ Data summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            summary = []
            summary.append(string.split('Rep|AA	Int|AA AA.v.AL pi.AA.v.AL || Rep|RA Int|RA RA.v.RL pi.RA.v.RL'))
            summary.append(string.split('Rep|AF	Int|AF AF.v.AL pi.AF.v.AL || Rep|RF	Int|RF RF.v.RL pi.RF.v.RL'))
            summary.append(string.split('Rep|AL	Int|AL	Rep|RL	Int|RL'))
            for thead in summary:
                droplist = []; fwidth = {}
                for field in main.fields():
                    if field not in thead: droplist.append(field)
                    else: fwidth[field] = 80
                fwidth['||'] = 5
                html += [self.resultTableHTML(thead,main.subset('NRID',gene),drop=droplist,width=fwidth)]
            ### ~ [4] ~ GO tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gtab = []
            for gtype in ['BP','CC','MF']:
                if gtype in godata:
                    gdict = {}
                    for gotup in godata[gtype]:
                        if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                        #gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                        gdict[gotup[1]] = goLink(gotup[1],gotup[0])
                    if gdict:
                        ghtml = []
                        for g in rje.sortKeys(gdict): ghtml.append(gdict[g])
                        gtab.append(('GO_%s' % gtype,string.join(ghtml,' ~ '),self.titleText('tab','go')))
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s' % gene,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % gene,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % gene,'']
        except: self.errorLog('seqDetailsHTML Error')
        ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #print string.join(html,'\n')
        return string.join(html,'\n')
#########################################################################################################################
    def frontPage(self,genelist):    ### Generates HTML front page for analysis
        '''
        Generates HTML front page for analysis.
        >> genelist:list of analysed genes
        >> domainlist:list of analysed domains
        >> slimlist:list of returned slims
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = htmlHead('Integrin Tail Proteomics Analysis',self.list['StyleSheets'],frontpage=True)
            html += '\n<table border="0" width="100%"><tr>\n<td width="60%"><h1>Integrin Tail Proteomics Analysis</h1>\n'
            html += '<br><p><i>This data has not yet been published and should not be used without permission.</i></p></td>\n'
            html += '<td width="40%" valign="top" align="right">\n'
            html += '<a href="http://www.rcsi.ie/"><img src="./resources/rcsi.png" height="100" alt="RCSI" title="The Royal College of Surgeons in Ireland"></a>\n'
            html += '<a href="http://www.personal.soton.ac.uk/re1u06/"><img src="./resources/SBS_100.png" height="100" alt="RJE Homepage" title="RJE Homepage"></a>\n'
            html += '</td></tr></table>\n\n'
            head_html = html[0:]
            main = self.db().getTable('Main')
            ppi = self.db().getTable('PPI')
            dbxref = self.db().getTable('DBXRef')
            jtxt = ' ~ '
            ### ~ [2] ~ Class Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            droplist = []
            for field in main.fields():
                if field not in ['NRID','Identifier','HGNC','Entrez','EnsG','Desc','AA','FF','enrAA','enrFF','Class','PPI']: droplist.append(field)
            fwidth = {'NRID':90,'Identifier':100,'HGNC':90,'Entrez':90,'EnsG':120,'Desc':250,'AA':50,'FF':50,'enrAA':50,'enrFF':50,'Class':50,'PPI':50}
            classtab = []          
            for c in rje.sortKeys(main.index('Class')):
                tablist = []
                chtml = '<H2>%s</H2>\n' % self.classDef(c)
                tablist.append(('Data',string.replace(self.resultTableHTML(main.fields(),main.subset('Class',c),drop=droplist,width=fwidth),'href="../','target="_blank" href="./'),'%s Summary' % c))
                ## ~ [2b] ~ Interactome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                intpng = './interactome_png/%s.interactome.png' % c
                alt = '%s interactome image not yet made' % c
                title = 'Graphical representation of %s interactome. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % c
                tablist.append(('Interactome','<h2>%s Interactome</h2>\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (c,intpng,intpng,alt,title),title))
                ## ~ [2c] ~ Finish Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                classtab.append((c,'%s%s' % (chtml,tabberHTML('%s-tabs' % c,tablist,level=1)),self.classDef(c)))
            ### ~ [3] ~ Front page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add a little text ? #!#
            html += tabberHTML('Main',classtab,level=0)
            html += htmlTail()
        except: self.errorLog('frontPage Error')
        ### ~ [3] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return {'front':html}
#########################################################################################################################
    def classDef(self,c):   ### Returns definition of Class
        '''Returns definition of Class.'''
        classdef = {'AA':'Genes strongly enriched for AApep but not FFpep binding',
                    'FF':'Genes strongly enriched for FFpep but not AApep binding',
                    'AA-FF':'Genes strongly enriched for AApep and FFpep binding',
                    'LYS':'Genes identified in experiments but not strongly enriched for AApep or FFpep binding',
                    'ITGA2B':'Known ITGA2B interactors not returned',
                    'ITGAL':'Known ITGAL interactors not returned or interacting with ITGA2B/5',
                    'ITGA5':'Known ITGA5 interactors not returned or interacting with ITGA2B',
                    'Integrin':'Known ITGA2B/ITGAL/ITGA5 interactors not returned by pulldowns'
                    }
        try: return classdef[c]
        except: self.errorLog('ClassDef error - "%s"' % c); return 'Unknown Class'
#########################################################################################################################
    def classCol(self,c,strict=True):   ### Returns colour of Class for PNG
        '''Returns definition of Class.'''
        classcol = {'AA':15,'FF':19,'AA-FF':8,'LYS':6,'ITGA2B':3,'ITGAL':3,'ITGA5':3,'Integrin':3}
        try: return classcol[c]
        except:
            if strict: self.errorLog('ClassDef error - "%s"' % c)
            return 12
#########################################################################################################################
    def genePage(self,gene):     ### Generates HTML for gene page
        '''
        Generates HTML for gene page.
        >> gene:str = Gene for page construction
        '''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dbxref = {}; seqid = '-'; desc = '<i>No sequence mapping.</i>'; godata = {}
            if gene in self.data('DBXRef'): 
                dbxref = self.data('DBXRef')[gene]
                seqid = self.data('DBXRef')[gene]['EnsLoci']
                desc = self.data('DBXRef')[gene]['EnsDesc']
                ensg = self.data('DBXRef')[gene]['EnsEMBL']
                try: godata = self.dict['GO'][ensg]
                except: pass
            html = htmlHead(gene,self.list['StyleSheets']) + self.seqDetailsHTML(gene,dbxref) #gene,seqid,dbxref,desc,godata)
            tablist = []
            ## ~ [1a] ~ Interactors tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('Interactors',self.interactorTableHTML(gene),'Table of returned %s interactors' % gene))
            tablist.append(('Full-PPI',self.interactorTableHTML(gene,full=True),'Full list of %s PPI' % gene))
            ## ~ [1b] ~ Raw Data tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rawdb = self.db().getTable('Raw')
            rawdata = rawdb.subset('NRID',gene)
            if rawdata: tablist.append(('Raw',self.resultTableHTML(rawdb.fields(),rawdata,drop=['NRID','HGNC','Description']),'Raw results data for %s' % gene))
            else: tablist.append(('Raw','<i>No proteomics identifications mapped to %s</i>' % gene,'Raw results data for %s' % gene))
            ## ~ [1c] ~ Interactome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = self.db().getTable('PPI')
            hdb = pdb.subset('Hub',gene)
            main = self.db().getTable('Main')
            ilist = rje.sortUnique(pdb.dataList(hdb.values(),'Spoke'))
            for i in ilist[0:]:
                if i not in main.index('NRID'): ilist.remove(i)
            if len(ilist) > 2 or (gene not in ilist and len(ilist) > 1):
                intpng = '../interactome_png/%s.interactome.png' % gene
                alt = '%s interactome image not yet made' % gene
                title = 'Graphical representation of %s interactome (returned proteins only). Click to open in new window. Larger interactomes may not be clear. See tables for details.' % gene
                tablist.append(('Interactome','<h2>%s Interactome</h2>\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (gene,intpng,intpng,alt,title),title))
            ## ~ [1d] ~ Specificity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # - Nested interactome: ITGA2B, AA, FF (inside), spokes in genelist outside
            if ilist: 
                intpng = '../nested_png/%s.nested.png' % gene
                alt = '%s nested interactome image not yet made' % gene
                title = 'Graphical representation of returned %s interactors and their interactions with AA, FF and ITGA2B/ITGAL/ITGA5. Interactions between spokes are not shown - see "Interactome". Click to open in new window. Larger interactomes may not be clear. See tables for details.' % gene
                tablist.append(('Specificity','<h2>%s Interactome</h2>\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (gene,intpng,intpng,alt,title),title))
            ## ~ [1e] ~ Complex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            intpng = '../complex_png/%s.complex.png' % gene
            alt = '%s complex interactome image not yet made' % gene
            title = 'Graphical representation of possible proteins complexed to %s. Click to open in new window. Larger interactomes may not be clear.' % gene
            tablist.append(('Complex','<h2>%s Complex</h2>\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (gene,intpng,intpng,alt,title),title))
            ### ~ [2] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
            html += tabberHTML(gene,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('genePage Error')
#########################################################################################################################
    def goPage(self,go,goid):   ### Make GO Pages
        '''Make GO Pages.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dbxref = self.db().getTable('DBXRef')
            godb = self.db().getTable('GO')
            domdb = self.db().getTable('PFam')
            odb = self.db().getTable('Occ')
            genes = []      # Genes
            sighubs = []    # Genes + domains
            sigspokes = []  # EnsLoci
            domains = []    # Domains
            slims = []      # Patterns
            ## ~ [0a] ~ Populate lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ensg in godb.dataList(godb.indexEntries('GO_Desc',go),'EnsG'):
                for gene in dbxref.dataList(dbxref.indexEntries('EnsEMBL',ensg),'Gene'):
                    if gene not in genes: genes.append(gene)
                    if gene not in sighubs and gene in odb.index('Hub'): sighubs.append(gene)
                for prot in dbxref.dataList(dbxref.indexEntries('EnsEMBL',ensg),'EnsLoci'):
                    if prot not in sigspokes and prot in odb.index('Seq'): sigspokes.append(prot)
                    for domain in domdb.dataList(domdb.indexEntries('Name',prot),'Type'):
                        if domain not in domains: domains.append(domain)
                        if domain not in sighubs and domain in odb.index('Hub'): sighubs.append(domain)
            for hub in sighubs:
                for pattern in odb.dataList(odb.indexEntries('Hub',hub),'Pattern'):
                    if pattern not in slims: slims.append(pattern)
            for spoke in sigspokes:
                for pattern in odb.dataList(odb.indexEntries('Seq',spoke),'Pattern'):
                    if pattern not in slims: slims.append(pattern)
            ### ~ [1] ~ Make Page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html = [htmlHead(go,self.list['StyleSheets']),
                    '','<!-- ~~~~ %s GO Summary details table ~~~~ -->' % goid,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="40%%"><a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s<b>GO:%s</b></FONT></a></td>' % (goid,gfont,goid),
                    '<td width="40%%">%s%s</FONT></td>' % (ifont,godb.dataList(godb.indexEntries('GO_Desc',go),'GO_Type')[0]),
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,go)]
            html += ['</table>','<!-- ~~~~ End %s GO Summary details table ~~~~ -->' % goid,'']
            ### ~ [2] Domain info Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = []; jtxt = ' ~ \n'
            ## ~ [2a] ~ Gene List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            genes.sort()
            for gene in genes: tabhtml.append(geneLink(gene))
            if tabhtml: tabhtml[0] = '<h2>%s Genes</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No genes in analysis mapped to %s</i>' % go]
            tablist.append(('Genes',string.join(tabhtml,jtxt),'Genes mapped to %s' % go))
            ## ~ [2b] ~ Domain List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            domains.sort()
            for hub in domains: tabhtml.append(domainLink(hub))
            if tabhtml: tabhtml[0] = '<h2>%s Domains</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No domains in analysis mapped to %s</i>' % go]
            tablist.append(('Domains',string.join(tabhtml,jtxt),'Domains mapped to %s' % go))
            ## ~ [2c] ~ SigHub List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            sighubs.sort()
            for hub in sighubs:
                if hub in genes: tabhtml.append(geneLink(hub))
                else: tabhtml.append(domainLink(hub))
            if tabhtml: tabhtml[0] = '<h2>%s Hubs with significant results</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s hubs returned significant results</i>' % go]
            tablist.append(('Hubs',string.join(tabhtml,jtxt),'Significant Hubs mapped to %s' % go))
            ## ~ [2d] ~ SigSpoke List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            sigspokes.sort()
            for prot in sigspokes:
                gene = self.gene(prot)
                if gene != prot: tabhtml.append(geneLink(gene))
                else: tabhtml.append(prot)
            if tabhtml: tabhtml[0] = '<h2>%s Spokes with significant results</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s spokes returned significant results</i>' % go]
            tablist.append(('Spokes',string.join(tabhtml,jtxt),'Significant Spokes mapped to %s' % go))
            ## ~ [2e] ~ SLiM List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            slims.sort()
            for slim in slims: tabhtml.append(slimLink(slim))
            if tabhtml: tabhtml[0] = '<h2>%s SLiM predictions</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s significant SLiM predictions</i>' % go]
            tablist.append(('SLiMs',string.join(tabhtml,jtxt),'Significant SLiM predictions mapped to %s' % go))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
            html = string.join(html,'\n') + tabberHTML(goid,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('goPage Error')
#########################################################################################################################
    ### <4> ### Tables pages                                                                                            #
#########################################################################################################################
    def interactorTableHTML(self,gene,full=False): ### Table of interactors with evidence & SLiMs returned for 4 datasets types
        '''
        Table of interactors in sample: Spoke, Evidence, Description, Class, AA, FF      
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = []
            ppi = self.db().getTable('PPI')
            hdb = ppi.subset('Hub',gene)
            main = self.db().getTable('Main')
            ilist = rje.sortUnique(ppi.dataList(hdb.values(),'Spoke'))
            ### ~ [2] Generate table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Table Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = ['','<h2>Protein-Protein Interactions</h2>','',
                    '<table width=1000 border=%d><tr><th colspan=3>%s Interactions</th><th colspan=5>Peptide Binding Results</th></tr>' % (self.stat['Border'],gene)]
            html += ['<tr><th align=left title="%s">Gene</th>' % self.titleText('ppi','gene'),
                     '<th align=left title="%s">Evidence</th>' % self.titleText('ppi','evidence'),
                     '<th align=left title="%s">Description</th>' % self.titleText('ppi','description'),
                     '<th align=left title="%s">Class</th>' % self.titleText('results','Class'),
                     '<th align=left title="%s">AA</th>' % self.titleText('results','AA'),
                     '<th align=left title="%s">enrAA</th>' % self.titleText('results','enrAA'),
                     '<th align=left title="%s">FF</th>' % self.titleText('results','FF'),
                     '<th align=left title="%s">enrFF</th>' % self.titleText('results','enrFF')]
            html.append('</tr>')
            ## ~ [2b] Table contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ix = 0
            for i in ilist[0:]:
                if i not in main.index('NRID') and not full: continue
                if i in ['AA','FF']: continue
                html += ['<tr valign="top"><td width=90>%s</td>' % self.geneLink(i)]
                hkey = '%s\t%s' % (gene,i)
                html += ['<td width=330>%s</td>' % string.join(string.split(hdb[hkey]['Evidence'],'|'),';<br>')]
                html += ['<td width=330>%s</td>' % self.geneDesc(i)]
                if i in main.index('NRID'):
                    html += ['<td width=50>%s</td>' % main.data()[i]['Class']]
                    html += ['<td width=50>%s</td>' % rje.expectString(main.data()[i]['AA'])]
                    html += ['<td width=50>%s</td>' % rje.expectString(main.data()[i]['enrAA'])]
                    html += ['<td width=50>%s</td>' % rje.expectString(main.data()[i]['FF'])]
                    html += ['<td width=50>%s</td>' % rje.expectString(main.data()[i]['enrFF'])]
                else: html += ['<td width=50>-</td>'] * 5
                html += ['</tr>','']
                ix += 1
            if not ix: html += ['<tr><td align=left colspan=6><i>No %s interactors returned by proteomics or known ITGA2B/ITGAl/ITGA5 interactors</i></td></tr>' % gene]
        except: self.errorLog('interactorTableHTML Error')
        ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += ['</table>','']
        return string.join(html,'\n')
#########################################################################################################################
    def resultTableHTML(self,headers,data,asdict=True,drop=[],width={}):     ### Write data to HTML table
        '''Write data to HTML table.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = headers[0:]   # Do not want to change original headers list!
            for h in drop:
                if h in headers: headers.remove(h)
            if asdict:
                datalist = []
                for key in rje.sortKeys(data): datalist.append(data[key])
            else: datalist = data
            ### ~ [2] ~ Write Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hcode = ''
            for h in headers:
                if h in width: hcode += '<th width=%s title="%s">%s</th>\n' % (width[h],self.titleText('results',h),h)
                else: hcode += '<th title="%s">%s</th>\n' % (self.titleText('results',h),h)
            rows = [hcode]
            #rows = ['<TD><B>%s</B></TD>' % string.join(headers,'</B></TD><TD><B>')]
            ## ~ [2b] ~ Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for dat in datalist:
                rdata = []
                for h in headers:
                    if h in dat:
                        if h[-2:] in ['AA','FF'] or h.find('.v.') > 1 or h[:3] in ['Rep','Int']:
                            try:
                                try: rdata.append(rje.expectString(dat[h]))
                                except: rdata.append(rje.expectString(string.atof(dat[h])))
                            except: rdata.append('-')
                        elif h in ['NRID','Hub','HGNC','Spoke']:
                            gene = self.gene(str(dat[h]))
                            rdata.append(self.geneLink(gene,frontpage=False,altgene=dat[h]))
                        else: rdata.append(str(dat[h]))
                    else: rdata.append('')
                    if h in width: rdata[-1] = '<TD WIDTH=%s>%s</TD>' % (width[h],rdata[-1])
                    else: rdata[-1] = '<TD>%s</TD>' % (rdata[-1])
                rows.append(string.join(rdata,''))
        except: self.errorLog('resultTableHTML error')
        return '<TABLE BORDER=%d><TR VALIGN="top">\n%s\n</TR></TABLE>\n' % (self.stat['Border'],string.join(rows,'\n</TR><TR VALIGN="top">\n'))
#########################################################################################################################
    ### <5> ### Misc HTML pages                                                                                         #
#########################################################################################################################
#########################################################################################################################
    ### <6> ### PNG Methods                                                                                             #
#########################################################################################################################
    def interactomePNG(self,gene):   ### Peforms the SLiMJim visualisation for a given Hub dataset               #V2.0
        '''Peforms the SLiMJim visualisation for a given Hub dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pngdir = '%sinteractome_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            basefile = rje.makePath('%s%s.interactome' % (pngdir,gene),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            pdb = self.db().getTable('PPI')
            hdb = pdb.subset('Hub',gene)
            main = self.db().getTable('Main')
            ilist = rje.sortUnique(pdb.dataList(hdb.values(),'Spoke'))
            for i in ilist[0:]:
                if i not in main.index('NRID'): ilist.remove(i)
            if len(ilist) < 2 or (gene in ilist and len(ilist) < 3): return    # No need!
            dis = self.obj['Dis']
            ### ~ [2] ~ Load distance matrix and make Tree/Heatmap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            nsftree = dis.upgma(objkeys=ilist,checksym=False)
            open('%s.nsf' % basefile,'w').write(nsftree)
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            upgma.rTree('%s.tree.csv' % basefile,seqname='short')
            reorder = upgma._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            dis.info['Name'] = '%s interactome GABLAM' % gene
            dis.saveMatrix(reorder,basefile+'.heatmap.tdt',delimit='\t')
            ### ~ [3] ~ Add PPI links between spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# New code #!#
            pobj = rje_ppi.PPI(self.log,self.cmd_list)
            ppi = ()
            for hub in ilist:
                ppi[hub] = []
                for spoke in pdb.indexDataList('Hub',hub,'Spoke'): ppi[hub].append(spoke)
            npos = pobj.rjeSpringLayout(ppi)
            classcol = {}   #{'AA':15,'FF':19,'AA-FF':8,'LYS':6,'ITGA2B':3,'ITGAL':3,'ITGA5':3,'Integrin':3}
            pobj.addCol(default=12,coldict=classcol,G=ppi,ckey='Class',edb=main)
            pobj.saveR(npos,basefile,ppi,cleantdt=False,backups=False)

            #return
            #names = reorder[0:]
            #if gene in names: names.remove(gene)
            #pfile = '%s.ppi.tdt' % basefile
            #if os.path.exists(pfile): os.unlink(pfile)
            #pdat = {}; sorted = {}
            #for p1 in names:
            #    datadict = {}
            #    psort = []
            #    for p2 in names:
            #        try: ppi = p2 in pdb.dataList(pdb.indexEntries('Hub',p1),'Spoke')    #self.dict['PPI'][p1]  #False
            #        except: ppi = False
            #        if ppi: datadict[p2] = 1
            #        else: datadict[p2] = 0
            #        psort.append('%s' % (1 - datadict[p2]))
            #    pdat[p1] = datadict
            #    psort = '%s-%s' % (rje.preZero(10000-sum(datadict.values()),10000),string.join(psort))
            #    if psort not in sorted: sorted[psort] = [p1]
            #    else: sorted[psort] += [p1]
            #porder = []
            #for psort in rje.sortKeys(sorted): porder += sorted[psort]
            #rje.delimitedFileOutput(self,pfile,porder)
            #for p1 in porder: rje.delimitedFileOutput(self,pfile,porder,datadict=pdat[p1])
            #datadict = {}
            #for p1 in porder: datadict[p1] = self.classCol(main.data()[p1]['Class'])
            #rje.delimitedFileOutput(self,pfile,porder,datadict=datadict)

            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "hminteractome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/hmhtml.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#self.deBug('%s.png' % basefile) 
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def nestedPNG(self,gene):     ### Performs SLiMJIM visualisation for a given motif                        #V2.0
        '''Performs SLiMJIM visualisation for a given motif.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pngdir = '%snested_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            basefile = '%s%s.nested' % (pngdir,gene)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            pdb = self.db().getTable('PPI')
            hdb = pdb.subset('Hub',gene)
            main = self.db().getTable('Main')
            ilist = rje.sortUnique(pdb.dataList(hdb.values(),'Spoke'))
            for i in ilist[0:]:
                if i not in main.index('NRID'): ilist.remove(i)
            if not ilist: return    # No need to make!
            ### ~ [2] ~ Generate nested network PPI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Inside = Hubs = Headers (& 'Spoke') / Outside = Spokes = rows
            if gene in ['ITGA2B','ITGAL','ITGA5']:
                headers = ['Spoke','ITGA2B','ITGAL','ITGA5']
                headers.remove(gene); headers.insert(1,gene); headers.insert(3,'AA'); headers.insert(4,'FF')
                hubs = headers[2:]
            else: headers = ['Spoke',gene,'ITGA2B','AA','ITGAL','FF','ITGA5']; hubs = headers[2:]
            spokes = ilist[0:]
            headers.append('Col')
            nfile = '%s.tdt' % basefile
            if os.path.exists(nfile): os.unlink(nfile)
            rje.delimitedFileOutput(self,nfile,headers)
            ndata = {}; sortdict = {'None':[]}
            for spoke in spokes:
                datadict = {'Spoke':spoke,'Col':self.classCol(main.data()[spoke]['Class']),gene:0}
                for hub in hubs:
                    h = hub     # humsf throwback
                    datadict[h] = 0
                    if spoke in pdb.dataList(pdb.indexEntries('Hub',hub),'Spoke'): datadict[h] = 1
                    if datadict[h] and 'Sort' not in datadict: datadict['Sort'] = h
                    elif datadict[h] and headers.index(h) < headers.index(datadict['Sort']): datadict['Sort'] = h
                try: skey = datadict.pop('Sort')
                except: skey = gene
                if skey not in sortdict: sortdict[skey] = []
                sortdict[skey].append(spoke)
                ndata[spoke] = datadict
            spokelist = []; altx = 0; #nonex = len(sortdict['None'])
            for h in headers[1:-1]:
                if h not in sortdict: continue
                for spoke in sortdict[h]: spokelist.append(spoke)
                if h == headers[1] and len(hubs) > 1: altx = len(sortdict[h]) / 2
                #sortdict['None'] = sortdict['None'][(nonex/3):]
            #if sortdict['None']: spokelist += sortdict['None']
            if altx: spokelist = spokelist[altx:] + spokelist[:altx]
            for spoke in spokelist: rje.delimitedFileOutput(self,nfile,headers,datadict=ndata[spoke])
            ### ~ [3] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "hmnested" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/hmhtml.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.deBug('%s.png' % basefile) 
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['tdt','r.tmp.txt']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def complexPNG(self,gene,maxsize=50):  ### Generates PNG file for predicted protein complex membership
        '''Generates PNG file for predicted protein complex membership.'''
        try:#### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pngdir = '%scomplex_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            basefile = '%s%s.complex' % (pngdir,gene)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            pdb = self.db().getTable('PPI')
            main = self.db().getTable('Main')
            nfile = '%s.tdt' % basefile
            ### ~ [1] ~ Establish complex using PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            core = [gene]   # List of core genes = query + genes with 2+ "Complex" ppi
            fuzz = []       # Peripheral proteins linked by one "complex" PPI only
            csize = 1       # Size of complex
            cppi = {gene:[]}       # PPI dictionary for complex only {hub:[spokes]}
            growth = [gene]
            while csize < maxsize and growth:
                oldsize = csize     # Store current size to check for growth?
                oldgrowth = growth[0:]; growth = []
                ## ~ [1a] ~ Add "Complex" PPI with growth genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for hub in oldgrowth:
                    for ppi in pdb.indexEntries('Hub',hub):
                        if ppi['Evidence'].lower().find('complex') < 1: continue
                        spoke = ppi['Spoke']
                        if spoke in core + fuzz: continue   # Already added
                        cppi[hub].append(spoke)
                        cppi[spoke] = []
                        for link in pdb.indexEntries('Hub',spoke):
                            if link['Evidence'].lower().find('complex') < 1: continue
                            if link['Spoke'] not in core + fuzz: continue
                            cppi[spoke].append(link['Spoke'])
                        if len(cppi[spoke]) == 0: raise ValueError
                        elif len(cppi[spoke]) == 1: fuzz.append(spoke)
                        else: core.append(spoke); growth.append(spoke)
                ## ~ [1b] ~ Check fuzz status and reduce growth to core genes ~~~~~~~~~~~~~~~~~~~~~ ##
                for spoke in fuzz:
                    if len(cppi[spoke]) == 0: raise ValueError
                    elif len(cppi[spoke]) == 1: continue
                    fuzz.remove(spoke); core.append(spoke); growth.append(spoke)
                csize = len(core) + len(fuzz)
            if csize > maxsize: fuzz = []

            #!# New code #!#
            pobj = rje_ppi.PPI(self.log,self.cmd_list)
            for hub in cppi:
                for spoke in cppi[hub]:
                    if hub not in cppi[spoke]: cppi[spoke].append(hub)
            for hub in cppi: cppi[hub].sort()
            npos = pobj.rjeSpringLayout(cppi)
            classcol = {}   #{'AA':15,'FF':19,'AA-FF':8,'LYS':6,'ITGA2B':3,'ITGAL':3,'ITGA5':3,'Integrin':3}
            pobj.addCol(default=12,coldict=classcol,G=cppi,ckey='Class',edb=main)
            pobj.saveR(npos,basefile,cppi,cleantdt=True,backups=False)
            xgmml = pobj.ppiToXGMML(cppi,'%s.complex' % gene)
            xgmml.dict['NodePos'] = npos
            xgmml.saveXGMML('%s.xgmml' % basefile)
            return       
            
            ### ~ [2] ~ Generate output with Core in centre and fuzz outside ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Make pure interactome. In future add two levels.
            core.sort(); fuzz.sort()
            names = core + fuzz; names.remove(gene); names.insert(0,gene)
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for p1 in names:
                datadict = {}
                for p2 in names:
                    if p2 in cppi[p1]: datadict[p2] = 1
                    else: datadict[p2] = 0
                rje.delimitedFileOutput(self,pfile,names,datadict=datadict)
            datadict = {}
            for p1 in names:
                try: datadict[p1] = self.classCol(main.data()[p1]['Class'])
                except: datadict[p1] = 12
            rje.delimitedFileOutput(self,pfile,names,datadict=datadict)

            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "interactome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/hmhtml.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#self.deBug('%s.png' % basefile) 
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['ppi.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: SLiMHTML Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
class SLiMFinder(slimfinder.SLiMFinder):
    '''For pickling only'''
def htmlHead(title,stylesheets=['../example.css','../redwards.css'],tabber=True,frontpage=False):    ### Returns text for top of HTML file
    '''
    Returns text for top of HTML file.
    >> title:str = Title of webpage.
    >> stylesheets:list = List of stylesheets to use.
    >> tabber:bool = whether page has tabber tabs
    '''
    html = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
            '<html lang="en">','','<!-- ~~~~~~~~~~~~~~~ HTML head data ~~~~~~~~~~~~~~~~~ -->','<head>',
            '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">', #!# Add additional metadata #!#
            '<title>%s</title>' % title,'']
    for stylesheet in stylesheets:
        if frontpage: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % string.replace(stylesheet,'../','./'))
        else: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % stylesheet)
    if tabber:
        html += ['','<!-- ~~~~~~~~~~~ Tabber Javascript ~~~~~~~~~~~ -->','<script type="text/javascript">',
                 'document.write(\'<style type="text/css">.tabber{display:none;}<\/style>\');',
                 'var tabberOptions = {',' /* Optional: start manually (run tabber) at end of file','*/',
                 '\'manualStartup\':true','};','</script>','<!-- Load the tabber code -->']
        if frontpage: html += ['<script type="text/javascript" src="./javascript/tabber.js"></script>','']
        else: html += ['<script type="text/javascript" src="../javascript/tabber.js"></script>','']
    html += ['</head>','<!-- ~~~~~~~~~~~~~~~ End of HTML head data ~~~~~~~~~~~~~~~~~ -->','','<body>','']
    #print string.join(html,'\n')
    return string.join(html,'\n')
#########################################################################################################################
def htmlTail(copyright='RJ Edwards 2010',tabber=True):  ### Returns text for bottom of HTML
    '''
    Returns text for bottom of HTML.
    >> copyright:str = copyright text'
    >> tabber:bool = whether page has tabber tabs
    '''
    t = string.split(time.asctime(time.localtime(time.time())))
    datetime = '%s %s %s' % (t[2],t[1],t[-1])
    html = ['','<!-- ~~~~~~~~~~~~~~ HTML tail data ~~~~~~~~~~~~~~~~~ -->',
            '<HR><FONT COLOR=#979E45 SIZE=2>&copy; %s. Last modified %s.</FONT></P>' % (copyright,datetime),'',
            '<script type="text/javascript">','/* manualStartup=true so need to run it now */',
            'tabberAutomatic(tabberOptions);','</script>','','</body>','</html>',
            '<!-- ~~~~~~~~~~~~~~ End of HTML tail data ~~~~~~~~~~~~~~~~~ -->']
    #print string.join(html,'\n')
    return string.join(html,'\n')
#########################################################################################################################
def tabberHTML(id,tablist,level=0):     ### Returns text for Tabber HTML
    '''
    Returns text for Tabber HTML.
    >> id:str = Identifier for Tabber object
    >> tablist:list = List of (tab_title, tab_html_text) tuples
    >> level:int = Level of Tabber object (base = level)
    '''
    jointxt = '\n' + '    ' * level 
    html = ['<!-- ~~~~~~~~~~~~~~~ %s Tabber Div ~~~~~~~~~~~~~~~ -->' % id,'','<div class="tabber" id="%s">' % id,'']
    #print html
    for tab in tablist:
        #print tab
        #print tab[0],tab[1]
        #print tabberTabHTML(tab[0],tab[1])
        if len(tab) > 2: html += string.split(tabberTabHTML(tab[0],tab[1],tab[2]),'\n')
        else: html += string.split(tabberTabHTML(tab[0],tab[1]),'\n')
    html += ['</div>','<!-- ~~~~~~~~~~~~~~~ End of %s Tabber Div ~~~~~~~~~~~~~~~ -->' % id,]
    #print string.join(html,jointxt)
    return string.join(html,jointxt)
#########################################################################################################################
def tabberTabHTML(id,text,title=''):          ### Returns text for TabberTab HTML
    '''
    Returns text for TabberTab HTML.
    >> title:str = Text for title of TabberTab
    >> text:str = HTML text for TabberTab content
    '''
    if not title: title = id
    html = ['','<!-- ~~~ %s TabberTab div ~~~ -->' % id,'<div class="tabbertab" title="%s" id="%s">' % (title,id),'']
    html += string.split(text,'\n')
    html += ['','</div>','<!-- ~~~ %s TabberTab end ~~~ -->' % id]
    #print string.join(html,'\n  ')
    if string.join(html).upper().find('<PRE>') >= 0: return string.join(html,'\n')       
    else: return string.join(html,'\n  ')
#########################################################################################################################
def geneLink(gene,frontpage=False):     ### Returns gene link text
    '''Returns gene link text.'''
    if frontpage: return '<a href="./gene/%s.html" target="_blank" title="%s results page">%s</a>' % (gene,gene,gene)
    else: return '<a href="../gene/%s.html" title="%s results page">%s</a>' % (gene,gene,gene)
#########################################################################################################################
def domainLink(domain,frontpage=False):     ### Returns gene link text
    '''Returns domain link text.'''
    if frontpage: return '<a href="./domain/%s.html" target="_blank" title="%s domain results page">%s</a>' % (domain,domain,domain)
    else: return '<a href="../domain/%s.html" title="%s domain results page">%s</a>' % (domain,domain,domain)
#########################################################################################################################
def randLink(dset,frontpage=False):     ### Returns gene link text
    '''Returns random dataset link text.'''
    if frontpage: return '<a href="./%s/%s.html" target="_blank" title="Random dataset results page">%s</a>' % (dset[:4],dset,dset)
    else: return '<a href="../%s/%s.html" title="Random dataset results page">%s</a>' % (dset[:4],dset,dset)
#########################################################################################################################
def slimLink(pattern,frontpage=False):  ### Returns SLiM link text
    '''Returns gene link text.'''
    if frontpage: return '<a href="./slim/%s.html" target="_blank" title="%s results page">%s</a>' % (rje_slim.slimFromPattern(pattern),pattern,pattern)
    else: return '<a href="../slim/%s.html" title="%s results page">%s</a>' % (rje_slim.slimFromPattern(pattern),pattern,pattern)
#########################################################################################################################
def goLink(go,goid,frontpage=False):  ### Returns GO link text
    '''Returns GO link text.'''
    return '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s" title="Browse AmiGO">%s</a>' % (goid,go)
    if frontpage: return '<a href="./go/%s.html" target="_blank" title="GO results page">%s</a>' % (goid,go)
    else: return '<a href="../go/%s.html" title="GO results page">%s</a>' % (goid,go)
#########################################################################################################################
### END OF SECTION III: MODULE METHODS                                                                                  #
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
    try: HMHTML(mainlog,cmd_list).run()

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
