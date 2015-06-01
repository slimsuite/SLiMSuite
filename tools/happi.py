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
Program:      HAPPI
Description:  Hyperlinked Analysis of Protein-Protein Interactions
Version:      1.2
Last Edit:    06/03/13
Citation:     Edwards et al. (2011), Molecular Biosystems DOI: 10.1039/C1MB05212H. [PMID: 21879107]
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module makes a set of linked webpages for the quick and dirty analysis of PPI networks around one or more
    sets of genes. A table of genes and their corresponding classes is used to make the initial front page tables.

    Additional data tables can also be loaded and linked via the "Gene" field for individual Gene Pages. These must
    be named in the format X.N.tdt, where N will be used as the tab name.

    Additional classes can be added using the makeclass and classorder lists. Any classes in classorder that are not
    found in indata will be added using the rules set by makeclass:
    - go = GO term description or ID
    - keyword = found in GO term description or protein description
    - desc = found in protein description
    - xref = adds PPI for protein identified from xrefdata table
    - add = will only add a gene to a class if it does not already have one (no multiple-class genes)

Commandline:
    ### ~ INPUT ~ ###
    indata=FILE     : Input data for front pages []
    pagehead=X      : Text for Front Page Header ['HAPPI Analysis of basefile']
    infotext=X      : Text to display under header ['This data has not yet been published and should not be used without permission.']
    genefield=X     : Field to be used for indentifying Genes ['Gene']
    geneclass=X     : Field used to identify class of Genes ['Class']
    classcol=X      : Table of "Class" and "Col" (soton$col indexes) to be used for PPI images [basefile.col.tdt]
    fillcol=T/F     : Fill in colour for missing class combinations [True]
    classorder=LIST : List of Class orders (otherwise alphabetical) []
    multiclass=T/F  : Whether to allow membership of multiple classes (joined by "-" [True]
    makeclass=LIST  : Generate classes from classorder LIST if missing from indata (go/keyword/desc/xref(+ppi)/add*) [keyword]
    genedata=LIST   : List of additional data tables for Gene-centric pages. Must have "Gene" or genefield field. []
    xrefdata=FILE   : File of Database cross-references. Must have "Gene" field. []
    xreftab=T/F     : Add database cross-reference tabs to the front page Class tabs [True]
    pairwise=FILE   : Pairwise protein-protein interation (PPI) file []
    addppi=LIST     : List of additional PPI pairwise files to add []
    addclass=T/F    : Whether to add the Classes themselves to the PPI networks as nodes [False]
    godata=FILE     : Delimited file of GO Data (Gene,GO_ID,GO_Type,GO_Desc) [basefile.go.tdt]
    gablam=FILE     : Delimited GABLAM results file for homology data [basefile.gablam.tdt]
    ppexpand=X      : Expand PPI by X levels for MCODE complex generation [1]
    ppcomplex=LIST  : List of different evidence codes for special MCODE analyses ['Complex']
    ppextra=T/F     : Make additional pages for genes added and returned in MCODE clusters [False]
    combine=LIST    : List of Classes to combine into Interactome for front page ('all' or '*' for all) []
    special=LIST    : Execute special analysis code for specific applications []

    ### ~ BASIC OUTPUT OPTIONS ~ ###
    htmlpath=PATH       : Path of parent html directory [./html]
    stylesheets=LIST    : List of CSS files to use ['../example.css','../redwards.css']
    border=X            : Border setting for tables [0]
    dropfields=LIST     : Fields to exclude from summary tables []
    makepng=T/F         : Whether to (look for and) make PNG files with R [True]
    pngmax=X            : Max number of genes for PNG construction [2000]
    svg=T/F             : Use SVG files instead of PNG files [True]
    xgmml=T/F           : Whether to also output XGMML files in PNG/SVG directories [False]
    makepages=LIST      : Types of pages to make [front,gene,go,mcode,class,interactome,expand]
    titletext=FILE      : File containing (Page,ID,Title) for mouseover text [titletext.tdt]
    gopages=X           : Create Pages of GO terms with X+ representative genes [5]
    maxgo=X             : Go terms above X genes will not have gene data tabs [500]
    nobots=T/F          : Whether to insert no-bot meta tag to pages [True]

    ### ~ FIGURE REMAKING OPTIONS ~ ###
    usepos=T/F          : Whether to use existing Node positions if found [True]
    updatepos=T/F       : Whether to run an additional round of the layout algorithm if npos found [False]

    *** See also RJE_PPI Fragment and MCODE options ***
    *** See also RJE_GO options ***

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_ppi, rje_uniprot, rje_zen
Other modules needed: None
"""
#    ### ~ ADVANCED OUTPUT OPTIONS ~ ###
#    fakehtml=T/F        : Whether to make UniFake HTML [True]
#    unifake=PATH        : Path to UniFake dat file(s) [./unifake/]
#    unipath=PATH        : Path to real UniProt dat file(s) [./uniprot/]
#    unireal=LIST        : Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]

#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_go, rje_html, rje_ppi, rje_sequence, rje_slim, rje_uniprot, rje_zen
import rje_dismatrix_V2, rje_tree
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation. Based on rje_hm_html.
    # 0.1 - Added special methods for specific analyses.
    # 0.2 - Added multiclass capability.
    # 0.3 - Added GO enrichment to MCODE complexes.
    # 0.4 - Added "expanded" MCODE complex and updated GO enrichment of complexes.
    # 0.5 - Added class colour compilation.
    # 0.6 - Replaced ppcomplex=T/F with list of different evidence codes for special MCODE analyses ['Complex']
    # 0.7 - Added ppextra option for adding additional pages for genes added and returned in MCODE clusters.
    # 0.8 - Added use of SVG PPI images in place of PNG.
    # 0.9 - Added GZipping of XGMML files to save space.
    # 1.0 - Fixed PNG/SVG XGMML directory issue. First release.
    # 1.1 - Fixed SVG class naming issue (Cannot have *.class.svg!)
    # 1.2 - Added addclass and refined output for Host-Pathogen PPI analysis.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Make and test code.
    # [ ] : Add advanced UniFake Options.
    # [Y] : Consider adding MCODE from "complex" evidence only
    # [Y] : Add optional loading of data types for additional data tables.
    # [Y] : Add Unique ranking of MCODE complexes. (RJE_DB.PY and RJE.PY)
    # [Y] : Add ordering of classes (rather than alphabetical)
    # [X] : Add option to subdivide classes into alphabetical tabs? # Added Gene List tab instead
    # [Y] : Add class fluff (singletons) to MCODE complexes?
    # [Y] : Read in MCODE tables if pre-existing.
    # [Y] : Add xgmml=T option for XGMML output into PNG directories (allows easy remaking of key images).
    # [X] : Handle genes in multiple classes. (string.split | ?)    # No need - indata removes
    # [Y] : Add XRefTab option for a database cross-reference tab on the front page.
    # [Y] : Add colourkey to front page.
    # [ ] : Add optional reading in of formats for data fields.
    # [ ] : Add "sideways" results table option for single dictionary?
    # [Y] : Add saving and reading of node positions to allow quick remaking of different size/colours.
    # [ ] : Add proper use of GOPHER and/or distance matrix for homology information.
    # [ ] : Add link to http://www.eurexpress.org/ee/databases/tdbsearch?mode=search&standardAtlas=on&getStandardAtlas=yes&multistageAtlas=on&getMultistageAtlas=yes&groupSelection=Group_All&findEverything=no&anatomySetNames=0&expressionStrengthOperator=g&expressionStrengthValue=&coverageOperator=g&coverageValue=&page=&searchField=All&searchValue=PAX6&Submit=Go
    # [ ] : Add use of CGI to reduce file volumes. (Maybe PHP with time?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('HAPPI', '1.1', 'July 2011', '2010')
    description = 'Hyperlinked Analysis of Protein-Protein Interactions'
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
            if rje.yesNo('Show RJE_PPI options?'): out.verbose(-1,4,text=rje_ppi.__doc__)
            if rje.yesNo('Show RJE_GO options?'): out.verbose(-1,4,text=rje_go.__doc__)
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
### SECTION II: HAPPI Class                                                                                             #
#########################################################################################################################
class HAPPI(rje.RJE_Object):     
    '''
    HAPPI Class. Author: Rich Edwards (2010).

    Info:str
    - InData = Input data for front pages []
    - GABLAM = Delimited GABLAM results file for homology data []
    - GeneField = Field to be used for indentifying Genes ['Gene']
    - GeneClass = Field used to identify class of Genes ['Class']
    - GOData = Delimited file of GO Data (Gene,GO_ID,GO_Type,GO_Desc) [basefile.go.tdt]
    - InfoText = Text to display under header ['This data has not yet been published and should not be used without permission.']
    - ClassCol = Table of "Class" and "Col" (soton$col indexes) to be used for PPI images []
    - XRefData = File of Database cross-references. Must have "Gene" field.
    - Pairwise = Pairwise protein-protein interation (PPI) file []
    - PageHead = Text for Front Page Header ['HAPPI Analysis of basefile']
    - HTMLPath = Path of parent html directory [./html]
    - TitleText = File containing (Page,ID,Title) [titletext.tdt]
    - UniFake = Path to UniFake dat file(s) [./unifake/]
    - UniPath = Path to real UniProt dat file(s) [./uniprot/]
    
    Opt:boolean
    - AddClass = Whether to add the Classes themselves to the PPI networks as nodes [False]
    - AlphTabs = Subdivide genes within each class alphabetically on front page [False]
    - FakeHTML = Whether to make UniFake HTML [True]
    - FillCol = Fill in colour for missing class combinations [True]
    - MakePNG = Whether to (look for and) make PNG files with R [True]
    - MultiClass = Whether to allow membership of multiple classes (joined by "-" [True]
    - NoBots = Whether to insert no-bot meta tag to pages [True]
    - PPExtra = Make additional pages for genes added and returned in MCODE clusters [False]
    - SVG = Use SVG files instead of PNG files [True]
    - UsePos = Whether to use existing Node positions if found [True]
    - UpdatePos = Whether to run an additional round of the layout algorithm if npos found [False]
    - XGMML = Whether to output XGMML into PNG directories [False]
    - XRefTab = Add database cross-reference tabs to the front page Class tabs [True]

    Stat:numeric
    - Border = Border setting for tables [0]
    - GOPages = Create Pages of GO terms with X+ representative genes [2]
    - MaxGO = Go terms above X genes will not have gene data tabs [500]
    - PNGMax = Max number of genes for PNG construction [2000]
    - PPExpand = Expand PPI by X levels for MCODE complex generation [1]

    List:list
    - AddPPI = List of additional PPI pairwise files to add []
    - Combine = List of Classes to combine into Interactome for front page []
    - ClassOrder = List of Class orders (otherwise alphabetical) []
    - DropFields = Fields to exclude from summary tables []
    - GeneData = List of additional data tables for Gene-centric pages. Must have "Gene" field. []
    - MakeClass = Generate classes from classorder LIST if missing from indata (go/keyword/xref(ppi)) []
    - MakePages = Types of pages to make [front,gene,go,mcode,class,interactome]
    - PPComplex = List of different evidence codes for special MCODE analyses ['Complex']
    - PPExtra = List of extra genes added from expanded MCODE clusters []
    - Special = Execute special analysis code for specific applications []
    - StyleSheets = List of CSS files to use ['../example.css','../redwards.css']
    - UniReal = Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]

    Dict:dictionary
    - GeneMap = Dictionary of mapped gene identifier to Gene
    - ClassGenes = Lists of genes for each class

    Obj:RJE_Objects
    - DB = RJE_DB Database object
    - PPI = RJE_PPI PPI Object
    - GO = RJE_GO GO Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['InData','GeneField','GeneClass','ClassCol','XRefData','Pairwise','HTMLPath','TitleText',
                         'UniFake','UniPath','GOData','GABLAM','PageHead']
        self.optlist = ['AddClass','FakeHTML','MakePNG','MultiClass','XGMML','AlphTabs','XRefTab','FillCol','NoBots',
                        'UsePos','UpdatePos']
        self.statlist = ['Border','GOPages','PPExpand','PNGMax']
        self.listlist = ['StyleSheets','UniReal','MakePages','GeneData','Combine','ClassOrder','MakeClass','Special',
                         'AddPPI','PPComplex','PPExtra']
        self.dictlist = ['GeneMap','ClassGenes']
        self.objlist = ['DB','GO','PPI','XRef']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'HTMLPath':rje.makePath('./html/'),'TitleText':'titletext.tdt','GeneField':'Gene','GeneClass':'Class',
                      'InfoText':'This data has not yet been published and should not be used without permission.'})
        self.setStat({'Border':2,'GOPages':5,'PPExpand':1,'PNGMax':2000,'MaxGO':500})
        self.setOpt({'FakeHTML':True,'MakePNG':True,'MultiClass':True,'Iridis':False,'XRefTab':True,'SVG':True,
                     'FillCol':True,'PPExtra':False,'UsePos':True,'UpdatePos':False,'AddClass':False})
        self.list['StyleSheets'] = ['../example.css','../slimhtml.css']
        self.list['MakePages'] = string.split('front,gene,go,mcode,class,interactome,expand',',')
        self.list['UniReal'] = string.split('AC,GN,RC,RX,CC,DR,PE,KW,FT',',')
        self.list['MakeClass'] = ['keyword']
        self.list['PPComplex'] = ['Complex']
        self.list['DropFields'] = ['EntrezCheck']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'info',['GeneField','GeneClass','PageHead'])
                self._cmdReadList(cmd,'file',['InData','ClassCol','XRefData','Pairwise','TitleText','GOData','GABLAM'])
                self._cmdReadList(cmd,'path',['DataPath','HTMLPath'])
                self._cmdReadList(cmd,'int',['Border','GOPages','PPExpand','PNGMax','MaxGO'])
                self._cmdReadList(cmd,'opt',['FakeHTML','MakePNG','MultiClass','XGMML','AlphTabs','XRefTab','SVG',
                                             'FillCol','NoBots','PPExtra','UsePos','UpdatePos','AddClass'])
                self._cmdReadList(cmd,'list',['UniReal','MakePages','Combine','ClassOrder','Special','PPComplex',
                                              'MakeClass','DropFields'])
                self._cmdReadList(cmd,'glist',['StyleSheets','GeneData','AddPPI'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.info['GABLAM'].lower() in ['','none']: self.info['GABLAM'] = '%s.gablam.tdt' % self.info['Basefile']
        if self.info['GOData'].lower() in ['','none']: self.info['GOData'] = '%s.go.tdt' % self.info['Basefile']
        if self.info['ClassCol'].lower() in ['','none']: self.info['ClassCol'] = '%s.col.tdt' % self.info['Basefile']
        if self.info['PageHead'].lower() in ['','none']: self.info['PageHead'] = 'HAPPI Analysis of %s' % self.info['Basefile']
        self.stat['PPExpand'] = max(0,self.stat['PPExpand'])
        self.list['MTypes'] = ['MCODE']
        if self.stat['PPExpand']: self.list['MTypes'] += ['Expanded']
        self.list['MTypes'] += self.list['PPComplex']
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            genelist = self.list['Genes']
            if not genelist: self.errorLog('No Gene List. Cannot run HAPPI.',printerror=False); return False
            try: self.opt['NoForks'] = self.opt['Win32'] or self.opt['NoForks'] or self.stat['Forks'] < 1
            except: self.opt['NoForks'] = True
            ### ~ [2] ~ Makes HTML Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forkpage = [] 
            ## ~ [2a] ~ Front Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'front' in self.list['MakePages']: self.frontPages()
            ## ~ [2b] ~ Make Gene Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'gene' in self.list['MakePages']:
                gx = 0.0; gtot = len(genelist + self.list['PPExtra'])
                for gene in genelist + self.list['PPExtra']:
                    self.progLog('\r#GENE','Making Gene Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    if 'slim' in self.list['Special']: hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'slim/'),rje_slim.slimFromPattern(gene))
                    else: hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
                    if not self.opt['Force'] and rje_html.checkHTML(hfile): continue
                    if self.opt['NoForks']: self.genePage(gene)
                    else: forkpage.append((gene,'gene'))
                if self.opt['NoForks']: self.printLog('\r#GENE','%s Gene HTML Pages made.' % rje.integerString(gtot))
                else: self.printLog('\r#GENE','%s Gene HTML Pages to fork.' % rje.iLen(forkpage))
            ## ~ [2c] GO Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'go' in self.list['MakePages']:
                gdb = self.db().getTable('GO')
                gx = 0.0; gtot = len(gdb.index('GO_Desc')); gpx = 0
                for go in rje.sortKeys(gdb.index('GO_Desc')):
                    self.progLog('\r#GO','Making GO Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    goid = gdb.dataList(gdb.indexEntries('GO_Desc',go),'GO_ID')[0]
                    if goid not in self.dict['GOPages']: continue
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
                    if not self.opt['Force'] and rje_html.checkHTML(hfile): continue
                    gpx += 1
                    if self.opt['NoForks']: self.goPage(go,goid)
                    else: forkpage.append((goid,go))
                if self.opt['NoForks']: self.printLog('\r#GO','%s GO HTML Pages made.' % rje.integerString(gpx))
                else: self.printLog('\r#GO','%s GO HTML Pages to fork.' % rje.integerString(gpx))
            try: self.forkPage(forkpage)
            except SystemExit: os._exit(0)
            except: raise
            ## ~ [2d] ~ Make cluster pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pobj = self.obj['PPI']
            for mtype in self.list['MTypes']:
                if 'mcode' in self.list['MakePages']:   # No longer each mtype.lower()
                    cdb = pobj.db(mtype)
                    for id in rje.sortKeys(cdb.index(mtype)): self.mcodePage(id,mtype)
                    self.printLog('\r#%s' % mtype.upper()[:5],'%s %s HTML Pages made.' % (rje.integerString(len(cdb.index(mtype))),mtype))
            if 'hm' in self.list['Special']:
                cdb = pobj.db('AvFnet')
                for id in rje.sortKeys(cdb.index('AvFnet')): self.mcodePage(id,'AvFnet')
                self.printLog('\r#AvF','%s AvFnet HTML Pages made.' % rje.integerString(len(cdb.index('AvFnet'))))
            ### ~ [3] ~ Make Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forkpng = [] 
            if self.opt['MakePNG']:
                pobj = self.obj['PPI']
                for mtype in self.list['MTypes']:
                    if 'mcode' in self.list['MakePages']:   # No longer each mtype.lower()
                        cdb = pobj.db(mtype)
                        for id in rje.sortKeys(cdb.index(mtype)):
                            if self.opt['NoForks']: self.interactomePNG(id,mtype.lower())
                            else: forkpng.append((id,mtype.lower()))
                if 'hm' in self.list['Special']:
                    cdb = pobj.db('AvFnet')
                    for id in rje.sortKeys(cdb.index('AvFnet')):
                        if self.opt['NoForks']: self.interactomePNG(id,'avfnet')
                        else: forkpng.append((id,'avfnet'))
                if 'class' in self.list['MakePages']:     # interactome_png
                    for pclass in rje.sortKeys(self.dict['ClassGenes']):
                        if self.opt['NoForks']: self.interactomePNG(pclass,'geneclass')
                        else: forkpng.append((pclass,'geneclass'))
                    if self.list['Combine']:
                        if self.opt['NoForks']: self.interactomePNG('combined','geneclass')
                        else: forkpng.append(('combined','geneclass'))
                if 'interactome' in self.list['MakePages']:     # interactome_png
                    gtot = len(genelist); gx = 0
                    for gene in genelist:
                        if self.opt['NoForks']: 
                            self.interactomePNG(gene,'interactome'); gx += 1
                            self.printLog('#IMG','Gene interactome %s of %s made' % (gx,gtot),log=False)
                        else: forkpng.append((gene,'interactome'))
                if 'go' in self.list['MakePages']:     # interactome_png
                    gtot = len(self.dict['GOPages']); gx = 0
                    for goid in rje.sortKeys(self.dict['GOPages']):
                        if self.opt['NoForks']:
                            self.interactomePNG(goid,'go'); gx += 1
                            self.printLog('#IMG','GO interactome %s of %s made' % (gx,gtot),log=False)
                        else: forkpng.append((goid,'go'))
                if 'interactome' in self.list['MakePages']:     # interactome_png
                    gtot = len(self.list['PPExtra']); gx = 0
                    for gene in self.list['PPExtra']:
                        if self.opt['NoForks']: 
                            self.interactomePNG(gene,'interactome'); gx += 1
                            self.printLog('#IMG','Extra Gene interactome %s of %s made' % (gx,gtot),log=False)
                        else: forkpng.append((gene,'interactome'))
                try: self.forkPNG(forkpng)
                except SystemExit: os._exit(0)
                except: raise
            ## ~ [3a] UniFake Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if 'fake' in self.list['MakePages']: self.uniFakeHTML()     #??#
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def forkPage(self,forkpng):  ### Fork out HTML page generation (lists of (id,type) tuples)
        '''Fork out HTML page generation (lists of (id,type) tuples).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not forkpng: return
            self.makeIndices()
            forkx = self.stat['Forks']      # Number of forks to have running at one time
            forks = {}      # Dictionary of active forks {pid:basefile}
            killforks = self.stat['KillForks']      # Time in seconds to wait after main thread has apparently finished
            killtime = time.time()
            pngx = len(forkpng)
            ### ~ [1] ~ Cycle through forkpng tuples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while forks or forkpng:
                ## ~ [1a] ~ Add more forks if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while forkpng and (len(forks) < forkx):     # Add more forks
                    (id,mtype) = forkpng.pop(0)
                    if mtype not in ['gene']: go = mtype; mtype = 'go'
                    basefile = 'tmp_%s-%s' % (mtype,id)
                    if mtype == 'gene' and 'slim' in self.list['Special']: basefile = 'tmp_%s-%s' % (mtype,rje_slim.slimFromPattern(id))
                    forkcmd = self.cmd_list + ['i=-1','log=%s.log' % basefile,'errorlog=None']
                    newpid = os.fork() 
                    if newpid == 0: # child
                        self.opt['Child'] = True
                        if mtype == 'gene': self.genePage(id)
                        if mtype == 'go': self.goPage(go,id)
                        os._exit(0)    # Exit process 
                    elif newpid == -1: self.errorLog('Problem forking %s.' % id,printerror=False)  
                    else: forks[newpid] = basefile            
                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(0.1)       # Sleep for 1s 
                forklist = self._activeForks(forks.keys())
                if len(forklist) != len(forks):
                    self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    for pid in rje.sortKeys(forks):    # Go through current forks
                        if pid not in forklist:
                            fork = forks.pop(pid)
                            killtime = time.time()  # Reset killtime - still doing stuff
                            rje.fileTransfer(fromfile='%s.log' % fork,tofile=self.log.info['LogFile'],deletefrom=True,append=True)
                    self.progLog('\r#FORK','%s of %s pages to fork     ' % (rje.iLen(forkpng),rje.iStr(pngx)))
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
            self.printLog('\r#FORK','End of %s HTML page forking.' % (rje.iStr(pngx)))
        except SystemExit: os._exit(0)
        except: self.errorLog('%s.forkPage error' % self)
#########################################################################################################################
    def ppi(self): return self.obj['PPI'].ppi()
    def png(self): return {True:'svg',False:'png'}[self.opt['SVG']]
#########################################################################################################################
    def classOrder(self,classes=[]):   ### Returns classes in order
        '''Returns classes in order.'''
        main = self.db('Main')
        if not classes: classes = rje.sortKeys(main.index(self.info['GeneClass']))
        if self.opt['MultiClass']: classes = string.split(string.join(classes,'-'),'-')
        order = []
        for gclass in self.list['ClassOrder']:
            if gclass in classes: order.append(gclass)
        for gclass in classes:
            if gclass not in order: order.append(gclass)
        return order
#########################################################################################################################
    def classGenes(self):   ### Generate self.dict['ClassGenes'] from main table
        '''Generate self.dict['ClassGenes'] from main table.'''
        main = self.db('Main')
        self.dict['ClassGenes'] = {}
        for iclass in rje.sortKeys(main.index(self.info['GeneClass'])):
            if self.opt['MultiClass']: gclasses = string.split(iclass,'-')
            else: gclasses = [iclass]
            igenes = main.indexDataList(self.info['GeneClass'],iclass,self.info['GeneField'])
            for gclass in gclasses:
                if gclass not in self.dict['ClassGenes']: self.dict['ClassGenes'][gclass] = igenes[0:]
                else: self.dict['ClassGenes'][gclass] = rje.listUnion(self.dict['ClassGenes'][gclass],igenes)
        for gclass in self.dict['ClassGenes']: self.dict['ClassGenes'][gclass].sort()
#########################################################################################################################
    def goIDfromName(self,go):  ### Returns GO_ID from GO_Desc
        '''Returns GO_ID from GO_Desc.'''
        try: return self.db('GO').indexEntries('GO_Desc',go)[0]['GO_ID']
        except:
            if self.opt['DeBug']: self.errorLog('')
            return None
#########################################################################################################################
    def makeClasses(self):    ### Makes missing classes, if appropriate
        '''Makes missing classes, if appropriate.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not self.list['MakeClass'] or not self.list['ClassOrder']: return
            main = self.db('Main')
            if self.opt['MultiClass']: indata = self.db('InData')
            else: indata = main
            indata.index(self.info['GeneClass'])
            dbxref = self.db('DBXRef'); tdb = self.db('TitleText')
            newclass = []
            self.list['MakeClass'] = rje.sortUnique(string.split(string.join(self.list['MakeClass']).lower()))
            explanations = {'keyword':'keyword in gene description and/or associated GO terms',
                            'desc':'keyword in gene description',
                            'gene':'gene identifier plus known interactors',
                            'xref':'database identifier cross-reference plus known interactors',
                            'go_desc':'GO term description','go_id':'GO term identifier'}
            ### ~ [1] ~ Make extra classes if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug(self.list['ClassOrder'])
            #self.deBug(rje.sortKeys(main.index(self.info['GeneClass'])))
            for gclass in self.list['ClassOrder']:
                if gclass in indata.index(self.info['GeneClass']): self.titleText('class',gclass); continue
                self.deBug(gclass)
                ## ~ [1a] ~ Get members for new classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cgenes = []; csource = []
                if 'keyword' in self.list['MakeClass']:
                    for go in self.db('GO').index('GO_Desc'):
                        if go.lower().find(gclass.lower()) >= 0:
                            cgenes += self.db('GO').indexDataList('GO_Desc',go,'Gene',sortunique=True)
                    for desc in dbxref.index('EnsDesc'):
                        if desc.lower().find(gclass.lower()) >= 0:
                            cgenes += dbxref.indexDataList('EnsDesc',desc,'Gene',sortunique=True)
                            if 'desc' in self.list['MakeClass']: csource = ['desc']
                    if cgenes: csource.append('keyword')
                elif 'desc' in self.list['MakeClass']:
                    for desc in dbxref.index('EnsDesc'):
                        if desc.lower().find(gclass.lower()) >= 0:
                            cgenes += dbxref.indexDataList('EnsDesc',desc,'Gene',sortunique=True)
                    if cgenes: csource.append('desc')
                if 'xref' in self.list['MakeClass']:
                    gene = self.gene(gclass)
                    self.deBug(gene)
                    if gene in dbxref.data():
                        cgenes = [gene]
                        if gene in self.ppi():
                            cgenes += rje.sortKeys(self.ppi()[gene])
                        if gene == gclass: csource.append('gene')
                        else: csource.append('xref')
                        self.deBug(cgenes)
                if 'go' in self.list['MakeClass']:
                    if go in self.db('GO').index('GO_Desc'):
                        cgenes += self.db('GO').indexDataList('GO_Desc',go,'Gene',sortunique=True)
                        csource.append('go_desc')
                    elif go in self.db('GO').index('GO_ID'):
                        cgenes += self.db('GO').indexDataList('GO_ID',go,'Gene',sortunique=True)
                        csource.append('go_id')
                cgenes = rje.sortUnique(cgenes)
                for gene in cgenes[0:]:
                    if gene not in main.dataKeys(): main.addEntry({self.info['GeneField']:gene})
                    if main.data()[gene][self.info['GeneClass']]:
                        if self.opt['MultiClass'] or 'add' in self.list['MakeClass']: main.data()[gene][self.info['GeneClass']] += '-%s' % gclass
                        else: cgenes.remove(gene)
                    else: main.data()[gene][self.info['GeneClass']] = gclass
                if not cgenes: self.printLog('#CLASS','No genes added for %s (%s)' % (gclass,string.join(self.list['MakeClass'],',')))
                else:
                    self.printLog('#CLASS','%s genes added for %s (%s)' % (rje.integerString(len(cgenes)),gclass,string.join(csource,',')))
                    newclass.append(gclass)
                    tkey = 'class\t%s' % gclass
                    if tkey not in tdb.data():
                        tdb.data()[tkey] = {'Page':'class','ID':gclass,'Title':'Genes with "%s" %s' % (gclass,explanations[csource[0]])}
                        for ctype in csource[1:]: tdb.data()[tkey]['Title'] += ' or %s' % explanations[ctype]
                    continue
            if newclass:
                self.printLog('#CLASS','%d new classes: %s' % (len(newclass),string.join(newclass,', ')))
                main.fillBlanks()
                main.index(self.info['GeneClass'],force=True)
            for gclass in self.list['ClassOrder']:
                tkey = 'mcode\t%s' % gclass; ckey = 'class\t%s' % gclass
                if tkey not in tdb.data(): tdb.data()[tkey] = {'Page':'mcode','ID':gclass,'Title':'No. of %s' % self.titleText('class',gclass)}
            tdb.saveToFile(tdb.info['Source'],delimit='\t',backup=self.opt['Backups'])
            self.classGenes()
        except: self.errorLog('HAPPI.makeClass() error')            
#########################################################################################################################
    def fillCol(self):  ### Fills out ClassCol dictionary for multiclass
        '''Fills out ClassCol dictionary for multiclass.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db('Main')
            ### ~ [1] ~ Fill out colour dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gclass in main.index(self.info['GeneClass']):
                if gclass in self.dict['ClassCol']: continue
                gcol = []
                for cclass in self.classOrder(string.split(gclass,'-')):
                    if cclass in self.dict['ClassCol']: gcol.append(self.dict['ClassCol'][cclass])
                if not gcol: continue
                self.dict['ClassCol'][gclass] = gcol[0]
                for col in gcol[1:]:
                    if gcol.count(col) > gcol.count(self.dict['ClassCol'][gclass]): self.dict['ClassCol'][gclass] = col
            #self.deBug(self.dict['ClassCol'])
        except: self.errorLog('HAPPI.fillCol() error')            
#########################################################################################################################
    def makeIndices(self):  ### Makes all the database indices used for HTML page generation (for forking)
        '''Makes all the database indices used for HTML page generation (for forking).'''
        main = self.db('Main')
        main.index(self.info['GeneClass'])
        main.index(self.info['GeneField'])
        indata = self.db().getTable('InData')
        if not indata: indata = main
        indata.index(self.info['GeneField'])
        try: self.db().getTable('full').index('Pattern')
        except: pass
        try:
            gdb = self.db().getTable('GO')
            gdb.index('GO_Desc'); gdb.index('GO_ID'); gdb.index('GO_Type'); gdb.index('GO_ID'); gdb.index('Gene')
        except: pass
        pobj = self.obj['PPI']
        for mtype in self.list['MTypes']: pobj.db(mtype).index(mtype)
        for dfile in self.list['GeneData']:
            try:
                name = string.split(dfile,'.')[-2]
                rawdb = self.db(name)
                if self.info['GeneField'] in rawdb.fields(): rawdb.index(self.info['GeneField'])
                elif 'Gene' in rawdb.fields(): rawdb.index('Gene')
                if 'slim' in self.list['Special'] and name == 'SlimPPI': rawdb.index('Hub')
            except: pass
        try:
            dbxref = self.db('DBXRef')
            for field in dbxref.fields(): dbxref.index(field)
        except: pass
        try:
            xdb = self.obj['PPI'].db('Node')
            for xfield in ['Uni','Seq']: xdb.index(xfield)
        except: pass
        try: self.db('ClassCol').index('Class')
        except: pass
#########################################################################################################################
    def setupTitleText(self):   ### Sets up titletext file and/or database table
        '''Sets up titletext file and/or database table.'''
        try:### ~ [0] Setup File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(self.info['TitleText']): tdb = self.db().addTable(self.info['TitleText'],['Page','ID'],name='TitleText')
            else:
                tdb = self.db().addEmptyTable('TitleText',['Page','ID','Title'])
                tdb.info['Source'] = self.info['TitleText']
            ### ~ [1] Populate with standard descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for standard in ['class	Description	Description of class/tab contents',
                             'class	Tab	Tab header',
                             'front	gene	Gene lists with source classes',
                             'tab go Gene Ontology Categories',
                             'gosig GO_Type Type of GO category',
                             'gosig GO_Desc GO category description',
                             'gosig GO_ID GO Category identifier',
                             'gosig N Number of genes in GO category',
                             'gosig p Uncorrected probability of N+ genes being in GO category',
                             'gosig bonf Bonferroni adjusted probability of N+ genes being in GO category',
                             'mcode	Complex	MCODE-predicted protein complex ID',
                             'mcode	Density	Node density of MCODE cluster',
                             'mcode	Edges	No. edges in MCODE cluster',
                             'mcode	Expanded	Expanded PPI MCODE subnetwork ID',
                             'mcode	MCODE	MCODE PPI subnetwork ID',
                             'mcode	Members	Proteins contained in MCODE cluster',
                             'mcode	Nodes	No. nodes in MCODE cluster',
                             'mcode	Rank	Rank of MCODE cluster',
                             'mcode	Score	MCODE Score (Nodes x Density) used for Ranking',
                             'mcode	Seed	Seed protein for MCODE cluster',
                             'ppi	description	Decription of interacting protein',
                             'ppi	evidence	Evidence (Database:experiment type) for protein interaction',
                             'ppi	gene	Gene Symbol for interacting protein',
                             'results	#	AutoID',
                             'results	%s	Source class for Gene' % self.info['GeneClass'],
                             'results	Gene	HGNC Gene Symbol',
                             'results	EnsEMBL	Cross-reference to Ensembl database',
                             'results	EnsLoci	Identifier for single-protein-per-gene Ensembl dataset',
                             'results	Entrez	Entrez Gene ID',
                             'results	EntrezCheck	Entrez Gene ID',
                             'results	HGNC	HGNC Gene Symbol',
                             'results	HPRD	Cross-reference to HPRD database',
                             'results	NRID	Non-redundant gene identifier mapped from Protein ID',
                             'results	OMIM	Cross-reference to OMIM database',
                             'results	Symbol	Gene symbol',
                             'results	UniProt	Cross-reference to UniProt database',
                             'results	Species	Species of gene sequence',
                             'xref	ensembl	Open EnsEMBL entry',
                             'xref	gene	Open Genecards entry',
                             'xref	hprd	Open HPRD entry',
                             'xref	omim	Open OMIM entry',
                             'xref	uniprot	Open UniProt entry',
                             'xreftab	Alias	Gene symbol aliases',
                             'xreftab	Desc	Description of identified protein',
                             'xreftab	EnsDesc	Ensembl gene description',
                             'xreftab	EnsEMBL	Cross-reference to Ensembl database',
                             'xreftab	EnsLoci	Identifier for single-protein-per-gene Ensembl dataset',
                             'xreftab	Entrez	Entrez Gene ID',
                             'xreftab	EntrezCheck	Entrez Gene ID',
                             'xreftab	Gene	Gene symbol for protein',
                             'xreftab	HGNC	HGNC Gene Symbol',
                             'xreftab	HPRD	Cross-reference to HPRD database',
                             'xreftab	NRID	Non-redundant gene identifier mapped from Protein ID',
                             'xreftab	OMIM	Cross-reference to OMIM database',
                             'xreftab	Symbol	Gene symbol',
                             'xreftab	UniProt	Cross-reference to UniProt database',
                             'xreftab	Species	Species of gene sequence']:
                title = string.split(standard)
                page = title[0]; id = title[1]; title = string.join(title[2:])
                tkey = '%s\t%s' % (page,id)
                if tkey not in tdb.data(): tdb.data()[tkey] = {'Page':page,'ID':id,'Title':title}
        except: self.errorLog('%s.setupTitleText error' % self)
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Main Database Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.setupTitleText()
            ## ~ [1a] PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppi = self.obj['PPI'] = rje_ppi.PPI(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
            if 'slim' in self.list['Special']: ppi.loadPairwisePPI(self.info['Pairwise'],clear=True,asdict=True,sym=True,evidence={'co-occurrence':'ignore'})   #,resave='%s.newppi.tdt' % self.info['Basefile'])
            else: ppi.loadPairwisePPI(self.info['Pairwise'],clear=True,asdict=True,sym=True)
            for pfile in self.list['AddPPI']:
                if 'slim' in self.list['Special']: ppi.loadPairwisePPI(pfile,clear=False,asdict=True,sym=True,evidence={'co-occurrence':'ignore'})  #,resave='%s.newppi.tdt' % self.info['Basefile'])
                elif pfile.lower() not in ['','none']: ppi.loadPairwisePPI(pfile,clear=False,asdict=True,sym=True)
            # This will populate PPI Table and create Node and Edge tables too
            ## ~ [1b] XRef table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = self.info['XRefData']
            if dfile.lower() not in ['none','']: dbxref = db.addTable(dfile,['Gene'],'All',name='DBXRef')
            #Gene    Entrez  HPRD    OMIM    UniProt EnsEMBL EnsLoci EnsDesc 
            #A1BG    1       00726   138670  P04217  ENSG00000121410 A1BG_HUMAN__ENSP00000263100     Alpha-1B-glycoprotein Precursor (Alpha-1-B glycoprotein)
            #x#for xref in string.split('A1BG    P04217  ENSG00000121410 A1BG_HUMAN__ENSP00000263100'): self.deBug(self.gene(xref))
            ## ~ [1c] GO Data and Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            go = self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            go.readGO()#; go.mapEnsGO()
            if not rje.exists(self.info['GOData']): self.makeGOData()
            #Gene    GO_ID   GO_Type GO_Desc
            if 'slim' in self.list['Special']:
                gdb = db.addTable(self.info['GOData'],['Pattern','GO_ID'],headers=string.split('Pattern GO_ID   Count   GO_Type GO_Desc'),name='GO')
                gdb.renameField('Pattern','Gene')
            else: gdb = db.addTable(self.info['GOData'],['Gene','GO_ID'],'All',name='GO')
            ## ~ [1d] Classes Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = self.info['InData']
            if self.opt['MultiClass']:
                if rje.exists(dfile):
                    indata = db.addTable(dfile,'All','All',name='InData')
                    for field in indata.fields(): self.titleText('results',field)
                else:
                    indata = db.addEmptyTable('InData',[self.info['GeneField'],self.info['GeneClass']],[self.info['GeneField']])
                    self.printLog('#DATA','No Input Data given. Will try to make ALL gene classes.')
                try:
                    main = db.addEmptyTable('Main',[self.info['GeneField'],self.info['GeneClass']],[self.info['GeneField']])
                    for entry in indata.entries():
                        gene = entry[self.info['GeneField']]
                        gclass = entry[self.info['GeneClass']]
                        if gene in main.dataKeys(): main.data()[gene][self.info['GeneClass']] = string.join(rje.sortUnique(string.split(main.data()[gene][self.info['GeneClass']],'-')+[gclass]),'-')
                        else: main.addEntry({self.info['GeneField']:gene,self.info['GeneClass']:gclass})
                        #x#self.deBug(main.data())
                    self.printLog('#MAIN','%s Gene entries from %s InData entries' % (rje.iStr(main.entryNum()),rje.iStr(indata.entryNum())))
                except KeyError:
                    self.errorLog('Problem loading main input data. Check field headers match GeneField ("%s") and GeneClass ("%s") settings' % (self.info['GeneField'],self.info['GeneClass']))
                    raise
                except: self.errorLog('Problem loading main input data.'); raise
            else:
                if rje.exists(dfile):
                    main = db.addTable(dfile,[self.info['GeneField']],'All',name='Main')
                    for field in main.fields(): self.titleText('results',field)
                else:
                    main = db.addEmptyTable('Main',[self.info['GeneField'],self.info['GeneClass']],[self.info['GeneField']])
                    self.printLog('#DATA','No Input Data given. Will try to make ALL gene classes.')
            if not main: main = db.addEmptyTable('Main',[self.info['GeneField'],self.info['GeneClass']],[self.info['GeneField']])
            self.makeClasses()
            if 'combined' in main.index(self.info['GeneClass']):
                self.errorLog('Cannot have Gene Class "combined". Rename and retry.', printerror=False); return False
            self.specialInData()    #!#
            genelist = self.list['Genes'] = self.geneList()
            # Add classes to PPI? (e.g. Host-Pathogen PPI)
            if self.opt['AddClass']:
                for entry in main.entries():
                    for gclass in string.split(entry[self.info['GeneClass']],'-'):
                        ppi.addPPI(gclass,entry[self.info['GeneField']],self.info['InData'])
                    if gclass not in genelist: genelist.append(gclass)
                genelist.sort()
            if not genelist: self.errorLog('No Gene List. Check input.',printerror=False); raise ValueError
            self.list['PPExtra'] = []
            if 'expand' in self.list['MakePages'] and self.opt['PPExtra']:
                px = 0; nx = 0
                for hub in genelist:
                    try:
                        for spoke in rje.sortKeys(self.ppi()[hub]):
                            if spoke not in self.list['PPExtra']: self.list['PPExtra'].append(spoke); px += 1
                    except: nx += 1
                self.printLog('#PPI','Added %s PPI genes to list for Gene Pages. (%d Input genes without PPI.)' % (rje.iStr(px),nx))
            if len(self.list['Combine']) == 1 and self.list['Combine'][0].lower() in  ['all','*']:
                self.list['Combine'] = self.classOrder()
            ## ~ [1e] Class Colours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = self.info['ClassCol']
            if rje.exists(dfile): cdb = db.addTable(dfile,['Class'],'All',name='ClassCol')
            else: cdb = None
            self.dict['ClassCol'] = {}
            if cdb:
                for entry in cdb.entries(): self.dict['ClassCol'][entry['Class']] = entry['Col']
                if self.opt['MultiClass'] and self.opt['FillCol']: self.fillCol()
            else:
                firstchoice = [5,8,10,13,15,17,19,20,27,28,39,43,44,45,47]
                for i in range(1,47):
                    if i not in firstchoice: firstchoice.append(i)
                while len(firstchoice) < len(main.index(self.info['GeneClass'])): firstchoice = firstchoice * 2
                self.dict['ClassCol']['Default'] = 12
                for gclass in self.classOrder(): self.dict['ClassCol'][gclass] = firstchoice.pop(0)
            ppi.addCol(default=12,coldict=self.dict['ClassCol'],ckey=self.info['GeneClass'],edb=main,addcol=48)
            ## ~ [1f] Additional Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for dfile in self.list['GeneData']:
                name = string.split(dfile,'.')[-2]
                try:
                    for field in db.addTable(dfile,'All','All',name=name).fields(): self.titleText('raw',field)
                except: self.errorLog('Probelm with %s gene data' % name)
            if self.opt['MultiClass']: self.list['GeneData'].append('InData.tdt')

            ### ~ [2] MCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            reformat = {}
            for f in string.split('Nodes Edges Rank'): reformat[f] = 'int'
            for f in string.split('Density Score'): reformat[f] = 'num'
            ## ~ [2a] MCODE clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for mtype in self.list['MTypes']:   #['MCODE','Complex','Expanded']:
                mfile = '%s.%s.tdt' % (self.info['Basefile'],mtype.lower())
                if os.path.exists(mfile) and self.opt['Force'] and self.yesNo('Remake %s?' % mfile): rje.backup(self,mfile,unlink=True,appendable=False)
                if os.path.exists(mfile):
                    cdb = ppi.db().addTable(mfile,mainkeys=['Seed'],datakeys='All',name=mtype)
                    cdb.dataFormat(reformat)
                    ppi.dict[mtype] = {}
                    for entry in cdb.entries():
                        try: ppi.dict[mtype][entry['Seed']] = string.split(entry['Members'],'|')
                        except: ppi.dict[mtype][entry['Seed']] = string.split(entry['Complex'],'|')
                else:
                    G = rje_ppi.subGraph(ppi.ppi(),self.geneList())
                    #if 'slim' in self.list['Special']: comtext = 'cloud'
                    #else: comtext = 'complex'
                    if mtype in self.list['PPComplex']: G = ppi.complexOnly(G,comtext=mtype.lower())
                    if mtype != 'MCODE' and self.stat['PPExpand']:
                        G = ppi.expandPPI(G,self.stat['PPExpand'],ppcomplex=mtype!='Expanded',comtext=mtype.lower())
                        G = rje_ppi.subGraph(ppi.ppi(),rje.sortKeys(G))
                    ppi.opt['Fluff'] = ppi.stat['Fluff'] > 0 and mtype == 'MCODE'
                    ppi.cmdMCODE(G,mtype)     # MCODE complexes saved to file and also stored in ppi.db('Complex') and ppi.dict['Complex'] = {seed:[]}
                    #ppi.dict[mtype] = ppi.dict.pop('Complex')
                    cdb = ppi.db(mtype)
                    #self.deBug(cdb)
                if 'Members' not in cdb.fields():
                    cdb.renameField('Complex','Members',log=True)
                    cdb.info['Name'] = mtype
                    cdb.list['Fields'].insert(0,mtype)
                    cdb.list['Fields'] += self.classOrder()
                    for entry in cdb.entries():
                        entry[mtype] = '%s_N%d' % (rje.preZero(entry['Rank'],cdb.entryNum()),entry['Nodes'])
                        for gclass in self.classOrder(): entry[gclass] = 0
                        for gene in ppi.dict[mtype][entry['Seed']]:
                            for gclass in self.geneClass(gene,split=True):
                                try: entry[gclass] += 1
                                except: pass
                    cdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],cdb.info['Name'].lower()))
                if self.opt['PPExtra']:
                    for entry in cdb.entries():
                        for gene in ppi.dict[mtype][entry['Seed']]:
                            if gene not in self.geneList() and gene not in self.list['PPExtra']: self.list['PPExtra'].append(gene)
            self.list['PPExtra'].sort()
            ## ~ [2d] Special complex analyses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'hm' in self.list['Special']: self.hmComplexes()
           
            ### ~ [3] GO Data and Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['GOPages'] = {}
            for id in rje.sortKeys(gdb.index('GO_ID')):
                self.dict['GOPages'][id] = []
                for gentry in gdb.indexEntries('GO_ID',id):
                    if gentry['Gene'] in genelist: self.dict['GOPages'][id].append(gentry['Gene'])
                if len(self.dict['GOPages'][id]) < self.stat['GOPages']: self.dict['GOPages'][id] = []
                if not self.dict['GOPages'][id]: self.dict['GOPages'].pop(id)
            self.printLog('#GO','%s GO IDs with %d+ genes' % (rje.integerString(len(self.dict['GOPages'])),self.stat['GOPages']))

            ### ~ [4] GABLAM Distance Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dfile = self.info['GABLAM']
            if rje.exists(dfile):
                self.obj['Dis'] = dis = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
                dis.loadFromDataTable(dfile)
                dis.forceSymmetry(method='mean')
                ## ~ [4a] Convert and add ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newmatrix = {}
                for obj1 in rje.sortKeys(dis.dict['Matrix']):
                    newmatrix[self.gene(obj1)] = {}
                    for obj2 in dis.dict['Matrix'][obj1]: newmatrix[self.gene(obj1)][self.gene(obj2)] = dis.dict['Matrix'][obj1][obj2]
                dis.dict['Matrix'] = newmatrix
                for gene in main.dataKeys():
                    if gene not in dis.dict['Matrix']: dis.dict['Matrix'][gene] = {}
                self.printLog('#DIS','GABLAM matrix loaded and converted to Gene Identifiers')
            else: self.obj['Dis'] = None

            ### ~ [5] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['HTMLPath'])
            subdirlist = ['gene','go','interactome_png','class_png','go_png']
            for mtype in self.list['MTypes']: subdirlist += [mtype.lower(),'%s_png' % mtype.lower()]
            if 'hm' in self.list['Special']: subdirlist += ['avfnet','avfnet_png']
            if 'slim' in self.list['Special']: subdirlist[0] = 'slim'
            for subdir in subdirlist:
                if self.opt['SVG']: rje.mkDir(self,rje.makePath(self.info['HTMLPath']+string.replace(subdir,'_png','_svg')))
                else: rje.mkDir(self,rje.makePath(self.info['HTMLPath']+subdir))
            self.printLog('#SETUP','HAPPI setup complete.'); return True
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def makeGOData(self):   ### Generates delimited GO File
        '''Generates delimited GO File.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfile = self.info['GOData']
            ghead = string.split('Gene    GO_ID   GO_Type GO_Desc')
            rje.delimitedFileOutput(self,gfile,ghead,rje_backup=True)
            godata = {}
            go = self.obj['GO']
            go.mapEnsGO(spec='HUMAN',gokey='EnsGO',fixhead=True)
            ### ~ [1] ~ Map and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ex,etot) = (0.0,len(go.dict['EnsGO']))
            for ensg in rje.sortKeys(go.dict['EnsGO']):
                self.progLog('\r#GO','Making delimited GO file (1/2): %.2f%%  ' % (ex/etot)); ex += 100.0
                gene = self.gene(ensg)  #; self.deBug('%s -> %s = %s' % (ensg,gene,go.dict['EnsGO'][ensg]))
                if gene not in godata: godata[gene] = {}
                for id in go.dict['EnsGO'][ensg]:
                    if id in godata[gene]: continue
                    godata[gene][id] = {'Gene':gene,'GO_ID':id,'GO_Type':go.go(id)['type'],'GO_Desc':go.go(id)['name']}
                #self.deBug('... %s' % godata)
            (ex,etot) = (0.0,len(godata))
            GFILE = open(gfile,'a'); delimit = rje.delimitFromExt(filename=gfile)
            for gene in rje.sortKeys(godata):
                self.progLog('\r#GO','Making delimited GO file (2/2): %.2f%%  ' % (ex/etot)); ex += 100.0
                for id in rje.sortKeys(godata[gene]):
                    outlist = []
                    for head in ghead: outlist.append(godata[gene][id][head])
                    rje.writeDelimit(GFILE,outlist,delimit)
            self.printLog('\r#GO','Making delimited GO file (2/2) finished.',log=False)
            GFILE.close()
            self.printLog('#GO','GO Table output (%s) complete' % gfile)
        except: self.errorLog('Problem during %s makeGOData()' % self)
#########################################################################################################################
    def geneList(self): return self.db('Main').dataKeys()
#########################################################################################################################
    def gene(self,xref):    ### Returns Gene for Given cross-reference
        '''Returns Gene for Given cross-reference.'''
        xdb = self.db('DBXRef')
        if xdb and xref in xdb.data(): self.dict['GeneMap'][xref] = xref
        if xref not in self.dict['GeneMap']:
            for xfield in ['UniProt','EnsEMBL','EnsLoci']:
                if not xdb or xfield not in xdb.fields(): continue
                try:
                    self.dict['GeneMap'][xref] = xdb.indexDataList(xfield,xref,'Gene')[0]; break
                except: pass
        if xref not in self.dict['GeneMap']:
            xdb = self.obj['PPI'].db('Node')
            for xfield in ['Uni','Seq']:
                if xfield not in xdb.fields(): continue
                try: self.dict['GeneMap'][xref] = xdb.indexDataList(xfield,xref,'Node')[0]; break
                except: pass
        if xref not in self.dict['GeneMap']:
            details = string.split(xref,'_')
            if details[0] != details[0].upper(): self.dict['GeneMap'][xref] = details[-1]
            else: self.dict['GeneMap'][xref] = string.join(details[:2],'_')
        if self.dict['GeneMap'][xref] in ['-','']: self.dict['GeneMap'][xref] = xref
        return self.dict['GeneMap'][xref]
#########################################################################################################################
    def geneClass(self,gene,missing='-',split=False):   ### Returns Class for Given Gene
        '''Returns Class for Given Gene.'''
        try:
            if split:
                if self.opt['MultiClass']: return string.split(self.db('Main').data()[gene][self.info['GeneClass']],'-')
                else: return [self.db('Main').data()[gene][self.info['GeneClass']]]
            return self.db('Main').data()[gene][self.info['GeneClass']]
        except: return missing
#########################################################################################################################
    def goLink(self,go,goid,frontpage=False):  ### Returns GO link text
        '''Returns GO link text.'''
        if goid not in self.dict['GOPages']: return '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s" title="Browse AmiGO">%s</a>' % (goid,go)
        if frontpage: return '<a href="./go/%s.html" target="_blank" title="GO results page">%s</a>' % (goid,go)
        else: return '<a href="../go/%s.html" title="GO results page">%s</a>' % (goid,go)
#########################################################################################################################
    def geneLink(self,gene,frontpage=False,altgene=None):    ### Returns linked Gene text
        if not altgene: altgene = gene
        if gene in self.list['Genes']:
            if 'slim' in self.list['Special']: return slimLink(gene,frontpage)
            if frontpage: return '<a href="./gene/%s.html" target="_blank" title="%s results page">%s</a>' % (gene,gene,altgene)
            else: return '<a href="../gene/%s.html" title="%s results page">%s</a>' % (gene,gene,altgene)
        elif gene in self.list['PPExtra']:
            if 'slim' in self.list['Special']: return '<i>%s</i>' % slimLink(gene,frontpage)
            if frontpage: return '<a href="./gene/%s.html" target="_blank" title="%s information page"><i>%s</i></a>' % (gene,gene,altgene)
            else: return '<a href="../gene/%s.html" title="%s information page"><i>%s</i></a>' % (gene,gene,altgene)
        else: return altgene
#########################################################################################################################
    ### <3> ### Front pages & general HTML code                                                                         #
#########################################################################################################################
    def titleText(self,page,id):    ### Returns title (mouseover text) for Page/ID combo
        '''Returns title (mouseover text) for Page/ID combo.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tdb = self.db().getTable('TitleText')
            if not tdb:
                if os.path.exists(self.info['TitleText']): tdb = self.db().addTable(self.info['TitleText'],['Page','ID'],name='TitleText')
                else:
                    tdb = self.db().addEmptyTable('TitleText',['Page','ID','Title'])
                    tdb.info['Source'] = self.info['TitleText']
                    open(self.info['TitleText'],'w').write('%s\n' % string.join(['Page','ID','Title'],'\t'))
            tkey = '%s\t%s' % (page,id)
            ### ~ [1] ~ Easy case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if tkey in tdb.data(): return tdb.data()[tkey]['Title']
            ### ~ [2] ~ Add new case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if page == 'gosig' and id in self.classOrder(): title = self.classDef(id) 
            elif page == 'mcode' and id in self.classOrder(): title = 'No. of "%s"' % self.classDef(id)
            elif page == 'raw' and tdb.data().has_key('results\t%s' % id): title = tdb.data()['results\t%s' % id]['Title']
            elif self.opt['Test']: title = '%s|%s' % (page,id)
            else: title = rje.choice('New mouseover (title) text for %s|%s title?' % (page,id),confirm=True)
            tdb.data()[tkey] = {'Page':page,'ID':id,'Title':title}
            open(self.info['TitleText'],'a').write('%s\n' % string.join([page,id,title],'\t'))
            return title
        except: self.errorLog('TitleText problem (%s|%s)' % (page,id)); return ''
#########################################################################################################################
    def classDef(self,c): return self.titleText('class',c)
#########################################################################################################################
    def frontPages(self): self.frontPage(); self.fullPage(); self.keyPNG()
#########################################################################################################################
    def frontPage(self):    ### Generates HTML front page for analysis
        '''Generates HTML front pages for analysis.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = self.info['HTMLPath'] + 'index.htm'
            if not self.opt['Force'] and rje_html.checkHTML(hfile): return self.printLog('#HTML','Found main %s page.' % hfile)
            html = htmlHead(self.info['PageHead'],self.list['StyleSheets'],frontpage=True,nobots=self.opt['NoBots'])
            html += '\n<table border="0" width="100%%"><tr>\n<td width="60%%"><h1>%s</h1>\n' % self.info['PageHead']
            html += '<br><p><i>This data has not yet been published and should not be used without permission.</i></p></td>\n'
            html += '<td width="40%" valign="top" align="right">\n'
            #html += '<a href="http://www.rcsi.ie/"><img src="./resources/rcsi.png" height="100" alt="RCSI" title="The Royal College of Surgeons in Ireland"></a>\n'
            html += '<a href="http://www.personal.soton.ac.uk/re1u06/"><img src="./resources/SBS_100.png" height="100" alt="RJE Homepage" title="RJE Homepage"></a>\n'
            html += '</td></tr></table>\n\n'
            head_html = html[0:]
            main = self.db().getTable('Main'); indb = self.db().getTable('InData')
            if not indb: indb = main
            dbxref = self.db().getTable('DBXRef')
            jtxt = ' ~ '
            fronttab = []          
            ### ~ [2] ~ Classes Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            classdata = []
            for c in self.classOrder(): classdata.append({'Tab':c,'Description':self.classDef(c)})
            if self.list['Combine']: classdata.append({'Tab':'Combined','Description':'Combined %s PPI' % string.join(self.list['Combine'],' + ')})
            for mtype in self.list['MTypes']:
                if mtype == 'MCODE': classdata.append({'Tab':mtype,'Description':'MCODE clusters'})
                else: classdata.append({'Tab':mtype,'Description':'%s MCODE clusters' % mtype})
            if 'hm' in self.list['Special']: classdata.append({'Tab':'AvFnet','Description':'Subnetworks of AApep or FFpep binders.'})
            chtml = '<H2>Gene Class/Tab overview for %s</H2>\n' % self.info['PageHead']
            chtml += '<table><tr valign="top"><td width="30%" border=0>'
            chtml += '<p>Click <a href="./full.htm">here</a> to access data by gene/GO class.</p>\n'
            chtml += self.resultTableHTML(['Tab','Description'],classdata,asdict=False,ttype='class')
            if self.db('ClassCol'): keyheight = int((len(self.db('ClassCol').index('Class'))+11)/2.0) * 50
            else: keyheight = int((len(self.dict['ClassCol'])+11)/2.0) * 50
            if self.opt['SVG']: chtml += '</td><td width="70%%">%s\n</td></tr></table>' % self.svgCode('./class_key.svg','Colour key for gene classes in images')
            else: chtml += '</td><td width="70%%"><img src="./class_key.png" alt="Class Key not yet made" title="Colour key for gene classes in images" height=%d>\n</td></tr></table>' % keyheight
            fronttab.append(('Gene Classes',chtml,'Gene Class/Tab overview for %s' % self.info['PageHead']))
            ### ~ [3] ~ Class Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for c in self.classOrder():
                tablist = []
                chtml = '<H2>%s</H2>\n' % self.classDef(c)
                genehtml = []; xref = []; cdata = []
                for gene in self.dict['ClassGenes'][c]:
                    genehtml.append(self.geneLink(gene,True))
                    if dbxref and gene in dbxref.data(): xref.append(dbxref.data()[gene])
                    for entry in indb.indexEntries(self.info['GeneField'],gene):
                        if entry[self.info['GeneClass']] == c: cdata.append(entry)
                tablist.append(('Genes',string.join(genehtml,jtxt),'%s Genes' % c))
                if dbxref and self.opt['XRefTab']: tablist.append(('XRef',string.replace(self.resultTableHTML(dbxref.fields(),xref,asdict=False,ttype='xreftab'),'href="../','target="_blank" href="./'),'%s database cross-references' % c))
                tablist.append(('Data',string.replace(self.resultTableHTML(indb.fields(),cdata,asdict=False),'href="../','target="_blank" href="./'),'%s Summary' % c))
                ## ~ [2b] ~ Interactome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                intpng = './geneclass_%s/%s.geneclass.png' % (self.png(),c)
                if self.opt['SVG']:
                    svglink = './geneclass_svg/%s.geneclass.svg' % c
                    title = 'Graphical representation of %s PPI. Larger interactomes may not be clear. See tables for details.' % c
                    ctab = '%s\n' % self.svgCode(svglink,title)
                else:
                    alt = '%s interactome image not yet made' % c
                    title = 'Graphical representation of %s PPI. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % c
                    ctab = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (intpng,intpng,alt,title)
                if self.opt['XGMML']: ctab = '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n%s' % (intpng[:-4],ctab)
                tablist.append(('Interactome','<h2>%s Interactome</h2>\n%s' % (c,ctab),title))
                #!# Do we want MCODE Complexes within each Class? #!#
                ## ~ [2c] ~ Finish Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fronttab.append((c,'%s%s' % (chtml,tabberHTML('%s-tabs' % c,tablist,level=1)),self.classDef(c)))
            if self.list['Combine']:
                chtml = '<H2>Combined %s PPI</H2>\n' % string.join(self.list['Combine'],' + ')
                intpng = './geneclass_%s/combined.geneclass.png' % self.png()
                if self.opt['XGMML']: chtml += '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                if self.opt['SVG']:
                    svglink = './geneclass_svg/combined.geneclass.svg'
                    title = 'Graphical representation of combined PPI. Larger interactomes may not be clear. See tables for details.'
                    fronttab.append(('Combined','%s\n%s\n' % (chtml,self.svgCode(svglink,title)),title))
                else:
                    alt = 'Combined interactome image not yet made'
                    title = 'Graphical representation of combined PPI. Click to open in new window. Larger interactomes may not be clear. See tables for details.'
                    fronttab.append(('Combined','%s\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (chtml,intpng,intpng,alt,title),title))
            ### ~ [4] ~ MCODE Complex Page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # MCODE complexes saved to file and also stored in ppi.db('Complex') and ppi.dict['Complex'] = {seed:[]}
            # ['ID','Seed','Nodes','Edges','Density','Score','Complex','Rank']
            for mtype in self.list['MTypes']:
                cdb = self.obj['PPI'].db(mtype)
                cdata = {}
                for entry in cdb.entries(): cdata[entry[mtype]] = entry
                if mtype == 'MCODE': mtxt = 'MCODE predicted clusters'
                else: mtxt = '%s MCODE predicted clusters' % mtype
                chtml = '<H2>%s</H2>\n' % mtxt
                chtml += string.replace(self.resultTableHTML(cdb.fields(),cdata,ttype='mcode'),'href="../','target="_blank" href="./')
                fronttab.append((mtype,chtml,mtxt))
            if 'hm' in self.list['Special']:
                cdb = self.obj['PPI'].db('AvFnet')
                cdata = {}
                for entry in cdb.entries(): cdata[entry['AvFnet']] = entry
                chtml = '<H2>AvF subnetworks</H2>\n'
                chtml += string.replace(self.resultTableHTML(cdb.fields(),cdata,ttype='mcode'),'href="../','target="_blank" href="./')
                fronttab.append(('AvFnet',chtml,'Subnetworks of AApep or FFpep binders.'))
            ### ~ [5] ~ Front page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html += tabberHTML('Main',fronttab,level=0)
            html += htmlTail()
            open(hfile,'w').write(html); open('%sl' % hfile,'w').write(html)
            self.printLog('#HTML','Made main %s page.' % hfile)            
        except: self.errorLog('frontPage Error')
#########################################################################################################################
    def geneTabs(self,genelist):    ### Returns A-Z tablisting for gene list (for tabberHTML('Gene',genetab,level=1))
        '''Returns A-Z tablisting for gene list (for tabberHTML('Gene',genetab,level=1)).'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not genelist: self.errorLog('No Genelist for Gene Tabs!',printerror=False); raise ValueError
            genetitle = {}; genedict = {}; genetab = []; jtxt = '~'
            for x in '%s_' % string.ascii_uppercase:
                if x == '_':
                    genetitle[x] = 'Genes not beginning A-Y'
                    genedict[x] = ['<h2>Genes not beginning A-Y</h2><p>\n']
                else:
                    genetitle[x] = 'Genes beginning with %s' % x
                    genedict[x] = ['<h2>Genes beginning %s...</h2><p id="nonsig">\n' % x]
            ### ~ [2] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gene in genelist:
                x = gene[0]
                gclass = self.geneClass(gene,'')
                if gclass: gtxt = '%s (%s)' % (self.geneLink(gene,frontpage=True),gclass)
                else: gtxt = gene
                if x not in string.ascii_uppercase: x = '_'
                gtitle = '%s summary page' % gene
                genedict[x].append('<a href="./gene/%s.html" target="_blank" title="%s">%s</a>' % (gene,gtitle,gtxt))
            ### ~ [3] ~ Make Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for x in '%s_' % string.ascii_uppercase:
                if len(genedict[x]) > 1: 
                    if x == '_': genetab.append(('*','%s\n</p>' % string.join(genedict[x],jtxt),genetitle[x]))
                    else: genetab.append((x,'%s\n</p>' % string.join(genedict[x],jtxt),genetitle[x]))
            return genetab
        except: self.errorLog('geneTabs Error')
#########################################################################################################################
    def fullPage(self):     ### Generates HTML full gene/GO lists for analysis
        '''Generates HTML full gene/GO lists for analysis.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = self.info['HTMLPath'] + 'full.htm'
            if not self.opt['Force'] and rje_html.checkHTML(hfile): return self.printLog('#HTML','Found main %s page.' % hfile)
            html = htmlHead(self.info['PageHead'],self.list['StyleSheets'],frontpage=True,nobots=self.opt['NoBots'])
            html += '\n<table border="0" width="100%%"><tr>\n<td width="60%%"><h1>%s</h1>\n' % self.info['PageHead']
            html += '<br><p><i>This data has not yet been published and should not be used without permission.</i></p></td>\n'
            html += '<td width="40%" valign="top" align="right">\n'
            #html += '<a href="http://www.rcsi.ie/"><img src="./resources/rcsi.png" height="100" alt="RCSI" title="The Royal College of Surgeons in Ireland"></a>\n'
            html += '<a href="http://www.personal.soton.ac.uk/re1u06/"><img src="./resources/SBS_100.png" height="100" alt="RJE Homepage" title="RJE Homepage"></a>\n'
            html += '</td></tr></table>\n\n'
            head_html = html[0:]
            main = self.db().getTable('Main')
            dbxref = self.db().getTable('DBXRef')
            jtxt = ' ~ '
            ## ~ [1a] ~ Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            genetab = self.geneTabs(self.list['Genes'])
            ## ~ [1b] ~ GO Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb = self.db().getTable('GO'); gotab = []
            got = {'bp':'Biological Process','cc':'Cellular Component','mf':'Molecular Function'}
            for type in ['bp','cc','mf']:
                gosub = {}
                for x in string.ascii_uppercase: gosub[x] = ['<h3>%s GO terms starting %s</h3>' % (got[type],x)]
                for id in gdb.indexDataList('GO_Type',type,'GO_ID'):
                    if id not in self.dict['GOPages']: continue
                    go = gdb.indexDataList('GO_ID',id,'GO_Desc')[0]
                    if not go: continue
                    ix = len(self.dict['GOPages'][id])
                    i = 0
                    while go[i].upper() not in gosub: i += 1
                    x = go[i].upper()
                    gosub[x].append('%s (%d)' % (self.goLink(go,id,frontpage=True),ix))
                gohtml = []
                for x in string.ascii_uppercase:
                    gosub[x].sort()
                    if len(gosub[x]) > 1: gohtml.append((x,string.join(gosub[x],'\n ~ '),'%s GO terms starting %s' % (got[type],x)))
                #gotab.append(('GO_%s' % type.upper(),string.join(gohtml,'\n ~ ')))
                if gohtml: gotab.append(('GO_%s' % type.upper(),tabberHTML('GO_%s' % type.upper(),gohtml,level=1),'%s GO terms' % got[type]))
            ### ~ [2] ~ Main code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fulltab = [('Genes',tabberHTML('Gene',genetab,level=1),self.titleText('front','gene'))] + gotab
            fullhtml = head_html
            fullhtml += '<p>Full results data can be accessed via the links on these pages.\n'
            fullhtml += 'Click <a href="./index.htm">here</a> to return to the summary data page.</p>'
            fullhtml += tabberHTML('Full',fulltab,level=0)
            fullhtml += htmlTail()
            open(hfile,'w').write(fullhtml)
            self.printLog('#HTML','Made main %s page.' % hfile)
        except: self.errorLog('fullPage Error')
#########################################################################################################################
    ### <4> ### Tables pages                                                                                            #
#########################################################################################################################
    def specialGeneData(self,name,gene):   ### Returns field to be used for extracting gene data
        '''Returns field to be used for extracting gene data.'''
        rawdb = self.db(name)
        if self.info['GeneField'] in rawdb.fields(): gfield = self.info['GeneField']
        elif 'Gene' in rawdb.fields(): gfield = 'Gene'
        else: gfield = None
        if 'slim' in self.list['Special'] and name == 'SlimPPI': gfield = 'Hub'
        if 'hm' in self.list['Special'] and name == 'raw':
            rawdata = {}
            if gene in self.db('Main').data():
                for ensp in string.split(self.db('Main').data()[gene]['Source'],'; '):
                    esub = rawdb.subset('Identifier',ensp)
                    for dkey in esub: rawdata[dkey] = esub[dkey]
        else:
            if gfield: rawdata = rawdb.subset(gfield,gene)
            else: rawdata = None
        return rawdata
#########################################################################################################################
    def specialGeneSummary(self,gene):  ### Returns HTML for gene summary data table on Gene Page
        '''Returns HTML for gene summary data table on Gene Page.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db().getTable('Main')
            html = [self.resultTableHTML(main.fields(),main.subset(self.info['GeneField'],gene))]
            ### ~ [1] ~ Application specific data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ H&M Summary Data Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'hm' in self.list['Special']:
                droplist = []
                for field in main.fields():
                    if field not in ['NRID','Identifier','AA','enrAA','FF','enrFF','Class','PPI','Source','AvF','MeanEnr']: droplist.append(field)
                fwidth = {'NRID':120,'Identifier':150,'Desc':250,'AA':50,'FF':50,'enrAA':50,'enrFF':50,'Class':50,'PPI':50}
                html = [self.resultTableHTML(main.fields(),main.subset('NRID',gene),drop=droplist,width=fwidth),'<br>']
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
                    html += [self.resultTableHTML(thead,main.subset('NRID',gene),drop=droplist,width=fwidth,ttype='raw'),'<br>']
            ## ~ [1b] ~ SLiM Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'slim' in self.list['Special']:
                slim = gene
                mdb = self.db().getTable('full')
                odb = self.db().getTable('Occurrence')
                genelist = []; domlist = []; randlist = []
                for entry in mdb.indexEntries('Pattern',slim):
                    if entry['Dataset'][-3:] == 'dom': domlist.append('%s (%s)' % (entry['Hub'],rje_slim.expectString(string.atof(entry['Sig']))))
                    elif entry['Dataset'][:4] in ['rseq','rupc']: randlist.append('%s (%s)' % (entry['Dataset'],rje_slim.expectString(string.atof(entry['Sig']))))
                    else: genelist.append('%s (%s)' % (entry['Hub'],rje_slim.expectString(string.atof(entry['Sig']))))
                genelist.sort(); domlist.sort()
                dfont = '<FONT SIZE=4 COLOR=#014359>'
                html = ['<table width="100%">']
                if genelist: html += ['<tr><td><p>%sGene Hubs (%d): </FONT>%s</p></td></tr>' % (dfont,len(genelist),string.join(genelist,' ~ '))]
                else: html += ['<tr><td><p>%sGene Hubs: </FONT><i>None.</i></p></td></tr>' % (dfont)]
                if domlist: html += ['<tr><td><p>%sDomain Hubs (%d): </FONT>%s</p></td></tr>' % (dfont,len(domlist),string.join(domlist,' ~ '))]
                else: html += ['<tr><td><p>%sDomain Hubs: </FONT><i>None.</i></p></td></tr>' % (dfont)]
                if randlist: html += ['<tr><td><p>%sRandom Datasets (%d): </FONT>%s</p></td></tr>' % (dfont,len(randlist),string.join(randlist,' ~ '))]
                else: html += ['<tr><td><p>%sRandom Datasets: </FONT><i>None.</i></p></td></tr>' % (dfont)]
                html += ['</table>','<!-- ~~~~ End %s SLiM Summary details table ~~~~ -->' % slim,'']

            return string.join(html, '\n')
        except: self.errorLog('specialGeneSummary() error'); return '<p><i>HAPPI code error!</i></p>'
#########################################################################################################################
    def specialDrop(self,drop=[],ttype='results'):  ### Modify dropped columns based on special analysis codes
        '''Modify dropped columns based on special analysis codes.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'hm' in self.list['Special']:
                if ttype == 'results':
                    drop += ['Identifier','HGNC','Entrez','EnsG','Sample']
                    drop += string.split('Rep|AA	Int|AA	PepNorm|AA	Rep|AF	Int|AF	PepNorm|AF	Rep|AL	Int|AL	PepNorm|AL')
                    drop += string.split('Rep|RA	Int|RA	PepNorm|RA	Rep|RF	Int|RF	PepNorm|RF	Rep|RL	Int|RL	PepNorm|RL')
                    drop += string.split('Rep|XL	Int|XL	PepNorm|XL	AA.v.AL	AF.v.AL	RF.v.RL	RA.v.RL	pj.AA.v.AL	pj.AF.v.AL')
                    drop += string.split('pj.RF.v.RL	pj.RA.v.RL	pi.AA.v.AL	pi.AF.v.AL	pi.RF.v.RL	pi.RA.v.RL')
                #if ttype == 'class':
                if ttype == 'xreftab': drop += ['Alias','Species','EntrezCheck','EnsDesc']
                #if ttype == 'mcode':
                if ttype == 'raw': drop += string.split('Identifier      HGNC    EnsG    Desc #')
            return self.list['DropFields'] + drop
        except: self.errorLog('specialDrop() Error')
#########################################################################################################################
    def resultTableHTML(self,headers,data,asdict=True,drop=[],width={},ttype='results'):     ### Write data to HTML table
        '''Write data to HTML table.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = headers[0:]   # Do not want to change original headers list!
            drop = self.specialDrop(drop,ttype)
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
                if h not in width and h in ['Members']: width[h] = 750
                if h in width: hcode += '<th width=%s title="%s">%s</th>\n' % (width[h],self.titleText(ttype,h),h)
                else: hcode += '<th title="%s">%s</th>\n' % (self.titleText(ttype,h),h)
            rows = [hcode]
            #rows = ['<TD><B>%s</B></TD>' % string.join(headers,'</B></TD><TD><B>')]
            ## ~ [2b] ~ Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for dat in datalist:
                rdata = []
                for h in headers:
                    if h in dat:
                        if h in ['Members']:
                            hdat = []
                            for gene in string.split(dat[h],'|'): hdat.append(self.geneLink(gene,frontpage=False))
                            rdata.append(string.join(hdat,'| '))
                        elif h in ['NRID','Hub','HGNC','Spoke','Gene',self.info['GeneField']]:
                            gene = self.gene(str(dat[h]))
                            rdata.append(self.geneLink(gene,frontpage=False,altgene=dat[h]))
                        elif h in ['Density','Score']:
                            try: rdata.append(rje.expectString(dat[h]))
                            except: rdata.append(dat[h])
                        elif h in ['AvF']:
                            try: rdata.append(float(rje.expectString(dat[h])))
                            except: rdata.append(dat[h])
                        elif h in self.list['MTypes']:
                            if h == 'MCODE': ctitle = 'MCODE cluster %s' % dat[h]
                            else: ctitle = '%s MCODE cluster %s' % (h,dat[h])
                            href = '../%s/%s.html' % (h.lower(),dat[h])
                            rdata += ['<a href="%s" title="%s" target="_blank">%s</a>' % (href,ctitle,dat[h])] 
                        elif h in ['AvFnet']:
                            ctitle = 'AvF subnetwork %s' % dat[h]
                            href = '../%s/%s.html' % (h.lower(),dat[h])
                            rdata += ['<a href="%s" title="%s" target="_blank">%s</a>' % (href,ctitle,dat[h])] 
                        elif h == 'XRef_Gene':
                            href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' % dat[h]
                            title = self.titleText('xref','gene')
                            rdata += ['<a href="%s" title="%s">%s</a>' % (href,title,dat[h])] 
                        elif h == 'UniProt':
                            href = 'http://www.uniprot.org/uniprot/%s' % dat[h]
                            title = self.titleText('xref','uniprot')
                            rdata += ['<a href="%s" title="%s">%s</a>' % (href,title,dat[h])] 
                        elif h == 'EnsEMBL':
                            href = 'http://www.ensembl.org/Homo_sapiens/geneview?gene=%s' % dat[h]
                            title = self.titleText('xref','ensembl')
                            rdata += ['<a href="%s" title="%s">%s</a>' % (href,title,dat[h])] 
                        elif h == 'HPRD':
                            href = 'http://www.hprd.org/summary?hprd_id=%s&isoform_id=%s_1&isoform_name=Isoform_1' % (dat[h],dat[h])
                            title = self.titleText('xref','hprd')
                            rdata += ['<a href="%s" title="%s">%s</a>' % (href,title,dat[h])] 
                        elif h == 'OMIM':
                            href = 'http://www.ncbi.nlm.nih.gov/omim/%s' % dat[h]
                            title = self.titleText('xref','omim')                        
                            rdata += ['<a href="%s" title="%s">%s</a>' % (href,title,dat[h])] 
                        else: rdata.append(str(dat[h]))
                    else: rdata.append('')
                    if h in width: rdata[-1] = '<TD WIDTH=%s>%s</TD>' % (width[h],rdata[-1])
                    else: rdata[-1] = '<TD>%s</TD>' % (rdata[-1])
                rows.append(string.join(rdata,''))
        except:
            self.errorLog('resultTableHTML error')
            self.deBug(headers)
            self.deBug(datalist)
            self.deBug(rdata)
        return '<TABLE BORDER=%d><TR VALIGN="top">\n%s\n</TR></TABLE>\n' % (self.stat['Border'],string.join(rows,'\n</TR><TR VALIGN="top">\n'))
#########################################################################################################################
    def interactorTableHTML(self,gene,full=False): ### Table of interactors with evidence & SLiMs returned for 4 datasets types
        '''
        Table of interactors in sample: Spoke, Evidence, Description, Class
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = []
            ppi = self.obj['PPI'].ppi()
            try: ilist = rje.sortUnique(ppi[gene])
            except: ilist = []
            ### ~ [2] Generate table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Table Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = ['','<h2>Protein-Protein Interactions</h2>','',
                    '<table width=1000 border=%d>' % self.stat['Border']]
            html += ['<tr><th align=left title="%s">Gene</th>' % self.titleText('ppi','gene'),
                     '<th align=left title="%s">Evidence</th>' % self.titleText('ppi','evidence'),
                     '<th align=left title="%s">Description</th>' % self.titleText('ppi','description'),
                     '<th align=left title="%s">Class</th>' % self.titleText('results','Class')]
            html.append('</tr>')
            ## ~ [2b] Table contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ix = 0
            for i in ilist[0:]:
                if i not in self.list['Genes'] and not full: continue
                html += ['<tr valign="top"><td width=100>%s</td>' % self.geneLink(i)]
                hkey = '%s\t%s' % (gene,i)
                html += ['<td width=400>%s</td>' % string.join(string.split(ppi[gene][i],'|'),';<br>')]
                html += ['<td width=400>%s</td>' % self.geneDesc(i)]
                html += ['<td width=100>%s</td>' % self.geneClass(i)]
                html += ['</tr>','']
                ix += 1
            if not ix: html += ['<tr><td align=left colspan=4><i>No %s interactors in data.</i></td></tr>' % gene]
        except: self.errorLog('interactorTableHTML Error')
        ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += ['</table>','']
        return string.join(html,'\n')
#########################################################################################################################
    ### <5> ### Gene pages & general HTML code                                                                          #
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
            seqid = '-'; desc = '<i>No sequence mapping.</i>'
            if 'slim' in self.list['Special']: slim = gene; code = rje_slim.slimFromPattern(slim)
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
            xrefhtml = []
            for db in ['Gene','UniProt','EnsEMBL','HPRD','OMIM']:
                if not dbxref or db not in dbxref: continue
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
                xrefhtml += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            if 'humsf' in self.list['Special']:
                if 'slim' in self.list['Special']: href = 'http://www.southampton.ac.uk/~re1u06/research/humsf09/slim/%s.html' % code
                else: href = 'http://www.southampton.ac.uk/~re1u06/research/humsf09/gene/%s.html' % gene
                title = alt = 'Human SLiMFinder %s results' % gene
                src = '../resources/humsf.png'
                xrefhtml += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            if 'slim' in self.list['Special']:
                href = 'http://bioware.ucd.ie/~compass/cgi-bin/formParser.py?name_server=slimsearch2&motif_str=%s' % slim
                title = alt = 'Human proteome SLiMSearch2 with %s' % slim
                src = '../resources/ucd.gif'
                xrefhtml += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            if xrefhtml: html += ['<tr><td colspan=2>'] + xrefhtml + ['</td></tr>']
            if not 'slim' in self.list['Special']: html += ['<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,desc)]
            html += ['</table>','<!-- ~~~~ End %s Gene Summary details table ~~~~ -->' % gene,'']
            ### ~ [3] ~ Summary Data Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db().getTable('Main')
            html += [self.specialGeneSummary(gene)]    #resultTableHTML(main.fields(),main.subset(self.info['GeneField'],gene))]
            ### ~ [4] ~ GO tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gtab = []; godata = {}
            for gentry in self.db('GO').indexEntries('Gene',gene):
                gtype = gentry['GO_Type'].upper()
                if gtype not in godata: godata[gtype] = []
                godata[gtype].append((gentry['GO_ID'],gentry['GO_Desc']))
            for gtype in ['BP','CC','MF']:
                if gtype in godata:
                    gdict = {}
                    for gotup in godata[gtype]:
                        if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                        #gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                        gdict[gotup[1]] = self.goLink(gotup[1],gotup[0])
                    if gdict:
                        ghtml = []
                        for g in rje.sortKeys(gdict): ghtml.append(gdict[g])
                        gotit = '%s Gene Ontology categories' % {'CC':'Cellular Component','BP':'Biological Process','MF':'Molecular Function'}[gtype]
                        gtab.append(('GO_%s' % gtype,string.join(ghtml,' ~ '),gotit))
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s' % gene,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % gene,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % gene,'']
        except: self.errorLog('seqDetailsHTML Error')
        ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #print string.join(html,'\n')
        return string.join(html,'\n')
#########################################################################################################################
    def geneDesc(self,gene):    ### Returns gene description
        '''Returns gene description.'''
        dbxref = self.db().getTable('DBXRef')
        for field in ['Gene','UniProt','EnsEMBL','EnsLoci']:
            if dbxref and gene in dbxref.index(field): return dbxref.indexEntries(field,gene)[0]['EnsDesc']
        return ''
#########################################################################################################################
    def genePage(self,gene):     ### Generates HTML for gene page
        '''
        Generates HTML for gene page.
        >> gene:str = Gene for page construction
        '''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'slim' in self.list['Special']: hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'slim/'),rje_slim.slimFromPattern(gene))
            else: hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
            dbxref = {}; seqid = '-'; desc = '<i>No sequence mapping.</i>'; godata = {}
            if gene in self.data('DBXRef'): 
                dbxref = self.data('DBXRef')[gene]
                seqid = self.data('DBXRef')[gene]['EnsLoci']
                desc = self.data('DBXRef')[gene]['EnsDesc']
                ensg = self.data('DBXRef')[gene]['EnsEMBL']
            html = htmlHead(gene,self.list['StyleSheets'],nobots=self.opt['NoBots']) + self.seqDetailsHTML(gene,dbxref) #gene,seqid,dbxref,desc,godata)
            tablist = []
            ## ~ [1a] ~ Interactors tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('Interactors',self.interactorTableHTML(gene),'Table of returned %s interactors' % gene))
            tablist.append(('Full-PPI',self.interactorTableHTML(gene,full=True),'Full list of %s PPI' % gene))
            ## ~ [1b] ~ UniProt tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uhtml = self.uniProtHTML(gene)
            if uhtml:
                uhtml = uhtml[uhtml.find('<H2>UniProt Format'):]
                uhtml = uhtml[:uhtml.find('</PRE>')] + '</PRE>\n'
                tablist.append(('UniProt',uhtml,'%s UniProt data' % gene))
            ## ~ [1b] ~ Raw Data tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if gene in self.geneList():
                for dfile in self.list['GeneData']:
                    name = string.split(dfile,'.')[-2]
                    rawdb = self.db(name)
                    rawdata = self.specialGeneData(name,gene)
                    if rawdata: tablist.append((name,self.resultTableHTML(rawdb.fields(),rawdata,drop=[self.info['GeneField'],'Gene'],ttype='raw'),'%s data for %s' % (name,gene)))
                    elif rawdata == None: tablist.append((name,'<i>No %s data mapping via "Gene" or "%s"</i>' % (name,self.info['GeneField']),'%s data for %s' % (name,gene)))
                    else: tablist.append((name,'<i>No %s data mapped to %s</i>' % (name,gene),'%s data for %s' % (name,gene)))
            ## ~ [1c] ~ Interactome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try: ilist = rje.sortKeys(self.ppi()[gene])
            except: ilist = []
            if len(ilist) > 2 or (gene not in ilist and len(ilist) > 1):
                if 'slim' in self.list['Special']: intpng = '../interactome_%s/%s.interactome.png' % (self.png(),rje_slim.slimFromPattern(gene))
                else: intpng = '../interactome_%s/%s.interactome.png' % (self.png(),gene)
                svglink = string.replace(intpng,'png','svg')
                alt = '%s interactome image not yet made' % gene
                title = 'Graphical representation of %s interactome (returned proteins only). Click to open in new window. Larger interactomes may not be clear. See tables for details.' % gene
                if self.opt['SVG']:
                    img = self.svgCode(svglink,title)
                    title = 'Graphical representation of %s interactome (returned proteins only). Larger interactomes may not be clear. See tables for details.' % gene
                else: img = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)
                if self.opt['XGMML']:
                    x = '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                    tablist.append(('Interactome','<h2>%s Interactome</h2>\n%s%s\n' % (gene,x,img),title))
                else: tablist.append(('Interactome','<h2>%s Interactome</h2>\n%s\n' % (gene,img),title))
            ## ~ [1d] ~ Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for mtype in self.list['MTypes']:
                cdb = self.obj['PPI'].db(mtype)
                mtitle = '%s clusters containing %s\n' % (mtype,gene)
                chtml = '<h2>%s clusters containing %s</h2>\n' % (mtype,gene)
                for complex in cdb.index(mtype):
                    seed = cdb.data()[cdb.index(mtype)[complex][0]]['Seed']
                    if gene not in self.obj['PPI'].dict[mtype][seed]: continue
                    ctitle = 'Predicted %s cluster %s (Seed = %s)' % (mtype,complex,seed)
                    chtml += '<h3>%s</h3>\n' % ctitle
                    chtml += self.resultTableHTML(cdb.fields(),[cdb.data()[seed]],asdict=False,ttype='mcode')
                    intpng = '../%s_%s/%s.%s.png' % (mtype.lower(),self.png(),complex,mtype.lower())
                    if self.opt['XGMML']: chtml += '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                    svglink = '../%s_svg/%s.%s.svg' % (mtype.lower(),complex,mtype.lower())
                    alt = '%s interactome image not yet made' % ctitle
                    title = 'Graphical representation of %s. Click to open in new window. Larger interactomes may not be clear.' % ctitle
                    if self.opt['SVG']:
                        img = self.svgCode(svglink,title)
                        title = 'Graphical representation of %s. Larger interactomes may not be clear. See tables for details.' % ctitle
                    else: img = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)
                    chtml += '%s\n' % img
                tablist.append((mtype,chtml,mtitle))
            ## ~ [1e] ~ Special ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'hm' in self.list['Special']:
                cdb = self.obj['PPI'].db('AvFnet')
                for complex in cdb.index('AvFnet'):
                    seed = cdb.data()[cdb.index('AvFnet')[complex][0]]['Seed']
                    if gene not in self.obj['PPI'].dict['AvFnet'][seed]: continue
                    ctitle = 'AvF subnetwork %s (Seed = %s)' % (complex,seed)
                    chtml = '<H2>%s</H2>\n' % ctitle
                    chtml += self.resultTableHTML(cdb.fields(),[cdb.data()[seed]],asdict=False,ttype='mcode')
                    #!#if self.opt['SVG']
                    intpng = '../avfnet_%s/%s.avfnet.png' % (self.png(),complex)
                    if self.opt['XGMML']: chtml += '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                    alt = '%s interactome image not yet made' % ctitle
                    title = 'Graphical representation of %s. Click to open in new window. Larger interactomes may not be clear.' % ctitle
                    tablist.append((complex,'%s\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (chtml,intpng,intpng,alt,title),title))
            ### ~ [2] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html += tabberHTML(gene,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('genePage Error')
#########################################################################################################################
    def uniProtHTML(self,gene):  ### Generate HTML-linked DAT file            
        '''Generate HTML-linked DAT file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xref = xdb = self.db('DBXRef')
            try: acc = xdb.data(gene)['UniProt']
            except:
                if self.opt['DeBug']: self.errorLog('XRef "%s"' % gene,quitchoice=False)
                return ''
            hdir = rje.makePath('%sunifake/' % self.info['HTMLPath'])
            rje.mkDir(self,hdir)
            ### ~ [2] ~ Generate DAT File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Generate general attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hfile = '%s%s.html' % (hdir,gene)
            if not self.opt['Force'] and rje_html.checkHTML(hfile): return open(hfile,'r').read()
            ## ~ [2b] ~ Compile UniProt data and add Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                shhh = self.log.opt['Silent']
                self.log.opt['Silent'] = True
                uni = rje_uniprot.UniProt(self.log,self.cmd_list)
                uni.readUniProt(clear=True,acclist=[acc],cleardata=False)
                self.log.opt['Silent'] = shhh
                if not uni.list['Entry']: return ''
                dat = uni.list['Entry'][0]
                sequence = dat.obj['Sequence'].info['Sequence'][0:]
            except:
                self.errorLog('Problem making HTML-linked DAT file for %s' % acc)
                return ''
            ### ~ [3] ~ Save data to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = htmlHead(gene,self.list['StyleSheets'])
            ## ~ [3a] ~ Gene summary as with hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html += self.seqDetailsHTML(gene,xref.data()[gene]) #gene,seqid,dbxref,desc,godata)
            ## ~ [3b] ~ UniFake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uhtml = ['<H2>UniProt Format</H2>\n<PRE>']
            entry = dat
            eseq = entry.obj['Sequence']
            ## Standard info ##
            for key in ['ID','AC','DT','DE','GN','OS']:
                if entry.dict['Data'].has_key(key):
                    for rest in entry.dict['Data'][key]:
                        if key == 'ID':
                            id = string.split(rest)
                            rest = string.join(['<a href="http://www.uniprot.org/uniprot/%s" title="Link to UniProt entry">%s</a>' % (id[0],id[0])] + id[1:])
                        #if key == 'GN':
                        #    id = rje.matchExp('^(GN\s+Name=)(\S+)(;.*)',rest)
                        #    if id: rest = string.join([id[0],self.geneLink(id[1]),id[2]],'')
                        uhtml.append('%s   %s\n' % (key,rje.chomp(rest)))
            ## Other data, except Features and sequence ##
            for key in rje.sortKeys(entry.dict['Data']):
                if key not in ['ID','AC','DT','DE','GN','OS','FT','SQ','SEQ','//']:
                    for rest in entry.dict['Data'][key]: uhtml.append('%s   %s\n' % (key,rje.chomp(rest)))
            ## Features ##
            entry.orderFT()
            for ftdict in entry.list['Feature']:
                (p1,p2) = (ftdict['Start'],ftdict['End'])
                ftxt = 'FT   %s' % ftdict['Type']
                while len(ftxt) < 14 or ftxt[-1] != ' ': ftxt += ' '
                ftxt += '%6s' % ('%d' % p1)
                while len(ftxt) > 20 and ftxt[-(len('%d' % p1)+2):-len('%d' % p1)] == '  ': ftxt = ftxt[:-(len('%d' % p1)+1)] + ftxt[-len('%d' % p1):]
                ftxt += '%7s' % ('%d' % p2)
                while len(ftxt) > 27 and ftxt[-(len('%d' % p2)+2):-len('%d' % p2)] == '  ': ftxt = ftxt[:-(len('%d' % p2)+1)] + ftxt[-len('%d' % p2):]
                if 'PosLink' in ftdict:
                    fsplit = string.split(ftxt)
                    for i in range(1,len(fsplit)): fsplit[i] = '%s%s</a>' % (ftdict['PosLink'],ftsplit[i])
                    ftxt = string.join(fsplit)
                ftxt += ' %s\n' % ftdict['Desc']
                uhtml.append(ftxt)
            ## Sequence/End ##
            uhtml.append('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % eseq.aaLen())),eseq.aaLen(),rje_sequence.MWt(eseq.info['Sequence'])))
            uniseq = eseq.info['Sequence'][0:]
            while len(uniseq) > 0:
                uhtml.append('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                uniseq = uniseq[60:]
            uhtml.append('//\n</PRE>\n')
            html += string.join(uhtml,'') + htmlTail()
            open(hfile,'w').write(html)
            return html
        except: self.errorLog('UniFake problem')
#########################################################################################################################
    def goPage(self,go,goid):   ### Make GO Pages
        '''Make GO Pages.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not go and goid: return False
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
            godata = self.obj['GO'].go(goid)
            ### ~ [1] ~ Make Page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html = [htmlHead(go,self.list['StyleSheets'],nobots=self.opt['NoBots']),
                    '','<!-- ~~~~ %s GO Summary details table ~~~~ -->' % goid,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="40%%"><a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s<b>GO:%s</b></FONT></a></td>' % (goid,gfont,goid),
                    '<td width="40%%">%s%s</FONT></td>' % (ifont,godata['type']),
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,go)]
            html += ['</table>','<!-- ~~~~ End %s GO Summary details table ~~~~ -->' % goid,'']
            ### ~ [2] Info Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = []
            ## ~ [2a] ~ GO parents  & children List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            for rtype in ['is_a','part_of']:
                try:
                    for rel in godata[rtype]: tabhtml.append('<li>%s %s %s</li>' % (go,rtype,self.goLink(self.obj['GO'].go(rel)['name'],rel)))                       
                except: pass
            if tabhtml: tabhtml.append('<br>')
            for id in godata['child_terms']:
                for rtype in ['is_a','part_of']:
                    try:
                        for rel in self.obj['GO'].go(id)[rtype]:
                            if rel == goid: tabhtml.append('<li>%s %s %s</li>' % (self.goLink(self.obj['GO'].go(id)['name'],id),rtype,self.obj['GO'].go(rel)['name']))                       
                    except: pass
            if tabhtml:
                tabhtml = ['<H2>GO Relationships</H2>','<ul>'] + tabhtml + ['</ul>']
                tablist.append(('GO',string.join(tabhtml,'\n'),'GO Term relationships for %s' % go))
            ## ~ [2b] ~ Gene List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []; jtxt = ' ~ '
            genes = self.dict['GOPages'][goid]
            genes.sort()
            for gene in genes:
                gclass = self.geneClass(gene,'')
                if gclass: gtxt = '%s (%s)' % (self.geneLink(gene),gclass)
                else: gtxt = gene
                tabhtml.append(gtxt)
            if tabhtml: tabhtml[0] = '<h2>%s Genes</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No genes in analysis mapped to %s</i>' % go]
            tablist.append(('Genes',string.join(tabhtml,jtxt),'Genes mapped to %s' % go))
            if len(genes) <= self.stat['MaxGO']:
                title = 'Input data for %s genes' % go
                datlist = []; main = self.db('Main')
                for gene in genes: datlist.append(main.data()[gene])
                tablist.append(('Data',self.resultTableHTML(main.fields(),datlist,asdict=False),title))
            #!# Consider adding Complexes? Parents? Children? #!#
            ## ~ [2c] ~ Interactome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            intpng = '../go_%s/%s.go.png' % (self.png(),goid)
            alt = '%s interactome image not made' % goid
            title = 'Graphical representation of %s interactome (returned proteins only). Click to open in new window. Larger interactomes may not be clear. See tables for details.' % go            
            if self.opt['SVG']:
                svglink = '../go_svg/%s.go.svg' % goid
                img = self.svgCode(svglink,title)
                title = 'Graphical representation of %s interactome (returned proteins only). Larger interactomes may not be clear. See tables for details.' % go
            else: img = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)
            if self.opt['XGMML']: x = '\n<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>' % intpng[:-4]
            else: x = ''
            tablist.append(('Interactome','<h2>%s Interactome</h2>%s\n%s</a>\n' % (go,x,img),title))
            ## ~ [2d] ~ Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for mtype in self.list['MTypes']:
                cdb = self.obj['PPI'].db(mtype)
                mtitle = '%s clusters containing %d+ %s genes' % (mtype,self.stat['GOPages'],go)
                chtml = '<h2>%s</h2>\n' % mtitle
                for complex in cdb.index(mtype):
                    seed = cdb.data()[cdb.index(mtype)[complex][0]]['Seed']
                    gox = 0
                    for gene in self.obj['PPI'].dict[mtype][seed]:
                        if gene in self.dict['GOPages'][goid]: gox += 1
                    if gox < self.stat['GOPages']: continue
                    ctitle = 'Predicted %s cluster %s (Seed = %s); %d x %s genes.' % (mtype,complex,seed,gox,go)
                    chtml += '<h3>%s</h3>\n' % ctitle
                    chtml += self.resultTableHTML(cdb.fields(),[cdb.data()[seed]],asdict=False,ttype='mcode')
                    intpng = '../%s_%s/%s.%s.png' % (mtype.lower(),self.png(),complex,mtype.lower())
                    if self.opt['XGMML']: chtml += '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                    svglink = '../%s_svg/%s.%s.svg' % (mtype.lower(),complex,mtype.lower())
                    alt = '%s interactome image not yet made' % ctitle
                    title = 'Graphical representation of %s. Click to open in new window. Larger interactomes may not be clear.' % ctitle
                    if self.opt['SVG']:
                        img = self.svgCode(svglink,title)
                        title = 'Graphical representation of %s. Larger interactomes may not be clear.' % gene
                    else: img = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)
                    #tablist.append(('%s:%s' % (mtype,complex),'%s\n%s\n' % (chtml,img),title))
                    chtml += '%s\n' % img
                tablist.append((mtype,chtml,mtitle))
            ## ~ [2e] ~ Special ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'hm' in self.list['Special']:
                cdb = self.obj['PPI'].db('AvFnet')
                for complex in cdb.index('AvFnet'):
                    seed = cdb.data()[cdb.index('AvFnet')[complex][0]]['Seed']
                    gox = 0
                    for gene in self.obj['PPI'].dict['AvFnet'][seed]:
                        if gene in self.dict['GOPages'][goid]: gox += 1
                    if gox < self.stat['GOPages']: continue
                    ctitle = 'AvF subnetwork %s (Seed = %s); %d x %s genes.' % (complex,seed,gox,go)
                    chtml = '<H2>%s</H2>\n' % ctitle
                    chtml += self.resultTableHTML(cdb.fields(),[cdb.data()[seed]],asdict=False,ttype='mcode')
                    #!#if self.opt['SVG']
                    intpng = '../avfnet_%s/%s.avfnet.png' % (self.png(),complex)
                    if self.opt['XGMML']: chtml += '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n' % intpng[:-4]
                    alt = '%s interactome image not yet made' % ctitle
                    title = 'Graphical representation of %s. Click to open in new window. Larger interactomes may not be clear.' % ctitle
                    tablist.append((complex,'%s\n<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (chtml,intpng,intpng,alt,title),title))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = string.join(html,'\n') + tabberHTML(goid,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('goPage(%s:%s) Error' % (goid,go))
#########################################################################################################################
    def mcodePage(self,id,mtype):   ### Generates HTML for MCODE and Complex pages
        '''Generates HTML for MCODE and Complex pages.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pobj = self.obj['PPI']
            #self.deBug(pobj)
            cdb = pobj.db(mtype)
            #self.deBug(cdb)
            seed = cdb.index(mtype)[id][0]
            #self.deBug(seed)
            cdata = cdb.data()[seed]
            #self.deBug(cdata)
            cgenes = string.split(cdata['Members'],'|')
            html = [htmlHead('%s_%s' % (mtype.upper(),id),self.list['StyleSheets'],nobots=self.opt['NoBots'])]
            tablist = []
            mname = '%s %s' % (mtype,id)
            ## ~ [1a] ~ Page header & summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html += ['','<!-- ~~~~ %s Cluster Summary details table ~~~~ -->' % id,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="30%%">%s<b>%s %s</b></FONT></a></td>' % (gfont,mtype,id),
                    '<td width="50%%">%s<b>(Seed = %s)</b></FONT></td>' % (ifont,seed),
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>']
            if mtype == 'MCODE': html.append('<tr><td colspan=2>%s<p>MCODE predicted cluster (%d genes)</p></FONT></td></tr>' % (dfont,len(cgenes)))
            elif mtype == 'AvFnet': html.append('<tr><td colspan=2>%s<p>AvF subnetwork (%d genes)</p></FONT></td></tr>' % (dfont,len(cgenes)))
            else: html.append('<tr><td colspan=2>%s<p>%s MCODE predicted cluster (%d genes)</p></FONT></td></tr>' % (dfont,mtype,len(cgenes)))
            gosig_insert = len(html)
            html += ['</table>',self.resultTableHTML(cdb.fields(),[cdata],asdict=False,ttype='mcode')]
           ## ~ [1b] ~ GO Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            go = self.obj['GO']
            gosig = {}; gosiglist = []
            gtab = []; godata = {}; gox = {}; gogenes = {}
            for gene in cgenes:
                for gentry in self.db('GO').indexEntries('Gene',gene):
                    gtype = gentry['GO_Type'].upper()
                    if gtype not in godata: godata[gtype] = []
                    if (gentry['GO_ID'],gentry['GO_Desc']) not in godata[gtype]: godata[gtype].append((gentry['GO_ID'],gentry['GO_Desc']))
                    if gentry['GO_ID'] not in gox: gox[gentry['GO_ID']] = 0
                    if gentry['GO_ID'] not in gogenes: gogenes[gentry['GO_ID']] = []
                    gogenes[gentry['GO_ID']].append(gene)
                    gox[gentry['GO_ID']] += 1
            for gtype in ['BP','CC','MF']:
                gosig[gtype] = {}
                gobonf = {}
                if gtype in godata:
                    gdict = {}
                    for gotup in godata[gtype]:
                        if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                        if gox[gotup[0]] < self.stat['GOPages']: continue
                        #gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                        try:
                            k = gox[gotup[0]]; N = len(cgenes); p = len(self.dict['GOPages'][gotup[0]]) / float(len(self.list['Genes']))
                            sig = rje.logBinomial(k,N,p,callobj=self)
                        except: sig = 1.0
                        if sig <= 0.05:
                            if sig not in gosig[gtype]: gosig[gtype][sig] = []
                            gosig[gtype][sig].append({'GO_ID':gotup[0],'GO_Desc':gotup[1],'GO_Type':gtype.lower(),'N':k,'p':rje.expectString(sig),'bonf':rje.expectString(sig*len(self.dict['GOPages']))})
                            gdict[gotup[1]] = '%s (%d; p=%s)' % (self.goLink(gotup[1],gotup[0]),gox[gotup[0]],rje.expectString(sig))
                            for gclass in self.dict['ClassGenes']:
                                gosig[gtype][sig][-1][gclass] = rje.listIntersect(self.dict['ClassGenes'][gclass],gogenes[gotup[0]])
                                gosig[gtype][sig][-1][gclass].sort()
                                gtxt = []
                                for gene in gosig[gtype][sig][-1][gclass]: gtxt.append(self.geneLink(gene,False))
                                gosig[gtype][sig][-1][gclass] = string.join(gtxt,'; ')
                            gobonf[gotup[0]] = sig
                        else: gdict[gotup[1]] = '%s (%d)' % (self.goLink(gotup[1],gotup[0]),gox[gotup[0]])
                    if gdict:
                        ghtml = []
                        for g in rje.sortKeys(gdict): ghtml.append(gdict[g])
                        gtab.append(('GO_%s' % gtype,string.join(ghtml,' ~ '),self.titleText('tab','go')))
                ## ~ Reduce if more sig parents or single children ~ ##
                subsig = []
                for goid in rje.sortKeys(gobonf):
                    for parent in go.parents(goid,all_levels=False):
                        if parent in gobonf and gobonf[parent] < gobonf[goid]: subsig.append(goid); break
                #if subsig: self.deBug('%s -> Chuck %s' % (gobonf,subsig))
                for goid in subsig: gobonf.pop(goid)
                for goid in rje.sortKeys(gobonf):
                    csigx = 0; mostsig = False
                    for child in go.children(goid):
                        if child in gobonf:
                            csigx += 1
                            if gobonf[child] > gobonf[goid]: mostsig = True
                    if 0 < csigx < 2 and not mostsig: subsig.append(goid)
                #if subsig: self.deBug('%s -> Chuck %s' % (gobonf,subsig))
                ## ~ Add HTML ~ ##
                gosig_added = False
                for sig in rje.sortKeys(gosig[gtype]):
                    for entry in gosig[gtype][sig]:
                        if entry['GO_ID'] not in subsig: gosiglist.append(entry)
                    if not gosig_added:
                        gsdata = gosig[gtype][sig][0]
                        gosig_sum = '%s: %s (%d; p=%s)' % (gtype,self.goLink(gsdata['GO_Desc'],gsdata['GO_ID']),gsdata['N'],gsdata['p'])
                        html.insert(gosig_insert,'<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,gosig_sum))
                        gosig_insert += 1; gosig_added = True
                #self.deBug('%s' % gosiglist)
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s' % mname,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % mname,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % mname,'']
            ### ~ [2] ~ Gene Data and PNG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db('Main'); gentries = []
            for gene in cgenes:
                try: gentries.append(main.data()[gene])
                except: pass
            tablist.append(('Data',self.resultTableHTML(main.fields(),gentries,asdict=False),'%s %s Summary' % (mtype,id)))
            if gosiglist: tablist.append(('GO_Sig',self.resultTableHTML(['GO_Type','GO_Desc','GO_ID','N','p','bonf']+self.classOrder(),gosiglist,asdict=False,ttype='gosig'),'GO Enrichment (versus all input genes)'))
            else: tablist.append(('GO_Sig','<i>No enriched GO terms</i>','GO Enrichment (versus all input genes)'))
            intpng = '../%s_%s/%s.%s.png' % (mtype.lower(),self.png(),id,mtype.lower())
            ctitle = '%s cluster %s (Seed = %s)' % (mtype,id,seed)
            alt = '%s interactome image not yet made' % ctitle
            title = 'Graphical representation of %s. Click to open in new window. Larger interactomes may not be clear.' % ctitle
            if self.opt['SVG']:
                svglink = '../%s_svg/%s.%s.svg' % (mtype.lower(),id,mtype.lower())
                ctitle = '%s cluster %s (Seed = %s)' % (mtype,id,seed)
                title = 'Graphical representation of %s. Larger interactomes may not be clear.' % ctitle
                img = self.svgCode(svglink,title)
            else: img = '<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>\n' % (intpng,intpng,alt,title)
            if self.opt['XGMML']: img = '<p>Download XGMML <a href="%s.xgmml.gz">here</a>.</p>\n%s' % (intpng[:-4],img)
            tablist.append(('Interactome',img,title))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + '%s/' % mtype.lower()),id)
            html = string.join(html,'\n') + tabberHTML(mname,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('%s %s page error' % (mtype,id))
#########################################################################################################################
    ### <6> ### PNG Methods                                                                                             #
#########################################################################################################################
    def classCol(self,c,strict=True):   ### Returns colour of Class for PNG
        '''Returns definition of Class.'''
        try: return self.dict['ClassCol'][c]
        except:
            if strict: self.errorLog('ClassDef error - "%s"' % c)
            if 'Default' in self.dict['ClassCol']: return self.dict['ClassCol']['Default']
            return 35
#########################################################################################################################
    def forkPNG(self,forkpng):  ### Fork out PNG generation (lists of (id,type) tuples)
        '''Fork out PNG generation (lists of (id,type) tuples).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not forkpng: return
            forkx = self.stat['Forks']      # Number of forks to have running at one time
            forks = {}      # Dictionary of active forks {pid:basefile}
            killforks = self.stat['KillForks']      # Time in seconds to wait after main thread has apparently finished
            killtime = time.time()
            pngx = len(forkpng)
            self.obj['PPI'].opt['ProgLog'] = False
            ### ~ [1] ~ Cycle through forkpng tuples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while forks or forkpng:
                ## ~ [1a] ~ Add more forks if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while forkpng and (len(forks) < forkx):     # Add more forks
                    (id,mtype) = forkpng.pop(0)
                    if mtype == 'interactome' and 'slim' in self.list['Special']: basefile = 'tmp_%s-%s' % (mtype,rje_slim.slimFromPattern(id))
                    else: basefile = 'tmp_%s-%s' % (mtype,id)
                    forkcmd = self.cmd_list + ['i=-1','log=%s.log' % basefile,'errorlog=None']
                    newpid = os.fork() 
                    if newpid == 0: # child
                        self.opt['Child'] = True
                        forkobj = HAPPI(cmd_list=forkcmd)
                        forkobj.log.info['LogFile'] = '%s.log' % basefile
                        forkobj.obj = self.obj
                        forkobj.info = self.info
                        forkobj.stat = self.stat
                        forkobj.opt = self.opt
                        forkobj.list = self.list
                        forkobj.dict = self.dict
                        forkobj.interactomePNG(id,mtype)
                        os._exit(0)    # Exit process 
                    elif newpid == -1: self.errorLog('Problem forking %s.' % id,printerror=False)  
                    else: forks[newpid] = basefile            
                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(1)       # Sleep for 1s 
                forklist = self._activeForks(forks.keys())
                if len(forklist) != len(forks):
                    self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    for pid in rje.sortKeys(forks):    # Go through current forks
                        if pid not in forklist:
                            killtime = time.time()  # Reset killtime - still doing stuff
                            try:
                                fork = forks.pop(pid)
                                self.printLog('#FORK','HAPPI Fork %s Finished! Transfering log details to %s.' % (fork,self.log.info['LogFile']),log=False)
                                rje.fileTransfer(fromfile='%s.log' % fork,tofile=self.log.info['LogFile'],deletefrom=True,append=True)
                            except: self.errorLog('Fork Error')
                            #self.deBug('%s.log? (%s)' % (basefile,not os.path.exists('%s.log' % basefile)))
                ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if time.time() - killtime > killforks:
                    self.printLog('#KILL','%d seconds of main thread inactivity. %d forks still active!' % (killforks,len(forks)))
                    for fork in forks: self.printLog('#FORK','Fork %s, PID %d still Active!' % (forks[fork],fork))
                    if self.stat['Interactive'] < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                    elif rje.yesNo('Kill hanging forks?'):
                        for fork in rje.sortKeys(forks):
                            self.printLog('#KILL','Killing Fork %s, PID %d.' % (forks[fork],fork))
                            os.system('kill %d' % fork)
                    else: killtime = time.time()
            self.printLog('#FORK','End of %s interactomes image forking.' % (rje.iStr(pngx)))
            self.obj['PPI'].opt['ProgLog'] = self.opt['ProgLog']
        except SystemExit: os._exit(0)
        except: self.errorLog('%s.forkPNG error' % self)
#########################################################################################################################
    def interactomePNG(self,id,type='interactome'):   ### Peforms the R visualisation for a given Gene's PPI Network
        '''Peforms the R visualisation for a given Gene's PPI Network.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SVG']: pngdir = '%s%s_svg/' % (self.info['HTMLPath'],type); rje.mkDir(self,pngdir)
            else: pngdir = '%s%s_png/' % (self.info['HTMLPath'],type); rje.mkDir(self,pngdir)
            if type == 'interactome' and 'slim' in self.list['Special']: basefile = rje.makePath('%s%s.%s' % (pngdir,rje_slim.slimFromPattern(id),type),wholepath=True)
            else: basefile = rje.makePath('%s%s.%s' % (pngdir,id,type),wholepath=True)
            if self.opt['SVG'] and os.path.exists('%s.svg' % basefile) and not self.opt['Force']: return False
            elif os.path.exists('%s.png' % basefile) and not self.opt['Force']: return False
            pobj = self.obj['PPI']
            main = self.db().getTable('Main')
            dis = self.obj['Dis']
            mtypes = string.split(string.join(self.list['MTypes']).lower())
            ## ~ [1a] ~ Type-specific code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if type == 'interactome':
                gene = id
                ilist = pobj.genePPI(gene)
                #!# Add optional full/partial? (Remove genes not in input data) #!#
                if gene in ilist: ilist.remove(gene)
            elif type == 'geneclass':
                if id == 'combined':
                    ilist = []
                    for gclass in self.list['Combine']: ilist += self.dict['ClassGenes'][gclass]
                else: ilist = self.dict['ClassGenes'][id]
            elif type in mtypes:
                i = mtypes.index(type)
                mtype = self.list['MTypes'][i]
                cdb = pobj.db(mtype)
                seed = cdb.data()[cdb.index(mtype)[id][0]]['Seed']
                ilist = pobj.dict[mtype][seed]                
            elif type == 'avfnet':
                cdb = pobj.db('AvFnet')
                seed = cdb.data()[cdb.index('AvFnet')[id][0]]['Seed']
                ilist = pobj.dict['AvFnet'][seed]
            elif type == 'go':
                ilist = self.dict['GOPages'][id]
            ## ~ [1b] ~ Universal code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not ilist: return False
            pobj = rje_ppi.PPI(self.log,self.cmd_list)
            pobj.obj = self.obj['PPI'].obj
            pobj.dict = self.obj['PPI'].dict
            pobj.opt['ProgLog'] = self.opt['NoForks']
            if len(ilist) > self.stat['PNGMax']:
                self.printLog('#MAX','Too many genes for %s.png (%s genes)' % (basefile,rje.iStr(len(ilist))))
                G = {'%s %s Genes' % (rje.iStr(len(ilist)),id):[]}
            else:
                if self.opt['SVG']: self.printLog('#%sSVG' % type.upper()[0],'Making %s.svg (%s genes)' % (basefile,rje.iStr(len(ilist))))
                else: self.printLog('#%sPNG' % type.upper()[0],'Making %s.png (%s genes)' % (basefile,rje.iStr(len(ilist))))
                #self.deBug('%s: %s' % (id,ilist))
                G = rje_ppi.subGraph(pobj.ppi(),ilist)      # Reduced PPI Graph
                G = rje_ppi.removeOrphans(G,self)
            if not G: G = {'No PPI':[]}
            nfile = '%s.npos.tdt' % basefile
            npos = {}
            if self.opt['UsePos'] and os.path.exists(nfile): npos = pobj.loadNPos(nfile)
            if rje.sortKeys(npos) != rje.sortKeys(G): npos = pobj.rjeSpringLayout(G); pobj.saveNPos(npos,nfile)
            elif self.opt['UpdatePos']:
                npos = rje_ppi.rjeSpringLayout(G,damping=pobj.stat['Damping'],callobj=pobj,walltime=pobj.stat['Walltime'],nudge=pobj.stat['NudgeCyc'],prepos=npos)
                if os.path.exists(nfile): os.unlink(nfile)
                pobj.saveNPos(npos,nfile)
            classcol = self.dict['ClassCol']   #{'AA':15,'FF':19,'AA-FF':8,'LYS':6,'ITGA2B':3,'ITGAL':3,'ITGA5':3,'Integrin':3}
            #!# Add Class Colour from DB Table #!#
            if pobj.opt['ColByDeg']: pobj.addCol(default=12,coldict=classcol,G=G,ckey=self.info['GeneClass'],edb=main,addcol=48)
            if self.opt['XGMML']:
                xgmml = pobj.ppiToXGMML(G,id)
                xgmml.dict['NodePos'] = npos
                xgmml.saveXGMML('%s.xgmml' % basefile)
                self.gZip('%s.xgmml' % basefile)
            if not dis:
                spokenum = len(npos)
                if 'slim' in self.list['Special']: font = 20 - int(min(80,spokenum)/10.0)
                if self.opt['SVG']: pobj.saveSVG(npos,basefile,G,width=1600,ntype='ellipse',backups=False)
                else: pobj.saveR(npos,basefile,G,cleantdt=True,backups=False)
                return True
            ### ~ [2] ~ Load distance matrix and make Tree/Heatmap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            nsftree = dis.upgma(objkeys=ilist,checksym=False)
            open('%s.nsf' % basefile,'w').write(nsftree)
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            upgma.rTree('%s.tree.csv' % basefile,seqname='short')
            reorder = upgma._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            dis.info['Name'] = '%s interactome GABLAM' % id
            dis.saveMatrix(reorder,basefile+'.heatmap.tdt',delimit='\t')
            pobj.saveTDT(npos,basefile,G,backups=False)
            ### ~ [3] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "hminteractome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/hmhtml.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
        return True
#########################################################################################################################
    def keyPNG(self):   ### Peforms the R visualisation for the class colour key
        '''Peforms the R visualisation for the class colour key.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['MakePNG']: return
            pngdir = self.info['HTMLPath']
            basefile = pngdir + 'class_key'
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            pobj = rje_ppi.PPI(self.log,self.cmd_list)
            if self.db('ClassCol'):
                ilist = rje.sortKeys(self.db('ClassCol').index('Class'))
                pobj.db().addEmptyTable('Node',self.db('ClassCol').fields(),keys=['Class']) ### Adds empty table
                pobj.db('Node').dict['Data'] = self.db('ClassCol').data()
            else:
                ilist = rje.sortKeys(self.dict['ClassCol'])
                pobj.db().addEmptyTable('Node',['Class','Col'],keys=['Class']) ### Adds empty table
                for gclass in self.dict['ClassCol']:
                    pobj.db('Node').dict['Data'][gclass] = {'Class':gclass,'Col':self.dict['ClassCol'][gclass]}
            ## ~ [1a] ~ Positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not ilist: return
            clist = self.classOrder()[0:]
            for c in clist:
                if c not in ilist: clist.remove(c)
            for c in ilist:
                if c not in clist: clist.append(c)
            ilist = clist
            npos = {}; G ={}
            for c in ilist:
                npos[c] = rje_ppi.getSpokeXY(ilist.index(c)+1,len(ilist),start_ang=0,scale=1.0); G[c] = []
                npos[c][0] += 1.0; npos[c][1] += 1.0
            #self.deBug(G); self.deBug(npos)
            if self.opt['SVG']: pobj.saveSVG(npos,basefile,G,width=1000,ntype='ellipse',cutspace=False,backups=False)
            else: pobj.saveR(npos,basefile,G,cleantdt=not self.opt['DeBug'],backups=False)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def svgCode(self,svglink,title):    ### Returns HTML Code for SVG file embedding.
        '''Returns HTML Code for SVG file embedding.'''
        try:
            svgfile = self.info['HTMLPath'] + string.join(string.split(svglink,'/')[1:],'/')
            svghtm = '<p title="%s">\n' % (title)
            try: (width,height) = rje.matchExp('<svg width="(\d+)" height="(\d+)" version="1.1"',open(svgfile,'r').read())
            except: (width,height) = (1600,1600)
            svghtm += '<embed src="%s" width="%s" height="%s" type="image/svg+xml"' % (svglink,width,height)
            svghtm += ' pluginspage="http://www.adobe.com/svg/viewer/install/" /></p>'
            return svghtm
        except: self.errorLog(rje_zen.Zen().wisdom()); return '<i>SVG code error!</i>'
#########################################################################################################################
    ### <7> ### "Special" Analysis Methods                                                                              #
#########################################################################################################################
    def specialInData(self):    ### Special run data manipulations
        '''Special run data manipulations.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db('Main')
            ### ~ [1] ~ H&M analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'hm' in self.list['Special']:
                for entry in main.entries():
                    if entry[self.info['GeneClass']] == 'LYS':
                        if float(entry['enrAA']) > 1.0 and float(entry['enrFF']) > 1.0: entry['Class'] = 'LYS-AA-FF'
                        elif float(entry['enrAA']) > 1.0: entry['Class'] = 'LYS-AA'
                        elif float(entry['enrFF']) > 1.0: entry['Class'] = 'LYS-FF'
        except: self.errorLog('Problem during %s specialInData()' % self)
#########################################################################################################################
    def hmComplexes(self):  ### Calculates Peptide enrichment for complexes
        '''Calculates Peptide enrichment for complexes.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            main = self.db('Main'); rnum = 1000
            ppi = self.obj['PPI']

            ### ~ [1] Add AvF Subnetworks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '%s.avfnet.tdt' % self.info['Basefile']
            reformat = {}
            for f in string.split('Nodes Edges Rank'): reformat[f] = 'int'
            for f in string.split('Density Score AvF'): reformat[f] = 'num'
            if os.path.exists(mfile) and self.opt['Force'] and self.yesNo('Remake %s?' % mfile): rje.backup(self,mfile,unlink=True,appendable=False)
            if os.path.exists(mfile):
                cdb = ppi.db().addTable(mfile,mainkeys=['Seed'],datakeys='All',name='AvFnet')
                cdb.dataFormat(reformat)
                ppi.dict['AvFnet'] = {}
                for entry in cdb.entries():
                    try: ppi.dict['AvFnet'][entry['Seed']] = string.split(entry['Members'],'|')
                    except: ppi.dict['AvFnet'][entry['Seed']] = string.split(entry['AvFnet'],'|')
            else:
                G = rje_ppi.subGraph(ppi.ppi(),self.geneList())
                if self.stat['PPExpand']: G = ppi.expandPPI(G,self.stat['PPExpand'],ppcomplex=self.opt['PPComplex'])
                ndb = ppi.db('Node')
                if 'AvF' not in ndb.fields(): ndb.addField('AvF')
                for gene in rje.sortKeys(G):
                    if gene not in ndb.data(): G.pop(gene); continue
                    try: ndb.data()[gene]['AvF'] = float(main.data()[gene]['AvF'])
                    except: ndb.data()[gene]['AvF'] = 0.0
                ppi.valueClusters(ppi=G,vfield='AvF',addlinks=False,save='avfnet')
                cdb = ppi.db('AvFnet')
            if 'AvFnet' not in cdb.fields():
                cdb.list['Fields'].insert(0,'AvFnet')
                cdb.list['Fields'] += self.classOrder()
                for entry in cdb.entries():
                    entry['AvFnet'] = '%s_N%d' % (rje.preZero(entry['Rank'],cdb.entryNum()),entry['Nodes'])
                    for gclass in self.classOrder(): entry[gclass] = 0
                    for gene in ppi.dict['AvFnet'][entry['Seed']]:
                        for gclass in self.geneClass(gene,split=True):
                            try: entry[gclass] += 1
                            except: pass
            # ['Seed','Nodes','Edges','Density','Score','Complex']

            ### ~ [2] Rand Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            needrand = self.opt['Force']
            for mtype in self.list['MTypes']:
                mdb  = ppi.db(mtype)
                for field in ['AvF','pAvF','MeanEnr','pMeanEnr']:
                    needrand = needrand or field not in mdb.fields()
                    if field not in mdb.fields(): mdb.addField(field)
            if not needrand: return
            ### ~ [3] Calculate enrichment values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mtype in self.list['MTypes']:
                mdb  = ppi.db(mtype)
                for entry in mdb.entries():
                    mgenes = string.split(entry['Members'],'|')
                    for field in ['AvF','pAvF','MeanEnr','pMeanEnr']: entry[field] = 0.0
                    meanenr = []
                    for gene in mgenes:
                        if gene not in main.data(): meanenr.append(1.0); continue
                        entry['AvF'] += float(main.data()[gene]['AvF']) / len(mgenes)
                        entry['MeanEnr'] += float(main.data()[gene]['MeanEnr']) / len(mgenes)
                        meanenr.append(float(main.data()[gene]['MeanEnr']))
                    entry['MeanEnr'] = rje.geoMean(meanenr)
            ### ~ [4] Assess using randomised gene-data links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rx = 0.0
            for r in range(rnum):
                self.progLog('\r#RAND','Complex enrichment randomisations: %.2f%%' % (rx/rnum))
                ## ~ [3a] Generate randomisations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rdata = {}
                ridlist = rje.randomList(main.index('NRID').keys())
                mkeys = rje.randomList(main.dataKeys())
                while mkeys: rdata[ridlist.pop(0)] = main.data()[mkeys.pop(0)]
                ## ~ [3b] Compare to observations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for mtype in self.list['MTypes']:
                    mdb  = ppi.db(mtype)
                    for entry in mdb.entries():
                        rval = {}
                        for field in ['AvF','MeanEnr']: rval[field] = 0.0
                        meanenr = []
                        mgenes = string.split(entry['Members'],'|')
                        for gene in mgenes:
                            if gene not in rdata: meanenr.append(1.0); continue
                            rval['AvF'] += float(rdata[gene]['AvF']) / len(mgenes)
                            rval['MeanEnr'] += float(rdata[gene]['MeanEnr']) / len(mgenes)
                            meanenr.append(float(rdata[gene]['MeanEnr']))
                        rval['MeanEnr'] = rje.geoMean(meanenr)
                        for field in ['AvF','MeanEnr']:
                            if rval[field] >= entry[field] > 0: entry['p%s' % field] += (1.0 / rnum)
                            elif rval[field] <= entry[field] < 0: entry['p%s' % field] += (1.0 / rnum)
                            elif entry[field] == 0: entry['p%s' % field] += (1.0 / rnum)
            self.printLog('\r#RAND','Complex enrichment randomisations complete.')
            ### ~ [5] Finish formatting and resave ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mtype in self.list['MTypes']:
                mdb  = ppi.db(mtype)
                for entry in mdb.entries():
                    for field in ['AvF','pAvF','MeanEnr','pMeanEnr']: entry[field] = rje.expectString(entry[field])
                mdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],mdb.info['Name'].lower()))
        except: self.errorLog('HM Complexes error')
#########################################################################################################################
### End of SECTION II: HAPPI Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def htmlHead(title,stylesheets=['../example.css','../redwards.css'],tabber=True,frontpage=False,nobots=True):    ### Returns text for top of HTML file
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
    if nobots: html.insert(6,'<meta name="robots" content="index,nofollow">')
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
    try: HAPPI(mainlog,cmd_list).run()

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
