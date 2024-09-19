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
Module:       RJE_SLiMHTML
Description:  Module for generating HTML for Human SLiMFinder study
Version:      0.9
Last Edit:    11/03/11
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    datapath=PATH       : Path to parent data directory [./]
    htmlpath=PATH       : Path of parent html directory [./html]
    stylesheets=LIST    : List of CSS files to use ['../example.css','../redwards.css']
    border=X            : Border setting for tables [0]
    dropfields=LIST     : Fields to exclude from summary tables []
    fakehtml=T/F        : Whether to make UniFake HTML [True]
    unifake=PATH        : Path to UniFake dat file(s) [./unifake/]
    unipath=PATH        : Path to real UniProt dat file(s) [./uniprot/]
    unireal=LIST        : Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]
    makepng=T/F         : Whether to (look for and) make PNG files with R [True]
    svg=T/F             : Make SVG files rather than PNG files [True]
    addrand=T/F         : Whether to add pages for random data [True]
    makepages=LIST      : Types of pages to make [front,gene,domain,rand,slim,fake,occaln,interactome,slimaln,go,nested]
    titletext=FILE      : File containing (Page,ID,Title) [titletext.tdt]
    slimdesc=FILE       : File containing descriptions for SLiMs [slimdesc.txt]
    baddset=LIST        : List of Bad Datasets to filter out of main and occ database tables []
    pround=T/F          : Whether to round off occurrence positions for PNGs [True]
    xgcut=X             : Significance cut-off for XGMML [0.01]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_ppi, rje_slim, rje_uniprot, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time, math
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_ensembl, rje_go, rje_ppi, rje_seq, rje_sequence, rje_slim, rje_svg, rje_uniprot, rje_zen
import rje_dismatrix_V2, rje_tree, rje_xgmml
import happi, slimfinder
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.3 - Added code for making Random Dataset pages
    # 0.4 - Updated UPC pages and added additional front pages.
    # 0.5 - Split front page into front and full. Added GO tabs/pages.
    # 0.6 - Added XGMML output.
    # 0.7 - Modified output for HumSF10 and HAPPI analysis.
    # 0.8 - Added SVG output. Integrated better with HAPPI code.
    # 0.9 - Added SLiM Descriptions.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add GO Tabs/Pages.
    # [ ] : Integrate with HAPPI.
    # [ ] : Replace all graphics with SVG.
    # [ ] : http://bioware.ucd.ie/~compass/cgi-bin/formParser.py?name_server=slimsearch2&motif_str=PP.Y
    # [?] : http://bioware.ucd.ie/~compass/cgi-bin/formParser.py?name_server=slimsearch2&motif_str=PP.Y&proteinlist_list=P51170
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_SLiMHTML', '0.9', 'March 2011', '2010')
    description = 'RJE SLiM HTML Module'
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
class SLiMHTML(rje.RJE_Object):     
    '''
    SLiMHTML Class. Author: Rich Edwards (2010).

    Info:str
    - DataPath = Path to parent data directory [./]
    - HTMLPath = Path of parent html directory [./html]
    - SlimDesc = File containing descriptions for SLiMs [slimdesc.txt]
    - TitleText = File containing (Page,ID,Title) [titletext.tdt]
    - UniFake = Path to UniFake dat file(s) [./unifake/]
    - UniPath = Path to real UniProt dat file(s) [./uniprot/]
    
    Opt:boolean
    - AddRand = Whether to add pages for random data [True]
    - FakeHTML = Whether to make UniFake HTML [True]
    - MakePNG = Whether to (look for and) make PNG files with R [True]
    - PRound = Whether to round off occurrence positions for PNGs [True]
    - SVG = Make SVG files rather than PNG files [True]

    Stat:numeric
    - Border = Border setting for tables [0]
    - XGCut = Significance cut-off for XGMML [0.01]

    List:list
    - BadDSet = List of Bad Datasets to filter out of main and occ database tables []
    - DropFields = Fields to exclude from summary tables []
    - MakePages = Types of pages to make [front,gene,domain,rand,slim,fake]
    - UniReal = Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]

    Dict:dictionary    
    - SlimDesc = Dictionary of {Pattern:Description text (can be HTML)}

    Obj:RJE_Objects
    - DB = RJE_DB Database object
    - PPI = RJE_PPI PPI Object
    - GO = RJE_GO GO Object
    - HAPPI = HAPPI Object for visualisations.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['DataPath','HTMLPath','TitleText']
        self.optlist = ['AddRand','FakeHTML','MakePNG','SVG']
        self.statlist = ['Border','XGCut']
        self.listlist = ['StyleSheets','UniReal','MakePages','BadDSet']
        self.dictlist = ['GeneMap','SlimDesc']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'DataPath':rje.makePath('./'),'HTMLPath':rje.makePath('./html/'),'Basefile':'humsf_html',
                      'UniFake':rje.makePath('./HumSF09_UniFake/'),'UniPath':rje.makePath('./HumSF09_UniProt/'),
                      'TitleText':'titletext.tdt','SlimDesc':'slimdesc.txt'})
        self.setStat({'Border':0,'XGCut':0.01})
        self.setOpt({'FakeHTML':True,'MakePNG':True,'AddRand':True,'Iridis':False,'PRound':True,'SVG':True})
        self.list['StyleSheets'] = ['../example.css','../slimhtml.css']
        self.list['DropFields'] = ['Original','DType','Hub','PPISet','Annotation','Seq','KCloud','DCloud','Analysis',
                                   'MisMatch','Variant','wtFDR','SeqWt','HomNum_mean','GlobID_mean','LocID_mean',
                                   'BestELM','ELMSim','ELMPattern','Score','MatchIC','dConn','sConn','WtFDR',
                                   'AutoAnn']
        self.list['MakePages'] = rje.split('front,gene,domain,rand,slim,fake,occaln,interactome,slimaln,go',',')
        self.list['UniReal'] = rje.split('AC,GN,RC,RX,CC,DR,PE,KW,FT',',')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['PPI'] = rje_ppi.PPI(self.log,self.cmd_list)
        self.obj['SVG'] = rje_svg.SVG(self.log,self.cmd_list); self.obj['SVG'].setup()
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
                self._cmdReadList(cmd,'file',['TitleText','SlimDesc'])
                self._cmdReadList(cmd,'path',['DataPath','HTMLPath','UniProt','UniFake'])
                self._cmdReadList(cmd,'int',['Border'])
                self._cmdReadList(cmd,'stat',['XGCut'])
                self._cmdReadList(cmd,'opt',['FakeHTML','MakePNG','PRound','SVG'])
                self._cmdReadList(cmd,'list',['StyleSheets','DropFields','UniReal','MakePages','BadDSet'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.opt['Test']: self.sigMappingJIM(); return
            genelist = rje.sortUnique(list(self.db().getTable('PPI').index('Hub').keys()) + list(self.db().getTable('PPI').index('Spoke').keys()))
            if self.opt['AddRand']:
                randgene = rje.sortUnique(list(self.db().getTable('Occ_rseq').index('Spoke').keys()) + list(self.db().getTable('Occ_rupc').index('Spoke').keys()))
                genelist = rje.sortUnique(genelist+randgene)
            self.list['Genes'] = genelist[0:]
            rseqlist = self.db().getTable('Main_rseq').index('Dataset')
            rupclist = self.db().getTable('Main_rupc').index('Dataset')
            domainlist = self.db().getTable('DomPPI').index('Domain').keys()
            if self.opt['AddRand']: slimlist = self.db().getTable('Occ').index('Pattern').keys()
            else: slimlist = self.db().getTable('Occ_real').index('Pattern').keys()
            hpath = self.info['HTMLPath']
            if 'xgmml' in self.list['MakePages']: self.xgmmlOut()
            if 'happi' in self.list['MakePages']: self.happiOut()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#if self.opt['FakeHTML']:
            if 'fake' in self.list['MakePages']: self.uniFakeHTML()
            if 'front' in self.list['MakePages']:
                html = self.frontPage(genelist,domainlist,slimlist)
                hfile = hpath + 'index.htm'
                open(hfile,'w').write(html['front'])
                hfile = hpath + 'index.html'
                open(hfile,'w').write(html['front'])
                hfile = hpath + 'full.htm'
                open(hfile,'w').write(html['full'])
            ## ~ [2a] Gene Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'gene' in self.list['MakePages']:
                gx = 0.0; gtot = len(genelist)
                for gene in genelist:
                    self.progLog('\r#GENE','Making Gene Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.genePage(gene)
                self.printLog('\r#GENE','%s Gene HTML Pages made.' % rje.integerString(gtot))
            ## ~ [2b] Domain Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'domain' in self.list['MakePages']:
                gx = 0.0; gtot = len(domainlist)
                for domain in domainlist:
                    self.progLog('\r#DOM','Making Domain Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'domain/'),domain)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.domainPage(domain)
                self.printLog('\r#DOM','%s Domain HTML Pages made.' % rje.integerString(gtot))
            ## ~ [2c] SLiM Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'slim' in self.list['MakePages']:
                gx = 0.0; gtot = len(slimlist)
                for slim in slimlist:
                    self.progLog('\r#SLIM','Making SLiM Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'slim/'),slim)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.slimPage(slim)
                self.printLog('\r#SLIM','%s SLiM HTML Pages made.' % rje.integerString(gtot))
                self.slimFTPage(genelist,domainlist,slimlist)
            elif 'slimft' in self.list['MakePages']: self.slimFTPage(genelist,domainlist,slimlist)
            ## ~ [2d] Make Significant Random Dataset Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'rand' in self.list['MakePages']:
                randlist = self.db().getTable('Sig_rand').index('Dataset')
                gx = 0.0; gtot = len(randlist)
                for rdset in randlist:
                    self.progLog('\r#RAND','Making Random Dataset Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + '%s/' % rdset[:4]),rdset)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.randPage(rdset)
                self.printLog('\r#RAND','%s Random Dataset HTML Pages made.' % rje.integerString(gtot))
            ## ~ [2e] GO Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'go' in self.list['MakePages']:
                gdb = self.db().getTable('GO')
                gx = 0.0; gtot = len(gdb.index('GO_Desc'))
                for go in rje.sortKeys(gdb.index('GO_Desc')):
                    self.progLog('\r#GO','Making GO Pages: %.2f%%' % (gx/gtot)); gx += 100.0
                    goid = gdb.dataList(gdb.indexEntries('GO_Desc',go),'GO_ID')[0]
                    hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
                    if not self.opt['Force'] and self.checkHTML(hfile): continue
                    self.goPage(go,goid)
                self.printLog('\r#GO','%s GO HTML Pages made.' % rje.integerString(gtot))
            ### ~ [3] Visualisations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forkpng = [] 
            if self.opt['MakePNG']:
                self.setupHAPPI(); hobj = self.obj['HAPPI']
                tdb = self.db().getTable('TP')  # Known Hub Type
                kdb = self.db().getTable('Knowns')  # Known Pattern
                inthub = tdb.dataList(tdb.indexEntries('Type','TP'),'Hub') + tdb.dataList(tdb.indexEntries('Known','!'),'Hub') + tdb.dataList(tdb.indexEntries('Known','!!'),'Hub')
                intpat = []
                for known in tdb.dataList(tdb.indexEntries('Type','TP'),'Known') + ['!','!!','?','?!']:
                    intpat += kdb.dataList(kdb.indexEntries('Known',known),'Pattern')
                intpat = rje.sortUnique(intpat)
                if 'nested' in self.list['MakePages']:
                    for pattern in rje.sortKeys(self.db().getTable('Sig').index('Pattern')):
                        if pattern not in intpat: continue
                        if hobj.opt['NoForks']: self.slimJimMotif(pattern)
                        else: forkpng.append((pattern,'nested'))                       
                    for pattern in rje.sortKeys(self.db().getTable('Sig').index('Pattern')):
                        if pattern in intpat: continue
                        if hobj.opt['NoForks']: self.slimJimMotif(pattern)
                        else: forkpng.append((pattern,'nested'))                       
                if 'occaln' in self.list['MakePages']: self.sigMappingJIM()
                if 'interactome' in self.list['MakePages']:
                    for dataset in rje.sortKeys(self.db().getTable('Main_real').index('Dataset')):
                        if rje.split(dataset,'.')[0] in inthub:
                            if hobj.opt['NoForks']: self.slimJimHub(dataset)
                            else: forkpng.append((dataset,'interactome'))                       
                    for dataset in rje.sortKeys(self.db().getTable('Main_real').index('Dataset')):
                        if rje.split(dataset,'.')[0] not in inthub:
                            if hobj.opt['NoForks']: self.slimJimHub(dataset)
                            else: forkpng.append((dataset,'interactome'))                       
                    for dataset in rje.sortKeys(self.db().getTable('Sig_rand').index('Dataset')): 
                            if hobj.opt['NoForks']: self.slimJimHub(dataset)
                            else: forkpng.append((dataset,'interactome'))                       
                if 'slimaln' in self.list['MakePages']:
                    for pattern in rje.sortKeys(self.db().getTable('Sig').index('Pattern')):
                        if pattern not in intpat: continue
                        for entry in self.db().getTable('Sig').indexEntries('Pattern',pattern):
                            if hobj.opt['NoForks']: self.slimJimHubMotif(entry)
                            else: forkpng.append((entry,'slimaln'))                                                  
                    for pattern in rje.sortKeys(self.db().getTable('Sig').index('Pattern')):
                        if pattern in intpat: continue
                        for entry in self.db().getTable('Sig').indexEntries('Pattern',pattern): 
                            if hobj.opt['NoForks']: self.slimJimHubMotif(entry)
                            else: forkpng.append((entry,'slimaln'))
                try: self.forkPNG(forkpng)
                except SystemExit: os._exit(0)
                except: raise
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def forkPNG(self,forkpng):  ### Fork out PNG generation (lists of (id,type) tuples)
        '''Fork out PNG generation (lists of (id,type) tuples).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not forkpng: return
            hobj = self.obj['HAPPI']
            forkx = hobj.stat['Forks']      # Number of forks to have running at one time
            forks = {}      # Dictionary of active forks {pid:basefile}
            killforks = hobj.stat['KillForks']      # Time in seconds to wait after main thread has apparently finished
            killtime = time.time(); sleepsec = 1
            pngx = len(forkpng)
            self.obj['PPI'].opt['ProgLog'] = False
            self.db('Occ').index('Dataset'); self.db('Occ').index('Pattern'); self.db('Occ').index('Hub')
            self.db('PPI').index('SpokeSeq'); self.db('PPI').index('Hub')
            self.db('DomPPI').index('Domain')
            
            ### ~ [1] ~ Cycle through forkpng tuples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while forks or forkpng:
                ## ~ [1a] ~ Add more forks if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while forkpng and (len(forks) < forkx):     # Add more forks
                    (id,mtype) = forkpng.pop(0)
                    if mtype == 'nested': basefile = 'tmp_%s-%s' % (mtype,rje_slim.slimFromPattern(id))
                    elif mtype == 'slimaln': basefile = 'tmp_%s-%s-%s' % (mtype,rje_slim.slimFromPattern(id['Pattern']),id['Spoke'])
                    else: basefile = 'tmp_%s-%s' % (mtype,id)
                    forkcmd = self.cmd_list + ['i=-1','log=%s.log' % basefile,'errorlog=None']
                    try: newpid = os.fork()
                    except:
                        self.errorLog('Problem forking %s.' % id,printerror=False); 
                        forkpng.append((id,mtype))
                        sleepsec *= 2; break
                    sleepsec = 1
                    if newpid == 0: # child
                        self.opt['Child'] = True
                        forkobj = SLiMHTML(cmd_list=forkcmd)
                        forkobj.log.info['LogFile'] = '%s.log' % basefile
                        forkobj.obj = self.obj
                        forkobj.info = self.info
                        forkobj.stat = self.stat
                        forkobj.opt = self.opt
                        forkobj.list = self.list
                        forkobj.dict = self.dict
                        if mtype == 'nested': forkobj.slimJimMotif(id)
                        elif mtype == 'interactome': forkobj.slimJimHub(id)
                        elif mtype == 'slimaln': forkobj.slimJimHubMotif(id)
                        os._exit(0)    # Exit process 
                    elif newpid == -1: self.errorLog('Problem forking %s.' % id,printerror=False)  
                    else: forks[newpid] = basefile
                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(sleepsec)       # Sleep for 1s 
                forklist = hobj._activeForks(forks.keys())
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
    def geneLink(self,gene,frontpage=False,altgene=None):     ### Returns gene link text
        '''Returns gene link text.'''
        if not altgene: altgene = gene
        if gene in self.list['Genes']: 
            if frontpage: return '<a href="./gene/%s.html" target="_blank" title="%s">%s</a>' % (gene,self.geneTitle(gene),altgene)
            else: return '<a href="../gene/%s.html" title="%s">%s</a>' % (gene,self.geneTitle(gene),altgene)
        return '<a title="%s">%s</a>' % (self.geneTitle(gene),altgene)
#########################################################################################################################
    def geneTitle(self,gene):   ### Returns mouseover text for gene
        '''Returns mouseover text for gene.'''
        desc = self.geneDesc(gene)
        if desc and gene in self.list['Genes']: return '%s (Click for %s results page)' % (desc,gene)
        elif gene in self.list['Genes']: return '(Click for %s results page)' % (gene)
        elif desc: return '%s (No results page for %s)' % (desc,gene)
        else: return '(No results page for %s)' % (gene)
#########################################################################################################################
    def geneDesc(self,gene):    ### Returns gene description
        '''Returns gene description.'''
        dbxref = self.db().getTable('DBXRef')
        for field in ['Gene','UniProt','EnsEMBL','EnsLoci']:
            if gene in dbxref.index(field): return dbxref.indexEntries(field,gene)[0]['EnsDesc']
        return ''
#########################################################################################################################
    def setupHAPPI(self):   ### Sets up HAPPI object for visualisations etc.
        '''Sets up HAPPI object for visualisations etc.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hobj = self.obj['HAPPI'] = happi.HAPPI(self.log,self.cmd_list)
            hobj.obj['DB'] = self.obj['DB']
            ### ~ [1] ~ Setup PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pmap = self.dict['PPIMap'] = {}     # Mapping of {ppi_node:hub}
            pobj = self.obj['PPI']; ppi = pobj.dict['PPI'] = {}
            hobj.obj['PPI'] = pobj
            pdb = self.db('PPI')    #Hub     HubUni  Spoke   SpokeUni        SpokeSeq        Evidence        ppi     bin     com     y2h
            ddb = self.db('DomPPI') #Domain  Spoke   SpokeUni        SpokeSeq        Hubs    Evidence        ppi     bin     com     y2h
            rdb = self.db('DSetSeq')#Dataset   Seq    (self.gene(seq))
            fdb = self.db('PFam')   #Type    Name    Start   End     Eval    Score
            for entry in pdb.entries():
                (hub,spoke,evidence) = (entry['Hub'],entry['Spoke'],entry['Evidence'])
                pmap[hub] = hub
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
                (hub,spoke) = (spoke,hub)
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
            for entry in ddb.entries():
                (hub,spoke,evidence) = (entry['Domain'],entry['Spoke'],entry['Evidence'])
                hub = 'd:%s' % hub; pmap[hub] = entry['Domain']
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
                (hub,spoke) = (spoke,hub)
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
            for entry in rdb.entries():
                try: (hub,spoke,evidence) = (entry['Dataset'],self.gene(entry['Seq']),'Random')
                except: continue
                hub = rje.join(rje.split(hub,'_')[:2],'_'); pmap[hub] = entry['Dataset']
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
                #self.deBug('%s: %s' % (hub,ppi[hub]))
                (hub,spoke) = (spoke,hub)
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
            for entry in fdb.entries():
                try: (hub,spoke,evidence) = (entry['Type'],self.gene(entry['Name']),'Domain')
                except: continue
                hub = 'd:%s' % hub; pmap[hub] = entry['Type']
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence
                (hub,spoke) = (spoke,hub)
                if hub not in ppi: ppi[hub] = {}
                if spoke not in ppi[hub]: ppi[hub][spoke] = evidence            
            try: hobj.opt['NoForks'] = hobj.opt['Win32'] or hobj.opt['NoForks'] or hobj.stat['Forks'] < 1
            except: hobj.opt['NoForks'] = True
            self.printLog('#HAPPI','HAPPI Object setup (Forking=%s)' % (not hobj.opt['NoForks']))
        except: self.errorLog('Fook')
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pick = self.unpickleMe()
            if pick: self = pick; return
            dpath = rje.makePath(self.info['DataPath'] + 'HumSF09_MainData')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)

            dfile = dpath + 'humsf09.pairwise_ppi.090505.tdt'
            ppi = self.db().addTable(dfile,['Hub','Spoke'],'All',name='PPI')
            #ppi.index('Hub'); ppi.index('Spoke')
            #Hub     HubUni  Spoke   SpokeUni        SpokeSeq        Evidence        ppi     bin     com     y2h
            #1082356 -       GSTK1   Q9Y2Q3  GSTK1_HUMAN__Q9Y2Q3     IntAct:anti bait coip   Y       Y       -       -

            dfile = dpath + 'humsf09.genemap.0905050.tdt'
            dbxref = self.db().addTable(dfile,['Gene'],'All',name='DBXRef')
            #Gene    Entrez  HPRD    OMIM    UniProt EnsEMBL EnsLoci EnsDesc ppi     bin     com     y2h
            #A1BG    1       00726   138670  P04217  ENSG00000121410 A1BG_HUMAN__ENSP00000263100     Alpha-1B-glycoprotein Precursor (Alpha-1-B glycoprotein)        1       1       0       0

            ## ~ [1a] Occurrence table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dfile = dpath + 'humsf09_slimfinder.occ.csv'
            occ = self.db().addTable(dfile,['Dataset','Rank','Pattern','Seq','Start_Pos'],'All',name='Occ')
            self.printLog('#BAD','%d Bad datasets to drop from results' % len(self.list['BadDSet']))
            occ.dropEntriesDirect('Dataset',self.list['BadDSet']) #?# Move and expand #?#
            if 'Spoke' not in occ.fields():
                occ.addField('Spoke'); occ.addField('DType')
                occ.list['Fields'].remove('Spoke')
                occ.list['Fields'].insert(occ.list['Fields'].index('Sig')+1,'Spoke')
                for entry in occ.entries():
                    if entry['Dataset'][:4] in ['rseq','rupc']: entry['DType'] = entry['Dataset'][:4]
                    else: entry['DType'] = 'real'
                    entry['Spoke'] = self.gene(entry['Seq'])
                    #try: entry['Spoke'] = dbxref.index('EnsLoci')[entry['Seq']][0]
                    #except:
                    #    if entry['Dataset'][:4] in ['rseq','rupc']: entry['Spoke'] = '-'; continue
                    #    self.errorLog('Cannot find gene for spoke "%s"' % entry['Seq'])
                    #    entry['Spoke'] = '!ERR!'
                occ.index('Spoke')
            for entry in occ.entries():
                entry['Rank'] = rje.preZero(int(entry['Rank']),1000)
                if entry['Spoke'] == '-': entry['Spoke'] = entry['Seq']
            occ.remakeKeys()
            #occ.index('Hub'); occ.index('Dataset'); occ.index('Pattern');
            if 'Spoke' in occ.fields():
                #occ.index('Spoke')
                occ.list['Fields'].remove('Spoke')
                occ.list['Fields'].insert(occ.list['Fields'].index('Sig')+1,'Spoke')
            #Hub,Dataset,Original,Rank,Pattern,Sig,Spoke,Seq,Start_Pos,End_Pos,Prot_Len,Match,Variant,MisMatch,Desc,Cons,HomNum,GlobID,LocID,Hyd,IUP,SA,PepSeq,PepDesign
            #12962935,12962935.com,12962935.comppi,1,G..G.GKT,0.0107820687891,DHX37,DHX37_HUMAN__Q8IY37,274,281,1157,GETGSGKT,G..G.GKT,0,Probable ATP-dependent RNA helicase DHX37 (EC 3.6.1.-)(DEAH box protein 37) [acc:Q8IY37 pep:ENSP00000311135 gene:ENSG00000150990],0.786344651024,25,0.675989628349,1.0,-0.135,0.2788875,1.03218724224,GETGSGKT,OK
            self.db().splitTable(occ,'DType')
            occ_rand = self.db().copyTable('Occ_rseq','Occ_rand')
            occ_temp = self.db().copyTable('Occ_rupc','Occ_temp')
            occ_rand = self.db().mergeTables(occ_rand,occ_temp)
    
            dfile = dpath + 'humsf09.domain_ppi.090505.tdt'
            domppi = self.db().addTable(dfile,['Domain','Spoke'],'All',name='DomPPI')
            #domppi.index('Domain'); domppi.index('Spoke')
            #Domain  Spoke   SpokeUni        SpokeSeq        Hubs    Evidence        ppi     bin     com     y2h
            #1-cysPrx_C      PPP2R1A P30153  2AAA_HUMAN__P30153      PRDX1,PRDX2     IntAct:anti tag coip    -       N       -       N

            dfile = dpath + 'humsf09.go.090505.tdt'
            godb = self.db().addTable(dfile,['EnsG','GO_ID'],'All',name='GO')
            self.dict['GO'] = {}
            for gkey in godb.data():
                (ensg,go) = rje.split(gkey)
                gtype = godb.data()[gkey]['GO_Type'].upper()
                if ensg not in self.dict['GO']: self.dict['GO'][ensg] = {}
                if gtype not in self.dict['GO'][ensg]: self.dict['GO'][ensg][gtype] = []
                self.dict['GO'][ensg][gtype].append((go,godb.data()[gkey]['GO_Desc']))
            #EnsG    GO_ID   GO_Type GO_Desc
            #ENSG00000000003 0007186 bp      G-protein coupled receptor protein signaling pathway
            go = self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            go.readGO()#; go.mapEnsGO()
            
            dfile = dpath + 'humsf09_slimfinder.full.csv'
            main = self.db().addTable(dfile,['Dataset','Rank','Pattern'],'All',name='Main')
            replacements = []
            for dset in self.list['BadDSet']: replacements += main.indexEntries('Dataset',dset)[:1]
            self.printLog('#BAD','%d Bad datasets to drop from results' % len(self.list['BadDSet']))
            main.dropEntriesDirect('Dataset',self.list['BadDSet']) #?# Move and expand #?#
            for rentry in replacements:
                for f in ['IC','ExpUP','Prob','Sig','Cons_mean','HomNum_mean','GlobID_mean','LocID_mean','Hyd_mean','IUP_mean','SA_mean','Occ','Support','UP','CloudSeq','CloudUP']: rentry[f] = 0.0
                for f in ['MotNum','Rank']: rentry[f] = 0
                rentry['Pattern'] = '-'
                rentry['Sig'] = 0.05
                rentry['KCloud'] = '%s|' % rentry['Dataset']
                rkey = main.makeKey(rentry)
                main.dict['Data'][rkey] = rentry
            main.dict['Index'] = {}
            for entry in main.entries(): entry['Rank'] = rje.preZero(int(entry['Rank']),1000)
            main.remakeKeys()
            main.renameField('SeqNum','SeqN'); main.renameField('UPNum','UPN'); main.renameField('AANum','AA')
            main.renameField('ELM','CM')
            sig = self.db().copyTable('Main','Sig')
            sig.dropEntries(['Rank<1'])
            #main.index('Pattern'); main.index('Hub'); main.index('Dataset'); main.index('Rank')
            #sig.index('Pattern'); sig.index('Hub'); sig.index('Dataset'); sig.index('Rank')
            #Dataset,RunID,Masking,Build,Chance,RunTime,SeqNum,UPNum,AANum,MotNum,Rank,Sig,FDR,Pattern,IC,Occ,Support,UP,ExpUP,Prob,Cloud,CloudSeq,CloudUP,Cons_mean,HomNum_mean,GlobID_mean,LocID_mean,Hyd_mean,IUP_mean,SA_mean,Original,DType,Hub,PPISet,Annotation,Class,ELM,BestELM,ELMSim,ELMPattern,Score,Match,MatchIC,KCloud,DCloud,Analysis
            #1-cysPrx_C.bindom,newbin,FreqConsDisComp-5-8,l5w2o2a1.FreqConsDisComp-5-8,Sig,00:01:28,23,16,3421,0,0,0.05,,-,,,,,,,,,,,,,,,,,1-cysPrx_C.bindom,real,1-cysPrx_C,bindom,,,0,,,,0,,0,1-cysPrx_C.bindom|,1-cysPrx_C|,real.bindom
            self.db().splitTable(main,'DType'); self.db().splitTable(sig,'DType')
            sig_rand = self.db().copyTable('Sig_rseq','Sig_rand')
            sig_temp = self.db().copyTable('Sig_rupc','Sig_temp')
            sig_rand = self.db().mergeTables(sig_rand,sig_temp)

            #dfile = dpath + 'humsf09.seqwt.tdt'
            #seqwt = self.db().addTable(dfile,['Seq','RType','PType'],'All',name='SeqWt')
            #self.db().splitTable(seqwt,'RType')

            dfile = dpath + 'humsf09.dsetseq.tdt'
            dsetseq = self.db().addTable(dfile,['Dataset','Seq'],'All',name='DSetSeq')
            dsetseq.dropEntriesDirect('RType',['real'])
            #x#self.db().splitTable(dsetseq,'RType')

            dfile = self.info['UniFake'] + 'ens_HUMAN.unifake.pfam.tdt'
            pfam = self.db().addTable(dfile,['Type','Name','Start','End'],'All',name='PFam')
            #pfam.index('Type'); pfam.index('Name')
            #Type    Name    Start   End     Eval    Score
            #Dysbindin       DBND1_HUMAN__ENSP00000002501    14      153     2.30e-92        317.8
            #http://pfam.sanger.ac.uk/family/Dysbindin

            dfile = dpath + 'humsf09.tp.tdt'
            dtp = self.db().addTable(dfile,['Known','Hub'],'All',name='TP')
            dfile = dpath + 'humsf09.knowns.tdt'
            dtp = self.db().addTable(dfile,['Pattern','Known'],'All',name='Knowns')

            comfile = dpath + 'humsf09.elmcomp.compare.tdt'
            cdb = self.db().addTable(comfile,mainkeys=['Name1','Name2'],name='CompariMotif')
            #for bad in self.list['ScreenELM']:
            #    self.printLog('#ELM','Screened out hits from ELM "%s"' % bad)
            #    cdb.dropEntries(['Name2==%s' % bad])
            #cdb.index('Motif1'); cdb.index('Motif2'); cdb.index('Name2')
            cdb.rankFieldByIndex('Motif1','Score','Rank',rev=True,absolute=True,lowest=True)
            cdb.deleteField('MotifFile'); cdb.deleteField('SearchDB'); cdb.deleteField('Motif1'); cdb.deleteField('Sim1')
            cdb.renameField('Name1','Pattern')
            cdb.renameField('Name2','SLiM')
            cdb.renameField('Sim2','Relationship')
            cdb.renameField('Desc2','Description')
            cdb.renameField('Motif2','Motif')
            cdb.newKey(['Pattern','Rank','SLiM'],startfields=True)

            afile = dpath + 'humsf09.annotation_links.tdt'
            adb = self.db().addTable(afile,mainkeys=['Annotation'],name='Annotation')
            if not adb:
                self.deBug(rje.dataDict(self,'humsf09.annotation_links.tdt',mainkeys=['Annotation'],getheaders=True))
                raise ValueError

            if self.info['SlimDesc'].lower() not in ['false','f']: self.setupSlimDesc()

            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['HTMLPath'])
            subdirlist = ['gene','domain','slim']
            if self.opt['AddRand']: subdirlist += ['rseq','rupc']
            for subdir in subdirlist: rje.mkDir(self,rje.makePath(self.info['HTMLPath']+subdir))
            #self.pickleMe()          
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def setupSlimDesc(self):    ### Setup SLiM Descriptions
        '''Setup SLiM Descriptions.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sd = self.dict['SlimDesc']
            ### ~ [1] Load existing descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(self.getStr('SlimCheck')):
                for line in open(self.getStr('SlimCheck'),'r').readlines():
                    self.progLog('\r#DESC','Loading %s SLiM descriptions' % rje.iLen(sd))
                    try:
                        pattern = rje.split(line)[0]
                        desc = rje.join(rje.split(rje.chomp(line))[1:])
                    except: continue
                    while desc[-1:] == ' ': desc = desc[:-1]
                    if not (pattern and desc): continue
                    if desc[-1] != '.': desc = '%s.' % desc
                    if pattern not in sd: sd[pattern] = desc
                    else: sd[pattern] = '%s %s' % (sd[pattern],desc)
                self.printLog('\r#DESC','Loaded %s SLiM descriptions.' % rje.iLen(sd))
            sx = len(sd)
            ### ~ [2] Auto-generate remaining descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            kdb = self.db('Knowns'); tdb = self.db('TP'); adb = self.db('Annotation')
            ktext = {'leak?':'BLAST leakage artefact','Nterm':'N-terminal motif','Cterm':'C-terminal motif'}
            for k in ktext: adb.addEntry({'Annotation':k,'Description':ktext[k],'Link':''})
            ## ~ [2a] ~ Known motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in kdb.entries():
                pattern = entry['Pattern']
                known = entry['Known']
                desc = ''
                if known in tdb.index('Known'): type = tdb.indexDataList('Known',known,'Type')
                else: type = []
                if '?' in known: desc = 'Possible '; known = known.replace('?','')
                if '!' in known: desc = '%sInteresting ' % desc; known = known.replace('!','')
                dtype = ''
                if 'TP' in type and 'OT' in type: dtype = 'True Positive/Off-target'
                elif 'TP' in type: dtype = 'True Positive'
                elif 'OT' in type: dtype = 'Off-target'
                elif 'FP' in type: dtype = 'False Positive'
                elif 'RAND' in type: dtype = 'Random dataset'
                elif 'Generic' not in type: dtype = 'Novel/Unknown'
                if dtype: desc = '%s%s ' % (desc,dtype)
                if known == 'CAAXbox': known = 'MOD_CAAXbox'
                if known in adb.index('Annotation'):
                    anentry = adb.data(known)
                    desc = '%s%s' % (desc,anentry['Description'])
                    if anentry['Link'] and anentry['Link'] != '-': desc = '%s [<a href="%s">linkout</a>]' % (desc,anentry['Link'])
                elif desc: desc = '%s%smotif'
                if desc:
                    if pattern not in sd: sd[pattern] = ''
                    if sd[pattern]: sd[pattern] += '\n<br>\n'
                    sd[pattern] += desc
            ## ~ [2b] ~ Novel and random motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            unslimlist = self.loadFromFile('humsf09.unannotated.txt',chomplines=True)
            if '-' in unslimlist: unslimlist.remove('-')
            for pattern in unslimlist:
                if pattern not in sd: sd[pattern] = 'Auto-annotated possible novel/unknown motif'
        except: self.errorLog('Setup SLiM Description error')            
#########################################################################################################################
    def logSig(self,sig,cap=-5):   ### Returns log10(Sig) capped at -5.
        sig = float(sig)
        if not sig: return cap
        return max(cap,math.log(sig,10))
#########################################################################################################################
    def xgmmlOut(self): ### Summary XGMML output
        '''Summary XGMML output.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pcut = self.stat['XGCut']
            if pcut > 0: self.printLog('#PCUT','XGMML output limited to results p <= %s' % pcut)
            xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            xgmml.setInfo({'Name':'HumSF09','Description':'SLiMFinder results from 2009 human proteome analysis.',
                           'Type':'SF-PPI'})
            edges = xgmml.dict['Edge'] = {}
            #- Edge = Dictionary of edges between nodes {Type:{(source,target):Attributes}}
            #- EdgeAtt = Dictionary of edge attributes {Att:Type}
            nodes = xgmml.dict['Node'] = {}
            #- Node = Dictionary of Nodes to be output {Node:Attributes}
            #- NodeAtt = Dictionary of node attributes {Att:Type}
            ppi = self.db().getTable('PPI')
            #Hub     HubUni  Spoke   SpokeUni        SpokeSeq        Evidence        ppi     bin     com     y2h
            #1082356 -       GSTK1   Q9Y2Q3  GSTK1_HUMAN__Q9Y2Q3     IntAct:anti bait coip   Y       Y       -       -
            dbxref = self.db().getTable('DBXRef')
            #Gene    Entrez  HPRD    OMIM    UniProt EnsEMBL EnsLoci EnsDesc ppi     bin     com     y2h
            occ = self.db().getTable('Occ')
            #Hub,Dataset,Original,Rank,Pattern,Sig,Spoke,Seq,Start_Pos,End_Pos,Prot_Len,Match,Variant,MisMatch,Desc,Cons,HomNum,GlobID,LocID,Hyd,IUP,SA,PepSeq,PepDesign
            #12962935,12962935.com,12962935.comppi,1,G..G.GKT,0.0107820687891,DHX37,DHX37_HUMAN__Q8IY37,274,281,1157,GETGSGKT,G..G.GKT,0,Probable ATP-dependent RNA helicase DHX37 (EC 3.6.1.-)(DEAH box protein 37) [acc:Q8IY37 pep:ENSP00000311135 gene:ENSG00000150990],0.786344651024,25,0.675989628349,1.0,-0.135,0.2788875,1.03218724224,GETGSGKT,OK
            sig = self.db().getTable('Sig')
            #Dataset,RunID,Masking,Build,Chance,RunTime,SeqNum,UPNum,AANum,MotNum,Rank,Sig,FDR,Pattern,IC,Occ,Support,UP,ExpUP,Prob,Cloud,CloudSeq,CloudUP,Cons_mean,HomNum_mean,GlobID_mean,LocID_mean,Hyd_mean,IUP_mean,SA_mean,Original,DType,Hub,PPISet,Annotation,Class,ELM,BestELM,ELMSim,ELMPattern,Score,Match,MatchIC,KCloud,DCloud,Analysis
            #1-cysPrx_C.bindom,newbin,FreqConsDisComp-5-8,l5w2o2a1.FreqConsDisComp-5-8,Sig,00:01:28,23,16,3421,0,0,0.05,,-,,,,,,,,,,,,,,,,,1-cysPrx_C.bindom,real,1-cysPrx_C,bindom,,,0,,,,0,,0,1-cysPrx_C.bindom|,1-cysPrx_C|,real.bindom
            ## ~ [0a] ~ Set up attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for x in ['Gene','Entrez','HPRD','OMIM','UniProt','EnsEMBL','EnsLoci','EnsDesc','Real','RUPC','RSeq','Hub','Spoke','Type','TP','OT','NVL']: xgmml.dict['NodeAtt'][x] = 'string'
            for x in ['Num','PPI','Occ']: xgmml.dict['NodeAtt'][x] = 'real'
            for x in ['Evidence','Start_Pos']: xgmml.dict['EdgeAtt'][x] = 'string'
            # Edge Evidence = PPI Evidence or List of DClouds (Occ/Clouds) or mediating proteins (DPI)
            # Edge Type = com/y2h/com-y2h/ppi/occ/rseq/rupc/slim
            for x in ['Sig']: xgmml.dict['EdgeAtt'][x] = 'real'
            ### ~ [1] ~ Add Proteins, Domains, SLiMs and Random Datasets as Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = occ.entryNum()
            for entry in occ.entries():
                self.progLog('\r#NODE','Generating XGMML nodes: %.2f%%' % (ex/etot)); ex += 100.0
                if float(entry['Sig']) > pcut > 0: continue
                hub = entry['Hub']; spoke = entry['Spoke']; pattern = entry['Pattern']
                ## Hub ##
                if hub not in nodes:
                    nodes[hub] = {'Spoke':'N','Occ':len(occ.indexDataList('Spoke',hub,'Pattern',sortunique=True)),
                                  'Real':'N','RUPC':'N','RSeq':'N','Type':'Protein','OT':'N','TP':'N','NVL':'N',
                                  'Col':39}
                    if hub in self.data('DBXRef'): nodes[hub] = rje.combineDict(nodes[hub],self.data('DBXRef')[hub])
                    nodes[hub]['PPI'] = len(ppi.indexEntries('Hub',hub))
                    if hub[:4].lower() in ['rupc','rseq']:
                        nodes[hub]['PPI'] = int(rje.split(hub,'-')[-1])
                        nodes[hub]['Type'] = hub[:4]
                        nodes[hub]['Col'] = {'rupc':44,'rseq':45}[hub[:4].lower()]
                    elif entry['Dataset'][-3:] == 'dom': nodes[hub]['Type'] = 'Domain'; nodes[hub]['Col'] = 38
                nodes[hub]['Hub'] = 'Y'
                ## Spoke ##
                if spoke not in nodes:
                    nodes[spoke] = {'Hub':'N','Occ':len(occ.indexDataList('Spoke',spoke,'Pattern',sortunique=True)),
                                    'Real':'N','RUPC':'N','RSeq':'N','Type':'Protein','OT':'N','TP':'N','NVL':'N',
                                    'Col':39}
                    if hub in self.data('DBXRef'): nodes[hub] = rje.combineDict(nodes[hub],self.data('DBXRef')[hub])
                    nodes[spoke]['PPI'] = len(ppi.indexEntries('Hub',spoke))
                ## SLiMs ##
                if pattern not in nodes:
                    nodes[pattern] = {'Spoke':'N','Hub':'N','Col':47,
                                      'Occ':len(occ.indexDataList('Pattern',pattern,'Spoke',sortunique=True)),
                                      'PPI':len(occ.indexDataList('Pattern',pattern,'Hub',sortunique=True)),
                                      'Real':'N','RUPC':'N','RSeq':'N','Type':'SLiM','OT':'N','TP':'N','NVL':'N'}
                ## RTypes ##
                if hub[:4].lower() == 'rupc': nodes[hub]['RUPC'] = 'Y'; nodes[spoke]['RUPC'] = 'Y'; nodes[pattern]['RUPC'] = 'Y';
                elif hub[:4].lower() == 'rseq': nodes[hub]['RSeq'] = 'Y'; nodes[spoke]['RSeq'] = 'Y'; nodes[pattern]['RSeq'] = 'Y';
                else: nodes[hub]['Real'] = 'Y'; nodes[spoke]['Real'] = 'Y'; nodes[pattern]['Real'] = 'Y';
                ## Class ##
                for sclass in sig.indexDataList('Pattern',pattern,'Class'):
                    if sclass in ['OT','TP']: nodes[hub][sclass] = 'Y'; nodes[spoke][sclass] = 'Y'; nodes[pattern][sclass] = 'Y';
                    elif len(sclass) < 2: nodes[hub]['NVL'] = 'Y'; nodes[spoke]['NVL'] = 'Y'; nodes[pattern]['NVL'] = 'Y';
            self.printLog('\r#NODE','Generated %s XGMML nodes.  ' % rje.integerString(len(nodes)))
            ## ?? Should non-significant protein hubs and domains be added too ?? ##
            ### ~ [2] ~ Add PPI, Occurrences and Cloud relationships as Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for etype in rje.split('com/y2h/com-y2h/ppi/dmi/occ/rseq/rupc/slim/pfam/cloud','/'): edges[etype] = {}
            ## ~ [2a] ~ Add PPI as Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            px = 0.0; ptot = ppi.entryNum()
            for pp in ppi.entries():
                self.progLog('\r#EDGE','Making PPI edges: %.2f%%' % (px/ptot)); px += 100.0
                hub = pp['Hub']; spoke = pp['Spoke']
                if pp['Evidence'].lower().find('complex') > 0 and pp['Evidence'].lower().find('yeast') > 0: etype = 'com-y2h'
                elif pp['Evidence'].lower().find('complex') > 0: etype = 'com'
                elif pp['Evidence'].lower().find('yeast') > 0: etype = 'y2h'
                else: etype = 'ppi'
                if spoke not in nodes or hub not in nodes: continue
                if (spoke,hub) in edges[etype]: continue
                edges[etype][(hub,spoke)] = rje.combineDict({'Sig':-len(rje.split(pp['Evidence'],'|'))},pp)
            ## ~ [2a] ~ Add Occurrences as Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ex = 0.0; etot = occ.entryNum()
            for entry in occ.entries():
                self.progLog('\r#EDGE','Making Occurrence edges: %.2f%%' % (ex/etot)); ex += 100.0
                if entry['Sig'] > pcut > 0: continue
                hub = entry['Hub']; spoke = entry['Spoke']; pattern = entry['Pattern']
                ## Hub-Pattern ##
                if (hub,pattern) in edges['slim']:
                    edges['slim'][(hub,pattern)]['Sig'] = min(edges['slim'][(hub,pattern)]['Sig'], self.logSig(entry['Sig']))
                    edges['slim'][(hub,pattern)]['Evidence'] = rje.join(rje.sortUnique(rje.split(edges['slim'][(hub,pattern)]['Evidence'],'|')+[entry['Dataset']]),'|')
                else:
                    edges['slim'][(hub,pattern)] = rje.combineDict({},entry)
                    edges['slim'][(hub,pattern)]['Sig'] = self.logSig(entry['Sig'])
                    edges['slim'][(hub,pattern)]['Evidence'] = edges['slim'][(hub,pattern)].pop('Dataset')
                    edges['slim'][(hub,pattern)]['Start_Pos'] = 'NA'
                ## Pattern-Spoke ##
                if (pattern,spoke) in edges['occ']:
                    edges['occ'][(pattern,spoke)]['Sig'] = min(edges['occ'][(pattern,spoke)]['Sig'], self.logSig(entry['Sig']))
                    edges['occ'][(pattern,spoke)]['Evidence'] = rje.join(rje.sortUnique(rje.split(edges['occ'][(pattern,spoke)]['Evidence'],'|')+[entry['Dataset']]),'|')
                    edges['occ'][(pattern,spoke)]['Start_Pos'] = rje.join(rje.sortUnique(rje.split(edges['occ'][(pattern,spoke)]['Start_Pos'],',')+[entry['Start_Pos']]),',')
                else:                    
                    edges['occ'][(pattern,spoke)] = rje.combineDict({},entry)
                    edges['occ'][(pattern,spoke)]['Sig'] = self.logSig(entry['Sig'])
                    edges['occ'][(pattern,spoke)]['Evidence'] = edges['occ'][(pattern,spoke)].pop('Dataset')
                ## Random/Domain-Spoke ##
                etype = ''
                if hub[:4] in ['rupc','rseq']: etype = hub[:4]
                elif entry['Dataset'][-3:] == 'dom': etype = 'dmi'
                else: continue
                if (hub,spoke) in edges[etype]:
                    edges[etype][(hub,spoke)]['Evidence'] = rje.join(rje.sortUnique(rje.split(edges[etype][(hub,spoke)]['Evidence'],'|')+[entry['Pattern']]),'|')
                else:
                    edges[etype][(hub,spoke)] = rje.combineDict({},entry)
                    edges[etype][(hub,spoke)]['Evidence'] = edges[etype][(hub,spoke)].pop('Pattern')
                    edges[etype][(hub,spoke)]['Sig'] = -1
                #if 'P4HB' in (hub,spoke): self.deBug('\n%s :: (%s,%s) :: %s // %s' % (etype,hub,spoke,edges[etype][(hub,spoke)],entry))
            ## ~ [2b] ~ Add Cloud relationships as Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ex = 0.0; etot = len(sig.index('KCloud'))
            for cloud in sig.index('KCloud'):
                self.progLog('\r#EDGE','Making Cloud Linkage edges: %.2f%%' % (ex/etot)); ex += 100.0
                cpat = []
                for pattern in sig.indexDataList('KCloud',cloud,'Pattern',sortunique=False):
                    if pattern not in nodes: continue
                    for prevpat in cpat:
                        if (pattern,prevpat) in edges['cloud']: edges['cloud'][(pattern,prevpat)]['Evidence'] = rje.join(rje.sortUnique(rje.split(edges['cloud'][(pattern,prevpat)]['Evidence'],',')+[cloud]),',')
                        else: edges['cloud'][(pattern,prevpat)] = {'Sig':-len(rje.listIntersect(occ.indexDataList('Pattern',pattern,'Spoke'),occ.indexDataList('Pattern',pattern,'Spoke'))),'Evidence':cloud}
                    cpat.append(pattern)
            ## ~ [2c] ~ Domain-Protein Links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pfam = self.db().getTable('PFam')
            ex = 0.0; etot = len(nodes)
            for gene in rje.sortKeys(nodes):
                self.progLog('\r#EDGE','Making Domain-Protein edges: %.2f%%' % (ex/etot)); ex += 100.0
                try:
                    seq = self.data('DBXRef')[gene]['EnsLoci']
                    pdom = pfam.subset('Name',seq)
                    for entry in pdom.values():
                        domain = entry['Type'] 
                        if domain not in nodes: continue
                        if (gene,domain) in edges['pfam']:
                            edges['pfam'][(gene,domain)]['Sig'] -= l
                            edges['pfam'][(gene,domain)]['Evidence'] = rje.join(rje.split(edges['pfam'][(gene,domain)]['Evidence'],'|')+['%s:%s' % (entry['Start'],entry['Eval'])],'|')
                        else:
                            edges['pfam'][(gene,domain)] = rje.combineDict({},entry)
                            edges['pfam'][(gene,domain)]['Sig'] = -l
                            edges['pfam'][(gene,domain)]['Evidence'] = '%s:%s' % (entry['Start'],entry['Eval'])            
                except: continue
            self.printLog('\r#EDGE','Made %d PPI, Cloud, Occurrence and Domain-Protein edge types.' % len(edges))
            ## ~ [2x] ~ Debug ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #for etype in rje.sortKeys(edges): self.deBug('%s = *** %s ***' % (rje.sortKeys(edges[etype]),etype))
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            p = rje_ppi.PPI(self.log,['name=humsf09']+self.cmd_list)
            p.xgmmlToPPI(xgmml)
            p.db('Node').saveToFile('%s.node.tdt' % self.info['Basefile'])
            p.db('Edge').saveToFile('%s.edge.tdt' % self.info['Basefile'])
            fullppi = p.dict['PPI']
            p.info['Basefile'] = '%s.cloud' % self.info['Basefile']
            p.dict['PPI'] = rje_ppi.subGraph(fullppi,rje.sortKeys(sig.index('Pattern')),addv=False)
            x = p.ppiToXGMML(p.dict['PPI'],name=p.info['Basefile'])
            x.saveXGMML()
            #x#p.main()
            p.info['Basefile'] = '%s.hub' % self.info['Basefile']
            p.dict['PPI'] = rje_ppi.subGraph(fullppi,rje.sortKeys(sig.index('Hub')),addv=False)
            x = p.ppiToXGMML(p.dict['PPI'],name=p.info['Basefile'])
            x.saveXGMML()
            #x#p.main()
            p.info['Basefile'] = '%s.spoke' % self.info['Basefile']
            p.dict['PPI'] = rje_ppi.subGraph(fullppi,rje.sortKeys(occ.index('Spoke')),addv=False)
            x = p.ppiToXGMML(p.dict['PPI'],name=p.info['Basefile'])
            x.saveXGMML()
            #x#p.main()
            #xgmml.saveXGMML('humsf09.xgmml')
        except: self.errorLog('Problem with XGMML output')            
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
            if gene in self.db('Occ').index('Spoke'):
                href = 'http://bioware.soton.ac.uk/slimdb/humsf10_spoke/gene/%s.html' % gene
                title = alt = '%s Spoke Network results' % gene
                src = '../resources/humsf.png'
                html += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            html += ['</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,desc)]
            html += ['</table>','<!-- ~~~~ End %s Gene Summary details table ~~~~ -->' % gene,'']
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
                        gtab.append(('GO_%s' % gtype,rje.join(ghtml,' ~ '),self.titleText('tab','go')))
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s' % gene,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % gene,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % gene,'']
        except: self.errorLog('seqDetailsHTML Error')
        ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #print rje.join(html,'\n')
        return rje.join(html,'\n')
#########################################################################################################################
    def frontPage(self,genelist,domainlist,slimlist):    ### Generates HTML front page for analysis
        '''
        Generates HTML front page for analysis.
        >> genelist:list of analysed genes
        >> domainlist:list of analysed domains
        >> slimlist:list of returned slims
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = htmlHead('Human SLiMFinder Analysis',self.list['StyleSheets'],frontpage=True)
            html += '\n<table border="0" width="100%"><tr>\n<td width="60%"><h1>Human PPI SLiMFinder Analysis</h1>\n'
            html += '<br><p><i>This data has not yet been published and should not be used without permission.</i></p></td>\n'
            html += '<td width="40%" valign="top" align="right">\n'
            html += '<a href="http://bioware.ucd.ie"><img src="./resources/ucd.gif" height="100" alt="Bioware server" title="Bioware server"></a>\n'
            html += '<a href="http://www.personal.soton.ac.uk/re1u06/"><img src="./resources/SBS_100.png" height="100" alt="RJE Homepage" title="RJE Homepage"></a>\n'
            html += '</td></tr></table>\n\n'
            head_html = html[0:]
            main = self.db().getTable('Main')
            sig = self.db().getTable('Sig_real')
            sig_rand = self.db().getTable('Sig_rand')
            occ = self.db().getTable('Occ')
            occ_real = self.db().getTable('Occ_real')
            occ_rand = self.db().getTable('Occ_rand')
            ppi = self.db().getTable('PPI')
            dbxref = self.db().getTable('DBXRef')
            jtxt = ' ~ '
            ## ~ [1a] ~ Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            genelist.sort(); genedict = {}; genetab = []; genetitle = {}
            for x in '%s_' % string.ascii_uppercase:
                if x == '_':
                    genetitle[x] = 'Genes not beginning A-Y'
                    genedict[x] = ['<h2>Genes not beginning A-Y</h2><p>\n']
                else:
                    genetitle[x] = 'Genes beginning with %s' % x
                    genedict[x] = ['<h2>Genes beginning %s...</h2><p id="nonsig">\n' % x]
            for gene in genelist:
                x = gene[0]
                gtxt = gene
                if gene in sig.index('Hub'): gtxt += ' (%d|' % len(sig.index('Hub')[gene])
                elif gene in main.index('Hub'): gtxt += ' (0|'
                else: gtxt += ' (-|'
                if gene in occ_real.index('Spoke'): gtxt += '%d|' % len(occ_real.index('Spoke')[gene])
                elif gene in ppi.index('Spoke'): gtxt += '0|'
                else: gtxt += '-|'
                if gene in occ_rand.index('Spoke'): gtxt += '%d)' % len(occ_rand.index('Spoke')[gene])
                elif gene in dbxref.index('Gene'): gtxt += '0)'
                else: gtxt += '-)'
                if x not in string.ascii_uppercase: x = '_'
                #!# Update mouseover title text to explain numbers #!#
                sigtitle = '%s results page (Sig. results Hub|Spoke|Random)' % gene
                nosigtitle = '%s summary page' % gene
                if gene in sig.index('Hub') or gene in occ.index('Spoke'): genedict[x].append('<strong><a href="./gene/%s.html" target="_blank" title="%s">%s</a></strong>' % (gene,sigtitle,gtxt))
                else: genedict[x].append('<a href="./gene/%s.html" target="_blank" title="%s">%s</a>' % (gene,nosigtitle,gtxt))
            for x in '%s_' % string.ascii_uppercase:
                if len(genedict[x]) > 1: 
                    if x == '_': genetab.append(('*','%s\n</p>' % rje.join(genedict[x],jtxt),genetitle[x]))
                    else: genetab.append((x,'%s\n</p>' % rje.join(genedict[x],jtxt),genetitle[x]))
            ## ~ [1b] ~ Domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            domainlist.sort(); domaindict = {}; domaintab = []; domtitle = {}
            for x in '%s_' % string.ascii_uppercase:
                if x == '_':
                    domtitle[x] = 'Domains not beginning A-Z'
                    domaindict[x] = ['<h2>Domains not beginning A-Z</h2><p>\n']
                else:
                    domtitle[x] = 'Domains beginning with %s' % x
                    domaindict[x] = ['<h2>Domains beginning %s...</h2><p>\n' % x]
            for domain in domainlist:
                x = domain[0].upper()
                if x not in string.ascii_uppercase: x = '_'
                gtxt = domain
                if domain in sig.index('Hub'): gtxt += ' (%d)' % len(sig.index('Hub')[domain])
                elif domain in main.index('Hub'): gtxt += ' (0)'
                else: gtxt += ' (-)'
                domaindict[x].append('<a href="./domain/%s.html" target="_blank" title="%s">%s</a>' % (domain,self.titleText('link','domain'),gtxt))
            for x in '%s_' % string.ascii_uppercase:
                if len(domaindict[x]) > 1:
                    if x == '_': domaintab.append(('*','%s\n</p>' % rje.join(domaindict[x],jtxt),domtitle[x]))
                    else: domaintab.append((x,'%s\n</p>' % rje.join(domaindict[x],jtxt),domtitle[x]))
            ## ~ [1c] ~ SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slimlist.sort(); slimdict = {}; slimtab = []; slim2dict = {}; slimtitle = {}
            for x in '%s_^' % string.ascii_uppercase:
                if x == '_': 
                    slimtitle[x] = 'SLiMs beginning with ambiguity'
                    slimdict[x] = ['<h2>SLiMs beginning with ambiguity...</h2><p>\n']
                    slim2dict[x] = ['<h2>SLiMs containing ambiguity</h2><p>\n']
                elif x == '^': 
                    slimtitle[x] = 'N-terminal and C-terminal SLiMs'
                    slimdict[x] = ['<h2>N-terminal SLiMs</h2><p>\n']
                    slim2dict[x] = ['<h2>C-terminal SLiMs</h2><p>\n']
                else:
                    slimtitle[x] = 'SLiMs beginning with %s' % x
                    slimdict[x] = ['<h2>SLiMs beginning %s...</h2><p>\n' % x]
                    slim2dict[x] = ['<h2>SLiMs containing %s...</h2><p>\n' % x]
            for slim in slimlist:
                code = rje_slim.slimFromPattern(slim)   
                x = slim[0]
                gtxt = slim
                try: gtxt += ' (%d|' % len(sig.index('Pattern')[slim])
                except: gtxt += ' (0|'
                try: gtxt += '%d|' % len(occ.index('Pattern')[slim])
                except: gtxt += '0|'
                if slim in sig_rand.index('Pattern'): gtxt += '%d)' % len(sig_rand.index('Pattern')[slim])
                else: gtxt += '0)'
                if x not in '%s^' % string.ascii_uppercase: x = '_'
                slimdict[x].append('<a href="./slim/%s.html" target="_blank" title="%s">%s</a>' % (code,self.titleText('link','slim'),gtxt))
                for x in '%s_^' % string.ascii_uppercase:
                    if x == '_' and '[' in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank" title="%s">%s</a>' % (code,self.titleText('link','slim'),gtxt))
                    elif x == '^':
                        if '$' in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank" title="%s">%s</a>' % (code,self.titleText('link','slim'),gtxt))
                    elif x in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank" title="%s">%s</a>' % (code,self.titleText('link','slim'),gtxt))
            for x in '%s_^' % string.ascii_uppercase:
                if len(slimdict[x]) > 1:
                    if x == '_': slimtab.append(('*','%s\n</p>' % rje.join(slimdict[x],jtxt),slimtitle[x]))
                    elif x == '^': slimtab.append(('Terminal','%s\n</p>' % rje.join(slimdict[x]+slim2dict[x],jtxt),slimtitle[x]))
                    else: slimtab.append((x,'%s\n</p>' % rje.join(slimdict[x]+slim2dict[x],jtxt),slimtitle[x]))
            ## ~ [1d] ~ Random data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rseqlist = rje.sortKeys(self.db().getTable('Main_rseq').index('Dataset'))
            rupclist = rje.sortKeys(self.db().getTable('Main_rupc').index('Dataset'))
            pots = [(3,3),(4,4),(5,5),(6,10),(11,15),(16,25),(26,40),(41,60),(61,100),(101,200)]
            pdict = {}; rseqdict = {}; rupcdict = {}; plist = []; rantitle = {}
            for pot in pots:
                pkey = '%d-%d' % pot
                for x in range(pot[0],pot[1]+1): pdict['%d' % x] = pkey
                plist.append(pkey)
                if pot[0] == pot[1]:
                    rantitle[pkey] = 'Significant RAND datasets with %s UPC' % pot[0]
                    rseqdict[pkey] = ['<h2>Significant RSeq datasets with %s UPC...</h2><p>\n' % pot[0]]
                    rupcdict[pkey] = ['<h2>Significant RUPC datasets with %s UPC...</h2><p>\n' % pot[0]]
                else:
                    rantitle[pkey] = 'Significant RAND datasets with %s UPC' % pkey
                    rseqdict[pkey] = ['<h2>Significant RSeq datasets with %s UPC...</h2><p>\n' % pkey]
                    rupcdict[pkey] = ['<h2>Significant RUPC datasets with %s UPC...</h2><p>\n' % pkey]
            plist.append('201+')
            rantitle['201+'] = 'Significant RAND datasets with 201+ UPC'
            rseqdict['201+'] = ['<h2>Significant RSeq datasets with 201+ UPC...</h2><p>\n']
            rupcdict['201+'] = ['<h2>Significant RUPC datasets with 201+ UPC...</h2><p>\n']
            for dset in rseqlist:
                if dset not in sig_rand.index('Dataset'): continue
                try: pot = pdict[rje.split(dset,'-')[-1]]
                except: pot = '201+'
                rseqdict[pot].append('<a href="./rseq/%s.html" target="_blank">%s (%d)</a>' % (dset,dset,len(sig_rand.index('Dataset')[dset])))
            for dset in rupclist:
                if dset not in sig_rand.index('Dataset'): continue
                try: pot = pdict[rje.split(dset,'-')[-1]]
                except: pot = '201+'
                rupcdict[pot].append('<a href="./rupc/%s.html" target="_blank">%s (%d)</a>' % (dset,dset,len(sig_rand.index('Dataset')[dset])))
            rseqtab = []
            rupctab = []
            for pot in plist:
                rseqtab.append((pot,'%s\n</p>' % rje.join(rseqdict[pot],jtxt),rje.replace(rantitle[pot],'RAND','random sequence (RSeq)')))
                rupctab.append((pot,'%s\n</p>' % rje.join(rupcdict[pot],jtxt),rje.replace(rantitle[pot],'RAND','random UPC (RUPC)')))
            ## ~ [1e] ~ Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datahtml = ['<!-- ~~~~~~~~ TARRED and Zipped Data Archives ~~~~~~~ -->',
                        '<h2>SLiMFinder Analysis Data Archive</h2>','','<dl>','',
                        ## Main ##
                        '<dt><a href="./data/humsf09_maindata.tar.gz" title="Download gzipped tar file"><strong>Main Data</strong></a>',
                        '<dd><p><i>Main results tables, with a few additional columns.</i></p><ul type="square">',
                        '<li> Summary results for all runs, including random data = humsf09_slimfinder.full.csv</li>',
                        '<li> All real data occurrences = humsf09_slimfinder.occ.csv</li>',
                        '<li> Pairwise interaction data (May 2009), including evidence = humsf09.pairwise_ppi.090505.tdt',
                        '<li> Domain-Protein interaction data (May 2009), including evidence = humsf09.domain_ppi.090505.tdt',
                        '<li> Sequence ID mapping file (May 2009) = humsf09.genemap.0905050.tdt',
                        '<li> Human sequence file = ens_HUMAN.loci.fas (May 2009)</li>','</ul><br>',
                        ## Sig ##
                        '<dt><a href="./data/humsf09_sig.tar.gz" title="Download gzipped tar file"><strong>Significant Datasets</strong></a>',
                        '<dd><p><i>Additional results files (*.tdt, *.txt, *.csv *.upc & *.fas) for Significant runs:</i></p><ul type="square">',
                        '<li> cloud.txt', '<li> dis.tdt', '<li> mapping.fas', '<li> maskaln.fas', '<li> motifaln.fas', '<li> upc</ul><br>',
                        # !! NB. The original *.occ.csv had a bug that incorrectly positions the motifs !! #
                        ## UniFake ##
                        '<dt><a href="./data/humsf09_unifake.tar.gz" title="Download gzipped tar file"><strong>UniFake</strong></a>',
                        '<dd><p><i>Uniprot-format DAT file of input sequences (May 2009) and PFam domain predictions for all human sequences (May 2009).</i></p>',
                        ## Non-Sig ##
                        '<dt><a href="./data/humsf09_nosig.tar.gz" title="Download gzipped tar file"><strong>Non-Significant Datasets</strong></a>',
                        '<dd><p><i>Additional results files (*.dis.tdt, *.upc & *.masked.fas) for Non-Significant runs.</i></p>',
                        ## GOPHER ##
                        '<dt><a href="./data/humsf09_gopheraln.tar.gz" title="Download gzipped tar file"><strong>GOPHER Orthologue Alignments</strong></a>',
                        '<dd><p><i>*.orthaln.fas for each human protein (May 2009).</i></p>',
                        ## PPIFas ##
                        '<dt><a href="./data/humsf09_ppifas.tar.gz" title="Download gzipped tar file"><strong>PPI Fasta files (input)</strong></a>',
                        '<dd><p><i>Directories containing input fasta files of PPI datasets (May 2009 data)</i></p>',
                        '</dl>','','<!-- ~~~~~~~~ End of Data Archives ~~~~~~~ -->','']
            ## ~ [1f] ~ True Positives Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tphtml = ['<!-- ~~~~~~~~ True Positives Table ~~~~~~~~~~~~~~ --!>',
                      '<h2>True Positive SLiMFinder Results</h2>','']
            tpfile = 'humsf09.tp_table_for_paper.tdt'
            if os.path.exists(tpfile):
                tpdata = rje.dataDict(self,tpfile,['Motif','Hub'],'all')
                tphtml += ['<table width=1000 border=%d><tr>' % self.stat['Border'],
                           '<th title="Known SLiM rediscovered">Motif</th>',
                           '<th title="Hub gene or domain returning known SLiM">Hub</th>',
                           '<th title="Top matching pattern for full PPI dataset">ppi</th>',
                           '<th title="Top matching pattern for Yeast-2-Hybrid PPI dataset">y2h</th>',
                           '<th title="Top matching pattern for binary-enriched PPI dataset">bin</th>',
                           '<th title="Top matching pattern for complex-enriched PPI dataset">com</th>',
                           '</tr>']
                for tpkey in rje.sortKeys(tpdata):
                    entry = tpdata[tpkey]
                    #self.deBug(entry)
                    tphtml.append('<tr>%s' % self.knownTableLink(entry['Motif']))
                    if entry['DType'] == 'gene': tphtml.append('<td>%s</td>' % self.geneLink(entry['Hub'],frontpage=True))
                    else: tphtml.append('<td>%s</td>' % domainLink(entry['Hub'],frontpage=True))
                    for ptype in ['ppi','y2h','bin','com']:
                        if entry[ptype] == 'n/a': tphtml.append('<td><i>n/a</i></td>')
                        elif entry[ptype]:
                            (pattern,sig) = rje.split(entry[ptype])
                            tphtml.append('<td>%s %s</td>' % (slimLink(pattern,frontpage=True),sig))
                        else: tphtml.append('<td></td>')
                    tphtml.append('</tr>')
                tphtml += ['</table>']
            else: tphtml += ['<i>TP Table not found.</i>']       
            ## ~ [1g] ~ Annotated SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            kdb = self.db().getTable('Knowns'); tdb = self.db().getTable('TP');
            ppi = self.db().getTable('PPI')
            ppi.index('Hub'); ppi.index('Spoke')
            skipknown = ['!','!!','CAAXbox','PY','THW','leak?','AKAP10','LIG_AP2?','Nterm','Cterm']
            annhtml = ['<!-- ~~~~~~~~ Annotated Motifs Table ~~~~~~~~~~~~~~ --!>',
                      '<h2>Annotated SLiMFinder Results</h2>','','<table border=%d><tr>' % self.stat['Border'],
                       '<th title="Known SLiM rediscovered" width=12%>Motif</th>',
                       '<th title="Returned patterns matching known SLiM" width=20%>Patterns</th>',
                       '<th title="Hub genes returning known SLiM" width=17%>Gene Hubs</th>',
                       '<th title="Hub domains returning known SLiM" width=17%>Domain Hubs</th>',
                       '<th title="Random sequence (RSeq) datasets returning known SLiM" width=17%>RSeq</th>',
                       '<th title="Random UPC (RUPC) datasets returning known SLiM" width=17%>RUPC</th>',
                       '</tr>']
            anntab = []; kgenes = []; kdoms = []; krseq = []; krupc = []; kslim = []; klist = []
            for known in rje.sortKeys(kdb.index('Known')):
                if known in skipknown or '?' in known or '!' in known: continue
                if known[:3] in ['LIG','MOD']: klist.append(known)
            for known in rje.sortKeys(kdb.index('Known')):
                if known in skipknown or '?' in known or '!' in known: continue
                if known[:3] not in ['LIG','MOD']: klist.append(known)
            i = 5
            khtml = {'Patterns':['<th title="Returned patterns matching known SLiM" width=80%>Patterns</th>'],
                     'Genes':['<th title="Hub genes returning known SLiM" width=80%>Gene Hubs</th>'],
                     'Domains':['<th title="Hub domains returning known SLiM" width=80%>Domain Hubs</th>'],
                     'RSeq':['<th title="Random sequence (RSeq) datasets returning known SLiM" width=80%>RSeq</th>'],
                     'RUPC':['<th title="Random UPC (RUPC) datasets returning known SLiM" width=80%>RUPC</th>']}
            ktitle = {}
            for k in khtml:
                ktitle[k] = rje.split(khtml[k][0],'"')[1]
                khtml[k] = [rje.replace(rje.join(annhtml[1:i],'\n'),'Results','Results (%s)' % k)] + khtml[k] + ['</tr>']
            for known in klist:
                annhtml += ['<tr>%s' % self.knownTableLink(known)]
                for k in khtml: khtml[k] += [annhtml[-1]] 
                kpat = []
                for entry in kdb.indexEntries('Known',known):
                    kpat.append(slimLink(entry['Pattern'],frontpage=True))
                    if entry['Pattern'] not in kslim: kslim.append(entry['Pattern'])
                kpat.sort()
                annhtml += ['<td>%s</td>' % rje.join(kpat,jtxt)]
                khtml['Patterns'].append('%s</tr>' % annhtml[-1])
                kpat = []
                for entry in tdb.indexEntries('Known',known):
                    if entry['Hub'][:4] in ['rseq','rupc']: continue
                    if entry['Hub'] in ppi.index('Hub') or entry['Hub'] in ppi.index('Spoke'): kpat.append(self.geneLink(entry['Hub'],frontpage=True))
                    if entry['Hub'] not in kgenes: kgenes.append(entry['Hub'])
                kpat.sort()
                annhtml += ['<td>%s</td>' % rje.join(kpat,jtxt)]
                khtml['Genes'].append('%s</tr>' % annhtml[-1])
                kpat = []
                for entry in tdb.indexEntries('Known',known):
                    if entry['Hub'][:4] in ['rseq','rupc']: continue
                    if entry['Hub'] not in ppi.index('Hub') and entry['Hub'] not in ppi.index('Spoke'): kpat.append(domainLink(entry['Hub'],frontpage=True))
                    if entry['Hub'] not in kdoms: kdoms.append(entry['Hub'])
                kpat.sort()
                annhtml += ['<td>%s</td>' % rje.join(kpat,jtxt)]
                khtml['Domains'].append('%s</tr>' % annhtml[-1])
                kpat = []
                for entry in tdb.indexEntries('Known',known):
                    if entry['Hub'][:4] != 'rseq': continue
                    if entry['Hub'] in main.index('Hub'):
                        dataset = main.indexEntries('Hub',entry['Hub'])[0]['Dataset']
                        kpat.append(randLink(dataset,frontpage=True))
                        if dataset not in krseq: krseq.append(dataset)
                    else:
                        continue    # Identifiers from old SF run?
                        kpat.append(entry['Hub'])
                        if entry['Hub'] not in krseq: krseq.append(entry['Hub'])
                kpat.sort()
                annhtml += ['<td>%s</td>' % rje.join(kpat,jtxt)]
                khtml['RSeq'].append('%s</tr>' % annhtml[-1])
                kpat = []
                for entry in tdb.indexEntries('Known',known):
                    if entry['Hub'][:4] != 'rupc': continue
                    if entry['Hub'] in main.index('Hub'):
                        dataset = main.indexEntries('Hub',entry['Hub'])[0]['Dataset']
                        kpat.append(randLink(dataset,frontpage=True))
                        if dataset not in krupc: krupc.append(dataset)
                    else:
                        continue    # Identifiers from old SF run?
                        kpat.append(entry['Hub'])
                        if entry['Hub'] not in krupc: krupc.append(entry['Hub'])
                kpat.sort()
                annhtml += ['<td>%s</td></tr>' % rje.join(kpat,jtxt)]
                khtml['RUPC'].append('%s</tr>' % annhtml[-1])
            for k in khtml: khtml[k] += ['<tr align=center><th title="Total counts of annotated SLiMs">Total:</th>']
            khtml['Patterns'] += ['<th title="Total number of different Patterns annotated as known">%d</th>' % len(kslim)]
            khtml['Genes'] += ['<th title="Total number of different Gene Hubs returning significant motifs annotated as known">%d</th>' % len(kgenes)]
            khtml['Domains'] += ['<th title="Total number of different Domain Hubs returning significant motifs annotated as known">%d</th>' % len(kdoms)]
            khtml['RSeq'] += ['<th title="Total number of different RSeq Hubs returning significant motifs annotated as known">%d (%s dom)</th>' % (len(krseq),rje.count(rje.join(krseq),'dom'))]
            khtml['RUPC'] += ['<th title="Total number of different RUPC Hubs returning significant motifs annotated as known">%d (%d dom)</th>' % (len(krupc),rje.count(rje.join(krupc),'dom'))]
            for k in ['Patterns','Genes','Domains','RSeq','RUPC']: 
                khtml[k] += ['</tr></table>']
                anntab.append((k,rje.join(khtml[k],'\n'),ktitle[k]))
                # Totals #
            annhtml += ['<tr align=center><th title="Total counts of annotated SLiMs">Total:</th>',
                        '<th title="Total number of different Patterns annotated as known">%d</th>' % len(kslim),
                        '<th title="Total number of different Gene Hubs returning significant motifs annotated as known">%d</th>' % len(kgenes),
                        '<th title="Total number of different Domain Hubs returning significant motifs annotated as known">%d</th>' % len(kdoms),
                        '<th title="Total number of different RSeq Hubs returning significant motifs annotated as known">%d (%s dom)</th>' % (len(krseq),rje.count(rje.join(krseq),'dom')),
                        '<th title="Total number of different RUPC Hubs returning significant motifs annotated as known">%d (%d dom)</th>' % (len(krupc),rje.count(rje.join(krupc),'dom')),
                        '</tr>']
            annhtml += ['</table>']
            annhtml = annhtml[:1] + [tabberHTML('Known',anntab,level=1)]
            ## ~ [1h] ~ SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            unslimlist = self.loadFromFile('humsf09.unannotated.txt',chomplines=True)
            if '-' in unslimlist: unslimlist.remove('-')
            ndb = self.db().copyTable('Main','Novel',replace=True)
            ndb.dropEntriesDirect('Pattern',unslimlist,inverse=True,log=True)
            ndb.dropEntriesDirect('Class',['Generic','OT','TP','RAND'],inverse=False,log=True)
            ndb.dataFormat({'Sig':'num'})
            ndb.rankField('Sig','SigRank',absolute=True,lowest=True)
            for entry in ndb.entries(): entry['SigRank'] = rje.preZero(int(entry['SigRank']),ndb.entryNum())
            ndb.newKey(['SigRank','Dataset','Pattern'],startfields=False)
            for entry in ndb.entries():
                analysis = rje.split(entry['Analysis'],'.')
                if analysis[0] == 'real': entry['Analysis'] = analysis[1]
                elif analysis[1][-3:] == 'dom': entry['Analysis'] = '%sdom' % analysis[0]
                else: entry['Analysis'] = analysis[0]
            unsplit = self.db().splitTable(ndb,'Analysis',asdict=True)
            unslimtab = []
            for type in ['ppi','y2h','bin','com','ppidom','y2hdom','bindom','comdom','rseq','rseqdom','rupc','rupcdom']:
                unslimhtml = '<h2>Potential Novel %s SLiMs</h2>\n' % type
                if type in unsplit: unslimhtml += self.resultTableHTML(unsplit[type].fields(),unsplit[type].data(),True,drop=self.list['DropFields']+['SigRank','Analysis','Class','Match'])
                else: unslimhtml += '<i>No novel SLiMs</i>'
                unslimtab.append((type,rje.replace(unslimhtml,'href="../','target="_blank" href="./'),'Potential Novel %s SLiMs' % type))
            if not 'this_is_the_old_code':            
                slimdict = {}; unslimtab = []; slim2dict = {}
                for x in '%s_^' % string.ascii_uppercase:
                    if x == '_': 
                        slimdict[x] = ['<h2>SLiMs beginning with ambiguity...</h2><p>\n']
                        slim2dict[x] = ['<h2>SLiMs containing ambiguity</h2><p>\n']
                    elif x == '^': 
                        slimdict[x] = ['<h2>N-terminal SLiMs</h2><p>\n']
                        slim2dict[x] = ['<h2>C-terminal SLiMs</h2><p>\n']
                    else:
                        slimdict[x] = ['<h2>SLiMs beginning %s...</h2><p>\n' % x]
                        slim2dict[x] = ['<h2>SLiMs containing %s...</h2><p>\n' % x]
                for slim in unslimlist:
                    if not slim or slim == '-': continue
                    code = rje_slim.slimFromPattern(slim)   
                    x = slim[0]
                    gtxt = slim
                    try: gtxt += ' (%d|' % len(sig.index('Pattern')[slim])
                    except: gtxt += ' (0|'
                    try: gtxt += '%d|' % len(occ.index('Pattern')[slim])
                    except: gtxt += '0|'
                    if slim in sig_rand.index('Pattern'): gtxt += '%d)' % len(sig_rand.index('Pattern')[slim])
                    else: gtxt += '0)'
                    if x not in '%s^' % string.ascii_uppercase: x = '_'
                    slimdict[x].append('<a href="./slim/%s.html" target="_blank">%s</a>' % (code,gtxt))
                    for x in '%s_^' % string.ascii_uppercase:
                        if x == '_' and '[' in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank">%s</a>' % (code,gtxt))
                        elif x == '^' and '$' in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank">%s</a>' % (code,gtxt))
                        elif x in slim: slim2dict[x].append('<a href="./slim/%s.html" target="_blank">%s</a>' % (code,gtxt))
                for x in '%s_^' % string.ascii_uppercase:
                    if len(slimdict[x]) > 1:
                        if x == '_': unslimtab.append(('*','%s\n\</p>' % rje.join(slimdict[x],jtxt)))
                        elif x == '^': unslimtab.append(('Terminal','%s\n</p>' % rje.join(slimdict[x]+slim2dict[x],jtxt)))
                        else: unslimtab.append(('%s' % x,'%s\n</p>' % rje.join(slimdict[x]+slim2dict[x],jtxt)))
            ## ~ [1i] ~ GO Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb = self.db().getTable('GO'); gotab = []
            got = {'bp':'Biological Process','cc':'Cellular Component','mf':'Molecular Function'}
            for type in ['bp','cc','mf']:
                #gohtml = ['<h2>%s GO terms</h2>' % got[type]]
                gosub = {}
                for x in string.ascii_uppercase: gosub[x] = ['<h3>%s GO terms starting %s</h3>' % (got[type],x)]
                for go in gdb.dataList(gdb.indexEntries('GO_Type',type),'GO_Desc'):
                    if not go: continue
                    i = 0
                    while go[i].upper() not in gosub: i += 1
                    x = go[i].upper()
                    #gohtml.append(goLink(go,gdb.dataList(gdb.indexEntries('GO_Desc',go),'GO_ID')[0],frontpage=True))
                    gosub[x].append(goLink(go,gdb.dataList(gdb.indexEntries('GO_Desc',go),'GO_ID')[0],frontpage=True))
                    #!# Add numbers of genes, domains & slims #!#
                gohtml = []
                for x in string.ascii_uppercase: gohtml.append((x,rje.join(gosub[x],'\n ~ '),'%s GO terms starting %s' % (got[type],x)))
                #gotab.append(('GO_%s' % type.upper(),rje.join(gohtml,'\n ~ ')))
                gotab.append(('GO_%s' % type.upper(),tabberHTML('GO_%s' % type.upper(),gohtml,level=1),'%s GO terms' % got[type]))
            ### ~ [2] ~ Main code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            moretext = '<p>Click <a href="./full.htm" target="_blank">here</a> to access data by gene/domain/pattern.</p>\n'
            moretext += '<p>Click to view overlapping UniProt features for\n'
            moretext += '<a href="./slimft_tp.htm" target="_blank">TP</a>, '
            moretext += '<a href="./slimft_known.htm" target="_blank">Annotated</a> and '
            moretext += '<a href="./slimft_novel.htm" target="_blank">Novel</a> '
            moretext += 'motifs.</p>\n\n'
            moretext += '<table><tr><td><a href="http://bioware.soton.ac.uk/slimdb/humsf10_spoke/">\n'
            moretext += '<img alt="Significant Spoke cluster analysis" src="./resources/humsf.png" align="MIDDLE" border="0" height="50" title="Significant Spoke cluster analysis">'
            moretext += 'Significant Spoke cluster analysis</a></td></tr>\n'
            moretext += '<tr><td><a href="http://bioware.soton.ac.uk/slimdb/humsf10_slim/">\n'
            moretext += '<img alt="Significant SLiM cluster analysis" src="./resources/humsf.png" align="MIDDLE" border="0" height="50" title="Significant SLiM cluster analysis">'
            moretext += 'Significant SLiM cluster analysis</a></td></tr></table>\n\n'
            tablist = [('Data',rje.join(datahtml,'\n'),self.titleText('front','data')),
                       ('TP',rje.join(tphtml,'\n'),self.titleText('front','tp')),
                       ('Annotated',rje.join(annhtml,'\n'),self.titleText('front','annotated')),
                       ('Novel',tabberHTML('Novel',unslimtab,level=1),self.titleText('front','novel')),
                       ('More...',moretext,self.titleText('front','more'))]
            fulltab = [('Genes',tabberHTML('Gene',genetab,level=1),self.titleText('front','gene')),
                       ('Domains',tabberHTML('Domain',domaintab,level=1),self.titleText('front','domain')),
                       ('SLiMs',tabberHTML('SLiM',slimtab,level=1),self.titleText('front','slim')),   #tabberHTML('Data',datatab,level=1))]
                       ('RSeq',tabberHTML('RSeq',rseqtab,level=1),self.titleText('front','rseq')),
                       ('RUPC',tabberHTML('RUPC',rupctab,level=1),self.titleText('front','rupc'))] + gotab
                       #('GO',tabberHTML('GO',gotab,level=1),self.titleText('front','go'))]
            ## ~ [2a] ~ Front page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add a little text ? #!#
            html += tabberHTML('Main',tablist,level=0)
            html += htmlTail()
            ## ~ [2b] ~ Full page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fullhtml = head_html
            fullhtml += '<p>Full results data can be accessed via the links on these pages.\n'
            fullhtml += 'Click <a href="./index.htm">here</a> to return to the summary data and download page.</p>\n\n'
            fullhtml += tabberHTML('Full',fulltab,level=0)
            fullhtml += htmlTail()
        except: self.errorLog('frontPage Error')
        ### ~ [3] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return {'front':html,'full':fullhtml}
#########################################################################################################################
    def slimFTPage(self,genelist,domainlist,slimlist):   ### Generates HTML front page of overlapping SLiM Features
        '''
        Generates HTML front page of overlapping SLiM Features.
        >> genelist:list of analysed genes
        >> domainlist:list of analysed domains
        >> slimlist:list of returned slims
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fthfile = self.info['HTMLPath'] + 'slimft.htm'
            if self.checkHTML(fthfile) and not self.opt['Force']: return
            self.progLog('\r#FTHTM','Making SLiMFT HTML Page: 0.00%')
            html = htmlHead('Human SLiMFinder Analysis',self.list['StyleSheets'],frontpage=True)
            html += '\n<table border="0" width="100%"><tr>\n<td width="60%"><h1>Human PPI SLiMFinder Analysis</h1>\n'
            html += '<br><p><i>This data has not yet been published and should not be used without permission.</i></p></td>\n'
            html += '<td width="40%" valign="top" align="right">\n'
            html += '<a href="http://bioware.ucd.ie"><img src="./resources/ucd.gif" height="100" alt="Bioware server" title="Bioware server"></a>\n'
            html += '<a href="http://www.personal.soton.ac.uk/re1u06/"><img src="./resources/SBS_100.png" height="100" alt="RJE Homepage" title="RJE Homepage"></a>\n'
            html += '</td></tr></table>\n\n'
            head_html = html[0:]
            main = self.db().getTable('Main')
            sig = self.db().getTable('Sig_real')
            sig_rand = self.db().getTable('Sig_rand')
            occ = self.db().getTable('Occ')
            occ_real = self.db().getTable('Occ_real')
            occ_rand = self.db().getTable('Occ_rand')
            ppi = self.db().getTable('PPI')
            dbxref = self.db().getTable('DBXRef')
            jtxt = ' ~ '
            hx = 0.0
            ### ~ [2] ~ HTML Feature listings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tphtml = ['<!-- ~~~~~~~~ True Positives ~~~~~~~~~~~~~~ --!>',
                      '<h2>UniFake features overlapping SLiM Occurrences (Patterns matching TP SLiMs)</h2>','']
            tpfile = 'humsf09.tp_table_for_paper.tdt'
            tpslims = []; ttot= 0
            if os.path.exists(tpfile):
                tpdata = rje.dataDict(self,tpfile,['Motif','Hub'],'all')
                for tpkey in rje.sortKeys(tpdata): 
                    entry = tpdata[tpkey]
                    for ptype in ['ppi','y2h','bin','com']:
                        if entry[ptype] == 'n/a': continue
                        elif entry[ptype]:
                            (pattern,sig) = rje.split(entry[ptype])
                            if pattern not in tpslims: tpslims.append(pattern)
                tpslims.sort(); tx = 0.0; ttot = len(tpslims)
                for slim in tpslims:
                    fthtml = self.frontPageFT(slim)
                    interest = rje.matchExp('FT   REGION.+(\Snteraction)',fthtml)
                    interest = interest or rje.matchExp('FT   VARIANT.+(\Sutation)',fthtml)
                    if interest: tphtml += [fthtml]
                    self.progLog('\r#FTHTM','Making SLiMFT HTML Page: %.2f%%' % (tx/ttot)); tx += 33.333
            else: tphtml += ['<i>TP Table not found.</i>']       
            ## ~ [1g] ~ Annotated SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            kdb = self.db().getTable('Knowns'); tdb = self.db().getTable('TP');
            ppi = self.db().getTable('PPI')
            ppi.index('Hub'); ppi.index('Spoke')
            skipknown = ['!','!!','CAAXbox','PY','THW','leak?','AKAP10','LIG_AP2?','Nterm','Cterm']
            annhtml = ['<!-- ~~~~~~~~ Annotated SLiMs ~~~~~~~~~~~~~~ --!>',
                      '<h2>UniFake features overlapping SLiM Occurrences (Annotated SLiMs)</h2>','']
            kslim = []
            for known in rje.sortKeys(kdb.index('Known')):
                if known in skipknown or '?' in known or '!' in known: continue
                for entry in kdb.indexEntries('Known',known):
                    if entry['Pattern'] not in tpslims + kslim: kslim.append(entry['Pattern'])
            kslim.sort(); kx = 0.0; ktot = len(kslim)
            for slim in kslim:
                fthtml = self.frontPageFT(slim)
                interest = rje.matchExp('FT   REGION.+(\Snteraction)',fthtml)
                interest = interest or rje.matchExp('FT   VARIANT.+(\Sutation)',fthtml)
                if interest: annhtml += [fthtml]
                self.progLog('\r#FTHTM','Making SLiMFT HTML Page: %.2f%%' % (33.333 + kx/ktot)); kx += 33.333
            ## ~ [1h] ~ Novel SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            unslimlist = self.loadFromFile('humsf09.unannotated.txt',chomplines=True)
            if '-' in unslimlist: unslimlist.remove('-')
            uhtml = ['<!-- ~~~~~~~~ Unknown/Novel SLiMs ~~~~~~~~~~~~~~ --!>',
                     '<h2>UniFake features overlapping SLiM Occurrences (Candidate Novel SLiMs)</h2>','']
            unslimlist.sort(); utot = ktot + ttot + len(unslimlist); ux = (ktot+ttot) * 100.0
            for slim in unslimlist:
                fthtml = self.frontPageFT(slim)
                interest = rje.matchExp('FT   REGION.+(\Snteraction)',fthtml)
                interest = interest or rje.matchExp('FT   VARIANT.+(\Sutation)',fthtml)
                if interest: uhtml += [fthtml]
                self.progLog('\r#FTHTM','Making SLiMFT HTML Page: %.2f%%' % (ux/utot)); ux += 100.0
            ### ~ [2] ~ Main code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = [('TP',rje.join(tphtml,'\n'),self.titleText('front','tp')),
                       ('Annotated',rje.join(annhtml,'\n'),self.titleText('front','annotated')),
                       ('Novel',rje.join(uhtml,'\n'),self.titleText('front','novel'))]
            tablist.append(('More...',
                            '<p>Click <a href="./index.htm">here</a> to return to the summary data and download page.</p>\n<p>Click <a href="./full.htm">here</a> to access data by gene/domain/pattern.</p>',
                            self.titleText('front','more')))
            ## ~ [2a] ~ Front page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add a little text ? #!#
            open(fthfile,'w').write(html + tabberHTML('SLiMFT',tablist,level=0) + htmlTail())
            self.printLog('\r#FTHTM','SLiMFT HTML Page data output to %s' % fthfile)
            fthfile = self.info['HTMLPath'] + 'slimft_tp.htm'
            open(fthfile,'w').write(rje.join([html]+tphtml+[htmlTail()],'\n'))
            self.printLog('\r#FTHTM','SLiMFT HTML Page data output to %s' % fthfile)
            fthfile = self.info['HTMLPath'] + 'slimft_known.htm'
            open(fthfile,'w').write(rje.join([html]+annhtml+[htmlTail()],'\n'))
            self.printLog('\r#FTHTM','SLiMFT HTML Page data output to %s' % fthfile)
            fthfile = self.info['HTMLPath'] + 'slimft_novel.htm'
            open(fthfile,'w').write(rje.join([html]+uhtml+[htmlTail()],'\n'))
            self.printLog('\r#FTHTM','SLiMFT HTML Page data output to %s' % fthfile)
        except: self.errorLog('slimFTPage Error')
#########################################################################################################################
    def knownTableLink(self,known):  ### Returns linked and titles Annotated motif name
        '''Returns linked and titles Annotated motif name.'''
        adb = self.db().getTable('Annotation')
        if known not in adb.data(): return '<td>%s</td>' % known
        entry = adb.data()[known]
        if entry['Link'] not in ['-','','.']: return '<td title="%s"><a href="%s" target="_blank">%s</a></td>' % (entry['Description'],entry['Link'],known)
        else: return '<td title="%s">%s</td>' % (entry['Description'],known)
#########################################################################################################################
    def tableClassLink(self,data):  ### Returns linked and titles Annotated motif name
        '''Returns linked and titles Annotated motif name.'''
        if 'Class' not in data or 'Annotation' not in data: return ''
        adb = self.db().getTable('Annotation')
        title = []
        for known in rje.split(data['Annotation'],','):
            if known not in adb.data(): title.append(known)
            else:
                entry = adb.data()[known]
                title.append('%s: %s' % (known,entry['Link']))
        return '<b title="%s">%s</b>' % (rje.join(title,'; '),data['Class'])
#########################################################################################################################
    def svgCode(self,svglink,title):    ### Returns HTML Code for SVG file embedding.
        '''Returns HTML Code for SVG file embedding.'''
        try:
            svgfile = self.info['HTMLPath'] + rje.join(rje.split(svglink,'/')[1:],'/')
            svghtm = '<p title="%s">\n' % (title)
            try: (width,height) = rje.matchExp('<svg width="(\d+)" height="(\d+)" version="1.1"',open(svgfile,'r').read())
            except: (width,height) = (1600,1600)
            svghtm += '<embed src="%s" width="%s" height="%s" type="image/svg+xml"' % (svglink,width,height)
            svghtm += ' pluginspage="http://www.adobe.com/svg/viewer/install/" /></p>'
            return svghtm
        except: self.errorLog(rje_zen.Zen().wisdom()); return '<i>SVG code error!</i>'
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
            tablist.append(('Interactors',self.interactorTableHTML(gene),'Table of %s interactors and SLiM predictions' % gene))
            ## ~ [1a2] ~ Domain Interactors tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('Domains',self.domainTableHTML(gene),'%s domains and table of domains interacting with %s' % (gene,gene)))
            ## ~ [1b] ~ Interactome tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            inthtml = []
            for ptype in ['ppi','y2h','bin','com']:
                dataset = '%s.%s' % (gene,ptype)
                if dataset not in self.db().getTable('Main').index('Dataset'): continue
                #intpng = rje.makePath('%sinteractome_png/%s.interactome.png' % (self.info['HTMLPath'],dataset),wholepath=True)
                if self.opt['SVG']:
                    svglink = '../interactome_svg/%s.interactome.svg' % dataset
                    title = 'Graphical representation of %s interactome. Larger interactomes may not be clear. See tables for details.' % dataset
                    ctab = '%s\n' % self.svgCode(svglink,title)
                    inthtml += ['<h2>%s Interactome</h2>' % dataset,ctab]
                else:
                    intpng = '../interactome_png/%s.interactome.png' % dataset
                    alt = '%s interactome image not yet made' % dataset
                    title = 'Graphical representation of %s interactome. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % dataset
                    inthtml += ['<h2>%s Interactome</h2>' % dataset,'<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)]
            if inthtml: tablist.append(('Interactome',rje.join(inthtml,'\n'),'Interactome graphics'))
            ## ~ [1c] ~ UPC tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('UPC',self.UPCHTML(gene),'UPC clustering of spoke proteins'))
            ## ~ [1c] ~ Hub Results Summary tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hubhtml = ['<h2>%s Hub results summary</h2>' % gene]
            maindb = self.db().getTable('Main_real')
            occdb = self.db().getTable('Occ_real')
            if gene in maindb.index('Hub'):
                headers = maindb.fields()
                data = maindb.subset('Hub',gene)
                hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Match']))
                if gene in occdb.index('Hub'):
                    hubhtml += ['<h2>SLiM occurrences for %s hub results</h2>' % gene]
                    headers = occdb.fields()
                    data = occdb.subset('Hub',gene)
                    hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']))
            else: hubhtml.append('<i>No Hub Datasets.</i>')
            tablist.insert(0,('Hub',rje.join(hubhtml,'\n'),'%s Hub results summary tables' % gene))
            #-- Summary graphics?
            ## ~ [1d] ~ Spoke Summary Tab (SLiM Occurrences in "Gene") ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sphtml = ['<h2>%s Spoke results summary</h2>' % gene]
            if gene in occdb.index('Spoke'):
                headers = occdb.fields()
                data = occdb.subset('Spoke',gene)
                sphtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Spoke','Desc']))
            elif gene in self.data('DBXRef'): sphtml.append('<i>No spoke results.</i>')
            else: sphtml.append('<i>No sequence mapping.</i>')
            tablist.insert(1,('Spoke',rje.join(sphtml,'\n'),'%s Spoke results summary tables' % gene))
            ## ~ [1e] ~ UniFake tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if gene in self.data('DBXRef'):
                ufile = rje.makePath(self.info['HTMLPath']+'unifake') + '%s.html' % seqid
                if os.path.exists(ufile):
                    uhtml = open(ufile,'r').read()
                    uhtml = uhtml[uhtml.find('<H2>UniProt Format'):]
                    uhtml = uhtml[:uhtml.find('</PRE>')] + '</PRE>\n'
                else: uhtml = '<i>UniFake file to be generated.</i>'
            else: uhtml = '<i>No sequence mapping.</i>'
            tablist.append(('UniFake',uhtml,'Modified %s UniProt entry' % gene))
            ## ~ [1f] ~ RSeq Results tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if gene in self.data('DBXRef'): 
                rseqhtml = ['<h2>%s RSeq results summary</h2>' % gene]
                maindb = self.db().getTable('Main_rseq')
                occdb = self.db().getTable('Occ_rseq')
                if gene in occdb.index('Spoke'):
                    headers = occdb.fields()
                    data = occdb.subset('Spoke',gene)
                    rseqhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Spoke','Desc']))
                else: rseqhtml.append('<i>No signficiant RSeq results.</i>')
                tablist.append(('RSeq',rje.join(rseqhtml,'\n'),'Randomised sequence (RSeq) datasets returning significant %s Spoke results' % gene))
            ## ~ [1g] ~ RUPC Results tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if gene in self.data('DBXRef'): 
                rseqhtml = ['<h2>%s RUPC results summary</h2>' % gene]
                maindb = self.db().getTable('Main_rupc')
                occdb = self.db().getTable('Occ_rupc')
                if gene in occdb.index('Spoke'):
                    headers = occdb.fields()
                    data = occdb.subset('Spoke',gene)
                    rseqhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Spoke','Desc']))
                else: rseqhtml.append('<i>No signficiant RUPC results.</i>')
                tablist.append(('RUPC',rje.join(rseqhtml,'\n'),'Randomised UPC (RUPC) datasets returning significant %s Spoke results' % gene))
            #- Modify the old code (below)

            ### ~ [2] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'gene/'),gene)
            html += tabberHTML(gene,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('genePage Error')
#########################################################################################################################
    def domainPage(self,domain):    ### Generates HTML for domain page
        '''
        Generates HTML for domain page.
        >> domain:str = domain for page construction
        '''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html = [htmlHead(domain,self.list['StyleSheets']),
                    '','<!-- ~~~~ %s Domain Summary details table ~~~~ -->' % domain,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="70%%"><a href="../domain/%s.html">%s<b>%s</b></FONT></a></td>' % (domain,gfont,domain),
                    '<td width="10%%"></td>',
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>']
            html += ['<tr><td colspan=2><a href="http://pfam.sanger.ac.uk/family/%s"><img src="../resources/pfam.png" title="Link to PFam entry"></a></td></tr>' % (domain)]
            html += ['<tr><td colspan=2>%s<p>%s PFam family HMM</p></FONT></td></tr>' % (dfont,domain)]
            html += ['</table>','<!-- ~~~~ End %s Domain Summary details table ~~~~ -->' % domain,'']
            ## ~ [1a] Add GO annotation for spoke proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pfam = self.db().getTable('PFam')
            pseq = pfam.dataList(pfam.indexEntries('Type',domain),'Name')
            dbxref = self.db().getTable('DBXRef')
            gtab = []
            for gtype in ['BP','CC','MF']:
                gdict = {}; gcount = {}
                for spoke in pseq:
                    try:
                        ensg = dbxref.indexEntries('EnsLoci',spoke)[0]['EnsEMBL']
                        godata = self.dict['GO'][ensg]
                    except: continue
                    if gtype not in godata: continue
                    for gotup in godata[gtype]:
                        if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                        #gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                        gdict[gotup[1]] = goLink(gotup[1],gotup[0])
                        if gotup[1] not in gcount: gcount[gotup[1]] = 0
                        gcount[gotup[1]] += 1
                if gdict:
                    ghtml = []
                    for g in rje.sortKeys(gdict): ghtml.append('%s (%d)' % (gdict[g],gcount[g]))
                    gtab.append(('GO_%s' % gtype,rje.join(ghtml,' ~ '),self.titleText('go',gtype)))
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s.' % domain,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % domain,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % domain,'']
            ### ~ [2] Domain info Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = []
            ## ~ [2a] Gene Domain stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pgene = []
            for seq in pseq:
                if seq in dbxref.index('EnsLoci'):
                    for entry in dbxref.indexEntries('EnsLoci',seq):
                        if entry['Gene'] not in self.list['Genes']: continue
                        #pgene.append('<a href="../gene/%s.html">%s</a>' % (entry['Gene'],entry['Gene']))
                        pgene.append(self.geneLink(entry['Gene']))
                pgene.sort()
            ghtml = ['','<h2>Genes containing PFam Domains %s</h2>' % domain,rje.join(pgene,' ~ ')]
            if not pgene: ghtml += ['<i>None.</i>']
            tablist.append(('Genes',rje.join(ghtml,'\n'),'Genes containing PFam Domains %s' % domain))
            ## ~ [2b] ~ Interactors tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('Interactors',self.domainInteractorsTableHTML(domain),'Table of %s interactors' % domain))
            ## ~ [2c] ~ Interactome tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            inthtml = []
            for ptype in ['ppi','y2h','bin','com']:
                dataset = '%s.%sdom' % (domain,ptype)
                if dataset not in self.db().getTable('Main').index('Dataset'): continue
                #intpng = rje.makePath('%sinteractome_png/%s.interactome.png' % (self.info['HTMLPath'],dataset),wholepath=True)
                if self.opt['SVG']:
                    svglink = '../interactome_svg/%s.interactome.svg' % dataset
                    title = 'Graphical representation of %s interactome. Larger interactomes may not be clear. See tables for details.' % dataset
                    ctab = '%s\n' % self.svgCode(svglink,title)
                    inthtml += ['<h2>%s Interactome</h2>' % dataset,ctab]
                else:
                    intpng = '../interactome_png/%s.interactome.png' % dataset
                    alt = '%s interactome image not yet made' % dataset
                    title = 'Graphical representation of %s interactome. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % dataset
                    inthtml += ['<h2>%s Interactome</h2>' % dataset,'<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)]
            if inthtml: tablist.append(('Interactome',rje.join(inthtml,'\n'),'Interactome graphics'))
            ## ~ [2c] ~ UPC tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tablist.append(('UPC',self.UPCHTML(domain,domain=True),'UPC clustering of spoke proteins'))
            ## ~ [2d] ~ Hub Results Summary tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hubhtml = ['<h2>%s Hub results summary</h2>' % domain]
            maindb = self.db().getTable('Main')
            occdb = self.db().getTable('Occ')
            if domain in maindb.index('Hub'):
                headers = maindb.fields()
                data = maindb.subset('Hub',domain)
                hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Match']))
                if domain in occdb.index('Hub'):
                    hubhtml += ['<h2>SLiM occurrences for %s hub results</h2>' % domain]
                    headers = occdb.fields()
                    data = occdb.subset('Hub',domain)
                    hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']))
            else: hubhtml.append('<i>No Hub Datasets.</i>')
            tablist.insert(0,('Hub',rje.join(hubhtml,'\n'),'%s Hub results summary tables' % domain))
            #-- Summary graphics?
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'domain/'),domain)
            #x#self.deBug(tablist)
            html = rje.join(html,'\n') + self.tabberHTML(domain,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('domainPage Error')
#########################################################################################################################
    def frontPageFT(self,slim): ### Extract SLiM FT HTML from relevant pages.
        '''Extract SLiM FT HTML from relevant pages.'''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            code = rje_slim.slimFromPattern(slim)
            mdb = self.db().getTable('Main')
            odb = self.db().getTable('Occ')
            html = ['\n<!-- ~~~~ %s SLiM Feature Summary table ~~~~ -->' % slim,'',
                    '<h3><a href="./slim/%s.html">%s</a></h3>' % (code,slim)]
            headers = mdb.fields()
            data = mdb.subset('Pattern',slim)
            html.append(rje.replace(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Match']),'../','./'))
            spokelist = []; rawspoke = []
            for entry in odb.indexEntries('Pattern',slim):
                if entry['Spoke'] != '-': spokelist.append(self.geneLink(entry['Spoke'],frontpage=True));
                else: spokelist.append(self.geneLink(self.gene(entry['Seq'],frontpage=True)));
                rawspoke.append(entry['Seq'])
            uspokelist = rje.sortUnique(spokelist); rawspoke = rje.sortUnique(rawspoke)
            for i in range(len(uspokelist)): uspokelist[i] = '%s (%d)' % (uspokelist[i],spokelist.count(uspokelist[i]))
            html += ['<p><b>Spokes (%d): </b>%s</p></hr width="50%%">' % (len(uspokelist),rje.join(uspokelist,' ~ '))]
            ### ~ [2] ~ UniFake Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'slim/'),code)
            try:
                uhtml = rje.split(open(hfile,'r').read(),'UniFake TabberTab')[1]
                uhtml = rje.split(uhtml,'PRE>')[1]
                html += ['<PRE>%sPRE>' % rje.replace(uhtml,'../','./')]
            except: html += ['<i>%s UniFake problem</i>' % slim]
            return rje.join(html,'\n')
        except: self.errorLog('%s frontPageFT Error' % slim)
#########################################################################################################################
    def slimPage(self,slim):    ### Generates HTML for domain page
        '''
        Generates HTML for domain page.
        >> domain:str = domain for page construction
        '''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            code = rje_slim.slimFromPattern(slim)
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            html = [htmlHead(slim,self.list['StyleSheets']),
                    '','<!-- ~~~~ %s SLiM Summary details table ~~~~ -->' % slim,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="40%%"><a href="../slim/%s.html">%s<b>%s</b></FONT></a></td>' % (code,gfont,slim),
                    '<td width="40%%">%s%s</FONT></td>' % (ifont,code),
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>']
            xrefhtml = []
            href = 'http://bioware.soton.ac.uk/slimdb/humsf10_slim/slim/%s.html' % code
            title = alt = '%s SLiM Network results' % slim
            src = '../resources/humsf.png'
            xrefhtml += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            href = 'http://bioware.ucd.ie/~compass/cgi-bin/formParser.py?name_server=slimsearch2&motif_str=%s' % (slim)
            title = alt = 'SLiMSearch human proteome with %s' % slim
            src = '../resources/ucd.gif'
            xrefhtml += ['<a href="%s"><img alt="%s" src="%s" align="BOTTOM" border="0" height="50" title="%s"></a>' % (href,alt,src,title)]
            if xrefhtml: html += ['<tr><td colspan=2>'] + xrefhtml + ['</td></tr>']
            if slim in self.dict['SlimDesc']: desc = self.dict['SlimDesc'][slim]
            else: desc = 'Unannotated/unknown motif.'
            html += ['</td></tr>','<tr><td colspan=2>%s<p>%s</p></FONT></td></tr>' % (dfont,desc)]
            mdb = self.db().getTable('Main')
            odb = self.db().getTable('Occ')
            #hublist = []
            #for entry in mdb.indexEntries('Pattern',slim): hublist.append('%s (%s)' % (self.geneLink(entry['Hub']),rje_slim.expectString(rje.atof(entry['Sig']))))
            #hublist.sort()
            #hublist = rje.join(hublist,' ~ ')
            #html += ['<tr><td colspan=2><p>%sHubs: </FONT>%s</p></td></tr>' % (dfont,hublist)]
            genelist = []; domlist = []; randlist = []
            for entry in mdb.indexEntries('Pattern',slim):
                if entry['Dataset'][-3:] == 'dom': domlist.append('%s (%s)' % (domainLink(entry['Hub']),rje_slim.expectString(rje.atof(entry['Sig']))))
                elif entry['Dataset'][:4] in ['rseq','rupc']: randlist.append('%s (%s)' % (randLink(entry['Dataset']),rje_slim.expectString(rje.atof(entry['Sig']))))
                else: genelist.append('%s (%s)' % (self.geneLink(entry['Hub']),rje_slim.expectString(rje.atof(entry['Sig']))))
            genelist.sort(); domlist.sort()
            if genelist: html += ['<tr><td colspan=2><p>%sGene Hubs (%d): </FONT>%s</p></td></tr>' % (dfont,len(genelist),rje.join(genelist,' ~ '))]
            else: html += ['<tr><td colspan=2><p>%sGene Hubs: </FONT><i>None.</i></p></td></tr>' % (dfont)]
            if domlist: html += ['<tr><td colspan=2><p>%sDomain Hubs (%d): </FONT>%s</p></td></tr>' % (dfont,len(domlist),rje.join(domlist,' ~ '))]
            else: html += ['<tr><td colspan=2><p>%sDomain Hubs: </FONT><i>None.</i></p></td></tr>' % (dfont)]
            if randlist: html += ['<tr><td colspan=2><p>%sRandom Datasets (%d): </FONT>%s</p></td></tr>' % (dfont,len(randlist),rje.join(randlist,' ~ '))]
            else: html += ['<tr><td colspan=2><p>%sRandom Datasets: </FONT><i>None.</i></p></td></tr>' % (dfont)]
            spokelist = []; rawspoke = []
            for entry in odb.indexEntries('Pattern',slim):
                if entry['Spoke'] != '-': spokelist.append(self.geneLink(entry['Spoke']));
                else: spokelist.append(self.geneLink(self.gene(entry['Seq'])));
                rawspoke.append(entry['Seq'])
            uspokelist = rje.sortUnique(spokelist); rawspoke = rje.sortUnique(rawspoke)
            for i in range(len(uspokelist)): uspokelist[i] = '%s (%d)' % (uspokelist[i],spokelist.count(uspokelist[i]))
            html += ['<tr><td colspan=2><p>%sSpokes (%d): </FONT>%s</p></td></tr>' % (dfont,len(uspokelist),rje.join(uspokelist,' ~ '))]
            html += ['</table>','<!-- ~~~~ End %s SLiM Summary details table ~~~~ -->' % slim,'']
            ## ~ [1a] Add GO annotation for spoke proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dbxref = self.db().getTable('DBXRef')
            gtab = []
            #self.deBug(rawspoke)
            for gtype in ['BP','CC','MF']:
                gdict = {}; gcount = {}
                for spoke in rawspoke:
                    try:
                        #self.deBug(spoke)
                        ensg = dbxref.indexEntries('EnsLoci',spoke)[0]['EnsEMBL']#; self.deBug(ensg)
                        godata = self.dict['GO'][ensg]#; self.deBug(godata)
                    except: continue
                    if gtype not in godata: continue
                    for gotup in godata[gtype]:
                        if gotup[1] in ['cellular_component','biological_process','molecular_function']: continue
                        #gdict[gotup[1]] = '<a href="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=GO:%s">%s</a>' % (gotup[0],gotup[1])
                        gdict[gotup[1]] = goLink(gotup[1],gotup[0])
                        if gotup[1] not in gcount: gcount[gotup[1]] = 0
                        gcount[gotup[1]] += 1
                if gdict:
                    ghtml = []
                    for g in rje.sortKeys(gdict):
                        if gcount[g] > 1: ghtml.append('%s (%d)' % (gdict[g],gcount[g]))
                    gtab.append(('GO_%s' % gtype,rje.join(ghtml,' ~ '),self.titleText('go',gtype)))
            if gtab:
                gtab.insert(0,('^','Click on GO tabs for Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) GO terms associated with %s (2+ hubs/spokes).' % slim,'Compress GO tabs'))
                html += ['','<!-- ~~~~ %s GO tabs ~~~~ -->' % slim,tabberHTML('GO',gtab),
                         '<!-- ~~~~ End %s GO tabs ~~~~ -->' % slim,'']
            #else:
            #    try: self.deBug(dbxref.indexEntries('EnsLoci',rawspoke[0]))
            #    except: self.deBug(rawspoke)
            ### ~ [2] Domain info Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = []
            ## ~ [2a] ~ Hub Results Summary tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hubhtml = ['<h2>%s Hub results summary</h2>' % slim]
            headers = mdb.fields()
            data = mdb.subset('Pattern',slim)
            hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Match']))
            tablist.append(('Hub',rje.join(hubhtml,'\n'),'%s Hub results summary' % slim))
            #-- Summary graphics?
            ## ~ [2b] ~ Spoke Summary Tab (SLiM Occurrences in "Gene") ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sphtml = ['<h2>%s Spoke occurrences</h2>' % slim]
            headers = odb.fields()
            data = odb.subset('Pattern',slim)
            sphtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']))
            tablist.append(('Spoke',rje.join(sphtml,'\n'),'%s Spoke results summary' % slim))
            #-- Summary graphics?
            ## ~ [2c] ~ CompariMotif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chtml = ['<h2>%s CompariMotif search results</h2>' % slim]
            cdb = self.db().getTable('CompariMotif')
            if slim in cdb.index('Pattern'):
                chtml.append(self.resultTableHTML(cdb.fields(),cdb.indexEntries('Pattern',slim),asdict=False,drop=['Pattern','Rank']))
            else: chtml.append('<i>No CompariMotif hits.</i>')
            tablist.append(('CompariMotif',rje.join(chtml,'\n'),'%s CompariMotif search results' % slim))
            ## ~ [2d] ~ Interactome tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            inthtml = []
            if self.opt['SVG']:
                svglink = '../nested_svg/%s.svg' % code
                title = 'Graphical representation of %s interactome. Larger interactomes may not be clear. See tables for details.' % slim
                ctab = '%s\n' % self.svgCode(svglink,title)
                inthtml += ['<h2>%s Interactome</h2>' % slim,ctab]
            else:
                intpng = '../nested_png/%s.png' % code
                alt = '%s interactome image not yet made' % code
                title = 'Graphical representation of %s interactome, Hubs in Centre and Spokes around edge. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % code
                inthtml += ['<h2>%s Interactome</h2>' % slim,'<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)]
            tablist.append(('Interactome',rje.join(inthtml,'\n'),'Interactome graphics'))
            ## ~ [2e] ~ UniFake tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ufhtml = ['<h2>UniFake features overlapping motif occurrences.</h2>','<PRE>']
            for seqid in rawspoke:
                ufile = rje.makePath(self.info['HTMLPath']+'unifake') + '%s.html' % seqid
                occpos = []
                for occ in data.values():
                    if occ['Seq'] != seqid: continue
                    occpos.append((int(occ['Start_Pos']),int(occ['End_Pos'])))
                dx = len(ufhtml)
                if os.path.exists(ufile):
                    uhtml = open(ufile,'r').read()
                    uhtml = uhtml[uhtml.find('<PRE>'):][5:]
                    uhtml = uhtml[:uhtml.find('</PRE>')]
                    uhtml = rje.split(uhtml,'\n')
                    for ft in uhtml:
                        if not ft: continue
                        ftsplit = rje.split(ft)
                        if not ftsplit: continue
                        if ftsplit[0] not in ['ID','DE','GN','FT']: continue
                        if ftsplit[0] == 'GN' and 'href=' not in ft:
                            id = rje.matchExp('^(GN\s+Name=)(\S+)(;.*)',ft)
                            if id: ft = rje.join([id[0],self.geneLink(id[1]),id[2]],'')
                        if ftsplit[0] != 'FT': ufhtml.append(ft); continue
                        if ftsplit[1] in ['CHAIN']: continue
                        if ftsplit[1] == 'SLIM':
                            try: start = rje.atoi(ftsplit[2]); end = rje.atoi(ftsplit[3])
                            except:
                                try:
                                    psplit = rje.split(rje.join(rje.split(ft,'>')[1:],'<'),'<')
                                    start = rje.atoi(psplit[0])
                                    end = rje.atoi(psplit[4])
                                except: continue
                        else:
                            try: start = rje.atoi(ftsplit[2]); end = rje.atoi(ftsplit[3])
                            except: continue
                        overlap = False
                        for (pstart,pend) in occpos:
                            if pstart <= end and pend >= end: overlap = True
                            elif pend >= start and pstart <= start: overlap = True
                            elif start <= pend and end >= pend: overlap = True
                            elif end >= pstart and start <= pstart: overlap = True
                        if overlap: ufhtml.append(ft)
                    ufhtml.append('//')
                else: ufhtml.append('<i>%s UniFake file to be generated.</i>' % seqid)
                #self.deBug(rje.join(ufhtml[dx:],'\n'))
            ufhtml.append('</PRE>')
            tablist.append(('UniFake',rje.join(ufhtml,'\n'),'Features from modified %s spoke UniProt entries' % slim))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'slim/'),code)
            html = rje.join(html,'\n') + tabberHTML(code,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('slimPage Error')
#########################################################################################################################
    def randPage(self,dataset):    ### Generates HTML for domain page
        '''
        Generates HTML for random dataset page.
        >> dataset:str = domain for page construction
        '''
        try:### ~ [1] ~ Setup & Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfont = '<FONT SIZE=6 FACE="Verdana" COLOR=#014359>'
            ifont = '<FONT SIZE=5 FACE="Verdana" COLOR=#979E45>'
            dfont = '<FONT SIZE=4 COLOR=#014359>'
            rtype = dataset[:4]
            html = [htmlHead(dataset,self.list['StyleSheets']),
                    '','<!-- ~~~~ %s Random Dataset Summary details table ~~~~ -->' % dataset,'<table width="100%">',
                    '<tr valign="top">'  #<td width="80%">', '<table>','<tr>',
                    '<td width="70%%"><a href="../%s/%s.html">%s<b>%s</b></FONT></a></td>' % (rtype,dataset,gfont,dataset),
                    '<td width="10%%"></td>',
                    '<td width="20%" align="right" rowspan="3">',
                    '<a href="../index.htm"><img src="../resources/SBS_100.png" height="100" align="RIGHT" border="0" alt="Home"></a>',
                    '</td></tr>']
            html += ['</table>','<!-- ~~~~ End %s Summary details table ~~~~ -->' % dataset,'']
            #!# Consider adding GO annotation for random dataset proteins? #!#
            ### ~ [2] Domain info Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tablist = []
            ## ~ [2a] Gene Domain stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rseq = self.db().getTable('DSetSeq')
            dbxref = self.db().getTable('DBXRef'); dbxref.index('EnsLoci')
            pseq = rseq.dataList(rseq.indexEntries('Dataset',dataset),'Seq')
            pgene = []
            for seq in pseq:
                #if seq in dbxref.index('EnsLoci'):
                #    for entry in dbxref.indexEntries('EnsLoci',seq):
                #        if entry['Gene'] not in self.list['Genes']: pgene.append(entry['Gene'])
                #        else: pgene.append('<a href="../gene/%s.html">%s</a>' % (entry['Gene'],entry['Gene']))
                pgene.append(self.geneLink(self.gene(seq)))
            pgene.sort()
            ghtml = ['','<h2>Genes in %s</h2>' % dataset,rje.join(pgene,' ~ ')]
            if not pgene: ghtml += ['<i>None.</i>']
            tablist.append(('Genes',rje.join(ghtml,'\n'),'Genes in %s' % dataset))
            ## ~ [2b] ~ Interactors tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Consider adding indicator of shared interactors?
            ## ~ [2c] ~ Interactome tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#tablist.append(('Interactome','<i>Add interactome image</i>'))
            #- Graphical view of interactome with simple list of interactors (linking to pages)
            ## ~ [2c] ~ UPC tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #dset = '%s.domppi' % domain; upc = None
            tablist.append(('UPC',self.UPCHTML(dataset),'UPC clustering of spoke proteins'))
            ## ~ [2d] ~ Hub Results Summary tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hubhtml = ['<h2>%s results summary</h2>' % dataset]
            maindb = self.db().getTable('Main')
            occdb = self.db().getTable('Occ')
            headers = maindb.fields()
            data = maindb.subset('Dataset',dataset)
            hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']+['Match']))
            if dataset in occdb.index('Dataset'):
                hubhtml += ['<h2>SLiM occurrences for %s results</h2>' % dataset]
                headers = occdb.fields()
                data = occdb.subset('Dataset',dataset)
                hubhtml.append(self.resultTableHTML(headers,data,True,drop=self.list['DropFields']))
            tablist.insert(0,('Sig',rje.join(hubhtml,'\n'),'SLiM occurrences for %s results' % dataset))
            ## ~ [2d] ~ Interactome tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            inthtml = []
            if self.opt['SVG']:
                svglink = '../interactome_svg/%s.interactome.svg' % dataset
                title = 'Graphical representation of %s interactome. Larger interactomes may not be clear. See tables for details.' % dataset
                ctab = '%s\n' % self.svgCode(svglink,title)
                inthtml += ['<h2>%s Interactome</h2>' % dataset,ctab]
            else:
                intpng = '../interactome_png/%s.interactome.png' % dataset
                alt = '%s interactome image not yet made' % dataset
                title = 'Graphical representation of %s interactome. Click to open in new window. Larger interactomes may not be clear. See tables for details.' % dataset
                inthtml += ['<h2>%s Interactome</h2>' % dataset,'<a href="%s" target="_blank"><img src="%s" alt="%s" title="%s" width=1000></a>' % (intpng,intpng,alt,title)]
            tablist.append(('Interactome',rje.join(inthtml,'\n'),'Interactome graphics'))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + '%s/' % rtype),dataset)
            html = rje.join(html,'\n') + self.tabberHTML(dataset,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
        except: self.errorLog('randPage Error')
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
            godata = self.obj['GO'].go(goid)
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
            ## ~ [2a] ~ GO parents  & children List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            for rtype in ['is_a','part_of']:
                if rtype not in godata: continue
                try:
                    self.deBug(godata[rtype])
                    for rel in godata[rtype]: tabhtml.append('<li>%s %s %s</li>' % (go,rtype,goLink(self.obj['GO'].go(rel)['name'],rel)))                       
                except:
                    self.errorLog('??')
                    pass
            if tabhtml: tabhtml.append('<br>')
            self.deBug(godata)
            for id in godata['child_terms']:
                for rtype in ['is_a','part_of']:
                    if rtype not in self.obj['GO'].go(id): continue
                    try:
                        for rel in self.obj['GO'].go(id)[rtype]:
                            if rel == goid: tabhtml.append('<li>%s %s %s</li>' % (goLink(self.obj['GO'].go(id)['name'],id),rtype,self.obj['GO'].go(rel)['name']))                       
                    except: pass
            self.deBug(tabhtml)
            if tabhtml:
                tabhtml = ['<H2>GO Relationships</H2>','<ul>'] + tabhtml + ['</ul>']
                tablist.append(('GO',rje.join(tabhtml,'\n'),'GO Term relationships for %s' % go))
            ## ~ [2a] ~ Gene List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            genes.sort()
            for gene in genes: tabhtml.append(self.geneLink(gene))
            if tabhtml: tabhtml[0] = '<h2>%s Genes</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No genes in analysis mapped to %s</i>' % go]
            tablist.append(('Genes',rje.join(tabhtml,jtxt),'Genes mapped to %s' % go))
            ## ~ [2b] ~ Domain List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            domains.sort()
            for hub in domains: tabhtml.append(domainLink(hub))
            if tabhtml: tabhtml[0] = '<h2>%s Domains</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No domains in analysis mapped to %s</i>' % go]
            tablist.append(('Domains',rje.join(tabhtml,jtxt),'Domains mapped to %s' % go))
            ## ~ [2c] ~ SigHub List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            sighubs.sort()
            for hub in sighubs:
                if hub in genes: tabhtml.append(self.geneLink(hub))
                else: tabhtml.append(domainLink(hub))
            if tabhtml: tabhtml[0] = '<h2>%s Hubs with significant results</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s hubs returned significant results</i>' % go]
            tablist.append(('Hubs',rje.join(tabhtml,jtxt),'Significant Hubs mapped to %s' % go))
            ## ~ [2d] ~ SigSpoke List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            sigspokes.sort()
            for prot in sigspokes:
                gene = self.gene(prot)
                if gene != prot: tabhtml.append(self.geneLink(gene))
                else: tabhtml.append(prot)
            if tabhtml: tabhtml[0] = '<h2>%s Spokes with significant results</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s spokes returned significant results</i>' % go]
            tablist.append(('Spokes',rje.join(tabhtml,jtxt),'Significant Spokes mapped to %s' % go))
            ## ~ [2e] ~ SLiM List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabhtml = []
            slims.sort()
            for slim in slims: tabhtml.append(slimLink(slim))
            if tabhtml: tabhtml[0] = '<h2>%s SLiM predictions</h2><p>\n%s' % (go,tabhtml[0])
            else: tabhtml = ['<i>No %s significant SLiM predictions</i>' % go]
            tablist.append(('SLiMs',rje.join(tabhtml,jtxt),'Significant SLiM predictions mapped to %s' % go))
            ### ~ [3] ~ Generate HTML page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s%s.html' % (rje.makePath(self.info['HTMLPath'] + 'go/'),goid)
            html = rje.join(html,'\n') + tabberHTML(goid,tablist)
            html += htmlTail()
            open(hfile,'w').write(html)
            self.deBug(html)
        except: self.errorLog('goPage Error')
#########################################################################################################################
    ### <4> ### Tables pages                                                                                            #
#########################################################################################################################
    def interactorTableHTML(self,gene): ### Table of interactors with evidence & SLiMs returned for 4 datasets types
        '''
        Table of interactors with evidence & SLiMs returned for 4 datasets types ("MeGene Hub SLiMs", "MeGene Spoke SLiMs")
        - Genes link to corresponding gene pages
        - MeGene headers - colspan=4
        - "-" = No dataset, "None" = No SLiMs, else SLiM (sig); SLiM (sig) etc.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = []
            ppi = self.db().getTable('PPI')
            hdb = ppi.subset('Hub',gene)
            sdb = ppi.subset('Spoke',gene)
            occ = self.db().getTable('Occ')
            hocc = occ.subset('Hub',gene)
            socc = occ.subset('Spoke',gene)
            ilist = rje.sortUnique(ppi.dataList(hdb.values(),'Spoke') + ppi.dataList(sdb.values(),'Hub'))   # Full list of interactors
            ### ~ [2] Generate table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Table Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = ['','<h2>Protein-Protein Interactions</h2>','',
                    '<table width=1000 border=%d><tr><th colspan=2>%s Interactions</th><th colspan=4>%s Hub SLiMs</th><th colspan=4>%s Spoke SLiMs</th></tr>' % (self.stat['Border'],gene,gene,gene)]
            html += ['<tr><th align=left title="%s">Gene</th>' % self.titleText('ppi','gene'),
                     '<th align=left title="%s">Evidence</th>' % self.titleText('ppi','evidence')]
            for i in range(2):
                for ptype in ['ppi','y2h','bin','com']: html.append('<th title="%s">%s</th>' %  (self.titleText('ppi',ptype),ptype))
            html.append('</tr>')
                    #rje.join(['<tr><th align=left>Gene</th><th align=left>Evidence'] + ['ppi','bin','com','y2h'] * 2 ,'</th><th>' )] + ['</th></tr>']
            ## ~ [2b] Table contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in ilist[0:]:
                html += ['<tr valign="top"><td width=90>%s</td>' % self.geneLink(i)]
                hkey = '%s\t%s' % (gene,i); skey = '%s\t%s' % (i,gene)
                if hkey in hdb: html += ['<td width=350>%s</td>' % rje.join(rje.split(hdb[hkey]['Evidence'],'|'),';<br>')]
                else: html += ['<td width=350>%s</td>' % rje.join(rje.split(sdb[skey]['Evidence'],'|'),';<br>')]
                ## [i] Hub SLiMs ##
                for type in ['ppi','y2h','bin','com']:
                    if hkey in hdb and hdb[hkey][type] == 'Y':
                        dset = '%s.%s' % (gene,type)
                        slist = []
                        for entry in occ.entryList(rje.sortKeys(hocc)):
                            if entry['Dataset'] == dset and entry['Spoke'] == i:
                                sres = '%s (%s)' % (slimLink(entry['Pattern']),rje_slim.expectString(rje.atof(entry['Sig'])))
                                if sres not in slist: slist.append(sres)
                        if slist: html += ['<td align=center width=70>%s</td>' % rje.join(slist,'; ')]
                        else: html += ['<td align=center width=70>None</td>']
                    else: html += ['<td align=center width=70>-</td>']
                ## [ii] Spoke SLiMs ##
                for type in ['ppi','y2h','bin','com']:
                    if skey in sdb and sdb[skey][type] == 'Y':
                        dset = '%s.%s' % (i,type)
                        slist = []
                        for entry in occ.entryList(rje.sortKeys(socc)):
                            if entry['Dataset'] == dset and entry['Spoke'] == gene:
                                slist.append('%s (%s)' % (slimLink(entry['Pattern']),rje_slim.expectString(rje.atof(entry['Sig']))))
                        if slist: html += ['<td align=center width=70>%s</td>' % rje.join(slist,'; ')]
                        else: html += ['<td align=center width=70>None</td>']
                    else: html += ['<td align=center width=70>-</td>']
                html += ['</tr>','']
        except: self.errorLog('interactorTableHTML Error')
        ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += ['</table>','']
        return rje.join(html,'\n')
#########################################################################################################################
    def domainTableHTML(self,gene): ### Table of interactors with evidence & SLiMs returned for 4 datasets types
        '''
        Table of interactors with evidence & SLiMs returned for 4 datasets types ("MeGene Hub SLiMs", "MeGene Spoke SLiMs")
        - Genes link to corresponding gene pages
        - MeGene headers - colspan=4
        - "-" = No dataset, "None" = No SLiMs, else SLiM (sig); SLiM (sig) etc.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = []
            ppi = self.db().getTable('DomPPI')
            sdb = ppi.subset('Spoke',gene)
            occ = self.db().getTable('Occ')
            socc = occ.subset('Spoke',gene)
            ilist = rje.sortUnique(ppi.dataList(sdb.values(),'Domain'))   # Full list of interactors
            ## ~ [1a] Gene Domain stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ddb = {}
            docc = {}
            pfam = self.db().getTable('PFam')
            try:
                seq = self.data('DBXRef')[gene]['EnsLoci']
                pdom = pfam.subset('Name',seq)
                for domain in pfam.dataList(pdom.values(),'Type'):
                    ddb[domain] = pfam.subset('Type',domain)    # Subset of original dictionary {key:datadict}
                    docc[domain] = occ.subset('Hub',domain)
            except: seq = None
            ## ~ [1b] Gene domain HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = ['','<h2>PFam Domains in %s</h2>' % gene]
            if ddb:
                dlist = []
                for domain in rje.sortKeys(ddb): dlist.append(domainLink(domain))
                html += ['<p>%s' % rje.join(dlist+['</p>'],'; ')]
            elif gene in self.db().getTable('DBXRef').index('Gene'): html += ['<p><i>None.</i></p>']
            else: html += ['<p><i>None. (No sequence mapped.)</i></p>']
            #html += ['<hr width=60%>','']
            ### ~ [2] Generate table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Table Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html += ['','<h2>Domain-Protein Interactions</h2>','',
                     '<table border=%d width=1000><tr><th colspan=3>%s Interactions</th><th colspan=4>%s Spoke SLiMs</th></tr>' % (self.stat['Border'],gene,gene)]
            html += ['<tr><th align=left title="%s">Domain</th>' % self.titleText('ppi','domain'),
                     '<th align=left title="%s">Interacting Hubs</th>' % self.titleText('ppi','domhubs'),
                     '<th align=left title="%s">Evidence</th>' % self.titleText('ppi','evidence')]
            for ptype in ['ppi','y2h','bin','com']: html.append('<th title="%s">%s</th>' %  (self.titleText('ppi',ptype),ptype))
            html.append('</tr>')
                     #rje.join(['<tr><th align=left>Domain</th><th align=left>Interacting Hubs</th><th align=left>Evidence'] + ['ppi','bin','com','y2h'],'</th><th>' )] + ['</th></tr>']
            ## ~ [2b] Table contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in ilist[0:]:
                skey = '%s\t%s' % (i,gene)
                html += ['<tr valign="top"><td width=150>%s</td>' % domainLink(i)]
                hlist = []
                for hub in rje.split(sdb[skey]['Hubs'],','): hlist.append(self.geneLink(hub))
                html += ['<td width=220>%s;</td>' % rje.join(hlist,'; ')]
                html += ['<td width=350>%s</td>' % rje.join(rje.split(sdb[skey]['Evidence'],'|'),'; ')]
                ## [i] Spoke SLiMs ##
                for type in ['ppi','y2h','bin','com']:
                    if skey in sdb and sdb[skey][type] == 'Y':
                        dset = '%s.%sdom' % (i,type)
                        slist = []
                        for entry in occ.entryList(rje.sortKeys(socc)):
                            if entry['Dataset'] == dset and entry['Spoke'] == gene:
                                sres = '%s (%s)' % (slimLink(entry['Pattern']),rje_slim.expectString(rje.atof(entry['Sig'])))
                                if sres not in slist: slist.append(sres)
                        if slist: html += ['<td width=70>%s</td>' % rje.join(slist,'; ')]
                        else: html += ['<td width=70>None</td>']
                    else: html += ['<td align=center width=70>-</td>']
                html += ['</tr>','']
        except: self.errorLog('domainTableHTML Error')
        ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += ['</table>','']
        return rje.join(html,'\n')
#########################################################################################################################
    def domainInteractorsTableHTML(self,domain):    ### Table of interactors with evidence & SLiMs returned for 4 datasets types
        '''
        Table of interactors with evidence & SLiMs returned for 4 datasets types ("MeGene Hub SLiMs")
        - Genes link to corresponding gene pages
        - MeGene headers - colspan=4
        - "-" = No dataset, "None" = No SLiMs, else SLiM (sig); SLiM (sig) etc.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = []
            pdb = self.db().getTable('DomPPI')
            ppi = pdb.subset('Domain',domain)
            occ = self.db().getTable('Occ') #.subset('Hub',domain)
            ilist = rje.sortUnique(pdb.dataList(ppi.values(),'Spoke'))   # Full list of interactors
            ### ~ [2] Generate table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Table Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html += ['','<h2>Domain-Protein Interactions</h2>','',
                     '<table border=%d width=1000><tr><th colspan=3>%s Interactions</th><th colspan=4>%s Hub SLiMs</th></tr>' % (self.stat['Border'],domain,domain)]
            html += ['<tr><th align=left title="%s">Gene</th>' % self.titleText('ppi','gene'),
                     '<th align=left title="%s">Interacting Hubs</th>' % self.titleText('ppi','domhubs'),
                     '<th align=left title="%s">Evidence</th>' % self.titleText('ppi','evidence')]
            for ptype in ['ppi','y2h','bin','com']: html.append('<th title="%s">%s</th>' %  (self.titleText('ppi',ptype),ptype))
            html.append('</tr>')
                     #rje.join(['<tr><th align=left>Gene</th><th align=left>Interacting Hubs</th><th align=left>Evidence'] + ['ppi','bin','com','y2h'],'</th><th>' )] + ['</th></tr>']
            ## ~ [2b] Table contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in ilist[0:]:
                skey = '%s\t%s' % (domain,i)
                html += ['<tr valign="top"><td width=150>%s</td>' % self.geneLink(i)]
                hlist = []
                for hub in rje.split(ppi[skey]['Hubs'],','): hlist.append(self.geneLink(hub))
                html += ['<td width=220>%s;</td>' % rje.join(hlist,'; ')]
                html += ['<td width=350>%s</td>' % rje.join(rje.split(ppi[skey]['Evidence'],'|'),'; ')]
                ## [i] Spoke SLiMs ##
                for type in ['ppi','y2h','bin','com']:
                    if skey in ppi and ppi[skey][type] == 'Y':
                        dset = '%s.%sdom' % (domain,type)
                        slist = []
                        for entry in occ.indexEntries('Hub',domain):
                            if entry['Dataset'] == dset and entry['Spoke'] == i:
                                sres = '%s (%s)' % (slimLink(entry['Pattern']),rje_slim.expectString(rje.atof(entry['Sig'])))
                                if sres not in slist: slist.append(sres)
                        if slist: html += ['<td width=70>%s</td>' % rje.join(slist,'; ')]
                        else: html += ['<td width=70>None</td>']
                    else: html += ['<td align=center width=70>-</td>']
                html += ['</tr>','']
        except: self.errorLog('domainInteractorsTableHTML Error')
        ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        html += ['</table>','']
        return rje.join(html,'\n')
#########################################################################################################################
    def pround(self,x):
        if self.opt['PRound']: p = 10 * int(x/10.0) + 5
        else: p = x
        return rje.preZero(p,10000)
#########################################################################################################################
    def resultTableHTML(self,headers,data,asdict=True,drop=[]):     ### Write data to HTML table
        '''Write data to HTML table.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = headers[0:]   # Do not want to change original headers list!
            for h in ['RunID','Masking','Build','RunTime','Chance','MotifNum','MotNum','Occ','UP','SigRank'] + drop:
                if h in headers: headers.remove(h)
            if asdict:
                datalist = []
                for key in rje.sortKeys(data): datalist.append(data[key])
            else: datalist = data
            ### ~ [2] ~ Write Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hcode = ''
            for h in headers: hcode += '<th title="%s">%s</th>\n' % (self.titleText('results',h),h)
            rows = [hcode]
            #rows = ['<TD><B>%s</B></TD>' % rje.join(headers,'</B></TD><TD><B>')]
            ## ~ [2b] ~ Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for dat in datalist:
                rdata = []
                for h in headers:
                    if h in dat:
                        if h in ['Cons','GlobID','LocID','Hyd','IUP','SA','IC'] or h[-4:] == 'mean':
                            try: rdata.append('%.2f' % rje.atof(dat[h]))
                            except: rdata.append('-')
                        elif h in ['Sig','FDR','dConn','sConn','ExpUP']:
                            try:
                                rdata.append(rje_slim.expectString(rje.atof(dat[h])))
                            except: rdata.append('-')
                        elif h == 'Pattern' and dat['Pattern'] != '-':
                            rdata.append('<a href="../slim/%s.html" title="%s Results page">%s</A>' % (rje_slim.slimFromPattern(dat[h]),dat[h],dat[h]))
                        elif h in ['Dataset','Hub','Seq','Spoke']:
                            if h == 'Dataset': gene = rje.join(rje.split(dat['Dataset'],'.')[:-1],'.')
                            elif h == 'Seq' and 'Spoke' in dat: gene = dat['Spoke']
                            else: gene = str(dat[h])
                            if h not in ['Seq','Spoke'] and 'Dataset' in dat and dat['Dataset'][-3:] == 'dom': rdata.append('<a href="../domain/%s.html" title="%s Results page">%s</A>' % (gene,gene,dat[h]))
                            elif h in ['Dataset','Hub'] and dat[h][:4] in ['rseq','rupc']:
                                if 'Dataset' in dat: rdata.append('<a href="../%s/%s.html" title="%s Results page">%s</A>' % (dat[h][:4],dat['Dataset'],dat['Dataset'],dat[h]))
                                else: rdata.append(str(dat[h]))
                            else: rdata.append(self.geneLink(gene,frontpage=False,altgene=dat[h]))
                        elif h == 'Start_Pos':
                            p = rje.atoi(dat[h])
                            pngdir = '../map_png/'
                            try: basefile = '%s%s.%s' % (pngdir,dat['Dataset'],dat['Seq'])
                            except: rdata.append('%d' % p); continue
                            ofile = '%s.%s.png' % (basefile,self.pround(p))
                            rdata.append('<a href="%s" target="_blank" title="Alignment of motif occurrence">%d</a>' % (ofile,p))
                        elif h in ['Rank','End_Pos']:
                            try: rdata.append('%d' % rje.atoi(dat[h]))
                            except: rdata.append(str(dat[h]))
                            if h == 'Rank':
                                if 'MotifNum' in dat: rdata[-1] = '%s / %s' % (rdata[-1],dat['MotifNum'])
                                elif 'MotNum' in dat: rdata[-1] = '%s / %s' % (rdata[-1],dat['MotNum'])
                        elif h == 'Support':
                            pngdir = '../motifaln_png/'
                            rdata.append('%s / %s / %s' % (dat['Occ'],dat[h],dat['UP']))
                            try: pfile = rje.makePath('%s%s.%d.%s.png' % (pngdir,dat['Dataset'],int(dat['Rank']),rje_slim.slimFromPattern(dat['Pattern'])),wholepath=True)
                            except: continue
                            rdata[-1] = '<a href="%s" title="Graphic of motif occurrences" target="_blank">%s</a>' % (pfile,rdata[-1])
                        elif h == 'Class' and dat[h]:
                            tdat = self.tableClassLink(dat)
                            if tdat: rdata.append(tdat)
                            else: rdata.append(dat[h])
                        else: rdata.append(str(dat[h]))
                        if h in ['Sig','FDR','Support'] and 'Pattern' in dat and dat['Pattern'] == '-': rdata[-1] = '-'
                    else: rdata.append('')
                    #if h in ['Rank','Sig']:
                    #    try: # Add link to PNG for Rank and Sig
                    #        if dat['Rank'] in ['0',0]: continue
                    #        pngdir = '%smotifaln_png/' % self.info['HTMLPath']
                    #        pfile = '%s.%s.%s.png' % (dat['Dataset'],sdata['Rank'],rje_slim.slimFromPattern(sdata['Pattern']))
                    #        if not os.path.exists('%s%s' % (pngdir,pfile)): continue
                    #        rdata[-1] = '<a href="../motifaln_png/%s" title="Graphic of motif occurrences">%s</a>' % (pfile,rdata[-1])
                    #    except: continue
                rows.append('<TD>%s</TD>' % rje.join(rdata,'</TD><TD>'))
        except: self.errorLog('resultTableHTML error')
        return '<TABLE BORDER=%d><TR VALIGN="top">\n%s\n</TR></TABLE>\n' % (self.stat['Border'],rje.join(rows,'\n</TR><TR VALIGN="top">\n'))
#########################################################################################################################
    ### <5> ### Misc HTML pages                                                                                         #
#########################################################################################################################
    def checkHTML(self,hpage):  ### Checks for existence of complete HTML page
        '''Checks for existence of complete HTML page.'''
        if not os.path.exists(hpage): return False
        html = open(hpage,'r').read().upper()
        if html.find('<HTML>') < 0 or html.find('</HTML>') < 0: return False
        return True
#########################################################################################################################
    def UPCHTML(self,hub,domain=False):  ### Makes HTML for UPC file, linking to genes
        '''Makes HTML for UPC file, linking to genes.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            upc = []; ulist = ['Full']; uhtml = []
            dtype = 'Real'
            if hub[:4] in ['rseq','rupc']: dtype = {'rupc':'RUPC','rseq':'RSeq'}[hub[:4]]
            for ptype in ['ppi','y2h','bin','com']:
                if dtype == 'Real':
                    dset = '%s.%s' % (hub,ptype)
                    if domain: dset = '%sdom' % dset
                else: dset = hub
                if os.path.exists('../2010-06-HumSF09_Data_Archive/HumSF09_%s_Sig/%s.upc' % (dtype,dset)): upc.append('../2010-06-HumSF09_Data_Archive/HumSF09_%s_Sig/%s.upc' % (dtype,dset))
                elif os.path.exists('../2010-06-HumSF09_Data_Archive/HumSF09_%s_NonSig/%s.upc' % (dtype,dset)): upc.append('../2010-06-HumSF09_Data_Archive/HumSF09_%s_NonSig/%s.upc' % (dtype,dset))
                if dtype != 'Real': break
                if ptype != 'ppi': ulist.append(ptype)
            if not upc: return '<p><i>UPC missing</i></p>'
            dbxref = self.db().getTable('DBXRef'); dbxref.index('EnsLoci')
            ### ~ [2] ~ Content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for upc in upc[0:]:
                utxt = open(upc,'r').readlines()[0:]
                uhtml += ['<H2>%s PPI UPC</H2>' % ulist.pop(0),'<p>%s</p>' % rje.chomp(utxt[0]),
                         '<table border=%d width=1000><tr>' % self.stat['Border'],
                         '<th title="UP Number">UP</th><th title="No. Proteins in UP">N</th><th title="Corrected UP Size (Min. Spanning Tree)">MST</th><th title="Spoke Proteins in UP">Spokes</th>']
                for up in utxt[2:]:
                    udet = rje.split(up)
                    uhtml += ['</tr><tr>','  <td width=30>%s</td><td width=50>%s</td><td width=70>%s</td><td width=850>' % (udet[0],udet[1],udet[2])]
                    spokes = []
                    for seq in udet[3:]:
                        gene = self.gene(seq)
                        spokes.append(self.geneLink(gene))
                        #if seq in dbxref.index('EnsLoci'):
                        #    gene = dbxref.data()[dbxref.index('EnsLoci')[seq][0]]['Gene']
                        #    spokes.append(self.geneLink(gene))
                        #else: spokes.append(seq)
                    spokes.sort()
                    uhtml += ['    %s' % rje.join(spokes,';\n    '),'  </td>']
                ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                uhtml += ['</tr></table>']
            return rje.join(uhtml,'\n')
        except: self.errorLog('UPCHTML (%s) error' % upc); return '<p><i>UPC file error</i></p>'
#########################################################################################################################
    def uniFakeHTML(self):  ### Generate composite HTML-linked DAT file            #774#                            #V2.0
        '''Generate composite HTML-linked DAT file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xref = self.obj['DB'].getTable('DBXRef')
            mdb = self.obj['DB'].getTable('Main')
            odb = self.obj['DB'].getTable('Occ')
            hdir = rje.makePath('%sunifake/' % self.info['HTMLPath'])
            rje.mkDir(self,hdir)
            ### ~ [2] ~ Generate DAT File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx, stot) = (0.0, xref.entryNum())
            for spoke in rje.sortKeys(xref.index('Gene')):
                ## ~ [2a] ~ Generate general attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#UNI','Generating HTML-linked DAT files: %.2f%%' % (sx/stot)); sx += 100.0
                if self.opt['Test'] and spoke[0] != 'Y': continue
                try: spokeseq = xref.data()[spoke]['EnsLoci']
                except:
                    self.errorLog('%s UniFake' % spoke)
                    #self.deBug(xref.data()[spoke])
                    continue
                hfile = '%s%s.html' % (hdir,spokeseq)
                if not self.opt['Force'] and self.checkHTML(hfile): continue
                try: fakeacc = rje.split(spokeseq,'__')[1]
                except: continue
                ## ~ [2b] ~ Compile UniProt data and add Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    shhh = self.log.opt['Silent']
                    self.log.opt['Silent'] = True
                    fakeuni = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=%s' % self.info['UniFake']])
                    fakeuni.readUniProt(clear=True,acclist=[fakeacc],cleardata=False)
                    self.log.opt['Silent'] = shhh
                    dat = fakeuni.list['Entry'][0]
                    #self.deBug(dat.dict)
                    sequence = dat.obj['Sequence'].info['Sequence'][0:]
                except:
                    self.errorLog('Problem making HTML-linked DAT file for %s' % spokeseq)
                    continue
                ## ~ [2c] ~ Add occurrences as UniProt features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if spokeseq in odb.index('Seq'):
                    mapdir = '../map_png/'
                    motdir = '../motifaln_png/'
                    maptit = 'Alignment of motif occurrence'
                    mottit = 'Summary graphic for motif'
                    for occ in odb.entryList(odb.index('Seq')[spokeseq]):
                        mappng = '%s.%s.%s.png' % (mapdir,occ['Dataset'],occ['Seq'],self.pround(rje.atoi(occ['Start_Pos'])))
                        motpng = '%s.%s.%d.%s.png' % (motdir,occ['Dataset'],int(occ['Rank']),rje_slim.slimFromPattern(dat['Pattern']))
                        ftdic = {'Type':'SLIM','Start':rje.atoi(occ['Start_Pos']),'End':rje.atoi(occ['End_Pos']),
                                 'PosLink':'<a href="%s" target="_blank" title="%s">' % (mappng,maptit), 
                                 'Desc':'%s Rank <a href="%s" target="_blank" title="%s">%s</a> SLiM %s (p=%s)' % (occ['Dataset'],motpng,mottit,int(occ['Rank']),occ['Pattern'],rje_slim.expectString(rje.atof(occ['Sig'])))}
                        dat.list['Feature'].append(ftdic)
                ## ~ [2d] ~ Add PPI to comments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for ppi in self.db().getTable('PPI').indexEntries('Hub',spoke):
                    plink = self.geneLink(ppi['Spoke'])
                    dat.dict['Data']['CC'].append('-!- PPI: %s <-> %s' % (plink,ppi['Evidence']))
                ### ~ [3] ~ Save data to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                html = htmlHead(spoke,self.list['StyleSheets'])
                ## ~ [3a] ~ Gene summary as with hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                html += self.seqDetailsHTML(spoke,xref.data()[spoke]) #gene,seqid,dbxref,desc,godata)
                ## ~ [3b] ~ UniFake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                uhtml = ['<H2>UniProt Format</H2>\n<PRE>']
                entry = dat
                eseq = entry.obj['Sequence']
                #self.deBug(entry.dict['Data'])
                ## Standard info ##
                for key in ['ID','AC','DT','DE','GN','OS']:
                    if key in entry.dict['Data']:
                        for rest in entry.dict['Data'][key]:
                            if key == 'ID':
                                id = rje.split(rest)
                                rest = rje.join(['<a href="http://www.uniprot.org/uniprot/%s" title="Link to UniProt entry">%s</a>' % (id[0],id[0])] + id[1:])
                            if key == 'GN':
                                id = rje.matchExp('^(GN\s+Name=)(\S+)(;.*)',rest)
                                if id: rest = rje.join([id[0],self.geneLink(id[1]),id[2]],'')
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
                        fsplit = rje.split(ftxt)
                        for i in range(1,len(fsplit)): fsplit[i] = '%s%s</a>' % (ftdict['PosLink'],ftsplit[i])
                        ftxt = rje.join(fsplit)
                    ftxt += ' %s\n' % self.linkFT(ftdict['Type'],ftdict['Desc'],spokeseq)
                    uhtml.append(ftxt)
                ## Sequence/End ##
                uhtml.append('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % eseq.aaLen())),eseq.aaLen(),rje_sequence.MWt(eseq.info['Sequence'])))
                uniseq = eseq.info['Sequence'][0:]
                while len(uniseq) > 0:
                    uhtml.append('     %s\n' % rje.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                    uniseq = uniseq[60:]
                uhtml.append('//\n</PRE>\n')
                html += rje.join(uhtml,'') + htmlTail()
                open(hfile,'w').write(html)
                self.printLog('#UNI','Generation of %s HTML-linked DAT file complete.' % spoke,screen=False)
            self.printLog('\r#UNI','Generation of HTML-linked DAT files complete.')
        except: self.errorLog('UniFake problem')
#########################################################################################################################
    def linkFT(self,type,desc,seq):  ### Adds hyperlinks to feature text
        '''Adds hyperlinks to feature text.'''
        try:### ~ [1] ~ Hyperlink SLiMs and datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if type == 'SLIM':
                dsplit = rje.split(desc)
                dset = dsplit[0]
                hub = rje.join(rje.split(dset,'.')[:-1],'.')
                pattern = dsplit[4]
                slim = rje_slim.slimFromPattern(pattern)
                if hub in self.db().getTable('PPI').index('Hub'): dsplit[0] = '<A HREF="../gene/%s.html" title="%s Results page">%s</A>' % (hub,hub,dsplit[0])
                elif dset[:4] in ['rseq','rupc']: dsplit[0] = '<A HREF="../%s/%s.html" title="%s Results page">%s</A>' % (dset[:4],dset,dset,dsplit[0])
                else: dsplit[0] = '<A HREF="../domain/%s.html" title="%s Results page">%s</A>' % (hub,hub,dsplit[0])
                dsplit[4] = '<A HREF="../slim/%s.html" title="%s Results page">%s</A>' % (slim,dsplit[4],dsplit[4])
                return rje.join(dsplit)
            if type in ['PFAM','DOMAIN']:
                hub = rje.split(desc)[0]
                if hub in self.db().getTable('PFam').index('Type'):
                    dsplit = ['<A HREF="../domain/%s.html" title="%s Results page">%s</A>' % (hub,hub,hub)] + rje.split(desc)[1:]
                    return rje.join(dsplit)
        except: self.errorLog(rje_zen.Zen().wisdom())
        return desc
#########################################################################################################################
    def tabberHTML(self,id,tablist,level=0):     ### Returns text for Tabber HTML
        '''Returns text for Tabber HTML.'''
        try: return tabberHTML(id,tablist,level)
        except: self.errorLog('TabberHTML Error')
        return '!Error!'
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
                    open(self.info['TitleText'],'w').write('%s\n' % rje.join(['Page','ID','Title'],'\t'))
            tkey = '%s\t%s' % (page,id)
            ### ~ [1] ~ Easy case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if tkey in tdb.data(): return tdb.data()[tkey]['Title']
            ### ~ [2] ~ Add new case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not rje.yesNo('Add new mouseover (title) text for %s|%s?' % (page,id)): return ''
            title = rje.choice('New mouseover (title) text for %s|%s title?' % (page,id),confirm=True)
            tdb.data()[tkey] = {'Page':page,'ID':id,'Title':title}
            open(self.info['TitleText'],'a').write('%s\n' % rje.join([page,id,title],'\t'))
            return title
        except: self.errorLog('TitleText problem (%s|%s)' % (page,id)); return ''
#########################################################################################################################
    ### <6> ### SLiMJIM Methods                                                                                         #
#########################################################################################################################
    def makeMapFas(self,mfile): ### Make *.mapping.fas from GOPHER alignment and occ data
        '''Make *.mapping.fas from GOPHER alignment and occ data.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdir = os.path.dirname(mfile)
            bdir = rje.join(rje.split(mdir,'_')[:2]+['Build/'],'_')
            dataset = rje.replace(os.path.basename(mfile),'.mapping.fas','')
            mdb = self.db().getTable('Main')
            patlist = mdb.dataList(mdb.indexEntries('Dataset',dataset),'Pattern')
            occdb = self.db().getTable('Occ')
            occlist = occdb.indexEntries('Dataset',dataset)
            ## ~ [0a] Setup SLiMFinder Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            basesf = slimfinder.SLiMFinder(self.log,['ResDir=%s' % mdir]+self.cmd_list)
            #sf.setInfo({'BuildPath':bdir})
            #sf = sf.pickleMe(load=True)
            pfile = '%s%s.l5w2o2a1.FreqConsDisComp-4-6FT' % (bdir,dataset)
            #self.deBug('%s.pickle.gz: %s' % (pfile,os.path.exists('%s.pickle.gz' % pfile)))
            sf = self.unpickleMe(pfile,process=True)
            sf.setStat({'WallTime':0,'Extras':1,'Basefile':'%s/%s' % (basesf.info['ResDir'],dataset)})
            sf.setInfo({'ResDir':basesf.info['ResDir']})
            sf.obj['SlimList'] = basesf.obj['SlimList']    # Should take SLiMCalc with it
            sf.setLog(self.log)
            sf.list['SigSlim'] = []
            ### ~ [1] Add Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pattern in patlist:
                slim = rje_slim.slimFromPattern(pattern)
                sf.list['SigSlim'].append(slim)
                sf.addSLiMToList(slim)
            ### ~ [2] Generate alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sf.extraOutput()
        except: self.errorLog('Problem during makeMapFas(%s)' % mfile,quitchoice=False)
        #self.deBug('How did that go?')
#########################################################################################################################
    def sigMappingJIM(self):    ### Generate PNG files from *.mapping.fas files                                     |0.1|
        '''Generate PNG files from *.mapping.fas files.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['MakePNG']: return
            self.dict['Forks'] = {}; self.list['Forks'] = {}
            self.opt['NoForks'] = self.opt['NoForks'] or self.stat['Forks'] < 2
            self.stat['KillTime'] = time.time()
            sigdir = '../2010-06-HumSF09_Data_Archive/HumSF09_Sig/'
            pngdir = '%smap_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            mfiles = glob.glob('%s*.mapping.fas' % sigdir)
            mtot = len(mfiles); mx = 0
            occdb = self.db().getTable('Occ')
            winsize = 30
            self.printLog('#MAP','%s mapping fasta files for PNG generation' % rje.integerString(mtot))
            ### ~ [1] Process each file in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dset in rje.sortKeys(occdb.index('Dataset')):
                if dset[:4] == 'rseq': rsigdir = rje.replace(sigdir,'Sig','RSeq_Sig')
                elif dset[:4] == 'rupc': rsigdir = rje.replace(sigdir,'Sig','RUPC_Sig')
                else: rsigdir = rje.replace(sigdir,'Sig','Real_Sig')
                mfile = '%s%s.mapping.fas' % (rsigdir,dset)
                #if mfile not in mfiles:
                #    mfile = rje.replace(mfile,'.com.','.comppi.')
                #    mfile = rje.replace(mfile,'.ppidom.','.domppi.')
                #    mfile = rje.replace(mfile,'.comdom.','.comdomppi.')
                #if mfile not in mfiles: self.errorLog('Where is %s?!' % mfile); continue
                if not os.path.exists(mfile):
                    targz = rje.replace(mfile,'mapping.fas','l5w2o2a1.FreqConsDisComp-4-6FT.tar.gz')
                    os.system('tar -xzf %s' % targz)
                #if not os.path.exists(mfile): self.makeMapFas(mfile)
                if not os.path.exists(mfile): self.errorLog('Where is %s?!' % mfile); continue
                ## ~ [1a] Load file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                mx += 1; self.printLog('#MAP','%s: %s of %s' % (mfile,rje.integerString(mx),rje.integerString(len(occdb.index('Dataset')))),log=False)
                seqcmd = ['accnr=F','seqnr=F','gnspacc=F','autoload=T','autofilter=F','replacechar=F','seqin=%s' % mfile]
                mseq = rje_seq.SeqList(self.log,self.cmd_list+seqcmd)
                mseqs = mseq.seqs()
                ## ~ [1b] Find spoke seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqsplit = {}; subseq = []; seqid = None
                for seq in mseqs:
                    if seq.info['Name'].find('Motifs') > 0:
                        if subseq: seqsplit[subseq[2].shortName()] = subseq[0:]
                        subseq = [seq]; seqid = seq.shortName()
                    else: subseq.append(seq)
                if subseq: seqsplit[subseq[2].shortName()] = subseq[0:]
                #self.deBug(rje.sortKeys(seqsplit))
                ## ~ [1c] Make Occ PNGs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqocc = {}; ox = 0; oxx = 0; px = []
                for entry in occdb.indexEntries('Dataset',dset):
                    enseq = entry['Seq']; ox += 1
                    pos = rje.atoi(entry['Start_Pos'])
                    basefile = '%s%s.%s' % (pngdir,dset,enseq)
                    ofile = '%s.%s.png' % (basefile,self.pround(pos))#seq.aaLen()))
                    if os.path.exists(ofile) and not self.opt['Force'] or self.pround(pos) in px: continue
                    if enseq not in seqocc: seqocc[enseq] = []
                    seqocc[enseq].append(entry); oxx += 1; px.append(self.pround(pos))
                self.printLog('\r#PNG','Make %d PNGs for %d of %d occ' % (len(px),oxx,ox))
                for enseq in rje.sortKeys(seqocc):
                    if enseq not in seqsplit: self.errorLog('Where is %s %s?' % (mfile,enseq)); continue
                    seq = seqsplit[enseq][2]
                    spoke = seqocc[enseq][0]['Spoke']
                    basefile = '%s%s.%s' % (pngdir,dset,enseq)
                    pseq = seqsplit[enseq][0:]
                    pseq[0].info['R'] = '%s' % dset
                    pseq[1].info['R'] = 'Masking'
                    q = 2   # Index of Query
                    pseq[q].info['R'] = spoke
                    for seq in pseq[(q+1):]: seq.info['R'] = rje.getFromDict(rje_ensembl.enscom,seq.info['SpecCode'])
                    ## ~ [3b] ~ Setup new SeqList, strip Query gaps, calculate RelCons ~~~~~~~~~~~~~~~~ ##
                    seqfile = '%s.aln.tdt' % basefile
                    if os.path.exists(seqfile): os.unlink(seqfile)
                    scmd = ['seqnr=F','accnr=F','replacechar=F']
                    rseq = rje_seq.SeqList(self.log,['minregion=3']+self.cmd_list+scmd+['autoload=F'])  #!# Check minregion #!#
                    rseq.seq = pseq
                    rseq.obj['QuerySeq'] = pseq[q]
                    rseq.tidyQueryGaps()
                    rseq.saveR(rseq.seq,seqfile,name='R')
                    #self.deBug(rseq.obj['QuerySeq'].info['Sequence'])
                    rseq.seq = pseq[q:]
                    relfile = '%s.rel.tdt' % basefile
                    if os.path.exists(relfile): os.unlink(relfile)
                    rseq.opt['MakePNG'] = False
                    rseq.relCons(relfile)
                    ## ~ [3c] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for entry in seqocc[enseq]:
                        pos = rje.atoi(entry['Start_Pos'])
                        ofile = '%s.%s.png' % (basefile,rje.preZero(self.pround(pos),10000))#seq.aaLen()))
                        if os.path.exists(ofile) and not self.opt['Force']: continue
                        rcmd = '%s --no-restore --no-save --args "occaln" "%s" "%s"' % (self.info['RPath'],basefile,self.pround(pos))    #!# Update R code #!#
                        rslimjim = '%s/rje.r' % rje.makePath(self.info['Path'])
                        rcmd += ' < "%s" > "%s.%s.r.tmp.txt"' % (rslimjim,basefile,self.pround(pos))
                        problems = self.rCall(rcmd,basefile)
                        if os.path.exists(ofile) and not self.opt['Test'] and not self.opt['Iridis'] and self.opt['NoForks']: 
                            for ext in ['%s.r.tmp.txt' % self.pround(pos)]:
                                if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
                        #self.deBug('...')
                    if not self.opt['Test'] and not self.opt['Iridis'] and self.opt['NoForks']:  
                        for ext in ['aln.tdt','rel.tdt']:
                            if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
                    #self.deBug('...')
                #x#if self.opt['Test']: break
            if not self.opt['NoForks']:
                forks = self.dict['Forks']
                while len(forks) >= self.stat['Forks']:
                    ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    time.sleep(1)       # Sleep for 1s 
                    forklist = self._activeForks(forks.keys())
                    if len(forklist) != len(forks):
                        self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                        for pid in rje.sortKeys(forks):    # Go through current forks
                            if pid not in forklist:
                                basefile = forks.pop(pid)
                                pngx = len(glob.glob('%s*png' % basefile))
                                self.printLog('#PNG','%d PNG files made for %s' % (pngx,basefile))
                                #!#if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
                                self.stat['KillTime'] = time.time()  # Reset killtime - still doing stuff
                    ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if time.time() - self.stat['KillTime'] > self.stat['KillForks']:
                        self.printLog('#KILL','%d seconds of main thread inactivity. %d forks still active!' % (self.stat['KillForks'],len(forks)))
                        for fork in forks: self.printLog('#FORK','Fork %s, PID %d still Active!' % (forks[fork],fork))
                        if self.stat['Interactive'] < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                        elif rje.yesNo('Kill hanging forks?'):
                            for fork in forks:
                                self.printLog('#KILL','Killing Fork %s, PID %d.' % (forks[fork],fork))
                                os.system('kill %d' % fork)
                        else: self.stat['KillTime'] = time.time()
            self.printLog('#MAP','%s of %s mapping files done' % (rje.integerString(mx),rje.integerString(mtot)))
        except: self.errorLog('Problem during sigMappingJIM()')
#########################################################################################################################
    def rCall(self,rcmd,basefile):   ### Performs given R call for PNG generation
        '''Performs given R call for PNG generation.'''
        ### ~ [0] ~ Check for PNG and skip if found and not Force ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #if glob.glob('%s*png' % basefile) and not self.opt['Force']:
        #    return self.printLog('#PNG','%s*png exists - skipping (Force=F)' % basefile)
        ### ~ [1] ~ "Standard" R Call without using IRIDIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.opt['Iridis'] = False  #!# Make this part of main RJE Object code #!#
        #if not self.opt['Iridis']: rcmd = rcmd + ' 2>&1'
        self.printLog('#RJIM',rcmd)
        if not self.opt['Iridis'] and self.opt['NoForks']:
            problems = os.popen(rcmd).read()
            if problems: self.errorLog(problems,printerror=False)
            pngx = len(glob.glob('%s*png' % basefile))
            self.printLog('#PNG','%d PNG files made for %s' % (pngx,basefile))
            if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
            return problems
        elif not self.opt['Iridis']:
            forks = self.dict['Forks']
            while len(forks) >= self.stat['Forks']:
                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(1)       # Sleep for 1s 
                forklist = self._activeForks(forks.keys())
                if len(forklist) != len(forks):
                    self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    for pid in rje.sortKeys(forks):    # Go through current forks
                        if pid not in forklist:
                            basefile = forks.pop(pid)
                            pngx = len(glob.glob('%s*png' % basefile))
                            self.printLog('#PNG','%d PNG files made for %s' % (pngx,basefile))
                            #!#if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
                            self.stat['KillTime'] = time.time()  # Reset killtime - still doing stuff
                ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if time.time() - self.stat['KillTime'] > self.stat['KillForks']:
                    self.printLog('#KILL','%d seconds of main thread inactivity. %d forks still active!' % (self.stat['KillForks'],len(forks)))
                    for fork in forks: self.printLog('#FORK','Fork %s, PID %d still Active!' % (forks[fork],fork))
                    if self.stat['Interactive'] < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                    elif rje.yesNo('Kill hanging forks?'):
                        for fork in forks:
                            self.printLog('#KILL','Killing Fork %s, PID %d.' % (forks[fork],fork))
                            os.system('kill %d' % fork)
                    else: self.stat['KillTime'] = time.time()
            ### Fork out RCMD ###
            newpid = os.fork() 
            if newpid == 0: # child
                self.opt['Child'] = True
                problems = os.popen(rcmd).read()
                if problems: self.errorLog(problems,printerror=False)            
                os._exit(0)    # Exit process 
            elif newpid == -1: self.errorLog('Problem forking %s.' % id,printerror=False)  
            else: forks[newpid] = basefile            
            return False
        ### ~ [2] ~ Complex R call using IRIDIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [2a] ~ Wait for free node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if len(self.dict['Running']) >= (self.nprocs()-1): host_id = self.runJobs()   # Wait for free node
        else:
            for host_id in range(1,self.nprocs()):
                if host_id not in self.dict['Running']: break
        try: self.list['Hosts'][host_id]
        except: self.errorLog('Problem with host_id "%s"' % host_id)
        ## ~ [2b] ~ Call R on remote node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        try:
            jdict = self.dict['Running'][host_id] = {'Basefile':basefile}      # Setup dictionary to fill
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as subjob %s on `hostname` ; module load R/2.7.1/gcc-3.4.4 ;' % (basefile,host_id)
            job = '%s %s'  % (initial_cmds, rcmd)
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            self.printLog('#RSH',rsh)
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#JOB','Running job as %s: %d running' % (cpid,len(self.dict['Running'])))
                return []
            else:                   # child process
                os.system(rsh)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('SLiMJIM IRIDIS rCall error')
#########################################################################################################################
    def gene(self,prot):    ### Return Gene for Protein
        pdb = self.db().getTable('PPI')
        if prot not in self.dict['GeneMap']:
            try:
                self.dict['GeneMap'][prot] = pdb.dataList(pdb.indexEntries('SpokeSeq',prot),'Spoke')[0]
                if self.dict['GeneMap'][prot] in ['','-']: raise ValueError
            except:
                dbxref = self.db().getTable('DBXRef')
                if prot not in dbxref.index('EnsLoci'):
                    details = rje.split(prot,'_')
                    if details[0] != details[0].upper(): self.dict['GeneMap'][prot] = details[-1]
                    else: self.dict['GeneMap'][prot] = rje.join(details[:2],'_')
                else: self.dict['GeneMap'][prot] = dbxref.dataList(dbxref.indexEntries('EnsLoci',prot),'Gene')[0]
            if self.dict['GeneMap'][prot] in ['-','']: self.dict['GeneMap'][prot] = prot
            for x in '!"%&\'()*+,/:;<=>?@[\\]': self.dict['GeneMap'][prot] = rje.replace(self.dict['GeneMap'][prot],x,'')
        elif self.dict['GeneMap'][prot] in ['-','']:
            self.dict['GeneMap'][prot] = prot                  
            for x in '!"%&\'()*+,/:;<=>?@[\\]': self.dict['GeneMap'][prot] = rje.replace(self.dict['GeneMap'][prot],x,'')
        return self.dict['GeneMap'][prot]
#########################################################################################################################
    def hubInteractomeSVG(self,dataset):    ### Generates Hub interactome tree & ppi visualisation
        '''Generates Hub interactome tree & ppi visualisation.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hobj = self.obj['HAPPI']
            svg = self.obj['SVG']
            pngdir = '%sinteractome_svg/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            basefile = rje.makePath('%s%s.interactome' % (pngdir,dataset),wholepath=True)
            if os.path.exists('%s.svg' % basefile) and not self.opt['Force']: return
            sigdir = '../2010-06-HumSF09_Data_Archive/HumSF09_Sig/'
            hub = dataset
            if dataset[:4] == 'rseq': rsigdir = rje.replace(sigdir,'Sig','RSeq_Sig')
            elif dataset[:4] == 'rupc': rsigdir = rje.replace(sigdir,'Sig','RUPC_Sig')
            else: rsigdir = rje.replace(sigdir,'Sig','Real_Sig'); hub = rje.split(dataset,'.')[0]
            if dataset not in self.db().getTable('Occ').index('Dataset'): rsigdir = rje.replace(rsigdir,'Sig','NonSig')
            disfile = '%s%s.dis.tdt' % (rsigdir,dataset)        # Make tree from this
            ## ~ [1a] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s%s.upc' % (rsigdir,dataset)            # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = rje.split(uline)[3:]
                upc.sort()
                uplist.append(rje.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = rje.split(upc)
                reorder += uprots     #!# Then reorder sequences! #!#
                if len(uprots) > 1:
                    for u in uprots: upid[self.gene(u)] = upx
                    upx += 1
                else: upid[self.gene(uprots[0])] = 0

            ### ~ [2] ~ Load distance matrix and make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nsf_file = '%s.nsf' % basefile
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            try:
                if not os.path.exists(nsf_file) or self.opt['Force']: raise ValueError
                upgma.loadTree(nsf_file,type='nsf',postprocess=False)
            except:                
                disdat = rje.dataDict(self,disfile,['SEQ'])
                dismat = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
                dismat.opt['ProgLog'] = hobj.opt['NoForks']
                for obj1 in reorder:
                    for obj2 in reorder:
                        dis = (rje.atof(disdat[obj1][obj2]) + rje.atof(disdat[obj2][obj1])) / 2.0
                        dismat.addDis(self.gene(obj1),self.gene(obj2),dis)
                nsftree = dismat.upgma()
                open(nsf_file,'w').write(nsftree)
                upgma.buildTree(nsftree,type='nsf',postprocess=False)
            svgdata = upgma.dictTree(seqname='short')
            for n in rje.sortKeys(svgdata):
                d = svgdata[n]
                if d['name'] in upid: d['col'] = upid[d['name']] + 1
            svgheight = max(1000,100 * int(0.081*len(reorder)+0.5))
            svgtext = svg.svgTree(basefile,svgdata,treesplit=0.75,font=8,maxfont=20,width=400,height=svgheight,save=False,xoffset=0,yoffset=0)

            ### ~ [3] ~ Add PPI links between spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pobj = rje_ppi.PPI(self.log,self.cmd_list+['backups=F'])
            pobj.obj = self.obj['PPI'].obj
            pobj.dict = self.obj['PPI'].dict
            pobj.opt['ProgLog'] = hobj.opt['NoForks']
            names = []; G = {}; pdb = self.db().getTable('PPI')
            for prot in reorder:
                if prot in pdb.index('SpokeSeq'): names.append(self.gene(prot))
            if hub in names: names.remove(hub)    # Will interact with everybody!
            for p1 in names:
                G[p1] = []
                for p2 in names:
                    try: ppi = p2 in pdb.dataList(pdb.indexEntries('Hub',p1),'Spoke')    #self.dict['PPI'][p1]  #False
                    except: ppi = False
                    if ppi: G[p1].append(p2)
            ncol = {}
            for p1 in names: ncol[p1] = upid[p1] + 1
            nfile = '%s.npos.tdt' % basefile
            npos = {}
            if hobj.opt['UsePos'] and os.path.exists(nfile): npos = pobj.loadNPos(nfile)
            if rje.sortKeys(npos) != rje.sortKeys(G): npos = pobj.rjeSpringLayout(G); pobj.saveNPos(npos,nfile)
            elif hobj.opt['UpdatePos']:
                npos = rje_ppi.rjeSpringLayout(G,damping=pobj.stat['Damping'],callobj=pobj,walltime=pobj.stat['Walltime'],nudge=pobj.stat['NudgeCyc'],prepos=npos)
                if os.path.exists(nfile): os.unlink(nfile)
                pobj.saveNPos(npos,nfile)
            self.deBug(ncol)
            svgwidth = 1000 + ((svgheight - 1000) / 2)
            svgtext += svg.networkPlot(npos,G,ncol,width=svgwidth,height=svgheight,xoffset=400)
            svg.svgFile(svgtext,'%s.svg' % basefile,width=svgwidth+400,height=svgheight)
            self.printLog('#SVG','Saved %s.svg' % basefile)
            return True
        except: self.errorLog(rje_zen.Zen().wisdom()); return False
#########################################################################################################################
    def slimJimHub(self,dataset):   ### Peforms the SLiMJim visualisation for a given Hub dataset               #V2.0
        '''Peforms the SLiMJim visualisation for a given Hub dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SVG']: return self.hubInteractomeSVG(dataset)
            pngdir = '%sinteractome_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            basefile = rje.makePath('%s%s.interactome' % (pngdir,dataset),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            sigdir = '../2010-06-HumSF09_Data_Archive/HumSF09_Sig/'
            hub = dataset
            if dataset[:4] == 'rseq': rsigdir = rje.replace(sigdir,'Sig','RSeq_Sig')
            elif dataset[:4] == 'rupc': rsigdir = rje.replace(sigdir,'Sig','RUPC_Sig')
            else: rsigdir = rje.replace(sigdir,'Sig','Real_Sig'); hub = rje.split(dataset,'.')[0]
            if dataset not in self.db().getTable('Occ').index('Dataset'): rsigdir = rje.replace(rsigdir,'Sig','NonSig')
            disfile = '%s%s.dis.tdt' % (rsigdir,dataset)        # Make tree from this
            ## ~ [1a] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s%s.upc' % (rsigdir,dataset)            # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = rje.split(uline)[3:]
                upc.sort()
                uplist.append(rje.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = rje.split(upc)
                reorder += uprots     #!# Then reorder sequences! #!#
                if len(uprots) > 1:
                    for u in uprots: upid[self.gene(u)] = upx
                    upx += 1
                else: upid[self.gene(uprots[0])] = 0

            ### ~ [2] ~ Load distance matrix and make Tree/Heatmap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disdat = rje.dataDict(self,disfile,['SEQ'])
            dismat = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
            for obj1 in reorder:
                for obj2 in reorder:
                    dis = (rje.atof(disdat[obj1][obj2]) + rje.atof(disdat[obj2][obj1])) / 2.0
                    dismat.addDis(self.gene(obj1),self.gene(obj2),dis)
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            nsftree = dismat.upgma()
            open('%s.nsf' % basefile,'w').write(nsftree)
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            upgma.rTree('%s.tree.csv' % basefile,seqname='short')
            reduced = upgma._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            dismat.info['Name'] = '%s interactome GABLAM' % hub
            dismat.saveMatrix(reduced,basefile+'.heatmap.tdt',delimit='\t')

            ### ~ [3] ~ Add PPI links between spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            names = []; pdb = self.db().getTable('PPI')
            for prot in reorder:
                if prot in pdb.index('SpokeSeq'): names.append(self.gene(prot))
            if hub in names: names.remove(hub)    # Will interact with everybody!
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for p1 in names:
                datadict = {}
                for p2 in names:
                    try: ppi = p2 in pdb.dataList(pdb.indexEntries('Hub',p1),'Spoke')    #self.dict['PPI'][p1]  #False
                    except: ppi = False
                    if ppi: datadict[p2] = 1
                    else: datadict[p2] = 0
                rje.delimitedFileOutput(self,pfile,names,datadict=datadict)
            datadict = {}
            for p1 in names: datadict[p1] = upid[p1]
            rje.delimitedFileOutput(self,pfile,names,datadict=datadict)

            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "interactome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/rje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <6...> ### OLD JIM Methods                                                                                      #
#########################################################################################################################
    def motifInteractomeSVG(self,pattern):  ### Generates interactome SVG for hubs and spokes returning pattern
        '''Generates interactome SVG for hubs and spokes returning pattern.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if pattern == '-': return
            pngdir = '%snested_svg/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            #?data = self.dict['Data']
            slim = rje_slim.slimFromPattern(pattern)
            basefile = rje.makePath('%s%s' % (pngdir,slim),wholepath=True)
            if os.path.exists('%s.svg' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            ### ~ [2] ~ Generate nested network PPI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            odb = self.db().getTable('Occ')
            ddb = self.db().getTable('DomPPI')
            pdb = self.db().getTable('PPI')
            hubs = odb.dataList(odb.indexEntries('Pattern',pattern),'Dataset')  # Replace with datasets in case domain and gene have same name?
            spokes = odb.dataList(odb.indexEntries('Pattern',pattern),'Spoke')
            for seq in odb.dataList(odb.indexEntries('Pattern',pattern),'Seq'):
                if self.gene(seq) not in spokes: spokes.append(self.gene(seq))
            if '-' in spokes: spokes.remove('-')
            spokes.sort()
            hublist = []; ncol = {}
            for i in range(len(hubs)):
                if hubs[i][:4] not in ['rseq','rupc']: hubs[i] = rje.split(hubs[i],'.')[0]
                if hubs[i] in pdb.index('Hub'):
                    hublist.append(hubs[i])  # Gene dataset
                    if hubs[i] in spokes: ncol[hublist[-1]] = rje_svg.soton[18]
                    else: ncol[hublist[-1]] = rje_svg.soton[4]
                elif hubs[i][:4] in ['rseq','rupc']:
                    hublist.append(rje.join(rje.split(hubs[i],'_')[:2],'_'))
                    if hubs[i][:4] in ['rseq']: ncol[hublist[-1]] = rje_svg.soton[7]
                    else: ncol[hublist[-1]] = rje_svg.soton[10]
                    try: self.deBug('%s >> %s: %s' % (hubs[i],hublist[-1],self.obj['PPI'].ppi()[hublist[-1]]))
                    except: self.deBug('%s >> %s: No PPI' % (hubs[i],hublist[-1]))
                else:
                    hublist.append('d:%s' % hubs[i])     # Domains dataset
                    ncol[hublist[-1]] = rje_svg.soton[2]
                    try: self.deBug('%s >> %s: %s' % (hubs[i],hublist[-1],self.obj['PPI'].ppi()[hublist[-1]]))
                    except: self.deBug('%s >> %s: No PPI' % (hubs[i],hublist[-1]))
            for spoke in spokes:
                if spoke not in ncol: ncol[spoke] = rje_svg.soton[5]
            ppistyle = 'style="stroke:"#9BA3A6";stroke-width:2"'
            domstyle = 'style="stroke:#007275;stroke-width:2"'
            ranstyle = 'style="stroke:#A67891;stroke-width:2"'
            sfstyle = 'style="stroke:rgb(0,0,0);stroke-width:2"'
            ilist = hublist + spokes

            ### HAPPI Code, adapted ###            
            if not ilist: return False
            hobj = self.obj['HAPPI']
            pobj = rje_ppi.PPI(self.log,self.cmd_list+['backups=F'])
            pobj.obj = self.obj['PPI'].obj
            pobj.dict = self.obj['PPI'].dict
            pobj.opt['ProgLog'] = hobj.opt['NoForks']

            #if len(ilist) > self.stat['PNGMax']:
            #    self.printLog('#MAX','Too many genes for %s.png (%s genes)' % (basefile,rje.iStr(len(ilist))))
            #    G = {'%s %s Genes' % (rje.iStr(len(ilist)),id):[]}
            #else:
            self.printLog('#SVG','Making %s.svg (%s genes)' % (basefile,rje.iStr(len(ilist))))
            self.deBug(ilist)
            G = rje_ppi.subGraph(pobj.ppi(),ilist)      # Reduced PPI Graph
            self.deBug(G)
            self.deBug(ncol)
            if not G: return False
            nfile = '%s.npos.tdt' % basefile
            npos = {}
            if hobj.opt['UsePos'] and os.path.exists(nfile): npos = pobj.loadNPos(nfile)
            if rje.sortKeys(npos) != rje.sortKeys(G): npos = pobj.rjeSpringLayout(G); pobj.saveNPos(npos,nfile)
            elif hobj.opt['UpdatePos']:
                npos = rje_ppi.rjeSpringLayout(G,damping=pobj.stat['Damping'],callobj=pobj,walltime=pobj.stat['Walltime'],nudge=pobj.stat['NudgeCyc'],prepos=npos)
                if os.path.exists(nfile): os.unlink(nfile)
                pobj.saveNPos(npos,nfile)
            newG = {}
            for vi in G:
                newG[vi] = {}
                for vj in G[vi]:
                    newG[vi][vj] = ppistyle 
                    if pobj.dict['PPI'][vi][vj] == 'Domain': newG[vi][vj] = domstyle
                    elif pobj.dict['PPI'][vi][vj] == 'Random': newG[vi][vj] = ranstyle
                    elif vi in hubs or vj in hubs: newG[vi][vj] = sfstyle
            G = newG                                        
            #if self.opt['XGMML']:
            #    xgmml = pobj.ppiToXGMML(G,id)
            #    xgmml.dict['NodePos'] = npos
            #    xgmml.saveXGMML('%s.xgmml' % basefile)
            #    self.gZip('%s.xgmml' % basefile)
            self.deBug(npos); self.deBug(G)
            pobj.saveSVG(npos,basefile,G,width=1600,ntype='ellipse',backups=False,nodecol=ncol)
            return True
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimMotif(self,pattern):     ### Performs SLiMJIM visualisation for a given motif                        #V2.0
        '''Performs SLiMJIM visualisation for a given motif.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if pattern == '-': return
            if self.opt['SVG']: return self.motifInteractomeSVG(pattern)
            pngdir = '%snested_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            #?data = self.dict['Data']
            slim = rje_slim.slimFromPattern(pattern)
            basefile = rje.makePath('%s%s' % (pngdir,slim),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            ### ~ [2] ~ Generate nested network PPI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            odb = self.db().getTable('Occ')
            ddb = self.db().getTable('DomPPI')
            pdb = self.db().getTable('PPI')
            hubs = odb.dataList(odb.indexEntries('Pattern',pattern),'Dataset')  # Replace with datasets in case domain and gene have same name?
            spokes = odb.dataList(odb.indexEntries('Pattern',pattern),'Spoke')
            remake = False
            for seq in odb.dataList(odb.indexEntries('Pattern',pattern),'Seq'):
                if self.gene(seq) not in spokes: spokes.append(self.gene(seq)); remake = True
            if '-' in spokes: spokes.remove('-'); remake = True
            spokes.sort()
            nfile = '%s.nested.tdt' % basefile
            if os.path.exists(nfile): os.unlink(nfile)
            headers = ['Spoke']
            for i in range(len(hubs)):
                if hubs[i][:4] not in ['rseq','rupc']: hubs[i] = rje.split(hubs[i],'.')[0]
            hubs = rje.sortUnique(hubs)
            for hub in hubs:
                if hub in pdb.index('Hub'): headers.append(hub)  # Gene dataset
            for hub in hubs:
                if hub in pdb.index('Hub') or hub[:4] in ['rseq','rupc']: continue
                headers.append('d%s' % hub)     # Domains dataset
            for hub in hubs:
                if hub[:4] in ['rseq','rupc']: headers.append(rje.join(rje.split(hub,'_')[:2],'_'))
            rje.delimitedFileOutput(self,nfile,headers)
            ndata = {}; sortdict = {}
            for spoke in spokes:
                datadict = {'Spoke':spoke}
                for hub in hubs:
                    if hub in pdb.index('Hub'): h = hub  # Gene dataset
                    elif hub[:4] in ['rseq','rupc']: h = rje.join(rje.split(hub,'_')[:2],'_')
                    else: h = 'd%s' % hub               # Domains dataset
                    datadict[h] = 0
                    if spoke in odb.dataList(odb.indexEntries('Hub',hub),'Spoke'): datadict[h] = 1
                    elif spoke in odb.dataList(odb.indexEntries('Dataset',hub),'Spoke'): datadict[h] = 1
                    #elif spoke in odb.dataList(odb.indexEntries('Hub',hub),'Seq'): datadict[h] = 1
                    #elif spoke in odb.dataList(odb.indexEntries('Dataset',hub),'Seq'): datadict[h] = 1
                    elif hub[:4] not in ['rseq','rupc']:    # Random dataset
                        if hub in pdb.index('Hub'):                 # Gene dataset
                            if spoke in pdb.dataList(pdb.indexEntries('Hub',hub),'Spoke'): datadict[h] = 1
                            #elif spoke in pdb.dataList(pdb.indexEntries('Hub',hub),'Seq'): datadict[h] = 1
                        if hub in ddb.index('Domain'):
                            if spoke in ddb.dataList(ddb.indexEntries('Domain',hub),'Spoke'): datadict[h] = 1
                            #elif spoke in pdb.dataList(pdb.indexEntries('Hub',hub),'Seq'): datadict[h] = 1
                    if datadict[h] and 'Sort' not in datadict: datadict['Sort'] = h
                    elif datadict[h] and headers.index(h) < headers.index(datadict['Sort']): datadict['Sort'] = h
                try:
                    skey = datadict.pop('Sort')
                    if skey not in sortdict: sortdict[skey] = []
                    sortdict[skey].append(spoke)
                except: self.errorLog('Problem with %s spoke %s' % (pattern,spoke))
                ndata[spoke] = datadict
            spokelist = []; altx = 0
            for h in headers[1:]:
                if h not in sortdict: continue
                for spoke in sortdict[h]: spokelist.append(spoke)
                if h == headers[1] and len(hubs) > 1: altx = len(sortdict[h]) / 2
            if altx: spokelist = spokelist[altx:] + spokelist[:altx]
            for spoke in spokelist: rje.delimitedFileOutput(self,nfile,headers,datadict=ndata[spoke])
            ### ~ [3] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "nested" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/rje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['nested.tdt','r.tmp.txt']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
            #self.deBug('%s.png' % basefile) 
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimHubMotif(self,sdata):    ### Peforms the SLiMJim visualisation for a given hub-motif data dictionary #V2.0
        '''Peforms the SLiMJim visualisation for a given hub-motif data dictionary.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Still needs standardising use of hub, spoke, dataset, spokeseq, pattern, slim. #!#
            if sdata['Pattern'] == '-': return
            pngdir = '%smotifaln_png/' % self.info['HTMLPath']
            rje.mkDir(self,pngdir)
            dset = dataset = sdata['Dataset']
            basefile = rje.makePath('%s%s.%s.%s' % (pngdir,dataset,sdata['Rank'],rje_slim.slimFromPattern(sdata['Pattern'])),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            sigdir = '../2010-06-HumSF09_Data_Archive/HumSF09_Sig/'
            if dset[:4] == 'rseq': rsigdir = rje.replace(sigdir,'Sig','RSeq_Sig')
            elif dset[:4] == 'rupc': rsigdir = rje.replace(sigdir,'Sig','RUPC_Sig')
            else: rsigdir = rje.replace(sigdir,'Sig','Real_Sig')
            if dataset not in self.db().getTable('Occ').index('Dataset'): rsigdir = rje.replace(rsigdir,'Sig','NonSig')
            gene = sdata['Hub'] #rje.replace(sdata['Dataset'],self.info['Suffix'],'')
            slim = sdata['Pattern']
            rank = rje.atoi(sdata['Rank'])
            basefile = rje.makePath('%s%s.%s.%s' % (pngdir,dataset,rank,rje_slim.slimFromPattern(slim)),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            disfile = '%s%s.dis.tdt' % (rsigdir,dataset)        # Make tree from this
            occfile = '%s%s.occ.csv' % (rsigdir,dataset)        # Get hit data from this
            motifaln = '%s%s.motifaln.fas' % (rsigdir,dataset)  # Or this?!
            ## ~ [1a] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s%s.upc' % (rsigdir,dataset)            # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = rje.split(uline)[3:]
                upc.sort()
                uplist.append(rje.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = rje.split(upc)
                reorder += uprots     #!# Then reorder sequences! #!#
                if len(uprots) > 1:
                    for u in uprots: upid[self.gene(u)] = upx
                    upx += 1
                else: upid[self.gene(uprots[0])] = 0

            ### ~ [2] ~ Load and tidy alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = ['seqin=%s' % motifaln,'autoload=T','seqnr=F','accnr=F','replacechar=F','maxx=0','maxgap=0']
            mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            #i = rank
            #findslim = rje.replace(rje.replace(slim,'^',''),'$','')
            slimsplit = rje.split(mseq.seq[0].info['Sequence'],'-XXXXXXXXXX-')
            #while (rje.replace(slimsplit[i-1][10:],'-','') not in [slim,findslim]): i += 1   # What does this do?! #
            for seq in mseq.seq[0:]:
                seq.info['Sequence'] = rje.split(seq.info['Sequence'],'-XXXXXXXXXX-')[rank-1]
                if not rje.replace(seq.info['Sequence'],'-',''): mseq.seq.remove(seq)
            self.printLog('#SEQ','%d seq from motifaln with Rank %d motif %s' % (mseq.seqNum()-1,rank,slim))
            ## ~ [2a] ~ Reorder sequences to group by UPC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ordseq = mseq.seq[:1]
            mseq.seq[0].info['Gene'] = mseq.seq[0].shortName()
            for seq in mseq.seq[1:]:
                seq.info['Gene'] = self.gene(seq.shortName())
                i = 1
                try:
                    while i < len(ordseq):
                        if reorder.index(seq.shortName()) < reorder.index(ordseq[i].shortName()): ordseq.insert(i,seq); break
                        i += 1
                except: self.errorLog('%s vs %s >> %s' % (seq.shortName(),ordseq[i].shortName(),reorder))
                if seq not in ordseq: ordseq.append(seq)
            mseq.seq = ordseq[0:]
            ## ~ [2b] ~ Modify Motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            msequence = mseq.seq[0].info['Sequence']
            newm = []
            r = 0
            while r < len(msequence):
                if msequence[r] in ['^','$']: r += 1; continue
                elif msequence[r] != '[': newm.append(msequence[r]); r += 1
                else:
                    amb = rje.matchExp('^\[([A-Z]+)\]',msequence[r:])[0]
                    newm.append('+')    #(amb)
                    r += len(amb)+2
            while len(newm) < len(msequence): newm.append('-')
            mseq.seq[0].info['Sequence'] = rje.join(newm,'')
            ## ~ [2c] ~ Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.aln.tdt' % basefile): os.unlink('%s.aln.tdt' % basefile)
            mseq.saveR(seqfile='%s.aln.tdt' % basefile,name='Gene')

            ### ~ [3] ~ Load full sequences and make profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdbdir = rje.join(rje.split(rsigdir,'_')[:2] + ['SLiMDB'],'_')
            seqfile = '%s/%s.slimdb' % (sdbdir,dataset)         # Get extra seq data from this
            scmd = ['seqin=%s' % seqfile,'autoload=T','replacechar=F','maxx=0','maxgap=0']
            dseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            seqdic = dseq.seqNameDic()
            occdic = self.loadOccData(occfile,['Dataset','Rank','Pattern','Seq','Start_Pos','End_Pos'])
            occdic.pop('Headers')
            #!# Add UPC weighting of profile? #!#
            ## ~ [3a] ~ Make a list of profile sequences +/- 5 of motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            proseq = []
            for okey in rje.sortKeys(occdic):
                odata = occdic.pop(okey)
                if odata['Pattern'] != slim or int(sdata['Rank']) != int(odata['Rank']): continue
                seq = seqdic[odata['Seq']]
                oseq = seq.info['Sequence'][max(0,rje.atoi(odata['Start_Pos'])-6):rje.atoi(odata['End_Pos'])+5]
                c = max(0,6 - rje.atoi(odata['Start_Pos']))
                oseq = '-' * c + oseq
                n = max(0,rje.atoi(odata['End_Pos'])+5-seq.aaLen())
                oseq = oseq + '-' * n
                proseq.append(oseq)
            #self.deBug(proseq)
            ### ~ [3b] ~ Convert to profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            poslist = ['-5','-4','-3','-2','-1']
            for pos in rje_slim.prestoFromCode(rje_slim.slimFromPattern(slim)):
                if pos in ['$','^']: continue
                elif len(pos) > 1: poslist.append('[%s]' % pos)
                elif pos == 'X': poslist.append('x')
                else: poslist.append(pos)
            poslist += ['1','2','3','4','5']
            ignorepos = ['-5','-4','-3','-2','-1','x','1','2','3','4','5','$','^']
            outfile = basefile + '.profile.tdt'
            if os.path.exists(outfile): os.unlink(outfile)
            ### ~ [3c] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = ['Pos'] + rje_seq.alph_protx[:-1]
            rje.delimitedFileOutput(self,outfile,headers,'\t')
            #self.deBug(poslist)
            try:
                for r in range(len(poslist)):
                    datadict = {}
                    for a in headers[1:]: datadict[a] = 0.0
                    for pseq in proseq:
                        while len(pseq) < len(poslist) and poslist[r] not in ignorepos and poslist[r].find(pseq[r]) < 0: pseq = pseq[:r] + '-' + pseq[r:]
                        try:
                            a = pseq[r]
                            datadict[a] += 1
                        except:
                            self.deBug('%s %s passed' % (pseq,r))
                            pass    # Not counting Xs
                    datadict['Pos'] = poslist[r]
                    rje.delimitedFileOutput(self,outfile,headers,'\t',datadict)
            except:
                self.errorLog('Problem during profile output for %s' % os.path.basename(basefile))
                for pseq in proseq: print pseq
                
            ### ~ [4] ~ Load distance matrix and make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disdat = rje.dataDict(self,disfile,['SEQ'])
            dismat = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
            reduced = []    # Reduced sequence set, ordered according to UPC
            for seq1 in mseq.seq[1:]:
                obj1 = seq1.shortName()
                reduced.append(obj1)
                for seq2 in mseq.seq[1:]:
                    obj2 = seq2.shortName()   
                    dis = (rje.atof(disdat[obj1][obj2]) + rje.atof(disdat[obj2][obj1])) / 2.0
                    dismat.addDis(seq1.info['Gene'],seq2.info['Gene'],dis)
            upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
            nsftree = dismat.upgma()
            open('%s.nsf' % basefile,'w').write(nsftree)
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            upgma.rTree('%s.tree.csv' % basefile,seqname='short')
            reduced = upgma._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            dismat.saveMatrix(reduced,basefile+'.heatmap.tdt',delimit='\t')

            ### ~ [5] ~ Add PPI links between spokes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #names = reduced[0:]
            names = []; pdb = self.db().getTable('PPI')
            for prot in reduced: names.append(self.gene(prot))
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for hub in names:
                datadict = {}
                for spoke in names:
                    try: ppi = spoke in pdb.dataList(pdb.indexEntries('Hub',hub),'Spoke') #self.dict['PPI'][hub]  #False
                    except: ppi = False
                    #for g1 in self.dict['SeqMap'][hub]:
                    #    for g2 in self.dict['SeqMap'][spoke]:
                    #        if g2 in self.dict['PPI'][g1]: ppi = True
                    if ppi: datadict[spoke] = 1
                    else: datadict[spoke] = 0
                rje.delimitedFileOutput(self,pfile,names,datadict=datadict)
            datadict = {}
            for spoke in names: datadict[spoke] = upid[spoke]
            rje.delimitedFileOutput(self,pfile,names,datadict=datadict)

            ### ~ [6] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "motif" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%s/rje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            #self.deBug('Done')
            ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','aln.tdt','profile.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadOccData(self,occfile,mainkeys): ### Loads and resorts OccData from OccFile
        '''Loads and resorts OccData from OccFile.'''
        delimit = rje.delimitFromExt(filename=occfile)
        return self.sortOccDic(rje.dataDict(self,occfile,mainkeys,'All',getheaders=True),mainkeys,delimit)
#########################################################################################################################
    def sortOccDic(self,occdic,mainkeys,delimit):   ### Converts Ranks and Positions to preZero strings
        '''Converts Ranks and Positions to preZero strings.'''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        maxval = {}     # Dictionary of maximum values
        for mkey in ['Rank','Start_Pos','End_Pos']: maxval[mkey] = 0
        for okey in occdic:
            if okey == 'Headers': continue
            for mkey in ['Rank','Start_Pos','End_Pos']:
                if mkey in occdic[okey]:
                    try: maxval[mkey] = max(maxval[mkey],rje.atoi(occdic[okey][mkey]))
                    except: self.errorLog('%s::%s' % (mkey,occdic[okey][mkey]))
        ### ~ [2] ~ Remake occdic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for okey in rje.sortKeys(occdic):
            if okey == 'Headers': continue
            odata = occdic.pop(okey)
            newkey = []
            for mkey in mainkeys:
                if mkey in maxval:
                    try: newkey.append(rje.preZero(rje.atoi(odata[mkey]),maxval[mkey]))
                    except: newkey.append(odata[mkey])
                else: newkey.append(odata[mkey])
            occdic[rje.join(newkey,delimit)] = odata
        return occdic
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
        if frontpage: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % rje.replace(stylesheet,'../','./'))
        else: html.append('<link rel="stylesheet" href="%s" TYPE="text/css" MEDIA="screen">' % stylesheet)
    if tabber:
        html += ['','<!-- ~~~~~~~~~~~ Tabber Javascript ~~~~~~~~~~~ -->','<script type="text/javascript">',
                 'document.write(\'<style type="text/css">.tabber{display:none;}<\/style>\');',
                 'var tabberOptions = {',' /* Optional: start manually (run tabber) at end of file','*/',
                 '\'manualStartup\':true','};','</script>','<!-- Load the tabber code -->']
        if frontpage: html += ['<script type="text/javascript" src="./javascript/tabber.js"></script>','']
        else: html += ['<script type="text/javascript" src="../javascript/tabber.js"></script>','']
    html += ['</head>','<!-- ~~~~~~~~~~~~~~~ End of HTML head data ~~~~~~~~~~~~~~~~~ -->','','<body>','']
    #print rje.join(html,'\n')
    return rje.join(html,'\n')
#########################################################################################################################
def htmlTail(copyright='RJ Edwards 2010',tabber=True):  ### Returns text for bottom of HTML
    '''
    Returns text for bottom of HTML.
    >> copyright:str = copyright text'
    >> tabber:bool = whether page has tabber tabs
    '''
    t = rje.split(time.asctime(time.localtime(time.time())))
    datetime = '%s %s %s' % (t[2],t[1],t[-1])
    html = ['','<!-- ~~~~~~~~~~~~~~ HTML tail data ~~~~~~~~~~~~~~~~~ -->',
            '<HR><FONT COLOR=#979E45 SIZE=2>&copy; %s. Last modified %s.</FONT></P>' % (copyright,datetime),'',
            '<script type="text/javascript">','/* manualStartup=true so need to run it now */',
            'tabberAutomatic(tabberOptions);','</script>','','</body>','</html>',
            '<!-- ~~~~~~~~~~~~~~ End of HTML tail data ~~~~~~~~~~~~~~~~~ -->']
    #print rje.join(html,'\n')
    return rje.join(html,'\n')
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
        if len(tab) > 2: html += rje.split(tabberTabHTML(tab[0],tab[1],tab[2]),'\n')
        else: html += rje.split(tabberTabHTML(tab[0],tab[1]),'\n')
    html += ['</div>','<!-- ~~~~~~~~~~~~~~~ End of %s Tabber Div ~~~~~~~~~~~~~~~ -->' % id,]
    #print rje.join(html,jointxt)
    return rje.join(html,jointxt)
#########################################################################################################################
def tabberTabHTML(id,text,title=''):          ### Returns text for TabberTab HTML
    '''
    Returns text for TabberTab HTML.
    >> title:str = Text for title of TabberTab
    >> text:str = HTML text for TabberTab content
    '''
    if not title: title = id
    html = ['','<!-- ~~~ %s TabberTab div ~~~ -->' % id,'<div class="tabbertab" title="%s" id="%s">' % (title,id),'']
    html += rje.split(text,'\n')
    html += ['','</div>','<!-- ~~~ %s TabberTab end ~~~ -->' % id]
    #print rje.join(html,'\n  ')
    if rje.join(html).upper().find('<PRE>') >= 0: return rje.join(html,'\n')
    else: return rje.join(html,'\n  ')
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
    try: SLiMHTML(mainlog,cmd_list).run()

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
