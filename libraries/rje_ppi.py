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
Module:       RJE_PPI
Description:  RJE Protein-Protein Interaction Module
Version:      2.9.0
Last Edit:    14/05/19
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is contains methods for manipulating protein-protein interaction (PPI) dictionaries. Database-specific
    classes and methods can be found in rje_hprd and rje_biogrid, while broader-scale functionality is found in the
    PINGU application. With time, generic PINGU functions should be migrated to this module and used by PINGU.
    
    The main purpose of rje_ppi is for use within other applications but there is also some standalone functionality
    for reading in a pairwise PPI file of delimited data with "Hub", "Spoke" and "Evidence" columns.
    
    Data is stored as an unannotated graph {Hub:[Spokes]} or an annotated graph {Hub:{Spoke:Evidence}} and rje_db 
    Database tables that store additional Node and Edge data. Edge data is read directly from PPI Pairwise input, with
    'Hub' and 'Edge' fields that combine to make the key. The Edge table will have an "Evidence" field and optional 
    additional fields. The Node table will have a 'Node' field as the key (that matches the 'Hub' or 'Spoke' from the
    Edge table). 

Commandline:
    ### Input Options: ###
    pairwise=FILE   : Input PPI pairwise file, containing Hub, Spoke and optional Evidence columns []
    nodelist=LIST   : Reduce input PPI to given Node list []
    nodemap=X       : Try to map Nodes first using Node table field X []
    expandppi=X     : Expand reduced Node list by X PPI levels [0]
    combineppi=LIST : List of pairwise PPI files (using same ID set) to import. Uses Evidence field or filename []
    nodefields=LIST : List of alternative A/B fields to replace Hub and Spoke fields [SYMBOL_,Gene_]
    ppisym=T/F      : Whether to enforce Hub/Spoke symmetry [True]

    ### Output/Processing Options: ###
    ppiout=FILE     : Save pairwise PPI file following processing (if rest=None) (T=basefile.pairwise.tdt) [None]
    tabout=T/F      : Output PPI data as Node and Edge tables [False]
    fragment=T/F    : Perform PPI fragmentation [False]
    fragsize=X      : Combine smaller fragments upto fragsize [200]
    minfrag=X       : Minimum fragment size to keep [3]
    mcode=T/F       : Perform MCODE clustering [False]
    xgmml=T/F       : Output network in XGMML style [False]
    xdir=PATH       : Directory for XGMML output [./XGMML]

    ### Layout Options: ###
    layout=X        : Layout to be used for XGMML output [spring]
    walltime=X      : Walltime (hours) for layouts [0.02]
    damping=X       : Force Directed Layout, damping parameter [0.9]
    colbydeg=T/F    : Whether to colour PNG output by node degree [False]
    nudgecyc=X      : Number of cycles between node nudges (try to bump out of unstable equilibria) [1000]

    ### MCODE Options: ###
    haircut=T/F     : Whether to perform "haircut" on MCODE complexes [False]
    multicut=T/F    : Whether to perform "haircut" on MCODE complexes for the purposes of looking at nodes 2+ times [True]
    fluff=X         : MCODE "fluff" threshold. <0 = No Fluff [0.5]
    vwp=X           : MCODE vertex weighting percentage [0.2]
    mink=X          : MCODE min k-core values [2]
    mindeg=X        : MCODE min degree for node scoring [2]

See also rje.py generic commandline options.

Uses general modules: copy, glob, math, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_xgmml, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, random, string, sys, time, math
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_svg, rje_xgmml, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 2.0 - Initial Compilation from earlier version. Complete reworking, including generation of PPI Object.
    # 2.1 - Added method for generating subnetworks based on +/- scoring.
    # 2.2 - Minor misc modifications.
    # 2.3 - Added extra Complex options.
    # 2.4 - Added SVG output.
    # 2.5 - Minor tweaks to R output and walltime default.
    # 2.6 - Added addPPI(hub,spoke,evidence) method. Added nodelist option.
    # 2.7 - Added tabout=T/F Output PPI data as Node and Edge tables [False]
    # 2.8 - Tweaked Spring Layout. Stores original Hub and Spoke Field.
    # 2.8.1 - Fixed bug with Spring Layout interruption message.
    # 2.9.0 - Added ppiout=FILE : Save pairwise PPI file following processing (if rest=None) [None]
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Create working PPI object.
    # [Y] : Add MCODE clustering algorithm.
    # [ ] : Add methods for saving Pairwise PPI tables.
    # [Y] : Add XGMML output.
    # [Y] : Add Spring-Embedded and/or Force-Directed layouts.
    # [ ] : Tidy so that all methods simply add to data and then add data clearing and purging methods.
    # [ ] : Improve discrimination (and conversion) between "flat" list and dictionary ppi dictionaries.
    # [ ] : Improve handling of unconnected parts of PPI graphs in rjeSpring Layout. (Decrease separation distance.)
    # |----> Fragment, layout each fragment and then reassemble next to each other? Or add gentle pull towards centre?
    # [ ] : Add iterative MCODE clustering (clusters of clusters)
    # [ ] : Add more clustering methods, such as Edge clustering.
    # [Y] : Adapt xgmml output to deal with non-symmetrical PPI.
    # [X] : Add option to output (evidence) instead of (pp) for XGMML etc. (Already uses optional Type field.)
    # [ ] : Rename walltime=X to be layout-specific. (Issues when running on supercomputer!)
    # [ ] : Upgrade to rje_obj object and add REST output.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('RJE_PPI', '2.9.0', 'May 2019', '2010')
    description = 'RJE Protein-Protein Interaction Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PPI Class                                                                                               #
#########################################################################################################################
class PPI(rje.RJE_Object):     
    '''
    Protein-Protein Interaction Class. Author: Rich Edwards (2010).

    Info:str
    - HubField = Original hub field heading read in from Pairwise PPI []
    - SpokeField = Original spoke field heading read in from Pairwise PPI []
    - Layout = Layout to be used for XGMML output [spring]
    - NodeMap = Try to map Nodes first using Node table field X []
    - Pairwise = Input PPI pairwise file, containing Hub, Spoke and optional Evidence columns []
    - PPIOut=FILE : Save pairwise PPI file following processing (if rest=None) [None]
    - XDir = Directory for XGMML output [./XGMML]
    
    Opt:boolean
    - ColByDeg = Whether to colour PNG output by node degree [False]
    - Fragment = Perform PPI fragmentation [False]
    - Haircut = Whether to perform "haircut" on MCODE complexes [False]
    - MCODE = Perform MCODE clustering [False]
    - MultiCut = Whether to perform "haircut" on MCODE complexes for the purposes of looking at nodes 2+ times [True]
    - PPISym = Whether to enforce Hub/Spoke symmetry [True]
    - TabOut = Output PPI data as Node and Edge tables [False]
    - XGMML = Output network in XGMML style [True]

    Stat:numeric
    - Damping = Force Directed Layout, damping parameter [0.9]
    - ExpandPPI = Expand reduced Node list by X PPI levels [0]
    - FragSize = Combine smaller fragments upto fragsize [200]
    - MinFrag = Minimum fragment size to keep [3]
    - Fluff = MCODE "fluff" threshold. <0 = No Fluff [0.5]
    - VWP = MCODE vertex weighting percentage [0.2]
    - MinK = MCODE min k-core values [2]
    - MinDeg = MCODE min degree for node scoring [2]
    - NudgeCyc = Number of cycles between node nudges (try to bump out of unstable equilibria) [1000]
    - Walltime = Walltime (hours) for layouts [1.0]
    
    List:list
    - Fragment = Generated if network fragmented
    - NodeList = Reduce input PPI to given Node list []
    - CombinePPI = List of pairwise PPI files (using same ID set) to import. Uses Evidence field or filename []
    - NodeFields = List of alternative A/B fields to replace Hub and Spoke fields [SYMBOL_,Gene_]

    Dict:dictionary
    - EdgeAtt = Dictionary of edge attributes {Att:Type}
    - NodeAtt = Dictionary of node attributes {Att:Type}
    - PPI = Main PPI dictionary: {hub:{spoke:evidence}}
    - Complex = Dictionary of {seed:[nodes]}

    Obj:RJE_Objects
    - DB = rje_db.Database object. Used to store Node and Edge data in 'Node' and 'Edge' tables.
    - SVG = rje_svg.SVG object. Used to create SVG output.
    '''
#########################################################################################################################
    def ppi(self,hub=None):
        if not hub: return self.dict['PPI']
        return self.dict['PPI'][hub]
    def hubs(self): return rje.sortKeys(self.ppi())
    def nodeNum(self,istr=False):
        if istr: return rje.integerString(len(self.dict['PPI']))
        else: return len(self.dict['PPI'])
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Layout','NodeMap','Pairwise','XDir','HubField','SpokeField','PPIOut']
        self.optlist = ['Fragment','MCODE','XGMML','Haircut','Fluff','ColByDeg','MultiCut','TabOut','PPISym']
        self.statlist = ['ExpandPPI','FragSize','MinFrag','Fluff','VWP','MinK','MinDeg','Walltime','Damping','NudgeCyc']
        self.listlist = ['CombinePPI','NodeFields','NodeList']
        self.dictlist = ['EdgeAtt','NodeAtt','PPI']
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Layout':'spring','XDir':rje.makePath('XGMML/')})
        self.setStat({'FragSize':200,'MinFrag':3,'Fluff':0.5,'VWP':0.2,'MinK':2,'MinDeg':2,'Walltime':0.02,'Damping':0.9,
                      'NudgeCyc':1000,'ExpandPPI':0})
        self.setOpt({'Fluff':True,'ColByDeg':False,'MultiCut':True,'PPISym':True})
        self.list['NodeFields'] = ['SYMBOL_','Gene_']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['SVG'] = None
        self.obj['DB'] = rje_db.Database(self.log, self.cmd_list)
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
                self._cmdReadList(cmd,'file',['Pairwise','PPIOut'])
                self._cmdReadList(cmd,'path',['XGMML'])
                self._cmdReadList(cmd,'info',['Name','NodeMap','Layout'])
                self._cmdReadList(cmd,'int',['ExpandPPI','FragSize','MinFrag','MinK','MinDeg','NudgeCyc'])
                self._cmdReadList(cmd,'stat',['Fluff','VWP','Walltime','Damping'])
                self._cmdReadList(cmd,'opt',['Fragment','MCODE','PPISym','XGMML','ColByDeg','MultiCut','TabOut'])
                self._cmdReadList(cmd,'list',['NodeFields','NodeList'])
                self._cmdReadList(cmd,'glist',['CombinePPI'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.opt['Fluff'] = self.stat['Fluff'] >= 0.0
        if self.info['PPIOut'].lower() in ['t','true']:
            ptxt = 'ppiout=%s' % self.info['PPIOut'].upper()
            self.info['PPIOut'] = '%s.pairwise.tdt' % self.baseFile()
            self.printLog('#PPIOUT','%s: set ppiout=%s' % (ptxt,self.info['PPIOut']))
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.main()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def test(self): ### Test code
        xfile = '%s.xgmml' % self.info['Basefile']
        if os.path.exists(xfile):
            print('Code up reading of XGMML!')
        G = randomGraph(200,400,name='Gene')
        npos = self.rjeSpringLayout(G)
        self.opt['ColByDeg'] = True
        svghtm = self.saveSVG(npos,basefile=self.info['Basefile'],G=G,font=0,width=1600,ntype='ellipse',backups=True)
        self.saveR(npos,basefile=self.info['Basefile'],G=G,cleantdt=False,backups=True)
        xgmml = self.ppiToXGMML(G,self.info['Basefile'])
        xgmml.saveXGMML(xfile)
        import rje_html        
        html = rje_html.htmlHead('PPI SVG Test',stylesheets=[],tabber=False) + svghtm
        html += '\n<hr>\n<img src=%s.png width=1600>\n' % self.info['Basefile']
        html += rje_html.htmlTail(copyright='RJ Edwards 2010',tabber=False)
        open('%s.htm' % self.info['Basefile'],'w').write(html)
#########################################################################################################################
    def main(self):  ### Main run method following loading etc.
        '''Main run method following loading etc..'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xbase = os.path.split(self.info['Basefile'])
            xbase = rje.makePath('%s%s%s' % (rje.makePath(xbase[0]),self.info['XDir'],xbase[1]),wholepath=True)
            if self.opt['XGMML']: rje.mkDir(self,xbase)
            if self.opt['Fragment']:
                self.fragment()
                #!# Add output in Grid layout using Frag Size etc. #!#
                if self.opt['MCODE']:
                    for partial in self.list['Fragment']: self.MCODE(subGraph(self.ppi(),partial),mink=self.stat['MinK'],vwp=self.stat['VWP'],haircut=self.opt['Haircut'],fluff=self.opt['Fluff'],fluffd=self.stat['Fluff'],mindeg=self.stat['Fluff'])
                elif self.opt['XGMML']:
                    if self.info['Layout'] == 'spring':
                        i = 1
                        for frag in self.list['Fragment']: # Add orphans option
                            name = '%s.%s' % (self.info['Basefile'],rje.preZero(i,len(self.list['Fragment']))); i+=1
                            self.springXGMML('%s.%s.xgmml' % (xbase,rje.preZero(i,len(self.list['Fragment']))),G=subGraph(self.ppiG(),frag),name=name)
            elif self.opt['MCODE']: self.MCODE(mink=self.stat['MinK'],vwp=self.stat['VWP'],haircut=self.opt['Haircut'],fluff=self.opt['Fluff'],fluffd=self.stat['Fluff'],mindeg=self.stat['Fluff'])
            elif self.opt['XGMML']:
                if self.info['Layout'] == 'spring':
                    self.springXGMML('%s.xgmml' % (self.info['Basefile']),G=self.ppiG(),name=self.info['Basefile'])
            #if self.opt['MCODE']: self.debug(self.dict['Complex'])
            if self.opt['MCODE'] and self.opt['XGMML']:
                if self.info['Layout'] == 'spring':
                    for seed in rje.sortKeys(self.dict['Complex']):
                        rank = self.db('Complex').data()[seed]['Rank']; N = self.db('Complex').data()[seed]['Nodes']
                        xfile = '%s.%s-N%d.xgmml' % (xbase,rje.preZero(rank,len(self.dict['Complex'])),N)
                        self.springXGMML(xfile,G=subGraph(self.ppi(),self.dict['Complex'][seed]),name=seed)
            ### ~ [2] Output of tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('PPIOut'): self.db('Edge').saveToFile(filename=self.getStr('PPIOut'))
            if self.getBool('TabOut'):
                self.db('Node').saveToFile()
                self.db('Edge').saveToFile(savefields=self.list['EdgeFields'])
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            load = self.opt['Fragment'] or self.opt['MCODE'] or self.opt['XGMML']  or self.getBool('TabOut') or self.getStrLC('PPIOut')
            self.list['EdgeFields'] = ['Hub','Spoke']
            if rje.exists(self.info['Pairwise']):
                if self.info['Basefile'].lower() in ['','none']: self.baseFile(rje.baseFile(self.info['Pairwise']))
                if self.info['Pairwise'] in self.list['CombinePPI']: self.list['CombinePPI'].remove(self.info['Pairwise'])
                self.list['CombinePPI'].insert(0,self.info['Pairwise'])
            if load and self.list['CombinePPI']:
                for ppifile in self.list['CombinePPI']: self.loadPairwisePPI(ppifile,clear=False,asdict=True,sym=self.getBool('PPISym'))
            ### ~ [2] NodeList Reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                    
                if self.list['NodeList']:
                    #self.deBug(self.list['NodeList'])
                    if self.info['NodeMap'] not in ['','none']:
                        nx = len(self.list['NodeList'])
                        for node in self.list['NodeList'][0:]:
                            #self.deBug(node)
                            #try: self.deBug(self.db('Node').data()[node])
                            #except: pass
                            if node in self.db('Node').index(self.info['NodeMap']):
                                #self.deBug(self.db('Node').indexDataList(self.info['NodeMap'],node,'Node'))
                                for mapped in self.db('Node').indexDataList(self.info['NodeMap'],node,'Node'):
                                    if mapped not in self.list['NodeList']: self.list['NodeList'].append(mapped)
                            else: self.printLog('#FAIL','Could not map %s "%s"' % (self.info['NodeMap'],node))
                        self.list['NodeList'].sort()
                        for node in self.list['NodeList'][0:]:
                            if node not in self.db('Node').data(): self.list['NodeList'].remove(node)
                        self.printLog('#MAP','%s nodes -> %s nodes after %s mapping' % (rje.iStr(nx),rje.iLen(self.list['NodeList']),self.info['NodeMap']))
                    nodelist = self.list['NodeList'][0:]
                    for i in range(self.getInt('ExpandPPI')):
                        for hub in nodelist[0:]:
                            if hub not in self.ppi(): continue
                            for spoke in self.ppi(hub):
                                if spoke not in nodelist: nodelist.append(spoke)
                        self.printLog('#ADD','Expand PPI NodeList (Level %d): %s -> %s nodes.' % (i+1,rje.iLen(self.list['NodeList']),rje.iLen(nodelist)))
                    #self.deBug(nodelist)
                    self.dict['PPI'] = subGraph(self.ppi(),nodelist,loops=True,addv=True)  # Addv as option?
                    if self.db('Node'): self.db('Node').dropEntriesDirect('Node',nodelist,inverse=True,log=True)
                    if self.db('Edge'):
                        self.db('Edge').dropEntriesDirect('Hub',nodelist,inverse=True,log=True)
                        self.db('Edge').dropEntriesDirect('Spoke',nodelist,inverse=True,log=True)                    
            if self.info['Name'].lower() in ['','none']: self.info['Name'] = self.info['Basefile']
            if self.info['Name'].lower() in ['','none']: self.info['Name'] = 'RJE_PPI'
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def splitComplex(self,filename=None,mcode=True):   ### Splits PPI network then predicts complexes and outputs
        '''Splits PPI network then predicts complexes and outputs.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename: self.loadPairwisePPI(filename,asdict=False)
            ### ~ [2] Split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.fragment()
            ### ~ [3] Complexes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for partial in self.list['Fragment']: self.MCODE(subGraph(self.ppi(),partial))
            ### ~ [4] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seed in rje.sortKeys(self.dict['Complex']):
                self.springXGMML('%s.%s.xgmml' % (self.info['Basefile'],seed),G=subGraph(self.ppi(),self.dict['Complex'][seed]),name=seed)
        except: self.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
    def ppiFromEdges(self,sym=None,clear=True): ### Generates self.dict['PPI'] from Edge table data
        '''
        Generates self.dict['PPI'] from Edge table data.
        >> sym:bool [None] = Whether to enforce symmetry. (Use self.getBool('PPISym') if None)
        >> clear:bool [True] = Whether to clear existing PPI.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear: G = self.dict['PPI'] = {}
            else: G = self.dict['PPI']
            if sym == None: sym = self.getBool('PPISym')
            for hub in self.db('Node').dataKeys():
                if hub not in G: G[hub] = {}
            ### ~ [1] Populate dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in self.db('Edge').entries():
                hub = entry['Hub']; spoke = entry['Spoke']
                if hub not in G: G[hub] = {}; self.warnLog('Hub "%s" missing from Node table?' % (hub),'missing_hubnode',suppress=True)
                if spoke in G[hub]: self.warnLog('Spoke "%s" over-written for Hub "%s"' % (spoke,hub),'spoke_overwrite',suppress=True)
                G[hub][spoke] = entry['Evidence']
            if not sym: return G
            for entry in self.db('Edge').entries():
                spoke = entry['Hub']; hub = entry['Spoke']
                if hub not in G: G[hub] = {}; self.warnLog('Spoke "%s" missing from Node table?' % (hub),'missing_hubnode',suppress=True)
                if spoke in G[hub]: continue
                G[hub][spoke] = entry['Evidence']
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <3> ### General PPI Methods                                                                                     #
#########################################################################################################################
    def checkIntegrity(self,mode='full',sym=True,orphans=True): ### Method to check certain conditions are met
        '''
        Method to check certain conditions are met and data has full integrity
        >> mode:str ['full'] = Assessment criteria
        # 'full' = Add nodes and edges from tables into self.ppi() and visa versa
        # 'ppi' = Add ppi nodes/edges to tables but remove missing ppi nodes/edges from tables
        # 'db' = Add nodes/edges to ppi dict but remove nodes/edges not in tables
        # 'min' = Only keep nodes/edges in both ppi dict and tables
        # 'comp' = Make sure all nodes in complexes are in ppi dict and tables, removing any that are not
        >> sym:bool [True] = Enforce symmetry before integrity check
        >> orphans:bool [True] = Allow orphan (unconnected) nodes
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sym: self.symmetry()
            if not orphans: self.purgeOrphans()
            ### ~ [2] Enforce Integrity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add rest of checking code here. #!#
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def fragment(self):     ### Fragments the PPI network into List of unconnected Nodes
        '''Fragments the PPI network into List of unconnected Nodes.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Fragment'] = fragment(self.ppi(),minfrag=self.stat['MinFrag'],callobj=self)
            #!# Output table of nodes that are not in a Fragment #!#
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def purgeOrphans(self): ### Removes all nodes without edges
        '''Removes all nodes without edges.'''
        try:### ~ [1] Purge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nx = self.nodeNum(); self.progLog('#PURGE','Purging unconnected nodes...')
            for node in self.hubs():
                if not self.ppi()[node]: self.ppi().pop(node)
            self.printLog('\r#PURGE','Purged unconnected nodes: %s -> %s nodes' % (rje.integerString(nx),self.nodeNum(True)))
            ndb = self.db('Node')
            if ndb: ndb.dropEntriesDirect('Node',self.hubs(),inverse=True,log=True)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def ppiG(self,loops=False): # Generate PPI lists, removing any self interactors (if loops=False)
        '''Generate PPI Graph (lists), removing any self interactors (if loops=False).'''
        ppi = {}
        for hub in self.hubs():
            try: ppi[hub] = rje.sortKeys(self.ppi()[hub])
            except: ppi[hub] = self.ppi()[hub][0:]; ppi[hub].sort()
            if hub in ppi[hub] and not loops: ppi[hub].remove(hub)
            if not ppi[hub]: ppi.pop(hub)
        return ppi
#########################################################################################################################
    def genePPI(self,gene):     ### Returns list of interactors with given gene
        '''Returns flattened graph of interactors with given gene.'''
        if gene not in self.ppi(): return []
        return rje.sortKeys(self.ppi()[gene])
#########################################################################################################################
    def symmetry(self,ppi=None,inplace=True):    ### Enforces symmetry
        '''
        Enforces symmetry.
        >> ppi:dict [None] = PPI graph dictionary. If None, will use self.ppi()
        >> inplace:bool [True] = Whether to update PPI graph in place (True) or return new dictionary (False)
        '''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not ppi: ppi = self.ppi()
        if not ppi: return ppi
        edb = self.db('Edge')
        asdict = True
        try: rje.sortKeys(ppi[ppi.keys()[0]])
        except: asdict = False
        if not inplace:
            newppi = {}
            for hub in rje.sortKeys(ppi):
                if asdict: 
                    newppi[hub] = {}
                    for spoke in ppi[hub]: newppi[hub][spoke] = ppi[hub][spoke]
                else: newppi[hub] = ppi[hub][0:]
            ppi = newppi
        ### ~ [2] ~ Symmetry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        sx = 0; hx = len(ppi)
        self.progLog('\r#SYM','Imposing symmetry on %s of %s input hubs -> %s hubs.' % (rje.iStr(sx),rje.iStr(hx),rje.iLen(ppi)))
        for spoke in rje.sortKeys(ppi):
            sx += 1
            for hub in ppi[spoke]:
                if hub not in ppi:
                    if asdict: ppi[hub] = {}
                    else: ppi[hub] = []
                if spoke not in ppi[hub]:
                    self.progLog('\r#SYM','Imposing symmetry on %s of %s input hubs -> %s hubs.' % (rje.iStr(sx),rje.iStr(hx),rje.iLen(ppi)))
                    if asdict: ppi[hub][spoke] = self.evidence(hub,spoke)
                    else: ppi[hub].append(spoke)
                    if edb:
                        sentry = edb.data(edb.makeKey({'Hub':spoke,'Spoke':hub}))
                        hentry = {'Hub':hub,'Spoke':spoke}
                        if edb.data(edb.makeKey(hentry)): hentry = edb.data(edb.makeKey(hentry))
                        for field in sentry:
                            if 'Hub' in field: symfield = rje.replace(field,'Hub','Spoke')
                            else: symfield = rje.replace(field,'Spoke','Hub')
                            if symfield not in hentry or not hentry[symfield]: hentry[symfield] = sentry[field]
                        if not edb.data(edb.makeKey(hentry)): edb.addEntry(hentry)
        self.printLog('\r#SYM','Imposed PPI symmetry on %s of %s input hubs -> %s hubs.' % (rje.iStr(sx),rje.iStr(hx),rje.iLen(ppi)))
                                          
        ### ~ [2] ~ Tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not asdict:
            for hub in ppi: ppi[hub].sort()
        return ppi
#########################################################################################################################
    def complexOnly(self,G,comtext='complex'):    ### Reduces interactions to those with "complex" support
        '''Reduces interactions to those with "complex" support.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not G: return G
            V = rje.sortKeys(G)
            asdict = True
            try: G[V[0]].keys()
            except: asdict = False
            for vi in V:
                if asdict: VJ = rje.sortKeys(G[vi])
                else: VJ = G[vi][0:]
                for vj in VJ:
                    if self.evidence(vi,vj).lower().find(comtext.lower()) >= 0: continue
                    if asdict: G[vi].pop(vj)
                    else: G[vi].remove(vj)
                    if vi in G[vj]:
                        if asdict: G[vj].pop(vi)
                        else: G[vj].remove(vi)
            return G    # NB. Edited in place anyway
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def expandPPI(self,G,ppexpand=1,flatten=False,ppcomplex=False,comtext='complex'):   ### Expands PPI Graph
        '''Expands PPI Graph.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if flatten: G = subGraph(G,rje.sortKeys(G))
            next = xnodes = rje.sortKeys(G)
            for ex in range(ppexpand):
                prev = next[0:]
                next = []
                gx = 0.0; gtot = len(prev)
                for g1 in prev[0:]:
                    self.progLog('\r#PPI','Expanding PPI network (Lvl %d): %.2f%%' % (ex+1,gx/gtot)); gx += 100.0
                    if g1 not in self.ppi(): continue
                    for g2 in self.ppi()[g1]:
                        if ppcomplex and self.evidence(g1,g2).lower().find(comtext.lower()) < 0: continue
                        if g2 not in xnodes:     # Add this gene
                            next.append(g2)
                            xnodes.append(g2)
                            if g2 not in G: G[g2] = []
                            if g1 not in G[g2]: G[g2].append(g1)
                            if g2 not in G[g1]: G[g1].append(g2)
                self.printLog('\r#PPI','Expanded PPI network (Lvl %d): %s nodes added' % (ex+1,len(next)))
            return G
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### Clustering Methods                                                                                      #
#########################################################################################################################
    def cmdMCODE(self,G,save='Complex'): return self.MCODE(G,mink=self.stat['MinK'],vwp=self.stat['VWP'],haircut=self.opt['Haircut'],fluff=self.opt['Fluff'],fluffd=self.stat['Fluff'],mindeg=self.stat['Fluff'],save=save)
    def MCODE(self,ppi=None,mink=2,vwp=0.2,haircut=True,fluff=False,fluffd=0.5,mindeg=2,save='Complex'):    ### Generates self.list['Clusters'] from self.ppi() using MCODE algorithm.
        '''Generates self.list['Clusters'] from self.ppi() using MCODE algorithm.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#MCODE','MCODE clustering: mink=%s, vwp=%s, haircut=%s, fluff=%s, fluffd=%s, mindeg=%s, multicut=%s' % (mink,vwp,haircut,fluff,fluffd,mindeg,self.opt['MultiCut']))
            #mink = 2        # Minimum k-core value.
            #mindeg = 2      # Minimum degree for node scoring.
            #haircut = True  # Perform haircut removal of singletons.
            #fluff = False   # Perform fluff addition analysis.
            #fluffd = 0.5    # Fluff density threshold
            #maxdepth = 100  # Maximum distance from seed to search for clusters.
            #vwp = 0.2   # Vertex weight percentage
            ## ~ [0a] ~ Generate PPI lits, removing any self interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fluff = fluff and fluffd >= 0.0
            if not ppi: ppi = self.ppiG()
            ### ~ [1] ~ Stage 1: Vertex Weighting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vweight = {}    # Dictionary of node weights
            vx = 0.0; vtot = len(ppi)
            for v in rje.sortKeys(ppi):
                self.progLog('\r#VWT','MCODE vertex weighting: %.2f%%' % (vx/vtot)); vx += 100.0
                N = [v] + ppi[v]              # NB. v itself cannot be included here, or the density of K would have to be 1!
                if len(N) < mindeg: vweight[v] = 0.0; continue
                GN = subGraph(ppi,N)    # Graph of N
                (K,k) = kCoreMax(GN,mink)  # K = highest k-core graph from N; k = highest k-core number from N
                GK = subGraph(ppi,K)   # Graph of K + v
                d = density(GK)         # Get density of K + v
                vweight[v] = k * d      # MCODE vertex weighting
                #print v, N, GK, K, k, d, vweight[v]
            self.printLog('\r#VWT','MCODE vertex weighting for %s nodes complete.' % rje.integerString(vtot))
            ### ~ [2] ~ Stage 2: Molecular Complex Prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seen = []   # List of seen vertexes (no multiple membership at this stage)
            vx = 0.0; complexes = []
            for s in rje.dictKeysSortedByValues(vweight,revsort=True):  # Take each seed vertex in turn
                self.progLog('\r#COM','MCODE finding clusters: %.2f%%' % (vx/vtot)); vx += 100.0
                if s in seen: continue
                tried = []          # List of vertexes assessed for this complex
                complex = [s]       # Current complex being predicted
                seen.append(s)
                growth = ppi[s][0:] # List of vertexes to check for inclusion
                while growth:
                    v = growth.pop(0)
                    if v in tried + seen: continue
                    if vweight[v] >= vweight[s]*(1 - vwp):   # NB. In paper this is just > but I am not sure that would work
                        complex.append(v); seen.append(v)
                        for n in ppi[v]:
                            if n not in complex+growth+tried: growth.append(n)
                    else: tried.append(v)
                if self.opt['MultiCut']:
                    G = subGraph(ppi,complex)
                    kcomp = kCore(G,2)
                    for loner in rje.listDifference(kcomp,complex): seen.remove(loner)
                complexes.append(complex)
            self.printLog('\r#COM','MCODE found %s clusters for %s nodes.' % (rje.integerString(len(complexes)),rje.integerString(vtot)))
            ### ~ [3] ~ Stage 3: Post-Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.db(save); cdbackup = not cdb
            if not cdb:
                cdb = self.db().addEmptyTable(save,['Seed','Nodes','Edges','Density','Score','Complex'],keys=['Seed'])
                self.dict[save] = {}
            cx = 0.0; ctot = len(complexes)
            ndensity = {}   # Dictionary of neighbourhood density scores (if fluff used)
            for complex in complexes[0:]:
                self.progLog('\r#POST','Post-processing MCODE clusters: %.f%%' % (cx/ctot)); cx += 100.0
                seed = complex[0]
                ## ~ [3a] ~ Remove complexes without 2-core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                G = subGraph(ppi,complex)
                #self.deBug(G); self.deBug(kCore(G,2))
                if not kCore(G,2): continue
                ## ~ [3b] ~ Optional addition of fluff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if fluff:
                    fluffed = complex[0:]
                    for v in complex[0:]:
                        for n in ppi[v]:
                            if n in fluffed: continue
                            if n not in ndensity: ndensity[n] = density(subGraph(ppi,ppi[n]+[n]))
                            if ndensity[n] >= fluffd: complex.append(n)
                            fluffed.append(n)
                    if len(complex) > len(G): G = subGraph(ppi,complex)
                # ~ [3c] ~ Optional haircut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if haircut:
                    complex = kCore(G,2)
                    G = subGraph(G,complex)
                complex.sort()
                E = 0; N = len(complex)
                for v in G: E += len(G[v])
                if seed in self.dict[save]:
                    self.printLog('\r#WARN','Cluster "%s" being overwritten! |%s| >> |%s|' % (seed,rje.join(self.dict[save][seed],'|'),rje.join(complex,'|')))
                self.dict[save][seed] = complex
                cdb.data()[seed] = {'Seed':seed,'Complex':rje.join(complex,'|'),'Nodes':N,'Edges':E/2,'Density':E/float(N*(N-1))}
                cdb.data()[seed]['Score'] = cdb.data()[seed]['Nodes'] * cdb.data()[seed]['Density'] 
            self.printLog('\r#POST','Post-processing MCODE: %s total clusters.' % rje.integerString(cdb.entryNum()))
            cdb.rankField('Score',newfield='Rank',rev=True,absolute=True,lowest=True,unique=True)
            if save.lower() not in ['','none']: cdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],save.lower()),backup=cdbackup)
        except: self.errorLog('Problem with MCODE implementation')
#########################################################################################################################       
    def valueClusters(self,ppi=None,vfield='Score',addlinks=False,save='valnet'):   ### Generate subnetworks based on node values for vfield
        '''
        Generate subnetworks based on node values for vfield.
        >> ppi:dict = PPI dictionary. Uses self.ppiG() if None.
        >> vfield:str ['Score'] = Field in 'Node' table used for scores. Also used for DB table name.
        >> addlinks:bool [False] = Whether to add links between pairs of subnetwork nodes.
        >> save:str ['valnet'] = Name used for saving to table.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ndb = self.db('Node')
            if ppi == None: ppi = self.ppiG()
            if not ppi: return self.printLog('#VAL','No PPI data for value clustering')
            if vfield not in ndb.fields(): return self.printLog('#VAL','Value field "%s" not found for value clustering' % vfield)
            V = rje.sortKeys(ppi)
            ## ~ [1a] ~ Prepare data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scores = {}         # Dictionary of node:score
            vnodes = {0.0:[]}   # Dictionary of score:[nodes]
            for vi in V:
                #scores[vi] = 0.0; vx = 0
                #for vj in ppi[vi]:
                #    if vi == vj: continue
                #    scores[vi] += self.nodeData(vj,vfield,0.0); vx += 1
                #if vx: scores[vi] = scores[vi] / vx
                #scores[vi] = (scores[vi] + self.nodeData(vi,vfield,0.0)) / 2.0
                #if scores[vi] > 0 and self.nodeData(vi,vfield,0.0) < 0: scores[vi] = 0.0    # Both scores and mean scores must be + or -
                #elif scores[vi] < 0 and self.nodeData(vi,vfield,0.0) > 0: scores[vi] = 0.0
                scores[vi] = self.nodeData(vi,vfield,0.0)
                if scores[vi] not in vnodes: vnodes[scores[vi]] = []
                vnodes[scores[vi]].append(vi)
            self.printLog('#VAL','%s of %s nodes with score 0.0' % (rje.iStr(len(vnodes[0.0])),rje.iStr(len(V))))
            for v in vnodes.pop(0.0): V.remove(v)
            ### ~ [2] Generate networks from mean vfield values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            subnetworks = []    # List of subnetwork lists of nodes
            while V:
                if -rje.sortKeys(vnodes)[0] > rje.sortKeys(vnodes)[-1]: seedval = rje.sortKeys(vnodes)[0]
                else: seedval = rje.sortKeys(vnodes)[-1]
                seed = vnodes[seedval][0]
                subnet = []             # Nodes in subnetwork
                growth = [seed]         # List of nodes added in previous iteration
                checked = []
                while growth:
                    vi = growth.pop(0); subnet.append(vi)
                    for vj in ppi[vi]:
                        if vj in subnet + growth + checked: continue
                        if seedval < 0 and scores[vj] < 0: growth.append(vj)
                        elif seedval > 0 and scores[vj] > 0: growth.append(vj)
                        else: checked.append(vj)
                    #self.deBug('%s: %s' % (vi,growth))
                for v in subnet:
                    V.remove(v)
                    vnodes[scores[v]].remove(v)
                    if not vnodes[scores[v]]: vnodes.pop(scores[v])
                subnetworks.append(subnet[0:])
                #if len(subnet) > 2: self.deBug(subnet)
            ### ~ [3] ~ Post-Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.db(vfield); cdbackup = not cdb
            if not cdb:
                cdb = self.db().addEmptyTable('%snet' % vfield,['Seed','Nodes','Edges','Density','Score',vfield,'Members'],keys=['Seed'])
                self.dict['%snet' % vfield] = {}
            cx = 0.0; ctot = len(subnetworks)
            ndensity = {}   # Dictionary of neighbourhood density scores (if fluff used)
            for complex in subnetworks[0:]:
                self.progLog('\r#POST','Post-processing %s subnetworks: %.f%%' % (vfield,cx/ctot)); cx += 100.0
                ## ~ [3a] ~ Remove complexes without 3 members ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if len(complex) < 3: continue
                seed = complex[0]
                ## ~ [3b] ~ Optional addition of linking nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if addlinks:
                    core = complex[0:]; core.sort()
                    for vi in core:
                        for vj in ppi[vi]:
                            if vj in complex: continue
                            for vk in ppi[vj]:
                                if vj in complex: break
                                if vk in [vi,vj]: continue
                                if vk in core: complex.append(vj)
                complex.sort() 
                #if addlinks: self.deBug(core); self.deBug(complex)
                ## ~ [3c] ~ Scoring and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                G = subGraph(ppi,complex)
                meanval = 0.0
                for v in complex: meanval += (self.nodeData(v,vfield,0.0) / float(len(complex)))
                E = 0; N = len(G)
                for v in G: E += len(G[v])
                if seed in self.dict['%snet' % vfield]:
                    self.printLog('\r#WARN','Subnetwork "%s" being overwritten! |%s| >> |%s|' % (seed,rje.join(self.dict[vfield][seed],'|'),rje.join(complex,'|')))
                self.dict['%snet' % vfield][seed] = complex
                cdb.data()[seed] = {'Seed':seed,'Members':rje.join(complex,'|'),'Nodes':N,'Edges':E/2,
                                    'Density':E/float(N*(N-1)),vfield:meanval}
                cdb.data()[seed]['Score'] = cdb.data()[seed]['Nodes'] * cdb.data()[seed]['Density'] 
            self.printLog('\r#POST','Post-processing %s: %s total clusters.' % (vfield,rje.integerString(cdb.entryNum())))
            cdb.rankField('Score',newfield='Rank',rev=True,absolute=True,lowest=True,unique=True)
            if save.lower() not in ['','none']: cdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],save),backup=cdbackup)
        except: self.errorLog('Problem with ppi.valueClusters()')            
#########################################################################################################################       
    ### <5> ### I/O Methods                                                                                             #
#########################################################################################################################
    def loadPairwisePPI(self,filename,clear=True,asdict=None,sym=True,evidence={},resave=None):   ### Loads PPI from pairwise delimited file
        '''
        Loads PPI from pairwise delimited file.
        >> filename:str = Filename of delimited Pairwise PPI file (Hub, Spoke, [Evidence])
        >> clear:bool [True] = whether to clear the Edge table prior to loading
        >> asdict:bool [None] = whether to load PPI data as a dictionary of dictionaries (with evidence)
        >> sym:bool [True] = whether to load PPI data symmetrically
        >> evidence:dict = dictionary of {evidence type:filter type} (including partial matches) to:
        - 'ignore' = expunge from PPI. Pairs without other support will be deleted.
        - 'exclude' = pairs with listed evidence codes will be skipped, regardless of other evidence
        - 'only' = pairs without listed evidence codes will be skipped
        - '*_exact' = only exact matches to evidence type will count.
        >> resave:str = Filename for PPI database to be re-saved following evidence filtering.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Extract file header data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = rje.readDelimit(rje.chomp(open(filename,'r').readline()),rje.delimitFromExt(filename=filename))
            #self.deBug(headers)
            hubfield = None; spokefield = None
            if 'Hub' in headers: hubfield = 'Hub'
            else:
                for field in self.list['NodeFields']:
                    nfield = '%sA' % field
                    if nfield in headers: hubfield = nfield; break
            if 'Spoke' in headers: spokefield = 'Spoke'
            else:
                for field in self.list['NodeFields']:
                    nfield = '%sB' % field
                    if nfield in headers: spokefield = nfield; break
            if not (hubfield and spokefield): self.printLog('#ERR','Cannot find Hub and/or Spoke fields for %s' % filename); raise ValueError
            self.setStr({'HubField':hubfield,'SpokeField':spokefield})
            ## ~ [0b] ~ Load PPI data from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb = self.db('Edge')
            if clear or not edb: pdb = self.db().addTable(filename,mainkeys=[hubfield,spokefield],name='Edge')
            else: pdb = self.db().addTable(filename,mainkeys=[hubfield,spokefield],name='PPITemp')
            if hubfield != 'Hub': pdb.renameField(hubfield,'Hub')
            if spokefield != 'Spoke': pdb.renameField(spokefield,'Spoke')
            pdb.list['Fields'].remove('Hub')
            pdb.list['Fields'].remove('Spoke')
            pdb.list['Fields'] = ['Hub','Spoke'] + pdb.list['Fields']
            if asdict == None: asdict = 'Evidence' in pdb.fields()
            if clear: self.dict['PPI'] = {}
            ppi = self.dict['PPI']
            edgedata = []
            nodedata = []
            ## ~ [0a] ~ Evidence filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if evidence and 'Evidence' in pdb.fields():
                (ex,enum) = (0.0,pdb.entryNum())
                for entry in pdb.entries():
                    self.progLog('\r#PPI','Filtering PPI according to evidence: %.2f%%' % (ex/enum)); ex += 100.0
                    ekeep = True; okeep = None
                    ekey = pdb.makeKey(entry)
                    evlist = rje.split(entry['Evidence'],'|')
                    for ev in evlist[0:]:
                        (spoke,link) = rje.split(ev,':')
                        for etype in rje.sortKeys(evidence):
                            ftype = evidence[etype]
                            if (ftype[-5:] == 'exact' and link != etype) or (etype not in link):
                                if ftype.find('only') == 0: okeep = not not okeep
                                continue
                            if ftype.find('exclude') == 0: ekeep = False; break
                            if ftype.find('ignore') == 0: evlist.remove(ev)
                            if ftype.find('only') == 0: okeep = True
                        if not ekeep: break
                    if okeep == False: ekeep = False
                    if not ekeep or not evlist: pdb.data().pop(ekey)
                    else: entry['Evidence'] = rje.join(evlist ,'|')
                self.printLog('\r#PPI','Filtered PPI according to evidence: %s of %s PPI remain' % (rje.iStr(pdb.entryNum()),rje.iStr(enum)))
                if resave: pdb.saveToFile(resave,backup=self.opt['Backups'])
            ### ~ [1] ~ Update PPI dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppix = 0
            for entry in pdb.entries():
                (hub,spoke) = (entry['Hub'],entry['Spoke'])
                if hub not in ppi:
                    if asdict: ppi[hub] = {}
                    else: ppi[hub] = []
                if spoke in ppi[hub]: continue
                if asdict and 'Evidence' in pdb.fields(): ppi[hub][spoke] = entry['Evidence']
                elif asdict: ppi[hub][spoke] = rje.baseFile(filename,strip_path=True)
                else: ppi[hub].append(spoke)
                if edb and pdb != edb:
                    ekey = edb.makeKey(entry); edb.data()[ekey] = {}
                    for field in edb.fields():
                        if field in pdb.fields(): edb.data()[ekey][field] = entry[field]
                        else: edb.data()[ekey][field] = ''
                ppix += 1
            asdictxt = {True:'(as dictionaries)',False:'(as lists)'}[asdict]
            self.printLog('#PPI','%s PPI added from %s %s' % (rje.iStr(ppix),filename,asdictxt))
            ### ~ [2] ~ Tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not asdict:
                for hub in ppi: ppi[hub].sort()
            if sym: self.symmetry()
            ## ~ [2a] ~ Update tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ndb = self.db('Node')
            if clear or not ndb: ndb = self.db().addEmptyTable('Node',['Node'],['Node'])
            ndict = ndb.data()
            for field in pdb.fields():
                if field in ['Hub','Spoke']: edgedata.append(field)
                elif field.find('Hub') == 0:
                    nodedata.append(field)
                    if field[3:] not in ndb.fields(): ndb.addField(field[3:])
                elif field.find('Spoke') == 0:
                    nodedata.append(field)
                    if field[5:] not in ndb.fields(): ndb.addField(field[5:])
                else: edgedata.append(field)
                if field in edgedata and field not in self.list['EdgeFields']: self.list['EdgeFields'].append(field)
                if self.getStrLC('PPIOut') and field not in edgedata: edgedata.append(field)
            for entry in pdb.entries():
                #self.deBug(entry)
                if entry['Hub'] not in ndict: ndict[entry['Hub']] = {'Node':entry['Hub']}
                if entry['Spoke'] not in ndict: ndict[entry['Spoke']] = {'Node':entry['Spoke']}
                #self.deBug(ndict[entry['Hub']]); self.deBug(ndict[entry['Spoke']])
                for field in nodedata:
                    if field[:3] == 'Hub': ndict[entry['Hub']][field[3:]] = entry[field]
                    if field[:5] == 'Spoke': ndict[entry['Spoke']][field[5:]] = entry[field]
                #self.deBug(ndict[entry['Hub']]); self.deBug(ndict[entry['Spoke']])
            if pdb.info['Name'] == 'PPITemp': self.db().list['Tables'].remove(pdb)
            else: pdb.dropFields(['Hub','Spoke'] + edgedata,inverse=True)
            ndb.fillBlanks(); self.db('Edge').fillBlanks()
        except: self.errorLog('%s.loadPairwisePPI error' % self)
#########################################################################################################################
    def addPPI(self,hub,spoke,evidence=None,asdict=None,sym=True):    ### Adds a single PPI to the tables
        '''Adds a single PPI to the tables.'''
        try:### ~ [1] Add Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ndb = self.db('Node')
            if hub not in ndb.data(): ndb.addEntry({'Node':hub})
            if spoke not in ndb.data(): ndb.addEntry({'Node':spoke})
            ### ~ [2] Add Edge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppi = self.dict['PPI']
            edb = self.db('Edge')
            if asdict == None: asdict = edb and 'Evidence' in edb.fields()
            if edb:
                ekey = edb.makeKey({'Hub':hub,'Spoke':spoke})
                if ekey in edb.data():
                    entry = edb.data(ekey)
                    if evidence and evidence not in entry['Evidence']: entry['Evidence'] += '|%s' % evidence
                else:
                    if evidence: edb.addEntry({'Hub':hub,'Spoke':spoke,'Evidence':evidence})
                    else: edb.addEntry({'Hub':hub,'Spoke':spoke})
            ### ~ [3] Add to PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if hub not in ppi:
                if asdict: ppi[hub] = {}
                else: ppi[hub] = []
            if spoke not in ppi[hub]:
                if asdict: ppi[hub][spoke] = evidence
                else: ppi[hub].append(spoke)
            elif asdict and evidence and evidence not in ppi[hub][spoke]: ppi[hub][spoke] += '|%s' % evidence
            ### ~ [4] Symmetry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sym: self.addPPI(spoke,hub,evidence,asdict,False)
        except: self.errorLog('%s.addPPI error' % self)
#########################################################################################################################
    def evidence(self,hub,spoke,default='NA'):  ### Returns Evidence, if it can find it
        '''Returns Evidence, if it can find it.'''
        try: return self.ppi()[hub][spoke]
        except: pass
        edb = self.db('Edge')
        if not edb or not 'Evidence' in edb.fields(): return default
        ekey = '%s%s%s' % (hub,edb.info['Delimit'],spoke)
        try: return edb.data()[ekey]['Evidence']
        except: return default
#########################################################################################################################
    def edgeType(self,hub,spoke):
        try:
            edb = self.db('Edge')
            return edb.data()['%s%s%s' % (hub,edb.info['Delimit'],spoke)]['Type']
        except: return 'pp'   
#########################################################################################################################
    def nodeData(self,node,field,default=''):   ### Returns node data if found in Node table
        '''Returns node data if found in Node table.'''
        try: return self.db('Node').data()[node][field]
        except: return default
#########################################################################################################################
    def addCol(self,default=0,coldict={},G={},ckey='SAMPLE',edb=None,addcol=0):       ### Adds auto colours to Node Table
        '''Adds auto colours to Node Table.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ndb = self.db('Node')
            if not ndb: ndb = self.db().addEmptyTable('Node',['Node'],['Node'])
            if not edb: edb = ndb
            if not G: G = self.ppi()
            needkey = not (self.opt['ColByDeg'] or 'Col' in ndb.fields() or coldict)
            if needkey and ckey not in edb.fields():
                self.printLog('#COL','No key "%s" in Table %s' % (ckey,edb.info['Name']))
                return False
            if 'Col' not in ndb.fields(): ndb.addField('Col')
            ndb.fillBlanks()
            ## ~ [0a] ~ Add missing nodes to Node Table, with addcol if desired ~~~~~~~~~~~~~~~~~~~ ##
            for node in G:
                if node not in ndb.data():
                    ndb.data()[node] = {'Node':node}
                    if addcol: ndb.data()[node]['Col'] = addcol
            ### ~ [1] Colour by Degree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not coldict and self.opt['ColByDeg']:
                for node in G:
                    entry = ndb.data()[node]
                    d = len(G[entry['Node']])
                    if d == 0: entry['Col'] = 6
                    elif d == 1: entry['Col'] = 7
                    elif d == 2: entry['Col'] = 19
                    elif 5 > d > 2: entry['Col'] = 5
                    elif 10 > d > 4: entry['Col'] = 3
                    elif 50 > d > 9: entry['Col'] = 2
                    else: entry['Col'] = 1
                return
            ### ~ [2] ~ Update colours from dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Default' in coldict: default = coldict['Default']
            for entry in ndb.entries():
                try: entry['Col'] = coldict[edb.data()[entry['Node']][ckey]]
                except:
                    if not 'Col' in entry or not entry['Col']: entry['Col'] = default
        except: self.errorLog('rje_ppi.addCol() error')
#########################################################################################################################
    def ppiToXGMML(self,ppi={},name=None,npos={}):   ### Returns populated XGMML object
        '''Returns populated XGMML object.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not ppi: ppi = self.ppiG()            
            xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            if name: xgmml.info['Name'] = name
            ndb = self.db('Node')
            edb = self.db('Edge')
            #- Edge = Dictionary of edges between nodes {Type:{(source,target):Attributes}}
            #- EdgeAtt = Dictionary of edge attributes {Att:Type}
            #- Node = Dictionary of Nodes to be output {Node:Attributes}
            #- NodeAtt = Dictionary of node attributes {Att:Type}
            if npos: xgmml.dict['NodePos'] = npos    # Dictionary of node positions {Node:[x,y]}
            if ndb:
                for field in ndb.fields():
                    try:
                        for entry in ndb.entries(): entry[field]/1
                        xgmml.dict['NodeAtt'][field] = 'real'
                    except: xgmml.dict['NodeAtt'][field] = 'string'      
            else: xgmml.dict['NodeAtt'] = {'Degree':'real'}
            #self.deBug(xgmml.dict['NodeAtt'])
            if edb:
                for field in edb.fields():
                    try:
                        for entry in edb.entries(): entry[field]/1
                        xgmml.dict['EdgeAtt'][field] = 'real'
                    except: xgmml.dict['EdgeAtt'][field] = 'string'      
            else: xgmml.dict['EdgeAtt'] = {'Evidence':'string'}
            #self.deBug(xgmml.dict['EdgeAtt'])
            symG = self.symmetry(ppi,inplace=False)
            xgmml.dict['Node'] = {}
            for hub in symG:
                xgmml.dict['Node'][hub] = {'Degree':len(symG[hub])}  
                try: xgmml.dict['Node'][hub] = rje.combineDict(xgmml.dict['Node'][hub],ndb.data()[hub],overwrite=False)
                except: pass
            xgmml.dict['Edge'] = {}
            #self.deBug(ppi==symG)
            for hub in rje.sortKeys(ppi):
                for spoke in ppi[hub]:
                    etype = self.edgeType(hub,spoke)
                    if etype not in xgmml.dict['Edge']: xgmml.dict['Edge'][etype] = {}
                    if (spoke, hub) in xgmml.dict['Edge'][etype]: continue
                    xgmml.dict['Edge'][etype][(hub,spoke)] = {}
                    if not edb: continue
                    ekey = '%s%s%s' % (hub,edb.info['Delimit'],spoke)
                    try: xgmml.dict['Edge'][etype][(hub,spoke)] = edb.data()[ekey]
                    except: self.deBug('%s :: %s,%s?' % (edb.dataKeys(),hub,spoke))
            return xgmml
        except: self.errorLog('%s.ppiToXGMML() problem' % self)
#########################################################################################################################
    def xgmmlToPPI(self,xgmml): ### Reads in PPI data from an RJE XGMML object (not a file (yet)).
        '''Reads in PPI data from an RJE XGMML object (not a file (yet)).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['PPI'] = {}
            if not xgmml: return
            ndb = self.db('Node')
            if not ndb: ndb = self.db().addEmptyTable('Node',['Node'],['Node'])
            edb = self.db('Edge')
            if not edb: edb = self.db().addEmptyTable('Edge',['Hub','Spoke'],['Hub','Spoke','Type'])
            ### ~ [1] ~ Copy Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'].lower() in ['','none']: self.info['Name'] = xgmml.info['Name']
            self.dict['NodeAtt'] = xgmml.dict['NodeAtt']
            self.dict['EdgeAtt'] = xgmml.dict['EdgeAtt']
            nfields = ['Node']
            for node in xgmml.dict['Node']:
                self.dict['PPI'][node] = {}
                ndb.dict['Data'][node] = rje.combineDict({'Node':node},xgmml.dict['Node'][node])
                nfields = rje.listUnion(nfields,xgmml.dict['Node'][node].keys())
            nfields.sort()
            for field in nfields:
                if field not in ndb.fields(): ndb.list['Fields'].append(field)
            efields = ['Hub','Spoke']
            for etype in xgmml.dict['Edge']:
                for (hub,spoke) in xgmml.dict['Edge'][etype]:
                    if 'Evidence' in xgmml.dict['Edge'][etype][(hub,spoke)]: self.dict['PPI'][hub][spoke] = xgmml.dict['Edge'][etype][(hub,spoke)]['Evidence']
                    else: self.dict['PPI'][hub][spoke] = etype
                    self.dict['PPI'][spoke][hub] = self.dict['PPI'][hub][spoke] 
                    edb.dict['Data']['%s%s%s' % (hub,edb.info['Delimit'],spoke)] = rje.combineDict({'Hub':hub,'Spoke':spoke,'Type':etype},xgmml.dict['Edge'][etype][(hub,spoke)])
                    efields = rje.listUnion(efields,xgmml.dict['Edge'][etype][(hub,spoke)].keys())
            efields.sort()
            for field in efields:
                if field not in edb.fields(): edb.list['Fields'].append(field)
            ndb.fillBlanks(); edb.fillBlanks()
        except: self.errorLog('%s.xgmmlToPPI() problem' % self)          
#########################################################################################################################
    def springXGMML(self,outfile=None,G={},name=None,rjespring=True,rpng=True,prejuggle=(10,10),scaleforce=False):    ### Outputs spring-embedded layout to XGMML
        '''Outputs spring-embedded layout to XGMML.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Basefile') == 'None': self.warnLog('No Basefile set for %s' % self)
            if not outfile: outfile = '%s.xgmml' % self.info['Basefile']
            if not G: G = self.ppiG()
            #self.deBug(G)
            xgmml = self.ppiToXGMML(G,name)
            G = self.symmetry(G)
            #self.deBug(self.stat)
            if rjespring: xgmml.dict['NodePos'] = self.rjeSpringLayout(G,prejuggle=prejuggle,scaleforce=scaleforce)
            else: xgmml.dict['NodePos'] = forceDirectedLayout(G,damping=self.stat['Damping'],callobj=self,walltime=self.stat['Walltime'])
            self.addCol(G)
            xgmml.saveXGMML(outfile)
            if rpng: self.saveR(xgmml.dict['NodePos'],rje.baseFile(outfile),G)
        except: self.errorLog('%s.springXGMML() problem' % self)
#########################################################################################################################
    def rjeSpringLayout(self,G,prejuggle=(10,10),scaleforce=False): return rjeSpringLayout(G,damping=self.stat['Damping'],callobj=self,walltime=self.stat['Walltime'],nudge=self.stat['NudgeCyc'],prejuggle=prejuggle,scaleforce=scaleforce)
#########################################################################################################################
    def loadNPos(self,nfile): ### Loads node positions from file
        '''Loads node positions from file.'''
        try:### ~ [0] ~ Load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            inpos = rje.dataDict(self,nfile,['Node'],['X','Y'])
            npos = {}
            for node in rje.sortKeys(inpos):
                try: npos[node] = (float(inpos[node]['X']),float(inpos[node]['Y']))
                except: pass
            return npos
        except:  self.errorLog('%s.saveNPos() problem' % self); return {}
#########################################################################################################################
    def saveNPos(self,npos,nfile): ### Outputs node positions to file
        '''Outputs node positions to file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nhead = ['Node','X','Y']
            rje.delimitedFileOutput(self,nfile,nhead,rje_backup=True)
            ### ~ [1] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in rje.sortKeys(npos): rje.delimitedFileOutput(self,nfile,nhead,datadict={'Node':node,'X':npos[node][0],'Y':npos[node][1]})
        except:  self.errorLog('%s.saveNPos() problem' % self)
#########################################################################################################################
    def saveSVG(self,npos,basefile=None,G={},font=0,width=1600,ntype='ellipse',cutspace=True,backups=True,nodecol={}): ### Outputs interaction network to SVG file using rje_svg.
        '''Outputs interaction network to SVG file using rje_svg.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Need: npos,G,nodecol={},font=font,width=1600,height=1200,ntype=ntype / filename
            if not self.obj['SVG']: self.obj['SVG'] = rje_svg.SVG(self.log,self.cmd_list); self.obj['SVG'].setup()
            svg = self.obj['SVG']
            if not basefile: basefile = self.info['Basefile']
            if not G: G = self.ppiG()
            G = subGraph(G,G.keys())
            ndb = self.db('Node')
            ## ~ [0a] ~ SVG file/backup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pfile = '%s.svg' % basefile
            if backups: rje.backup(self,pfile)
            if os.path.exists(pfile): os.unlink(pfile)
            ## ~ [0b] ~ Node colours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not nodecol:
                if self.opt['ColByDeg']:
                    self.addCol(G=G)
                    ndb = self.db('Node')
                nodecol = {}
                #self.deBug(ndb); self.deBug(ndb.fields()); self.deBug(rje.sortKeys(ndb.data()))
                if ndb and 'Col' in ndb.fields():
                    for vj in rje.sortKeys(npos):
                        #self.deBug(vj); self.deBug(ndb.data()[vj])
                        try: nodecol[vj] = svg.col(int(ndb.data()[vj]['Col']))
                        except: nodecol[vj] = svg.col(48)
            ## ~ [0c] ~ Image width and height ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            spokenum = len(npos)
            if font <= 0: font = 20 - int(min(100,spokenum)/10.0)
            xmax = 0.0; ymax = 0.0
            for xy in npos.values(): xmax = max(xmax,xy[0]); ymax = max(ymax,xy[1])
            if xmax: yx = ymax / xmax
            else: yx = width
            if yx > 1:  # Invert XY if long and thin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for node in npos: npos[node] = (npos[node][1],npos[node][0])
                tmax = xmax; xmax = ymax; ymax = tmax
                yx = 1.0 / yx
            width = min(width,max(400,spokenum * font * 4))
            height = max(400,int(width * yx))
            ### ~ [1] ~ Save PPI SVG File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug(npos); self.deBug(G);
            svgtext = svg.networkPlot(npos,G,nodecol,font,width,height,ntype,cutspace)
            svg.svgFile(svgtext,pfile,width,height)
            svghtm = '<embed src="%s" width="%d" height="%d" type="image/svg+xml"' % (pfile,width,height)
            svghtm += ' pluginspage="http://www.adobe.com/svg/viewer/install/" />'
            return svghtm
            
        except: self.errorLog('%s.saveSVG() problem' % self); return False
#########################################################################################################################
    def saveTDT(self,npos,basefile=None,G={},backups=True):  ### Outputs interaction network to TDT file for R etc.
        '''Outputs interaction network to TDT file for R etc.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basefile: basefile = self.info['Basefile']
            if not G: G = self.ppiG()
            ndb = self.db('Node')
            ### ~ [1] ~ Save PPI File for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfile = '%s.ppi.tdt' % basefile
            phead = rje.sortKeys(G)
            if backups: rje.backup(self,pfile)
            if os.path.exists(pfile): os.unlink(pfile)
            headers = {}
            for h in phead: headers[h] = rje.replace(rje.replace(h,'[','-'),']','-')
            rje.delimitedFileOutput(self,pfile,phead,datadict=headers)
            for vi in rje.sortKeys(G):
                datadict = {}
                for vj in rje.sortKeys(G):
                    if vj in G[vi]: datadict[vj] = 1    #!# Add Edge strength (thickness) #!#
                    else: datadict[vj] = 0
                rje.delimitedFileOutput(self,pfile,phead,datadict=datadict)
            for xy in [0,1]:
                datadict = {}
                for vj in rje.sortKeys(G): datadict[vj] = npos[vj][xy]
                rje.delimitedFileOutput(self,pfile,phead,datadict=datadict)
            if self.opt['ColByDeg']:
                #self.deBug('ColByDeg')
                self.addCol(G=G)
                ndb = self.db('Node')
            if ndb and 'Col' in ndb.fields():
                datadict = {}
                for vj in rje.sortKeys(G):
                    try: datadict[vj] = ndb.data()[vj]['Col']
                    except: datadict[vj] = 0
                rje.delimitedFileOutput(self,pfile,phead,datadict=datadict)
            return True
        except: self.errorLog('%s.saveTDT() problem' % self); return False
#########################################################################################################################
    def saveR(self,npos,basefile=None,G={},cleantdt=True,backups=True):  ### Outputs interaction network to PNG file
        '''Outputs interaction network to PNG file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.saveTDT(npos,basefile,G,backups): raise ValueError
            ### ~ [3] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "pureppi" "%s"' % (self.info['RPath'],basefile)
            self.info['PathR'] = rje.makePath(os.path.abspath(rje.join(rje.split(sys.argv[0],os.sep)[:-1]+[''],os.sep)))
            rcall = rje.makePath('%srje.r' % self.info['PathR'],wholepath=True)
            if not os.path.exists(rcall):
                self.info['PathR'] = rje.makePath(os.path.abspath(rje.join(rje.split(sys.argv[0],os.sep)[:-2]+['libraries','r',''],os.sep)))
                rcall = rje.makePath('%srje.r' % self.info['PathR'],wholepath=True)
            #rcall = '%s/rje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rcall,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanext = ['r.tmp.txt']
            if cleantdt: cleanext.append('ppi.tdt')
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in cleanext:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog('%s.saveR() problem' % self)
#########################################################################################################################
    def rCall(self,rcmd,basefile):   ### Performs given R call for PNG generation
        '''Performs given R call for PNG generation.'''
        ### ~ [0] ~ Check for PNG and skip if found and not Force ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #if glob.glob('%s*png' % basefile) and not self.opt['Force']:
        #    return self.printLog('#PNG','%s*png exists - skipping (Force=F)' % basefile)
        ### ~ [1] ~ "Standard" R Call without using IRIDIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.opt['Iridis'] = False  #!# Make this part of main RJE Object code #!#
        #if not self.opt['Iridis']: rcmd = rcmd + ' 2>&1'
        self.printLog('#RCMD',rcmd)
        if not self.opt['Iridis']:
            problems = os.popen(rcmd).read()
            if problems: self.errorLog(problems,printerror=False)
            pngx = len(glob.glob('%s*png' % basefile))
            self.printLog('#PNG','%d PNG files made for %s' % (pngx,basefile))
            if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
            return problems
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
### End of SECTION II: PPI Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: Graph Theory Methods and other generic stuff added by me                                               #
#########################################################################################################################
def kCoreMax(G,mink=2,loops=False):    ### Returns highest k-core graph and k for graph G
    '''
    Returns highest k-core graph and k for graph G.
    >> G:dictionary of graph {node:[nodes]}
    >> mink:int = Min. k-core value
    '''
    bestK = []; bestk = 0; k = mink
    K = kCore(G,k)
    while K:
        bestK = K; bestk = k; k += 1
        K = kCore(subGraph(G,K,loops),k)
    return (bestK,bestk)
#########################################################################################################################       
def kCore(G,k=2,loops=False):   ### Returns k-core graph for k and graph G
    '''
    Returns k-core graph for k and graph G.
    >> G:dictionary of graph {node:[nodes]}
    >> k:int = k-core value
    '''
    V = rje.sortKeys(G) # Vertex list
    GK = subGraph(G,V,loops)
    K = []              # k-Core vertex list
    while V != K:
        V = rje.sortKeys(GK)
        K = V[0:]
        for v in V:
            if len(GK[v]) < k: K.remove(v)
        GK = subGraph(G,K,loops)
    return K
#########################################################################################################################       
def density(G,loops=False): ### Return density for graph G
    '''
    Return density for graph G. *Every* node must have an entry in G and no links outside G.
    >> G:dictionary of graph {node:[nodes]}
    >> loops:bool [False] = whether self-connecting loops are allowed
    '''
    if loops: Emax = len(G) * len(G)
    else: Emax = len(G) * (len(G) - 1)
    E = 0.0
    for v in G: E += len(G[v])
    if not E or not Emax: return 0.0
    return E / Emax
#########################################################################################################################       
def densitySubGraph(G,loops=False): ### Return density for graph G
    '''
    Return density for graph G. *Every* node must have an entry in G. Links to nodes outside G allowed
    >> G:dictionary of graph {node:[nodes]}
    >> loops:bool [False] = whether self-connecting loops are allowed
    '''
    if loops: Emax = len(G) * len(G)
    else: Emax = len(G) * (len(G) - 1)
    E = 0.0
    for v in G: E += len(rje.listIntersect(G[v],G.keys()))
    if not E or not Emax: return 0.0
    return E / Emax
#########################################################################################################################       
def subGraph(G,V,loops=False,addv=True):  ### Returns subgraph of G for given vertices V:
    '''
    Returns subgraph of G for given vertices V.
    >> G:dict of {node:[nodes]}
    >> V:list of nodes
    >> loops:bool = whether to allow loop [False]
    >> addv:bool = whether to add vertices from V that are missing from G [True]
    '''
    subgraph = {}
    for v in V:
        if v not in G: subgraph[v] = []
        else:
            subgraph[v] = rje.listIntersect(G[v],V)
            if not loops and v in subgraph[v]: subgraph[v].remove(v)
            subgraph[v].sort()
    return subgraph
#########################################################################################################################       
def fragment(G,sorted=True,callobj=None,minfrag=3):    ### Fragments the Graph into Lists of unconnected Nodes
    '''Fragments the Graph into Lists of unconnected Nodes.'''
    fragments = []; seen = []; vx = 0.0; vtot = len(G); remx = 0
    for v in rje.sortKeys(G):
        if callobj: callobj.progLog('\r#FRAG','Fragmenting Graph: %.2f%%' % (vx/vtot)); vx += 100.0
        if v in seen: continue
        frag = [v]; growth = [v]
        while growth:
            for n in G[growth.pop(0)]:
                if n in frag: continue
                frag.append(n); growth.append(n)
        seen += frag
        if len(frag) >= minfrag: fragments.append(frag)
        else: remx += len(frag)
    if callobj: callobj.printLog('\r#FRAG','Fragmented %s nodes into %s unconnected Graphs.' % (rje.integerString(vtot),rje.integerString(len(fragments))))
    if callobj: callobj.printLog('\r#REM','%s nodes removed during fragmentation (Fragments < %d members).' % (rje.integerString(remx),minfrag))
    if sorted: fragments = rje.sortListsByLen(fragments,rev=True)
    return fragments
#########################################################################################################################
def getSpokeXY(i,inum,start_ang=0,scale=1.0): ### Calculate X and Y from Angle ###
	spangle = 360 / inum
	A1 = spangle * (i-1) + start_ang
	if A1 > 360: A1 -= 360 
	A = A1
	while A >= 90: A -= 90 
	A = (A / 360.0) * (2 * math.pi)	# Convert to rads
	o = math.sin(A)
	a = math.cos(A)
	if A1 < 90: x = o; y = a
	elif 180 > A1 >= 90: x = a; y = -o
	elif 270 > A1 >= 180: x = -o; y = -a
	else: x = -a; y = o
	return [x*scale,y*scale]
#########################################################################################################################
def rjeStartPosNodeOrder(G,callobj=None):   ### Returns the order of the nodes to place for layouts
    '''Returns the order of the nodes to place for layouts.'''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        vdeg = {}
        for v in rje.sortKeys(G):
            if len(G[v]) not in vdeg: vdeg[len(G[v])] = []
            vdeg[len(G[v])].append(v)
        ### ~ [1] ~ Order nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        toplace = []; place1 = []; place2 = []; place3 = []
        for deg in rje.sortKeys(vdeg,True):
            if deg < 3: place3 += vdeg[deg]
            elif deg == 2: place2 += vdeg[deg]
            else: place1 += vdeg[deg]
        for placelist in [place1,place2,place3]:
            if not placelist: continue
            if not toplace: toplace.append(placelist.pop(0))
            while placelist:
                next = None; nx = 0
                for v in placelist:
                    vx = len(rje.listIntersect(G[v],toplace))
                    if vx > nx: next = v; nx = vx
                    #if v in neighbours: next = v; break
                if not next: next = placelist[0]
                toplace.append(next)
                placelist.remove(next)
        return toplace        
    except:
        if callobj: callobj.errorLog('rje_ppi.rjeStartPosNodeOrder() Error')
#########################################################################################################################
def OLDrjeStartPosNodeOrder(G,callobj=None):   ### Returns the order of the nodes to place for layouts
    '''Returns the order of the nodes to place for layouts.'''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        vdeg = {}
        for v in rje.sortKeys(G):
            if len(G[v]) not in vdeg: vdeg[len(G[v])] = []
            vdeg[len(G[v])].append(v)
        ### ~ [1] ~ Order nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        toplace = []; neighbours = []
        for deg in rje.sortKeys(vdeg,True):
            if not toplace: toplace.append(vdeg[deg].pop(0)); neighbours = G[toplace[0]][0:]
            while vdeg[deg]:
                next = None; nx = 0
                for v in vdeg[deg]:
                    vx = len(rje.listIntersect(G[v],toplace))
                    if vx > nx: next = v; nx = vx
                    #if v in neighbours: next = v; break
                if not next: next = vdeg[deg][0]
                toplace.append(next); neighbours += G[next][0:]
                vdeg[deg].remove(next)
        return toplace        
    except:
        if callobj: callobj.errorLog('rje_ppi.rjeStartPosNodeOrder() Error')
#########################################################################################################################
def rjeStartPos(v,G,npos,callobj=None,nudge=0.1,ignoresingles=False):     ### Returns starting position of new node, v
    '''
    Returns starting position of new node, v.
    >> v:str = Node identifier
    >> G:dict = PPI network as dictionary
    >> npos:dict = Dictionary of node positions for G that has been placed
    >> callobj:Object [None] = calling object for log messages etc.
    >> nudge:float [0.1] = Give node a random nudge with that stderr
    >> ignoresingles:bool [False] = Whether to ignore singleton spokes when determining placed spokes
    << new [x,y] positions of node v
    '''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        vpos = [0.0,0.0]
        placed = rje.listIntersect(npos.keys(),G[v])    # v interactors with position info
        if ignoresingles:
            for vj in placed[0:]:
                if len(G[vj]) < 2: placed.remove(vj)
        ### ~ [1] ~ Place node ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [1a] ~ No interactors -> Line up just beyond existing nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not G[v]:
            (miny,maxy) = (0.0,0.0)
            for (x,y) in npos.values():
                vpos[0] = max(vpos[0],x+1.0)
                miny = min(miny,y); maxy = max(maxy,y)
            if len(G) < 2: vpos[1] = (maxy-miny) / 2.0
            else: vpos[1] = miny + (rje.sortKeys(G).index(v) / float(len(G)-1)) * (maxy-miny)
        ## ~ [1b] ~ No PPI placed -> place in circle with radius 3.0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        elif not placed: vpos = getSpokeXY(rje.sortKeys(G).index(v),len(G),scale=3.0)
        ## ~ [1c] ~ One PPI placed -> place in circle around it with radius 0.5 ~~~~~~~~~~~~~~~~~~~ ##
        elif len(placed) == 1:
            vj = placed[0]
            vpos = getSpokeXY(G[vj].index(v),len(G[vj]),scale=0.5)
            vpos[0] += npos[vj][0]
            vpos[1] += npos[vj][1]
        ## ~ [1d] ~ Multiple PPI placed -> place in mean position of PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        else:
            for vj in placed:   
                vpos[0] += (npos[vj][0] / float(len(placed)))
                vpos[1] += (npos[vj][1] / float(len(placed)))
        ### ~ [2] ~ Add random nudging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if nudge > 0.0:
            vpos[0] += random.gauss(0.0,nudge)    # Nudge
            vpos[1] += random.gauss(0.0,nudge)    # Nudge
        return vpos
    except:
        if callobj: callobj.errorLog('rje_ppi.rjeStartPos Error')
#########################################################################################################################
def delimitedNPos(npos,callobj):
    ###TEMP###
    nfile = 'c:\\development\\npos.tdt'; nhead = ['Node','x','y']
    rje.delimitedFileOutput(callobj,nfile,nhead,rje_backup=True)
    for v in rje.sortKeys(npos): rje.delimitedFileOutput(callobj,nfile,nhead,datadict={'Node':v,'x':npos[v][0],'y':npos[v][1]})
#########################################################################################################################
def xyDis(ipos,jpos,mindis=0.1):    ### Calculate distance from two coordinates
    '''Calculate distance from two coordinates.'''
    dx = ipos[0] - jpos[0]
    dy = ipos[1] - jpos[1]
    if not dx and not dy: dx = mindis  # Cannot directly sit on top of each other - force some distance
    try: dij = math.sqrt(dx*dx + dy*dy)
    except: dij = max(rje.modulus(dx),rje.modulus(dy))   # Must be big!
    return dij
#########################################################################################################################
def rjeSpringLayout(G,damping=0.9,klim=0.01,callobj=None,cycles=1e6,walltime=0.1,nudge=1000,prepos={},prejuggle=(10,50),scaleforce=False):    ### Returns a dictionary of node positions for Graph G
    '''
    Returns a dictionary of node positions for Graph G using a simple spring-embedded force directed layout.
    >> G:dict = Graph as dictionary. (Will be replaced with dictionary of lists if not already.)
    >> damping:float [0.9] = Proportion of Force translated into Movement (1-D)
    >> klim:float [0.1] = Stop once kinetics (summed distance moved) below this threshold
    >> cycles:int [1000] = Number of cycles to perform (max)
    >> walltime:float [0.1] = Walltime in hours before killing
    >> nudge:int [1000] = Number of cycles before nudging random node
    >> prepos:dict {} = Starting positions for nodes.
    >> prejuggle:(int,int) = Tuple of (times,stepsize) for spring layout cycles by replacement of nodes at spoke average
    >> scaleforce:bool [False] = Whether to scale spring force inversely to node degree
    << npos: dict = Node positions
    '''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        npos = {}               # Node positions.
        if not G: return npos
        N = len(G)              # Number of nodes (vertices)
        V = rje.sortKeys(G)     # List of Nodes
        if N == 1: return {V[0]:[0.5,0.5]}
        G = subGraph(G,V)       # Make sure it's a dictionary of lists
        maxforce = 5.0          # Max repulsion and maximum distance beyond which there is no force (e.g. Nodes have no influence)
        start = time.time()     # Setup walltime
        userandom = True
        placement = [False,False]
        k = -1
        prejuggle = [prejuggle[0],prejuggle[1]] # Needs to be able to reassign.
        ## ~ [0a] ~ Setup initial node positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        toplace = []
        nudged = 'PostPlacement'
        if prepos: npos = prepos; nudged = []
        else:
            toplace = rjeStartPosNodeOrder(G,callobj)
            i = 3
            while len(toplace) > i and len(G[toplace[i]]) > 2: i += 1   #!# Test circle #!#
            for v in toplace[:i]: npos[v] = getSpokeXY(toplace.index(v),i,scale=max(2.0,math.sqrt(float(i))))
            toplace = toplace[i:]
        ## ~ [1a] Unconnected nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        extras = []
        for v in V:
            if not G[v]: extras.append(v)
        if len(extras) == len(G):   # All unconnected!
            npos = {}
            for v in extras: npos[v] = getSpokeXY(extras.index(v),N,scale=max(2.0,math.sqrt(float(N))))
            if callobj:
                callobj.printLog('\r#LAY','No connected nodes. Returning circle layout.')
                #callobj.deBug(G); callobj.deBug(npos)
            return npos
        for v in extras:
            G.pop(v)
            if v in toplace: toplace.remove(v)
            if v in npos: npos.pop(v)
        ### ~ [1] ~ Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding %s Nodes...' % (rje.iStr(N)),log=False)
            i = len(npos) - len(G)
            if not nudge: nudge = -1
            nloop = nudge; prevpull = []
            if prejuggle[0]: nloop = prejuggle[1]
            while i < cycles or toplace:
                ## ~ [1^] ~ Place additonal node if nodes waiting to be placed ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if toplace:
                    vi = toplace.pop(0)
                    if vi in npos: i+=1; continue
                    if not G[vi]: i+=1; continue
                    v = G[vi][0]
                    if len(G[vi]) == 1 and v in npos:
                        singles = []
                        for vj in G[v]:
                            if len(G[vj]) == 1: singles.append(vj)
                        if len(singles) == 1:
                            vj = singles[0]
                            npos[vj] = rjeStartPos(vj,{vj:[v],v:[vj]},{v:npos[v]},callobj)
                        else:
                            for vj in singles:
                                npos[vj] = getSpokeXY(singles.index(vj),len(singles),scale=random.gauss(1.0,0.1))
                                npos[vj][0] += npos[v][0]; npos[vj][1] += npos[v][1]
                    else: npos[vi] = rjeStartPos(vi,G,npos,callobj)
                elif nudged == 'PostPlacement':     # Redistribute singletons
                    nudged = []
                    for v in V:
                        if v not in G or len(G[v]) < 2: continue
                        singles = []
                        for vj in G[v]:
                            if len(G[vj]) == 1: singles.append(vj)
                        if len(singles) == 1:
                            vj = singles[0]
                            npos[vj] = rjeStartPos(vj,{vj:[v],v:[vj]},{v:npos[v]},callobj)
                        else:
                            for vj in singles:
                                npos[vj] = getSpokeXY(singles.index(vj),len(singles),scale=random.gauss(1.0,0.1))
                                npos[vj][0] += npos[v][0]; npos[vj][1] += npos[v][1]
                    #if callobj and callobj.opt['Test'] and rje.yesNo('Stop here?'):  break
                    continue
                elif prejuggle[0]:
                    if callobj: callobj.progLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; k = %s) [%d juggles] ' % (rje.iStr(i+1),rje.expectString(k),prejuggle[0]))
                    if nloop: nloop -= 1
                    else:   # End of prejuggle cycle
                        prejuggle[0] -= 1
                        if prejuggle[0]: nloop = prejuggle[1]
                        else: nloop = nudge
                        ## ~ Juggle nodes to mean of PPI ~ ##
                        juggle = []
                        juggleG = {}
                        for v in V:
                            if v in G:
                                juggleG[v] = []
                                for vj in G[v]:
                                    if len(G[vj]) > 1: juggleG[v].append(vj)
                                if len(juggleG[v]) > 1: juggle.append(v)
                        #if callobj: callobj.deBug(juggleG)
                        for j in range(prejuggle[1]):
                            juggle = rje.randomList(juggle)
                            #if callobj: callobj.deBug(juggle)
                            if not juggle: continue
                            for v in juggle:
                                npos[v] = rjeStartPos(v,G,npos,callobj,ignoresingles=True)
                                #if callobj: callobj.deBug('%s -> %s' % (v,npos[v]))
                            meandis = 0.0
                            for v in juggle:
                                for vj in juggleG[v]:
                                    meandis += (xyDis(npos[v],npos[vj])/(len(juggleG[v])*len(juggle)))
                            if meandis < 1.0: scalePos(npos,1/meandis)
                        nudged = 'PostPlacement'
                        continue
                elif not nloop: nloop = nudge
                else: nloop -= 1
                ## ~ [1a] ~ Calculate Forces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                i += 1      # Cycle
                k = 0.0     # Distance moved by nodes this round
                F = {}      # Dictionary of {node:(X,Y)} forces
                maxpull = (0.0,None)    # Node being pulled the most
                #if callobj: callobj.debug(npos)
                for vi in rje.sortKeys(npos):    # Loop through nodes
                    F[vi] = [0.0,0.0]
                    P = [0.0,0.0]
                    centroid = [0.0,0.0]
                    for vj in rje.sortKeys(npos):    # Loop through other nodes
                        if vi == vj: continue
                        ## ~ Distances ~ ##
                        dx = npos[vi][0] - npos[vj][0]
                        dy = npos[vi][1] - npos[vj][1]
                        if not dx and not dy: dx = 0.1  # Cannot directly sit on top of each other - force some distance
                        try: dij = math.sqrt(dx*dx + dy*dy)
                        except: dij = max(rje.modulus(dx),rje.modulus(dy))   # Must be big!
                        ## ~ Attraction ~ ##
                        if vj in G[vi]: # NB. Bigger distance = bigger pull!
                            if scaleforce:
                                F[vi][0] -= dx/len(G[vj])  # Proportion of Force in X (in direction of vi->vj)
                                F[vi][1] -= dy/len(G[vj])  # Proportion of Force in Y (in direction of vi->vj)
                                P[0] -= dx/len(G[vj])  # Proportion of Force in X (in direction of vi->vj)
                                P[1] -= dy/len(G[vj])  # Proportion of Force in Y (in direction of vi->vj)
                            else:
                                F[vi][0] -= dx  # Proportion of Force in X (in direction of vi->vj)
                                F[vi][1] -= dy  # Proportion of Force in Y (in direction of vi->vj)
                                P[0] -= dx  # Proportion of Force in X (in direction of vi->vj)
                                P[1] -= dy  # Proportion of Force in Y (in direction of vi->vj)
                            centroid[0] += (npos[vj][0] / float(len(G[vi])))
                            centroid[1] += (npos[vj][1] / float(len(G[vi])))
                        ## ~ Repulsion ~ ##
                        if dij > maxforce: continue             # Nodes too remote to influence each other
                        elif not G[vj] and dij > 1.0: continue  # Singleton nodes do not repel beyond distance 1.0
                        Fij = min(maxforce,1.0 / (dij ** 2))    # Avoid really close nodes flinging each other miles
                        F[vi][0] += (dx * Fij / dij)    # Proportion of Force in X (in direction of vj->vi)
                        F[vi][1] += (dy * Fij / dij)    # Proportion of Force in Y (in direction of vj->vi)
                    #if callobj: callobj.debug(F[vi])
                    if not nloop:
                        #try: P = math.sqrt(P[0] ** 2 + P[1] ** 2)
                        #except: P = 0.0
                        dx = npos[vi][0] - centroid[0]
                        dy = npos[vi][1] - centroid[1]
                        try: P = math.sqrt(dx ** 2 + dy ** 2)    # Use distance from centroid as indicator, not Pull
                        except: P = 0.0
                        if len(G[vi]) < 2: P -= 1.0
                        if P < 1.0: P = 0.0
                        if P > maxpull[0] and vi not in prevpull: maxpull = (P,vi)
    #            print 'F:', F
                ## ~ [1b] ~ Update positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    #            for vi in V:    # Loop through nodes
                    nbudge = [(1 - damping) * F[vi][0], (1 - damping) * F[vi][1]]
                    jump = 1.0
                    if userandom: jump = max(0.1,random.gauss(1.0,0.5))
                    for xy in [0,1]:
                        if nbudge[xy] >= 0: nbudge[xy] = min(maxforce,nbudge[xy])
                        else: nbudge[xy] = max(-maxforce,nbudge[xy])
                        npos[vi][xy] += jump * nbudge[xy]
                        k += rje.modulus(nbudge[xy])
    #            print 'F:', F
    #            print '>>>\n', i, k, npos, '\n<<<'
                #print i, k
                #if callobj: callobj.debug(npos)
                if k <= klim and not toplace: break
                if walltime > 0.0 and (time.time() - start) / 3600.0 > walltime and not toplace:
                    if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding: Walltime (%.2f hours) reached! ' % (walltime))
                    break
                if not toplace and not nloop and maxpull[1]:
                    v = maxpull[1]
                    nudged.append(v)
                    prevpull = nudged[-nudged.count(v):]
                    for v in prevpull:
                        singles = []; ipos = {}
                        for vj in G[v]:   # Move singleton interactors
                            if len(G[vj]) > 1: ipos[vj] = npos[vj]
                            singles.append(vj)
                        if len(G[v]) == 1: npos[v] = rjeStartPos(v,G,npos,callobj)   # Place in midpoint of current interactors
                        else: npos[v] = rjeStartPos(v,G,ipos,callobj)   # Place in midpoint of current interactors
                        if len(singles) == 1:
                            vj = singles[0]
                            npos[vj] = rjeStartPos(vj,{vj:[v],v:[vj]},{v:npos[v]},callobj)
                        else:
                            for vj in singles:
                                npos[vj] = getSpokeXY(singles.index(vj),len(singles),scale=random.gauss(1.5,0.1))
                                npos[vj][0] += npos[v][0]; npos[vj][1] += npos[v][1]
                        if nudged.count(v) / 2 != nudged.count(v) / 2.0: continue
                    if callobj: callobj.progLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; k = %s) >>> Nudged %s <<<     ' % (rje.iStr(i+1),rje.expectString(k),rje.join(rje.sortUnique(prevpull),'|')))
                    #if callobj and callobj.opt['Test'] and rje.yesNo('Stop here?'):  break
                    maxpull = (0.0,None)
                elif callobj: callobj.progLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; k = %s) ' % (rje.iStr(i+1),rje.expectString(k)))
        except KeyboardInterrupt:
            if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding interrupted. (%s Cyc; k = %s) ' % (rje.iStr(i+1),rje.expectString(k)))
        except: raise
        ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        elimit = max(20,int(math.sqrt(float(len(extras))))+1); odd = True
        while extras:
            if odd: addsingles = extras[:elimit+1]
            else: addsingles = extras[:elimit]
            odd = not odd
            eG = {}; epos = {}; 
            for v in addsingles: G[v] = []; eG[v] = []; extras.remove(v)
            for v in addsingles: epos[v] = rjeStartPos(v,eG,npos,callobj,nudge=0.0)
            for v in addsingles: npos[v] = epos.pop(v)
        minpos = npos[V[0]][0:]; maxpos = npos[V[0]][0:]; 
        for v in V: minpos[0] = min(minpos[0],npos[v][0]); minpos[1] = min(minpos[1],npos[v][1])
        for v in V:
            npos[v][0] -= minpos[0]; npos[v][1] -= minpos[1];
            maxpos[0] = max(maxpos[0],npos[v][0]); maxpos[1] = max(maxpos[1],npos[v][1])
        if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; %.1f X, %.1f Y)' % (rje.iStr(i+1),maxpos[0],maxpos[1]))
        #if callobj: callobj.debug(npos)
    except:
        if callobj: callobj.errorLog('rje_ppi.forceDirectedLayout Error')
    return npos
#########################################################################################################################
def OLDrjeSpringLayout(G,damping=0.9,klim=0.01,callobj=None,cycles=1e6,walltime=1.0,nudge=1000,placement={}):    ### Returns a dictionary of node positions for Graph G
    '''Returns a dictionary of node positions for Graph G using a simple spring-embedded force directed layout.'''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        npos = {}               # Node positions.
        N = len(G)              # Number of nodes (vertices)
        V = rje.sortKeys(G)     # List of Nodes
        G = subGraph(G,V)       # Make sure it's a dictionary of lists
        #damping = 0.5           # Proportion of Force translated into Movement (1-D)
        maxforce = 5.0          # Max repulsion and maximum distance beyond which there is no force (e.g. Nodes have no influence)
        #cycles = 1000           # Number of cycles to perform (max)
        #klim = 0.1              # Stop once kinetics (summed distance moved) below this threshold
        start = time.time()     # Setup walltime
        userandom = True
        placement = [False,False]
        ## ~ [0a] ~ Setup initial node positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        #for v in V:
        #    x = V.index(v)      # set up initial node positions randomly // make sure no 2 nodes are in exactly the same position
        #    npos[v] = [(x/int(math.sqrt(N))),(x - x/int(math.sqrt(N)) * int(math.sqrt(N)))]
        #if callobj: callobj.debug(npos)
        toplace = []
        if placement: npos = placement
        else:
            vdeg = {}
            for v in V:
                if len(G[v]) not in vdeg: vdeg[len(G[v])] = []
                vdeg[len(G[v])].append(v)
            #if callobj: callobj.deBug(vdeg)
            for deg in rje.sortKeys(vdeg,True): toplace += vdeg[deg]    # List of nodes ordered by degree (big to small)
            #if callobj: callobj.deBug(toplace)
            v = toplace.pop(0)
            npos[v] = [0.0,0.0]
            try: placing = rje.sortKeys(G[v])
            except: placing = G[v][0:]
            if v in placing: placing.remove(v)
            #if callobj: callobj.deBug(placing)
            placefirst = placing[:3]
            px = 1
            while placefirst:
                v = placefirst.pop(0); toplace.remove(v)
                npos[v] = getSpokeXY(px,3.0); px += 1
                for vj in G[v]:
                    if vj not in npos and vj not in toplace: toplace.append(vj)
        #print npos, '...', V
        ### ~ [1] ~ Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding %s Nodes...' % (rje.iStr(N)),log=False)
        i = len(npos) - len(G)
        if not nudge: nudge = -1
        nloop = nudge; prevpull = []
        nudged = 'PostPlacement'
        while i < cycles or toplace:
            ## ~ [1^] ~ Place additonal node if nodes waiting to be placed ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if toplace:
                if placing:
                    v = None
                    for vi in toplace:
                        if vi in placing: v = vi; break
                    if not v: v = toplace[0]
                    if len(G[v]) < 3 and not placement[1]: scalePos(npos,2.0)
                    if len(G[v]) < 2 and not placement[0]: scalePos(npos,3.0)
                    npos[v] = [0.0,0.0]
                    placed = rje.listIntersect(npos.keys(),G[v])
                    if len(placed) == 1:
                        vj = placed[0]
                        npos[v] = getSpokeXY(G[vj].index(v),len(G[vj]),scale=1.5)
                        npos[v][0] += npos[vj][0]
                        npos[v][1] += npos[vj][1]
                    else:
                        for vj in placed:   # Place in midpoint of current interactors then budge once
                            npos[v][0] += (npos[vj][0] / float(len(placed)))
                            npos[v][1] += (npos[vj][1] / float(len(placed)))
                    toplace.remove(v)
                else:
                    v = toplace.pop(0)    # No interactors to place
                    npos[v] = getSpokeXY(V.index(v),len(V),scale=3.0)
                for vj in G[v]:
                    if vj not in npos and vj not in toplace: toplace.append(vj)
            elif nudged == 'PostPlacement':     # Redistribute singletons
                nudged = []
                for vj in G:
                    if len(G[vj]) > 1 or not G[vj]: continue
                    v = G[vj][0]
                    vsingle = []
                    for vi in G[v]:
                        if len(G[vi]) == 1: vsingle.append(vi)
                    vsingle.sort()
                    npos[vj] = getSpokeXY(vsingle.index(vj),len(vsingle),scale=1.5)
                    npos[vj][0] += npos[v][0]
                    npos[vj][1] += npos[v][1]
                #if callobj and callobj.opt['Test'] and rje.yesNo('Stop here?'):  break
                continue
            elif not nloop: nloop = nudge
            else: nloop -= 1
            ## ~ [1a] ~ Calculate Forces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            i += 1      # Cycle
            k = 0.0     # Distance moved by nodes this round
            F = {}      # Dictionary of {node:(X,Y)} forces
            maxpull = (0.0,None)    # Node being pulled the most
            #if callobj: callobj.debug(npos)
            for vi in rje.sortKeys(npos):    # Loop through nodes
                F[vi] = [0.0,0.0]
                P = [0.0,0.0]
                for vj in rje.sortKeys(npos):    # Loop through other nodes
                    if vi == vj: continue
                    ## ~ Distances ~ ##
                    dx = npos[vi][0] - npos[vj][0]
                    dy = npos[vi][1] - npos[vj][1]
                    if not dx and not dy: dx = 0.1  # Cannot directly sit on top of each other - force some distance
                    try: dij = math.sqrt(dx*dx + dy*dy)
                    except: dij = max(rje.modulus(dx),rje.modulus(dy))   # Must be big!
                    ## ~ Attraction ~ ##
                    if vj in G[vi]:
                        F[vi][0] -= dx  # Proportion of Force in X (in direction of vi->vj)
                        F[vi][1] -= dy  # Proportion of Force in Y (in direction of vi->vj)
                        P[0] -= dx  # Proportion of Force in X (in direction of vi->vj)
                        P[1] -= dy  # Proportion of Force in Y (in direction of vi->vj)
                    ## ~ Repulsion ~ ##
                    if dij > maxforce: continue             # Nodes too remote to influence each other
                    elif not G[vj] and dij > 1.0: continue  # Singleton nodes do not repel beyond distance 1.0
                    Fij = min(maxforce,1.0 / (dij ** 2))    # Avoid really close nodes flinging each other miles
                    F[vi][0] += (dx * Fij / dij)    # Proportion of Force in X (in direction of vj->vi)
                    F[vi][1] += (dy * Fij / dij)    # Proportion of Force in Y (in direction of vj->vi)
                #if callobj: callobj.debug(F[vi])
                if not nloop: 
                    try: P = math.sqrt(P[0] ** 2 + P[1] ** 2)
                    except: P = 0.0
                    if P > maxpull[0] and vi not in prevpull: maxpull = (P,vi)
#            print 'F:', F
            ## ~ [1b] ~ Update positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
#            for vi in V:    # Loop through nodes
                nbudge = [(1 - damping) * F[vi][0], (1 - damping) * F[vi][1]]
                jump = 1.0
                if userandom: jump = max(0.1,random.gauss(1.0,0.5))
                for xy in [0,1]:
                    if nbudge[xy] >= 0: nbudge[xy] = min(maxforce,nbudge[xy])
                    else: nbudge[xy] = max(-maxforce,nbudge[xy])
                    npos[vi][xy] += jump * nbudge[xy]
                    k += rje.modulus(nbudge[xy])
#            print 'F:', F
#            print '>>>\n', i, k, npos, '\n<<<'
            #print i, k
            #if callobj: callobj.debug(npos)
            if k <= klim and not toplace: break
            if walltime > 0.0 and (time.time() - start) / 3600.0 > walltime and not toplace:
                if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding: Walltime (%.2f hours) reached! ' % (walltime))
                break
            if not toplace and not nloop and maxpull[1]:
                v = maxpull[1]
                nudged.append(v)
                prevpull = nudged[-nudged.count(v):]
                for v in prevpull:
                    if len(G[v]) == 1:
                        vj = G[v][0]
                        try:
                            dij = math.sqrt(npos[vj][0] ** 2 + npos[vj][1] ** 2)
                            npos[v][0] = (1.5 + dij) * npos[vj][0] / dij
                            npos[v][1] = (1.5 + dij) * npos[vj][1] / dij
                        except: pass
                        continue
                    else:
                        npos[v] = [0.0,0.0]
                        for vj in G[v]:   # Place in midpoint of current interactors 
                            npos[v][0] += (npos[vj][0] / float(len(G[v])))
                            npos[v][1] += (npos[vj][1] / float(len(G[v])))
                    if nudged.count(v) / 2 != nudged.count(v) / 2.0: continue
                    for vj in G[v]:   # Move singleton interactors
                        if len(G[vj]) > 1: continue
                        vsingle = []
                        for vi in G[v]:
                            if len(G[vi]) == 1: vsingle.append(vi)
                        vsingle.sort()
                        npos[vj] = getSpokeXY(vsingle.index(vj),len(vsingle),scale=1.5)
                        npos[vj][0] += npos[v][0]
                        npos[vj][1] += npos[v][1]
                if callobj: callobj.progLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; k = %s) >>> Nudged %s <<<     ' % (rje.iStr(i+1),rje.expectString(k),rje.join(rje.sortUnique(prevpull),'|')))
                #if callobj and callobj.opt['Test'] and rje.yesNo('Stop here?'):  break
                maxpull = (0.0,None)
            elif callobj: callobj.progLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; k = %s) ' % (rje.iStr(i+1),rje.expectString(k)))
        ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        minpos = npos[V[0]][0:]; maxpos = npos[V[0]][0:]; 
        for v in V: minpos[0] = min(minpos[0],npos[v][0]); minpos[1] = min(minpos[1],npos[v][1])
        for v in V:
            npos[v][0] -= minpos[0]; npos[v][1] -= minpos[1];
            maxpos[0] = max(maxpos[0],npos[v][0]); maxpos[1] = max(maxpos[1],npos[v][1])
        if callobj: callobj.printLog('\r#LAY','RJE Force directed spring embedding (%s Cyc; %.1f X, %.1f Y)' % (rje.iStr(i+1),maxpos[0],maxpos[1]))
        #if callobj: callobj.debug(npos)
    except:
        if callobj: callobj.errorLog('rje_ppi.forceDirectedLayout Error')
    return npos
#########################################################################################################################
def forceDirectedLayout(G,damping=0.5,klim=0.0001,callobj=None,walltime=1.0,killzone=True):     ### Returns a dictionary of node positions for Graph G
    '''Returns a dictionary of node positions for Graph G.'''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #if callobj: callobj.debug(G);
        start = time.time(); looplen = 10; looping = ['0.0'] * looplen; looped = ''
        velocity = {}; npos = {}; N = len(G); V = rje.sortKeys(G); COULOMB = 8.987e9
        maxsize = 100.0; minvel = klim
        for n in range(N):
            try: math.sqrt(3 * ((2 * maxsize) ** 2)); maxsize *= 10.0
            except: break
        maxsize = max(1e6,maxsize)
        for v in V:
            velocity[v] = [0.0,0.0] # set up initial node velocities to (0,0)
            x = V.index(v)          # set up initial node positions randomly // make sure no 2 nodes are in exactly the same position
            npos[v] = [10.0*(x/int(math.sqrt(N))),10*(x - x/int(math.sqrt(N)) * int(math.sqrt(N)))]
        #if callobj: callobj.debug(npos)
        total_kinetic_energy = klim + 1; t = 0.0
        ### ~ [1] ~ Loop to Equilibrium ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        spring_constant = 1.0; timestep = 0.5; nodemass = 1.0; spring_len = 1.0; max_limit = False; dead = []
        while total_kinetic_energy > klim:  # Loop until the simulation has stopped moving
            #x#print total_kinetic_energy, npos
            if walltime > 0.0 and (time.time() - start) / 3600.0 > walltime:
                if callobj: callobj.printLog('\r#LAY','Force directed spring embedding: Walltime (%.2f hours) reached! ' % (walltime))
                break
            looped = '%s:%s' % (looped,looping.pop(0)); looping.append('%.3e' % total_kinetic_energy)
            if rje.join(looping,':') in looped:
                if callobj: callobj.printLog('\r#LAY','Force directed spring embedding: Cycle (len %d) detected! ' % (looplen))
                break
            if callobj:
                if max_limit: callobj.progLog('\r#LAY','Force directed spring embedding: t = %.2e; Kinetic E = %.3e [!]' % (t,total_kinetic_energy))
                else: callobj.progLog('\r#LAY','Force directed spring embedding: t = %.2e; Kinetic E = %.3e ' % (t,total_kinetic_energy))
            total_kinetic_energy = 0.0      # running sum of total kinetic energy over all particles
            max_limit = False; t += timestep
            F = {}
            for vi in V:
                if vi in dead: continue
                F[vi] = [0.0, 0.0]     # running sum of total force on this particular node
                for vj in V:           # Coulomb_repulsion( this_node, other_node )
                    if vj == vi: continue
                    dx = npos[vi][0] - npos[vj][0]
                    dy = npos[vi][1] - npos[vj][1]
                    try: dij = math.sqrt(dx*dx + dy*dy)
                    except: dij = dx + dy
                    dij = max(dij,0.001)
                    #if callobj: callobj.debug('%s %s <-> %s %s = (%s,%s) = %s' % (vi,npos[vi],vj,npos[vj],dx,dy,dij)) 
                    try: Fij = spring_constant / (dij ** 2) #COULOMB / (dij ** 2)
                    except: Fij = 0.0
                    F[vi][0] += (dx * Fij / dij)
                    F[vi][1] += (dy * Fij / dij)
                    if vj in G[vi]:     # Hooke_attraction
                        #if callobj: callobj.debug('%s %s <-> %s %s = %.2f = %.2fN' % (vi,npos[vi],vj,npos[vj],dij,Fij)) 
                        disp = dij - spring_len
                        Fij = -spring_constant * disp
                        #if callobj: callobj.debug('%s <-> %s = Disp %.2f = %.2fN' % (vi,vj,disp,Fij))
                        F[vi][0] += (dx * Fij / dij)
                        F[vi][1] += (dy * Fij / dij)
                # Adjust velocity // without damping, it moves forever
                #if callobj: callobj.debug(F[vi])
                #if callobj: callobj.debug(npos[vi])
                velocity[vi][0] = (velocity[vi][0] + timestep * F[vi][0]) * (1 - damping)
                if rje.modulus(velocity[vi][0]) < minvel: velocity[vi][0] = 0.0
                npos[vi][0] = npos[vi][0] + (timestep * velocity[vi][0])
                #if callobj: callobj.debug(npos[vi])
                if npos[vi][0] >= maxsize: npos[vi][0] = maxsize; velocity[vi][0] = 0.0; max_limit = True
                if npos[vi][0] <= -maxsize: npos[vi][0] = -maxsize; velocity[vi][0] = 0.0; max_limit = True
                velocity[vi][1] = (velocity[vi][1] + timestep * F[vi][1]) * (1 - damping)
                if rje.modulus(velocity[vi][1]) < minvel: velocity[vi][1] = 0.0
                #if callobj: callobj.debug(npos[vi])
                npos[vi][1] = npos[vi][1] + (timestep * velocity[vi][1])
                if npos[vi][1] >= maxsize: npos[vi][1] = maxsize; velocity[vi][1] = 0.0; max_limit = True
                if npos[vi][1] <= -maxsize: npos[vi][1] = -maxsize; velocity[vi][1] = 0.0; max_limit = True
                if max(rje.modulus(npos[vi][0]),rje.modulus(npos[vi][1])) >= maxsize and killzone: dead.append(vi); continue
                try: total_kinetic_energy += nodemass * (velocity[vi][0]**2 + velocity[vi][1]**2)    # Mass x velocity
                except: pass
                #if callobj: callobj.debug(npos[vi])
            #x#if not rje.yesNo('Continue'): break
        if callobj and max_limit:
            callobj.printLog('\r#LAY','Force directed spring embedding: t = %.2e; Kinetic E = %.3e [!]' % (t,total_kinetic_energy))
            callobj.printLog('\r#MAX','Warning: maximum size limit reached. Check Layout')
        elif callobj: callobj.printLog('\r#LAY','Force directed spring embedding: t = %.2e; Kinetic E = %.3e. ' % (t,total_kinetic_energy))
        ### ~ [2] ~ Adjust positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        minpos = npos[V[0]][0:]; maxpos = npos[V[0]][0:]; 
        for v in V: minpos[0] = min(minpos[0],npos[v][0]); minpos[1] = min(minpos[1],npos[v][1])
        for v in V:
            npos[v][0] -= minpos[0]; npos[v][1] -= minpos[1];
            maxpos[0] = max(maxpos[0],npos[v][0]); maxpos[1] = max(maxpos[1],npos[v][1])
        ###TEMP###
        #nfile = 'c:\\development\\npos.tdt'; nhead = ['Node','x','y']
        #rje.delimitedFileOutput(callobj,nfile,nhead,rje_backup=True)
        #for v in V: rje.delimitedFileOutput(callobj,nfile,nhead,datadict={'Node':v,'x':npos[v][0],'y':npos[v][1]})
    except:
        if callobj: callobj.errorLog('rje_ppi.forceDirectedLayout Error')
    return npos
#########################################################################################################################       
def scalePos(npos,scale):   ### Multiplies all npos coordinates by scale
    '''Multiplies all npos coordinates by scale.'''
    for v in npos: npos[v][0] *= scale; npos[v][1] *= scale
#########################################################################################################################       
def randomGraph(N,E,name='N'):   ### Generate a random graph with N nodes and E edges
    '''Generate a random graph with N nodes and E edges.'''
    G = {}
    for n in range(N): G['%s%s' % (name,rje.preZero(n+1,N))] = []
    V = rje.sortKeys(G)
    maxE = N * (N - 1) / 2
    E = min(maxE,E)
    invert = E > (maxE / 2.0)
    if invert:
        for v in V:
            G[v] = V[0:]
            G[v].remove(v)
        E = maxE - E
    #print G
    ex = 0
    #print N, maxE, E, invert
    while ex < E:
        [vi,vj] = random.sample(range(N),2)
        vi = V[vi]; vj = V[vj]
        if invert:
            if vj not in G[vi]: continue
            G[vi].remove(vj); G[vj].remove(vi); ex += 1
        else:
            if vj in G[vi]: continue
            G[vi].append(vj); G[vj].append(vi); ex += 1
    return G
#########################################################################################################################       
### End of SECTION III: Graph Theory Methods                                                                            #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: OLD PPI Dictionary METHODS                                                                              #
#########################################################################################################################
def removeOrphans(G,callobj=None):  ### Removes all nodes without edges
    '''Removes all nodes without edges.'''
    try:### ~ [1] Purge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        N = len(G); V = rje.sortKeys(G)
        for v in V:
            if not G[v]: G.pop(v)
        if callobj: callobj.printLog('\r#PURGE','Purged unconnected nodes: %s -> %s nodes' % (rje.integerString(N),rje.integerString(len(G))))
        return G
    except:
        if callobj: callobj.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
def reciprocate(callobj=None,ppi={}):   ### This will go through and fill in any missing reciprocal interactions
    '''
    This will go through and fill in any missing reciprocal interactions.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    << ppi:returned dictionary with extra interactions added.
    '''
    try:
        ###~Work through each gene and check/add its reciprocal interactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppx = len(ppi)
        ppi = copy.copy(ppi)
        for hub in rje.sortKeys(ppi):
            for spoke in rje.sortKeys(ppi[hub]):
                if spoke not in ppi: ppi[spoke] = {}
                if hub not in ppi[spoke]: ppi[spoke][hub] = []
                for ppi_type in ppi[hub][spoke]:
                    if ppi_type not in ppi[spoke][hub]: ppi[spoke][hub].append(ppi_type)
            if callobj: callobj.log.printLog('\r#HUB','%d hubs => %d reciprocal hubs' % (ppx,len(ppi)),newline=False,log=False)
        if callobj: callobj.log.printLog('\r#HUB','%d hubs => %d reciprocal hubs' % (ppx,len(ppi)))
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.reciprocate()')
        else: raise
    return ppi
#########################################################################################################################
def combineTypes(callobj=None,ppi={},newtype=None): ### Combine PPI of different types with same hub-spoke pair
    '''
    Combine PPI of different types with same hub-spoke pair.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    >> newtype:str [None] = if given, all interactions will be combined into this new type
    << ppi:returned dictionary with extra interactions added.
    '''
    try:
        ###~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppi = copy.copy(ppi)        # Leave original dictionary unchanged
        ###~Merge Synonyms~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for hub in rje.sortKeys(ppi):
            for spoke in ppi[hub]:
                if newtype: ppi[hub][spoke] = [newtype]   # Only one interaction type
                else: ppi[hub][spoke] = [rje.join(rje.sortUnique(ppi[hub][spoke]),'')]
        if callobj: callobj.log.printLog('\r#COMBINE','PPI types combined for %d hubs' % (len(ppi)))
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.combineTypes()')
        else: raise
    return ppi
#########################################################################################################################
def mergeSynonyms(callobj=None,ppi={},synonyms={},logtext='synonyms'): ### This will go through and merge ppi for synonyms
    '''
    This will go through and merge ppi for synonyms.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    >> synonyms:dictionary of {gene:[synonyms]}
    >> logtext:str ['synonyms'] = description of synonym dictionary (could be matched by topology etc.)
    << ppi:returned dictionary with extra interactions added.
    '''
    try:
        ###~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppi = copy.copy(ppi)        # Leave original dictionary unchanged
        syn = copy.copy(synonyms)   # Leave original dictionary unchanged
        mx = 0
        ###~Merge Synonyms~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for hub in rje.sortKeys(syn):
            if hub not in syn:  continue        # Removed by earlier alias
            if hub not in ppi: ppi[hub] = {}
            #x#print hub, ppi[hub]
            for alias in syn[hub]:
                if alias == hub: continue
                if alias in syn: syn.pop(alias)     # Remove from potential hub list
                for replaceme in rje.sortKeys(ppi): # Replace alias with hub in ppi dictionary
                    if alias in ppi[replaceme]:
                        if hub not in ppi[replaceme]: ppi[replaceme][hub] = []
                        for itype in ppi[replaceme].pop(alias):
                            if itype not in ppi[replaceme][hub]: ppi[replaceme][hub].append(itype)
                if alias not in ppi: continue       # No interactions to merge
                mx += 1
                for spoke in ppi[alias]:
                    #x#print alias, spoke, hub, ppi[hub]
                    if spoke not in ppi[hub]: ppi[hub][spoke] = []
                    for itype in ppi[alias][spoke]:
                        if itype not in ppi[hub][spoke]: ppi[hub][spoke].append(itype)
                ppi.pop(alias)
            if callobj: callobj.log.printLog('\r#MERGE','%d %s merged to %d entries' % (len(synonyms),logtext,mx),newline=False,log=False)
        if callobj: callobj.log.printLog('\r#MERGE','%d %s merged to %d entries' % (len(synonyms),logtext,mx))
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.mergeSynonyms()')
        else: raise
    return ppi
#########################################################################################################################
def removeSingletons(callobj=None,ppi={},safelist=[]):  ### This will go through and remove singleton interactors
    '''
    This will go through and remove singleton interactors.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    >> safelist:list of gene IDs to keep singleton interactors of
    << ppi:returned dictionary with extra interactions added.
    '''
    try:
        ###~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppi = copy.copy(ppi)    # Leave original dictionary unchanged
        sx = 0
        ###~Remove Singletons~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for hub in rje.sortKeys(ppi):
            if hub not in ppi: continue
            if not ppi[hub]: ppi.pop(hub); continue
            if len(ppi[hub]) == 1:
                spoke = ppi[hub].keys()[0]
                if spoke not in safelist:
                    sx += 1
                    ppi.pop(hub)
                    if spoke in ppi and hub in ppi[spoke]:
                        try: ppi[spoke].pop(hub)
                        except: ppi[spoke].remove(hub)
                        if not ppi[spoke]: ppi.pop(spoke)
            if callobj: callobj.log.printLog('\r#REM','%d boring singletons removed' % (sx),newline=False,log=False)
        if callobj: callobj.log.printLog('\r#REM','%d boring singletons removed. %s hubs remain.' % (sx,rje.integerString(len(ppi))))
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.removeSingletons()')
        else: raise
    return ppi
#########################################################################################################################
def compressTopology(callobj=None,ppi={},safelist=[],matchtype=True):  ### Compresses nodes with same interactions into single entry
    '''
    This will go through and remove singleton interactors.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    >> safelist:list of gene IDs to ignore for compression
    >> matchtype:bool [True] = whether to check types of interaction in addition to just genes
    << ppi:returned dictionary with compressed genes (separated by /).
    '''
    try:
        ###~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppi = copy.copy(ppi)    # Leave original dictionary unchanged
        sharedtopology = {}     # Dictionary of {gene:[list of genes sharing topology]}
        ###~Identify shared topology~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for hub in rje.sortKeys(ppi):
            if hub in safelist: continue    # Do not check the safelist
            for check in rje.sortKeys(ppi):
                if check == hub: continue   # Don't merge with self, that's just silly!
                if check in safelist or rje.sortKeys(ppi[hub]) != rje.sortKeys(ppi[check]): continue    # Not the same
                mymatch = True
                if matchtype:
                    for spoke in ppi[hub]:
                        if rje.sortUnique(ppi[hub][spoke]) != rje.sortUnique(ppi[check][spoke]): mymatch = False
                if mymatch:
                    if hub not in sharedtopology: sharedtopology[hub] = []
                    if check not in sharedtopology: sharedtopology[check] = []
                    if hub not in sharedtopology[check]: sharedtopology[check].append(hub)
                    if check not in sharedtopology[hub]: sharedtopology[hub].append(check)
            if callobj: callobj.log.printLog('\r#TOP','%d shared topologies' % len(sharedtopology),newline=False,log=False)
        if callobj: callobj.log.printLog('\r#TOP','%d shared topologies' % len(sharedtopology))
        ###~Merge by shared topology~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ppi = mergeSynonyms(callobj,ppi,synonyms=sharedtopology,logtext='shared topologies')
        for hub in rje.sortKeys(ppi):   ## Now combine names using /
            for spoke in rje.sortKeys(ppi[hub]):
                if spoke in sharedtopology:
                    combo = rje.join(rje.sortUnique([spoke]+sharedtopology[spoke]),'/')
                    ppi[hub][combo] = ppi[hub].pop(spoke)
            if hub in sharedtopology: 
                combo = rje.join(rje.sortUnique([spoke]+sharedtopology[hub]),'/')
                ppi[combo] = ppi.pop(hub)
        if callobj: callobj.log.printLog('\r#TOP','%d PPI groups by shared topologies' % len(ppi))
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.compressTopology()')
        else: raise
    return ppi
#########################################################################################################################
def cytoscapeSave(callobj=None,basefile='ppi',ppi={},node={},edge={},bipolar=[]):  ### Outputs data in cytoscape format
    '''
    This will go through and remove singleton interactors.
    >> callobj:RJE Object calling method. (Used for error handling)
    >> ppi:dictionary of {hub:{spoke:[List of interaction types]}}
    >> basefile:str ['ppi'] = leader for files 
    >> node:dict = {Type:{gene:annotation}}
    >> edge:dict = {Type:{gene1:{gene2:{intType:annotation}}}}
    >> bipolar:list = List of genes to be output in both directions, else trim later PPI from earlier genes.
    '''
    try:
        ###~Output main SIF file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        sif = '%s.sif' % basefile
        SIF = open(sif,'w')
        for hub in rje.sortKeys(ppi):
            if hub not in ppi: continue    # Removed from earlier
            tmp_ppi = {}
            for spoke in rje.sortKeys(ppi[hub]):
                if spoke in bipolar and hub not in bipolar: continue    # Will be handled later
                ## Add to lines of output ##
                for itype in ppi[hub][spoke]:
                    if itype not in tmp_ppi: tmp_ppi[itype] = []
                    tmp_ppi[itype].append(spoke)
                ## Remove from future output ##
                if spoke in bipolar: continue
                if spoke in ppi and hub in ppi[spoke]: ppi[spoke].pop(hub)
                if spoke in ppi and not ppi[spoke]: ppi.pop(spoke)
            for itype in rje.sortKeys(tmp_ppi): SIF.write('%s %s %s\n' % (hub,itype,rje.join(tmp_ppi[itype])))
        SIF.close()
        if callobj: callobj.log.printLog('#SIF','%s output for %s hubs' % (sif,rje.integerString(len(ppi))))

        ###~Output Node Attributes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for type in node:
            nx = 0
            NODE = open('%s.%s.noa' % (basefile,type.lower()),'w')
            NODE.write('%s\n' % type)
            for gene in rje.sortKeys(node[type]):
                nx += 1
                NODE.write('%s = %s\n' % (gene,node[type][gene]))
            NODE.close()
            if callobj: callobj.log.printLog('#NOA','%s node attributes output to %s.%s.noa' % (rje.integerString(nx),basefile,type.lower()))

        ###~Output Node Attributes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        for type in edge:
            EDGE = open('%s.%s.eoa' % (basefile,type.lower()),'w')
            EDGE.write('%s\n' % type)
            for hub in rje.sortKeys(edge[type]):
                for spoke in rje.sortKeys(edge[type][hub]):
                    for itype in rje.sortKeys(edge[type][hub][spoke]): EDGE.write('%s (%s) %s = %s\n' % (hub,itype,spoke,edge[type][hub][spoke][itype]))
            EDGE.close()
    except:
        if callobj: callobj.log.errorLog('Major problem during rje_ppi.cytoscapeSave()')
        else: raise
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        PPI(mainlog,cmd_list).run()
        #PPI(mainlog,cmd_list).test()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
