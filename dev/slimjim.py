#!/usr/bin/python

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
Module:       SLiMJIM
Description:  Short Linear Motif Jolly Image Maker
Version:      2.0
Last Edit:    08/09/10
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program is designed to be used in conjunction with the slimjim.r R script to generate useful images for helping
    the interpretation of SLiMFinder runs. These can be embedded into a series of linked html pages for ease of navigat-
    ion and data exploration.

    Alternatively, SLiMJIM can be used to generate single images for single SLiMFinder runs, without the extra linking
    and PPI information.

Commandline:
    ### ~ SLiMFinder Run Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    resfile=FILE    : Main SLiMFinder results file []
    resdir=PATH     : Path to extra SLiMFinder results []
    alndir=PATH     : Path to orthologue alignments for those missing mapping.fas []
    pairwise=FILE   : Pingu Pairwise output file []
    name=X          : Name of analysis ['SLiMFinder Human PPI Analysis']
    suffix=X        : Suffix added to gene symbols for dataset name [.ppi]
    ensgopath=PATH  : Path to EnsGO files   (!!! Restricted to Humans Currently !!!)
    singlesf=FILE   : Path to *.slimdb file for which to generate graphics for single run []
    ### ~ Specific Analysis Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hublist=LIST    : List of hub datasets to look at (*.ppi). Can use wildcards. []
    runlist=LIST    : Optional list to limit results to certain run IDs []
    slimlist=LIST   : Optional list to limit results to certain SLiMs []
    runtypes=LIST   : List of run types to generate output for (if multiple runs) []
    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hubs=T/F        : Whether to generate output for Hubs [True]
    domains=T/F     : Whether to generate output for Domains [True]
    spokes=T/F      : Whether to generate output for Spokes [True]
    motifs=T/F      : Whether to generate output for Motifs [True]
    mains=T/F       : Whether to generate main analysis pages [True]
    makejim=T/F     : Whether to generate data and PNG files [False]
    jimdir=PATH     : Output directory for generated PNG files [./]
    mapdir=PATH     : Output directory for intermediate composite *.mapping.fas files [sjmap/]
    iridis=T/F      : Whether running on IRIDIS - fork out R jobs [False]
    ### ~ HTML Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    makehtml=T/F    : Whether to generate linked HTML output files [False]
    htmlall=T/F     : Whether to make HTML for all datasets or only those identified with hublist [False]
    htmldir=PATH    : Output path for HTML files [html/]
    htmlcss=X       : File to be used for cascading style sheet (in */resource/)
    htmlpng=X       : Graphic for corner of header frames - links back to Home ['SBS_100.png']
    unihtml=T/F     : Generate missing linked HTML DAT files [True]
    unipath=PATH    : Path to indexed real UniProt DAT files []
    unifake=PATH    : Path to indexed UniFake DAT file []
    unireal=LIST    : Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]
    tableonly=T/F   : Only make delimited text versions of tables, not all HTML files [False]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_dismatrix_V2, rje_ensembl, rje_genemap, rje_go, rje_iridis, rje_slim, rje_tree, rje_seq,
    rje_sequence, rje_uniprot, rje_zen
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_dismatrix_V2, rje_ensembl, rje_genemap, rje_go, rje_slim, rje_tree, rje_seq, rje_sequence, rje_uniprot, rje_zen
import slimjim_V1
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Altered basis of all protein naming to be gene-based. Added more visualisations and UniProt data.
    # 0.2 - Tidied up management of Hubs and Spokes.
    # 0.3 - Added SingleSF mode and tidied some of R code.
    # 0.4 - Added use of rje_iridis and code from sfmap2png.
    # 1.0 - Basic functionality all present. Fixed ordering bug.
    # 1.1 - Added tableonly=T/F option.
    # 2.0 - Updated to handle multiple runs and domain datasets.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [?] : Add extraction of specific motifs only
    # [ ] : Add <META HTTP-EQUIV="Description" NAME="Description" CONTENT="xxx">
    # [ ] : Add <META HTTP-EQUIV="Keywords" NAME="Keywords" CONTENT="xxx">
    # [ ] : Remove norobot instructions once published
    # [ ] : Add domain tables to front page in addition to genes and motifs.
    # [ ] : Add different PPI types as <H3> sets following main gene/domain interactors.
    # [ ] : Add domain pages, similar in general design to gene hub tables.
    # [ ] : Add domains contained by gene to gene pages.
    # [ ] : Add genes containing domain to domain pages.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('SLiMJIM', '2.0', 'October 2009', '2008')
    description = 'Short Linear Motif Jolly Image Maker'
    author = 'Dr Richard J. Edwards.'
    comments = ['An installation of R is necessary for this program to work.',
                'This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
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
            if rje.yesNo('Show GeneMap commandline options?'): out.verbose(-1,4,text=rje_genemap.__doc__)
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMJIM Class                                                                                           #
#########################################################################################################################
class SLiMJIM(slimjim_V1.SLiMJIM):     
    '''
    SLiMJIM Class. Author: Rich Edwards (2008).

    Info:str
    - Name = Name of analysis ['SLiMFinder Human PPI Analysis']
    - ResFile = Main SLiMFinder results file []
    - ResDir = Path to extra SLiMFinder results []
    - AlnDir = Path to orthologue alignments for those missing mapping.fas []
    - Pairwise = Pingu Pairwise output file []
    - EnsGOPath = Path to EnsGO files 
    - JIMDir = Output directory for generated PNG files [./]
    - HTMLDir = Output path for HTML files [html/]
    - HTMLCSS = File to be used for cascading style sheet (in */resource/)
    - HTMLPNG = Graphic for corner of header frames - links back to Home ['SBS_100.png']
    - MapDir = Output directory for intermediate composite *.mapping.fas files [sjmap/]
    - SingleSF = Path to *.slimdb file for which to generate graphics for single run []
    - Suffix = Suffix added to gene symbols for dataset name [.ppi]
    - UniFake = Path to indexed UniFake DAT file []
    
    Opt:boolean
    - HTMLAll = Whether to make HTML for all datasets or only those identified with hublist [False]
    - MakeHTML = Whether to generate linked HTML output files [False]
    - MakeJIM = Whether to generate data and PNG files [False]
    - Domains = Whether to generate output for Domains [True]
    - Hubs = Whether to generate output for Hubs [True]
    - Iridis = Whether running on IRIDIS - fork out R jobs [False]
    - Spokes = Whether to generate output for Spokes [True]
    - Motifs = Whether to generate output for Motifs [True]
    - Mains = Whether to generate main analysis pages [True]
    - TableOnly = Only make delimited text versions of tables, not all HTML files [False]

    Stat:numeric

    List:list
    - Datasets = List of datasets to analyse (read from ResDir using HubList)
    - HubList = List of hub datasets to look at (*.ppi). Can use wildcards. []
    - RunList = Optional list to limit results to certain run IDs []
    - SlimList = Optional list to limit results to certain SLiMs []
    - OccHead = List of headers read from Occ Results []
    - RunTypes = List of run types to generate output for (if multiple runs) []
    - SumHead = List of headers read from Summary Results []
    - UniReal = Real UniProt data to add to UniFake output (Empty=None) [AC,GN,RC,RX,CC,DR,PE,KW,FT]

    Dict:dictionary
    - EnsLoci = Links for genes to EnsLoci
    - GO = Dictionary of Gene-GO Mappings {Gene:{[GO]}}
    - PPI = Stores PPI data loaded from Pairwise file
    - SeqMap = Maps EnsLoci to genes for PPI mapping
    - SLiMFinder = Stores summary results for SLiMFinder run (selected datasets/runs only)
    - SlimHubs = List of hubs containing each SLiM
    - SlimSpokes = List of spokes containing each SLiM
    - SpokeSlims = List of slims contained by each spoke {Gene:[Patterns]}

    Obj:RJE_Objects
    - GeneMap = rje_genemap object for links to other databases for HTML making    
    - GO = rje_go.GO Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['ResFile','ResDir','Pairwise','JIMDir','HTMLDir','HTMLCSS','UniFake','AlnDir','HTMLPNG',
                         'Suffix','EnsGOPath','SingleSF','PyPath','iIni','MapDir']
        self.optlist = ['MakeHTML','MakeJIM','HTMLAll','UniHTML','Hubs','Spokes','Motifs','Mains','Iridis','RjePy','TableOnly','Domains']
        self.statlist = ['SubSleep']
        self.listlist = ['Datasets','HubList','RunList','SlimList','OccHead','SumHead','Hosts','SubJobs','UniReal','RunTypes']
        self.dictlist = ['PPI','EnsLoci','SlimHubs','SlimSpokes','SpokeSlims','GO','Running']
        self.objlist = ['GeneMap','GO']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'HTMLDir':rje.makePath('html/'),'HTMLCSS':'bioinf.css','Name':'SLiMFinder Human PPI Analysis',
                      'JIMDir':rje.makePath('html/'),'HTMLPNG':'SBS_100.png','Suffix':'.ppi',
                      'MapDir':rje.makePath('sjmap/')})
        for o in ['Hubs','Spokes','Motifs','Mains','UniHTML','Domains']: self.opt[o] = True
        self.list['UniReal'] = string.split('AC,GN,RC,RX,CC,DR,PE,KW,FT',',')
        ### Iridis ###
        self.setInfo({'PyPath':rje.makePath('/rhome/re1u06/Serpentry/')})
        self.setStat({'SubSleep':1.0})
        self.setOpt({'RjePy':True})
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
                self._cmdReadList(cmd,'info',['HTMLCSS','Name','HTMLPNG','Suffix'])
                self._cmdReadList(cmd,'file',['ResFile','Pairwise','SingleSF'])
                self._cmdReadList(cmd,'path',['ResDir','JIMDir','HTMLDir','UniFake','AlnDir','EnsGOPath','MapDir'])
                self._cmdReadList(cmd,'list',['HubList','RunList','SlimList','UniReal','RunTypes'])
                self._cmdReadList(cmd,'opt',['MakeHTML','MakeJIM','HTMLAll','UniHTML','Hubs','Spokes','Motifs','Mains',
                                             'Iridis','TableOnly','Domains'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.opt['Iridis']: self._iridisCmdList()               # Allows easy incorporation when inheriting class
#########################################################################################################################
    ### <2> ### Main Pipeline Class Methods                                                                             #
#########################################################################################################################
    def run(self):  ### Main run method.
        '''Main run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile(self.info['SingleSF']): return self.singleSF()
            self.setup()
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Iridis']: self.iridisRun()
            elif self.opt['TableOnly']: self.makeTables()
            else:
                if self.opt['MakeJIM']: self.makeJIM()
                if self.opt['MakeHTML']: self.makeHTML()
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):  ### Loads and organises relevant data.                                                        #V2.0
        '''Loads and organises relevant data.'''
        try:### ~ [0] ~ Setup Data Storage Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ General Data dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['SLiMFinder'] = {}    # *All* results for specified RunIDs
            self.dict['PPI'] = {}           # Dictionary of {HubGene:{SpokeGene:{Data}}}
            self.dict['SeqMap'] = {}        # Dictionary of {SpokeSeq:[Genes]}
            self.dict['EnsLoci'] = {}       # Dictionary of {SpokeGene:SpokeSeq}
            ## ~ [0b] ~ Results data lists and dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Hubs'] = []          # List of all hub genes with datasets for processing
            self.list['Spokes'] = []        # List of all spoke genes contained in one or more processed datasets
            self.dict['SlimHubs'] = {}      # Dictionary of {Pattern:{Hub:[SpokeSeq]}}
            self.dict['SlimSpokes'] = {}    # Dictionary of {Pattern:{SpokeSeq:[Hubs]}}
            self.dict['HubSlims'] = {}      # Dictionary of {Hub:[SlimList]}
            self.dict['SpokeSlims'] = {}    # Dictionary of {SpokeSeq:[SlimList]}

            ### ~ [1] ~ Generic Data Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Read Mapping Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['GeneMap'] = rje_genemap.GeneMap(self.log,self.cmd_list)
            self.obj['GeneMap'].run()
            ## ~ [1b] Read GO Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            self.obj['GO'].readGO()
            self.obj['GO'].mapEnsGO()
            self.dict['GO'] = self.obj['GO'].dict['EnsGO']
            ## ~ [1c] Read PPI from Pairwise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#PPI','Processing %s: 0.00%%' % (self.info['Pairwise']))
            pdata = rje.dataDict(self,self.info['Pairwise'],['Hub','Spoke'],datakeys=['Hub','Spoke','SpokeUni','SpokeSeq','Evidence'],headers=['Hub','HubUni','Spoke','SpokeUni','SpokeSeq','Evidence'])
            (px,ptot) = (0.0,len(pdata))
            for pkey in pdata.keys()[0:]:
                self.progLog('\r#PPI','Processing %s: %.2f%%' % (self.info['Pairwise'],(px/ptot))); px += 100.0
                idata = pdata.pop(pkey)
                p1 = idata.pop('Hub')
                if p1 not in self.dict['PPI']: self.dict['PPI'][p1] = {}
                p2 = idata.pop('Spoke')
                self.dict['PPI'][p1][p2] = idata
                spokeseq = idata['SpokeSeq']
                if spokeseq in ['','-']: self.dict['EnsLoci'][p2] = '-'; continue
                if spokeseq not in self.dict['SeqMap'].keys(): self.dict['SeqMap'][spokeseq] = []
                if p2 not in self.dict['SeqMap'][spokeseq]: self.dict['SeqMap'][spokeseq].append(p2)
                self.dict['EnsLoci'][p2] = spokeseq
            self.printLog('\r#PPI','Processing %s: %s hubs; %s PPI read.' % (self.info['Pairwise'],rje.integerString(len(self.dict['PPI'])),rje.integerString(ptot/2)))
            for spokeseq in rje.sortKeys(self.dict['SeqMap']):
                if len(self.dict['SeqMap'][spokeseq]) > 1: self.printLog('#DUP','%s => %s' % (spokeseq,self.dict['SeqMap'][spokeseq]))

            ### ~ [2] Load SLiMFinder results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdata = rje.dataDict(self,self.info['ResFile'],['Dataset','RunID','Rank','Pattern'],'All',getheaders=True)
            self.list['SumHead'] = sdata.pop('Headers')
            (sx,stot,mx) = (0.0,len(sdata),0)
            for skey in rje.sortKeys(sdata):
                self.progLog('\r#SLIM','Reading SLiMFinder results: %.2f%%' % (sx/stot)); sx += 100.0
                mdata = sdata.pop(skey)
                dataset = mdata['Dataset']
                pattern = mdata['Pattern']
                if self.info['Suffix'] and string.split(dataset,'.')[-1] != self.info['Suffix']: continue
                hub = string.replace(dataset,self.info['Suffix'],'')
                if self.list['RunList'] and mdata['RunID'] not in self.list['RunList']: continue
                if self.list['HubList'] and hub not in self.list['HubList']: continue
                if hub not in self.list['Hubs']: self.list['Hubs'].append(hub)
                if self.list['SlimList'] and pattern not in self.list['SlimList']: continue
                if pattern != '-':
                    if pattern not in self.dict['SlimHubs']: self.dict['SlimHubs'][pattern] = {}
                    if hub not in self.dict['SlimHubs'][pattern]: self.dict['SlimHubs'][pattern][hub] = mdata['Sig']
                if dataset not in self.dict['SLiMFinder']: self.dict['SLiMFinder'][dataset] = []
                self.dict['SLiMFinder'][dataset].append(mdata); mx += 1
            self.list['Hubs'].sort()
            self.printLog('\r#SLIM','Read SLiMFinder results: %s results from %s lines.' % (rje.integerString(mx),rje.integerString(stot)))
            ## ~ [2a] ~ Sort SLiMFinder results by Rank ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (sx,stot) = (0.0, len(self.dict['SLiMFinder']))
            for dataset in rje.sortKeys(self.dict['SLiMFinder']):
                self.progLog('\r#SLIM','Sorting SLiMFinder results by Rank: %.2f%%' % (sx/stot))
                sdata = []
                for mdata in self.dict['SLiMFinder'][dataset]:
                    i = 0
                    while i < len(sdata) and string.atoi(sdata[i]['Rank']) < string.atoi(mdata['Rank']): i += 1
                    sdata.insert(i,mdata)
                self.dict['SLiMFinder'][dataset] = sdata[0:]
            self.printLog('\r#SLIM','Sorting SLiMFinder results by Rank complete.')
            ## ~ [2b] ~ Make SpokeList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (sx,stot) = (0.0, len(self.list['Hubs']))
            for hub in self.list['Hubs']:
                self.progLog('\r#SPOKE','Identifying spokes: %.2f%%' % (sx/stot)); sx += 100.0
                if hub not in self.dict['PPI']: continue
                for spoke in self.dict['PPI'][hub]:
                    if self.dict['EnsLoci'][spoke] != '-' and spoke not in self.list['Spokes']: self.list['Spokes'].append(spoke)
            self.list['Spokes'].sort()
            self.printLog('\r#SPOKE','Identified %s spoke genes.' % (rje.integerString(len(self.list['Spokes']))))

            ### ~ [3] ~ Load Occurrence Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ox,otot) = (0.0,len(self.list['Hubs']))
            for hub in self.list['Hubs']:
                dataset = '%s%s' % (hub,self.info['Suffix'])
                occfile = '%s%s.occ.csv' % (self.info['ResDir'],dataset)
                self.progLog('\r#OCC','Loading Occurrences for %s hubs: %.2f%%' % (rje.integerString(otot),(ox/otot))); ox += 100.0
                if not os.path.exists(occfile): continue
                occ = rje.dataDict(self,occfile,['Pattern'],['Seq'],lists=True)
                self.dict['HubSlims'][hub] = []
                for pattern in rje.sortKeys(occ):     
                    if self.list['SlimList'] and pattern not in self.list['SlimList']: continue
                    self.dict['HubSlims'][hub].append(pattern)
                    if pattern not in self.dict['SlimSpokes']: self.dict['SlimSpokes'][pattern] = {}
                    #if pattern not in self.dict['SlimHubs']: self.dict['SlimHubs'][pattern] = {}
                    #self.dict['SlimHubs'][pattern][hub] = occ[pattern]['Seq'][0:]
                    for spokeseq in occ[pattern]['Seq']:
                        if spokeseq not in self.dict['SlimSpokes'][pattern]: self.dict['SlimSpokes'][pattern][spokeseq] = []
                        if hub not in self.dict['SlimSpokes'][pattern][spokeseq]: self.dict['SlimSpokes'][pattern][spokeseq].append(hub)
                        if spokeseq not in self.dict['SpokeSlims']: self.dict['SpokeSlims'][spokeseq] = []
                        if pattern not in self.dict['SpokeSlims'][spokeseq]: self.dict['SpokeSlims'][spokeseq].append(pattern)
            self.printLog('\r#OCC','Loaded Occurrences for %s hubs: %s SLiMs; %s Spokes.' % (rje.integerString(otot),rje.integerString(len(self.dict['SlimSpokes'])),rje.integerString(len(self.dict['SpokeSlims']))))

            ### ~ [4] ~ Special extra data correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'] == 'PTPRD SLiMFinder':
                for hub in self.list['Hubs']:
                    if 'PTPRD' not in self.dict['PPI'][hub]: self.dict['PPI'][hub]['PTPRD'] = {'SpokeUni':'P23468','SpokeSeq':'PTPRD_HUMAN__P23468','Evidence':'Laavanya:Misc'}  
                    if hub not in self.dict['PPI']['PTPRD']:
                        try: self.dict['PPI']['PTPRD'][hub] = {'SpokeSeq':self.dict['EnsLoci'][hub],'Evidence':'Laavanya:Misc'}
                        except: self.dict['PPI']['PTPRD'][hub] = {'Evidence':'Laavanya:Misc'}
                self.printLog('#PTPRD','Extra PTPRD interaction data added.')
            if self.info['Name'] == 'PTPRD Hub':
                for spokeseq in self.dict['SlimSpokes']:
                    spoke = self.gene(spokeseq)
                    if 'PTPRD' not in self.dict['PPI'][spoke]: self.dict['PPI'][spoke]['PTPRD'] = {'SpokeUni':'P23468','SpokeSeq':'PTPRD_HUMAN__P23468','Evidence':'Laavanya:Misc'}  
                    if spoke not in self.dict['PPI']['PTPRD']:
                        try: self.dict['PPI']['PTPRD'][spoke] = {'SpokeSeq':self.dict['EnsLoci'][spoke],'Evidence':'Laavanya:Misc'}
                        except: self.dict['PPI']['PTPRD'][spoke] = {'Evidence':'Laavanya:Misc'}
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def go(self,id=None): return self.obj['GO'].go(id)
#########################################################################################################################
    def addGeneGO(self,gene,go):    ### Adds gene-go relationship (and all parents) to self.dict['GO']
        '''Adds gene-go relationship (and all parents) to self.dict['GO'].'''
        if gene not in self.dict['GO']: self.dict['GO'][gene] = []
        if go in self.dict['GO'][gene]: return
        self.dict['GO'][gene].append(go)
        for parent in self.obj['GO'].parents(go): self.addGeneGO(gene,parent)
#########################################################################################################################
    def gene(self,prot):    ### Returns Gene associated with given protein ID
        '''Returns Gene associated with given protein ID.'''
        try:
            for gene in self.dict['SeqMap'][prot]:
                if gene in self.list['Spokes']: return gene
            for gene in self.dict['SeqMap'][prot]:
                if gene in self.list['Hubs']: return gene
            return self.dict['SeqMap'][prot][0]
        except: return prot
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
                    try: maxval[mkey] = max(maxval[mkey],string.atoi(occdic[okey][mkey]))
                    except: self.errorLog('%s::%s' % (mkey,occdic[okey][mkey]))
        ### ~ [2] ~ Remake occdic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for okey in rje.sortKeys(occdic):
            if okey == 'Headers': continue
            odata = occdic.pop(okey)
            newkey = []
            for mkey in mainkeys:
                if mkey in maxval:
                    try: newkey.append(rje.preZero(string.atoi(odata[mkey]),maxval[mkey]))
                    except: newkey.append(odata[mkey])
                else: newkey.append(odata[mkey])
            occdic[string.join(newkey,delimit)] = odata
        return occdic
#########################################################################################################################
    ### <3> ### JIM Methods                                                                                             #
#########################################################################################################################
    def makeJIM(self):  ### Generate Graphics                                                                       #V2.0
        '''Generate Graphics.'''
        ### ~ [1] ~ Hub Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        data = self.dict['Data']
        mdb = self.obj['DB'].getTable('Main')
        if self.opt['Hubs']:
            for hub in self.list['Hubs']:
                for ppitype in self.list['RunTypes']:
                    if data[hub]['type'] == 'domain': suffix = ppitype + 'dom'
                    else: suffix = ppitype
                    self.slimJimHub(hub,suffix)
                    dataset = '%s%s' % (hub,suffix)
                    for sdata in mdb.entryList(mdb.index('Dataset')[dataset]):
                        if self.list['SlimList'] and sdata['Pattern'] not in self.list['SlimList']: continue
                        self.slimJimHubMotif(sdata)
        if self.opt['Spokes']: self.slimJimSpokes()
        if self.opt['Motifs']:
            for pattern in rje.sortKeys(self.dict['SlimHubs']):
                if self.list['SlimList'] and pattern not in self.list['SlimList']: continue
                self.slimJimMotif(pattern)        
#########################################################################################################################
    def slimJimHub(self,hub,ppitype):   ### Peforms the SLiMJim visualisation for a given Hub dataset               #V2.0
        '''Peforms the SLiMJim visualisation for a given Hub dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['JIMDir'])
            dataset = '%s.%s' % (hub,ppitype)
            basefile = rje.makePath('%shub/%s/%s.interactome' % (self.info['JIMDir'],hub,hub),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            disfile = '%s%s.dis.tdt' % (self.info['ResDir'],dataset)        # Make tree from this
            ## ~ [1a] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s%s.upc' % (self.info['ResDir'],dataset)            # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = string.split(uline)[3:]
                upc.sort()
                uplist.append(string.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = string.split(upc)
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
                    dis = (string.atof(disdat[obj1][obj2]) + string.atof(disdat[obj2][obj1])) / 2.0
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
            names = []
            for prot in reorder: names.append(self.gene(prot))
            if hub in names: names.remove(hub)    # Will interact with everybody!
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for p1 in names:
                datadict = {}
                for p2 in names:
                    try: ppi = p2 in self.dict['PPI'][p1]  #False
                    except: ppi = False
                    if ppi: datadict[p2] = 1
                    else: datadict[p2] = 0
                rje.delimitedFileOutput(self,pfile,names,datadict=datadict)
            datadict = {}
            for p1 in names: datadict[p1] = upid[p1]
            rje.delimitedFileOutput(self,pfile,names,datadict=datadict)

            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "interactome" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%srje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt','.nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimHubMotif(self,sdata):    ### Peforms the SLiMJim visualisation for a given hub-motif data dictionary #V2.0
        '''Peforms the SLiMJim visualisation for a given hub-motif data dictionary.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Still needs standardising use of hub, spoke, dataset, spokeseq, pattern, slim. #!#
            if sdata['Pattern'] == '-': return
            rje.mkDir(self,self.info['JIMDir'])
            dataset = sdata['Dataset']
            gene = sdata['Hub'] #string.replace(sdata['Dataset'],self.info['Suffix'],'')
            slim = sdata['Pattern']
            rank = string.atoi(sdata['Rank'])
            basefile = rje.makePath('%shub/%s/%s.%s.%s' % (self.info['JIMDir'],gene,gene,rank,rje_slim.slimFromPattern(slim)),wholepath=True)
            rje.mkDir(self,basefile)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            disfile = '%s%s.dis.tdt' % (self.info['ResDir'],dataset)        # Make tree from this
            occfile = '%s%s.occ.csv' % (self.info['ResDir'],dataset)        # Get hit data from this
            seqfile = '%s%s.slimdb' % (self.info['ResDir'],dataset)         # Get extra seq data from this
            motifaln = '%s%s.motifaln.fas' % (self.info['ResDir'],dataset)  # Or this?!
            ## ~ [1a] ~ Load UPC groupings for sequence ordering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upcfile = '%s%s.upc' % (self.info['ResDir'],dataset)            # Get UPC from this
            uplist = []
            for uline in open(upcfile,'r').readlines()[2:]:
                upc = string.split(uline)[3:]
                upc.sort()
                uplist.append(string.join(upc))
            uplist.sort()
            reorder = []
            upx = 1             # UP identifier counter
            upid = {}           # Dictionary of prot:UP_ID
            for upc in uplist:
                uprots = string.split(upc)
                reorder += uprots     #!# Then reorder sequences! #!#
                if len(uprots) > 1:
                    for u in uprots: upid[self.gene(u)] = upx
                    upx += 1
                else: upid[self.gene(uprots[0])] = 0

            ### ~ [2] ~ Load and tidy alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = ['seqin=%s' % motifaln,'autoload=T','seqnr=F','accnr=F','replacechar=F']
            mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            i = rank
            while (string.replace(string.split(mseq.seq[0].info['Sequence'],'-XXXXXXXXXX-')[i-1][10:],'-','') != slim): i += 1
            for seq in mseq.seq[0:]:
                seq.info['Sequence'] = string.split(seq.info['Sequence'],'-XXXXXXXXXX-')[rank-1]
                if not string.replace(seq.info['Sequence'],'-',''): mseq.seq.remove(seq)
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
                if msequence[r] != '[': newm.append(msequence[r]); r += 1
                else:
                    amb = rje.matchExp('^\[([A-Z]+)\]',msequence[r:])[0]
                    newm.append(amb)
                    r += len(amb)+2
            while len(newm) < len(msequence): newm.append('-')
            mseq.seq[0].info['Sequence'] = newm
            ## ~ [2c] ~ Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.aln.tdt' % basefile): os.unlink('%s.aln.tdt' % basefile)
            mseq.saveR(seqfile='%s.aln.tdt' % basefile,name='Gene')

            ### ~ [3] ~ Load full sequences and make profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            scmd = ['seqin=%s' % seqfile,'autoload=T','replacechar=F']
            dseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            seqdic = dseq.seqNameDic()
            occdic = self.loadOccData(occfile,['Dataset','Rank','Pattern','Seq','Start_Pos','End_Pos'])
            occdic.pop('Headers')
            #!# Add UPC weighting of profile? #!#
            ## ~ [3a] ~ Make a list of profile sequences +/- 5 of motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            proseq = []
            for okey in rje.sortKeys(occdic):
                odata = occdic.pop(okey)
                if odata['Pattern'] != slim or sdata['Rank'] != odata['Rank']: continue
                seq = seqdic[odata['Seq']]
                oseq = seq.info['Sequence'][max(0,string.atoi(odata['Start_Pos'])-6):string.atoi(odata['End_Pos'])+5]
                c = max(0,6 - string.atoi(odata['Start_Pos']))
                oseq = '-' * c + oseq
                n = max(0,string.atoi(odata['End_Pos'])+5-seq.aaLen())
                oseq = oseq + '-' * n
                proseq.append(oseq)
            ### ~ [3b] ~ Convert to profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            poslist = ['-5','-4','-3','-2','-1']
            for pos in rje_slim.prestoFromCode(rje_slim.slimFromPattern(slim)):
                if len(pos) > 1: poslist.append('[%s]' % pos)
                elif pos == 'X': poslist.append('x')
                else: poslist.append(pos)
            poslist += ['1','2','3','4','5']
            ignorepos = ['-5','-4','-3','-2','-1','x','1','2','3','4','5']
            outfile = basefile + '.profile.tdt'
            if os.path.exists(outfile): os.unlink(outfile)
            ### ~ [3c] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = ['Pos'] + rje_seq.alph_protx[:-1]
            rje.delimitedFileOutput(self,outfile,headers,'\t')
            try:
                for r in range(len(poslist)):
                    datadict = {}
                    for a in headers[1:]: datadict[a] = 0.0
                    for pseq in proseq:
                        while len(pseq) < len(poslist) and poslist[r] not in ignorepos and poslist[r].find(pseq[r]) < 0: pseq = pseq[:r] + '-' + pseq[r:]
                        try:
                            a = pseq[r]
                            datadict[a] += 1
                        except: pass    # Not counting Xs
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
                    dis = (string.atof(disdat[obj1][obj2]) + string.atof(disdat[obj2][obj1])) / 2.0
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
            names = []
            for prot in reduced: names.append(self.gene(prot))
            pfile = '%s.ppi.tdt' % basefile
            if os.path.exists(pfile): os.unlink(pfile)
            rje.delimitedFileOutput(self,pfile,names)
            for hub in names:
                datadict = {}
                for spoke in names:
                    try: ppi = spoke in self.dict['PPI'][hub]  #False
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
            rslimjim = '%srje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','aln.tdt','profile.tdt','r.tmp.txt','nsf']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimMotif(self,pattern):     ### Performs SLiMJIM visualisation for a given motif                        #V2.0
        '''Performs SLiMJIM visualisation for a given motif.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if pattern == '-': return
            if not pattern in self.dict['SlimSpokes']: return self.printLog('#SLIM','Pattern "%s" not in SLiMSPokes' % pattern)
            data = self.dict['Data']
            slim = rje_slim.slimFromPattern(pattern)
            rje.mkDir(self,self.info['JIMDir'])
            basefile = rje.makePath('%smotif/%s/%s' % (self.info['JIMDir'],slim,slim),wholepath=True)
            if os.path.exists('%s.png' % basefile) and not self.opt['Force']: return
            rje.mkDir(self,basefile)
            ### ~ [2] ~ Generate nested network PPI file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hubs = rje.sortKeys(self.dict['SlimHubs'][pattern]['ALL'])
            spokes = []
            for sseq in self.dict['SlimSpokes'][pattern]['ALL']: spokes.append(self.gene(sseq))
            spokes.sort()
            nfile = '%s.nested.tdt' % basefile
            if os.path.exists(nfile): os.unlink(nfile)
            headers = ['Spoke'] + hubs
            rje.delimitedFileOutput(self,nfile,headers)
            for spoke in spokes:
                datadict = {'Spoke':spoke}
                for hub in hubs:
                    if spoke in data[hub]['ppi'] or hub in data[spoke]['ppi']: datadict[hub] = 1
                    else: datadict[hub] = 0
                rje.delimitedFileOutput(self,nfile,headers,datadict=datadict)
            ### ~ [3] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "nested" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%srje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
            ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                for ext in ['nested.tdt','r.tmp.txt']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimSpokeMapping(self): ### Generates SLiMJIM composite mapping files                                    #V2.0
        '''Generates SLiMJIM composite mapping files.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapdir = self.info['MapDir']
            mx = len(glob.glob('%s*sjmap.fas' % mapdir))
            if os.path.exists(mapdir) and mx and not self.opt['Force']:
                return self.printLog('#SFMAP','%s %s*sjmap.fas found. Skipping slimJimSpokeMapping()' % (rje.integerString(mx),mapdir))
            rje.mkDir(self,mapdir)
            spokeseqlist = {}       # Dictionary of {spokeGene:[Sequence objects]}
            ### ~ [2] ~ Read in Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in self.list['Hubs']:
                for ppitype in self.list['RunTypes']:
                    dataset = '%s%s' % (hub,ppitype)
                    mapaln =  '%s%s.mapping.fas' % (self.info['ResDir'],dataset)        
                    if not os.path.exists(mapaln): continue
                    scmd = ['seqin=%s' % mapaln,'autoload=T','seqnr=F','accnr=F','replacechar=F']
                    mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
                    while mseq.seq:
                        ## ~ [2a] ~ Read in all sequences for one spoke ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        pseq = [mseq.seq.pop(0)]  # Pseq = list of sequences for this protein
                        while mseq.seq:
                            if mseq.seq[0].info['Name'].find('Motifs') > 0 and string.split(mseq.seq[0].info['Name'])[1] == 'Motifs': break   # Next protein
                            pseq.append(mseq.seq.pop(0))
                        ## ~ [2b] ~ Update relevant sequence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        ensloci = pseq[2].shortName()
                        if ensloci in spokeseqlist:         # Update Motifs
                            slims = spokeseqlist[ensloci][0]
                            if pseq[0].aaLen():
                                ppi = pseq[0].info['Sequence']
                                for r in range(len(ppi)):
                                    if ppi[r] != '-': slims.info['Sequence'] = rje.strSub(slims.info['Sequence'],r,r,ppi[r])
                        else: spokeseqlist[ensloci] = pseq[0:]
            ### ~ [3] ~ Save sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ensloci in rje.sortKeys(spokeseqlist):
                mseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                mseq.info['Name'] = '%s%s.sjmap.fas' % (mapdir,ensloci)
                mseq.seq = spokeseqlist[ensloci]
                mseq.saveFasta()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimJimSpoke(self,spoke):   ### Generate SpokeAln PNGs for spoke using SLiMJIM composite mapping file       #V2.0
        '''Generate SpokeAln PNGs for spoke using SLiMJIM composite mapping file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapdir = self.info['MapDir']
            self.printLog('#SPOKE','%s: %d of %d' % (spoke,self.list['Spokes'].index(spoke),len(self.list['Spokes'])),log=False)
            try: ensloci = self.dict['EnsLoci'][spoke]
            except: return self.printLog('#ENSLOCI','No EnsLoci for spoke "%s"' % spoke)
            if ensloci == '-': return self.printLog('#ENSLOCI','No EnsLoci for spoke "%s"' % spoke)
            jdir = rje.makePath('%sspoke/%s/' % (self.info['JIMDir'],spoke))
            rje.mkDir(self,jdir)
            basefile = '%s%s' % (jdir,ensloci)
            if os.path.exists('%s.00001-00110.png' % basefile) and not self.opt['Force']: return
            mapfile = '%s%s.sjmap.fas' % (mapdir,ensloci)
            ### ~ [2] ~ Read in Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ### ~ [3] ~ Make SLiMJIM visualisations for spoke  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Rename sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists(mapfile):
                scmd = ['seqin=%s' % mapfile,'autoload=T','seqnr=F','accnr=F','replacechar=F']
                mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
                pseq = mseq.seq[0:]
                pseq[0].info['R'] = 'PPI SLiMs'
                pseq[1].info['R'] = 'Masked'
                q = 2   # Index of Query 
            else: 
                try: alnfile = '%s%s.orthaln.fas' % (self.info['AlnDir'],string.split(ensloci,'__')[1])
                except: self.printLog('#SPOKE','No alignments found for %s' % ensloci); return
                if not os.path.exists(alnfile): self.printLog('#SPOKE','No alignments found for %s' % ensloci); return
                scmd = ['seqin=%s' % alnfile,'autoload=T','seqnr=F','accnr=F','replacechar=F']
                mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
                pseq = mseq.seq[0:]
                #!# Need to add masked sequence at some point #!#
                q = 0   # Index of Query
            pseq[q].info['R'] = spoke
            for seq in pseq[(q+1):]: seq.info['R'] = rje.getFromDict(rje_ensembl.enscom,seq.info['SpecCode'])
            ## ~ [3b] ~ Setup new SeqList, strip Query gaps, calculate RelCons ~~~~~~~~~~~~~~~~ ##
            seqfile = '%s.aln.tdt' % basefile
            if os.path.exists(seqfile): os.unlink(seqfile)
            rseq = rje_seq.SeqList(self.log,['minregion=3']+self.cmd_list+scmd+['autoload=F'])
            rseq.seq = pseq
            rseq.obj['QuerySeq'] = pseq[q]
            rseq.tidyQueryGaps()
            rseq.saveR(rseq.seq,seqfile,name='R')
            #self.deBug(rseq.obj['QuerySeq'].info['Sequence'])
            rseq.seq = pseq[q:]
            relfile = '%s.rel.tdt' % basefile
            if os.path.exists(relfile): os.unlink(relfile)
            rseq.relCons(relfile)
            ## ~ [3c] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rcmd = '%s --no-restore --no-save --args "spokealn" "%s"' % (self.info['RPath'],basefile)
            rslimjim = '%srje.r' % self.info['Path']
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            problems = self.rCall(rcmd,basefile)
        except: self.errorLog('SLiMJIM visualisation error for spoke "%s"' % spoke)
#########################################################################################################################
    def slimJimSpokes(self):    ### Generate SpokeAln PNGs for all spokes                                           #V2.0
        '''Generate SpokeAln PNGs for all spokes.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map = self.obj['GeneMap'].dict['Data']
            spokeseqlist = {}       # Dictionary of {spokeGene:[Sequence objects]}
            ### ~ [2] ~ Read in Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in self.list['Hubs']:
                for ppitype in self.list['RunTypes']:
                    dataset = '%s%s' % (hub,ppitype)
                    mapaln =  '%s%s.mapping.fas' % (self.info['ResDir'],dataset)        
                    if not os.path.exists(mapaln): continue
                    scmd = ['seqin=%s' % mapaln,'autoload=T','seqnr=F','accnr=F','replacechar=F']
                    mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
                    while mseq.seq:
                        ## ~ [2a] ~ Read in all sequences for one spoke ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        pseq = [mseq.seq.pop(0)]  # Pseq = list of sequences for this protein
                        while mseq.seq:
                            if mseq.seq[0].info['Name'].find('Motifs') > 0 and string.split(mseq.seq[0].info['Name'])[1] == 'Motifs': break   # Next protein
                            pseq.append(mseq.seq.pop(0))
                        ## ~ [2b] ~ Update relevant sequence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        ensloci = pseq[2].shortName()
                        if ensloci in spokeseqlist:         # Update Motifs
                            slims = spokeseqlist[ensloci][0]
                            if pseq[0].aaLen():
                                ppi = pseq[0].info['Sequence']
                                for r in range(len(ppi)):
                                    if ppi[r] != '-': slims.info['Sequence'] = rje.strSub(slims.info['Sequence'],r,r,ppi[r])
                        else: spokeseqlist[ensloci] = pseq[0:]
            ### ~ [3] ~ Make SLiMJIM visualisations for spoke  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0  # Number of errors
            for spoke in self.list['Spokes']:
                try:
                    self.printLog('#SPOKE','%s: %d of %d' % (spoke,self.list['Spokes'].index(spoke),len(self.list['Spokes'])),log=False)
                    try: ensloci = self.dict['EnsLoci'][spoke]
                    except: continue
                    if ensloci == '-': continue
                    jdir = rje.makePath('%sspoke/%s/' % (self.info['JIMDir'],spoke))
                    rje.mkDir(self,jdir)
                    basefile = '%s%s' % (jdir,ensloci)
                    if os.path.exists('%s.00001-00110.png' % basefile) and not self.opt['Force']: continue
                    ## ~ [3a] ~ Rename sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if ensloci in spokeseqlist and spokeseqlist[ensloci]:
                        pseq = spokeseqlist[ensloci][0:]
                        pseq[0].info['R'] = 'PPI SLiMs'
                        pseq[1].info['R'] = 'Masked'
                        q = 2   # Index of Query 
                    else: 
                        try: alnfile = '%s%s.orthaln.fas' % (self.info['AlnDir'],string.split(ensloci,'__')[1])
                        except: self.printLog('#SPOKE','No alignments found for %s' % ensloci); continue
                        if not os.path.exists(alnfile): self.printLog('#SPOKE','No alignments found for %s' % ensloci); continue
                        scmd = ['seqin=%s' % alnfile,'autoload=T','seqnr=F','accnr=F','replacechar=F']
                        mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
                        pseq = mseq.seq[0:]
                        #!# Need to add masked sequence at some point #!#
                        q = 0   # Index of Query
                    pseq[q].info['R'] = spoke
                    for seq in pseq[(q+1):]: seq.info['R'] = rje.getFromDict(rje_ensembl.enscom,seq.info['SpecCode'])
                    ## ~ [3b] ~ Setup new SeqList, strip Query gaps, calculate RelCons ~~~~~~~~~~~~~~~~ ##
                    seqfile = '%s.aln.tdt' % basefile
                    if os.path.exists(seqfile): os.unlink(seqfile)
                    rseq = rje_seq.SeqList(self.log,['minregion=3']+self.cmd_list+scmd+['autoload=F'])
                    rseq.seq = pseq
                    rseq.obj['QuerySeq'] = pseq[q]
                    rseq.tidyQueryGaps()
                    rseq.saveR(rseq.seq,seqfile,name='R')
                    #self.deBug(rseq.obj['QuerySeq'].info['Sequence'])
                    rseq.seq = pseq[q:]
                    relfile = '%s.rel.tdt' % basefile
                    if os.path.exists(relfile): os.unlink(relfile)
                    rseq.relCons(relfile)
                    ## ~ [3c] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    rcmd = '%s --no-restore --no-save --args "spokealn" "%s"' % (self.info['RPath'],basefile)
                    rslimjim = '%srje.r' % self.info['Path']
                    rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
                    problems = self.rCall(rcmd,basefile)
                    if os.path.exists('%s.00001-00110.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                        for ext in ['r.tmp.txt']:   # ['aln.tdt','rel.tdt','r.tmp.txt']:
                            if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
                except: self.errorLog('SLiMJIM visualisation error for spoke "%s"' % spoke); ex += 1
            self.printLog('#SPOKE','Generation of Spoke SLiMJIMs complete. %d Problems.' % ex)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### SLiMJIM HTML generation Methods                                                                         #
#########################################################################################################################
    def makeHTML(self): ### Makes a series of linked HTML pages from read in data                                   #V2.0
        '''Makes a series of linked HTML pages from read in data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Check directories etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hdir = self.info['HTMLDir']
            rje.mkDir(self,hdir)
            css = self.info['HTMLCSS']  #!# Check exists? #!#
            ### ~ [2] ~ Make Front Pages and Hub/Spoke/Motif link lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Mains']: self.mainHTML()
            ### ~ [3] ~ Make Hub pages, including spokes and SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Spokes']:
                (tx,tt) = (0.0,len(self.list['Spokes']))
                for spoke in self.list['Spokes']:
                    self.progLog('\r#SPOKE','Making Spoke HTML: %.2f%%' % (tx/tt)); tx += 100.0
                    self.spokeHTML(spoke)
                self.printLog('\r#SPOKE','Made HTML for %s Spokes.' % rje.integerString(tt))
            ### ~ [4] ~ Make Hub pages, including spokes and SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Hubs']:
                (tx,tt) = (0.0,len(self.dict['Data']))
                for hub in rje.sortKeys(self.dict['Data']):
                    self.progLog('\r#HUB','Making Hub HTML: %.2f%%' % (tx/tt)); tx += 100.0
                    try: self.hubHTML(hub)
                    except: break
                self.printLog('\r#HUB','Made HTML for %s Hubs.' % rje.integerString(tt))
            ### ~ [5] ~ Make SLiM pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Motifs']:
                odb = self.obj['DB'].getTable('Occ')
                odb.index('Pattern')
                (px,ptot) = (0.0,len(odb.dict['Index']['Pattern']))
                for pattern in rje.sortKeys(odb.dict['Index']['Pattern']):
                    self.progLog('\r#SLIM','Generating Motif HTML: %.2f%%' % (px/ptot)); px += 100.0
                    self.motifHTML(pattern)
                self.printLog('\r#SLIM','Generation of HTML for %s Motifs complete.' % rje.integerString(ptot))
            ### ~ [6] ~ Make linked UniFake pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['UniHTML']: self.uniFakeHTML()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def checkHTML(self,hpage):  ### Checks for existence of complete HTML page
        '''Checks for existence of complete HTML page.'''
        if not os.path.exists(hpage): return False
        html = open(hpage,'r').read()
        if html.find('<HTML>') < 0 or html.find('</HTML>') < 0: return False
        return True
#########################################################################################################################
    def mainHTML(self): ### Makes main HTML front pages                                                             #V2.0
        '''Makes main HTML front pages.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdir = self.info['HTMLDir']
            rje.mkDir(self,hdir)
            rje.mkDir(self,'%shub/' % hdir)
            #rje.mkDir(self,'%sdomain/' % hdir)
            rje.mkDir(self,'%sspoke/' % hdir)
            rje.mkDir(self,'%smotif/' % hdir)
            css = self.info['HTMLCSS']
            map = self.obj['GeneMap']
            db = self.obj['DB']
            odb = db.getTable('Occ')
            data = self.dict['Data']
            ## ~ [1a] ~ Gene/Motif counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            (dx,gx,mx,_x) = ({'_':0},{'_':0},{'_':0},{'_':0})
            for x in alph: (dx[x],gx[x],mx[x],_x[x]) = (0,0,0,0)
            for hub in self.dict['Data']:
                x = hub[0]
                if x not in alph: x = '_'
                if data[hub]['type'] == 'gene': gx[x] += 1
                else: dx[x] += 1
            for pattern in rje.sortKeys(odb.dict['Index']['Pattern']):
                x = pattern[0]
                if x not in alph: x = '_'
                mx[x] += 1
                for x in alph:
                    if pattern.count(x) > 0: _x[x] += 1                    
            
            ### ~ [2] ~ Main index page with frames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + 'index.html' 
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % self.info['Name'])
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('<meta name="robots" content="noindex" />\n')
            HTML.write('<meta name="robots" content="nofollow" />\n')           
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET ROWS=150,*><FRAME SRC="main.head.html" NAME="HEAD_FRAME" NORESIZE><FRAME SRC="main.main.html" NAME="MAIN_FRAME">\n')
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()
            ## ~ [2a] ~ Main page split into links and body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + 'main.main.html' 
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % self.info['Name'])
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET COLS=300,*><FRAME SRC="main.link.html" NAME="LINK_FRAME"><FRAME SRC="hub/A.html" NAME="BODY_FRAME">\n')
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()

            ### ~ [3] ~ Header frame with summary data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + 'main.head.html'
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % self.info['Name'])
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><TABLE WIDTH=100%% BORDER=0>\n<TR><TD>\n')
            HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="index.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (self.info['Name']))
            HTML.write('<TD><A HREF="index.html" TARGET=_top><IMG SRC="resources/%s" ALIGN=TOP BORDER=0></A></TD></TR>\n</TABLE>\n' % self.info['HTMLPNG'])
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [4] ~ Links bar to Genes and Motifs with given letters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + 'main.link.html'
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % self.info['Name'])
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><H2>A-Z Results</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n')
            HTML.write('<TD><FONT COLOR=#014359><B>Gene</B></FONT></TD>\n')
            HTML.write('<TD><FONT COLOR=#014359><B>Domain</B></FONT></TD>\n')
            HTML.write('<TD><FONT COLOR=#014359><B>Motif</B></FONT></TD>\n')
            for x in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if not (gx[x] or mx[x] or _x[x] or dx[x]): continue
                HTML.write('</TR><TR><TD>')
                if gx[x]: HTML.write('<A HREF="hub/%s.html" TARGET=BODY_FRAME>%s.. (%d)</A></TD><TD>' % (x,x,gx[x]))
                else: HTML.write('</TD>%s.. (0)<TD>' % x)
                if dx[x]: HTML.write('<A HREF="hub/dom%s.html" TARGET=BODY_FRAME>%s.. (%d)</A></TD><TD>' % (x,x,dx[x]))
                else: HTML.write('</TD>%s.. (0)<TD>' % x)
                if mx[x]: HTML.write('<A HREF="motif/%s_.html" TARGET=BODY_FRAME>%s.. (%d)</A> | ' % (x,x,mx[x]))
                elif x in 'BJOUXZ': HTML.write('-')
                else: HTML.write('%s.. (0) | ' % (x))
                if _x[x]: HTML.write('<A HREF="motif/_%s.html" TARGET=BODY_FRAME>*%s* (%d)</A></TD>\n' % (x,x,_x[x]))
                elif x in 'BJOUXZ': HTML.write('</TD>\n')
                else: HTML.write('*%s* (0)</TD>\n' % (x))
            x = '_'
            if (gx[x] or mx[x] or dx[x]): 
                HTML.write('</TR><TR><TD>')
                if gx[x]: HTML.write('<A HREF="hub/_.html" TARGET=BODY_FRAME>*.. (%d)</A></TD><TD>' % (gx[x]))
                else: HTML.write('</TD>*.. (0)<TD>')
                if dx[x]: HTML.write('<A HREF="hub/dom_.html" TARGET=BODY_FRAME>*.. (%d)</A></TD><TD>' % (dx[x]))
                else: HTML.write('</TD>%s.. (0)<TD>' % x)
                if mx[x]: HTML.write('<A HREF="motif/_.html" TARGET=BODY_FRAME>*.. (%d)</A></TD>\n' % (mx[x]))
                else: HTML.write('*%s* (0)</TD>\n' % (x))
            HTML.write('</TR></TABLE>\n')
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [5] ~ Main tables of genes & domains with given letters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            (hx,htot) = (0.0,len(self.dict['Data']))
            for x in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_':
                ## ~ [5a] ~ Setup gene page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gpage = hdir + 'hub/%s.html' % x
                GHTML = open(gpage,'w')
                if x != '_': GHTML.write('<HTML><HEAD><TITLE>%s Genes</TITLE>\n' % x)
                else: GHTML.write('<HTML><HEAD><TITLE>* Genes</TITLE>\n')
                GHTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                GHTML.write('</HEAD>\n')
                GHTML.write('<BODY><H2>%s Genes</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n' % x)
                GHTML.write('<TD><FONT COLOR=#014359><B>Gene</B></FONT></TD>\n')
                for ppitype in self.list['RunTypes']:
                    GHTML.write('<TD><FONT COLOR=#014359><B>%s</B></FONT></TD>\n' % ppitype)
                GHTML.write('</TR>\n')
                ## ~ [5b] ~ Setup domain page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                dpage = hdir + 'hub/dom%s.html' % x
                DHTML = open(gpage,'w')
                if x != '_': DHTML.write('<HTML><HEAD><TITLE>%s Domains</TITLE>\n' % x)
                else: DHTML.write('<HTML><HEAD><TITLE>* Domains</TITLE>\n')
                DHTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                DHTML.write('</HEAD>\n')
                DHTML.write('<BODY><H2>%s Domains</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n' % x)
                DHTML.write('<TD><FONT COLOR=#014359><B>Domain</B></FONT></TD>\n')
                for ppitype in self.list['RunTypes']:
                    DHTML.write('<TD><FONT COLOR=#014359><B>%s</B></FONT></TD>\n' % ppitype)
                DHTML.write('</TR>\n')
                ## ~ [5c] ~ Output hubs and interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for hub in rje.sortKeys(data):
                    #x#if hub not in self.list['Hubs'] and hub not in self.list['Spokes']: continue
                    if data[hub]['type'] == 'gene': HTML = GHTML
                    else: HTML = DHTML
                    if (x != '_' and hub[:1] != x) or (x == '_' and 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.find(hub[:1]) >=0): continue
                    self.progLog('\r#MAIN','Generating main HTML A-Z Gene/Domain pages: %.1f%%' % (hx/htot)); hx += 100.0
                    if hub in self.list['Hubs']: 
                        slimx = data[hub]['slimx']
                        HTML.write('<TR><TD WIDTH=10%%><A HREF="%s/%s.html" TARGET=_top>%s</A> [%d' % (hub,hub,hub,slimx))
                    else: HTML.write('<TR><TD WIDTH=10%%>%s [-' % (hub))
                    if data[hub]['type'] == 'domain': HTML.write(']</TD>\n') 
                    elif self.worthyHTMLSpoke(hub): 
                        hseq = self.dict['EnsLoci'][hub]    
                        try: spokex = data[hub]['spokeslimx']
                        except: spokex = 0
                        HTML.write(' / <A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%d</A>]</TD>\n' % (hub,hub,spokex))
                    else: HTML.write(' / -]</TD>\n')
                    for ppitype in self.list['RunTypes']:
                        if data[hub]['type'] == 'domain': dataset = '%s.%sdom' % (hub,ppitype)
                        else: dataset = '%s.%s' % (hub,ppitype)
                        occseq = []
                        if dataset in odb.index('Dataset'):
                            for entry in odb.entryList(odb.index('Dataset')[dataset]): occseq.append(entry['Seq'])
                        ilist = []
                        self.deBug(data[hub][ppitype])
                        for spoke in rje.sortKeys(data[hub][ppitype]):
                            try:
                                sseq = self.dict['EnsLoci'][spoke]
                                if sseq in occseq: ilist.append('<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A> [<A HREF="%s/%s.%s.html" TARGET=BODY_FRAME>%d</A>]' % (spoke,spoke,spoke,hub,hub,sseq,occseq.count(sseq)))
                                elif self.worthyHTMLSpoke(spoke): ilist.append('<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A> [0]' % (spoke,spoke,spoke))
                                else: ilist.append('%s [0]' % (spoke))
                            except: ilist.append('%s [-]' % (spoke))
                        HTML.write('<TD>%s;</TD>\n' % string.join(ilist,'; '))
                    HTML.write('</TR>\n')
                for HTML in [GHTML, DHTML]:
                    HTML.write('</TABLE>\n')
                    HTML.write('</BODY></HTML>')
                    HTML.close()
            self.printLog('\r#MAIN','Generating main HTML A-Z Gene/Domain pages complete.')
                
            ### ~ [6] ~ Main tables of Motifs with given letters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            (mx,mtot) = (0.0,41)
            ## ~ [6a] ~ Starting letters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            for x in 'ACDEFGHIKLMNPQRSTVWY_':
                self.progLog('\r#MAIN','Generating main HTML A-Z Motif pages: %.1f%%' % (mx/mtot)); mx += 100.0
                if x == '_': hpage = hdir + 'motif/_.html'
                else: hpage = hdir + 'motif/%s_.html' % x
                HTML = open(hpage,'w')
                if x != '_': HTML.write('<HTML><HEAD><TITLE>%s.. Motifs</TITLE>\n' % x)
                else: HTML.write('<HTML><HEAD><TITLE>* Motifs</TITLE>\n')
                HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                HTML.write('</HEAD>\n')
                if x != '_': HTML.write('<BODY><H2>%s.. Motifs</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n' % x)
                else: HTML.write('<BODY><H2>* Motifs</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n')
                HTML.write('<TD WIDTH=15%%><FONT COLOR=#014359><B>Motif</B></FONT></TD>\n')
                for ppitype in self.list['RunTypes']:
                    HTML.write('<TD><FONT COLOR=#014359><B>%s</B></FONT></TD>\n' % ppitype)
                HTML.write('<TD WIDTH=30%%><FONT COLOR=#014359><B>Spokes</B></FONT></TD>\n')
                HTML.write('</TR>\n')
                for pattern in rje.sortKeys(odb.index('Pattern')):
                    if (x != '_' and pattern[:1] != x) or (x == '_' and 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.find(pattern[:1]) >=0): continue
                    slim = rje_slim.slimFromPattern(pattern)
                    HTML.write('<TR><TD WIDTH=10%%><A HREF="../motif/%s/%s.html" TARGET=_top>%s</A></TD>' % (slim,slim,pattern))
                    # Hubs (sig)
                    for ppitype in self.list['RunTypes']:
                        ilist = []
                        if ppitype in self.dict['SlimHubs'][pattern]:
                            for hub in rje.sortKeys(self.dict['SlimHubs'][pattern][ppitype]):
                                try:
                                    if data[hub]['type'] == 'gene': ilist.append('<A HREF="../hub/%s/%s.html" TARGET=_top>%s</A> (%s)' % (hub,hub,hub,self.dict['SlimHubs'][pattern][ppitype][hub]))
                                    else: ilist.append('<A HREF="../domain/%s/%s.html" TARGET=_top>%s</A> (%s)' % (hub,hub,hub,self.dict['SlimHubs'][pattern][ppitype][hub]))
                                except: self.errorLog('%s missing from data (%s)' % (hub,pattern))
                            ilist.sort()
                            HTML.write('<TD>%s;</TD>\n' % string.join(ilist,'; '))
                        else: HTML.write('<TD>-</TD>\n')
                        # Spokes
                    if pattern in self.dict['SlimSpokes']:
                        ilist = []
                        for spokeseq in rje.sortKeys(self.dict['SlimSpokes'][pattern]):
                            spoke = self.gene(spokeseq)
                            if spoke in self.list['Hubs']: ilist.append('<A HREF="../hub/%s/%s.html" TARGET=_top>%s</A> (<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A>)' % (spoke,spoke,spoke,spoke,spoke,len(self.dict['SlimSpokes'][pattern][spokeseq])))
                            else: ilist.append('%s (<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A>)' % (spoke,spoke,spoke,len(self.dict['SlimSpokes'][pattern][spokeseq])))
                        ilist.sort()
                        HTML.write('<TD WIDTH=30%%>%s;</TD>\n' % string.join(ilist,'; '))
                    else: HTML.write('<TD WIDTH=50%%><I>Unpack Occurrence Tables</I></TD></TR>\n')
                    HTML.write('</TR>\n')
                HTML.write('</TABLE>\n')
                HTML.write('</BODY></HTML>')
                HTML.close()
            ## ~ [6b] ~ Motifs containing x ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            for x in 'ACDEFGHIKLMNPQRSTVWY':
                self.progLog('\r#MAIN','Generating main HTML A-Z Motif pages: %.1f%%' % (mx/mtot)); mx += 100.0
                hpage = hdir + 'motif/_%s.html' % x
                HTML = open(hpage,'w')
                HTML.write('<HTML><HEAD><TITLE>Motifs containing %s</TITLE>\n' % x)
                HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                HTML.write('</HEAD>\n')
                HTML.write('<BODY><H2>Motifs containing %s</H2>\n<TABLE WIDTH=100%% BORDER=2>\n<TR>\n' % x)
                HTML.write('<TD WIDTH=15%%><FONT COLOR=#014359><B>Motif</B></FONT></TD>\n')
                for ppitype in self.list['RunTypes']:
                    HTML.write('<TD><FONT COLOR=#014359><B>%s</B></FONT></TD>\n' % ppitype)
                HTML.write('<TD WIDTH=30%%><FONT COLOR=#014359><B>Spokes</B></FONT></TD>\n')
                HTML.write('</TR>\n')
                for pattern in rje.sortKeys(odb.index('Pattern')):
                    if pattern.find(x) < 0: continue
                    slim = rje_slim.slimFromPattern(pattern)
                    HTML.write('<TR><TD WIDTH=10%%><A HREF="../motif/%s/%s.html" TARGET=_top>%s</A></TD>' % (slim,slim,pattern))
                    # Hubs (sig)
                    for ppitype in self.list['RunTypes']:
                        if ppitype in self.dict['SlimHubs'][pattern]:
                            if ppitype == 'ALL': continue
                            ilist = []
                            for hub in rje.sortKeys(self.dict['SlimHubs'][pattern][ppitype]):
                                if data[hub]['type'] == 'gene': ilist.append('<A HREF="../hub/%s/%s.html" TARGET=_top>%s</A> (%s)' % (hub,hub,hub,self.dict['SlimHubs'][pattern][ppitype][hub]))
                                else: ilist.append('<A HREF="../domain/%s/%s.html" TARGET=_top>%s</A> (%s)' % (hub,hub,hub,self.dict['SlimHubs'][pattern][ppitype][hub]))
                            ilist.sort()
                            HTML.write('<TD>%s;</TD>\n' % string.join(ilist,'; '))
                        else: HTML.write('<TD>-</TD>\n')
                    # Spokes
                    if pattern in self.dict['SlimSpokes']:
                        ilist = []
                        for spokeseq in rje.sortKeys(self.dict['SlimSpokes'][pattern]):
                            spoke = self.gene(spokeseq)
                            if spoke in self.list['Hubs']: ilist.append('<A HREF="../hub/%s/%s.html" TARGET=_top>%s</A> (<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A>)' % (spoke,spoke,spoke,spoke,spoke,len(self.dict['SlimSpokes'][pattern][spokeseq])))
                            else: ilist.append('%s (<A HREF="../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A>)' % (spoke,spoke,spoke,len(self.dict['SlimSpokes'][pattern][spokeseq])))
                        ilist.sort()
                        HTML.write('<TD WIDTH=50%%>%s;</TD></TR>\n' % string.join(ilist,'; '))
                    else: HTML.write('<TD WIDTH=50%%><I>Unpack Occurrence Tables</I></TD></TR>\n')
                    #!# Eventually make a nested network graphic of hubs and spokes #!#
                HTML.write('</TABLE>\n')
                HTML.write('</BODY></HTML>')
                HTML.close()
            self.printLog('\r#MAIN','Generating main HTML A-Z Motif pages complete.')
        except: self.errorLog(rje_zen.Zen().wisdom())                    
#########################################################################################################################
    def worthyHTMLSpoke(self,spoke,allowhubs=False):    ### Returns True if it "deserves" HTML made, else False
        '''Returns True if it "deserves" HTML made, else False.'''
        if spoke not in self.list['Spokes'] and not (allowhubs and spoke in self.list['Hubs']): return False
        try: ensloci = self.dict['EnsLoci'][spoke]
        except: return False
        if ensloci == '-': return False
        return True
#########################################################################################################################
    def motifGoHTML(self,genelist): ### Returns HTML of GO children for gene list
        '''Returns HTML of GO children for gene list.'''
        try:### ~ [1] ~ Look for gene, then try symbol, then try EnsEMBL mapped to symbol ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['GO'].dict['GO']: return '<I>No GO Data read by SLiMJIM.</I>'
            go = []
            for mygene in genelist:
                gene = self.gene(mygene)                                          # Try mapping to gene in case protein
                if gene in self.dict['GO']: go += self.dict['GO'][gene]          # Mapped gene found in dictionary
                elif not self.obj['GeneMap'].getGeneData(gene): continue                 # No mapping to EnsEMBL via GeneMap
                else:
                    ens = rje.dictValues(self.obj['GeneMap'].getGeneData(gene),'EnsEMBL',valtype='str') # Try to map EnsEMBL
                    if ens and ens in self.dict['GO']: go += self.dict['GO'][ens]  # Found EnsEMBL mapping in GO dictionary
            ### ~ [2] ~ Convert list to HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lengo = 0
            while lengo != len(go):     # Remove parents from list
                lengo = len(go)
                for id in go[0:]:
                    for child in self.obj['GO'].go(id)['child_terms']:
                        if child != id and child in go: go.remove(id); break
            if not go: return '<I>No GO categories associated with %s.</I>' % gene
            glist = []
            for id in go:
                if id[:3] != 'GO:': id = 'GO:%s' % id
                glist.append('<A HREF="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=%s" TARGET=BODY_FRAME>%s</A>' % (id,self.obj['GO'].go(id)['name']))
            return string.join(glist,'; ')
        except: self.errorLog('SLiMJIM.getGO(%s) error' % gene)
        return '<I>SLiMJIM goTHML Error!.</I>'                                                           # Return empty list if failed
#########################################################################################################################
    def goHTML(self,gene):   ### Returns HTML of GO children for a given gene (or protein) using EnsGO and GeneMap
        '''Returns HTML of GO children for a given gene (or protein) using EnsGO and GeneMap.'''
        try:### ~ [1] ~ Look for gene, then try symbol, then try EnsEMBL mapped to symbol ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['GO'].dict['GO']: return '<I>No GO Data read by SLiMJIM.</I>'
            gene = self.gene(gene)                                          # Try mapping to gene in case protein
            if gene in self.dict['GO']: go = self.dict['GO'][gene]          # Mapped gene found in dictionary
            elif not self.obj['GeneMap'].getGeneData(gene):                 # No mapping to EnsEMBL via GeneMap
                return '<I>No Gene Mapping available for processing GO Data.</I>'    
            else:
                ens = rje.dictValues(self.obj['GeneMap'].getGeneData(gene),'EnsEMBL',valtype='str') # Try to map EnsEMBL
                if ens and ens in self.dict['GO']: go = self.dict['GO'][ens]  # Found EnsEMBL mapping in GO dictionary
                elif ens: return '<I>No GO categories associated with %s.</I>' % ens
                else: return '<I>EnsEMBL mapping failed. Could not extract GO categories.</I>'
            ### ~ [2] ~ Convert list to HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug('%s: %s' % (gene,go))
            lengo = 0
            while lengo != len(go):     # Remove parents from list
                lengo = len(go)
                for id in go[0:]:
                    self.deBug('%s: %s' % (id,self.obj['GO'].go(id)['child_terms']))
                    for child in self.obj['GO'].go(id)['child_terms']:
                        if child != id and child in go: go.remove(id); break
            if not go: return '<I>No GO categories associated with %s.</I>' % gene
            glist = []
            for id in go:
                if id[:3] != 'GO:': id = 'GO:%s' % id
                glist.append('<A HREF="http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=%s" TARGET=BODY_FRAME>%s</A>' % (id,self.obj['GO'].go(id)['name']))
            return string.join(glist,'; ')
        except: self.errorLog('SLiMJIM.getGO(%s) error' % gene)
        return '<I>SLiMJIM goTHML Error!.</I>'                                                           # Return empty list if failed
#########################################################################################################################
    def spokeHTML(self,spoke):  ### Makes HTML for a given spoke gene                                               #V2.0
        '''Makes HTML for a given spoke gene.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Neaten with spokeseq etc. #!#
            if not self.worthyHTMLSpoke(spoke): return
            hdir = self.info['HTMLDir']
            css = self.info['HTMLCSS']
            map = self.obj['GeneMap'].dict['Data']
            try: seq = map[spoke]['EnsLoci']
            except: seq = None
            hdir = rje.makePath('%sspoke/%s/' % (hdir,spoke))
            rje.mkDir(self,hdir)
            hpage = hdir + '%s.html' % spoke
            if not self.opt['Force'] and self.checkHTML(hpage): return
            data = self.dict['Data']
            ## ~ [1a] ~ Occurrence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            odb = self.obj['DB'].getTable('Occ')
            hubocc = {}     # {hub:{type:no. occ}}
            for ppitype in self.list['RunTypes']: hubocc[ppitype] = {}
            for hub in rje.sortKeys(data[spoke]['ppi']):
                for ppitype in self.list['RunTypes']: hubocc[ppitype][hub] = 0
                if hub not in odb.index('Hub'): continue
                hx = 0
                for entry in odb.entryList(odb.index('Hub')[hub]):
                    if entry['Seq'] == seq: hx += 1
                for ppitype in self.list['RunTypes']:
                    if spoke in data[hub][ppitype]: hubocc[ppitype][hub] += hx
            spocc = {}      # {key:{occdata}}
            if seq in odb.index('Seq'):
                for entry in odb.entryList(odb.index('Seq')[seq]):
                    spocc['%s_%s' % (entry['Dataset'],rje.preZero(entry['Rank'],1000))] = entry
            ### ~ [2] ~ HTML setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % spoke)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            ## ~ [2a] ~ Gene summary as with hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<BODY>\n')
            HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="../../hub/%s/%s.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (spoke,spoke,spoke))
            if seq: HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B><A HREF="../../unifake/%s.html" TARGET=BODY_FRAME>%s</A></B></FONT></TD>\n' % (seq,seq))
            else: HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B>-</B></FONT></TD>\n')
            try:
                HTML.write('<TD>\n')
                gdata = map[spoke]
                if 'Symbol' in gdata and gdata['Symbol']: HTML.write('<A HREF="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" TARGET=BODY_FRAME><IMG ALT="GeneCards" SRC="../../resources/genecard.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['Symbol'])
                if 'UniProt' in gdata and gdata['UniProt']: HTML.write('<A HREF="http://www.uniprot.org/uniprot/%s" TARGET=BODY_FRAME><IMG ALT="UniProt" SRC="../../resources/logo_ebi.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['UniProt'])
                if 'EnsEMBL' in gdata and gdata['EnsEMBL']: HTML.write('<A HREF="http://www.ensembl.org/Homo_sapiens/geneview?gene=%s" TARGET=BODY_FRAME><IMG ALT="EnsEMBL" SRC="../../resources/e-bang.gif" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['EnsEMBL'])
                if 'HPRD' in gdata and gdata['HPRD']:  HTML.write('<A HREF="http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1" TARGET=BODY_FRAME><IMG ALT="HPRD" SRC="../../resources/hprd.png" ALIGN=BOTTOM BORDER=0></A>\n' % (gdata['HPRD'],gdata['HPRD']))
                if 'OMIM' in gdata and gdata['OMIM']:  HTML.write('<A HREF="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s" TARGET=BODY_FRAME><IMG ALT="OMIM" SRC="../../resources/omim.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['OMIM'])
            except: pass    #self.errorLog('GeneMap data problem for %s' % spoke)
            HTML.write('</TD></TR>\n')
            try:
                HTML.write('<TD>\n')
                gdata = map[spoke]
                if 'Desc' in gdata: HTML.write('<TR><TD COLSPAN=3><FONT SIZE=4 COLOR=#014359><P>%s</FONT></P></TD></TR>\n' % gdata['Desc'])
            except: pass
            HTML.write('</TD></TR>\n')
            HTML.write('<TR><TD COLSPAN=3><P>%s</P></TD></TR>\n' % self.goHTML(spoke)) 
            HTML.write('</TABLE>\n')
            ## ~ [2b] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<H2>%s Interactors</H2>\n' % spoke)
            for ppitype in self.list['RunTypes']:
                HTML.write('<H3>%s Interactors</H3>\n' % ppitype)
                ilist = []
                for hub in rje.sortKeys(hubocc[ppitype]):
                    if hubocc[ppitype][hub]: ilist.append('<A HREF="../../hub/%s/%s.html" TARGET=_top><U>%s</U></A> [<A HREF="../../hub/%s/%s.%s.html">%d</A>]' % (hub,hub,hub,hub,hub,seq,hubocc[ppitype][hub]))
                    elif hub in self.list['Hubs']: ilist.append('<A HREF="../../hub/%s/%s.html" TARGET=_top>%s</A>' % (hub,hub,hub))
                    else: ilist.append('%s' % hub)
                HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
            ## ~ [2c] ~ Occurrence Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<H2>SLiM Occurrences in %s</H2>\n' % (spoke))
            if spocc:
                self.writeHTMLTable(HTML,self.list['OccHead'],spocc,asdict=True,gene=None,drop=['Seq','Desc'])
                seqlen = string.atoi(spocc.values()[0]['Prot_Len'])
            else:
                HTML.write('<P><I>No predicted SLiMs in %s</I></P>\n' % (spoke))
                seqlen = 0
            ## ~ [2d] ~ Graphics and End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Intersperse with tables at some point? #!#
            HTML.write('<H2>Summary Graphics</H2>\n')
            vis = glob.glob('%s*.png' % (hdir))
            if not vis and seqlen:
                i = 1
                while i < seqlen:
                    vis.append('%s.%s-%s.png' % (seq,rje.preZero(i,99999),rje.preZero(i+109,99999)))
                    if (i+109) >= seqlen: break
                    i += 100
            vis.sort()
            if vis:
                for pfile in vis: 
                    png = os.path.basename(pfile)
                    HTML.write('<A HREF="%s" TARGET=_blank><IMG SRC="%s" ALT="%s" WIDTH=100%%></A><BR>\n' % (png,png,png))
            else: HTML.write('<P><I>No Summary Graphics made</I></P>\n')
            HTML.write('</BODY></HTML>')
            HTML.close()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def motifHTML(self,pattern):    ### Makes HTML for a given SLiM pattern                                         #V2.0
        '''Makes HTML for a given SLiM pattern.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Check and neaten #!#
            hdir = self.info['HTMLDir']
            css = self.info['HTMLCSS']
            if pattern == '-': return
            slim = rje_slim.slimFromPattern(pattern)
            hdir = rje.makePath('%smotif/%s/' % (hdir,slim))
            rje.mkDir(self,hdir)
            if not self.opt['Force'] and self.checkHTML(hdir + '%s.interactome.html' % slim): return
            odb = self.obj['DB'].getTable('Occ')
            spokes = {}; allspokes = []
            hubs = {}; allhubs = []
            allgenes = []
            for ppitype in self.list['RunTypes']: spokes[ppitype] = []; hubs[ppitype] = []
            for entry in odb.entryList(odb.index('Pattern')[pattern]):
                ppitype = string.split(entry['Dataset'],'.')[-1][:3]
                if entry['Hub'] not in hubs[ppitype]: hubs[ppitype].append(entry['Hub'])
                if entry['Seq'] not in spokes[ppitype]: spokes[ppitype].append(self.gene(entry['Seq']))
                if entry['Hub'] not in allgenes: allgenes.append(entry['Hub'])
                if entry['Seq'] not in allgenes: allgenes.append(self.gene(entry['Seq']))
                spokes[ppitype].sort(); hubs[ppitype].sort()
                allspokes += spokes[ppitype]
                allhubs += hubs[ppitype]
            allspokes = rje.sortUnique(allspokes); allhubs = rje.sortUnique(allhubs); 

            ### ~ [2] ~ Make main frame pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Outer frame page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + '%s.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET ROWS=150,*><FRAME SRC="%s.head.html" NAME="HEAD_FRAME" NORESIZE><FRAME SRC="%s.main.html" NAME="MAIN_FRAME">\n' % (slim,slim))
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()
            ## ~ [2b] ~ Inner Frame page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + '%s.main.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET COLS=300,*><FRAME SRC="%s.link.html" NAME="LINK_FRAME"><FRAME SRC="%s.results.html" NAME="BODY_FRAME">\n' % (slim,slim))
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()

            ### ~ [3] ~ Main Header page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.head.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><TABLE WIDTH=100%% BORDER=0>\n<TR><TD>\n')
            HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="%s.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (slim,pattern))
            HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45>%s</FONT></TD></TR>\n' % (slim))
            HTML.write('<TR><TD COLSPAN=3><P>%s</P></TD></TR>\n</TABLE>' % self.motifGoHTML(allgenes))
            HTML.write('</TD><TD><A HREF="../../index.html" TARGET=_top><IMG SRC="../../resources/%s" ALIGN=TOP BORDER=0></A></TD></TR>\n</TABLE>\n' % self.info['HTMLPNG'])
            HTML.write('</BODY></HTML>')
            HTML.close()
            ### ~ [3] ~ Hub/Spoke links page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.link.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><H2>%s Results</H2><UL>\n' % pattern)
            HTML.write('<A HREF="%s.results.html" TARGET=BODY_FRAME><LI>Summary Tables</A>\n' % slim)
            HTML.write('<A HREF="%s.interactome.html" TARGET=BODY_FRAME><LI>Interactome</A>\n' % slim)
            HTML.write('</UL>\n<H2>Hubs</H2>\n')
            for ppitype in self.list['RunTypes']:
                if not hubs[ppitype]: continue
                HTML.write('</UL>\n<H3>%s hubs</H3>\n<OL>\n' % ppitype)
                for hub in hubs[ppitype]:
                    HTML.write('<A HREF="../../hub/%s/%s.html" TARGET=_top><LI>%s</A> (<A HREF="../../hub/%s/%s.%s.html" TARGET=BODY_FRAME>%s</A>)\n' % (hub,hub,hub,hub,hub,slim,self.dict['SlimHubs'][pattern][ppitype][hub]))
                HTML.write('</OL>\n')
            HTML.write('\n<H2>Spokes</H2>\n')
            for ppitype in self.list['RunTypes']:
                if not spokes[ppitype]: continue
                HTML.write('<H3>%s spokes</H3>\n<OL>\n' % ppitype)
                for spoke in spokes[ppitype]:
                    try: HTML.write('<A HREF="../../spoke/%s/%s.html" TARGET=BODY_FRAME><LI>%s</A> [%s]\n' % (spoke,spoke,spoke,len(self.dict['SlimSpokes'][pattern][self.dict['EnsLoci'][spoke]])))
                    except: HTML.write('<A HREF="../../spoke/%s/%s.html" TARGET=BODY_FRAME><LI>%s</A> [?]\n' % (spoke,spoke,spoke))
            HTML.write('</OL>\n</BODY></HTML>')
            HTML.close()

            ### ~ [4] ~ Main Results page = Results table(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.results.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s Results</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n<BODY>')
            ## ~ [4a] ~ Summary Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mdb = self.obj['DB'].getTable('Main')
            HTML.write('<H2>%s SLiMFinder Summary Results</H2>\n' % (pattern))
            for ppitype in self.list['RunTypes']:
                HTML.write('<H3>%s</H3>\n' % (ppitype))
                slimdat = {}
                for entry in mdb.entryList(mdb.index('Pattern')[pattern]):
                    if string.split(entry['Dataset'],'.')[-1][:3] == ppitype: slimdat[entry['Dataset']] = entry
                self.writeHTMLTable(HTML,self.list['SumHead'],slimdat,asdict=True,drop=['Pattern'])    
            HTML.write('<H2>%s SLiM Occurrence Results</H2>\n' % (pattern))
            for ppitype in self.list['RunTypes']:
                HTML.write('<H3>%s</H3>\n' % (ppitype))
                slimdat = {}
                for entry in odb.entryList(odb.index('Pattern')[pattern]):
                    if string.split(entry['Dataset'],'.')[-1][:3] == ppitype: slimdat[odb.makeKey(entry)] = entry
                self.writeHTMLTable(HTML,self.list['SumHead'],slimdat,asdict=True,drop=['Pattern'])    
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [5] ~ Interactome summary page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.interactome.html' % slim
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s Interactome</TITLE>\n' % pattern)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n<BODY>')
            ## ~ [5a] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<H2>%s Hubs</H2>\n' % pattern)
            ilist = []
            for hub in allhubs:
                sigtxt = []
                for ppitype in self.dict['SlimHubs'][pattern]:
                    if ppitype == 'ALL': continue
                    if hub in self.dict['SlimHubs'][pattern][ppitype]:
                        sigtxt.append('%s:%s' % (ppitype,self.dict['SlimHubs'][pattern][ppitype][hub]))
                sigtxt = string.join(sigtxt,'; ')
                ilist.append('<A HREF="../../hub/%s/%s.html" TARGET=_top>%s</A> (<A HREF="../../hub/%s/%s.%s.html" TARGET=BODY_FRAME>%s</A>)' % (hub,hub,hub,hub,hub,slim,sigtxt))
            HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
            HTML.write('<H2>%s Spokes</H2>\n' % pattern)
            ilist = []
            for spoke in allspokes:
                sigtxt = []
                for ppitype in self.dict['SlimSpokes'][pattern]:
                    sigtxt.append('%s:%s' % (ppitype,len(self.dict['SlimSpokes'][pattern][ppitype])))
                sigtxt = string.join(sigtxt,'; ')
                ilist.append('<A HREF="../../spoke/%s/%s.html" TARGET=BODY_FRAME>%s</A> [%s]' % (spoke,spoke,spoke,sigtxt))
            HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
            HTML.write('<H2>%s Interactome</H2>\n' % pattern)
            HTML.write('<A HREF="%s.png" TARGET=_blank"><IMG SRC="%s.png" WIDTH=100%% ALT="%s Interactome graphics"></A>\n' % (slim,slim,slim))
            HTML.write('</BODY></HTML>')
            HTML.close()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def hubHTML(self,hub):  ### Generates hub HTML for dataset                                                      #V2.0
        '''Generates hub HTML for dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = self.dict['Data']
            hdir = self.info['HTMLDir']
            css = self.info['HTMLCSS']
            map = self.obj['GeneMap']
            gene = hub
            if gene in self.dict['EnsLoci'] and self.dict['EnsLoci'][gene] != '-': seq = self.dict['EnsLoci'][gene]
            else: seq = None
            hdir = rje.makePath('%shub/%s/' % (hdir,gene))
            rje.mkDir(self,hdir)
            mdb = self.obj['DB'].getTable('Main')
            odb = self.obj['DB'].getTable('Occ')
            ### ~ [2] ~ Make main frame pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Outer frame page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + '%s.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET ROWS=150,*><FRAME SRC="%s.head.html" NAME="HEAD_FRAME" NORESIZE><FRAME SRC="%s.main.html" NAME="MAIN_FRAME">\n' % (gene,gene))
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()
            ## ~ [2b] ~ Inner Frame page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + '%s.main.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<FRAMESET COLS=300,*><FRAME SRC="%s.link.html" NAME="LINK_FRAME"><FRAME SRC="%s.results.html" NAME="BODY_FRAME">\n' % (gene,gene))
            HTML.write('<NOFRAMES><BODY><H1>You don\'t support frames!</H1></BODY></NOFRAMES></FRAMESET></HTML>')
            HTML.close()
            ### ~ [3] ~ Main Header page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.head.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><TABLE WIDTH=100%% BORDER=0>\n<TR><TD>\n')
            HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="%s.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (gene,gene))
            try: HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B><A HREF="../../unifake/%s.html" TARGET=BODY_FRAME>%s</A></B></TD>\n' % (self.dict['EnsLoci'][gene],self.dict['EnsLoci'][gene]))
            except: HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B>-</B></TD>\n')
            try:
                HTML.write('<TD>\n')
                gdata = self.obj['GeneMap'].dict['Data'][gene]  #1500#
                if 'Symbol' in gdata and gdata['Symbol']: HTML.write('<A HREF="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" TARGET=BODY_FRAME><IMG ALT="GeneCards" SRC="../../resources/genecard.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['Symbol'])
                if 'UniProt' in gdata and gdata['UniProt']: HTML.write('<A HREF="http://www.uniprot.org/uniprot/%s" TARGET=BODY_FRAME><IMG ALT="UniProt" SRC="../../resources/logo_ebi.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['UniProt'])
                if 'EnsEMBL' in gdata and gdata['EnsEMBL']: HTML.write('<A HREF="http://www.ensembl.org/Homo_sapiens/geneview?gene=%s" TARGET=BODY_FRAME><IMG ALT="EnsEMBL" SRC="../../resources/e-bang.gif" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['EnsEMBL'])
                if 'HPRD' in gdata and gdata['HPRD']:  HTML.write('<A HREF="http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1" TARGET=BODY_FRAME><IMG ALT="HPRD" SRC="../../resources/hprd.png" ALIGN=BOTTOM BORDER=0></A>\n' % (gdata['HPRD'],gdata['HPRD']))
                if 'OMIM' in gdata and gdata['OMIM']:  HTML.write('<A HREF="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s" TARGET=BODY_FRAME><IMG ALT="OMIM" SRC="../../resources/omim.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['OMIM'])
            except: self.printLog('#ERR','GeneMap data problem for %s' % gene,screen=False)
            HTML.write('</TD></TR>\n')
            try:
                HTML.write('<TD>\n')
                gdata = self.obj['GeneMap'].dict['Data'][gene]
                if 'Desc' in gdata: HTML.write('<TR><TD COLSPAN=3><FONT SIZE=4 COLOR=#014359><P>%s</FONT></P></TD></TR>\n' % gdata['Desc'])
            except: pass
            HTML.write('</TD></TR>\n')
            HTML.write('<TR><TD COLSPAN=3><P>%s</P></TD></TR>\n</TABLE>' % self.goHTML(gene))
            HTML.write('</TD><TD><A HREF="../../index.html" TARGET=_top><IMG SRC="../../resources/%s" ALIGN=TOP BORDER=0></A></TD></TR>\n</TABLE>\n' % self.info['HTMLPNG'])
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [4] ~ Spoke/SLiM links page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] ~ Setup data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            occx = {}
            for ppitype in self.list['RunTypes']:
                occx[ppitype] = {}
                for p1 in data[gene][ppitype]: occx[ppitype][p1] = 0
            if hub in odb.index('Hub'):
                for entry in odb.entryList(odb.index('Hub')[hub]):
                    ppitype = string.split(entry['Dataset'],'.')[-1][:3]
                    for spoke in self.obj['DB'].getTable('SeqMap').data()[entry['Seq']]['Spoke']:
                        if spoke in occx[ppitype]: occx[ppitype][spoke] += 1
            ## ~ [4b] ~ Make page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpage = hdir + '%s.link.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n')
            HTML.write('<BODY><H2>%s Results</H2><UL>\n' % gene)
            HTML.write('<A HREF="%s.results.html" TARGET=BODY_FRAME><LI>Summary Tables</A>\n' % gene)
            HTML.write('<A HREF="%s.interactome.html" TARGET=BODY_FRAME><LI>Interactome</A>\n' % gene)
            if seq and seq in odb.index('Seq'):
                HTML.write('<A HREF="../../spoke/%s/%s.html" TARGET=BODY_FRAME><LI>Spoke Results</A>\n' % (gene,gene))
            if seq:
                HTML.write('<A HREF="../../unifake/%s.html" TARGET=BODY_FRAME><LI>UniFake</A>\n' % (seq))
            HTML.write('</UL>\n<H2>SLiMS</H2>\n')
            if hub in mdb.index('Hub'):
                ekeys = mdb.index('Hub')[hub]
                for ppitype in self.list['RunTypes']:
                    HTML.write('<H3>%s SLiMS</H3>\n' % ppitype)
                    sdata = {}
                    for e in ekeys:
                        entry = mdb.data()[e]
                        if entry['PPISet'][:3] == ppitype or entry['Rank'] < 1: continue
                        slim = rje_slim.slimFromPattern(entry['Pattern'])
                        try: sdata[rje.preZero(entry['Rank'],len(ekeys))] = '<A HREF="%s.%s.html" TARGET=BODY_FRAME><LI>%s</A> [%s] (%s)\n' % (gene,slim,entry['Pattern'],entry['Occ'],rje_slim.expectString(entry['Sig']))
                        except: sdata[rje.preZero(entry['Rank'],len(ekeys))] = '<A HREF="%s.%s.html" TARGET=BODY_FRAME><LI>%s</A> [%s] (%s)\n' % (gene,slim,entry['Pattern'],entry['Occ'],entry['Sig'])
                    if sdata:
                        HTML.write('<OL>\n')
                        for rkey in rje.sortKeys(sdata): HTML.write(sdata[rkey])
                        HTML.write('</OL>\n\n')
                    else: HTML.write('<UL><LI>No SLiMs returned</UL>\n\n')
            else: HTML.write('<UL><LI>No SLiMFinder runs.</UL>\n\n')
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [5] ~ Main Results page = Results table(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.results.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s Results</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n<BODY>')
            ## ~ [5a] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<H2>%s Interactors</H2>\n' % gene)
            for ppitype in self.list['RunTypes']:
                if not occx[ppitype]: continue
                HTML.write('<H3>%s</H3>\n' % ppitype)
                ilist = []
                for spoke in rje.sortKeys(occx[ppitype]):
                    if occx[ppitype][spoke]:
                        seq = self.dict['EnsLoci'][spoke]
                        ilist.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A> [<A HREF="%s.%s.html">%d</A>]' % (spoke,spoke,spoke,gene,seq,occx[ppitype][spoke]))
                    else: ilist.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A> [-]' % (spoke,spoke,spoke))
                HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
            ## ~ [5b] ~ Summary Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<H2>%s SLiMFinder Summary Results</H2>\n' % (gene))
            self.writeHTMLTable(HTML,self.list['SumHead'],mdb.subset('Hub',hub),asdict=False,gene=gene,drop=['Hub'])   
            if gene in odb.index('Hub'): 
                HTML.write('<H2>%s SLiM Occurrence Results</H2>\n' % (gene))
                self.writeHTMLTable(HTML,self.list['OccHead'],odb.subset('Hub',hub),asdict=True,gene=gene,drop=['Hub'])   
            else: HTML.write('<P><I>No SLiMs returned.</I></P>/n')
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [6] ~ Interactome summary page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpage = hdir + '%s.interactome.html' % gene
            HTML = open(hpage,'w')
            HTML.write('<HTML><HEAD><TITLE>%s Interactome</TITLE>\n' % gene)
            HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
            HTML.write('</HEAD>\n<BODY>')
            ## ~ [6a] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ppitype in self.list['RunTypes']:
                if not occx[ppitype]: continue
                HTML.write('<H2>%s Interactome</H2>\n' % ppitype)
                ilist = []
                for spoke in rje.sortKeys(occx[ppitype]):
                    if occx[ppitype][spoke]:
                        seq = self.dict['EnsLoci'][spoke]
                        ilist.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A> [<A HREF="%s.%s.html">%d</A>]' % (spoke,spoke,spoke,gene,seq,occx[ppitype][spoke]))
                    else: ilist.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A> [-]' % (spoke,spoke,spoke))
                HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
            #!# Add a figure for each ppitype? #!#
            HTML.write('<A HREF="%s.interactome.png" TARGET=_blank"><IMG SRC="%s.interactome.png" WIDTH=100%% ALT="%s Interactome graphics"></A>\n' % (gene,gene,gene))
            HTML.write('</BODY></HTML>')
            HTML.close()

            ### ~ [7] ~ Hub-Spoke Summary pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spoke in data[gene]['ppi']:
                hpage = hdir + '%s.%s.html' % (gene,spoke)
                if not self.opt['Force'] and self.checkHTML(hpage): continue
                HTML = open(hpage,'w')
                HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % gene)
                HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                HTML.write('</HEAD>\n')
                ## ~ [6a] ~ Gene summary as with hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sgene = self.gene(spoke)
                HTML.write('<BODY>\n')
                HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="../%s/%s.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (sgene,sgene,sgene))
                HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B><A HREF="../../unifake/%s.html" TARGET=BODY_FRAME>%s</A></B></TD>\n' % (spoke,spoke))
                try:
                    HTML.write('<TD>\n')
                    gdata = self.obj['GeneMap'].dict['Data'][sgene]
                    if 'Symbol' in gdata and gdata['Symbol']: HTML.write('<A HREF="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" TARGET=BODY_FRAME><IMG ALT="GeneCards" SRC="../../resources/genecard.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['Symbol'])
                    if 'UniProt' in gdata and gdata['UniProt']: HTML.write('<A HREF="http://www.uniprot.org/uniprot/%s" TARGET=BODY_FRAME><IMG ALT="UniProt" SRC="../../resources/logo_ebi.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['UniProt'])
                    if 'EnsEMBL' in gdata and gdata['EnsEMBL']: HTML.write('<A HREF="http://www.ensembl.org/Homo_sapiens/geneview?gene=%s" TARGET=BODY_FRAME><IMG ALT="EnsEMBL" SRC="../../resources/e-bang.gif" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['EnsEMBL'])
                    if 'HPRD' in gdata and gdata['HPRD']:  HTML.write('<A HREF="http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1" TARGET=BODY_FRAME><IMG ALT="HPRD" SRC="../../resources/hprd.png" ALIGN=BOTTOM BORDER=0></A>\n' % (gdata['HPRD'],gdata['HPRD']))
                    if 'OMIM' in gdata and gdata['OMIM']:  HTML.write('<A HREF="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s" TARGET=BODY_FRAME><IMG ALT="OMIM" SRC="../../resources/omim.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['OMIM'])
                except: pass    #self.errorLog('GeneMap data problem for %s' % sgene)
                HTML.write('</TD></TR>\n')
                try:
                    HTML.write('<TD>\n')
                    gdata = self.obj['GeneMap'].dict['Data'][sgene]
                    if 'Desc' in gdata: HTML.write('<TR><TD COLSPAN=3><FONT SIZE=4 COLOR=#014359><P>%s</FONT></P></TD></TR>\n' % gdata['Desc'])
                except: pass
                HTML.write('</TD></TR>\n')
                HTML.write('<TR><TD COLSPAN=3><P>%s</P></TD></TR>\n' % self.goHTML(sgene)) 
                try: HTML.write('<TR><TD COLSPAN=3><P><B>%s<->%s Evidence:</B> %s</P></TD></TR>\n' % (gene,sgene,string.join(data[gene]['ppi'][sgene],'; ')))
                except: HTML.write('<TR><TD COLSPAN=3><P><B>%s<->%s Evidence:</B> Not in DB</P></TD></TR>\n' % (gene,sgene))
                HTML.write('</TABLE>\n')
                ## ~ [6b] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>%s Interactors</H2>\n' % sgene)
                for ppitype in self.list['RunTypes']:
                    if not ppitype in data[sgene]: continue
                    if not data[sgene][ppitype]: continue
                    HTML.write('<H3>%s</H3>\n' % ppitype)
                    ilist = []
                    for igene in rje.sortKeys(data[sgene][ppitype]):
                        ilist.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A>' % (igene,igene,igene))
                    HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))

                ## ~ [6c] ~ Occurrence Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>%s SLiM Occurrences in %s</H2>\n' % (gene,spoke))
                socc = {}
                if spoke in odb.index('Seq'):
                    for entry in odb.entryList(odb.index('Seq')[spoke]):
                        if entry['Hub'] == gene: socc[odb.makeKey(entry)] = entry
                if socc:
                    self.writeHTMLTable(HTML,self.list['OccHead'],socc,asdict=True,gene=gene,drop=['Hub','Seq','Desc'])
                    #!# In future, intersperse #!#
                    vis = []
                    for okey in rje.sortKeys(socc):
                        (start,end) = (string.atoi(socc[okey]['Start_Pos']),string.atoi(socc[okey]['End_Pos']))
                        i = rje.preZero(int(start/100) * 100 + 1,99999)
                        j = rje.preZero(int(end/100) * 100 + 110,99999)
                        png = '../../spoke/%s/%s.%s-%s.png' % (sgene,spoke,i,j)
                        if png not in vis: vis.append(png)
                    vis.sort()
                    HTML.write('<H2>Summary Graphics</H2>\n')
                    for png in vis: HTML.write('<A HREF="%s" TARGET=_blank><IMG SRC="%s" ALT="%s" WIDTH=100%%></A>\n' % (png,png,os.path.basename(png)))
                else:
                    HTML.write('<P><I>No predicted %s-interacting SLiMs in %s</I></P>\n' % (gene,spoke))
                ## ~ [6d] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('</BODY></HTML>')
                HTML.close()

            ### ~ [7] ~ Hub-Motif Summary Page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gene not in odb.index('Hub'): return
            for sdata in odb.entryList(odb.index('Hub')[gene]):
                slim = rje_slim.slimFromPattern(sdata['Pattern'])
                hpage = hdir + '%s.%s.html' % (gene,slim)
                if not self.opt['Force'] and self.checkHTML(hpage): continue
                HTML = open(hpage,'w')
                HTML.write('<HTML><HEAD><TITLE>%s %s</TITLE>\n' % (gene,sdata['Pattern']))
                HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../../resources/%s"/>\n' % css)
                HTML.write('</HEAD>\n')
                ## ~ [7a] ~ Summary information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<BODY><H1>%s in %s</H1>\n' % (sdata['Pattern'],sdata['Dataset']))
                self.writeHTMLTable(HTML,self.list['SumHead'],[sdata],asdict=False,gene=gene,drop=['Dataset'])
                ## ~ [7b] ~ HubList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>Hubs returning %s</H2>\n' % (sdata['Pattern']))
                shubs = []
                for shub in rje.sortKeys(self.dict['SlimHubs'][sdata['Pattern']]['ALL']):
                    shared = 0
                    for p2 in data[shub]['ppi']:
                        if self.dict['EnsLoci'][p2] and self.dict['EnsLoci'][p2] in self.dict['SlimSpokes'][sdata['Pattern']]: shared += 1
                    shubs.append('<A HREF="../%s/%s.html" TARGET=_top>%s</A> [%s] (%s)' % (shub,shub,shub,shared,self.dict['SlimHubs'][sdata['Pattern']]['ALL'][shub]))
                HTML.write('<P>%s;</P>\n' % string.join(shubs,'; '))
                ## ~ [7c] ~ Graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>Summary Graphics</H2>\n')
                HTML.write('<A HREF="%s.%s.%s.png" TARGET=_blank><IMG SRC="%s.%s.%s.png" ALT="Summary Graphics" WIDTH=100%%></A>\n' % (gene,sdata['Rank'],slim,gene,sdata['Rank'],slim))
                ## ~ [7d] ~ Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                socc = {}
                for entry in odb.entryList(odb.index('Pattern')[sdata['Pattern']]):
                    if entry['Hub'] == gene: socc[odb.makeKey(entry)] = entry
                HTML.write('<H2>%s Occurrences in %s-interactors</H2>\n' % (sdata['Pattern'],gene))    
                self.writeHTMLTable(HTML,self.list['OccHead'],socc,asdict=True,gene=gene,drop=['Dataset','Pattern'])
                HTML.write('</BODY></HTML>')
                HTML.close()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def writeHTMLTable(self,HTML,headers,data,asdict=True,gene=None,drop=[]):     ### Write data to HTML table
        '''Write data to HTML table.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = headers[0:]   # Do not want to change original headers list!
            for h in ['RunID','Masking','Build','RunTime'] + drop:
                if h in headers: headers.remove(h)
            if asdict:
                datalist = []
                for key in rje.sortKeys(data): datalist.append(data[key])
            else: datalist = data
            ### ~ [2] ~ Write Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rows = ['<TD><B>%s</B></TD>' % string.join(headers,'</B></TD><TD><B>')]
            ## ~ [2b] ~ Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for dat in datalist:
                rdata = []
                for h in headers:
                    if h in dat:
                        if h in ['Cons','GlobID','LocID','Hyd','IUP','SA','IC','ExpUP'] or h[-4:] == 'mean':
                            try: rdata.append('%.3f' % string.atof(dat[h]))
                            except: rdata.append('-')
                        elif gene and h == 'Seq':
                            #rdata.append('<A HREF="%s.%s.html" TARGET=BODY_FRAME>%s (%s)</A>' % (gene,dat[h],string.split(dat[h],'__')[0],string.split(dat[h],'__')[1]))
                            #!# Replaced with Gene in table #!#
                            rdata.append('<A HREF="%s.%s.html" TARGET=BODY_FRAME>%s</A>' % (gene,dat[h],self.dict['SeqMap'][dat[h]][0]))
                        elif h == 'Seq':
                            rdata.append('<A HREF="../../unifake/%s.html" TARGET=BODY_FRAME>%s</A>' % (dat[h],self.dict['SeqMap'][dat[h]][0]))
                        elif h == 'Pattern' and dat['Pattern'] == '-': rdata.append('-')
                        elif gene and h == 'Pattern':
                            slim = rje_slim.slimFromPattern(dat['Pattern'])
                            rdata.append('<A HREF="%s.%s.html" TARGET=BODY_FRAME>%s</A>' % (gene,slim,dat[h]))
                        elif h == 'Pattern' and not gene and 'Dataset' in dat:
                            hgene = string.replace(dat['Dataset'],self.info['Suffix'],'')
                            slim = rje_slim.slimFromPattern(dat['Pattern'])
                            rdata.append('<A HREF="../../hub/%s/%s.%s.html" TARGET=BODY_FRAME>%s</A>' % (hgene,hgene,slim,dat[h]))
                        elif h == 'Dataset' and not gene:
                            hgene = string.replace(dat['Dataset'],self.info['Suffix'],'')
                            rdata.append('<A HREF="../../hub/%s/%s.html" TARGET=_top>%s</A>' % (hgene,hgene,dat[h]))
                        elif h in ['Rank','Start_Pos','End_Pos']: rdata.append('%d' % string.atoi(dat[h]))
                        else: rdata.append(dat[h])
                    else: rdata.append('')
                rows.append('<TD>%s</TD>' % string.join(rdata,'</TD><TD>'))
            ## ~ [2c] ~ Write to HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            HTML.write('<TABLE BORDER=2><TR>\n%s\n</TR></TABLE>\n' % string.join(rows,'\n</TR><TR>\n'))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def uniFakeHTML(self):  ### Generate composite HTML-linked DAT file                                             #V2.0
        '''Generate composite HTML-linked DAT file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = self.dict['Data']
            mdb = self.obj['DB'].getTable('Main')
            odb = self.obj['DB'].getTable('Occ')
            hdir = rje.makePath('%sunifake/' % self.info['HTMLDir'])
            rje.mkDir(self,hdir)
            css = self.info['HTMLCSS']
            ### ~ [2] ~ Generate DAT File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx, stot) = (0.0, len(data))
            for spoke in rje.sortKeys(data):
                ## ~ [2a] ~ Generate general attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#UNI','Generating HTML-linked DAT files: %.2f%%' % (sx/stot)); sx += 100.0
                try: spokeseq = self.dict['EnsLoci'][spoke]
                except: continue
                hfile = '%s%s.html' % (hdir,spokeseq)
                if not self.opt['Force'] and self.checkHTML(hfile): continue
                fakeacc = string.split(spokeseq,'__')[1]
                realacc = None
                try:
                    gdata = self.obj['GeneMap'].dict['Data'][spoke]
                    if 'UniProt' in gdata and gdata['UniProt']: realacc = gdata['UniProt']
                except: pass
                ## ~ [2b] ~ Compile UniProt data and add Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    if not self.list['UniReal']: realacc = None
                    elif not realacc: self.printLog('#UNI','No Real AccNum for %s' % seq)
                    shhh = self.log.opt['Silent']
                    self.log.opt['Silent'] = True
                    realuni = rje_uniprot.UniProt(self.log,self.cmd_list)
                    if realacc: realuni.readUniProt(clear=True,acclist=[realacc],cleardata=False)
                    fakeuni = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=%s' % self.info['UniFake']])
                    fakeuni.readUniProt(clear=True,acclist=[fakeacc],cleardata=False)
                    self.log.opt['Silent'] = shhh
                    dat = fakeuni.list['Entry'][0]
                    #self.deBug(dat.dict)
                    sequence = dat.obj['Sequence'].info['Sequence'][0:]
                    ## Map and Add Features from actual UniProt entry ##
                    for entry in realuni.list['Entry']:
                        #self.deBug(entry.dict['Data'])
                        for key in self.list['UniReal']:
                            if key != 'FT' and entry.dict['Data'].has_key(key):
                                if dat.dict['Data'].has_key(key): dat.dict['Data'][key] = entry.dict['Data'][key] + dat.dict['Data'][key]
                                else: dat.dict['Data'][key] = entry.dict['Data'][key][0:]
                        for ft in entry.list['Feature'][0:]:
                            if 'FT' not in self.list['UniReal']: break
                            ft_start = ft['Start']
                            ft_end = ft['End']
                            ft_seq = entry.obj['Sequence'].info['Sequence'][ft_start-1:ft_end]
                            if ft_seq == sequence[ft_start-1:ft_end]: dat.list['Feature'].append(ft); continue
                            fudge = 1
                            while fudge:
                                if ft_start - fudge < 1 and ft_end + fudge > len(sequence): fudge = 0; break
                                if ft_start - fudge >= 1 and ft_seq == sequence[ft_start-1-fudge:ft_end-fudge]: fudge = -fudge; break
                                if ft_end + fudge <= len(sequence) and ft_seq == sequence[ft_start-1+fudge:ft_end+fudge]: break
                                fudge += 1
                            if fudge:
                                ft['Start'] = ft_start + fudge
                                ft['End'] = ft_end + fudge
                                dat.list['Feature'].append(ft)
                except:
                    self.errorLog('Problem making HTML-linked DAT file for %s' % spokeseq)
                    continue
                ## ~ [2c] ~ Add occurrences as UniProt features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for occ in odb.entryList(odb.index('Seq')[spokeseq]):
                    ftdic = {'Type':'SLIM','Start': string.atoi(occ['Start_Pos']),'End':string.atoi(occ['End_Pos']),
                            'Desc':'%s Rank %s SLiM %s (p=%s)' % (string.replace(occ['Dataset'],self.info['Suffix'],''),occ['Rank'],occ['Pattern'],occ['Sig'])}
                    dat.list['Feature'].append(ftdic)
                ## ~ [2d] ~ Save data to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML = open(hfile,'w')
                HTML.write('<HTML><HEAD><TITLE>%s</TITLE>\n' % spoke)
                HTML.write('<LINK REL="stylesheet" TYPE="text/css" HREF="../resources/%s"/>\n' % css)
                HTML.write('</HEAD>\n')
                ## ~ [6a] ~ Gene summary as with hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<BODY>\n')
                HTML.write('<TABLE WIDTH=100%% BORDER=0>\n<TR><TD><A HREF="../hub/%s/%s.html" TARGET=_top><FONT SIZE=6 FACE="Verdana" COLOR=#014359><B>%s</B></FONT></A></TD>\n' % (spoke,spoke,spoke))
                HTML.write('<TD><FONT SIZE=5 FACE="Verdana" COLOR=#979E45><B>%s</P></TD>\n' % (spokeseq))  #!# Add UniFake link #!#
                try:
                    HTML.write('<TD>\n')
                    gdata = self.obj['GeneMap'].dict['Data'][spoke]
                    if 'Symbol' in gdata and gdata['Symbol']: HTML.write('<A HREF="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" TARGET=BODY_FRAME><IMG ALT="GeneCards" SRC="../resources/genecard.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['Symbol'])
                    if 'UniProt' in gdata and gdata['UniProt']: HTML.write('<A HREF="http://www.uniprot.org/uniprot/%s" TARGET=BODY_FRAME><IMG ALT="UniProt" SRC="../resources/logo_ebi.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['UniProt'])
                    if 'EnsEMBL' in gdata and gdata['EnsEMBL']: HTML.write('<A HREF="http://www.ensembl.org/Homo_sapiens/geneview?gene=%s" TARGET=BODY_FRAME><IMG ALT="EnsEMBL" SRC="../resources/e-bang.gif" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['EnsEMBL'])
                    if 'HPRD' in gdata and gdata['HPRD']:  HTML.write('<A HREF="http://www.hprd.org/summary?protein=%s&isoform_id=%s_1&isoform_name=Isoform_1" TARGET=BODY_FRAME><IMG ALT="HPRD" SRC="../resources/hprd.png" ALIGN=BOTTOM BORDER=0></A>\n' % (gdata['HPRD'],gdata['HPRD']))
                    if 'OMIM' in gdata and gdata['OMIM']:  HTML.write('<A HREF="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s" TARGET=BODY_FRAME><IMG ALT="OMIM" SRC="../resources/omim.png" ALIGN=BOTTOM BORDER=0></A>\n' % gdata['OMIM'])
                except: self.errorLog('GeneMap data problem for %s' % spoke)
                HTML.write('</TD></TR>\n')
                try:
                    HTML.write('<TD>\n')
                    gdata = self.obj['GeneMap'].dict['Data'][spoke]
                    if 'Desc' in gdata: HTML.write('<TR><TD COLSPAN=3><FONT SIZE=4 COLOR=#014359><P>%s</FONT></P></TD></TR>\n' % gdata['Desc'])
                except: pass
                HTML.write('</TD></TR>\n')
                HTML.write('</TABLE>\n')
                ## ~ [6b] ~ Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>%s Interactors</H2>\n' % spoke)
                ilist = []
                for p2 in rje.sortKeys(data[spoke]['ppi']):
                    if 'seq' in data[p2]: ilist.append('<A HREF="../hub/%s/%s.html" TARGET=_top>%s</A>' % (p2,p2,p2))
                    else: ilist.append('%s' % (p2))
                HTML.write('<P>%s;</P>\n' % string.join(ilist,'; '))
                ## ~ [6c] ~ Occurrence Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                HTML.write('<H2>UniProt Format</H2><PRE>\n')
                entry = dat
                eseq = entry.obj['Sequence']
                #self.deBug(entry.dict['Data'])
                ## Standard info ##
                for key in ['ID','AC','DT','DE','GN','OS']:
                    if entry.dict['Data'].has_key(key):
                        for rest in entry.dict['Data'][key]: HTML.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## Other data, except Features and sequence ##
                for key in rje.sortKeys(entry.dict['Data']):
                    if key not in ['ID','AC','DT','DE','GN','OS','FT','SQ','SEQ','//']:
                        for rest in entry.dict['Data'][key]: HTML.write('%s   %s\n' % (key,rje.chomp(rest)))
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
                    ftxt += ' %s\n' % self.linkFT(ftdict['Type'],ftdict['Desc'],spokeseq)
                    HTML.write(ftxt)
                ## Sequence/End ##
                HTML.write('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % eseq.aaLen())),eseq.aaLen(),rje_sequence.MWt(eseq.info['Sequence'])))
                uniseq = eseq.info['Sequence'][0:]
                while len(uniseq) > 0:
                    HTML.write('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                    uniseq = uniseq[60:]
                HTML.write('//\n</PRE></BODY></HTML>')
                HTML.close()
                self.printLog('#UNI','Generation of %s HTML-linked DAT file complete.' % spoke,screen=False)
            self.printLog('\r#UNI','Generation of HTML-linked DAT files complete.')
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def linkFT(self,type,desc,seq):  ### Adds hyperlinks to feature text
        '''Adds hyperlinks to feature text.'''
        try:### ~ [1] ~ Hyperlink SLiMs and datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if type == 'SLIM':
                dsplit = string.split(desc)
                dset = dsplit[0]
                pattern = dsplit[4]
                slim = rje_slim.slimFromPattern(pattern)
                dsplit[0] = '<A HREF="../hub/%s/%s.%s.html" TARGET=BODY_FRAME>%s</A>' % (dset,dset,seq,dsplit[0])
                dsplit[4] = '<A HREF="../hub/%s/%s.%s.html" TARGET=BODY_FRAME>%s</A>' % (dset,dset,slim,dsplit[4])
                return string.join(dsplit)
        except: self.errorLog(rje_zen.Zen().wisdom())
        return desc
#########################################################################################################################
    ### <5> ### Single Data SLiMJIM PNG generation Methods                                                              #
#########################################################################################################################
    def singleSF(self):     ### Generates graphics for single SLiMFinder dataset
        '''Generates graphics for single SLiMFinder dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.info['JIMDir'])
            if self.opt['Iridis']: os.system('module load R/2.7.1/gcc-3.4.4')
            self.opt['Iridis'] = False
            self.info['Basefile'] = rje.baseFile(self.info['SingleSF'])
            ### ~ [2] ~ Generate PNGs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.singleSFMotifJIM()     # Generates motif alignment, profile and tree
            self.singleSFMappingJIM()   # Generates SLiMJIM spoke alignment PDFs from *.mapping.fas
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def singleSFMotifJIM(self): ### Peforms the SLiMJim motif visualisations for a given SLiMFinder dataset
        '''Peforms the SLiMJim motif visualisations for a given SLiMFinder dataset.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disfile = '%s.dis.tdt' % (self.info['Basefile'])        # Make tree from this
            occfile = '%s.occ.csv' % (self.info['Basefile'])        # Get hit data from this
            seqfile = '%s.slimdb' % (self.info['Basefile'])         # Get extra seq data from this
            motifaln = '%s.motifaln.fas' % (self.info['Basefile'])  # Get regional alignment from this
            basedir = rje.makePath('%s/%s' % (self.info['JIMDir'],rje.baseFile(seqfile,True)))
            rje.mkDir(self,basedir)
            for file in [disfile,occfile,seqfile,motifaln]:
                if not os.path.exists(file):
                    self.errorLog('Cannot generate motif PNG without %s' % file)
                    return False
            ## ~ [1a] ~ Load general data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scmd = ['seqin=%s' % seqfile,'autoload=T','replacechar=F']
            dseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            seqdic = dseq.seqNameDic()
            occdic = self.loadOccData(occfile,['Dataset','Rank','Pattern','Seq','Start_Pos','End_Pos'])
            occdic.pop('Headers')
            disdat = rje.dataDict(self,disfile,['SEQ'])
            ## ~ [1b] ~ Individual Motif alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            motdic = {}     # Dictionary of {slimseq:[seqs]}
            scmd = ['seqin=%s' % motifaln,'autoload=T','seqnr=F','accnr=F','replacechar=F']
            mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)
            fullseq = mseq.seq[0:]
            rank = 1
            malnlist = string.split(mseq.seq[0].info['Sequence'],'-XXXXXXXXXX-')
            for maln in malnlist:
                motdic[maln] = [mseq._addSeq(mseq.seq[0].shortName(),maln)]
                slim = string.replace(maln[10:],'-','')
                for seq in fullseq[1:]:
                    sequence = string.split(seq.info['Sequence'],'-XXXXXXXXXX-')[rank-1]
                    if string.replace(sequence,'-',''): motdic[maln].append(mseq._addSeq(seq.shortName(),sequence))
                rank += 1
            mseq.seq = fullseq[0:]
            self.printLog('#SEQ','%d motif alignment read from motifaln' % (len(motdic)))

            ### ~ [2] ~ Work through Motif Alignments in turn
            for maln in malnlist:
                mseq.seq = motdic[maln][0:]
                pattern = string.replace(maln[10:],'-','')
                slim = rje_slim.slimFromPattern(pattern)
                rank = malnlist.index(maln) + 1
                basefile = '%s%s.%s.%s' % (basedir,rje.baseFile(seqfile,True),rje.preZero(rank,len(malnlist)),slim)
                ## ~ [2a] ~ Modify Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                msequence = mseq.seq[0].info['Sequence']
                newm = []
                r = 0
                while r < len(msequence):
                    if msequence[r] != '[': newm.append(msequence[r]); r += 1
                    else:
                        amb = rje.matchExp('^\[([A-Z]+)\]',msequence[r:])[0]
                        newm.append(amb)
                        r += len(amb)+2
                while len(newm) < len(msequence): newm.append('-')
                mseq.seq[0].info['Sequence'] = newm
                ## ~ [2b] ~ Draw tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                dismat = rje_dismatrix_V2.DisMatrix(self.log,self.cmd_list)
                reduced = []    # Reduced sequence set, ordered according to UPC
                for seq1 in mseq.seq[1:]:
                    obj1 = seq1.shortName()
                    for seq2 in mseq.seq[1:]:
                        obj2 = seq2.shortName()   
                        dis = (string.atof(disdat[obj1][obj2]) + string.atof(disdat[obj2][obj1])) / 2.0
                        dismat.addDis(obj1,obj2,dis)
                upgma = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
                nsftree = dismat.upgma()
                open('%s.nsf' % basefile,'w').write(nsftree)
                upgma.buildTree(nsftree,type='nsf',postprocess=False)
                if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
                upgma.rTree('%s.tree.csv' % basefile,seqname='short')
                reduced = upgma._vertOrder(internal=False,namelist=True)
                #if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
                #dismat.saveMatrix(reduced,basefile+'.heatmap.tdt',delimit='\t')
                ## ~ [2c] ~ Reorder sequences allowing for duplications ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newseq = [mseq.seq[0]]
                for seq in mseq.seq[1:]:
                    i = 1
                    while i < len(newseq):
                        if reduced.index(seq.shortName()) < reduced.index(newseq[i].shortName()): break
                        i += 1
                    newseq.insert(i,seq)
                mseq.seq = newseq[0:]
                ## ~ [2d] ~ Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if os.path.exists('%s.aln.tdt' % basefile): os.unlink('%s.aln.tdt' % basefile)
                mseq.saveR(seqfile='%s.aln.tdt' % basefile,name='AccNum')
                ## ~ [2e] ~ Make a list of profile sequences +/- 5 of motif ~~~~~~~~~~~~~~~~~~~~~~~ ##
                proseq = []
                for okey in rje.sortKeys(occdic):
                    odata = occdic[okey]
                    if odata['Pattern'] != pattern or rank != string.atoi(odata['Rank']): continue
                    seq = seqdic[odata['Seq']]
                    oseq = seq.info['Sequence'][max(0,string.atoi(odata['Start_Pos'])-6):string.atoi(odata['End_Pos'])+5]
                    c = max(0,6 - string.atoi(odata['Start_Pos']))
                    oseq = '-' * c + oseq
                    n = max(0,string.atoi(odata['End_Pos'])+5-seq.aaLen())
                    oseq = oseq + '-' * n
                    proseq.append(oseq)
                ### ~ [2f] ~ Convert to profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                poslist = ['-5','-4','-3','-2','-1']
                for pos in rje_slim.prestoFromCode(rje_slim.slimFromPattern(slim)):
                    if len(pos) > 1: poslist.append('[%s]' % pos)
                    elif pos == 'X': poslist.append('x')
                    else: poslist.append(pos)
                poslist += ['1','2','3','4','5']
                ignorepos = ['-5','-4','-3','-2','-1','x','1','2','3','4','5']
                outfile = basefile + '.profile.tdt'
                if os.path.exists(outfile): os.unlink(outfile)
                ### ~ [2g] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                headers = ['Pos'] + rje_seq.alph_protx[:-1]
                rje.delimitedFileOutput(self,outfile,headers,'\t')
                try:
                    for r in range(len(poslist)):
                        datadict = {}
                        for a in headers[1:]: datadict[a] = 0.0
                        for pseq in proseq:
                            while len(pseq) < len(poslist) and poslist[r] not in ignorepos and poslist[r].find(pseq[r]) < 0: pseq = pseq[:r] + '-' + pseq[r:]
                            try:
                                a = pseq[r]
                                datadict[a] += 1
                            except: pass    # Not counting Xs
                        datadict['Pos'] = poslist[r]
                        rje.delimitedFileOutput(self,outfile,headers,'\t',datadict)
                except:
                    self.errorLog('Problem during profile output for %s' % os.path.basename(basefile))
                    for pseq in proseq: print pseq
                ## ~ [2h] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rcmd = '%s --no-restore --no-save --args "sfslim2png" "%s"' % (self.info['RPath'],basefile)
                rslimjim = '%srje.r' % self.info['Path']
                rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
                problems = self.rCall(rcmd,basefile)
                ## ~ [6a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not self.opt['Iridis']: 
                    for ext in ['tree.csv','aln.tdt','profile.tdt','r.tmp.txt','nsf']:
                        if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def singleSFMappingJIM(self):   ### Generate SLiMJIM spoke alignment PDFs from *.mapping.fas
        '''Generate SLiMJIM spoke alignment PDFs from *.mapping.fas.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basedir = rje.makePath('%s/%s' % (self.info['JIMDir'],os.path.basename(self.info['Basefile'])))
            mapseq = {}     # Dictionary of {dataset:[seqs]}
            scmd = ['autoload=T','seqnr=F','accnr=F','replacechar=F','seqin=%s.mapping.fas' % self.info['Basefile']]
            mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)     #!# Removed ['minregion=3']+ #!#
            while mseq.seq:
                ## ~ [1a] ~ Read in all sequences for one spoke ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pseq = [mseq.seq.pop(0)]  # Pseq = list of sequences for this protein
                while mseq.seq:
                    if mseq.seq[0].info['Name'].find('Motifs') > 0 and string.split(mseq.seq[0].info['Name'])[1] == 'Motifs': break   # Next protein
                    pseq.append(mseq.seq.pop(0))
                ## ~ [1b] ~ Update relevant sequence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                mapseq[pseq[0].shortName()] = pseq[0:]
            self.printLog('#ALN','%d distinct alignments identified' % len(mapseq))
            ### ~ [2] ~ Make SLiMJIM visualisations for each protein  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0  # Number of errors
            for mapping in rje.sortKeys(mapseq):
                try:
                    ## ~ [3a] ~ Rename sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    pseq = mapseq[mapping][0:]
                    basefile = basedir + pseq[0].shortName()
                    if self.interactive() > 0 and not rje.yesNo(pseq[0].shortName()): continue
                    qryname = pseq[2].shortName()
                    pseq[0].info['R'] = pseq[0].shortName()[len(qryname)+1:]
                    pseq[1].info['R'] = 'Masked'
                    for seq in pseq[2:]: seq.info['R'] = seq.info['ID']
                    ## ~ [3b] ~ Setup new SeqList, strip Query gaps, calculate RelCons ~~~~~~~~~~~~~~~~ ##
                    seqfile = '%s.aln.tdt' % basefile
                    if os.path.exists(seqfile): os.unlink(seqfile)
                    rseq = rje_seq.SeqList(self.log,self.cmd_list+scmd+['autoload=F'])
                    rseq.seq = pseq
                    rseq.obj['QuerySeq'] = pseq[2]
                    rseq.tidyQueryGaps()
                    rseq.saveR(rseq.seq,seqfile,name='R')
                    rseq.seq = pseq[2:]
                    relfile = '%s.rel.tdt' % basefile
                    if os.path.exists(relfile): os.unlink(relfile)
                    rseq.relCons(relfile)
                    self.deBug(rseq.obj['QuerySeq'].cmd_list)
                    ## ~ [3c] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    rcmd = '%s --no-restore --no-save --args "sfmap2png" "%s"' % (self.info['RPath'],basefile)
                    rslimjim = '%srje.r' % self.info['Path']
                    rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
                    problems = self.rCall(rcmd,basefile)
                except: self.errorLog('SLiMJIM visualisation error for "%s"' % mapping); ex += 1
            self.printLog('#SLIMJIM','Generation of SLiMJIMs complete. %d Problems.' % ex)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <6> ### IRIDIS PNG generation Methods                                                                           #
#########################################################################################################################
                                        ## >>> See slimjim_V1 <<< ##
#########################################################################################################################
    ### <7> ### SLiMJIM Table generation Methods                                                                        #
#########################################################################################################################
                                        ## >>> See slimjim_V1 <<< ##
#########################################################################################################################
### End of SECTION II: SLiMJIM Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
class GeneMap(rje_genemap.GeneMap):     # For unpickling?
    '''No methods to be added here. See rje_genemap.'''
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
    try: SLiMJIM(mainlog,cmd_list).run()
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
