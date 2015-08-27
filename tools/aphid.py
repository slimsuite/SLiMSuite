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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      APHID
Description:  Automated Processing of High-resolution Intensity Data
Version:      2.2
Last Edit:    10/07/14
Citation:     Raab et al. (2010), Proteomics 10: 2790-2800. [PMID: 20486118]
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This modules takes for input the partially processed results of MS analysis, with intensity data, filters based on
    scores thresholds, removes redundancy (using PINGU) and calculates relative intensity scores. PINGU is then used to
    generate outputs for use with Cytoscape and other visualisation tools.

    Input takes the form of a delimited text file with the following column headers: Expt, Subpop, Identifier,
    logInt, Score. In addition, other columns may be present (and may be used to filter data). The "unique" column headers allow
    individual identifications to be isolated, which is important for data filtering and intensity mapping.

    Intensities are converted into relative intensities with two options: (1) the redundancy level determines which
    results are combined. This may be simply at the identified peptide level, the gene level, or even at the protein
    family level (as determined through BLAST homology). If "fam" is used then GABLAM will be used to generate families
    for a given level of identity.
    
    **Note.** Core functionality has not been checked/developed since 2010.

Commandline:
    ### ~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    data=FILE       : Delimited file containing input data [None]
    seqin=FILE      : Sequence file containing hits (if pepnorm != count) [None]
    basefile=X      : Text base for output files and resdir [default based on data file name]
    unique=LIST     : Headers that, with Identifier, Treatment & Replicate, constitute unique entries [Slice]
    identifier=X    : Column containing protein indentifications [Identifier]
    treatment=X     : Column heading to identify different experimental treatment samples (e.g. condition/control) [Treatment]
    replicate=X     : Column heading to identify experimntal replicates [Replicate]
    intensity=X     : Column containing intensity values [logInt]
    logint=T/F      : Whether intensity value is a log intensity [True]
    pepcount=X      : Column containing peptide counts [rI]
    statfilter=LIST : List of filters for data [Score>-3]
    force=T/F       : Regenerate intermediate files even if found [False]

    ### ~ Intensity Combination Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    nr=X            : Level for redundancy removal (gene/pep/fam) [gene]
    famcut=X        : Percentage identity to be used by GABLAM for clustering [0.0]
    normalise=X     : Scoring strategy for normalising intensity (ppm/shared/none/mwt/frag/count) [ppm]
    flatout=T/F     : Whether to output reduced GeneMap flatfiles (*.data.tdt & *.aliases.tdt) [False]

    ### ~ Enrichment options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    absence=X       : Min normalised combined score (arbitrarily assigned to proteins totally absent from samples) [0.5]
    combine=X       : Methods for combining different replicates (max,mean,min,geo,full) [mean]
    convert=LIST    : List of X:Y where Treatment X will be renamed Y []
    enrpairs=LIST   : List of X:Y where enrichment will be restricted to X:Y ratios []
    jackknife=T/F   : Whether to perform jack-knifing tests on enrpairs [True]
    enrcut=X        : Enrichment > X for jack-knifing test [1.0]
    blanks=T/F      : Whether to include most blank enrichment/jacknife rows for missing NRID [True]

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import math, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genemap, rje_scoring, rje_seq, rje_slimcore
import gablam, pingu_V3
import rje_dismatrix_V2 as rje_dismatrix
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added fam NR using GABLAM and SLiMCore.
    # 0.2 - Added replicates and pepnum abundance.
    # 0.3 - Improvements incorporating GeneMap version of Pingu.
    # 0.4 - Allowed mutliple runs from one call. Added fullmonty option. Moved and corrected compare mode.
    # 1.0 - Full initial working version.
    # 1.1 - Added PNG visualisations.
    # 1.2 - Added additional filtering according to replicate counts (minrep=X).
    # 1.3 - Added additional filtering according to min peptide occurrences (minpep=X,Y)
    # 2.0 - Reduced options and tidied code substantially with additional intermediate data outputs.
    # 2.1 - Reduced import commands.
    # 2.2 - Updated for revised SLiMCore.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add fam NR - GABLAM and use SLiMCore UPC maker.
    # [Y] : Check additional replicates
    # [Y] : Check peptide number abundance measures
    # [X] : Add filtering of Identifications not in Sequence File
    # [ ] : Check redundancy by family is working
    # [ ] : Fix PNG output?
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('APHID', '2.2', 'July 2014', '2007')
    description = 'Automated Processing of High-resolution Intensity Data'
    author = 'Dr Richard J. Edwards.'
    comments = [#'Note: R is required for PNG visualisations.',
                'This program is still in development and has not been published.'] #,rje_zen.Zen().wisdom()]
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
            if rje.yesNo('Show PINGU commandline options?'): out.verbose(-1,4,text=pingu.__doc__)
            if rje.yesNo('Show BLAST commandline options (for GABLAM)?'): out.verbose(-1,4,text=pingu.__doc__)
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
### SECTION II: APHID Class                                                                                             #
#########################################################################################################################
class APHID(rje.RJE_Object):     
    '''
    APHID Class. Author: Rich Edwards (2007).

    Info:str
    - Data = Delimited file containing input data [None]
    - Identifier = Column containing protein indentifications [Identifier]
    - Treatment = Column identifier to be used to identify different samples [Treatment]
    - Replicate = Column heading used to identify replicates [Replicate]
    - Intensity = Column containing intensity values [logInt]
    - PepCount = Column containing peptide counts [rI]
    - Combine = Methods for combining different replicates (max,mean,min,geo) [mean]
    - NR = Level for redundancy removal (gene/pep/fam) [gene]
    - Normalise = Scoring strategy for normalising intensity (ppm/shared/none/mwt/frag/count) [ppm]
    
    Opt:boolean
    - Blanks = Whether to include most blank enrichment/jacknife rows for missing NRID [True]
    - FlatOut = Whether to output reduced GeneMap flatfiles (*.data.tdt & *.aliases.tdt) [False]
    - JackKnife = Whether to perform jack-knifing tests on enrpairs [True]
    - LogInt = Whether intensity value is a log intensity [True]

    Stat:numeric
    - Absence = Log intensity score to be arbitrarily assigned to proteins/genes absent from Treatments [4.0]
    - EnrCut = Enrichment > X for jack-knifing test [1.0]
    - FamCut = Percentage identity to be used by GABLAM for clustering [0.0]

    List:list
    - StatFilter = list of filters for data [Score>-3]
    - Unique = List of column headers that, together with Identifier, constitute unique entries [Treatment,Band,HCT_No]
    - Convert = List of X:Y where Treatment X will be renamed Y []
    - EnrPairs = List of X:Y where enrichment will be restricted to X:Y ratios []

    Dict:dictionary
    - StatFilter = Manipulated dictionary from Statfilter list
    
    Obj:RJE_Objects
    - DB = rje_db.Database object used for loading and manipulating data
    - GeneMap = rje_genemap object used for ID mapping
    - SeqList = Sequence Object 
    - PINGU = Pingu object containing hits (if pepnorm != none), i.e. used for MS search
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Data','Identifier','Treatment','Replicate','Intensity','PepCount','Combine','NR',
                         'Normalise']
        self.optlist = ['LogInt','JackKnife','FlatOut','Blanks']
        self.statlist = ['Absence','FamCut','EnrCut']
        self.listlist = ['StatFilter','Unique','Convert','EnrPairs']
        self.dictlist = ['StatFilter']
        self.objlist = ['DB','GeneMap','SeqList','PINGU']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'PepCount':'rI','Intensity':'logInt','Combine':'mean','NR':'gene','Normalise':'ppm'})
        for i in ['Identifier','Treatment','Replicate']: self.info[i] = i
        self.setOpt({'LogInt':True,'JackKnife':True,'FlatOut':True,'Blanks':True})
        self.setStat({'Absence':0.5,'FamCut':0.0,'EnrCut':1.0})
        self.setList({'Unique':['Slice'],'StatFilter':['Score>-3']})
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
                self._cmdReadList(cmd,'file',['Data'])
                self._cmdReadList(cmd,'info',['Identifier','Treatment','Replicate','Intensity','PepCount','Combine',
                                              'NR','Normalise'])
                self._cmdReadList(cmd,'stat',['Absence','FamCut','EnrCut'])
                self._cmdReadList(cmd,'opt',['LogInt','JackKnife','FlatOut','Blanks'])
                self._cmdReadList(cmd,'list',['StatFilter','Unique','Convert','EnrPairs'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Run Methods                                                                                        #
#########################################################################################################################
    def run(self):  ### Main run method.
        '''Main run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] Convert to NR identifiers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.makeNR()
            ### ~ [3] Normalise intensities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.normalise()
            self.combine()
            ### ~ [4] Calculate Enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.enrichment()
            ### ~ [5] Jack-knifing assessment of enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.jackKnife()
            self.bootstrap()
        except: self.errorLog('APHID run error')
#########################################################################################################################
    def setup(self):    ### Checks input data and reads into database and filters                                   |2.0|
        '''Checks input data and reads into self attributes. Returns True/False depending on whether OK.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list+['autoload=T'])
            self.checkInputFiles(['Data'],ask=True)
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = rje.baseFile(self.info['Data'])
            self.obj['GeneMap'] = rje_genemap.GeneMap(self.log,self.cmd_list+['flatout=F','pickleout=F'])
            self.obj['GeneMap'].run()
            ### ~ [1] ~ Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in ['Identifier','Treatment','Replicate']:
                if self.info[i].lower() not in ['','none'] + self.list['Unique']: self.list['Unique'].append(self.info[i])
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            maindb = db.addTable(self.info['Data'],self.list['Unique'],name='Main')
            maindb.info['Delimit'] = '\t'
            self.printLog('#DATA','%s datapoints read from %s' % (rje.integerString(maindb.entryNum()),self.info['Data']))
            ## ~ [1a] ~ StatFilter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['StatFilter'] = rje_scoring.setupStatFilter(self,maindb.fields(),filterlist=self.list['StatFilter'])
            data = rje_scoring.statFilter(self,self.data('Main'),self.dict['StatFilter'])
            self.printLog('#DATA','%s datapoints remain after filtering.' % (rje.integerString(len(data))))
            maindb.dict['Data'] = data
            ## ~ [1b] ~ Convert Treatments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Convert']:
                self.printLog('#CONV','%d Conversions: %s' % (len(self.list['Convert']),string.join(self.list['Convert'],', ')))
                treatment = self.info['Treatment']
                convert = {}; cx = 0
                for cpair in self.list['Convert']:
                    if len(string.split(cpair,':')) == 2: convert[string.split(cpair,':')[0]] = string.split(cpair,':')[1]
                for entry in maindb.entries():
                    if entry[treatment] in convert: entry[treatment] = convert[entry[treatment]]; cx += 1
                self.printLog('#CONV','%d %s conversions' % (cx,treatment))
            return True
        except: self.log.errorLog('Problem with APHID.setup()')
        return False
#########################################################################################################################
    ### <3> ### Gene Mapping/NR Methods                                                                                 #
#########################################################################################################################
    def geneMap(self,id):   ### Gets best gene symbol for gene and updates GeneMap dictionary                       |2.0|
        '''Gets best gene symbol for gene and updates GeneMap dictionary.'''
        try: return self.obj['GeneMap'].bestMap(id)
        except: self.errorLog('APHID.GeneMap error'); raise
#########################################################################################################################
    def makeNR(self):   ### Converts to NR identifiers and outputs expanded conversion table                        |2.0|
        '''Converts to NR identifiers and outputs expanded conversion table.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map = self.obj['GeneMap']
            map.dict['BestMap'] = {}
            id = self.info['Identifier']
            maindb = self.db().getTable('Main')
            fields = maindb.fields()
            for field in ['NRID','HGNC','EnsG']:
                fields.insert(fields.index(id)+1,field)
                maindb.makeField(fieldname=field)
            maindb.list['Fields'] = fields
            maindb.makeField(fieldname='Desc')
            self.info['NR'] = self.info['NR'].lower()
            if self.info['NR'] not in ['pep','gene']:
                self.printLog('#NR','NR by "%s" not supported - will use NR by gene' % self.info['NR'])
                self.info['NR'] = 'gene'
            ### ~ [1] Convert to NR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['UpdateMap'] = {'ARP11':'ACTR3B',
                                      'ENSG00000198618':'PPIA',
                                      'ENSG00000206505':'HLA-A',
                                      'ENSG00000206441':'HLA-A',
                                      'ENSP00000376189':'RPS27A',
                                      'ENSP00000376685':'SH3BGRL2',
                                      'ENSP00000378021':'PRDX5',
                                      'ENSP00000378455':'ABHD14B',
                                      'ENSP00000381136':'PNP',
                                      'ENSP00000381540':'RAB1A',
                                      'ENSP00000382839':'HSPA1B',
                                      'ENSP00000383494':'HLA-A',
                                      'ENSP00000383556':'ARHGDIA'}
            ## ~ [1a] Add NR IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            unmapped = []
            ex = 0.0; etot = maindb.entryNum()
            for entry in maindb.entries():
                self.progLog('\r#NR','NR mapping: %.2f%%' % (ex/etot)); ex += 100.0
                pep = entry[id]
                nrgene = self.geneMap(pep)
                if nrgene in self.dict['UpdateMap']:
                    map.addAlias(self.dict['UpdateMap'][nrgene],nrgene)
                    nrgene = self.dict['UpdateMap'][nrgene]
                if nrgene not in map.dict['Data'].keys() + unmapped:
                    unmapped.append(nrgene)
                    self.deBug('%s >> %s (%s)' % (pep,nrgene,map.aliases(pep,full=False)))
                try: entry['EnsG'] = map.dict['Data'][nrgene]['EnsEMBL']
                except: entry['EnsG'] = '-'
                try: entry['HGNC'] = map.dict['Data'][nrgene]['Symbol']
                except: entry['HGNC'] = '-'
                try:
                    entry['Desc'] = map.dict['Data'][nrgene]['Desc']
                    if not entry['Desc']: entry['Desc'] = nrgene
                except: entry['Desc'] = 'Unmapped peptide'
                if self.info['NR'] == 'gene': entry['NRID'] = nrgene
                elif self.info['NR'] == 'pep': entry['NRID'] = pep
                else: raise ValueError
                if nrgene in unmapped: self.deBug(entry)
            self.printLog('\r#NR','NR mapping %s entries complete.' % rje.integerString(ex))
            self.printLog('#MAP','%s IDs remain unmapped' % rje.integerString(len(unmapped)))
            self.deBug('...')
            if unmapped:
                unmapped.sort()
                UNMAP = open('%s.unmapped.txt' % self.info['Basefile'],'w')
                for unmap in unmapped: UNMAP.write(string.join([unmap]+map.aliases(unmap,full=False)+['\n']))
                UNMAP.close()
            ## ~ [1b] Save intermediate files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maindb.saveToFile('%s.data_nr_map.tdt' % self.info['Basefile'])
            pinhead = ['Identifier','Sample']; pinfile = '%s.pingu_genes.tdt' % self.info['Basefile']
            rje.delimitedFileOutput(self,pinfile,pinhead,rje_backup=True)
            for gene in rje.sortKeys(maindb.index('NRID',force=True)): rje.delimitedFileOutput(self,pinfile,pinhead,datadict={'Identifier':gene,'Sample':self.info['Basefile']})
            self.printLog('#PINGU','%s genes output to %s for PINGU' % (rje.integerString(len(maindb.index('NRID'))),pinfile))
            ## ~ [1c] Compress database table and save NR data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for field in self.list['Unique']:
                if field in [self.info['Treatment'],self.info['Replicate'],self.info['Identifier']]: continue
                maindb.deleteField(field)
            newkeys = [self.info['Treatment'],self.info['Replicate'],'NRID']
            rules = {self.info['Intensity']:'max', self.info['PepCount']:'max',
                     'Score':'min',id:'list','EnsG':'list','HGNC':'list'}
            maindb.compress(newkeys,rules,default='mean',best=[])
            maindb.saveToFile('%s.nr.tdt' % self.info['Basefile'])
            ## ~ [1d] Additional GeneMap output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['FlatOut']:
                self.obj['GeneMap'].info['Basefile'] = self.info['Basefile']
                self.obj['GeneMap'].reduceGeneData(maindb.index('NRID').keys())
                self.obj['GeneMap'].flatOut()
        except: self.errorLog('Major error with APHID.makeNR()')
#########################################################################################################################
    ### <4> ### Normalise/Combine Methods                                                                               #
#########################################################################################################################
    def normalise(self):   ### Normalises abundance within sample                                                   |2.0|
        '''Normalises abundance within sample.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            intensity = self.info['Intensity']; pepcount = self.info['PepCount']
            self.info['Normalise'] = self.info['Normalise'].lower()
            if self.info['Normalise'] not in ['ppm','none']:
                self.printLog('#NORM','"%s" normalisation not supported. Will use ppm.' % self.info['Normalise'])
                self.info['Normalise'] = 'ppm'
            maindb = self.db().getTable('Main')
            if 'Intensity' not in maindb.fields(): maindb.list['Fields'].append('Intensity')
            for entry in maindb.entries():
                if self.opt['LogInt']: entry['Intensity'] = math.pow(10,float(entry[intensity]))
                else: entry['Intensity'] = float(entry[intensity])
            sampcomb = '#%s#%s#%s#' % (self.info['Treatment'],maindb.info['Delimit'],self.info['Replicate'])
            samples = maindb.index(sampcomb,make=True)
            sampdb = self.db().copyTable(maindb,'Samples')
            rules = {'NRID':'list','Intensity':'sum', self.info['PepCount']:'sum'}
            newkeys = [self.info['Treatment'],self.info['Replicate']]
            for field in sampdb.fields():
                if field not in newkeys + rules.keys(): sampdb.deleteField(field)
            sampdb.compress(newkeys,rules)
            sampdb.makeField(fieldname='NRCount')
            sampdb.list['Fields'].remove('NRID')
            sampdb.list['Fields'].append('NRID')
            for key in sampdb.datakeys():
                try: sampdb.data()[key]['NRCount'] = len(samples[key])
                except: self.errorLog('Problem with NRCount for "%s"' % key)
            sampdb.saveToFile('%s.samples.tdt' % self.info['Basefile'])
            ### ~ [1] Calculate normalised intensity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maindb.makeField(fieldname='IntNorm'); maindb.makeField(fieldname='PepNorm');
            for entry in maindb.entries():
                try:
                    sample = '%s%s%s' % (entry[self.info['Treatment']],maindb.info['Delimit'],entry[self.info['Replicate']])
                    if sample not in sampdb.datakeys(): raise ValueError
                    if self.info['Normalise'] == 'ppm':
                        entry['IntNorm'] = 1e6 * float(entry['Intensity']) / float(sampdb.data()[sample]['Intensity'])
                        entry['PepNorm'] = 1e6 * float(entry[pepcount]) / float(sampdb.data()[sample][pepcount])
                    elif self.info['Normalise'] == 'none':
                        entry['IntNorm'] = entry['Intensity']
                        entry['PepNorm'] = entry[pepcount]
                    else: raise ValueError
                except: self.errorLog('Problem processing entry %s' % (maindb.makeKey(entry)))
            maindb.saveToFile('%s.norm.tdt' % self.info['Basefile'])
        except: self.errorLog('Major error with APHID.normalise()')
#########################################################################################################################
    def combine(self,normtable='Main',comtable='Combined',save=True):   ### Combines replicates into single value per Treatment-NRID combo  |2.0|
        '''
        Combines replicates into single value per Treatment-NRID combo.
        >> normtable:str ['Main'] = Name of Table containing normalised data for combining.
        >> comtable:str ['Combined'] = Name of Table containing combined replicate data.
        >> save:bool [True] = Whether to save combined data table 
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            intensity = 'Intensity'; pepcount = self.info['PepCount']
            maindb = self.db().getTable(normtable)
            combo = '#%s#%s#%s#' % (self.info['Treatment'],maindb.info['Delimit'],'NRID')
            combo = maindb.index(combo,make=True)
            combdb = self.db().copyTable(maindb,comtable)
            combdb.makeField(fieldname='Rep')
            rules = {self.info['Identifier']:'list'}
            for field in [intensity, 'IntNorm', pepcount, 'PepNorm']: rules[field] = self.info['Combine'].lower()
            newfields = [self.info['Treatment'],'NRID','Rep', intensity, 'IntNorm', pepcount, 'PepNorm',
                         self.info['Identifier'], 'HGNC','EnsG','Desc']
            for field in combdb.fields():
                if field not in newfields: combdb.deleteField(field)
            combdb.list['Fields'] = newfields
            ### ~ [1] Compress data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            combdb.compress([self.info['Treatment'],'NRID'],rules)
            combdb.newKey([self.info['Treatment'],'NRID'])
            for key in combdb.datakeys():
                try: combdb.data()[key]['Rep'] = len(combo[key])
                except: self.errorLog('Combine problem with Reps for "%s"' % key)
            if save: combdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],comtable.lower()))
            return combdb
        except: self.errorLog('Major error with APHID.combine()')
#########################################################################################################################
    ### <5> ### Enrichment Methods                                                                                      #
#########################################################################################################################
    def enrichment(self,comtable='Combined',enrtable='Enrichment',save=True):   ### Pairwise enrichment values      |2.0|
        '''
        Calculates pairwise enrichment values.
        >> comtable:str ['Combined'] = Name of Table containing combined replicate data.
        >> enrtable:str ['Enrichment'] = Name of Table containing enrichment data.
        >> save:bool = Whether to save Enrichment table [True]
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Enrichment Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            treatment = self.info['Treatment']
            maindb = self.db().getTable('Main')
            combdb = self.db().getTable(comtable)
            self.deBug(combdb.datakeys())
            treatments = combdb.index(treatment)
            nrcomb = combdb.index('NRID')
            enrdb = self.db().copyTable(combdb,enrtable)
            enrdb.deleteField('Intensity'); enrdb.deleteField('rI')
            enrdb.reshapeWide(treatment,reshape=['Rep','IntNorm','PepNorm'])
            enrdb.newKey('NRID')
            ## ~ [0b] Treatment combos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            enrpairs = []
            for t1 in treatments:
                for t2 in treatments:
                    if self.list['EnrPairs'] and '%s:%s' % (t1,t2) not in self.list['EnrPairs']: continue
                    enr = '%s.v.%s' % (t1,t2)
                    enrdb.makeField(fieldname=enr)
                    enrpairs.append(enr)
            ### ~ [1] Go through each gene in turn calculating enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nx = 0.0; ntot = len(nrcomb)
            self.printLog('\r#NR','%s NR ID with *any* combined data.' % rje.integerString(ntot))
            for nrid in rje.sortKeys(nrcomb):
                self.progLog('\r#ENR','Calculating enrichment: %.2f%%' % (nx/ntot)); nx += 100.0
                enrdict = {}
                for enr in enrpairs:
                    (t1,t2) = string.split(enr,'.v.') 
                    c1 = '%s%s%s' % (t1,maindb.info['Delimit'],nrid)
                    c2 = '%s%s%s' % (t2,maindb.info['Delimit'],nrid)
                    #c1 = '#%s#%s#%s#' % (t1,maindb.info['Delimit'],nrid)
                    #c2 = '#%s#%s#%s#' % (t2,maindb.info['Delimit'],nrid)
                    if c1 not in combdb.datakeys() and c2 not in combdb.datakeys(): continue
                    if c1 in combdb.datakeys(): enrdict[enr] = float(combdb.data()[c1]['IntNorm'])
                    else: enrdict[enr] = self.stat['Absence']
                    if c2 in combdb.datakeys(): enrdict[enr] /= float(combdb.data()[c2]['IntNorm'])
                    else: enrdict[enr] /= self.stat['Absence']
                    enrdb.data()[nrid][enr] = enrdict[enr]
                try:
                    if not enrdict:
                        if not self.opt['Blanks']: enrdb.data().pop(nrid)
                        ntot -= 1; nx -= 100.0
                except:
                    #self.errorLog('Enrichment problem for %s' % nrid)
                    self.bugPrint(enrpairs); self.bugPrint(nrid)
                    try: self.deBug(enrdb.data()[nrid])
                    except: self.deBug(nrid in enrdb.data())
            self.printLog('\r#ENR','Calculated enrichment for %s NR ID.' % rje.integerString(ntot))
            ### ~ [2] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: enrdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],enrtable.lower()))
            return enrdb
        except: self.errorLog('Major error with APHID.enrichment()')
#########################################################################################################################
    ### <6> ### Jack-Knifing  Methods                                                                                   #
#########################################################################################################################
    def jackKnife(self):    ### Assess Pairwise enrichment values using (possible) jack-knife approach              |2.0|
        '''Assess Pairwise enrichment values using (possible) jack-knife approach.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['JackKnife']: return
            treatment = self.info['Treatment']; replicate = self.info['Replicate']
            maindb = self.db().getTable('Main')
            enrdb = self.db().getTable('Enrichment')
            treatments = self.db().getTable('Combined').index(treatment)
            ## ~ [0a] Enrichment Pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            enrpairs = []
            for t1 in treatments:
                for t2 in treatments:
                    if self.list['EnrPairs'] and '%s:%s' % (t1,t2) not in self.list['EnrPairs']: continue
                    enr = '%s.v.%s' % (t1,t2)
                    enrpairs.append(enr)
            ## ~ [0b] Pairwise Replicate combinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            samples = {}
            sampcomb = '#%s#%s#%s#' % (self.info['Treatment'],maindb.info['Delimit'],self.info['Replicate'])
            for sample in maindb.index(sampcomb,make=True):
                (treat,rep) = string.split(sample,maindb.info['Delimit'])
                if treat not in samples: samples[treat] = {}
                samples[treat][rep] = sample
            ### ~ [1] Perform jack-knifing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for enr in enrpairs:
                (t1,t2) = string.split(enr,'.v.')
                ctot = len(samples[t1]) * len(samples[t2]); cx = 0
                ## ~ [1a] Generate new data tables with reduced data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                jackenr = {}     # Dictionaries of tables
                for r1 in rje.sortKeys(samples[t1]):
                    for r2 in rje.sortKeys(samples[t2]):
                        cx += 1; self.printLog('\r#JACK','Jacknifing %s %d of %d' % (enr,cx,ctot),log=False)
                        j = '%s-%s.%s-%s' % (t1,r1,t2,r2)
                        jdb = self.db().copyTable(maindb,j)
                        jx = jdb.entryNum()
                        #jdb.dropEntries(['%s==%s' % (sampcomb,samples[t1][r1]),'%s==%s' % (sampcomb,samples[t2][r2])])
                        for jkey in jdb.datakeys():
                            entry = jdb.data()[jkey]
                            if entry[treatment] == t1 and entry[replicate] == r1: jdb.data().pop(jkey)
                            elif entry[treatment] == t2 and entry[replicate] == r2: jdb.data().pop(jkey)
                            elif entry[treatment] not in [t1,t2]: jdb.data().pop(jkey)
                        self.printLog('#DROP','Reduced %s entries to %s for Jackknife %s' % (rje.integerString(jx),rje.integerString(jdb.entryNum()),j))
                        #for t3 in treatments:
                        #    if t3 not in [t1,t2]: jdb.dropEntries(['%s==%s' % (treatment,t3)])
                        jcomb = self.combine(normtable=j,comtable='%s.Combined' % j,save=self.opt['Test'])
                        jackenr[j] = self.enrichment(comtable='%s.Combined' % j,enrtable='%s.Enrichment' % j,save=self.opt['Test'])
                        self.db().list['Tables'].remove(jdb); self.db().list['Tables'].remove(jcomb)
                ## ~ [1b] Assess observed full-data enrichment with jackknives ~~~~~~~~~~~~~~~~~~~~ ##
                if not jackenr: self.printLog('#JACK','No jackknife combos for %s' % enr); continue
                penr = 'pj.%s' % enr; enrdb.makeField(fieldname=penr)
                ex = 0.0; etot = enrdb.entryNum(); eex = 0
                for entry in enrdb.entries():
                    self.progLog('\r#JACK','Jackknifing %s > %s: %.2f%%' % (enr,self.stat['EnrCut'],(ex/etot))); ex += 100.0
                    if not entry[enr]: entry[penr] = 'n/a'; continue
                    #x#elif float(entry[enr]) <= self.stat['EnrCut']: entry[penr] = '-'; continue
                    else: entry[penr] = 0.0; eex += 1
                    ekey = enrdb.makeKey(entry)
                    for jdb in jackenr.values():
                        try: jentry = jdb.data()[ekey]
                        except: continue    #self.errorLog('Jackknife DB error'); 
                        if jentry[enr] and float(jentry[enr]) > self.stat['EnrCut']: entry[penr] += 1
                    entry[penr] /= len(jackenr)
                self.printLog('\r#JACK','Jackknifing %s > %s: %s entries' % (enr,self.stat['EnrCut'],rje.integerString(eex)))
            ### ~ [2] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            enrdb.saveToFile('%s.jackknife.tdt' % self.info['Basefile'])
        except: self.errorLog('Major error with APHID.enrichment()')
#########################################################################################################################
    def bootstrap(self):    ### Assess Pairwise enrichment values using (possible?) bootstrap approach              |2.0|
        '''Assess Pairwise enrichment values using (possible?) bootstrap approach.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['JackKnife']: return
            treatment = self.info['Treatment']; replicate = self.info['Replicate']
            maindb = self.db().getTable('Main')
            enrdb = self.db().getTable('Enrichment')
            treatments = self.db().getTable('Combined').index(treatment)
            ## ~ [0a] Enrichment Pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            enrpairs = []
            for t1 in treatments:
                for t2 in treatments:
                    if self.list['EnrPairs'] and '%s:%s' % (t1,t2) not in self.list['EnrPairs']: continue
                    enr = '%s.v.%s' % (t1,t2)
                    enrpairs.append(enr)
            ## ~ [0b] Pairwise Replicate combinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            samples = {}
            sampcomb = '#%s#%s#%s#' % (self.info['Treatment'],maindb.info['Delimit'],self.info['Replicate'])
            for sample in maindb.index(sampcomb,make=True):
                (treat,rep) = string.split(sample,maindb.info['Delimit'])
                if treat not in samples: samples[treat] = {}
                samples[treat][rep] = sample
            ### ~ [1] Perform bootstrapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for enr in enrpairs:
                (t1,t2) = string.split(enr,'.v.')
                ctot = len(samples[t1]) * len(samples[t2]); cx = 0
                ## ~ [1a] Generate new data tables with reduced data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                jackenr = {}     # Dictionaries of tables
                for r1 in rje.sortKeys(samples[t1]):
                    for r2 in rje.sortKeys(samples[t2]):
                        cx += 1; self.printLog('\r#BOOT','Bootstrapping %s %d of %d' % (enr,cx,ctot),log=False)
                        j = '%s-%s.%s-%s' % (t1,r1,t2,r2)
                        jdb = self.db().copyTable(maindb,j)
                        jx = jdb.entryNum()
                        #jdb.dropEntries(['%s==%s' % (sampcomb,samples[t1][r1]),'%s==%s' % (sampcomb,samples[t2][r2])])
                        for jkey in jdb.datakeys():
                            entry = jdb.data()[jkey]
                            if entry[treatment] == t1 and entry[replicate] == r1: pass      #jdb.data().pop(jkey)
                            elif entry[treatment] == t2 and entry[replicate] == r2: pass    #jdb.data().pop(jkey)
                            else: jdb.data().pop(jkey)
                            #elif entry[treatment] not in [t1,t2]: jdb.data().pop(jkey)
                        self.printLog('#DROP','Reduced %s entries to %s for Bootstrap %s' % (rje.integerString(jx),rje.integerString(jdb.entryNum()),j))
                        #for t3 in treatments:
                        #    if t3 not in [t1,t2]: jdb.dropEntries(['%s==%s' % (treatment,t3)])
                        jcomb = self.combine(normtable=j,comtable='%s.Combined' % j,save=self.opt['Test'])
                        jackenr[j] = self.enrichment(comtable='%s.Combined' % j,enrtable='%s.Enrichment' % j,save=self.opt['Test'])
                        self.db().list['Tables'].remove(jdb); self.db().list['Tables'].remove(jcomb)
                ## ~ [1b] Assess observed full-data enrichment with bootstraps ~~~~~~~~~~~~~~~~~~~~ ##
                if not jackenr: self.printLog('#BOOT','No bootstrap combos for %s' % enr); continue
                penr = 'pi.%s' % enr; enrdb.makeField(fieldname=penr)
                ex = 0.0; etot = enrdb.entryNum(); eex = 0
                for entry in enrdb.entries():
                    self.progLog('\r#BOOT','Bootstrapping %s > %s: %.2f%%' % (enr,self.stat['EnrCut'],(ex/etot))); ex += 100.0
                    if not entry[enr]: entry[penr] = 'n/a'; continue
                    #x#elif float(entry[enr]) <= self.stat['EnrCut']: entry[penr] = '-'; continue
                    else: entry[penr] = 0.0; eex += 1
                    ekey = enrdb.makeKey(entry)
                    for jdb in jackenr.values():
                        try: jentry = jdb.data()[ekey]
                        except: continue    #self.errorLog('Jackknife DB error'); 
                        if jentry[enr] and float(jentry[enr]) > self.stat['EnrCut']: entry[penr] += 1
                    entry[penr] /= len(jackenr)
                self.printLog('\r#BOOT','Bootknifing %s > %s: %s entries' % (enr,self.stat['EnrCut'],rje.integerString(eex)))
            ### ~ [2] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            enrdb.saveToFile('%s.bootstrap.tdt' % self.info['Basefile'])
        except: self.errorLog('Major error with APHID.bootstrap()')
#########################################################################################################################
    ### <X> ### OLD Methods                                                                                             #
#########################################################################################################################
    def makeFamilies(self): ### Uses GABLAM to organise genes into families based on sequence identity.
        '''Calculates sample-specific enrichments given normalised intensities.'''
        try:### ~ [1] Setup objects etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pingu = self.obj['PINGU']
            seqfile = pingu.info['ResDir'] + '%s.Combined.fas' % self.info['Basefile']
            famgablam = gablam.GABLAM(self.log,self.cmd_list)
            famgablam.opt['QryAcc'] = False
            gfile = famgablam.info['GablamOut'] = pingu.info['ResDir'] + '%s.gablam.tdt' % self.info['Basefile']
            famgablam.info['HitSumOut'] = pingu.info['ResDir'] + '%s.hitsum.tdt' % self.info['Basefile']
            famgablam.info['QueryDB'] = famgablam.info['SearchDB'] = seqfile
            slimcore = rje_slimcore.SLiMCore(self.log,self.cmd_list)
            slimcore.setOpt({'EFilter':True})
            slimcore.setInfo({'GablamDis':gfile})

            ### ~ [2] Generate Families ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Use GABLAM and SLiMCore to generate UPC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(gfile) or (self.interactive() >= 2 and not rje.yesNo('Use existing %s?' % gfile)):
                famgablam.gablam()
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','seqin=%s' % seqfile]
            slimcore.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            slimcore.setupBasefile()
            slimcore.setNum({'StartTime':time.time()})
            slimcore.makeUPC()
            ## ~ [2b] Map genes onto UPC for conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upgenes = {}
            seqdict = slimcore.obj['SeqList'].seqNameDic()
            ensloci = {}    # Mapping of genes to sequence short names
            ensname = {}    # Mapping of sequence short names to gene names
            for seq in slimcore.obj['SeqList'].seq:
                gene = rje.matchExp(' gene:(\S+)\]',seq.info['Description'])[0]
                ensloci[gene] = seq.shortName()
                upc = slimcore.getUP(seq)
                if upc not in upgenes: upgenes[upc] = []
                upgenes[upc].append(gene)
                symbol = pingu.geneMap(gene)
                ensname[seq.shortName()] = symbol
            ## ~ [2c] Convert UPC to family dictionary for each gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Family'] = {}
            for upc in upgenes:
                symblist = []
                for gene in upgenes[upc]:
                    symbol = pingu.geneMap(gene)
                    if symbol not in symblist: symblist.append(symbol)
                symblist.sort()
                for i in range(len(symblist)):
                    if symblist[0][:4] == 'ENSG': symblist.append(symblist.pop(0))
                    else: break
                fx = len(upgenes[upc]) - 1
                fam = '%s+%d' % (symblist[0],fx)
                for gene in upgenes[upc]: self.dict['Family'][gene] = fam
            #x#self.deBug(self.dict['Family'])
    
            ### ~ [3] Generate Cytoscape network of families ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            matrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            matrix.loadFromDataTable(gfile,normalise=100.0,inverse=True,checksym=False)  ### Loads matrix from GABLAM
            matrix.forceSymmetry()  ### Use min distance for each pair of values
            matrix.rename(ensname)  ### Goes through matrix and renames objects using given dictionary
            matrix.saveCytoscape(self.info['Basefile'],type='gablam',cutoff=1.0,inverse=False)  ### Output for Cytoscape

            return True
        except:
            self.log.errorLog('APHID.makeFamilies() failure')
            raise
        return False
#########################################################################################################################
    ### <6> ### Visualisation Methods                                                                                   #
#########################################################################################################################
    def visualisations(self):   ### Generate Enrinchment visualisations
        '''Generate Enrinchment visualisations.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pingu = self.obj['PINGU']
            pingu.genericData()
            self.printLog('#PNG','Visualise %s.*.enrichment.tdt files' % self.info['Absolute Base'])
            tdtfiles = rje.getFileList(self,filelist=['%s.*.enrichment.tdt' % self.info['Absolute Base']],subfolders=False,summary=False)
            tdtfiles.sort()
            for efile in tdtfiles: self.enrichPNG(efile)
        except: self.errorLog('APHID.visualisations() failure')
#########################################################################################################################
    def enrichPNG(self,efile):  ### Generate visualisations for one enrichment file
        '''Generate visualisations for one enrichment file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pingu = self.obj['PINGU']
            edata = rje.dataDict(self,efile,['Identifier'],getheaders=True)
            samples = edata.pop('Headers')[4:]
            print efile
            print '>>>', samples
            newfile = string.replace(efile,'enrichment','enr_pingu')
            newhead = ['Sample','Identifier','Score']
            rje.delimitedFileOutput(self,newfile,newhead)
            for e in rje.sortKeys(edata):
                for sample in samples:
                    if edata[e][sample]:
                        try: score = string.atof(edata[e][sample])
                        except:
                            self.errorLog('%s: %s' % (sample, edata[e][sample]))
                            score = -1
                        if sample.find('~v~') > 0 and score <= self.stat['EnrCut']: continue
                        rje.delimitedFileOutput(self,newfile,newhead,datadict={'Sample':sample,'Identifier':e,'Score':edata[e][sample]})
            self.printLog('#TDT','Enrichment data output to %s for Pingu analysis' % newfile)
            self.obj['PINGU'].info['ResDir'] = rje.makePath(rje.baseFile(newfile))
            self.obj['PINGU'].list['Data'] = [newfile]
            self.obj['PINGU'].run()

            ### ~ [2] ~ Visualisations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Interactome for each sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for sample in pingu.dict['Datasets']: pingu.interactomePNG(sample)

            ## ~ [2b] ~ GO enrichment for each sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            basefile = '%s%s' % (self.obj['PINGU'].info['ResDir'],rje.baseFile(newfile,strip_path=True))
            basefile = string.join(string.split(basefile,'.')[:-1],'.')
            for g in ['go','goslim_power','goslim_generic']:
                rcmd = '%s --no-restore --no-save --args "%s" "%s"' % (self.info['RPath'],g,basefile)
                rslimjim = '%srje_call.r' % self.info['Path']
                rcmd += ' < "%s" > "%s.r.tmp.txt" 2>&1' % (rslimjim,basefile)
                self.printLog('#RCALL',rcmd)
                problems = os.popen(rcmd).read()
                if problems: self.errorLog(problems,printerror=False)
                #if not os.path.exists('%s.png' % basefile): self.errorLog('%s.png not created' % basefile,printerror=False)
                if os.path.exists('%s.png' % basefile) and not self.opt['Test'] and not False: 
                    for ext in ['r.tmp.txt']:
                        if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))

            ## ~ [2c] ~ Venn comparisons for samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # ... Combos of original samples
            # ... Enrichments vs other samples
        except: self.errorLog('APHID.enrichPNG(%s) failure' % efile)
#########################################################################################################################
### End of SECTION II: APHID Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
class GeneMap(rje_genemap.GeneMap):
    '''Just for pickling.'''
#########################################################################################################################
class PINGU(pingu.PINGU):
    '''Just for pickling.'''
#########################################################################################################################
def eString(value):  ### Returns formatted string for _expect value
    '''Returns formatted string for _expect value.'''
    try:
        if value >= 1000: return '%.2e' % value
        if value >= 10: return '%.1f' % value
        elif value >= 0.1: return '%.2f' % value
        elif value >= 0.01: return '%.3f' % value
        else: return '%.2e' % value
    except:
        print expect
        raise
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
    try: APHID(mainlog,cmd_list).run()
        
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
