#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_taxmap
Description:  Mapping taxonomic groups onto sequences based on Newick trees.
Version:      0.4.0
Last Edit:    06/06/16
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    ### ~ Input Options ~ ###
    nwklist=FILES       : List of Newick format files. (*.nwk names must match accnum.) Wildcards allowed. [*.nwk]
    protdesc=FILE       : Delimited file of protein, description [None]
    classify=FILELIST   : List of files containing protein IDs (first column) for additional taxsum outputs [*.class]
    basefile=X          : Base for output files. Will reuse unless force=T ['taxmap']
    taxbase=X           : Base of previous run to load results from, if not basefile []

    ### ~ Taxon Assignment Options ~ ###
    minboot=X       : Minimum bootstrap value for an "in-clade" with query protein [0.5]
    noneboot=X      : Bootstrap support to give "None" ratings [1.0]
    minscore=X      : Filter out species codes with score < minscore [1.0]
    minsum=X        : Convert taxa with score < minsum to "Other" for *.taxsum.tdt output [10.0]
    minclass=X      : Convert taxa with score < minclass to "Other" for *.CLASS.tdt output [1.0]
    bootweight=T/F  : Whether to weight taxon scores for each protein by clade bootstrap support for filtering [True]
    bootfilter=X    : Filter bootstrap support < X to "Uncertain" [0.0]
    monophyly=T/F   : Enforce strict monophyly and change any multiple assignments to "Uncertain" [False]
    taxfilter=LIST  : List of taxonomic assignments to filter from rank/summary (e.g. "Uncertain") []

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_tree, rje_taxonomy
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial Working Version.
    # 0.2.0 - Added bootstraps and minscore.
    # 0.2.1 - Fixed missing species names in *.taxa.tdt output.
    # 0.3.0 - Added monophyly=T/F   : Enforce strict monophyly and change any multiple assignments to "Uncertain" [False]
    # 0.4.0 - Added TaxFilter, minsum=X and classify=FILELIST files containing protein IDs (first column) for additional taxsum outputs [*.class]
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add parsing of gene names.
    # [Y] : Add parsing and processing of bootstraps.
    # [Y] : Add filtering of low-abundance taxa. (minscore=X)
    # [ ] : Add inference/propagation of high-abundance taxa. (WinnerTakesAll)
    # [Y] : Calculate both weighted and unweighted. (Use option for filter) Move noneboot to TaxaMap/Filter (replace).
    # [Y] : Add meanboot to taxasum.
    # [Y] : Add table of taxonomic ranks to output.
    # [Y] : Add taxa to screen from taxasum table and calculation.
    # [X] : Add split output for hypothetical (and NOT) proteins only.
    # [X] : Add function to restrict analysis to a set of gene descriptions only.
    # [Y] : classify=FILELIST files containing protein IDs (first column) for additional taxsum outputs [*.class]
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('TAXMAP', '0.4.0', 'June 2016', '2016')
    description = 'Mapping taxonomic groups onto sequences based on Newick trees'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
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
### SECTION II: TaxMap Class                                                                                            #
#########################################################################################################################
class TaxMap(rje_obj.RJE_Object):
    '''
    Class. Author: Rich Edwards (2015).

    Str:str
    - ProtDesc=FILE   : Delimited file of protein, description [None]

    Bool:boolean
    - BootWeight=T/F  : Whether to weight taxon scores for each protein by clade bootstrap support [True]
    - Monophyly=T/F   : Enforce strict monophyly and change any multiple assignments to "Uncertain" [False]
    - TaxBase=X           : Base of previous run to load results from, if not basefile []

    Int:integer

    Num:float
    - BootFilter=X    : Filter bootstrap support < X to "Uncertain" [0.0]
    - MinBoot=X       : Minimum bootstrap value for an "in-clade" with query protein [0.5]
    - MinClass=X      : Convert taxa with score < minclass to "Other" for *.CLASS.tdt output [1.0]
    - MinScore=X      : Filter out taxonomic groups with score < minscore [1.0]
    - MinSum=X        : Convert taxa with score < sumscore to "Other" for *.taxsum.tdt output [10.0]
    - NoneBoot=X      : Bootstrap support to give "None" ratings [1.0]

    File:file handles with matching str filenames
    
    List:list
    - Classify=FILELIST   : List of files containing protein IDs (first column) for additional taxsum outputs [*.class]
    - NwkList=FILES   : List of Newick format files. (*.nwk names must match accnum.) Wildcards allowed. [*.nwk]
    - TaxFilter=LIST  : List of taxonomic assignments to filter from rank/summary (e.g. "Uncertain") []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['ProtDesc','TaxBase']
        self.boollist = ['BootWeight','Monophyly']
        self.intlist = []
        self.numlist = ['BootFilter','MinBoot','MinScore','MinClass','MinSum','NoneBoot']
        self.filelist = []
        self.listlist = ['Classify','NwkList','TaxFilter']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'BootWeight':True,'Monophyly':False})
        self.setInt({})
        self.setNum({'BootFilter':0.0,'MinBoot':0.5,'MinClass':1.0,'MinScore':1.0,'MinSum':10.0,'NoneBoot':1.0})
        self.obj['DB'] = rje_db.Database(self.log,['tuplekeys=T','basefile=taxmap']+self.cmd_list)
        self.baseFile(self.obj['DB'].baseFile())
        self.list['Classify'] = glob.glob('*.class') # List of files using wildcards and glob
        self.list['NwkList'] = glob.glob('*.nwk') # List of files using wildcards and glob
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
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['TaxBase'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['ProtDesc'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BootWeight','Monophyly'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                self._cmdReadList(cmd,'float',['BootFilter','MinBoot','MinClass','MinScore','MinSum','NoneBoot']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['TaxFilter'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Classify','NwkList']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not (self.db('taxa') and self.db('taxamap')):
                self.treeListSPCode()
                self.filterSPCode()
                self.taxaMap()
            self.summaryScores()
            #!# Idenitfy duplicates with non-overlapping taxa (quantify overlap at each rank!)
            self.classify()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Protein descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['ProtDesc'] = {}
            if self.getStrLC('ProtDesc'):
                for fline in open(self.getStr('ProtDesc'),'r').readlines():
                    [prot,desc] = string.split(rje.chomp(fline),maxsplit=1)
                    self.dict['ProtDesc'][prot] = desc
                #self.db().addTable(self.getStr('ProtDesc'),mainkeys=['protein'],datakeys='All',headers=['protein','description'],ignore=['#'],name='protdesc',expect=True)
            ## ~ [0b] Look for previous run results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            taxdb = self.db('taxa',add=True,forcecheck=True,mainkeys=['spcode'])
            if not taxdb and self.getStrLC('TaxBase') and not self.force():
                spfile = '%s.taxa.tdt' % self.getStr('TaxBase')
                taxdb = db.addTable(spfile,mainkeys=['spcode'],name='taxa',expect=False)
            mapdb = self.db('taxamap',add=True,forcecheck=True,mainkeys=['protein'])
            if not mapdb and self.getStrLC('TaxBase') and not self.force():
                spfile = '%s.taxamap.tdt' % self.getStr('TaxBase')
                mapdb = db.addTable(spfile,mainkeys=['protein'],name='taxamap',expect=False)
            if taxdb and mapdb:
                taxdb.dataFormat({'boot':'num'})
                mapdb.dataFormat({'boot':'num'})
                return True
            ## ~ [0c] Taxonomy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['Taxonomy'] = rje_taxonomy.Taxonomy(self.log,self.cmd_list)
            self.obj['Taxonomy'].setup(force=False)
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def treeListSPCode(self):  ### Main taxa mapping from list of tree files
        '''Main taxa mapping from list of tree files.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            specdb = self.db('spcode',add=True,forcecheck=True,mainkeys=['protein'])
            if not specdb and self.getStrLC('TaxBase') and not self.force():
                spfile = '%s.spcode.tdt' % self.getStr('TaxBase')
                specdb = db.addTable(spfile,mainkeys=['protein'],name='spcode',expect=False)
            if specdb: specdb.dataFormat({'boot':'num'}); return True
            specdb = db.addEmptyTable('spcode',['protein','boot','spcode','inpara','paralogues'],['protein'])
            #dupdb = db.addEmptyTable('para',['protein','paralogues'],['protein'])
            self.dict['Duplicates'] = {}    # {prot1:[dups]}
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for nwkfile in self.list['NwkList']:
                tree = rje_tree.Tree(self.log,self.cmd_list)
                tree.loadTree(nwkfile,seqlist=None,postprocess=False)
                seqacc = rje.baseFile(nwkfile,strip_path=True)
                # Identify node corresponding to query sequence
                seqnode = None
                for node in tree.nodes():
                    try:
                        if string.split(node.shortName(),'__')[1] == seqacc: seqnode = node
                    except: pass    # Internal node or bad sequence format
                if not seqnode:
                    self.warnLog('Could not find %s in %s nodes!' % (seqacc,nwkfile))
                    continue
                # Get species code for query sequence
                seqspec = tree.cladeSpec(seqnode)
                if len(seqspec) != 1: self.warnLog('Could not find species in %s node!' % (seqacc)); continue
                seqspec = seqspec.keys()[0]
                if seqspec != string.split(seqnode.shortName(),'_')[1]: raise ValueError('Species mismatch for %s & %s' % (seqacc,seqnode.shortName()))
                # Find ancestor with closest orthologue outgroup
                rootnode = tree._getRootNode()
                if not rootnode: self.warnLog('Could not find root node in %s!' % (nwkfile)); continue
                ancnode = seqnode.ancNode()
                try: bootx = float(ancnode.ancBranch().stat['Bootstrap'])/tree.stat['Bootstraps']
                except: bootx = 1.0
                inparanode = None    # Node to define in-paralogues
                ancspec = tree.cladeSpec(ancnode)
                while len(ancspec) < 2 or bootx < self.getNum('MinBoot'):
                    inparanode = ancnode    # All same species
                    if ancnode == rootnode: break
                    ancnode = ancnode.ancNode(); ancspec = tree.cladeSpec(ancnode)
                    try: bootx = float(ancnode.ancBranch().stat['Bootstrap'])/tree.stat['Bootstraps']
                    except: bootx = 1.0
                ancspec.pop(seqspec)    # Now only have counts of closest other species
                # Update table, replacing species codes with genera?
                sentry = {'protein':seqacc,'spcode':rje.sortUnique(ancspec.keys())}
                sentry['boot'] = bootx
                if not ancspec: sentry['spcode'] = ['None']; sentry['boot'] = self.getNum('NoneBoot')
                sentry['spcode'] = string.join(sentry['spcode'],'|')
                # Establish list of duplicate proteins
                inpara = []     # List of in-paralogue nodes
                inparacc = []   # List of in-paralogue accnum
                if inparanode: inpara = tree._nodeClade(inparanode,internal=False)
                self.dict['Duplicates'][seqacc] = []
                for node in tree._nodeClade(rootnode,internal=False):
                    if node == seqnode: continue
                    if len(string.split(node.shortName(),'_')) < 2: continue
                    if string.split(node.shortName(),'_')[1] == seqspec:
                        paracc = string.split(node.shortName(),'__')[1]
                        if node in inpara: inparacc.append(paracc)
                        else: self.dict['Duplicates'][seqacc].append(paracc)
                sentry['inpara'] = string.join(inparacc,'|')
                sentry['paralogues'] = string.join(self.dict['Duplicates'][seqacc],'|')
                specdb.addEntry(sentry)
            ## Update specdb and save
            specdb.saveToFile()
            #dupdb.saveToFile()
            return True
        except:
            self.errorLog(self.zen())
            return False
#########################################################################################################################
    def taxaMap(self):      ### Maps species codes onto different taxonomic ranks.
        '''Maps species codes onto different taxonomic ranks.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            tax = self.obj['Taxonomy']
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            specdb = self.db('spcode')
            #descdb = self.db('protdesc')
            ranks = ['genus','family','order','class','phylum']
            rankmap = {}    # SPCODE to Taxon dictionary
            rankfields = ['protein']+ranks+specdb.fields()[1:]
            #if descdb: rankfields.append('desc')
            if self.getStrLC('ProtDesc'):
                rankfields.append('desc'); px = 0
                for prot in self.dict['ProtDesc']:
                    if prot.lower() in ['','protein','gene']: continue
                    pentry = {'protein':prot,'spcode':'None','boot':self.getNum('NoneBoot')}
                    pkey = specdb.makeKey(pentry)
                    if pkey not in specdb.dataKeys(): specdb.addEntry(pentry); px += 1
                self.printLog('#PROT','Added %s proteins from %s without trees.' % (rje.iStr(px),self.getStr('ProtDesc')))
            rankdb = db.addEmptyTable('taxamap',rankfields,['protein'])
            for rank in ranks: rankmap[rank] = {'None':'None','Unmapped':'Unmapped','Uncertain':'Uncertain'}
            taxdb = db.addEmptyTable('taxa',['spcode','taxid','name']+ranks,['spcode'])

            sx = 0.0; stot = specdb.entryNum()
            for entry in specdb.entries():
                self.progLog('\r#SPEC','Processing species: %.2f%%' % (sx/stot)); sx += 100.0
                #if descdb:
                    #try: entry['desc'] = descdb.data(descdb.makeKey(entry))['description']
                try: entry['desc'] = self.dict['ProtDesc'][entry['protein']]
                except: entry['desc'] = ''
                for spcode in string.split(entry['spcode'],'|'):
                    if spcode in rankmap['genus']: continue
                    tentry = {'spcode':spcode}
                    try:
                        taxid = tax.mapToTaxID(spcode,nodeonly=True,warn=False)[0]
                        rank = tax.dict['Rank'][taxid]
                        tentry['taxid'] = taxid
                        tentry['name'] = tax.getSpecies(taxid)
                    except:
                        self.warnLog('Unable to map species code "%s" to TaxID -> "Unmapped"' % spcode)
                        taxid = 'Unmapped'
                        rank = 'genus'
                    # Loop through different ranks
                    for ri in range(len(ranks)):
                        nextrank = ranks[ri]
                        while rank not in ranks[ri:] and taxid in tax.dict['Parent']:
                            taxid = tax.dict['Parent'][taxid]
                            rank = tax.dict['Rank'][taxid]
                            #self.debug('%s: %s' % (tax.dict['Rank'][taxid],tax.getSpecies(taxid)))
                        if taxid in tax.dict['Parent']: taxon = tax.getSpecies(taxid)
                        else: taxon = 'Unmapped'
                        if rank != nextrank:
                            if self.getBool('Monophyly'): taxon = 'Uncertain'
                            else: taxon = '%s %s.' % (taxon,nextrank[:3])
                        rankmap[nextrank][spcode] = taxon
                        tentry[nextrank] = taxon
                    taxdb.addEntry(tentry)
                rentry = {}
                for nextrank in ranks:
                    taxa = []
                    for spcode in string.split(entry['spcode'],'|'): taxa.append(rankmap[nextrank][spcode])
                    if len(taxa) > 1 and 'None' in taxa:
                        self.warnLog('None in: %s' % string.join(rje.sortUnique(taxa),'|'))
                        taxa.remove('None')
                    if len(taxa) > 1 and 'Unmapped' in taxa: taxa.remove('Unmapped')
                    if len(taxa) > 1 and self.getBool('Monophyly'): rentry[nextrank] = 'Uncertain'
                    else: rentry[nextrank] = string.join(rje.sortUnique(taxa),'|')
                rankdb.addEntry(rje.combineDict(rentry,entry))
            self.printLog('\r#SPEC','%s proteins with species codes processed.' % rje.iStr(stot))
            rankdb.saveToFile()
            taxdb.saveToFile()
        except: self.errorLog('%s.taxaMap error' % self.prog())
#########################################################################################################################
    def summaryScores(self,rankdb=None,sumstr='taxasum',minsum='MinSum'):   ### Generates summary scores from rank table.
        '''Generates summary scores from rank table.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not rankdb: rankdb = self.db('taxamap')
            sumdb = db.addEmptyTable(sumstr,['rank','taxon','count','bootwt','meanboot','perc','wtperc'],['rank','taxon'])
            ranks = ['genus','family','order','class','phylum']
            ### ~ [2] Normalise to reduced levels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for rank in ranks:
                self.printLog('\r#RANK','Normalising %s data.' % rank)
                taxsum = {}; ranksum = 0.0  # Summed counts for taxa and rank total
                taxwt = {}; wtsum = 0.0     # Bootstrap-weighted summed counts for taxa and rank total
                bootsum = {}; bootx = {}    # Sum and count of bootstrap values for mean boot numbers
                for entry in rankdb.entries():
                    taxa = string.split(entry[rank],'|')
                    for taxon in taxa:
                        if taxon in self.list['TaxFilter']: continue
                        if taxon not in taxsum:
                            taxsum[taxon] = 0.0; taxwt[taxon] = 0.0
                            bootsum[taxon] = 0.0; bootx[taxon] = 0
                        taxsum[taxon] += 1.0 / len(taxa)
                        ranksum += 1.0 / len(taxa)
                        taxweight = entry['boot']
                        bootsum[taxon] += entry['boot']; bootx[taxon] += 1
                        taxwt[taxon] += taxweight / len(taxa)
                        wtsum += taxweight / len(taxa)
                otherx = 0
                for taxon in rje.sortKeys(taxsum):
                    if taxon == 'Other': continue
                    if taxsum[taxon] < self.getNum(minsum):
                        if 'Other' not in taxsum:
                            taxsum['Other'] = 0.0
                            taxwt['Other'] = 0.0
                            bootsum['Other'] = 0.0
                            bootx['Other'] = 0.0
                        taxsum['Other'] += taxsum.pop(taxon)
                        taxwt['Other'] += taxwt.pop(taxon)
                        bootsum['Other'] += bootsum.pop(taxon)
                        bootx['Other'] += bootx.pop(taxon)
                        otherx += 1
                self.printLog('#MINSUM','%s %s taxa converted to "Other" (count < minsum=%s)' % (rje.iStr(otherx),rank,self.getNum(minsum)))
                for taxon in taxsum: sumdb.addEntry({'rank':rank,'taxon':taxon,'count':rje.dp(taxsum[taxon],1),
                                                     'perc':rje.sf(100.0*taxsum[taxon]/ranksum),
                                                     'bootwt':rje.dp(taxwt[taxon],1),'meanboot':rje.dp(bootsum[taxon]/bootx[taxon],3),
                                                     'wtperc':rje.sf(100.0*taxwt[taxon]/wtsum)})
            ## ~ [2a] Rank taxa by counts such that highest is Rank 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sumdb.rankFieldByIndex('rank','count',rev=True,absolute=True,lowest=True)
            sumdb.rankFieldByIndex('rank','bootwt',rev=True,absolute=True,lowest=True)
            ## ~ [2b] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sumdb.saveToFile()
        except: self.errorLog('%s.summaryScores error' % self.prog())
#########################################################################################################################
    def classify(self): ### Generate summary tables for each protein class
        '''Generate summary tables for each protein class.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            rankdb = self.db('taxamap')
            for cfile in self.list['Classify']:
                pclass = rje.baseFile(cfile,strip_path=True)
                clist = []
                for fline in open(cfile,'r').readlines():
                    prot = string.split(rje.chomp(fline),maxsplit=1)[0]
                    if prot: clist.append(prot)
                self.printLog('#CLASS','%s "%s" class proteins read from %s' % (rje.iLen(clist),pclass,cfile))
                if not clist:
                    self.warnLog('No proteins read from %s' % (cfile))
                    continue
                classdb = db.copyTable(rankdb,pclass)
                classdb.dropEntriesDirect('protein',clist,inverse=True)
                if not classdb.entries():
                    self.warnLog('No "%s" proteins found in TaxaMap table' % (pclass))
                    continue
                self.summaryScores(classdb,pclass,'MinClass')
        except: self.errorLog('%s.classify() error' % self.prog())
#########################################################################################################################
    def filterSPCode(self):     ### Filters species codes according to mincount and shared taxa at different levels.
        '''Filters species codes according to mincount and shared taxa at different levels.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            tax = self.obj['Taxonomy']
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            specdb = self.db('spcode')
            parents = {}    # Dictionary of {spcode:parents}
            taxsum = {}
            # Reduced according to low abundance and/or higher level taxa of species in clade
            fx = 0; bfx = 0; ufx = 0
            for ekey in specdb.dataKeys():
                entry = specdb.data(ekey)
                if entry['spcode'] == 'None':
                    entry['boot'] = self.getNum('NoneBoot')
                    continue
                if entry['boot'] < self.getNum('BootFilter'):
                    self.printLog('#FILT','%s: filtered -> "Uncertain" (bootstrap %s < bootfilter=%s).' % (entry['protein'],entry['boot'],self.getNum('BootFilter')))
                    entry['spcode'] = 'Uncertain'; bfx += 1
                    continue
                #self.debug(entry)
                spcodes = string.split(entry['spcode'],'|')
                for spcode in spcodes[0:]:
                    if spcode not in parents:
                        parents[spcode] = []
                        try: taxid = tax.mapToTaxID(spcode,nodeonly=True,warn=False)[0]
                        except: continue
                        while taxid in tax.dict['Parent']:
                            taxid = tax.dict['Parent'][taxid]
                            parsp = tax.getSpCode(taxid,invent=False,warn=False)
                            if parsp: parents[spcode].append(parsp)
                    if not parents[spcode] and len(spcodes) > 1:
                        self.printLog('#FILT','%s: filtered unmapped spcode %s.' % (entry['protein'],spcode))
                        spcodes.remove(spcode); ufx += 1
                    for parsp in parents[spcode]:
                        if parsp in spcodes:
                            self.printLog('#FILT','%s: filtered %s as parent of %s.' % (entry['protein'],parsp,spcode))
                            spcodes.remove(parsp); fx += 1
                for taxon in spcodes[0:]:
                    if taxon not in taxsum: taxsum[taxon] = 0.0
                    if self.getBool('BootWeight'): taxweight = entry['boot']
                    else: taxweight = 1.0
                    taxsum[taxon] += taxweight / len(spcodes)
                entry['spcode'] = string.join(spcodes,'|')
            self.printLog('#FILT','Filtered %s species codes with co-occurring child taxa' % rje.iStr(fx))
            self.printLog('#FILT','Filtered %s unmapped species codes with co-occurring mapped taxa' % rje.iStr(ufx))
            if self.getNum('BootFilter') > 0.0: self.printLog('#FILT','Filtered %s proteins with bootstrap < bootfilter=%s' % (rje.iStr(bfx),self.getNum('BootFilter')))
                #self.debug(entry)

            fx = 0
            for ekey in specdb.dataKeys():
                entry = specdb.data(ekey)
                if entry['spcode'] in ['None','Uncertain']: continue
                #self.debug(entry)
                spcodes = string.split(entry['spcode'],'|')
                for spcode in spcodes[0:]:
                    if self.getNum('MinScore') > 0 and self.getNum('MinScore') > taxsum[spcode]:
                        self.printLog('#FILT','%s: filtered %s < minscore=%s.' % (entry['protein'],spcode,self.getNum('MinScore')))
                        spcodes.remove(spcode); fx += 1
                if spcodes: entry['spcode'] = string.join(spcodes,'|')
                else: self.printLog('#FILT','%s filter aborted: no spcode left!' % (entry['protein']))
                #self.debug(entry)
            self.printLog('#FILT','Filtered %s species codes failing to meet minscore=%s.' % (rje.iStr(fx),self.getNum('MinScore')))

        except: self.errorLog('%s.filterSPCode error' % self.prog())
#########################################################################################################################
    def winnerTakesAll(self):   ### Cleans up taxamap to remove other taxa when dominant taxon present.
        '''Cleans up taxamap to remove other taxa when dominant taxon present.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxdb = self.db('taxa')
            sumdb = self.db('taxsum')
            rankdb = self.db('taxamap')
            ## ~ [1a] Rank taxa by counts such that highest is Rank 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Now in summaryScores
            if self.getBool('BootWeight'): rankfield = 'bootwt'
            else: rankfield = 'count'
            sumdb.rankFieldByIndex('rank',rankfield,'ranked',rev=True,absolute=True,lowest=True)
            ### ~ [2]
        except: self.errorLog('%s.winnerTakesAll error' % self.prog())
#########################################################################################################################
### End of SECTION II: TaxMap Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

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
    try: TaxMap(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################