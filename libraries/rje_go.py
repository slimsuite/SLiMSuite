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
Module:       rje_go
Description:  Gene Ontology Parsing/Manipulation Module
Version:      1.2
Last Edit:    05/05/10
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    Parse the OBO V1.2 text download from GO into a data structure for querying by other programs, making GO Slims and
    mapping parent terms etc. Download input file from http://www.geneontology.org/GO.downloads.ontology.shtml. If no
    file is given, will attempt to download and parse http://www.geneontology.org/ontology/gene_ontology_edit.obo.

Commandline:
    obofile=FILE    : Input GO OBO V1.2 download [None]
    webobo=T/F      : Whether to download from GO website if file not given [True]
    goslim=LIST     : List of GO IDs to form basis of GO Slim []
    parentterms=LIST: Terms relating to parent relationships ['is_a','part_of']
    ensgopath=PATH  : Path to EnsGO files   (!!! Restricted to Humans Currently !!!)

Uses general modules: copy, glob, os, string, sys, time
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
    # 1.0 - Working version with RJE modules. Added EnsGO Mapping from SLiMJIM.
    # 1.1 - Updated EnsGO mapping for new BioMart structure.
    # 1.2 - And again.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add mapping of other sequences to GO terms.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_GO', '1.2', 'May 2010', '2008')
    description = 'Gene Ontology Parsing/Manipulation Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
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
### Module Attributes
gotypes = {'bp':'biological_process','cc':'cellular_component','mf':'molecular_function'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GO Class                                                                                                #
#########################################################################################################################
class GO(rje.RJE_Object):     
    '''
    GO Class. Author: Rich Edwards (2008).

    Info:str
    - EnsGOPath = Path to EnsGO files   (!!! Restricted to Humans Currently !!!)
    - OBOFile = Input GO OBO V1.2 download [None]
    
    Opt:boolean
    - WebOBO = Whether to download from GO website if file not given [True]

    Stat:numeric

    List:list
    - GOSlim = List of GO IDs to form basis of GO Slim []
    - ParentTerms = List of terms relating to parent relationships ['is_a','part_of']

    Dict:dictionary
    - AltID = Dictionary of {Alternate ID:GO ID}
    - GO = Main dictionary of {GO ID:{'name':str,'is_a':[GO IDs],'part_of':[GO_IDs],'type':bp/mf/cc}}
    - Subset = GO Slim subsets {ID:{'name':str,'terms':[GO IDs]}}

    Obj:RJE_Objects
    '''
#########################################################################################################################
    def go(self,id=None):   ### Returns specific GO data, or whole GO dictionary if no ID given, or empty if missing
        '''Returns specific GO data, or whole GO dictionary if no ID given, or empty if missing.'''
        if not id: return self.dict['GO']
        if id[:3] == 'GO:': id = id[3:]
        if id in self.dict['GO']: return self.dict['GO'][id]
        return {}
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['OBOFile','EnsGOPath']
        self.optlist = ['WebOBO']
        self.statlist = []
        self.listlist = ['GOSlim','ParentTerms']
        self.dictlist = ['AltID','GO','Subset']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.list['ParentTerms'] = ['is_a','part_of']
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
                self._cmdReadList(cmd,'file',['OBOFile'])
                self._cmdReadList(cmd,'opt',['WebOBO'])
                self._cmdReadList(cmd,'list',['GOSlim','ParentTerms'])
                self._cmdReadList(cmd,'path',['EnsGOPath'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def setup(self,read=True,children=True,map=True):   ### Sets up GO Object
        '''Sets up GO Object.'''
        if read: self.readGO()
        if children: self.makeChildren()
        if map: self.mapEnsGO()
#########################################################################################################################
    def test(self): ### Development method
        '''Development method.'''
        self.readGO()
        self.mapEnsGO()
        gohead = ['EnsG','GO_ID','GO_Type','GO_Desc']
        gofile = 'test.go.tdt'
        rje.delimitedFileOutput(self,gofile,gohead,rje_backup=True)
        gx = 0.0; gtot = len(self.dict['EnsGO'])
        for gene in rje.sortKeys(self.dict['EnsGO']):
            self.progLog('\r#ENSGO','Compiling %s: %.2f%%' % (gofile,gx/gtot)); gx += 100.0
            for goid in self.dict['EnsGO'][gene]:
                godata = {'EnsG':gene, 'GO_ID':goid}
                godata['GO_Type'] = self.dict['GO'][goid]['type']
                godata['GO_Desc'] = self.dict['GO'][goid]['name']
                rje.delimitedFileOutput(self,gofile,gohead,datadict=godata)
        self.printLog('\r#ENSGO','Compiling %s all done: %s genes.' % (gofile,rje.integerString(gtot)))
#########################################################################################################################
    ### <2> ### Loading/Parsing Methods                                                                                 #
#########################################################################################################################
    def readGO(self):   ### Reads from self.info['OBOFile'] if present, or tries web download
        '''Reads from self.info['OBOFile'] if present, or tries web download.'''
        try:### ~ [1] Try to read file lines and parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['OBOFile'].lower() not in ['','none']:
                if os.path.exists(self.info['OBOFile']): return self.parseGO(self.loadFromFile(self.info['OBOFile']))
                self.errorLog('OBO V1.2 file "%s" missing!' % self.info['OBOFile'],printerror=False)
            ### ~ [2] Try to parse from website ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#OBO','Downloading http://www.geneontology.org/ontology/gene_ontology_edit.obo ...',newline=False,log=False)
            obolines = urllib2.urlopen('http://www.geneontology.org/ontology/gene_ontology_edit.obo').readlines()
            self.printLog('\r#OBO','Downloaded http://www.geneontology.org/ontology/gene_ontology_edit.obo: %s lines.' % rje.integerString(len(obolines)))
            return self.parseGO(obolines)                    
        except: self.log.errorLog('GO.readGO() failed')
        return False
#########################################################################################################################
    def parseGO(self,glines,clear=True,obselete=False):   ### Parses GO Data from list of glines from OBO file
        '''
        Parses GO Data from list of glines from OBO file.
        >> glines:list of text lines read from OBO file
        >> clear:opt [True] = Whether to clear self.dict before reading in data
        >> obselete:opt [False] = Whether to read in obselete terms
        << returns True/False depending on success
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear:
                self.dict['AltID'] = {}
                self.dict['GO'] = {}
                self.dict['Subset'] = {}
            id = 'subsets'         # Current term being parsed
            ### ~ [2] ~ Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (gx,gtot) = (0.0,len(glines))
            while glines:
                self.printLog('\r#PARSE','Parsing %s GO terms: %.1f%%' % (rje.integerString(len(self.dict['GO'])),gx/gtot),newline=False,log=False)
                gx += 100.0
                ## ~ [2a] ~ Establish ID of current GO terms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gline = rje.chomp(glines.pop(0))
                if rje.matchExp('^id:\s+GO:(\d+)',gline):
                    id = rje.matchExp('^id:\s+GO:(\d+)',gline)[0]
                    self.dict['GO'][id] = {}
                    continue
                elif not id: continue
                elif rje.matchExp('^(\S+):\s+(\S.+)$',gline): (type,data) = rje.matchExp('^(\S+):\s+(\S.+)$',gline)
                elif gline[:1] in ['[','']: id = ''; continue
                ## ~ [2b] ~ Parse details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    if type == 'is_obsolete' and data.lower()[:4] == 'true':
                        self.dict['GO'].pop(id)
                        id = ''
                    elif type in ['name','def']: self.dict['GO'][id][type] = data
                    elif rje.matchExp('^subsetdef: (\S+) \"(\S.+)\"',gline):
                        (subset,desc) = rje.matchExp('^subsetdef: (\S+) \"(\S.+)\"',gline)
                        self.dict['Subset'][subset] = {'name':desc,'terms':[]}
                    elif type == 'namespace':
                        g = string.split(data,'_')
                        self.dict['GO'][id]['type'] = '%s%s' % (g[0][0],g[1][0])
                    elif type in ['is_a','relationship']:
                        parent = rje.matchExp('GO:(\d+)',data)[0]
                        if type != 'is_a': type = string.split(data)[0]
                        if type not in self.list['ParentTerms']: self.list['ParentTerms'].append(type)
                        if type not in self.dict['GO'][id]: self.dict['GO'][id][type] = []
                        self.dict['GO'][id][type].append(parent)
                    elif type == 'subset': self.dict['Subset'][string.split(gline)[1]]['terms'].append(id)
                    elif type == 'alt_id':
                        alt_id = rje.matchExp('GO:(\d+)',data)[0]
                        if alt_id in self.dict['AltID']: self.dict['AltID'][alt_id].append(id)
                        else: self.dict['AltID'][alt_id] = [id]
                    elif type in ['xref','synonym']:
                        if type not in self.dict['GO'][id]: self.dict['GO'][id][type] = []
                        self.dict['GO'][id][type].append(data)
                except: self.errorLog('GO.parseGO(%s) error' % gline)
            self.printLog('\r#PARSE','Parsed %s GO terms and %d subsets.' % (rje.integerString(len(self.dict['GO'])),len(self.dict['Subset'])))
            ### ~ [3] ~ Tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.makeChildren()
            self.makeGOSlim()
            for subset in self.dict['Subset']: self.dict['Subset'][subset]['terms'].sort()
            self.list['ParentTerms'].sort()
            return True
        except: self.log.errorLog('GO.parseGO() failed')
        return False
#########################################################################################################################
    ### <3> ### Additional Linking Methods                                                                              #
#########################################################################################################################
    def makeChildren(self): ### Goes through GO dictionary and adds 'child_terms' to dictionary
        '''Goes through GO dictionary and adds 'child_terms' to dictionary.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            go = self.go()
            for id in rje.sortKeys(go): go[id]['child_terms'] = []
            (gx,gtot) = (0.0,len(go)*2)
            ### ~ [2] ~ Add children ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for id in rje.sortKeys(go):
                self.progLog('#GO','Making GO children: %.2f%%' % (gx/gtot)); gx += 100.0
                for term in self.list['ParentTerms']:
                    if term in go[id]:
                        for parent in go[id][term]:
                            if parent not in go: self.errorLog('%s %s ID "%s" missing!' % (id,term,parent),printerror=False)
                            else: go[parent]['child_terms'].append(id)
            for id in rje.sortKeys(go):
                self.progLog('#GO','Making GO children: %.2f%%' % (gx/gtot)); gx += 100.0
                go[id]['child_terms'].sort()
            self.printLog('\r#GO','Making GO children complete.')
        except:
            self.errorLog('Major problem with GO.makeChildren()')
            raise
#########################################################################################################################
    def children(self,id): return self.go()[id]['child_terms']
#########################################################################################################################
    def parents(self,id,all_levels=True):  ### Returns all 'parents' for given ID
        '''Returns all 'parents' for given ID.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            parents = []
            go = self.go(id)
            ### ~ [2] ~ Add parent terms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for term in self.list['ParentTerms']:
                if term in go:
                    for parent in go[term]:
                        if parent not in parents and parent in self.go():
                            parents.append(parent)
                            if not all_levels: continue
                            for gramp in self.parents(parent):
                                if gramp not in parents: parents.append(gramp)
            parents.sort()
            return parents
        except: self.errorLog('Problem with GO.parents(%s)' % id)
#########################################################################################################################
    ### <4> ### GO Slim Methods                                                                                         #
#########################################################################################################################
    def makeGOSlim(self):   ### Adds a subset using self.list['GOSlim']
        '''Adds a subset using self.list['GOSlim'].'''
        self.dict['Subset']['goslim_user'] = {'name':'User-defined GO Slim','terms':[]}
        for id in self.list['GOSlim']:
            if id[:3] == 'GO:': id = id[3:]
            if id in self.go(): self.dict['Subset']['goslim_user']['terms'].append(id)
            else: self.errorLog('Unrecognised/Obselete GO Slim ID: %s' % id,printerror=False)
        if not self.dict['Subset']['goslim_user']['terms']:
            self.dict['Subset'].pop('goslim_user')
            return
        self.dict['Subset']['goslim_user']['terms'].sort()
        self.printLog('#GOSLIM','%d GO Terms added to goslim_user Subset' % len(self.dict['Subset']['goslim_user']['terms']))
#########################################################################################################################
    def powerGO(self,numbers,sig=0.01,samples='all',total='Total',countkey='counts',ignore=[]):  ### Special GO power calculation for GO slim set
        '''
        Special GO power calculation for GO slim set.
        >> numbers:dictionary of {Sample:Count}
        >> sig:float [0.01] = Desired significance level to achieve. Currently uncorrected. Add Bonf/FDR with time.
        >> samples:str ['all'] = Whether sig must be achievable for 'any' or 'all' samples.
        >> total:str ['Total'] = Sample containing Total counts to compare against
        >> countkey:str ['counts'] = Key identifying count dictionary for each GO term and 'total' count sample
        - self.go(id)[countkey] = {Sample:count}
        >> ignore:list of Samples to ignore from calculation
        << returns a list of GO IDs that meet criteria
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            N = numbers[total]        # Total count for calculating expectations/probabilities
            nlist = []                  # List of counts for subsamples to be assessed
            for sample in numbers:
                if sample not in ignore + [total]: nlist.append(numbers[sample])
            nlist = rje.sortUnique(nlist,xreplace=False,num=True)
            ### ~ [2] ~ Generate Power Range ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            plist = []                  # List of acceptable Total counts for subset
            nx = 0.0
            for i in range(1,N+1):      # Look at all possible levels of occurrence
                self.progLog('#POW','Calculating GO term power: %.1f%%' % (nx/N))
                nx += 100.0
                ok = 0
                p = float(i) / N        # Probability of each gene having this term
                for n in nlist:         # Look at each subset
                    k1 = min(i,n)       # Want to look at largest possible count for sample-term pairing
                    k2 = max(0,n-(N-i)) # Also want to look at the likelihood of under-representation
                    if rje.binomial(k1,n,p,callobj=self) <= sig: ok += 1
                    elif (1 - rje.binomial(k2+1,n,p,callobj=self)) <= sig: ok += 1
                    #!# Add under-representation too! #!#
                    if ok and samples == 'any': break
                if (ok and samples == 'any') or ok == len(nlist): plist.append(i)
            self.printLog('\r#POW','Calculation of GO term power complete.',log=False)
            self.deBug(nlist)
            ### ~ [3] ~ Generate GO Slim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            terms = []
            (ix,itot) = (0.0,len(self.go()))
            for id in rje.sortKeys(self.go()):
                self.progLog('#POW','Assessing terms for power: %.1f%% (%s terms)' % (ix/itot,rje.iLen(terms)))
                ix += 100.0
                if self.go(id)[countkey][total] in plist: terms.append(id)
            self.printLog('\r#POW','Assessed terms for statistical power, p <= %s: %s GO terms' % (sig,rje.iLen(terms)))
            #!# Add correction terms #!#
            self.deBug(terms)
            return terms
        except: self.errorLog('Major problem with GO.powerGO()')
        return []
#########################################################################################################################
    def topTerms(self,slimx=20,parents=False,total='Total',countkey='counts'):  ### Selects top terms for GO slim set
        '''
        Selects top terms for GO slim set.
        >> slimx:int [20] = Desired min. number of terms for each GO domain.
        >> parents:bool [False] = Whether parents and children both allowed in list
        >> total:str ['Total'] = Sample containing Total counts for assessment
        >> countkey:str ['counts'] = Key identifying count dictionary for each GO term and 'total' count sample
        - self.go(id)[countkey] = {Sample:count}
        << returns a list of GO IDs that meet criteria
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#self.opt['DeBug'] = True
            terms = []                          # List of terms
            dom = {'cc':{},'bp':{},'mf':{}}     # Dictionary of {domain:{count:[IDs]}}
            for id in self.go():
                n = self.go(id)[countkey][total]
                type = self.go(id)['type']
                if n not in dom[type]: dom[type][n] = [id]
                else: dom[type][n].append(id)
            ### ~ [2] ~ Generate Top Terms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(dom)
            for type in dom:
                dterms = []                     # Terms for this domain only
                dkeys = rje.sortKeys(dom[type]) # Counts, low to high
                dkeys.reverse()                 # Counts, high to low
                (dx,dtot) = (0.0,len(dkeys))
                while dkeys and len(dterms) < slimx: # Keep looping
                    self.deBug('%s: %s' % (type,dterms))
                    self.progLog('#TOP','Generating top %d %s terms: %.1f%%' % (slimx,type,dx/dtot))
                    dx += 100.0
                    n = dkeys.pop(0)            # Remove from list
                    dterms += dom[type][n]      # Add terms to term list
                    if parents: continue        # Don't care if parents and children all mixed up
                    for id in dterms[0:]:
                        if id not in dterms: continue               # Previously-removed parent
                        for par in self.parents(id):                # Check all parents
                            if par in dterms: dterms.remove(par)    # Remove parent term
                self.printLog('\r#TOP','Identified %s top %s terms: >= %s genes' % (rje.iLen(dterms),type,rje.iStr(n)))
                terms += dterms                 # Found a stable list of terms
            self.deBug(terms)
            return terms
        except: self.errorLog('Major problem with GO.topTerms()')
        return []
#########################################################################################################################
    ### <5> ### GO Mapping Methods                                                                                      #
#########################################################################################################################
    def mapEnsGO(self,spec='HUMAN',gokey='EnsGO',fixhead=True):   ### Extracts EnsEMBL GO mapping data from a BioMart download
        '''Extracts EnsEMBL GO mapping data from a BioMart download.'''
        ### ~ [1] ~ Setup paths and files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if gokey not in self.dict: self.dict[gokey] = {}
        ensmap = []
        for gtype in ['GO','GO.BP','GO.CC','GO.MF']:
            gfile = self.info['EnsGOPath'] + 'ens_%s.%s.tdt' % (spec,gtype)
            if os.path.exists(gfile): ensmap.append(gfile)
        if not ensmap:
            self.errorLog('EnsEMBL-GO mapping file (%s) missing' % self.info['EnsGOPath'],printerror=False)
            return False             
        ### ~ [2] ~ Parse Gene-GO Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        mainkeys = ['Ensembl Gene ID','GO ID']
        for gfile in ensmap:
            if fixhead:
                headers = string.split(rje.chomp(open(gfile,'r').readlines()[0]),'\t')
                if 'Ensembl Gene ID' in headers: mainkeys = ['Ensembl Gene ID']
                else: mainkeys = headers[:1]
                if 'GO Term Accession' in headers: mainkeys.append('GO Term Accession')
                elif 'GO Term Accession (bp)' in headers: mainkeys.append('GO Term Accession (bp)')
                elif 'GO Term Accession (mf)' in headers: mainkeys.append('GO Term Accession (mf)')
                elif 'GO Term Accession (cc)' in headers: mainkeys.append('GO Term Accession (cc)')
                elif 'GO ID' in headers: mainkeys.append('GO ID')
                else: mainkeys.append(headers[2])
                self.printLog('#HEAD','%s' % (string.join(mainkeys,' / ')))
            self.progLog('\r#GO','Mapping EnsEMBL GO...')
            ensdata = rje.dataDict(self,gfile,mainkeys)
            (mx,mtot) = (0.0,len(ensdata))
            obselete_go = []
            for map in ensdata:
                self.progLog('\r#GO','Mapping EnsEMBL GO: %.2f%%' % (mx/mtot)); mx += 100.0
                try: (gene,go) = string.split(map)
                except: continue    # no GO!
                ## Update dictionaries ##
                if go[:3] == 'GO:': go = go[3:]
                if go in self.go(): self.addGeneGO(gene,go,gokey)
                elif go in self.dict['AltID']:
                    for id in self.dict['AltID'][go]: self.addGeneGO(gene,id,gokey)
                elif go not in obselete_go: obselete_go.append(go)
            self.printLog('\r#GO','Mapping EnsEMBL GO from %s complete.' % os.path.basename(gfile))
            #self.verbose(0,0,self.dict['GO'],1)
#########################################################################################################################
    def addGeneGO(self,gene,go,gokey):    ### Adds gene-go relationship (and all parents) to self.dict['GO']
        '''Adds gene-go relationship (and all parents) to self.dict[gokey].'''
        if gene not in self.dict[gokey]: self.dict[gokey][gene] = []
        if go in self.dict[gokey][gene]: return
        self.dict[gokey][gene].append(go)
        for parent in self.parents(go): self.addGeneGO(gene,parent,gokey)
#########################################################################################################################    
    def makeGOGenes(self,gokey='EnsGO'):  ### Makes a dictionary of {GO:[Genes]}
        '''Makes a dictionary of {GO:[Genes]}.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gokey in ['name','is_a','part_of','type','child_terms']:
                self.errorLog('Cannot have "%s" as GOGenes key - reserved for GO' % gokey); raise ValueError
            if gokey not in self.dict:
                self.errorLog('"%s" mappings missing!' % gokey); raise ValueError
            ### ~ [2] ~ Process GO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (gx,gtot,ix) = (0.0,len(self.dict[gokey]),0)
            for gene in rje.sortKeys(self.dict[gokey]):
                self.progLog('\r#GENES','Making GO Gene lists: %.1f%%' % (gx/gtot)); gx += 100.0
                for go in self.dict[gokey][gene]:
                    if gokey in self.dict['GO'][go]: self.dict['GO'][go][gokey].append(gene)
                    else: self.dict['GO'][go][gokey] = [gene]; ix += 1
            self.printLog('\r#GENES','Making GO Gene lists complete: %s GO IDs with genes' % rje.integerString(ix))
        except: self.errorLog('Major problem with GO.makeGOGenes()')
#########################################################################################################################    
    def getGOGenes(self,id,gokey='EnsGO'):  ### Returns genes for given GO ID
        '''Returns genes for given GO ID.'''
        if gokey not in self.dict['GO'][id]: return []
        else: return self.dict['GO'][id][gokey]
#########################################################################################################################    
### End of SECTION II: GO Class                                                                                         #
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
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        GO(mainlog,cmd_list).test()
        print '\n\n *** No standalone functionality! *** \n\n'
        print rje_zen.Zen().wisdom()
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
