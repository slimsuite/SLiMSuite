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
Module:       rje_synteny
Description:  Gene/Protein Annotation Transfer and Synteny Tool
Version:      0.0.3
Last Edit:    26/06/17
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is based on the TopHits analysis of PAGSAT. It takes as input:
    - A Reference Genome sequence [refgenome=FILE]
    - A Reference Genome feature table with Gene and CDS features [ftfile=FILE]
    - GABLAM search results (e.g. from PAGSAT) of Genes/Proteins against assembly regions [ftgablam=FILE]




Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    refgenome=FILE  : Fasta/Genbank file of reference genome for assessment (also *.gb for full functionality) [None]
    ftfile=FILE     : File of reference genome features [*.Feature.tdt]
    ftgablam=FILE   : GABLAM search results (e.g. from PAGSAT) of Genes/Proteins against assembly regions [None]
    basefile=FILE   : Basename for output files [ftgablam basefile]
    ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    proteins=T/F    : Whether input is protein sequences (True) or genes (False) [False]
    tophitbuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genbank, rje_obj, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.0.1 - Altered problematic ValueError to warnLog()
    # 0.0.2 - Updated the synteny mappings to be m::n instead of m:n for Excel compatibility.
    # 0.0.3 - Added catching of the Feature locus/accnum mismatch issue.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [ ] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Make setupReference() method in rje_genbank and use for this (and PAGSAT?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_SYNTENY', '0.0.3', 'June 2017', '2016')
    description = 'Gene/Protein Annotation Transfer and Synteny Tool'
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
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
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
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
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
### SECTION II: Synteny Class                                                                                           #
#########################################################################################################################
class Synteny(rje_obj.RJE_Object):
    '''
    Synteny Class. Author: Rich Edwards (2015).

    Str:str
    - FTFile=FILE     : File of reference genome features [*.Feature.tdt]
    - FTGablam=FILE   : GABLAM search results (e.g. from PAGSAT) of Genes/Proteins against assembly regions [None]
    - RefBase=X       : Basefile for reference genome for assessment (*.gb) [None]
    - RefGenome=FILE  : Fasta file of reference genome for assessment (also *.gb for full functionality) [None]

    Bool:boolean
    - Proteins=T/F    : Whether input is protein sequences (True) or genes (False) [False]

    Int:integer

    Num:float
    - TopHitBuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FTFile','FTGablam','RefBase','RefGenome']
        self.boollist = ['Proteins']
        self.intlist = []
        self.numlist = ['TopHitBuffer']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({'TopHitBuffer':1.0})
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
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['FTFile','FTGablam','RefGenome'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['Proteins'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                self._cmdReadList(cmd,'float',['TopHitBuffer']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.topHitSynteny(self.db('features'),self.db('ftgablam'),proteins=self.getBool('Proteins'))
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup Reference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Set str['RefBase'], str['RefGenome'], bool['Features'] and str['FTFile']
            rje_genbank.setupRefGenome(self)
            fullfas = '%s.full.fas' % self.getStr('RefBase')
            if self.getStr('RefGenome') == fullfas or not rje.exists(self.getStr('RefGenome')):
                self.printLog('#NAMES','%s sequence names will not be suitable.' % fullfas)
                self.printLog('#NAMES','Please modify gene names in %s (e.g. ChrX) and save as %s.fas.' % (fullfas,self.getStr('RefBase')))
                self.printLog('#NAMES','Then re-run with refgenome=%s.fas.' % self.getStr('RefBase'))
                return False
            if not self.getBool('Features'):
                raise IOError('No Synteny analysis without Genbank features table.')
            self.obj['RefSeq'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('RefGenome'),'autoload=T','seqmode=file'])
            ### ~ [2] Setup Database Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            ftdb = db.addTable(self.getStr('FTFile'),mainkeys=['locus','feature','position'],name='features')
            if self.getBool('Proteins'): ftdb.dropEntriesDirect('feature',['CDS'],inverse=True)
            else: ftdb.dropEntriesDirect('feature',['gene'],inverse=True)
            ftdb.dataFormat({'start':'int','end':'int'})
            gdb = db.addTable(self.getStr('FTGablam'),mainkeys=['Qry','Hit'],name='ftgablam',expect=True)
            if not gdb: raise IOError('GABLAM file "%s" missing!' % (self.getStr('FTGablam')))
            ### ~ [3] Setup Basefile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.baseFile(return_none=None): self.baseFile(rje.baseFile(self.getStr('FTGablam')))
            db.baseFile(self.baseFile())
            ### ~ [4] Setup Run Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('TopHitBuffer') < 0: self.warnLog('Cannot have TopHitBuffer < 0.0!'); self.setNum({'TopHitBuffer': 0.0})
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
    def topHitSynteny(self,ftdb,gdb,proteins=False,tabname='synteny'): ### Performs TopHits and Synteny analysis on ftdb features.
        '''
        Performs TopHits and Synteny analysis on ftdb features.
        >> ftdb:Table = Database table of reference genome features for analysis (gene or CDS)
        >> gabdb:Table = Database table of features versus assembly GABLAM.
        >> proteins:bool [False] = Whether analysis is Proteins (else Genes)
        >> tabname:str ['synteny'] = Top Hits table name.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if proteins: gtype = 'Proteins'
            else: gtype = 'Genes'
            self.printLog('#~~#','## ~~~~~ %s TopHits/Synteny Analysis ~~~~~ ##' % gtype)
            # RefSeq is the reference genome sequence in a SeqList object
            refseq = self.obj['RefSeq']
            acc2chr = {}    # This maps sequence accession numbers onto reference chromosomes
            for seq in refseq.seqs(): acc2chr[rje.split(refseq.seqAcc(seq),'.')[0]] = refseq.seqGene(seq)
            # Check for data integrity
            lfield = 'locus'
            if lfield not in ftdb.fields(): lfield = 'Locus'
            if lfield not in ftdb.fields(): self.warnLog('Could not find "locus" or "Locus" field in features table! No AccNum/Locus check.')
            else:
                for locus in ftdb.indexKeys(lfield):
                    if locus not in acc2chr:
                        self.warnLog('Feature locus "%s" missing from Reference Sequence. AccNum/Locus mismatch in Reference?' % locus)
                        ftdb.dropEntriesDirect(lfield,[locus])

            ### ~ [1] Generate initial hits table from features and GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Reformat Qry sequence names to match feature locus tags
            for entry in gdb.entries():
                if '__' in entry['Qry']: entry['Qry'] = rje.split(entry['Qry'],'__')[-1]
                else: entry['Qry'] = rje.split(entry['Qry'],'_')[-1]
            # Join tables
            fdb = db.joinTables(name=tabname,join=[(ftdb,'locus_tag'),(gdb,'Qry')],newkey=['locus_tag','Hit'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
            # Update table fields
            gfield = 'Contig%s' % gtype[:4]     # Will rename to SPXXXXX.X.[P|G]Y
            fdb.addField(gfield)
            gfield = 'Hit%s' % gtype[:4]     # Will rename to SPXXXXX.X.[P|G]Y -> Gets reference gene when 1:1 or synteny
            fdb.addField(gfield)
            fdb.setFields(['locus_tag','Hit','Contig%s' % gtype[:4],gfield,'QryLen','Qry_AlnLen','Qry_AlnID','Qry_Start','Qry_End','Hit_Start','Hit_End','Hit_Dirn','locus','position','start','end','product','gene_synonym','db_xref','note'])
            fdb.dataFormat({'Qry_AlnID':'num','Qry_Start':'int'})
            fdb.addField('Chrom',after='locus')
            fdb.addField('Contig',after='Chrom')
            fdb.renameField('locus_tag','Gene')
            fdb.renameField('gene_synonym','Synonym')
            fdb.renameField('db_xref','XRef')
            fdb.renameField('Hit_Dirn','Dirn')
            for field in ['locus','position','start','end','product','note']: fdb.renameField(field,rje.strSentence(field))
            # Reformat assembly sequence names
            for entry in fdb.entries():
                # Assign strand data >>, ><, <> or << [ref dirn, hit dirn]
                if entry['Dirn'] == 'Bwd': entry['Dirn'] = '<'
                else: entry['Dirn'] = '>'
                if entry['Position'].startswith('complement'): entry['Dirn'] = '<' + entry['Dirn']
                else: entry['Dirn'] = '>' + entry['Dirn']
                # Chromosome
                entry['Chrom'] = acc2chr[entry['Locus']]
                if not entry['Hit']: continue
                # Pull out contig and position data
                [entry['Contig'],strain,na,entry[gfield]] = rje.split(entry['Hit'],'_')
                try: # Fragment names in form: hcq10_MBG11A__SP16495.10.011478-016127
                    if '-' in rje.split(entry['Hit'],'.')[-1]:
                        cdata = rje.split(entry['Hit'],'.')[-1]  # Position info
                        cdata = rje.split(cdata,'-')
                    else: # Old fragments in form: hcq10_MBG11A__SP16495.10-011478.016127
                        cdata = rje.split(entry['Hit'],'-')[-1]  # Position info
                        cdata = rje.split(cdata,'.')
                    entry['Hit_Start'] = int(cdata[0])
                    entry['Hit_End'] = int(cdata[1])
                except: self.errorLog('Problem with %s [%s]' % (entry['Hit'],entry),quitchoice=True)

            ### ~ [2] Report missing genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nullx = len(fdb.indexDataList('Hit','','Gene'))
            self.printLog('#HITS','%s of %s %s have no Assembly hit.' % (rje.iStr(nullx),rje.iLen(fdb.indexKeys('Gene')),gtype))

            ### ~ [3] Rename genes/proteins and check for overlaps (shouldn't be any) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prevhit = (None,None,-1,-1)                  # Tuples of (ctg,acc,start,end)
            gx = 0; hx = len(fdb.indexKeys('Hit'))  # Gene/Protein counter and total
            self.progLog('\r#HITS','%s...' % gfield)
            #self.debug(fdb.indexKeys('Hit')[:100])
            #self.debug(fdb.indexKeys('Hit')[-100:])
            for ghit in fdb.indexKeys('Hit'):
                # Compare start with prev end
                if not ghit: continue
                thishit = rje.matchExp('^(\S+)_\S+__(\S+)-(\d+)\.(\d+)$',ghit)
                if not thishit: thishit = rje.matchExp('^(\S+)_\S+__(\S+)\.(\d+)-(\d+)$',ghit)
                if prevhit and prevhit[0] == thishit[0]:
                    if int(thishit[-2]) < int(prevhit[-1]):
                        #raise ValueError('Overlapping fragment %s vs %s. Should not happen!' % (ghit,prevhit))
                        self.warnLog('Overlapping fragment %s vs %s. Should not happen!' % (ghit,prevhit),'overlap',quitchoice=True)
                prevhit = thishit
                # Rename protein/gene
                gx += 1
                for entry in fdb.indexEntries('Hit',ghit):
                    entry['Contig%s' % gtype[:4]] = entry[gfield] = '%s.%s%s' % (thishit[1],gtype[:1],rje.preZero(gx,hx))

            ### ~ [4] Generate TopHits assignments and filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # TopHitBuffer determines how similar to the TopHit a gene match can be. (Based on Qry_AlnID)
            # The bestentries dictionary will contain 1(+ if tying) best matches
            # All other entries that are not within the TopHitBuffer (for BOTH sides) will be deleted
            # The remaining entries are what defines the m:n mapping numbers
            # Any 1:n mappings have the reference gene identity transferred to gfield
            tophitbuffer = self.getNum('TopHitBuffer')
            ## ~ [4a] Identify best hits (each direction) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            besthit = {}            # Dictionary of {gene:best %ID} for ref and assembly
            bestentries = {}        # Dictionary of {gene:[best entries]} for ref and assembly
            self.progLog('\r#HITS','Best %s hits...' % gtype.lower()[:-1])
            for entry in fdb.entries():
                if not entry['Qry_AlnID']: continue
                for gid in (entry['Gene'],entry[gfield]):
                    if gid not in besthit or entry['Qry_AlnID'] > besthit[gid]:
                        besthit[gid] = entry['Qry_AlnID']
                        bestentries[gid] = [entry]
                    elif entry['Qry_AlnID'] == besthit[gid]: bestentries[gid].append(entry)
            ## ~ [4b] Filter anything beyond TopHitBuffer of best hits (both directions) ~~~~~~~~~~ ##
            ftotx = fdb.entryNum()
            self.progLog('\r#HITS','Top %s hit buffer...' % gtype.lower()[:-1])
            for entry in fdb.entries():
                if not entry['Qry_AlnID']: continue
                # Only keep matches where BOTH sequences are within tophitbuffer
                if entry['Qry_AlnID'] < (max(besthit[entry['Gene']],besthit[entry[gfield]]) - tophitbuffer):
                    fdb.dropEntry(entry)
            self.printLog('#HITS','%s %s reduced to %s within %.2f%% of best hit.' % (rje.iStr(ftotx),gtype,rje.iStr(fdb.entryNum()),tophitbuffer))
            ## ~ [4c] Classify mappings (1:1,1:n,n:1,n:n) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fdb.addField('Mapping')
            fdb.index('Gene'); fdb.index(gfield)
            self.progLog('\r#HITS','Classify m:n %s mappings...' % gtype.lower()[:-1])
            for entry in fdb.entries():
                if not entry['Qry_AlnID']: entry['Mapping'] = '1::0'; continue
                n = len(fdb.index('Gene')[entry['Gene']])
                m = len(fdb.index(gfield)[entry[gfield]])
                entry['Mapping'] = '%s::%s' % (m,n)
            ## ~ [4d] Identify 1:x mappings and map ref ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#HITS','Identify and map 1:x %s mappings...' % gtype.lower()[:-1])
            for m2n in fdb.index('Mapping'):
                if m2n[:2] != '1:': continue
                if m2n == '1::0': continue
                for gene in fdb.indexDataList('Mapping',m2n,'Gene'):
                    gx = 1; gtot = len(fdb.index('Gene')[gene])
                    for entry in fdb.indexEntries('Gene',gene):
                        if gtot == 1: entry[gfield] = gene      # 1:1 = easy straight mapping
                        elif entry in bestentries[gene] and len(bestentries[gene]) == 1: entry[gfield] = gene
                        elif gtot > 1: entry[gfield] = '%s|%d' % (gene,gx); gx += 1
            #fdb.indexReport('Mapping') # Report subset?

            ### ~ [5] Synteny Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Initially, look at direct synteny. Might want to extend at some point.
            # entry['Dirn'] = strand data >>, ><, <> or << [ref dirn, hit dirn]
            fdb.addFields(['Chr5','Chr3','Ctg5','Ctg3','Synteny'])
            # Chr5 & Chr3 = Genes 5' & 3' of Reference Gene (accounting for gene Dirn)
            # Ctg5 & Ctg3 = Genes 5' & 3' of Assembly Gene (accounting for gene Dirn)
            synteny_updated = True  # Boolean marker whether to keep looping
            geneord = {}            # Gene order list per chromosome/contig
            genedir = {}            # Gene direction list per chromosome/contig
            sloop = 0               # Number of Synteny loops
            while synteny_updated:
                synteny_updated = False
                ## ~ [5a] Reference gene order and direction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#SYN','Ref gene order...')
                for chrom in fdb.index('Chrom'): geneord[chrom] = [chrom]; genedir[chrom] = ['|']
                for entry in fdb.sortedEntries('Start'):
                    chrom = entry['Chrom']
                    if geneord[chrom][-1] != entry['Gene']:     # May have duplicates thanks to 1:n mapping
                        geneord[chrom].append(entry['Gene'])
                        genedir[chrom].append(entry['Dirn'][0])
                for chrom in fdb.index('Chrom'): geneord[chrom] += [chrom]; genedir[chrom] += ['|']
                ## ~ [5b] Assembly gene order and direction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#SYN','Assembly gene order...')
                for chrom in fdb.index('Contig'):
                    if not chrom: continue
                    geneord[chrom] = [chrom]; genedir[chrom] = ['|']
                #self.debug(fdb.sortedEntries('Hit')[:10])
                #self.debug(fdb.sortedEntries('Hit')[-10:])
                for entry in fdb.sortedEntries('Hit'):
                    chrom = entry['Contig']
                    if entry[gfield] and geneord[chrom][-1] != entry[gfield]:     # May have duplicates thanks to n:n mapping
                        geneord[chrom].append(entry[gfield])
                        genedir[chrom].append(entry['Dirn'][0])
                for chrom in fdb.index('Contig'):
                    if not chrom: continue
                    geneord[chrom] += [chrom]; genedir[chrom] += ['|']
                ## ~ [5c] Update synteny of assembly genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Syn = 5' and 3' match
                # Syn5 & Syn3 = only 5' or 3' match
                # Inv, Inv5 & Inv3 = As Syn, Syn5 and Syn3 but gene direction reversed
                # Tan, Tan5 & Tan3 = As Syn
                # *TR as above but 2+ copies of self-gene between flanking syntenic genes
                ex = 0.0; etot = fdb.entryNum()
                for entry in fdb.entries():
                    try:
                        self.progLog('\r#SYN','Updating synteny: %.2f%%' % (ex/etot)); ex += 100.0
                        ## Skip if classification will not change
                        if entry['Gene'] == 'YDL075W': self.debug(entry)
                        if entry['Synteny'] in ['Syn','Inv','SynTR','InvTR','NA']: continue
                        prevsyn = entry['Synteny']  # Check for updated synteny
                        ## Establish flanking genes
                        chrom = entry['Chrom']      # Set chromosome
                        entry['Chr5'] = geneord[chrom][geneord[chrom].index(entry['Gene'])-1]
                        entry['Chr3'] = geneord[chrom][geneord[chrom].index(entry['Gene'])+1]
                        if entry['Dirn'][0] == '<': (entry['Chr5'],entry['Chr3']) = (entry['Chr3'],entry['Chr5'])
                        if not entry[gfield]:   # This might get set by m:n management, when worked out
                            entry['Synteny'] = 'NA'
                            synteny_updated = True
                            continue
                        contig = entry['Contig']
                        entry['Tan5'] = entry['Ctg5'] = geneord[contig][geneord[contig].index(entry[gfield])-1]
                        entry['Tan3'] = entry['Ctg3'] = geneord[contig][geneord[contig].index(entry[gfield])+1]
                        ## Check for tandem repeat
                        tr = False
                        ti = 1
                        while rje.split(entry['Tan5'],'|')[0] == rje.split(entry[gfield],'|')[0]:  #!# Check initial names!
                            tr = True; ti += 1
                            entry['Tan5'] = geneord[contig][geneord[contig].index(entry[gfield])-ti]
                        ti = 1
                        while rje.split(entry['Tan3'],'|')[0] == rje.split(entry[gfield],'|')[0]:
                            tr = True; ti += 1
                            entry['Tan3'] = geneord[contig][geneord[contig].index(entry[gfield])+ti]
                        ## Update for direction
                        if entry['Dirn'][1] == '<':
                            (entry['Ctg5'],entry['Ctg3']) = (entry['Ctg3'],entry['Ctg5'])
                            (entry['Tan5'],entry['Tan3']) = (entry['Tan3'],entry['Tan5'])
                        if tr:
                            ctg5 = rje.split(entry['Tan5'],'|')[0]   # Adjust for gene|x 1:n mapping
                            ctg3 = rje.split(entry['Tan3'],'|')[0]   # Adjust for gene|x 1:n mapping
                        else:
                            ctg5 = rje.split(entry['Ctg5'],'|')[0]   # Adjust for gene|x 1:n mapping
                            ctg3 = rje.split(entry['Ctg3'],'|')[0]   # Adjust for gene|x 1:n mapping
                        ## Establish synteny
                        if ctg5 == entry['Chr5'] and ctg3 == entry['Chr3']: entry['Synteny'] = 'Syn'
                        elif ctg5 == entry['Chr5'] and ctg3 == chrom and entry['Chr3'] == contig: entry['Synteny'] = 'Syn'
                        elif ctg3 == entry['Chr3'] and ctg5 == chrom and entry['Chr5'] == contig: entry['Synteny'] = 'Syn'
                        elif ctg5 == entry['Chr3'] and ctg3 == entry['Chr5']: entry['Synteny'] = 'Inv'
                        elif ctg5 == entry['Chr5']: entry['Synteny'] = 'Syn5'
                        elif ctg5 == entry['Chr3']: entry['Synteny'] = 'Inv3'
                        elif ctg3 == entry['Chr5']: entry['Synteny'] = 'Inv5'
                        elif ctg3 == entry['Chr3']: entry['Synteny'] = 'Syn3'
                        else: entry['Synteny'] = 'Non'
                        if tr: entry['Synteny'] = entry['Synteny'][:3] + 'TR' + entry['Synteny'][3:]
                        ## Update n:1 gene mapping based on Synteny
                        # NB. bestentries = Dictionary of {gene:[best entries]} for ref and assembly
                        # where gene is gid in (entry['Gene'],entry[gfield]):
                        gid = entry[gfield]     # Old gene prior to update
                        synteny_updated = synteny_updated or prevsyn != entry['Synteny']    # Changed synteny
                        if entry['Synteny'] not in ['Syn','Inv']: continue
                        if entry['Mapping'] == '1::1': continue  # Nothing else to update
                        if not entry['Mapping'].endswith(':1'):
                            # Check for a single Syn rating for assembly gene and update, then delete rest of entries
                            synx = 0
                            for sentry in fdb.indexEntries(gfield,gid):
                                if sentry['Synteny'] in ['Syn','Inv']: synx += 1
                            if synx > 1: continue
                        entry[gfield] = entry['Gene']
                        gi = geneord[contig].index(gid)
                        #self.bugPrint(geneord[contig][gi-1:][:3])
                        geneord[contig][geneord[contig].index(gid)] = entry['Gene']
                        #self.bugPrint(geneord[contig][gi-1:][:3])
                        #self.debug(entry)
                        for gentry in fdb.indexEntries(gfield,gid):
                            if gentry[gfield] != gid: continue
                            if gentry[gfield] == gentry['Gene']: self.debug('Scrubbing!: %s' % gentry)
                            gentry[gfield] = gentry['Ctg5'] = gentry['Ctg3'] = ''
                    except:
                        self.errorLog('Oops')
                        self.debug(entry)
                        self.warnLog('Synteny problem for %s="%s"' % (entry['Gene'],entry[gfield]),warntype='Synteny',quitchoice=True,suppress=True,dev=True)
                fdb.indexReport('Synteny',force=True)
                if synteny_updated: sloop += 1; self.printLog('#LOOP','Loop %d complete: synteny updated.' % sloop)
                else: sloop += 1; self.printLog('#LOOP','Convergence after synteny loop %d.' % sloop)

            ### ~ [6] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [6a] Clean up m:n proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            remx = 0
            for gene in fdb.index('Gene',force=True):
                gentries = fdb.indexEntries('Gene',gene)
                if len(gentries) < 2: continue
                mapping = fdb.indexDataList('Gene',gene,gfield)
                if '' in mapping and len(mapping) > 1:
                    for gentry in gentries:
                        if not gentry[gfield]: fdb.dropEntry(gentry); remx +=1
            self.printLog('#REM','%s redundant m:n entries dropped following synteny mapping.' % rje.iStr(remx))

            ## ~ [6b] TopHits/Synteny Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in fdb.entries(): entry['Hit'] = rje.split(entry['Hit'],'-')[0]   # Update entry['Hit'] for R Script
            fdb.saveToFile()
            ## ~ [6c] GeneOrder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ordtxt = '%s.%s.order.txt' % (self.baseFile(),gtype)
            ORD = open(ordtxt,'w'); cx = 0
            for chrom in rje.sortKeys(geneord):
                chromord = []; cx += 1
                for i in range(len(geneord[chrom])):
                    chromord.append('%s%s%s' % (genedir[chrom][i],geneord[chrom][i],genedir[chrom][i]))
                ORD.write('%s\n' % rje.join(chromord))
                self.debug('%s\n' % rje.join(chromord))
            ORD.close()
            self.printLog('#ORD','%s order output to %s' % (gtype,ordtxt))

        except: self.errorLog('%s.topHitSynteny error' % self.prog())
#########################################################################################################################
    def PAGSATtopHits(self):    ### Generates output for Gene/Protein TopHits analysis
        '''Returns the Reference Features table, if given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT TopHits/Synteny Analysis ~~~~~ ##')
            db = self.db()
            if not self.force() and rje.exists('%s.Genes.TopHits.tdt' % db.baseFile()) and rje.exists('%s.Proteins.TopHits.tdt' % db.baseFile()):
                self.printLog('#HITS','Genes.TopHits and Proteins.TopHits tables found (force=F).'); return True
            tophitbuffer = self.getNum('TopHitBuffer')
            refseq = self.obj['RefSeq']
            ftdict = db.splitTable(self.ftdb(),'feature',asdict=True,keepfield=False,splitchar=None,values=['gene','CDS'])    # Reference features tables
            ftdict['Genes'] = ftdict.pop('gene')
            #self.debug(ftdict['Genes'].entries()[:10])
            ftdict['Proteins'] = ftdict.pop('CDS')
            #self.debug(ftdict['Proteins'].index('note').keys()[:10])
            acc2chr = {}
            for seq in refseq.seqs(): acc2chr[rje.split(refseq.seqAcc(seq),'.')[0]] = refseq.seqGene(seq)
            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gtype in ['Genes','Proteins']:
                ## ~ [1a] Load GABLAM data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.getBool('%sSummary' % gtype[:4]): continue
                gdb = db.addTable(filename='%s.%s.Fragments.gablam.tdt' % (self.fileBase('GABLAM','Base'),gtype),mainkeys=['Qry','Hit'],name='%s.gablam' % gtype,expect=True)
                if not gdb: self.warnLog('%s.%s.Fragments.gablam.tdt missing!' % (self.fileBase('GABLAM','Base'),gtype)); continue
                ## ~ [1b] Filter, join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #gdb.dropEntries(['Rank>1']) #!# No more! Filter on
                for entry in gdb.entries():
                    if '__' in entry['Qry']: entry['Qry'] = rje.split(entry['Qry'],'__')[-1]
                    else: entry['Qry'] = rje.split(entry['Qry'],'_')[-1]
                fdb = db.joinTables(name='%s.TopHits' % gtype,join=[(ftdict[gtype],'locus_tag'),(gdb,'Qry')],newkey=['locus_tag','Hit'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
                self.debug(fdb.datakeys()[:10])
                gfield = 'Contig%s' % gtype[:4]     # Rename to SPXXXXX.X.[P|G]Y
                fdb.addField(gfield)
                fdb.setFields(['locus_tag','Hit',gfield,'QryLen','Qry_AlnLen','Qry_AlnID','Qry_Start','Qry_End','HitLen','Hit_Start','Hit_End','locus','position','start','end','product','gene_synonym','db_xref','note'])
                fdb.dataFormat({'Qry_AlnID':'num','Qry_Start':'int'})
                fdb.addField('Chrom',after='locus')
                fdb.addField('Contig',after='Chrom')
                for entry in fdb.entries():
                    entry['Chrom'] = acc2chr[entry['locus']]
                    if not entry['Hit']: continue
                    [entry['Contig'],strain,na,entry[gfield]] = rje.split(entry['Hit'],'_')    # hcq10_MBG11A__SP16495.10-011478.016127
                    #x#cdata = rje.split(entry['Hit'],'-')[-1]  # Position info
                    #x#cdata = rje.split(cdata,'.')
                    #Now: # hcq10_MBG11A__SP16495.10.011478-016127
                    try:
                        if '-' in rje.split(entry['Hit'],'.')[-1]:
                            cdata = rje.split(entry['Hit'],'.')[-1]  # Position info
                            cdata = rje.split(cdata,'-')
                        else:
                            cdata = rje.split(entry['Hit'],'-')[-1]  # Position info
                            cdata = rje.split(cdata,'.')
                        entry['Hit_Start'] = int(cdata[0])
                        entry['Hit_End'] = int(cdata[1])
                    except: self.errorLog('Problem with %s [%s]' % (entry['Hit'],entry),quitchoice=True)
                fdb.renameField('locus_tag','Gene')
                fdb.renameField('gene_synonym','Synonym')
                fdb.renameField('db_xref','XRef')
                for field in ['locus','position','start','end','product','note']: fdb.renameField(field,rje.strSentence(field))
                #self.debug(fdb.index('Note').keys()[:10])
                ## ~ [1c] Report missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nullx = 0
                #for entry in fdb.entries():
                #    if not entry['Hit']: nullx += 1
                #self.printLog('#HITS','%s of %s %s have no Assembly hit.' % (rje.iStr(nullx),rje.iStr(fdb.entryNum()),gtype))
                nullx = len(fdb.indexDataList('Hit','','Gene'))
                self.printLog('#HITS','%s of %s %s have no Assembly hit.' % (rje.iStr(nullx),rje.iLen(fdb.indexKeys('Gene')),gtype))
                ## ~ [1d] Check overlapping hits & rename (should be any, right?) ~~~~~~~~~~~~~~~~~ ##
                prevhit = (None,-1,-1)  #?# Make tuple of (ctg,start,end)?
                gx = 0; hx = len(fdb.indexKeys('Hit'))
                self.progLog('\r#HITS','%s...' % gfield)
                for ghit in fdb.indexKeys('Hit'):
                    # Compare start with prev end
                    if not ghit: continue
                    self.bugPrint(ghit)
                    thishit = rje.matchExp('^(\S+)_\S+__(\S+)-(\d+)\.(\d+)$',ghit)
                    if not thishit: thishit = rje.matchExp('^(\S+)_\S+__(\S+)\.(\d+)-(\d+)$',ghit)
                    if prevhit and prevhit[0] == thishit[0]:
                        if int(thishit[-1]) < int(prevhit[-2]):
                            #raise ValueError('Overlapping fragment %s vs %s. Should not happen!' % (ghit,prevhit))
                            self.warnLog('Overlapping fragment %s vs %s. Should not happen!' % (ghit,prevhit),'overlap',quitchoice=True)
                    prevhit = thishit
                    # Rename protein/gene
                    gx += 1
                    for entry in fdb.indexEntries('Hit',ghit): entry[gfield] = '%s.%s%s' % (thishit[1],gtype[:1],rje.preZero(gx,hx))
                ## ~ [1e] Filter according to best hits (each direction) ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                besthit = {}            # Dictionary of {gene:best %ID}
                bestentries = {}        # Dictionary of {gene:[best entries]}
                self.progLog('\r#HITS','Best Hits...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: continue
                    for gid in (entry['Gene'],entry[gfield]):
                        #!# Make this a dictionary to store tied best hits! Add Best field
                        if gid not in besthit or entry['Qry_AlnID'] > besthit[gid]:
                            besthit[gid] = entry['Qry_AlnID']
                            bestentries[gid] = [entry]
                        elif entry['Qry_AlnID'] == besthit[gid]: bestentries[gid].append(entry)
                ftotx = fdb.entryNum()
                self.progLog('\r#HITS','Top Hit Buffer...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: continue
                    # Only keep matches where BOTH sequences are within tophitbuffer
                    if entry['Qry_AlnID'] < (max(besthit[entry['Gene']],besthit[entry[gfield]]) - tophitbuffer):
                        fdb.dropEntry(entry)
                self.printLog('#HITS','%s %s reduced to %s within %.2f%% of best hit.' % (rje.iStr(ftotx),gtype,rje.iStr(fdb.entryNum()),tophitbuffer))
                ## ~ [1f] Classify mappings (1:1,1:n,n:1,n:n) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fdb.addField('Mapping')
                fdb.index('Gene'); fdb.index(gfield)
                self.progLog('\r#HITS','Classify m:n mappings...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: entry['Mapping'] = '1::0'; continue
                    n = len(fdb.index('Gene')[entry['Gene']])
                    m = len(fdb.index(gfield)[entry[gfield]])
                    entry['Mapping'] = '%s::%s' % (m,n)
                ## ~ [1g] Identify 1:x mappings and map ref ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#HITS','Identify and map 1:x mappings...')
                for m2n in fdb.index('Mapping'):
                    if m2n[:2] != '1:': continue
                    if m2n == '1::0': continue
                    for gene in fdb.indexDataList('Mapping',m2n,'Gene'):
                        gx = 1; gtot = len(fdb.index('Gene')[gene])
                        for entry in fdb.indexEntries('Gene',gene):
                            if gtot == 1: entry[gfield] = gene      # 1:1 = easy straight mapping
                            elif entry in bestentries[gene] and len(bestentries[gene]) == 1: entry[gfield] = gene
                            elif gtot > 1: entry[gfield] = '%s|%d' % (gene,gx); gx += 1
                #fdb.indexReport('Mapping') # Report subset?
                ## ~ [1h] Synteny ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fdb.addFields(['Chr5','Chr3','Ctg5','Ctg3','Synteny'])
                synteny_updated = True; sloop = 0
                while synteny_updated:
                    synteny_updated = False
                    #i# Initially, look at direct synteny. Might want to extend at some point.
                    # Reference genes
                    synteny = {}
                    self.progLog('\r#SYN','Ref synteny...')
                    for chrom in fdb.index('Chrom'): synteny[chrom] = ['NA']
                    for entry in fdb.sortedEntries('Start'):
                        chrom = entry['Chrom']
                        if synteny[chrom][-1] != entry['Gene']: synteny[chrom].append(entry['Gene'])
                    for chrom in fdb.index('Chrom'): synteny[chrom] += ['NA']
                    # Assembly genes
                    self.progLog('\r#SYN','Assembly synteny...')
                    for chrom in fdb.index('Contig'): synteny[chrom] = ['NA']
                    for entry in fdb.sortedEntries('Hit'):
                        chrom = entry['Contig']
                        if synteny[chrom][-1] != entry[gfield]: synteny[chrom].append(entry[gfield])
                    for chrom in fdb.index('Contig'): synteny[chrom] += ['NA']
                    # Update synteny
                    self.progLog('\r#SYN','Updating synteny...')
                    for entry in fdb.entries():
                        if entry['Synteny'] in ['Full','Tandem']: continue
                        prevsyn = entry['Synteny']
                        chrom = entry['Chrom']
                        entry['Chr5'] = synteny[chrom][synteny[chrom].index(entry['Gene'])-1]
                        entry['Chr3'] = synteny[chrom][synteny[chrom].index(entry['Gene'])+1]
                        contig = entry['Contig']
                        entry['Ctg5'] = synteny[contig][synteny[contig].index(entry[gfield])-1]
                        entry['Ctg3'] = synteny[contig][synteny[contig].index(entry[gfield])+1]
                        ctg5 = rje.split(entry['Ctg5'],'|')[0]
                        ctg3 = rje.split(entry['Ctg3'],'|')[0]
                        if ctg5 == entry['Chr5'] and ctg3 == entry['Chr3']: entry['Synteny'] = 'Full'
                        elif ctg5 == entry['Chr3'] and ctg3 == entry['Chr5']: entry['Synteny'] = 'Full' # Reverse?
                        elif ctg5 in [entry['Chr5'],entry['Chr3']]: entry['Synteny'] = 'Partial'
                        elif ctg3 in [entry['Chr5'],entry['Chr3']]: entry['Synteny'] = 'Partial'
                        else: entry['Synteny'] = 'Orphan'
                        # Special identification of 2:1 "tandem compression" and 1:2 "tandem duplication"
                        # Special 2:1 or 1:2
                        if (entry['Mapping'].startswith('1:') or entry['Mapping'].endswith(':1')) and entry['Synteny'] == 'Partial':
                            chrneighbours = [entry['Chr5'],entry['Chr3']]
                            ctgneighbours = [ctg5,ctg3]
                            missing = rje.listDifference(chrneighbours,ctgneighbours)
                            if entry['Mapping'].endswith(':1') and missing:
                                missing = missing[0]
                                chrneighbours.append(synteny[chrom][synteny[chrom].index(missing)-1])
                                chrneighbours.append(synteny[chrom][synteny[chrom].index(missing)+1])
                                if len(rje.listIntersect(chrneighbours,ctgneighbours)) == 2: entry['Synteny'] = 'Tandem'
                            missing = rje.listDifference(ctgneighbours,chrneighbours)
                            if entry['Mapping'].startswith('1:') and missing:
                                missing = missing[0]
                                if missing not in synteny[contig]:
                                    variant = None
                                    #self.bugPrint('%s vs "%s|"' % (synteny[contig],missing))
                                    for gvar in synteny[contig]:
                                        if gvar.startswith('%s|' % missing):
                                            #self.bugPrint(gvar)
                                            if variant: variant = None; break
                                            variant = gvar
                                    missing = variant
                                    #self.deBug(missing)
                                if missing in synteny[contig]:
                                    ctgneighbours.append(synteny[contig][synteny[contig].index(missing)-1])
                                    ctgneighbours.append(synteny[contig][synteny[contig].index(missing)+1])
                                #else:
                                    #self.bugPrint(entry)
                                    #self.bugPrint(chrneighbours)
                                    #self.bugPrint(ctgneighbours)
                                    #self.bugPrint(missing)
                                    #self.debug(synteny[contig])
                                if len(rje.listIntersect(chrneighbours,ctgneighbours)) == 2: entry['Synteny'] = 'Tandem'
                        #if prevsyn != entry['Synteny'] and sloop: self.deBug('%s: %s => %s' % (entry,prevsyn,entry['Synteny']))
                        synteny_updated = synteny_updated or prevsyn != entry['Synteny']    # Change
                    # if x:1 keep Full synteny, else keep best hit
                    hitchanged = {}      # Will need to remake keys if Gene/Hit changes
                    for xn in fdb.indexKeys('Mapping'):
                        if not xn.endswith(':1'): continue
                        for gid in fdb.indexDataList('Mapping',xn,gfield):
                            syntenic = 'Full' in fdb.indexDataList(gfield,gid,'Synteny')   # Keep only syntenic
                            for entry in fdb.indexEntries(gfield,gid):
                                if syntenic and entry['Synteny'] != 'Full': fdb.dropEntry(entry)
                                elif not syntenic and besthit[gid] != entry['Qry_AlnID']: fdb.dropEntry(entry)
                                else: hitchanged[gid] = entry['Gene'] # Use (one of) best gene name
                    if hitchanged:
                        for gid in hitchanged:
                            for entry in fdb.indexEntries(gfield,gid): entry[gfield] = hitchanged[gid]
                        fdb.index(gfield,force=True)
                    # if x:y ...?
                    fdb.indexReport('Synteny')
                    if synteny_updated: sloop += 1; self.printLog('#LOOP','Loop %d complete: synteny updated.' % sloop)
                    else: sloop += 1; self.printLog('#LOOP','Convergence after synteny loop %d.' % sloop)

                ## ~ [1x] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries(): entry['Hit'] = rje.split(entry['Hit'],'-')[0]   # Update entry['Hit'] for R Script
                fdb.saveToFile()


        except: self.errorLog('%s.topHits error' % self.prog())
#########################################################################################################################
### End of SECTION II: Synteny Class                                                                                    #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Synteny(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
