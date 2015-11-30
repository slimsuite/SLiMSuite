#!/usr/bin/python

# See below for name and description
# Copyright (C) 2013 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SNP_Mapper
Description:  SNP consensus sequence to CDS mapping 
Version:      0.4.0
Last Edit:    17/11/15
Copyright (C) 2013  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is in development. The following documentation should be considered aspirational rather than accurate!

    This module is for mapping SNPs onto a genbank file (or similar feature annotation) and reporting the possible effect
    of SNPs in addition to more detailed output for individual genes and proteins. The main SNP output will be restricted
    to feature reporting and coding changes.

    NOTE: At present, this does not handle indels properly, nor multiple SNPS affecting the same amino acid. It cannot
    deal with introns.

    Primary input for this program is:

    1. Reference genome sequence. This should be a DNA sequence file in which the accession numbers match the `Locus`
    field for the SNP file (below). If a Genbank file is given, the `*.full.fas` file will be used.

    2. A feature table with the headers:
    locus	feature	position	start	end	product	gene_synonym	note	db_xref	locus_tag	details
    If this is not given, a `*.Feature.tdt` file based on `seqin=FILE` will be sought and loaded.

    3. The SNP file (`snpfile=FILE`) should have a `Locus`, `Pos`, `REF` and `ALT` fields, which will be used to map onto
    the features file and uniquely mark variants in the reference genome. If running in (default) altpos=T mode, this
    file will represent the mapping of a single pair of genomes: there will be no multiple-allele entries and the file
    headers should include AltPos and AltLocus for making unique entries. (Reference sequences can map to multiple query
    sequences.) If running on BCF output, alleles will be split on commas but there will be no compiled sequence output.
    (** Not yet implemented! **)

To be added:
    Sequence output and recognition of BCF files to be added.

Old Function:
    This module reads in an alignment of a coding sequence with a consensus sequence and a list of polymorphic sites
    relative to consensus. Polymorphisms are mapped onto the coding sequence and the desired output produced.

Commandline:
    ### ~ Basic SNP mapping functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Sequence input file with accession numbers matching Locus IDs, or Genbank file.  []
    spcode=X        : Overwrite species read from file (if any!) with X if generating sequence file from genbank [None]
    ftfile=FILE     : Input feature file (locus,feature,position,start,end) [*.Feature.tdt]
    ftskip=LIST     : List of feature types to exclude from analysis [source,telomere]
    ftbest=LIST     : List of features to exclude if earlier feature in list overlaps position [CDS,mRNA,gene,mobile_element]
    snpfile=FILE    : Input table of SNPs to map and output (should have locus and pos info, see above) []
    snphead=LIST    : List of SNP file headers []
    snpdrop=LIST    : List of SNP fields to drop []
    altpos=T/F      : Whether SNP file is a single mapping (with AltPos) (False=BCF) [True]
    basefile=FILE   : Root of output file names (same as SNP input file by default) []
    
    ### ~ Old Options (need reviving) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    batch=LIST  : List of alignment files to read X.aln - must have X_polymorphisms.txt too [*.aln]
    screenmatch=LIST: List of genotypes for which to screen out matching SNPs []

    ### ~ Obsolete Options (roll back to pre v0.3.0 if required) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genotypes=LIST  : List of snpfile headers corresponding to genotypes to map SNPs []
    genbase=FILE    : Basefile for Genbank output. Will use base of seqin if None []
    snpkeys=LIST    : Additional headers to use as keys for SNP file (e.g. if mapping done by chromosome) []

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_seqlist, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genbank, rje_obj, rje_seqlist, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation. Batch mode for mapping SNPs needs updating.
    # 0.1 - SNP mapping against a GenBank file.
    # 0.2 - Fixed complement strand bug.
    # 0.3.0 - Updated to work with RATT(/Mummer?) snp output file. Improved docs.
    # 0.4.0 - Major reworking for easier updates and added functionality. (Convert to 1.0.0 when complete.)
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Standard SNP table and add converters to this standard format. (Base on mpileup?)
    # [ ] : Fix indel handling.
    # [ ] : Break down into more general/reusable methods that can be applied to other mapping. (e.g. PAGSAT)
    # [ ] : Deal with introns!
    # [ ] : Add syn/NS/nonsense rating to aa changes (extra field)
    # [ ] : Update input to be able to handle a sequence file and features file.
    # [ ] : Add a summary of the SNP types for each Locus[/SNPKeys].
    # [ ] : If required:     snpalt=X        : SNP table alternative genotype field [ALT]
    # [ ] : In addition to any annotation extracted using rje_genbank, additional features can be given to be mapped on to the genomic sequence?
    # [ ] : Or an option to merge results from several runs (e.g. with different features?)
    # [ ] : Add sequence output for different features.
    # [ ] : Add ratings for CDS based on Ka/Ks and no. deletions etc.
    # [ ] : Add output of possible duplications (and total deletions?) - May need to compare to self-analysis to be sure.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('SNP_MAPPER', '0.4.0', 'November 2015', '2013')
    description = 'SNP consensus sequence to CDS mapping'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class SNPMap(rje_obj.RJE_Object):     
    '''
    SNPMap Class. Author: Rich Edwards (2013).

    Str:str
    - FTFile=FILE     : Input feature file (locus,feature,position,start,end) [*.Feature.tdt]
    - SeqIn = Sequence input file with accession numbers matching Locus IDs, or Genbank file.  []
    - SNPFile = Input table of SNPs to map and output (should have locus and pos info, see above) []
  
    Bool:boolean
    - AltPos=T/F      : Whether SNP file is a single mapping (with AltPos) (False=BCF) [True]

    Int:integer

    Num:float
    
    List:list
    - OLD: Batch = List of alignment files to read X.aln - must have X_polymorphisms.txt too [*.aln]
    - FTBest=LIST     : List of features to exclude if earlier feature in list overlaps position [CDS,mRNA,gene,mobile_element]
    - FTSkip=LIST     : List of feature types to exclude from analysis [source]
    - SNPDrop=LIST    : List of SNP fields to drop []
    - SNPHead=LIST    : List of SNP file headers []

    Dict:dictionary

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - GenBank = rje_genbank.GenBank object
    - SeqList = rje_seqlist.SeqList object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FTFile','SeqIn','SNPFile']
        self.boollist = ['AltPos']
        self.intlist = []
        self.numlist = []
        self.listlist = ['FTBest','FTSkip','SNPDrop','SNPHead']
        self.dictlist = []
        self.objlist = ['DB','GenBank','SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({})
        self.setBool({'AltPos':True})
        self.setInt({})
        self.setNum({})
        self.list['FTSkip'] = ['source','telomere']
        self.list['FTBest'] = ['CDS','mRNA','gene','mobile_element']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['TupleKeys=T'])
        self.obj['GenBank'] = rje_genbank.GenBank(self.log,self.cmd_list+['tabout=T'])
        self.obj['GenBank'].obj['DB'] = self.obj['DB']
        self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F'])
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
                self._cmdReadList(cmd,'file',['GenBase','SeqIn','SNPFile'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['AltPos'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['FTBest','FTSkip','SNPDrop','SNPHead'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Batch']) # List of files using wildcards and glob
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
            #!# Add method types #!#
            self.featureSNPs()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Load Sequences and Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupReference(): raise ValueError('Reference setup failed.')
            if not self.setupSNPTable(): raise ValueError('SNP Table setup failed.')
            self.db().baseFile(self.baseFile())
            return True
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupReference(self):   ### Loads input sequences and reference features.
        '''Loads input sequences and reference features.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gb = self.obj['GenBank']
            seqlist = self.obj['SeqList']
            seqin = self.getStr('SeqIn')
            gbftfile = None
            ### ~ [1] Load Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(seqin): raise IOError('Cannot find seqin=%s!' % seqin)
            if open(seqin,'r').readline()[:1] != '>':   # Not fasta file: assume Genbank. Should have features file.
                gbase = rje.baseFile(seqin)
                if gbase.endswith('.full'): gbase = rje.baseFile(gbase)
                gb.baseFile(gbase)
                if 'full' not in gb.list['FasOut']: gb.list['FasOut'].append('full')
                if not gb.loadFeatures() or not os.path.exists('%s.full.fas' % gb.baseFile()): gb.run()
                seqin = '%s.full.fas' % gb.baseFile()
                gbftfile = '%s.Feature.tdt' % gbase
                gb.db().deleteTable('Feature')   # Want to reload from table
                self.db().addTable(gbftfile,name='Feature',expect=True,mainkeys=['locus','feature','position'])
            seqlist.loadSeq(seqin)
            seqdict = seqlist.makeSeqNameDic('accnum')
            ### ~ [2] Load Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.db('Feature') and (self.getStr('FTFile') == gbftfile or not self.getStrLC('FTFile')):
                self.printLog('#FT','Genbank FTFile "%s" already loaded.' % gbftfile)
            elif not rje.exists(self.getStr('FTFile')): raise IOError('Cannot find ftfile=%s!' % self.getStr('FTFile'))
            elif self.db('Feature'):    # Add features
                ftdb = self.db().addTable(self.getStr('FTFile'),name='FTFile',expect=True,mainkeys=['locus','feature','position'])
                self.db().mergeTables(self.db('Feature'),ftdb)
            else: self.db().addTable(self.getStr('FTFile'),name='Feature',expect=True,mainkeys=['locus','feature','position'])
            ftdb = self.db('Feature')
            ftdb.dropEntriesDirect('feature',self.list['FTSkip'])
            ftdb.dataFormat({'start':'int','end':'int'})
            ftdb.newKey(['locus','start','feature','position'])    # Reorder for sorting
            #self.debug(ftdb.dataKeys()[:10])
            #self.debug(ftdb.entries()[:10])
            ### ~ [3] Check sequences and loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqloci = rje.sortKeys(seqdict)
            ftloci = ftdb.indexKeys('locus')
            okloci = rje.listIntersect(seqloci,ftloci)
            self.printLog('#LOCUS','%s sequence loci; %s feature loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(ftloci),rje.iLen(okloci)))
            self.printLog('#SETUP','Reference setup complete!')
            return True
        except: self.errorLog('Problem during %s setupReference.' % self); return False  # Setup failed
#########################################################################################################################
    def setupSNPTable(self):    ### Loads, reformats and checks SNP Table
        '''Loads, reformats and checks SNP Table.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('SNPFile')): raise IOError('Cannot find snpfile=%s!' % self.getStr('SNPFile'))
            if not self.getStrLC('Basefile'):
                self.baseFile(rje.baseFile(self.getStr('SNPFile'),strip_path=True))
            if self.getBool('AltPos'):
                skeys = ['Locus','Pos','AltLocus','AltPos']
            else:
                skeys = ['Locus','Pos','REF']
            ### ~ [1] Load SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db().addTable(self.getStr('SNPFile'),mainkeys=skeys,headers=self.list['SNPHead'],name='SNP')
            snpdb.dropFields(self.list['SNPDrop'])
            if 'ALT' not in snpdb.fields(): raise ValueError('Cannot find `ALT` field!')
            snpdb.dataFormat({'Pos':'int'})
            ### ~ [2] Reformat Indels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AltPos'):
                snpdb.dataFormat({'AltPos':'int'})
                snpdb.remakeKeys()
                #?# Need to convert indels at some point to BCF-style #?# Try without first #?#
            else:
                #!# Add reformatting to split multiple alleles #!#
                snpdb.remakeKeys()
            ### ~ [3] Add default SNPType ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb.addField('SNPType')
            for sentry in snpdb.entries():
                # NOTE: SNP Data should now be formatted either:
                # (a) like RATT output with Pos, AltLocus, AltPos, REF and ALT
                # (b) Like BCF with Pos, REF and ALT
                ## ~ [2a] First define type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if sentry['REF'] == sentry['ALT']: stype = 'REF'   # No mismatch!
                elif sentry['REF'] in ['.','-'] or len(sentry['REF']) < len(sentry['ALT']): stype = 'INS'
                elif sentry['ALT'] in ['.','-'] or len(sentry['REF']) > len(sentry['ALT']): stype = 'DEL'
                else: stype = 'SNP'
                sentry['SNPType'] = stype
            snpdb.dropEntriesDirect('SNPType',['REF'])
            ### ~ [4] Check and rename loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdb = self.db('Feature')
            seqloci = self.obj['SeqList'].seqNameDic().keys()
            snploci = snpdb.indexKeys('Locus')
            okloci = rje.listIntersect(seqloci,snploci)
            self.printLog('#LOCUS','%s sequence loci; %s SNP loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(snploci),rje.iLen(okloci)))
            if not okloci:   # Check for reformatting
                locupdate = False
                for locus in snploci:
                    newlocus = string.join(string.split(locus,'.')[:-1])
                    if newlocus and newlocus in seqloci:
                        locupdate = True
                        for entry in snpdb.indexEntries('Locus',locus): entry['Locus'] = newlocus
                if locupdate:
                    snpdb.remakeKeys()
                    snploci = snpdb.indexKeys('Locus')
                    okloci = rje.listIntersect(seqloci,snploci)
                    self.printLog('#LOCUS','%s sequence loci; %s SNP loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(snploci),rje.iLen(okloci)))
            ftloci = ftdb.indexKeys('locus')
            okloci = rje.listIntersect(snploci,ftloci)
            self.printLog('#LOCUS','%s feature loci; %s SNP loci; %s in common.' % (rje.iLen(ftloci),rje.iLen(snploci),rje.iLen(okloci)))
            return True
        except: self.errorLog('Problem during %s setupSNPTable.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def featureSNPs(self):  ### Map loaded SNP and Feature tables
        '''Map loaded SNP and Feature tables'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db('SNP'); snpdb.index('Locus')
            ftdb = self.db('Feature'); ftdb.index('locus')
            seqlist = self.obj['SeqList']
            seqdict = seqlist.seqNameDic()
            #self.debug(seqdict.keys())
            mapdb = self.db().addEmptyTable('snpmap',snpdb.fields()+['Strand','GB']+ftdb.fields()+['SNPType','SNPEffect'],snpdb.keys()+['feature','start','end']) #['feature','start','end','protein_id','details',]
            snpdb.addField('Mapped',evalue=False)
            gb = self.obj['GenBank']
            ### ~ [2] Process one locus at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for locus in snpdb.indexKeys('Locus'):
                try: fullseq = seqlist.getSeq(seqdict[locus],format='tuple')[1]
                except: self.warnLog('Failed to get sequence for locus "%s": no SNP mapping.' % locus); continue
                if locus not in gb.dict['Sequence']: gb.dict['Sequence'][locus] = fullseq
                self.mapLocusSNPs(fullseq,locus)
            ### ~ [3] Summary files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AltPos'):
                lockeys = ['Locus','AltLocus']
                poskeys = ['Pos','AltPos']
            else:
                lockeys = ['Locus']
                poskeys = ['Pos']
            sumdb = self.db().copyTable(mapdb,'summary')
            sumfields = lockeys + ['feature','position','locus_tag','product']
            snptypes = ['NON','NS','SYN','SNP','INS','DEL']
            for snptype in snptypes:
                sumdb.addField(snptype,evalue=0)
                sumfields.append(snptype)
            for sentry in sumdb.entries(): sentry[sentry['SNPType']] = 1
            sumdb = self.db().copyTable(sumdb,'loci')
            sumdb.compress(sumfields[:4],default='sum')
            sumdb.keepFields(sumfields)
            cdsdb = self.db().copyTable(sumdb,'cds')
            cdsdb.remakeKeys()
            cdsdb.dropEntriesDirect('feature',['CDS'],inverse=True)
            cdsdb.compress(lockeys+['position'])
            cdsdb.dropFields(['feature','SNP'])
            cdsdb.addField('Ka',evalue=0)
            cdsdb.addField('Ks',evalue=0)
            for centry  in cdsdb.entries():
                centry['locus'] = centry['Locus']   # Needed for gb.featureSequence()
                # Frequency of synonymous subs
                fsyn = rje_sequence.sequenceKs(gb.featureSequence(centry),self)
                Ns = fsyn * len(gb.featureSequence(centry))
                Na = len(gb.featureSequence(centry)) - Ns
                centry['Ka'] = centry['NS'] / Na
                centry['Ks'] = centry['SYN'] / Ns
            cdsdb.makeField('Ka/Ks')
            cdsdb.saveToFile()
            sumdb.dropFields(['product'])
            sumdb.saveToFile()
            sumdb.setStr({'Name':'features'})
            sumdb.compress(sumfields[:3],default='sum')
            sumdb.dropFields(['position','locus_tag','product'])
            sumdb.saveToFile()
            sumdb = self.db('summary')
            sumdb.compress(lockeys+poskeys,default='max')
            for sentry in sumdb.entries():
                for i in range(len(snptypes)):
                    if sentry[snptypes[i]]:
                        for stype in snptypes[i+1:]: sentry[stype] = 0
                        break
            sumdb.compress(lockeys,default='sum')
            sumdb.keepFields(lockeys+snptypes)
            sumdb.saveToFile()
            ### ~ [4] Remove gene SNPs mapped to other features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Generate dictionary of mapped SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpft = {}
            for mentry in mapdb.entries():
                if mentry['feature'] not in self.list['FTBest']: continue   #!# Make case-insensitive? #!#
                snp = (mentry['Locus'],mentry['Pos'])
                if snp not in snpft: snpft[snp] = {}
                if mentry['feature'] not in snpft[snp]: snpft[snp][mentry['feature']] = []
                snpft[snp][mentry['feature']].append(mentry)
            ## ~ [4b] Process dictionary of mapped SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ftdrop = self.list['FTBest'][0:]
            dx = {}
            for ft in self.list['FTBest'][1:]: dx[ft] = 0
            while len(ftdrop) > 1:
                master = ftdrop.pop(0)
                for snp in snpft.keys():
                    if master not in snpft[snp]: continue
                    for ft in ftdrop:
                        if ft not in snpft[snp]: continue
                        for mentry in snpft[snp].pop(ft): mapdb.dropEntry(mentry); dx[ft] += 1
            for ft in self.list['FTBest'][1:]:
                self.printLog('#%s' % ft.upper()[:5],'%s mapped %s SNPs removed. (Mapped to other features).' % (rje.iStr(dx[ft]),ft))
            ### ~ [5] Finish and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapdb.saveToFile()
            #!# Add code for counts etc. Summarise by gene/feature. Do this before Reduction! #!#
        except: self.errorLog('%s.featureSNPs error' % self)
#########################################################################################################################
    def mapLocusSNPs(self,fullseq,locus):  ### Map loaded SNP and Feature tables
        '''Map loaded SNP and Feature tables'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gb = self.obj['GenBank']
            snpdb = self.db('SNP')
            #self.debug(snpdb.fields())
            skeys = snpdb.index('Locus')[locus][0:]
            ftdb = self.db('Feature')
            mapdb = self.db('snpmap')
            # MapDB Fields: snpdb.fields()+['Strand','GB']+ftdb.fields()+['SNPType','SNPEffect']
            # MapBD Keys: snpdb.keys()+['feature','start','end'])
            ### ~ [2] Map SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fx = 0  # No. failures
            skipx = 0   # No. skipped features
            for fkey in ftdb.index('locus')[locus]: # These should be in a consistent order: ['locus','start','feature','position']
                fentry = ftdb.data(fkey)
                self.debug(fentry)
                ## ~ [2a] Check for features to skip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if fentry['feature'] in self.list['FTSkip']: skipx += 1; continue
                #Fixed!# if 'join' in fentry['position']: self.warnLog('%s %s %s: Introns not handled' % (locus,fentry['feature'],fentry['position']))
                ## ~ [2b] Skip SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while skeys:
                    sentry = snpdb.data(skeys[0])
                    if (sentry['Pos'] + len(sentry['REF'])) <= fentry['start']: skeys.pop(0)  # No overlap (gone past)
                    else: break     # Caught up with features!
                ## ~ [2c] Process overlapping SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ftsnps = []     # Tuple of (Pos,REF,ALT)
                sx = 0
                for skey in skeys:
                    sentry = snpdb.data(skey)
                    if sentry['Pos'] > fentry['end']: break     # Gone beyond feature
                    ftsnps.append((sentry['Pos'],sentry['REF'],sentry['ALT']))
                    self.bugPrint('---\nFT: %s' % fentry)
                    self.bugPrint('SNP: %s' % sentry)
                    # NOTE: SNP Data should now be formatted either:
                    # (a) like RATT output with Pos, AltLocus, AltPos, REF and ALT
                    # (b) Like BCF with Pos, REF and ALT
                    try:
                        mentry = rje.combineDict({'Strand':'+','GB':fullseq[sentry['Pos']]},sentry)
                        mentry = rje.combineDict(mentry,fentry)
                        cseq = gb.featureSequence(fentry)
                        cpos = gb.featurePos(fentry,sentry['Pos'])  # Position 1-L
                        if 'join' in fentry['position']: self.debug(cpos)
                        if not cpos: raise ValueError('Feature outside position data!')
                        cpos -= 1   # Change to 0<L
                        mentry['length'] = len(cseq)
                        if sentry['REF'] == '.': mentry['GB'] = '-'
                        if cpos > 0:
                            try: mentry['GB'] = cseq[cpos]
                            except: print sentry, fentry, len(cseq), cpos; raise
                            if 'complement' in fentry['position']: mentry['Strand'] = '-'
                            elif sentry['REF'] == '.': cpos += 1
                        ### >>> OLD >>>
                        if not 'dealingwithintrons':
                            cseq = fullseq[fentry['start']-1:fentry['end']]
                            cpos = sentry['Pos'] - fentry['start']
                            #!# Add code to return cseq and cpos for join(...) loci #!#
                            #!# Will need to change feature to intron too if required #!#
                            mentry['length'] = len(cseq)
                            #?#rje_sequence.sequenceKs(sequence,self)
                            try: mentry['GB'] = cseq[cpos]
                            except: print sentry, fentry, len(cseq), cpos; raise
                            if sentry['REF'] == '.': cpos += 1; mentry['GB'] = '-'
                            if 'complement' in fentry['position']:
                                cseq = rje_sequence.reverseComplement(cseq)
                                cpos = len(cseq) - cpos - 1
                                mentry['Strand'] = '-'
                        ### <<< OLD <<<<
                        # Convert seq and cpos based on join()
                        if cpos < 1: mentry['SNPEffect'] = 'intron'
                        elif fentry['feature'] == 'CDS':
                            transl = rje.matchExp('/transl_table="?(\d+)"?',fentry['details'])
                            if transl: transl = transl[0]
                            else: transl = '1'
                            i = 3
                            while i <= cpos: i += 3
                            codon = cseq[i-3:i]
                            snp = mentry['ALT']
                            if snp in ['.','-']: snp = ''     # Deletion
                            if sentry['REF'] in ['.','-']:  # Insertion
                                if mentry['Strand'] == '-': snpseq = cseq[:cpos] + rje_sequence.reverseComplement(snp) + cseq[cpos:]
                                else: snpseq = cseq[:cpos] + snp + cseq[cpos:]
                            else:
                                if mentry['Strand'] == '-': snpseq = cseq[:cpos] + rje_sequence.reverseComplement(snp) + cseq[cpos+1:]
                                else: snpseq = cseq[:cpos] + snp + cseq[cpos+1:]
                            snpcodon = snpseq[i-3:i]
                            #print '\t', gtype, mentry[gtype], codon, '->', snpcodon
                            if mentry['SNPType'] in ['INS','DEL']:
                                mentry['SNPEffect'] = '%d:%s->%s' % ((cpos+2)/3,rje_sequence.dna2prot(cseq[i-3:],transl=transl),string.split(rje_sequence.dna2prot(snpseq[i-3:],transl=transl),'*')[0])
                            else:
                                mentry['SNPEffect'] = '%s%d%s' % (rje_sequence.dna2prot(codon,transl=transl),(cpos+2)/3,rje_sequence.dna2prot(snpcodon,transl=transl))
                                # Update SNPType
                                if mentry['SNPEffect'][0] == mentry['SNPEffect'][-1]: mentry['SNPType'] = 'SYN'
                                elif mentry['SNPEffect'][-1] == '*': mentry['SNPType'] = 'NON'
                                else: mentry['SNPType'] = 'NS'
                        else: mentry['SNPEffect'] = fentry['feature']
                        #self.deBug(mentry)
                        mapdb.addEntry(mentry); sx += 1
                        sentry['Mapped'] = True
                    except:
                        self.errorLog('Bugger')
                        self.warnLog('Failed: %s' % sentry); fx += 1
                        mentry = rje.combineDict({'feature':'FAILED'},sentry)
                        mapdb.addEntry(mentry)
                        self.debug(mentry)
                ## ~ [2d] Combine all SNPs in ftsnps for new sequence output ~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Add code #!#

            ### ~ [3] Process intergenic SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ix = 0
            for sentry in snpdb.indexEntries('Locus',locus):
                if not sentry['Mapped']:
                    mentry = rje.combineDict({'Strand':'na'},sentry)
                    mentry = rje.combineDict(mentry,{'GB':fullseq[sentry['Pos']-1],'feature':'intergenic'})
                    mapdb.addEntry(mentry); ix += 1
            self.printLog('#SNPMAP','%s of %s %s SNPs mapped to intergenic regions only.' % (rje.iStr(ix),rje.iLen(snpdb.indexEntries('Locus',locus)),locus))
            if fx: self.printLog('#FAIL','%s %s SNPs failed to be parsed/mapped.' % (rje.iStr(fx),locus))
            else: self.printLog('#SNPMAP','%s %s SNPs failed to be parsed/mapped.' % (rje.iStr(fx),locus))
        except: self.errorLog('%s.mapLocusSNPs(%s) error' % (self,locus))
#########################################################################################################################
    ### <4> ### Old SNPMap Methods to revise                                                                            #
#########################################################################################################################
    def oldRun(self):  ### Main controlling run method
        '''Main controlling run method'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batchfiles = []
            for batch in self.list['Batch']: batchfiles += glob.glob(batch)
            self.log.printLog('#FILES','%d files identified for batch run' % len(batchfiles))

            ### ~ [2] Run each batch file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for file in batchfiles: self.processGene(file)
            self.log.printLog('#ZEN',rje_zen.Zen().wisdom())
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def processGene(self,file):  ### Main controlling run method for a single gene
        '''Main controlling run method for a single gene.'''
        try:### ~ [1] Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if string.count(file,'.') > 1: return
            seqlist = self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % file])
            #!# Replace with unaligned sequences and list of exon positions #!#
            gene = seqlist.info['Basefile']
            self.log.printLog('#GENE','Processing files for %s gene' % gene)
            ## ~ [1b] Regenerate alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            exons = self.loadFromFile('%s_exons.txt' % gene,chomplines=True)
            if not exons: return self.log.errorLog('Missing exon data for %s' % gene,printerror=False)
            cds = string.replace(seqlist.seq[1].info['Sequence'],'-','')
            aln = '-' * len(seqlist.seq[0].info['Sequence'])
            for exon in exons:
                exon = string.split(exon)
                (start,end) = (int(exon[0]),int(exon[1]))
                elen = end - start + 1
                eseq = cds[:elen]
                if len(eseq) != elen: return self.log.errorLog('Not enough CDS for exon (%s,%s) - maybe be indels in sequence?' % (exon[0],exon[1]),printerror=False)
                cds = cds[elen:]
                aln = aln[:start-1] + eseq + aln[end:]
                if len(aln) != len(seqlist.seq[0].info['Sequence']): raise ValueError
            seqlist.seq[1].info['Sequence'] = aln[0:]
            seqlist.saveFasta(seqfile='%s.aln.fas' % gene)

            ### ~ [2] Map SNPs onto CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdata = rje.dataDict(self,'%s_polymorphisms.txt' % gene,mainkeys=['SNP'],delimit='\t')
            self.log.printLog('#SNP','%d SNPs read from %s_polymorphisms.txt' % (len(snpdata),gene))
            ## ~ [2a] Find ATG start site in consensus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            consensus = seqlist.seq[0].info['Sequence']
            cds = seqlist.seq[1].info['Sequence']
            dna = string.replace(cds,'-','')
            open('%s.protein.txt' % gene,'w').write(rje_sequence.dna2prot(dna))     #!# Tranlsation tables! #!#
            atg = 0     # Position in consensus that matches ATG of CDS  (0<L)
            while cds[atg] == '-': atg += 1
            ## ~ [2b] Map SNPs relative to ATG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdict = {}    # List of new SNP data dictionaries
            for snp in snpdata:
                rje.combineDict(snpdata[snp],{'SNP':int(snp),'CDS':int(snp)-atg,'ConSeq':consensus[int(snp)-1],'CDSSeq':cds[int(snp)-1]})
                snpdict[int(snp)] = snpdata[snp]
            ## ~ [2c] Extract protein information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for snpkey in rje.sortKeys(snpdict):
                snp = snpdict[snpkey]
                if snp['CDSSeq'] == '-': continue
                if snp['CDSSeq'] not in [snp['MinAllele'],snp['MajAllele']]: self.log.errorLog('SNP %s CDS "%s" matches neither allele' % (snp['SNP'],snp['CDSSeq']),printerror=False)
                myseq = string.replace(cds[atg:snp['SNP']],'-','')
                if myseq[-1] != snp['CDSSeq']: print myseq, '!=', snp['CDSSeq']
                snp['ProteinPos'] = ((len(myseq)-1) / 3) + 1
                frame = snp['CodonPos'] = [3,1,2][len(myseq) % 3]
                codon = snp['Codon'] = dna[(snp['ProteinPos']-1)*3:snp['ProteinPos']*3]
                snp['MajAA'] = rje_sequence.dna2prot(codon[:frame-1]+snp['MajAllele']+codon[frame:])
                snp['MinAA'] = rje_sequence.dna2prot(codon[:frame-1]+snp['MinAllele']+codon[frame:])
                snp['Nonsyn'] = snp['MajAA'] != snp['MinAA']
                self.deBug(snp)
                self.deBug(myseq)
                self.deBug(rje_sequence.dna2prot(myseq))
                self.deBug(dna[:snp['ProteinPos']*3])
                self.deBug(rje_sequence.dna2prot(dna[:snp['ProteinPos']*3]))

            ### ~ [3] Output new data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.snpmap.tdt' % gene
            headers = ['SNP','MajAllele','MinAllele','CDS','ConSeq','CDSSeq','ProteinPos','Codon','CodonPos','MajAA','MinAA','Nonsyn']
            rje.delimitedFileOutput(self,outfile,headers,'\t',rje_backup=True)
            for snp in rje.sortKeys(snpdict): rje.delimitedFileOutput(self,outfile,headers,'\t',snpdict[snp])
                
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def OLDsetup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Load SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            delim = rje.delimitFromExt(filename=self.getStr('SNPFile'))
            skeys = self.list['SNPKeys']
            if self.list['SNPHead'] and 'Locus' in self.list['SNPHead']: skeys = ['Locus']
            elif not self.list['SNPHead'] and 'Locus' in rje.readDelimit(open(self.getStr('SNPFile'),'r').readline(),delim): skeys = ['Locus']
            skeys.append('Pos')
            snpdb = self.db().addTable(self.getStr('SNPFile'),mainkeys=skeys,headers=self.list['SNPHead'],name='SNP')
            for field in snpdb.fields():
                if field.lower() == 'ref' and field != 'REF': snpdb.renameField(field,'REF')
            if '' in snpdb.fields(): snpdb.dropField('')
            if 'EXTRA' not in snpdb.fields(): snpdb.addField('EXTRA')
            snpdb.dataFormat({'Pos':'int'})
            for gtype in self.list['Genotypes'][0:]:
                if gtype not in snpdb.fields():
                    self.printLog('#ERR','Could not find "%s" in %s' % (gtype,self.getStr('SNPFile')))
                    self.list['Genotypes'].remove(gtype)
            self.printLog('#GTYPE','%d genotypes to map SNPs' % len(self.list['Genotypes']))
            ## ~ [1a] Check ScreenMatch versus Genotypes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for g1 in self.list['ScreenMatch'][0:]:
                if g1 not in self.list['Genotypes']:
                    self.printLog('#GTYPE','ScreenMatch genotype "%s" not recognised.' % g1)
                    self.list['ScreenMatch'].remove(g1)
            if self.list['ScreenMatch']:
                sx = snpdb.entryNum()
                for sentry in snpdb.entries()[0:]:
                    screened = True
                    for g1 in self.list['ScreenMatch'][:-1]:
                        for g2 in self.list['ScreenMatch'][1:]:
                            if sentry[g1] != sentry[g2]: screened = False
                    if screened: snpdb.dropEntry(sentry)
                if snpdb.entryNum() != sx: self.printLog('#SNP','%d SNPs screened (matching %s SNPs)' % (sx - snpdb.entryNum(),string.join(self.list['ScreenMatch'],',')))
            ### ~ [2] Load GenBank ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gb = self.obj['GenBank']
            if 'full' not in gb.list['FasOut']: gb.list['FasOut'].append('full')
            if not gb.loadFeatures() or not os.path.exists('%s.full.fas' % gb.baseFile()): gb.run()
            self.obj['SeqList'].loadSeq('%s.full.fas' % gb.baseFile())
            self.obj['SeqList'].nextSeq()
            ftdb = self.db('Feature')
            ftdb.dataFormat({'start':'int','end':'int'})
            self.db().baseFile(self.baseFile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def OLDfeatureSNPs(self):  ### Use loaded SNP and Feature table
        '''
        Use loaded SNP and Feature table.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Change this to work on one locus at a time and give method split SNP and FT tables?
            #!# Or just index SNP and FT tables on locus... (Add locus if not given when loading.)
            snpdb = self.db('SNP')
            for sentry in snpdb.entries():
                for gtype in self.list['Genotypes']:
                    if gtype != 'ref' and sentry[gtype] == '-': sentry[gtype] = sentry['REF']
            ftdb = self.db('Feature')
            mapdb = self.db().addEmptyTable('snpmap',snpdb.fields()+['Strand','GB']+ftdb.fields(),snpdb.keys()+['feature','start','end']) #['feature','start','end','protein_id','details',]
            if 'Locus' in mapdb.fields() and 'locus' in mapdb.fields(): mapdb.dropField('locus')
            for gtype in self.list['Genotypes']: mapdb.addField('%s_Effect' % gtype)
            #fullseq = self.obj['SeqList'].getSeq(format='tuple')[1]     # Only 1 sequence ! #
            seqdict = self.obj['SeqList'].makeSeqNameDic('accnum')
            ### ~ [2] Map SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fx = 0  # No. failures
            for sentry in snpdb.entries():
                if 'Locus' in sentry: fullseq = self.obj['SeqList'].getSeq(seqdict[sentry['Locus']],format='tuple')[1]
                else: fullseq = self.obj['SeqList'].getSeq(format='tuple')[1]
                try:
                    #??? What is this REF stuff about?  ???#
                    if 'REF' in sentry and 'ALT' in sentry:
                        if len(sentry['REF']) > len(sentry['ALT']):
                            if sentry['REF'] == fullseq[sentry['Pos']-1:][:len(sentry['REF'])]:
                                sentry['EXTRA'] = 'Deletion'
                            elif sentry['ALT'] == fullseq[sentry['Pos']-1:][:len(sentry['ALT'])]:
                                sentry['EXTRA'] = 'WT-Insertion'
                            else:
                                sentry['EXTRA'] = 'REF/ALT Mismatch'
                                self.deBug(sentry['REF'])
                                self.deBug(fullseq[sentry['Pos']-1:][:len(sentry['REF'])])
                                self.deBug(sentry['ALT'])
                                self.deBug(fullseq[sentry['Pos']-1:].find(sentry['REF']))
                        elif len(sentry['REF']) < len(sentry['ALT']):
                            if sentry['ALT'] == fullseq[sentry['Pos']-1:][:len(sentry['ALT'])]:
                                sentry['EXTRA'] = 'WT-Deletion'
                            elif sentry['REF'] == fullseq[sentry['Pos']-1:][:len(sentry['REF'])]:
                                sentry['EXTRA'] = 'Insertion'
                            else:
                                sentry['EXTRA'] = 'ALT/REF Mismatch'
                                self.deBug(sentry['ALT'])
                                self.deBug(fullseq[sentry['Pos']-1:][:len(sentry['ALT'])])
                                self.deBug(sentry['REF'])
                                self.deBug(fullseq[sentry['Pos']-1:].find(sentry['REF']))
                    mapped = False
                    for fentry in ftdb.entries():
                        if fentry['feature'] in self.list['FTSkip']: continue
                        if 'Locus' in sentry and fentry['locus'] != sentry['Locus']: continue
                        if fentry['start'] <= sentry['Pos'] <= fentry['end']:
                            mentry = rje.combineDict({'Strand':'+'},sentry)
                            mentry = rje.combineDict(mentry,fentry)
                            cseq = fullseq[fentry['start']-1:fentry['end']]
                            cpos = sentry['Pos'] - fentry['start']
                            if 'complement' in fentry['position']:
                                cseq = rje_sequence.reverseComplement(cseq)
                                cpos = len(cseq) - cpos - 1
                                mentry['Strand'] = '-'
                            mentry['GB'] = cseq[cpos]
                            if fentry['feature'] == 'CDS':
                                i = 3
                                while i <= cpos: i += 3
                                codon = cseq[i-3:i]
                                #print sentry['ref'], mentry['GB'], codon
                                for gtype in self.list['Genotypes']:
                                    if mentry[gtype] in ['ins','del']:
                                        mentry['%s_Effect' % gtype] = 'frameshift'
                                        #!# Not sure if this is OK #!#
                                        if mentry[gtype] == 'ins': mentry[gtype] = mentry['GB'] + 'N'
                                        if mentry[gtype] == 'del': mentry[gtype] = ''
                                        snpseq = cseq[:cpos] + mentry[gtype] + cseq[cpos+1:]
                                        mentry['%s_Effect' % gtype] = '%d:%s->%s' % ((cpos+2)/3,rje_sequence.dna2prot(cseq[i-3:]),string.split(rje_sequence.dna2prot(snpseq[i-3:]),'*')[0])
                                        #!# ^^^^^^^^^^^^^^^^^^^^^^ #!#
                                    elif len(mentry[gtype]) > 1 and 'REF' in mentry:    # REF vs ALT
                                        #if len(mentry['REF']) > len(mentry['ALT']):
                                        #!# Modify to match above #!#
                                        mentry['%s_Effect' % gtype] = 'indel%d' % rje.modulus(len(mentry['REF']) - len(mentry['ALT']))
                                    elif mentry[gtype] == '.' and 'REF' in mentry:    # Deletion
                                        mentry[gtype] = ''
                                        #!# Not sure if this is right? #!#
                                        snpseq = cseq[:cpos] + mentry[gtype] + cseq[cpos+1:]
                                        mentry['%s_Effect' % gtype] = '%d:%s->%s' % ((cpos+2)/3,rje_sequence.dna2prot(cseq[i-3:]),string.split(rje_sequence.dna2prot(snpseq[i-3:]),'*')[0])
                                    elif 'REF' in mentry and mentry['REF'] == '.':    # Insertion
                                        mentry[gtype] = mentry['GB'] + mentry[gtype]
                                        mentry['REF'] = mentry['GB']
                                        #!# Not sure if this is right? #!#
                                        snpseq = cseq[:cpos] + mentry[gtype] + cseq[cpos+1:]
                                        mentry['%s_Effect' % gtype] = '%d:%s->%s' % ((cpos+2)/3,rje_sequence.dna2prot(cseq[i-3:]),string.split(rje_sequence.dna2prot(snpseq[i-3:]),'*')[0])
                                    else:
                                        if mentry['Strand'] == '-': snpseq = cseq[:cpos] + rje_sequence.reverseComplement(mentry[gtype]) + cseq[cpos+1:]
                                        else: snpseq = cseq[:cpos] + mentry[gtype] + cseq[cpos+1:]
                                        snpcodon = snpseq[i-3:i]
                                        #print '\t', gtype, mentry[gtype], codon, '->', snpcodon
                                        mentry['%s_Effect' % gtype] = '%s%d%s' % (rje_sequence.dna2prot(codon),(cpos+2)/3,rje_sequence.dna2prot(snpcodon))
                            else:
                                for gtype in self.list['Genotypes']: mentry['%s_Effect' % gtype] = fentry['feature']
                            #self.deBug(mentry)
                            mapdb.addEntry(mentry)
                            mapped = True
                    if not mapped:
                        mentry = rje.combineDict({'Strand':'na'},sentry)
                        mentry = rje.combineDict(mentry,{'GB':fullseq[sentry['Pos']-1],'feature':'intergenic'})
                        mapdb.addEntry(mentry); fx += 1
                except:
                    self.warnLog('Failed: %s' % sentry); fx += 1
                    mentry = rje.combineDict({'feature':'FAILED'},sentry)
                    mapdb.addEntry(mentry)
            ### ~ [3] Remove gene SNPs mapped to other features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'gene' in mapdb.index('feature'):
                #!# Incorporate locus properly #!#
                #!# Then repeat with mRNA to remove CDS annotation #!#
                #!# Ditto mobile_element #!#
                gx = 0
                for mkey in mapdb.index('feature')['gene'][0:]:
                    mentry = mapdb.data(mkey)
                    if len(mapdb.indexDataList('Pos',mentry['Pos'],'feature')) > 1: mapdb.dropEntry(mentry); gx += 1
                self.printLog('#GENE','%d mapped gene features removed. (Mapped to other features.)' % gx)
            mapdb.saveToFile()
            self.printLog('#FAIL','%s SNPs failed to be parsed/mapped.' % rje.iStr(fx))
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################
### End of SECTION II: SNPMap Class                                                                                     #
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
    try: SNPMap(mainlog,cmd_list).run()

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
