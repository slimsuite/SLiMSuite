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
Version:      1.2.1
Last Edit:    03/05/21
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

    The default FTBest hierarchy for `*.ftypes.tdt` output is:
        CDS,mRNA,tRNA,rRNA,ncRNA,misc_RNA,gene,mobile_element,LTR,rep_origin,telomere,centromere,misc_feature,intergenic

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
    ftskip=LIST     : List of feature types to exclude from analysis [source]
    ftbest=LIST     : List of features to exclude if earlier feature in list overlaps position [(see above)]
    snpbyftype=T/F  : Whether to output mapped SNPs by feature type (before FTBest filtering) [False]
    snpfile=FILE    : Input table of SNPs to map and output (should have locus and pos info, see above) []
    snphead=LIST    : List of SNP file headers []
    snpdrop=LIST    : List of SNP fields to drop []
    altpos=T/F      : Whether SNP file is a single mapping (with AltPos) (False=BCF) [True]
    altft=T/F       : Use AltLocus and AltPos for feature mapping (if altpos=T) [False]
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
    # 0.5.0 - Added CDS rating.
    # 0.6.0 - Added AltFT mapping mode (map features to AltLocus and AltPos)
    # 0.7.0 - Added additional fields for processing Snapper output. (Hopefully will still work for SAMTools etc.)
    # 0.8.0 - Added parsing of GFF file from Prokka.
    # 0.8.1 - Corrected "intron" classification for first position of features. Updated FTBest defaults.
    # 1.0.0 - Version that works with Snapper V1.0.0. Not really designed for standalone running any more.
    # 1.1.0 - Added pNS and modified the "Positive" CDS rating to be pNS < 0.05.
    # 1.1.1 - Updated pNS calculation to include EXT mutations and substitution frequency.
    # 1.2.0 - SNPByFType=T/F  : Whether to output mapped SNPs by feature type (before FTBest filtering) [False]
    # 1.2.1 - Fixed GFF parsing bug.
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
    # [Y] : Deal with introns!
    # [ ] : Add syn/NS/nonsense rating to aa changes (extra field)
    # [ ] : Update input to be able to handle a sequence file and features file.
    # [ ] : Add a summary of the SNP types for each Locus[/SNPKeys].
    # [ ] : If required:     snpalt=X        : SNP table alternative genotype field [ALT]
    # [ ] : In addition to any annotation extracted using rje_genbank, additional features can be given to be mapped on to the genomic sequence?
    # [ ] : Or an option to merge results from several runs (e.g. with different features?)
    # [ ] : Add sequence output for different features.
    # [ ] : Add ratings for CDS based on Ka/Ks and no. deletions etc.
    # [ ] : Add output of possible duplications (and total deletions?) - May need to compare to self-analysis to be sure.
    # [ ] : Add blank entries for CDS without any SNPs. (Autofill the altlocus if flanks found?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('SNP_MAPPER', '1.2.3', 'May 2021', '2015')
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class SNPMap(rje_obj.RJE_Object):     
    '''
    SNPMap Class. Author: Rich Edwards (2013).

    Str:str
    - FTFile=FILE     : Input feature file (locus,feature,position,start,end) [*.Feature.tdt]
    - GenBank = Genbank file from which Features/Sequences generated [None]
    - SeqIn = Sequence input file with accession numbers matching Locus IDs, or Genbank file.  []
    - SNPFile = Input table of SNPs to map and output (should have locus and pos info, see above) []
  
    Bool:boolean
    - AltFT=T/F       : Use AltLocus and AltPos for feature mapping (if altpos=T) [False]
    - AltPos=T/F      : Whether SNP file is a single mapping (with AltPos) (False=BCF) [True]
    - SNPByFType=T/F  : Whether to output mapped SNPs by feature type (before FTBest filtering) [False]

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
        self.strlist = ['FTFile','GenBank','SeqIn','SNPFile']
        self.boollist = ['AltFT','AltPos','SNPByFType']
        self.intlist = []
        self.numlist = []
        self.listlist = ['FTBest','FTSkip','SNPDrop','SNPHead']
        self.dictlist = []
        self.objlist = ['DB','GenBank','SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({})
        self.setBool({'AltPos':True,'AltFT':False,'SNPByFType':False})
        self.setInt({})
        self.setNum({})
        self.list['FTSkip'] = ['source']
        self.list['FTBest'] = rje.split('CDS,mRNA,tRNA,rRNA,ncRNA,misc_RNA,gene,mobile_element,LTR,rep_origin,telomere,centromere,misc_feature,intergenic',',')
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
                self._cmdReadList(cmd,'file',['FTFile','GenBase','SeqIn','SNPFile'])  # String representing file path
                self._cmdReadList(cmd,'bool',['AltFT','AltPos','SNPByFType'])  # True/False Booleans
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
        if self.getBool('AltFT') and not self.getBool('AltPos'):
            self.printLog('#CMD','Cannot have altft=True if altpos=False. AltFT=False.')
            self.setBool({'AltFT':False})
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
            if not self.setupSNPTable(): return False
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
            #!# Replace with rje_genbank code? #!#
            if open(seqin,'r').readline()[:1] != '>': gbfile = seqin
            elif rje.exists('%s.gb' % seqin): gbfile = '%s.gb' % seqin
            elif rje.exists('%s.gbk' % seqin): gbfile = '%s.gbk' % seqin
            else: gbfile = '%s.gb' %  rje.baseFile(seqin)
            #x#if open(seqin,'r').readline()[:1] != '>':   # Not fasta file: assume Genbank. Should have features file.
            if rje.exists(gbfile):   # Not fasta file: assume Genbank. Should have features file.
                self.setStr({'GenBank':gbfile})
                gb.setStr({'SeqIn':gbfile})
                gbase = rje.baseFile(gbfile)
                #if gbase.endswith('.full'): gbase = rje.baseFile(gbase)
                gb.baseFile(gbase)
                if 'full' not in gb.list['FasOut']: gb.list['FasOut'].append('full')
                #x#if not gb.loadFeatures() or not os.path.exists('%s.full.fas' % gb.baseFile()): gb.run()
                if not gb.loadFeatures() or not os.path.exists(seqin): gb.run()
                if not rje.exists(seqin): seqin = '%s.full.fas' % gb.baseFile()
                gbftfile = '%s.Feature.tdt' % gbase
                gb.db().deleteTable('Feature')   # Want to reload from table
                self.db().addTable(gbftfile,name='Feature',expect=True,mainkeys=['locus','feature','position'])
            else: self.printLog('#GB','%s not found: will process with seqin=FILE and ftfile=FILE' % gbfile)
            self.setStr({'SeqIn':seqin})
            seqlist.loadSeq(seqin)
            seqdict = seqlist.makeSeqNameDic('max')
            ### ~ [2] Load Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.db('Feature') and (self.getStr('FTFile') == gbftfile or not self.getStrLC('FTFile')):
                self.printLog('#FT','Genbank FTFile "%s" already loaded.' % gbftfile)
            elif not rje.exists(self.getStr('FTFile')):
                self.warnLog('Cannot find ftfile=%s!' % self.getStr('FTFile'))
            elif self.db('Feature'):    # Add features
                ftdb = self.db().addTable(self.getStr('FTFile'),name='FTFile',expect=True,mainkeys=['locus','feature','position'])
                self.db().mergeTables(self.db('Feature'),ftdb)
            elif open(self.getStr('FTFile'),'r').readline().startswith('##gff-'):
                baddetails = 0
                protected = ['locus','feature','position','start','end']
                ftdb = self.db().addEmptyTable('Feature',['locus','feature','position','start','end','product','gene_synonym','note','db_xref','locus_tag','details'],['locus','feature','position'])
                for gline in open(self.getStr('FTFile'),'r').readlines():
                    if gline.startswith('##FASTA'): break
                    elif not gline or gline.startswith('#'): continue
                    gdata = rje.readDelimit(gline)
                    gentry = {'locus':gdata[0],'feature':gdata[2],'start':gdata[3],'end':gdata[4],
                              'product':'','gene_synonym':'','note':'','db_xref':'','locus_tag':'','details':''}
                    if gdata[6] == '+': gentry['position'] = '%s..%s' % (gdata[3],gdata[4])
                    else: gentry['position'] = 'complement(%s..%s)' % (gdata[3],gdata[4])
                    ginfo = rje.split(gdata[8],';')
                    for detail in ginfo[0:]:
                        try:
                            [dkey,dvalue] = rje.split(detail,'=',maxsplit=1)
                            if dkey in protected: continue
                            if dkey in gentry:
                                gentry[dkey] = dvalue
                                ginfo.remove(detail)
                            elif dkey in ['Parent','ID','Name']:
                                gentry[dkey] = dvalue
                        except:
                            baddetails += 1
                    if not gentry['locus_tag']:
                        for att in ['Parent','ID','Name']:
                            if att in gentry and gentry[att]:
                                gentry['locus_tag'] = gentry[att]
                                break
                    gentry['details'] = rje.join(ginfo,'; ')
                    ftdb.addEntry(gentry)
                self.printLog('#FTFILE','%s features parsed from GFF %s' % (rje.iStr(ftdb.entryNum()),self.getStr('FTFile')))
                if baddetails: self.warnLog('{0} details had issues with parsing. (Lacking "=" or ";" inside quotes.)'.format(rje.iStr(baddetails)))
                newftfile = '%s.Feature.tdt' % rje.baseFile(seqin,strip_path=True)
                ftdb.saveToFile(newftfile)
            else: self.db().addTable(self.getStr('FTFile'),name='Feature',expect=True,mainkeys=['locus','feature','position'])
            ftdb = self.db('Feature')
            if ftdb:
                ftdb.dropEntriesDirect('feature',self.list['FTSkip'])
                ftdb.dataFormat({'start':'int','end':'int'})
                ftdb.newKey(['locus','start','feature','position'])    # Reorder for sorting
            #self.debug(ftdb.dataKeys()[:10])
            #self.debug(ftdb.entries()[:10])
            ### ~ [3] Check sequences and loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                seqloci = rje.sortKeys(seqdict)
                ftloci = ftdb.indexKeys('locus')
                okloci = rje.listIntersect(seqloci,ftloci)
                self.printLog('#LOCUS','%s sequence keys; %s feature loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(ftloci),rje.iLen(okloci)))
            else: self.printLog('#LOCUS','No feature file: no Locus features to map.')
            self.printLog('#SETUP','Reference setup complete!')
            return True
        except: self.errorLog('Problem during %s setupReference.' % self); return False  # Setup failed
#########################################################################################################################
    def setupSNPTable(self,keepmatches=False):    ### Loads, reformats and checks SNP Table
        '''
        Loads, reformats and checks SNP Table.
        >> keepmatches:bool [False] = Whether to keep REF/ALT matches. (Want to do this if FT ends added by Snapper.)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('SNPFile')): raise IOError('Cannot find snpfile=%s!' % self.getStr('SNPFile'))
            if not self.getStrLC('Basefile'):
                self.baseFile(rje.baseFile(self.getStr('SNPFile'),strip_path=True))
            if self.getBool('AltPos'):
                skeys = ['Locus','Pos','AltLocus','AltPos']
            else:
                skeys = ['Locus','Pos','ALT']
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
                if sentry['REF'] == sentry['ALT']: stype = 'ID'   # No mismatch!
                elif sentry['REF'] in ['.','-'] or len(sentry['REF']) < len(sentry['ALT']): stype = 'INS'
                elif sentry['ALT'] in ['.','-'] or len(sentry['REF']) > len(sentry['ALT']): stype = 'DEL'
                else: stype = 'SNP'
                sentry['SNPType'] = stype
            if not keepmatches: snpdb.dropEntriesDirect('SNPType',['ID'])
            if not snpdb.entries(): self.printLog('#SNP','No SNPs for analysis!'); return False
            ### ~ [4] Check and rename loci ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdb = self.db('Feature')
            seqloci = self.obj['SeqList'].seqNameDic().keys()
            if self.getBool('AltFT'): locfield = 'AltLocus'
            else: locfield = 'Locus'
            snploci = snpdb.indexKeys(locfield)
            okloci = rje.listIntersect(seqloci,snploci)
            self.printLog('#LOCUS','%s sequence loci; %s SNP loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(snploci),rje.iLen(okloci)))
            if not okloci:   # Check for reformatting
                # Add extra field with [Alt]LocusName if Locus changed?
                snpdb.addField('%sName' % locfield)
                locupdate = False
                for locus in snploci:
                    newlocus = rje.join(rje.split(locus,'.')[:-1])
                    if newlocus and newlocus in seqloci:
                        locupdate = True
                        for entry in snpdb.indexEntries(locfield,locus):
                            entry['%sName' % locfield] = locus
                            entry[locfield] = newlocus
                if locupdate:
                    snpdb.remakeKeys()
                    snploci = snpdb.indexKeys('Locus')
                    okloci = rje.listIntersect(seqloci,snploci)
                    self.printLog('#LOCUS','%s sequence loci; %s SNP loci; %s in common.' % (rje.iLen(seqloci),rje.iLen(snploci),rje.iLen(okloci)))
            if ftdb:
                ftloci = ftdb.indexKeys('locus')
                okloci = rje.listIntersect(snploci,ftloci)
                self.printLog('#LOCUS','%s feature loci; %s SNP loci; %s in common.' % (rje.iLen(ftloci),rje.iLen(snploci),rje.iLen(okloci)))
            else: self.printLog('#LOCUS','No Feature table: no Locus features to map.')
            return True
        except: self.errorLog('Problem during %s setupSNPTable.' % self); raise  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def featureSNPs(self):  ### Map loaded SNP and Feature tables
        '''Map loaded SNP and Feature tables'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db('SNP'); snpdb.index('Locus')
            ftdb = self.db('Feature')
            if not ftdb:
                self.printLog('#LOCUS','No Feature table: no featuresSNP output.')
                return False
            ftdb.index('locus')
            seqlist = self.obj['SeqList']
            seqdict = seqlist.seqNameDic()
            #self.debug(seqdict.keys())
            snpfields = snpdb.fields()
            if 'SNPType' in snpfields: snpfields.remove('SNPType')
            mapdb = self.db().addEmptyTable('snpmap',snpfields+['Strand','GB']+ftdb.fields()+['SNPType','SNPEffect'],list(snpdb.keys())+['feature','start','end']) #['feature','start','end','protein_id','details',]
            snpdb.addField('Mapped',evalue=False)
            gb = self.obj['GenBank']
            #!# Need to improve the fsyn calculation based on mutation frequencies
            mutdict = {}   #i# (n1,n2):count
            self.progLog('#SUBDIC','Generating substitution dictionary...')
            for entry in snpdb.entries():
                n1 = entry['REF']
                n2 = entry['ALT']
                nkey = (n1,n2)
                if nkey in mutdict: mutdict[nkey] += 1
                else: mutdict[nkey] = 1
            self.printLog('\r#SUBDIC','Generation of substitution dictionary complete.')
            mutations = []
            for nkey in mutdict:
                mutations.append((nkey[0],nkey[1],mutdict[nkey]))
            mutdict = rje_sequence.mutationDict(mutations,gaps=False,nosub=False,callobj=self)   #i# (n1,n2,count)
            self.deBug(mutdict)
            ### ~ [2] Process one locus at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AltFT'): snploci = snpdb.indexKeys('AltLocus')
            else: snploci = snpdb.indexKeys('Locus')
            #self.debug(snploci)
            for locus in snploci:
                if locus not in ftdb.index('locus'): ftlocus = rje.split(locus,'_')[-1]
                else: ftlocus = locus
                try: fullseq = seqlist.getSeq(seqdict[locus],format='tuple')[1]
                except: self.warnLog('Failed to get sequence for locus "%s": no SNP mapping.' % locus); continue
                if ftlocus not in gb.dict['Sequence']: gb.dict['Sequence'][ftlocus] = fullseq
                self.mapLocusSNPs(fullseq,locus)
            ### ~ [3] Summary files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AltPos'):
                lockeys = ['Locus','AltLocus']
                poskeys = ['Pos','AltPos']
            else:
                lockeys = ['Locus']
                poskeys = ['Pos']
            ## ~ [3a] *.summary.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # NB. This is set up here but finished and saved later
            sumdb = self.db().copyTable(mapdb,'summary')
            sumfields = lockeys + ['feature','position','locus_tag','product','details','RefSNP','AltQV','AltCN','RefCN']
            snptypes = ['ID','NON','EXT','NS','SYN','SNP','INS','DEL']
            for snptype in snptypes:
                sumdb.addField(snptype,evalue=0)
                sumfields.append(snptype)
            for sentry in sumdb.entries(): sentry[sentry['SNPType']] = 1
            ## ~ [3b] *.loci.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Summary per locus
            # NB. This is set up here but finished and saved later
            sumdb = self.db().copyTable(sumdb,'features')
            compfields = sumfields[:4]
            if 'RefSNP' in sumdb.fields(): compfields.append('RefSNP')
            sumdb.compress(compfields,default='sum',rules={'RefSNP':'mean','AltCN':'mean','RefCN':'mean','MeanCNV':'mean','MaxCNV':'max','AltQV':'mean'})
            sumdb.keepFields(sumfields)
            ## ~ [3c] *.cds.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Summary per CDS features
            cdsdb = self.db().copyTable(sumdb,'cds')
            cdsdb.remakeKeys()
            cdsdb.dropEntriesDirect('feature',['CDS'],inverse=True)
            compfields = lockeys+['position']
            if 'RefSNP' in cdsdb.fields(): compfields.append('RefSNP')
            cdsdb.compress(compfields)
            #!# Make sure that SNPs in the same position are not covered twice?!
            #i# Will now be once for the "Reference" match and one mean for all others
            #!# Actually using the mean for a CDS, not the sum!
            cdsdb.dropFields(['feature','SNP'])
            cdsdb.addField('length',evalue=0)
            cdsdb.addField('Ka',evalue=0)
            cdsdb.addField('Ks',evalue=0)
            cdsdb.addField('Ka/Ks',evalue=0.0)
            cdsdb.addField('pSYN',evalue=1.0)
            cdsdb.addField('pNS',evalue=1.0)
            matchesallowed = False
            for centry  in cdsdb.entries():
                centry['locus'] = centry['Locus']   # Needed for gb.featureSequence()
                matchesallowed = matchesallowed or centry['ID'] > 0
                # Frequency of synonymous subs
                ftseq = gb.featureSequence(centry)
                centry['length'] = len(ftseq)
                fsyn = rje_sequence.sequenceKs(ftseq,self,mutdict=mutdict)
                Ns = fsyn * len(ftseq)
                Na = len(ftseq) - Ns
                self.bugPrint('FTSeq: %s' % ftseq)
                self.bugPrint('FSyn: %s' % fsyn)
                self.bugPrint('Ns: %s' % Ns)
                self.bugPrint('Na: %s' % Na)
                centry['Ka'] = centry['NS'] / Na
                centry['Ks'] = centry['SYN'] / Ns
                if centry['Ks']: centry['Ka/Ks'] = centry['Ka'] / centry['Ks']
                else: centry['Ka/Ks'] = centry['Ka'] / (0.5 / Ns)
                nonsynx = centry['NS'] + centry['NON'] + centry['EXT']
                totalx = nonsynx + centry['SYN']
                if totalx > 0:
                    try: centry['pSYN'] = rje.logBinomial(centry['SYN'],totalx,fsyn,callobj=self)
                    except:
                        if self.dev(): self.warnLog('%s' % centry)
                        self.errorLog('pSYN Error! %s => pSYN = -1.' % str(cdsdb.makeKey(centry))); centry['pSYN'] = -1
                    try:
                        centry['pNS'] = rje.logBinomial(nonsynx, totalx, (1-fsyn),callobj=self)
                    except:
                        if self.dev(): self.warnLog('%s' % centry)
                        self.errorLog('pNS Error! %s => pNS = -1.' % str(cdsdb.makeKey(centry)));
                        centry['pNS'] = -1
                self.debug(centry)
            # Add CDS Rating
            cdsdb.addField('Rating')
            cdsdb.addField('Integrity')
            for centry  in cdsdb.entries(): # ['NON','NS','SYN','SNP','INS','DEL']
                if centry['NON']: centry['Rating'] = 'Truncation'
                elif centry['EXT']: centry['Rating'] = 'Extension'
                elif not centry['SYN'] and not centry['NS']: centry['Rating'] = 'Perfect'
                elif 0 <= centry['pNS'] < 0.05: centry['Rating'] = 'Positive'
                #elif centry['NS'] > centry['SYN']: centry['Rating'] = 'Positive'
                elif not centry['NS']: centry['Rating'] = 'Conserved'
                elif 0 <= centry['pSYN'] < 0.05: centry['Rating'] = 'Constrained'
                else: centry['Rating'] = 'Neutral'
                if not (centry['INS'] + centry['DEL']):
                    if centry['ID'] == 2 or not matchesallowed: centry['Integrity'] = 'Complete'
                    else: centry['Integrity'] = 'Partial'
                elif (centry['INS'] - centry['DEL']) % 3: centry['Integrity'] = 'Indels'
                else: centry['Integrity'] = 'Indelx3'
            cdsdb.saveToFile()
            # Make CDS Summary of counts for each rating (in log)
            cdsdb.makeField('#Rating#|#Integrity#','Summary')
            cdsdb.indexReport('Summary')
            ## ~ [3d] *.features.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Overall summary per feature
            sumdb.dropFields(['product','details'])
            sumdb.saveToFile()
            ## ~ [3e] *.ftypes.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Overall summary per feature type per contig pair?
            sumdb.setStr({'Name':'ftypes'})
            sumdb.compress(sumfields[:3],default='sum',rules={'RefSNP':'mean','AltCN':'mean','RefCN':'mean','MeanCNV':'mean','MaxCNV':'max','AltQV':'mean'})
            sumdb.dropFields(['position','locus_tag','product'])
            sumdb.saveToFile()
            ## ~ [3f] *.summary.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Overall summary per contig pair
            sumdb = self.db('summary')
            sumdb.compress(lockeys+poskeys,default='max',rules={'RefSNP':'mean','AltCN':'mean','RefCN':'mean','MeanCNV':'mean','MaxCNV':'max','AltQV':'mean'})
            #?# What was this doing?! #?#
            #for sentry in sumdb.entries():
            #    for i in range(len(snptypes)):
            #        if sentry[snptypes[i]]:
            #            for stype in snptypes[i+1:]: sentry[stype] = 0
            #            break
            sumdb.compress(lockeys,default='sum',rules={'RefSNP':'mean','AltCN':'mean','RefCN':'mean','AltQV':'mean'})
            sumdb.keepFields(lockeys+snptypes)
            sumdb.addField('SUBTOT',evalue=0)
            sumdb.addField('SNPTOT',evalue=0)
            for sentry in sumdb.entries():
                for stype in snptypes[1:]:  # Not ID
                    sentry['SNPTOT'] += sentry[stype]
                    if stype in ['INS','DEL']: continue
                    sentry['SUBTOT'] += sentry[stype]
            sumdb.saveToFile()
            totdb = self.db().copyTable(sumdb,'total')
            totdb.addField('TEMP',evalue='TOTAL')
            totdb.compress(['TEMP'],default='sum',rules={'RefSNP':'mean','AltCN':'mean','RefCN':'mean','AltQV':'mean'})
            tentry = totdb.entries()[0]
            tentry['Locus'] = tentry['AltLocus'] = tentry.pop('TEMP')
            sumdb.addEntry(tentry)
            tkey = sumdb.makeKey(tentry)
            sumdb.saveToFile(append=True,savekeys=[tkey])
            ### ~ [4] Remove gene SNPs mapped to other features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SNPByFType'):
                for typedb in self.db().splitTable(mapdb,'feature',asdict=False,keepfield=True):
                    typedb.saveToFile()
                    self.db().deleteTable(typedb)
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
            if self.getBool('AltFT'): lockey = 'AltLocus'; posfield = 'AltPos'
            else: lockey = 'Locus'; posfield = 'Pos'
            gb = self.obj['GenBank']
            snpdb = self.db('SNP')
            #self.debug(snpdb.fields())
            skeys = snpdb.index(lockey)[locus][0:]
            ftdb = self.db('Feature')
            if locus not in ftdb.index('locus'): ftlocus = rje.split(locus,'_')[-1]
            else: ftlocus = locus
            mapdb = self.db('snpmap')
            # MapDB Fields: snpdb.fields()+['Strand','GB']+ftdb.fields()+['SNPType','SNPEffect']
            # MapBD Keys: snpdb.keys()+['feature','start','end'])
            ### ~ [2] Map SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fx = 0  # No. failures
            skipx = 0   # No. skipped features
            if ftlocus not in ftdb.index('locus'):
                self.warnLog('Locus %s not found in Feature Table' % ftlocus)
                return False
            for fkey in ftdb.index('locus')[ftlocus]: # These should be in a consistent order: ['locus','start','feature','position']
                fentry = ftdb.data(fkey)
                #self.debug(fentry)
                ## ~ [2a] Check for features to skip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if fentry['feature'] in self.list['FTSkip']: skipx += 1; continue
                #Fixed!# if 'join' in fentry['position']: self.warnLog('%s %s %s: Introns not handled' % (locus,fentry['feature'],fentry['position']))
                ## ~ [2b] Skip SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while skeys:
                    sentry = snpdb.data(skeys[0])
                    if (sentry[posfield] + len(sentry['REF'])) <= fentry['start']: skeys.pop(0)  # No overlap (gone past)
                    else: break     # Caught up with features!
                ## ~ [2c] Process overlapping SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ftsnps = []     # Tuple of (Pos,REF,ALT)
                sx = 0
                for skey in skeys:
                    sentry = snpdb.data(skey)
                    if sentry[posfield] > fentry['end']: break     # Gone beyond feature
                    ftsnps.append((sentry[posfield],sentry['REF'],sentry['ALT']))
                    self.bugPrint('---\nFT: %s' % fentry)
                    self.bugPrint('SNP: %s' % sentry)
                    # NOTE: SNP Data should now be formatted either:
                    # (a) like RATT output with Pos, AltLocus, AltPos, REF and ALT
                    # (b) Like BCF with Pos, REF and ALT
                    try:
                        try: mentry = rje.combineDict({'Strand':'+','GB':fullseq[sentry[posfield]-1]},sentry)
                        except:
                            self.warnLog('Problem adding GB %s %s %s (%s nt sequence)' % (locus,posfield,rje.iStr(sentry[posfield]-1),rje.iLen(fullseq)))
                            mentry = rje.combineDict({'Strand':'+','GB':'?'},sentry)
                        mentry = rje.combineDict(mentry,fentry)
                        cseq = gb.featureSequence(fentry)
                        cpos = gb.featurePos(fentry,sentry[posfield])  # Position 1-L
                        #if 'join' in fentry['position']: self.debug(cpos)
                        if not cpos: raise ValueError('Feature outside position data!')
                        cpos -= 1   # Change to 0<L
                        mentry['length'] = len(cseq)
                        if sentry['REF'] == '.': mentry['GB'] = '-'
                        if cpos >= 0:
                            try: mentry['GB'] = cseq[cpos]
                            except:
                                #x#print sentry, fentry, len(cseq), cpos;
                                raise
                            if 'complement' in fentry['position']: mentry['Strand'] = '-'
                            elif sentry['REF'] == '.': cpos += 1    # Want to compare alt nt to next reference nt
                        ### >>> OLD >>>
                        if not 'dealingwithintrons':
                            cseq = fullseq[fentry['start']-1:fentry['end']]
                            cpos = sentry[posfield] - fentry['start']
                            #!# Add code to return cseq and cpos for join(...) loci #!#
                            #!# Will need to change feature to intron too if required #!#
                            mentry['length'] = len(cseq)
                            #?#rje_sequence.sequenceKs(sequence,self)
                            try: mentry['GB'] = cseq[cpos]
                            except:
                                #x#print sentry, fentry, len(cseq), cpos
                                raise
                            if sentry['REF'] == '.': cpos += 1; mentry['GB'] = '-'
                            if 'complement' in fentry['position']:
                                cseq = rje_sequence.reverseComplement(cseq)
                                cpos = len(cseq) - cpos - 1
                                mentry['Strand'] = '-'
                        ### <<< OLD <<<<
                        # Convert seq and cpos based on join()
                        # Remember that cpos is now 0<L
                        if cpos < 0: mentry['SNPEffect'] = 'intron'
                        elif fentry['feature'] == 'CDS':
                            transl = rje.matchExp('/transl_table="?(\d+)"?',fentry['details'])
                            if transl: transl = transl[0]
                            else: transl = '1'
                            i = 3   # End of codon
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
                                mentry['SNPEffect'] = '%d:%s->%s' % ((cpos+3)/3,rje_sequence.dna2prot(cseq[i-3:],transl=transl),rje.split(rje_sequence.dna2prot(snpseq[i-3:],transl=transl),'*')[0])
                            else:
                                mentry['SNPEffect'] = '%s%d%s' % (rje_sequence.dna2prot(codon,transl=transl),(cpos+3)/3,rje_sequence.dna2prot(snpcodon,transl=transl))
                                # Update SNPType
                                if sentry['REF'] == sentry['ALT']: mentry['SNPType'] = 'ID'
                                elif mentry['SNPEffect'][0] == mentry['SNPEffect'][-1]: mentry['SNPType'] = 'SYN'
                                elif mentry['SNPEffect'][-1] == '*': mentry['SNPType'] = 'NON'
                                elif mentry['SNPEffect'][0] == '*': mentry['SNPType'] = 'EXT'
                                else: mentry['SNPType'] = 'NS'
                        else: mentry['SNPEffect'] = fentry['feature']
                        # Do not map internal identical nucleotides
                        if mentry['SNPType'] == 'ID' and mentry['Pos'] not in [fentry['start'],fentry['end']]: continue
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
            for sentry in snpdb.indexEntries(lockey,locus):
                if not sentry['Mapped']:
                    mentry = rje.combineDict({'Strand':'na'},sentry)
                    mentry = rje.combineDict(mentry,{'GB':fullseq[sentry[posfield]-1],'feature':'intergenic'})
                    mapdb.addEntry(mentry); ix += 1
            self.printLog('#SNPMAP','%s of %s %s SNPs mapped to intergenic regions only.' % (rje.iStr(ix),rje.iLen(snpdb.indexEntries(lockey,locus)),locus))
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
            for bfile in batchfiles: self.processGene(bfile)
            self.log.printLog('#ZEN',rje_zen.Zen().wisdom())
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def processGene(self,gfile):  ### Main controlling run method for a single gene
        '''Main controlling run method for a single gene.'''
        try:### ~ [1] Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.count(gfile,'.') > 1: return
            seqlist = self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % gfile])
            #!# Replace with unaligned sequences and list of exon positions #!#
            gene = seqlist.info['Basefile']
            self.log.printLog('#GENE','Processing files for %s gene' % gene)
            ## ~ [1b] Regenerate alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            exons = self.loadFromFile('%s_exons.txt' % gene,chomplines=True)
            if not exons: return self.log.errorLog('Missing exon data for %s' % gene,printerror=False)
            cds = rje.replace(seqlist.seq[1].info['Sequence'],'-','')
            aln = '-' * len(seqlist.seq[0].info['Sequence'])
            for exon in exons:
                exon = rje.split(exon)
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
            dna = rje.replace(cds,'-','')
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
                myseq = rje.replace(cds[atg:snp['SNP']],'-','')
                if myseq[-1] != snp['CDSSeq']:
                    self.warnLog('{0} != {1}'.format(myseq, snp['CDSSeq']))
                snp['ProteinPos'] = ((len(myseq)-1) / 3) + 1
                frame = snp['CodonPos'] = [3,1,2][len(myseq) % 3]
                codon = snp['Codon'] = dna[(snp['ProteinPos']-1)*3:snp['ProteinPos']*3]
                snp['MajAA'] = rje_sequence.dna2prot(codon[:frame-1]+snp['MajAllele']+codon[frame:])
                snp['MinAA'] = rje_sequence.dna2prot(codon[:frame-1]+snp['MinAllele']+codon[frame:])
                snp['Nonsyn'] = snp['MajAA'] != snp['MinAA']
                #self.deBug(snp)
                #self.deBug(myseq)
                #self.deBug(rje_sequence.dna2prot(myseq))
                #self.deBug(dna[:snp['ProteinPos']*3])
                #self.deBug(rje_sequence.dna2prot(dna[:snp['ProteinPos']*3]))

            ### ~ [3] Output new data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.snpmap.tdt' % gene
            headers = ['SNP','MajAllele','MinAllele','CDS','ConSeq','CDSSeq','ProteinPos','Codon','CodonPos','MajAA','MinAA','Nonsyn']
            rje.delimitedFileOutput(self,outfile,headers,'\t',rje_backup=True)
            for snp in rje.sortKeys(snpdict): rje.delimitedFileOutput(self,outfile,headers,'\t',snpdict[snp])
                
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

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
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
