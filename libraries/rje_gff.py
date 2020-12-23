#!/usr/bin/python

# See below for name and description
# Copyright (C) 2018 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_gff
Description:  GFF File Parser and Manipulator
Version:      0.2.1
Last Edit:    20/11/20
Webserver:    http://www.slimsuite.unsw.edu.au/servers/gff.php
Copyright (C) 2018  Richard J. Edwards - See source code for GNU License Notice

Function:
    The GFF file given by gffin=FILE will be parsed and the components optionally output to tables, a text file of
    comment lines (starting `#`) and fasta format sequences if given. The GFF filename sets the output prefix,
    which can be over-ridden with basefile=FILE.

    The default fields parsed from the GFF are: `locus, source, feature, start, end, score, strand, phase, attributes`.
    Additional fields can be extracted from the attributes field, using `attributes=LIST`. Setting `attributes="*"` or
    `attributes=all` will extract all attributes into additional fields. Note that the `attributes` field itself will be
    kept unless `attfield=F` is used to remove it.

    `integrity=T` will perform checks that the features do not go outside the range of the parsed sequence-region and/or
    fasta sequences.

    `indelwarn=T` and `stopwarn=T` will identify adjacent CDS features that may have sequencing and/or translation
    errors. `indelwarn=T` looks for adjacent CDS features with the same (or hyplist=LIST) "product" (warnfield=X)
    annotation that are within 3 nt of each other (generally overlapping) and might thus represent a fragmented ORF due
    to a frameshift error. `stopwarn=T` identifies similar features that have exactly one codon between them, which
    could represent an atypical genetic code being mis-translated as a stop codon.

    `joinseq=T` will output joined sequences to `*.joined.gff` and, if sequences are parsed, `*.joined.aa.fas` and
    `*.joined.nt.fas`. For protein sequence translations, `stopwarn` sequences are joined with a `*`. `indelwarn`
    sequences are joined with flanking and internal `xx` pairs that delineate the overlapping parts of each
    annotated protein sequence.

    NOTE: Only GFF3 is currently supported.

Commandline:
    ### ~ Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    gffin=FILE      : Input GFF file to parse [None]
    seqin=FILE      : Optional fasta file of reference sequences [None]
    gfftab=T/F      : Whether to output parsed GFF file as a delimited table with headers [True]
    gffloci=T/F     : Whether to parse sequence-region GFF comments to `*.loci.tdt` [True]
    gffcomment=T/F  : Whether to output parsed GFF comments to `*.comments.txt` [False]
    gfffasta=T/F    : Whether to output parsed GFF sequences to `*.fasta` [False]
    attributes=LIST : List of attributes (X=Y;) to pull out into own fields ("*" or "all" for all) [*]
    attfield=T/F    : Whether to keep the full attribute field as parsed from the GFF file [False]
    gffout=FILE     : Save updated GFF format to FILE [None]
    gffseq=T/F      : Whether to include sequences in updated GFF file [False]

    ### ~ GFF Processing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    integrity=T/F   : Perform GFF integrity check based on parsed sequence-region comments and/or fasta [True]
    indelwarn=T/F   : Perform check for possible indels based on overlapping/close common features [True]
    hypindel=INT    : Number of hypothetical proteins that can be involved in a possible indel (0-2) [1]
    stopwarn=T/F    : Perform check for possible codon table stop codon errors based on close common features [True]
    warnfield=X     : Attribute field to use for generating indel or stop codon warnings [product]
    idfield=X       : Attribute field to use for CDS gene ID [ID]
    hyplist=LIST    : List of warnfield values to identify as hypothetical protein ['hypothetical protein']
    cdsfeatures=LIST: List of feature types to count as CDS for warning checks [CDS]
    joinseq=T/F     : Whether to join sequences possible affected by stop codons or frameshifts [False]

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, random, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj
import rje_db, rje_seqlist, rje_sequence
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Basic functional version.
    # 0.1.1 - Modified for splice isoform handling
    # 0.1.2 - Fixed parsing of GFFs with sequence-region information interspersed with features.
    # 0.1.3 - Added option to parseGFF to switch off the attribute parsing.
    # 0.2.0 - Added gff output with ability to fix GFF of tab delimit errors
    # 0.2.1 - Added restricted feature parsing from GFF.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add splitByFeature option that will output GFF split by different features
    # [ ] : Add good and bad feature lists for filtering features.
    # [Y] : More nuanced handling of indels with unknown proteins.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_GFF', '0.2.1', 'November 2020', '2018')
    description = 'GFF File Parser and Manipulator'
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class GFF(rje_obj.RJE_Object):
    '''
    GFF Class. Author: Rich Edwards (2018).

    Str:str
    - GFFIn=FILE      : Input GFF file to parse [None]
    - GFFOut=FILE     : Save updated GFF format to FILE [None]
    - IDField=X       : Attribute field to use for CDS gene ID [ID]
    - WarnField=X     : Attribute field to use for generating indel or stop codon warnings [product]

    Bool:boolean
    - AttField=T/F    : Whether to keep the full attribute field as parsed from the GFF file [True]
    - GFFTab=T/F      : Whether to output parsed GFF file as a delimited table with headers [False]
    - GFFComment=T/F  : Whether to output parsed GFF comments to `*.comments.txt` [False]
    - GFFLoci=T/F     : Whether to parse sequence-region GFF comments to `*.loci.tdt` [False]
    - GFFFasta=T/F    : Whether to output parsed GFF sequences to `*.fasta` [False]
    - GFFSeq=T/F      : Whether to include sequences in updated GFF file [False]
    - IndelWarn=T/F   : Perform check for possible indels based on overlapping/close common features [False]
    - Integrity=T/F   : Perform GFF integrity check based on parsed sequence-region comments and/or fasta [True]
    - JoinSeq=T/F     : Whether to join sequences possible affected by stop codons or frameshifts [False]
    - StopWarn=T/F    : Perform check for possible codon table stop codon errors based on close common features [False]

    Int:integer
    - HypIndel=INT    : Number of hypothetical proteins that can be involved in a possible indel (0-2) [1]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - Attributes=LIST : List of attributes (X=Y;) to pull out into own fields ("*" or "all" for all) []
    - CDSFeatures=LIST: List of feature types to count as CDS for warning checks [CDS]
    - Comments = List of strings storing the comments parsed from the GFF file
    - HypList=LIST   : List of warnfield values to identify as hypothetical protein ['hypothetical protein']

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object, saving parsed features and sequence regions
    - SeqList = rje_seqlist.SeqList object, saving data as tuples
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['GFFIn','GFFOut','IDField','SeqIn','WarnField']
        self.boollist = ['AttField','GFFTab','GFFComment','GFFLoci','GFFFasta','IndelWarn','Integrity','JoinSeq','StopWarn']
        self.intlist = ['HypIndel']
        self.numlist = []
        self.filelist = []
        self.listlist = ['Attributes','CDSFeatures','Comments','HypList']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GFFIn':'None','IDField':'ID','WarnField':'product'})
        self.setBool({'AttField':False,'GFFTab':True,'GFFComment':False,'GFFLoci':True,'GFFFasta':False,'GFFSeq':False,
                      'IndelWarn':True,'Integrity':True,'JoinSeq':False,'StopWarn':True})
        self.setInt({'HypIndel':1})
        self.setNum({})
        self.list['Attributes'] = ['*']
        self.list['HypList'] = ['hypothetical protein']
        self.list['CDSFeatures'] = ['CDS']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
        self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=F'])
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
                self._cmdReadList(cmd,'str',['IDField','WarnField'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['GFFIn','GFFOut'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['AttField','GFFTab','GFFComment','GFFLoci','GFFFasta','GFFSeq',
                                              'IndelWarn','Integrity','JoinSeq','StopWarn'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['HypIndel'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Attributes','CDSFeatures','HypList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        The GFF file given by gffin=FILE will be parsed and the components optionally output to tables, a text file of
        comment lines (starting `#`) and fasta format sequences if given. The GFF filename sets the output prefix,
        which can be over-ridden with basefile=FILE.

        The default fields parsed from the GFF are: `locus, source, feature, start, end, score, strand, phase, attributes`.
        Additional fields can be extracted from the attributes field, using `attributes=LIST`. Setting `attributes="*"` or
        `attributes=all` will extract all attributes into additional fields. Note that the `attributes` field itself will be
        kept unless `attfield=F` is used to remove it.

        `integrity=T` will perform checks that the features do not go outside the range of the parsed sequence-region and/or
        fasta sequences.

        `indelwarn=T` and `stopwarn=T` will identify adjacent CDS features that may have sequencing and/or translation
        errors. `indelwarn=T` looks for adjacent CDS features with the same (or unknown=LIST) "product" (warnfield=X)
        annotation that are within 3 nt of each other (generally overlapping) and might thus represent a fragmented ORF due
        to a frameshift error. `stopwarn=T` identifies similar features that have exactly one codon between them, which
        could represent an atypical genetic code being mis-translated as a stop codon.

        `joinseq=T` will output joined sequences to `*.joined.gff` and, if sequences are parsed, `*.joined.aa.fas` and
        `*.joined.nt.fas`. For protein sequence translations, `stopwarn` sequences are joined with a `*`. `indelwarn`
        sequences are joined with flanking and internal `XX` pairs that delineate the overlapping parts of each
        annotated protein sequence.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.parseGFF(self.getStr('GFFIn'))
            self.checkIntegrity()
            self.cdsWarnings()
            self.saveGFFData()
            self.joinSeq()
            self.restSetup()    #i# Trying to fix REST output issue
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('GFFIn'): raise IOError('No GFF file given! (Set with gffin=FILE)')
            if not rje.exists(self.getStr('GFFIn')): raise IOError('GFF file "%s" not found! (Set with gffin=FILE)' % self.getStr('GFFIn'))
            if rje.exists(self.getStr('SeqIn')):
                 self.obj['SeqList'].loadSeq(self.getStr('SeqIn'))
            elif self.getStrLC('SeqIn'): raise IOError('SeqIn Fasta file "%s" not found! (Set with seqin=FASFILE)' % self.getStr('SeqIn'))
            else: self.printLog('#SEQIN','No seqin=FASFILE sequences provided. Will use sequences from GFF file.')
            if not self.getStrLC('Basefile'):
                self.setBasefile(rje.baseFile(self.getStr('GFFIn')))
                self.printLog('#BASE','Output basefile set: %s' % self.baseFile())
                self.restSetup()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        gffin = GFF input file [txt]
        seqin = optional reference sequences fasta file [fas]
        comments = parsed GFF comments [txt]
        features = parsed GFF features [tdt]
        loci = parsed GFF reference sequence loci [tdt]
        fasta = parsed reference sequences from GFF file [fas]
        joined = table of joined CDS ORF features [tdt]
        joined.aa = protein sequence of joined CDS ORF features [fas]
        joined.nt = nucleotide sequence of joined CDS ORF features [fas]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            for skey in ['GFFIn','SeqIn']: self.dict['Output'][skey.lower()] = self.getStr(skey)
            if self.getBool('GFFTab'): self.dict['Output']['features'] = '%s.features.tdt' % self.baseFile()
            if self.getBool('GFFComment'): self.dict['Output']['comments'] = '%s.comments.txt' % self.baseFile()
            if self.getBool('GFFLoci'): self.dict['Output']['loci'] = '%s.loci.tdt' % self.baseFile()
            if self.getBool('GFFFasta'):
                self.dict['Output']['fasta'] = '%s.fasta' % self.baseFile()
                self.dict['Output']['fasta'] = '%s.fasta not yet implemented' % self.baseFile()
            if self.getBool('JoinSeq'):
                self.dict['Output']['joined'] = '%s.joined.tdt' % self.baseFile()
                self.dict['Output']['joined.aa'] = '%s.joined.aa.fas' % self.baseFile()
                self.dict['Output']['joined.nt'] = '%s.joined.nt.fas' % self.baseFile()
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['gffin','seqin','comments','features','loci','fasta','joined','joined.aa','joined.nt']
#########################################################################################################################
    ### <3> ### GFF Parsing Methods                                                                                     #
#########################################################################################################################
    def parseGFF(self,gfile,parseattributes=True,attfields=None,fix=True,ftypes=()):      ### Generic method
        '''
        Main GFF Parsing method. Parses comments, sequences and features from gfile.
        >> gfile:str = GFF file to parse
        >> parseattributes:bool [True] = Whether to parse attributes
        >> attfields:dict [None] = translation of attributes to field names. Defaults to lower case attribute.
        >> fix:bool [True] = Whether to try to fix entries without tab separation.
        >> ftypes:tuple [()] = Tuple or list of restricted feature types to parse (case sensitive)

        ##gff-version 3
        ##sequence-region fca0000601_BEN4355A1__BEN4355A1F4A.0000601 1 28423
        ##sequence-region fca0001101_BEN4355A1__BEN4355A1F4A.0001101 1 41383
        ...
        pad39_BEN4355A1__BEN4355A1P4P.39        Prodigal:2.6    CDS     120961  122751  .       +       0       ID=BEN4355A1P4_39384;inference=ab initio prediction:Prodigal:2.6,protein motif:Cdd:COG4654;locus_tag=BEN4355A1P4_39384;product=Cytochrome c551/c552
        pad39_BEN4355A1__BEN4355A1P4P.39        Prodigal:2.6    CDS     123180  123419  .       +       0       ID=BEN4355A1P4_39385;inference=ab initio prediction:Prodigal:2.6;locus_tag=BEN4355A1P4_39385;product=hypothetical protein
        pad39_BEN4355A1__BEN4355A1P4P.39        Prodigal:2.6    CDS     123541  123864  .       +       0       ID=BEN4355A1P4_39386;inference=ab initio prediction:Prodigal:2.6,protein motif:CLUSTERS:PHA2517;locus_tag=BEN4355A1P4_39386;product=putative transposase OrfB
        ##FASTA
        >fca0000601_BEN4355A1__BEN4355A1F4A.0000601
        CGTGTTTCATCGGGATTGCTTTCCAGTTTTGATGCAGTTGCATCTATCTAGATGCGGATG
        CATTGATTGTCAATCTTGATGCAGATGCATCTATCCCTATGATGCCGGCATGAATGACAA
        ...
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            sdb = db.addEmptyTable('loci',['locus','start','end','sequence'],['locus'])
            gfields = ['#','locus', 'source', 'ftype', 'start', 'end', 'score', 'strand', 'phase']
            #i# attfield=F does not stop parsing of attributes but stops storage of attributes field
            if self.getBool('AttField'): gfields.append('attributes')
            ## ~ [1a] Setup parsing of attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            addatt = False
            attfieldlist = []
            afields = []
            if not attfields: attfields = {}
            if parseattributes:
                if 'all' in rje.listLower(self.list['Attributes']) or '*' in self.list['Attributes']: addatt = True
                else:
                    #attfieldlist = rje.listLower(self.list['Attributes'])
                    attfieldlist = self.list['Attributes']
                    for fstr in ['WarnField','IDField']:
                        #ffield = self.getStrLC(fstr)
                        #if ffield and ffield not in gfields:
                        #    gfields.append(ffield)
                        #    self.printLog('#FIELD','Added %s=%s field to Attributes parsing list' % (fstr,ffield))
                        if self.getStrLC(fstr) and self.getStr(fstr) not in attfieldlist:
                            attfieldlist.append(self.getStr(fstr))
                            self.printLog('#FIELD','Added %s=%s to Attributes parsing list' % (fstr,self.getStr(fstr)))
                    for field in attfieldlist:
                        if field not in attfields:
                            if field.lower() in gfields:
                                attfields[field] = 'att-{0}'.format(field.lower())
                                self.warnLog('Attribute "%s" will be parsed as "att-%s"' % (field,field.lower()))
                            else: attfields[field] = field.lower()
                        elif attfields[field] in gfields: self.warnLog('Attribute "%s" will replace GFF field "%s"' % (field,field.lower()))
                        if attfields[field] not in gfields:
                            afields.append(attfields[field])
            #gdb = db.addEmptyTable('features',gfields,['locus','strand','start','end','source','ftype'])
            gdb = db.addEmptyTable('features',gfields+afields,['#'])
            gffdata = gdb.dict['Data'] = {}   # Entries being parsed
            gx = 0  # Feature counter
            nx = 0  # Skipped feature counter
            comments = self.list['Comments'] = []
            ### ~ [2] Parse file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# NOTE: This is not designed to be memory efficient. May want a memsaver version in future.
            glines = open(gfile,'r').readlines()
            ## ~ [2a] Comments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('#GFF','Parsing GFF comments...')
            while glines and glines[0].startswith('#'):
                gtext = rje.chomp(glines.pop(0))
                comments.append(gtext)
                seqdata = rje.matchExp('##sequence-region (\S+) (\d+) (\d+)',gtext)
                if seqdata: sdb.addEntry({'locus':seqdata[0],'start':int(seqdata[1]),'end':int(seqdata[2]),'sequence':''})
            self.printLog('\r#COMM','%s comments parsed; %s sequence-regions' % (rje.iLen(comments),rje.iStr(sdb.entryNum())))
            sregx = sdb.entryNum()
            ## ~ [2b] Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ftypes:
                self.debug('%s' % ftypes)
                ftnr = rje.sortUnique(ftypes)
                if '' in ftnr: ftnr.remove('')
                self.printLog('#GFF','Parsing %d feature types: %s' % (len(ftnr),','.join(ftnr)))
            else: self.printLog('#GFF','Parsing all feature types.')
            self.progLog('#GFF','Parsing GFF features...')
            prand = random.randint(1,200)
            while glines:
                if self.dev() or not prand:
                    self.progLog('\r#GFF','Parsing GFF features... %s features; %s sequence-regions; %s features skipped' % (rje.iStr(gx),rje.iStr(sregx),rje.iStr(nx)))
                    prand = random.randint(1,200)
                else: prand -= 1
                if glines[0].startswith('##FASTA'): break
                gtext = rje.chomp(glines.pop(0))
                if gtext.startswith('#'):
                    comments.append(gtext)
                    seqdata = rje.matchExp('##sequence-region (\S+) (\d+) (\d+)',gtext)
                    if seqdata: sdb.addEntry({'locus':seqdata[0],'start':int(seqdata[1]),'end':int(seqdata[2]),'sequence':''})
                    sregx = sdb.entryNum()
                    continue
                gdata = string.split(gtext,'\t')
                gentry = {}
                try:
                    if ftypes and gdata[2] not in ftypes: nx +=1; continue
                    for col in ['locus', 'source', 'ftype', 'start', 'end', 'score', 'strand', 'phase', 'attributes']:
                        gentry[col] = gdata.pop(0)
                    for col in ['start', 'end']:
                        gentry[col] = int(gentry[col])
                    if gdata: raise ValueError('Too many fields in GFF line!')
                except:
                    try:
                        if not fix: raise
                        #i# Add attempt to fix
                        gdata = gtext.split()
                        for col in ['locus', 'source', 'ftype', 'start', 'end', 'score', 'strand', 'phase']:
                            gentry[col] = gdata.pop(0)
                        for col in ['start', 'end']:
                            gentry[col] = int(gentry[col])
                        gentry['attributes'] = ' '.join(data)
                    except: self.errorLog('Problem parsing GFF line: %s' % gtext)

                if parseattributes:
                    if addatt:
                        attlist = rje.longCmd(string.split(gentry['attributes'],';'))
                        #self.debug('%s' % attlist)
                        for attdata in attlist:
                            if not attdata: continue
                            try:
                                [att,val] = string.split(attdata,'=')
                                if att not in attfields:
                                    if att.lower() in gfields:
                                        attfields[att] = 'att-{0}'.format(att.lower())
                                        self.warnLog('Attribute "%s" will be parsed as "att-%s"' % (att,att.lower()))
                                    else: attfields[att] = att.lower()
                                gentry[attfields[att]] = val
                                if attfields[att] not in gdb.fields(): gdb.addField(attfields[att])
                            except:
                                self.warnLog('Problem with GFF attribute: {0}'.format(attdata))
                    else:
                        attstr = ';{0};'.format(gentry['attributes'])
                        for att in attfieldlist:
                            if rje.matchExp(';{0}=([^;]+);'.format(att),attstr):
                                gentry[attfields[att]] = rje.matchExp(';{0}=([^;]+);'.format(att),attstr)[0]
                gx += 1
                gentry['#'] = gx
                gffdata[gentry['#']] = gentry
            #if self.dev(): self.printLog('\r#GFF','%s features and %s sequence-regions parsed from %s' % (rje.iStr(gx),rje.iStr(sregx),gfile))
            self.printLog('\r#GFF','Parsing GFF features complete: %s features; %s sequence-regions; %s features skipped.' % (rje.iStr(gx),rje.iStr(sregx),rje.iStr(nx)))
            self.printLog('\r#GFF','%s features and %s sequence-regions parsed from %s' % (rje.iStr(gdb.entryNum()),rje.iStr(sregx),gfile))
            gdb.indexReport('ftype')
            ## ~ [2c] Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fx = 0
            if glines and glines[0].startswith('##FASTA'):
                self.progLog('#GFF','Parsing GFF fasta...')
                glines.pop(0)
                while glines and glines[0].startswith('>'):
                    locus = string.split(glines.pop(0))[0][1:]
                    sequence = ''
                    while glines and not glines[0].startswith('>'): sequence += rje.chomp(glines.pop(0))
                    if not sequence: self.warnLog('Fasta for %s missing sequence' % locus)
                    fx += 1
                    try: sdb.data(locus)['sequence'] = sequence
                    except:
                        self.warnLog('Parsed fasta for locus not found in sequence-regions: %s' % locus)
                        sdb.addEntry({'locus':locus,'start':1,'end':len(sequence),'sequence':sequence})
            self.printLog('\r#FASTA','Fasta sequences read for %s of %s loci' % (rje.iStr(fx),rje.iStr(sdb.entryNum())))
            self.printLog('#PARSE','Parsing of %s complete.' % gfile)
            if glines: self.warnLog('%s GFF lines remain unparsed' % rje.iLen(glines))
            return True
        except: self.errorLog('%s.parseGFF(%s) error' % (self.prog(),gfile))
#########################################################################################################################
    def checkIntegrity(self):      ### Add/compare the sequences from the fasta file and check positions versus sequence.
        '''
        Add/compare the sequences from the fasta file and check positions versus sequence.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# ['locus','start','end','sequence']
            sdb = self.db('loci')
            #i# ['locus', 'source', 'ftype', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
            gdb = self.db('features')
            seqlist = self.obj['SeqList']
            seqdict = seqlist.makeSeqNameDic()
            ### ~ [2] Check sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            errx = 0
            for entry in sdb.entries():
                locus = entry['locus']
                if seqdict:
                    try: sequence = seqdict[locus]
                    except:
                        self.printLog('#FAIL','Locus "%s" not found in seqin=FILE' % locus)
                        errx += 1; continue
                    if entry['sequence'] and entry['sequence'] != sequence:
                        self.printLog('#FAIL','Locus "%s" sequence mismatch!' % locus)
                        errx += 1; continue
                    else: entry['sequence'] = sequence
                if entry['end'] > len(entry['sequence']):
                    self.printLog('#FAIL','Locus "%s" sequence too short versus sequence-region!' % locus)
                    errx += 1; continue
            ### ~ [3] Check features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in gdb.entries():
                locus = entry['locus']
                if entry['start'] < sdb.data(locus)['start']:
                    self.printLog('#FAIL','%s feature start outside sequence-region: %s' % (locus,entry['attributes']))
                    errx += 1; continue
                if entry['end'] > sdb.data(locus)['end']:
                    self.printLog('#FAIL','%s feature end outside sequence-region: %s' % (locus,entry['attributes']))
                    errx += 1; continue
                #!# Add checking of start and stop #!#

            ### ~ [4] Report error count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if errx: raise ValueError('%s integrity errors. See #FAIL messages in log.' % rje.iStr(errx))
            self.printLog('#CHECK','GFF integrity check complete: no errors spotted.')
            return True
        except: self.errorLog('%s.checkIntegrity error' % self.prog()); return False
#########################################################################################################################
    def cdsWarnings(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('features')
            sdb = self.db('loci')
            if self.getBool('StopWarn') or self.getBool('IndelWarn'): fdb.addField('warn')
            else: return False
            pfield = self.getStr('WarnField')
            idfield = self.getStr('IDField')
            stopx = 0
            longstopx = 0
            indelx = 0
            overx = 0
            seqdict = {}; mx = 0
            for entry in sdb.entries():
                if entry['sequence']:
                    seqdict[entry['locus']] = entry['sequence']
                else: mx += 1
            if self.getBool('JoinSeq') and mx:
                self.errorLog('Set joinseq=T but %s of %s loci missing sequences: joinseq=F.' % (rje.iStr(mx),rje.iStr(sdb.entryNum())),printerror=False)
                self.setBool({'JoinSeq':False})
            joinseq = self.getBool('JoinSeq')
            if joinseq:
                jdb = self.db().addEmptyTable('joined',fdb.fields()+['hypothetical','idlist','aaseq'],fdb.keys())
            ### ~ [2] Scan CDS features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# CDS entries, sorted by: 'locus','start','end','source','feature'
            cdslist = fdb.indexEntries('ftype',self.list['CDSFeatures'])
            cdsx = len(cdslist)
            #!# Need to check the order these come back in for multiple CDS types.
            self.printLog('#CDS','%s CDS features to scan for CDS warnings' % rje.iLen(cdslist))
            jentry = {}
            while len(cdslist) > 1:
                cdsi = cdslist.pop(0)
                cdsj = cdslist[0]
                self.bugPrint('%s\nvs\n%s' % (cdsi,cdsj))
                skip = cdsi['locus'] != cdsj['locus'] or cdsi['strand'] != cdsj['strand']
                #i# Check for matching product
                unki = cdsi[pfield] in self.list['HypList']
                unkj = cdsj[pfield] in self.list['HypList']
                if not skip: self.bugPrint('prod:%s unki:%s unkj:%s' % (cdsi[pfield] == cdsj[pfield],unki,unkj))
                skip = skip or not (cdsi[pfield] == cdsj[pfield] or unki or unkj)
                if not skip: self.bugPrint('End:%d after Start:%d?' % (cdsi['end'],cdsj['start'] - 4))
                unkx = 0
                if unki: unkx += 1
                if unkj: unkx += 1
                skip = skip or unkx > self.getInt('HypIndel')
                #?# Check for matching source?
                if skip and unkx > 1 and self.getBool('IndelWarn') and cdsi['end'] > cdsj['start'] - 4:
                    cdsi['warn'] = 'overlap:%s' % cdsj[idfield]
                    self.printLog('#OVRLAP','Overlapping hypothetical CDS (%s strand): %s (%s) and %s (%s)' % (cdsi['strand'],cdsi[idfield],cdsi[pfield],cdsj[idfield],cdsj[pfield]))
                    overx += 1
                if skip:
                    if joinseq and jentry:    # Process join
                        jdb.addEntry(jentry)
                        jentry = {}
                    continue
                #i# Note that Prokka includes stop codon in feature. This might be common!
                if self.getBool('StopWarn') and not (cdsj['start'] - cdsi['end'] - 1) % 3 and cdsi['end'] < (cdsj['start'] - 4):
                    if cdsi['strand'] == '-':
                        cdsi['warn'] = 'longstop:%d-%d:%s' % (cdsi['end'],cdsj['start']+2,cdsj[idfield])
                        self.printLog('#TRUNC','Possible readthrough to in-frame downstream ORF (%s strand): %s (%s) to %s (%s)' % (cdsi['strand'],cdsj[idfield],cdsj[pfield],cdsi[idfield],cdsi[pfield]))
                    else:
                        cdsi['warn'] = 'longstop:%d-%d:%s' % (cdsi['end']-2,cdsj['start'],cdsj[idfield])
                        self.printLog('#TRUNC','Possible readthrough to in-frame downstream ORF (%s strand): %s (%s) to %s (%s)' % (cdsi['strand'],cdsi[idfield],cdsi[pfield],cdsj[idfield],cdsj[pfield]))
                    longstopx += 1
                if self.getBool('StopWarn') and (cdsi['end'] == cdsj['start'] - 4 or cdsi['end'] == cdsj['start'] - 1):
                    if cdsi['strand'] == '-':
                        cdsi['warn'] = 'stop:%d-%d:%s' % (cdsi['end'],cdsi['end']+2,cdsj[idfield])
                        self.printLog('#STOP','Possible readthrough (%s strand): %s (%s) to %s (%s)' % (cdsi['strand'],cdsj[idfield],cdsj[pfield],cdsi[idfield],cdsi[pfield]))
                    else:
                        cdsi['warn'] = 'stop:%d-%d:%s' % (cdsi['end']-2,cdsi['end'],cdsj[idfield])
                        self.printLog('#STOP','Possible readthrough (%s strand): %s (%s) to %s (%s)' % (cdsi['strand'],cdsi[idfield],cdsi[pfield],cdsj[idfield],cdsj[pfield]))
                    stopx += 1
                    if joinseq:
                        if not jentry:
                            jentry = rje.combineDict({'idlist':cdsi[idfield],'hypothetical':unkx,'shift':0},cdsi)
                            jentry['aaseq'] = seqdict[jentry['locus']][jentry['start']-1:jentry['end']]
                            if jentry['strand'] == '-': jentry['aaseq'] = rje_sequence.reverseComplement(jentry['aaseq'])
                            self.bugPrint(jentry['aaseq'])
                            jentry['aaseq'] = rje_sequence.dna2prot(jentry['aaseq'])
                            if '*' in jentry['aaseq'][:-1]: self.warnLog('%s has internal stop codons!' % cdsi[idfield])
                            self.bugPrint(jentry['aaseq'])
                        else:
                            jentry['warn'] += ',%s' % cdsi['warn']
                            if unkj: jentry['hypothetical'] += 1
                        jentry['end'] = cdsj['end']
                        jentry[idfield] = '%s+' % jentry[idfield]
                        jentry['idlist'] += ',%s' % cdsj[idfield]
                        if jentry[pfield] in self.list['HypList']: jentry[pfield] = cdsj[pfield]
                        aaj = seqdict[cdsj['locus']][cdsj['start']-1:cdsj['end']]
                        if jentry['strand'] == '-': aaj = rje_sequence.reverseComplement(aaj)
                        aaj = rje_sequence.dna2prot(aaj)
                        if '*' in aaj[:-1]: self.warnLog('%s has internal stop codons!' % cdsj[idfield])
                        if jentry['strand'] == '-':
                            if cdsi['end'] == cdsj['start'] - 4:
                                jentry['aaseq'] = jentry['aaseq'] + '*' + aaj
                            else:
                                jentry['aaseq'] = jentry['aaseq'] + aaj
                        else:
                            if cdsi['end'] == cdsj['start'] - 4:
                                jentry['aaseq'] = aaj + '*' + jentry['aaseq']
                            else:
                                jentry['aaseq'] = aaj + jentry['aaseq']
                    self.bugPrint(jentry['aaseq'])
                elif self.getBool('IndelWarn') and cdsi['end'] > cdsj['start'] - 4:
                    x = min(cdsi['end'],cdsj['start'])
                    y = max(cdsi['end'],cdsj['start'])
                    cdsi['warn'] = 'indel:%d-%d:%s' % (x,y,cdsj[idfield])
                    self.printLog('#INDEL','Possible frameshift indel (%s strand): %s (%s) and %s (%s)' % (cdsi['strand'],cdsi[idfield],cdsi[pfield],cdsj[idfield],cdsj[pfield]))
                    indelx += 1
                    if joinseq:
                        if not jentry:
                            jentry = rje.combineDict({'idlist':cdsi[idfield],'hypothetical':unkx,'shift':0},cdsi)
                            jentry['aaseq'] = seqdict[jentry['locus']][jentry['start']-1:jentry['end']]
                            if jentry['strand'] == '-': jentry['aaseq'] = rje_sequence.reverseComplement(jentry['aaseq'])
                            self.bugPrint(jentry['aaseq'])
                            jentry['aaseq'] = rje_sequence.dna2prot(jentry['aaseq'])
                            if '*' in jentry['aaseq'][:-1]: self.warnLog('%s has internal stop codons!' % cdsi[idfield])
                            self.bugPrint(jentry['aaseq'])
                        else:
                            jentry['warn'] += ',%s' % cdsi['warn']
                            if unkj: jentry['hypothetical'] += 1
                        jentry['idlist'] += ',%s' % cdsj[idfield]
                        jentry[idfield] = '%s+' % jentry[idfield]
                        jentry['end'] = cdsj['end']
                        if jentry[pfield] in self.list['HypList']: jentry[pfield] = cdsj[pfield]
                        aaj = seqdict[cdsj['locus']][cdsj['start']-1:cdsj['end']]
                        if jentry['strand'] == '-': aaj = rje_sequence.reverseComplement(aaj)
                        aaj = rje_sequence.dna2prot(aaj)
                        if '*' in aaj[:-1]: self.warnLog('%s has internal stop codons!' % cdsj[idfield])
                        #i# Shift is the change in reading frame from the current cds to the next one
                        shift = (cdsi['start'] % 3) - (cdsj['start'] % 3)
                        #i# Convert shift into an insertion (need to bump cdsj forwards one) or deletion (bring cdsj back)
                        if shift in [1,-2]: indel = 1
                        elif shift in [-1,2]: indel = -1
                        else: indel = 0     # Might be a double indel!
                        #i# jentry['shift'] contains sum of previous indels
                        jlen = jentry['end'] - jentry['start'] + 1 + indel + jentry['shift']
                        jentry['shift'] += indel
                        #i# Need to remove previous joinseq overlaps for calculating new overlap
                        jxx = string.split(jentry['aaseq'],'xx')
                        self.bugPrint(jentry['aaseq'])
                        while len(jxx) > 1:
                            jxx = [jxx[0]+jxx[1]+jxx[3]] + jxx[4:]
                        jxx = jxx[0]
                        self.bugPrint(jxx)
                        if jlen % 3: self.warnLog('%s Frameshift joinseq is not an even number of codons!' % jentry[idfield])
                        if jentry['strand'] == '-':
                            if aaj.endswith('*'): aaj = aaj[:-1]
                            overlap =  len(jxx) + len(aaj) - jlen/3
                            if overlap > 0:
                                jentry['aaseq'] = aaj[:-overlap] + 'xx' + aaj[-overlap:] + 'xx' + jentry['aaseq'][:overlap] + 'xx' + jentry['aaseq'][overlap:]
                            else:
                                if overlap < 0: aaj += '*'; overlap += 1
                                if overlap == -1: aaj += 'X'
                                elif overlap < -1:
                                    self.warnLog('%s Frameshift joinseq has negative overlap (%d)!' % (jentry[idfield],overlap))
                                    self.bugPrint(aaj)
                                    self.debug(jentry)
                                jentry['aaseq'] = aaj + jentry['aaseq']
                        else:
                            if jentry['aaseq'].endswith('*'): jentry['aaseq'] = jentry['aaseq'][:-1]
                            overlap =  len(jxx) + len(aaj) - jlen/3
                            if overlap > 0:
                                jentry['aaseq'] = jentry['aaseq'][:-overlap] + 'xx' + jentry['aaseq'][-overlap:] + 'xx' + aaj[:overlap] + 'xx' + aaj[overlap:]
                            else:
                                if overlap < 0: jentry['aaseq'] += '*'; overlap += 1
                                if overlap == -1: jentry['aaseq'] += 'X'
                                elif overlap < -1:
                                    self.warnLog('%s Frameshift joinseq has negative overlap (%d)!' % (jentry[idfield],overlap))
                                    self.bugPrint(aaj)
                                    self.debug(jentry)
                                jentry['aaseq'] = jentry['aaseq'] + aaj
                        self.bugPrint('%s\n\n' % jentry['aaseq'])
                elif joinseq and jentry:    # Process join
                    jdb.addEntry(jentry)
                    jentry = {}
            if joinseq and jentry: jdb.addEntry(jentry)
            ### ~ [3] Report warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodtxt = 'of %s CDS features without warnings' % rje.iStr(cdsx)
            if self.getBool('StopWarn'):
                self.printLog('#STOP','%s Possible readthroughs' % rje.iStr(stopx))
                self.printLog('#TRUNC','%s Possible long-range readthroughs. (Not joined if joinseq=T)' % rje.iStr(longstopx))
                cdsx = cdsx - stopx
                goodtxt += ' (excluding long-range readthroughs)'
            else:
                self.printLog('#STOP','stopwarn=F')
            if self.getBool('IndelWarn'):
                self.printLog('#INDEL','%s Possible frameshift indels' % rje.iStr(indelx))
                self.printLog('#OVRLAP','%s Overlapping hypothetical CDS' % rje.iStr(overx))
                cdsx = cdsx - indelx - overx
            else:
                self.printLog('#INDEL','indelwarn=F')
            self.printLog('#CDS','%s %s.' % (rje.iStr(cdsx),goodtxt))
            return True     # Method successful
        except: self.errorLog('Problem during %s cdsWarnings.' % self.prog()); return False  # Method failed
#########################################################################################################################
    def saveGFFData(self):    ### Save GFF data in various formats.
        '''
        Save GFF data in various formats:
        - GFFTab=T/F      : Whether to output parsed GFF file as a delimited table with headers [False]
        - GFFComment=T/F  : Whether to output parsed GFF comments to `*.comments.txt` [False]
        - GFFLoci=T/F     : Whether to parse sequence-region GFF comments to `*.loci.tdt` [False]
        - GFFFasta=T/F    : Whether to output parsed GFF sequences to `*.fasta` [False]
        - GFFOut=FILE     : Save updated GFF format to FILE [None]
        - GFFSeq=T/F      : Whether to include sequences in updated GFF file [False]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('GFFTab'):
                self.db('features').saveToFile()
            if self.getBool('GFFLoci'):
                #?# Add summarise data from SeqList?
                self.db('loci').saveToFile(savefields=['locus','start','end'])
            if self.getBool('GFFComment'):
                cfile = '%s.comments.txt' % self.baseFile()
                rje.backup(self,cfile,appendable=False)
                open(cfile,'w').write(string.join(self.list['Comments'],'\n'))
                self.printLog('#OUT','GFF Comments output to %s.' % cfile)
                self.dict['Output']['comments'] = cfile    #i# Trying to fix REST output issue
            if self.getBool('GFFFasta'):
                self.printLog('#DEV','GFFFasta output not yet implemented.')
            if self.getStrLC('GFFOut'):
                gffdb = self.db('features')
                if not gffdb or not gffdb.entryNum():
                    self.printLog('#OUT','Cannot save to GFFOut: no GFF features parsed.')
                else:
                    gfields = ['locus', 'source', 'ftype', 'start', 'end', 'score', 'strand', 'phase']
                    #i# Save comments first then, table with selected headers, minus the headers!
                    gffdb.saveToFile(self.getStr('GFFOut'),delimit='\t',backup=True,append=False,savefields=gfields,log=True,headers=False,comments=self.list['Comments'])

            return True     # Setup successful
        except: self.errorLog('Problem during %s saveGFFData.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def joinSeq(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdb = self.db('joined')
            sdb = self.db('loci')
            if not (self.getBool('JoinSeq') and jdb): return False
            ### ~ [2] Save JoinSeq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jfields = jdb.fields()
            if not self.dev(): jfields.remove('aaseq')
            if jdb.entryNum(): self.db('joined').saveToFile(savefields=jfields)
            else: self.printLog('#JOIN','No sequence joins to output'); return False
            jfas = '%s.joined.aa.fas' % self.baseFile()
            jnt = '%s.joined.nt.fas' % self.baseFile()
            rje.backup(self,jfas)
            rje.backup(self,jnt)
            JFAS = open(jfas,'a')
            JNT = open(jnt,'a')
            for jkey in jdb.dataKeys():
                jentry = jdb.data(jkey)
                jid = jentry[self.getStr('IDField')]
                jdesc = '%s %s%s (%s) hypothetical:%d' % (jentry[self.getStr('WarnField')],jid,jentry['warn'],jentry['idlist'],jentry['hypothetical'])
                i = 1
                while jid[-i] == '+': i += 1
                i -= 1
                jid = '%s.%d' % (jid[:-i],i)
                JFAS.write('>%s %s\n%s\n' % (jid,jdesc,jentry['aaseq']))
                ntseq = sdb.data(jentry['locus'])['sequence'][jentry['start']-1:jentry['end']]
                if jentry['strand'] == '-': ntseq = rje_sequence.reverseComplement(ntseq)
                JNT.write('>%s %s\n%s\n' % (jid,jdesc,ntseq))
            JFAS.close()
            JNT.close()
            self.printLog('#JFAS','%s joinseq output to %s' % (rje.iStr(jdb.entryNum()),jfas))
            self.warnLog('Multiple frameshift joins (3+ ORFs) may have messed up overlaps')
            self.printLog('#DEV','joinSeq GFF and nt fasta output not yet implemented.')
            #i# Check that loci have sequences and make seqdict
            return True     # Setup successful
        except: self.errorLog('Problem during %s joinSeq.' % self.prog()); return False  # Setup failed
#########################################################################################################################
### End of SECTION II: GFF Class                                                                                        #
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
    try: GFF(mainlog,cmd_list).run()

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
