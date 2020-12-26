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
Module:       Snapper
Description:  Genome-wide SNP Mapper
Version:      1.7.0
Last Edit:    30/01/19
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    Snapper is designed to generate a table of SNPs from a BLAST comparison of two genomes, map those SNPs onto genome
    features, predict effects and generate a series of output tables to aid exploration of genomic differences.

    A basic overview of the Snapper workflow is as follows:

    1. Read/parse input sequences and reference features.

    2. All-by-all BLAST of query "Alt" genome against reference using GABLAM.

    3. Reduction of BLAST hits to Unique BLAST hits in which each region of a genome is mapped onto only a single region
    of the other genome. This is not bidirectional at this stage, so multiple unique regions of one genome may map onto
    the same region of the other.

    4. Determine Copy Number Variation (CNV) for each region of the genome based on the unique BLAST hits. This is
    determined at the nucleotide level as the number of times that nucleotide maps to unique regions in the other genome,
    thus establishing the copy number of that nucleotide in the other genome.

    5. Generate SNP Tables based on the unique local BLAST hits. Each mismatch or indel in a local BLAST alignment is
    recorded as a SNP.

    6. Mapping of SNPs onto reference features based on SNP reference locus and position.

    7. SNP Type Classification based on the type of SNP (insertion/deletion/substitution) and the feature in which it
    falls. CDS SNPs are further classified according to codon changes.

    8. SNP Effect Classification for CDS features predicting their effects (in isolation) on the protein product.

    9. SNP Summary Tables for the whole genome vs genome comparison. This includes a table of CDS Ratings based on the
    numbers and types of SNPs. For the `*.summary.tdt` output is, each SNP is only mapped to a single feature according
    to the FTBest hierarchy, removing SNPs mapping to one feature type from feature types lower in the list:
    - CDS,mRNA,tRNA,rRNA,ncRNA,misc_RNA,gene,mobile_element,LTR,rep_origin,telomere,centromere,misc_feature,intergenic

    Version 1.1.0 introduced additional fasta output of the genome regions with zero coverage in the other genome, i.e.
    the regions in the *.cnv.tdt file with CNV=0. Regions smaller than `nocopylen=X` [default=100] are deleted and then
    those within `nocopymerge=X` [default=20] of each other will be merged for output. This can be switched off with
    `nocopyfas=F`.

    Version 1.6.0 added filterself=T/F to filter out self-hits prior to Snapper pipeline. seqin=FILE sequences that are
    found in the Reference (matched by name) will be renamed with the prefix `alt` and output to `*.alt.fasta`. This is
    designed for identifying unique and best-matching homologous contigs from whole genome assemblies, where seqin=FILE
    and reference=FILE are the same. In this case, it is recommended to increase the `localmin=X` cutoff.

    Version 1.7.0 add the option to use minimap2 instead of BLAST+ for speed, using `mapper=minimap`.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input genome to identify variants in []
    reference=FILE  : Fasta (with accession numbers matching Locus IDs) or genbank file of reference genome. []
    basefile=FILE   : Root of output file names (same as SNP input file by default) [<SNPFILE> or <SEQIN>.vs.<REFERENCE>]
    nocopyfas=T/F   : Whether to output CNV=0 fragments to *.nocopy.fas fasta file [True]
    nocopylen=X     : Minimum length for CNV=0 fragments to be output [100]
    nocopymerge=X   : CNV=0 fragments within X nt of each other will be merged prior to output [20]
    makesnp=T/F     : Whether or not to generate Query vs Reference SNP tables [True]
    localsAM=T/F    : Save local (and unique) hits data as SAM files in addition to TDT [False]
    filterself=T/F  : Filter out self-hits prior to Snapper pipeline (e.g for assembly all-by-all) [False]
    mapper=X        : Program to use for mapping files against each other (blast/minimap) [blast]
    ### ~ Reference Feature Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    spcode=X        : Overwrite species read from file (if any!) with X if generating sequence file from genbank [None]
    ftfile=FILE     : Input feature file (locus,feature,position,start,end) [*.Feature.tdt]
    ftskip=LIST     : List of feature types to exclude from analysis [source]
    ftbest=LIST     : List of features to exclude if earlier feature in list overlaps position [(see above)]
    ### ~ SNP Mapping Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    snpmap=FILE     : Input table of SNPs for standalone mapping and output (should have locus and pos info) [None]
    snphead=LIST    : List of SNP file headers []
    snpdrop=LIST    : List of SNP fields to drop []
    altpos=T/F      : Whether SNP file is a single mapping (with AltPos) (False=BCF) [True]
    altft=T/F       : Use AltLocus and AltPos for feature mapping (if altpos=T) [False]
    localsort=X     : Local hit field used to sort local alignments for localunique reduction [Identity]
    localmin=X      : Minimum length of local alignment to output to local stats table [10]
    localidmin=PERC : Minimum local %identity of local alignment to output to local stats table [0.0]
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
import rje, rje_db, rje_obj, rje_genbank, rje_seqlist, rje_sequence, snp_mapper
import rje_blast_V2 as rje_blast
import gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Tidied up with improved run pickup.
    # 0.2.0 - Added FASTQ and improved CNV output along with all features.
    # 0.2.1 - Fixed local output error. (Query/Qry issue - need to fix this and make consistent!) Fixed snp local table revcomp bug.
    # 0.2.2 - Corrected excess CNV table output (accnum AND shortname).
    # 0.2.3 - Corrected "intron" classification for first position of features. Updated FTBest defaults.
    # 1.0.0 - Working version with completed draft manual. Added to SeqSuite.
    # 1.0.1 - Fixed issues when features missing.
    # 1.1.0 - NoCopy fasta output
    # 1.2.0 - makesnp=T/F : Whether or not to generate Query vs Reference SNP tables [True]
    # 1.3.0 - localsAM=T/F : Save local (and unique) hits data as SAM files in addition to TDT [False] - via GABLAM
    # 1.4.0 - localidmin=PERC : Minimum local %identity of local alignment to output to local stats table [0.0]
    # 1.4.1 - Modified warning for AccNum/Locus mismatch in Reference.
    # 1.5.0 - Added pNS and modified the "Positive" CDS rating to be pNS < 0.05.
    # 1.6.0 - filterself=T/F  : Filter out self-hits prior to Snapper pipeline (e.g for assembly all-by-all) [False]
    # 1.6.0 - Added renaming of alt sequences that are found in the Reference for self-comparisons.
    # 1.6.1 - Fixed bug for reducing to unique-unique pairings that was over-filtering.
    # 1.7.0 - Added mapper=minimap setting, compatible with GABLAM v2.30.0 and rje_paf v0.1.0.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add processing of reference genbank file.
    # [Y] : Add GABLAM SNPTable analysis. (Run, generate inverted local table and repeat.)
    # [Y] : Add SNP_Mapper mapping of genbank features to SNP Tables.
    # [Y] : Add reading of QV from FASTQ.
    # [Y] : Add CNV directly to feature table.
    # [Y] : Add feature Start/End positions to SNP table.
    # [ ] : Add splitting of local hits if gap string too long.
    # [ ] : Add compilation of adjacent gaps into longer indels.
    # [ ] : Consolidate CNV of regions using localmin. (Currently determined on a per-nucleotide level)
    # [ ] : Consider making some of the outputs dependent on extras=X setting.
    # [ ] : Change reference to refgenome and use rje_genbank SetupReference() method.
    # [ ] : Add some summary R plots of data.
    # [ ] : Add option to run without some of the SNP-Ft output for compatibility with PAGSAT.
    # [ ] : Add method to generate SNP Frequency plot from table.
    # [ ] : Tidy up the run() method comments and documentation.
    # [ ] : Add CDS GFF output: see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for gaps
    # [ ] : Check that the snp.tdt CN fields are giving the correct output.
    # [Y] : Update to use rje_paf.py once in place.
    # [ ] : Update to use --cs output from minimap2 when mapper=minimap, rather than rje_paf conversion.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Snapper', '1.7.0', 'January 2019', '2016')
    description = 'Genome-wide SNP Mapper'
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
def setupProgram(extra_cmd=[]): ### Basic Setup of Program when called from commandline.
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
        cmd_list = rje.getCmdList(sys.argv[1:]+extra_cmd,info=info)   # Reads arguments and load defaults from program.ini
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
### SECTION II: Snapper Class                                                                                           #
#########################################################################################################################
class Snapper(rje_obj.RJE_Object):
    '''
    Snapper Class. Author: Rich Edwards (2015).

    Str:str
    - SeqIn=FASFILE   : Input genome to identify variants in []
    - Reference=FILE  : Fasta (with accession numbers matching Locus IDs) or genbank file of reference genome. []
    - SNPMap=FILE     : Input table of SNPs for standalone mapping and output (should have locus and pos info) [None]

    Bool:boolean
    - FilterSelf=T/F  : Filter out self-hits prior to Snapper pipeline (e.g for assembly all-by-all) [False]
    - MakeSNP=T/F     : Whether or not to generate Query vs Reference SNP tables [True]
    - NoCopyFas=T/F   : Whether to output CNV=0 fragments to *.nocopy.fas fasta file [True]

    Int:integer
    - NoCopyLen=X     : Minimum length for CNV=0 fragments to be output [100]
    - NoCopyMerge=X   : CNV=0 fragments within X nt of each other will be merged prior to output [20]

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database
    - GenBank = rje_genbank.GenBank for Reference (SNPMap.obj['GenBank'])
    - SNPMap = snp_mapper.SNPMap
    - RefSeqList = rje_seqlist.SeqList for Reference (SNPMap.obj['SeqList'])
    - SeqList = rje_seqlist.SeqList for SeqIn
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['SeqIn','Reference','SNPMap']
        self.boollist = ['FilterSelf','MakeSNP','NoCopyFas']
        self.intlist = ['NoCopyLen','NoCopyMerge']
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'FilterSelf':False,'MakeSNP':True,'NoCopyFas':True})
        self.setInt({'NoCopyLen':100,'NoCopyMerge':20})
        self.setNum({})
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
                self._cmdRead(cmd,type='str',att='Reference',arg='refgenome')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['SeqIn','Reference','SNPMap'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['FilterSelf','MakeSNP','NoCopyFas'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['NoCopyLen','NoCopyMerge'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.debugging(): self.printLog('#CMD','Snapper Commandlist: %s' % string.join(self.cmd_list))
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Setup objects and process Genbank File
            if not self.setup(): raise ValueError('Sequence setup failed.')
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpmap = self.obj['SNPMap']
            ## ~ [2a] ~ Standalone Feature mapping for existing SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Update this to be able to add CNV data if desired (e.g. SNP Table from another source) #!#
            if self.getStrLC('SNPMap'): return self.standaloneSNPMap()
            ## ~ [2b] ~ Generate SNP Tables with GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add capacity to pick up #!#
            # Setup GABLAM commandline to have correct input/output files
            gcmd = ['seqin=%s' % snpmap.getStr('SeqIn'),'searchdb=%s' % self.getStr('SeqIn'),'basefile=%s' % self.baseFile()]
            # Add commands for generating correct data tables
            gcmd += ['local=T','localunique=T','snptable=F','nrseq=F','fullblast=T','localaln=T','dismat=F','distrees=F','reftype=both']
            #if self.dev(): gcmd += ['debug=F']
            # Run GABLAM to generate local and unique tables
            # Recognise and load intermediate tables
            rungablam = False
            ufile = '%s.unique.tdt' % self.baseFile()
            locdb = None; db = None; blast = None
            if rje.exists(ufile) and not self.force():
                self.printLog('#UNIQUE','%s found! (force=F).' % ufile)
                try:
                    blast = rje_blast.BLASTRun(self.log,['blastf=F']+gcmd)
                    db = self.obj['DB'] = blast.obj['DB'] = rje_db.Database(self.log,gcmd)
                    db.addTable(name='hitsum',mainkeys=['Qry'],delimit='\t')
                    db.addTable(name='gablam',mainkeys=['Qry','Hit'],delimit='\t')
                    locdb = db.addTable(name='local',mainkeys=['Qry','Hit','AlnNum'],delimit='\t')
                    if not locdb or not locdb.hasField('SbjSeq'):
                        rungablam = True
                        db.list['Tables'] = []
                        self.printLog('#GABLAM','Need QrySeq and AlnSeq for SNP Table: re-running GABLAM.')
                    else:
                        #!# Update to recongnise Qry, Hit, AlnNum keys
                        db.addTable(name='unique',mainkeys=['Query','Hit','AlnID'],delimit='\t')
                    self.warnLog('Using %s (force=F): check filterself=T/F settings.' % ufile)
                except: self.errorLog('Unique Table present but processing error!'); rungablam = True
            else: rungablam = True
            # Setup GABLAM object - handles some parameters as well as running GABLAM if needed
            gabobj = self.obj['GABLAM'] = gablam.GABLAM(self.log,['qryacc=F','dna=T','blastp=blastn','localmin=10','selfhit=T']+self.cmd_list+gcmd)
            if self.i() >= 0 and gabobj.getInt('LocalMin') < 10:
                gabobj.setInt({'LocalMin':rje.getInt('Recommended min local length >=10 for LocalUnique SNP Table',default='0',confirm=True)})
            #i# FilterSelf currently works by switching selfhit=F. Need to check whether this works sufficiently.
            gabobj.setBool({'SelfHit':not self.getBool('FilterSelf')})

            ## Run GABLAM if required
            if rungablam:
                gabobj.run()
                self.debug('GABLAM has run!')
                blast = gabobj.obj['BLAST']
                db = self.obj['DB'] = blast.obj['DB']    # NB. GABLAM does not have a DB object! Handled by BLAST.
                locdb = blast.db('Local')
                if not locdb: locdb = blast.db('local')
            ftdb = snpmap.db('Feature')
            if ftdb: db.list['Tables'].append(ftdb)
            #self.bugPrint(blast.db().tables())
            #self.debug(blast.db().tableNames())

            ## Make sure data formatting is correct
            dformat = {'AlnID':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int',
            'SbjStart':'int','SbjEnd':'int','Length':'int','Positives':'int','AlnNum':'int','QryLen':'int','HitLen':'int'}
            for table in db.tables(): table.dataFormat(dformat)

            ### ~ [3] ~ Rename Alt sequence found in the Reference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Need to rename awkward sequence names: adding "alt" prefix
            #i# Reference is the Query and the "alt" genome assembly (SeqIn) is the hit
            ## ~ [3a] ~ Identify duplicate names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            refseq = snpmap.obj['SeqList']; refloci = []
            self.bugPrint(refseq.names())
            altseq = self.obj['SeqList']; altloci = []
            self.bugPrint(altseq.names())
            newname = {}
            for seq in refseq.seqs():
                sname = refseq.shortName(seq)
                refloci.append(sname)
            for seq in altseq.seqs():
                sname = altseq.shortName(seq)
                if sname in refloci:
                    newname[sname] = 'alt%s' % sname
                    altloci.append('alt%s' % sname)
                    self.warnLog('AltSeq %s in Reference. Renamed alt%s to avoid issues.' % (sname,sname))
                else: altloci.append(sname)
            self.printLog('#SNAME','%s of %s Alt sequence names also found in Reference' % (rje.iLen(newname),rje.iLen(altloci)))
            ## ~ [3b] ~ Replace duplicate names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for table in db.list['Tables']:
                if not newname or 'Hit' not in table.fields(): continue
                hx = 0
                for entry in table.entries():
                    if entry['Hit'] in newname: entry['Hit'] = newname[entry['Hit']]; hx += 1
                self.printLog('#RENAME','Hit renamed in %s of %s %s entries.' % (rje.iStr(hx),rje.iStr(table.entryNum()),table.name()))
                table.remakeKeys()
            ## ~ [3c] ~ Update SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if newname:
                newfile = '%s.alt.fasta' % self.basefile()
                rje.backup(self,newfile)
                NEWFAS = open(newfile,'w')
                for seq in altseq.seqs():
                    sname = altseq.shortName(seq)
                    if sname in newname: NEWFAS.write('>alt%s\n%s\n' % altseq.getSeq(seq))
                    else: NEWFAS.write('>%s\n%s\n' % altseq.getSeq(seq))
                NEWFAS.close()
                self.printLog('#SAVE','%s sequences saved to %s' % (rje.iStr(altseq.seqNum()),newfile))
                altseq.loadSeq(newfile)

            ### ~ [4] ~ Main Snapper Pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            udb = blast.db('unique')        #!# Need to rationalise all of these into the same DB object!
            if not udb: raise ValueError('BLAST.reduceLocal() failed: No SNP Table output.')
            # Create reference unique tables by inverting Qry and Sbj
            locdb.renameField('Qry','Query')
            locdb.renameField('AlnNum','AlnID')
            refdb = self.db().copyTable(locdb,'ref.local')
            refdb.dataFormat({'AlnID':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'})
            # Invert the Query and Hit data
            for entry in refdb.entries():
                [entry['Query'],entry['Hit'],entry['QryStart'],entry['QryEnd'],entry['QrySeq'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq']] = [entry['Hit'],entry['Query'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq'],entry['QryStart'],entry['QryEnd'],entry['QrySeq']]
                if entry['QryStart'] > entry['QryEnd']:
                    [entry['QryStart'],entry['QryEnd'],entry['SbjStart'],entry['SbjEnd']] = [entry['QryEnd'],entry['QryStart'],entry['SbjEnd'],entry['SbjStart']]
                    entry['QrySeq'] = rje_sequence.reverseComplement(entry['QrySeq'])
                    entry['SbjSeq'] = rje_sequence.reverseComplement(entry['SbjSeq'])
            refudb = blast.reduceLocal(refdb,sortfield=gabobj.getStr('LocalSort'),minloclen=gabobj.getInt('LocalMin'))
            blast.db().deleteTable(refdb)
            if not refudb: raise ValueError('Reference BLAST.reduceLocal() failed.')
            # Re-invert the Query and Hit data
            for entry in refudb.entries():
                [entry['Query'],entry['Hit'],entry['QryStart'],entry['QryEnd'],entry['QrySeq'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq']] = [entry['Hit'],entry['Query'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq'],entry['QryStart'],entry['QryEnd'],entry['QrySeq']]
                if entry['QryStart'] > entry['QryEnd']:
                    [entry['QryStart'],entry['QryEnd'],entry['SbjStart'],entry['SbjEnd']] = [entry['QryEnd'],entry['QryStart'],entry['SbjEnd'],entry['SbjStart']]
                    entry['QrySeq'] = rje_sequence.reverseComplement(entry['QrySeq'])
                    entry['SbjSeq'] = rje_sequence.reverseComplement(entry['SbjSeq'])
            refudb.setStr({'Name':'refunique'})
            refudb.saveToFile(savefields=['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd'])
            # Generate and combine nocoverage tables
            nocovdb = blast.noCoverage(udb,save=False)
            nocovdb.dropEntriesDirect('QryHit',['Qry'])
            for entry in nocovdb.entries(): entry['QryHit'] = 'Alt'     # Or 'SeqIn'?
            refnocovdb = blast.noCoverage(refudb,save=False)
            refnocovdb.dropEntriesDirect('QryHit',['Qry'])
            for entry in refnocovdb.entries(): entry['QryHit'] = 'Ref'; nocovdb.addEntry(entry)
            nocovdb.renameField('QryHit','Source')
            #!# Is this used? #!#
            # Calculate duplicate regions for SeqIn and Reference

            #!# This needs to be updated in response to altseq renaming

            #refseq = snpmap.obj['SeqList']; refloci = []
            #self.bugPrint(refseq.names())
            #altseq = self.obj['SeqList']; altloci = []
            #self.bugPrint(altseq.names())
            self.dict['DupRegions'] = dupregions = {}     # {locus:[depth list]}
            for seq in refseq.seqs():
                sname = refseq.shortName(seq)
                #refloci.append(sname)
                dupregions[sname] = [0] * refseq.seqLen(seq)
            #!# Catch if SeqIn = Reference (avoid doing this with same names! Check!!)
            for seq in altseq.seqs():
                sname = altseq.shortName(seq)
                #altloci.append(sname)
                if sname in dupregions: self.warnLog('AltSeq %s in Reference. Avoid duplicated names: will cause issues.' % sname)
                dupregions[sname] = [0] * altseq.seqLen(seq)
            # Add AccNum for Feature mapping -> Point to the same list objects
            duploci = rje.sortKeys(dupregions)
            for sname in duploci:
                sacc = string.split(sname,'_')[-1]
                if sacc not in dupregions: dupregions[sacc] = dupregions[sname]
            #self.debug(rje.sortKeys(dupregions))
            #self.debug(udb.indexKeys('Query'))
            #self.debug(refudb.indexKeys('Hit'))
            for entry in udb.entries():
                locus = entry['Query']
                for i in range(entry['QryStart']-1,entry['QryEnd']): dupregions[locus][i] += 1
            for entry in refudb.entries():
                locus = entry['Hit']
                try:
                    if entry['SbjStart'] > entry['SbjEnd']:
                        for i in range(entry['SbjEnd']-1,entry['SbjStart']): dupregions[locus][i] += 1
                    else:
                        for i in range(entry['SbjStart']-1,entry['SbjEnd']): dupregions[locus][i] += 1
                except:
                    self.bugPrint('%s dupregions: %d' % (locus,len(dupregions[locus])))
                    self.bugPrint(entry)
                    self.errorLog('DupRegions error!',quitchoice=True)
            # Generate table of duplicated regions from dupregions (includes nocoverage - could drop this later?)
            # NOTE: These are regions that are duplicated IN THE OTHER GENOME
            cnvdb = db.addEmptyTable('cnv',['Source','Locus','Start','End','Copy'],['Source','Locus','Start'])
            for locus in duploci:
                x = 0; y = 1
                while x < len(dupregions[locus]):
                    while y < len(dupregions[locus]) and dupregions[locus][y] == dupregions[locus][x]: y += 1
                    cnvdb.addEntry({'Source':{True:'Ref',False:'Alt'}[locus in refloci],'Locus':locus,'Start':x+1,'End':y,
                                    'Copy':dupregions[locus][x]})
                    x = y; y += 1
            cnvdb.saveToFile()
            self.featureCNVTable()
            #!# Add a more sophisticated analysis that maps these regions BACK onto the other genome: should be in
            #!# the local table.
            # Generate *.nocopy.fas output
            if self.getBool('NoCopyFas'):
                nocopyfas = '%s.nocopy.fas' % db.baseFile()
                self.progLog('#NOCOPY','Generating %s output...' % (nocopyfas))
                NOCOPY = open(nocopyfas,'w'); nx = 0
                #i# refseq = snpmap.obj['SeqList']; refloci = []
                #i # altseq = self.obj['SeqList']; altloci = []
                for seq in refseq.seqs():
                    (sname,sequence) = refseq.getSeq(seq)
                    sname = refseq.shortName(seq)
                    reglist = []
                    for centry in cnvdb.indexEntries('Locus',sname):
                        if centry['Copy'] > 0: continue
                        if centry['Source'] != 'Ref': self.warnLog('Ref/Alt name conflict! - %s' % sname)
                        regx = int(centry['Start'])
                        regy = int(centry['End'])
                        if (regy - regx + 1) < self.getInt('NoCopyLen'): continue
                        reglist.append((regx,regy))
                    reglist = rje.collapseTupleList(reglist,joindistance=self.getInt('NoCopyMerge'))
                    for (regx,regy) in reglist:
                        NOCOPY.write('>%s.%s-%s\n%s\n' % (sname,rje.preZero(regx,len(sequence)),rje.preZero(regy,len(sequence)),sequence[regx-1:regy]))
                        nx += 1
                for seq in altseq.seqs():
                    (sname,sequence) = altseq.getSeq(seq)
                    sname = altseq.shortName(seq)
                    reglist = []
                    for centry in cnvdb.indexEntries('Locus',sname):
                        if centry['Copy'] > 0: continue
                        if centry['Source'] != 'Alt': self.warnLog('Ref/Alt name conflict! - %s' % sname)
                        regx = int(centry['Start'])
                        regy = int(centry['End'])
                        if (regy - regx + 1) < self.getInt('NoCopyLen'): continue
                        reglist.append((regx,regy))
                    reglist = rje.collapseTupleList(reglist,joindistance=self.getInt('NoCopyMerge'))
                    for (regx,regy) in reglist:
                        NOCOPY.write('>%s.%s-%s\n%s\n' % (sname,rje.preZero(regx,len(sequence)),rje.preZero(regy,len(sequence)),sequence[regx-1:regy]))
                        nx += 1
                self.printLog('\r#NOCOPY','%s sequences output to %s.' % (rje.iStr(nx),nocopyfas))

            if not self.getBool('MakeSNP'):
                self.printLog('#NOSNP','Snapper SNP generation switched off (makesnp=F).')
                return

            # Generate SNP Tables for SeqIn and Reference
            self.printLog('#~~#','## ~~~~~ Query Genome SNP Table ~~~~~ ##')
            snpdb = self.snpTableFromLocal(udb,save=False,addft=True)         # Keys: ['Locus','Pos','AltLocus','AltPos']
            snpdb.setStr({'Name':'alt.snp'})
            self.printLog('#~~#','## ~~~~~ Reference Genome SNP Table ~~~~~ ##')
            refsnpdb = self.snpTableFromLocal(refudb,save=False,addft=True)
            refsnpdb.setStr({'Name':'snp'})
            if self.dev(): refsnpdb.saveToFile('%s.snp.full.tdt' % self.baseFile())
            #refsnpdb.newKey(['AltLocus','AltPos','Locus','Pos'])    # Should now match snpdb
            # - No longer required as refudb re-inverted
            snpdb.addFields(['RefSNP','RefCN','AltCN'],evalue=0)
            refsnpdb.addFields(['RefSNP','RefCN','AltCN'],evalue=0)
            self.progLog('#REFSNP','Generating Reference SNP list...')
            refsnps = rje.listIntersect(snpdb.dataKeys(),refsnpdb.dataKeys())
            self.debug('%s SNP vs %s Ref -> %s RefSNP' % (rje.iLen(snpdb.dataKeys()),rje.iLen(refsnpdb.dataKeys()),rje.iLen(refsnps)))
            ex = 0.0; etot = len(refsnps)
            for ekey in refsnps:
                self.progLog('\r#REFSNP','Generating Reference SNP list: %.1f%%' % (ex/etot)); ex += 100.0
                snpdb.data(ekey)['RefSNP'] = 1
                refsnpdb.data(ekey)['RefSNP'] = 1
            self.printLog('#REFSNP','Generated Reference SNP list: %s SNPs' % rje.iLen(refsnps))
            ex = 0.0; etot = snpdb.entryNum()
            for entry in snpdb.entries():
                self.progLog('\r#CNV','Adding copy number data to SNP Table: %.2f%%' % (ex/etot)); ex += 100.0
                # AltCN uses the number of times the Reference region appears = other Alt loci hitting the same Ref
                entry['AltCN'] = dupregions[entry['Locus']][entry['Pos']-1]
                # RefCN uses the number of times the SeqIn region appears = other Alt loci hitting the same Ref
                entry['RefCN'] = dupregions[entry['AltLocus']][entry['AltPos']-1]
            if etot: self.printLog('\r#CNV','Adding copy number data to SNP Table: %.2f%%' % (ex/etot)); ex += 100.0
            else: self.printLog('#CNV','No SNP table entries for copy number data!')
            snpdb.saveToFile()
            #?# Make optional:
            refsnpdb.dropEntriesDirect('RefSNP',[0])
            ex = 0.0; etot = refsnpdb.entryNum()
            for entry in refsnpdb.entries():
                self.progLog('\r#CNV','Adding copy number data to Reference SNP Table: %.2f%%' % (ex/etot)); ex += 100.0
                # AltCN uses the number of times the Reference region appears = other Alt loci hitting the same Ref
                entry['AltCN'] = dupregions[entry['Locus']][entry['Pos']-1]
                # RefCN uses the number of times the SeqIn region appears = other Alt loci hitting the same Ref
                entry['RefCN'] = dupregions[entry['AltLocus']][entry['AltPos']-1]
            if etot: self.printLog('\r#CNV','Adding copy number data to Reference SNP Table: %.2f%%' % (ex/etot)); ex += 100.0
            else: self.printLog('#CNV','No SNP table entries for copy number data!')
            refsnpdb.saveToFile()
            self.debug('SNP Tables saved')
            #?# Convert reference names to locus (if required) for SNP_Mapper? (Updated SNP_Mapper to handle?)
            #!# Optionally compress long indels. (Check recognition by snp_mapper)
            # Output snp tables and ref.snp tables
            # - OR: flag refsnps in snp table and just output = better [What does this mean? Have a RefSNP Y/N output? Or a number? 1/0?]

            ## ~ [2c] Cross-Reference SNP Tables to Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Make optional?:
            snpmap.setStr({'SNPFile':'%s.snp.tdt' % self.baseFile()})
            #x#snpmap.setBool({'AltFT':True})
            if snpmap.setupSNPTable(keepmatches=True):
                snpmap.db().baseFile(self.baseFile())
                snpmap.featureSNPs()
                if self.dev(): self.printLog('#DEV','Tables: %s.' % string.join(snpmap.db().tableNames(),'; '))
                for table in snpmap.db().tables():
                    if table.name() in ['SNP','cds','ftypes','features','summary','snpmap']:
                        #self.printLog('#NAME','SNPMap %s table renamed "ref.%s"' % (table.name(),table.name()))
                        table.setStr({'Name':'ref.%s' % table.name()})
            # "Regular" snpmap file
            #!# Add extras=X setting to control alt output
            snpmap.setStr({'SNPFile':'%s.alt.snp.tdt' % self.baseFile()})
            if snpmap.setupSNPTable(keepmatches=True):
                snpmap.db().baseFile('%s.alt' % self.baseFile())
                snpmap.featureSNPs()
                #for table in snpmap.db().tables():
                #    if table.name() in ['SNP','cds','ftypes','features','summary','snpmap']:
                #        self.printLog('#NAME','SNPMap %s table renamed "alt.%s"' % (table.name(),table.name()))
                #        table.setStr({'Name':'alt.%s' % table.name()})
            #snpmap.db().getTable('SNP').setStr({'Name':'snp'})



            #?# Cross-reference snp and refsnp summaries to identify duplications? (different SNP numbers)
            #?# Optional output of altered protein sequences
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Reference'): raise IOError('No reference=FILE given!')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['TupleKeys=True'])
            self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
            ## ~ [1a] Setup FastQ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fastqfile = '%s.fastq' % rje.baseFile(self.getStr('SeqIn'))
            self.dict['QScore'] = {}
            if rje.exists(fastqfile):
                #!# Improve this!
                qlines = open(fastqfile,'r').readlines()
                while len(qlines) > 3:
                    name = string.split(qlines.pop(0))[0][1:]
                    sequence = rje.chomp(qlines.pop(0))
                    qlines.pop(0)
                    qscore = rje.chomp(qlines.pop(0))
                    if len(sequence) != len(qscore): self.warnLog('%s Fastq sequence and Qscores of different lengths!' % name)
                    else: self.dict['QScore'][name] = qscore
                self.printLog('#FASTQ','Loaded QV scores for %s sequences from %s.' % (rje.iLen(self.dict['QScore']),fastqfile))
            else: self.printLog('#FASTQ','%s not found: no QV scores' % fastqfile)
            ## ~ [1b] Setup Reference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # This will take a fasta or genbank file. If genbank, will check/generate feature table and use *.full.fas
            # for sequence mapping for SNP effects. If fasta file but <SEQBASE>.gb or <SEQIN>.gb is present, this genbank
            # file will be used if required. Otherwise ftfile=FILE must be present.
            if self.getStrLC('SNPMap'): self.cmd_list.append('snpfile=%s' % self.getStrLC('SNPFile'))
            self.obj['SNPMap'] = snp_mapper.SNPMap(self.log,self.cmd_list)
            self.obj['SNPMap'].setStr({'SeqIn':self.getStr('Reference')})
            self.obj['SNPMap'].setupReference()
            ## ~ [1c] Setup basefile for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.baseFile(return_none=None):
                if self.obj['SNPMap'].getStrLC('SNPFile') and rje.exists(self.obj['SNPMap'].getStr('SNPFile')):
                    newbase = rje.baseFile(self.obj['SNPMap'].getStr('SNPFile'))
                else:
                    newbase = '%s.vs.%s' % (rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(self.getStr('Reference'),strip_path=True))
                self.baseFile(newbase)
            self.printLog('#BASE','Output filename base: %s' % self.baseFile())
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
    def standaloneSNPMap(self): ### Standalone SNP File Feature Mapping
        '''
        Standalone SNP File Feature Mapping.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ Standalone SNP File Feature Mapping ~~~~~ ##')
            snpmap = self.obj['SNPMap']
            mapfile = '%s.snpmap.tdt' % snpmap.baseFile()
            if self.dev(): self.debug(mapfile)
            if not self.force() and rje.exists(mapfile):
                self.printLog('#MAP','SNPMap file found (force=F): %s' % mapfile)
                return True
            snpfile = snpmap.getStr('SNPFile')
            if not rje.exists(snpfile): raise IOError('SNPFile "%s" not found!' % snpfile)
            if not snpmap.setupSNPTable(): raise ValueError('SNP Table setup failed.')
            snpmap.db().baseFile(snpmap.baseFile())
            ### ~ [2] ~ Standalone Feature mapping for existing SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpmap.featureSNPs()
            return True
        except: self.errorLog('%s.standaloneSNPMap() error' % self.prog()); return False
#########################################################################################################################
    def snpTableFromLocal(self,locdb=None,queries=[],hits=[],save=True,addft=True):    ### Outputs a SNP table (similar format to nucmer) from local alignments
        '''
        Outputs a SNP table (similar format to nucmer) from local alignments.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param save:bool = Whether to save SNP table once generated.
        @param addft:bool = Whether to add Feature Start/End positions to SNP Table
        @return: SNP table
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NB. This is based on BLAST SNP Table. Might want to add to this table rather than create separately to
            # enable use of other SNP tables sources!
            db = self.db()
            snpdb = db.addEmptyTable('SNP',['Locus','Pos','REF','ALT','AltLocus','AltPos'],keys=['Locus','Pos','AltLocus','AltPos'])
            if self.dict['QScore']: snpdb.addField('AltQV')
            snpdb.setBool({'TupleKeys':True})
            if not locdb: locdb = self.db('Local')
            if rje.listDifference(queries,locdb.index('Query')):
                self.warnLog('Queries for snpTableFromLocal() not found in local table','missing_queries',suppress=True)
            if rje.listDifference(hits,locdb.index('Hit')):
                self.warnLog('Hits for snpTableFromLocal() not found in local table','missing_hits',suppress=True)
            revcomp = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
            ## ~ [0a] Feature Start/End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add ProgLog messages.
            fex = 0
            ftdb = self.db('Feature')   # Add start and stops to SNP Table
            ftends = {}
            if addft and ftdb:
                for locus in ftdb.index('locus'):
                    ftends[locus] = ftdb.indexDataList('locus',locus,'start') + ftdb.indexDataList('locus',locus,'end')
                self.printLog('#ADDFT','Start and End feature positions compiled for %s loci' % rje.iLen(ftends))
            ### ~ [1] Generate SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lx = 0; ltot = locdb.entryNum()
            for lentry in locdb.entries():
                if queries and lentry['Query'] not in queries: continue
                if hits and lentry['Hit'] not in hits: continue
                self.progLog('\r#SNP','Generating SNP Table: %.1f%%' % (lx/ltot)); lx += 100.0
                fwd = True
                sdir = 1
                if lentry['SbjStart'] > lentry['SbjEnd']: sdir = -1; fwd = False
                qpos = lentry['QryStart'] - 1
                spos = lentry['SbjStart'] - sdir
                for ai in range(len(lentry['AlnSeq'])):
                    if lentry['QrySeq'][ai] != '-': qpos += 1
                    if lentry['SbjSeq'][ai] != '-': spos += sdir
                    ftpos = False
                    # Check for FT end
                    if lentry['Query'] in ftends and qpos in ftends[lentry['Query']]: ftpos = True
                    elif string.split(lentry['Query'],'_')[-1] in ftends and qpos in ftends[string.split(lentry['Query'],'_')[-1]]: ftpos = True
                    elif lentry['Hit'] in ftends and spos in ftends[lentry['Hit']]: ftpos = True
                    elif string.split(lentry['Hit'],'_')[-1] in ftends and spos in ftends[string.split(lentry['Hit'],'_')[-1]]: ftpos = True
                    # Dump unwanted identities
                    if lentry['QrySeq'][ai].upper() == lentry['SbjSeq'][ai].upper() and not ftpos: continue   # Identity
                    sentry = {'Locus':lentry['Query'],'Pos':qpos,'REF':lentry['QrySeq'][ai].upper(),
                              'ALT':lentry['SbjSeq'][ai].upper(),'AltLocus':lentry['Hit'],'AltPos':spos}
                    #???# if not fwd: sentry['ALT'] = revcomp[sentry['ALT']]
                    if lentry['Hit'] in self.dict['QScore']:
                        try:
                            if lentry['SbjSeq'][ai] == '-': # Mean of flanking scores
                                sentry['AltQV'] = phredScore(self.dict['QScore'][lentry['Hit']][spos-1])
                                sentry['AltQV'] += phredScore(self.dict['QScore'][lentry['Hit']][spos])
                                sentry['AltQV'] /= 2.0
                            else: sentry['AltQV'] = phredScore(self.dict['QScore'][lentry['Hit']][spos-1])
                        except: self.warnLog('PhredScore problem for %s:%s' % (lentry['Hit'],spos),'phredscore',suppress=True)
                    if lentry['QrySeq'][ai].upper() == lentry['SbjSeq'][ai].upper() and ftpos: fex += 1 #self.bugPrint(sentry)
                    snpdb.addEntry(sentry)
            if addft: self.printLog('#ADDFT','%s Start and End feature positions added to SNP table' % rje.iStr(fex))
            self.printLog('#SNPDB','%s SNP table entries.    ' % rje.iStr(snpdb.entryNum()))
            ### ~ [2] Output SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: snpdb.saveToFile()
            return snpdb
        except: self.errorLog('Snapper.snpTableFromLocal() error'); return None
#########################################################################################################################
    def featureCNVTable(self):  ### Output feature table with CNV data
        '''Output feature table with CNV data.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdb = self.db('Feature')               # Add start and stops to SNP Table
            if not ftdb: return False
            dupregions = self.dict['DupRegions']    # {locus:[depth list]}
            ftdb.addFields(['MeanCNV','MaxCNV'])    # Mean and Max copy number
            for locus in ftdb.index('locus'):
                if locus not in dupregions:
                    self.warnLog('Locus "%s" not found in DupRegions list. AccNum/Locus mismatch in Reference?' % locus)
                    continue
                for ftentry in ftdb.indexEntries('locus',locus):
                    cnv = dupregions[locus][ftentry['start']-1:ftentry['end']]
                    if not cnv:
                        self.warnLog('Zero length DupRegion for Locus "%s".' % locus)
                        continue
                    ftentry['MeanCNV'] = float(sum(cnv)) / len(cnv)
                    ftentry['MaxCNV'] = max(cnv)
            ftdb.saveToFile('%s.ftcnv.tdt' % self.baseFile())
        except: self.errorLog('Snapper.featureCNVTable() error'); return None
#########################################################################################################################
### End of SECTION II: Snapper Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def phredScore(qchar,qscale=33): ### Returns the Phred quality score for a given character.
    '''Returns the Phred quality score for a given character.'''
    return ord(qchar) - qscale
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
    try: Snapper(mainlog,['tuplekeys=T']+cmd_list).run()

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
