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
Module:       rje_paf
Description:  Minimap2 PAF parser and converter
Version:      0.13.2
Last Edit:    07/08/24
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to use Minimap alignments in PAF format (with the --cs flag) to replace blastn for all-by-all
    genome assembly comparisons for PAGSAT etc.

    Input is a PAF file with --cs flag, produced by minimap:

    `minimap2 --cs -N 100 -p 0.0001 -x asm20 $REFERENCE $ASSEMBLY > $PAFILE`

    Setting `pafin=minimap` will run minimap2 and generate $BASEFILE.paf output, which will then be processed.

    Output is a series of tables that emulate GABLAM output:

    * hitsum.tdt = summary of hits per query sequence.
    * gablam.tdt = summary stats for each query-hit pair.
    * local.tdt = local hits table. Will include alignments if `localaln=T`
    * qryunique.tdt = local hits table reduced to unique coverage of query ($ASSEMBLY) sequences (qryunique=T).
    * hitunique.tdt = local hits table reduced to unique coverage of hit ($REFERENCE) sequences (hitunique=T).

    NOTE: When Snapper runs rje_paf it inverts the reference and assembly so that the query is the reference.

    The default functionality is to emulate GABLAM (and thereby replace blast+ in GABLAM with the mapper=minimap flag).
    However, excess fields can be reduced with `mockblast=T`. Additional fields are output either way.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    pafin=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping []
    basefile=STR    : Base for file outputs [PAFIN basefile]
    seqin=FASFILE   : Input genome to identify variants in ($ASSEMBLY) []
    reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
    assembly=FASFILE: As seqin=FASFILE
    uniqueout=T/F   : Whether to output *.qryunique.tdt and *.hitunique.tdt tables of unique coverage [True]
    uniquehit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [False]
    alnseq=T/F      : Whether to use alnseq-based processing (True) or CS-Gstring processing [False]
    localaln=T/F    : Whether to output local alignments in Local Table [False]
    mockblast=T/F   : Whether to output mock BLAST headers even when not appropriate [True]
    ### ~ Minimap2 run/mapping options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [0]
    endextend=X     : Extend minimap2 hits to end of sequence if query region with X bp of end [0]
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    mapopt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm20]
    mapsplice=T/F   : Switch default minimap2 options to `-x splice -uf -C5` [False]
    reads=FILELIST  : List of fasta/fastq read files for minimap2 mapping to BAM. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    ### ~ Coverage checking mode (dev only) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    checkpos=TDTFILE: File of Locus, Start, End positions for read coverage checking [None]
    checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
    spanid=X        : Generate sets of read IDs that span checkpos regions, grouped by values of field X []
    ### ~ Variant mode (dev only) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    snptableout=T/F : Generated output of filtered variants to SNP Table [False]
    qcut=X          : Min. quality score for a mapped read to be included [0]
    minqn=X         : Min. number of non-N reads meeting qcut for output (after indel filtering) [10]
    rid=T/F         : Whether to include Read ID (number) lists for each allele [True]
    readnames=T/F   : Output the read names to the RID file [False]
    indels=T/F      : Whether to include indels in "SNP" parsing [True]
    skiploci=LIST   : List of loci to exclude from pileup parsing (e.g. mitochondria) []
    mincut=X        : Minimum read count for minor allele (proportion if <1) [0.05]
    absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    biallelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    ignoren=T/F     : Whether to exclude "N" calls from alleles [True]
    ignoreref=T/F   : If False will always keep Reference allele and call fixed change as SNP [True]
    ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
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
import rje, rje_forker, rje_obj, rje_db, rje_seqlist, rje_sequence
#import numpy as np
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version. Compatible with GABLAM v2.30.0 and Snapper v1.7.0.
    # 0.2.0 - Added endextend=X : Extend minimap2 hits to end of sequence if with X bp [10]
    # 0.3.0 - Added mapsplice mode for dealing with transcript mapping.
    # 0.3.1 - Correct PAF splicing bug.
    # 0.4.0 - Added TmpDir and forking for GABLAM conversion.
    # 0.5.0 - Added uniquehit=T/F : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [False]
    # 0.6.0 - Added CS alignment manipulation methods.
    # 0.6.1 - Added additional error-handling for CS parsing errors.
    # 0.7.0 - Added alnseq=T/F : Whether to use alnseq-based processing (True) or CS-Gstring processing (dev only) [False]
    # 0.7.1 - Disabled endextend due to bug.
    # 0.7.2 - Fixed alnlen bug - minimap2 ignore Ns for alignment length calculation. Fixed endextend bug.
    # 0.7.3 - Added minlocid and minloclen filtering to PAF parsing to speed up processing. Set default alnseq=F.
    # 0.8.0 - Added developmental variant calling options.
    # 0.9.0 - Added long-read mapping to BAM option.
    # 0.10.0 - Added longreadMinimapPAF() and checkpos=TDTFILE options.
    # 0.10.1 - Added spanid=X: Generate sets of read IDs that span checkpos regions, based on values of field X []
    # 0.10.2 - Fixed formatting for Python 2.6 back compatibility for servers.
    # 0.10.3 - Fixing issues of PAF files not being generated.
    # 0.11.0 - Added HiFi read type.
    # 0.12.0 - Added readnames=T/F : Output the read names to the RID file [False]
    # 0.13.0 - Updated default -x to asm20 (5% divergence). Use minimap2 map-hifi setting for hifi read mapping.
    # 0.13.1 - Py3 updates.
    # 0.13.2 - Added fix to pafin reading issue.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [?] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [Y] : Read in PAF format
    # [Y] : Convert to local.tdt, unique.tdt and refunique.tdt
    # [Y] : Convert to gablam.tdt (from unique.tdt)
    # [Y] : Convert to hitsum.tdt
    # [ ] : Generate CNV table.
    # [ ] : Generate SNP tables.
    # [N] : Add snapper mode to reverse reference and assembly? (And output SNP tables.)
    # [Y] : Make the GABLAM reduction quiet.
    # [Y] : Add localaln=T : Output alignment sequences to local hits table [False]
    # [Y] : Modify refunique output name - in Snapper, it is reciprocal unique hits (I think).
    # [ ] : Table('Hit',['Query','Rank','Hit','Description','BitScore','E-Value','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Query','Hit'])
    # [ ] : Update to use existing results if found and not force=T.
    # [ ] : Add GFF output from local hits table.
    # [ ] : Add fasta output of hit sequences. (Account for introns.)
    # [Y] : Possible update of cs string when using endextend
    # [ ] : Add minlocid etc. cutoffs to speed up processing.
    # [ ] : Fix endextend bug
    # [Y] : Add reads=LIST etc. and move the long read BAM generation here?
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_PAF', '0.13.1', 'May 2022', '2019')
    description = 'Minimap2 PAF parser and converter'
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
#pafhead = rje.split('Qry QryLen QryStart QryEnd Strand Hit HitLen HitStart HitEnd Identity Length Quality')
pafhead = 'Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality'.split()
# Minimap2 outputs mapping positions in the Pairwise mApping Format (PAF) by default. PAF is a TAB-delimited text format
# with each line consisting of at least 12 fields as are described in the following table:

# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start coordinate (0-based)
# 4	int	Query end coordinate (0-based)
# 5	char	"+" if query/target on the same strand; "-" if opposite
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start coordinate on the original strand
# 9	int	Target end coordinate on the original strand
# 10	int	Number of matching bases in the mapping
# 11	int	Number bases, including gaps, in the mapping
# 12	int	Mapping quality (0-255 with 255 for missing)

pafaln = 'tp cm s1 s2 NM MD AS ms nn ts cg cs dv'.split()
# When alignment is available, column 11 gives the total number of sequence matches, mismatches and gaps in the
# alignment; column 10 divided by column 11 gives the BLAST-like alignment identity. When alignment is unavailable,
# these two columns are approximate. PAF may optionally have additional fields in the SAM-like typed key-value format.
# Minimap2 may output the following tags:

# Tag	Type	Description
# tp	A	Type of aln: P/primary, S/secondary and I,i/inversion
# cm	i	Number of minimizers on the chain
# s1	i	Chaining score
# s2	i	Chaining score of the best secondary chain
# NM	i	Total number of mismatches and gaps in the alignment
# MD	Z	To generate the ref sequence in the alignment
# AS	i	DP alignment score
# ms	i	DP score of the max scoring segment in the alignment
# nn	i	Number of ambiguous bases in the alignment
# ts	A	Transcript strand (splice mode only)
# cg	Z	CIGAR string (only in PAF)
# cs	Z	Difference string
# dv	f	Approximate per-base sequence divergence

# The cs tag encodes difference sequences in the short form or the entire query AND reference sequences in the long form.
# It consists of a series of operations:
#
# Op	Regex	Description
# =	[ACGTN]+	Identical sequence (long form)
# :	[0-9]+	Identical sequence length
# *	[acgtn][acgtn]	Substitution: ref to query
# +	[acgtn]+	Insertion to the reference
# -	[acgtn]+	Deletion from the reference
# ~	[acgtn]{2}[0-9]+[acgtn]{2}	Intron length and splice signal

# This module also adds ! [0-9]+ for extended matches of unknown type
# NOTE: The cs tag (as parsed) is based on the target sequence, such that a -ve strand hit will be the reverse complement
#       of the true hit and always has target start < end.

## NOTE: The following settings will be interesting to play with for different applications:
# -D	If query sequence name/length are identical to the target name/length, ignore diagonal anchors. This option also reduces DP-based extension along the diagonal.
# -P	Retain all chains and don't attempt to set primary chains. Options -p and -N have no effect when this option is in use.
# -p FLOAT	Minimal secondary-to-primary score ratio to output secondary mappings [0.8]. Between two chains overlaping over half of the shorter chain (controlled by -M), the chain with a lower score is secondary to the chain with a higher score. If the ratio of the scores is below FLOAT, the secondary chain will not be outputted or extended with DP alignment later. This option has no effect when -X is applied.
# -N INT	Output at most INT secondary alignments [5]. This option has no effect when -X is applied.
#########################################################################################################################
locformats = {'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int','Positives':'int',
              'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PAF Class                                                                                               #
#########################################################################################################################
class PAF(rje_obj.RJE_Object):
    '''
    PAF Class. Author: Rich Edwards (2019).

    Str:str
    - CheckPos=TDTFILE    : File of Locus, Start, End positions for read coverage checking [None]
    - Minimap2=PROG   : Full path to run minimap2 [minimap2]
    - PAFIn=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping []
    - Reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
    - SeqIn=FASFILE   : Input genome to identify variants in ($ASSEMBLY) []
    - SpanID=X        : Generate sets of read IDs that span checkpos regions, based on values of field X []
    - TmpDir=PATH     : Temporary directory to use when forking [./tmp/]

    Bool:boolean
    - AlnSeq = Whether to use alnseq-based processing (True) or CS-Gstring processing (dev only) [False]
    - BiAllelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    - IgnoreN=T/F     : Whether to exclude "N" calls from alleles [True]
    - IgnoreRef=T/F   : If False will always keep Reference allele and call fixed change as SNP [True]
    - Indels=T/F      : Whether to include indels in "SNP" parsing [True]
    - LocalAln = Whether to keep local alignments in Local Table [False]
    - MapSplice=T/F   : Switch default minimap2 options to `-x splice -uf -C5` [False]
    - MockBLAST = Whether to output mock BLAST headers even when not appropriate [True]
    - ReadNames=T/F   : Output the read names to the RID file [False]
    - RID=T/F         : Whether to include Read ID (number) lists for each allele [True]
    - SNPTableOut=T/F : Generated output of filtered variants to SNP Table [False]
    - UniqueHit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [False]
    - UniqueOut=T/F   : Whether to output *.qryunique.tdt and *.hitunique.tdt tables of unique coverage [True]

    Int:integer
    - AbsMinCut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    - EndExtend=X         : Extend minimap2 hits to end of sequence if with X bp [0]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [0]
    - MinQN=X         : Min. number of non-N reads meeting qcut for output (after indel filtering) [10]
    - QCut=X          : Min. quality score for a mapped read to be included [0]

    Num:float
    - MinCut=X        : Minimum read count for minor allele (proportion if <1) [0.05]
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]

    File:file handles with matching str filenames
    
    List:list
    - CheckFields=LIST    : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    - CheckFlanks=LIST       : List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
    - Reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - SkipLoci=LIST   : List of loci to exclude from pileup parsing (e.g. mitochondria) []

    Dict:dictionary    
    MapOpt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm20]

    Obj:RJE_Objects
    - DB: Database object for main IO
    - Reference: SeqList Object of $REFERENCE
    - Query: SeqList Object of $ASSEMBLY
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['CheckPos','Minimap2','PAFIn','SeqIn','SpanID','Reference','TmpDir']
        self.boollist = ['AlnSeq','BiAllelic','IgnoreN','IgnoreRef','Indels','LocalAln','MapSplice','MockBLAST','ReadNames','RID','SNPTableOut','UniqueHit','UniqueOut']
        self.intlist = ['AbsMinCut','EndExtend','MinLocLen','MinQN','QCut']
        self.numlist = ['MinCut','MinLocID']
        self.filelist = []
        self.listlist = ['CheckFields','CheckFlanks','Reads','ReadType','SkipLoci']
        self.dictlist = ['MapOpt']
        self.objlist = ['Assembly','Reference']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'Minimap2':'minimap2','TmpDir':'./tmp/'})
        self.setBool({'AlnSeq':False,'LocalAln':False,'MapSplice':False,'MockBLAST':True,'UniqueHit':False,'UniqueOut':True,
                      'BiAllelic':False,'IgnoreN':True,'IgnoreRef':True,'Indels':True,'ReadNames':False,'RID':True,'SNPTableOut':False})
        self.setInt({'AbsMinCut':2,'EndExtend':0,'MinLocLen':1,'MinQN':10,'QCut':0})
        self.setNum({'MinCut':0.05,'MinLocID':0.0})
        self.dict['MapOpt'] = {} #'N':'100','p':'0.0001','x':'asm20'}
        self.list['CheckFlanks'] = [0,100,1000,5000]
        self.list['CheckFields'] = ['Locus','Start','End']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['ReadType'] = ['ont']
        self.obj['Forker']  = rje_forker.Forker(self.log,['logfork=F']+self.cmd_list)
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
                self._cmdRead(cmd,type='file',att='SeqIn',arg='assembly')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='file',att='Reference',arg='searchdb')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='file',att='Reference',arg='refgenome')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['Minimap2','SpanID'])   # Normal strings
                self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['CheckPos','PAFIn','SeqIn','Reference'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['AlnSeq','BiAllelic','IgnoreN','IgnoreRef','Indels','LocalAln','MapSplice','MockBLAST','ReadNames','RID','SNPTableOut','UniqueHit','UniqueOut'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['AbsMinCut','EndExtend','MinLocLen','MinQN','QCut'])   # Integers
                self._cmdReadList(cmd,'perc',['MinLocID'])   # 0-100 percentage, converted x100 if <=1
                self._cmdReadList(cmd,'float',['MinCut']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CheckFields','ReadType','SkipLoci'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['CheckFlanks']) # Comma separated list as integers
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Reads']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['MapOpt']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.list['CheckFlanks'] = rje.intList(self.list['CheckFlanks'])
        self.list['CheckFlanks'].sort()
        #if not self.getBool('AlnSeq'): self.warnLog('alnseq=F mode is developmental - please report odd behaviour.')
        #if self.getInt('EndExtend') > 0:
        #    self.warnLog('Endextend>0 bug may cause some incorrect trimming. Watch for alignment positions warnings and consider running with endextend=0')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('CheckPos'): return self.checkPos()
            self.setup()
            if self.getStrLC('PAFIn') in ['minimap','minimap2']:
                self.setStr({'PAFIn':'%s.paf' % self.baseFile()})
                self.printLog('#PAFIN','Minimap2 PAF file set: %s' % self.getStr('PAFIn'))
                self.minimap2()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.parsePAF(): return False
            if self.dev() and self.getBool('SNPTableOut'): return self.pafToVariants()
            if not self.pafToLocal(minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen')): raise ValueError
            self.localToUnique(save=self.getBool('UniqueOut'))   #!# Sort out output of unique tables!
            self.uniqueToHitSum()
            self.localToGABLAM()
            #? Optional deletion of empty self.getStr('TmpDir') ?#
            return self.db()
        except SystemExit: raise    # Fork child
        except:
            self.errorLog('RJE_PAF run error')
            return None  # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Basefile'): self.baseFile(rje.baseFile(self.getStr('PAFIn')))
            if not self.getStrLC('Basefile'): self.baseFile('rje_paf')
            self.printLog('#BASE',self.baseFile())
            if self.getBool('MapSplice'): defaults = {'x':'splice'} #,'uf':'','C5':''}
            else: defaults = {'N':'100','p':'0.0001','x':'asm20'}
            self.dict['MapOpt'] = rje.combineDict(defaults,self.dict['MapOpt'],overwrite=True)
            #self.devLog('#MAPOPT','%s' % self.dict['MapOpt'])
            if self.getInt('Forks') > 0:
                if 't' not in self.dict['MapOpt']: self.dict['MapOpt']['t'] = self.getInt('Forks')
            elif self.i() >=0 and rje.yesNo('Forking recommended to speed up operation? Set forks>0?'):
                self.setInt({'Forks':max(0,rje.getInt('Set number of forks:',confirm=True))})
                self.printLog('#FORKS','Forks set to forks=%d.' % self.getInt('Forks'))
            else: self.printLog('#FORKS','Forking recommended to speed up operation: set forks>0 to use multiple threads.')
            ## ~ [1a] Database Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T','basefile=%s' % self.baseFile()])
            ## ~ [1b] Reference sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('AlnSeq') and not self.getStrLC('Reference'):
                raise ValueError('Reference missing: %s' % rje.argString(self.cmd_list))
            if self.getStrLC('Reference'):
                self.printLog('#REFIN','Loading reference (target) sequences.')
                self.obj['Reference'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Reference'),'autoload=T','seqmode=file'])
            ## ~ [1c] Reference sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('AlnSeq') and not self.getStrLC('SeqIn'):
                raise ValueError('SeqIn missing: %s' % rje.argString(self.cmd_list))
            if self.getStrLC('SeqIn'):
                self.printLog('#SEQIN','Loading assembly (query) sequences.')
                self.obj['Assembly'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('SeqIn'),'autoload=T','seqmode=file'])
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
        There is currently no specific help available on REST output for this program.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### PAF Parsing Methods                                                                                     #
#########################################################################################################################
    def parsePAF(self,pafin=None,table='paf',filter=True,debugstr='NOTDEBUGGING'):  ### Parse PAF into Database table.
        '''
        Parse PAF into Database table.
        >> pafin:file [None] = PAF file to load, or use
        >> table:str ['paf'] = Name of table to enter data into
        >> filter:bool [True] = Whether to impose length/identity filtering during parsing.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not pafin:
                if not self.getStrLC('PAFIn'): raise ValueError('No PAFIn given!')
                pafin = self.getStr('PAFIn')
            if not rje.exists(pafin): raise IOError('PAFIn file "%s" not found!' % pafin)
            self.setStr({'PAFIn':pafin})
            tmpshush = pafin.endswith('.tmp') and not self.debugging()
            #self.printLog('#~~#','## ~~~~~ Parsing PAF Alignments ~~~~~ ##')
            if not tmpshush:
                self.headLog('Parsing PAF Alignments')
                self.printLog('\r#PAF','Parsing %s' % pafin)

            ### ~ [2] Load PAF file with auto-counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lfiltx = 0; ifiltx = 0
            minloclen = self.getInt('MinLocLen')
            minlocid = self.getPerc('MinLocID')
            headx = len(pafhead)
            pafdb = db.addEmptyTable(table,['#']+pafhead+pafaln,['#'])
            px = 0
            with open(pafin,'r') as PAF:
                for line in PAF:
                    if not tmpshush: self.progLog('\r#PAF','Parsing %s lines' % (rje.iStr(px)))
                    px += 1
                    data = rje.chomp(line).split('\t')
                    #i# Added for checkpos
                    if len(data) < len(pafhead): continue
                    # Generate entry
                    pentry = {'#':px}
                    for i in range(headx):
                        field =  pafhead[i]
                        if field in ['Qry','Hit','Strand']: pentry[field] = data[i]
                        else: pentry[field] = rje.atoi(data[i])
                    for pdat in data[headx:]:
                        pentry[pdat[:2]] = pdat[3:]
                    if pentry['Qry'].endswith(debugstr): self.bugPrint('%s' % pentry)
                    if filter:
                        nn = 0
                        try:
                            if 'nn' in pentry: nn = rje.atoi(pentry['nn'][2:])
                        except:
                            if self.debugging(): self.errorLog('Unexpected error!')
                        if (pentry['Length'] + nn) < minloclen:
                            lfiltx += 1
                            if pentry['Qry'].endswith(debugstr): self.debug('MinLen: %s' % pentry)
                            continue
                        if float(pentry['Identity']+nn)/(pentry['Length']+nn) < minlocid:
                            if pentry['Qry'].endswith(debugstr): self.debug('MinID: %s' % pentry)
                            ifiltx += 1; continue
                    # Reformat data to 1->L positions
                    for field in 'QryStart SbjStart'.split(): pentry[field] += 1
                    # Tidy the 'cs' field
                    if 'cs' in pentry and pentry['cs'].startswith('Z:'): pentry['cs'] = pentry['cs'][2:]
                    pafdb.addEntry(pentry)
                    # Debugging temp code
                    if self.debugging() == 'DISABLE':
                        cstats = self.statsFromCS(pentry['cs'])
                        nn = rje.atoi(pentry['nn'][2:])
                        self.bugPrint(cstats)
                        self.bugPrint(pafdb.entrySummary(pentry,collapse=True))
                        self.deBug('Length=%s; Identity=%s; nn=%s; :%s; Alnlen=%s' % (pentry['Length'],pentry['Identity'],nn,cstats[':'],cstats['AlnLen']))


            if not tmpshush:
                self.printLog('\r#PAF','Parsed %s lines from %s' % (rje.iStr(px),pafin))
                #if not pafdb.entryNum(): self.warnLog('No hits parsed from %s: check minimap2=PROG setting' % pafin)
                if lfiltx: self.printLog('#FILT','%s PAF lines filtered with (Length+nn) < minloclen=%d' % (rje.iStr(lfiltx),minloclen))
                if ifiltx: self.printLog('#FILT','%s PAF lines filtered with (Identity+nn)/(Length+nn) < minlocid=%.2f%%' % (rje.iStr(ifiltx),minlocid*100.0))

            ### ~ [3] Return PAF Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if self.dev(): pafdb.saveToFile()
            return pafdb
        except: self.errorLog('%s.parsePAF error' % self.prog()); return False
#########################################################################################################################
    def pafToLocal(self,pafdb=None,save=True,alnseq=False,minlocid=0,minloclen=1):  ### Parse PAF table into Local alignment table.
        '''
        Parse PAF table into Local alignment table.
        >> pafin:file [None] = PAF file to load, or use
        >> save:bool [True] = Save local hits table
        >> alnseq:bool [False] = Whether to use pafToLocalAlnSeq method
        >> minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        >> minlocid:int [1] = Minimum local length to keep.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alnseq or self.getBool('AlnSeq') or self.getBool('LocalAln'): return self.pafToLocalAlnSeq(pafdb,save,True,minlocid,minloclen)
            self.headLog('PAF CS to Local Alignment Conversion')
            if not pafdb: pafdb = self.db('paf')
            locdb = self.db().copyTable(pafdb,'local')
            locfields = 'Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd'.split()
            locfields += ['QryLen','SbjLen']

            ### ~ [2] Add AlnNum and re-index data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locdb.joinFields('Pair',['Qry','Hit'],join='|',replace=False)
            locdb.rankFieldByIndex(index='Pair',rankfield='Identity',newfield='AlnNum',rev=True,absolute=True,unique=True)
            locdb.newKey(['Qry','Hit','AlnNum'])
            locdb.makeField(formula='#Identity#',fieldname='Positives',evalue=-1)
            locdb.makeField(formula='#Identity#',fieldname='BitScore',evalue=-1)
            locdb.makeField(fieldname='Expect',evalue=-1)
            locdb.dataFormat({'Positives':'int','BitScore':'int','Expect':'num','AlnNum':'int'})
            locdb.keepFields(locfields+['Strand','nn','cs'])
            locdb.list['Fields'] = locfields+['Strand','nn','cs']

            ### ~ [3] Process CS strings into Local dictionary table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lx = 0.0; ltot = locdb.entryNum(); px = 0
            badlenx = 0; badidx = 0     # Length and Identity mismatches
            extx = 0;   # Count of end extensions
            for lentry in locdb.entries():
                self.progLog('#ALN','Converting PAF cs alignments: %.1f%%' % (lx/ltot)); lx += 100.0
                if not lentry['cs']: continue
                cstats = self.statsFromCS(lentry['cs'])
                #if cstats['n'] > 0 and cstats['AlnLen'] == (lentry['Length'] + cstats['n']):
                #    self.debug('AlnLen should not include the %d n' % cstats['n'])
                nn = rje.atoi(lentry['nn'][2:])
                if (cstats['AlnLen'] + nn) != lentry['Length']:
                    badlenx += 1
                    self.deBug('%s\n-> %s -> %s' % (locdb.entrySummary(lentry,collapse=True),lentry['cs'],cstats))
                if cstats[':'] != lentry['Identity']:
                    badidx += 1
                    #self.deBug('%s\n-> %s -> %s' % (lentry,lentry['cs'],cstats))
                ## ~ [3a] EndExtend ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                aln = lentry['cs']
                hang5 = lentry['QryStart'] - 1
                sbj5 = lentry['SbjStart'] - 1
                hang3 = lentry['QryLen'] - lentry['QryEnd']
                sbj3 = lentry['SbjLen'] - lentry['SbjEnd']
                rev = lentry['Strand'] == '-'
                if rev: (sbj5,sbj3) = (sbj3,sbj5)
                if 0 < hang5 <= self.getInt('EndExtend'):
                    extx += 1
                    if hang5 > sbj5:
                        endgap = hang5 - sbj5
                        #if rev: aln = '%s?%d!%d' % (aln,sbj5,endgap)
                        #else:
                        aln = '!%d?%d%s' % (endgap,sbj5,aln)
                    else:
                        #if rev: aln = '%s?%d' % (aln,hang5)
                        #else:
                        aln = '?%d%s' % (hang5,aln)
                    lentry['QryStart'] -= hang5
                    if rev: lentry['SbjEnd'] += min(sbj5,hang5)
                    else: lentry['SbjStart'] -= min(sbj5,hang5)
                if 0 < hang3 <= self.getInt('EndExtend'):
                    extx += 1
                    if hang3 > sbj3:
                        endgap = hang3 - sbj3
                        #if rev: aln = '!%d?%d%s' % (endgap,sbj3,aln)
                        #else:
                        aln = '%s?%d!%d' % (aln,sbj3,endgap)
                    else:
                        #if rev: aln = '?%d%s' % (hang3,aln)
                        #else:
                        aln = '%s?%d' % (aln,hang3)
                    lentry['QryEnd'] += hang3
                    if lentry['Strand'] == '+': lentry['SbjEnd'] += min(sbj3,hang3)
                    else: lentry['SbjStart'] -= min(sbj3,hang3)
                lentry['cs'] = aln
                ## ~ [3b] Generate stats from CS string ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cstats = self.statsFromCS(lentry['cs'])
                lentry['Length'] = cstats['AlnLen']
                lentry['Identity'] = cstats[':']
                px += 1
            self.printLog('#ALN','Conversion of %s PAF cs alignments complete.' % rje.iStr(px))
            if extx:
                self.warnLog('%s end extensions (endextend=%d): may be some minor identity underestimates.' % (rje.iStr(extx),self.getInt('EndExtend')))
            if badlenx:
                self.warnLog('%s Length/AlnLen mismatches from PAF file. (May be disparity with N treatment.)' % rje.iStr(badlenx))
            if badidx:
                self.warnLog('%s Identity mismatches from PAF file (PAF identity ignores Ns).' % rje.iStr(badidx))
            if px != ltot:
                raise ValueError('%s cs alignments missing from PAF file: check --cs flag was used' % rje.iStr(ltot-px))

            ### ~ [4] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if minloclen > 1: locdb.dropEntries('Length<%d' % minloclen)
            if minlocid > 0:
                locdb.makeField('100.0*Identity/Length','LocID')
                locdb.dropEntries('LocID<%s' % minlocid)
                locdb.dropField('LocID')
            #self.deBug('?')

            ### ~ [5] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: locdb.saveToFile(savefields=locfields+['Strand','cs'])

            return locdb
        except: self.errorLog('%s.pafToLocal error' % self.prog()); return False
#########################################################################################################################
    def localToUnique(self,locdb=None,save=True):  ### Convert local alignments to unique coverage of queries and hits.
        '''Convert local alignments to unique coverage of queries and hits.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AlnSeq'): return self.localToUniqueAlnSeq(locdb,save)
            self.headLog('PAF Local Alignment to Unique Hits')
            if not locdb: locdb = self.db('local')
            locfields = 'Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd'.split()

            ### ~ [2] Add AlnNum and re-index data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniqdb = self.reduceLocalCS(locdb,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
            uniqdb.rename('hitunique')
            qryuniqdb = self.reduceLocalCS(locdb,byqry=True,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
            qryuniqdb.rename('qryunique')

            ### ~ [3] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save:
                uniqdb.saveToFile(savefields=locfields)
                qryuniqdb.saveToFile(savefields=locfields)

            return uniqdb
        except: self.errorLog('%s.localToUnique error' % self.prog()); return False
#########################################################################################################################
    def uniqueToHitSum(self,save=True):  ### Parse PAF Local alignment table into HitSum table.
        '''Parse PAF Local alignment table into HitSum table.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('PAF Local and Unique Hits to HitSum')
            self.printLog('#INFO','Converting query-unique reduced local hits into hit summary statistics')
            locdb = self.db('local')
            if self.getBool('UniqueHit'):
                self.printLog('#INFO','Using hit-unique reduced local hits for GABLAM summary statistics')
                locdb = self.db('hitunique')
            locdb.index('Qry')
            uniqdb = self.db('qryunique')
            sumdb = self.db().copyTable(uniqdb,'hitsum')
            sumfields = 'Qry HitNum MaxScore EVal Description Length Coverage Identity Positives'.split()
            rules={'Coverage':'sum','Identity':'sum','Positives':'sum','QryEnd':'max','QryStart':'min','CovHit':'sum'}
            sumfields += ['CovHit','QryStart','QryEnd','Qry_AlnLen','Qry_AlnID']

            ### ~ [2] Add AlnNum and re-index data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sumdb.renameField('BitScore','MaxScore')
            sumdb.renameField('Expect','EVal')
            sumdb.renameField('QryLen','Length')
            sumdb.makeField('QryEnd-QryStart+1','Coverage')
            sumdb.dataFormat({'Coverage':'int'})
            sumdb.addField('CovHit',evalue=1)
            sumdb.keepFields(sumfields+['Hit','AlnNum'])
            sumdb.compress(newkeys=['Qry'],rules=rules,default='max',best=[],joinchar='|')
            sumdb.makeField('100.0*Coverage/Length','Qry_AlnLen')
            sumdb.makeField('100.0*Identity/Length','Qry_AlnID')
            for field in sumfields:
                if field not in sumdb.fields(): sumdb.addField(field)
            sumdb.list['Fields'] = sumfields
            if not self.getBool('MockBLAST'):
                sumdb.dropField('EVal')
                sumdb.renameField('MaxScore','MaxID')
            for sentry in sumdb.entries():
                sentry['HitNum'] = len(locdb.indexDataList('Qry',sentry['Qry'],'Hit'))
                sentry['CovHit'] = len(uniqdb.indexDataList('Qry',sentry['Qry'],'Hit'))
                if not self.getBool('MockBLAST'): sentry['MaxID'] = int(sentry['MaxID'])

            ### ~ [3] Add missing Queries from SeqList file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['Assembly']
            if not seqlist or not seqlist.seqNum(): self.warnLog('No assembly/seqin file loaded - cannot add missing queries to hitsum')
            else:
                for seq in seqlist.seqs():
                    sname = seqlist.shortName(seq)
                    sdesc = seqlist.seqDesc(seq)
                    if sname not in sumdb.dataKeys():
                        slen = seqlist.seqLen(seq)
                        sentry = {'Qry':sname,'HitNum':0,'MaxScore':0,'EVal':1000,'Description':sdesc,'Length':slen,
                                  'Coverage':0,'Identity':0,'Positives':0,
                                  'CovHit':0,'QryStart':0,'QryEnd':0,'Qry_AlnLen':0.0,'Qry_AlnID':0.0}
                        sumdb.addEntry(sentry)
                    else: sumdb.data(sname)['Description'] = sdesc

            ### ~ [4] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: sumdb.saveToFile(sfdict={'Qry_AlnLen':4,'Qry_AlnID':4})

            return sumdb
        except: self.errorLog('%s.uniqueToHitSum error' % self.prog()); return False
#########################################################################################################################
    def localToGABLAM(self,save=True):  ### Parse PAF Unique Local alignment tables into GABLAM table.
        '''Parse PAF Unique Local alignment tables into GABLAM table..'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('Forks') > 1: return self.forkLocalToGABLAM(save)
            #self.printLog('#~~#','## ~~~~~ PAF Local Alignment to GABLAM Stats ~~~~~ ##')
            self.headLog('PAF Local Alignment to GABLAM Stats')
            #gabfields = rje.split('Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID       Qry_AlnSim      Qry_Dirn        Qry_Start       Qry_End Qry_OrderedAlnLen       Qry_OrderedAlnID        Qry_OrderedAlnSim       Qry_OrderedDirn Qry_OrderedStart        Qry_OrderedEnd  Hit_AlnLen      Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End Hit_OrderedAlnLen       Hit_OrderedAlnID        Hit_OrderedAlnSim       Hit_OrderedDirn Hit_OrderedStart        Hit_OrderedEnd')
            gabfields = 'Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID Qry_AlnSim Qry_Dirn Qry_Start Qry_End Hit_AlnLen Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End'.split()
            locdb = self.db('local')
            gabdb = self.db().copyTable(locdb,'gablam')

            ### ~ [2] Compress local table copy for general stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gabdb.renameField('SbjLen','HitLen')
            gabdb.renameField('BitScore','Score')
            gabdb.renameField('Expect','Eval')
            gabdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Score':'sum'},default='max',best=[],joinchar='|')
            #gabdb.index('Qry',force=True)
            gabdb.rankFieldByIndex(index='Qry',rankfield='Identity',newfield='Rank',rev=True,absolute=True,unique=True)
            gabdb.dataFormat({'Rank':'int'})
            gabdb.newKey(['Qry','Rank'])
            gabdb.keepFields( 'Qry Hit Rank Score EVal QryLen HitLen'.split() )
            gabdb.list['Fields'] = gabfields

            ### ~ [3] Add GABLAM stats from Qry-Hit specifc unique mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = gabdb.entryNum(); filtx = 0
            for gentry in gabdb.entries():
                self.progLog('#GABLAM','Generating GABLAM table from local hits: %.2f%%' % (gx/gtot)); gx += 100
                #if not self.debugging():
                self.quiet()
                #self.log.opt['Silent'] = True
                qry = gentry['Qry']
                hit = gentry['Hit']
                if self.getBool('AlnSeq'):
                    refdb = self.reduceLocal(locdb,byqry=True,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
                else:
                    refdb = self.reduceLocalCS(locdb,byqry=True,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
                if not refdb.entryNum():    # Filtered by MinLocID
                    filtx += 1
                    self.talk()
                    gabdb.dropEntry(gentry)
                    continue

                #i# FUTURE UPDATE:
                #i# Also need to reduce unique hits to ordered ones for the ordered stats!
                #i# This will need to be a future addition
                refdb.makeField('QryEnd-QryStart+1','Coverage')
                refdb.dataFormat({'Coverage':'int'})
                refdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Coverage':'sum','QryStart':'min','Strand':'str'},default='max',best=[],joinchar='|')
                gentry['Qry_AlnLen'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Coverage'] / gentry['QryLen'] )
                gentry['Qry_AlnID'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Identity'] / gentry['QryLen'] )
                gentry['Qry_AlnSim'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Identity'] / gentry['QryLen'] )
                gentry['Qry_Dirn'] = 'Fwd'
                gentry['Qry_Start'] = refdb.data((qry,hit))['QryStart']
                gentry['Qry_End'] = refdb.data((qry,hit))['QryEnd']
                self.db().deleteTable(refdb)

                if self.getBool('AlnSeq'):
                    hitdb = self.reduceLocal(locdb,byqry=False,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
                else:
                    hitdb = self.reduceLocalCS(locdb,byqry=False,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
                if not hitdb.entryNum(): self.warnLog('Asymmetric loss of %s-%s alignments based on MinLocID' % (qry,hit),quitchoice=True,warntype='locfilter',suppress=True)
                for hentry in hitdb.entries():
                    if hentry['Strand'] == '-' and self.getBool('AlnSeq'): [ hentry['SbjStart'], hentry['SbjEnd'] ] = [ hentry['SbjEnd'], hentry['SbjStart'] ]
                hitdb.makeField('SbjEnd-SbjStart+1','Coverage')
                hitdb.dataFormat({'Coverage':'int'})
                hitdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Coverage':'sum','QryStart':'min','Strand':'str'},default='max',best=[],joinchar='|')
                gentry['Hit_AlnLen'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Coverage'] / gentry['HitLen'] )
                gentry['Hit_AlnID'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Identity'] / gentry['HitLen'] )
                gentry['Hit_AlnSim'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Identity'] / gentry['HitLen'] )
                gentry['Hit_Dirn'] = {'+':'Fwd','-':'Bwd'}[hitdb.data((qry,hit))['Strand']]
                gentry['Hit_Start'] = hitdb.data((qry,hit))['SbjStart']
                gentry['Hit_End'] = hitdb.data((qry,hit))['SbjEnd']
                self.db().deleteTable(hitdb)

                #self.log.opt['Silent'] = silent
                self.talk()
            self.printLog('#GABLAM','Generated GABLAM table from %s Qry-Hit pairs (%s MinLocID filtered).' % (rje.iStr(gx),rje.iStr(filtx)))

            ### ~ [4] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: gabdb.saveToFile()

            return gabdb
        except:
            #self.log.opt['Silent'] = silent
            self.talk()
            self.errorLog('%s.localToGABLAM error' % self.prog()); return False
#########################################################################################################################
    ### <3b> ### AlnSeq-based PAF Parsing Methods                                                                       #
#########################################################################################################################
    def pafToLocalAlnSeq(self,pafdb=None,save=True,alnseq=True,minlocid=0,minloclen=1):  ### Parse PAF table into Local alignment table.
        '''
        Parse PAF table into Local alignment table.
        >> pafin:file [None] = PAF file to load, or use
        >> minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','## ~~~~~ PAF to Local Alignment Conversion ~~~~~ ##')
            self.headLog('PAF to Local Alignment Conversion')
            if not pafdb: pafdb = self.db('paf')
            locdb = self.db().copyTable(pafdb,'local')
            #i# Was filtering here but PAF AlnLen/Length mismatches make this dangerous
            # if minloclen > 1: locdb.dropEntries('Length<%d' % minloclen)
            # if minlocid > 0:
            #     locdb.makeField('100.0*Identity/Length','LocID')
            #     locdb.dropEntries('LocID<%s' % minlocid)
            # self.deBug('?')
            locfields = 'Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd'.split()
            locfields += ['QryLen','SbjLen']

            ### ~ [2] Add AlnNum and re-index data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locdb.joinFields('Pair',['Qry','Hit'],join='|',replace=False)
            locdb.rankFieldByIndex(index='Pair',rankfield='Identity',newfield='AlnNum',rev=True,absolute=True,unique=True)
            locdb.newKey(['Qry','Hit','AlnNum'])
            locdb.makeField(formula='#Identity#',fieldname='Positives',evalue=-1)
            locdb.makeField(formula='#Identity#',fieldname='BitScore',evalue=-1)
            locdb.makeField(fieldname='Expect',evalue=-1)
            locdb.dataFormat({'Positives':'int','BitScore':'int','Expect':'num','AlnNum':'int'})
            locdb.keepFields(locfields+['Strand','cs'])
            locdb.list['Fields'] = locfields+['Strand','cs']

            self.printLog('#SEQ','Making sequence name dictionaries...',log=False)
            qryseqs = self.obj['Assembly']
            qrydict = qryseqs.makeSeqNameDic('short')
            sbjseqs = self.obj['Reference']
            sbjdict = sbjseqs.makeSeqNameDic('short')

            #!# NEED TO MAKE QrySeq, SbjSeq and AlnSeq from cs #!#
            for field in ['QrySeq','AlnSeq','SbjSeq']: locdb.addField(field)
            if self.getBool('LocalAln'): locfields += ['QrySeq','AlnSeq','SbjSeq']
            lx = 0.0; ltot = locdb.entryNum(); px = 0
            alncorrx = 0  # Number of AlnLen -> Length adjustments
            for lentry in locdb.entries():
                self.progLog('#ALN','Converting PAF cs alignments: %.1f%%' % (lx/ltot)); lx += 100.0
                aln = lentry['cs']
                if not aln: continue
                if aln.startswith('Z:'): aln = aln[2:]
                #self.bugPrint('|===|')
                #self.bugPrint(lentry)
                # EndExtend
                hang5 = lentry['QryStart'] - 1
                sbj5 = lentry['SbjStart'] - 1
                hang3 = lentry['QryLen'] - lentry['QryEnd']
                sbj3 = lentry['SbjLen'] - lentry['SbjEnd']

                rev = lentry['Strand'] == '-'
                #!# Move ? to other end if rev

                if lentry['Strand'] == '-': (sbj5,sbj3) = (sbj3,sbj5)
                if 0 < hang5 <= self.getInt('EndExtend'):
                    if hang5 > sbj5:
                        endgap = hang5 - sbj5
                        if rev: aln = '%s?%d!%d' % (aln,sbj5,endgap)
                        else: aln = '!%d?%d%s' % (endgap,sbj5,aln)
                    else:
                        if rev: aln = '%s?%d' % (aln,hang5)
                        else: aln = '?%d%s' % (hang5,aln)
                    lentry['QryStart'] -= hang5
                    if rev: lentry['SbjEnd'] += min(sbj5,hang5)
                    else: lentry['SbjStart'] -= min(sbj5,hang5)
                if 0 < hang3 <= self.getInt('EndExtend'):
                    if hang3 > sbj3:
                        endgap = hang3 - sbj3
                        if rev: aln = '!%d?%d%s' % (endgap,sbj3,aln)
                        else: aln = '%s?%d!%d' % (aln,sbj3,endgap)
                    else:
                        if rev: aln = '?%d%s' % (hang3,aln)
                        else: aln = '%s?%d' % (aln,hang3)
                    lentry['QryEnd'] += hang3
                    if lentry['Strand'] == '+': lentry['SbjEnd'] += min(sbj3,hang3)
                    else: lentry['SbjStart'] -= min(sbj3,hang3)
                lentry['cs'] = aln

                #if self.dev(): self.debug('%s => %s' % (aln,self.mapCS(lentry)))

                # Process
                for x in ':-+*~?!': aln = aln.replace(x,' %s' % x)
                #self.bugPrint('|===|')
                #self.bugPrint(lentry)
                aln = aln.split()
                #self.debug(aln)
                for field in ['QrySeq','AlnSeq','SbjSeq']: lentry[field] = ''
                qseq = sseq = aseq = ''
                qfull = qryseqs.seqSequence(qrydict[lentry['Qry']])[lentry['QryStart']-1:lentry['QryEnd']]
                sfull = sbjseqs.seqSequence(sbjdict[lentry['Hit']])[lentry['SbjStart']-1:lentry['SbjEnd']]
                if lentry['Strand'] == '-': qfull = rje_sequence.reverseComplement(qfull)
                qwarnx = 0; swarnx = 0
                for cs in aln:

                    #if '?' in lentry['cs'] or cs[0]=='*':
                    #    self.bugPrint('')
                    #    self.bugPrint(qfull[:100])
                    #    self.bugPrint(sfull[:100])
                    #    self.bugPrint(cs)
                    try:
                        if cs[0] == ':':    # Identity
                            ilen = rje.atoi(cs[1:])
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += '|' * ilen
                        elif cs[0] == '*':    # Mismatch
                            ilen = 1
                            try:
                                if cs[2] != qfull[0].lower(): qwarnx += 1
                                if cs[1] != sfull[0].lower(): swarnx += 1
                            except:
                                if not qfull: self.errorLog('CS parsing error: query sequence (%s) consumed before CS string (hit:%s)' % (lentry['Qry'],lentry['Hit']),printerror=False)
                                elif not sfull: self.errorLog('CS parsing error: hit sequence (%s) consumed before CS string (query:%s)' % (lentry['Hit'],lentry['Qry']),printerror=False)
                                else: self.errorLog('Problem parsing CS element "%s" (query:%s, hit:%s)' % (cs,lentry['Qry'],lentry['Hit']))
                                self.warnLog('CS parsing prematurely ended (query:%s, hit:%s)' % (lentry['Qry'],lentry['Hit']))
                                break
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += ' '
                        elif cs[0] == '-':    # Insertion in subject
                            ilen = len(cs[1:])
                            qseq += '-' * ilen
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += ' ' * ilen
                        elif cs[0] == '+':    # Deletion in subject
                            ilen = len(cs[1:])
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += '-' * ilen
                            aseq += ' ' * ilen
                        elif cs[0] == '~':    # Intron in subject
                            ilen = rje.atoi(cs[3:-2])
                            qseq += '-' * ilen
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += '~' * ilen
                        elif cs[0] == '!':    # Missing in subject
                            ilen = rje.atoi(cs[1:])
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += '^' * ilen
                            aseq += '^' * ilen
                        elif cs[0] == '?':    # Uncertain match
                            ilen = rje.atoi(cs[1:])
                            lentry['Length'] += ilen
                            qadd = qfull[:ilen]; qfull = qfull[ilen:]; qseq += qadd
                            sadd = sfull[:ilen]; sfull = sfull[ilen:]; sseq += sadd
                            for i,qnt in enumerate(qadd):
                                snt = sadd[i]
                                if qnt == snt: aseq += '|'; lentry['Identity'] += 1; lentry['Positives'] += 1
                                else: aseq += ' '
                        else: raise ValueError('Unexpected cs element: %s' % cs)
                    except:
                        self.bugPrint(qfull[:100])
                        self.bugPrint(sfull[:100])
                        self.bugPrint(cs)
                        self.debug(aln)
                        self.deBug(lentry)
                        self.errorLog('Unanticipated CS parsing error (query:%s, hit:%s)' % (lentry['Qry'],lentry['Hit']))
                        self.warnLog('CS parsing prematurely ended (query:%s, hit:%s)' % (lentry['Qry'],lentry['Hit']))
                        break

                #self.bugPrint(lentry)
                if qwarnx: self.warnLog('%s CS mismatch does not match loaded query %s sequence' % (qwarnx,lentry['Qry']),dev=self.debugging())
                if swarnx: self.warnLog('%s CS mismatch does not match loaded subject %s sequence' % (swarnx,lentry['Hit']),dev=self.debugging())
                if qfull: self.warnLog('Query sequence not consumed by alignment: %s' % qfull,dev=self.debugging())
                if sfull: self.warnLog('Subject sequence not consumed by alignment: %s' % sfull,dev=self.debugging())

                #if qwarnx or swarnx or qfull or sfull: self.debug(lentry)
                #else: self.debug('OK: %s' % lentry)

                #:37+ac:83+a:52+c:73*cg*gc:90-g:20-c:54-c:33+t:52+t:53-c:22-a:28*at*ta:15-at*at:25-aa:21+c:12-a:9-a:40+t:28-c:9*ac:64-ca:90-t:54+c:13-t:5-t:7-g:19-a:11-a:90+a:13+aa:14-c:52-g:7-c:7+a:51-t:22+t:11+t:11-t:44-a:12-ct*gt:66-c:62-g:23*ac:32+t:15+a:83-a:29+a:29+c:122+a:17+g:113-a:35+t:137-t:31+t:49+c:76+t:152+c:25+g:43+t:257+a:23+g:591+a:361-a:89+t:310+c:36*ca*at*tg:295+g:503+a:173-g:216+a:50-c:1332-a:1462-a:569-a:801-g:1716+g:34+a:657+a:55+a:30+a:240+a:2737-a:20062+a:685+c:10563+a:143986

                #i# Check for mismatched lengths
                alen = max(len(qseq),len(aseq),len(sseq))
                lentry['QrySeq'] = qseq + '?' * (len(qseq) - alen)
                lentry['AlnSeq'] = aseq + ' ' * (len(aseq) - alen)
                lentry['SbjSeq'] = sseq + '?' * (len(sseq) - alen)

                ## Sort out direction
                if lentry['Strand'] == '-':
                    [ lentry['SbjStart'], lentry['SbjEnd'] ] = [ lentry['SbjEnd'], lentry['SbjStart'] ]
                    lentry['QrySeq'] = rje_sequence.reverseComplement(lentry['QrySeq'])
                    lentry['SbjSeq'] = rje_sequence.reverseComplement(lentry['SbjSeq'])
                    lentry['AlnSeq'] = rje.strReverse(lentry['AlnSeq'])

                if not lentry['SbjSeq']:
                    self.debug(aln)
                    self.debug(lentry)
                px += 1

                alnlen = len(lentry['AlnSeq']) - lentry['AlnSeq'].count('^')  - lentry['AlnSeq'].count('~')
                if lentry['Length'] != alnlen:
                    self.debugging() and self.warnLog('Length/AlnLen mismatch for %s vs %s (Aln %s): %s vs %s [cs:%s]' % (lentry['Qry'],lentry['Hit'],lentry['AlnNum'],lentry['Length'],alnlen,lentry['cs']),
                                                      'alnlen',suppress=True,quitchoice=False)
                    #i# Check the start and end positions
                    qlen = len(lentry['QrySeq']) - lentry['QrySeq'].count('-')  - lentry['QrySeq'].count('^')
                    qend = min(lentry['QryStart'],lentry['QryEnd']) + qlen - 1
                    if qend != max(lentry['QryStart'],lentry['QryEnd']):
                        self.warnLog('QryAln length mismatch: Start:%s; End:%s QrySeq:%s' % (lentry['QryStart'],lentry['QryEnd'],qlen))
                    hlen = len(lentry['SbjSeq']) - lentry['SbjSeq'].count('-')  - lentry['SbjSeq'].count('^')
                    hend = min(lentry['SbjStart'],lentry['SbjEnd']) + hlen - 1
                    if hend != max(lentry['SbjStart'],lentry['SbjEnd']):
                        self.warnLog('SbjAln length mismatch: Start:%s; End:%s SbjSeq:%s' % (lentry['SbjStart'],lentry['SbjEnd'],hlen))
                    lentry['Length'] = alnlen
                    alncorrx += 1

                #self.debug('%d (paf) vs %d (aln)' % (lentry['Length'],alnlen))

            self.printLog('#ALN','Conversion of %s PAF cs alignments complete.' % rje.iStr(px))
            if alncorrx:
                #self.warnLog('%s Length/AlnLen mismatches from PAF file. (Seems to be Mimimap2 bug. Run debug=T for details.)' % rje.iStr(alncorrx))
                self.warnLog('%s minor Length/AlnLen mismatches from PAF file. (Run debug=T for details.)' % rje.iStr(alncorrx))
            if px != ltot and (alnseq or self.getBool('LocalAln')):
                raise ValueError('%s cs alignments missing from PAF file: check --cs flag was used' % rje.iStr(ltot-px))
            if minloclen > 1: locdb.dropEntries('Length<%d' % minloclen)
            if minlocid > 0:
                locdb.makeField('100.0*Identity/Length','LocID')
                locdb.dropEntries('LocID<%s' % minlocid)
                locdb.dropField('LocID')
            #self.deBug('?')

            #self.debug(locdb.data(('ctgIII_MBG11A__SP16533T7.03','ctgIII_MBG11A__SP16533T7.03',1)))

            ### ~ [4] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Add dropping of Positives at some point - will need to fix code elsewhere
            #if not self.getBool('MockBLAST'):
            #    locdb.dropField('Positives')
            #if self.debugging(): locdb.saveToFile()
            #locdb.keepFields(locfields+['Strand','cs'])
            if save: locdb.saveToFile(savefields=locfields+['Strand','cs'])

            return locdb
        except: self.errorLog('%s.pafToLocal error' % self.prog()); return False
#########################################################################################################################
    def localToUniqueAlnSeq(self,locdb=None,save=True):  ### Convert local alignments to unique coverage of queries and hits.
        '''Convert local alignments to unique coverage of queries and hits.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','## ~~~~~ PAF Local Alignment to Unique Hits ~~~~~ ##')
            self.headLog('PAF Local Alignment to Unique Hits')
            if not locdb: locdb = self.db('local')
            locfields = 'Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd'.split()

            ### ~ [2] Add AlnNum and re-index data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniqdb = self.reduceLocal(locdb,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
            uniqdb.rename('hitunique')
            qryuniqdb = self.reduceLocal(locdb,byqry=True,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
            qryuniqdb.rename('qryunique')

            ### ~ [3] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save:
                uniqdb.saveToFile(savefields=locfields)
                qryuniqdb.saveToFile(savefields=locfields)

            return uniqdb
        except: self.errorLog('%s.localToUnique error' % self.prog()); return False
#########################################################################################################################
    def forkLocalToGABLAM(self,save=True):  ### Parse PAF Local alignment tables into GABLAM table.
        '''
        Parse PAF Local alignment tables into GABLAM table, forking out each qry-hit pair (i.e. a single GABLAM entry)
        as a pair of forks - one reducing the hits to unique Query coverage for query stats, and one reducing hits to
        unique Hit coverage for Hit stats.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','## ~~~~~ Forking PAF Local Alignment to GABLAM Stats ~~~~~ ##')
            self.headLog('Forking PAF Local Alignment to GABLAM Stats')
            locdb = self.db('local')
            tmpdir = rje.makePath('%s.tmp/' % self.baseFile()) #x#self.getStr('TmpDir')
            rje.mkDir(self,tmpdir,log=self.dev() or self.debugging())
            self.bugPrint('%s: %s' % (tmpdir,rje.exists(tmpdir)))
            forkx = max(1,self.getInt('Forks'))
            ## ~ [1a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['GABTmp'] = []  # List of temporary files made during forking
            self.list['Forked'] = [None] * forkx     # List of (type,qry,hit) tuples forked
            self.list['ToFork'] = []  # List of (type,qry,hit) tuples to fork
            self.silence()
            qhdb = self.db().splitTable(locdb,'Qry',asdict=True,keepfield=True,add=False)
            #x#for qry in locdb.index('Qry'):
            for qry in rje.sortKeys(qhdb):
                qhdb[qry] = self.db().splitTable(qhdb[qry],'Hit',asdict=True,keepfield=True,add=False)
                #x#for hit in locdb.indexDataList('Qry',qry,'Hit',sortunique=True):
                for hit in rje.sortKeys(qhdb[qry]):
                    for qh in ['qry','hit']:
                        self.list['ToFork'].append((qh,qry,hit))
            self.talk()
            self.list['ToFork'].sort()
            self.setNum({'KillTime':time.time()})
            self.printLog('#FORK','%s Qry-Hit pairs for forking' % rje.iStr(len(self.list['ToFork'])/2))

            ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forktot = len(self.list['ToFork'])
            forked = 0
            finished = 0
            skipped = 0
            pidcheck = None
            if self.dev() or self.debugging(): pidcheck = '%s.pid' % self.baseFile()
            while self.list['Forked']:
                ## ~ [2a] Process finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if pidcheck: PIDCHECK = open(pidcheck,'w')
                for fdict in self.list['Forked'][0:]:
                    if fdict:
                        try:
                            pid = fdict['PID']
                            if pidcheck: PIDCHECK.write('%s: %s\n' % (self.list['Forked'].index(fdict),pid))
                            if str(pid).split()[0] == 'WAIT': status = 1
                            else: (status,exit_stat) = os.waitpid(pid,os.WNOHANG)
                        except:
                            self.errorLog('!')
                            status = 1
                    else: status = 1

                    if status > 0:
                        self.list['Forked'].remove(fdict)
                        self.setNum({'KillTime':time.time()})
                        if fdict:
                            if not rje.exists(fdict['File']):
                                self.warnLog('Fork finished but output file "%s" missing!' % fdict['File'])
                            finished += 1
                            self.progLog('#GABLAM','Forking %s local qry-hit GABLAM reductions: %.1f%% (%s forked; %s skipped).' % (rje.iStr(forktot),100.0*finished/forktot,rje.iStr(forked),rje.iStr(skipped)))
                        if not self.list['ToFork']: continue      # No more runs to fork
                        # Fork out GABLAM conversion
                        (qh,qry,hit) = self.list['ToFork'].pop(0)
                        gfile = '%s%s.%s.%s.tdt' % (tmpdir,rje.fileSafeString(qry),rje.fileSafeString(hit),qh)
                        self.list['GABTmp'].append(gfile)
                        if rje.exists(gfile) and not self.force():
                            finished += 1
                            skipped += 1
                            self.list['Forked'].append(None)
                        else:
                            forked += 1
                            fdict = {'File':gfile}
                            self.list['Forked'].append(fdict)
                            fdict['ID'] = 'Fork %d' % self.list['Forked'].index(fdict)
                            ### ~ [2] ~ Add Fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                            cpid = os.fork()        # Fork child process
                            if cpid:                # parent process records pid of child rsh process
                                fdict['PID'] = cpid
                            else:                   # child process
                                try:
                                    if self.debugging(): self.quiet()
                                    else: self.silence()
                                    #x#forkdb = self.reduceLocal(locdb,byqry=qh=='qry',queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'))
                                    if self.getBool('AlnSeq'):
                                        forkdb = self.reduceLocal(qhdb[qry][hit],byqry=qh=='qry',queries=[qry],hits=[hit],quiet=False,minloclen=self.getInt('MinLocLen'),minlocid=self.getNum('MinLocID'))
                                    else:
                                        forkdb = self.reduceLocalCS(qhdb[qry][hit],byqry=qh=='qry',queries=[qry],hits=[hit],quiet=False,minloclen=self.getInt('MinLocLen'),minlocid=self.getNum('MinLocID'))
                                    forkdb.saveToFile(gfile,backup=False)
                                except:
                                    self.errorLog('Something went wrong with forked GABLAM local hit reduction -> %s' % gfile)
                                self.talk()
                                os._exit(0)
                if pidcheck:
                    PIDCHECK.close()
                    #self.deBug(open(pidcheck,'r').read())

                ## ~ [2b] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if time.time() - self.getNum('KillTime') > self.getNum('KillForks'):
                    self.verbose(0,1,'\n%d seconds of main thread inactivity. %d forks still active!' % (self.getNum('KillForks'),len(self.list['Forked'])),1)
                    for fdict in self.list['Forked']:
                        self.verbose(0,2,' => Fork %s, PID %d still Active!' % (fdict['ID'],fdict['PID']),1)
                    if self.i() < 0 or rje.yesNo('Kill Main Thread?'):
                        raise ValueError('%d seconds of main thread inactivity. %d forks still active!' % (self.getNum('KillForks'),len(self.list['Forked'])))
                    elif rje.yesNo('Kill hanging forks?'):
                        for fdict in self.list['Forked']:
                            self.printLog('#KILL','Killing Fork %s, PID %d.' % (fdict['ID'],fdict['PID']))
                            os.system('kill %d' % fdict['PID'])
                    else: self.setNum({'KillTime':time.time()})
                ## ~ [2c] Sleep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(self.getNum('ForkSleep'))
            self.printLog('#GABLAM','Forked %s local qry-hit GABLAM reductions: %s complete (%s forked; %s skipped).' % (rje.iStr(forktot),rje.iStr(finished),rje.iStr(forked),rje.iStr(skipped)))

            ### ~ [3] ~ Load forked out files and process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gabfields = 'Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID Qry_AlnSim Qry_Dirn Qry_Start Qry_End Hit_AlnLen Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End'.split()
            gabdb = self.db().copyTable(locdb,'gablam')
            ## ~ [3b] Compress local table copy for general stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gabdb.renameField('SbjLen','HitLen')
            gabdb.renameField('BitScore','Score')
            gabdb.renameField('Expect','Eval')
            gabdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Score':'sum'},default='max',best=[],joinchar='|')
            #gabdb.index('Qry',force=True)
            gabdb.rankFieldByIndex(index='Qry',rankfield='Identity',newfield='Rank',rev=True,absolute=True,unique=True)
            gabdb.dataFormat({'Rank':'int'})
            gabdb.newKey(['Qry','Rank'])
            gabdb.keepFields( 'Qry Hit Rank Score EVal QryLen HitLen'.split() )
            gabdb.list['Fields'] = gabfields

            ### ~ [4] Add GABLAM stats from Qry-Hit specifc unique mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = gabdb.entryNum(); filtx = 0
            for gentry in gabdb.entries():
                self.progLog('#GABLAM','Generating GABLAM table from local hits: %.2f%%' % (gx/gtot)); gx += 100
                self.log.silence()
                #self.log.opt['Silent'] = True
                qry = gentry['Qry']
                hit = gentry['Hit']
                gfile = '%s%s.%s.qry.tdt' % (tmpdir,rje.fileSafeString(qry),rje.fileSafeString(hit))
                refdb = self.db().addTable(gfile,mainkeys=['Qry','Hit'],expect=True,name='tempqry')
                if not refdb or not refdb.entryNum():    # Filtered by MinLocID
                    self.talk()
                    if not refdb:
                        self.errorLog('Problem loading %s: no %s,%s GABLAM entry' % (gfile,qry,hit))
                    filtx += 1
                    gabdb.dropEntry(gentry)
                    self.db().deleteTable(refdb)
                    continue

                #i# FUTURE UPDATE:
                #i# Also need to reduce unique hits to ordered ones for the ordered stats!
                #i# This will need to be a future addition
                refdb.dataFormat(locformats)
                refdb.makeField('QryEnd-QryStart+1','Coverage')
                refdb.dataFormat({'Coverage':'int'})
                refdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Coverage':'sum','QryStart':'min','Strand':'str'},default='max',best=[],joinchar='|')
                gentry['Qry_AlnLen'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Coverage'] / gentry['QryLen'] )
                gentry['Qry_AlnID'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Identity'] / gentry['QryLen'] )
                gentry['Qry_AlnSim'] = '%.2f' % (100.0 * refdb.data((qry,hit))['Identity'] / gentry['QryLen'] )
                gentry['Qry_Dirn'] = 'Fwd'
                gentry['Qry_Start'] = refdb.data((qry,hit))['QryStart']
                gentry['Qry_End'] = refdb.data((qry,hit))['QryEnd']
                self.db().deleteTable(refdb)

                gfile = '%s%s.%s.hit.tdt' % (tmpdir,rje.fileSafeString(qry),rje.fileSafeString(hit))
                hitdb = self.db().addTable(gfile,mainkeys=['Qry','Hit'],expect=True,name='temphit')
                if not hitdb or not hitdb.entryNum():
                    self.talk()
                    self.warnLog('Asymmetric loss of %s-%s alignments based on MinLocID' % (qry,hit),quitchoice=True,warntype='locfilter',suppress=True)
                    if not hitdb:
                        self.errorLog('Problem loading %s: no %s,%s GABLAM entry' % (gfile,qry,hit))
                        self.db().deleteTable(hitdb)
                        continue
                hitdb.dataFormat(locformats)
                for hentry in hitdb.entries():
                    if hentry['Strand'] == '-' and self.getBool('AlnSeq'): [ hentry['SbjStart'], hentry['SbjEnd'] ] = [ hentry['SbjEnd'], hentry['SbjStart'] ]
                hitdb.makeField('SbjEnd-SbjStart+1','Coverage')
                hitdb.dataFormat({'Coverage':'int'})
                hitdb.compress(newkeys=['Qry','Hit'],rules={'Identity':'sum','Coverage':'sum','QryStart':'min','Strand':'str'},default='max',best=[],joinchar='|')
                gentry['Hit_AlnLen'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Coverage'] / gentry['HitLen'] )
                gentry['Hit_AlnID'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Identity'] / gentry['HitLen'] )
                gentry['Hit_AlnSim'] = '%.2f' % (100.0 * hitdb.data((qry,hit))['Identity'] / gentry['HitLen'] )
                gentry['Hit_Dirn'] = {'+':'Fwd','-':'Bwd'}[hitdb.data((qry,hit))['Strand']]
                gentry['Hit_Start'] = hitdb.data((qry,hit))['SbjStart']
                gentry['Hit_End'] = hitdb.data((qry,hit))['SbjEnd']
                self.db().deleteTable(hitdb)

                #self.log.opt['Silent'] = silent
                self.talk()
            self.talk()
            self.printLog('#GABLAM','Generated GABLAM table from %s Qry-Hit pairs (%s MinLocID filtered).' % (rje.iStr(gtot),rje.iStr(filtx)))
            #self.debug(gabdb.entryNum())

            ### ~ [4] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: gabdb.saveToFile()
            if not self.dev() and not self.debugging():
                missx = 0
                for gfile in self.list['GABTmp']:
                    if rje.exists(gfile): os.unlink(gfile)
                    else: missx += 1
                try: os.rmdir(tmpdir)
                except: self.errorLog('Problem deleting %s' % tmpdir)
                self.printLog('#DEL','%s temporary GABLAM files deleted from %s (%s missing)' % (rje.iLen(self.list['GABTmp']),tmpdir,rje.iStr(missx)))

            return gabdb
        except SystemExit: raise    # Child
        except:
            #self.log.opt['Silent'] = silent
            self.talk()
            self.errorLog('%s.localToGABLAM error' % self.prog()); return False
#########################################################################################################################
    ### <4> ### BLAST Trimming Methods                                                                                  #
#########################################################################################################################
    def reduceLocal(self,locdb=None,queries=[],hits=[],sortfield='Identity',keepself=False,minloclen=1,minlocid=0,byqry=False,quiet=False):    ### Reduces local BLAST alignments to cover each hit region only once
        '''
        Reduces local BLAST alignments to cover each hit region only once. Local alignments are sorted by identity
        (unless sortfield=X changed) and processed in order. Other local alignments that overlap are truncated and
        updated. Any alignments completely overlapped are removed. Processes and returns a COPY of the table.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param sortfield:str ['Identity'] = LocalDB field used to sort local alignments.
        @param keepself:bool [False] = Whether to include self query-hit pairs in assessment.
        @param minloclen:int [0] = Minimum local length to keep.
        @param minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        @param byqry:bool [False] = Whether to reduce Local table by Query rather than hit
        @param quiet:bool [False] = Whether to suppress log/stdout progress messages.
        @return: copy of local table, filtered and reduced.
        '''
        silent = self.log.opt['Silent']
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.opt['Silent'] = silent or quiet
            if not locdb: locdb = self.db('Local')
            ## ~ [0a] Invert and run if byqry=True ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if byqry:
                #self.printLog('#~~#','## ~~~~~ Inverted PAF Local Alignment Reduction Preparation ~~~~~ ##')
                self.headLog('Inverted PAF Local Alignment Reduction Preparation')
                # Create Query-unique tables by inverting Qry and Sbj
                refdb = self.db().copyTable(locdb,'qry.%s' % locdb.name())
                refdb.dataFormat({'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'})
                # Invert the Query and Hit data
                self.printLog('#INVERT','Inverting local BLAST hits for qryunique=T reduction')
                for entry in refdb.entries():
                    #self.deBug(entry)
                    [entry['Qry'],entry['Hit'],entry['QryStart'],entry['QryEnd'],entry['QrySeq'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq']] = [entry['Hit'],entry['Qry'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq'],entry['QryStart'],entry['QryEnd'],entry['QrySeq']]
                    if entry['QryStart'] > entry['QryEnd']:
                        [entry['QryStart'],entry['QryEnd'],entry['SbjStart'],entry['SbjEnd']] = [entry['QryEnd'],entry['QryStart'],entry['SbjEnd'],entry['SbjStart']]
                        entry['AlnSeq'] = entry['AlnSeq'][::-1]
                        entry['QrySeq'] = rje_sequence.reverseComplement(entry['QrySeq'])
                        entry['SbjSeq'] = rje_sequence.reverseComplement(entry['SbjSeq'])
                #self.debug(refdb.indexKeys('Qry'))
                refdb.remakeKeys()
                refudb = self.reduceLocal(refdb,hits,queries,sortfield,keepself,minloclen,minlocid)
                self.db().deleteTable(refdb)
                if not refudb: raise ValueError('Reference BLAST.reduceLocal() failed.')
                #self.debug(refudb.indexKeys('Qry'))
                # Re-invert the Qry and Hit data
                self.printLog('#INVERT','Inverting local BLAST hits following qryunique=T reduction')
                for entry in refudb.entries():
                    [entry['Qry'],entry['Hit'],entry['QryStart'],entry['QryEnd'],entry['QrySeq'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq']] = [entry['Hit'],entry['Qry'],entry['SbjStart'],entry['SbjEnd'],entry['SbjSeq'],entry['QryStart'],entry['QryEnd'],entry['QrySeq']]
                    if entry['QryStart'] > entry['QryEnd']:
                        [entry['QryStart'],entry['QryEnd'],entry['SbjStart'],entry['SbjEnd']] = [entry['QryEnd'],entry['QryStart'],entry['SbjEnd'],entry['SbjStart']]
                        entry['QrySeq'] = rje_sequence.reverseComplement(entry['QrySeq'])
                        entry['SbjSeq'] = rje_sequence.reverseComplement(entry['SbjSeq'])
                        entry['AlnSeq'] = entry['AlnSeq'][::-1]
                refudb.remakeKeys()
                refudb.setStr({'Name':refudb.name()[4:]})
                #self.debug(refudb.indexKeys('Qry'))
                self.log.opt['Silent'] = silent
                return refudb
            ## ~ [0b] Setup local hit reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.printLog('#~~#','## ~~~~~ PAF Local Alignment Reduction ~~~~~ ##')
            self.headLog('PAF Local Alignment Reduction')
            if sortfield not in locdb.fields():
                raise ValueError('"%s" is not a valid Local hit table field (%s)' % (sortfield,rje.join(locdb.fields(),'|')))
            # Create a copy to protect initial data
            bdb = self.db().copyTable(locdb,'%s.reduced' % locdb.name())
            # ['Qry','Hit','AlnNum','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            # Filter if required
            if queries: bdb.dropEntriesDirect('Qry',queries,inverse=True)
            if hits: bdb.dropEntriesDirect('Hit',hits,inverse=True)
            if not keepself: bdb.dropEntries(['Qry=Hit'])
            bdb.dataFormat({'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int','QryLen':'int','HitLen':'int'})
            btot = bdb.entryNum()
            bdb.dropEntries(['Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing short local hits')

            if minlocid > 0.0:
                badidx = 0
                for entry in bdb.entries()[0:]:
                    if 100.0 * float(entry['Identity']) / int(entry['Length']) < minlocid:
                        badidx += 1
                        bdb.dropEntry(entry)
                self.printLog('#MINID','Dropped %s entries < LocalIDMin=%s%%' % (rje.iStr(badidx),rje.sf(minlocid)))

            mx = btot - bdb.entryNum()
            ### ~ [1] Cycle through local alignments, reducing as required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bentries = bdb.sortedEntries(sortfield,reverse=sortfield in ['Identity','Positives','Length','BitScore'])   # List of all entries (sorted) to process
            alignpos = {}; ax = 0   # Dictionary of {Hit:[(start,stop) list of positions included in local aln]}
            while bentries:
                ## ~ [1a] Grab next best remaining hit from bentries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = bentries.pop(0)     # This is best remaining hit
                #self.debug('QrySeq: %d; AlnSeq: %d; SbjSeq: %d' % (len(entry['QrySeq']),len(entry['AlnSeq']),len(entry['SbjSeq'])))
                ax += 1
                self.progLog('\r#LOCALN','Processing local alignments: %s -> %s' % (rje.iLen(bentries),rje.iStr(ax)))
                hit = entry['Hit']
                region = (min(entry['SbjStart'],entry['SbjEnd']),max(entry['SbjStart'],entry['SbjEnd']))
                ## ~ [1b] Update alignpos dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if hit not in alignpos: alignpos[hit] = []
                alignpos[hit].append(region)
                alignpos[hit] = rje.collapseTupleList(alignpos[hit])
                ## ~ [1c] Adjust/Filter remaining entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ex = 0
                while ex < len(bentries):
                    xentry = bentries[ex]
                    yentry = None   # Will be created if splitting
                    # Skip if different Hit
                    if entry['Hit'] != xentry['Hit']: ex += 1; continue
                    # Check for overlapping regions and remove
                    xregion = (min(xentry['SbjStart'],xentry['SbjEnd']),max(xentry['SbjStart'],xentry['SbjEnd']))
                    # No overlap
                    if xregion[1] < region[0] or xregion[0] > region[1]: ex += 1; continue
                    # Completely overlapped: remove
                    elif xregion[0] >= region[0] and xregion[1] <= region[1]:
                        bdb.dropEntry(xentry)
                        bentries.pop(ex)
                        continue
                    # Middle covered: split
                    elif region[0] > xregion[0] and region[1] < xregion[1]:
                        #self.bugPrint('\nEntry splitting: %s vs %s' % (xregion,region))
                        xalnx = max(bdb.indexDataList('Hit',hit,'AlnNum'))
                        yentry = rje.combineDict({'AlnNum':xalnx+1},xentry,overwrite=False)
                        self.bugLog('#ALNID','Splitting %s vs %s Aln %d -> %d & %d' % (xentry['Qry'],xentry['Hit'],xentry['AlnNum'],xentry['AlnNum'],yentry['AlnNum']))
                        #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},xentry,overwrite=False,replaceblanks=False))
                        self.trimLocal(xentry,trimend='End',trimto=region[0],sortends=True)   # Trim the end back to region[0]
                        #self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},xentry,overwrite=False,replaceblanks=False))
                        #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},yentry,overwrite=False,replaceblanks=False))
                        self.trimLocal(yentry,trimend='Start',trimto=region[1],sortends=True)   # Trim the start back to region[1]
                        #self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},yentry,overwrite=False,replaceblanks=False))
                    # Overlap at one end
                    elif region[0] <= xregion[1] <= region[1]:  # End overlaps with focal entry
                        #self.bugPrint('\nEnd overlap: %s vs %s' % (xregion,region))
                        self.trimLocal(xentry,trimend='End',trimto=region[0],sortends=True)   # Trim the end back to region[0]
                    elif region[0] <= xregion[0] <= region[1]:  # Start overlaps with focal entry
                        #self.bugPrint('\nStart overlap: %s vs %s' % (xregion,region))
                        self.trimLocal(xentry,trimend='Start',trimto=region[1],sortends=True)   # Trim the start back to region[1]
                    else: raise ValueError('Entry filtering has gone wrong: %s vs %s' % (xregion,region))
                    ## Check lengths
                    if xentry['Length'] >= minloclen: ex += 1
                    else: bdb.dropEntry(xentry); bentries.pop(ex); mx += 1
                    if yentry and yentry['Length'] >= minloclen:
                        bdb.addEntry(yentry)
                        bentries.insert(ex,yentry)
                        ex += 1
                    elif yentry: mx += 1
            ### ~ [2] Check and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#LOCALN','Reduced local alignments to unique coverage: %s -> %s (%s failed to meet minloclen=%d)' % (rje.iStr(btot),rje.iStr(ax),rje.iStr(mx),minloclen))
            if ax != bdb.entryNum(): raise ValueError('EntryNum mismatch following reduceLocal(): %s best entries but %s local alignments' % (ax,bdb.entryNum()))
            self.log.opt['Silent'] = silent
            return bdb
        except:
            self.log.opt['Silent'] = silent
            self.errorLog('Problem during BLASTRun.reduceLocal()')
            return None
#########################################################################################################################
    def reduceLocalCS(self,locdb=None,queries=[],hits=[],sortfield='Identity',keepself=False,minloclen=1,minlocid=0,byqry=False,quiet=False):    ### Reduces local BLAST alignments to cover each hit region only once
        '''
        Reduces local BLAST alignments to cover each hit region only once. Local alignments are sorted by identity
        (unless sortfield=X changed) and processed in order. Other local alignments that overlap are truncated and
        updated. Any alignments completely overlapped are removed. Processes and returns a COPY of the table.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param sortfield:str ['Identity'] = LocalDB field used to sort local alignments.
        @param keepself:bool [False] = Whether to include self query-hit pairs in assessment.
        @param minloclen:int [0] = Minimum local length to keep.
        @param minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        @param byqry:bool [False] = Whether to reduce Local table by Query rather than hit
        @param quiet:bool [False] = Whether to suppress log/stdout progress messages.
        @return: copy of local table, filtered and reduced.
        '''
        silent = self.log.opt['Silent']
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.opt['Silent'] = silent or quiet
            if not locdb: locdb = self.db('Local')
            if byqry: target = 'Qry'
            else: target = 'Hit'
            ## ~ [0a] Setup local hit reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog('PAF Local Alignment Reduction')
            if sortfield not in locdb.fields():
                raise ValueError('"%s" is not a valid Local hit table field (%s)' % (sortfield,rje.join(locdb.fields(),'|')))
            # Create a copy to protect initial data
            bdb = self.db().copyTable(locdb,'%s.reduced' % locdb.name())
            bdb.dataFormat({'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int','QryLen':'int','HitLen':'int'})
            # ['Qry','Hit','AlnNum','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            ## ~ [0b] Filter if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if queries: bdb.dropEntriesDirect('Qry',queries,inverse=True)
            if hits: bdb.dropEntriesDirect('Hit',hits,inverse=True)
            if not keepself: bdb.dropEntries(['Qry=Hit'])
            btot = bdb.entryNum()
            # MinLocLen
            bdb.dropEntries(['Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing short local hits')
            # MinLocID
            if minlocid > 0.0:
                badidx = 0
                for entry in bdb.entries()[0:]:
                    if 100.0 * float(entry['Identity']) / int(entry['Length']) < minlocid:
                        badidx += 1
                        bdb.dropEntry(entry)
                self.printLog('#MINID','Dropped %s entries < LocalIDMin=%s%%' % (rje.iStr(badidx),rje.sf(minlocid)))
            # Removed count
            mx = btot - bdb.entryNum()

            ### ~ [1] Cycle through local alignments, reducing as required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bentries = bdb.sortedEntries(sortfield,reverse=sortfield in ['Identity','Positives','Length','BitScore'])   # List of all entries (sorted) to process
            alignpos = {}; ax = 0   # Dictionary of {Hit:[(start,stop) list of positions included in local aln]}
            while bentries:
                ## ~ [1a] Grab next best remaining hit from bentries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = bentries.pop(0)     # This is best remaining hit
                ax += 1
                self.progLog('\r#LOCALN','Processing local alignment CS: %s -> %s' % (rje.iLen(bentries),rje.iStr(ax)))
                hit = entry[target]
                if byqry:
                    region = (min(entry['QryStart'],entry['QryEnd']),max(entry['QryStart'],entry['QryEnd']))
                else:
                    region = (min(entry['SbjStart'],entry['SbjEnd']),max(entry['SbjStart'],entry['SbjEnd']))
                ## ~ [1b] Update alignpos dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if hit not in alignpos: alignpos[hit] = []
                alignpos[hit].append(region)
                alignpos[hit] = rje.collapseTupleList(alignpos[hit])
                ## ~ [1c] Adjust/Filter remaining entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ex = 0
                while ex < len(bentries):
                    xentry = bentries[ex]
                    if xentry and xentry['SbjStart'] > xentry['SbjEnd']: self.debug('\n\nX: %s' % xentry)
                    yentry = None   # Will be created if splitting
                    # Skip if different Hit
                    if entry[target] != xentry[target]: ex += 1; continue
                    # Check for overlapping regions and remove
                    if byqry: xregion = (min(xentry['QryStart'],xentry['QryEnd']),max(xentry['QryStart'],xentry['QryEnd']))
                    else: xregion = (min(xentry['SbjStart'],xentry['SbjEnd']),max(xentry['SbjStart'],xentry['SbjEnd']))
                    # No overlap
                    if xregion[1] < region[0] or xregion[0] > region[1]: ex += 1; continue
                    # Completely overlapped: remove
                    elif xregion[0] >= region[0] and xregion[1] <= region[1]:
                        bdb.dropEntry(xentry)
                        bentries.pop(ex)
                        continue
                    # Middle covered: split
                    elif region[0] > xregion[0] and region[1] < xregion[1]:
                        #self.bugPrint('\nEntry splitting: %s vs %s' % (xregion,region))
                        if not hit in bdb.index(target):
                            self.warnLog('Cannot find %s=%s entries!' % (target,hit))
                            self.bugPrint('%s' % bdb.index(target))
                            self.debug(bdb.entrySummary(xentry,collapse=True))
                        xalnx = max(bdb.indexDataList(target,hit,'AlnNum'))
                        yentry = rje.combineDict({'AlnNum':xalnx+1},xentry,overwrite=False)
                        self.bugLog('#ALNID','Splitting %s vs %s Aln %d -> %d & %d' % (xentry['Qry'],xentry['Hit'],xentry['AlnNum'],xentry['AlnNum'],yentry['AlnNum']))
                        xentry = self.trimCS(xentry,newend=region[0],target=target)
                        yentry = self.trimCS(yentry,newstart=region[1],target=target)
                    # Overlap at one end
                    elif region[0] <= xregion[1] <= region[1]:  # End overlaps with focal entry
                        xentry = self.trimCS(xentry,newend=region[0],target=target)
                    elif region[0] <= xregion[0] <= region[1]:  # Start overlaps with focal entry
                        xentry = self.trimCS(xentry,newstart=region[1],target=target)
                    else: raise ValueError('Entry filtering has gone wrong: %s vs %s' % (xregion,region))
                    ## Check lengths
                    if xentry and xentry['Length'] >= minloclen: ex += 1
                    else: bdb.dropEntry(bentries[ex]); bentries.pop(ex); mx += 1
                    if yentry and yentry['Length'] >= minloclen:
                        bdb.addEntry(yentry)
                        bentries.insert(ex,yentry)
                        ex += 1
                    elif yentry: mx += 1
                    # Debug
                    if xentry and xentry['SbjStart'] > xentry['SbjEnd']: self.debug('\n\nX: %s' % xentry)
                    if yentry and yentry['SbjStart'] > yentry['SbjEnd']: self.debug('Y: %s' % yentry)
            ### ~ [2] Check and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in bdb.indexEntries('Strand','-'):
                [entry['SbjStart'], entry['SbjEnd']] = [entry['SbjEnd'], entry['SbjStart']]
            self.printLog('\r#LOCALN','Reduced local alignment CS to unique coverage: %s -> %s (%s failed to meet minloclen=%d)' % (rje.iStr(btot),rje.iStr(ax),rje.iStr(mx),minloclen))
            if ax != bdb.entryNum(): raise ValueError('EntryNum mismatch following reduceLocal(): %s best entries but %s local alignments' % (ax,bdb.entryNum()))
            self.log.opt['Silent'] = silent
            return bdb
        except:
            self.log.opt['Silent'] = silent
            self.errorLog('Problem during BLASTRun.reduceLocalCS()')
            return None
#########################################################################################################################
    def trimLocal(self,lentry,trimend,trimto,sortends=True,debug=False):  # Trims local alignment entry data to hit coordinates
        '''
        Trims local alignment entry data.
        @param lentry: local alignment entry
        @param trimend: which end to trim (Start/End)
        @param trimto: position to trim up to (and including)
        @param sortends: whether to switch Start/End if reversed match is
        @return: modified lentry. (NOTE: Modified in place.)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            for field in ['QrySeq','SbjSeq','AlnSeq']:
                if not lentry[field]:
                    self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
                    raise ValueError('Empty %s!' % field)
            fwd = lentry['SbjStart'] <= lentry['SbjEnd']
            if sortends and not fwd: trimend = {'Start':'End','End':'Start'}[trimend]
            #if debug: self.bugPrint('Fwd:%s => %s' % (fwd,trimend))
            # Check need to trim
            if fwd and not lentry['SbjStart'] <= trimto <= lentry['SbjEnd']:
                self.warnLog('Fwd trimLocal() called but Start <= trimpos <= End not met.','trimlocal',suppress=True)
                return lentry
            if not fwd and not lentry['SbjStart'] >= trimto >= lentry['SbjEnd']:
                self.warnLog('Bwd trimLocal() called but End <= trimpos <= Start not met.','trimlocal',suppress=True)
                return lentry
            # Sort out starting positions in Query/Subject
            #qpos = lentry['QryStart'] #lentry['Qry%s' % trimend]
            spos = lentry['SbjStart']   #lentry['Sbj%s' % trimend]
            # Sort out starting position in alignment
            ai = 0
            if fwd: sdir = 1
            else: sdir = -1
            if trimend == 'Start': trimto += sdir  # For Start trim, want to go beyond trimto for alignment split
            ### ~ [1] Trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Locate trimto ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            #self.bugPrint(lentry)
            #self.bugPrint('Fwd:%s Scan %d pos until %d reaches %s' % (fwd,len(lentry['SbjSeq']),spos,trimto))
            while (fwd and spos < trimto) or (not fwd and spos > trimto):
                ai += 1
                try:
                    if lentry['SbjSeq'][ai] != '-': spos += sdir
                except:
                    wdict = rje.combineDict({'QrySeq':len(lentry['QrySeq']),'SbjSeq':len(lentry['SbjSeq']),'AlnSeq':len(lentry['AlnSeq'])},lentry,overwrite=False,replaceblanks=False)
                    if fwd:
                        self.warnLog('Unexpectedly reached end of SbjSeq during fwd local reduction trimming to %d: %s; Dropped alignment' % (trimto,wdict))
                    else:
                        self.warnLog('Unexpectedly reached end of SbjSeq during rev local reduction trimming to %d: %s; Dropped alignment' % (trimto,wdict))
                    lentry['Length'] = -1
                    return lentry
            #self.bugPrint('Reached AlnPos %d = SbjPos %d => Trim %s' % (ai,spos,trimend))
            ## ~ [1b] Split sequence and update stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if trimend == 'Start':  # Removing Start
                lentry['QryStart'] += len(lentry['QrySeq'][:ai]) - lentry['QrySeq'][:ai].count('-')
                lentry['QrySeq'] = lentry['QrySeq'][ai:]
                lentry['SbjStart'] += sdir * (len(lentry['SbjSeq'][:ai]) - lentry['SbjSeq'][:ai].count('-'))
                lentry['SbjSeq'] = lentry['SbjSeq'][ai:]
                if lentry['Positives']:
                    lentry['Positives'] -= lentry['AlnSeq'][:ai].count('+')
                    lentry['Positives'] -= lentry['AlnSeq'][:ai].count('|')
                lentry['Identity'] -= lentry['AlnSeq'][:ai].count('|')
                #lentry['Length'] -= len(lentry['AlnSeq'][:ai])
                lentry['AlnSeq'] = lentry['AlnSeq'][ai:]
                lentry['Length'] = len(lentry['AlnSeq']) - lentry['AlnSeq'].count('^')  - lentry['AlnSeq'].count('~')
            else:   # Removing end
                lentry['QryEnd'] -= (len(lentry['QrySeq'][ai:]) - lentry['QrySeq'][ai:].count('-'))
                lentry['QrySeq'] = lentry['QrySeq'][:ai]
                lentry['SbjEnd'] -= sdir * (len(lentry['SbjSeq'][ai:]) - lentry['SbjSeq'][ai:].count('-'))
                lentry['SbjSeq'] = lentry['SbjSeq'][:ai]
                if lentry['Positives']:
                    lentry['Positives'] -= lentry['AlnSeq'][ai:].count('+')
                    lentry['Positives'] -= lentry['AlnSeq'][ai:].count('|')
                lentry['Identity'] -= lentry['AlnSeq'][ai:].count('|')
                #lentry['Length'] -= len(lentry['AlnSeq'][ai:])
                lentry['AlnSeq'] = lentry['AlnSeq'][:ai]
                lentry['Length'] = len(lentry['AlnSeq']) - lentry['AlnSeq'].count('^')  - lentry['AlnSeq'].count('~')
            ### ~ [2] Return trimmed entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            return lentry
        except:
            self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            self.errorLog('BLASTRun.trimLocal error'); raise
#########################################################################################################################
    ### <5> ### Running Minimap2                                                                                        #
#########################################################################################################################
    def minimap2(self,mapopt=None,save=True,invert=False):  ### Run minimap2 and return or save PAF.
        '''Parse PAF Local alignment table into HitSum table.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SeqIn'):
                raise ValueError('SeqIn missing: %s' % rje.argString(self.cmd_list))
            if not self.getStrLC('Reference'):
                raise ValueError('Reference missing: %s' % rje.argString(self.cmd_list))
            if not mapopt: mapopt = self.dict['MapOpt']
            maprun = self.getStr('Minimap2') + ' --cs'
            for arg in mapopt:
                if arg == '-cs': continue
                if mapopt[arg] not in ['',None]: maprun += ' -%s %s' % (arg,mapopt[arg])
                else: maprun += ' -%s' % (arg)
            if invert: maprun += ' %s %s' % (self.getStr('SeqIn'),self.getStr('Reference'))
            else: maprun += ' %s %s' % (self.getStr('Reference'),self.getStr('SeqIn'))
            ### ~ [2] Run Minimap2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save:
                if not self.getStrLC('PAFIn'): self.setStr({'PAFIn':'%s.paf' % self.baseFile()})
                if rje.exists(self.getStr('PAFIn')) and not self.force():
                    self.printLog('#PAF','%s found and force=F: minimap2 run skipped.' % self.getStr('PAFIn'))
                else:
                    try:
                        mmv = os.popen('%s --version' % self.getStr('Minimap2')).read()
                        if not mmv: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting')
                        open('%s.cmd' % self.getStr('PAFIn'),'w').write('#minimap2 v%s' % mmv)
                    except: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % self.getStr('Minimap2'))
                    self.printLog('#SYS','%s > %s' % (maprun,self.getStr('PAFIn')))
                    open('%s.cmd' % self.getStr('PAFIn'),'a').write('%s\n' % maprun)
                    os.system('%s > %s' % (maprun,self.getStr('PAFIn')))
                return rje.exists(self.getStr('PAFIn'))
            else:
                return os.popen(maprun).read()
        except: self.errorLog('%s.minimap2 error' % self.prog()); raise
#########################################################################################################################
    ### <6> ### PAF CS string manipulations                                                                             #
#########################################################################################################################
# Op	Regex	Description
# =	[ACGTN]+	Identical sequence (long form)
# :	[0-9]+	Identical sequence length
# *	[acgtn][acgtn]	Substitution: ref to query
# +	[acgtn]+	Insertion to the reference
# -	[acgtn]+	Deletion from the reference
# ~	[acgtn]{2}[0-9]+[acgtn]{2}	Intron length and splice signal

# This module also adds ! [0-9]+ for extended matches of unknown type
# NOTE: The cs tag (as parsed) is based on the target sequence, such that a -ve strand hit will be the reverse complement
#       of the true hit and always has target start < end.
    #i# alnSeqToCS => convert QrySeq, AlnSeq, SbjSeq to a CS string
    def alnSeqToCS(self):
        return
    #i# alnSeqToGstring
    #i# CSToGstring -> produces pair
    #># generate a new GABLAM string of "Gstring" of start:[|+*-]\d+ list:end  >> paf2GString (Qry & Sbj pair)
    #># this is actually generated as part statsFromCS(cs)
    #i# GstringToCS -> cs from pair (is this needed?)
#########################################################################################################################
    def combineGstring(self,gstr1,gstr2):   ### Combine two Gstrings into the best GABLAM combination
        '''
        Combine two Gstrings into the best GABLAM combination.
        :param gstr1 [str]: First GString to combine
        :param gstr2 [str]: Second GString to combine
        :return: gstr [str]: Combined GString with the best combination
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdata1 = gstr1.split(':')
            if not len(gdata1) == 3: raise ValueError('Problem with GString1 format: %s' % gstr1)
            gdata2 = gstr2.split(':')
            if not len(gdata2) == 3: raise ValueError('Problem with GString2 format: %s' % gstr1)
            gstr1 = gdata1[1]
            gstr2 = gdata2[1]
            for x in '|+*-':
                gstr1 = gstr1.replace(x,' %s' % x)
                gstr2 = gstr2.replace(x,' %s' % x)
            gstr1 = gstr1.split()
            gstr2 = gstr2.split()
            # x=start; y=end; t=type; n=length
            x1 = rje.atoi(gdata1[0])
            x2 = rje.atoi(gdata2[0])
            gstr = '%d:' % min(x1,x2)
            t1 = gstr1[0][0]; n1 = rje.atoi(gstr1[0][1:]); y1 = x1 + n1 - 1
            t2 = gstr2[0][0]; n2 = rje.atoi(gstr2[0][1:]); y2 = x2 + n2 - 1
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while gstr1 and gstr2:
                # No overlap
                if y1 < x2:
                    gstr += gstr1.pop(0)
                    x1 = y1 + 1
                    if gstr1:
                        t1 = gstr1[0][0]; n1 = rje.atoi(gstr1[0][1:]); y1 = x1 + n1 - 1
                    continue
                if y2 < x1:
                    gstr += gstr2.pop(0)
                    x2 = y2 + 1
                    if gstr2:
                        t2 = gstr2[0][0]; n2 = rje.atoi(gstr2[0][1:]); y2 = x2 + n2 - 1
                    continue
                # Perfect match
                if x1 == x2 and y1 == y2:
                    if '|+*-'.index(t2) < '|+*-'.index(t1): gstr += gstr1.pop(0); gstr2.pop(0)
                    else: gstr += gstr2.pop(0); gstr1.pop(0)
                    x1 = y1 + 1
                    if gstr1:
                        t1 = gstr1[0][0]; n1 = rje.atoi(gstr1[0][1:]); y1 = x1 + n1 - 1
                    x2 = y2 + 1
                    if gstr2:
                        t2 = gstr2[0][0]; n2 = rje.atoi(gstr2[0][1:]); y2 = x2 + n2 - 1
                    continue
                # Split off 5' if x1 != x2
                if x1 < x2:
                    trimx = x2 - x1
                    gstr += '%s%d' % (t1,trimx)
                    x1 += trimx
                    gstr1[0] = '%s%d' % (t1,n1 - trimx)
                elif x2 < x1:
                    trimx = x1 - x2
                    gstr += '%s%d' % (t2,trimx)
                    x2 += trimx
                    gstr2[0] = '%s%d' % (t2,n2 - trimx)
                # Split off 3' if y1 != y2
                if y1 < y2:
                    trimy = y2 - y1
                    gstr2[0] = '%s%d' % (t2,n2 - trimy)
                    gstr2.insert(1,'%s%d' % (t2,trimy))
                elif y2 < y1:
                    trimy = y1 - y2
                    gstr1[0] = '%s%d' % (t1,n1 - trimy)
                    gstr1.insert(1,'%s%d' % (t1,trimy))
            if gstr1: gstr += rje.join(gstr1) + ':%s' % gdata1[2]
            else: gstr += rje.join(gstr2) + ':%s' % gdata2[2]
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return gstr
        except: self.errorLog('%s.combineGstring error' % self.prog()); raise
#########################################################################################################################
    def mapCS(self,lentry):  ### Generates a list of (qrypos,sbjpos,cs) tuples for a deconstructed cs string
        '''
        Generates a list of (qrypos,sbjpos,cs) tuples for a deconstructed cs string. Positions are the first base of any
        multibase elements.
        >> lentry:locdb entry = PAF file line loaded into locdb. Needs QryStart QryEnd SbjStart SbjEnd Strand cs
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Based on 'PAF to Local Alignment Conversion'
            #i# locfields = rje.split('Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd')
            #i# locfields += ['QryLen','SbjLen']
            #i# locdb.list['Fields'] = locfields+['Strand','cs']
            cspos = []  # List of tuples
            aln = lentry['cs']
            if not aln: return []
            if aln.startswith('Z:'): aln = aln[2:]
            #!# EndExtension has already been performed
            for x in ':-+*~?!_': aln = aln.replace(x,' %s' % x)
            aln = aln.split()
            qpos = lentry['QryStart']
            spos = lentry['SbjStart']; sadd = +1
            if lentry['Strand'] == '-': spos = lentry['SbjEnd']; sadd = -1
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cs in aln:
                #self.bugPrint(cs)
                if cspos and cs[0] in '+_': cspos.append((qpos,spos-1,cs))
                elif cspos and cs[0] in '-~': cspos.append((qpos-1,spos,cs))
                else: cspos.append((qpos,spos,cs))
                try:
                    if cs[0] == ':':    # Identity
                        ilen = rje.atoi(cs[1:])
                        qpos += ilen
                        spos += (ilen * sadd)
                    elif cs[0] == '*':    # Mismatch
                        ilen = 1
                        qpos += ilen
                        spos += (ilen * sadd)
                    elif cs[0] == '-':    # Insertion in subject
                        ilen = len(cs[1:])
                        spos += (ilen * sadd)
                    elif cs[0] == '+':    # Deletion in subject
                        ilen = len(cs[1:])
                        qpos += ilen
                    elif cs[0] == '~':    # Intron in subject
                        ilen = rje.atoi(cs[3:-2])
                        spos += (ilen * sadd)
                    elif cs[0] == '!':    # Missing in subject
                        ilen = rje.atoi(cs[1:])
                        qpos += ilen
                    elif cs[0] == '?':    # Uncertain match
                        ilen = rje.atoi(cs[1:])
                        lentry['Length'] += ilen
                        qpos += ilen
                        spos += (ilen * sadd)
                    else: raise ValueError('Unexpected cs element: %s' % cs)
                except:
                    self.bugPrint(qfull[:100])
                    self.bugPrint(sfull[:100])
                    self.bugPrint(cs)
                    self.debug(aln)
                    self.deBug(lentry)
                    self.errorLog('Unanticipated CS parsing error (query:%s, hit:%s)' % (lentry['Qry'],lentry['Hit']))
                    self.warnLog('CS parsing prematurely ended (query:%s, hit:%s)' % (lentry['Qry'],lentry['Hit']))
                    break
            cspos.append((qpos,spos,''))    #i# Final entry is next position in sequences.
            if qpos != lentry['QryEnd'] + 1:
                self.warnLog('"%s" should end at QryPos %d, not QryEnd %d' % (lentry['cs'],qpos-1,lentry['QryEnd']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],cspos))
            if lentry['Strand'] == '-' and spos != lentry['SbjStart'] - 1:
                self.warnLog('"%s" should end at SbjPos %d, not SbjStart %d' % (lentry['cs'],spos+1,lentry['SbjStart']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],cspos))
            elif lentry['Strand'] != '-' and spos != lentry['SbjEnd'] + 1:
                self.warnLog('"%s" should end at SbjPos %d, not SbjEnd %d' % (lentry['cs'],spos-1,lentry['SbjEnd']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],cspos))
            return cspos
        except: self.errorLog('%s.mapCS error' % self.prog()); return False
#########################################################################################################################
    def trimCS(self,lentry,newstart=0,newend=0,target='Qry'):  ### Trims lentry to revised start and end positions.
        '''
        Trims lentry to revised start and end positions.
        >> lentry:locdb entry = PAF file line loaded into locdb. Needs QryStart QryEnd SbjStart SbjEnd Strand cs
        >> newstart:int [0] = New start position - trim anything before this.
        >> newend:int [0] = New end position - trim anything after this.
        >> target:str ['Qry'] = Target of trimming (Qry or Hit/Sbj)
        << nentry:dict = Returns trimmed lentry. Also modified IN PLACE. Will return empty dictionary if totally removed.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if target == 'Hit': target = 'Sbj'
            if not newstart: newstart = lentry['%sStart' % target]
            if not newend: newend = lentry['%sEnd' % target]
            trim = False
            if newstart > lentry['%sEnd' % target]: return {}
            elif newstart > lentry['%sStart' % target]: trim = True
            if newend < lentry['%sStart' % target]: return {}
            elif newend < lentry['%sEnd' % target]: trim = True
            if not trim: return lentry
            nentry = rje.combineDict({},lentry)
            ### ~ [2] Trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Generate a list of (qrypos,sbjpos,cs) tuples
            #i# Remember that if Strand=='-' then Sbj positions are from the end, not the start.
            csmap = self.mapCS(lentry)
            (qpos,spos,null) = csmap[-1]
            newmap = []
            self.bugPrint('%s %s:%s-%s' % (lentry,target,newstart,newend))

            if qpos != lentry['QryEnd'] + 1:
                self.warnLog('"%s" should end at QryPos %d, not QryEnd %d' % (lentry['cs'],qpos-1,lentry['QryEnd']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],csmap))
            if lentry['Strand'] == '-' and spos != lentry['SbjStart'] - 1:
                self.warnLog('"%s" should end at SbjPos %d, not SbjStart %d' % (lentry['cs'],spos+1,lentry['SbjStart']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],csmap))
            elif lentry['Strand'] != '-' and spos != lentry['SbjEnd'] + 1:
                self.warnLog('"%s" should end at SbjPos %d, not SbjEnd %d' % (lentry['cs'],spos-1,lentry['SbjEnd']))
                self.bugPrint('\n%s' % lentry)
                self.debug('%s -> %s' % (lentry['cs'],csmap))

            csmap = csmap[:-1]

            #i# First, exclude all outside of trimmed region
            for (qpos,spos,cs) in csmap:
                (qlen,slen) = self.lenPartCS(cs)
                #self.bugPrint('%s:%d,%d' % (cs,qlen,slen))
                (qend,send) = (qpos,spos)
                #self.bugPrint('%s,%s %s' % (qpos,spos,cs))
                if target == 'Qry':
                    if qlen: qend = qpos + qlen - 1
                    if newend < qpos: break
                    elif cs[0] in '~_' and newend < qend: break    #i# do not trim partway into an intron
                    if newstart > qend: continue
                elif lentry['Strand'] == '-':
                    if slen: send = spos - slen + 1
                    if newend < send: continue
                    if newstart > send: break
                else:
                    if slen: send = spos + slen - 1
                    if newend < spos: break
                    elif cs[0] in '~_' and newend < send: break     #i# do not trim partway into an intron
                    if newstart > send: continue
                #if cs[0] in '~_':
                #    self.bugPrint('%s,%s %s -> +%s,%s -> %s,%s' % (qpos,spos,cs,qlen,slen,qend,send))
                #    continue #i# introns are always skipped: do not trim partway into an intron
                #self.debug('add')
                newmap.append((qpos,spos,cs))
            #if target == 'Qry': notrim = '~_-'
            #else:
            notrim = '~_+-'
            while newmap and newmap[0][2][0] in notrim: newmap.pop(0)   #i# do not trim partway into an intron/deletion
            while newmap and newmap[-1][2][0] in notrim: (qpos,spos,cs) = newmap.pop(-1)   #i# do not trim partway into an intron/deletion
            if not newmap: return {}
            prevcs = None
            nextcs = (qpos,spos,cs)
            starti = csmap.index(newmap[0])
            if starti: prevcs = csmap[starti-1]
            self.bugPrint('\n%s %s:%s-%s (%s)' % (self.db('local').entrySummary(lentry,collapse=True),target,newstart,newend,lentry['Strand']))
            self.bugPrint('|\nFull: %s' % csmap)
            self.bugPrint('|\nTrim: %s' % newmap)
            #i# Now check ends for trimming
            if target == 'Qry':
                #i# Look at start first
                if newmap[0][0] < newstart:  # Trim start
                    trimx = newstart - newmap[0][0]
                    qpos = newmap[0][0] + trimx
                    cs = newmap[0][2]
                    if cs[0] in ':?':    # Identity
                        ilen = rje.atoi(cs[1:])
                        ilen -= trimx
                        spos = newmap[0][1] + trimx
                        if lentry['Strand'] == '-': spos = newmap[0][1] - trimx
                        newmap[0] = (qpos,spos,'%s%d' % (cs[0],ilen))
                    elif cs[0] in '+!':    # Deletion in subject
                        newmap[0] = (qpos,newmap[0][1],cs[0] + cs[trimx+1:])
                    #i# Nothing else should change qpos
                    else: raise ValueError('Query pos/newstart mismatch for "%s"!' % cs)
                #i# Then check end
                cs = newmap[-1][2]
                endpos = newmap[-1][0]
                qlen = self.lenPartCS(cs)[0]
                if qlen: endpos += (qlen - 1)
                if newend < endpos:  # Trim end
                    trimx = endpos - newend
                    if cs[0] in ':?': cs = '%s%d' % (cs[0],rje.atoi(cs[1:]) - trimx)
                    elif cs[0] in '+!': cs = cs[:-trimx]
                    else: raise ValueError('Query pos/newend mismatch for "%s"!' % cs)
                    newmap[-1] = (newmap[-1][0],newmap[-1][1],cs)
            elif lentry['Strand'] == '-':
                #i# Look at start first
                if newmap[0][1] > newend:  # Trim start
                    trimx = newend - newmap[0][1]
                    spos = newmap[0][1] - trimx
                    cs = newmap[0][2]
                    if cs[0] in ':?':    # Identity
                        ilen = rje.atoi(cs[1:])
                        ilen -= trimx
                        qpos = newmap[0][0] + trimx
                        newmap[0] = (qpos,spos,'%s%d' % (cs[0],ilen))
                    elif cs[0] in '-~':    # Deletion/intron in subject
                        newmap[0] = (newmap[0][0],spos,cs[0] + cs[trimx+1:])
                    #i# Nothing else should change qpos
                    else: raise ValueError('Query pos/newstart mismatch for "%s"!' % cs)
                #i# Then check end
                cs = newmap[-1][2]
                endpos = newmap[-1][1]
                slen = self.lenPartCS(cs)[1]
                if slen: endpos -= (slen - 1)
                if newstart > endpos:  # Trim end
                    trimx = endpos - newstart
                    if cs[0] in ':?': cs = '%s%d' % (cs[0],rje.atoi(cs[1:]) - trimx)
                    elif cs[0] in '+!': cs = cs[:-trimx]
                    else:
                        self.debug('%s = %s -> endpos=%s' % (newmap[-1],slen,endpos))
                        raise ValueError('Query pos/newend mismatch for "%s"!' % cs)
                    newmap[-1] = (newmap[-1][0],newmap[-1][1],cs)
            else:
                #i# Look at start first
                if newmap[0][1] < newstart:  # Trim start
                    trimx = newstart - newmap[0][1]
                    spos = newmap[0][1] + trimx
                    cs = newmap[0][2]
                    if cs[0] in ':?':    # Identity
                        ilen = rje.atoi(cs[1:])
                        ilen -= trimx
                        qpos = newmap[0][0] + trimx
                        newmap[0] = (qpos,spos,'%s%d' % (cs[0],ilen))
                    elif cs[0] in '-~':    # Deletion/intron in subject
                        newmap[0] = (newmap[0][0],spos,cs[0] + cs[trimx+1:])
                    #i# Nothing else should change qpos
                    else: raise ValueError('Query pos/newstart mismatch for "%s"!' % cs)
                #i# Then check end
                cs = newmap[-1][2]
                endpos = newmap[-1][1]
                slen = self.lenPartCS(cs)[1]
                if slen: endpos += (slen - 1)
                if newend < endpos:  # Trim end
                    trimx = endpos - newend
                    if cs[0] in ':?': cs = '%s%d' % (cs[0],rje.atoi(cs[1:]) - trimx)
                    elif cs[0] in '+!': cs = cs[:-trimx]
                    else: raise ValueError('Query pos/newend mismatch for "%s"!' % cs)
                    newmap[-1] = (newmap[-1][0],newmap[-1][1],cs)
            self.bugPrint('|\nTrimmed: %s' % newmap)
            ### ~ [3] Regenerate CS string from newmap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newcs = []
            for (qpos,spos,cs) in newmap: newcs.append(cs)
            endstat = self.statsFromCS(newcs[-1])
            #if nentry['Qry'] == 'EOG092E076U':
            #    self.debug('\n%s' % nentry)
            #    self.debug('%s: %s' % (newcs[-1],endstat))
            newcs = rje.join(newcs,'')
            nentry['cs'] = newcs
            #i# Updated Start and End positions in nentry
            qstart = sstart = 0
            while newmap[qstart][2][0] == '-': qstart += 1
            while newmap[sstart][2][0] == '+': sstart += 1
            nentry['QryStart'] = newmap[qstart][0]
            nentry['QryEnd'] = newmap[-1][0] + max(0,endstat['QryCov']-1)
            if lentry['Strand'] == '-':
                nentry['SbjEnd'] = newmap[sstart][1]
                nentry['SbjStart'] = newmap[-1][1] - max(0,endstat['SbjCov']-1)
            else:
                nentry['SbjStart'] = newmap[sstart][1]
                nentry['SbjEnd'] = newmap[-1][1] + max(0,endstat['SbjCov']-1)
            self.bugPrint('\n%s\n|-- %s | %s | %s\n%s:%s-%s -> Qry:%s-%s | Sbj:%s-%s' % (self.db('local').entrySummary(nentry,collapse=True),prevcs,nentry['cs'],nextcs,target,newstart,newend,nentry['QryStart'],nentry['QryEnd'],nentry['SbjStart'],nentry['SbjEnd']))
            #i# Update general stats and check for consistency
            csdict = self.statsFromCS(newcs)
            nentry['Length'] = csdict['AlnLen']
            nentry['QryGS'] = '%d:%s:%d' % (nentry['QryStart'],csdict['QryGS'],nentry['QryEnd'])
            nentry['SbjGS'] = '%d:%s:%d' % (nentry['SbjStart'],csdict['SbjGS'],nentry['SbjEnd'])
            self.mapCS(nentry)

            #if nentry['Qry'] == 'EOG092E076U': self.debug('\n%s' % nentry)

            nentry['Identity'] = csdict[':']
            for field in nentry: lentry[field] = nentry[field]
            #self.debug('%s' % nentry)
            return nentry
        except: self.errorLog('%s.trimCS error' % self.prog()); return False
#########################################################################################################################
    def lenPartCS(self,cs): ### Returns (qry,sbj) lengths for CS element
        '''
        ### Returns (qry,sbj) lengths for CS element
        :param cs:str []
        :return:
        '''
        try:### ~ [1] Calculate Length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not cs: return (0,0)
            if len(cs) < 2:
                self.debug('"%s"' % aln)
                raise ValueError('cs element length < 2: "%s"' % cs)
            if cs[0] in ':?': ilen = rje.atoi(cs[1:])
            elif cs[0] in '~_':    # Intron in subject
                ilen = rje.atoi(cs[3:-2])
                #self.debug('%s: %d' % (cs,ilen))
            elif cs[0] == '*': ilen = 1   # Mismatch
            else: ilen = len(cs[1:])
            ### ~ [2] Return lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (qlen,slen) = (0,0)
            if cs[0] in ':*+_!': qlen = ilen
            if cs[0] in ':*-~': slen = ilen
            return (qlen,slen)
        except: self.errorLog('%s.lenPartCS error' % self.prog()); raise
#########################################################################################################################
    def statsFromCS(self,cs):  ### Generates a disctionary of stats from a cs string
        '''
        Generates a disctionary of stats from a cs string. Other stats can be extracted from these numbers.
        Note: the Gstrings produced do not have the flanking start and end positions:
        # start:[|+*-]\d+ list:end
        >> cs:str [] = CS string with ':' (identity), '*' (mismatch), '-' (sbj insertion), '+' (sbj deletion), '~' (sbj intron), '_' (qry intron), '!' (sbj missing), '?' (uncertain match)
        << dist: dictionary of elements, plus 'AlnLen', 'QryGS', 'SbjGS', 'QryCov', 'SbjCov'
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ncount = cs.count('n')
            csdict = {}
            aln = cs
            if aln.startswith('Z:'): aln = aln[2:]
            for x in ':-+*~?!_':
                csdict[x] = 0
                aln = aln.replace(x,' %s' % x)
            aln = aln.split()
            qrygs = []
            sbjgs = []
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cs in aln:
                if not cs: continue
                if len(cs) < 2:
                    self.debug('"%s"' % aln)
                    raise ValueError('cs element length < 2: "%s"' % cs)
                if cs[0] in ':?': ilen = rje.atoi(cs[1:])
                elif cs[0] in '~_':    # Intron in subject
                    ilen = rje.atoi(cs[3:-2])
                elif cs[0] == '*': ilen = 1   # Mismatch
                else: ilen = len(cs[1:])
                csdict[cs[0]] += ilen
                if cs[0] == ':':
                    if qrygs and qrygs[-1][0] == '|': qrygs[-1] = '|%d' % (rje.atoi(qrygs[-1][1:])+ilen)
                    else: qrygs.append('|%d' % ilen);
                    if sbjgs and sbjgs[-1][0] == '|': sbjgs[-1] = '|%d' % (rje.atoi(sbjgs[-1][1:])+ilen)
                    else: sbjgs.append('|%d' % ilen)
                elif cs[0] in '*':
                    if qrygs and qrygs[-1][0] == '*': qrygs[-1] = '*%d' % (rje.atoi(qrygs[-1][1:])+1)
                    else: qrygs.append('*%d' % ilen)
                    if sbjgs and sbjgs[-1][0] == '*': sbjgs[-1] = '*%d' % (rje.atoi(sbjgs[-1][1:])+1)
                    else: sbjgs.append('*%d' % ilen)
                elif cs[0] in '?': qrygs.append('*%d' % ilen); sbjgs.append('*%d' % ilen)
                elif cs[0] in '-~!': sbjgs.append('-%d' % ilen)
                elif cs[0] in '+_': qrygs.append('-%d' % ilen)
                else: raise ValueError('Unexpected cs element: %s' % cs)
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            csdict['AlnLen'] = sum(csdict.values()) - csdict['~'] - csdict['_'] - ncount  #i# Ignore introns and n's
            csdict['n'] = ncount
            csdict['QryGS'] = rje.join(qrygs,'')
            csdict['SbjGS'] = rje.join(sbjgs,'')
            #i# Ignore introns
            csdict['QryCov'] = csdict[':'] + csdict['*'] + csdict['+'] + csdict['!'] + csdict['?'] #x# + csdict['_']
            csdict['SbjCov'] = csdict[':'] + csdict['*'] + csdict['-'] + csdict['?'] #x# + csdict['~']
            return csdict
        except: self.errorLog('%s.statsFromCS error' % self.prog()); return False
#########################################################################################################################
    def revGstring(self,gstr):  ### Reverses Gstring, e.g. for Subject on -ve strand
        '''
        Reverses Gstring, e.g. for Subject on -ve strand
        :param gstr [str]: Gstring
        :return revgstr [str]: reversed GSstring [str]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# This is not actually required, as a CS string and Gstring is always 5' to 3'
            return gstr
        except: self.errorLog('%s.revGstring error' % self.prog()); return False
#########################################################################################################################
# Op	Regex	Description
# =	[ACGTN]+	Identical sequence (long form)
# :	[0-9]+	Identical sequence length
# *	[acgtn][acgtn]	Substitution: ref to query
# +	[acgtn]+	Insertion to the reference
# -	[acgtn]+	Deletion from the reference
# ~	[acgtn]{2}[0-9]+[acgtn]{2}	Intron length and splice signal
#########################################################################################################################
    def invertCS(self,lentry):  ### Inverts a CS string local entry so that the query is now the hit anf vice versa
        '''
        Inverts a CS string local entry so that the query is now the hit anf vice versa.
        >> lentry:locdb entry = PAF file line loaded into locdb. Needs QryStart QryEnd SbjStart SbjEnd Strand cs
        << nentry:locdb entry with Qry and Hit reversed.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Based on 'PAF to Local Alignment Conversion'
            #i# locfields = rje.split('Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd')
            #i# locfields += ['QryLen','SbjLen']
            #i# locdb.list['Fields'] = locfields+['Strand','cs']
            cspos = []  # List of tuples
            aln = lentry['cs']
            if not aln: return []
            if aln.startswith('Z:'): aln = aln[2:]
            #!# EndExtension has already been performed
            for x in ':-+*~?!_': aln = aln.replace(x,' %s' % x)
            aln = aln.split()
            nentry = rje.combineDict({},lentry)
            naln = []
            replace = {'-':'+','+':'-','!':'-','~':'_','_':'~'}
            #i# Reverse complement if switching strands
            revcomp = lentry['Strand'] == '-'
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cs in aln:
                # Replacements
                if cs[0] in replace:
                    if revcomp: cs = '%s%s' % (cs[0],rje_sequence.reverseComplement(cs[1:]))
                    naln.append('%s%s' % (replace[cs[0]],cs[1:]))
                elif cs[0] == '*':
                    if revcomp: cs = '%s%s' % (cs[0],rje_sequence.reverseComplement(cs[1:]))
                    naln.append('*%s%s' % (cs[2],cs[1]))
                else: naln.append(cs)
            ## ~ [2a] Reverse strand ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if lentry['Strand'] == '-':
                naln.reverse()
            nentry['cs'] = rje.join(naln,'')
            for x in ['Len','Start','End']:
                nentry['Qry%s' % x] = lentry['Sbj%s' % x]
                nentry['Sbj%s' % x] = lentry['Qry%s' % x]
            nentry['Qry'] = lentry['Hit']
            nentry['Hit'] = lentry['Qry']
            return nentry
        except: self.errorLog('%s.invertCS error' % self.prog()); return False
#########################################################################################################################
    ### <7> ### CS string variant calling functions                                                                     #
#########################################################################################################################
# The old samtools settings:
    # minqn=X         : Min. number of reads meeting qcut (QN) for output [10]
    # rid=T/F         : Whether to include Read ID (number) lists for each allele [True]
    # snponly=T/F     : Whether to restrict parsing output to SNP positions (will use mincut settings below) [False]
    # indels=T/F      : Whether to include indels in "SNP" parsing [True]
    # skiploci=LIST   : List of loci to exclude from pileup parsing (e.g. mitochondria) []
    # snptableout=T/F : Output filtered alleles to SNP Table [False]
    # ### ~ SNP Frequency Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    # mincut=X        : Minimum read count for minor allele (proportion if <1) [1]
    # absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    # biallelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    # ignoren=T/F     : Whether to exclude "N" calls from alleles [True]
    # ignoreref=T/F   : If False will always keep Reference allele and call fixed change as SNP [False]
#i# pafhead = Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality .. cs
# New PAFTools variant call settings
    # qcut=X          : Min. mapping quality score (0-255) to include a read [0]


# Want to generate/output:
    # snpdb = self.db('snp')  #Locus	Pos	Ref	N	QN	Seq	Dep	RID
    #chrI_YEAST__BK006935	26	A	148	148	A:107|C:2	0.74	A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
    #==> S288C-ISH.chrI.Q30.10.tdt <==

    #==> S288C-ISH.chrI.rid.tdt <==
    #RID	Locus	Start	End     Name
    #176 chrI_YEAST__BK006935	61	351


#########################################################################################################################
    def pafToVariants(self):    ### Performs variant calling on parsed PAF (CS) data
        '''
        Performs variant calling on parsed PAF (CS) data.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            #i# PAF data hase already been loaded into PAFDB
            #i# pafhead = Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality .. cs
            pafdb = self.db('paf')
            rdb = db.addEmptyTable('rid',['Qry','QryStart','QryEnd','SbjStart','SbjEnd','Hit','Strand','cs'],['Qry','QryStart','QryEnd','SbjStart','Hit'])
            sbjseqs = self.obj['Reference']
            sbjdict = {}
            if sbjseqs: sbjdict = sbjseqs.makeSeqNameDic('short')

            ### ~ [2] Convert PAF to Reads (RID) table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# The first step is to invert all of the PAF hits so that the Reference is the Query. This will mean that
            #i# when self.mapCS(lentry) tuples are created, they will be sorted according to the reference position,
            #i# which will enable them to be processed in parallel.
            #i# During this process, reads are filtered according to the QCut mapping quality.
            qcut = 0 # self.getInt('QCut')
            strand = '+-'   # self.getStr('Strand')   # Only keep reads on this strand
            #!# Check revcomp is right with this setting!
            pafdb.dropEntries(['Quality<%d' % qcut],inverse=False,log=True,logtxt='Mapping Quality < %d' % qcut,purelist=False,keylist=False)
            pafdb.dropEntriesDirect('Strand',strand,inverse=True,log=True,force=False)
            for entry in pafdb.entries():
                rentry = rdb.addEntry(self.invertCS(entry))
                #self.bugPrint(pafdb.entrySummary(entry,fields=['Qry','QryStart','QryEnd','SbjStart','SbjEnd','Hit','Strand'],collapse=True))
                #self.bugPrint(pafdb.entrySummary(entry,fields=['cs'],collapse=True))
                #self.mapCS(entry)
                #self.bugPrint('>>>>>>> Invert >>>>>>>>>')
                #self.bugPrint(rdb.entrySummary(rentry,fields=['Qry','QryStart','QryEnd','SbjStart','SbjEnd','Hit','Strand'],collapse=True))
                #self.bugPrint(rdb.entrySummary(rentry,fields=['cs'],collapse=True))
                #self.mapCS(rentry)

            #i# Reads are given a RID based on revised QryStart
            rdb.newKey(['Qry','QryStart','QryEnd','SbjStart','Hit'])    # This is just to get the right sorting
            rdb.addField('RID'); rid = 0
            for entry in rdb.entries(sorted=True): rid += 1; entry['RID'] = rid
            rdb.newKey(['RID'],startfields=True)

            #i# Cut down and rename fields
            #RID Locus Start End Name Strand cs
            rdb.renameField('Qry','Locus')
            rdb.renameField('QryStart','Start')
            rdb.renameField('QryEnd','End')
            rdb.renameField('Hit','Name')

            #i# Save cutdown fields to RID table
            #RID	Locus	Start	End     Name
            rdb.saveTable(savefields=['RID','Locus','Start','End','Name','Strand','cs'])

            ### ~ [3] ~ Generate RID SNP table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mindepth = self.getInt('MinQN')   # minqn from samtools
            indels = self.getBool('Indels') #Whether to include indels in "SNP" parsing [True]
            skiploci = self.list['SkipLoci']  # List of loci to exclude from pileup parsing (e.g. mitochondria) []
            snptableout = self.getBool('SNPTableOut') # Output filtered alleles to SNP Table [False]
            mincut = self.getNum('MinCut')        # Minimum read count for minor allele (proportion if <1) [1]
            absmincut = self.getInt('AbsMinCut')     # Absolute minimum read count for minor allele (used if mincut<1) [2]
            biallelic = self.getBool('BiAllelic')   # Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
            ignoren = self.getBool('IgnoreN')     # Whether to exclude "N" calls from alleles [True]
            ignoreref = self.getBool('IgnoreRef')   # If False will always keep Reference allele and call fixed change as SNP [False]
            snpdb = db.addEmptyTable('snp',['Locus','Pos','Ref','N','QN','Seq','Dep','VarFreq','RID'],['Locus','Pos'])
            # Locus	Pos	Ref	N	QN	Seq	Dep	RID
            #chrI_YEAST__BK006935	26	A	148	148	A:107|C:2	0.74	A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            varlist = []    # List of (qrypos,rid,cs) built from read mapCS lists of (qrypos,sbjpos,cs) tuples for a deconstructed cs string
            for locus in rdb.indexKeys('Locus'):
                if locus in skiploci: continue
                refseq = ''
                if locus in sbjdict: refseq = sbjseqs.seqSequence(sbjdict[locus])
                ## ~ [3a] ~ Build RID CS Mapping list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rx = 0.0; rtot = rdb.entryNum()
                for rentry in rdb.indexEntries('Locus',locus):
                    self.progLog('\r#MAPCS','Parsing CS Mapping for %s: %.1f%%' % (locus,rx/rtot)); rx += 100.0
                    rid = rentry['RID']
                    rentry['QryStart'] = rentry['Start']
                    rentry['QryEnd'] = rentry['End']
                    for (qrypos,sbjpos,cs) in self.mapCS(rentry):
                        if cs: varlist.append((qrypos,cs,rid))
                    varlist.append((rentry['End']+1,'END',rid))
                self.progLog('\r#MAPCS','Parsing CS Mapping for %s: sorting...' % locus)
                varlist.sort()
                baddep = []
                ## ~ [3b] ~ Parse out SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Need to make sure that the matching RIDs are stored
                matchrids = []  # List of currently matching RIDs
                delrids = []    # List of current deletion RIDs
                rx = 0.0; rtot = len(varlist); vx = 0
                while varlist:
                    pos = varlist[0][0]
                    varn = 0    # Number of variant reads
                    ## Generate dictionary of c
                    posvar = {'N':[]}     # cs:[RIDs]
                    while varlist and varlist[0][0] == pos:
                        (pos,cs,rid) = varlist.pop(0)
                        self.progLog('\r#CSVAR','Extracting variants from CS Mapping for %s: %.2f%%' % (locus,rx/rtot)); rx += 100.0
                        #self.bugPrint('%s %s "%s"\n|-- :%s\n|-- -%s' % (pos,rid,cs,matchrids,delrids))
                        # Update deletion and match lists
                        if cs[0] != '-' and rid in delrids: delrids.remove(rid)
                        if cs[0] not in ':=' and rid in matchrids: matchrids.remove(rid)
                        if cs == 'END': continue
                        if cs[0] == '-':  # Deletion from the reference
                            if rid not in delrids: delrids.append(rid)
                            continue
                        if cs[0] in ':=':
                            if rid not in matchrids: matchrids.append(rid)
                            continue
                        # Introns are absent of read coverage, so skip
                        if cs[0] == '~':  # Intron length and splice signal
                            continue
                        # Update posvar dictionary with cs element
                        if cs[0] == '+' and not indels:  # Insertion to the reference
                            cs = 'N'
                        if cs not in posvar: posvar[cs] = []
                        if cs[0] == '*':  # Substitution: ref to query
                            pass
                        posvar[cs].append(rid)
                        varn += 1
                    if not varlist: break
                    # Update deletions and matches
                    if delrids and indels: posvar['-'] = delrids[0:]
                    elif delrids: posvar['N'] += delrids[0:]
                    if matchrids: posvar[':'] = matchrids[0:]
                    elif not ignoreref: posvar[':'] = []
                    varn += len(delrids)
                    varn += len(matchrids)
                    refdep = len(matchrids)
                    #i# Filter variants here
                    if varn < mindepth:
                        continue
                    # Drop low abundance and ambiguous variants
                    vardep = int(varn * mincut + 0.5)        # Minimum read count for minor allele (proportion if <1) [1]
                    vardep = max(vardep,absmincut)
                    qn = varn
                    if ignoren and 'N' in posvar:
                        qn = varn - len(posvar.pop('N'))
                    ref = '.'
                    if refseq: ref = refseq[pos-1].upper()
                    for cs in posvar.keys():
                        if cs == ':' and not ignoreref: continue
                        elif len(posvar[cs]) < vardep: qn -= len(posvar.pop(cs))
                        elif cs[0] == '*': ref = cs[2].upper()
                    if len(posvar) < 2: continue   #i# No variants escape filtering
                    # Apply biallelic filter
                    if biallelic and len(posvar) != 2: continue
                    # Add SNP to table
                    seq = []
                    rid = []
                    #if refseq: self.bugPrint('%s: %s' % (refseq[max(0,pos-10):pos-1].lower()+refseq[pos-1].upper()+refseq[pos:pos+10].lower(),posvar))
                    #if refseq: self.bugPrint(refseq[max(0,pos-10):pos-1].lower()+refseq[pos-1].upper()+refseq[pos:pos+10].lower())
                    #else: self.bugPrint('?: %s' % (posvar))
                    for cs in posvar:
                        posvar[cs].sort()
                        if cs == ':': alt = ref.upper()
                        elif cs == '-': alt = '-'
                        elif cs[0] == '*': alt = cs[1].upper()
                        elif cs[0] == '+':
                            alt = cs.upper()
                            #if ref == '.': alt = cs.upper()
                            #else: alt = ref + cs[1:].upper()
                        else: alt = cs[1:].upper()
                        seq.append('%s:%s' % (alt,len(posvar[cs])))
                        rid.append('%s:%s' % (alt,str(posvar[cs])[1:-1].replace(' ','')))
                    seq.sort(); seq = rje.join(seq,'|')
                    rid.sort(); rid = rje.join(rid,'|')
                    ventry = {'Locus':locus,'Pos':pos,'Ref':ref,'N':varn,'QN':qn,'Seq':seq,'Dep':varn,'RID':rid,'VarFreq':float(qn-refdep)/qn}
                    ventry = snpdb.addEntry(ventry); vx += 1
                    #self.debug(snpdb.entrySummary(ventry))
                self.printLog('\r#CSVAR','Extracted %s variants from CS Mapping for %s.   ' % (rje.iStr(vx),locus))


# Op	Regex	Description
# =	[ACGTN]+	Identical sequence (long form)
# :	[0-9]+	Identical sequence length
# *	[acgtn][acgtn]	Substitution: ref to query
# +	[acgtn]+	Insertion to the reference
# -	[acgtn]+	Deletion from the reference
# ~	[acgtn]{2}[0-9]+[acgtn]{2}	Intron length and splice signal


            ### ~ [4] ~ Save SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Include headers containing key settings to check when re-loading data
            if snptableout: snpdb.saveToFile(sfdict={'VarFreq':3})


        except: self.errorLog('%s.pafToVariants error' % self.prog()); return False
#########################################################################################################################
    ### <8> ### Long read mapping to PAF/BAM files                                                                      #
#########################################################################################################################
    def longreadMinimapPAF(self,paffile=None):  ### Performs long read versus assembly minimap2 and saves to PAF
        '''
        Performs long read versus assembly minimap2 and saves to PAF
        :return: paffile
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            paf = self
            if not paffile:
                if self.getStrLC('PAFIn') and not self.getStrLC('PAFIn') in ['minimap','minimap2']:
                    paffile = self.getStr('PAFIn')
                else:
                    paffile = self.baseFile() + '.paf'

            ### ~ [2] ~ Generate individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.force() or not rje.exists(paffile):
                if not self.list['Reads']: raise IOError('Cannot generate long read PAF file without reads=LIST')
                mmv = os.popen('%s --version' % paf.getStr('Minimap2')).read()
                if not mmv: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % paf.getStr('Minimap2'))
                rje.backup(self,paffile)
                #!# Check these BAM files have headers! #!#
                paflist = []; rx = 0
                if not self.list['ReadType']:
                    self.warnLog('Read Type not given (pb/ont/hifi): check readtype=LIST. Will use "ont".')
                elif len(self.list['ReadType']) == 1 and len(self.list['Reads']) != 1:
                    self.printLog('#READS','Using "%s" as read type for all long reads' % self.list['ReadType'][0])
                elif len(self.list['ReadType']) != len(self.list['Reads']):
                    self.warnLog('reads=FILELIST vs readtype=LIST length mismatch: check readtype=LIST. Will cycle if needed.')
                for readfile in self.list['Reads']:
                    if not rje.exists(readfile): raise IOError('Read file "{0}" not found!'.format(readfile))
                    if self.list['ReadType']:
                        try: rtype = self.list['ReadType'][rx]; rx +=1
                        except: rtype = self.list['ReadType'][0]; rx = 1
                        if rtype in ['pacbio','pac']: rtype = 'pb'
                        if rtype in ['hifi','ccs']: rtype = 'hifi'
                        if rtype not in ['ont','pb','hifi']:
                            self.warnLog('Read Type "%s" not recognised (pb/ont/hifi): check readtype=LIST. Will use "ont".' % rtype)
                            rtype = 'ont'
                    else: rtype = 'ont'
                    prefix = '{0}.{1}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(readfile,strip_path=True))
                    maplog = '{0}.log'.format(prefix)
                    ## ~ [2a] Make SAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    maprun = '{0} -t {1} --secondary=no -o {2}.paf -L -x map-{3} {4} {5}'.format(paf.getStr('Minimap2'),self.threads(),prefix,rtype,self.getStr('SeqIn'),readfile)
                    #i# Minimap2 now has a HiFi mapping mode
                    #if rtype in ['hifi']:
                    #    maprun = '{0} -t {1} --secondary=no -o {2}.paf -L -x asm20 {4} {5}'.format(paf.getStr('Minimap2'),self.threads(),prefix,rtype,self.getStr('SeqIn'),readfile)
                    logline = self.loggedSysCall(maprun,maplog,append=False)
                    #!# Add check that run has finished #!#
                    rpfile = '{0}.paf'.format(prefix)
                    if not rje.exists(rpfile): raise IOError('Minimap2 output not found: {0}'.format(rpfile))
                    paflist.append(rpfile)

            ### ~ [3] ~ Merge individual PAF files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if len(paflist) > 1:
                    # samtools merge - merges multiple sorted input files into a single output.
                    bammerge = 'cat {0} > {1}'.format(' '.join(paflist),paffile)
                    logline = self.loggedSysCall(bammerge,append=True)
                    if not rje.exists(paffile): raise IOError('Merged PAF file "%s" not generated' % paffile)
                    for sortbam in paflist: os.unlink(sortbam)
                else: os.rename(paflist[0],paffile)

            return paffile
        except:
            self.errorLog('PAF.longreadMinimapPAF() error'); raise
            #return None
#########################################################################################################################
    def longreadMPileup(self): return self.longreadMinimapBAM()
    def longreadMinimapBAM(self):  ### Performs long read versus assembly minimap2 and to converts BAM
        '''
        Performs long read versus assembly minimap2 and converts to BAM file
        :return: bamfile/None
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            paf = self
            bamfile = self.baseFile() + '.bam'

            ### ~ [2] ~ Generate individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.force() or not rje.exists(bamfile):
                mmv = os.popen('%s --version' % paf.getStr('Minimap2')).read()
                if not mmv: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % paf.getStr('Minimap2'))
                rje.backup(self,bamfile)
                #!# Check these BAM files have headers! #!#
                bamlist = []; rx = 0
                if not self.list['ReadType']:
                    self.warnLog('Read Type not given (pb/ont): check readtype=LIST. Will use "ont".')
                elif len(self.list['ReadType']) == 1 and len(self.list['Reads']) != 1:
                    self.printLog('#READS','Using "%s" as read type for all long reads' % self.list['ReadType'][0])
                elif len(self.list['ReadType']) != len(self.list['Reads']):
                    self.warnLog('reads=FILELIST vs readtype=LIST length mismatch: check readtype=LIST. Will cycle if needed.')
                for readfile in self.list['Reads']:
                    if not rje.exists(readfile): raise IOError('Read file "{0}" not found!'.format(readfile))
                    if self.list['ReadType']:
                        try: rtype = self.list['ReadType'][rx]; rx +=1
                        except: rtype = self.list['ReadType'][0]; rx = 1
                        if rtype in ['pacbio','pac']: rtype = 'pb'
                        if rtype in ['hifi','ccs']: rtype = 'hifi'
                        if rtype not in ['ont','pb','hifi']:
                            self.warnLog('Read Type "%s" not recognised (pb/ont): check readtype=LIST. Will use "ont".' % rtype)
                            rtype = 'ont'
                    else: rtype = 'ont'
                    prefix = '{0}.{1}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(readfile,strip_path=True))
                    maplog = '{0}.log'.format(prefix)
                    ## ~ [2a] Make SAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    maprun = '{0} -t {1} --secondary=no -o {2}.sam -L -ax map-{3} {4} {5}'.format(paf.getStr('Minimap2'),self.threads(),prefix,rtype,self.getStr('SeqIn'),readfile)
                    if rtype in ['hifi']:
                        maprun = '{0} -t {1} --secondary=no -o {2}.sam -L -ax asm20 {4} {5}'.format(paf.getStr('Minimap2'),self.threads(),prefix,rtype,self.getStr('SeqIn'),readfile)
                    logline = self.loggedSysCall(maprun,maplog,append=False)
                    #!# Add check that run has finished #!#
                    ## ~ [2b] Converting SAM to BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #sam2bam = 'samtools view -bo {0}.tmp.bam -@ {1} -S {2}.sam'.format(prefix,self.threads()-1,prefix)
                    self.printLog('#BAM','Converting SAM to BAM. Using a single thread due to past issues of missing data.')
                    sam2bam = 'samtools view -bo {0}.tmp.bam -S {1}.sam'.format(prefix,prefix)
                    logline = self.loggedSysCall(sam2bam,maplog,append=True,nologline='No stdout from sam2bam')
                    #!# Add check that run has finished #!#
                    ## ~ [2c] Sorting BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.printLog('#BAM','Sorting BAM file.')
                    sortbam = '{0}.bam'.format(prefix)
                    bamsort = 'samtools sort -@ {0} -o {1}.bam -m 6G {2}.tmp.bam'.format(self.threads()-1,prefix,prefix)
                    logline = self.loggedSysCall(bamsort,maplog,append=True)
                    #!# Add check that run has finished #!#
                    if not rje.exists(sortbam): raise IOError('Sorted BAM file "%s" not generated' % sortbam)
                    os.unlink('{0}.sam'.format(prefix))
                    os.unlink('{0}.tmp.bam'.format(prefix))
                    bamlist.append(sortbam)

            ### ~ [3] ~ Merge individual BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if len(bamlist) > 1:
                    # samtools merge - merges multiple sorted input files into a single output.
                    bammerge = 'samtools merge -@ {0} {1} {2}'.format(self.threads()-1,bamfile,' '.join(bamlist))
                    logline = self.loggedSysCall(bammerge,append=True)
                    if not rje.exists(bamfile): raise IOError('Merged BAM file "%s" not generated' % bamfile)
                    for sortbam in bamlist: os.unlink(sortbam)
                else: os.rename(bamlist[0],bamfile)

            ## ~ [3a] Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            baifile = '{0}.bai'.format(bamfile)
            rje.checkForFiles(filelist=[bamfile,baifile],basename='',log=self.log,cutshort=False,ioerror=False)
            if self.needToRemake(baifile,bamfile):
                makebai = 'samtools index -b {0} {1}.bai'.format(bamfile,bamfile)
                logline = self.loggedSysCall(makebai,append=True)
                #os.system('samtools index -b {0} {1}.bai'.format(bamfile,bamfile))

            return bamfile
        except:
            self.errorLog('PAF.longreadMinimapBAM() error')
            return None
#########################################################################################################################
    ### <9> ### Read coverage position checking                                                                         #
#########################################################################################################################
    def checkPos(self,save=True,cdb=None):  ### Checks read coverage spanning given positions
        '''
        Checks read coverage spanning given positions.
        >> save:bool [True] = Whether to save table
        >> cdb:Table [None] = Existing table with self.list['CheckFields'] fields for checking.
        << returns table of checked positions
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forker = self.obj['Forker']
            ## ~ [1a] ~ Setup basefile and PAF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('Basefile'): self.baseFile(rje.baseFile(self.getStr('PAFIn'),strip_path=True))
            if not self.getStrLC('Basefile'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
            basefile = self.baseFile()
            self.printLog('#BASE',self.baseFile())
            db = self.db()
            if not db:
                db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T','basefile=%s' % self.baseFile()])
            if self.getStrLC('PAFIn') in ['minimap','minimap2']:
                self.setStr({'PAFIn':'%s.paf' % self.baseFile()})
                self.printLog('#PAFIN','Minimap2 PAF file set: %s' % self.getStr('PAFIn'))
            ## ~ [1b] ~ Setup coverage check table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            covflanks = self.list['CheckFlanks']
            if not len(self.list['CheckFields']) == 3:
                raise ValueError('checkfields=LIST must have exactly 3 elements: Locus, Start, End. %d found!' % len(self.list['CheckFields']))
            [locusfield,startfield,endfield] = self.list['CheckFields']
            #cdb = db.addTable(self.getStr('CheckPos'),mainkeys=self.list['CheckFields'],name='check',expect=True)
            if not cdb:
                cdb = db.addTable(self.getStr('CheckPos'),mainkeys='auto',name='check',expect=True)
            cdb.dataFormat({startfield:'int',endfield:'int'})
            cdb.compress(self.list['CheckFields'],rules={self.getStr('SpanID'):'list'})
            cdb.setStr({'Delimit':'\t'})
            cdb.addField('MaxFlank5',evalue=-1)
            cdb.addField('MaxFlank3',evalue=-1)
            for covx in covflanks: cdb.addField('Span%d' % covx,evalue=0)
            ## ~ [1c] Temp directory for forked depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)
            ## ~ [1d] Output directory for SpanID reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            iddir = '%s_spanid/' % basefile
            if self.getStrLC('SpanID'):
                if not rje.exists(iddir): rje.mkDir(self,iddir)


            ### ~ [2] ~ Check/Generate PAF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            paffile = self.longreadMinimapPAF()
            if not paffile: raise IOError('PAF file "{0}" not found'.format(paffile))

            ### ~ [3] ~ Fork out checking of positions using PAF file and awk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for centry in cdb.entries():
                locus = centry[locusfield]
                cstart = centry[startfield]
                cend = centry[endfield]
                tmpfile = '{0}{1}.{2}.{3}.{4}.tmp'.format(tmpdir,basefile,locus,cstart,cend)
                if rje.exists(tmpfile):
                    if not self.force() and len(open(tmpfile,'r').readline().split()) > 1: skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                #i# NOTE: PAF positions are 0-based, so need to subtract 1
                forker.list['ToFork'].append("awk -F '\\t' '$6==\"{0}\" && $8<={1} && $9>={2}' {3} > {4} && echo 'awk complete' >> {5}".format(locus,int(cstart)-1,int(cend)-1,paffile,tmpfile,tmpfile))
                self.bugPrint(forker.list['ToFork'][-1])
            self.printLog('#CHECK','{0} Coverage check regions queued for forking ({1} existing files deleted); {2} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [3b] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 1:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of coverage analysis completed.')
                else:
                    try:
                        self.errorLog('Coverage forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Coverage forking did not complete')

            ### ~ [4] Read in and process coverage calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # pafhead = rje.split('Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality')
            #!# Need to add read lists for each contaminant #!#
            ## ~ [4a] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cx = 0.0; ctot = cdb.entryNum(); spanid = {}
            for centry in cdb.entries():
                self.progLog('\r#CHECK','Processing checkpos coverage: %.2f%%' % (cx/ctot)); cx += 100.0
                locus = centry[locusfield]
                cstart = centry[startfield]
                cend = centry[endfield]
                centry['MaxFlank5'] = cstart - 1
                tmpfile = '{0}{1}.{2}.{3}.{4}.tmp'.format(tmpdir,basefile,locus,cstart,cend)
                try:
                    if not rje.exists(tmpfile):
                        raise IOError('Cannot find {0}'.format(tmpfile))
                    self.quiet()
                    pafdb = self.parsePAF(tmpfile,'tmp')
                    self.talk()
                    for pentry in pafdb.entries():
                        if pentry['Hit'] != locus: raise ValueError('Hit sequence mismatch for "%s"' % tmpfile)
                        centry['MaxFlank3'] = pentry['SbjLen'] - centry[endfield]
                        for covx in covflanks:
                            if pentry['SbjStart'] <= max(1,centry[startfield]-covx) and pentry['SbjEnd'] >= min(pentry['SbjLen'],centry[endfield]+covx):
                                centry['Span%d' % covx] += 1
                    if self.getStrLC('SpanID'):
                        idlist = pafdb.dataList(pafdb.entries(),'Qry',sortunique=False,empties=False)
                        for spanner in centry[self.getStr('SpanID')].split('|'):
                            if spanner not in spanid: spanid[spanner] = []
                            spanid[spanner] += idlist
                    db.deleteTable(pafdb)
                except:
                    self.talk()
                    self.errorLog('Coverage result processing error',quitchoice=True)
                    continue
            self.printLog('\r#CHECK','Processing checkpos coverage complete!')
            ## ~ [4b] Save and tidy data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if save: cdb.saveToFile()
            if not self.debugging() and not self.dev():
                tx = 0
                for centry in cdb.entries():
                    locus = centry[locusfield]
                    cstart = centry[startfield]
                    cend = centry[endfield]
                    tmpfile = '{0}{1}.{2}.{3}.{4}.tmp'.format(tmpdir,basefile,locus,cstart,cend)
                    os.unlink(tmpfile); tx += 1
                self.printLog('#TMP','%s temp files deleted' % rje.iStr(tx))
            ## ~ [4c] Save spanning read IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('SpanID'):
                nonx = 0; idx = 0
                for spanner in rje.sortKeys(spanid):
                    if not spanid[spanner]:
                        nonx += 1
                        if self.v() > 0: self.printLog('#SPANID','No %s spanning read IDs to output.' % (spanner))
                        continue
                    spout = '%s%s.%s.span.id' % (iddir,basefile,spanner)
                    spids = rje.sortUnique(spanid[spanner])
                    open(spout,'w').write('\n'.join(spids+['']))
                    self.printLog('#SPANID','%s %s spanning read IDs output to %s' % (rje.iLen(spids),spanner,spout))
                    idx += 1
                self.printLog('#SPANID','%s regions with spanning read IDs output to %s.' % (rje.iStr(idx),iddir))
                self.printLog('#SPANID','%s regions without spanning read IDs to output.' % (rje.iStr(nonx)))
            return cdb

        except: self.errorLog('Problem during %s checkPos().' % self.prog()); return False  # Setup failed
#########################################################################################################################
### End of SECTION II: PAF Class                                                                                        #
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
    try: PAF(mainlog,cmd_list).run()

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
