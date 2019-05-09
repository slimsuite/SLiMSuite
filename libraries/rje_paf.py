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
Version:      0.5.0
Last Edit:    02/05/19
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to use Minimap alignments in PAF format (with the --cs flag) to replace blastn for all-by-all
    genome assembly comparisons for PAGSAT etc.

    Input is a PAF file with --cs flag, produced by minimap:

    `minimap2 --cs -N 100 -p 0.0001 -x asm5 $REFERENCE $ASSEMBLY > $PAFILE`

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
    localaln=T/F    : Whether to output local alignments in Local Table [False]
    mockblast=T/F   : Whether to output mock BLAST headers even when not appropriate [True]
    ### ~ Minimap2 run/mapping options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [0]
    endextend=X     : Extend minimap2 hits to end of sequence if query region with X bp of end [10]
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    mapopt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm5]
    mapsplice=T/F   : Switch default minimap2 options to `-x splice -uf -C5` [False]
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
import rje, rje_obj, rje_db, rje_seqlist, rje_sequence
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
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_PAF', '0.5.0', 'May 2019', '2019')
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
#pafhead = string.split('Qry QryLen QryStart QryEnd Strand Hit HitLen HitStart HitEnd Identity Length Quality')
pafhead = string.split('Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality')
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

pafaln = string.split('tp cm s1 s2 NM MD AS ms nn ts cg cs dv')
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
    - Minimap2=PROG   : Full path to run minimap2 [minimap2]
    - PAFIn=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping []
    - Reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
    - SeqIn=FASFILE   : Input genome to identify variants in ($ASSEMBLY) []
    - TmpDir=PATH     : Temporary directory to use when forking [./tmp/]

    Bool:boolean
    - LocalAln = Whether to keep local alignments in Local Table [False]
    - MapSplice=T/F   : Switch default minimap2 options to `-x splice -uf -C5` [False]
    - MockBLAST = Whether to output mock BLAST headers even when not appropriate [True]
    - UniqueHit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [False]
    - UniqueOut=T/F   : Whether to output *.qryunique.tdt and *.hitunique.tdt tables of unique coverage [True]

    Int:integer
    - EndExtend=X         : Extend minimap2 hits to end of sequence if with X bp [10]
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [0]

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    
    MapOpt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm5]

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
        self.strlist = ['Minimap2','PAFIn','SeqIn','Reference','TmpDir']
        self.boollist = ['LocalAln','MapSplice','MockBLAST','UniqueHit','UniqueOut']
        self.intlist = ['EndExtend','MinLocLen']
        self.numlist = ['MinLocID']
        self.filelist = []
        self.listlist = []
        self.dictlist = ['MapOpt']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'Minimap2':'minimap2','TmpDir':'./tmp/'})
        self.setBool({'LocalAln':False,'MapSplice':False,'MockBLAST':True,'UniqueHit':False,'UniqueOut':True})
        self.setInt({'EndExtend':10,'MinLocLen':1})
        self.setNum({'MinLocID':0.0})
        self.dict['MapOpt'] = {} #'N':'100','p':'0.0001','x':'asm5'}
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
                self._cmdRead(cmd,type='file',att='SeqIn',arg='assembly')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='file',att='Reference',arg='searchdb')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='file',att='Reference',arg='refgenome')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['Minimap2'])   # Normal strings
                self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['PAFIn','SeqIn','Reference'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['LocalAln','MapSplice','MockBLAST','UniqueHit','UniqueOut'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['EndExtend','MinLocLen'])   # Integers
                self._cmdReadList(cmd,'perc',['MinLocID'])   # 0-100 percentage, converted x100 if <=1
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['MapOpt']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.getStrLC('PAFIn') in ['minimap','minimap2']:
                self.setStr({'PAFIn':'%s.paf' % self.baseFile()})
                self.printLog('#PAFIN','Minimap2 PAF file set: %s' % self.getStr('PAFIn'))
                self.minimap2()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.parsePAF(): return False
            if not self.pafToLocal(minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen')): raise ValueError
            self.localToUnique(save=self.getBool('UniqueOut'))   #!# Sort out output of unique tables!
            self.uniqueToHitSum()
            self.localToGABLAM()
            #? Optional deletion of empty self.getStr('TmpDir') ?#
            return self.db()
        except SystemExit: raise    # Fork child
        except:
            self.errorLog(self.zen())
            return None  # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Basefile'): self.baseFile(rje.baseFile(self.getStr('PAFIn')))
            if not self.getStrLC('Basefile'): self.baseFile('rje_paf')
            self.printLog('#BASE',self.baseFile())
            if self.getBool('MapSplice'): defaults = {'x':'splice'} #,'uf':'','C5':''}
            else: defaults = {'N':'100','p':'0.0001','x':'asm5'}
            self.dict['MapOpt'] = rje.combineDict(defaults,self.dict['MapOpt'],overwrite=True)
            self.devLog('#MAPOPT','%s' % self.dict['MapOpt'])
            if self.getInt('Forks') > 0:
                if 't' not in self.dict['MapOpt']: self.dict['MapOpt']['t'] = self.getInt('Forks')
            elif self.i() >=0 and rje.yesNo('Forking recommended to speed up operation? Set forks>0?'):
                self.setInt({'Forks':max(0,rje.getInt('Set number of forks:',confirm=True))})
                self.printLog('#FORKS','Forks set to forks=%d.' % self.getInt('Forks'))
            else: self.printLog('#FORKS','Forking recommended to speed up operation: set forks>0 to use multiple threads.')
            ## ~ [1a] Database Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T','basefile=%s' % self.baseFile()])
            ## ~ [1b] Reference sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('Reference'):
                raise ValueError('Reference missing: %s' % rje.argString(self.cmd_list))
            self.printLog('#REFIN','Loading reference (target) sequences.')
            self.obj['Reference'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Reference'),'autoload=T','seqmode=file'])
            ## ~ [1c] Reference sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('SeqIn'):
                raise ValueError('SeqIn missing: %s' % rje.argString(self.cmd_list))
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
    def parsePAF(self,pafin=None,table='paf'):  ### Parse PAF into Database table.
        '''
        Parse PAF into Database table.
        >> pafin:file [None] = PAF file to load, or use
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not pafin:
                if not self.getStrLC('PAFIn'): raise ValueError('No PAFIn given!')
                pafin = self.getStr('PAFIn')
            if not rje.exists(pafin): raise IOError('PAFIn file "%s" not found!' % pafin)
            self.setStr({'PAFIn':pafin})
            #self.printLog('#~~#','## ~~~~~ Parsing PAF Alignments ~~~~~ ##')
            self.headLog('Parsing PAF Alignments')
            self.printLog('\r#PAF','Parsing %s' % pafin)

            ### ~ [2] Load PAF file with auto-counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headx = len(pafhead)
            pafdb = db.addEmptyTable(table,['#']+pafhead+pafaln,['#'])
            px = 0
            with open(pafin,'r') as PAF:
                for line in PAF:
                    self.progLog('\r#PAF','Parsing %s lines' % (rje.iStr(px)))
                    px += 1
                    data = string.split(rje.chomp(line),'\t')
                    # Generate entry
                    pentry = {'#':px}
                    for i in range(headx):
                        field =  pafhead[i]
                        if field in ['Qry','Hit','Strand']: pentry[field] = data[i]
                        else: pentry[field] = string.atoi(data[i])
                    for pdat in data[headx:]:
                        pentry[pdat[:2]] = pdat[3:]
                    # Reformat data to 1->L positions
                    for field in string.split('QryStart SbjStart'): pentry[field] += 1
                    pafdb.addEntry(pentry)
            self.printLog('\r#PAF','Parsed %s lines from %s' % (rje.iStr(px),pafin))
            #if not pafdb.entryNum(): self.warnLog('No hits parsed from %s: check minimap2=PROG setting' % pafin)

            ### ~ [3] Return PAF Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dev(): pafdb.saveToFile()
            return pafdb
        except: self.errorLog('%s.parsePAF error' % self.prog()); return False
#########################################################################################################################
    def pafToLocal(self,pafdb=None,save=True,alnseq=True,minlocid=0,minloclen=1):  ### Parse PAF table into Local alignment table.
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
            locfields = string.split('Qry Hit AlnNum BitScore Expect Length Identity Positives QryStart QryEnd SbjStart SbjEnd')
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
                # Process
                for x in ':-+*~?!': aln = aln.replace(x,' %s' % x)
                #self.bugPrint('|===|')
                #self.bugPrint(lentry)
                aln = string.split(aln)
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
                            ilen = string.atoi(cs[1:])
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += '|' * ilen
                        elif cs[0] == '*':    # Mismatch
                            ilen = 1
                            if cs[2] != qfull[0].lower(): qwarnx += 1
                            if cs[1] != sfull[0].lower(): swarnx += 1
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
                            ilen = string.atoi(cs[3:-2])
                            qseq += '-' * ilen
                            sseq += sfull[:ilen]; sfull = sfull[ilen:]
                            aseq += '~' * ilen
                        elif cs[0] == '!':    # Missing in subject
                            ilen = string.atoi(cs[1:])
                            qseq += qfull[:ilen]; qfull = qfull[ilen:]
                            sseq += '^' * ilen
                            aseq += '^' * ilen
                        elif cs[0] == '?':    # Uncertain match
                            ilen = string.atoi(cs[1:])
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
                        raise

                #self.bugPrint(lentry)
                if qwarnx: self.warnLog('%s CS mismatch does not match loaded query %s sequence' % (qwarnx,lentry['Qry']),dev=self.debugging())
                if swarnx: self.warnLog('%s CS mismatch does not match loaded subject %s sequence' % (swarnx,lentry['Hit']),dev=self.debugging())
                if qfull: self.warnLog('Query sequence not consumed by alignment: %s' % qfull,dev=self.debugging())
                if sfull: self.warnLog('Subject sequence not consumed by alignment: %s' % sfull,dev=self.debugging())

                #if qwarnx or swarnx or qfull or sfull: self.debug(lentry)
                #else: self.debug('OK: %s' % lentry)

                #:37+ac:83+a:52+c:73*cg*gc:90-g:20-c:54-c:33+t:52+t:53-c:22-a:28*at*ta:15-at*at:25-aa:21+c:12-a:9-a:40+t:28-c:9*ac:64-ca:90-t:54+c:13-t:5-t:7-g:19-a:11-a:90+a:13+aa:14-c:52-g:7-c:7+a:51-t:22+t:11+t:11-t:44-a:12-ct*gt:66-c:62-g:23*ac:32+t:15+a:83-a:29+a:29+c:122+a:17+g:113-a:35+t:137-t:31+t:49+c:76+t:152+c:25+g:43+t:257+a:23+g:591+a:361-a:89+t:310+c:36*ca*at*tg:295+g:503+a:173-g:216+a:50-c:1332-a:1462-a:569-a:801-g:1716+g:34+a:657+a:55+a:30+a:240+a:2737-a:20062+a:685+c:10563+a:143986

                lentry['QrySeq'] = qseq
                lentry['AlnSeq'] = aseq
                lentry['SbjSeq'] = sseq

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
                self.warnLog('%s Length/AlnLen mismatches from PAF file. (Seems to be Mimimap2 bug. Run debug=T for details.)' % rje.iStr(alncorrx))
            if px != ltot and (alnseq or self.getBool('LocalAln')):
                raise ValueError('%s cs alignments missing from PAF file: check --cs flag was used' % rje.iStr(ltot-px))
            if minloclen > 1: locdb.dropEntries('Length<%d' % minloclen)
            if minlocid > 0:
                locdb.makeField('100.0*Identity/Length','LocID')
                locdb.dropEntries('LocID<%s' % minlocid)
                locdb.dropField('LocID')
            self.deBug('?')

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
    def localToUnique(self,locdb=None,save=True):  ### Convert local alignments to unique coverage of queries and hits.
        '''Convert local alignments to unique coverage of queries and hits.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','## ~~~~~ PAF Local Alignment to Unique Hits ~~~~~ ##')
            self.headLog('PAF Local Alignment to Unique Hits')
            if not locdb: locdb = self.db('local')
            locfields = string.split('Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd')

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
    def uniqueToHitSum(self,save=True):  ### Parse PAF Local alignment table into HitSum table.
        '''Parse PAF Local alignment table into HitSum table.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','## ~~~~~ PAF Local and Unique Hits to HitSum ~~~~~ ##')
            self.headLog('PAF Local and Unique Hits to HitSum')
            self.printLog('#INFO','Converting query-unique reduced local hits into hit summary statistics')
            locdb = self.db('local')
            if self.getBool('UniqueHit'):
                self.printLog('#INFO','Using hit-unique reduced local hits for GABLAM summary statistics')
                locdb = self.db('hitunique')
            locdb.index('Qry')
            uniqdb = self.db('qryunique')
            sumdb = self.db().copyTable(uniqdb,'hitsum')
            sumfields = string.split('Qry HitNum MaxScore EVal Description Length Coverage Identity Positives')
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
            #gabfields = string.split('Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID       Qry_AlnSim      Qry_Dirn        Qry_Start       Qry_End Qry_OrderedAlnLen       Qry_OrderedAlnID        Qry_OrderedAlnSim       Qry_OrderedDirn Qry_OrderedStart        Qry_OrderedEnd  Hit_AlnLen      Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End Hit_OrderedAlnLen       Hit_OrderedAlnID        Hit_OrderedAlnSim       Hit_OrderedDirn Hit_OrderedStart        Hit_OrderedEnd')
            gabfields = string.split('Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID Qry_AlnSim Qry_Dirn Qry_Start Qry_End Hit_AlnLen Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End')
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
            gabdb.keepFields( string.split('Qry Hit Rank Score EVal QryLen HitLen') )
            gabdb.list['Fields'] = gabfields

            ### ~ [3] Add GABLAM stats from Qry-Hit specifc unique mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = gabdb.entryNum(); filtx = 0
            for gentry in gabdb.entries():
                self.progLog('#GABLAM','Generating GABLAM table from local hits: %.2f%%' % (gx/gtot)); gx += 100
                self.quiet()
                #self.log.opt['Silent'] = True
                qry = gentry['Qry']
                hit = gentry['Hit']
                refdb = self.reduceLocal(locdb,byqry=True,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
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

                hitdb = self.reduceLocal(locdb,byqry=False,queries=[qry],hits=[hit],quiet=False,minlocid=self.getNum('MinLocID'),minloclen=self.getInt('MinLocLen'))
                if not hitdb.entryNum(): self.warnLog('Asymmetric loss of %s-%s alignments based on MinLocID' % (qry,hit),quitchoice=True,warntype='locfilter',suppress=True)
                for hentry in hitdb.entries():
                    if hentry['Strand'] == '-': [ hentry['SbjStart'], hentry['SbjEnd'] ] = [ hentry['SbjEnd'], hentry['SbjStart'] ]
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
                            if string.split('%s' % pid)[0] == 'WAIT': status = 1
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
                                    forkdb = self.reduceLocal(qhdb[qry][hit],byqry=qh=='qry',queries=[qry],hits=[hit],quiet=False,minloclen=self.getInt('MinLocLen'),minlocid=self.getNum('MinLocID'))
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
            gabfields = string.split('Qry Hit Rank Score EVal QryLen HitLen Qry_AlnLen Qry_AlnID Qry_AlnSim Qry_Dirn Qry_Start Qry_End Hit_AlnLen Hit_AlnID        Hit_AlnSim      Hit_Dirn        Hit_Start       Hit_End')
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
            gabdb.keepFields( string.split('Qry Hit Rank Score EVal QryLen HitLen') )
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
                    if hentry['Strand'] == '-': [ hentry['SbjStart'], hentry['SbjEnd'] ] = [ hentry['SbjEnd'], hentry['SbjStart'] ]
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
            self.debug(gabdb.entryNum())

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
                raise ValueError('"%s" is not a valid Local hit table field (%s)' % (sortfield,string.join(locdb.fields(),'|')))
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
                    except: raise ValueError('Could not detect minimap2 version: check minimap2=PROG setting (%s)' % self.getStr('Minimap2'))
                    self.printLog('#SYS','%s > %s' % (maprun,self.getStr('PAFIn')))
                    open('%s.cmd' % self.getStr('PAFIn'),'w').write('%s\n' % maprun)
                    os.system('%s > %s' % (maprun,self.getStr('PAFIn')))
                return rje.exists(self.getStr('PAFIn'))
            else:
                return os.popen(maprun).read()
        except: self.errorLog('%s.minimap2 error' % self.prog()); raise
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
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
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
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
