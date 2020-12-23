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
Module:       rje_genomics
Description:  Genomics data reformatting module
Version:      0.8.1
Last Edit:    07/08/20
Copyright (C) 2018  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

    reformat: convert TDT file into GFF or SAM
    ncbi: combine NCBI accession numbers and annotation with local data
    gffmap: convert GFF files from one ID set to another
    samfilt: filter read alignments from SAM file

Commandline:
    ### ~ General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    runmode=X       : Run mode (reformat/ftgff/ncbi/makemap/gffmap/diphap/fqreads/fas2bed/ncbinr/gapgff/locgff/samfilt) [reformat]
    basefile=X      : Base for output files, including log
    ### ~ Reformat Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    tdtfile=FILE    : Input delimited text file with data to convert [None]
    tdtkeys=LIST    : Input fields that define unique entries [Qry,Hit,AlnNum]
    queryfield=X    : Field defining the Query ("Read") name [Qry]
    targetfield=X   : Field defining the Target (Genome contig) name [Hit]
    begfield=X      : Field for beginning position [QryStart]
    endfield=X      : Field for end position [HitStart]
    reformat=X      : Output format (GFF3/SAM) [GFF3]
    ### ~ NCBI Annotation Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input assembly fasta file []
    seqstyle=X      : Sequence naming format for seqin=FASFILE sequences [dipnr]
    ncbifas=FASFILE : NCBI assembly fasta file []
    ncbigff=GFFFILE : NCBI annotation GFF file []
    ### ~ NCBINR Protein filtering Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input protein fasta file - accnum should match CDS feature Name []
    ncbigff=GFFFILE : NCBI annotation GFF file (locus naming format not important) []
    ### ~ GFF Map Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    gffs=FILELIST   : List of GFF files to convert - will be renamed BASEFILE.*.gff3 [*.gff,*.gff3]
    mapping=FILE    : File of old -> new ID mapping (e.g. from ncbi formatting) [mapping.csv]
    seqin=FASFILE   : Input assembly fasta file (makemap mode) []
    mapfas=FASFILE  : Alternative fasta file for GFF ID mapping (makemap mode) []
    ### ~ FASTQ Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    fqfiles=FILELIST: List of fastq files (may be gzipped) to process [*.fq,*.fq.gz,*.fastq,*.fastq.gz]
    ### ~ GAP GFF Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input assembly fasta file []
    mingap=INT      : Minimum gap length to annotation [10]
    ### ~ SAM Filter Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sam=FILE        : Input SAM file to filter []
    minmaplen=INT   : Minimum number of matching template positions to keep SAM hit [0]
    outsam=FILE     : Output SAM file [$BASEFILE.filtered.sam]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time, glob
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_samtools, rje_seqlist
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added ncbi annotation reformatting.
    # 0.2.0 - Added diphap renaming of pseduodiploid assembly sequences.
    # 0.3.0 - Added fqreads counting of reads from fastq.gz.
    # 0.4.0 - Added GFFMap function for mapping IDs onto others in GFF files. And mapfas mode to make mapping table.
    # 0.5.0 - Added Fas2Bed function for making a BED file of all contigs for bedtools coverage etc.
    # 0.6.0 - Added GapGFF option to generate GFF of assembly gaps.
    # 0.6.1 - Fixed TDTKeys bug.
    # 0.7.0 - Added loc2gff mode for converting local hits table to GFF3 output.
    # 0.8.0 - Added samfilt: filter read alignments from SAM file
    # 0.8.1 - Fixed gapgff bug that had first gap in header.
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
    # [ ] : Add parallel unzipping (pigz?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_GENOMICS', '0.8.1', 'August 2020', '2016')
    description = 'Misc genomics module'
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
### SECTION II: Genomics Class                                                                                          #
#########################################################################################################################
class Genomics(rje_obj.RJE_Object):
    '''
    Genomics Class. Author: Rich Edwards (2016).

    Str:str
    - BegField=X      : Field for beginning position [QryStart]
    - EndField=X      : Field for end position [HitStart]
    - MapFas=FASFILE  : Alternative fasta file for GFF ID mapping []
    - Mapping=FILE    : File of old -> new ID mapping (e.g. from ncbi formatting) [mapping.csv]
    - NCBIFas=FASFILE : NCBI assembly fasta file []
    - NCBIGFF=GFFFILE : NCBI annotation GFF file []
    - OutSAM=FILE     : Output SAM file [$BASEFILE.filtered.sam]
    - QueryField=X    : Field defining the Query ("Read") name [Qry]
    - Reformat=X      : Output format (GFF3/SAM) [GFF3]
    - RunMode=X       : Run mode (reformat/ncbi) [reformat]
    - SAM=FILE        : Input SAM file to filter []
    - SeqIn=FASFILE   : Input assembly fasta file []
    - SeqStyle=X      : Sequence naming format for seqin=FASFILE sequences [dipnr]
    - TargetField=X   : Field defining the Target (Genome contig) name [Hit]
    - TDTFile=FILE    : Input delimited text file with data to convert [None]

    Bool:boolean

    Int:integer
    - MinGap=INT      : Minimum gap length to annotation [10]
    - MinMapLen=INT   : Minimum number of matching template positions to keep SAM hit [0]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - FQFiles=FILELIST: List of fastq files (may be gzipped) to process [*.fq,*.fq.gz,*.fastq,*.fastq.gz]
    - GFFs=FILELIST   : List of GFF files to convert - will be renamed BASEFILE.*.gff3 [*.gff,*.gff3]
    - TDTKeys=LIST    : Input fields that define unique entries [Qry,Hit,AlnNum]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB
    - SeqList
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BegField','EndField','MapFas','Mapping','NCBIFas','NCBIGFF','QueryField','Reformat','RunMode',
                        'SeqIn','SeqStyle','TargetField','TDTFile','SAM','OutSAM']
        self.boollist = []
        self.intlist = ['MinGap','MinMapLen']
        self.numlist = []
        self.filelist = []
        self.listlist = ['FQFiles','GFFs','TDTKeys']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'BegField':'QryStart','EndField':'HitStart','QueryField':'Qry','Reformat':'GFF3',
                     'RunMode':'reformat','SeqStyle':'dipnr','Mapping':'mapping.csv',
                     'TargetField':'Hit'})
        self.setBool({})
        self.setInt({'MinGap':10})
        self.setNum({})
        self.list['TDTKeys'] = ['Qry','Hit','AlnNum']
        setfq = self.getStrLC('RunMode') == 'fqreads'
        for cmd in self.cmd_list: setfq = setfq and not cmd.lower().startswith('fqfiles=')
        if setfq: self._cmdReadList('fqfiles=*.fq,*.fq.gz,*.fastq,*.fastq.gz','glist',['FQFiles'])
        setgff = self.getStrLC('RunMode') == 'gffmap'
        for cmd in self.cmd_list: setgff = setgff and not cmd.lower().startswith('gffs=')
        if setgff: self._cmdReadList('gffs=*.gff,*.gff3','glist',['GFFs'])
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
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
                self._cmdReadList(cmd,'str',['BegField','EndField','NCBIFas','NCBIGFF','QueryField','Reformat','RunMode',
                        'SeqIn','SeqStyle','TargetField','TDTFile'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['MapFas','Mapping','TDTFile','SAM','OutSAM'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinGap','MinMapLen'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['TDTKeys'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['FQFiles','GFFs']) # List of files using wildcards and glob
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
            if self.getStrLC('RunMode') == 'loc2gff': self.loc2GFF3()
            if self.getStrLC('RunMode') == 'reformat':
                if self.getStrLC('Reformat') == 'gff3': self.saveGFF3()
                elif self.getStrLC('Reformat') == 'sam': self.saveSAM()
                else: raise ValueError('Reformat="%s" not recognised!' % self.getStr('Reformat'))
            if self.getStrLC('RunMode') == 'ncbi':
                self.snakeNCBI()    #!# Make this more generic!
            if self.getStrLC('RunMode') == 'diphap':
                self.dipHap()
            if self.getStrLC('RunMode') == 'fqreads':
                self.fastqcreads()
            if self.getStrLC('RunMode') == 'gffmap':
                self.gffMap()
            if self.getStrLC('RunMode') == 'makemap':
                self.makeMap()
            if self.getStrLC('RunMode') == 'fas2bed':
                self.fas2Bed()
            if self.getStrLC('RunMode') == 'ncbinr':
                self.NCBINR()
            if self.getStrLC('RunMode') == 'gapgff':
                self.gapGFF()
            if self.getStrLC('RunMode') == 'samfilt':
                self.samFilt()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('RunMode') in ['reformat','loc2gff']:
                self.db().addTable(filename=self.getStr('TDTFile'),mainkeys=self.list['TDTKeys'],name='input',expect=True)
                if not self.baseFile(return_none=None): self.baseFile(rje.baseFile(self.getStr('TDTFile')))
            #if self.getStrLC('RunMode') == 'ncbi':
            #    self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file'])
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
    ### <3> ### Data Saving Methods                                                                                     #
#########################################################################################################################
    def loc2GFF3(self,table='input'): ### Use BLAST to convert a loaded local table to GFF3
        '''Use BLAST to convert a loaded local table to GFF3.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #Qry	Hit	AlnNum	BitScore	Expect	Length	Identity	Positives	QryStart	QryEnd	SbjStart	SbjEnd
            #YARCTy1-1	ctgIA_MBGISH__MBGISH002.01	1	10619	0.0	5858	5822	0	68	5925	166826	160969
            blast = rje_blast.BLASTRun(log=self.log,cmd_list=['blastf=F']+self.cmd_list+['checkblast=F'])
            blast.obj['DB'] = self.obj['DB']
            gfffile = '%s.gff3' % rje.baseFile(self.getStr('TDTFile'))
            blast.saveGFF(gfffile,locdb=self.db(table),cdsmode=False)
        except: self.errorLog('%s.saveGFF3 error' % self.prog())
#########################################################################################################################
    def saveGFF3(self,table='input'):   ### Save data in GFF3 format
        '''
        Save data in GFF3 format.
        ### ~ Reformat Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        tdtfile=FILE    : Input delimited text file with data to convert [None]
        tdtkeys=LIST    : Input fields that define unique entries [Qry,Hit,AlnNum]
        queryfield=X    : Field defining the Query ("Read") name [Qry]
        targetfield=X   : Field defining the Target (Genome contig) name [Hit]
        begfield=X      : Field for beginning position [QryStart]
        endfield=X      : Field for end position [HitStart]
        reformat=X      : Output format (GFF3/SAM) [GFF3]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #Qry	Hit	AlnNum	BitScore	Expect	Length	Identity	Positives	QryStart	QryEnd	SbjStart	SbjEnd
            #YARCTy1-1	ctgIA_MBGISH__MBGISH002.01	1	10619	0.0	5858	5822	0	68	5925	166826	160969

            gffdb = self.db().copyTable(table,'%s.gff' % self.db(table).name())
            gfields = string.split('seqid source type start end score strand phase attributes')

            gffdb.newKey(['seqid','start','AlnID','end','type','source','attributes'])
            gffdb.keepFields(gfields+['AlnID'])
            ### ~ [2] Output GFF File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gcomments = ['##gff-version 3','#Generated by %s' % self.log.runDetails()]
            gcomments.append('#Full Command List: %s' % rje.argString(rje.tidyArgs(self.log.cmd_list)))
            gffdb.saveToFile(filename,delimit='\t',backup=True,append=append,savefields=gfields,log=True,headers=False,comments=gcomments)

            return ValueError('Sorry: save as GFF not yet implemented!')
        except: self.errorLog('%s.saveGFF3 error' % self.prog())
#########################################################################################################################
    def gapGFF(self):
        '''
        gapgff = Generating a GFF file with all the runs of gaps above a certain size.
        seqin=FASFILE   : Input assembly fasta file []
        mingap=INT      : Minimum gap length to annotation [10]

        pri49_NOTSC__NW_020716611.1     GapGFF  gap    7781    130001  .       -       .       ID=gap0;Note="xN gap"
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.getStr('SeqIn')
            mingap = self.getInt('MinGap')
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % seqin])
            ## ~ [0a] Output GFF File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gff = '%s.gap.gff3' % rje.baseFile(seqin)
            rje.backup(self,gff)
            GFF = open(gff,'w')
            gcomments = ['##gff-version 3','#Generated by %s' % self.log.runDetails()]
            gcomments.append('#Full Command List: %s\n' % rje.argString(rje.tidyArgs(self.log.cmd_list)))
            GFF.write(string.join(gcomments,'\n'))

            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0; gx = 0
            while seqlist.nextSeq():
                sx += 1; self.progLog('\r#SEQ','Processed %s sequences: %s gaps' % (rje.iStr(sx),rje.iStr(gx)))
                (name,sequence) = seqlist.getSeq(format='tuple')
                sname = seqlist.shortName()
                # Generate gaps
                sequence = string.replace(sequence.upper(),'N','|N|')
                #sequence = string.replace(sequence,'N||N','NN')
                sequence = string.join(string.split(sequence, '||'),'')
                gapseq = string.split(sequence,'|')
                # Process gaps
                pos = 0     # Position in sequence
                for chunk in gapseq:
                    startx = pos + 1
                    pos += len(chunk)
                    if chunk.startswith('N') and len(chunk) >= mingap:  # Add a gap
                        gstr = string.join([sname,'GapGFF','gap','%d' % startx,'%d' % pos,'.','+','.','ID=gap%d;Note=%sN gap' % (gx,rje.iStr(len(chunk)))],'\t')
                        GFF.write('%s\n' % gstr)
                        gx += 1
            self.printLog('\r#SEQ','Processed %s sequences: %s gaps' % (rje.iStr(sx),rje.iStr(gx)))
            GFF.close()
            self.printLog('#GFF','%s gaps output to %s' % (rje.iStr(gx),gff))

        except: self.errorLog('%s.gapGFF() error' % self.prog())
#########################################################################################################################
    def snakeNCBI(self):
        '''snakencbi = Reformatting the fasta and GFF files for snake genomes.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            babsfas = self.getStr('SeqIn')
            babsseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % babsfas])
            ncbifas = self.getStr('NCBIFas')
            ncbiseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % ncbifas])
            ncbigff = self.getStr('NCBIGFF')
            seqstyle = self.getStrLC('SeqStyle')
            mapping = {}    # Dictionary of GFF sequence ID mapping
            seqout = '%s.fasta' % self.baseFile()
            SEQOUT = open(seqout,'w')
            ### ~ [1] Read sequences and establish mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while ncbiseq.nextSeq():
                self.progLog('\r#SEQ','Processed %s primary sequences' % rje.iLen(mapping))
                babsseq.nextSeq()
                nseq = ncbiseq.currSeq()
                nacc = ncbiseq.seqAcc()
                bacc = babsseq.seqAcc()
                bdesc = babsseq.seqDesc()
                btype = string.split(bdesc)[0]
                bspec = babsseq.seqSpec()
                if seqstyle in ['dipnr']:
                    hapid = rje.matchExp('HAP(\d+)$',babsseq.shortName())[0]
                    newname = 'pri%s_%s__%s %s %s %s' % (hapid,bspec,nacc,btype,bacc,ncbiseq.seqDesc())
                else: raise ValueError('SeqStyle not recognised!')
                mapping[nacc] = string.split(newname)[0]
                SEQOUT.write('>%s\n%s\n' % (newname,nseq[1]))
            SEQOUT.close()
            self.printLog('\r#FAS','Processed %s primary sequences -> %s' % (rje.iLen(mapping),seqout))
            ### ~ [2] Output identity mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '%s.mapping.csv' % self.baseFile()
            MAP = open(mfile,'w')
            MAP.write('ncbi,mapname\n')
            for nacc in rje.sortKeys(mapping):
                MAP.write('%s,%s\n' % (nacc,mapping[nacc]))
            MAP.close()
            self.printLog('#MAP','%s ID mappings output to %s' % (rje.iLen(mapping),mfile))
            ### ~ [3] Update GFF file and split into sources ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gfile = ncbigff
            gout = {}   # Dictionary of source:FILE
            gx = {}; zx = 0
            headlines = []  # Header lines for GFF
            if gfile.endswith('.gz'):
                GFF = os.popen('zcat %s' % gfile)
            else:
                GFF = open(gfile,'r')
            gline = GFF.readline()
            while gline:
                if gline.startswith('#'): headlines.append(gline)
                else:
                    gdata = string.split(gline,'\t')
                    if len(gdata) < 2: continue
                    gsource = gdata[1]
                    if gdata[0] in mapping: gdata[0] = mapping[gdata[0]]
                    if gsource not in gout:
                        gout[gsource] = open('%s.ncbi.%s.gff3' % (self.baseFile(),gsource),'w')
                        gout[gsource].writelines(headlines)
                        gx[gsource] = 0
                    gout[gsource].write(string.join(gdata,'\t')); gx[gsource] += 1
                    zx += 1
                self.progLog('\r#GFF','Parsing %d GFF sources: %s lines' % (len(gx),rje.iStr(zx)),rand=0.01)
                gline = GFF.readline()
            self.printLog('\r#GFF','Parsing %d GFF sources: %s lines' % (len(gx),rje.iStr(zx)))
            for gsource in rje.sortKeys(gout):
                self.printLog('\r#GFF','Parsed %s %s GFF lines.' % (gsource,rje.iStr(gx[gsource])))
                gout[gsource].close()

        except:
            self.errorLog('Error in snakeNCBI',printerror=True,quitchoice=True)
#########################################################################################################################
    def makeMap(self):
        '''makemap = Make a mapping file from two fasta files.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infas = self.getStr('SeqIn')
            inseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % infas])
            mapfas = self.getStr('MapFas')
            mapseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s' % mapfas])
            mapping = {}    # Dictionary of GFF sequence ID mapping
            ### ~ [1] Read sequences and establish mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while inseq.nextSeq():
                self.progLog('\r#SEQ','Processed %s sequences' % rje.iLen(mapping))
                mapseq.nextSeq()
                mapping[inseq.shortName()] = mapseq.shortName()
            self.printLog('\r#MAP','Processed %s sequences.' % (rje.iLen(mapping)))
            ### ~ [2] Output identity mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '%s.mapping.csv' % self.baseFile()
            MAP = open(mfile,'w')
            MAP.write('original,mapped\n')
            for nacc in rje.sortKeys(mapping):
                MAP.write('%s,%s\n' % (nacc,mapping[nacc]))
            MAP.close()
            self.printLog('#MAP','%s ID mappings output to %s' % (rje.iLen(mapping),mfile))
        except:
            self.errorLog('Error in makeMap',printerror=True,quitchoice=True)
#########################################################################################################################
    def gffMap(self):
        '''gffmap = Map IDs from GFF files to make new ones.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapping = {}    # Dictionary of GFF sequence ID mapping
            for mline in open(self.getStr('Mapping'),'r').readlines()[1:]:
                mdata = string.split(rje.chomp(mline),',')
                mapping[mdata[0]] = mdata[1]
            self.printLog('#MAP','%s mapping IDs' % rje.iLen(mapping))
            ### ~ [1] Update GFF files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gff in self.list['GFFs']:
                GFF = open(gff,'r')
                newgff = string.split(rje.baseFile(gff,strip_path=True),'.')
                newgff[0] = self.baseFile()
                newgff = '%s.gff3' % string.join(newgff,'.')
                if newgff == gff:
                    try: raise ValueError(gff)
                    except: self.errorLog('New GFF file cannot have same name as original!')
                    continue
                rje.backup(self,newgff)
                NEWGFF = open(newgff,'w')
                gline = 'Go!'; gx = 0
                while gline:
                    gline = GFF.readline()
                    if gline.startswith('#'): NEWGFF.write(gline); gx += 1
                    else:
                        gdata = string.split(gline,'\t')
                        if len(gdata) < 2: continue
                        if gdata[0] in mapping: gdata[0] = mapping[gdata[0]]
                        else: continue
                        NEWGFF.write(string.join(gdata,'\t'))
                        gx += 1
                    self.progLog('\r#GFF','Parsing %s: %s lines' % (gff,rje.iStr(gx)),rand=0.001)
                NEWGFF.close()
                self.printLog('\r#GFF','Parsing %s %s lines -> %s' % (rje.iStr(gx),gff,newgff))

        except:
            self.errorLog('Error in snakeNCBI',printerror=True,quitchoice=True)
#########################################################################################################################
    def dipHap(self):
        '''Update sequence names.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list)
            SEQOUT = open(self.baseFile() + '.dipnr.fasta','w')
            diplist = []
            haplist = []
            for seq in seqlist.seqs():
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in diplist: haplist.append(hap)
                else: diplist.append(hap)
            for seq in seqlist.seqs():
                sname = string.split(seqlist.shortName(seq),'_')
                hap = rje.matchExp('HAP(\d+)',seqlist.shortName(seq))
                if hap in haplist:
                    if hap in diplist: haptxt = 'haploidA'; diplist.remove(hap); sname[0] = 'pri%s' % hap
                    else: haptxt = 'haploidB'; sname[0] = 'alt%s' % hap
                else: haptxt = 'diploid'; sname[0] = 'pri%s' % hap
                sname = string.join(sname,'_')
                SEQOUT.write('>%s %s %s\n%s\n' % (sname,haptxt,seqlist.seqDesc(seq),seqlist.seqSequence(seq)))
                self.printLog('#DIP','%s: %s\n' % (sname,haptxt))
            SEQOUT.close()
        except:
            self.log.errorLog('Error in rje_genomics.dipHap()',printerror=True,quitchoice=True)
#########################################################################################################################
    def NCBINR(self):
        '''
        Generate a BED file of contigs, which can then be used with bedtools, e.g.

        awk '$3 == "tRNA"' ts10xv2.ncbipri.ncbi.tRNAscan-SE.gff3 | bedtools coverage -a test.bed -b stdin
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list)
            seqdict = seqlist.makeSeqNameDic('accnum')
            GFF = open(self.getStr('NCBIGFF'), 'r')
            gffprot = []  # List of (locus, strand, start, end, GeneID)  gdata[0,6,3,4,<8>]
            gene2prot = {} # Dictionary of {GeneID: accnum list}
            parents = {}  # CDS to RNA

            ### ~ [1] Parse GFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gline = 'Go!'; px = 0
            while gline:
                gline = rje.chomp(GFF.readline())
                if gline.startswith('#') or not gline: continue
                gdata = string.split(gline, '\t')
                if len(gdata) < 9: continue
                gtype = gdata[2]
                if gtype not in ['gene','CDS']: continue
                self.progLog('\r#GFF','Parsing %s GFF genes; %s proteins' % (rje.iLen(gffprot),rje.iStr(px)))
                try: geneid = rje.matchExp('GeneID:(\d+)[,;]',gdata[8])[0]
                except: raise ValueError(gdata[8])
                if gtype == 'CDS':
                    try: acc = rje.matchExp('Name=(XP_\d+\.\d+);',gdata[8])[0]
                    except: self.warnLog(gdata[8]); continue
                    parent = rje.matchExp('Parent=(rna\d+);',gdata[8])[0]
                    parents[acc] = parent
                    if geneid not in gene2prot: gene2prot[geneid] = []
                    if acc not in gene2prot[geneid]: gene2prot[geneid].append(acc); px += 1
                else:
                    (glocus, gstrand, gstart, gend) = (gdata[0], gdata[6], int(gdata[3]), int(gdata[4]))
                    gffprot.append( (glocus, gstrand, gstart, gend, geneid) )
            self.printLog('\r#GFF','Parsed %s GFF genes; %s proteins.' % (rje.iLen(gffprot),rje.iStr(px)))

            ### ~ [2] Sort and reduce ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Longest isoform per gene -> update gene2prot to single sequence
            for geneid in gene2prot:
                if len(geneid) == 1: gene2prot[geneid] = gene2prot[geneid][0]
                else:
                    longx = 0; best = 0
                    for acc in gene2prot[geneid]:
                        if seqlist.seqLen(seqdict[acc]) > longx:
                            longx = seqlist.seqLen(seqdict[acc])
                            best = acc
                    if not best: raise ValueError('GeneID %s has no best protein sequence!' % geneid)
                    gene2prot[geneid] = best
            #i# Remove identical genes in identical positions
            badgenes = []
            gffprot.sort()
            i = 0
            while i < (len(gffprot) - 1):
                j = i + 1
                while j < len(gffprot) and (gffprot[i][0],gffprot[i][1]) == (gffprot[j][0],gffprot[j][1]):
                    if (gffprot[i][3] < gffprot[j][2]): break # Not overlapping
                    try:
                        if seqlist.seqSequence(seqdict[gene2prot[gffprot[i][4]]]) == seqlist.seqSequence(seqdict[gene2prot[gffprot[j][4]]]):
                            badgenes.append(gene2prot[gffprot[j][4]])
                            self.printLog('#NR','%s and %s are overlapping in GFF with identical sequence!' % (seqdict[gene2prot[gffprot[i][4]]],seqdict[gene2prot[gffprot[j][4]]]))
                    except: pass #!# Add warning
                    j += 1
                i += 1
            #i# Drop bad genes
            badgenes = rje.sortUnique(badgenes)
            self.printLog('\r#NR','%s redundant genes to remove' % rje.iLen(badgenes))
            for geneid in rje.sortKeys(gene2prot):
                if geneid in badgenes: gene2prot.pop(geneid)
            goodacc = gene2prot.values()
            goodpar = []
            for acc in goodacc: goodpar.append(parents[acc])

            ### ~ [3] Cleanup fasta and GFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fasout = '%s.prot.nr.fas' % self.baseFile()
            gffout = '%s.nr.gff3' % self.baseFile()
            FASOUT = open(fasout,'w'); sx = 0
            for seq in seqlist.seqs():
                if seqlist.seqAcc(seq) in goodacc:
                    FASOUT.write('>%s\n%s\n' % seqlist.getSeq(seq)); sx += 1
            self.printLog('#FAS','%s proteins output to %s' % (rje.iStr(sx),fasout))
            FASOUT.close()

            GFFOUT = open(gffout,'w')
            GFF.seek(0)
            gline = 'Go!'; gx = 0
            while gline:
                gline = GFF.readline()
                if not gline: continue
                if gline.startswith('#'): GFFOUT.write(gline)
                gdata = string.split(gline, '\t')
                if len(gdata) < 9: continue
                geneid = rje.matchExp('GeneID:(\d+)[,;]',gdata[8])[0]
                if geneid not in gene2prot: continue
                self.progLog('\r#GFF','%s feature lines output to %s' % (rje.iStr(gx),gffout),rand=0.001)
                gtype = gdata[2]
                if gtype == 'gene': GFFOUT.write(gline); gx += 1
                elif gtype == 'CDS':
                    acc = rje.matchExp('Name=(XP_\d+\.\d+);',gdata[8])[0]
                    if acc in goodacc: GFFOUT.write(gline); gx += 1
                elif gtype in ['mRNA','transcript']:
                    parent = rje.matchExp('ID=(rna\d+);',gdata[8])[0]
                    if parent in goodpar: GFFOUT.write(gline); gx += 1
                else:
                    try: parent = rje.matchExp('Parent=(rna\d+);',gdata[8])[0]
                    except: raise ValueError(gdata[8])
                    if parent in goodpar: GFFOUT.write(gline); gx += 1
            self.printLog('\r#GFF','%s feature lines output to %s' % (rje.iStr(gx),gffout))
            GFFOUT.close()

        except:
            self.log.errorLog('Error in rje_genomics.NCBINR()',printerror=True,quitchoice=True)
#########################################################################################################################
    def fas2Bed(self):
        '''
        Reduce protein fasta file and GFF file based on redundant annotations.

        >prot_PSETE__XP_026558790.1 trans-acting T-cell-specific transcription factor GATA-3, partial [Pseudonaja textilis]
        >prot_PSETE__XP_026562786.1 CUGBP Elav-like family member 2 isoform X1 [Pseudonaja textilis]
        >prot_PSETE__XP_026562874.1 CUGBP Elav-like family member 2 isoform X2 [Pseudonaja textilis]

        pri49_NOTSC__NW_020716611.1     Gnomon  gene    7781    130001  .       -       .       ID=gene0;Dbxref=GeneID:113411864;Name=LOC113411864;gbkey=Gene;gene=LOC113411864;gene_biotype=protein_coding;partial=true;start_range=.,7781
        pri49_NOTSC__NW_020716611.1     Gnomon  mRNA    7781    130001  .       -       .       ID=rna0;Parent=gene0;Dbxref=GeneID:113411864,Genbank:XM_026667177.1;Name=XM_026667177.1;gbkey=mRNA;gene=LOC113411864;model_evidence=Supporting evidence includes similarity to: 6 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 11 samples with support for all annotated introns;partial=true;product=dedicator of cytokinesis protein 1;start_range=.,7781;transcript_id=XM_026667177.1
        pri49_NOTSC__NW_020716611.1     Gnomon  exon    129898  130001  .       -       .       ID=id1;Parent=rna0;Dbxref=GeneID:113411864,Genbank:XM_026667177.1;gbkey=mRNA;gene=LOC113411864;partial=true;product=dedicator of cytokinesis protein 1;transcript_id=XM_026667177.1
        pri49_NOTSC__NW_020716611.1     Gnomon  CDS     129898  129983  .       -       0       ID=cds0;Parent=rna0;Dbxref=GeneID:113411864,Genbank:XP_026522962.1;Name=XP_026522962.1;gbkey=CDS;gene=LOC113411864;partial=true;product=dedicator of cytokinesis protein 1;protein_id=XP_026522962.1


        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list)   # This should be the protein sequence
            BEDOUT = open(self.baseFile() + '.bed','w')
            GENOUT = open(self.baseFile() + '.genome','w')
            self.progLog('#SEQ','Processing...')
            for seq in seqlist.seqs():
                BEDOUT.write('%s\t0\t%d\n' % (seqlist.shortName(seq),seqlist.seqLen(seq)))
                GENOUT.write('%s\t%d\n' % (seqlist.shortName(seq),seqlist.seqLen(seq)))
            BEDOUT.close()
            self.progLog('\r#BED',self.baseFile() + '.bed saved.')
            GENOUT.close()
            self.progLog('\r#SEQLEN',self.baseFile() + '.genome saved.')
        except:
            self.log.errorLog('Error in rje_genomics.fas2Bed()',printerror=True,quitchoice=True)
#########################################################################################################################
    def fastqcreads(self):
        '''Update sequence names.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #for fqgz in glob.glob('*.fq.gz')+glob.glob('*.fastq.gz')+glob.glob('*.fq')+glob.glob('*.fastq'):
            for fqgz in self.list['FQFiles']:
                x = 0; bp = 0; n = 0
                if fqgz.endswith('.gz'): FQ = os.popen('zcat %s' % fqgz)
                else: FQ = open(fqgz,'r')
                fline = None
                self.progLog('\r#FASTQ','Parsing %s reads: %dk' % (fqgz,n/1000))
                while fline or x == 0:
                    fline = rje.chomp(FQ.readline()); x += 1
                    if x % 4 == 2:
                        bp += len(fline); n += 1
                        if not n % 10000: self.progLog('\r#FASTQ','Parsing %s reads: %.2fM' % (fqgz,n/1e6))
                self.printLog('\r#FASTQ','Parsed %s reads: %s -> %s bp' % (fqgz,rje.iStr(n),rje.iStr(bp)))
                self.printLog('\r#RLEN','Mean %s read length: %.2f' % (fqgz,1.0*bp/n))
                FQ.close()
        except:
            self.errorLog('Error in rje_genomics.fastqcreads()',printerror=True,quitchoice=True)
#########################################################################################################################
    def samFilt(self):  ### Filters SAM file into another SAM file
        '''
        Filters SAM file into another SAM file.
        >> filename:str = Pileup file name
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SAM'): raise ValueError('Need to set sam=FILE for runmode=samfilt')
            filename = self.getStr('SAM')
            if not rje.exists(filename): raise IOError('sam={0} file missing'.format(filename))
            if filename.endswith('.bam'): SAM = os.popen('samtools view %s' % filename)
            else: SAM = open(filename,'r')
            if self.getStrLC('OutSAM'): outfile = self.getStr('OutSAM')
            else: outfile = '{0}.filtered.sam'.format(self.baseFile())
            rje.backup(self,outfile,appendable=False)
            OUT = open(outfile,'w')
            rid = 0         # Read counter (ID counter)
            unmappedx = 0   # Counter of unmapped reads
            filtx = 0       # Counter of filtered reads
            oldblasrwarn = None
            ## ~ [1a] Filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            minmaplen = self.getInt('MinMapLen')
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in SAM:
                ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.startswith('@'): OUT.write(line); continue
                samdata = string.split(line)
                if len(samdata) < 11: OUT.write(line); continue
                self.progLog('\r#SAM','Parsing %s: %s reads (+ %s unmapped); %s filtered...' % (filename,rje.iStr(rid),rje.iStr(unmappedx),rje.iStr(filtx)),rand=0.1)
                cigstr = samdata[5]
                if cigstr == '*': unmappedx += 1; continue
                rid += 1
                rname = samdata[0]
                #!# Add additional filtering options here
                try:
                    cigdata = rje_samtools.parseCigar(cigstr)
                except:
                    self.errorLog('Problem with cigar string (%s): "%s"' % (rname,cigstr))
                    continue
                if rje_samtools.cigarAlnLen(cigdata) < minmaplen: filtx += 1; continue
                OUT.write(line)
            self.printLog('#SAM','Parsed %s: %s filtered reads output to %s; %s filtered.' % (filename,rje.iStr(rid),outfile,rje.iStr(filtx)))
            if unmappedx: self.printLog('#SAM','Parsed and rejected %s unmapped reads.' % (rje.iStr(unmappedx)))
            if oldblasrwarn: self.warnLog(oldblasrwarn)
            SAM.close()
            OUT.close()
            return True
        except: self.errorLog('%s.parseSAM() error' % (self.prog())); return False
#########################################################################################################################
 ### End of SECTION II: Genomics Class                                                                                   #
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
    try: Genomics(mainlog,cmd_list).run()

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
