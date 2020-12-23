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
Module:       rje_exonerate
Description:  Runs Exonerate and parses output
Version:      0.5.1
Last Edit:    14/09/18
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    exopt=X         : Exonerate options to add. Will split on whitespace. No semicolons allowed. []
    memsaver=T/F    : Whether to run in memsaver mode that stores exonerate output in intermediate text file [False]
    cleanup=T/F     : Remove intermediate files after run if memsaver=T [False]
    gzip=T/F        : Whether to gzip (and gunzip) exonerate results files (not Windows) [True]

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, subprocess, sys, time
from subprocess import Popen, PIPE
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_sequence
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial Working Version
    # 0.1.1 - Fixed no model bug for DNA searches.
    # 0.2.0 - Added SeqList objects to store compiled aligned sequences.
    # 0.3.0 - Added Exonerate Options.
    # 0.4.0 - Added MemSaver mode.
    # 0.4.1 - Fixed bug in GFF output that had ID and Name uniqueness swapped.
    # 0.5.0 - Added QryDesc to GFF output.
    # 0.5.1 - Fixed bug in hitsum combined score.
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
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_EXONERATE', '0.5.1', 'Sept 2018', '2018')
    description = 'Runs Exonerate and parses output'
    author = 'Dr Richard J. Edwards & Timothy G. Amos.'
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
### SECTION II: Exonerate Class                                                                                               #
#########################################################################################################################
class Exonerate(rje_obj.RJE_Object):
    '''
    Exonerate Class. Author: Rich Edwards (2018).

    Str:str
    - ExOpt = Additional exonerate options
    
    Bool:boolean
    - MemSaver=T/F    : Whether to run in memsaver mode that stores exonerate output in intermediate text file [False]
    - Cleanup=T/F     : Remove intermediate files after run if memsaver=T [False]
    - GZip=T/F        : Whether to gzip (and gunzip) exonerate results files (not Windows) [True]

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DNAHits = rje_seqlist.SeqList object storing sequence tuples
    - ProtHits = rje_seqlist.SeqList object storing sequence tuples
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['ExOpt']
        self.boollist = ['MemSaver','Cleanup','GZip']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'MemSaver':False,'Cleanup':False,'GZip':True})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
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
                self._cmdReadList(cmd,'str',['ExOpt'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['MemSaver','Cleanup','GZip'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
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
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def exonerate(self,qryfas, genome, model,exonerate='exonerate',bestn=0):
        '''
        Runs exonerate and parses output into lists for processing.
        { query: {'gff':[outputlines], 'cigar':[outputlines], 'alignment':[outputlines], 'vulgar':[[headerlist], {header:value}, {header:value}, ...] }
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            EXFILE = None
            exfile = '%s.%s' % (self.baseFile(),model)  # Used in memsaver mode
            query_dic = {}
            header_list = ['query_id', 'query_start', 'query_end', 'query_strand', 'target_id', 'target_start', 'target_end', 'target_strand', 'score', '<label, query_length, target_length> triplets']
            excmd = [exonerate, qryfas, genome, '--showtargetgff', '--showcigar']
            if model: excmd += ['--model', model]
            if bestn: excmd += ['--bestn', '%d' % bestn]
            if self.getStrLC('ExOpt'): excmd += string.split(self.getStr('ExOpt'))
            self.printLog('#RUN',string.join(excmd))
            extext = []
            if self.getBool('MemSaver'):
                gzfile = '%s.gz' % exfile
                if rje.exists(gzfile): self.gUnzip(gzfile)
                if rje.exists(exfile) and not self.force():
                    self.printLog('#EXFILE','Found %s (force=F). Assuming complete.' % exfile)
                else:
                    rje.backup(self,exfile)
                    self.printLog('#SAVER','memsaver=T: Exonerate output directed to %s.' % exfile)
                    EXFILE = open(exfile,'w')
                    if subprocess.call(excmd, stdout=EXFILE): raise IOError('Exonerate call did not complete!')
                    EXFILE.close()
                    self.printLog('#EXFILE','%s generated.' % exfile)
                EXFILE = open(exfile,'r')
            else:
                extext = Popen(excmd, stdout=PIPE).stdout.readlines()
            output_format = ''
            while extext or EXFILE:
                #line = process.stdout.readline().rstrip()
                if EXFILE:
                    line = EXFILE.readline()
                    if not line: break
                    line = rje.chomp(line)
                else: line = rje.chomp(extext.pop(0))
                if line:
                    if line.startswith('         Query:'):
                        query = line.split(':', 1)[1].split(' ')[1]
                        #for q in rje.sortKeys(query_dic):
                        #    self.bugPrint('%s: %s' % (q,rje.sortKeys(query_dic[q])))
                        #self.debug(query)
                    if line == 'C4 Alignment:':
                        output_format = 'alignment'
                    elif line == '# --- START OF GFF DUMP ---':
                        output_format = 'gff'
                    elif line.startswith('vulgar:'):
                        output_format = 'vulgar'
                        fields = line.split(' ', 10)[1:]
                        if output_format in query_dic[query]:
                            query_dic[query][output_format].append({})
                        else:
                            query_dic[query][output_format] = [header_list, {}]
                        for header, field in zip(header_list, fields):
                            query_dic[query][output_format][-1][header] = field
                        #self.debug(query_dic[query][output_format])
                    elif line.startswith('cigar:'):
                        output_format = 'cigar'
                        if output_format in query_dic[query]:
                            query_dic[query][output_format].append(line.replace('cigar: ', ''))
                        else:
                            query_dic[query][output_format] = [line.replace('cigar: ', '')]
                    elif line == '------------' or line.startswith('Command line:') or line.startswith('Hostname:') or line == '# --- END OF GFF DUMP ---' or line == '#' or line.startswith('-- completed exonerate analysis'):
                        pass
                    elif output_format:
                        if query in query_dic:
                            if output_format in query_dic[query]:
                                query_dic[query][output_format].append(line)
                            else:
                                query_dic[query][output_format] = [line]
                        else:
                            query_dic[query] = {output_format:[line]}
                #elif process.poll() is not None:
                #    break
                elif output_format == 'alignment':
                    try: query_dic[query][output_format].append(line)
                    except: pass
                self.vPrint(line,v=1)
            if EXFILE:
                EXFILE.close()
                if self.getBool('Cleanup'):
                    os.unlink(exfile)
                    self.printLog('#CLEAN','%s deleted.' % exfile)
                elif self.getBool('GZip'): self.gZip(exfile)
            return query_dic
        except: self.errorLog('%s.exonerate error' % self.prog()); raise
#########################################################################################################################
    def alignmentToLocal(self,alignment=[],protqry=False):    ### Converts alignment into local hits table
        '''
        Converts alignment into local hits table.
        >> alignment:list of alignment text strings parsed from exonerate output.
        >> protqry:bool[False] = Whether query is protein
        << returns local database table.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vfields = ['Qry','Hit','AlnID','Score','Expect','Length','Identity','Positives','QryStart','QryEnd','HitStart','HitEnd','QrySeq','HitSeq','AlnSeq','Rank','Phase','HitStrand']
            vdb = self.db().addEmptyTable('local',vfields,['Qry','Hit','AlnID'])

            ### ~ [2] Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            '''
                     Query: FAXD1_NOTSC (P82807) Venom prothrombin activator notecarin-D1 [Notechis scutatus scutatus]
                    Target: ahap_PSETE__EBS10XV2AHAP187 haploidB edges=694320..157489 left=833615 right=281503 ver=1.9 style=4:[revcomp]
                     Model: protein2genome:local
                 Raw score: 1170
               Query range: 19 -> 295
              Target range: 12312786 -> 12307250
            
                   20 : AlaGluSerAsnValPheLeuLysSerLysValAlaAsnArgPheLeuGlnArg :       37
                        ..!...|||   ||||||||||||||||||||||||||||||||||||||||||
                        CysSerSerLeuValPheLeuLysSerLysValAlaAsnArgPheLeuGlnArg
             12312786 : TGTTCTTCTTTAGTATTCTTAAAAAGCAAAGTGGCAAATAGATTTTTGCAAAGA : 12312735
            
                  264 : {G}  >>>> Target Intron 7 >>>>  {ly}GluIleAspIleSerArg :      270
                        {|}           1304 bp           {||}|||||||||||||||!!!
                        {G}++                         ++{ly}GluIleAspIleSerSer
             12308652 : {G}gt.........................ag{GG}GAAATAGACATATCAAGC : 12307328
            
                  289 : ValProProAsnTyrTyrTyr :      295
                        |||||| !!!..||| !!|||
                        ValProAlaThrTyrAspTyr
             12307273 : GTTCCTGCCACGTATGACTAT : 12307251
            '''
            qry = None
            hit = None
            alnx = {}
            ventry = {}
            parsing = alignment[0:]
            rank = 1

            while parsing:
                line = parsing.pop(0)
                #self.bugPrint(line)
                # Query
                if rje.matchExp('Query: (\S+)',line):
                    if ventry: vdb.addEntry(ventry)
                    ventry = {'Qry':rje.matchExp('Query: (\S+)',line)[0],'QrySeq':'','HitSeq':'','AlnSeq':'','Rank':rank}
                    rank += 1
                # Hit
                if rje.matchExp('Target: (\S+)',line):
                    ventry['Hit'] = rje.matchExp('Target: (\S+)',line)[0]
                    qh = (ventry['Qry'],ventry['Hit'])
                    if qh in alnx: alnx[qh] += 1
                    else: alnx[qh] = 1
                    ventry['AlnID'] = alnx[qh]
                # Score
                if rje.matchExp('core: (\S+)',line):
                    ventry['Score'] = int(rje.matchExp('core: (\S+)',line)[0])
                # Alignment
                if rje.matchExp('^\s+(\d+) : (.+) :\s+(\d+)',line):
                    adata = rje.matchExp('^\s+(\d+) : (.+) :\s+(\d+)',line)
                    #self.bugPrint('= new aln: %s ->  %s' % (adata[0],adata[2]))
                    start = int(adata[0])
                    end = int(adata[2])
                    aln = adata[1]
                    x = line.find(aln)
                    if 'QryStart' not in ventry: ventry['QryStart'] = start
                    ventry['QryEnd'] = end
                    ventry['QrySeq'] += aln
                    #self.bugPrint('^%s$' % ventry['QrySeq'])

                    line = parsing.pop(0)
                    #self.bugPrint(line)
                    #self.bugPrint(']%s[' % aln)
                    #self.bugPrint(']%s[' % line[x:x+len(aln)])
                    ventry['AlnSeq'] += line[x:x+len(aln)]
                    #self.debug('^%s$' % ventry['AlnSeq'])

                    #self.bugPrint(parsing[0])
                    adata = rje.matchExp('^\s+(\d+) : (.+) :\s+(\d+)',parsing.pop(0))
                    if not adata:
                        #self.deBug(parsing[0])
                        adata = rje.matchExp('^\s+(\d+) : (.+) :\s+(\d+)',parsing.pop(0))
                    if not adata: raise ValueError('Partial alignment! Truncated output?')
                    #self.bugPrint('+ hit aln: %s ->  %s' % (adata[0],adata[2]))
                    start = int(adata[0])
                    end = int(adata[2])
                    aln = adata[1]
                    if 'HitStart' not in ventry: ventry['HitStart'] = start
                    ventry['HitEnd'] = end
                    ventry['HitSeq'] += aln
            if ventry: vdb.addEntry(ventry)
            ## Seq Check
            for ventry in vdb.entries():
                #self.bugPrint('^%s$' % ventry['QrySeq'])
                #self.bugPrint('^%s$' % ventry['AlnSeq'])
                #self.bugPrint('^%s$' % ventry['HitSeq'])
                if len(ventry['QrySeq']) != len(ventry['AlnSeq']) or len(ventry['QrySeq']) != len(ventry['HitSeq']):
                    self.debug(ventry)
                    raise ValueError('Alignment sequence length mismatch! Qry:%d ; Aln:%d ; Hit:%d' % (len(ventry['QrySeq']),len(ventry['AlnSeq']),len(ventry['HitSeq'])))

            ### ~ [3] Split on introns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DNAHits'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=None','seqmode=tuple','autoload=F','dna=T'])
            self.obj['ProtHits'] = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=None','seqmode=tuple','autoload=F'])

            #i# Protein Position Conversion
            if protqry:
                for ventry in vdb.entries():
                    # 1->1, 2->4, 3->7 = 1+3*(n-1)
                    ventry['QryStart'] = 1+3*(ventry['QryStart']-1)
                    if ventry['QrySeq'].startswith('{'):
                        codend = ventry['QrySeq'].find('}')
                        # {X} = phase 2, find = 2
                        if codend == 2: ventry['QryStart'] += 2
                        # {XX} = phase 1, find = 3
                        elif codend == 3: ventry['QryStart'] += 1
                        else: raise ValueError('QrySeq {} bracket mismatch!: %s' % ventry)
                    ventry['QryEnd'] = ventry['QryStart'] + len(ventry['QrySeq']) - string.count(ventry['QrySeq'],'-') - 1

            vdb.newKey(['Qry','Rank','Hit','AlnID'])
            for vkey in vdb.dataKeys():
                ventry = vdb.data(vkey)
                #i# Make a combined hitseq to output to fasta
                #># phap_PSETE__EBS10XV2PHAP187.FAXD1_NOTSC.XXX
                hitname = '%s.ex%s %s %s-%s' % (ventry['Qry'],ventry['Rank'],ventry['Hit'],rje.iStr(ventry['HitStart']),rje.iStr(ventry['HitEnd']))
                hitseq = ''
                phase = (ventry['QryStart'] + 2) % 3
                alnx = 1
                vkeyentries = [ventry]
                dirn = 1
                if ventry['HitEnd'] < ventry['HitStart']:
                    dirn = -1
                    ventry['HitStrand'] = '-'
                else: ventry['HitStrand'] = '+'
                for seq in ['HitSeq','QrySeq','AlnSeq']:
                    ventry[seq] = string.replace(ventry[seq],'}','')
                    ventry[seq] = string.replace(ventry[seq],'{','')
                while rje.matchExp('(\s+>>>> Target Intron \d+ >>>>\s+)',ventry['QrySeq']):
                    intron = rje.matchExp('(\s+>>>> Target Intron \d+ >>>>\s+)',ventry['QrySeq'])[0]
                    x = ventry['QrySeq'].find(intron)
                    y = x + len(intron)
                    intronlen = int(rje.matchExp('(\d+) bp',ventry['AlnSeq'][x:y])[0])
                    #i# Create a new entry of the first exon
                    newentry = rje.combineDict({},ventry)
                    for seq in ['HitSeq','QrySeq','AlnSeq']:
                        newentry[seq] = newentry[seq][:x]
                    newentry['AlnID'] = '%s.%d' % (ventry['AlnID'],alnx); alnx += 1
                    newentry['QryEnd'] = newentry['QryStart'] + len(newentry['QrySeq']) - string.count(newentry['QrySeq'],'-') - 1
                    newentry['HitEnd'] = newentry['HitStart'] + (len(newentry['HitSeq']) - string.count(newentry['HitSeq'],'-') - 1) * dirn
                    newentry['Length'] = x
                    newentry['Identity'] = string.count(newentry['AlnSeq'],'|')
                    vkeyentries.append(vdb.addEntry(newentry))
                    hitseq += newentry['HitSeq']
                    #i# Update ventry to be the rest of the hit
                    for seq in ['HitSeq','QrySeq','AlnSeq']:
                        ventry[seq] = ventry[seq][y:]
                    ventry['QryStart'] = newentry['QryEnd'] + 1
                    if protqry: ventry['QryEnd'] = ventry['QryStart'] + len(ventry['QrySeq']) - string.count(ventry['QrySeq'],'-') - 1
                    ventry['HitStart'] = newentry['HitEnd'] + intronlen * dirn
                #i# Calculate length and identity of final exon
                ventry['AlnID'] = '%s.%d' % (ventry['AlnID'],alnx)
                ventry['Length'] = len(ventry['AlnSeq'])
                ventry['Identity'] = string.count(ventry['AlnSeq'],'|')
                #i# Add sequence hits
                hitname += ' (%d alignment blocks)' % alnx
                hitseq += ventry['HitSeq']
                hitseq = string.replace(hitseq,'-','')
                protseq = rje_sequence.dna2prot('%s%s' % ('N' * phase,hitseq))
                self.obj['ProtHits']._addSeq(hitname,protseq)
                if ventry['HitStart'] > ventry['HitEnd']: hitseq = rje_sequence.reverseComplement(hitseq)
                self.obj['DNAHits']._addSeq(hitname,hitseq)

                #i# Update AlnID for proper float sorting
                for ventry in vkeyentries:
                    (vcore,vx) = string.split(ventry['AlnID'],'.')
                    ventry['AlnID'] = '%s.%s' % (vcore,rje.preZero(int(vx),alnx))
                    #self.debug(ventry)
            vdb.dataFormat({'AlnID':'string'})
            vdb.remakeKeys()
            self.debug(vdb.dataKeys())

            ## Seq Check
            for ventry in vdb.entries():
                #self.bugPrint('^%s$' % ventry['QrySeq'])
                #self.bugPrint('^%s$' % ventry['AlnSeq'])
                #self.bugPrint('^%s$\n' % ventry['HitSeq'])
                if len(ventry['QrySeq']) != len(ventry['AlnSeq']) or len(ventry['QrySeq']) != len(ventry['HitSeq']):
                    self.debug(ventry)
                    raise ValueError('Alignment sequence length mismatch! Qry:%d ; Aln:%d ; Hit:%d' % (len(ventry['QrySeq']),len(ventry['AlnSeq']),len(ventry['HitSeq'])))

            udb = self.reduceLocal(byqry=True)
            udb.rename('unique')
            udb.newKey(['Qry','Rank','Hit','AlnID'])
            self.debug(vdb.dataKeys())

            #i# Calculate exon phase
            for ventry in vdb.entries() + udb.entries(): ventry['Phase'] = (ventry['QryStart'] - 1) % 3

            #i# Protein Position Conversion
            if protqry:
                for ventry in vdb.entries():
                    ventry['QryStart'] = (ventry['QryStart']+2)/3
                    ventry['QryEnd'] = (ventry['QryEnd']+2)/3
                for ventry in udb.entries():
                    ventry['QryStart'] = (ventry['QryStart']+2)/3
                    ventry['QryEnd'] = (ventry['QryEnd']+2)/3

            #vdb.remakeKeys()
            return vdb

        except: self.errorLog('%s.alignmentToLocal error' % self.prog()); raise
#########################################################################################################################
    def reduceLocal(self,minloclen=0,minlocid=0,byqry=False):    ### Reduces local alignments to cover each hit region only once
        '''
        Reduces local BLAST alignments to cover each hit region only once. Uses BLAST object to do this, so first converts
        table to BLAST local table format, reduces, then converts back.
        @param minloclen:int [0] = Minimum local length to keep.
        @param minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        @param byqry:bool [False] = Whether to reduce Local table by Query rather than hit
        @return: copy of local table, filtered and reduced.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Generate a BLAST-compatible format; Rename fields and make AlnID an integer
            fieldswap = [('Qry','Query'), ('HitStart','SbjStart'), ('HitEnd','SbjEnd'), ('HitSeq','SbjSeq'), ('Score','BitScore')]
            bdb = self.db().copyTable(self.db('local'),'blast')
            for (ef,bf) in fieldswap:
                if ef in bdb.fields(): bdb.renameField(ef,bf)
            bdb.renameField('AlnID','ExAln')
            bdb.rankField('ExAln',newfield='AlnID',rev=False,absolute=True,lowest=True,unique=False)
            bdb.newKey(['Query','Hit','AlnID'])
            self.debug(bdb.dataKeys())
            blast = rje_blast.BLASTRun(self.log,['tuplekeys=T']+self.cmd_list)
            blast.obj['DB'] = self.db()
            ### ~ [1] Reduce Local ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add minloclen and minlocid?
            udb = blast.reduceLocal(bdb,minloclen=minloclen,minlocid=minlocid,byqry=byqry,sortfield='BitScore')
            for entry in udb.entries():
                entry['AlnID'] = entry['ExAln']
                self.bugPrint(entry)
            udb.remakeKeys()
            self.debug(udb.dataKeys())
            udb.dropField('ExAln')
            for (ef,bf) in fieldswap:
                if bf in udb.fields():
                    udb.renameField(bf,ef)
            return udb
        except: self.errorLog('%s.reduceLocal error' % self.prog()); raise
#########################################################################################################################
    def localToHitSum(self,locdb,sumname=''):  ### Compress a local hit table to a summary table for each alignment set.
        '''
        Compress a local hit table to a summary table for each alignment set.
        :param locdb: Local hits database Table
        :return: Hitsum database Table
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not sumname: sumname = '%s.%s' % (locdb.name(),'hitsum')
            vdb = self.db().copyTable(locdb,sumname)
            sfields = ['Qry','Hit','Rank','Score','AlnNum','Length','Identity','QryLen','QryStart','QryEnd','HitStart','HitEnd','HitStrand','Phase']
            if 'QryLen' not in locdb.fields(): sfields.remove('QryLen')
            ### ~ [1] Compress and Tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vdb.addField('AlnNum',evalue=1)
            for entry in vdb.entries():
                if entry['HitStart'] > entry['HitEnd']:
                    [entry['HitStart'], entry['HitEnd']] = [entry['HitEnd'],entry['HitStart']]
            rules = {'Score':'max','QryLen':'max','QryStart':'min','QryEnd':'max','HitStart':'min','HitEnd':'max','AlnID':'list','HitStrand':'list','Length':'sum','Identity':'sum','AlnNum':'sum'}
            vdb.compress(['Qry','Rank','Hit'],rules,'text')
            vdb.setFields(sfields)
            vdb.dataFormat({'Score':'int','QryLen':'int','QryStart':'int','QryEnd':'int','HitStart':'int','HitEnd':'int','Length':'int','Identity':'int'})
            #i# NOTE: Phase only makes sense for CDS searches. Proteins will be protein positions!!
            for ventry in vdb.entries(): ventry['Phase'] = (ventry['QryStart'] - 1) % 3
            return vdb
        except: self.errorLog('%s.localToHitSum error' % self.prog()); raise
#########################################################################################################################
    def gff3(self,locdb,sumdb,exmodel,filename=''):  ### Converts local and summary tables into GFF table.
        '''
        Converts local and summary tables into GFF table.
        :param locdb: Local hits database Table
        :param sumdb: Hit summary database Table
        :param exmodel: Exonerate model
        :param filename: Save to filename, if given.
        :return: GFF database Table
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not exmodel: exmodel = 'dna2dna'
            lochits = {}    # Dictionary of {hitsumkey:[local entries]}
            for hkey in sumdb.dataKeys(): lochits[hkey] = []
            for lkey in locdb.dataKeys():
                entry = locdb.data(lkey)
                hkey = (entry['Qry'],entry['Rank'],entry['Hit'])
                if hkey in lochits: lochits[hkey].append(entry)
                else: self.warnLog('Local hit (%s) without summary (%s)' % (lkey,hkey))

            gffdb = self.db().copyTable(locdb,'%sgff' % locdb.name()[:3])
            gffdb.addField('type',evalue='similarity')  #!# Add cds and gene types
            gffdb.addField('attributes',evalue='')  #!# Add cds and gene types

            e2gfields = [('Hit','seqid'), ('HitStart','start'), ('HitEnd','end'), ('Phase','phase'), ('Score','score'), ('HitStrand','strand')]
            gfields = string.split('seqid source type start end score strand phase attributes')

            # ID is unique
            # Name is not unique
            # Parent links exons or cds to gene
            geneid = {}

            for hkey in sumdb.dataKeys():
                entry = sumdb.data(hkey)
                entry['Phase'] = '.'
                entry['AlnID'] = entry['Rank']
                geneid[hkey] = '%s-%d' % (entry['Qry'],entry['Rank'])
                entry['attributes'] = 'ID=%s;Name=%s;' % (geneid[hkey],geneid[hkey])
                if 'QryDesc' in sumdb.fields(): entry['attributes'] += 'Note=Similar to %s;' % (entry['QryDesc'])
                if exmodel in ['est2genome','protein2genome']: entry['type'] = 'gene'
                else: entry['type'] = 'similarity'
                gffdb.addEntry(entry)

            for lkey in locdb.dataKeys():
                entry = gffdb.data(lkey)
                hkey = (entry['Qry'],entry['Rank'],entry['Hit'])
                if exmodel in ['est2genome','protein2genome']: entry['type'] = 'CDS'
                entry['attributes'] = 'ID=%s:%s;Parent=%s;Name=%s-%s;' % (geneid[hkey],entry['type'],geneid[hkey],geneid[hkey],entry['AlnID'])

            gffdb.addField('source',evalue='rje_apollo:exonerate:%s' % exmodel)

            qrytype = 'Qry'
            for entry in gffdb.entries():
                if entry['HitStart'] > entry['HitEnd']:
                    [entry['HitStart'], entry['HitEnd']] = [entry['HitEnd'],entry['HitStart']]
                entry['attributes'] += 'qstart=%s;' % min(entry['%sStart' % qrytype],entry['%sEnd' % qrytype])
                entry['attributes'] += 'qend=%s;' % max(entry['%sStart' % qrytype],entry['%sEnd' % qrytype])
                if 'QryLen' in entry and entry['QryLen']:
                    field = 'QryLen'
                    entry['attributes'] += 'qrylen=%s;' % (entry[field])
                entry['attributes'] += 'identity=%s;' % (entry['Identity'])
                entry['attributes'] += 'score=%s;' % (entry['Score'])
                entry['attributes'] = entry['attributes'][:-1]


            for (ef,gf) in e2gfields:
                if ef in gffdb.fields():
                    gffdb.renameField(ef,gf)
                else: gffdb.addField(gf,evalue='.')

            gffdb.newKey(['seqid','start','AlnID','end','type','source','attributes'])
            gffdb.keepFields(gfields+['AlnID'])
            ### ~ [2] Output GFF File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename:
                gcomments = ['##gff-version 3','#Generated by %s' % self.log.runDetails()]
                if not self.getStrLC('Rest'): gcomments.append('#Full Command List: %s' % rje.argString(rje.tidyArgs(self.log.cmd_list)))
                gffdb.saveToFile(filename,delimit='\t',backup=True,append=False,savefields=gfields,log=True,headers=False,comments=gcomments)

            return gffdb
        except: self.errorLog('%s.gff3 error' % self.prog()); raise
#########################################################################################################################
### End of SECTION II: Exonerate Class                                                                                  #
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
    try:#NewClass(mainlog,cmd_list).run()
        print rje_obj.zen(), '\n\n *** No standalone functionality! *** \n\n'

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
