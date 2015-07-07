#!/usr/bin/python

# See below for name and description
# Copyright (C) 2012 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SeqMapper
Description:  Sequence Mapping Program
Version:      2.1
Last Edit:    16/04/14
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for mapping one set of protein sequences onto a different sequence database, using Accession Numbers
    etc where possible and then using GABLAM when no direct match is possible. The program gives the following outputs:
    - *.*.mapped.fas = Fasta file of successfully mapped sequences
    - *.*.missing.fas = Fasta file of sequences that could not be mapped
    - *.*.mapping.tdt = Delimited file giving details of mapping (Seq, MapSeq, Method)
    If combine=T then the *.missing.fas file will not be created and unmapped sequences will be output in *.mapped.fas.
    Note that the possible mappings are all identified through BLAST and so a protein with matching IDs etc. but not
    hitting with BLAST will NOT be mapped. Currently only mapping of protein or nucleotides onto a protein database is
    supported.

    Unless the interactivity setting is set to 2 or more (i=2), sequences that are mapped using Name, AccNum, Sequence
    (100% identical sequences), ID or DescAcc will be mapped onto the first appropriate sequence. If automap > 0, then
    the best sequence according to the mapstat will be mapped automatically. If two sequences tie, the other two possible
    stats will also be used to rank the hits. If still tied and mapfocus is not "both" then the sequences will be ranked
    using both query and hit stats. If still tied, the first sequence will be selected.

    Any sequences that fall below automap (or i>1) but meet the minmap criteria will be ranked according to their BLAST
    rankings and then presented for a user decision. Presentation will be in reverse order, so that in the case of many
    possible mappings, the best options remain clear and on screen. The default choice (selected by hitting ENTER) will
    be the best ranked according to GABLAM stats, which will have been moved to position 1 if not already there. (BLAST
    rankings and GABLAM rankings will not always agree.)

    SeqMapper will enter a user menu if i>1 or seqin and/or mapdb are missing. If i=0 and one of these is missing, a
    simple prompt will ask for the missing files. If i<0 and one of these is missing, the program will exit. 

Commandline:
    ### Input Options ###
    seqin=FILE      : File of sequences to be mapped [None]
    mapdb=FILE      : File of sequences to map sequences onto [None]
    startfrom=X     : Shortname or AccNum of seqin file to startfrom (will append results) (memsaver=T only) [None]

    ### Output Options ###
    resfile=FILE    : Base of output filenames (*.mapped.fas, *.missing.fas & *.mapping.tdt)  [seqin.mapdb]
    combine=T/F     : Combine both fasta files in one (e.g. include unmapped sequences in *.mapping.fas) [False]
    gablamout=T/F   : Output GABLAM statistics for mapped sequences, including "straight" matches [True]
    append=T/F      : Append rather than overwrite results files [False]
    delimit=X       : Delimiter for *.mapping.* file (will set extension) [tab]

    ### Mapping Options ###
    mapspec=X       : Maps sequences onto given species code. "Self" = same species as query. "None" = any. [None]
    mapping=LIST    : Possible ways of mapping sequences (in pref order) [Name,AccNum,Sequence,ID,DescAcc,GABLAM,grep]
        - Name = First word of sequence name
        - Sequence = Identical sequence
        - grep = grep-based searching of sequence if no hits
        - ID = SwissProt style ID of GENE_SPECIES (note that the species may be changed according to mapspec)
        - AccNum = Primary Accession Number
        - DescAcc = Accession Number featured in description line in form "\WAccNum\W", where \W is non-
    skipgene=LIST   : List of "genes" in protein IDs to ignore [ens,nvl,ref,p,hyp,frag]
    mapstat=X       : GABLAM Stat to use for mapping assessment (if GABLAM in mapping list) (ID/Sim/Len) [ID]
    minmap=X        : Minimum value of mapstat for any mapping to occur [90.0]
    automap=X       : Minimum value of mapstat for automatic mapping to occur (if i<1) [99.5]
    ordered=T/F     : Whether to use GABLAMO rather than GABLAM stat [True]
    mapfocus=X      : Focus for mapping statistic, i.e. which sequence must meet requirements [query]
        - query = Best if query is ultimate focus and maximises closeness of mapped sequence)
        - hit = Best if lots of sequence fragments are in mapdb and should be allowed as mappings
        - either = Best if both above conditions are true
        - both = Gets most similar sequences in terms of length but can be quite strict where length errors exist

    ### Advanced BLAST Options ###
    blaste=X    : E-Value cut-off for BLAST searches (BLAST -e X) [1e-4]
    blastv=X    : Number of BLAST hits to return per query (BLAST -v X) [20]
    blastf=T/F  : Complexity Filter (BLAST -F X) [False]

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_menu, rje_obj, rje_seq, rje_seqlist, rje_zen
#########################################################################################################################
import rje_blast_V2 as rje_blast
import rje_sequence
#########################################################################################################################
### History
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Basic working version for protein databases.
    # 1.1 - Modified run() method to be called from other programs
    # 1.2 - Added grep method
    # 2.0 - Reworked with new Object format, new BLAST(+) module and new seqlist module.
    # 2.1 - Added catching of failure to read input sequences. Removed 'Run' from GABLAM table.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Update with new modules.
    # [ ] : Add DNA-DNA and DNA-Protein functionality?
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SeqMapper', '2.0', 'February 2013', '2006')
    description = 'Sequence Mapping Program'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
mapping_methods = ['name', 'accnum', 'sequence', 'id', 'descacc', 'gablam','grep']
method_info = {'name':'Name', 'accnum':'AccNum', 'sequence':'Sequence', 'id':'ID', 'descacc':'Description', 'gablam':'',
               'grep':'Grep Search'}
gstat_type = {'sim':'Sim','id':'ID','len':'Len'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SeqMapper Class                                                                                         #
#########################################################################################################################
class SeqMapper(rje_obj.RJE_Object):     
    '''
    SeqMapper Class. Author: Rich Edwards (2006).

    Str:str
    - SeqIn = File of sequences to be mapped [None]
    - MapDB = File of sequences to map sequences onto [None]
    - ResFile = Base of output filenames [seqin.mapdb]
    - MapSpec = Maps sequences onto given species code. "Self" = same species as query. "None" = any. [Self]
    - MapStat = GABLAM Stat to use for mapping assessment (if GABLAM in mapping list) (ID/Sim/Len) [ID]
    - MapFocus = Focus for mapping statistic, i.e. which sequence must meet requirements [query]
        - query = Best if query is ultimate focus and maximises closeness of mapped sequence)
        - hit = Best if lots of sequence fragments are in mapdb and should be allowed as mappings
        - either = Best if both above conditions are true
        - both = Gets most similar sequences in terms of length but can be quite strict where length errors exist
    - StartFrom = Shortname or AccNum of seqin file to startfrom (will append results) [None]
    - MapFas = mapped fasta file
    - MissFas = failed fasta file
    - MapRes = delimited file of mappings
    
    Bool:boolean
    - Ordered = Whether to use GABLAMO rather than GABLAM stat [True]
    - GablamOut = Output GABLAM statistics for mapped sequences, including "straight" matches [True]
    - Combine = Combine both fasta files in one (e.g. include unmapped sequences in *.mapping.fas) [False]

    Int:integer

    Num:float
    - MinMap = Minimum value of mapstat for any mapping to occur [90.0]
    - AutoMap = Minimum value of mapstat for automatic mapping to occur (if i<1) [99.5]
    
    List:list

    Dict:dictionary    
    - Mapping = Possible ways of mapping sequences (in pref order) [Name,AccNum,Sequence,ID,DescAcc,GABLAM,grep]
        - Name = First word of sequence name
        - ID = SwissProt style ID of GENE_SPECIES (note that the species may be changed according to mapspec)
        - AccNum = Primary Accession Number
        - DescAcc = Accession Number featured in description line in form "\WAccNum\W", where \W is non-
    - SkipGene = List of "genes" in protein IDs to ignore [ens,nvl,ref,p,hyp,frag]
    - Headers = list of headers for delimited file output
    - Mapped = List of mapped sequence IDs to check redundancy

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - MapDB = rje_seqlist.SeqList object of MapDB
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['SeqIn','MapDB','ResFile','MapSpec','MapStat','MapFocus','StartFrom','MapFas','MissFas','MapRes']
        self.boollist = ['Ordered','GablamOut','Combine']
        self.intlist = []
        self.numlist = ['MinMap','AutoMap']
        self.listlist = ['Mapping','SkipGene','Headers','Mapped']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=True,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'MapSpec':'None','MapStat':'ID','MapFocus':'Query'})
        self.setBool({'Combine':False})
        self.setInt({})
        self.setNum({'MinMap':90.0,'AutoMap':99.5})
        self.list['Mapping'] = mapping_methods[0:]
        self.list['SkipGene'] = string.split('ens,nvl,ref,p,hyp,frag',',')
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
                self._cmdReadList(cmd,'file',['SeqIn','MapDB','ResFile'])  
                self._cmdReadList(cmd,'str',['MapSpec','MapStat','MapFocus','StartFrom'])
                self._cmdReadList(cmd,'bool',['Ordered','GablamOut','Combine'])
                self._cmdReadList(cmd,'float',['MinMap','AutoMap'])
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Mapping'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.str['MapFocus'].lower() == 'qry': self.str['MapFocus'] = 'query'
        if self.str['StartFrom'].lower() == 'none': self.str['StartFrom'] = ''
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,imenu=False,outputmap=True,returndict=False):      ### Main controlling run Method
        '''
        Main controlling run Method.
        >> imenu:boolean = Whether to initiate interactive menu if appropriate [False].
        >> outputmap:boolean = Whether to output mapping into a file [True]
        >> returndict:boolean = Whether to return a dictionary of {searchname:mappedname} (no previous mapping) [False]
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(imenu): raise ValueError
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
            if not seqlist.seqNum(): self.warnLog('No sequences loaded for mapping.'); return {}
            ## ~ [0a] Setup BLAST Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = rje_blast.BLASTRun(self.log,['blaste=1e-4','blastv=20','blastf=F']+self.cmd_list+['v=-1'])
            blast.setStr({'DBase':self.getStr('MapDB'),'Type':'blastp','InFile':self.getStr('SeqIn'),
                         'Name':'%s-%s.blast' % (rje.baseFile(self.str['SeqIn'],True),rje.baseFile(self.str['MapDB'],True))})  
            blast.setStat({'HitAln':blast.getStat('OneLine')})
            blast.list['ResTab'] = ['Search','Hit','GABLAM']
            if seqlist.nt(): blast.str['Type'] = 'blastx'
            ## ~ [0b] Setup Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if outputmap: self._setupOutput()                           ## Output Files ##
            if returndict: mapdict = {}
            else: self._setupMapped()                                   ## Previously Mapped Sequences ##
            seqx = seqlist.seqNum()             ## Number of sequences ##
            ### ~ [1] BLAST Search Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#BLAST','BLASTing %s vs %s.\n *** This could take some time if files are large. Please be patient! ***' % (self.str['SeqIn'],self.str['MapDB']),log=False)
            ## ~ [1a] Perform BLAST Unless it exists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast.run(format=True)
            self.obj['DB'] = blast.obj['DB']
            ## ~ [1b] Mapping from searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.debug(self.getStr('MapDB'))
            self.obj['MapDB'] = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqmode=file','seqin=%s' % self.str['MapDB']])
            self.obj['MapDB'].loadSeq(self.getStr('MapDB'))
            self.debug('%s' % self.obj['MapDB'].list['Seq'])
            sx = 0
            while seqlist.nextSeq() != None:
                search = seqlist.getSeq(format='short')
                sx += 1
                ## Check StartFrom ##
                if self.str['StartFrom']:
                    if self.str['StartFrom'] != search:
                        self.progLog('\r#SKIP','Looking for %s: skipping %d seqs' % (self.str['StartFrom'],sx))
                        continue
                    self.str['StartFrom'] = ''
                    self.printLog('\r#SKIP','Starting from %s: skipped %d seqs' % (self.str['StartFrom'],sx))
                ## Check if in Mapped ##
                if search in self.list['Mapped']:
                    resdict = {'Query':search,'Hit':search,'Method':'Already Mapped!'}
                    self.printLog('#FAS','%s already in output - not duplicating in %s' % (search,self.str['MapFas']))
                    rje.delimitedFileOutput(self,self.str['MapRes'],self.list['Headers'],rje.getDelimit(self.cmd_list),resdict)
                    continue
                ### Map Sequence ###
                self.printLog('#MAP','Mapping %s seqs: %s of %s' % (self.str['SeqIn'],rje.integerString(sx),rje.integerString(seqx)))
                mapname = self.mapSeq(seqlist,blast,search)
                if returndict: mapdict[search] = mapname
            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#MAP','Mapping of %s (%s seqs) complete.' % (self.str['SeqIn'],rje.integerString(seqx)))           
            if os.path.exists(blast.str['Name']) and not (self.getBool('DeBug') or self.test()): os.unlink(blast.str['Name'])     #!# Add option to keep BLAST! #!#
            if returndict: return mapdict
        except: self.errorLog('Error in SeqMapper.run()',printerror=True,quitchoice=True); raise   
#########################################################################################################################
    def setup(self,imenu=False):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if imenu and self.i() > 0: self._setupMenu()     ## Menu ##
            needblast = False
            for method in self.list['Mapping'][0:]:         ## Mapping Methods ##
                if method.lower() not in mapping_methods:
                    self.errorLog('Mapping method "%s" not recognised!' % method,printerror=False)
                    self.list['Mapping'].remove(method)
            if not self.list['Mapping']:
                self.errorLog('No mapping methods selected!' % method,printerror=False)
                raise ValueError
            self.printLog('#MAP','Mapping: %s' % string.join(self.list['Mapping'],', '))
            self.checkInputFiles(['SeqIn','MapDB'],imenu)   ## Sequence Files ##
            return True
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def _setupMenu(self):   ### Interactive Menu for setup
        '''Interactive Menu for setup.'''
        try:### ~ [1] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seqfile in ['SeqIn','MapDB']: self.setStr({seqfile: rje.getFileName('%s Sequence File?' % seqfile,self.getStr(seqfile))})
            ### ~ [2] Setup Menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            menuhead = 'SeqMapper Options Menu'
            # (code,desc,type,key)
            menulist = [('','### INPUT OPTIONS ###','',''),
                        ('I','Sequence File to Map','info','SeqIn'),
                        ('DB','Mapping Database File','info','MapDB'),
                        ('S','StartFrom Name/AccNum','info','StartFrom'),
                        ('','\n### OUTPUT OPTIONS ###','',''),
                        ('R','Name of results file','info','ResFile'),
                        ('G','Output GABLAM Stats','opt','GablamOut'),
                        ('A','Append output','opt','Append'),
                        ('R','Name of results file','info','ResFile'),
                        ('D','Delimiter','info','Delimit'),
                        ('','\n### MAPPING OPTIONS ###','',''),
                        ('SP','Species for mapping (Code/Self/None)','info','MapSpec'),
                        ('M','Mapping option list','list','Mapping'),
                        ('K','SkipGene list','list','SkipGene'),
                        ('ST','Statistic for GABLAM mapping','info','MapStat'),
                        ('MM','Min value for mapping','stat','MinMap'),
                        ('AM','Min value for auto-mapping','stat','AutoMap'),
                        ('O','Use Ordered GABLAM','opt','Ordered'),
                        ('F','Focus of GABLAM mapping','info','MapFocus'),
                        ('','','return',True),
                        ('X','Exit menu and proceed.','return',True)]
            rje_menu.menu(self,menuhead,menulist,choicetext='Please select (Blank to proceed):',changecase=True)
            return
        except: self.errorLog('Noooooo! Error in SeqMapper._setupMenu()',quitchoice=True)
#########################################################################################################################
    def _setupOutput(self): ### Sets up output files self.str['MapFas','MissFas','MapRes']
        '''Sets up output files self.str['MapFas','MissFas','MapRes'].'''
        ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        delimit = rje.getDelimit(self.cmd_list)
        if self.str['StartFrom'].lower() in ['','none']: self.str['StartFrom'] = ''
        else:
            self.bool['Append'] = True
            self.printLog('#CMD','StartFrom = "%s" so Append=T' % self.str['StartFrom'])
        ### ~ [1] General ResFile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        files = {'MapFas':'mapping.fas','MissFas':'missing.fas','MapRes':'mapping.%s' % rje.delimitExt(delimit)}
        if self.getBool('Combine'): files.pop('MissFas')
        if self.str['ResFile'].lower() in ['','none']:
            self.str['ResFile'] = '%s.%s' % (rje.baseFile(self.str['SeqIn']),rje.baseFile(self.str['MapDB'],strip_path=True))
        for file in files.keys():
            self.setStr({file: self.getStr('ResFile') + '.' + files[file]})
            rje.backup(self,self.getStr(file))
        ### ~ [2] Headers for MapRes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #!# Consider replacing with rje_db object? #!#
        self.list['Headers'] = ['Query','Hit','Method','MapRank','BlastRank','EVal','Score']
        for qh in ['Query','Hit']:
            self.list['Headers'] += ['%s_Species' % qh]
            if self.bool['GablamOut']:
                for st in ['Len','Sim','ID']:
                    self.list['Headers'] += ['%s_%s' % (qh,st)]
        rje.delimitedFileOutput(self,self.str['MapRes'],self.list['Headers'],delimit)
#########################################################################################################################
    def _setupMapped(self):     ### Sets up list of Previously Mapped Sequences
        '''Sets up list of Previously Mapped Sequences.'''
        ### Setup ###
        self.list['Mapped'] = []    # List of mapped sequence names
        if not self.bool['Append'] or not os.path.exists(self.str['MapFas']): return
        ### Previous Sequences ###
        seqlist = rje_seq.SeqList(None,['i=-1','v=-1','autoload=F','seqin=%s' % self.str['MapFas']])
        SEQFILE = open(filename,'r')
        lastline = ''
        sx = 0
        ### Count ###
        while 1:
            (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
            seqlist.seq = []
            if nextseq:
                sx += 1
                self.list['Mapped'].append(nextseq.shortName())
            else:
                break
        SEQFILE.close()
        self.printLog('#MAP','Read names of %s previously mapped sequences for redundancy checking' % rje.integerString(sx))
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def mapSeq(self,seqlist,blast,search,outputmap=True): ### Performs actual mapping of sequence
        '''
        Performs actual mapping of sequence.
        >> seq:SeqList object containing Sequence Object to be mapped
        >> blast:BLAST_Run object to perform BLAST and GABLAM
        >> search:Current BLAST search object for mapping
        >> outputmap:boolean = Whether to output mapping into a file [True]
        << returns shortName() of mapped sequence (or None if none)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seq = seqlist.getSeq(format='tuple')
            mapseq = self.obj['MapDB']
            hits = blast.db('Hit').indexEntries('Query',search)
            self.printLog('#HITS','%s vs %s = %d hits' % (search,blast.str['DBase'],len(hits)))
            hitseq = {}; hitdata = {}
            for entry in hits:
                hitseq[entry['Hit']] = mapseq.getDictSeq(entry['Hit'],format='tuple')
                hitdata[entry['Hit']] = entry
            resdict = {'Query':search,'Hit':None,'Method':'Failed','Query_Species':rje_sequence.specCodeFromName(seq[0])}
            ### ~ [1] Order Hits and Check Species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hits,hitdict) = self.orderHits(seq,hits,hitseq)
            self.debug(hits)
            self.debug(hitdict)
            ### ~ [2] Attempt mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for method in self.list['Mapping']:
                resdict['Hit'] = self.mapHit(seq,hits,hitdict,method.lower())
                if resdict['Hit']:
                    resdict['Method'] = method[:1].upper() + method[1:].lower()
                    break
                elif method == 'gablam' and (len(hits) > 0):
                    resdict['Method'] = 'Rejected'
            self.debug(resdict)
            ### ~[3] Output! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if resdict['Hit']:  #hitdict[hit]['Data']['ShortName']
                hit = resdict['Hit']['Hit']     # resdict['Hit'] is the BLAST table entry for Hit
                shortname = hitdict[hit]['Data']['ShortName']   # This is just hit!
                self.printLog('#MAP','%s mapped to %s (by %s)' % (string.split(seq[0])[0],shortname,resdict['Method']))
                ## Update Stats ##
                self.debug('')
                resdict['BlastRank'] = hitdata[hit]['Rank']
                for key in hitdict[hit]: resdict[key] = hitdict[hit][key]
                ## Fasta and Redundancy ##
                if shortname in self.list['Mapped']: self.printLog('#MAP','%s already mapped before - not duplicating in %s' % (shortname,self.str['MapFas']))
                else:
                    self.list['Mapped'].append(shortname)
                    if outputmap:
                        open(self.str['MapFas'],'a').write('>%s\n%s\n' % (hitseq[hit][0],hitseq[hit][1]))
                resdict['Hit_Species'] = hitdict[hit]['Data']['SpecCode']
                resdict['Hit'] = shortname
            else:
                ### ~ [2] GREP-based search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if 'grep' in self.list['Mapping']:
                    greplist = []; hitseq = ''
                    self.printLog('#GREP','grep %s %s -B 1' % (seq[1],blast.str['DBase']),log=False)
                    for line in os.popen('grep %s %s -B 1' % (seq[1],blast.str['DBase'])).readlines():
                        if line[:1] == '>': greplist.append(string.split(line[1:])[0])
                        elif not hitseq: hitseq = rje.chomp(line)
                    if greplist:
                        shortname = greplist.pop(0)
                        resdict['Hit'] = shortname
                        resdict['Method'] = 'Grep'
                        resdict['Qry_ID'] = '100.0'
                        resdict['Qry_Len'] = len(seq[1])
                        resdict['Hit_Len'] = len(hitseq)
                        resdict['Hit_ID'] = 100.0 * len(hitseq) / len(seq[1])
                        try: resdict['Hit_Species'] = string.split(shortname,'_')[1]
                        except: pass
                        if shortname in self.list['Mapped']:
                            self.printLog('#MAP','%s already mapped before - not duplicating in %s' % (shortname,self.str['MapFas']))
                        else:
                            self.list['Mapped'].append(shortname)
                            if outputmap: open(self.str['MapFas'],'a').write('>%s\n%s\n' % (shortname,hitseq))
                    for extra in greplist: self.printLog('#GREP','Warning! Query "%s" also hit "%s" with grep!' % (string.split(seq[0])[0],extra))
                if not resdict['Hit'] and self.bool['Combine']:
                    ## Fasta and Redundancy ##
                    shortname = string.split(seq[0])[0]
                    if shortname in self.list['Mapped']:
                        self.printLog('#FAS','%s already in output - not duplicating in %s' % (shortname,self.str['MapFas']))
                    else:
                        self.list['Mapped'].append(shortname)
                        if outputmap:
                            open(self.str['MapFas'],'a').write('>%s\n%s\n' % (seq[0],seq[1]))
                elif outputmap:
                    open(self.str['MissFas'],'a').write('>%s\n%s\n' % (seq[0],seq[1]))
                self.printLog('#MISS','%s mapping %s' % (resdict['Query'],resdict['Method']))
            if outputmap:
                rje.delimitedFileOutput(self,self.str['MapRes'],self.list['Headers'],rje.getDelimit(self.cmd_list),resdict)
            return resdict['Hit']

        except:
            self.errorLog('Fudgesticks! SeqMapper.mapSeq(%s) has died!' % seq[0],quitchoice=True)
            return False
#########################################################################################################################
    def orderHits(self,seq,hits,hitseq):    ### Returns orderd hits and fuller hitdict with GABLAM Stats etc.
        '''
        Returns orderd hits and fuller hitdict with GABLAM Stats etc.
        >> seq:Query Sequence Object (name,sequence)
        >> hits:List of hit entries 
        >> hitseq:Dictionary of {Hit:Sequence}
        << Tuple of (hits,hitdict)
        '''
        try:### ~ [0] Order hits by BLAST Rank ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sorted = hits[0:]
            for hit in hits: sorted[hit['Rank']-1] = hit            
            ### ~ [1] Screen Species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hits = sorted
            hits.reverse()
            mapspec = None
            qryspec = rje_sequence.specCodeFromName(seq[0])
            if self.getStr('MapSpec').lower() == 'self': mapspec = qryspec
            elif self.getStr('MapSpec').lower() not in ['any','none']: mapspec = self.getStr('MapSpec').upper()
            for hit in hits[0:]:
                #self.debug(hit)
                #self.debug('%s vs %s' % (rje_sequence.specCodeFromName(hitseq[hit['Hit']][0]), mapspec))
                if mapspec and rje_sequence.specCodeFromName(hitseq[hit['Hit']][0]) != mapspec:
                    hits.remove(hit)
                    hitseq.pop(hit['Hit'])
                elif not mapspec and rje_sequence.specCodeFromName(hitseq[hit['Hit']][0]) == qryspec:
                    hits.remove(hit)
                    hits.append(hit)    # => Move hits from same species to top of list
            hits.reverse()
            ### ~ [2] Convert GABLAM Stats to Percentages and move to hitdict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.debug(self.db())
            gstat = 'GABLAMO'
            if not self.bool['Ordered']: gstat = 'GABLAM'
            hitdict = {}
            for hit in hits[0:]:
                hitname = hit['Hit']
                hitdict[hitname] = {'Seq':hitseq[hitname],'MapRank':hits.index(hit)+1,'EVal':hit['E-Value'],'Score':hit['BitScore']}
                hitgab = self.db('GABLAM').data(self.db('GABLAM').makeKey({'Query':hit['Query'],'Hit':hitname,'QryHit':'Hit'}))
                qrygab = self.db('GABLAM').data(self.db('GABLAM').makeKey({'Query':hit['Query'],'Hit':hitname,'QryHit':'Query'}))
                for st in ['ID','Len','Sim']:
                    hitdict[hitname]['Query_%s' % st] = 100 * float(qrygab['%s %s' % (gstat,st)]) / len(seq[1])
                    hitdict[hitname]['Hit_%s' % st] = 100 * float(hitgab['%s %s' % (gstat,st)]) / hit['Length']
                    hitdict[hitname]['Both_%s' % st] = (hitdict[hitname]['Query_%s' % st] + hitdict[hitname]['Hit_%s' % st]) / 2.0
                    if hitdict[hitname]['Query_%s' % st] >= hitdict[hitname]['Hit_%s' % st]:
                        hitdict[hitname]['Either_%s' % st] = hitdict[hitname]['Query_%s' % st]
                    else:
                        hitdict[hitname]['Either_%s' % st] = hitdict[hitname]['Hit_%s' % st]
            ### ~ [3] Rank Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            focus = self.str['MapFocus'][:1].upper() + self.str['MapFocus'][1:].lower()
            scorelist = []
            multorder = ['Both_Len','Both_Sim','Both_ID','%s_Len' % focus,'%s_Sim' % focus,'%s_ID' % focus]
            for hit in hits:
                hitname = hit['Hit']
                s = 0
                m = 1           
                for mult in multorder:
                    s += (hitdict[hitname][mult] * m)
                    m *= 1000
                scorelist.append(s)
            ranks = rje.rankList(scorelist,rev=True,absolute=True,lowest=True)  ### Returns rank of scores as list
            rankdict = {}
            for h in range(len(hits)):
                hit = hits[h]
                hitname = hit['Hit']
                r = ranks[h]
                hitdict[hitname]['HitRank'] = r
                if rankdict.has_key(r): rankdict[r].append(hit)
                else: rankdict[r] = [hit]
            newlist = []
            for r in rje.sortKeys(rankdict): newlist += rankdict[r]
            return (newlist,hitdict)
        except: self.errorLog('Bad goings on in SeqMapper.orderHits()',quitchoice=True)
#########################################################################################################################
    def mapHit(self,seq,hits,hitdict,method):     ### Tries to map seq onto hitseq and returns hit if successful
        '''
        Tries to map seq onto hitseq and returns hit if successful.
        >> seq:Query Sequence Object
        >> hits:List of hits in rough order of goodness
        >> hitdict:Dictionary of {hitname:stats}
        >> method:Mapping method to use
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (name,sequence) = seq
            data = rje_sequence.extractNameDetails(name,self)
            data['Sequence'] = seq[1]
            data['ShortName'] = string.split(seq[0])[0]
            for hit in hitdict:
                hitdict[hit]['Data'] = rje_sequence.extractNameDetails(hitdict[hit]['Seq'][0],self)
                hitdict[hit]['Data']['Sequence'] = hitdict[hit]['Seq'][1]
                hitdict[hit]['Data']['ShortName'] = string.split(hitdict[hit]['Seq'][0])[0]
            ### SkipGene ###
            if method == 'id' and rje.matchExp('^(\S+)_\S+',data['ID']):
                gene = rje.matchExp('^(\S+)_\S+',data['ID'])
                if gene in self.list['SkipGene']:
                    return None
            ### Name, AccNum, Sequence and ID ###
            if method_info[method] in ['Name', 'AccNum', 'Sequence', 'ID']:
                for hit in hits:
                    hitdata = hitdict[hit['Hit']]['Data']
                    if hitdata[method_info[method]] == data[method_info[method]]:
                        if self.i() < 2 or rje.yesNo('Map %s to %s?' % (data['ShortName'],hitdata['ShortName'])):
                            return hit
            ### DescAcc ###
            if method == 'descacc':
                for hit in hits:
                    hitdata = hitdict[hit['Hit']]['Data']
                    if rje.matchExp('\W(%s)\W' % data['AccNum'],hitdata['Name']):
                        if self.i() < 2 or rje.yesNo('Map %s to %s?' % (data['ShortName'],hitdata['ShortName'])):
                            return hit
            ### GABLAM ###
            if method != 'gablam': return None
            focus = self.str['MapFocus'][:1].upper() + self.str['MapFocus'][1:].lower()
            gstat = gstat_type[self.str['MapStat'].lower()]
            possibles = []  # List of Hits that meet MinMap criterion
            for hit in hits:
                hitname = hit['Hit']
                hitdata = hitdict[hit['Hit']]['Data']
                if self.getNum('AutoMap') > 0.0 and hitdict[hitname]['%s_%s' % (focus,gstat)] >= self.getNum('AutoMap'):
                    if self.i() < 2 or rje.yesNo('Map %s to %s?' % (data['ShortName'],hitdata['ShortName'])):
                        return hit
                elif hitdict[hitname]['%s_%s' % (focus,gstat)] >= self.getNum('MinMap'):
                    possibles.append(hit)
            ### Manual GABLAM Choice ###
            if self.i() < 0 or not possibles: return None
            possibles.reverse()
            print '\nMapping options for %s:\n' % data['ShortName']
            for p in range(len(possibles)):
                hit = possibles[p]
                hitname = hit['Hit']
                hitdata = hitdict[hit['Hit']]['Data']
                print '<%d> %s (%d aa) =\t' % (len(possibles)-p,hitdata['Name'],hit['Length']),
                print '%.1f%% Qry Len,' % (100.0 * hit['Length'] / len(seq[1])),
                print '%.1f%% ID (%.1f%% Sim, %.1f%% Cov.)' % (hitdict[hitname]['Hit_ID'],hitdict[hitname]['Hit_Sim'],hitdict[hitname]['Hit_Len']),
                print '(Qry: %.1f%% ID (%.1f%% Sim, %.1f%% Cov.)' % (hitdict[hitname]['Query_ID'],hitdict[hitname]['Query_Sim'],hitdict[hitname]['Query_Len'])
            choice = -1
            print '<0> No mapping.\n'
            ## Choice ##
            while 1:
                choice = rje.getInt('Select sequence to replace %s?' % data['ShortName'],default=1,confirm=True)
                i = len(possibles) - choice
                if choice == 0: # No mapping
                    if self.i() < 2 or rje.yesNo('No GABLAM mapping for %s?' % (data['ShortName'])): return None
                elif choice > 0 and choice <= len(possibles):    
                    hit = possibles[i]
                    hitdata = hitdict[hit['Hit']]['Data']
                    if self.i() < 2 or rje.yesNo('Map %s to %s?' % (data['ShortName'],hitdata['ShortName'])): return hit
        except:
            self.errorLog('Problem during SeqMapper.mapHit(%s)' % method,quitchoice=True)
            return None
#########################################################################################################################
### End of SECTION II: SeqMapper Class                                                                                  #
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
    try: SeqMapper(mainlog,['i=1','v=1']+cmd_list).run(imenu=True)

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
