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
Module:       rje_apollo
Description:  WebApollo genome search program
Version:      0.6.1
Last Edit:    02/07/18
Webserver:    http://www.slimsuite.unsw.edu.au/servers/apollo.php
Copyright (C) 2018  Richard J. Edwards - See source code for GNU License Notice

Function:
    RJE_APOLLO is a wrapper for GABLAM that peforms a limited number of focused GABLAM functions and adds some additional
    data cleanup and outputs. It is primarily designed for use as part of the SLiMSuite REST servers as a way to search
    the WebApollo genome servers of the Edwards Lab. See the REST output and Apollo webserver for more details.

    Version 0.2.0 added the option of running exonerate rather than BLAST+, which will probably become the default option
    on the server.

Commandline:
    ### ~ Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Fasta file of input sequences to search against genome []
    genome=FASFILE  : Fasta file of genome to search against
    qrytype=X       : Type of seqin query file (cds/transcript/dna/protein) [cds]
    basefile=FILE   : Root of output files [apollo]
    maxqry=INT      : Maximum number of sequences allowed for seqin [0]

    ### ~ GABLAM Search Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    blaste=NUM      : E-Value cut-off for BLAST searches (BLAST -e X) [1e-10]
    tophits=INT     : Sets max number of BLAST hits returned (blastb and blastv) [10]
    localmin=INT    : Minimum length of local alignment to output to local stats table [30]
    localidmin=PERC : Minimum local %identity of local alignment to output to local stats table [75.0]
    keepblast=T/F   : Whether to keep the blast results files rather than delete them [True]
    blasttask=X     : Flavour of blast to use (BLAST -task X) (NOTE: blastn default is megablast) [blast+ default]
    exonerate=T/F   : Whether to use exonerate in place of GABLAM. Can also set blasttask=exonerate [False]
    exopt=X         : Additional exonerate options []

    ### ~ System Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    webapollo=URL   : URL base for Web Apollo genome browser [https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/]
    genomeid=X      : Genome identifier code for Web Apollo genome browser (with optional .public), or SPEC:ID list []
    tracks=LIST     : Comma separated list of tracks to display for WebApollo server []
    exoneratebin=PATH : Path to exonerate program bin (can be '' if added to environment path) ['']

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
import rje, rje_db, rje_exonerate, rje_obj, rje_html, rje_seqlist
import tga_exonerate
import rje_blast_V2 as rje_blast
import gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial Working Version.
    # 0.1.1 - Modified HTML output.
    # 0.2.0 - Added basic exonerate output.
    # 0.3.0 - Added transcript query type.
    # 0.4.0 - Added MaxQry=INT : Maximum number of sequences allowed for seqin [0]
    # 0.5.0 - Added QryDesc to exonerate GFF output.
    # 0.6.0 - Added cdict parsing to genomeid.
    # 0.6.1 - Debugging genomeID dict and fixed cdict parsing bug.
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
    # [ ] : Consider adding samtools conversion to BAM format.
    # [ ] : Add proper GFF output
    # [ ] : Add QBLAST equivalent output from unique hits table.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_APOLLO', '0.6.1', 'July 2018', '2018')
    description = 'WebApollo genome search program'
    author = 'Dr Richard J. Edwards & Timothy G. Amos'
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
### SECTION II: Apollo Class                                                                                               #
#########################################################################################################################
class Apollo(rje_obj.RJE_Object):
    '''
    Apollo Class. Author: Rich Edwards (2015).

    Str:str
    - BlastTask=X     : Flavour of blast to use (BLAST -task X) (NOTE: blastn default is megablast) [blast+ default]
    - ExonerateBin=PATH : Path to exonerate program bin (can be '' if added to environment path) ['']
    - SeqIn=FASFILE   : Fasta file of input sequences to search against genome []
    - Genome=FASFILE  : Fasta file of genome to search against
    - GenomeID=X      : Genome identifier code for Web Apollo genome browser []
    - QryType=X       : Type of seqin query file (cds/transcript/dna/protein) [cds]
    - Tracks          : Comma separated list of tracks to display for WebApollo server []
    - WebApollo=URL   : URL base for Web Apollo genome browser [https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/]

    Bool:boolean
    - Exonerate=T/F   : Whether to use exonerate in place of GABLAM [False]

    Int:integer
    - MaxQry=INT : Maximum number of sequences allowed for seqin [0]
    - TopHits=INT     : Sets max number of BLAST hits returned (blastb and blastv) [10]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - Tracks=LIST     : Comma separated list of tracks to display for WebApollo server []

    Dict:dictionary
    - GenomeID  : Dictionary of {Species code: GenomeID}. '*' will be used for any unassigned species codes.

    Obj:RJE_Objects
    - BLAST = rje_blast.BLASTRun object that performs qassemble run for single queries
    - GABLAM = gablam.GABLAM object that performs most of the actual processing.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BlastTask','ExonerateBin','Genome','GenomeID','QryType','SeqIn','Tracks','WebApollo']
        self.boollist = ['Exonerate']
        self.intlist = ['MaxQry','TopHits']
        self.numlist = []
        self.filelist = []
        self.listlist = ['Tracks']
        self.dictlist = ['GenomeID']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'QryType':'cds','WebApollo':'https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/',
                     'ExonerateBin':''})
        self.setBool({'Exonerate':False})
        self.setInt({'TopHits':10})
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
                self._cmdReadList(cmd,'str',['BlastTask','GenomeID','QryType','WebApollo'])   # Normal strings
                self._cmdReadList(cmd,'path',['ExonerateBin'])  # String representing directory path
                self._cmdReadList(cmd,'file',['Genome','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['Exonerate'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxQry','TopHits'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Tracks'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['GenomeID']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if not self.dict['GenomeID']: self.dict['GenomeID'] = {'*':self.getStr('GenomeID')}
        self.printLog('#GENID','Apollo genome IDs: %s' % self.dict['GenomeID'])
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Exonerate'):
                self.exonerate()
            else:
                #i# Run GABLAM to generate main outputs and cleanup
                gab = self.GABLAM()
                if not gab: raise ValueError('GABLAM failed')
                #i# Run QAssemble on single sequence
                seqlist = gab.obj['SeqList']
                if seqlist.seqNum() == 1: self.QBLAST()

            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Rest'): self.restSetup()
            self.setStr({'Tracks':string.join(self.list['Tracks'],',')})
            ### ~ [2] Exonerate run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('BlastTask') in ['exonerate','est2genome','protein2genome']:
                self.setBool({'Exonerate':True})
                self.printLog('#BLAST','blasttask=%s -> switched exonerate=T.' % self.getStrLC('BlastTask'))
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

        ### General REST Outputs:
        html = Reduced summary table with links to WebApollo genome
        seqin = Input search sequence [fasta]
        hitsum = One line summary of hits per query [tdt]
        local = Local alignments table [tdt]
        locgff = Local alignments table in GFF3 format [gff3]
        unique = Local alignments table, reduced to non-overlapping query regions [tdt]

        ### Exonerate REST Outputs:
        alignment = Exonerate alignment [txt]
        exgff = Exonerate GFF2 output [gff2]
        cigar = Exonerate cigar output [txt]
        dnahits = Fasta file of hit sequences [fasta]
        prothits = Fasta file of translated hits (CDS or protein only) [fasta]

        ### BLAST+ REST Outputs:
        gablam = Full GABLAM summary table of each Query-Hit pair [tdt]
        blast = Raw BLAST results [txt]
        unigff = Local BLAST+ alignments table, reduced to non-overlapping query regions, in GFF3 format [gff]
        sam = Local BLAST alignments table in SAM alignment format (blastn only) [sam]
        qblast = BLAST results formatted in query assembly `-outfmt 4` [txt]
        qaln = Query sequence aligned to consensus for each hit from QBLAST [fasta]
        qhits = Degapped consensus hits from QBLAST [fasta]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        outfmt = ['html','seqin','hitsum','local','locgff','unique']
        if self.getBool('Exonerate'):
            outfmt += ['alignment','exgff','cigar','dnahits','prothits']
        else:
            outfmt += ['gablam','blast','unigff','sam','qblast','qaln','qhits']
        return outfmt
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
# Public:
# https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/120/jbrowse/index.html?loc=chr_WON710A1__DCMF160804:1232840..1233371&tracks=
# Logged in:
# https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/annotator/loadLink?loc=chr_WON710A1__DCMF160804:1232840..1233371&organism=120&tracks=

#     python ~/Dropbox/_Repository_/slimsuite/tools/seqsuite.py blast -blasti tlr1.fas -blastd ../../data/2017-08-18.Genome/Nscutatus.k81L500.2017-08-18.fasta blastp=tblastn -blasto Nscutatus.k81L500.tlr1.qblast -qassemblefas -qcomplete -basefile Nscutatus.k81L500.tlr1 blasta=8 qconsensus=hit
#

#########################################################################################################################
    def exonerate(self):   ### Runs Exonerate to generate main outputs.
        '''
        Runs Exonerate to generate main outputs.

        Outputs:
        *.exonerate.txt = replaces BLAST output
        *.hitsum.tdt = hitsum = One line summary of hits per query [tdt]
        *.local.tdt = local = Local BLAST alignments table [tdt]
        *.local.gff = locgff = Local BLAST alignments table in GFF3 format [gff3]
        *.local.sam = sam = Local BLAST alignments table in SAM alignment format (blastn only) [sam]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            exobj = rje_exonerate.Exonerate(self.log,self.cmd_list)
            exmodel = ''     # This is for DNA versus DNA.
            if self.getStrLC('BlastTask') in ['est2genome','protein2genome']:
                exmodel = self.getStrLC('BlastTask')
            elif self.getStrLC('QryType') in ['cds','transcript','trans']: exmodel = 'est2genome'
            elif self.getStrLC('QryType') == 'prot': exmodel = 'protein2genome'
            elif self.getStrLC('QryType') != 'dna': raise ValueError('Unrecognised qrytype: %s' % self.getStrLC('QryType'))
            if self.getInt('MaxQry') > 0:
                seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
                #!# Add max length and max combined length too?
                if self.getInt('MaxQry') < seqlist.seqNum():
                    if self.getStr('Rest'): raise ValueError('Server Apollo query limit (%d sequences) exceeded. Please reduce and rerun.' % self.getInt('MaxQry'))
                    else: raise ValueError('Apollo query limit (maxqry=%d sequences) exceeded. Please reduce and rerun.' % self.getInt('MaxQry'))

            ### ~ [2] Run Exonerate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Output']['seqin'] = self.getStr('SeqIn')
            self.printLog('#RUN','Running exonerate %s' % exmodel)
            #exodata = tga_exonerate.exonerate(self.getStr('SeqIn'), self.getStr('Genome'), exmodel, '%sexonerate' % self.getStr('ExonerateBin'),self.getInt('TopHits'))
            exodata = exobj.exonerate(self.getStr('SeqIn'), self.getStr('Genome'), exmodel, '%sexonerate' % self.getStr('ExonerateBin'),self.getInt('TopHits'))
            if not exodata: raise ValueError('Exonerate run failed!')
            #self.debug(exodata)
            exout = {'gff':[],'cigar':[],'alignment':[],'vulgar':[]}
            excmd = ['exonerate --model', exmodel, '--showtargetgff', '--showcigar']
            if self.getInt('TopHits'): excmd += ['--bestn', '%d' % self.getInt('TopHits')]
            exout['alignment'] = [self.log.runDetails(),string.join(excmd),'']

            ## ~ [2a] Process output text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=tuple'])
            queries = []
            seqlen = {}
            for seq in seqlist.seqs(): seqlen[seqlist.shortName(seq)] = seqlist.seqLen(seq)
            for qry in seqlist.names():
                if qry not in exodata: self.warnLog('Could not find query %s in exonerate data' % qry)
                else: queries.append(qry)
            for qry in rje.sortKeys(exodata):
                if qry not in queries:
                    self.warnLog('Could not find exonerate query %s in seqin sequences' % qry)
                    queries.append(qry)
            vfields = ['query_id', 'query_start', 'query_end', 'query_strand', 'target_id', 'target_start', 'target_end', 'target_strand', 'score', '<label, query_length, target_length> triplets','Rank','AlnID','alnx']
            vdb = self.db().addEmptyTable('local',vfields,['query_id', 'AlnID', 'target_id'])
            alnid = {}
            #self.bugPrint(queries)
            for query in queries:
                hits = []
                alnx = 0
                exdata = exodata[query]
                ## GFF output - self.dict['Output']['gff']
                if 'gff' in exdata:
                    if exout['gff']:
                        while exdata['gff'] and exdata['gff'][0].startswith('#'): exdata['gff'].pop(0)
                    commenting = True
                    while commenting:
                        commenting = exdata['gff'][0].startswith('#')
                        exout['gff'].append(exdata['gff'].pop(0))
                    while exdata['gff']:
                        if exdata['gff'][0].startswith('#'): exdata['gff'].pop(0)
                        else: exout['gff'].append(exdata['gff'].pop(0))
                else: self.warnLog('No GFF data for %s' % query)

                ## Cigar strings - self.dict['Output']['sam']
                #i# eogA_NOTSC__TS10XV2PRIEOG0907027J 0 2241 + phap_NOTSC__TS10XV2PHAP49 1692572 1640045 - 11103  M 83 D 14119 M 851 D 5363 M 37 D 460 M 35 D 2995 M 188 D 8449 M 43 D 1365 M 72 D 937 M 199 D 8990 M 94 D 421 M 113 D 3443 M 77 D 3067 M 230 D 677 M 219
                #!# Convert to SAM format
                if 'cigar' in exdata:
                    exout['cigar'] += exdata['cigar']
                else: self.warnLog('No cigar data for %s' % query)

                ## Alignment - self.dict['Output']['blast']
                if 'alignment' in exdata:
                    exout['alignment'] += exdata['alignment']
                else: self.warnLog('No alignment data for %s' % query)

                ## Vulgar - self.dict['Output']['local'] and self.dict['Output']['html']
                #self.debug('%s -> %s' % (query,exdata))
                if 'vulgar' in exdata:
                    #self.debug(exdata['vulgar'])
                    for vdata in exdata['vulgar'][1:]:
                        if vdata['target_id'] not in hits: hits.append(vdata['target_id'])
                        vdata['Rank'] = hits.index(vdata['target_id']) + 1
                        qh = (vdata['query_id'],vdata['target_id'])
                        if qh not in alnid: alnid[qh] = 1
                        else: alnid[qh] += 1
                        alnx += 1
                        vdata['alnx'] = alnx
                        vdata['AlnID'] = alnid[qh]
                        vdb.addEntry(vdata)
                        #self.bugPrint('\nAlnX: %s' % alnid)
                        #self.debug('%s -> %d' % (vdata,vdb.entryNum()))
                else: self.warnLog('No vulgar data for %s' % query)

            ## ~ [2b] Cleanup and save outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for (etype,otype) in [('cigar','cigar'),('alignment','alignment')]:
                self.dict['Output'][otype] = '%s.%s.txt' % (self.baseFile(),etype)
                open(self.dict['Output'][otype],'w').write(string.join(exout[etype],'\n'))

            self.dict['Output']['exgff'] = '%s.gff2' % (self.baseFile())
            open(self.dict['Output']['exgff'],'w').write(string.join(exout['gff'],'\n'))

            #i# Replace with parsed alignment if possible
            ldb = None
            if exout['alignment']:
                ldb = exobj.alignmentToLocal(exout['alignment'],protqry=self.getStrLC('QryType') == 'prot')
            if ldb:
                lfields = ['Qry','Rank','Hit','AlnID','Score','Length','Identity','QryLen','QryStart','QryEnd','HitStart','HitEnd','HitStrand']
                if self.getStrLC('QryType') in ['prot','cds']:
                    lfields.append('Phase')
                    exobj.obj['ProtHits'].saveSeq(seqfile='%s.prot.fas' % self.baseFile())
                    self.dict['Output']['prothits'] = '%s.prot.fas' % self.baseFile()
                exobj.obj['DNAHits'].saveSeq(seqfile='%s.dna.fas' % self.baseFile())
                self.dict['Output']['dnahits'] = '%s.dna.fas' % self.baseFile()
                ldb.addField('QryLen',evalue=0)
                for entry in ldb.entries(): entry['QryLen'] = seqlen[entry['Qry']]
                vdb = ldb
                self.debug(ldb.dataKeys())
                udb = exobj.db('unique')
                udb.addField('QryLen',evalue=0)
                for entry in udb.entries(): entry['QryLen'] = seqlen[entry['Qry']]
                self.debug(udb.dataKeys())
                udb.saveToFile('%s.unique.tdt' % (self.baseFile()),savefields=lfields)
                self.dict['Output']['unique'] = '%s.unique.tdt' % (self.baseFile())
            else:
                vdb.renameField('query_id','Qry')
                vdb.renameField('query_start','QryStart')
                vdb.renameField('query_end','QryEnd')
                vdb.renameField('query_strand','QryStrand')
                vdb.renameField('target_id','Hit')
                vdb.renameField('target_start','HitStart')
                vdb.renameField('target_end','HitEnd')
                vdb.renameField('target_strand','HitStrand')
                vdb.renameField('score','Score')
                vdb.remakeKeys()
                lfields = ['Qry','Hit','AlnID','Score','QryStart','QryEnd','QryStrand','HitStart','HitEnd','HitStrand']
            vdb.saveToFile('%s.local.tdt' % (self.baseFile()),savefields=lfields)
            self.dict['Output']['local'] = '%s.local.tdt' % (self.baseFile())

            #!# Compress
            if ldb:
                sfields = ['Qry','Hit','Rank','Score','AlnNum','Length','Identity','QryLen','QryStart','QryEnd','HitStart','HitEnd','HitStrand']
                vdb = exobj.localToHitSum(ldb,'hitsum')
                #i# Add QryDesc to hitsum table
                vdb.addField('QryDesc',after='Identity')
                seqdict = seqlist.seqNameDic()
                for entry in vdb.entries():
                    try:
                        seqdesc = seqlist.seqDesc(seqdict[entry['Qry']])
                        if not seqdesc: seqdesc = entry['Qry']
                    except: seqdesc = entry['Qry']
                    entry['QryDesc'] = seqdesc
                #i# Local phase makes sense for protein but not hitsums, which will be protein positions?
                if self.getStrLC('QryType') in ['cds']: sfields.append('Phase')
                hdb = self.db().copyTable(vdb,'html')
                hdb.setFields(sfields)
                vdb.saveToFile('%s.hitsum.tdt' % (self.baseFile()),savefields=sfields)
                self.dict['Output']['locgff'] = '%s.gff3' % (self.baseFile())
                exobj.gff3(ldb,vdb,exmodel,self.dict['Output']['locgff'])
            else:
                sfields = ['Qry','Hit','Rank','Score','QryStart','QryEnd','HitStart','HitEnd','HitStrand']
                rules = {'Score':'sum','QryStart':'min','QryEnd':'max','HitStart':'min','HitEnd':'max','AlnID':'list','HitStrand':'list','Length':'sum','Identity':'sum'}
                vdb.compress(['Qry','Rank','Hit'],rules,'text')
                vdb.dataFormat({'Score':'int','QryStart':'int','QryEnd':'int','HitStart':'int','HitEnd':'int','Length':'int','Identity':'int'})
                vdb.saveToFile('%s.hitsum.tdt' % (self.baseFile()),savefields=sfields)
                hdb = vdb
            self.dict['Output']['hitsum'] = '%s.hitsum.tdt' % (self.baseFile())
            self.dict['Output']['gablam'] = 'No GABLAM table generated for Exonerate searches.'

            #!# Get rid of Sequences; Get fields in right order for HTML. Dump data-sort: not working #!#

            ### ~ [3] Generate HTML table with links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('GenomeID'):
                #if self.dev():
                for entry in hdb.entries():
                    try:
                        spcode = entry['Hit'].split('__')[0].split('_')[-1]
                        if spcode not in self.dict['GenomeID']: spcode = '*'
                    except:
                        spcode = '*'
                    if spcode in self.dict['GenomeID']:
                        self.printLog('#GID','%s -> %s (%s)' % (entry['Hit'],spcode,self.dict['GenomeID'][spcode]))
                        gid = string.split(self.dict['GenomeID'][spcode],'.')
                        public = len(gid) > 1 and gid[1] == 'public'
                        private = len(gid) < 2 or gid[1] == 'private'
                        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/120/jbrowse/index.html?loc=chr_WON710A1__DCMF160804:1232840..1233371&tracks=
                        if public:
                            entry['Hit'] = '<a href="%s%s/jbrowse/index.html?loc=%s:%s..%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),gid[0],entry['Hit'],entry['HitStart'],entry['HitEnd'],self.getStr('Tracks'),entry['Hit'])
                            #x#if private: entry['WebApollo'] += ' | '
                        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/annotator/loadLink?loc=chr_WON710A1__DCMF160804:1232840..1233371&organism=120&tracks=
                        else:
                            entry['Hit'] = '<a href="%sannotator/loadLink?loc=%s:%s..%s&organsim=%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),entry['Hit'],entry['HitStart'],entry['HitEnd'],gid[0],self.getStr('Tracks'),entry['Hit'])
                    else: self.printLog('#GID','%s -> %s (No apollo)' % (entry['Hit'],spcode))
                #else:
                #    gid = string.split(self.getStr('GenomeID'),'.')
                #    public = len(gid) > 1 and gid[1] == 'public'
                #    private = len(gid) < 2 or gid[1] == 'private'
                #    for entry in hdb.entries():
                #        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/120/jbrowse/index.html?loc=chr_WON710A1__DCMF160804:1232840..1233371&tracks=
                #        if public:
                #            entry['Hit'] = '<a href="%s%s/jbrowse/index.html?loc=%s:%s..%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),gid[0],entry['Hit'],entry['HitStart'],entry['HitEnd'],self.getStr('Tracks'),entry['Hit'])
                #            #x#if private: entry['WebApollo'] += ' | '
                #        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/annotator/loadLink?loc=chr_WON710A1__DCMF160804:1232840..1233371&organism=120&tracks=
                #        else:
                #            entry['Hit'] = '<a href="%sannotator/loadLink?loc=%s:%s..%s&organsim=%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),entry['Hit'],entry['HitStart'],entry['HitEnd'],gid[0],self.getStr('Tracks'),entry['Hit'])
                hfile = '%s.html' % self.baseFile()
                #!# Add header and footer if not REST?
                tdtitle = {}    #'Qry,Hit,Rank,Score,EVal,QryLen,HitLen,Qry_AlnLen,Qry_AlnID,Hit_Start,Hit_End'
                html = rje_html.HTML(self.log,self.cmd_list)
                HTML = open(hfile,'w')
                HTML.write(html.htmlHead(title=self.baseFile(),tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
                #HTML.write('<a name="head"><h1>%s Hits</h1></a>\n\n' % self.baseFile())
                HTML.write('<p>Click on the <code>Hit</code> field entries to be taken to the WebApollo genome browser. Login may be required.</p>')
                #!# Add a brief description of what should happen when clicking the links
                fieldsort = {'*':'string','Score':'int','Rank':'int','QryStart':'int','QryEnd':'int','HitStart':'int','HitEnd':'int','AlnID':'string','HitStrand':'int','Length':'int','Identity':'int'}
                #x#fieldsort = {}
                #x#HTML.write(string.replace(rje_html.dbTableToHTML(hdb,tdtitles=tdtitle,datasort=fieldsort),'<table','<table id="parse"'))
                HTML.write(rje_html.dbTableToHTML(hdb,tdtitles=tdtitle,datasort=fieldsort,tabid='parse'))
                HTML.write(html.htmlTail(tabber=False,stupidtable=True))
                HTML.close()
                self.printLog('#HTML','HTML Apollo links output: %s' % hfile)
                self.dict['Output']['html'] = hfile
            else:
                self.printLog('#HTML','No GenomeID given: no HTML Apollo links output.')
                self.dict['Output']['html'] = 'No GenomeID given: no HTML Apollo links output.'
        except: self.errorLog('%s.Exonerate error' % self.prog()); return None
#########################################################################################################################
    def GABLAM(self):   ### Runs GABLAM to generate main outputs.
        '''
        Runs GABLAM to generate main outputs.

        Outputs:
        *.gablam.tdt = gablam = Full GABLAM summary table of each Query-Hit pair [tdt]
        *.hitsum.tdt = hitsum = One line summary of hits per query [tdt]
        *.local.tdt = local = Local BLAST alignments table [tdt]
        *.unique.tdt = unique = Local BLAST alignments table, reduced to non-overlapping query regions [tdt]
        *.local.gff = locgff = Local BLAST alignments table in GFF3 format [gff3]
        *.unique.gff = unigff = Local BLAST alignments table, reduced to non-overlapping query regions, in GFF3 format [gff3]
        *.local.sam = sam = Local BLAST alignments table in SAM alignment format (blastn only) [sam]
        *.unique.sam -> Deleted
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup GABLAM commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdefaults = ['keepblast=T','localmin=30','localidmin=75','blaste=1e-10','tophits=10']
            gcmd = gdefaults + self.cmd_list + ['fullblast=T','localsam=T','localgff=T','qryunique=T','localunique=T','qassemble=T','searchdb=%s' % self.getStr('Genome')]
            if self.getStrLC('Rest'): gcmd.append('blastdir=%s/' % os.path.split(self.baseFile())[0])
            if self.getStrLC('QryType') in ['cds','prot']: gcmd.append('cdsgff=T')
            if self.getStrLC('QryType') in ['cds','dna','transcript','trans']: gcmd.append('blastp=blastn')
            elif self.getStrLC('QryType') in ['prot']: gcmd.append('blastp=tblastn')
            else: raise ValueError('Unrecognised qrytype: %s' % self.getStrLC('QryType'))
            if self.dev(): self.printLog('#DEV','GABLAM Cmd: %s' % string.join(gcmd))
            ### ~ [2] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gab = self.obj['GABLAM'] = gablam.GABLAM(self.log,gcmd)
            gab.gablam()
            for otype in ['gablam','hitsum','local','unique']:
                self.dict['Output'][otype] = '%s.%s.tdt' % (self.baseFile(),otype)
            for otype in ['local','unique']:
                self.dict['Output']['%sgff' % otype[:3]] = '%s.%s.gff' % (self.baseFile(),otype)
            self.dict['Output']['sam'] = '%s.local.sam' % (self.baseFile())
            self.dict['Output']['blast'] = '%s.blast' % (self.baseFile())
            self.dict['Output']['seqin'] = 'SeqIn'
            usam = '%s.unique.sam' % (self.baseFile())
            if rje.exists(usam): os.unlink(usam)
            ### ~ [3] Generate HTML table with links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('GenomeID'):
                gfile = '%s.gablam.tdt' % (self.baseFile())
                hdb = self.db().addTable(gfile,name='html',mainkeys=['Qry','Hit'])
                hdb.keepFields(string.split('Qry,Hit,Rank,Score,EVal,QryLen,HitLen,Qry_AlnLen,Qry_AlnID,Hit_Start,Hit_End',','))
                hdb.dataFormat({'Rank':'int'})
                hdb.newKey(['Qry','Rank','Hit'])
                #hdb.addField('WebApollo')
                #gid = string.split(self.getStr('GenomeID'),'.')
                #public = len(gid) > 1 and gid[1] == 'public'
                #private = len(gid) < 2 or gid[1] == 'private'
                for entry in hdb.entries():
                    try:
                        spcode = entry['Hit'].split('__')[0].split('_')[-1]
                        if spcode not in self.dict['GenomeID']: spcode = '*'
                    except:
                        spcode = '*'
                    if spcode in self.dict['GenomeID']:
                        gid = string.split(self.dict['GenomeID'][spcode],'.')
                        public = len(gid) > 1 and gid[1] == 'public'
                        private = len(gid) < 2 or gid[1] == 'private'
                        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/120/jbrowse/index.html?loc=chr_WON710A1__DCMF160804:1232840..1233371&tracks=
                        if public:
                            entry['Hit'] = '<a href="%s%s/jbrowse/index.html?loc=%s:%s..%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),gid[0],entry['Hit'],entry['Hit_Start'],entry['Hit_End'],self.getStr('Tracks'),entry['Hit'])
                            #x#if private: entry['WebApollo'] += ' | '
                        # https://edwapollo.babs.unsw.edu.au/apollo-2.0.8/annotator/loadLink?loc=chr_WON710A1__DCMF160804:1232840..1233371&organism=120&tracks=
                        else:
                            entry['Hit'] = '<a href="%sannotator/loadLink?loc=%s:%s..%s&organsim=%s&tracks=%s" target="_blank">%s</a>' % (self.getStr('WebApollo'),entry['Hit'],entry['Hit_Start'],entry['Hit_End'],gid[0],self.getStr('Tracks'),entry['Hit'])
                hfile = '%s.html' % self.baseFile()
                #!# Add header and footer if not REST?
                tdtitle = {}    #'Qry,Hit,Rank,Score,EVal,QryLen,HitLen,Qry_AlnLen,Qry_AlnID,Hit_Start,Hit_End'
                html = rje_html.HTML(self.log,self.cmd_list)
                HTML = open(hfile,'w')
                HTML.write(html.htmlHead(title=self.baseFile(),tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
                #HTML.write('<a name="head"><h1>%s Hits</h1></a>\n\n' % self.baseFile())
                HTML.write('<p>Click on the <code>Hit</code> field entries to be taken to the WebApollo genome browser. Login may be required.</p>')
                #!# Add a brief description of what should happen when clicking the links
                HTML.write(rje_html.dbTableToHTML(hdb,tdtitles=tdtitle))
                HTML.write(html.htmlTail(tabber=False))
                HTML.close()
                self.printLog('#HTML','HTML Apollo links output: %s' % hfile)
                self.dict['Output']['html'] = hfile
            else:
                self.printLog('#HTML','No GenomeID given: no HTML Apollo links output.')
                self.dict['Output']['html'] = 'No GenomeID given: no HTML Apollo links output.'
            return gab
        except: self.errorLog('%s.GABLAM error' % self.prog()); return None
#########################################################################################################################
    def QBLAST(self):   ### Runs GABLAM to generate main outputs.
        '''
        ### ~ QAssemble Options (single query only) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        qassemblefas=T/F: Special mode for running with outfmt=4 and then converting to fasta file [True]
        qcomplete=T/F   : Whether the query sequence should be full-length in qassemblefas output [True]
        qconsensus=X    : Whether to convert QAssemble alignments to consensus sequences (None/Hit/Full) [Hit]
        qfasdir=PATH    : Output directory for QAssemble alignments [./]

        Outputs:
        *.qblast = qblast = BLAST results formatted in query assembly `-outfmt 4`
        *.<QRY>.qaln.hitcons.fas -> *.qaln.fas = qaln = Query sequence aligned to consensus for each hit from QBLAST
        *.qhits.fas = qhits = Degapped consensus hits from QBLAST

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup QBLAST commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gab = self.obj['GABLAM']
            seqlist = gab.obj['SeqList']
            seq = seqlist.nextSeq()
            bcmd = ['qassemblefas=T','qcomplete=T','qconsensus=Hit','qfasdir=%s/' % os.path.split(self.baseFile())[0],'blasto=%s.qblast' % self.baseFile(),'blastd=%s' % self.getStr('Genome'),'blasti=%s' % self.getStr('SeqIn')]
            if self.getStrLC('QryType') in ['cds','dna','transcript','trans']: bcmd.append('blastp=blastn')
            elif self.getStrLC('QryType') in ['prot']: bcmd.append('blastp=tblastn')
            else: raise ValueError('Unrecognised qrytype: %s' % self.getStrLC('QryType'))
            self.printLog('#DEV','BLAST Cmd: %s' % string.join(bcmd))
            ### ~ [2] Run QBLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje_blast.BLASTRun(self.log,['tuplekeys=T']+gab.cmd_list+bcmd).run()
            qfas = '%s.%s.qaln.hitcons.fas' % (self.baseFile(),seqlist.shortName())
            self.printLog('#QFAS','%s: %s' % (qfas,rje.exists(qfas)))
            qout = '%s.qaln.fas' % (self.baseFile())
            os.rename(qfas,qout)
            hout = '%s.qhits.fas' % (self.baseFile())
            HOUT = open(hout,'w'); hx =0
            for qline in open(qout,'r').readlines()[2:]:
                if qline.startswith('>'): HOUT.write(qline); hx += 1
                else: HOUT.write(string.replace(qline,'-',''))
            HOUT.close()
            self.printLog('#QHITS','%s degapped consensus hit sequences saved to %s' % (rje.iStr(hx),hout))
            self.dict['Output']['qblast'] = '%s.qblast' % self.baseFile()
            self.dict['Output']['qaln'] = qout
            self.dict['Output']['qhits'] = hout
            #?# Degap
            # seqsuite.py seq -seqin QFAS/Nscutatus.k81L500.tlr1.TLR1_HUMAN__Q15399.qaln.hitcons.fas degap basefile=tlr1 -seqout tlr1hits.fas
            return True
        except: self.errorLog('%s.QBLAST error' % self.prog()); return None
#########################################################################################################################
### End of SECTION II: Apollo Class                                                                                     #
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
    try: Apollo(mainlog,['basefile=apollo']+cmd_list).run()

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
