#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       PAGSAT
Description:  Pairwise Assembled Genome Sequence Analysis Tool
Version:      1.6.1
Last Edit:    30/10/15
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for the assessment of an assembled genome versus a suitable reference. For optimal results, the
    reference genome will be close to identical to that which should be assembled. However, comparative analyses should
    still be useful when different assemblies are run against a related genome - although there will not be the same
    expectation for 100% coverage and accuracy, inaccuracies would still be expected to make an assembly less similar
    to the reference.

    Main input for PAGSAT is an assembled genome in fasta format (`assembly=FILE`) and a reference genome in fasta format
    with corresponding `*.gb` genbank download for feature extraction. The

Output:
    Main output is a number of delimited text files and PNG graphics made with R. Details to follow.

Commandline:
    ### ~ Input/Setup Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    assembly=FILE   : Fasta file of assembled contigs to assess [None]
    refgenome=FILE  : Fasta file of reference genome for assessment (also *.gb for full functionality) [None]
    spcode=X        : Species code for reference genome (if not already processed by rje_genbank) [None]
    minqv=X         : Minimum mean QV score for assembly contigs (read from *.qv.csv) [20]
    mincontiglen=X  : Minimum contig length to retain in assembly (QV filtering only) [1000]
    casefilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    ### ~ Reference vs Assembly Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minlocid=X      : Minimum percentage identity for local hits mapping to chromosome coverage [99.0]
    minloclen=X     : Mininum length for local hits mapping to chromosome coverage [1000]
    genesummary=T/F : Whether to include reference gene searches in summary data [True]
    protsummary=T/F : Whether to include reference protein searches in summary data [True]
    tophitbuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]
    diploid=T/F     : Whether to treat assembly as a diploid [False]
    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : Basename for output files and directories. [assembly+ref]
    chromalign=T/F  : Whether to align chromosomes with contigs [True]
    rgraphics=T/F   : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
    dotplots=T/F    : Whether to use gablam.r to output dotplots for all ref vs assembly. [False]
    report=T/F      : Whether to generate HTML report [True]
    ### ~ Comparison Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    compare=FILES   : Compare assemblies selected using a list of *.Summary.tdt files (wildcards allowed). []
    fragcov=LIST    : List of coverage thresholds to count min. local BLAST hits (checks integrity) [50,90,95,99]
    chromcov=LIST   : Report no. of chromosomes covered by a single contig at different %globID (GABLAM table) [95,98,99]
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
import rje, rje_db, rje_genbank, rje_html, rje_obj, rje_seqlist, rje_sequence, rje_tree, rje_tree_group, rje_xref
import rje_blast_V2 as rje_blast
import rje_dismatrix_V3 as rje_dismatrix
import gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0.0 - Initial working version for based on rje_pacbio assessment=T.
    # 1.1.0 - Fixed bug with gene and protein summary data. Removed gene/protein reciprocal searches. Added compare mode.
    # 1.1.1 - Added PAGSAT output directory for tidiness!
    # 1.1.2 - Renamed the PacBio class PAGSAT.
    # 1.2.0 - Tidied up output directories. Added QV filter and Top Gene/Protein hits output.
    # 1.2.1 - Added casefilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    # 1.3.0 - Added tophitbuffer=X and initial synteny analysis for keeping best reference hits.
    # 1.4.0 - Added chrom-v-contig alignment files along with *.ordered.fas.
    # 1.4.1 - Made default chromalign=T.
    # 1.4.2 - Fixed casefilter=F.
    # 1.5.0 - diploid=T/F     : Whether to treat assembly as a diploid [False]
    # 1.6.0 - mincontiglen=X  : Minimum contig length to retain in assembly [1000]
    # 1.6.1 - Added diploid=T/F to R PNG call.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add reduced functionality if genbank file not given (i.e. no features).
    # [ ] : Add HTML report output.
    # [ ] : Add indels to QAssemble
    # [ ] : Have a min minloclen (=localcut) but then try increasing by 1kb chunks without increasing "N" periods?
    # [Y] : Calculate the difference from Reference as a number for comparisons.
    # [ ] : Pull out "N" regions of the reference - QAssemble back against the pre-assembly and subreads.
    # [ ] : Add (interactive?) reformatting of *.full.fas for refgenome input.
    # [Y] : Include minloclen in relevant file names - GABLAM and default basefile.
    # [Y] : Consider using minloclen for localcut=X GABLAM Cut-off length for local alignments contributing to global stats.
    # [Y] : Option to switch off Gene and Protein searches for increased speed. (Need to edit R summGraph() too)
    # [ ] : Improve gene search to extend hits to full length of genes/proteins.
    # [ ] : Separate minloclen and localcut?
    # [ ] : Add thumbnails=T/F for report and R graphics
    # [ ] : Add maxcontig=X for a single summary page. (Ask to proceed if more contigs?)
    # [ ] : Add summary stats for assembly and reference? (Could load them into SeqList objects)
    # [Y] : Add reading of unitig coverage and quality scores from quiver (*.qv.csv files)
    # [Y] : Need to tidy up outputs. Generate a *.PAGSAT/ directory for most outputs: update R script accordingly.
    # [N] : blastdir=PATH   : Path for blast results file (unless keepblast=F) [./assembly.BLAST/]
    # [Y] : Add separate GABLAM directory for feature searches? (./assembly.GABLAM/) (These ignore cutoffs, right?)
    # [Y] : Add (and distinguish) minlocid to output file names as well as minloclen. (Make integer? *.LXXX.IDXX.*)
    # [X] : Contemplate setting softmask=F for BLAST searches. (Why are some genes missing?!)
    # [ ] : Consider replacing GABLAM fragfas with own gene extraction algorithm that extends ORFs & combines exons.
    # [ ] : Move R graphics such that it can be re-run in isolation (if summary.png missing).
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('PAGSAT', '1.6.1', 'October 2015', '2015')
    description = 'Pairwise Assembled Genome Sequence Analysis Tool'
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
### SECTION II: PAGSAT Class                                                                                            #
#########################################################################################################################
class PAGSAT(rje_obj.RJE_Object):
    '''
    PAGSAT Class. Author: Rich Edwards (2015).

    Str:str
    - Assembly=FILE   : Fasta file of assembled contigs to assess [None]
    - BaseBase        : Path-trimmed Basefile for outputs that do not involve data filtering etc.
    - CutBase         : Path-trimmed Basefile for outputs that include data filtering etc.
    - GABLAMDir       : Parent directory for all BLAST and GABLAM searches.
    - RefBase=X       : Basefile for reference genome for assessment (*.gb) [None]
    - RefGenome=FILE  : Fasta file of reference genome for assessment (also *.gb for full functionality) [None]
    - ResDir          : Results directory = BASEFILE.PAGSAT/

    Bool:boolean
    - CaseFilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    - ChromAlign=T/F  : Whether to align chromosomes with contigs (slow!) [False]
    - Diploid=T/F     : Whether to treat assembly as a diploid [False]
    - DotPlots=T/F    : Whether to use gablam.r to output dotplots for all ref vs assembly. [False]
    - GeneSummary=T/F : Whether to include reference gene searches in summary data [True]
    - ProtSummary=T/F : Whether to include reference protein searches in summary data [True]
    - RGraphics=T/F   : Whether to generate PNG graphics using R [True]
    - Report=T/F      : Whether to generate HTML report [True]

    Int:integer
    - MinContigLen=X  : Minimum contig length to retain in assembly [1000]
    - MinLocLen=X     : Mininum length for local hits mapping to chromosome coverage [100]
    - MinQV=X         : Minimum mean QV score for assembly contigs (read from *.qv.csv) [20]

    Num:float
    - MinLocID=X      : Minimum percentage identity for local hits mapping to chromosome coverage [0.99]
    - TopHitBuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]

    File:file handles with matching str filenames
    
    List:list
    - ChromCov=LIST   : Report no. of chromosomes covered by a single contig at different %globID (GABLAM table) [95,98,99]
    - Compare=FILES   : Special mode to compare a list of *.Summary.tdt files (wildcards allowed). []
    - FragCov=LIST    : List of coverage thresholds to count min. local BLAST hits (checks integrity) [50,90,95,99]

    Dict:dictionary
    - QV              : Dictionary of contig:QV score read from *.qv.csv

    Obj:RJE_Objects
    - Features  : Reference Features Database table (reused for Compare=T).
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Assembly','BaseBase','CaseFilter','CutBase','GABLAMDir','RefBase','RefGenome','ResDir']
        self.boollist = ['ChromAlign','Diploid','DotPlots','GeneSummary','ProtSummary','RGraphics','Report']
        self.intlist = ['MinContigLen','MinLocLen','MinQV']
        self.numlist = ['MinLocID','TopHitBuffer']
        self.filelist = []
        self.listlist = ['ChromCov','Compare','FragCov']
        self.dictlist = ['QV']
        self.objlist = ['DB','Features']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GABLAMDir':rje.makePath('GABLAM/'),'ResDir':rje.makePath('PAGSAT/')})
        self.setBool({'CaseFilter':True,'ChromAlign':True,'Diploid':False,'DotPlots':False,
                      'GeneSummary':True,'ProtSummary':True,'RGraphics':True,'Report':True})
        self.setInt({'MinLocLen':1000,'MinQV':20,'MinContigLen':1000})
        self.setNum({'MinLocID':99.0,'TopHitBuffer':1.0})
        self.list['ChromCov'] = [95,98,99]
        self.list['FragCov'] = [50,90,95,99]
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
                #self._cmdReadList(cmd,'str',[])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Assembly','RefGenome'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['CaseFilter','ChromAlign','Diploid','DotPlots','GeneSummary','ProtSummary','RGraphics','Report'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinContigLen','MinLocLen','MinQV'])   # Integers
                self._cmdReadList(cmd,'float',['TopHitBuffer']) # Floats
                self._cmdReadList(cmd,'perc',['MinLocID'])
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',[])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['ChromCov','FragCov'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Compare']) # List of files using wildcards and glob
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
            if self.list['Compare']: return self.compare()
            if self.getBool('Report'): return self.report()
            else: return self.assessment()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Setup ~~~~~ ##')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            #if self.getStrLC('REST') and not self.basefile(return_none=''): self.basefile('pagsat')
            if self.list['Compare']:
                if not self.basefile(return_none=''): self.basefile('pagsat')
                return True
            checkfiles = []
            ## ~ [1a] Reference Genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.exists('%s.gb' % self.getStr('RefGenome')): self.setStr({'RefBase':self.getStr('RefGenome')})
            else: self.setStr({'RefBase':rje.baseFile(self.getStr('RefGenome'))})
            rungb = False   # Whether to run rje_genbank on RefGenome
            for rfile in ['full.fas','gene.fas','prot.fas','Feature.tdt']:
                gfile = '%s.%s' % (self.getStr('RefBase'),rfile)
                self.printLog('#CHECK','%s: %s' % (gfile,{True:'Found.',False:'Missing!'}[os.path.exists(gfile)]))
                rungb = rungb or not os.path.exists(gfile)
                checkfiles.append(gfile)
            #self.debug('Run Genbank: %s' % rungb)
            if rungb:
                gcmd = ['protacc=locus_tag','details=product,gene_synonym,note,db_xref']   # Defaults
                gcmd += self.cmd_list   # Can over-ride/add. This include spcode=X
                gcmd += ['seqin=%s.gb' % self.getStr('RefBase'),'taxdir=','tabout=T','fasout=full,gene,cds,prot']
                rje_genbank.GenBank(self.log,gcmd).run()
                for cfile in checkfiles:
                    if not rje.exists(cfile): raise IOError('Cannot find %s!' % cfile)
            if not self.getStr('RefGenome').endswith('.fas') or not rje.exists(self.getStr('RefGenome')):
                self.printLog('#NAMES','%s.full.fas sequence names will not be suitable.' % self.getStr('RefBase'))
                self.printLog('#NAMES','Please modify gene names in %s.full.fas (e.g. ChrX) and save as %s.fas re-run.' % (self.getStr('RefBase'),self.getStr('RefBase')))
                self.printLog('#NAMES','Then re-run with refgenome=%s.fas.' % self.getStr('RefBase'))
                return False
            self.printLog('#REF','Reference genome fasta: %s' % self.getStr('RefGenome'))
            if not rje.exists(self.getStr('RefGenome')): raise IOError('Cannot find RefGenome: %s!' % self.getStr('RefGenome'))
            ## ~ [1b] Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.exists(self.getStr('Assembly')): raise IOError('Cannot find Assembly: %s!' % self.getStr('Assembly'))
            qvfile = '%s.qv.csv' % rje.baseFile(self.getStr('Assembly'))
            if rje.exists(qvfile):
                qdb = self.db().addTable(qvfile,['contig_id'],name='QV')
                qdb.dataFormat({'mean_coverage':'num','mean_qv':'num'})
                qdb.addField('Contig'); qdb.addField('Seq'); qdb.addField('Len')
                assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=file','usecase=T'])
                qx = 0
                while assembly.nextSeq():
                    seqname = assembly.seqName()
                    sname = assembly.shortName()
                    for qentry in qdb.entries():
                        if '%s|' % qentry['contig_id'] in seqname:
                            if qentry['Seq']: raise ValueError('Multiple QV sequence mapping error! (%s)' % seqname)
                            qentry['Seq'] = assembly.obj['Current']; qentry['Contig'] = sname
                            qentry['Len'] = assembly.seqLen()
                            if sname in self.dict['QV']: raise ValueError('Multiple QV sequence mapping error! (%s)' % seqname)
                            self.dict['QV'][sname] = qentry['mean_qv']; qx += 1
                if qx != qdb.entryNum(): raise ValueError('Only %d of %d QV values mapped to sequences!' % (qx,qdb.entryNum()))
                if qx != assembly.seqNum(): raise ValueError('%d QV values but %d sequences!' % (qx,assembly.seqNum()))
                # Filter qdb on XCov and/or QV and see if entries lost
                qvcut = self.getInt('MinQV')
                lencut = self.getInt('MinContigLen')
                qvcutlen = 0
                #qdb.dropEntries(['mean_qv<%d' % qvcut],inverse=False,log=True,logtxt='QV Filtering')
                for qentry in qdb.sortedEntries('Contig'):
                    if qentry['mean_qv'] < qvcut:
                        self.printLog('#MINQV','%s failed to meet QV>=%d: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],qvcut,qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                        qvcutlen += qentry['Len']
                        qdb.dropEntry(qentry)
                    elif qentry['Len'] < lencut:
                        self.printLog('#MINLEN','%s failed to meet Len>=%d: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],lencut,qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                        qvcutlen += qentry['Len']
                        qdb.dropEntry(qentry)
                    else:
                        self.printLog('#QV','%s: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                if qx != qdb.entryNum() or self.getBool('CaseFilter'):
                    self.printLog('#QV','%d of %d contigs (%.2f kb) fail to meet QV>=%d' % (qx-qdb.entryNum(),qx,qvcutlen/1000.0,qvcut))
                    abase = rje.baseFile(self.getStr('Assembly'),strip_path=True)
                    qvbase = '%s.QV/%s.qv%d' % (abase,abase,qvcut)
                    rje.mkDir(self,qvbase)
                    qvfas = '%s.fas' % (qvbase)
                    if rje.exists(qvfas) and not self.force():
                        self.printLog('#QV','QVFas file found (force=F): %s' % qvfas)
                        QVFAS = None
                    else:
                        # Make new Assembly file in new directory. BASEFILE.Assembly.QVX
                        rje.backup(self,qvfas)
                        QVFAS = open(qvfas,'w')
                    sx = 0; trimx = 0
                    for seq in qdb.indexKeys('Seq'):
                        (name,sequence) = assembly.getSeq(seq)
                        if self.getBool('CaseFilter'):
                            i = 0; slen = j = len(sequence)
                            while i < len(sequence) and sequence[i] == sequence[i].lower(): i += 1
                            while j and sequence[j-1] == sequence[j-1].lower(): j -= 1
                            sequence = sequence[i:j]
                            trimx += slen - len(sequence)
                            if len(sequence) != slen: self.printLog('#TRIM','%s: %s bp QV trimmed.' % (assembly.shortName(seq),rje.iStr(slen-len(sequence))))
                            if len(sequence) < lencut: self.printLog('#TRIM','%s: %s bp QV trimmed => now too short!' % (assembly.shortName(seq),rje.iStr(slen-len(sequence)))); continue
                            if not sequence: self.warnLog('No quality (upper case) sequence for %s!' % name); continue
                        if QVFAS: QVFAS.write('>%s\n%s\n' % (name, sequence)); sx += 1
                    if self.getBool('CaseFilter'):
                        self.printLog('#TRIM','Total assembly QV trimmed: %.3f kb' % (trimx/1000.0))
                    if QVFAS:
                        QVFAS.close()
                        self.printLog('\r#OUT','%s of %s sequences output to %s.' % (sx,assembly.seqNum(),qvfas))
                        # Save qdb too!
                        qdb.saveToFile('%s.qv.csv' % qvbase)
                    self.setStr({'Assembly':qvfas})
                else: self.printLog('#QV','%d of %d contigs meet QV>=%d' % (qdb.entryNum(),qx,qvcut))
            else: self.warnLog('%s not found: no QV filtering/reporting/' % qvfile)
            ### ~ [2] Report on Key settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#INPUT','Assembly file: %s' % self.getStr('Assembly'))
            self.printLog('#INPUT','Reference genome file: %s' % self.getStr('RefGenome'))
            self.printLog('#INPUT','Reference genbank file: %s.gb' % self.getStr('RefBase'))
            minloclen = self.getInt('MinLocLen')
            self.printLog('#PARAM','Min. local BLAST alignment length: %sbp ("L%d")' % (rje.iStr(minloclen),minloclen))
            minlocid = self.getNum('MinLocID')
            self.printLog('#PARAM','Min. local BLAST alignment %%identity: %s%% ("ID%d")' % (rje.sf(minlocid,3),int(minlocid)))
            if not self.baseFile(return_none=None):
                self.baseFile('%s.%s' % (rje.baseFile(self.getStr('Assembly'),strip_path=True),rje.baseFile(self.getStr('RefGenome'),strip_path=True)))
            basedir = rje.makePath('%s.PAGSAT/' % self.baseFile())
            gabdir = rje.makePath('%s.GABLAM/' % self.baseFile())
            rje.mkDir(self,basedir)
            rje.mkDir(self,gabdir)
            self.setStr({'GABLAMDir':gabdir,'ResDir':basedir,'BaseBase':os.path.basename(self.baseFile())})
            self.printLog('#GABDIR',self.getStr('GABLAMDir'))   # Directory for GABLAM and BLAST output
            self.printLog('#PAGDIR',self.getStr('ResDir'))      # Directory for PAGSAT output
            # Primary basefile is PAGSAT directory with additional cut-off information added
            self.baseFile('%s%s.L%dID%d' % (basedir,os.path.basename(self.baseFile()),minloclen,int(minlocid)))
            self.setStr({'CutBase':os.path.basename(self.baseFile())})
            self.printLog('#BASE',self.getStr('BaseBase'))      # Root basefile
            self.printLog('#PAGOUT','%s.*' % self.baseFile())   # Primary output basefile
            if self.baseFile() != self.fileBase(): raise ValueError('Dev problem! Please report. Avoid use of basefile=X.')
            self.db().basefile(self.basefile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def fileBase(self,resdir='Res',base='Cut',extra=None):   ### Returns appropriate file output basename
        '''
        Returns appropriate file output basename. Default should return same as self.baseFile()
        @param dir:str ['Res'] = Whether output directory is 'GABLAM' or 'Res'.
        @param base:str ['Cut'] = Whether file basename included cutoffs ('Cut') or not ('Base')
        @return: path constructed from GABLAMDir/ResDir and CutBase/BaseBase.
        '''
        filebase = '%s%s' % (self.getStr('%sDir' % resdir),self.getStr('%sBase' % base))
        if extra: filebase = '%s.%s' % (filebase,extra)
        return filebase
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        coverage = main results table
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['coverage']
#########################################################################################################################
    ### <3> ### PAGSAT GABLAM Methods                                                                                   #
#########################################################################################################################
    def runGABLAM(self,gabcmd,gtype='Reference',blastgz=True): ### Runs GABLAM with given commands and basefile, managing *blast file.
        '''
        Runs GABLAM with given commands and basefile, managing *blast file. Adds minloclen for actual run.
        >> gabcmd:list = List of GABLAM run commands. Will extract basefile from this list.
        >> blastgz:bool [True] = whether to make/process a general BLAST file if possible. (e.g. rename and g(un)zip)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = None
            for gcmd in gabcmd:
                if gcmd.startswith('basefile='): basefile = string.split(gcmd,'=',1)[1]
            self.printLog('#~~#','## ~~~~~ %s GABLAM ~~~~~ ##' % basefile)
            # Set blastbase = basefile name for blast file
            if gtype == 'Self': blastbase = string.join(string.split(basefile,'.')[:-1],'.')
            else: blastbase = self.fileBase('GABLAM','Base',gtype)
            #if basefile.endswith('.%d' % self.getInt('MinLocLen')): blastbase = string.join(string.split(basefile,'.')[:-1],'.')
            #else: blastbase = basefile
            ### ~ [1] Run GABLAM, processing BLAST file as required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blastfile = '%s.blast' % blastbase
            # Unzip file if required
            if rje.exists('%s.gz' % blastfile) and blastgz:
                os.system('gunzip %s.gz' % blastfile)
                self.printLog('#GUNZIP','%s unzipped.' % blastfile)
            # Rename BLAST file if found
            if rje.exists(blastfile) and blastbase != basefile:
                os.rename(blastfile,'%s.blast' % basefile)
                self.printLog('#BLAST','%s -> %s.blast' % (blastfile,basefile))
            elif not rje.exists(blastfile): self.printLog('#BLAST','%s not found.' % (blastfile))
            # Run GABLAM
            gablam.GABLAM(self.log,gabcmd).run()
            # Rename BLAST file if kept
            if rje.exists('%s.blast' % basefile) and blastbase != basefile:
                os.rename('%s.blast' % basefile,blastfile)
                self.printLog('#BLAST','%s.blast -> %s' % (basefile,blastfile))
            # (Re)zip BLAST file if required
            if rje.exists(blastfile) and blastgz:
                os.system('gzip %s' % blastfile)
                self.printLog('#GZIP','%s (re)zipped.' % blastfile)
            return True
        except: self.errorLog('%s.runGABLAM() error' % self.prog()); return False
#########################################################################################################################
    def qAssembleGABLAM(self,gtype='Reference'):   ### Generate GABLAM analyses for assembly assessment.
        '''
        Generate GABLAM analyses for assembly assessment.
        >> gtype:str ['qassemble'] = type of GABLAM to run (Reference/Assembly/Genes/Proteins/Genes.Reciprocal/Proteins.Reciprocal)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Process reference gb file if required and find genome, genes and proteins
            refbase = self.getStr('RefBase')    # Used for searches of genes and proteins
            gdir = self.getStr('GABLAMDir')     # GABLAMDir. Also used as BLAST directory
            if gtype in ['Reference','Assembly']:  # These searches used cutoff data
                gbase = self.fileBase('GABLAM','Cut',gtype)
            else: gbase = self.fileBase('GABLAM','Base',gtype)     # Other searches do not need cutoffs.

            #basefile = self.baseFile()
            #gdir = '%s.GABLAM/' % basefile
            #if gtype not in ['Assembly','Reciprocal'] and basefile.endswith('.%d' % self.getInt('MinLocLen')):
            #    gdir = '%s.GABLAM/' % string.join(string.split(basefile,'.')[:-1],'.')
            #    gbase = '%s%s' % (gdir,string.join(string.split(self.baseFile(strip_path=True),'.')[:-1],'.'))
            #else: gbase = '%s%s' % (gdir,self.baseFile(strip_path=True))

            self.printLog('#GDIR',gdir)
            rje.mkDir(self,gdir,log=True)
            if not self.force():
                runfound = True
                if gtype == 'Reference' and not rje.exists('%s.hitsum.tdt' % gbase): runfound = False
                if gtype == 'Reference' and not rje.exists('%s.gablam.tdt' % gbase): runfound = False
                if gtype == 'Reference' and not rje.exists('%s.local.tdt' % gbase): runfound = False
                if gtype != 'Reference' and not rje.exists('%s.hitsum.tdt' % (gbase)): runfound = False
                self.printLog('#GABLAM','%s GABLAM %s.* found: %s' % (gtype,gbase,runfound))
                if runfound: return True
            ## Might want to set forks for speed up
            if self.getInt('Forks') < 2: self.printLog('#INFO','Consider setting forks=X to speed up multiple BLASTs.')
            ## GABLAM defaults that can be over-ridden by commandline
            gabdefault = ['keepblast=T','fullblast=T','blastdir=%s' % gdir,'dismat=F','distrees=F','disgraph=F','dotplots=F']
            ## GABLAM options that must be set
            gabcmd = gabdefault + self.cmd_list + ['qassemble=T','outstats=GABLAM','qryacc=F','percres=T','localcut=%d' % self.getInt('MinLocLen'),'basefile=%s' % gbase]
            fascmd = ['fasout=T','fragfas=T','combinedfas=T','localcut=0']

            ### ~ [1] Perform different GABLAM searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # 1. QAssemble GABLAM of the reference genome against the assembly to get full genome coverage.
            if gtype == 'Reference':
                gcmd = ['seqin=%s' % self.getStr('RefGenome'),'searchdb=%s' % self.getStr('Assembly'),'dna=T','blastp=blastn','dotplots=%s' % self.getBool('DotPlots'),'dotlocalmin=%s' % self.getInt('MinLocLen')]
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 2. QAssemble GABLAM of the assembly versus reference genome to assess assembly accuracy in terms of excess sequence.
            if gtype == 'Assembly':
                gcmd = ['seqin=%s' % self.getStr('Assembly'),'searchdb=%s' % self.getStr('RefGenome'),'dna=T','blastp=blastn']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 3. QAssemble GABLAM of the reference genes (from Genbank annotation) to assess accuracy in terms of annotated features.
            if gtype == 'Genes':
                gcmd = fascmd + ['seqin=%s.gene.fas' % refbase,'searchdb=%s' % self.getStr('Assembly'),'dna=T','blastp=blastn','fasdir=%s.GeneHits/'% gbase]
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Fragments search of hit genes
            if gtype == 'Genes.Fragments':
                gcmd = ['localcut=0','seqin=%s.gene.fas' % refbase,'searchdb=%s.fas' % gbase[:-10],'dna=T','blastp=blastn']
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Reciprocal search of hit genes
            if gtype == 'Genes.Reciprocal':
                gcmd = ['localcut=0','searchdb=%s.gene.fas' % refbase,'seqin=%s.fas' % gbase[:-11],'dna=T','blastp=blastn']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 4. QAssemble GABLAM of the reference proteins (from Genbank annotation) to assess accuracy in terms of proteome coverage.
            if gtype == 'Proteins':
                gcmd = fascmd + ['seqin=%s.prot.fas' % refbase,'searchdb=%s' % self.getStr('Assembly'),'blastp=tblastn','fasdir=%s.ProtHits/'% gbase]
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Fragments search of hit protein-coding sequences by proteins
            if gtype == 'Proteins.Fragments':
                gcmd = ['localcut=0','seqin=%s.prot.fas' % refbase,'searchdb=%s.fas' % gbase[:-10],'blastp=tblastn']
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Reciprocal search of hit protein-coding sequences versus proteins
            #Q# Should these be translated first?
            if gtype == 'Proteins.Reciprocal':
                gcmd = ['localcut=0','searchdb=%s.prot.fas' % refbase,'seqin=%s.fas' % gbase[:-11],'blastp=blastx']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
        except: self.errorLog('%s.qAssembleGABLAM error' % self.prog())
#########################################################################################################################
    def selfQAssemble(self,genome=None):    ### Runs QAssemble GABLAM against self.
        '''
        Runs QAssemble GABLAM against self.
        >> genome:str = Genome to self QAssemble. Uses self.getStr('RefGenome') if None.
        << selfdb = Database Object with HitSum and Local tables loaded.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not genome: genome = self.getStr('RefGenome')
            genbase = '%s.L%dID%d' % (rje.baseFile(genome),self.getInt('MinLocLen'),self.getInt('MinLocID'))
            selfdb = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % genbase])
            ## ~ [1] Check for existing run data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.force():
                hdb = selfdb.addTable(mainkeys=['Qry'],name='hitsum',expect=False)
                gdb = selfdb.addTable(mainkeys=['Qry','Hit'],name='gablam',expect=False)
                ldb = selfdb.addTable(mainkeys=['Qry','Hit','AlnNum'],name='local',expect=False)
                if hdb and ldb and gdb: return selfdb
            ## ~ [2] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Might want to set forks for speed up
            if self.getInt('Forks') < 2: self.printLog('#INFO','Consider setting forks=X to speed up multiple BLASTs.')
            ## Process reference gb file if required and find genome, genes and proteins
            blastdir = rje.makePath(os.path.dirname(genbase))
            gabdefault = ['keepblast=T','fullblast=T','blastdir=%s' % blastdir,'dismat=F','distrees=F','disgraph=F','dotplots=F']
            ## GABLAM options that must be set
            gabcmd = gabdefault + self.cmd_list + ['qassemble=T','outstats=GABLAM','qryacc=F','percres=T']
            gabcmd += ['seqin=%s' % genome,'dna=T','blastp=blastn','basefile=%s' % genbase,'selfhit=T','selfsum=T']
            #gablam.GABLAM(self.log,gabcmd).run()
            self.runGABLAM(gabcmd,'Self')
            hdb = selfdb.addTable(mainkeys=['Qry'],name='hitsum',expect=False)
            gdb = selfdb.addTable(mainkeys=['Qry','Hit'],name='gablam',expect=False)
            ldb = selfdb.addTable(mainkeys=['Qry','Hit','AlnNum'],name='local',expect=False)
            if hdb and ldb and gdb: return selfdb
            else: return None
        except: self.errorLog('%s.selfQAssemble error' % self.prog()); return None
#########################################################################################################################
    def disMatrixOut(self,gtables=[]): ### Output BLAST GABLAM distance matrix and Tree.
        '''
        Output BLAST GABLAM distance matrix and Tree.
        >> gtables:list [] = List of GABLAM summary tables to combine into dismatrix and tree.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gablam = rje_dismatrix.DisMatrix(self.log,self.cmd_list+['basefile=%s' % self.baseFile()])
            gablam.setStr({'Name':'%s GABLAM' % self.baseFile()})
            diskey = 'Qry_AlnID'
            namedict = {}
            spcode = None

            ### ~ [1] ~ Generate and output DisMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for table in gtables:
                #self.debug(self.printLog('#DIS','Adding distances from %s' % table.name()))
                if not table: raise ValueError('Problem with PacBio dismatrix tables: %s' % gtables)
                table.dataFormat({diskey:'float','Hit_AlnID':'float'})
                for entry in table.entries():
                    if not spcode: spcode = string.split(entry['Qry'],'_')[1]
                    gablam.addDis(entry['Qry'],entry['Hit'],100.0-entry[diskey])
                    # Ref vs Assembly needs to be added both ways round
                    if not gablam.getDis(entry['Hit'],entry['Qry']): gablam.addDis(entry['Hit'],entry['Qry'],100.0-entry['Hit_AlnID'])
                    #self.debug('%s vs %s (%s & %s)' % (entry['Qry'],entry['Hit'],gablam.getDis(entry['Qry'],entry['Hit']),gablam.getDis(entry['Hit'],entry['Qry'])))
                    # Add to namedict
                    if entry['Qry'] not in namedict:
                        slen = float(entry['QryLen'])
                        if slen > 1e6: namedict[entry['Qry']] = '%s (%.2f Mb)' % (entry['Qry'],slen/1e6)
                        elif slen > 1e3: namedict[entry['Qry']] = '%s (%.2f kb)' % (entry['Qry'],slen/1e3)
                        else: namedict[entry['Qry']] = '%s (%.3f kb)' % (entry['Qry'],slen/1e3)
                        if entry['Qry'] in self.dict['QV']: namedict[entry['Qry']] = '%s; QV=%s)' % (namedict[entry['Qry']][:-1],rje.sf(self.dict['QV'][entry['Qry']],3))
                    if entry['Hit'] not in namedict:
                        slen = float(entry['HitLen'])
                        if slen > 1e6: namedict[entry['Hit']] = '%s (%.2f Mb)' % (entry['Hit'],slen/1e6)
                        elif slen > 1e3: namedict[entry['Hit']] = '%s (%.2f kb)' % (entry['Hit'],slen/1e3)
                        else: namedict[entry['Hit']] = '%s (%.3f kb)' % (entry['Hit'],slen/1e3)
                        if entry['Hit'] in self.dict['QV']: namedict[entry['Hit']] = '%s; QV=%s)' % (namedict[entry['Hit']][:-1],rje.sf(self.dict['QV'][entry['Hit']],3))
            gablam.saveMatrix(filename='%s.dismatrix.csv' % self.baseFile())

            ### ~ [2] ~ Make DisMatrix symmetrical on mean difference and generate tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gablam.forceSymmetry(method='mean',missing=100.0)
            self.setNum({'MST':gablam.MST()})
            ## ~ [2a] ~ Generate tree outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upgma = rje_tree.Tree(self.log,['savetype=none','treeformats=nwk,text,png,r']+self.cmd_list+['autoload=F'])
            nsftree = gablam.upgma()
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            self.printLog('#TREE','Total TreeLen=%s; MST=%s' % (rje.sf(upgma.treeLen(),3),rje.sf(gablam.MST(),3)))
            for node in upgma.node:
                if node.info['Name'] in namedict: node.info['Name'] = namedict[node.info['Name']]
            upgma.basefile(self.basefile())
            upgma.info['GroupSpecies'] = spcode
            self.printLog('#SPEC','Looking for tree duplications based on %s' % spcode)
            upgma.opt['QueryGroup'] = True
            rje_tree_group._dupGroup(upgma,useseq=False) #.findDuplicationsNoSeq(spcode)
            upgma.saveTrees()
            #gablam.savePNG(upgma)
        except: self.errorLog('Major problem with %s.disMatrixOut' % self)
#########################################################################################################################
    ### <4> ### PAGSAT Assessment Methods                                                                               #
#########################################################################################################################
    def assessment(self,qassemble=True):   ### Generate GABLAM analyses and summarise assembly assessment.
        '''Generate GABLAM analyses and summarise assembly assessment.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Assessment ~~~~~ ##')
            db = self.db()
            if not self.force():
                wanted = ['Summary.tdt',    # Summary delimited text file
                          'png']            # Summary tree of assembly vs reference chromosomes
                complete = True
                for wext in wanted:
                    wfile = '%s.%s' % (self.baseFile(),wext)
                    self.printLog('#CHECK','%s: %s.' % (wfile,os.path.exists(wfile)))
                    complete = complete and os.path.exists(wfile)
                if complete:
                    self.printLog('#SKIP','Assessment run found (force=F).')
                    if self.getBool('RGraphics'): return self.rGraphics()
                    return True
            minloclen = self.getInt('MinLocLen')    #!# Add to basefile #!#
            minid = self.getNum('MinLocID') / 100.0
            ## ~ [0a] Load Sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Read in reference and assembly to seqlists
            self.obj['RefSeq'] = refseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('RefGenome'),'autoload=T','seqmode=file'])
            self.obj['Assembly'] = assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=file'])
            ## ~ [0b] Reference All-by-all (including self-assessment benchmark) ~~~~~~~~~~~~~~~~~~ ##
            # NOTE: Cannot do a combined reference+assembly self-GABLAM as QAssemble stats would be messed up.
            # >> This is actually better in some ways as reference self-search can be re-used!
            selfdb = self.selfQAssemble()
            if not selfdb: raise ValueError('Reference self-GABLAM failure.')
            for table in ['hitsum','local','gablam']:
                tdb = selfdb.getTable(table)
                tdb.setStr({'Name':'Self.%s' % table})
                tdb.baseFile(self.baseFile())
                db.list['Tables'].append(tdb)
            ## ~ [0c] Assembly All-by-all for dismatrix and tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aselfdb = self.selfQAssemble(self.getStr('Assembly'))
            if not aselfdb: raise ValueError('Assembly self-GABLAM failure.')
            ## ~ [0d] Load/Generate GABLAM Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #gdir = '%s.GABLAM/' % self.baseFile()
            gdir = self.getStr('GABLAMDir')
            gbase = self.fileBase('GABLAM','Cut')
            db.baseFile(gbase)
            # Look to add QAssemble GABLAM tables: run qAssembleGABLAM() if force/missing.
            self.qAssembleGABLAM('Reference')
            qrefdb = self.db(table='Reference.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
            if not qrefdb: raise IOError('Failed to generate %s.Reference.hitsum.tdt' % self.db().baseFile())
            gdb = self.db(table='Reference.gablam',add=True,forcecheck=True,mainkeys=['Qry','Hit'])
            # Read in local genome alignment results
            locdb = self.db(table='Reference.local',add=True,forcecheck=True,mainkeys=['Qry','Hit','AlnNum'])
            #if qassemble and not locdb: self.qAssembleGABLAM('Assembly'); return self.assessment(False)
            #if qassemble and not (qrefdb and gdb): self.qAssembleGABLAM('Assembly'); return self.assessment(False)
            self.qAssembleGABLAM('Assembly')
            qassdb = self.db(table='Assembly.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
            #if qassemble and not qassdb: self.qAssembleGABLAM('Reciprocal'); return self.assessment(False)
            #qassdb.setStr({'Name':'Assembly.hitsum'})

            # Read in gene and protein results
            ftbase = self.fileBase('GABLAM','Base')
            db.baseFile(ftbase)
            #i# NOTE: Reciprocal Gene and Protein searches have been disabled for now as not that useful.
            genedb = protdb = grepdb = prepdb = None
            if self.getBool('GeneSummary'):
                self.qAssembleGABLAM('Genes')
                genedb = self.db(table='Genes.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not genedb: self.qAssembleGABLAM('Genes'); return self.assessment(False)
                self.qAssembleGABLAM('Genes.Fragments')    # Used only for TopHits Analysis
                #self.qAssembleGABLAM('Genes.Reciprocal')
                #grepdb = self.db(table='Genes.Reciprocal.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not grepdb: self.qAssembleGABLAM('Genes.Reciprocal'); return self.assessment(False)
            if self.getBool('ProtSummary'):
                self.qAssembleGABLAM('Proteins')            # Used only for TopHits Analysis
                protdb = self.db(table='Proteins.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not protdb: self.qAssembleGABLAM('Proteins'); return self.assessment(False)
                self.qAssembleGABLAM('Proteins.Fragments')  # Used only for TopHits Analysis
                #self.qAssembleGABLAM('Proteins.Reciprocal')
                #prepdb = self.db(table='Proteins.Reciprocal.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not prepdb: self.qAssembleGABLAM('Proteins.Reciprocal'); return self.assessment(False)
            db.baseFile(self.baseFile())
            for table in db.tables(): table.baseFile(self.baseFile())

            # Plot both chr-contig and chr-chr in same table
            # Generate stats on duplication (contigs > self) and fragmentation (no local hits?) for each chromosome
            # ? Can we also map genes and proteins to chromosomes and make the summary per chromosome as well as total?
            # ? Should we also perform the gene/protein searches against both Reference and Assembly for benchmark?


            ### ~ [1] Order assembly contigs and reference by matches (based on biggest local alignments) ~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSAT Contig Ordering ~~~~~ ##')
            locdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            locdb.dropField('Positives')
            locdb.makeField('Identity/Length','Local')
            locdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor local hits')
            # Create lists of contigs and chromosomes from seqlists
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            #self.debug(chrdict)
            chrom = refseq.names()
            contigs = assembly.names()
            mapping = []
            # Work through local alignments in identity (count) order for primary pairings
            for lentry in locdb.sortedEntries('Identity',reverse=True):
                if not contigs: break
                if lentry['Hit'] not in contigs: continue
                # If contig OR chromosome still in list:
                # a. Add to assembly list: (chr,start,end,contig,start,end)
                mapping.append((lentry['Qry'],lentry['QryStart'],lentry['QryEnd'],lentry['Hit'],lentry['SbjStart'],lentry['SbjEnd']))
                # b. Remove chr and contig from lists
                if lentry['Qry'] in chrom: chrom.remove(lentry['Qry'])
                if lentry['Hit'] in contigs: contigs.remove(lentry['Hit'])
            if chrom: self.printLog('#CHROM','%d reference chromosomes without primary contig hits: %s' % (len(chrom),string.join(chrom,', ')))
            else: self.printLog('#CHROM','No reference chromosomes without primary contig hits.')
            if contigs: self.printLog('#CHROM','%d assembly contigs without primary reference hits: %s' % (len(contigs),string.join(contigs,', ')))
            else: self.printLog('#CHROM','No assembly contigs without primary reference hits.')
            # Sort assembly list
            mapping.sort()
            c2cmap = mapping[0:]
            # Cycle through chromosomes in order and output matching contigs to *.ordered.fas
            ordfile = '%s.ordered.fas' % self.baseFile()
            if self.force() or not rje.exists(ordfile):
                ordered = []
                for chrom in refseq.names():
                    for pair in mapping[0:]:
                        if pair[0] != chrom: continue
                        ordered.append(seqdict[pair[3]])
                        mapping.remove(pair)
                # Add (and log) any unmatched contigs to *.ordered.fas
                for seq in contigs: ordered.append(seqdict[seq])
                assembly.saveSeq(ordered,seqfile=ordfile)

            ### ~ [2] Chromosome-Contig Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb = self.db(table='ChromAlign',add=True,forcecheck=True,mainkeys=['Chrom'])
            #!# Need to add a diploid=T mode: duplicate queries to W&C copies, run, then recombine #!#
            if self.getBool('ChromAlign') and not adb and self.getBool('Diploid'): self.chromAlign(c2cmap,True)
            if self.getBool('ChromAlign') and not adb: adb = self.chromAlign(c2cmap)

            ### ~ [3] Generate Summary Tables, Charts and Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Summary Tables ~~~~~ ##')
            ## ~ [1a] Generate Distance Matrix and Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Think about how best to make the distance matrix and tree. Add an assembly self-QAssemble and combine?
            self.disMatrixOut([gdb,self.db('Self.gablam'),aselfdb.getTable('gablam')])
            ## ~ [1b] Generate Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Add summary columns etc.
            sumtab = []
            for table in [self.db('Self.hitsum'),qrefdb,qassdb,genedb,protdb,grepdb,prepdb]:
                if not table: continue
                self.printLog('#TABLE','### ~~~~~~~~ %s ~~~~~~~~ ###' % table.name())
                table.setStr({'Name':string.replace(table.name(),'hitsum','Coverage')})
                table.dataFormat({'HitNum':'int','Positives':'int','EVal':'num','Length':'int','Coverage':'int','Identity':'int'})
                table.makeField('Length-Coverage','Missing',evalue=0)
                table.makeField('Coverage-Identity','Errors',evalue=0)
                table.dataFormat({'Missing':'int','Errors':'int'})
                table.addField('Perfect',evalue=0)
                for entry in table.entries():
                    if not entry['Missing'] and not entry['Errors']: entry['Perfect'] = 1
                table.saveToFile()
                sumtab.append(db.copyTable(table,string.replace(table.name(),'.Coverage','')))
            # Output summary statistics
            covdb = None
            if adb: fullsum = sumtab + db.splitTable(adb,'Type')
            else: fullsum = sumtab
            for table in fullsum:
                if table not in sumtab:
                    table.dropFields(['Insertions','Deletions']); table.renameField('Chrom','Qry')
                    table.setStr({'Name':string.split(table.name(),'_')[-1]+'Align'})
                else: table.dropFields(['MaxScore','EVal'])
                table.addField('N',evalue=1)
                table.addField('Summary',evalue=table.name())
                table.compress(['Summary'],default='sum')
                if table not in sumtab: table.list['Fields'] = covdb.fields()
                else: table.list['Fields'] = ['Summary'] + table.fields()[:-1]
                table.dropField('Qry')
                if covdb: db.mergeTables(covdb,table)
                else: covdb = table; table.setStr({'Name':'Summary'})
            covdb.saveToFile()
            # Generate Excel file
            if self.dev() and False:
                try:
                    import XlsxWriter
                    #!# Add code here
                except: self.errorLog('XlsxWriter not on system: no Excel output')



            ### ~ [4] Generate coverage plotting data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Coverage Data ~~~~~ ##')
            ## ~ [4a] Coverage of Reference Chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chrdb = db.addEmptyTable('covplot.chrom',['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom'],['Chrom','Pos'])
            selfdb = self.db('Self.local')  # This is the chromosome versus chromosome plot
            selfdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            selfdb.dropField('Positives')
            selfdb.makeField('Identity/Length','Local')
            selfdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor (self) local hits')
            for chrom in locdb.index('Qry'):    # Should be the same for both.
                chromlen = refseq.seqLen(chrdict[chrom])
                posdict = {chromlen:{}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for pos in locdb.indexDataList('Qry',chrom,'QryStart') + selfdb.indexDataList('Qry',chrom,'QryStart'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos-1) not in posdict: posdict[pos-1] = {}
                for pos in locdb.indexDataList('Qry',chrom,'QryEnd') + selfdb.indexDataList('Qry',chrom,'QryEnd'):
                    if pos not in posdict: posdict[pos] = {}
                    if pos < chromlen and (pos+1) not in posdict: posdict[pos+1] = {}
                poslist = rje.sortKeys(posdict)
                for lentry in locdb.indexEntries('Qry',chrom) + selfdb.indexEntries('Qry',chrom):
                    hit = lentry['Hit']
                    for pos in poslist:
                        if pos < lentry['QryStart']: continue
                        if pos > lentry['QryEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['SbjStart'])
                for pos in poslist:
                    pentry = {'Chrom':chrom,'Pos':pos,
                              'ContigNum':0,'HitNum':0,'Class':'N','Contigs':[],
                              'ChromHit':0,'ChromNum':0,'RefChrom':[]}
                    for contig in posdict[pos]:
                        if contig in seqdict:   # Contig
                            pentry['ContigNum'] += 1
                            pentry['HitNum'] += len(posdict[pos][contig])
                            pentry['Contigs'].append(string.split(contig,'_')[0])
                        else:
                            pentry['ChromNum'] += 1
                            pentry['ChromHit'] += len(posdict[pos][contig])
                            pentry['RefChrom'].append(string.split(contig,'_')[0])
                    pentry['Contigs'].sort()
                    pentry['Contigs'] = string.join(pentry['Contigs'],';')
                    pentry['RefChrom'].sort()
                    pentry['RefChrom'] = string.join(pentry['RefChrom'],';')
                    # Class
                    if pentry['HitNum'] == 1: pentry['Class'] = 'U'
                    elif pentry['ContigNum'] == 1: pentry['Class'] = 'C'
                    elif pentry['ContigNum'] == pentry['HitNum'] == 2: pentry['Class'] = 'D'
                    elif pentry['ContigNum'] > 1: pentry['Class'] = 'M'
                    chrdb.addEntry(pentry)
            chrdb.saveToFile()
            ## ~ [4b] Coverage of Assembled Chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Rationalise titles of these tables?
            # Make sure that SbjStart < SbjEnd. Store direction for 3c below.
            locdb.addField('Dirn',evalue='Fwd')
            for lentry in locdb.entries():
                if lentry['SbjStart'] > lentry['SbjEnd']:
                    (lentry['SbjStart'],lentry['SbjEnd']) = (lentry['SbjEnd'],lentry['SbjStart'])
                    lentry['Dirn'] = 'Rev'
            # This is now basically the same but swapping chromosomes and contigs
            chrdb = db.addEmptyTable('covplot.contig',['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom'],['Chrom','Pos'])
            chrdict = refseq.seqNameDic()   # Chromosomes to sequences
            ctgdict = seqdict               # Contigs to sequences
            selfdb = aselfdb.getTable('local')  # Assembly versus self
            selfdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            selfdb.dropField('Positives')
            selfdb.makeField('Identity/Length','Local')
            selfdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor (self) local hits')
            for chrom in locdb.index('Hit'):    # Should be the same for both.
                posdict = {assembly.seqLen(ctgdict[chrom]):{}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for pos in locdb.indexDataList('Hit',chrom,'SbjStart') + selfdb.indexDataList('Qry',chrom,'QryStart'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos-1) not in posdict: posdict[pos-1] = {}
                for pos in locdb.indexDataList('Hit',chrom,'SbjEnd') + selfdb.indexDataList('Qry',chrom,'QryEnd'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos+1) not in posdict: posdict[pos+1] = {}
                poslist = rje.sortKeys(posdict)
                for lentry in locdb.indexEntries('Hit',chrom):
                    hit = lentry['Qry']
                    for pos in poslist:
                        if pos < lentry['SbjStart']: continue
                        if pos > lentry['SbjEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['QryStart'])
                for lentry in selfdb.indexEntries('Qry',chrom):
                    hit = lentry['Hit']
                    for pos in poslist:
                        if pos < lentry['QryStart']: continue
                        if pos > lentry['QryEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['SbjStart'])
                for pos in poslist:
                    pentry = {'Chrom':chrom,'Pos':pos,
                              'ContigNum':0,'HitNum':0,'Class':'N','Contigs':[],
                              'ChromHit':0,'ChromNum':0,'RefChrom':[]}
                    for contig in posdict[pos]:
                        if contig in chrdict:   # Chromosome
                            pentry['ContigNum'] += 1
                            pentry['HitNum'] += len(posdict[pos][contig])
                            pentry['Contigs'].append(string.split(contig,'_')[0])
                        else:
                            pentry['ChromNum'] += 1
                            pentry['ChromHit'] += len(posdict[pos][contig])
                            pentry['RefChrom'].append(string.split(contig,'_')[0])
                    pentry['Contigs'].sort()
                    pentry['Contigs'] = string.join(pentry['Contigs'],';')
                    pentry['RefChrom'].sort()
                    pentry['RefChrom'] = string.join(pentry['RefChrom'],';')
                    # Class
                    if pentry['HitNum'] == 1: pentry['Class'] = 'U'
                    elif pentry['ContigNum'] == 1: pentry['Class'] = 'C'
                    elif pentry['ContigNum'] == pentry['HitNum'] == 2: pentry['Class'] = 'D'
                    elif pentry['ContigNum'] > 1: pentry['Class'] = 'M'
                    chrdb.addEntry(pentry)
            chrdb.saveToFile()
            ## ~ [4c] Directional coverage of Assembly Contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ctgdb = db.addEmptyTable('dirplot.contig',['Contig','Pos','FwdHit','FwdChrom','RevHit','RevChrom','Chroms','Class'],['Contig','Pos'])
            for chrom in locdb.index('Hit'):
                posdict = {'Fwd':{assembly.seqLen(seqdict[chrom]):{}},'Rev':{assembly.seqLen(seqdict[chrom]):{}}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for lentry in locdb.indexEntries('Hit',chrom):
                    for dirn in ['Fwd','Rev']:
                        pos = lentry['SbjStart']
                        if pos not in posdict[dirn]: posdict[dirn][pos] = {}
                        if (pos-1) not in posdict[dirn]: posdict[dirn][pos-1] = {}
                        pos = lentry['SbjEnd']
                        if pos not in posdict[dirn]: posdict[dirn][pos] = {}
                        if (pos+1) not in posdict[dirn]: posdict[dirn][pos+1] = {}
                poslist = rje.sortUnique(rje.sortKeys(posdict['Fwd']) + rje.sortKeys(posdict['Rev']),num=True)
                for lentry in locdb.indexEntries('Hit',chrom):
                    dirn = lentry['Dirn']
                    hit = lentry['Qry']
                    for pos in poslist:
                        if pos < lentry['SbjStart']: continue
                        if pos > lentry['SbjEnd']: break    # Next entry
                        if hit not in posdict[dirn][pos]: posdict[dirn][pos][hit] = []
                        posdict[dirn][pos][hit].append(lentry['QryStart'])
                for pos in poslist:
                    pentry = {'Contig':chrom,'Pos':pos,'Class':'N','Chroms':[]}
                    for dirn in ['Fwd','Rev']:
                        pentry['%sHit' % dirn] = 0; pentry['%sChrom' % dirn] = 0
                        if pos in posdict[dirn]:
                            for contig in posdict[dirn][pos]:
                                pentry['%sHit' % dirn] += len(posdict[dirn][pos][contig])
                                pentry['%sChrom' % dirn] += 1
                                pchrom = string.split(contig,'_')[0]
                                if pchrom not in pentry['Chroms']: pentry['Chroms'].append(pchrom)
                    pentry['Chroms'].sort()
                    # Class
                    if (pentry['FwdHit']+pentry['RevHit']) == 1: pentry['Class'] = 'U'
                    elif len(pentry['Chroms']) == 1: pentry['Class'] = 'C'
                    elif len(pentry['Chroms']) > 1: pentry['Class'] = 'M'
                    pentry['Chroms'] = string.join(pentry['Chroms'],';')
                    ctgdb.addEntry(pentry)
            ctgdb.saveToFile()

            ### ~ [6] Generate additional summary of Gene and Protein Top Hits using Ref features ~~~~~~~~~~~~~~~~~~~ ###
            self.topHits()

            ### ~ [5] Generate detailed chromosome to contig mappings from local alignments ~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('RGraphics'): return self.rGraphics()
        except: self.errorLog('%s.assessment error' % self.prog())
#########################################################################################################################
    def chromAlign(self,c2cmap,diploid=False): ### Generates PAGSAT Chrom Alignment tables
        '''Returns whether *Plots/*.summary.png has been made.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSAT Chromosome Alignments (Diploid: %s) ~~~~~ ##' % diploid)
            db = self.db()
            alndir = rje.makePath('%s.ALN/' % self.baseFile())
            rje.mkDir(self,alndir)
            minloclen = self.getInt('MinLocLen')
            minid = self.getNum('MinLocID') / 100.0
            refseq = self.obj['RefSeq']
            assembly = self.obj['Assembly']
            ### ~ [1] Reference vs Assembly local BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Read BLAST into local database table including alignments
            blastfile = '%s.blast' % self.fileBase('GABLAM','Base','Reference')
            # Unzip file if required
            blastgz = rje.exists('%s.gz' % blastfile)
            if blastgz:
                os.system('gunzip %s.gz' % blastfile)
                self.printLog('#GUNZIP','%s unzipped.' % blastfile)
            if not rje.exists(blastfile): raise IOError('#BLAST','%s not found.' % (blastfile))
            # Better to re-read in the original BLAST
            blast = rje_blast.blastObj(self.log,self.cmd_list+['blastp=blastn'])
            blast.readBLAST(resfile=blastfile,clear=False,gablam=False,unlink=False,local=True,screen=True,log=False,keepaln=True)
            # (Re)zip BLAST file if required
            if blastgz:
                os.system('gzip %s' % blastfile)
                self.printLog('#GZIP','%s (re)zipped.' % blastfile)
            ## Use data from:
            bdb = blast.db('Local')
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            bdb.dataFormat({'Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'})
            bdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor local hits')
            # If diploid mode, duplicate queries here >> -A and -B
            if self.getBool('Diploid') and diploid:
                for entry in bdb.entries():
                    newentry = rje.combineDict({},entry)
                    newentry['Query'] = '%s-B' % entry['Query']
                    bdb.addEntry(newentry)
                    entry['Query'] = '%s-A' % entry['Query']
                bdb.remakeKeys()
            ### ~ [2] Chromosome-Contig Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bentries = bdb.sortedEntries('Identity',reverse=True)   # List of all entries (sorted) to process
            aentries = {}                                           # Dictionary of final entries for alignments
            alignpos = {}; ax = 0   # Dictionary of {chr/ctg:[(start,stop) list of positions included in local aln]}
            ## Process local alignments, keeping good ones and modifying the rest
            while bentries:
                entry = bentries.pop(0)     # This is best remaining hit
                ax += 1
                self.progLog('\r#LOCALN','Processing local alignments: %s -> %s' % (rje.iLen(bentries),rje.iLen(aentries)))
                entry['Qry'] = entry['Query']; entry['Sbj'] = entry['Hit']
                self.bugPrint('%s & %s' % (entry['Query'],entry['Sbj']))
                qh = 'Qry'
                if entry[qh] not in aentries: aentries[entry[qh]] = {}
                aentries[entry[qh]][entry['%sStart' % qh]] = entry
                # Update alignpos
                for qh in ['Qry','Sbj']:
                    if entry[qh] not in alignpos: alignpos[entry[qh]] = []
                    alignpos[entry[qh]].append((min(entry['%sStart' % qh],entry['%sEnd' % qh]),max(entry['%sStart' % qh],entry['%sEnd' % qh])))
                    alignpos[entry[qh]].sort()
                    i = 1
                    while i < len(alignpos[entry[qh]]):
                        if alignpos[entry[qh]][i][0] <= alignpos[entry[qh]][i-1][1]+1:
                            alignpos[entry[qh]][i-1] = (alignpos[entry[qh]][i-1][0],max(alignpos[entry[qh]][i-1][1],alignpos[entry[qh]][i][1]))
                            alignpos[entry[qh]].pop(i)
                        else: i += 1
                # Adjust/Filter remaining entries
                ex = 0
                while ex < len(bentries):
                    aentry = bentries[ex]
                    # Skip or delete if Query/Hit pair do not match
                    if entry['Query'] != aentry['Query'] and entry['Hit'] != aentry['Hit']: ex += 1; continue
                    if entry['Query'] != aentry['Query'] and entry['Hit'] == aentry['Hit']: bentries.pop(ex); continue
                    # Check for overlapping regions in Query
                    #??# What if in middle! :-(
                    #!# Need to fix this! #!#
                    dump = False
                    for (qstart,qend) in alignpos[entry['Query']]:
                        if qend < aentry['QryStart'] or qstart > aentry['QryEnd']: continue
                        if qstart <= aentry['QryStart'] and qend >= aentry['QryEnd']: dump = True; break
                        # Trim according to query positions
                        apos = 0                        # Position in alignment
                        qpos = aentry['QryStart'] - 1   # Position in query (prior to apos)
                        sdir = {True:1,False:-1}[aentry['SbjStart']<aentry['SbjEnd']]
                        spos = aentry['SbjStart'] - sdir  # Position in subject (prior to apos)
                        trimstart = aentry['QryEnd'] >= qend >= aentry['QryStart']  # and aentry['QryEnd'] >= qstart
                        trimend = aentry['QryStart'] <= qstart <= aentry['QryEnd'] #and aentry['QryStart'] <= qend
                        if (trimstart and trimend): self.bugPrint('Q: (%s,%s) vs (%s,%s)!' % (qstart,qend,aentry['QryStart'],aentry['QryEnd']))
                        while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                            if aentry['QrySeq'][apos] != '-': qpos += 1
                            if aentry['SbjSeq'][apos] != '-': spos += sdir
                            if trimstart and qpos > qend:   # Trim!
                                for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                aentry['QryStart'] = qpos
                                aentry['SbjStart'] = spos
                                apos = 1; trimstart = False
                            elif trimend and qpos == qstart - 1:
                                for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                aentry['QryEnd'] = qpos
                                aentry['SbjEnd'] = spos
                                break
                            else: apos += 1
                    if dump: bentries.pop(ex); continue
                    if entry['Query'] == aentry['Query'] and entry['Hit'] != aentry['Hit']: ex += 1; continue   # Can have 2+ contigs
                    # Check for overlapping regions in Hit
                    for (hstart,hend) in alignpos[entry['Hit']]:
                        if hend < min(aentry['SbjStart'],aentry['SbjEnd']) or hstart > max(aentry['SbjStart'],aentry['SbjEnd']): continue
                        if hstart <= min(aentry['SbjStart'],aentry['SbjEnd']) and hend >= max(aentry['SbjStart'],aentry['SbjEnd']): dump = True; continue
                        # Trim according to subject positions
                        apos = 0                        # Position in alignment
                        qpos = aentry['QryStart'] - 1   # Position in query (prior to apos)
                        sdir = {True:1,False:-1}[aentry['SbjStart']<aentry['SbjEnd']]
                        spos = aentry['SbjStart'] - sdir  # Position in subject (prior to apos)
                        if sdir > 0:    # Fwd hit
                            trimstart = aentry['SbjEnd'] >= hend >= aentry['SbjStart'] #and aentry['SbjEnd'] >= hstart
                            trimend = aentry['SbjStart'] <= hstart <= aentry['SbjEnd'] #and aentry['SbjStart'] <= hend
                            #if (trimstart and trimend): self.warnLog('(%s,%s) vs %s!' % (hstart,hend,aentry))
                            if (trimstart and trimend): self.bugPrint('Fwd: (%s,%s) vs (%s,%s)!' % (hstart,hend,aentry['SbjStart'],aentry['SbjEnd']))
                            while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                                if aentry['QrySeq'][apos] != '-': qpos += 1
                                if aentry['SbjSeq'][apos] != '-': spos += sdir
                                if trimstart and spos > hend:   # Trim!
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                    aentry['QryStart'] = qpos
                                    aentry['SbjStart'] = spos
                                    apos = 1; trimstart = False
                                elif trimend and spos == hstart - 1:
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                    aentry['QryEnd'] = qpos
                                    aentry['SbjEnd'] = spos
                                    break
                                else: apos += 1
                        else:
                            trimstart = aentry['SbjEnd'] <= hstart <= aentry['SbjStart'] #and aentry['SbjEnd'] <= hend
                            trimend = aentry['SbjStart'] >= hend >= aentry['SbjEnd'] #and aentry['SbjStart'] >= hstart
                            #if (trimstart and trimend): self.warnLog('(%s,%s) vs %s!' % (hstart,hend,aentry))
                            if (trimstart and trimend): self.bugPrint('Bwd: (%s,%s) vs (%s,%s)!' % (hstart,hend,aentry['SbjStart'],aentry['SbjEnd']))
                            while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                                if aentry['QrySeq'][apos] != '-': qpos += 1
                                if aentry['SbjSeq'][apos] != '-': spos += sdir
                                if trimstart and spos == hstart - 1:   # Trim!
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                    aentry['QryStart'] = qpos
                                    aentry['SbjStart'] = spos
                                    apos = 1; trimstart = False
                                elif trimend and spos == hend + 1:
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                    aentry['QryEnd'] = qpos
                                    aentry['SbjEnd'] = spos
                                    break
                                else: apos += 1
                    if dump: bentries.pop(ex); continue
                    ex += 1
            #self.debug(alignpos)
            self.printLog('\r#LOCALN','Processed local alignments -> %s pairwise aln.' % (rje.iStr(ax)))
            #bdb.dict['Data'] = {}
            #for entry in aentries: bdb.addEntry(entry)
            #bdb.newKey(['Query','QryStart'])
            #self.debug(bdb.dataKeys())
            ## Generate localn files:
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            adb = self.db('ChromAlign')
            if not adb: adb = db.addEmptyTable('ChromAlign',['Chrom','Ctid','Type','Length','Identity','Coverage','Insertions','Deletions'],['Chrom','Ctid'])
            ddb = db.addEmptyTable('ChromAlignLoc', ['Query','Hit','Ctid','AlnID','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd'],['Query','Hit','Ctid','AlnID'])
            for chrom in rje.sortKeys(aentries):
                self.debug(chrom)
                if self.getBool('Diploid') and diploid: dentry = {'Query':chrom[:-2],'Ctid':chrom[-1]}
                else: dentry = {'Query':chrom,'Ctid':'H'}
                for entry in aentries[chrom].values():
                    self.debug(entry)
                    #if self.getBool('Diploid') and diploid: dentry = {'Query':chrom[:-2],'Ctid':entry['Query'][-1]}
                    #else: dentry = {'Query':chrom,'Ctid':'H'}
                    ddentry = rje.combineDict({},dentry,overwrite=False)
                    ddb.addEntry(rje.combineDict(ddentry,entry,overwrite=False))
                if self.getBool('Diploid') and diploid: chrseq = refseq.getSeq(chrdict[chrom[:-2]])
                else: chrseq = refseq.getSeq(chrdict[chrom])
                contigs = []
                poslist = rje.sortKeys(aentries[chrom])
                entry = aentries[chrom][poslist[0]]
                qname = string.split(entry['Query'],'__')[0]
                if self.getBool('Diploid') and diploid: qname += chrom[-2:]
                #self.debug(qname)
                aname = string.split(qname,'_')[0] + '_ALIGN'
                hname = string.split(qname,'_')[0] + '_' + string.split(entry['Hit'],'_')[1]
                qentry = {'Chrom':qname,'Type':'Reference','Ctid':dentry['Ctid']}; hentry = {'Chrom':hname,'Type':'Assembly','Ctid':dentry['Ctid']}
                qpos = 0    # Last position covered
                qseq = aseq = hseq = ''
                for spos in poslist:
                    entry = aentries[chrom][spos]
                    if entry['Hit'] not in contigs: contigs.append(entry['Hit'])
                    qname += ' %s(%s..%s)' % (string.split(entry['Qry'],'_')[-1],entry['QryStart'],entry['QryEnd'])
                    hname += ' %s(%s..%s)' % (string.split(entry['Hit'],'_')[-1],entry['SbjStart'],entry['SbjEnd'])
                    addseq = chrseq[1][qpos:entry['QryStart']]
                    if addseq:
                        qseq += addseq
                        aseq += '-' * len(addseq)
                        hseq += '.' * len(addseq)
                    if qseq: qseq += '.....'; aseq += 'XXXXX'; hseq += '.....'  # Small spacer between alignments
                    qseq += entry['QrySeq']
                    aseq += entry['AlnSeq']
                    hseq += entry['SbjSeq']
                    qpos = entry['QryEnd']
                addseq = chrseq[1][qpos:]
                if addseq:
                    qseq += addseq
                    aseq += '-' * len(addseq)
                    hseq += '.' * len(addseq)
                aseq = string.replace(aseq,' ','x')
                qentry['Length'] = len(chrseq[1])
                hentry['Length'] = 0
                for ctg in contigs: hentry['Length'] += assembly.seqLen(seqdict[ctg])
                qentry['Identity'] = hentry['Identity'] = string.count(aseq,'|')
                qentry['Insertions'] = hentry['Deletions'] = string.count(hseq,'-')
                qentry['Deletions'] = hentry['Insertions'] = string.count(qseq,'-')
                qentry['Coverage'] = hentry['Coverage'] = string.count(aseq,'x') + qentry['Identity'] - qentry['Insertions'] - qentry['Deletions']
                if not diploid:
                    ALNFAS = open('%s%s.aln.fas' % (alndir,chrom),'w')
                    ALNFAS.write('>%s\n%s\n' % (qname,qseq))
                    ALNFAS.write('>%s\n%s\n' % (aname,aseq))
                    ALNFAS.write('>%s\n%s\n' % (hname,hseq))
                    ALNFAS.close()
                    adb.addEntry(qentry)
                if diploid == self.getBool('Diploid'):
                    adb.addEntry(hentry)
            if diploid == self.getBool('Diploid'):
                ddb.saveToFile()
            if not diploid:
                adb.makeField('Length-Coverage','Missing',evalue=0)
                adb.makeField('Coverage-Identity','Errors',evalue=0)
                adb.addField('Perfect',evalue=0)
                for entry in adb.entries():
                    if not entry['Missing'] and not entry['Errors']: entry['Perfect'] = 1
                adb.compress(['Chrom'],default='sum')
                adb.dropField('Ctid')
                adb.saveToFile()
            #else: return adb

            # Build up and combine local alignments
            # >> Build up best -> Worst (BitScore), keeping if new sequences involved and trimming ends where appropriate
            # >> Fill in the gaps by splitting the shortest sequence in the middle and inserting gaps
            # Q. How to deal with inversions? RevComp and lower case?

            # Create lists of contigs and chromosomes from seqlists
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            # Cycle through chromosomes in order and output matching contigs to *.ordered.fas
            for chrom in rje.sortKeys(chrdict): # refseq.names():
                #self.debug('%s -> %s' % (chrom,chrdict[chrom]))
                ALNFAS = open('%s%s.fas' % (alndir,chrom),'w')
                ALNFAS.write(refseq.fasta(chrdict[chrom]))
                for pair in c2cmap[0:]:
                    if pair[0] != chrom: continue
                    c2cmap.remove(pair)
                    if pair[-2] > pair[-1]:     # reverse!
                        (name,sequence) = assembly.getSeq(seqdict[pair[3]])
                        name += 'Reverse Complement'
                        #self.bugPrint(sequence[:100])
                        sequence = rje_sequence.reverseComplement(sequence)
                        #self.deBug(sequence[-100:])
                        ALNFAS.write('>%s\n%s\n' % (name,sequence))
                    else: ALNFAS.write(assembly.fasta(seqdict[pair[3]]))
                ALNFAS.close()
                #!# The following takes a long time and should be replaced with a better genome alignment program! (Exonerate?)
                #rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s%s.fas' % (alndir,chrom),'autoload=T','autofilter=F','dna=T','align=F']).align(outfile='%s%s.aln.fas' % (alndir,chrom))
            return adb
        except: self.errorLog('%s.chromAlign error' % self.prog()); raise
#########################################################################################################################
    def rGraphics(self):    ### Generates PAGSAT R Graphics
        '''Returns whether *Plots/*.summary.png has been made.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT R Graphics ~~~~~ ##')
            basename = self.baseFile(strip_path=True)
            rsfile ='%s.Plots/%s.summary.png' % (self.baseFile(),basename) # Summary reference figure.
            if not self.force() and os.path.exists(rsfile):
                self.printLog('#CHECK','%s: Found. (Force=F)' % rsfile)
                return True
            rcmd = '%s --no-restore --no-save --args "pagsat" "%s"' % (self.getStr('RPath'),self.baseFile())
            rcmd += ' "refbase=%s"' % rje.baseFile(self.getStr('RefGenome'))
            rcmd += ' "minloclen=%d"' % self.getInt('MinLocLen')
            if not self.getBool('Diploid'):
                rcmd += ' "diploid=F"'
            rdir = '%slibraries/r/' % slimsuitepath
            rcmd += ' "rdir=%s" < "%srje.r" > "%s.r.tmp.txt"' % (rdir,rdir,self.baseFile())
            self.printLog('#RPNG',rcmd)
            problems = os.popen(rcmd).read()
            if problems:
                for ptxt in problems: self.warnLog(ptxt)
            # Optional cleanup of *.r.tmp.txt ?
            return os.path.exists(rsfile)
        except: self.errorLog('%s.rGraphics error' % self.prog()); return False
#########################################################################################################################
    def topHits(self):    ### Generates output for Gene/Protein TopHits analysis
        '''Returns the Reference Features table, if given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT TopHits/Synteny Analysis ~~~~~ ##')
            db = self.db()
            if not self.force() and rje.exists('%s.Genes.TopHits.tdt' % db.baseFile()) and rje.exists('%s.Proteins.TopHits.tdt' % db.baseFile()):
                self.printLog('#HITS','Genes.TopHits and Proteins.TopHits tables found (force=F).'); return True
            tophitbuffer = self.getNum('TopHitBuffer')
            if tophitbuffer < 0: self.warnLog('Cannot have TopHitBuffer < 0.0!'); tophitbuffer = 0.0
            refseq = self.obj['RefSeq']
            ftdict = db.splitTable(self.ftdb(),'feature',asdict=True,keepfield=False,splitchar=None,values=['gene','CDS'])    # Reference features tables
            ftdict['Genes'] = ftdict.pop('gene')
            #self.debug(ftdict['Genes'].entries()[:10])
            ftdict['Proteins'] = ftdict.pop('CDS')
            #self.debug(ftdict['Proteins'].index('note').keys()[:10])
            acc2chr = {}
            for seq in refseq.seqs(): acc2chr[string.split(refseq.seqAcc(seq),'.')[0]] = refseq.seqGene(seq)
            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gtype in ['Genes','Proteins']:
                ## ~ [1a] Load GABLAM data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.getBool('%sSummary' % gtype[:4]): continue
                gdb = db.addTable(filename='%s.%s.Fragments.gablam.tdt' % (self.fileBase('GABLAM','Base'),gtype),mainkeys=['Qry','Hit'],name='%s.gablam' % gtype,expect=True)
                if not gdb: self.warnLog('%s.%s.Fragments.gablam.tdt missing!' % (self.fileBase('GABLAM','Base'),gtype)); continue
                ## ~ [1b] Filter, join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #gdb.dropEntries(['Rank>1']) #!# No more! Filter on
                for entry in gdb.entries(): entry['Qry'] = string.split(entry['Qry'],'_')[-1]
                fdb = db.joinTables(name='%s.TopHits' % gtype,join=[(ftdict[gtype],'locus_tag'),(gdb,'Qry')],newkey=['locus_tag','Hit'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
                #self.debug(fdb.index('note').keys()[:10])
                gfield = 'Contig%s' % gtype[:4]     # Rename to SPXXXXX.X.[P|G]Y
                fdb.addField(gfield)
                fdb.setFields(['locus_tag','Hit',gfield,'QryLen','Qry_AlnLen','Qry_AlnID','Qry_Start','Qry_End','Hit_Start','Hit_End','locus','position','start','end','product','gene_synonym','db_xref','note'])
                fdb.dataFormat({'Qry_AlnID':'num','Qry_Start':'int'})
                fdb.addField('Chrom',after='locus')
                fdb.addField('Contig',after='Chrom')
                for entry in fdb.entries():
                    entry['Chrom'] = acc2chr[entry['locus']]
                    if not entry['Hit']: continue
                    [entry['Contig'],strain,na,entry[gfield]] = string.split(entry['Hit'],'_')    # hcq10_MBG11A__SP16495.10-011478.016127
                    cdata = string.split(entry['Hit'],'-')[-1]  # Position info
                    cdata = string.split(cdata,'.')
                    entry['Hit_Start'] = int(cdata[0])
                    entry['Hit_End'] = int(cdata[1])
                fdb.renameField('locus_tag','Gene')
                fdb.renameField('gene_synonym','Synonym')
                fdb.renameField('db_xref','XRef')
                for field in ['locus','position','start','end','product','note']: fdb.renameField(field,rje.strSentence(field))
                #self.debug(fdb.index('Note').keys()[:10])
                ## ~ [1c] Report missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nullx = 0
                #for entry in fdb.entries():
                #    if not entry['Hit']: nullx += 1
                #self.printLog('#HITS','%s of %s %s have no Assembly hit.' % (rje.iStr(nullx),rje.iStr(fdb.entryNum()),gtype))
                nullx = len(fdb.indexDataList('Hit','','Gene'))
                self.printLog('#HITS','%s of %s %s have no Assembly hit.' % (rje.iStr(nullx),rje.iLen(fdb.indexKeys('Gene')),gtype))
                ## ~ [1d] Check overlapping hits & rename (should be any, right?) ~~~~~~~~~~~~~~~~~ ##
                prevhit = (None,-1,-1)  #?# Make tuple of (ctg,start,end)?
                gx = 0; hx = len(fdb.indexKeys('Hit'))
                self.progLog('\r#HITS','%s...' % gfield)
                for ghit in fdb.indexKeys('Hit'):
                    # Compare start with prev end
                    if not ghit: continue
                    #self.bugPrint(ghit)
                    thishit = rje.matchExp('^(\S+)_\S+__(\S+)-(\d+)\.(\d+)$',ghit)
                    if prevhit and prevhit[0] == thishit[0]:
                        if int(thishit[-1]) < int(prevhit[-2]): raise ValueError('Overlapping fragment %s. Should not happen!' % ghit)
                    prevhit = thishit
                    # Rename protein/gene
                    gx += 1
                    for entry in fdb.indexEntries('Hit',ghit): entry[gfield] = '%s.%s%s' % (thishit[1],gtype[:1],rje.preZero(gx,hx))
                ## ~ [1e] Filter according to best hits (each direction) ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                besthit = {}            # Dictionary of {gene:best %ID}
                bestentries = {}        # Dictionary of {gene:[best entries]}
                self.progLog('\r#HITS','Best Hits...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: continue
                    for gid in (entry['Gene'],entry[gfield]):
                        #!# Make this a dictionary to store tied best hits! Add Best field
                        if gid not in besthit or entry['Qry_AlnID'] > besthit[gid]:
                            besthit[gid] = entry['Qry_AlnID']
                            bestentries[gid] = [entry]
                        elif entry['Qry_AlnID'] == besthit[gid]: bestentries[gid].append(entry)
                ftotx = fdb.entryNum()
                self.progLog('\r#HITS','Top Hit Buffer...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: continue
                    # Only keep matches where BOTH sequences are within tophitbuffer
                    if entry['Qry_AlnID'] < (max(besthit[entry['Gene']],besthit[entry[gfield]]) - tophitbuffer):
                        fdb.dropEntry(entry)
                self.printLog('#HITS','%s %s reduced to %s within %.2f%% of best hit.' % (rje.iStr(ftotx),gtype,rje.iStr(fdb.entryNum()),tophitbuffer))
                ## ~ [1f] Classify mappings (1:1,1:n,n:1,n:n) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fdb.addField('Mapping')
                fdb.index('Gene'); fdb.index(gfield)
                self.progLog('\r#HITS','Classify m:n mappings...')
                for entry in fdb.entries():
                    if not entry['Qry_AlnID']: entry['Mapping'] = '1:0'; continue
                    n = len(fdb.index('Gene')[entry['Gene']])
                    m = len(fdb.index(gfield)[entry[gfield]])
                    entry['Mapping'] = '%s:%s' % (m,n)
                ## ~ [1g] Identify 1:x mappings and map ref ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#HITS','Identify and map 1:x mappings...')
                for m2n in fdb.index('Mapping'):
                    if m2n[:2] != '1:': continue
                    if m2n == '1:0': continue
                    for gene in fdb.indexDataList('Mapping',m2n,'Gene'):
                        gx = 1; gtot = len(fdb.index('Gene')[gene])
                        for entry in fdb.indexEntries('Gene',gene):
                            if gtot == 1: entry[gfield] = gene      # 1:1 = easy straight mapping
                            elif entry in bestentries[gene] and len(bestentries[gene]) == 1: entry[gfield] = gene
                            elif gtot > 1: entry[gfield] = '%s|%d' % (gene,gx); gx += 1
                #fdb.indexReport('Mapping') # Report subset?
                ## ~ [1h] Synteny ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fdb.addFields(['Chr5','Chr3','Ctg5','Ctg3','Synteny'])
                synteny_updated = True; sloop = 0
                while synteny_updated:
                    synteny_updated = False
                    #i# Initially, look at direct synteny. Might want to extend at some point.
                    # Reference genes
                    synteny = {}
                    self.progLog('\r#SYN','Ref synteny...')
                    for chrom in fdb.index('Chrom'): synteny[chrom] = ['NA']
                    for entry in fdb.sortedEntries('Start'):
                        chrom = entry['Chrom']
                        if synteny[chrom][-1] != entry['Gene']: synteny[chrom].append(entry['Gene'])
                    for chrom in fdb.index('Chrom'): synteny[chrom] += ['NA']
                    # Assembly genes
                    self.progLog('\r#SYN','Assembly synteny...')
                    for chrom in fdb.index('Contig'): synteny[chrom] = ['NA']
                    for entry in fdb.sortedEntries('Hit'):
                        chrom = entry['Contig']
                        if synteny[chrom][-1] != entry[gfield]: synteny[chrom].append(entry[gfield])
                    for chrom in fdb.index('Contig'): synteny[chrom] += ['NA']
                    # Update synteny
                    self.progLog('\r#SYN','Updating synteny...')
                    for entry in fdb.entries():
                        if entry['Synteny'] in ['Full','Tandem']: continue
                        prevsyn = entry['Synteny']
                        chrom = entry['Chrom']
                        entry['Chr5'] = synteny[chrom][synteny[chrom].index(entry['Gene'])-1]
                        entry['Chr3'] = synteny[chrom][synteny[chrom].index(entry['Gene'])+1]
                        contig = entry['Contig']
                        entry['Ctg5'] = synteny[contig][synteny[contig].index(entry[gfield])-1]
                        entry['Ctg3'] = synteny[contig][synteny[contig].index(entry[gfield])+1]
                        ctg5 = string.split(entry['Ctg5'],'|')[0]
                        ctg3 = string.split(entry['Ctg3'],'|')[0]
                        if ctg5 == entry['Chr5'] and ctg3 == entry['Chr3']: entry['Synteny'] = 'Full'
                        elif ctg5 == entry['Chr3'] and ctg3 == entry['Chr5']: entry['Synteny'] = 'Full' # Reverse?
                        elif ctg5 in [entry['Chr5'],entry['Chr3']]: entry['Synteny'] = 'Partial'
                        elif ctg3 in [entry['Chr5'],entry['Chr3']]: entry['Synteny'] = 'Partial'
                        else: entry['Synteny'] = 'Orphan'
                        # Special identification of 2:1 "tandem compression" and 1:2 "tandem duplication"
                        # Special 2:1 or 1:2
                        if (entry['Mapping'].startswith('1:') or entry['Mapping'].endswith(':1')) and entry['Synteny'] == 'Partial':
                            chrneighbours = [entry['Chr5'],entry['Chr3']]
                            ctgneighbours = [ctg5,ctg3]
                            missing = rje.listDifference(chrneighbours,ctgneighbours)
                            if entry['Mapping'].endswith(':1') and missing:
                                missing = missing[0]
                                chrneighbours.append(synteny[chrom][synteny[chrom].index(missing)-1])
                                chrneighbours.append(synteny[chrom][synteny[chrom].index(missing)+1])
                                if len(rje.listIntersect(chrneighbours,ctgneighbours)) == 2: entry['Synteny'] = 'Tandem'
                            missing = rje.listDifference(ctgneighbours,chrneighbours)
                            if entry['Mapping'].startswith('1:') and missing:
                                missing = missing[0]
                                if missing not in synteny[contig]:
                                    variant = None
                                    #self.bugPrint('%s vs "%s|"' % (synteny[contig],missing))
                                    for gvar in synteny[contig]:
                                        if gvar.startswith('%s|' % missing):
                                            #self.bugPrint(gvar)
                                            if variant: variant = None; break
                                            variant = gvar
                                    missing = variant
                                    #self.deBug(missing)
                                if missing in synteny[contig]:
                                    ctgneighbours.append(synteny[contig][synteny[contig].index(missing)-1])
                                    ctgneighbours.append(synteny[contig][synteny[contig].index(missing)+1])
                                #else:
                                    #self.bugPrint(entry)
                                    #self.bugPrint(chrneighbours)
                                    #self.bugPrint(ctgneighbours)
                                    #self.bugPrint(missing)
                                    #self.debug(synteny[contig])
                                if len(rje.listIntersect(chrneighbours,ctgneighbours)) == 2: entry['Synteny'] = 'Tandem'
                        #if prevsyn != entry['Synteny'] and sloop: self.deBug('%s: %s => %s' % (entry,prevsyn,entry['Synteny']))
                        synteny_updated = synteny_updated or prevsyn != entry['Synteny']    # Change
                    # if x:1 keep Full synteny, else keep best hit
                    hitchanged = {}      # Will need to remake keys if Gene/Hit changes
                    for xn in fdb.indexKeys('Mapping'):
                        if not xn.endswith(':1'): continue
                        for gid in fdb.indexDataList('Mapping',xn,gfield):
                            syntenic = 'Full' in fdb.indexDataList(gfield,gid,'Synteny')   # Keep only syntenic
                            for entry in fdb.indexEntries(gfield,gid):
                                if syntenic and entry['Synteny'] != 'Full': fdb.dropEntry(entry)
                                elif not syntenic and besthit[gid] != entry['Qry_AlnID']: fdb.dropEntry(entry)
                                else: hitchanged[gid] = entry['Gene'] # Use (one of) best gene name
                    if hitchanged:
                        for gid in hitchanged:
                            for entry in fdb.indexEntries(gfield,gid): entry[gfield] = hitchanged[gid]
                        fdb.index(gfield,force=True)
                    # if x:y ...?
                    fdb.indexReport('Synteny')
                    if synteny_updated: sloop += 1; self.printLog('#LOOP','Loop %d complete: synteny updated.' % sloop)
                    else: sloop += 1; self.printLog('#LOOP','Convergence after synteny loop %d.' % sloop)

                ## ~ [1x] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries(): entry['Hit'] = string.split(entry['Hit'],'-')[0]   # Update entry['Hit'] for R Script
                fdb.saveToFile()


        except: self.errorLog('%s.topHits error' % self.prog())
#########################################################################################################################
    ### <5> ### PAGSAT Report Methods                                                                                   #
#########################################################################################################################
    def ftdb(self): ### Returns the Reference Features table, if given
        '''Returns the Reference Features table, if given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['Features']: return self.obj['Features']
            db = self.db()
            ### ~ [1] Load Reference Feature Table (if refgenome given) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('RefGenome'):
                ftdb = db.addTable('%s.Feature.tdt' % rje.baseFile(self.getStr('RefGenome')),mainkeys=['locus','feature','position'],name='features')
                ftdb.dataFormat({'start':'int','end':'int'})
            else: self.printLog('#REFFT','Cannot load/report/assess features without refgenome=FILE.'); ftdb = None
            self.obj['Features'] = ftdb
            return self.obj['Features']
        except: self.errorLog('%s.ftdb error' % self.prog()); return None
#########################################################################################################################
    def report(self):   ### Generates HTML reports of PAGSAT assessment
        '''Generates HTML reports of PAGSAT assessment.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basename = self.baseFile(strip_path=True)
            wanted = ['Summary.tdt',    # Summary delimited text file
                      'png',            # Summary tree of assembly vs reference chromosomes
                      'Plots/%s.summary.png' % basename] # Summary reference figure.
            for wext in wanted:
                wfile = '%s.%s' % (self.baseFile(),wext)
                if not os.path.exists(wfile):
                    self.printLog('#CHECK','%s: Missing! Will generate.' % wfile)
                    self.assessment(); break
                self.printLog('#CHECK','%s: Found.' % wfile)
            for wext in wanted:
                wfile = '%s.%s' % (self.baseFile(),wext)
                if not os.path.exists(wfile): raise IOError('Cannot find %s!' % wfile)
            ## ~ [0a] Setup HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = rje_html.HTML(self.log,self.cmd_list)
            hfile = '%s.report.html' % self.baseFile()
            ## ~ [0b] Dictionary of links and descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            desc = {'tree':'Summary tree of chromosomes vs contigs (% global identity)',
                    'summary':'Summary plot of assembly against reference chromosomes',
                    'chromalign':'Summary plot of aligned contigs against reference chromosomes',
                    'assembly':'Summary plot of reference chromosomes against assembly contigs'}





            ### ~ [1] Initial Summary and quick links to sections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hbody = ['<h1>%s Report</h1>' % basename,'']
            hbody += ['<p>','Quick Links:']
            for id in ['Tree','Summary PNG','ChromAlign PNG','Assembly PNG']:
                link = string.split(id)[0].lower()
                hbody.append('~ <a href="#%s" title="%s">%s</a>' % (link,desc[link],id))
            hbody += ['</p>','']


            hbody += ['<h2>%s Summary Table</h1>' % basename,'']
            sumtable = open('%s.Summary.tdt' % self.basefile()).read()
            hbody.append(rje_html.tableToHTML(sumtable,'\t',tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''))

            hbody += ['<a name="tree"><h2>%s Summary Tree</h1><a>' % basename,'']
            hbody.append('<a href="./%s.png"><img src="./%s.png" width="100%%" title="%s"></a>' % (basename,basename,desc['tree']))

            for png in ['summary','chromalign','assembly']:
                pngfile = '%s.Plots/%s.%s.png' % (self.baseFile(),basename,png)
                pnglink = './%s.Plots/%s.%s.png' % (basename,basename,png)
                hbody += ['<a name="%s"><h2>%s</h1><a>' % (png,desc[png]),'']
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,desc[png]))

            # SECTIONS:
            # Reference-based Summary

            self.warnLog('PAGSAT.Report() only partially implemented. Sorry!')

            # Include features overlapping:
            # (a) the missing regions
            # (b) repeated regions (where assembly > reference) in *.covplot.chrom.tdt

            '''
            head MBG8150.SP16481.hcq.sgd.srt.1000.covplot.chrom.tdt
            Chrom	Pos	HitNum	ContigNum	Contigs	Class	ChromHit	ChromNum	RefChrom
            chrIII_S288C__BK006937.2	0	0	0		N	0	0
            chrIII_S288C__BK006937.2	1	2	2	hcq1;hcq11	D	1	1	chrIII
            chrIII_S288C__BK006937.2	10	3	3	hcq1;hcq11;hcq6	M	1	1	chrIII
            chrIII_S288C__BK006937.2	11	4	4	hcq0;hcq1;hcq11;hcq6	M	1	1	chrIII
            chrIII_S288C__BK006937.2	11225	2	2	hcq16;hcq7	D	2	2	chrIII;chrXI
            chrIII_S288C__BK006937.2	11226	1	1	hcq16	U	1	1	chrIII
            chrIII_S288C__BK006937.2	114	6	6	hcq0;hcq1;hcq11;hcq4;hcq5;hcq6	M	2	2	chrIII;chrXIV
            chrIII_S288C__BK006937.2	1149	9	9	hcq0;hcq1;hcq11;hcq16;hcq2;hcq4;hcq5;hcq6;hcq7	M	2	2	chrIII;chrXIV
            chrIII_S288C__BK006937.2	115	7	7	hcq0;hcq1;hcq11;hcq4;hcq5;hcq6;hcq7	M	2	2	chrIII;chrXIV
            '''



            ### ~ [X] Output HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            HTML = open(hfile,'w')
            HTML.write(html.htmlHead(title=basename,tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
            HTML.write(string.join(hbody,'\n'))
            HTML.write(html.htmlTail(tabber=False))
            HTML.close()
            self.printLog('#HTML','HTML report output: %s' % hfile)
        except: self.errorLog('%s.report error' % self.prog())
#########################################################################################################################
    ### <6> ### PAGSAT Comparison Methods                                                                               #
#########################################################################################################################
    def compare(self):  ### Generates summary of statistics across multiple PAGSAT runs.
        '''Generates summary of statistics across multiple PAGSAT runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Compare mode (%d files) ~~~~~ ##' % len(self.list['Compare']))
            db = self.db()
            # This method essentially wants to read and combine the summary data from several runs, and then extract the
            # most useful information for assessing assessment quality, including:
            # - %coverage and %accuracy for assembly and reference.
            # - no. contigs
            # - optional gene/protein data if present
            ## ~ [0a] Load Reference Feature Table (if refgenome given) ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            repeatft = ['rRNA','mobile','LTR','centromere','telomere']   # self.list['RepeatFT']
            repeatft.sort()
            if repeatft: ftdb = self.ftdb()
            else: ftdb = None
            if ftdb:
                ftdb.dropEntriesDirect('feature',repeatft,inverse=True)
                self.printLog('#RPTFT','%s repeat features to exclude in "Uniq" outputs.' % rje.iStr(ftdb.entryNum()))
                if not ftdb.entryNum(): ftdb = None
            else:
                if not self.getStrLC('RefGenome'): self.printLog('#RPTFT','Cannot filter repeat features without refgenome=FILE.')
                if not repeatft: self.printLog('#RPTFT','Cannot filter repeat features without repeatft=LIST.')
            loc2chr = {}                # Will be a locus -> chromosome dictionary
            fragcov = self.list['FragCov']   # = [50,90,95,99] List of coverage thresholds to count (local table)
            chromcov = self.list['ChromCov'] # = [95,98,99] No. of chromosomes covered by a single contig (GABLAM table)
            fragcov.sort(); chromcov.sort()

            ### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            compfields = ['Assembly','N','%AssCov','%AssAcc','%AssAlnCov','%AssAlnAcc','Multiplicity','Parsimony','%RefCov','%RefAcc','%RefAlnCov','%RefAlnAcc','Missing','Errors','Extra','Duplicate','TreeLen']
            if ftdb: compfields += ['UniqCov','UniqCtg','UniqDup','RepeatFT']
            for chromx in chromcov: compfields.append('Chrom%d' % chromx)
            for fragx in fragcov:
                compfields.append('Frag%d' % fragx)
                if ftdb: compfields.append('UniqFrag%d' % fragx)
            compdb = db.addEmptyTable('compare',compfields,['Assembly'])   # Table of final comparison data
            if self.getBool('GeneSummary'): compdb.addFields(['%GeneCov','%GeneAcc'])#,'%GeneIntegrity'])
            if self.getBool('ProtSummary'): compdb.addFields(['%ProtCov','%ProtAcc'])#,'%ProtIntegrity'])
            for pfile in self.list['Compare']:
                #!# Can/should this whole process be moved into a function that can be run on a single dataset? #!#
                #!# Can then simply compile the datasets if found with the right headers.
                ## ~ [1a] Load and Process Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pbase = rje.baseFile(pfile,strip_path=True)
                if pbase.endswith('.Summary'): pbase = string.join(string.split(pbase,'.')[:-1],'.')     # Strip Summary
                basedir = rje.makePath('%s.PAGSAT/' % string.join(string.split(pbase,'.')[:-1],'.'))
                gabdir = rje.makePath('%s.GABLAM/' % string.join(string.split(pbase,'.')[:-1],'.'))
                self.setStr({'GABLAMDir':gabdir,'ResDir':basedir,'BaseBase':string.join(string.split(pbase,'.')[:-1],'.'),'CutBase':pbase})
                try: pdb = db.addTable(pfile,['Summary'],name=pbase,expect=True)
                except: self.errorLog('Cannot load PAGSAT Summary table "%s": check format' % pfile); continue
                pdb.dataFormat({'Length':'int','Coverage':'int','Identity':'int','Missing':'int','Errors':'int'})
                centry = {'Assembly':pbase}
                for entry in pdb.entries():
                    if entry['Summary'] == 'Reference':
                        centry['%RefCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Coverage
                        centry['%RefAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['Missing'] = entry['Missing']
                        centry['Errors'] = entry['Errors']
                        centry['Duplicate'] = 0
                    if entry['Summary'] == 'ReferenceAlign':
                        centry['%RefAlnCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Coverage
                        centry['%RefAlnAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['Duplicate'] = 0
                    if entry['Summary'] == 'Assembly':
                        centry['%AssCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Validity
                        centry['%AssAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['Multiplicity'] = float(entry['Coverage']) / pdb.data('Reference')['Coverage']
                        centry['Parsimony'] = float(entry['Length']) / pdb.data('Reference')['Coverage']
                        centry['Extra'] = entry['Missing']
                        centry['N'] = entry['N']
                    if entry['Summary'] == 'AssemblyAlign':
                        centry['%AssAlnCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Validity
                        centry['%AssAlnAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    if entry['Summary'] == 'Genes':
                        centry['%GeneCov'] = 100.0 * entry['Coverage'] / entry['Length']
                        centry['%GeneAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    if entry['Summary'] == 'Proteins':
                        centry['%ProtCov'] = 100.0 * entry['Coverage'] / entry['Length']
                        centry['%ProtAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    #i# The Reciprocal Gene searches do not seem to very useful as summary data.
                    #if entry['Summary'] == 'Genes.Reciprocal':
                    #    centry['%GeneIntegrity'] = 100.0 * entry['Identity'] / entry['Length']
                    #if entry['Summary'] == 'Proteins.Reciprocal':
                    #    centry['%ProtIntegrity'] = 100.0 * entry['Identity'] / entry['Length']
                centry = compdb.addEntry(centry)
                #self.debug(centry)
                #self.debug('%s' % compdb.data(pbase))
                ## ~ [1b] Load and process CovPlot Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #MBG8150.SP16481.hcq.sgd.srt.1000.covplot.chrom.tdt
                cfile = string.replace(pfile,'Summary','covplot.chrom')
                if not rje.exists(cfile): self.warnLog('Could not locate %s' % cfile); continue
                try: cdb = db.addTable(cfile,['Chrom','Pos'],name='covplot',expect=True)
                except: self.errorLog('Cannot load PAGSAT chromosome coverage table "%s": check format' % cfile); continue
                cdb.dataFormat({'Pos':'int','ContigNum':'int','ChromNum':'int'})
                # ['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom']
                covdata = {}    # Dict of {chrom:{pos:excess}}
                for entry in cdb.entries():
                    if entry['Chrom'] not in covdata: covdata[entry['Chrom']] = {}
                    covdata[entry['Chrom']][entry['Pos']] = entry['ContigNum'] - entry['ChromNum']
                centry['Duplicate'] = 0
                for chrom in covdata:
                    cpos = rje.sortKeys(covdata[chrom])
                    (x,i) = (0,0)
                    while i < len(cpos):
                        cx = covdata[chrom][cpos[i]]    # Contig hits - chrom hits
                        if cx > 0: centry['Duplicate'] +=  cx * (cpos[i] - x)
                        x = cpos[i]
                        i += 1
                ## ~ [1c] Unique coverage and duplication ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                bad_loc2chr = False
                if ftdb:        #!# NOTE: Unique stuff is not working!
                    centry['RepeatFT'] = string.join(repeatft,';')
                    # Note: locus=BK006937; Chrom=chrIII_S288C__BK006937.2
                    if not loc2chr:     # Only need to make once
                        for locus in ftdb.index('locus'):
                            loc2chr[locus] = None
                            for chrom in cdb.index('Chrom'):
                                if string.split(chrom,'__')[1].startswith(locus): loc2chr[locus] = chrom
                            if not loc2chr[locus]: self.warnLog('CovPlot missing %s?' % locus); bad_loc2chr = True
                        #self.debug(loc2chr)
                    # First, modify cdb to include R entries where Contig hits = Chrom hits = 0
                    for ft in ftdb.entries():
                        chrom = loc2chr[ft['locus']]
                        if not chrom or not cdb.indexEntries('Chrom',chrom): continue
                        fmin = ()   # (pos,centry) closest to 5' end of feature (to be at Pos=start-1)
                        fmax = ()   # (pos,centry) closest to 3' end of feature (to be at Pos=end+1)
                        for chrentry in cdb.indexEntries('Chrom',chrom)[0:]:
                            if chrentry['Pos'] < (ft['start']-1):
                                if not fmin or chrentry['Pos'] > fmin[0]: fmin = (chrentry['Pos'],chrentry)
                                continue    # No overlap, so keep
                            if chrentry['Pos'] > (ft['end']+1):
                                if not fmax or chrentry['Pos'] < fmax[0]: fmax = (chrentry['Pos'],chrentry)
                                continue    # No overlap, so keep
                            if not fmin or chrentry['Pos'] < fmin[0]: fmin = (chrentry['Pos'],chrentry)
                            if not fmax or chrentry['Pos'] > fmax[0]: fmax = (chrentry['Pos'],chrentry)
                            cdb.dropEntry(chrentry)
                        chrentry = fmin[1]
                        cdb.addEntry(rje.combineDict({'Chrom':chrom,'Pos':ft['start']-1},chrentry,overwrite=False))
                        cdb.addEntry({'Chrom':chrom,'Pos':ft['start'],'ContigNum':0,'Class':'R','ChromNum':0})
                        cdb.addEntry({'Chrom':chrom,'Pos':ft['end'],'ContigNum':0,'Class':'R','ChromNum':0})
                        chrentry = fmax[1]
                        cdb.addEntry(rje.combineDict({'Chrom':chrom,'Pos':ft['end']+1},chrentry,overwrite=False))
                    # Then calculate UniqDup as before
                    covdata = {}    # Dict of {chrom:{pos:excess}}
                    uniqdata = {}   # Dict of {chrom:{pos:class}}
                    uniqlen = totlen = 0; uniqcov = uniqctg = 0
                    for entry in cdb.entries():
                        if entry['Chrom'] not in covdata: covdata[entry['Chrom']] = {}; uniqdata[entry['Chrom']] = {}
                        try: covdata[entry['Chrom']][entry['Pos']] = entry['ContigNum'] - entry['ChromNum']
                        except: self.debug(entry)
                        uniqdata[entry['Chrom']][entry['Pos']] = entry['Class']
                    centry['UniqDup'] = 0
                    for chrom in covdata:
                        cpos = rje.sortKeys(covdata[chrom])
                        totlen += cpos[-1] - 1
                        uniqlen += cpos[-1] - 1
                        (x,i) = (0,0)
                        while i < len(cpos):
                            if uniqdata[chrom][cpos[i]] == 'R': uniqlen -= cpos[i] - x
                            elif uniqdata[chrom][cpos[i]] != 'N': uniqcov += cpos[i] - x    # Was C/U but I think
                            if uniqdata[chrom][cpos[i]] in ['C','U']: uniqctg += cpos[i] - x    # Was C/U but I think
                            cx = covdata[chrom][cpos[i]]    # Contig hits - chrom hits
                            if cx > 0: centry['UniqDup'] +=  cx * (cpos[i] - x)
                            x = cpos[i]
                            i += 1
                    centry['UniqCov'] = 100.0 * uniqcov / uniqlen
                    centry['UniqCtg'] = 100.0 * uniqctg / uniqlen
                db.deleteTable(cdb)
                ## ~ [1d] Load an process Tree file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                tfile = string.replace(pfile,'Summary.tdt','nwk')
                tree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
                tree.loadTree(tfile,postprocess=False)
                centry['TreeLen'] = tree.treeLen()
                ## ~ [1e] FragX and ChrX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# NOTE: This needs to be improved with respect to checking options etc. #!#
                #lfile = '%s/%s.local.tdt' % (string.replace(pfile,'Summary.tdt','GABLAM'),pbase)
                # Not: MBG479.SP16499.hcq.qv20.sgd.srt.PAGSAT/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.GABLAM/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.local.tdt
                # -> MBG479.SP16499.hcq.qv20.sgd.srt.GABLAM/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.Reference.local.tdt
                lfile = '%s.local.tdt' % self.fileBase('GABLAM','Cut','Reference')
                if not rje.exists(lfile): self.warnLog('Could not locate %s' % lfile); continue
                try: locdb = db.addTable(lfile,['Qry','Hit','AlnNum'],name='local',expect=True)
                except: self.errorLog('Cannot load local hit table "%s": check format' % lfile); continue
                locdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                for lentry in locdb.entries():
                    if lentry['SbjStart'] > lentry['SbjEnd']:
                        (lentry['SbjStart'],lentry['SbjEnd']) = (lentry['SbjEnd'],lentry['SbjStart'])
                # Use uniqlen and totlen calculated above to count number of local BLAST hits needed to exceed length thresholds
                #fragcov = [50,90,95,99]     # List of coverage thresholds to count (local table)
                for fragx in fragcov: centry['Frag%d' % fragx] = centry['UniqFrag%d' % fragx] = 0
                # Make lists of coverage (start, end), merging as required, sum up and compare to totlen
                # For uniqlen, start with a list of features before adding local hits
                covdict = {'Frag':{},'UniqFrag':{}}    # Dictionary of {chromosome:[(start,end)]
                covtot = {'Frag':{},'UniqFrag':{}}  # Dictionary of {Chromosome:total coverage}
                covlen = {'Frag':totlen,'UniqFrag':uniqlen}  # Dictionary of {Chromosome:total coverage}
                for qry in locdb.index('Qry'):
                    covdict['Frag'][qry] = []
                    covdict['UniqFrag'][qry] = []
                    covtot['Frag'][qry] = 0
                ucovdict = covdict['UniqFrag']
                for ft in ftdb.entries():
                    qry = loc2chr[ft['locus']]
                    if not qry: continue
                    ucovdict[qry].append((ft['start'],ft['end']))
                for qry in locdb.index('Qry'):
                    ucovdict[qry].sort()
                    x = 1
                    while x < len(ucovdict[qry]):
                        if ucovdict[qry][x][0] <= (ucovdict[qry][x-1][1] + 1):    # Merge
                            ucovdict[qry][x-1] = (ucovdict[qry][x-1][0],max(ucovdict[qry][x-1][1],ucovdict[qry][x][1]))
                            ucovdict[qry].pop(x)
                        else: x += 1
                    covtot['UniqFrag'][qry] = 0
                    for (i,j) in ucovdict[qry]: covtot['UniqFrag'][qry] += (j - i + 1)
                    #self.debug(ucovdict[qry])
                    #self.debug(covtot['UniqFrag'][qry])
                uniqlen = sum(covtot['UniqFrag'].values())
                covlen['UniqFrag'] = totlen - uniqlen
                # Add local hits in size order.
                hitx = 0    # Hit counter
                for lentry in locdb.sortedEntries('Identity',reverse=True):
                    hitx += 1
                    qry = lentry['Qry']
                    for c in ['Frag','UniqFrag']:
                        qfrag = covdict[c][qry]
                        qfrag.append((lentry['QryStart'],lentry['QryEnd']))
                        qfrag.sort()
                        x = 1
                        while x < len(qfrag):
                            if qfrag[x][0] <= (qfrag[x-1][1] + 1):    # Merge
                                qfrag[x-1] = (qfrag[x-1][0],max(qfrag[x-1][1],qfrag[x][1]))
                                qfrag.pop(x)
                            else: x += 1
                        covtot[c][qry] = 0
                        for (i,j) in qfrag: covtot[c][qry] += (j - i + 1)
                        # Assess coverage:
                        for fragx in fragcov:
                            if centry['%s%d' % (c,fragx)]: continue
                            if c == 'UniqFrag':
                                if 100.0 * (sum(covtot[c].values()) - uniqlen) / covlen[c] >= fragx: centry['%s%d' % (c,fragx)] = hitx
                            elif 100.0 * sum(covtot[c].values()) / covlen[c] >= fragx: centry['%s%d' % (c,fragx)] = hitx
                    if centry['Frag%d' % fragcov[-1]] and centry['UniqFrag%d' % fragcov[-1]]: break
                db.deleteTable(locdb)
                if bad_loc2chr: loc2chr = {}

                #chromcov = [50,95,98,99]    # No. of chromosomes covered by a single contig (GABLAM table)
                #gfile = string.replace(lfile,'local','gablam')
                gfile = '%s.gablam.tdt' % self.fileBase('GABLAM','Cut','Reference')
                if not rje.exists(gfile): self.warnLog('Could not locate %s' % gfile); continue
                try: gdb = db.addTable(gfile,['Qry','Hit'],name='gablam',expect=True)
                except: self.errorLog('Cannot load GABLAM table "%s": check format' % gfile); continue
                gdb.dropEntriesDirect('Rank',['1'],inverse=True)
                gxfield = 'Qry_AlnID'   # Could also use Qry_AlnLen
                gdb.dataFormat({gxfield:'float'})
                for chromx in chromcov:
                    centry['Chrom%d' % chromx] = 0
                    for gentry in gdb.entries():
                        if gentry[gxfield] >= chromx: centry['Chrom%d' % chromx] += 1

                db.deleteTable(gdb)

                self.debug(centry)
            compdb.saveToFile()



            pagfiles = []
            pagbase = []    # List of basefiles for PAGSAT results

            # Load in summary table, add assembly name and then combine with others
            # Reshape wide and then reshape long again!

            #Summary	HitNum	Length	Coverage	Identity	Positives	Missing	Errors	Perfect	N
            #Assembly	1173	13235834	13233590	13232765	13232765	2244	825	28	120
            #Reference	1190	12157104	12124470	12123837	12123837	32634	633	0	17
            #Self	268	12157104	12157104	12157104	12157104	0	0	17	17


            datatypes = ['Genes','Genes.Reciprocal','Proteins','Proteins.Reciprocal','Reference','Self']


        except: self.errorLog('%s.compare error' % self.prog())
#########################################################################################################################
### End of SECTION II: PAGSAT Class                                                                                     #
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
    try: PAGSAT(mainlog,cmd_list).run()

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
