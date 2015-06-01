#!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      GFESSA
Description:  Genome-Free EST SuperSAGE Analysis
Version:      1.4
Last Edit:    20/08/13
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program is for the automated processing, mapping and identification-by-homology for SuperSAGE tag data for
    organisms without genome sequences, relying predominantly on EST libraries etc. Although designed for genome-free
    analysis, there is no reason why transcriptome data from genome projects cannot be used in the pipeline. 

    GFESSA aims to take care of the following main issues:
    1. Removal of unreliable tag identification/quantification based on limited count numbers.
    2. Converting raw count values into enrichment in one condition versus another.
    3. Calculating mean quantification for genes based on all the tags mapping to the same sequence.
    4. The redundancy of EST libraries, by mapping tags to multiple sequences where necessary and clustering sequences
    on shared tags.

    The final output is a list of the sequences identified by the SAGE experiment along with enrichment data and
    clustering based on shared tags.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    tagfile=FILE    : File containing SuperSAGE tags and counts [None]
    tagfield=X      : Field in tagfile containing tag sequence. (All others should be counts) ['Tag sequence']
    expconvert=FILE : File containing 'Header', 'Experiment' conversion data [None]
    experiments=LIST: List of (converted) experiment names to use []
    seqin=FILE      : File containing EST/cDNA data to search for tags within [None]
    tagindex=FILE   : File containing possible tags and sequence names from seqin file [*.tag.index]
    tagmap=FILE     : Tag to sequence mapping file to over-ride auto-generated file based on Seqin and Mismatch [None]
    tagstart=X      : Sequence starting tags ['CATG']
    taglen=X        : Length of sequence tags [26]
    ### ~ PROCESS OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    mintag=X        : Minimum total number of counts for a tag to be included (summed replicates) [0]
    minabstag=X     : Minimum individual number of counts for a tag to be included (ANY one replicate) [5]
    minexptag=X     : Minimum number of experiments for a tag to be included (no. replicates) [3]
    allreptag=X     : Filter out any Tags that are not returned by ALL replicates of X experiments [0]
    minenrtag=X     : Minimum number of counts for a tag to be retained for enrichment etc. (summed replicates) [15]
    enrcut=X        : Minimum mean fold change between experiments [2.5]
    pwenr=X         : Minimum fold change between pairwise experiment comparisons [1.0]
    expand=T/F      : Whether to expand from enriched TAGs to unenriched TAGs through shared sequence hits [True]
    mismatch=X      : No. mismatches to allow. -1 = Exact matching w/o BLAST [-1]
    bestmatch=T/F   : Whether to stop looking for more mismatches once hits of a given stringency found [True]
    normalise=X     : Method for normalising tag counts within replicate (None/ppm) [ppm]
    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : "Base" name for all results files, e.g. X.gfessa.tdt [TAG file basename]
    longtdt=T/F     : Whether to output "Long" format file needed for R analysis [True]

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import budapest 
import rje_seq, rje_seqlist, rje_sequence
import rje_blast_V2 as rje_blast
import rje, rje_db, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation using exact matches only.
    # 0.1 - BLAST-based inexact search method.
    # 0.2 - Removed sequence annotation and clustering. Added extra enrichment clustering.
    # 1.0 - Updated to fix basefile issue and improve documentation, including manual. Add mean cluster enrichment.
    # 1.1 - Added  minabstag and minexptag to give more control over low abundance tag filtering
    # 1.2 - Added longtdt to output "Long" format file needed for R analysis.
    # 1.3 - Tidied module imports.
    # 1.4 - Switched to rje_blast_V2. More work needed for BLAST+.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Store sequence IDs with mismatches - look at no mismatches first.
    # [Y] : Pre-filter based on enrichment.
    # [?] : Possible optional FIESTA Assembly of mapped ESTs, followed by remapping, during mapped EST processing.
    # [ ] : Add reading of existing results from previous normalisation/enrichment/enrichment filtering.
    # [ ] : Add rje_menu interactive mode.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('GFESSA', '1.4', 'July 2011', '2011')
    description = 'Genome-Free EST SuperSAGE Analysis'
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GFESSA Class                                                                                            #
#########################################################################################################################
class GFESSA(budapest.Budapest):     
    '''
    GFESSA Class. Author: Rich Edwards (2011).

    Info:str
    - ExpConvert = File containing 'Header', 'Experiment' conversion data [None]
    - Name = Name of input TAG file
    - Normalise = Method for normalising tag counts within replicate (None/ppm) [ppm]
    - TagField = Field in tagfile containing tag sequence. (All others should be counts) ['Tag sequence']
    - TagIndex = File containing possible tags and sequence names from seqin file [*.tag.index]
    - TagMap = Tag to sequence mapping file to over-ride auto-generated file based on Seqin and Mismatch [None]
    - TagStart = Sequence starting tags ['GATC']
    
    Opt:boolean
    - BestMatch = Whether to stop looking for more mismatches onces hits of a given stringency found [True]
    - Expand = Whether to expand from enriched TAGs to unenriched TAGs through shared sequence hits [True]
    - LongTDT = Whether to output "Long" format file needed for R analysis [True]

    Stat:numeric
    - EnrCut = Minimum mean fold change between experiments [2.5]
    - MinAbsTag = Minimum individual number of counts for a tag to be included (ANY one replicate) [5]
    - MinExpTag = Minimum number of experiments for a tag to be included (no. replicates) [3]
    - AllRepTag = Filter out any Tags that are not returned by ALL replicates of X experiments [0]
    - MinEnrTag = Minimum number of counts for a tag to be retained for enrichment etc. (summed replicates) [15]
    - MinTag = Minimum number of counts for a tag to be included (summed replicates) [0]
    - Mismatch = No. mismatches to allow. -1 = Exact matching w/o BLAST [-1]
    - PWEnr = Minimum fold change between pairwise experiment comparisons [1.0]
    - TagLen = Length of sequence tags [26]
    
    List:list
    - Experiments = List of (converted) experiment names to use []

    Dict:dictionary
    # Replaced with database objects? #
    - Details = For each Hit, a list of text details to be output to *.details.txt
    - Peptides = Dictionary of {Translation shortName:[Peptide list]}
    - PepSeq = Dictionary of {peptide:[Sequence objects]}
    - PepTypes = Dictionary of {'Common':[peplist],'Cluster':[peplist],'Unique':[peplist]}
    - RF = Dictionary of {EST ID:[Acceptable RF translations]}
    - Support = Dictionary of {no. peptides:[List of IDs with this many peptides]}

    Obj:RJE_Objects
    - DB = Database object storing relevant data
    - SeqList = SeqList object containing input ESTs
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['ExpConvert','Name','Normalise','SeqIn','TagField','TagIndex','TagStart','TagMap']
        self.optlist = ['Expand','BestMatch','LongTDT']
        self.statlist = ['EnrCut','MinEnrTag','MinTag','Mismatch','PWEnr','TagLen','MinAbsTag','MinExpTag','AllRepTag']
        self.listlist = ['Experiments']
        self.dictlist = []
        self.objlist = ['DB','SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'TagStart':'CATG','Normalise':'ppm','TagField':'Tag sequence'})
        self.setStat({'EnrCut':2.5,'MinEnrTag':15,'MinTag':0,'Mismatch':-1,'PWEnr':1.0,'TagLen':26,'MinAbsTag':5,'MinExpTag':3,'AllRepTag':0})
        self.setOpt({'Expand':True,'BestMatch':True,'LongTDT':True})
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
                ### Class Options ### 
                self._cmdRead(cmd,type='info',att='Name',arg='tagfile')     # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'info',['Normalise','TagField','TagStart'])
                self._cmdReadList(cmd,'file',['ExpConvert','SeqIn','TagIndex','TagMap'])
                self._cmdReadList(cmd,'int',['MinTag','Mismatch','MinEnrTag','TagLen','MinAbsTag','MinExpTag','AllRepTag'])
                self._cmdReadList(cmd,'opt',['BestMatch','Expand','LongTDT'])
                self._cmdReadList(cmd,'num',['EnrCut','PWEnr'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False    # Loads basic input data
            ### ~ [2] ~ Map Tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mdb,normalise,enrich) = self.mapTags()
            ### ~ [3] ~ Filter and Normalise Tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if normalise:
                self.setBasefile('%s.%da%de%dr%d' % (self.basefile(),self.stat['MinTag'],self.stat['MinAbsTag'],self.stat['MinExpTag'],self.stat['AllRepTag']))
                self.normalise(mdb)
            if self.getBool('LongTDT'): self.longTDT(mdb)
            ### ~ [4] ~ Tag Enrichment Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setBasefile('%s.enr%d' % (self.basefile(),self.stat['MinEnrTag']))
            if enrich: self.enrichment(mdb)
            ### ~ [5] ~ Process mapped ESTs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.enrFilter(mdb)                
            self.processHits(mdb)
            return True
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def seqlist(self):  ### Loads and/or returns SeqList object
        if not self.obj['SeqList']: self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F'])
        return self.obj['SeqList']
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                                #V1.0
        '''Main class setup method.'''
        try:### ~ [0] General Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Check files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.checkForFile(self.info['Name']): self.printLog('#FILE','Tag File "%s": found' % self.info['Name'])
            else: self.printLog('#FILE','Tag File "%s": Missing!' % self.info['Name']); raise IOError
            if rje.checkForFile(self.info['SeqIn']): self.printLog('#FILE','Sequence File "%s": found' % self.info['SeqIn'])
            else: self.printLog('#FILE','Sequence File "%s": Missing!' % self.info['SeqIn']); raise IOError
            ## ~ [0b] Setup Basefile etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setStr({'InfileBase': rje.baseFile(self.info['Name'],strip_path=True)})
            self.setStr({'GFESSABase': '%s.%s' % (rje.baseFile(self.info['SeqIn'],strip_path=True),self.getStr('InfileBase'))})
            self.info['CoreBaseFile'] = self.basefile()
            if self.info['Basefile'].lower() in ['','none']:
                self.setBasefile(self.getStr('InfileBase'))
                self.printLog('#BASE','Output file names: %s.*' % self.basefile())
            else: self.printLog('#BASE','Output file names: %s.* and %s.*' % (self.getStr('InfileBase'),self.info['Basefile']))
            ### ~ [1] Load Input Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tagdb = self.setupTagFile()
            ### ~ [2] Generate Tag Index for SeqIn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['Mismatch'] < 0 and self.info['TagStart']: idb = self.setupTagIndex()
            ### ~ [3] Setup Basefile for rest of analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = '%s.%s' % (rje.baseFile(self.info['SeqIn'],strip_path=True),rje.baseFile(self.info['Name'],strip_path=True))
            self.printLog('#BASE','Additional output file names: %s.*' % self.info['Basefile'])
            if self.stat['Mismatch'] >= 0: basefile += '.%d' % self.stat['Mismatch']
            self.setBasefile(basefile)
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def setupTagFile(self): ### Setup main tag file                                                                 #V1.0
        '''Setup main tag file. This contains the input counts for each tag from the SuperSAGE experiment.'''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Check for existing Tag table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tagdb = None
            for base in (self.basefile(),self.getStr('InfileBase')):
                tfile = '%s.Tag.tdt' % base
                if rje.exists(tfile) and (not self.opt['Force'] or self.yesNo('Use existing %s?' % tfile,default='N')):
                    try: tagdb = self.db().addTable(tfile,['Tag'],name='Tag')
                    except: tagdb = None
                if tagdb or self.basefile() == self.getStr('InfileBase'): break
            ## ~ [1b] Load Tag table if missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not tagdb:
                try:
                    tagdb = self.db().addTable(self.info['Name'],[self.getStr('TagField')],name='Tag')
                    tagdb.renameField(self.getStr('TagField'),'Tag')
                except: tagdb = self.db().addTable(self.info['Name'],['Tag'],name='Tag')
                if rje.exists(self.info['ExpConvert']):
                    edb = self.db().addTable(self.info['ExpConvert'],['Header'],['Experiment'],headers=['Header','Experiment'],name='ExpConvert')
                    for header in edb.dataKeys():
                        if header in tagdb.fields(): tagdb.renameField(header,edb.data()[header]['Experiment'])
            ## ~ [1c] Set experiment list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            all_experiments = []; conditions = {}; condition_warn = 0
            for field in tagdb.fields():
                if field == 'Tag': continue
                if field == 'Count': break
                all_experiments.append(field)
                try: conditions[field] = rje.matchExp('^(\S+)\d+',field)[0]
                except: condition_warn += 1; conditions[field] = field
            if condition_warn: self.errorLog('Could not pull out Condition and Replicate for %d fields' % condition_warn,printerror=False)
            self.dict['Conditions'] = conditions
            ## ~ [1d] ~ Format Tag File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            reformat = {'Count':'int','MaxCount':'int','ECount':'int','AllRep':'int'}
            for field in all_experiments: reformat[field] = 'int'
            tagdb.dataFormat(reformat)
            ### ~ [2] ~ Generate Sequence Counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'Count' not in tagdb.fields():
                tagdb.addField('Count',evalue=0); tagdb.addField('ECount',evalue=0); tagdb.addField('MaxCount',evalue=0)
                tagdb.addField('AllRep',evalue=0)
                for entry in tagdb.entries():
                    entry['AllRep'] = rje.sortUnique(conditions.keys())
                    for field in all_experiments:
                        entry['Count'] += entry[field]
                        entry['MaxCount'] = max(entry['MaxCount'],entry[field])
                        if entry[field] > 0: entry['ECount'] += 1
                        elif conditions[field] in entry['AllRep']: entry['AllRep'].remove(conditions[field])
                    entry['AllRep'] = len(entry['AllRep'])
                self.printLog('#TAG','Tag counts for %d input experiments calculated' % len(all_experiments))
                tagdb.saveToFile()  # Full counts for all experiments saved to file
            ## ~ [2a] ~ Defined possible reduced experiment set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Experiments']:
                process = []
                for experiment in self.list['Experiments']:
                    if experiment in all_experiments: process.append(experiment)
                    else: self.printLog('#ERR','Experiment "%s" not found in Input' % experiment)
                self.list['Experiments'] = process
                for experiment in process: conditions.pop(experiment)
            else: self.list['Experiments'] = all_experiments
            self.printLog('#EXP','%d experiments to process: %s.' % (len(self.list['Experiments']),string.join(self.list['Experiments'],'; ')))
            ## ~ [2b] ~ Recalculate for reduced experiment set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if len(self.list['Experiments']) != len(all_experiments):
                for entry in tagdb.entries():
                    entry['Count'] = 0; entry['ECount'] = 0; entry['MaxCount'] = 0
                    entry['AllRep'] = rje.sortUnique(conditions.keys())
                    for field in self.list['Experiments']:
                        entry['Count'] += entry[field]
                        entry['MaxCount'] = max(entry['MaxCount'],entry[field])
                        if entry[field] > 0: entry['ECount'] += 1
                        elif conditions[field] in entry['AllRep']: entry['AllRep'].remove(conditions[field])
                    entry['AllRep'] = len(entry['AllRep'])
                self.printLog('#TAG','Tag counts recalculated for %d selected experiments' % len(self.list['Experiments']))
            ### ~ [3] ~ Finish and return Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return tagdb
        except: self.errorLog('%s.setupTagFile error' % self)
#########################################################################################################################
    def setupTagIndex(self):    ### Sets up index for sequence file of all possible SuperSAGE tags.                 #V1.0
        '''Sets up index for sequence file of all possible SuperSAGE tags. Only used for BLAST-free method.'''
        try:### ~ [0] General Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.info['TagStart']: raise ValueError
            self.setInt({'TagLen':max(len(self.info['TagStart']),self.getInt('TagLen'))})
            ## ~ [0a] Check and load pre-existing tag file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tagindex = '%s.%s.%d.tag.index' % (rje.baseFile(self.info['SeqIn']),self.info['TagStart'],self.getInt('TagLen'))
            if not self.needToRemake(tagindex,self.info['SeqIn']) or (os.path.exists(tagindex) and self.yesNo('Use existing %s?' % tagindex,default='N')):
                return self.db().addTable(tagindex,['#'],['Tag','Seq'],name='TagSeq')
            ### ~ [1] Generate Tag Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','usecase=False'])
            tagdict = {}; nx = 0; sx = 0
            (fwd,rev) = (self.info['TagStart'],rje_sequence.reverseComplement(self.info['TagStart']))
            while seqlist.nextSeq():
                self.progLog('\r#TAG','Indexing possible tags in %s: %s -> %s' % (seqlist.getStr('Name'),rje.iStr(sx),rje.iStr(nx))); sx += 1
                (name,fseq) = seqlist.currSeq()
                tags = []; name = string.split(name)[0]
                if fwd in fseq:
                    frag = fseq[0:]
                    while fwd in frag:
                        frag = string.join(string.split(frag,fwd)[1:],fwd)
                        tag = self.info['TagStart'] + frag
                        if len(tag) > self.getInt('TagLen'): tags.append(tag[:self.getInt('TagLen')])
                if rev in fseq:
                    frag = rje_sequence.reverseComplement(fseq)
                    while fwd in frag:
                        frag = string.join(string.split(frag,fwd)[1:],fwd)
                        tag = self.info['TagStart'] + frag
                        if len(tag) > self.getInt('TagLen'): tags.append(tag[:self.getInt('TagLen')])
                for tag in tags:
                    if 'N' in tag: continue
                    if tag not in tagdict: tagdict[tag] = []
                    if name not in tagdict[tag]: tagdict[tag].append(name); nx += 1
            self.printLog('\r#TAG','Found %s occurrences of %s different possible tags in %s %s seqs.' % (rje.iStr(nx),rje.iLen(tagdict.keys()),rje.iStr(sx),seqlist.getStr('Name')))
            seqlist.tidy()
            ### ~ [2] Generate Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            idb = self.db().addEmptyTable('TagSeq',['#','Tag','Seq'],['#']); tx = 0
            for tag in rje.sortKeys(tagdict):
                tagdict[tag].sort()
                for name in tagdict[tag]:
                    tx += 1; i = rje.preZero(tx,nx)
                    idb.dict['Data'][i] = {'#':i,'Tag':tag,'Seq':name}
            if not self.opt['Test']: idb.saveToFile(tagindex,'\t')
            return idb     # Setup successful
        except:
            self.errorLog('Problem during %s setupTagIndex.' % self)
            if self.yesNo('Use BLAST-based method instead?',default='N'): return None
            raise   # Setup failed
#########################################################################################################################
    ### <3> ### Tag Mapping Methods                                                                                     #
#########################################################################################################################
    def mapTags(self):  ### Map SuperSAGE tags onto EST/cDNA sequences                                              #V1.0
        '''Map SuperSAGE tags onto EST/cDNA sequences.'''
        try:### ~ [1] Check for existing file to load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            normbase = '%s.%da%de%dr%d' % (self.basefile(),self.stat['MinTag'],self.stat['MinAbsTag'],self.stat['MinExpTag'],self.stat['AllRepTag'])
            if self.info['Normalise'].lower() == 'ppm': normbase += '.ppm'
            reformat = {'Count':'int','ECount':'int','SeqN':'int','Tag':'str','MaxCount':'int','AllRep':'int'}
            longfile = '%s.Long.tdt' % normbase
            ## ~ [1a] Enriched, normalised file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            enrfile = '%s.enr%d.Tag.tdt' % (normbase,self.stat['MinEnrTag'])
            if rje.exists(enrfile) and (rje.exists(longfile) or not self.getBool('LongTDT')) and not self.opt['Force']:
                tdb = self.db().addTable(enrfile,['Tag'],name='Tag')
                for field in tdb.fields():
                    if field not in reformat: reformat[field] = 'num'
                tdb.dataFormat(reformat)
                if self.info['Normalise'].lower() == 'ppm': self.printLog('#ENR','Loaded enriched and normalised data from %s' % enrfile)
                else: self.printLog('#ENR','Loaded enriched data from %s' % enrfile)
                self.setBasefile(normbase)
                self.setComparisons()
                return (tdb,False,False)
            ## ~ [1b] Normalised file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            normfile = '%s.Tag.tdt' % normbase
            if rje.exists(normfile) and not self.opt['Force']:
                tdb = self.db().addTable(normfile,['Tag'],name='Tag')
                for field in tdb.fields():
                    if field not in reformat: reformat[field] = 'num'
                tdb.dataFormat(reformat)
                if self.info['Normalise'].lower() == 'ppm': self.printLog('#NORM','Loaded normalised data from %s' % enrfile)
                else: self.printLog('#TAG','Loaded Tag data from %s' % enrfile)
                self.setBasefile(normbase)
                return (tdb,False,True)
            ## ~ [1c] Basic TagMap file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            if self.info['TagMap'].lower() not in ['','none']: tsfile = self.info['TagMap']
            else: tsfile = '%s.TagMap.tdt' % self.basefile()
            if rje.exists(tsfile) and not self.opt['Force']:
                tdb = self.db().addTable(tsfile,['Tag'],name='TagMap')
                for field in self.list['Experiments']: reformat[field] = 'float'    # Input might be normalised already!
                for mm in range(self.stat['Mismatch']): reformat['SeqN%d' % mm] = 'int'
                tdb.dataFormat(reformat); tdb.info['Basefile'] = self.basefile()
                return (tdb,True,True)
            ### ~ [2] Perform Tag Mapping using BLAST if index not used ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.db('TagSeq'): return (self.blastMapTags(),True,True)
            ### ~ [3] Join sequence-tag mapping with experimental results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tdb = self.db().joinTables('TagMap',join=[('Tag','Tag'),('TagSeq','Tag')])
            tdb.compress(['Tag'],rules={'Seq':'list'},default='max',best=[])
            tdb.deleteField('#'); tdb.deleteField('AutoID'); tdb.fillBlanks()
            ### ~[4] Count Hit Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tdb.addField('SeqN'); seqn = 0
            (tx,ttot) = (0.0,tdb.entryNum())
            for entry in tdb.entries():
                self.progLog('\r#SEQN','Tag sequence counts: %.2f%%' % (tx/ttot)); tx += 100.0
                if entry['Seq']: entry['SeqN'] = len(string.split(entry['Seq'],';')); seqn += entry['SeqN']
                else: entry['SeqN'] = 0
            self.printLog('\r#SEQN','Tag sequence counts: %s tags; %s Seq-Tag pairs.' % (rje.iStr(tdb.entryNum()),rje.iStr(seqn)))
            tdb.saveToFile()
            return (tdb,True,True)
        except: self.errorLog('%s.mapTags error' % self); raise
#########################################################################################################################
    def blastMapTags(self): ### Map SuperSAGE tags onto EST/cDNA sequences using BLAST                              #V1.0
        '''Map SuperSAGE tags onto EST/cDNA sequences using BLAST.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.seqlist()
            tdb = self.db().copyTable('Tag','TagMap')
            #tdb.info['Basefile'] = self.db().info['Basefile']
            tdb.addField('SeqN',evalue=0); tdb.addField('Seq',evalue=[])
            for x in range(self.stat['Mismatch']+1): tdb.addField('SeqN%d' % x,evalue=0); tdb.addField('Seq%d' % x,evalue='')
            ## ~ [0a] ~ Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.stat['Forks'] > 1 and not self.opt['NoForks']: self.cmd_list = ['blasta=%d' % self.stat['Forks']] + self.cmd_list
            blast = rje_blast.blastObj(self.log,['blastv=1000','blastb=1000','blaste=0.0001']+self.cmd_list+['blastf=T','blastg=F','blastf=F','blastcf=F','blastp=blastn'],type='Dev')
            ## ~ [0a] ~ Setup fasta file of Tag sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqfile = rje.baseFile(self.info['Name']) + '.tag.fas'
            if self.opt['Force'] or not rje.exists(seqfile):
                output = []
                for tag in self.db('Tag').datakeys(): output.append('>%s\n%s' % (tag,tag))
                open(seqfile,'w').write(string.join(output,'\n'))
            ## ~ [0b] ~ Format Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast.setInfo({'InFile':seqfile,'DBase':self.info['SeqIn'],
                           'Name':'%s.%s.blast' % (rje.baseFile(self.info['Name'],True),rje.baseFile(self.info['SeqIn'],True))})
            blast.formatDB(fasfile=self.info['SeqIn'],force=self.opt['Force'],protein=False)
            if not rje_blast.checkForDB(dbfile=self.info['SeqIn'],checkage=False,log=self.log,protein=False,oldblast=blast.getBool('OldBLAST')):
                self.errorLog('FormatDB failed for unknown reasons. Check blastpath=X and permissions.',printerror=False)
                raise IOError
            ### ~ [1] ~ Perform and read BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#BLAST','Performing BLAST of Tags versus SeqIn. May take some time!')
            blast.blast(cleandb=False,use_existing=not self.opt['Force'])
            RESFILE = open(blast.getStr('Name'),'r'); RESFILE.seek(0,2); fend = RESFILE.tell(); RESFILE.close()
            ## ~ [1a] ~ Add to TagMap Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hx = 0; mx = 0; fpos = 0; maxmm = 0
            while fpos >= 0:
                self.progLog('\r#HIT','Extracting %s hits from BLAST: %.2f%% (%s rejected > %d mismatch)' % (rje.iStr(hx),(100.0*fpos/fend),rje.iStr(mx),self.stat['Mismatch']))
                (search,fpos) = blast.readNextBLASTSearch(fpos,gablam=True,local=False)
                if not search: break
                if blast.oldBLAST():
                    tag = search.info['Name']; entry = tdb.data()[tag]; entry['Seq'] = []
                    for x in range(self.stat['Mismatch']+1): entry['Seq%d' % x] = []
                    for hit in search.hit:
                        try:
                            seq = hit.info['Name']; hx += 1
                            mm = len(tag) - hit.dict['GABLAM']['Query']['GABLAMO ID']
                            maxmm = max(maxmm,mm)
                            if mm > self.stat['Mismatch']: mx += 1; continue
                            entry['Seq'].append(seq)
                            entry['Seq%d' % mm].append(seq)
                            entry['SeqN%d' % mm] += 1
                        except: self.errorLog('%s Hit Error' % hit.info['Name'])
                else:
                    tag = search['Query']; entry = tdb.data()[tag]; entry['Seq'] = []
                    for x in range(self.stat['Mismatch']+1): entry['Seq%d' % x] = []
                    for seq in blast.queryHits(tag):
                        try:
                            hx += 1
                            gdict = blast.gablamData(seq,tag)
                            mm = len(tag) - gdict['Query']['GABLAMO ID']
                            maxmm = max(maxmm,mm)
                            if mm > self.stat['Mismatch']: mx += 1; continue
                            entry['Seq'].append(seq)
                            entry['Seq%d' % mm].append(seq)
                            entry['SeqN%d' % mm] += 1
                        except: self.errorLog('%s Hit Error' % hit.info['Name'])
            self.printLog('\r#HIT','Extracted %s hits from BLAST vs %s (%s <= %d mismatch)' % (rje.iStr(hx),self.info['SeqIn'],rje.iStr(hx-mx),self.stat['Mismatch']))
            self.printLog('#MM','These BLAST settings are giving hits of up to %d mismatches.' % (maxmm))
            tx = 0.0; ttot = tdb.entryNum(); seqn = 0
            for entry in tdb.entries():
                self.progLog('\r#SEQN','Tag sequence counts: %.2f%%' % (tx/ttot)); tx += 100.0
                entry['SeqN'] = len(entry['Seq']); entry['Seq'] = string.join(entry['Seq'],';'); seqn += entry['SeqN']
                for x in range(self.stat['Mismatch']+1): entry['Seq%d' % x] = string.join(entry['Seq%d' % x],';'); 
            self.printLog('\r#SEQN','Tag sequence counts: %s tags; %s Seq-Tag pairs.' % (rje.iStr(tdb.entryNum()),rje.iStr(seqn)))
            tdb.saveToFile()
            return tdb
        except: self.errorLog('%s.blastMapTags error' % self); return None
#########################################################################################################################
    ### <4> ### Tag Normalisation Methods                                                                               #
#########################################################################################################################
    def normalise(self,mdb): ### Normalise DB Table Evidence counts                                                 #V1.0
        '''Normalise DB Table Evidence counts.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.info['Name'] = 'Tag'
            mdb.dropEntries(['Count<%d' % self.stat['MinTag']])
            mdb.dropEntries(['MaxCount<%d' % self.stat['MinAbsTag']])
            mdb.dropEntries(['ECount<%d' % self.stat['MinExpTag']])
            mdb.dropEntries(['AllRep<%d' % self.stat['AllRepTag']])
            if self.info['Normalise'].lower() not in ['ppm','','none']:
                self.errorLog('"%s" normalisation not recognised' % self.info['Normalise'],printerror=False)
            ### ~ [1] ~ PPM Normalisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Normalise'].lower() == 'ppm':
                self.setBasefile('%s.ppm' % (self.basefile()))
                etotals = {}
                for e in self.list['Experiments']: etotals[e] = 0
                ex = 0.0; etot = mdb.entryNum()
                for entry in mdb.entries():
                    self.progLog('\r#PPM','PPM Normalisation: %.2f%%' % (ex/etot)); ex += 50.0
                    for e in self.list['Experiments']: etotals[e] += entry[e]
                for entry in mdb.entries():
                    self.progLog('\r#PPM','PPM Normalisation: %.2f%%' % (ex/etot)); ex += 50.0
                    for e in self.list['Experiments']: entry[e] = 1e6 * entry[e] / float(etotals[e])
                for e in self.list['Experiments']: self.printLog('\r#PPM','PPM Normalisation: %s = %s' % (e,rje.iStr(etotals[e])))
            ### ~ [X] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.saveToFile()            
        except: self.errorLog('%s.normalise error' % self)
#########################################################################################################################
    def longTDT(self,mdb):  ### Outputs "Long" format for R analysis
        '''Outputs "Long" format for R analysis.'''
        try:### ~ [0] ~ Setup new file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            longfile = '%s.Long.tdt' % self.basefile()
            if os.path.exists(longfile) and not self.force(): return
            LONG = open(longfile,'w')
            if 'WTL1' in mdb.fields(): LONG.write('%s\n' % string.join(['Tag','Experiment','Expression','Strain','Light'],'\t'))
            else: LONG.write('%s\n' % string.join(['Tag','Experiment','Expression'],'\t'))
            ### ~ [1] ~ Process "Wide" data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = mdb.entryNum(); lx = 0
            for mkey in mdb.dataKeys():
                self.progLog('\r#LONG','Converting to "Long" format: %.2f%%' % (ex/etot)); ex += 100.0
                entry = mdb.data(mkey)
                for field in mdb.fields():
                    if field in ['Tag','Count','ECount','MaxCount','AllRep','Seq','SeqN']: continue
                    elif field[:3] == 'Seq': continue
                    outlist = [entry['Tag'],field,'%s' % entry[field]]
                    if 'WTL1' in mdb.fields():
                        outlist.append(field[:2])
                        outlist.append({'L':'Low','M':'Med','H':'High'}[field[2]])
                    LONG.write('%s\n' % string.join(outlist,'\t')); lx += 1
            self.printLog('\r#LONG','Converted %s "Short" rows to %s "Long" rows.' % (rje.iStr(int(ex)),rje.iStr(lx)))
            LONG.close()
            self.printLog('\r#LONG','Saved %s "Long" data to %s.' % (rje.iStr(lx),longfile))
        except: self.errorLog('%s.longTDT error' % self); self.deBug('!')
#########################################################################################################################
    ### <5> ### Tag Enrichment Methods                                                                                  #
#########################################################################################################################
    def setComparisons(self):   ### Generates self.list['Comparisons'] from self.list['Experiments']
        '''Generates self.list['Comparisons'] from self.list['Experiments'].'''
        try:### ~ [0] ~ Setup list of different conditions from experiments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            conditions = {}
            for e in self.list['Experiments']:
                try:
                    c = rje.matchExp('^(\S+)\d+',e)[0]
                    if c not in conditions: conditions[c] = []
                    conditions[c].append(e)
                except: self.errorLog('Need experiments in form ExpX to calculate enrichment',printerror=False); return False
            self.dict['Conditions'] = conditions
            clist = rje.sortKeys(conditions); comp = []
            ### ~ [1] ~ Generate list of pairs to compare for enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for c1 in clist:
                for c2 in clist:
                    if clist.index(c2) <= clist.index(c1): continue
                    comp.append('%s/%s' % (c1,c2))
            self.list['Comparisons'] = comp
            self.printLog('#COMP','%d pairwise enrichment comparison(s) to make' % len(comp))
        except: self.errorLog('%s.setComparisons error' % self)
#########################################################################################################################
    def enrichment(self,mdb):   ### Calculates enrichment of individual tags
        '''Calculates enrichment of individual tags.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.dropEntries(['Count<%d' % self.stat['MinEnrTag']])
            ## ~ [0a] ~ Regenerate reordered field list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mfields = mdb.fields()
            mdb.list['Fields'] = []
            for field in mfields:
                if field[:3] != 'Seq': mdb.list['Fields'].append(field)
            for field in mfields:
                if field[:4] == 'SeqN': mdb.list['Fields'].append(field)
            for field in mfields:
                if field[:4] != 'SeqN' and field[:3] == 'Seq': mdb.list['Fields'].append(field)
            self.printLog('#DB','TagDB fields: %s' % string.join(mdb.fields()))
            ## ~ [0b] ~ Generate a list of different conditions from experiments ~~~~~~~~~~~~~~~~~~ ##
            self.setComparisons()
            ## ~ [0c] ~ Generate list of pairs to compare for enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for comp in self.list['Comparisons']:
                mdb.addField(comp)
                mdb.addField('minabs%s' % comp)
            comp = self.list['Comparisons']
            conditions = self.dict['Conditions']
            ## ~ [0d] ~ Setup minimum tag count value for absent data (cannot divide by zero) ~~~~~ ##
            minval = 1000.0
            ex = 0.0; etot = mdb.entryNum()
            for entry in mdb.entries():
                self.progLog('\r#ENR','Calculating min Tag level: %.2f%%' % (ex/etot)); ex += 100.0
                for e in self.list['Experiments']:
                    if entry[e]: minval = min(entry[e],minval)
            minval /= max(2.0,self.getNum('EnrCut'))
            ### ~ [1] ~ Calculate Enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = mdb.entryNum()
            for entry in mdb.entries():
                self.progLog('\r#ENR','Calculating Tag enrichment: %.2f%% ' % (ex/etot)); ex += 100.0
                for c in comp:
                    [c1,c2] = string.split(c,'/')
                    num = 0.0; den = 0.0; v1 = []; v2 = []
                    for e in conditions[c1]: num += entry[e]; v1.append(entry[e])
                    for e in conditions[c2]: den += entry[e]; v2.append(entry[e])
                    num /= len(conditions[c1])
                    num = max(num,minval)
                    den /= len(conditions[c2])
                    den = max(den,minval)
                    entry[c] = num / den
                    if entry[c] >= 1: entry['minabs%s' % c] = max(min(v1),minval) / max(max(v2),minval)
                    else: entry['minabs%s' % c] = max(max(v1),minval) / max(min(v2),minval)
            self.printLog('\r#ENR','Calculation of Tag enrichment complete.')
            mdb.saveToFile()
        except: self.errorLog('%s.enrichment error' % self); raise
#########################################################################################################################
    ### <5> ### EST Hit processing Methods                                                                              #
#########################################################################################################################
    def enrFilter(self,mdb):    ### Filters results according to enrichment
        '''Filters results according to enrichment.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            enrcut = self.getNum('EnrCut'); pwenr = self.getNum('PWEnr')
            if enrcut <= 1.0 and pwenr <= 1.0: return
            self.setBasefile('%s.filtered' % (self.basefile()))
            goodtag = []; goodseq = []; seq2tag = {}; tag2seq = {}
            ### ~ [1] ~ Identify enriched Tags and Sequences from initial data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = mdb.entryNum()
            for entry in mdb.entries():
                self.progLog('\r#ENR','Filtering on Tag enrichment: %.2f%% ' % (ex/etot)); ex += 100.0
                enriched = False
                for c in self.list['Comparisons']:
                    if pwenr > 0:
                        if entry[c] >= 1.0 and entry['minabs%s' % c] < pwenr: continue
                        elif entry[c] < 1.0 and (1.0 / entry['minabs%s' % c]) < pwenr: continue
                    if entry[c] >= enrcut or entry[c] <= 1.0 / enrcut: enriched = True; break
                if enriched: goodtag.append(entry['Tag'])
                if entry['Seq']: tag2seq[entry['Tag']] = string.split(entry['Seq'],';')
                else: tag2seq[entry['Tag']] = []
                for seq in tag2seq[entry['Tag']]:
                    if seq not in seq2tag: seq2tag[seq] = []
                    seq2tag[seq].append(entry['Tag'])
                    if enriched and seq not in goodseq: goodseq.append(seq)
            ### ~ [2] ~ Expand using Seq-Tag links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            expandseq = goodseq[0:]
            expand = self.opt['Expand']; xx = 1
            while expand and expandseq:
                self.progLog('\r#ENR','Filtering on Tag enrichment: %s     ' % (rje.iLen(expandseq)))
                seq = expandseq.pop(0)
                for tag in seq2tag[seq]:
                    if tag in goodtag: continue
                    goodtag.append(tag)
                    for newseq in string.split(mdb.data()[tag]['Seq'],';'):
                        if newseq in goodseq: continue
                        goodseq.append(newseq); expandseq.append(newseq)
            self.printLog('\r#ENR','Filtered on Tag enrichment: %s Tags; %s Seqs' % (rje.iLen(goodtag),rje.iLen(goodseq)))
            ### ~ [3] ~ Reduce dataset and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.dropEntriesDirect('Tag',goodtag,inverse=True)
            mdb.saveToFile()
        except: self.errorLog('%s.enrFilter error' % self)
#########################################################################################################################                    
    def processHits(self,tdb=None):  ### Processes Hit ESTs into a new table
        '''Processes Hit ESTs into a new table.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['CoreBaseFile'].lower() not in ['','none']: self.info['Basefile'] = self.info['CoreBaseFile']
            pdb = self.db().addEmptyTable('HitSeq',['Seq','Tags'],['Seq'])
            pdata = pdb.dict['Data']
            seqlist = None
            if rje.checkForFile(self.info['SeqIn']):
                pdb.addField('Desc','Seq')
                seqlist = self.seqlist(); seqlist.seq = []
                estfile = '%s.gfessa.hits.fas' % self.info['Basefile']
                blast = rje_blast.blastObj(self.log,self.cmd_list,type='Dev')
                blast.formatDB(fasfile=self.info['SeqIn'],force=False,protein=False)
            ### ~ [1] ~ Identify unique Sequence hits and pull out sequences if file given ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tx = 0.0; ttot = tdb.entryNum()
            for tag in tdb.datakeys():
                entry = tdb.data()[tag]
                self.progLog('\r#HIT','Extracting Hit data from Tag matches: %.2f%%' % (tx/ttot)); tx += 100.0
                taghits = []
                if self.stat['Mismatch'] >= 0 and self.opt['BestMatch']:
                    mm = 0
                    while mm <= self.stat['Mismatch']:
                        if entry['Seq']: taghits = string.split(entry['Seq'],';'); break
                        else: mm += 1
                else: taghits = string.split(entry['Seq'],';')
                for hit in taghits:
                    if not hit: continue
                    if hit not in pdata:
                        pdata[hit] = {'Seq':hit,'Tags':[]}
                        if seqlist:
                            if blast.oldBLAST(): hitseq = seqlist.seqFromFastaCmd(hit)
                            else: hitseq = seqlist.seqFromBlastDBCmd(hit)
                            pdata[hit]['Desc'] = hitseq.info['Description']
                    if tag not in pdata[hit]['Tags']: pdata[hit]['Tags'].append(tag)
            self.printLog('\r#HIT','Extracting Hit data from Tag matches: %s hits' % rje.iLen(pdata))
            if seqlist: seqlist.saveFasta(seqfile=estfile)  # Considering adding back the options for n-trimming etc.
            ### ~ [2] ~ Group Sequences based on Tag/Seq clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            all_seq = rje.sortKeys(pdata); sx = 0.0; stot = len(all_seq)
            all_tags = tdb.datakeys()[0:]
            groups = []
            while all_seq:
                self.progLog('\r#GRP','Clustering Seq/Tag combos: %.2f%%' % (sx/stot))
                nextgroup = []; addme = [all_seq[0]]
                while addme:
                    add = addme.pop(0)
                    nextgroup.append(add)
                    if add in all_seq:
                        all_seq.remove(add); sx += 100.0
                        for tag in pdata[add]['Tags']:
                            if tag not in nextgroup and tag not in addme: addme.append(tag)
                    elif add in all_tags:
                        all_tags.remove(add)
                        for seq in string.split(tdb.data()[add]['Seq'],';'):
                            if seq not in all_seq: continue     # Done already, or not to be done
                            if seq not in nextgroup + addme: addme.append(seq)
                groups.append(nextgroup[0:])
            self.printLog('\r#GRP','Clustered %s Seq & %s Tag: %s combos; %s Tags w/o hits' % (rje.iLen(pdata),tdb.entryNum(),rje.iLen(groups),rje.iLen(all_tags)))
            ## ~ [2a] ~ Make Cluster Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb = self.db().addEmptyTable('Cluster',['Cluster','CSeqN','CTagN','CSeq','CTag'],['Cluster'])
            tdb.addField('Cluster',evalue=0); pdb.addField('Cluster'); gx = 0
            tdb.addField('CSeqN',evalue=0); pdb.addField('CSeqN')
            tdb.addField('CTagN',evalue=0); pdb.addField('CTagN')
            for group in groups:
                gx += 1; tags = []; seqs = []
                for member in group:
                    if member in tdb.data(): tags.append(member)
                    elif member in pdata: seqs.append(member)
                    else: self.errorLog('Cannot find %s in Seq or Tag?!' % member,printerror=False)
                gdb.addEntry({'Cluster':gx,'CSeqN':len(seqs),'CTagN':len(tags),'CSeq':string.join(seqs,'; '),'CTag':string.join(tags,'; ')})
                for tag in tags: tdb.data()[tag]['Cluster'] = gx; tdb.data()[tag]['CSeqN'] = len(seqs); tdb.data()[tag]['CTagN'] = len(tags)
                for tag in seqs: pdb.data()[tag]['Cluster'] = gx; pdb.data()[tag]['CSeqN'] = len(seqs); pdb.data()[tag]['CTagN'] = len(tags)
            tdb.saveToFile('%s.gfessa.tags.tdt' % self.info['Basefile'])
            ### ~ [3] ~ Calculate Mean Values for Enrichment etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(pdata)
            for e in self.list['Experiments'] + self.list['Comparisons']: pdb.addField(e)
            cluster_enr = {}
            for comp in self.list['Comparisons']:
                pdb.addField('C|%s' % comp); cluster_enr[comp] = {}
                for gx in range(len(groups)): cluster_enr[comp][gx+1] = []
            ## ~ [3a] ~ Enrichment for sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for seq in pdata:
                self.progLog('\r#MEAN','Calculating mean values for hits: %.2f%%' % (sx/stot)); sx += 100.0
                entry = pdata[seq]
                entry['Tags'].sort()
                for e in self.list['Experiments'] + self.list['Comparisons']:
                    if len(entry['Tags']) > 1:
                        entry[e] = []
                        for tag in entry['Tags']: entry[e].append(tdb.data()[tag][e])
                        if e in self.list['Experiments']: entry[e] = rje.meansd(entry[e])[0]
                        else: entry[e] = rje.geoMean(entry[e])
                    else:
                        tag = entry['Tags'][0]
                        entry[e] = tdb.data()[tag][e]
                cluster = entry['Cluster']
                for comp in self.list['Comparisons']:
                    for tag in entry['Tags']: cluster_enr[comp][cluster].append(tdb.data()[tag][comp])
                entry['Tags'] = string.join(entry['Tags'],'; ')
            self.printLog('\r#MEAN','Calculating mean values for hits complete.')
            ## ~ [3b] ~ Enrichment for clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#MEAN','Calculating cluster means...')
            for comp in self.list['Comparisons']:
                gdb.addField('C|%s' % comp)
                for gx in range(len(groups)):
                    cluster_enr[comp][gx+1] = rje.geoMean(cluster_enr[comp][gx+1])
                    gdb.data()['%d' % (gx+1)]['C|%s' % comp] = cluster_enr[comp][gx+1]
            gdb.saveToFile('%s.gfessa.clusters.tdt' % self.info['Basefile'])
            ## ~ [3c] ~ Update hit table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for seq in pdata:
                entry = pdata[seq]
                cluster = entry['Cluster']
                for comp in self.list['Comparisons']: entry['C|%s' % comp] = cluster_enr[comp][cluster]
            self.printLog('\r#MEAN','Added mean cluster values to hit table.')
            pdb.saveToFile('%s.gfessa.hits.tdt' % self.info['Basefile'])
        except: self.errorLog('%s.processHits error' % self)
#########################################################################################################################
### End of SECTION II: GFESSA Class                                                                                     #
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
    try: GFESSA(mainlog,cmd_list).run()

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
