#!/usr/bin/python

# GOPHER - Generation of Orthologous Proteins from High-Throughput Estimation of Relationships
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / Biological Sciences, University of Southamtpon, SO16 7PX, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

'''
Program:      GOPHER
Description:  Generation of Orthologous Proteins from Homology-based Estimation of Relationships
Version:      3.4.1
Last Edit:    31/01/15
Citation:     Davey, Edwards & Shields (2007), Nucleic Acids Res. 35(Web Server issue):W455-9. [PMID: 17576682]
Copyright (C) 2005 Richard J. Edwards - See source code for GNU License Notice

Function:
    This script is designed to take in two sequences files and generate datasets of orthologous sequence alignments.
    The first `seqin` sequence set is the 'queries' around which orthologous datasets are to be assembled. This is now
    optimised for a dataset consisting of one protein per protein-coding gene, although splice variants should be dealt
    with OK and treated as paralogues. This will only cause problems if the postdup=T option is used, which restricts
    orthologues returned to be within the last post-duplication clade for the sequence.

    The second `orthdb` is the list of proteins from which the orthologues will be extracted. The `seqin` sequences are
    then BLASTed against the orthdb and processed (see below) to retain putative orthologues using an estimation of the
    phylogenetic relationships based on pairwise sequences similarities.

    NB. As of version 2.0, gopher=FILE has been replaced with seqin=FILE for greater rje python consistency. The `allqry`
    option has been removed. Please cleanup the input data into a desired non-redundant dataset before running GOPHER.
    (In many ways, GOPHER's strength is it's capacity to be run for a single sequence of interest rather than a whole
    genome, and it is this functionality that has been concentrated on for use with PRESTO and SLiM Pickings etc.) The
    output of statistics for each GOPHER run has also been discontinued for now but may be reintroduced with future
    versions. The `phosalign` command (to produce a table of potential phosphorylation sites (e.g. S,T,Y) across
    orthologues for special conservation of phosphorylation prediction analyses) has also been discontinued for now.

    Version 2.1 has tightened up on the use of rje_seq parameters that were causing trouble otherwise. It is now the
    responsibility of the user to make sure that the orthologue database meets the desired criteria. Duplicate accession
    numbers will not be tolerated by GOPHER and (arbitrary) duplicates will be deleted if the sequences are the same, or
    renamed otherwise. Renaming may cause problems later. It is highly desirable not to have two proteins with the same
    accession number but different amino acid sequences. The following commands are added to the rje_seq object when input
    is read: accnr=T unkspec=F specnr=F gnspacc=T. Note that unknown species are also not permitted.

    Version 3.0 has improved directory organisation for multi-species and multi-orthdb GOPHER runs on the same system, in
    line with the bioware.ucd.ie webserver. Additional data cleanup has been added too. (NB. Experimental Sticky GOPHER
    runs will not work with Organise=T.) The output directory is now set at the highest level by gopherdir=PATH. If
    organise=T then a subdirectory will be created within this directory named after the `orthdb` (path and extension
    stripped) and each species will have its own subdirectory within this in turn. By default, these are named using the
    UniProt species codes, which are read from the input sequences. To use TaxIDs, the unispec=FILE option must be used.
    This should point to a file that has one line per UniProt species code, starting with the species code and ending
    with :TaxID: where TaxID is a number. The file should be sorted by species code.

Processing:
    The process for dataset assembly is as follows for each protein :

    1. BLAST against orthdb [`orthblast`]
    * BLASTs saved in BLAST/AccNum.blast

    2. Work through BLAST hits, indentifying paralogues (query species duplicates) and the closest homologue from each
    other species. This involves a second BLAST of the query versus original BLAST hits (blaste=10, no complexity
    filter). The best sequence from each species is kept, i.e. the one with the best similarity to the query and not part
    of a clade with any paralogue that excludes the query. (If postdup=T, the hit must be in the query's post duplication
    clade.) In addition hits:  [`orthfas`]
    * Must have minimum identity level with Query
    * Must be one of the 'good species' [goodspec=LIST]
    * Save reduced sequences as ORTH/AccNum.orth.fas
    * Save paralogues identified (and meeting minsim settings) in PARA/AccNum.para.fas

    3. Align sequences with MUSCLE  [`orthalign`]
    * ALN/AccNum.orthaln.fas

    4. Generate an unrooted tree with (ClustalW or PHYLIP)  [`orthtree`]
    * TREE/AccNum.orth.nsf

    Optional paralogue/subfamily output:  (These are best not used with Force=T or FullForce=T)
    2a. Alignment of query protein and any paralogues >minsim threshold (paralign=T/F). The parasplice=T/F controls
    whether splice variants are in these paralogue alignments (where identified using AccNum-X notation).
    * PARALN/AccNum.paraln.fas

    2b. Pairwise combinations of paralogues and their orthologues aligned, with "common" orthologues removed from the
    dataset, with a rooted tree and group data for BADASP analysis etc. (parafam=T)
    * PARAFAM/AccNum+ParaAccNum.parafam.fas
    * PARAFAM/AccNum+ParaAccNum.parafam.nsf
    * PARAFAM/AccNum+ParaAccNum.parafam.grp

    2c. Combined protein families consisting of a protein, all the paralogues > minsim and all orthologues for each in a
    single dataset. Unaligned. (gopherfam=T)
    * SUBFAM/AccNum.subfam.fas

    NB. The subfamily outputs involve Gopher calling itself to ensure the paralogues have gone through the Gopher
    process themselves. This could potentially cause conflict if forking is used.

Commandline:
    ### Basic Input/Output ###
    seqin=FILE      : Fasta file of 'query' sequences for orthology discovery (over-rides uniprotid=LIST) []
    uniprotid=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.fas and use as seqin=FILE. []
    orthdb=FILE     : Fasta file with pool of sequences for orthology discovery []. Should contain query sequences.
    startfrom=X     : Accession Number / ID to start from. (Enables restart after crash.) [None]
    dna=T/F         : Whether to analyse DNA sequences (not optimised) [False]

    ### GOPHER run control parameters ###
    orthblast   : Run to blasting versus orthdb (Stage 1).
    orthfas     : Run to output of orthologues (Stage 2). 
    orthalign   : Run to alignment of orthologues (Stage 3).
    orthtree    : Run to tree-generation (Stage 4). [default!]

    ### GOPHER Orthologue identifcation Parameters ###
    postdup=T/F     : Whether to align only post-duplication sequences [False]
    reciprocal=T/F  : Use Reciprocal Best Hit method instead of standard GOPHER approach [False]
    fullrbh=T/F     : Whether RBH method should run BLAST searches for all potential RBH orthologues [False]
    compfilter=T/F  : Whether to use complexity filter and composition statistics for *initial* BLAST [True]
    minsim=X        : Minimum %similarity of Query for each "orthologue" [40.0]
    simfocus=X      : Style of similarity comparison used for MinSim and "Best" sequence identification [query]
        - query = %query must > minsim (Best if query is ultimate focus and maximises closeness of returned orthologues)
        - hit = %hit must > minsim (Best if lots of sequence fragments are in searchdb and should be retained)
        - either = %query > minsim OR %hit > minsim (Best if both above conditions are true)
        - both = %query > minsim AND %hit > minsim (Most stringent setting)
    gablamo=X       : GABLAMO measure to use for similarity measures [Sim]
        - ID = %Identity (from BLAST)
        - Sim = %Similarity (from BLAST)
        - Len = %Coverage (from BLAST)
        - Score = BLAST BitScore (Reciprocal Match only)
    goodX=LIST      : Filters where only sequences meeting the requirement of LIST are kept.
                      LIST may be a list X,Y,..,Z or a FILE which contains a list [None]
                        - goodacc  = list of accession numbers
                        - goodseq  = list of sequence names
                        - goodspec = list of species codes
                        - gooddb   = list of source databases
                        - gooddesc = list of terms that, at least one of which must be in description line
    badX=LIST       : As goodX but excludes rather than retains filtered sequences

    ### Additional run control options ###
    repair=T/F      : Repair mode - replace previous files if date mismatches or files missing.
                      (Skip missing files if False) [True]
    force=T/F       : Whether to force execution at current level even if results are new enough [False]
    fullforce=T/F   : Whether to force current and previous execution even if results are new enough [False]
    dropout=T/F     : Whether to "drop out" at earlier phases, or continue with single sequence [False]
    ignoredate=T/F  : Ignores the age of files and only replaces if missing [False]
    savespace=T/F   : Save space by deleting intermediate blast files during orthfas [True]
    maxpara=X       : Maximum number of paralogues to consider (large gene families can cause problems) [50]
    oldblast=T/F    : Run with old BLAST rather than BLAST+ [False]

    ### Additional Output Options ###
    runpath=PATH    : Directory from which to run GOPHER. (NB. Will look for input here unless full paths given) [./]
    gopherdir=PATH  : Parent directory for output of files [./]
    organise=T/F    : Output files according to orthdb an species (code or TaxaID - need conversion) [True]
    unispec=FILE    : UniProt species file for TaxaID conversion. Will use TaxaID instead of Species Code if given [None]
    alnprog=X       : Choice of alignment program to use (clustalw/clustalo/muscle/mafft/fsa) [clustalo]
    orthid=X        : File identifier (Lower case) for orthology files [orth]
    paralign=T/F    : Whether to produce paralogue alignments (>minsim) in PARALN/ (assuming run to orthfas+) [False]
    parasplice=T/F  : Whether splice variants (where identified) are counted as paralogues [False]
    parafam=T/F     : Whether to paralogue paired subfamily alignments (>minsim) (assuming run to orthfas+) [False]
    gopherfam=T/F   : Whether to combined paralogous gopher orthologues into protein families (>minsim) (assuming run to orthfas+) [False]
    sticky=T/F      : Switch on "Sticky Orthologous Group generation" [False] *** Experimental ***
    stiggid=X       : Base for Stigg ID numbers [STIGG]
'''
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, gc, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_sequence, rje_tree, rje_uniprot
import rje_blast_V2 as rje_blast
import rje_dismatrix_V2 as rje_dismatrix    #!# Check that this is OK! #!#
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 3.0 - See archived GOPHER 1.9 and gopher_V2 2.9 for history and obselete options.
    # 3.0 - Added organise=T/F and gopherdir=PATH for improved file organisation. Tightened savespace. 
    # 3.0 - Added compfilter=T/F for improved complexity filter and composition statistics control for *initial* BLAST.
    # 3.0 - Changed default tree extension to *.nwk for compatibility with MEGA. Deleted _phosAlign() method.
    # 3.0 - Added orthology ID option and alignment program to customise output further.
    # 3.1 - Added full reciprocal best hit method. (fullrbh=T/F)
    # 3.2 - Minor tweak to prevent unwanted directory generation for programs using existing GOPHER alignments.
    # 3.3 - Added rje_blast_V2 to use BLAST+. Run with oldblast=T to stick with old NCBI BLAST.
    # 3.4 - Fixed FullRBH paralogue duplication issue.
    # 3.4.1 - Fixed stripXGap issue. (Why was this being implemented anyway?). Added REST output.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add a reciprocal best hit method
    # -- [Y] : Modify to use any stat including BLAST score for comparison
    # -- [Y] : Add a orthid=X option for GOPHER vs Sticky vs PostDup vs MBH.
    # -- [ ] : Add option to read in distances from file.
    # -- [Y] : Update Reciprocal best hit method to be intensive, i.e. read/run the BLAST stage of GOPHER for each species.
    # [ ] : Consider re-introduction of GOPHER Statistics
    # [ ] : Consider addition of a menu front-end
    # [ ] : Add or scrap orthanc (GASP) run mode. (Replaced by gopherfam etc?)
    # [ ] : Consider automated BADASP analysis?
    # [ ] : Add automated GABLAM table generation for GOPHER alignments?
    # [Y] : Output files according to species (code or TaxaID - need conversion) (separate run directories)
    # [Y] : Output files into directory named after orthdb.fas
    # [Y] : Organise=T/F and GopherDir=PATH/ settings for organising output files. (V3.0)
    # [ ] : Improve handling and cleanup of unnecessary files. (Document too.)
    # [ ] : Give main run() method a tidy at some point.
    # [ ] : Check/Fix Sticky GOPHER runs: will not work with Organise=T.
    # [Y] : Should composition statistics and complexity filter always be on for the BLAST?
    # [Y] : Add and check rje_blast_V2 and BLAST+
    # [Y?]: Tidy up forking use of RunPath so that log files do not end up in RunPath within RunPath!
    # [ ] : Add generation of a summary table of pairwise orthologues following (OrthFas) run.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('GOPHER', '3.4.1', 'January 2015', '2005')
    description = 'Generation of Orthologous Proteins from Homology-based Estimation of Relationships'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please cite SLiMDisc webserver paper. (Davey, Edwards & Shields 2007)']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Show RJE_SEQ commandline options (redundancy etc.)?'): out.verbose(-1,4,text=rje_seq.__doc__)
            if rje.yesNo('Show RJE_TREE commandline options (used for tree-making)?'): out.verbose(-1,4,text=rje_tree.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Gopher Class:                                                                                           #
#########################################################################################################################
class Gopher(rje.RJE_Object):     
    '''
    Gopher main controller class. Author: Rich Edwards (2005). This class contains a lot of forking options for use on
    certain machines with multiple processors.    

    Info:str
    - Name = Name of Input sequence file (gopher=FILE)
    - OrthDB = Name of sequence database to fnd orthologues in (orthdb=FILE)
    - StartFrom = AccNum/ID to start running from
    - StiggID = Base for Stigg ID numbers [STIGG]
    - Gopher = Alternative input file name - use only if seqin in ['','None']
    - GopherDir = Parent directory for output of files [./]
    - UniSpec = UniProt species file for TaxaID conversion. Will use TaxaID instead of Species Code if given [None]
    
    Opt:boolean
    - DNA = Whether to analyse DNA sequences (not optimised) [False]
    - DropOut = Whether to "drop out" at earlier phases, or continue with single sequence [False]
    - FullForce = Whether to force current and previous execution even if results are new enough [False]
    - IgnoreDate = Ignores the age of files and only replaces if missing [False]
    - Organise = Output files according to orthdb an species (code or TaxaID - need conversion) [True]
    - Sticky = Switch on "Sticky Orthologous Group generation" [False]

    Stat:numeric

    List:list
    - UniprotID=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.fas and use as seqin=FILE. []

    Dict:dictionary
    - Output = REST server output. Only works for single sequence runs.
    - TaxaCode = Dictionary of {Species Code: TaxaID} if Organise and UniSpec are being used.

    Obj:RJE_Objects
    - SeqIn = SeqList object for handling input sequences
    - BLAST = BlastRun object for handling BLAST searches
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','OrthDB','StartFrom','Gopher','StiggID','GopherDir','UniSpec']
        self.statlist = []
        self.optlist = ['FullForce','IgnoreDate','DNA','Child','Sticky','DropOut','Organise','OldBLAST']
        self.listlist = ['UniprotID']
        self.dictlist = ['Output','TaxaCode']
        self.objlist = ['SeqIn','BLAST']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'GopherDir':rje.makePath('./')})  #'Name':'seqin.fas','OrthDB':'nr.fas'})
        self.setOpt({'Organise':True,'OldBLAST':False})
        ### Other Attributes ###
        self._setForkAttributes()   
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### ~ [1] General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
                self._generalCmd(cmd)
                self._forkCmd(cmd)  
                ### ~ [2] Specific Class Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
                self._cmdRead(cmd,type='file',att='Name',arg='seqin')
                self._cmdReadList(cmd,'opt',['FullForce','IgnoreDate','DNA','Sticky','DropOut','Organise','OldBLAST'])
                self._cmdReadList(cmd,'file',['OrthDB','Gopher','UniSpec'])  
                self._cmdReadList(cmd,'path',['GopherDir'])  
                self._cmdReadList(cmd,'list',['UniprotID'])
                self._cmdReadList(cmd,'info',['StartFrom','StiggID'])
            except: self.log.errorLog('Problem with cmd [%s]:' % cmd)
        self.cmd_list += ['seqin=None']     # Other objects should not use seqin!
#########################################################################################################################
    ### <2> ### Main Forker Method                                                                                      #
#########################################################################################################################
    def setupBlast(self):     ### Sets up and self.obj['BLAST'] object                                              #V3.3
        '''Sets up a self.obj['BLAST'] object.'''
        try:### ~ [1] Setup BLAST object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gopherblast = rje_blast.blastObj(log=copy.deepcopy(self.log),cmd_list=['blastv=1000','v=0']+self.cmd_list,type='Dev')
            if self.dev(): self.deBug(gopherblast)
            gopherblast.setStr({'Type':'blastp','DBase':self.getStr('OrthDB')})
            gopherblast.setInt({'HitAln':0})
            #!# Temporary fix for hitaln:0 problem crashing BLAST. Why? Possible change to default gopher=T/F?
            #gopherblast.setInt({'HitAln':gopherblast.getInt('OneLine')})
            if self.opt['DNA']: gopherblast.setStr({'Type':'blastn'})
            ### ~ [2] FormatDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(self.getStr('OrthDB')):
                self.errorLog('BLAST Database file %s missing!' % self.getStr('OrthDB'),printerror=False)
                sys.exit(1)
            gopherblast.formatDB(protein=not self.getBool('DNA'),force=False,checkage=not self.getBool('IgnoreDate'),details=True)
            self.obj['BLAST'] = gopherblast
            return gopherblast
        except:
            self.errorLog('Disaster with Gopher.setupBlast(%s)' % self.getStr('OrthDB'))
            raise
#########################################################################################################################
    def setMode(self):  ### Returns mode to run from self.cmd_list                                                  #V3.0
        if self.opt['SoapLab'] or self.opt['Sticky']: return 'orthfas'
        mode = gophermodes[-1]    # By default, run all
        modes = gophermodes[1:] + ['phosalign']
        for cmd in self.cmd_list:
            if cmd in modes: mode = cmd  # Command found
            elif cmd == 'orthaln': mode = 'orthalign'
        return mode
#########################################################################################################################
    def setupSeqIn(self):   ### Sets up SeqIn, using UniprotID list if appropriate
        '''Sets up SeqIn, using UniprotID list if appropriate.'''
        if self.list['UniprotID'] and not rje.exists(self.getStr('Name')):
            if self.getStrLC('Basefile'): fasfile = '%s.fas' % self.baseFile(runpath=True)
            else: fasfile = '%s%s.fas' % (self.getStr('RunPath'),self.prog().lower())
            datfile = fasfile[:-3] + 'dat'
            unicmd = self.cmd_list + ['datout=%s' % datfile]
            rje_uniprot.UniProt(self.log,unicmd).run()
            if not rje.exists(datfile): return False
            seqcmd = ['gnspacc=T','usecase=T','accnr=F','seqnr=F'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            seqcmd += ['seqin=%s' % datfile,'seqout=%s' % fasfile,'reformat=fas']
            self.setStr({'Name':fasfile})
            rje_seq.SeqList(self.log,seqcmd).run()
        #if rje.exists(seqlist.name()): self.dict['Output']['seqin'] = seqlist.name()
        #else: self.dict['Output']['seqin'] = seqlist.fasta()
        return True
#########################################################################################################################
    def run(self,runpath=None):      ### Main Gopher run, including forking                                         #V3.0
        '''Main Gopher run, including forking.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            startpath = os.path.abspath(os.curdir)
            if not runpath: runpath = self.info['RunPath']
            if os.path.abspath(runpath) != startpath:
                self.log.printLog('#PATH','Running from %s' % os.path.abspath(runpath))
                rje.mkDir(self,runpath)
                os.chdir(runpath)
            ## ~ [1a] Forking Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            forkx = self.stat['Forks']      # Number of forks to have running at one time
            if self.opt['Win32'] or forkx < 1: self.opt['NoForks'] = True
            forks = []      # List of active forks
            forkaccs = {}   # Dictionary of Accession numbers for fork PIDs {pid:accnum}
            killforks = self.stat['KillForks']      # Time in seconds to wait after main thread has apparently finished
            killtime = time.time()
            ## ~ [1b] Gopher Variables Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.setupSeqIn(): raise ValueError('Unable to generate input from UniprotID list.')
            seqfile = self.info['Name']
            if seqfile.lower() in ['','none']: seqfile = self.info['Gopher']
            if not os.path.exists(seqfile):
                self.log.errorLog('Sequence file %s missing!' % seqfile,printerror=False)
                sys.exit(1)                        
            basefile = rje.baseFile(seqfile)
            startfrom = self.info['StartFrom']   # First sequence to run through gopher
            if startfrom == 'None': startfrom = None
            rejects = []    # Rejected AccNum (happens when multiple isoforms present in a file)
            mode = self.setMode()
            ## ~ [1c] Gopher BLAST Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gopherblast = self.setupBlast()
            gopherblast.setStr({'DBase':self.info['OrthDB']})      # Altered by formatDB?
            ## ~ [1d] Sticky Mode Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.stickySetup(seqfile,startfrom) #!# Experimental #!#

            ### ~ [1] Forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
            gseqx = 0   # Counter for GOPHER
            gopherseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=T','unkspec=F','specnr=F','gnspacc=T'])
            SEQFILE = open(seqfile, 'r')
            readingseq = True   # Still reading in sequences from file
            lastline = ''
            while readingseq or len(forks):
                if self.opt['Sticky']: readingseq = self.list['StiGGing']
                ## ~ [1a] Read sequence and start new fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while readingseq and (len(forks) < forkx or self.opt['NoForks']):     # Add more forks
                    killtime = time.time()  # Reset killtime - still doing stuff
                    nextseq = None
                    if readingseq:  # Get next sequence from file
                        if self.opt['Sticky']:
                            if not self.list['StiGGing']: readingseq = False; break
                            stigg = self.list['StiGGing'].pop(0)
                            self.list['StiGGed'].append(stigg)
                            nextseq = gopherseq.seqFromFastaCmd(id=stigg,dbase=gopherblast.getStr('DBase'))
                            if not nextseq: continue    #!# What are consequences of this error?! #!#
                        else: (nextseq, lastline) = gopherseq.nextFasSeq(fileobject=SEQFILE,lastline=lastline)
                    else: break  # Not ready!
                    if nextseq:
                        gseqx +=1 
                        gopherseq.seq = [nextseq]
                        # Skip to start if appropriate
                        self._logTransfer(accnum=nextseq.info['AccNum'])    # Cleans up logs generated by previous killings!
                        if startfrom:
                            if nextseq.shortName().find('%s ' % startfrom) >= 0: startfrom = None
                            else: self.stiGGMe(seq=nextseq); continue
                        # Check rejects
                        if nextseq.shortName() in self.gopherRejects():  # Skip
                            self.verbose(0,3,'\n *** Ignoring %s. (Rejected) ***' % nextseq.shortName(),2)
                            continue
                        # Check for filtering (especially goodspec)
                        gopherseq.autoFilter(cmd_list=self.cmd_list+['filterout=None','seqout=None'])
                        if nextseq not in gopherseq.seq: continue   # Filtered out. Message should be in Log #
                        # Prep Sequence for run
                        gopherdir = self.gopherDir(nextseq)
                        # Check BLAST object not messed up #
                        if self.opt['DNA']: gopherblast.setStr({'Type':'blastn'})
                        else: gopherblast.setStr({'Type':'blastp'})
                        # Add new fork (or not!)
                        vtext = 'Sequence %s: %s' % (rje.integerString(gseqx), nextseq.shortName())
                        if self.opt['Sticky']: vtext = '%s (%s StiGGing)' % (vtext,rje.integerString(len(self.list['StiGGing'])))
                        if self.opt['NoForks']:
                            self.verbose(0,2,'\n\n *** GOPHER %s ***' % vtext,1)
                            if self.getStrLC('Rest'): self.restSetup()
                            solo_gopher = GopherFork(log=self.log,cmd_list=self.cmd_list+['i=-1','errorlog=None'])
                            solo_gopher.dict['Output'] = self.dict['Output']
                            solo_gopher.setInfo({'Name':nextseq.shortName(),'OutPath':gopherdir})
                            solo_gopher.obj['Sequence'] = nextseq
                            solo_gopher.obj['BLAST'] = gopherblast
                            solo_gopher.obj['Parent'] = self
                            gopherblast.log = self.log
                            solo_gopher.run(mode)
                            solo_gopher.restOutputs()
                            self.stiGGMe(seq=nextseq)
                        else:
                            self.verbose(0,2,'GOPHER Fork - %s...' % vtext,1)
                            forkcmd = self.cmd_list + ['i=-1','log=%s.log' % nextseq.info['AccNum'],'errorlog=None','runpath=./']
                            newpid = os.fork() 
                            if newpid == 0: # child
                                self.opt['Child'] = True
                                forked_gopher = GopherFork(cmd_list=forkcmd)
                                forked_gopher.setInfo({'Name':nextseq.shortName(),'OutPath':gopherdir})
                                forked_gopher.obj['Sequence'] = nextseq
                                forked_gopher.obj['BLAST'] = gopherblast
                                forked_gopher.obj['Parent'] = self
                                forked_gopher.run(mode)
                                sys.exit()    # Exit process 
                            elif newpid == -1: self.log.errorLog('Problem forking %s.' % new_fork_id,printerror=False)  # Error!
                            else:
                                forkaccs[newpid] = nextseq.info['AccNum']
                                forks.append(newpid)    # Add fork to list 
                    else: readingseq = False    # No more sequences to process - just deal with existing forks
            
                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.opt['NoForks']: time.sleep(1)       # Sleep for 1s 
                forklist = self._activeForks(forks)
                if len(forklist) != len(forks):
                    self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    forks = forklist[0:]
                    for pid in forkaccs.keys():   # Go through current forks
                        if pid not in forks:
                            fork = forkaccs.pop(pid)
                            killtime = time.time()  # Reset killtime - still doing stuff
                            self.verbose(0,2,' => GOPHER Fork %s Finished! Transfering log details to %s.' % (fork,self.log.info['LogFile']),1)
                            self._logTransfer(accnum=fork)
                            self.stiGGMe(acc=fork)
                self.verbose(3,3,'End of a Cycle.',2)

                ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if (time.time() - killtime) > killforks:
                    self.verbose(0,1,'\n%d seconds of main thread inactivity. %d forks still active!' % (killforks,len(forks)),1)
                    for fork in forks:
                        self.verbose(0,2,' => GopherFork %s, PID %d still Active!' % (forkaccs[fork],fork),1)
                    if self.stat['Interactive'] < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                    elif rje.yesNo('Kill hanging forks?'):
                        for fork in forks:
                            self.log.printLog('#KILL','Killing GopherFork %s, PID %d.' % (forkaccs[fork],fork))
                            os.system('kill %d' % fork)
                    else: killtime = time.time()

            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SEQFILE.close()
            if len(forks) > 0: self.log.errorLog('%d Gopher Forks still active after %d seconds of mainthread inactivity' % (len(forks),killforks),quitchoice=True,printerror=False)
            elif not self.opt['NoForks']: self.verbose(0,1,'Gopher Forks have finished.',2)
            if os.path.abspath(runpath) != startpath:
                self.log.printLog('#PATH','Reverting to %s' % startpath)
                os.chdir(startpath)

            ### ~ [3] Sticky ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
            self.stiGG(gopherseq,gopherblast)
        except SystemExit:
            if self.opt['Child']: sys.exit()
            else: raise ValueError
        except: self.log.errorLog('Error in Gopher run()'); raise
#########################################################################################################################
    def _activeForks(self,pidlist=[]):   ### Checks Process IDs of list and returns list of those still running.
        '''
        Checks Process IDs of list and returns list of those still running.
        >> pidlist:list of integers = Process IDs
        '''
        try:
            ### Checks PIDs and removes any that have finished running ###
            oldpids = pidlist[0:]
            for pid in oldpids[0:]:
                (newpid,exitcode) = os.waitpid(pid,os.WNOHANG)
                while exitcode > 255: exitcode -= 256
                if newpid == pid and exitcode == 0: oldpids.remove(pid)
                elif newpid == pid:
                    oldpids.remove(pid)
                    self.log.errorLog('WARNING!: PID %d returned with exit code %d.' % (pid,exitcode),printerror=False)
            return oldpids
        except:
            self.log.errorLog('Error in _activeForks(%s)' % pidlist)
            raise   
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _logTransfer(self,accnum=''):   ### Transfers details from accnum.log to self.log
        '''
        Transfers details from accnum.log to self.log.
        accnum:str = leader for accnum.log
        '''
        try:
            if os.path.exists('%s.log' % accnum):   
                open(self.log.info['LogFile'],'a').write(open('%s.log' % accnum,'r').read())
                os.unlink('%s.log' % accnum)
        except: self.errorLog('Problem transfering %s log details to %s.' % (accnum,self.log.log_file))       
#########################################################################################################################
    def gopherRejects(self):  ### Returns list of rejects from gopher.rejects
        '''Returns list of rejects from gopher.rejects.'''
        if os.path.exists('gopher.rejects'): return self.loadFromFile('gopher.rejects', v=1, chomplines=True)
        return []
#########################################################################################################################
    def gopherDir(self,nextseq,mkdir=True):    ### Returns Parent directory for GopherFork outputs                  #V3.1
        '''Returns Parent directory for GopherFork outputs.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gopherdir = rje.makePath('./')
            if self.info['GopherDir'].lower() not in ['','none']: gopherdir = self.info['GopherDir']
            spec = nextseq.info['SpecCode']
            ### ~ [1] ~ Special Organisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Organise']:
                orthdb = rje.baseFile(self.info['OrthDB'],strip_path=True)
                taxid = spec
                try:
                    if spec in self.dict['TaxaCode']: taxid = self.dict['TaxaCode'][spec]
                    elif rje.checkForFile(self.info['UniSpec']):
                        taxaline = rje.lineFromIndex(spec,self.info['UniSpec'],re_index='^(\S+)\s',sortunique=False,xreplace=True)
                        taxid = rje.matchExp(':(\d+):\s*$',taxaline)[0]
                        self.dict['TaxaCode'][spec] = taxid
                except: self.errorLog('Problem getting TaxaID for %s. Will use %s' % (spec,taxid))
                gopherdir = rje.makePath('%s%s/%s/' % (gopherdir,orthdb,taxid))
            if mkdir: rje.mkDir(self,gopherdir,log=True)
            return gopherdir
        except: self.errorLog('Problem with gopherDir setup (%s).' % (nextseq))       
#########################################################################################################################
    ### <4> ### Sticky Orthologous Groups Methods                                                                       #
#########################################################################################################################
    def stickySetup(self,seqfile,startfrom):  ### Sets up special StiGG (Sticky GOPHER Groups) options
        '''Sets up special StiGG (Sticky GOPHER Groups) options.'''
        if not self.opt['Sticky']: return
        self.list['StiGGMe'] = rje_seq.SeqInfoListFromFile(self,seqfile,'short')    #,startfrom)  # To StiGG
        if self.info['StiggID'].lower() == 'none': self.info['StiggID'] = ''
        while '' in self.list['StiGGMe']: self.list['StiGGMe'].remove('')
        self.list['StiGGing'] = self.list['StiGGMe'][0:]                                    # In line to be StiGGed
        if self.opt['Test']: self.list['StiGGing'] = []
        self.list['StiGGed'] = []                                                           # Already StiGGed
        self.list['StiGGPrime'] = []    # Not part of StiggMe but StiGGed directly by StiggMe sequences.
        self.list['StiGGHead'] = ['AccNum','Spec','StiGG','Orth','OrthSpec','OrthAcc']
        self.info['StiGGFile'] = '%s.stigg.tdt' % rje.baseFile(seqfile) 
#########################################################################################################################
    def stiGGMe(self,seq=None,acc=None):  ### Process StGG sticky list.
        '''Process StGG dictionary etc.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Sticky']: return
            if seq: acc = seq.info['AccNum']
            else: self.deBug('Fork.StiggMe - %s' % acc)
            seqfile = 'ORTH' + os.sep + '%s.orth.fas' % acc
            stickylist = rje_seq.SeqInfoListFromFile(self,seqfile,'short')  # To StiGG
            try: stickylist += string.split(open('ORTH' + os.sep + '%s.sticky.id' % acc,'r').read(),'\n')
            except: self.errorLog('No *.sticky.id file made for %s' % acc)
            self.deBug(stickylist)
            if seq: stigg = seq.shortName()
            else: stigg = stickylist[0]
            self.deBug('%s:%s' % (stigg,stigg in self.list['StiGGed']))
            #x#self.list['StiGGed'].append(stigg)
            ### ~ [1] Check StiGGMe List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            stiggprime = stigg in self.list['StiGGMe']
            stiggme = stigg in self.list['StiGGMe'] + self.list['StiGGPrime']
            for stiggorth in stickylist:
                if not stiggorth: continue
                if stiggprime and stiggorth not in self.list['StiGGMe'] + self.list['StiGGPrime']: self.list['StiGGPrime'].append(stiggorth)
                if stiggorth in self.list['StiGGMe'] + self.list['StiGGPrime']: stiggme = True
            if not stiggme:
                self.printLog('#SKIP','No %s sticky sequences of interest - skipping' % stigg)
                return      # Do not add extras if none of the sequences are of interest
            ### ~ [2] Extend StiGGing List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for stiggorth in stickylist:
                if not stiggorth: continue
                if stiggorth not in self.list['StiGGed'] + self.list['StiGGing']: self.list['StiGGing'].append(stiggorth)
        except: self.errorLog('Problem with GOPHER.stiGGMe()')
#########################################################################################################################
    def getSticky(self,stigg,gopherseq,gopherblast,accdict):
        if stigg in self.dict['Sticky']: return self.dict['Sticky'][stigg]
        stickylist = self.getOrth(stigg,gopherseq,gopherblast,accdict)
        acc = accdict[stigg]
        try: stickylist += string.split(open('ORTH' + os.sep + '%s.sticky.id' % acc,'r').read(),'\n')
        except: self.errorLog('No *.sticky.id file made for %s' % acc)
        while '' in stickylist: stickylist.remove('')
        self.dict['Sticky'][stigg] = stickylist
        return stickylist[0:]
#########################################################################################################################
    def getOrth(self,stigg,gopherseq,gopherblast,accdict):
        if stigg in self.dict['Orth']: return self.dict['Orth'][stigg]
        seq = gopherseq.seqFromFastaCmd(id=stigg,dbase=gopherblast.getStr('DBase'))
        if not seq:
            self.errorLog('Problem with %s - will try workaround.' % stigg)
            acc = string.split(stigg,'_')[-1]
        else: acc = seq.info['AccNum']
        accdict[stigg] = acc
        seqfile = 'ORTH' + os.sep + '%s.orth.fas' % acc
        if os.path.exists(seqfile): self.dict['Orth'][stigg] = rje_seq.SeqInfoListFromFile(self,seqfile,'short')
        else: self.dict['Orth'][stigg] = [stigg]
        return self.dict['Orth'][stigg][0:]
#########################################################################################################################
    def stiGG(self,gopherseq,gopherblast):    ### Generate and Process StGG dictionary from orthology files
        '''Generate and Process StGG dictionary from orthology files.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Sticky']: return
            self.dict['Sticky'] = {}    # Dictionary used to make StiGG sets
            self.dict['Orth'] = {}    # Dictionary used to make StiGG sets
            self.dict['StiGG'] = {}     # Dictionary used to make StiGG sets
            self.list['StiGG'] = []     # List of StiGG sets
            stiggprime = self.list['StiGGMe'] + self.list['StiGGPrime'] 
            accdict = {}
            rje.delimitedFileOutput(self,self.info['StiGGFile'],self.list['StiGGHead'],rje_backup=True)
            gopherseq.opt['Silent'] = True
            ### ~ [1] Generate Full StiGG Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(self.list['StiGGMe']); gx = 0; ox = 0
            for stigg in self.list['StiGGMe']:
                sx += 100.0
                if stigg in self.dict['StiGG']: continue
                self.dict['StiGG'][stigg] = self.getOrth(stigg,gopherseq,gopherblast,accdict)
                ## ~ [1a] Expand StiGG group ~ ##
                expand = True
                while expand:
                    expand = False
                    for orth in self.dict['StiGG'][stigg][0:]:
                        for stick in self.getSticky(orth,gopherseq,gopherblast,accdict):    # Consider adding this to orth?
                            self.progLog('\r#STIGG','Generating StiGG dictionary: %.2f%% (%s StiGG|%s|%d)     ' % (sx/stot,rje.integerString(gx),stigg,len(self.dict['StiGG'][stigg])))
                            for stickorth in self.getOrth(stick,gopherseq,gopherblast,accdict):
                                if stickorth in self.dict['StiGG'][stigg]:  # Shares an Orthologue with StiGG
                                    if stick not in self.dict['StiGG'][stigg]:
                                        if stick != stigg and stick in self.dict['StiGG']:
                                            self.errorLog('%s [<< %s] already in StiGG dictionary?!' % (stick,stigg))
                                            self.deBug('%s: %s' % (stigg,string.join(self.dict['StiGG'][stigg],', ')))
                                            self.deBug('%s: %s' % (stick,string.join(self.dict['StiGG'][stick],', ')))
                                            self.list['StiGG'].remove(self.dict['StiGG'][stick])
                                            ox -= len(self.dict['StiGG'][stick]); gx -= 1
                                            for x in self.dict['StiGG'][stick]: self.dict['StiGG'][stigg].append(x)
                                        else: self.dict['StiGG'][stigg].append(stick)
                                        expand = True
                ## ~ [1b] Reciprocate Stigg Group ~ ##
                sgrp = self.dict['StiGG'][stigg][0:]; ox += len(sgrp)
                for orth in sgrp: self.dict['StiGG'][orth] = sgrp
                self.list['StiGG'].append(sgrp)
                gx = len(self.list['StiGG'])
            self.printLog('\r#STIGG','Making sticky Orthologous groups complete: %s Seq -> %s Orth -> %s StiGG' % (rje.integerString(stot),rje.integerString(ox),rje.integerString(gx)))
            gx = 0
            for sgrp in self.list['StiGG']:
                agrp = []; gx += 1
                sid = '%s%s' % (self.info['StiggID'],rje.preZero(gx,stot))
                for orth in sgrp[0:]:
                    try: agrp.append(accdict[orth])
                    except: self.errorLog('Jeepers! %s not found in StiGG.accdict.' % stigg)
                for orth in sgrp[0:]:
                    acc = accdict[orth]
                    open('ORTH' + os.sep + '%s.stigg.id' % acc,'w').write(string.join(sgrp,'\n'))
                    open('ORTH' + os.sep + '%s.stigg.acc' % acc,'w').write(string.join(agrp,'\n'))
                    seqfile = rje.makePath('ORTH/%s.orth.fas' % acc,wholepath=True)
                    if os.path.exists(seqfile): gopherseq.loadSeqs(seqfile=seqfile)
                    else: continue
                    qry = gopherseq.querySeq(acc)
                    for seq in gopherseq.seq:
                        datadict = {'AccNum':acc,'StiGG':sid,'Spec':qry.info['SpecCode'],'Orth':seq.shortName(),
                                    'OrthSpec':seq.info['SpecCode'],'OrthAcc':seq.info['AccNum']}
                        rje.delimitedFileOutput(self,self.info['StiGGFile'],self.list['StiGGHead'],datadict=datadict)
            self.printLog('\r#STIGG','Output of sticky Orthologous groups complete: %s Seq -> %s Orth -> %s StiGG' % (rje.integerString(stot),rje.integerString(ox),rje.integerString(gx)))
            
        except: self.errorLog('Problem with GOPHER.stiGG()')
#########################################################################################################################
    ### <5> ### REST Setup Methods                                                                                      #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        ## ~ OrthBLAST outputs ~ ##
        qry = Query file
        blast = BLAST results file
        blastid = BLAST hit ID file

        ## ~ Orthologue Identification ~ ##
        orthfas = Unaligned fasta file of GOPHER orthologues
        orthid = ID list of GOPHER orthologues

        ## ~ Orthologue alignment ~ ##
        alnfas = Aligned fasta file of GOPHER orthologues

        ## ~ Orthologue tree ~ ##
        nwk = Tree of GOPHER orthologues in Newick format
        tree = Tree of GOPHER orthologues in plain text format
        png = Tree of GOPHER orthologues in PNG format

        ## ~ Paralogues ~ ##
        parafas = Unaligned fasta file of GOPHER paralogues
        paraid = ID list of GOPHER paralogues
        paralnfas = Aligned fasta file of GOPHER paralogues

        Note that individual outputs can be requested for single proteins using the special REST call:

            http://rest.slimsuite.unsw.edu.au/gopher&acc=X&spcode=X

        If this protein has already been processed by GOPHER, the relevant output will be returned directly as plain
        text. If not, a jobid will be returned, which will have the desired output once run.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Output'] = {}
            # Only generate tabs for output that is made.
            #for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['qry','blast','blastid','orthfas','orthid','alnfas','nwk','tree','png','parafas','paraid','paralnfas']
#########################################################################################################################
    def restFileReturn(self,acc,spec):   ### Returns contents of desired output file, else None
        '''Returns contents of desired output file, else None.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfmt = self.getStrLC('Rest')
            if not outfmt: outfmt = 'alnfas'
            gopherdir = self.getStr('GopherDir')
            orthdb = rje.baseFile(self.getStr('OrthDB'),strip_path=True)
            taxid = spec
            try:
                if spec in self.dict['TaxaCode']: taxid = self.dict['TaxaCode'][spec]
                elif rje.checkForFile(self.getStr('UniSpec')):
                    taxaline = rje.lineFromIndex(spec,self.getStr('UniSpec'),re_index='^(\S+)\s',sortunique=False,xreplace=True)
                    taxid = rje.matchExp(':(\d+):\s*$',taxaline)[0]
                    self.dict['TaxaCode'][spec] = taxid
            except: self.errorLog('Problem getting TaxaID for %s. Will use %s' % (spec,taxid))
            gopherdir = rje.makePath('%s%s/%s/' % (gopherdir,orthdb,taxid))
            seq = rje_sequence.Sequence(self.log,self.cmd_list)
            seq.info['AccNum'] = acc
            gopherfork = GopherFork(self.log,self.cmd_list)
            gopherfork.obj['Sequence'] = seq
            gopherfork.info['OutPath'] = gopherdir
            ### ~ [1] ~ Outfmt output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfmt == 'rlc':
                rlcfile = gopherfork.gFile('ALN','orthaln.rlc.txt')
                self.printLog('#GRLC','%s: %s' % (rlcfile,os.path.exists(rlcfile)))
                if os.path.exists(rlcfile): return open(rlcfile,'r').read()
                else: return None
            self.printLog('#GDIR',gopherdir)
            gopherfork.restOutputs()
            if outfmt in gopherfork.dict['Output']: return open(gopherfork.dict['Output'][outfmt],'r').read()
            return None
        except: return 'ERROR: %s' % self.errorLog('restFileReturn(%s) failed' % outfmt)
#########################################################################################################################
### End of Section II: Gopher Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: GopherFork Class:                                                                                      #
#########################################################################################################################
class GopherFork(rje.RJE_Object):     
    '''
    GopherFork Class. Author: Rich Edwards (2005).

    This class is designed to handle each forking of the main Gopher program. This will carry out the various stages
    associated with each sequence read by the main gopher method.

    Info:str
    - AlnProg = Choice of alignment program to use (clustalw/clustalo/muscle/mafft/fsa) [clustalo]
    - Name = Sequence ShortName
    - OrthDB = Fasta file with pool of sequences for orthology discovery [orthdb.fas]
    - OrthID = File identifier (Lower case) for orthology files [orth]
    - SimFocus = Style of MinSim used [query]
        - query = %query must > minsim (Best if query is ultimate focus and maximises closeness of returned orthologues)
        - hit = %hit must > minsim (Best if lots of sequence fragments are in searchdb and should be retained)
        - either = %query > minsim OR %hit > minsim (Best if both above conditions are true)
        - both = %query > minsim AND %hit > minsim (Gets most similar sequences in terms of length)
    - GABLAMO Key = GABLAMO measure to use for similarity measures [Sim]
        - ID = %Identity (from BLAST)
        - Sim = %Similarity (from BLAST)
        - Len = %Coverage (from BLAST)
    - Mode = Gopher mode being run for fork [None]
    - OutPath = Parent path for output files, incorporating Organise and UniSpec settings [./]
    
    Opt:boolean
    - CompFilter = Whether to use complexity filter and composition statistics for *initial* BLAST [True]
    - DropOut = Whether to "drop out" at earlier phases, or continue with single sequence [False]
    - Force = Whether to force execution even if results are new enough [False]
    - FullForce = Whether to force current and previous execution even if results are new enough [False]
    - OldBLAST = Whether to use old BLAST in place of BLAST+ [False]
    - NoExec = Whether to simply report what *would* be executed [False]
    - Organise = Output files according to orthdb an species (code or TaxaID - need conversion) [True]
    - PostDup = Whether to align only post-duplication sequences [True]
    - IgnoreDate = Ignores the age of files and only replaces if missing [False]
    - SaveSpace = Save space by deleting intermediate blast files during orthfas [True]
    - Paralign = Whether to produce paralogue alignments (>minsim) in PARALN/ [True]
    - ParaSplice = Whether to allow splice variants in paralogue alignments (where identified) [False]
    - ParaFam = Whether to paralogue paired subfamily alignments (>minsim) (assuming run to orthfas+) [False]
    - GopherFam = Whether to combined paralogous gopher orthologues into protein families (>minsim) (assuming run to orthfas+) [False]
    - FullRBH = Whether RBH method should run BLAST searches for all potential RBH orthologues [False]
    - Reciprocal = Use Reciprocal Best Hit method instead of standard GOPHER approach [False]
    - Sticky = Switch on "Sticky Orthologous Group generation" [False]

    Stat:numeric
    - MinSim = Min %Sim for orth vs qry or qry vs orth (i.e. shorter versus longer)
    - MaxPara = Maximum number of paralogues to consider (large gene families can cause problems) [50]

    List:list

    Dict:dictionary    
    - Output = REST server output. Only works for single sequence runs.

    Obj:RJE_Objects
    - Parent:gopher.Gopher object controlling this run.
    - Sequence:rje_seq.Sequence = 'Query' Sequence
    - BLAST:rje_blast.BLASTRun = Blast Run. (If None, invoke 'gopher' mode)

    Run modes:
    orthblast   :   Run to blasting versus orthdb (Stage 1).
    orthfas     :   Run to output of orthologues (Stage 2). 
    orthalign   :   Run to alignment of orthologues (Stage 3). 
    orthtree    :   Run to generation of trees (Stage 4). 
    orthanc     :   Run to generation of AncSeqs (Stage 5, GASP).
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AlnProg','Name','OrthID','OrthDB','GABLAMO Key','SimFocus','Mode','OutPath']
        self.statlist = ['MinSim','MaxPara']
        self.optlist = ['CompFilter','PostDup','Force','FullForce','NoExec','Repair','IgnoreDate','SaveSpace','FullRBH',
                        'Reciprocal','Paralign','ParaSplice','ParaFam','GopherFam','DNA','Sticky','DropOut','Organise']
        self.listlist = []
        self.dictlist = ['Output']
        self.objlist = ['Sequence','BLAST']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'GABLAMO Key':'Sim','SimFocus':'query','Mode':None,'OutPath':rje.makePath('./'),'AlnProg':'clustalo','OrthID':'orth'})
        self.setStat({'MinSim':40.0,'MaxPara':50})
        self.setOpt({'Repair':True,'SaveSpace':True,'CompFilter':True,'Organise':True,'OldBLAST':False})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdRead(cmd,type='info',att='GABLAMO Key',arg='gablamo')
                self._cmdReadList(cmd.lower(),'info',['SimFocus','AlnProg','OrthID'])
                self._cmdReadList(cmd,'file',['OrthDB'])
                self._cmdReadList(cmd,'stat',['MinSim'])
                self._cmdReadList(cmd,'int',['MaxPara'])
                self._cmdReadList(cmd,'opt',['PostDup','Force','FullForce','NoExec','Repair','IgnoreDate','SaveSpace',
                                             'Paralign','ParaFam','GopherFam','DNA','Sticky','DropOut','Reciprocal',
                                             'CompFilter','Organise','FullRBH','OldBLAST'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        if self.stat['MinSim'] < 1: self.stat['MinSim'] *= 10
        if self.opt['FullForce']: self.opt['Force'] = True
        if self.getBool('FullRBH'): self.setBool({'Reciprocal':True})
#########################################################################################################################
    def gFile(self,subdir,extension,acc=None):   ### Returns full path including OUTPATH and AccNum                          #V3.0
        '''Returns full path including OUTPATH and AccNum: self.info['OutPath']/subdir/AccNum.extension.'''
        if not acc: acc = self.obj['Sequence'].info['AccNum']
        if 'ALN' in subdir and self.opt['Organise']: subdir = self.info['AlnProg'] + subdir
        #subdir = string.replace(rje.makePath(subdir),'ORTH',self.info['OrthID'].upper())
        subdir = rje.makePath(subdir)
        extension = string.replace(extension,'orth',self.info['OrthID'].lower())
        return rje.makePath('%s%s%s.%s' % (self.info['OutPath'],subdir,acc,extension),wholepath=True)
#########################################################################################################################
    def restOutputs(self):   ### Updates REST outputs that have been generated.
        '''Updates REST outputs that have been generated.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfmts = {'qry':self.gFile('BLAST','qry'),       # Query file
                        'blast':self.gFile('BLAST','blast'),     # BLAST results file
                        'blastid':self.gFile('BLAST','blast.id'),  # BLAST hit ID file
                        'orthfas':self.gFile('ORTH','orth.fas'),   # Generated now
                        'orthid':self.gFile('ORTH','orth.id'),    # Generated now
                        'parafas':self.gFile('PARA','para.fas'),   # Define here as will check to delete later
                        'paraid':self.gFile('PARA','closepara.id'),   # Used by self._parAlign()
                        'alnfas':self.gFile('ALN','orthaln.fas'),
                        'nwk':self.gFile('TREE','orth.nwk'),
                        'tree':self.gFile('TREE','orth.tree.txt'),
                        'png':self.gFile('TREE','orth.png'),
                        'paralnfas':self.gFile('PARALN','paraln.fas')} # Will generate this now
            for outfmt in outfmts:
                #print outfmt, outfmts[outfmt], rje.exists(outfmts[outfmt])
                if rje.exists(outfmts[outfmt]): self.dict['Output'][outfmt] = outfmts[outfmt]
        except: return 'ERROR: %s' % self.errorLog('restOutputs() failed')
#########################################################################################################################
    ### <2> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def run(self,mode):      ### Main GopherFork Run method                                                         #V3.0
        '''Main GopherFork Run method. Calls _runStage(mode) as appropriate.'''
        ### ~ [1] ~ Perform usual run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.info['Mode'] = mode
        if mode in gophermodes:
            self._runStage(mode)
            if gophermodes.index(mode) >= gophermodes.index('orthfas'):     ### Run far enough
                self._parAlign()
                self._paraFam()
                self._gopherFam()
        else: self.errorLog('GopherFork(%s) Problem: %s mode not recognised' % (self.info['Name'],mode))
        ### ~ [2] ~ Special SoapLab adjustment and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['SoapLab']:
            ## ~ [2a] ~ SoapLab output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #orthfas = 'ORTH' + os.sep + '%s.orth.fas' % self.obj['Sequence'].info['AccNum']
            orthfas = self.gFile('ORTH','orth.fas')
            if os.path.exists(orthfas): open('gopher.fas','w').writelines(open(orthfas,'r').readlines())
            else: open('gopher.fas','w').write('>%s\n%s' % (self.obj['Sequence'].info['Name'],self.obj['Sequence'].info['Sequence']))
            ## ~ [2b] ~ Clean up unwanted files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gpath in ['BLAST','ORTH','PARA']:
                for gfile in glob.glob(self.gFile(gpath,'*')): os.unlink(gfile)
                cleandir = rje.makePath('%s%s' % (self.info['OutPath'],gpath))
                if os.path.exists(cleandir) and not glob.glob(self.gFile(gpath,'*')): os.rmdir(cleandir)
            for gfile in ['i_seqin','i_orthdb']: rje_blast.cleanupDB(callobj=self,dbfile=gfile,deletesource=False)   
#########################################################################################################################
    def _runStage(self,stage):  ### Generic method for checking a stage has been completed OK                       #V3.0
        '''
        Generic method for checking a stage has been completed OK.
        >> stage:str = current stage
        << returns False if ending here due to failure etc. or True if stage completed and ready for next stage.
        '''
        ### ~ [1] ~ Setup assessment lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #!# Possibly make an rje_pipeline module to handle this kind of thing? #!#
        try:## ~ [1a] ~ Generic variable setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = self.obj['BLAST']
            ## ~ [1b] ~ Setup Stages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if stage == 'dbases':  return False   # Handled in main thread - only here as it is a pre-requisite for other stages
            stagelist = ['dbases'] + gophermodes 
            stageneed = {}
            # Master Database
            stageneed['dbases'] = [blast.getStr('DBase')]
            # OrthBlast
            stageneed['orthblast'] = [self.gFile('BLAST','blast.id')]   # This is the only file actually used next
            # OrthFas
            stageneed['orthfas'] = [self.gFile('ORTH','orth.fas')]
            if self.opt['Sticky']: stageneed['orthfas'] += [self.gFile('ORTH','sticky.id')]
            # OrthAln
            stageneed['orthalign'] = [self.gFile('ALN','orthaln.fas')]  #!# Got up to here. Replace all paths with self.gFile() #!#
            # OrthTree
            stageneed['orthtree'] = [self.gFile('TREE','orth.nwk')]
            ## ~ [1c] ~ Lists for self._needStage() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            needfiles = stageneed[stage]
            prevstage = stagelist[stagelist.index(stage)-1]
            prevfile = stageneed[stagelist[stagelist.index(stage)-2]][-1]
        except:
            self.log.errorLog('Major Problem with GopherFork(%s)._runStage(%s) Setup.' % (self.info['Name'],stage))
            sys.exit(1)
        ### ~ [2] ~ Checks for run needs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:## ~ [2a] ~ If Repairmode off, only check current level: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.opt['Repair'] and not self._needStage(stage,needfiles,prevfile): return False
            ## ~ [2b] ~ Check previous stage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self._runStage(prevstage): pass      # Earlier stage ran, so this needs to as well
            elif prevstage != 'dbases' and self._needStage(prevstage,stageneed[prevstage],prevfile):
                # Running the previous stage failed and yet in still needs to be run - exit!
                self.printLog('#QRY','Sequence %s has dropped out of GOPHER at %s.' % (self.info['Name'],prevstage))
                return False
            elif not self._needStage(stage,needfiles,prevfile): return False   # Check own status: No need 
            ## ~ [2c] ~ NoExec - report intent only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['NoExec']: return self.printLog('#RUN','Run GopherFork(%s) %s.' % (self.info['Name'],stage))
        ### ~ [3] ~ Execute Specific Stage Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if stage == 'orthblast': return self._orthBlast(self.opt['CompFilter'])
            elif stage == 'orthfas': return self._orthFas()
            elif stage == 'orthalign': return self._orthAln()
            elif stage == 'orthtree': return self._orthTree()
            return True
        except:
            self.errorLog('Major Problem with GopherFork(%s)._runStage(%s).' % (self.info['Name'],stage))
            sys.exit(1)
#########################################################################################################################
    def _needStage(self,stage,needfiles,prevfile):  ### Generic method for checking a stage has been completed OK   #V3.0
        '''
        Generic method for checking a stage has been completed OK.
        >> stage:str = current stage
        >> needfiles:list = list of filenames produced by stage in order created
        >> prevfile:str = filename for last (i.e. stats?) output of previous stage
        '''
        try:### ~ [1] ~ Check for missing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for file in needfiles:
                if not os.path.exists(file): return True                    # File missing
            ### ~ [2] ~ Repairmode and old files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Repair'] and not self.opt['IgnoreDate']:
                for file in needfiles:
                    if rje.isYounger(file,prevfile) != file: return True    # File not newer than prevfile
            ### ~ [3] ~ Forced running ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['FullForce'] or (self.opt['Force'] and self.info['Mode'] == stage): return True          
            return False
        except: self.errorLog('Major Problem with GopherFork(%s)._needStage(%s).' % (self.info['Name'],stage)); raise
#########################################################################################################################
    def _orthBlast(self,compf=True):   ### Blasts Query against OrthDB and populates BLAST object                   #V3.3
        '''
        Blasts Query against OrthDB and populates BLAST object.
        1. BLAST against orthdb.
        2. Read BLAST and save sequence IDs as BLAST/AccNum.blast.id
        << Returns True if Query OK and BLAST OK. Otherwise, returns False.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('OldBLAST'): return self._orthBlastLegacy(compf)
            ## ~ [1a] Query Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = self.obj['Sequence']  # Query Sequence Object
            acc = qry.info['AccNum']
            spec = qry.info['SpecCode']
            ## ~ [1b] Files and Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            qfile = self.gFile('BLAST','qry')       # Query file
            bfile = self.gFile('BLAST','blast')     # BLAST results file
            ifile = self.gFile('BLAST','blast.id')  # BLAST hit ID file
            rje.mkDir(self,bfile,log=True)
            ## ~ [1c] BLAST Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = self.obj['BLAST']
            blast.setStr({'Name':bfile,'InFile':qfile})
            blast.setBool({'Complexity Filter':compf,'Composition Statistics': compf})
            ## ~ [1d] Check for files etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            needtoblast = self.opt['Force']
            if self.opt['Force'] or not os.path.exists(qfile) or (rje.isYounger(qfile,blast.getStr('DBase')) != qfile and not self.opt['IgnoreDate']):   # Either missing, else database younger than query
                open(qfile,'w').write('>%s\n%s\n' % (qry.info['Name'], qry.info['Sequence']))
            ### ~ [2] BLAST against orthdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.blast(wait=True,use_existing=os.path.exists(bfile) and not self.force())    # Will check bfile and qfile etc.
            ### ~ [3] Check BLAST results for query sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.readBLAST(clear=True,gablam=False)
            if not blast.hitNum():
                if compf:
                    self.printLog('#COMPF','%s yielded no BLAST hits. Low complexity sequence? Will try with -F F' % self.info['Name'])
                    if os.path.exists(bfile): os.unlink(bfile)
                    if os.path.exists(ifile): os.unlink(ifile)
                    return self._orthBlast(compf=False)
                else:
                    self.printLog('#QRY','%s yielded no BLAST hits. Short sequence?' % self.info['Name'])
                    if not self.opt['DropOut']: open(ifile,'w').write(qry.shortName()); return True
                    return False
            ## ~ [3a] Check for Self Hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qryhit = qry.shortName() in blast.db('Hit').index('Hit')
            if not qryhit:
                self.printLog('#QRY','WARNING: Query protein %s not (hit) in search database!' % (self.info['Name']))
            ### ~ [4] Save hit IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.saveHitIDs(outfile=ifile)
            ### ~ [5] Evaluate need for next stage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if blast.hitNum() > 1 or (not qryhit and blast.hitNum() > 0): return True           # May be orthologues
            elif not self.opt['DropOut']: open(ifile,'w').write(qry.shortName()); return True   # Query only hits itself
            else: return False      # Query only hits itself               
        except: self.errorLog('Major Problem with GopherFork(%s).blast()' % self.info['Name']); raise
#########################################################################################################################
    def _orthBlastLegacy(self,compf=True):   ### Blasts Query against OrthDB and populates BLAST object                   #V3.0
        '''
        Blasts Query against OrthDB and populates BLAST object.
        1. BLAST against orthdb.
        2. Read BLAST and save sequence IDs as BLAST/AccNum.blast.id
        << Returns True if Query OK and BLAST OK. Otherwise, returns False.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Query Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = self.obj['Sequence']  # Query Sequence Object
            acc = qry.info['AccNum']
            spec = qry.info['SpecCode']
            ## ~ [1b] Files and Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            qfile = self.gFile('BLAST','qry')
            bfile = self.gFile('BLAST','blast')
            ifile = self.gFile('BLAST','blast.id')
            rje.mkDir(self,bfile,log=True)
            ## ~ [1c] BLAST Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = self.obj['BLAST']
            blast.setStr({'Name':bfile,'InFile':qfile})
            blast.setBool({'Complexity Filter':compf,'Composition Statistics': compf})
            ## ~ [1d] Check for files etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            needtoblast = self.opt['Force']
            if self.opt['Force'] or not os.path.exists(qfile) or (rje.isYounger(qfile,blast.getStr('DBase')) != qfile and not self.opt['IgnoreDate']):   # Either missing, else database younger than query
                open(qfile,'w').write('>%s\n%s\n' % (qry.info['Name'], qry.info['Sequence']))
            if not os.path.exists(bfile) or (not self.opt['IgnoreDate'] and rje.isYounger(bfile,blast.getStr('DBase')) != bfile):   # Either missing, else database younger than results
                needtoblast = True
            elif not os.path.exists(qfile) or (not self.opt['IgnoreDate'] and rje.isYounger(bfile,qfile) != bfile): # Either missing, else query younger than results
                needtoblast = True
            else: needtoblast = not blast.checkBLAST()
            ### ~ [2] BLAST against orthdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if needtoblast: blast.blast(wait=True)           
            ### ~ [3] Check BLAST results for query sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast.readBLAST(clear=True)
            try:
                search = blast.search[0]
                if not search.hit: raise ValueError    # No BLAST Hits - probably a short sequence - chuck out at gopher stage!
            except:
                if compf:
                    self.log.printLog('#COMPF','%s yielded no BLAST hits. Low complexity sequence? Will try with -F F' % self.info['Name'])
                    if os.path.exists(bfile): os.unlink(bfile)
                    if os.path.exists(ifile): os.unlink(ifile)
                    return self._orthBlast(compf=False)
                else:
                    self.log.printLog('#QRY','%s yielded no BLAST hits. Short sequence?' % self.info['Name'])
                    if not self.opt['DropOut']: open(ifile,'w').write(qry.shortName()); return True
                    return False
            ## ~ [3a] Check for Self Hit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qryhit = None
            for hit in search.hit:
                if hit.info['Name'] == qry.shortName():        
                    qryhit = hit
                    break
            if not qryhit: self.log.printLog('#QRY','WARNING: Query protein %s not (hit) in search database!' % (self.info['Name']))
            ### ~ [4] Save hit IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            search.saveHitIDs(outfile=ifile)
            ### ~ [5] Evaluate need for next stage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if search.hitNum() > 1 or (not qryhit and search.hitNum() > 0): return True         # May be orthologues
            elif not self.opt['DropOut']: open(ifile,'w').write(qry.shortName()); return True   # Query only hits itself
            else: return False      # Query only hits itself               
        except: self.errorLog('Major Problem with GopherFork(%s).blast()' % self.info['Name']); raise
#########################################################################################################################
    def _orthFas(self):     ### Makes an unaligned Fasta file of potential orthologues.                             #V3.3
        '''
        Makes an unaligned Fasta file of potential orthologues.
        > Save reduced sequences as ORTH/AccNum.orth.fas
        '''
        try:### ~ [1] Setup Variables needed for all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Setup'
            if self.getBool('OldBLAST'): return self._orthFasLegacy()
            ## ~ [1a] Setup variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = self.obj['Sequence']
            spec = qry.info['SpecCode']
            stickhits = []
            ssx = 0     # Counter for save space deletions
            ## ~ [1b] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ifile = self.gFile('BLAST','blast.id')  # Generated by _orthBlast()
            qfile = self.gFile('BLAST','qry')       # Generated by _orthBlast()
            bfile = self.gFile('ORTH','orth.blast') # Will generate this now
            ffile = self.gFile('ORTH','full.fas')   # Generated now but cleared by savespace
            ofile = self.gFile('ORTH','orth.fas')   # Generated now
            afile = self.gFile('ORTH','orth.id')    # Generated now
            mfile = self.gFile('ORTH','minsim.fas') # Clean this file up
            pfile = self.gFile('PARA','para.fas')   # Define here as will check to delete later
            pbfile = self.gFile('PARA','para.blast')# Clean this file up
            rje.mkDir(self,ofile,log=True)
            rje.mkDir(self,pfile,log=True)
            ## ~ [1c] Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = self.obj['BLAST']
            blast.setStr({'Name':bfile})
            blast.clear()
            ## ~ [1d] Setup GABLAMO Matrix and stats. Used for comparing Queries and Hits ~~~~~~~~~ ##
            gablamo_matrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            gablamo_matrix.opt['Symmetric'] = False
            if self.info['GABLAMO Key'] not in ['Sim','ID','Len','Score']:
                self.errorLog('GABLAMO Key (%s) not Sim/ID/Len. Changed to Sim.' % self.info['GABLAMO Key'],printerror=False,quitchoice=True)
                self.info['GABLAMO Key'] = 'Sim'
            if self.info['GABLAMO Key'] != 'Score': gkey = 'GABLAMO %s' % self.info['GABLAMO Key']
            else: gkey = self.info['GABLAMO Key']

            ### ~ [2] Read in original BLAST Hits and filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Read BLAST Hits'
            ## ~ [2a] Check for ID list file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(ifile):
                if self._orthBlast(self.opt['CompFilter']): return self._orthFas()    # Try _orthBlast again and then try running self
                elif not self.opt['DropOut']:
                    open(ofile,'w').write('>%s\n%s\n' % (qry.info['Name'],qry.info['Sequence']))
                    return True
                else:                   # BLAST Unsuccessful
                    self.opt['OrthAlign'] = False   
                    return False
            ## ~ [2b] Read in proteins hit from original BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=T','unkspec=F','specnr=F','gnspacc=T'])
            idlist = string.split(open(ifile,'r').read())
            idseq.seq = idseq.seqFromBlastDBCmd(id=idlist,dbase=blast.getStr('DBase'),expect=False) #Qry might not be in OrthDB
            ## ~ [2c] Filter sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if idseq.opt['AutoFilter']: idseq.autoFilter(cmd_list=self.cmd_list+['filterout=None','seqout=None'])
            ## ~ [2d] Check Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if idseq.querySeq(query=qry.shortName()):   # Remove and re-add for safety/consistency
                idseq.removeSeq(text='Query sequence in BLAST database',seq=idseq.obj['QuerySeq'])
                self.printLog('#QRY','Query sequence %s replaced with Input sequence.' % qry.shortName())
            else: self.printLog('#QRY','Query sequence %s added to BLAST Hits' % qry.shortName())
            idseq._addSeq(qry.info['Name'],qry.info['Sequence'])    ### Adds a new Sequence Object to list
            idseq.seq = idseq.seq[-1:] + idseq.seq[:-1]
            idseq.querySeq(query=qry.shortName())       # Make Query 
            qseq = idseq.obj['QuerySeq']
            ## ~ [2e] Save filtered sequences for second BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idseq.saveFasta(seqfile=ffile)
            if not rje.exists(qfile): open(qfile,'w').write('>%s\n%s\n' % (qseq.info['Name'],qseq.info['Sequence']))

            ### ~ [3] New BLAST to get full GABLAM stats for GOPHER calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'New BLAST'
            ## ~ [3a] Make new BLAST database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newblast = rje_blast.blastObj(log=self.log,cmd_list=self.cmd_list,type='New')
            newblast.setStr(blast.getStr())
            newblast.setBool(blast.getBool())
            newblast.setInt(blast.getInt())
            newblast.setNum(blast.getNum())
            newblast.setStr({'InFile':qfile,'DBase':ffile})     # BLAST Query against first round hits
            newblast.setInt({'OneLine':idseq.seqNum(),'HitAln':idseq.seqNum()})
            newblast.setNum({'E-Value':10})
            newblast.setBool({'Complexity Filter':False,'Composition Statistics':False})
            newblast.formatDB(protein=not self.opt['DNA'])
            newblast.blast(cleandb=True)    # BLAST vs homologues and remove database files
            ## ~ [3b] Read BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newblast.readBLAST(clear=True,gablam=True)
            if self.opt['SaveSpace']:
                for ssfile in [ffile,blast.getStr('Name'),newblast.getStr('Name')]:
                    if os.path.exists(ssfile): os.unlink(ssfile); ssx += 1
            if not newblast.hitNum():   # No hits, including self
                self.errorLog('%s has no BLAST hits (including self!).' % self.info['Name'],printerror=False)
                self.opt['OrthAlign'] = False   # Nothing to align!
                if self.opt['SaveSpace']: self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
                return False
            ### ~ [3c] Map sequences onto hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'MapSeq'
            hitdict = {}        # This will start as a dictionary of {search:hitdict} and have paralogue searches added
            search = qry.shortName()
            hitdict[search] = newblast.hitToSeq(idseq)            # Dictionary of {Hit:Sequence Object}
            if not qseq or qseq not in hitdict[search].values():    # No self hit
                self.errorLog('Cannot find self hit for %s in reduced BLAST hits.' % self.info['Name'],printerror=False)
                if self.opt['SaveSpace']: self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
                if not self.opt['DropOut']:
                    open(ofile,'w').write('>%s\n%s\n' % (qry.info['Name'],qry.info['Sequence']))
                    return True
                self.opt['OrthAlign'] = False   # Nothing to align!
                return False
            
            ### ~ [4] Similarity Filter, Best Species Hit & Paralogues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Similarity Filter'
            ## Best spec dictionaries store the Sim scores (and Reciprocal score) for the best hit of that sequence:
            ## - best_qspec = {seq:best hit of query species}
            ## - best_hspec = {hspec:best hit for query of hspec sequences}
            best_qspec = {}     # Dictionary of {seq:{'Sim':X,'Rec':X,'Seq':seq}}
            best_hspec = {}     # Dictionary of {hspec:{'Sim':X,'Rec':X,'Seq':seq}}
            paraseq = []        # Paralogous sequences
            baseacc = []        # Base AccNums for spotting splice variants
            if not self.opt['ParaSplice']:  # Add query base accnum to baseacc to filter splice variants of query
                if rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum']): baseacc.append(rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum'])[0])
                else: baseacc.append(qseq.info['AccNum'])
            simlost = 0
            qlen = qseq.aaLen()
            ## ~ [4a] Go through each hit in turn and establish similarity and best hits ~~~~~~~~~~ ##
            for hentry in newblast.db('Hit').entries():
                ## ~~~ GABLAMO ~~~ ##
                hit = hentry['Hit']
                search = hentry['Query']
                hseq = hitdict[search][hit]
                if hseq == qseq:
                    if gkey == 'Score': gablamo_matrix.addDis(qseq,qseq,hentry['BitScore'])
                    else: gablamo_matrix.addDis(qseq,qseq,100.0)
                    continue
                if not hseq:
                    self.errorLog('Unable to match idseq sequence for %s hit %s' % (self.info['Name'],hit),False,False)
                    self.printLog('#CHECK','Check %s files for problem hit %s' % (self.info['Name'],hit),screen=False)
                    newblast.db('Hit').dropEntry(hentry)
                gdict = {'Hit':newblast.db('GABLAM').data(makekey={'Hit':hit,'Query':search,'QryHit':'Hit'}),
                         'Query':newblast.db('GABLAM').data(makekey={'Hit':hit,'Query':search,'QryHit':'Query'})}
                hspec = hseq.info['SpecCode']
                if gkey == 'Score': qvh = hentry['BitScore']
                else: qvh = float(100 * gdict['Query'][gkey]) / float(qlen)
                hvq = 0
                if gkey == 'Score': hvq = hentry['BitScore']
                elif hentry['Length'] > 0: hvq = float(100 * gdict['Hit'][gkey]) / float(hentry['Length'])
                else: self.errorLog('%s Hit %s has zero length!' % (self.info['Name'],hit),False,False)
                ## ~~~ Paralogue ~~~ ##
                if hspec == spec:   # Paralogue or splice variant
                    if not self.opt['ParaSplice']:
                        pacc = hseq.info['AccNum']
                        if rje.matchExp('^(\S+)\-(\d+)$',pacc):
                            pacc = rje.matchExp('^(\S+)\-(\d+)$',pacc)[0]
                        if pacc in baseacc:     # Splice variant => Toss Hit
                            newblast.db('Hit').dropEntry(hentry)
                            idseq.removeSeq(text='Splice variant (%s)' % pacc,seq=hseq)
                            continue                   
                    paraseq.append(hseq)
                ## ~~~ MinSim Filter: Only MinSim potential orthlogues so that all paralogues are used for comparison ~~~ ##
                badsim = False
                if gkey != 'Score': 
                    if self.info['SimFocus'] in ['query','both'] and qvh < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] in ['hit','both'] and hvq < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] == 'either' and qvh < self.stat['MinSim'] and qvh < self.stat['MinSim']: badsim = True
                if badsim and hspec != spec:
                    simlost += 1
                    newblast.db('Hit').dropEntry(hentry)
                    while hseq in idseq.seq: idseq.seq.remove(hseq)     #?# Can a sequence be in here twice?! #?#
                    continue
                ## ~~~ Update best_hit ~~~ ##
                gablamo_matrix.addDis(qseq,hseq,qvh)
                gablamo_matrix.addDis(hseq,qseq,hvq)
                best_qspec[hseq] = {'Sim':hvq,'Rec':qvh,'Seq':qseq}     # At this point, the Query is the best qspec sequence for hseq
                if best_hspec.has_key(hspec) and (qvh < best_hspec[hspec]['Sim'] or (qvh == best_hspec[hspec]['Sim'] and hvq <= best_hspec[hspec]['Rec'] and best_hspec[hspec]['Seq'])):
                    continue    # Equal or better already exists
                # Best so far for that species...
                best_hspec[hspec] = {'Sim':qvh,'Rec':hvq,'Seq':hseq}
            if gkey == 'Score': self.printLog('#SIM','%s: %d hits. Using BitScore, so no Min. Similarity threshold.' % (self.info['Name'],newblast.hitNum(search)))
            else: self.printLog('#SIM','%s: %d of %d hits failed to meet Min. Similarity threshold (%s %.2f%%) - removed.' % (self.info['Name'],simlost,(newblast.hitNum(search)+simlost),gkey,self.stat['MinSim']))
            ## ~ [4b] Remove all but best hits of each species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orthseq = []    # List of sequences to compare paralogues to
            for hit in newblast.queryHits(search):
                hseq = hitdict[search][hit]
                hspec = hseq.info['SpecCode']
                if hspec == spec: continue  # Paralogue 
                if best_hspec.has_key(hspec) and best_hspec[hspec]['Seq'] == hseq:    # Best of species
                    orthseq.append(hseq)
                else:   # Bye bye...!
                    stickhits.append(hseq.shortName())
                    try:
                        hentry = newblast.hitData(hit,search)
                        newblast.db('Hit').dropEntry(hentry)
                        idseq.seq.remove(hseq)
                    except:
                        self.log.errorLog('Problem during remove (%s/%s).' % (hit,hseq.shortName()),printerror=False)
            self.printLog('#SIM','%s: %d best hits for species.' % (self.info['Name'],len(orthseq)))
            ## ~ [4c] Special Reciprocal Best Hit Method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry_paraseq = paraseq[0:]
            if self.getBool('Reciprocal') and self.getBool('FullRBH'):
                px = len(paraseq)
                paradict = {}
                for pseq in paraseq: paradict[pseq.shortName()] = pseq
                if len(paradict) != len(paraseq): raise ValueError
                self.printLog('#RBH','Performing GopherBLAST searches for %d putative orthologues' % len(orthseq))
                for orth in orthseq:
                    orthdir = self.obj['Parent'].gopherDir(orth)
                    orthfork = GopherFork(log=self.log,cmd_list=self.cmd_list+['i=-1','errorlog=None'])
                    orthfork.setInfo({'Name':orth.shortName(),'OutPath':orthdir})
                    orthfork.obj['Sequence'] = orth
                    orthfork.obj['BLAST'] = self.obj['BLAST']
                    orthfork.obj['Parent'] = self.obj['Parent']
                    orthfork.log = self.log
                    orthfork.run('orthblast')
                    orthidseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=T','unkspec=F','specnr=F','gnspacc=T'])
                    orthidseq.seq = []
                    for id in open(orthfork.gFile('BLAST','blast.id'),'r').readlines():
                        id = rje.chomp(id)
                        if id: orthidseq.seqFromBlastDBCmd(id=id,dbase=blast.getStr('DBase'))
                    for oseq in orthidseq.seq:
                        if oseq.shortName() in paradict: continue
                        if qseq.sameSpec(oseq) and oseq not in paraseq: paraseq.append(oseq); paradict[oseq.shortName()] = oseq
                self.printLog('#RBH','%d additional paralogues added.' % (len(paraseq)-px))
                if len(paradict) != len(paraseq): raise ValueError

            ### ~ [5] ParaSeq BLAST to establish Orthologue-Paralogue relationships ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# When updating to make full reciprocal BLAST version, this will need to call full GOPHER BLAST for minsim sequences #!#
            _stage = 'ParaSeq BLAST'
            self.printLog('#PARA','%s: %d paralogues to consider.' % (self.info['Name'],len(paraseq)))
            if len(paraseq) > self.stat['MaxPara']:
                self.printLog('#PARA','%s: %d paralogues reduced to %d (MaxPara).' % (self.info['Name'],len(paraseq),self.stat['MaxPara']))
                paraseq = paraseq[:self.stat['MaxPara']]
            if paraseq and orthseq:
                idseq.saveFasta(seqs=orthseq,seqfile=mfile)
                idseq.saveFasta(seqs=paraseq,seqfile=pfile)
                newblast.setStr({'DBase':mfile,'Name':pbfile,    # Will generate this now
                                 'InFile':pfile})
                newblast.formatDB(protein=not self.opt['DNA'])
                newblast.blast(cleandb=True)                    # BLAST and cleanup paralogue vs minsim
                newblast.readBLAST(clear=True,gablam=True)     # Add paralogue BLASTs onto Query BLAST
                for search in newblast.queries():
                    hitdict[search] = newblast.hitToSeq(idseq,query=search)    # Dictionary of {Hit Object:Sequence Object}
                ## ~ [5a] GABLAMO Matrix and Best Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                _stage = 'GABLAMO Matrix'
                for pseq in paraseq:
                    psearch = pseq.shortName()   # newblast.search[paraseq.index(pseq)+1]
                    plen = pseq.aaLen()
                    for hit in newblast.queryHits(psearch):
                        hentry = newblast.hitData(hit,psearch)
                        hseq = hitdict[psearch][hit]
                        if hseq == pseq:
                            if gkey == 'Score': gablamo_matrix.addDis(pseq,pseq,hentry['BitScore'])
                            else: gablamo_matrix.addDis(pseq,pseq,100.0)
                            continue
                        gdict = newblast.gablamData(hit,psearch)    # hit.globalFromLocal(plen)
                        qvh = 0
                        hvq = 0
                        if gkey == 'Score': qvh = hentry['BitScore']
                        elif plen > 0: qvh = float(100 * gdict['Query'][gkey]) / float(plen)
                        else: self.errorLog('%s Paralogue %s has zero length!' % (self.info['Name'],pseq.shortName()),False,False)
                        if gkey == 'Score': hvq = hentry['BitScore']
                        elif hentry['Length'] > 0: hvq = float(100 * gdict['Hit'][gkey]) / float(hentry['Length'])
                        else: self.log.errorLog('%s Paralogue %s Hit %s has zero length!' % (self.info['Name'],pseq.shortName(),hit),False,False)
                        gablamo_matrix.addDis(pseq,hseq,qvh)
                        gablamo_matrix.addDis(hseq,pseq,hvq)
                        if hseq.info['SpecCode'] != spec:       # Hit dif species
                            if hvq > best_qspec[hseq]['Sim']:   # New best hit of query species
                                best_qspec[hseq] = {'Sim':hvq,'Rec':qvh,'Seq':pseq}
            #i# Should now have an all-by-all GABLAMO Matrix and best hit data

            ### ~ [6] Paralign Option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Paralign'
            if not (qry_paraseq and self.opt['Paralign']):   # No need to perform
                ## ~ [6a] MinSim Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for seq in qry_paraseq[0:]:
                    qvh = gablamo_matrix.getDis(qseq,seq)
                    hvq = gablamo_matrix.getDis(seq,qseq)
                    badsim = False
                    if self.info['SimFocus'] in ['query','both'] and qvh < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] in ['hit','both'] and hvq < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] == 'either' and qvh < self.stat['MinSim'] and qvh < self.stat['MinSim']: badsim = True
                    if badsim: qry_paraseq.remove(seq)
                ## ~ [6b] ReOrder by similarity to Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                neworder = []
                rankdict = {}           
                for pseq in qry_paraseq:
                    rkey = '%s-%s-%s' % (rje.preZero(int(gablamo_matrix.getDis(qseq,pseq) * 1000),100000),rje.preZero(int(gablamo_matrix.getDis(pseq,qseq) * 1000),100000),pseq.shortName())
                    rankdict[rkey] = pseq
                for rkey in rje.sortKeys(rankdict,revsort=True):
                    neworder.append(rankdict[rkey])
                qry_paraseq = neworder[0:]
                ## ~ [6c] Splice variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.opt['ParaSplice']:
                    baseacc = []
                    if rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum']): baseacc.append(rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum'])[0])
                    else: baseacc.append(qseq.info['AccNum'])
                    for pseq in qry_paraseq:
                        pacc = pseq.info['AccNum']
                        if rje.matchExp('^(\S+)\-(\d+)$',pacc): pacc = rje.matchExp('^(\S+)\-(\d+)$',pacc)[0]
                        if pacc in baseacc: # Splice variant
                            qry_paraseq.remove(pseq)
                        else: baseacc.append(pacc)
                ## ~ [6d] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cfile = self.gFile('PARA','closepara.id')   # Used by self._parAlign()
                open(cfile,'w').write(string.join(rje_seq.seqInfoList(qry_paraseq),'\n'))
                self._parAlign(qseq,qry_paraseq)

            ### ~ [7] Orthologue assignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Orthologues'
            ## To be considered an orthologue:
            ## -1- QvH must be >= QvH for any other sequence of same species
            ## -2- HvQ must be > HvP
            ## -3- QvH must be > QvP
            ## ~ [7a] Establish Query's closest paralogue ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if best_hspec.has_key(spec):        # Paralogues to consider
                cseq = best_hspec[spec]['Seq']  # This is the Query's closest paralogue - used for PostDup option
            else: cseq = None
            ## ~ [7b] Work through each potential orthologue ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for hseq in idseq.seq:
                hspec = hseq.info['SpecCode']   # Species code for hit sequence
                if hspec == spec: continue      # Paralogue - do not consider here                    
                qvh = gablamo_matrix.getDis(qseq,hseq)  # Query vs Hit similarity
                hvq = gablamo_matrix.getDis(hseq,qseq)  # Hit vs Query similarity
                hvp = -1        # Hit vs Paralogue Similarity
                pseq = None
                if not best_qspec.has_key(hseq):  #!# Why not??!! Looking at all idseq sequences: should all be in best_qseq #!#
                    self.log.errorLog('Hit "%s" behaving strangely for %s: excluded.' % (hseq.shortName(),self.info['Name']),printerror=False,quitchoice=False)
                    best_hspec[hspec] = {'Seq':None}
                    stickhits.append(hseq.shortName())   
                    continue
                if best_qspec[hseq]['Seq'] != qseq:     # This should be the closest Query species seq to Hit
                    if self.opt['Reciprocal']:          # Simple Reciprocal Best Hit Method: rejected!
                        best_hspec[hspec]['Seq'] = None  
                        continue
                    pseq = best_qspec[hseq]['Seq']
                    hvp = gablamo_matrix.getDis(hseq,pseq)
                    pvh = gablamo_matrix.getDis(pseq,hseq)
                    pvq = gablamo_matrix.getDis(pseq,qseq)
                    qvp = gablamo_matrix.getDis(qseq,pseq)
                else: continue                          # Reciprocal best hits
                ### ~~~ Strict within orthologous clade? ~~~ ###
                if self.opt['PostDup'] and cseq:
                    cvq = gablamo_matrix.getDis(cseq,qseq)
                    qvc = gablamo_matrix.getDis(qseq,cseq)
                    cvh = gablamo_matrix.getDis(cseq,hseq)
                    hvc = gablamo_matrix.getDis(hseq,cseq)
                    if (hvc >= hvq or qvc >= qvh):    # Query or Hit closer to paralogue
                        best_hspec[hspec]['Seq'] = None    
                        continue
                ### ~~~ May be outside of clade but not within a different paralogous clade ~~~ ###
                ## >> Hit closer to query than to paralogue. Query closer to hit than to paralogue.
                if hvq > hvp and qvh > qvp: continue    # Orthologue within query clade
                ## >> Query closer to paralogue than to hit. Paralogue closer to query than to hit.
                elif pvq > pvh and qvp > qvh:           # Orthologue outside of query/paralogue duplication clade
                    stickhits.append(pseq.shortName())                    
                    continue
                ## >> Otherwise, may be in clade with paralogue to exclusion of Query
                else: best_hspec[hspec]['Seq'] = None
            ## ~ [7c] Recreate list, ordered by similarity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            simseq = {} # Dictionary of {sim_spec:seq}
            for hspec in best_hspec.keys():
                if hspec != spec and best_hspec[hspec]['Seq']:    #!# Why use Sim and Seq, why not just store seq and use gablamo_matrix? #!#
                    newkey = '%.6f_%s' % ((string.atof(best_hspec[hspec]['Sim'])/100.0),best_hspec[hspec]['Seq'])
                    simseq[newkey] = best_hspec[hspec]['Seq']
            ## ~ [7d] Generate sequence list and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orths = [qseq]
            for key in rje.sortKeys(simseq,revsort=True): orths.append(simseq[key])
            idseq.info['Name'] = ofile
            idseq.seq = orths[0:]
            idseq.saveFasta()
            ## ~ [7e] OrthID List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orthid = []
            for seq in orths: orthid.append(seq.shortName())
            open(afile,'w').write(string.join(rje.sortUnique(orthid),'\n'))
            ## ~ [7f] Optional Sticky list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['Sticky']: open(self.gFile('ORTH','sticky.id'),'w').write(string.join(rje.sortUnique(stickhits),'\n'))

            ### ~ [8] SaveSpace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['DNA'] and os.path.exists(pfile): os.unlink(pfile)
            if self.opt['SaveSpace']:
                for ssfile in [pbfile,mfile]:
                    if os.path.exists(ssfile): os.unlink(ssfile); ssx += 1
                self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
            return True
        except:
            self.errorLog('Major Problem with GopherFork(%s).orthFas() %s' % (self.info['Name'],_stage))
            self.deBug('Kill me now, or I will go on!')
            self.opt['OrthAlign'] = False   # Nothing to align!
            return False
#########################################################################################################################
    def _orthFasLegacy(self):     ### Makes an unaligned Fasta file of potential orthologues.                             #V3.0
        '''
        Makes an unaligned Fasta file of potential orthologues.
        > Save reduced sequences as ORTH/AccNum.orth.fas
        '''
        try:### ~ [1] Setup Variables needed for all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Setup'
            ## ~ [1a] Setup variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry = self.obj['Sequence']
            spec = qry.info['SpecCode']
            stickhits = []
            ssx = 0     # Counter for save space deletions
            ## ~ [1b] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ifile = self.gFile('BLAST','blast.id')  # Generated by _orthBlast()
            qfile = self.gFile('BLAST','qry')       # Generated by _orthBlast()
            bfile = self.gFile('ORTH','orth.blast') # Will generate this now
            ffile = self.gFile('ORTH','full.fas')   # Generated now but cleared by savespace
            ofile = self.gFile('ORTH','orth.fas')   # Generated now
            afile = self.gFile('ORTH','orth.id')    # Generated now
            mfile = self.gFile('ORTH','minsim.fas') # Clean this file up
            pfile = self.gFile('PARA','para.fas')   # Define here as will check to delete later
            pbfile = self.gFile('PARA','para.blast')# Clean this file up
            rje.mkDir(self,ofile,log=True)
            rje.mkDir(self,pfile,log=True)
            ## ~ [1c] Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blast = self.obj['BLAST']
            blast.setStr({'Name':bfile})
            blast.search = []
            ## ~ [1d] Setup GABLAMO Matrix and stats. Used for comparing Queries and Hits ~~~~~~~~~ ##
            gablamo_matrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            gablamo_matrix.opt['Symmetric'] = False
            if self.info['GABLAMO Key'] not in ['Sim','ID','Len','Score']:
                self.errorLog('GABLAMO Key (%s) not Sim/ID/Len. Changed to Sim.' % self.info['GABLAMO Key'],printerror=False,quitchoice=True)
                self.info['GABLAMO Key'] = 'Sim'
            if self.info['GABLAMO Key'] != 'Score': gkey = 'GABLAMO %s' % self.info['GABLAMO Key']
            else: gkey = self.info['GABLAMO Key']

            ### ~ [2] Read in original BLAST Hits and filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Read BLAST Hits'
            ## ~ [2a] Check for ID list file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(ifile):
                if self._orthBlast(self.opt['CompFilter']): return self._orthFas()    # Try _orthBlast again and then try running self
                elif not self.opt['DropOut']:
                    open(ofile,'w').write('>%s\n%s\n' % (qry.info['Name'],qry.info['Sequence']))
                    return True
                else:                   # BLAST Unsuccessful
                    self.opt['OrthAlign'] = False   
                    return False
            ## ~ [2b] Read in proteins hit from original BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=T','unkspec=F','specnr=F','gnspacc=T'])
            idseq.seq = []
            for id in open(ifile,'r').readlines():
                id = rje.chomp(id)
                if id: idseq.seqFromFastaCmd(id=id,dbase=blast.getStr('DBase'))
            ## ~ [2c] Filter sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if idseq.opt['AutoFilter']: idseq.autoFilter(cmd_list=self.cmd_list+['filterout=None','seqout=None'])
            ## ~ [2d] Check Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if idseq.querySeq(query=qry.shortName()):  # Query not in BLAST IDs
                idseq.removeSeq(text='Query sequence in BLAST database',seq=idseq.obj['QuerySeq'])
                self.printLog('#QRY','Query sequence %s replaced with Input sequence.' % qry.shortName())
            else: self.printLog('#QRY','Query sequence %s added to BLAST Hits' % qry.shortName())
            idseq._addSeq(qry.info['Name'],qry.info['Sequence'])    ### Adds a new Sequence Object to list
            idseq.seq = idseq.seq[-1:] + idseq.seq[:-1]
            idseq.querySeq(query=qry.shortName())       # Make Query 
            qseq = idseq.obj['QuerySeq']
            ## ~ [2e] Save filtered sequences for second BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idseq.saveFasta(seqfile=ffile)
            if not rje.exists(qfile): open(qfile,'w').write('>%s\n%s\n' % (qseq.info['Name'],qseq.info['Sequence']))

            ### ~ [3] New BLAST to get full GABLAM stats for GOPHER calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'New BLAST'
            ## ~ [3a] Make new BLAST database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newblast = rje_blast.blastObj(log=self.log,cmd_list=self.cmd_list,type='Dev')
            newblast.setInfo(blast.getStr())
            newblast.setOpt(blast.getBool())
            newblast.setInt(blast.getInt())
            newblast.setNum(blast.getNum())
            newblast.setStr({'InFile':qfile,'DBase':ffile})
            newblast.setInt({'OneLine':idseq.seqNum(),'HitAln':idseq.seqNum()})
            newblast.setNum({'E-Value':10})
            newblast.setBool({'Complexity Filter':False,'Composition Statistics':False})
            newblast.formatDB(protein=not self.opt['DNA'])
            newblast.blast(cleandb=True)    # BLAST vs homologues and remove database files
            ## ~ [3b] Read BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newblast.readBLAST(clear=True,gablam=True)
            search = newblast.search[0]
            if self.opt['SaveSpace']:
                for ssfile in [ffile,blast.getStr('Name'),newblast.getStr('Name')]:
                    if os.path.exists(ssfile): os.unlink(ssfile); ssx += 1
            if search.hitNum() < 1: # No hits, including self
                self.errorLog('%s has no BLAST hits (including self!).' % self.info['Name'],printerror=False)
                self.opt['OrthAlign'] = False   # Nothing to align!
                if self.opt['SaveSpace']: self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
                return False

            ### ~ [3c] Map sequences onto hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'MapSeq'
            hitdict = {}        # This will start as a dictionary of {search:hitdict} and have paralogue searches added
            hitdict[search] = search.hitSeq(idseq,proglog=False)    # Dictionary of {Hit Object:Sequence Object}
            if not qseq or qseq not in hitdict[search].values():    # No self hit
                self.errorLog('Cannot find self hit for %s in reduced BLAST hits.' % self.info['Name'],printerror=False)
                if self.opt['SaveSpace']: self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
                if not self.opt['DropOut']:
                    open(ofile,'w').write('>%s\n%s\n' % (qry.info['Name'],qry.info['Sequence']))
                    return True
                self.opt['OrthAlign'] = False   # Nothing to align!
                return False
            
            ### ~ [4] Similarity Filter, Best Species Hit & Paralogues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Similarity Filter'
            ## Best spec dictionaries store the Sim scores (and Reciprocal score) for the best hit of that sequence:
            ## - best_qspec = {seq:best hit of query species}
            ## - best_hspec = {hspec:best hit for query of hspec sequences}
            best_qspec = {}     # Dictionary of {seq:{'Sim':X,'Rec':X,'Seq':seq}}
            best_hspec = {}     # Dictionary of {hspec:{'Sim':X,'Rec':X,'Seq':seq}}
            paraseq = []        # Paralogous sequences
            baseacc = []        # Base AccNums for spotting splice variants
            if not self.opt['ParaSplice']:  # Add query base accnum to baseacc to filter splice variants of query
                if rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum']): baseacc.append(rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum'])[0])
                else: baseacc.append(qseq.info['AccNum'])
            simlost = 0
            qlen = qseq.aaLen()
            ## ~ [4a] Go through each hit in turn and establish similarity and best hits ~~~~~~~~~~ ##
            for hit in search.hit[0:]:
                ## ~~~ GABLAMO ~~~ ##
                hseq = hitdict[search][hit]
                if hseq == qseq:
                    if gkey == 'Score': gablamo_matrix.addDis(qseq,qseq,hit.stat['BitScore'])
                    else: gablamo_matrix.addDis(qseq,qseq,100.0)
                    continue
                if not hseq:
                    self.errorLog('Unable to match idseq sequence for %s hit %s' % (self.info['Name'],hit.info['Name']),False,False)
                    self.printLog('#CHECK','Check %s files for problem hit %s' % (self.info['Name'],hit.info['Name']),screen=False)
                    search.hit.remove(hit)
                hspec = hseq.info['SpecCode']
                gdict = hit.globalFromLocal(qlen)
                hit.aln = []    #!# Clear Aln objects to save memory - move to rje_blast GABLAM at some point? #!#
                if gkey == 'Score': qvh = hit.stat['BitScore']
                else: qvh = float(100 * gdict['Query'][gkey]) / float(qlen)
                hvq = 0
                if gkey == 'Score': hvq = hit.stat['BitScore']
                elif hit.stat['Length'] > 0: hvq = float(100 * gdict['Hit'][gkey]) / float(hit.stat['Length'])
                else: self.errorLog('%s Hit %s has zero length!' % (self.info['Name'],hit.info['Name']),False,False)
                ## ~~~ Paralogue ~~~ ##
                if hspec == spec:   # Paralogue or splice variant
                    if not self.opt['ParaSplice']:
                        pacc = hseq.info['AccNum']
                        if rje.matchExp('^(\S+)\-(\d+)$',pacc):
                            pacc = rje.matchExp('^(\S+)\-(\d+)$',pacc)[0]
                        if pacc in baseacc:     # Splice variant => Toss Hit
                            search.hit.remove(hit)
                            idseq.removeSeq(text='Splice variant (%s)' % pacc,seq=hseq)
                            continue                   
                    paraseq.append(hseq)
                ## ~~~ MinSim Filter: Only MinSim potential orthlogues so that all paralogues are used for comparison ~~~ ##
                badsim = False
                if gkey != 'Score': 
                    if self.info['SimFocus'] in ['query','both'] and qvh < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] in ['hit','both'] and hvq < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] == 'either' and qvh < self.stat['MinSim'] and qvh < self.stat['MinSim']: badsim = True
                if badsim and hspec != spec:
                    simlost += 1
                    search.hit.remove(hit)
                    while hseq in idseq.seq: idseq.seq.remove(hseq)     #?# Can a sequence be in here twice?! #?#
                    continue
                ## ~~~ Update best_hit ~~~ ##
                gablamo_matrix.addDis(qseq,hseq,qvh)
                gablamo_matrix.addDis(hseq,qseq,hvq)
                best_qspec[hseq] = {'Sim':hvq,'Rec':qvh,'Seq':qseq}     # At this point, the Query is the best qspec sequence for hseq
                if best_hspec.has_key(hspec) and (qvh < best_hspec[hspec]['Sim'] or (qvh == best_hspec[hspec]['Sim'] and hvq <= best_hspec[hspec]['Rec'] and best_hspec[hspec]['Seq'])):
                    continue    # Equal or better already exists
                # Best so far for that species...
                best_hspec[hspec] = {'Sim':qvh,'Rec':hvq,'Seq':hseq}
            if gkey == 'Score': self.printLog('#SIM','%s: %d hits. Using BitScore, so no Min. Similarity threshold.' % (self.info['Name'],search.hitNum()))
            else: self.printLog('#SIM','%s: %d of %d hits failed to meet Min. Similarity threshold (%s %.2f%%) - removed.' % (self.info['Name'],simlost,(search.hitNum()+simlost),gkey,self.stat['MinSim']))
            ## ~ [4b] Remove all but best hits of each species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orthseq = []    # List of sequences to compare paralogues to
            for hit in search.hit[0:]:
                hseq = hitdict[search][hit]
                hspec = hseq.info['SpecCode']
                if hspec == spec: continue  # Paralogue 
                if best_hspec.has_key(hspec) and best_hspec[hspec]['Seq'] == hseq:    # Best of species
                    orthseq.append(hseq)
                else:   # Bye bye...!
                    stickhits.append(hseq.shortName())
                    try:
                        search.hit.remove(hit)
                        idseq.seq.remove(hseq)
                    except:
                        self.log.errorLog('Problem during remove (%s/%s).' % (hit.info['Name'],hseq.shortName()),printerror=False)
            self.printLog('#SIM','%s: %d best hits for species.' % (self.info['Name'],len(orthseq)))
            ## ~ [4c] Special Reciprocal Best Hit Method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qry_paraseq = paraseq[0:]
            if self.getBool('Reciprocal') and self.getBool('FullRBH'):
                px = len(paraseq)
                self.printLog('#RBH','Performing GopherBLAST searches for %d putative orthologues' % len(orthseq))
                for orth in orthseq:
                    orthdir = self.obj['Parent'].gopherDir(orth)
                    orthfork = GopherFork(log=self.log,cmd_list=self.cmd_list+['i=-1','errorlog=None'])
                    orthfork.setInfo({'Name':orth.shortName(),'OutPath':orthdir})
                    orthfork.obj['Sequence'] = orth
                    orthfork.obj['BLAST'] = self.obj['BLAST']
                    orthfork.obj['Parent'] = self.obj['Parent']
                    orthfork.log = self.log
                    orthfork.run('orthblast')
                    orthidseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None','accnr=T','unkspec=F','specnr=F','gnspacc=T'])
                    orthidseq.seq = []
                    for id in open(orthfork.gFile('BLAST','blast.id'),'r').readlines():
                        id = rje.chomp(id)
                        if id: orthidseq.seqFromFastaCmd(id=id,dbase=blast.getStr('DBase'))
                    for oseq in orthidseq.seq:
                        if qseq.sameSpec(oseq) and oseq not in paraseq: paraseq.append(oseq)
                self.printLog('#RBH','%d additional paralogues added.' % (len(paraseq)-px))

            ### ~ [5] ParaSeq BLAST to establish Orthologue-Paralogue relationships ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# When updating to make full reciprocal BLAST version, this will need to call full GOPHER BLAST for minsim sequences #!#
            _stage = 'ParaSeq BLAST'
            self.printLog('#PARA','%s: %d paralogues to consider.' % (self.info['Name'],len(paraseq)))
            if len(paraseq) > self.stat['MaxPara']:
                self.printLog('#PARA','%s: %d paralogues reduced to %d (MaxPara).' % (self.info['Name'],len(paraseq),self.stat['MaxPara']))
                paraseq = paraseq[:self.stat['MaxPara']]
            if paraseq and orthseq:
                idseq.saveFasta(seqs=orthseq,seqfile=mfile)
                newblast.info['DBase'] = mfile
                idseq.saveFasta(seqs=paraseq,seqfile=pfile)
                newblast.info['Name'] = pbfile    # Will generate this now
                newblast.info['InFile'] = pfile
                newblast.formatDB(protein=not self.opt['DNA'])
                newblast.blast(cleandb=True)        # BLAST and cleanup paralogue vs minsim
                newblast.readBLAST(clear=False,gablam=True)    # Add paralogue BLASTs onto Query BLAST
                for search in newblast.search:
                    hitdict[search] = search.hitSeq(idseq,proglog=False)    # Dictionary of {Hit Object:Sequence Object}
                ## ~ [5a] GABLAMO Matrix and Best Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                _stage = 'GABLAMO Matrix'
                for pseq in paraseq:
                    psearch = newblast.search[paraseq.index(pseq)+1]
                    plen = pseq.aaLen()
                    for hit in psearch.hit:
                        hseq = hitdict[psearch][hit]
                        if hseq == pseq:
                            if gkey == 'Score': gablamo_matrix.addDis(pseq,pseq,hit.stat['BitScore'])
                            else: gablamo_matrix.addDis(pseq,pseq,100.0)
                            continue
                        gdict = hit.globalFromLocal(plen)
                        hit.aln = []    #!# Clear Aln objects to save memory - move to rje_blast GABLAM at some point? #!#
                        qvh = 0
                        hvq = 0
                        if gkey == 'Score': qvh = hit.stat['BitScore']
                        elif plen > 0: qvh = float(100 * gdict['Query'][gkey]) / float(plen)
                        else: self.log.errorLog('%s Paralogue %s has zero length!' % (self.info['Name'],pseq.shortName()),False,False)
                        if gkey == 'Score': hvq = hit.stat['BitScore']
                        elif hit.stat['Length'] > 0: hvq = float(100 * gdict['Hit'][gkey]) / float(hit.stat['Length'])
                        else: self.log.errorLog('%s Paralogue %s Hit %s has zero length!' % (self.info['Name'],pseq.shortName(),hit.info['Name']),False,False)
                        gablamo_matrix.addDis(pseq,hseq,qvh)
                        gablamo_matrix.addDis(hseq,pseq,hvq)
                        if hseq.info['SpecCode'] != spec:       # Hit dif species
                            if hvq > best_qspec[hseq]['Sim']:   # New best hit of query species
                                best_qspec[hseq] = {'Sim':hvq,'Rec':qvh,'Seq':pseq}
            #i# Should now have an all-by-all GABLAMO Matrix and best hit data

            ### ~ [6] Paralign Option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Paralign'
            if not (qry_paraseq and self.opt['Paralign']):   # No need to perform
                ## ~ [6a] MinSim Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for seq in qry_paraseq[0:]:
                    qvh = gablamo_matrix.getDis(qseq,seq)
                    hvq = gablamo_matrix.getDis(seq,qseq)
                    badsim = False
                    if self.info['SimFocus'] in ['query','both'] and qvh < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] in ['hit','both'] and hvq < self.stat['MinSim']: badsim = True
                    if self.info['SimFocus'] == 'either' and qvh < self.stat['MinSim'] and qvh < self.stat['MinSim']: badsim = True
                    if badsim: qry_paraseq.remove(seq)
                ## ~ [6b] ReOrder by similarity to Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                neworder = []
                rankdict = {}           
                for pseq in qry_paraseq:
                    rkey = '%s-%s-%s' % (rje.preZero(int(gablamo_matrix.getDis(qseq,pseq) * 1000),100000),rje.preZero(int(gablamo_matrix.getDis(pseq,qseq) * 1000),100000),pseq.shortName())
                    rankdict[rkey] = pseq
                for rkey in rje.sortKeys(rankdict,revsort=True):
                    neworder.append(rankdict[rkey])
                qry_paraseq = neworder[0:]
                ## ~ [6c] Splice variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.opt['ParaSplice']:
                    baseacc = []
                    if rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum']): baseacc.append(rje.matchExp('^(\S+)\-(\d+)$',qseq.info['AccNum'])[0])
                    else: baseacc.append(qseq.info['AccNum'])
                    for pseq in qry_paraseq:
                        pacc = pseq.info['AccNum']
                        if rje.matchExp('^(\S+)\-(\d+)$',pacc): pacc = rje.matchExp('^(\S+)\-(\d+)$',pacc)[0]
                        if pacc in baseacc: # Splice variant
                            qry_paraseq.remove(pseq)
                        else: baseacc.append(pacc)
                ## ~ [6d] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cfile = self.gFile('PARA','closepara.id')   # Used by self._parAlign()
                open(cfile,'w').write(string.join(rje_seq.seqInfoList(qry_paraseq),'\n'))
                self._parAlign(qseq,qry_paraseq)

            ### ~ [7] Orthologue assignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Orthologues'
            ## To be considered an orthologue:
            ## -1- QvH must be >= QvH for any other sequence of same species
            ## -2- HvQ must be > HvP
            ## -3- QvH must be > QvP
            ## ~ [7a] Establish Query's closest paralogue ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if best_hspec.has_key(spec):        # Paralogues to consider
                cseq = best_hspec[spec]['Seq']  # This is the Query's closest paralogue - used for PostDup option
            else: cseq = None
            ## ~ [7b] Work through each potential orthologue ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for hseq in idseq.seq:
                hspec = hseq.info['SpecCode']   # Species code for hit sequence
                if hspec == spec: continue      # Paralogue - do not consider here                    
                qvh = gablamo_matrix.getDis(qseq,hseq)  # Query vs Hit similarity
                hvq = gablamo_matrix.getDis(hseq,qseq)  # Hit vs Query similarity
                hvp = -1        # Hit vs Paralogue Similarity
                pseq = None
                if not best_qspec.has_key(hseq):  #!# Why not??!! Looking at all idseq sequences: should all be in best_qseq #!#
                    self.log.errorLog('Hit "%s" behaving strangely for %s: excluded.' % (hseq.shortName(),self.info['Name']),printerror=False,quitchoice=False)
                    best_hspec[hspec] = {'Seq':None}
                    stickhits.append(hseq.shortName())   
                    continue
                if best_qspec[hseq]['Seq'] != qseq:     # This should be the closest Query species seq to Hit
                    if self.opt['Reciprocal']:          # Simple Reciprocal Best Hit Method: rejected!
                        best_hspec[hspec]['Seq'] = None  
                        continue
                    pseq = best_qspec[hseq]['Seq']
                    hvp = gablamo_matrix.getDis(hseq,pseq)
                    pvh = gablamo_matrix.getDis(pseq,hseq)
                    pvq = gablamo_matrix.getDis(pseq,qseq)
                    qvp = gablamo_matrix.getDis(qseq,pseq)
                else: continue                          # Reciprocal best hits
                ### ~~~ Strict within orthologous clade? ~~~ ###
                if self.opt['PostDup'] and cseq:
                    cvq = gablamo_matrix.getDis(cseq,qseq)
                    qvc = gablamo_matrix.getDis(qseq,cseq)
                    cvh = gablamo_matrix.getDis(cseq,hseq)
                    hvc = gablamo_matrix.getDis(hseq,cseq)
                    if (hvc >= hvq or qvc >= qvh):    # Query or Hit closer to paralogue
                        best_hspec[hspec]['Seq'] = None    
                        continue
                ### ~~~ May be outside of clade but not within a different paralogous clade ~~~ ###
                ## >> Hit closer to query than to paralogue. Query closer to hit than to paralogue.
                if hvq > hvp and qvh > qvp: continue    # Orthologue within query clade
                ## >> Query closer to paralogue than to hit. Paralogue closer to query than to hit.
                elif pvq > pvh and qvp > qvh:           # Orthologue outside of query/paralogue duplication clade
                    stickhits.append(pseq.shortName())                    
                    continue
                ## >> Otherwise, may be in clade with paralogue to exclusion of Query
                else: best_hspec[hspec]['Seq'] = None
            ## ~ [7c] Recreate list, ordered by similarity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            simseq = {} # Dictionary of {sim_spec:seq}
            for hspec in best_hspec.keys():
                if hspec != spec and best_hspec[hspec]['Seq']:    #!# Why use Sim and Seq, why not just store seq and use gablamo_matrix? #!#
                    newkey = '%.6f_%s' % ((string.atof(best_hspec[hspec]['Sim'])/100.0),best_hspec[hspec]['Seq'])
                    simseq[newkey] = best_hspec[hspec]['Seq']
            ## ~ [7d] Generate sequence list and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orths = [qseq]
            for key in rje.sortKeys(simseq,revsort=True): orths.append(simseq[key])
            idseq.info['Name'] = ofile
            idseq.seq = orths[0:]
            idseq.saveFasta()
            ## ~ [7e] OrthID List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            orthid = []
            for seq in orths: orthid.append(seq.shortName())
            open(afile,'w').write(string.join(rje.sortUnique(orthid),'\n'))
            ## ~ [7f] Optional Sticky list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['Sticky']: open(self.gFile('ORTH','sticky.id'),'w').write(string.join(rje.sortUnique(stickhits),'\n'))

            ### ~ [8] SaveSpace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['DNA'] and os.path.exists(pfile): os.unlink(pfile)
            if self.opt['SaveSpace']:
                for ssfile in [pbfile,mfile]:
                    if os.path.exists(ssfile): os.unlink(ssfile); ssx += 1
                self.printLog('#MEM','%d files deleted to save disc space.' % (ssx))
            return True
        except:
            self.errorLog('Major Problem with GopherFork(%s).orthFas() %s' % (self.info['Name'],_stage))
            self.deBug('Kill me now, or I will go on!')
            self.opt['OrthAlign'] = False   # Nothing to align!
            return False
#########################################################################################################################
    def _orthAln(self):     ### Makes an aligned Fasta file of potential orthologues.                               #V3.0
        '''
        Makes an aligned Fasta file of potential orthologues.
        5. Align sequences with MUSCLE
            > ALN/AccNum.orthaln.fas
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# What are sfile and stat used for? Don't seem to do anything!
            #?# Add stats and stats output for REST output?
            ## ~ [1a] ~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fasname = self.gFile('ORTH','orth.fas')
            sfile = self.gFile('ALN','gopher_alnfas')
            afile = self.gFile('ALN','orthaln.fas')
            rje.mkDir(self,afile,log=True)
            ## ~ [1b] ~ Stats wanted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            stat = {}
            stat['FILE'] = sfile
            stat['TDT'] = 'gopher_orthaln'
            stats = ['gopher_id','aln_len']
            stat['gopher_id'] = self.obj['Sequence'].info['AccNum']     
            stat['aln_len'] = 0         
            ### ~ [2] ~ Read Fasta file & Align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            orths = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['alnprog=%s' % self.info['AlnProg'],'accnr=T','unkspec=F','specnr=F','gnspacc=T','seqin=%s' % fasname,'memsaver=F'])
            orthaln = orths.align(self.info['AlnProg'],afile)
            #orthaln = orths.muscleAln(outfile=afile)
            #if orthaln == None: orthaln = orths.clustalAln(outfile=afile)
            orths.mapSeq(orthaln)
            #orths.tidyXGaps() # ??? Why do this?! If needing to remove 100% gaps: orths.stripGap(orths.seqNum(),gaps=['-'])
            orths.saveFasta(seqfile=afile)
            stat['aln_len'] = orthaln.seq[0].seqLen()
            return True
        except: self.errorLog('Major Problem with GopherFork(%s).orthAln()' % (self.info['Name'])); raise
#########################################################################################################################
    def _orthTree(self):    ### Makes an unrooted (by default) tree using rje_tree.py                               #V3.0
        '''Makes an unrooted (by default) tree using rje_tree.py.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            acc = self.obj['Sequence'].info['AccNum']
            afile = self.gFile('ALN','orthaln.fas')
            if not os.path.exists(afile):   # No alignment file to make tree from
                return False
            tfile = self.gFile('TREE','orth.nwk')
            rje.mkDir(self,tfile,log=True)
            ## ~ [1a] Check need to run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists(tfile) and not self.opt['Force'] and not self.opt['FullForce']:   # May be no need to run
                if rje.isYounger(tfile,afile) == tfile or self.opt['IgnoreDate']: return True   # Already done!

            ### ~ [2] Read Fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list + ['seqin=%s' % afile,'autoload=T','autofilter=F','memsaver=F','accnr=T','unkspec=F','specnr=F','gnspacc=T']
            orths = rje_seq.SeqList(log=self.log,cmd_list=seqcmd)
            if orths.seqNum() < 4:  ### Too small
                self.printLog('#TREE','%s has too few sequence for tree construction.' % afile)
                return False

            ### ~ [3] Make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tree = rje_tree.Tree(log=self.log,cmd_list=['root=mid']+self.cmd_list+['autoload=F','savetree=%s' % tfile])
            tree.makeTree(orths,keepfile=False)
            if tree.opt['OutputBranchLen']: withbranchlengths = 'Length'
            else: withbranchlengths = 'none'
            outnames = tree.info['OutNames']
            maxnamelen = tree.stat['TruncNames']
            tree.info['Basefile'] = self.gFile('TREE','orth')
            try: tree.saveTrees(seqname=outnames,blen=withbranchlengths)
            except: tree.saveTree(seqnum=True,seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            return True
        except: self.errorLog('Major Problem with GopherFork(%s)._orthTree()' % self.info['Name']); raise
#########################################################################################################################
    def _parAlign(self,qseq=None,paraseq=None,gablamo_matrix=None): ### Makes an alignment of close paralogues.     #V3.0
        '''
        Makes an aligned Fasta file of close paralogues: PARALN/accnum.paraln.fas
        >> qseq:Sequence Object (Query)
        >> paraseq:list of Sequence Objects (Paralogues)
        >> gablamo_matrix:dictionary of GABLAMO distances (seq1:seq2)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfile = self.gFile('PARA','para.fas')       # Should exist already
            parout = self.gFile('PARALN','paraln.fas')  # Will generate this now
            ## ~ [1a] Check need to run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.opt['Paralign'] or paraseq == []: return False      # No need to perform
            if os.path.exists(parout) and not self.opt['Force'] and not self.opt['FullForce']:  # May be no need to run
                if rje.isYounger(parout,pfile) == parout or self.opt['IgnoreDate']: return True     # Already done!
            rje.mkDir(self,parout,log=True)
            ## ~ [1b] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ocmd = self.cmd_list + ['alnprog=%s' % self.info['AlnProg'],'seqin=%s' % pfile,'autoload=F','accnr=T','unkspec=F','specnr=F','gnspacc=T']
            outseq = rje_seq.SeqList(log=self.log,cmd_list=ocmd)
            if not qseq: qseq = self.obj['Sequence']
            if paraseq: outseq.seq = [qseq] + paraseq
            else:
                outseq.loadSeqs(aln=False)
                outseq.seq = [qseq] + outseq.seq[0:]
            outseq.info['Name'] = parout
            ### ~ [2] Align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outseq.seqNum() < 2: return False    # No paralogues
            outseq._checkAln(aln=True,realign=True)
            outseq.saveFasta(seqfile=parout)
            return True
        except:
            self.errorLog('Major Problem with GopherFork(%s).parAlign()' % self.info['Name'])
            self.deBug('Kill me now, or I will go on!')
            return False
#########################################################################################################################
    def _paraFam(self):     ### Makes a rooted tree of two paralogous gopher subfams using rje_tree.py              #V3.0
        '''Makes a rooted tree of two paralogous gopher subfams using rje_tree.py.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['DNA'] or not self.opt['ParaFam']: return False
            acc = self.obj['Sequence'].info['AccNum']
            ofile = self.gFile('ORTH','orth.fas')   # Should exist
            pfile = self.gFile('PARA','para.fas')   # Should exist
            for file in [ofile,pfile]:
                if not os.path.exists(file):   # No orthologues and/or paralogues
                    self.printLog('#PARA','No %s - No ParaFam for me!' % file,log=False)
                    return False
            rje.mkDir(self,self.gFile('PARAFAM/','tmp'),log=True)
            ## ~ [1a] Setup Paralogue SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pcmd = self.cmd_list + ['seqin=%s' % pfile,'fasdb=%s' % self.info['OrthDB'],'autoload=F','accnr=T','unkspec=F','specnr=F','gnspacc=T']
            paraseq = rje_seq.SeqList(log=self.log,cmd_list=pcmd)
            paraseq.loadSeqs(aln=False)

            ### ~ [2] Generate paralogue Gophers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tmpfile = 'tmp%s.fas' % rje.randomString(8)
            paraseq.saveFasta(seqfile=tmpfile)
            gcmd = self.cmd_list + ['orthfas','paralign=F','parafam=F','gopherfam=F','seqin=%s' % tmpfile]    #!# Could have (full)force issues #!#
            pgopher = Gopher(self.log,gcmd)
            #self.deBug('ParaFam: Run Gopher on Paralogues:\n%s' % gcmd)
            pgopher.run()
            #self.deBug('ParaFam run of Gopher finished!')
            rje_blast.cleanupDB(self,tmpfile,deletesource=True)

            ### ~ [3] Load query orthologues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qcmd = self.cmd_list + ['seqin=%s' % ofile,'autoload=T','autofilter=F','accnr=T','unkspec=F','specnr=F','gnspacc=T']
            qseqs = rje_seq.SeqList(log=self.log,cmd_list=qcmd)
            orths = qseqs.seq[0:]

            ### ~ [4] Make combined datasets and trees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for para in paraseq.seq[0:]:
                ## ~ [4a] Setup pair ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pacc = para.info['AccNum']
                accsort = [acc,pacc]
                accsort.sort()
                pfile = self.gFile('ORTH','orth.fas',pacc)
                sfile = self.gFile('PARAFAM','parafam.fas','%s+%s' % (accsort[0],accsort[1]))
                if os.path.exists(sfile) and not self.opt['Force'] and not self.opt['FullForce']:  # May be no need to run
                    if (rje.isYounger(sfile,pfile) == sfile and rje.isYounger(sfile,ofile) == sfile) or self.opt['IgnoreDate']:
                        continue
                ## ~ [4b] Query orthologues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                dropx = 0       # No. shared orthologues removed
                qseqs = orths[0:]
                qids = rje_seq.seqInfoList(qseqs)
                ## ~ [4c] Load paralogous subfam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pcmd = self.cmd_list + ['seqin=%s' % pfile,'autoload=T','autofilter=F']     # Already filtered.
                comblist = rje_seq.SeqList(log=self.log,cmd_list=pcmd)
                pseqs = comblist.seq[0:]
                ## ~ [4d] Drop shared orthologues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for seq in pseqs[0:]:
                    if seq.shortName() in qids:
                        dropx += 1
                        qi = qids.index(seq.shortName())
                        qseqs = qseqs[:qi] + qseqs[(qi+1):]
                        qids.remove(seq.shortName())
                        pseqs.remove(seq)
                ## ~ [4e] Combine and align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                comblist.seq = qseqs + pseqs
                comblist.info['Name'] = sfile
                comblist._checkAln(aln=True,realign=True)
                comblist.saveFasta()
                ## ~ [4f] Check Numbers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if comblist.seqNum() < 4:  ### Too small
                    self.printLog('#PARAFAM','%s has too few sequence for tree construction.' % sfile)
                    continue
                ## ~ [4g] Save Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gfile = self.gFile('PARAFAM','parafam.grp','%s+%s' % (acc,pacc))
                GROUPS = open(gfile, 'w')
                GROUPS.write('Group 1: %s query orthologues\n' % acc)
                GROUPS.write('%s\n' % string.join(qids,'\n'))
                GROUPS.write('Group 2: Paralogue %s & orthologues\n' % pacc)
                GROUPS.write('%s\n' % string.join(rje_seq.seqInfoList(pseqs),'\n'))
                GROUPS.close()
                ## ~ [4h] Make Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ofile = self.gFile('PARAFAM','parafam.out','%s+%s' % (acc,pacc))
                OUT = open(ofile,'w')
                OUT.write(string.join(qids,'\n'))
                OUT.close()
                tfile = self.gFile('PARAFAM','parafam.nwk','%s+%s' % (acc,pacc))
                tcmd = self.cmd_list + ['autoload=F','savetree=%s' % tfile,'root=%s' % ofile,'group=%s' % gfile]
                tree = rje_tree.Tree(log=self.log,cmd_list=tcmd)
                tree.info['Basefile'] = rje.baseFile(tfile)
                tree.makeTree(comblist,keepfile=False)
                os.unlink(ofile)
                if tree.opt['OutputBranchLen']: withbranchlengths = 'Length'
                else: withbranchlengths = 'none'
                outnames = tree.info['OutNames']
                maxnamelen = tree.stat['TruncNames'] 
                tree.saveTrees(seqname=outnames,blen=withbranchlengths)  ### Generates all tree file formats selected in normal format
                #tree.saveTree(seqnum=True,seqname=outnames,maxnamelen=maxnamelen,blen=withbranchlengths)
            return True
        except:
            self.errorLog('Major Problem with GopherFork(%s)._paraFam()' % self.info['Name'])
            raise
#########################################################################################################################
    def _gopherFam(self):   ### Makes a combined protein family for protein & paralogues > minsim + orthologues     #V3.0
        '''
        Makes a combined protein family dataset (unaligned) consisting of a protein, all the paralogues > minsim and all
        orthologues for each in a single dataset.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['GopherFam']: return False
            acc = self.obj['Sequence'].info['AccNum']
            ofile = self.gFile('ORTH','orth.fas')       # Should exist
            pfile = self.gFile('PARA','para.fas')       # Should exist
            ffile = self.gFile('SUBFAM','subfam.fas')   # Make now
            rje.mkDir(self,ffile,log=True)
            ## ~ [1a] Check need to run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(ofile) or not os.path.exists(pfile): # No orthologues and/or paralogues
                self.printLog('#PARA','No %s or %s - No gopherFam for me!' % (ofile,pfile),log=False)
                return False
            if os.path.exists(ffile) and not self.opt['Force'] and not self.opt['FullForce']:  # May be no need to run
                if (rje.isYounger(ffile,pfile) == ffile and rje.isYounger(ffile,ofile) == ffile) or self.opt['IgnoreDate']:
                    return True     # Already done!
            ## ~ [1b] Setup Paralogue SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pcmd = self.cmd_list + ['seqin=%s' % pfile,'fasdb=%s' % self.info['OrthDB'],'autoload=F']   # Already filtered
            paraseq = rje_seq.SeqList(log=self.log,cmd_list=pcmd)
            paraseq.loadSeqs(aln=False)

            ### ~ [2] Generate paralogue Gophers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tmpfile = 'tmp%s.fas' % rje.randomString(8)
            paraseq.saveFasta(seqfile=tmpfile)
            gcmd = self.cmd_list + ['orthfas','paralign=F','parafam=F','gopherfam=F','seqin=%s' % tmpfile]    #!# Could have (full)force issues #!#
            pgopher = Gopher(self.log,gcmd)
            #self.deBug('GopherFam: Run Gopher on Paralogues:\n%s' % gcmd)
            pgopher.run()
            #self.deBug('GopherFam run of Gopher finished!')
            rje_blast.cleanupDB(self,tmpfile,deletesource=True)

            ### ~ [3] Load query orthologues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qcmd = self.cmd_list + ['seqin=%s' % ofile,'autoload=T','autofilter=F','accnr=T','unkspec=F','specnr=F','gnspacc=T']
            qseqs = rje_seq.SeqList(log=self.log,cmd_list=qcmd)

            ### ~ [4] Make combined dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            idlist = rje_seq.seqInfoList(qseqs.seq[0:])
            for para in paraseq.seq[0:]:
                ## Load paralogous subfam ##
                pacc = para.info['AccNum']
                pfile = self.gFile('ORTH','orth.fas',pacc)
                pcmd = self.cmd_list + ['seqin=%s' % pfile,'autoload=T','autofilter=F'] # Already filtered
                comblist = rje_seq.SeqList(log=self.log,cmd_list=pcmd)
                ## Drop shared orthologues ##
                for seq in comblist.seq[0:]:
                    if seq.shortName() not in idlist:
                        idlist.append(seq.shortName())
                        qseqs.seqFromFastaCmd(seq.shortName(),self.info['OrthDB'])

            ### ~ [5] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qseqs.saveFasta(seqfile=ffile)
            return True
        except:
            self.log.errorLog('Major Problem with GopherFork(%s)._gopherFam()' % self.info['Name'])
            raise
#########################################################################################################################
### End of SECTION II: GopherFork Class                                                                                 #
#########################################################################################################################
    
                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: SPECIFIC METHODS                                                                                        #
#########################################################################################################################
gophermodes = ['dbases','orthblast','orthfas','orthalign','orthtree']        
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MAIN PROGRAM                                                                                             #
#########################################################################################################################
def runMain():
    ### ~ Basic Setup of Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return 
    ### ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Gopher(log=mainlog,cmd_list=cmd_list).run()
    ### ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
        mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
        sys.exit(1)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except SystemExit: sys.exit(1)
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV
#########################################################################################################################
