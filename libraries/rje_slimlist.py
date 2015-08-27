#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_slimlist
Description:  SLiM dataset manager
Version:      1.7.2
Last Edit:    28/03/15
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a replace for the rje_motiflist module and contains the SLiMList class, a replacement for the
    MotifList class. The primary function of this class is to load and store a list of SLiMs and control generic SLiM
    outputs for such programs as SLiMSearch. This class also controls motif filtering according to features of the motifs
    and/or their occurrences.

    Although not actually designed with standalone functionality in mind, as of V1.3 it is possible to load motifs (inc.
    downloading ELM), filter/process motifs and output data using motifout=FILE and motinfo=FILE if desired.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ###
    motifs=FILE     : File of input motifs/peptides [None]
                      Single line per motif format = 'Name Sequence #Comments' (Comments are optional and ignored)
                      Alternative formats include fasta, SLiMDisc output and raw motif lists.
                      "elm" will download, process and load ELM classes using the ELM API.
    reverse=T/F     : Reverse the motifs - good for generating a test comparison data set [False]
    wildscram=T/F   : Perform a wildcard spacer scrambling - good for generating a test comparison data set [False]
    motifout=FILE   : Name of output file for reformatted/filtered SLiMs (PRESTO format) [None]
    ftout=T/F       : Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    mismatch=LIST   : List of X:Y pairs for mismatch dictionary, where X mismatches allowed for Y+ defined positions []    
    motinfo=FILE    : Filename for output of motif summary table (if desired) [None]
    dna=T/F         : Whether motifs should be considered as DNA motifs [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Input I: Motif Filtering Options ###
    minpos=X        : Min number of defined positions [0]
    minfix=X        : Min number of fixed positions for a motif to contain [0]
    minic=X         : Min information content for a motif (1 fixed position = 1.0) [0.0]
    varlength=T/F   : Whether to motifs can have flexible-length elements [True]
    goodmotif=LIST  : List of text to match in Motif names to keep (can have wildcards) []
    nrmotif=T/F     : Whether to remove redundancy in input motifs [False]

    ### Advanced Input II: Motif reformatting options ###
    trimx=T/F       : Trims Xs from the ends of a motif [False]
    minimotif=T/F   : Input file is in minimotif format and will be reformatted (PRESTO File format only) [False]
    ambcut=X        : Cut-off for max number of choices in ambiguous position to be shown as variant [10]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Output I: Motif Occurrence Statistics Options ###
    slimcalc=LIST   : List of additional statistics to calculate for occurrences - Cons,SA,Hyd,Fold,IUP,Chg []
    winsize=X       : Used to define flanking regions for stats. If negative, will use flanks *only* [0]
    iupath=PATH     : The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    iucut=X         : Cut-off for IUPred results (0.0 will report mean IUPred score) [0.0]
    iumethod=X      : IUPred method to use (long/short) [short]
    percentile=X    : Percentile steps to return in addition to mean [0]
    peptides=T/F    : Whether to output peptide sequences based on motif and winsize [False]

    ### Advanced Output II: Alignment Settings ###
    usealn=T/F      : Whether to search for and use alignemnts where present. [False]
    gopher=T/F      : Use GOPHER to generate missing orthologue alignments in alndir - see gopher.py options [False]
    alndir=PATH     : Path to alignments of proteins containing motifs [./] * Use forward slashes (/)
    alnext=X        : File extension of alignment files, accnum.X [aln.fas]
    protalndir=PATH : Output path for Protein Alignments [ProteinAln/]
    motalndir=PATH  : Output path for Motif Alignments []
    flanksize=X     : Size of sequence flanks for motifs [30]
    xdivide=X       : Size of dividing Xs between motifs [10]
    fullforce=T/F   : Whether to force regeneration of alignments using GOPHER

    ### Advanced Output III: Conservation Parameters ###   
    * see rje_slimcalc.py documentation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seq, rje_slim, rje_slimcalc, rje_zen
import rje_motif_V3 as rje_motif    # Still used for reformatting minimotif
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added SLiMCalc object to replace rje_motif_stats.
    # 0.2 - Added DNA option.
    # 0.3 - Added wildcard spacer scrambling - good for generating a test comparison data set.
    # 0.4 - Added motif-splitting based on ("|") options.
    # 0.5 - Fixed minor typo bug.
    # 0.6 - Added reading of ELM classes download for motifs. Modified splitting functions to cope with non-wildcard runs.
    # 1.0 - Functional module with lower case motif splitting fixed and ? -> .{0,1} replacement.
    # 1.1 - Modified to work with GOPHER V3.0 for alignments.
    # 1.2 - Added some extra functions for CompariMotif Memsaver mode
    # 1.3 - Added auto-download of ELM data.
    # 1.4 - Modified code to be compatible with SLiMCore V2.x objects.
    # 1.5 - Added run() method for slimsuite.py compatibility. Improved split motif handling.
    # 1.6 - Modified to read in new ELM class download file with extra header information. Added varlength=T/F filter.
    # 1.6 - Modified so that filtering one element of a split motif removes all.
    # 1.7.0 - Added direct feeding of motif file content for loading (for REST servers).
    # 1.7.1 - Modified input to allow motif=X in additon to motifs=X.
    # 1.7.2 - Fixed bug that could not accept variable length motifs from commandline. Improved error message.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add defined ambiguities to automatically include. (e.g. MSMS mode)
    # [ ] : Make sure that (A|B) would be replaced with [AB] and not split.
    # [Y] : Replace ? with {0,1}
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_SLiMList', '1.7.2', 'March 2015', '2007')
    description = 'SLiM List Management Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['In development',rje_zen.Zen().wisdom()]
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
### SECTION II: SLiMList Class                                                                                          #
#########################################################################################################################
class SLiMList(rje.RJE_Object):     
    '''
    SLiMList Class. Author: Rich Edwards (2007).

    Info:str
    - AlnDir = Path to alignment files [./]
    - AlnExt = File extensions of alignments: AccNum.X [aln.fas]
    - ConScore = Type of conservation score used:  [abs]
        - abs = absolute conservation of motif: reports percentage of homologues in which conserved
        - prop = conservation of amino acid properties
        - pos = positional conservation: each position treated independently 
        - all = all three methods for comparison purposes
    - MotAlnDir = Directory name for output of protein aligments []
    - Motifs = Name of input motif file [None]
    - MotInfo = Filename for output of motif summary table (if desired) [None]
    - MotifOut = Filename for output of reformatted (and filtered?) motifs in PRESTO format [None]
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix)
    - ProtAlnDir = Directory name for output of protein aligments [ProteinAln/]
    
    Opt:boolean
    - AlnGap = Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore as putative sequence fragments [True]
    - Compare = whether being called by CompariMotif (has some special requirements!)
    - ConsAmb = Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    - ConsInfo = Weight positions by information content [True]
    - DNA = Whether motifs should be considered as DNA motifs [False]
    - FTOut = Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    - FullForce = Whether to force regeneration of alignments using GOPHER
    - Gopher = Use GOPHER to generate missing orthologue alignments in outdir/Gopher - see gopher.py options [False]
    - MiniMotif = Input file is in minimotif format and will be reformatted [False]
    - NRMotif = Whether to remove redundancy in input motifs [False]
    - Peptides = Whether to output peptide sequences based on motif [False]
    - Reverse = Reverse the motifs - good for generating a test comparison data set [False]
    - TrimX = Trims Xs from the ends of a motif
    - UseAln = Whether to look for conservation in alignments
    - WildScram = Perform a wildcard spacer scrambling - good for generating a test comparison data set [False]
    - VarLength = Whether to motifs can have flexible-length elements [True]

    Stat:numeric
    - AmbCut = Cut-off for max number of choices in ambiguous position to be shown as variant [10]
        For mismatches, this is the max number of choices for an ambiguity to be replaced with a mismatch wildcard
    - ConsWeight = Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
    - FlankSize = Size of sequence flanks for motifs in MotifAln [30]
    - MinPos = Minimum number of defined positions in motif [0]
    - MinFix = Min number of fixed positions for a motif to contain [0]
    - MinIC = Min information content for a motif (1 fixed position = 1.0) [0.0]
    - Percentile = Percentile steps to return in addition to mean [25]
    - WinSize = Used to define flanking regions for stats. If negative, will use flanks *only* [0]
    - XDivide = Size of dividing Xs between motifs [10]

    List:list
    - Alphabet = List of letters in alphabet of interest
    - GoodMotif = List of text to match in Motif names to keep (can have wildcards) []
    - Motif = List of rje_slim.SLiM objects
    - OccStats - List of occurrence statistics to calculate []
    - Rejects - List if rejected SLiM core names, used for checking for variant filtering.

    Dict:dictionary
    - ConsSpecLists = Dictionary of {BaseName:List} lists of species codes for special conservation analyses
    - ElementIC = Dictionary of {Position Element:IC}
    - MisMatch = Dictionary of {X mismatches:Y+ defined positions}
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) {}

    Obj:RJE_Objects
    - AAPropMatrix = rje_aaprop.AAPropMatrix object
    - INPUT = Input file handle used for MemSaver mode
    - SeqList = SeqList object of search database
    - SLiMCalc = SLiMCalc object for SLiM/Occurrence attribute calculations
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AlnDir','AlnExt','ConScore','MotAlnDir','Motifs','MotifOut','PosMatrix','ProtAlnDir','MotInfo']
        self.optlist = ['AlnGap','Compare','ConsAmb','ConsInfo','FTOut','FullForce','Gopher','MiniMotif','NRMotif',
                        'Reverse','TrimX','UseAln','Peptides','DNA','WildScram','VarLength']
        self.statlist = ['AmbCut','ConsWeight','FlankSize','MinPos','MinFix','MinIC','Percentile','WinSize','XDivide']
        self.listlist = ['Alphabet','GoodMotif','Motif','OccStats','Rejects']
        self.dictlist = ['ConsSpecLists','ElementIC','MisMatch','PosMatrix']
        self.objlist = ['AAPropMatrix','SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'AlnDir':rje.makePath('',return_blank=False),'AlnExt':'aln.fas','ConScore':'abs',
                      'MotAlnDir':rje.makePath('MotifAln/'),'ProtAlnDir':rje.makePath('ProteinAln/')})
        self.setOpt({'AlnGap':True,'ConsAmb':True,'ConsInfo':True,'VarLength':True})
        self.setStat({'AmbCut':10,'FlankSize':30,'Percentile':25,'XDivide':10})
        self.obj['SLiMCalc'] = rje_slimcalc.SLiMCalc(self.log,self.cmd_list)
        self.obj['SLiMCalc'].setupHeaders()
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] Read in commands from cmd_list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._cmdReadList(cmd,'info',['AlnExt','ConScore'])
                self._cmdReadList(cmd,'file',['Motifs','MotifOut','PosMatrix','MotInfo'])
                self._cmdRead(cmd,'file','Motifs','motif')
                self._cmdReadList(cmd,'path',['AlnDir','MotAlnDir','ProtAlnDir'])
                self._cmdReadList(cmd,'opt',['AlnGap','Compare','ConsAmb','ConsInfo','FTOut','FullForce','Gopher','DNA',
                                             'MiniMotif','NRMotif','Reverse','TrimX','UseAln','Peptides','WildScram',
                                             'VarLength'])
                self._cmdReadList(cmd,'int',['AmbCut','FlankSize','MinPos','MinFix','WinSize','XDivide'])
                self._cmdReadList(cmd,'stat',['ConsWeight','MinIC','Percentile'])
                self._cmdReadList(cmd,'list',['Alphabet','GoodMotif','OccStats'])
                self._cmdReadList(cmd,'cdict',['MisMatch'])
                ### Special ###
                if rje.matchExp('^conspec=(.+)',cmd):
                    conspec = rje.matchExp('^conspec=(.+)',cmd)[0]
                    consfiles = glob.glob(conspec)
                    if len(consfiles) > 0:     # File
                        for cfile in consfiles:
                            self.dict['ConsSpecLists'][rje.baseFile(cfile,True)] = rje.listFromCommand(cfile)
                    else: self.dict['ConsSpecLists']['SPEC'] = rje.listFromCommand(conspec,checkfile=False)
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        ## ~ [1a] Alphabet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not self.list['Alphabet']:
            if self.opt['DNA']: self.list['Alphabet'] = rje_slim.default_nts
            else: self.list['Alphabet'] = rje_slim.default_aas

        ### ~ [2] Make MisMatch dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try: 
            ## ~ [2a] Convert to numbers and get max mm number ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maxmm = 0
            for m in self.dict['MisMatch'].keys()[0:]:
                mx = int(m)
                ax = int(self.dict['MisMatch'].pop(m))
                if ax < mx: (ax,mx) = (mx,ax)   # Entered backwards!
                if mx > 0: self.dict['MisMatch'][mx] = ax
                maxmm = max(maxmm,mx)
            ## ~ [2b] Fill in blanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if maxmm:
                if not 1 in self.dict['MisMatch']: self.dict['MisMatch'][1] = 0
                for mx in range(1,maxmm+1):
                    if not mx in self.dict['MisMatch']: self.dict['MisMatch'][mx] = self.dict['MisMatch'][mx-1]
        except:
            self.log.errorLog('Problem with MisMatch dictionary')
            self.dict['MisMatch'] = {}
#########################################################################################################################
    ### <2> ### General Attribute Methods                                                                               #
#########################################################################################################################
    def motifNum(self): return len(self.list['Motif'])
    def slimNum(self): return len(self.list['Motif'])
    def slims(self): return self.list['Motif']
    def motifs(self): return self.list['Motif']
    def name(self): return self.getStr('Name')
#########################################################################################################################
    def slimCoreName(self,mname): ### Returns core name of motif
        """Returns core name of motif."""
        nsplit = string.split(mname,'_')
        if len(nsplit[-1]) == 1 and nsplit[-1] in string.ascii_letters: mname = string.join(nsplit[:-1],'_')
        return mname
#########################################################################################################################
    def nameList(self,remsplit=False): ### Returns list of motif names
        nlist = []
        for motif in self.slims():
            mname = motif.getStr('Name')
            if remsplit: mname = self.slimCoreName(mname)
            if mname not in nlist: nlist.append(mname)
        nlist.sort()
        return nlist
#########################################################################################################################
    def slimDict(self,corelist=True): ### Returns {corename:[SLiMs]} or {name:slim}
        '''Returns {corename:[SLiMs]} or {name:slim}.'''
        #!# Add self.dict object to store?
        slimdict = {}
        for slim in self.slims():
            name = slim.getStr('Name')
            if corelist:
                if string.split(name,'_')[-1] in string.ascii_letters: name = string.join(string.split(name,'_')[:-1],'_')
                if name in slimdict: slimdict[name].append(slim)
                else: slimdict[name] = [slim]
            else:
                self.warnLog('SLiM "%s" being over-written in slimdict' % name)
                slimdict[name] = [slim]
#########################################################################################################################
    ### <3> ### Motif loading/reformatting methods                                                                      #
#########################################################################################################################
    def run(self): return self.loadMotifs() #!# Add a fuller method with text returned for REST.
#########################################################################################################################
    def loadMotifs(self,motfile=None,clear=True,mlines=[]):      ### Loads motifs and populates self.list['Motif']
        '''
        Loads motifs and populates self.list['Motif'].
        >> file:str = filename or self.info['Motif'] if None
        >> clear:boolean = whether to clear self.list['Motif'] before loading [True]
        >> mlines:list [] = motif file content to over-ride motfile.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear:
                self.list['Motif'] = []
                self.list['Rejects'] = []   # List of motif (cores) rejected to enable rejection of split variants.
            preloadx = self.motifNum()
            if not mlines:
                if not motfile: motfile = self.info['Motifs']
                if motfile.lower() == 'elm': motfile = self.downloadELM()
                ## ~ [1a] Check for split file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                splitsfile = rje.baseFile(motfile) + '.split.motifs'
                if rje.exists(splitsfile) and not self.force() and (self.i() < 1 or rje.yesNo('%s found? Load instead?' % splitsfile)):
                    self.printLog('#SPLIT','Found %s. Will load instead.' % splitsfile)
                    motfile = splitsfile

            ### ~ [2] Read data from file (or direct from list) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if os.path.exists(motfile): mlines = self.loadFromFile(motfile,chomplines=True)
                if mlines: self.info['Name'] = motfile
                elif not os.path.exists(motfile):
                    if motfile.find(',') > 0:
                        mlines = string.split(motfile,',')  # Motif List
                        i = 0
                        while i < (len(mlines)-1):
                            if string.count(mlines[0],'{') == string.count(mlines[0],'}'): i += 1
                            else: mlines[0] += mlines.pop(1)
                    elif motfile.lower() not in ['','none'] and (self.stat['Interactive'] < 1 or rje.yesNo('No lines read from "%s". Use as motif?' % motfile)):
                        #self.warnLog('File "%s" not found!' % motfile)
                        self.log.printLog('#MOT','No lines read from "%s": using as motif.' % motfile)
                        mlines = [motfile]
                    else: return False
                    self.info['Name'] = 'Commandline'

            ### ~ [3] Process motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0  # Number read
            if mlines[0].find('Dataset') == 0:
                try: delimit = mlines[0][7]
                except: delimit = rje.delimitFromExt(filename=motfile)
            else: delimit = rje.delimitFromExt(filename=motfile)
            fhead = rje.readDelimit(mlines[0],delimit)
            if fhead[0].startswith('#ELM_Classes'):
                while mlines[0].startswith('#'): mlines = mlines[1:]
                fhead = rje.readDelimit(mlines[0],delimit)
            #self.deBug(fhead)
            ## ~ [3a] MiniMotif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['MiniMotif']:
                self.reformatMiniMotif(mlines)
                self.log.printLog('#MOT','Motifs read from %s (%d lines): %d retained.' % (motfile,len(mlines),self.motifNum()-preloadx))
                return True
            ## ~ [3b] Fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            if mlines[0][0] == '>':     # Motifs from Fasta file
                mtotx = len(mlines)
                (name,seq) = ('','')
                while mlines or name:
                    if mlines: line = mlines.pop(0)
                    if line[0] == '>':
                        if name:
                            mx += 1
                            self._addMotif(name=name,seq=seq,reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                            if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,mtotx-len(mlines),self.motifNum()-preloadx),log=False,newline=False)
                            name = ''
                        name = line[1:]
                        seq = ''
                    else: seq += line
            ## ~ [3b] TEIRESIAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            elif mlines[0][:2] == '##' and rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',mlines[5]):  # TEIRESIAS
                for line in mlines:
                    motif = rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',line)
                    if motif:
                        mx += 1
                        self._addMotif(name=motif[2],seq=motif[2],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            ## ~ [3c] SLiMDisc rank file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            elif mlines[0][:2] == '#-' and rje.matchExp('Input dataset:\s+(\S+)',mlines[1]):  # SLiMDisc
                name_root = os.path.splitext(os.path.basename(rje.matchExp('Input dataset:\s+(\S+)',mlines[1])[0]))[0]
                for line in mlines:
                    motif = rje.matchExp('^\((\d+)\)\s+\S+\s+(\S+)',line)
                    if motif:
                        mx += 1
                        self._addMotif(name='%s#%s' % (name_root,motif[0]),seq=motif[1],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            ## ~ [3d] Delimited file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif 'Pattern' in fhead:
                if 'Name' in fhead: namekeys = ['Name']
                else:
                    namekeys = []
                    for key in ['Dataset','RunID','Rank']:
                        if key in fhead: namekeys.append(key)
                namekeys += ['Pattern']
                sfdata = rje.dataDict(self,motfile,mainkeys=namekeys,datakeys=['Pattern'],delimit=delimit)
                for returned in rje.sortKeys(sfdata):
                    data = sfdata[returned]
                    if data['Pattern'] and data['Pattern'] not in ['.','-','X','!']:
                        name = string.replace(string.join(string.split(returned,delimit),'#'),' ','_')
                        mx += 1
                        self._addMotif(name=name,seq=data['Pattern'],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            ## ~ [3e] ELM classes download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif 'Regex' in fhead:
                sfdata = rje.dataDict(self,motfile,mainkeys=['ELMIdentifier'],datakeys=['ELMIdentifier','Regex','Description'],delimit=delimit,ignore=['#'])
                for returned in rje.sortKeys(sfdata):
                    data = sfdata[returned]
                    if data['ELMIdentifier'] == 'ELMIdentifier': continue
                    name = '%s %s' % (data['ELMIdentifier'],data['Description'])
                    #self.deBug('%s -> %s' %(name,data['Regex']))
                    mx += 1
                    self._addMotif(name=name,seq=data['Regex'],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                    if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            ## ~ [3f] PRESTO format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:                
                for line in mlines:
                    if not line: continue   # Skip blank lines
                    mx += len(self.motifsFromPRESTOLine(line))
                    if self.stat['Verbose'] < 2: self.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            ## ~ [3x] Motif read error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.slims() and not os.path.exists(motfile):
                self.printLog('#ERR','File %s not found and not recognised as motif(s). Check for typos.' % motfile)
                raise IOError('%s not found' % motfile)

            ### ~ [4] Additional filtering of split motif variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rx = 0
            for slim in self.slims()[0:]:
                if self.slimCoreName(slim.name()) in self.list['Rejects']:
                    self.removeMotif(slim,'%s (%s) removed: %s variant rejected.  ' % (slim.name(),slim.pattern(),self.slimCoreName(slim.name())))
                    rx += 1
            if rx: self.printLog('\r#REJ','%s motifs rejected following filtering of split motif variants.' % rje.iStr(rx))

            ### ~ [5] Summarise and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sumtxt = '%d motifs read from %s (%d lines): %d retained.' % (mx,motfile,len(mlines),self.motifNum()-preloadx)
            if self.opt['Reverse']: sumtxt += ' (Reversed)'
            if self.opt['WildScram']: sumtxt += ' (WildScram)'
            self.printLog('\r#MOT',sumtxt)
            self.motifOut(self.info['MotifOut'])
            if self.getStrLC('MotInfo') not in ['','none']: self.motifInfo()
            return True
        except: self.log.errorLog('Error in loadMotifs()')     
        return False
#########################################################################################################################
    def motifsFromPRESTOLine(self,line):    ### Returns a list of motif objects from one PRESTO-format line
        '''Returns a list of motif objects from one PRESTO-format line.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not line: return []
            prex = self.motifNum()
            ### ~ [2] Parse Motif(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            desc = string.join(string.split(line,'#')[1:],'#')
            while desc[:1] == ' ': desc = desc[1:]
            motif = rje.matchExp('^(\S+)\s+(\S+)',string.split(line,'#')[0])
            if motif:
                self._addMotif(name='%s %s' % (motif[0],desc),seq=motif[1],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
            elif rje.matchExp('^(\S+)',string.split(line,'#')[0]):   # Pure patterns: name = pattern
                motif = rje.matchExp('^(\S+)',string.split(line,'#')[0])[0]
                self._addMotif(name='%s %s' % (motif,desc),seq=motif,reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
            return self.motifs()[prex:]        
        except: self.log.errorLog('Error in motifsFromPRESTOLine()')     
#########################################################################################################################
    def reformatMiniMotif(self,mlines):   ### Reformats MiniMotif file, compressing motifs as appropriate
        '''
        Reformats MiniMotif file, compressing motifs as appropriate.
        >> mlines:list of lines read from input file
        '''
        try:### ~ [1] Make dictionaries of {name:[patterns]} and {name:desc} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            patdict = {}    # Dictionary of name:patterns
            desdict = {}    # Dictionary of name:description
            adddict = {}    # Dictionary of name and extra letters to add
            for m in mlines:
                ## ~ [1a] Get details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                details = rje.matchExp('^(\S+)\s+(\S+)\s+#\s*(.+)',m)
                if details:
                    (name,pattern,desc) = details
                    pattern = rje_motif.reformatMiniMotif(self,pattern)
                ## ~ [1b] Update dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for mname in patdict.keys():
                    if mname[:1] == name[:1] and desdict[mname] == desc:    # Add to this motif
                        patdict[mname].append(pattern)
                        if name[-1] not in adddict[mname]: adddict[mname].append(name[-1])
                        name = ''
                        break
                if name:
                    patdict[name] = [pattern]
                    adddict[name] = [name[-1]]
                    desdict[name] = desc

            ### ~ [2] Compress MiniMotifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for name in rje.sortKeys(patdict):
                newpat = rje_motif.defineMotif(self,patdict[name],profile=False,minfreq=0.05,minocc=1,ambcut=19)
                add = adddict[name]
                add.sort()
                mname = name[:-1] + string.join(add,'')
                self.log.printLog('#MOT','Motif "%s" => %s' % (mname,newpat))
                if len(newpat) > 1:     ### Make several
                    ext = 'abcdefghijklmnopqrstuvwxyz'
                    for i in range(len(newpat)): self._addMotif(name='%s%s %s' % (mname,ext[i],desdict[name]),seq=newpat[i],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
                else: self._addMotif(name='%s %s' % (mname,desdict[name]),seq=newpat[0],reverse=self.opt['Reverse'],check=True,wildscram=self.opt['WildScram'])
        except:
            self.log.errorLog('Error in _reformatMiniMotif()')     
            raise   
#########################################################################################################################
    def _reject(self,corename): self.list['Rejects'].append(corename); return None
#########################################################################################################################
    def _addMotif(self,name,seq,reverse=False,check=False,logrem=True,wildscram=False):  ### Adds new motif to self.list['Motif'].
        '''
        Adds new motif to self.list['Motif']. Checks redundancy etc.
        >> name:str = Motif Name
        >> seq:str = Motif Sequence read from file
        >> reverse:boolean [False] = whether to reverse sequence
        >> check:boolean [False] = whether to check redundancy and sequence length
        >> logrem:boolean [True] = whether to log removal of motifs
        >> wildscram:boolean [False] = whether to perform wildcard scrambling on motif
        << returns Motif object or None if failed
        '''
        try:### ~ [1] Setup name and description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            desc = string.join(string.split(name)[1:])
            name = string.split(name)[0]
            if string.split(name,'_')[-1] in string.ascii_letters: corename = string.join(string.split(name,'_')[:-1],'_')
            else: corename = name
            if not rje.matchExp('([A-Za-z])',seq):
                self.progLog('\r#REM',' ' * 100)
                self.printLog('\r#REM','Motif %s ("%s") has no discernable positions' % (name,seq))
                return self._reject(corename)
            ## ~ [1a] Check Name for GoodMotif filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['GoodMotif']:
                dump = True
                for good in self.list['GoodMotif']:
                    if rje.matchExp('^(%s)$' % string.replace(good,'*','\S*'),name) or rje.matchExp('^(%s)$' % string.replace(good,'*','\S*'),corename):
                        dump = False
                        break
                if dump:
                    self.progLog('\r#REM',' ' * 100)
                    self.log.printLog('\r#REM','Motif "%s" not in goodmotif=LIST.' % (name))
                    return self._reject(corename)
            ## ~ [1b] Reverse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if reverse: name = '%s_rev' % name
            if wildscram: name = '%s_scram' % name
            ## ~ [1c] Motif splitting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seq = string.replace(seq,'?','{0,1}')
            #self.bugPrint('%s:%d' % (seq,seq.find('|')))
            seq = string.replace(seq,'?','{0,1}')
            if check and not self.getBool('VarLength') and '{' in seq:
                self.progLog('\r#REM',' ' * 100)
                self.log.printLog('\r#REM','Motif %s (%s) has variable length positions (varlength=F).' % (name, seq),screen=logrem,log=logrem)
                return self._reject(corename)
            if rje_slim.needToSplitPattern(seq): return self._splitMotif(name,seq,desc)
            #if seq.find('|') > 1: return self._splitMotif(name,seq,desc)

            ### ~ [2] Setup SLiM object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newmotif = rje_slim.SLiM(log=self.log,cmd_list=self.cmd_list)
            newmotif.info['Name'] = name
            if desc: newmotif.info['Description'] = desc
            else: newmotif.info['Description'] = name
            newmotif.info['Sequence'] = seq
            newmotif.obj['SLiMList'] = self
            newmotif.opt['TrimX'] = self.opt['TrimX']
            if not newmotif.format(reverse=reverse):
                self.progLog('\r#REM',' ' * 100)
                self.printLog('\r#REJECT','Formatting problem with %s sequence %s: Rejected.' % (name,seq))
                return self._reject(corename)
            if wildscram: newmotif.wildScram()

            ### ~ [3] Check Length filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if check and not self.getBool('VarLength') and newmotif.varLength():
                self.progLog('\r#REM',' ' * 100)
                self.log.printLog('\r#REM','Motif %s (%s) has variable length positions (varlength=F).' % (newmotif.info['Name'], newmotif.info['Sequence']),screen=logrem,log=logrem)
                return self._reject(corename)
            if check and newmotif.slimFix() < self.stat['MinFix']:
                self.progLog('\r#REM',' ' * 100)
                self.log.printLog('\r#REM','Motif %s (%s) has < %d fixed positions (%d).' % (newmotif.info['Name'], newmotif.info['Sequence'],self.stat['MinFix'],newmotif.slimFix()),screen=logrem,log=logrem)
                return self._reject(corename)
            elif check and newmotif.slimPos() < self.stat['MinPos']:
                self.progLog('\r#REM',' ' * 100)
                self.log.printLog('\r#REM','Motif %s (%s) has < %d defined positions (%d).' % (newmotif.info['Name'], newmotif.info['Sequence'],self.stat['MinFix'],newmotif.slimPos()),screen=logrem,log=logrem)
                return self._reject(corename)
            elif check and newmotif.stat['IC'] < self.stat['MinIC']:
                self.progLog('\r#REM',' ' * 100)
                self.log.printLog('\r#REM','Motif %s (%s) has insufficient IC (%.1f < %.1f).' % (newmotif.info['Name'], newmotif.info['Sequence'],newmotif.stat['IC'],self.stat['MinIC']),screen=logrem,log=logrem)
                return self._reject(corename)
            else: self.verbose(2,4,'%s: %s (%d non-X aa)' % (newmotif.info['Name'], newmotif.info['Sequence'],newmotif.slimPos()),1)

            ### ~ [4] Check Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if check and self.opt['NRMotif']:
                for motif in self.slims():
                    if newmotif.info['Slim'] == motif.info['Slim']:
                        self.progLog('\r#REM',' ' * 100)
                        self.log.printLog('\r#REM','Motif %s is identical to %s and has been removed.' % (newmotif.info['Name'],motif.info['Name']),screen=logrem,log=logrem)
                        return self._reject(corename)
                    elif motif.info['Slim'].find(newmotif.info['Slim']) >= 0:
                        self.progLog('\r#REM',' ' * 100)
                        self.log.printLog('\r#REM','Motif %s is a subsequence of %s and has been removed.' % (newmotif.info['Name'],motif.info['Name']),screen=logrem,log=logrem)
                        return self._reject(corename)

            ### ~ [5] Complete Formatting and Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newmotif.misMatch(mismatch=self.dict['MisMatch'])
            self.list['Motif'].append(newmotif)
            return newmotif
        except:
            self.log.errorLog('Error in MotifList._addMotif()',quitchoice=True)     
            return None
#########################################################################################################################
    def _splitMotif(self,name,seq,desc):     ### Splits complex motifs on "|" and adds each separately
        '''Splits complex motifs on "|" and adds each separately.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bases = []          # This is a list of basic strings to be added to with each variation
            newmotifs = [seq]   # New motif variants
            splitting = True

            ### ~ [1] Cycle and split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while splitting:
                bases = newmotifs[0:]
                newmotifs = []
                splitting = False
                for base in bases:
                    #self.bugPrint('%s >> %s' % (bases,newmotifs))
                    i = base.find('|')
                    if i < 0:
                        vardef = rje.matchExp('([\[A-Za-z\^\]]+)\{(\d+),(\d+)\}',base)
                        if not vardef: newmotifs.append(base); continue
                        splitting = True
                        (nonwild,m,n) = vardef
                        m = int(m)
                        n = int(n)
                        if nonwild[-1] == ']': nonwild = '[%s' % string.split(nonwild,'[')[-1]
                        else: nonwild = nonwild[-1]
                        #self.deBug('%s -> %s{%d,%d}' % (base,nonwild,m,n))
                        if m == n: newmotifs.append(string.replace(base,'%s{%d,%d}' % (nonwild,m,n),nonwild*m))
                        else:
                            for x in range(m,n+1): newmotifs.append(string.replace(base,'%s{%d,%d}' % (nonwild,m,n),nonwild*x))
                        continue
                    splitting = True
                    pre = 0 # Number of pre parentheses - check X|Y split has surrounding parentheses (or add)
                    x = i   # x = start of selection to vary
                    while pre < 1:
                        x -= 1
                        if x < 0:   # The whole pattern is an X|Y split with no outside brackets.
                            base = '(%s)' % base
                            i += 1
                            #self.bugPrint('%s -> %s' % (base[1:-1],base))
                            break
                        if base[x] == '(': pre += 1
                        if base[x] == ')': pre -= 1
                    ## ~ [1a] Find selection spot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    pre = post = 0  # Number of pre/post parentheses
                    x = y = i       # x & y = ends of selection to vary
                    #self.deBug(i)
                    while pre < 1:
                        x -= 1
                        if x < 0:   # The whole pattern is an X|Y split with no outside brackets.
                            raise ValueError
                        if base[x] == '(': pre += 1
                        if base[x] == ')': pre -= 1
                        #self.bugPrint('%s:%s:%s - %s>>%s' % (x,i,y,base[x],base[y]))
                    while post < 1:
                        y += 1
                        if y >= len(base): raise ValueError
                        if base[y] == '(': post -= 1
                        if base[y] == ')': post += 1
                        #self.bugPrint('%s:%s:%s - %s>>%s' % (x,i,y,base[x],base[y]))
                    ## ~ [1b] Replace swap region with variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #self.bugPrint('%s >> %s & %s' % (base,base[x+1:i],base[i+1:y]))
                    #self.bugPrint('%s:%s:%s - %s>>%s' % (x,i,y,base[x],base[y]))
                    for var in (base[x+1:i],base[i+1:y])[0:]:
                        if var[0] == '(' and var[-1] == ')': var = var[1:-1]
                        newmotifs.append(base[:x] + var + base[y+1:])
                #self.deBug('%s >> %s' % (bases,newmotifs))

            ### ~ [2] Clean up motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newmotifs = rje.sortUnique(newmotifs)
            if len(newmotifs) == 1: return self._addMotif(name,newmotifs[0])
            elif len(newmotifs) > len(string.ascii_letters):
                self.printLog('#SPLIT','Too many variants (%d) for motif split: redefine %s!' % (len(newmotifs),name))
                raise ValueError
            if self.info['MotifOut'].lower() in ['','none']: self.info['MotifOut'] = rje.baseFile(self.info['Name']) + '.split.motifs'
            self.printLog('\r#SPLIT','Input motif "%s" - %s - split into %d variants' % (name,seq,len(newmotifs)))
            ## ~ [2a] Add splits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            splitmotifs = []
            for m in range(len(newmotifs)):
                newseq = newmotifs[m]
                newname = '%s_%s %s variant %d/%d' % (name,string.ascii_letters[m],desc,m+1,len(newmotifs))
                newmotif = self._addMotif(newname,newseq)
                if newmotif: splitmotifs.append(newmotif)
            ## ~ [2b] Check for cleaned up variants: remove all if any part is rejected? ~~~~~~~~~~ ##
            self.printLog('#SPLIT','%d of %d %s variants added (not filtered)' % (len(splitmotifs),len(newmotifs),name))
            if len(splitmotifs) != len(newmotifs):
                for newmotif in splitmotifs: self.removeMotif(newmotif)
                self.printLog('#REM','Removed %d %s variants: other variants filtered.' % (len(splitmotifs),name))

        except ValueError:
            self.printLog('#REJECT','Formatting problem with %s sequence %s: Rejected.' % (name,seq))
            return None
        except: self.errorLog('MotifList._splitMotif Error')
#########################################################################################################################
    def removeMotif(self,Motif,remtxt=''):    ### Removes motif and occurrences from self.
        '''
        Removes motif and occurrences from self.
        >> Motif:Motif object to remove
        >> remtxt:str = Text to output to log
        '''
        try:### ~ [1] Remove Motif from list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if Motif in self.list['Motif']:
                self.list['Motif'].remove(Motif)
                if remtxt: self.log.printLog('\r#REM',remtxt)
            elif Motif: self.log.errorLog('Cannot remove SLiM "%s" - not in SLiMList!' % Motif.info['Name'],printerror=False)
            else: return self.log.errorLog('No Motif to remove')
        except: self.log.errorLog('Error in SLiMList.removeMotif()')     
#########################################################################################################################
    def downloadELM(self):  ### Downloads ELM from website (if missing) and returns processed ELM classes file name
        '''Downloads ELM from website (if missing) and returns processed ELM classes file name.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            downloader = rje_obj.RJE_Object(self.log,self.cmd_list)
            db = downloader.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            downloader.dict['SourceURL'] = {'ELMClass':'http://elm.eu.org/elms/browse_elms.html?q=&submit=tsv',
                                            'ELMInstance':'http://elm.eu.org/elms/browse_instances.tsv?q=*&taxon=&instance_logic=',
                                            'ELMInteractors':'http://elm.eu.org/interactions/as_tsv'}
            downloader.setStr({'SourcePath':'./',
                               'ELMClass':'elm_classes.tsv',
                               'ELMInstance':'elm_instances.tsv',
                               'ELMInteractors':'elm_interaction_domains.tsv'})

            ### ~ [2] Generate ELM motif file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elmc = db.addTable(downloader.sourceDataFile('ELMClass',download=True),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
            elmc.dataFormat({'#Instances':'int'})
            elmbase = '%s' % rje.baseFile(downloader.getStr('ELMClass'),strip_path=True)
            motif_file = '%s.motifs' % elmbase
            if rje.checkForFile(motif_file) and not self.force():
                self.printLog('#ELM','%s file found (force=F).' % motif_file)
            else:
                rje.backup(self,motif_file)
                motif_out = []
                for elm in elmc.dataKeys():
                    entry = elmc.data(elm)
                    try: motif_out.append('%s  %s  # %s [%d ELM instances]' % (entry['ELMIdentifier'],entry['Regex'],entry['Description'],entry['#Instances']))
                    except: self.debug(entry); raise
                open(motif_file,'w').write(string.join(motif_out,'\n'))
                self.printLog('#ELM','%s motif patterns output to %s' % (elmc.entryNum(),motif_file))

            ### ~ [3] Optional download of other ELM data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.i() < 0 or rje.yesNo('Check/Download ELM instance and interactor data?'):
                for etype in ['ELMInstance','ELMInteractors']:
                    self.printLog('#ELM','%s: %s' % (etype,downloader.sourceDataFile(etype,download=True)))

            return motif_file
        except: self.errorLog('Error in SLiMList.downloadELM()')
#########################################################################################################################
    ### <4> ### Motif occurrence methods                                                                                #
#########################################################################################################################
    def motifOcc(self,byseq=False,justdata=None,fastacmd=False,nested=True,maxocc=0):    ### Returns MotifOcc list/dictionary 
        '''
        Returns a list or dictionary of MotifOccurrences. Partially converts occurrences to old rje_motifocc style with
        Seq and Motif as part of occurrence.
        >> byseq:bool [False] = return a dictionary of {Sequence:{Motif:OccList}}, else {Motif:{Sequence:OccList}} 
        >> justdata:str [None] = if given a value, will return this data entry for each occ rather than the object itself
        >> fastacmd:bool [False] = whether to return FastaCmd instead of Sequence if Sequence missing
        >> nested:bool [True] = whether to return a nested dictionary or just the occurrences per Motif/Seq
        << returns dictionary or plain list of MotifOcc if byseq and bymotif are both False        
        '''
        try:### ~ [1] Setup return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motifocc = {}

            ### ~ [2] Populate dictionary according to settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for Motif in self.slims():
                if maxocc > 0 and Motif.occNum() > maxocc: continue
                for Seq in Motif.dict['Occ']:
                    for occ in Motif.dict['Occ'][Seq]:
                        occ['Seq'] = Seq
                        occ['Motif'] = Motif
                        if fastacmd and not Seq: Seq = occ['FastaCmd']
                        if byseq: (key1,key2) = (Seq,Motif)
                        else: (key2,key1) = (Seq,Motif)
                        if not motifocc.has_key(key1):
                            if nested: motifocc[key1] = {}
                            else: motifocc[key1] = []
                        if nested and not motifocc[key1].has_key(key2): motifocc[key1][key2] = []
                        ## Update ##
                        if justdata:
                            try: val = occ[justdata]
                            except: continue    # Skip this occurrence as desired data is missing
                        else: val = occ
                        if nested: motifocc[key1][key2].append(val)
                        else: motifocc[key1].append(val)
            return motifocc
        except:
            self.log.errorLog('Problem with MotifList.motifOcc()')
            raise
#########################################################################################################################
    ### <5> ### Motif attribute methods                                                                                 #
#########################################################################################################################
    def calculateOccAttributes(self,silent=False,wallobj=None):   ### Executes rje_slimcalc calculations via rje_slimlist object
        '''Executes rje_slimcalc calculations via rje_slimlist object.'''
        try:### ~ [1] Call rje_slimcalc for the occurrences on each sequence in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if wallobj:
                slimocc = self.motifOcc(byseq=True,nested=False,maxocc=wallobj.getInt('MaxOcc'))
                if wallobj.getInt('MaxOcc') > 0:
                    for Motif in self.slims():
                        mx = Motif.occNum()
                        if mx > wallobj.getInt('MaxOcc'):
                            self.printLog('#MAX','Maximum occurrence threshold (%d) exceeded for Motif "%s" (%d occ)' % (wallobj.getInt('MaxOcc'),Motif.info['Name'],mx))
            else: slimocc = self.motifOcc(byseq=True,nested=False)
            (sx,stot) = (0.0,len(slimocc))
            for seq in slimocc:
                if not ('cons' in self.obj['SLiMCalc'].list['SLiMCalc'] or silent):
                    self.log.printLog('\r#ATT','Calculating occurrence attributes: %.1f%%' % (sx/stot),newline=False,log=False)
                    sx += 100.0
                self.obj['SLiMCalc'].occStats(slimocc[seq],xpad=0,progress=False,silent=silent)
                if wallobj: wallobj.wallTime()
            if not ('cons' in self.obj['SLiMCalc'].list['SLiMCalc'] or silent):
                self.log.printLog('\r#ATT','Calculating occurrence attributes complete.',log=False)
        except: self.log.errorLog('Problem during rje_slimlist.calculateOccAttributes')
#########################################################################################################################
    def combMotifOccStats(self,revlist=['Hyd'],progress=True,motiflist=[]):    ### Combines mean and percentile stats for the Occurrences of a Motif
        '''
        Combines mean and percentile stats for the Occurrences of each Motif. Combines all attributes stored in
        self.obj['SLiMCalc'].list['Headers']
        >> statlist:list of stats to combine from occurrences
        >> revlist:list of stats that should be ordered from low(best) to high(worst) rather than the other way round
        >> progress:boolean [True] = whether to log progress of stat combination
        >> motiflist:list of motifs to combine (will use all if none given)
        '''
        try:### ~ [1] Setup Motif OccList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motifocc = self.motifOcc(byseq=False,nested=False)
            for Motif in motifocc.keys():
                if motiflist and Motif not in motiflist: motifocc.pop(Motif)
            
            ### ~ [2] Combine attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
            (mx,mtot) = (0.0,len(motifocc))
            for Motif in motifocc:
                if progress: self.log.printLog('\r#COMB','Combined Occurrence Attributes: %.1f%%' % (mx/mtot),log=False,newline=False)
                mx += 100.0
                try: self.obj['SLiMCalc'].combMotifOccStats(motifocc[Motif],revlist)
                except: self.log.errorLog('Major problem with rje_slimCalc.combMotifOccStats(%s)' % Motif.info['Name'])
            if progress: self.log.printLog('\r#COMB','Combined Occurrence Attributes complete.')
        except:
            self.log.errorLog('Major problem with rje_motiflist.combMotifOccStats()')
            raise
#########################################################################################################################
    def setupFilters(self,slimheaders=[],occheaders=[]): self.obj['SLiMCalc'].setupFilters(slimheaders,occheaders)
#########################################################################################################################
    ### <6> ### Motif output methods                                                                                    #
#########################################################################################################################
    def motifOut(self,filename='None',motlist=[],pattern=False):      ### Outputs motifs in SLiMSuite format
        '''
        Outputs motifs in SLiMSuite format.
        >> filename:str [None] = Name for output file. Will not output if '' or 'None' - returns as text instead.
        >> motlist:list of Motif objects to output. If [], will use self.list['Motifs']
        >> pattern:bool [False] = Whether to only output/return motif patterns (no names/descriptions).
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not motlist: motlist = self.slims()
            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename and filename.lower() != 'none':
                rje.backup(self,filename)
                OUT = open(filename,'a')
                for Motif in motlist:
                    if pattern: OUT.write('%s\n' % (Motif.info['Sequence']))
                    else: OUT.write('%s  %s  # %s\n' % (Motif.info['Name'],Motif.info['Sequence'],Motif.info['Description']))
                OUT.close()
                self.printLog('#OUT','%s motifs output to %s in SLiMSuite format.' % (rje.iLen(motlist),filename))
            else:
                motiftxt = ''
                for Motif in motlist:
                    if pattern: motiftxt += '%s\n' % (Motif.info['Sequence'])
                    else: motiftxt += '%s  %s  # %s\n' % (Motif.info['Name'],Motif.info['Sequence'],Motif.info['Description'])
                return motiftxt
        except: self.log.errorLog('Problem during SLiMList.motifOut()')
#########################################################################################################################
    def proteinAlignments(self,alndir='',hitname='AccNum'):    ### Generates copies of protein alignments, with motif hits marked.
        '''Generates copies of protein alignments, with motif hits marked.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alndir and not os.path.exists(alndir): rje.mkDir(self,alndir)
            ### ~ [2] Generate Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motif_occ = self.motifOcc(byseq=True,nested=False)
            for seq in motif_occ.keys():
                self.obj['SLiMCalc'].singleProteinAlignment(seq,motif_occ[seq],alndir,hitname,usegopher=True)                
        except:            
            self.log.errorLog('Major problem with SlimPicker.proteinAlignments',quitchoice=True)
#########################################################################################################################
    def motifAlignments(self,resfile='motifaln.fas'):      ### Makes motif alignments from occurrences
        '''
        Makes motif alignments from occurrences. MotifOcc objects should have Sequence objects associated with them. If
        necessary, add a method to go through and generate Sequence objects using MotifOcc.info['FastaCmd'].
        >> resfile:str = Name of output file
        '''
        try:### ~ [0] New long output for easier conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['XDivide'] < 1:
                motif_occ = self.motifOcc()
                for Motif in self.slims(): self.motifAlnLong(Motif,motif_occ[Motif],append=True,memsaver=False,resfile=resfile)
                return
            
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motif_occ = self.motifOcc(byseq=True,justdata='Pos')
            for Seq in motif_occ.keys():
                for Motif in motif_occ[Seq].keys():
                    newpos = []
                    for pos in motif_occ[Seq][Motif]: newpos.append(int(pos) - 1)
                    motif_occ[Seq][Motif] = newpos
            ltxt = 'Constructing motif alignments: %s motifs & %s seqs.' % (rje.integerString(self.motifNum()),rje.integerString(len(motif_occ)))
            self.log.printLog('#ALN',ltxt,log=False,newline=False)
            extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
            extract_seq.opt['ReplaceChar'] = False
            ## ~ [1a] Setup positions in sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Work off motif_occ = Dictionary of {Sequence:{Motif:[positions]}} (Pos is now from 0 to L-1)
            max_occ = {}            # Dictionary of Sequence:Max occurrences for any motif
            max_pos = 0             # Max motif position
            occ_seq = {'Motif':{}}  # Dictionary of {Sequence:{Motif:[aligned sequence fragments]}}
            motiflist = []          # List of motifs with 1+ occurrences
            for motif in self.slims(): occ_seq['Motif'][motif] = []
            for seq in motif_occ.keys():
                occ_seq[seq] = {}
                max_occ[seq] = 0
                for motif in self.slims():
                    occ_seq[seq][motif] = []
                    if motif not in motif_occ[seq].keys(): continue    #!# Improve way missing sequences/motifs are dealt with #!#
                    if len(motif_occ[seq][motif]) > max_occ[seq]: max_occ[seq] = len(motif_occ[seq][motif])
                    if len(motif_occ[seq][motif]) > 0 and motif not in motiflist: motiflist.append(motif)
                    for pos in motif_occ[seq][motif]:
                        if pos >= max_pos: max_pos = pos + 1
            self.log.printLog('\r#ALN','%s: setup complete.' % ltxt,log=False)
            neworder = self.slims()[0:]
            for Motif in self.slims()[0:]:
                if Motif not in motiflist: neworder.remove(Motif)
            motiflist = neworder[0:]

            ### ~ [2] Make sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for motif in motiflist:
                pattern = motif.info['Sequence']
                self.log.printLog('\r#ALN','%s: %.f%%.' % (ltxt,100.0*motiflist.index(motif)/len(motiflist)),log=False,newline=False)
                ## Pattern info ##
                patseq = '-%s-' % string.split(motif.info['Name'])[0]   # Name of motif
                while len(patseq) < (len(rje.preZero(max_pos,max_pos)) + 2): patseq += '-'
                namelen = len(patseq)
                overlap = len(pattern) - motif.slimLen()    # Extra length of pattern vs. longest occurrence
                patseq += '-' * self.stat['FlankSize']
                patseq += pattern       
                patseq += '-' * (self.stat['FlankSize'] - overlap)
                occ_seq['Motif'][motif].append(patseq[0:])
                ## Occurrences ##
                for seq in motif_occ.keys():
                    if not motif_occ[seq].has_key(motif): motif_occ[seq][motif] = {}
                    for x in range(max_occ[seq]):  # Need an entry for each potential occurrence
                        if x < len(motif_occ[seq][motif]):  # Actual entry
                            r = motif_occ[seq][motif][x]
                            patseq = '-%s-' % rje.preZero(r+1,max_pos)
                            while len(patseq) < namelen: patseq += '-'
                            (left,right) = (r-self.stat['FlankSize'],r+motif.slimLen()-1+self.stat['FlankSize'])
                            if left < 0: patseq += '-' * -left + seq.info['Sequence'][:right]
                            else: patseq += seq.info['Sequence'][left:right]
                            patseq += '-' * (len(occ_seq['Motif'][motif][0]) - len(patseq))
                        else: patseq = '-' * len(occ_seq['Motif'][motif][0])    # Add blank one                            
                        occ_seq[seq][motif].append(patseq[0:])
            self.log.printLog('\r#ALN','%s: 100.0%%.' % (ltxt))

            ### ~ [3] Output sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Motif Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            name = '%s motifs' % self.info['Name']
            outlist = []
            for motif in motiflist:
                pattern = motif.info['Sequence']
                outlist.append(occ_seq['Motif'][motif][0])
            extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ## ~ [3b] Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for seq in motif_occ.keys():
                name = seq.info['Name']
                for x in range(max_occ[seq]):
                    outlist = []
                    for motif in motiflist:
                        pattern = motif.info['Sequence']
                        outlist.append(occ_seq[seq][motif][x])
                    extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ## ~ [3c] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            extract_seq.info['Name'] = resfile
            extract_seq.saveFasta()

        except: self.log.errorLog('Major problem with SLiMList.motifAlignments()',quitchoice=True)
#########################################################################################################################
    def motifAlnLong(self,Motif,seq_occ,append=False,memsaver=False,resfile=''):     ### Makes a single MotifAln output
        '''
        Makes a single MotifAln output.
        >> Motif:Motif Object
        >> seq_occ:dictionary of {Seq:[MotifOcc]}
        >> append:bool [False] = whether to append file or create new
        >> memsaver:bool [False] = whether output is to be treated as a single sequence or all occurrences
        >> resfile:str [''] = name for motif alignment output file
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not resfile: resfile = Motif.info['Name'] + '.motifaln.fas'
            extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
            extract_seq.opt['ReplaceChar'] = False
            ## ~ [1a] Motif Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            max_pos = 10000
            patseq = '-' * (len(rje.preZero(max_pos,max_pos)) + 2)
            overlap = len(Motif.info['Sequence']) - Motif.stat['FullLength']       # Extra length of pattern vs. longest occurrence
            patseq += '-' * self.stat['FlankSize']
            patseq += Motif.info['Sequence']
            patseq += '-' * (self.stat['FlankSize'] - overlap)
            if memsaver or not append: extract_seq._addSeq(Motif.info['Name'],patseq)      ### Output Motif itself
            motlen = len(patseq)
                
            ### ~ [2] Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for Seq in seq_occ:
                for Occ in seq_occ[Seq]:
                    if Occ.obj['Motif'] != Motif: continue
                    patseq = '-%s-' % rje.preZero(Occ.stat['Pos'],max_pos)
                    r = Occ.stat['Pos'] - 1
                    (left,right) = (r-self.stat['FlankSize'],r+Motif.stat['FullLength']-1+self.stat['FlankSize'])
                    if left < 0: patseq += '-' * -left + Seq.info['Sequence'][:right]
                    else: patseq += Seq.info['Sequence'][left:right]
                    patseq += '-' * (motlen - len(patseq))
                    extract_seq._addSeq(Seq.shortName(),patseq)

            ### ~ [3] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            extract_seq.saveFasta(seqfile=resfile,append=append,log=False)
        except: self.log.errorLog('Major problem with MotifList.MotifAlnLong',quitchoice=True)            
#########################################################################################################################
### End of SECTION II: SLiMList Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SLiMList methods to add                                                                                #
#########################################################################################################################
    def checkForOcc(self,Occ,merge=True):   ### Returns existing MotifOcc if one is found, else given Occ
        '''
        Returns existing MotifOcc if one is found, else given Occ.
        >> Occ:MotifOcc object to check
        >> merge:bool [True] = whether to merge Occ with found occurrence, if there is one
        '''
        try:
            ### Setup ###
            Seq = Occ.obj['Seq']
            Motif = Occ.obj['Motif']    #!# Given Occ should already have mapped Motif #!#
            if not Occ.getData('Hit') and not seq:
                return Occ
            if not Motif:
                Motif = self.mapMotif(Occ,update=True)
                Occ.obj['Motif'] = Motif
            ### Check for Motif ###
            mymotoccs = self.motifOcc()
            if not mymotoccs.has_key(Motif):     ## No occs for this motif!
                return Occ
            ### MapSeq and generate possible occlist ###
            occlist = []
            if Seq and mymotoccs[Motif].has_key(Seq):
                occlist = mymotoccs[Motif][Seq]
            else:
                for key in mymotoccs[Motif].keys():
                    if key:     # Seq object
                        if Occ.getData('Hit') in [key.shortName(),key.info['AccNum'],key.info['Name']]:
                            occlist = mymotoccs[Motif][key]
                            Occ.obj['Seq'] = key
                            break
                        continue
                    else:
                        for myocc in mymotoccs[Motif][key]:
                            if myocc.getData('Hit') and myocc.getData('Hit') == Occ.getData('Hit'):
                                occlist.append(myocc)
            ### CheckForOcc ###
            for myocc in occlist:
                if myocc.getData('Pos') and myocc.getData('Pos') == Occ.getData('Pos') and myocc.getData('Variant') == Occ.getData('Variant'):
                    if merge:
                        self.mergeOcc([myocc,Occ])
                    return myocc
            return Occ
        except:
            self.log.errorLog('Problem with MotifList.checkForOcc()',quitchoice=True)
            return Occ
#########################################################################################################################
    def mapMotif(self,Occ,update=True):     ### Returns Motif Object based on self.slims() and Occ data
        '''
        Returns Motif Object and updates Occ, based on self.slims() and Occ data.
        >> Occ:MotifOcc object to check
        >> update:bool [True] = whether to update own list['Motifs'] and/or Motif objects
        '''
        try:
            ### Check for exact Motif ###
            if Occ.obj['Motif'] in self.slims():
                return Occ.obj['Motif']
            ### Check for Motif Name with compatible varlist ###
            for Motif in self.slims():
                if Motif.info['Name'] == Occ.getData('Motif') and Occ.getData('Variant') in Motif.list['Variants']:
                    return Motif
                elif Motif.info['Name'] == Occ.getData('Motif') and update:
                    varlist = Motif.list['Variants'] + [Occ.info['Variant']]
                    newdef = rje_motif.defineMotif(self,occlist=varlist,minfreq=0,minocc=1)
                    if len(newdef) == 1:    # OK to integrate
                        Motif.info['Sequence'] = newdef[0]
                        Motif.format(msmode=self.opt['MSMS'])
                        Motif.makeVariants(msmode=self.opt['MSMS'],ambvar=not self.opt['Compare'])   
                        Motif.misMatches(mismatch=self.dict['MisMatch'],msmode=self.opt['MSMS'])
                        return Motif
            ### Check for Motif Name matching variant ###
            for Motif in self.slims():
                if Motif.info['Name'] == Occ.getData('Variant'):
                    return Motif
            ### Check varlist only ###
            if not Occ.getData('Motif') and Occ.getData('Variant'):
                for Motif in self.slims():
                    if Occ.getData('Variant') in Motif.list['Variants']:
                        return Motif
            ### Make new Motif ###
            if update:
                return self._addMotif(Occ.getData('Variant'),Occ.getData('Variant'))
            ### Finish ###
            return None
        except:
            self.log.errorLog('Problem with MotifList.mapMotif()',quitchoice=True)
            return None
#########################################################################################################################
    def mapPattern(self,pattern,update=True):     ### Returns Motif Object if pattern matches original sequence 
        '''
        Returns Motif Object based on self.slims(). Will add if update=True and pattern missing.
        >> pattern:str = motif pattern to check
        >> update:bool [True] = whether to update own list['Motifs'] and/or Motif objects
        '''
        try:
            ### Check for Motif Name with compatible varlist ###
            for Motif in self.slims()[0:]:
                if not Motif: self.list['Motifs'].remove(Motif)
                elif Motif.info['Sequence'] == pattern: return Motif
            ### Make new Motif ###
            if update:
                return self._addMotif(pattern,pattern)
            ### Finish ###
            return None
        except:
            self.log.errorLog('Problem with MotifList.mapMotif()',quitchoice=True)
            return None
#########################################################################################################################
    ### <3> ### Methods for Loading and Reformatting motifs                                                             #
#########################################################################################################################
    def addOcc(self,Seq=None,Motif=None,data={},merge=False):     ### Adds a MotifOcc object with the data given
        '''
        Adds a MotifOcc object with the data given.
        >> Seq:Sequence object against in which the occurrence lies
        >> Motif:Motif object
        >> data:dictionary of data to add {'Info':{},'Stat':{},'Data':{}}
        >> merge:bool [False] = whether to check for existing occurrence and merge if found
        '''
        try:
            ### Make new MotifOcc object ###
            newocc = rje_motifocc.MotifOcc(self.log,self.cmd_list)
            ### Add Data ###
            newocc.obj = {'Seq':Seq,'Motif':Motif}
            if data.has_key('Info'): newocc.setInfo(data['Info'])
            if data.has_key('Stat'): newocc.setStat(data['Stat'])
            if data.has_key('Data'): newocc.dict['Data'] = data['Data']
            ### Add and Return ###
            if merge:
                oldocc = self.checkForOcc(newocc,merge=True)
                if oldocc != newocc: return oldocc   ## Added to existing Occ ##
            self.list['MotifOcc'].append(newocc)
            return newocc
        except:
            self.log.errorLog('Problem adding MotifOcc')
            return None
#########################################################################################################################
    def mergeOcc(self,occlist,overwrite=False):     ### Adds Info, Stat and Data from occlist[1:] to occlist[0]
        '''
        Adds Info, Stat and Data from occlist[1:] to occlist[0]. If overwrite=False, then the occurrences should be in
        order of preferred data, as only missing values will be added. If overwrite=True, then later MotifOcc in the list
        will overwrite the values of the earlier ones. (Note that it is always the first MotifOcc object that is changed.
        Note also that no checks are made that the objects *should* be merged. (See self.checkForOcc())
        >> occlist:list of MotifOcc to merge. Will add to occlist[0] from occlist[1:]
        >> overwrite:bool [False] = If False, will only add missing Info/Stat/Data entries. If True, will add all.
        '''
        try:
            ### Setup focal Occ ###
            MainOcc = occlist[0]
            ### Merge ##
            for Occ in occlist[1:]:
                if overwrite:
                    MainOcc.setInfo(Occ.info)
                    MainOcc.setStat(Occ.stat)
                    for key in Occ.dict['Data']:
                        MainOcc.dict['Data'][key] = Occ.dict['Data'][key]
                else:
                    for key in Occ.info:
                        if not MainOcc.info.has_key(key):
                            MainOcc.info[key] = Occ.info[key]
                    for key in Occ.stat:
                        if not MainOcc.stat.has_key(key):
                            MainOcc.stat[key] = Occ.stat[key]
                    for key in Occ.dict['Data']:
                        if not MainOcc.dict['Data'].has_key(key):
                            MainOcc.dict['Data'][key] = Occ.dict['Data'][key]
        except:
            self.log.errorLog('Problem with MotifList.mergeOcc()')
            raise
#########################################################################################################################
    ### <4> ### Methods for Motif stats and filters                                                                     #
#########################################################################################################################
    def rankMotifs(self,stat,cutoff=0): ### Reranks Motifs using stat and reduces to cutoff if given
        '''
        Reranks Motifs using stat and reduces to cutoff if given.
        >> stat:str = Stat to use for ranking
        >> cutoff:int [0] = number of top ranks to keep
        '''
        rev = True
        if stat in ['Rank']: rev = False     # Low is good!
        premotifs = self.slims()[0:]
        newmotifs = rje_scoring.rankObj(self,self.slims()[0:],stat,cutoff=cutoff,rev=rev)   ### Ranks objects using numerical data
        for Motif in premotifs:
            if Motif not in newmotifs: self.removeMotif(Motif)
        self.list['Motifs'] = newmotifs[0:]
#########################################################################################################################
    def patternStats(self,log=False):      ### Performs calculations all motifs based on basic pattern (info['Sequence'])
        '''Performs calculations all motifs based on basic pattern (info['Sequence']), adding to Motif.stat/info/opt.'''
        mx = 0.0
        for Motif in self.slims():
            if log: self.log.printLog('\r#STAT','SLiM Stats: %.f%%' % (mx/self.motifNum()),log=False,newline=False)
            mx += 100.0
            Motif.patternStats()
        if log: self.log.printLog('\r#STAT','SLiM Stats complete!')
#########################################################################################################################
    def setupDomFilter(self):   ### Sets up self.dict['DomFilter']
        '''Sets up self.dict['DomFilter'].'''
        try:
            ### Read in data ###
            if self.info['DomFilter'].lower() in ['','none']:
                return False
            if not os.path.exists(self.info['DomFilter']):
                self.log.errorLog('Cannot find DomFilter file "%s" - will not use' % self.info['DomFilter'],printerror=False)
                return False
            domdata = rje.dataDict(self,self.info['DomFilter'],mainkeys=['Name','Start','End'],datakeys=['Type'],delimit='\t')
            ### Process Stage 1 ###
            self.dict['DomFilter'] = {}
            for keyset in domdata.keys():
                [name,start,end] = string.split(keyset,'\t')
                if not self.dict['DomFilter'].has_key(name):
                    self.dict['DomFilter'][name] = []
                self.dict['DomFilter'][name].append((string.atoi(start),string.atoi(end)))
            ### Process Stage 2 ###
            for name in self.dict['DomFilter'].keys():
                domlist = self.dict['DomFilter'][name][0:]
                lenlist = []
                for dom in domlist:
                    lenlist.append(dom[1]-dom[0])
                ranklist = rje.rankList(lenlist,absolute=True,lowest=True)
                newlist = []
                for r in range(1,len(ranklist)+1):
                    for i in range(len(ranklist)):
                        if ranklist[i] == r:
                            newlist.append(domlist[i])
                if len(newlist) != len(domlist):
                    raise ValueError
                self.dict['DomFilter'][name] = newlist[0:]
                
            return True
        except:
            self.log.errorLog('Error in setupDomFilter()')
            return False
#########################################################################################################################
    def statFilterMotifs(self,statfilter):   ### Filters motifs using statfilter
        '''Filters motifs using statfilter.'''
        premotifs = self.slims()[0:]
        newmotifs = rje_scoring.statFilterObj(self,self.slims()[0:],statfilter)
        for Motif in premotifs:
            if Motif not in newmotifs:
                self.removeMotif(Motif)
#########################################################################################################################

    def statFilter(self,occdata={},statfilter={}):  ### Filters using rje_scoring and removes filtered MotifOcc
        '''
        Filters using rje_scoring and removes filtered MotifOcc.
        >> occdata:dictionary of {Occ:Statdict}
        >> statfilter:dictionary to statfilters
        << occdata: returns reduced occdata dictionary (same object - pre-reassign if original required!)
        '''
        try:    #!# Add Restrict/Exclude? #!#
            occlist = occdata.keys()[0:]
            occdata = rje_scoring.statFilter(self,occdata,statfilter)
            for Occ in occlist[0:]:
                if Occ in self.list['MotifOcc'] and not occdata.has_key(Occ):
                    self.list['MotifOcc'].remove(Occ)
            return occdata
        except:
            self.log.errorLog('Problem with MotifList.statFilter()')
            return occdata
        
#########################################################################################################################


    def OLDstatFilter(self,datadict):  ### Returns True if motif should be kept according to statfilter
        #!# Replace with rje_scoring #!#
        '''
        Returns True if motif should be kept according to self.obj['Presto'].list['StatFilter']. 
        >> datadict:dictionary of hit stats
        << True/False if accepted/filtered
        '''
        try:
            ### Setup ###
            presto = self.obj['Presto']

            ### Restricted and Exclusive Masking of Motifs ###
            motif = datadict['MOTIF']
            vmotif = string.replace(string.replace(motif,'_fix',''),'_var','')
            vmatch = '=%s' % datadict['MATCHSEQ']
            if presto.dict['Restrict']:
                if presto.dict['Restrict'].has_key(vmotif):
                    #X#presto.deBug('\n(%s,%s,%s)\n' % (datadict['HIT'],datadict['START_POS'],datadict['END_POS']))
                    #X#presto.deBug('%s: %s' % (motif,presto.dict['Restrict'][motif]))
                    if [datadict['HIT'],datadict['START_POS'],datadict['END_POS']] not in presto.dict['Restrict'][vmotif]:
                        return False
                if presto.dict['Restrict'].has_key(vmatch):
                    #X#presto.deBug('\n(%s,%s,%s)\n' % (datadict['HIT'],datadict['START_POS'],datadict['END_POS']))
                    #X#presto.deBug('%s: %s' % (motif,presto.dict['Restrict'][motif]))
                    if [datadict['HIT'],datadict['START_POS'],datadict['END_POS']] not in presto.dict['Restrict'][vmatch]:
                        return False
                if not presto.dict['Restrict'].has_key(vmatch) and not presto.dict['Restrict'].has_key(vmotif):
                    return False
            if presto.dict['Exclude']:
                if presto.dict['Exclude'].has_key(vmotif):
                    if [datadict['HIT'],datadict['START_POS'],datadict['END_POS']] in presto.dict['Exclude'][vmotif]:
                        return False            
                if presto.dict['Exclude'].has_key(vmatch):
                    if [datadict['HIT'],datadict['START_POS'],datadict['END_POS']] in presto.dict['Exclude'][vmatch]:
                        return False            

            ### Filter patterns ###
            #!# Add presto.list['StatFilter'] when adding NewScore=X commands #!#
            filtdata = statFilter(self,data={motif:datadict},statfilter=presto.dict['StatFilter'])
            if filtdata.has_key(motif):
                return True

            ### Finish ###
            return False
        except:
            self.log.errorLog('Error in SlimPicker.reRank()',printerror=True,quitchoice=True)
            return False
#########################################################################################################################            
    ### <5> ### Motif output methods                                                                                    #
#########################################################################################################################
## Some of these outputs were originally designed for slim_pickings but have been generalised for MotifOcc results:
## - motifOut outputs motifs in PRESTO format

## - MotifAln produces a single file for the dataset of alignments of each motif occurrence [Memsaver=F only]
## - ProteinAln produces a single file per sequence of all the motifs positioned on the sequence alignmnet
## - FTOut produces a single file for the dataset of extracted UniProt features with the motif positions added
## - MotifInfo produces a summary table of motif information
#########################################################################################################################
    def outputs(self):  ### Processes addition outputs for motif occurrences (Alignments and Features)
        '''Processes addition outputs for motif occurrences (Alignments and Features).'''
        try: ###???###
            if self.opt['FTOut']:
                self.FTOut()
            
        except:
            self.log.errorLog('Problem with rje_motiflist.outputs()',quitchoice=True)
#########################################################################################################################
    def FTOut(self,acc_occ):    ### Produces a single for of extracted UniProt features with the motif positions added
        '''
        Produces a single file for the dataset of extracted UniProt features with the motif positions added.
        >> acc_occ:dictionary of {AccNum:{Motif:[(position:match)]}
        '''
        try:
            ### Setup Features Dictionary ###
            slim_ft = {}
            for acc in acc_occ.keys()[0:]:
                slim_ft[acc] = []
                svacc = None
                if rje.matchExp('^(\S+)\-(\d+)',acc):   # Splice variant
                    svacc = rje.matchExp('^(\S+)\-(\d+)',acc)[0]
                    slim_ft[svacc] = []
                mx = 0  # Count of motif features
                for motif in acc_occ[acc].keys():
                    for (pos,match) in acc_occ[acc][motif]:
                        mx += 1
                        my_ft = {'Type':'PRESTO','Start':pos,'End':(pos+len(match)-1),
                                 'Desc':'PRESTO match (%s) to %s (%s).' % (match,motif.info['Name'],motif.info['Sequence'])}
                        slim_ft[acc].append(my_ft)
                        if svacc:   # Splice variant
                            sv_ft = {'Type':my_ft['Type'],'Start':my_ft['Start'],'End':my_ft['End'],
                                     'Desc':'%s (Splicevar %s)' % (my_ft['Desc'],acc)}
                            slim_ft[svacc].append(sv_ft)
                if not mx:
                    acc_occ.pop(acc)

            ### Extract uni_extract from UniProt ###
            ## Setup UniProt File ##
            my_entries = []
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
            uniprot.list['Extract'] = rje.sortKeys(acc_occ)
            uniprot.opt['MemSaver'] = False
            uniprot.run()
            my_entries += uniprot.list['Entry'][0:]
            unipaths = []
            if self.info['UniPaths'] not in ['','None']:
                unipaths = string.split(self.info['UniPaths'],',')
            for path in unipaths:
                uniprot.opt['Append'] = True
                uniprot.info['UniPath'] = rje.makePath(path,return_blank=False)
                uniprot.run()   #!# Add DomTable = Makes a table of domains from uniprot file [False] #!#
                my_entries += uniprot.list['Entry'][0:]     #!# Add features to existing entries at some point #!#

            ### Features Table incorporating SLiMs ###
            # Extract list is AccNums from slim_seq seqlist.
            uniprot.opt['Append'] = self.opt['Append']
            accout = []
            for entry in my_entries:
                acc = entry.obj['Sequence'].info['AccNum']
                if slim_ft.has_key(acc):
                    entry.list['Feature'] += slim_ft[acc]
                accout.append(acc)
            for acc in slim_ft.keys():
                svacc = None
                if rje.matchExp('^(\S+)\-(\d+)',acc):
                    svacc = rje.matchExp('^(\S+)\-(\d+)',acc)[0]
                if acc not in accout and (not svacc or svacc not in accout):
                    _entry = rje_uniprot.UniProtEntry(log=self.log,cmd_list=self.cmd_list)
                    _entry.obj['Sequence'].info['AccNum'] = acc
                    _entry.list['Feature'] = slim_ft[acc]
                    my_entries.append(_entry)
            uniprot.list['Entry'] = my_entries[0:]
            uniprot.ftTable('%s.features.tdt' % self.info['ResFile'])
        except:
            self.log.errorLog('Major problem during PRESTO.FTOut()')
#########################################################################################################################
    def motifInfo(self):   ### Produces summary table for motifs, including expected values and information content.
        '''Produces summary table for motifs, including expected values and information content.'''
        try:
            #!# Needs tidying. Does not currently output expected values (from what?) #!#
            #!# Add option to return a database table instead of output to file #!#
            ### Setup Results ###
            outfile = self.info['MotInfo']
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            rje.backup(self,outfile)    # Backs up if existing and not self.opt['Append']
            headers = ['Motif','Pattern','Description','PosLength','MinLength','MaxLength','FixLength']
            if not self.opt.has_key('MotifIC') or self.opt['MotifIC']: headers.append('IC')
            rje.delimitedFileOutput(self,outfile,headers,delimit,datadict={})   # Output headers

            ### Output ###
            for motif in self.motifs():
                datadict = {'Motif':motif.info['Name'],
                            'Pattern':motif.pattern(),
                            'Description':motif.info['Description'],
                            'IC':'%.3f' % motif.stat['IC'],
                            'MaxLength':motif.slimLen(),
                            'MinLength':motif.slimMinLen(),'FixLength':motif.slimFix(),'PosLength':motif.slimPos()
                            }
                rje.delimitedFileOutput(self,outfile,headers,delimit,datadict)
            self.printLog('#INFO','Motif info output to %s.' % outfile)
            
        except: self.log.errorLog('Error in SLiMList.motifInfo()',quitchoice=True)
#########################################################################################################################
    ### <5> ### Expectation methods                                                                                     #
#########################################################################################################################
    def seqExp(self,seq=None):   ### Populates motif.dict[Expect(MM)] for single Sequence object
        '''
        Populates self.dict[Expect].
        >> seq:Sequence object to consider
        '''
        try:
            if not seq:
                return
            aafreq = rje.dictFreq(seq.aaFreq(),total=True)
            for Motif in self.list['Motifs']:
                Motif.dict['ExpectMM'][seq] = Motif.expectation(aafreq,aanum=aafreq['Total'],seqnum=1)
                Motif.dict['Expect'][seq] = Motif.dict['ExpectMM'][seq][0]
        except:
            self.log.errorLog('Problem with Motif._seqExp(%s)' % seq.shortName(),quitchoice=True)
#########################################################################################################################
    def seqListExp(self,seqlist=None,filename='',cutoff=0):   ### Populates self.dict[Expect(MM)] for SeqList or seqfile
        '''
        Populates self.dict[Expect].
        >> seqlist:SeqList object to consider
        >> filename:Sequence files to consider (must be fasta)
        >> cutoff:float [0] = Expectation cutoff, will remove motif if exceeded.
        '''
        try:
            ### Setup ###
            if seqlist:
                aafreq = seqlist.aaFreq(alphabet=rje_seq.alph_protx,fromfile=None,total=True)
                seqnum = seqlist.seqNum()
            elif filename.lower() not in ['','none'] and os.path.exists(filename):
                aafreq = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F','seqin=None']).aaFreq(alphabet=rje_seq.alph_protx,fromfile=filename,total=True)
                seqnum = rje_seq.SeqCount(self,filename)
            else:
                self.log.errorLog('No seqlist given for seqListExp and file "%s" not found' % filename,printerror=False)
                return
            if cutoff < 0:
                cutoff = float(seqnum) / -cutoff
            ### Expect and Cutoff ###            
            for Motif in self.list['Motifs'][0:]:
                Motif.dict['ExpectMM'][filename] = Motif.expectation(aafreq,aanum=aafreq['Total'],seqnum=seqnum)
                Motif.dict['Expect'][filename] = Motif.dict['ExpectMM'][filename][0]
                if cutoff > 0 and Motif.dict['Expect'][filename] > cutoff:
                    self.removeMotif(Motif,remtxt='Motif "%s" discarded: ExpOcc = %s > ExpCut (%s)' % (Motif.info['Name'],rje_motif.expectString(Motif.dict['Expect'][filename]),rje_motif.expectString(cutoff)))
        except:
            self.log.errorLog('Problem with Motif._seqListExp()',quitchoice=True)
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
        
    ### Rest of Functionality... ###
    try: SLiMList(mainlog,cmd_list).run()

    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
