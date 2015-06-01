#!/usr/local/bin/python

# MotifList - RJE MotifList Object (based on PRESTO)
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_motiflist
Description:  RJE Motif List Module
Version:      1.0
Last Edit:    03/04/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the MotifList Class, which is designed to replace many of the functions that previously formed
    part of the Presto Class. This class will then be used by PRESTO, SLiMPickings and CompariMotif (and others?) to
    control Motif loading, redundancy and storage. MotifOcc objects will replace the previous PrestoSeqHit objects and
    contain improved data commenting and retrieval methods. The MotifList class will contain methods for filtering motifs
    according to individual or combined MotifOcc data.

    The options below should be read in by the MotifList object when it is instanced with a cmd_list and therefore do not
    need to be part of any class that makes use of this object unless it has conflicting settings.

    The Motif Stats options are used by MotifList to calculate statistics for motif occurrences, though this data will
    actually be stored in the MotifOcc objects themselves. This includes conservation statistics.

    Note. Additional output parameters, such as motifaln and proteinaln settings, and stat filtering/novel scores are not
    stored in this object, as they will be largely dependent on the main programs using the class, and the output from
    those programs. (This also enables statfilters etc. to be used with stats not related to motifs and their occurrences
    if desired.)

MotifList Commands:
    ## Basic Motif Input/Formatting Parameters ##
    motifs=FILE     : File of input motifs/peptides [None]
                      Single line per motif format = 'Name Sequence #Comments' (Comments are optional and ignored)
                      Alternative formats include fasta, SLiMDisc output and raw motif lists.
    minpep=X        : Min length of motif/peptide X aa [2]
    minfix=X        : Min number of fixed positions for a motif to contain [0]
    minic=X         : Min information content for a motif (1 fixed position = 1.0) [2.0]
    trimx=T/F       : Trims Xs from the ends of a motif [False]
    nrmotif=T/F     : Whether to remove redundancy in input motifs [False]
    minimotif=T/F   : Input file is in minimotif format and will be reformatted (PRESTO File format only) [False]
    goodmotif=LIST  : List of text to match in Motif names to keep (can have wildcards) []
    ambcut=X        : Cut-off for max number of choices in ambiguous position to be shown as variant [10]
    reverse=T/F     : Reverse the motifs - good for generating a test comparison data set [False]
    msms=T/F        : Whether to include MSMS ambiguities when formatting motifs [False]

    ## Motif Occurrence Statistics Options ##
    winsa=X         : Number of aa to extend Surface Accessibility calculation either side of motif [0]
    winhyd=X        : Number of aa to extend Eisenberg Hydrophobicity calculation either side of motif [0]
    windis=X        : Extend disorder statistic X aa either side of motif (use flanks *only* if negative) [0]
    winchg=X        : Extend charge calculations (if any) to X aa either side of motif [0]
    winsize=X       : Sets all of the above window sizes (use flanks *only* if negative) [0]
    slimchg=T/F     : Calculate Asolute, Net and Balance charge statistics (above) for occurrences [False]
    iupred=T/F      : Run IUPred disorder prediction [False]
    foldindex=T/F   : Run FoldIndex disorder prediction [False]
    iucut=X         : Cut-off for IUPred results (0.0 will report mean IUPred score) [0.0]
    iumethod=X      : IUPred method to use (long/short) [short]
    domfilter=FILE  : Use the DomFilter options, reading domains from FILE [None] ?? Check how this works ??
    ftout=T/F       : Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    percentile=X    : Percentile steps to return in addition to mean [0]

    ## Conservation Parameters ##   ??? Add separate SlimCons option ???
    usealn=T/F      : Whether to search for and use alignemnts where present. [False]
    gopher=T/F      : Use GOPHER to generate missing orthologue alignments in alndir - see gopher.py options [False]
    alndir=PATH     : Path to alignments of proteins containing motifs [./] * Use forward slashes (/)
    alnext=X        : File extension of alignment files, accnum.X [aln.fas]
    alngap=T/F      : Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore
                      as putative sequence fragments [False]  (NB. All X regions are ignored as sequence errors.)
    conspec=LIST    : List of species codes for conservation analysis. Can be name of file containing list. [None]
    conscore=X      : Type of conservation score used:  [pos]
                        - abs = absolute conservation of motif using RegExp over matched region
                        - pos = positional conservation: each position treated independently 
                        - prop = conservation of amino acid properties
                        - all = all three methods for comparison purposes
    consamb=T/F     : Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    consinfo=T/F    : Weight positions by information content (does nothing for conscore=abs) [True]
    consweight=X    : Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
                        - 0 gives equal weighting to all. Negative values will upweight distant sequences.
    posmatrix=FILE  : Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) [None]
    aaprop=FILE     : Amino Acid property matrix file. [aaprop.txt]

    ## Alignment Settings ##
    protalndir=PATH : Output path for Protein Alignments [ProteinAln/]
    motalndir=PATH  : Output path for Motif Alignments []
    flanksize=X     : Size of sequence flanks for motifs [30]
    xdivide=X       : Size of dividing Xs between motifs [10]

    ## System Settings ##
    iupath=PATH     : The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    ?? memsaver=T/F    : Whether to store all results in Objects (False) or clear as search proceeds (True) [True] ??
    ?- should this be controlled purely by the calling program? Probably!
    fullforce=T/F   : Whether to force regeneration of alignments using GOPHER

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_aaprop, rje_disorder, rje_motif_V3, rje_motif_cons, rje_scoring, rje_seq, rje_sequence,
    rje_blast, rje_uniprot
Other modules needed: rje_dismatrix, 
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_aaprop, rje_motif_stats, rje_scoring, rje_seq
import rje_motif_V3 as rje_motif
import rje_uniprot, rje_motifocc
#########################################################################################################################
### History
# 0.0 - Initial construction of module based on PRESTO for use with PRESTO V5.0 and SLiMPickings V3.0.
#########################################################################################################################
### Major Functionality to Add
# [ ] : Finish initial compilation based on PRESTO.
# [ ] : Check that all the imported modules are needed.
# [ ] : Check use of Alphabet and assess whether it can be dumped.
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit) = ('RJE_MotifList', '1.0', 'January 2007')  
    description = 'RJE Motif List Module'
    author = 'Dr Richard J. Edwards.'
    return rje.Info(program,version,last_edit,description,author,time.time())
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program,info.version,time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,-1,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
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
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: MotifList Class:                                                                                        #
#########################################################################################################################
class MotifList(rje.RJE_Object):     
    '''
    Motif List Class for use with PRESTO, SLiMPickings etc. Author: Rich Edwards (2007). Based on PRESTO 4.2.

    Info:str
    - Name = Name of input motif/peptide file
    - MotifOut = Filename for output of reformatted (and filtered?) motifs in PRESTO format [None]
    - DomFilter = Use the DomFilter options, reading domains from FILE [None]
    - AlnDir = Path to alignment files
    - AlnExt = File extensions of alignments: AccNum.X
    - ConScore = Type of conservation score used:  [abs]
        - abs = absolute conservation of motif: reports percentage of homologues in which conserved
        - prop = conservation of amino acid properties
        - pos = positional conservation: each position treated independently 
        - all = all three methods for comparison purposes
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix)
    - ProtAlnDir = Directory name for output of protein aligments [ProteinAln/]
    - MotAlnDir = Directory name for output of protein aligments []
    
    Opt:boolean
    - NRMotif = Whether to remove redundancy in input motifs [False]
    - Expect = Whether to calculate crude 'expected' values based on AA composition.
    - MSMS = Whether to run in MSMS mode
    - UseAln = Whether to look for conservation in alignments
    - Reverse = Reverse the motifs - good for generating a test comparison data set [False]
    - IUPred = Run IUPred disorder prediction [False]
    - FoldIndex = Run FoldIndex disorder prediction [False]
    - ConsInfo = Weight positions by information content [True]
    - AlnGap = Whether to count proteins in alignments that have 100% gaps over motif (True) or (False) ignore as putative sequence fragments [True]
    - Gopher = Use GOPHER to generate missing orthologue alignments in outdir/Gopher - see gopher.py options [False]
    - ConsAmb = Whether to calculate conservation allowing for degeneracy of motif (True) or of fixed variant (False) [True]
    - TrimX = Trims Xs from the ends of a motif
    - SlimChg = Calculate Asolute, Net and Balance charge statistics (above) for occurrences [False]
    - MiniMotif = Input file is in minimotif format and will be reformatted [False]
    - Compare = whether being called by CompariMotif (has some special requirements!)
    - FullForce = Whether to force regeneration of alignments using GOPHER
    
    Stat:numeric
    - AmbCut = Cut-off for max number of choices in ambiguous position to be shown as variant [10]
        For mismatches, this is the max number of choices for an ambiguity to be replaced with a mismatch wildcard
    - MinPep = Minimum length of motif/peptide (non-X characters)
    - MinFix = Min number of fixed positions for a motif to contain [0]
    - MinIC = Min information content for a motif (1 fixed position = 1.0) [2.0]
    - WinSA = Number of aa to extend Surface Accessibility calculation either side of motif [0]
    - WinHyd = Number of aa to extend Eisenberg Hydrophobicity calculation either side of motif [0]
    - WinDis = Extend disorder statistic X aa either side of motif [0]
    - WinChg = Extend charge calculations (if any) to X aa either side of motif [0]
    - WinSize = Used for peptide design and also to set all of above [0]
    - ConsWeight = Weight given to global percentage identity for conservation, given more weight to closer sequences [0]
    - FlankSize = Size of sequence flanks for motifs in MotifAln [30]
    - XDivide = Size of dividing Xs between motifs [10]
    - FTOut = Make a file of UniProt features for extracted parent proteins, where possible, incoroprating SLIMs [*.features.tdt]
    - Percentile = Percentile steps to return in addition to mean [25]

    List:list
    - Alphabet = List of letters in alphabet of interest
    - Motifs = List of rje_Motif.Motif objects
    - MotifOcc = List of MotifOcc objects 
    - GoodMotif = List of text to match in Motif names to keep (can have wildcards) []

    Dict:dictionary
    - ConsSpecLists = Dictionary of {BaseName:List} lists of species codes for special conservation analyses
    - ElementIC = Dictionary of {Position Element:IC}
    - PosMatrix = Score matrix for amino acid combinations used in pos weighting. (conscore=pos builds from propmatrix) {}

    Obj:RJE_Objects
    - AAPropMatrix = rje_aaprop.AAPropMatrix object
    - SeqList = SeqList object of search database
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics 
        self.infolist = ['Name','AlnDir','AlnExt','ConScore','PosMatrix','ProtAlnDir','DomFilter','MotAlnDir']
        self.statlist = ['AmbCut','MinPep','MinFix','WinSA','WinHyd','WinDis','WinChg','ConsWeight','WinSize','MinIC',
                         'FlankSize','XDivide','Percentile']
        self.optlist = ['Expect','MSMS','UseAln','Reverse','AlnGap','NRMotif','IUPred','FoldIndex','Compare',
                        'ConsInfo','ConsAmb','Gopher','TrimX','SlimChg','MiniMotif','FTOut','FullForce']
        self.listlist = ['Motifs','Alphabet','MotifOcc','GoodMotif']
        self.dictlist = ['ConsSpecLists','ElementIC','PosMatrix','MisMatch']
        self.objlist = ['SeqList','AAPropMatrix']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'MinFix':0,'MinPep':2,'FlankSize':30,'XDivide':10,'MinShare':2,'AmbCut':10,'MotDesc':3,
                      'XPad':0,'XPadDB':0,'MinIC':2.0,'Percentile':0})
        for w in ['SA','Hyd','Dis','Chg','Size']:
            self.stat['Win%s' % w] = 0
        self.setOpt({'Expect':True,'MatchIC':True,'MemSaver':True,'ConsAmb':True,'ConsInfo':True})
        self.setInfo({'AlnDir':'Gopher/ALN/','AlnExt':'orthaln.fas','ConScore':'pos',
                      'ProtAlnDir':rje.makePath('ProteinAln/'),'MotifAlnDir':''})
        self.list['Alphabet'] = string.split('A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y',',')
#########################################################################################################################
    def clear(self):    ### Clears objects and lists
        '''Clears objects and lists.'''
        self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
        self.list['Motifs'] = []
        self.list['MotifOcc'] = []
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

                ### Basic Options ###
                self._cmdRead(cmd,type='file',att='Name',arg='motifs')  
                self._cmdReadList(cmd,'info',['AlnExt','ConScore'])
                self._cmdReadList(cmd,'path',['AlnDir','ProtAlnDir','MotAlnDir'])
                self._cmdReadList(cmd,'file',['PosMatrix','DomFilter'])
                self._cmdReadList(cmd,'int',['AmbCut','MinPep','MinFix','WinSA','WinHyd','WinDis','WinChg','ConsWeight',
                                             'WinSize','FlankSize','XDivide'])
                self._cmdReadList(cmd,'stat',['MinIC'])
                self._cmdReadList(cmd,'list',['GoodMotif','Alphabet'])
                self._cmdReadList(cmd,'opt',['Expect','MSMS','UseAln','Reverse','AlnGap','NRMotif','IUPred','FoldIndex',
                                             'ConsInfo','ConsAmb','Gopher','TrimX','SlimChg','MiniMotif','FTOut',
                                             'FullForce'])
                self._cmdRead(cmd,type='opt',att='UseAln',arg='slimcons')

                ### Extra Options for SLiMPickings Compatibility ###
                self._cmdRead(cmd,'opt','IUPred','slimiup')
                self._cmdRead(cmd,'opt','FoldIndex','slimfold')

                ### Special ###
                if rje.matchExp('^conspec=(.+)',cmd):
                    conspec = rje.matchExp('^conspec=(.+)',cmd)[0]
                    consfiles = glob.glob(conspec)
                    if len(consfiles) > 0:     # File
                        for cfile in consfiles:
                            self.dict['ConsSpecLists'][rje.baseFile(cfile,True)] = rje.listFromCommand(cfile)
                    else:
                        self.dict['ConsSpecLists']['SPEC'] = rje.listFromCommand(conspec,checkfile=False)
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)

        ### Setup GopherDir and AlnDir ##
        self.setupGopher()
                
        ### Conservation Attributes ###
        if self.opt['UseAln'] and self.info['ConScore'].lower() in ['all','prop']:
            self.obj['AAPropMatrix'] = rje_aaprop.AAPropMatrix(self.log,self.cmd_list)
        if self.opt['UseAln'] and self.info['ConScore'].lower() in ['all','prop','pos']: self.posMatrix()    
#########################################################################################################################
    ### <2> ### Methods for returning Occurrences and MotifList attributes                                              #
#########################################################################################################################
    def motifNum(self): return len(self.list['Motifs'])
#########################################################################################################################
    def motifs(self): return self.list['Motifs']
#########################################################################################################################
    def motifOcc(self,byseq=False,justdata=None,fastacmd=False,nested=True):    ### Returns MotifOcc list/dictionary 
        '''
        Returns a list or dictionary of MotifOccurrences.
        >> byseq:bool [False] = return a dictionary of {Sequence:{Motif:OccList}}, else {Motif:{Sequence:OccList}} 
        >> justdata:str [None] = if given a value, will return this data entry for each occ rather than the object itself
        >> fastacmd:bool [False] = whether to return FastaCmd instead of Sequence if Sequence missing
        >> nested:bool [True] = whether to return a nested dictionary or just the occurrences per Motif/Seq
        << returns dictionary or plain list of MotifOcc if byseq and bymotif are both False        
        '''
        try:
            motifocc = {}
            for occ in self.list['MotifOcc']:
                ## Setup Occ ##
                (Seq,Motif) = (occ.obj['Seq'],occ.obj['Motif'])
                if fastacmd and not Seq: Seq = occ.getData('FastaCmd')
                if byseq: (key1,key2) = (Seq,Motif)
                else: (key2,key1) = (Seq,Motif)
                if not motifocc.has_key(key1):
                    if nested: motifocc[key1] = {}
                    else: motifocc[key1] = []
                if nested and not motifocc[key1].has_key(key2): motifocc[key1][key2] = []
                ## Update ##
                if justdata: val = occ.getData(justdata)
                else: val = occ
                if nested: motifocc[key1][key2].append(val)
                else: motifocc[key1].append(val)
            return motifocc
        except:
            self.log.errorLog('Problem with MotifList.motifOcc()')
            raise
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
    def mapMotif(self,Occ,update=True):     ### Returns Motif Object based on self.motifs() and Occ data
        '''
        Returns Motif Object and updates Occ, based on self.motifs() and Occ data.
        >> Occ:MotifOcc object to check
        >> update:bool [True] = whether to update own list['Motifs'] and/or Motif objects
        '''
        try:
            ### Check for exact Motif ###
            if Occ.obj['Motif'] in self.motifs():
                return Occ.obj['Motif']
            ### Check for Motif Name with compatible varlist ###
            for Motif in self.motifs():
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
            for Motif in self.motifs():
                if Motif.info['Name'] == Occ.getData('Variant'):
                    return Motif
            ### Check varlist only ###
            if not Occ.getData('Motif') and Occ.getData('Variant'):
                for Motif in self.motifs():
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
        Returns Motif Object based on self.motifs(). Will add if update=True and pattern missing.
        >> pattern:str = motif pattern to check
        >> update:bool [True] = whether to update own list['Motifs'] and/or Motif objects
        '''
        try:
            ### Check for Motif Name with compatible varlist ###
            for Motif in self.motifs()[0:]:
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
    def loadMotifs(self,file=None,clear=True):      ### Loads motifs and populates self.list['Motifs']
        '''
        Loads motifs and populates self.list['Motifs'].
        >> file:str = filename or self.info['Name'] if None
        >> clear:boolean = whether to clear self.list['Motifs'] before loading [True]
        '''
        try:
            ### Setup ###
            if file: self.info['Name'] = file
            else: file = self.info['Name']
            if clear: self.list['Motifs'] = []
            preloadx = self.motifNum()

            ### Read File ###
            mlines = []
            if os.path.exists(file): mlines = self.loadFromFile(file)
            if not mlines:
                if file.find(',') > 0:   # Motif List
                    mlines = string.split(file,',')
                elif file.lower() not in ['','none'] and (self.stat['Interactive'] < 1 or rje.yesNo('No lines read from "%s". Use as motif?' % file)):
                    self.log.printLog('#MOT','No lines read from "%s": using as motif.' % file)
                    mlines = [file]
                else: return  False

            ### Process MiniMotif ###
            if self.opt['MiniMotif']:
                self._reformatMiniMotif(mlines)
                self.log.printLog('#MOT','Motifs read from %s (%d lines): %d retained.' % (file,len(mlines),self.motifNum()-preloadx))
                return False

            ## Process other formats ##            
            mx = 0  # Number read
            if mlines[0][0] == '>':     # Motifs from Fasta file
                name = None
                seq = None
                while mlines:
                    line = rje.chomp(mlines.pop(0))
                    if line[0] == '>':
                        if name:
                            mx += 1
                            self._addMotif(name=name,seq=seq,reverse=self.opt['Reverse'],check=True)
                            if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
                        name = line[1:]
                        seq = ''
                    else:
                        seq += line
            elif mlines[0][:2] == '##' and rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',mlines[5]):  # TEIRESIAS
                for line in mlines:
                    motif = rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',line)
                    if motif:
                        mx += 1
                        self._addMotif(name=motif[2],seq=motif[2],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            elif mlines[0][:2] == '#-' and rje.matchExp('Input dataset:\s+(\S+)',mlines[1]):  # SLiMDisc
                name_root = os.path.splitext(os.path.basename(rje.matchExp('Input dataset:\s+(\S+)',mlines[1])[0]))[0]
                for line in mlines:
                    motif = rje.matchExp('^\((\d+)\)\s+\S+\s+(\S+)',line)
                    if motif:
                        mx += 1
                        self._addMotif(name='%s#%s' % (name_root,motif[0]),seq=motif[1],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            elif mlines[0].find('Dataset,RunID,') == 0:     # SLiMFinder
                delimit = rje.delimitFromExt(filename=file)
                sfdata = rje.dataDict(self,file,mainkeys=['Dataset','RunID','Rank','Pattern'],datakeys=['Dataset','RunID','Rank','Pattern'])
                for returned in rje.sortKeys(sfdata):
                    data = sfdata[returned]
                    if data['Pattern']:
                        name = string.replace(string.join([data['Dataset'],data['RunID'],data['Rank'],data['Pattern']],'#'),' ','_')
                        mx += 1
                        self._addMotif(name=name,seq=data['Pattern'],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            elif mlines[0].find('Dataset,SeqNum,TotalAA,Rank,Score,Pattern') == 0:      # SlimPickings (SlimDisc 1.3)
                for line in mlines[1:]:
                    motif = string.split(line,',')
                    if len(motif) > 5:
                        mx += 1
                        self._addMotif(name='%s#%s' % (motif[0],motif[3]),seq=motif[5],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            elif mlines[0].find('Dataset,SeqNum,FullMST,TotalAA,Rank,Score,Pattern') == 0:      # SlimPickings (SlimDisc 1.4)
                for line in mlines[1:]:
                    motif = string.split(line,',')
                    if len(motif) > 6:
                        mx += 1
                        self._addMotif(name='%s#%s' % (motif[0],motif[4]),seq=motif[6],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
            else:                
                for line in mlines:
                    desc = ''
                    if line.find('#') >= 0:
                        desc = rje.chomp(line[line.find('#')+1:])
                        while desc[:1] == ' ':
                            desc = desc[1:]
                        line = line[:line.find('#')]
                    motif = rje.matchExp('^(\S+)\s+(\S+)',line)
                    if motif:
                        mx += 1
                        self._addMotif(name='%s %s' % (motif[0],desc),seq=motif[1],reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
                    elif rje.matchExp('^(\S+)',line):   # Pure patterns: name = pattern
                        mx += 1
                        motif = rje.matchExp('^(\S+)',line)[0]
                        self._addMotif(name='%s %s' % (motif,desc),seq=motif,reverse=self.opt['Reverse'],check=True)
                        if self.stat['Verbose'] < 2: self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx),log=False,newline=False)
                
            ### Summarise ###                            
            self.log.printLog('\r#MOT','%d motifs read from %s (%d lines): %d retained.' % (mx,file,len(mlines),self.motifNum()-preloadx))
            return True

        except:
            self.log.errorLog('Error in loadMotifs()')     
            raise   
#########################################################################################################################
    def _reformatMiniMotif(self,mlines):   ### Reformats MiniMotif file, compressing motifs as appropriate
        '''
        Reformats MiniMotif file, compressing motifs as appropriate.
        >> mlines:list of lines read from input file
        '''
        try:
            ### Make dictionaries of {name:[patterns]} and {name:desc} ###
            patdict = {}    # Dictionary of name:patterns
            desdict = {}    # Dictionary of name:description
            adddict = {}    # Dictionary of name and extra letters to add
            for m in mlines:
                ## Get details ##
                details = rje.matchExp('^(\S+)\s+(\S+)\s+#\s*(.+)',m)
                if details:
                    (name,pattern,desc) = details
                    pattern = rje_motif.reformatMiniMotif(self,pattern)
                ## Update dictionaries ##
                for mname in patdict.keys():
                    if mname[:1] == name[:1] and desdict[mname] == desc:    # Add to this motif
                        patdict[mname].append(pattern)
                        if name[-1] not in adddict[mname]:
                            adddict[mname].append(name[-1])
                        name = ''
                        break
                if name:
                    patdict[name] = [pattern]
                    adddict[name] = [name[-1]]
                    desdict[name] = desc

            ### Compress MiniMotifs ###
            for name in rje.sortKeys(patdict):
                newpat = rje_motif.defineMotif(self,patdict[name],profile=False,minfreq=0.05,minocc=1,ambcut=19)
                add = adddict[name]
                add.sort()
                mname = name[:-1] + string.join(add,'')
                self.log.printLog('#MOT','Motif "%s" => %s' % (mname,newpat))
                if len(newpat) > 1:     ### Make several
                    ext = 'abcdefghijklmnopqrstuvwxyz'
                    for i in range(len(newpat)):
                        self._addMotif(name='%s%s %s' % (mname,ext[i],desdict[name]),seq=newpat[i],reverse=self.opt['Reverse'],check=True)
                else:
                    self._addMotif(name='%s %s' % (mname,desdict[name]),seq=newpat[0],reverse=self.opt['Reverse'],check=True)
                
        except:
            self.log.errorLog('Error in _reformatMiniMotif()')     
            raise   
#########################################################################################################################
    def _addMotif(self,name,seq,reverse=False,check=False,logrem=True):  ### Adds new motif to self.list['Motifs']. 
        '''
        Adds new motif to self.list['Motifs']. Checks redundancy etc.
        >> name:str = Motif Name
        >> seq:str = Motif Sequence read from file
        >> reverse:boolean [False] = whether to reverse sequence
        >> check:boolean [False] = whether to check redundancy and sequence length
        >> logrem:boolean [True] = whether to log removal of motifs
        << returns Motif object or None if failed
        '''
        try:
            ### Reverse ##
            desc = string.join(string.split(name)[1:])
            name = string.split(name)[0]
            if reverse:
                name = '%s_rev' % name

            ### Check Name ###
            if self.list['GoodMotif']:
                dump = True
                for good in self.list['GoodMotif']:
                    if rje.matchExp('^(%s)$' % string.replace(good,'*','\S*'),name):
                        dump = False
                        break
                if dump:
                   self.log.printLog('\r#REM','Motif "%s" not in goodmotif=LIST.' % (name)) 

            ### Setup ###
            newmotif = rje_motif.Motif(log=self.log,cmd_list=self.cmd_list)
            newmotif.info['Name'] = name
            if desc:
                newmotif.info['Description'] = desc
            else:
                newmotif.info['Description'] = name
            newmotif.info['Sequence'] = seq
            newmotif.obj['MotifList'] = self
            newmotif.opt['TrimX'] = self.opt['TrimX']
            if not newmotif.format(msmode=self.opt['MSMS'],reverse=reverse):
                self.log.errorLog('Formatting problem with %s sequence %s: Rejected.' % (name,seq),printerror=False,quitchoice=False)
                return None

            ### Check Length ###
            if check and newmotif.stat['FixLength'] < self.stat['MinFix']:
                self.log.printLog('\r#REM','Motif %s (%s) has < %d fixed positions (%d).' % (newmotif.info['Name'], newmotif.info['Sequence'],self.stat['MinFix'],newmotif.stat['FixLength']),screen=logrem,log=logrem)
                return None
            elif check and newmotif.stat['Length'] < self.stat['MinPep']:
                self.log.printLog('\r#REM','Motif %s (%s) is too short (%d non-X aas < %d).' % (newmotif.info['Name'], newmotif.info['Sequence'],newmotif.stat['Length'],self.stat['MinPep']),screen=logrem,log=logrem)
                return None
            elif check and newmotif.stat['IC'] < self.stat['MinIC']:
                self.log.printLog('\r#REM','Motif %s (%s) has insufficient IC (%.1f < %.1f).' % (newmotif.info['Name'], newmotif.info['Sequence'],newmotif.stat['IC'],self.stat['MinIC']),screen=logrem,log=logrem)
                return None
            else:
                self.verbose(2,4,'%s: %s (%d non-X aa)' % (newmotif.info['Name'], newmotif.info['Sequence'],newmotif.stat['Length']),1)

            ### Check Redundancy ###
            if check and self.opt['NRMotif']:
                for motif in self.list['Motifs'][0:]:
                    if string.join(newmotif.list['PRESTO'],'-') == string.join(motif.list['PRESTO'],'-'):
                        self.log.printLog('\r#REM','Motif %s is identical to %s and has been removed.' % (newmotif.info['Name'],motif.info['Name']),screen=logrem,log=logrem)
                        return None
                    elif string.join(motif.list['PRESTO'],'-').find(string.join(newmotif.list['PRESTO'],'-')) >= 0:
                        self.log.printLog('\r#REM','Motif %s is a subsequence of %s and has been removed.' % (newmotif.info['Name'],motif.info['Name']),screen=logrem,log=logrem)
                        return None

            ### Complete Formatting ###
            newmotif.makeVariants(msmode=self.opt['MSMS'],ambvar=not self.opt['Compare'])   # Makes self.list['Variants'] of length variants and basic self.dict['Search'] with no mismatches
            newmotif.misMatches(mismatch=self.dict['MisMatch'],msmode=self.opt['MSMS'])
                

            #!# Change Reversal process - move to method that reverses all motif at end and no output #!#
            
            ### Reverse Output ###
            if reverse and self.opt['Reverse']:
                REVMOT = open(self.info['Name'],'a')
                REVMOT.write('%s %s # Reversed by MotifList\n' % (newmotif.info['Name'],newmotif.info['Sequence']))
                REVMOT.close()

            ### Finish ###
            self.list['Motifs'].append(newmotif)
            #X#self.deBug(newmotif.info)
            #X#self.deBug(newmotif.list['Variants'])
            #X#self.deBug(newmotif.dict['Search'])
            return newmotif
        except:
            self.log.errorLog('Error in MotifList._addMotif()',quitchoice=True)     
            return None
#########################################################################################################################
    def removeMotif(self,Motif,remtxt=''):    ### Removes motif and occurrences from self.
        '''
        Removes motif and occurrences from self.
        >> Motif:Motif object to remove
        >> remtxt:str = Text to output to log
        '''
        try:
            ### Remove Motif ###
            if Motif in self.list['Motifs']:
                self.list['Motifs'].remove(Motif)
                if remtxt: self.log.printLog('\r#REM',remtxt)
            elif Motif:
                self.log.errorLog('Cannot remove motif "%s" - not in MotifList!' % Motif.info['Name'],printerror=False)
            else: return self.log.errorLog('Cannot remove motif "%s" - not in MotifList!' % Motif,printerror=False)

            ### Remove Occurrences ###
            for Occ in self.list['MotifOcc'][0:]:
                if Occ.obj['Motif'] == Motif: self.list['MotifOcc'].remove(Occ)
        except:
            self.log.errorLog('Error in MotifList.removeMotif()')     
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
        premotifs = self.motifs()[0:]
        newmotifs = rje_scoring.rankObj(self,self.motifs()[0:],stat,cutoff=cutoff,rev=rev)   ### Ranks objects using numerical data
        for Motif in premotifs:
            if Motif not in newmotifs: self.removeMotif(Motif)
        self.list['Motifs'] = newmotifs[0:]
#########################################################################################################################
    def patternStats(self,log=False):      ### Performs calculations all motifs based on basic pattern (info['Sequence'])
        '''Performs calculations all motifs based on basic pattern (info['Sequence']), adding to Motif.stat/info/opt.'''
        mx = 0.0
        for Motif in self.motifs():
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
        premotifs = self.motifs()[0:]
        newmotifs = rje_scoring.statFilterObj(self,self.motifs()[0:],statfilter)
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
    def combMotifOccStats(self,statlist=[],revlist=[],log=True,motiflist=[]):    ### Combines mean and percentile stats for the Occurrences of a Motif
        '''
        Combines mean and percentile stats for the Occurrences of a Motif.
        >> statlist:list of stats to combine from occurrences
        >> revlist:list of stats that should be ordered from low(best) to high(worst) rather than the other way round
        >> log:boolean [True] = whether to log progress of stat combination
        '''
        try:
            ### Setup Motif OccList ###
            motifocc = self.motifOcc(byseq=False,nested=False)
            mx = 0.0
            for Motif in motifocc.keys():
                if log: self.log.printLog('\r#COMB','Combined Occurence Stats: %.1f%%' % (mx/len(motifocc)),log=False,newline=False)
                mx += 100.0
                if motiflist and Motif not in motiflist: continue
                for stat in statlist:
                    ## Values for each occurrence ##
                    occval = []     # List of given stat for each occurrence
                    for Occ in motifocc[Motif]:
                        if stat == 'Cons':  # Special case
                            if Occ.stat['ALL_HOM'] > 0:
                                occval.append(Occ.stat['ALL_CONS'])
                        else:
                            try:
                                occval.append(Occ.stat[stat])
                            except:
                                print statlist
                                print Occ.stat
                                raise
                    if not occval:
                        continue
                    ## Basic Mean ##
                    try:
                        Motif.stat['%s_mean' % stat] = float(sum(occval)) / len(occval)
                    except:
                        print occval
                        raise
                    ## Percentiles ##
                    if self.stat['Percentile'] <= 0:
                        continue
                    if stat in revstat:
                        occval.sort(reverse=True)
                    else:
                        occval.sort()
                    pc = 100.0
                    ovx = len(occval)
                    for i in range(ovx):
                        ipc = 100.0 * float((ovx-1) - i) / (ovx-1)  # This position
                        jpc = 100.0 * float(ovx - i) / (ovx-1)      # Next position
                        if ipc == pc:   # Exact percentile
                            Motif.stat['%s_pc%d' % (stat,int(pc+0.5))] = occval[i]
                        elif ipc > pc and jpc < pc and (i+1) < ovx:
                            Motif.stat['%s_pc%d' % (stat,int(pc+0.5))] = float(occval[i] + occval[i+1]) / 2
                        elif (i+1) == ovx:
                            Motif.stat['%s_pc0' % stat] = occval[i]
                        else:
                            continue
                        pc -= self.stat['Percentile']
            if log:
                self.log.printLog('\r#COMB','Combined Occurence Stats complete.')
        except:
            self.log.errorLog('Major problem with rje_motiflist.combMotifOccStats()')
            raise
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
    ### <5> ### Motif conservation methods                                                                              #
#########################################################################################################################
    def setupGopher(self):  ### Sets up GOPHER directory etc.
        '''Sets up GOPHER directory etc.'''
        try:
            if not self.opt['Gopher']:
                return
            self.opt['UseAln'] = True
            if not self.info.has_key('GopherDir'):
                #X#self.deBug(string.split(self.info['AlnDir'],os.sep))
                #X#self.deBug(string.split(self.info['AlnDir'],os.sep)[:1])
                #X#self.deBug(string.split(self.info['AlnDir'],os.sep)[-2:-1])
                if self.info['AlnDir'][:3] == 'ALN' or string.split(self.info['AlnDir'],os.sep)[:1] == ['ALN']:   # Current Directory
                    self.info['GopherDir'] = rje.makePath('./',return_blank=False)
                    self.deBug('1. GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                elif string.split(self.info['AlnDir'],os.sep)[-1] == 'ALN':   # Gopher is one directory up
                    self.info['GopherDir'] = rje.makePath(string.join(string.split(self.info['AlnDir'],os.sep)[:-1],os.sep),return_blank=False)
                    self.deBug('2. GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                elif string.split(self.info['AlnDir'],os.sep)[-2:-1] == ['ALN']:
                    self.info['GopherDir'] = rje.makePath(string.join(string.split(self.info['AlnDir'],os.sep)[:-2],os.sep),return_blank=False)
                    self.deBug('3. GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                else:
                    self.info['GopherDir'] = self.info['AlnDir']
                    self.info['AlnDir'] = rje.makePath(self.info['AlnDir'] + 'ALN')
                    self.deBug('4. GopherDir = %s; AlnDir = %s' % (self.info['GopherDir'],self.info['AlnDir']))
                #X#self.deBug('GopherDir = %s' % self.info['GopherDir'])
                if self.info['GopherDir'] and not os.path.exists(self.info['GopherDir']):
                    rje.mkDir(self,self.info['GopherDir'])
                if not os.path.exists(self.info['AlnDir']):
                    rje.mkDir(self,self.info['AlnDir'])
                self.info['AlnExt'] = 'orthaln.fas'
        except:
            self.log.errorLog('Problem with MotifList.setupGopher(). UseAln cancelled')
            self.opt['UseAln'] = False
#########################################################################################################################
    def posMatrix(self):    ### Loads and builds PosMatrix for Conservation Scoring
        '''Loads and builds PosMatrix for Conservation Scoring.'''
        try:
            ### Setup PosMatrix ###
            self.dict['PosMatrix'] = {}

            ### Generate from property matrix if appropriate ###
            if self.opt['UseAln'] and self.info['ConScore'].lower() in ['all','prop']:
                propx = len(self.obj['AAPropMatrix'].prop)
                for a1 in self.list['Alphabet']:
                    for a2 in self.list['Alphabet']:
                        self.dict['PosMatrix']['%s%s' % (a1,a2)] = float(propx - self.obj['AAPropMatrix'].pdif['%s%s' % (a1,a2)]) / float(propx)
                self.dict['PropPosMatrix'] = self.dict.pop('PosMatrix')
                self.dict['PosMatrix'] = {}
                if self.info['ConScore'].lower() in ['prop']:
                    return

            ### Make default PosMatrix ###
            for aa in self.list['Alphabet']:
                self.dict['PosMatrix']['%s%s' % (aa,aa)] = 1.0

            ### Look for PosMatrix File ###
            posfile = self.info['PosMatrix']
            if posfile.lower() in ['','none']:
                return
            if not os.path.exists(posfile):
                self.log.errorLog('PosMatrix file "%s" not found!' % posfile,printerror=False)
                return

            ### Load Matrix ###
            plines = self.loadFromFile(posfile)
            _alphabet = plines[0].split()
            if len(_alphabet) > 1:  # Treat as matrix
                for a in _alphabet:
                    line = plines[_alphabet.index(a)+1].split()
                    if len(line) != (len(self.alphabet)+1):
                        self.log.errorLog('%s has wrong format! Does not match %s' % (line, _alphabet))
                        self.info['PosMatrix'] = 'None'
                        self.posMatrix()
                        raise ValueError
                    for i in range(1,len(line)):
                        score = float(line[i])
                        if score:
                            self.dict['PosMatrix']['%s%s' % (a,_alphabet[i-1])] = score
                self.log.printLog('#CONS','Position Scoring Matrix set from %s as matrix.' % posfile)
            else:   # Treat as equivalence files
                for line in plines:
                    match = rje.matchExp('^(\S+)',line)
                    if match:
                        for a1 in match[0]:
                            for a2 in match[0]:
                                self.dict['PosMatrix']['%s%s' % (a1,a2)] = 1.0
                self.log.printLog('#CONS','Position Scoring Matrix set from %s as equivalences.' % posfile)
        except:
            self.log.errorLog('Major problem in Presto.posMatrix(%s)' % self.info['PosMatrix'],quitchoice=True)
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
    def motifOut(self,filename='None',motlist=[]):      ### Outputs motifs in PRESTO format
        '''
        Outputs motifs in PRESTO format. #!# Check exactly how different formats are stored etc. #!#
        >> filename:str [None] = Name for output file. Will not output if '' or 'None'.
        >> motlist:list of Motif objects to output. If [], will use self.list['Motifs']
        '''
        try:
            ### Setup ###
            if not motlist:
                motlist = self.list['Motifs'][0:]
            if filename and filename.lower() != 'none':
                rje.backup(self,filename)
                OUT = open(filename,'a')
                for Motif in motlist:
                    OUT.write('%s  %s  # %s\n' % (Motif.info['Name'],Motif.info['Sequence'],Motif.info['Description']))
                OUT.close()
        except:
            self.log.errorLog('Problem during MotifList.motifOut()')
#########################################################################################################################
    def proteinAlignments(self,alndir='',hitname='AccNum'):    ### Generates copies of protein alignments, with motif hits marked.
        '''Generates copies of protein alignments, with motif hits marked.'''
        try:
            ### Setup ###
            if alndir and not os.path.exists(alndir):
                os.mkdir(alndir)
            ### Generate Alignments ###
            motif_occ = self.motifOcc(byseq=True)
            for seq in motif_occ.keys():
                occlist = []
                for Motif in motif_occ[seq].keys():
                    occlist += motif_occ[seq][Motif]
                self.singleProteinAlignment(seq,occlist,alndir,hitname,usegopher=True)                
        except:            
            self.log.errorLog('Major problem with SlimPicker.proteinAlignments',quitchoice=True)
#########################################################################################################################
    def singleProteinAlignment(self,seq,occlist,alndir='',hitname='AccNum',usegopher=True,savefasta=True):    ### Generates copies of protein alignments, with motif hits marked.
        '''
        Generates copies of protein alignment, with motif hits marked. Return SeqList object.
        >> seq:Sequence object for alignment
        >> occlist:list of MotifOcc objects for sequence
        >> alndir:str = Alignment directory for output
        >> hitname:str = Format of hitnames, used for naming files
        >> usegopher:boolean = whether to look to use Gopher if settings correct (set this False if done already)
        >> savefasta:boolean [True] = whether to save fasta file or simply return SeqList object alone.
        '''
        try:
            ### ~ Try to find alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Alignment'
            if self.opt['UseAln']: aln = rje_motif_stats.loadOrthAln(self,seq,usegopher=usegopher)
            if not self.opt['UseAln'] or not aln:
                alncmd = ['seqin=None','query=%s' % seq.shortName(),'accnr=F','seqnr=F','autofilter=F','align=T','gnspacc=F'] 
                aln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
                aln.seq = [seq]
                aln.obj['QuerySeq'] = seq

            ### ~ Check Query and File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Check SeqList'
            if not aln.obj['QuerySeq'] and not aln.querySeq(query=seq.info['AccNum']):
                self.log.printLog('#ERR','Problem finding %s in %s.' % (self.info['Name'],file))
                return None
            qry = aln.obj['QuerySeq']

            ### ~ Create new sequence containing motifs, mapped onto Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'MotifSeq'
            motifseq = ['-'] * qry.seqLen()
            for Occ in occlist:
                #x#print Occ.info, Occ.stat
                #x#print qry.info['Sequence']
                #x#print qry.info['MaskSeq']
                (start,end) = rje_motif_stats.findOccPos(self,Occ,qry)
                i = 0
                for a in qry.info['Sequence'][start:end]:
                    if a == '-':
                        continue
                    motifseq[start+i] = Occ.info['Variant'][i]
                    i += 1
            aln._addSeq('Motifs',string.join(motifseq,''))
            aln.seq = aln.seq[-1:] + aln.seq[:-1]
            if savefasta:
                if qry.info.has_key(hitname):
                    aln.saveFasta(seqfile='%s%s.proteinaln.fas' % (alndir,qry.info[hitname]),log=False)
                else:
                    aln.saveFasta(seqfile='%s%s.proteinaln.fas' % (alndir,qry.shortName()),log=False)
            return aln

        except:            
            self.log.errorLog('Major problem with MotifList.proteinAlignments(%s)' % _stage)
            return None
#########################################################################################################################
    def motifAlignments(self,resfile='motifaln.fas'):      ### Makes motif alignments from occurrences
        '''
        Makes motif alignments from occurrences. MotifOcc objects should have Sequence objects associated with them. If
        necessary, add a method to go through and generate Sequence objects using MotifOcc.info['FastaCmd'].
        >> resfile:str = Name of output file
        '''
        try:
            ### Long ###
            if self.stat['XDivide'] < 1:
                motif_occ = self.motifOcc()
                for Motif in self.motifs():     #!# Add log output #!#
                    self.motifAlnLong(Motif,motif_occ[Motif],append=True,memsaver=False,resfile=resfile)
            
            ### Setup ###
            _stage = 'Setup'
            motif_occ = self.motifOcc(byseq=True,justdata='Pos')
            for Seq in motif_occ.keys():
                for Motif in motif_occ[Seq].keys():
                    newpos = []
                    for pos in motif_occ[Seq][Motif]:
                        newpos.append(int(pos) - 1)
                    motif_occ[Seq][Motif] = newpos
            ltxt = 'Constructing motif alignments: %s motifs & %s seqs.' % (rje.integerString(self.motifNum()),
                                                                            rje.integerString(len(motif_occ)))
            self.log.printLog('#ALN',ltxt,log=False,newline=False)
            extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
            extract_seq.opt['ReplaceChar'] = False

            ### Setup positions in sequences ###
            ## Work off motif_occ = Dictionary of {Sequence:{Motif:[positions]}} (Pos is from 0 to L-1)
            _stage = 'SeqPos'
            max_occ = {}            # Dictionary of Sequence:Max occurrences for any motif
            max_pos = 0             # Max motif position
            occ_seq = {'Motif':{}}  # Dictionary of {Sequence:{Motif:[aligned sequence fragments]}}
            motiflist = []          # List of motifs with 1+ occurrences
            for motif in self.list['Motifs']:
                occ_seq['Motif'][motif] = []
            for seq in motif_occ.keys():
                occ_seq[seq] = {}
                max_occ[seq] = 0
                for motif in self.list['Motifs']:
                    occ_seq[seq][motif] = []
                    if motif not in motif_occ[seq].keys():
                        continue    #!# Improve way missing sequences/motifs are dealt with #!#
                    if len(motif_occ[seq][motif]) > max_occ[seq]:
                        max_occ[seq] = len(motif_occ[seq][motif])
                    if len(motif_occ[seq][motif]) > 0 and motif not in motiflist:
                        motiflist.append(motif)
                    for pos in motif_occ[seq][motif]:
                        if pos >= max_pos:
                            max_pos = pos + 1
            self.log.printLog('\r#ALN','%s: setup complete.' % ltxt,log=False)
            neworder = self.list['Motifs'][0:]
            for Motif in self.list['Motifs'][0:]:
                if Motif not in motiflist: neworder.remove(Motif)
            motiflist = neworder[0:]

            ### Make sequences ###
            _stage = 'Sequences'
            for motif in motiflist:
                pattern = motif.info['Sequence']
                self.log.printLog('\r#ALN','%s: %.f%%.' % (ltxt,100.0*self.list['Motifs'].index(motif)/len(self.list['Motifs'])),log=False,newline=False)
                ## Pattern info ##
                patseq = '-%s-' % string.split(motif.info['Name'])[0]   # Name of motif
                while len(patseq) < (len(rje.preZero(max_pos,max_pos)) + 2):
                    patseq += '-'
                namelen = len(patseq)
                overlap = len(pattern) - motif.stat['FullLength']       # Extra length of pattern vs. longest occurrence
                #x#print pattern, motif.stat['FullLength'], overlap, int(overlap/2), int((overlap+1)/2)
                patseq += '-' * self.stat['FlankSize']
                patseq += pattern       
                patseq += '-' * (self.stat['FlankSize'] - overlap)
                occ_seq['Motif'][motif].append(patseq[0:])
                ## Occurrences ##
                for seq in motif_occ.keys():
                    if not motif_occ[seq].has_key(motif):
                        motif_occ[seq][motif] = {}
                    for x in range(max_occ[seq]):  # Need an entry for each potential occurrence
                        if x < len(motif_occ[seq][motif]):  # Actual entry
                            r = motif_occ[seq][motif][x]
                            patseq = '-%s-' % rje.preZero(r+1,max_pos)
                            while len(patseq) < namelen:
                                patseq += '-'
                            (left,right) = (r-self.stat['FlankSize'],r+motif.stat['FullLength']-1+self.stat['FlankSize'])
                            if left < 0:
                                patseq += '-' * -left + seq.info['Sequence'][:right]
                            else:
                                patseq += seq.info['Sequence'][left:right]
                            patseq += '-' * (len(occ_seq['Motif'][motif][0]) - len(patseq))
                        else:   # Add blank one
                            patseq = '-' * len(occ_seq['Motif'][motif][0])
                        occ_seq[seq][motif].append(patseq[0:])
            self.log.printLog('\r#ALN','%s: 100.0%%.' % (ltxt))

            ### Output sequence file ###
            _stage = 'Output'
            ## Motif Data ##
            name = '%s motifs' % self.info['Name']
            outlist = []
            for motif in motiflist:
                pattern = motif.info['Sequence']
                outlist.append(occ_seq['Motif'][motif][0])
            extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ## Occurrences ##
            for seq in motif_occ.keys():
                name = seq.info['Name']
                for x in range(max_occ[seq]):
                    outlist = []
                    for motif in motiflist:
                        pattern = motif.info['Sequence']
                        outlist.append(occ_seq[seq][motif][x])
                    extract_seq._addSeq(name,string.join(outlist,sep='X' * self.stat['XDivide']))
            ### Output ###
            extract_seq.info['Name'] = resfile
            extract_seq.saveFasta()

        except:            
            self.log.errorLog('Major problem with MotifList.motifAlignments(%s)' % _stage,quitchoice=True)
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
        try:
            ### Setup ###
            if not resfile:
                resfile = Motif.info['Name'] + '.motifaln.fas'
            extract_seq = rje_seq.SeqList(self.log,['logrem=F']+self.cmd_list+['seqin=None','autofilter=T','newlog=F'])
            extract_seq.opt['ReplaceChar'] = False

            ### Motif Setup ###
            max_pos = 10000
            patseq = '-' * (len(rje.preZero(max_pos,max_pos)) + 2)
            overlap = len(Motif.info['Sequence']) - Motif.stat['FullLength']       # Extra length of pattern vs. longest occurrence
            patseq += '-' * self.stat['FlankSize']
            patseq += Motif.info['Sequence']
            patseq += '-' * (self.stat['FlankSize'] - overlap)
            if memsaver or not append:      ### Output Motif itself
                extract_seq._addSeq(Motif.info['Name'],patseq)
            motlen = len(patseq)
                
            ### Occurrences ###
            for Seq in seq_occ:
                for Occ in seq_occ[Seq]:
                    if Occ.obj['Motif'] != Motif:
                        continue
                    patseq = '-%s-' % rje.preZero(Occ.stat['Pos'],max_pos)
                    r = Occ.stat['Pos'] - 1
                    (left,right) = (r-self.stat['FlankSize'],r+Motif.stat['FullLength']-1+self.stat['FlankSize'])
                    if left < 0:
                        patseq += '-' * -left + Seq.info['Sequence'][:right]
                    else:
                        patseq += Seq.info['Sequence'][left:right]
                    patseq += '-' * (motlen - len(patseq))
                    extract_seq._addSeq(Seq.shortName(),patseq)

            ### Output ###
            extract_seq.saveFasta(seqfile=resfile,append=append,log=False)
        except:            
            self.log.errorLog('Major problem with MotifList.MotifAlnLong',quitchoice=True)            
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
    def motifInfo(self,expfile=None,statlist=[]):   ### Produces summary table for motifs, including expected values and information content.
        '''
        Produces summary table for motifs, including expected values and information content.
        >> expfile:str = Filename from which expectations have been calculated and should be output
        >> statlist:list of extra statistics to output
        '''
        try:
            #!# Check these things! #!#

            #!# Either cut-down on output, else reduce what it does and make sure objects have relevant data in them #!#

            
            ### Setup Results ###
            outfile = self.info['MotInfo']
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            rje.backup(self,outfile)    # Backs up if existing and not self.opt['Append']
            headers = ['Motif','Pattern','Description','MaxLength','MinLength','FixLength','FullLength']
            if self.opt['Expect'] and not self.opt['Compare']:
                headers.append('Expect')
            if not self.opt.has_key('MotifIC') or self.opt['MotifIC']:
                headers.append('IC')
            if self.opt.has_key('Searched') and self.opt['Searched']:
                headers += ['OccNum','OccSeq']
            headers += statlist
            rje.delimitedFileOutput(self,outfile,headers,delimit,datadict={})   # Output headers

            ### Output ###
            for motif in self.list['Motifs']:
                datadict = {'Motif':motif.info['Name'],
                            'Pattern':motif.info['Sequence'],
                            'Description':motif.info['Description'],
                            'IC':'%.3f' % motif.stat['IC'],
                            'MaxLength':motif.stat['Length']
                            }
                if expfile:
                    datadict['Expect'] = rje_motif.expectString(motif.dict['Expect'][expfile])
                for s in ['MinLength','FixLength','FullLength','OccNum','OccSeq'] + statlist:
                    datadict[s] = motif.getData(s)
                rje.delimitedFileOutput(self,outfile,headers,delimit,datadict)
            
        except:
            self.log.errorLog('Error in presto.motifInfo()',quitchoice=True)
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
### End of SECTION II: MotifList                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:        
        presto = Presto(log=mainlog,cmd_list=cmd_list)
        presto.run()

    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
