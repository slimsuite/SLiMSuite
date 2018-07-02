#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMMaker
Description:  SLiM generator from aligned peptide sequences
Version:      1.7.0
Last Edit:    01/02/17
Citation:     Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]
Webserver:    http://www.slimsuite.unsw.edu.au/servers/slimmaker.php
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program has a fairly simple function of reading in a set of sequences and generating a regular expression motif
    from them. It is designed with protein sequences in mind but should work for DNA sequences too. Input sequences can
    be in fasta format or just plain text (with no sequence headers) and should be aligned already. If varlength=F then
    gapped positions will be ignored (treated as Xs) and variable length wildcards are not returned. If varlength=T, any
    gapped positions will be assessed based on the ungapped peptides at that position and a variable length inserted.
    This variable-length position may be a wildcard or it may be a defined position if there is sufficient signal in the
    peptides with amino acids at that position.

    SLiMMaker considers each column of the input in turn and compresses it into a regular expression element according to
    some simple rules, screening out rare amino acids and converting particularly degenerate positions into wildcards.
    Each amino acid in the column that occurs at least X times (as defined by minseq=X) is considered for the regular
    expression definition for that position. The full set of amino acids meeting this criterion is then assessed for
    whether to keep it as a defined position, or convert into a wildcard. First, if the number of different amino acids
    meeting this criterion is zero or above a second threshold (maxaa=X), the position is defined as a wildcard. Second,
    the proportion of input sequences matching the amino acid set is compared to a minimum frequency criterion
    (minfreq=X). Failing to meet this minimum frequency will again result in a wildcard. Otherwise, the amino acid set is
    added to the SLiM definition as either a fixed position (if only one amino acid met the `minseq` criterion) or as a
    degenerate position. Finally, leading and trailing wildcards are removed.

    By default, each defined position in a motif will contain amino acids that (a) occur in at least three sequences
    each, (b) have a combined frequency of >=75%, and (c) have 5 or fewer different amino acids (that occur in 3+
    sequences). The same `minseq=X` threshold is also used to determine whether flexible length *defined* positions are
    generated (if `varlength=T`), i.e. to have a flexible-length non-wildcard position, at least minseq sequences must
    have a gap at that position. This does not apply to flexible-length wildcards.

    Note. Unless the "iterate" function is used, the final motif only contains defined positions that match a given 
    frequency of the input (75% by default). Because positions are considered independently, however, the final motif
    might occur in fewer than 75% of the input sequences. SLiMSearch can be used to check the occurrence stats.

    Version 1.5.0 incorporates a new peptide alignment mode to deal with unaligned peptides. This is controlled by the
    `peptalign=T/F/X` option, which is set to True by default. If given a regular expression, this will be used to guide
    the alignment. Otherwise, the longest peptides will be used as a guide and the minimum number of gaps added to
    shorter peptides. PeptCluster peptide distance measures are used to assess different variants, starting with simple
    sequence identity, then amino acid properties (if ties) and finally PAM distances. One of the latter can be set as
    the priority using `peptdis=X`. Peptide alignment assumes that peptides have termini (^ & $) or flanking wildcards
    added. If not, set `termini=F`.

    Version 1.6.0 added the option to incorporate amino acid equivalencies to extend motif sites beyond the top X% of
    amino acids. This works by identifying a degenerate set of amino acids as normal using `minseq=X` and then checking
    whether these form a subset of an equivalence group prior to the `minfreq=X` filter. If so, it will try extending the
    degenerate position to incorporate additional members of the equivalence group. For example, `IL` could incorporate
    additional `MVF` amino acids of an `FILMV` group. Only amino acids represented in the peptides will be added. Single
    amino acids will also be extended, e.g. `S` could be extended to `ST`. This mode is switched on with `extendaa=T`.
    The `equiv=LIST` option sets the equivalence groups.

    If two or more equivalence groups could be extended, the one with the most members will be chosen. If tied, the one
    with fewest possible amino acids (from `equiv=LIST`) will be chosen. If still tied, the first group in the list will
    take precedence.

Commandline:
    ### ~ SLiMMaker Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    peptides=LIST   : These can be entered as a list or a file. If a file, lines following '#' or '>' are ignored
    maxlen=INT      : Maximum length for peptide [50]
    peptalign=T/F/X : Align peptides. Will use as guide regular expression, else T/True for regex-free alignment. [True]
    minseq=X        : Min. no. of sequences for an aa to be in [3]
    minfreq=X       : Min. combined freq of accepted aa to avoid wildcard [0.75]
    maxaa=X         : Max. no. different amino acids for one position [5]
    ignore=X        : Amino acid(s) to ignore. (If nucleotide, would be N-) ['X-']
    dna=T/F         : Whether "peptides" are actually DNA fragments [False]
    iterate=T/F     : Whether to perform iterative SLiMMaker, re-running on matched peptides with each iteration [False]
    varlength=T/F   : Whether to identifies gaps in aligned peptides and generate variable length motif [True]
    extendaa=T/F    : Whether to extend ambiguous aa using equivalence list [False]
    equiv=LIST      : List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,KRH,DE,ST]

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_slim, rje_zen
import peptcluster
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial Working Version. Some minor modifications for SLiMBench including iterative SLiMMaker.
    # 1.1 - Modified to work with end of line characters.
    # 1.2.0 - Modified to work with REST servers.
    # 1.3.0 - Added varlength option to identify gaps in aligned peptides and generate variable length motif.
    # 1.3.1 - Fixed varlength option to work with end of peptide gaps. (Gaps ignored completely - should not be there!)
    # 1.4.0 - Add iteration REST output.
    # 1.4.1 - Add unmatched peptides REST output.
    # 1.4.2 - Fixed bug with variable length wildcards at start of sequence.
    # 1.5.0 - Added peptalign=X functionality, using PeptCluster peptide alignment.
    # 1.6.0 - Added equiv=LIST : List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,FYH,KRH,DE,ST]
    # 1.6.1 - Fixed peptide case bug.
    # 1.7.0 - Added maxlen parameter.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add a calculation of how many input sequences are matched by the motif.
    # [Y] : Add an iterative SLiMMaker mode.
    # [ ] : Add normalisation for UPC relationships.
    # [Y] : Add peptide alignment using PeptCluster.
    # [ ] : |-- Deal with Regex variants? (Winner takes all? Ignore?)
    # [ ] : Add equivalence file expansion of positions to include similar AA.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('SLiMMaker', '1.7.0', 'January 2017', '2012')
    description = 'SLiM generator from aligned peptide sequences'
    author = 'Dr Richard J. Edwards.'
    comments = ['Cite: Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]',
                'Please report bugs to Richard.Edwards@UNSW.edu.au']
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
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
### SECTION II: SLiMMaker Class                                                                                         #
#########################################################################################################################
class SLiMMaker(rje_obj.RJE_Object):     
    '''
    SLiMMaker Class. Author: Rich Edwards (2012).

    Str:str
    - Ignore = Amino acid(s) to ignore. (If nucleotide, would be N-) ['X-']
    - PeptAlign = Align peptides. Will use as guide regular expression, else T/True for regex-free alignment. []

    Bool:boolean
    - DNA = Whether "peptides" are actually DNA fragments [False]
    - ExtendAA = Whether to extend ambiguous aa using equivalence list [False]
    - Iterate = Whether to perform iterative SLiMMaker, re-running on matched peptides with each iteration [False]
    - VarLength=T/F   : Whether to identifies gaps in aligned peptides and generate variable length motif [True]

    Int:integer
    - MinSeq = Min. no. of sequences for an aa to be in [3]
    - MaxAA = Max. no. different amino acids for one position [5]
    - MaxLen=INT      : Maximum length for peptide [50]

    Num:float
    - MinFreq = Min. combined freq of accepted aa to avoid wildcard [0.75]
    
    List:list
    - Equiv = List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,KRH,DE,ST]
    - Input  = List of original input peptides
    - Peptides  = These can be entered as a list or a file. If a file, lines following '#' or '>' are ignored

    Dict:dictionary    

    Obj:RJE_Objects
    - PeptCluster = peptcluster.PeptCluster object for peptide distances
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Ignore','PeptAlign']
        self.boollist = ['DNA','ExtendAA','Iterate','VarLength']
        self.intlist = ['MinSeq','MaxAA','MaxLen']
        self.numlist = ['MinFreq']
        self.listlist = ['Equiv','Input','Peptides']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'Ignore':'X-','Basefile':'slimmaker','PeptAlign':'True'})
        self.setBool({'VarLength':True})
        self.setInt({'MinSeq':3,'MaxAA':5,'MaxLen':50})
        self.setNum({'MinFreq':0.75})
        self.list['Equiv'] = string.split('AGS,ILMVF,FYW,KRH,DE,ST',',')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['PeptCluster'] = peptcluster.PeptCluster(self.log,['peptdis=None']+self.cmd_list)
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
                self._cmdReadList(cmd,'str',['Ignore','PeptAlign'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['DNA','ExtendAA','Iterate','VarLength'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinSeq','MaxAA','MaxLen'])   # Integers
                self._cmdReadList(cmd,'float',['MinFreq']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Equiv','Peptides'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getBool('DNA'): self.str['Ignore'] = string.replace(self.str['Ignore'].upper(),'X','N')
        if self.getStrUC('PeptAlign') in ['F','FALSE']: self.setStr({'PeptAlign':''})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,iterate=None,log=True):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slim = ''
            if iterate == None: iterate = self.getBool('Iterate')
            elif iterate: self.setBool({'Iterate':True})
            if not self.setup(log=log): return ('','SLiMMaker setup failed. Check log.')
            if not self.list['Input']: self.list['Input'] = self.list['Peptides'][0:]
            equiv = []
            if self.getBool('ExtendAA'):
                #self.warnLog('Equivalence mode (extendaa=T) not yet implemented! Please contact author.')
                self.printLog('#EQUIV','[%s]' % string.join(self.list['Equiv'],'] ['))
                equiv = self.list['Equiv'][0:]
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slim = rje_slim.makeSlim(self.list['Peptides'],self.getInt('MinSeq'),self.getNum('MinFreq'),self.getInt('MaxAA'),self,self.getStr('Ignore'),self.getBool('VarLength'),equiv)
            self.dict['Output']['slim'] = slim
            if log: self.printLog('#SLIM','SLiM generated: "%s"' % slim)
            if not slim: return (slim,'Unable to make a SLiM with these settings and peptides')
            ## ~ [2a] ~ Assess matches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            matched = []
            if self.getBool('DNA'): regexp = string.replace(slim,'N','.')
            else: regexp = string.replace(slim,'X','.')
            for peptide in self.list['Peptides']:
                searchpep = string.replace('X%sX' % peptide,'$X','')
                searchpep = string.replace(searchpep,'X^','')
                searchpep = string.replace(searchpep,'-','')
                try:
                    if rje.matchExp('(%s)' % regexp,searchpep): matched.append(peptide)
                except: self.errorLog('Error with SLiM/peptide match, %s vs %s' % (regexp,searchpep))
            sx = len(matched)
            matchstr = 'SLiM matches %d of %d sequences (%.1f%%).' % (sx,len(self.list['Peptides']),(100.0*sx)/len(self.list['Peptides']))
            if log: self.printLog('#FREQ',matchstr)
            if iterate:
                self.dict['Output']['iterate'] += '%s: %s\n' % (slim,matchstr)
                self.dict['Output']['iterate'] += '-> %s\n' % string.join(matched,',')
            if iterate and (len(matched) != len(self.list['Peptides'])):
                if not matched: return (slim,'Unable to make an interative SLiM with these settings and peptides')
                if self.getStrLC('PeptAlign'):
                    self.list['Peptides'] = string.split(string.replace(string.join(matched),'-',''))
                else: self.list['Peptides'] = matched
                return self.run(iterate=True)
            ### ~ [3] REST Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if iterate:
                (matchstr,matched) = self.inputMatches(regexp)
                if log: self.printLog('#FREQ',matchstr)
            self.dict['Output']['match'] = matchstr
            self.dict['Output']['matches'] = string.join(matched,'\n')
            try:
                unmatched = self.list['Input'][0:]
                for pep in matched: unmatched.remove(pep)
                self.dict['Output']['unmatched'] = string.join(unmatched,'\n')
            except: self.dict['Output']['unmatched'] = self.errorLog('SLiMMaker Umatched Error')
            return (slim,matchstr)
        except: return (slim,self.errorLog('SLiMMaker Error'))
#########################################################################################################################
    def inputMatches(self,regexp):  ### Returns the matches for the original peptides
        '''Returns the matches for the original peptides.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            matched = []
            for peptide in self.list['Input']:
                searchpep = string.replace('X%sX' % peptide,'$X','')
                searchpep = string.replace(searchpep,'X^','')
                searchpep = string.replace(searchpep,'-','')
                try:
                    if rje.matchExp('(%s)' % regexp,searchpep): matched.append(peptide)
                except: self.errorLog('Error with SLiM/peptide match, %s vs %s' % (regexp,searchpep))
            sx = len(matched)
            matchstr = 'SLiM matches %d of %d input sequences (%.1f%%).' % (sx,len(self.list['Input']),(100.0*sx)/len(self.list['Input']))
        except: self.errorLog('Error with inputMatches()'); matchstr = 'Error with inputMatches()'; matched = []
        return (matchstr,matched)
#########################################################################################################################
    def setup(self,log=True):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Peptides'] = rje.listUpper(self.list['Peptides'])
            mypep = []; peplen = 0; pepaln = True
            for inpep in self.list['Peptides']:
                pep = string.split(inpep,'>')[0]
                pep = string.split(pep,'#')[0]
                for pep in string.split(pep):
                    if len(pep) > self.getInt('MaxLen'):
                        if log: self.warnLog('Peptide ignored - exceeds maxlen (%d): %s' % (self.getInt('MaxLen'),pep))
                    elif pep:
                        mypep.append(pep)
                        if peplen and len(pep) != peplen: pepaln = False
                        else: peplen = len(pep)
            inx = len(self.list['Peptides'])
            self.list['Peptides'] = mypep
            if 'peptides' not in self.dict['Output']: self.dict['Output']['peptides'] = string.join(self.list['Peptides'],'\n')
            peptxt = {True:'fragments',False:'peptides'}[self.getBool('DNA')]
            if log: self.printLog('#PEP','%s %s read from %s input %s' % (len(mypep),peptxt,inx,peptxt))
            self.setInt({'MinSeq':max(1,self.getInt('MinSeq'))})
            if len(self.list['Peptides']) < self.getInt('MinSeq'):
                if log: self.warnLog('Too few peptides (%d) for minseq=%d' % (len(self.list['Peptides']),self.getInt('MinSeq')))
                return False
            if not pepaln and self.getStrLC('PeptAlign'):
                if log: self.printLog('#ALN','Unaligned peptides: aligning with PeptCluster')
                pclust = self.obj['PeptCluster']
                self.list['Peptides'] = pclust.peptAlign(self.getStrUC('PeptAlign'),self.list['Peptides'])
                self.dict['Output']['aligned'] = string.join(self.list['Peptides'],'\n')
            elif not pepaln: self.printLog('#WARN','Warning: %s are not all same length!' % peptxt)
            self.restSetup()
            return True
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### REST Output Methods                                                                                     #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        ### Running SLiMMaker:
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        ### Available REST Outputs:
        slim = Short Linear Motif pattern returned
        match = Number of input peptides matched by the SLiM
        peptides = Original input peptides
        aligned = aligned peptides
        matches = Peptides matching the SLiM
        unmatched = Peptides not matching the SLiM
        iterate = SLiM/Peptide iterations. [&iterate=T]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Iterate') and 'iterate' not in self.dict['Output']: self.dict['Output']['iterate'] = ''
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['slim','match','peptides','aligned','matches','unmatched','iterate']
#########################################################################################################################
### End of SECTION II: SLiMMaker Class                                                                                  #
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
    try: SLiMMaker(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
    #mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
