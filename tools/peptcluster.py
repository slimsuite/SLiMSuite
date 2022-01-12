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
Module:       PeptCluster
Description:  Peptide Clustering Module
Version:      1.5.4
Last Edit:    01/11/17
Webserver:    http://www.slimsuite.unsw.edu.au/servers/peptcluster.php
Copyright (C) 2012  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program is for simple sequence-based clustering of short (aligned) peptide sequences. First, a pairwise distance
    matrix is generated from the peptides. This distance matrix is then used to generate a tree using a distance method
    such as Neighbour-Joining or UPGMA.

    Default distances are amino acid property differences loaded from an amino acid property matrix file.

    Version 1.5.0 incorporates a new peptide alignment mode to deal with unaligned peptides. This is controlled by the
    `peptalign=T/F/X` option, which is set to True by default. If given a regular expression, this will be used to guide
    the alignment. Otherwise, the longest peptides will be used as a guide and the minimum number of gaps added to
    shorter peptides. Pairwise peptide distance measures are used to assess different variants, starting with amino acid
    properties, then simple sequence identity (if ties) and finally PAM distances. One of the latter can be set as
    the priority using `peptdis=X`. Peptide alignment assumes that peptides have termini (^ & $) or flanking wildcards
    added. If not, set `termini=F`.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    peptides=LIST   : These can be entered as a list or file. Lines following '#' or '>' ignored [peptides.txt]
    aaprop=FILE     : File of amino acid properties [aaprop.txt]
    aadis=FILE      : Alternative amino acid distance matrix [None]
    peptalign=T/F/X : Align peptides. Will use as guide regular expression, else T/True for regex-free alignment. []
    termini=T/F     : Whether peptides for alignment have termini (^ & $) or X flanking regex match [True]
    maxgapvar=X     : Maximum number of consecutive gaps to allow for peptide alignment without Regex guide [3]
    maxgapx=X       : Maximum total number of gaps to allow for peptide alignment without Regex guide [5]

    ### ~ CLUSTER OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    peptdis=X       : Method for generating pairwise peptide distances (id/prop/pam) [prop]
    peptcluster=X   : Clustering mode (upgma/wpgma/neighbor) [upgma]
    
    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    savedis=FILE    : Output distance matrix to file [peptides.*]
    outmatrix=X     : Type for output matrix - tdt / csv / mysql / phylip / png [tdt]
    savetree=FILE   : Save generated tree as FILE [peptides.peptdis.peptcluster.*]
    treeformats=LIST: List of output formats for generated trees (nsf/nwk/text/r/png/bud/qspec/cairo/te/svg/html) [nwk,text]
    
See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_aaprop, rje_obj, rje_pam, rje_slim, rje_tree, rje_zen
import rje_dismatrix_V2 as rje_dismatrix
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial working version.
    # 1.1 - Modified output a little.
    # 1.2 - CGI Output tuple.
    # 1.3 - PAM clustering and WPGMA.
    # 1.4 - Bug fixes for end of sequence characters and different length peptides.
    # 1.5.0 - Added peptalign=T/F/X function for aligning peptides using regex or minimal gap addition. Added REST.
    # 1.5.1 - Updated REST output. Removed peptide redundancy.
    # 1.5.2 - Improved clarity of warning message.
    # 1.5.3 - Added catching of splitPattern() error during alignment.
    # 1.5.4 - Fixed error message error.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Generate basic functional version using simple property distances and UPGMA.
    # [Y] : Add PAM-based distance.
    # [Y] : Add phylip tree methods.
    # [Y] : Add a method to return a tuple for a CGI script.
    # [ ] : Add BLOSUM-based distance method.
    # [ ] : Add proper handling of alphabets and gaps.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('PeptCluster', '1.5.4', 'November 2017', '2012')
    description = 'Peptide Clustering Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
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
### SECTION II: PeptCluster Class                                                                                       #
#########################################################################################################################
class PeptCluster(rje_obj.RJE_Object):     
    '''
    PeptCluster Class. Author: Rich Edwards (2012).

    Str:str
    - AAProp = File of amino acid properties [aaprop.txt]
    - AADis = Alternative amino acid distance matrix [None]
    - PeptAlign = Align peptides. Will use as guide regular expression, else T/True for regex-free alignment. []
    - PeptDis = Method for generating pairwise peptide distances (id/prop/pam) [prop]
    - PeptCluster = Clustering mode (upgma/wpgma/neighbor) [upgma]
    - Peptides = Name of file containing peptides (read in as list but used for basefile) [peptides.txt]
    - SaveDis = Output distance matrix to file [peptides.*]
    - OutMatrix = Type for output matrix - tdt / csv / mysql / phylip / png [tdt]
    - SaveTree = Save generated tree as FILE [peptides.peptdis.peptcluster.*]
    
    Bool:boolean
    - Termini = Whether peptides for alignment have termini (^ & $) or X flanking regex match [True]

    Int:integer
    - MaxGapVar = Maximum number of consecutive gaps to allow for peptide alignment [3]
    - MaxGapX = Maximum total number of gaps to allow for peptide alignment without Regex guide [5]

    Num:float
    
    List:list
    - Peptides = These can be entered as a list or a file. If a file, lines following '#' or '>' are ignored [peptides.txt]
    - TreeFormats = List of output formats for generated trees (nsf/nwk/text/r/png/bud/qspec/cairo/te/svg/html) [nwk,text]

    Dict:dictionary

    Obj:RJE_Objects
    - AADis = Distance Matrix containing pairwise AA distances
    - AAProp = AA Property object for making AA distances from properties
    - PAM = PAM Control object
    - PeptDis = Distance Matrix containing pairwise peptide distances
    - Tree = rje_tree.Tree object containing clustered peptides.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['AAProp','AADis','PeptAlign','PeptDis','PeptCluster','Peptides','SaveDis','OutMatrix','SaveTree']
        self.boollist = ['Termini']
        self.intlist = ['MaxGapVar','MaxGapX']
        self.numlist = []
        self.listlist = ['Peptides','TreeFormats']
        self.dictlist = []
        self.objlist = ['AADis','AAProp','PAM','PeptDis','Tree']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'AAProp':'aaprop.txt','AADis':'None','PeptDis':'prop','PeptCluster':'upgma','SaveDis':'',
                     'OutMatrix':'tdt','SaveTree':'','Peptides':'peptides.txt'})
        self.setBool({'Termini':True})
        self.setInt({'MaxGapVar':3,'MaxGapX':5})
        self.setNum({})
        self.list['Peptides'] = []
        self.list['TreeFormats'] = ['text','nwk']
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
                self._cmdReadList(cmd.lower(),'str',['PeptAlign','PeptDis','PeptCluster','OutMatrix'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['AAProp','AADis','Peptides','SaveDis','SaveTree'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['Termini'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxGapVar','MaxGapX'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Peptides','TreeFormats'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if not self.list['Peptides']: self._cmdRead('peptides=peptides.txt','list','Peptides')
        if self.getStrUC('PeptAlign') in ['F','FALSE']: self.setStr({'PeptAlign':''})
        self.list['TreeFormats'] = rje.listLower(self.list['TreeFormats'])
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Generate Peptide Distances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Optional peptide alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Peptides'] = self.peptAlign(save=True)
            self.dict['Output']['aligned'] = string.join(self.list['Peptides'],'\n')
            if not self.getStrLC('PeptDis'):
                self.printLog('#END','No peptide distance method: clustering cancelled.')
                return
            ## ~ [2b] ~ Generate Peptide Distances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.preparePeptides()
            for pep1 in self.list['Peptides']:
                for pep2 in self.list['Peptides']: self.peptDist(pep1,pep2)
            self.dict['Output']['matrix'] = self.str['SaveDis']
            self.obj['PeptDis'].saveMatrix(filename=self.str['SaveDis'])
            ### ~ [3] ~ Perform Peptide Clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.peptCluster()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [0] Setup File names etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('SaveDis').lower() in ['','none']:
                base = 'peptides'
                if rje.checkForFile(self.getStr('Peptides')): base = rje.baseFile(self.getStr('Peptides'))
                if self.baseFile().lower() not in ['','none']: base = self.baseFile()
                self.baseFile(base)
                self.setStr({'SaveDis':'%s.%s.%s' % (base,self.getStr('PeptDis'),self.getStr('PeptCluster'))})
            if self.getStr('OutMatrix') in ['tdt','csv','png','phylip']: self.str['SaveDis'] += '.%s' % self.getStr('OutMatrix')[:3]
            else: self.str['SaveDis'] += '.txt'
            self.dict['Output']['peptides'] = string.join(self.list['Peptides'],'\n')
            ### ~ [1] Setup Distance Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['AADis'] = rje_dismatrix.DisMatrix(self.log,['nsf2nwk=T']+self.cmd_list)
            self.obj['AADis'].info['Name'] = 'Pairwise AA distances'
            self.obj['PeptDis'] = rje_dismatrix.DisMatrix(self.log,['nsf2nwk=T']+self.cmd_list)
            self.obj['PeptDis'].info['Name'] = 'Pairwise peptide distances'
            ### ~ [2] Optional loading of AA Distance Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('AADis').lower() not in ['','none']: self.obj['AADis'].loadMatrix(self.getStr('AADis'))
            else:
                self.obj['AAProp'] = aaprop = rje_aaprop.AAPropMatrix(self.log,self.cmd_list)
                #aaprop.readAAProp()    # Does this on loading!
                for aa in aaprop.pdif: self.obj['AADis'].addDis(aa[0],aa[1],aaprop.pdif[aa])
            return True
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def preparePeptides(self):  ### Prepares peptides for clustering
        '''Prepares peptides for clustering.'''
        try:### ~ [1] Check Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mypep = []; peplen = 0; pepaln = True
            for inpep in self.list['Peptides']:
                pep = string.split(inpep,'>')[0]
                pep = string.split(pep,'#')[0]
                for pep in string.split(pep)[0:]:
                    pep = string.replace(pep,'^','-')
                    pep = string.replace(pep,'$','-')
                    if pep:
                        if pep in mypep: self.warnLog('Redundant peptide "%s"' % pep); continue
                        mypep.append(pep)
                        if peplen and len(pep) != peplen: pepaln = False
                        else: peplen = len(pep)
            inx = len(self.list['Peptides'])
            self.list['Peptides'] = mypep
            peptxt = {True:'fragments',False:'peptides'}[self.getBool('DNA')]
            self.printLog('#PEP','%s %s read from %s input %s' % (len(mypep),peptxt,inx,peptxt))
            if not pepaln: self.printLog('#WARN','Warning: %s are not all same length!' % peptxt)
            #!# Consider adding interactive trigger of alignment
        except: self.errorLog('Problem during %s preparePeptides.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def peptAlign(self,regex=None,peptides=[],peptdis=None,termini=None,save=False): ### Align peptides using regular expression
        '''
        Align peptides using regular expression.
        >> regex:str [None] = Regular expression to use for alignment of peptides.
        >> peptides:list [] = List of peptides to align using regex.
        >> peptdis:str [None] = Peptide distance method to use first.
        >> termini:bool [None] = Whether peptides for alignment have termini (^ & $) or X flanking regex match.
        >> save:bool [True] = Whether to save to
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['PeptDis']: self.setup()
            failx = 0       # Number of failures
            ## ~ [0a] Setup method attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not peptides: peptides = self.list['Peptides']
            if not regex:
                if self.getStrLC('PeptAlign'): regex = self.getStrUC('PeptAlign')
                else: return peptides[0:]
            if termini == None: termini = self.getBool('Termini')
            if not peptdis: peptdis = self.getStrLC('PeptDis')
            ## ~ [0b] Setup SLiM and alignment attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if regex in ['T','TRUE']:   # SLiM-free alignment
                maxlen = 0; maxgapx = 0
                slimvar = {}    # Dictionary of {gapx:gap pos}
                for pept in peptides: maxlen = max(maxlen,len(pept))
                for pept in peptides: maxgapx = max(maxgapx,maxlen-len(pept))
                maxgapx = min(maxgapx,self.getInt('MaxGapX'))
                for gapx in range(1,maxgapx+1):
                    slimvar[gapx] = []
                    peptlen = maxlen - gapx
                    if termini: gapcombos = rje.listCombos([range(peptlen)[1:-1]] * gapx,maxrep=self.getInt('MaxGapVar'))
                    else: gapcombos = rje.listCombos([range(peptlen)] * gapx,maxrep=self.getInt('MaxGapVar'))
                    for gapvar in gapcombos[0:]:
                        gapvar.sort(); gapvar.reverse()
                        if gapvar not in slimvar[gapx]: slimvar[gapx].append(gapvar)
                    self.printLog('#GAPX','PeptLen: %d; MaxLen: %d; Termini: %s => %s x %d gap variants.' % (peptlen,maxlen,termini,rje.iLen(slimvar[gapx]),gapx))
                #slimvar = range(maxlen)
                #if termini: slimvar = slimvar[1:-1]     # All possible positions for gaps. Generate combos as needed.
            else:
                #!# Need to deal with multiple regex?! (Use one with most matches and only keep that one?!)
                if rje_slim.needToSplitPattern(regex):
                    try: splits = rje_slim.splitPattern(regex)
                    except:
                        self.errorLog('SLiM "%s" splitPattern failure' % regex)
                        self.warnLog('Will try SLiM-free alignment following splitPattern() error.')
                        return self.peptAlign('TRUE',peptides,peptdis,termini,save)
                    self.printLog('#SPLIT','%s => %s' % (regex,string.join(splits,' | ')))
                    newregex = ''; bestpep = []
                    for regsplit in splits:
                        regexpep = []
                        for pept in peptides:
                            if termini and rje.matchExp('(%s)' % regsplit,pept[1:-1]): regexpep.append(pept)
                            elif not termini and rje.matchExp('(%s)' % regsplit,pept): regexpep.append(pept)
                        if len(regexpep) > len(bestpep): bestpep = regexpep[0:]; newregex = regsplit
                    self.printLog('#REGEX','%s => %s (%d/%d peptides)' % (regex,newregex,len(bestpep),len(peptides)))
                    regex = newregex
                    for pept in peptides[0:]:
                        if pept not in bestpep: self.warnLog('%s does not match %s!' % (pept,regex)); peptides.remove(pept); failx += 1
                slim = rje_slim.slimFromPattern(regex)
                self.printLog('#GUIDE','SLiM Guide: %s' % slim)
                slimpos = string.split(slim,'-')
                maxlen = rje_slim.slimLen(slim)
                if regex.startswith('^'): maxlen -= 1
                if regex.endswith('$'): maxlen -= 1
                slimvar = []    # Make variants of SLiMs (wildcard spacers only)
                w = 1
                while len(slimpos) > w: slimvar.append(slimpos[w]); w += 2  # Add wildvar spacers only
                maxvar = []
                for var in slimvar: maxvar.append(var[-1])   # Smallest number of wildcards: used to assess - to add
                slimvar = rje.listCombos(slimvar)   # Returns all possible combinations, used for building variants
                if termini: maxlen += 2
            ## ~ [0c] Setup Peptide Distance methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dismethods = ['id','prop','pam']
            if peptdis:
                try: dismethods.insert(0,dismethods.pop(dismethods.index(peptdis)))
                except: self.warnLog('PeptDis method "%s" not recognised.' % peptdis)

            ### ~ [1] Cycle through peptides and make variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            variants = {}   # Dictionary of {peptide:[variants]}
            singletons = [] # Peptides with single variants
            for pept in peptides[0:]:
                self.progLog('\r#VAR','%s peptides: %s singletons; %s with possible variants.' % (rje.iLen(peptides),rje.iLen(singletons),rje.iLen(variants)))
                variants[pept] = []
                ## ~ [1a] ~ Make list of peptide length variants adding - at all possible positions ~ ##
                if regex in ['T','TRUE']:   # SLiM-free alignment
                    gapx = maxlen - len(pept)
                    if gapx > self.getInt('MaxGapX'):
                        self.warnLog('Peptide %s exceeds MapGapX=%d; rejected.' % (pept,self.getInt('MaxGapX')))
                        peptides.remove(pept)
                        continue
                    self.bugPrint(slimvar)
                    self.debug('%s: %s vs %s = %d' % (pept,len(pept),maxlen,gapx))
                    if gapx:    # Try all gap combinations
                        for gapvar in slimvar[gapx]:    #rje.listCombos([slimvar] * gapx):
                            peptvar = pept
                            for gap in gapvar: peptvar = peptvar[:gap] + '-' + peptvar[gap:]
                            if peptvar not in variants[pept]: variants[pept].append(peptvar)
                    else: variants[pept] = [pept]
                ## ~ [1b] ~ Make list of peptide length variants, adding - to regex wildvar positions ##
                else:
                    self.debug(slimvar)
                    for var in slimvar:
                        self.debug(var)
                        peptvar = ''    # Add new variant
                        i = 0
                        if termini: peptvar += pept[i]; i += 1
                        if regex[0] != '^': peptvar += pept[i]; i += 1
                        for wi in range(len(var)):
                            wy = int(maxvar[wi]); wx = int(var[wi])
                            if wx: peptvar += pept[i:i+wx]; i += wx
                            peptvar += '-' * (wy - wx)  # Add a number of gaps equal to maxvar for same position minus slimvar
                            if i >= len(pept): break
                            peptvar += pept[i]; i += 1
                        # Keep variants that match regex and maxlen
                        if termini:
                            if regex[-1] != '$':
                                if i < len(pept): peptvar += pept[i]
                                i += 1
                            self.debug('%s >> %s' % (pept,peptvar))
                            rmatch = rje.matchExp('(%s)' % regex,peptvar[1:-1])
                            self.bugPrint('%s vs %s: %s' % (regex,peptvar[1:-1],rmatch))
                            self.debug('%s vs %s and %s vs %s' % (i,len(pept),len(peptvar),maxlen))
                            keepvar = i == len(pept) and rmatch and len(peptvar) == maxlen
                        else: keepvar = i == len(pept) and rje.matchExp('(%s)' % regex,peptvar) and len(peptvar) == maxlen
                        self.bugPrint('%s %s: %s x %s => %s = %s' % (regex,slim,var,pept,peptvar,keepvar))
                        if keepvar and peptvar not in variants[pept]: variants[pept].append(peptvar)
                ## ~ [1c] ~ Check Peptide variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if len(variants[pept]) == 1: singletons.append(variants.pop(pept)[0])
                elif not variants[pept]: self.warnLog('No %s gap variants completely match %s!' % (pept,regex)); variants.pop(pept); failx += 1
            self.printLog('#VAR','%s peptides: %s singletons; %s with possible variants.' % (rje.iLen(peptides),rje.iLen(singletons),rje.iLen(variants)))
            self.debug(singletons)
            self.debug(variants)

            ### ~ [2] ~ Sort peptides by increasing numbers of variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # For 2+ variants, rank by mean PeptDist versus single variants
            # Keep best (including ties) and cycle
            # Iterate until no more variants filtered
            # If variants remain, switch score method and iterate again
            # If variants remain after all score methods, keep first variant
            peptdis = dismethods.pop(0)
            while variants:
                self.progLog('\r#VAR','%s peptide tidy: %s singletons; %s with variants.' % (rje.iLen(peptides),rje.iLen(singletons),rje.iLen(variants)))
                if self.obj['PeptDis']: self.obj['PeptDis'].dict['Matrix'] = {}
                comppept = singletons[0:]
                prevarx = 0; postvarx = 0
                if not comppept: comppept = rje.listJoin(variants.values(),sortunique=True)
                for pept in rje.sortKeys(variants):
                    self.progLog('\r#VAR','%s peptide tidy: %s singletons; %s with variants.' % (rje.iLen(peptides),rje.iLen(singletons),rje.iLen(variants)))
                    scores = {}; prevarx += len(variants[pept])
                    for peptvar in variants[pept]:
                        dis = 0.0
                        for pep2 in comppept:
                            if termini: dis += self.peptDist(peptvar[1:-1],pep2[1:-1],peptdis)
                            else: dis += self.peptDist(peptvar,pep2,peptdis)
                        dis /= len(comppept)
                        if dis not in scores: scores[dis] = []
                        scores[dis].append(peptvar)
                    variants[pept] = scores.pop(rje.sortKeys(scores)[0])    # Keep lowest scoring variant(s)
                    if len(variants[pept]) == 1: singletons.append(variants.pop(pept)[0])
                    else: postvarx += len(variants[pept])
                self.printLog('#PDIS','%s distances: %s => %s variants.' % (peptdis,prevarx,postvarx))
                if prevarx == postvarx:
                    if dismethods: peptdis = dismethods.pop(0)
                    else: break
            self.printLog('#VAR','%s peptides tidied: %s singletons; %s with variants.' % (rje.iLen(peptides),rje.iLen(singletons),rje.iLen(variants)))
            if variants:
                self.warnLog('Unable to select all variants using distances.')
                self.printLog('#VAR','Arbitrary variants picked for %s peptides' % rje.iLen(variants))
                for pept in rje.sortKeys(variants): singletons.append(variants.pop(pept)[0])

            ### ~ [3] Remove 100% gapped positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(maxlen-1,-1,-1):
                degap = True
                for pept in singletons:
                    if pept[i] != '-': degap = False; break
                if degap:
                    for p in range(len(singletons)): singletons[p] = singletons[p][:i] + singletons[p][i+1:]

            ### ~ [4] Save and/or return peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save:
                open('%s.aligned.txt' % self.baseFile(),'w').write(string.join(singletons,'\n'))
                self.printLog('#OUT','%s aligned peptides output to %s.aligned.txt' % (rje.iLen(singletons),self.baseFile()))
            return singletons

        except: self.errorLog('%s.peptAlign error' % self); raise
#########################################################################################################################
    def peptDist(self,pep1,pep2,dismethod=None,force=False):### Calculate peptide distance and updates self.obj['PeptDis']
        '''
        Calculate peptide distance and updates self.obj['PeptDis'].
        >> pep1:str = Peptide sequence 1
        >> pep2:str = Peptide sequence 2
        >> dismethod:str [None] = Distance method to use. If None, will use self.str['PeptDis']
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if pep1 == pep2: self.obj['PeptDis'].addDis(pep1,pep2,0.0)
            elif not force and self.obj['PeptDis'].getDis(pep1,pep2,None): return self.obj['PeptDis'].getDis(pep1,pep2,None)
            if not dismethod: dismethod = self.getStr('PeptDis')
            ### ~ [2] Calculate Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Simple Identity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if dismethod == 'id':
                dis = 0.0; res = 0.0
                for i in range(len(pep1)):
                    if i >= len(pep2): dis += 1; res += 1; continue
                    if pep1[i] != pep2[i] or pep1[i] == 'X': dis += 1
                    if pep1[i] != '-' or pep2[i] != '-': res += 1
                dis = dis/res
            ## ~ [2b] AA Properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif dismethod == 'prop':
                dis = 0.0; res = 0.0
                for i in range(len(pep1)):
                    if i >= len(pep2): continue
                    if pep1[i] == pep2[i] == '-': continue   # Replace with positive check for alphabet?
                    if self.obj['AADis'].getDis(pep1[i],pep2[i],None) == None:
                        self.warnLog('Cannot find AA Distance for %s vs %s. (Will use %s)' % (pep1[i],pep2[i],self.obj['AADis'].maxDis()))
                    dis += self.obj['AADis'].getDis(pep1[i],pep2[i],self.obj['AADis'].maxDis())
                #self.obj['PeptDis'].addDis(pep1,pep2,dis)
            ## ~ [2a] PAM Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif dismethod == 'pam':
                if not self.obj['PAM']: self.obj['PAM'] = rje_pam.PamCtrl(self.log,self.cmd_list)
                dis = self.obj['PAM'].pamML(desc='pamML',ancseq=pep1,descseq=pep2)
                #self.errorLog('Peptide Distance Method "%s" not ready!' % dismethod,printerror=False)
                #raise ValueError
            ## ~ [2x] Unrecognised method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                self.errorLog('Peptide Distance Method "%s" not recognised' % dismethod,printerror=False)
                raise ValueError
            ### ~ [3] Finish and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['PeptDis'].addDis(pep1,pep2,dis)
            return dis
        except ValueError: raise 
        except: self.errorLog('%s.peptDis error' % self); 
#########################################################################################################################
    def peptCluster(self):  ### Performs actual peptide clustering and stores results in self.obj['Tree']
        '''Performs actual peptide clustering and stores results in self.obj['Tree'].'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            base = rje.baseFile(self.getStr('SaveDis'))
            pretree = ['treeformats=nwk,text','basefile=%s' % base]
            ### ~ [1] ~ Phylip Neighbor method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('PeptCluster') == 'neighbor':
                disfile = '%s.phy' % base
                fasfile = '%s.fas' % base
                treecmd = ['autoload=T','maketree=neighbor','disin=%s' % disfile,'seqin=%s' % fasfile]
                pretree += ['root=mid']
                if disfile != self.getStr('SaveDis'):
                    rje.backup(self,disfile)
                    self.obj['PeptDis'].saveMatrix(filename=disfile,format='phylip')   ### Saves matrix
                if 'peptides=%s' % fasfile not in self.cmd_list:
                    rje.backup(self,fasfile)
                    FAS = open(fasfile,'w')
                    for pep in self.list['Peptides']: FAS.write('>%s\n%s\n' % (pep,pep))
                    FAS.close()
                tree = self.obj['Tree'] = rje_tree.Tree(self.log,pretree+self.cmd_list+treecmd)
            ### ~ [2] ~ UPGMA method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                if self.getStr('PeptCluster') not in ['wpgma','upgma']:
                    self.errorLog('PeptCluster method "%s" not recognised. Will use UPGMA' % self.getStr('PeptCluster'),printerror=False)
                    base = string.replace(base,self.getStr('PeptCluster'),'upgma')
                    pretree += ['basefile=%s' % base]
                if self.getStr('PeptCluster') == 'upgma': nsftree = self.obj['PeptDis'].upgma()
                elif self.getStr('PeptCluster') == 'wpgma': nsftree = self.obj['PeptDis'].wpgma()
                #nwkfile = '%s.nwk' % base
                #treecmd += ['nsfin=%s' % nwkfile]
                #rje.backup(self,nwkfile)
                #open(nwkfile,'w').write(nsftree)
                treecmd = ['autoload=F']
                tree = self.obj['Tree'] = rje_tree.Tree(self.log,pretree+self.cmd_list+treecmd)
                tree.buildTree(nsftree)
            ### ~ [3] ~ Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for node in tree.node:
                if node.info['Name'] in self.list['Peptides']: node.stat['ID'] = self.list['Peptides'].index(node.info['Name']) + 1
            tree.saveTrees()
            for outfmt in tree.list['TreeFormats']:
                treefile = '%s.%s' % (tree.info['Basefile'],rje_tree.formatext[outfmt])
                self.dict['Output'][outfmt] = treefile
        except: self.errorLog('%s.peptDis error' % self);
#########################################################################################################################
    def cgiTuple(self):  ### Returns a tuple of distance matrix,Text tree and Newick Tree.
        '''Returns a tuple of distance matrix,Text tree and Newick Tree.'''
        try:### ~ [0] ~ Return tuple ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('AADis').lower() not in ['','none']: aaprop = {}
            else: aaprop = self.obj['AAProp'].prop    #[property][aa]
            return (self.obj['PeptDis'].dict['Matrix'],self.obj['Tree'].textTree(),self.obj['Tree']._makeNSFTree(),aaprop,self.obj['AADis'].dict['Matrix'])
        except: self.errorLog('%s.cgiTuple error' % self); 
#########################################################################################################################
    ### <3> ### REST Output Methods                                                                                     #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        peptides = original input peptides
        aligned = aligned peptides
        matrix = peptide distance matrix
        X = peptide tree(s), where X is nsf/nwk/text/png/svg (`treeformats=LIST`: only `nwk` & `text` made by default.)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Iterate') and 'iterate' not in self.dict['Output']: self.dict['Output']['iterate'] = ''
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['peptides','aligned','matrix'] + self.list['TreeFormats']
#########################################################################################################################
### End of SECTION II: PeptCluster Class                                                                                #
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
    try: PeptCluster(mainlog,cmd_list).run()
        #print rje_zen.Zen().wisdom(), '\n\n *** No standalone functionality! *** \n\n'

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
