#!/usr/bin/python

# See below for name and description
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_disorder
Description:  Disorder Prediction Module
Version:      1.5.1
Last Edit:    28/05/19
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module currently has limited function and no standalone capability, though this may be added with time. It is
    designed for use with other modules. The disorder Class can be given a sequence and will run the appropriate
    disorder prediction software and store disorder prediction results for use in other programs. The sequence will have
    any gaps removed.

    Currently six disorder prediction methods are implemented:

    * `IUPred2` : From V1.5.0, IUPred2 API calls can be made to https://iupred2a.elte.hu/iupred2a/. The IUPred method
    needs to be given as `iumethod=X`, where `X` is `short`, `long`, `anchor`, or `redox`. This can also be triggered by
    setting `disorder=X` to `iushort2` or `iulong2`, `iuanchor2` (or simply `anchor2`), or `iuredox2`.

    * `IUPred`: Dosztanyi Z, Csizmok V, Tompa P & Simon I (2005). J. Mol. Biol. 347, 827-839. This has to be installed
    locally. It is available on request from the IUPred website and any use of results should cite the method. (See
    http://iupred.enzim.hu/index.html for more details.) IUPred returns a value for each residue, which by default,
    is determined to be disordered if > 0.5. The IUPred method needs to be given as `iumethod=short` or `iumethod=long`.
    This can also be triggered by setting `disorder=iushort` or `disorder=iulong`.

    * `FoldIndex`: This is run directly from the website (http://bioportal.weizmann.ac.il/fldbin/findex) and more simply
    returns a list of disordered regions. You must have a live web connection and `curl` on your system to use this method!

    * `ANCHOR`: Meszaros B, Simon I & Dosztanyi Z (2009). PLoS Comput Biol 5(5): e1000376. This has to be installed
    locally. It is available on request from the ANCHOR website and any use of results should cite the method. (See
    http://anchor.enzim.hu/ for more details.) ANCHOR returns a probability value for each residue, which by default,
    is determined to be disordered if > 0.5.

    * `Parse`: Parsed disorder from protein sequence name, e.g. DisProt download.
    #X-Y = disordered region; &X-Y = ordered region [0.0]

    * `IUScoreDir`: From V1.2.0, pre-calculated disorder scores can be loaded from a file (see below). To activate this mode,
    the `disorder=X` setting should match the `<DISORDER>` part of the filename (below). If this is one of the other
    recognised disorder predictors above, missing files will be generated if possible. Otherwise, files must be present
    for disorder prediction to occur.

    For IUPred, the individual residue results are stored in Disorder.list['ResidueDisorder']. For all methods, the
    disordered regions are stored in Disorder.list['RegionDisorder'] as (start,stop) tuples.

    ### IUScoreDir:

    V1.2.0 introduced the optional use of `IUScoreDir/` to save and/or (re)load lists of disorder scores. These are files
    named <ACC>.<DISORDER>.txt where <ACC> is the accession number of the protein. If `md5acc=T` then an md5 hash of the
    sequence is used instead: `hashlib.md5(<SEQUENCE>).hexdigest()`.
    
Commandline:
    ### General Options: ###
    disorder=X      : Disorder method to use (iupred2/iupred/foldindex/anchor/anchor2/parse) [iushort2]
    strict=T/F      : Whether to exit with error if disorder method not found [False]
    iucut=X         : Cut-off for score-based method (e.g. IUPred/ANCHOR) results [0.2]
    iumethod=X      : IUPred method to use (long/short/redox/anchor) [short]
    sequence=X      : Sequence to predict disorder for (autorun) []
    name=X          : Name of sequence to predict disorder for []
    minregion=INT   : Minimum length of an ordered/disordered region [0]
    minorder=INT    : Minimum length of an ordered region; over-rides minregion if >-1 [-1]
    smoothing=X     : Smoothing mode for minregion=X and minorder=X smoothing (foldfirst/sequence) [foldfirst]
    iuscoredir=PATH : Path in which to save protein acc.DISORDER.txt score files for re-use []
    discalculate=T/F: Whether to try to calculate disorder if existing score not loaded [True]
    md5acc=T/F      : Whether to use md5sum hexdigest hashing of sequence in place of accession numbers [True]

    ### System Settings: ###
    iupath=PATH     : The full path to the IUPred executable [c:/bioware/iupred/iupred.exe]
    anchor=PATH     : Full path to ANCHOR executable []
    filoop=INT      : Number of times to try connecting to FoldIndex server [10]
    fisleep=INT     : Number of seconds to sleep between attempts [2]
    iuchdir=T/F     : Whether to change to IUPred directory and run (True) or rely on IUPred_PATH env variable [False]

Uses general modules: copy, os, string, sys, time, urllib2
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import hashlib, os, random, string, sys, time, urllib2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added parsing of disorder from name as an option instead of disorder prediction
    # 0.2 - Added Folded tuple as well as disordered
    # 0.3 - Added PrintLog opt attribute
    # 0.4 - Added option for correct use of IUPred_PATH environment variable
    # 0.5 - Added Minimum length of an ordered/disordered region
    # 0.6 - Added ANCHOR prediction.
    # 0.7 - Added globProportion calculation.
    # 0.8 - Added makeRegions() method.
    # 1.0.0 - Added random disorder function and elevated to v1.x as fully functional for SLiMSuite
    # 1.1.0 - Added strict option for disorder method selection. Added minorder=X.
    # 1.2.0 - Added saving and loading scores to IUScoreDir/.
    # 1.3.0 - Switched default behaviour to be md5acc=T.
    # 1.4.0 - Fixed up disorder=parse and disorder=foldindex.
    # 1.5.0 - Added iupred2 and anchor2 parsing from URL using accnum. Made default disorder=iushort2.
    # 1.5.1 - Fixed iupred2 URL generation error.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add region buffer - min length of region allowed before ignoring
    # [Y] : Add a PrintLog option that controls whether printing to Log or not.
    # [ ] : Neaten and tidy
    # [ ] : Add domain-based disorder stuff?
    # [Y] : Add option to store/re-read prediction data for later use.
    # [Y] : Add reading data from IUScore directory.
    # [ ] : Add IUPred2A and ANCHOR2 URL calls, if iupred=url or anchor=url - will send accession number.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_DISORDER', '1.5.1', 'May 2019', '2008')
    description = 'Disorder Prediction Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### SECTION II: Disorder Class                                                                                          #
#########################################################################################################################
class Disorder(rje.RJE_Object):     
    '''
    Disorder Prediction Class. Author: Rich Edwards (2005).

    Info:str
    - ANCHOR = Full path to ANCHOR executable
    - Sequence = sequence to give to disorder prediction (will have gaps removed)
    - Disorder = Disorder method to use (iupred/foldindex) [None]
    - IUPath = The full path to the IUPred exectuable [c:/bioware/iupred/iupred.exe]
    - IUMethod = IUPred method to use (long/short) [short]
    - IUScoreDir=PATH : Path in which to save protein acc.DISORDER.txt score files for re-use [None]
    - Smoothing = Smoothing mode for minregion=X and minorder=X smoothing (foldfirst/sequence) [foldfirst]
    
    Opt:boolean
    - DisCalculate=T/F: Whether to try to calculate disorder if existing score not loaded [True]
    - Flat = whether ResidueDisorder is a "flat" 1/0 or graded (e.g. raw IUPred)
    - IUChDir = Whether to change to IUPred directory and run (True) or rely on IUPred_PATH env variable [False]
    - MD5Acc=T/F      : Whether to use md5sum hexdigest hashing of sequence in place of accession numbers [True]
    - PrintLog = whether to print disorder prediction status to Log [False]
    - Strict=T/F  : Whether to exit with error if disorder method not found [False]

    Stat:numeric
    - IUCut = Cut-off for IUPred results [0.2]
    - FILoop = Number of times to try connecting to FoldIndex server [10]
    - FISleep = Number of seconds to sleep between attempts [2]
    - MinOrder=X  : Minimum length of an ordered region; over-rides minregion if >-1 [-1]
    - MinRegion = Minimum length of an ordered/disordered region [0]
    
    List:list
    - ResidueDisorder = individual (IUPRed) residue results 
    - RegionDisorder = disordered regions as (start,stop) tuples (1->L)
    - RegionFold = folded (i.e. not disordered regions) as (start,stop) tuples (1->L)
    
    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Sequence','Disorder','IUPath','IUMethod','ANCHOR','IUScoreDir']
        self.optlist = ['DisCalculate','Flat','PrintLog','IUChDir','Strict','MD5Acc']
        self.statlist = ['IUCut','FILoop','FISleep','MinOrder','MinRegion','Smoothing']
        self.listlist = ['ResidueDisorder','RegionDisorder','RegionFold']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='',opt=False,stat=0.2,obj=None,setlist=True,setdict=True)
        self.setInfo({'IUPath':rje.makePath('c:/bioware/iupred/iupred.exe',wholepath=True),'IUMethod':'short',
                      'Disorder':'iushort2','Smoothing':'foldfirst','IUScoreDir':''})
        self.setStat({'FILoop':10,'FISleep':2,'MinOrder':-1,'MinRegion':0,'IUCut':0.2})
        self.setOpt({'MD5Acc':True,'DisCalculate':True})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'info',['Disorder','Name','Sequence','IUMethod','Smoothing'])
                self._cmdRead(cmd,type='stat',att='IUCut')
                self._cmdRead(cmd,type='fullpath',att='IUPath')
                self._cmdRead(cmd,type='fullpath',att='ANCHOR')
                self._cmdReadList(cmd,'int',['FILoop','FISleep','MinOrder','MinRegion'])
                self._cmdReadList(cmd,'opt',['IUChDir','Strict','MD5Acc'])
                self._cmdReadList(cmd,'path',['IUScoreDir'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        ## ~ Setup Disorder method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        disorder = self.getStrLC('Disorder')
        if disorder in ['iupred','iupred2']:
            self.info['Disorder'] = 'iu%s' % self.getStrLC('IUMethod')
            if disorder.endswith('2'): self.info['Disorder'] = '%s2' % self.getStrLC('Disorder')
        ### AutoRun ###
        if self.info['Sequence'].lower() not in ['none','']: self.disorder(self.info['Sequence'])
#########################################################################################################################
    def md5hash(self,sequence=''): return hashlib.md5(sequence).hexdigest()
#########################################################################################################################
    ### <3> ### Disorder Prediction Methods                                                                             #
#########################################################################################################################
    def disorder(self,sequence='',name=''):     ### Takes a sequence, degaps, and runs disorder prediction
        '''
        Takes a sequence, degaps, and runs disorder prediction.
        >> sequence:str = protein sequence
        >> name:str = (optional) name for sequence - goes in self.info['Name']
        '''
        try:### ~ [1] ~ Setup sequence and name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sequence: self.info['Sequence'] = string.join(string.split(sequence,'-'),'')
            if name: self.info['Name'] = name
            sname = string.split(self.info['Name'])[0]
            if not self.info['Sequence']:
                self.log.errorLog('Cannot calculate disorder: no AAs given in sequence!',printerror=False)
                return False
            ### ~ [2] ~ Load Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            loaded = self.loadDisorder()
            #self.debug(loaded)
            if loaded: return loaded
            if not self.getBool('DisCalculate'):
                self.warnLog('Failed to load %s %s disorder (discalculate=F)' % (sname,self.getStrLC('Disorder')))
                return self.noDisorder()
            ### ~ [3] ~ Run Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disorder = self.getStrLC('Disorder')
            if disorder in ['anchor2','iuanchor2','iushort','iulong','iushort2','iulong2','iuredox2']: return self.iuPred()
            elif disorder == 'foldindex': return self.foldIndex()
            elif disorder == 'anchor': return self.ANCHOR()
            elif disorder == 'parse': return self.parseDisorder()
            elif disorder == 'random': return self.randomDisorder()
            else:
                if self.getBool('Strict'):
                    self.errorLog('Cannot calculate %s disorder: disorder method "%s" not recognised (strict=T)!' % (sname,disorder),printerror=False)
                    return False
                else:
                    self.warnLog('No disorder calculation for %s: disorder method "%s" not recognised (strict=F)!' % (sname,disorder))
                    return self.noDisorder()
        except:
            self.log.errorLog('Error in Disorder.disorder(%s)' % name,quitchoice=True)
            return False
#########################################################################################################################
    def loadDisorder(self): ### Looks for existing disorder prediction results and loads if found.
        '''Looks for existing disorder prediction results and loads if found.'''
        try:### ~ [1] ~ Setup sequence and name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sname = string.split(self.getStr('Name'))[0]
            if not self.getStrLC('IUScoreDir'): return False
            sequence = self.getStrUC('Sequence')
            if self.getBool('MD5Acc'): acc = self.md5hash(sequence)
            else: acc = string.split(sname,'__',maxsplit=1)[-1]
            ### ~ [2] ~ Look for file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ifile = '%s%s.%s.txt' % (self.getStr('IUScoreDir'),acc,self.getStrLC('Disorder'))
            #self.debug(ifile)
            if not os.path.exists(ifile): return False
            ### ~ [3] ~ Parse file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['ResidueDisorder'] = []
            fline = open(ifile,'r').readline()
            for dstr in string.split(fline)[1:]: self.list['ResidueDisorder'].append(string.atof(dstr))
            if len(self.list['ResidueDisorder']) != len(sequence):
                self.errorLog('%s Disorder score length mismatch (%d score vs %d pos)' % (sname,len(self.list['ResidueDisorder']),len(sequence)),printerror=False)
                self.list['ResidueDisorder'] = []; return False
            return True
        except:
            self.log.errorLog('Error in Disorder.loadDisorder(%s)' % sname,quitchoice=True)
            return False
#########################################################################################################################
    def saveDisorder(self): ### Saves disorder prediction results to file.
        '''Saves disorder prediction results to file.'''
        try:### ~ [1] ~ Setup sequence and name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sname = string.split(self.getStr('Name'))[0]
            if not self.getStrLC('IUScoreDir'): return False
            sequence = self.getStrUC('Sequence')
            if self.getBool('MD5Acc'):
                acc = self.md5hash(sequence)
                sname = sequence
            else: acc = string.split(sname,'__',maxsplit=1)[-1]
            ### ~ [2] ~ Set up file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ifile = '%s%s.%s.txt' % (self.getStr('IUScoreDir'),acc,self.getStrLC('Disorder'))
            ### ~ [3] ~ Save to File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['ResidueDisorder']:  # Should have scores
                dlist = []
                for x in self.list['ResidueDisorder']:
                    dlist.append('%f' % x)
                rje.mkDir(self,ifile)
                open(ifile,'w').write('%s\t%s\n' % (sname,string.join(dlist)))
        except:
            self.log.errorLog('Error in Disorder.saveDisorder(%s)' % sname,quitchoice=True)
            return False
#########################################################################################################################
    def iuPred(self,retry=2):     ### Runs IUPred disorder prediction
        '''Runs IUPred disorder prediction.'''
        mydir = os.path.abspath(os.curdir)
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            disorder = self.getStrLC('Disorder')
            iumethod = self.getStrLC('IUMethod')
            if disorder.startswith('iushort'): iumethod = 'short'
            elif disorder.startswith('iulong'): iumethod = 'long'
            elif disorder.startswith('iuredox'): iumethod = 'redox'
            elif 'anchor' in disorder: iumethod = 'anchor'
            ## ~ [1b] Setup sequence and temp file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sequence = self.info['Sequence'].upper()
            name = self.info['Name'][:4] + rje.randomString(8)
            tmp = name + '.tmp'
            sname = string.split(self.getStr('Name'))[0]
            acc = string.split(sname,'__',maxsplit=1)[-1]

            #!# Temp shunt to old method #!#
            if retry < 2 and disorder.endswith('2'): disorder = disorder[:-1]

            ### ~ [2] Run Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dlines = []
            ## ~ [2a] Run IUPred2 server online using accnum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if disorder.endswith('2'):
                url = 'https://iupred2a.elte.hu/iupred2a/%s/%s' % (iumethod,acc)
                try:
                    dlines = urllib2.urlopen(url).readlines()
                    if 'not found' in dlines[0]: raise ValueError(dlines[0])
                except:
                    self.errorLog(url)
                    self.warnLog('%s disorder failure for %s: trying %s' % (disorder,self.info['Name'],disorder[:-1]))
                    disorder = disorder[:-1]

            ## ~ [2b] Run IUPred2 on the commandline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Add an option to try online then switch to local?
            if not disorder.endswith('2'):
                iupath = string.join(string.split(self.info['IUPath'],os.sep)[:-1],os.sep)
                iupred = string.split(self.info['IUPath'],os.sep)[-1]
                if self.opt['IUChDir']: os.chdir(string.join(string.split(self.info['IUPath'],os.sep)[:-1],os.sep))
                open(tmp,'w').write('>%s\n%s\n' % (name,sequence))
                if self.opt['IUChDir'] and self.opt['Win32']: iucmd = '%s %s %s' % (iupred,tmp,iumethod)
                elif self.opt['IUChDir']: iucmd = './%s %s %s' % (iupred,tmp,iumethod)
                else: iucmd = '%s %s %s' % (self.info['IUPath'],tmp,iumethod)
                dlines = os.popen(iucmd).readlines()
                try: os.unlink(tmp)
                except: self.errorLog('Cannot delete %s!' % tmp)
                if self.opt['IUChDir']: os.chdir(mydir)
            if self.info['Name'] not in ['','None']: name = self.info['Name']
            self.list['ResidueDisorder'] = []

            ### ~ [3] Parse Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            iuregex = '^\s*(\d+)\s+(\S)\s+(\S+)'
            if iumethod in ['anchor','redox']: iuregex = '^\s*(\d+)\s+(\S)\s+\S+\s+(\S+)'
            for d in dlines:
                if rje.matchExp(iuregex,d):
                    dm = rje.matchExp(iuregex,d)
                    pos = string.atoi(dm[0])
                    aa = dm[1]
                    score = string.atof(dm[2])
                    i = len(self.list['ResidueDisorder'])
                    if sequence[i] != aa:
                        self.log.errorLog('%s: Position %d is %s in sequence but %s in IUPred output!' % (name,pos,sequence[i],aa),printerror=False)
                        raise ValueError
                    if pos != (i + 1):
                        self.log.errorLog('%s: Position %d reached in IUPred output but previous results missing!' % (name,pos),printerror=False)
                        raise ValueError
                    self.list['ResidueDisorder'].append(score)
            if len(self.list['ResidueDisorder']) != len(sequence):
                self.log.errorLog('%s: Sequence = %d aa but IUPred results stop at %s!' % (name,len(sequence),len(self.list['ResidueDisorder'])),printerror=False)
                raise ValueError
            self.saveDisorder()

            ### ~ [4] Make Regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            start = 0
            fstart = 0
            i = 0
            dx = 0
            while i < len(sequence):
                score = self.list['ResidueDisorder'][i]
                i += 1
                if not start and score > self.stat['IUCut']:    ### Start new disorder ###
                    start = i
                elif start and score <= self.stat['IUCut']:     ### End!
                    self.list['RegionDisorder'].append((start,i-1))
                    dx += i - start
                    start = 0
                if not fstart and score <= self.stat['IUCut']:    ### Start new fold ###
                    fstart = i
                elif fstart and score > self.stat['IUCut']:     ### End!
                    self.list['RegionFold'].append((fstart,i-1))
                    fstart = 0
            if start:
                self.list['RegionDisorder'].append((start,len(sequence)))
                dx += len(sequence) + 1 - start
            if fstart: self.list['RegionFold'].append((fstart,len(sequence)))
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','IUPred (%s) Disorder prediction complete: %d disorder regions, %d disordered aa' % (self.info['IUMethod'].lower(),len(self.list['RegionDisorder']),dx))
            return True
        except:
            if self.opt['IUChDir']: os.chdir(mydir)            
            if retry:
                self.printLog('#RETRY','Trying %s again...' % name)
                return self.iuPred(retry-1)
            self.log.errorLog('Error in Disorder.iuPred(%s). Disorder prediction failed. Check (setenv?) IUPred_PATH environment variable.' % name)
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            #try: os.system('rm %s*tmp' % (rje.makePath(os.path.split(self.info['IUPath'])[0])))
            #except: pass
            return False
#########################################################################################################################
    def makeRegions(self,sequence=None): ### Make regions from ResidueDisorder List                                 #V0.8
        '''Make regions from ResidueDisorder List.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not sequence: sequence = self.info['Sequence'].upper()
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            start = 0
            fstart = 0
            i = 0
            dx = 0
            while i < len(sequence):
                if i >= len(self.list['ResidueDisorder']):
                    raise ValueError('Sequence length (%d aa) exceeds %d disorder scores!' % (len(sequence),len(self.list['ResidueDisorder'])))
                else: score = self.list['ResidueDisorder'][i]
                i += 1
                if not start and score > self.stat['IUCut']:    ### Start new disorder ###
                    start = i
                elif start and score <= self.stat['IUCut']:     ### End!
                    self.list['RegionDisorder'].append((start,i-1))
                    dx += i - start
                    start = 0
                if not fstart and score <= self.stat['IUCut']:    ### Start new fold ###
                    fstart = i
                elif fstart and score > self.stat['IUCut']:     ### End!
                    self.list['RegionFold'].append((fstart,i-1))
                    fstart = 0
            if start:
                self.list['RegionDisorder'].append((start,len(sequence)))
                dx += len(sequence) + 1 - start
            if fstart: self.list['RegionFold'].append((fstart,len(sequence)))
            self.minRegion()
            if self.opt['PrintLog']: self.printLog('\r#DIS','%s Disorder prediction complete: %d disorder regions, %d disordered aa' % (self.getStr('Disorder'),len(self.list['RegionDisorder']),dx))
            return True
        except: self.errorLog('Disorder.makeRegions(%s) error' % self.info['Name']); raise
#########################################################################################################################
    def ANCHOR(self,retry=2):     ### Runs ANCHOR disorder prediction
        '''Runs ANCHOR disorder prediction.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Setup sequence and temp file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sequence = self.info['Sequence'].upper()
            name = self.info['Name'][:4] + rje.randomString(8)
            tmp = name + '.tmp'
            ## ~ [0b] ~ Setup ANCHOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            apath = self.info['ANCHOR']
            if os.path.basename(apath) == 'anchor': apath = os.path.dirname(apath)
            anchor = rje.makePath(apath) + 'anchor'
            if not os.path.exists(anchor):
                self.errorLog('Path "%s" not found!' % anchor,printerror=False)
                retry = 0; raise IOError
            ### ~ [1] Run ANCHOR Disorder prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            open(tmp,'w').write('>%s\n%s\n' % (name,sequence))
            acmd = '%s %s -d %s' % (anchor,tmp,apath)
            dlines = os.popen(acmd).readlines()
            try: os.unlink(tmp)
            except: self.errorLog('Cannot delete %s!' % tmp)
            ### ~ [2] Read in results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'] not in ['','None']: name = self.info['Name']
            self.list['ResidueDisorder'] = []
            for d in dlines:
                if d[:1] == '#': continue
                if rje.matchExp('^(\d+)\s+(\S)\s+(\S+)',d):
                    dm = rje.matchExp('^(\d+)\s+(\S)\s+(\S+)',d)
                    pos = string.atoi(dm[0])
                    aa = dm[1]
                    score = string.atof(dm[2])
                    i = len(self.list['ResidueDisorder'])
                    if sequence[i] != aa:
                        self.log.errorLog('%s: Position %d is %s in sequence but %s in ANCHOR output!' % (name,pos,sequence[i],aa),printerror=False)
                        raise ValueError
                    if pos != (i + 1):
                        self.log.errorLog('%s: Position %d reached in ANCHOR output but previous results missing!' % (name,pos),printerror=False)
                        raise ValueError
                    self.list['ResidueDisorder'].append(score)
            if len(self.list['ResidueDisorder']) != len(sequence):
                self.log.errorLog('%s: Sequence = %d aa but ANCHOR results stop at %s!' % (name,len(sequence),len(self.list['ResidueDisorder'])),printerror=False)
                raise ValueError
            self.saveDisorder()
            ### ~ [3] ~ Make Regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            start = 0
            fstart = 0
            i = 0
            dx = 0
            while i < len(sequence):
                score = self.list['ResidueDisorder'][i]
                i += 1
                if not start and score > self.stat['IUCut']:    ### Start new disorder ###
                    start = i
                elif start and score <= self.stat['IUCut']:     ### End!
                    self.list['RegionDisorder'].append((start,i-1))
                    dx += i - start
                    start = 0
                if not fstart and score <= self.stat['IUCut']:    ### Start new fold ###
                    fstart = i
                elif fstart and score > self.stat['IUCut']:     ### End!
                    self.list['RegionFold'].append((fstart,i-1))
                    fstart = 0
            if start:
                self.list['RegionDisorder'].append((start,len(sequence)))
                dx += len(sequence) + 1 - start
            if fstart: self.list['RegionFold'].append((fstart,len(sequence)))
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','ANCHOR Disorder prediction complete: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),dx))
            return True
        except:
            if retry:
                self.printLog('#RETRY','Trying %s again...' % name)
                return self.ANCHOR(retry-1)
            self.log.errorLog('Error in Disorder.ANCHOR(%s). Disorder prediction failed.' % name)
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            return False
#########################################################################################################################
    def summary(self):  ### Returns a string summary of the disorder.
        '''Returns a string summary of the disorder.'''
        #!# Add this method
        if self.opt['PrintLog']: self.log.printLog('\r#DIS','IUPred (%s) Disorder prediction complete: %d disorder regions, %d disordered aa' % (self.info['IUMethod'].lower(),len(self.list['RegionDisorder']),dx))
        return
#########################################################################################################################
    def foldIndex(self):     ### Runs FoldIndex disorder prediction
        '''Runs FoldIndex disorder prediction.'''
        try:
            ### Setup sequence and name ###
            sequence = self.info['Sequence']
            name = self.info['Name'][:4] + rje.randomString(8)
            tmp = name + '.tmp'

            ### Run Disorder ###
            retry = self.stat['FILoop']
            #X#url = "http://bioportal.weizmann.ac.il/fldbin/findex"
            url = "https://fold.weizmann.ac.il/fldbin/findex"
            params = "m=xml&sq=" + sequence # + "  "
            #url = '%s?%s' % (url,params)
            while retry:
                # try:
                #     self.bugLog('#URL','%s?%s' % (url,params))
                #     flines = urllib2.urlopen(url, params).readlines()
                #     #flines = urllib2.urlopen(url).readlines()
                # except:
                #     #!# Catch urllib2.HTTPError: HTTP Error 404: Not Found
                #     flines = []

                flines = []
                syscmd = 'curl -o %s "%s?%s"' % (tmp,url,params)
                self.bugLog('#SYS','%s?%s' % (url,params))
                os.system(syscmd)
                if rje.exists(tmp):
                    flines = open(tmp,'r').readlines()
                    os.unlink(tmp)
                if flines:
                    break
                retry -= 1
                time.sleep(self.stat['FISleep'])
            if not flines:
                self.log.errorLog('FoldIndex run for "%s" failed.' % self.info['Name'],printerror=False)
                self.list['ResidueDisorder'] = []
                self.list['RegionDisorder'] = []
                return False
            ### Process ###
            self.list['ResidueDisorder'] = [0.0] * len(sequence)
            self.list['RegionDisorder'] = []
            for f in flines:
                if rje.matchExp('<segment start="(\d+)" end="(\d+)" len="(\d+)"',f):
                    fm = rje.matchExp('<segment start="(\d+)" end="(\d+)" len="(\d+)"',f)
                    self.list['RegionDisorder'].append((string.atoi(fm[0]),string.atoi(fm[1])))
                    for i in range(string.atoi(fm[0])-1,string.atoi(fm[1])):
                        self.list['ResidueDisorder'][i] = 1.0
            ### Update Regions
            self.updateRegions()
            self.saveDisorder()
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','FoldIndex Disorder prediction complete: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),sum(self.list['ResidueDisorder'])))
            self.opt['Flat'] = True
            return True
        except:
            self.log.errorLog('Error in Disorder.foldIndex(%s)' % self.info['Name'],quitchoice=True)
            return False
#########################################################################################################################
    def parseDisorder(self):    ### Parses disordered regions from sequence name (e.g. DisProt download)
        '''
        Parses disordered regions from sequence name (e.g. DisProt download).
        #X-Y = disordered region [1.0]; &X-Y = ordered region [0.0]; All else neutral [0.5];

        If no disordered regions are found, background will be set to disordered [1.0]. Else, if no ordered regions are
        found, the background will be set to ordered [0.0]. If both are found, missing regions are neutral [0.5].
        '''
        try:
            ### Setup sequence and name ###
            sequence = self.info['Sequence']
            name = self.info['Name']
            self.list['RegionDisorder'] = []
            scoredict = {'#':1.0,'&':0.0}

            ### Set background
            disreg = False; ordreg = False
            for region in string.split(name)[1:]:
                if rje.matchExp('^[#&](\d+)-(\d+)',region):
                    if rje.matchExp('^([#&])(\d+)-(\d+)',region)[0] == '#': disreg = True
                    elif rje.matchExp('^([#&])(\d+)-(\d+)',region)[0] == '&': outreg = True
            if not disreg: bgscore = 1.0
            elif not ordreg: bgscore = 0.0
            else: bgscore = 0.5
            self.list['ResidueDisorder'] = [bgscore] * len(sequence)

            ### Process ###
            for region in string.split(name)[1:]:
                if rje.matchExp('^[#&](\d+)-(\d+)',region):
                    (i,x,y) = rje.matchExp('^([#&])(\d+)-(\d+)',region)
                    score = scoredict[i]
                    start = string.atoi(x)
                    end = string.atoi(y)
                    for r in range(start,end): self.list['ResidueDisorder'][r] = score
                    if i == '#': self.list['RegionDisorder'].append((start,end))
                    elif i == '&': self.list['RegionFold'].append((start,end))
            ### Update Regions
            self.updateRegions()

            #self.bugPrint(name)
            #self.debug(self.list['ResidueDisorder'])
            #self.debug(self.list['RegionDisorder'])
            self.saveDisorder()
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','DisProt Disorder parsing complete: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),self.list['ResidueDisorder'].count(1.0)))
            return True
        except:
            self.log.errorLog('Error in Disorder.foldIndex(%s)' % self.info['Name'],quitchoice=True)
            return False
#########################################################################################################################
    def updateRegions(self):    ### Updates missing RegionFold or RegionDisorder if only one added.
        '''
        Updates missing RegionFold or RegionDisorder if only one added.
        :return:
        '''
        sequence = self.info['Sequence']
        self.list['RegionDisorder'].sort()
        disreg = self.list['RegionDisorder']
        self.list['RegionFold'].sort()
        ordreg = self.list['RegionFold']
        if not disreg and not ordreg: self.list['RegionDisorder'].append((1,len(sequence)))
        elif not self.list['RegionFold']:
            if self.list['RegionDisorder'][0][0] > 1:
                self.list['RegionFold'].append((1,self.list['RegionDisorder'][0][0]-1))
            for n in range(len(self.list['RegionDisorder'])-1):
                if self.list['RegionDisorder'][n][1] < self.list['RegionDisorder'][n+1][0]:
                    self.list['RegionFold'].append(( self.list['RegionDisorder'][n][1]+1, self.list['RegionDisorder'][n+1][0]-1 ))
            if self.list['RegionDisorder'][-1][1] < len(sequence):
                self.list['RegionFold'].append((self.list['RegionDisorder'][-1][1]+1,len(sequence)))
        else:
            if self.list['RegionFold'][0][0] > 1:
                self.list['RegionDisorder'].append((1,self.list['RegionFold'][0][0]-1))
            for n in range(len(self.list['RegionFold'])-1):
                if self.list['RegionFold'][n][1] < self.list['RegionFold'][n+1][0]:
                    self.list['RegionDisorder'].append(( self.list['RegionFold'][n][1]+1, self.list['RegionFold'][n+1][0]-1 ))
            if self.list['RegionFold'][-1][1] < len(sequence):
                self.list['RegionDisorder'].append((self.list['RegionFold'][-1][1]+1,len(sequence)))
#########################################################################################################################
    def noDisorder(self,score=1.0):    ### Generates list of 1.0 disorder scores, defaulting to disorder.
        '''Generates list of 1.0 disorder scores, defaulting to disorder.'''
        try:
            ### Setup sequence ###
            sequence = self.info['Sequence'].upper()

            ### Generate Random Disorder ###
            self.list['ResidueDisorder'] = [score] * len(sequence)

            ### Make Regions ###
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []

            if score >= self.getStat('IUCut'):
                self.list['RegionDisorder'] = [(1,len(sequence))]
            else: self.list['RegionFold'] = [(1,len(sequence))]

            #if self.opt['PrintLog']: self.log.printLog('\r#DIS','No Disorder prediction: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),dx))
            return True
        except:
            self.log.errorLog('Error in Disorder.randomDisorder(%s)' % self.info['Name'],quitchoice=True)
            return False
#########################################################################################################################
    def randomDisorder(self):    ### Generates random disorder scores.
        '''Generates random disorder scores.'''
        try:
            ### Setup sequence ###
            sequence = self.info['Sequence'].upper()

            ### Generate Random Disorder ###
            self.list['ResidueDisorder'] = []
            for i in sequence:
                score = random.random()
                self.list['ResidueDisorder'].append(score)
            self.saveDisorder()

            ### Make Regions ###
            self.list['RegionDisorder'] = []
            self.list['RegionFold'] = []
            start = 0
            fstart = 0
            i = 0
            dx = 0
            while i < len(sequence):
                score = self.list['ResidueDisorder'][i]
                i += 1
                if not start and score > self.stat['IUCut']:    ### Start new disorder ###
                    start = i
                elif start and score <= self.stat['IUCut']:     ### End!
                    self.list['RegionDisorder'].append((start,i-1))
                    dx += i - start
                    start = 0
                if not fstart and score <= self.stat['IUCut']:    ### Start new fold ###
                    fstart = i
                elif fstart and score > self.stat['IUCut']:     ### End!
                    self.list['RegionFold'].append((fstart,i-1))
                    fstart = 0
            if start:
                self.list['RegionDisorder'].append((start,len(sequence)))
                dx += len(sequence) + 1 - start
            if fstart: self.list['RegionFold'].append((fstart,len(sequence)))
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','Random Disorder prediction complete: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),dx))
            return True
        except:
            self.log.errorLog('Error in Disorder.randomDisorder(%s)' % self.info['Name'],quitchoice=True)
            return False
#########################################################################################################################
    def globProportion(self,absolute=False):    ### Returns the proportion that is globular
        '''Returns the proportion that is globular.'''
        temp = self.list['ResidueDisorder'][0:]
        self.flatten()
        globsum = len(self.list['ResidueDisorder']) - sum(self.list['ResidueDisorder'])
        self.list['ResidueDisorder'] = temp
        if absolute: return globsum
        return float(globsum) / len(self.list['ResidueDisorder'])
#########################################################################################################################
    def flatten(self):  ### Converts RegionDisorder into "flat" 0.0 or 1.0
        '''Converts ResidueDisorder into "flat" 0.0 or 1.0.'''
        if self.opt['Flat']: return
        newlist = []
        for d in self.list['ResidueDisorder']:
            if d > self.stat['IUCut']: newlist.append(1.0)
            else: newlist.append(0.0)
        self.list['ResidueDisorder'] = newlist
        self.opt['Flat'] = True
#########################################################################################################################
    def minRegion(self):    ### Reduced self.list['RegionDisorder']/self.list['RegionFold'] using stat['MinRegion']
        '''Reduced self.list['RegionDisorder']/self.list['RegionFold'] using stat['MinRegion'].'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Smoothing') != 'sequence': return self.foldMinRegion()
            if self.stat['MinOrder'] < 0: self.stat['MinOrder'] = self.stat['MinRegion']
            if not self.list['RegionFold'] or not self.list['RegionDisorder'] or max(self.stat['MinRegion'],self.stat['MinOrder']) < 1: return
            foldbackup = self.list['RegionFold'][0:]
            disorderbackup = self.list['RegionDisorder'][0:]
            nterm = 'Disorder'                              # State of nterminus
            swap = {'Disorder':'Fold','Fold':'Disorder'}    # State swapper
            if self.list['RegionFold'][0][0] == 1: nterm = 'Fold'
            edges = []
            for region in self.list['RegionFold'] + self.list['RegionDisorder']: edges += [region[0],region[1]]
            edges.sort()
            #x#self.deBug('Ord:%s\nDis:%s' % (self.list['RegionFold'],self.list['RegionDisorder']))
            ### ~ [2] ~ Remove small regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for n in range(1,self.stat['MinRegion']):   # Work up in size
                i = 1
                while i < len(edges):
                    if (edges[i] - edges[i-1] + 1) == n:    # This region is too small
                        if i == 1: nterm = swap[nterm]; edges = [1] + edges[i+2:]
                        elif (i+2) > len(edges): edges = edges[:i-2] + edges[-1:]
                        else: edges = edges[:i-2] + edges[i+2:]
                    else: i += 2
            #x#self.deBug(' => %s' % (edges))
            ### ~ [3] ~ Regenerate Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['RegionFold'] = []
            self.list['RegionDisorder'] = []
            state = nterm
            while edges:
                self.list['Region%s' % state].append((edges[0],edges[1]))
                edges = edges[2:]
                state = swap[state]            
        except:
            self.errorLog('Problem applying MinRegion(%d)' % self.stat['MinRegion'],quitchoice=self.dev() or self.debugging())
            self.list['RegionFold'] = foldbackup
            self.list['RegionDisorder'] = disorderbackup
        #x#self.deBug('::\nOrd:%s\nDis:%s\n' % (self.list['RegionFold'],self.list['RegionDisorder']))
#########################################################################################################################
    def foldMinRegion(self):    ### Reduced self.list['RegionDisorder']/self.list['RegionFold'] using stat['MinRegion']
        '''Reduced self.list['RegionDisorder']/self.list['RegionFold'] using stat['MinRegion'].'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MinOrder'] < 0: self.stat['MinOrder'] = self.stat['MinRegion']
            if not self.list['RegionFold'] or not self.list['RegionDisorder'] or max(self.stat['MinRegion'],self.stat['MinOrder']) < 1: return
            foldbackup = self.list['RegionFold'][0:]
            disorderbackup = self.list['RegionDisorder'][0:]
            nterm = 'Disorder'                              # State of nterminus
            swap = {'Disorder':'Fold','Fold':'Disorder'}    # State swapper
            if self.list['RegionFold'][0][0] == 1: nterm = 'Fold'
            edges = []
            for state in ['Fold','Disorder']:
                for region in self.list['Region%s' % state]: edges += [(region[0],region[1],state)]
            edges.sort()
            #x#self.deBug('Ord:%s\nDis:%s' % (self.list['RegionFold'],self.list['RegionDisorder']))
            ### ~ [2] ~ Remove small regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.bugPrint('%s' % edges)
            for n in range(1,max(self.stat['MinRegion'],self.stat['MinOrder'])):   # Work up in size
                if len(edges) == 1: break
                #i# First, look at ordered regions
                if n < self.stat['MinOrder']:
                    self.bugPrint('>>> Fold:%d' % n)
                    i = 0
                    while i < len(edges):
                        if edges[i][2] == 'Fold' and (edges[i][1] - edges[i][0] + 1) == n:    # This region is too small
                            if i == 0:
                                nterm = swap[nterm]
                                edges[i+1] = (edges[i][0],edges[i+1][1],edges[i+1][2])
                                edges = edges[i+1:]
                                self.bugPrint('%s' % edges)
                            elif (i+1) == len(edges):
                                edges[i-1] = (edges[i-1][0],edges[i][1],edges[i-1][2])
                                edges = edges[:i]
                                self.bugPrint('%s' % edges)
                            else:
                                self.bugPrint('%d vs %d (%s)' % (i,len(edges),(i+1) == len(edges)))
                                self.bugPrint(edges[i-1:i+2])
                                edges[i-1] = (edges[i-1][0],edges[i+1][1],edges[i-1][2])
                                edges = edges[:i] + edges[i+2:]
                                self.bugPrint('%s' % edges)
                        else: i += 1
                #i# Next, look at disordered regions
                if n < self.stat['MinRegion']:
                    self.bugPrint('>>> Disorder:%d' % n)
                    i = 0
                    while i < len(edges):
                        if edges[i][2] == 'Disorder' and (edges[i][1] - edges[i][0] + 1) == n:    # This region is too small
                            if i == 0:
                                nterm = swap[nterm]
                                edges[i+1] = (edges[i][0],edges[i+1][1],edges[i+1][2])
                                edges = edges[i+1:]
                                self.bugPrint('%s' % edges)
                            elif (i+1) == len(edges):
                                edges[i-1] = (edges[i-1][0],edges[i][1],edges[i-1][2])
                                edges = edges[:i]
                                self.bugPrint('%s' % edges)
                            else:
                                edges[i-1] = (edges[i-1][0],edges[i+1][1],edges[i-1][2])
                                edges = edges[:i] + edges[i+2:]
                                self.bugPrint('%s' % edges)
                        else: i += 1
            #x#self.deBug(' => %s' % (edges))
            ### ~ [3] ~ Regenerate Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['RegionFold'] = []
            self.list['RegionDisorder'] = []
            state = nterm
            while edges:
                reg = edges.pop(0)
                self.list['Region%s' % state].append((reg[0],reg[1]))
                if reg[2] != state: raise ValueError
                state = swap[state]
        except:
            self.errorLog('Problem applying MinRegion(%d)' % self.stat['MinRegion'],quitchoice=self.dev() or self.debugging())
            self.list['RegionFold'] = foldbackup
            self.list['RegionDisorder'] = disorderbackup
        #x#self.deBug('::\nOrd:%s\nDis:%s\n' % (self.list['RegionFold'],self.list['RegionDisorder']))
#########################################################################################################################
### End of SECTION II: Disorder Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Disorder(mainlog,cmd_list)
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
### END OF SECTION III                                                                                                  #
#########################################################################################################################
