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
Version:      0.8
Last Edit:    06/08/14
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module currently has limited function and no standalone capability, though this may be added with time. It is
    designed for use with other modules. The disorder Class can be given a sequence and will run the appropriate
    disorder prediction software and store disorder prediction results for use in other programs. The sequence will have
    any gaps removed.

    Currently four disorder prediction methods are implemented:
    * IUPred : Dosztanyi Z, Csizmok V, Tompa P & Simon I (2005). J. Mol. Biol. 347, 827-839. This has to be installed
    locally. It is available on request from the IUPred website and any use of results should cite the method. (See
    http://iupred.enzim.hu/index.html for more details.) IUPred returns a value for each residue, which by default,
    is determined to be disordered if > 0.5.
    * FoldIndex : This is run directly from the website (http://bioportal.weizmann.ac.il/fldbin/findex) and more simply
    returns a list of disordered regions. You must have a live web connection to use this method!
    * ANCHOR : Meszaros B, Simon I & Dosztanyi Z (2009). PLoS Comput Biol 5(5): e1000376. This has to be installed
    locally. It is available on request from the ANCHOR website and any use of results should cite the method. (See
    http://anchor.enzim.hu/ for more details.) ANCHOR returns a probability value for each residue, which by default,
    is determined to be disordered if > 0.5.
    * Parse: Parsed disorder from protein sequence name, e.g. DisProt download.
    #X-Y = disordered region; &X-Y = ordered region [0.0]

    For IUPred, the individual residue results are stored in Disorder.list['ResidueDisorder']. For both methods, the
    disordered regions are stored in Disorder.list['RegionDisorder'] as (start,stop) tuples.
    
Commandline:
    ### General Options ###
    disorder=X  : Disorder method to use (iupred/foldindex/anchor/parse) [iupred]
    iucut=X     : Cut-off for IUPred/ANCHOR results [0.2]
    iumethod=X  : IUPred method to use (long/short) [short]
    sequence=X  : Sequence to predict disorder for (autorun) []
    name=X      : Name of sequence to predict disorder for []
    minregion=X : Minimum length of an ordered/disordered region [0]

    ### System Settings ###
    iupath=PATH : The full path to the IUPred executable [c:/bioware/iupred/iupred.exe]
    anchor=PATH : Full path to ANCHOR executable []
    filoop=X    : Number of times to try connecting to FoldIndex server [10]
    fisleep=X   : Number of seconds to sleep between attempts [2]
    iuchdir=T/F : Whether to change to IUPred directory and run (True) or rely on IUPred_PATH env variable [False]

Uses general modules: copy, os, string, sys, time, urllib2
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time, urllib2
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
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add region buffer - min length of region allowed before ignoring
    # [Y] : Add a PrintLog option that controls whether printing to Log or not.
    # [ ] : Neaten and tidy
    # [ ] : Add domain-based disorder stuff?
    # [ ] : Add option to store/re-read prediction data for later use.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_DISORDER', '0.8', 'August 2014', '2008')
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
    
    Opt:boolean
    - Flat = whether ResidueDisorder is a "flat" 1/0 or graded (e.g. raw IUPred)
    - IUChDir = Whether to change to IUPred directory and run (True) or rely on IUPred_PATH env variable [False]
    - PrintLog = whether to print disorder prediction status to Log [False]

    Stat:numeric
    - IUCut = Cut-off for IUPred results [0.2]
    - FILoop = Number of times to try connecting to FoldIndex server [10]
    - FISleep = Number of seconds to sleep between attempts [2]
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
        self.infolist = ['Sequence','Disorder','IUPath','IUMethod','ANCHOR']
        self.optlist = ['Flat','PrintLog','IUChDir']
        self.statlist = ['IUCut','FILoop','FISleep','MinRegion']
        self.listlist = ['ResidueDisorder','RegionDisorder','RegionFold']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='',opt=False,stat=0.2,obj=None,setlist=True,setdict=True)
        self.setInfo({'IUPath':rje.makePath('c:/bioware/iupred/iupred.exe',wholepath=True),'IUMethod':'short',
                      'Disorder':'iupred'})
        self.setStat({'FILoop':10,'FISleep':2,'MinRegion':0,'IUCut':0.2})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'info',['Disorder','Name','Sequence','IUMethod'])
                self._cmdRead(cmd,type='stat',att='IUCut')
                self._cmdRead(cmd,type='fullpath',att='IUPath')
                self._cmdRead(cmd,type='fullpath',att='ANCHOR')
                self._cmdReadList(cmd,'int',['FILoop','FISleep','MinRegion'])
                self._cmdReadList(cmd,'opt',['IUChDir'])                                        
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        ### AutoRun ###
        if self.info['Sequence'].lower() not in ['none','']: self.disorder(self.info['Sequence'])
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
            self.info['Sequence'] = string.join(string.split(sequence,'-'),'')
            if name: self.info['Name'] = name
            if not self.info['Sequence']:
                self.log.errorLog('Cannot calculate disorder: no AAs given in sequence!',printerror=False)
                return False
            ### ~ [2] ~ Run Disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Disorder'].lower() == 'iupred': return self.iuPred()
            elif self.info['Disorder'].lower() == 'foldindex': return self.foldIndex()
            elif self.info['Disorder'].lower() == 'anchor': return self.ANCHOR()
            elif self.info['Disorder'].lower() == 'parse': return self.parseDisorder()
            else:
                self.log.errorLog('Cannot calculate disorder: no disorder method given!',printerror=False)
                return False
        except:
            self.log.errorLog('Error in Disorder.disorder(%s)' % name,quitchoice=True)
            return False
#########################################################################################################################
    def iuPred(self,retry=2):     ### Runs IUPred disorder prediction
        '''Runs IUPred disorder prediction.'''
        mydir = os.path.abspath(os.curdir)
        try:
            ### Setup sequence and temp file ###
            sequence = self.info['Sequence'].upper()
            name = self.info['Name'][:4] + rje.randomString(8)
            tmp = name + '.tmp'
            
            ### Run Disorder ###
            iupath = string.join(string.split(self.info['IUPath'],os.sep)[:-1],os.sep)
            iupred = string.split(self.info['IUPath'],os.sep)[-1]
            if self.opt['IUChDir']: os.chdir(string.join(string.split(self.info['IUPath'],os.sep)[:-1],os.sep))
            open(tmp,'w').write('>%s\n%s\n' % (name,sequence))
            if self.opt['IUChDir'] and self.opt['Win32']: iucmd = '%s %s %s' % (iupred,tmp,self.info['IUMethod'].lower())
            elif self.opt['IUChDir']: iucmd = './%s %s %s' % (iupred,tmp,self.info['IUMethod'].lower())
            else: iucmd = '%s %s %s' % (self.info['IUPath'],tmp,self.info['IUMethod'].lower())
            dlines = os.popen(iucmd).readlines()
            try: os.unlink(tmp)
            except: self.errorLog('Cannot delete %s!' % tmp)
            if self.opt['IUChDir']: os.chdir(mydir)            
            if self.info['Name'] not in ['','None']: name = self.info['Name']
            self.list['ResidueDisorder'] = []
            for d in dlines:
                if rje.matchExp('^\s*(\d+)\s+(\S)\s+(\S+)',d):
                    dm = rje.matchExp('^\s*(\d+)\s+(\S)\s+(\S+)',d)
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

            ### Run Disorder ###
            retry = self.stat['FILoop']
            url = "http://bioportal.weizmann.ac.il/fldbin/findex"
            params = "m=xml&sq=" + sequence  + "  " 
            while retry:
                try:
                    flines = urllib2.urlopen(url, params).readlines()
                except:
                    flines = []
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
        '''
        try:
            ### Setup sequence and name ###
            sequence = self.info['Sequence']
            name = self.info['Name']
            self.list['ResidueDisorder'] = [0.5] * len(sequence)
            self.list['RegionDisorder'] = []
            scoredict = {'#':1.0,'&':0.0}

            ### Process ###
            for region in string.split(name)[1:]:
                if rje.matchExp('^[#&](\d+)-(\d+)',region):
                    (i,x,y) = rje.matchExp('^([#&])(\d+)-(\d+)',region)
                    score = scoredict[i]
                    start = string.atoi(x) - 1
                    end = string.atoi(y)
                    for r in range(start,end): self.list['ResidueDisorder'][r] = score
                    if i == '#': self.list['RegionDisorder'].append((start,end))
            self.minRegion()
            if self.opt['PrintLog']: self.log.printLog('\r#DIS','DisProt Disorder parsing complete: %d disorder regions, %d disordered aa' % (len(self.list['RegionDisorder']),self.list['ResidueDisorder'].count(1.0)))
            return True
        except:
            self.log.errorLog('Error in Disorder.foldIndex(%s)' % self.info['Name'],quitchoice=True)
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
            if not self.list['RegionFold'] or not self.list['RegionDisorder'] or self.stat['MinRegion'] < 1: return 
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
