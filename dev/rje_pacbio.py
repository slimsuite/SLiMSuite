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
Module:       rje_pacbio
Description:  Miscellaneous Utilities for PacBio Sequencing
Version:      1.1.0
Last Edit:    27/05/15
Webserver:    http://www.slimsuite.unsw.edu.au/servers/pacbio.php
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module estimates the % genome coverage and accuracy for different X coverage of a genome using PacBio
    sequencing, i.e. assuming a non-biased error distribution. Calculations use binomial/poisson distributions, assuming
    independence of sites. Accuracy is based on >50% reads covering a particular base
    having the correct call. Assuming random calls at the other positions, 25% of the "wrong" positions will be correct
    by chance. In reality, it will be even higher than this, assuming majority calls are used. Wrong calls will be split
    between three possible incorrect bases. Accuracy is therefore a conservative estimate.

    All calculations are based on *assembled* reads, and therefore using the full `smrtreads=X` value for SMRT cells will
    overestimate coverage. Note that `smrtreads=X` can be used to input sequence capacity in Gb (or Mb) rather than read
    counts by changing `smrtunits=X`.

Output:
    Main output is a results table containing the following fields:
    * XCoverage = estimated average X genome coverage.
    * SMRT = estimated number of SMRT cells.
    * %Coverage = estimated percentage genome coverage.
    * %Accuracy = estimated percentage of covered bases with correct base calls.
    * %Xn = 0+ columns giving % sites with coverage >= Xn (`xnlist=LIST`).

Commandline:
    ### ~ Genome Coverage Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    coverage=T/F    : Whether to generate coverage report [True]
    avread=X        : Average read length (bp) [20000]
    errperbase=X    : Error-rate per base [0.14]
    genomesize=X    : Genome size (bp) [4e9]
    maxcov=X        : Maximmum X coverage to calculate [100]
    smrtreads=X     : Average assemble output of a SMRT cell [50000]
    smrtunits=X     : Units for smrtreads=X (reads/Gb/Mb) [reads]
    bysmrt=T/F      : Whether to output estimated  coverage by SMRT cell rather than X coverage [False]
    xnlist=LIST     : Additional columns giving % sites with coverage >= Xn [10,25,50,100]
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
import rje, rje_db, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 1.0.0 - Initial working version for server.
    # 1.1.0 - Added xnlist=LIST : Additional columns giving % sites with coverage >= Xn [10,25,50,100].
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [ ] : Improved error estimation.
    # [Y] : Option to do a per-SMRT cell analysis.
    # [ ] : Add costing.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('PacBio', '1.1.0', 'May 2015', '2015')
    description = 'Generic RJE Module'
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
### SECTION II: PacBio Class                                                                                            #
#########################################################################################################################
class PacBio(rje_obj.RJE_Object):
    '''
    PacBio Class. Author: Rich Edwards (2015).

    Str:str
    - SMRTUnits=X     : Units for smrtreads=X (reads/Gb/Mb) [reads]

    Bool:boolean
    - BySMRT=T/F      : Whether to output estimated by SMRT cell rather than X coverage [False]
    - Coverage=T/F    : Whether to generate coverage report [True]

    Int:integer
    - MaxCov=X        : Maximmum X coverage to calculate [100]

    Num:float
    - AvRead=X        : Average read length (bp) [20000]
    - ErrPerBase=X    : Error-rate per base [0.14]
    - GenomeSize=X    : Genome size (bp) [4e9]
    - SMRTReads=X     : Average number of reads of a SMRT cell [50000]

    File:file handles with matching str filenames
    
    List:list
    - Accuracy : list of %accuracy for each xcoverage level
    - XnList=LIST     : Additional columns giving % sites with coverage >= Xn [10,25,50,100]

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['SMRTUnits']
        self.boollist = ['BySMRT','Coverage']
        self.intlist = ['MaxCov']
        self.numlist = ['AvRead','ErrPerBase','GenomeSize','SMRTReads']
        self.filelist = []
        self.listlist = ['Accuracy','XnList']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'SMRTUnits':'reads'})
        self.setBool({'BySMRT':False,'Coverage':True})
        self.setInt({'MaxCov':100})
        self.setNum({'AvRead':20000,'ErrPerBase':0.14,'GenomeSize':4e9,'SMRTReads':50000})
        self.basefile('pacbio')
        self.list['XnList'] = [10,25,50,100]
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
                self._cmdReadList(cmd,'str',['SMRTUnits'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BySMRT','Coverage'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxCov'])   # Integers
                self._cmdReadList(cmd,'float',['AvRead','ErrPerBase','GenomeSize','SMRTReads']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['XnList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Coverage'): self.coverage()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.db().basefile(self.basefile())
            self.list['Accuracy'] = [0,1.0 - self.getNum('ErrPerBase')]
            ## ~ [1a] SMRTReads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while self.getStrLC('SMRTUnits') not in ['reads','gb','mb']:
                txt = 'SMRTUnits "%s" not recognised'
                if self.getNum('SMRTReads') < 10: smrtunits = 'Gb'
                elif self.getNum('SMRTReads') > 10000: smrtunits = 'reads'
                else: smrtunits = 'Mb'
                if self.i() < 0 or rje.yesNo('%s: switch to (%s) %s?' % (txt,self.getNum('SMRTReads'),smrtunits)):
                    self.setStr({'SMRTUnits':smrtunits})
                elif self.i() >0: self.setStr({'SMRTUnits':rje.choice('SMRTUnits (reads/Gb/Mb)?')})
                self.printLog('#UNITS','%s => %s' % (txt,self.getStr('SMRTUnits')))
            if self.getStrLC('SMRTUnits') in ['gb','mb']:
                smrttotal = self.getNum('SMRTReads') * {'gb':1e9,'mb':1e6}[self.getStrLC('SMRTUnits')]
                txt =  '%s %s @ %.3f kb/read' % (self.getNum('SMRTReads'),self.getStr('SMRTUnits'),self.getNum('AvRead')/1000.0)
                self.setNum({'SMRTReads':smrttotal/self.getNum('AvRead')})
                txt += ' => %s reads' % rje.iStr(int(self.getNum('SMRTReads')))
                self.printLog('#READS',txt)
            ## ~ [1b] XnList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xnlist = []
            for xn in self.list['XnList']:
                if xn == '': continue
                try:
                    ixn = int(xn)
                    if xn not in [ixn,'%d' % ixn]: self.printLog('#XN','"%s" -> %dX' % (xn,ixn))
                    if ixn == 0: self.printLog('#XN','No point in 0X output: use 1-%Coverage.')
                    elif ixn == 1: self.printLog('#XN','No point in 1X output: use %Coverage.')
                    else: xnlist.append(ixn)
                except: self.errorLog('Could not process %s as part of XnList. (Integers only.)' % xn)
            xnlist.sort()
            if xnlist: self.printLog('#XN','XnList: %sX.' % string.join(string.split('%s' % xnlist,','),'X, ')[1:-1])
            self.list['XnList'] = xnlist
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
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
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def coverage(self): ### Calculates estimated % coverage and accuracy of genome sequencing.
        '''Calculates estimated % coverage and accuracy of genome sequencing.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # XCoverage, SMRT, %Cov, Accuracy
            if self.getBool('BySMRT'): ckey = 'SMRT'
            else: ckey = 'XCoverage'
            cfields = ['XCoverage','SMRT','%Coverage','Accuracy']
            for xn in self.list['XnList']: cfields.append('%%X%d' % xn)
            cdb = self.db().addEmptyTable('coverage',cfields,[ckey])

            ### ~ [2] Calculate stats for one round ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('\r#XCOV','Calculating coverage stats...')
            cov_per_base_per_read = self.getNum('AvRead') / self.getNum('GenomeSize')
            if self.getBool('BySMRT'): reads = self.getInt('SMRTReads')    # If going per SMRT cell
            else: reads = int(0.5 + self.getNum('GenomeSize') / self.getNum('AvRead')) # if going per X coverage
            # Calculate X coverage counts using binomial
            bases = int(self.getNum('GenomeSize'))
            xcov = []   # List where index is X coverage and number is proportion of reads
            while bases > 1:
                try: xcov.append(rje.logBinomial(len(xcov),reads,cov_per_base_per_read,exact=True,callobj=self))
                except: xcov.append(rje.logPoisson(len(xcov),reads*cov_per_base_per_read,exact=True,callobj=self))
                bases -= self.getNum('GenomeSize') * xcov[-1]
                if len(xcov) > reads: raise ValueError('XCoverage cannot exceed read count!')
            cyccov = xcov[0:]
            self.debug(xcov)

            ### ~ [3] Cycle through rounds, multiplying by X coverage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cx = 0.0; ctot = self.getInt('MaxCov'); xcoverage = 0.0
            while xcoverage < self.getInt('MaxCov'):
                self.progLog('\r#XCOV','Calculating coverage stats: %.1f%% (%d|%d)' % ((cx/ctot,cdb.entryNum()+1,len(cyccov)))); cx += 100.0
                # Update xcov: calculate %bases at different X coverage
                if cdb.entryNum():  # Equivalent of starting with [1.0] (!00% 0 @ 0X)
                    prevcov = cyccov[0:]
                    cyccov = [0.0] * (len(prevcov)*2 - 1)
                    for xi in range(len(prevcov)):
                        for xj in range(len(xcov)):
                            x = xi + xj
                            cyccov[x] += (prevcov[xi] * xcov[xj])
                while(cyccov[-1]) < 1.0 / self.getNum('GenomeSize'): cyccov.pop(-1)
                # Calculate accuracy: For each X coverage, calculate % bases with >50% correct
                accuracy = 0.0
                for x in range(len(cyccov[1:])): accuracy += cyccov[x] * self.accuracy(x)
                accuracy = 100.0 * accuracy / sum(cyccov[1:])
                # SMRT cells versus coverage
                xcoverage += self.getNum('AvRead') * reads / self.getNum('GenomeSize')
                smrt = (self.getNum('GenomeSize') * xcoverage) / (self.getNum('AvRead') * self.getNum('SMRTReads'))
                # Update cdb
                #centry = {'XCoverage':'%.3f' % xcoverage,'SMRT':'%.2f' % smrt,'%Coverage':100.0 * (1.0-cyccov[0]),'Accuracy':accuracy}
                centry = {'XCoverage':rje.sf(xcoverage,3),'SMRT':rje.sf(smrt,3),'%Coverage':100.0 * (1.0-cyccov[0]),'Accuracy':accuracy}
                for xn in self.list['XnList']:
                    if xn <= len(cyccov): centry['%%X%d' % xn] = rje.sf(100.0*sum(cyccov[xn:]),4)
                    else: centry['%%X%d' % xn] = 0.000
                cdb.addEntry(centry)
            self.progLog('\r#XCOV','Calculated coverage stats upto %dX coverage.' % self.getInt('MaxCov'))

            ### ~ [4] Save results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for xkey in cdb.dataKeys():
                cdb.dict['Data'][float(xkey)] = cdb.dict['Data'].pop(xkey)
            cdb.saveToFile()

            return
        except: self.errorLog('%s.coverage error' % self.prog())
#########################################################################################################################
    def accuracy(self,xcoverage):   ### Calculate accuracy (if required) at xcoverage and returns
        '''Calculate accuracy (if required) at xcoverage and returns.'''
        try:### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while len(self.list['Accuracy']) <= xcoverage:
                if not int((1.0 - self.list['Accuracy'][-1]) * self.getNum('GenomeSize')):
                    self.list['Accuracy'].append(1.0)
                    continue
                xcov = len(self.list['Accuracy'])
                majority = int(xcov/2.0) + 1        # Number of correct reads needed for majority
                try: self.list['Accuracy'].append(rje.logBinomial(majority,xcov,1.0 - self.getNum('ErrPerBase'),exact=False,callobj=self))
                except: self.list['Accuracy'].append(rje.logPoisson(majority,xcov*(1.0 - self.getNum('ErrPerBase')),exact=False,callobj=self))
                self.debug(self.list['Accuracy'])
            return self.list['Accuracy'][xcoverage]
        except: self.errorLog('%s.accuracy error' % self.prog())
#########################################################################################################################
### End of SECTION II: PacBio Class                                                                                     #
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
    try: PacBio(mainlog,cmd_list).run()

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
