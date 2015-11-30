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
Module:       scap
Description:  Sequence Composition Assessment of Peptides
Version:      0.0
Last Edit:    05/01/09
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    seqin=FILE      : File with input sequences for testing [None]
    xmerback=FILE   : File with sequences for background Markov Chain calculations (can be pickle) [seqin]    
    scapback=FILE   : File with sequences for background SCAP calculations [None]    
    sorted=T/F      : Whether to use sorted xmer fragments to reduce memory requirements [False]
    xmers=X[,Y]     : Deal with Xmers [upto Ymers] [1]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_markov, rje_seq, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add expectation calculated from actual numbers, rather than simply the mean of 0.05.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SCAP', '0.0', 'January 2009', '2009')
    description = 'Sequence Composition Assessment of Peptides'
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
### SECTION II: SCAP Class                                                                                              #
#########################################################################################################################
class SCAP(rje.RJE_Object):     
    '''
    SCAP Class. Author: Rich Edwards (2009).

    Info:str
    - ScapBack = File with sequences for background SCAP calculations [None]    
    - XmerBack = File with sequences for background Markov Chain calculations (can be pickle) [seqin]    
    
    Opt:boolean

    Stat:numeric

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - Markov = rje_markov.Markov object for background Markov Chain calculations     
    - ScapBack = SeqList object for background SCAP calculations [xmerback or seqin]    
    - SeqList = SeqList object with input sequences for testing [None]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['ScapBack','XmerBack']
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['Markov','ScapBack','SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['ScapBack','XmerBack'])  # No need for arg if arg = att.lower()
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.scap()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = self.obj['SeqList'] = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list)   #!# Add code for memsaver/autoload=F #!#
            self.printLog('#SCAP','%s sequences loaded for SCAP analysis' % rje.integerString(seqlist.seqNum()))
            ## ~ [1b] ~ Xmer background file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mseqfile = self.info['XmerBack']
            if mseqfile.lower() in ['','none']: mseqfile = self.info['XmerBack'] = seqlist.info['Name']
            markov = self.obj['Markov'] = rje_markov.Markov(self.log,['autoload=T','accnr=F','seqnr=F']+self.cmd_list+['seqin=%s' % mseqfile,'direction=both','markov=F','scap=T'])
            markov.setup()
            maxx = markov.stat['MaxXmer']
            if self.info['Basefile'].lower() in ['','none']:
                self.info['Basefile'] = '%s.scap' % rje.baseFile(seqlist.info['Name'],True)
                if markov.opt['Sorted']: self.info['Basefile'] = '%s.sorted' % self.info['Basefile']
            basefile = self.info['Basefile']
            self.printLog('#MARKOV','Markov setup complete')
            ## ~ [1c] ~ SCAP Background file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scapfile = self.info['ScapBack']
            if scapfile.lower() in ['','none',seqlist.info['Name'].lower()]: self.obj['ScapBack'] = self.obj['SeqList']
            elif scapfile == mseqfile: self.obj['ScapBack'] = markov.obj['SeqList'] 
            else: self.obj['ScapBack'] = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list+['seqin=%s' % scapfile])
            self.printLog('#SCAP','%s sequences for SCAP Background' % rje.integerString(seqlist.seqNum()))

            ### ~ [2] Markov Chains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if mseqfile == seqlist.info['Name']: markov.obj['SeqList'] = seqlist
            elif mseqfile == self.obj['ScapBack'].info['Name']: markov.obj['SeqList'] = self.obj['ScapBack']
            mpickle = markov.unpickleMe()
            if mpickle: markov = self.obj['Markov'] = mpickle
            if not markov.suftree() or not markov.pretree() or maxx > markov.stat['MaxXmer']:
                markov.run()
                markov.pickleMe()
            markov.opt['DeBug'] = self.opt['DeBug']
            self.deBug(markov.opt)
            self.deBug(markov.stat)
            #self.deBug(markov.suftree())
            #self.deBug(markov.pretree())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### SCAP Methods                                                                                            #
#########################################################################################################################
    def scap(self):     ### Full SCAP method
        '''Full SCAP method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            markov = self.obj['Markov']
            minx = markov.stat['MinXmer']
            maxx = markov.stat['MaxXmer']
            headers = ['seq','type','sorted']
            for x in range(minx,maxx+1): headers.append('X%d' % x)
            delimit = rje.getDelimit(self.cmd_list,'\t')
            scapfile = '%s.%s' % (self.info['Basefile'],rje.delimitExt(delimit))
            rje.delimitedFileOutput(self,scapfile,headers,delimit,rje_backup=True)
            ### ~ [2] SCAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (sx,stot) = (0.0,self.obj['SeqList'].seqNum())
            for seq in self.obj['SeqList'].seq:
                self.progLog('\r#SCAP','SCAP processing Query to %s: %.2f%%' % (scapfile,(sx/stot))); sx += 100.0
                datadict = {'seq':seq.shortName(),'type':'qry','sorted':markov.opt['Sorted']}
                for x in range(minx,maxx+1): 
                    datadict['X%d' % x] = self.scapSeq(seq.info['Sequence'],x)
                    if datadict['X%d' % x] > 0.001: datadict['X%d' % x] = '%.4f' % datadict['X%d' % x]
                    else: datadict['X%d' % x] = '%.3e' % datadict['X%d' % x]
                rje.delimitedFileOutput(self,scapfile,headers,delimit,datadict)
            self.printLog('\r#SCAP','SCAP processed Query to %s for %s sequences.' % (scapfile,rje.integerString(stot)))
            ## ~ [2b] Background ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.obj['ScapBack'] != self.obj['SeqList']:
                (sx,stot) = (0.0,self.obj['ScapBack'].seqNum())
                for seq in self.obj['ScapBack'].seq:
                    self.progLog('\r#SCAP','SCAP processing Background to %s: %.2f%%' % (scapfile,(sx/stot))); sx += 100.0
                    datadict = {'seq':seq.shortName(),'type':'bg','sorted':markov.opt['Sorted']}
                    for x in range(minx,maxx+1):
                        datadict['X%d' % x] = self.scapSeq(seq.info['Sequence'],x)
                        if datadict['X%d' % x] > 0.001: datadict['X%d' % x] = '%.4f' % datadict['X%d' % x]
                        else: datadict['X%d' % x] = '%.3e' % datadict['X%d' % x]
                    rje.delimitedFileOutput(self,scapfile,headers,delimit,datadict)
                self.printLog('\r#SCAP','SCAP processed Background to %s for %s sequences.' % (scapfile,rje.integerString(stot)))
            if markov.opt['Sorted']: self.printLog('#SCAP','Sorted SCAP run complete')
            else: self.printLog('#SCAP','UnSorted SCAP run complete')
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def scapSeq(self,sequence,maxx=0):    ### SCAP method for a single sequence
        '''SCAP method for a single sequence'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            markov = self.obj['Markov']
            if not maxx: maxx = markov.stat['MaxXmer']
            plist = []
            ### ~ [2] Calculate SCAP probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for r in range(len(sequence)):
                pref = sequence[max(0,r+1-maxx):r+1]    # Should include residue r!
                suff = sequence[r:r+maxx]               # Should include residue r!
                #X#print maxx, pref, suff
                p = self.xmerProb(pref)
                if maxx > 1: p *= self.xmerProb(suff,prefix=True) 
                self.deBug(p)
                plist.append(p)
            self.deBug('%s :: %s' % (plist,rje.geoMean(plist)))
            return rje.geoMean(plist)
            return rje.meansd(plist)[0]
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def xmerProb(self,xmer,prefix=False):   ### Returns SCAP probability for count for given xmer from markov.dict tree
        '''
        Returns count for given xmer from self.dict tree.
        >> xmer:str = Xmer of interest
        >> prefix:bool [False] = Use Prefix tree rather than suffix tree
        '''
        ### ~ [1] ~ Choose tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        markov = self.obj['Markov']
        if markov.opt['Sorted'] and len(xmer) < markov.stat['MaxXmer']: return 1.0
        if prefix: _sufdic = markov.pretree(); xmer = rje.strReverse(xmer)
        else: _sufdic = markov.suftree()
        if markov.opt['Sorted']: xmer = rje.strSort(xmer[:-1]) + xmer[-1]
        prex = 0
        self.deBug('%s :: %s' % (xmer,rje.sortKeys(_sufdic)))
        ### ~ [2] ~ Find subtree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for x in range(len(xmer)):
            if xmer[x] in _sufdic.keys():
                if not _sufdic.has_key('e'):
                    _sufdic['e'] = 0.0
                    for a in markov.list['Alphabet']:
                        if _sufdic.has_key(a):
                            fa = _sufdic[a]['='] / float(_sufdic['='])
                            _sufdic['e'] += (fa * fa)
                prex = _sufdic['='] * _sufdic['e']; _sufdic = _sufdic[xmer[x]]
            elif xmer[x] not in markov.list['Alphabet'] or not prex: return 1.0
            else: return 0.5 / (0.05 * prex)    # Arbitrary small number!
            self.deBug('%s %d [%s] :: %d = %s' % (xmer,x,xmer[x],prex,_sufdic))
        ### ~ [3] ~ Calculate SCAP value ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if prex: return _sufdic['='] / prex
        else: return 1.0
#########################################################################################################################
### End of SECTION II: SCAP Class                                                                                       #
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
    try: SCAP(mainlog,cmd_list).run()

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
