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
# Author contact: <redwards@cabbagesofdoom.co.uk>
# School of Biotechnology and Biomolecular Sciences, University of New South Wales, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       SLiMSuite
Description:  Short Linear Motif analysis Suite
Version:      1.8.1
Last Edit:    09/05/19
Citation:     Edwards RJ & Palopoli N (2015): Methods Mol Biol. 1268:89-141. [PMID: 25555723]
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMSuite is designed to be a front end for the SLiMSuite set of sequence analysis tools. The relevant tool is given
    by the first system command, or selected using `prog=X` (or `program=X`). As much as possible, SLiMSuite will emulate
    running that tool from the commandline, adding any matching `X.ini` file to the default commandline options read in
    (*before* settings read from slimsuite.ini itself). By default, the SLiMCore tool will be called
    (libraries/rje_slimcore.py) and read in commands from slimcore.ini.

    Help for the selected tool can be accessed using the `help=T` option. Note that `-h`, `-help` or `help` alone will
    trigger the SLiMSuite help (this!). As `-help` or `help` will also set `help=T`, these commands will trigger both the
    SLiMSuite help and the selected program help (unless over-ruled by `help=F`). An explicit `help=T` command will only
    trigger the selected program help.

    Running SLiMSuite should also try importing all the main SLiMSuite modules, testing for download errors etc.

SLiMSuite tools:
    The list of tools recognised by `prog=X` will be added here as the relevant code is added:
    - SLiMCore = rje_slimcore.SLiMCore. SLiMSuite core module with MegaSLiM and UPC functions.
    - SLiMFarmer = slimfarmer.SLiMFarmer. SLiMSuite job forking/HPC controller.
    - QSLiMFinder = qslimfinder.QSLiMFinder. Query-based Short Linear Motif Finder - de novo SLiM prediction.
    - SLiMFinder = slimfinder.SLiMFinder. Short Linear Motif Finder - de novo SLiM prediction.
    - SLiMList = rje_slimlist.SLiMList. Short Linear Motif manipulation/filtering module.
    - SLiMProb = slimprob.SLiMProb. Short Linear Motif Probability - known SLiM prediction.
    - SLiMBench = slimbench.SLiMBench. SLiM discovery benchmarking module.
    - SLiMMaker = slimmaker.SLiMMaker. Simple SLiM generation from aligned peptide sequences.
    - PeptCluster = peptcluster.PeptCluster. Peptide alignment, pairwise distance and clustering tool.

Example use to run SLiMFinder:
    python SLiMSuitePATH/tools/slimsuite.py slimfinder

Please also see the SeqSuite documentation for additional utilities, which can be run from SLiMSuite or SeqSuite.

Commandline:
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    prog=X      # Identifies the tool to be used. Will load defaults from X.ini (before slimsuite.ini) [help]
    test=T/F    # Trigger test function for program X, if implemented. If `prog=test` a general test will run. [False]
    help=T/F    # Return the help documentation for program X. [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
### NOTE: INSTRUCTIONS TO ADD MODULE...
##  1. Add text, class and description to docstring above.
##  2. Add import statement for module.
##  3. Update mod dictionary
##  4. Add object creation to self.setup().
##  5. Check that the added class has a run() method.
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'extras/'))
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
sys.path.append(os.path.join(slimsuitepath,'slimsuite/extras/'))
sys.path.append(os.path.join(slimsuitepath,'slimsuite/libraries/'))
sys.path.append(os.path.join(slimsuitepath,'slimsuite/tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_slimcore, rje_slimlist
import comparimotif_V3 as comparimotif
import seqsuite
import qslimfinder, slimbench, slimfarmer, slimfinder, slimmaker, slimprob, slimmutant
import gablam, gopher, haqesac, multihaq, peptcluster
import pingu_V4 as pingu
import presto_V5 as presto
import badasp, budapest, fiesta, gasp, gfessa, picsi, unifake, seqmapper, seqforker, rje_seqgen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on SeqSuite.
    # 1.0 - Moved to tools/ for general release. Added reading and using of SeqSuite programs.
    # 1.1 - Added slimlist.
    # 1.2 - Added SLiMBench.
    # 1.3.0 - Added SLiMMaker and modified code to work with REST services.
    # 1.4.0 - Added RLC and Disorder progs to call SLiMCore. Added CompariMotif.
    # 1.5.0 - Added peptcluster and peptalign calls.
    # 1.5.1 - Changed disorder to iuscore to avoid module conflict.
    # 1.5.2 - Updated XRef REST call.
    # 1.6.0 - Removed SLiMCore as default. Default will now show help.
    # 1.7.0 - Updated to work with symbolic link in main slimsuite/ path.
    # 1.7.1 - Added error raising for protected REST alias data.
    # 1.8.0 - Added BUSCOMP and basic test function.
    # 1.8.1 - Updated documentation and added IUPred2. General tidy up and new example data for protocols paper.
    # 1.9.0 - Added SAAGA, Diploidocus, SynBad and Genomics. Slight module tidy for GitHub updates.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add menu-driven functions and menu=T/F option.
    # [ ] : Add proper test function.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SLiMSuite', '1.9.0', 'December 2020', '2014')
    description = 'Short Linear Motif analysis Suite'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',
                'Please see individual programs for citation details.',rje_obj.zen()]
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
            if rje.yesNo('Show SeqSuite commandline options?'): out.verbose(-1,4,text=seqsuite.__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(argcmd,fullcmd=True): ### Basic Setup of Program when called from commandline.
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
        cmd_list = rje.getCmdList(argcmd,info=info)         # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list,fullcmd=fullcmd) # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
mod = {'slimcore':rje_slimcore,'core':rje_slimcore,'rje_slimcore':rje_slimcore,'slimlist':rje_slimlist,'rje_slimlist':rje_slimlist,
       'slimfinder':slimfinder,'qslimfinder':qslimfinder,'slimprob':slimprob,'slimmaker':slimmaker,
       'slimfarmer':slimfarmer,'farm':slimfarmer,'slimbench':slimbench,'rlc':rje_slimcore,'iuscore':rje_slimcore,
       'comparimotif':comparimotif,'peptcluster':peptcluster,'peptalign':peptcluster}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMSuite Class                                                                                         #
#########################################################################################################################
class SLiMSuite(rje_obj.RJE_Object):
    '''
    SLiMSuite Class. Author: Rich Edwards (2014).

    Str:str
    - Name = Identifies the tool to be used. Will load defaults from X.ini (after seqsuite.ini) [seqlist]

    Bool:boolean
    - Help = Show the help documentation for program X. (Note that  [False]
    - Test = Trigger test function for program X, if implemented. If `prog=test` a general test will run. [False]

    Int:integer

    Num:float

    List:list

    Dict:dictionary

    Obj:RJE_Objects
    - Prog = Main program Object
    - ProgInfo = Info object for Program being run.
    '''
#########################################################################################################################
    def prog(self): ### Returns program name
        if self.obj['ProgInfo']: return self.obj['ProgInfo'].program
        else: return self.getStr('Name')
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Name']
        self.boollist = ['Help','Test']
        self.intlist = []
        self.numlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['Prog','ProgInfo']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        #self.setStr({'Name':'slimcore'})
        if sys.argv[1:] and sys.argv[1] in mod.keys() + seqsuite.mod.keys() + ['test']: self.setStr({'Name':sys.argv[1]})
        self.setBool({})
        self.setInt({})
        self.setNum({})
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
                self._cmdRead(cmd,type='str',att='Name',arg='prog')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='Name',arg='program')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='Name',arg='tool')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['Help','Test'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if not self.getStrLC('Name'):
            if sys.argv[1:] and sys.argv[1] not in mod.keys() + seqsuite.mod.keys() + ['help','-h','help=T']:
                self.printLog('#CORE','Argument "%s" not recognised as SLiMSuite program.' % sys.argv[1])
            #else: self.printLog('#CORE','No SLiMSuite program given: default to SLiMCore.')
            self.setStr({'Name':'help'})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,rest=False,expectpickle=False):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(rest): return False
            if self.getStrLC('Name') == 'test': return self.test()
            slimobj = self.obj['Prog']
            self.obj['Info'] = self.log.obj['Info']
            slimobj.log.obj['Info'] = info = self.obj['ProgInfo']
            slimobj.log.info['Name'] = info.program
            ### ~ [2] ~ Special REST server run code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rest:
                ## ~ [2a] Look for existing pickle: should only exist if previous JobID called ~~~~ ##
                objpickle = slimobj.unpickleMe('%s' % slimobj.basefile(runpath=True))
                if expectpickle: return objpickle
                if objpickle: slimobj = objpickle
                ## ~ [2b] Run program unless rest call is return simple information ~~~~~~~~~~~~~~~ ##
                if slimobj.getStrLC('Rest') not in ['help','version','check','outfmt']:
                    try:
                        for cmd in self.cmd_list:
                            if cmd.endswith('=!PROTECTED!'): raise ValueError('%s - check for password-protected data aliases' % cmd)
                        if self.getStrLC('Name') == 'haqesac': slimobj.run(setobjects=True)
                        elif self.getStrLC('Name') == 'xref': slimobj.run(rest=True)
                        else: slimobj.run()       #!# Note that SLiMFinder raises SystemExit upon premature ending #!#
                    except: self.errorLog('Problem with %s run' % info.program)
                    slimobj.log.endLog(info)
                    #slimobj.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
                #slimobj.log.obj['Info'] = info
                return slimobj     # Calling program should call restOutput()
            ## ~ [2a] Special test run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Test'):
                try: return slimobj.test()
                except:
                    self.printLog('#TEST','No test mode implemented for %s' % info.program)
                    return False
            ### ~ [3] ~ Regular run code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Name') == 'haqesac': slimobj.run(setobjects=True)
            else: slimobj.run()
            self.printLog('#RUN','%s V%s run finished.' % (info.program,info.version))
            slimobj.log.obj['Info'] = self.obj['Info']
            slimobj.log.info['Name'] = self.obj['Info'].program
            return slimobj
        except SystemExit:
            if rest: return slimobj
            else: raise
        except: self.errorLog(self.zen()); return False
#########################################################################################################################
    def setup(self,rest=False):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] ~ Setup Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['Prog'] = None
            prog = self.getStrLC('Name')
            if prog in mod:
                i = self.obj['ProgInfo'] = mod[prog].makeInfo()
                self.printLog('#PROG','%s V%s: %s' % (i.program,i.version,i.description))
                progcmd = rje.getCmdList(['basefile=%s' % prog],info=i) + self.cmd_list + ['newlog=F']
                out = rje.Out(cmd_list=progcmd)
                out.printIntro(i)
                if self.getBool('Help'): progcmd = mod[prog].cmdHelp(i,out,['help']+progcmd)
                purgelist = seqsuite.purgelist
                self.printLog('#CMD','Full %s CmdList: %s' % (i.program,rje.argString(rje.tidyArgs(progcmd,nopath=self.getStrLC('Rest') and not self.dev(),purgelist=purgelist))),screen=False)
                #self.debug(prog)
            ## ~ [1a] ~ Make self.obj['Prog'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if prog in ['slimcore','rje_slimcore','core']: self.obj['Prog'] = rje_slimcore.SLiMCore(self.log,progcmd)
                elif prog in ['rlc','iuscore']: self.obj['Prog'] = rje_slimcore.SLiMCore(self.log,progcmd+['prog=%s' % prog])
                elif prog in ['slimlist','rje_slimlist']: self.obj['Prog'] = rje_slimlist.SLiMList(self.log,progcmd)
                elif prog in ['slimfinder']: self.obj['Prog'] = slimfinder.SLiMFinder(self.log,progcmd)
                elif prog in ['qslimfinder']: self.obj['Prog'] = qslimfinder.QSLiMFinder(self.log,progcmd)
                elif prog in ['slimprob']: self.obj['Prog'] = slimprob.SLiMProb(self.log,progcmd)
                elif prog in ['slimmaker']: self.obj['Prog'] = slimmaker.SLiMMaker(self.log,progcmd)
                elif prog in ['slimfarmer','farm']: self.obj['Prog'] = slimfarmer.SLiMFarmer(self.log,progcmd)
                elif prog in ['slimbench']: self.obj['Prog'] = slimbench.SLiMBench(self.log,progcmd)
                elif prog in ['comparimotif']: self.obj['Prog'] = comparimotif.CompariMotif(self.log,progcmd)
                elif prog in ['peptcluster']: self.obj['Prog'] = peptcluster.PeptCluster(self.log,progcmd)
                elif prog in ['peptalign']: self.obj['Prog'] = peptcluster.PeptCluster(self.log,['peptalign=T']+progcmd+['peptdis=None'])
                self.obj['Prog'].dict['Output']['help'] = mod[prog].__doc__
            elif prog in seqsuite.mod:
                seqsuiteobj = seqsuite.SeqSuite(self.log,self.cmd_list)
                self.obj['Prog'] = seqsuiteobj.setup()
                self.obj['ProgInfo'] = seqsuiteobj.obj['ProgInfo']
                self.obj['Prog'].dict['Output']['help'] = seqsuite.mod[prog].__doc__
            elif prog == 'test':
                self.obj['Prog'] = self
                self.obj['ProgInfo'] = self.obj['ProgInfo']
                self.obj['Prog'].dict['Output']['help'] = __doc__

            ### ~ [2] ~ Failure to recognise program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['Prog']:
                if self.getStrLC('Name') != 'help':
                    if not rest: self.printLog('#ERR','Program "%s" not recognised.' % self.getStr('Name'))
                    if self.i() < 0 or rest: return False
                #!# Try SeqSuite? #!#
                if self.getStrLC('Name') == 'help' or rje.yesNo('Show SLiMSuite help with program options?'):
                    extracmd = cmdHelp(cmd_list=['help'])[1:]
                    if extracmd:
                        self.cmd_list += extracmd
                        self._cmdList()
                        if prog != self.getStrLC('Name'): return self.setup()
                self.setStr({'Name':rje.choice('Give program name (Blank or CTRL+C to quit)')})
                if self.getStrLC('Name'): return self.setup()
                else: return False
            return self.obj['Prog']     # Setup successful
        except KeyboardInterrupt: return False
        except SystemExit: raise
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### SLiMSuite Testing Methods                                                                               #
#########################################################################################################################
    def test(self):      ### Generic method
        '''
        SLiMSuite Testing method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#TEST','SLiMSuite modules imported OK.')
            #import comparimotif_V3 as comparimotif
            #import qslimfinder, slimbench, slimfarmer, slimfinder, slimmaker, slimprob, slimmutant
            #import gablam, gopher, haqesac, multihaq, peptcluster
            #import pingu_V4 as pingu
            #import presto_V5 as presto
            #import badasp, budapest, fiesta, gasp, gfessa, picsi, unifake, seqmapper, seqforker, rje_seqgen
            return
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################
### End of SECTION II: SLiMSuite Class                                                                                  #
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
def checkProg(prog): return prog in mod or prog in seqsuite.mod
#########################################################################################################################
def rest(argcmd,expectpickle=False):
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram(argcmd,fullcmd=False) # Commands coming from system arguments or REST call
    except: return 'Unexpected error during program setup: %s' % sys.exc_info()[0]

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: return SLiMSuite(mainlog,cmd_list).run(rest=True,expectpickle=expectpickle)

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except: return 'Fatal error in main %s run.' % info.program
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram(sys.argv[1:])
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: SLiMSuite(mainlog,cmd_list).run()

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
