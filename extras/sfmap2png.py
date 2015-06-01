#!/usr/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       sfmap2png
Description:  Converts SLiMFinder Mapping files to PNGs
Version:      0.0
Last Edit:    01/09/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module takes a SLiMFinder *.mapping.fas file and uses some R visualisations to convert it into relative
    conservation/disorder PNG files.

Commandline:
    seqin=FILE  : *.mapping.fas file to use as input data []

See also rje.py generic commandline options and rje_disorder.py commands.

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
import rje, rje_disorder, rje_seq, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SFMap2PNG', '0.0', 'September 2008', '2008')
    description = 'Converts SLiMFinder Mapping files to PNGs'
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
            if rje.yesNo('Show disorder commandline options?'): out.verbose(-1,4,text=rje_disorder.__doc__)
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
### SECTION II: SFMap2PNG Class                                                                                         #
#########################################################################################################################
class SFMap2PNG(rje.RJE_Object):     
    '''
    SFMap2PNG Class. Author: Rich Edwards (2008).

    Info:str
    
    Opt:boolean

    Stat:numeric

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
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
                #self._cmdRead(cmd,type='info',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.slimJimMapping()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def slimJimMapping(self):   ### Generate SLiMJIM PNGs for all sequences
        '''Generate SpokeAln PNGs for all spokes.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapseq = {}     # Dictionary of {dataset:[seqs]}
            scmd = ['autoload=T','seqnr=F','accnr=F','replacechar=F']
            mseq = rje_seq.SeqList(self.log,self.cmd_list+scmd)     #!# Removed ['minregion=3']+ #!#
            while mseq.seq:
                ## ~ [1a] ~ Read in all sequences for one spoke ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pseq = [mseq.seq.pop(0)]  # Pseq = list of sequences for this protein
                while mseq.seq:
                    if mseq.seq[0].info['Name'].find('Motifs') > 0 and string.split(mseq.seq[0].info['Name'])[1] == 'Motifs': break   # Next protein
                    pseq.append(mseq.seq.pop(0))
                ## ~ [1b] ~ Update relevant sequence dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                mapseq[pseq[0].shortName()] = pseq[0:]
            self.printLog('#ALN','%d distinct alignments identified' % len(mapseq))
            ### ~ [2] ~ Make SLiMJIM visualisations for each protein  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0  # Number of errors
            for mapping in rje.sortKeys(mapseq):
                try:
                    ## ~ [3a] ~ Rename sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    basefile = pseq[0].shortName()
                    if self.interactive() > 0 and not rje.yesNo(basefile): continue
                    qryname = pseq[2].shortName()
                    pseq = mapseq[mapping][0:]
                    pseq[0].info['R'] = pseq[0].shortName()[len(qryname)+1:]
                    pseq[1].info['R'] = 'Masked'
                    for seq in pseq[2:]: seq.info['R'] = seq.info['ID']
                    ## ~ [3b] ~ Setup new SeqList, strip Query gaps, calculate RelCons ~~~~~~~~~~~~~~~~ ##
                    seqfile = '%s.aln.tdt' % basefile
                    if os.path.exists(seqfile): os.unlink(seqfile)
                    rseq = rje_seq.SeqList(self.log,self.cmd_list+scmd+['autoload=F'])
                    rseq.seq = pseq
                    rseq.obj['QuerySeq'] = pseq[2]
                    rseq.tidyQueryGaps()
                    rseq.saveR(rseq.seq,seqfile,name='R')
                    rseq.seq = pseq[2:]
                    relfile = '%s.rel.tdt' % basefile
                    if os.path.exists(relfile): os.unlink(relfile)
                    rseq.relCons(relfile)
                    self.deBug(rseq.obj['QuerySeq'].cmd_list)
                    ## ~ [3c] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    rcmd = '%s --no-restore --no-save --args "sfmap2png" "%s"' % (self.info['RPath'],basefile)
                    rslimjim = '%srje.r' % self.info['Path']
                    rcmd += ' < "%s" > "%s.r.tmp.txt" 2>&1' % (rslimjim,basefile)
                    self.printLog('#RSLIM',rcmd)
                    problems = os.popen(rcmd).read()
                    if problems: self.errorLog(problems,printerror=False)
                    pngx = len(glob.glob('%s*png' % basefile))
                    self.printLog('#PNG','%d PNG files made for %s' % (pngx,basefile))
                    if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
                except: self.errorLog('SLiMJIM visualisation error for "%s"' % mapping); ex += 1
            self.printLog('#SLIMJIM','Generation of SLiMJIMs complete. %d Problems.' % ex)

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: SFMap2PNG Class                                                                                  #
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
    try:SFMap2PNG(mainlog,cmd_list).run()

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
