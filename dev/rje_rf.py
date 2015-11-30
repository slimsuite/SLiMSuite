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
Module:       rje_rf
Description:  RF prediction play module
Version:      0.0
Last Edit:    16/11/10
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is not for public consumption.

Commandline:

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
import rje, rje_seq, rje_sequence, rje_zen
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
    (program, version, last_edit, copyright) = ('RJE_RF', '0.0', 'November 2010', '2010')
    description = 'Generic RJE Module'
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class RF(rje.RJE_Object):     
    '''
    RF Class. Author: Rich Edwards (2010).

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
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list)
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
            return
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
    def rfAtt(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rfhead = ['Att','RF1','RF2','RF3','RF-1','RF-2','RF-3','ObsRF1','ObsRF2','ObsRF3','ObsRF-1','ObsRF-2','ObsRF-3','ExpRF1','ExpRF2','ExpRF3','ExpRF-1','ExpRF-2','ExpRF-3']
            rfdata = {}; rfobs = {}; rfexp = {}; ntfreq = {}
            for rf in ['RF1','RF2','RF3','RF-1','RF-2','RF-3']:
                rfdata[rf] = {}; rfobs[rf] = {}; rfexp[rf] = {}
                for x in rje_seq.alph_protx[:-1] + ['*']: rfdata[rf][x] = 0; rfobs[rf][x] = 0; rfexp[rf][x] = 0
                for a1 in rje_seq.alph_protx[:-1] + ['*']:
                    for a2 in rje_seq.alph_protx[:-1] + ['*']: rfdata[rf]['%s%s' % (a1,a2)] = 0; rfobs[rf]['%s%s' % (a1,a2)] = 0; rfexp[rf]['%s%s' % (a1,a2)] = 0
            for x in rje_seq.alph_dna[:-1]: ntfreq[x] = 0
            seqlist = self.obj['SeqList'] 
            ### ~ [2] Count sequence attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,seqlist.seqNum())
            for seq in seqlist.seq:
                self.progLog('\r#ATT','Counting sequence attributes: %.2f%%' % (sx/stot)); sx += 100.0
                for x in seq.info['Sequence']:
                    if x in ntfreq: ntfreq[x] += 1
                rf6 = rje_sequence.sixFrameTranslation(seq.info['Sequence'])
                for r in rf6:
                    rseq = rf6[r]
                    rf = 'RF%d' % r
                    for i in range(len(rseq)):
                        a = rseq[i]; dia = rseq[i:i+2]
                        if a in rfdata[rf]: rfdata[rf][a] += 1
                        if dia in rfdata[rf]: rfdata[rf][dia] += 1
            self.printLog('\r#ATT','Counting sequence attributes complete.')
            ### ~ [3] Calculate Observed & Expected ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ntobs = rje.dictFreq(ntfreq,total=True,newdict=True)
            ntcomp = {'Total':ntobs['Total']}
            for xy in ['AT','GC']: ntcomp[xy[0]] = ntobs[xy[1]]; ntcomp[xy[1]] = ntobs[xy[0]]
            for rf in ['RF1','RF2','RF3','RF-1','RF-2','RF-3']:
                aafreq = {}
                for a in rje_seq.alph_protx[:-1] + ['*']: aafreq[a] = rfdata[rf][a]
                aafreq = rje.dictFreq(aafreq,total=True,newdict=True)
                for a in rje_seq.alph_protx[:-1] + ['*']: rfobs[rf][a] = rfdata[rf][a]; rfexp[rf][a] = 0
                for n1 in 'GATC':
                    for n2 in 'GATC':
                        for n3 in 'GATC':
                            codon = '%s%s%s' % (n1, n2, n3)
                            aa = rje_sequence.dna2prot(codon)
                            if rf[-2] == '-': rfexp[rf][aa] += (int(ntobs['Total']/3.0) * ntcomp[n1] * ntcomp[n2] * ntcomp[n3])
                            else: rfexp[rf][aa] += (int(ntobs['Total']/3.0) * ntobs[n1] * ntobs[n2] * ntobs[n3])
                            #self.deBug('%s: %s x %s x %s x %s' % (aa,(ntobs['Total'] - 2), rfobs[rf][n1], rfobs[rf][n2], rfobs[rf][n3]))
                            #self.deBug('%s: %s' % (aa,rfexp[rf][aa]))
                for a1 in rje_seq.alph_protx[:-1] + ['*']:
                    for a2 in rje_seq.alph_protx[:-1] + ['*']:
                        rfexp[rf]['%s%s' % (a1,a2)] = (aafreq['Total'] - 1) * aafreq[a1] * aafreq[a2]
                        rfobs[rf]['%s%s' % (a1,a2)] = rfdata[rf]['%s%s' % (a1,a2)] 
            ### ~ [4] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rfile = rje.baseFile(seqlist.info['Name']) + '.rf.tdt'
            rje.delimitedFileOutput(self,rfile,rfhead,rje_backup=True)
            for a in rje_seq.alph_protx[:-1] + ['*']:
                data = {'Att':a}
                for rf in ['RF1','RF2','RF3','RF-1','RF-2','RF-3']:
                    data['Obs%s' % rf] = rfobs[rf][a]
                    data['Exp%s' % rf] = '%.2f' % rfexp[rf][a]
                    data[rf] = rje.expectString(rfobs[rf][a] / rfexp[rf][a])
                rje.delimitedFileOutput(self,rfile,rfhead,datadict=data)
            for a1 in rje_seq.alph_protx[:-1] + ['*']:
                for a2 in rje_seq.alph_protx[:-1] + ['*']:
                    a = '%s%s' % (a1,a2)
                    data = {'Att':a}
                    for rf in ['RF1','RF2','RF3','RF-1','RF-2','RF-3']:
                        data['Obs%s' % rf] = rfobs[rf][a]
                        data['Exp%s' % rf] = '%.2f' % rfexp[rf][a]
                        data[rf] = rje.expectString(rfobs[rf][a] / rfexp[rf][a])
                    rje.delimitedFileOutput(self,rfile,rfhead,datadict=data)
            self.printLog('#TDT','TDT output complete.')
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: RF Class                                                                                         #
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
    try: RF(mainlog,cmd_list).rfAtt()

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
