#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_sleeper
Description:  RJE Sleeper module to keep Putty active on Bioinformatics/Cerberus
Version:      0.0
Last Edit:    10/10/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Occasionally prints to screen to keep putty active.

Commandline:
    sleep=X : Seconds to sleep for between prints [600 (10 mins)]
    wake=X  : Total seconds until finishing [864,000 (10 days)]
    tofile=T/F  : Whether to print to file [True]
    toscreen=T/F  : Whether to print to file [False]

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import os
#import random
import re
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
#########################################################################################################################
### Major Functionality to Add
# [ ] List here
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_SLEEPER'
        version = "0.0"
        last_edit = "October 05"  
        description = "RJE Sleeper Module"
        author = "Dr Richard J. Edwards."
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
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
            print "\n\nHelp for " + info.program + " V:" + info.version + ": " + time.asctime(time.localtime(info.start_time)) + "\n"
            out.verbose(-1,-1,text=__doc__)
            sys.exit()
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
#############################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    ### <0> ### Objects setup
    try:
        ## <a> ## Initial Command Setup & Info
        cmd_list = sys.argv[1:]
        info = makeInfo()
        cmd_list = rje.getCmdList(cmd_list,info=info)      ### Load defaults from program.ini
        ## <b> ## Out object
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except:
        print "Problem during initial setup."
        raise
#############################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CLASSES                                                                                                 #
#########################################################################################################################

#########################################################################################################################
### Sleeper Class:                                                                                                      #
#########################################################################################################################
class Sleeper(rje.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2005).

    Info:str
    
    Opt:boolean
    - ToFile = Whether to print to file [True]
    - ToScreen = Whether to print to file [False]

    Stat:numeric
    - Sleep = Seconds to sleep between prints [600]
    - Wake = Total seconds until finishing [864,000 (10 days)]

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str []
        - Opt:boolean ['ToFile','ToScreen']
        - Stats:float ['Sleep','Wake']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = []
        self.optlist = ['ToFile','ToScreen']
        self.statlist = ['Sleep','Wake']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        ### <c> ### Other Attributes
        self.stat['Sleep'] = 600
        self.stat['Wake'] = 864000
        self.stat['Interactive'] = -1
        self.opt['ToFile'] = True
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='int',att='Sleep')  
                self._cmdRead(cmd,type='int',att='Wake')  
                self._cmdRead(cmd,type='opt',att='ToFile')  
                self._cmdRead(cmd,type='opt',att='ToScreen')  
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Method                                                                                             #
#########################################################################################################################
    def _sleeper(self):      ### Sleeper method
        '''
        Sleeps, prints, sleeps etc.
        '''
        try:
            total_sleep = 0
            while total_sleep < self.stat['Wake']:
                time.sleep(self.stat['Sleep'])
                self.log.printLog('#TIME','rje_sleep keeping window active!',screen=self.opt['ToScreen'],log=self.opt['ToFile'])
                total_sleep += self.stat['Sleep']
        except:
            self.log.errorLog('Error in _sleeper',printerror=True,quitchoice=False)
            raise   
#########################################################################################################################
### End of Sleeper Template                                                                                            #
#########################################################################################################################
       
#########################################################################################################################
## End of SECTION II                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    try:
        ### <0>  ### Basic Setup of Program
        [info,out,mainlog,cmd_list] = setupProgram()
        
        ### <1> ### Rest...
        sleeper = Sleeper(mainlog,cmd_list+['i=-1'])
        sleeper._sleeper()

        ### <X> ### End
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog("User terminated.\n")
    except:
        print "Unexpected error:", sys.exc_info()[0]
    mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
