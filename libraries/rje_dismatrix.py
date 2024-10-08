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
Module:       rje_dismatrix
Description:  Distance Matrix Module 
Version:      1.0
Last Edit:    05/12/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    DisMatrix Class. Stores distance matrix data.

Commandline:
    outmatrix=X : Type for output matrix - text / mysql / phylip

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
# [ ] : List here
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_DISMATRIX'
        version = '1.0'
        last_edit = 'December 05'  
        description = 'Distance Matrix Module'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print( 'Problem making Info object.')
        raise
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
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
    except: rje.printf('Problem during initial setup.'); raise
#############################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CLASSES                                                                                                 #
#########################################################################################################################

#########################################################################################################################
### DisMatrix Class: 
#########################################################################################################################
class DisMatrix(rje.RJE_Object):     
    '''
    Sequence Distance Matrix Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of matrix
    - Type = Identifying type (e.g. 'PWAln ID'). Should be same as key in SeqList.obj
    - Description = Extra information if desired
    - OutMatrix = Type for output matrix - text / mysql / phylip
    
    Opt:boolean
    - Symmetric = Whether Seq1->Seq2 = Seq2->Seq1

    Stat:numeric

    Obj:RJE_Objects

    Other:
    - matrix = dictionary of Sequence Object pairs and their distance {(Seq1,Seq2):Dis}
    '''
    ### Attributes
    matrix = {}
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','Type','Description','OutMatrix']
        - Stats:float []
        - Opt:boolean ['Symmetric']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','Description','OutMatrix']
        self.statlist = []
        self.optlist = ['Symmetric']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        ### <c> ### Other Attributes
        matrix = {}
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
                self._cmdRead(cmd,type='info',att='OutMatrix')  
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        return
#########################################################################################################################
    def addDis(self,seq1,seq2,dis):     ### Adds a distance to matrix
        '''
        Adds a distance to matrix.
        >> seq1 and seq2:Sequence Objects
        >> dis:Float = 'Distance' measure
        '''
        self.matrix[(seq1,seq2)] = dis
#########################################################################################################################
    def getDis(self,seq1,seq2):     ### Returns distance from matrix or None if comparison not made.
        '''
        Returns distance from matrix or None if comparison not made.
        >> seq1 and seq2:Sequence Objects
        << Float = 'Distance' measure
        '''
        if (seq1,seq2) in self.matrix:
            return self.matrix[(seq1,seq2)]
        elif self.opt['Symmetric'] and (seq2,seq1) in self.matrix:
            return self.matrix[(seq2,seq1)]
        else:
            return None
#########################################################################################################################
    def saveMatrix(self,sequences,filename='matrix.txt',delimit=',',format=None,log=True):   ### Saves matrix in file
        '''
        Saves matrix in file. Uses None for missing values.
        >> sequences:list of Sequence Objects that form keys of matrix in output order
        >> filename:str = output file name
        >> delimit:str = separator between columns
        >> format:str = 'text' [Default], 'mysql' = lower case header, 'phylip' = phylip
        >> log:Boolean = whether to print report to log [True]
        '''
        if not format:
            format = self.info['OutMatrix']
        if format == 'None' and self.opt['MySQL']:
            format = 'mysql'
        elif format == 'None':
            format = 'text'
            
        MATRIX = open(filename, 'w')
        ### <a> ### Header
        header = ['SEQ']
        for seq in sequences:
            try:
                header.append(seq.shortName())
            except:
                header.append(seq)
        if format == 'phylip':
            MATRIX.write('%d\n' % len(sequences))
            delimit = ' '
        else:
            MATRIX.write('%s\n' % rje.join(header,delimit))
        ### <b> ### Matrix
        for seq1 in sequences:
            if format == 'phylip':
                try:
                    if len(seq1.shortName()) < 10:
                        sname = seq1.shortName()
                    else:
                        sname = '%d' % (sequences.index(seq1)+1)
                except:
                    sname = seq1[:10]
                while len(sname) < 10:
                    sname = '%s ' % sname
                matline = [sname]
            else:
                matline = [seq1.shortName()]
            for seq2 in sequences:
                matline.append(str(self.getDis(seq1,seq2)))
            MATRIX.write('%s\n' % rje.join(matline,delimit))
        MATRIX.close()
        if log:
            self.log.printLog('#MAT','%s Distance matrix out to %s (%s format).' % (self.info['Name'],filename,format))
#########################################################################################################################
### End of DisMatrix Class
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
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:        
        print( '\n\n *** No standalone functionality! *** \n\n')
        
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
