#!/usr/local/bin/python

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
Module:       rje_scansite
Description:  Converter for Scansite to MySQL
Version:      0.0
Last Edit:    11/01/06
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    Converts multiple scansite output files into a single table for MySQL upload. (Use rje_mysql.py to generate Build
    Statement.)

Commandline:
    filelist=LIST   : List of files for upload. Can use wildcard. [*.scansite.txt]
    outfile=FILE    : Output file name [scansite.tdt]

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import glob
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
        program = 'RJE_SCANSITE'
        version = '0.0'
        last_edit = 'January 06'  
        description = 'Scansite Converter Module'
        author = 'Dr Richard J. Edwards.'
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
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
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
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#############################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
re_scansite = '^(\S+)\s+(\S+)\s+(\S+)\s+(\S)(\d+)\s+(\S+)\s+(\d+\.\d+)\S\s+(\S+)\s+(\S+)'  # 8 fields
#########################################################################################################################
### SECTION II: CLASSES                                                                                                 #
#########################################################################################################################

#########################################################################################################################
### Scansite Class:                                                                                                     #
#########################################################################################################################
class Scansite(rje.RJE_Object):     
    '''
    Scansite Converter Class. Author: Rich Edwards (2006).

    Info:str
    - Name = Name of output file
    
    Opt:boolean

    Stat:numeric

    List:list
    - FileList = List of input files

    Dict:dictionary    

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Opt:boolean []
        - Stats:float []
        - List:list ['FileList']
        - Dict:dictionary []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name']
        self.optlist = []
        self.statlist = []
        self.listlist = ['FileList']
        self.dictlist = []
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.info['Name'] = 'scansite.tdt'
        self.list['FileList'] = glob.glob('*.scansite.txt')
        ### <c> ### Other Attributes
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
                self._cmdRead(cmd,type='info',att='Name',arg='outfile')
                self._cmdRead(cmd,type='glist',att='FileList')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Converter Method                                                                                   #
#########################################################################################################################
    def convert(self,filelist=[],outfile=None):      ### Converts scansite output files in FileList to Outfile
        '''
        Converts scansite output files in FileList to Outfile.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            if len(filelist) < 1:
                filelist = self.list['FileList']
            if not outfile:
                outfile = self.info['Name']          
            if len(filelist) < 1:
                self.log.errorLog('No scansite files to convert! %s unchanged/not made.' % outfile,printerror=False)
                return False
            delimit = rje.getDelimit(self.cmd_list)
            ext = rje.delimitExt(delimit)
            if ext != outfile[-3:]:
                newfile = outfile[:-3] + ext
                if rje.yesNo('Change file name from %s to %s?' % (outfile, newfile)):
                    outfile = newfile
            self.log.printLog('#OUT','Converting %d file(s), output to %s.' % (len(filelist),outfile))

            ### Output File ###
            _stage = 'Output File'
            if not self.opt['Append'] or not os.path.exists(outfile):   # Create with header
                OUTFILE = open(outfile,'w')
                headers = ['seq_id','enzyme','enz_group','aa','pos','score','percentile','matchseq','sa']
                rje.writeDelimit(OUTFILE,headers,delimit)
            else:
                OUTFILE = open(outfile,'a')

            ### Conversion ###
            _stage = 'Conversion'
            sx = 0
            for infile in filelist:
                if not os.path.exists(infile):
                    self.log.errorLog('Input file %s does not exist! :o(' % infile,False,False)
                    continue
                fx = 0
                INFILE = open(infile,'r')
                inline = rje.nextLine(INFILE)
                while inline != None:
                    if rje.matchExp(re_scansite,inline):
                        scanlist = rje.matchExp(re_scansite,inline)
                    rje.writeDelimit(OUTFILE,scanlist,delimit)
                    sx += 1
                    fx += 1
                    rje.progressPrint(self,sx)
                    inline = rje.nextLine(INFILE)
                self.log.printLog('#OUT','%s scansite results from %s. (%s Total.)' % (rje.integerString(fx),infile,rje.integerString(sx)))
                INFILE.close()

            ### End ###
            _stage = 'End'
            OUTFILE.close()
            self.log.printLog('#OUT','%s scansite results output to %s.' % (rje.integerString(sx),outfile))
            return True            
        except:
            self.log.errorLog('Error in convert(%s)' % _stage,printerror=True,quitchoice=False)
            raise   
#########################################################################################################################
### End of Scansite Template                                                                                            #
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
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:        
        scansite = Scansite(mainlog,cmd_list)
        scansite.convert()
        
    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
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
