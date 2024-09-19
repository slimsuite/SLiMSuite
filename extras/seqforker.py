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
Program:      SeqForker
Description:  Generic Sequence Analysis Forking Script
Version:      1.0
Last Edit:    23/09/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to take a large input sequence dataset, split it into chunks, fork out the chunks to some
    process and then stick all the results back together at the end. It is designed to be flexible and work for any
    analyses where each sequence is looked at independently and filenames can be given to the program performing the
    analysis.

Commandline:
    seqin=FILE          : Input sequence file [None]
    split=X             : Number of sequences per split file [0]
    startfrom=X         : Will pick up program at this sequence, where X is the name or accession number [None]
                          * Remember to set append=T if picking up a crashed run *
    append=T/F          : Will append to results files rather than overwrite [False]
    forkprog=PATH       : Common program system call for all forked splits of sequence file (e.g. program) [None]
    forkcmd="blah blah" : Common system commands for all forked splits of sequence file (e.g. program options) [None]
    outcmd=X            : Command line option for giving output file name. Will be altered to match forked splits (*.*) [None]
    seqincmd=X          : Command given to program for input file name [seqin=]
                          * SeqForker will stitch together "forkprog seqincmd=FILE forkcmd outcmd" *
    
General Commandline:
    v=X         : Sets verbosity (-1 for silent) [0]
    i=X         : Sets interactivity (-1 for full auto) [0]
    log=FILE    : Redirect log to FILE [Default = calling_program.log]
    newlog=T/F  : Create new log file. [Default = False: append log file]

Forking Commandline:
    noforks=T/F : Whether to avoid forks [False]
    forks=X     : Number of parallel sequences to process at once [0]
    killforks=X : Number of seconds of no activity before killing all remaining forks. [3600]

Uses general modules: copy, glob, os, re, string, sys, time
Uses RJE modules: rje, rje_seq
Other modules needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
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
import rje_seq
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 1.0 - Intital Complete Working Version
#########################################################################################################################
### Major Functionality to Add
# [ ] Move splitting of file to within forking routine? Should speed up and reduce disk usage.
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'SeqForker'
        version = "1.0"
        last_edit = "September 05"  
        description = "Generic Sequence Forking Module"
        author = "Dr Richard J. Edwards."
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
### SeqForker Class                                                                                                     #
#########################################################################################################################
class SeqForker(rje.RJE_Object):     
    '''
    SeqForker Class. Author: Rich Edwards (2005).

    Info:str
    - StartFrom = Will pick up program at this sequence, where X is the name or accession number [None]
    - ForkProg = Common program system call for all forked splits of sequence file (e.g. program) [None]
    - ForkCmd = Common system commands for all forked splits of sequence file (e.g. program options) [None]
    - SeqInCmd = Command given to program for input file name [seqin=]
    - OutCmd = Command line option for giving output file name. Will be altered to match forked splits [None]

    Opt:boolean

    Stat:numeric
    - Split = Number of sequences per split file [0]

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['StartFrom','ForkProg','SeqInCmd','ForkCmd','OutCmd']
        - Opt:boolean []
        - Stats:float []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['StartFrom','ForkProg','SeqInCmd','ForkCmd','OutCmd']
        self.optlist = []
        self.statlist = ['Split']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        self.info['SeqInCmd'] = 'seqin='
        ### <c> ### Other Attributes
        self._setForkAttributes()   # Delete if no forking
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
                self._forkCmd(cmd)  # Delete if no forking
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='StartFrom')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='info',att='ForkProg')
                self._cmdRead(cmd,type='info',att='SeqInCmd')
                self._cmdRead(cmd,type='info',att='ForkCmd')
                self._cmdRead(cmd,type='info',att='OutCmd')
                self._cmdRead(cmd,type='stat',att='Split')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Forker Method                                                                  #
#########################################################################################################################
    def forker(self):      ### Generic forking method
        '''
        Generic method for forking (without threads).
        Add description here (and arguments.)
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Fork Setup'
            forkx = int(self.stat['Forks'])   # Number of forks to have running at one time
            if self.opt['Win32'] or forkx < 1:
                self.opt['NoForks'] = True
            forks = []  # List of active fork PIDs
            killforks = int(self.stat['KillForks']) # Time in seconds to wait after main thread has apparently finished
            forking_condition = True    # Condition to keep forking

            ### Sequence List setup ###
            _stage = '<1> Forking'
            seqx = 0    # Sequence Counter
            subx = 0    # Subset sequence counter
            outfile = None  # Output file name
            randlist = []   # List of random strings for split sequence files
            filedict = {}   # Dictionary of input files for each random string
            seqlist = rje_seq.SeqList(log=self.log,cmd_list=['autoload=F']+self.cmd_list)
            seqlist.makeBaseFile()
            SEQFILE = open(seqlist.info['Name'], 'r')
            (seq,lastline) = seqlist.nextFasSeq(SEQFILE,'Starting')
            while seq:
                seqlist.seq = [seq]
                if self.info['StartFrom'] != 'None':    # Not yet reached wanted sequence
                    if self.info['StartFrom'] in [seq.info['Name'], seq.info['ID'], seq.info['AccNum'], seq.shortName()]:
                        self.info['StartFrom'] = 'None'
                if self.info['StartFrom'] == 'None':    # Wanted sequence
                    if outfile:   # Create new file
                        SEQOUT = open(outfile,'a')
                    else:
                        rs = rje.randomString(6)
                        while rs in randlist:
                            rs = rje.randomString(6)
                        outfile = '%s.%s.fas' % (seqlist.info['Basefile'],rs)
                        SEQOUT = open(outfile,'w')
                        randlist.append(rs)
                        filedict[rs] = outfile
                    SEQOUT.write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
                    SEQOUT.close()
                    seqx += 1
                    subx += 1
                    if subx == self.stat['Split']:  # Finished split
                        self.log.printLog('#SEQ','%s sequences output to %s.' % (rje.integerString(subx),outfile))
                        outfile = None
                        subx = 0
                (seq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
            if subx > 0:
                self.log.printLog('#SEQ','%s sequences output to %s.' % (rje.integerString(subx),outfile))
            self.log.printLog('#SEQ','%s sequences output in total to %d files.' % (rje.integerString(seqx),len(randlist)))
            # Now have the list of random strings in randlist (in order) and filenames in filedict

            ### <1> ### Forking                            
            killtime = time.time()
            dealt_with = 0      # Split files dealt with
            while dealt_with < len(randlist) or len(forks):

                ## <a> ## forks
                _stage = '<1a> New Forks'
                while dealt_with < len(randlist) and (len(forks) < forkx or self.opt['NoForks']):     # Add more forks
                    _stage = '<1a-i> Fork: Get stuff for fork'
                    killtime = time.time()  # Reset killtime - still doing stuff
                    # Add new fork
                    _stage = '<1a-ii> Fork: New Fork'
                    new_fork_id = randlist[dealt_with]
                    dealt_with += 1
                    
                    outcmd = rje.split(self.info['OutCmd'],'.')
                    if len(outcmd) > 1:
                        outcmd = outcmd[:-1] + [new_fork_id] + outcmd[-1:]
                    else:
                        outcmd = outcmd + [new_fork_id] + ['resfile']
                    outcmd = rje.join(outcmd,'.')
                    forkcmd = '%s %s%s %s %s log=%s.log newlog=T i=-1' % (self.info['ForkProg'],self.info['SeqInCmd'],filedict[new_fork_id],outcmd,self.info['ForkCmd'],new_fork_id)
                    if self.opt['NoForks']:
                        os.system(forkcmd)                        
                    else:   # Forks
                        newpid = os.fork() 
                        if newpid == 0: # child 
                            os.system(forkcmd)
                            sys.exit()    # Exit process 
                        elif newpid == -1: # error
                            self.log.errorLog('Problem forking %s.' % new_fork_id)
                        else: 
                            forks.append(newpid)    # Add fork to list 
            
                ## <b> ## Monitor and remove finished forks
                _stage = '<1b> Finished Forks'
                forklist = self._activeForks(forks)
                if len(forklist) != len(forks):
                    self.verbose(0,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    forks = forklist[0:]
                        
                self.verbose(3,3,'End of a Cycle.',2)

                ## <c> ## Look for eternal hanging of forks
                _stage = '<1c> Hanging'
                if time.time() - killtime > killforks:
                    self.verbose(0,1,'\n%d seconds of main program inactivity. %d forks still active!' % (killforks,len(forks)),1)
                    for fork in forks:
                        self.verbose(0,2,' => Fork PID %d still Active!' % (fork),1)
                    if rje.yesNo('Kill?'):
                        break   #!# killing options
                    else:
                        killtime = time.time()

            ### <3> ### Finish
            _stage = '<3> Finish'
            if len(forks) > 0:
                self.log.errorLog('%d Forks still active after %d seconds of main program inactivity' % (len(forks),killforks),True)
            else:
                self.verbose(0,1,'Forks have finished.',2)

            ### <4> ### Recompile results
            for randstr in randlist:
                os.unlink(filedict[randstr])
                rje.fileTransfer(fromfile='%s.log' % randstr,tofile=self.log.info['Name'],deletefrom=True)            
                outfiles = glob.glob('*.%s.*' % randstr)
                for outfile in outfiles:
                    compfile = outfile.split('.')
                    compfile.remove(randstr)
                    compfile = rje.join(compfile,'.')
                    if randstr == randlist[0] and os.path.exists(compfile) and not self.opt['Append']:
                        os.unlink(compfile)
                    rje.fileTransfer(fromfile=outfile,tofile=compfile,deletefrom=True)
                    self.verbose(1,2,'Copying results data from %s to %s...' % (outfile,compfile),0)
                self.verbose(0,1,'%d results files copied for Split %d.' % (len(outfiles),(randlist.index(randstr)+1)),1)
            self.log.printLog('#OUT','Results for %d splits compiled.' % len(randlist))
            
        except SystemExit:  # Don't want forks raising an Exception upon exiting
            sys.exit()
        except:
            self.log.errorLog('Error in forker(%s):' % _stage,printerror=True,quitchoice=False)     
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _method(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
        except:
            self.log.errorLog('Error in _method(%s)' % _stage,printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SeqForker                                                                                                    #
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
        seqforker = SeqForker(mainlog,cmd_list)
        seqforker.forker()        

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
