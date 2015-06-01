#!/usr/local/bin/python

"""
Module:       rem_parser
Description:  Module to Parse removed sequence links from Log Files
Version:      0.0
Last Edit:    26/04/05

Function:
    To be added.

Commandline:
    remlog=FILE : File containing removal data

Uses general modules: os, sys, threading, time
Uses RJE modules: rje
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
#import copy
import os
#import random
#import re
import sys
import threading    # Threading & Forking bit only
import time
#############################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#############################################################################################################################
### History
# 0.0 - Initial Compilation.
#############################################################################################################################
### Major Functionality to Add
# [ ] List here
#############################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'REM_PARSER'
        version = "0.0"
        last_edit = "Apr 05"  
        description = "Removed Sequence Log parser"
        author = "Dr Richard J. Edwards."
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#############################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out([])
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print "\n\nHelp for " + info.program + " V:" + info.version + ": " + time.asctime(time.localtime(info.start_time)) + "\n"
            out.verbose(out.verbosity,out.interactive,__doc__,1)
            out.verbose(out.verbosity,out.interactive,rje.__doc__,1)
            # ! # Add additional user modules and objects here!
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.interactive > 0:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except KeyboardInterrupt:
        raise
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
        out = rje.Out(cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except:
        print "Problem during initial setup."
        raise
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#############################################################################################################################

#############################################################################################################################
### RemParser: 
#############################################################################################################################
class RemParser(rje.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of File
    
    Opt:boolean

    Stat:numeric

    Obj:RJE_Objects
    '''
    ### Attributes
    acc_links = {}  # Dictionary of removed accession numbers and unltimate links to others.
    rejects = {}    # Dictionary of dictionaries of rejects
#############################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#############################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Stats:float []
        - Opt:boolean []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name']
        for info in self.infolist:
            self.info[info] = 'None'
        self.statlist = []
        for stat in self.statlist:
            self.stat[stat] = 0.0
        self.optlist = []
        for opt in self.optlist:
            self.opt[opt] = False
        self.objlist = []
        for obj in self.objlist:
            self.obj[obj] = None
        ### <b> ### Defaults
        ### <c> ### Other Attributes
        self.acc_links = {}
        self.rejects = {}
#############################################################################################################################
#    def _cmdRead(self,type='info',att=None,arg=None,cmd=None):     ### Sets self.type[att] from commandline command cmd
#############################################################################################################################
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
                self._cmdRead(type='info',att='Name',arg='remlog',cmd=cmd)
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#############################################################################################################################
    def run(self):  ### Loads and parses file
        '''Loads and parses file.'''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            self.log.printLog('#INF','Reading data from %s.' % self.info['Name'])
            REMLOG = open('%s' % self.info['Name'],'r')
            self.acc_links = {} # Links to other sequences
            self.rejects = {}   # Dictionary of dictionary of rejects: [kept,reason,decision]
            tmplink = {}   # Dictionary of temporary links
            
            ### <1> ### Build tmplinks and remreason from file
            while 1:
                _stage = '<1a> Readline'
                line = REMLOG.readline()
                if line in ('', None):
                    break
                _stage = '<1b> Process'
                details = rje.matchExp('Deleted\s+(\S+):\s+(\S.+)$',line)
                #REM    01:15:01        nvl_HUMAN/ENSP00000299720
                if details == [] or details == None:
                    details = rje.matchExp('^#REM\s+\S+\s+(\S+)\s+(\S.+)$',line)
                print details
                if details:
                    name = rje.matchExp('\S+\/(\S+)$', details[0])
                    #print name
                    if name:
                        acc = name[0]
                    else:
                        acc = details[0]
                    self.rejects[acc] = {}
                    _stage = '<1b-i> Extract details'
                    # kept = acc_num of kept sequence 
                    # reason = Redundancy / Duplication / Too short
                    # decision = SwissProt over TrEMBL / Arbitrary etc.
                    if details[1].find('Too short') >= 0:   # Short Sequence Filter
                        self.rejects[acc]['Reason'] = 'Too Short'
                        self.rejects[acc]['Kept'] = 'None'
                        self.rejects[acc]['Decision'] = 'N/A'
                    elif details[1].find('no BLAST hits') >= 0:   # Probable short sequence
                        self.rejects[acc]['Reason'] = 'No (self) BLAST.'
                        self.rejects[acc]['Kept'] = 'None'
                        self.rejects[acc]['Decision'] = 'N/A'
                    elif details[1].find(' ID ') >= 0:   # %ID Redundancy Filter Probable short sequence
                        remdetail = rje.matchExp('^(\S.+\S)\s+ID vs (\S+): (\S.+)$', details[1])
                        if remdetail:
                            self.rejects[acc]['Reason'] = rje.matchExp('^(\S+)',remdetail[0])[0]
                            name = rje.matchExp('\S+\/(\S+)$', remdetail[1])
                            if name:
                                self.rejects[acc]['Kept'] = name[0]
                            else:
                                self.rejects[acc]['Kept'] = remdetail[1]
                            self.rejects[acc]['Decision'] = remdetail[2]
                        else:
                            print 'ID: Problem matching line:\n%s\n' % line
                    elif details[1].find('duplication') >= 0:   # Duplication/Redundancy from GOPHER
                        remdetail = rje.matchExp('vs (\S+): (\S.+)$', details[1])
                        #print remdetail
                        if remdetail:
                            self.rejects[acc]['Reason'] = 'In-paralogue/Redundancy'
                            name = rje.matchExp('\S+\/(\S+)$', remdetail[0])
                            if name:
                                self.rejects[acc]['Kept'] = name[0]
                            else:
                                self.rejects[acc]['Kept'] = remdetail[0]
                            self.rejects[acc]['Decision'] = remdetail[1]
                        else:
                            print 'duplication: Problem matching line:\n%s\n' % line
                    elif details[1].find('edundancy') >= 0 and details[1].find('with') >= 0:   # Duplication/Redundancy from GOPHER
                        remdetail = rje.matchExp('with (\S+)[:\.]', details[1])
                        if remdetail:
                            self.rejects[acc]['Reason'] = 'In-paralogue/Redundancy'
                            name = rje.matchExp('\S+\/(\S+)$', remdetail[0])
                            if name:
                                self.rejects[acc]['Kept'] = name[0]
                            else:
                                self.rejects[acc]['Kept'] = remdetail[0]
                            self.rejects[acc]['Decision'] = 'Better HUMAN BLAST hit.'
                        else:
                            print 'Redundancy (with): Problem matching line:\n%s\n' % line
                    elif details[1].find('edundancy') >= 0 and details[1].find('vs') >= 0:   # Duplication/Redundancy from GOPHER
                        remdetail = rje.matchExp('vs (\S+)\: (\S.+)$', details[1])
                        if remdetail:
                            self.rejects[acc]['Reason'] = 'In-paralogue/Redundancy'
                            name = rje.matchExp('\S+\/(\S+)$', remdetail[0])
                            if name:
                                self.rejects[acc]['Kept'] = name[0]
                            else:
                                self.rejects[acc]['Kept'] = remdetail[0]
                            self.rejects[acc]['Decision'] = remdetail[1]
                        else:
                            print 'Redundancy (vs): Problem matching line:\n%s\n' % line
                    else:
                        print 'Deleted: Problem matching line:\n%s\n(%s)\n' % (line,details)
                    tmplink[acc] = self.rejects[acc]['Kept']
                else:
                    print 'Absolute Problem matching line to /Deleted\s+(\S+):\s+(\S.+)$/:\n%s\n' % line
            print '\n\n*** Reading finished ***\n\n'
                
            ### <2> ### Build links from tmplinks
            _stage = '<2> Building links'
            REMLOG.close()
            ACCLINKS = open('acc_links.tdt','a')
            ACCLINKS.write('rej_acc\tacc_num\n')
            for acc in tmplink.keys()[0:]:
                print '...%d' % len(tmplink),
                link = tmplink.pop(acc)
                links = [acc,link]
                print '-> %s' %link,
                while link in tmplink.keys() and link not in links:  # Link itself removed - change and recycle
                    link = tmplink[link]
                    links.append(link)
                    print '-> %s' %link,
                if link in self.acc_links.keys():
                    self.acc_links[acc] = self.acc_links[link]
                else:
                    self.acc_links[acc] = link
                ACCLINKS.write('%s\t%s\n' % (acc, self.acc_links[acc]))
            ACCLINKS.close()

            ### <3> ### Output rejects
            _stage = '<3> Output Rejects'
            REJECTS = open('rejects.tdt','a')
            REJECTS.write('rej_acc\tkept_acc\treason\tdecision\n')
            for acc in self.rejects.keys():
                REJECTS.write('%s\t%s\t%s\t%s\n' % (acc,self.rejects[acc]['Kept'],self.rejects[acc]['Reason'],self.rejects[acc]['Decision']))
            REJECTS.close()
            
        except:
            self.log.errorLog('Cataclysmic Error with RemParser.run() stage %s.' % _stage)
#############################################################################################################################
### End of RemParser
#############################################################################################################################
#############################################################################################################################
## End of SECTION II
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                           #
#############################################################################################################################

#############################################################################################################################
### END OF SECTION III
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                               #
#############################################################################################################################
def runMain():
    try:
        ### <0>  ### Basic Setup of Program
        [info,out,mainlog,cmd_list] = setupProgram()
        
        ### <1> ### Rest...
        remparser = RemParser(mainlog,cmd_list)
        #print remparser
        remparser.run()

        ### <X> ### End
    except KeyboardInterrupt:
        mainlog.errorLog("User terminated.\n")
    except:
        print "Unexpected error:", sys.exc_info()[0]
    mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#############################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#############################################################################################################################
### END OF SECTION IV
#############################################################################################################################
