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
Module:       peptide_stats
Description:  HRB Peptide Statistics Generator
Version:      1.0
Last Edit:    09/11/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Generates property outputs for HRB peptides.

Commandline:
    seqin=FILE  : Peptide Sequence File
    aaprop=FILE : Amino Acid property matrix file. [aaprop.txt]
    delimit=X   : Text delimiter for text or mysql file

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje, rje_aaprop, rje_pam, rje_seq
Other modules needed: rje_blast, rje_dismatrix, rje_sequence, rje_uniprot
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
import rje_aaprop
#import rje_pam
import rje_seq
#########################################################################################################################
### History
# 1.0 - Working version based on peptide_dismatrix
#########################################################################################################################
### Major Functionality to Add
# [ ] : List here
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'PEPTIDE_STATS'
        version = '1.0'
        last_edit = 'November 05'  
        description = 'HRB Peptide Statistics Generator'
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
### PepStats Class:                                                                                                     #
#########################################################################################################################
class PepStats(rje.RJE_Object):     
    '''
    PepStats Class. Author: Rich Edwards (2005).

    Info:str
    
    Opt:boolean

    Stat:numeric

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
        - Opt:boolean []
        - Stats:float []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        
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
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Class Methods                                                                                           #
#########################################################################################################################
    def _pepStats(self):      ### Peptide Distance
        '''
        Peptide Distance.
        '''
        try:
            ### Setup ###
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T'])
            aaprop = rje_aaprop.AAPropMatrix(self.log,self.cmd_list)
            aaprop.makePropDif()
            delimit = rje.getDelimit(self.cmd_list)

            ### Output File Setup ###
            OUTFILE = open('hrb.pepstats.%s' % rje.delimitExt(delimit),'w')
            headlist = ['peptide']
            ## 10 Dimensional Peptide Property Output ##
            for property in rje.sortKeys(aaprop.prop):
                headlist.append(property.lower())
                for aa in aaprop.prop[property].keys():
                    try:
                        if aa not in ['-','X']:
                            aaprop.prop[property][aa] = string.atoi(aaprop.prop[property][aa])
                    except:
                        print aaprop.prop, property, aa, aaprop.prop[property][aa]
                        raise
            ## Additional Stats ##
            headlist.append('net_charge')
            #headlist.append('hydrophobicity')
            headlist.append('charge_balance')
            headlist.append('hydrophobic_balance')
            #headlist.append('hydrophobicity_balance')
            ## Output
            rje.writeDelimit(OUTFILE,headlist,delimit)
            
            ### Calculate stats ###
            for pep in seqlist.seq:
                pepname = pep.shortName()
                if rje.matchExp('^(\S+_\d[CQ])',pepname):
                    pepname = rje.matchExp('^(\S+_\d[CQ])',pepname)[0]
                outlist = [pepname]
                pepseq = pep.info['Sequence']
                ## 10 Dimensional Peptide Property Output ##
                for property in rje.sortKeys(aaprop.prop):
                    px = 0
                    for aa in pepseq:
                        px += aaprop.prop[property][aa]
                    outlist.append('%d' % px)
                ## Additional Stats ##
                net_charge = 0
                for aa in pepseq:
                    net_charge += (aaprop.prop['Positive'][aa] - aaprop.prop['Negative'][aa])
                outlist.append('%d' % net_charge)
                charge_balance = 0
                hydrophobic_balance = 0
                for r in range(len(pepseq)):
                    charge_balance += aaprop.prop['Charged'][pepseq[r]] * (1.0 / (r+1))
                    charge_balance -= aaprop.prop['Charged'][pepseq[r]] * (1.0 / (10-r))
                    hydrophobic_balance += aaprop.prop['Hydrophobic'][pepseq[r]] * (1.0 / (r+1))
                    hydrophobic_balance -= aaprop.prop['Hydrophobic'][pepseq[r]] * (1.0 / (10-r))
                outlist.append('%.3f' % charge_balance)
                outlist.append('%.3f' % hydrophobic_balance)
                rje.writeDelimit(OUTFILE,outlist,delimit)
                        
            ### Finish ###
            OUTFILE.close()
                    
        except:
            self.log.errorLog('Error in _pepStats',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of PepStats Class                                                                                               #
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
        pepstats = PepStats(mainlog,cmd_list)
        pepstats._pepStats()

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
