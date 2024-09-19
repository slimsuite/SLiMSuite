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
Module:       peptide_dismatrix
Description:  Peptide Distance Matrix Generator
Version:      1.1
Last Edit:    27/10/05
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Generates distance matrix for input peptides.

Commandline:
    seqin=FILE  : Peptide Sequence File
    aaprop=FILE : Amino Acid property matrix file. [aaprop.txt]
    outmatrix=X : Type for output matrix - text / mysql / phylip [Phylip]
    delimit=X   : Text delimiter for text or mysql file
    method=X    : Peptide Distance Method to use [ds_prop]
        - ds_prop = Denis Shields distance using AA properties
        - ds_id = Denis Shields distance using AA identity
        - pam = ML PAM distance
        - tot_prop = Summed difference across all properties in all residues
        - best_prop = Minimum aligned property difference by circularising one peptide and sliding versus other

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje, rje_aaprop, rje_dismatrix, rje_pam, rje_seq
Other modules needed: rje_blast, rje_sequence, rje_uniprot
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
import rje_dismatrix
import rje_pam
import rje_seq
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 1.0 - Working version for DS Method
# 1.1 - Altered ds_id and ds_prop to subtract mean self-distance. Add best_prop method.
#########################################################################################################################
### Major Functionality to Add
# [ ] : List here
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'PEPTIDE_DISMATRIX'
        version = '1.1'
        last_edit = 'October 05'  
        description = 'Peptide DisMatrix Module'
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
### PepDis Class:                                                                                                   #
#########################################################################################################################
class PepDis(rje.RJE_Object):     
    '''
    PepDis Class. Author: Rich Edwards (2005).

    Info:str
    - OutMatrix = Type for output matrix - text / mysql / phylip
    - Method = Peptide Distance Method to use [ds_prop]
    
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
        - Info:str ['OutMatrix','Method']
        - Opt:boolean []
        - Stats:float []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['OutMatrix','Method']
        self.optlist = []
        self.statlist = []
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        self.info['OutMatrix'] = 'phylip'
        self.info['Method'] = 'ds_prop'
        
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
                self._cmdRead(cmd,type='info',att='OutMatrix')
                self._cmdRead(cmd,type='info',att='Method')  
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Class Methods                                                                                #
#########################################################################################################################
    def _pepDis(self):      ### Peptide Distance
        '''
        Peptide Distance.
        '''
        try:
            ### <0> ### Setup
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T'])
            dismatrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            dismatrix.info['Name'] = self.info['Method']
            dismatrix.opt['Symmetric'] = True
            if self.info['Method'] in ['ds_prop','tot_prop','best_prop']:
                aaprop = rje_aaprop.AAPropMatrix(self.log,self.cmd_list)
                #aaprop.readAAProp()
                aaprop.makePropDif()
            elif self.info['Method'] == 'pam':
                pam = rje_pam.PamCtrl(log=self.log,cmd_list=self.cmd_list)
            ### <1> ### Make DisMatrix
            for seq1 in seqlist.seq:
                for seq2 in seqlist.seq:
                    if seqlist.seq.index(seq1) > seqlist.seq.index(seq2):   # No need to calculate - symmetrical!
                        continue
                    dis = 0
                    if seq1 == seq2 and self.info['OutMatrix'] == 'phylip':
                        dis = 0
                    elif self.info['Method'] in ['ds_prop','ds_id']:
                        (self_dis1,self_dis2) = (0,0)
                        for r1 in range(seq1.seqLen()):
                            for r2 in range(r1,seq2.seqLen()):
                                (a1,a2) = (seq1.info['Sequence'][r1],seq2.info['Sequence'][r2])
                                (s1,s2) = (seq1.info['Sequence'][r2],seq2.info['Sequence'][r1])
                                phys_dis = r2 - r1
                                if self.info['Method'] == 'ds_prop':
                                    dis += (aaprop.pdif['%s%s' % (a1,a2)] * (seq1.seqLen() - phys_dis))
                                    self_dis1 += (aaprop.pdif['%s%s' % (a1,s1)] * (seq1.seqLen() - phys_dis))
                                    self_dis2 += (aaprop.pdif['%s%s' % (a2,s2)] * (seq1.seqLen() - phys_dis))
                                elif self.info['Method'] == 'ds_id' and a1 != a2:
                                    dis += (seq1.seqLen() - phys_dis)
                                if self.info['Method'] == 'ds_id' and a1 != s1:
                                    self_dis1 += (seq1.seqLen() - phys_dis)
                                if self.info['Method'] == 'ds_id' and a2 != s2:
                                    self_dis2 += (seq1.seqLen() - phys_dis)
                        dis -= (self_dis1 + self_dis2) / 2.0
                    elif self.info['Method'] == 'tot_prop':
                        proptot = {}
                        for property in aaprop.prop.keys():
                            proptot[property] = {seq1:0.0,seq2:0.0}
                        for seq in [seq1,seq2]:
                            for r in range(seq.seqLen()):
                                aa = seq.info['Sequence'][r]
                                for property in aaprop.prop.keys():
                                    proptot[property][seq] += string.atof(aaprop.prop[property][aa])
                        for property in aaprop.prop.keys():
                            if proptot[property][seq1] > proptot[property][seq2]:
                                dis += (proptot[property][seq1] - proptot[property][seq2])
                            else:
                                dis += (proptot[property][seq2] - proptot[property][seq1])
                    elif self.info['Method'] == 'pam':
                        dis = pam.pamML(ancseq=seq1.info['Sequence'],descseq=seq2.info['Sequence'])
                    elif self.info['Method'] == 'best_prop':
                        min_dis = seq1.seqLen() * len(aaprop.prop)
                        pepseq1 = seq1.info['Sequence']
                        for c in range(seq1.seqLen()):  # Circular start
                            dis = 0
                            pepseq2 = seq2.info['Sequence'][c:] + seq2.info['Sequence'][:c]
                            for r in range(seq1.seqLen()):
                                (a1,a2) = (pepseq1[r],pepseq2[r])
                                dis += aaprop.pdif['%s%s' % (a1,a2)]
                            if dis < min_dis:
                                min_dis = dis
                        dis = min_dis
                    dismatrix.addDis(seq1,seq2,dis)
            ### <2> ### Output
            if self.info['OutMatrix'] == 'phylip':
                delimit = ' '
                format = 'phylip'
            else:
                delimit = rje.getDelimit(self.cmd_list,',')
                format = 'None'
            outfile = '%s.%s.%s' % (rje.baseFile(seqlist.info['Name'],True),self.info['Method'],rje.delimitExt(delimit))
            dismatrix.saveMatrix(seqlist.seq,outfile,delimit,format=format)
                    
        except:
            self.log.errorLog('Error in _pepDis',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of PepDis Template                                                                                            #
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
        pepdis = PepDis(mainlog,cmd_list)
        pepdis._pepDis()
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
########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
