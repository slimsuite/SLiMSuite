#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_motifocc
Description:  Motif Occurrence Module
Version:      0.0
Last Edit:    29/01/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the MotifOcc class. This class if for storing methods and attributes pertinent to an individual
    occurrence of a motif, i.e. one Motif instance in one sequence at one position. This class is loosely based on (and
    should replace) the old PRESTO PrestoHit object. (And, to some extent, the PrestoSeqHit object.) This class is
    designed to be flexible for use with PRESTO, SLiMPickings and CompariMotif, among others.

    In addition to storing the standard info and stat dictionaries, this object will store a "Data" dictionary, which
    contains the (program-specific) data to be output for a given motif. All data will be in string format. The
    getData() and getStat() methods will automatically convert from string to numerics as needed.

Commandline:
    This module has no standalone functionality and should not be called from the commandline.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_seq, rje_sequence
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
#########################################################################################################################
### Major Functionality to Add
# [ ] : List here
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit) = ('RJE_MOTIFOCC', '0.0', 'January 2007')
    description = 'Motif Occurrence Module'
    author = 'Dr Richard J. Edwards.'
    return rje.Info(program,version,last_edit,description,author,time.time())
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info:
            info = makeInfo()
        if not out:
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
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: MotifOcc Class                                                                                          #
#########################################################################################################################
class MotifOcc(rje.RJE_Object):     
    '''
    Motif Occurrence Class. Author: Rich Edwards (2007).

    Note that, unlike other RJE_Object classes, this list of attributes will not be created for all runs with default
    values. Instead, defaults are return for missing values using the getStat() and getData() methods. 

    The list of attributes here is for reference when using other programs. Please try to be consistent across calling
    programs (e.g. PRESTO and SLiM Pickings) for management issues but this is not enforced anywhere. (Some conservation
    and disorder etc. functions will, however, need and return set attributes.

    Note also that some of the original attributes of the PrestoHit object are no longer stored as they can be
    retrieved from the Motif and Sequence object associated with the MotifOcc.

    Info:str
    - Name   : The name of the regular expression given in the input file.
    - SearchVar : The actual RegExp making the match
    - Variant: The sequence of the regular expression variant that is closest to the match.
    - Match  : The sequence of the matched region in the target sequence.
    - Expect : A crude measure of how often this MATCH would be expected by chance in this database. Formatted.
    - PepSeq : Peptide sequence if presto.opt['Peptides'] = True
    
    Stat:numeric
    - ID    : The percentage identity between VARIANT and MATCH. (Always 1 unless mismatches allowed!)
    - Cons  : Conservation of motif across alignment (usealn=T).
    - Hom   : Number of homologues used in conservation score.
    - GlobID: Mean Global %ID between Query and Homologues
    - LocID : Mean Local %ID between Query and Homologues across motif
    - Pos   : The start position of the MATCH in the HIT. (1-n)
    - SA    : The mean Surface Accessibility score (Janin et al) for the matched region
    - Hyd   : The mean Eisenberg hydrophobicity score for the matched region
    *In addition, when MS-MS mode is used, there are the following additional stats:
    - SNT      : The minimum number of single nucleotide substitutions needed to change a stretch of DNA encoding the
                 MATCH into one encoding the VARIANT. (Always 0 unless mismatches allowed! This is most useful when
                 searching a database that may contain sequencing errors or, for MSMS peptides, proteins from a different
                 species.)
    - ORFMwt   : The predicted molecular weight of HIT.
    - M-ORFMWt : The predicted molecular weight of HIT starting with the first Methionine.
    - FragMWt  : The predicted molecular weight of the tryptic fragment of HIT containing MATCH.
    - N[KR]    : The number of amino acids N-terminal of MATCH before reaching a K or R. (A value of 1 indicates that
                 MATCH includes the start of a tryptic fragment.)
    - C[KR]    : The number of amino acids C-terminal of MATCH before reaching a K or R. (A value of 0 indicates that
                 MATCH includes the end of a tryptic fragment.)
    - Rating   : A crude score to help rank hits. The higher the better. [Details to be added later.]

    Dict:
    - Data = String values of output data (can also be in Stat or Info)
    
    Obj:RJE_Objects
    - Seq = rje_sequence.Sequence object
    - Motif = rje_motif.Motif object
    '''
    def hit(self):
        if self.obj['Seq']:
            return self.obj['Seq'].shortName()
        return None
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object. Unlike other RJE classes, most of these are optional.'''
        ### Basics ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = ['Data']
        self.objlist = ['Seq','Motif']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ### 
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <3> ### Methods for populating attributes                                                                       #
#########################################################################################################################
    def _addSeq(self,seqname,fasdb):    ### Fishes seqname from fasdb (using fastacmd) and sets as self.obj['Seq']
        '''Fishes seqname from fasdb (using fastacmd) and sets as self.obj['Seq'].'''
        try:
            ### Setup ###
            scmd = self.cmd_list[0:] + ['seqin=None','autoload=F']
            self.obj['Seq'] = rje_seq.SeqList(self.log,scmd).seqFromFastaCmd(seqname,fasdb)
            if not self.obj['Seq']:
                self.log.errorLog('Fastacmd failure (%s from %s). Check win32=T/F!' % (seqname,fasdb),printerror=False)
        except:
            self.log.errorLog('Error in MotifOcc._addSeq()',printerror=True,quitchoice=True)
#########################################################################################################################
    def memSaver(self):     ### Cuts down attributes to bare minimum
        '''Cuts down attributes to bare minimum. Stores sequence name in self.info['FastaCmd']'''
        try:
            ### Basic Motif Details ###
            basic_info = {}
            for i in ['Match','SearchVar']:     # ['Match','SearchVar'] enough for Motif._hitStats(self,match,searchvar)
                basic_info[i] = self.getData(i,dlist=['info'])
            basic_stat = {}
            for s in ['Pos']:
                basic_stat[s] = self.getData(s,dlist=['stat'],str=False)
            ### Update attributes ##
            self.info = basic_info
            self.stat = basic_stat
            self.opt = {}
            self.dict = {'Data':{}}
            ### Clear Sequence ###
            if self.obj['Seq']:
                self.info['FastaCmd'] = self.obj['Seq'].shortName()
                self.obj['Seq'] = None
        except:
            self.log.errorLog('Error in MotifOcc.memSaver()',printerror=True,quitchoice=True)
#########################################################################################################################
    def occData(self):  ### Expands data dictionary for output
        '''Expands data dictionary for output.'''
        if self.obj['Motif']:
            self.dict['Data']['Motif'] = self.obj['Motif'].info['Name']
            self.stat['Len'] = len(self.obj['Motif'].info['Sequence'])
        else:
            self.stat['Len'] = len(self.getData('Match',default=''))
        if self.obj['Seq']:
            self.dict['Data']['Hit'] = self.obj['Seq'].shortName()
        else:
            self.dict['Data']['Hit'] = self.getData('FastaCmd')
        self.stat['Start_Pos'] = self.stat['Pos']
        self.stat['End_Pos'] = self.stat['Pos'] + self.stat['Len'] - 1
#########################################################################################################################
### End of SECTION II: MotifOcc Class                                                                                   #
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
        print '\n\n *** No standalone functionality! *** \n\n'
        
    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
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
