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
Module:       rje_slimfungo
Description:  SLiM Functional GO classification module
Version:      0.0
Last Edit:    28/03/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module sticks together GO categories with output from SLiMSearch to try and assess functional significance
    within categories.

Commandline:
    occdata=FILE    : File of motif occurrences. Must have  [None]
    gomap=FILE      : Mapping of sequences to GO categories [None]
    minocc=X        : Min occurrences for a motif in a GO category [6]

Optional:
    seqin=FILE      : File of input sequences [None]

RJE_GO Commandline:    
    obofile=FILE    : Input GO OBO V1.2 download [None]
    webobo=T/F      : Whether to download from GO website if file not given [True]
    goslim=LIST     : List of GO IDs to form basis of GO Slim []

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_seq, rje_go
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_go, rje_zen
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
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_SLiMFunGo', '0.0', 'June 2008', '2008')
    description = 'SLiM Functional GO classification module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
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
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMFunGO Class                                                                                     #
#########################################################################################################################
class SLiMFunGO(rje.RJE_Object):     
    '''
    SLiMFunGO Class. Author: Rich Edwards (2008).

    Info:str
    - OccData = File of motif occurrences. Must have  [None]
    - GOMap = Mapping of sequences to GO categories [None]

    Opt:boolean

    Stat:numeric
    - MinOcc = Min occurrences for a motif in a GO category [6]
    
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
        self.infolist = ['OccData','GOMap']
        self.optlist = []
        self.statlist = ['MinOcc']
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'MinOcc':6})
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
                self._cmdReadList(cmd,'file',['OccData','GOMap'])  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'int',['MinOcc'])  # No need for arg if arg = att.lower()
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mygo = rje_go.GO(self.log,self.cmd_list)
            mygo.readGO()
            gomap = rje.dataDict(self,self.info['GOMap'],mainkeys=['Ensembl Gene ID'],datakeys=['GO ID'],lists=True)
            self.deBug(rje.sortKeys(gomap)[:100])
            #!# Replace 'Ensembl Gene ID' with commandline parameter at some point #!#
            self.printLog('#GOMAP','Loaded GO mappings for %s sequence IDs' % (rje.integerString(len(gomap))))
            slimocc = rje.dataDict(self,self.info['OccData'],mainkeys=['Motif','Seq','Start_Pos','End_Pos'],datakeys=['Motif','Seq','Start_Pos','End_Pos','Cons','HomNum'])
            self.printLog('#OCC','Loaded Data for %s motif occurrences.' % (rje.integerString(len(slimocc))))
            ## ~ [1a] ~ Sequence mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list)
            seqmap = {}
            (sx,stot) = (0.0,seqlist.seqNum())
            for seq in seqlist.seq:
                self.progLog('#SEQMAP','Mappings sequence IDs: %.1f%%' % (sx/stot)); sx += 100.0
                if rje.matchExp('gene:(\S+)\]',seq.info['Name']): seqmap[seq.shortName()] = rje.matchExp('gene:(\S+)\]',seq.info['Name'])[0]
            self.printLog('\r#SEQMAP','Mappings %s sequence IDs complete: %s mapped' % (rje.integerString(stot),rje.integerString(len(seqmap))))
            self.deBug(rje.sortKeys(seqmap)[:100])

            ### ~ [2] ~ Output new data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goocc = {}
            outfile = string.join(string.split(self.info['OccData'],'.')[:-1] + ['slimfungo','tdt'],'.')
            headers = ['GO','Motif','Type','Seq','Start_Pos','End_Pos','Cons','HomNum']
            for okey in slimocc.keys():
                self.progLog('#NEW','Making new GO occurrences: %s    ' % (rje.integerString(len(slimocc))))
                data = slimocc.pop(okey)
                gene = seq = data['Seq']
                type = 'fwd'
                if string.split(data['Motif'],'_')[-1] in ['rev','scram']:
                    type = string.split(data['Motif'],'_')[-1]
                    data['Motif'] = string.join(string.split(data['Motif'],'_')[:-1],'_')
                if gene not in gomap and gene in seqmap: gene = seqmap[gene]
                golist = []
                if gene in gomap:
                    for id in gomap[gene]: golist += mygo.parents(id)
                else: golist = ['NoGo']
                self.deBug('%s:%s::%s' % (seq,gene,golist))
                for id in rje.sortUnique(golist,False,False):
                    if id not in goocc: goocc[id] = {}
                    if motif not in goocc[id]: goocc[id][motif] = {'fwd':[],'rev':[],'scram':[]}
                    goocc[id][motif][type].append(rje.combineDict({'GO':id,'Type':type},data))
            self.printLog('\r#NEW','Making new GO occurrences complete.    ' % (rje.integerString(len(slimocc))))

            rje.delimitedFileOutput(self,outfile,headers,rje_backup=True)
            (mx,ox,ix,itot) = (0,0,0.0,len(goocc))
            for id in rje.sortKeys(goocc):
                for motif in rje.sortKeys(goocc[id]):
                    for type in rje.sortKeys(goocc[id][motif]):
                        if len(goocc[id][motif][type] < self.stat['MinOcc']): goocc[id][motif].pop(type)
                    if len(goocc[id][motif]) < 2 or 'fwd' not in goocc[id][motif]: continue
                    mx += 1
                    for type in goocc[id][motif]:
                        for occ in goocc[id][motif][type]: rje.delimitedFileOutput(self,outfile,headers,datadict=occ); ox += 1
                self.progLog('#OUT','Output to %s: %.2f%% :: %s motifs; %s occ.' % (outfile,ix/itot,rje.integerString(mx),rje.integerString(ox)))
            self.printLog('\r#OUT','Output of occurrences to %s is now complete: %s motifs; %s occ.' % (outfile,rje.integerString(mx),rje.integerString(ox)))

        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: SLiMFunGO Class                                                                                  #
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
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:SLiMFunGO(mainlog,cmd_list).run()
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
