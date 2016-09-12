#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_spf
Description:  SPF Level Extraction Tool
Version:      0.1.0
Last Edit:    10/04/16
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module converts an SPF file into a file per taxonomic level, avoiding issues of blank taxa.

Commandline:
    infile=FILE : The input SPF file [input.spf]
    log=FILE    : The name of the output log (*.log) [rje_spf.log]

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.1.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [ ] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_SPF', '0.1.0', 'April 2016', '2016')
    description = 'SPF Level Extraction Tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
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
### SECTION II: SPF Class                                                                                               #
#########################################################################################################################
class SPF(rje_obj.RJE_Object):
    '''
    Class. Author: Rich Edwards (2015).

    Str:str
    - InFile=FILE : The input SPF file [input.spf]

    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
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
        self.strlist = ['InFile']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'InFile':'input.spf'})
        self.setBool({})
        self.setInt({})
        self.setNum({})
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
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['InFile'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infile = self.getStr('InFile')
            while not rje.exists(infile):
                infile = rje.choice('File "%s" not found. Input file name? (Blank to quit):' % infile)
                if not infile: return self.printLog('#QUIT','Execution terminated!')
            db = rje_db.Database(self.log,self.cmd_list)
            db.basefile(rje.baseFile(infile))
            sdb = db.addTable(infile,mainkeys='#',delimit='\t',name='SPF.Mod')
            levels = {'Level_1':'k','Level_2':'p','Level_3':'c','Level_4':'o','Level_5':'f','Level_6':'g','Level_7':'s'}
            # k__Bacteria	p__Proteobacteria	c__Alphaproteobacteria	o__Rhodospirillales	f__Rhodospirillaceae	g__	s__	denovo44
            # Unassigned	unclassified	unclassified	unclassified	unclassified	unclassified	unclassified	denovo49
            ### ~ [1] Modify Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dupnames = []
            parents = {}    # Parent for each term
            renamed = []
            ex = 0.0; etot = sdb.entryNum()
            for entry in sdb.entries():
                self.progLog('\r#SPF','Modifying SPF content: %.1f%%' % (ex/etot)); ex += 100.0
                taxon = ''
                parent = ''
                #self.debug(entry)
                for lvl in ['Level_1','Level_2','Level_3','Level_4','Level_5','Level_6','Level_7']:
                    entry[lvl] = string.replace(entry[lvl],'unidentified','unclassified')
                    #entry[lvl] = string.replace(entry[lvl],'Incertae_sedis','Incertae_sedis-%s' % levels[lvl])
                    null = '%s__' % levels[lvl]
                    #self.bugPrint(null)
                    #self.bugPrint(entry[lvl])
                    if entry[lvl] in [null,'Unassigned','unclassified','%sunclassified' % null,'%sunidentified' % null,'%sunculturedfungus' % null,'%sIncertae_sedis' % null,'%sunclassified_sp.' % null]:
                        if not taxon or taxon.endswith('unclassified'): entry[lvl] = '%sunclassified' % null
                        #elif taxon.endswith('unassigned)'): entry[lvl] = '%s%s' % (null,taxon[3:])
                        #elif taxon.endswith('unassigned)'): entry[lvl] = '%s(%s;%s-unassigned)' % (null,string.split(taxon,'(')[1][:-1],levels[lvl])
                        elif taxon.endswith('unassigned)'): entry[lvl] = '%s%s;%s-unassigned)' % (null,taxon[3:][:-1],levels[lvl])
                        else: entry[lvl] = '%s%s(%s-unassigned)' % (null,taxon[3:],levels[lvl])
                    if entry[lvl] in parents:
                        #self.debug(parents[entry[lvl]])
                        if parent in parents[entry[lvl]]: entry[lvl] = parents[entry[lvl]][parent]
                        else:
                            self.bugPrint(entry[lvl])
                            self.bugPrint(parents[entry[lvl]])
                            renamed.append(entry[lvl])
                            newtax = '%s%d' % (entry[lvl],renamed.count(entry[lvl]))
                            self.warnLog('%s had multiple parents (%s & %s) -> %s' % (entry[lvl],string.join(parents[entry[lvl]],'|'),parent,newtax))
                            parents[newtax] = {parent:newtax}
                            parents[entry[lvl]][parent] = newtax
                            entry[lvl] = newtax
                            self.deBug(parents[entry[lvl]])
                    elif parent: parents[entry[lvl]] = {parent:entry[lvl]}
                    parent = entry[lvl]
                    if entry[lvl][3:] == taxon[3:]:
                        if (entry[lvl],taxon) not in dupnames: dupnames.append((entry[lvl],taxon))
                    #self.bugPrint(entry[lvl])
                    taxon = entry[lvl]
                #self.debug(entry)
                #self.debug(parents)
            self.printLog('\r#SPF','Modifying SPF content complete.')
            dupnames.sort()
            for (dupA,dupB) in dupnames: self.warnLog('Duplicate taxa names: %s & %s' % (dupA,dupB))
            ### ~ [2] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.saveToFile(savefields=sdb.list['Fields'][1:])
            ### ~ [3] Compress to different taxonomic levels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            compress = ['Level_1','Level_2','Level_3','Level_4','Level_5','Level_6','Level_7','#']
            dump = compress.pop(-1)
            rules = {'Observation Ids':'list',dump:'str'}
            sdb.dropField('Observation Ids')
            while compress:
                sdb.compress(compress,rules=rules,default='sum',best=[],joinchar='|')
                #if dump == '#':
                sdb.dropField(dump)
                sdb.saveToFile('%s.SPF.%s.%s.spf' % (rje.baseFile(infile),compress[-1],levels[compress[-1]]))
                dump = compress.pop(-1); rules[dump] = 'list'
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: SPF Class                                                                                        #
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
    try: SPF(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
