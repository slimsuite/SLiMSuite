#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMDIP
Description:  Short Linear Motif-Domain Interaction Prediction
Version:      0.0.0
Last Edit:    15/12/15
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program combines the outputs from several other programs to produce a filtered list of predicted SLiM-mediated
    interactions in which:
    - hub proteins contain a SLiM-binding domain
    - spoke proteins contain a predicted SLiM
    - PPI data supports an interaction between hub and spoke for matching SLiMs and SLiM-binding domains

    The basic functionality of this program is a four-way join of data:


    A gene field is generated from the Seq field in occfile (the first element, split on "_"), which is used to cross-
    reference to the Hub field of the PPIFile.

Commandline:
    ### ~ SLiMDIP Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    occfile=DSVFILE : SLiM prediction (e.g. SLiMProb) occurrence file with Motif and Seq fields [None]
    ppifile=PPIFILE : Pairwise PPI file with HubUni and SpokeUni fields [None]
    domfile=DSVFILE : Protein-domain linkage file [None]
    dmifile=DSVFILE : Delimited text file linking Motifs with Proteins via SLiM-interaction domains [None]

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
import rje, rje_db, rje_obj, rje_slimlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
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
    # [ ] : Make more generic that HubUni and SpokeUni -> Seq mapping.
    # [ ] : Generate the DMI file
    # [ ] : Make sure DOMFILE can have Domain and Protein fields or recognise Uniprot xref output.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SLiMDIP', '0.0.0', 'December 2015', '2015')
    description = 'Generic RJE Module'
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
### SECTION II: SLiMDIP Class                                                                                           #
#########################################################################################################################
class SLiMDIP(rje_obj.RJE_Object):
    '''
    SLiMDIP Class. Author: Rich Edwards (2015).

    Str:str
    DomFile=DSVFILE : Protein-domain linkage file [None]
    DMIFile=DSVFILE : Delimited text file linking Motifs with Proteins via SLiM-interaction domains [None]
    OccFile=DSVFILE : SLiM prediction (e.g. SLiMProb) occurrence file with Motif and Seq fields [None]
    PPIFile=PPIFILE : Pairwise PPI file with HubUni and SpokeUni fields [None]

    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['DomFile','DMIFile','OccFile','PPIFile']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['SLiMList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        #self._setForkAttributes()   # Delete if no forking
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
                self._cmdReadList(cmd,'file',['DomFile','DMIFile','OccFile','PPIFile'])  # String representing file path
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
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.joinData()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## Set Basefile
            if not self.basefile(return_none=None): self.basefile(rje.baseFile(self.getStr('OccFile')))
            tabkeys = {'OccFile':['dataset','runid','motif','seq','start_pos','end_pos','variant'],
                       'DomFile':['domain','uniprot'],
                       'DMIFile':['motif','domain'],
                       'PPIFile':['hub','spoke']}
            ## Load Tables
            for dfile in ['DomFile','DMIFile','OccFile','PPIFile']:
                dbtable = db.addTable(self.getStr(dfile),mainkeys=tabkeys[dfile],name=dfile,expect=True,replace=False,uselower=True)
                self.tidyMotifNames(dbtable)
                if dfile == 'OccFile':
                    #dbtable.addField('uniprot')
                    dbtable.addField('gene')
                    for entry in dbtable.entries():
                        #entry['uniprot'] = string.split(entry['seq'],'_')[-1]  # Don't want this: uniprot is spoke!
                        entry['gene'] = string.split(entry['seq'],'_')[0]
                elif dfile == 'DomFile':
                    dbtable.compress(['domain','uniprot'],default='str')
                    dbtable.keepFields(['domain','uniprot'])
                elif dfile == 'DMIFile':
                    dbtable.compress(['motif','domain'],default='str')
                    dbtable.keepFields(['motif','domain'])
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def tidyMotifNames(self,dbtable):    ### Tidy the motif names in given dbtable
        '''Tidy the motif names in given dbtable.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slist = self.obj['SLiMList']
            if 'motif' not in dbtable.fields(): return 0
            mx = 0
            for entry in dbtable.entries():
                newname = slist.slimCoreName(entry['motif'])
                if newname != entry['motif']: entry['motif'] = newname; mx += 1
            self.printLog('#MOTIF','%s motif names corrected for SLiMList splitting.' % rje.iStr(mx))
            if mx: dbtable.remakeKeys()
            return mx
        except: self.errorLog('Problem during %s tidyMotifNames.' % self.prog()); raise
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def joinData(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            odb = self.db('OccFile')
            ### ~ [2] Join Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #>> name:str [''] = Name for new table. If not given will become "TableX"
            #>> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
                #formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
                #new table. If Field does not exist, it will be added. (Field may be a Formula.)
            #>> newkey:list [] = If None, will make a new "AutoID" key field. Else, will use the given Fields.
            #>> cleanup:bool [True] = If True will delete any Fields generated just to make the join
            #>> delimit:str ['\t'] = Delimiter to be used to join the key fields
            #>> empties:bool [True] = Whether to keep entries that do not link to 1+ tables with empty values or delete them.
            #>> check:bool [False] = Whether to check for entries that are not being joined.
            #>> keeptable:bool [True] = Whether to add new table to self.list['Tables']
            okeys = odb.keys()
            okeys += ['domain']
            self.printLog('#JOIN','Joining occurrence table to DMI via motif')
            self.bugPrint(odb.fields())
            self.bugPrint(self.db('DMIFile').fields())
            self.debug(odb.index('motif')['LIG_CtBP_PxDLS_1'])
            odb.addField('motupper')
            for entry in odb.entries(): entry['motupper'] = entry['motif'].upper()
            odb = db.joinTables(name='occ-dmi',join=[(odb,'motupper'),(self.db('DMIFile'),'motif',['domain'])],newkey=okeys,empties=False)
            odb.dropField('motupper')
            self.debug(odb.index('motif')['LIG_CtBP_PxDLS_1'])
            okeys.insert(-1,'uniprot')
            self.bugPrint(okeys)
            self.printLog('#JOIN','Joining table to DMI uniprot via domain')
            self.bugPrint(odb.fields())
            self.bugPrint(self.db('DomFile').fields())
            ddb = db.joinTables(name='occ-dmi-dom',join=[(odb,'domain'),(self.db('DomFile'),'domain',['uniprot'])],newkey=okeys,empties=False)
            self.debug(ddb.index('motif')['LIG_CtBP_PxDLS_1'])
            okeys = okeys[:-1]
            ddb.compress(okeys,rules={'domain':'list'},default='str')
            self.printLog('#JOIN','Joining table to PPI via gene|uniprot pairs')
            self.bugPrint(ddb.fields())
            self.bugPrint(self.db('PPIFile').fields())
            dip = db.joinTables(name='slimdip',join=[(ddb,'#gene#|#uniprot#'),(self.db('PPIFile'),'#hub#|#spokeuni#')],newkey=okeys,empties=False,cleanup=True)
            self.debug(dip.index('motif')['LIG_CtBP_PxDLS_1'])
            dip.saveToFile()
            ### ~ [3] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Reduce to unique sets of data of interest #!#
            unique = ['dataset','runid','motif','domain','hub','spoke']
            extra = self.db('PPIFile').fields()[2:]
            dip.compress(unique)
            dip.keepFields(unique+extra)
            dip.setStr({'Name':'dmi'})
            dip.saveToFile()
            #dip.dropEntries(['gene!=hub'])
            return dip
        except: self.errorLog('%s.method error' % self.prog())
#########################################################################################################################
### End of SECTION II: SLiMDIP Class                                                                                    #
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
    try: SLiMDIP(mainlog,cmd_list).run()

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
