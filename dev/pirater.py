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
Module:       Pirater
Description:  Protein-protein Interaction Rater
Version:      0.0.0
Last Edit:    01/01/15
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program is a one-off rater of PPI data for picking up different types of interaction.

Commandline:

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
import rje, rje_db, rje_obj, rje_xref
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
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
    (program, version, last_edit, copy_right) = ('Pirater', '0.0', 'April 2015', '2015')
    description = 'Protein-protein Interaction Rater'
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class Pirater(rje_obj.RJE_Object):
    '''
    Pirater Class. Author: Rich Edwards (2015).

    Str:str
    
    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Pairwise']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
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
                self._cmdReadList(cmd,'file',['Pairwise'])  # String representing file path
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
            self.pirater()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ppidb = self.db('PPITypes',add=True,mainkeys=['UniA','UniB'])
            if not ppidb: ppidb = self.ppiTypes()
            if not ppidb: return False
            ### ~ [2] Load Pairwise PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Pairwise'): raise IOError('No pairwise PPI file given.')
            ppidb = db.addTable(self.getStr('Pairwise'),['Hub','Spoke','HubUni','SpokeUni','HubTaxID','SpokeTaxID','Evidence','IType'],['Hub','Spoke'])

            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def ppiTypes(self):  ### Mapping PPI types to experiment types & databases.
        '''Mapping PPI types to experiment types & databases.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()      # Might want to move some of this to the setup() method.
            xref = rje_xref.XRef(self.log,self.cmd_list)
            xref.setup()
            xdb = xref.db('xref')
            xdb.dropEntriesDirect('Status',['APPROVED'],inverse=True)
            ## Generate a list of Uniprot IDs to serve as background = all HGNC-mapped Uniprot entries.
            unilist = xdb.indexKeys('Uniprot')  # Make self.str['UniField']?
            ## MasterData directory
            mdir = rje.makePath('../../../MasterData-Apr15/data/')


            #?# Add ELM Class table with IC and filter #?#


            ### ~ [1] Table of Uniprot:Pfam links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            u2pfile = '%s2015-04-15.SLiMBenchSource/slimbench.Pfam.HUMAN.tdt' % mdir
            u2pdb = db.addTable(u2pfile,mainkeys=["Uniprot","Pfam"],name='Uni2Pfam')
            # Uniprot	Pfam
            # A0A087WTH1	PF04505
            u2pdb.renameField('Uniprot','PfamUni')

            ### ~ [2] ELM Interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            e2dfile = '%s2015-04-15.SLiMBenchSource/elm_interaction_domains.2015-04-15.tsv' % mdir
            e2ddb = db.addTable(e2dfile,mainkeys=["ELM identifier","Interaction Domain Id"],datakeys=["ELM identifier","Interaction Domain Id"],name='ELM2Pfam')
            e2ddb.renameField("ELM identifier",'ELM'); e2ddb.renameField("Interaction Domain Id",'Pfam')
            # "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
            # "CLV_NRD_NRD_1"	"PF00675"	"Peptidase_M16"	"Insulinase (Peptidase family M16)"
            # "CLV_PCSK_FUR_1"	"PF00082"	"Peptidase_S8"	"Subtilase family"
            e2ifile = '%s2015-04-15.SLiMBenchSource/elm_interactions.2015-04-15.tsv' % mdir
            e2idb = db.addTable(e2ifile,mainkeys=["Elm","Domain","interactorElm","interactorDomain"],name='ELM2PPI')
            e2idb.renameField("Elm",'ELM'); e2idb.renameField("Domain",'Pfam')
            e2idb.renameField("interactorElm",'ELMUni'); e2idb.renameField("interactorDomain",'PfamUni')
            # Elm	Domain	interactorElm	interactorDomain	StartElm	StopElm	StartDomain	StopDomain	AffinityMin	AffinityMax	PMID	taxonomyElm	taxonomyDomain
            # CLV_Separin_Fungi	PF03568	Q12158	Q03018	175	181	1171	1571	None	None	10403247,14585836	"559292"(Saccharomyces cerevisiae S288c)	"559292"(Saccharomyces cerevisiae S288c)
            # CLV_Separin_Fungi	PF03568	Q12158	Q03018	263	269	1171	1571	None	None	10403247	"559292"(Saccharomyces cerevisiae S288c)	"559292"(Saccharomyces cerevisiae S288c)
            e2pdb = db.copyTable(e2idb,'ELM2Domain')
            e2pdb.compress(['ELM','Pfam'])
            e2pdb.keepFields(['ELM','Pfam'])
            db.mergeTables(e2ddb,e2pdb)
            e2pdb = db.copyTable(e2idb,'ELM2PfamUni')
            # Make elmi and elmd - Need to add ELM instances
            # elmi = known ElmUni PfamUni PPI
            e2idb.compress(['ELMUni','PfamUni'])
            e2idb.newField('PPIType',evalue='elmi')
            e2idb.autoID()
            e2idb.setFields(['#','ELMUni','PfamUni','PPIType'])
            # elmp = known ElmUni Pfam
            e2pdb.compress(['ELMUni','Pfam'])
            e2pdb.dropFields(['ELM','PfamUni'])
            e2pdb.newField('PPIType',evalue='elmp')
            elmpdb = db.joinTables(name='elmp',join=[(e2pdb,'Pfam'),(u2pdb,'Pfam')],newkey=[],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            db.deleteTable(e2pdb)
            elmpdb.dropField('Pfam')
            elmpdb.compress(['ELMUni','PfamUni'])
            elmpdb.renameField('AutoID','#')
            elmpdb.setFields(['#','ELMUni','PfamUni','PPIType'])
            # elmd (== edmi) = ELM -> Known ElmUni - Pfam -> PfamUni
            # ndmi = ELM -> NoMask ElmUni - Pfam -> PfamUni
            # ddmi = ELM -> DisMask ElmUni - Pfam -> PfamUni
            # bdmi = ELM -> BothMask ElmUni - Pfam -> PfamUni
            for table in db.tables(): self.bugPrint('%s >> %s' % (table.name(),string.join(table.fields(),', ')))
            self.debug('>>>>')

            ### ~ [3] SLiMProb results of ELM search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdir = rje.makePath('../2015-04-25.HGNC-ELM-SLiMProb/')
            sdb = None
            for runid in ['NoMask','DisMask','BothMask']:
                dname = runid
                if not sdb: dname = 'slimprob'
                rdb = db.addTable('%shgnc-elm.%s.occ.csv' % (sdir,runid.lower()),mainkeys=['RunID','Motif','Seq'],datakeys=['RunID','Motif','Seq'],name=dname)
                # Dataset,RunID,Masking,Motif,Seq,Start_Pos,End_Pos,Prot_Len,Pattern,Match,Variant,MisMatch,Desc,UPC
                # CLV_C14_Caspase3-7.P42574.ppi,SLiMBench,FreqFT,CLV_C14_Caspase3-7,BID_HUMAN__P55957,27,31,195,[DEST][^P][^CDEFHWY]D[AGNS],SCSDN,[DEST][^P][^CDEFHWY]D[AGNS],0,RecName: Full=BH3-interacting domain death agonist; AltName: Full=p22 BID; Short=BID; Contains: RecName: Full=BH3-interacting domain death agonist p15; AltName: Full=p15 BID; Contains: RecName: Full=BH3-interacting domain death agonist p13; AltName: Full=p13 BID; Contains: RecName: Full=BH3-interacting domain death agonist p11; AltName: Full=p11 BID;,6
                # nomask, dismask and bothmask SLiMProb runs => *.occ.csv file
                # reduce to ELM-Uniprot pairs, keeping runID as list
                rdb.renameField('Motif','ELM'); rdb.renameField('Seq','ELMUni')
                rdb.compress(['ELM','ELMUni'],rules={'RunID':'list'})
                for entry in rdb.entries(): entry['ELMUni'] = string.split(entry['ELMUni'],'_')[-1]
                for entry in rdb.entries(): entry['RunID'] = '%sdmi' % entry['RunID'][:1].lower()
                rdb.renameField('RunID','PPIType')
                if sdb: db.mergeTables(sdb,rdb)
                else: sdb = rdb
            # Expand ELM and ELMUni to ELM, ELMUni, Pfam, then PfamUni
            tmpdb = db.joinTables(name='tmp',join=[(sdb,'ELM'),(e2ddb,'ELM',['Pfam'])],newkey=[],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            db.deleteTable(sdb)
            tmpdb.compress(['ELMUni','Pfam'])
            tmpdb.dropField('ELM')
            dmidb = db.joinTables(name='dmi',join=[(tmpdb,'Pfam'),(u2pdb,'Pfam')],newkey=[],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            dmidb.compress(['ELMUni','PfamUni'])
            dmidb.renameField('AutoID','#')
            dmidb.setFields(['#','ELMUni','PfamUni','PPIType'])
            # Update DMI fields
            for table in db.tables(): self.bugPrint('%s >> %s' % (table.name(),string.join(table.fields(),', ')))
            self.bugPrint('>>>>')
            for table in [elmpdb,e2idb,dmidb]:
                self.debug('%s >> %s' % (table.name(),string.join(table.fields(),', ')))
                table.newKey('#')
                for entry in table.entries():
                    if entry['ELMUni'] > entry['PfamUni']: (entry['ELMUni'],entry['PfamUni']) = (entry['PfamUni'],entry['ELMUni'])
            db.mergeTables(elmpdb,e2idb)
            db.mergeTables(dmidb,elmpdb)
            dmidb.renameField('ELMUni','UniA')
            dmidb.renameField('PfamUni','UniB')
            dmidb.setFields(['#','UniA','UniB','PPIType'])

            ### ~ [4] Table of DDI from 3DID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ddfile = '%s2015-04-16.PPI/3DID/3did_flat' % mdir
            # grep PF 3did_flat | head
            # #=ID	1-cysPrx_C	1-cysPrx_C	 (PF10417.4@Pfam	PF10417.4@Pfam)
            # #=ID	1-cysPrx_C	AhpC-TSA	 (PF10417.4@Pfam	PF00578.16@Pfam)
            dddb = db.addEmptyTable('DDI',['HubPfam','SpokePfam'],['HubPfam','SpokePfam'])
            for dline in open(ddfile,'r').readlines():
                if dline.startswith('#=ID'):
                    if rje.matchExp('^#=ID\s+\S+\s+\S+\s+\((PF\d+)\.\S+\s+(PF\d+)\.',dline):
                        (hubpfam,spokepfam) = rje.matchExp('^#=ID\s+\S+\s+\S+\s+\((PF\d+)\.\S+\s+(PF\d+)\.',dline)
                        dddb.addEntry({'HubPfam':hubpfam,'SpokePfam':spokepfam})
                    else: self.warnLog(rje.chomp(dline))
            dddb.renameField('HubPfam','Pfam')
            tdb = db.joinTables(name='temp',join=[(dddb,'Pfam'),(u2pdb,'Pfam')],newkey=[],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            tdb.dropField('Pfam')
            tdb.renameField('PfamUni','UniA')
            tdb.renameField('SpokePfam','Pfam')
            tdb.renameField('AutoID','OldID')   #!# Improve join with AutoID/# key!
            ddiddb = db.joinTables(name='ddi',join=[(tdb,'Pfam'),(u2pdb,'Pfam')],newkey=[],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            ddiddb.dropField('Pfam')
            ddiddb.renameField('AutoID','#')
            ddiddb.renameField('PfamUni','UniB')
            ddiddb.newField('PPIType',evalue='ddi')
            ddiddb.setFields(['#','UniA','UniB','PPIType'])

            ### ~ [5] Compile all data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Compress to HubUni-SpokeUni-PPIType data (where HubUni<SpokeUni to purge symmetry)
            # Keep "best" DMI type: e > b > d > n
            db.mergeTables(dmidb,ddiddb)
            dmidb.setStr({'Name':'PPITypes'})
            dmidb.compress(['UniA','UniB'])
            dmidb.setFields(['UniA','UniB','PPIType'])
            dmidb.saveToFile()
            return dmidb
        except: self.errorLog('%s.ppiTypes() error' % self.prog()); return None
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        There is currently no specific help available on REST output for this program. Run with &rest=help for general
        options. Run with &rest=full to get full server output. Individual outputs can be identified/parsed:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ...

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Main Pirater Methods                                                                                    #
#########################################################################################################################
    def pirater(self):  ### Mapping PPI types to experiment types & databases.
        '''Mapping PPI types to experiment types & databases.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppidb = self.db('PPITypes')

            # Join Tables of PPITypes and parsed PPI
            # Make index on PPITypes and Evidence
            # Look for overlaps

            return
        except: self.errorLog('%s.pirater() error' % self.prog())
#########################################################################################################################
### End of SECTION II: Pirater Class                                                                                    #
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
    try: Pirater(mainlog,['basefile=pirater']+cmd_list).run()

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
