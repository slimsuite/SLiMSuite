#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_hivqsf
Description:  HIV QSLiMFinder analysis pipeline
Version:      0.0
Last Edit:    19/07/11
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    hhpid=FILE      : Delimited text file containing HHPID interactions []
    hivseq=FILE     : File containing HIV sequences []
    genemap=FILE    : Delimited human gene mapping data []
    ppidir=PATH     : Path to human protein PPI files [] 

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seq, rje_zen
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
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_HIVQSF', '0.0', 'July 2011', '2011')
    description = 'HIV QSLiMFinder analysis pipeline'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
class HIVQSF(rje_obj.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2011).

    Str:str
    - GeneMap = Delimited human gene mapping data []
    - HHPID = Delimited text file containing HHPID interactions []
    - HIVSeq = File containing HIV sequences []
    - Pairwise = Pairwise Human PPI Table
    - PPIDir = Path to human protein PPI files [] 
    
    Bool:boolean

    Int:integer

    Num:float
    
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
        self.strlist = ['GeneMap','HHPID','Pairwise','HIVSeq','PPIDir']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                #self._forkCmd(cmd)  # Delete if no forking
                ### Class Options ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'file',['GeneMap','HHPID','Pairwise','HIVSeq'])
                self._cmdReadList(cmd,'path',['PPIDir'])  
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.processHHPID() and os.path.exists(self.getStr('PPIDir')): self.makePPI()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.db().addTable(self.getStr('GeneMap'),mainkeys=['Gene'],datakeys='All',name='GeneMap')
            self.db().addTable(self.getStr('Pairwise'),mainkeys=['Hub','Spoke'],datakeys='All',name='PPI')
            self.loadHHPID()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### HHPID Methods                                                                                           #
#########################################################################################################################
    def loadHHPID(self):    ### Load HHPID interactions
        '''Load HHPID interactions.'''
        try:### ~ [1] Setup HHPID Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('HHPID').lower() in ['','none']: return self.printLog('#HHPID','No HHPID file to load')
            hdb = self.db().addTable(self.getStr('HHPID'),mainkeys='auto',datakeys='All',name='HHPID')
            for field in ['#Tax ID 1','Tax ID 2','product accession.version 2','last update timestamp']: hdb.dropField(field)
            hdb.renameField('Gene ID 1','EntrezHIV')
            hdb.renameField('product accession.version 1','AccHIV')
            hdb.renameField('product name 1','HIV')
            hdb.renameField('Interaction short phrase','Interaction')
            hdb.renameField('Gene ID 2','Entrez')
            hdb.renameField('product name 2','Description')
            hdb.renameField('PubMed ID (PMID) list','PMID')
            for itype in rje.sortKeys(hdb.index('Interaction')): self.printLog('#HHPID','%s => %s entries' % (itype, len(hdb.index('Interaction')[itype])))
            hdb.dropEntriesDirect('Interaction',['binds','complexes with','interacts with'],inverse=True)
            return True
        except: self.errorLog('%s.loadHHPID error' % self)
#########################################################################################################################
    def processHHPID(self): ### Process HHPID interactions
        '''Process HHPID interactions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile('%s.HHPIDMap.tdt' % self.basefile()):
                mdb = self.db().addTable('%s.HHPIDMap.tdt' % self.basefile(),['HIV','Gene'],'All',name='HHPIDMap')
                return mdb
            hdb = self.db('HHPID')
            gdb = self.db('GeneMap')
            pdb = self.db('PPI')
            mdb = self.db().joinTables(name='HHPIDMap',join=[(hdb,'Entrez'),(gdb,'Entrez')],newkey=['#'],empties=False,keeptable=True)
            for field in mdb.fields()[0:]:
                if field not in ['#','AccHIV','EntrezHIV','HIV','Entrez','Gene','Symbol','UniProt','EnsEMBL','EnsLoci']: mdb.dropField(field)
            mdb.compress(['HIV','Gene'],default='str'); mdb.dropField('#')
            mdb.saveToFile()
            ### ~ [2] Save viral accession numbers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            open('%s.hivacc' % self.getStr('Basefile'),'w').write('%s\n' % string.join(rje.sortKeys(mdb.index('AccHIV')),'\n'))
            return mdb
        except: self.errorLog('%s.processHHPID error' % self); return False
#########################################################################################################################
    def makePPI(self):  ### Generates files for Human-HIV PPI analysis
        '''Generates files for Human-HIV PPI analysis.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('HIVSeq'),'autoload=T'])
            if not seqlist.seqs(): return False
            seqmap = seqlist.seqNameDic('Max')
            mdb = self.db('HHPIDMap')
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hivacc in mdb.index('AccHIV'):
                # map HIV accession numbers on to sequences seqNameDic
                accnum = string.split(hivacc,'.')[0]
                hivseq = seqmap[accnum]              
                # extract short HIV name from sequence ID
                hivgene = string.split(hivseq.shortName(),'_')[0].upper()
                # create directory named after HIV gene
                #self.progLog('\r#PPI','Generating human-HIV PPI fasta files for %s' % (hivgene))
                rje.mkDir(self,'%s/' % hivgene,log=True)
                # copy human PPI files into directories, adding HIV gene
                ex = 0.0; etot = len(mdb.index('AccHIV')[hivacc])
                for entry in mdb.indexEntries('AccHIV',hivacc):
                    self.progLog('\r#PPI','Generating human-HIV PPI fasta files for %s %s PPI' % (rje.iStr(etot),hivgene))
                    pfile = self.getStr('PPIDir') + entry['Symbol'] + '.ppi.fas'
                    if rje.exists(pfile):
                        FAS = open('%s/%s.%s.ppi.fas' % (hivgene,hivgene.lower(),entry['Symbol']),'w')
                        FAS.write('>%s\n%s\n' % (hivseq.info['Name'],hivseq.getSequence()))
                        FAS.write(open(pfile,'r').read())
                        FAS.close()
                    else: self.errorLog('Cannot find human PPI file for %s interactor "%s"' % (entry['HIV'],entry['Symbol']),printerror=False)
                self.printLog('\r#PPI','Generated human-HIV PPI fasta files for %s %s (%s) PPI.' % (rje.iStr(etot),entry['HIV'],hivgene))                                      
        except: self.errorLog('%s.makePPI error' % self); return False
#########################################################################################################################
### End of SECTION II: HIVQSF Class                                                                                     #
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
    try: HIVQSF(mainlog,['basefile=hiv_qsf']+cmd_list).run()

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
