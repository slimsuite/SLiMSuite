#!/usr/bin/python

# rje_blast - BLAST Control Module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_blast
Description:  BLAST Control Module
Version:      1.16.0
Last Edit:    21/08/20
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Performs BLAST searches and loads results into objects. Peforms GABLAM conversion of local alignments into global
    alignment statistics. Remember to set blastp=X:  blastx for DNA vs prot; tblastn for Prot vs DNA)

Objects:
    BLASTRun = Full BLAST run
    BLASTSearch = Information for a single Query search within a BLASTRun
    BLASTHit = Detailed Information for a single Query-Hit pair within BLASTRun
    PWAln = Detailed Information for each aligned section of a Query-Hit Pair

Commandline:
    blastpath=X     : path for blast files [''] (Use fwd slashes)
    
    blastp=X        : BLAST program (BLAST -p X) [blastp]
    blasti=FILE     : Input file (BLAST -i FILE) [None]
    blastd=FILE     : BLAST database (BLAST -d FILE) [None]
    formatdb=T/F    : Whether to (re)format BLAST database [False]
    blasto=FILE     : Output file (BLAST -o FILE) [*.blast]

    blaste=X        : E-Value cut-off for BLAST searches (BLAST -e X) [1e-4]
    blastv=X        : Number of one-line hits per query (BLAST -v X) [500]
    blastb=X        : Number of hit alignments per query (BLAST -b X) [250]  

    blastf=T/F      : Complexity Filter (BLAST -F X) [True]
    blastcf=T/F     : Use BLAST Composition-based statistics (BLAST -C X) [False]
    blastg=T/F      : Gapped BLAST (BLAST -g X) [True]

    blasta=X        : Number of processors to use (BLAST -a X) [1]
    blastopt=FILE   : File containing raw BLAST options (applied after all others) []
    ignoredate=T/F  : Ignore date stamps when deciding whether to regenerate files [False]

    gablamfrag=X    : Length of gaps between mapped residue for fragmenting local hits [100]
    localcut=X      : Cut-off length for local alignments contributing to global GABLAM stats) [0]

Uses general modules: os, re, string, sys, time
Uses RJE modules: rje
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Working Compilation.
    # 0.1 - No Out Object in Objects
    # 1.0 - Corrected to work with blastn (and blastp)
    # 1.1 - Added special calling for Cerberus
    # 1.2 - Added GABLAM and GABLAMO to BlastHit
    # 1.3 - Added GABLAM calculation upon reading BLAST results and clearing Alignment sequences to save memory
    # 1.4 - Tidied up the module with improved logging and progress reporting. Added dbCleanup.
    # 1.5 - Added checking for multiple hits with same name and modified BLAST_Run.hitToSeq()
    # 1.6 - Added nucleotide vs protein searches to GABLAM
    # 1.7 - Added nucleotide vs nucleotide searches to GABLAM
    # 1.8 - Added local alignment summary output to ReadBLAST()
    # 1.9 - Added BLAST -C
    # 1.10- Added BLAST -g
    # 1.11- Added gablamfrag=X : Length of gaps between mapped residue for fragmenting local hits [100]
    # 1.12- Altered checkDB and cleanupDB to spot index files split over multiple files (*.00.p* etc.)
    # 1.13- Added localcut=X : Cut-off length for local alignments contributing to global GABLAM stats) [0]
    # 1.14- Added blast.checkProg(qtype,stype) to check whether blastp setting matches sequence formats.
    # 1.15- Added OldBLAST/Legacy option to Object for compatibility with rje_blast_V2. (Always True!)
    # 1.16.0 - Initial Python3 code conversion.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Output of alignment with options for line lengths and numbers +-
    # [ ] : OrderAln by any stat ['BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd']
    # [Y] : Optionfile and commandline options blast-e=X etc.
    # [ ] : Test with other blast programs
    # [ ] : Add standalone running for blast searching?
    # [ ] : Replace blast.search with list['Search'] etc?
    # [ ] : Add documentation for implementation details - each method?
    # [ ] : Locate which classes/methods call BLASTRun.hitToSeq and look to improve reporting etc.
    # [ ] : Fix DNA implementation of GABLAM to allow Ordered GABLAM in either direction.
    # [ ] : Add positional information to GABLAM dictionary - start and end of aligned portions
    # [ ] : Add oritentation of Query and Hit for DNA GABLAM
    # [ ] : Check/fix the database format checking of DNA databases
    # [Y] : Update to new Module Structre (V2.0)
    # [ ] : Check and Tidy this To Do list!
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_BLAST', '1.16.0', 'August 2018', '2005')
    description = 'RJE BLAST Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please contact seqsuite@gmail.com if this breaks - BLAST changes sometimes!']
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
    except: print('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: BLASTRun Class                                                                                          #
#########################################################################################################################
class BLASTRun(rje.RJE_Object):     
    '''
    BLASTRun Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Output File Name (BLAST -o X)
    - Type = Program (blastp etc.) (BLAST -p X)
    - DBase = Database to search (BLAST -d X)
    - InFile = Input file name (BLAST -i X)
    - OptionFile = file containing a string of BLAST options to append to commandline
    - BLAST Path = path to blast programs
    - BLASTCmd = system command used to generate BLAST in self.blast()
    - BLASTOpt = Additional BLAST options
    
    Opt:boolean
    - Composition Statistics
    - Complexity Filter = (BLAST -F) [True]
    - FormatDB = whether to (re)format database before blasting
    - GappedBLAST = Gapped BLAST (BLAST -g X) [True]
    - IgnoreDate = Ignore date stamps when deciding whether to regenerate files [False]

    Stat:numeric
    - E-Value = e-value cut-off (BLAST -e X) [1e-4]
    - OneLine = Number of one-line hits per query (BLAST -v X) [500]
    - HitAln  = Number of hit alignments per query (BLAST -b X) [250]
    - DBLen = Length of Database (letters)
    - DBNum = Number of Sequences in Database
    - BLASTa = Number of processors to use (BLAST -a X) [1]

    List:list

    Dict:dictionary    

    Obj:RJE_Objects

    Other:
    search:list = list of BLASTSearch Objects
    '''
    ### Attributes
    search = [] # List of BLASTSearch Objects   #!# Also change rje_hmm if changing this #!#
    def searchNum(self): return len(self.search)
    def legacy(self): return True
    def oldBLAST(self): return True
    def hitNum(self):
        hx = 0
        for search in self.search: hx += search.hitNum()
        return hx
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Name','Type','DBase','InFile','OptionFile','BLAST Path','BLASTCmd','BLASTOpt']
        self.statlist = ['E-Value','OneLine','HitAln','DBLen','DBNum','BLASTa']
        self.optlist = ['Complexity Filter','FormatDB','Composition Statistics','GappedBLAST','IgnoreDate']
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'BLAST Path':'','Type':'blastp'})
        self.setOpt({'FormatDB':False,'Composition Statistics':False,'Complexity Filter':True,'IgnoreDate':False,
                     'OldBLAST':True})
        self.setStat({'E-Value':0.0001,'OneLine':500.0,'HitAln':250.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.search = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] ~ Standard Command Line processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdRead(type='path',att='BLAST Path',arg='blastpath',cmd=cmd)
                self._cmdRead(type='info',att='Name',arg='blasto',cmd=cmd)
                self._cmdRead(type='info',att='Type',arg='blastp',cmd=cmd.lower())
                self._cmdRead(type='info',att='DBase',arg='blastd',cmd=cmd)
                self._cmdRead(type='info',att='InFile',arg='blasti',cmd=cmd)
                self._cmdRead(type='opt',att='Complexity Filter',arg='blastf',cmd=cmd)
                self._cmdRead(type='opt',att='Composition Statistics',arg='blastcf',cmd=cmd)
                self._cmdReadList(cmd,'opt',['FormatDB','IgnoreDate'])
                self._cmdReadList(cmd,'int',['BLASTa'])
                self._cmdRead(type='opt',att='GappedBLAST',arg='blastg',cmd=cmd)
                self._cmdRead(type='stat',att='E-Value',arg='blaste',cmd=cmd)
                self._cmdRead(type='int',att='OneLine',arg='blastv',cmd=cmd)
                self._cmdRead(type='int',att='HitAln',arg='blastb',cmd=cmd)
                self._cmdRead(type='file',att='OptionFile',arg='blastopt',cmd=cmd)
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [2] ~ Special ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['FormatDB']:
            if self.info['Type'] in ['blastn','tblastn','tblastx']: self.formatDB(protein=False)
            else: self.formatDB()
        return
#########################################################################################################################
    ### <2> ### BLAST Search                                                                                            #
#########################################################################################################################
    def formatDB(self,fasfile=None,protein=True,force=True,log=True,checkage=None):    ### BLAST formats database given
        '''
        BLAST formats database given.
        >> fasfile:str = Name of file to form database [None]
        >> protein:boolean = whether protein [True]
        >> force:boolean = whether to overwrite an existing formatted DB [True]
        '''
        ### ~ [1] Check and update DBase as necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if checkage == None: checkage = not self.getAttribute('opt','IgnoreDate',default=False)
        if fasfile: self.info['DBase'] = fasfile
        elif self.info['DBase'] == 'None': self.printLog('#ERR','Cannot format database "None"!'); raise ValueError
        ### ~ [2] Call formatdb if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if force or not checkForDB(dbfile=self.info['DBase'],checkage=checkage,log=self.log,protein=protein):
            if log: formatDB(self.info['DBase'],self.info['BLAST Path'],protein,self.log)
            else: formatDB(self.info['DBase'],self.info['BLAST Path'],protein)
        else: self.printLog('#DB ','%s already formatted. (Force = False).' % self.info['DBase'],log=log)
#########################################################################################################################
    def checkProg(self,qtype='Unknown',stype='Unknown',log=True): ### Checks sequence types against blastp
        '''
        Checks sequence types against blastp.
        >> qtype:str [None] = Query sequence type
        >> stype:str [None] = SearchDB sequence type
        >> log:bool [True] = Whether to output result of check to log
        '''
        try:### ~ [1] ~ Check BLASTP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qtype = qtype.lower(); stype = stype.lower()
            problem = False
            if self.info['Type'] in ['blastp','blastx'] and stype[1:] == 'na': problem = True
            elif self.info['Type'] in ['blastp','tblastn'] and qtype[1:] == 'na': problem = True
            elif self.info['Type'] in ['blastn','tblastn','tblastx'] and stype[:4] in ['prot','pept']: problem = True
            elif self.info['Type'] in ['blastn','blastx','tblastx'] and qtype[:4] in ['prot','pept']: problem = True
            if not problem:
                if log: self.printLog('#PROG','BLAST Program OK (%s): %s vs %s' % (self.info['Type'],qtype,stype))
                return True
            ### ~ [2] ~ Predict BLASTP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blastp = ''
            if qtype[:4] in ['prot','pept']:
                if stype[:4] in ['prot','pept']: blastp = 'blastp'
                elif stype[1:] == 'na': blastp = 'tblastn'
            elif qtype[1:] == 'na':
                if stype[:4] in ['prot','pept']: blastp = 'blastx'
                elif stype[1:] == 'na': blastp = 'blastn'
            if blastp:
                if self.i() < 0 or rje.yesNo('%s vs %s: switch BLAST program to %s?' % (qtype,stype,blastp)):
                    self.info['Type'] = blastp
                    if log: self.printLog('#PROG','BLAST Program changed (%s): %s vs %s' % (self.info['Type'],qtype,stype))
                    return True
            else: blastp = 'Unknown'
            self.printLog('#ERR','WARNING: %s vs %s and blastp=%s (Recommended: %s)' % (qtype,stype,self.info['Type'],blastp))
            return False
        except: self.errorLog('Something went wrong with blast.checkProg()'); return False
#########################################################################################################################
    def blast(self,wait=True,type=None,cleandb=False,use_existing=False,log=False):    ### Performs BLAST using object attributes
        '''
        Performs BLAST using object attributes.
        >> wait:boolean  = whether to wait for BLAST. [True]
        >> type:str = type of BLAST search [None]
        >> cleandb:bool = whether to cleanup (delete) searchDB files after search [False]
        >> use_existing:bool = if True, will check for existing result and use if newer than files
        >> log:bool = Whether to log BLAST run
        '''
        try:### ~ [1] ~ Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Check for Existing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if use_existing and self.checkBLAST(logcheck=False):                ### BLAST Results exist
                if self.getAttribute('opt','IgnoreDate',default=False): return  ### Don't check age!
                if rje.isYounger(self.info['DBase'],self.info['Name']) == self.info['Name'] and rje.isYounger(self.info['InFile'],self.info['Name']) == self.info['Name']: return
            ## ~ [1b] ~ Setup BLAST Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blastpath = self.info['BLAST Path']
            if self.stat['OneLine'] < self.stat['HitAln']:
                self.printLog('#CMD','Reduced reported alignments from %s to match one-line reporting number of %s.' % (rje.integerString(self.stat['HitAln']),rje.integerString(self.stat['OneLine'])))
                self.stat['HitAln'] = self.stat['OneLine']
            ### ~ [2] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Make system command for BLAST call ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            command = blastpath + 'blastall'
            if type: self.info['Type'] = type
            command += ' -p %s' % self.info['Type']
            command += ' -i %s' % self.info['InFile']
            command += ' -d %s' % self.info['DBase']
            command += ' -o %s' % self.info['Name']
            command += ' -F %s' % str(self.opt['Complexity Filter'])[0] 
            try: command += ' -C %s' % str(self.opt['Composition Statistics'])[0]
            except: pass    # For pickled BLASTs
            command += ' -e %e' % self.stat['E-Value']
            command += ' -v %d' % self.stat['OneLine']
            command += ' -b %d' % min(self.stat['HitAln'],self.stat['OneLine'])
            if not self.opt['GappedBLAST']: command += ' -g F'
            if self.stat['BLASTa'] > 1: command += ' -a %d' % self.stat['BLASTa']
            if self.info['BLASTOpt'].lower() not in ['','none']: command = '%s %s' % (command,self.info['BLASTOpt'])
            if rje.exists(self.info['OptionFile']):
                for line in open(self.info['OptionFile'],'r').readlines(): command = '%s %s' % (command,rje.chomp(line))
            if not wait: command += ' &'
            ##~ [2b] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.info['BLASTCmd'] = command
            if log: self.printLog('\r#SYS',command)
            elif self.stat['Verbose'] > 0: self.log.printLog('\r#SYS',command,log=False)
            #x#print command
            os.system(command)
            ## ~ [2c] ~ Cleanup database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if cleandb: cleanupDB(self,self.info['DBase'])
        except: self.errorLog('Fatal Error during BLASTRun.blast()'); raise
#########################################################################################################################
    def readBLAST(self,resfile=None,clear=False,gablam=False,unlink=False,local=False,screen=True,log=False):  ### Reads BLAST Results into objects
        '''
        Reads BLAST Results into objects.
        >> resfile:str = Results File (set as self.info['Name'])
        >> clear:Boolean = whether to clear current searches (True) or just append (False) [False]
        >> gablam:Boolean = whether to calculate gablam statistics and clear alignments to save memory [False]
        >> unlink:Boolean = whether to delete BLAST results file after reading [False]
        >> local:Boolean = whether to store Local alignment dictionary with basic alignment data [False]
        >> screen:Bool [False] = whether to output reading details to screen
        >> log:Bool [False] = whether to output reading details to log
        << returns True if (apparently) read to completion OK, else False
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if resfile != None: self.info['Name'] = resfile
            resline = self.loadFromFile(filename=self.info['Name'],v=1,checkpath=False,chomplines=True)
            try:
                newtype = string.split(resline[0])[0].lower()
                if newtype in ['blastn','blastp','blastx','tblastn','tblastx','rpsblast','rpstblastn']: self.info['Type'] = newtype
            except:
                if not resline: self.errorLog('No lines read in from %s.' % self.info['Name'],printerror=False)
                else: self.errorLog('Error with lines read in from %s.' % self.info['Name'],printerror=False)
                raise IOError
            if clear: self.search = []
            ### ~ [2] ~ Read in Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtxt = 'Reading %s %s BLAST results' % (self.info['Name'],self.info['Type']); (sx,hx) = (0,0)
            self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
            search = None; readhits = False
            i = 0
            hitaln = 0
            while i < len(resline):
                line = resline[i]
                ## ~ [2a] ~ Basic Search Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.find('Query=') == 0:    # New Search
                    search = self._addSearch(); readhits = False
                    search.info['Name'] = rje.matchExp('^Query=\s+(\S+)', line)[0]
                    sx += 1
                elif re.search('^\s+\(\S+ letters\)',line):
                    len_match = string.replace(rje.matchExp('^\s+\((\S+) letters\)',line)[0],',','')
                    search.stat['Length'] = string.atoi(len_match)
                elif line.find('Number of letters in database:') >= 0:
                    dblen = rje.matchExp('Number of letters in database:\s+(\d\S*)', line)[0]
                    dblen = re.sub('\D','',dblen)
                    self.stat['DBLen'] = string.atoi(dblen)
                elif line.find('Number of sequences in database:') >= 0:
                    dbnum = rje.matchExp('Number of sequences in database:\s+(\d\S*)', line)[0]
                    dbnum = re.sub('\D','',dbnum)
                    self.stat['DBNum'] = string.atoi(dbnum)
                elif rje.matchExp('(\S+) sequences; (\S+) total letters',line):
                    (dbnum,dblen) = rje.matchExp('(\S+) sequences; (\S+) total letters',line)
                    self.stat['DBNum'] = string.atoi(re.sub('\D','',dbnum))
                    self.stat['DBLen'] = string.atoi(re.sub('\D','',dblen))
                elif line.find('Number of sequences better than') >= 0:
                    self.stat['E-Value'] = string.atof(rje.matchExp('Number of sequences better than\s+(\S+):', line)[0])
                ## ~ [2b] ~ One-line hit data (BLASTHit) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if string.split(line)[-1] == 'N': (si,ei) = (-3,-2)
                    else: (si,ei) = (-2,-1)
                    i += 2  # Skip blank line
                    while rje.matchExp('^(\S+)\s.*\s(\S*\d)\s+(\S*\d)\s*$',resline[i]):
                        #match = rje.matchExp('^(\S+)\s.*\s(\d\S*)\s+(\S+\d)\s*$',resline[i])
                        match = string.split(resline[i])
                        hit = search._addHit()
                        hit.setInfo({'Name':match[0],'Type':self.info['Type']})
                        hit.stat['BitScore'] = string.atof(match[si])
                        eval = match[ei]
                        if eval.find('e') == 0: eval = '1' + eval
                        hit.stat['E-Value'] = string.atof(eval)
                        i += 1
                        hx += 1
                    line = resline[i]   # End of one-lines (blank line)
                    self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
                    search.checkHitNames()
                    hitaln = 0; readhits = True
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    search.hit = []
                    self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
                    hitaln = 0; readhits = True
                ## ~ [2c] ~ Aln Hit data (PWAln) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('>') == 0:
                    if not readhits: i += 1; continue   # '>' Character within name happened to occur at start of line
                    hitname = rje.matchExp('^>(\S+)\s*',line)[0]
                    if hitaln >= len(search.hit):
                        self.errorLog('Apparent hits for %s have exceeded the %d found. (Hit %s.) BLAST read-through or sequence name format error?' % (search.info['Name'],hitname,len(search.hit)),False,False)
                        raise ValueError
                    #if hitaln not in search.hit:
                    #    self.errorLog('Cannot find hitaln %s for %s hits. Sequence name format error?' % (hitaln,search.info['Name']))
                    #    raise ValueError
                    if hitname != search.hit[hitaln].info['Name']:      # Identify hit object
                        for hit in search.hit:
                            if hit.info['Name'] == hitname: hitaln = search.hit.index(hit)
                        if hitname != search.hit[hitaln].info['Name']:
                            self.errorLog('Problem with BLAST results - %s single-line hits and alignments do not match' % search.info['Name'],printerror=False,quitchoice=True)
                            i += 1; continue
                    hit = search.hit[hitaln]
                    hitaln += 1
                    aln = None
                    while i < (len(resline) - 1):
                        i += 1
                        line = resline[i]
                        ## Hit Length ##
                        if re.search('^\s+Length\s+=\s+(\d+)',line):
                            hit.stat['Length'] = string.atoi(rje.matchExp('^\s+Length = (\d+)',line)[0])
                        ## New Aln Block ##
                        elif hit.stat['Length'] and rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line):
                            aln = hit._addAln()
                            scores = rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line)
                            aln.stat['BitScore'] = string.atof(scores[0])
                            eval = scores[1]
                            if eval.find('e') == 0: eval = '1' + eval
                            aln.stat['Expect'] = string.atof(eval)
                            i += 1
                            if re.search('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives = (\d+)/(\d+)\s?',resline[i]):
                                sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives\s+=\s+(\d+)/(\d+)\s?',resline[i])
                                aln.stat['Positives'] = string.atoi(sim[2])
                            else: sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?',resline[i])
                            try:
                                aln.stat['Length'] = string.atoi(sim[1])
                                aln.stat['Identity'] = string.atoi(sim[0])
                            except: self.errorLog('Problem reading line "%s"' % resline[i])
                            if resline[i+1].find('Frame') >= 0: i += 1    # BLASTX or TBLASTN
                        ## Alignment ##
                        elif line.find('Query:') == 0 and aln != None:
                            # Query Line
                            alnstructure = rje.matchExp('^(Query:\s+\d+\s+)(\S+)\s+(\d+)',line)
                            if not alnstructure:    # Something's wrong!
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                self.deBug(line)
                                aln._clear()
                                break
                            leader = len(alnstructure[0])
                            aln.info['QrySeq'] += alnstructure[1]
                            aln.stat['QryEnd'] = string.atoi(alnstructure[2])
                            if aln.stat['QryStart'] == 0:
                                aln.stat['QryStart'] = string.atoi(rje.matchExp('^Query:\s+(\d+)\s',alnstructure[0])[0])
                            # Alignment Line
                            i += 1; aln.info['AlnSeq'] += resline[i][leader:(leader+len(alnstructure[1]))]
                            # Subject Line
                            i += 1; subject = rje.matchExp('^Sbjct:\s+(\d+)\s*(\D+)\s+(\d+)',resline[i])
                            if not subject:
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                self.deBug(resline[i])
                                aln._clear()
                                break
                            aln.info['SbjSeq'] += subject[1]
                            aln.stat['SbjEnd'] = string.atoi(subject[2])
                            if aln.stat['SbjStart'] == 0: aln.stat['SbjStart'] = string.atoi(subject[0])
                        ## End of Aln ##
                        elif line.find('>') == 0:
                            i -= 1  # i Advanced anyway, so subtract 1
                            break   # so '>' will still be caught as next Aln
                        ## End of All Alns ##
                        elif line.find('Database:') >= 0 and hit.alnNum() > 0: break   # hit.alnNum() to avoid matching 'Database:' in hit description
                        elif line.find('BLAST') == 0: break
                        elif line.find('TBLASTN') == 0: break
                ## ~ [2d] ~ Continue reading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##      
                i += 1
            ptxt = '%s complete: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx))
            self.printLog('\r#BLAST','%s (vs %s Sequences of total %s letters)' % (ptxt,rje.integerString(self.stat['DBNum']),rje.integerString(self.stat['DBLen'])),log=log,screen=screen)
            ### ~ [3] ~ Local alignment options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if local:
                self.progLog('#LOCAL','Generating local aln dictionaries...',screen=screen)
                for search in self.search:
                    for hit in search.hit: hit.makeLocalDict()
                self.progLog('\r#LOCAL','Generation of local aln dictionaries complete!',screen=screen)
            ### ~ [4] ~ GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gablam:  
                logtxt = 'Calculating GABLAM statistics from %s BLAST results' % self.info['Name']
                _hitx = 0
                sx = 0.0
                for search in self.search:
                    self.progLog('\r#GABLAM','%s: %.1f%%' % (logtxt,sx/len(self.search)),screen=screen)
                    _hitx += search.hitNum()
                    search.gablam()
                    sx += 100.0
                logtxt = 'Calculation of GABLAM statistics from %s BLAST results (%s hits) complete.' % (self.info['Name'],rje.integerString(_hitx))
                self.printLog('\r#GABLAM',logtxt,log=log,screen=screen)
            ### ~ [5] ~ Unlink if necessary and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if unlink and os.path.exists(self.info['Name']): os.unlink(self.info['Name'])
            return True
        except:
            if self.checkBLAST():
                self.errorLog('Fatal Error during BLASTRun.readBLAST() despite intact BLAST results.')
                raise
            else:
                self.errorLog('BLASTRun.readBLAST() failed due to incomplete BLAST results file.')
                return False
#########################################################################################################################
    def readNextBLASTSearch(self,fpos=0,resfile=None,gablam=False,local=False):  ### Reads BLAST Result into objects
        '''
        Reads BLAST Results into objects.
        >> resfile:str = Results File 
        >> gablam:Boolean = whether to calculate gablam statistics and clear alignments to save memory [False]
        >> local:Boolean = whether to store Local alignment dictionary with basic alignment data [False]
        >> fpos:Integer = position in results file to start reading from.
        << returns tuple: (Search object or None if no more results to be read,fpos).
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fpos < 0: return (None,-1)
            if not resfile: resfile = self.info['Name']
            RESFILE = open(resfile,'r'); RESFILE.seek(fpos)
            ### ~ [2] ~ Read in Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            line = RESFILE.readline()
            search = None; readhits = False; hitaln = 0; i = 0
            while line:
                while i > 0 and line: fpos = RESFILE.tell(); line = RESFILE.readline(); i -= 1
                if not line: break
                ## ~ [2a] ~ Basic Search Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.find('Query=') == 0:    # New Search
                    if search: break            # Read one whole search, so stop and return
                    search = BLASTSearch(log=self.log); readhits = False
                    search.info['Name'] = rje.matchExp('^Query=\s+(\S+)', line)[0]
                elif re.search('^\s+\(\S+ letters\)',line):
                    len_match = string.replace(rje.matchExp('^\s+\((\S+) letters\)',line)[0],',','')
                    search.stat['Length'] = string.atoi(len_match)
                elif line.find('Number of letters in database:') >= 0:
                    dblen = rje.matchExp('Number of letters in database:\s+(\d\S*)', line)[0]
                    dblen = re.sub('\D','',dblen)
                    self.stat['DBLen'] = string.atoi(dblen)
                elif line.find('Number of sequences in database:') >= 0:
                    dbnum = rje.matchExp('Number of sequences in database:\s+(\d\S*)', line)[0]
                    dbnum = re.sub('\D','',dbnum)
                    self.stat['DBNum'] = string.atoi(dbnum)
                elif rje.matchExp('(\S+) sequences; (\S+) total letters',line):
                    (dbnum,dblen) = rje.matchExp('(\S+) sequences; (\S+) total letters',line)
                    self.stat['DBNum'] = string.atoi(re.sub('\D','',dbnum))
                    self.stat['DBLen'] = string.atoi(re.sub('\D','',dblen))
                elif line.find('Number of sequences better than') >= 0:
                    self.stat['E-Value'] = string.atof(rje.matchExp('Number of sequences better than\s+(\S+):', line)[0])
                ## ~ [2b] ~ One-line hit data (BLASTHit) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if string.split(line)[-1] == 'N': (si,ei) = (-3,-2)
                    else: (si,ei) = (-2,-1)
                    RESFILE.readline(); fpos = RESFILE.tell(); line = RESFILE.readline()   # Skip blank line
                    while rje.matchExp('^(\S+)\s.*\s(\S*\d)\s+(\S*\d)\s*$',line):
                        match = string.split(line)
                        hit = search._addHit()
                        hit.setInfo({'Name':match[0],'Type':self.info['Type']})
                        hit.stat['BitScore'] = string.atof(match[si])
                        eval = match[ei]
                        if eval.find('e') == 0: eval = '1' + eval
                        hit.stat['E-Value'] = string.atof(eval)
                        fpos = RESFILE.tell(); line = RESFILE.readline()
                    search.checkHitNames()
                    hitaln = 0; readhits = True
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    search.hit = []
                    hitaln = 0; readhits = True
                ## ~ [2c] ~ Aln Hit data (PWAln) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('>') == 0:
                    if not readhits: i += 1; continue   # '>' Character within name happened to occur at start of line
                    hitname = rje.matchExp('^>(\S+)\s*',line)[0]
                    if hitaln >= len(search.hit):
                        self.errorLog('Apparent hits for %s have exceeded the %d found. (Hit %s.) BLAST read-through or sequence name format error?' % (search.info['Name'],hitname,len(search.hit)),False,False)
                        raise ValueError
                    if hitname != search.hit[hitaln].info['Name']:      # Identify hit object
                        for hit in search.hit:
                            if hit.info['Name'] == hitname: hitaln = search.hit.index(hit)
                        if hitname != search.hit[hitaln].info['Name']:
                            self.errorLog('Problem with BLAST results - %s single-line hits and alignments do not match' % search.info['Name'],printerror=False,quitchoice=True)
                            i += 1; continue
                    hit = search.hit[hitaln]
                    hitaln += 1
                    aln = None
                    while line:
                        fpos = RESFILE.tell(); line = RESFILE.readline()
                        if not line: break
                        ## Hit Length ##
                        if re.search('^\s+Length\s+=\s+(\d+)',line):
                            hit.stat['Length'] = string.atoi(rje.matchExp('^\s+Length = (\d+)',line)[0])
                        ## New Aln Block ##
                        elif hit.stat['Length'] and rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line):
                            aln = hit._addAln()
                            scores = rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line)
                            aln.stat['BitScore'] = string.atof(scores[0])
                            eval = scores[1]
                            if eval.find('e') == 0: eval = '1' + eval
                            aln.stat['Expect'] = string.atof(eval)
                            fpos = RESFILE.tell(); line = RESFILE.readline()
                            if re.search('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives = (\d+)/(\d+)\s?',line):
                                sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives\s+=\s+(\d+)/(\d+)\s?',line)
                                aln.stat['Positives'] = string.atoi(sim[2])
                            else: sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?',line)
                            try:
                                aln.stat['Length'] = string.atoi(sim[1])
                                aln.stat['Identity'] = string.atoi(sim[0])
                            except: self.errorLog('Problem reading line "%s"' % line)
                            fpos = RESFILE.tell(); line = RESFILE.readline()
                            if line.find('Frame') < 0: RESFILE.seek(fpos)    # Not BLASTX or TBLASTN
                        ## Alignment ##
                        elif line.find('Query:') == 0 and aln != None:
                            # Query Line
                            alnstructure = rje.matchExp('^(Query:\s+\d+\s+)(\S+)\s+(\d+)',line)
                            if not alnstructure:    # Something's wrong!
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                aln._clear()
                                break
                            leader = len(alnstructure[0])
                            aln.info['QrySeq'] += alnstructure[1]
                            aln.stat['QryEnd'] = string.atoi(alnstructure[2])
                            if aln.stat['QryStart'] == 0:
                                aln.stat['QryStart'] = string.atoi(rje.matchExp('^Query:\s+(\d+)\s',alnstructure[0])[0])
                            # Alignment Line
                            fpos = RESFILE.tell(); line = RESFILE.readline(); aln.info['AlnSeq'] += line[leader:(leader+len(alnstructure[1]))]
                            # Subject Line
                            fpos = RESFILE.tell(); line = RESFILE.readline(); subject = rje.matchExp('^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)',line)
                            if not subject:
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                aln._clear()
                                break
                            aln.info['SbjSeq'] += subject[1]
                            aln.stat['SbjEnd'] = string.atoi(subject[2])
                            if aln.stat['SbjStart'] == 0: aln.stat['SbjStart'] = string.atoi(subject[0])
                        ## End of Aln ##
                        elif line.find('>') == 0:
                            RESFILE.seek(fpos)  # Return to beginning of this line
                            break   # so '>' will still be caught as next Aln
                        ## End of All Alns ##
                        elif line.find('Database:') >= 0 and hit.alnNum() > 0: break   # hit.alnNum() to avoid matching 'Database:' in hit description
                        elif line.find('BLAST') == 0: break
                        elif line.find('TBLASTN') == 0: break
                ## ~ [2d] ~ Continue reading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##      
                i += 1
            RESFILE.close()
            ### ~ [3] ~ Local alignment options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if local and search:
                for hit in search.hit: hit.makeLocalDict()
            ### ~ [4] ~ GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gablam and search: search.gablam()
            return (search,fpos)
        except:
            try: RESFILE.close()
            except: pass
            self.errorLog('Error during BLASTRun.readNextBLASTSearch()'); raise
#########################################################################################################################
    def checkBLAST(self,resfile=None,logcheck=True):  ### Checks that each BLAST started has an end
        '''
        Checks that each BLAST started has an end, thus identifying BLAST runs that have been terminated prematurely and
        need to be re-run.
        >> resfile:str = Results File (set as self.info['Name'])
        >> logcheck:boolean = Whether to print findings to log [True]
        << True if OK, False if not.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if resfile != None: self.info['Name'] = resfile
            if not os.path.exists(self.info['Name']):
                if logcheck: self.errorLog('BLAST results file "%s" missing!' % self.info['Name'],printerror=False)
                return False
            RESFILE = open(self.info['Name'],'r');
            RESFILE.seek(0,2); fend = RESFILE.tell(); RESFILE.seek(0)
            line = RESFILE.readline()
            if not line:
                if logcheck: self.errorLog('No lines read in from %s.' % self.info['Name'],printerror=False)
                return False
            ### ~ [2] ~ Check Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtxt = 'Checking %s BLAST results' % self.info['Name']
            px = 0.0; self.progLog('\r#BLAST','%s: %.2f%%' % (logtxt,px/fend))
            last = 'Top'    # Need at least one search to be intact
            while line:
                self.progLog('\r#BLAST','%s: %.2f%%' % (logtxt,px/fend))
                if line.find('Query=') == 0 and last in [None,'Top']:
                    last = line   # New Search
                    px = 100.0 * RESFILE.tell(); self.progLog('\r#BLAST','%s: %.2f%%' % (logtxt,px/fend))
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if last.find('Query=') != 0:    # Wrong place!
                        if logcheck: self.errorLog('%s check: Found "Query=" in wrong place!' % self.info['Name'],printerror=False)
                        return False
                    last = line
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    if last.find('Query=') != 0:    # Wrong place!
                        if logcheck: self.errorLog('%s check: Found "Query=" in wrong place!' % self.info['Name'],printerror=False)
                        return False
                    last = line
                elif line.find('Effective length') >= 0 or line.find('Gap Penalties:') >= 0 or line[:5] in ['BLAST','Matri']:
                    if not last: pass
                    elif last.find('Sequences producing significant alignments:') >= 0: last = None     # One-line hits
                    elif last.find('***** No hits found ******') >= 0: last = None                      # No Hits
                line = RESFILE.readline()
            ### ~ [3] ~ Finish check and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            RESFILE.close()
            if last: self.printLog('\r#BLAST','%s: BLAST incomplete?!' % logtxt,log=logcheck); return False
            else: self.printLog('\r#BLAST','%s: BLAST intact.' % logtxt,log=logcheck); return True
        except:
            try: RESFILE.close()
            except: pass
            self.errorLog('Error during BLASTRun.checkBLAST().'); return False
#########################################################################################################################
    def _addSearch(self):   ### Adds and returns a new search object
        '''Adds and returns a new search object.'''
        newsearch = BLASTSearch(log=self.log,cmd_list=self.cmd_list)#; newsearch.setOpt({'DeBug':self.getBool('DeBug')})
        self.search.append(newsearch)
        return newsearch
#########################################################################################################################
    ### <3> ### Output                                                                                                  #
#########################################################################################################################
    def saveCutBLAST(self,outfile=None):   ### Saves a cutdown version of current BLAST (no alignments)
        '''
        Saves a cutdown version of current BLAST (no alignments), which has enough lines to be successfully read in by
        the BLASTRun class. (Obviously, no GABLAM statistics can be calculated!
        >> outfile:str = outfile to use if different from self.info['Name']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfile: self.info['Name'] = outfile
            CUT = open(self.info['Name'], 'w')
            ### ~ [2] ~ Output cut-down search data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for search in self.search:
                CUT.write('Query= %s\n\n' % search.info['Name'])
                CUT.write('         (%d letters)\n' % search.stat['Length'])
                CUT.write('Number of letters in database: %d\n' % self.stat['DBLen'])
                CUT.write('Number of sequences in database: %d\n' % self.stat['DBNum'])
                CUT.write('Number of sequences better than %e:\n' % self.stat['E-Value'])
                CUT.write('\nSequences producing significant alignments:\n\n')
                if search.hitNum() > 0:
                    for hit in search.hit:
                        CUT.write('%s  %.1f  %e\n' % (hit.info['Name'],hit.stat['BitScore'],hit.stat['E-Value']))
                else: CUT.write('***** No hits found ******\n')
                CUT.write('\n\n')
            CUT.close()
        except: self.errorLog('Major error during BLASTRun.saveCutBLAST().',quitchoice=True)
#########################################################################################################################
    def hitToSeq(self,seqlist,searchlist=[],filename=None,appendfile=False):   ### Saves hits from given searches to sequence object/file
        '''
        Saves hits from given searches to sequence object and, if given, a file.
        >> seqlist:rje_seq.SeqList Object *Necessary!* 
        >> seachlist:list of BLASTSearch objects [self.search if none]
        >> filename:str = Name of fasta output file - no save if None [None]
        >> appendfile:bool = Whether to append file
        << returns dictionary of {Hit:Sequence}
        '''
        try:### ~ [1] ~ Map hits to seqlist sequences and return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitseq = {}
            if searchlist == []: searchlist = self.search[0:]
            for search in searchlist:            
                for hit in search.hit:
                    newseq = seqlist.seqFromFastaCmd(hit.info['Name'],self.info['DBase'])
                    if newseq:
                        hitseq[hit] = newseq
                        self.verbose(2,4,'Seq %d retrieved from %s: %s.' % (seqlist.seqNum(),self.info['DBase'],newseq.shortName()),1)
                    else: self.errorLog('No Seq retrieval for %s using BLASTRun.hitToSeq().' % hit.info['Name'],printerror=False)
            if hitseq and filename and filename.lower() != 'none': seqlist.saveFasta(seqfile=filename,append=appendfile)
            return hitseq
        except: self.errorLog('Major error during BLASTRun.hitToSeq().');  raise
#########################################################################################################################
    def searchSeq(self,seqlist,proglog=True,inverse=False):   ### Returns dictionary of searches as sequences from seqlist (or None if missing)
        '''
        Returns dictionary of searchesas sequences from seqlist (or None if missing).
        >> seqlist:SeqList object
        >> proglog:bool [True] = whether to log progress 
        >> inverse:bool [False] = whether to reverse dictionary to sequence:search
        << {Search Object: Sequence Object}
        '''
        try:### ~ [1] ~ Setup (Clear objects) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            searchdic = {}
            for search in self.search: searchdic[search] = None
            ### ~ [2] ~ Map SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist:
                seqdic = seqlist.seqNameDic(proglog=False)  #proglog)
                (mx,hx) = (0,0.0)
                for search in self.search:
                    hx += 100.0
                    if proglog:
                        ltxt = 'Mapping %s searches onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
                        self.printLog('\r#SEARCH','%s: %.2f%%' % (ltxt,(hx/self.searchNum())),newline=False,log=False)
                    if seqdic.has_key(search.info['Name']):
                        searchdic[search] = seqdic[search.info['Name']]
                        mx += 1
                    else:
                        for key in seqdic.keys():
                            if key.find('|%s' % search.info['Name']) > 0:
                                searchdic[search] = seqdic[key]
                                mx += 1
                                break   #!# Add report_none=False?
                if proglog:
                    ltxt = 'Mapping %s searches onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
                    self.printLog('\r#SEARCH','%s: %s of %s mapped.' % (ltxt,rje.integerString(mx),rje.integerString(self.searchNum())))
            ### ~ [3] ~ Finish and return (empty if no seqlist) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if inverse:
                invdic = {}
                for search in searchdic:
                    if searchdic[search]: invdic[searchdic[search]] = search
                return invdic
            return searchdic                        
        except: self.errorLog('Major Problem with BLASTSearch.searchSeq().'); raise       
#########################################################################################################################
### END OF SECTION II: BLASTRun Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: BLASTSearch Class                                                                                      #
#########################################################################################################################
class BLASTSearch(rje.RJE_Object):     
    '''
    BLAST Search Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Short Name of Search Query (First word of description line)
    
    Opt:boolean

    Stat:numeric
    - Length = Sequence Length

    List:list

    Dict:dictionary    

    Obj:RJE_Objects

    Other:
    hit : list of BLASTHit Objects
    '''
    ### Attributes
    hit = []    #!# Also change rje_hmm if changing this #!#
    def hitNum(self): return len(self.hit)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Name']
        self.statlist = ['Length']
        self.optlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.hit = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try: self._generalCmd(cmd)   ### General Options ###
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def _addHit(self):   ### Adds and returns a new search object
        '''Adds and returns a new search object.'''
        newhit = BLASTHit(log=self.log,cmd_list=self.cmd_list)#; newhit.setOpt({'DeBug':self.getBool('DeBug')})
        newhit._cmdList()
        self.hit.append(newhit)
        return newhit
#########################################################################################################################
    def checkHitNames(self): ### Checks for multiple hits of same name
        '''Checks for multiple hits of same name.'''
        try:### ~ [1] ~ Make Name List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            names = []
            for hit in self.hit: names.append(hit.info['Name'])                
            ### ~ [2] ~ Dictionary of multiple occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            multiples = {}  # Dictionary of multiple occurrences
            for name in names:
                if names.count(name) > 1: multiples[name] = names.count(name)
            ### ~ [3] ~ Warnings! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for name in rje.sortKeys(multiples): self.errorLog('Sequence "%s" hit %d times! Check that sequence names are unique in database.' % (name,multiples[name]),printerror=False)
        except: self.errorLog('Problem with %s checkHitNames()' % self.info['Name'])
#########################################################################################################################
    def hitSeq(self,seqlist,proglog=True,inverse=False):   ### Returns dictionary of hits as sequences from seqlist (or None if missing)
        '''
        Returns dictionary of hits as sequences from seqlist (or None if missing).
        >> seqlist:SeqList object
        >> proglog:bool [True] = whether to log progress 
        >> inverse:bool [False] = whether to reverse dictionary to sequence:hit
        << {Hit Object: Sequence Object}
        '''
        try:### ~ [1] ~ Setup (Clear) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitdic = {}
            for hit in self.hit: hitdic[hit] = None
            ltxt = 'Mapping %s hits onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
            ### ~ [2] ~ Map SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist:
                seqdic = seqlist.seqNameDic(proglog=proglog)
                (mx,hx) = (0,0.0)
                for hit in self.hit:
                    hx += 100.0
                    if proglog: self.progLog('\r#HIT','%s: %.2f%%' % (ltxt,(hx/self.hitNum())))
                    if seqdic.has_key(hit.info['Name']):
                        hitdic[hit] = seqdic[hit.info['Name']]
                        mx += 1
                    else:
                        for key in seqdic.keys():
                            if key.find('|%s' % hit.info['Name']) > 0:
                                hitdic[hit] = seqdic[key]
                                mx += 1
                                break   #!# Add report_none=False?
                if proglog: self.printLog('\r#HIT','%s: %s of %s mapped.' % (ltxt,rje.integerString(mx),rje.integerString(self.hitNum())))
            ### Return (empty if no seqlist) ###
            if inverse:
                invdic = {}
                for hit in hitdic:
                    if hitdic[hit]: invdic[hitdic[hit]] = hit
                    return invdic
            return hitdic                        
        except: self.errorLog('Major Problem with BLASTSearch.hitSeq().'); raise       
#########################################################################################################################
    def saveHitIDs(self,outfile=None):   ### Saves IDs of hits in file (e.g. for later fastacmd extraction)
        '''
        Saves IDs of hits in file (e.g. for later fastacmd extraction).
        >> outfile:str = outfile to use.
        '''
        try:### ~ [1] ~ Save simple list of IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.hitNum() < 1: return
            if outfile == None: outfile = '%s.id' % self.info['Name'] 
            IDS = open(outfile, 'w')
            for hit in self.hit: IDS.write('%s\n' % hit.info['Name'])
            IDS.close()
        except: self.errorLog('Major error during BLASTSearch.saveHitIDs().')
#########################################################################################################################
    def gablam(self,keepaln=False):   ### Performs GABLAM analysis on all hits (unless done) and clears alignments
        '''
        Performs GABLAM analysis on all hits (unless done) and clears alignments.
        >> keepaln:bool [False] = Whether to store GABLAM alignments in Hit object
        '''
        try:### ~ [1] ~ Perform GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hit in self.hit:
                if not hit.opt['GABLAM']:
                    hit.globalFromLocal(self.stat['Length'],keepaln)
                    for aln in hit.aln:
                        aln.info['QrySeq'] = ''     #!# Saves memory #!#
                        aln.info['SbjSeq'] = ''     #!# Saves memory #!#
                        aln.info['AlnSeq'] = ''     #!# Saves memory #!#
        except: self.log.errorLog('Major error with BLASTSearch.gablam().')
#########################################################################################################################
### END OF SECTION III: BLASTSearch Class                                                                               #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: BLASTHit Class                                                                                          # 
#########################################################################################################################
class BLASTHit(rje.RJE_Object):     
    '''
    BLAST Hit Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Short Name of Hit Sequence (First Word of description line)
    
    Opt:boolean
    - GABLAM = Whether GABLAM search has been performed

    Stat:numeric
    - BitScore
    - E-Value
    - Length
    - GablamFrag = Length of gaps between mapped residue for fragmenting local hits [100]
    - LocalCut = Cut-off length for local alignments contributing to global GABLAM stats) [0]

    List:list

    Dict:Dictionary
    - GABLAM = GABLAM and GABLAMO stats
    - Local = Local alignment dictionary {alnID:{stats}}

    Obj:RJE_Objects

    Other:
    - aln:list of PWAln Objects
    '''
    ### Attributes
    aln = []    # List of PWAln Objects
    def alnNum(self): return len(self.aln)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Name','Type','Description']
        self.statlist = ['BitScore','E-Value','Length','GablamFrag','LocalCut']
        self.optlist = ['GABLAM']
        self.listlist = []
        self.objlist = []
        self.dictlist = ['GABLAM','Local']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setdict=True)
        self.setInt({'GablamFrag':0,'LocalCut':0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.aln = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'int',['GablamFrag','LocalCut'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def _addAln(self):  ### Adds and returns a PWAln Object
        '''Adds and returns a PWAln Object.'''
        newaln = PWAln(log=self.log)
        self.aln.append(newaln)
        return newaln
#########################################################################################################################
    def makeLocalDict(self):    ### Generates summary Local alignment dictionary
        '''Generates summary Local alignment dictionary.'''
        try:### ~ [1] Compile alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(self.aln)):
                aln = self.aln[i]
                self.dict['Local'][i+1] = rje.combineDict({},aln.stat) 
        except: self.errorLog('Major problem during BLASTHit.makeLocalDict()')
#########################################################################################################################
    def globalFromLocal(self,qrylen,keepaln=False):  ### Returns a dictionary of global alignment stats from Query-Hit local alignments
        '''
        Returns a dictionary of global alignment stats from Query-Hit local alignments.
        >> querylen:int = length of query sequence
        >> keepaln:bool [False] = Whether to store GABLAM alignments in Hit object
        << {Query:{},Hit:{}}, where each value is a dictionary of [GABLAM(O) ID, GABLAM(O) Sim, GABLAM(O) Len]
        '''
        try:
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['GABLAM']: return self.dict['GABLAM']       # Already performed: do not repeat on empty aln
            gablam = {} # Global Alignment from BLAST Local Alignment Matrix
            for aln in ['Qry','QryO']: gablam[aln] = ['-'] * int(qrylen)                # O = ordered
            for aln in ['Hit','HitO']: gablam[aln] = ['-'] * int(self.stat['Length'])   # O = ordered
            for qh in ['Qry','Hit']:
                gablam['%sDirn' % qh] = 'None'
                gablam['%sStart' % qh] = len(gablam[qh])
                gablam['%sEnd' % qh] = 0
                gablam['%sDirnO' % qh] = 'None'
                gablam['%sStartO' % qh] = len(gablam[qh])
                gablam['%sEndO' % qh] = 0
            lpairs = []     # List of paired query-hit residue tuples for ordered GABLAM
            orientation = {'None':{True:'Bwd',False:'Fwd'},'Fwd':{True:'Both',False:'Fwd'},
                           'Bwd':{True:'Bwd',False:'Both'},'Both':{True:'Both',False:'Both'}}

            ### ~ [2] Compile alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug('GABLAMO for Hit %s: %d alignments' % (self.info['Name'],self.alnNum()))
            lcutx = 0
            for aln in self.aln:
                if aln.stat['Length'] < self.stat['LocalCut']: lcutx += 1; continue
                ## ~ [2a] Assess for backwards hit as in nucleotide BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qbackwards = aln.stat['QryEnd'] < aln.stat['QryStart']       # Whether match reversed (e.g. BLASTX)
                sbackwards = aln.stat['SbjEnd'] < aln.stat['SbjStart']       # Whether match reversed (e.g. TBLASTN)
                gablam['QryDirn'] = orientation[gablam['QryDirn']][qbackwards]
                gablam['HitDirn'] = orientation[gablam['HitDirn']][sbackwards]
                gablam['QryStart'] = min(gablam['QryStart'],aln.stat['QryStart'],aln.stat['QryEnd'])
                gablam['HitStart'] = min(gablam['HitStart'],aln.stat['SbjStart'],aln.stat['SbjEnd'])
                gablam['QryEnd'] = max(gablam['QryEnd'],aln.stat['QryStart'],aln.stat['QryEnd'])
                gablam['HitEnd'] = max(gablam['HitEnd'],aln.stat['SbjStart'],aln.stat['SbjEnd'])
                #x#    self.log.printLog('#REV','Query has reverse orientation for alignment %s: Ignoring for GABLAM.' % self.aln.index(aln))
                #x#    continue
                ## ~ [2b] Assess aln for GABLAMO ordering with previous alns ~~~~~~~~~~~~~~~~~~~~~~ ##
                order_ok = True #x# not (qbackwards or sbackwards)      # GABLAMO is in forward direction only
                #!# Update method to allow GABLAMO in all directions #!#
                for pair in lpairs:
                    if qbackwards == sbackwards:
                        if aln.stat['QryStart'] < pair[0] and aln.stat['SbjStart'] > pair[1]: order_ok = False
                        elif aln.stat['QryStart'] > pair[0] and aln.stat['SbjStart'] < pair[1]: order_ok = False
                        elif aln.stat['QryEnd'] < pair[0] and aln.stat['SbjEnd'] > pair[1]: order_ok = False
                        elif aln.stat['QryEnd'] > pair[0] and aln.stat['SbjEnd'] < pair[1]: order_ok = False
                    else:
                        if aln.stat['QryStart'] < pair[0] and aln.stat['SbjStart'] < pair[1]: order_ok = False
                        elif aln.stat['QryStart'] > pair[0] and aln.stat['SbjStart'] > pair[1]: order_ok = False
                        elif aln.stat['QryEnd'] < pair[0] and aln.stat['SbjEnd'] < pair[1]: order_ok = False
                        elif aln.stat['QryEnd'] > pair[0] and aln.stat['SbjEnd'] > pair[1]: order_ok = False
                if order_ok:
                    lpairs.append((aln.stat['QryStart'],aln.stat['SbjStart']))
                    lpairs.append((aln.stat['QryEnd'],aln.stat['SbjEnd']))
                    gablam['QryDirnO'] = orientation[gablam['QryDirnO']][qbackwards]
                    gablam['HitDirnO'] = orientation[gablam['HitDirnO']][sbackwards]
                    gablam['QryStartO'] = min(gablam['QryStartO'],aln.stat['QryStart'],aln.stat['QryEnd'])
                    gablam['HitStartO'] = min(gablam['HitStartO'],aln.stat['SbjStart'],aln.stat['SbjEnd'])
                    gablam['QryEndO'] = max(gablam['QryEndO'],aln.stat['QryStart'],aln.stat['QryEnd'])
                    gablam['HitEndO'] = max(gablam['HitEndO'],aln.stat['SbjStart'],aln.stat['SbjEnd'])
                ## ~ [2c] Perform GABLAM calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qres = aln.stat['QryStart'] - 1     # Count from zero in list
                hres = aln.stat['SbjStart'] - 1
                shop = qhop = 1
                if qbackwards: qhop = -1
                if sbackwards: shop = -1
                iloop = 1       # Number of loops needed for each aln position
                if self.info['Type'] in ['blastx','tblastn']: iloop = 3
                try:
                    for r in range(aln.stat['Length']): # Sim
                        for i in range(iloop):     # Break to exit this loop, depending on nt/protein GABLAM
                            if aln.info['AlnSeq'][r] == '+':
                                if gablam['Qry'][qres] != '|': gablam['Qry'][qres] = '+'
                                if gablam['Hit'][hres] != '|': gablam['Hit'][hres] = '+'
                                if gablam['QryO'][qres] != '|' and order_ok: gablam['QryO'][qres] = '+'
                                if gablam['HitO'][hres] != '|' and order_ok: gablam['HitO'][hres] = '+'
                            elif aln.info['AlnSeq'][r] == '|' or re.search('[A-Za-z]',aln.info['AlnSeq'][r]):   # ID
                                gablam['Qry'][qres] = '|'
                                gablam['Hit'][hres] = '|'
                                if order_ok:
                                    gablam['QryO'][qres] = '|'
                                    gablam['HitO'][hres] = '|'
                            else:
                                if gablam['Qry'][qres] == '-': gablam['Qry'][qres] = 'X'
                                if order_ok and gablam['QryO'][qres] == '-': gablam['QryO'][qres] = 'X'
                                if gablam['Hit'][hres] == '-': gablam['Hit'][hres] = 'X'
                                if order_ok and gablam['HitO'][hres] == '-': gablam['HitO'][hres] = 'X'
                            if i < 2 and self.info['Type'] == 'blastx' and aln.info['QrySeq'][r] not in [' ','-']: qres += qhop    # DNA triplet
                            if i < 2 and self.info['Type'] == 'tblastn' and aln.info['SbjSeq'][r] not in [' ','-']: hres += shop    # DNA triplet
                        if re.search('[A-Za-z]',aln.info['QrySeq'][r]) or aln.info['QrySeq'][r] == '*': qres += qhop
                        if re.search('[A-Za-z]',aln.info['SbjSeq'][r]) or aln.info['SbjSeq'][r] == '*': hres += shop
                    if qbackwards: qres += 2
                    if sbackwards: hres += 2
                    if qres != aln.stat['QryEnd']:
                        #!# Could be DNA sequence #!#
                        dres = ((qres - aln.stat['QryStart']) * 3) + aln.stat['QryStart'] + 2
                        if dres != aln.stat['QryEnd']:
                            self.log.errorLog('Hit %s: Query end position should be %d but reached %d (or %d) in Aln process!' % (self.info['Name'],aln.stat['QryEnd'],qres,dres),False,False)
                            print(aln.info)
                            print(aln.stat)
                            print(gablam)
                            raw_input('Continue?')
                            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                                gablam[aln] = ['X'] * qrylen    #!# Hit!! #!#
                            break
                    if hres != aln.stat['SbjEnd']:
                        ## Check backwards match ##
                        bwd_hres = hres - aln.stat['SbjStart']
                        #!# Could be DNA sequence #!#  dres = ((qres - aln.stat['QryStart']) * 3) + aln.stat['QryStart'] + 2
                        if aln.stat['SbjStart'] != bwd_hres:
                            self.log.errorLog('Hit %s: Subject end position should be %d but reached %d in Aln process (Bwd: %d)!' % (self.info['Name'],aln.stat['SbjEnd'],hres,bwd_hres),False,False)
                            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                                gablam[aln] = ['X'] * qrylen    #!# Hit!! #!#
                            break
                        else: self.printLog('#GAB','Subject appears to have reverse orientation. GABLAMO Results will be wrong!')
                except: self.errorLog('GABLAM problems: GABLAM stats will be wrong for Hit %s' % self.info['Name'])

            ### ~ [3] Make GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.verbose(2,4,'\n*** Hit: %s ***' % self.info['Name'],1)
            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                gablam[aln] = string.join(gablam[aln],'')
                #self.verbose(2,4,'%s: %s' % (aln, gablam[aln]),1)
                #self.bugPrint('%s: %s' % (aln, gablam[aln]))
                #self.deBug('%s: Len = %d; Unaln = %d' % (aln,len(gablam[aln]),string.count(gablam[aln],'-')))
            gdict = { 'Query':{}, 'Hit':{} }
            gdict['Query']['GABLAM Len'] = qrylen - string.count(gablam['Qry'],'-')
            gdict['Query']['GABLAM ID'] = string.count(gablam['Qry'],'|')
            gdict['Query']['GABLAM Sim'] = gdict['Query']['GABLAM ID'] + string.count(gablam['Qry'],'+')
            gdict['Query']['GABLAMO Len'] = qrylen - string.count(gablam['QryO'],'-')
            gdict['Query']['GABLAMO ID'] = string.count(gablam['QryO'],'|')
            gdict['Query']['GABLAMO Sim'] = gdict['Query']['GABLAMO ID'] + string.count(gablam['QryO'],'+')
            gdict['Hit']['GABLAM Len'] = self.stat['Length'] - string.count(gablam['Hit'],'-')
            gdict['Hit']['GABLAM ID'] = string.count(gablam['Hit'],'|')
            gdict['Hit']['GABLAM Sim'] = gdict['Hit']['GABLAM ID'] + string.count(gablam['Hit'],'+')
            gdict['Hit']['GABLAMO Len'] = self.stat['Length'] - string.count(gablam['HitO'],'-')
            gdict['Hit']['GABLAMO ID'] = string.count(gablam['HitO'],'|')
            gdict['Hit']['GABLAMO Sim'] = gdict['Hit']['GABLAMO ID'] + string.count(gablam['HitO'],'+')
            for gscore in ['Start','End','Dirn']:
                gdict['Query']['GABLAM %s' % gscore] = gablam['Qry%s' % gscore]
                gdict['Query']['GABLAMO %s' % gscore] = gablam['Qry%sO' % gscore]
                gdict['Hit']['GABLAM %s' % gscore] = gablam['Hit%s' % gscore]
                gdict['Hit']['GABLAMO %s' % gscore] = gablam['Hit%sO' % gscore]
            for gab in ['','O']:
                gdict['Hit']['GABLAM%s Frag' % gab] = []
                if self.getInt('GablamFrag') < 1: break     # No GABLAM Fragging
                i = 0
                #self.deBug('HIT SEQUENCE: %s' % gablam['Hit%s' % gab])
                #self.bugPrint('GABLAM%s: %s' % (gab,min(max(-1,gablam['Hit%s' % gab].find('|',i)),max(-1,gablam['Hit'].find('+',i)))))
                #self.deBug('Hit%s: Len = %d; Unaln = %d' % (gab,len(gablam['Hit%s' % gab]),string.count(gablam['Hit%s' % gab],'-')))
                while i > -1 and i < len(gablam['Hit%s' % gab]):
                    self.progLog('\r#FRAG','GABLAM%s Fragging: %.2f%%' % (gab,i*100.0/len(gablam['Hit%s' % gab])))
                    i1 = max(-1,gablam['Hit%s' % gab].find('|',i))
                    i2 = max(-1,gablam['Hit%s' % gab].find('+',i))
                    i = max(i1,i2)
                    if i < 0: break
                    elif min(i1,i2) > -1: i = min(i1,i2)
                    fragend = fragstart = i
                    gx = 0
                    while gx < self.getInt('GablamFrag') and i < (len(gablam['Hit%s' % gab]) - 1):
                        self.progLog('\r#FRAK','GABLAM%s Fragging: %.2f%%' % (gab,i*100.0/len(gablam['Hit%s' % gab])))
                        i += 1
                        if gablam['Hit%s' % gab][i] == '-': gx += 1
                        else: gx = 0; fragend = i
                    gdict['Hit']['GABLAM%s Frag' % gab].append((fragstart,fragend))
                    if gablam['Hit%s' % gab][i] != '-': i += 1
                    #self.deBug('%s' % (gdict['Hit']['GABLAM%s Frag' % gab]))
                self.progLog('\r#FRAG','GABLAM%s Fragging done!     ' % gab)
            if keepaln: gdict['Aln'] = gablam
            self.dict['GABLAM'] = gdict
            self.opt['GABLAM'] = True
            if lcutx: self.printLog('#CUT','%d local alignments < length %d ignored for GABLAM (localcut=X)' % (lcutx,self.stat['LocalCut']))
            return gdict
        except:
            self.log.errorLog('Major problem during BLASTHit.globalFromLocal vs Hit %s (qylen=%d)' % (self.info['Name'],qrylen))
            gdict = { 'Query':{}, 'Hit':{} }
            for stat in ['GABLAM Len','GABLAM ID','GABLAM Sim','GABLAMO Len','GABLAMO ID','GABLAMO Sim']:
                gdict['Query'][stat] = -1
                gdict['Hit'][stat] = -1
            self.dict['GABLAM'] = gdict
            self.opt['GABLAM'] = True #!# Really?
            return gdict
#########################################################################################################################
### END OF SECTION IV: BLASTHit Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: PWAln Class                                                                                              # 
#########################################################################################################################
class PWAln(rje.RJE_Object):     
    '''
    Pairwise Aligned Sequence (Global/Local) Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name
    - QrySeq = Text of aligned portion (Query only)
    - SbjSeq = Text of aligned portion (Subject only)
    - AlnSeq = Text of aligned portion (Alignment line)
    
    Opt:boolean

    Stat:numeric
    - BitScore
    - Expect
    - Length (integer value)
    - Identity (integer value)
    - Positives (integer value)
    - QryStart (1-N)
    - QryEnd
    - SbjStart (1-N)
    - SbjEnd

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Name','QrySeq','SbjSeq','AlnSeq']
        self.statlist = ['BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd']
        self.optlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.info['Name'] = 'None'
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try: self._generalCmd(cmd)   ### General Options ###
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def _clear(self):   ### Clears all data - dodgy alignment
        '''Clears all data - dodgy alignment.'''
        self.setStat({'BitScore':0.0,'Expect':1000,'Length':1,'Identity':0,'Positives':0,'QryStart':1,'QryEnd':1,'SbjStart':1,'SbjEnd':1})
        self.setInfo({'QrySeq':'X','SbjSeq':'X','AlnSeq':' '})
#########################################################################################################################
### END OF SECTION V: PWAln Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION VI: GENERAL MODULE METHODS                                                                                  #
#########################################################################################################################
def formatDB(fasfile,blastpath,protein=True,log=None):  ### Formats a blastDB
    '''
    Formats a blastDB.
    >> fasfile:str = file to format
    >> blastpath:str = path to BLAST programs
    >> protein:boolean = whether protein sequences
    >> log:Log Object
    '''
    try:### ~ [1] Setup Command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        path = rje.makePath(blastpath)
        #!# -a F option is to remove SeqPortNew errors with -o T (??!!)
        #command = blastpath + 'formatdb -i %s -o T -a F' % fasfile
        if os.sep == '\\': command = blastpath + 'formatdb -i "%s" -o T' % fasfile
        else: command = blastpath + 'formatdb -i %s -o T -a F' % fasfile
        if protein: command += ' -p T'
        else: command += ' -p F'
        ### ~ [2] Execute command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if log != None: log.printLog('#DB ',command,True)
        os.system(command)
    except:
        if log: log.errorLog('Major Problem during rje_blast.formatDB(%s).' % fasfile)
        else: print('Major Problem during rje_blast.formatDB(%s).' % fasfile)
        raise
#########################################################################################################################
def checkForDB(dbfile=None,checkage=True,log=None,protein=True):     ### Checks for BLASTDB files and returns True or False as appropriate
    '''
    Checks for BLASTDB files and returns True or False as appropriate.
    >> dbfile:str = sequence file forming basis of database
    >> checkage:Boolean = also check that blastdb files are newer than dbfile
    >> log:Log Object
    >> protein:boolean = whether database is protein
    '''
    try:### ~ [1] ~ Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(dbfile):
            if log: log.errorLog('%s missing' % dbfile,False,False)
            return False
        if protein: suffix = ['phr','pin','psd','psi','psq']
        else: suffix = ['nhr','nin','nsd','nsi','nsq']
        missing = []
        for suf in suffix:
            if os.path.exists('%s.%s' % (dbfile,suf)):
                if checkage and rje.isYounger('%s.%s' % (dbfile,suf),dbfile) == dbfile:
                    if log: log.errorLog('%s.%s too old' % (dbfile,suf),False,False)
                    return False
            elif os.path.exists('%s.00.%s' % (dbfile,suf)):
                if checkage and rje.isYounger('%s.00.%s' % (dbfile,suf),dbfile) == dbfile:
                    if log: log.errorLog('%s.%s too old' % (dbfile,suf),False,False)
                    return False
            else: missing.append(suf)
        if missing:
            if log and len(missing) < 5: log.errorLog('%s.%s missing' % (dbfile,string.join(missing,'/')),False,False)
            return False
        return True
    except:
        if log: log.errorLog('Major Problem during rje_blast.checkForDB().')
        else: print('Major Problem during rje_blast.checkForDB().')
        raise
#########################################################################################################################
def cleanupDB(callobj=None,dbfile=None,deletesource=False):     ### Deletes files created by formatdb
    '''
    Deletes files created by formatdb.
    >> callobj:Object = object calling method
    >> dbfile:str = sequence file forming basis of database
    >> deletesource:boolean = whether to delete dbfile as well
    '''
    try:
        for suf in ['phr','pin','psd','psi','psq','pal','nhr','nin','nsd','nsi','nsq','pni','pnd']:
            if os.path.exists('%s.%s' % (dbfile,suf)): os.unlink('%s.%s' % (dbfile,suf))
            for x in range(100):
                xstr = rje.preZero(x,99)
                if os.path.exists('%s.%s.%s' % (dbfile,xstr,suf)): os.unlink('%s.%s.%s' % (dbfile,xstr,suf))
        if deletesource and os.path.exists(dbfile): os.unlink(dbfile)
    except:
        if callobj: callobj.log.errorLog('Major Problem during rje_blast.cleanupDB().')
        else: print('Major Problem during rje_blast.cleanupDB().')
        raise
#########################################################################################################################
def expectString(_expect):  ### Returns formatted string for _expect value
    '''Returns formatted string for _expect value.'''
    try:
        if _expect > 10: return '%.1f' % _expect
        elif _expect > 0.1: return '%.2f' % _expect
        elif _expect > 0.001: return '%.3f' % _expect
        else: return '%.2e' % _expect
    except: print(expect); raise
#########################################################################################################################
### END OF SECTION VI                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION VII: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print('Unexpected error during program setup:', sys.exc_info()[0])
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        me = BLASTRun(mainlog,cmd_list)     # Will formatdb 
        me.blast(use_existing=not me.opt['Force'],log=True)           
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################