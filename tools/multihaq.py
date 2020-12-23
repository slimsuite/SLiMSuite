#!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      MultiHAQ
Description:  Multi-Query HAQESAC controller
Version:      1.4.3
Last Edit:    31/07/20
Citation:     Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: 20924652]
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a wrapper for multiple HAQESAC runs where different query proteins are to be BLASTed against the same
    search database(s) and run through HAQESAC with the same settings. The default expectation is that some queries will
    be returned by the HAQESAC runs of other queries and may therefore be skipped as a result, although this can be
    switched off using screenqry=F. For large runs, the first phase of MulitHAQ will take a long time to run. In these
    cases, it may be desirable to set the second, interactive, phase running before it has finished. This is achieved
    using the "chaser" option, which will set the second phase in motion, "chasing" the progress of the first. To avoid
    jumbled log output, this should be given a different log file using log=FILE.

    Note: that all options will be output into a haqesac.ini file in the haqdir path, for both HAQESAC runs within the
    framework of MultiHAQ itself and also for later runs using the batch file produced. Any generic HAQESAC options
    should therefore be placed into a multihaq.ini file, not a haqesac.ini file and multiple runs with different settings
    using the same haqdir should be avoided.

    Note: Because HAQESAC makes use of RJE_SEQ filtering options, they will NOT be applied to the MultiHAQ query input
    file prior to analysis. To filter this input, run it through rje_seq.py separately in advance of running multihaq.

Commandline:
    ### ~~~ INPUT OPTIONS ~~~ ###
    seqin=FILE      : Input query sequences [None]
    blast2fas=LIST  : List of databases to BLAST queries against prior to HAQESAC []
    addqueries=T/F  : Whether to add query database to blast2fas list [True]
    keepblast=T/F   : Whether to keep BLAST results files [True]
    blastcut=X      : Restrict HAQESAC and MultiHAQ BLAST searches to top X BLAST hits [0]
    multicut=X      : Restrict MultiHAQ BLASTs to the top X hits from each database (over-rides blastcut) [0]

    ### ~~~ MULTIHAQ OPTIONS ~~~ ###
    haqesac=T/F     : Run HAQESAC (True) or limit to batch file output (False) [True]
    multihaq=T/F    : Whether to run HAQESAC in two-phase multihaq mode [True]
    screenqry=T/F   : Whether to look for queries in previous runs and give option to skip [True]
    autoskip=T/F    : Whether to automatically skip queries found in previous runs [False]
    chaser=T/F      : Option for running second phase of multihaq as second run whilst first run in progress [False]
    force=T/F       : Whether to force re-running of stages (True) or pick-up pre-existing runs (False) [False]

    ### ~~~ OUTPUT OPTIONS ~~~~ ###
    haqdir=PATH     : Directory in which to output HAQESAC files and perform run [seqin_HAQESAC]
    haqblastdir=PATH: Directory in which MultiHAQ BLAST2FAS BLAST runs will be performed [./HAQBLAST/]

See also haqesac.py commands and rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
### User modules - remember to add *.__doc__ to cmdHelp() below ###
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
import haqesac, rje, rje_seq, rje_zen
import slimfarmer
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added chaser and autoskip options.
    # 1.0 - Fully working version. Fixed minor basefile bug. Added blastcut filter.
    # 1.1 - Improved pickup of aborted run.
    # 1.2 - Changed defaults to autoskip=F.
    # 1.2.1 - Updated documentation to include the HAQESAC reference.
    # 1.2.2 - Switched default to keepblast=T. Added forking blasta=X command to BLAST.
    # 1.3.0 - MultiCut : Restrict BLAST to the top X hits from each database [100]
    # 1.4.0 - Added SLiMFarmer batch forking if autoskip=F and i=-1.
    # 1.4.1 - Added haqblastdir=PATH: Directory in which MultiHAQ BLAST2FAS BLAST runs will be performed [./HAQBLAST/]
    # 1.4.2 - Fixed issue with SLiMFarmer for i<0 runs.
    # 1.4.3 - Updated warnings if BLAST2FAS files not found.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add forking for Phase I multihaq runs?
    # [ ] : Upgrade to new object type.
    # [ ] : Add HTML output. Currently coded up in ProtHunter.
    # [Y] : Separate the MultiHAQ and HAQESAC blastcut values.
    # [ ] : Name the _HAQESAC output directory using basefile, not seqin, by default.
    # [ ] : Add an error in the BLAST2Fas does not find any BLAST databases.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('MULTIHAQ', '1.4.2', 'July 2020', '2009')
    description = 'Multi-Query HAQESAC controller'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please cite: Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504.',
                'For HAQESAC results, please cite: Edwards et al. (2007), Nature Chem. Biol. 3(2):108-112.']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show HAQESAC commandline options?'): out.verbose(-1,4,text=haqesac.__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
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
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: MultiHAQ Class                                                                                          #
#########################################################################################################################
class MultiHAQ(rje.RJE_Object):     
    '''
    MultiHAQ Class. Author: Rich Edwards (2009).

    Info:str
    - HaqDir = Directory in which to output HAQESAC files and perform run [seqin_HAQESAC]
    - HAQBLASTDir=PATH: Directory in which MultiHAQ BLAST2FAS BLAST runs will be performed [./HAQBLAST/]
    
    Opt:boolean
    - AddQueries = Whether to add query database to blast2fas list [True]
    - AutoSkip = Whether to automatically skip queries found in previous runs [False]
    - Chaser = Option for running second phase of multihaq as second run whilst first run in progress [False]
    - Force = Whether to force re-running of stages (True) or pick-up pre-existing runs (False) [False]
    - HAQESAC = Run HAQESAC (True) or limit to batch file output (False) [True]
    - MultiHAQ = Whether to run HAQESAC in two-phase multihaq mode [True]
    - ScreenQry = Whether to look for queries in previous runs and give option to skip [True]

    Stat:numeric
    - BlastCut = Maximum number of sequences to have in dataset (BLAST query against NR dataset.)
    - MultiCut = Restrict BLAST to the top X hits from each database [0]

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - SeqList = Main query sequence input SeqList object
    - SLiMFarmer = SLiMFarmer batch forking object.
    '''
#########################################################################################################################
    def name(self): return self.obj['SeqList'].info['Name']
    def seqNum(self): return self.obj['SeqList'].seqNum()
    def seqs(self): return self.obj['SeqList'].seqs()
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Blast2Fas','HaqDir','HAQBLASTDir']
        self.optlist = ['AddQueries','AutoSkip','Chaser','HAQESAC','MultiHAQ','ScreenQry']
        self.statlist = ['BlastCut','MultiCut']
        self.listlist = []
        self.dictlist = []
        self.objlist = ['SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.basefile('MultiHAQ')
        self.setOpt({'Chaser':False,'AutoSkip':False})
        self.setStr({'HAQBLASTDir':rje.makePath('./HAQBLAST/')})
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
                ### Class Options ### 
                self._cmdReadList(cmd,'info',['Blast2Fas'])
                self._cmdReadList(cmd,'path',['HaqDir','HAQBLASTDir'])
                self._cmdReadList(cmd,'int',['BlastCut','MultiCut'])
                self._cmdReadList(cmd,'opt',['AddQueries','AutoSkip','Chaser','HAQESAC','MultiHAQ','ScreenQry'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return self.printLog('#ERR','Setup failed: quitting MultiHAQ')
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blast2fas = not self.opt['Chaser'] and self.blast2fas()
            if not self.opt['Chaser']: self.haqBatch(force=blast2fas)
            if self.opt['HAQESAC']:
                if self.i() < 0 and not self.getBool('AutoSkip') and self.seqNum() > 1 and self.getInt('Forks') > 1:
                    self.farmHAQ()
                else:
                    self.multiHAQ(secondrun=self.opt['Chaser'])
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = rje_seq.SeqList(self.log,['keepblast=T']+self.cmd_list+['autofilter=F','align=F','haqbat=None'])
            self.obj['SeqList']._checkForDup(True)
            if not self.seqNum(): self.errorLog('No sequences loaded!',printerror=False); return False
            if not self.obj['SeqList'].list['Blast2Fas']:
                if not self.getStrLC('Blast2Fas'): self.warnLog('No blast2fas=FILES paths given. Check for rogue spaces etc. unless only self-searching seqin=FILE',quitchoice=True)
            if self.opt['AddQueries'] and self.name() not in self.obj['SeqList'].list['Blast2Fas']: self.obj['SeqList'].list['Blast2Fas'].append(self.name())
            ### ~ [2] Setup Results Directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['HaqDir'].lower() in ['','none']: self.info['HaqDir'] = '%s_HAQESAC/' % rje.baseFile(self.name(), strip_path=True)
            rje.mkDir(self,self.info['HaqDir'])
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def blast2fas(self):    ### Executes BLAST2FAS and copies results files
        '''Executes BLAST2FAS and copies results files.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            need2blast = self.opt['Force']
            null_file = '%s.blast2fas_null.txt' % self.baseFile(); nx = 0; null_list = []
            if os.path.exists(null_file): null_list = string.split(open(null_file,'r').read(),'\n')
            self.debug(null_file)
            for seq in self.seqs():
                if seq.info['AccNum'] in null_list: nx += 1; continue
                hfile = rje.makePath('%s%s.fas' % (self.info['HaqDir'],seq.info['AccNum']),wholepath=True)
                for db in self.obj['SeqList'].list['Blast2Fas']:
                    self.debug(rje.isYounger(hfile,db))
                    self.debug(rje.isYounger(hfile,db) == hfile)
                    need2blast = need2blast or not rje.isYounger(hfile,db) == hfile
            if not need2blast:
                self.printLog('#BLAST','All HAQESAC input files found (%s w/o BLAST hits) - no BLAST2Fas (force=F)' % nx)
                return False
            ### ~ [2] Execute ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.backup(self,null_file); nx = 0
            if self.getInt('MultiCut'): self.obj['SeqList'].cmd_list += ['blastb=%d' % self.getInt('MultiCut'),'blastv=%d' % self.getInt('MultiCut')]
            elif self.getInt('BlastCut'): self.obj['SeqList'].cmd_list += ['blastb=%d' % self.getInt('BlastCut'),'blastv=%d' % self.getInt('BlastCut')]
            if self.getInt('Forks'): self.obj['SeqList'].cmd_list += ['blasta=%d' % self.getInt('Forks')]
            rje_seq.Blast2Fas(self.obj['SeqList'],self.getStr('HAQBLASTDir'))
            for seq in self.seqs():
                sbfile = '%s%s.blast.fas' % (self.getStr('HAQBLASTDir'),seq.info['AccNum'])
                if os.path.exists(sbfile):
                    hfile = rje.makePath('%s%s.fas' % (self.info['HaqDir'],seq.info['AccNum']),wholepath=True)
                    os.rename(sbfile,hfile)
                    if os.path.exists('%s.pickle' % rje.baseFile(hfile)): os.unlink('%s.pickle' % rje.baseFile(hfile))
                    if os.path.exists('%s.pickle.gz' % rje.baseFile(hfile)): os.unlink('%s.pickle.gz' % rje.baseFile(hfile))
                else: open(null_file,'a').write('%s\n' % seq.info['AccNum']); nx += 1
            if nx: self.printLog('#BLAST','%s Accession Numbers without BLAST2Fas hits output to %s' % (nx,null_file))
            self.printLog('#BLAST','%s HAQESAC input files made using BLAST2Fas' % (self.seqNum()-nx))
            return True
        except: self.errorLog('Major problem with MultiHAQ.blast2fas'); raise
#########################################################################################################################
    def haqBatch(self,force=False): ### Generates Batch and INI files for HAQESAC runs
        '''Generates Batch and INI files for HAQESAC runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batfile = rje.makePath('%shaqesac.bat' % self.info['HaqDir'],wholepath=True)
            inifile = rje.makePath('%shaqesac.ini' % self.info['HaqDir'],wholepath=True)
            if force or self.force() or not rje.exists(batfile) or not rje.exists(inifile): rje.backup(self,batfile); rje.backup(self,inifile)
            else: return self.printLog('#HAQBAT','HAQESAC Batch files found.')
            ### ~ [1] Make INI File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            haqcmd = []
            for cmd in self.cmd_list:
                if cmd[:4].lower() != 'ini=': haqcmd.append(cmd)
            if self.opt['MultiHAQ']: haqcmd += ['multihaq=T','force=F']
            open(inifile,'w').write(string.join(haqcmd,'\n'))
            ### ~ [2] Make Batch file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                acc = seq.info['AccNum']
                haqcmd = ['seqin=%s.fas' % acc, 'query=%s' % acc, 'basefile=%s' % acc]
                open(batfile,'a').write('python %shaqesac.py %s\n' % (self.info['Path'],string.join(haqcmd)))
            self.printLog('#HAQBAT','HAQESAC Batch file output to %s' % batfile)
        except: self.errorLog('Major problem with MultiHAQ.haqBatch',quitchoice=True)
#########################################################################################################################
    def multiHAQ(self,secondrun=False):     ### Executes main HAQESAC runs
        '''Executes main HAQESAC runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            finalrun = secondrun == self.opt['MultiHAQ']    # Whether this is the manual HAQESAC phase
            qryacc = self.obj['SeqList'].accList()          # Full list of Query accession numbers
            processed = []                                  # List of processed sequence accession numbers
            ### ~ [1] Peform HAQESAC runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                ## ~ [1a] Check AutoSkip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                acc = seq.info['AccNum']
                if finalrun and acc in processed and (self.opt['AutoSkip'] or (self.i() >=0 and rje.yesNo('%s already covered by previous HAQESAC. Skip?' % seq.shortName()))):
                    self.printLog('#SKIP','%s already covered by previous HAQESAC: Skipped' % seq.shortName()); continue
                ## ~ [1b] Check Whether to run (re-runs and low sequence number) ~~~~~~~~~~~~~~~~~~ ##
                logfile = rje.makePath('%s%s.log' % (self.info['HaqDir'],acc),wholepath=True)
                infile = rje.makePath('%s%s.fas' % (self.info['HaqDir'],acc),wholepath=True)
                pkfile = rje.makePath('%s%s.pickle' % (self.info['HaqDir'],acc),wholepath=True)
                pkzfile = rje.makePath('%s%s.pickle.gz' % (self.info['HaqDir'],acc),wholepath=True)
                if not os.path.exists(infile): self.printLog('#SKIP','%s input file %s not found: Skipped' % (seq.shortName(),infile)); continue
                if not finalrun and not self.opt['Force'] and rje.isYounger(pkzfile,infile) == pkzfile:
                    self.printLog('#SKIP','%s run detected: Skipped' % seq.shortName()); continue
                if not finalrun and not self.opt['Force'] and rje.isYounger(pkfile,infile) == pkfile:
                    self.printLog('#SKIP','%s run detected: Skipped' % seq.shortName()); continue
                inseqx = rje_seq.SeqCount(self,infile)
                if inseqx < 2: self.printLog('#SKIP','Only one sequence found in %s: Skipped' % (infile)); continue
                ## ~ [1c] Pause if running in Chaser Mode and no Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pickled = os.path.exists(pkfile) or os.path.exists('%s.gz' % pkfile); tm = 0
                while secondrun and self.opt['Chaser'] and not pickled:
                    self.progLog('#WAIT','No %s pickle. Sleeping for %d min.' % (acc,tm))
                    time.sleep(60*tm); tm += 1
                    pickled = os.path.exists(pkfile) or os.path.exists('%s.gz' % pkfile)
                    if not pickled:
                        try: rje.choice('Press <ENTER> to try again, or <CTRL+C> to Quit')
                        except:
                            self.printLog('#PICKLE','No %s pickle.' % (acc,tm))
                            self.printLog('\r#MULTI','Exiting multiHAQ "Chaser" run.'); return
                ## ~ [1d] Run HAQESAC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                runhaqesac = True
                pngfile = rje.makePath('%s%s.png' % (self.info['HaqDir'],acc),wholepath=True)
                if not self.force() and rje.exists(pngfile):
                    self.printLog('#SKIP','Found evidence of completed run: %s (force=F). Skipping.' % pngfile)
                    runhaqesac = False
                ancfile = rje.makePath('%s%s.anc.fas' % (self.info['HaqDir'],acc),wholepath=True)
                if not self.force() and rje.exists(ancfile):
                    self.printLog('#SKIP','Found evidence of completed run: %s (force=F). Skipping.' % ancfile)
                    runhaqesac = False
                #if not finalrun or self.opt['Force'] or rje.isYounger(logfile,nsfile) != logfile:
                if runhaqesac:
                    haqcmd = ['ini=haqesac.ini','seqin=%s.fas' % acc, 'query=%s' % acc, 'basefile=%s' % acc, 'newlog=F']
                    self.printLog('#HAQ','Running HAQESAC for %s - will have own log etc.' % seq.shortName(),log=False)
                    os.chdir(self.info['HaqDir'])
                    info = haqesac.makeInfo()
                    haqcmd = rje.getCmdList(haqcmd,info=info)
                    out = rje.Out(cmd_list=haqcmd)    # Sets up Out object for controlling output to screen
                    out.printIntro(info)                                # Prints intro text using details from Info object
                    haqlog = rje.setLog(info,out,haqcmd)                 # Sets up Log object for controlling log file output
                    try: haqesac.HAQESAC(log=haqlog, cmd_list=haqcmd).run(setobjects=True)
                    except:
                        os.chdir(self.info['RunPath'])
                        if self.i() >= 0 and rje.yesNo('Problem with %s HAQESAC run. Abort?' % seq.shortName()): raise KeyboardInterrupt
                    os.chdir(self.info['RunPath'])
                    if finalrun: self.printLog('#HAQ','HAQESAC final round run for %s' % seq.shortName())
                    else: self.printLog('#HAQ','HAQESAC first round run for %s' % seq.shortName())
                ## ~ [1e] Update ScreenQry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.opt['ScreenQry'] or not finalrun: continue
                qacclist = []
                for qacc in rje_seq.SeqList(self.log,['seqin=%s' % infile,'autoload=T','autofilter=F']).accList():
                    if qacc in qryacc and qacc != acc: qacclist.append(qacc)
                    if qacc in qryacc and qacc not in processed: processed.append(qacc)
                self.printLog('#QRY','%d other queries found in %s: [%s]' % (len(qacclist),infile,string.join(qacclist,'; ')))
                self.printLog('#QRY','%d of %d queries processed' % (len(processed),self.seqNum()))
            ### ~ [2] MultiHAQ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not finalrun: self.printLog('#MULTI','Executing second round of multiHAQ'); self.multiHAQ(True)
        except: self.errorLog('Major problem with MultiHAQ.multiHAQ',quitchoice=True)
#########################################################################################################################
    def farmHAQ(self):  ### Uses SLiMFarmer to farm out the HAQESAC runs
        '''Uses SLiMFarmer to farm out the HAQESAC runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batfile = os.path.abspath(rje.makePath('%shaqesac.bat' % self.info['HaqDir'],wholepath=True))
            self.printLog('#FARM',batfile)
            if not rje.exists(batfile): raise IOError('Cannot find %s' % batfile)
            farmcmd = ['subjobs=%s' % batfile,'farm=batch','qsub=F','i=-1','runpath=%s' % os.path.abspath(self.info['HaqDir'])]
            if self.opt['MultiHAQ']:
                haqfarm = ['First round','Second round']
            else: haqfarm = ['Complete run']

            ### ~ [1] Peform HAQESAC runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for farmrun in haqfarm:
                self.printLog('#CHDIR','Changing directory for %s farming: %s' % (farmrun,self.info['HaqDir']))
                os.chdir(self.info['HaqDir'])
                farmer = slimfarmer.SLiMFarmer(self.log,self.cmd_list+farmcmd)
                farmer.slimFarm()
                os.chdir(self.info['RunPath'])
                self.printLog('#CHDIR','Changed directory post-farming: %s' % self.info['RunPath'])
                self.printLog('#FARM','HAQESAC %s farming complete.' % farmrun)
            return True

            #!# Add identifying and skipping of partial runs.

            for seq in self.seqs():
                ## ~ [1a] Check AutoSkip ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                acc = seq.info['AccNum']
                if finalrun and acc in processed and (self.opt['AutoSkip'] or (self.i() >=0 and rje.yesNo('%s already covered by previous HAQESAC. Skip?' % seq.shortName()))):
                    self.printLog('#SKIP','%s already covered by previous HAQESAC: Skipped' % seq.shortName()); continue
                ## ~ [1b] Check Whether to run (re-runs and low sequence number) ~~~~~~~~~~~~~~~~~~ ##
                logfile = rje.makePath('%s%s.log' % (self.info['HaqDir'],acc),wholepath=True)
                infile = rje.makePath('%s%s.fas' % (self.info['HaqDir'],acc),wholepath=True)
                pkfile = rje.makePath('%s%s.pickle' % (self.info['HaqDir'],acc),wholepath=True)
                pkzfile = rje.makePath('%s%s.pickle.gz' % (self.info['HaqDir'],acc),wholepath=True)
                if not os.path.exists(infile): self.printLog('#SKIP','%s input file %s not found: Skipped' % (seq.shortName(),infile)); continue
                if not finalrun and not self.opt['Force'] and rje.isYounger(pkzfile,infile) == pkzfile:
                    self.printLog('#SKIP','%s run detected: Skipped' % seq.shortName()); continue
                if not finalrun and not self.opt['Force'] and rje.isYounger(pkfile,infile) == pkfile:
                    self.printLog('#SKIP','%s run detected: Skipped' % seq.shortName()); continue
                inseqx = rje_seq.SeqCount(self,infile)
                if inseqx < 2: self.printLog('#SKIP','Only one sequence found in %s: Skipped' % (infile)); continue
                ## ~ [1c] Pause if running in Chaser Mode and no Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pickled = os.path.exists(pkfile) or os.path.exists('%s.gz' % pkfile); tm = 0
                while secondrun and self.opt['Chaser'] and not pickled:
                    self.progLog('#WAIT','No %s pickle. Sleeping for %d min.' % (acc,tm))
                    time.sleep(60*tm); tm += 1
                    pickled = os.path.exists(pkfile) or os.path.exists('%s.gz' % pkfile)
                    if not pickled:
                        try: rje.choice('Press <ENTER> to try again, or <CTRL+C> to Quit')
                        except:
                            self.printLog('#PICKLE','No %s pickle.' % (acc,tm))
                            self.printLog('\r#MULTI','Exiting multiHAQ "Chaser" run.'); return
                ## ~ [1d] Run HAQESAC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                runhaqesac = True
                pngfile = rje.makePath('%s%s.png' % (self.info['HaqDir'],acc),wholepath=True)
                if not self.force() and rje.exists(pngfile):
                    self.printLog('#SKIP','Found evidence of completed run: %s (force=F). Skipping.' % pngfile)
                    runhaqesac = False
                ancfile = rje.makePath('%s%s.anc.fas' % (self.info['HaqDir'],acc),wholepath=True)
                if not self.force() and rje.exists(ancfile):
                    self.printLog('#SKIP','Found evidence of completed run: %s (force=F). Skipping.' % ancfile)
                    runhaqesac = False

        except:
            os.chdir(self.info['RunPath'])
            self.errorLog('Major problem with MultiHAQ.farmHAQ',quitchoice=True)
#########################################################################################################################
### End of SECTION II: MultiHAQ Class                                                                                   #
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
    try: MultiHAQ(mainlog,cmd_list).run()

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
