#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_qsub
Description:  QSub Generating module
Version:      1.12.0
Last Edit:    13/02/23
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to make a job file and call it with qsub. Walltime and depend commands are included. 

Commandline:
    - program=X     : Program call for Qsub (and options) [None]
    - job=X         : Name of job file (.job added) [qsub]
    - qpath=PATH    : Path to change directory too [current path]
    - rjepy=T/F     : Whether program is an RJE *.py script (adds python PyPath/) [True]
    - pypath=PATH   : Path for RJE Python scripts [/rhome/re1u06/Serpentry/]
    - nodes=X       : Number of nodes to run on [4]
    - ppn=X         : Processors per node [12]
    - walltime=X    : Walltime for qsub job (hours) [60]
    - depend=LIST   : List of job ids to wait for before starting job (dependhpc=X added) []
    - pause=X       : Wait X seconds before attempting showstart [5]
    - report=T/F    : Pull out running job IDs and run showstart [False]
    - email=X       : Email address to email job stats to at end ['']
    - mailstart=T/F : Whether to email user at start of run [False]
    - hpc=X         : Name of HPC system for depend ['IRIDIS4']
    - dependhpc=X   : Name of HPC system for depend (will parse from hostname if blank) ['kman.restech.unsw.edu.au']
    - vmem=X        : Virtual Memory limit for run (GB) [48]
    - modules=LIST  : List of modules to add in job file []
    - modpurge=T/F  : Whether to purge loaded modules in qsub job file prior to loading [True]
    - precall=LIST  : List of additional commands to run between module loading and program call []
    - jobwait=T/F   : Whether to wait for the job to finish before exiting [False]
    - monitor=T/F   : Whether to add "-k oed" to enable instant job monitoring. [False]
    - startbash=T/F : Whether to add "-S /bin/bash" to qsub command [False]

Uses general modules: copy, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Initial working script.
    # 1.1 - Reworked for running on Iridis.
    # 1.2 - Added showstart.
    # 1.3 - Converted for IRIDIS3.
    # 1.4 - Added showstart pause.
    # 1.5 - Added emailing of job stats after run. Added vmem limit.
    # 1.6 - Added modules=LIST : List of modules to add in job file [clustalo,mafft]
    # 1.6.1 - Added R/3.1.1 to modules.
    # 1.6.2 - Updated module list: blast+/2.2.30,clustalw,clustalo,fsa,mafft,muscle,pagan,R/3.1.1
    # 1.6.3 - Tweaked the showstart command for katana.
    # 1.7.0 - Added option for email when job started
    # 1.8.0 - Added modpurge=T/F : Whether to purge loaded modules in qsub job file prior to loading [True]
    # 1.9.0 - Added precall=LIST  : List of additional commands to run between module loading and program call []
    # 1.9.1 - Removed default module list: causing conflicts. Better to have in INI file.
    # 1.9.2 - Modified qsub() to return job ID.
    # 1.9.3 - Updates the order of the qsub -S /bin/bash flag.
    # 1.9.4 - Updated qsub to use mem flag rather than vmem. Replaced showstart with qsub -T.
    # 1.9.5 - Added "-k oed" to enable instant job monitoring and option to drop -S /bin/bash/. Tidied walltime.
    # 1.9.6 - Updated qsub to use both mem and vmem flags.
    # 1.10.0 - Added jobwait to wait until the job has finished.
    # 1.10.1 - New default dependhpc=X  ['kman.restech.unsw.edu.au']
    # 1.10.2 - Added notification that job is running when jobwait=T.
    # 1.11.0 - Added output of date and time at end of job script too. (Gives a record of total time running.)
    # 1.11.1 - Added job run summary output to end of stdout.
    # 1.12.0 - Updated for new katana OS and qsub. Dropped vmem again. Switched nodes/ppn for select/ncpus.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Check out PBS options within job file
    # [Y] : Modify to work on HPC other than IRIDIS.
    # [ ] : Add i=-1 to call if rjepy=T.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copy_right) = ('RJE_QSUB', '1.12.0', 'February 2023', '2006')
    description = 'QSub Generating module'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
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
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #                                                                                                          #
#########################################################################################################################
class QSub(rje.RJE_Object):     
    '''
    QSub Generating Class. Author: Rich Edwards (2005).

    Info:str
    - Email = Email address to email job stats to at end ['']
    - HPC = Name of HPC system for depend ['IRIDIS4']
    - DependHPC = Name of HPC system for depend ['kman.restech.unsw.edu.au']
    - Program = Program call for Qsub (with options) [None]
    - PyPath = Path for Python scripts [python /rhome/re1u06/Serpentry/]
    - QPath = Path to change directory too [current path]
    - Job = Name of job file (.job added) [qsub]
    
    Opt:boolean
    - JobWait : Whether to wait for the job to finish before exiting [False]
    - MailStart = Whether to email user at start of run [False]
    - ModPurge=T/F : Whether to purge loaded modules in qsub job file prior to loading [True]
    - Monitor=T/F   : Whether to add "-k oed" to enable instant job monitoring. [False]
    - RjePy = Whether program is an RJE *.py script (adds python /home/richard/Python_Modules) [True]
    - Report = Pull out running job IDs and run showstart [False]
    - StartBash=T/F : Whether to add "-S /bin/bash" to qsub command [False]

    Stat:numeric
    - Nodes = Number of nodes to run on
    - Pause = No. seconds to pause prior to showstart [5]
    - PPN = Processors per node [12]
    - VMem = Virtual Memory limit for run (GB) [48]
    - Walltime = Walltime for qsub job (hours) [60]

    List:list
    - Depend = List of job ids to wait for before starting job []
    - Modules = List of modules to add in job file [blast+/2.2.30,clustalw,clustalo,fsa,mafft,muscle,pagan,R/3.1.1]
    - PreCall=LIST  : List of additional commands to run between module loading and program call []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Program','QPath','Job']
        - Opt:boolean ['ModPurge','RjePy']
        - Stats:float ['Nodes','Walltime','VMem']
        - List:list ['Depend','Modules']
        - Dict:dictionary []
        - Obj:RJE_Object []
        '''
        ### Basics ###
        self.infolist = ['Program','QPath','Job','PyPath','Email','HPC','DependHPC']
        self.optlist = ['JobWait','ModPurge','Monitor','RjePy','Report','MailStart','StartBash']
        self.statlist = ['Nodes','Walltime','PPN','Pause']
        self.listlist = ['Depend','PreCall','Modules']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
        self.setInfo({'QPath':os.path.abspath(os.curdir),'Job':'rje_%s' % rje.randomString(4),
                      'PyPath':'/home/re1u06/Serpentry/','Email':'',
                      'HPC':'IRIDIS4','DependHPC':'kman.restech.unsw.edu.au'})
        self.setStat({'Walltime':60,'Nodes':1,'PPN':12,'Pause':5,'VMem':48})
        self.setOpt({'JobWait':False,'Report':False,'MailStart':False,'ModPurge':True,'Monitor':False,'StartBash':False})
        #self.list['Modules'] = rje.split('blast+/2.2.30,clustalw,clustalo,fsa,mafft,muscle,pagan,R/3.1.1,fasttree,phylip',',')
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
                self._cmdReadList(cmd,'info',['Program','Job','Email','HPC','DependHPC'])
                self._cmdReadList(cmd,'path',['QPath','PyPath'])
                self._cmdReadList(cmd,'int',['Nodes','PPN','Pause','VMem'])
                self._cmdReadList(cmd,'stat',['Walltime'])
                self._cmdReadList(cmd,'opt',['JobWait','RjePy','Report','MailStart','ModPurge','Monitor','StartBash'])
                self._cmdRead(cmd,type='opt',att='JobWait',arg='wait')
                self._cmdReadList(cmd,'list',['Depend','Modules','PreCall'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStr('Email').lower() in ['none','']: self.info['Email'] = ''
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def qsub(self):      ### Creates job and calls with qsub
        '''Creates job and calls with qsub. Returns qsub job ID or 0 if jobwait=True and job completed.'''
        try:### Basics ###
            hr = int(self.stat['Walltime'])
            min = int((0.5+(self.stat['Walltime'] - hr)*60.0))
            if self.opt['Report']: return self.report()
            jobstr = rje.replace('%s.job' % self.info['Job'],'.job','')
            jlist = ['#!/bin/bash',
                     '#PBS -N %s' % jobstr,  #,'#PBS -q batch',
                     '#PBS -l select=%d:ncpus=%d' % (self.stat['Nodes'],self.stat['PPN']),
                     '#PBS -l walltime=%d:%s:00' % (hr,rje.preZero(min,60)),
                     '#PBS -l mem=%dgb' % self.getInt('VMem'),
                     '']     #10
            #if not os.popen('hostname').read().startswith('katana.science.unsw.edu.au'):
            #    jlist[-2] = '#PBS -l mem=%dgb' % self.getInt('VMem')
            if self.getBool('Monitor'):
                if self.getBool('JobWait'):
                    self.warnLog('Cannot run with wait=T and monitor=T: switched monitor=F')
                    self.setBool({'Monitor':False})
                else:
                    jlist += ['#PBS -k oed']
            if self.getStr('Email'):
                jlist += ['#PBS -M %s' % self.getStr('Email'),'#PBS -m ae']
                if self.getBool('MailStart'): jlist[-1] = '#PBS -m bae'
            jlist += ['### Define number of processors','NPROCS=`wc -l < $PBS_NODEFILE`',
                      'echo Running on host `hostname`','echo Time is `date`','echo Directory is `pwd`', #2
                      'echo This jobs runs on the following processors:','echo `cat $PBS_NODEFILE`','',                #5
                      'echo This job has allocated $NPROCS cpus','']
            self.printLog('#PPN','%d Node(s) requested: %d PPN.' % (self.getInt('Nodes'),self.getInt('PPN')))
            self.printLog('#MEM','%s GB Mem requested.' % (self.getStat('VMem')))
            if self.getBool('ModPurge'):
                jlist.append('module purge')
                self.printLog('#MOD','Modules purged (modpurge=T)')
            for mod in self.list['Modules']:
                if mod.lower() not in ['','none']: jlist.append('module add %s' % mod)
            if self.list['Modules']: self.printLog('#MOD','Modules added: %s' % rje.join(self.list['Modules'],'; '))
            for pcall in self.list['PreCall']:
                self.printLog('#PCALL',pcall)
                jlist.append(pcall)
            #x#jlist = ['#!/bin/sh']   # New Iridis shell script method!
            ### Directory & Program ###
            jlist.append('cd %s' % self.info['QPath'])
            pcall = self.info['Program']
            if self.opt['RjePy']: pcall = 'python ' + self.info['PyPath'] + pcall
            jlist.append(pcall)
            ### Completion message
            jlist += ['','echo ---','qstat -f $PBS_JOBID','echo ---']
            jlist += ['','echo','echo Time is `date`','echo Job complete']
            ### Output and call ###
            job = '{0}.job'.format(jobstr) #rje.replace('%s.job' % self.info['Job'],'.job.job','.job')
            open(job,'w').write(rje.join(jlist,'\n'))
            self.printLog('#DIR',self.info['QPath'])
            self.printLog('#RUN',pcall)
            #qsub = 'qsub %s -S /bin/sh -l walltime=%d:%d:00,nodes=%d:ppn=2' % (job,hr,min,self.stat['Nodes'])
            qsub = 'qsub'
            if self.getBool('StartBash'): qsub += ' -S /bin/bash'
            if self.list['Depend']:
                qsub += ' -W depend=afterany'
                #for id in self.list['Depend']: qsub += ':%s.bio-server' % id
                myhost = self.getStr('DependHPC')
                if not self.getStrLC('DependHPC'):
                    myhost = rje.split(os.popen('hostname').read())[0]
                for id in self.list['Depend']: qsub += ':%s.%s' % (id,myhost)
            qsub += ' %s' % (job)
            self.printLog('#JOB',qsub)
            if self.test():
                self.printLog('#TEST','Test mode: will not place job in queue.')
                self.verbose(0,1,rje.join(['>>>>>']+jlist+['<<<<<',''],'\n'))
                return False
            qrun = os.popen(qsub).read()
            self.printLog('#QSUB',qrun)
            qid = rje.split(qrun,'.')[0]
            showstart = 'qstat -T'
            if os.popen('hostname').read().startswith('katana.science.unsw.edu.au'):
                showstart = 'showstart'
            self.printLog('#SHOW','Attempt %s %s in %s sec' % (showstart,qrun,self.stat['Pause']),log=False)
            time.sleep(self.stat['Pause'])
            for qline in os.popen('%s %s' % (showstart,qrun)):   #qid):
                if rje.chomp(qline): self.printLog('#INFO', qline, timeout=False)

            ### Wait for job to be completed
            if self.getBool('JobWait'):
                if self.getBool('Monitor'): raise ValueError('Cannot run with wait=T and monitor=T')
                self.printLog('#WAIT','Waiting for job {0} to finish'.format(qid))
                ofile = '{0}.o{1}'.format(rje.replace('%s.job' % self.info['Job'],'.job',''),qid)
                running = False
                while not rje.exists(ofile):
                    qstat = rje.atoi( os.popen("qstat | grep '^{0}' -c".format(qid)).read().split()[0] )
                    if not qstat:
                        self.printLog('#QSTAT','Job {0} disappeared from qstat'.format(qid))
                        break
                    elif not running:
                        try:
                            qstat = rje.split( os.popen("qstat | grep '^{0}'".format(qid)).read().split()[4] )
                            if qstat == 'R':
                                running = True
                                self.printLog('#QSTAT','Job {0} running...'.format(qid))
                        except: pass
                    time.sleep( max(1,self.getInt('Pause')) )
                owait = 300
                while owait and not rje.exists(ofile):
                    owait -= 1
                    time.sleep(1)
                if rje.exists(ofile):
                    if 'Job complete' in os.popen('tail -n 1 {0}'.format(ofile)).read():
                        self.printLog('#DONE','{0} job ({1}) complete.'.format(jobstr, qid))
                        return 0
                    else:
                        self.printLog('#FAIL','{0} job ({1}) failed to finish.'.format(jobstr, qid))
                        return qid
                else:
                    self.printLog('#FAIL','{0} job ({1}) failed to generate {2}.'.format(jobstr, qid, ofile))

            return qid
        except: self.errorLog('Error in qsub()'); return False
#########################################################################################################################
    def report(self):   ### Run qstat to get job list then showstart on each job
        '''Run qstat to get job list then showstart on each job .'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qidlist = []
            qidjob = {}
            ### ~ [2] ~ Read in List of IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qline in os.popen('qstat'):
                try:
                    (qid,job) = rje.matchExp('^(\d+)\.\S+\s+(\S+)',qline)
                    qidlist.append(qid)
                    qidjob[qid] = job
                except: continue
            ### ~ [3] ~ Report ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#QSTAT','%d jobs in queue.' % len(qidlist))
            for qid in qidlist:
                self.printLog('#JOB', '%s = %s' % (qid,qidjob[qid]), timeout=False)
                for qline in os.popen('showstart %s' % qid):
                    if rje.chomp(qline): self.printLog('#INFO', qline, timeout=False)
            self.printLog('#ZEN',rje_zen.Zen().wisdom())
        except: self.errorLog('QSub.report problem')            
#########################################################################################################################
### End of SECTION II: QSub Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: QSub(mainlog,cmd_list).qsub()
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
