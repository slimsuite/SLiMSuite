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
Module:       Diploidocus
Description:  In silico diploid data generator.
Version:      0.1.0
Last Edit:    13/10/17
Copyright (C) 2017  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
    parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
    parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
    identifier table.)

    A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
    unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
    selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
    two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
    This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
    no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
    subreads will be added from the other parent if they reduce the difference in cumulative output for each parent, or
    until `lenfilter=X` is reached.

    Final output will be a `*.LXXXRQXX.fasta` file in which each parent has a similar total sequence content and for
    which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
    assemblies, where one parent has higher quality data than the other.

    NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
    higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
    minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
    relaxed. Similarly, only sequences above `lenfilter=X` in length will be output. These are the figures given in the
    `LXXXRQXX` part of the output file, e.g. defaults of RQ>=0.84 and Len>=500 generates `*.L500RQ84.fas`.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []
    basefile=FILE   : Root of output file names [diploidocus]
    genomesize=INT  : Haploid genome size (bp) [13.1e6]
    ### ~ Filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rqfilter=X      : Minimum RQ for output subreads [0.84]
    lenfilter=X     : Min read length for filtered subreads [500]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    See also SMRTSCAPE `summarise=T` options if `*.unique.tdt`/`*.smrt.tdt` have not been pre-generated with SMRTSCAPE.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist
import smrtscape
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Fixed bugs with parent basefile, genome size default and Sequel data parsing.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Diploidocus', '0.1.0', 'October 2017', '2017')
    description = 'In silico diploid data generator.'
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
### SECTION II: Diploidocus Class                                                                                       #
#########################################################################################################################
class Diploidocus(rje_obj.RJE_Object):
    '''
    Diploidocus Class. Author: Rich Edwards (2015).

    Str:str
    - Parent1=FOFN    : File of file names for subreads fasta files on Parent 1. []
    - Parent2=FOFN    : File of file names for subreads fasta files on Parent 2. []

    Bool:boolean

    Int:integer
    - GenomeSize=INT  : Haploid genome size (bp) [13.1e6]
    - LenFilter=X     : Min read length for filtered subreads [500]

    Num:float
    - RQFilter=X      : Minimum RQ for output subreads [0.84]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Parent1','Parent2']
        self.boollist = []
        self.intlist = ['GenomeSize','LenFilter']
        self.numlist = ['RQFilter']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({'LenFilter':500,'GenomeSize':13.1e6})
        self.setNum({'RQFilter':0.84})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'file',['Parent1','Parent2'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['GenomeSize','LenFilter'])   # Integers
                self._cmdReadList(cmd,'float',['RQFilter']) # Floats
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
            return self.inSilicoHybrid()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''): self.baseFile('diploidocus')
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Main Run Method                                                                                         #
#########################################################################################################################
    def inSilicoHybrid(self):  ### Filter and combine subreads from parent and output to fasta file.
        '''
        Filter and combine subreads from parent and output to fasta file.

        This module generates balanced "in silico diploid" PacBio subread data from two sequenced haploid parents. Each
        parent must first be run through SMRTSCAPE to generate subread summary data. (This will be performed if missing. Each
        parent needs a `*.fofn` file of subread file names, `*.unique.tdt` unique subreads table and `*.smrt.tdt` SMRT cell
        identifier table.)

        A new set of subreads is then generated from the combined set of parent subreads. This is done by first ranking the
        unique subreads from each parent by length. First, the longest subread from each parent are compared and the shortest
        selected to be the first subread of the diploid. (The shortest is taken to minimise length differences between the
        two parents.) Next, the longest subread from the next parent that is no longer than the previous subread is added.
        This cycles, picking a read from the the parent with fewest cumulative bases each cycle. The longest subread that is
        no longer than the previous subread is selected. This continues until one parent runs out of subreads. Additional
        subreads will be added from the other parent if they reduce the difference in cumulative output for each parent.

        Final output will be a `*.subreads.fasta` file in which each parent has a similar total sequence content and for
        which the subread length distributions should also be similar. This is to overcome biases in resulting diploid
        assemblies, where one parent has higher quality data than the other.

        NOTE: If performing downstream filtering by Read Quality (RQ), this might reintroduce a bias if one parent has much
        higher RQ values than the other. The `rqfilter=X` setting can therefore be used to restrict output to  reads with a
        minimum RQ value. By default this is 0.84. If you do not get enough sequence output, this setting may need to be
        relaxed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Parent 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 1 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent1: %s' % self.getStr('Parent1'))
            base1 = rje.baseFile(self.getStr('Parent1'))
            parent1 = smrtscape.SMRTSCAPE(self.log,['genomesize=13.1e6']+self.cmd_list+['batch=%s' % self.getStr('Parent1'),'basefile=%s' % base1])
            parent1.setup()
            udb1 = parent1.udb()
            cdb = parent1.db('smrt',add=True,mainkeys=['Name'])
            cdb.dataFormat({'SMRT':'int'})
            cx = cdb.entryNum()
            ## ~ [0a] Parent 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ SETUP PARENT 2 ~~~~~~~~~~~~~~~~~~~~ #')
            self.printLog('#FOFN','Parent2: %s' % self.getStr('Parent2'))
            base2 = rje.baseFile(self.getStr('Parent2'))
            parent2 = smrtscape.SMRTSCAPE(self.log,['genomesize=13.1e6']+self.cmd_list+['batch=%s' % self.getStr('Parent2'),'basefile=%s' % base2])
            parent2.setup()
            udb2 = parent2.udb()
            cdb2 = parent2.db('smrt',add=True,mainkeys=['Name'])
            cdb2.dataFormat({'SMRT':'int'})
            # Shift all of the Parent2 SMRT IDs to avoid conflict with Parent1
            for entry in cdb2.entries() + udb2.entries(): entry['SMRT'] = entry['SMRT'] + cx
            cdb = parent1.db().mergeTables(cdb,cdb2)
            ## ~ [0c] Output Sequence File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~ DIPLOIDOCUS SUBREADS ~~~~~~~~~~~~~~~~~~~~ #')
            minlen = self.getInt('LenFilter')
            minrq = self.getNum('RQFilter')
            rqstr = '%s' % minrq
            filtfile = '%s.L%sRQ%s.fasta' % (self.baseFile(),minlen,rqstr[2:])
            ## ~ [0d] Input Sequence Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqbatch = []   # List of SeqList objects
            self.printLog('#BATCH','%s sequence files to process.' % rje.iLen(parent1.list['Batch']+parent2.list['Batch']))
            for seqfile in parent1.list['Batch']+parent2.list['Batch']:
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=F','seqin=%s' % seqfile,'autofilter=F']
                seqbatch.append(rje_seqlist.SeqList(self.log,seqcmd))
            self.printLog('#BATCH','%s sequence files to summarise.' % rje.iLen(seqbatch))
            if not seqbatch: raise IOError('No batch input fasta files found! Make sure parentN=FILE settings given *.fofn.')
            ## ~ [0e] Setup subread lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elists = [udb1.sortedEntries('Len',reverse=True),udb2.sortedEntries('Len',reverse=True)]
            plen = [0,0]    # Summed lengths for each parent
            pseq = [0,0]    # Total sequence number for each parent
            prq = [0,0]     # Total sequence RQ for each parent (convert to mean)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            lastlen = max(elists[0][0]['Len'],elists[1][0]['Len'])    # Length of last selected read
            for elist in elists:
                while elist and elist[0]['RQ'] < minrq: elist.pop(0)
            if not elists[0] or not elists[1]: raise ValueError('No Unique ZMW subreads for one or both parents!')
            nextp = 0       # Index of next parent to use
            if elists[0][0]['Len'] < elists[1][0]['Len']: nextp = 1

            ### ~ [1] Filter and Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Filter Unique Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            zmwlist = []    # List of (smrt,zmw) meeting filtering criteria
            ux = 0.0; utot = len(elists[0])+len(elists[1])
            while lastlen:
                self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq: elist.pop(0); ux += 100.0
                if elist and elist[0]['Len'] < minlen: ux += 100.0 * len(elist); elist = []
                if not elist: nextp = 1 - nextp; break  # Finish
                entry = elist.pop(0); ux += 100.0
                zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                plen[nextp] += entry['Len']
                prq[nextp] += entry['RQ']
                pseq[nextp] += 1
                if plen[1-nextp] <= plen[nextp]: nextp = 1 - nextp
                lastlen = entry['Len']
            ## ~ [1b] Final processing of last reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while elists[nextp]:
                elist = elists[nextp]
                while elist and elist[0]['RQ'] < minrq:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    elist.pop(0); ux += 100.0
                while elist and elist[0]['Len'] >= minlen:
                    self.progLog('\r#DIP','Diploidising subreads: %.2f%%' % (ux/utot))
                    entry = elist.pop(0); ux += 100.0
                    pdiff = rje.modulus(plen[0]-plen[1])
                    ediff = rje.modulus(plen[nextp]+entry['Len']-plen[1-nextp])
                    if ediff >= pdiff: elists[nextp] = []; break    #Finish!
                    zmwlist.append((entry['SMRT'],entry['ZMW'],entry['Pos']))
                    plen[nextp] += entry['Len']
                    prq[nextp] += entry['RQ']
                    pseq[nextp] += 1
            self.printLog('\r#DIP','Diploidising subreads complete: %s subreads to output.' % rje.iLen(zmwlist))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent1'),rje.iStr(pseq[0]),rje.iStr(plen[0]),1.0*plen[0]/self.getInt('GenomeSize'),prq[0]/pseq[0]))
            self.printLog('\r#DIP','%s: %s seq; %s bp (%.1fX); %.3f mean RQ.' % (self.getStr('Parent2'),rje.iStr(pseq[1]),rje.iStr(plen[1]),1.0*plen[1]/self.getInt('GenomeSize'),prq[1]/pseq[1]))
            ## ~ [1b] Extract Filtered Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,filtfile)
            SEQOUT = open(filtfile,'w')
            sx = 0.0; stot = 0; sn = len(seqbatch); fx = 0
            for seqlist in seqbatch:
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/9/0_3967 RQ=0.784
                si = 100.0/seqlist.seqNum(); stot += seqlist.seqNum()
                for seq in seqlist.seqs():
                    self.progLog('\r#OUT','Extracting subreads: %.2f%%' % (sx/sn)); sx += si
                    (name,sequence) = seqlist.getSeq(seq)
                    try: [smrt,zmw,pos,rq] = string.split(string.replace(name,'/',' '))
                    except:
                        [smrt,zmw,pos] = string.split(string.replace(name,'/',' '))
                        rq = minrq
                    if (cdb.data(smrt)['SMRT'],int(zmw),pos) not in zmwlist: continue
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fx += 1
            self.printLog('\r#OUT','Saved %s filtered subreads to %s.' % (rje.iStr(fx),filtfile))

            ### ~ [2] Summarise Filtered File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=T','seqin=%s' % filtfile,'autofilter=F']
            rje_seqlist.SeqList(self.log,seqcmd)

            return True
        except: self.errorLog('%s.run error' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: Diploidocus Class                                                                                #
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
    try: Diploidocus(mainlog,cmd_list).run()

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
