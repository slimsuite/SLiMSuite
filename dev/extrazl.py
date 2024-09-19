#!/usr/bin/python

# See below for name and description
# Copyright (C) 2023 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       Extrazl
Description:  Extraction of ZMW longest subreads
Version:      0.1.0
Last Edit:    26/05/23
Copyright (C) 2023  Richard J. Edwards - See source code for GNU License Notice

Function:
    Extrazl is a simple tool for extracting the longest subread of each ZMW from a fastq file, with an optional minimum
    length filter. By default, ZMWs with over three subreads will not be output, as these should be used for CCS reads
    and those not producing CCS reads probably have quality issues.

    For the bamstats function with BAM input, samtools must be installed on the system. This will parse a collection of
    statistics from a mapped BAM file to assist with differentiating good and bad subreads with a view to adding some
    additional filters in future.

Commandline:
    ### ~ Main Extrazl run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence fastq file (may be zipped) [None]
    basefile=FILE   : Root of output file names [$SEQINBASE]
    seqout=FILE     : Output fastq file name [$BASEFILE.extrazl.fastq]
    minlen=INT      : Optional minimum length filter [0]
    maxpass=INT     : Maximumn number of passes detected in a ZMW to keep (0 = no filter) [3]
    zstats=T/F      : Whether to collect and output some stats per ZMW during processing [False]
    bamstats=BAM    : Extract mapping stats from a SAM/BAM file of mapped reads to use for filtering [None]
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
import rje, rje_db, rje_obj, rje_samtools, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.0.1 - Fixed minlen filter bug.
    # 0.1.0 - Added maxpass=INT : Maximumn number of passes detected in a ZMW to keep [3]; Add zstats=T/F.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add summarise and gzipping options for output.
    # [ ] : Add bamstats parsing of the mapped subreads BAM file.
    # [ ] : Add input file checks.
    # [ ] : Add recongition of BAM input for seqin and include full processing - check speed versus python. (pyparse=T/F)
    # [ ] : Add multithreaded processing - process chunks in parallel into temp directory.
    # [ ] : Enable fasta input/output (fastx=fasta or auto-recognition)
    # [ ] : Rework to use two passes? (Identify the subreads to output, then output the subreads.)
    # [ ] : Add warning that the zstats processing is very slow (due to the compression scoring).
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Extrazl', '0.1.0', 'May 2023', '2023')
    description = 'Extraction of ZMW longest subreads'
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
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
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
### SECTION II: Extrazl Class                                                                                           #
#########################################################################################################################
class Extrazl(rje_obj.RJE_Object):
    '''
    Extrazl Class. Author: Rich Edwards (2021).

    Str:str
    - BAMStats=BAM    : Extract mapping stats from a SAM/BAM file of mapped reads to use for filtering [None]
    - SeqIn=FILE      : Input sequence fastq file (may be zipped) [None]
    - SeqOut=FILE     : Output fastq file name [$BASEFILE.extrazl.fastq]

    Bool:boolean
    - ZStats=T/F      : Whether to collect and output some stats per ZMW during processing [False]

    Int:integer
    - MaxPass=INT     : Maximumn number of passes detected in a ZMW to keep [3]
    - MinLen=INT      : Optional minimum length filter [0]

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
        self.strlist = ['BAMStats','SeqIn','SeqOut']
        self.boollist = ['ZStats']
        self.intlist = ['MaxPass','MinLen']
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'ZStats':False})
        self.setInt({'MaxPass':3})
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
                self._cmdReadList(cmd,'file',['BAMStats','SeqIn','SeqOut'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['ZStats'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxPass','MinLen'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'perc',['Att'])  # Percentage, converts to 1-100 scale.
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
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.extrazl()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log, self.cmd_list + ['tuplekeys=T'])
            if self.getStrLC('BAMStats'):
                raise ValueError('BAMStats mode not yet implemented.')
            if not self.getStrLC('SeqIn'): raise IOError('Need to provide seqin=FILE')
            self.printLog('#SEQIN', 'Input fastq file: %s' % self.getStr('SeqIn'))
            if not self.baseFile(return_none=''): self.baseFile(rje.baseFile(self.getStr('SeqIn')))
            self.printLog('#BASE', 'Output file basename: %s' % self.baseFile())
            if not self.getStrLC('SeqOut'): self.setStr({'SeqOut':'{0}.extrazl.fastq'.format(self.baseFile())})
            self.printLog('#SEQOUT', 'Output fastq file: %s' % self.getStr('SeqOut'))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Subread Extraction Methods                                                                              #
#########################################################################################################################
    def extrazl(self):  ### Parses data from SeqIn and outputs longest subreads into SeqOut.
        '''
        Parses data from SeqIn and outputs longest subreads into SeqOut.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.getStr('SeqIn')
            seqout = self.getStr('SeqOut')
            rje.backup(self,seqout)
            SEQOUT = open(seqout,'w')
            if seqin.endswith('.gz'): SEQIN = os.popen('cat {0} | zcat'.format(seqin))
            else: SEQIN = open(seqin,'r')
            zdb = self.getBool('ZStats')
            if zdb:
                zdb = self.db().addEmptyTable('zstats',['zmw','rn','minlen','maxlen','a','c','g','t','zlib'],['zmw'],log=True)

            ### ~ [2] Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            zmw = ''
            z = {'zmw':'','rn':0,'minlen':0,'maxlen':0,'a':0,'c':0,'g':0,'t':0,'zlib':0}
            zlen = []
            inx = 0; outx = 0; zx = 0; rx = 0; passfiltx = 0
            subread = ()  # (name,sequence,quality)
            self.progLog('\r#FASTQ','Parsing ZMW from fastq: {0} subreads ({1} ZMW) -> {2} output'.format(rje.iStr(inx),rje.iStr(zx),rje.iStr(outx)))
            while SEQIN:
                self.progLog('\r#FASTQ','Parsing ZMW from fastq: {0} subreads ({1} ZMW) -> {2} output'.format(rje.iStr(inx),rje.iStr(zx),rje.iStr(outx)),rand=0.01)
                line = rje.chomp(SEQIN.readline())
                if line == "":
                    name = ""
                    sequence = ""
                    quality = ""
                    zmw = ''
                    zlen.append(0)
                elif line[:1] != '@': raise ValueError('Problem with fastq file format!')
                else:
                    #i# Read sequence
                    name = line[1:]
                    sequence = rje.chomp(SEQIN.readline())
                    SEQIN.readline()
                    quality = rje.chomp(SEQIN.readline())
                    if len(sequence) != len(quality): raise ValueError('Problem with fastq file format!')
                    inx += 1; rx += 1
                    zlen.append(len(sequence))
                #i# Compare to previous ZMW and output or update
                prevzmw = zmw
                if name:
                    try:
                        [zmw,zi,zj] = rje.matchExp('^(\S+)\/(\d+)_(\d+)',name)
                        if len(sequence) != (int(zj) - int(zi)):
                            self.warnLog('{0} -> {1} vs {2} length mismatch!'.format(name,rje.iStr(int(zj) - int(zi)),rje.iStr(len(sequence))))
                    except:
                        self.warnLog('Problem with sequence name: {0}'.format(name))

                if prevzmw and zmw != prevzmw:
                    zx += 1
                    if zdb:
                        z = {'zmw': zmw, 'rn': rx, 'minlen': min(zlen[:-1]), 'maxlen': len(subread[1]), 'zlib': self.compressionScore(subread[1])}
                        for n in 'acgt':
                            z[n] = subread[1].lower().count(n)
                        zdb.addEntry(z)
                    if len(subread[1]) >= self.getInt('MinLen'):
                        if rx > self.getInt('MaxPass') > 0:
                            passfiltx += 1
                            self.debug('Skipping {0} ({1} passes)'.format(subread[0], rx))
                        else:
                            self.debug('{0} -> {1}'.format(subread[0],rje_seqlist.dnaLen(len(subread[1]))))
                            SEQOUT.write('@{0}\n{1}\n+\n{2}\n'.format(subread[0],subread[1],subread[2]))
                            outx += 1
                    rx = 0
                if line == "": break
                if not prevzmw or prevzmw != zmw or len(sequence) > len(subread[1]):
                    subread = (name, sequence, quality)
                    self.bugPrint('{0} = {1}'.format(name, rje_seqlist.dnaLen(len(sequence))))

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#FASTQ','Parsed ZMW subreads from fastq: {0} subreads ({1} ZMW) -> {2} output'.format(rje.iStr(inx),rje.iStr(zx),rje.iStr(outx)))
            if self.getInt('MaxPass') >= 0:
                self.printLog('\r#FILT','Filtered {0} ZMW subreads that met length criteria but exceeded maxpass={1}'.format(rje.iStr(passfiltx),self.getInt('MaxPass')))
            else:
                self.printLog('\r#FILT','No maxpass filter (maxpass={0})'.format(self.getInt('MaxPass')))
            SEQIN.close()
            SEQOUT.close()
            self.printLog('\r#SEQOUT','{0} longest ZMW reads output to: {1}'.format(rje.iStr(outx),seqout))
            if zdb:
                zdb.saveToFile()
            return True
        except:
            try: SEQIN.close()
            except: pass
            try: SEQOUT.close()
            except: pass
            self.errorLog('%s.extrazl error' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: Extrazl Class                                                                                    #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: Extrazl(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
