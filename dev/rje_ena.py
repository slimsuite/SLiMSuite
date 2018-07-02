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
Module:       rje_ena
Description:  ENA Submission XML Generator
Version:      0.0.0
Last Edit:    15/02/18
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    This script is designed to read directory contents and ask questions to generate XML files required from programmatic
    ENA uploads.

Commandline:
    datatype=X      : Type of data in directory [PacBio]
    dirlist=LIST    : List of directories to process (blank for all) []

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
import rje, rje_obj
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
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('ENA', '0.0.0', 'February 2018', '2016')
    description = 'ENA Submission XML Generator'
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
### SECTION II: ENA Class                                                                                               #
#########################################################################################################################
class ENA(rje_obj.RJE_Object):
    '''
    ENA Class. Author: Rich Edwards (2018).

    Str:str
    - DataType=X   : Type of data in directory [PacBio]

    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - DirList=LIST    : List of directories to process (blank for all) []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['DataType']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = ['DirList']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'DataType':'pacbio'})
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
                self._cmdReadList(cmd,'str',['DataType'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['DirList']) # List of files using wildcards and glob
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
            self.runXML()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        There is currently no specific help available on REST output for this program.
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
    def runXML(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datatype = self.getStr('DataType').lower()
            exps = {}        # Experiment alias: run directory list
            runs = {}        # Run alias: file list
            run2run = {}    # Convert runs to run aliases

            ## ~ [1a] Get Files and Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dirlist = rje.listDir(self,subfolders=False,folders=True,files=False,summary=True)
            filelist = rje.listDir(self,folder=os.getcwd(),subfolders=True,folders=False,files=True,summary=True)

            current = os.getcwd()
            curlen = len(current) + 1
            for pathlist in [dirlist,filelist]:
                for i in range(len(pathlist)):
                    path = pathlist[i]
                    if path.startswith(current): pathlist[i] = path[curlen:]
                    else: raise ValueError(path)

            ## ~ [1b] Clean up files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dx = len(dirlist)
            self.debug(filelist)
            self.debug(dirlist)
            self.debug(self.list['DirList'])
            if self.list['DirList']:
                dirlist = rje.listIntersect(dirlist,self.list['DirList'])
            self.printLog('#DIR','Process %d of %d directories' % (len(dirlist),dx))
            keepext = []
            if datatype == 'pacbio': keepext = ['.h5','.xml']
            for filename in filelist[0:]:
                ext = os.path.splitext(filename)[1]
                if len(string.split(filename,os.sep)) < 2 or ext not in keepext or string.split(filename,os.sep)[0] not in dirlist:
                    filelist.remove(filename)
            self.printLog('#FILES','%s files kept from %s directories' % (rje.iLen(filelist),rje.iLen(dirlist)))
            self.debug(filelist[:10])
            self.debug(filelist[-10:])

            ### ~ [2] Parse runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if datatype == 'pacbio':
                for filename in filelist[0:]:
                    self.printLog('#FILE',filename)
                    filedata = string.split(filename,os.sep)
                    parent = filedata[0]    # This should be a directory containing runs
                    experiment = filedata[1]
                    expalias = string.join(string.split(experiment,'.')[:2],'.')
                    run = string.join(filedata[1:3],os.sep)
                    if expalias not in exps: exps[expalias] = []
                    if run not in exps[expalias]: exps[expalias].append(run)
                    runalias = '%s-%d' % (expalias,len(exps[expalias]))
                    run2run[run] = runalias
                    runfile = filedata[-1]
                    if runalias not in runs: runs[runalias] = []
                    runs[runalias].append(filename)
                    self.printLog('#PARSE','%s - %s: (%d) %s' % (expalias,runalias,len(runs[runalias]),filename))

            ### ~ [3] Generate XML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Experiment XML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            efile = '%s.exp.xml' % self.baseFile()
            elines = ['<?xml version="1.0" encoding="UTF-8"?>','<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.experiment.xsd">']
            for experiment in rje.sortKeys(exps):
                ex = experiment[-1:]
                elines += ['    <EXPERIMENT alias="%s" center_name="">' % experiment,
                            '        <TITLE>Cane toad whole genome sequencing - PacBio library %s</TITLE>' % ex,
                            '        <STUDY_REF accession="ERP106543"/>',
                            '        <DESIGN>',
                            '            <DESIGN_DESCRIPTION/>',
                            '            <SAMPLE_DESCRIPTOR accession="ERS2169570"/>',
                '            <LIBRARY_DESCRIPTOR>',
                '                <LIBRARY_NAME>%s</LIBRARY_NAME>' % experiment,
                '                <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>',
                '                <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>',
                '                <LIBRARY_SELECTION>size fractionation</LIBRARY_SELECTION>',
                '                <LIBRARY_LAYOUT>',
                '                    <SINGLE/>',
                '                </LIBRARY_LAYOUT>',
                '                <LIBRARY_CONSTRUCTION_PROTOCOL></LIBRARY_CONSTRUCTION_PROTOCOL>',
                '               </LIBRARY_DESCRIPTOR>',
                '        </DESIGN>',
                '        <PLATFORM>',
                '            <PACBIO_SMRT>',
                '                <INSTRUMENT_MODEL>PacBio RS II</INSTRUMENT_MODEL>',
                '            </PACBIO_SMRT>',
                '        </PLATFORM>',
                '        <EXPERIMENT_ATTRIBUTES>',
                '            <EXPERIMENT_ATTRIBUTE>',
                '                <TAG>Size selection</TAG>',
                '                <VALUE>15-50 kb</VALUE>',
                '            </EXPERIMENT_ATTRIBUTE>',
                '            <EXPERIMENT_ATTRIBUTE>',
                '                <TAG>Sequencing Chemistry</TAG>',
                '                <VALUE>P6C4</VALUE>',
                '            </EXPERIMENT_ATTRIBUTE>',
                '        </EXPERIMENT_ATTRIBUTES>',
                '    </EXPERIMENT>']
            elines += ['</EXPERIMENT_SET>']
            open(efile,'w').write(string.join(elines,'\n'))
            self.printLog('#EXP','Experiment data saved to %s' % efile)
            ## ~ [3b] Run XML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rfile = '%s.run.xml' % self.baseFile()
            rlines = ['<?xml version="1.0" encoding="UTF-8"?>','<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.run.xsd">']
            for experiment in rje.sortKeys(exps):
                for run in exps[experiment]:
                    runalias = run2run[run]
                    rlines += ['   <RUN alias="%s" center_name="">' % runalias,'    <EXPERIMENT_REF refname="%s"/>' % experiment,
                            '      <DATA_BLOCK>','        <FILES>']
                    for filename in runs[runalias]:
                        rlines += ['             <FILE filename="%s" filetype="PacBio_HDF5">' % filename,'             </FILE>']
                    rlines += ['        </FILES>','      </DATA_BLOCK>','  </RUN>']
            rlines += ['</RUN_SET>']
            open(rfile,'w').write(string.join(rlines,'\n'))
            self.printLog('#RUN','Run data saved to %s' % rfile)
            ## ~ [3c] Submission XML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xfile = '%s.xml' % self.baseFile()
            xlines = ['<?xml version="1.0" encoding="UTF-8"?>','<SUBMISSION alias="%s" center_name="">' % self.baseFile(),
                        '   <ACTIONS>','      <ACTION>','         <ADD source="%s" schema="experiment"/>' % efile,
                        '      </ACTION>','      <ACTION>','         <ADD source="%s" schema="run"/>' % rfile,
                        '      </ACTION>','   </ACTIONS>','</SUBMISSION>']
            open(xfile,'w').write(string.join(xlines,'\n'))
            self.printLog('#SUBXML','Submission XML saved to %s' % xfile)
            return
        except: self.errorLog('%s.method error' % self.prog())
#########################################################################################################################
### End of SECTION II: ENA Class                                                                                        #
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
    try: ENA(mainlog,cmd_list).run()

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
