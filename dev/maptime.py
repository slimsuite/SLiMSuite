#!/usr/bin/python

# See below for name and description
# Copyright (C) 2012 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       MapTime
Description:  MapTime TimeLine and TimePoint format conversion
Version:      0.1
Last Edit:    12/11/12
Copyright (C) 2012  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a very simple MapTime file conversion program. Currently, three different file formats can be inter-
    converted:
    - Delimited text, either tab (*.tdt) or comma (*.csv) separated.
    - Plain descriptive text (*.txt) that could be used as input for rje_glossary.py (TimePoints only)
    - Database input lists for direct entry to MapTime database (*.db).
    
    Delimited text files have the following headers:    
    - TimePoint Name    TimePoint Description   Source URL  Year    yearUnit    month   day 
    - keyword1  keyword2    keyword3    keyword4    keyword5

Commandline:

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Basic file conversion functionality.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add TimeLines and Keywords to Database tables.
    # [ ] : Add glossary text output.
    # [ ] : Add glossary HTML output using rje_glossary.py. (Will need to move to libraries.)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('MapTime', '0.1', 'November 2012', '2012')
    description = 'MapTime TimeLine and TimePoint format conversion'
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
### SECTION II: MapTime Class                                                                                           #
#########################################################################################################################
class MapTime(rje_obj.RJE_Object):     
    '''
    MapTime Class. Author: Rich Edwards (2012).

    Str:str
    - FileIn = Name of input file with MapTime data (auto-recognise format) [None]
    - FileOut = Name of output file [None]
    - Format = Format of output file [None]
    
    Bool:boolean

    Int:integer

    Num:float
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database object for storing TimeLine and TimePoints etc.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FileIn','FileOut','Format']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({})
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
                self._cmdReadList(cmd,'str',['Format'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['FileIn','FileOut'])  # String representing file path 
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
            self.loadTimePoints(self.getStr('FileIn'))
            self.saveTimePoints(self.getStr('FileOut'),self.getStr('Format'))
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = db = rje_db.Database(self.log,self.cmd_list)
            db.addEmptyTable('TimePoints',['TimePoint Name','TimePoint Description','Source URL','Year','yearUnit','month','day','keyword1','keyword2','keyword3','keyword4','keyword5'],keys=['TimePoint Name'])         
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Loading and Saving Methods                                                                              #
#########################################################################################################################
    def loadTimePoints(self,filename):  ### Load TimePoints from file of various formats
        '''Load TimePoints from file of various formats.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(filename): return self.errorLog('File %s missing!' % filename)
            data = open(filename,'r').readlines()
            db = self.db('TimePoints')
            
            ### ~ [2] Load from File Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Delimited File Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.split(data[0])[0] == 'TimePoint Name':    #
                ftype = 'delimited text file'
                temp = self.db().addTable(filename,mainkeys=['TimePoint Name'],name='temp')
                for entry in temp.entries(): db.addEntry(entry)
                db.deleteTable(temp)
            ## ~ [2b] File of Database Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif data[0][0] == '(':
                ftype = 'database string'
                for line in data:
                    line = rje.chomp(line)
                    while line[-1:] == ' ': line = line[:-1]
                    pdata = rje.split(rje.replace(line[2:-3],', ',','),"','")
                    if not pdata: continue
                    if rje.matchExp('^(\d+)$',pdata[0]): pdata.pop(0)   # Database output with key ID numbers
                    entry = {}
                    for field in db.fields(): entry[field] = pdata[db.fields().index(field)]
                    db.addEntry(entry)
            ## ~ [2c] Glossary Text File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                ftype = 'glossary text file'
                for line in data:
                    if '(TimePoint)' not in line: continue
                    # American Independence. (TimePoint) 1776 AD, 4 July. The US declared independence from the British Empire. Source: <http://en.wikipedia.org/wiki/United_States_Declaration_of_Independence>[Wikipedia]. (Keywords: history)
                    pdata = rje.split(line,'. ')
                    if pdata[2][-2:] == 'ya':
                        pdata[1] = '%s. %s' % (pdata[1],pdata.pop(2))
                    entry = {'TimePoint Name':pdata[0]}
                    try: entry['Source URL'] = rje.matchExp('Source: <(\S+)>',line)[0]
                    except: self.errorLog('Cannot read Source URL')
                    try: entry['TimePoint Description'] = rje.matchExp('^(\S.+\S) Source: <',rje.join(pdata[2:],'. '))[0]
                    except: self.errorLog('Cannot read TimePoint Description: %s' % line)
                    if pdata[1][-2:] == 'ya':
                        [entry['Year'],entry['yearUnit']] = rje.split(pdata[1])[-2:]
                    else:
                        try:
                            ydata = rje.matchExp('(\d+) (\S+), (\d+) (\S+)$',pdata[1])
                            if ydata:
                                for i in range(4): entry[['Year','yearUnit','month','day'][i]] = ydata[i]   
                            else: (entry['Year'],entry['yearUnit']) = rje.matchExp('(\d+) (\S+)$',pdata[1])
                        except: self.errorLog('Cannot parse time from %s' % pdata[1])
                    kfield = ['keyword1','keyword2','keyword3','keyword4','keyword5']
                    try: 
                        keywords = rje.split(rje.matchExp('\(Keywords: (\S.+)\)',pdata[-1])[0],', ')
                        while keywords and kfield:
                            entry[kfield.pop(0)] = keywords.pop(0)
                        while kfield: entry[kfield.pop(0)] = 'blank'
                        if keywords: self.printLog('#ERR','%d extra Keywords (%s)!' % (len(keywords),rje.join(keywords,', ')))
                    except: self.errorLog('Cannot read Keywords (%s)' % pdata[-1])
                    db.addEntry(entry)
            ### ~ [3] Summarise Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#TP','Timepoints read from %s: %s TimePoints total.' % (ftype,db.entryNum()))
            return True
        except: self.errorLog('%s.loadTimePoints(%s) error' % (self,filename)); return False
#########################################################################################################################
    def saveTimePoints(self,filename='',format='tdt',entries=[]):   ### Saves TimePoints to a file
        '''
        Saves TimePoints to a file from main TimePoints table.
        >> filename:str [''] = Output filename. Will use basefile if none given.
        >> format:str ['tdt'] = Output file format (csv/tsv/txt/db)
        >> entries:list [] = Entries from main table to output. (All if none given).
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db('TimePoints')
            if format.lower() in ['','none']: format = rje.split(filename.lower(),'.')[-1]
            if not filename: filename = '%s.%s' % (self.basefile(),format)
            if not entries: entries = db.entries()
            ### ~ [2] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Simple delimited file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if format in ['csv','tdt']: 
                self.blanksToEmpty()
                rje.delimitedFileOutput(self,filename,db.fields(),rje_backup=True)
                for entry in entries: rje.delimitedFileOutput(self,filename,db.fields(),datadict=entry)
            ## ~ [2b] Text file output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                self.emptyToBlank()
                rje.backup(self,filename)
                OUT = open(filename,'a')
                for entry in entries:
                    if format == 'db':
                        outlist = []
                        for field in db.fields(): outlist.append(entry[field])
                        out_txt = '%s' % outlist
                        OUT.write('(%s);\n' % out_txt[1:-1])
                    else:
                        # American Independence. (TimePoint) 1776 AD, 4 July. The US declared independence from the British Empire. Source: <http://en.wikipedia.org/wiki/United_States_Declaration_of_Independence>[Wikipedia]. (Keywords: history)
                        out_text = '%s. (TimePoint) ' % entry['TimePoint Name']
                        if entry['month'] in ['','blank']: out_text += '%s %s.' % (entry['Year'],entry['yearUnit'])
                        else: out_text += '%s %s, %s %s.' % (entry['Year'],entry['yearUnit'],entry['month'],entry['day'])
                        out_text = '%s %s Source: <%s>[%s].' % (out_text,entry['TimePoint Description'],entry['Source URL'],entry['Source URL'])
                        klist = []
                        for i in range(1,6):
                            if entry['keyword%d' % i] not in ['','blank']: klist.append(entry['keyword%d' % i])
                        out_text = '%s (Keywords: %s)' % (out_text,rje.join(klist,', '))
                        OUT.write('%s\n' % out_text)
            self.printLog('#OUT','%d entries output to %s' % (len(entries),filename))
        except: self.errorLog('%s.saveTimePoints(%s) error' % (self,filename)); return False
#########################################################################################################################
    ### <4> ### TimePoint Data Editing Methods                                                                          #
#########################################################################################################################
    def blanksToEmpty(self):    ### Replace 'blank' values with empty values
        '''Replace 'blank' values with empty values.'''
        db = self.db('TimePoints'); bx = 0
        for entry in db.entries():
            for field in db.fields():
                if entry[field] == 'blank': entry[field] = ''; bx += 1
        self.printLog('#DB','%s blank values represented with empty values' % rje.iStr(bx))
#########################################################################################################################
    def emptyToBlank(self):     ### Replace empty values with 'blank' values
        '''Replace empty values with 'blank' values.'''
        db = self.db('TimePoints'); bx = 0
        for entry in db.entries():
            for field in db.fields():
                if entry[field] == '': entry[field] = 'blank'; bx += 1
        self.printLog('#DB','%s empty values represented with blank values' % rje.iStr(bx))
#########################################################################################################################
### End of SECTION II: MapTime Class                                                                                    #
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
    try: MapTime(mainlog,cmd_list).run()

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
