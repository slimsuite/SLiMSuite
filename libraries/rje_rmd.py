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
Module:       rje_rmd
Description:  R Markdown generation and execution module
Version:      0.1.0
Last Edit:    01/02/21
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

    The R markdown itself is run using:

    Rscript -e 'library(rmarkdown); rmarkdown::render("/path/to/test.Rmd", "html_document")'

    NOTE: The "html_document" over-rules the content of the file itself, e.g. "pdf_document" can turn it into a PDF
    rather than HTML.

    NOTE: Running the above generates a `*.html` file in the same place as the `*.Rmd` file (not the run directory).

    NOTE: For HTML output, R must be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

Commandline:

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added docHTML.
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
    # [Y] : Rscript -e 'library(rmarkdown); rmarkdown::render("/path/to/test.Rmd", "html_document")'
    # [ ] : Add logging of Rmd outputs
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_Rmd', '0.1.0', 'January 2021', '2019')
    description = 'R Markdown generation and execution module'
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
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
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
### SECTION II: Rmd Class                                                                                               #
#########################################################################################################################
class Rmd(rje_obj.RJE_Object):
    '''
    Rmd Class. Author: Rich Edwards (2019).

    Str:str
    
    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - CodeChunks = List of code chunk names to avoid duplication

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = []
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = ['CodeChunks']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        self.list['CodeChunks'] = []
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
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
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
            self.rmdTest()
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
    ### <3> ### Rmd Output Methods                                                                                      #
#########################################################################################################################
    def rmdKnit(self,rmdfile,document='html',stdout=False):  ### Knit Rmd to HTML/PDF file
        '''
        Knit Rmd to HTML/PDF file.
        >> rmdfile:str = R markdown file to knit
        >> document:str ['html'] = type of document to knit into
        << success:bool = whether output is generated
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.%s' % (rje.baseFile(rmdfile),document)
            rcmd = 'Rscript -e \'library(rmarkdown); rmarkdown::render("%s", "%s_document")\'' % (rmdfile,document)
            self.printLog('#RCMD',rcmd)
            rcmd += ' 2>&1'
            if self.v() < 2 and not stdout: os.popen(rcmd).read()
            else:
                self.progLog('#RCMD','Knitting %s...' % (rmdfile))
                os.system(rcmd)
            success = rje.exists(outfile)
            if success: self.printLog('#RCMD','%s generated from %s' % (outfile,rmdfile))
            else:
                self.printLog('#SYS','If pandoc error, try setting global variable: export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc')
                self.printLog('#SYS','If no pandoc error, check that required libraries in %s are installed' % rmdfile)
                raise IOError('%s not created' % outfile)
            return True
        except: self.errorLog('%s.rmdKnit error: check R installation' % self.prog()); return False
#########################################################################################################################
    def rmdOutput(self,rmdfile=None,header={},elements=[]):    ### Generate Rmd output file
        '''
        Generate Rmd output file. Call self.rmdKnit(rmdfile) to convert to another format.
        >> rmdfile:str [self.str['RmdFile']] = Full/relative path to Rmd output file.
        >> header:dict {} = Dictionary of Rmd header elements. Will default to self.log.obj['Info'] and HTML.
        >> elements:list [] = List of tuples (type,content) to output into file. For R code, content will be a dictionary

        This method puts together the basic elements of an R markdown file into a text document that can be knitted into
        HTML or PDF using self.rmdKnit().
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return

        except: self.errorLog('%s.rmdOutput error' % self.prog())
#########################################################################################################################
    def rmdHead(self,title=None,author=None,date=None,extra=[],setup=True,toc=True):    ### Generate Rmd output file
        '''
        Generate Rmd output file. Call self.rmdKnit(rmdfile) to convert to another format.
        >> rmdfile:str [self.str['RmdFile']] = Full/relative path to Rmd output file.
        >> header:dict {} = Dictionary of Rmd header elements. Will default to self.log.obj['Info'] and HTML.
        >> elements:list [] = List of tuples (type,content) to output into file. For R code, content will be a dictionary

        This method puts together the basic elements of an R markdown file into a text document that can be knitted into
        HTML or PDF using self.rmdKnit().
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcode = '---\n'
            if not title: title = '%s - %s' % (self.prog(),self.basefile())
            rcode += 'title: "%s"\n' % title
            if not author: author = self.log.obj['Info'].author
            rcode += 'author: "%s"\n' % author
            if not date: date = rje.dateTime(dateonly=True)
            rcode += 'date: "%s"\n' % date   #07/02/2019
            for (key,value) in extra:
                rcode += '%s: "%s"\n' % (key,value)
            rcode += 'output:\n  html_document:\n    css: http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css\n'
            if toc: rcode += rje.join(['    toc: true','    toc_float: true','    toc_collapsed: false','    toc_depth: 3','    number_sections: true',''],'\n')
            rcode += '---\n\n'

            if setup:
                rcode += '%s\n\n<a name="Top" />\n\n' % setupTest

            self.debug(rcode)
            return rcode
        except: self.errorLog('%s.rmdOutput error' % self.prog())
#########################################################################################################################
    def rmdTable(self,delimfile=None,name='dbtable',codesection=True,loadtable=True,showtable=True,delim='tab',kable=None,rows=10,cols=10):  ### Output table
        '''
        Output table. If the table is larger than rows tall, or cols wide, paged_table will be used. Otherwise, kable
        will be used. This can be over-ridden by setting kable=True, or kable=False (for paged_table).
        :param delimfile: delimited text file
        :param name: name for data.frame (and code section if needed)
        :param codesection: give R code wrapping text
        :param loadtable: load table into R object
        :param showtable: display the table
        :param delimit: tab/csv
        :return: text of Rmd code chunk
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcode = '# Load and display %s\n' % name
            # Read table
            if loadtable:
                rcode = '# Load and display %s\n' % delimfile
                if delim == 'csv':
                    rcode += '%s <- read.csv("%s", header = TRUE, stringsAsFactors = FALSE, comment.char = "")\n' % (name,delimfile)
                else:
                    rcode += '%s <- read.delim("%s", header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)\n' % (name,delimfile)
            # Show table
            if showtable:
                if kable == None:
                    rcode += 'if(nrow(%s) > %d | ncol(%s) > %d){\n' % (name,rows,name,cols)
                    rcode += '    rmarkdown::paged_table(%s)\n' % (name)
                    rcode += '}else{\n'
                    rcode += '    knitr::kable(%s, row.names = FALSE)\n' % name
                    rcode += '}\n'
                elif kable:
                    rcode += 'knitr::kable(%s, row.names = FALSE)\n' % name
                else:
                    #rcode += 'rmarkdown::paged_table(%s, options = list(rows.print = %d, max.print = %d, cols.print = %d, rownames.print = FALSE))\n' % (name,rows,max,cols)
                    #rcode += 'rmarkdown::paged_table(%s, options = list(rows.print = %d))\n' % (name,rows)
                    rcode += 'paged_table(%s)\n' % (name)
            # Code section
            if codesection:
                codename = name
                if codename in self.list['CodeChunks']:
                    i = 1
                    while '%s%d' % (name,i) in self.list['CodeChunks']: i += 1
                    codename = '%s%d' % (name,i)
                self.list['CodeChunks'].append(codename)
                rcode = '```{r %s, echo=FALSE}\n%s\n```\n\n' % (codename,rcode)
            # Return text
            return rcode
        except: self.errorLog('%s.rmdOutput error' % self.prog())
#########################################################################################################################
    def rmdTest(self):  ### Generates a test Rmd File
        '''Generates a test Rmd File.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rmdfile = self.basefile() + '.Rmd'
            RMD = open(rmdfile,'w')
            RMD.write(self.rmdHead())
            #RMD.write(setupTest)
            RMD.write(mdTest)
            RMD.write(rcodeTest)
            RMD.write('## Tables\n\n')
            for tdtfile in glob.glob('*.tdt'):
                RMD.write('```\n%s\n```\n\n%s\n\n' % (tdtfile,self.rmdTable(tdtfile)))
            RMD.write(htmlTest)
            RMD.close()
            self.rmdKnit(rmdfile)
        except: self.errorLog('%s.rmdOutput error' % self.prog())
#########################################################################################################################
### End of SECTION II: Rmd Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def docHTML(self):  ### Generate Rmd and HTML documents from main run() method docstring.
    '''Generate Rmd and HTML documents from main run() method docstring.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = self.log.obj['Info']
        if not self.getStrLC('Basefile'): self.baseFile(info.program.lower())
        prog = '%s V%s' % (info.program,info.version)
        rmd = Rmd(self.log,self.cmd_list)
        rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
        #!# Replace this with documentation text?
        rtxt += rje.replace(self.run.__doc__,'\n        ','\n')
        rtxt += '\n\n<br>\n<small>&copy; 2023 Richard Edwards | rich.edwards@uwa.edu.au</small>\n'
        rmdfile = '%s.docs.Rmd' % self.baseFile()
        open(rmdfile,'w').write(rtxt)
        self.printLog('#RMD','RMarkdown %s documentation output to %s' % (prog,rmdfile))
        rmd.rmdKnit(rmdfile)
    except:
        self.errorLog(self.zen())
        raise   # Delete this if method error not terrible
#########################################################################################################################
headTest = '''---
title: "RJE_RMD"
author: "Rich Edwards"
date: "07/02/2019"
output: html_document
---
'''
#########################################################################################################################
setupTest = '''
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
'''
#########################################################################################################################
mdTest = '''
## R Markdown

<a name="aname">?</a>

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

Trying some internal links:

* [cars](#cars)
* [head2](#rmarkdown)
* [aname](#aname)
'''
#########################################################################################################################
rcodeTest = '''
```{r cars}
summary(cars)
```
'''
#########################################################################################################################
htmlTest = '''
<hr>
<small>&copy; Richard Edwards 2019</small>

'''
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
    try: Rmd(mainlog,['basefile=test']+cmd_list).run()

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
