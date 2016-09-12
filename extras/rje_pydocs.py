#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <seqsuite@gmail.com>
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
# Author contact: <seqsuite@gmail.com> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_pydocs
Description:  Python Module Documentation & Distribution
Version:      2.16.3
Last Edit:    03/12/15
Copyright (C) 2011 Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to aid handling and maintenance of linked Python Modules by creating limited documentation
    for the module, its objects and their methods. Possible interactions between modules and methods can also be
    identified. In addition, this module can be used to make the distribution directory for a python project and form the
    basis of linked webpages.    

Commandline:
    ### ~~~ Python Module Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    name=X          : Name for PyDoc run. Used for file naming and within documentation files. ['pydocs']
    pylist=LIST     : List of python modules to upload. Can have * wildcards. Will add '.py' if missing. ['*']
    pypath=PATH     : Path to python modules. Will also look in listed sourcedir subfolders. [../]
    sourcedir=LIST  : List of subdirectories in which to look for modules [tools,extras,libraries,legacy]
    addimports=T/F  : Whether to add identified imported modules to python module list [True]
    docsource=PATH  : Input path for Python Module documentation (manuals etc.) ['../docs/']

    ### ~~~ Python Module Processing and Documentation Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    html=T/F        : Whether to make a basic HTML page of module docstrings (will make linked fun later) [False]
    fulldoc=T/F     : Whether to generate full docstring output including Classes and Methods [False]
    methodskip=LIST : List of methods to skip documentation [makeInfo,cmdHelp,setupProgram,Main]
    methodcap=X     : Maximum number of method calls before collapsed to single line [0]
    calls=T/F       : Whether to output Method Calls [False]
    self=T/F        : Whether to include 'self' calls of methods if calls=T [False]
    docdir=PATH     : Output path for Python Module documentation ['../docs/']
    stylesheets=LIST: List of style sheets to use for HTML ['rje_tabber.css','re1u06.css']
    stylepath=PATH  : Path to style sheets and javascript code folder etc. [http://www.southampton.ac.uk/~re1u06/]
    author=X        : Author name to put at bottom of webpages [RJ Edwards]
    keywords=LIST   : List of keywords to add in additon to module names/descriptions []
    modlinks=LIST   : List of rje_ppi formats for module import links (xgmml,tdt,svg,r,png) [xgmml]
    makepages=T/F   : Special run to generate default cmd/ and docs/ pages for commandline option docs [False]
    logourl=URL     : URL to SLiMSuite program logos ['http://www.slimsuite.unsw.edu.au/graphics/']
    manualurl=URL   : URL to SLiMSuite program manuals ['http://docs.slimsuite.unsw.edu.au/software/slimsuite/docs/manuals/']
    resturl=URL     : URL for REST server to pull outfmt docs ['http://rest.slimsuite.unsw.edu.au/']

    ### ~~~ Distribution and Webpage options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###    
    distribute=LIST : Names of distributions - gets details from distributions.txt ['SLiMSuite']
    distdir=PATH    : Output directory for distribution directories ['../packages/']
    webdir=PATH     : Directory containing webpage resources [../html/']
    email=X         : E-Mail address for general contact [seqsuite@gmail.com]
    release=X       : Release for packages [YYYY-MM-DD]
    
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, re, shutil, string, sys, time, urllib2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_html, rje_obj, rje_ppi, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 2.0 - Initial Compilation based on rje_pydocs 1.7 (now rje_pydocs_V1.py), converted to new module framework.
    # 2.1 - Added R code and INI directory.
    # 2.2 - Minor bug fixes for distribution making code.
    # 2.3 - Added Release text (default YYMMDD) to be inserted into downloads.
    # 2.4 - Implemented reading of version numbers and outputting updates to log.
    # 2.5 - Converted to make use of settings/ directory in place of libraries/ini/
    # 2.6 - Minor bug fixes for readme etc. Changed to read ini for packages from ./defaults/.
    # 2.7 - Added rje_ppi output for module links.
    # 2.8 - Added parsing of commandline options from docstring and cmdRead calls.
    # 2.8 - Added docsource=PATH  : Input path for Python Module documentation (manuals etc.) ['../docs/']
    # 2.9 - Attempts to fix some broken links and sort out manuals confusion.
    # 2.10- Added legacy/ to default drives for processing.
    # 2.11- Added development option for multifile HTML readme output.
    # 2.12- Continued modification of the HTML output for SLiMSuite package generation.
    # 2.13- Centred logo and altered restricted width (300) to height (120).
    # 2.14.0 - Added handling of X.Y.Z versions.
    # 2.15.0 - Added parsing and generation of "pages" for new rest server docs functions.
    # 2.15.1 - Tweaked formatting of outfmt and docstring documentation.
    # 2.15.2 - Tweaked formatting of docstring documentation.
    # 2.15.3 - Fixed URL formatting of docstring documentation.
    # 2.16.0 - Added Webserver tab to doc parsing from settings/*.form.
    # 2.16.1 - Added parsing of imports within a try/except block. (Cannot be on same line as try: or except:)
    # 2.16.2 - Tweaked makePages() output.
    # 2.16.3 - Fixed docstring REST parsing to work with _V* modules.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Reading in of Python modules to new Classes using rje_db Database Objects and Tables.
    # [N] : Add an ignore list of subfolders to ignore
    # [Y] : Improved output of documentation to html and text files.
    # [Y] : Distributions with new file setup (seqsuite, slimsuite and rjesuite)
    # [Y] : Add manuals to lists of entries to read from distributions.txt
    # [Y] : Add GNU license etc. automatically to output files. (Also style sheets and javascript.)
    # [?] : Update web defaults
    # [Y] : New webpage output with tabs and improved sidebar read from other file.
    # [?] : Add a parameter to ignore listed methods (verbose etc.) for methodcalls
    # [Y] : Add author to webpage generation
    # [X] : Automatically copy distribution files
    # [Y] : Add a general docs folder to the packages directory (for distributions) to cleanup web copying.
    # [Y] : Add reading in of old distributions and output report of updated modules etc.
    # [Y] : Add the settings/ directory to packages and include an rje.ini file with all defaults commented.
    # [Y] : Test ppi links output and remove dev.
    # [ ] : Add a table of descriptions for the "other files"
    # [ ] : Add one line summaries to either files themselves or another file to read in.
    # [ ] : Add reading of a Publications tab, much like the package summary tab. ("X_pub.html")
    # [ ] : Check that reading old histories is working OK. Seems to be returning different things for each package.
    # [Y] : Add parsing and compiling of commandline options.
    # [ ] : Add reading of markdown text for generating readme files - also change how docstring is handled.
    # [ ] : Tidy up docs in general: docs_pdf and docs_logos etc. Add a DocSource attribute, separate from DocDir.
    # [ ] : Modify/Control tab usage in HTML. (Multiple tiers?)
    # [ ] : Tidy up use of logos. (HAQESAC is giant!)
    # [Y] : Make a tabber=F option with alternative readme output. -> onefile option in readme output.
    # [X] : Add update tab to individual modules.
    # [ ] : Add exclusion list for module parsing?
    # [Y] : Replace the rje.ini file from the packages with defaults.ini.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('RJE_PYDOCS', '2.16.3', 'December 2015', '2011')
    description = 'Python Module Documentation & Distribution'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
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
web_defaults = {'Logo':'../_Documentation/binflogo.jpg',                    # Logo image file
                'Style':'http://www.southampton.ac.uk/~re1u06/software/redwards.css',  # Style file for HTML pages
                'EMail':'../_Documentation/email.png',                      # E-mail image file
                'Background':'',                                            # Website background image file
                'Authors':'Richard J. Edwards','Copyright':'2013',
                'Address':'Centre for Biological Sciences, University of Southampton, UK.',
                'Contact':'software@cabbagesofdoom.co.uk'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PyDoc Class                                                                                             #
#########################################################################################################################
class PyDoc(rje_obj.RJE_Object):     
    '''
    PyDoc python documentation & distribution Class. Author: Rich Edwards (2011).

    Str:str
    - Author = Author name to put at bottom of webpages [RJ Edwards]
    - DocDir = Output path for Python Module documentation ['../docs/']
    - DocSource = Input path for Python Module documentation (manuals etc.) ['../docs/']
    - DistDir = Output directory for distribution directories ['../packages/']
    - EMail = E-Mail address for general contact [seqsuite@gmail.com]
    - LogoURL = URL to SLiMSuite program logos ['http://www.slimsuite.unsw.edu.au/graphics/']
    - ManualURL = URL to SLiMSuite program manuals ['http://docs.slimsuite.unsw.edu.au/software/slimsuite/docs/manuals/']
    - Name = Name for PyDoc run. Used for file naming and within documentation files. ['pydocs']
    - PyPath = Path to python modules. Will also look in listed sourcedir subfolders. [../]
    - Release = Release for packages [YYMMDD]
    - RestURL = URL for REST server to pull outfmt docs ['http://rest.slimsuite.unsw.edu.au/']
    - StylePath = Path to style sheets and javascript code folder etc. [http://www.southampton.ac.uk/~re1u06/]
    - WebDir = Output directory for webpages [../html/']
    
    Bool:boolean
    - AddImports = Whether to add identified imported modules to python module list [True]
    - Calls = Whether to output Method Calls [False]
    - FullDoc = Whether to generate full docstring output including Classes and Methods [False]
    - HTML = Whether to make a basic HTML page of module docstrings (will make linked fun later) [False]
    - MakePages = Special run to generate default cmd/ and docs/ pages for commandline option docs [False]
    - Self = Whether to include 'self' calls of methods if calls=T [False]

    Int:integer
    - MethodCap = Maximum number of method calls before collapsed to single line [0]

    Num:float
    
    List:list
    - Distribute = Names of distributions - gets details from distributions.txt ['seqsuite','slimsuite']
    - Keywords = List of keywords to add in additon to module names/descriptions []
    - MethodSkip = List of methods to skip documentation ['makeInfo','cmdHelp','setupProgram','Main']
    - ModLinks = List of rje_ppi formats for module import links (xgmml,tdt,svg,r) (dev=T) [xgmml]
    - PyList = List of modules to upload. Can have * wildcards. Will add '.py' if missing. ['*']
    - SourceDir = List of subdirectories in which to look for modules [tools,extras,libraries]
    - StyleSheets = List of style sheets to use for HTML ['rje_tabber.css','re1u06.css']
    - Updates = List of Updates to be stored for output into documentation. 

    Dict:dictionary
    - PyModules = Dictionary of {path:Object} for python modules

    Obj:RJE_Objects
    - DB = Main Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Author','DocDir','DocSource','DistDir','Name','PyPath','RestURL','StylePath','WebDir','EMail','Release','LogoURL','ManualURL']
        self.boollist = ['AddImports','Calls','FullDoc','HTML','MakePages','Self']
        self.intlist = ['MethodCap']
        self.numlist = []
        self.listlist = ['Distribute','Keywords','MethodSkip','OtherFiles','PyList','SourceDir','StyleSheets','Updates','ModLinks']
        self.dictlist = ['PyModules']
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'DocDir':rje.makePath('../docs/'),'DocSource':rje.makePath('../docs/'),
                      'DistDir':rje.makePath('../packages/'),'Name':'pydocs',
                      'PyPath':rje.makePath('../'),'StylePath':rje.makePath('http://www.southampton.ac.uk/~re1u06/'),
                      'WebDir':rje.makePath('../html/'),'Author':'RJ Edwards','EMail':'seqsuite@gmail.com',
                      'LogoURL':'http://www.slimsuite.unsw.edu.au/graphics/',
                      'ManualURL':'http://docs.slimsuite.unsw.edu.au/software/slimsuite/docs/manuals/',
                      'RestURL':'http://rest.slimsuite.unsw.edu.au/'})
        self.setBool({'AddImports':True})
        self.setInt({})
        self.setNum({})
        ymd = time.localtime(time.time())[0:3]
        self.str['Release'] = '%s-%s-%s' % (ymd[0],rje.preZero(ymd[1],12),rje.preZero(ymd[2],31))
        self.list['SourceDir'] = ['tools','extras','libraries','legacy']
        self.list['PyList'] = ['*']
        self.list['MethodSkip'] = ['makeInfo','cmdHelp','setupProgram','Main']
        self.list['StyleSheets'] = ['rje_tabber.css','re1u06.css']
        self.list['Distribute'] = ['SLiMSuite']
        self.list['ModLinks'] = ['xgmml']
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
                self._cmdReadList(cmd,'str',['Author','Name','EMail','Release','LogoURL','ManualURL','RestURL'])
                self._cmdReadList(cmd,'path',['DocDir','DocSource','DistDir','PyPath','StylePath','WebDir'])
                self._cmdReadList(cmd,'bool',['AddImports','Calls','FullDoc','HTML','MakePages','Self'])
                self._cmdReadList(cmd,'int',['MethodCap'])
                self._cmdReadList(cmd,'list',['Distribute','Keywords','MethodSkip','PyList','SourceDir','StyleSheets'])  
                self._cmdReadList(cmd.lower(),'list',['ModLinks'])  
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            if self.dev(): print self.list['PyList'][0:]; return
            self.parseModules()
            if self.getBool('MakePages'): return self.makePages()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('FullDoc'):
                self.saveDocs('%s%s.txt' % (self.getStr('DocDir'),self.getStr('Name')))
                if self.getBool('HTML'):
                    self.htmlDocs('%s%s.html' % (self.getStr('DocDir'),self.getStr('Name')))
                    #if self.dev():
                    self.htmlReadMe(self.getStr('DocDir')+'readme/',modlist=[],zipfile='',level=0,onefile=False)
                    #else: self.htmlReadMe(self.getStr('DocDir'),modlist=[],zipfile='',level=0)
            ## ~ [2a] PPI Module import links output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['ModLinks']: self.ppiModLinks()
            ### ~ [3] ~ Generate distributions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for software in self.list['Distribute']: self.distribute(software)
            self.makePackageIndex()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                                #V2.0
        '''Main class setup method.'''
        try:### ~ [1] Expand PyList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pylist = []     # Generate a new list containing actual modules
            for subdir in [''] + self.list['SourceDir']:
                pypath = '%s%s' % (self.getStr('PyPath'),rje.makePath(subdir))
                if not os.path.exists(pypath): self.printLog('#WARN','Warning! Path "%s" not found!' % pypath); continue
                for pymod in self.list['PyList']:
                    if pymod[-3:] != '.py': pymod = '%s.py' % pymod
                    pylist += rje.getFileList(self,pypath,[pymod],subfolders=False,filecount=len(pylist))
            self.printLog('\r#PYMOD','%s Python Modules found' % rje.iLen(pylist))
            ## ~ [1a] ~ Sanity check that enough found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if len(pylist) < len(self.list['PyList']):
                self.printLog('#WARN','Warning! %d PyList elements given but only %d modules found!' % (len(self.list['PyList']),len(pylist)))
                if self.i() >= 0 and not rje.yesNo('Continue?'): return False
            ### ~ [2] Setup Database Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db.addEmptyTable('Module',['File','SourceDir','Module','Program','Description','Version','Last Edit','Author','Imports','Imported_By','Classes','Methods'],['File'])
            db.addEmptyTable('Class',['File','Class','Methods'],['File','Class'])
            db.addEmptyTable('Method',['File','Class','Method','CallAttributes'],['File','Class','Method'])
            db.addEmptyTable('Options',['File','Class','Type','Attribute','Argument','ArgType','Description','Default'],['File','Class','Argument'])
            ### ~ [3] Create directories if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('fulldoc'): rje.mkDir(self,self.getStr('DocDir'))
            if self.basefile().lower() in ['','none']: self.basefile('%s%s' % (self.getStr('DocDir'),self.getStr('Name')))
            ### ~ [4] Finish and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['PyList'] = pylist               
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def findModule(self,pymod,addlib=False):    ### Looks for file corresponding to module and returns if found     #V2.0
        '''Looks for file corresponding to module and returns if found.'''
        try:### ~ [1] Check each directory in turn and return file if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sourcedir = self.list['SourceDir'][0:]
            if addlib: sourcedir += ['tools','libraries']
            if pymod[-3:] != '.py': pymod = '%s.py' % pymod
            for subdir in [''] + sourcedir:
                pyfile = '%s%s%s' % (self.getStr('PyPath'),rje.makePath(subdir),pymod)
                if os.path.exists(pyfile): return pyfile
            return None
        except: self.errorLog('Problem during %s findModule(%s).' % (self,module)); return None
#########################################################################################################################
    ### <3> ### Reading and Parsing Methods                                                                             #
#########################################################################################################################
    def parseModules(self):  ### Reads in and parses modules from PyList                                            #V2.0
        '''Reads in and parses modules from PyList.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pyfile in self.list['PyList'][0:]: self.addModule(pyfile)
            for x in ['Module','Class','Method']:
                self.printLog('#%s' % x.upper(),'%s %s parsed' % (rje.iStr(self.db(x).entryNum()),x))
            ### ~ [2] Update Imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.updateImports()
        except: self.errorLog('%s.parseModules error' % self)
#########################################################################################################################
    def updateImports(self,modlist=[]): ### Updates the Imported_By lists for entries in modlist only
        '''
        Updates the Imported_By lists for entries in modlist only. If modlist is empty then all modules are updated.
        >> modlist:list = List of module names for updated import data
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not modlist: modlist = self.db('Module').indexKeys('Module')
            pyfiles = self.sortedPyFiles(modlist)
            for py in pyfiles: self.db('Module').data(py)['Imported_By'] = 'None'
            ### ~ [2] Update Imported_By ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pyfile in pyfiles:
                entry = self.db('Module').data(pyfile)
                for imod in string.split(entry['Imports'],', '):
                    if imod in self.db('Module').index('Module'):
                        for mentry in self.db('Module').indexEntries('Module',imod):
                            if mentry['Imported_By'] != 'None': mentry['Imported_By'] += ', %s' % entry['Module']
                            else: mentry['Imported_By'] = entry['Module']
        except: self.errorLog('%s.updateImports error' % self)
#########################################################################################################################
    def addModule(self,pyfile,docsonly=False):  ### Reads in and parses module from file                         #V2.15.0
        '''
        Reads in and parses module from file.
        >> pyfile:str = Python file to parse.
        >> docsonly:bool [False] = Whether parsing is limited for REST doc pages only.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pyfile = os.path.abspath(pyfile)
            if not docsonly and pyfile in self.db('Module').datakeys(): return   # Already parsed
            fulltext = self.loadFromFile(pyfile,v=0,checkpath=True,chomplines=False)
            data = {'File':pyfile,'Module':rje.baseFile(pyfile,strip_path=True),'FullText':fulltext,'DocString':'','Description':'',
                    'Version':'','Program':'','Last Edit':'','Author':'Unknown','SourceDir':'Error',
                    'Imports':[],'Imported_By':'None','Classes':[],'Methods':[]}
            classes = []; methods = []; imports = []
            modargs = {}    # Dictionary of argument:(type,default) parsed from docstring
            if not docsonly: odb = self.db('Options')
            ## ~ [1a] ~ Source directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for sourcedir in self.list['SourceDir']:
                module = data['Module']
                if pyfile.find(sourcedir) >= 0 and os.path.exists('%s%s%s.py' % (self.getStr('PyPath'),rje.makePath(sourcedir),module)): 
                    data['SourceDir'] = sourcedir
                    break
            ### ~ [2] Parse Full Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            content = fulltext[0:]      # This is will be eaten away and processed
            ptot = len(content); px = 0.0
            i = 0
            fulltext = ''
            while i < len(content):
                self.progLog('\r#MOD','Parsing module %s: %.1f%%' % (pyfile,px/ptot)); px += 100.0
                line = content[i]
                ## ~ [2a] ~ Parse DocString ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line[0:3] in ['\'\'\'','"""']:   # Started doc
                    if rje.matchExp('\'\'\'(.+)\'\'\'',line):
                        data['DocString'] = rje.matchExp('\'\'\'(.+)\'\'\'',line)[0] + '\n'
                        i += 1
                        continue
                    if rje.matchExp('"""(.+)"""',line):
                        data['DocString'] = rje.matchExp('"""(.+)"""',line)[0] + '\n'
                        i += 1
                        continue
                    docstring = ''
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if line[0:3] in ['\'\'\'','"""']:   # Ended docstring
                            break
                        else:
                            if docsonly: pass
                            elif rje.matchExp('\s+(\S+)=(\S+)\s*[:#]\s+(.+)\s*\[(.*)\]',line):
                                argdata = rje.matchExp('\s+(\S+)=(\S+)\s*[:#]\s+(.+)\s*\[(.*)\]',line)
                                modargs[argdata[0]] = (argdata[1],argdata[2],argdata[3])
                                odb.addEntry({'File':pyfile,'Class':'__doc__','Type':'cmd','Argument':argdata[0],'ArgType':argdata[1],'Default':argdata[3],'Description':argdata[2]})
                            elif rje.matchExp('\s+(\S+)=(\S+)\s*[:#].+\[(.*)\]',line):
                                self.warnLog('Cannot correctly parse command option: %s' % line,'cmdparse',quitchoice=self.test(),suppress=True,dev=True)
                            docstring += line
                            i += 1
                            for key in ['Description','Version','Last Edit','Program','Module']:
                                if rje.matchExp('^%s:\s+(\S.+)\s*$' % key,line):
                                    if key == 'Module':
                                        if not data['Program']: data['Program'] = rje.matchExp('^%s:\s+(\S.+)\s*$' % key,rje.chomp(line))[0]
                                    else: data[key] = rje.matchExp('^%s:\s+(\S.+)\s*$' % key,rje.chomp(line))[0]
                    data['DocString'] = docstring
                    i += 1
                ## ~ [2b] Parse Class Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('class') == 0 and docsonly: break
                elif line.find('class') == 0 and not data['DocString']:
                    self.printLog('#ERR','Cannot parse data from %s - format not recognised' % pyfile)
                    break
                    return False
                elif rje.matchExp('^class (\S+)\(',line):   # Newclass
                    line = string.replace(line,'(',' (',1)
                    classname = rje.matchExp('^class (\S+) \(',line)[0]
                    classtext = [line]
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if re.search('^def ',line) or re.search('^class ',line):    # Reached next thing
                            break
                        else:
                            classtext.append(line) 
                            i += 1
                    data['Classes'].append((classname,classtext))
                ## ~ [2c] Parse Method Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif rje.matchExp('^def (\S+)\((.*)\):',line):   # New method
                    line = string.replace(line,'(',' (',1)
                    method_detail = rje.matchExp('^def (\S+) \((.*)\):',line)
                    methodtext = [line]
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if re.search('^def ',line) or re.search('^class ',line) or re.search('^if ',line): break   # Reached next thing                            
                        else:
                            methodtext.append(line) 
                            i += 1
                        if rje.matchExp("author = \'(\S.+)\'",line) and data['Author'] == 'Unknown': data['Author'] = rje.matchExp("author = \'(\S.+)\'",line)[0]
                    if method_detail[0] not in self.list['MethodSkip']:
                        data['Methods'].append((method_detail[0],method_detail[1],methodtext))
                ## ~ [2d] Parse "Main" method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif re.search('^if __name__',line):   # "Main"
                    methodtext = [line]
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if re.search('^def ',line) or re.search('^class ',line):    # Reached next thing
                            break
                        else:
                            methodtext.append(line) 
                            i += 1
                    if 'Main' not in self.list['MethodSkip']:
                        data['Methods'].append(('Main','n/a',methodtext))
                ## ~ [2e] Other line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:                
                    fulltext += line
                    i += 1
            if not data['Program']: data['Program'] = data['Module']
            if data['Author'] == 'Unknown' and '_' in data['Module']: data['Author'] = string.split(data['Module'],'_')[0]
            ## ~ [2f] Extract imported modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for line in string.split(fulltext,'\n'):
                if rje.matchExp('^\s*import (\S.+)$',line) or rje.matchExp('^\s*from (\S+)',line):
                    try: implist = string.split(string.replace(rje.matchExp('^\s*import (\S.+)$',line)[0],', ',','),',')
                    except: implist = string.split(rje.matchExp('^\s*from (\S+)',line)[0])
                    for impmod in implist:
                        newmod = string.split(impmod)[0]
                        newfile = self.findModule(newmod,addlib=True)
                        if not newfile: continue
                        data['Imports'].append(newmod)
            importlist = data['Imports'][0:]
            ## ~ [2g] Update fulltext or return if docsonly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if docsonly: return data
            data['FullText'] = fulltext     # This now contains everything not extracted for other purposes!
            ptot = len(data['Classes']) + len(data['Methods']) + 1; px = 100.0
            self.progLog('\r#MOD','Parsing module %s: %.1f%%' % (pyfile,px/ptot))
            ### ~ [3] Parse Classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            classlist = []
            for (classname,classtext) in data['Classes']:
                self.progLog('\r#MOD','Parsing module %s: %.1f%%' % (pyfile,px/ptot)); px += 1
                entry = {'File':pyfile,'Class':classname,'DocString':'','Methods':[],'FullText':classtext}
                self.addClass(entry,modargs)
                classlist.append(classname)
            ### ~ [4] Parse Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for (methodname,methodarg,methodtext) in data['Methods']:
                self.progLog('\r#MOD','Parsing module %s: %.1f%%' % (pyfile,px/ptot)); px += 1
                entry = {'File':pyfile,'Class':'Module','Method':methodname,'DocString':'',
                         'CallAttributes':methodarg,'FullText':methodtext}
                self.addMethod(entry)
            ### ~ [5] Update main database entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data['Classes'] = string.join(classlist,'; ')
            data['Methods'] = len(data['Methods'])
            data['Imports'] = string.join(data['Imports'],', ')
            self.printLog('\r#MOD','Parsing module %s complete: %d Classes; %d Methods' % (pyfile,len(classlist),data['Methods']))
            self.db('Module').addEntry(data)
            #for h in self.db('Module').fields(): self.deBug('%s: %s' % (h,data[h]))
            ### ~ [6] Look for Imported Modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for newmod in importlist:
                newfile = self.findModule(newmod)
                if not newfile or newfile in self.list['PyList']: continue
                self.list['PyList'].append(newfile)
                self.addModule(newfile)
                self.printLog('#MOD','Added module %s to list (imported by %s)' % (newmod,data['Module']))
        except: self.errorLog('%s.parseModule(%s) error' % (self,pyfile))
#########################################################################################################################
    def addClass(self,entry,modargs={}):   ### Reads in and parses class from entry['FullText']                     #V2.7
        '''
        Reads in and parses class from entry['FullText'].
        >> entry:dict = Class data parsed from module
        >> modargs:dict = dictionary of argument:(type,default) parsed from module docstring.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            odb = self.db('Options')
            content = entry['FullText'][1:]
            fulltext = ''
            i = 0
            #self.deBug(modargs)
            ### ~ [2] Parse Content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while i < len(content):
                line = content[i]
                ## ~ [2a] Docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if re.search('^\s+\'\'\'',line) or re.search('^\s+"""',line):   # Started doc
                    if rje.matchExp('\'\'\'(.+)\'\'\'',line):
                        entry['DocString'] += rje.matchExp('\'\'\'(.+)\'\'\'',line)[0] + '\n'
                        i += 1
                        continue
                    if rje.matchExp('"""(.+)"""',line):
                        entry['DocString'] += rje.matchExp('"""(.+)"""',line)[0] + '\n'
                        i += 1
                        continue
                    docstring = ''
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if re.search('^\s+\'\'\'',line) or re.search('^\s+"""',line): break  # Ended docstring
                        else: docstring += line; i += 1
                    entry['DocString'] = docstring
                    i += 1
                ## ~ [2b] Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif re.search('\s+def (\S+)\((.*)\):',line):   # New method
                    line = string.replace(line,'(',' (',1)
                    method_detail = rje.matchExp('\s+def (\S+) \((.*)\):',line)
                    methodtext = [line]
                    i += 1
                    cmdread = False
                    while i < len(content):
                        line = content[i]
                        if re.search('^\s+def ',line): break
                        else: methodtext.append(line); i += 1
                        #self.deBug(line)
                        ## Special commandline option parsing
                        if '#' in line: line = line[:line.find('#')]
                        line = string.join(string.split(line),'')
                        line = string.replace(line,'.lower()','')
                        line = string.replace(line,'.upper()','')
                        if cmdread and not 'self' in line:
                            line = 'self._cmdReadList(cmd,%s,[%s' % (ctype,line)
                            if line[-1] == ')': line = line[:-1]    # Remove trailing ) for last line of argument
                        if rje.matchExp('self\._cmdRead(\S+)',line) and '%s' not in line:   # New Commandline argument(s)
                            try: cmdata = string.split(rje.matchExp('self\._cmd\S+\(([^)]+)',line)[0],',')
                            except:
                                self.warnLog('Troubling parsing cmdRead: %s' % line,'cmdreadparse',quitchoice=True,suppress=True,dev=True)
                                continue
                            for cx in range(len(cmdata)):
                                cmdata[cx] = string.replace(string.split(cmdata[cx],'=')[-1],'"','')
                                cmdata[cx] = string.replace(cmdata[cx],"'","")
                                cmdata[cx] = string.replace(cmdata[cx]," ","")
                            while '' in cmdata: cmdata.remove('')
                            try: cmd = cmdata[0]; ctype = cmdata[1]
                            except: self.deBug('%s: %s' % (entry['File'],cmdata))
                            if cmdata[-1] == 'cmd':
                                cmdata = [cmdata[-1]] + cmdata[:-1]
                            #self.deBug(cmdata)
                            cmdargs = [cmdata[-1]]
                            cmdatt = None
                            if cmdata[2][0] == '[':    # self._cmdReadList(cmd,type,optionlist)
                                cmdargs = cmdata[2:]
                                cmdargs[0] = cmdargs[0][1:]
                                if cmdargs[-1][-1] == ']': cmdargs[-1] = cmdargs[-1][:-1]; cmdread = False
                                else: cmdRead = True
                            elif len(cmdata) > 3: cmdatt = cmdata[2]
                            if len(cmdargs) > 1 and cmdatt: self.warnLog('Conflicting argument/attribute numbers from %s: %s' % (entry['File'],line),quitchoice=True,suppress=True,dev=True)
                            for arg in cmdargs:
                                #self.deBug(line)
                                argentry = {'File':entry['File'],'Class':entry['Class'],'Type':ctype,'Argument':arg.lower(),'ArgType':'?','Default':'?'}
                                if cmdatt: argentry['Attribute'] = cmdatt
                                else: argentry['Attribute'] = arg
                                if arg.lower() in modargs: [argentry['ArgType'],argentry['Description'],argentry['Default']] = modargs[arg.lower()]
                                #self.deBug(argentry)
                                odb.addEntry(argentry)
                                docentry = odb.makeKey({'File':entry['File'],'Class':'__doc__','Argument':arg.lower()})
                                if docentry in odb.data(): odb.dict['Data'].pop(docentry)
                            #!# Can we also read defaults? How?! setStr etc. ?? Too complex, I think. docstring?!
                            #!# Read docstring, parse options and defaults and match where possible.
                            #!# ArgType and Default parsed from docstring, matching on 'Argument'
                    ## Add Method
                    entry['Methods'].append((method_detail[0],method_detail[1],methodtext))
                ## ~ [2c] Remaining text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else: fulltext += line; i += 1
            entry['FullText'] = fulltext
            ### ~ [3] Add Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for (methodname,methodarg,methodtext) in entry['Methods']:
                data = {'File':entry['File'],'Class':entry['Class'],'Method':methodname,'DocString':'',
                         'CallAttributes':methodarg,'FullText':methodtext}
                self.addMethod(data)
            ### ~ [4] Finish and Update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            entry['Methods'] = len(entry['Methods'])
            self.db('Class').addEntry(entry)
        except: self.errorLog('%s.parseClass() error' % self)
#########################################################################################################################
    def addMethod(self,entry):   ### Reads in and parses class from entry['FullText']                                #V2.0
        '''Reads in and parses class from entry['FullText'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            i = 0
            content = entry['FullText'][1:]
            fulltext = ''
            self.verbose(2,2,'%s(%s)...' % (entry['Method'],entry['CallAttributes']),0) #v,i
            ### ~ [2] Parse Content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while i < len(content):
                line = content[i]
                ## ~ [2a] Docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if re.search('^\s+\'\'\'',line) or re.search('^\s+"""',line):   # Started doc
                    if rje.matchExp('\'\'\'(.+)\'\'\'',line):
                        entry['DocString'] += rje.matchExp('\'\'\'(.+)\'\'\'',line)[0] + '\n'
                        i += 1
                        continue
                    if rje.matchExp('"""(.+)"""',line):
                        entry['DocString'] += rje.matchExp('"""(.+)"""',line)[0] + '\n'
                        i += 1
                        continue
                    docstring = ''
                    i += 1
                    while i < len(content):
                        line = content[i]
                        if re.search('^\s+\'\'\'',line) or re.search('^\s+"""',line):   # Ended docstring
                            i += 1
                            break
                        else:
                            docstring += line
                            i += 1
                    entry['DocString'] = docstring
                ## ~ [2b] Method Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:
                    fulltext += line
                    i += 1
            ### ~ [3] Finish and Update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.db('Method').addEntry(entry)
        except: self.errorLog('%s.parseMethod() error' % self)
#########################################################################################################################
    def parseDocFunction(self,docstring,html=False):    ### Parses the "function" part of the docstring 
        '''
        Parses the "function" part of the docstring and returns as text or html.
        >> docstring:str = the string from which to parse to the relevant text.
        >> html:boolean [False] = whether to return as HTML
        '''
        try:### ~ [1] Parse out function text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lines = string.split(self.info['DocString'],'\n')
            parse = False
            function = []
            for line in lines:
                if line.find('Function:') == 0: parse = True                    # Found Function section
                elif line[:1] not in ['',' ']: parse = False                    # Reached next section
                elif parse: function.append(string.join(string.split(line)))    # Add re-formed line to function
            ### ~ [2] Tidy function text and add HTML code if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for f in range(len(function)):
                if function[f][:1] == '*':
                    if html:
                        function[f] = '<LI>' + function[f][2:]
                        if function[f-1] == '': function[f-1] = '<UL>'
                        if f > len(function): function.append('</UL>')
                        elif function[f+1] == '': function[f+1] = '</UL>'
                elif html and f == 0: function[0] = '<P>' + function[0]
                elif not function[f]:
                    if html:
                        if f > len(function): function[f] = '</P>'
                        else: function[f] = '</P><P>'
                    else: function[f] = '\n\n'
            return string.join(function)
        except: self.errorLog('%s.parseDocFunction error' % self)
#########################################################################################################################
    ### <4> ### Documentation output Methods                                                                            #
#########################################################################################################################
    def sortedPyFiles(self,modlist=[]):    ### Returns Python Files in order to be dealt with
        '''
        Returns Python Files in order to be dealt with.
        >> modlist:list [] = List of modules to return files for. (If empty, will return all.)
        '''
        try:### ~ [1] Order by Source Directory and then Name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pylist = []
            db = self.db('Module')
            for sourcedir in self.list['SourceDir']:
                for pyfile in db.dataKeys():
                    entry = db.data(pyfile)
                    module = entry['Module']
                    if modlist and module not in modlist: continue
                    if entry['SourceDir'] == sourcedir: pylist.append(pyfile)
                    #if pyfile.find(sourcedir) >= 0 and os.path.exists('%s%s%s.py' % (self.getStr('PyPath'),rje.makePath(sourcedir),module)): 
                    #    pylist.append(pyfile)
            return pylist
        except:
            self.errorLog('Error in %s.sortedPyFiles()' % self); raise
#########################################################################################################################
    def saveDocs(self,filename='pydocs.txt',append=False):      ### Prints docs for modules, classes and methods to file
        '''
        Prints docs for modules, classes and methods to file.
        >> filename:str = output file name
        >> append:boolean
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###           
            if append:
                self.printLog('#DOC','Appending docstrings to %s' % filename)
                PYDOC = open(filename,'a')
            else:
                rje.mkDir(self,filename)
                self.printLog('#DOC','Writing docstrings to %s' % filename)
                PYDOC = open(filename,'w')
            db = self.db('Module'); cdb = self.db('Class'); mdb = self.db('Method'); odb = self.db('Options')
            db.saveToFile(); cdb.saveToFile(); mdb.saveToFile(); odb.saveToFile()
            dx = cx = mx = 0
            ### ~ [2] Output Docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sourcedir in self.list['SourceDir']:
                PYDOC.write('-%s:\n\n' % sourcedir)
                for pyfile in db.dataKeys():
                    entry = db.data(pyfile)
                    module = entry['Module']
                    if not pyfile.find(sourcedir) >= 0 or not os.path.exists('%s%s%s.py' % (self.getStr('PyPath'),rje.makePath(sourcedir),module)): continue
                    ## ~ [2a] ~ Module docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    mtxt = '### ~~~ Module %s ~ [%s] ~~~ ###' % (module,pyfile)
                    while len(mtxt) < 122: mtxt = mtxt[:5] + '~' + mtxt[5:-5] + '~' + mtxt[-5:]
                    PYDOC.write('%s\n\n%s\n' % (mtxt,entry['DocString'])); dx += 1
                    ## ~ [2b] ~ Class docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for centry in cdb.indexEntries('File',pyfile):
                        PYDOC.write('\n** Class %s **\n%s\n' % (centry['Class'],centry['DocString'])); cx += 1
                        for mentry in mdb.indexEntries('File',pyfile):
                            if mentry['Class'] != centry['Class']: continue
                            PYDOC.write('%s.%s(%s)\n%s\n' % (mentry['Class'],mentry['Method'],mentry['CallAttributes'],mentry['DocString'])); mx += 1
                    ## ~ [2c] ~ Method docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    PYDOC.write('\n** Module Methods **\n\n')
                    for mentry in mdb.indexEntries('File',pyfile):
                        if mentry['Class'] != 'Module': continue
                        PYDOC.write('%s.%s(%s)\n%s\n' % (module,mentry['Method'],mentry['CallAttributes'],mentry['DocString'])); mx += 1
                PYDOC.write('\n\n\n')
            PYDOC.close()            
            self.printLog('#DOC','Output to %s complete: %s modules; %s classes; %s methods' % (filename,rje.iStr(dx),rje.iStr(cx),rje.iStr(mx)))
        except: self.errorLog('Error in %s.saveDocs()' % self)     
#########################################################################################################################
    def htmlDocs(self,filename='pydocs.html'):      ### Prints docs for modules, classes and methods to file
        '''
        Prints docs for modules to HTML file.
        >> filename:str = output file name
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#HTML','HTML documentation output...')
            PYDOC = open(filename,'w')
            db = self.db('Module')
            ## ~ [1a] HTML Header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            title = 'RJE_PYDOCS %s documentation' % self.getStr('Name')   #!# Get better title
            stylesheets = []
            for css in self.list['StyleSheets']: stylesheets.append(self.getStr('StylePath')+css)
            htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=rje.sortKeys(db.index('Module')),javascript='%sjavascript/' % self.getStr('StylePath'))
            PYDOC.write(htmlhead)
            ### ~ [2] Output Module List, linking to relevant code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            PYDOC.write('<A NAME="Top"><H2>%s</H2></A>\n' % title)
            ### ~ [3] General Information Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headtabs = []
            ## ~ [3a] Module Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for sourcedir in self.list['SourceDir']:
                htm = ['<h3>Python Modules in %s:</h3>\n' % (rje.makePath(sourcedir))]
                htm.append('<UL>')
                for pyfile in db.indexDataList('SourceDir',sourcedir,'File'):             #dataKeys():
                    entry = self.db('Module').data(pyfile)
                    module = entry['Module']
                    program = entry['Program']
                    if not pyfile.find('%s/%s.py' % (sourcedir,module)) >= 0: continue
                    if entry['Version']: htm.append('<LI><B><A HREF="#%s-%s">%s [version %s]</A></B> %s' % (sourcedir,module,program,entry['Version'],entry['Description']))
                    else: htm.append('<LI><B><A HREF="#%s-%s">%s</A></B> %s' % (sourcedir,module,program,entry['Description']))
                    if program.lower() != module.lower(): htm[-1] += ' (%s.py)' % module
                htm.append('</UL>')
                headtabs.append((sourcedir,string.join(htm,'\n'),'Full list of python modules from %s/' % sourcedir))
            ## ~ [3b] Other Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['OtherFiles']:
                htm = ['<H3>Other Files</H3>','<UL>']
                for file in self.list['OtherFiles']: htm.append('<LI>%s' % os.path.basename(file))
                htm.append('</UL>')
                headtabs.append(('Other files',string.join(htm,'\n'),'Other files used by python modules'))
            ## ~ [3c] Updates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Updates']:
                htm = ['<H3>Updates since last release:</H3>'] + self.list['Updates']
                headtabs.append(('Updates',string.join(htm,'\n'),'Updates since last release'))
            ## ~ [3d] GNU License ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htm = ['<HR>','<H2>GNU License</H2>',
                   '<P>Copyright (C) %s %s &lt;%s&gt;' % (web_defaults['Copyright'],self.getStr('Author'),self.getStr('EMail')),'</P><P>',
                   'This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.',
                   '</P><P>','This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.',
                   '</P><P>','You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA'
                   '</P><P>','Author contact: &lt;%s&gt; / %s' % (self.getStr('EMail'),web_defaults['Address']),
                   '</P><P>','To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py']
            headtabs.append(('GNU',string.join(htm,'\n'),'GNU License agreement for use'))
            ## ~ [3x] Output Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            PYDOC.write(rje_html.tabberHTML('readme',headtabs,level=0))
            ### ~ [3] Output Module Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            PYDOC.write(self.pydocHTML(db.indexKeys('Module')))
            ### ~ [4] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            PYDOC.write(rje_html.htmlTail('%s %s' % (self.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1])))
            PYDOC.close()
            self.printLog('\r#HTML','HTML documentation output to %s' % filename)
        except: self.errorLog('Error in %s.htmlDocs()' % self)     
#########################################################################################################################
    def pyModLinkRef(self,pyfile,module,onefile=True):  ### Returns the <A NAME> cross-reference for module
        '''
        Returns the <A NAME> cross-reference for module.
        >> pyfile:str = File name of the current python module (where to start looking)
        >> module:str = Name of module to link to
        >> onefile:bool [True] = whether to output full readme as one file (True) or one file per module (False)
        '''
        aname = ''
        for sourcedir in self.list['SourceDir']:
            for otherfile in self.db('Module').dataKeys():
                if not otherfile.find('%s/%s.py' % (sourcedir,module)) >= 0: continue
                if pyfile.find('%s/' % sourcedir) >= 0 or not aname:
                    if onefile: aname = "#%s-%s" % (sourcedir,module)
                    else: aname = "../%s/%s.html" % (sourcedir,module)
                #self.deBug('%s:%s > %s' % (pyfile,module,aname))
        return aname
#########################################################################################################################
    def pydocHTML(self,modlist,level=0,onefile=True,outdir=''):    ### Returns HTML code for listed modules
        '''
        Returns HTML code for listed modules.
        >> modlist:list = List of modules needing to be output
        >> level:int [0] = Tabber tab level
        >> onefile:bool [True] = whether to output full readme as one file (True) or one file per module (False)
        >> outdir:string [''] = output directory for readme files
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = ''
            pyfiles = self.sortedPyFiles(modlist)
            ### ~ [2] Output Module Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sourcedir in self.list['SourceDir']:
                hdir = '%s%s/' % (outdir,sourcedir)
                if not onefile: rje.mkDir(self,hdir)
                html += '\n\n<!-- ~~~ Modules in %s ~~~ -->\n<h3>%s/:</h3>\n' % (sourcedir,sourcedir)
                for pyfile in pyfiles:
                    entry = self.db('Module').data(pyfile)
                    module = entry['Module']
                    if not pyfile.find('%s/%s.py' % (sourcedir,module)) >= 0: continue
                    html += '\n\n<A NAME="%s-%s"></A>' % (sourcedir,module)
                    if onefile: html += self.modHTML(pyfile)
                    else:
                        hfile = '%s%s.html' % (hdir,module)
                        HTML = open(hfile,'w')
                        title = '%s/%s ReadMe' % (sourcedir,module)
                        keywords = self.list['Keywords'] + [module]
                        stylepath = self.getStr('StylePath')
                        if stylepath.startswith('.'): stylepath = '../%s' % self.getStr('StylePath')
                        stylesheets = []
                        for css in self.list['StyleSheets']: stylesheets.append(stylepath+css)
                        htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=keywords,javascript=stylepath)
                        HTML.write(htmlhead)
                        #HTML.write('<p><a href="%sreadme.html">Return to main %s readme.html</a></p>' % (outdir,self.getStr('Name')))
                        HTML.write(self.modHTML(pyfile,onefile=False))
                        HTML.write(rje_html.htmlTail('%s %s' % (self.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1])))
                        HTML.close()
                        html = '<p>See individual module readme pages for details.</p>'
            return html
        except: self.errorLog('Error in %s.pydocHTML()' % self,quitchoice=True)
#########################################################################################################################
    def modHTML(self,pyfile,level=0,modlist=[],reduced=False,onefile=True):     ### Returns HTML code for pymodule documentation
        '''
        Returns HTML code for pymodule documentation.
        >> pyfile:str = The python file identifying the module for output
        >> level:int [0] = The tabber tab level (increment if nested inside other tabs)
        >> modlist:list [] = List of python modules on page for linking purposes
        >> reduced:bool [False] = Whether to return a reduced docstring with limited data
        >> onefile:bool [True] = whether to output full readme as one file (True) or one file per module (False)
        '''
        ### ~ [1] Output header and reference point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cdb = self.db('Class'); mdb = self.db('Method')
        try:
            entry = self.db('Module').data(pyfile)
            if entry['Version']: pyhtml = '<h2>%s [version %s] %s ~ ' % (entry['Program'],entry['Version'],entry['Description'])
            else: pyhtml = '<h2>%s %s ~ ' % (entry['Program'],entry['Description'])
            if entry['Program'].lower() != entry['Module'].lower(): pyhtml += ' (%s.py)' % entry['Module']
        except:
            self.errorLog('%s cannot be processed with modHTML' % pyfile)
            return ''
        if onefile: pyhtml += '<font size=-1>[<A HREF="#Top">Top</A>]</font></h2>\n\n'
        else: pyhtml += '<font size=-1>[<A HREF="../readme.html">%s Readme</A>]</font></h2>\n\n' % self.getStr('Name')
        #PYDOC.write('source: %s</p><p>' % (pyfile))
        ### ~ [2] Output and format docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        docstring = string.replace(entry['DocString'],'<','&gt;')
        docstring = string.replace(docstring,'>','&lt;')
        docstring = string.split(docstring,'\n')
        htmldoc = ['<p>']
        code = False    # Whether currently in a code block
        for line in docstring:
            #if re.search('^(\S.+:)(.+)$',line):
            #    match = rje.matchExp('^(\S.+:)(.*)$',line)
            if code and not (re.search('^(.+\s)(\S+)=(\S+)(\s*:\s*)(\S.+)$',line) or rje.matchExp('^\s+(-)',line)):
                htmldoc[-1] += '</code></pre>\n'; code = False
            if ':' in line:
                match = (line[:line.find(':')+1],line[line.find(':')+1:])
                #self.deBug('%s = %s' % (entry['Module'],match))
                if string.split(match[0])[-1] in ['Function:','Functions:','Commands:','Commandline:','Options:']:
                    htmldoc[-1] += '</p>\n\n'
                    if reduced and string.split(match[0])[-1] in ['Function:','Functions:']: htmldoc[-1] += '<p>'
                    else: htmldoc.append('<p>')
                    htmldoc[-1] += '<B>%s</B>%s<BR>' % match
                elif string.split(match[0])[-1] in ['Cite:','Citation:']:
                    if rje.matchExp('^%s\s+(\S.+)\s+\[PMID:\s*(\d+)\]' % match[0],line):
                        pmid = rje.matchExp('^%s\s+(\S.+)\s+\[PMID:\s*(\d+)\]' % match[0],line)
                        match = (match[0],' <a href="http://www.ncbi.nlm.nih.gov/pubmed/%s?dopt=Abstract" TARGET="_blank">%s</a>' % (pmid[1],pmid[0]))
                    #else: match = (string.split(match[0])[0],'%s%s' % (string.join(string.split(match[0])[1:],match[1])))
                    htmldoc[-1] += '<B>%s</B>%s<BR>' % match
                elif match[0] in ['Module:','Program:','Description:','Version:','Last Edit:']:
                    htmldoc[-1] += '<B>%s</B>%s<BR>' % match
                elif rje.matchExp('^(.+\s)(\S+)=(\S+)(\s*:\s*)(\S.+)$',line):
                    match = rje.matchExp('^(.+\s)(\S+)=(\S+)(\s*:\s*)(\S.+)$',line)
                    if code: htmldoc[-1] += '%s<B><FONT COLOR=#AB1210>%s=</FONT><FONT COLOR=#653A28>%s</FONT>%s</B>%s\n' % match
                    else: htmldoc[-1] += '</p>\n<pre><code>%s<B><FONT COLOR=#AB1210>%s=</FONT><FONT COLOR=#653A28>%s</FONT>%s</B>%s\n' % match; code = True
                elif rje.matchExp('^\s*(\S+)=\S+',match[0]):
                    htmldoc[-1] += '<B>%s</B>%s<BR>' % match
                else: htmldoc[-1] += '%s<BR>' % line
                if match[0] == 'Last Edit:':    # Add Imports and Imported By
                    if reduced: continue
                    htmldoc[-1] += '<B>Imports:</B> '
                    for imod in string.split(entry['Imports'],', '):
                        aname = self.pyModLinkRef(pyfile,imod,onefile)
                        if aname: htmldoc[-1] += '<A HREF="%s">%s</A>, ' % (aname,imod)
                        else: htmldoc[-1] += '%s, ' % (imod)
                    if htmldoc[-1][-2:] == ', ': htmldoc[-1] = htmldoc[-1][:-2]
                    htmldoc[-1] += '<BR>'
                    htmldoc[-1] += '<B>Imported By:</B> '
                    for imod in string.split(entry['Imported_By'],', '):
                        aname = self.pyModLinkRef(pyfile,imod,onefile)
                        if aname: htmldoc[-1] += '<A HREF="%s">%s</A>, ' % (aname,imod)
                        else: htmldoc[-1] += '%s, ' % (imod)
                    if htmldoc[-1][-2:] == ', ': htmldoc[-1] = htmldoc[-1][:-2]
                    htmldoc[-1] += '<BR>'
            elif re.search('^\S',line) or re.search('^\s*#',line):
                htmldoc[-1] += '<B>%s</B><BR>' % string.replace(rje.chomp(line),'(C)','&copy;')
            elif re.search('^(.+\s)(\S+)=(\S+)(\s*:\s*)(\S.+)$',line):
                match = rje.matchExp('^(.+\s)(\S+)=(\S+)(\s*:\s*)(\S.+)$',line)
                #htmldoc[-1] += '%s<LI><B><FONT COLOR=RED>%s=</FONT><FONT COLOR=DARKRED>%s</FONT>%s</B>%s<BR>' % match
                #htmldoc[-1] += '%s<B><FONT COLOR=#AB1210>* %s=</FONT><FONT COLOR=#653A28>%s</FONT>%s</B>%s<BR>' % match
                if code: htmldoc[-1] += '%s<B><FONT COLOR=#AB1210>%s=</FONT><FONT COLOR=#653A28>%s</FONT>%s</B>%s\n' % match
                else:
                    htmldoc[-1] += '</p>'
                    htmldoc[-1] += '<pre><code>%s<B><FONT COLOR=#AB1210>%s=</FONT><FONT COLOR=#653A28>%s</FONT>%s</B>%s\n' % match; code = True
            elif re.search('\S',line):
                if code: htmldoc[-1] += '%s\n' % line
                else: htmldoc[-1] += '%s<BR>' % line
            else: htmldoc[-1] += '<BR>'
        if not docstring or htmldoc[-1] == '<p><BR>': htmldoc[-1] = '<p><i><FONT COLOR=#AB1210>This module was not recognised by rje_pydocs.</FONT></i>'
        elif not docstring[0]: htmldoc[-1] = '<p><i><FONT COLOR=#AB1210>This module was not recognised by rje_pydocs.</FONT></i>'
        elif docstring[0][:4] == '####': htmldoc[-1] = '<p><i><FONT COLOR=#AB1210>This module was not recognised by rje_pydocs.</FONT></i>'
        htmldoc[-1] += '</p>\n\n'
        doctab = [('Summary',htmldoc[0],'Summary of %s module' % entry['Module'])]
        for i in range(1,len(htmldoc)):
            tabid = string.split(htmldoc[i],':')[0][6:]
            if string.split(tabid)[0] == 'Commandline': tabid = 'Commandline'
            doctab.append((tabid,htmldoc[i],'Docstring %s for %s module' % (tabid,entry['Module'])))
        if reduced: return doctab
        ### ~ [3] ~ History and ToDo docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        todox = 0
        for mentry in mdb.indexEntries('File',pyfile):
            if mentry['Class'] != 'Module' or mentry['Method'] not in ['history','todo']: continue
            if mentry['Method'] == 'history':
                htm = '<h3>%s Module Version History</h3>\n\n' % (entry['Module'])
                htm += '<pre><code>%s</code></pre>\n\n' % (mentry['DocString'])
                doctab.append(('History',htm,'%s Module version history' % entry['Module']))
            else:
                htm = '<h3>%s Module ToDo Wishlist</h3>\n\n' % (entry['Module'])
                htm += '<pre><code>%s</code></pre>\n\n' % (mentry['DocString'])
                todox = len(doctab)
                doctab.append(('ToDo',htm,'%s Module ToDo Wishlist' % entry['Module']))
        ### ~ [4] Output Class and Method Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [4a] Class Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        for centry in cdb.indexEntries('File',pyfile):
            htm = '<h3>%s Class</h3>\n<pre>%s</pre>\n\n' % (centry['Class'],centry['DocString'])
            for mentry in mdb.indexEntries('File',pyfile):
                if mentry['Class'] != centry['Class']: continue
                docstring = string.replace(string.replace(mentry['DocString'],'>','&gt;'),'<','&lt;')
                if docstring and docstring[:4] != '    ':
                    htm += '<b>%s.%s(<i>%s</i>)</b>\n<pre><code>    %s</code></pre>\n\n' % (mentry['Class'],mentry['Method'],mentry['CallAttributes'],docstring)
                else: htm += '<b>%s.%s(<i>%s</i>)</b>\n<pre><code>%s</code></pre>\n\n' % (mentry['Class'],mentry['Method'],mentry['CallAttributes'],string.replace(docstring,'        ','    '))
            doctab.append(('%s Class' % centry['Class'],htm,'Documentation for %s Class and associated methods' % centry['Class']))
        ## ~ [4b] ~ Method docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        htm = '<h3>%s Module Methods</h3>\n\n' % entry['Module']; mx = 0
        for mentry in mdb.indexEntries('File',pyfile):
            if mentry['Class'] != 'Module' or mentry['Method'] in ['history','todo']: continue
            docstring = string.replace(string.replace(mentry['DocString'],'>','&gt;'),'<','&lt;'); mx += 1
            if docstring and docstring[:4] != '    ':
                htm += '<b>%s.%s(<i>%s</i>)</b>\n<pre><code>    %s</code></pre>\n\n' % (entry['Module'],mentry['Method'],mentry['CallAttributes'],docstring)
            else: htm += '<b>%s.%s(<i>%s</i>)</b>\n<pre><code>%s</code></pre>\n\n' % (entry['Module'],mentry['Method'],mentry['CallAttributes'],docstring)
        if mx: doctab.append(('Methods',htm,'Documentation for %s module-level methods' % entry['Module']))
        ### ~ [5] Write to PYDOC html file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if todox and not reduced: doctab.append(doctab.pop(todox))
        pyhtml += rje_html.tabberHTML(pyfile,doctab,level)
        return pyhtml
#########################################################################################################################
    def htmlReadMe(self,outdir='',modlist=[],zipfile='',level=0,linkothers=False,onefile=True):    ### Makes readme documentation from Modules
        '''
        Makes readme documentation for set of modules read in.
        >> outdir:string [''] = output directory for readme files
        >> modlist:list [] = list of modules to make readme for. If empty, use all read in.
        >> zipfile:string [''] = zipfile name. If none, will not provide installation instructions.
        >> level:int [0] = Tabber level for code
        >> linkothers:bool [False] = whether to add links to other files
        >> onefile:bool [True] = whether to output full readme as one file (True) or one file per module (False)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db('Module'); htmlcode = ''
            rje.mkDir(self,outdir)
            if not level: HTML = open(outdir+'readme.html','w')
            if not modlist: modlist = db.indexKeys('Module')
            pyfiles = self.sortedPyFiles(modlist)
            keywords = self.list['Keywords'] + modlist
            for entry in db.entryList(pyfiles): keywords.append(entry['Description'])
            stylesheets = []
            for css in self.list['StyleSheets']: stylesheets.append(self.getStr('StylePath')+css)
            ### ~ [1] HTML Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            title = '%s Release ReadMe' % self.getStr('Name')
            htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=keywords,javascript=self.getStr('StylePath'))   #%sjavascript/' % self.getStr('StylePath'))
            if not level: HTML.write(htmlhead)
            ### ~ [2] ReadMe header information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headerhtml = ['<A NAME="Top"><H2>ReadMe documentation for release of %s software</H2></A>' % self.getStr('Name'),
                          '<P><B>Distribution compiled:</B> %s' % time.asctime(time.localtime(time.time())),
                          '<P><B>Questions/Comments?:</B> please contact ',
                          '<A HREF="mailto:%s">%s</A></P>' % (self.getStr('EMail'),self.getStr('EMail'))]
            if level: htmlcode += string.join(headerhtml,'\n')
            else: HTML.write(string.join(headerhtml,'\n'))
            ### ~ [3] General Information Tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headtabs = []
            ## ~ [3a] Module Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for sourcedir in self.list['SourceDir']:
                htm = ['<h3>Python Modules in %s:</h3>\n' % (rje.makePath(sourcedir))]
                htm.append('<UL>')
                for pyfile in db.indexDataList('SourceDir',sourcedir,'File'):   #for pyfile in db.dataKeys():
                    entry = self.db('Module').data(pyfile)
                    module = entry['Module']
                    if module not in modlist: continue
                    #self.deBug('%s - %s - %s' % (sourcedir,module,pyfile))
                    program = entry['Program']
                    if onefile: href = '#%s-%s' % (sourcedir,module)    # Use placement within file
                    elif level: href = './readme/%s/%s.html' % (sourcedir,module)         # Use dir structure
                    else: href = './%s/%s.html' % (sourcedir,module)         # Use dir structure
                    if entry['Version']: htm.append('<LI><B><A HREF="%s">%s [version %s]</A></B> %s' % (href,program,entry['Version'],entry['Description']))
                    else: htm.append('<LI><B><A HREF="%s">%s</A></B> %s' % (href,program,entry['Description']))
                    if program.lower() != module.lower(): htm[-1] += ' (%s.py)' % module
                htm.append('</UL>')
                headtabs.append((sourcedir,string.join(htm,'\n'),'Full list of python modules in %s/' % sourcedir))
            ## ~ [3b] Other Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['OtherFiles']:
                docs = ['<H3>Manuals</H3>','<UL>']
                htm = ['<H3>Other Files</H3>','<UL>']
                for file in self.list['OtherFiles']:
                    if linkothers: htm.append('<LI><a href="%s">%s</a>' % (file,file))
                    else: htm.append('<LI>%s' % os.path.basename(file))
                    if file[-4:].lower() == '.pdf': docs.append(htm.pop(-1))
                htm.append('</UL>'); docs.append('</UL>')
                if len(docs) > 3: headtabs.append(('docs',string.join(docs,'\n'),'Documentation (PDF Manuals)'))
                if len(htm) > 3: headtabs.append(('files',string.join(htm,'\n'),'Other files used by python modules'))
            ## ~ [3c] Installation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if zipfile:
                htm = ['<HR>','<H2>Installation Instructions</H2>','<OL>',
                       '<LI>Place the %s.%s.tar.gz file in chosen directory (e.g. c:\\bioware\\) and unpack.' % (zipfile,self.getStr('Release')),
                       '<LI>A subdirectory %s will be created containing all the files necessary to run.' % zipfile,
                       '</OL>','<P>The software should run on any system that has <A HREF="http://www.python.org">Python</A> installed.',
                       'Additional software may be necessary for full functionality. Further details can be found in the manuals supplied.</P>']
                headtabs.append(('installation',string.join(htm,'\n'),'Installation instructions for %s' % self.getStr('Name')))
            ## ~ [3d] Updates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Updates']:
                htm = ['<H3>Updates since last release:</H3>'] + self.list['Updates']
                headtabs.append(('Updates',string.join(htm,'\n'),'Updates since last release'))
            ## ~ [3e] GNU License ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htm = ['<HR>','<H2>GNU License</H2>',
                   '<P>Copyright (C) %s %s &lt;%s&gt;' % (web_defaults['Copyright'],self.getStr('Author'),self.getStr('EMail')),'</P><P>',
                   'This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.',
                   '</P><P>','This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.',
                   '</P><P>','You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA'
                   '</P><P>','Author contact: &lt;%s&gt; / %s' % (self.getStr('EMail'),web_defaults['Address']),
                   '</P><P>','To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py']
            headtabs.append(('GNU',string.join(htm,'\n'),'GNU License agreement for use'))
            #!# Add citation and webserver details #!#
            ## ~ [3x] Output Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if level:
                htmlcode += rje_html.tabberHTML('readme',headtabs,level)
                if onefile: htmlcode += self.pydocHTML(modlist,level)
            else:
                HTML.write(rje_html.tabberHTML('readme',headtabs,level))
                HTML.write(self.pydocHTML(modlist,level,onefile,outdir))
            ### ~ [4] End of ReadMe ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if level: return htmlcode
            HTML.write(rje_html.htmlTail('%s %s' % (self.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1])))
            HTML.close()        
            self.printLog('#HTML','Generation of readme.html for %s distribution complete.' % self.getStr('Name'))
        except: self.errorLog('Error in PyDoc.makeReadMe(%s)' % self.getStr('Name'))     
#########################################################################################################################
    def ppiModLinks(self): ### Outputs module import relationships using rje_ppi.
        '''Outputs module import relationships using rje_ppi.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppi = rje_ppi.PPI(self.log,self.cmd_list+['PPISym=F'])
            ppi.basefile(self.basefile())
            db = ppi.obj['DB'] = self.db()
            ## ~ [0a] Setup PPI tables from Module table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ndb = db.copyTable(self.db('Module'),'Node')
            for mod in ndb.index('Module'):
                if len(ndb.index('Module')[mod]) > 1:
                    self.warnLog('Module "%s" duplicated' % entry['Module'],warntype='ModLink.ModDup',quitchoice=True,suppress=True,dev=True)
            ndb.renameField('Module','Node')
            ndb.makeField('#SourceDir#','SAMPLE')
            ndb.newKey('Node',True)
            ## ~ [0b] Generate Edge table from Module table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb = db.addEmptyTable('Edge',['Hub','Spoke','Evidence','ParentType','Type','Links','Strength'],['Hub','Spoke'])
            for entry in self.db('Module').entries():
                ptype = entry['SourceDir']
                importlist = string.split(entry['Imports'],', ')
                for imported in importlist:
                    if imported in ['','None']: continue                    
                    try:
                        idata = ndb.data(imported)
                        ilinks = rje.listIntersect(importlist,string.split(idata['Imported_By'],', '))
                        edb.addEntry({'Hub':entry['Module'],'Spoke':imported,'Evidence':'Imports','Type':'import','ParentType':ptype,'Links':len(ilinks),'Strength':1.0/(len(ilinks)+1)})                        
                    except:
                        self.deBug(idata)
                        self.errorLog('Problem finding imported "%s" module information' % imported)
            ppi.ppiFromEdges()
            ### ~ [1] Use PPI object output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'xgmml' in self.list['ModLinks']:              
                #xfile = '%s.xgmml' % self.getStr('Basefile')
                #xgmml = ppi.ppiToXGMML(name=self.getStr('Basefile'))
                #xgmml.saveXGMML(xfile)
                #!# Might want/need to actually create an XGMML object directly?
                #- Edge = Dictionary of edges between nodes {Type:{(source,target):Attributes}}
                #- EdgeAtt = Dictionary of edge attributes {Att:Type}
                #- Node = Dictionary of Nodes to be output {Node:Attributes}
                #- NodeAtt = Dictionary of node attributes {Att:Type}
                ppi.springXGMML(rpng=False,prejuggle=(10,10),scaleforce=True)
            #!# Add other outputs.            

        except: self.errorLog('Error during PyDocs.ppiModLinks()')
#########################################################################################################################
    ### <5> ### Software distribution Methods                                                                           #
#########################################################################################################################
    def distribute(self,distribution=None,targz=True):  ### Reads in modules and copies relevant files and readmes.
        '''
        Reads in modules and copies relevant files and readmes. Gets key data from documentation.txt and reads as
        commands into new self attributes. Note that, unlike pydocs 1.x, all the modules are now read in to database
        objects in advance of any distributions being made. Objects are therefore not created here. This method will,
        however, re-write the "Imported By" lists for the module subset in the distributions.
        >> distribution:str = Name to look for in file (uses self.info['Name'] if None)
        #distribution_name#
        modules=LIST    : List of modules to document and distribute - adds imported modules.
        otherfiles=LIST : List of other files to copy into distribution directory. (Remember " " if spaces.)
        keywords=LIST   : List of keywords for readme pages.
        >> targz:bool [True] = Whether to tar and zip distribution once made. (Not in windows.)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if distribution: self.setStr({'Name':distribution})
            else: distribution = self.getStr('Name')
            distribution = distribution.lower()
            outdir = rje.makePath(self.getStr('DistDir') + distribution)    # Directory into which to copy/make files.
            images = {}
            ### ~ [1] Read in distribution.txt and make new attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['TabSets'] = []       # List of keys for the WebTabs dictionary of module sets for tabs
            self.list['WebTabs'] = []       # List of programs/modules to have tabs on webpage
            self.dict['WebTabs'] = {}       # Dictionary of programs/module sets to have tabs on webpage
            self.list['Modules'] = []       # List of modules to include in distribution (*.py)
            self.list['Manuals'] = []       # List of manuals include in distribution
            self.list['OtherFiles'] = []    # List of other files to include in distribution
            self.list['Keywords'] = []      # List of keywords for HTML
            self.list['Include'] = []       # List of other packages that are included in download
            self.list['Updates'] = []       # List of updates since last distribution
            dislines = self.loadFromFile('distribution.txt',v=0,checkpath=True,chomplines=True)
            process = False
            for line in dislines:
                if line.find('#%s#' % distribution) == 0: process = True       # Current distribution of choice
                elif line.find('#') == 0: process = False
                elif process: self._cmdReadList(line,'list',['Modules','WebTabs','Include','OtherFiles','Keywords','Manuals'])
                if self.list['WebTabs']:
                    self.debug(self.list['WebTabs'])
                    if '=' in self.list['WebTabs'][0]: (tabset,self.list['WebTabs'][0]) = string.split(self.list['WebTabs'][0],'=')
                    else: tabset = 'Main %s' % self.getStr('Name')
                    self.list['TabSets'].append(tabset)
                    self.dict['WebTabs'][tabset] = self.list['WebTabs'][0:]
                    self.list['WebTabs'] = []
            for tabset in self.list['TabSets']: self.list['WebTabs'] += self.dict['WebTabs'][tabset]
            self.debug(self.dict['WebTabs'])
            self.debug(self.list['WebTabs'])
            #self.deBug(self.list['Modules'])
            for prog in self.list['WebTabs']: self.list['Modules'] += self.db('Module').dataList(self.db('Module').indexEntries('Program',prog),'Module')
            self.list['Modules'] = rje.sortUnique(self.list['Modules'])
            if '*' in self.list['Modules']: self.list['Modules'] = self.db('Module').indexKeys('Module')
            self.deBug(self.list['WebTabs'])
            #self.deBug(self.list['Modules'])
            self.deBug(self.list['Manuals'])
            for pymod in self.list['WebTabs'] + self.list['Modules']:
                if pymod not in self.list['Keywords']: self.list['Keywords'].append(pymod)
                for logobase in ['%s%s%s' % (self.getStr('PyPath'),rje.makePath('docs/'),pymod.lower()),'%s%s%s' % (self.getStr('PyPath'),rje.makePath('docs/'),string.split(pymod.lower(),'_')[0])]:
                    for imgext in ['png','gif','jpg']:
                        for ifile in ['%s_logo.%s' % (logobase,imgext),'%s.%s' % (logobase,imgext)]:
                            if os.path.exists(ifile) and ifile not in self.list['OtherFiles']: self.list['OtherFiles'].append(ifile)
                            if os.path.exists(ifile) and pymod not in images: images[pymod] = ifile
            if self.list['Manuals'] == ['*']:
                self.list['OtherFiles'] += rje.getFileList(self,'%s%s' % (self.getStr('PyPath'),rje.makePath('docs/manuals/')),['*anual.pdf'],subfolders=False,filecount=len(self.list['OtherFiles']))
                self.list['OtherFiles'] += rje.getFileList(self,'%s%s' % (self.getStr('PyPath'),rje.makePath('docs/manuals/')),['*Appendices.pdf'],subfolders=False,filecount=len(self.list['OtherFiles']))
            else:
                for manual in self.list['Manuals'][0:]:
                    if os.path.exists('%smanuals/%s_manual.pdf' % (self.getStr('DocSource'),manual.lower())): self.list['OtherFiles'].append('docs/manuals/%s_manual.pdf' % manual.lower())
                    elif os.path.exists('%smanuals/%s Manual.pdf' % (self.getStr('DocSource'),manual)): self.list['OtherFiles'].append('docs/manuals/%s Manual.pdf' % manual)
                    elif os.path.exists('%smanuals/%s.pdf' % (self.getStr('DocSource'),manual)): self.list['OtherFiles'].append('docs/manuals/%s.pdf' % manual)
                    elif os.path.exists('%s%s Manual.pdf' % (self.getStr('DocSource'),manual)): self.list['OtherFiles'].append('docs/%s Manual.pdf' % manual)
                    elif os.path.exists('%s%s.pdf' % (self.getStr('DocSource'),manual)): self.list['OtherFiles'].append('docs/%s.pdf' % manual)
                    else: self.warnLog('Cannot find PDF Manual for "%s" in %s!' % (manual,self.getStr('DocSource')),'missing_manual'); self.list['Manuals'].remove(manual)
            self.printLog('#DIST','%d modules and %d other files read from distribution.txt for %s' % (len(self.list['Modules']),len(self.list['OtherFiles']),distribution))
            self.deBug(self.list['OtherFiles'])
            ### ~ [2] Update module list list as necessary from import commands. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            modscan = self.list['Modules'][0:]
            while modscan:
                pymod = modscan.pop(0)
                for entry in self.db('Module').indexEntries('Module',pymod):
                    self.printLog('#IMP','Module %s imports: %s' % (pymod,entry['Imports']))
                    for imod in string.split(entry['Imports'],', '):
                        if imod and imod not in self.list['Modules']:
                            self.printLog('#ADD','%s -> %s' % (pymod,imod))
                            if imod in self.db('Module').index('Module'): self.list['Modules'].append(imod); modscan.append(imod)
                            else: self.printLog('#WARN','Module "%s" wants to import module "%s" but not found!' % (pymod,imod))
            ## ~ [2a] ~ Update list of other files to import ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for file in ['docs/gnu_general_public_license.txt','docs/gnu_lesser_general_public_license.txt']:
                if file not in self.list['OtherFiles']: self.list['OtherFiles'].append(file)
            for css in self.list['StyleSheets'] + ['tabber.js']:
                file = 'docs/%s' % css
                if file not in self.list['OtherFiles']: self.list['OtherFiles'].append(file)
            self.printLog('#DIST','%d modules and %d other files compiled for %s' % (len(self.list['Modules']),len(self.list['OtherFiles']),distribution))
            ## ~ [2b] ~ INI Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for pymod in self.list['Modules'][0:]:
                #if os.path.exists('%s%s%s.ini' % (self.getStr('PyPath'),rje.makePath('libraries/ini/'),pymod.lower())):
                #    self.list['OtherFiles'].append('%s%s.ini' % (rje.makePath('ini/'),pymod.lower()))
                if pymod.lower() == 'rje': pymod = 'defaults'
                if os.path.exists('%s%s%s.ini' % (self.getStr('PyPath'),rje.makePath('defaults/'),pymod.lower())):
                    self.list['OtherFiles'].append('%s%s.ini' % (rje.makePath('defaults/'),pymod.lower()))
            ## ~ [2c] ~ Add R files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for rfile in glob.glob('%s%s*.r' % (self.getStr('PyPath'),rje.makePath('libraries/r/'))):
                self.list['OtherFiles'].append('%s%s' % (rje.makePath('r/'),os.path.basename(rfile)))
            #self.deBug(self.list['Modules'])
            self.deBug(self.list['OtherFiles'])
            if not self.list['Modules']: return self.errorLog('Cannot process %s' % distribution,printerror=False)
            ### ~ [3] Generate distributions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.test(): self.deBug('Kill me now if you don\'t want the distribution.')
            ## ~ [3a] Make directory and delete current contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            prevmod = {}
            if os.path.exists(outdir):
                modfile = '%s%s.modules.tdt' % (outdir,distribution)
                if os.path.exists(modfile): prevmod = rje.dataDict(self,modfile,['Module'],['Module','Version','Last Edit','Description'])
                self.deBug(prevmod)
                if not self.test(): rje.deleteDir(self,outdir,contentsonly=True,confirm=self.i()>0 or self.getBool('DeBug'))    # Clear contents
            else:
                if not self.test(): rje.mkDir(self,outdir)
            mx = fx = 0
            ## ~ [3b] Copy modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for pyfile in self.sortedPyFiles(self.list['Modules']):
                entry = self.db('Module').data(pyfile)
                newfile = rje.makePath(outdir + entry['SourceDir']) + entry['Module'] + '.py'
                if self.test(): self.printLog('#COPY','%s -> %s' % (pyfile,newfile)); mx += 1; continue
                rje.mkDir(self,newfile)
                try: shutil.copy(pyfile, newfile); self.printLog('#COPY','%s -> %s' % (pyfile,newfile)); mx += 1
                except: self.errorLog('Cannot copy "%s" to "%s"!' % (pyfile, newfile))
            ## ~ [3c] Copy other files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['OtherFiles'] = rje.sortUnique(string.split(string.join(self.list['OtherFiles'],','),','))
            readme_others = []  # List of other files for readme file
            for file in self.list['OtherFiles'][0:]:
                #self.deBug(file)
                destination_dir = 'libraries/'
                for altdir in ['docs/','data/','settings/','docs/manuals/']:
                    if file.find(altdir) >= 0: destination_dir = altdir
                for altdir in ['defaults/']:    
                    if file.find(altdir) >= 0: destination_dir = 'settings/'
                for altdir in ['ini/','r/']:    
                    if file.find(altdir) == 0:
                        destination_dir = 'libraries/%s' % altdir
                        file = 'libraries/%s' % file
                (fromfile,tofile) = (rje.makePath(file,wholepath=True),rje.makePath('%s%s%s' % (outdir,destination_dir,os.path.basename(file)),wholepath=True))
                #self.deBug('#COPY','%s -> %s' % (fromfile,tofile))
                if self.test(): self.printLog('#COPY','%s -> %s' % (fromfile,tofile)); continue
                if os.path.exists(tofile): self.printLog('#DUP','File "%s" being created twice! %s not copied' % (tofile,fromfile)); continue
                if not os.path.exists(fromfile): fromfile = self.getStr('PyPath') + fromfile
                if not os.path.exists(fromfile): self.printLog('#ERR','File "%s" not found! Not copied' % (fromfile)); continue
                readme_others.append(rje.makePath('./%s%s' % (destination_dir,os.path.basename(file)),wholepath=True))
                rje.mkDir(self,tofile)
                try: shutil.copy(fromfile, tofile); self.printLog('#COPY','%s -> %s' % (fromfile,tofile)); fx += 1
                except: self.errorLog('Cannot copy "%s" to "%s"!' % (fromfile, tofile),quitchoice=False)
            self.printLog('#COPY','%s modules and %d other files copied into %s' % (mx,fx,outdir))
            ## ~ [3d] Save cut-down module table containing version numbers & identify updates ~~~~ ##
            updates = {}
            modout = '%s%s.modules.tdt' % (outdir,distribution)
            modhead = ['Module','Version','Last Edit','Description']
            if not self.test(): rje.delimitedFileOutput(self,modout,modhead)
            for pyfile in self.sortedPyFiles(self.list['Modules']):
                entry = self.db('Module').data(pyfile)
                if not self.test(): rje.delimitedFileOutput(self,modout,modhead,datadict=entry)
                mod = entry['Module']
                history = []
                for mentry in self.db('Method').indexEntries('File',pyfile):
                    if mentry['Method'] == 'history': history = string.split(mentry['DocString'],'\n'); break
                if mod not in prevmod or prevmod[mod]['Version'] != entry['Version']:
                    if mod in prevmod:
                        updates[mod] = [prevmod[mod]['Version'],entry['Version']]
                        self.printLog('#VNUM','%s: Version %s -> Version %s' % (mod,updates[mod][0],updates[mod][1]))
                        self.list['Updates'].append('<p><b>&bull; %s:</b> <i>Updated from Version %s.</i>' % (mod,updates[mod][0]))
                    else:
                        updates[mod] = ['-',entry['Version']]
                        self.printLog('#VNUM','%s: Creation -> Version %s' % (mod,updates[mod][1]))
                        self.list['Updates'].append('<p><b>&bull; %s:</b> <i>Created/Renamed.</i>' % mod)
                    #!# Extract from module history the relevant update information                    
                    updated = not mod in prevmod
                    lastv = ''
                    for dline in history:
                        if rje.matchExp('# (\d+\.\d+\.?\d*)\s?-\s(\S.+)$',dline): lastv = rje.matchExp('# (\d+\.\d+\.?\d*)\s?-\s(\S.+)$',dline)[0]
                        self.deBug(rje.matchExp('# (\d\.\d+\.?\d*)\s?-\s(\S.+)$',dline))
                        if not updated and rje.matchExp('# %s\s?-\s(\S.+)$' % prevmod[mod]['Version'],dline): updated = True
                        elif updated and rje.matchExp('# (\d\.\d+\.?\d*)\s?-\s(\S.+)$',dline):
                            (v,text) = rje.matchExp('# (\d\.\d+\.?\d*)\s?-\s(\S.+)$',dline)
                            self.printLog('#V%s' % v,text)
                            updates[mod].append('Version %s: %s' % (v,text))                    
                            self.list['Updates'].append('<br>&rarr; Version %s: %s' % (v,text))
                    self.list['Updates'].append('</p>')
                    if lastv != entry['Version']: self.errorLog('Module %s Version %s but history() ends at %s' % (mod,entry['Version'],lastv),printerror=False)
                if self.test(): self.deBug('>>>')
            self.printLog('#VNUM','%d modules identify version updates' % len(updates))
            self.deBug('Update End')
            if self.test(): return
            ## ~ [3e] Make readme.html and save in distribution directory ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            readme_others.sort()
            self.list['OtherFiles'] = readme_others
            csspath = self.getStr('StylePath')
            #if self.dev():
            self.setStr({'StylePath':'../docs/'})
            self.htmlReadMe(rje.makePath(outdir)+'readme/',modlist=[],zipfile='',level=0,onefile=False)
            #else:
            #    self.setStr({'StylePath':'./docs/'})
            #    self.htmlReadMe(rje.makePath(outdir),self.list['Modules'],distribution,linkothers=True)
            self.setStr({'StylePath':csspath})
            ### ~ [4] Tar and Zip the distribution if possible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Win32'): self.printLog('#TARGZ','Cannot tar and gzip in Windows!')
            else:
                mydir = os.path.abspath(os.curdir)
                self.printLog('#TAR','tar -czf %s.%s.tar.gz %s/' % (distribution,self.getStr('Release'),distribution))
                try:
                    os.chdir(self.getStr('DistDir'))
                    if os.system('tar -czf %s.%s.tar.gz %s/' % (distribution,self.getStr('Release'),distribution)): raise ValueError
                    else: os.chdir(mydir); self.printLog('#DIST','Distribution package for %s made and tarred/zipped.' % distribution)
                except: os.chdir(mydir); self.errorLog('Something went wrong tarring and zipping %s package' % distribution)
                os.chdir(mydir)
            ### ~ [5] Make webpage for online version of package ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [5a] Setup general HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #htmlfile = rje.makePath('%s%s.html' % (outdir,distribution),True)
            htmlfile = '%sindex.html' % outdir
            linkhtml = ''
            if os.path.exists('%slinks.html' % self.getStr('WebDir')): linkhtml += open('%slinks.html' % self.getStr('WebDir'),'r').read()
            gnu = '<A HREF="./docs/gnu_general_public_license.txt">GNU General Public License</A>'
            ## ~ [5b] HTML headers, including keywords and tabber code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            title = '%s software package' % self.getStr('Name')   #!# Get better title
            stylesheets = []
            for css in self.list['StyleSheets']: stylesheets.append(self.getStr('StylePath')+css)
            htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=self.list['Keywords'],javascript=self.getStr('StylePath'))
            ## ~ [5c] Output sidebar links, read in from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if linkhtml:
                htmlhead += '<TABLE width=100% border=0><TR width=100%><TD width=20% valign=top align=LEFT>\n\n'
                htmlhead += linkhtml
                htmlhead += '\n\n</TD><TD width=3% valign=top></TD><TD width=77%  valign=top>\n\n'
            ## ~ [5d] Generate general text, linking to output and help ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htmlbody = '<a name="Top"></a><h1>%s</h1>\n\n' % title
            htmlbody += '<p>This webpage contains the latest <a href="../%s.%s.tar.gz">download</a> and documentation for the <scap>%s</scap> package.' % (distribution,self.getStr('Release'),self.getStr('Name'))
            if self.list['WebTabs']:
                proglist = string.join(self.list['WebTabs'][:-1],', ')
                if len(self.list['WebTabs']) > 1: proglist += ' and %s' % self.list['WebTabs'][-1]
                else: proglist = self.list['WebTabs'][-1]
                #htmlbody += 'The main programs in %s are: %s. Please see the tabs below for more details. Release notes can be found in the <b>ReadMe</b> tab.' % (self.getStr('Name'),proglist)
                htmlbody += 'Please see the tabs below for more details. Release notes can be found in the <b>ReadMe</b> tab.'
                if distribution in ['SLiMSuite','SeqSuite']: htmlbody += '\nMore information can also be found at the <a href="http://slimsuite.blogspot.co.uk">SLiMSuite blog</a>.'
                htmlbody += '</p>\n\n'
            else: htmlbody += ' %s contains a number of functional Python modules, which are listed in the <b>ReadMe</b> tab below along with release notes.</p>\n\n' % self.getStr('Name')
            if len(self.list['Include']) > 1:
                htmlbody += '<p>%s also includes the programs in the ' % self.getStr('Name')
                htmlbody += '%s and %s packages. \nSee ' % (string.join(self.list['Include'][:-1],', '),self.list['Include'][-1])
                for include in self.list['Include'][:-1]: htmlbody += '<a href="../%s/">%s page</a>, ' % (include.lower(),include)
                include = self.list['Include'][-1]
                htmlbody += 'and <a href="../%s/">%s page</a> for more details.</p>\n\n' % (include.lower(),include)
            elif self.list['Include']:
                include = self.list['Include'][0]
                htmlbody += '<p>%s also includes the programs in the %s package. See <a href="../%s/">%s page</a> for more details.</p>\n\n' % (self.getStr('Name'),include,include.lower(),include)
            htmltabs = [(self.getStr('Name'), 'The current %s release is %s. Click on the tabs for details.' % (self.getStr('Name'),self.getStr('Release')), self.getStr('Name'))]
            tabtext = '<h3>Availability</h3>\n'
            tabtext += '<p>The tools in <scap>%s</scap> are freely available for local installation under a %s as part of the ' % (self.getStr('Name'),gnu)
            tabtext += '<a href="../%s.%s.tar.gz">%s package</a> (%s release). Please see the Manuals and <a href="./readme/readme.html">ReadMe</a> for more details.</p>\n\n' % (distribution,self.getStr('Release'),self.getStr('Name'),self.getStr('Release'))
            if self.getStr('EMail'): tabtext += '<P>To contact the author, e-mail: <A HREF="mailto:%s">%s</A></P>' % (self.getStr('EMail'),self.getStr('EMail'))
            htmltabs.append(('Availability',tabtext,'%s availability' % self.getStr('Name')))
            tabtext = '<h3>Citing %s</h3>\n' % self.getStr('Name')
            tabtext += '<p>When publishing analyses performed with this software, please cite the individual papers listed for the relevant program/module. '
            tabtext += 'If no paper or URL is listed, please cite this website.</p>\n\n'
            htmltabs.append(('How to cite',tabtext,'Citing %s' % self.getStr('Name')))
            ## ~ [5e] ~ Generate Intro tab, read in from file (if found) ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s%s.html' % (self.getStr('WebDir'),distribution)):
                htmltabs.append(('Summary',open('%s%s.html' % (self.getStr('WebDir'),distribution)).read(),'Introduction to the %s package' % self.getStr('Name')))
            ## ~ [5f] ~ Generate ReadMe tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if self.dev():
            htmltabs.append(('ReadMe',self.htmlReadMe(rje.makePath(outdir)+'readme/',self.list['Modules'],distribution,level=1,linkothers=True,onefile=False),'Full ReadMe for the %s package' % self.getStr('Name')))
            #else:
            #    htmltabs.append(('ReadMe',self.htmlReadMe(rje.makePath(outdir),self.list['Modules'],distribution,level=1,linkothers=True),'Full ReadMe for the %s package' % self.getStr('Name')))
            if '*' in self.list['Manuals']: self.list['Manuals'].remove('*')
            for file in self.list['OtherFiles']:
                if file[-4:].lower() != '.pdf' or file in self.list['Manuals']: continue
                for manpath in ['docs\/manuals','docs']:
                    manual = rje.matchExp('%s\/(\.+) Manual.pdf' % manpath, file)
                    if manual and manual[0] not in self.list['Manuals']: self.list['Manuals'].append(manual[0])
                #if rje.matchExp('docs\/(\.+) Manual.pdf', file) and rje.matchExp('docs\/(\.+) Manual.pdf', file)[0] not in self.list['Manuals']: self.list['Manuals'].append(file)
                #if rje.matchExp('docs\/manuals\/(\.+) Manual.pdf', file) and rje.matchExp('docs\/manuals\/(\.+) Manual.pdf', file)[0] not in self.list['Manuals']: self.list['Manuals'].append(file)
            if self.list['Manuals']:
                docs = ['<H3>PDF Manuals</H3>','<p>Please note that not all programs have manuals and some manuals are out of date.',
                        'If in doubt, check the readme information for the latest options and default settings.',
                        'Please report any anomalous behaviour. Suggestions for improvements to programs and documentation are also appreciated.</p>',
                        '','<UL>']
                for manual in self.list['Manuals'][0:]:
                    try: desc = self.db('Module').dataList(self.db('Module').indexEntries('Program',manual),'Description')[0]
                    except: desc = 'Additional information'
                    if os.path.exists('%smanuals/%s_manual.pdf' % (self.getStr('DocSource'),manual.lower())): file = './docs/manuals/%s Manual.pdf' % manual.lower()
                    elif os.path.exists('%smanuals/%s Manual.pdf' % (self.getStr('DocSource'),manual)): file = './docs/manuals/%s Manual.pdf' % manual
                    elif os.path.exists('%smanuals/%s.pdf' % (self.getStr('DocSource'),manual)): file = './docs/manuals/%s.pdf' % manual
                    elif os.path.exists('%s%s Manual.pdf' % (self.getStr('DocSource'),manual)): file = './docs/%s Manual.pdf' % manual
                    else: file = './docs/%s.pdf' % manual
                    docs.append('<LI><a href="%s" TARGET="_blank">%s</a> - %s' % (file,file,desc))
                docs.append('</UL>')
                htmltabs.append(('Manuals',string.join(docs,'\n'),'Documentation (PDF Manuals)'))
            htmlbody += rje_html.tabberHTML(self.getStr('Name'),htmltabs,level=0)
            ## ~ [5g] ~ Generate Tab per named module in WebTabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for tabset in self.list['TabSets']:
                htmlbody += '\n<h3>%s programs/modules</h3>\n' % tabset
                tabtext = 'The %s programs/modules in %s are: %s. Click on the tabs, read the manuals, or see the <a href="./readme/readme.html">ReadMe</a> for more details.' % (tabset,self.getStr('Name'),string.join(self.dict['WebTabs'][tabset],', '))
                if self.getStr('Name') in ['SLiMSuite','SeqSuite']: htmlbody += '\nMore information can also be found at the <a href="http://slimsuite.blogspot.co.uk">SLiMSuite blog</a>.'
                htmltabs = [('^', tabtext, '%s programs/modules' % tabset)]
                for pymod in self.dict['WebTabs'][tabset][0:]:
                    tabname = pymod
                    try: pymod = self.db('Module').dataList(self.db('Module').indexEntries('Program',tabname),'Module')[0]
                    except: self.printLog('#WARN','Cannot find module for "%s"' % tabname)
                    webtab = []
                    for hfile in rje.getFileList(self,self.getStr('WebDir'),['%s*.html' % tabname.lower()],subfolders=False,summary=False):
                        hfilename = rje.baseFile(hfile,True)
                        if hfilename.startswith('slimsuite'): continue  # Interference from webserver files!
                        webtab.append((hfilename,open(hfile,'r').read(),open(hfile,'r').readlines()[1]))    # Make sure line 2 is the title
                    ix = len(webtab)
                    pyfile = self.findModule(pymod)
                    if pyfile:
                        webmod = self.modHTML(pyfile,level=1,modlist=self.list['Modules'],reduced=True)
                        if webmod:
                            webtab += webmod
                            if tabname in images:
                                ifile = './docs/%s' % os.path.basename(images[tabname])
                                newtop = webtab[ix][1]
                                newtop = string.split(newtop,'Copyright &copy;')
                                newtop = '<table><tr align="top"><td width="300"><center><img height="120" src="%s"></center></td><td>\n%s</td></tr></table>\n\nCopyright &copy;%s' % (ifile,newtop[0],string.join(newtop[1:],'Copyright &copy;'))
                                webtab[ix] = (webtab[ix][0],newtop,webtab[ix][2])
                    elif pymod != '*': self.printLog('#WARN','Cannot find pyfile for "%s" (%s)' % (pymod,tabname))
                    if webtab: htmltabs.append((tabname,rje_html.tabberHTML(tabname,webtab,level=1),'Details for %s program of %s package' % (tabname,self.getStr('Name'))))
                htmlbody += rje_html.tabberHTML('%s programs/modules' % tabset,htmltabs,level=0)
            ## ~ [5h] End of HTML code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htmltail = rje_html.htmlTail('%s %s' % (self.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1]))
            if linkhtml: htmltail = '</TD></TR></TABLE>\n\n' + htmltail
            ## ~ [5i] Save and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            open(htmlfile,'w').write(htmlhead+htmlbody+htmltail)
            self.printLog('#HTML','HTML summary for %s output to %s' % (self.getStr('Name'),htmlfile))
        except: self.errorLog('Error in PyDoc.distribute(%s)' % self.getStr('Name'),printerror=True,quitchoice=True)
#########################################################################################################################
    ### <6> ### Webpage generation methods                                                                              #
#########################################################################################################################
    def makePackageIndex(self): ### Makes a summary info page for the distributions listed
        '''Makes a summary info page for the distributions listed.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self._cmdList()     # Re-read general commandline options
            if not self.list['Distribute']: return False
            db = self.db('Module')
            for pymod in db.indexKeys('Module') + db.indexKeys('Program'):
                if pymod not in self.list['Keywords']: self.list['Keywords'].append(pymod)
            for str in ['Author','Name','EMail']:
                if self.getStr(str) not in self.list['Keywords']: self.list['Keywords'].append(self.getStr(str))
            outdir = self.getStr('DistDir')
            htmlfile = '%sindex.html' % outdir
            ### ~ [1] Create HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            linkhtml = ''
            if os.path.exists('%slinks.html' % self.getStr('WebDir')): linkhtml += open('%slinks.html' % self.getStr('WebDir'),'r').read()
            title = 'Distributed software packages'     #!# Get better title
            stylesheets = []
            for css in self.list['StyleSheets']: stylesheets.append(self.getStr('StylePath')+css)
            htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=self.list['Keywords'],javascript=self.getStr('StylePath'))
            ## ~ [1a] Output sidebar links, read in from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if linkhtml:
                htmlhead += '<TABLE width=100% border=0><TR width=100%><TD width=20% valign=top align=LEFT>\n\n'
                htmlhead += linkhtml
                htmlhead += '\n\n</TD><TD width=3% valign=top></TD><TD width=77%  valign=top>\n\n'
            ## ~ [1b] Generate general text, linking to output and help ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htmlbody = '<h1>%s</h1>\n\n' % title
            htmlbody += '<p>This webpage contains links, downloads and documentation for the '
            if len(self.list['Distribute']) > 1:
                for distribution in self.list['Distribute'][:-1]:
                    htmlbody += '<a href="./%s/">%s</a> (<a href="./%s.%s.tar.gz">download</a>), ' % (distribution.lower(),distribution,distribution.lower(),self.getStr('Release'))
                distribution = self.list['Distribute'][-1]
                htmlbody += 'and <a href="./%s/">%s</a> (<a href="./%s.%s.tar.gz">download</a>) packages.</p>' % (distribution.lower(),distribution,distribution.lower(),self.getStr('Release'))
            else: 
                distribution = self.list['Distribute'][-1]
                htmlbody += '<a href="./%s/">%s</a> (<a href="./%s.%s.tar.gz">download</a>) package.</p>' % (distribution.lower(),distribution,distribution.lower(),self.getStr('Release'))
            gnu = '<A HREF="./%s/docs/gnu_general_public_license.txt">GNU General Public License</A>' % distribution.lower()
            htmlbody += '<p>These tools are freely available for local installation under a %s. ' % (gnu)
            htmlbody += 'Please see the Manuals and ReadMe for more details.</p>\n\n'
            if self.getStr('EMail'): htmlbody += '<P>To contact the author, e-mail: <A HREF="mailto:%s">%s</A></P>' % (self.getStr('EMail'),self.getStr('EMail'))
            htmlbody += '\n<h3>Citing %s</h3>\n' % self.getStr('Name')
            htmlbody += '<p>When publishing analyses performed with this software, please cite the individual papers listed for the relevant program/module.'
            htmlbody += 'If no program is listed, please cite this website.</p>\n\n'          
            ### ~ [2] ~ Generate Summary tabs, read in from files (if found) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            htmltabs = []
            for distribution in self.list['Distribute']:
                if os.path.exists('%s%s.html' % (self.getStr('WebDir'),distribution.lower())):
                    htmltabs.append((distribution,open('%s%s.html' % (self.getStr('WebDir'),distribution.lower())).read(),'Introduction to the %s package' % distribution))
            htmlbody += rje_html.tabberHTML(self.getStr('Name'),htmltabs,level=0)
            ## ~ [5h] End of HTML code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htmltail = rje_html.htmlTail('%s %s' % (self.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1]))
            if linkhtml: htmltail = '</TD></TR></TABLE>\n\n' + htmltail
            ## ~ [5i] Save and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            open(htmlfile,'w').write(htmlhead+htmlbody+htmltail)
            self.printLog('#HTML','HTML summary for packages output to %s' % (htmlfile))
        except: self.errorLog('Error in PyDoc.makePackageIndex(%s)' % self.getStr('Name'),printerror=True,quitchoice=True)
#########################################################################################################################
    ### <7> ### REST Server documentation page methods                                                                  #
#########################################################################################################################
    def formatDocString(self,docstring,outfmt=False):    ### Adds <code> formatting and page links to docstring text
        '''Adds <code> formatting and page links to docstring text.'''
        try:### ~ [1] Format docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            docstring = string.replace(docstring,'<','&lt;')
            docstring = string.replace(docstring,'>','&gt;')
            docstring = string.replace(docstring,'&lt;=','&lt;!!*!!')
            docstring = string.replace(docstring,'&gt;=','&gt;!!*!!')
            while rje.matchExp('(http:\S+)=(\S+)',docstring):
                http = rje.matchExp('(http:\S+)=(\S+)',docstring)
                docstring = string.replace(docstring,'%s=%s' % (http[0],http[1]),'%s!!*!!%s' % (http[0],http[1]))
            while rje.matchExp('(http:\S+)',docstring):
                http = rje.matchExp('(http:\S+)',docstring)[0]
                newhttp = string.replace(http,'http','!!HTTP!!')
                docstring = string.replace(docstring,http,'<a href!!*!!"%s">%s</a>' % (newhttp,newhttp))
            if outfmt:
                while rje.matchExp('(\S+) =',docstring):    # opt=value
                    cmd = rje.matchExp('(\S+) =',docstring)[0]
                    docstring = string.replace(docstring,'%s =' % cmd,'<code>%s</code> !!*!!' % (cmd),1)
                while rje.matchExp('(&\S+=\S*\w)',docstring):    # opt=value
                    cmdstr = rje.matchExp('(&\S+=\S*\w)',docstring)[0]
                    [cmd,opt] = string.split(cmdstr,'=',maxsplit=1)
                    opt = string.replace(opt,'=','!!*!!')   # For REST calls etc.
                    docstring = string.replace(docstring,cmdstr,'<code>%s!!*!!%s</code>' % (cmd,opt),1)
                pretext = rje.matchExp('(###[~]+###.*\n# OUTFMT\:.*\n\.\.\..*)',docstring)
                if pretext: docstring = string.replace(docstring,pretext[0],'<pre>%s</pre>' % string.replace(pretext[0],'\n','<PREND>'))
            while rje.matchExp('(`\w\S*=\S*\w`)',docstring):    # `opt=value`
                cmdstr = rje.matchExp('(`\w\S*=\S*\w`)',docstring)[0]
                [cmd,opt] = string.split(string.replace(cmdstr,'`',''),'=')
                docstring = string.replace(docstring,cmdstr,'<code>[%s!!*!!%s]{cmd:%s}</code>' % (cmd,opt,cmd),1)
            while rje.matchExp('(\w\S*=\S*\w)',docstring):    # opt=value
                cmdstr = rje.matchExp('(\w\S*=\S*\w)',docstring)[0]
                try:
                    [cmd,opt] = string.split(cmdstr,'=')
                    docstring = string.replace(docstring,cmdstr,'<code>[%s!!*!!%s]{cmd:%s}</code>' % (cmd,opt,cmd),1)
                except: docstring = string.replace(docstring,cmdstr,'<code>%s</code>' % string.replace(cmdstr,'=','!!*!!'),1)
            while rje.matchExp('(({cmd[^\n]+)\[([^<\r\n]+)\][\r\n])',docstring):
                codestr = rje.matchExp('(({cmd[^\n]+)\[([^<\r\n]+)\][\r\n])',docstring)
                docstring = string.replace(docstring,codestr[0],'%s[<code>%s</code>]\n' % (codestr[1],codestr[2]),1)
            while rje.matchExp('(`[^`]+`)',docstring):    # `code`
                codestr = rje.matchExp('(`[^`]+`)',docstring)[0]
                docstring = string.replace(docstring,codestr,'<code>%s</code>' % (string.replace(codestr,'`','')),1)
            while rje.matchExp('(\S+)\.py',docstring):
                docmod = rje.matchExp('(\S+)\.py',docstring)[0]
                if self.findModule(docmod): docstring = string.replace(docstring,'%s.py' % docmod,'[%s!!.!!py]{module:%s}' % (docmod,docmod),1)
                else: docstring = string.replace(docstring,'%s.py' % docmod,'%s!!.!!py' % docmod,1)
        except:
            self.errorLog('formatDocString error',quitchoice=False)
            docstring += '<p><i>Something went wrong with docstring formatting!</i></p>'
        ### ~ [2] Tidy and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        docstring = string.replace(docstring,'!!HTTP!!','http')
        docstring = string.replace(docstring,'!!*!!','=')
        docstring = string.replace(docstring,'!!.!!','.')
        docstring = string.replace(docstring,'</code><code>','')
        if outfmt:
            docstring = docstring.replace('\n','<br>\n')
            docstring = docstring.replace('<PREND>','\n')
        return docstring
#########################################################################################################################
    def parseToDocTabs(self,pymod):    ### Parses the docstring etc. for module and returns HTML for REST docs page
        '''
        Parses the docstring etc. for module and returns HTML for REST docs page.
        '''
        data = {}
        try:### [1] Parse key module contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pyfile = self.findModule(pymod)
            if not pyfile: pyfile = self.findModule('rje_%s' % pymod)
            if not pyfile: return {'DocString':'Module %s not found.' % pymod}
            data = self.addModule(pyfile,docsonly=True)
            for dkey in ['FullText','Classes','Imported_By']: data.pop(dkey)
            # ['Author', 'Description', 'DocString', 'File', 'Imports', 'Last Edit', 'Methods', 'Module', 'Program', 'SourceDir', 'Version']
            ### [2] Format commandline options for {cmd:page} links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug(data['DocString'])
            docstring = self.formatDocString(data['DocString'])
            if 'Uses general modules:' in docstring: docstring = docstring[:docstring.find('Uses general modules:')]
            ### [3] Convert docstring into tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug(docstring)
            docstring = string.split(docstring,'\n')
            #self.debug(docstring)
            ## [3a] Summary Information Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            infotable = '<table border=0 width="100%" id="linktab">\n'
            tablist = []
            tabid = data['Program']
            tabhtml = 'Empty Docstring!'
            tabdesc = 'Summary info for %s V%s' % (data['Program'],data['Version'])
            while docstring:
                nextline = rje.chomp(docstring.pop(0))
                if not nextline: continue
                if ':' in nextline:
                    match = [nextline[:nextline.find(':')+1],nextline[nextline.find(':')+1:]]
                    if string.split(match[0])[-1] in ['Cite:','Citation:']:
                        if rje.matchExp('^%s\s+(\S.+)\s+\[PMID:\s*(\d+)\]' % match[0],nextline):
                            pmid = rje.matchExp('^%s\s+(\S.+)\s+\[PMID:\s*(\d+)\]' % match[0],nextline)
                            match = (match[0],' <a href="http://www.ncbi.nlm.nih.gov/pubmed/%s?dopt=Abstract" TARGET="_blank">%s</a>' % (pmid[1],pmid[0]))
                        infotable += '<tr><th id="linktab">%s</th><td id="doctab2">%s</td></tr>\n' % (match[0],match[1])
                    elif string.split(match[0])[-1] in ['Program:']:
                        infotable += '<tr><th id="linktab">%s</th><td id="doctab3">%s</td></tr>\n' % (match[0],match[1])
                    else:
                        # http now managed by formatDocString.
                        #if rje.matchExp('(http:\S+)',match[1]):
                        #    http = rje.matchExp('(http:\S+)',match[1])[0]
                        #    match[1] = string.replace(match[1],http,'<a href="%s" TARGET="_blank">%s</a>' % (http,http))
                        infotable += '<tr><th id="linktab">%s</th><td id="doctab">%s</td></tr>\n' % (match[0],match[1])
                elif 'Copyright' in nextline:
                    infotable += '</table>\n'
                    logo = '%sEdwardsLabLogo_300px.jpg' % self.getStr('LogoURL')
                    for img in ['png','gif','jpg']:
                        logofile = '%s_logo.%s' % (data['Program'].lower(),img)
                        if rje.exists('%slogos/%s' % (self.getStr('DocSource'),logofile)):
                            logo = '%s%s' % (self.getStr('LogoURL'),logofile); break
                        #else: print 'No >> %slogos/%s' % (self.getStr('DocSource'),logofile)
                    tabhtml = '<table width="100%" border=0>\n'
                    tabhtml += '<tr><td width=150><img src="%s" width=150 border=0></td>\n' % logo
                    tabhtml += '<td width=50><p></p></td>\n'
                    tabhtml += '<td>%s</td>\n' % infotable
                    tabhtml += '</tr></table>\n<!-- End of Summary Table -->\n'
                    tabhtml += '<p>%s' % string.replace(nextline,'(C)','&copy;')
                    if data['Imports']:
                        tabhtml += '</p>\n<hr>\n'
                        tabhtml += '<p><b>Imported modules:</b> \n'
                        for imod in data['Imports']: tabhtml += '<code>[%s]{module:%s}</code> \n' % (imod,imod)
                    tabhtml += '</p>\n<hr>\n<p>See <a href="http://slimsuite.blogspot.com/" target="_blank">SLiMSuite Blog</a>'
                    manfile = '%smanuals/%s_manual.pdf' % (self.getStr('DocSource'),data['Program'].lower())
                    if rje.exists(manfile): tabhtml += ' and <a href="%s%s_manual.pdf" target="_blank"><code>%s Manual</code></a>' % (self.getStr('ManualURL'),data['Program'].lower(),data['Program'])
                    tabhtml += ' for further documentation.'
                    break
            ## [3b] Parse rest of docstring into tabs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while docstring:
                #self.debug(docstring[:3])
                nextline = rje.chomp(docstring.pop(0))
                newtab = rje.matchExp('^(\S.+):',nextline)
                if newtab and not newtab[0].endswith('{module') and not newtab[0].startswith('#'):
                    tabhtml = self.tabHTMLCleanup('%s</p>' % tabhtml)
                    tablist.append((tabid,tabhtml,tabdesc))
                    tabid = newtab[0]
                    tabhtml = '<h2>%s</h2>\n<p>' % (tabid)
                    tabdesc = 'Docstring %s for %s module' % (tabid,pymod)
                else:
                    #nextline = string.join(string.split(nextline))
                    if rje.matchExp('^\s*(#+)\s~?\s?(\S[^~]+)(:|\s~)',nextline):
                        subhead = rje.matchExp('^\s*(#+)\s~?\s?(\S[^~]+)(:|\s~)',nextline)
                        nextline = '<h%d>%s</h%d>' % (len(subhead[0]),subhead[1],len(subhead[0]))
                    elif rje.matchExp('^\s*(#+\s~+\s#+)',nextline): nextline = '<hr>'
                    if nextline.startswith('    '): nextline = nextline[4:]
                    #!# Need to add ability to deal with wrapped bullets
                    if nextline.startswith('* '):
                        try:
                            if not string.split(tabhtml)[-1].endswith('</li>') and not string.split(tabhtml)[-1].endswith('</li><br>'): tabhtml += '<ul>'
                        except: tabhtml += '<ul>'
                        nextline = '<li>%s</li>' % nextline
                        try:
                            if string.split(docstring[0])[0].startswith('* '): nextline += '</ul>'
                        except: nextline += '</ul>'
                        if rje.matchExp('(\* (\S+)\s+=)',nextline):
                            bullet = rje.matchExp('(\* (\S+)\s+=)',nextline)
                            nextline = string.replace(nextline,bullet[0],'<code>%s</code> =' % bullet[1])
                    if not rje.matchExp('(\S)',nextline): nextline = ''
                    elif nextline.startswith('<h'): nextline += '<p>'
                    elif nextline.endswith(']'):
                        nextline += '<br>'
                        if docstring and not docstring[0]: docstring.pop(0)
                    elif string.split(tabhtml)[-1].endswith('</li>'): pass
                    elif tabid.lower() not in ['function','summary','introduction']: nextline += '<br>'
                    tabhtml += '%s\n' % nextline
                    #self.debug(nextline)
            tabhtml = self.tabHTMLCleanup('%s</p>' % tabhtml)
            tablist.append((tabid,tabhtml,tabdesc))
            ## [3c] Add History tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabid = 'History'
            tabhtml = '<h3>%s Module Version History</h3>\n' % tabid
            try:
                for (method,mdesc,mdetail) in data['Methods']:
                    if method != 'history': continue
                    tabhtml += '<pre>%s</pre>\n\n' % string.join(mdetail[2:-2],'')
                    tablist.append((tabid,tabhtml,tabdesc))
            except: self.errorLog('History')
            ## [3d] Add REST Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            restmod = pymod
            if restmod[-3:-1] == '_V': restmod = restmod[:-3]
            try:
                resturl = '%s%s&rest=outfmt' % (self.getStr('RestURL'),restmod)
                outfmt = self.formatDocString(urllib2.urlopen(resturl).read(),outfmt=True)
            except: outfmt = self.errorLog('Failure to parse %s outfmt' % restmod)
            if not outfmt.startswith('ERROR') and not outfmt.startswith('#ERR'):
                outfmt = '<h2>%s REST Output formats</h2>\n\n%s' % (data['Program'],outfmt)
                tablist.append(('REST',outfmt,'%s REST Output formats' % data['Program']))
            ## [3e] Add Server Tab ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try:
                # Need *.form file in servers/ directory with code!
                formfile = '%sservers/%s.form' % (self.getStr('PyPath'),data['Program'].lower())
                if rje.exists(formfile):
                    formhtml = open(formfile,'r').read()
                    tablist.append(('Server',formhtml,'Webserver form for running REST server.'))
            except: self.errorLog('parseToDocTabs:Server',quitchoice=False)
            data['DocString'] = rje_html.tabberHTML(pymod,tablist)
        except: data['DocString'] = self.errorLog('Error in PyDoc.parseToDocTabs(%s)' % pymod)
        return data
#########################################################################################################################
    def tabHTMLCleanup(self,tabhtml):   ### Performs final processing/cleanup of documentation tabhtml for parseToDocTabs
        '''Performs final processing/cleanup of documentation tabhtml for parseToDocTabs.'''
        try:### ~ [1] Cleanup rogue paragraphs, headings and line endings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug(tabhtml)
            tabhtml = string.replace(tabhtml,'\n\n','</p>\n\n<p>')
            tabhtml = string.replace(tabhtml,']<br></p>',']<br>')
            tabhtml = string.replace(tabhtml,'<br>\n  \n  <p>','<br>\n')
            tabhtml = string.replace(tabhtml,'<br></p>','</p>')
            tabhtml = string.replace(tabhtml,'<p><h','<h')
            tabhtml = string.replace(tabhtml,'<p></p>','')
            tabhtml = string.replace(tabhtml,'\n</p>','</p>\n\n')
            tabhtml = string.replace(tabhtml,'<br>\n\n<p>','<br>\n')
            #self.debug([tabhtml])
            #self.debug(tabhtml.find('<br>\n  \n  <h'))
            tabhtml = string.replace(tabhtml,'<br>\n\n<h','</p>\n\n<h')
            tabhtml = string.replace(tabhtml,'<br>\nSee also','</p>\n\n<p>See also')
            tabhtml = string.replace(tabhtml,'<p></p>','')
            while '\n\n\n' in tabhtml: tabhtml = string.replace(tabhtml,'\n\n\n','\n\n')
            while tabhtml.endswith('\n'): tabhtml = tabhtml[:-1]
            #self.debug([tabhtml])
            #self.bugPrint(tabhtml)
            #self.debug('<p></p>' in tabhtml)
        except: self.errorLog('PyDocs.tabHTMLCleanup() error')
        return tabhtml
#########################################################################################################################
    def makePages(self):    ### Generates docs/pages/cmd/ content from populated database tables
        '''Generates docs/pages/cmd/ content from populated database tables.'''
        try:### [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pagedir = '%spages/' % self.getStr('DocDir')
            prevdir = '%spages/' % self.getStr('DocSource')
            odb = self.db('Options')    # File	Class	Type	Attribute	Argument	ArgType	Description	Default
            odb.addField('page')
            for entry in odb.entries():
                arg = entry['ArgType']
                entry['page'] = arg
                if arg in ['"X"','X\\t','0','N(,X)', 'T/F/X']: entry['page'] = 'X'
                if arg in [ 'FILE(s)', 'FILE/LIST', 'FILES']: entry['page'] = 'FILELIST'
                if arg in ['T', 'T/F', 'T/F\\t']: entry['page'] = 'boolean'
                if arg in ['X,Y', 'X[,Y]']: entry['page'] = 'minmax'
                if arg in ['X,Y,..','X,Y,..,Z']: entry['page'] = 'LIST'
                if arg in ['PATH/']: entry['page'] = 'PATH'
            ### [1] Cmd Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'%scmd/' % pagedir)
            for arg in odb.index('Argument'):
                if arg == '?': continue
                apage = '%scmd/%s.page' % (pagedir,arg)
                pagetxt = '<h1>SLiMSuite Documentation</h1>\n'
                pagetxt += '<p><i>Please note that command documentation sometimes lags behind the latest module versions.\n'
                pagetxt += 'Click on module names to check current documentation.</i></p>\n'
                pagetxt += '<h2>cmd: %s</h2>\n' % (arg)
                pagetxt += '<p><b>Option Type:</b>'
                for otype in odb.indexDataList('Argument',arg,'ArgType'):
                    if otype == 'T/F': pagetxt += ' <code>[%s]{docs:boolean}</code>' % (otype)
                    else: pagetxt += ' <code>[%s]{docs:%s}</code>' % (otype,otype)
                pagetxt += '</p>\n'
                pagetxt += '<p><b>Modules:</b>'
                for ofile in odb.indexDataList('Argument',arg,'File'):
                    omod = string.split(ofile,'/')[-1][:-3]
                    pagetxt += ' <code>[%s]{module:%s}</code>' % (omod,omod)
                pagetxt += '</p>\n\n'
                #pagetxt += '<!-- INSERT CUSTOM TEXT HERE -->'
                prevpage = '%scmd/%s.page' % (prevdir,arg)
                if os.path.exists(prevpage):
                    prevtext = open(prevpage,'r').read()
                    prevtext = prevtext[prevtext.find('<!-- INSERT CUSTOM TEXT HERE -->'):prevtext.find('<!-- END CUSTOM TEXT -->')]
                    if not prevtext: prevtext = '<!-- INSERT CUSTOM TEXT HERE -->\n\n'
                else: prevtext = '<!-- INSERT CUSTOM TEXT HERE -->\n\n'
#                    prevtext = string.replace(open(prevpage,'r').read(),'\n','!!!N!!!')
#                    prevtext = rje.matchExp('<!-- INSERT CUSTOM TEXT HERE -->(.+)<!-- END CUSTOM TEXT -->',prevtext)
#                    if prevtext: prevtext = string.replace(prevtext[0],'!!!N!!!','\n')
#                    else: prevtext = '\n\n'
#                else: prevtext = '\n\n'
                pagetxt += prevtext
                pagetxt += '<!-- END CUSTOM TEXT -->\n\n'
                pagetxt += '<hr>'
                pagetxt += '<table width="100%" border=0>\n'
                pagetxt += '<tr id="linktab"><th>Module</th><th>Option</th><th>Description</th><th>Default</th></tr>\n'
                modrows = []
                for entry in odb.indexEntries('Argument',arg):
                    omod = string.split(entry['File'],'/')[-1][:-3]
                    modtxt = '<tr>\n  <td id="linktab2" align=left>[%s]{module:%s}</td>\n' % (omod,omod)
                    modtxt += '  <td><code>[%s=%s]{docs:%s}</code></td>\n' % (arg,entry['ArgType'],entry['page'])
                    modtxt += '  <td>%s</td>\n' % entry['Description']
                    modtxt += '  <td><code>%s</code></td>\n</tr>\n' % entry['Default']
                    if modtxt not in modrows: modrows.append(modtxt)
                modrows.sort()
                pagetxt += string.join(modrows,'')
                pagetxt += '</tr></table>\n'
                open(apage,'w').write(pagetxt)
            ### [2] Option Type Pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'%sdocs/' % pagedir)
            self.debug(rje.sortKeys(odb.index('page')))
            for arg in odb.index('page'):
                if arg == '?': continue
                if arg == 'T/F': apage = '%sdocs/boolean.page' % pagedir
                else: apage = '%sdocs/%s.page' % (pagedir,arg)
                self.debug(apage)
                pagetxt = '<h2>Option Type: %s</h2>\n' % (arg)
                pagetxt += '<p><b>Argument:</b>'
                for cmd in odb.indexDataList('page',arg,'Argument'):
                    pagetxt += ' <code>[%s]{cmd:%s}</code>' % (cmd,cmd)
                pagetxt += '</p>\n'
                pagetxt += '<p><b>Modules:</b>'
                for ofile in odb.indexDataList('page',arg,'File'):
                    omod = string.split(ofile,'/')[-1][:-3]
                    pagetxt += ' <code>[%s]{module:%s}</code>' % (omod,omod)
                pagetxt += '</p>\n\n'
                #pagetxt += '<!-- INSERT CUSTOM TEXT HERE -->'
                prevpage = string.replace(apage,pagedir,prevdir)
                if os.path.exists(prevpage):
                    prevtext = open(prevpage,'r').read()
                    prevtext = prevtext[prevtext.find('<!-- INSERT CUSTOM TEXT HERE -->'):prevtext.find('<!-- END CUSTOM TEXT -->')]
                    if prevtext: self.bugPrint(prevtext)
                    #prevtext = rje.matchExp('<!-- INSERT CUSTOM TEXT HERE -->(.+)<!-- END CUSTOM TEXT -->',open(prevpage,'r').read())
                    #if prevtext: prevtext = prevtext[0]
                    else: prevtext = '<!-- INSERT CUSTOM TEXT HERE -->\n\n'
                else: prevtext = '<!-- INSERT CUSTOM TEXT HERE -->\n\n'
                pagetxt += prevtext
                pagetxt += '<!-- END CUSTOM TEXT -->\n\n'
                pagetxt += '<table width="100%" border=0>\n'
                pagetxt += '<tr id="linktab"><th>Module</th><th>Option</th><th>Description</th><th>Default</th></tr>\n'
                modrows = []
                for entry in odb.indexEntries('page',arg):
                    omod = string.split(entry['File'],'/')[-1][:-3]
                    otype = entry['Argument']
                    modtxt = '<tr>\n  <td id="linktab2" align=left>[%s]{module:%s}</td>\n' % (omod,omod)
                    modtxt += '  <td><code>[%s=%s]{cmd:%s}</code></td>\n' % (otype,entry['ArgType'],otype)
                    modtxt += '  <td>%s</td>\n' % entry['Description']
                    modtxt += '  <td><code>%s</code></td>\n</tr>\n' % entry['Default']
                    modrows.append(modtxt)
                modrows.sort()
                pagetxt += string.join(modrows,'')
                pagetxt += '</tr></table>\n'
                open(apage,'w').write(pagetxt)
            ### [3] End of Page Generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except: return self.errorLog('Error in PyDoc.makePages()')
#########################################################################################################################
### End of SECTION II: PyDoc Class                                                                                      #
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
    try: PyDoc(mainlog,cmd_list).run()

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
