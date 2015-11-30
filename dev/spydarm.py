#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SPyDARM
Description:  SLiMSuite Python Documentation And ReadMe
Version:      0.0.0
Last Edit:    01/01/15
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    ### ~ SPyDARM Release Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X  : Prefix for output files [slimsuite]
    prevbase=X  : Prefix of output files for previous release (to update) [slimsuite]
    backbase=X  : Prefix of backed up output files for previous release [slimsuite.#DATE]
    ### ~~~ Python Module Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    pylist=LIST     : List of python modules to upload. Can have * wildcards. Will add '.py' if missing. ['*']
    pypath=PATH     : Path to python modules. Will also look in listed sourcedir subfolders. [../]
    sourcedir=LIST  : List of subdirectories in which to look for modules [tools,extras,libraries,legacy]
    addimports=T/F  : Whether to add identified imported modules to python module list [True]
    ### ~ SPyDARM Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    output=LIST : List of outputs to generate [release,history,updates,readme,dependencies]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
sys.path.append(os.path.join(slimsuitepath,'extras/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_html, rje_obj, rje_pydocs
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [?] : Add REST outputs to restSetup() and restOutputOrder()
    # [?] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SPyDARM', '0.0.0', 'May 2015', '2015')
    description = 'Generic RJE Module'
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
### SECTION II: SPyDARM Class                                                                                           #
#########################################################################################################################
class SPyDARM(rje_obj.RJE_Object):
    '''
    SPyDARM Class. Author: Rich Edwards (2015).

    Str:str
    - BackBase = Prefix of backed up output files for previous release [slimsuite.#DATE]
    - PrevBase = Prefix of output files for previous release (to update) [slimsuite]

    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - Output = List of outputs to generate [release,readme,dependancies]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - PyDoc = rje_pydocs.PyDoc object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BackBase','PrevBase']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = ['Output']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.baseFile('slimsuite')
        self.setStr({'BackBase':'slimsuite.%s' % rje.dateTime(dateonly=True),'PrevBase':'slimsuite'})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        self.list['Output'] = string.split('release,history,updates,readme,dependencies',',')
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
                self._cmdReadList(cmd,'str',['BackBase','PrevBase'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'lclist',['Output'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with %s.cmd:%s' % (self.me(),cmd))
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dependencies(): return False
            if not self.release(): return False
            # Generate a readme text file
            return True
        except:
            self.errorLog(self.zen())
            return False
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pydoc = self.obj['PyDoc'] = rje_pydocs.PyDoc(self.log,self.cmd_list+['name=%s' % self.basefile(strip_path=True)])
            if not pydoc.setup(): return False
            pydoc.parseModules()
            db = self.obj['DB'] = pydoc.db()
            db.basefile(self.basefile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        There is currently no specific help available on REST output for this program. Run with &rest=help for general
        options. Run with &rest=full to get full server output. Individual outputs can be identified/parsed:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ...

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Main SPyDARM Methods                                                                                    #
#########################################################################################################################
    def release(self):  ### Generate the release information tables.
        '''Generate the release information tables.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            hdb = None      # History table: dir, module, version, update, release
            basefile = self.basefile()
            prevbase = self.getStr('PrevBase')
            backbase = self.getStr('BackBase')
            if backbase == prevbase: raise ValueError('BackBase cannot match PrevBase ("%s")' % prevbase)
            if backbase == basefile: raise ValueError('BackBase cannot match BaseFile ("%s")' % basefile)
            ## ~ [1a] Load & Backup previous release ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for sfile in ['release.tdt','history.tdt','readme.txt','updates.html']:
                pfile = '%s.%s' % (prevbase,sfile)
                bfile = '%s.%s' % (backbase,sfile)
                if os.path.exists(pfile):
                    if os.path.exists(bfile): rje.backup(self,bfile)
                    open(bfile,'w').write(open(pfile,'r').read())
                    self.printLog('#BACK','%s => %s' % (pfile,bfile))
                    if sfile == 'history.tdt': hdb = db.addTable(filename=pfile,mainkeys=['Dir','Module','Version'],name='history',expect=True)
            ### ~ [2] Generate slimsuite.release.tdt based on *.Module.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = db.copyTable(self.db('Module'),'release')
            rdb.renameField('SourceDir','Dir')
            rdb.newKey(['Dir','Module'])
            rdb.dropFields(['File','Classes','Methods'])
            if 'release' in self.list['Output']: rdb.saveToFile(backup=False)
            ### ~ [3] Generate slimsuite.history.tdt, parsed from docstrings, based on pydoc.distribute() ~~~~~~~~~~~ ###
            if not hdb:
                hdb = db.addEmptyTable('history',['Dir','Module','Version','Update','Release'],['Dir','Module','Version'])
            self.makeHistory()
            # Generate slimsuite.readme.txt based on pydoc.saveDocs()
            if 'readme' in self.list['Output']: self.saveReadMe('%s.readme.txt' % basefile)
            return True
        except: self.errorLog('%s.release error' % self.prog()); return False
#########################################################################################################################
    def readMeHeader(self):
        rtxt = '# ReadMe documentation for SLiMSuite software\n\n'
        rtxt += '> Copyright (C) 2015 Richard J. Edwards\n\n'
        rtxt += '> Distribution compiled: %s\n\n' % time.asctime(time.localtime(time.time()))
        rtxt += '> Please visit http://www.slimsuite.unsw.edu.au for additional documentation\n'
        rtxt += '\n\n# GNU License\n\n'
        rtxt += '> This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License\n'
        rtxt += '> as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.\n'
        rtxt += '> \n'
        rtxt += '> This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied\n'
        rtxt += '> warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.\n'
        rtxt += '> \n'
        rtxt += '> You should have received a copy of the GNU General Public License along with this program; if not, write to\n'
        rtxt += '> the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.\n'
        rtxt += '> \n'
        rtxt += '> Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.\n'
        rtxt += '> \n'
        rtxt += '> To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py.\n'
        rtxt += '\n\n# Module Docstrings\n\n'
        return rtxt
#########################################################################################################################
    def saveReadMe(self,filename='pydocs.txt',append=False):      ### Prints docs for modules to file
        '''
        Prints docs for modules to file.
        >> filename:str = output file name
        >> append:boolean
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pydoc = self.obj['PyDoc']
            if append:
                self.printLog('#DOC','Appending docstrings to %s' % filename)
                PYDOC = open(filename,'a')
            else:
                rje.mkDir(self,filename)
                self.printLog('#DOC','Writing docstrings to %s' % filename)
                PYDOC = open(filename,'w')
                PYDOC.write(self.readMeHeader())
            db = self.db('Module')
            dx = 0
            ### ~ [2] Output Docstrings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for sourcedir in pydoc.list['SourceDir']:
                PYDOC.write('-%s:\n\n' % sourcedir)
                for pyfile in db.dataKeys():
                    entry = db.data(pyfile)
                    module = entry['Module']
                    if not pyfile.find(sourcedir) >= 0 or not os.path.exists('%s%s%s.py' % (pydoc.getStr('PyPath'),rje.makePath(sourcedir),module)): continue
                    ## ~ [2a] ~ Module docstring ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    mtxt = '### ~~~ Module %s ~ [%s] ~~~ ###' % (module,pyfile)
                    while len(mtxt) < 122: mtxt = mtxt[:5] + '~' + mtxt[5:-5] + '~' + mtxt[-5:]
                    PYDOC.write('%s\n\n%s\n' % (mtxt,entry['DocString'])); dx += 1
                PYDOC.write('\n\n\n')
            PYDOC.close()
            self.printLog('#DOC','Output to %s complete: %s modules.' % (filename,rje.iStr(dx)))
        except: self.errorLog('Error in %s.saveDocs()' % self.prog())
#########################################################################################################################
    def makeHistory(self):  ### Extracts history information from docstrings.
        '''Extracts history information from docstrings.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdb = self.db('history')
            mdb = self.db('Module')
            pydoc = self.obj['PyDoc']
            uhtml = []  # Update HTML text
            udir = ''

            ### ~ [1] Work through python modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mkey in mdb.dataKeys():
                entry = mdb.data(mkey)
                pyfile = entry['File']
                mod = entry['Module']
                prev = '-'
                lastv = ''
                ## ~ [1a] Parse out history text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                history = []    # History lines
                updates = []    # Version number entries
                for mentry in self.db('Method').indexEntries('File',pyfile):
                    if mentry['Method'] == 'history': history = string.split(mentry['DocString'],'\n'); break
                for dline in history:
                    if rje.matchExp('# (\d+\.\d+\.?\d*)\s?-\s(\S.+)$',dline):
                        (v,text) = rje.matchExp('# (\d\.\d+\.?\d*)\s?-\s(\S.+)$',dline)
                        if v == lastv:  # Continuation of update
                            if prev == lastv: continue  # In previous release
                            updates[-1]['Update'] += ' %s' % text
                            continue
                        lastv = v
                        ventry = {'Dir':entry['SourceDir'],'Module':mod,'Version':v,'Update':text,'Release':rje.dateTime(dateonly=True)}
                        vkey = hdb.makeKey(ventry)
                        if hdb.data(vkey): prev = v
                        else: updates.append(ventry)
                ## ~ [1b] Assess/report updates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if prev != lastv:
                    if entry['SourceDir'] != udir: uhtml.append('<h2>Updates in %s/:</h2>\n' % entry['SourceDir']); udir = entry['SourceDir']
                if prev == '-' and lastv:
                    self.printLog('#VNUM','%s: Creation -> Version %s' % (mod,entry['Version']))
                    uhtml.append('<p><b>&bull; %s:</b> <i>Created/Renamed/moved.</i>' % mod)
                elif prev not in ['-',lastv]:
                    self.printLog('#VNUM','%s: Version %s -> Version %s' % (mod,prev,entry['Version']))
                    uhtml.append('<p><b>&bull; %s:</b> <i>Updated from Version %s.</i>' % (mod,prev))
                for ventry in updates:
                    self.printLog('#V%s' % ventry['Version'],ventry['Update'])
                    uhtml.append('<br>&rarr; Version %s: %s' % (ventry['Version'],ventry['Update']))
                    hdb.addEntry(ventry)
                if uhtml and uhtml[-1] != '</p>': uhtml.append('</p>')
                if lastv != entry['Version']: self.warnLog('Module %s Version %s but history() ends at %s' % (mod,entry['Version'],lastv))
                self.deBug('>>>')
            if 'history' in self.list['Output']: hdb.saveToFile(backup=False)

            ### ~ [2] Make updates.html file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'updates' in self.list['Output']:
                htmlfile = '%s.updates.html' % self.basefile()
                title = 'SLiMSuite updates'
                stylesheets = []
                for css in pydoc.list['StyleSheets']: stylesheets.append(pydoc.getStr('StylePath')+css)
                htmlhead = rje_html.htmlHead(title,stylesheets,tabber=True,frontpage=False,nobots=False,keywords=pydoc.list['Keywords'],javascript=pydoc.getStr('StylePath'))
                htmlbody = string.join(['<h1>SLiMSuite updates</h1>'] + uhtml,'\n')
                htmltail = rje_html.htmlTail('%s %s' % (pydoc.getStr('Author'),string.split(time.asctime(time.localtime(time.time())))[-1]))
                open(htmlfile,'w').write(htmlhead+htmlbody+htmltail)
                self.printLog('#HTML','HTML update summary output to %s' % (htmlfile))
        except: self.errorLog('Error in %s.makeHistory()' % self.prog())
#########################################################################################################################
    def dependencies(self): ### Check dependencies of loaded modules.
        '''Check dependencies of loaded modules.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb = self.db('Module')     # Use Module, Imports, Imported_By
            for mfield in ['Module', 'Imports', 'Imported_By']:
                mdb.index(mfield,splitchar=', ')
                mdb.indexReport(mfield)
            ### ~ [1] Add/Check dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.addField('Dependencies')
            for mkey in mdb.dataKeys():
                mentry = mdb.data(mkey)
                mod = mentry['Module']
                mdir = mentry['SourceDir']
                # Cycle through imports and add their imports as dependencies
                imports = string.split(mentry['Imports'],', ')
                mentry['Dependencies'] = []
                ilist = imports[0:]
                while ilist:
                    imod = ilist.pop(0)
                    if imod == 'None': continue
                    dlist = mdb.indexDataList('Module',imod,'Imports',sortunique=False)
                    dlist = string.join(dlist,', ')
                    dlist = string.split(dlist,', ')
                    for dmod in dlist:
                        if dmod in ['None',''] + imports + mentry['Dependencies']: continue
                        mentry['Dependencies'].append(dmod)
                        ilist.append(dmod)
                mentry['Dependencies'].sort()
                if not mentry['Dependencies']: mentry['Dependencies'] = ['None']
                mentry['Dependencies'] = string.join(mentry['Dependencies'],', ')
                # Cycle through imported_by modules and check for orphans (e.g. Nones and libraries/ only)
                if mentry['SourceDir'] != 'libraries': continue
                self.debug(mentry)
                if mentry['Imported_By'] == 'None':
                    self.warnLog('%s/%s is not imported by anything!' % (mdir,mod))
                    continue
                idir = []   # List of sourcedir for importing modules
                imported = string.split(mentry['Imported_By'],', ')
                ilist = imported[0:]
                while ilist:
                    imod = ilist.pop(0)
                    if imod == 'None': continue
                    for sdir in mdb.indexDataList('Module',imod,'SourceDir'):
                        if sdir not in idir: idir.append(sdir)
                    dlist = mdb.indexDataList('Module',imod,'Imported_By',sortunique=False)
                    dlist = string.join(dlist,', ')
                    dlist = string.split(dlist,', ')
                    for dmod in dlist:
                        if dmod in ['None'] + imported: continue
                        imported.append(dmod)
                        ilist.append(dmod)
                required = []
                for sdir in idir[0:]:
                    if sdir not in ['legacy','libraries','servers','dev']:
                        required.append(sdir)
                        idir.remove(sdir)
                if required: continue
                if idir: self.warnLog('%s/%s is a dependency for %s only.' % (mdir,mod,string.join(idir,'/')))
            ### ~ [2] Output dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'dependencies' in self.list['Output']:
                all_dependencies = []
                dfile = '%s.dependencies.txt' % self.basefile()
                DFILE = open(dfile,'w')
                DFILE.write('# The following modules have dependencies:\n')
                for mkey in mdb.dataKeys():
                    mentry = mdb.data(mkey)
                    if mentry['SourceDir'] == 'libraries': continue
                    DFILE.write('%s/%s: ' % (mentry['SourceDir'],mentry['Module']))
                    if mentry['Imports'] == mentry['Dependencies'] == 'None': DFILE.write('None\n')
                    else:
                        if mentry['Imports'] != 'None':
                            DFILE.write(mentry['Imports'])
                            if mentry['Dependencies'] != 'None': DFILE.write(', ')
                            all_dependencies += string.split(mentry['Imports'],', ')
                        if mentry['Dependencies'] != 'None':
                            DFILE.write(mentry['Dependencies'])
                            all_dependencies += string.split(mentry['Dependencies'],', ')
                        DFILE.write('\n')
                DFILE.write('\n')
                DFILE.write('# Full module list required, including dependencies:\n')
                for mkey in mdb.dataKeys():
                    mentry = mdb.data(mkey)
                    if mentry['SourceDir'] == 'libraries' and mentry['Module'] not in all_dependencies: continue
                    DFILE.write('%s/%s\n' % (mentry['SourceDir'],mentry['Module']))
                DFILE.write('\n')
                DFILE.write('# Orphan modules, not required:\n')
                for mkey in mdb.dataKeys():
                    mentry = mdb.data(mkey)
                    if mentry['SourceDir'] == 'libraries' and mentry['Module'] not in all_dependencies:
                        DFILE.write('%s/%s\n' % (mentry['SourceDir'],mentry['Module']))
                DFILE.close()
                self.printLog('#OUT','Module dependencies output to %s.' % dfile)
                return True
        except: self.errorLog('Error in %s.makeHistory()' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: SPyDARM Class                                                                                    #
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
    try: SPyDARM(mainlog,cmd_list).run()

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
