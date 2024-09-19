#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_glossary
Description:  RJE Glossary HTML Maker
Version:      1.4
Last Edit:    17/12/12
Copyright (C) 2012  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module is to convert a plain text file of glossary terms and definitions into a webpage complete
    with optional hyperlinks between terms.

Commandline:
    ### Input Options ###
    infile=FILE   : Input file of Term:Definition [glossary.tdt]
    termsplit=X   : Text to use for splitting term from description (if file extension not recognised) [tab]
    name=X        : Title for HTML output ['Glossary']
    terms=LIST    : List of terms to extract from glossary (all by default) []

    ### Term Linking Options ###
    aname=T/F     : Whether to hyperlink terms using a name refs [True]
    plurals=T/F   : Whether to map plurals onto singular terms [True]
    href=T/F      : Whether to added external hyperlinks for <url> and <url>[text] in descriptions [True]
        
    ### Output Options ###
    htmlstyle=X   : Output HTML style for splits (tab/table/h3) [h3]
    splits=X      : Number of sets to divide terms into [6]
    outfile=FILE  : Output file name (input.html by default) []
    copyright=X   : Copyright statement for page ['RJ Edwards 2012']
    keeporder=T/F : Keep output order the same as input order (unless terms=LIST given). Uses termsplit=X. [False]

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
import rje, rje_html, rje_obj, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Working version, including text setup for webserver.
    # 1.1 - Added href=T option to add external hyperlinks for <url> and <url>[text] in descriptions [True]
    # 1.2 - Added recognition of _italics_ markup.
    # 1.3 - Fixed minor italicising bug.
    # 1.4 - Added keeporder=T/F to maintain input order (e.g. for MapTime).
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Deal with plurals.
    # [ ] : Add better treatment of aliases and abbreviations. (Can parse from "See X")
    # [Y] : Add recognition of italics.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_GLOSSARY', '1.4', 'December 2012', '2012')
    description = 'RJE Glossary HTML Maker'
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
### SECTION II: Glossary Class                                                                                          #
#########################################################################################################################
class Glossary(rje_obj.RJE_Object):     
    '''
    Glossary Class. Author: Rich Edwards (2012).

    Str:str
    - InFile = Input file of Term:Definition [glossary.tdt]
    - HTMLStyle = Output HTML style for splits (tab/table/h3) [h3]
    - Name = Title for HTML output 
    - OutFile = Output file name (*.html)
    - TermSplit = Text to use for splitting term from description (if file extension not recognised) [tab]
    
    Bool:boolean
    - AName = Whether to hyperlink terms using a name refs [True]
    - HRef = Whether to added external hyperlinks for <url> and <url>[text] in descriptions [True]
    - KeepOrder = Keep output order the same as input order (unless terms=LIST given) [False]
    - Plurals = Whether to map plurals onto singular terms [True]

    Int:integer
    - Splits = Number of sets to divide terms into [6]

    Num:float
    
    List:list
    - Terms = List of input terms in input case

    Dict:dictionary
    - Glossary = Nested dictionary of lower case term words and '=' definition    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['InFile','HTMLStyle','OutFile','Name','TermSplit']
        self.boollist = ['AName','HRef','KeepOrder','Plurals']
        self.intlist = ['Splits']
        self.numlist = []
        self.listlist = ['Terms']
        self.dictlist = ['Glossary']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'InFile':'glossary.tdt','HTMLStyle':'h3','Name':'Glossary','TermSplit':'tab'})
        self.setBool({'AName':True,'HRef':True,'Plurals':True})
        self.setInt({'Splits':6})
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
                self._cmdReadList(cmd,'str',['HTMLStyle','Name','TermSplit'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['InFile','OutFile'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['AName','HRef','KeepOrder','Plurals'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['Splits'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Terms'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,gtext=''):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup(gtext)
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = self.glossaryHTML()
            hobj = self.obj['HTML']
            date = rje.split(time.asctime(time.localtime(time.time())))
            date = '%s %s %s' % (date[2],date[1],date[-1])
            hobj.info['Copyright'] += '. Generated by rje_glossary.py'
            title = '%s' % self.getStr('Name')
            tabber = self.getStr('HTMLStyle').lower() == 'tab'
            frontpage = True
            html = '%s\n\n%s\n\n%s' % (hobj.htmlHead(title,tabber,frontpage),html,hobj.htmlTail(tabber))
            if not gtext:   # Replace with CGI option
                rje.backup(self,self.getStr('OutFile'),appendable=False)
                open(self.getStr('OutFile'),'w').write(html)
                self.printLog('#HTML','%s HTML output to %s' % (title,self.getStr('OutFile')))
            return html
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self,gtext=''):    ### Main class setup method. gtext will over-ride input file.
        '''Main class setup method. gtext will over-ride input file.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['HTML'] = rje_html.HTML(self.log,self.cmd_list)
            ## ~ [1a] File names etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.basefile().lower() in ['','none']: self.basefile(rje.baseFile(self.getStr('InFile')))
            if self.getStr('OutFile').lower() in ['','none']: self.str['OutFile'] = '%s.html' % self.basefile()
            ## ~ [1b] Read in Glossary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            interms = []
            if gtext:
                delimit = self.getStr('TermSplit')
                if delimit.lower() == 'tab': delimit = '\t'
                if delimit.lower() == 'space': delimit = ' '
                if delimit.lower() == 'comma': delimit = ','
                if delimit.lower() == 'period (.)': delimit = '.'
                if delimit.lower() == 'colon': delimit = ':'
                glossary = {}
                for line in rje.split(gtext,'\n'):
                    splitline = rje.split(line,delimit)
                    if delimit == '.' and (splitline[-1] in ['',' ']): splitline = splitline[:-1]
                    if not splitline: continue
                    (term,definition) = (splitline[0],rje.join(splitline[1:],delimit))
                    if term == 'Term' and not glossary: continue
                    if term:
                        glossary[term] = {'Definition':definition}
                        interms.append(term)
            else: 
                try:
                    if not self.getBool('KeepOrder') and open(self.getStr('InFile'),'r').readline()[:4] == 'Term': 
                        glossary = rje.dataDict(self,self.getStr('InFile'),mainkeys=['Term'],datakeys=['Term','Definition'])
                    else: return self.setup(open(self.getStr('InFile'),'r').read())
                except: 
                    self.errorLog('Problem reading input as dataDict(). Will try as text.')
                    return self.setup(open(self.getStr('InFile'),'r').read())
            if self.list['Terms']:
                for term in glossary:
                    if term not in self.list['Terms']: glossary.pop(term)
            elif self.getBool('KeepOrder'): self.list['Terms'] = interms
            else: self.list['Terms'] = rje.sortKeys(glossary)
            for term in glossary: glossary[term] = glossary[term]['Definition']
            ### ~ [2] Create Full Glossary Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nested = {}
            for term in glossary:
                tdict = nested
                for word in rje.split(term.lower()):
                    if word not in tdict: tdict[word] = {}
                    tdict = tdict[word]
                tdict['='] = glossary[term]
            self.dict['Glossary'] = nested
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def glossaryHTML(self): ### Main glossary generator. Establishes links and outputs HTML.
        '''Main glossary generator. Establishes links and returns HTML.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Add HTML links to definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('AName'): self.addLinks(self.dict['Glossary'])
            ## ~ [1b] Divide terms into sets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            termsets = []
            self.setInt({'Splits': min(self.getInt('Splits'),len(self.list['Terms']))})
            if self.getInt('Splits') > 1 and self.getStr('HTMLStyle') == 'tab': 
                termx = int(len(self.list['Terms'])/float(self.getInt('Splits')))
                while termx * self.getInt('Splits') < len(self.list['Terms']): termx += 1
                i = 0
                while i < len(self.list['Terms']):
                    termsets.append(self.list['Terms'][i:i+termx])
                    i += termx
                ## ~ [1c] Calculate number of letters needed for a-b set identifiers ~~~~~~~~~~~~~~~~~~ ##
                x = 0; unique = False
                while not unique:
                    x += 1; unique = True
                    setends = []
                    for termset in termsets:
                        #if termset[0][:x] in setends: unique = False
                        #else:
                        setends.append(termset[0][:x])
                        if termset[-1][:x] in setends: unique = False
                        else: setends.append(termset[-1][:x])
                x = 3   #!# Over-riding for the moment. Add as option?
                ## ~ [1d] Generate termset references ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                termref = []
                for termset in termsets: termref.append('%s-%s' % (termset[0][:x],termset[-1][:x]))
                self.deBug(termref)
            if not self.getBool('KeepOrder') and self.getStr('HTMLStyle') != 'tab':
                quicklink = []
                for a in string.ascii_uppercase:
                    for term in self.list['Terms']:
                        if term[0].upper() == a: quicklink.append(a); break
                for a in quicklink:
                    termsets.append([])
                    for term in self.list['Terms']:
                        if term[0].upper() == a: termsets[-1].append(term)
            ### ~ [2] Generate HTML Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hobj = self.obj['HTML']
            html = []   # List of HMTL lines to join and output
            ## ~ [2a] HTML Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            title = '%s' % self.getStr('Name')
            ## ~ [2b] Glossary Title and shortcuts to sets (unless tabbed) ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html.append('<h1>%s</h1>' % title)
            #if termsets and self.getStr('HTMLStyle') != 'tab':
            #    html.append('<ul>')
            #    for ref in termref: html.append('<li><a href="#%s">%s</a>' % (ref.lower(),ref))
            #    html.append('</ul>')
            if not self.getBool('KeepOrder') and self.getStr('HTMLStyle') != 'tab':
                html.append('<h3><a href="#%s">%s</a>' % (quicklink[0],quicklink[0]))
                for a in quicklink[1:]: html.append('~ <a href="#%s">%s</a>' % (a,a))
                html[-1] += '</h3>'
            if not termsets:
                termsets = [self.list['Terms']]
                #termref = [termref.append('%s-%s' % (termsets[0][0][:x],termsets[0][-1][:x]))]
                termref = [title]
                quicklink = ['Terms']
                self.deBug(termref)
            ## ~ [2c] Glossary sets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            htmlsets = []
            for termset in termsets:
                sethtml = []    #'<ul>']
                for term in termset:
                    sethtml.append('<!-- %s -->' % term)
                    akey = []; glossary = self.dict['Glossary']
                    for word in rje.split(term.lower()): akey.append(word); glossary = glossary[word]
                    #sethtml.append('<a name="%s"><li></a><p><scaps><b>%s.</b></scaps>' % (rje.join(akey,'_'),term))
                    sethtml.append('<a name="%s"><p></a><scaps><b>%s.</b></scaps>' % (rje.join(akey,'_'),term))
                    dkey = {True:'+',False:'='}['+' in glossary]
                    self.deBug(dkey)
                    if not glossary[dkey]: glossary[dkey] = '<i>No definition</i>'
                    if glossary[dkey][-1] != '.': glossary[dkey] += '.'
                    sethtml.append('%s</p>' % glossary[dkey])
                #sethtml.append('</ul>')
                htmlsets.append(rje.join(sethtml,'\n'))
            ## ~ [2d] HTML Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('HTMLStyle') == 'tab':
                tablist = []
                for hset in htmlsets:
                    termset = termsets[htmlsets.index(hset)]
                    if len(termsets) > 1: tref = termref[htmlsets.index(hset)]
                    else: tref = 'Glossary'
                    tablist.append((tref,hset,"Terms '%s' to '%s'" % (termset[0],termset[-1])))
                html.append(rje_html.tabberHTML('GlossaryTabs',tablist))
            elif self.getStr('HTMLStyle') == 'table':
                html += ['<table border=%d align="top">' % hobj.stat['Border'],'<tr><th>#</th><th>Definitions</th></tr>']
                for i in range(len(htmlsets)):
                    html.append('<tr><td valign=TOP align=CENTER><a name="%s"></a><h3>%s</h3></td>' % (quicklink[i],quicklink[i]))
                    html.append('<td valign=TOP>%s</td></tr>' % htmlsets[i])
                html.append('</table>')
            else:
                for i in range(len(htmlsets)):
                    html.append('<a name="%s"><%s></a>%s</%s>' % (quicklink[i],self.getStr('HTMLStyle'),quicklink[i],self.getStr('HTMLStyle')))
                    html.append(htmlsets[i])
            return rje.join(html,'\n')
        except: return self.errorLog('%s.glossary error' % self)
#########################################################################################################################
    def addLinks(self,nested): ### Adds href aname links to definitions.
        '''Adds href aname links to definitions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            endstrip = [')','.',',',':',';','!']
            if self.getBool('Plurals'): endstrip.append('s')
            for term in rje.sortKeys(nested):
                if term == '=':
                    linkdef = []
                    rawdef = rje.split(rje.replace(nested['='],'(','( '))
                    while rawdef:
                        glossary = self.dict['Glossary']
                        if self.getBool('HRef') and rje.matchExp('<(\S+)>',rawdef[0]):
                            safetynet = rawdef[0:]
                            url = rje.matchExp('<(\S+)>',rawdef[0])[0]
                            if rje.matchExp('<(\S+)>\[(\S+)',rawdef[0]): rawdef[0] = '[%s' % rje.matchExp('<(\S+)>\[(\S+)',rawdef[0])[1]
                            elif rje.matchExp('<(\S+)>(\S+)',rawdef[0]): rawdef[0] = '[%s]%s' % (url,rje.matchExp('<(\S+)>(\S+)',rawdef[0])[1])
                            else: rawdef[0] = '[%s]' % url
                            try:
                                while ']' not in rawdef[0]: rawdef[0] = '%s %s' % (rawdef[0],rawdef.pop(1))
                                (linktext,linkextra) = rje.matchExp('\[(.+)\](\S*)',rawdef.pop(0))
                                if url[:3] not in ['htt','ftp']: url = 'http://%s' % url
                                linkdef.append('<a href="%s">%s</a>%s' % (url,linktext,linkextra))
                                continue
                            except:
                                self.errorLog('Problem parsing URL from "%s"' % nested['='])
                                rawdef = safetynet
                        if rawdef[0].lower() not in glossary:
                            if rawdef[0].lower()[:-1] not in glossary or rawdef[0].lower()[-1] not in endstrip:
                                linkdef.append(rawdef.pop(0)); continue
                        akey = []; alink = []
                        while rawdef and (rawdef[0].lower() in glossary or rawdef[0].lower()[:-1] in glossary):
                            if rawdef[0].lower() in glossary and '=' in glossary[rawdef[0].lower()]: rterm = rawdef[0].lower()
                            elif len(rawdef) > 1 and rawdef[0].lower() in glossary and (rawdef[1].lower() in glossary[rawdef[0].lower()] or rawdef[1].lower()[:-1] in glossary[rawdef[0].lower()]): rterm = rawdef[0].lower()
                            elif rawdef[0].lower()[-1] in endstrip and rawdef[0].lower()[:-1] in glossary: rterm = rawdef[0].lower()[:-1]
                            elif rawdef[0].lower() in glossary: rterm = rawdef[0].lower()
                            else: break
                            glossary = glossary[rterm]
                            akey.append(rterm)
                            alink.append(rawdef.pop(0))
                        akey = rje.join(akey,'_')
                        if '=' in glossary:
                            alink = rje.join(alink)
                            if nested == glossary: linkdef.append(alink)
                            elif self.getStr('HTMLStyle') != 'tab':
                                if alink[-1] in endstrip and alink[-1] != 's': linkdef.append('<a href="#%s">%s</a>%s' % (akey,alink[:-1],alink[-1]))
                                else: linkdef.append('<a href="#%s">%s</a>' % (akey,alink))
                            else:
                                if alink[-1] in endstrip and alink[-1] != 's': linkdef.append('<scaps>%s</scaps>%s' % (alink[:-1],alink[-1]))
                                else: linkdef.append('<scaps>%s</scaps>' % (alink))
                        else:
                            linkdef.append(alink[0])
                            rawdef = alink[1:] + rawdef
                    nested['+'] = rje.replace(rje.join(linkdef),'( ','(')
                    while rje.matchExp(' _([^_]+)_',nested['+']):
                        italics = rje.matchExp(' _([^_]+)_',nested['+'])[0]
                        nested['+'] = rje.replace(nested['+'],' _%s_' % italics,' <i>%s</i>' % italics)
                    #self.deBug(nested)
                elif term != '+': self.addLinks(nested[term])
        except: self.errorLog('%s.addLinks error' % self)
#########################################################################################################################
### End of SECTION II: Glossary Class                                                                                   #
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
    try: Glossary(mainlog,cmd_list).run()

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
