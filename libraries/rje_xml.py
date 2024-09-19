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
Module:       rje_xml.py
Description:  XML (Text) Parsing Module
Version:      0.2
Last Edit:    06/08/13
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the XML class for parsing XML files. These are parsed into a generic set of nested dictionaries
    and XML objects stored in the root XML object. The main attributes of interest in the XML object for retrieving data
    are:
    * info['Name'] = element name
    * info['Content'] = text content of element (if any)
    * dict['Attributes'] = dictionary of element attributes (if any)
    * list['XML'] = list of XML objects containing nested elements (if any)

    For example, the following short XML:    
    < ?xml version="1.0" encoding="ISO-8859-1"? >
    < database name="EnsEMBL" ftproot="ftp://ftp.ensembl.org/pub/" outdir="EnsEMBL" >
          < file path="current_aedes_aegypti/data/fasta/pep/*.gz" >Yellow Fever Mosquito< /file >
          < file path="current_anopheles_gambiae/data/fasta/pep/*.gz" >Malaria Mosquito< /file >
    < /database >

    The following XML objects would be created:
    * XML root object:
        XML.info['Name'] = filename
        XML.list['XML'] = [XML1]
    * XML1: 
        XML1.info['Name'] = 'database'
        XML1.info['Content'] = ''
        XML1.dict['Atrributes'] = {'name':"EnsEMBL",'ftproot':"ftp://ftp.ensembl.org/pub/",'outdir':"EnsEMBL"}
        XML1.list['XML'] = [XML2,XML3]
    * XML2:
        XML2.info['Name'] = 'file'
        XML2.info['Content'] = 'Yellow Fever Mosquito'
        XML2.dict['Atrributes'] = {'path':"current_aedes_aegypti/data/fasta/pep/*.gz"}
        XML2.list['XML'] = []
    * XML3:
        XML3.info['Name'] = 'file'
        XML3.info['Content'] = 'Malaria Mosquito'
        XML3.dict['Atrributes'] = {'path':"current_anopheles_gambiae/data/fasta/pep/*.gz"}
        XML3.list['XML'] = []
    
    The top level XML list (XML.list['XML']) is returned by the parseXML() method of the class, which populates all the
    objects.

Commandline:
    - parse=FILE        : Source file for reading XML file [None]
    - attributes=LIST   : List of Attributes to exclusively extract []
    - elements=LIST     : List of Elements to exclusively extract []

Uses general modules: copy, os, string, sys, time, xml.sax
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
try: import urllib2 as urllib2
except: import urllib.request as urllib2
import xml.sax
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added xml.sax functions.
    # 0.2 - Added parsing from URL.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add options for only extracting partial data, i.e. certain element types and/or attributes.
    # [X] : Convert nested objects into dictionary - cannot because of multiple elements of same type!
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_XML', '0.2', 'January 2013', '2006')
    description = 'XML Parsing Module'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
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
### SECTION II: XML Class                                                                                               #
#########################################################################################################################
class XML(rje.RJE_Object):     
    '''
    XML Class. Author: Rich Edwards (2005).

    Info:str
    - Name = element name
    - Content = text content of element (if any)
    - Source = File name for top level Element
    
    Opt:boolean

    Stat:numeric

    List:list
    - Attributes = List of acceptable Attributes
    - Elements = List of acceptable Elements
    - ParentTags = List of Text Tags for parents, all the way up
    - XML = list of XML objects containing nested elements (if any)

    Dict:dictionary    
     - Attributes = dictionary of element attributes (if any) as strings
     - Data = Conversion into data dictionary
     - Schema = dictionary of element type schema for main Element

    Obj:RJE_Objects
    - Parent = XML object that corresponds to parent Element.
    - Processor = Object to process each level one XML object, which is then deleted.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Content','Source']
        self.optlist = []
        self.statlist = []
        self.listlist = ['Attributes','Elements','ParentTags','XML']
        self.dictlist = ['Attributes','Data','Schema']
        self.objlist = ['Parent','Processor']
        ### Defaults ###
        self._setDefaults(info='',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
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
                self._cmdRead(cmd,type='info',att='Source',arg='parse')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'list',['Attributes','Elements'])     
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### XML Parsing Methods                                                                                     #
#########################################################################################################################
    def parseXML(self,xmlfile=None,url=None):      ### XML Parsing method. Populates XML List and returns.
        '''
        XML Parsing method. Populates XML List and returns.
        >> xmlfile:str = XML text file to parse
        >> url:str = URL to parse instead of file
        << xml_list:list of XML objects to return
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xmlfile and xmlfile.lower() != 'none': self.info['Source'] = xmlfile     # Source file for XML data
            elif url: self.info['Source'] = url; xmlfile = ''
            self.obj['Handler'] = XMLHandler()                  # XML Handler Object
            p = xml.sax.make_parser()                           # XML Parser
            self.obj['Handler'].xml = self                      # Make sure Handler can refer to and populate self
            p.setContentHandler(self.obj['Handler'])            # Set Handler for Parser
            ### ~ [2] Read and Parse data, one line at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Source'): X = open(self.info['Source'],'r')
            else: X = urllib2.urlopen(url)
            while True:
                r = X.readline()
                if not r: break
                p.feed(r)
                if self.getAttribute('opt','Test') and len(self.list['XML']) == 20: raise ValueError
            p.close()
            X.close()
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self
        except:
            self.log.errorLog('Error in XML.parseXML(%s)' % self.info['Name'],printerror=True,quitchoice=True)
            self.obj['Handler'].makeSchema()
            print( self.dict['Schema'])
            self.outputSchema()
            return None
#########################################################################################################################
    def getXML(self,taglist):   ### Returns XML object matching taglist
        '''Returns XML object matching taglist.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            intags = taglist[0:]
            myxml = self
            ### ~ [2] Delve into Elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while intags:
                tag = intags.pop(0)
                for child in myxml.list['XML']:             #!# Should make a list! #!#
                    #print tag, child.info['Name'], '?'
                    if child.info['Name'] == tag:
                        myxml = child
                        break
                if myxml.info['Name'] == tag: continue
                self.log.errorLog('Cannot find "%s" object from %s' % (tag,rje.join(taglist,' > ')),printerror=False)
                #print self.info
                #print self.list
                #print self.dict
                return None
            return myxml
        except: self.log.errorLog('Problem with XML.getXML()')            
#########################################################################################################################
    ### <3> ### XML Detail output                                                                                       #
#########################################################################################################################
    def outputSchema(self,format='txt',filename='schema.txt'):  ### Formats and outputs shema
        '''
        Formats and outputs shema.
        >> format:str [txt] = Type of output format
        >> filename:str [schema.txt] = Name for output file
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            schema = self.dict['Schema']
            rje.backup(self,filename)
            level = 1
            ### ~ [2] Process Schema ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.schedict(schema,filename,level)            

        except:
            self.log.errorLog('Problem outputting schema')
            print( schema)
#########################################################################################################################
    def schedict(self,schema,filename,level,tag='xml'):
        children = schema['<>']
        atts = schema[':']
        print( '%s%s [%d] : %s' % ('|--' * level,tag,level,rje.join(atts,',')))
        open(filename,'a').write('%s%s [%d] : %s\n' % ('|--' * level,tag,level,rje.join(atts,',')))
        for child in children: self.schedict(schema[child],filename,level+1,child)
#########################################################################################################################
### End of SECTION II: XML Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: XMLHandler Class                                                                                       #
#########################################################################################################################
class XMLHandler(xml.sax.ContentHandler):   ### Adapted from Python in a Nutshell
    #xml = None     # XML Object that ultimately stores all the data
    #parsing = []   # List XML objects currently being parsed (in order of parentage)
    #schatts = {}   # Dictionary of schema attributes
    def printLog(self,type,text,newline=True,log=True,clear=0): self.xml.log.printLog(type,text,newline=newline,log=log,clear=clear)
    def errorLog(self,errortext,printerror=True): self.xml.log.errorLog(errortext,printerror=printerror)
#########################################################################################################################
    def startDocument(self):    ### Called when document begins
        '''Called when document begins.'''
        #x#self.printLog('#PARSE','Parsing %s' % self.xml.info['Name'],newline=False,log=False)
        self.parsing = []
        self.schemalist = []
        self.schatts = {}
        self.x = 0  # Higher level Elements
        self.e = 0  # Number of elements
        self.r = 0  # Number of Elements retained
#########################################################################################################################
    def endDocument(self):    ### Called when document ends
        '''Called when document begins.'''
        ### ~ [1] Finish parsing (and check complete) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.printLog('\r#PARSE','Parsing %s complete!' % self.xml.info['Name'])
        if self.parsing: raise ValueError
#########################################################################################################################
    def makeSchema(self):   ### Generate schema dictionary from schemalist 
        #print self.schatts
        schema = {'<>':[],':':[]}
        for slist in self.schemalist:
            sdict = schema
            schjoin = rje.join(slist,':')
            while slist:
                tag = slist.pop(0)
                if tag not in sdict['<>']: sdict['<>'].append(tag)
                if tag not in sdict: sdict[tag] = {'<>':[]}
                sdict = sdict[tag]
            sdict[':'] = self.schatts[schjoin][0:]
        self.xml.dict['Schema'] = schema
#########################################################################################################################
    def startElement(self,tag,attributes):  ### Called when a new element begins
        ### ~ [1] Generate XML object for element ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.printLog('\r#PARSE','Parsing %s: %s elements (%s level 1; %s retained)' % (self.xml.info['Name'],rje.integerString(self.e),rje.integerString(self.x),rje.integerString(self.r)),False,False)
        self.e += 1
        if not self.parsing:
            myxml = self.xml   # Very first Element goes in main XML Object
            myxml.list['ParentTags'] = [tag]
        else:
            myxml = XML(log=self.xml.log,cmd_list=self.xml.cmd_list)
            myxml.obj['Parent'] = self.parsing[-1]
            myxml.list['ParentTags'] = myxml.obj['Parent'].list['ParentTags'] + [tag]
            if self.parsing[-1] == self.xml: self.x += 1
        if tag in self.xml.list['Elements'] or not self.xml.list['Elements']:
            if myxml.obj['Parent']: self.parsing[-1].list['XML'].append(myxml)
            if myxml.list['ParentTags'] not in self.schemalist: self.schemalist.append(myxml.list['ParentTags'])
        myxml.info['Name'] = tag
        self.parsing.append(myxml)
        myxml.stat['Level'] = len(self.parsing)
        ### ~ [2] Update Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        schjoin = rje.join(myxml.list['ParentTags'],':')
        if schjoin not in self.schatts: self.schatts[schjoin] = []
        for q in attributes.getQNames():
            if q in self.xml.list['Attributes'] or not self.xml.list['Attributes']:   # Only add if wanted
                myxml.dict['Attributes'][q] = attributes.getValueByQName(q)
                if q not in self.schatts[schjoin]: self.schatts[schjoin].append(q)
#########################################################################################################################
    def endElement(self,tag):   ### Called when Element ends. Check last XML object and update.
        '''Called when Element ends. Check last XML object and update.'''
        ### ~ [1] Check Last Tag matches this one ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        myxml = self.parsing.pop(-1)    # Remove from current list
        if tag != myxml.info['Name']: raise ValueError      # Check it is the last Element
        self.r += 1
        if len(self.parsing) == 1 and self.xml.obj['Processor']:   # Try to process Level 1 Element and Delete
            try:
                self.xml.obj['Processor'].processXML(myxml)
                self.xml.list['XML'].remove(myxml)
                self.r -= 1
            except: self.errorLog('Problem processing XML with %s.processXML()' % self.xml.obj['Processor'].info['Name'])
        #x#print 'End', tag, myxml.dict['Attributes'], myxml.info
#########################################################################################################################
    def characters(self,data):  ### Process text content
        '''Process text content.'''
        myxml = self.parsing[-1]
        content = rje.join(rje.split(data))
        if content and myxml.info['Content']: myxml.info['Content'] = rje.join([myxml.info['Content'],content],'\n')
        elif content: myxml.info['Content'] = content
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
    try: XML(mainlog,cmd_list).parseXML()
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
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
