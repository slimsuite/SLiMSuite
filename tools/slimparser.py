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
Module:       SLiMParser
Description:  SLiMSuite REST output parsing tool.
Version:      0.5.0
Last Edit:    27/09/17
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for parsing the full REST output of a program into a couple of dictionaries, with options to output
    the data to files or convert to/from json format.

    If `restin=FILE` is a URL, this will be interpreted as a REST command for API access. Use with `rest=X` and
    `pureapi=T` to print the output to STDOUT once the run is complete. Use in conjunction with `v=-1` to avoid
    additional STDOUT output and `silent=T` to avoid log generation.

    REST URLs can include files to be uploaded. These must be prefixed with `file:`, e.g. `&seqin=file:input.fas`. If the
    specified file exists then the content will replace the file name in the REST call.

Commandline:
    ### ~ INPUT/OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    restin=FILE     # Full REST output file, REST jobid (number), or full REST call URL (http://...) []
    password=X      # Optional password for REST jobid retrieval [None]
    restout=T/F     # Whether to save extracted elements to individual files [False]
    rest=X          # Return text for just the REST output element X []
    restbase=X      # Basefile for parsed REST output that lacks defined filename [jobid]
    restoutdir=PATH # Path for output of parsed REST output [./]
    pureapi=T/F     # Whether to return the text returned from the REST call. Needs rest=X. [False]
    ### ~ REST API Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    resturl=X       # URL of rest server ['http://rest.slimsuite.unsw.edu.au/']
    refresh=X       # Initial number of seconds for between checks for job status (will double) [5]
    maxrefresh=X    # Maximum number of seconds for incrementing check refresh [600]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import json, os, string, sys, time, urllib2
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_html
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.0.1 - Fixed RestKeys bug.
    # 0.1.0 - Added retrieval and parsing of existing server job. Added password.
    # 0.2.0 - Added API access to REST server if restin is REST call (i.e. starts with http:)
    # 0.2.1 - Added PureAPI output of API REST call returned text.
    # 0.3.0 - Added parsing of input files to give to rest calls.
    # 0.3.1 - Fixed issue that had broken REST server full output.
    # 0.3.2 - Fixed issue reading files for full output.
    # 0.3.3 - Tidied output names when restbase=jobid.
    # 0.3.4 - Tweaked error messages.
    # 0.4.0 - Added simple json format output.
    # 0.5.0 - Added restkeys/outputs output of REST output keys.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program: parse saved output from REST server.
    # [Y] : Add capacity to call rest server and parse results directly.
    # [ ] : Add json conversion.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SLiMParser', '0.5.0', 'September 2017', '2014')
    description = 'SLiMSuite REST output parsing tool'
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
### SECTION II: SLiMParser Class                                                                                        #
#########################################################################################################################
class SLiMParser(rje_obj.RJE_Object):
    '''
    SLiMParser Class. Author: Rich Edwards (2014).

    Str:str
    Password=X      : Optional password for REST jobid retrieval [None]
    RestIn=FILE     : Full REST output file or REST jobid (number). []
    Rest=X          : Return text for just the REST output element X []
    RestBase=X      : Basefile for parsed REST output that lacks defined filename [rest]
    RestOutDir=PATH : Path for output of parsed REST output [./]
    RestURL         : URL of rest server ['http://rest.slimsuite.unsw.edu.au/']

    Bool:boolean
    PureAPI=T/F     # Whether to return the text returned from the REST call. Needs rest=X. [False]
    RestOut=T/F     : Whether to save extracted elements to individual files [False]

    Int:integer
    MaxRefresh=X    # Maximum number of seconds for incrementing check refresh [600]
    Refresh=X       # Initial number of seconds for between checks for job status (will double) [5]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - RestKeys = List of REST outfmt keys in order read in.

    Dict:dictionary
    - Output = REST output text parsed from RestIn
    - Outfile = REST output files parsed from RestIn

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Password','RestIn','Rest','RestBase','RestOutDir','RestURL']
        self.boollist = ['PureAPI','RestOut']
        self.intlist = ['MaxRefresh','Refresh']
        self.numlist = []
        self.filelist = []
        self.listlist = ['RestKeys']
        self.dictlist = ['Output','Outfile']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'RestOutDir':rje.makePath('./'),'RestURL':'http://rest.slimsuite.unsw.edu.au/'})
        self.setBool({'PureAPI':False,'RestOut':False})
        self.setInt({'MaxRefresh':600,'Refresh':5})
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
                self._cmdRead(cmd,type='str',att='RestIn',arg='jobid')
                self._cmdReadList(cmd,'str',['Password','Rest','RestIn','RestURL'])   # Normal strings
                self._cmdReadList(cmd,'path',['RestOutDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['RestBase'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['PureAPI','RestOut'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxRefresh','Refresh'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with SLiMParser cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,maxparsesize=0):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('RestIn').startswith('http:') and not self.restAPI(): return 'ERROR: REST call failed.'
            restparse = self.parse()
            if not restparse: return 'ERROR: Rest Parsing error.'
            if self.getBool('PureAPI'): return restparse
            elif self.getBool('RestOut'): self.save()
        except: self.errorLog('Problem with main SLiMParser Run()')
        return self.restOutput(maxparsesize=maxparsesize)
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Check and modify URL if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('RestIn').startswith('http:'):
                #!# Check for rest URL and add if missing
                #!# Split on &
                restcmd = string.split(self.getStr('RestIn'),'&')
                for i in range(len(restcmd)):
                    if '=' not in restcmd[i]: continue
                    (opt,value) = string.split(restcmd[i],'=',1)
                    if value.startswith('file:'):   # Conversion of cmd=file:FILE into cmd=CONTENT
                        rfile = string.split(value,':',1)[1]
                        #!# Consider adding max size constraint. Probably a URL size limit.
                        if rje.exists(rfile):
                            restcmd[i] = '%s=%s' % (opt,rje.chomp(string.join(open(rfile,'r').readlines(),'\\n')))
                            if '&' in restcmd[i]:
                                self.warnLog('%s "&" => "+" conversions for %s.' % (rje.iStr(restcmd[i].count('&')),rfile))
                                restcmd[i] = string.replace(restcmd[i],'&','+')
                        else: self.warnLog('File "%s" not found.' % rfile,quitchoice=True)
                self.setStr({'RestIn':string.join(restcmd,'&')})
            ## ~ [1b] Direct Parsing of output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:   # Convert to file
                self.setStr({'RestIn':rje.makePath(self.getStr('RestIn'),True)})
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### REST Output parsing methods                                                                             #
#########################################################################################################################
    def parse(self):    ### Parse REST file into dictionaries
        '''Parse REST file into dictionaries.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['RestKeys'] = []
            rbase = '%s%s' % (self.getStr('RestOutDir'),rje.baseFile(self.getStr('RestBase'),strip_path=True,keepext=True))
            if rje.exists(self.getStr('RestIn')): restin = open(self.getStr('RestIn'),'r').read()
            elif rje.matchExp('^(\d+)$',self.getStr('RestIn')):
                url = '%sretrieve&jobid=%s&password=%s' % (self.getStr('RestURL'),self.getStr('RestIn'),self.getStr('Password'))
                if self.getBool('PureAPI') and self.getStrLC('Rest'): url += '&rest=%s' % (self.getStr('Rest'))
                else: url += '&rest=full'
                restin = urllib2.urlopen(url).read()
                if self.getBool('PureAPI'): return restin
            else: raise IOError('%s not found!' % self.getStr('RestIn'))
            jobid = None
            ### ~ [2] Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for restdata in string.split(restin,'###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'):
                if not jobid:
                    self.dict['Output']['intro'] = restdata
                    prog = rje.matchExp('Output for (\S+)',restdata)[0]
                    self.dict['Output']['prog'] = prog
                    jobid = rje.matchExp('JobID: (\d+)',restdata)[0]
                    self.dict['Output']['jobid'] = jobid
                    if not self.getStrLC('RestBase'): rbase = '%s%s' % (self.getStr('RestOutDir'),jobid)
                    self.dict['Outfile']['jobid'] =  '%s.jobid' % (rbase)
                    continue
                restlines = string.split(restdata,'\n')
                rparse = string.split(restlines.pop(0))
                if rparse[0] != '#': self.errorLog('REST output format error: %s' % string.join(rparse),printerror=False); continue
                if rparse[1][-1] != ':': self.errorLog('REST output format error: %s' % string.join(rparse),printerror=False); continue
                rkey = rparse[1][:-1]
                try:
                    rfile = '%s.%s' % (rbase,rje.baseFile(rparse[2],strip_path=True,keepext=True))
                except: rfile = ''
                if not rfile: rfile = '%s.%s' % (rbase,rkey)
                rfile = string.replace(rfile,'%s.%s.' % (jobid,jobid),'%s.' % jobid)
                self.dict['Output'][rkey] = string.join(restlines,'\n')
                self.dict['Outfile'][rkey] = rfile
                self.list['RestKeys'].append(rkey)
            self.printLog('#PARSE','Parsed %s: %d REST outputs.' % (self.getStr('RestIn'),len(self.dict['Output'])))
            return True
        except: self.errorLog('%s.parse error' % self); return False
#########################################################################################################################
    def save(self):     ### Saves parsed REST output to files
        '''Saves parsed REST output to files.'''
        rbase = '%s%s' % (self.getStr('RestOutDir'),rje.baseFile(self.getStr('RestBase'),strip_path=True,keepext=True))
        rje.mkDir(self,self.getStr('RestOutDir'))
        outputs = rje.sortKeys(self.dict['Output'])
        if self.getStrLC('Rest') in outputs: outputs = [self.getStrLC('Rest')]
        elif self.getStrLC('Rest') in ['full','text']:
            outfile = '%s.rest' % rbase
            open(outfile,'w').write(self.restFullOutput())
            self.printLog('#OUT','%s: %s' % (self.getStrLC('Rest'),outfile))
            return True
        elif self.getStrLC('Rest'):
            self.printLog('#OUTFMT','REST output format "%s" not recognised.' % self.getStrLC('Rest'))
            if self.i() < 0 or not rje.yesNo('Output all parsed outputs?'): return False
            outfile = '%s.rest' % rbase
            open(outfile,'w').write(self.restFullOutput())
            self.printLog('#OUT','full: %s' % (outfile))
            return True
        for rkey in outputs:
            if rkey in self.dict['Outfile']:
                rje.backup(self,self.dict['Outfile'][rkey])
                open(self.dict['Outfile'][rkey],'w').write(self.dict['Output'][rkey])
                self.printLog('#OUT','%s: %s' % (rkey,self.dict['Outfile'][rkey]))
            elif rkey not in ['intro']: self.warnLog('No outfile parsed/generated for %s output' % rkey)
#########################################################################################################################
    def jsonOutput(self,outfmt=None,maxparsesize=0):    ### Returns json output for outfmt
        '''Returns json output for outfmt.'''
        return json.dumps(self.restOutput(outfmt,maxparsesize))     #!# Need to improve this!
#########################################################################################################################
    def jsonText(self,text,asjson=False):
        if asjson: return json.dumps(text)
        else: return text
#########################################################################################################################
    def restOutput(self,outfmt=None,maxparsesize=0,asjson=False):    ### Returns rest output for outfmt
        '''Returns rest output for outfmt.'''
        if not outfmt: outfmt = self.getStrLC('Rest')
        if not outfmt: self.jsonText('No REST output',asjson)
        if outfmt in self.dict['Output']:
            rfile = string.split(self.dict['Output'][outfmt],'\n')[0]
            if rje.exists(rfile):
                fext = string.split(rfile,'.')[-1]
                if fext in ['png']:
                    self.debug(rfile)
                    self.jsonText(rfile,asjson)
                nbytes = os.path.getsize(rfile)
                if nbytes > maxparsesize > 0:   # Too large to parse
                    otext = '%s is too large to return (%s > %s)' % (os.path.basename(rfile),rje.humanByteSize(nbytes),rje.humanByteSize(maxparsesize))
                    try: jobid = self.dict['Output']['jobid']
                    except: jobid = None
                    resturl = '%sretrieve&jobid=%s&rest=%s[&password=X]' % (self.getStr('RestURL'),jobid,outfmt)
                    if not jobid or outfmt == self.getStrLC('Rest'): return self.jsonText('ERROR: %s' % (otext),asjson)
                    else: return self.jsonText('%s in full output. Try %s.' % (otext,resturl),asjson)
                else:
                    delimit = rje.delimitFromExt(filename=rfile,write=False)
                    if asjson and delimit in [',','\t']:
                        jtext = []
                        for rline in open(rfile,'r').readlines():
                            jtext.append(json.dumps(rje.readDelimit(rline,delimit)))
                        return '[%s]' % string.join(jtext,',\n        ')
                    #!# Add json parsing of fasta files?
                    else:
                        outtxt = open(rfile,'r').read()
                        if not outtxt.endswith('\n'): outtxt += '\n'
                        return self.jsonText(outtxt,asjson)
            elif asjson and outfmt in self.dict['Outfile']:
                pass    #!# Sort out json formatting here based on file extension!
            return self.dict['Output'][outfmt]
        elif outfmt in ['parse','format']:
            intro = '<pre>%s</pre>\n\n' % self.restOutput('intro')
            return self.jsonText(intro,asjson)
        elif outfmt in ['default','full']: return self.jsonText(self.restFullOutput(maxparsesize),asjson)
        elif outfmt in ['restkeys','outputs']: return string.join(self.list['RestKeys']+[''],'\n')
        return self.jsonText('No %s output generated.' % outfmt,asjson)
#########################################################################################################################
    def restFullOutput(self,maxparsesize=0):   ### Returns full REST output from file
        '''Returns full REST output from file.'''
        if rje.exists(self.getStr('RestIn')) and not self.force(): return open(self.getStr('RestIn'),'r').read()
        try: jobid = self.dict['Output']['jobid']
        except: jobid = None
        rtxt = '%s\n' % self.dict['Output']['intro']
        for rkey in self.list['RestKeys']:
            rtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            rtxt += '# %s: %s\n' % (rkey,self.dict['Outfile'][rkey])
            #?# Q. Why only if jobid? Is it not always good to replace with content?
            if jobid and rje.exists(string.split(self.dict['Output'][rkey],'\n')[0]): ### File given instead of content
                rfile = string.split(self.dict['Output'][rkey],'\n')[0]
                fext = string.split(rfile,'.')[-1]
                nbytes = os.path.getsize(rfile)
                if nbytes > maxparsesize > 0:   # Too large to parse
                    otext = '%s is too large to return (%s > %s)' % (os.path.basename(rfile),rje.humanByteSize(nbytes),rje.humanByteSize(maxparsesize))
                    resturl = '%sretrieve&jobid=%s&rest=%s[&password=X]' % (self.getStr('RestURL'),jobid,rkey)
                    rtxt += '%s in full output. Try %s.' % (otext,resturl)
                elif rfile.endswith('.png'):
                    rtxt += '%s\n' % rfile
                    #rtxt += 'Cannot return graphic in full output\n'
                #elif fext in ['htm','html']:
                #    rtxt += 'Cannot return HTML in full output\n'
                else:
                    outtxt = open(rfile,'r').read()
                    if not outtxt.endswith('\n'): outtxt += '\n'
                    rtxt += outtxt
            else: rtxt += '%s\n' % self.dict['Output'][rkey]
        return rtxt
#########################################################################################################################
    ### <3> ### REST API calling methods                                                                                #
#########################################################################################################################
    def restAPI(self):  ### Make a rest call and update RestIn with JobID if successful
        '''Make a rest call and update RestIn with JobID if successful.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            restcall = '%s&rest=jobid' % self.getStr('RestIn')
            self.printLog('#REST',restcall)
            refresh = self.getInt('Refresh')
            ### ~ [1] Set job running ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jobid = rje.chomp(urllib2.urlopen(restcall).read())
            self.printLog('#JOBID',jobid)
            ### ~ [2] Wait for completion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            checkurl = '%scheck&jobid=%s' % (self.getStr('RestURL'),jobid)
            self.printLog('#CHECK',checkurl)
            check = rje.chomp(urllib2.urlopen(checkurl).read())
            while check in ['Queued','Running']:
                self.progLog('\r#RUN',check)
                time.sleep(refresh)
                refresh = min(self.getInt('MaxRefresh'),refresh*2)
                check = rje.chomp(urllib2.urlopen(checkurl).read())
            ### ~ [3] Return JobID if finished ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if check == 'Finished':
                self.printLog('\r#RUN','REST call complete: restin=%s' % jobid)
                self.setStr({'RestIn':jobid})
                if not self.getStrLC('RestBase'): self.setStr({'RestBase':jobid})
                return jobid
            else: self.printLog('#FAIL','REST check error: %s' % check)
        except: self.errorLog('%s.restAPI error' % self)
        return False
#########################################################################################################################
### End of SECTION II: SLiMParser Class                                                                                 #
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
    try:
        slimparser = SLiMParser(mainlog,cmd_list)
        restparse = slimparser.run()
        if slimparser.getBool('PureAPI'):
            mainlog.endLog(info)
            print restparse
            return

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
