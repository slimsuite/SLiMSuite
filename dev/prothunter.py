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
Module:       ProtHunter
Description:  Protein Hunter
Version:      0.0.0
Last Edit:    29/01/15
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    ProtHunter is essentially a wrapper for plugging together four processing steps into a single pipeline:

    1. Reformatting and/or translating a set of protein/DNA target sequences ready for searching.
    2. GABLAM or HMM search of candidate proteins/domains in target sequences.
    3. MultiHAQ annotation of candidate targets using additional protein databases.
    4. Generation of summary HTML files for easy data exploration.

    During development, not all of these functions will be available.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=SEQFILE       # File containing target protein/DNA sequences
    candidates=SEQFILE  # File containing candidate protein sequences to go hunting for []

    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X  # This will control output file names and also be used to pick up later analyses for HTML etc. []

    ### ~ TARGET PROCESSING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    ### ~ HTML OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    stylesheets=LIST    : List of CSS files to use ['http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css']


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
import rje, rje_db, rje_html, rje_obj, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
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
    (program, version, last_edit, copy_right) = ('ProtHunter', '0.0.0', 'January 2015', '2015')
    description = 'Protein Hunter'
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
### SECTION II: ProtHunter Class                                                                                        #
#########################################################################################################################
class ProtHunter(rje_obj.RJE_Object):
    '''
    ProtHunter Class. Author: Rich Edwards (2015).

    Str:str
    
    Bool:boolean

    Int:integer

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
        self.strlist = ['Candidates','SeqIn']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['HTML'] = rje_html.HTML(self.log,self.cmd_list)
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
                self._cmdReadList(cmd,'file',['Candidates','SeqIn'])  # String representing file path
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
            # Reformat sequences
            # GABLAM search or HMM search to identify domains
            # MultiHAQ Search
            # HTML generation
            self.makeHTML()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Basefile'): self.basefile(rje.baseFile(self.getStr('SeqIn')))
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
    ### <3> ### HTML Generation Methods                                                                                 #
#########################################################################################################################
    def makeHTML(self): ### Generates HTML pages for interactive navigation.
        '''Generates HTML pages for interactive navigation.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = self.basefile()
            scmd = self.cmd_list + ['seqin=%s' % self.getStr('Candidates'),'autoload=T','autofilter=F','seqmode=file']
            candseq = rje_seqlist.SeqList(self.log,scmd)
            # All files and directories are named after basefile:
            # *.fas = original target PROTEIN sequences (with original descriptions)
            scmd = self.cmd_list + ['seqin=%s' % self.getStr('SeqIn'),'autoload=T','autofilter=F','seqmode=file']
            seqlist = rje_seqlist.SeqList(self.log,scmd)
            # *.gablam.tdt = GABLAM results with match details. (Might have *.hmmer.tdt instead.)
            gdb = self.db().addTable('%s.gablam.tdt' % basefile,mainkeys=['Qry','Hit'],name='gablam',expect=False)
            # - Contains candidate proteins as Queries and Target proteins as hits
            # *.HAQESAC/ = directory containing individual HAQESAC runs, named after Hit accnum
            haqdir = rje.makePath('./%s.HAQESAC/' % basefile)

            ### ~ [2] Generate front page ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s.html' % basefile
            hobj = self.obj['HTML']
            hobj.list['StyleSheets'] = ['http://www.slimsuite.unsw.edu.au/stylesheets/rje_tabber.css',
                                        'http://www.slimsuite.unsw.edu.au/stylesheets/slimhtml.css']
            html = hobj.htmlHead(basefile)
            # Front page should have:
            html += '<h1>%s</h1>\n\n' % basefile
            htabs = []      # (tab_id, tab_html_text[, tab_title])
            # Target protein list (with links to HAQ HTML)
            ctext = '%s\n' % string.join(['Name','Descripton','Length'],'\t')
            seqdict = seqlist.makeSeqNameDic('short')
            if gdb: hitlist = gdb.indexKeys('Hit')
            else: hitlist = rje.sortKeys(seqdict)
            for name in hitlist:
                seq = seqdict[name]
                cseq = [name,seqlist.seqDesc(seq),'%s aa' % seqlist.seqLen(seq)]
                acc = seqlist.seqAcc(seq)
                if os.path.exists('%s%s.log' % (haqdir,acc)):
                    cseq[0] = '<a href="%s%s.html">%s</a>' % (haqdir,acc,cseq[0])
                ctext += '%s\n' % string.join(cseq,'\t')
            htabs.append(('Hits',rje_html.tableToHTML(ctext,'\t',tabid='parse'),'Target sequences hit by candidates.'))
            # GABLAM/HMM table (with above links)
            if gdb:
                ctext = '%s\n' % string.join(gdb.fields(),'\t')
                for gline in open('%s.gablam.tdt' % basefile,'r').readlines()[1:]:
                    gdata = string.split(gline,'\t')
                    acc = string.split(gdata[0],'__')[-1]
                    gdata[0] = '<a href="http://www.uniprot.org/uniprot/%s" target="_blank">%s</a>' % (acc,gdata[0])
                    acc = string.split(gdata[1],'__')[-1]
                    gdata[1] = '<a href="%s%s.html">%s</a>' % (haqdir,acc,gdata[1])
                    ctext += '%s\n' % string.join(gdata,'\t')
                htabs.append(('GABLAM',rje_html.tableToHTML(ctext,'\t',tabid='parse'),'GABLAM hit table.'))
            # Candidate list (with DB links)
            if candseq.seqNum():
                ctext = '%s\n' % string.join(['AccNum','ID','Descripton','Length'],'\t')
                accdict = candseq.makeSeqNameDic('accnum')
                for acc in rje.sortKeys(accdict):
                    seq = accdict[acc]
                    cseq = [acc,candseq.seqID(seq),candseq.seqDesc(seq),'%s aa' % candseq.seqLen(seq)]
                    cseq[0] = '<a href="http://www.uniprot.org/uniprot/%s" target="_blank">%s</a>' % (acc,acc)
                    ctext += '%s\n' % string.join(cseq,'\t')
                htabs.append(('Candidates',rje_html.tableToHTML(ctext,'\t',tabid='parse'),'Candidate sequences to search.'))
            html += hobj.tabberHTML('GABLAM',htabs)
            html += hobj.htmlTail()
            open(hfile,'w').write(html)

            ### ~ [3] Generate sequence-specific pages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #?# Move this to HAQESAC or MultiHAQ
            for i in range(len(hitlist)):
                hit = string.split(hitlist[i],'__')[-1]
                logfile = '%s%s.log' % (haqdir,hit)
                seqbase = logfile[:-4]
                hfile = '%s.html' % seqbase
                html = hobj.htmlHead(seqbase)
                # Front page should have:
                html += '<h1>%s</h1>\n\n' % seqbase
                html += '<p>Click <a href="../%s.html">here</a> to return to results summary. \n' % basefile
                if i: html += 'Previous: <a href="./%s.html"><code>%s</code></a>. \n' % (string.split(hitlist[i-1],'__')[-1],hitlist[i-1])
                if i < len(hitlist)-1: html += 'Next: <a href="./%s.html"><code>%s</code></a>. \n' % (string.split(hitlist[i+1],'__')[-1],hitlist[i+1])
                html += '</p>\n'
                htabs = []      # (tab_id, tab_html_text[, tab_title])
                for ftype in ['png','tree.txt','fas','nwk','log']:
                    seqfile = '%s.%s' % (seqbase,ftype)
                    if not os.path.exists(seqfile): continue
                    tabtext = '<p><a href="./%s">./%s</a></p>\n' % (os.path.basename(seqfile),os.path.basename(seqfile))
                    if ftype == 'png':
                        tabtext += '<a href="./%s"><img src="%s" width="100%%"></a>\n' % (os.path.basename(seqfile),os.path.basename(seqfile))
                        tabdesc = 'PNG of %s tree.' % seqbase
                    else:
                        tabtext += '<pre>%s</pre>\n' % open(seqfile,'r').read()
                        if ftype == 'tree.txt':
                            for xref in hitlist:
                                reptext = '<a href="./%s.html">%s</a>' % (string.split(xref,'__')[-1],xref)
                                tabtext = string.replace(tabtext,': %s ' % xref,': %s ' % reptext)
                            while rje.matchExp('(: \S+_(\S+)__(\S+) )',tabtext):
                                (oldtext,sid,spec,spacc) = rje.matchExp('(: (\S+)_(\S+)__(\S+) )',tabtext)
                                newtext = ': %s_<a href="http://www.uniprot.org/taxonomy/?query=%s&sort=score" target="_blank">%s</a>__<a href="http://www.uniprot.org/uniprot/%s" target="_blank">%s</a> ' % (sid,spec,spec,spacc,spacc)
                                tabtext = string.replace(tabtext,oldtext,newtext)
                        tabdesc = '%s output' % seqfile
                    htabs.append((ftype,tabtext,tabdesc))
                if htabs: html += hobj.tabberHTML(os.path.basename(seqbase),htabs)
                else: html += '<p><i>No output found for <code>%s</code>!</i></p>\n' % hit
                html += hobj.htmlTail()
                open(hfile,'w').write(html)
        except: self.errorLog('Problem with %s.makeHTML()' % self.prog())
#########################################################################################################################
### End of SECTION II: ProtHunter Class                                                                                 #
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
    try: ProtHunter(mainlog,cmd_list).run()

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
