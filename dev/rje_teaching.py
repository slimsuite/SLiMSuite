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
Module:       rje_teaching
Description:  Miscellaneous utilities used for teaching and marking
Version:      0.0
Last Edit:    08/05/12
Copyright (C) 2012  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_html, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_sequence, rje_zen
import rje_html
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_TEACHING', '0.0', 'March 2012', '2012')
    description = 'Miscellaneous utilities used for teaching and marking'
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
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### SECTION II: Teaching Class                                                                                          #
#########################################################################################################################
class Teaching(rje_obj.RJE_Object):     
    '''
    Teaching Class. Author: Rich Edwards (2012).

    Str:str
    - InFile = Input filename
    - Job = Identifier determining which methods to run
    
    Bool:boolean

    Int:integer

    Num:float
    
    List:list
    - Batch = Batch list of input filenames
    - StyleSheets = List of stylesheets for HTML output

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Job','InFile','OutFile']
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.listlist = ['Batch','StyleSheets','Students']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'OutFile':'../BIOL2013.marks.tdt'})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        self.list['StyleSheets'] = ['http://www.southampton.ac.uk/~re1u06/re1u06.css']
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
                self._cmdReadList(cmd,'str',['Job'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['InFile','OutFile'])  # String representing file path 
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'inr',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['StyleSheets','Students'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Batch']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Select Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            job = self.getStr('Job')
            if job.lower() == 'test2013': self.test2013()
            elif job.lower() == 'mark2013': self.mark2013()
            elif job.lower() == 'summary2013': self.summary2013()
            elif job.lower() == 'results2013': self.mark2013(results=True)
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else: self.printLog('#JOB','Do not recognise job "%s"!' % job)
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### BIOL2013 Project Methods                                                                                #
#########################################################################################################################
    def test2013(self): ### Extracts data from BIOL2013 data file and generates formatted HTML report
        '''Extracts data from BIOL2013 data file and generates formatted HTML report.'''
        try:### ~ [1] Extract Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ekeys,edata) = self.extract2013()
            ### ~ [2] Make HTML file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '%s.htm' % rje.baseFile(self.getStr('InFile'))
            rje.backup(self,hfile)
            html = [rje_html.htmlHead(self.getStr('InFile'),self.list['StyleSheets'],tabber=False,frontpage=True,nobots=True,keywords=[],javascript='')]
            html.append('<h1>Auto-extraction of data from %s</h1>' % self.getStr('InFile'))
            ### Main content ###
            sections = ['Student Name','Species','EnsFasta','5\' UTR Length','Fusion ORF data','Annotated fusion protein sequence']
            for ekey in ekeys:
                if ekey in sections: html.append('<h2>Section %d</h2>' % (sections.index(ekey) + 1))
                if '\n' in edata[ekey]:
                    html.append('<p><b>%s</b>:</p>' % (ekey))
                    html+= ['<font face="Courier New"><pre>',edata[ekey],'</font></pre>','']
                else: html.append('<li><b>%s</b>: %s' % (ekey,edata[ekey]))

            html += [rje_html.htmlTail(tabber=False)]
            open(hfile,'w').write(string.join(html,'\n'))
            self.printLog('#HTML','HTML output: %s' % hfile)
        except: self.errorLog('%s.test2013 error' % self)
#########################################################################################################################
    def mark2013(self,results=False): ### Extracts data from BIOL2013 data files, marks and returns output
        '''Extracts data from BIOL2013 data files, marks and returns output.'''
        try:### ~ [1] Load in general data and setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edatakeys = ['Student Name','Student Number','Student Username',
                         'Species','Common Name','Ensembl Transcript ID','Ensembl Gene ID','Gene Name (MGI Symbol)','Gene Description','Transcript Length (bp)','Protein Length (aa)','No. Exons',
                         'EnsFasta',
                         '5\' UTR Length','ORF Start','ORF End','Reading Frame','Annotated Context','Annotated Strength','AIC position (start of codon)','AIC residue (amino acid) position','AIC Context','AIC Strength',
                         'Fusion ORF data','400nt from CMV30 primer','800nt from the CMV24 primer',
                         'Annotated fusion protein sequence','Annotated fusion pI','Annotated fusion Mw','Alternative fusion protein sequence','Alternative fusion pI','Alternative fusion Mw','Truncation length difference (aa)','Mw Difference (kDa 3sf)','Alignment'
                         ]
            if results:
                rje.mkDir(self,'Results/')
                student_summary = self.summary2013(save=False)
                full_hfile = '../BIOL2013.datafiles.htm'
                rje.backup(self,full_hfile)
                full_html = [rje_html.htmlHead('BIOL2013 Coursework Mark Summary',self.list['StyleSheets'],tabber=False,frontpage=True,nobots=True,keywords=[],javascript='')]
            #self.deBug(results)
            ## ~ [1a] Setup database object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            #Gene	ENSG	ENST	Use?	Student	Protein	CDS	UTR5	Codon	Core	Context	Class	GC	Strength	uSTOP	StartPos	GoodPos	eORF	eType	eCore	eCount	eLen	tORF	tCore	tCount	tType	tLen	RF1	RF2	BIOL3050
            tdb = db.addTable('../Mus_musculus.NCBIM37.64.BIOL2013.tdt',mainkeys=['ENST'],name='trans',expect=True)
            #self.deBug(tdb.fields())
            #University Username	First Name	Last Name	Transcript
            assignfile = '../BIOL2013.Bioinformatics.Assignments.tdt'
            udb = db.addTable(assignfile,['University Username'],name='students',expect=True)
            udb.renameField('University Username','Student')
            udb.renameField('Transcript','ENST')
            #self.deBug(udb.fields())
            tfields = ['Gene','ENSG','ENST','Student','Protein','CDS','UTR5','Codon','Core','Strength','uSTOP','StartPos','GoodPos','tORF','tCore','tLen','tType']
            rdb = db.joinTables('results',join=[(udb,'ENST'),(tdb,'ENST',tfields)],newkey=['Student'],cleanup=True,delimit='\t',empties=False,check=False,keeptable=True)
            self.deBug(rdb.data('ens1g11'))
            self.deBug(rdb.data('sen1g10'))
            for field in edatakeys: rdb.addField(field,evalue=0)
            for sx in range(7): rdb.addField('Section%d' % sx,evalue=0)
            rdb.addField('Total',evalue=0)
            for field in ['Name','Date','Submission','Comments','File']: rdb.addField(field);
            edb = db.addTable('../BIOL2013.mart_export.tdt',['Ensembl Transcript ID','Ensembl Exon ID'],name='ensembl')
            #Ensembl Gene ID	Ensembl Transcript ID	Ensembl Protein ID	Ensembl Exon ID	Description
            edb.addField('Exons',evalue=1)
            edb.compress(['Ensembl Transcript ID'],default='sum')
            ## ~ [1b] Setup sequence data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            epath = '/home/re1u06/researchfiles/SBSBINF/Projects/ActiveMain/2009-11-Initiation/2012-01-BIOL2013/'
            #cds = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%sMus_musculus.NCBIM37.64.cds.fas' % epath,'seqmode=index'])
            #utr = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%sMus_musculus.NCBIM37.64.utr5.fas' % epath,'seqmode=index'])
            cdna = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=../BIOL2013.cDNA.fas','seqmode=index'])
            ## ~ [1c] Get file list of student data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sfiles = glob.glob('BIOL20132011-1220Project20Data20File_*')    # 2011-12
            sfiles = glob.glob('BIOL2013*File_*')    # 2012-13
            self.printLog('#GLOB','%d submission text files identified' % len(sfiles))
            self.list['Students'] = string.split(string.join(self.list['Students']).lower())
            trouble = ['erc1g11']
            for si in range(len(sfiles)-1,-1,-1):
                for moveme in trouble:
                    if moveme in sfiles[si]: sfiles.append(sfiles.pop(si)); break
            self.deBug(sfiles)
            genedesc = {}
            if os.path.exists('genedesc.tdt'):
                for line in self.loadFromFile('genedesc.tdt',chomplines=True):
                    desc = string.split(line,'\t')
                    if desc[0] not in genedesc: genedesc[desc[0]] = {}
                    genedesc[desc[0]][desc[1]] = desc[2]

            ### ~ [2] Process files, marking as they are read in ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datafiles = []; failed = []; marked = []
            for sfile in sfiles:
                restxt = []
                self.printLog('#NEXT',sfile)
                if open(sfile,'r').readline()[:4] != 'Name': datafiles.append(sfile); continue
                ## ~ [2a] Process submission info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ikey = None; data = {}
                for line in self.loadFromFile(sfile,chomplines=True):
                    #self.debug(line)
                    if rje.matchExp('^Name:\s*(\S.+) \((\S+)\)',line): (data['Name'],data['Student']) = rje.matchExp('^Name:\s*(\S.+) \((\S+)\)',line)
                    elif rje.matchExp('^Date Submitted:\s*(\S.+)$',line): data['Date'] = rje.matchExp('^Date Submitted:\s*(\S.+)$',line)[0]
                    elif not line: continue
                    elif not string.split(line): continue
                    elif line[:9] == 'Comments:': ikey = 'Comments'
                    elif string.split(line)[0] == 'Submission': ikey = 'Submission'
                    elif string.split(line)[0] == 'Files:': ikey = 'File'
                    elif line and ikey:
                        if ikey not in data: data[ikey] = line
                        else: data[ikey] = data[ikey] + '\n' + line
                    #self.deBug(data)
                if 'Submission' in data:
                    while rje.matchExp('(<p.+</p>)',data['Submission']): data['Submission'] = string.replace(data['Submission'],rje.matchExp('(<p.+</p>)',data['Submission'])[0],'')
                    while rje.matchExp('(###.+</div>)',data['Submission']): data['Submission'] = string.replace(data['Submission'],rje.matchExp('(###.+</div>)',data['Submission'])[0],'')
                    while rje.matchExp('(<a.+</a>)',data['Submission']): data['Submission'] = string.replace(data['Submission'],rje.matchExp('(<a.+</a>)',data['Submission'])[0],'')
                    while rje.matchExp('(<span.+</span>)',data['Submission']): data['Submission'] = string.replace(data['Submission'],rje.matchExp('(<span.+</span>)',data['Submission'])[0],'')
                    data['Submission'] = string.replace(data['Submission'],'\n','; ')
                    if data['Submission'] and data['Submission'] != 'There is no student submission text data for this assignment.': print data['Submission']; #self.deBug(data)
                else: data['Submission'] = 'There is no student submission text data for this assignment.'
                if 'Comments' in data: data['Comments'] = string.replace(data['Comments'],'\n','; ')
                try: data['File'] = rje.matchExp('.+Filename:\s*(\S.+)$',data['File'])[0]
                except: self.printLog('#NULL','No submission file recognised for %s' % sfile)
                self.deBug(data)
                ## ~ [2b] Try to process data file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    self.printLog('#>>>','%s <<<#' % data['Name'],timeout=False)
                    restxt.append(('#>>>','%s <<<#' % data['Name']))
                    if self.list['Students'] and data['Student'] not in self.list['Students']: self.printLog('#SKIP','Student %s not in list' % data['Student']); continue
                    for dkey in ['Name','Student','Date','File']:
                        self.printLog('#%s' % dkey.upper()[:4],data[dkey])
                        restxt.append(('#%s' % dkey.upper()[:4],data[dkey]))
                    self.str['InFile'] = data['File']
                    #self.deBug(self.str['InFile'])
                    if self.str['InFile'][-4:] in ['.odt','docx','.rtf','.bad']:
                        if os.path.exists('%s.txt' % self.str['InFile']): self.str['InFile'] = '%s.txt' % self.str['InFile']
                        else: raise ValueError
                    elif os.path.exists('%s.bad' % self.str['InFile']): self.str['InFile'] = '%s.bad' % self.str['InFile']
                    if os.path.exists(self.str['InFile']): (ekeys,edata) = self.extract2013(logdata=False)
                    elif os.path.exists(string.replace(self.str['InFile'],'Datafile','Data File')):
                        self.str['InFile'] = string.replace(self.str['InFile'],'Datafile','Data File')
                        (ekeys,edata) = self.extract2013(logdata=False)
                    else: (ekeys,edata) = (edatakeys,{})
                except:
                    self.errorLog('%s Extraction Problem' % data['File'])
                    (ekeys,edata) = (edatakeys,{})
                if not edata:
                    self.printLog('#FAIL','Failed to process "%s"' % self.str['InFile'])
                    self.printLog('#FILE','No data file processed for %s' % sfile)
                    failed.append(data)
                    self.printLog('#>>>','<<<#',timeout=False)
                    continue
                #self.deBug(edata)
                ## ~ [2c] Reformat data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try:
                    self.printLog('#MARK','Marking %s (%s) entry' % (data['Name'],data['Student']))
                    for ekey in edata.keys():
                        edata[ekey] = string.replace(edata[ekey],'"','')
                    for ekey in ['Annotated fusion Mw','Alternative fusion Mw','Mw Difference (kDa 3sf)']:
                        edata[ekey] = string.replace(edata[ekey],',','')
                        edata[ekey] = string.replace(edata[ekey],'kDa','')
                        edata[ekey] = string.replace(edata[ekey],'Da','')
                        try: edata[ekey] = float(edata[ekey])
                        except:
                            self.printLog('#ERR','Format error: %s = %s -> 0.0' % (ekey,edata[ekey]))
                            restxt.append(('#ERR','Format error: %s = %s -> 0.0' % (ekey,edata[ekey])))
                            edata[ekey] = 0.0
                    edata['Reading Frame'] = string.replace(edata['Reading Frame'],'Frame ','')
                    for ekey in ['Transcript Length (bp)','Protein Length (aa)','No. Exons','5\' UTR Length','ORF Start','ORF End','Reading Frame','AIC position (start of codon)','AIC residue (amino acid) position','Truncation length difference (aa)']:
                        edata[ekey] = string.replace(edata[ekey],',','')
                        try: edata[ekey] = int(edata[ekey])
                        except:
                            self.printLog('#ERR','Format error: %s = %s -> 0' % (ekey,edata[ekey]))
                            restxt.append(('#ERR','Format error: %s = %s -> 0' % (ekey,edata[ekey])))
                            edata[ekey] = 0
                    for ekey in edatakeys:
                        if ekey not in edata:
                            self.printLog('#ERR','%s data missing' % (ekey))
                            restxt.append(('#ERR','%s data missing' % (ekey)))
                            edata[ekey] = ''
                ## ~ [2d] Mark processed file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    entry = rdb.data(data['Student'])
                    self.debug(entry)
                    if not entry:
                        self.deBug(rdb.dataKeys())
                        raise ValueError('Cannot find student "%s" in results table. Check assignments file (%s).' % (data['Student'],assignfile))
                    if entry['Total']:
                        self.printLog('#DUP','%s already marked once!' % data['Student'])
                        x = 1
                        while '%s-%d' % (data['Student'],x) in rdb.data(): x += 1
                        dupentry = rje.combineDict({'Student':'%s-%d' % (data['Student'],x)},entry,overwrite=False,copyblanks=True)
                        rdb.addEntry(dupentry)
                    for sx in range(7): entry['Section%d' % sx] = 0
                    if entry['tType'] == 'Bad': entry['tType'] = 'Weak'
                    if entry['Strength'] == 'Bad': entry['Strength'] = 'Weak'
                    ens = edb.data(entry['ENST'])
                    seq = cdna.getSeq(entry['ENST'],format='tuple')[1].upper()
                    #self.deBug(seq)
                    for dkey in data: entry[dkey] = data[dkey]
                    # ~ Data File naming and format ~ #
                    self.printLog('#FILE','FILE TYPE AND NAMING',timeout=False)
                    restxt.append(('#FILE','File format (plain text) = 5%; File name (*_Surname_Initial.biol2013.txt) = 5%.'))
                    self.deBug(data['Student'])
                    if len(string.split(data['Name'])) > 1 and string.split(data['Name'])[-2] in ['Mohd']: surname = string.join(string.split(data['Name'])[-2:])
                    else: surname = string.split(data['Name'])[-1]
                    if '.docx.' in self.str['InFile']: entry['Section0'] = 0; restxt.append(('#FILE','Incorrect file type (docx not plain text) and name.'))
                    elif '.odt.' in self.str['InFile']: entry['Section0'] = 0; restxt.append(('#FILE','Incorrect file type (odt not plain text) and name.'))
                    elif '.rtf.' in self.str['InFile']: entry['Section0'] = 0; restxt.append(('#FILE','Incorrect file type (rtf not plain text) and name.'))
                    elif '.error.' in self.str['InFile']: entry['Section0'] = 0; restxt.append(('#FILE','Incorrect file type (rtf not plain text) and name.'))
                    elif '.bad' in self.str['InFile']:
                        if rje.matchExp('(\S*\d_%s_\S+\.biol2013\.txt.bad$)' % surname.lower(),self.str['InFile'].lower()): entry['Section0'] = 5; restxt.append(('#FILE','File format error (not plain text).'))
                        else: entry['Section0'] = 0; restxt.append(('#FILE','Incorrect file type (not plain text) and name.'))
                    elif rje.matchExp('(\S*\d_%s_\S+\.biol2013\.txt$)' % surname.lower(),self.str['InFile'].lower()): entry['Section0'] = 10; restxt.append(('#FILE','Correct file naming and type'))
                    else: entry['Section0'] = 5; restxt.append(('#FILE','File naming error'))
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 1','STUDENT INFORMATION',timeout=False)
                    restxt.append(('#SECTION 1','STUDENT INFORMATION'))
                    for ekey in edatakeys[0:3]:
                        self.printLog('#DATA','%s = %s' % (ekey,edata[ekey]))
                        restxt.append(('#DATA','%s = %s' % (ekey,edata[ekey])))
                    entry['Section1'] = 1
                    student = '%s %s' % (string.split(edata['Student Name'])[0],string.split(edata['Student Name'])[-1])
                    if not student: student = 'None'
                    #if not edata['Student Name'] or (edata['Student Name'] != data['Name'] and student != data['Name'] and not rje.yesNo('%s = %s?' % (edata['Student Name'],data['Name']),default='N')): entry['Section1'] = 0
                    if not edata['Student Name'] or (edata['Student Name'] != data['Name'] and student != data['Name']): entry['Section1'] = 0
                    entry['Student Name'] = edata['Student Name']
                    if not edata['Student Number'] or edata['Student Number'] == '00000000': entry['Section1'] = 0
                    entry['Student Number'] = edata['Student Number']
                    if edata['Student Username'].lower() != entry['Student'].lower(): entry['Section1'] = 0
                    entry['Student Username'] = edata['Student Username']
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 2','GENE/TRANSCRIPT INFORMATION',timeout=False)
                    restxt.append(('#SECTION 2','GENE/TRANSCRIPT INFORMATION'))
                    for ekey in edatakeys[3:12]:
                        self.printLog('#DATA','%s = %s' % (ekey,edata[ekey]))
                        restxt.append(('#DATA','%s = %s' % (ekey,edata[ekey])))
                    entry['Section2'] = 0
                    if edata['Species'] == 'Mus musculus': entry['Species'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Species (Mus musculus).')
                        restxt.append(('#ERR','ERROR. Species (Mus musculus).'))
                    if edata['Common Name'] == 'Mouse': entry['Common Name'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Common Name (Mouse).')
                        restxt.append(('#ERR','ERROR. Common Name (Mouse).'))
                    if edata['Ensembl Transcript ID'] == entry['ENST']: entry['Ensembl Transcript ID'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Transcript ID (%s).' % entry['ENST'])
                        restxt.append(('#ERR','ERROR. Transcript ID (%s).' % entry['ENST']))
                    if edata['Ensembl Gene ID'] == entry['ENSG']: entry['Ensembl Gene ID'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Gene ID (%s).' % entry['ENSG'])
                        restxt.append(('#ERR','ERROR. Gene ID (%s).' % entry['ENSG']))
                    if edata['Gene Name (MGI Symbol)'].upper() == entry['Gene'].upper(): entry['Gene Name (MGI Symbol)'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Gene Name (%s).' % entry['Gene'])
                        restxt.append(('#ERR','ERROR. Gene Name (%s).' % entry['Gene']))
                    ensdesc = string.replace(ens['Description'].lower(),' ','')
                    edatdesc = string.replace(edata['Gene Description'].lower(),' ','')
                    edatdesc2 = string.replace(edatdesc,'protein','')
                    if edata['Gene Description'] and (ensdesc.find(edatdesc) == 0 or ensdesc.find(edatdesc2) == 0): entry['Gene Description'] = 1; entry['Section2'] += 1
                    elif edata['Gene Description'] and edata['Gene Description'] in genedesc and ens['Description'] in genedesc[edata['Gene Description']]:
                        if genedesc[edata['Gene Description']] == 'True':
                            self.printLog('#DESC','Gene Description: "%s"; Accepted: "%s".' % (ens['Description'],edata['Gene Description']))
                        else:
                            self.printLog('#ERR','ERROR. Gene Description (%s).' % ens['Description'])
                            restxt.append(('#ERR','ERROR. Gene Description (%s).' % ens['Description']))
                    elif edata['Gene Description'] and rje.yesNo('%s\n = \n%s\n' % (edata['Gene Description'],ens['Description'])):
                        self.printLog('#DESC','Gene Description: "%s"; Accepted: "%s".' % (ens['Description'],edata['Gene Description']))
                        entry['Gene Description'] = 1; entry['Section2'] += 1
                        open('genedesc.tdt','a').write('%s\t%s\tTrue\n' % (edata['Gene Description'],ens['Description']))
                    else:
                        if edata['Gene Description']: open('genedesc.tdt','a').write('%s\t%s\tFalse\n' % (edata['Gene Description'],ens['Description']))
                        self.printLog('#ERR','ERROR. Gene Description (%s).' % ens['Description'])
                        restxt.append(('#ERR','ERROR. Gene Description (%s).' % ens['Description']))
                    if int(edata['Transcript Length (bp)']) == len(seq): entry['Transcript Length (bp)'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Transcript Length (%s).' % len(seq))
                        restxt.append(('#ERR','ERROR. Transcript Length (%s).' % len(seq)))
                    self.deBug(entry['Protein'])
                    cds_seq = ''; skey = ''
                    for sline in string.split(edata['EnsFasta'],'\n'):
                        if rje.matchExp('(\S+):KNOWN',sline): skey = rje.matchExp('(\S+):KNOWN',sline)[0]
                        elif rje.matchExp('(\S+):PUTATIVE',sline): skey = rje.matchExp('(\S+):PUTATIVE',sline)[0]
                        elif rje.matchExp('(\S+):NOVEL',sline): skey = rje.matchExp('(\S+):NOVEL',sline)[0]
                        elif sline[:1] == '>': skey = 'err'
                        elif skey == 'cds': cds_seq += string.join(string.split(sline.upper()),'')
                    if cds_seq[-3:].lower() in ['taa','tga','tag']: protlen = int(entry['Protein']) - 1
                    else: protlen = int(entry['Protein']); self.printLog('#STOP','Student %s had no STOP codon!' % entry['Student'])
                    if int(edata['Protein Length (aa)']) == protlen: entry['Protein Length (aa)'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. Protein Length (%s).' % (protlen))
                        restxt.append(('#ERR','ERROR. Protein Length (%s).' % (protlen)))
                    if int(edata['No. Exons']) == ens['Exons']: entry['No. Exons'] = 1; entry['Section2'] += 1
                    else:
                        self.printLog('#ERR','ERROR. No. Exons (%s).' % ens['Exons'])
                        restxt.append(('#ERR','ERROR. No. Exons (%s).' % ens['Exons']))
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 3','SEQUENCES',timeout=False)
                    restxt.append(('#SECTION 3','SEQUENCES'))
                    entry['EnsFasta'] = ''; entry['Section3'] = 2   # 2 marks for correct four sequences only. 2 marks per sequence.
                    eseq = {'pep':'', 'cdna':'', 'utr5':'', 'cds':'', 'err':''}
                    for sline in string.split(edata['EnsFasta'],'\n'):
                        if rje.matchExp('(\S+):KNOWN',sline): skey = rje.matchExp('(\S+):KNOWN',sline)[0]
                        elif rje.matchExp('(\S+):PUTATIVE',sline): skey = rje.matchExp('(\S+):PUTATIVE',sline)[0]
                        elif rje.matchExp('(\S+):NOVEL',sline): skey = rje.matchExp('(\S+):NOVEL',sline)[0]
                        elif sline[:1] == '>':
                            skey = 'err'; self.printLog('#ERR',sline[1:]); entry['Section3'] = 0
                            restxt.append(('#ERR',sline[1:]))
                        elif skey in eseq: eseq[skey] += string.join(string.split(sline.upper()),'')
                        else:
                            self.printLog('#ERR','Unexpected sequence type "%s"' % skey); entry['Section3'] = 0; skey = ''
                            restxt.append(('#ERR','Unexpected sequence type "%s"' % skey))
                    for skey in eseq:
                        if skey != 'err' and not eseq[skey]: entry['Section3'] = 0
                    if eseq['cdna'] == seq: entry['EnsFasta'] += 'T'; entry['Section3'] += 2
                    if seq.find(eseq['utr5']) == 0: entry['EnsFasta'] += 'U'; entry['Section3'] += 2
                    if eseq['cds'] in seq: entry['EnsFasta'] += 'C'; entry['Section3'] += 2
                    if eseq['cds'][-3:].lower() in ['taa','tga','tag']: 
                        if eseq['pep'] == rje_sequence.dna2prot(eseq['cds'])[:-1]: entry['EnsFasta'] += 'P'; entry['Section3'] += 2
                    else:
                        if eseq['pep'] == rje_sequence.dna2prot(eseq['cds']): entry['EnsFasta'] += 'P'; entry['Section3'] += 2
                    scode = {'cdna':('T','cDNA'),'cds':('C','CDS'),'pep':('P','Protein'),'utr5':('U','5\' UTR')}
                    for stype in ['cdna','cds','pep','utr5']:
                        if scode[stype][0] in entry['EnsFasta']:
                            self.printLog('#%s' % stype.upper(),'Correct %s sequence' % scode[stype][1])
                            restxt.append(('#%s' % stype.upper(),'Correct %s sequence' % scode[stype][1]))
                        else:
                            self.printLog('#%s' % stype.upper(),'ERROR. Missing or incorrect %s sequence' % scode[stype][1])
                            restxt.append(('#%s' % stype.upper(),'ERROR. Missing or incorrect %s sequence' % scode[stype][1]))
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 4','TRANSLATION AND INITIATION',timeout=False)
                    restxt.append(('#SECTION 4','TRANSLATION AND INITIATION'))
                    for ekey in edatakeys[13:23]:
                        self.printLog('#DATA','%s = %s' % (ekey,edata[ekey]))
                        restxt.append(('#DATA','%s = %s' % (ekey,edata[ekey])))
                    if int(edata['5\' UTR Length']) == len(eseq['utr5']): entry['5\' UTR Length'] = 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#ERR','ERROR. 5\' UTR length %d vs expected %d' % (int(edata['5\' UTR Length']),len(eseq['utr5'])))
                        restxt.append(('#ERR','ERROR. 5\' UTR length %d vs expected %d' % (int(edata['5\' UTR Length']),len(eseq['utr5']))))
                    if int(edata['ORF Start']) == (1+len(eseq['utr5'])): entry['ORF Start'] = 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#ERR','ERROR. ORF Start %d vs expected %d' % (int(edata['ORF Start']),1+len(eseq['utr5'])))
                        restxt.append(('#ERR','ERROR. ORF Start %d vs expected %d' % (int(edata['ORF Start']),1+len(eseq['utr5']))))
                    #self.deBug('ORF End: %d vs expected %d' % (int(edata['ORF End']),len(eseq['utr5']+eseq['cds'])-3))
                    if eseq['cds'][-3:].lower() in ['taa','tga','tag']:
                        if int(edata['ORF End']) == len(eseq['utr5']+eseq['cds'])-3: entry['ORF End'] = 1;  entry['Section4'] += 2
                        else:
                            self.printLog('#ERR','ERROR. ORF End %d vs expected %d' % (int(edata['ORF End']),len(eseq['utr5']+eseq['cds'])-3))
                            restxt.append(('#ERR','ERROR. ORF End %d vs expected %d' % (int(edata['ORF End']),len(eseq['utr5']+eseq['cds'])-3)))
                    elif int(edata['ORF End']) == len(eseq['utr5']+eseq['cds']): entry['ORF End'] = 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#ERR','ERROR. ORF End %d vs expected %d' % (int(edata['ORF End']),len(eseq['utr5']+eseq['cds'])))
                        restxt.append(('#ERR','ERROR. ORF End %d vs expected %d' % (int(edata['ORF End']),len(eseq['utr5']+eseq['cds']))))
                    rf = len(eseq['utr5']) + 1
                    while rf > 3: rf -= 3
                    if int(edata['Reading Frame']) == rf: entry['Reading Frame'] = 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#ERR','ERROR. RF %d vs expected %d' % (int(edata['Reading Frame']),rf))
                        restxt.append(('#ERR','ERROR. RF %d vs expected %d' % (int(edata['Reading Frame']),rf)))
                    if entry['tCore'] != eseq['cds'][int(entry['tORF'])-len(eseq['utr5'])-3:][:7]:
                        self.printLog('#WARN','Predicted AIC Context error: %s vs expected %s' % (entry['tCore'],eseq['cds'][int(entry['tORF'])-len(eseq['utr5'])-3:][:7]))
                        restxt.append(('#WARN','Predicted AIC Context error: %s vs expected %s' % (entry['tCore'],eseq['cds'][int(entry['tORF'])-len(eseq['utr5'])-3:][:7])))
                    entry['tType'] = 'Weak'
                    if edata['AIC Context'][:1] in 'GA':
                        if edata['AIC Context'][-1:] == 'G': entry['tType'] = 'Strong'
                        else: entry['tType'] = 'MidR'
                    elif edata['AIC Context'][-1:] == 'G': entry['tType'] = 'MidG'                      
                    entry['Strength'] = 'Weak'
                    if edata['Annotated Context'][:1] in 'GA':
                        if edata['Annotated Context'][-1:] == 'G': entry['Strength'] = 'Strong'
                        else: entry['Strength'] = 'MidR'
                    elif edata['Annotated Context'][-1:] == 'G': entry['Strength'] = 'MidG'                      
                    for epair in [('Annotated Context','Core'),('Annotated Strength','Strength'),('AIC Context','tCore'),('AIC Strength','tType')]:
                        if edata[epair[0]].lower() == entry[epair[1]].lower(): entry[epair[0]] = 1; entry['Section4'] += 2
                        else:
                            self.printLog('#ERR','ERROR. %s %s vs expected %s' % (epair[0],edata[epair[0]],entry[epair[1]]))
                            restxt.append(('#ERR','ERROR. %s %s vs expected %s' % (epair[0],edata[epair[0]],entry[epair[1]])))
                    #if edata['Annotated Context'] == entry['Core']: entry['Annotated Context'] = 1;  entry['Section4'] += 2
                    #if edata['Annotated Strength'] == entry['Strength']: entry['Annotated Strength'] = 1;  entry['Section4'] += 2
                    #if edata['AIC Context'] == entry['tCore']: entry['AIC Context'] = 1;  entry['Section4'] += 2
                    #if edata['AIC Strength'] == entry['tType']: entry['AIC Strength'] = 1;  entry['Section4'] += 2
                    if int(edata['AIC position (start of codon)']) == int(entry['tORF'])-len(eseq['utr5'])+1: entry['AIC position (start of codon)'] += 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#AIC','ERROR. AIC position (start of codon): %d vs expected %d' % (int(edata['AIC position (start of codon)']),int(entry['tORF'])-len(eseq['utr5'])+1))
                        restxt.append(('#AIC','ERROR. AIC position (start of codon): %d vs expected %d' % (int(edata['AIC position (start of codon)']),int(entry['tORF'])-len(eseq['utr5'])+1)))
                    if int(edata['AIC residue (amino acid) position']) == int(entry['tLen'])+1: entry['AIC residue (amino acid) position'] += 1;  entry['Section4'] += 2
                    else:
                        self.printLog('#AIC','ERROR. AIC residue (amino acid) position: %d vs expected %d' % (int(edata['AIC residue (amino acid) position']),int(entry['tLen'])+1))
                        restxt.append(('#AIC','ERROR. AIC residue (amino acid) position: %d vs expected %d' % (int(edata['AIC residue (amino acid) position']),int(entry['tLen'])+1)))
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 5','PLASMID ANALYSIS',timeout=False)
                    restxt.append(('#SECTION 5','PLASMID ANALYSIS'))
                    orfdata = ['Orf#','Frame','Start','End','Length','50','100','150','200','DYKDDDDK*']
                    orflines = string.split(edata['Fusion ORF data'],'\n')
                    entry['Fusion ORF data'] = 0
                    try:
                        for i in range(10):
                            if orfdata[i] in orflines[i]: entry['Fusion ORF data'] += 1; entry['Section5'] += 1
                    except:
                        self.printLog('#ORF','ERROR. Problem with ORF data "%s"...' % orfdata[i])
                        restxt.append(('#ORF','ERROR. Problem with ORF data "%s"...' % orfdata[i]))
                        for oline in orflines:
                            self.printLog('#ERR',oline)
                            restxt.append(('#ERR',oline))
                    self.printLog('#ORF','Identified %s of 10 expected ORF lines' % entry['Fusion ORF data'])
                    restxt.append(('#ORF','Identified %s of 10 expected ORF lines' % entry['Fusion ORF data']))
                    cmv30 = 'AATGTCGTAATAACCCCGCCCCGTTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCGTTTAGTGAACCGTCAGAATTAAGCTTGCGGCCGCGAATTCATCGATAGATCTGAT' + seq
                    cmv30 = cmv30[:400]
                    ecmv30 = string.replace(edata['400nt from CMV30 primer'],'\n','').upper()
                    if string.replace(edata['400nt from CMV30 primer'],'\n','').upper() == cmv30:
                        entry['400nt from CMV30 primer'] = 5; self.printLog('#CMV30','Sequence matches plasmid + cDNA')
                        restxt.append(('#CMV30','Sequence matches plasmid + cDNA'))
                    else:
                        self.printLog('#CMV30','ERROR. Sequence does not match plasmid + cDNA')
                        restxt.append(('#CMV30','ERROR. Sequence does not match plasmid + cDNA'))
                    if not entry['400nt from CMV30 primer']:
                        self.deBug('CCMV30: %s (%d)\nECMV30: %s (%d)' % (cmv30,len(cmv30),ecmv30,len(ecmv30)))

                    cmv24 = 'TATTAGGACAAGGCTGGTGGGCACTGGAGTGGCAACTTCCAGGGCCAGGAGAGGCACTGGGGAGGGGTCACAGGGATGCCACCCGGGATCACTACTTGTCATCGTCATCCTTGTAGTCGATGTCATGATCTTTATAATCACCGTCATGGTCTTTGTAGTCAGCCCGGGATCCTCTAGAGTCGACTGGTACCGAT'
                    ecmv24 = string.replace(edata['800nt from the CMV24 primer'],'\n','').upper()
                    if ecmv24.find(cmv24) == 0:
                        entry['800nt from the CMV24 primer'] = 5; self.printLog('#CMV24','First %d nt matches plasmid' % len(cmv24))
                        restxt.append(('#CMV24','First %d nt matches plasmid' % len(cmv24)))
                    else:
                        self.printLog('#CMV24','ERROR. First %d nt does not match plasmid' % len(cmv24))
                        restxt.append(('#CMV24','ERROR. First %d nt does not match plasmid' % len(cmv24)))
                    extra = rje_sequence.reverseComplement(ecmv24[len(cmv24):])
                    if extra not in seq:
                        entry['800nt from the CMV24 primer'] = 0;
                        self.printLog('#CMV24','ERROR. Sequence does not match cDNA')
                        restxt.append(('#CMV24','ERROR. Sequence does not match cDNA'))
                    else:
                        self.printLog('#CMV24','Sequence matches cDNA')
                        restxt.append(('#CMV24','Sequence matches cDNA'))
                    if len(ecmv24) != 800:
                        entry['800nt from the CMV24 primer'] = 0; self.printLog('#CMV24','ERROR. Wrong length (%d)' % len(ecmv24))
                        restxt.append(('#CMV24','ERROR. Wrong length (%d)' % len(ecmv24)))
                    else:
                        self.printLog('#CMV24','Correct length (%d)' % len(ecmv24))
                        restxt.append(('#CMV24','Correct length (%d)' % len(ecmv24)))
                    if not entry['800nt from the CMV24 primer']:
                        self.deBug('CCMV24: %s\nECMV24: %s (%d)' % (cmv24,ecmv24,len(ecmv24)))
                    entry['Section5'] += entry['800nt from the CMV24 primer'] + entry['400nt from CMV30 primer']
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Section 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    self.printLog('#SECTION 6','ORF ANALYSIS',timeout=False)
                    restxt.append(('#SECTION 6','ORF ANALYSIS'))
                    for ekey in edatakeys[-9:-1]:
                        self.printLog('#DATA','%s = %s' % (ekey,edata[ekey]))
                        restxt.append(('#DATA','%s = %s' % (ekey,edata[ekey])))
                    flag = 'SVPVDSRGSRADYKDHDGDYKDHDIDYKDDDDK'
                    fusion = string.replace(edata['Annotated fusion protein sequence'],'\n','').upper()
                    if fusion[:100] == eseq['pep'][:100] and flag in fusion:
                        entry['Annotated fusion protein sequence'] = 2; entry['Section6'] += 2
                        self.printLog('#PROT','Fusion Protein sequence OK')
                        restxt.append(('#PROT','Fusion Protein sequence OK'))
                    else:
                        self.printLog('#PROT','ERROR. Fusion Protein sequence.')
                        restxt.append(('#PROT','ERROR. Fusion Protein sequence.'))
                    #self.deBug('Fusion Mw: %s vs %s (%s)' % (edata['Annotated fusion Mw'],rje_sequence.MWt(fusion,'mono'),edata['Annotated fusion Mw']-rje_sequence.MWt(fusion,'mono')))
                    #self.deBug('Fusion Mw: %s vs %s (%s)' % (edata['Annotated fusion Mw'],rje_sequence.MWt(fusion,'ave'),edata['Annotated fusion Mw']-rje_sequence.MWt(fusion,'ave')))
                    if rje.modulus(rje.dp(edata['Annotated fusion Mw'],1) - rje.dp(rje_sequence.MWt(fusion,'ave'),1)) <= 0.11:
                        if edata['Annotated fusion pI']: entry['Annotated fusion pI'] = 2; entry['Section6'] += 2
                        entry['Annotated fusion Mw'] = 2; entry['Section6'] += 2
                    else:
                        self.printLog('#MWT','ERROR. Ann Mwt: %s vs Calc. %s -> %s vs Calc. %s' % (edata['Annotated fusion Mw'],rje_sequence.MWt(fusion,'ave'),rje.dp(edata['Annotated fusion Mw'],1),rje.dp(rje_sequence.MWt(fusion,'ave'),1)))
                        restxt.append(('#MWT','ERROR. Ann Mwt: %s vs Calc. %s -> %s vs Calc. %s' % (edata['Annotated fusion Mw'],rje_sequence.MWt(fusion,'ave'),rje.dp(edata['Annotated fusion Mw'],1),rje.dp(rje_sequence.MWt(fusion,'ave'),1))))
                    trunc = string.replace(edata['Alternative fusion protein sequence'],'\n','').upper()
                    if trunc[:1] == 'M' and trunc in fusion and flag in trunc and trunc != fusion:
                        entry['Alternative fusion protein sequence'] = 4; entry['Section6'] += 4
                        self.printLog('#PROT','Truncated Fusion Protein sequence OK.')
                        restxt.append(('#PROT','Truncated Fusion Protein sequence OK.'))
                    else:
                        self.printLog('#PROT','ERROR. Truncated Fusion Protein sequence.')
                        restxt.append(('#PROT','ERROR. Truncated Fusion Protein sequence.'))
                    if not trunc and 'M' in fusion:
                        trunc = fusion[fusion[1:].find('M')+1:]
                    #self.deBug('Alt Mw: %s vs %s' % (edata['Alternative fusion Mw'],rje_sequence.MWt(trunc,'ave')))
                    if rje.modulus(rje.dp(edata['Alternative fusion Mw'],1) - rje.dp(rje_sequence.MWt(trunc,'ave'),1)) <= 0.11:
                        if edata['Alternative fusion pI']: entry['Alternative fusion pI'] = 2; entry['Section6'] += 2
                        entry['Alternative fusion Mw'] = 2; entry['Section6'] += 2
                    else:
                        self.printLog('#MWT','ERROR. Alt Mwt: %s vs calc. %s -> %s vs calc. %s ' % (edata['Alternative fusion Mw'],rje_sequence.MWt(trunc,'ave'),rje.dp(edata['Alternative fusion Mw'],1),rje.dp(rje_sequence.MWt(trunc,'ave'),1)))
                        restxt.append(('#MWT','ERROR. Alt Mwt: %s vs calc. %s -> %s vs calc. %s ' % (edata['Alternative fusion Mw'],rje_sequence.MWt(trunc,'ave'),rje.dp(edata['Alternative fusion Mw'],1),rje.dp(rje_sequence.MWt(trunc,'ave'),1))))
                    if edata['Truncation length difference (aa)'] ==  edata['AIC residue (amino acid) position'] - 1: entry['Truncation length difference (aa)'] = 3; entry['Section6'] += 3
                    elif len(fusion) - len(trunc) == int(edata['Truncation length difference (aa)']): entry['Truncation length difference (aa)'] = 2; entry['Section6'] += 2
                    else:
                        self.printLog('#LEN','ERROR. Len Dif: Calc. %s vs %s' % ((len(fusion) - len(trunc)),int(edata['Truncation length difference (aa)'])))
                        restxt.append(('#LEN','ERROR. Len Dif: Calc. %s vs %s' % ((len(fusion) - len(trunc)),int(edata['Truncation length difference (aa)']))))
                        self.deBug(fusion)
                        self.deBug(trunc)
                    if rje.sf((edata['Annotated fusion Mw']-edata['Alternative fusion Mw'])/1000.0,3) == float(edata['Mw Difference (kDa 3sf)']):
                        entry['Mw Difference (kDa 3sf)'] = 3; entry['Section6'] += 3
                    elif '%s' % rje.sf((edata['Annotated fusion Mw']-edata['Alternative fusion Mw'])/1000.0,3) == '%s' % float(edata['Mw Difference (kDa 3sf)']):
                        entry['Mw Difference (kDa 3sf)'] = 3; entry['Section6'] += 3
                    else:
                        self.printLog('#MWT','ERROR. Mwt Dif: Calc. %s vs %s' % (rje.sf((edata['Annotated fusion Mw']-edata['Alternative fusion Mw'])/1000.0,3),float(edata['Mw Difference (kDa 3sf)'])))
                        restxt.append(('#MWT','ERROR. Mwt Dif: Calc. %s vs %s' % (rje.sf((edata['Annotated fusion Mw']-edata['Alternative fusion Mw'])/1000.0,3),float(edata['Mw Difference (kDa 3sf)']))))
                    aseqs = []
                    for aline in string.split(edata['Alignment'],'\n'):
                        if aline[:1] == '>': aseqs.append('')
                        else:
                            try: aseqs[-1] += string.replace(aline.upper(),' ','')
                            except:
                                self.printLog('#ALNERR',aline)
                                restxt.append(('#ALNERR',aline))
                    aligned = True
                    for aseq in aseqs:
                        if string.replace(aseq,'-','') == eseq['pep'] and len(aseqs) > 3: entry['Alignment'] = 5
                        if not aseq or len(aseq) != len(aseqs[0]): aligned = False
                    if entry['Alignment']:
                        self.printLog('#ALN','Query protein sequence found in Alignment')
                        restxt.append(('#ALN','Query protein sequence found in Alignment'))
                    elif aseqs:
                        self.printLog('#ALN','ERROR. Query protein sequence not found in Alignment')
                        restxt.append(('#ALN','ERROR. Query protein sequence not found in Alignment'))
                    if aligned and aseqs:
                        entry['Alignment'] += 5;  self.printLog('#ALN','Alignment aligned. (Sequences all same length.)')
                        restxt.append(('#ALN','Alignment aligned. (Sequences all same length.)'))
                    elif aseqs:
                        self.printLog('#ALN','ERROR. Alignment not aligned. (Sequences not all same length.)')
                        restxt.append(('#ALN','ERROR. Alignment not aligned. (Sequences not all same length.)'))
                    else:
                        self.printLog('#ALN','ERROR. Alignment missing or incorrect format.')
                        restxt.append(('#ALN','ERROR. Alignment missing or incorrect format.'))
                    entry['Section6'] += entry['Alignment']
                    # ~ Total ~ #
                    restxt.append('Marks')
                    entry['Total'] = 0
                    maxmark = [10,1,9,10,20,20,30]
                    for sx in range(7):
                        entry['Total'] += entry['Section%d' % sx]
                        if sx:
                            self.printLog('#SECT%d' % sx,'%d / %d' % (entry['Section%d' % sx],maxmark[sx]))
                            restxt.append(('#SECTION %d' % sx,'%d / %d' % (entry['Section%d' % sx],maxmark[sx])))
                        else:
                            self.printLog('#FILE','%d / %d' % (entry['Section%d' % sx],maxmark[sx]))
                            restxt.append(('#FILE','%d / %d' % (entry['Section%d' % sx],maxmark[sx])))
                    self.printLog('#TOTAL','%d' % entry['Total'])
                    restxt.append(('#TOTAL','%d' % entry['Total']))
                    self.printLog('#>>>','End: %s <<<#' % data['Name'],timeout=False)
                    restxt.append(('#>>>','End: %s <<<#' % data['Name']))
                    marked.append(data['Name'])
                    ### ~ [3] Generate Results HTML File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                    if results:
                        (hkeys,hdata) = self.extract2013()
                        hfile = 'Results/BIOL2013.%s.%s.htm' % (entry['Student'].lower(),entry['ENST'])
                        #self.deBug(hfile)
                        html = [rje_html.htmlHead(self.getStr('InFile'),self.list['StyleSheets'],tabber=False,frontpage=True,nobots=True,keywords=[],javascript='')]
                        html.append('<h1>BIOL2013 Results for %s (%s)<h1>' % (entry['Student'].lower(),entry['ENST']))
                        ### Main content ###
                        html.append('<h2>Auto-extraction of data from %s</h2>' % self.getStr('InFile'))
                        sections = ['Student Name','Species','EnsFasta','5\' UTR Length','Fusion ORF data','Annotated fusion protein sequence']
                        self.deBug(hkeys)
                        self.deBug(hdata['Alignment'])
                        for hkey in hkeys:
                            hdata[hkey] = string.replace(hdata[hkey],'<','&lt;')
                            if hkey in sections: html.append('<h3>Section %d</h3>' % (sections.index(hkey) + 1))
                            if '\n' in hdata[hkey] or hkey == 'Alignment':
                                self.deBug(hkey)
                                html.append('<p><b>%s</b>:</p>' % (hkey))
                                html+= ['<font face="Courier New"><pre>',hdata[hkey],'</pre></font>','']
                            else: html.append('<li><b>%s</b>: %s' % (hkey,hdata[hkey]))
                        self.deBug(html[-10:])
                        ### Assessment & Marks ###
                        html.append('<h2>Assessment notes</h2>'); marks = False
                        for result in restxt:
                            if result == 'Marks': html.append('<h2>Marks</h2><table border=1>'); marks = True
                            elif marks and 'TOTAL' in result[0]: html.append('<tr><td><b>%s</b></td><td><b>%s%%</b></td></tr>' % (result[0][1:],result[1]))
                            elif marks and '>>>' in result[0]: html += ['</table><hr>%s %s' % result]
                            elif marks: html.append('<tr><td><b>%s</b></td><td>%s</td></tr>' % (result[0][1:],result[1]))
                            elif 'SECTION' in result[0]: html.append('<h3>%s: %s</h3>' % (result[0][1:],result[1]))
                            else: html.append('<b><li>%s:</b> %s' % (result[0][1:],result[1]))
                        if student_summary:
                            html.append('<h2>Lab Report Marks</h2>')
                            html.append(student_summary[entry['Student'].lower()])
                        full_html += html[1:]
                        html += [rje_html.htmlTail(tabber=False)]
                        open(hfile,'w').write(string.join(html,'\n'))
                        self.printLog('#HTML','HTML output: %s' % hfile)
                        open(full_hfile,'w').write(string.join(full_html,'\n'))
                        self.printLog('#HTML','HTML output: %s' % full_hfile)

                #self.deBug(entry)
                except:
                    try:
                        self.errorLog('Problem processing %s file' % data['Name'])
                        failed.append(data)
                        entry['Total'] = 'ERROR'
                        #self.deBug(entry)
                    except: self.errorLog('Real problem with %s submission' % sfile)
            for failure in failed:
                if failure['Name'] not in marked: self.printLog('#FAIL',failure['Name'])
            #if not results:
            rdb.saveToFile(self.getStr('OutFile'))
            if results: self.summary2013()
        except: self.errorLog('%s.mark2013 error' % self)
#########################################################################################################################
    def summary2013(self,save=True):  ### Generates HTML of Summary Lab Report and Data File Marks
        '''Generates HTML of Summary Lab Report and Data File Marks.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            student_summary = {}
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            # Username	First Name	Last Name
            # Title	Contents	Abbrev	Rationale
            # Table of features	Plasmid Image	General Methods
            # Image of clone fragment	Image of cloned plasmid with RE sites marked	Gel images and restriction analysis	Discussion: sequencing vs RE	Discussion: Mol weights	Domain organisation	Alignment	Table of Sequences
            # Further Work, Extra Analyses
            # Clear figures, legends, subheadings, cross-refs etc.
            # FINAL MARK (UNROUNDED)	FINAL MARK
            proj = db.addTable('../biol2013.labreports.tdt',['Username'],name='mark_summary')
            sectionheads = ['Student','Introduction','Methods','Results and Discussion','Understanding/Extras','Presentation','TOTAL']
            sections = {'Student':proj.fields()[:4],'Introduction':proj.fields()[4:7],'Methods':proj.fields()[7:10],
                        'Results and Discussion':proj.fields()[10:-3],'Understanding/Extras':proj.fields()[-3:-2],'Presentation':proj.fields()[-2:-1],'TOTAL':proj.fields()[-1:]}
            for shead in sectionheads[1:-1]: proj.addField(shead,evalue=0)
            for entry in proj.entries():
                for shead in sectionheads[1:-1]: 
                    for field in sections[shead]:
                        try: entry[shead] += float(entry[field])
                        except: pass
            proj.renameField('FINAL MARK','Lab Report Total')
            marks = rje.combineDict(proj.data().pop('-'),{'File':10,'Section1':1,'Section2':9,'Section3':10,'Section4':20,'Section5':20,'Section6':30,'Data File Total':100})
            ## ~ [0b] ~ Data File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dat = db.addTable('../BIOL2013.marks.tdt',['Student'],name='data_file')
            dat.renameField('Total','Data File Total')
            for field in ['Student Name','Student Number','Section0','Section1','Section2','Section3','Section4','Section5','Section6','Data File Total']: proj.addField(field,evalue=0)
            for entry in proj.entries():
                stud = entry['Username']
                data = dat.data(stud)
                for field in ['Student Name','Student Number','Section0','Section1','Section2','Section3','Section4','Section5','Section6','Data File Total']: entry[field] = data[field]
            proj.renameField('Section0','File')
            ### ~ [1] ~ Generate HTML of Summary Lab Report and Data File Marks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = '../mark_summary.htm'
            rje.backup(self,hfile)
            html = [rje_html.htmlHead('BIOL2013 Coursework Mark Summary',self.list['StyleSheets'],tabber=False,frontpage=True,nobots=True,keywords=[],javascript='')]
            for entry in proj.entries():
                hi = len(html)
                html.append('<h1>%s</h1>' % (entry['Student Number']))
                html.append('<h2>%s (%s)</h2>' % (entry['Student Name'],entry['Username']))
                ### Main content ###
                html.append('<table border=1>')
                html.append('<tr><th colspan=2>Lab Report (15%)</th></tr>')
                for section in sectionheads[1:-1]:
                    for field in sections[section]: html.append('<tr><td>%s</td><td width="50" align="center">%s / %s</td></tr>' % (field,entry[field],marks[field]))
                    html.append('<tr><th><b>%s</b></th><th width="100" align="center">%s / %s</th></tr>' % (section,entry[section],int(marks[section])))
                section = 'Lab Report Total'
                html.append('<tr><th><b>%s:</b></th><th width="100" align="center">%s / %s</th></tr>' % (section,entry[section],int(marks[section])))
                html += ['</table>']
                student_summary[entry['Username']] = string.join(html[hi:],'\n')
                #for field in ['Introduction','Methods','Results and Discussion','Understanding/Extras','Presentation','Lab Report Total']:
                #    html.append('<tr><td width="180"><b>%s:</b></td><td width="50" align="center">%s</td></tr>' % (field,entry[field]))
                html.append('<br><table border=1>')
                html.append('<tr><th colspan=2>Data File (10%)</th></tr>')
                for field in ['File','Section1','Section2','Section3','Section4','Section5','Section6']:
                    html.append('<tr><td width="180"><b>%s:</b></td><td width="100" align="center">%s / %s</td></tr>' % (field,entry[field],int(marks[field])))
                for field in ['Data File Total']:
                    html.append('<tr><td width="180"><b>%s:</b></td><td width="100" align="center"><b>%s / %s</b></td></tr>' % (field,entry[field],int(marks[field])))
                html.append('</table>')
                html.append('<br><table border=1>')
                if entry['Data File Total'] == 'ERROR':
                    entry['Data File Total'] = 0.0
                    for field in ['File','Section1','Section2','Section3','Section4','Section5','Section6']:
                        entry['Data File Total'] += float(entry[field])
                #self.deBug(entry)
                if not entry['Data File Total']: entry['Data File Total'] = 0.0
                if not entry['Lab Report Total']: entry['Lab Report Total'] = 0.0
                html.append('<tr><th width="180">TOTAL (25%%)</th><th width="100">%d%%</th></tr>' % rje.dp((0.4 * float(entry['Data File Total'])) + (0.6 * float(entry['Lab Report Total'])),0))
                html.append('</table>')
                html.append('<p><font size=-1 family="Georgia"><i>NB. Final marks may be subject to additional moderation.</i></font></p>')
                html.append('<hr>')

            html += [rje_html.htmlTail(tabber=False)]
            if save:
                open(hfile,'w').write(string.join(html,'\n'))
                self.printLog('#HTML','HTML output: %s' % hfile)
            
            return student_summary
        except: self.errorLog('%s.summary2013 error' % self)
#########################################################################################################################
    def OLDsummary2013(self):  ### Generates HTML of Summary Lab Report and Data File Marks
        '''Generates HTML of Summary Lab Report and Data File Marks.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            # Username	First Name	Last Name
            # Title	Contents	Abbrev	Rationale
            # Table of features	Plasmid Image	General Methods
            # Image of clone fragment	Image of cloned plasmid with RE sites marked	Gel images and restriction analysis	Discussion: sequencing vs RE	Discussion: Mol weights	Domain organisation	Alignment	Table of Sequences
            # Further Work, Extra Analyses
            # Clear figures, legends, subheadings, cross-refs etc.
            # FINAL MARK (UNROUNDED)	FINAL MARK
            proj = db.addTable('biol2013.labreports.tdt',['Username'],name='mark_summary')
            sectionheads = ['Student','Introduction','Methods','Results and Discussion','Understanding/Extras','Presentation','TOTAL']
            sections = {'Student':proj.fields()[:3],'Introduction':proj.fields()[3:6],'Methods':proj.fields()[6:9],
                        'Results and Discussion':proj.fields()[9:-4],'Understanding/Extras':proj.fields()[-4:-3],'Presentation':proj.fields()[-3:-2],'TOTAL':proj.fields()[-1:]}
            for shead in sectionheads[1:-1]: proj.addField(shead,evalue=0)
            for entry in proj.entries():
                for shead in sectionheads[1:-1]: 
                    for field in sections[shead]:
                        try: entry[shead] += float(entry[field])
                        except: pass
            proj.renameField('FINAL MARK','Lab Report Total')
            ## ~ [0b] ~ Data File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dat = db.addTable('mark_full.tdt',['Student'],name='data_file')
            dat.renameField('Total','Data File Total')
            for field in ['Student Name','Student Number','Section0','Section1','Section2','Section3','Section4','Section5','Section6','Data File Total']: proj.addField(field,evalue=0)
            for entry in proj.entries():
                stud = entry['Username']
                data = dat.data(stud)
                for field in ['Student Name','Student Number','Section0','Section1','Section2','Section3','Section4','Section5','Section6','Data File Total']: entry[field] = data[field]
            proj.renameField('Section0','File')
            ### ~ [1] ~ Generate HTML of Summary Lab Report and Data File Marks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hfile = 'mark_summary.htm'
            rje.backup(self,hfile)
            html = [rje_html.htmlHead('BIOL2013 Coursework Mark Summary',self.list['StyleSheets'],tabber=False,frontpage=True,nobots=True,keywords=[],javascript='')]
            for entry in proj.entries():
                html.append('<h1>%s</h1>' % (entry['Student Number']))
                html.append('<h2>%s (%s)</h2>' % (entry['Student Name'],entry['Username']))
                ### Main content ###
                html.append('<table border=1>')
                html.append('<tr><th colspan=2>Data File (10%)</th></tr>')
                for field in ['File','Section1','Section2','Section3','Section4','Section5','Section6','Data File Total']:
                    html.append('<tr><td width="180"><b>%s:</b></td><td width="50" align="center">%s</td></tr>' % (field,entry[field]))
                html.append('<tr><td></td><td></td></tr>')
                html.append('<tr><th colspan=2>Lab Report (15%)</th></tr>')
                for field in ['Introduction','Methods','Results and Discussion','Understanding/Extras','Presentation','Lab Report Total']:
                    html.append('<tr><td width="180"><b>%s:</b></td><td width="50" align="center">%s</td></tr>' % (field,entry[field]))
                html.append('<tr><td></td><td></td></tr>')
                if not entry['Data File Total']: entry['Data File Total'] = 0.0
                if not entry['Lab Report Total']: entry['Lab Report Total'] = 0.0
                html.append('<tr><th>TOTAL (25%%)</th><th>%d</th></tr>' % rje.dp((0.4 * float(entry['Data File Total'])) + (0.6 * float(entry['Lab Report Total'])),0))
                html.append('</table><br><hr><br>')
                html.append('<font size=-1 family="Georgia"><i>NB. Final marks may differ from draft Blackboard results following additional checks. Marks will <b>not</b> have gone down.</i></font>')

            html += [rje_html.htmlTail(tabber=False)]
            open(hfile,'w').write(string.join(html,'\n'))
            self.printLog('#HTML','HTML output: %s' % hfile)
            

        except: self.errorLog('%s.summary2013 error' % self)
#########################################################################################################################
    def extract2013(self,logdata=True):  ### Extracts data from BIOL2013 data file
        '''Extracts data from BIOL2013 data file.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edatakeys = ['Student Name','Student Number','Student Username',
                         'Species','Common Name','Ensembl Transcript ID','Ensembl Gene ID','Gene Name (MGI Symbol)','Gene Description','Transcript Length (bp)','Protein Length (aa)','No. Exons',
                         'EnsFasta',
                         '5\' UTR Length','ORF Start','ORF End','Reading Frame','Annotated Context','Annotated Strength','AIC position (start of codon)','AIC residue (amino acid) position','AIC Context','AIC Strength',
                         'Fusion ORF data','400nt from CMV30 primer','800nt from the CMV24 primer',
                         'Annotated fusion protein sequence','Annotated fusion pI','Annotated fusion Mw','Alternative fusion protein sequence','Alternative fusion pI','Alternative fusion Mw','Truncation length difference (aa)','Mw Difference (kDa 3sf)','Alignment'
                         ]
            edata = {}
            for ekey in edatakeys: edata[ekey] = ''
            ### ~ [2] Read in Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FILE','Extract data from %s' % self.getStr('InFile'))
            section = '0'; multiline = None
            for eline in self.loadFromFile(self.getStr('InFile'),chomplines=True):
                if not eline: continue
                ## ~ [2a] ~ Mutli-line data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if multiline and eline[:1] != '#':
                    if eline[:2] == 'pI': multiline = None
                    else:
                        if multiline == 'Alignment':
                            if '>' not in eline: eline = string.join(string.split(eline),'')
                            edata[multiline] = '%s\n%s' % (edata[multiline],eline)
                            #self.deBug('%s' % [eline])
                            #self.deBug([edata[multiline]])
                        else:
                            edata[multiline] = '%s\n%s' % (edata[multiline],eline)
                        continue
                elif multiline:
                    multiline = None    # Reached end of multline
                    if '\n' in edata[ekey]: self.printLog('#DATA','%s = %s...%s' % (ekey,string.split(edata[ekey],'\n')[0][:10],string.split(edata[ekey],'\n')[-1][-10:]))
                    continue
                ## ~ [2b] ~ Simple Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if string.split(eline,':')[0] in edata:
                    ekey = string.split(eline,':')[0]
                    if edata[ekey]: self.errorLog('Data for %s read twice!' % ekey,printerror=False)
                    edata[ekey] = string.join(string.split(string.join(string.split(eline,':')[1:],':')))
                    if logdata: self.printLog('#DATA','%s = %s' % (ekey,edata[ekey]))
                    continue
                ## ~ [2c] ~ File markers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if rje.matchExp('^### SECTION (\d+):',eline):
                    section = rje.matchExp('^### SECTION (\d+):',eline)[0]
                    continue
                if rje.matchExp('### END SECTION (\d+):',eline):
                    section = '%sx' % rje.matchExp('### END SECTION (\d+):',eline)[0]
                    continue
                ## ~ [2d] ~ Special data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if eline[:1] == '>' and section == '3':
                    edata['EnsFasta'] = eline
                    multiline = 'EnsFasta'
                    continue
                if eline[:1] != '#' and section == '5':
                    edata['Fusion ORF data'] = eline
                    multiline = 'Fusion ORF data'
                    section = '5a'
                    continue
                if eline[:1] != '#' and section == '5a':
                    edata['400nt from CMV30 primer'] = eline
                    multiline = '400nt from CMV30 primer'
                    section = '5b'
                    continue
                if eline[:1] != '#' and section == '5b':
                    edata['800nt from the CMV24 primer'] = eline
                    multiline = '800nt from the CMV24 primer'
                    section = '5c'
                    continue
                if eline[:1] != '#' and section == '6':
                    if eline[:2] == 'pI': edata['Annotated fusion pI'] = string.join(string.split(string.join(string.split(eline,':')[1:],':')))
                    elif eline[:2] == 'Mw': edata['Annotated fusion Mw'] = string.join(string.split(string.join(string.split(eline,':')[1:],':'))); section = '6a'
                    else: 
                        multiline = 'Annotated fusion protein sequence'
                        edata[multiline] = eline
                    continue
                if eline[:1] != '#' and section == '6a':
                    if eline[:2] == 'pI': edata['Alternative fusion pI'] = string.join(string.split(string.join(string.split(eline,':')[1:],':')))
                    elif eline[:2] == 'Mw': edata['Alternative fusion Mw'] = string.join(string.split(string.join(string.split(eline,':')[1:],':'))); section = '6b'
                    else: 
                        multiline = 'Alternative fusion protein sequence'
                        edata[multiline] = eline
                    continue
                if eline[:1] != '#' and section == '6b':
                    multiline = 'Alignment'
                    edata[multiline] = eline
                    section = '6c'
                    continue
            ### ~ [3] ~ Strip spaces from some data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(edata)
            for ekey in ['400nt from CMV30 primer','800nt from the CMV24 primer','Annotated fusion protein sequence','Alternative fusion protein sequence']: #,'Alignment']:
                edata[ekey] = string.join(string.split(edata[ekey]),'')
            return (edatakeys,edata)
        except: self.errorLog('%s.extract2013 error' % self,quitchoice=False); return (edatakeys,{})
#########################################################################################################################
### End of SECTION II: Teaching Class                                                                                   #
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
    try: Teaching(mainlog,cmd_list).run()

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
