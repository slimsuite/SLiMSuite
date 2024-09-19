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
Module:       rje_mascot
Description:  Module for loading, manipulating and saving MASCOT data for BUDAPEST and PICSI
Version:      1.2
Last Edit:    27/02/13
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The main function of this module is to read in CSV-exported data from MASCOT and then save it in formats that can be
    more readily utilised and interrogated by other programs.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    mascot=FILE     : Name of MASCOT csv file [None]
    itraq=T/F       : Whether data is from an iTRAQ experiment [False]
    empai=T/F       : Whether emPAI data is present in MASCOT file [True]
    samples=LIST    : List of X:Y, where X is an iTRAQ isotag and Y is a sample []
    
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
    # 0.0 - Initial Compilation based on BUDAPEST 1.9 for use with BUDAPEST 2.0.
    # 1.0 - Working version including reading in of iTRAQ data.
    # 1.1 - Fixed bugs for reading in data with unmatched peptides and iTRAQ data.
    # 1.2 - Added 
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add reading in of non-quantitative MASCOT data.
    # [Y] : Add reading in of iTRAQ data.
    # [Y] : Reading in of unmatched peptide data.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_MASCOT', '1.2', 'February 2013', '2012')
    description = 'Module for loading, manipulating and saving MASCOT data for BUDAPEST and PICSI'
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
### SECTION II: MASCOT Class                                                                                            #
#########################################################################################################################
class MASCOT(rje_obj.RJE_Object):     
    '''
    MASCOT Class. Author: Rich Edwards (2012).

    Str:str
    - MASCOT = Name of input MASCOT file
    
    Bool:boolean
    - emPAI = Whether emPAI data is present in MASCOT file [True]
    - iTRAQ = Whether data is from an iTRAQ experiment [False]

    Int:integer

    Num:float
    
    List:list

    Dict:dictionary    
    - Samples = Dictionary of Tag:Sample {}

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['MASCOT']
        self.boollist = ['emPAI','iTRAQ']
        self.intlist = []
        self.numlist = []
        self.listlist = []
        self.dictlist = ['Samples']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({})
        self.setBool({'emPAI':True})
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
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['MASCOT'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['emPAI','iTRAQ'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['Samples']) # Splits comma separated X:Y pairs into dictionary
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
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['delimit=,'])
            self.splitMascot()
            if self.getBool('iTRAQ') and self.dict['Samples']: self.iTRAQSamples()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### General MASCOT Reading Methods                                                                          #
#########################################################################################################################
    def splitMascot(self):  ### Reads the MASCOT file and splits into header, hits and unmatched files.
        '''Reads the MASCOT file and splits into header, hits and unmatched files.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            infile = self.getStr('MASCOT')
            if self.basefile().lower() in ['','none']: self.basefile(rje.baseFile(self.getStr('MASCOT')))
            #x#self.deBug(self.basefile())
            headfile = '%s.header.txt' % self.basefile()
            hitsfile = '%s.mascot.csv' % self.basefile()
            peptfile = '%s.nohits.csv' % self.basefile()
            if rje.isYounger(self.getStr('MASCOT'),hitsfile) == hitsfile and not self.force():
                return self.printLog('#FILE','%s file found (force=F)' % hitsfile)
            ### ~ [1] Split MASCOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headlines = []
            csvhead = []
            mdb = None
            mx = 0
            itraq = []
            prot_data = {}
            for mline in open(self.getStr('MASCOT'),'r').readlines():
                mx += 1     # Index of next line in case needed for iTRAQ reading!
                ## ~ [1a] Skip down until Header found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not headlines and mline.find('Header') < 0: continue
                ## ~ [1b] Add Header lines to headlines until results headers found ~~~~~~~~~~~~~~~ ##
                if not csvhead and mline.find('prot_hit_num') < 0: headlines.append(mline); continue
                ## ~ [1c] Sort out MASCOT results headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if mline.find('prot_hit_num') >= 0:
                    ## ~ Read Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    open(headfile,'w').writelines(headlines)
                    csvhead = rje.readDelimit(rje.join(rje.split(rje.chomp(mline))),',')
                    while '' in csvhead: csvhead.remove('')
                    ## ~ Sort out iTRAQ headers (missing) ~~~~~~~~~ ##
                    if self.getBool('iTRAQ'):
                        iline = open(self.getStr('MASCOT'),'r').readlines()[mx]
                        for isplit in rje.readDelimit(iline,',')[len(csvhead):]:  # Should be start of iTRAQ data
                            if '/' in isplit: itraq.append(isplit)
                        self.printLog('#ITRAQ',rje.join(itraq))
                        csvhead += itraq
                        idb = db.addEmptyTable('itraq',['prot_hit_num','prot_acc','prot_desc','itraq','ratio','n','geomean','summary'],keys=['prot_hit_num','itraq'])
                        idb.info['Delimit'] = ','
                    ## ~ Add emPAI header (also missing) ~~~~~~~~~~ ##
                    if self.getBool('emPAI'): csvhead.append('empai')
                    ## ~ Set up Database Table ~~~~~~~~~~~~~~~~~~~~ ##
                    self.printLog('#HEAD',rje.join(csvhead,'; '))
                    mdb = db.addEmptyTable('mascot',csvhead,keys=['prot_hit_num','pep_query'])
                    mdb.info['Delimit'] = ','
                elif mline.find('Peptide matches') >= 0:
                    mdb.saveToFile()
                    if self.getBool('emPAI'): csvhead.remove('empai')
                    mdb = db.addEmptyTable('nohits',csvhead,keys=['pep_query'])
                    for field in mdb.fields():
                        if field[:4] == 'prot': mdb.dropField(field)
                    mdb.info['Delimit'] = ','
                    continue
                elif rje.chomp(mline):
                    #self.deBug('%s ... %s' % (mline[:20],mline.find('Peptide matches')))
                    data = rje.readDelimit(mline,',')
                    entry = {}; pretraq = True
                    #self.deBug(csvhead); self.deBug(itraq);
                    for d in range(len(csvhead)+len(itraq)):
                        if d >= len(data): break
                        if data[d] in itraq: dhead = data[d]; pretraq = False
                        elif data[d] == 'emPAI': entry['empai'] = data[d+1]; pretraq = False
                        elif pretraq and d < len(csvhead): dhead = csvhead[d]
                        elif pretraq: continue      # Unmatched peptides will not have emPAI or iTRAQ data
                        #self.deBug('%s > %s' % (data[d],dhead))
                        if d and data[d-1] == 'emPAI': continue
                        elif data[d] in itraq + ['emPAI']: continue
                        elif dhead not in entry: entry[dhead] = data[d]
                        #self.deBug('%s = %s' % (dhead,entry[dhead]))
                    if entry['prot_acc']: prot_data[entry['prot_hit_num']] = {'prot_acc':entry['prot_acc'],'prot_desc':entry['prot_desc']}
                    if self.getBool('iTRAQ') and 'Quantitation summary for protein' in data:
                        d = data.index('Quantitation summary for protein') + 1
                        if entry['prot_hit_num'] in prot_data:
                            pacc = prot_data[entry['prot_hit_num']]['prot_acc']
                            pdesc = prot_data[entry['prot_hit_num']]['prot_desc']
                        else:
                            pacc = entry['prot_acc']
                            pdesc = entry['prot_desc']
                        while d < len(data):
                            if data[d] in itraq:
                                idb.addEntry({'prot_hit_num':entry['prot_hit_num'],'prot_acc':pacc,'prot_desc':pdesc,
                                              'itraq':data[d],'ratio':data[d+1],'n':data[d+2],'geomean':data[d+3],'summary':data[d+4]})
                            d += 1
                    #self.deBug(entry)
                    if entry['prot_hit_num'] or entry['pep_query']: mdb.addEntry(entry)
            mdb.saveToFile()
            if self.getBool('iTRAQ'): idb.saveToFile()
            self.deBug('')
            return True
        except: self.errorLog('Error reading MASCOT file'); return False
#########################################################################################################################
    ### <4> ### iTRAQ Summary Methods                                                                                   #
#########################################################################################################################
    def iTRAQSamples(self): ### Uses self.dict['Samples'] and self.db('itraq') to summarise hit data
        '''Uses self.dict['Samples'] and self.db('itraq') to summarise hit data.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db(); idb = self.db('itraq')
            mdb = db.copyTable(idb,'itraq_summary')
            gdb = db.copyTable(idb,'itraq_geomean')
            ### ~ [1] Reformat Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.dropField('geomean'); gdb.dropField('ratio'); gdb.renameField('geomean','ratio')
            for sdb in [mdb,gdb]:
                sdb.dropField('summary');
                sdb.dropEntriesDirect('ratio','---')
                sdb.dropEntriesDirect('ratio','NN')
                sdb.dataFormat({'ratio':'num','n':'int'})
                ## ~ [1a] Drop tags with Samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (ex,etot) = (0.0,sdb.entryNum())
                for entry in sdb.entries():
                    self.progLog('\r#ITRAQ','Drop isotags without Sample info: %.2f%%' % (ex/etot)); ex += 100.0
                    tags = rje.split(entry['itraq'],'/')
                    if tags[0] not in self.dict['Samples'] or tags[1] not in self.dict['Samples']: sdb.dropEntry(entry)
                self.printLog('\r#ITRAQ','Dropped all isotags without Sample info: %s of %s entries remain' % (rje.iStr(sdb.entryNum()),rje.iStr(etot)))
                ## ~ [1b] Reshape, rename, invert and remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sdb.reshapeWide('itraq',['ratio','n'])
                samples = rje.sortUnique(self.dict['Samples'].values())
                ratios = []
                self.printLog('#SAMP',rje.join(samples,', '))
                for s1 in samples:
                    for s2 in samples[samples.index(s1):]:
                        newfield = '%s/%s' % (s1,s2)
                        sdb.addField(newfield)
                        sdb.addField('%s_Min' % newfield)
                        sdb.addField('%s_Max' % newfield)
                        sdb.addField('%s_Dirn' % newfield)
                        ratios.append(newfield)
                        for entry in sdb.entries(): entry[newfield] = []
                for field in sdb.fields():
                    if '|' in field:
                        (score,tags) = rje.split(field,'|')
                        tag = rje.split(tags,'/')
                        if int(tag[0]) > int(tag[1]):   ### Invert
                            newfield = '%s|%s/%s' % (score,tag[1],tag[0])
                            if newfield in sdb.fields(): sdb.dropField(newfield); continue
                            sdb.renameField(field,newfield)
                            if score == 'ratio':
                                for entry in sdb.entries():
                                    if entry[newfield]: entry[newfield] = 1.0 / entry[newfield]
                            tag = (tag[1],tag[0])
                            field = newfield
                        s1 = self.dict['Samples'][tag[0]]
                        s2 = self.dict['Samples'][tag[1]]
                        newname = '%s|%s%s/%s%s' % (score,s1,tag[0],s2,tag[1])
                        sdb.renameField(field,newname)
                        if score == 'n': continue
                        newfield = '%s/%s' % (s1,s2)
                        invfield = '%s/%s' % (s2,s1)
                        for entry in sdb.entries():
                            if entry[newname] and newfield in sdb.fields(): entry[newfield].append(entry[newname])
                            elif entry[newname]: entry[invfield].append(1.0/entry[newname])
                ## ~ [1c] Calculate Geometric mean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (ex,etot) = (0.0,sdb.entryNum())
                for entry in sdb.entries():
                    self.progLog('\r#GEO','Calculating Geometric means: %.2f%%' % (ex/etot)); ex += 100.0
                    for ratio in ratios:
                        if entry[ratio]:
                            entry['%s_Min' % ratio] = min(entry[ratio])
                            entry['%s_Max' % ratio] = max(entry[ratio])
                            try: entry[ratio] = rje.geoMean(entry[ratio])
                            except: self.deBug(entry)
                            if entry[ratio] > 1 and entry['%s_Min' % ratio] > 1: entry['%s_Dirn' % ratio] = 'UP'
                            elif entry[ratio] < 1 and entry['%s_Max' % ratio] < 1: entry['%s_Dirn' % ratio] = 'DOWN'
                        else: entry['%s_Dirn' % ratio] = entry['%s_Min' % ratio] = entry['%s_Max' % ratio] = entry[ratio] = ''
                self.printLog('\r#GEO','Geometric mean calculations complete')
                sdb.saveToFile()
        except: self.errorLog('iTRAQSamples error')
#########################################################################################################################
### End of SECTION II: MASCOT Class                                                                                     #
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
    try: MASCOT(mainlog,cmd_list).run()

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
