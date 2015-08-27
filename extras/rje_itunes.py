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
Module:       rje_itunes
Description:  iTunes PlayList Processor
Version:      0.1
Last Edit:    01/01/14
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module is to read in tab delimited exported iTunes playlists and process data about Albums,
    Artists etc. Files should be name *.YYMMDD.tdt so that the dates can be extracted and used to identify changes.

Commandline:
    itunes=LIST     : List of itunes playlist files *.YYMMDD.tdt [itunes.*.tdt]
    mintracks=X     : Min. number of tracks for averages [1]
    addscore=T/F    : Adds a score of 100x Rating for each track [True]
    tophtml=T/F     : Whether to output top tunes HTML summary [True]
    toplist=X       : Number of entries to feature in HTML top lists [20]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_html, rje_obj, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added Plays/Track, default Album Artist and topHTML() method.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add additional options
    # [Y] : Strip all non-standard characters from names etc.
    # [ ] : Find out what is going wrong with Green Day album names.
    # [ ] : Handle accented letters
    # [ ] : Correct mean Rating calculations to exclude unrated tracks.
    # [ ] : Add number of tracks to Top Albums table.
    # [ ] : Rename "Top 20 Artists" to "Top 20 Tracks": Artists table -> Tracks.
    # [ ] : Allow longer album names. (Or is this an iTunes cutoff?)
    # [ ] : Add table styles (see slimhtml.css) and/or better annotate CSS selection.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_ITUNES', '0.1', 'January 2014', '2011')
    description = 'iTunes PlayList Processor'
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
# iTunes fields and formatting
iformat = {'Name':'str',
           'Artist':'str',
           'Composer':'str',
           'Album':'str',
           'Grouping':'del',
           'Genre':'str',
           'Size':'int',
           'Time':'int',
           'Disc Number':'int',
           'Disc Count':'int',
           'Track Number':'int',
           'Track Count':'int',
           'Year':'int',
           'Date Modified':'str',#'datetime',
           'Date Added':'str',#'datetime',
           'Bit Rate':'int',
           'Sample Rate':'int',
           'Volume Adjustment':'del',
           'Kind':'str',
           'Equaliser':'del',
           'Comments':'del',
           'Plays':'int',
           'Last Played':'str',#'datetime',
           'Skips':'int',
           'Last Skipped':'str',#'datetime',
           'My Rating':'int',
           'Location':'str'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class iTunes(rje_obj.RJE_Object):     
    '''
    iTunes Class. Author: Rich Edwards (2011).

    Str:str
    - TopHTML = HTML to be output if tophtml=T

    Bool:boolean
    - AddScore = Adds a score of 100x Rating for each track [True]
    - TopHTML = Whether to output top tunes HTML summary [True]

    Int:integer
    - MinTracks = Min. number of tracks for averages [1]
    - TopList = Number of entries to feature in HTML top lists [20]

    Num:float
    
    List:list
    - iTunes = List of itunes playlist files *.YYMMDD.tdt [itunes.*.tdt]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Main Databse object for loading and manipulating playlist data 
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['TopHTML']
        self.boollist = ['AddScore','TopHTML']
        self.intlist = ['MinTracks','TopList']
        self.numlist = []
        self.listlist = ['iTunes']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'Basefile':'iTunes','TopHTML':''})
        self.setBool({'AddScore':True,'TopHTML':True})
        self.setInt({'MinTracks':1,'TopList':20})
        self.setNum({})
        self.list['iTunes'] = glob.glob('itunes.*.tdt')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.cmd_list = ['basefile=iTunes'] + self.cmd_list
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
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'glist',['iTunes'])  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'int',['MinTracks'])
                self._cmdReadList(cmd,'bool',['AddScore','TopHTML'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tables = self.db().tables()[0:]
            ## ~ [2a] ~ Calculate Differences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for table1 in tables:
                for table2 in tables[tables.index(table1)+1:]: 
                    self.difference(table1,table2)
            ## ~ [2b] ~ Calculate Averages & Generate HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for table in self.db().tables()[0:]: self.average(table)
            ## ~ [2c] ~ Output HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('TopHTML'):
                html = rje_html.HTML(self.log,self.cmd_list)
                hfile = '%s.html' % self.basefile()
                rje.backup(self,hfile)
                open(hfile,'w').write(html.htmlHead(title=self.basefile(),tabber=False)+self.getStr('TopHTML')+html.htmlTail(False))
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fixfields = ['Location','Name','Artist','Composer','Album']
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            #self.deBug(self.list['iTunes'])
            for ifile in self.list['iTunes']:
                #self.deBug(string.split(open(ifile,'r').readline(),'\t'))
                idb = db.addTable(ifile,mainkeys=['Location'],name=rje.baseFile(ifile,True))
                for field in iformat:
                    if iformat[field] == 'del' and field in idb.fields(): idb.dropField(field)
                idb.dataFormat(iformat)
                idb.addField('Album_Artist','Album')
                idb.addField('Tracks',evalue=1)
                if self.getBool('AddScore'): idb.addField('Score',evalue=0)
                for entry in idb.entries():
                    for field in fixfields:
                        newval = ''
                        for x in entry[field]:
                            if x.isalnum() or x in '\\/: -_()[].~$': newval += x
                        entry[field] = newval
                    entry['Album_Artist'] = entry['Artist']
                    try:
                        for divider in ['\\\\','\\',':','/']:
                            if len(string.split(entry['Location'],divider)) > 2: 
                                entry['Album_Artist'] = string.split(entry['Location'],divider)[-3]
                                break
                    except:
                        self.errorLog('!')
                        self.deBug(entry['Location'])
                    if not entry['Plays']: entry['Plays'] = 0
                    if not entry['Skips']: entry['Skips'] = 0
                    if self.getBool('AddScore'): 
                        if entry['My Rating']: entry['Score'] = (entry['My Rating'] - 60) / 20.0
                idb.remakeKeys()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def average(self,table,groups=['Artist','Album','Album_Artist'],fields=['Plays','My Rating','Skips','Score']):   ### Generates average stats
        '''
        Generates average stats.
        >> table:Table = iTunes database table
        >> groups:list = fields to be focus of compression for summaries
        >> fields:list = fields to be compressed and summarised.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            irules = {'Tracks':'sum','Score':'sum','Plays':'sum','Skips':'sum'}
            for field in fields: 
                if field not in irules: irules[field] = 'mean'
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for group in groups:
                #if group not in table.fields():
                new = db.copyTable(table,'%s.%s' % (table.name(),group))
                newfields = [group,'Tracks'] + fields 
                if group == 'Album': newfields.insert(0,'Album_Artist')
                new.dropFields(newfields+['Location'],inverse=True)
                new.compress(newfields[:newfields.index(group)+1],irules,default='str')
                new.dropField('Location')
                new.dropEntries(['Tracks<%d' % self.getInt('MinTracks')])
                if group == 'Album': new.makeField(formula='Plays/Tracks',fieldname='Plays/Track')
                new.saveToFile('%s.tdt' % new.name())
            ### ~ [3] HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if not self.getBool('TopHTML'): continue
                if group == 'Artist':
                    sortfield = 'Plays'
                    outfields = ['Name','Artist','Album','Plays']
                    out = table
                if group == 'Album':
                    sortfield = 'Plays/Track'
                    outfields = ['Album','Album_Artist','Plays/Track']
                    out = new
                if group == 'Album_Artist':
                    sortfield = 'Plays'
                    outfields = ['Album_Artist','Plays','Tracks']
                    out = new
                html = ['<h2>%s - Top %d %ss</h2>' % (new.name(),self.getInt('TopList'),string.replace(group,'_',' ')),
                        '<table border=1 width="100%"><tr>']
                for field in ['#'] + outfields: html.append('<th>%s</th>' % field)
                html.append('</tr>')
                toplist = []
                lastx = -1      # Last sortfield value
                for ikey in rje.sortKeys(out.index(sortfield),revsort=True):
                    self.debug(ikey)
                    self.debug(out.index(sortfield)[ikey])
                    for entry in out.indexEntries(sortfield,ikey):
                        self.debug(entry)
                        if entry[sortfield] == lastx:
                            html.append('<tr><td>=</td>')
                        elif len(toplist) >= self.getInt('TopList'): break
                        else:
                            html.append('<tr><td>%d</td>' % (len(toplist)+1))
                        toplist.append(entry)
                        for field in outfields: html.append('<td>%s</td>' % entry[field])
                        html.append('</tr>')
                        lastx = entry[sortfield]
                html.append('</table>\n<!-- End of %s Top %d %ss -->\n\n' % (out.name(),self.getInt('TopList'),string.replace(group,'_',' ')))
                self.str['TopHTML'] += string.join(html,'\n')

        except: self.errorLog('%s.average error' % self)
#########################################################################################################################
    def difference(self,table1,table2): ### Generates differences as new table
        '''
        Generates differences as new table.
        >> table1:Table = iTunes database table to compare
        >> table2:Table = iTunes database table to compare
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dfields = ['Name','Artist','Composer','Album','Album_Artist','Genre','Time','Disc Number','Disc Count','Track Number','Track Count','Year','Date Added','Plays','Last Played','Skips','Last Skipped','My Rating','Location','Tracks','Score']
            db = self.db()
            tabindex = '#Artist#|#Album#|#Track Number#|#Name#'
            try:
                age1 = string.atoi(string.split(table1.name(),'.')[-1])
                age2 = string.atoi(string.split(table2.name(),'.')[-1])
                table1.index(tabindex,make=True)
                table2.index(tabindex,make=True)
                if age1 < age2: oldtable = table1; newtable = table2; newdate = age2
                else: newtable = table1; oldtable = table2; newdate = age1              
                diftable = db.copyTable(newtable,'%s-%s' % (oldtable.name(),string.split(newtable.name(),'.')[-1]))
                diftable.keepFields(dfields+[tabindex])
                diftable.addField('Status')
            except: self.errorLog('Cannot generate differences for %s and %s' % (table1,table2))
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#NEW','%s tracks in new iTunes export.' % rje.iStr(newtable.entryNum()))
            self.printLog('#OLD','%s tracks in old iTunes export.' % rje.iStr(oldtable.entryNum()))
            oldfiles = oldtable.datakeys()[0:]
            for entry in diftable.entries():
                ## ~ [2a] Find pair of entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if entry['Location'] in oldfiles: oldentry = oldtable.data(entry['Location'])
                elif entry[tabindex] in oldtable.index(tabindex): 
                    oldentry = oldtable.indexEntries(tabindex,entry[tabindex])[0]
                    if len(oldtable.indexEntries(tabindex,entry[tabindex])) == 1: pass
                    else:
                        self.printLog('#DUP','Duplicate entries for %s' % entry[tabindex])
                        for ientry in oldtable.indexEntries(tabindex,entry[tabindex]):
                            if ientry['Location'] in oldfiles: oldentry = ientry; break
                else: oldentry = None
                #self.deBug(entry)
                #self.deBug(oldentry)
                ## ~ [2b] Generate Differences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not oldentry:
                    entry['Status'] = 'New'
                    continue
                #self.deBug(oldentry['Location'] in oldfiles)
                if oldentry['Location'] in oldfiles: oldfiles.remove(oldentry['Location'])
                #self.deBug(len(oldfiles))
                changed = False
                for field in ['Plays','Skips','My Rating']:
                    if entry[field] != oldentry[field]: changed = True
                    try: entry[field] -= oldentry[field]
                    except: pass    # Keep new value - probably empty in old entry
                if changed: entry['Status'] = 'Changed'
                else: entry['Status'] = 'Unchanged'            
            ### ~ [3] Add missing old entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            reportdel = rje.yesNo('Report deleted %s tracks?' % diftable.name())
            for old in oldfiles:
                entry = diftable.addEntry(oldtable.data(old))
                entry['Status'] = 'Deleted'
                if reportdel: self.printLog('#DEL','%s: %s [%s]' % (entry['Artist'],entry['Name'],entry['Album']))
            ### ~ [4] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for status in rje.sortKeys(diftable.index('Status')):
                self.printLog('#STAT','%s: %d tracks' % (status.upper(),len(diftable.index('Status')[status])))
            self.printLog('#TRACK','%s tracks in total' % rje.iStr(diftable.entryNum()))
            self.deBug('?')
            for table in [table1,table2,diftable]: table.dropField(tabindex)
            diftable.saveToFile('%s.tdt' % diftable.name())
        except: self.errorLog('%s.difference() error' % self)
#########################################################################################################################
### End of SECTION II: iTunes Class                                                                                     #
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
    try: iTunes(mainlog,cmd_list).run()

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
