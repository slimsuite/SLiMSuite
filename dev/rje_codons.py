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
Module:       rje_codons
Description:  RJE Codon Usage module
Version:      0.1
Last Edit:    16/08/11
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    Codon Usage Bias module

Commandline:
    seqin=FILE  : Sequence file for Codon Usage analysis

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_seq, rje_sequence, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_sequence, rje_obj, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added Entropy-based calculation of CU Bias
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_CODONS', '0.1', 'August 2011', '2011')
    description = 'RJE Codon Usage module'
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class Codons(rje_obj.RJE_Object):     
    '''
    Codons Class. Author: Rich Edwards (2011).

    Str:str
    
    Bool:boolean

    Int:integer

    Num:float
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object
    - SeqList = rje_seq.SeqList object
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
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        #self.setInfo({})
        #self.setBool({})
        #self.setInt({})
        #self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['SeqList'] = rje_seq.SeqList(self.log,['autoload=T','dna=T']+self.cmd_list)
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
                #self._cmdReadList(cmd,'str',['Att'])  # No need for arg if arg = att.lower()
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.expectedCodonUsage()
            self.codonUsageBias()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList']
            if self.getStr('Basefile').lower() in ['','none']:
                self.str['Basefile'] = rje.baseFile(seqlist.getStr('Name'))
                self.obj['DB'].setInfo({'Basefile':self.str['Basefile']})
            ## ~ [1a] Genetic Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb = self.db().addEmptyTable('Code',['Codon','AA'],['Codon'])
            for codon in rje_sequence.genetic_code: cdb.addEntry({'Codon':codon,'AA':rje_sequence.genetic_code[codon]})
            cdb.index('AA')
            ### ~ [2] Calculate Codon Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            codons = rje.sortKeys(rje_sequence.genetic_code)
            db = self.db().addEmptyTable('Codons',['Seq','Len']+codons,['Seq'])
            sx = 0.0; seqx = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#COD','Calculating codon usage: %.2f%%' % (sx/seqx)); sx += 100.0
                entry = rje_sequence.codons(seq.getSequence(),{})
                #self.deBug(entry); self.deBug(entry.values())
                entry['Len'] = sum(entry.values())
                entry['Seq'] = seq.getStr('AccNum')
                db.addEntry(entry)
            self.printLog('\r#COD','Codon usage calculated for %s sequences' % rje.iStr(seqx))
            db.fillBlanks(blank=0,fillempty=True)
            db.saveToFile()
            ### ~ [3] Calculate NT Count Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nt = ['C','A','G','U']
            for i in [1,2,3]:
                for n in ['C','A','G','U']: nt.append('%s|%d' % (n,i))
            ndb = self.db().addEmptyTable('NT',['Seq','Len']+nt,['Seq'])
            sx = 0.0; seqx = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#NT','Calculating NT Counts: %.2f%%' % (sx/seqx)); sx += 100.0
                entry = rje_sequence.aaFreq(rje.replace(seq.getSequence(),'T','U'),{'C':0,'A':0,'G':0,'U':0},False)
                entry['Len'] = sum(entry.values())
                entry['Seq'] = seq.getStr('AccNum')
                centry = db.data(entry['Seq'])
                for i in [1,2,3]:
                    for n in ['C','A','G','U']: entry['%s|%d' % (n,i)] = 0
                for codon in codons:
                    for i in [1,2,3]:
                        n = codon[i-1]
                        entry['%s|%d' % (n,i)] += centry[codon]
                ndb.addEntry(entry)
            self.printLog('\r#NT','NT Counts calculated for %s sequences' % rje.iStr(seqx))
            ndb.saveToFile()
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def expectedCodonUsage(self):     ### Calculate expected codon usage from Frequency data
        '''Calculate expected codon usage from Frequency data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aacode = self.db('Code').index('AA')
            nt = ['C','A','G','U']; codons = rje.sortKeys(rje_sequence.genetic_code)
            cdb = self.db('Codons'); ndb = self.db('NT')
            nsumdb = self.db().copyTable(ndb,'NTPos',replace=True)
            nsumdb.dropField('Len')
            for n in ['C','A','G','U']: nsumdb.renameField(n,'%s|All' % n)
            nsumdb.reshapeLong('Pos',reshape=['C','A','G','U'])
            nsumdb.compress(['Pos'],{'Pos':'str','Seq':'str'},default='sum')
            nsumdb.dropField('Seq'); nsumdb.addField('Total')
            for entry in nsumdb.entries():
                pos = entry.pop('Pos'); entry.pop('Total')
                rje.dictFreq(entry)
                entry['Pos'] = pos
            nsumdb.saveToFile()
            nexentry = nsumdb.data('3')
            fdb = self.db().addEmptyTable('Freq',['Seq','Len']+nt+codons+['Total'],['Seq'])
            edb = self.db().copyTable(cdb,'Expected',replace=True)
            ### ~ [2] Calculate Frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            x = 0.0; etot = cdb.entryNum()
            for oldentry in cdb.entries():
                self.progLog('\r#FREQ','Calculating Frequencies: %.2f%%' % (x/etot)); x += 100.0
                entry = rje.combineDict({},oldentry)
                seq = entry['Seq']; entry['Total'] = entry.pop('Len')
                exentry = edb.data(seq)
                ntentry = rje.combineDict({},ndb.data()[seq])
                ntentry.pop('Seq'); ntentry.pop('Len')
                rje.dictFreq(ntentry)
                ntentry['Len'] = ntentry.pop('Total')
                for aa in aacode:
                    ax = 0.0; ex = 0.0
                    for codon in aacode[aa]:
                        ax += entry[codon]
                        exentry[codon] = nexentry[codon[0]] * nexentry[codon[1]] * nexentry[codon[2]]
                        ex += exentry[codon]
                    for codon in aacode[aa]:
                        if ax: entry[codon] = len(aacode[aa]) * entry[codon] / ax
                        else: entry[codon] = 0.0
                        exentry[codon] = ax * (exentry[codon] / ex)
                fdb.addEntry(rje.combineDict(entry,ntentry))
            self.printLog('\r#Freq','Frequencies calculated for %s entries' % rje.iStr(etot))
            fdb.saveToFile(); edb.saveToFile()
        except: self.errorLog('%s.expectedCodonUsage error' % self)
#########################################################################################################################
    def codonUsageBias(self):   ### Calculate bias in Codon Usage using Entropy-based measure
        '''Calculate bias in Codon Usage using Entropy-based measure.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aacode = self.db('Code').index('AA')
            nt = ['C','A','G','U']; codons = rje.sortKeys(rje_sequence.genetic_code)
            cdb = self.db('Codons'); edb = self.db('Expected')
            ## ~ [1a] Setup bias table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb = self.db().addEmptyTable('Bias',['Seq','Len','Bias','WtBias','AbsBias'],['Seq'])
            ### ~ [2] Calculate Frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            x = 0.0; etot = cdb.entryNum()
            for codentry in cdb.entries():
                self.progLog('\r#BIAS','Calculating Bias: %.2f%%' % (x/etot)); x += 100.0
                expentry = edb.data(codentry['Seq'])
                entry = {'Seq':codentry['Seq'],'Len':codentry['Len'],'Bias':0.0,'WtBias':0.0,'AbsBias':0.0}
                aafreq = {}
                for aa in aacode:
                    aafreq[aa] = 0.0
                    for code in aacode[aa]: aafreq[aa] += codentry[code]
                for aa in aacode:
                    if not aafreq[aa]: continue
                    for code in aacode[aa]: 
                        entry['Bias'] += rje.modulus(codentry[code] - expentry[code]) / aafreq[aa]
                        entry['WtBias'] += rje.modulus(codentry[code] - expentry[code]) / codentry['Len']
                        entry['AbsBias'] += rje.modulus(codentry[code] - (aafreq[aa]/len(aacode[aa]))) / codentry['Len']
                bdb.addEntry(entry)
            self.printLog('\r#BIAS','Codon Usage entropy bias calculated for %s entries' % rje.iStr(etot))
            bdb.saveToFile()
        except: self.errorLog('%s.expectedCodonUsage error' % self)
#########################################################################################################################
    def codonUsageEntropyBias(self):   ### Calculate bias in Codon Usage using Entropy-based measure
        '''Calculate bias in Codon Usage using Entropy-based measure.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aacode = self.db('Code').index('AA')
            nt = ['C','A','G','U']; codons = rje.sortKeys(rje_sequence.genetic_code)
            cdb = self.db('Codons'); edb = self.db('Expected')
            ## ~ [1a] Setup bias table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb = self.db().addEmptyTable('Bias',['Seq','Len','Bias','ExpBias','WtBias','ExpWtBias'],['Seq'])
            ### ~ [2] Calculate Frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            x = 0.0; etot = cdb.entryNum()
            for codentry in cdb.entries():
                self.progLog('\r#BIAS','Calculating Bias: %.2f%%' % (x/etot)); x += 100.0
                expentry = edb.data(codentry['Seq'])
                entry = {'Seq':codentry['Seq'],'Len':codentry['Len'],'Bias':0.0,'ExpBias':0.0,'WtBias':0.0,'ExpWtBias':0.0}
                aafreq = {}
                for aa in aacode:
                    aafreq[aa] = 0.0
                    for code in aacode[aa]: aafreq[aa] += codentry[code]
                rje.dictFreq(aafreq,total=False)
                for aa in aacode:
                    entry['Bias'] += rje.entropyDict(codentry,aacode[aa])
                    entry['ExpBias'] += rje.entropyDict(expentry,aacode[aa])
                    entry['WtBias'] += (aafreq[aa] * rje.entropyDict(codentry,aacode[aa]))
                    entry['ExpWtBias'] += (aafreq[aa] * rje.entropyDict(expentry,aacode[aa]))
                bdb.addEntry(entry)
            self.printLog('\r#BIAS','Codon Usage entropy bias calculated for %s entries' % rje.iStr(etot))
            bdb.saveToFile()
        except: self.errorLog('%s.expectedCodonUsage error' % self)
#########################################################################################################################
### End of SECTION II: Codons Class                                                                                     #
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
    try: Codons(mainlog,cmd_list).run()

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
