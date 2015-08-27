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
Module:       ProDigIS
Description:  Protein Digestion In Silico
Version:      0.2
Last Edit:    21/06/11
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module is to take one or more lists of proteins, perform in silico digestions using different
    proteases (and combinations) and then output some predicted stats about peptide numbers etc. It is hoped that this
    data can be used in time to predict which proteins can be potentially identified by Mass Spec and, for a given
    proteome, which combination of proteases should maximise coverage.

Commandline:
    ## Basic Input Parameters ##
    seqfiles=LIST   : Sequence files for input. Wildcards permitted. RJE_SEQ filtering WILL be applied. [*.fas]
    source=FILE     : File containing source protein data (including UniProt AccNum column) [None]
    proteases=LIST  : List of proteases to use [tryp,aspn,gluc,lysc]
    peptides=FILE   : File containing list of identified peptides [None]
    pepcut=X        : Protease used to generate list of identified peptides [tryp]
    positives=FILE  : File of positively identified proteins matching peptide lists [None]

    ## Basic Processing Parameters ##
    combcut=T/F     : Whether to peform combined digestions with pairs of proteases [True]
    nterm=T/F       : Whether to include N-terminal peptides [False]
    nrpep=T/F       : Whether to only include the non-redundant (unique) peptides [False]
    cysweight=T/F   : Whether to weight peptide probabilities according to cysteine count [True]

    ## Basic Output Parameters ##
    minpeplen=X     : Minimum peptide length to consider [0]
    maxpeplen=X     : Maximum peptide length to individually return [40]
    pepmwt=T/F      : Whether to output peptide mol weights in addition to lengths [True]
    cyscount=T/F    : Whether to perform peptide count with Cysteine numbers [True]

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
import rje, rje_db, rje_obj, rje_seq, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added probability calculations based on hydrophobicity, serine and cysteine.
    # 0.2 - Added cysteine count and cysteine weighting.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('ProDigIS', '0.2', 'June 2011', '2011')
    description = 'Protein Digestion In Silico'
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
proteases = {'tryp':['K:','R:'],
             'aspn':[':N'], # AspN (cleaves N-terminal side of Asp bonds)
             'gluc':['E:'], # GluC (cleaves C-terminal side of Glu)
             'lysc':['K:'], # Lys C (cuts the C-terminal side of Lys)
            }
bad_combo = [('lysc','tryp')]   # Combinations of proteases that are not worth trying
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: ProDigIS Class                                                                                          #
#########################################################################################################################
class ProDigIS(rje_obj.RJE_Object):     
    '''
    ProDigIS Class. Author: Rich Edwards (2011).

    Str:str
    - PepCut = Protease used to generate list of identified peptides [tryp]
    - Peptides = File containing list of identified peptides [None]
    - Positives = File of positively identified proteins matching peptide lists [None]
    - Source = File containing source protein data (including UniProt AccNum column) [None]
    
    Bool:boolean
    - CombCut = Whether to peform combined digestions with pairs of proteases [True]
    - CysWeight = Whether to weight peptide probabilities according to cysteine count [True]
    - CysCount = Whether to perform peptide count with Cysteine numbers [True]
    - NRPep = Whether to only include the non-redundant (unique) peptides [False]
    - NTerm = Whether to include N-terminal peptides [False]
    - PepMWt = Whether to output peptide mol weights in addition to lengths [True]

    Int:integer
    - MinPepLen = Minimum peptide length to consider [0]
    - MaxPepLen = Maximum peptide length to individually return [41]

    Num:float
    
    List:list
    - SeqFiles = Sequence files for input. Wildcards permitted. RJE_SEQ filtering WILL be applied. [*.fas]
    - NegPep = List of missing peptides from positive protein set
    - Peptides = List of peptides identified by MS
    - PosPep = List of positively identified peptides from positive protein set
    - Proteases = List of proteases to use [tryp,aspn,gluc,lysc]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db object.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['PepCut','Peptides','Positives','Source']
        self.boollist = ['CombCut','NRPep','NTerm','PepMWt','CysWeight','CysCount']
        self.intlist = ['MaxPepLen','MinPepLen']
        self.numlist = []
        self.listlist = ['Proteases','SeqFiles']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'PepCut':'tryp'})
        self.setBool({'CombCut':True,'NTerm':False,'PepMWt':True,'NRPep':False,'CysWeight':True,'CysCount':True})
        self.setInt({'MaxPepLen':41})
        self.setNum({})
        self._cmdRead('seqfiles=*.fas',type='glist',att='SeqFiles')  # Default = all *.fas
        self.list['Proteases'] = string.split('tryp,aspn,gluc,lysc',',')
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
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'file',['Peptides','Positives','Source'])
                self._cmdReadList(cmd,'bool',['CombCut','NRPep','NTerm','PepMWt','CysWeight','CysCount'])  
                self._cmdReadList(cmd.lower(),'list',['Proteases'])  
                self._cmdReadList(cmd,'glist',['SeqFiles'])
                self._cmdReadList(cmd,'int',['MaxPepLen','MinPepLen'])
                self._cmdReadList(cmd,'info',['PepCut'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self._digest()
            self._output()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db = self.db().addEmptyTable('ProDigIS',['AccNum','Protease','PepCount'],['AccNum','Protease'])
            if self.getInt('MinPepLen') > 0: db.addField('MinPepLen')
            if self.getBool('NRPep'): db.addField('NRPep')
            if rje.exists(self.getStr('Source')):
                fdb = self.db().addTable(self.getStr('Source'),mainkeys=['AccNum'],name='Source')
                fdb.addField('File')
                fdb.addField('ProtMWt')
            else: fdb = self.db().addEmptyTable('Source',['AccNum','File','ProtMWt'],['AccNum'])
            for i in range(1,self.getInt('MaxPepLen')+1): db.addField(i)
            if self.getBool('PepMWt'):
                for i in range(1,self.getInt('MaxPepLen')+1): db.addField(i*100.0)
            ### ~ [2] Load Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F'])
            self.obj['SeqList'].seq = fullseq = []
            for seqfile in self.list['SeqFiles']:
                file = rje.baseFile(seqfile,True)
                seqlist = rje_seq.SeqList(self.log,['autofilter=T','gnspacc=T','seqnr=F']+self.cmd_list+['seqin=%s' % seqfile,'autoload=T'])
                fullseq += seqlist.seqs()
                for seq in seqlist.seqs():
                    accnum = seq.getStr('AccNum')
                    try:
                        entry = fdb.data()[accnum]
                        if 'File' in entry and entry['File']: self.errorLog('%s found in %s AND %s!' % (accnum,entry['File'],file),printerror=False)
                        entry['File'] = file
                        entry['ProtMWt'] = seq.MWt()
                    except:
                        entry = {'AccNum':accnum,'File':file,'ProtMWt':seq.MWt()}
                        fdb.addEntry(entry)
                    self.deBug(fdb.dict['Data'][seq.getStr('AccNum')])
            self.printLog('#SEQ','%s sequences to analyse in total' % rje.iLen(fullseq))
            fdb.fillBlanks()
            ### ~ [3] Setup Peptide Probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self._peptideProbabilities():
                db.addField('LenExp','PepCount');
                if self.getBool('PepMWt'): db.addField('MWtExp','LenExp'); db.addField('Len7Exp','MWtExp')
                else: db.addField('Len7Exp','LenExp')
                db.addField('Len37','Len7Exp')
                if self.getBool('PepMWt'):
                    db.addField('Len5','MWtExp'); db.addField('MWt5','Len5')
                    db.addField('Len3','MWtExp'); db.addField('MWt3','Len3')
                else: db.addField('Len5','LenExp'); db.addField('Len3','LenExp')
            return
            ### ~ [4] Temp GABLAM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db().addTable('Chlam_Pos.vs.embl_bacteria.hitsum.tdt',['Qry'],name='GABLAM')
            ndb = self.db().addTable('Chlam_Neg.vs.embl_bacteria.hitsum.tdt',['Qry'],name='GNeg')
            self.db().mergeTables(gdb,ndb,overwrite=True,matchfields=True)
            gdb.renameField('Qry','AccNum')
            tmp = self.db().joinTables(name='blast',join=[('Source','AccNum'),('GABLAM','AccNum')],newkey=['AccNum','File'],keeptable=False)
            tmp.saveToFile()
            tmp.compress(['File'],default='mean')
            tmp.dropFields(['AccNum'])
            tmp.info['Name'] = 'blastsum'
            tmp.saveToFile()
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def _positiveAndNegativePeptides(self): ### Populates PosPep and NegPep Lists
        '''Populates PosPep and NegPep Lists.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfile = '%s.peptides.tdt' % self.basefile()
            #if rje.exists(pfile) and not self.getBool('Force'):
            #    try:
            #        pdb = self.db().addTable(pfile,['Peptide'],name='Peptides')
            #        pdb.dataFormat(reformat={'Len':'int','MWt':'num','Cys':'int','Ser':'int','Hyd':'num'})
            #        self.list['Peptides'] = self.list['PosPep'] = pdb.index('Pos')['Y']
            #        self.list['NegPep'] = pdb.index('Positive')['Neg']
            #        return pdb
            #    except: pass
            if not rje.exists(self.getStr('Peptides')) or not rje.exists(self.getStr('Positives')): return False
            self.list['Peptides'] = peplist = self.loadFromFile(self.getStr('Peptides'),chomplines=True)
            seqlist = rje_seq.SeqList(self.log,['autofilter=T','gnspacc=T','seqnr=F']+self.cmd_list+['seqin=%s' % self.getStr('Positives'),'autoload=T'])
            pdb = self.db().addEmptyTable('Peptides',['Peptide','NR','Pos','Len','MWt','C','HPW','DENQ','M','Hyd'],['Peptide'])
            ### ~ [1] ~ Digest Positives ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            protease = self.getStr('PepCut')
            self.list['PosPep'] = poslist = []; self.list['NegPep'] = neglist = []; sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#PEP','Processing positive proteins (%s): %.2f%%' % (protease,sx/stot)); sx += 100.0
                sequence = seq.getSequence()
                for cut in proteases[protease]: sequence = string.join(string.split(sequence,string.replace(cut,':','')),cut)
                frag = string.split(sequence,':')
                while '' in frag: frag.remove('')
                if not self.getBool('NTerm'): frag = frag[1:]
                for pep in frag[0:]:
                    if pep not in poslist: poslist.append(pep)
            self.printLog('\r#PEP','Processed positive proteins (%s): %s peptides' % (protease,rje.iLen(poslist)))
            ## ~ [1b] ~ Peptide Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            allpep = []; self.list['Redundant'] = redundant = []
            sx = 0.0; stot = self.obj['SeqList'].seqNum() 
            for seq in self.obj['SeqList'].seqs():
                self.progLog('\r#DIG','%s Digesting sequences: %.2f%%' % (protease,sx/stot)); sx += 100.0
                sequence = seq.getSequence()
                for cut in proteases[protease]:
                    sequence = string.join(string.split(sequence,string.replace(cut,':','')),cut)
                for frag in string.split(sequence,':'):
                    if frag in allpep: redundant.append(frag)
                    else: allpep.append(frag)
            self.printLog('\r#DIG','%s Digesting %s sequences complete.' % (protease,rje.iStr(stot)))   
            ## ~ [1c] ~ Process fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            px = 0.0; ptot = len(poslist)
            for pep in poslist[0:]:
                self.progLog('\r#PEP','Processing positive peptides (%s): %.2f%%' % (protease,px/ptot)); px += 100.0
                entry = {'Peptide':pep,'MWt':rje_sequence.MWt(pep),'Hyd':rje_sequence.eisenbergHydropathy(pep,returnlist=False),
                         'Len':len(pep),'NR':'Y','Pos':'Y'}
                if pep not in peplist: poslist.remove(pep); neglist.append(pep); entry['Pos'] = 'N'
                if pep in redundant: entry['NR'] = 'N'
                for aacomb in ['C','HPW','DENQ','M']:
                    x = 0
                    for a in aacomb: x += pep.count(a)
                    entry[aacomb] = x
                pdb.addEntry(entry)
            self.printLog('\r#PEP','Processing positive peptides (%s) complete: %s Pos; %s Neg.' % (protease,rje.iLen(poslist),rje.iLen(neglist)))
            ### ~ [2] ~ Save Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb.saveToFile(pfile)
            POS = open('%s.positives.fas' % self.basefile(),'w'); NEG = open('%s.negatives.fas' % self.basefile(),'w')
            for pep in poslist: POS.write('>%s\n%s\n' % (pep,pep))
            for pep in neglist: NEG.write('>%s\n%s\n' % (pep,pep))
            POS.close(); self.printLog('#FAS','%s peptides output to %s.positives.fas' % (rje.iLen(poslist),self.basefile()))
            NEG.close(); self.printLog('#FAS','%s peptides output to %s.negatives.fas' % (rje.iLen(neglist),self.basefile()))
            return pdb
        except: self.errorLog('Problem during %s._positiveAndNegativePeptides().' % self); return None  # Setup failed
#########################################################################################################################
    def _peptideProbabilities(self):    ### Read in peptides and positives and calculate probability of return
        '''Read in peptides and positives and calculate probability of return.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('CysWeight'): return self._cysteinePeptideProbabilities()
            self._positiveAndNegativePeptides()
            #return self.printLog('#NOPROB','Probability calculation temporarily suspended')
            pfile = '%s.pep_prob.tdt' % self.basefile()
            if rje.exists(pfile) and not self.getBool('Force'):
                try:
                    pdb = self.db().addTable(pfile,['PepSize'],name='PepProb')
                    pdb.dataFormat(reformat={'PepSize':'num','Positive':'int','Negative':'int','Prob':'num'})
                    for entry in pdb.entries():
                        if entry['PepSize'] < 100: entry['PepSize'] = int(entry['PepSize'])
                    return pdb
                except: pass
            pdb = self.db().addEmptyTable('PepProb',['PepSize','Positive','Negative','Prob'],['PepSize'])
            if not rje.exists(self.getStr('Peptides')) or not rje.exists(self.getStr('Positives')): return False
            ## ~ [0a] ~ Load Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            peplist = self.loadFromFile(self.getStr('Peptides'),chomplines=True)
            ## ~ [0b] ~ Load Positives ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = rje_seq.SeqList(self.log,['autofilter=T','gnspacc=T','seqnr=F']+self.cmd_list+['seqin=%s' % self.getStr('Positives'),'autoload=T'])
            ### ~ [1] ~ Digest Positives and Update PepProb Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            protease = self.getStr('PepCut')
            ## ~ [1a] ~ Create new database entry to fill with data ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edict = {}
            for i in range(1,self.getInt('MaxPepLen')+1):
                edict[i] = pdb.addEntry({'PepSize':i,'Positive':0,'Negative':0,'Prob':1.0})
                if self.getBool('PepMWt'): edict[i*100.0] = pdb.addEntry({'PepSize':i*100.0,'Positive':0,'Negative':0,'Prob':1.0})
            ## ~ [1b] ~ For each recognition site of each protease, mark cuts with ":" ~~~~~~~~ ##
            poslist = []; neglist = []; sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#PEP','Processing positive proteins (%s): %.2f%%' % (protease,sx/stot)); sx += 100.0
                sequence = seq.getSequence()
                for cut in proteases[protease]: sequence = string.join(string.split(sequence,string.replace(cut,':','')),cut)
                frag = string.split(sequence,':')
                while '' in frag: frag.remove('')
                if not self.getBool('NTerm'): frag = frag[1:]
                for pep in frag[0:]:
                    if self.getBool('NRPep') and pep in self.list['Redundant']: continue
                    if pep not in poslist: poslist.append(pep)
            self.printLog('\r#PEP','Processed positive proteins (%s): %s peptides' % (protease,rje.iLen(poslist)))
            ## ~ [1c] ~ Process fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            px = 0.0; ptot = len(poslist)
            for pep in poslist[0:]:
                self.progLog('\r#PEP','Processing positive peptides (%s): %.2f%%' % (protease,px/ptot)); px += 100.0
                plen = min(len(pep),self.getInt('MaxPepLen'))
                if pep in peplist: edict[plen]['Positive'] += 1
                else: edict[plen]['Negative'] += 1; poslist.remove(pep); neglist.append(pep)
                if self.getBool('PepMWt'):
                    pwt = 100.0 * min(int((rje_sequence.MWt(pep)+99)/100.0),self.getInt('MaxPepLen'))
                    if pep in peplist: edict[pwt]['Positive'] += 1
                    else: edict[pwt]['Negative'] += 1
            self.printLog('\r#PEP','Processing positive peptides (%s) complete.' % protease)
            ## ~ [1d] # Calculate peptide probabilities for protease combo ~~~~~~~~~~~~~~~~~~~~ ##
            for entry in edict.values():
                try: entry['Prob'] = float(entry['Positive']) / float(entry['Positive']+entry['Negative'])
                except: entry['Prob'] = 0.0
            ### ~ [2] ~ Save File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb.saveToFile(pfile)
            return pdb
        except: self.errorLog('Problem during %s._peptideProbabilities().' % self); return None  # Setup failed
#########################################################################################################################
    def _cysteinePeptideProbabilities(self):    ### Read in peptides and positives and calculate probability of return
        '''Read in peptides and positives and calculate probability of return.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pepdb = self._positiveAndNegativePeptides()
            maxcys = max(pepdb.index('C').keys())
            pdb = self.db().addEmptyTable('PepProb',['PepSize','CysCount','Positive','Negative','Prob'],['PepSize','CysCount'])
            edict = {}
            for i in range(1,self.getInt('MaxPepLen')+1):
                edict[i] = {}
                for c in range(-1,maxcys+1):
                    edict[i][c] = pdb.addEntry({'PepSize':i,'CysCount':c,'Positive':0,'Negative':0,'Prob':1.0})
                    if self.getBool('PepMWt'): edict[i*100.0][c] = pdb.addEntry({'PepSize':i*100.0,'CysCount':c,'Positive':0,'Negative':0,'Prob':1.0})
            ## ~ [1b] ~ For each recognition site of each protease, mark cuts with ":" ~~~~~~~~ ##
            ex = 0.0; etot = pepdb.entryNum()
            for entry in pepdb.entries():    #addEmptyTable('Peptides',['Peptide','NR','Pos','Len','MWt','C','HPW','DENQ','M','Hyd'],['Peptide'])
                self.progLog('\r#PEP','Processing peptides probabilities: %.2f%%' % (ex/etot)); ex += 100.0
                pep = entry['Peptide']
                if self.getBool('NRPep') and entry['NR'] == 'N': continue
                plen = min(entry['Len'],self.getInt('MaxPepLen'))
                if entry['Pos'] == 'Y': edict[plen][entry['C']]['Positive'] += 1; edict[plen][-1]['Positive'] += 1
                else: edict[plen][entry['C']]['Negative'] += 1; edict[plen][-1]['Negative'] += 1
                if self.getBool('PepMWt'):
                    pwt = 100.0 * min(int((entry['MWt']+99)/100.0),self.getInt('MaxPepLen'))
                    if entry['Pos'] == 'Y': edict[pwt][entry['C']]['Positive'] += 1
                    else: edict[pwt][entry['C']]['Negative'] += 1
            ## ~ [1d] # Calculate peptide probabilities for protease combo ~~~~~~~~~~~~~~~~~~~~ ##
            for entry in pdb.entries():
                try: entry['Prob'] = float(entry['Positive']) / float(entry['Positive']+entry['Negative'])
                except: entry['Prob'] = 0.0
                self.deBug(entry)
            self.printLog('\r#PEP','Processing peptides probabilities is complete.')
            ### ~ [2] ~ Save File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb.saveToFile()
            return pdb
        except: self.errorLog('Problem during %s._peptideProbabilities().' % self); return None  # Setup failed
#########################################################################################################################
    ### <3> ### Main Digestion Methods                                                                                  #
#########################################################################################################################
    def protCombo(self):    ### Return proteases combinations
        '''Return proteases combinations.'''
        self.list['Proteases'].sort()
        prot_combo = self.list['Proteases'][0:]
        if self.getBool('CombCut'):
            for p1 in self.list['Proteases'][0:]:
                for p2 in self.list['Proteases'][0:]:
                    if self.list['Proteases'].index(p2) <= self.list['Proteases'].index(p1): continue
                    if (p1,p2) not in bad_combo: prot_combo.append('%s+%s' % (p1,p2))
        return prot_combo
#########################################################################################################################
    def _digest(self): ### Main digestion of sequences and population of results database
        '''Main digestion of sequences and population of results database.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db('ProDigIS')
            prot_combo = self.protCombo()
            ## ~ [1] ~ Peptide Probability Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = self.db('PepProb'); pdict = {}
            if pdb:
                if self.getBool('CysWeight'):
                    for plen in pdb.index('PepSize').keys(): pdict[plen] = {}
                    for entry in pdb.entries(): pdict[entry['PepSize']][entry['CysCount']] = entry
                else:
                    for entry in pdb.entries(): pdict[entry['PepSize']] = entry
            ### ~ [2] Process each sequence in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(self.int)
            for prot in prot_combo:
                allpep = []; redundant = []; maxcys = 0
                sx = 0.0; stot = self.obj['SeqList'].seqNum() 
                for seq in self.obj['SeqList'].seqs():
                    self.progLog('\r#DIG','%s Digesting sequences: %.2f%%' % (prot,sx/stot)); sx += 100.0
                    sequence = seq.getSequence()
                    for protease in string.split(prot,'+'):
                        for cut in proteases[protease]:
                            sequence = string.join(string.split(sequence,string.replace(cut,':','')),cut)
                    for frag in string.split(sequence,':'):
                        if frag in allpep: redundant.append(frag)
                        else: allpep.append(frag); maxcys = max(maxcys,frag.count('C'))
                self.printLog('\r#DIG','%s Digesting %s sequences complete.' % (prot,rje.iStr(stot)))
                if self.getBool('CysCount'):
                    for c in range(maxcys+1): db.addField('Cys%d' % c)
            ### ~ [3] Process each sequence in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                sx = 0.0; stot = self.obj['SeqList'].seqNum() 
                for seq in self.obj['SeqList'].seqs():
                    self.progLog('\r#DIG','%s Digesting sequences: %.2f%%' % (prot,sx/stot)); sx += 100.0
                    acc = seq.getStr('AccNum')
                    ## ~ [2a] ~ Create new database entry to fill with data ~~~~~~~~~~~~~~~~~~~~~~~ ##
                    entry = {'AccNum':acc,'Protease':prot}
                    for i in range(1,self.getInt('MaxPepLen')+1):
                        entry[i] = 0
                        if self.getBool('PepMWt'): entry[i*100.0] = 0
                    sequence = seq.getSequence()
                    ## ~ [2b] ~ For each recognition site of each protease, mark cuts with ":" ~~~~ ##
                    for protease in string.split(prot,'+'):
                        for cut in proteases[protease]:
                            sequence = string.join(string.split(sequence,string.replace(cut,':','')),cut)
                    ## ~ [2c] ~ Cut into fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    frag = string.split(sequence,':')
                    while '' in frag: frag.remove('')
                    self.deBug(frag)
                    entry['PepCount'] = len(frag)
                    if not self.getBool('NTerm'): frag = frag[1:]
                    if self.getInt('MinPepLen') > 0: 
                        for pep in frag[0:]:
                            if len(pep) < self.getInt('MinPepLen'): frag.remove(pep)
                    entry['MinPepLen'] = len(frag)
                    if self.getBool('NRPep'):
                        for pep in frag[0:]:
                            if pep in redundant: frag.remove(pep)
                        entry['NRPep'] = len(frag)
                    if self.getBool('CysCount'):
                        for c in range(maxcys+1): entry['Cys%d' % c] = 0
                        for pep in frag: entry['Cys%d' % pep.count('C')] += 1
                    if pdict: entry['LenExp'] = 0.0; entry['MWtExp'] = 0.0; entry['Len7Exp'] = 0.0
                    ## ~ [2d] ~ Process fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for pep in frag[0:]:
                        plen = min(len(pep),self.getInt('MaxPepLen'))
                        self.deBug('"%s" -> %d' % (pep,plen))
                        entry[plen] += 1
                        if pdict:
                            if self.getBool('CysWeight'):
                                try: pprob = pdict[plen][pep.count('C')]['Prob']
                                except: pprob = 0.0
                            else: pprob = pdict[plen]['Prob']
                        if pdict: entry['LenExp'] += pprob
                        if pdict and 7 <= plen: entry['Len7Exp'] += pprob
                        if self.getBool('PepMWt'):
                            pwt = 100.0 * min(int((rje_sequence.MWt(pep)+99)/100.0),self.getInt('MaxPepLen'))
                            entry[pwt] += 1
                            if pdict: entry['MWtExp'] += pprob
                    entry['Len3'] = rje.logPoisson(3,entry['LenExp'],callobj=self)
                    if self.getBool('PepMWt'): entry['MWt3'] = rje.logPoisson(3,entry['MWtExp'],callobj=self)
                    entry['Len5'] = rje.logPoisson(5,entry['LenExp'],callobj=self)
                    if self.getBool('PepMWt'): entry['MWt5'] = rje.logPoisson(5,entry['MWtExp'],callobj=self)
                    entry['Len37'] = rje.logPoisson(3,entry['Len7Exp'],callobj=self)
                    db.addEntry(entry)
                self.printLog('\r#DIG','%s Digesting %s sequences complete.' % (prot,rje.iStr(stot)))
        except: self.errorLog('%s._digest error' % self)
#########################################################################################################################
    def _output(self):  ### Output results to protein and summary tables
        '''Ouput results to protein and summary tables.'''
        try:### ~ [1] Full protein table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tmp = self.db().joinTables(name='protein',join=[('Source','AccNum'),('ProDigIS','AccNum')],newkey=['AccNum','File','Protease'],keeptable=False)
            tmp.saveToFile()
            ### ~ [2] Summary table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tmp.compress(['File','Protease'],default='mean')
            tmp.dropFields(['AccNum'])
            tmp.info['Name'] = 'summary'
            tmp.saveToFile()
        except: self.errorLog('%s._output error' % self)
#########################################################################################################################
### End of SECTION II: ProDigIS Class                                                                                   #
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
    try:ProDigIS(mainlog,['basefile=prodigis']+cmd_list,newstyle=True).run()
        #print rje_zen.Zen().wisdom(), '\n\n *** No standalone functionality! *** \n\n'

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
