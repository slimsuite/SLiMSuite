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
Module:       rje_chlamydia
Description:  Bespoke Chlamydia analysis pipeline
Version:      0.2.0
Last Edit:    13/11/15
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    winsize=X   : Window size for NT composition plots [5000]
    winstep=X   : Window step size for NT composition plots [1000]

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
import rje, rje_codons, rje_db, rje_obj, rje_seqlist, rje_zen
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
    (program, version, last_edit, copyright) = ('RJE_CHLAMYDIA', '0.0', 'August 2011', '2011')
    description = 'Bespoke Chlamydia analysis pipeline'
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
### SECTION II: Chlamydia Class                                                                                         #
#########################################################################################################################
class Chlamydia(rje_obj.RJE_Object):     
    '''
    Chlamydia Class. Author: Rich Edwards (2011).

    Str:str
    
    Bool:boolean

    Int:integer
    - WinSize=X   : Window size for NT composition plots [5000]
    - WinStep=X   : Window step size for NT composition plots [1000]

    Num:float
    
    List:list

    Dict:dictionary
    - SeqList = {type:rje_seqlist.SeqList object}

    Obj:RJE_Objects
    - DB = Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = []
        self.boollist = []
        self.intlist = ['WinSize','WinStep']
        self.numlist = []
        self.listlist = []
        self.dictlist = ['SeqList']
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Basefile':'NC_010287','Mode':'cont'})
        self.setBool({})
        self.setInt({'WinSize':5000,'WinStep':1000})
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
                ### Class Options ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['Mode'])  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'int',['WinSize','WinStep'])  # No need for arg if arg = att.lower()
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Mode') == 'cont': self.contamination()
            else:
                self.setup()
                self.genomeScan()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            seqcmd = self.cmd_list + ['autoload=T','seqmode=file','seqindex=T']
            dfile = '%s.data.tdt' % self.basefile()
            ### ~ [2] Load Sequence Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['SeqList']['full'] = rje_seqlist.SeqList(self.log,seqcmd + ['seqmode=list'])
            self.debug(self.dict['SeqList']['full'].seqNum())
            if self.dict['SeqList']['full'].seqNum(): return
            self.dict['SeqList']['full'] = rje_seqlist.SeqList(self.log,seqcmd + ['seqin=%s.full.fas' % (self.basefile()),'seqmode=list'])
            for stype in ['CDS','gene','prot']:
                seq = self.dict['SeqList'][stype] = rje_seqlist.SeqList(self.log,seqcmd + ['seqin=%s.%s.fas' % (self.basefile(),stype)])
                seq.dict['SeqDict'] = {}
                for s in seq.list['Seq']:
                    (name,sequence) = seq.getSeq(s)
                    seq.dict['SeqDict'][string.split(string.split(name)[0],'_')[-1]] = s
            ### ~ [3] Database Compilation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.exists(dfile) and not self.getBool('Force'): db.addTable(dfile,name='data',mainkeys=['tag'])
            else:
                ## ~ [3a] ~ Load part tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fdb = db.addTable('%s.function.tdt' % self.basefile(),name='function',mainkeys=['tag'])
                fdb.dropField('description')
                edb = db.addTable('%s.expression.tdt' % self.basefile(),name='expression',mainkeys=['key'])
                nx = 0
                edb.fillBlanks(blank='0',fillempty=True)
                for ekey in rje.sortKeys(edb.data()):
                    entry = edb.data(ekey)
                    for field in edb.fields():
                        if entry[field] == 'na': entry[field] = '0.0'; nx += 1
                self.printLog('#TDT','Updated %s entries for expression table' % rje.iStr(nx))
                kdb = db.addTable('%s.proteinkey.tdt' % self.basefile(),name='proteinkey',mainkeys=['key'])
                xdb = db.addTable('%s.dbxref.tdt' % self.basefile(),name='dbxref',mainkeys=['tag'])
                xdb.dropField('gene')   # Pull from genbank instead
                #pdb = db.addTable('%s.cysweight.tdt' % self.basefile(),name='cysweight',mainkeys=['AccNum'])
                pdb = db.addTable('%s.protein.tdt' % self.basefile(),name='prodigis',mainkeys=['AccNum'])
                pdb.addField('NRPep5','NRPep',0); pdb.addField('NRPep7','NRPep5',0)
                for x in range(5,51):
                    xfield = '%d' % x
                    if xfield not in pdb.fields(): continue
                    for entry in pdb.entries():
                        entry['NRPep5'] += int(entry[xfield])
                        if x >= 7: entry['NRPep7'] += int(entry[xfield])
                for field in pdb.fields()[0:]:
                    if field not in ['AccNum','File','ProtMWt','PepCount','LenExp','Len3','Len5','Len7Exp','Len37','NRPep','NRPep5','NRPep7','Cys0']: pdb.dropField(field)
                #pdb.renameField('AccNum','uniprot')
                #pdb.newKey(['uniprot'])
                pdb.renameField('AccNum','tag')
                pdb.newKey(['tag'])
                mdb = db.addTable('%s.PNASmaintable.tdt' % self.basefile(),name='main',mainkeys=['tag'])
                tdb = db.addTable('%s.tmhmm.tdt' % self.basefile(),name='TMHMM',mainkeys=['acc_num'])
                ## ~ [3b] ~ Load and process features table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gdb = db.addTable('%s.Feature.tdt' % self.basefile(),name='feature',mainkeys=['locus','feature','position'])
                gdb.dropEntriesDirect('feature',['CDS'],inverse=True)
                gdb.list['Fields'] += ['tag','start','end','gene','product']
                for entry in gdb.entries():
                    pos = rje.matchExp('(\d+)\.\.(\d+)',entry['position'])
                    if entry['position'][:4] == 'comp': entry['start'] = pos[1]; entry['end'] = pos[0]
                    else: entry['start'] = pos[0]; entry['end'] = pos[1]
                    try: entry['tag'] = rje.matchExp('locus_tag="(\S+)"',entry['details'])[0]
                    except: entry['tag'] = '-'
                    try: entry['gene'] = rje.matchExp('gene="(\S+)"',entry['details'])[0]
                    except: entry['gene'] = ''
                    try: entry['product'] = string.split(string.split(entry['details'],'/product="')[1],'"')[0]
                    except: entry['product'] = ''
                gdb.dropEntriesDirect('tag',['-'])
                gdb.newKey(['tag'])
                for field in ['locus','feature','position','details']: gdb.dropField(field)
                ## ~ [3c] ~ Codon Bias Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cfile = '%s.CDS.Bias.tdt' % self.basefile()
                if not rje.exists(cfile) or self.getBool('Force'):
                    rje_codons.Codons(self.log,self.cmd_list+['seqin=%s.CDS.fas' % self.basefile(),'backups=F']).run()
                bdb = db.addTable(cfile,name='Bias',mainkeys=['Seq'])
                bdb.renameField('Len','AALen')
                ndb = db.addTable('%s.CDS.NT.tdt' % self.basefile(),name='NT',mainkeys=['Seq'])
                ndb.renameField('Len','NTLen')
                for field in ndb.fields():
                    if field != string.replace(field,'U','T'): ndb.renameField(field,string.replace(field,'U','T'))
                ## ~ [3d] ~ Join tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                temp = db.joinTables(name='temp',join=[(edb,'key'),(kdb,'key')],newkey=['key'],cleanup=True,keeptable=True)
                #pfields = pdb.fields()[0:]
                #pfields.remove('uniprot')
                #temp2 = db.joinTables(name='temp2',join=[(xdb,'uniprot'),(pdb,'uniprot',pfields)],newkey=['tag'],cleanup=True,keeptable=True)
                #data = db.joinTables(name='data',join=[(temp2,'tag'),(fdb,'tag'),(gdb,'tag'),(bdb,'Seq'),(ndb,'Seq'),(temp,'tag'),(mdb,'tag')],newkey=['tag'],cleanup=True,keeptable=True)
                data = db.joinTables(name='data',join=[(pdb,'tag'),(xdb,'tag'),(fdb,'tag'),(tdb,'acc_num'),(gdb,'tag'),(bdb,'Seq'),(ndb,'Seq'),(temp,'tag'),(mdb,'tag')],newkey=['tag'],cleanup=True,keeptable=True)
                data.dropField('Seq')
                ## ~ [3e] ~ Fill out data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                data.fillBlanks(blank='0.0',fields=['eb','rb'],fillempty=True)
                #for entry in data.entries():
                #    if entry['tag'] not in self.dict['SeqList']['CDS'].dict['SeqDict']: entry['function'] = 'Non-CDS'
                data.fillBlanks(blank='Unassigned',fields=['function'],fillempty=True)
                data.fillBlanks()
                data.fillBlanks(blank='no mapping',fields=['description'],fillempty=True)
                data.saveToFile(dfile)
                allfields = data.list['Fields'][0:]
                data.list['Fields'] = ["tag","File","PepCount","LenExp","Len3","Len5","Len7Exp","Len37","NRPep",'NRPep5','NRPep7',"Cys0",
                                       "pi","mass","function","new_function","tm","start","end","AALen","Bias",
                                       "WtBias","AbsBias",'NTLen','C','A','G','T','C|3','A|3','G|3','T|3',
                                       'eb_1.1','eb_1.2','eb_2.1','eb_2.2','rb_1.1','rb_1.2','rb_2.1','rb_2.2','eb','rb']
                data.saveToFile('%s.cutdata.tdt' % self.basefile())
                data.list['Fields'] = allfields
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def genomeScan(self):   ### Produces tables of nucleotide composition for R graphical display
        '''Produces tables of nucleotide composition for R graphical display.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','#~~~~~~~~~ Genome Scan ~~~~~~~~~~#')
            cfile = '%s.gene_pos.tdt' % self.basefile()
            gfile = '%s.genome_nt.tdt' % self.basefile()
            if os.path.exists(gfile) and not self.getBool('Force'): return
            GENSCAN = open(gfile,'w')
            GENSCAN.write(string.join(['Pos','G','C','A','T\n'],'\t'))
            genseq = self.dict['SeqList']['full'].nextSeq()[1].upper()
            ### ~ [2] Calculate sliding windows of nt composition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pos = 0; step = self.getInt('WinStep'); win = self.getInt('WinSize')
            while pos < len(genseq):
                winseq = genseq[pos-(win/2):][:win]
                if len(winseq) < win:
                    if (pos-(win/2)) < 0: winseq = genseq[pos-(win/2):] + winseq
                    else: winseq = winseq + genseq[:win-len(winseq)]
                GENSCAN.write('%s' % pos)
                for n in 'GCAT': GENSCAN.write('\t%.3f' % (float(string.count(winseq,n))/len(winseq)))
                GENSCAN.write('\n')
                pos += step
            GENSCAN.close()
            self.printLog('#SCAN','Genomescan output to %s' % gfile)
            ### ~ [3] Calculate protein expression density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.db('data'): return
            self.db().info['Basefile'] = self.basefile()
            cutdb = self.db().copyTable('data','gene_pos')
            cutdb.list['Fields'] = ['tag','start','end','NTLen','C','A','G','T','C|3','A|3','G|3','T|3','eb_1.1','eb_1.2','eb_2.1','eb_2.2','rb_1.1','rb_1.2','rb_2.1','rb_2.2','eb','rb']
            cutdb.saveToFile(cfile)
            pdb = self.db().addEmptyTable('protein_density',['Pos','eb_p','eb_n','rb_p','rb_n'],['Pos'])
            cutdb.dataFormat({'start':'int','end':'int','eb':'num','rb':'num'})
            pos = 0; win = 10000
            while pos < len(genseq):
                pos += win
                pentry = {'Pos':pos,'eb_p':0.0,'eb_n':0.0,'rb_p':0.0,'rb_n':0.0}
                for entry in cutdb.entries():
                    if min(entry['start'],entry['end']) <= pos and max(entry['start'],entry['end']) > (pos - win):
                        strand = {True:'p',False:'n'}[entry['start'] < entry['end']]
                        outside = max(0,(pos-win) - min(entry['start'],entry['end']))
                        outside += max(0,max(entry['start'],entry['end']) - pos)
                        prop = max(entry['start'],entry['end']) - min(entry['start'],entry['end'])
                        prop = float(prop - outside) / prop
                        pentry['eb_%s' % strand] += prop * entry['eb']
                        pentry['rb_%s' % strand] += prop * entry['rb']
                pdb.addEntry(pentry)
            pdb.saveToFile()
        except: self.errorLog('%s.genomeScan error' % self)
#########################################################################################################################
    def contamination(self):    ### Compares peptides from Chlamydia and human and outputs summaries
        '''Compares peptides from Chlamydia and human and outputs summaries.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            mods = ['none']
            
            ### >>>> Shortcut reanalysis without modifications >>>> ###
            pepfile = '%s.chlam_peptides.tdt' % self.basefile()
            if not self.force() and os.path.exists(pepfile):
                pepdb = db.addTable(pepfile,mainkeys=['key','seqmod'],name='chlam_nomod')
                pepdb.dropFields(['pass','modification'])
                pepdb.compress(['key','seq'],default='max')
                pepdb.dropFields(['seqmod'])
                for entry in pepdb.entries():
                    for field in pepdb.fields():
                        if 'len' not in field: continue
                        try:
                            if entry[field] and int(entry[field]): entry[field] = len(entry['seq'])
                            else: entry[field] = ''
                        except: self.errorLog('%s >> %s' % (entry,field),quitchoice=True)
                tdb = pepdb
                comprules = {'key':'str','pi':'str','mass':'str'}
                shapefields = []
                for field in pepdb.fields():
                    if 'len' in field: comprules[field] = 'mean'
                    if len(string.split(field,'|')) > 1 and string.split(field,'|')[0] not in shapefields: shapefields.append(string.split(field,'|')[0])
                print shapefields
                tdb.compress(['protein'],rules=comprules,default='sum')
                tdb.dropFields(['seq'])
                tdb.saveToFile()

                tdb.info['Name'] = 'chlam_nomod_summary'
                tdb.addField('temp',evalue=1)
                tdb.compress(['temp'],rules=comprules,default='sum')
                tdb.reshapeLong('exp',shapefields)
                tdb.newKey(['exp'])
                tdb.dropFields(['exp']+shapefields,inverse=True)
                tdb.saveToFile()

                return
            ### <<<< End Shortcut reanalysis without modifications <<<< ###

            ## ~ [0a] ~ Load EB and RB human peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','## ~ [0a] ~ Load EB and RB human peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##',log=False)
            #protein.key	protein.Entry	protein.Accession	protein.Description	protein.dataBaseType	protein.score	protein.falsePositiveRate	protein.avgMass	protein.MatchedProducts	protein.matchedPeptides	protein.digestPeps	protein.seqCover(%)	protein.MatchedPeptideIntenSum	protein.top3MatchedPeptideIntenSum	protein.MatchedProductIntenSum	protein.fmolOnColumn	protein.ngramOnColumn	protein.AutoCurate	protein.Key_ForHomologs	protein.SumForTotalProteins	peptide.Rank	peptide.Pass	peptide.matchType	peptide.modification	peptide.mhp	peptide.seq	peptide.OriginatingSeq	peptide.seqStart	peptide.seqLength	peptide.pI	peptide.componentID	peptide.MatchedProducts	peptide.UniqueProducts	peptide.ConsectiveMatchedProducts	peptide.ComplementaryMatchedProducts	peptide.rawScore	peptide.score	peptide.(X)-P Bond	peptide.MatchedProductsSumInten	peptide.MatchedProductsTheoretical	peptide.MatchedProductsString	peptide.ModelRT	peptide.Volume	peptide.CSA	peptide.ModelDrift	peptide.RelIntensity	peptide.AutoCurate	precursor.leID	precursor.mhp	precursor.mhpCal	precursor.retT	precursor.inten	precursor.calcInten	precursor.charge	precursor.z	precursor.mz	precursor.fraction	precursor.numFrac	precursor.fwhm	precursor.liftOffRT	precursor.infUpRT	precursor.infDownRT	precursor.touchDownRT	prec.rmsFWHMDelta	peptidePrecursor.deltaMhpPPM
            humedb = db.addTable('EB_IA_final_peptide.csv',mainkeys=['protein.Accession','peptide.Rank','peptide.seq','peptide.modification'],datakeys=['protein.key','protein.Accession','protein.Entry','peptide.Rank','peptide.Pass','peptide.seq','peptide.modification','peptide.OriginatingSeq'],name='humaneb')
            humrdb = db.addTable('RB_IA_final_peptide.csv',mainkeys=['protein.Accession','peptide.Rank','peptide.seq','peptide.modification'],datakeys=['protein.key','protein.Accession','protein.Entry','peptide.Rank','peptide.Pass','peptide.seq','peptide.modification','peptide.OriginatingSeq'],name='humanrb')
            for humdb in [humedb,humrdb]:
                humdb.info['Delimit'] = '\t'
                humdb.addField('exp',evalue=humdb.info['Name'][-2:])
                humdb.renameField('protein.Accession','Protein')
                humdb.renameField('protein.Entry','Species')
                for entry in humdb.entries(): entry['Species'] = string.split(entry['Species'],'_')[-1]
                humdb.dropEntriesDirect('Species',['HUMAN'],inverse=True)
                for field in ['Rank','Pass','seq','OriginatingSeq','modification']: humdb.renameField('peptide.%s' % field,field)
                humdb.dataFormat({'Rank':'int'})
                for mod in humdb.index('modification'):
                    if mod.lower() and mod.lower() not in mods: mods.append(mod.lower())
                humdb.addField('seqmod')
                for entry in humdb.entries():
                    if entry['modification'] and mods.index(entry['modification'].lower()): entry['seqmod'] = '%s-%d' % (entry['seq'],mods.index(entry['modification'].lower()))
                    else: entry['seqmod'] = entry['seq']
            humtdb = db.copyTable(humedb,'humantot')
            humtdb.newKey(['Protein','Rank','seq','modification','exp'])
            db.mergeTables(humtdb,db.copyTable(humrdb,'temp',add=False))
            humtdb.compress(['Protein','seq','Pass'],rules={'Rank':'max'})
            ## ~ [0b] ~ Load Proteomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','## ~ [0b] ~ Load Proteomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##',log=False)
            # Load human proteome
            hseqfile = '/home/re1u06/researchfiles/SBSBINF/Databases/DBase_120225/EnsEMBL/ens_HUMAN.loci.fas'
            hseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % hseqfile])
            # Load Chlamydia proteome
            cseqfile = '../2011-07-18-Genome/NC_010287.proteome.fas'
            cseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % cseqfile])
            # Load matched protein list
            rbpep = rje.listFromCommand('../2011-05-ProDigIS/soton_rb_peptides.txt')
            ebpep = rje.listFromCommand('../2011-05-ProDigIS/soton_rb_peptides.txt')
            ## ~ [0c] ~ Load EB and RB Chlamydia peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','## ~ [0c] ~ Load EB and RB Chlamydia peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##',log=False)
            chlamdb = {'EB':[],'RB':[]}
            for pfile in glob.glob('./Soton*/*peptide.csv'):
                (er,uniq) = rje.matchExp('\./Soton(\S\S)\S+_(\d+)/',pfile)
                chlamdb[er].append(db.addTable(pfile,mainkeys=['protein.key','protein.name','peptide.Rank','peptide.seq','peptide.modification'],datakeys=['protein.name','peptide.Rank','peptide.Pass','peptide.seq','peptide.modification','peptide.OriginatingSeq'],name=uniq))
            edb = chlamdb['EB'].pop(0); edb.info['Name'] = 'chlam_eb'
            while chlamdb['EB']: db.mergeTables(edb,chlamdb['EB'].pop(0))
            rdb = chlamdb['RB'].pop(0); rdb.info['Name'] = 'chlam_rb'
            while chlamdb['RB']: db.mergeTables(rdb,chlamdb['RB'].pop(0))
            # Load EB and RB matching peptide file
            #edb = db.addTable('../2011-05-ProDigIS/SotonEB_peptide_pjss.csv',mainkeys=['protein.name','peptide.Rank','peptide.seq','peptide.modification'],datakeys=['protein.name','peptide.Rank','peptide.Pass','peptide.seq','peptide.OriginatingSeq'],name='chlam_eb')
            #rdb = db.addTable('../2011-05-ProDigIS/SotonRB_peptide_pjss.csv',mainkeys=['protein.name','peptide.Rank','peptide.seq','peptide.modification'],datakeys=['protein.name','peptide.Rank','peptide.Pass','peptide.seq','peptide.OriginatingSeq'],name='chlam_rb')
            for chlamdb in [edb,rdb]:
                chlamdb.info['Delimit'] = '\t'
                chlamdb.addField('exp',evalue=chlamdb.info['Name'][-2:])
                chlamdb.renameField('protein.name','Protein'); chlamdb.renameField('protein.key','key')
                for field in ['Rank','Pass','seq','OriginatingSeq','modification']: chlamdb.renameField('peptide.%s' % field,field)
                chlamdb.dataFormat({'Rank':'int'})
                for mod in chlamdb.index('modification'):
                    if mod.lower() and mod.lower() not in mods: mods.append(mod.lower())
                chlamdb.addField('seqmod')
                chlamdb.addField('Species',evalue='UNKNOWN')
                for entry in chlamdb.entries():
                    if 'Chlamydia trachomatis' in entry['Protein'] or '_CHLT2' in entry['Protein']: entry['Species'] = 'CHLT2'
                    if entry['modification'] and mods.index(entry['modification'].lower()): entry['seqmod'] = '%s-%d' % (entry['seq'],mods.index(entry['modification'].lower()))
                    else: entry['seqmod'] = entry['seq']
                    if not entry['OriginatingSeq']: entry['OriginatingSeq'] = entry['seq']
                chlamdb.dropEntriesDirect('Species',['CHLT2'],inverse=True)
                chlamdb.remakeKeys()
            ## ~ Load Protein Key Mapping ~ ##
            kdb = db.addTable('NC_010287.proteinkey.tdt',mainkeys=['key'],name='keys')
            xdb = db.addTable('NC_010287.dbxref.tdt',mainkeys=['tag'],name='xref')
            tdb = db.copyTable(edb,'chlam_temp')
            self.deBug(tdb.entries()[0])
            tdb.newKey(['Protein','Rank','Pass','seq','modification','exp'])
            db.mergeTables(tdb,db.copyTable(rdb,'temp',add=False))
            kdb = db.joinTables(name='full_xref',join=[(kdb,'tag'),(xdb,'tag')],newkey=kdb.keys(),keeptable=True)
            tdb = db.joinTables(name='chlam_tot',join=[(tdb,'key'),(kdb,'key')],newkey=tdb.keys(),keeptable=True)
            self.deBug(tdb.keys())
            self.deBug(tdb.entries()[0])

            ### ~ [1] ~ Add Human Data to combined Chlamydia Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','### ~ [1] ~ Add Human Data to combined Chlamydia Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###',log=False)
            tdb.renameField('Pass','pass'); tdb.renameField('Protein','protein');
            for entry in tdb.entries():
                entry['pass'] = string.atoi(entry['pass'][-1])
            keep = ['key'] + xdb.fields() + ['description','protein','exp','pass','seq','seqmod','modification']
            tdb.newKey(['tag','exp','pass','seqmod','Rank'])
            tdb.compress(['tag','exp','pass','seqmod'],rules={'pass':'min'})
            tdb.dropFields(keep,inverse=True)
            self.deBug(tdb.keys())
            self.deBug(tdb.entries()[0])

            ## ~ [1a] ~ Map ID'd peptides onto Chlamydia, human and EB/RB hits ~~~~~~~~~~~~~~~~~~~~ ##
            bothpep = {'eb':[],'rb':[]}     # Peptides found in both species
            chlampep = {'eb':[],'rb':[]}    # Peptides found only in Chlamydia
            uniqpep = {'eb':[],'rb':[]}     # Peptides found only in a single protein in Chlamydia
            fx = len(tdb.fields())
            tdb.addField('pep',evalue=1)
            for field in ['pass1','pass2','hsap1','hsap2','uniq1','uniq2']: tdb.addField(field,evalue=0)
            comprules = {'pass':'min','key':'min'}
            for field in ['pep','pass1','pass2','hsap','uniq']:
                tdb.addField('%s_len' % field)
                comprules[tdb.fields()[-1]] = 'mean'
            shapefields = tdb.fields()[fx:]

            for entry in tdb.entries():
                epass = 'Pass%d' % entry['pass']
                entry[epass.lower()] = 1
                plen = entry['pep_len'] = len(entry['seq'])
                plen = entry['pass%d_len' % entry['pass']] = len(entry['seq'])
                hsap = False
                if entry['exp'] == 'eb':
                    if 'Pass1' in humedb.indexDataList('seqmod',entry['seqmod'],'Pass'): entry['hsap1'] += 1; hsap = True
                    if 'Pass2' in humedb.indexDataList('seqmod',entry['seqmod'],'Pass'): entry['hsap2'] += 1; hsap = True
                    if entry['seqmod'] not in humedb.index('seqmod'):
                        if entry['seq'] in humedb.index('seq'): self.errorLog('EB mod peptide %s not found in Human EB but unmod *is* found in Human EB!' % entry['seqmod'],printerror=False)
                        if entry['seqmod'] in humrdb.index('seqmod'): self.errorLog('EB peptide %s not found in Human EB but found in Human RB!' % entry['seqmod'],printerror=False)
                else:
                    if 'Pass1' in humrdb.indexDataList('seqmod',entry['seqmod'],'Pass'): entry['hsap1'] += 1; hsap = True
                    if 'Pass2' in humrdb.indexDataList('seqmod',entry['seqmod'],'Pass'): entry['hsap2'] += 1; hsap = True
                    if entry['seqmod'] not in humrdb.index('seqmod'):
                        if entry['seq'] in humrdb.index('seq'): self.errorLog('RB mod peptide %s not found in Human RB but unmod *is* found in Human RB!' % entry['seqmod'],printerror=False)
                        if entry['seqmod'] in humedb.index('seqmod'): self.errorLog('RB peptide %s not found in Human RB but found in Human EB!' % entry['seqmod'],printerror=False)
                if hsap: entry['hsap_len'] = plen; bothpep[entry['exp']].append(entry['seq']); continue
                chlampep[entry['exp']].append(entry['seq'])
                entry['uniq1'] = entry['pass1']
                entry['uniq2'] = entry['pass2']
                entry['uniq_len'] = plen
                for altentry in tdb.indexEntries('seqmod',entry['seqmod']):
                    if altentry['tag'] == entry['tag']: continue
                    entry['uniq1'] = entry['uniq2'] = entry['uniq_len'] = 0
                if entry['uniq1'] or entry['uniq2']: uniqpep[entry['exp']].append(entry['seq'])

            tdb.reshapeWide('exp',shapefields)
            fillfields = tdb.fields()[13:]
            for field in fillfields[0:]:
                if 'len' in field: fillfields.remove(field)
            tdb.fillBlanks(0,fillfields,fillempty=True)
            for entry in tdb.entries():
                if entry['modification'] == 0: entry['modification'] = ''

            for field in shapefields:
                tdb.addField('%s|tot' % field)
                if field[-3:] == 'len':
                    comprules['%s|eb' % field] = 'mean'
                    comprules['%s|rb' % field] = 'mean'
                    comprules['%s|tot' % field] = 'mean'
            for entry in tdb.entries():
                for field in shapefields:
                    if entry['%s|eb' % field] and not entry['%s|rb' % field]: entry['%s|tot' % field] = entry['%s|eb' % field]
                    elif entry['%s|rb' % field] and not entry['%s|eb' % field]: entry['%s|tot' % field] = entry['%s|rb' % field]
                    else: entry['%s|tot' % field] = max(entry['%s|eb' % field],entry['%s|rb' % field])

            tdb.info['Name'] = 'chlam_peptides'
            tdb.saveToFile()

            tdb.info['Name'] = 'chlam_proteins'
            tdb.compress(['protein'],rules=comprules,default='sum')
            tdb.dropFields(['pass','seq','modification','seqmod'])
            tdb.saveToFile()

            tdb.info['Name'] = 'chlam_summary'
            tdb.addField('temp',evalue=1)
            tdb.compress(['temp'],rules=comprules,default='sum')
            tdb.reshapeLong('exp',shapefields)
            tdb.newKey(['exp'])
            tdb.dropFields(['exp']+shapefields,inverse=True)
            tdb.saveToFile()

            bothpep['tot'] = bothpep['eb'] + bothpep['rb']     # Peptides found in both species
            chlampep['tot'] = chlampep['eb'] + chlampep['rb']# Peptides found only in Chlamydia
            uniqpep['tot'] = uniqpep['eb'] + uniqpep['rb']
            for er in bothpep:
                open('%s.%s.bothpep.txt' % (self.basefile(),er),'w').write(string.join(rje.sortUnique(bothpep[er]),'\n'))
                open('%s.%s.chlampep.txt' % (self.basefile(),er),'w').write(string.join(rje.sortUnique(chlampep[er]),'\n'))
                open('%s.%s.uniqpep.txt' % (self.basefile(),er),'w').write(string.join(rje.sortUnique(uniqpep[er]),'\n'))
            



            return



            #Peptide numbers for C. trachomatis/human
            #1.	Number of  chlamydial peptides assigned for each protein from RBs
            #2.	Number of  chlamydial peptides assigned for each protein from EBs
            #3.	Number of  chlamydial peptides assigned from both EB and RB combined, with redundancy removed
            #4.	Number of  unique chlamydial peptides assigned for each protein from RBs
            #5.	Number of unique chlamydial peptides assigned for each protein from EBs
            #6.	Number of unique chlamydial peptides assigned for EBs and RBs combined with redundancy removed
            #7.	Total number of human peptides identified in EB (Length would be useful)
            #8.	Total number of human peptides identified in RB (Length would be useful)
            #9.	Total number of human peptides identified in EB and RB
            #10.	Human peptides matching pass 1 chlamydia peptides for RB (sequence would be useful)
            #11.	Human peptides matching pass 2 chlamydia peptides for EB (sequence would be useful)

            #An accession number and protein description would be useful where possible, i.e.,  the number of chlamydial peptides for each protein.





            tdb.compress(['Protein','seq'],rules={'Rank':'max','Pass':'list'})
            for entry in tdb.entries(): entry['exp'] = 'tot'

            
            ### ~ [1] ~ Map ID'd peptides onto Chlamydia, human and EB/RB hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapkey = 'seqmod'
            self.deBug(rje.sortKeys(humdb.index(mapkey)))
            self.printLog('#~~#','## ~ [1] ~ Map ID\'d peptides onto Chlamydia, human and EB/RB hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###',log=False)
            for chlamdb in [edb,rdb,tdb]:
                bothpep = []
                chlampep = []
                uniqpep = []
                if chlamdb == edb: humdb = humedb
                elif chlamdb == rdb: humdb = humrdb
                else: humdb = humtdb
                chlamdb.addField('Pep',evalue=1)
                chlamdb.addField('Pass1',evalue=0)
                chlamdb.addField('Pass2',evalue=0)
                chlamdb.addField('Hsap1',evalue=0)
                chlamdb.addField('Hsap2',evalue=0)
                chlamdb.addField('Uniq1',evalue=0)
                chlamdb.addField('Uniq2',evalue=0)
                comprules = {'Rank':'max'}
                for field in ['Pep','Pass1','Pass2','Hsap','Uniq']:
                    chlamdb.addField('%s_len' % field)
                    comprules[chlamdb.fields()[-1]] = 'mean'
                for entry in chlamdb.entries():
                    if 'Pass1' in entry['Pass']: entry['Pass'] = 'Pass1'
                    else: entry['Pass'] = 'Pass2'
                    entry[entry['Pass']] += 1
                    entry['Pep_len'] = plen = len(entry['seq'])
                    entry['%s_len' % entry['Pass']] = plen
                    hsap = False
                    self.deBug(entry[mapkey])
                    self.deBug(entry[mapkey] in humdb.index(mapkey))
                    if 'Pass1' in humdb.indexDataList(mapkey,entry[mapkey],'Pass'): entry['Hsap1'] += 1; bothpep.append(entry[mapkey]); hsap = True
                    if 'Pass2' in humdb.indexDataList(mapkey,entry[mapkey],'Pass'): entry['Hsap2'] += 1; bothpep.append(entry[mapkey]); hsap = True
                    if hsap: entry['Hsap_len'] = plen; continue
                    chlampep.append(entry[mapkey])
                    entry['Uniq1'] = entry['Pass1']
                    entry['Uniq2'] = entry['Pass2']
                    entry['Uniq_len'] = plen
                    for altentry in chlamdb.indexEntries(mapkey,entry[mapkey]):
                        if altentry['Protein'] == entry['Protein']: continue
                        entry['Uniq1'] = entry['Uniq2'] = 0
                    if entry['Uniq1'] or entry['Uniq2']: uniqpep.append(entry[mapkey])
                chlamdb.dropFields(['Pass','Rank','seq','OriginatingSeq','modification'])
                chlamdb.compress(['Protein'],rules=comprules,default='sum')
                #chlamdb.dropField('Rank')
                chlamdb.saveToFile()
                open('%s.%s.bothpep.txt' % (self.basefile(),chlamdb.info['Name']),'w').write(string.join(rje.sortUnique(bothpep),'\n'))
                open('%s.%s.chlampep.txt' % (self.basefile(),chlamdb.info['Name']),'w').write(string.join(rje.sortUnique(chlampep),'\n'))
                open('%s.%s.uniqpep.txt' % (self.basefile(),chlamdb.info['Name']),'w').write(string.join(rje.sortUnique(uniqpep),'\n'))
                chlamdb.newKey(['Protein','exp'])
            db.mergeTables(edb,rdb)
            db.mergeTables(edb,tdb)
            cdb = db.copyTable(edb,'chlam_summary')
            edb.info['Name'] = 'chlam_pep'
            edb.reshapeWide('exp',edb.fields()[-7:])
            edb.saveToFile()
            cdb.compress(['exp'],rules=comprules,default='sum')
            cdb.dropField('Protein')
            cdb.saveToFile()
            


            # - twice maybe, once using EnsEMBL sequences directly, once using EB/RB search
            # - Numbers of unique Pass1/2 human peptides, and numbers matching Chlam
            # - Numbers of matched peptides per Chlam gene: total, eb, rb, human (e/r), unique (e/r), ens (e/r)
            ## Do complete digest of Chlam and search against Human
        except: self.errorLog('%s.contamination error' % self)
#########################################################################################################################
### End of SECTION II: Chlamydia Class                                                                                  #
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
    try: Chlamydia(mainlog,cmd_list).run()

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
