#!/usr/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_yeast
Description:  Yeast PPI & Sequence Module
Version:      0.0
Last Edit:    22/12/10
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This Yeast SLiMFinder Analysis Module is designed to read, assess, store and combine the different yeast information
    for SLiMFinder analysis, using the YGOB for generating alignments for conservation masking etc.

    The main (default) inputs are:
    - Y2H_union.txt = High quality binary interactions from Vidal lab
    - Pillars.tab = YGOB orthology pillars for yeast species
    - Proteins.fas = Fasta file of protein sequences to match pillars.tab
    
Commandline:
    ### ~ INPUT DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence file containing yeast sequences [Proteins.fas]
    pillars=FILE    : YGOB pillars file [Pillars.tab]
    ppifile=FILE    : Input PPI data (two columns of binary interactors) [Y2H_union.txt]
    xref=FILE       : Yeast identifier XRef file (e.g. BioMart download) [yeast_xref.20101222.tdt]
    sgd2sp=T/F      : Convert SGD identifiers into SwissProt identifiers [False]
    gopher=T/F      : Whether to feed Pillars into GOPHER Orthology (at BLAST ID stage) [False]
    
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_sequence, rje_uniprot, rje_zen
import rje_blast_V1 as rje_blast
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
    (program, version, last_edit, copyright) = ('YEAST', '0.0', 'December 2010', '2010')
    description = 'Yeast PPI & Sequence Module'
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
### SECTION II: Yeast Class                                                                                             #
#########################################################################################################################
class Yeast(rje.RJE_Object):     
    '''
    Yeast Class. Author: Rich Edwards (2008).

    Info:str
    - Pillars = YGOB pillars file [Pillars.tab]
    - PPIFile = Input PPI data (two columns of binary interactors) [Y2H_union.txt]
    - XRef = Yeast identifier XRef file (e.g. BioMart download) [yeast_xref.20101222.tdt]
    
    Opt:boolean
    - Gopher = Whether to feed Pillars into GOPHER Orthology (at BLAST ID stage) [False]
    - SGD2SP = Convert SGD identifiers into SwissProt identifiers [False]

    Stat:numeric

    List:list
    - Pillars = List of lists of YGOB gene pillars
    - YeastSeq = AccNum of Scer sequences

    Dict:dictionary
    - PPI = Dictionary of {hub:[spokes]}
    - Rename = Dictionary for renaming to UniProt IDs

    Obj:RJE_Objects
    - DB = Database object
    - GOPHER = GOPHER Object
    - SeqList = Sequence list storing yeast sequences [Proteins.fsa]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Pillars','PPIFile','XRef']
        self.optlist = ['SGD2SP','Gopher']
        self.statlist = []
        self.listlist = ['Pillars','YeastSeq']
        self.dictlist = ['PPI','Rename']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Pillars':'Pillars.tab','PPIFile':'Y2H_union.txt','XRef':'yeast_xref.20101222.tdt'})
        ### Other Attributes ###
        self.obj['SeqList'] = rje_seq.SeqList(self.log,['accnr=F','seqnr=F','autoload=T','seqin=Proteins.fsa']+self.cmd_list)
        self.dict['SeqDict'] = self.obj['SeqList'].seqNameDic(proglog=self.stat['Verbose']>0)
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['Pillars','PPIFile','XRef'])
                self._cmdReadList(cmd,'opt',['SGD2SP','Gopher'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.reformatSeq()
            if self.opt['Gopher']: self.gopher()
            self.makePPIDatasets()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.loadPPI()          # Load pairwise interaction data
            self.loadPillars()      # Load YGOB Pillar data
            self.loadXRef()
            #x#self.checkSeqList()     # Check sequence integrity
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Data Loading Class Methods                                                                                #
#########################################################################################################################
    def loadPPI(self):  ### Load pairwise interaction data
        '''Load pairwise interaction data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(self.info['PPIFile']): return False
            ### ~ [2] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in open(self.info['PPIFile'],'r').readlines():
                try: [pa,pb] = string.split(rje.chomp(line))[:2]
                except: continue
                for ppi in [(pa,pb),(pb,pa)]:
                    if ppi[0] not in self.dict['PPI']: self.dict['PPI'][ppi[0]] = []
                    if ppi[1] not in self.dict['PPI'][ppi[0]]: self.dict['PPI'][ppi[0]].append(ppi[1])
                self.progLog('\r#PPI','Loading PPI data: %s proteins' % rje.integerString(len(self.dict['PPI'])))
            self.printLog('\r#PPI','Loaded PPI data for %s proteins' % rje.integerString(len(self.dict['PPI'])))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def loadPillars(self):  ### Load YGOB Pillar data
        '''Load YGOB Pillar data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(self.info['Pillars']): return False
            ### ~ [2] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in self.loadFromFile(filename=self.info['Pillars'],chomplines=True):
                pillars = string.split(line)
                #self.deBug('%s = %d' % (pillars,len(pillars)))
                if len(pillars) < 17: continue
                pillars = pillars[:5] + pillars[6:]     # Remove ancestral gene
                while '---' in pillars: pillars.remove('---')
                #self.deBug('%s = %d' % (pillars,len(pillars)))
                if pillars: self.list['Pillars'].append(pillars)
                self.progLog('\r#YGOB','Loading Pillar data: %s loci' % rje.integerString(len(self.list['Pillars'])))
            self.printLog('\r#YGOB','Loaded Pillar data for %s loci' % rje.integerString(len(self.list['Pillars'])))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def loadXRef(self):     ### Load Identifier XRef Data
        '''Load Identifier XRef Data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.exists('%s.xref.tdt' % self.info['Basefile']) and not self.opt['Force']: 
                return self.db().addTable('%s.xref.tdt' % self.info['Basefile'],mainkeys=['#'],datakeys='All',name='XRef')
            if not rje.checkForFile(self.info['XRef']): return False
            changehead = {'Ensembl Gene ID':'EnsG','Ensembl Protein ID':'EnsP','Associated Gene Name':'Gene',
                          'Associated Gene DB':'GeneDB','UniProt/SwissProt ID':'UniprotID',
                          'UniProt/SwissProt Accession':'UniProt','SGD Gene':'SGD'}
            ### ~ [2] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xref = self.db().addTable(self.info['XRef'],mainkeys='All',datakeys='All',name='XRef')
            for field in changehead:
                if field in xref.fields(): xref.renameField(field,changehead[field])
            xref.saveToFile('%s.xref.tdt' % self.info['Basefile']); return xref
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def checkSeqList(self):  ### Check sequence integrity
        '''Check sequence integrity.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqdict = self.dict['SeqDict']
            ### ~ [2] Check PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []
            self.progLog('\r#MISS','Checking PPI IDs: %d missing' % len(missing))
            for p in rje.sortKeys(self.dict['PPI']):
                if p not in seqdict: missing.append(p); self.progLog('\r#MISS','Checking PPI IDs: %d missing' % len(missing))
            self.printLog('\r#MISS','Checking PPI IDs complete: %d missing' % len(missing))
            open('yeast.ppi.missing.txt','w').write(string.join(missing,'\n'))
            ### ~ [3] Check Pillar data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = []
            self.progLog('\r#MISS','Checking Pillar IDs: %d missing' % len(missing))
            for pillar in self.list['Pillars']:
                for p in pillar:
                    if p not in seqdict: missing.append(p); self.progLog('\r#MISS','Checking Pillar IDs: %d missing' % len(missing))
            self.printLog('\r#MISS','Checking Pillar IDs complete: %d missing' % len(missing))
            open('yeast.pillar.missing.txt','w').write(string.join(missing,'\n'))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <4> ### PPI Dataset generation methods                                                                          #
#########################################################################################################################
    def reformatSeq(self):  ### Reformats yeast sequence names
        '''Reformats yeast sequence names.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SGD2SP']: return self.sgd2sp()
            if self.obj['SeqList'].info['Name'] != 'Proteins.fsa': return
            self.obj['SeqList'].info['Name'] = 'ygob.fas'
            ### ~ [2] ~ Reformat sequences and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SEQ','Reformatting sequence names...',log=False)
            for seq in self.obj['SeqList'].seq:
                seq.opt['Yeast'] = True
                seq.extractDetails(gnspacc=True)
            self.obj['SeqList'].saveFasta()
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def sgd2sp(self):   ### Reformats yeast sequence names and outputs new data for GOPHER
        '''Reformats yeast sequence names and outputs new data for GOPHER.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            inseq = self.obj['SeqList']
            uni = rje_uniprot.UniProt(self.log,self.cmd_list+['datout=None'])
            xref = self.db('XRef')
            self.dict['Rename'] = {}
            ## ~ [1a] ~ Check or Make UniProt extraction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ufile = '%s.dat' % self.info['Basefile']
            if os.path.exists(ufile) and not self.opt['Force']: uni.readUniProt(ufile,clear=True,cleardata=False)
            else:
                uni.readUniProt(clear=True,acclist=rje.sortKeys(xref.index('UniProt')),cleardata=False)
                uni.saveUniProt(ufile)
            ## ~ [1b] ~ Make dictionary of UniProt sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uniseq = {}
            for entry in uni.entries():
                seq = entry.obj['Sequence']
                uniseq[seq.info['AccNum']] = seq
            self.printLog('\r#USEQ','%s UniProt Sequences extracted (%s Ensembl AccNum)' % (rje.iStr(len(uniseq)), rje.iStr(len(xref.index('UniProt')))))
            ### ~ [2] ~ Reformat sequences and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            yseq = []       # List of YEAST sequence objects
            (sx,stot) = (0.0,inseq.seqNum())
            for seq in inseq.seqs():
                self.progLog('\r#SEQ','Reformatting sequence names: %.2f%%' % (sx/stot)); sx += 100.0
                if seq.info['SpecCode'] != 'YEAST': continue
                yseq.append(seq)
                sgd = seq.info['AccNum']; newname = seq.info['Name']
                try:
                    for x in xref.indexEntries('EnsG',sgd):
                        acc = x['UniProt']
                        if acc: newname = '%s [Gene:%s EnsG:%s SGD:%s AccNum:%s]' % (seq.info['Name'],x['Gene'],x['EnsG'],x['SGD'],acc)
                        else: newname = '%s [Gene:%s EnsG:%s SGD:%s AccNum:-]' % (seq.info['Name'],x['Gene'],x['EnsG'],x['SGD']); continue
                        if acc not in uniseq: self.printLog('\r#UNIERR','Unable to find UniProt sequence %s (%s)' % (acc,sgd)); continue
                        useq = uniseq[acc]
                        if useq.info['Sequence'] != seq.info['Sequence']: self.printLog('\r#SEQERR','%s sequence <> %s sequence' % (sgd,acc)); continue
                        nsplit = string.split(newname)
                        nsplit[0] = '%s__%s' % (x['UniprotID'],acc)
                        newname = string.join(nsplit)
                        self.dict['Rename'][sgd] = acc
                        break
                except: self.errorLog('%s problem' % sgd)
                seq.info['Name'] = newname
                seq.extractDetails(gnspacc=True)
            self.printLog('\r#SEQ','Reformatting sequence names complete.')
            ## ~ [2a] ~ Save renamed sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.exists('%s.ygob.fas' % self.info['Basefile']):
                inseq.saveFasta(seqfile='%s.ygob.fas' % self.info['Basefile'])
            if not rje.exists('%s.yeast.fas' % self.info['Basefile']):
                inseq.saveFasta(seqs=yseq,seqfile='%s.yeast.fas' % self.info['Basefile'])
            self.list['YeastSeq'] = inseq.accList(yseq)
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def gopher(self):  ### Sets up data for GOPHER run
        '''Sets up data for GOPHER run.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'BLAST/')
            rje_blast.BLASTRun(self.log,self.cmd_list).formatDB(fasfile='%s.ygob.fas' % self.info['Basefile'],protein=True,force=False)
            rje_blast.BLASTRun(self.log,self.cmd_list).formatDB(fasfile='%s.yeast.fas' % self.info['Basefile'],protein=True,force=False)
            seqdict = self.obj['SeqList'].seqNameDic('AccNum')
            ymap = self.dict['PillarMap'] = {}
            ### ~ [2] Convert Pillars to BLAST IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (px,ptot) = (0.0,len(self.list['Pillars'])); ox = 0
            for pillar in self.list['Pillars']:
                self.progLog('\r#YGOB','Converting YGOB Pillars for GOPHER: %.2f%%' % (px/ptot)); px += 100
                newpillar = []
                for yid in pillar:
                    seq = rje_sequence.Sequence(self.log,self.cmd_list)
                    seq.opt['Yeast'] = True
                    #self.deBug(yid)
                    seq.info['Name'] = yid
                    seq.extractDetails(gnspacc=True)
                    #self.deBug(seq.info)
                    ygob = seq.info['AccNum']
                    if ygob in self.dict['Rename']: acc = self.dict['Rename'][ygob]
                    else: acc = ygob
                    ymap[yid] = acc
                    if acc not in seqdict: self.printLog('\r#GENE','Non-coding gene %s (%s)? Cannot find in fasta file' % (acc,yid)); continue
                    try:
                        newpillar.append(seqdict[acc].shortName())
                    except:
                        print yid, ygob, acc
                        self.errorLog(rje_zen.Zen().wisdom())
                if not newpillar: continue
                for ygob in pillar:
                    acc = ymap[ygob]
                    if acc not in seqdict: continue
                    if acc in self.list['YeastSeq'] or (not self.list['YeastSeq'] and seqdict[acc].info['SpecCode'] == 'YEAST'):
                        open(rje.makePath('BLAST/%s.blast.id' % acc,wholepath=True),'w').write(string.join(newpillar,'\n'))
                        ox += 1
            self.progLog('\r#YGOB','Converted YGOB Pillars for GOPHER: %s BLAST ID files.' % rje.iStr(ox))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
    def makePPIDatasets(self):  ### Generate PPI datasets from pairwise data
        '''Generate PPI datasets from pairwise data.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'YeastPPI/')
            seqdict = self.dict['SeqDict']
            ### ~ [2] Parse data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot,fx) = (0.0,len(self.dict['PPI']),0)
            for hub in rje.sortKeys(self.dict['PPI']):
                self.progLog('\r#FAS','Generating %s PPI fasta files: %.2f' % (rje.integerString(fx),hx/htot)); hx += 100.0
                if len(self.dict['PPI'][hub]) < 3: continue
                seqs = []
                for spoke in self.dict['PPI'][hub]:
                    if spoke not in seqdict: continue
                    seqs.append(seqdict[spoke])
                if len(seqs) < 3: continue
                self.obj['SeqList'].saveFasta(seqs,rje.makePath('YeastPPI/%s.fas' % hub,wholepath=True),log=False); fx+=1
            self.printLog('\r#FAS','Generation of %s PPI fasta files from %s hubs complete.' % (rje.integerString(fx),rje.integerString(htot)))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: Yeast Class                                                                                      #
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
    try:Yeast(mainlog,cmd_list).run()

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
