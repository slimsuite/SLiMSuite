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
Module:       rje_ingolia
Description:  Custom script for Ingolia SI compilation
Version:      0.0
Last Edit:    24/03/12
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genemap, rje_obj, rje_seqlist, rje_sequence, rje_zen
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
    (program, version, last_edit, copyright) = ('RJE_INGOLIA', '0.0', 'March 2012', '2012')
    description = 'Custom script for Ingolia SI compilation'
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
### SECTION II: Ingolia Class                                                                                           #
#########################################################################################################################
class Ingolia(rje_obj.RJE_Object):     
    '''
    Ingolia Class. Author: Rich Edwards (2012).

    Str:str
    
    Bool:boolean

    Int:integer

    Num:float
    
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
        self.strlist = []
        self.boollist = []
        self.intlist = []
        self.numlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({})
        self.setBool({})
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
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'inr',['Att'])   # Integers
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
            self.seqSubset2()
            #self.patisMap()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                                |0.0|
        '''Main class setup method.'''
        try:### ~ [1] Load tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.baseFile('ingolia')
            names = {'s1a':'synthesis','s1c':'efficiency',  # Not sure what Table S1 data is!
                     's2a':'genes','s2b':'pauses','s2c':'stops','s3':'starts','s4':'lincRNA'}
            for tab in ['s1a','s1b','s1c','s1d','s2a','s2b','s2c','s3','s4']:
                if tab in names: name = names[tab]
                else: name = tab
                dfile = '%s.%s.csv' % (self.baseFile(),tab)
                if name == 'pauses': self.db().addTable(dfile,mainkeys=['UCSC ID','Codon'],name=name)
                elif name == 'starts': self.db().addTable(dfile,mainkeys=['knownGene','Init Codon [nt]'],name=name)
                else: self.db().addTable(dfile,name=name)
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def seqSubset(self):    ### Extracts sequence subset from MOUSE cDNA and Peptide libraries
        '''Extracts sequence subset from MOUSE cDNA and Peptide libraries.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists('%s.map.tdt' % self.baseFile()):
                mdb = self.db().addTable('%s.map.tdt' % self.baseFile(),mainkeys=['Ingolia'],name='map')
            else:
                ### ~ [2] Load Mouse Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                xfile = '../../../../../Databases/DBase_120225/MGI/mousemap.120324.data.tdt'
                xref = db.addTable(xfile,mainkeys=['Gene'],name='xref')
                afile = '../../../../../Databases/DBase_120225/MGI/mousemap.120324.alias.tdt'
                self.obj['Map'] = rje_genemap.GeneMap(self.log,self.cmd_list)
                #self.obj['Map'].loadPickle('../../../../../Databases/DBase_120225/MGI/mousemap.120324.pickle')
                self.obj['Map'].loadData(['sourcedata=%s' % xfile,'aliases=%s' % afile])
                ing_genes = string.split(string.join(self.db('starts').index('Gene').keys()).upper())
                map = self.obj['Map']
                ing_map = {}
                for gene in ing_genes: ing_map[gene] = map.bestMap(gene)
                ing_mgi = rje.sortUnique(ing_map.values())
                self.printLog('#MUSG','%s Ingolia genes mapped onto %s MGI genes' % (rje.iLen(ing_genes),rje.iLen(ing_mgi)))
                xdb = self.db('xref')
                bad_genes = []
                for gene in ing_mgi[0:]:
                    if gene not in xdb.data():
                        self.printLog('#MAP','Cannot map gene "%s" from Ingolia data!' % gene)
                        bad_genes.append(gene); ing_mgi.remove(gene)
                self.printLog('#BAD','Failed to map %s genes from Ignolia' % rje.iLen(bad_genes))
                open('ingolia.bad.txt','w').write(string.join(bad_genes))
                ### ~ [2] EnsEMBL subset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                ing_musg = xdb.dataList(xdb.entryList(ing_mgi),'EnsEMBL',sortunique=True)
                if '' in ing_musg: ing_musg.remove('')
                self.printLog('#MUSG','%s Ingolia genes mapped onto %s EnsEMBL genes' % (rje.iLen(ing_genes),rje.iLen(ing_musg)))
                if not ing_musg: raise ValueError
                self.deBug(ing_musg[:10])
                for stype in ['cdna','pep']:
                    seqfile = '../MOUSE/Mus_musculus.NCBIM37.66.%s.all.fa' % stype
                    if self.getBool('Force') or not os.path.exists(seqfile):
                        seqout = 'Ingolia.%s.all.fa' % stype
                        seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'seqout=%s' % seqout,'autofilter=T','autload=T','seqmode=file','gooddesc=%s' % string.join(ing_musg,',')]
                        rje_seqlist.SeqList(self.log,seqcmd)
                mdb = self.db().addEmptyTable('map',['Ingolia','Gene','EnsEMBL'],['Ignolia'])
                for gene in ing_map:
                    entry = {'Ingolia':gene,'Gene':ing_map[gene]}
                    if entry['Gene'] in bad_genes: entry['EnsEMBL'] = ''
                    else: entry['EnsEMBL'] = xdb.data()[ing_map[gene]]['EnsEMBL']
                    mdb.addEntry(entry)
            seqfile = 'Ingolia.cdna.all.fa'
            seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'autofilter=F','autload=T','seqmode=file']
            iseq = rje_seqlist.SeqList(self.log,seqcmd)
            if 'ENST' not in mdb.fields():
                mdb.addField('ENST',evalue='')
                while iseq.nextSeq():
                    (iname,icdna) = iseq.getSeq()
                    musg = rje.matchExp('gene:(\S+)',iname)[0]
                    for entry in mdb.indexEntries('EnsEMBL',musg):
                        if entry['ENST']: entry['ENST'] += ',%s' % string.split(iname)[0]
                        else: entry['ENST'] = string.split(iname)[0]
                mdb.saveToFile()
            ### ~ [3] Generate new start sites from Ignolia Harrington data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('starts')
            sdb.dataFormat({'Init Codon [nt]':'int'})
            icod = 'Init Codon [nt]'
            icon = 'Init Context [-3 to +4]'
            sdb.info['Name'] = 'mapped_start'
            sdb.addField('ENST'); sdb.addField('ENSP');
            ENST = open('Ingolia_harr.cdna.all.fa','w')
            ENSP = open('Ingolia_harr.pep.all.fa','w')
            ex = 0.0; etot = sdb.entryNum(); sx = 0; fx = 0
            minpep = 20
            for entry in sdb.entries():
                self.progLog('\r#ING','Mapping Ignolia Harrington Starts: %.2f%%' % (ex/etot)); ex += 100.0
                self.deBug(entry)
                entry[icon] = entry[icon].upper()
                gene = entry['Gene'].upper()
                mentry = mdb.data(gene)
                entry['ENST'] = entry['ENSI'] = ''
                cdnaseq = peptseq = ''
                if not mentry or not mentry['ENST']: fx += 1; continue
                self.deBug(mentry)
                mtype = 'fail'
                for trans in string.split(mentry['ENST'],','):
                    (tname,tseq) = iseq.getDictSeq(trans,format='tuple')
                    if tseq[entry[icod]-4:][:7] == entry[icon]:
                        ipept = string.split(rje_sequence.dna2prot(tseq[entry[icod]-1:]),'*')[0]
                        if len(ipept) > len(peptseq):
                            entry['ENST'] = trans
                            cdnaseq = tseq
                            peptseq = ipept
                            mtype = 'exact'
                if not entry['ENST']:
                    for trans in string.split(mentry['ENST'],','):
                        (tname,tseq) = iseq.getDictSeq(trans,format='tuple')
                        if entry[icon] in tseq:
                            i = tseq.find(entry[icon])
                            ipept = string.split(rje_sequence.dna2prot(tseq[i+3:]),'*')[0]
                            self.deBug('%s: %s' % (entry[icon],ipept))
                            if len(ipept) > len(peptseq):
                                entry['ENST'] = trans; 
                                cdnaseq = tseq
                                peptseq = string.split(rje_sequence.dna2prot(tseq[i+3:]),'*')[0]
                                mtype = 'find'
                if not entry['ENST']:
                    self.printLog('\r#ING','Unable to find Harrington start for %s %s (%s)' % (gene,entry[icod],entry[icon]))
                    fx += 1; continue
                elif len(peptseq) < minpep:
                    self.printLog('\r#ING','Peptide from mapped Harrington start for %s %s (%s) too short!' % (gene,entry[icod],entry[icon]),screen=False)
                    fx += 1; continue
                id = rje.preZero(int(ex/100),etot)
                ENST.write('>ENSINGT%s mtype:%s enst:%s gene:%s ingolia:%s mgi:%s\n%s\n' % (id,mtype,entry['ENST'],mentry['EnsEMBL'],entry['Gene'],mentry['Gene'],cdnaseq))
                ENSP.write('>ENSINGP%s mtype:%s enst:%s gene:%s transcript:ENSINGT%s ingolia:%s mgi:%s\n%s\n' % (id,mtype,entry['ENST'],mentry['EnsEMBL'],id,entry['Gene'],mentry['Gene'],peptseq))
                sx += 1
            sdb.saveToFile('%s.mapped_start.tdt' % self.baseFile())
            ENST.close(); ENSP.close()
            self.printLog('\r#ING','Output %s Ingolia peptides and transcripts. %s failed.' % (rje.iStr(sx),rje.iStr(fx)))
            return
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################
    def seqSubset2(self):    ### Extracts sequence subset from MOUSE cDNA and Peptide libraries
        '''Extracts sequence subset from MOUSE cDNA and Peptide libraries.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists('%s.map.tdt' % self.baseFile()):
                mdb = self.db().addTable('%s.map.tdt' % self.baseFile(),mainkeys=['Ingolia'],name='map')
            else:
                ### ~ [2] Load Mouse Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                xfile = '../../../../../Databases/DBase_120225/MGI/mousemap.120324.data.tdt'
                xref = db.addTable(xfile,mainkeys=['Gene'],name='xref')
                afile = '../../../../../Databases/DBase_120225/MGI/mousemap.120324.alias.tdt'
                self.obj['Map'] = rje_genemap.GeneMap(self.log,self.cmd_list)
                #self.obj['Map'].loadPickle('../../../../../Databases/DBase_120225/MGI/mousemap.120324.pickle')
                self.obj['Map'].loadData(['sourcedata=%s' % xfile,'aliases=%s' % afile])
                ing_genes = string.split(string.join(self.db('starts').index('Gene').keys()).upper())
                map = self.obj['Map']
                ing_map = {}
                for gene in ing_genes: ing_map[gene] = map.bestMap(gene)
                ing_mgi = rje.sortUnique(ing_map.values())
                self.printLog('#MUSG','%s Ingolia genes mapped onto %s MGI genes' % (rje.iLen(ing_genes),rje.iLen(ing_mgi)))
                xdb = self.db('xref')
                bad_genes = []
                for gene in ing_mgi[0:]:
                    if gene not in xdb.data():
                        self.printLog('#MAP','Cannot map gene "%s" from Ingolia data!' % gene)
                        bad_genes.append(gene); ing_mgi.remove(gene)
                self.printLog('#BAD','Failed to map %s genes from Ignolia' % rje.iLen(bad_genes))
                open('ingolia.bad.txt','w').write(string.join(bad_genes))
                ### ~ [2] EnsEMBL subset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                ing_musg = xdb.dataList(xdb.entryList(ing_mgi),'EnsEMBL',sortunique=True)
                if '' in ing_musg: ing_musg.remove('')
                self.printLog('#MUSG','%s Ingolia genes mapped onto %s EnsEMBL genes' % (rje.iLen(ing_genes),rje.iLen(ing_musg)))
                if not ing_musg: raise ValueError
                self.deBug(ing_musg[:10])
                for stype in ['cdna','pep']:
                    seqfile = '../MOUSE/Mus_musculus.NCBIM37.66.%s.all.fa' % stype
                    if self.getBool('Force') or not os.path.exists(seqfile):
                        seqout = 'Ingolia.%s.all.fa' % stype
                        seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'seqout=%s' % seqout,'autofilter=T','autload=T','seqmode=file','gooddesc=%s' % string.join(ing_musg,',')]
                        rje_seqlist.SeqList(self.log,seqcmd)
                mdb = self.db().addEmptyTable('map',['Ingolia','Gene','EnsEMBL'],['Ignolia'])
                for gene in ing_map:
                    entry = {'Ingolia':gene,'Gene':ing_map[gene]}
                    if entry['Gene'] in bad_genes: entry['EnsEMBL'] = ''
                    else: entry['EnsEMBL'] = xdb.data()[ing_map[gene]]['EnsEMBL']
                    mdb.addEntry(entry)
            seqfile = 'Ingolia.cdna.all.fa'
            seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'autofilter=F','autload=T','seqmode=file']
            iseq = rje_seqlist.SeqList(self.log,seqcmd)
            if 'ENST' not in mdb.fields():
                mdb.addField('ENST',evalue='')
                while iseq.nextSeq():
                    (iname,icdna) = iseq.getSeq()
                    musg = rje.matchExp('gene:(\S+)',iname)[0]
                    for entry in mdb.indexEntries('EnsEMBL',musg):
                        if entry['ENST']: entry['ENST'] += ',%s' % string.split(iname)[0]
                        else: entry['ENST'] = string.split(iname)[0]
                mdb.saveToFile()
            ### ~ [3] Generate new start sites from Ignolia Harrington data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db('starts')
            sdb.dataFormat({'Init Codon [nt]':'int'})
            icod = 'Init Codon [nt]'
            icon = 'Init Context [-3 to +4]'
            sdb.info['Name'] = 'mapped_start'
            sdb.addField('ENST'); sdb.addField('ENSP'); sdb.addField('ENSI');
            ENST = open('IngExact.cdna.all.fa','w')
            ENSP = open('IngExact.pep.all.fa','w')
            ex = 0.0; etot = sdb.entryNum(); sx = 0; fx = 0
            minpep = 20
            for entry in sdb.entries():
                self.progLog('\r#ING','Mapping Ignolia Harrington Starts: %.2f%%' % (ex/etot)); ex += 100.0
                #self.deBug(entry)
                entry[icon] = entry[icon].upper()
                gene = entry['Gene'].upper()
                mentry = mdb.data(gene)
                entry['ENST'] = entry['ENSI'] = ''
                cdnaseq = peptseq = ''
                if not mentry or not mentry['ENST']: fx += 1; continue
                #self.deBug(mentry)
                mtype = 'fail'
                for trans in string.split(mentry['ENST'],','):
                    (tname,tseq) = iseq.getDictSeq(trans,format='tuple')
                    self.deBug('%s vs %s' % (tseq[entry[icod]-3:][:7],entry[icon]))
                    if tseq[entry[icod]-3:][:7] == entry[icon]:
                        ipept = string.split(rje_sequence.dna2prot(tseq[entry[icod]:]),'*')[0]
                        self.deBug(ipept)
                        if len(ipept) > len(peptseq):
                            entry['ENST'] = trans
                            cdnaseq = tseq
                            peptseq = ipept
                            mtype = 'exact'
                if not entry['ENST']:
                    self.printLog('\r#ING','Unable to find Harrington start for %s %s (%s)' % (gene,entry[icod],entry[icon]),screen=False)
                    fx += 1; continue
                elif len(peptseq) < minpep:
                    self.printLog('\r#ING','Peptide from mapped Harrington start for %s %s (%s) too short!' % (gene,entry[icod],entry[icon]),screen=False)
                    fx += 1; continue
                id = rje.preZero(int(ex/100),etot)
                entry['ENSI'] = 'ENSINGT%s' % id
                entry['ENSP'] = 'ENSINGP%s' % id
                ENST.write('>ENSINGT%s mtype:%s enst:%s gene:%s ingolia:%s mgi:%s\n%s\n' % (id,mtype,entry['ENST'],mentry['EnsEMBL'],entry['Gene'],mentry['Gene'],cdnaseq))
                ENSP.write('>ENSINGP%s mtype:%s enst:%s gene:%s transcript:ENSINGT%s ingolia:%s mgi:%s\n%s\n' % (id,mtype,entry['ENST'],mentry['EnsEMBL'],id,entry['Gene'],mentry['Gene'],peptseq))
                sx += 1
            sdb.saveToFile('%s.mapped_exact.tdt' % self.baseFile())
            ENST.close(); ENSP.close()
            self.printLog('\r#ING','Output %s Ingolia peptides and transcripts. %s failed.' % (rje.iStr(sx),rje.iStr(fx)))
            return
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################
    def patisMap(self): ### Maps PATIS run onto Ingloia starts
        '''Maps PATIS run onto Ingloia starts.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #pdb = self.db().addTable('../../2011-09-BIOL3050-11/Mus_musculus.NCBIM37.64.patis.tdt',mainkeys=['ENST'],name='patis')
            pdb = self.db().addTable('../Ingolia.gatic.aic.tdt',mainkeys=['transcript'],name='patis')
            #for entry in pdb.entries():
            #    for tups in [('Core','Strength'),('eCore','eType'),('tCore','tType')]:
            #        core = entry[tups[0]]
            #        sfield = tups[1]
            #        entry[sfield] = 'Strong'
            #        if not core: entry[sfield] = 'None'
            #        elif core[3:-1] != 'ATG': entry[sfield] = 'NonAUG'
            #        elif core[0] not in 'GA' and core[-1] != 'G': entry[sfield] = 'Weak'
            #        elif core[0] and core[-1] != 'G': entry[sfield] = 'MidR'
            #        elif core[0] not in 'GA': entry[sfield] = 'MidG'
            mdb = self.db().addTable(mainkeys=['ENST','Init Codon [nt]'],name='mapped_start')
            pdb.renameField('transcript','ENST')
            jdb = self.db().joinTables(name='patis_map',join=[(pdb,'ENST'),(mdb,'ENST')],newkey=['ENST','Init Codon [nt]'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
            jdb.addField('PType'); jdb.dropField('ENSP')
            for entry in jdb.entries():
                if not entry['Init Context [-3 to +4]']: entry['PType'] = ''
                elif entry['eContext'] == entry['Init Context [-3 to +4]']:
                    entry['PType'] = 'eORF'
                    if int(entry['Init Codon [nt]']) - int(entry['eORF']) in [-1]: entry['PType'] = 'eORF+'
                elif entry['tContext'] == entry['Init Context [-3 to +4]']:
                    entry['PType'] = 'tORF'
                    if int(entry['Init Codon [nt]']) - int(entry['tORF']) in [-1]: entry['PType'] = 'tORF+'
                elif entry['context'] == entry['Init Context [-3 to +4]']:
                    entry['PType'] = 'Annotated'
                    if int(entry['Init Codon [nt]']) - int(entry['start']) in [-1]: entry['PType'] = 'Annotated+'
                else: entry['PType'] = '?'
            jdb.saveToFile()
            ### Annotated ###
            adb = self.db().copyTable(jdb,'annotated')
            adb.dropEntriesDirect('PType',''); myrules = {}
            for field in ['All','Annotated','Annotated+']: adb.addField(field,evalue=0); myrules[field] = 'sum'
            for entry in adb.entries():
                entry['All'] = 1
                if 'Annotated' in entry['PType']: entry['Annotated'] = 1
                if 'Annotated+' in entry['PType']: entry['Annotated+'] = 1
            adb.dropFields(['ENST','Init Codon [nt]','strength','All','Annotated','Annotated+'],inverse=True)
            adb.compress(['ENST','strength'],default='max'); adb.dropField('Init Codon [nt]')
            adb.compress(['strength'],rules=myrules,default='sum'); adb.dropField('ENST')
            for entry in adb.entries(): self.deBug(entry)
            adb.saveToFile()
            ### Ratings ###
            jdb.info['Name'] = 'patis_rating'
            jdb.addField('Count',evalue=1)
            jdb.dropEntriesDirect('PType','')
            jdb.dropFields(['ENST','Init Codon [nt]','Product','PType','Count'],inverse=True)
            jdb.compress(['Product','PType'],default='sum'); jdb.dropFields(['Init Codon [nt]','ENST'])
            jdb.saveToFile()
            
            
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################
### End of SECTION II: New Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
class GeneMap(rje_genemap.GeneMap):
    '''Just for pickling.'''
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

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
    try: Ingolia(mainlog,cmd_list).run()

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
