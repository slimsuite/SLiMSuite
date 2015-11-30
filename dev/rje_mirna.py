#!/usr/local/bin/python

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
Module:       rje_miRNA
Description:  Miscellaneous miRNA Analysis Module
Version:      0.3
Last Edit:    15/07/13
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module controls various custom miRNA analyses. It is envisioned that it will be replaced in time with other
    modules.

Commandline:
    ### ~ INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE  : Main file for analysis, containing UTRs []
    maxseq=X    : Maximum number of sequences to process [500]
    seedfreq=T/F: Whether to use overall Seed Frequencies per sequence for slimprob [False]
    
    ### ~ ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    run=LIST    : List of analyses/manipulations (nrutr/enspita/cleal) []
    enspath=PATH: Path to EnsLoci files [/scratch/Databases/NewDB/EnsEMBL/]

    ### ~ CLEAL Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minmax=X        : Minimum value for Maximum expression value [1.0]
    arraydata=FILE  : File of array data for Cleal analysis []
    miranda=FILE    : Miranda prediction file []   
    validated=FILE  : Validate mRNA target file []   

Uses general modules: glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_seq, rje_slim, rje_zen, slimfinder
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, sets, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_slim, rje_zen
import slimfinder
#import rje_blast, rje_sequence, rje_scoring, rje_xgmml
#import rje_slimcalc, rje_slimcore, rje_slimlist
#import rje_motif_V3 as rje_motif            # Used for expect method only 
#import rje_dismatrix_V2 as rje_dismatrix
#import comparimotif_V3 as comparimotif
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Shifted focus of analysis to use PITA instead.
    # 0.2 - Added Cleal Analysis of Placental versus miRanda downloads (http://www.microrna.org/microrna/getDownloads.do)
    # 0.3 - Tidied up some import calls.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [X] : Calculate observed vs. expected for 8mers - seem to be many over-represented.
    # [ ] : Needs full testing and documented functionality.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_miRNA', '0.3', 'July 2013', '2008')
    description = 'Miscellaneous miRNA Analysis Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: miRNA Class                                                                                             #
#########################################################################################################################
class miRNA(rje.RJE_Object):     
    '''
    miRNA Class. Author: Rich Edwards (2008).

    Info:str
    - ArrayData = File of array data for Cleal analysis []
    - EnsPath = Path to EnsLoci files [/scratch/Databases/NewDB/EnsEMBL/]
    - Miranda = Miranda prediction file []   
    - Validated = Validate mRNA target file []   
    
    Opt:boolean
    - SeedFreq = Whether to use overall Seed Frequencies per sequence for slimprob [False]

    Stat:numeric
    - MaxSeq = Maximum number of sequences to process [500]
    - MinMax = Minimum value for Maximum expression value [1.0]

    List:list
    - Run = List of analyses/manipulations (nrutr) []    

    Dict:dictionary    

    Obj:RJE_Objects
    - SeqList = UTR Sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['EnsPath','Seqin','Validated','ArrayData','Miranda']
        self.optlist = ['SeedFreq']
        self.statlist = ['MaxSeq','MinMax']
        self.listlist = ['Run']
        self.dictlist = []
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'EnsPath':rje.makePath('/scratch/Databases/NewDB/EnsEMBL/')})
        self.setStat({'MaxSeq':500,'MinMax':1.0})
        ### Other Attributes ###
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
                self._cmdReadList(cmd,'file',['Seqin','Validated','ArrayData','Miranda'])
                self._cmdReadList(cmd,'path',['EnsPath'])
                self._cmdReadList(cmd,'int',['MaxSeq'])
                self._cmdReadList(cmd,'stat',['MinMax'])
                self._cmdReadList(cmd,'opt',['SeedFreq'])
                self._cmdReadList(cmd,'list',['Run'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def run(self):  ### Main class run method.
        '''Main class run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Run'] = string.split(string.join(self.list['Run']).lower())
            if 'nrutr' in self.list['Run']:
                self.obj['SeqList'] = seqlist = rje_seq.SeqList(self.log,['accnr=F','seqnr=F','gnspacc=F','dna=T']+self.cmd_list) # Load UTR sequences
                self.nrUTR(seqlist)     # Convert UTRs into single UTR per gene
            if 'enspita' in self.list['Run']: self.ensPITA()
            if 'cleal' in self.list['Run']: self.cleal()
                
            ### ~ [X] Old Stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            if not self.list['Run'] and rje.yesNo('Load UTR sequences?'):
                self.obj['SeqList'] = seqlist = rje_seq.SeqList(self.log,self.cmd_list+['accnr=F','seqnr=F','gnspacc=F','dna=T']) # Load UTR sequences
                if rje.yesNo('NR UTR?'): self.nrUTR(seqlist)     # Convert UTRs into single UTR per gene
                if rje.yesNo('Generate seed sequences?'): self.seedSeq(seqlist)     # Generate seed sequence lists
                if rje.yesNo('Assess 8mer representation?'): self.xmer(seqlist)     # Assess xmers occurrences using AAFreq
            if not self.list['Run'] and rje.yesNo('Run SLiMFinder?'): self.slimFinder()
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <2> ### Data Preparation Methods                                                                                #
#########################################################################################################################
    def nrUTR(self,seqlist):    ### Reduces and reformats UTR sequences. (Must be EnsEMBL download.)
        '''Reduces and reformats UTR sequences. (Must be EnsEMBL download.)'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            oldname = seqlist.info['Name']
            newname = string.replace(seqlist.info['Name'],'.fasta','.fas')
            if oldname == newname and rje.yesNo('Sequence file %s already reformatted?' % oldname):
                return self.printLog('#SEQ','Sequence file %s already reformatted' % oldname)
            #>ENSG00000177757|1|protein_coding|ENST00000326734|Description|742614|745077

            ### ~ [2] ~ Reduce UTRs to one per gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            utr = {}   # Dictionary of {EnsGene:Sequence object}
            for seq in seqlist.seq:
                if string.count(seq.info['Name'],'|') >= 5:
                    try: [gene,csome,type,tran,desc,start,end] = string.split(seq.info['Name'],'|')
                    except: [gene,csome,type,tran,start,end,desc] = string.split(seq.info['Name'],'|') + ['']
                elif string.count(seq.info['Name'],'|') == 2:
                    [gene,type,tran,pep,desc] = string.split(seq.info['Name'],'|') + ['','']
                else:
                    try: [gene,type,tran,pep,desc] = string.split(seq.info['Name'],'|')
                    except:
                        try: [gene,type,tran,pep,desc] = string.split(seq.info['Name'],'|') + ['']
                        except: self.errorLog(seq.info['Name'])
                if gene in utr and utr[gene].aaLen() >= seq.aaLen(): continue
                if type != 'protein_coding': continue   #!# protein_coding Only for now! #!#
                utr[gene] = seq
                if string.count(seq.info['Name'],'|') >= 5: seq.info['Name'] = '%s [%s; %s] %s:%s-%s %s' % (gene,tran,type,csome,start,end,desc)
                else: seq.info['Name'] = '%s [%s; %s; %s] %s' % (gene,tran,pep,type,desc)

            ### ~ [3] ~ Reduce seqlist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist.info['Name'] = newname
            seqlist.seq = []
            for gene in rje.sortKeys(utr): seqlist.seq.append(utr[gene])
            seqlist.saveFasta()
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise
#########################################################################################################################
    def ensPITA(self):  ### Generates data for running PITA on EnsEMBL data
        '''Generates data for running PITA on EnsEMBL data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ufiles = glob.glob('*UTR.fas')      # Files to process
            open('mir.bat','w')                 # Batch file run
            ### ~ [2] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for u in ufiles:
                spec = u[:5]
                ## ~ [2a] ~ Load sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ensloci = rje_seq.SeqList(self.log,['accnr=F','seqnr=F','gnspacc=F']+self.cmd_list+['seqin=%sens_%s.loci.fas' % (self.info['EnsPath'],spec)])
                eseq = {}
                for seq in ensloci.seqs()[0:]:
                    try: gene = rje.matchExp('gene:(\S+)\]',seq.info['Name'])[0]
                    except: ensloci.seq.remove(seq); continue
                    eseq[gene] = seq
                    seq.info['Name'] = '%s %s' % (gene,seq.info['Name'])
                utrs = rje_seq.SeqList(self.log,['accnr=F','seqnr=F','gnspacc=F','dna=T']+self.cmd_list+['seqin=%s' % u])
                useq = {}
                for seq in utrs.seqs(): useq[seq.shortName()] = seq
                ## ~ [2b] ~ Make new files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                genes = list(sets.Set(eseq.keys()).intersection(sets.Set(useq.keys())))
                self.printLog('#GENE','%s => %s common EnsLoci/UTR genes' % (spec,len(genes)))
                if not genes:
                    print rje.sortKeys(useq)[:100]
                    print rje.sortKeys(eseq)[:100]
                    continue
                for ens in eseq:
                    if ens not in genes: ensloci.seq.remove(eseq[ens])
                ensloci.info['Name'] = 'mir_%s.pep.fas' % spec
                ensloci.saveFasta()
                for utr in useq:
                    if utr not in genes: utrs.seq.remove(useq[utr])
                utrs.info['Name'] = 'mir_%s.utr.fas' % spec
                utrs.saveFasta()
                ## ~ [2c] ~ Add batch commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try: mir = {'CAEEL':'worm','DROME':'fly'}[spec]
                except: mir = spec.lower()
                open('mir.bat','a').write('perl /home/re1u06/Bioware/PITA/pita_prediction.pl -utr %s -mir /home/re1u06/Bioware/PITA/known_mirs/%s_mirs.fasta -prefix mir_%s\n' % (utrs.info['Name'],mir,spec))
                open('mir.bat','a').write('python /home/re1u06/Serpentry/gablam.py gnspacc=F seqin=%s log=mir_%s\n' % (ensloci.info['Name'],spec))
            ### ~ [3] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for b in open('mir.bat','r').readlines(): self.printLog('#BAT',b)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def newStuff(self): ### Make something out of this
        '''Make something out of this.'''
        try: return
        ### TRANSLATIONS ###
        # AccNum :: FBtr0071764
        # Part of Description :: length=XXXX; parent=FBgn0000008; => GeneID
        ### 3' UTRs ###
        # ShortName: FBtr0071764
        # Part of Description :: length=XXXX; parent=FBgn0000008; 
        ### EnsEMBL ###
        # Can the same be done for EnsEMBL sequences using EnsLoci: add UTRs to downloads and process as in nrUTR?
        # *** Try to download old versions of UTRs *** #
        # Caenorhabditis_elegans  CAEEL   WS180.48        26903   20140   0       0       1736    0       14330   4074    0
        # Drosophila_melanogaster DROME   BDGP4.3.48      8127    4713    11664   9326    2389    0       11641   6       3
        # Homo_sapiens    HUMAN   NCBI36.48       43371   21798   3222    2753    17206   0       3011    1581    1199
        # Mus_musculus    MOUSE   NCBIM37.48
        # ... Try to use nrUTR above in conjunction with EnsLoci sequences <<<
        
        # 1. Reduce to longest translation and longest UTR for each gene.
        # 2. GABLAM translations and use to estimate total UPC
        # 3. Run PITA on UTRs :: perl /home/re1u06/Bioware/PITA/pita_prediction.pl -utr FILE -mir /home/re1u06/Bioware/PITA/known_mirs/fly_mirs.fasta -prefix X
        # ...  Also  human_mirs.fasta  mouse_mirs.fasta  worm_mirs.fasta
        # ..? Should PITA be run on one sequence at a time and results parsed as we go?
        # ..? Or run one miRNA at a time?
        # 4. Use translations to generate UPC for each mir
        # 5. Look at overlap between mirs and assess using 2.
        except: self.errorLog('Shucks')
#########################################################################################################################
    ### <3> ### Jane Cleal custom analysis pipeline                                                                     #
#########################################################################################################################
    def cleal(self):    ### Custom analysis pipeline for mapping array data onto validated and predicted targets
        '''Custom analysis pipeline for mapping array data onto validated and predicted targets.'''
        try:### ~ [0] ~ Set up database object and load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.basefile('cleal_mirna')
            vdb = db.addTable(self.info['Validated'],mainkeys=['mir','gene'],name='ValidatedIn')
            #mir	gene	NM
            #143	MAPK7	NM_139034
            for entry in vdb.entries(): entry['mir'] = 'HSA-MIR-%s' % entry['mir'].upper()
            vdb.renameField('mir','mirna_name')
            vdb.renameField('gene','Gene')
            vdb.list['Fields'] = ['Gene','mirna_name','NM']
            adb = db.addTable(self.info['ArrayData'],mainkeys=['gene_id'],name='ArrayData')
            #gene_id	miRNA	mature_acc	gene_acc	219_pl	112_pl	170_pl
            #35816	LET-7B: HSA-LET-7B; MMU-LET-7B; RNO-LET-7B	MIMAT0000775 MIMAT0000522 MIMAT0000063	MI0000063 MI0000558 MI0000829	0.34	0.31	0.22
            adb.renameField('miRNA','mirna_name')
            for entry in adb.entries():
                if entry['mirna_name'][:3] in ['LET','MIR']: entry['mirna_name'] = 'HSA-%s' % (string.split(entry['mirna_name'].upper(),':')[0])
                else: entry['mirna_name'] = string.split(entry['mirna_name'].upper(),':')[0]
            adb.newKey(['mirna_name'])
            adb.dropFields(['gene_id','mature_acc','gene_acc',''])
            adb.dataFormat({'219_pl':'float','112_pl':'float','170_pl':'float'})
            adb.fillBlanks(0.0,fillempty=True)
            adb.addField('max_pl'); adb.addField('mean_pl')
            self.deBug(adb.fields())
            for entry in adb.entries():
                entry['max_pl'] = max(entry['219_pl'],entry['112_pl'],entry['170_pl'])
                entry['mean_pl'] = sum([entry['219_pl'],entry['112_pl'],entry['170_pl']])/3.0
            adb.dropEntries(['max_pl==0'])
            adb.dropEntries(['max_pl<%f' % self.stat['MinMax']],logtxt='MinMax expression filter')
            adb.saveToFile(); #vdb.saveToFile(); mdb.saveToFile()
            mdb = db.addTable(self.info['Miranda'],mainkeys=['#mirbase_acc','gene_symbol'],name='MirandaIn')
            #mirbase_acc	mirna_name	gene_id	gene_symbol	transcript_id	ext_transcript_id	mirna_alignment	alignment	gene_alignment	mirna_start	mirna_end	gene_start	gene_end	genome_coordinates	conservation	align_score	seed_cat	energy	mirsvr_score
            #MIMAT0000062	hsa-let-7a	5270	SERPINE2	uc002vnu.2	NM_006216	uuGAUAUGUUGGAUGAU-GGAGu	  | :|: ||:|| ||| |||| 	aaCGGUGAAAUCU-CUAGCCUCu	2	21	495	516	[hg19:2:224840068-224840089:-]	0.5684	122	0	-14.73	-0.7269
            for entry in mdb.entries(): entry['mirna_name'] = entry['mirna_name'].upper()
            mdb.renameField('gene_symbol','Gene')
            mdb.newKey(['Gene','mirna_name'],True)
            ### ~ [1] ~ Compile new tables for HAPPI Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vdb2 = db.joinTables(name='Validated',join=[(vdb,'mirna_name'),(adb,'mirna_name')],newkey=['Gene','mirna_name'],empties=False,keeptable=True)
            vdb2.saveToFile()
            mdb2 = db.joinTables(name='miRanda',join=[(mdb,'mirna_name'),(adb,'mirna_name')],newkey=['Gene','mirna_name'],empties=False,keeptable=True)
            mdb2.saveToFile()
            ### ~ [2] ~ Generate Gene-Centred Summary Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = db.copyTable('Validated','Genes')
            gdb.dropFields(['Gene','mirna_name','219_pl','112_pl','170_pl','max_pl','mean_pl'],inverse=True)
            gdb.compress(['Gene'],rules={'max_pl':'max','mirna_name':'list'},default='sum')
            gdb.addField('Class',evalue='Validated')
            gdb.newKey(['Gene','Class'])
            mgdb = db.copyTable('miRanda','MGene')
            mgdb.dropFields(['Gene','mirna_name','219_pl','112_pl','170_pl','max_pl','mean_pl'],inverse=True)
            mgdb.compress(['Gene'],rules={'max_pl':'max','mirna_name':'list'},default='sum')
            mgdb.addField('Class',evalue='miRanda')
            mgdb.newKey(['Gene','Class'])
            db.mergeTables(gdb,mgdb,overwrite=True,matchfields=True)
            gdb.addField('mirna_count',after='mirna_name')
            for entry in gdb.entries(): entry['mirna_count'] = len(string.split(entry['mirna_name'],';'))
            #gdb = db.joinTables(name='Genes',join=[(vgdb,'gene'),(mgdb,'gene')],newkey=['gene'])
            #gdb.fillBlanks()
            self.deBug(gdb.fields())
            #for entry in gdb.entries():
            #    if entry['miRanda_Class'] == 'miRanda':
            #        if entry['Class'] == 'Validated': entry['Class'] = 'Validated-miRanda'
            #        else: entry['Class'] = 'miRanda'
            #gdb.dropFields(['miRanda_Class'])
            gdb.saveToFile()
            #self.deBug(gdb.index('mirna_count').keys())
            for dkey in gdb.dataKeys()[0:]:
                entry = gdb.data(dkey)
                if entry['Class'] == 'Validated': continue
                if entry['mirna_count'] < 10: gdb.data().pop(dkey)
            gdb.saveToFile('%s.Genes10.tdt' % self.basefile())
                

        except: self.errorLog('Oh no!')        
#########################################################################################################################
    ### <X> ### Old Co-occurrence analysis Methods                                                                      #
#########################################################################################################################
    def xmer(self,seqlist): ### Use AAFreq to assess general over/under-representation of 8mers
        '''Use AAFreq to assess general over/under-representation of 8mers.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nfreq = seqlist.aaFreq(alphabet=['A','C','G','T'],fromfile=self.info['Seqin'],total=True)
            total = nfreq.pop('Total')  # Total no. bps
            snum = seqlist.seqNum()     # No. sequences
            slen = float(total) / snum  # Mean sequence length
            xfile = 'mirna.seedprob.tdt'
            head = ['Seed','Count','Exp','OverP','OverC','UnderP','UnderC']
            bonf = 4 ** 8
            seqfiles = glob.glob('miUTR/*.acc')
            seqfiles.sort()
            (sx,stot) = (0,len(seqfiles))
            self.list['CommonSeed'] = []    # List of seeds where Obs > Exp
            if os.path.exists('mirna.commonseed.txt') and rje.yesNo('Use mirna.commonseed.txt data?'):
                self.list['CommonSeed'] = self.loadFromFile('mirna.commonseed.txt',chomplines=True)
                self.printLog('\r#PROB','%s over-represented seeds.' % rje.integerString(len(self.list['CommonSeed'])))
                return
            rje.delimitedFileOutput(self,xfile,head,'\t',rje_backup=True)
            ### ~ [2] ~ Crude method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for file in seqfiles:
                self.printLog('\r#PROB','Calculating crude seed probabilities: %.1f%%' % (sx/stot),newline=False,log=False)
                sx += 100.0
                seed = string.split(os.path.basename(file),'.')[0]  # 8mer seed sequence
                ax = len(open(file,'r').readlines())                # No. of sequences with seed
                p = 1.0
                for a in seed: p *= nfreq[a]
                p = rje.logPoisson(1,slen*p,callobj=self)   # Mean prob of 1+ occurrence per sequence
                ex = p * stot                                # Expected no.
                over = rje.logPoisson(ax,ex,callobj=self)
                under = 1.0 - rje.logPoisson(ax+1,ex,callobj=self)
                overc = rje.binomial(1,bonf,over,callobj=self)
                underc = rje.binomial(1,bonf,under,callobj=self)
                data = {'Seed':seed,'Count':ax,'Exp':'%.1f' % ex,'OverP':over,'UnderP':under,'OverC':overc,'UnderC':underc}
                rje.delimitedFileOutput(self,xfile,head,'\t',data)
                if ax > ex: self.list['CommonSeed'].append(seed)
            self.printLog('\r#PROB','Calculating crude seed probabilities complete: %s over-represented' % rje.integerString(len(self.list['CommonSeed'])))
            open('mirna.commonseed.txt','w').write(string.join(self.list['CommonSeed'],'\n'))
        except: self.log.errorLog('%s: xmer() error.' % rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimFinder(self):   ### Mask miRNA seeds from sequences containing them and runs SLiMFinder.
        '''Mask miRNA seeds from sequences containing them and runs SLiMFinder.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfcmd = 'batch= slimlen=8 maxwild=0 preamb=F dna=T masking=F termini=F absminocc=3 extras=F savespace=2 append=T resfile=mirna_slimfinder.csv runid=miRNA efilter=F'
            sqcmd = string.split(string.replace(string.join(self.cmd_list),'seqin','fasdb')) + ['autoload=T','dna=T']
            seqfiles = glob.glob('miUTR/*.acc')
            seqfiles.sort()
            (sx,stot) = (0,len(seqfiles))
            if not self.obj['SeqList']: self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['accnr=F','seqnr=F','gnspacc=F','dna=T'])
            accdict = self.obj['SeqList'].seqNameDic('AccNum')
            #!# Add skiplist from mirna_slimfinder.csv #!#
            ### ~ [2] Calculate frequencies of each 8mer for revised SLiMFinder slimProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seedfreq = {}
            if os.path.exists('seedfreq.tdt'):
                for line in self.loadFromFile('seedfreq.tdt',chomplines=True)[1:]:
                    [seed,freq] = string.split(line)
                    seedfreq[seed] = string.atof(freq)
            else:
                fx = 0.0
                FREQ = open('seedfreq.tdt','w')
                FREQ.write('Seed     Freq\n')
                for file in seqfiles:
                    self.printLog('\r#FREQ','Calculating seed freq: %.1f%%' % (fx/stot),newline=False,log=False)
                    fx += 100.0
                    seed = string.split(os.path.basename(file),'.')[0]
                    ax = float(len(open(file,'r').readlines()))
                    seedfreq[seed] = ax / stot
                    FREQ.write('%s %s\n' % (seed,seedfreq[seed]))
                self.printLog('\r#FREQ','Calculating seed freq complete: saved to seedfreq.tdt')
                FREQ.close()
                
            ### ~ [3] Process each file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for file in seqfiles:
                sx += 1
                try:
                    ## ~ [2a] ~ Mask out seed sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    seed = string.split(os.path.basename(file),'.')[0]
                    ax = len(open(file,'r').readlines())
                    if ax > self.stat['MaxSeq']:
                        self.printLog('#MAX','miRNA seed %s has %s > %s UTRs: skipped' % (seed,rje.integerString(ax),rje.integerString(self.stat['MaxSeq'])))
                        continue
                    #try: seqlist = rje_seq.SeqList(self.log,sqcmd+['accnr=F','seqnr=F','gnspacc=F','seqin=%s' % file])   # Load UTR sequences
                    #except:
                    seqlist = rje_seq.SeqList(self.log,sqcmd+['accnr=F','seqnr=F','gnspacc=F','seqin=None'])   # Load UTR sequences
                    seqlist.seq = []
                    for acc in self.loadFromFile(file,chomplines=True): seqlist.seq.append(accdict[acc])
                    for seq in seqlist.seq: seq.info['Sequence'] = string.replace(seq.info['Sequence'],seed,'NNNNNNNN')
                    ## ~ [2b] ~ Save file and run SLiMFinder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    seqlist.saveFasta(seqfile='%s.fas' % seed)
                    del seqlist
                    msf = miSLiMFinder(self.log,string.split(sfcmd)+self.cmd_list+['seqin=%s.fas' % seed])
                    msf.opt['SeedFreq'] = self.opt['SeedFreq']
                    try: msf.list['CommonSeed'] = self.list['CommonSeed']
                    except: msf.list['CommonSeed'] = []
                    msf.dict['SeedFreq'] = seedfreq
                    msf.run()
                    os.unlink('%s.fas' % seed)
                    self.printLog('#RUN','%s of %s seed UTR files analysed' % (rje.integerString(sx),rje.integerString(stot)))
                except: self.log.errorLog(seed)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def OLDslimFinder(self):     ### Mask miRNA seeds from sequences containing them and runs SLiMFinder.
        '''Mask miRNA seeds from sequences containing them and runs SLiMFinder.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfcmd = 'batch= slimlen=8 maxwild=0 preamb=F dna=T masking=F termini=F absminocc=3 extras=F savespace=2 append=T resfile=mirna_slimfinder.csv runid=miRNA efilter=F'
            sqcmd = string.split(string.replace(string.join(self.cmd_list),'seqin','fasdb')) + ['autoload=T','dna=T']
            seqfiles = glob.glob('miUTR/*.acc')
            seqfiles.sort()
            (sx,stot) = (0,len(seqfiles))
            #!# Add skiplist from mirna_slimfinder.csv #!#
            ### ~ [2] Process each file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for file in seqfiles:
                sx += 1
                ## ~ [2a] ~ Mask out seed sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seed = string.split(os.path.basename(file),'.')[0]
                ax = len(open(file,'r').readlines())
                if ax > self.stat['MaxSeq']:
                    self.printLog('#MAX','miRNA seed %s has %s > %s UTRs: skipped' % (seed,rje.integerString(ax),rje.integerString(self.stat['MaxSeq'])))
                    continue
                seqlist = rje_seq.SeqList(self.log,sqcmd+['accnr=F','seqnr=F','gnspacc=F','seqin=%s' % file])   # Load UTR sequences
                for seq in seqlist.seq: seq.info['Sequence'] = string.replace(seq.info['Sequence'],seed,'NNNNNNNN')
                ## ~ [2b] ~ Save file and run SLiMFinder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqlist.saveFasta(seqfile='%s.fas' % seed)
                os.system('python /home/re1u06/Serpentry/slimfinder.py %s %s seqin=%s.fas' % (sfcmd,string.join(self.cmd_list),seed))
                os.unlink('%s.fas' % seed)
                self.printLog('#RUN','%s of %s seed UTR files analysed' % (rje.integerString(sx),rje.integerString(stot)))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def seedSeq(self,seqlist):    ### Generate sequence lists for potential seed patterns.
        '''Generate sequence lists for potential seed patterns .'''
        try:### ~ [1] Generate sequence lists for potential seed patterns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seedseq = {'':seqlist.seq[0:]}
            for n in range(8):
                seeds = rje.sortKeys(seedseq)
                (sx,stot) = (0.0,len(seeds))
                for seed in seeds:
                    self.printLog('\r#SEED','Growing seeds: %.2f%% (%s seeds)' % (12.5*n+sx/stot,rje.integerString(len(seedseq))),newline=False,log=False)
                    sx += 12.5
                    baseseq = seedseq.pop(seed)
                    for b in 'ACGT':
                        newseed = seed + b
                        if len(newseed) < 5: newseq = baseseq[0:]     # Assume found in all!
                        else:
                            newseq = []
                            for seq in baseseq:
                                if seq.info['Sequence'].find(newseed) >= 0: newseq.append(seq)
                        if len(newseq) >= 3: seedseq[newseed] = newseq[0:]      # Must be in 3+ UTRs
                        self.deBug('%s > %s >> %d' % (b,newseed,len(newseq)))
            self.printLog('\r#SEED','Growing seeds complete: %s seeds have 3+ UTRs' % (rje.integerString(len(seedseq))))
            ### ~ [2] Output to sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'miUTR/')
            seeds = rje.sortKeys(seedseq)
            (sx,stot) = (0.0,len(seeds))
            for seed in seeds:
                self.printLog('\r#FAS','Outputting acclist files for search: %.1f%%' % (sx/stot),newline=False,log=False)
                sx += 100.0
                seqlist.saveAcc(seqs=seedseq[seed],accfile='miUTR/%s.utr.acc' % seed,log=False)
            self.printLog('\r#FAS','%s acclist files for searching output to miUTR/*.utr.acc' % rje.integerString(stot))                
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: miRNA Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
class miSLiMFinder(slimfinder.SLiMFinder):
    '''
    Modified SLiMFinder class to jump straight to the chase to get 8nt SLiMs.
    '''
#########################################################################################################################
    def run(self):  ### Modified Main SLiMFinder Run Method
        '''
        Main SLiMFinder Run Method:
        1. Input:
            - Read sequences into SeqList
        2. SLiMBuild:
            - Check for existing Pickle and load if found. Check appropriate parameter settings and re-run if desired.
            - or - Save sequences as fasta and mask sequences in SeqList
            -  Perform BLAST and generate UPC based on saved fasta.
            - Calculate AAFreq for each sequence and UPC.
            - Find all dimer motifs in dataset using MinWild/MaxWild parameters.
            - Extend to SLiMs and add ambiguity
        5. Identify significant SLiMs.
        6. Output results and tidy files.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:###~INPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            self.setupBasefile()
            self.loadAADimerFreq()
            ## Check whether to bother running dataset at all - Check Input versus Min and Max Seq ##
            if self.stat['MaxSeq'] > 0 and self.stat['MaxSeq'] < self.obj['SeqList'].seqNum():
                self.log.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.integerString(self.obj['SeqList'].seqNum()),rje.integerString(self.stat['MaxSeq'])))
                return False
            if self.seqNum() < self.stat['MinOcc']:
                self.log.printLog('#SEQ','Insufficient Sequences (%d) for MinOcc setting (%d). Run aborted.' % (self.seqNum(),self.stat['MinOcc']))
                return

            ###~SLiMBuild~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.stat['StartTime'] = time.time()
            ## UPC and MinOcc settings: needed to identify correct pickle so must be done first ##
            if not self.setupFocus(): return    ## Setup Focus before UPC & MST - also check against Pickle ##
            if not self.makeUPC():
                return self.log.errorLog('Error during makeUPC(). Abandoning %s run' % self.dataset(),printerror=False)
            if not self.setupMinOcc(): 
                self.log.printLog('#UPC','Insufficient UPC (%d) for MinOcc setting (%d). Run aborted.' % (self.UPNum(),self.stat['MinOcc']))
                return
            if self.stat['MaxUPC'] >= self.stat['MinOcc'] and self.stat['MaxUPC'] < self.UPNum():
                self.log.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). Run aborted.' % (self.UPNum(),self.stat['MaxUPC']))
                return

            ## Check for existing pickle to replace SLiMBuild portion ##
            self.setupResults()                 ## Sets up OccStats filter etc. - check against Pickle ##
            
            ## Setup Main Results File early in case of user intervention ##
            self.backupOrCreateResFile()

            ## AA Frequency Calculations made early as needed superficially in SLiMBuild ##
            self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
            if self.opt['MaskFreq']: self.makeAAFreq()
            else: 
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                self.makeAAFreq()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
            self.adjustAATotals()

            ## Execute SLiMBuild if pickle not loaded, else recalculate Bonferroni ##
            if self.opt['SlimBuild'] or self.opt['SlimChance'] or self.opt['SlimDisc']:
                self.makeBonferroni()   # Estimates and reports total no. motifs in dataset
                self.dict['Extremf.'][8] = (4 ** 8) - len(self.list['CommonSeed'])  # Only looking at rarer subset
                self.printLog('#SPACE','Reduced 8mers using CommonSeed to %s' % rje.integerString(self.dict['Extremf.'][8]))
                #>>> Replace here <<<#
                self.makeDimers()       # Makes all ai.{0,x}aj dimers
                #X#self.reduceDimers()     # Reduces to interesting subset
                self.makeSLiMs()        # Makes all SLiMs with sufficient support
                #>>> Replace here <<<#
                self.deBug(self.dict['AAFreq'])
                self.deBug(self.dict['DimFreq'])

            ### Special MotifSeq Output ###
            # This must occur after Input masking but needs no AA Frequencies or SLiMBuild #
            if self.motifSeq() and not self.opt['SlimBuild']: return     

            ###~Post-SLiMBuild Processing/Filtering before SLiMChance and Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
            ## Non-SLiMChance filtering of motifs ##
            self.dict['ElementIC'] = {}
            self.filterSLiMs()
            ## TEIRESIAS Output ##
            if self.opt['Teiresias'] or self.opt['SlimDisc']: self.teiresias()
            if not self.opt['SlimChance'] and not self.opt['SlimDisc']:
                self.log.printLog('#PROB','SlimChance=F and SlimDisc=F : No SLiM probability calculations')
                return

            ###~SLiMChance Probability and Significance Calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ptxt = 'Calculating SLiM Probabilities (p<=%s)' % rje_slim.expectString(self.stat['ProbCut'])
            sx = 0
            self.dict['PCut'] = {}
            for slim in self.dict['Slim']:
                self.log.printLog('\r#PROB','%s %.2f%% (%s Sig SLiMs)' % (ptxt,sx/self.slimNum(),rje.integerString(len(self.list['SigSlim']))),newline=False,log=False)
                sx += 100.0
                self.sigSlim(slim)
                self.wallTime()
            self.log.printLog('\r#PROB','%s complete: %s Sig SLiMs.' % (ptxt,rje.integerString(len(self.list['SigSlim']))))
            ## Reduce to Significant SLiMs ##
            self.obj['SlimList'].list['Motif'] = []
            for slim in self.dict['Slim'].keys()[0:]:
                if slim not in self.list['SigSlim']: self.dict['Slim'].pop(slim)
            for slim in self.list['SigSlim']: self.obj['SlimList'].list['Motif'].append(self.addSLiMToList(slim))

            ###~SLiMFinder Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ### Output results and tidy files. ###
            if not (self.occFilter() or self.statFilter()): self.calculateSLiMOccStats()
            self.obj['SlimList'].combMotifOccStats()
            self.tidyMotifObjects()     # Temporary solution to problem with unknown cause
            self.makeClouds()           # Identifies "clouds" of motifs - similar patterns that overlap            
            self.rankScore()            # Converts rankings into Numeric
            self.results()              # Controls SLiMFinder results output
            self.slimCheck()            # Additional SlimCheck Motifs 
            if self.opt['Extras']: self.extraOutput()   # MotifList Outputs 
            self.tarZipSaveSpace()      # Tarring, Zipping and Saving Space 

            ### End ###
            self.log.printLog('#RES','SLiMFinder results output to %s and %s*' % (self.info['ResFile'],self.info['Basefile']))
            if self.opt['Win32'] and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.stat['WallTime'] <= 0 or (time.time() - self.stat['StartTime']) < (self.stat['WallTime']*3600): raise
            if self.list['Headers']: self.results(aborted=True)
            return False # Walltime reached
        except:
            self.log.errorLog('Error in SLiMFinder.run()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <5> ### SLiMBuild Generation Methods                                                                            #
#########################################################################################################################
    def makeDimers(self):   ### Finds all possible dimers with wildcards, using MaxWild stat
        '''Finds all possible dimers with wildcards, using MinWild/MaxWild stat.'''
        try:### ~ [1] ~ Setup Dimer frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Dimers'] = {}
            self.dict['DimFreq'] = {}
            uplist = self.list['UP'][0:]
            nonx = {}   # Total count of non-X positions in UPC
            for seq in self.seqs():
                try: self.dict['DimFreq'][seq] = [(seq.nonN()-2.0)/seq.nonN()]
                except: self.dict['DimFreq'][seq] = [0.0] 
            for upc in uplist:
                try:
                    df = []
                    for seq in upc: df.append(self.dict['DimFreq'][seq][0])
                    self.dict['DimFreq'][upc] = [rje.meanse(df)[0]]
                except: self.dict['DimFreq'][upc] = [0.0] 
        except:
            self.log.errorLog('Major problem during makeDimers()')
            raise
#########################################################################################################################
    def makeSLiMs(self):    ### Makes SLiMs with enough support from Dimers
        '''Makes SLiMs with enough support from Dimers.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Slim'] = {}
            
            ### ~ [2] ~ Read in miRNA seeds to  ###
            (sx,stot) = (0.0,self.seqNum())
            for seq in self.seqs():
                ## Setup Sequence ##
                self.printLog('\r#SEED','Reading %s seeds: %.1f%%' % (rje.integerString(self.slimNum()),(sx/stot)),newline=False,log=False)
                sx += 100.0
                sequence = seq.info['Sequence'].upper()
                if self.opt['DNA']: sequence = string.replace(sequence,'N','X')
                if self.opt['Termini']: sequence = '^%s$' % sequence
                (ss,st) = (0.0,len(sequence))
                ## Setup UPC and DimFreq ##
                upc = self.getUP(seq)
                ## Find Seeds ##
                for i in range(len(sequence)):
                    ## Choose first position and check for wildcard ##
                    r = i
                    if self.opt['Termini']: r = i - 1
                    seed = sequence[i:i+8]
                    if seed.count('X') + seed.count('N') + seed.count('.') > 1: continue    # Ignore wildcard seeds
                    if seed in self.list['CommonSeed']: continue    # Skip common 8mers
                    if len(seed) < 8: continue
                    newslim = string.join(rje.strList(seed),'-0-')
                    if newslim not in self.dict['Slim']: self.dict['Slim'][newslim] = {'Occ':[],'UP':[]}
                    self.dict['Slim'][newslim]['Occ'].append((seq,r))
                    if upc not in self.dict['Slim'][newslim]['UP']: self.dict['Slim'][newslim]['UP'].append(upc)
                    if not self.dict['SeqOcc'].has_key(newslim): self.dict['SeqOcc'][newslim] = {seq:1}
                    elif not self.dict['SeqOcc'][newslim].has_key(seq): self.dict['SeqOcc'][newslim][seq] = 1
                    else: self.dict['SeqOcc'][newslim][seq] += 1
            self.log.printLog('\r#SEED','Read %s 8mer seeds from %d seq.' % (rje.integerString(self.slimNum()),self.seqNum()))
        except:
            self.log.errorLog('Major problem during makeSLiMs()')
            raise
#########################################################################################################################
    def slimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM.'''
        try:
            ###~Calculate prob of 1+ occ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            if self.opt['SeedFreq']: return self.seedFreqSlimProb(slim)
            p1 = {}         # Dictionary of {upc:chance of 1+ occ in upc}
            ##~~Setup pattern and variable-lenght multiplier~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            poslist = []    # List of AA positions in SLiM
            wildlist = []   # List of wildcard lengths in SLiM
            wild = False    # Whether next part is a wildcard length
            mult = 1        # Variable-length multiplier
            for part in string.split(slim,'-'):      # Split SLiM code in components
                ## Update lists ##
                if wild: wildlist.append(part)
                else: poslist.append(part)
                ## Calculate multiplier ##
                if wild:
                    (minx,maxx) = (self.stat['MaxWild'],0)
                    for x in part:
                        minx = min(minx,int(x))
                        maxx = max(maxx,int(x))
                    mult *= (int(maxx) - int(minx) + 1)
                wild = not wild
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.list['UP']:
                if self.dict['AADimerFreq']: (k,N,p) = self.aaDp1(slim,upc)
                else:
                    ## Setup  parameters for binomial ##
                    N = self.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
                    p = 1.0                                 # Probability of SLiM at each position
                    k = 1                                   # Number of successful trials (occurrences)
                    if self.opt['SeqOcc'] and self.slimOccNum(slim,upc) > 1: k = self.slimOccNum(slim,upc)
                    ## Calculate p and N from AAFreq and DimFreq ##
                    for pos in poslist:     # AA position
                        posfreq = 0.0
                        for aa in pos: posfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        p *= posfreq
                    for dim in wildlist:    # DimerFreq
                        dimfreq = 0.0
                        for x in dim:
                            try: dimfreq += self.dict['DimFreq'][upc][int(x)]   # Options for wildcard length
                            except: pass
                        N *= (dimfreq / len(dim))       # Mutliply by mean dimer frequency
                    N *= mult       # Each length variant is effectively another position the SLiM could occur
                    if p > 1: p = 1.0   # Cannot in reality have p > 1!
                    ## Calculate binomial ##
                    p1[upc] = rje.binomial(k,N,p,usepoisson=False,callobj=self)
            ## Extra verbosity. Remove at some point? ##
            self.deBug('%s: %s' % (patternFromCode(slim),p1.values()))
            self.deBug('%s: %s vs %s' % (patternFromCode(slim),self.slimUP(slim),sum(p1.values())))

            ###~Calculate overall probability of observed support~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## All observed occurrences ##
            self.dict['Slim'][slim]['ExpUP'] = sum(p1.values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimUP(slim), self.UPNum(), self.dict['Slim'][slim]['ExpUP']/self.UPNum()) # Use mean p1+
            if k <= 0: self.dict['Slim'][slim]['Prob'] = 1.0
            else: self.dict['Slim'][slim]['Prob'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            ## Catch for binomial problems. Should no longer happen. ##
            if self.dict['Slim'][slim]['Prob'] <= 0:    # Shouldn't happen now! #
                self.log.errorLog('Probability for %s <= 0.0 due to numerical limitations: Given arbitrary 1e-16.!' % (patternFromCode(slim)),printerror=False)
                self.dict['Slim'][slim]['Prob'] = 1e-16
                self.deBug(p1)
            ## Correction for restricted focal sequences ##
            self.focusAdjustment(slim)
        except:
            self.log.errorLog('Error with slimProb(%s)' % slim)
            self.dict['Slim'][slim]['Prob'] = 1.0
#########################################################################################################################
    def seedFreqSlimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM.'''
        try:
            ###~Calculate prob of 1+ occ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            p1 = {}         # Dictionary of {upc:chance of 1+ occ in upc}
            seed = patternFromCode(slim)
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.list['UP']: p1[upc] = min(self.dict['MST'][upc] * self.dict['SeedFreq'][seed],1.0)
            ## Extra verbosity. Remove at some point? ##
            self.deBug('%s: %s' % (patternFromCode(slim),p1.values()))
            self.deBug('%s: %s vs %s' % (patternFromCode(slim),self.slimUP(slim),sum(p1.values())))

            ###~Calculate overall probability of observed support~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## All observed occurrences ##
            self.dict['Slim'][slim]['ExpUP'] = sum(p1.values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimUP(slim), self.UPNum(), self.dict['Slim'][slim]['ExpUP']/self.UPNum()) # Use mean p1+
            if k <= 0: self.dict['Slim'][slim]['Prob'] = 1.0
            else: self.dict['Slim'][slim]['Prob'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            ## Catch for binomial problems. Should no longer happen. ##
            if self.dict['Slim'][slim]['Prob'] <= 0:    # Shouldn't happen now! #
                self.log.errorLog('Probability for %s <= 0.0 due to numerical limitations: Given arbitrary 1e-16.!' % (patternFromCode(slim)),printerror=False)
                self.dict['Slim'][slim]['Prob'] = 1e-16
                self.deBug(p1)
            ## Correction for restricted focal sequences ##
            self.focusAdjustment(slim)
        except:
            self.log.errorLog('Error with slimProb(%s)' % slim)
            self.dict['Slim'][slim]['Prob'] = 1.0
#########################################################################################################################
def patternFromCode(slim): return slimfinder.patternFromCode(slim)
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: miRNA(mainlog,cmd_list).run()
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
