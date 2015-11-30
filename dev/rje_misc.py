#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_misc
Description:  Miscellenous script storage module
Version:      0.50.0
Last Edit:    12/04/15
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for "one-off" scripts for odd-jobs that often need to be run multiple times and then forgotten
    forever:

    - aicpaper = Conversion of AIC PATIS tables for Jan 2014 AIC paper.
    - crisp = Analysis of Crisp et al paper data
    - testing = SWC boot camp testing
    - holger = reformat fasta
    - yunan = combine and normalised Yunan data
    - sageshape = reshape the SuperSAGE data for ANOVA and FDR analysis
    - 3did = generate 3DID benchmarking datasets
    - sfstick = Stick together SF files
    - dtfas = Reformat DT_Sequence data    
    - laavanya = Prepares PTPRD PPI datasets for SLiMFinder
    - elm_presto_qsub = Makes a list of python command for qsub job
    - elm_presto_compile = Compiles all the PRESTO results from ELM conservation analyses into a single file
    - ppi_go_datafarm = Copies PPI and GO datasets into the current MotifAln folder
    - minimotif = Reformats mini_motif.txt into standard motif RegExps
    - wag = reconstruct a PAM1 matrix from the Goldman WAG matrix
    - hprd_elm_copy = Copies HPRD ELM datasets for SLiMDisc analysis
    - hprd_elm_comp = Extracts rows of *.compare.tdt files where the HPRD dataset matches its ELM
    - neduvamotif = Reformats raw motif file from Neduva & Russell work
    - disordertest = Tests the two disorder prediction servers and compares
    - ygob = Makes orthologue alignments from YGOB files
    - ygob2 = Makes orthologue alignments from YGOB files following the original ygob!
    - ensloci = Counts EnsEMBL loci
    - go_cc = Extracts a list of GO cellular component IDs into a file corresponding to SLiMDisc input
    - elm_gcut = Cuts down the ELM GABLAM results as they *should* have been output!
    - ensgo = Makes GO Datasets from EnsEMBL Loci data and BioMart download (& ensgo2)
    - elmdup = Identifies "duplicate" proteins from different species using orthaln.fas files
    - elmaln = Sorts ELM alignments
    - elmFT = Makes a Table of SwissProt Annotations from mapping.tdt and UniProt
    - termX = Reformats all fasta files in the directory to be UC for X aa at each end of the sequence and LC for rest
    - 143ppi = Make a table of 14-3-3 interacting proteins from HPRD
    - allslim = runs slimfinder on all datasets (*.dat & *.fas) in directory
    - slimtest = Makes a bunch of random datasets to run SLiMFinder on
    - slimtestsum = Makes a summary table of random datasets run with SLiMFinder
    - sfvst = SLiMBuild vs. TEIRESIAS time test w/o ambiguity
    - ensdat = EnsDat cleanup following rje_ensembl boo-boo
    - slimcons = Custom assembly of EnsDatHuman SLiMFinder runs and CompariMotif results for data analysis in R
    - locustdb = Temp analysis of Locust EST database
    - phomo = PhosphoMotif Reformatter
    - flyseq = Extract sequences from Flybase chado_xml
    - arath2go = Convert EMBL_ARATH sequence to GO using Arabidopsis download
    - intenrich = Reformat enriched integrin proteomics results for Pingu
    - sftargz = TarGZ cleanup of SLiMFinder results
    - taxadbsum = E hux Taxa DB summary cleanup
    - bencog = Makes Ben COG links
    - humsf09 = Human SF09 cleanup
    - biol3050 = BIOL3050 clone vs ENST cDNA matchup
    - jrjspf = Reformatting an SPF file for Joe Jenkins

Commandline:
    job=X       : Identifier for the job to be performed [None]
    infile=FILE : Name of input file for relevant task [None]

Uses general modules: glob, os, string, sys, time
Uses RJE modules: rje, rje_blast, rje_db, rje_disorder, rje_seq, rje_uniprot, rje_zen
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_db, rje_seq, rje_seqlist, rje_sequence, rje_uniprot, rje_xml, rje_zen
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.1 - Initial Compilation with elm_presto_compile job.
    # 0.2 - Added ppi_go_datafarm job
    # 0.3 - minimotif = Reformats mini_motif.txt into standard motif RegExps
    # 0.4 - Added elm_presto_qsub
    # 0.5 - wag = reconstruct a PAM1 matrix from the Goldman WAG matrix
    # 0.6 - hprd_elm_copy = Copies HPRD ELM datasets for SLiMDisc analysis
    # 0.7 - hprd_elm_comp = Extracts rows of *.compare.tdt files where the HPRD dataset matches its ELM
    # 0.8 - neduvamotif = Reformats raw motif file from Neduva & Russell work
    # 0.9 - disordertest = Tests the two disorder prediction servers and compares
    # 0.10 - ygob = Makes orthologue alignments from YGOB files
    # 0.11 - ensloci = Counts EnsEMBL loci
    # 0.12 - go_cc = Extracts a list of GO cellular component IDs into a file corresponding to SLiMDisc input
    # 0.13 - elm_gcut = Cuts down the ELM GABLAM results as they *should* have been output!
    # 0.14 - ensgo = Makes GO Datasets from EnsEMBL Loci data and BioMart download
    # 0.15 - elmdup = Identifies "duplicate" proteins from different species using orthaln.fas files
    # 0.16 - elmaln = Sorts ELM alignments
    # 0.17 - elmFT = Makes a Table of SwissProt Annotations from mapping.tdt and UniProt
    # 0.18 - termX = Reformats all fasta files in the directory to be UC for X aa at each end of the sequence and LC for rest
    # 0.19 - allslim = runs slimfinder on all datasets (*.dat & *.fas) in directory
    # 0.20 - slimtest = Makes a bunch of random datasets to run SLiMFinder on
    # 0.21 - slimtestsum = Makes a summary table of random datasets run with SLiMFinder
    # 0.22 - sfvst = SLiMBuild vs. TEIRESIAS time test w/o ambiguity
    # 0.23 - mailer = E-Mail test
    # 0.24 - ensdat = EnsDat cleanup following rje_ensembl boo-boo
    # 0.25 - slimcons = Custom assembly of EnsDatHuman SLiMFinder runs and CompariMotif results for data analysis in R
    # 0.26 - locustdb = Temp analysis of Locust EST database
    # 0.27 - phomo = PhosphoMotif Finder reformatting
    # 0.28 - flyseq = Extract sequences from Flybase chado_xml
    # 0.29 - arath2go = Convert EMBL_ARATH sequence to GO using Arabidopsis download
    # 0.30 - intenrich = Reformat enriched integrin proteomics results for Pingu
    # 0.31 - laavanya = Prepares PTPRD PPI datasets for SLiMFinder
    # 0.32 - sftargz = TarGZ cleanup of SLiMFinder results
    # 0.33 - taxadbsum = E hux Taxa DB summary cleanup
    # 0.34 - bencog = Makes Ben COG links
    # 0.35 - sftidy = SLiMFinder grand results tidy up.
    # 0.36 - prpsi = Make SI table for PRP project
    # 0.37 - humsf09 = Human SF09 cleanup
    # 0.38 - biol3050 = BIOL3050 clone vs ENST cDNA matchup
    # 0.39 - biol2018 = BIOL2018 domain identification.
    # 0.40 - dtfas = Reformat DT_Sequence data
    # 0.41 - sfstick = Stick together SF files
    # 0.42 - 3did = generate 3DID benchmarking datasets
    # 0.43 - sageshape = reshape the SuperSAGE data for ANOVA and FDR analysis
    # 0.44 - yunan = combine and normalised Yunan data
    # 0.45 - kroc = Kieren ROC analysis
    # 0.46 - holger = reformat fasta
    # 0.47 - testing = SWC boot camp misc testing code
    # 0.48 - aicpaper = Conversion of AIC PATIS tables for Jan 2014 AIC paper.
    # 0.49.0 - jrjspf = Reformatting an SPF file for Joe Jenkins
    # 0.50.0 - crisp = Analysis of Crisp et al paper data
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : General Tidy up if bothered. (Not really that kind of module!)
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_MISC', '0.50.0', 'April 2015', '2007')
    description = 'Miscellaneous Odd-jobs Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Wise %s always says:' % rje_zen.Zen()._noun(),'\t"%s"' % rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        print cmd_list, help
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?','N'): out.verbose(-1,4,text=rje.__doc__)
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
### SECTION II: OddJob CLASS                                                                                            #
#########################################################################################################################
class OddJob(rje.RJE_Object):     
    '''
    Misc OddJobs Class. Author: Rich Edwards (2005).

    Info:str
    - Name = name of job (job=X command)
    - InFile = name of input file for relevant task [None]
    
    Opt:boolean

    Stat:numeric

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics 
        self.infolist = ['InFile']
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ### 
                self._cmdRead(cmd,type='info',att='Name',arg='job')
                self._cmdReadList(cmd,'file',['InFile'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### OddJob Run Methods                                                                                      #
#########################################################################################################################
    def run(self):      ### Runs odd-job specified by self.info['Name']
        '''Runs odd-job specified by self.info['Name'].'''
        try:### Laavanya ###
            if self.info['Name'] == 'laavanya': self.laavanya()
            elif self.info['Name'] == 'crisp': self.Crisp()
            elif self.info['Name'] == 'jrjspf': self.JRJSPF()
            elif self.info['Name'] == 'aicpaper': self.aicPaper() # Conversion of AIC PATIS tables for Jan 2014 AIC paper.
            elif self.info['Name'] == 'testing': self.testing()
            elif self.info['Name'] == 'holger': self.holger()
            elif self.info['Name'] == 'kroc': self.kroc()
            elif self.info['Name'] == 'yunan': self.yunan()
            elif self.info['Name'] in ['shapesage','sageshape']: self.shapeSAGE()
            elif self.info['Name'] == '3did': self.dmi()
            elif self.info['Name'] == 'patfas': self.patFas()
            elif self.info['Name'] == 'sfstick': self.sfStick()
            elif self.info['Name'] == 'sfstick2': self.sfStick2()
            elif self.info['Name'] == 'dtfas': self.dtFas()
            elif self.info['Name'] == 'biol3050': self.biol3050()
            elif self.info['Name'] == 'biol2018': self.biol2018()
            ### PRP SI V.36 ###
            elif self.info['Name'] == 'prpsi': self.prpSI()
            ### SLiMFinder grand results tidy up. V35
            elif self.info['Name'] == 'sftidy': self.sfTidy()
            elif self.info['Name'] == 'humsf09': self.SF09()
            ### Ben COG analysis ###
            elif self.info['Name'] == 'bencog': self.benCOG()
            ### E hux TaxaDB ###
            elif self.info['Name'] == 'taxadbsum': self.taxaDBSum()
            ### IntEnrich ###
            elif self.info['Name'] == 'intenrich': self.intEnrich()
            ### SLiMFinder TarGZ ###
            elif self.info['Name'] == 'sftargz': self.sfTarGZ()
            ### FlyBase ###
            elif self.info['Name'] == 'flyseq': self.flySeq()     # Extract sequences from Flybase chado_xml
            ### SlimCons / RedRed Play ###
            elif self.info['Name'] == 'slimcons': self.slimCons()
            ### arath2go = Convert EMBL_ARATH sequence to GO using Arabidopsis download ###
            elif self.info['Name'] == 'arath2go': self.arath2GO()
            ### LocustDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elif self.info['Name'] == 'locustdb':
                seqs = rje_seq.SeqList(self.log,self.cmd_list)
                newseq = {} # rje_seq.SeqList(self.log,self.cmd_list+['seqin=','autoload=F'])
                bestrf = {}
                sx = 0
                for seq in seqs.seq:
                    sx += 1
                    seq.sixFrameTranslation()
                    text = '%d|%s: ' % (sx,seq.shortName())
                    best = (0,0)
                    ox = 0
                    for rf in rje.sortKeys(seq.dict['Translation']):
                        for orf in string.split(seq.dict['Translation'][rf],'*'):
                            if len(orf) < seqs.stat['MinLen']: continue
                            ox += 1
                            newname = '%s_ORF%s RF%d; %d aa.' % (seq.shortName(),rje.preZero(ox,99),rf,len(orf))
                            newseq[newname] = orf
                            #x#newseq._addSeq(newname,orf)
                        aa = len(rje_sequence.bestORF(seq.dict['Translation'][rf]))
                        if aa > best[0]: best = (aa,rf)
                        text += '%d=%d; ' % (rf,aa)
                    self.log.printLog('#ORF','%s %d ORFs' % (text,ox))
                    if best[1] not in bestrf: bestrf[best[1]] = 0
                    bestrf[best[1]] += 1
                self.log.printLog('#BEST','%s' % bestrf)
                #newseq.saveFasta(seqfile='locust_orfs.fas')
                FAS = open('locust_orfs.fas','w')
                (sx,stot) = (0,len(newseq))
                for name in rje.sortKeys(newseq):
                    FAS.write('>%s\n%s\n' % (name,newseq[name]))
                    sx += 1
                    self.log.printLog('\r#FAS','Saving translated ORFs: %.1f%% (%s)' % (sx*100.0/stot,rje.integerString(sx)),newline=False,log=False)
                FAS.close()
                self.log.printLog('\r#FAS','Saving translated ORFs complete: %s sequences.' % (rje.integerString(sx)))

            ### PhosphoMotif ###
            elif self.info['Name'] == 'phomo': self.phoMo()     ### PhosphoMotif Finder reformatting
            ### EnsDat cleanup
            elif self.info['Name'] == 'ensdat': self.ensDat() ### EnsDat cleanup following rje_ensembl boo-boo
            ### SLiMFinder ###
            elif self.info['Name'] == 'sfvst':    ### SLiMBuild vs. TEIRESIAS time test w/o ambiguity
                self.slimfinderVsTeiresias()
            elif self.info['Name'] == 'slimtest': ### Makes a bunch of random datasets to run SLiMFinder on
                self.slimTest()
            elif self.info['Name'] == 'slimtestsum': ### Makes a bunch of random datasets to run SLiMFinder on
                self.slimTestSum()
            elif self.info['Name'] == 'allslim': ### runs slimfinder on all datasets (*.dat & *.fas) in directory
                ifiles = glob.glob('*.dat') + glob.glob('*.fas')
                import slimfinder
                for i in ifiles:
                    if i == ifiles[0]:
                        slimfinder.SLiMFinder(self.log,['resfile=allslim.csv']+self.cmd_list+['seqin=%s' % i]).run()
                    else:
                        slimfinder.SLiMFinder(self.log,['resfile=allslim.csv']+self.cmd_list+['seqin=%s' % i,'append=T']).run()
            ### Make a table of 14-3-3 interacting proteins from HPRD ###
            elif self.info['Name'] == '143ppi':
                self.ppi1433()
            ### Terminal Case ###
            elif rje.matchExp('^term(\d+)',self.info['Name']):
                t = string.atoi(rje.matchExp('^term(\d+)',self.info['Name'])[0])
                files = glob.glob('*.fas')
                for fas in files:
                    seqlist = rje_seq.SeqList(self.log,['seqin=%s' % fas])
                    uniprot = rje_uniprot.UniProt(self.log,['ucft=IM','lcft=EM'])
                    for seq in seqlist.seq:
                        if seq.aaLen() > (2*t):
                            seq.info['Sequence'] = seq.info['Sequence'][:t] + seq.info['Sequence'][t:-t].lower() + seq.info['Sequence'][-t:]
                        seq.addSequence(seq.info['Sequence'])
                        uniprot.addFromSeq(seq)
                    seqlist.saveFasta(case=True)
                    uniprot.saveUniProt(string.replace(fas,'.fas','.dat'),append=False)
            ### ELM Protein Stuff ###
            elif self.info['Name'] == 'elmdup':
                self.elmDup()
            elif self.info['Name'] == 'elmaln':
                self.elmAln()
            elif self.info['Name'] == 'elmdom':
                self.elmDom()
            ### Count EnsEMBL loci ###
            elif self.info['Name'] == 'ensloci':
                self.ensLoci()
            ### YGOB Orthologue alignments ###
            elif self.info['Name'] == 'ygob': # Makes orthologue alignments from YGOB files
                self.ygobOrthologues()
            elif self.info['Name'] == 'ygob2': # Makes orthologue alignments from YGOB files
                self.ygobOrthologues2()
            ### Disorder Test ###
            elif self.info['Name'] == 'disordertest': # Tests the two disorder prediction servers and compares
                self.disorderTest() 
            ### Neduva & Russell significant motifs ###
            elif self.info['Name'] == 'neduvamotif':
                raw = 'C:\\Documents and Settings\\redwards\\My Documents\\redwards\\My Projects\\Motifs\\Literature Motifs\\N&R_Raw.txt'
                rlines = self.loadFromFile(raw,chomplines=True)
                ## Read motifs ##
                motifs = [] ### List of motif entries
                ppi = ''    ### PPI Dataset
                data = ''   ### Dataset returning motif
                for line in rlines:
                    if rje.matchExp('\*\*\*(\S+)\*\*\*',line):  ### PPI dataset
                        ppi = rje.matchExp('\*\*\*(\S+)\*\*\*',line)[0]
                    elif rje.matchExp('^(\S+)\s+(\d\S+)\s+(\d+)\s+(\d+)\s+(\d\S+)',line):   ### Motif
                        dx += 1
                        (motif,scons,occ,seqnum,pvalue) = rje.matchExp('^(\S+)\s+(\d\S+)\s+(\d+)\s+(\d+)\s+(\d\S+)',line)
                        motifs.append('%s_%s_%d %s\t# (%s/%s); Scons=%s; p=%s;' % (ppi,data,dx,motif,occ,seqnum,scons,pvalue))
                        self.log.printLog('\r#MOT','Processing: %d motifs.' % len(motifs),newline=False,log=False)
                    elif line:   ### Prot/Domain dataset
                        info = string.split(line)
                        if len(info) > 1:
                            data = '%s_%s' % (info[0],info[1])
                        else:
                            data = info[0]
                        dx = 0
                ## Output ##
                motifs.sort()
                outfile = string.replace(raw,'Raw','SigMotif')
                OUT = open(outfile,'w')
                OUT.write(string.join(motifs,'\n'))
                OUT.close()
                self.log.printLog('\r#MOT','Processing: %d motifs output.' % len(motifs))                
            
            ### hprd_elm_copy = Copies HPRD ELM datasets for SLiMDisc analysis ###
            elif self.info['Name'] == 'hprd_elm_copy':
                hfile = glob.glob('Datasets/*.fas')
                for h in hfile:
                    d = os.path.splitext(os.path.basename(h))[0]
                    for v in ['amb','uhs','mts']:
                        print 'cp %s SlimDisc/%s_%s.fas' % (h,d,v)
                        os.system('cp %s SlimDisc/%s_%s.fas' % (h,d,v))

            ### hprd_elm_comp = Extracts rows of *.compare.tdt files where the HPRD dataset matches its ELM ###
            elif self.info['Name'] == 'hprd_elm_comp':
                cfile = glob.glob('*.compare.tdt')
                OUT = open('hprd_elm.matches.tdt','w')
                rje.writeDelimit(OUT,['file','HPRD','ELM','Motif1','Motif2','Sim1','Sim2','MatchPos','Code'],delimit='\t')
                mx = 0
                results = {}
                for comp in cfile:
                    file = string.replace(comp,'.compare.tdt','')
                    results[file] = {}
                    for cline in self.loadFromFile(comp)[1:]:
                        data = string.split(rje.chomp(cline),'\t')
                        if len(data) < 8:
                            continue
                        try:
                            hprd = rje.matchExp('HPRD0*([1-9]\d+)_',data[0])[0]
                            if rje.matchExp('(\D%s\D)' % hprd,data[1]):
                                try:
                                    (dataset,best) = rje.matchExp('^(\S+)_(\d+)$',data[0])
                                except:
                                    self.errorLog(rje.matchExp('^(\S+)_(\d+)$',data[0]))
                                if not results[file].has_key(dataset):
                                    results[file][dataset] = {'Best':best,'Total':0}
                                results[file][dataset]['Total'] += 1
                                rje.writeDelimit(OUT,[file] + data,'\t')
                                mx += 1
                                print '\r', mx,
                            else:
                                self.deBug('\D%s\D not in %s' % (hprd,data[1]))
                        except:
                            self.log.errorLog('Balls!')
                            print '\n', file, data
                print
                OUT.close()
                OUT = open('hprd_elm.performance.tdt','w')
                rje.writeDelimit(OUT,['file','dataset','best','total'],delimit='\t')
                for file in results.keys():
                    for dataset in results[file].keys():
                        rje.writeDelimit(OUT,[file,dataset,results[file][dataset]['Best'],results[file][dataset]['Total']],delimit='\t')
                OUT.close()
                
            ### ELM PRESTO QSub ###
            elif self.info['Name'] == 'elm_presto_qsub':
                Q = open('elm_presto_qsub.tmp','w')
                methods = ['abs','pam','prop']
                weights = [-3,3]    # [2,-1,-2]
                ambs = 'TF'
                revs = 'TF'
                extra = True
                for m in methods:
                    call = 'python /home/richard/Python_Modules/presto.py ini=newcons.ini conscore=%s resfile=%s' % (m,m)
                    for w in weights:
                        for a in ambs:
                            for r in revs:
                                if extra:
                                    Q.write('%s_w%s_a%s_r%s consweight=%s consamb=%s reverse=%s winsa=3 winhyd=3\n' % (call,w,a,r,w,a,r))
                                    extra = False
                                else:
                                    Q.write('%s_w%s_a%s_r%s consweight=%s consamb=%s reverse=%s\n' % (call,w,a,r,w,a,r))
                Q.close()
            
            ### ELM PRESTO COMPILE ###
            elif self.info['Name'] == 'elm_presto_compile':
                csvlist = glob.glob('*.presto.csv')
                COMP = open('newcons_presto.csv','w')
                head = True
                for csv in csvlist:
                    conscore = csv[:csv.find('.presto')]
                    SUB = open(csv,'r')
                    lines = SUB.readlines()
                    if head:
                        COMP.write('conscore,%s' % lines[0])
                        head = False
                    for line in lines[1:]:
                        COMP.write('%s,%s' % (conscore,line))
                    SUB.close()
                COMP.close()

            ### SlimPicks PPI/GO Datasets ###
            elif self.info['Name'] == 'ppi_go_datafarm':
                mlist = glob.glob('*.slimpicks.fas')
                for m in mlist:
                    d = m[:-len('.slimpicks.fas')]
                    if d[:2] == 'GO':
                        dpath = '/home/richard/SlimDisc_Projects/IPI_GO/'
                    else:
                        dpath = '/home/richard/SlimDisc_Projects/PPI_Datasets/'
                    print 'cp %s%s.fas .' % (dpath,d)
                    os.system('cp %s%s.fas .' % (dpath,d))
                print 'Done!'

            ### New ELM Reformat (Jan '08) ###
            elif self.info['Name'] == 'newelm':
                mlines = self.loadFromFile(filename='D:\\Databases\\PhosphoMotif Finder\\elms.dat')
                OUT = open('D:\\Databases\\PhosphoMotif Finder\\ELM.motifs','w')
                for m in mlines:
                    details = rje.matchExp('^(\S+)\s+(\S+)\s+(\S.+)\s+$',m)
                    OUT.write('%s\t%s\t# %s [ELM]\n' % details)
                OUT.close()
            ### MiniMotif Reformat ###
            elif self.info['Name'] == 'minimotif':
                mlines = self.loadFromFile(filename='D:\\Databases\\PhosphoMotif Finder\\mini_motif_raw.txt')
                OUT = open('D:\\Databases\\PhosphoMotif Finder\\MnM.motifs','w')
                for m in mlines:
                    details = rje.matchExp('^(\S+)\s+(\S+)\s+(\S.+)',m)
                    if details:
                        (name,pattern,desc) = details
                        ## ~ Modify Pattern ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        old = pattern
                        if pattern[0] == '<': pattern = string.replace(pattern,'<','^')
                        else: pattern = '-%s' % pattern
                        pattern = string.replace(pattern,'>',' $')
                        pattern = string.replace(pattern,'?','.')
                        pattern = string.replace(pattern,'-',' -')
                        while rje.matchExp('(\-(\S+) )',pattern):
                            match = rje.matchExp('(\-(\S+) )',pattern)
                            self.deBug(match)
                            if string.count(match[1],'/') > 0:
                                amb = string.replace(match[1],'/','')
                                pattern = string.replace(pattern,match[0],'[%s]' % amb)
                            else: pattern = string.replace(pattern,match[0],match[1])
                            #self.log.printLog('\r#MOT','%s => %s         |<=' % (old,pattern),newline=False,log=False)
                        while rje.matchExp('(\-(\S+))$',pattern):
                            match = rje.matchExp('(\-(\S+))$',pattern)
                            if string.count(match[1],'/') > 0:
                                amb = string.replace(match[1],'/','')
                                pattern = string.replace(pattern,match[0],'[%s]' % amb)
                            else: pattern = string.replace(pattern,match[0],match[1])
                        pattern = string.replace(pattern,' ','')
                        self.log.printLog('\r#MOT','%s => %s' % (old,pattern))
                        ## ~ Modify Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                        desc = string.join(string.split(string.replace(desc,'"','')))
                        ## >> Stupid MnM << ##
                        desc = string.replace(desc,'&dopt=Abstract; http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=',',')
                        desc = string.replace(desc,'&dopt=Abstract; http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=retrieve&db=pubmed&list_uids=',',')
                        desc = string.replace(desc,'; http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Retrieve&dopt=AbstractPlus&list_uids=',',')
                        desc = string.replace(desc,', http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=',',')
                        desc = string.replace(desc,';http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?itool=abstractplus&db=pubmed&cmd=Retrieve&dopt=abstractplus&list_uids=',',')
                        desc = string.replace(desc,';http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=',',')
                        desc = string.replace(desc,' ; http://','&')
                        desc = string.replace(desc,'; http://','&')
                        ## >> Stupid MnM End << ##
                        try: (words,subcell) = rje.matchExp('^(\S.+)\s+(\S+)\s+http',desc)
                        except:
                            print '>>>\n%s\n<<<' % m
                            continue
                        try: pmid = string.split(rje.matchExp("list_uids=(\S+)",desc)[0],'&')[0]
                        except: pmid = ''
                        source = string.split(desc)
                        s = -1
                        while source[s-1][:4] != 'http': s -= 1
                        source = string.join(source[s:])
                        desc = '# %s [CC=%s] [MnM]' % (words,subcell)
                        if pmid: desc += ' [PMID:%s]' % pmid
                        if source != 'N/A': desc += '[%s]' % source
                        OUT.write('%s\t%s\t%s\n' % (name,pattern,desc))
                        print '%s\t%s\t%s' % (name,pattern,desc)
                OUT.close()

            ### WAG ###
            elif self.info['Name'] == 'wag':
                self.wag()

            ### IPI GO CC ###
            elif self.info['Name'] == 'go_cc':  # Extracts a list of GO cellular component IDs into a file corresponding to SLiMDisc input
                gotxt = '/home/richard/Motifs/IPI_GO_Source_Data/GO.txt'
                DO = open('go_map.do','w')
                cc = []
                for line in self.loadFromFile(gotxt):
                    match = rje.matchExp('^GO:(\d+)\s+(\S.+\S)\s+C\s*$',line)
                    if match:
                        (go,desc) = match
                        DO.write('replace desc = "%s" if dataset == "GO_%s"\n' % (desc,go))
                        cc.append('GO_%s' % go)
                        self.log.printLog('\r#GO','%s CC GO terms.' % rje.integerString(len(cc)),newline=False,log=False)
                self.log.printLog('\r#GO','%s CC GO terms.' % rje.integerString(len(cc)))
                DO.close()
                open('go_cc.txt','w').write(string.join(cc,'\n'))

            ### IPI GO CC II ###
            elif self.info['Name'] == 'go_cc2':     # Extracts GO cellular component Datasets
                cc = []
                cx = 0
                for go in self.loadFromFile('go_cc.txt',chomplines=True):
                    cx += 1
                    if os.path.exists('/home/richard/SlimDisc_Projects/IPI_GO/%s.fas' % go):
                        os.mkdir('GO_CC_SLiMDisc/%s/' % go)
                        os.system('cp /home/richard/SlimDisc_Projects/IPI_GO/%s.fas GO_CC_SLiMDisc/.' % go)
                        os.system('cp /home/richard/SlimDisc_Projects/IPI_GO/%s/*.out /home/richard/SlimDisc_Projects/IPI_GO/%s/*.fasta GO_CC_SLiMDisc/%s/.' % (go,go,go))
                        cc.append(go)
                        self.log.printLog('\r#GO','%s of %s CC GO datasets.' % (rje.integerString(len(cc)),rje.integerString(cx)),newline=False,log=False)
                self.log.printLog('\r#GO','%s of %s CC GO datasets.' % (rje.integerString(len(cc)),rje.integerString(cx)))
                open('go_fas.txt','w').write(string.join(cc,'\n'))

            ### IPI GO CC III ###
            elif self.info['Name'] == 'go_cc3':     # Cleanup Files
                for go in glob.glob('*.fas.bak'):
                    gcmd = self.cmd_list + ['seqin=%s' % go,'seqout=%s' % go[:-4],'memsaver=F','accnr=T','reformat=fasta','autofilter=True']
                    rje_seq.SeqList(self.log,gcmd).saveFasta(seqfile=go[:-4])

            ### EnsLoci GO Datasets ###
            elif self.info['Name'] == 'ensgo':  # Makes GO Datasets from EnsEMBL Loci data and BioMart download
                self.ensLociGo()

            ### ELM GABLAM CLEANUP ###
            elif self.info['Name'] == 'elm_gcut':   # Cuts down the ELM GABLAM results as they *should* have been output!
                gabs = glob.glob('*.gablam.tdt')[0:]
                for gab in gabs:
                    os.rename(gab,'%s.old' % gab)
                    GABIN = open('%s.old' % gab,'r')
                    GABOUT = open(gab,'w')
                    for line in GABIN.readlines():
                        gdata = string.split(line)
                        if not gdata:
                            continue
                        if gdata[0] == 'Qry':   # Header
                            GABOUT.write(line)
                        elif len(gdata) < 19:
                            continue
                        else:
                            q = string.atof(gdata[11])
                            h = string.atof(gdata[17])
                            if q > 90 or h > 90:
                                GABOUT.write(line)
                    GABIN.close()
                    GABOUT.close()

            ### E-Mail ###
            elif self.info['Name'] == 'mailer':
                ### ~ Setup message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                send_from = 'mailer@cabbagesofdoom.co.uk'
                send_to = 'richard.edwards@ucd.ie'
                subject = 'Test'
                text = 'This is a test random mail'
                server = "mail.ucd.ie"

                ### ~ Setup EMail ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
                import smtplib, email
                from email.MIMEMultipart import MIMEMultipart
                from email.MIMEBase import MIMEBase
                from email.MIMEText import MIMEText
                from email.Utils import COMMASPACE, formatdate
                from email import Encoders
                msg = MIMEMultipart()
                msg['From'] = send_from
                msg['To'] = COMMASPACE.join(send_to)
                msg['Date'] = formatdate(localtime=True)
                msg['Subject'] = subject

                msg.attach(MIMEText(text))

                ### ~ Send EMail ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                #x#smtp = smtplib.SMTP(server)
                #x#smtp.sendmail(send_from, send_to, msg.as_string())
                #x#smtp.close()
                '''
                    os.system('/usr/lib/sendmail -oi -t -odq")
                    or die "Can't fork for sendmail: $!\n";
                    print SENDMAIL <<"EOF";
                    From: User Originating Mail <me\@host>
                    To: Final Destination <you\@otherhost>
                    Subject: A relevant subject line

                    Body of the message goes here, in as many lines as you like.
                    EOF
                    close(SENDMAIL)     or warn "sendmail didn't close nicely";
                '''

            ### Nothing! ###
            else:
                self.log.errorLog('Job "%s" not recognised.' % self.info['Name'],printerror=False)
        except:
            self.log.errorLog('Error in rje_misc.run(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def Crisp(self):   # SPF file conversion
        '''SPF file conversion.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infile = self.getStr('InFile')
            db = rje_db.Database(self.log,self.cmd_list)
            db.basefile(rje.baseFile(infile,strip_path=True))
            hdb = db.addTable(infile,mainkeys=['gene_ID'],name='hgt.nr')
            ### ~ [1] Compress data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdb.keepFields(string.split('species,gene_ID,desc,transcript_ID,h_index,blastx_taxon,og,metazoa_uniprot_best_hit_accession,archaea_uniprot_best_hit_accession,bacteria_uniprot_best_hit_accession,fungi_uniprot_best_hit_accession,protist_uniprot_best_hit_accession,plants_uniprot_best_hit_accession',','))
            taxa = []; rules = {'n_taxon':'sum','nr_taxon':'sum'}
            for field in hdb.fields()[0:]:
                if field.endswith('_accession'):
                    taxon = string.split(field,'_')[0]
                    if taxon == 'plants': taxon = 'plant'
                    taxa.append(taxon)
                    hdb.renameField(field,taxon)
            newfields = ['species','blastx_taxon','n_taxon','nr_taxon']
            taxa.sort()
            taxa.insert(0,'og')
            for taxon in taxa:
                hdb.addField('n_%s' % taxon,evalue=1); rules['n_%s' % taxon] = 'sum'; rules['nr_%s' % taxon] = 'sum'
                newfields += ['n_%s' % taxon,'nr_%s' % taxon]
                for entry in hdb.entries():
                    if entry[taxon] == 'NA': entry[taxon] = ''; entry['n_%s' % taxon] = 0
            hdb.compress(['species','blastx_taxon'],rules,default='list',joinchar='|')
            for taxon in taxa: hdb.addField('nr_%s' % taxon)
            for entry in hdb.entries():
                for taxon in taxa:
                    entry['nr_%s' % taxon] = len(string.split(entry[taxon],'|'))
                try: entry['n_taxon'] = entry['n_%s' % entry['blastx_taxon']]
                except: entry['n_taxon'] = 0
                try: entry['nr_taxon'] = entry['nr_%s' % entry['blastx_taxon']]
                except: entry['nr_taxon'] = 0
            ### ~ [2] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for field in hdb.fields():
                if field not in newfields: newfields.append(field)
            hdb.list['Fields'] = newfields
            hdb.saveToFile()
            ### ~ [3] Totals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hdb.compress(['species'],rules,default='list',joinchar='|')
            for entry in hdb.entries(): entry['blastx_taxon'] = 'total'
            hdb.info['Name'] = 'hgt.totals'
            hdb.saveToFile()
        except: self.errorLog('Error in rje_misc.JRJSPF()')
#########################################################################################################################
    def JRJSPF(self):   # SPF file conversion
        '''SPF file conversion.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infile = self.getStr('InFile')
            db = rje_db.Database(self.log,self.cmd_list)
            db.basefile(rje.baseFile(infile))
            sdb = db.addTable(infile,mainkeys='#',delimit='\t',name='SPF.Mod')
            levels = {'Level_1':'k','Level_2':'p','Level_3':'c','Level_4':'o','Level_5':'f','Level_6':'g','Level_7':'s'}
            # k__Bacteria	p__Proteobacteria	c__Alphaproteobacteria	o__Rhodospirillales	f__Rhodospirillaceae	g__	s__	denovo44
            # Unassigned	unclassified	unclassified	unclassified	unclassified	unclassified	unclassified	denovo49
            ### ~ [1] Modify Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dupnames = []
            parents = {}    # Parent for each term
            renamed = []
            ex = 0.0; etot = sdb.entryNum()
            for entry in sdb.entries():
                self.progLog('\r#SPF','Modifying SPF content: %.1f%%' % (ex/etot)); ex += 100.0
                taxon = ''
                parent = ''
                #self.debug(entry)
                for lvl in ['Level_1','Level_2','Level_3','Level_4','Level_5','Level_6','Level_7']:
                    entry[lvl] = string.replace(entry[lvl],'unidentified','unclassified')
                    #entry[lvl] = string.replace(entry[lvl],'Incertae_sedis','Incertae_sedis-%s' % levels[lvl])
                    null = '%s__' % levels[lvl]
                    #self.bugPrint(null)
                    #self.bugPrint(entry[lvl])
                    if entry[lvl] in [null,'Unassigned','unclassified','%sunclassified' % null,'%sunidentified' % null,'%sunculturedfungus' % null,'%sIncertae_sedis' % null,'%sunclassified_sp.' % null]:
                        if not taxon or taxon.endswith('unclassified'): entry[lvl] = '%sunclassified' % null
                        #elif taxon.endswith('unassigned)'): entry[lvl] = '%s%s' % (null,taxon[3:])
                        #elif taxon.endswith('unassigned)'): entry[lvl] = '%s(%s;%s-unassigned)' % (null,string.split(taxon,'(')[1][:-1],levels[lvl])
                        elif taxon.endswith('unassigned)'): entry[lvl] = '%s%s;%s-unassigned)' % (null,taxon[3:][:-1],levels[lvl])
                        else: entry[lvl] = '%s%s(%s-unassigned)' % (null,taxon[3:],levels[lvl])
                    if entry[lvl] in parents:
                        #self.debug(parents[entry[lvl]])
                        if parent in parents[entry[lvl]]: entry[lvl] = parents[entry[lvl]][parent]
                        else:
                            self.bugPrint(entry[lvl])
                            self.bugPrint(parents[entry[lvl]])
                            renamed.append(entry[lvl])
                            newtax = '%s%d' % (entry[lvl],renamed.count(entry[lvl]))
                            self.warnLog('%s had multiple parents (%s & %s) -> %s' % (entry[lvl],string.join(parents[entry[lvl]],'|'),parent,newtax))
                            parents[newtax] = {parent:newtax}
                            parents[entry[lvl]][parent] = newtax
                            entry[lvl] = newtax
                            self.deBug(parents[entry[lvl]])
                    elif parent: parents[entry[lvl]] = {parent:entry[lvl]}
                    parent = entry[lvl]
                    if entry[lvl][3:] == taxon[3:]:
                        if (entry[lvl],taxon) not in dupnames: dupnames.append((entry[lvl],taxon))
                    #self.bugPrint(entry[lvl])
                    taxon = entry[lvl]
                #self.debug(entry)
                #self.debug(parents)
            self.printLog('\r#SPF','Modifying SPF content complete.')
            dupnames.sort()
            for (dupA,dupB) in dupnames: self.warnLog('Duplicate taxa names: %s & %s' % (dupA,dupB))
            ### ~ [2] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.saveToFile(savefields=sdb.list['Fields'][1:])
            ### ~ [3] Compress to different taxonomic levels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            compress = ['Level_1','Level_2','Level_3','Level_4','Level_5','Level_6','Level_7','#']
            dump = compress.pop(-1)
            rules = {'Observation Ids':'list',dump:'str'}
            sdb.dropField('Observation Ids')
            while compress:
                sdb.compress(compress,rules=rules,default='sum',best=[],joinchar='|')
                #if dump == '#':
                sdb.dropField(dump)
                sdb.saveToFile('%s.SPF.%s.%s.spf' % (rje.baseFile(infile),compress[-1],levels[compress[-1]]))
                dump = compress.pop(-1); rules[dump] = 'list'

        except: self.errorLog('Error in rje_misc.JRJSPF()')
#########################################################################################################################
    def aicPaper(self): # Conversion of AIC PATIS tables for Jan 2014 AIC paper.
        '''Conversion of AIC PATIS tables for Jan 2014 AIC paper.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = 'Homo_sapiens.GRCh37.p13'
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db.basefile('aic_paper')
            fields = ['Gene',       # Gene Identifier
                      'ENSG',       # EnsEMBL gene ID
                      'ENST',       # EnsEMBL transcript ID
                      'CDS',        # Length of CDS
                      'UTR5',       # Position of annotated Start Codon (length of 5' UTR)
                      # Annotated start codon data
                      'Codon',      # Annotated Start Codon
                      'Core',       # Core sequence nnnNNNn
                      'Context',    # Start codon plus 15nt flanks
                      'Class',      # Class of transcript start site (Unique, Redundant, Internal)
                      'GC',         # Percentage GC in window +15 - +85
                      'Strength',   # Strength of annotated Start Codon (Good/Bad/NonAUG/Fail)
                      # eORF and tORF data
                      'uSTOP',      # Position of upstream in-frame stop codon
                      'StartPos',   # Position of annotated start codon
                      'GoodPos',    # Position of most 5' in-frame Good start codon (RnnAUGn or nnnAUGG)
                      'eORF',       # Position of most 5' AIC
                      'eType',      # Strength of most 5' AIC
                      'eCore',      # Core sequence of most 5' AIC
                      'eCount',     # Number of 5' AIC (Good/Bad/NonAUG)
                      'eLen',       # Size of max eORF extension (aa)
                      'tORF',       # Position of 3' AIC
                      'tCore',      # Strength of 3' AIC
                      'tCount',     # Number of 3' AIC (Good/Bad/NonAUG)
                      'tType',      # Strength of 3' AIC
                      'tLen',       # Size of min tORF truncation (aa)
                      'RF1',        # Number of start sites in 5' UTR RF1 that would span annotated start
                      'RF2',        # Number of start sites in 5' UTR RF2 that would span annotated start
                      # Verdict
                      'BIOL3050']   # Whether this transcript is potentially suitable for BIOL3050

            ### Fig 1B: Count and % all AUG IC -3/+4 ###
            if False:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ FIG 1B ~~~~~~~~~~~~~~~~~~~~~ #')
                pdb = db.addTable('%s.patis.tdt' % basefile,mainkeys=['Gene','ENST'],name='Fig1B')
                pdb.keepFields(['Gene','ENST','Core'],log=False)
                pdb.addField('-3'); pdb.addField('+4'); pdb.addField('N',evalue=1)
                for entry in pdb.entries():
                    entry['-3'] = entry['Core'][0]
                    entry['+4'] = entry['Core'][-1]
                    if '-' in [entry['-3'],entry['+4']]:
                        self.warnLog('%s' % entry)
                pdb.dropField('Core')
                pdb.compress(['+4','-3'],default='sum')
                pdb.dropFields(['Gene','ENST'])
                pdb.saveToFile()

            ### Fig 1C: counts by efficiency
            if False:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ FIG 1C ~~~~~~~~~~~~~~~~~~~~~ #')
                pdb = db.addTable('%s.patis.tdt' % basefile,mainkeys=['Gene','ENST'],name='Fig1C')
                pdb.keepFields(['Gene','ENST','Strength'],log=False); rules={}
                for cstr in pdb.index('Strength'): pdb.addField('T|%s' % cstr,evalue=0); rules[pdb.fields()[-1]] = 'sum'
                for cstr in pdb.index('Strength'): pdb.addField('G|%s' % cstr,evalue=0); rules[pdb.fields()[-1]] = 'max'
                for entry in pdb.entries(): entry['T|%s' % entry['Strength']] = 1; entry['G|%s' % entry['Strength']] = 1
                pdb.dropField('Strength')
                pdb.compress(['Gene'],rules,'str')
                pdb.dropField('ENST'); pdb.addField('Temp',evalue='!')
                pdb.compress(['Temp'],default='sum')
                pdb.dropField('Gene')
                pdb.reshapeLong('Strength',['G','T'])
                pdb.newKey('Strength')
                pdb.dropField('Temp')
                pdb.saveToFile()

            ### Fig 1D: weak annotated IC with tORFs (1st Type/Length)
            if True:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ FIG 1D ~~~~~~~~~~~~~~~~~~~~~ #')
                pdb = db.addTable('%s.patis.tdt' % basefile,mainkeys=['ENST'],name='Fig1D')
                pdb.keepFields(['ENST','Strength','tType','tCore','tLen'],log=False )
                pdb.dataFormat({'tLen':'int'})
                pdb.dropEntriesDirect('Strength',['Weak','NonAUG'],inverse=True)
                pdb.dropEntriesDirect('tType',[''])
                pdb.dropField('Strength')
                bins = [-1,0,40,80,120,160,200]
                for i in range(len(bins)-1): pdb.addField('%d-%d' % (bins[i]+1,bins[i+1]),evalue=0)
                pdb.addField('%d+' % bins[-1])
                for entry in pdb.entries():
                    entry['tCore'] = entry['tCore'][-4:-1]
                    binned = False
                    for i in range(len(bins)-1):
                        if entry['tLen'] < bins[i+1]: entry[pdb.fields()[i-len(bins)]] = 1; binned = True; break
                    if not binned: entry[pdb.fields()[-1]] = 1
                    self.deBug(entry)
                pdb.compress(['tCore'],default='sum')
                pdb.keepFields(['tCore']+pdb.fields()[-len(bins):])
                pdb.saveToFile()

            ### Fig 2B: all eORFs (1st type/length)
            if True:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ FIG 2B ~~~~~~~~~~~~~~~~~~~~~ #')
                pdb = db.addTable('%s.patis.tdt' % basefile,mainkeys=['ENST'],name='Fig2B')
                pdb.keepFields(['ENST','Class','eType','eCore','eLen'],log=False )
                pdb.dataFormat({'eLen':'int'})
                pdb.dropEntriesDirect('Class',['Internal'])
                pdb.dropEntriesDirect('eType',[''])
                pdb.dropField('Class')
                bins = [-1,0,40,80,120,160,200]
                for i in range(len(bins)-1): pdb.addField('%d-%d' % (bins[i]+1,bins[i+1]),evalue=0)
                pdb.addField('%d+' % bins[-1])
                for entry in pdb.entries():
                    entry['eCore'] = entry['eCore'][-4:-1]
                    binned = False
                    for i in range(len(bins)-1):
                        if entry['eLen'] < bins[i+1]: entry[pdb.fields()[i-len(bins)]] = 1; binned = True; break
                    if not binned: entry[pdb.fields()[-1]] = 1
                    self.deBug(entry)
                pdb.compress(['eCore'],default='sum')
                pdb.keepFields(['eCore']+pdb.fields()[-len(bins):])
                pdb.saveToFile()


        except: self.errorLog('Error in rje_misc(%s)' % self.info['Name'],printerror=True)#,quitchoice=True)
#########################################################################################################################
    def testing(self):  ### Code testing method
        '''Code testing method.'''
        try:
            self.printLog('#PATH','SLiMSuite: %s' % rje.slimsuite)
            raise ValueError("What happens to his text?")
        except: self.errorLog('Error in rje_misc(%s)' % self.info['Name'],printerror=True)#,quitchoice=True)
#########################################################################################################################
    def holger(self):   ### Reformat EHux fasta and extract to TDT
        '''Reformat EHux fasta and extract to TDT.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            best = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=Emihu1_best_transcripts.fasta','autoload=T'])
            all = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=Emihu1_all_transcripts.fasta','autoload=T'])
            db = rje_db.Database(self.log,self.cmd_list+['delimit=,'])
            tdb = db.addEmptyTable('transcripts',['file','name','number','desc','sequence'],keys=['file','name'])
            for seq in best.seqs():
                (name,sequence) = best.getSeq(seq,'tuple')
                tdb.addEntry({'file':'best','name':name,'number':string.split(name,'|')[2],'desc':string.split(name,'|')[3],'sequence':sequence})
            for seq in all.seqs():
                (name,sequence) = all.getSeq(seq,'tuple')
                tdb.addEntry({'file':'all','name':name,'number':string.split(name,'|')[2],'desc':string.split(name,'|')[3],'sequence':sequence})
            tdb.saveToFile('holger_transcripts.csv')

        except: self.errorLog('Error in rje_misc(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def yunan(self):    ### Combine ATG numbers with ANOVA table data and calculate ratios based on normalised data
        '''Combine ATG numbers with ANOVA table data and calculate ratios based on normalised data.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            min_expr = 0.01     # Value to give zero data
            ### ~ [1] ~ Combine data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db.basefile('plantago_atmap')
            pdb = db.addTable('plantago_atg.pme.tdt',name='pme',mainkeys=['EST'])
            ndb = db.addTable('plantago_unatg.nvl.tdt',name='nvl',mainkeys=['EST'])
            db.mergeTables(pdb,ndb,overwrite=True)
            pdb.renameField('EST','Contig')
            adb = db.addTable('../plantago.111202.anova.fdr.tdt',name='anova',mainkeys=['Contig'])
            for field in adb.fields()[0:]:
                if field[0] in ['u','t']: adb.dropField(field)
            mdb = db.joinTables(name='atmap',join=[(pdb,'Contig'),(adb,'Contig')],newkey=['Contig'],cleanup=True,delimit='\t',empties=True,check=True,keeptable=True)
            mdb.saveToFile()
            ### ~ [2] ~ Normalise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            normdb = db.copyTable(mdb,'norm')
            for field in normdb.fields()[0:]:
                if field[0] not in 'ABCD' or field[:2] == 'Df': normdb.dropField(field)
            normdb.addField('Normalised',after='Contig',evalue='Sum')
            normdb.compress(['Normalised'],default='sum',best=[])
            normdb.dropField('Contig')
            normdb.saveToFile()
            ## ~ [2a] ~ Convert value ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.info['Name'] = 'ppm'
            sumdb = normdb.data('Sum')
            mdb.addField('meanA',evalue=0.0)
            mdb.addField('meanB',evalue=0.0)
            mdb.addField('meanC',evalue=0.0)
            mdb.addField('meanD',evalue=0.0)
            mdb.addField('B/A',evalue=0.0)
            mdb.addField('D/C',evalue=0.0)
            mdb.addField('C/A',evalue=0.0)
            mdb.addField('D/B',evalue=0.0)
            ex = 0.0; etot = mdb.entryNum(); absmin = 1e7
            for entry in mdb.entries():
                self.progLog('\r#NORM','Normalising to ppm %.2f%%' % (ex/etot)); ex += 100.0
                for field in normdb.fields()[1:]:
                    entry[field] = 1e6 * float(entry[field]) / sumdb[field]
                    if entry[field] > 0.0: absmin = min(absmin,entry[field])
                    else: entry[field] = 0.01
                    x = field[0]
                    entry['mean%s' % x] += entry[field]
                for x in 'ABCD': entry['mean%s' % x] /= 6
                #B/A, D/C, C/A and D/B.
                entry['B/A'] = entry['meanB'] / entry['meanA']
                entry['D/C'] = entry['meanD'] / entry['meanC']
                entry['C/A'] = entry['meanC'] / entry['meanA']
                entry['D/B'] = entry['meanD'] / entry['meanB']
            self.printLog('\r#NORM','Normalised to ppm: min. value = %s' % absmin)
            self.printLog('\r#ZERO','Zero PPM values replaced with %s' % min_expr)
            mdb.saveToFile()           
            ### ~ [3] ~ Compress on AGI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.info['Name'] = 'at_ppm'
            mdb.dropEntriesDirect('TopHit',['-'])
            mdb.dropFields(['RFX','Trans','GeneLen','HitID','HitDef','HitEVal','eTot','Df.Plant','SumSq.Plant','MeanSq.Plant','F.Plant','Df.CO2','SumSq.CO2','MeanSq.CO2','F.CO2','Df.Interaction','SumSq.Interaction','MeanSq.Interaction','F.Interaction','Df.Residuals','SumSq.Residuals','MeanSq.Residuals','rank'])
            mdb.compress(['TopHit'],{'Contig':'list','B/A':'geomean','D/C':'geomean','C/A':'geomean','D/B':'geomean','E':'max','p.Plant':'min','p.CO2':'min','p.Interaction':'min','p':'min','FDR':'min'},'sum')
            mdb.renameField('B/A','geoB/A')
            mdb.renameField('D/C','geoD/C')
            mdb.renameField('C/A','geoC/A')
            mdb.renameField('D/B','geoD/B')
            mdb.addField('B/A',evalue=0.0)
            mdb.addField('D/C',evalue=0.0)
            mdb.addField('C/A',evalue=0.0)
            mdb.addField('D/B',evalue=0.0)
            mdb.addField('ContigNum',evalue=0)
            mdb.addField('AT','TopHit',evalue=0)
            for entry in mdb.entries():
                entry['B/A'] = entry['meanB'] / entry['meanA']
                entry['D/C'] = entry['meanD'] / entry['meanC']
                entry['C/A'] = entry['meanC'] / entry['meanA']
                entry['D/B'] = entry['meanD'] / entry['meanB']
                entry['ContigNum'] = len(string.split(entry['Contig'],';'))
                entry['AT'] = string.split(entry['TopHit'],'__')[-1]
            ### ~ [5] ~ Reorganise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.list['Fields'].append(mdb.list['Fields'].pop(0))
            mdb.newKey(['AT'],startfields=True)
            mdb.dropFields(['TopHit'])
            mdb.dataFormat(reformat={'p':'num'})
            mdb.rankField('p',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False)
            for entry in mdb.entries(): entry['Rank'] = mdb.entryNum() - entry['Rank'] + 1
            mdb.addField('newFDR','FDR',evalue=0)
            for entry in mdb.entries():
                try: entry['newFDR'] = (mdb.entryNum() * entry['p']) / entry['Rank']
                except: entry['newFDR'] = '!'
            mdb.saveToFile()           

            
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def shapeSAGE(self):    ### Reshape SuperSAGE data for R analysis
        '''Reshape SuperSAGE data for R analysis.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            #self.basefile('DT.annotated.MT_tag_comp.14.ppm.Tag')
            # Tag	WTL1	WTL2	WTL3	WTM1	WTM2	WTM3	WTH1	WTH2	WTH3	L1L1	L1L2	L1L3	L1M1	L1M2	L1M3	L1H1	L1H2	L1H3	L2L1	L2L2	L2L3	L2M1	L2M2	L2M3	L2H1	L2H2	L2H3	Count	ECount	Seq	SeqN
            tdb = db.addTable('%s.tdt' % self.basefile(),name='Long',mainkeys=['Tag'])
            tdb.dropFields(['Count','ECount','Seq','SeqN','AllRep','MaxCount'])
            for field in tdb.fields()[1:]: tdb.renameField(field,'Expression|%s' % field)
            tdb.reshapeLong('Experiment',reshape=['Expression'])
            tdb.addField('Strain'); tdb.addField('Light')
            for entry in tdb.entries():
                entry['Strain'] = entry['Experiment'][:2]
                entry['Light'] = {'L':'Low','M':'Med','H':'High'}[entry['Experiment'][2]]
            tdb.saveToFile()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def dmi(self):   ### Make 3DID DMI benchmarks
        '''Make 3DID DMI benchmarks.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            base = db.info['Basefile'] = '3did.DMI'
            dmi = db.addTable('3did.DMI.csv',name='DMI',mainkeys=['ID'])
            for entry in dmi.entries(): entry['PDB'] = entry['PDB'].upper()
            dmi.makeField('#PDB#:#dch#','DomPDB'); dmi.makeField('#PDB#:#pch#','MotPDB')
            pdb = db.addTable('3did.DMI.pdb.ens_HUMAN.mapping.tdt',name='PDB',mainkeys=['Query'])
            #pdb.dropEntriesDirect('Hit',['None'])
            for entry in pdb.entries(): entry['Query'] = string.split(entry['Query'],'|')[0]
            xref = db.addTable('../../HGNC/genemap.111005.data.tdt',name='HGNC',mainkeys=['Gene'])
            ### ~ [1] ~ Join Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            temp = db.joinTables(name='Temp',join=[(dmi,'DomPDB'),(pdb,'Query',['Hit'])],newkey=['ID'],empties=False)
            temp.renameField('Hit','DomSeq')
            temp2 = db.joinTables(name='Temp2',join=[(temp,'DomSeq'),(xref,'EnsLoci',['Gene'])],newkey=['ID','Gene'],empties=False)
            map = db.joinTables(name='PDBMap',join=[(temp2,'MotPDB'),(pdb,'Query',['Hit'])],newkey=['ID','Gene'],empties=False)
            map.renameField('Hit','MotSeq'); map.renameField('Gene','Hub')
            map.saveToFile()
            ### ~ [2] ~ Make sequence query files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map.dropEntriesDirect('DomSeq',['None'])
            map.dropEntriesDirect('MotSeq',['None'])
            map.makeField('#MotSeq#','qry')
            for entry in map.entries(): entry['qry'] = string.split(entry['qry'],'_')[-1]
            map.makeField('#peptide#.#Hub#.#qry#','dataset')
            map.dataFormat({'pstart':'int','pend':'int'})
            ensloci = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=../../EnsEMBL/ens_HUMAN.loci.fas','autoload=T','seqmode=file'])
            ensdict = ensloci.seqNameDic()
            SEQOUT = open('%s.pseq.fas' % base,'w')
            flank = 5
            ppidir = '../../Pingu/PPIFas/'
            rje.mkDir(self,'3DID_PPIFas/',log=True)
            for dataset in map.index('dataset'):
                pseq = ''; qry = ''; hub = ''; ex = 0
                for entry in map.indexEntries('dataset',dataset):
                    #while len(seq) < entry['pend']: seq += 'x'
                    qry = entry['MotSeq']; hub = entry['Hub']; ex += 1
                pseq = string.join(map.indexDataList('dataset',dataset,'pseq'))
                self.printLog('#DSET','%s: %d entries' % (dataset,ex))
                SEQOUT.write('>%s %s\n%s\n' % (dataset,qry,pseq))
                ppifas = ppidir + hub + '.ppi.fas'
                outfas = '3DID_PPIFas/%s.ppi.fas' % dataset
                if not os.path.exists(ppifas):
                    self.printLog('#NOPPI','No PPI file for hub %s' % hub)
                    continue
                ppiseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % ppifas,'autoload=T','seqnr=F','accnr=F','usecase=T'])
                ppidict = ppiseq.seqNameDic()
                if qry in ppidict: qseq = ppidict[qry]
                elif qry in ensdict:
                    (name,sequence) = ensloci.getSeq(ensdict[qry],'tuple')
                    self.printLog('#QRY','Adding Query %s to %s' % (qry,dataset))
                    qseq = ppiseq._addSeq(name,sequence)                    
                else:
                    self.printLog('#NOQRY','Cannot find query sequence %s' % qry)
                    continue
                if qseq in ppiseq.seq: ppiseq.seq.remove(qseq)
                else: self.printLog('#QRY','Adding Query %s to %s' % (qry,dataset))
                newseq = qseq.info['Sequence'].lower()
                uppseq = qseq.info['Sequence'].upper()
                for frag in map.indexDataList('dataset',dataset,'pseq'):
                    x = 0
                    self.deBug(newseq)
                    while uppseq.find(frag,x) >= 0:
                        i = uppseq.find(frag,x)
                        newseq = newseq[:max(i-flank,0)] + uppseq[max(i-flank,0):i+len(frag)+flank] + newseq[i+len(frag)+flank:]
                        self.deBug(newseq)
                        x = i + 1
                PPI = open(outfas,'w')
                PPI.write('>%s\n%s\n' % (qseq.info['Name'],newseq))
                for seq in ppiseq.seq: PPI.write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence'].upper()))
                PPI.close()
                self.printLog('#FAS','%s: %d seq' % (outfas,ppiseq.seqNum()+1))
                    
            SEQOUT.close()

            
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def sfStick(self):   ### Stick together SF Files
        '''Stick together SF Files.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.stat['SigCut'] = 0.05
            for cmd in self.cmd_list: self._cmdReadList(cmd,'stat',['SigCut'])
            sigcut = self.stat['SigCut']
            db = rje_db.Database(self.log,['basefile=qsf_v_sf.%s' % sigcut]+self.cmd_list)
            ipath = './'
            table = db.addTable('combined.csv',name='Full',mainkeys=['Dataset','RunID','Rank','Pattern'])
            elm = db.addTable('ELM.090608.tdt',name='ELM',mainkeys=['ELM'])
            ### ~ [1] ~ Reduce ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table.dropFields(['Dataset','RunID','SeqNum','UPNum','Rank','Pattern','Sig','Cloud','CloudSeq','CloudUP'],inverse=True,log=True)
            table.addField('Species','Dataset','All')
            table.addField('ELM','Dataset','')
            table.addField('Query','ELM','')
            table.addField('Prog','ELM','SF')
            table.addField('Method','ELM','')
            table.addField('ELM_Pattern','ELM','')
            table.makeField('[#CloudSeq#|#CloudUP#]','Support')
            nsx = 0
            for dkey in table.dataKeys():
                entry = table.data(dkey)
                [entry['ELM'],entry['Query'],dset] = string.split(entry['Dataset'],'.')
                if 'QSF' in entry['RunID']: entry['Prog'] = 'QSF'
                if 'Subset' in entry['Dataset']: entry['Species'] = 'Subset'
                if '_' not in dset: entry['Dataset'] = '%s_whole' % entry['Dataset']
                entry['Method'] = '%s_%s' % (entry['Prog'],string.split(entry['Dataset'],'_')[-1])
                if string.atof(entry['Sig']) > sigcut:
                    if entry['Rank'] not in ['0','1']: table.dict['Data'].pop(dkey); continue
                    if entry['Rank'] == '1': entry['Rank'] = '0'; entry['Pattern'] = '-'; entry['Sig'] = '0.1'; nsx += 1
                if entry['Rank'] == '0': entry['Cloud'] = '1'; entry['Support'] = '-'
                try: entry['ELM_Pattern'] = elm.data(entry['ELM'])['Definition']
                except: entry['ELM_Pattern'] = '?!'
            self.printLog('#NS','%d NS entries cleared (p > %s)' % (nsx,sigcut))
            table.dropEntriesDirect('Query','NoQry')
            table.makeField('#Prog#|#Dataset#|#Cloud#','FCloud')
            table.compress(['FCloud'],rules={'Rank':'min','Sig':'min'},default='mean',best=['Sig','Rank'])
            table.renameField('Dataset','File')
            #table = db.joinTables('Full',join=[(table,'ELM'),(elm,'ELM')],newkey=['FCloud'],empties=True)
            table.saveToFile()
            ### ~ [2] ~ Reshape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table.info['Name'] = 'Cloud'
            table.makeField('#ELM#|#Query#','Dataset','File')
            table.makeField('#Dataset#|#Cloud#','DCloud','Dataset')
            table.makeField('#SeqNum# (#UPNum#)','N')
            table.newKey(['DCloud','Species','Method'],startfields=False)
            roc = db.copyTable(table,'Verdict')
            table.dropFields(['DCloud','ELM_Pattern','Species','Cloud','Method','N','Pattern','Sig','Support'],inverse=True,log=True)
            table.list['Fields'] = ['DCloud','ELM_Pattern','Cloud','Species','Method','N','Pattern','Sig','Support']
            summary = db.copyTable(table,'Summary')
            table.reshapeWide('Method',reshape=['N','Pattern','Sig','Support'])
            #table.renameField('N|All.QSF_flank5','N|All'); table.renameField('N|Subset.QSF_flank5','N|Subset')
            table.renameField('N|QSF_flank5','N')
            for field in table.fields()[0:]:
                if field[:2] == 'N|': table.dropField(field)
            newfields = ['DCloud','ELM_Pattern','Cloud','Species','N']
            for reg in ['whole','win100','flank5']:
                for method in ['SF','QSF']:
                    for score in ['Pattern','Sig','Support']: newfields.append('%s|%s_%s' % (score,method,reg))
            table.list['Fields'] = newfields
            table.saveToFile()
            
            table = summary
            table.makeField('#Pattern# (p=#Sig#) #Support#','SLiM')
            table.dropFields(['Pattern','Sig','Support'])
            for entry in table.entries():
                if entry['SLiM'] == '- (p=0.1) -': entry['SLiM'] = '-'
            table.reshapeWide('Method',reshape=['N','SLiM'])
            table.renameField('N|QSF_flank5','N')
            for field in table.fields()[0:]:
                if field[:2] == 'N|': table.dropField(field)
                if field[:4] == 'SLiM': table.renameField(field,field[5:])
            newfields = ['DCloud','ELM_Pattern','Cloud','Species','N']
            for reg in ['whole','win100','flank5']:
                for method in ['SF','QSF']: newfields.append('%s_%s' % (method,reg))
            table.list['Fields'] = newfields
            table.saveToFile()

            table.info['Name'] = 'Wide'
            table.reshapeWide('Species',reshape=table.fields()[4:])
            newfields = ['DCloud','ELM_Pattern','Cloud']
            for sub in ['All','Subset']:
                newfields.append('N|%s' % sub)
                for reg in ['whole','win100','flank5']:
                    for method in ['SF','QSF']: newfields.append('%s_%s|%s' % (method,reg,sub))
            table.list['Fields'] = newfields
            table.saveToFile()


            roc.dropFields(['DCloud','ELM','ELM_Pattern','Species','Cloud','Method','N','SeqNum','Pattern','Sig','Support'],inverse=True,log=True)
            roc.list['Fields'] = ['DCloud','ELM','ELM_Pattern','Cloud','Species','Method','N','SeqNum','Pattern','Sig','Support']
            tpfile = 'qsf_v_sf.TP.csv'
            if os.path.exists(tpfile):
                tp = db.addTable(tpfile,name='TP',mainkeys=['ELM','Pattern'])
                tp.dropEntriesDirect('Verdict',['OT','TP','FP','FN'],inverse=True)
            else: tp = db.addEmptyTable('TP',['ELM','Pattern','Verdict'],['ELM','Pattern'])
            roc.addField('Verdict')
            self.deBug(roc.fields())
            for entry in roc.entries():
                self.deBug(entry)
                vkey = tp.makeKey(entry)
                self.deBug(vkey)
                if vkey not in tp.data():
                    if entry['Pattern'] == '-': tp.addEntry({'ELM':entry['ELM'],'Pattern':entry['Pattern'],'Verdict':'FN'})
                    else:
                        try:
                            tp.addEntry({'ELM':entry['ELM'],'Pattern':entry['Pattern'],
                                         'Verdict':rje.choice('Classify %s vs %s - %s :' % (entry['Pattern'],entry['ELM'],entry['ELM_Pattern']),'?',True).upper()})
                        except: break
                entry['Verdict'] = tp.data()[vkey]['Verdict']
            tp.saveToFile(tpfile)
            perf = db.copyTable(roc,'QvS')
            roc.saveToFile()
            roc.info['Name'] = 'ROC'
            for v in ['OT','TP','FP','FN']: roc.addField(v,evalue=0)
            for entry in roc.entries(): entry[entry['Verdict']] += 1
            roc.dropFields(['Pattern','Sig','Support','Verdict','Cloud'])
            roc.compress(['ELM','Species','Method'],{'SeqNum':'mean'},default='sum')
            roc.addField('SN'); roc.addField('SP'); roc.addField('OTT'); roc.addField('OTF')
            for entry in roc.entries():
                entry['SN'] = 100.0 * float(entry['TP']) / max(int(entry['SeqNum']),(entry['TP'] + entry['FN']))
                if entry['TP'] or entry['FP']: entry['SP'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'])
                else: entry['SP'] = ''
                if entry['TP'] or entry['FP'] or entry['OT']:
                    entry['OTT'] = 100.0 * float(entry['TP']+entry['OT']) / (entry['TP'] + entry['FP'] + entry['OT'])
                    entry['OTF'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'] + entry['OT'])
                else: entry['OTT'] = ''; entry['OTF'] = ''           
            roc.addField('QRegion','Method')
            for entry in roc.entries():
                (entry['Method'],entry['QRegion']) = string.split(entry['Method'],'_')
            roc.newKey(['ELM','Species','Method','QRegion'])
            roc.saveToFile()
            
            roc.info['Name'] = 'ROCmean'
            roc.compress(['Species','Method','QRegion'],default='mean')
            roc.dropFields(['ELM','ELM_Pattern','N','DCloud'])
            roc.saveToFile()

            perf.dropEntriesDirect('Cloud',['1'],inverse=True)
            perf.dropFields(['Cloud','Pattern','Support'])
            perf.addField('QRegion','Method')
            for entry in perf.entries():
                (entry['Method'],entry['QRegion']) = string.split(entry['Method'],'_')
            perf.newKey(['DCloud','Species','Method','QRegion'],startfields=False)
            perf.reshapeWide('Method',reshape=['Sig','Verdict'])
            perf.addField('QvS',evalue=0); perf.addField('QSF',evalue=0); perf.addField('SF',evalue=0); perf.addField('Tie',evalue=0)
            qsf = 'Verdict|QSF'; sf = 'Verdict|SF'
            for entry in perf.entries():
                if entry[qsf] == 'TP' and float(entry['Sig|QSF']) < float(entry['Sig|SF']): entry['QvS'] = 1; entry['QSF'] = 1
                elif entry[qsf] == 'TP' and entry[sf] != 'TP': entry['QvS'] = 1; entry['QSF'] = 1
                elif entry[sf] == 'TP' and float(entry['Sig|QSF']) > float(entry['Sig|SF']): entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'TP' and entry[qsf] != 'TP': entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'TP' and entry[qsf] == 'TP': entry['Tie'] += 1
                elif entry[qsf] == 'FP' and entry[sf] != 'FP': entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'FP' and entry[qsf] != 'FP': entry['QvS'] = 1; entry['QSF'] = 1
            perf.dropFields(['Sig|SF','Sig|QSF',qsf,sf])
            perf.compress(['ELM','Species','QRegion'],{'SeqNum':'mean'},default='sum')
            perf.addField('Null')
            for entry in perf.entries():
                entry['QvS'] = 100.0 * entry['QvS'] / float(entry['SeqNum'])
                entry['QSF'] = 100.0 * entry['QSF'] / float(entry['SeqNum'])
                entry['SF'] = 100.0 * entry['SF'] / float(entry['SeqNum'])
                entry['Tie'] = 100.0 * entry['Tie'] / float(entry['SeqNum'])
                entry['Null'] = 100.0 - (entry['QSF'] + entry['SF'] + entry['Tie'])
            perf.saveToFile()
            perf.info['Name'] = 'QvSmean'
            perf.compress(['Species','QRegion'],default='mean')
            perf.dropFields(['Species','QRegion','QSF','SF','Tie','Null'],inverse=True)
            perf.saveToFile()

                
            
                             
            ##
            return
            
            
                
            
            return            
            table.dropField('Keep'); table.dropField('Rank');                
            ### ~ [2] ~ Join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sf.makeField('#ELM#|#Species#','Joiner'); sfq.makeField('#ELM#|#Species#','Joiner')
            self.deBug(rje.sortKeys(sf.index('Joiner')))
            self.deBug(rje.sortKeys(sfq.index('Joiner')))
            sf = db.joinTables(name='Temp',join=[(sfq,'Joiner',['Query','Dataset']),(sf,'Joiner')],newkey=['Dataset'])
            sf.dropField('SF_Dataset'); sf.dropField('Joiner')
            sf.saveToFile('tempcheck.tdt')
            top = db.joinTables(name='TopRanks',join=[(sf,'Dataset'),(sfq,'Dataset',['Pattern','Sig','Occ','Support','UP','Cloud','CloudSeq']),(qsf,'Dataset',['Pattern','Sig','Occ','Support','UP','Cloud','CloudSeq'])],newkey=['Dataset'])
            top.saveToFile('qsf_versus_sf.tdt')
            

        except:
            self.log.errorLog('Error in rje_misc.sfStick(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def sfStick2(self):   ### Stick together SF Files
        '''Stick together SF Files.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.stat['SigCut'] = 0.05
            for cmd in self.cmd_list: self._cmdReadList(cmd,'stat',['SigCut'])
            sigcut = self.stat['SigCut']
            db = rje_db.Database(self.log,['basefile=qsf_3did.%s' % sigcut]+self.cmd_list)
            ipath = './'
            table = db.addTable('qsf_3did.csv',name='Full',mainkeys=['Dataset','RunID','Rank','Pattern'])
            elm = db.addTable('ELM.090608.tdt',name='ELM',mainkeys=['ELM'])
            ### ~ [1] ~ Reduce ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table.dropFields(['Dataset','RunID','SeqNum','UPNum','Rank','Pattern','Sig','Cloud','CloudSeq','CloudUP'],inverse=True,log=True)
            table.addField('Hub','Dataset','All')
            table.addField('ELM','Dataset','')
            table.addField('Query','ELM','')
            table.addField('Prog','ELM','SF')
            table.addField('ELM_Pattern','ELM','')
            table.makeField('[#CloudSeq#|#CloudUP#]','Support')
            nsx = 0
            for dkey in table.dataKeys():
                entry = table.data(dkey)
                [entry['ELM'],entry['Hub'],entry['Query'],dset] = string.split(entry['Dataset'],'.')
                if 'QSF' in entry['RunID']: entry['Prog'] = 'QSF'
                if string.atof(entry['Sig']) > sigcut:
                    if entry['Rank'] not in ['0','1']: table.dict['Data'].pop(dkey); continue
                    if entry['Rank'] == '1': entry['Rank'] = '0'; entry['Pattern'] = '-'; entry['Sig'] = '0.1'; nsx += 1
                if entry['Rank'] == '0': entry['Cloud'] = '1'; entry['Support'] = '-'
                try: entry['ELM_Pattern'] = elm.data(entry['ELM'])['Definition']
                except: entry['ELM_Pattern'] = '?!'
            self.printLog('#NS','%d NS entries cleared (p > %s)' % (nsx,sigcut))
            table.makeField('#Prog#|#Dataset#|#Cloud#','FCloud')
            table.compress(['FCloud'],rules={'Rank':'min','Sig':'min'},default='mean',best=['Sig','Rank'])
            table.renameField('Dataset','File')
            #table = db.joinTables('Full',join=[(table,'ELM'),(elm,'ELM')],newkey=['FCloud'],empties=True)
            table.saveToFile()
            ### ~ [2] ~ Reshape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table.info['Name'] = 'Cloud'
            table.makeField('#ELM#|#Hub#|#Query#','Dataset','File')
            table.makeField('#Dataset#|#Cloud#','DCloud','Dataset')
            table.makeField('#SeqNum# (#UPNum#)','N')
            table.newKey(['DCloud','Prog'],startfields=False)
            roc = db.copyTable(table,'Verdict')
            table.dropFields(['DCloud','ELM_Pattern','Cloud','Prog','N','Pattern','Sig','Support'],inverse=True,log=True)
            table.list['Fields'] = ['DCloud','ELM_Pattern','Cloud','Prog','N','Pattern','Sig','Support']
            summary = db.copyTable(table,'Summary')
            table.reshapeWide('Prog',reshape=['N','Pattern','Sig','Support'])
            #table.renameField('N|All.QSF_flank5','N|All'); table.renameField('N|Subset.QSF_flank5','N|Subset')
            table.renameField('N|QSF','N')
            for field in table.fields()[0:]:
                if field[:2] == 'N|': table.dropField(field)
            newfields = ['DCloud','ELM_Pattern','Cloud','N']
            for prog in ['SF','QSF']:
                for score in ['Pattern','Sig','Support']: newfields.append('%s|%s' % (score,prog))
            table.list['Fields'] = newfields
            table.saveToFile()
            
            table = summary
            table.makeField('#Pattern# (p=#Sig#) #Support#','SLiM')
            table.dropFields(['Pattern','Sig','Support'])
            for entry in table.entries():
                if entry['SLiM'] == '- (p=0.1) -': entry['SLiM'] = '-'
            table.reshapeWide('Prog',reshape=['N','SLiM'])
            table.renameField('N|QSF','N')
            for field in table.fields()[0:]:
                if field[:2] == 'N|': table.dropField(field)
                if field[:4] == 'SLiM': table.renameField(field,field[5:])
            newfields = ['DCloud','ELM_Pattern','Cloud','N']
            for method in ['SF','QSF']: newfields.append('%s' % (method))
            table.list['Fields'] = newfields
            table.saveToFile()

            roc.dropFields(['DCloud','ELM','ELM_Pattern','Cloud','Prog','N','SeqNum','Pattern','Sig','Support'],inverse=True,log=True)
            roc.list['Fields'] = ['DCloud','ELM','ELM_Pattern','Cloud','Prog','N','SeqNum','Pattern','Sig','Support']
            tpfile = 'qsf_3did.TP.csv'
            if os.path.exists(tpfile):
                tp = db.addTable(tpfile,name='TP',mainkeys=['ELM','Pattern'])
                tp.dropEntriesDirect('Verdict',['OT','TP','FP','FN'],inverse=True)
            else: tp = db.addEmptyTable('TP',['ELM','Pattern','Verdict'],['ELM','Pattern'])
            roc.addField('Verdict')
            self.deBug(roc.fields())
            for entry in roc.entries():
                self.deBug(entry)
                vkey = tp.makeKey(entry)
                self.deBug(vkey)
                if vkey not in tp.data():
                    if entry['Pattern'] == '-': tp.addEntry({'ELM':entry['ELM'],'Pattern':entry['Pattern'],'Verdict':'FN'})
                    else:
                        try:
                            tp.addEntry({'ELM':entry['ELM'],'Pattern':entry['Pattern'],
                                         'Verdict':rje.choice('Classify %s vs %s - %s :' % (entry['Pattern'],entry['ELM'],entry['ELM_Pattern']),'?',True).upper()})
                        except: break
                entry['Verdict'] = tp.data()[vkey]['Verdict']
            tp.saveToFile(tpfile)
            perf = db.copyTable(roc,'QvS')
            roc.saveToFile()
            roc.info['Name'] = 'ROC'
            for v in ['OT','TP','FP','FN']: roc.addField(v,evalue=0)
            for entry in roc.entries(): entry[entry['Verdict']] += 1
            roc.dropFields(['Pattern','Sig','Support','Verdict','Cloud'])
            roc.compress(['ELM','Prog'],{'SeqNum':'mean'},default='sum')
            roc.addField('SN'); roc.addField('SP'); roc.addField('OTT'); roc.addField('OTF')
            for entry in roc.entries():
                entry['SN'] = 100.0 * float(entry['TP']) / max(int(entry['SeqNum']),(entry['TP'] + entry['FN']))
                if entry['TP'] or entry['FP']: entry['SP'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'])
                else: entry['SP'] = ''
                if entry['TP'] or entry['FP'] or entry['OT']:
                    entry['OTT'] = 100.0 * float(entry['TP']+entry['OT']) / (entry['TP'] + entry['FP'] + entry['OT'])
                    entry['OTF'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'] + entry['OT'])
                else: entry['OTT'] = ''; entry['OTF'] = ''           
            roc.saveToFile()
            
            roc.info['Name'] = 'ROCmean'
            roc.compress(['Prog'],default='mean')
            roc.dropFields(['ELM','ELM_Pattern','N','DCloud'])
            roc.saveToFile()

            perf.dropEntriesDirect('Cloud',['1'],inverse=True)
            perf.dropFields(['Cloud','Pattern','Support'])
            perf.newKey(['DCloud','Prog'],startfields=False)
            perf.reshapeWide('Prog',reshape=['Sig','Verdict'])
            perf.addField('QvS',evalue=0); perf.addField('QSF',evalue=0); perf.addField('SF',evalue=0); perf.addField('Tie',evalue=0)
            qsf = 'Verdict|QSF'; sf = 'Verdict|SF'
            for entry in perf.entries():
                if entry[qsf] == 'TP' and float(entry['Sig|QSF']) < float(entry['Sig|SF']): entry['QvS'] = 1; entry['QSF'] = 1
                elif entry[qsf] == 'TP' and entry[sf] != 'TP': entry['QvS'] = 1; entry['QSF'] = 1
                elif entry[sf] == 'TP' and float(entry['Sig|QSF']) > float(entry['Sig|SF']): entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'TP' and entry[qsf] != 'TP': entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'TP' and entry[qsf] == 'TP': entry['Tie'] += 1
                elif entry[qsf] == 'FP' and entry[sf] != 'FP': entry['QvS'] = -1; entry['SF'] = 1
                elif entry[sf] == 'FP' and entry[qsf] != 'FP': entry['QvS'] = 1; entry['QSF'] = 1
            perf.dropFields(['Sig|SF','Sig|QSF',qsf,sf])
            perf.compress(['ELM','Prog'],{'SeqNum':'mean'},default='sum')
            perf.addField('Null')
            for entry in perf.entries():
                entry['QvS'] = 100.0 * entry['QvS'] / float(entry['SeqNum'])
                entry['QSF'] = 100.0 * entry['QSF'] / float(entry['SeqNum'])
                entry['SF'] = 100.0 * entry['SF'] / float(entry['SeqNum'])
                entry['Tie'] = 100.0 * entry['Tie'] / float(entry['SeqNum'])
                entry['Null'] = 100.0 - (entry['QSF'] + entry['SF'] + entry['Tie'])
            perf.saveToFile()
            perf.info['Name'] = 'QvSmean'
            perf.compress(['Prog'],default='mean')
            perf.dropFields(['Prog','QSF','SF','Tie','Null'],inverse=True)
            perf.saveToFile()

                
            
                             
            ##
            return
            
            
                
            
            return            
            table.dropField('Keep'); table.dropField('Rank');                
            ### ~ [2] ~ Join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sf.makeField('#ELM#|#Species#','Joiner'); sfq.makeField('#ELM#|#Species#','Joiner')
            self.deBug(rje.sortKeys(sf.index('Joiner')))
            self.deBug(rje.sortKeys(sfq.index('Joiner')))
            sf = db.joinTables(name='Temp',join=[(sfq,'Joiner',['Query','Dataset']),(sf,'Joiner')],newkey=['Dataset'])
            sf.dropField('SF_Dataset'); sf.dropField('Joiner')
            sf.saveToFile('tempcheck.tdt')
            top = db.joinTables(name='TopRanks',join=[(sf,'Dataset'),(sfq,'Dataset',['Pattern','Sig','Occ','Support','UP','Cloud','CloudSeq']),(qsf,'Dataset',['Pattern','Sig','Occ','Support','UP','Cloud','CloudSeq'])],newkey=['Dataset'])
            top.saveToFile('qsf_versus_sf.tdt')
            

        except:
            self.log.errorLog('Error in rje_misc.sfStick(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def kroc(self): ### Generate Specificity and Sensitivity scores for ROC calculations
        basefile='kroc'
        kdb = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % basefile])
        skdb = None
        for sigcut in [0.05,0.01,0.005,0.001,0.0001]:
            if not os.path.exists('%s.roc_mean-%s.tdt' % (basefile,sigcut)):
                #self.rocData(sigcut,basefile)
                self.roc3DID(sigcut,basefile)
            if skdb: kdb.mergeTables(skdb,kdb.addTable('%s.roc_mean-%s.tdt' % (basefile,sigcut),mainkeys=['Method','Sig'],name='%s' % sigcut))
            else: skdb = kdb.addTable('%s.roc_mean-%s.tdt' % (basefile,sigcut),mainkeys=['Method','Sig'],name='roc')
        ##program <- c("QSLiM", "SLiM")
        #cut_off <- c(0.05, 0.01, 0.001)
        #masking <- c("disorder", "conservation", "both", "none")
        #flank <- c("whole", "flank5", "win100", "win300")
        skdb.addField('Flank','Method'); skdb.addField('Masking','Method'); skdb.addField('Program','Method')
        skdb.renameField('SN','Sens'); skdb.renameField('SP','Spec');
        skdb.addField('Sort','Flank')
        remsort = False
        for entry in skdb.entries():
            try:
                [entry['Masking'],entry['Flank'],entry['Program']] = string.split(entry['Method'],'.')
                if entry['Flank'] == 'All': entry['Flank'] = 'whole'
                entry['Sort'] = ['site',"flank5", "win100", "win300",'whole'].index(entry['Flank'])
                entry['Masking'] = {'both':'both','no':'none','dis':'disorder','cons':'conservation'}[entry['Masking'][:-4]]
            except:
                remsort = True
                [entry['Masking'],entry['Program']] = string.split(entry['Method'],'.')
                entry['Masking'] = {'both':'both','no':'none','dis':'disorder','cons':'conservation'}[entry['Masking'][:-4]]
            entry['Program'] = {'SF':'SLiM','QSF':'QSLiM'}[entry['Program']]
        if remsort: skdb.dropFields(['Sort','Flank'])
        skdb.dropField('Method')
        if 'Sort' in skdb.fields(): skdb.newKey(['Program','Masking','Sort','Sig'])
        else: skdb.newKey(['Program','Masking','Sig']); skdb.fillBlanks(blank=0,fields=['Spec','OTT','OTF'],fillempty=True)
        skdb.saveToFile()

        skdb.dropFields(['TP','OT','FP','FN'])
        skdb.dataFormat({'Sens':'num','Spec':'num','OTT':'num','OTF':'num'})
        skdb.reshapeWide('Program',['Sens','Spec','OTT','OTF'])
        skdb.info['Name'] = 'comp'

        for stat in ['Sens','Spec','OTT','OTF']:
            skdb.makeField('%s|QSLiM-%s|SLiM' % (stat,stat),stat)

            #skdb.addField(stat)
            #for entry in skdb.entries(): entry[stat] = entry['%s|QSLiM'] - entry['%s|SLiM' % (stat,stat),stat)


        skdb.saveToFile()
#########################################################################################################################
    def roc3DID(self,sigcut=0.05,basefile='test'):   ### Generate Specificity and Sensitivity scores for ROC calculations
        '''
        Generate Specificity and Sensitivity scores for ROC calculations.
        >> filename:str = Name of full data file
        >> sigcut [0.05] = Significance cut-off
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Setup rje_db.Database object
            db = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % basefile])

            if os.path.exists('%s.full.tdt' % basefile):
                resdb = db.addTable('%s.full.tdt' % basefile,mainkeys=['Method','Dataset','Rank','Pattern'],name='full')
            else:
                # Load results data into 'full' table
                uniquekey = ['dataset','runid','rank','pattern']    #Replace!
                resdb = db.addTable('All_3DID_SLiM.tdt',mainkeys=uniquekey,name='full')
                resdb.addField('Program',evalue='SF')
                resdb.newKey(['dataset','runid','rank','pattern','Program']) # Make new key for table  
                qdb = db.addTable('All_3DID_QSLiM.tdt',mainkeys=uniquekey,name='qsf')
                qdb.addField('Program',evalue='QSF')
                qdb.newKey(['dataset','runid','rank','pattern','Program']) # Make new key for table  
                db.mergeTables(resdb,qdb,overwrite=True,matchfields=True)

                resdb.dropFields(['masking','elm'])
                resdb.renameField('runid', 'Masking')
                resdb.renameField('pos', 'Verdict')
                resdb.renameField('seqnum', 'SeqNum')
                resdb.renameField('rank', 'Rank')
                resdb.renameField('pattern', 'Pattern')
                resdb.renameField('sig', 'Sig')
                resdb.renameField('cloud', 'Cloud')
                resdb.addField('ELM'); resdb.addField('Hub'); resdb.addField('Query'); resdb.addField('Flank')
                for entry in resdb.entries():
                    [entry['ELM'],entry['Hub'],entry['Query'],entry['Flank']] = string.split(entry['dataset'],'.')
                    entry['Flank'] = string.split(entry['Flank'],'_')[-1]

                # Reduce to useful data
                # See also code at SBSBINF/Tools/svn/docs/dev_pydocs.html#libraries-rje_db
                # Fields that we want - Use resdb.renameField(oldname, newname) to rename if needed
                resdb.makeField('#Masking#.#Program#','Method')
                resdb.makeField('#ELM#|#Hub#|#Query#','Dataset')
                resdb.newKey(['Method','Dataset','Rank','Pattern']) # Make new key for table  
                useful_fields = ['Dataset','Method','ELM','Query','SeqNum',       # Dataset fields
                                 'Rank','Pattern','Sig','Cloud',                        # Results fields
                                 'Verdict'                                              # TP/FP/OT/FN
                                 ]
                resdb.dropFields(useful_fields,inverse=True)    # Keep only above fields
                resdb.list['Fields'] = ['Method','Dataset','ELM','Query','SeqNum','Rank','Sig','Pattern','Cloud','Verdict']
                resdb.saveToFile()
                for field in resdb.list['Fields'][:2]: resdb.index(field)

            #Define the exact unique method
            # Method = Masking.flank.program (i.e Disorder.whole.QSLiMFinder)
            resdb.dropFields(['Masking','Flank','Program'])     # These only need to split back up for output (maybe)

            # Reformat data (i.e. numbers as numbers)
            resdb.dataFormat({'Sig':'num','SeqNum':'int','Rank':'int'})

            resdb.info['Name'] = 'roc_qry-%s' % sigcut

            # Add Significance filtering
            for entry in resdb.entries()[0:]:
                if entry['Sig'] > sigcut:
                    if entry['Rank'] > 1: resdb.dropEntry(entry)
                    else: entry['Pattern'] = '-'; entry['Verdict'] = 'FN'; entry['Rank'] = 0; entry['Cloud'] = ''
            resdb.dropField('Sig')

            # Convert ratings into numerical form for ease of handling
            verdict = 'Verdict'
            for vtype in ['TP','OT','FP','FN']: 
                resdb.addField(vtype,evalue=0)
            for entry in resdb.entries():
                entry[entry[verdict]] = 1
            resdb.dropField('Verdict')

            # Want to compress on cloud first
            resdb.makeField('#Method#|#Dataset#|#Cloud#','MCloud')
            resdb.index('MCloud')
            resdb.dropField('Cloud')
            resdb.addField('CloudNum',evalue=1)
            resdb.compress(['MCloud'],{'SeqNum':'mean'},default='sum')  # Add alternatives
            resdb.dropFields(['Rank','Pattern','Cloud'])

            # For now, use each pattern in cloud and normalise by cloudsize
            for entry in resdb.entries():
                for vtype in ['TP','OT','FP','FN']: 
                    entry[vtype] /= entry['CloudNum']
            resdb.dropField('CloudNum')

            # Compress on Dataset/Method
            resdb.makeField('#Method#|#Dataset#','QRun')
            resdb.index('QRun')
            resdb.compress(['QRun'],default='mean')
            resdb.dropField('MCloud')

            # Calculate SP, SN, OTT and OTF for each run
            resdb.addField('SN'); resdb.addField('SP'); resdb.addField('OTT'); resdb.addField('OTF')
            for entry in resdb.entries():
                entry['SN'] = 100.0 * float(entry['TP']) #/ max(int(entry['SeqNum']),(entry['TP'] + entry['FN']))
                if entry['TP'] or entry['FP']: entry['SP'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'])
                else: entry['SP'] = ''
                if entry['TP'] or entry['FP'] or entry['OT']:
                    entry['OTT'] = 100.0 * float(entry['TP']+entry['OT']) / (entry['TP'] + entry['FP'] + entry['OT'])
                    entry['OTF'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'] + entry['OT'])
                else: entry['OTT'] = ''; entry['OTF'] = ''           
            resdb.saveToFile()
            
            resdb.info['Name'] = 'roc_elmhub-%s' % sigcut
            resdb.makeField('#Method#|#ELM#|#Hub#','ERun')
            resdb.index('ERun')
            resdb.compress(['ERun'],default='mean')
            resdb.dropFields(['ELM','Dataset','QRun','Query','SeqNum'])
            resdb.saveToFile()

            resdb.info['Name'] = 'roc_mean-%s' % sigcut
            resdb.compress(['Method'],default='mean')
            resdb.dropFields(['ERun'])
            resdb.addField('Sig',after='Method',evalue=sigcut)
            resdb.saveToFile()

            return

        except:
            self.errorLog('Error in rocData(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def rocData(self,sigcut=0.05,basefile='test'):   ### Generate Specificity and Sensitivity scores for ROC calculations
        '''
        Generate Specificity and Sensitivity scores for ROC calculations.
        >> filename:str = Name of full data file
        >> sigcut [0.05] = Significance cut-off
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Setup rje_db.Database object
            db = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % basefile])

            if os.path.exists('%s.full.tdt' % basefile):
                resdb = db.addTable('%s.full.tdt' % basefile,mainkeys=['Method','Dataset','Rank','Pattern'],name='full')
            else:
                # Load results data into 'full' table
                uniquekey = ['dataset','runid','rank','pattern']    #Replace!
                resdb = db.addTable('All_SLiM.tdt',mainkeys=uniquekey,name='full')
                resdb.addField('Program',evalue='SF')
                resdb.newKey(['dataset','runid','rank','pattern','Program']) # Make new key for table  
                qdb = db.addTable('All_3DID_QSLiM.tdt',mainkeys=uniquekey,name='qsf')
                #qdb = db.addTable('All_QSLiM.tdt',mainkeys=uniquekey,name='qsf')
                qdb.addField('Program',evalue='QSF')
                qdb.dropField('id')
                qdb.newKey(['dataset','runid','rank','pattern','Program']) # Make new key for table  
                db.mergeTables(resdb,qdb,overwrite=True,matchfields=True)

                resdb.dropFields(['masking','elm'])
                resdb.renameField('runid', 'Masking')
                resdb.renameField('pos', 'Verdict')
                resdb.renameField('seqnum', 'SeqNum')
                resdb.renameField('rank', 'Rank')
                resdb.renameField('pattern', 'Pattern')
                resdb.renameField('sig', 'Sig')
                resdb.renameField('cloud', 'Cloud')
                resdb.addField('ELM'); resdb.addField('Query'); resdb.addField('Flank')
                for entry in resdb.entries():
                    [entry['ELM'],entry['Query'],entry['Flank']] = string.split(entry['dataset'],'.')
                    entry['Flank'] = string.split(entry['Flank'],'_')[-1]

                # Reduce to useful data
                # See also code at SBSBINF/Tools/svn/docs/dev_pydocs.html#libraries-rje_db
                # Fields that we want - Use resdb.renameField(oldname, newname) to rename if needed
                resdb.makeField('#Masking#.#Flank#.#Program#','Method')
                resdb.makeField('#ELM#|#Query#','Dataset')
                resdb.newKey(['Method','Dataset','Rank','Pattern']) # Make new key for table  
                useful_fields = ['Dataset','Method','ELM','Query','SeqNum',       # Dataset fields
                                 'Rank','Pattern','Sig','Cloud',                        # Results fields
                                 'Verdict'                                              # TP/FP/OT/FN
                                 ]
                resdb.dropFields(useful_fields,inverse=True)    # Keep only above fields
                resdb.list['Fields'] = ['Method','Dataset','ELM','Query','SeqNum','Rank','Sig','Pattern','Cloud','Verdict']
                resdb.saveToFile()
                for field in resdb.list['Fields'][:2]: resdb.index(field)

            #Define the exact unique method
            # Method = Masking.flank.program (i.e Disorder.whole.QSLiMFinder)
            resdb.dropFields(['Masking','Flank','Program'])     # These only need to split back up for output (maybe)

            # Reformat data (i.e. numbers as numbers)
            resdb.dataFormat({'Sig':'num','SeqNum':'int','Rank':'int'})

            resdb.info['Name'] = 'roc_qry-%s' % sigcut

            # Add Significance filtering
            for entry in resdb.entries()[0:]:
                if entry['Sig'] > sigcut:
                    if entry['Rank'] > 1: resdb.dropEntry(entry)
                    else: entry['Pattern'] = '-'; entry['Verdict'] = 'FN'; entry['Rank'] = 0; entry['Cloud'] = ''
            resdb.dropField('Sig')

            # Convert ratings into numerical form for ease of handling
            verdict = 'Verdict'
            for vtype in ['TP','OT','FP','FN']: 
                resdb.addField(vtype,evalue=0)
            for entry in resdb.entries():
                entry[entry[verdict]] = 1
            resdb.dropField('Verdict')

            # Want to compress on cloud first
            resdb.makeField('#Method#|#Dataset#|#Cloud#','MCloud')
            resdb.index('MCloud')
            resdb.dropField('Cloud')
            resdb.addField('CloudNum',evalue=1)
            resdb.compress(['MCloud'],{'SeqNum':'mean'},default='sum')  # Add alternatives
            resdb.dropFields(['Rank','Pattern','Cloud'])

            # For now, use each pattern in cloud and normalise by cloudsize
            for entry in resdb.entries():
                for vtype in ['TP','OT','FP','FN']: 
                    entry[vtype] /= entry['CloudNum']
            resdb.dropField('CloudNum')

            # Compress on Dataset/Method
            resdb.makeField('#Method#|#Dataset#','QRun')
            resdb.index('QRun')
            resdb.compress(['QRun'],default='mean')
            resdb.dropField('MCloud')

            # Calculate SP, SN, OTT and OTF for each run
            resdb.addField('SN'); resdb.addField('SP'); resdb.addField('OTT'); resdb.addField('OTF')
            for entry in resdb.entries():
                entry['SN'] = 100.0 * float(entry['TP']) #/ max(int(entry['SeqNum']),(entry['TP'] + entry['FN']))
                if entry['TP'] or entry['FP']: entry['SP'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'])
                else: entry['SP'] = ''
                if entry['TP'] or entry['FP'] or entry['OT']:
                    entry['OTT'] = 100.0 * float(entry['TP']+entry['OT']) / (entry['TP'] + entry['FP'] + entry['OT'])
                    entry['OTF'] = 100.0 * float(entry['TP']) / (entry['TP'] + entry['FP'] + entry['OT'])
                else: entry['OTT'] = ''; entry['OTF'] = ''           
            resdb.saveToFile()
            
            resdb.info['Name'] = 'roc_elm-%s' % sigcut
            resdb.makeField('#Method#|#ELM#','ERun')
            resdb.index('ERun')
            resdb.compress(['ERun'],default='mean')
            resdb.dropFields(['ELM','Dataset','QRun','Query','SeqNum'])
            resdb.saveToFile()

            resdb.info['Name'] = 'roc_mean-%s' % sigcut
            resdb.compress(['Method'],default='mean')
            resdb.dropFields(['ERun'])
            resdb.addField('Sig',after='Method',evalue=sigcut)
            resdb.saveToFile()

            return

        except:
            self.errorLog('Error in rocData(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def patFas(self):   ### Reformat Patrick sequence data
        '''Reformat Patrick sequence data.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            db = db.addTable('patrick.tdt',name='patrick',mainkeys=['Number'])
            self.deBug(db.fields())
            FAS = open('patrick.fas','w')
            for e in db.entries():
                FAS.write('>%s %s %s (%s)\n%s\n' % (e['mRNA Name'],e['Genbank Identifier'],e['Annotation'],e['Organism'],e['mRNA Sequence']))
            FAS.close()
        except:
            self.log.errorLog('Error in rje_misc.patFas(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def dtFas(self):    ### Reformat downloaded Algal sequence data
        '''Reformat downloaded Algal sequence data.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#ZEN',rje_zen.Zen().wisdom())
            db = rje_db.Database(self.log,self.cmd_list)
            descdb = db.addTable('DT_sequence_table.tdt',name='SeqDesc',mainkeys=['Seq. Name'])
            seq_cmd_list = self.cmd_list + ['seqin=DT_SeqClean454Isotigs_singletons.fna','seqmode=file','autoload=T']
            seqlist = rje_seqlist.SeqList(self.log,seq_cmd_list)
            ### ~ [1] ~ Make new sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ANN = open('DT.annotated.fas','w'); UNK = open('DT.unknown.fas','w'); ax = 0; ux = 0
            MISS = open('DT.missing.fas','w'); mx = 0
            while seqlist.nextSeq():
                self.progLog('\r#FAS','Making fasta files: %s' % seqlist.progress())
                cseq = seqlist.currSeq()
                cdata = string.split(cseq[0])
                try:
                    desc = descdb.data(cdata[0])['Seq. Description']
                    if desc == '---NA---': UNK.write('>%s %s\n%s\n' % (cdata[0],string.join(cdata[1:]),cseq[1])); ux+=1
                    else: ANN.write('>%s %s %s\n%s\n' % (cdata[0],desc,string.join(cdata[1:]),cseq[1])); ax+=1
                except: MISS.write('>%s %s\n%s\n' % (cdata[0],string.join(cdata[1:]),cseq[1])); mx+=1
            ANN.close(); UNK.close(); MISS.close()
            self.printLog('#FAS','Fasta files generated: %s annotated; %s unknown; %s missing.' % (rje.iStr(ax),rje.iStr(ux),rje.iStr(mx)))
            self.printLog('#ZEN',rje_zen.Zen().wisdom())
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def biol2018(self): ### Identify suitable domains for BIOL2018 projects
        '''Identify suitable domains for BIOL2018 projects.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            pdb = db.addTable('ens_HUMAN.unifake.pfam.tdt','all')
            pdb.index('Type'); pdb.index('Name')
            ## ~ [0a] ~ Set parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            minlen = 75
            minocc = 5
            maxocc = 20            
            ## ~ [0b] ~ Setup dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            data = {}
            dom2gene = {}
            gene2dom = {}
            ### ~ [1] ~ Parse data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in pdb.index('Name'):
                structure = pdb.indexDataList('Name',seq,'Type',sortunique=False)
                structure.sort()
                if len(structure) < 2: continue
                structure = string.join(structure,';')
                for entry in pdb.indexEntries('Name',seq):
                    dom = entry['Type']; seq = entry['Name']; x = int(entry['Start']); y = int(entry['End'])
                    dlen = y - x + 1
                    if dlen < minlen: continue
                    if len(pdb.indexDataList('Type',dom,'Name')) > maxocc: continue
                    if dom not in data: data[dom] = {}
            for seq in pdb.index('Name'):
                structure = []
                for dom in pdb.indexDataList('Name',seq,'Type',sortunique=False):
                    if dom in data: structure.append(dom)   #'%s (%d)' % (dom,len(pdb.indexDataList('Type',dom,'Name'))))
                structure.sort()
                if len(structure) < 2: continue
                structure = string.join(structure,'; ')
                for entry in pdb.indexEntries('Name',seq):
                    dom = entry['Type']; seq = entry['Name']; x = int(entry['Start']); y = int(entry['End'])
                    if dom not in data: continue
                    if structure not in data[dom]: data[dom][structure] = []
                    if seq not in data[dom][structure]: data[dom][structure].append(seq)
            reduction = True
            while reduction:
                reduction = False
                for dom in rje.sortKeys(data):
                    self.progLog('\r#DOM','%s dom  ' % rje.iStr(len(data)))
                    for structure in rje.sortKeys(data[dom]):
                        domlist = string.split(structure,'; ')
                        newstruct = []
                        for sdom in domlist:
                            if sdom in data and len(data[sdom]) > 1: newstruct.append(sdom)
                            else: reduction = True
                        if not newstruct: data[dom].pop(structure); continue
                        shared = {}
                        for sdom in newstruct[0:]:
                            if sdom == dom: continue
                            shared = 0
                            for sharestruct in data[dom]:
                                if sdom in string.split(sharestruct,'; '): shared += 1
                            if shared == len(data[dom]): newstruct = []; reduction = True; break    # Redundant domain
                        if len(newstruct) < 2: data[dom].pop(structure); reduction = True; continue
                        newstruct = string.join(newstruct,'; ')
                        data[dom][newstruct] = data[dom].pop(structure)
                    if not data[dom]: data.pop(dom) #?#
            for dom in data:
                for structure in rje.sortKeys(data[dom]):
                    domlist = string.split(structure,'; ')
                    newstruct = []
                    for sdom in domlist: newstruct.append('%s (%d/%d)' % (sdom,len(pdb.indexDataList('Type',sdom,'Name')),len(data[sdom])))
                    newstruct = string.join(newstruct,'; ')
                    data[dom][newstruct] = data[dom].pop(structure)
            BIOL = open('biol2018.tdt','w')
            BIOL.write('Domain\tStructure\tSeqN\tSeq\n')
            for dom in rje.sortKeys(data):
                if len(data[dom]) < 2: continue
                for structure in rje.sortKeys(data[dom]):
                    data[dom][structure].sort()
                    for seq in data[dom][structure]: BIOL.write('%s\t%s\t%d\t%s\n' % (dom,structure,len(data[dom][structure]),seq))
            BIOL.close()
        except: self.errorLog('Error in rje_misc.run(%s)' % self.info['Name'],printerror=True,quitchoice=True)
#########################################################################################################################
    def biol3050(self): ### Compare clones to ENST cDNAs
        '''Compare clones to ENST cDNAs.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            maindb = db.addEmptyTable('main',['clone','enst','hgnc','strand','cand'],keys=['clone','enst'])
            est = rje_seq.SeqList(self.log,self.cmd_list+['seqin=clones_est.fas'])
            nuc = rje_seq.SeqList(self.log,self.cmd_list+['seqin=clones_nuc.fas'])
            hum = rje_seq.SeqList(self.log,self.cmd_list+['seqin=Homo_sapiens.GRCh37.59.cdna.fas'])
            candidates = string.split(open('biol3050.more_genes.txt').read(),'\n')
            seqdict = {}; enst2gene = {}; sx = 0.0; stot = hum.seqNum()
            for seq in hum.seqs():
                self.progLog('\r#SEQ','Making sequence dictionary: %.2f%%' % (sx/stot)); sx += 100.0
                seq.trimPolyA()
                details = string.split(seq.info['Name'],'|')
                cdna = seq.info['Sequence']
                enst = details[0]; hgnc = details[-1]
                if cdna not in seqdict: seqdict[cdna] = []
                seqdict[cdna].append(enst); enst2gene[enst] = hgnc
            self.printLog('\r#SEQ','Making sequence dictionary: %.2f%%' % (sx/stot))
            fullens = string.join(seqdict.keys())

            cloneseq = {}; sx = 0.0; stot = est.seqNum() + nuc.seqNum()
            for clone in est.seqs() + nuc.seqs():
                self.progLog('\r#SEQ','Making sequence dictionary 2: %.2f%%' % (sx/stot)); sx += 100.0
                clone.trimPolyA()
                for cseq in [clone.info['Sequence'], rje_sequence.reverseComplement(clone.info['Sequence'])]:
                    if cseq not in cloneseq: cloneseq[cseq] = []
                    cloneseq[cseq].append(clone)
            self.printLog('\r#SEQ','Making sequence dictionary 2: %.2f%%' % (sx/stot))
            fullclone = string.join(cloneseq.keys())

            ### ~ [1] ~ Match clones to transcripts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #sx = 0.0; stot = len(cloneseq) + len(seqdict); cpairs = []
            #for cseq in cloneseq:
            #    self.progLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (len(cpairs),sx/stot)); sx += 100.0
            #    ix = 0
            #    eseq = rje.matchExp('(\S*%s\S*)' % cseq,fullens[ix:])
            #    self.deBug(eseq)
            #    while eseq:
           #         eseq = eseq[0]; ix = fullens.find(eseq) + len(eseq)
            #        if (cseq,eseq) not in cpairs: cpairs.append((cseq,eseq))
            #        eseq = rje.matchExp('(\S*%s\S*)' % cseq,fullens[ix:])
            #for eseq in seqdict:
            #    self.progLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (len(cpairs),sx/stot)); sx += 100.0
            #    ix = 0
            #    cseq = rje.matchExp('(\S*%s\S*)' % eseq,fullclone[ix:])
            #    self.deBug(cseq)
            #    while cseq:
            #        cseq = cseq[0]; ix = fullclone.find(cseq) + len(cseq)
            #        if (cseq,eseq) not in cpairs: cpairs.append((cseq,eseq))
            #        cseq = rje.matchExp('(\S*%s\S*)' % eseq,fullclone[ix:])
            #self.printLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (len(cpairs),sx/stot)); sx += 100.0
            sx = 0.0; stot = len(seqdict); cpairs = []
            for eseq in seqdict:
                self.progLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (len(cpairs),sx/stot)); sx += 100.0
                for cseq in cloneseq:
                    if (cseq,eseq) in cpairs: continue
                    if eseq in cseq or cseq in eseq: cpairs.append((cseq,eseq))
            self.printLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (len(cpairs),sx/stot))
            cx = 0
            for (cseq,eseq) in cpairs:
                for clone in cloneseq[cseq]:
                    for enst in seqdict[eseq]:
                        entry = {'clone':string.split(clone.shortName(),'|')[3],'enst':enst,'hgnc':enst2gene[enst],
                                 'strand':dirn,'cand':enst2gene[enst] in candidates}
                        dkey = maindb.makeKey(entry)
                        maindb.data()[dkey] = entry
                        cx += 1
            self.printLog('\r#MATCH','Matching %d clones to human cDNAs: %.2f%%' % (cx,sx/stot))
            maindb.saveToFile('clones_enst.tdt')
        except: self.errorLog('Error in rje_misc.biol3050')            
#########################################################################################################################
    def processXML(self,myxml):     ### Processes an XML object from rje_xml
        '''Processes an XML object from rje_xml.'''
        try:### ~ [1] FlyBase ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'] == 'flyseq':
                (organism,name,sequence,acc) = ('UNK','p','X','NA')
                # Organism
                try: organism = myxml.getXML(['organism_id','organism','abbreviation']).info['Content']
                except: self.log.errorLog(rje_zen.Zen().wisdom())
                # Gene
                try: name = myxml.getXML(['name']).info['Content']
                except: self.log.errorLog(rje_zen.Zen().wisdom())
                # Sequence
                try: sequence = myxml.getXML(['residues']).info['Content']
                except: self.log.errorLog(rje_zen.Zen().wisdom())
                # Acc
                #try: acc = myxml.getXML(['feature_synonym','synonym_id','name']).info['Content']
                #except: self.log.errorLog(rje_zen.Zen().wisdom())
                self.obj['SeqList']._addSeq('%s_%s__%s' % (name,organism,acc),sequence)
                open('tmp','a').write('>%s %s %s\n%s\n' % (name,organism,acc,sequence))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def benCOG(self):   ### Makes Ben COG files
        '''Makes Ben COG files.'''
        try:### ~ [0] DATIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            unidict = {}
            ex = 0; ix = 0; efile = 'ethanol_spec.dat'
            DAT = open(efile,'r')
            while efile:
                line = rje.chomp(DAT.readline())
                if line[:2] == '': break
                if line[:2] not in ['AC','ID']: continue
                if line[:5] == 'ID   ':
                    primacc = ''
                    id = string.split(line)[1]
                    [gene,spec] = string.split(id,'_')
                    ex += 1
                if line[:5] == 'AC   ':
                    acclist = string.split(string.replace(line,';',''))[1:]
                    if 'AC' in acclist: acclist.remove('AC')
                    if not primacc: primacc = acclist[0]
                    for acc in acclist: unidict[acc] = {'ID':acc,'AccNum':primacc,'Spec':spec}; ix += 1
                    if gene != primacc: unidict[id] = {'ID':id,'AccNum':primacc,'Spec':spec}; ix += 1
                    self.progLog('#ACC','%s acc/id read for %s entries' % (rje.integerString(ix),rje.integerString(ex)))
            self.printLog('\r#ACC','%s acc/id read for %s entries' % (rje.integerString(ix),rje.integerString(ex)))
            DAT.close()
            ### ~ [1] COG File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cog = rje.dataDict(self,'/home/bbi1g08/BenIent/UniProtAC2eggNOG.tsv',['COGAcc'],['COGAcc','COG'],'\t',['COGAcc','COG'])
            self.printLog('#COG','COG families read for %s acc/id' % (rje.integerString(len(cog))))
            ## ~ [1a] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for acc in cog.keys():
                if acc not in unidict: cog.pop(acc)
            self.printLog('#COG','COG families reduced to %s UniProt acc/id' % (rje.integerString(len(cog))))
            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] COG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cogtdt = 'ben_cog.tdt'
            coghead = ['COGAcc','COG','AccNum','Spec']
            rje.delimitedFileOutput(self,cogtdt,coghead)
            for acc in rje.sortKeys(cog):
                cog[acc]['AccNum'] = unidict[acc]['AccNum']
                cog[acc]['Spec'] = unidict[acc]['Spec']
                rje.delimitedFileOutput(self,cogtdt,coghead,datadict=cog[acc])
            self.printLog('#COG','COG data output for %s UniProt acc/id' % (rje.integerString(len(cog))))
            ## ~ [2b] UniMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cogtdt = 'unimap.tdt'
            coghead = ['ID','AccNum','Spec']
            rje.delimitedFileOutput(self,cogtdt,coghead)
            for acc in rje.sortKeys(unidict): rje.delimitedFileOutput(self,cogtdt,coghead,datadict=unidict[acc])
            self.printLog('#UNI','UniProt AccNum mapping output for %s UniProt acc/id' % (rje.integerString(len(unidict))))
            
            
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def taxaDBSum(self):    ### Cleans up E hux TaxaDB results
        '''Cleans up E hux TaxaDB results.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datafile = 'taxadbsum.csv'
            data = {}       # strain:{hit:{???}}
            pep2prot = {}   # strain:{peptide:[hits]}
            id2prot = {}    # strain:{id:hit}
            prot2desc = {}
            fullpeplist = {}    
            pepcon = {}     # Convert pep:longer pep
            speclist = []   # List of species codes
            ### ~ [1] Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            indata = rje.dataDict(self,datafile,['strain','prot_hit_num'],'All',lists=True)
            for ikey in rje.sortKeys(indata):
                (strain,id) = string.split(ikey,',')
                prot = indata[ikey]['prot_acc'][0]
                desc = string.replace(indata[ikey]['prot_desc'][0],'Full=','')
                if desc[3:7] == 'Name': desc = desc[9:]
                prot2desc[prot] = desc; self.printLog('#DESC','%s = %s' % (prot,desc))
                indata[ikey]['pep_seq'] = string.join(indata[ikey]['pep_seq'],'|')
                pepconv = string.replace(indata[ikey]['pep_seq'],'I','L')
                pepconv = string.replace(pepconv,'Q','K')
                peplist = rje.sortUnique(string.split(pepconv,'|'))
                indata[ikey]['pep_seq'] = string.join(rje.sortUnique(string.split(indata[ikey]['pep_seq'],'|')),'|')
                if strain not in data:
                    data[strain] = {}
                    pep2prot[strain] = {}
                    id2prot[strain] = {}
                    fullpeplist[strain] = []
                    pepcon[strain] = {}
                fullpeplist[strain] += peplist
                id2prot[strain][id] = prot
                spec = string.split(prot,'_')[1]
                if spec not in speclist: speclist.append(spec)
                data[strain][prot] = {'strain':strain,'pepcount':len(peplist),'hit':id,'desc':desc,
                                      'accnum':string.split(prot,'_')[-1],'spec':spec,
                                      'pep_uniq':0,'peplist':indata[ikey]['pep_seq'],'conpep':peplist[0:],
                                      'pep_rem':0}
                for pep in peplist:
                    if pep not in pep2prot[strain]:
                        pep2prot[strain][pep] = []
                    pep2prot[strain][pep].append(prot)
            ## ~ [1a] Convert peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for strain in fullpeplist:
                fullpeplist[strain] = rje.sortUnique(fullpeplist[strain])
                for pep in fullpeplist[strain][0:]:
                    for pep2 in fullpeplist[strain]:
                        if pep != pep2 and pep in pep2:
                            pepcon[strain][pep] = pep2
                            fullpeplist[strain].remove(pep)
                            break
                for pep in pepcon[strain]:
                    while pepcon[strain][pep] in pepcon[strain]: pepcon[strain][pep] = pepcon[strain][pepcon[pep]]
                self.printLog('#PEP','%s %s peptide conversions' % (len(pepcon[strain]),strain))
                #self.deBug(pepcon[strain])
                #self.deBug(rje.sortKeys(pep2prot[strain]))
                pp = 0; pm = 0
                for prot in data[strain]:
                    for pep in data[strain][prot]['conpep'][0:]:
                        if pep in pepcon[strain]:
                            newpep = pepcon[strain][pep]
                            if newpep not in data[strain][prot]['conpep']: data[strain][prot]['conpep'].append(newpep); pp += 1
                            data[strain][prot]['conpep'].remove(pep); pm += 0
                            if prot not in pep2prot[strain][newpep]: pep2prot[strain][newpep].append(prot)
                            if pep in pep2prot[strain]: pep2prot[strain].pop(pep)
                    data[strain][prot]['pep_con'] = len(data[strain][prot]['conpep'])
                self.printLog('#PEP','%s %s converted peptides added; %s removed' % (pp,strain,pm))
            ### ~ [2] Calculate Unique/Redundancy status ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for strain in pep2prot:
            ## ~ [2a] Species Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                remx = 0
                for prot in data[strain]:
                    if data[strain][prot]['spec'] != 'EMIHU': continue
                    for pep in data[strain][prot]['conpep']:
                        for prot2 in pep2prot[strain][pep][0:]:
                            if data[strain][prot2]['spec'] == 'EMIHU': continue
                            pep2prot[strain][pep].remove(prot2)
                            data[strain][prot2]['conpep'].remove(pep)
                            data[strain][prot2]['pep_rem'] += 1; remx += 1
                self.printLog('#REM','%s %s peptides removed from non-EMIHU hits' % (rje.integerString(remx),strain))
            ## ~ [2b] One-hit wonders ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for prot in data[strain]:
                    if len(data[strain][prot]['conpep']) < 2:
                        for pep in data[strain][prot]['conpep']:
                            #if pep in pep2prot[strain] and prot in pep2prot[strain][pep]:
                            pep2prot[strain][pep].remove(prot)
            ## ~ [2c] Unique peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ux = 0
                for pep in pep2prot[strain]:
                    #self.deBug(pep)
                    if len(pep2prot[strain][pep]) == 1: data[strain][pep2prot[strain][pep][0]]['pep_uniq'] += 1; ux += 1
                self.printLog('#UNIQ','%s unique %s peptides' % (rje.integerString(ux),strain))
            ## ~ [2d] Total Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                summary = {'HITS':len(data[strain]),'REJECT':0,'UNIQUE':0,'NR':0,'REDUNDANT':0}
                rx = 0
                for prot in data[strain]:
                    #if data[strain][prot]['unique']: data[strain][prot]['red'] = False; continue
                    data[strain][prot]['pep_red'] = 0   # Redundant peptides found in proteins with unique peptides
                    data[strain][prot]['pep_nr'] = 0    # Redundant peptides found only in proteins without unique peptides
                    for pep in data[strain][prot]['conpep']:
                        if pep2prot[strain][pep] == [prot]: continue
                        upep = False
                        for prot2 in pep2prot[strain][pep]:
                            if data[strain][prot2]['pep_uniq']: upep = True; break
                        if upep: data[strain][prot]['pep_nr'] += 1   # Redundant peptide not found in unique protein
                        else: data[strain][prot]['pep_red'] += 1
                    if len(data[strain][prot]['conpep']) < 2: data[strain][prot]['class'] = 'REJECT'; rx += 1
                    elif data[strain][prot]['pep_uniq']: data[strain][prot]['class'] = 'UNIQUE'
                    elif data[strain][prot]['pep_nr']: data[strain][prot]['class'] = 'NR'
                    else: data[strain][prot]['class'] = 'REDUNDANT'; rx += 1
                    summary[data[strain][prot]['class']] += 1
                self.printLog('#REJ','%s rejected %s hits' % (rje.integerString(rx),strain))
                for x in rje.sortKeys(summary): self.printLog('#%s' % strain,'%s %s' % (summary[x],x))

            ### ~ [3] Species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            speclist.sort()
            species = {}
            for spec in speclist:
                grep = os.popen('grep %s /scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/UniProt/uniprot.spec.tdt' % spec).read()
                species[spec] = string.split(grep,':')[-4]
                self.printLog('#SPEC','%s = %s' % (spec,species[spec]))

            ### ~ [END] Output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = 'taxadbsum.clean.tdt'
            headers = ['strain','hit','class','accnum','spec','species','desc','pepcount','pep_con','pep_rem','pep_uniq','pep_nr','pep_red','peplist','conpep']
            rje.delimitedFileOutput(self,outfile,headers,datadict={},rje_backup=True)
            for strain in rje.sortKeys(data):
                for prot in rje.sortKeys(data[strain]):
                    data[strain][prot]['species'] = species[data[strain][prot]['spec']]                                                                               
                    rje.delimitedFileOutput(self,outfile,headers,datadict=data[strain][prot])
                    
                        
                
            
        except: self.errorLog('Errg')
#########################################################################################################################
    def laavanya(self): ### Prepares PTPRD datasets for SLiMFinder
        '''Prepares PTPRD datasets for SLiMFinder.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,'PTPRD_PPI')
            ptprd = os.popen('grep -A 130 PTPRD_HUMAN /rhome/re1u06/Data/ensdat.dat').read()
            ppath = rje.makePath('/rhome/re1u06/Data/PinguScreenDDI/')
            ppi = rje.dataDict(self,'pingu.ppi.tdt',mainkeys=['EnsLoci'],datakeys=['Gene'],lists=True)
            ppi.pop('-')
            accgenes = {}
            for key in rje.sortKeys(ppi):
                try:
                    acc = string.split(key,'__')[1]
                    if acc not in accgenes: accgenes[acc] = []
                    accgenes[acc] += ppi[key]['Gene']
                except: self.errorLog('Error with %s' % key)
            ### ~ [2] ~ Make Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for acc in self.loadFromFile('ptprd_humap.acc',chomplines=True):
                if acc in accgenes:
                    for gene in accgenes[acc]:
                        pacc = self.loadFromFile('%s%s.ppi.acc' % (ppath,gene),chomplines=True)
                        px = len(pacc)
                        if 'P23468' not in pacc: px += 1
                        if px >= 3:
                            self.printLog('#PPI','Acc %s => %s => %d PPI data' % (acc,gene,px))
                            os.system('cp %s%s.ppi.dat PTPRD_PPI/%s.dat' % (ppath,gene,gene))
                            if 'P23468' not in pacc: open('PTPRD_PPI/%s.dat' % gene,'a').write(ptprd)
                        else: self.printLog('#ACC','Acc %s => %s => insufficient PPI data' % (acc,gene))
                else: self.printLog('#ACC','Acc %s has no PPI data' % acc)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def intEnrich(self): ### Converts APHID enrichement file for Pingu
        '''Converts APHID enrichement file for Pingu.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            infile = 'hm_comb.full.filtered.gene.count.enrichment.tdt'
            outfile = 'hm_comb.enriched.tdt'
            headers = ['Sample','Identifier','Enrichment']
            data = rje.dataDict(self,infile,getheaders=True)
            ihead = data.pop('Headers')
            samples = ['Activated','Lysate','Resting']
            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.delimitedFileOutput(self,outfile,headers,rje_backup=True)
            for id in rje.sortKeys(data):
                for s1 in samples:
                    for s2 in samples:
                        if s1 == s2:
                            datadict = {'Sample':s1,'Enrichment':data[id][s1],'Identifier':id}
                            if not datadict['Enrichment'] or float(datadict['Enrichment']) <= 0.5: continue
                        else:
                            datadict = {'Sample':'%s-%s' % (s1[:3],s2[:3]),'Enrichment':data[id]['%s/%s' % (s1,s2)],'Identifier':id}
                            if not datadict['Enrichment'] or float(datadict['Enrichment']) <= 2: continue
                        rje.delimitedFileOutput(self,outfile,headers,datadict=datadict)
            self.printLog('#ZEN',rje_zen.Zen().wisdom())
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def arath2GO(self): ### Converts ARATH sequence to GO (if possible)
        '''Converts ARATH sequence to GO (if possible).'''
        try:
            accfile = '/scratch/Projects-Proteomics/2008-01-Emiliania_huxleyi/HeadercorrectedF003975/HeadercorrectedF003975.ARATH.mapping.acc'
            acclist = self.loadFromFile(accfile,chomplines=True)
            mapfile = '/scratch/Databases/NewDB/Arabidopsis/UniProt2AGI.20070820'
            atmap = {}
            for line in self.loadFromFile(mapfile,chomplines=True):
                if rje.matchExp('^(\S+)\s+(AT\S+)',line):   # Mapping lines
                    (acc,at) = rje.matchExp('^(\S+)\s+(AT\S+)',line)
                    if acc in acclist: atmap[acc] = string.split(at,':')
            gofile = '/scratch/Databases/NewDB/Arabidopsis/ATH_GO_GOSLIM.20080216.txt'
            #Column headers :explanation
            #1. locus name: standard AGI convention name
            #2. TAIR accession:the unique identifier for an object in the TAIR database-the object type is the prefix, followed by a unique accession number(e.g. gene:12345).
            #3.object name : the name of the object (gene, protein, locus) being annotated.#
            #4. GO term: the actual string of letters corresponding to the GO ID
            #5. GO ID: the unique identifier for a GO term.
            #6. TAIR Keyword ID: the unique identifier for a keyword in the TAIR database.
            #7.  Aspect: F=molecular function, C=cellular component, P=biological process.
            #8. GOslim term: high level GO term helps in functional categorization.
            #9. Evidence code: three letter code for evidence types (see: http://www.geneontology.org/GO.evidence.html).
            #10. Reference: Either a TAIR accession for a reference (reference table: reference_id) or reference from PubMed (e.g. PMID:1234).
            #11. Annotating database: TAIR or TIGR
            #12. Date annotated: date the annotation was made.
            gomap = {}
            godata = {}
            gosdat = {}
            for line in self.loadFromFile(gofile,chomplines=True):
                data = string.split(line,'\t')
                if len(data) < 12: continue
                gene = data[0]
                if gene not in gomap: gomap[gene] = []
                go = data[3]
                goid = data[4]
                gotype = data[6]
                goslim = data[7]
                godata[go] = {'GO':go,'ID':goid,'Type':gotype,'Slim':goslim,'Acc':[]}
                gosdat[goslim] = {'Type':gotype,'Acc':[],'ID':'.','GO':'.','Slim':goslim}
                if go not in gomap[gene]: gomap[gene].append(go)
            acclist.sort()
            for acc in acclist:
                gomap[acc] = []
                try: genes = atmap[acc]
                except:
                    self.log.errorLog('Acc %s not mapped!' % acc)
                    continue
                for map in genes:
                    gene = string.split(map,'.')[0]
                    try: gomap[acc] += gomap[gene][0:]
                    except: self.log.errorLog('Problem with gomap for %s' % gene)
                gomap[acc] = rje.sortUnique(gomap[acc],False)
                for go in gomap[acc]:
                    godata[go]['Acc'].append(acc)
                    gosdat[godata[go]['Slim']]['Acc'].append(acc)  
            gfile = 'HeadercorrectedF003975.ARATH.gomap.tdt'
            ghead = ['Slim','Type','GO','ID','N','Acc']
            rje.delimitedFileOutput(self,gfile,ghead,'\t',rje_backup=True)
            for goslim in rje.sortKeys(gosdat):
                if not gosdat[goslim]['Acc']: continue
                gosdat[goslim]['N'] = len(gosdat[goslim]['Acc'])
                gosdat[goslim]['Acc'] = string.join(gosdat[goslim]['Acc'],',')
                rje.delimitedFileOutput(self,gfile,ghead,'\t',gosdat[goslim])
            for go in rje.sortKeys(godata):
                if not godata[go]['Acc']: continue
                godata[go]['N'] = len(godata[go]['Acc'])
                godata[go]['Acc'] = string.join(godata[go]['Acc'],',')
                rje.delimitedFileOutput(self,gfile,ghead,'\t',godata[go])
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def flySeq(self):   ### Extract sequences from FlyBase chado_xml
        '''Extract sequences from FlyBase chado_xml.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flyxml = rje_xml.XML(self.log,self.cmd_list)
            flyxml.obj['Processor'] = self
            self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F','autofilter=F'])
            open('tmp','w')
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flyxml.parseXML('chado_proteins.xml')
            self.obj['SeqList'].saveFasta(seqfile='chado.fas')
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def slimCons(self):     ### Manipulate and output slimcons results data for R analysis
        '''Manipulate and output slimcons results data for R analysis.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            projectdir = "D:\\Projects - Motifs\\2007-12 EnsDatHuman SLiMFinder"
            os.chdir(projectdir)		    # Set working directory
            rdb = rje_db.Database(self.log,self.cmd_list)
            outfile = 'slimcons_analysis.txt'

            ### ~ [1] ~ How many "significant" motifs are there? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['RunID','p05-Total','p05-Dataset','p05-Cloud','FDR05-Total','FDR05-Dataset','FDR05-Cloud']
            sigdata = {}
            sf = rdb.addTable('slimcons.csv',['Dataset','RunID','Pattern'],delimit=',')
            sf.dropFields(['Masking','Build','RunTime'])
            dx = len(sf.index('Dataset'))
            self.log.printLog('#DATASET','Results returned for %s Datasets' % rje.integerString(dx))
            for runid in rje.sortKeys(sf.index('RunID')): sigdata[runid] = {'RunID':runid}
            ## ~ [1a] Total Counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            p = 0.05
            sf.info['Name'] = 'Total'
            sf.dropEntries(['Rank<1','Sig>0.05'])
            rindex = sf.index('RunID')
            for runid in rindex: sigdata[runid]['p05-Total'] = len(rindex[runid])
            ## ~ [1b] Dataset Counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ds = rdb.copyTable(sf,'Dataset')
            ds.dropEntries(['Rank>1'])
            dindex = ds.index('RunID')
            for runid in dindex: sigdata[runid]['p05-Dataset'] = len(dindex[runid])
            ## ~ [1c] Cloud Counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cl = rdb.copyTable(sf,'Cloud')
            cl.compress(['Dataset','RunID','Cloud'],rules={'Sig':'min','Rank':'min'},best=['Sig'])
            cindex = cl.index('RunID')
            for runid in cindex: sigdata[runid]['p05-Cloud'] = len(cindex[runid])
            cl.saveToFile(backup=False)
            ## ~ [1d] FDR calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tables = {'Total':sf,'Dataset':ds,'Cloud':cl}
            for level in tables:
                runsplits = rdb.splitTable(tables[level],'RunID')
                for runid in sigdata:
                    focus = rdb.getTable('%s_%s' % (level,runid))
                    siglist = []
                    for entry in focus.entries():
                        sig = float(entry['Sig'])
                        if sig not in siglist: siglist.append(sig)
                    siglist.sort()
                    siglist.reverse()
                    for sig in siglist:
                        expected = dx * sig
                        focus.dropEntries(['Sig>%s' % sig],log=False)
                        fdr = expected / focus.entryNum()
                        if fdr <= p: break
                    sigdata[runid]['FDR05-%s' % level] = '%s (p=%s)' % (focus.entryNum(),sig)
            ## ~ [1e] Output summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            OUT = open(outfile,'w').write('Summary of "significant" motif numbers from different runs\n\n')
            rje.delimitedFileOutput(self,outfile,headers,delimit='\t')
            for runid in sigdata: rje.delimitedFileOutput(self,outfile,headers,delimit='\t',datadict=sigdata[runid])
                    
            ### ~ [2] ~ What is their relationship to known motifs? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cm = rdb.addTable('slimcons.compare.tdt',['MotifFile','SearchDB','Name1','Name2'])
            t1 = rdb.joinTables(name='AllData',join=[(sf,'#Dataset###RunID###Rank###Pattern#'),(cm,'Name1')],newkey=['Dataset','RunID','Pattern','Name2'],cleanup=True,delimit='\t',empties=True)
            t1.dropFields(['#Dataset###RunID###Rank###Pattern#','MotifFile','SearchDB','Name1','Motif1','Desc1'])
            t1.saveToFile(backup=False)
            ## ~ [2a] ~ Reduce to best match per motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            t1.info['Name'] = 'BestMatch'
            #x#t1.opt['DeBug'] = True
            t1.compress(['Dataset','RunID','Pattern'],best=['Score','NormIC','MatchIC','MatchPos'])
            t1.saveToFile(backup=False)
            t1.dropEntries(['MatchPos<3'])      # Actually reduces to just motifs with a match
            
        

            return
                                                                   








            ## ~ [1a] Temp replacement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ##for dkey in cm.data(): cm.data()[dkey]['Name1'] = string.join(string.split(cm.data()[dkey]['Name1'],'#')[:3] + [cm.data()[dkey]['Motif1']],'#')
            ##cm.saveToFile(backup=True)   # Temp
            
            ### ~ [2] Join data from Dataset summaries and motif comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #t1.index('Name1')
            ## ~ [2a] Compress to a single entry per cloud ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            t1.compress(['Dataset','RunID','Cloud'],rules={'MatchIC':'max'},default='text')
            t1.saveToFile(backup=False)

            return        
            #t1.saveToFile(filename='newtable.tdt',backup=False)
            #slist = red.splitTable(t1,'RunID')
            slist = red.splitTable(sf,'RunID')
            for s in slist:
                j = red.joinTables(name='Data',join=[(s,'#Dataset###RunID###Rank###Pattern#'),(cm,'Name1')],newkey=['Dataset','RunID','Pattern','SearchDB','Name2'],cleanup=True,delimit='\t',empties=True)
                j.index('Name1')
                try: red.list['Tables'].remove(s)
                except:
                    self.log.errorLog(rje_zen.Zen().wisdom())
                    print red.list
                j.info['Name'] = string.join(string.split(s.name(),'_')[:1],'_')
                j.saveToFile(backup=False)
                #s.saveToFile(backup=False)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def phoMo(self):    ### PhosphoMotif Finder reformatting
        '''
        PhosphoMotif Finder reformatting. The HTML source for the pages should be stripped down to just the table row
        entries with dividers establish ST vs Y and kinase vs binding.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            phomo = self.loadFromFile("D:\\Databases\\PhosphoMotif Finder\\phosphomotif_html.txt",chomplines=True)
            OUT = open("D:\\Databases\\PhosphoMotif Finder\\phosphomotif.motifs",'w')
            type = ['','']      # AA, kinase/binding
            x = {}              # Counter

            ### ~ [2] ~ Read and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in phomo:
                ## ~ [2a] Type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not rje.matchExp('(\S)',line): continue      # Blank
                elif line[:1] == 'Y': type = ['KIN','Y']
                elif line[:2] == 'ST': type = ['KIN','ST']
                elif line[:4] == 'Bind': type[0] = 'BIND'
                ptype = string.join(type,'_')
                if ptype not in x: x[ptype] = 1
                ## ~ [2b] Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line[:4] == '<tr>':      # New entry
                    line = string.replace(line,'span style="color: red;">','')
                    line = string.replace(line,'</span>','')
                    this = string.split(line,'</td>')[1:]
                    pattern = string.split(string.join(string.split(this[0],'<'),''),'>')[1]
                    pattern = string.replace(pattern,'p','')
                    pattern = string.replace(pattern,'*','')
                    pattern = string.replace(pattern,'/','')
                    pattern = string.replace(pattern,'X','.')
                    desc = string.split(string.join(string.split(this[1],'<'),''),'>')[1]
                    if desc[-1] != ' ': desc += ' '
                    desc += '[PMF] [PMID:%s]' % rje.matchExp('list_uids=(\S+)"',string.split(string.join(string.split(this[2],'<'),''),'>')[1])[0]
                    name = '%s%s' % (ptype,rje.preZero(x[ptype],len(phomo)))
                    x[ptype] += 1
                    OUT.write('%s\t%s\t# %s\n' % (name, pattern, desc))
                    print '%s\t%s\t# %s' % (name, pattern, desc)
            OUT.close()
            self.log.printLog('#OUT','Output complete!')
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def ensDat(self):   ### EnsDat cleanup following rje_ensembl boo-boo
        '''EnsDat cleanup following rje_ensembl boo-boo.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN = open('ens_HUMAN.ensdat.dat','r')
            OUT = open('new_ensdat.dat','w')
            dx = 0      # No. of proteins processed
            hgnc = rje.dataDict(self,'../../HGNC/HPRD.genecards.tdt',['EnsEMBL'],['Symbol','Alias'],lists=True)

            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            line = '//'
            while line == '//':
                ## ~ [2a] Read in next entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = []
                while line:
                    line = rje.chomp(IN.readline())
                    entry.append(line)
                    if line == '//': break
                ## ~ [2b] Correct entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Parse all acc numbers and add to AC
                if entry[1][:2] != 'AC': print string.join(['']+entry,'\n')
                acclist = [rje.matchExp('^AC\s+(\S+);',entry[1])[0]]
                accmap = rje.matchExp('^DE.+\[acc:(\S+) pep:(\S+) gene:(\S+)\]',entry[3])
                if not accmap:
                    print string.join(['']+entry,'\n')
                    print accmap
                for acc in accmap:
                    if acc not in acclist + ['-']: acclist.append(acc)
                entry[1] = 'AC   %s' % string.join(acclist,'; ')
                # Cross-reference GeneCards Gene and add it and EnsG to GN
                ens = accmap[-1]
                genelist = []
                if ens in hgnc:
                    for gene in hgnc[ens]['Symbol'] + [ens] + hgnc[ens]['Alias']:
                        if gene not in genelist: genelist.append(gene)
                else: genelist = [ens]
                entry.insert(4,'GN   %s' % string.join(genelist,'; '))
                # Modify the DOMAIN/PFAM assignment for new iucut.
                for i in range(len(entry)):
                    ft = entry[i]
                    if ft[:2] != 'FT': continue
                    ftx = string.split(ft)
                    if ftx[1] not in ['PFAM']: continue
                    if float(ftx[-1]) < 0.5: entry[i] = string.replace(entry[i],'PFAM','DOMAIN')
                # Re-output entry #
                OUT.write(string.join(entry,'\n'))
                dx += 1
                self.log.printLog('\r#DAT','Reformatted %s entries' % rje.integerString(dx),newline=False,log=False)

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
            self.log.printLog('\r#DAT','Reformatted %s entries' % rje.integerString(dx))
            IN.close()
            OUT.close()
                
        except: self.log.errorLog('EnsDat error')            
#########################################################################################################################
    def slimfinderVsTeiresias(self):    ### Runs SLiMBuild and TEIRESIAS on files and compares run times
        '''Runs SLiMBuild and TEIRESIAS on files and compares run times.'''
        try:### ~ Setup sequence files to search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ### ~ These files are based on the HPRD interaction dataset for AAA (AAA_hprd.upc)
            makefiles = False
            abcfam = string.split('ABCD2_HUMAN__Q9UBJ2 ABCD1_HUMAN__P33897 ABCD3_HUMAN__P28288 TAP1_HUMAN__NP_000584.2 TAP2_HUMAN__Q03519 CFTR_HUMAN__P13569 ABCB8_HUMAN__Q9NUT2 ABCG8_HUMAN__Q9H221 ABCG1_HUMAN__P45844 ABCG4_HUMAN__Q9H172 ABCG5_HUMAN__Q9H222')
            strfam = string.split('ESR1_HUMAN__P03372 RXRA_HUMAN__P19793 THRB_HUMAN__P10828 RORA_HUMAN__P35398 PPARG_HUMAN__P37231 PPARA_HUMAN__Q07869')
            mcmfam = string.split('MCM2_HUMAN__P49736 MCM4_HUMAN__P33991 MCM5_HUMAN__P33992 MCM6_HUMAN__Q14566 MCM3_HUMAN__P25205 MCM7_HUMAN__P33993')
            singles = string.split('AKAP8L_HUMAN__Q9ULX6 SSNA1_HUMAN__O43805 TBC1D4_HUMAN__O60343 ATP5C1_HUMAN__P36542 EEF1G_HUMAN__P26641 PSMA7_HUMAN__O14818 RPLP1_HUMAN__P05386 ORC6L_HUMAN__Q9Y5N6 BRCA2_HUMAN__P51587')
            ### Make files with different levels of homology ###
            rje_blast.BLASTRun(self.log,self.cmd_list).formatDB(fasfile='AAA_hprd.fas',protein=True,force=False)
            dataset = singles[0:] + [strfam[0],mcmfam[0]]   # Should be 11 sequences (add from abcfam)
            svtfiles = []
            for r in range(10):    ### Range from 1-10 sequences from abcfam
                n = r + 1
                if makefiles:
                    mylist = abcfam[:n] + dataset[r:]
                    myseq = rje_seq.SeqList(self.log,self.cmd_list)
                    for id in mylist: myseq.seqFromFastaCmd(id,'AAA_hprd.fas')
                    myseq.saveFasta(seqfile='svt_1_%d.fas' % n,name='Teiresias')
                svtfiles.append('svt_1_%d.fas' % n)
            singles.reverse()
            dataset = singles[0:] + [abcfam[0]]   # Should be 10 sequences now
            for r in range(5):      ### Range from 1-5 sequences from each of mcmfam and strfam
                n = r + 1
                if makefiles:
                    mylist = strfam[:n] + mcmfam[:n] + dataset[2*r:]
                    myseq = rje_seq.SeqList(self.log,self.cmd_list)
                    for id in mylist: myseq.seqFromFastaCmd(id,'AAA_hprd.fas')
                    myseq.saveFasta(seqfile='svt_2_%d.fas' % n,name='Teiresias')
                svtfiles.append('svt_2_%d.fas' % n)
            dataset = singles[0:]   # Should be 9 sequences now
            for r in range(4):      ### Range from 1-4 sequences from each of mcmfam and strfam
                n = r + 1
                if makefiles:
                    mylist = abcfam[:n] + strfam[:n] + mcmfam[:n] + dataset[3*r:]
                    myseq = rje_seq.SeqList(self.log,self.cmd_list)
                    for id in mylist: myseq.seqFromFastaCmd(id,'AAA_hprd.fas')
                    myseq.saveFasta(seqfile='svt_3_%d.fas' % n,name='Teiresias')
                svtfiles.append('svt_3_%d.fas' % n)

            ### ~ Run TEIRESIAS and SLiMBuild on files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tcmd = '/home/richard/Bioware/teiresias_char -l2 -w14 -c1 -k3 -v -p'     #-i -o
            scmd = 'python /home/richard/TestBed/TestPython/slimfinder.py v=-1 preamb=F wildvar=F combamb=F'
            for svt in svtfiles:
                ## Run TEIRESIAS ##
                self.log.printLog('#RUN','Running TEIRESIAS on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s -i%s -oTEIRESIAS/%s.amb.out' % (tcmd,svt,os.path.split(svt)[1]))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished TEIRESIAS run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s TEIRESIAS_Amb %.1f seconds' % (os.path.split(svt)[1],runtime))
                ## Run SLiMBuild ##
                self.log.printLog('#RUN','Running SLiMBuild on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s seqin=%s' % (scmd,svt))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished SLiMBuild run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s SLiMBuild_Amb %.1f seconds' % (os.path.split(svt)[1],runtime))

            tcmd = '/home/richard/Bioware/teiresias_char -l2 -w14 -c1 -k3 -v -p -bequiv.txt'     #-i -o
            scmd = 'python /home/richard/TestBed/TestPython/slimfinder.py v=-1 preamb=T wildvar=F combamb=F'
            for svt in svtfiles:
                ## Run TEIRESIAS ##
                self.log.printLog('#RUN','Running TEIRESIAS on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s -i%s -oTEIRESIAS/%s.amb.out' % (tcmd,svt,os.path.split(svt)[1]))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished TEIRESIAS run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s TEIRESIAS_Amb %.1f seconds' % (os.path.split(svt)[1],runtime))
                ## Run SLiMBuild ##
                self.log.printLog('#RUN','Running SLiMBuild on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s seqin=%s' % (scmd,svt))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished SLiMBuild run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s SLiMBuild_Amb %.1f seconds' % (os.path.split(svt)[1],runtime))
            
            for svt in svtfiles:
                ## Run SLiMBuild ##
                self.log.printLog('#RUN','Running SLiMBuild on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s seqin=%s wildvar=T' % (scmd,svt))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished SLiMBuild run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s SLiMBuild_WildAmb %.1f seconds' % (os.path.split(svt)[1],runtime))
                ## Run SLiMBuild ##
                self.log.printLog('#RUN','Running SLiMBuild on %s' % os.path.split(svt)[1])
                start_time = time.time()
                os.system('%s seqin=%s wildvar=T combamb=T' % (scmd,svt))
                runtime = time.time() - start_time
                self.log.printLog('#RUN','Finished SLiMBuild run on %s' % os.path.split(svt)[1])
                self.log.printLog('#RUN','%s SLiMBuild_CombAmb %.1f seconds' % (os.path.split(svt)[1],runtime))
                
        except: self.log.errorLog('slimfinderVsTeiresias is fecked.')
#########################################################################################################################
    def slimTestSum(self):  ### Reads in a log file (self.info['InFile']) and converts into data summary table
        '''Reads in a log file (self.info['InFile']) and converts into data summary table.'''
        try:
            ### Setup ###
            headers = ['Dataset','SeqNum','UPNum','AANum','MaskAA','DimNum','SlimNum','BonfNum','SigNum','RunTime']
            outfile = rje.baseFile(self.info['InFile']) + '.summary.csv'
            rje.delimitedFileOutput(self,outfile,headers,',',rje_backup=True)
            dx = 0

            ### Memory inefficient Way! :o) ###
            for line in self.loadFromFile(self.info['InFile']):
                if not rje.matchExp('^#(\S+)\s+(\d+):(\d+):(\d+)\s+(\S.+)$',line):
                    continue
                (type,hr,min,sec,details) = rje.matchExp('^#(\S+)\s+(\d+):(\d+):(\d+)\s+(\S.+)$',line)
                if type == 'SEQ':
                    (sx,dat) = rje.matchExp('^(\d+) sequences loaded from (\S+)',details)
                    datadict = {'Dataset':dat,'SeqNum':string.atoi(sx)}
                    start = (3600 * int(hr)) + (60 * int(min)) + int(sec)
                elif type == 'UP':
                    datadict['UPNum'] = rje.matchExp('^(\d+)',details)[0]
                elif type == 'ADJ':
                    (ax,mx) = rje.matchExp('(\d+) AA to (\d+) AA',string.replace(details,',',''))
                    datadict['AANum'] = int(ax)
                    datadict['MaskAA'] = int(mx)
                elif type == 'DIM':
                    datadict['DimNum'] = rje.matchExp('(\d+) dimers',string.replace(details,',',''))[0]
                elif type == 'SLIM' and rje.matchExp('^(\d+) SLiMs',string.replace(details,',','')):
                    datadict['SlimNum'] = rje.matchExp('^(\d+) SLiMs',string.replace(details,',',''))[0]
                elif type == 'BON':
                    datadict['BonfNum'] = rje.matchExp('(\d+) equivalent',string.replace(details,',',''))[0]
                elif type == 'PROB':
                    datadict['SigNum'] = rje.matchExp('(\d+) Sig',string.replace(details,',',''))[0]
                    datadict['RunTime'] = (3600 * int(hr)) + (60 * int(min)) + int(sec) - start
                    rje.delimitedFileOutput(self,outfile,headers,',',datadict)
                    dx += 1
                    self.log.printLog('\r#OUT','%s datasets summarised' % rje.integerString(dx),newline=False,log=False)
            self.log.printLog('\r#OUT','%s datasets summarised' % rje.integerString(dx))
                    
        except:
            self.log.errorLog('Poop')
#########################################################################################################################
    def slimTest(self):     ### Makes a bunch of random datasets to run SLiMFinder on
        '''Makes a bunch of random datasets to run SLiMFinder on.'''
        try:
            ### Setup Dataset sizes ###
            datsize = []    # List of dataset sizes to analyse
            dx = 3
            while dx <= 100:
                datsize.append(dx)
                if dx < 10:
                    dx += 1
                elif dx < 18:
                    dx += 2
                elif dx < 30:
                    dx += 3
                elif dx < 50:
                    dx += 5
                else:
                    dx += 10
            reps = range(1,11)   # No. of repetitions of each dataset size
            

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TOTALLY RANDOM ONE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            import rje_seqgen
            seqgencmd = self.cmd_list + ['blastgen=F','keepnames=F','seqlen=200,800','xmerseq=uniform_aa.fas','append=F']
            for rep in reps:
                for dat in datsize:
                    rcmd = seqgencmd + ['outfile=UniformRandSeq/Uniform_%d_%d.fas' % (dat,rep),'randname=Uniform_%d_%d_' % (dat,rep),'seqnum=%d' % dat]
                    rje_seqgen.SeqGen(self.log,rcmd).seqGen()
            return


            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OLD HUMAN ONE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            ### Input Dataset ###
            #inputseq = '/home/richard/Databases/EnsEMBL/ens_HUMAN.loci.fas'
            #seqx = rje_seq.SeqCount(self,inputseq)
            allseq = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','accnr=F','seqnr=F'])

            ### Output Datasets ###
            import rje_seqgen
            import random
            seqgencmd = self.cmd_list + ['blastgen=T','keepnames=F','append=F']
            for rep in reps:
                for dat in datsize:
                    ## Human Seqs ##
                    seqi = []
                    for i in range(dat):
                        r = random.randint(1,allseq.seqNum()) - 1
                        seqi.append(allseq.seq[r])
                    seqi = rje.sortUnique(seqi,False,num=True)
                    filename = 'HsRand_%d_%d.fas' % (dat,rep)
                    allseq.saveFasta(seqi,filename)
                    ## Random Seqs ##
                    froot = 'Rand_%d_%d' % (dat,rep)
                    scmd = seqgencmd + ['randname=%s_' % froot,'seqin=%s' % filename,'outfile=%s.fas' % froot]
                    rje_seqgen.SeqGen(self.log,scmd).blastGen()
            self.log.printLog('#OUT','%d random sequence files made' % (dat*rep*2))
                        
        except:
            self.log.errorLog('SlimTest is efooked.')
#########################################################################################################################
    def wag(self):  ### Reconstructs a PAM1 matrix from a Goldman WAG matrix
        '''Reconstructs a PAM1 matrix from a Goldman WAG matrix.'''
        try:
            ### Setup ###
            wlines = self.loadFromFile('goldman.wag')
            aas = string.split(wlines[0])
            codes = string.split(wlines[1])
            rawfreqs = string.split(wlines[2])
            freq = {}
            for i in range(len(rawfreqs)):
                freq[aas[i]] = string.atof(rawfreqs[i])
            prob = {}
            for r in range(3,22):
                subs = string.split(wlines[r])
                for i in range(len(subs)):
                    prob['%s%s' % (aas[i],aas[r-2])] = string.atof(subs[i])
                    prob['%s%s' % (aas[r-2],aas[i])] = string.atof(subs[i])
            
            ### Calculate s ###
            s = 0.01
            step = 0.000001
            solve = True
            bests = 1.000000
            bestdif = -1
            while solve and s >= step:
                ## Scaler ##
                s = s - step
                self.log.printLog('\r#WAG','Considering s = %.6f; Best s = %.6f (Dif = %.6f)' % (s,bests,bestdif),log=False,newline=False)
                ## Self Subs ##
                newprobs = rje.scaledict(dict=prob,scale=s)
                toobig = False
                for a in aas:
                    newprobs['%s%s' % (a,a)] = 1.0
                    for key in prob.keys():
                        if key[0] == a:
                            newprobs['%s%s' % (a,a)] -= newprobs[key]
                            if newprobs['%s%s' % (a,a)] < 0.0:  # Overshot possibility
                                toobig = True
                                break
                    if toobig:
                        break
                if toobig:
                    continue
                #print 'PAM!!', 
                ## PAM1 ##
                dsum = 0.0
                for a in aas:
                    dsum += freq[a] * newprobs['%s%s' % (a,a)]
                dif = 0.99 - dsum
                if dif < 0:
                    dif = -dif
                if dif < bestdif or bestdif < 0:
                    bestdif = dif
                    bests = s
            ### Output best s ###
            self.log.printLog('\r#WAG','Considered all s <= 0.010000; Best s = %.6f (Dif = %.6f)' % (bests,bestdif))
            PAM = open('wag.pam','w')
            rje.writeDelimit(PAM,aas,' ')
            newprobs = rje.scaledict(dict=prob,scale=bests)
            for a in aas:
                newprobs['%s%s' % (a,a)] = 1.0
                for key in prob.keys():
                    if key[0] == a:
                        newprobs['%s%s' % (a,a)] -= newprobs[key]
            for i in range(len(aas)):
                out = [codes[i]]
                a = aas[i]
                for b in aas:
                    out.append('%.6f' % newprobs['%s%s' % (a,b)])
                rje.writeDelimit(PAM,out,' ')
            PAM.close()

        except:
            self.log.errorLog('Knackered', quitchoice=True)
#########################################################################################################################
    def disorderTest(self): ### Tests the two disorder prediction servers and compares
        '''Tests the two disorder prediction servers and compares.'''
        try:
            ### Setup ###
            seqlist = rje_seq.SeqList(self.log,self.cmd_list)

            ### Disorder ###
            for seq in seqlist.seq:
                seq.cmd_list = self.cmd_list + ['disorder=iupred','iumethod=long']
                seq.disorder()
                (seq.obj['IUPred'],seq.obj['Disorder']) = (seq.obj['Disorder'],None)
                seq.cmd_list = self.cmd_list + ['disorder=iupred','iumethod=short']
                seq.disorder()
                (seq.obj['IUPredShort'],seq.obj['Disorder']) = (seq.obj['Disorder'],None)
                seq.cmd_list = self.cmd_list + ['disorder=foldindex']
                seq.disorder()
                (seq.obj['FoldIndex'],seq.obj['Disorder']) = (seq.obj['Disorder'],None)
                seq.cmd_list = self.cmd_list + ['disorder=parse']
                seq.disorder()
                (seq.obj['DisProt'],seq.obj['Disorder']) = (seq.obj['Disorder'],None)

            ### Full Out ###
            FULL = open('disorder_full.tdt','w')
            rje.writeDelimit(FULL,['Seq','Pos','IUPredLong','IUPredShort','FoldIndex','DisProt'])
            for seq in seqlist.seq:
                for i in range(len(seq.obj['FoldIndex'].list['ResidueDisorder'])):
                    outlist = [seq.shortName(),'%d' % (i+1),
                               '%s' % seq.obj['IUPred'].list['ResidueDisorder'][i],
                               '%s' % seq.obj['IUPredShort'].list['ResidueDisorder'][i],
                               '%s' % seq.obj['FoldIndex'].list['ResidueDisorder'][i],
                               '%s' % seq.obj['DisProt'].list['ResidueDisorder'][i]]
                    rje.writeDelimit(FULL,outlist)
            FULL.close()

            ### Regions ###
            REG = open('disorder_regions.tdt','w')
            rje.writeDelimit(REG,['Seq','Method','Start','End'])
            for seq in seqlist.seq:
                for o in ['IUPred','IUPredShort','FoldIndex','DisProt']:
                    for reg in seq.obj[o].list['RegionDisorder']:
                        outlist = [seq.shortName(),o,'%d' % reg[0],'%d' % reg[1]]
                        rje.writeDelimit(REG,outlist)
            REG.close()                    
                
        except:
            self.log.errorLog('Shite. disorderTest() is buggered', quitchoice=True)
#########################################################################################################################
    def ygobOrthologues(self):  ### Makes orthologue alignments from YGOB files
        '''Makes orthologue alignments from YGOB files.'''
        try:
            ### (Re)Format Database ###
            spec_re = {'^(CAGL\S+)':'CANGA','^(Scas\S+)':'SACCA','^(Y\S+)':'YEAST','^(A\S+)':'ASHGO',
                       '^(KLL\S+)':'KLULA','^(Kwal\S+)':'KLUWA','^(Sklu\S+)':'SACKL'}
            ygob = rje.makePath('C:/Documents and Settings/redwards/My Documents/redwards/Databases/YGOB/ygob.fas',True)
            if True:
                NEWFAS = open(ygob,'w')
                YGOB = open(rje.makePath('C:/Documents and Settings/redwards/My Documents/redwards/Databases/YGOB/Proteins.fsa',True),'r')
                line = YGOB.readline()
                yx = 0
                (name,seq) = ('','')
                ydict = {}
                while line:
                    if line[0] == '>':  # New sequence
                        if name:
                            ydict[name] = seq
                        match = ()
                        for spre in spec_re.keys():
                            match = rje.matchExp(spre,line[1:])
                            if match:
                                break
                        if not match:
                            self.log.errorLog('Cannot reformat %s' % line[1:],printerror=False)
                            raise ValueError
                        NEWFAS.write('>ygob_%s__%s\n' % (spec_re[spre],match[0]))
                        name = 'ygob_%s__%s' % (spec_re[spre],match[0])
                        seq = ''
                        yx += 1
                        self.log.printLog('\r#FAS','Reformatting YGOB fasta: %s' % rje.integerString(yx),newline=False,log=False)
                    else:
                        NEWFAS.write(string.replace(line,'*',''))
                        seq += rje.chomp(string.replace(line,'*',''))
                    line = YGOB.readline()
                ydict[name] = seq
                NEWFAS.close()
                YGOB.close()
                self.log.printLog('\r#FAS','Reformatting YGOB fasta: %s' % rje.integerString(yx))
                #X#rje_blast.formatDB(ygob,'c:/bioware/blast/',True,self.log)

            ### Make Orthologues ###
            alndir = rje.makePath('C:/Documents and Settings/redwards/My Documents/redwards/Databases/YGOB/ALN/')
            if not os.path.exists(alndir):
                os.mkdir(alndir)
            ylines = self.loadFromFile(rje.makePath('C:/Documents and Settings/redwards/My Documents/redwards/Databases/YGOB/Pillars.tab',True),chomplines=True)
            yx = 0.0
            for y in ylines:
                y = string.replace(y,'undef','---')
                yx += 100.0
                orths1 = string.split(y)[:-3]
                orths2 = string.split(y)[3:]
                orths2.reverse()
                if len(orths1) != 7 or len(orths2) != 7:
                    self.log.errorLog('Cannot parse "%s"' % y,printerror=False)
                    self.deBug(orths1)
                    self.deBug(orths2)
                    continue
                for olist in [orths1,orths2]:
                    if olist[2] != '---':   ### Yeast sequence in column
                        plist = [olist[2]] + olist[:2] + olist[3:]
                        while '---' in plist:
                            plist.remove('---')
                        fascmd = []
                        for p in plist:
                            match = ()
                            for spre in spec_re.keys():
                                match = rje.matchExp(spre,p)
                                if match:
                                    break
                            if not match:
                                self.log.errorLog('Cannot reformat %s' % p,printerror=False)
                                self.deBug(y)
                                self.deBug(plist)
                                raise ValueError
                            fascmd.append('ygob_%s__%s' % (spec_re[spre],match[0]))
                        if not fascmd:
                            self.deBug(plist)
                            continue
                        seqfile='%s%s.orthaln.fas' % (alndir,plist[0])
                        YSEQ = open(seqfile,'w')
                        for f in fascmd:
                            if ydict.has_key(f):
                                YSEQ.write('>%s\n%s\n' % (f,ydict[f]))
                            else:
                                self.log.errorLog('No sequence read for %s!' % f,printerror=False)
                            #X#yseq.seqFromFastaCmd(f,ygob)
                        YSEQ.close()
                        #X#self.deBug(fascmd)
                self.log.printLog('\r#ORTH','Generating orthologue datasets: %.2f%%' % (yx/len(ylines)),newline=False,log=False)
            self.log.printLog('\r#ORTH','Generating orthologue datasets: %.2f%%' % (yx/len(ylines)))

            #C. glabrata Position 1    #
            #S. castellii Position 1
            #S. cerevisiae Position 1
            #A. gossypii
            #K. lactis
            #K. waltii
            #S. kluyveri
            #S. cerevisiae Position 2
            #S. castellii Position 2
            #C. glabrata Position 2                
        except:
            self.log.errorLog('Shite. YGOB Orthologues() is buggered.', quitchoice=True)
#########################################################################################################################
    def ygobOrthologues2(self):  ### Makes orthologue alignments from YGOB files
        '''Makes orthologue alignments from YGOB files.'''
        try:
            ### Setup ###
            ygobdir = ''    #rje.makePath('C:/Documents and Settings/redwards/My Documents/redwards/Databases/YGOB')
            orthdir = rje.makePath(ygobdir + 'ORTH')
            alndir = rje.makePath(ygobdir + 'ALN')
            if not os.path.exists(orthdir):
                os.mkdir(orthdir)

            ### Rename ###
            if True:
                ygob1 = glob.glob('%s*.orthaln.fas' % alndir)
                ygob2 = glob.glob('%s*.orth.fas' % orthdir)
                for y in ygob1 + ygob2:
                    orth = string.replace(y,'orth','ygob')
                    os.rename(y,orth)
                self.log.printLog('#MOVE','Renamed %d files.' % len(ygob1))
                return

            ### Align ###
            ygob1 = glob.glob('%s*.orthaln.fas' % alndir)
            ygob2 = glob.glob('%s*.orth.fas' % orthdir)
            for orth in ygob2:
                file = os.path.basename(orth)
                aln = alndir + string.replace(file,'orth','orthaln')
                if aln not in ygob1:    # Already done
                    scmd = self.cmd_list + ['seqin=%s' % orth]
                    seqs = rje_seq.SeqList(log=self.log,cmd_list=scmd)
                    seqs.muscleAln(outfile=aln,mapseq=True)
            self.log.printLog('#ALN','Aligned %d files.' % len(ygob2))
                
        except:
            self.log.errorLog('Shite. YGOB Orthologues() is buggered.', quitchoice=True)
#########################################################################################################################
    def ensLoci(self):  ### Calculates stats on EnsLoci
        '''Calculates stats on EnsLoci.'''
        try:
            ### Setup ###
            ens_known = {}  # Dictionary of {species:{gene:count}}
            ens_novel = {}  # Dictionary of {species:{gene:count}}

            ### EnsEMBL Known ###
            for known in glob.glob('*.known.fas'):
                species = rje.matchExp('^ens_(\S+)\.known',known)[0]
                ens_known[species] = {}
                ens_novel[species] = {}
                ENS = open(known,'r')
                kx = 0
                while True:
                    line = ENS.readline()
                    if line in ['',None]:
                        break
                    match = rje.matchExp('^>(\S+)\s.+gene:(\S+)\s',line)
                    if match:
                        kx += 1
                        (seqname,gene) = match
                        if ens_known[species].has_key(gene):
                            ens_known[species][gene] += 1
                        else:
                            ens_known[species][gene] = 1
                        self.log.printLog('\r#LOCI','%s: %s known; %s loci' % (species,rje.integerString(kx),rje.integerString(len(ens_known[species]))),newline=False,log=False)
                    elif line[:1] == '>':
                        self.log.errorLog('No gene for: %s' % line[1:], printerror=False)
                self.log.printLog('\r#LOCI','%s: %s known; %s loci' % (species,rje.integerString(kx),rje.integerString(len(ens_known[species]))))
                ENS.close()

            ###  EnsEMBL Novel ###
            for novel in glob.glob('*.novel.fas'):
                species = rje.matchExp('^ens_(\S+)\.novel',novel)[0]
                ENS = open(novel,'r')
                kx = sum(ens_known[species].values())
                nx = 0
                while True:
                    line = ENS.readline()
                    if line in ['',None]:
                        break
                    match = rje.matchExp('^>(\S+)\s.+gene:(\S+)\s',line)
                    if match:
                        nx += 1
                        (seqname,gene) = match
                        if ens_known[species].has_key(gene):
                            ens_known[species][gene] += 1
                        elif ens_novel[species].has_key(gene):
                            ens_novel[species][gene] += 1
                        else:
                            ens_novel[species][gene] = 1
                        self.log.printLog('\r#LOCI','%s: %s known; %s known loci; %s novel; %s novel loci' % (species,rje.integerString(kx),rje.integerString(len(ens_known[species])),rje.integerString(nx),rje.integerString(len(ens_novel[species]))),newline=False,log=False)
                    elif line[:1] == '>':
                        self.log.errorLog('No gene for: %s' % line[1:], printerror=False)
                self.log.printLog('\r#LOCI','%s: %s known; %s known loci; %s novel; %s novel loci' % (species,rje.integerString(kx),rje.integerString(len(ens_known[species])),rje.integerString(nx),rje.integerString(len(ens_novel[species]))))
                ENS.close()
                
        except:
            self.log.errorLog('Wah wah wah! Cry like a baby - ensLoci() is busted.',quitchoice=True)
#########################################################################################################################
    def ensLociGo(self):    ### Partitions GO Datasets by Type and outputs data 
        '''Makes GO Datasets from EnsEMBL Loci data and BioMart download.'''
        try:
            ### Setup Files and Directories ###
            #os.mkdir('EnsGO_CC')
            #os.mkdir('EnsGO_BP')
            #os.mkdir('EnsGO_MF')
            DO = open('go_map.do','a')
            #DO.write('gen desc = dataset\n\n')
            
            ### Read in dictionary ###
            godesc = {}     # Dictionary of code:description
            golines = self.loadFromFile('go_ids_05-10-06.txt')
            (cx,px,fx) = (0,0,0)
            lx = 0.0
            for line in golines:
                lx += 100.0
                if line[:1] == '!':
                    continue
                go = string.split(line)
                if len(go) >= 3:
                    desc = string.join(go[1:-1])
                    id = string.replace(go[0],':','_')
                    type = go[-1]
                    if type not in ['P','C','F']:
                        self.log.errorLog('Fuck. What is "%s"?' % rje.chomp(line),printerror=False)
                        continue
                    gofile = 'GO_DATASETS/%s.fas' % id
                    if not os.path.exists(gofile):
                        self.log.errorLog('%s not found' % gofile,printerror=False)
                        continue
                    if type == 'C':
                        os.rename(gofile,'EnsGO_CC/%s.fas' % id)
                        cx += 1
                    elif type == 'P':
                        os.rename(gofile,'EnsGO_BP/%s.fas' % id)
                        px += 1
                    elif type == 'F':
                        os.rename(gofile,'EnsGO_MF/%s.fas' % id)
                        fx += 1
                    else:
                        continue
                    DO.write('replace desc = "%s" if dataset == "%s"\n' % (desc,id))
                self.log.printLog('\r#GO','%s CC; %s MF; %s BP; %.1f%%' % (rje.integerString(cx),rje.integerString(fx),rje.integerString(px),lx/len(golines)),log=False,newline=False)
            self.log.printLog('\r#GO','%s CC; %s MF; %s BP; %.1f%%' % (rje.integerString(cx),rje.integerString(fx),rje.integerString(px),lx/len(golines)))
            DO.close()

        except:
            self.log.errorLog('Wah wah wah! Cry like a baby - ensLociGo() is busted.',quitchoice=True)
#########################################################################################################################
    def OLDensLociGo(self):    ### Makes GO Datasets from EnsEMBL Loci data and BioMart download
        '''
        Makes GO Datasets from EnsEMBL Loci data and BioMart download.
        This method has been run and retired to make way for the next round.
        '''
        try:
            ### Read in dictionary ###
            godesc = {}     # Dictionary of code:description
            genego = {}     # Dictionary of gene:[gocodelist]
            protgo = {}     # Dictionary of prot:[gocodelist]
            golines = self.loadFromFile('ens_human_go.csv')[1:]
            lx = 0.0
            for line in golines:
                lx += 100.0
                go = string.split(line,',')
                if len(go) >= 5:
                    godesc[go[3]] = go[4]
                    if not genego.has_key(go[0]):
                        genego[go[0]] = []
                    if go[3] not in genego[go[0]]:
                        genego[go[0]].append(go[3])
                    if not protgo.has_key(go[2]):
                        protgo[go[2]] = []
                    if go[3] not in protgo[go[2]]:
                        protgo[go[2]].append(go[3])
                self.log.printLog('\r#GO','%s GO; %s Genes; %s Proteins; %.1f%%' % (rje.integerString(len(godesc)),rje.integerString(len(genego)),rje.integerString(len(protgo)),lx/len(golines)),log=False,newline=False)
            self.log.printLog('\r#GO','%s GO; %s Genes; %s Proteins; %.1f%%' % (rje.integerString(len(godesc)),rje.integerString(len(genego)),rje.integerString(len(protgo)),lx/len(golines)))

            ### Make GO Files ###
            genego = {} #!# Save memory #!#
            os.mkdir('GO_DATASETS')
            seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=ens_HUMAN.loci.fas'])
            sx = 0.0
            while seqlist.seq:
                seq = seqlist.seq.pop(0)
                pep = rje.matchExp('pep:(\S+)',seq.info['Name'])[0]
                if protgo.has_key(pep):
                    for go in protgo.pop(pep):
                        gofile = 'GO_DATASETS/%s.fas' % string.replace(go,':','_')
                        GO = open(gofile,'a')
                        GO.write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
                        GO.close()
                sx += 100.0
                self.log.printLog('\r#GO','Making GO datasets %.1f%%' % (sx/seqlist.seqNum()),log=False,newline=False)
            self.log.printLog('\r#GO','Making GO datasets complete!')
        
            ### Make GO Desc ###
            DO = open('go_map.do','w')
            gx = 0
            for go in godesc.keys():
                gx += 1
                DO.write('replace desc = "%s" if dataset == "GO_%s"\n' % (godesc[go],go))
                self.log.printLog('\r#GO','%s GO terms.' % rje.integerString(gx),newline=False,log=False)
            self.log.printLog('\r#GO','%s GO terms.' % rje.integerString(gx))
            DO.close()

        except:
            self.log.errorLog('Wah wah wah! Cry like a baby - ensLociGo() is busted.',quitchoice=True)
#########################################################################################################################
    def elmDup(self):   ### Identifies "duplicate" proteins from different species using orthaln.fas files
        '''Identifies "duplicate" proteins from different species using orthaln.fas files.'''
        try:
            gopherdir = rje.makePath('../Gopher/ALN')
            enseq = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T'])
            ensprot = []
            for seq in enseq.seq:
                ensprot.append(seq.shortName())

            dup = [['Prot1','Prot2']] 
            for seq in enseq.seq:
                name = seq.shortName()
                acc = seq.info['AccNum']
                fas = gopherdir + acc + '.orthaln.fas'
                if not os.path.exists(fas):
                    continue
                orths = rje_seq.SeqList(self.log ,self.cmd_list+['seqin=%s' % fas,'autofilter=F','autoload=T'])
                for oseq in orths.seq:
                    o = oseq.shortName()
                    if o != name and o in ensprot:
                        dup.append([name,o])

            DUP = open('elm_dup.tdt','w') 
            for pair in dup:
                DUP.write('%s\n' % string.join(pair,'\t'))
            DUP.close()

        except:
            self.log.errorLog('Wah wah wah! Cry like a baby - elmDup() is busted.',quitchoice=True)
#########################################################################################################################
    def elmAln(self):   ### Cleans up and copies relevant ELM alignments
        '''Cleans up and copies relevant ELM alignments'''
        try:
            ### Setup ###
            gopherdir = rje.makePath('../Gopher/ALN')
            ygobdir = rje.makePath('/home/richard/YGOB/YGOB_ALN/')
            alndir = rje.makePath('../Alignments')
            humprot = []    # List of human proteins in alignments - want to test for novel ELMs in rest

            ### AccNum ###
            enseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=elm_prot.fas','autoload=T'])
            elmacc = []
            for seq in enseq.seq:
                elmacc.append(seq.info['AccNum'])
            
            ### Alignments ###
            for acc in elmacc:
                ygob = ygobdir + acc + '.ygobaln.fas'
                goph = gopherdir + acc + '.orthaln.fas'
                new = alndir + acc + '.elmaln.fas'
                if os.path.exists(new):
                    scmd = self.cmd_list + ['seqin=%s' % new,'autoload=T','autofilter=F']
                    gseq = rje_seq.SeqList(self.log,scmd)
                    for seq in gseq.seq:
                        if seq.info['SpecCode'] == 'HUMAN':
                            humprot.append(seq.shortName())
                elif os.path.exists(ygob):    # Yeast sequences
                    os.system('cp %s %s' % (ygob,new))
                elif os.path.exists(goph):   # Not yeast
                    scmd = self.cmd_list + ['seqin=%s' % goph,'seqout=%s' % new,'badspec=YEAST','autoload=T','autofilter=T','degap=T']
                    gseq = rje_seq.SeqList(self.log,scmd)
                    gseq.info['Name'] = new
                    gseq._checkAln(aln=True,realign=True)
                    gseq.saveFasta()
                    for seq in gseq.seq:
                        if seq.info['SpecCode'] == 'HUMAN':
                            humprot.append(seq.shortName())
                else:
                    open('no_aln.txt','a').write('%s\n' % acc)

            ### HumProt ###
            open('elm_human.txt','w').write(string.join(humprot,'\n'))
            
        except:
            self.log.errorLog('Wah wah wah! Cry like a baby - elmAln() is busted.',quitchoice=True)
#########################################################################################################################
    def elmDom(self):   ### Uses mapping.tdt and UniProt to make file of features and proteins without information
        '''Uses mapping.tdt and UniProt to make file of features and proteins without information.'''
        try:
            ### Read in Mapping ###
            mlines = self.loadFromFile('mapping.tdt',chomplines=True)
            headers = string.split(mlines[0],'\t')
            mappings = {}
            for line in mlines:
                ## Read line ##
                data = string.split(line,'\t')
                if not data or data == headers:
                    continue
                if len(data) != len(headers):
                    self.log.errorLog('Data Mismatch:\n%s\n <> \n%s' % (data,headers),printerror=False,quitchoice=True)
                    continue
                ## Assign to temp dictionary ##
                datadict = {}
                for h in range(len(headers)):
                    datadict[headers[h]] = data[h]
                ## Add to mappings ##
                elmprot = datadict['Hit']
                if elmprot == 'None' or datadict['Query_Len'] != '100.0' or datadict['Hit_Len'] != '100.0':
                    continue
                (id,acc) = rje.matchExp('^(\S+_\S+)__(\S+)',datadict['Query'])
                mappings[acc] = elmprot
            open('elm_map100.txt','w').write(string.join(mappings.values(),'\n'))
            self.log.printLog('\r#MAP','Read %s 100%% len mappings' % (rje.integerString(len(mappings))))

            ### Get UniProt ###
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
            entrydict = uniprot.accDict(mappings.keys())
            self.log.printLog('\r#UNI','Extracted UniProt entries for %s AccNum' % (rje.integerString(len(entrydict))))

            ### Output Features ###
            headers = ['Name','AccNum','Type','Start','End','Desc']
            rje.delimitedFileOutput(self,'elm_ft.tdt',headers,'\t',rjebackup=True)
            for acc in rje.sortKeys(mappings):
                entry = entrydict[acc]
                ## Make dictionary of {start:{end:[features]}}
                ft_dict = {}
                for ft in entry.list['Feature']:
                    if ft['Type'] not in ['TRANSMEM','DOMAIN']:
                        continue
                    ft_start = ft['Start']
                    if not ft_dict.has_key(ft_start):
                        ft_dict[ft_start] = {}
                    ft_end = ft['End']
                    if not ft_dict[ft_start].has_key(ft_end):
                        ft_dict[ft_start][ft_end] = []
                    ft_dict[ft_start][ft_end].append(ft)
                ## Sort and output ##
                for ft_start in rje.sortKeys(ft_dict):
                    for ft_end in rje.sortKeys(ft_dict[ft_start]):
                        for ft in ft_dict[ft_start][ft_end]:
                            datadict = {'Name':mappings[acc],'AccNum':acc}
                            for fk in ['Type','Start','End','Desc']:
                                datadict[fk] = '%s' % ft[fk]
                            rje.delimitedFileOutput(self,'elm_ft.tdt',headers,'\t',datadict)    
            self.log.printLog('\r#FT','Extracted UniProt features output to elm_ft.tdt')

        except:
            self.log.errorLog('Shitsticks!')
#########################################################################################################################
    def ppi1433(self):   ### Makes a table of 14-3-3 interacting proteins from HPRD
        '''Makes a table of 14-3-3 interacting proteins from HPRD.'''
        try:### Setup ###
            isoforms = {'Beta':'YWHAB','Epsilon':'YWHAE','Eta':'YWHAH','Gamma':'YWHAG','Sigma':'SFN','Theta':'YWHAQ','Zeta':'YWHAZ'}
            headers = ['Gene','Species','Protein','HPRD_ID','UniProt','GenBank','Beta','Epsilon','Eta','Gamma','Sigma','Theta','Zeta']
            gene_data = {}
            ### Read in ###
            for iso in rje.sortKeys(isoforms):
                file = '%s_hprd.fas' % isoforms[iso]
                for seq in rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % file]).seq:
                    (gene,hprd,gb,sp) = rje.matchExp('\[Gene:(\S+) HPRD:(\S+); gb:(\S+); sp:(\S+)\]',seq.info['Name'])
                    if gene == '-': gene = 'HPRD_%s' % hprd
                    if gene_data.has_key(gene): gene_data[gene][iso] = '1'
                    else: gene_data[gene] = {'Gene':gene,'Species':'Human','Protein':seq.shortName(),'HPRD_ID':hprd,
                                             'UniProt':sp,'GenBank':gb,iso:'1'}
            ### Output ###
            out = '1433_ppi.tdt'
            rje.backup(self,out)
            rje.delimitedFileOutput(self,out,headers,'\t')
            for gene in rje.sortKeys(gene_data): rje.delimitedFileOutput(self,out,headers,'\t',gene_data[gene])
            self.log.printLog('#OUT','Output for %d interactors complete' % len(gene_data))
        except: self.log.errorLog('Wham! Bam! This is fucked. Jerk.')
#########################################################################################################################
    def sfTarGZ(self):  ### Cleanup SLiMFinder targz data
        '''Cleanup SLiMFinder targz data.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bases = []
            self.progLog('\r#BASE','Counting datasets...')
            for g in glob.glob('*.*'):
                self.progLog('\r#BASE','Counting datasets: %s' % rje.integerString(len(bases)))
                base = string.join(string.split(g,'.')[:2],'.')
                if base not in bases: bases.append(base)
            ### ~ [1] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (bx,btot) = (0.0,len(bases))
            for base in bases:
                self.progLog('\r#TARGZ','Cleaning TARGZ for %s datasets: %.3f%%' % (rje.integerString(len(bases)),(bx/btot))); bx += 100.0
                done = '%s.done' % base
                self.deBug('\n%s: %s' % (base,os.path.exists(done)))
                if os.path.exists(done): continue
                oldtar = '%s.tar.gz' % base
                newtar = '%s.l5w2o2a1.FreqConsDisComp-5-8.tar.gz' % base
                for gz in glob.glob('%s.*.gz' % base):
                    if string.split(gz,'.')[-2] in ['pickle','tar']: continue
                    os.popen('gunzip %s' % gz).read()
                if os.path.exists(oldtar): os.unlink(oldtar)
                if os.path.exists(newtar): x = os.popen('tar -xzf %s' % newtar).read()
                for b in glob.glob('SLiMFinder090817/*'):   #glob.glob('SLiMFinder090817/%s.*' % base):
                    if string.split(b,'.')[-1] not in ['csv','upc','fas','tdt','txt'] or string.split(b,'.')[-2] in ['phydis','compare']: os.unlink(b)
                    else:
                        n = os.path.basename(b)
                        if os.path.exists(n): os.unlink(n)
                        os.rename(b,n)
                        self.deBug('-> %s: %s' % (n,os.path.exists(n)))
                    self.deBug('<- %s: %s' % (b,os.path.exists(b)))
                #x#os.popen('mv SLiMFinder090817/%s.* .' % base).read()
                keepers = []
                for b in glob.glob('%s.*' % base):
                    if string.split(b,'.')[-1] not in ['csv','upc','fas','tdt','txt','gz','pickle'] or string.split(b,'.')[-2] in ['phydis','compare']: os.unlink(b)
                    elif string.split(b,'.')[-1] not in ['gz']: keepers.append(b)
                self.deBug('>> %s <<' % string.join(keepers,' | '))
                open('%s.done' % base,'w').write(string.join(keepers,'\n'))
                if os.path.exists(newtar): os.unlink(newtar)
                #x# Later? #x# x = os.popen('tar -czf %s %s' % (newtar,string.join(keepers))).read()
                self.deBug(glob.glob('SLiMFinder090817/*'))
            self.printLog('\r#TARGZ','Cleaning TARGZ for %s datasets complete!' % (rje.integerString(len(bases))))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def sfTidy(self):   ### Grand Cleanup of SLiMFinder results data
        '''
        Grand Cleanup of SLiMFinder results data.
        - *.tar.gz to be moved to Archive
        - *.upc, *.slimdb, *.pickle.gz to be moved to Build
        - *.fas, *.occ.csv, *.cloud.txt & *.dis.tdt to be moved to Sig (if occ.csv found) + copies of *.upc 
        - *.fas & *.dis.tdt to be moved to NoSig (if no occ.csv found) + copies of *.upc 
        - *.blast and *.phydis.txt to be deleted
        - rseq* and rupc* to be moved to Random
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outdir = '/scratch/RJE_Filestore/SBSBINF/Projects/Motifs/2010-01_Human-SLiMFinder/'
            rje.mkDir(self,outdir)
            moved = {}
            for subdir in ['Archive','Build','Sig','NoSig','Random','Delete','GO']:
                rje.mkDir(self,'%s/HumSF09_%s/' % (outdir,subdir))
                moved[subdir] = 0
            ## ~ [0b] Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            extdict = {'tar.gz':'Archive'}
            for ext in ['upc','slimdb','pickle.gz']: extdict[ext] = 'Build'     # UPC is special!
            for ext in ['occ.csv','cloud.txt']: extdict[ext] = 'Sig'    
            for ext in ['masked.fas','mapping.fas','motifaln.fas','maskaln.fas','dis.tdt']: extdict[ext] = 'NoSig'
            for ext in ['blast','phydis.txt','done']: extdict[ext] = 'Delete'  # to be deleted
            ## ~ [0c] Input data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            indir = ['/scratch/RJE_Filestore/SBSBINF/Projects/Motifs/2009-09_SLiMFRAP/SLiMFinder090817/',
                     '/scratch/RJE_Filestore/SBSBINF/Projects/Motifs/2009-10_BinSF/SLiMFinder091008/',
                     '%s/HumSF09_Unknown/' % outdir]
            glist = []      # List of files from indir
            slist = []      # List of files awaiting movement depending on Sig Status of base
            bases = []      # Dataset bases
            sig = []        # List of bases with occ.csv file
            self.progLog('\r#FILES','Getting files...')
            for i in indir: glist += rje.getFileList(self,i,['*.*'],subfolders=False,summary=True,filecount=len(glist))
            self.printLog('\r#FILES','Getting files: %5s' % rje.integerString(len(glist)))
            
            ### ~ [1] Process files 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0.0; gtot = len(glist)
            while glist:
                self.progLog('\r#PASS1','Processing files: %.3f%% (%d Sig)' % (gx/gtot, len(sig))); gx += 100.0
                g = glist.pop(0); file = os.path.basename(g)
                base = string.join(string.split(file,'.')[:2],'.')
                ext = string.join(string.split(file,'.')[2:][-2:],'.')
                if string.split(ext,'.')[0][-3:] in ['y2h','bin','ppi','dom']:
                    base = '%s.%s' % (base,string.split(ext,'.')[0])
                    ext = string.join(string.split(ext,'.')[1:],'.')
                if string.split(base,'.')[-1][-3:] not in ['y2h','bin','ppi','dom']:
                    ext = '%s.%s' % (string.split(base,'.')[-1],ext)
                    base = string.join(string.split(base,'.')[:-1],'.')
                newdir = ''
                ## ~ [1a] Classify data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if base[:4] in ['rseq','rupc']: newdir = 'Random'; ext = ''
                elif base[:3] in ['GO_']: newdir = 'GO'; ext = ''
                elif ext in extdict: newdir = extdict[ext]
                else:
                    newdir = extdict[ext] = rje.choice('Fate for unknown "%s" extension?' % ext,'Unknown',confirm=True)
                    rje.mkDir(self,'%s/HumSF09_%s/' % (outdir,extdict[ext]))
                    if newdir not in moved: moved[newdir] = 0
                if ext == 'occ.csv' and base not in sig: sig.append(base)
                ## ~ [1b] Move/copy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not newdir: continue
                newfile = '%s/HumSF09_%s/%s' % (outdir,newdir,file)
                if ext == 'upc':
                    rje.fileTransfer(fromfile=g,tofile=newfile,deletefrom=False,append=False)
                    slist.append(g)
                    moved[newdir] += 1
                elif newdir in ['Sig','NoSig']: slist.append(g)
                else: os.rename(g,newfile); moved[newdir] += 1

            ### ~ [2] Process files 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(slist)
            self.printLog('\r#PASS1','Processed %s files: %s Sig datasets; %s files remain' % (rje.integerString(gtot),rje.integerString(len(sig)),rje.integerString(stot)))
            while slist:
                self.progLog('\r#PASS2','Processing files: %.3f%%' % (sx/stot)); sx += 100.0
                g = slist.pop(0); file = os.path.basename(g)
                base = string.join(string.split(file,'.')[:2],'.')
                ## ~ [1a] Classify data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if base in sig: newdir = 'Sig'
                else: newdir = 'NoSig'
                ## ~ [1b] Move/copy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newfile = '%s/HumSF09_%s/%s' % (outdir,newdir,file)
                os.rename(g,newfile); moved[newdir] += 1
            self.printLog('\r#PASS2','Processed %s files: %s Sig datasets' % (rje.integerString(gtot),rje.integerString(len(sig))))
            for subdir in rje.sortKeys(moved):
                self.printLog('#%s' % subdir,'%s files moved to HumSF09_%s/' % (rje.integerString(moved[subdir]),subdir))

            ### ~ [3] Clean Random & GO files from before ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Random Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rlist = []
            self.progLog('\r#FILES','Getting files...')
            for subdir in ['Archive','Build','Sig','NoSig']:
                rlist += rje.getFileList(self,'%s/HumSF09_%s/' % (outdir,subdir),['rseq*','rupc*'],subfolders=False,summary=True,filecount=len(glist))
            self.printLog('\r#FILES','Getting random files: %5s' % rje.integerString(len(rlist)))
            rx = 0.0; rtot = len(rlist); mx = 0
            while rlist:
                self.progLog('\r#RAND','Cleaning random files: %.3f%%' % (rx/rtot)); rx += 100.0
                g = rlist.pop(0); file = os.path.basename(g)
                ## ~ [3a] Move/copy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newfile = '%s/HumSF09_Random/%s' % (outdir,file)
                os.rename(g,newfile); mx += 1
            self.printLog('\r#RAND','Cleaning random files: moved %s' % rje.integerString(mx))
            ## ~ [3b] GO Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rlist = []
            self.progLog('\r#FILES','Getting files...')
            for subdir in ['Archive','Build','Sig','NoSig']:
                rlist += rje.getFileList(self,'%s/HumSF09_%s/' % (outdir,subdir),['GO_*'],subfolders=False,summary=True,filecount=len(glist))
            self.printLog('\r#FILES','Getting GO files: %5s' % rje.integerString(len(rlist)))
            rx = 0.0; rtot = len(rlist); mx = 0
            while rlist:
                self.progLog('\r#RAND','Cleaning GO files: %.3f%%' % (rx/rtot)); rx += 100.0
                g = rlist.pop(0); file = os.path.basename(g)
                ## ~ [3a] Move/copy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newfile = '%s/HumSF09_GO/%s' % (outdir,file)
                os.rename(g,newfile); mx += 1
            self.printLog('\r#RAND','Cleaning GO files: moved %s' % rje.integerString(mx))
                
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def prpSI(self):    ### PRP SI maker
        '''PRP SI maker.'''
        try:### Load ###
            ass = rje.dataDict(self,'assembly.tdt',['CONS'],['EST'],lists=True)
            print len(ass)
            tlr = rje.dataDict(self,'tlrfam.tdt',getheaders=True)
            print rje.sortKeys(tlr)
            igs = rje.dataDict(self,'igsfam.tdt')
            print rje.sortKeys(igs)
            ### Setup output ###
            head = tlr.pop('Headers')[:-1]
            head.insert(1,'CONS')
            head += ['Fam','Candidate']
            data = {}
            rje.delimitedFileOutput(self,'prpsi.tdt',head,rje_backup=True)
            for cons in rje.sortKeys(tlr):
                if cons == 'EST': continue
                tlr[cons]['CONS'] = cons
                tlr[cons]['Candidate'] = tlr[cons]['TLR'] 
                for est in ass[cons]['EST']:
                    data[est] = rje.combineDict({'EST':est,'Fam':'TLR'},tlr[cons],overwrite=False)
                    self.deBug(data[est])
            for cons in rje.sortKeys(igs):
                if cons == 'EST': continue
                igs[cons]['CONS'] = cons
                igs[cons]['Candidate'] = igs[cons]['IGS'] 
                for est in ass[cons]['EST']:
                    data[est] = rje.combineDict({'EST':est,'Fam':'IGSF'},igs[cons],overwrite=False)
                    self.deBug(data[est])
            ### Output ###
            for est in rje.sortKeys(data): rje.delimitedFileOutput(self,'prpsi.tdt',head,datadict=data[est])
                

        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def SF09x1(self): ### HumSF09 Tidy
        '''HumSF09 Tidy.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seq.SeqList(self.log,self.cmd_list)
            seqdict = seqlist.seqNameDic()
            ofile = 'humsf09_slimfinder.occ.csv'
            occdict = rje.dataDict(self,ofile,['Hub','Dataset','Rank','Pattern','Seq','Match','Cons','LocID'],'All',getheaders=True)
            headers = occdict.pop('Headers')
            ## ~ [0a] ~ New File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newfile = 'humsf09_slimfinder.newocc.csv'
            rje.delimitedFileOutput(self,newfile,headers,rje_backup=True)
            done = []
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ox = 0.0; otot = len(occdict)
            for okey in rje.sortKeys(occdict):
                self.progLog('\r#OCC','Processing Occ: %.2f%%' % (ox/otot)); ox += 100.0
                occ = occdict[okey]
                try:
                    seq = seqdict[occ['Seq']]
                    occ['Start_Pos'] = seq.info['Sequence'].find(occ['Match'])
                    occ['End_Pos'] = occ['Start_Pos'] + string.atoi(occ['End_Pos']) - 1
                    check = string.join([occ['Hub'],occ['Dataset'],occ['Pattern'],occ['Seq'],'%d %d' % (occ['Start_Pos'],occ['End_Pos'])])
                    if check in done: self.printLog('#DUP','Oops. Double hit: %s' % check)
                    else: done.append(check)
                except: self.errorLog('Problem mapping %s occurrence %s %s' % (occ['Seq'],occ['Pattern'],occ['Match']))
                rje.delimitedFileOutput(self,newfile,headers,datadict=occ)
            self.printLog('\r#OCC','%s Occ processed -> %s' % (rje.integerString(otot),rje.integerString(len(done))))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def SF09x2(self): ### HumSF09 Tidy 2
        '''HumSF09 Tidy.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #seqlist = rje_seq.SeqList(self.log,self.cmd_list)
            #seqdict = seqlist.seqNameDic()
            ofile = 'humsf09_slimfinder.occ.csv'
            occdict = rje.dataDict(self,ofile,['Hub','Dataset','Rank','Pattern','Seq','Match','Start_Pos','End_Pos'],'All',getheaders=True)
            oheaders = occdict.pop('Headers')
            oheaders.insert(oheaders.index('Seq'),'Spoke')
            pfile = '../HumSF09_PPI/pingu.pairwise.tdt'
            ppidict = rje.dataDict(self,pfile,['Hub','Spoke'],'All',getheaders=True)
            pheaders = ppidict.pop('Headers') + ['ppi','bin','com','y2h']
            genemap = rje.dataDict(self,pfile,['SpokeSeq'],['Spoke'],getheaders=False)
            ## ~ [0a] ~ New Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newofile = 'humsf09_slimfinder.newocc.csv'
            rje.delimitedFileOutput(self,newofile,oheaders,rje_backup=True)
            newpfile = '../HumSF09_PPI/humsf09_ppi.pairwise.tdt'
            rje.delimitedFileOutput(self,newpfile,pheaders,rje_backup=True)
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ox = 0.0; otot = len(occdict); genex = 0
            for okey in rje.sortKeys(occdict):
                self.progLog('\r#OCC','Processing Occ: %.2f%%' % (ox/otot)); ox += 100.0
                occ = occdict[okey]
                if occ['Hub'][:3] == 'GO_': continue
                try:
                    occ['Spoke'] = genemap[occ['Seq']]['Spoke']
                    genex += 1
                except: self.errorLog('Problem mapping gene for %s' % (occ['Seq']))
                rje.delimitedFileOutput(self,newofile,oheaders,datadict=occ)
            self.printLog('\r#OCC','%s Occ processed -> %s' % (rje.integerString(otot),rje.integerString(genex)))
            ## ~ [1b] ~ PPIPairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppifas = {}     # Hub:{type:[spokeseq]}
            px = 0.0; ptot = len(ppidict)
            for hub in rje.sortKeys(ppidict):
                self.progLog('\r#PPI','Processing PPI: %.2f%%' % (px/ptot)); px += 100.0
                ppi = ppidict[hub] 
                seq = ppi['SpokeSeq']
                if hub not in ppifas:
                    ppifas[hub] = {}
                    for type in ['ppi','bin','com','y2h']:
                        ppifas[hub][type] = []
                        hfile = '../HumSF09_PPI/HumSF09_PPIFas-%s/%s.%s.fas' % (type,hub,type)
                        if os.path.exists(hfile): ppifas[hub][type] = SeqInfoListFromFile(self,hfile,key='short',startfrom=None)
                for type in ['ppi','bin','com','y2h']:
                    if seq in ppifas[hub][type]: ppi[type] = 'Y'
                    elif ppifas[hub][type]: ppi[type] = 'N'
                    else: ppi[type] = '-'
                rje.delimitedFileOutput(self,newpfile,pheaders,datadict=ppi)
            self.printLog('\r#PPI','%s PPI processed.' % (rje.integerString(ptot)))
                
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def SF09x3(self): ### HumSF09 Tidy 3
        '''HumSF09 Tidy.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mfile = '/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/HGNC/genemap_090505.data.tdt'
            mheaders = ['Gene','Entrez','HPRD','OMIM','UniProt','EnsEMBL','EnsLoci','EnsDesc']
            mapdict = rje.dataDict(self,mfile,['Gene'],mheaders)
            mheaders += ['ppi','bin','com','y2h']
            pfile = 'pingu.pairwise.tdt'
            ppidict = rje.dataDict(self,pfile,['Hub','Spoke'],'All',getheaders=True)
            pheaders = ppidict.pop('Headers') + ['ppi','bin','com','y2h']
            genemap = rje.dataDict(self,pfile,['SpokeSeq'],['Spoke'],getheaders=False)
            ## ~ [0a] ~ New Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newmfile = 'humsf09.genemap.0905050.tdt'
            rje.delimitedFileOutput(self,newmfile,mheaders,rje_backup=True)
            newpfile = 'humsf09_ppi.pairwise.tdt'
            rje.delimitedFileOutput(self,newpfile,pheaders,rje_backup=True)
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ PPIPairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppifas = {}     # Hub:{type:[spokeseq]}
            px = 0.0; ptot = len(ppidict); nomap = []; mapped = []; spokemap = []
            for pair in rje.sortKeys(ppidict):
                self.progLog('\r#PPI','Processing PPI: %.2f%%' % (px/ptot)); px += 100.0
                ppi = ppidict[pair]
                hub = ppi['Hub']
                seq = ppi['SpokeSeq']
                if hub not in ppifas:
                    ppifas[hub] = {}
                    for type in ['ppi','bin','com','y2h']:
                        ppifas[hub][type] = []
                        hfile = 'HumSF09_PPIFas-%s/%s.%s.fas' % (type,hub,type)
                        self.deBug(hfile)
                        if os.path.exists(hfile): ppifas[hub][type] = rje_seq.SeqInfoListFromFile(self,hfile,key='short',startfrom=None)
                for type in ['ppi','bin','com','y2h']:
                    if seq in ppifas[hub][type]: ppi[type] = 'Y'
                    elif ppifas[hub][type]: ppi[type] = 'N'
                    else: ppi[type] = '-'
                try:
                    if hub not in mapped:
                        for type in ['ppi','bin','com','y2h']: mapdict[hub][type] = len(ppifas[hub][type])
                        rje.delimitedFileOutput(self,newmfile,mheaders,datadict=mapdict[hub])
                        mapped.append(hub)
                except:
                    if hub not in nomap:
                        nomap.append(hub)
                        if rje.matchExp('(\D)',hub):
                            self.errorLog('Problem with mapping for %s' % hub,quitchoice=False)
                if ppi['Spoke'] not in nomap + mapped + spokemap: spokemap.append(ppi['Spoke'])
                self.deBug(ppi)
                rje.delimitedFileOutput(self,newpfile,pheaders,datadict=ppi)
            for spoke in spokemap:
                if spoke in mapped + nomap: continue
                for type in ['ppi','bin','com','y2h']: mapdict[spoke][type] = '-'
                rje.delimitedFileOutput(self,newmfile,mheaders,datadict=mapdict[spoke])
            self.printLog('\r#PPI','%s PPI processed.' % (rje.integerString(ptot)))
            nfile = 'humsf09.no_map.0905050.tdt'
            nomap.sort()
            open(nfile,'w').write(string.join(nomap,'\n'))
            self.printLog('#NOMAP','%s hubs without mapping -> %s' % (rje.integerString(len(nomap)),nfile))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def SF09x4(self): ### HumSF09 Tidy 4
        '''HumSF09 Tidy.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfile = 'humsf09.pairwise_ppi.090505.tdt'
            ppidict = rje.dataDict(self,pfile,['Hub','Spoke'],'All',getheaders=False)
            self.printLog('#PPI','%s pairwise PPI' % rje.integerString(len(ppidict)))
            hspokes = rje.dataDict(self,pfile,['Hub'],['Spoke'],getheaders=False,lists=True)
            self.printLog('#HUB','%s pairwise PPI Hubs' % rje.integerString(len(hspokes)))
            genemap = rje.dataDict(self,pfile,['SpokeSeq'],['Spoke'],getheaders=False)
            self.printLog('#MAP','%s pairwise genes mapped' % rje.integerString(len(genemap)))
            dfile = '../HumSF09_UniFake/ens_HUMAN.unifake.pfam.tdt'
            dommap = rje.dataDict(self,dfile,['Type'],['Name'],getheaders=False,lists=True)
            self.printLog('#DOM','%s domains mapped' % rje.integerString(len(dommap)))
            ## ~ [0a] ~ New Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newdfile = 'humsf09.domain_ppi.090505.tdt'
            dheaders = ['Domain','Spoke','SpokeUni','SpokeSeq','Hubs','Evidence','ppi','bin','com','y2h']
            rje.delimitedFileOutput(self,newdfile,dheaders,rje_backup=True)
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppifas = {}     # Dom:{type:[spokeseq]}
            px = 0.0; ptot = len(dommap)
            for domain in rje.sortKeys(dommap):
                self.progLog('\r#DOM','Processing Domains: %.2f%%' % (px/ptot)); px += 100.0
                ## Map domains onto hubs ##
                dhubs = []
                for seq in dommap[domain]['Name']:
                    if genemap.has_key(seq): dhubs.append(genemap[seq]['Spoke'])
                if not dhubs: continue
                ## Map domains onto spokes through hubs ##
                dspokes = {}    # spoke:[hubs]
                for hub in dhubs:
                    if hub not in hspokes: continue
                    for spoke in hspokes[hub]['Spoke']:
                        if spoke not in dspokes: dspokes[spoke] = []
                        dspokes[spoke].append(hub)
                ## Output domain-spoke pairs ##
                for spoke in dspokes:
                    sdata = {'Domain':domain,'Hubs':string.join(dspokes[spoke],',')}
                    evidence = []
                    for hub in dspokes[spoke]:
                        pkey = '%s\t%s' % (hub,spoke)
                        if pkey not in ppidict: continue
                        sdata = rje.combineDict(sdata,ppidict[pkey],overwrite=False)
                        for e in string.split(ppidict[pkey]['Evidence'],'|'):
                            if e not in evidence: evidence.append(e)
                    evidence.sort()
                    sdata['Evidence'] = string.join(evidence,'|')
                    try: seq = sdata['SpokeSeq']
                    except:
                        self.errorLog('No sequence for %s' % spoke)
                        continue
                    if domain not in ppifas:
                        ppifas[domain] = {}
                        for type in ['ppi','bin','com','y2h']:
                            ppifas[domain][type] = []
                            hfile = '../HumSF09_PPIFas/HumSF09_PPIFas-%sdom/%s.%sdom.fas' % (type,domain,type)
                            self.deBug(hfile)
                            if os.path.exists(hfile): ppifas[domain][type] = rje_seq.SeqInfoListFromFile(self,hfile,key='short',startfrom=None)
                    for type in ['ppi','bin','com','y2h']:
                        if seq in ppifas[domain][type]: sdata[type] = 'Y'
                        elif ppifas[domain][type]: sdata[type] = 'N'
                        else: sdata[type] = '-'
                    self.deBug(sdata)
                    rje.delimitedFileOutput(self,newdfile,dheaders,datadict=sdata)
            self.printLog('\r#DOM','%s domains processed.' % (rje.integerString(ptot)))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def SF09(self): ### HumSF09 Tidy 5
        '''HumSF09 Tidy.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfile = 'HumSF09_MainData/humsf09.pairwise_ppi.090505.tdt'
            ppidict = rje.dataDict(self,pfile,['Hub','Spoke'],'All',getheaders=True)
            self.printLog('#PPI','%s pairwise PPI' % rje.integerString(len(ppidict)))
            #dfile = 'HumSF09_MainData/humsf09.domain_ppi.090505.tdt'
            #domppi = rje.dataDict(self,dfile,['Domain','Spoke'],'All',getheaders=True)
            #self.printLog('#PPI','%s domain PPI' % rje.integerString(len(domppi)))
            mfile = 'HumSF09_MainData/humsf09.genemap.0905050.tdt'
            genemap = rje.dataDict(self,mfile,['Gene'],'All',getheaders=False)
            self.printLog('#MAP','%s pairwise genes mapped' % rje.integerString(len(genemap)))
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppifas = {}     # Hub:{type:[spokeseq]}
            badkeys = []
            px = 0.0; ptot = len(ppidict) #+ len(domppi)
            for pkey in rje.sortKeys(ppidict):
                self.progLog('\r#PPI','Processing PPI: %.2f%%' % (px/ptot)); px += 100.0
                if pkey == 'Headers': continue
                try: (hub,spoke) = string.split(pkey)
                except:
                    if pkey not in badkeys: badkeys.append(pkey)
                    continue
                seq = ppidict[pkey]['SpokeSeq']
                if hub not in ppifas:
                    ppifas[hub] = {}
                    for type in ['ppi','bin','com','y2h']:
                        ppifas[hub][type] = []
                        hfile = './HumSF09_PPIFas/HumSF09_PPIFas-%s/%s.%s.fas' % (type,hub,type)
                        hfile2 = './HumSF09_PPIFas/HumSF09_PPIFas-%s/%s.%sppi.fas' % (type,hub,type)
                        #hfile = '../HumSF09_PPIFas/HumSF09_PPIFas-%sdom/%s.%sdom.fas' % (type,domain,type)
                        self.deBug(hfile)
                        if os.path.exists(hfile2): os.rename(hfile2,hfile)
                        if os.path.exists(hfile): ppifas[hub][type] = rje_seq.SeqInfoListFromFile(self,hfile,key='short',startfrom=None)
                for type in ['ppi','bin','com','y2h']:
                    if seq in ppifas[hub][type]: ppidict[pkey][type] = 'Y'
                    elif ppifas[hub][type]: ppidict[pkey][type] = 'N'
                    else: ppidict[pkey][type] = '-'
            ## ~ [1b] ~ Domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ppifas = {}     # Hub:{type:[spokeseq]}
            for pkey in []: #rje.sortKeys(domppi):
                break # ! DONE OK !
                self.progLog('\r#PPI','Processing PPI: %.2f%%' % (px/ptot)); px += 100.0
                if pkey == 'Headers': continue
                try: (hub,spoke) = string.split(pkey)
                except:
                    if pkey not in badkeys: badkeys.append(pkey)
                    continue
                seq = domppi[pkey]['SpokeSeq']
                if hub not in ppifas:
                    ppifas[hub] = {}
                    for type in ['ppi','bin','com','y2h']:
                        ppifas[hub][type] = []
                        hfile = './HumSF09_PPIFas/HumSF09_PPIFas-%sdom/%s.%sdom.fas' % (type,hub,type)
                        hfile2 = './HumSF09_PPIFas/HumSF09_PPIFas-%sdom/%s.dom%s.fas' % (type,hub,type)
                        hfile3 = './HumSF09_PPIFas/HumSF09_PPIFas-%sdom/%s.%sdomppi.fas' % (type,hub,type)
                        self.deBug(hfile)
                        if os.path.exists(hfile2): os.rename(hfile2,hfile)
                        if os.path.exists(hfile3): os.rename(hfile3,hfile)
                        if os.path.exists(hfile): ppifas[hub][type] = rje_seq.SeqInfoListFromFile(self,hfile,key='short',startfrom=None)
                for type in ['ppi','bin','com','y2h']:
                    if seq in ppifas[hub][type]: domppi[pkey][type] = 'Y'
                    elif ppifas[hub][type]: domppi[pkey][type] = 'N'
                    else: domppi[pkey][type] = '-'
            self.printLog('\r#PPI','%s genes/domains processed.' % (rje.integerString(ptot)))
            ### ~ [2] Re-output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            phead = ppidict.pop('Headers')
            rje.delimitedFileOutput(self,pfile,phead,rje_backup=True)
            px = 0.0; ptot = len(ppidict)
            for pkey in rje.sortKeys(ppidict):
                self.progLog('\r#OUT','Output PPI: %.2f%%' % (px/ptot)); px += 100.0
                if pkey in badkeys: continue
                rje.delimitedFileOutput(self,pfile,phead,datadict=ppidict[pkey])
            self.printLog('\r#OUT','Processing PPI: %.2f%%' % (px/ptot))
            return # !DONE DOMAINS OK!
            dhead = domppi.pop('Headers')
            rje.delimitedFileOutput(self,dfile,dhead,rje_backup=True)
            px = 0.0; ptot = len(domppi)
            for pkey in rje.sortKeys(domppi):
                self.progLog('\r#OUT','Output Domain PPI: %.2f%%' % (px/ptot)); px += 100.0
                if pkey in badkeys: continue
                rje.delimitedFileOutput(self,dfile,dhead,datadict=domppi[pkey])
            self.printLog('\r#OUT','Output Domain PPI: %.2f%%' % (px/ptot))
            
#domppi            
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: OddJob Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: Misc methods                                                                                           #
#########################################################################################################################
def abNprob(a,b,N,overlap=0):     ### Returns probabilities of overlap in a given b/N and in b given a/N
    '''Returns probabilities of no overlap in a given b/N and in b given a/N.'''
    callobj = OddJob()
    pa = 0.0
    pb = 0.0
    for k in range(overlap+1):
        pa += rje.binomial(k,a,float(b)/N,exact=True,callobj=callobj)
        pb += rje.binomial(k,b,float(a)/N,exact=True,callobj=callobj)
    ea = a * float(b) / N
    eb = b * float(a) / N
    print 'Prob <=%d from %d a in %d b = %.4f (Exp %.3f)' % (overlap,a,b,pb,eb)
    print 'Prob <=%d from %d b in %d a = %.4f (Exp %.3f)' % (overlap,b,a,pa,ea)
#########################################################################################################################
### End of SECTION III: Misc methods                                                                                    #
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
    try: OddJob(log=mainlog,cmd_list=cmd_list).run()
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