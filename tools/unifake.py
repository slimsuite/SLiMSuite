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
Program:      UniFake
Description:  Fake UniProt DAT File Generator
Version:      1.3
Last Edit:    17/04/12
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This program runs a number of in silico predication programs and converts protein sequences into a fake UniProt DAT
    flat file. Additional features may be given as one or more tables, using the features=LIST option. Please see the
    UniFake Manual for more details. 

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence file. See rje_seq documentation for filtering options. [None]
    spcode=X        : Species code to use if it cannot be established from sequence name [None]
    features=LIST   : List of files of addtional features in delimited form []    
    aliases=FILE    : File of aliases to be added to Accession number list (for indexing) [None]
    pfam=FILE       : PFam HMM download [None]
    unipath=PATH    : Path to real UniProt Datafile (will look here for DB Index file made with rje_dbase)
    unireal=LIST    : Real UniProt data to add to UniFake output ['AC','GN','RC','RX','CC','DR','PE','KW']
    fudgeft=T/F     : Fudge the real features left/right until a sequence match is found [True]

    ### ~ PROCESSING/OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    unifake=LIST    : List of predictions to add to entries [tmhmm,signalp,disorder,pfam,uniprot]
    datout=FILE     : Name of output DAT file [Default input FILE.dat]
    disdom=X        : Disorder threshold below which to annotate PFam domain as "DOMAIN" [0.0]
    makeindex=T/F   : Whether to make a uniprot index file following run [False]
    ensdat=T/F      : Look for acc/pep/gene in sequence name [False]
    tmhmm=FILE      : Path to TMHMM program [None]
    cleanup=T/F     : Remove TMHMM files after run [True]
    signalp=FILE    : Path to SignalP program [None]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_tm, rje_uniprot, rje_zen
#import rje_hmm_V1 as rje_hmm
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Version with full functionality.
    # 1.1 - Made changes to temp file names for parallel running
    # 1.2 - Added option for reading in Features and limited other info from real UniProt entry
    # 1.3 - Removed searching in UniProt if not UniReal features given.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Expand automatic species identification.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('UNIFAKE', '1.3', 'April 2012', '2008')
    description = 'Fake UniProt DAT File Generator'
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
            if rje.yesNo('Show rje_disorder commandline options?'): out.verbose(-1,4,text=rje_disorder.__doc__)
            if rje.yesNo('Show rje_seq commandline options?'): out.verbose(-1,4,text=rje_seq.__doc__)
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
### SECTION II: UniFake Class                                                                                           #
#########################################################################################################################
class UniFake(rje.RJE_Object):     
    '''
    UniFake DAT Generator Class. Author: Rich Edwards (2008).

    Info:str
    - Aliases = File of aliases to be added to Accession number list (for indexing) [None]
    - DatOut = Name of output DAT file [Default input FILE.dat]
    - PFam = PFam HMM download [None]
    - SPCode = Species code to use if it cannot be established from sequence name [None]
    - TMHMM = Path to TMHMM program [None]
    - SignalP = Path to SignalP program [None]
    
    Opt:boolean
    - CleanUp = Remove TMHMM directories
    - EnsDat = Look for acc/pep/gene in sequence name [False]
    - FudgeFT = Fudge the real features left/right until a sequence match is found [True]
    - MakeIndex = Whether to make a uniprot index file following run [False]
    
    Stat:numeric
    - DisDom = Disorder threshold below which to annotate PFam domain as "DOMAIN" [0.0]

    List:list
    - Features = Files of addtional features in delimited form [None]    
    - UniFake = List of predictions to add to DAT file [tmhmm,signalp,disorder,pfam]
    - UniReal = Real UniProt data to add to UniFake output ['AC','GN','RC','RX','CC','DR','PE','KW']

    Dict:dictionary    
    - Aliases = Aliases to be added to Accession number list {ID:[Aliases]}
    - Features = Addtional features for given sequence ID {ID:[{Feature,Start,End,Description}]}

    Obj:RJE_Objects
    - SeqList = Main sequence list object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Aliases','DatOut','PFam','SPCode','SignalP','TMHMM']
        self.optlist = ['CleanUp','EnsDat','MakeIndex','FudgeFT']
        self.statlist = ['DisDom']
        self.listlist = ['Features','UniFake']
        self.dictlist = ['Aliases','Features','UniReal']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
        self.setOpt({'CleanUp':True})
        self.list['UniFake'] = ['tmhmm','signalp','disorder','pfam','uniprot']
        self.list['UniReal'] = ['AC','GN','RC','RX','CC','DR','PE','KW']
        self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','gnspacc=T','datout=F'])
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
                self._cmdReadList(cmd,'file',['Aliases','DatOut','PFam','SignalP','TMHMM'])
                self._cmdReadList(cmd,'info',['SPCode'])
                self._cmdReadList(cmd,'opt',['CleanUp','EnsDat','FudgeFT','MakeIndex'])
                self._cmdReadList(cmd,'stat',['DisDom'])
                self._cmdReadList(cmd,'list',['UniFake','UniReal'])
                self._cmdReadList(cmd,'glist',['Features'])  
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def run(self):  ### Performs main run method, including both setup and UniFake
        '''Performs main run method, including both setup and UniFake.'''
        ### ~ [1] ~ Setup aliases and features dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.setup()
        ### ~ [2] ~ Perform main UniFake file generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.uniFake()
        ### ~ [3] ~ Index files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['MakeIndex']:
            i = self.stat['Interactive']
            self.stat['Interactive'] = -1
            self.info['UniPath'] = rje.makePath(os.path.split(self.info['DatOut'])[0])
            rje_uniprot.processUniProt(self,makeindex=True,makespec=False,makefas=False)
            self.stat['Interactive'] = i
#########################################################################################################################
    def setup(self,clear=True):    ### Sets up the Aliases and Features dictionaries
        '''Sets up the Aliases and Features dictionaries.'''
        try:### ~ [1] ~ Setup Alias dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear: self.dict['Aliases'] = {}
            self.loadAlias(self.info['Aliases'])
            ### ~ [2] ~ Setup Features dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear: self.dict['Features'] = {}
            for ftfile in self.list['Features']: self.loadFeatures(ftfile)
        except: self.errorLog('Disaster during UniFake.setup(). Aliases/Features may not be loaded')
#########################################################################################################################
    def uniFake(self,seqs=[],store=False):  ### Main UniFake method. Runs on sequences in self.obj['SeqList'] if no seqs.
        '''Main UniFake method. Runs on sequences in self.obj['SeqList'] if no seqs given.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            unifake = string.split(string.join(self.list['UniFake']).lower())
            seqlist = self.obj['SeqList']
            if seqs: seqlist.seq = seqs
            else: seqs = seqlist.seq
            (sx,seqnum) = (0,seqlist.seqNum())
            ## ~ [1b] Setup UniProt object and output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)   # UniProt object for saving data
            if self.info['DatOut'].lower() in ['','none']: self.info['DatOut'] = rje.baseFile(seqlist.info['Name']) + '.dat'
            datfile = self.info['DatOut']
            if os.path.exists(datfile): rje.backup(self,datfile)
            if store: seqlist.obj['UniProt'] = uniprot
            ## ~ [1c] Setup RJE_HMM object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'pfam' in unifake:
                hmm = rje_hmm.HMMRun(self.log,self.cmd_list+['force=T'])
                hmmfile = '%s.pfam.tdt' % rje.baseFile(datfile)
                if os.path.exists(hmmfile): rje.backup(self,hmmfile)
                hmm.list['HMM'] = [self.info['PFam']]
                hmm.opt['HMMPFam'] = True
            else: hmm = None
            ## ~ [1d] Setup RJE_TM object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'signalp' in unifake: tm = rje_tm.TM(self.log,self.cmd_list)
            else: tm = None
            ### ~ [2] ~ Perform UniFake processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in seqs:
                sx += 1
                name = seq.shortName()                    
                self.printLog('#SEQ','Processing %s (%s aa) %s...' % (seq.shortName(),rje.integerString(seq.aaLen()),seq.info['Description'][:50]))
                try:
                    ## ~ [2a] ~ Basic data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    utmp = 'tmp%s.%s' % (rje.randomString(5),seq.info['AccNum'])
                    open('%s.fas' % utmp,'w').write('>%s\n%s\n' % (seq.shortName(),seq.info['Sequence']))
                    udata = {'CC':['-!- Features generated using unifake.py'],'AC':[]}
                    if seq.info['SpecCode'] in ['Unknown','UNK']: seq.info['SpecCode'] = self.info['SPCode']
                    #x#elif seq.info['Species'] != 'None': udata['OS'] = [seq.info['Species']]     #!# Check how well this works. Add spectable? #!#
                    ## ~ [2b] ~ Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if self.opt['EnsDat'] and rje.matchExp('\[acc:(\S+) pep:(\S+) gene:(\S+)\]',seq.info['Name']):
                        details = rje.matchExp('\[acc:(\S+) pep:(\S+) gene:(\S+)\]',seq.info['Name'])
                        self.addAlias(seq.info['AccNum'],details[0])
                        self.addAlias(seq.info['AccNum'],details[1])
                        self.addAlias(seq.info['AccNum'],details[2])
                        udata['GN'] = [details[2]]
                    for id in [seq.shortName(),seq.info['AccNum']]:
                        if id in self.dict['Aliases']: udata['AC'].append('%s;' % string.join(self.dict['Aliases'][id],'; '))
                    ## ~ [2c] ~ Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    ft = []     # List of features for sequence
                    for id in [seq.shortName(),seq.info['AccNum'],seq.info['ID']]:
                        if id in self.dict['Features']: ft += self.dict['Features'][id]                        
                    ## ~ [2d] IUPRED disorder prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if 'disorder' in self.list['UniFake']:
                        try:
                            seq.disorder()
                            dis = seq.obj['Disorder']
                            for disorder in seq.obj['Disorder'].list['RegionDisorder']:
                                ft.append({'Type':'DISORDER','Desc':'Predicted disorder: %s' % seq.obj['Disorder'].info['Disorder'],'Start':disorder[0],'End':disorder[1]})
                                if dis.info['Disorder'].lower() == 'iupred': ft[-1]['Desc'] = '%s > %.2f' % (ft[-1]['Desc'],dis.stat['IUCut'])
                            for fold in seq.obj['Disorder'].list['RegionFold']:
                                ft.append({'Type':'ORDER','Desc':'Predicted order: %s' % seq.obj['Disorder'].info['Disorder'],'Start':fold[0],'End':fold[1]})
                                if dis.info['Disorder'].lower() == 'iupred': ft[-1]['Desc'] = '%s <= %.2f' % (ft[-1]['Desc'],dis.stat['IUCut'])
                        except: self.log.errorLog('UniFake disorder problem for %s.' % name)
                    ## ~ [2e] PFam HMM domain prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if hmm:
                        try:
                            hmm.setInfo({'SearchDB':'%s.fas' % utmp,'HMMOut':'%s.hmm.out' % utmp})      # This will be made for each sequence                    
                            hmm.search = []
                            hmm.list['HMMRes'] = [hmm.hmmSearch(self.info['PFam'],outfile=hmm.info['HMMOut'])]   # Used in hmmTable
                            hmm.hmmTable(outfile=hmmfile,append=True)
                            if 'disorder' in self.list['UniFake']: disorder = seq.obj['Disorder'].list['ResidueDisorder']          # individual (IUPRed) residue results
                            else: disorder = []
                            if hmm.search: udata['CC'].append('PFam: HMMer PFam search vs %s (Modified %s)' % (self.info['PFam'],time.ctime(os.path.getmtime(self.info['PFam']))))
                            else:
                                udata['CC'].append('-!- ERROR: PFam HMMer Search failure!')
                                out = {'Type':'!ERROR!','Name':name}
                                rje.delimitedFileOutput(self,hmmfile,['Type','Name','Start','End','Eval','Score'],datadict=out)
                            for search in hmm.search:
                                for hit in search.hit:
                                    for aln in hit.aln:
                                        pfamft = {'Start':aln.stat['SbjStart'],'End':aln.stat['SbjEnd'],'Type':'PFAM',
                                                   'Desc':'%s PFam HMM Eval: %.2e; Score: %.1f' % (search.info['Name'],aln.stat['Expect'],aln.stat['BitScore'])}
                                        if disorder:
                                            region = disorder[aln.stat['SbjStart']-1:aln.stat['SbjEnd']]
                                            hmmdisorder = float(sum(region)) / len(region)
                                            pfamft['Desc'] = '%s; IUPRed: %.2f' % (pfamft['Desc'],hmmdisorder)
                                            if hmmdisorder < self.stat['DisDom']: pfamft['Type'] = 'DOMAIN'
                                        ft.append(pfamft)
                        except: self.log.errorLog('UniFake PFam HMM problem for %s.' % name)                  
                    ## ~ [2f] TMHMM transmembrane topology prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if 'tmhmm' in unifake:
                        try:
                            tmdat = os.popen('%s %s.fas -short' % (self.info['TMHMM'],utmp)).readlines()
                            domlist = rje_tm.domainList(rje_tm.parseTMHMM(tmdat[0]))
                            for tmdom in domlist:
                                ft.append(tmdom)
                                ft[-1]['Desc'] = 'TMHMM topology prediction'
                                ft[-1]['Start'] = string.atoi(ft[-1]['Start'])
                                ft[-1]['End'] = string.atoi(ft[-1]['End'])
                            if len(domlist) > 1: udata['CC'].append('TMHMM: %d TM domains; N-Term %s' % ((len(domlist)-1)/2,domlist[0]['Type']))
                            else: udata['CC'].append('TMHMM: 0 TM domains')
                        except: self.log.errorLog('UniFake TMHMM problem for %s.' % name)
                    ## ~ [2g] SIGNALP signal peptide prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if 'signalp' in unifake:
                        try:
                            os.system('%s -f short -t euk %s.fas > %s.signalp' % (self.info['SignalP'],utmp,utmp))
                            tm.signalp = {}
                            tm.parseSignalP('%s.signalp' % utmp)
                            sigp = tm.signalp.pop(seq.shortName())
                            cpos = 0
                            if sigp['nn_ymax?'] == 'Y':
                                cpos = string.atoi(sigp['nn_ymaxpos'])
                                desc = 'SignalP NN prediction'
                            if sigp['hmm_cmax?'] == 'Y':
                                hmm_c = string.atoi(sigp['hmm_cmaxpos'])
                                if cpos == 0:
                                    cpos = hmm_c
                                    desc = 'SignalP HMM prediction'
                                else:
                                    if hmm_c < cpos:
                                        cpos = hmm_c
                                        desc = 'SignalP HMM prediction (NN also Y)'
                                    else: desc += ' (HMM also Y)'
                            if cpos > 0: ft.append({'Type':'SIGNALP','Desc':desc,'Start':1,'End':cpos})
                        except: self.log.errorLog('UniFake SignalP problem for %s.' % name)
                    ## ~ [2h] Convert to UniProt and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.addRealUniProt(seq,udata,ft)
                    self.deBug(ft)
                    if not store: uniprot.list['Entry'] = []
                    if uniprot.addFromSeq(seq,data=udata,ft=ft):    ### Converts into UniProtEntry object 
                        if not store: uniprot.saveUniProt(datfile,append=True)
                        #x#open(self.info['DatPickup'],'a').write('%s\n' % seq.shortName())
                ## ~ [2f] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                except: self.log.errorLog('Problem during UniFake(%s)' % name)
                for tmp in glob.glob('%s*' % utmp): os.unlink(tmp)
                self.printLog('#UNIFAKE','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(sx),rje.integerString(seqnum-sx)),log=False)
            if store: uniprot.saveUniProt(datfile,append=False)
            if self.opt['CleanUp']:
                for tmp in glob.glob('TMHMM*'):
                    if os.path.isdir(tmp): os.rmdir(tmp)            
        except: self.errorLog('Oh, the shame of it! Trouble during UniFake.uniFake()')
#########################################################################################################################
    ### <3> ### Alias Methods                                                                                           #
#########################################################################################################################
    def addAlias(self,id,alias):  ### Add id to self.dict['Alias'][alias]
        '''Add alias to self.dict['Aliases'][id].'''
        try:
            badalias = ['!FAILED!','','NONE','-','N/A']     # Entries to most definitely ignore!
            if alias.upper() in badalias or id.upper() in badalias or id == alias: return
            if id not in self.dict['Aliases']: self.dict['Aliases'][id] = []
            if alias not in self.dict['Aliases'][id]: self.dict['Aliases'][id].append(alias)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadAlias(self,sourcefile):  ### Loads Alias data
        '''
        Loads Alias data.
        >> sourcefile:str = Source filename
        '''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if sourcefile.lower() in ['','none']: return 
            if not os.path.exists(sourcefile): return self.log.errorLog('Alias file "%s" not found' % (sourcefile),printerror=False)
            data = rje.dataDict(self,sourcefile,datakeys=['Aliases'],lists=True)
            ### ~ [2] Parse out Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (hx,htot) = (0.0,len(data))
            for id in data:
                self.log.printLog('\r#ALIAS','Processing %s: %.1f%%' % (sourcefile,hx/htot),newline=False,log=False)
                hx += 100.0
                ## ~ [2a] Update self.dict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for alist in data[id]['Aliases']:
                    for alias in string.split(alist,','): self.addAlias(id,alias)
                if id in self.dict['Aliases']: self.dict['Aliases'][id].sort()
            self.log.printLog('\r#ALIAS','Processed %s: %s IDs with aliases' % (sourcefile,rje.integerString(len(self.dict['Aliases']))))           
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### Add Features Methods                                                                                    #
#########################################################################################################################
    def loadFeatures(self,ftfile):  ### Loads features from given file
        '''Loads features from given file.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ftfile in ['','none']: return
            if not os.path.exists(ftfile): return self.printLog('#ERR','Features file "%s" missing')
            delimit = rje.delimitFromExt(filename=ftfile)
            ## ~ [1a] ~ Establish headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headers = rje.readDelimit(open(ftfile,'r').readline(),delimit)
            mainkeys = [headers[0]]
            hmap = {}
            for h in headers: hmap[h.lower()] = h
            pos = ''    # Leader for start/end positions
            if 'ft_start' in hmap or 'ft_end' in hmap: pos = 'ft_'
            for h in ['feature','%sstart' % pos,'%send' % pos,'description']:
                if h not in hmap: return self.printLog('#ERR','No %s field detected in "%s" features file' % (h,ftfile))
                mainkeys.append(hmap[h])
            mainkeys.remove(hmap['description'])
            ### ~ [2] ~ Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ftdata = rje.dataDict(self,ftfile,mainkeys,['description'],delimit,headers,lists=True)
            (mx,mtot,fx) = (0.0,len(ftdata),0)
            for mainkey in rje.sortKeys(ftdata):
                self.progLog('\r#FT','Loading features from %s: %.2f%%' % (ftfile,mx/mtot))
                mx += 100.0                                                                           
                (id,ft,start,end) = string.split(mainkey,delimit)
                if id == mainkeys[0]: continue
                if id not in self.dict['Features']: self.dict['Features'][id] = []
                for desc in ftdata[mainkey][hmap['description']]:
                    fx += 1
                    self.dict['Features'][id].append({'Type':ft,'Start':int(start),'End':int(end),'Desc':desc})
            self.printLog('\r#FT','Loaded %s features for %s IDs from %s' % (rje.integerString(fx),rje.integerString(len(self.dict['Features'])),ftfile))
        except: self.errorLog('UniFake.loadFeatures error ["%s"]' % ftfile)
#########################################################################################################################
    def addRealUniProt(self,seq,udata,ftlist):  ### Updates feature list ft using real UniProt where possible and makes NR
        '''
        Updates feature list ft using real UniProt where possible and makes NR.
        >> seq:Sequence object = target of UniFake
        >> udata:UniProt Entry Data dictionary *Modified in place*
        >> ftlist:list of feature dictionaries to add to (and make NR) *Modified in place*
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['UniReal']: return
            realuni = rje_uniprot.UniProt(self.log,self.cmd_list+['datout=None'])
            realuni.readUniProt(clear=True,acclist=[seq.shortName(),seq.info['ID'],seq.info['AccNum']],cleardata=False)
            if not realuni.list['Entry']: self.printLog('#UNI','No Real AccNum for %s' % seq.shortName())
            sequence = seq.info['Sequence'][0:]
            ### ~ [1] ~ Map and Add Features from actual UniProt entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in realuni.list['Entry']:
                if 'uniprot' not in self.list['UniFake']: break
                for key in self.list['UniReal']:    #['AC','GN','RC','RX','CC','DR','PE','KW']:
                    if entry.dict['Data'].has_key(key):
                        if udata.has_key(key): udata[key] = entry.dict['Data'][key] + udata[key]
                        else: udata[key] = entry.dict['Data'][key][0:]
                for ft in entry.list['Feature'][0:]:
                    ft_start = ft['Start']
                    ft_end = ft['End']
                    ft_seq = entry.obj['Sequence'].info['Sequence'][ft_start-1:ft_end]
                    if ft_seq == sequence[ft_start-1:ft_end]: ftlist.append(ft); continue
                    if not self.opt['FudgeFT']: continue
                    fudge = 1
                    while fudge:
                        if ft_start - fudge < 1 and ft_end + fudge > len(sequence): fudge = 0; break
                        if ft_start - fudge >= 1 and ft_seq == sequence[ft_start-1-fudge:ft_end-fudge]: fudge = -fudge; break
                        if ft_end + fudge <= len(sequence) and ft_seq == sequence[ft_start-1+fudge:ft_end+fudge]: break
                        fudge += 1
                    if fudge:
                        ft['Start'] = ft_start + fudge
                        ft['End'] = ft_end + fudge
                        ftlist.append(ft)
            ### ~ [2] ~ Make FT list NR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            i = 0
            while i < len(ftlist):
                if ftlist.count(ftlist[i]) > 1: ftlist.pop(i)
                else: i += 1
        except: self.errorLog('UniFake.addRealUniProt error [%s]' % seq.shortName())
#########################################################################################################################
### End of SECTION II: UniFake Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
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
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: UniFake(mainlog,cmd_list).run()
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
