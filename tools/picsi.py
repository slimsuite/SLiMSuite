#!/usr/bin/python

# See below for name and description
# Copyright (C) 2010 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      PICSI
Description:  Proteomics Identification from Cross-Species Inference
Version:      1.2
Last Edit:    26/03/14
Citation:     Jones, Edwards et al. (2011), Marine Biotechnology 13(3): 496-504. [PMID: 20924652]
Copyright (C) 2010  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for cross-species protein identifications using searches against NCBInr, for example. MASCOT
    processing uses BUDAPEST. Hits are then converted into peptides for redundancy removal. Hits from a given (known)
    query species are preferentially kept and any peptides belonging to those hits are purged from hits returned by
    other species. All hits are then classified:
    - UNIQUE    : Contains 2+ peptides, including 1+ unique peptides
    - NR        : Contains 2+ peptides; None unique but 1+ peptides only found in other "NR" proteins
    - REDUNDANT : Contains 2+ peptides but all found in proteins also containing UNIQUE peptides
    - REJECT    : Identified by <2 peptides once query-species peptides subtracted

Commandline:
    ### ~~~ INPUT OPTIONS ~~~ ###
    seqin=FILE      : File containing sequences hit during searches [None]
    resfiles=FILES  : List of CSV MASCOT output files []
    sumfile=FILE    : Delimited summary file containing search,prot_hit_num,prot_acc,prot_desc,pep_seq
    qryspec=FILE    : Species code for Query species [EMIHU]
    spectdt=FILE    : File of UniProt species translations ['/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/UniProt/uniprot.spec.tdt']
    ### ~~~ OUTPUT OPTIONS ~~~ ###
    basefile=FILE   : Base for output files [picsi]
    delimit=X       : Delimiter for output files [tab]

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import budapest
import rje, rje_seq, rje_zen
import rje_blast_V2 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Basic working version.
    # 1.1 - Updated to blast_V2 and BLAST+.
    # 1.2 - Updated to BUDAPEST 2.3 and rje_mascot.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add documentation!
    # [ ] : Update to rje_seqlist.
    # [ ] : Add minpep setting and check protein/peptide classification.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('PICSI', '1.2', 'March 2014', '2010')
    description = 'Proteomics Identification from Cross-Species Inference'
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
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
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PICSI Class                                                                                             #
#########################################################################################################################
class PICSI(rje.RJE_Object):     
    '''
    Main PICSI Class. Author: Rich Edwards (2010).

    Info:str
    - QrySpec = Species code for Query species [EMIHU]
    - SeqIn = File containing sequences hit during searches [None]
    - SpecTDT = File of UniProt species translations ['/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/UniProt/uniprot.spec.tdt']
    - SumFile = Delimited summary file containing search,prot_hit_num,prot_acc,prot_desc,pep_seq
    
    Opt:boolean

    Stat:numeric

    List:list
    - ResFiles = List of CSV MASCOT output files []
    
    Dict:dictionary
    - Acc2Seq = Dictionary of {Acc:Sequence Object}
    - MapGI = {gi:Sequence Object}
    - Searches = BUDAPEST data dictionaries with MASCOT data read and stored
    
    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['QrySpec','SeqIn','SpecTDT','SumFile']
        self.optlist = []
        self.statlist = []
        self.listlist = ['ResFiles']
        self.dictlist = ['Acc2Seq','Searches']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'SumFile':'picsi.tdt','QrySpec':'EMIHU',
                      'SpecTDT':'/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/UniProt/uniprot.spec.tdt'})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'file',['SeqIn','SpecTDT','SumFile'])
                self._cmdReadList(cmd,'info',['QrySpec'])
                self._cmdReadList(cmd,'glist',['ResFiles'])  
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.picsi()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method. Makes sumfile if necessary.
        '''Main class setup method. Makes sumfile if necessary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.debug(self.getStrLC('SumFile')); self.debug(self.getStr('SumFile'))
            if self.getStrLC('Basefile') in ['','none']: self.baseFile(rje.baseFile(self.info['SumFile']))
            if self.getStrLC('SumFile') in ['','none']: self.info['SumFile'] = '%s.tdt' % self.basefile()
            self.printLog('#SUM','Summary file: %s' % self.getStr('SumFile'))
            if os.path.exists(self.info['SumFile']) and not self.opt['Force']:
                if rje.yesNo('%s found. Use these results?' % self.info['SumFile']):
                    return self.printLog('#SUM','Summary results file found. No MASCOT processing.')
            mapgi = False
            ### ~ [2] Process MASCOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mfile in self.list['ResFiles']:
                bud = budapest.Budapest(self.log,self.cmd_list+['mascot=%s' % mfile])
                bud.info['Name'] = mfile
                bud.readMascot()
                self.dict['Searches'][mfile] = bud.dict['Hits']
                protacclist = rje.sortKeys(bud.dict['Hits'])
                for protacc in protacclist:
                    if rje.matchExp('gi\|(\d+)',protacc): mapgi = True
                accfile = '%s.%s.protacc' % (self.baseFile(),rje.baseFile(mfile))
                self.debug(accfile)
                open(accfile,'w').write(string.join(protacclist,'\n'))
                self.printLog('#MFILE','%s: %s proteins.' % (mfile,rje.iLen(protacclist)))
            ## ~ [2a] gi Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if mapgi:
            #    mapgi = self.dict['MapGI'] = seqlist.seqNameDic('NCBI')
            #    open('mapgi.tmp','w').write(string.join(rje.sortKeys(mapgi),'\n'))
            ### ~ [3] Setup seqlist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seq.SeqList(self.log,['gnspacc=T']+self.cmd_list)
            self.dict['Acc2Seq'] = seqlist.seqNameDic('Max')
            ### ~ [4] Generate Summary File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sumhead = string.split('search,prot_hit_num,prot_acc,prot_desc,pep_seq',',')
            rje.delimitedFileOutput(self,self.info['SumFile'],sumhead,rje_backup=True)
            for mfile in rje.sortKeys(self.dict['Searches']):
                bud = self.dict['Searches'][mfile]
                for protacc in rje.sortKeys(bud)[0:]:
                    protname = bud[protacc]['prot_acc']
                    protdesc = bud[protacc]['prot_desc']
                    if rje.matchExp('gi\|(\d+)',protacc):
                        gi = rje.matchExp('gi\|(\d+)',protacc)[0]
                        try:
                            protname = self.dict['Acc2Seq'][gi].shortName()
                            protdesc = self.dict['Acc2Seq'][gi].info['Description']
                        except: protname = 'gi_UNK__%s' % gi
                    #x#print protname, protdesc, bud[protacc]
                    for pep in bud[protacc]['Peptides']:
                        data = {'search':rje.baseFile(mfile,True),'prot_desc':protdesc,'prot_acc':protname,
                                'pep_seq':pep,'prot_hit_num':bud[protacc]['prot_hit_num']}
                        rje.delimitedFileOutput(self,self.info['SumFile'],sumhead,datadict=data)
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _method(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def picsi(self):    ### Cleans up cross-species search results
        '''Cleans up cross-species search results.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            datafile = self.info['SumFile']
            delimit = rje.delimitFromExt(filename=self.info['SumFile'])
            data = {}       # search:{hit:{???}}
            pep2prot = {}   # search:{peptide:[hits]}
            id2prot = {}    # search:{id:hit}
            prot2desc = {}
            fullpeplist = {}    
            pepcon = {}     # Convert pep:longer pep
            speclist = []   # List of species codes
            ### ~ [1] Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            indata = rje.dataDict(self,datafile,['search','prot_hit_num'],'All',lists=True)
            for ikey in rje.sortKeys(indata):
                (search,id) = string.split(ikey,delimit)
                prot = indata[ikey]['prot_acc'][0]
                desc = string.replace(indata[ikey]['prot_desc'][0],'Full=','')
                if desc[3:7] == 'Name': desc = desc[9:]
                prot2desc[prot] = desc; self.printLog('#DESC','%s = %s' % (prot,desc))
                indata[ikey]['pep_seq'] = string.join(indata[ikey]['pep_seq'],'|')
                pepconv = string.replace(indata[ikey]['pep_seq'],'I','L')
                pepconv = string.replace(pepconv,'Q','K')
                peplist = rje.sortUnique(string.split(pepconv,'|'))
                indata[ikey]['pep_seq'] = string.join(rje.sortUnique(string.split(indata[ikey]['pep_seq'],'|')),'|')
                if search not in data:
                    data[search] = {}
                    pep2prot[search] = {}
                    id2prot[search] = {}
                    fullpeplist[search] = []
                    pepcon[search] = {}
                fullpeplist[search] += peplist
                id2prot[search][id] = prot
                spec = string.split(prot,'_')[1]
                if spec not in speclist: speclist.append(spec)
                data[search][prot] = {'search':search,'pepcount':len(peplist),'hit':id,'desc':desc,'spec':spec,
                                      'pep_uniq':0,'peplist':indata[ikey]['pep_seq'],'conpep':peplist[0:],
                                      'pep_rem':0}
                try: data[search][prot]['accnum'] = self.dict['Acc2Seq'][prot].info['AccNum']
                except: data[search][prot]['accnum'] = string.split(prot,'__')[-1]
                for pep in peplist:
                    if pep not in pep2prot[search]:
                        pep2prot[search][pep] = []
                    pep2prot[search][pep].append(prot)
            ## ~ [1a] Convert peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for search in fullpeplist:
                fullpeplist[search] = rje.sortUnique(fullpeplist[search])
                for pep in fullpeplist[search][0:]:
                    for pep2 in fullpeplist[search]:
                        if pep != pep2 and pep in pep2:
                            pepcon[search][pep] = pep2
                            fullpeplist[search].remove(pep)
                            break
                for pep in pepcon[search]:
                    while pepcon[search][pep] in pepcon[search]: pepcon[search][pep] = pepcon[search][pepcon[pep]]
                self.printLog('#PEP','%s %s peptide conversions' % (len(pepcon[search]),search))
                #self.deBug(pepcon[search])
                #self.deBug(rje.sortKeys(pep2prot[search]))
                pp = 0; pm = 0
                for prot in data[search]:
                    for pep in data[search][prot]['conpep'][0:]:
                        if pep in pepcon[search]:
                            newpep = pepcon[search][pep]
                            if newpep not in data[search][prot]['conpep']: data[search][prot]['conpep'].append(newpep); pp += 1
                            data[search][prot]['conpep'].remove(pep); pm += 0
                            if prot not in pep2prot[search][newpep]: pep2prot[search][newpep].append(prot)
                            if pep in pep2prot[search]: pep2prot[search].pop(pep)
                    data[search][prot]['pep_con'] = len(data[search][prot]['conpep'])
                self.printLog('#PEP','%s %s converted peptides added; %s removed' % (pp,search,pm))
            ### ~ [2] Calculate Unique/Redundancy status ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for search in pep2prot:
            ## ~ [2a] Species Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                remx = 0
                for prot in data[search]:
                    if data[search][prot]['spec'] != self.info['QrySpec']: continue
                    for pep in data[search][prot]['conpep']:
                        for prot2 in pep2prot[search][pep][0:]:
                            if data[search][prot2]['spec'] == self.info['QrySpec']: continue
                            pep2prot[search][pep].remove(prot2)
                            data[search][prot2]['conpep'].remove(pep)
                            data[search][prot2]['pep_rem'] += 1; remx += 1
                self.printLog('#REM','%s %s peptides removed from non-%s hits' % (rje.integerString(remx),search,self.info['QrySpec']))
            ## ~ [2b] One-hit wonders ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for prot in data[search]:
                    if len(data[search][prot]['conpep']) < 2:
                        for pep in data[search][prot]['conpep']:
                            #if pep in pep2prot[search] and prot in pep2prot[search][pep]:
                            pep2prot[search][pep].remove(prot)
            ## ~ [2c] Unique peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ux = 0
                for pep in pep2prot[search]:
                    #self.deBug(pep)
                    if len(pep2prot[search][pep]) == 1: data[search][pep2prot[search][pep][0]]['pep_uniq'] += 1; ux += 1
                self.printLog('#UNIQ','%s unique %s peptides' % (rje.integerString(ux),search))
            ## ~ [2d] Total Redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                summary = {'HITS':len(data[search]),'REJECT':0,'UNIQUE':0,'NR':0,'REDUNDANT':0}
                rx = 0
                for prot in data[search]:
                    #if data[search][prot]['unique']: data[search][prot]['red'] = False; continue
                    data[search][prot]['pep_red'] = 0   # Redundant peptides found in proteins with unique peptides
                    data[search][prot]['pep_nr'] = 0    # Redundant peptides found only in proteins without unique peptides
                    for pep in data[search][prot]['conpep']:
                        if pep2prot[search][pep] == [prot]: continue
                        upep = False
                        for prot2 in pep2prot[search][pep]:
                            if data[search][prot2]['pep_uniq']: upep = True; break
                        if upep: data[search][prot]['pep_red'] += 1     # Redundant peptide found in unique protein
                        else: data[search][prot]['pep_nr'] += 1         # Redundant peptide NOT found in unique protein
                    if len(data[search][prot]['conpep']) < 2: data[search][prot]['class'] = 'REJECT'; rx += 1
                    elif data[search][prot]['pep_uniq']: data[search][prot]['class'] = 'UNIQUE'
                    elif data[search][prot]['pep_nr']: data[search][prot]['class'] = 'NR'
                    else: data[search][prot]['class'] = 'REDUNDANT'; rx += 1
                    summary[data[search][prot]['class']] += 1
                self.printLog('#REJ','%s rejected %s hits' % (rje.integerString(rx),search))
                for x in rje.sortKeys(summary): self.printLog('#%s' % search,'%s %s' % (summary[x],x))

            ### ~ [3] Species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            speclist.sort()
            species = {}
            for spec in speclist:
                try:
                    grep = os.popen('grep %s %s' % (spec,self.info['SpecTDT'])).read()
                    species[spec] = string.split(grep,':')[-4]
                    self.printLog('#SPEC','%s = %s' % (spec,species[spec]))
                except: species[spec] = '?'

            ### ~ [END] Output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.clean.tdt' % rje.baseFile(self.info['SumFile'])
            headers = ['search','hit','class','accnum','spec','species','desc','pepcount','pep_con','pep_rem','pep_uniq','pep_nr','pep_red','peplist','conpep']
            if self.dict['Acc2Seq']: headers.insert(3,'cluster')
            rje.delimitedFileOutput(self,outfile,headers,datadict={},rje_backup=True)
            for search in rje.sortKeys(data):
                if self.dict['Acc2Seq']: self.clusterGoodSeq(search,data[search])
                for prot in rje.sortKeys(data[search]):
                    if rje.matchExp('^gi:(\d+).+\[(\S.+\S)\]$',data[search][prot]['desc']):
                        data[search][prot]['species'] = rje.matchExp('^gi:(\d+).+\[(\S.+\S)\]$',data[search][prot]['desc'])[1]
                    else: data[search][prot]['species'] = species[data[search][prot]['spec']]                                                                               
                    rje.delimitedFileOutput(self,outfile,headers,datadict=data[search][prot])
                                
        except: self.errorLog('Errg')
#########################################################################################################################
    #!# Modify this BUDAPEST clustering algorithm? #!#
    def clusterGoodSeq(self,searchset,data):   ### Clusters good sequences returned by search and updates data dictionary
        '''Clusters good sequences returned by search and updates data dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Extract Non-rejected sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = rje_seq.SeqList(self.log,['gnspacc=T']+self.cmd_list+['autoload=F'])
            #self.deBug(rje.sortKeys(self.dict['Acc2Seq']))
            for prot in rje.sortKeys(data):
                if data[prot]['class'] != 'REJECT': seqlist.seq.append(self.dict['Acc2Seq'][data[prot]['accnum']])
            if not seqlist.seqNum():
                return self.printLog('#NULL','No %s sequences remain for clustering' % searchset)
            seqfile = '%s.%s.tmpdb' % (self.info['Basefile'],searchset)
            seqlist.saveFasta(seqfile=seqfile)
            seqdict = seqlist.seqNameDic()

            ### ~ [2] Perform BLAST and generate hit matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try:
                blast = rje_blast.blastObj(self.log,['blastf=T','blaste=1e-4']+self.cmd_list+['dna=F'],type='New')
                clusters = blast.blastClusters(seqfile,seqdict=seqdict,keepblast=False) 
            except:
                self.errorLog('Problem with new BLAST clustering')
                blast = rje_blast.blastObj(self.log,['blastf=T','blaste=1e-4']+self.cmd_list+['dna=F'],type='Old')
                blast.setInfo({'InFile':seqfile,'DBase':seqfile,'Name':'%s.tmp.blast' % self.info['Basefile'],'Type':'blastp'})
                blast.setStat({'OneLine':seqlist.seqNum(),'HitAln':0})
                blast.formatDB(fasfile=seqfile,force=True,protein=True)
                blast.blast(cleandb=False,use_existing=False,log=True)
                blast.readBLAST(gablam=False,unlink=True,log=True)
                rje_blast.cleanupDB(self,seqfile,deletesource=True)
                ## ~ [2a] Cluster by BLAST hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cluster = {}    # Dictionary of {seq:hit seqs} for clustering
                for search in blast.search:
                    seq = seqdict[search.info['Name']]
                    cluster[seq] = []
                    for hit in search.hit: cluster[seq].append(seqdict[hit.info['Name']])
                #self.deBug(cluster)
                ## ~ [2b] Combine clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                clusters = []   # List of [seqs] in clusters
                for seq in seqlist.seqs():
                    if seq not in cluster: continue
                    newcluster = [seq]
                    hits = cluster.pop(seq)
                    while hits:
                        hit = hits.pop(0)
                        if hit not in newcluster: newcluster.append(hit)
                        if hit in cluster: hits += cluster.pop(hit)
                    clusters.append(newcluster)
            self.printLog('#CLUSTER','%d clusters of %s proteins hits' % (len(clusters),searchset))
            #self.deBug(clusters)

            ### ~ [3] Assign peptides to consensi as "Common", "Cluster" or "Unique" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Match peptides to sequence lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pepcons = {}
            for seq in seqlist.seqs():
                prot = seq.shortName()  #.info['AccNum']
                for pep in data[prot]['conpep']:
                    if pep not in pepcons: pepcons[pep] = []
                    pepcons[pep].append(seq)
            self.dict['PepSeq'] = pepcons
            ## ~ [3b] Classify peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['PepTypes'] = {'Common':[],'Cluster':[],'Unique':[]}
            for pep in pepcons:
                if len(pepcons[pep]) == 1: self.dict['PepTypes']['Unique'].append(pep); continue
                pepclus = []
                for seq in pepcons[pep]:
                    for cluster in clusters:
                        if seq in cluster and cluster not in pepclus: pepclus.append(cluster)
                if len(pepclus) == 1: self.dict['PepTypes']['Cluster'].append(pep)
                else: self.dict['PepTypes']['Common'].append(pep)
            ## ~ [3c] Summarise Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#PEP','%d different %s Peptide sequences' % (len(pepcons),searchset))
            for ptype in ['Common','Cluster','Unique']: self.dict['PepTypes'][ptype].sort()
            self.printLog('#UNIQ','%d Unique to one consensus' % (len(self.dict['PepTypes']['Unique'])))
            self.printLog('#CLUS','%d Resticted to one cluster' % (len(self.dict['PepTypes']['Cluster'])))
            self.printLog('#COMM','%d Common to multiple clusters' % (len(self.dict['PepTypes']['Common'])))

            ### ~ [4] Update dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cx = 0
            for cluster in clusters:
                cx += 1
                for seq in cluster:
                    prot = seq.shortName()  #info['AccNum']
                    data[prot]['cluster'] = cx

            ### ~ [5] Peptide Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            peptdt = '%s.%s.peptides.tdt' % (self.info['Basefile'],searchset)
            pephead = ['Peptide','Classification','Hits']
            rje.delimitedFileOutput(self,peptdt,pephead,rje_backup=True)
            for ptype in ['Common','Cluster','Unique']:
                for pep in self.dict['PepTypes'][ptype]:
                    data = {'Peptide':pep,'Classification':ptype,'Hits':seqlist.accList(self.dict['PepSeq'][pep])}
                    data['Hits'].sort()
                    data['Hits'] = string.join(data['Hits'],'|')
                    rje.delimitedFileOutput(self,peptdt,pephead,datadict=data)
            self.printLog('#PEP','Peptide details output to %s' % peptdt)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: PICSI Class                                                                                      #
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
    try: PICSI(mainlog,cmd_list).run()

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
