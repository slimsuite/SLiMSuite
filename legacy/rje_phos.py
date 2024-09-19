#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_phos
Description:  RJE Phosphorylation Module
Version:      0.1
Last Edit:    26/10/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for the control and management of different phosphorylation sources and for mapping known phosphosites
    onto other proteins. Initially, the reading and processing of PhosphoELM is the only source implemented, although
    phosphorylation prediction may be added with time.

    The phosBLAST method will take an input dataset and BLAST them against the phosphorylation database(s), align hits
    above a certain GABLAM %ID and mark the phosphosites on the sequence. A second %ID cut-off will be used to
    determine which sites are marked as identities and which as homologies. Outputs will include an alignment with all
    the phosphoHomologues plus a simpler "alignment" of the query protein and marked sites. Phosphorylation predictions
    may be added to this second simpler "alignment". A UniProt-format file could also be produced with phosphoSites
    marked.

Commandline:
    ### ~ PhosphoELM Input Options (Creates Fasta file if pELM given) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    pelm=FILE       : Filename for phosphoELM download [None]
    filterseq=T/F   : Apply rje_seq sequence filters to phosphoELM data [False]
    pelmfas=FILE    : Filename for fasta file output of pELM sequences [pelm.fas]
    ### ~ PhosBLAST Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    phosblast=FILE  : Fasta file of sequences to perform phosBLAST method against pELM [None]
    usespec=T/F     : Use species codes for determing same species for ID matches [True]
    idsim=X         : Percentage identity (GABLAM; phosblast qry) for marking as identity [95.0]
    homsim=X        : Percentage identity (GABLAM; phosblast qry) for marking as homologue [40.0]
    phosdat=T/F     : Whether to produce a modified UniProt-format file with potential phosphoSites as features [False]
    phosres=FILE    : Delimited text file containing input sequence, position and evidence [*.phosres.tdt]

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_seq
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_uniprot, rje_zen
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation. Basic pELM parsing done.
    # 0.1 - Added phosBLAST method.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Add post-loading filtering of sequences. [Needs testing]
    # [Y] : Add UniProt cross-referencing for sequence details.
    # [ ] : Add method to BLAST sequences against phosphoDB, align hits and mark phosphosites (ID & Homology)
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_PHOS', '0.0', 'October 2007', '2007')
    description = 'RJE Phosphorylation Module'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
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
### SECTION II: PhosphoSeq Class                                                                                        #
#########################################################################################################################
class PhosphoSeq(rje.RJE_Object):     
    '''
    PhosphoSeq Class. Author: Rich Edwards (2007).

    Info:str
    - PELM = Filename for phosphoELM download [None]
    - PELMFas = Filename for fasta file output of pELM sequences [pelm.fas]
    - PhosBlast = Fasta file of sequences to perform phosBLAST method against pELM [None]
    - PhosRes = Delimited text file containing input sequence, position and evidence [*.phosres.tdt]

    Opt:boolean
    - FilterSeq = Apply rje_seq sequence filters to phosphoELM data [False]
    - UseSpec = Use species codes for determing same species for ID matches [True]
    - PhosDat = Whether to produce a modified UniProt-format file with potential phosphoSites as features [False]

    Stat:numeric
    - IDSim = Percentage identity (GABLAM; phosblast qry) for marking as identity [95.0]
    - HomSim = Percentage identity (GABLAM; phosblast qry) for marking as homologue [40.0]

    List:list

    Dict:dictionary
    - PhosphoSites = Dictionary of {Seq:{Pos:details}}

    Obj:RJE_Objects
    - SeqList = rje_seq.SeqList() object for storing sequences
    - UniProt = rje_uniprot.UniProt() object for storing UniProt data
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['PELM','PELMFas','PhosBlast','PhosRes']
        self.optlist = ['FilterSeq','UseSpec','PhosDat']
        self.statlist = ['IDSim','HomSim']
        self.listlist = []
        self.dictlist = ['PhosphoSites']
        self.objlist = ['SeqList','UniProt']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'PELMFas':'pelm.fas'})
        self.setStat({'IDSim':95.0,'HomSim':40.0})
        self.setOpt({'UseSpec':True})
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
                self._cmdReadList(cmd,'file',['PELM','PELMFas','PhosBlast','PhosRes'])
                self._cmdReadList(cmd,'opt',['FilterSeq','UseSpec','PhosDat'])
                self._cmdReadList(cmd,'stat',['IDSim','HomSim'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <3> ### Main Run Methods                                                                                        #
#########################################################################################################################
    def run(self):  ### Main method for standalone functionality
        '''Main method for standalone functionality.'''
        self.readPELM()
        if self.info['PhosBlast'].lower() not in ['','none']: self.mapPhosByBLAST(self.info['PhosBlast'])
#########################################################################################################################
    def readPELM(self): ### Reads phosphoELM into classes. Extracts UniProt data if available for Species etc.
        '''Reads phosphoELM into classes. Extracts UniProt data if available for Species etc.'''
        try:### ~ [1] Setup & Read File into Data Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = rje.dataDict(self,self.info['PELM'],mainkeys=['acc','position'])
            seqdict = {}    # Dictionary of Acc:Sequence

            ### ~ [2] Generate PhosphoSites dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdict = self.dict['PhosphoSites']
            for dkey in data:
                ## ~ [2a] Basic acc, seq and pos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (acc,pos) = rje.split(dkey)
                pos = string.atoi(pos)
                if acc not in pdict: pdict[acc] = {}
                if pos not in pdict[acc]: pdict[acc][pos] = {}
                ## ~ [2b] PhosphoELM data with checks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if acc not in seqdict: seqdict[acc] = data[dkey]['sequence']
                elif seqdict[acc] != data[dkey]['sequence']: self.log.printLog('#ERR','Warning. Sequence mismatch for %s' % acc)
                if 'aa' not in pdict[acc][pos]: pdict[acc][pos]['aa'] = data[dkey]['code']
                elif pdict[acc][pos]['aa'] != data[dkey]['code']: self.log.printLog('#ERR','Warning. PhosphoSite mismatch for %s at pos %d: %s not %s' % (acc,pos,data[dkey]['code'],pdict[acc][pos]['aa']))
                if data[dkey]['code'] != seqdict[acc][(pos-1):pos]: self.log.printLog('#ERR','Warning. PhosphoSeq mismatch for %s at pos %d: %s not %s' % (acc,pos,data[dkey]['code'],seqdict[acc][pos-1:pos]))

            ### ~ [3] Make sequence objects and update PhosphoSites keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Setup objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            acclist = rje.sortKeys(seqdict)
            pelmuni = rje_uniprot.UniProt(self.log,self.cmd_list)   # UniProt entry
            unidict = pelmuni.accDict(acclist)        # Dictionary of {acc:UniProtEntry}
            pelmseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])            # SeqList object
            ## ~ [3b] Add one sequence for each AccNum and update seqdict  ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Look out for splice variants! (There are some!) - Copy UniProt and change sequence & AccNum #!#
            for acc in acclist:     #!# Make accdict of {acc:Seq} using unidict and seqlist #!#
                sequence = seqdict[acc]
                try:
                    uni = unidict[rje.split(acc,'-')[0]]
                    desc = uni.obj['Sequence'].info['Description']
                    name = '%s__%s %s' % (uni.obj['Sequence'].info['ID'],acc,desc)
                    if sequence != uni.obj['Sequence'].info['Sequence']:
                        self.log.printLog('#WARNING','Sequence mismatch for UniProt entry %s' % acc)
                except:
                    self.log.errorLog('Problem with %s' % acc)
                    name = '%s_UNK__%s' % (acc,acc)             #!# Add sequences where UniProt missing #!#
                seqdict[acc] = pelmseq._addSeq(name,sequence)
            ## ~ [3c] Filtering of sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['FilterSeq']:
                pelmseq.autoFilter()
                for acc in acclist:
                    if seqdict[acc] not in pelmseq.seq: seqdict.pop(acc)
                acclist = rje.sortKeys(seqdict)
            ## ~ [3d] Save sequences for BLASTing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not os.path.exists(self.info['PELMFas']) or self.stat['Interactive'] < 0 or rje.yesNo('%s exists: overwrite?' % self.info['PELMFas']):
                pelmseq.saveFasta(seqfile=self.info['PELMFas'])
            self.obj['SeqList'] = pelmseq
            self.obj['UniProt'] = pelmuni
        except: self.log.errorLog('Problem during PhosphoSeq.readPELM')
#########################################################################################################################
    def mapPhosByBLAST(self,fasfile):   ### BLAST sequences against phosphoDB, align hits & mark sites (ID & Homology)
        '''BLAST sequences against phosphoDB, align hits and mark phosphosites (ID & Homology).'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup fasfile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scmd = self.cmd_list + ['seqin=%s' % fasfile,'autoload=T','autofilter=F']
            qseqlist = rje_seq.SeqList(self.log,scmd)
            qdict = qseqlist.seqNameDic()
            ## ~ [1b] Setup results files/directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            basefile = rje.baseFile(fasfile)
            if self.info['PhosRes'].lower() in ['','none']: self.info['PhosRes'] = '%s.phosres.tdt' % basefile
            headers = ['Name','Pos','AA','PELM','PELMPos','Evidence']
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.info['PhosRes']))
            rje.delimitedFileOutput(self,self.info['PhosRes'],headers,delimit,rje_backup=True)
            ppath = rje.makePath('PhosALN')
            rje.mkDir(self,ppath)
            ## ~ [1c] Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pblast = rje_blast.BLASTRun(self.log,self.cmd_list+['formatdb=F'])
            pblast.setInfo({'Name':'%s.p.blast' % rje.baseFile(fasfile),'DBase':self.info['PELMFas'],'InFile':fasfile})
            pblast.setStat({'HitAln':pblast.stat['OneLine']})
            pblast.opt['Complexity Filter'] = False
            pblast.formatDB(force=False)
            ## ~ [1d] Setup GABLAM Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gkey = 'GABLAMO ID' #x# % self.info['GABLAMO Key']
            for g in ['ID','Hom']:
                if self.stat['%sSim' % g] < 1.0: self.stat['%sSim' % g] *= 100.0
                self.stat['%sSim' % g] = max(0.0,self.stat['%sSim' % g])

            ### ~ [2] PhosphoBLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pblast.blast(use_existing=True,log=True)    # BLAST
            pblast.readBLAST(gablam=True)               # Read in
            while pblast.search:
                ## ~ [2a] Align relevant hits from each BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                search = pblast.search.pop(0)
                qseq = qdict[search.info['Name']]
                idlist = []
                qlen = qseq.aaLen()
                hitdict = search.hitSeq(self.obj['SeqList'])
                aln = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F','autofilter=F'])
                aln.seq = [qseq]
                pdict = {}      # Dictionary of {hseq:[poslist]}
                rdict = {qseq:0}      # Dictionary of {hseq:res}
                for hit in search.hit[0:]:
                    hseq = hitdict[hit]
                    pdict[hseq] = []
                    for pos in rje.sortKeys(self.dict['PhosphoSites'][hseq.info['AccNum']]): pdict[hseq].append(pos)
                    if hit.info['Name'] == search.info['Name']:
                        if qseq.getSequence(case=False,gaps=False) != hseq.getSequence(case=False,gaps=False):
                            self.log.errorLog('Major problem: Search/Hit sequence mismatch for same sequence "%s"' % hit.info['Name'])
                        idlist.append(qseq)
                        pdict[qseq] = pdict.pop(hseq)
                        continue
                    gdict = hit.globalFromLocal(qlen)
                    qvh = float(100 * gdict['Query'][gkey]) / float(qlen)
                    if qvh < self.stat['HomSim']:
                        pdict.pop(hseq)
                        continue
                    aln.seq.append(hseq)
                    if (qseq.sameSpec(hseq) or not self.opt['UseSpec']) and qvh >= self.stat['IDSim']: idlist.append(hseq)
                    rdict[hseq] = 0
                aln.muscleAln()   #x#outfile='%s%s.phosaln.fas' % (ppath,qseq.info['AccNum']))
                aln._addSeq('PhosAln','-' * qseq.seqLen())
                aln.info['Name'] = '%s%s.phosaln.fas' % (ppath,qseq.info['AccNum'])
                ## ~ [2b] Map phosphorylations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                print '>>>\n', aln.seq, pdict.keys(), rdict.keys()
                for a in range(qseq.seqLen()):
                    if qseq.info['Sequence'][a] != '-': rdict[qseq] += 1
                    for hseq in pdict:
                        if hseq.info['Sequence'][a] == '-': continue
                        if hseq != qseq: rdict[hseq] += 1
                        if rdict[hseq] in pdict[hseq] and qseq.info['Sequence'][a] == hseq.info['Sequence'][a]:  # Phosphosite
                            pdata = {'Name':search.info['Name'],'Pos':rdict[qseq],'AA':qseq.info['Sequence'][a],
                                     'PELM':hseq.shortName(),'PELMPos':rdict[hseq],'Evidence':'Hom'}
                            if hseq == qseq: pdata['Evidence'] = 'Self'
                            elif hseq in idlist: pdata['Evidence'] = 'ID'
                            rje.delimitedFileOutput(self,self.info['PhosRes'],headers,delimit,pdata)
                            self.addPhos(aln.seq[-1],a,pdata['Evidence'])
                ## ~ [2c] Add Scansite/NetPhos if made? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## ~ [2d] Save alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                aln.saveFasta()


            # Align hits for each > X %ID
            # Map phosphosites onto alignment and output #
            
            return
        except: self.log.errorLog('Problem during PhosphoSeq.mapPhosByBLAST')
#########################################################################################################################
    def addPhos(self,seq,pos,type): ### Adds relevant marker to given sequence
        '''Adds relevant marker to given sequence.'''
        hierarchy = ['*','?','H','M','L','-']
        markers = {'Hom':'?','ID':'*','Self':'*'}
        p = markers[type]
        x = seq.info['Sequence'][pos]
        if x not in hierarchy or hierarchy.index(x) > hierarchy.index(p):
            seq.info['Sequence'] = rje.strSub(seq.info['Sequence'],pos,pos,p)
#########################################################################################################################
### End of SECTION II: PhosphoSeq Class                                                                                 #
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
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
        
    ### Rest of Functionality... ###
    try: PhosphoSeq(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
