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
Module:       rje_gquad
Description:  Custom G-Quadruplex Analysis Tool
Version:      0.0
Last Edit:    19/03/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This is a customised analysis pipeline. Much of this might be moved into a FlyBase Class and module.

Commandline:
    ### ~ MAIN FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    makeflyseq=T/F  : Re-annotate D.mel gene sequences with CDS and Exon positions [False]
    codons=T/F      : Perform codon/triplet analysis [False]

    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ... None yet!
    
Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_sequence, rje_zen
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
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('G-QUAD', '0.0', 'March 2008', '2008')
    description = 'Custom G-Quadruplex Analysis Tool'
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
### SECTION II: GQuad Class                                                                                             #
#########################################################################################################################
class GQuad(rje.RJE_Object):     
    '''
    GQuad Class. Author: Rich Edwards (2008).

    Info:str
        
    Opt:boolean
    - Codons = Perform codon/triplet analysis [False]
    - MakeFlySeq = Re-annotate D.mel gene sequences with CDS and Exon positions [False]

    Stat:numeric

    List:list

    Dict:dictionary
    - NTFreq = Nucleotide frequencies for CDS {nt:freq}
    - CodonFreq = Dictionary of {aa:{codon:freq}} based on genetic code and NTFreq

    Obj:RJE_Objects
    - Exons = SeqList object containing Exon sequences
    - Genes = SeqList object containing Gene sequences
    - Transcripts = SeqList object containing transcript (CDS?) sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = []
        self.optlist = ['Codons','MakeFlySeq']
        self.statlist = []
        self.listlist = []
        self.dictlist = ['NTFreq','CodonFreq']
        self.objlist = ['SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
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
                self._cmdReadList(cmd,'opt',['Codons','MakeFlySeq'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Run Method                                                                                         #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['MakeFlySeq']: self.makeFlySeq()
            if self.opt['Codons']: self.codons()
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <3> ### Sequence Generation Methods                                                                             #
#########################################################################################################################
    def makeFlySeq(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flybase = rje.makePath('/scratch/Databases/NewDB/FlyBase/Fasta/')
            scmd = ['accnr=F','seqnr=F','gnspacc=F']
            genes = rje_seq.SeqList(self.log, self.cmd_list+['seqin=%sdmel-all-gene-r5.5.fasta' % flybase]+scmd)
            cds = rje_seq.SeqList(self.log, self.cmd_list+['seqin=%sdmel-all-CDS-r5.5.fasta' % flybase]+scmd)
            exons = rje_seq.SeqList(self.log, self.cmd_list+['seqin=%sdmel-all-exon-r5.5.fasta' % flybase]+scmd)

            ### ~ [1] ~	Read in full-length gene and note start and end positions in parent scaffold ~~~~~~~~~~~~~~~~ ###
            genedict = {}   # Dictionary of {ID:Sequence object}
            (gx,gtot) = (0.0,genes.seqNum())
            for gene in genes.seq:
                self.log.printLog('\r#GENE','Processing Gene Annotation: %.1f%%' % (gx/gtot),newline=False,log=False)
                gx += 100
                (id,scaffold,pos,name,glen) = rje.matchExp('^(\S+)\s.+loc=(\S+):(\S+);.+name=(\S+);.+length=(\d+);',gene.info['Name'])
                if string.atoi(glen) != gene.aaLen(): self.log.errorLog('%s Length mismatch!' % id, printerror=False)
                genedict[id] = gene
                gene.setInfo({'Scaffold':scaffold,'Gene':name})
                try: (end,start) = rje.matchExp('^complement\((\d+)\.\.(\d+)\)',pos)
                except: (start,end) = rje.matchExp('^(\d+)\.\.(\d+)',pos)
                (start,end) = (string.atoi(start),string.atoi(end))
                gene.opt['Complement'] = start > end        # Sequence on "lagging" strand
                gene.setStat({'Start':start,'End':end})
                gene.list['CDS'] = []       # Will add CDS sequences here
                gene.list['Exon'] = []      # Will add exon sequences here
            self.log.printLog('\r#GENE','Processing Gene Annotation complete!')
                           
            ### ~ [2] ~ Read in associated CDS sequences and note start and end positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (cx,ctot) = (0.0,cds.seqNum())
            for seq in cds.seq:
                self.log.printLog('\r#CDS','Processing CDS Annotation: %.1f%%' % (cx/ctot),newline=False,log=False)
                cx += 100
                try: (id,scaffold,pos,name,glen,parent) = rje.matchExp('^(\S+)\s.+loc=(\S+):(\S+);.+name=(\S+);.+length=(\d+);.+parent=(\S+),\S+;',seq.info['Name'])
                except:
                    self.log.errorLog(seq.info['Name'])
                    raise
                if string.atoi(glen) != seq.aaLen(): self.log.errorLog('%s Length mismatch!' % id, printerror=False)
                seq.obj['Parent'] = gene = genedict[parent]
                try: (end,start) = rje.matchExp('^complement\((\d+)\..*\.(\d+)\)',pos)
                except:
                    try: (start,end) = rje.matchExp('^join\((\d+)\..*\.(\d+)\)',pos)
                    except: (start,end) = rje.matchExp('^(\d+)\.\.(\d+)',pos)
                (start,end) = (string.atoi(start),string.atoi(end))
                seq.opt['Complement'] = start > end        # Sequence on "lagging" strand
                seq.setStat({'Start':start,'End':end})
                gene.list['CDS'].append(seq)
            self.log.printLog('\r#CDS','Processing CDS Annotation complete!')
                
            ### ~ [3] ~ Read in associated exons and note start and end positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ex,etot) = (0.0,exons.seqNum())
            for seq in exons.seq:
                self.log.printLog('\r#EXON','Processing Exon Annotation: %.1f%%' % (ex/etot),newline=False,log=False)
                ex += 100
                try: (id,scaffold,pos,name,parent) = rje.matchExp('^(\S+)\s.+loc=(\S+):(\S+);.+name=(\S+);.+parent=(\S+);',seq.info['Name'])
                except:
                    self.log.errorLog(seq.info['Name'])
                    raise
                seq.obj['Parent'] = gene = genedict[string.split(parent,',')[0]]
                try: (end,start) = rje.matchExp('^complement\((\d+)\..*\.(\d+)\)',pos)
                except:
                    try: (start,end) = rje.matchExp('^join\((\d+)\..*\.(\d+)\)',pos)
                    except: (start,end) = rje.matchExp('^(\d+)\.\.(\d+)',pos)
                (start,end) = (string.atoi(start),string.atoi(end))
                seq.opt['Complement'] = start > end        # Sequence on "lagging" strand
                seq.setStat({'Start':start,'End':end})
                gene.list['Exon'].append(seq)
            self.log.printLog('\r#EXON','Processing Exon Annotation complete!')
                
            ### ~ [4] ~ Regenerate output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] ~ Convert to relative positions and store ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (gx,gtot) = (0.0,genes.seqNum())
            for gene in genes.seq:
                glen = gene.aaLen()
                self.log.printLog('\r#GENE','Generating new Gene Annotation: %.1f%%' % (gx/gtot),newline=False,log=False)
                gx += 100
                clist = []
                for seq in gene.list['CDS']:
                    if gene.opt['Complement']:  # Must substract from "wrong" end and reverse
                        start = gene.stat['Start'] - seq.stat['Start']
                        end = gene.stat['Start'] - seq.stat['End']
                    else:
                        start = seq.stat['Start'] - gene.stat['Start']
                        end = seq.stat['End'] - gene.stat['Start']
                    pos = '%s-%s' % (rje.preZero(start,glen),rje.preZero(end,glen))
                    clist.append(pos)
                clist = rje.sortUnique(clist,xreplace=False)
                elist = []
                for seq in gene.list['Exon']:
                    if gene.opt['Complement']:  # Must substract from "wrong" end and reverse
                        start = gene.stat['Start'] - seq.stat['Start']
                        end = gene.stat['Start'] - seq.stat['End']
                    else:
                        start = seq.stat['Start'] - gene.stat['Start']
                        end = seq.stat['End'] - gene.stat['Start']
                    pos = '%s-%s' % (rje.preZero(start,glen),rje.preZero(end,glen))
                    elist.append(pos)
                elist = rje.sortUnique(elist,xreplace=False)
                gene.info['Name'] = '%s_%s__%s Length=%d; CDS=%s; Exons=%s;' % (gene.info['Gene'],gene.info['SpecCode'],gene.info['AccNum'],gene.aaLen(),string.join(clist,','),string.join(elist,','))
            self.log.printLog('\r#GENE','Generating new Gene Annotation complete!')
            ## ~ [4b] ~ Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            genes.saveFasta(seqfile='flybase_DROME.genes.fas')

        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### Codons Methods                                                                                          #
#########################################################################################################################
    def codons(self):  ### Main codons analysis method
        '''Main codons analysis method.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flybase = rje.makePath('/scratch/Databases/NewDB/FlyBase/Fasta/')
            scmd = ['accnr=F','seqnr=F','gnspacc=F']
            cds = rje_seq.SeqList(self.log, self.cmd_list+['seqin=%sdmel-all-CDS-r5.5.fasta' % flybase]+scmd)
            gcode = rje_sequence.genetic_code

            ### ~ [1] ~ Make codon frequency tables (a) Observed, (b) Based on NTFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nts = ['A','C','G','T']
            ntfreq = cds.aaFreq(alphabet=nts)
            codons = []     # List of codons
            obs_cfreq = {}  # Observed codon frequencies
            nts_cfreq = {}  # Codon frequencies from NT frequencies
            obs_tfreq = {}  # Observed triplet frequencies
            nts_tfreq = {}  # Predicted triplet frequencies from NT frequencies
            ocd_tfreq = {}  # Predicted triplet frequencies from observed codon frequencies
            ncd_tfreq = {}  # Predicted triplet frequencies from nt-predicted codon frequencies
            ## ~ [1a] ~ Setup dictionaries using nt freqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for n1 in nts:
                for n2 in nts:
                    for n3 in nts:
                        cod = '%s%s%s' % (n1,n2,n3)
                        codons.append(cod)
                        aa = gcode[string.replace(cod,'T','U')]
                        if aa not in obs_cfreq: obs_cfreq[aa] = {}
                        if aa not in nts_cfreq: nts_cfreq[aa] = {}
                        obs_cfreq[aa][cod] = 0.0
                        nts_cfreq[aa][cod] = ntfreq[n1] * ntfreq[n2] * ntfreq[n3]
                        obs_tfreq[cod] = 0.0
                        nts_tfreq[cod] = ntfreq[n1] * ntfreq[n2] * ntfreq[n3]
                        ocd_tfreq[cod] = 0.0
                        ncd_tfreq[cod] = 0.0
            nts_tfreq = rje.dictFreq(nts_tfreq,total=False)                                 # Normalise triplet freq.
            for aa in nts_cfreq: nts_cfreq[aa] = rje.dictFreq(nts_cfreq[aa],total=False)    # Normalise codon freq.
            self.log.printLog('#FREQ','Frequency dictionaries set up.')
            ## ~ [1b] ~ Observed codon freq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (sx,stot) = (0.0,cds.seqNum())
            for seq in cds.seq[0:]:
                self.log.printLog('\r#OBS','Calculating observed codon frequencies: %.1f%%' % (sx/stot),newline=False,log=False)
                sx += 100.0
                try: (id,scaffold,pos,name,glen,parent) = rje.matchExp('^(\S+)\s.+loc=(\S+):(\S+);.+name=(\S+);.+length=(\d+);.+parent=(\S+),\S+;',seq.info['Name'])
                except:
                    self.log.errorLog(seq.info['Name'])
                    raise
                try: exons = rje.matchExp('^complement\((\d+\..*\.\d+)\)',pos)[0]
                except:
                    try: exons = rje.matchExp('^join\((\d+\..*\.\d+)\)',pos)[0]
                    except: exons = rje.matchExp('^(\d+\.\.\d+)',pos)[0]
                self.deBug(exons)
                exons = string.split(exons,',')
                elen = []
                try:
                    for exon in exons:
                        (start,end) = string.split(exon,'..')
                        elen.append(string.atoi(end) - string.atoi(start) + 1)
                except:
                    self.log.errorLog(id)
                    cds.seq.remove(seq)
                    continue
                        
                if pos[:4] == 'comp': elen.reverse()
                seq.list['ExonLen'] = elen
                self.deBug(elen)
                if sum(elen) != seq.aaLen(): self.log.errorLog('%s exon length error' % id,printerror=False)
                if seq.aaLen()/3 != seq.aaLen()/3.0:
                    self.log.errorLog('%s not a multiple of 3nt long!' % id,printerror=False)
                    cds.seq.remove(seq)
                    continue
                #!# Add use exon option - single full-length exon if false (mature mRNA) #!#
                sequence = seq.info['Sequence'][0:]
                if string.count(sequence,'N') > 0:
                    self.log.errorLog('%s has 1+ Ns!' % id,printerror=False)
                    cds.seq.remove(seq)
                    continue
                while sequence:
                    cod = sequence[:3]
                    sequence = sequence[3:]
                    aa = gcode[string.replace(cod,'T','U')]
                    obs_cfreq[aa][cod] += 1
            for aa in obs_cfreq: obs_cfreq[aa] = rje.dictFreq(obs_cfreq[aa],total=False)    # Normalise codon freq.
            self.log.printLog('\r#OBS','Calculating observed codon frequencies complete.')

            ### ~ [2] ~ Generate Triplet freq. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,stot) = (0.0,cds.seqNum())
            for seq in cds.seq:
                self.log.printLog('\r#TRIP','Calculating triplet frequencies: %.1f%%' % (sx/stot),newline=False,log=False)
                sx += 100.0
                elen = seq.list['ExonLen'] 
                sequence = seq.info['Sequence'][0:]
                aa = ''
                cod = ''
                ax = 0      # Measure sequence length processed for exon boundary checks
                while sequence:
                    prevcod = cod
                    cod = sequence[:3]
                    prevaa = aa
                    sequence = sequence[3:]
                    aa = gcode[string.replace(cod,'T','U')]
                    ## ~ [2a] ~ Predicted Triplet Freq. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for cod2 in obs_cfreq[aa]:
                        if elen[0] > ax + 3:    # Exon boundary beyond this codon
                            ocd_tfreq[cod2] += obs_cfreq[aa][cod2]
                            ncd_tfreq[cod2] += nts_cfreq[aa][cod2]
                        if prevaa:              # Look at overlap with previous codon
                            for cod1 in obs_cfreq[prevaa]:
                                for i in range(1,3):
                                    if elen[0] > ax + i:    # Exon boundary beyond overlap
                                        acod = cod1[i:] + cod2[:i]
                                        ocd_tfreq[acod] += (obs_cfreq[prevaa][cod1] * obs_cfreq[aa][cod2])
                                        ncd_tfreq[acod] += (nts_cfreq[prevaa][cod1] * nts_cfreq[aa][cod2])
                    ## ~ [2b] ~ Observed Triplet Freq. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if elen[0] > ax + 3:    # Exon boundary beyond this codon
                        obs_tfreq[cod] += 1
                    if prevcod:              # Look at overlap with previous codon
                        for i in range(1,3):
                            if elen[0] > ax + i:    # Exon boundary beyond overlap
                                acod = prevcod[i:] + cod[:i]
                                obs_tfreq[acod] += 1
                    # Check exons #
                    ax += 3
                    if ax >= elen[0]: ax -= elen.pop(0)
            obs_tfreq = rje.dictFreq(obs_tfreq,total=False)
            ocd_tfreq = rje.dictFreq(ocd_tfreq,total=False)
            ncd_tfreq = rje.dictFreq(ncd_tfreq,total=False)    
            self.log.printLog('\r#TRIP','Calculating triplet frequencies complete.')

            ### ~ [3] ~ Output results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['Triplet','AA','Degen','Obs_Codon','NT_Codon','Obs_Trip','NT_Trip','ObCod_Trip','NTCod_Trip']
            tfile = 'quad_triplet.tdt'
            rje.delimitedFileOutput(self,tfile,headers,rje_backup=True)
            for cod in codons:
                aa = gcode[string.replace(cod,'T','U')]
                datadict = {'Triplet':cod,'AA':aa,'Degen':len(obs_cfreq[aa]),'Obs_Codon':obs_cfreq[aa][cod],
                            'NT_Codon':nts_cfreq[aa][cod],'Obs_Trip':obs_tfreq[cod],'NT_Trip':nts_tfreq[cod],
                            'ObCod_Trip':ocd_tfreq[cod],'NTCod_Trip':ncd_tfreq[cod]}
                rje.delimitedFileOutput(self,tfile,headers,datadict=datadict)
            self.log.printLog('#OUT','Triplet & codon data output to %s' % tfile)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: GQuad Class                                                                                      #
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
    try: GQuad(mainlog,cmd_list).run()
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
