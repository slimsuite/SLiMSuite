#!/usr/local/bin/python

# rje_tm.py - Transmembrane and Signal Peptide prediction module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#  
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#  
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_tm    
Description:  Tranmembrane and Signal Peptide Prediction Module
Version:      1.2
Last Edit:    16/08/07
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Will read in results from tmhmm and/or signalp files as appropriate and append output to:
    - tm.tdt        = TM domain counts and orientation
    - domains.tdt   = Domain table
    - singalp.tdt   = SingalP data (use to add signal peptide domains to domains table using mySQL

Commandline:
    tmhmm=FILE  : TMHMM output file [None]
    signalp=FILE: SignalP output file [None]
    mysql=T/F   : Output results in tdt files for mySQL import [True]

    seqin=FILE      : Sequence file for which predictions have been made [None]
    maskcleave=T/F  : Whether to output sequences with cleaved signal peptides masked. [False]
    source=X        : Source text for mySQL file ['tmhmm']

Uses general modules: os, string, sys, threading, time
Uses RJE modules: rje, rje_seq
Other modules required: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import os, re, string, sys, time
#############################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_seq, rje_zen
#############################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Basic functional version.
    # 1.1 - Output of 'truncated' proteins - cleaved signalp peptides replaced with Xs.
    # 1.2 - Added standalone parsing methods as first step in development (& for rje_ensembl.ensDat())
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Commandline option for running TMHMM and SignalP
    # [ ] : Forking for TMHMM and SignalP -> split input, fork, compile, delete
    # [ ] : Remodel entire TM class to be in line with new module structure
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_TM', '1.2', 'August 2007', '2007')
    description = 'Transmembrane and SignalP Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is under development and may contain bugs!',rje_zen.Zen().wisdom()]
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
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)
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
### SECTION II: TM Class                                                                                                # 
#########################################################################################################################
class TM(rje.RJE_Object):     
    '''
    TM Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Sequence file for which predictions have been made [None]
    - TMHMM = File for TMHMM output
    - SignalP = File for SignalP output
    - Source = Source text for mySQL file ['tmhmm']
    
    Opt:boolean
    - MySQL = Whether to output files to mySQL tables (*.tdt)
    - MaskCleave = Whether to output sequences with cleaved signal peptides masked. [False]

    Stat:numeric

    Obj:RJE_Objects
    - SeqList = rje_seq.SeqList object.
    '''
    ### Attributes
    tmhmm = {}      # Dictionary of AccNum:results (itself a dictionary)
    signalp = {}    # Dictionary of AccNum:results (itself a dictionary)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','TMHMM','SignalP','Source']
        - Stats:float []
        - Opt:boolean ['MySQL','MaskCleave']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','TMHMM','SignalP','Source']
        self.statlist = []
        self.optlist = ['MySQL','MaskCleave']
        self.objlist = ['SeqList']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        self.obj['SeqList'] = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=F','autofilter=F'])
        self.info['Source'] = 'tmhmm'
        ### <c> ### Other Attributes
        self.tmhmm = {}
        self.signalp = {}
#########################################################################################################################
#    def _cmdRead(self,type='info',att=None,arg=None,cmd=None):     ### Sets self.type[att] from commandline command cmd
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='Name',arg='seqin')
                self._cmdRead(cmd,type='info',att='TMHMM')
                self._cmdRead(cmd,type='info',att='SignalP')
                self._cmdRead(cmd,type='info',att='Source')
                self._cmdRead(cmd,type='opt',att='MySQL')
                self._cmdRead(cmd,type='opt',att='MaskCleave')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Text Parsing Methods
#########################################################################################################################
    def parseTMHMM(self,file=None):     ### Parses TMHMM into dictionary
        '''
        Parses TMHMM into dictionary self.tmhmm.
        >> file:str = will read from file if given, else self.info['TMHMM']
        '''
        try:
            ### <a> ### Setup
            _stage = '<a> Setup'
            tmhmm_pattern = '^(\S+)\s+(len.+)$'
            if file == None:
                file = self.info['TMHMM']
            if file == 'None':
                self.verbose(0,1,'No TMHMM file given.',2)
                return
            self.verbose(0,3,'Parsing TMHMM file %s...' % file,0)
            TMRES = open(file, 'r')
            ### <b> ### Read in
            _stage = '<b> Read in'
            while 1:
                tmline = re.sub('\t',' ',TMRES.readline())
                if tmline:
                    tmres = rje.matchExp(tmhmm_pattern,tmline)
                    if tmres:
                        acc = tmres[0]
                        if rje.matchExp('^\S+__(\S+)',acc):
                            acc = rje.matchExp('^\S+__(\S+)',acc)[0]
                        self.tmhmm[acc] = {}
                        reslist = string.split(tmres[1])
                        for res in reslist:
                            split = string.split(res,'=')
                            self.tmhmm[acc][split[0]] = split[1]
                else:
                    break
            TMRES.close()
            self.verbose(0,1,'Done!',2)
        except:
            self.log.errorLog('Problem with parseTMHMM() %s.' % _stage)
#########################################################################################################################
    def getDomains(self,acc):     ### Returns domains for a given sequence = list of dictionaries ['Type','Start','End']
        '''
        Returns domains for a given sequence = list of dictionaries ['Type','Start','End']
        >> file:str = will read from file if given, else self.info['TMHMM']
        '''
        try:
            ### <a> ### Setup
            _stage = '<a> Setup'
            domainlist = []
            domains = self.tmhmm[acc]['Topology']
            tm = False
            dom = 'CYTOPLASMIC'
            if domains[0] == 'o':
                dom = 'EXTRACELLULAR'
            domains = re.sub('o', '-', domains)
            domains = re.sub('i', '-', domains)
            domains = string.split('1' + domains + self.tmhmm[acc]['len'],'-')
            started = False
            while len(domains) > 1:
                _stage = '<b-i> TM Dom Write'
                start = domains.pop(0)
                end = domains[0]
                if tm:
                    type = 'TRANSMEMBRANE'
                else:
                    type = dom
                    if started:
                        #print start,
                        start = '%d' % (string.atoi(start) + 1)
                        #print start
                    else:
                        started = True
                    if len(domains) > 1:
                        #print end,
                        end = '%d' % (string.atoi(end) - 1)
                        #print end
                domainlist.append({'Type':type,'Start':start,'End':end})
                if tm:
                    tm = False
                    if dom == 'CYTOPLASMIC':
                        dom = 'EXTRACELLULAR'
                    else:
                        dom = 'CYTOPLASMIC'
                else:
                    tm = True
            _stage = '<b-ii> TM Dom Check'
            if tm == False:
                self.log.errorLog('Problem with %s TM domains - wrong number of domains!' % acc)
                return []
            return domainlist
        except:
            self.log.errorLog('Problem with getDomains() %s.' % _stage)
#########################################################################################################################
    def parseSignalP(self,file=None):     ### Parses SignalP into dictionary
        '''
        Parses SignalP into dictionary self.signalp. This takes heavily from Paul's SignalP module. SignalP should be
        run using the "-f short" option to give simple one-line predictions per protein.
        >> file:str = will read from file if given, else self.info['SignalP']
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nnNFields = 14 # number of fields in a neural networks result
            hmmNFields = 7 # number of fields in a hidden markov model result
            if file == None: file = self.info['SignalP']
            if file.lower() in ['none','']: return self.log.errorLog('Cannot parse SignalP - no results file given',printerror=False)
            elif not os.path.exists(file): return self.log.errorLog('Cannot parse SignalP - file "%s" not found' % file,printerror=False)
            self.log.printLog('#SIGP','Parsing SignalP file %s...' % file,log=False,newline=False)

            ### ~ [2] Read in results file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for iline in self.loadFromFile(file,v=2,chomplines=True):
                if not iline or re.match(r"^#", iline): continue
                line = re.sub("\t", " ", iline)
                
                # Clean up and build list of fields
                fields = []
                for field in line.split(' '):
                    field.strip()
                    if field not in ['', None]: fields.append(field)
                ncols = len(fields)
                results = {}
                
                # Based on the number of fields in each line of the result, we can 
                # determine whether a hidden markov model was used, a neural network
                # was used, or both. We build our list of results accordingly.
                # XXX - This is a bit messy. It would be a good idea to find a cleaner
                #       way of doing this. 
                if ncols == hmmNFields:
                    # Hidden Markov Model result
                    results['exclamation'] = fields[1]
                    results['cmax']    = fields[2]
                    results['cmaxpos'] = fields[3]
                    results['cmax?']   = fields[4]
                    results['sprob']   = fields[5]
                    results['sprob?']  = fields[6]
                elif ncols == nnNFields:
                    # Neural Networks result
                    results['cmax']    = fields[1]
                    results['cmaxpos'] = fields[2]
                    results['cmax?']   = fields[3]
                    results['ymax']    = fields[4]
                    results['ymaxpos'] = fields[5]
                    results['ymax?']   = fields[6]
                    results['smax']    = fields[7]
                    results['smaxpos'] = fields[8]
                    results['smax?']   = fields[9]
                    results['smean']   = fields[10]
                    results['smean?']  = fields[11]
                    results['d']       = fields[12]
                    results['d?']      = fields[13]
                elif ncols == hmmNFields + nnNFields:
                    # Hidden Markove Model + Neural Networks result
                    results['nn_cmax']    = fields[1]
                    results['nn_cmaxpos'] = fields[2]
                    results['nn_cmax?']   = fields[3]
                    results['nn_ymax']    = fields[4]
                    results['nn_ymaxpos'] = fields[5]
                    results['nn_ymax?']   = fields[6]
                    results['nn_smax']    = fields[7]
                    results['nn_smaxpos'] = fields[8]
                    results['nn_smax?']   = fields[9]
                    results['nn_smean']   = fields[10]
                    results['nn_smean?']  = fields[11]
                    results['nn_d']       = fields[12]
                    results['nn_d?']      = fields[13]
                    fields[0] = fields[14]
                    results['hmm_exclamation'] = fields[15]
                    results['hmm_cmax']   = fields[16]
                    results['hmm_cmaxpos']= fields[17]
                    results['hmm_cmax?']  = fields[18]
                    results['hmm_sprob']  = fields[19]
                    results['hmm_sprob?'] = fields[20]
                else: continue
                
                self.signalp[fields[0]] = results

            self.log.printLog('\r#SIGP','Parsed %d SignalP results from %s.' % (len(self.signalp),file))
        except:
            self.log.errorLog('Problem with parseSignalP() %s.' % _stage)
#########################################################################################################################
    ### <3> ### MySQL Output
#########################################################################################################################
    def mySQLOut(self,tmfile='tm.tdt',domfile='domains.tdt',sigfile='signalp.tdt',makenew=False):     ### Output to tdt files
        '''
        Output to tdt files.
        >> tmfile:str = File to save TM numbers in (no save if None)
        >> domfile:str = File to save domain data in (no save if None)
        >> sigfile:str = File to save choice SignalP data in (no save if None)
        >> makenew:boolean [False] = whether to make new files (True) or append (False)
        '''
        try:
            ### <a> ### Setup
            _stage = '<a> Setup'
            siglist = ['nn_cleavemax','nn_cleavepos','nn_cleave','nn_sig_mean','nn_sig','nn_dscore','nn_d',
                       'hmm_cmax','hmm_cpos','hmm_cleave','hmm_sigprob','hmm_sig']

            self.deBug(self.tmhmm)
            if tmfile and self.tmhmm.keys():
                _stage = '<a-i> TM File'
                if makenew or os.access(tmfile, os.F_OK) == False:
                    TMFILE = open(tmfile, 'w')
                    TMFILE.write('acc_num\ttm\tnterm\tcterm\n')
                else:
                    TMFILE = open(tmfile, 'a')
            if domfile and (self.tmhmm.keys()+self.signalp.keys()):
                _stage = '<a-ii> Dom File'
                if makenew or os.access(domfile, os.F_OK) == False:
                    DOMFILE = open(domfile, 'w')
                    DOMFILE.write('acc_num\tdomain\tdom_start\tdom_end\tsource\n')
                else:
                    DOMFILE = open(domfile, 'a')
            if sigfile and self.signalp.keys():
                _stage = '<a-iii> Sig File'
                if makenew or os.access(sigfile, os.F_OK) == False:
                    SIGFILE = open(sigfile, 'w')
                    sheader = string.join(siglist,'\t')
                    SIGFILE.write('acc_num\t%s\n' % sheader)
                else:
                    SIGFILE = open(sigfile, 'a')

            ### <b> ### TMHMM
            for acc in self.tmhmm.keys():
                _stage = '<b> TMHMM'
                TMFILE.write('%s\t%s\t%s\t%s\n' % (acc,self.tmhmm[acc]['PredHel'],self.tmhmm[acc]['Topology'][0],self.tmhmm[acc]['Topology'][-1]))
                domains = self.tmhmm[acc]['Topology']
                tm = False
                dom = 'CYTOPLASMIC'
                if domains[0] == 'o':
                    dom = 'EXTRACELLULAR'
                domains = re.sub('o', '-', domains)
                domains = re.sub('i', '-', domains)
                domains = string.split('1' + domains + self.tmhmm[acc]['len'],'-')
                started = False
                while len(domains) > 1:
                    _stage = '<b-i> TM Dom Write'
                    start = domains.pop(0)
                    end = domains[0]
                    if tm:
                        type = 'TRANSMEMBRANE'
                    else:
                        type = dom
                        if started:
                            #print start,
                            start = '%d' % (string.atoi(start) + 1)
                            #print start
                        else:
                            started = True
                        if len(domains) > 1:
                            #print end,
                            end = '%d' % (string.atoi(end) - 1)
                            #print end
                    DOMFILE.write('%s\n' % string.join([acc,type,start,end,self.info['Source']],'\t'))
                    if tm:
                        tm = False
                        if dom == 'CYTOPLASMIC':
                            dom = 'EXTRACELLULAR'
                        else:
                            dom = 'CYTOPLASMIC'
                    else:
                        tm = True
                _stage = '<b-ii> TM Dom Check'
                if tm == False:
                    self.log.errorLog('Problem with %s TM domains - wrong number of domains!' % acc)
                
            sigdic = {'nn_cleavemax':'nn_ymax','nn_cleavepos':'nn_ymaxpos','nn_cleave':'nn_ymax?','nn_sig_mean':'nn_smean','nn_sig':'nn_smean?',
                      'nn_dscore':'nn_d','nn_d':'nn_d?','hmm_cmax':'hmm_cmax','hmm_cpos':'hmm_cmaxpos','hmm_cleave':'hmm_cmax?',
                      'hmm_sigprob':'hmm_sprob','hmm_sig':'hmm_sprob?'}

            ### <c> ### SignalP
            for acc in self.signalp.keys():
                _stage = '<c> SignalP'
                accout = acc
                if re.search('_HUMAN_(\S+)$', acc):
                    accout = rje.matchExp('_HUMAN_(\S+)$', acc)[0]
                writelist = [accout]
                for stat in siglist:
                    #print accout, stat
                    writelist.append(self.signalp[acc][sigdic[stat]])
                #print '%s\n' % string.join(writelist,'\t')
                SIGFILE.write('%s\n' % string.join(writelist,'\t'))
                _stage = '<c-ii> SingalP domains'
                nn_y = string.atoi(self.signalp[acc]['nn_ymaxpos']) - 1
                hmm_c = string.atoi(self.signalp[acc]['hmm_cmaxpos']) - 1
                if self.signalp[acc]['nn_d?'] == 'Y':
                    DOMFILE.write('%s\n' % string.join([accout,'SIGNALP','1','%d' % nn_y,'signalp-NN'],'\t'))
                if self.signalp[acc]['nn_ymax?'] == 'Y':
                    DOMFILE.write('%s\n' % string.join([accout,'CLEAVAGE','%d' % nn_y,'%d' % (nn_y+1),'signalp-NN'],'\t'))
                if self.signalp[acc]['hmm_sprob?'] == 'Y':
                    DOMFILE.write('%s\n' % string.join([accout,'SIGNALP','1','%d' % hmm_c,'signalp-HMM'],'\t'))
                if self.signalp[acc]['hmm_cmax?'] == 'Y':
                    DOMFILE.write('%s\n' % string.join([accout,'CLEAVAGE','%d' % hmm_c,'%d' % (hmm_c+1),'signalp-HMM'],'\t'))
                    
            ### <d> ### Finish
            _stage = '<d> Finish'
            if tmfile and self.tmhmm.keys():
                TMFILE.close()
            if domfile and (self.tmhmm.keys()+self.signalp.keys()):
                DOMFILE.close()
            if sigfile and self.signalp.keys():
                SIGFILE.close()
            return
        except:
            self.log.errorLog('Problem with mySQLOut() %s.' % _stage)
#########################################################################################################################
    ### <4> ### Masked Cleavage Output
#########################################################################################################################
    def maskCleave(self):   ### Outputs masked cleavage sequences to file
        '''
        Outputs masked cleavage sequences to file.
        '''
        try:
            ### <a> ### Setup
            _stage = '<a> Setup'
            seqlist = self.obj['SeqList']
            seqlist.loadSeqs()
            seqlist.degapSeq()
            outfile = '%s.cleaved.fas' % seqlist.info['Basefile']

            ### <b> ### MaskSeqs
            _stage = '<b> Mask'
            self.verbose(0,3,'Masking cleaved signalp petides from %s into %s...' % (seqlist.info['Name'],outfile),0)
            cx = 0
            for seq in seqlist.seq:
                cpos = 0
                acc = seq.info['AccNum']
                if acc in self.signalp.keys():
                    sigp = self.signalp.pop(acc)
                else:
                    continue
                if sigp['nn_ymax?'] == 'Y':
                    cpos = string.atoi(sigp['nn_ymaxpos']) 
                if sigp['hmm_cmax?'] == 'Y':
                    hmm_c = string.atoi(sigp['hmm_cmaxpos'])
                    if cpos==0 or (cpos > 0 and hmm_c < cpos):
                        cpos = hmm_c
                if cpos > 0:
                    cx += 1
                    seq.info['Sequence'] = 'X' * cpos + seq.info['Sequence'][cpos:]
                    if cx/100 == cx/100.0:
                        self.verbose(0,4,'.',0)
            self.verbose(0,1,'Done! %d sequences masked.' % cx,1)

            ### <c> ### Save
            _stage = '<c> Save'
            seqlist.saveFasta(seqfile=outfile)

            ### <d> ### Run TMHMM
            _stage = '<d> Run TMHMM'
            if rje.yesNo('Run TMHMM on this now?'):
                os.system('tmhmm %s.cleaved.fas -short > %s.cleaved.tmhmm' % (seqlist.info['Basefile'],seqlist.info['Basefile']))           
            return
        except:
            self.log.errorLog('Problem with maskCleave() %s.' % _stage)
#########################################################################################################################
## End of SECTION II: TM Class                                                                                          #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def parseTMHMM(tmline):     ### Returns a dictionary of TMHMM data from a TMHMM line
    '''Returns a dictionary of TMHMM data from a TMHMM line.'''
    tmdata = string.split(rje.chomp(tmline))
    tmdict = {'Seq':tmdata.pop(0)}
    for tm in tmdata:
        (tkey,tval) = string.split(tm,'=')
        tmdict[tkey] = tval
    return tmdict
#########################################################################################################################
def domainList(tmdict):   ### Returns list of TMHMM domain dictionaries {Type/Start/End} from tmdictionary
    '''Returns list of TMHMM domain dictionaries {Type/Start/End} from Topology string.'''
    ### ~ [1] Setup variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    domainlist = []
    tm = False
    dom = 'CYTOPLASMIC'
    domains = tmdict['Topology']
    started = False

    ### ~ [2] Adjust domain listing to include ends of sequences and orientation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if domains[0] == 'o': dom = 'EXTRACELLULAR'
    domains = re.sub('o', '-', domains)
    domains = re.sub('i', '-', domains)
    domains = string.split('1' + domains + tmdict['len'],'-')

    ### ~ [3] Parse domains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    while len(domains) > 1:
        start = domains.pop(0)
        end = domains[0]
        if tm: type = 'TRANSMEMBRANE'
        else:
            type = dom
            if started: start = '%d' % (string.atoi(start) + 1)
            else: started = True
            if len(domains) > 1: end = '%d' % (string.atoi(end) - 1)
        domainlist.append({'Type':type,'Start':start,'End':end})
        if tm:
            if dom == 'CYTOPLASMIC': dom = 'EXTRACELLULAR'
            else: dom = 'CYTOPLASMIC'
        tm = not tm
    if tm == False: raise ValueError
    if len(domainlist) == 1: domainlist[-1]['Type'] = 'NONMEMBRANE'
    return domainlist
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
    try:
        tm = TM(log=mainlog,cmd_list=cmd_list)
        tm.parseSignalP()
        tm.parseTMHMM()
        if tm.opt['MySQL']: tm.mySQLOut(makenew=not tm.opt['Append'])
        if tm.opt['MaskCleave']: tm.maskCleave()

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
