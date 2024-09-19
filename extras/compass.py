#!/usr/local/bin/python

# COMPASS - Comparison Of Motif Predictions Across Sequences/Servers
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       compass
Description:  Comparison Of Motif Predictions Across Sequences/Servers
Version:      2.2
Last Edit:    19/05/06
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Outputs a table of Scansite predictions etc as arranged by alignment > *.compass.tdt. This file should then be
    opened using the COMPASS.xls Excel workbook for analysis using its Macros.

    This is designed to compare multiple aligned sequences - typically orthologues - and/or multiple servers for making
    the same kinds of predictions. Although explicitly designed with certain servers (eg. scansite, netphos) in mind,
    any data that is presented in the right format can be uploaded and compared. The 'Server' column of the results takes
    its values from the extension of the input files. E.g. *.scansite will have the server 'scansite'.

Input files:
    - sequence alignment (may be one sequence)
    - *.server or *.server.txt server results 

Commandline:
    seqin=FILE      : Input alignment file
    partial=T/F     : Allow partial scansite results (NULL values) [False]
    autorun=T/F     : Automated querying of web servers (within RCSI only) [False]
    results=T/F     : Whether to output summary file or simply check/generate results (Aligned seqs only) [True]
    recase=T/F      : Whether to look for accession numbers in case-independent fashion (scansite results) [True]
    txt=T/F         : Whether to look by deafult for *.server.txt files (True) or *.server files (False) [True]

    scansite=FILE   : Motif prediction files for sequences in scansite format [*.scansite, *.automotif & *.elm]
    netphos=FILE    : Motif prediction files for sequences in netphos format [*.netphos]
    uniprot=FILE    : Sequence features/details in UniProt download format [*.uniprot]
    tmhmm=FILE      : Transmembrane topology data in TMHMM single-line format [*.tmhmm]
    signalp=FILE    : SignalP results in single-line format [*.signalp]
    serverlist=X,Y,..,Z : List of servers for results to be read from
                        ['scansite','netphos','automotif','elm','uniprot','tmhmm','signalp','antigenic','disorder']
    
Uses general modules: copy, os, string, sys, time
Uses RJE modules: rje, rje_seq, rje_tm, rje_uniprot
Other modules needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import os
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_seq
import rje_tm
import rje_uniprot
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 0.1 - NetPhos parsing
# 1.0 - Operating version with Minor changes.
# 1.1 - Added uniprot download parsing
# 1.2 - Renamed COMPASS from COMP and altered dictionary indexing. Added 'Server'. TMHMM and SignalP not functioning.
# 2.0 - Overhaul of Module Structure.
# 2.1 - Allow results files to be *.server.txt or just *.server. Added serverlist=X,Y,Z option
# 2.2 - Updated use of general methods, lists and dictionaries. Better error handling and log output.
#########################################################################################################################
### Major Functionality to Add
# [ ] Add Psipred/PDB/DSSP
# [ ] Add output from BADASP?
# [ ] Incorporate automated searching [autorun=T/F]
# [ ] Add reading of sequence groups if available?!
# [ ] Alter uniprot filtering to have one entry per residue rather than start & stop - make this an option [endsonly=T/F]
# [ ] Incorporate SLiMDisc results.
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'COMPASS'
        version = '2.2'
        last_edit = 'May 06'
        description = 'Comparison Of Motif Predictions Across Sequences/Servers'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print( 'Problem making Info object.')
        raise
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
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
    except: rje.printf('Problem during initial setup.'); raise
#############################################################################################################################
### END OF SECTION I
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Compass Class:                                                                                          #
#########################################################################################################################
class Compass(rje.RJE_Object):     
    '''
    Comparison Of Motif Prediction Across Sequences/Servers Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of sequence group for output (=rje_seq.SeqList.info['Basefile']
   
    Opt:boolean
    - Partial = Allow partial scansite results (NULL values) [False]
    - AutoRun = Automated querying of web servers (within RCSI only)
    - Results = Whether to output summary file or simply check/generate results (Aligned seqs only)
    - ReCase = Whether to look for accession numbers in case-independent fashion (scansite results)
    - Txt = Whether to look by default for *.server.txt files (rather than *.server)

    Stat:numeric

    List:list
    - Format = List of acceptable file formats
    - Server = List of servers that should have results.
    - File = List of file dictionaries
        # 'Name':str = filename
        # 'Server':str = server name
        # 'Format':str = format of data

    Dict:dictionary    
    - Format = Dictionary of default formats for each server
    - Feature = Dictionary of features.
        # Key(s) = [accnum:[pos:[server,group,feature]]].
        # Value = dictionary of results keys = [score,percentile,sa]

    Obj:RJE_Objects
    - SeqList = rje_seq.SeqList Object   
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Stats:float []
        - Opt:boolean ['Partial','AutoRun','Results','ReCase','Txt']
        - List:list ['Format','Server','File']
        - Dict:dictionary ['Format','Feature']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name']
        self.statlist = []
        self.optlist = ['Partial','AutoRun','Results','ReCase','Txt']
        self.listlist = ['Format','Server','File']
        self.dictlist = ['Format','Feature']
        self.objlist = ['SeqList']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.opt['Results'] = True
        self.opt['ReCase'] = True
        self.opt['Txt'] = True
        ### <c> ### Other Attributes
        self.dict['Feature'] = {}
        self.list['Format'] = ['scansite','netphos','uniprot','tmhmm','signalp']
        self.list['Server'] = ['scansite','netphos','automotif','elm','uniprot','tmhmm','signalp','antigenic','disorder']
        self.dict['Format'] = {'scansite':'scansite',
                          'netphos':'netphos',
                          'automotif':'scansite',
                          'elm':'uniprot',
                          'uniprot':'uniprot',
                          'antigenic':'uniprot',
                          'disorder':'uniprot',
                          'tmhmm':'tmhmm',
                          'signalp':'signalp'}
        self.list['File'] = []
#########################################################################################################################
    def setSeqList(self,seqlist):   ### Sets given seqlist as self.obj['SeqList']
        '''Sets given seqlist as self.obj['SeqList'].'''
        try:
            _stage = '<1> Check Aln'
            self.obj['SeqList'] = seqlist
            self.info['Name'] = seqlist.info['Basefile']
            if self.opt['Results'] and self.stat['Interactive'] < 0:
                seqlist._checkAln(aln=True,realign=True)
                return
            else:
                seqlist._checkAln()

            _stage = '<2> Realign Options'
            if self.opt['Results'] and seqlist.opt['Aligned'] == False:
                if rje.yesNo('Sequences are not aligned but Summary Results output asked for. Realign?'):
                    seqlist.muscleAln(outfile=None,mapseq=True)
                    seqlist._checkAln(aln=True)
                else:
                    self.log.printLog('#ALN','Unaligned sequences: no Summary Results table.')
                    self.opt['Results'] = False
        except:
            self.log.errorLog('Problem with setSeqList(%s)' % _stage)
            raise
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
                self._cmdRead(type='opt',att='Partial',cmd=cmd)
                self._cmdRead(type='opt',att='AutoRun',cmd=cmd)
                self._cmdRead(type='opt',att='Results',cmd=cmd)
                self._cmdRead(type='opt',att='Txt',cmd=cmd)
                self._cmdRead(cmd,type='list',att='Server',arg='serverlist')
                #if rje.matchExp('^serverlist=(\S+)',cmd):
                #    self.list['Server'] = rje.split(rje.matchExp('^serverlist=(\S+)',cmd)[0],sep=',')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def setFileList(self):   ### Sets filelist dictionaries and warns of duplicates etc.
        '''Sets filelist dictionaries and warns of duplicates etc.'''
        try:
            ### Setup ###
            _stage = 'Setup'
            self.list['File'] = []
            servercount = {}
            for server in self.list['Server']:
                servercount[server] = 0

            ### Read File Names etc. ###
            _stage = 'Read File Names'
            self.verbose(0,2,'Reading input file names and servers...',1)
            for cmd in self.cmd_list:
                eq = cmd.find('=')
                if eq > 0 and cmd[:eq] in self.list['Format']:  # Command giving filename
                    _format = cmd[:eq]
                    _filename = cmd[eq+1:]
                    _server = _filename
                    if _server[-4:] == '.txt':
                       _server = _server[:-4] 
                    if _server.find('.') >= 0:
                        _server = rje.matchExp('\.([A-Za-z0-9_\-]+)$',_server)[0]
                    self.list['File'].append({'Name':_filename,'Server':_server,'Format':_format})
                    if _server in servercount.keys():
                        servercount[_server] += 1
                    else:
                        servercount[_server] = 1
                
            ### Check List and Add Defaults ###
            _stage = 'Check List'
            _addfiles = []
            for server in self.list['Server']:
                if servercount[server] == 0:
                    _filename = '%s.%s' % (self.info['Name'],server)
                    if self.opt['Txt']:
                        _filename = _filename + '.txt'
                    _addfiles.append({'Name':_filename,'Server':server,'Format':self.dict['Format'][server]})
                    servercount[server] += 1
                if servercount[server] > 1:
                    self.verbose(0,2,'WARNING: %d files for %s server. Later results will overwrite earlier duplicates.' % (servercount[server],server),1)
            self.list['File'] = _addfiles + self.list['File']

            ### Summarise ###
            _stage = 'Summarise'
            for file in self.list['File'][0:]:
                self.verbose(0,4,'- %s = %s server results in %s format.' % (file['Name'],file['Server'],file['Format']),1)
            self.verbose(0,1,'=> %d results files for %d servers.' % (len(self.list['File']),len(servercount)),2)

            ### Additional Options ###
            if self.stat['Interactive'] > 0 and not rje.yesNo('Use all %d files?' % len(self.list['File'])):
                for file in self.list['File'][0:]:
                    if not rje.yesNo('Use %s [%s server results in %s format]?' % (file['Name'],file['Server'],file['Format'])):
                        if rje.yesNo('Discard %s. Are you sure?' % file['Name']):
                            self.list['File'].remove(file)
            #!# Add option here at some point to add remove files to list?

        except:
            self.log.errorLog('Major Problem with setFileList(%s).' % _stage)
#########################################################################################################################
    def _updateResults(self,file,seqfilebase,append,missing,missingseq):    ### Prints prompt to update results and/or calls automated script.
        '''
        Prints prompt to update results and/or calls automated script.
        >> file:dictionary of results file details
        >> seqfilebase:str = base for name of sequence file for server upload
        >> append:str = whether appending or saving as...
        >> missing:list of missing AccNum
        >> missingseq:list of missing Sequence Objects
        << True (continue) or False
        '''
        try:
            _stage = '<1> Automated Updates'
            if file['Server'] in self.list['Server'] and self.opt['AutoRun']:
                print( 'AutoRun feature not yet implemented. Sorry!')
                if self.stat['Interactive'] < 0:
                    return False
            seqfile = None

            _stage = '<2> Save Missing Sequences'
            print( '\nNOTE: %d of %d %s results appear to be missing or in wrong format.' % (len(missing),self.obj['SeqList'].seqNum(),file['Server']))
            self.verbose(1,4,'Missing: %s' % missing,1)
            seqlist = self.obj['SeqList']
            if file['Server'] in ['scansite','automotif','elm']: ## AccNum and sequence on one line
                seqfile = '%s.scanseq' % seqfilebase
                seqlist.saveScanSeq(seqs=missingseq,seqfile=seqfile)
            elif file['Server'] in ['netphos']: ## All sequences, fasta format with numbers
                seqfile = '%s.netphos.fas' % self.info['Name']
                append = 'save as'
                missingseq = seqlist.seq
                seqlist.saveFasta(seqfile=seqfile,name='AccNum',id=True)
            elif file['Server'] in ['tmhmm']: ## Fasta format, accession numbers only
                seqfile = '%s.acc.fas' % seqfilebase
                seqlist.saveFasta(seqs=missingseq,seqfile=seqfile,name='AccNum')
            elif file['Server'] in ['signalp']: ## Truncated seqs in Fasta format, accession numbers only
                seqfile = '%s.signalp.fas' % seqfilebase
                sigseq = copy.deepcopy(seqlist)
                sigmiss = copy.deepcopy(missingseq)
                sigseq.seq = sigmiss
                sigseq.truncSeq(70)
                sigseq.saveFasta(seqfile=seqfile,name='AccNum')
            elif file['Server'] in ['uniprot']: ## Accession numbers only
                seqfile = '%s.acc' % seqfilebase
                seqlist.saveAcc(seqs=missingseq,accfile=seqfile,uniprot=True)
            if seqfile:
                print( ' => Missing sequences saved to %s.' % seqfile)
            else:
                print( ' => Unknown server details. Cannot output relevant file format for upload.')

            _stage = '<3> Upload instructions'
            upload = 'Unknown Server'
            if file['Server'] == 'scansite':
                upload = 'Scansite Parallel version as a File of Protein Sequences (not IDs):\nhttp://stjuderesearch.org/scansite/ (Low Stringency setting)'
            elif file['Server'] == 'netphos':
                upload = 'NetPhos:\nhttp://www.cbs.dtu.dk/services/NetPhos/ (no graphics)'               
            elif file['Server'] == 'uniprot':
                upload = 'UniProt:\nhttp://us.expasy.org/sprot/sprot-retrieve-list.html'
            elif file['Server'] == 'tmhmm':
                upload = 'TMHMM (one line per protein):\nhttp://www.cbs.dtu.dk/services/TMHMM/'
            elif file['Server'] == 'signalp':
                upload = 'SignalP (Short output (no graphics!)):\nhttp://www.cbs.dtu.dk/services/SignalP/'
                upload = upload + '\nNB. Sequences have been truncated to 70aa.'
            elif file['Server'] == 'automotif':
                upload = 'AutoMotif using Python script:\npython queryAutoMotif.py %s' % seqfile
            elif file['Server'] == 'elm':
                upload = 'ELM using Python script:\nCurrently unavailable, sorry!' #!# Get script from Norman!
            print( '\n1. Upload file to %s' % upload)
            print( '\n2. When results are ready, %s %s.' % (append, file['Name']))
            
            _stage = '<4> User decision'
            if rje.yesNo('\n3. Filed updated? Yes to continue, No to bypass or give alternative file.'):

                return True
            else:
                newfile = rje.choice('Enter new filename or blank to bypass: ')

                if newfile == '':
                    return False
                else:
                    file['Name'] = newfile
                    return True
        except KeyboardInterrupt:
            raise
        except:
            self.log.errorLog('Major Problem with Compass._updateResults(%s).' % _stage)
            raise
#########################################################################################################################
    ### <2> ### Make feature dictionary etc.                                                                            #
#########################################################################################################################
    def run(self,clear=True):  ### Main method calling primary COMPASS methods.
        '''
        Main method calling primary COMPASS methods.
        - setSeqList()
        - setFileList()
        - readFileData()
        - compResults()
        > clear:Boolean [True] = whether to clear self.dict['Feature'] before reading results.
        '''
        try:
            ### SeqList setup ###
            _stage = 'SeqList Setup'
            if self.obj['SeqList'] == None:
                seqlist = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list)
                if seqlist.seqNum() < 1:
                    self.log.errorLog('Must input sequence alignment using seqin=FILE',printerror=False,quitchoice=False)
                    raise ValueError
                self.setSeqList(seqlist)

            ### FileList Setup ###
            _stage = 'FileList Setup'
            self.setFileList()

            ### Read results into self.dict['Feature'] ###
            _stage = 'Read Results'
            if clear:
                self.dict['Feature'].clear()
            for file in self.list['File']:
                _stage = 'General Read Results'
                #filename = file['Name']
                #server = file['Server']
                #format = file['Format']
                print( '\nLoad %s server results from %s (%s format)...' % (file['Server'],file['Name'],file['Format']))
                while 1:
                    _stage = '<2a-ii> General Read Results'
                    missing = self.readFileData(filedic=file)
                    #print missing
                    if len(missing) == 0 or (self.stat['Interactive'] < 0 and self.opt['AutoRun'] == False):   # No Seqs missing or ignoring
                        break
                    elif self.opt['Partial'] and len(missing) < seqlist.seqNum(): # Missing some only - OK to continue
                        break
                    else:
                        if rje.checkForFile(file['Name']):  # File present but results missing - append results
                            append = 'append to'
                        else:   # No results file! Make new file.
                            append = 'save as'
                        if len(missing) < seqlist.seqNum(): # Missing some seqs only
                            seqfile = '%s.missing' % self.info['Name']
                            missingseq = []
                            for seq in seqlist.seq:
                                if seq.info['AccNum'] in missing:
                                    missingseq.append(seq)
                        else:   # Missing all seqs!
                            seqfile = self.info['Name']
                            missingseq = seqlist.seq
                        _stage = '<2a-iii> Update details'
                        if self._updateResults(file,seqfile,append,missing,missingseq):
                            continue
                        else:
                            break
    
            ### <3> ### Results Output
            _stage = '<3> Results Output'
            if self.opt['Results'] and self.dict['Feature']:
                self._compResults()
            elif self.opt['Results']:
                self.log.errorLog('No results read in! No output generated.',printerror=False,quitchoice=False)
            else:
                self.log.printLog('#INF','No output generated.')
        except KeyboardInterrupt:
            raise
        except SystemExit:
            raise
        except:
            self.log.errorLog('Major Problem with COMPASS run(%s).' % _stage,quitchoice=True)
#########################################################################################################################
    def readFileData(self,filedic):  ### Method that looks for files, prompts for creation and calls relevant processing.
        '''
        Method that looks for files, prompts for creation and calls relevant processing.
        >> filedic:dictionary object of Name, Server and Format.
        << Returns list of missing AccNum. This is empty if all present.

        This is where additional default servers should be added. Each server output format should have its own method.
        These methods, e.g. _processScansite() should return a list of missing accession numbers (empty if all found).
        The action taken when sequences are missing will then depend on the 'Partial' setting and the interactivity.
        Saving of sequences in the correct format for server upload is handled by methods in the SeqList class.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            filename = filedic['Name']
            server = filedic['Server']
            format = filedic['Format']
            acclist = self.obj['SeqList'].accList()
            
            ### Check for File ###
            _stage = 'Check for File'
            reslines = []
            if not rje.checkForFile(filename):
                self.log.errorLog('%s results file %s missing!' % (server, filename),printerror=False,quitchoice=False)
                return acclist
                
            ### Process Results ###
            _stage = 'Process Results'
            if format == 'scansite':
                return self._processScansite(filedic)
            elif format == 'netphos':
                return self._processNetphos(filedic)
            elif format == 'uniprot':
                return self._processUniprot(filedic)
            elif format == 'tmhmm':
                return self._processTMHMM(filedic)
            elif format == 'signalp':
                return self._processSignalP(filedic)
            else:
                self.log.errorLog('%s results file %s in unrecognised format (%s)!' % (server, filename, format),printerror=False,quitchoice=False)
                return []

        except SystemExit:
            raise
        except:
            self.log.errorLog('Major Problem with readFileData(%s).' % _stage)
#########################################################################################################################
    def _addFeature(self,accnum,pos,enz_info,feature):  ### Adds a feature to self.dict['Feature']
        '''
        Adds a feature to self.dict['Feature'].
        >> accnum:str = Accession number of sequence
        >> pos:str = position in sequence as string
        >> enz_info:tuple of (server,group,feature)
        >> feature:dictionary of {'score':score,'perc':percentile,'sa':sa,'aa':aa}
        << Returns 0 if problem and 1 if OK.
        '''
        try:
            if accnum not in self.dict['Feature'].keys():
                self.dict['Feature'][accnum] = {}
            if pos not in self.dict['Feature'][accnum].keys():
                self.dict['Feature'][accnum][pos] = {}
                #self.dict['Feature'][accnum][pos][enz_info] = {}
            self.dict['Feature'][accnum][pos][enz_info] = feature
            return 1
        except:
            self.log.errorLog('Problem with adding feature %s:%s %s %s.' % (accnum, pos, enz_info, feature),quitchoice=True)
            return 0
#########################################################################################################################
    def _processScansite(self,filedic):  ### Populates dictionaries
        '''
        Populates feature dictionary with scansite format results.
        >> filedic:dictionary of file info: Name, Format, Server
        << Returns list of missing AccNum. This is empty if all present
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            seq_accs = self.obj['SeqList'].accList()
            missing = seq_accs[0:]
            server = filedic['Server']
            
            ### Read lines ###
            _stage = 'Read'
            RES = open(filedic['Name'],'r')
            _linex = 0
            _ftx = 0
            prevacc = None
            while 1:
                line = RES.readline()
                if line in ['',None]:
                    break
                scan = rje.matchExp('^(\S+)\s+(\S+)\s+(\S+)\s+(\S)(\d+)\s+(\S+)\s+(\S+\%)\s+(\S+)\s+(\S+)\s*$',line)   # Results line
                ptuple = (rje.integerString(_linex),server,rje.integerString(_ftx),len(seq_accs)-len(missing),len(missing))
                self.log.printLog('\r#SCAN','%s %s lines read. %s features for %d sequences. %d missing.' % ptuple,newline=False,log=False)
                _linex += 1
                if scan:
                    [accnum, enzyme, group, aa, pos, score, percentile, sequence, sa] = scan
                    enz_info = (server,group,enzyme)
                    accnum = self._reCase(accnum,seq_accs)
                    if prevacc != accnum and not missing:   # Finished reading!
                        break
                    prevacc = accnum
                    if accnum not in seq_accs:
                        continue
                    if accnum in missing:
                        missing.remove(accnum)
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _ftx += self._addFeature(accnum,pos,enz_info,feature)
            RES.close()
            ptuple = (rje.integerString(_linex),server,rje.integerString(_ftx),len(seq_accs)-len(missing),len(missing))
            self.log.printLog('\r#SCAN','%s %s lines read. %s features for %d sequences. %d missing.' % ptuple)
                    
            ### Return Missing ###
            return missing       
        except:
            self.log.errorLog('Major Problem with Compass._processScansite(%s).' % _stage)
            return seq_accs
#########################################################################################################################
    def _reCase(self,accnum,seq_accs):  ### Looks for case-altered accnum in list and returns 'true' accnum
        '''
        Looks for case-altered accnum in list and returns 'true' accnum.
        >> accnum:str = Accession number as returned by website results
        >> seq_accs:list of real accnums
        '''
        try:
            if accnum in seq_accs or self.opt['ReCase'] == False:
                return accnum
            for acc in seq_accs:
                if acc.lower() == accnum.lower():
                    return acc
            return accnum
        except:
            self.log.errorLog('Major Problem with _reCase(%s).' % accnum)
            return accnum
#########################################################################################################################
    def _processNetphos(self,filedic):  ### Populates dictionaries
        '''
        Populates feature dictionary with netphos format results.
        >> filedic:dictionary of file info: Name, Format, Server
        << Returns list of missing AccNum. This is empty if all present
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            seq_accs = self.obj['SeqList'].accList()
            missing = seq_accs[0:]
            server = filedic['Server']
            _numbered = False
            
            ### Read lines ###
            _stage = 'Read'
            RES = open(filedic['Name'],'r')
            _linex = 0
            _ftx = 0
            while 1:
                line = RES.readline()
                if line in ['',None]:
                    break
                scan = rje.matchExp('^(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\*(\S)\*\s*$',line)   # Positive Results line
                if scan:
                    [accnum, pos, sequence, score, aa] = scan
                else:
                    scan = rje.matchExp('^(\S+)\s+(\d+)\s+(\S)\s+(\S+)\s+(\d\S+)',line)     # Pre-processed results
                    if scan:
                        [accnum, pos, aa, sequence, score] = scan
                # ENSCAFP0000      6   AKQPSDVSS  0.941  *S*
                ptuple = (rje.integerString(_linex),server,rje.integerString(_ftx),len(seq_accs)-len(missing),len(missing))
                self.log.printLog('\r#NET','%s %s lines read. %s features for %d sequences. %d missing.' % ptuple,newline=False,log=False)
                _linex += 1
                if scan:
                    enzyme = 'NP_%s_kin' % aa
                    group = '%s_kin' % aa
                    enz_info = (server,group,enzyme)
                    sa = 'NULL'
                    percentile = 1 - string.atof(score)
                    percentile = '%.3f' % percentile
                    if accnum not in seq_accs and rje.matchExp('^(\d+)$',accnum):
                        _numbered = True
                        accnum = seq_accs[string.atoi(accnum)-1]
                    if accnum not in seq_accs:
                        if _numbered:
                            return seq_accs    # Should have results for all seqs and seqs only if numbers given
                        else:
                            continue
                    if accnum in missing:
                        missing.remove(accnum)
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _ftx += self._addFeature(accnum,pos,enz_info,feature)
                #X#else:
                #X#    raw_input('\nCannot read:\n%s <ENTER>' % line)
            RES.close()
            ptuple = (rje.integerString(_linex),server,rje.integerString(_ftx),len(seq_accs)-len(missing),len(missing))
            self.log.printLog('\r#NET','%s %s lines read. %s features for %d sequences. %d missing.' % ptuple)

            ### <2> ### Re-save NetPhos if desirable        #!# Not updated yet #!#
            _stage = '<2> Resave Results'
            if _numbered and (self.stat['Interactive'] < 0 or rje.yesNo('Re-save NetPhos results with Accession Numbers?')):
                try:
                    try:
                        os.unlink('%s.bak' % filedic['Name'])
                        os.rename(filedic['Name'],'%s.bak' % filedic['Name'])
                        self.log.printLog('#INF','NetPhos file %s backed up to %s.bak' % (filedic['Name'],filedic['Name']))
                    except:
                        self.log.errorLog('Problem generating backup file, %s.bak.' % filedic['Name'])
                        if rje.yesNo('Proceed without backup?') == False:
                            raise

                    NP = open(filedic['Name'],'w')
                    _writex = 0
                    for line in reslines:
                        scan = rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$',line)   # Positive Results line
                        if scan:
                            [accnum, pos, sequence, score, aa] = scan
                            accnum = seq_accs[string.atoi(accnum)-1]
                            if aa == '.':
                                aa = sequence[4]
                            NP.write('%s\n' % rje.join([accnum, pos, sequence, score, aa],'\t'))
                            _writex += 1
                    NP.close()
                    self.log.printLog('#INF','%d reformatted result lines saved to %s.' % (_writex, filedic['Name']))
                except:
                    self.log.errorLog('Problem generating new results file, %s.' % filedic['Name'])
                    
            ### Return Missing ###
            return missing              
        except:
            self.log.errorLog('Major Problem with Compass._processNetPhos(%s).' % _stage)
            return seq_accs
#########################################################################################################################
    def _processUniprot(self,filedic):  ### Populates dictionaries
        '''
        Populates feature dictionary with uniprot format results.
        >> filedic:dictionary of file info: Name, Format, Server
        << Returns list of missing AccNum. This is empty if all present
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            seq_accs = self.obj['SeqList'].accList()    #!# Strip splice variants #!#
            missing = seq_accs[0:]
            server = filedic['Server']
            uniprot = rje_uniprot.UniProt(log=self.log,cmd_list=self.cmd_list)
            if not uniprot.readUniProt(filename=filedic['Name'],clear=True,acclist=seq_accs,use_index=False):
                uniprot.entry = []
                return seq_accs
            #self.log.printLog('#INF','%d Uniprot-style entries read in from %s.' % (len(uniprot.entry),filedic['Name']))

            ### <1> ### Read features
            _featurex = 0
            _linex = 0
            ### <1> ### Read Features
            _stage = '<1> Read features'
            for entry in uniprot.list['Entry']:
                _stage = '<1a> AccNum'
                accnum = entry.obj['Sequence'].info['AccNum']
                if accnum not in seq_accs:
                    continue    # Should have results for seqs only
                if accnum in missing:   # In results file (but may have no features)
                    missing.remove(accnum)
                featurelist = []
                _featurex += len(entry.list['Feature'])
                
                _stage = '<1b> FeatureList'
                for feature in entry.list['Feature']:
                    #print feature
                    if feature['Start'] == feature['End']:
                        featurelist.append(copy.deepcopy(feature))
                        featurelist[-1]['AA'] = entry.obj['Sequence'].info['Sequence'][featurelist[-1]['Start']-1]
                    else:
                        featurelist.append(copy.deepcopy(feature))
                        if feature['Type'] in ['VARSPLIC','CONFLICT']:
                            featurelist[-1]['Desc'] = '%s Start.' % feature['Type'].lower()
                        else:
                            featurelist[-1]['Desc'] = '%s Start' % featurelist[-1]['Desc']
                        featurelist[-1]['AA'] = entry.obj['Sequence'].info['Sequence'][featurelist[-1]['Start']-1]
                        featurelist.append(copy.deepcopy(feature))
                        if feature['Type'] in ['VARSPLIC','CONFLICT']:
                            featurelist[-1]['Desc'] = '%s End.' % feature['Type'].lower()
                        else:
                            featurelist[-1]['Desc'] = '%s End' % featurelist[-1]['Desc']
                        featurelist[-1]['AA'] = entry.obj['Sequence'].info['Sequence'][featurelist[-1]['End']-1]
                        featurelist[-1]['Start'] = featurelist[-1]['End']

                _stage = '<1c> Convert to Scansite'
                for feature in featurelist:                        
                    group = feature['Type']
                    enzyme = feature['Desc']
                    enz_info = (server,group,enzyme)
                    pos = '%d' % feature['Start']
                    sequence = feature['AA']
                    aa = sequence
                    score = '1.0'
                    if enzyme.lower().find('by similarity') >= 0 or server not in ['uniprot']:
                        percentile = '100.000%'
                    else:
                        percentile = '0.000%'
                    sa = 'NULL'
                    details = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,details)
                    
            ### <2> ### Return Missing
            _stage = '<2> Return Missing'
            self.log.printLog('#INF','%d %s results (%d features) read in for %d sequences. %d missing.' % (_linex,server,_featurex,self.obj['SeqList'].seqNum(),len(missing)))
            return missing       
        except:
            self.log.errorLog('Major Problem with _processUniprot(%s).' % _stage)
            return seq_accs
#########################################################################################################################
    def _processTMHMM(self,filedic):  ### Populates dictionaries
        '''
        Populates feature dictionary with tmhmm format results.
        >> filedic:dictionary of file info: Name, Format, Server
        << Returns list of missing AccNum. This is empty if all present
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            seq_accs = self.obj['SeqList'].accList()
            missing = seq_accs[0:]
            server = filedic['Server']
            tm = rje_tm.TM(log=self.log,cmd_list=[])
            tm.parseTMHMM(file=filedic['Name'])
            # Now have topologies in tm.tmhmm[accnum]['topology']
            
            ### <1> ### Read Domains
            _featurex = 0
            _linex = 0
            _stage = '<1> Read domains'
            for key in tm.tmhmm.keys():
                _stage = '<1a> AccNum'
                accnum = key
                if accnum not in seq_accs:
                    if accnum.find('_') >= 0:
                        accnum = rje.matchExp('_([A-Za-z0-9\.\-]+)$',accnum)[0]
                    if accnum not in seq_accs:
                        continue    # Should have results for seqs only
                if accnum in missing:
                    missing.remove(accnum)
                    
                _stage = '<1b> Convert to Scansite'
                counts = {}     # Counts of different domains
                #print key, tm.getDomains(key)
                for domain in tm.getDomains(key):
                    _featurex += 1
                    #print domain
                    if domain['Type'] == 'TRANSMEMBRANE':
                        group = 'TRANSMEM'
                    else:
                        group = 'TOP_DOM'
                    sequence = 'NULL'
                    aa = sequence
                    score = '1.0'
                    percentile = '100.000%'
                    sa = 'NULL'
                    if domain['Type'] in counts.keys():
                        counts[domain['Type']] += 1
                    else:
                        counts[domain['Type']] = 1

                    # Start
                    _stage = '<1c> Start as UniProt Feature'
                    enzyme = '%s %d. Start.' % (domain['Type'],counts[domain['Type']])
                    enz_info = (server,group,enzyme)
                    pos = domain['Start']
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)

                    # End
                    _stage = '<1d> End as UniProt Feature'
                    enzyme = '%s %d. End.' % (domain['Type'],counts[domain['Type']])
                    enz_info = (server,group,enzyme)
                    pos = domain['End']
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)
                    
            ### <2> ### Return Missing
            _stage = '<2> Return Missing'
            self.log.printLog('#INF','%d %s results (%d features) read in for %d sequences. %d missing.' % (_linex,server,_featurex,self.obj['SeqList'].seqNum(),len(missing)))
            return missing       
        except:
            self.log.errorLog('Major Problem with _processTMHMM(%s).' % _stage)
            return seq_accs
#########################################################################################################################
    def _processSignalP(self,filedic):  ### Populates dictionaries
        '''
        Populates feature dictionary with signalp format results.
        >> filedic:dictionary of file info: Name, Format, Server
        << Returns list of missing AccNum. This is empty if all present
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            seq_accs = self.obj['SeqList'].accList()
            missing = seq_accs[0:]
            server = filedic['Server']
            tm = rje_tm.TM(log=self.log,cmd_list=[])
            tm.parseSignalP(file=filedic['Name'])
            # Now have topologies in tm.signalp[accnum]['topology']
            
            ### <1> ### Read Features
            _linex = 0
            _stage = '<1> Read features'
            for key in tm.signalp.keys():
                _stage = '<1a> AccNum'
                accnum = key
                #print accnum, seq_accs
                if accnum not in seq_accs:
                    if accnum.find('_') >= 0:
                        accnum = rje.matchExp('_([A-Za-z0-9\.\-]+)$',accnum)[0]
                    #print accnum, seq_accs
                    if accnum not in seq_accs:
                        continue    # Should have results for seqs only
                if accnum in missing:
                    missing.remove(accnum)
                
                _stage = '<1b> Convert to Scansite'
                [sequence, aa, sa] = ['NULL'] * 3
                nn_y = string.atoi(tm.signalp[key]['nn_ymaxpos']) - 1
                hmm_c = string.atoi(tm.signalp[key]['hmm_cmaxpos']) - 1
                hmm_p = string.atof(tm.signalp[key]['hmm_sprob'])

                _stage = '<1b-i> NN Sig'
                if tm.signalp[key]['nn_d?'] == 'Y':
                    score = tm.signalp[key]['nn_d']
                    percentile = '100.000%'
                    group = 'SIGNALP'
                    enzyme = 'signalp-NN End.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % nn_y
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)

                _stage = '<1b-ii> NN Cleave'
                if tm.signalp[key]['nn_ymax?'] == 'Y':
                    score = tm.signalp[key]['nn_ymax']
                    percentile = '100.000%'
                    group = 'CLEAVAGE'
                    enzyme = 'signalp-NN Cleavage.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)

                _stage = '<1b-iii> HMM Sig'
                if tm.signalp[key]['hmm_sprob?'] == 'Y':
                    score = tm.signalp[key]['hmm_sprob']
                    percentile = '%f%%' % (100.0 * (1 - hmm_p))
                    group = 'SIGNALP'
                    enzyme = 'signalp-HMM End.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)

                _stage = '<1b-iv> HMM Cleave'
                if tm.signalp[key]['hmm_cmax?'] == 'Y':
                    score = tm.signalp[key]['hmm_cmax']
                    percentile = '%f%%' % (100.0 * (1 - hmm_p))
                    group = 'CLEAVAGE'
                    enzyme = 'signalp-HMM Cleavage.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    feature = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    _linex += self._addFeature(accnum,pos,enz_info,feature)
                    
            ### <2> ### Return Missing
            _stage = '<2> Return Missing'
            self.log.printLog('#INF','%d %s features read in for %d sequences. %d missing.' % (_linex,server,self.obj['SeqList'].seqNum(),len(missing)))
            return missing       
        except:
            self.log.errorLog('Major Problem with _processSignalP(%s).' % _stage)
            return seq_accs
#########################################################################################################################
    def _sigDic(self,sigfile=None,tm=None,check=True):  ### Populates dictionaries
        '''
        Populates dictionaries.
        >> sigfile:str = File from which SignalP-format results are to be loaded.
        >> tm:rje_tm.TM Object
        >> check:boolean = Whether to perform checks on data integrity. [True]
        << Returns True if all AccNum present or False if not.
        '''
        try:
 
            ### <1> ### Read Features
            _featurex = 0
            _stage = '<1> Read features'
            for key in tm.tmhmm.keys():
                _stage = '<1a> AccNum'
                accnum = key
                if accnum not in seq_accs:
                    if accnum.find('_') >= 0:
                        accnum = rje.matchExp('_([A-Za-z0-9\.\-]+)$',accnum)[0]
                    if accnum not in seq_accs:
                        continue    # Should have results for seqs only
                
                _stage = '<1b> Convert to Scansite'
                sequence = 'NULL'
                aa = sequence
                sa = 'NULL'
                    
                nn_y = string.atoi(tm.signalp[accnum]['nn_ymaxpos']) - 1
                hmm_c = string.atoi(tm.signalp[accnum]['hmm_cmaxpos']) - 1
                if tm.signalp[accnum]['nn_d?'] == 'Y':
                    score = tm.signalp[accnum]['nn_d']
                    percentile = '100.000%'
                    group = 'SIGNALP'
                    enzyme = 'signalp-NN End.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % nn_y
                    #!# This section is now the same for all and can be made a method?
                    if accnum not in self.dict['Feature'].keys():
                        self.dict['Feature'][accnum] = {}
                    if pos not in self.dict['Feature'][accnum].keys():
                        self.dict['Feature'][accnum][pos] = {}
                        self.dict['Feature'][accnum][pos][enz_info] = {}
                    self.dict['Feature'][accnum][pos][enz_info] = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    #!#

                if self.signalp[accnum]['nn_ymax?'] == 'Y':
                    score = tm.signalp[accnum]['nn_ymax']
                    percentile = '100.000%'
                    group = 'CLEAVAGE'
                    enzyme = 'signalp-NN Cleavage.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    #!# This section is now the same for all and can be made a method?
                    if accnum not in self.dict['Feature'].keys():
                        self.dict['Feature'][accnum] = {}
                    if pos not in self.dict['Feature'][accnum].keys():
                        self.dict['Feature'][accnum][pos] = {}
                        self.dict['Feature'][accnum][pos][enz_info] = {}
                    self.dict['Feature'][accnum][pos][enz_info] = {'score':score,'perc':percentile,'sa':sa,'aa':aa}

                if self.signalp[accnum]['hmm_sprob?'] == 'Y':
                    score = tm.signalp[accnum]['hmm_sprob']
                    percentile = '%f%%' % (100.0 * (1 - tm.signalp[accnum]['hmm_sprob']))
                    group = 'SIGNALP'
                    enzyme = 'signalp-HMM End.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    #!# This section is now the same for all and can be made a method?
                    if accnum not in self.dict['Feature'].keys():
                        self.dict['Feature'][accnum] = {}
                    if pos not in self.dict['Feature'][accnum].keys():
                        self.dict['Feature'][accnum][pos] = {}
                        self.dict['Feature'][accnum][pos][enz_info] = {}
                    self.dict['Feature'][accnum][pos][enz_info] = {'score':score,'perc':percentile,'sa':sa,'aa':aa}

                if self.signalp[accnum]['hmm_cmax?'] == 'Y':
                    score = tm.signalp[accnum]['hmm_cmax']
                    percentile = '%f%%' % (100.0 * (1 - tm.signalp[accnum]['hmm_sprob']))
                    group = 'CLEAVAGE'
                    enzyme = 'signalp-HMM Cleavage.'
                    enz_info = (server,group,enzyme)
                    pos = '%d' % (nn_y + 1)
                    #!# This section is now the same for all and can be made a method?
                    if accnum not in self.dict['Feature'].keys():
                        self.dict['Feature'][accnum] = {}
                    if pos not in self.dict['Feature'][accnum].keys():
                        self.dict['Feature'][accnum][pos] = {}
                        self.dict['Feature'][accnum][pos][enz_info] = {}
                    self.dict['Feature'][accnum][pos][enz_info] = {'score':score,'perc':percentile,'sa':sa,'aa':aa}
                    
            ### <2> ### Check
            _stage = '<2> Check'    #!# This checking routine needs improvement - it only works for first of type
            self.log.printLog('#INF','%d UniProt features match given accession numbers.' % _featurex)
            _missing = 0
            for seq in self.obj['SeqList'].seq:
                if seq.info['AccNum'] not in self.dict['Feature'].keys():
                    self.log.printLog('#ERR','No results for %s: Missing?!' % seq.info['AccNum'])
                    _missing += 1
            self.log.printLog('#INF','Results read in for %d sequences. %d missing.' % (self.obj['SeqList'].seqNum(),_missing))
            if check and _missing > 0:
                self.log.printLog('#ERR','WARNING: %d Missing!' % _missing)
                return False
            return True
        
        except:
            self.log.errorLog('Major Problem with sigDic(%s).' % _stage)
#########################################################################################################################
    ### <3> ### Generation and Output of Results                                                                        #
#########################################################################################################################
    def _compResults(self):  ### Populates objects as necessary and executes methods
        '''Populates objects as necessary and executes methods.'''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            seqpos = {}     # Dictionary of {seq:position}
            myscan = copy.deepcopy(self.dict['Feature'])   # Local copy of feature
            seqlist = self.obj['SeqList'].seq
            TDT = open('%s.compass.tdt' % self.info['Name'], 'w')  # Output file
            TDT.write('Aln\tAln\tAln\tAln')    # Position, Server, Group, Feature
            for seq in seqlist:
                seqpos[seq] = 0
                for i in range(5):  # Pos, AA, Score, Perc, SA
                    TDT.write('\t%s' % seq.shortName())
            TDT.write('\nPos\tServer\tGroup\tFeature')    # Position, Server, Group, Feature
            for seq in seqlist:
                TDT.write(rje.join(['','Pos','AA','Score','Perc','SA'],'\t'))

            ### <1> ### Work through residues
            _stage = '<1> Residues'
            for r in range(seqlist[0].seqLen()):
                ## <a> ## Sequences with residues
                _stage = '<1a> Sequence Residues'
                res_seq = []    # List of sequences with residue at pos r
                for seq in seqlist:
                    if seq.info['Sequence'][r] != '-':
                        seqpos[seq] += 1
                        res_seq.append(seq)
                #print res_seq
                ## <b> ## ResidueScan
                _stage = '<1b> Res Scan'
                res_scan = {}   # Dic of dics of enzymes from myscan [accnum:[enz_info:stats]]
                for seq in res_seq:
                    accnum = seq.info['AccNum']
                    pos = '%d' % seqpos[seq]
                    if myscan.has_key(accnum) and myscan[accnum].has_key(pos):   # Scansite entry for this position
                        res_scan[accnum] = myscan[accnum].pop(pos)  # Do not need position (only one dealt with at a time)
                #self.verbose(0,0,res_scan,1)
                ## <c> ## Enzyme Info Keys
                _stage = '<1c> List of Enzyme Info keys'
                res_enz = []
                for enzdic in res_scan.values():
                    for enz_info in enzdic.keys():
                        if enz_info not in res_enz:
                            res_enz.append(enz_info)
                if res_enz:
                    self.verbose(1,3,'%d features for AlnPos %d: %s' % (len(res_enz),r,res_enz),1)
                ## <d> ## Output from list
                _stage = '<1d> Output'
                for enz_info in res_enz:
                    TDT.write('\n%d\t%s\t%s\t%s' % (r+1, enz_info[0], enz_info[1], enz_info[2]))
                    for seq in seqlist:
                        accnum = seq.info['AccNum']
                        if seq in res_seq:  # Has residue
                            if res_scan.has_key(accnum) and res_scan[accnum].has_key(enz_info):  # Has this enzyme!
                                rsa = res_scan[accnum].pop(enz_info)
                                if rsa['aa'] == 'NULL' or rsa['aa'] == 'X':
                                    rsa['aa'] = seq.info['Sequence'][r]
                                if rsa['aa'] != seq.info['Sequence'][r]:
                                    self.log.printLog('#ERR','%s Alignment/Prediction Server Error! %s residue %d (aln %d) = %s. Prediction Result = %s.' % (enz_info,accnum,seqpos[seq],r+1,seq.info['Sequence'][r],rsa['aa']))

                                    #if res_scan[accnum]['aa'] == seq.info['Sequence'][r+1]:     # Shift one right - seems to be common Scansite Error!
                                    #    self.log.printLog('#ERR','=> shifted prediction one position right.')
                                    #    #print 'X\t\t%s\t\t\t' % seq.info['Sequence'][r]
                                    #    TDT.write('X\t\t%s\t\t\t' % seq.info['Sequence'][r])
                                    #    newpos = '%d' % (seqpos[seq] + 1)
                                    #    if newpos not in self.scansite[accnum].keys():
                                    #        myscan[accnum][newpos] = {}
                                    #        myscan[accnum][newpos]['aa'] = seq.info['Sequence'][r+1]
                                    #        myscan[accnum][newpos]['sa'] = 'NULL'
                                    #    myscan[accnum][newpos][enzyme] = copy.deepcopy(rsa)
                                    #    #print accnum, newpos, enzyme, self.scansite[accnum][newpos]
                                    #else:
                                    TDT.write(rje.join(['','%d' % seqpos[seq],'*%s*' % rsa['aa'],rsa['score'],rsa['perc'],rsa['sa']],'\t'))
                                #print ['','%d' % seqpos[seq],res_scan[accnum]['aa'],rsa['score'],rsa['perc'],res_scan[accnum]['sa']]
                                else:
                                    TDT.write(rje.join(['','%d' % seqpos[seq],rsa['aa'],rsa['score'],rsa['perc'],rsa['sa']],'\t'))
                            elif self.dict['Feature'].has_key(accnum):   # Residue but no enzyme
                                TDT.write('\t\t%s\t\t\t' % seq.info['Sequence'][r])
                            else:   # No scan
                                TDT.write(rje.join(['','NULL',seq.info['Sequence'][r],'NULL','NULL','NULL'],'\t'))
                        else:   # No residue
                            TDT.write('\t-' * 5)
                ## <e> ## Checks
                _stage = '<1e> Checks'
                scancheck = {}
                for accnum in res_scan.keys():
                    if res_scan[accnum]:
                        scancheck[accnum] = res_scan[accnum]
                if scancheck:
                    self.log.errorLog('res_scan still has entries remaining: %s' % scancheck,printerror=False,quitchoice=True)
                            
            ### <2> ### End Check
            _stage = '<2> End Check'
            TDT.close()
            scancheck = {}
            for accnum in myscan.keys():
                if myscan[accnum]:
                    scancheck[accnum] = myscan[accnum]
            if scancheck:
                self.log.errorLog('myscan still has entries remaining: %s' % scancheck,printerror=False,quitchoice=True)
        except:
            self.log.errorLog('Major Problem with compResults(%s): ' % _stage,False,True)
#########################################################################################################################
## End of SECTION II: Compass Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:        
        compass = Compass(log=mainlog,cmd_list=cmd_list)
        compass.run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
