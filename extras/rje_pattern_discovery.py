#!/usr/local/bin/python

# See below for name and description
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
Module: 	  rje_pattern_discovery
Description:  Pattern Discovery Module 
Version:	  1.3
Last Edit:	  07/10/06
Copyright (C) 2006	Richard J. Edwards - See source code for GNU License Notice

Function:
	Calls and tidies TEIRESIAS. Will add SLIMDISC in time. Will also read motifs and output information content.

SLiMDisc Searching (Tested with slimdisc_V1.3.py):
	seqin=FILE			: Sequence File to search [None]
	slimfiles=LIST		: List of files for SLiMDisc discovery. May have wildcards. (Over-ruled by seqin=FILE.) [*.fas] 
	minsup=X    		: Min. number of sequences to have in file [3]
	maxsup=X    		: Max. number of sequences to have in file [0]
	slimdisc=T/F		: Whether to run list of files through SLiMDisc [False]
	slimopt="X"			: Text string of additional SLiMDisc options [""]
	useres=T/F			: Whether to use existing results files or overwrite (-BT -TT) [True]
	remhub=X			: If X > 0.0, removes "hub" protien ("HUB_PPI") and any proteins >=X% identity to hub [0.0]
						: Renames datasets as -RemHub, -KeptHub or -NoHub
	slimsupport=X		: Min. SLiMDisc support (-S X). If < 1, it is a proportion of input dataset size. [0.1]
	slimranks=X			: Return top X SLiMDisc ranked sequences [1000]
	slimwall=X			: TEIRESIAS walltime (minutes) in SLiMDisc run (-W X) [60]
	slimquery=T/F		: Whether to pull out Query Protein from name of file qry_QUERY [False]
	memsaver=T/F		: Whether to run SLiMDisc in memory_saver mode (-X T) [True]
	bigfirst=T/F		: Whether to run the biggest datasets first (e.g. for ICHEC taskfarm) [True]
	slimversion=X		: Version of SLiMDisc to run (See slimcall=X for batch file jobs) [1.4]  

Batch File Output options:
	batchout=FILE		: Create a batch file containing individual seqin=FILE calls. [None]
	slimcall="X"		: Call for SLiMDisc in batch mode. May have leading commands. ["python slimdisc_V1.4.py"]

TEIRESIAS Searching: *Currently windows-tested only* (Pretty obselete with functional SLiMDisc)
	mysql=T/F			: MySQL mode [True]
	delimit=X			: Text delimiter [\\t]
	info=FILE			: Calculate information content of motifs in FILE (based on AA Freq from seqin, if given) [None]
	igap=X				: Information Content "Gap penalty" (wildcard penalisation) [0]
	outfile=FILE		: Output file name [seqin.teiresias.*]
	teiresias=T/F		: Whether to perform TEIRESIAS search on seqin=FILE [True]
	teiresiaspath=PATH	: Path to TEIRESIAS ['c:/bioware/Teiresias/teiresias_char.exe'] * Use forward slashes (/)
	teiresiasopt=X		: Options for TEIRESIAS Search (Remember "X Y Z" for spaced cmds) e.g. "-bc:\bioware\Teiresias\equiv.txt"
						["-l3 -w10 -c1 -k2 -p"]

Uses general modules: copy, glob, math, os, re, string, sys, time
Uses RJE modules: rje, rje_seq, presto
Other modules needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot, rje_motif
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS																			#
#########################################################################################################################
import copy
import glob
import math
import os
import re
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_seq
import presto
import rje_motif
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 1.0 - Added information content calculation from Norman Davey
# 1.1 - Added calling SLiMDisc
# 1.2 - Updated for use with slimdisc_V1.3 and added batchout mode.
# 1.3 - Added slimversion=X		: Version of SLiMDisc to run [1.4]
#########################################################################################################################
### Major Functionality to Add
# [Y] : Add additional SLIMDISC options
# [?] : Add compilation of SLiMDisc results from *.rank files (see slim_pickings.py)
# [ ] : Add use of orthologue conservation to assessment of SLiMDisc results
# [Y] : Add batchout=FILE to pipe slimdisc calling commands etc. into a file rather than running. (Still renames hubs)
# [ ] : Add rerun=T/F			: Re-run SLiMDisc even if *.rank results are newer than *.out resulst [True]
#########################################################################################################################
def makeInfo(): 	### Makes Info object
	'''Makes rje.Info object for program.'''
	try:
		start_time = time.time()
		program = 'RJE_PATTERN_DISCOVERY'
		version = '1.4'
		last_edit = 'October 06'  
		description = 'Pattern Discovery Module'
		author = 'Dr Richard J. Edwards.'
		info = rje.Info(program,version,last_edit,description,author,start_time)
		return info
	except:
		print 'Problem making Info object.'
		raise
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Show RJE_SEQ commandline options?'):
                out.verbose(-1,4,text=rje_seq.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    ### <0> ### Objects setup
    try:
        ## <a> ## Initial Command Setup & Info
        cmd_list = sys.argv[1:]
        info = makeInfo()
        cmd_list = rje.getCmdList(cmd_list,info=info)      ### Load defaults from program.ini
        ## <b> ## Out object
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#############################################################################################################################
### END OF SECTION I																									#
#########################################################################################################################

													### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Class PatternDiscovery: 																				#
#########################################################################################################################
class PatternDiscovery(rje.RJE_Object): 	
	'''
	PatternDiscovery Class. Author: Rich Edwards (2005).

	Info:str
	- Name = Output file name
	- SeqIn = Name of single input file
	- TeiresiasPath = Path to TEIRESIAS ['c:/bioware/Teiresias/teiresias_char.exe'] * Use forward slashes (/)
	- TeiresiasOpt = Options for TEIRESIAS Search (Remember "X Y Z" for spaced cmds) ["-l3 -w10 -c1 -k2 -p"]
	- Info = List of motifs for information content calculations
	- SlimOpt = Text string of additional SLiMDisc options [""]
	- BatchOut = Create a batch file containing individual seqin=FILE calls. [None]
	- SlimCall = Call for SLiMDisc in batch mode. May have leading commands. ["python slimdisc_V1.4.py"]
	- SlimVersion = Version of SLiMDisc to run [1.4]
	
	Opt:boolean
	- Teiresias = Whether to perform TEIRESIAS search on seqin=FILE [False]
	- SLiMDisc = Whether to run list of files through SLiMDisc [False]
	- UseRes = Whether to use existing results files or overwrite (-BT -TT) [True]
	- SlimQuery = Whether to pull out Query Protein from name of file qry_QUERY [False]
	- BigFirst = Whether to run the biggest datasets first (e.g. for ICHEC taskfarm) [True]

	Stat:numeric
	- InfoGapPen = Information Content "Gap penalty" (wildcard penalisation) [0]
	- MinSup = Min. number of sequences to have in file [3]
	- MaxSup = Max. number of sequences to have in file [200]
	- RemHub = If X > 0.0, removes "hub" protien ("HUB_PPI") and any proteins >=X% identity to hub [0.0]
	- SlimRanks = Return top X SLiMDisc ranked sequences [1000]
	- SlimWall = TEIRESIAS walltime (minutes) in SLiMDisc run (-W X) [60]
	- SlimSupport = Min. SLiMDisc support (-S X). If < 1, it is a proportion of input dataset size [0.1

    List:list
    - Pattern = List of pattern objects
    - SlimFiles = List of files for SLiMDisc discovery. May have wildcards [*.fas]

    Dict:dictionary    

	Obj:RJE_Objects
	- SeqList = Sequence List Object
	'''
#########################################################################################################################
	### <1> ### Class Initiation etc.: sets attributes																	#
#########################################################################################################################
	def _setAttributes(self):	### Sets Attributes of Object
		'''
		Sets Attributes of Object:
		- Info:str ['TeiresiasPath','TeiresiasOpt','Info','SlimOpt','SeqIn','BatchOut','SlimCall','SlimVersion']
		- Opt:boolean ['Teiresias','MySQL','SLiMDisc','UseRes','SlimQuery','BigFirst']
		- Stats:float ['InfoGapPen','MinSup','MaxSup','RemHub','SlimRanks','SlimWall','SlimSupport']
		- List:list ['Pattern','SlimFiles']
		- Dict:dictionary []
		- Obj:RJE_Object []
		'''
		### <a> ### Basics 
		self.infolist = ['TeiresiasPath','TeiresiasOpt','Info','SlimOpt','SeqIn','BatchOut','SlimCall','SlimVersion']
		self.optlist = ['Teiresias','MySQL','SLiMDisc','UseRes','SlimQuery','BigFirst']
		self.statlist = ['InfoGapPen','MinSup','MaxSup','RemHub','SlimRanks','SlimWall','SlimSupport']
		self.listlist = ['Pattern','SlimFiles']
		self.dictlist = []
		self.objlist = []
		### <b> ### Defaults
		self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
		self.info['TeiresiasPath'] = 'c:/bioware/Teiresias/teiresias_char.exe'
		self.info['TeiresiasOpt'] = '-l3 -w10 -c1 -k2 -p'
		self.info['SlimOpt'] = ''
		self.opt['UseRes'] = True
		self.opt['MemSaver'] = True
		self.stat['SlimRanks'] = 1000
		self.stat['SlimWall'] = 60
		self.stat['MinSup'] = 3
		self.stat['MaxSup'] = 0
		self.stat['SlimSupport'] = 0.1
		self.list['SlimFiles'] = ['*.fas']
		self.info['SlimCall'] = 'python slimdisc_V1.4.py'
		self.info['SlimVersion'] = '1.4'
		self.opt['BigFirst'] = True
#########################################################################################################################
	def _cmdList(self): 	### Sets Attributes from commandline
		'''
		Sets attributes according to commandline parameters:
		- see .__doc__ or run with 'help' option
		'''
		for cmd in self.cmd_list:
			try:
				### <a> ### General Options
				self._generalCmd(cmd)
				self._forkCmd(cmd)	# Delete if no forking
				### <b> ### Class Options
				self._cmdRead(cmd,type='info',att='Name',arg='outfile')
				self._cmdRead(cmd,type='info',att='BatchOut')
				self._cmdRead(cmd,type='opt',att='Teiresias')
				self._cmdRead(cmd,type='opt',att='SLiMDisc')
				self._cmdRead(cmd,type='opt',att='UseRes')
				self._cmdRead(cmd,type='int',att='MinSup')
				self._cmdRead(cmd,type='int',att='MaxSup')
				self._cmdRead(cmd,type='info',att='TeiresiasPath')
				self._cmdRead(cmd,type='info',att='TeiresiasOpt')
				self._cmdRead(cmd,type='info',att='SlimOpt')
				self._cmdRead(cmd,type='info',att='Info')
				self._cmdRead(cmd,type='stat',att='InfoGapPen',arg='igap')
				self._cmdRead(cmd,type='glist',att='SlimFiles')
				self._cmdRead(cmd,type='opt',att='SlimQuery')
				self._cmdRead(cmd,type='stat',att='RemHub')
				self._cmdRead(cmd,type='int',att='SlimRanks')
				self._cmdRead(cmd,type='stat',att='SlimWall')
				self._cmdRead(cmd,type='stat',att='SlimSupport')
				self._cmdRead(cmd,type='info',att='SeqIn')
				self._cmdRead(cmd,type='info',att='SlimCall')
				self._cmdRead(cmd,type='info',att='SlimVersion')
				self._cmdRead(cmd,type='opt',att='BigFirst')
			except:
				self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
	### <2> ### Main Class Methods																						#
#########################################################################################################################
	def run(self):		### Main Run method
		'''
		Main Run method.
		'''
		try:
			### SLiMDisc Run ###
			if self.opt['SLiMDisc']:
				return self.slimDisc()
			
			### TEIRESIAS ###
			if self.opt['Teiresias']:
				## Setup ##
				seqlist = rje_seq.SeqList(self.log,self.cmd_list)
				infile = '%s.teiresias.fas' % rje.baseFile(seqlist.info['Name'],True)
				outfile = '%s.teiresias.out' % rje.baseFile(seqlist.info['Name'],True)
				run_teiresias = True
				if rje.isYounger(outfile,infile) == outfile:
					if self.stat['Interactive'] < 1 or not rje.yesNo('%s and %s exist already. Regenerate?' % (infile,outfile),'N'):
						run_teiresias = False
				## Run TEIRESIAS ##
				if run_teiresias:
					seqlist.saveFasta(seqfile=infile,name='Teiresias')	### Saves sequences in fasta format
					command = rje.makePath(self.info['TeiresiasPath'],True)
					command += ' -i%s -o%s %s' % (infile,outfile,self.info['TeiresiasOpt'])
					self.log.printLog('#CMD',command)
					os.system(command)
				## Read Results ##
				self.verbose(0,2,'Reading TEIRESIAS output from %s...' % outfile,1)
				self.list['Pattern'] = []
				RESULTS = open(outfile,'r')
				line = RESULTS.readline()
				while line:
					if rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)\s+(\d.+\d)$',line): # New pattern
						self.addTeiresiasPattern(rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)\s+(\d.+\d)$',line))
					elif len(line) > 3 and line[0] != '#':
						self.log.errorLog('Did not recognise line: %s' % line,False,False)
					line = RESULTS.readline()
				RESULTS.close()
				patx = len(self.list['Pattern'])
				self.log.printLog('#PAT','%s TEIRESIAS patterns read from %s.' % (rje.integerString(patx),outfile))
				## Calculate Information Content ##
				aafreq = seqlist.aaFreq()
				self.verbose(0,3,'Calculating Information Content & Length stats...',0)
				occx = 0
				for pattern in self.list['Pattern']:
					pattern.stat['Info'] = self.calculateScore(pattern.info['Pattern'],aafreq)
					pattern._makeLength()
					occx += 1
					rje.progressPrint(self,occx,patx/100,patx/10)
				self.verbose(0,1,'...Done!',2)
				## Prepare Results ##
				delimit = rje.getDelimit(self.cmd_list)
				if self.info['Name'] == 'None':
					self.info['Name'] = '%s.teiresias.%s' % (rje.baseFile(seqlist.info['Name'],True),rje.delimitExt(delimit))
				if self.opt['MySQL']:	# Two tables
					patfile = os.path.splitext(self.info['Name'])
					occfile = '%s.occ%s' % (patfile[0],patfile[1])
					patfile = '%s.patterns%s' % (patfile[0],patfile[1])
					if self.opt['Append']:
						PATFILE = open(patfile,'a')
						OCCFILE = open(occfile,'a')
					else:
						PATFILE = open(patfile,'w')
						rje.writeDelimit(PATFILE,['pattern','tot_occ','seq_occ','info','len','fix','wild'],delimit)
						OCCFILE = open(occfile,'a')
						rje.writeDelimit(OCCFILE,['seq_id','pos','pattern','pat_match'],delimit)
				else:
					if self.opt['Append']:
						RESFILE = open(self.info['Name'],'a')
					else:
						RESFILE = open(patfile,'w')
						rje.writeDelimit(RESFILE,['Sequence Name','Position','Pattern','Match','Total Occurrences','Num Sequences','Information Content','Length','Fixed','Wildcard'],delimit)
				## Save Results ##
				occx = 0
				for pattern in self.list['Pattern']:
					patstats = []
					for stat in ['OccCount','SeqCount','Info','Length','Fixed','Wildcards']:
						patstats.append('%d' % pattern.stat[stat])
					patstats[2] = '%.3f' % pattern.stat['Info']
					if self.opt['MySQL']:	# Two tables
						rje.writeDelimit(PATFILE,[pattern.info['Pattern']] + patstats,delimit)
					for occ in rje.sortKeys(pattern.occ):
						seq = seqlist.seq[occ]
						for pos in pattern.occ[occ]:
							match = seq.info['Sequence'][pos:(pos+pattern.stat['Length'])]
							outlist = [seq.shortName(),'%d' % pos,pattern.info['Pattern'],match]
							if self.opt['MySQL']:	# Two tables
								rje.writeDelimit(OCCFILE,outlist,delimit)
							else:
								rje.writeDelimit(RESFILE,outlist+patstats,delimit)
							occx += 1
				if self.opt['MySQL']:	# Two tables
					PATFILE.close()
					OCCFILE.close()
					self.log.printLog('#OUT','%s patterns output to %s.' % (rje.integerString(patx),patfile))
					self.log.printLog('#OUT','%s pattern occurrences output to %s.' % (rje.integerString(occx),occfile))
				else:
					RESFILE.close()
					self.log.printLog('#OUT','%s occurrences of %s patterns output to %s.' %
									  (rje.integerString(occx),rje.integerString(patx),self.info['Name']))

			### InfoContent ###
			elif self.info['Info'] != 'None':
				## Setup ##
				alphabet = rje_seq.alph_protx 
				if not os.path.exists(self.info['Info']):
					self.log.errorLog('Input file %s missing!' % self.info['Info'],False,False)
					return False
				else:
					mypresto = presto.Presto(self.log,self.cmd_list)
					mypresto.loadMotifs(file=self.info['Info'],clear=True)
				seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T'])
				if seqlist.seqNum() > 0:
					aafreq = seqlist.aaFreq(alphabet=None,fromfile=None,loadfile=None,total=False)  ### Returns dictionary of AA (& gap etc.) frequencies
				else:
					aafreq = {}
					for aa in alphabet:
						aafreq[aa] = 1.0 / len(alphabet)
				alphabet = aafreq.keys()
				maxinfo = 0 
				for aa in alphabet:
					maxinfo +=  (aafreq[aa] * math.log(aafreq[aa],2))
				## Output ##
				delimit = rje.getDelimit(self.cmd_list)
				ext = rje.delimitExt(delimit)
				outfile = '%s.info.%s' % (rje.baseFile(self.info['Info'],True,['.txt','.%s' % ext]),ext)
				if self.opt['Append']:
					OUTFILE = open(outfile,'a')
				else:
					OUTFILE = open(outfile,'w')
					rje.writeDelimit(OUTFILE,['motif','pattern','info'],delimit)
				
				## Calculate Information Scores ##
				for motif in mypresto.motif:
					self.verbose(2,4,motif.info['Sequence'],0)
					pattern = string.replace(motif.info['Sequence'],'X','.')
					elements = string.split(pattern,'-')
					pattern = ''
					for el in elements:
						if el.find('.{') == 0:	# Ambiguous spacer length - compress
							pattern += '.'
						else:
							pattern += el
					self.verbose(2,2,'=> %s' % pattern,1)
					motif.stat['Info'] = self.calculateInformationContent(pattern,aafreq,maxinfo,self.stat['InfoGapPen'])
					self.verbose(0,3,'%s (%s) = %.2f' % (motif.info['Name'],pattern,motif.stat['Info']),1)
					## Output ##
					rje.writeDelimit(OUTFILE,[motif.info['Name'],pattern,'%.2f' % motif.stat['Info']],delimit)
				
				## Finish ##
				OUTFILE.close()
		except:
			self.log.errorLog('Error in run().',printerror=True,quitchoice=False)
			raise	# Delete this if method error not terrible
#########################################################################################################################
	def addTeiresiasPattern(self,patstats): ### Adds a pattern using given patstats
		'''
		Adds a pattern using given patstats.
		>> patsats:list of strings [tot_occ,seq_occ,pattern,occ_list]
		'''
		### Setup ###
		[tot_occ,seq_occ,pattern,occ_list] = patstats
		occ_list = string.split(occ_list)
		if len(occ_list) < 1 or (len(occ_list)/2 != len(occ_list)/2.0):
			self.log.errorLog('Problem with pattern %s: wrong length for occ_list.' % pattern,False,False)
			return
		newpat = Pattern(self.log,self.cmd_list)
		### Details ###
		newpat.info['Pattern'] = pattern
		newpat.stat['OccCount'] = string.atoi(tot_occ)
		newpat.stat['SeqCount'] = string.atoi(seq_occ)
		#newpat._makeLength()
		### Occurrences ###
		newpat.occ = {}
		while occ_list:
			seq = string.atoi(occ_list.pop(0))
			pos = string.atoi(occ_list.pop(0))
			if seq in newpat.occ.keys():
				newpat.occ[seq].append(pos)
			else:
				newpat.occ[seq] = [pos]
		self.list['Pattern'].append(newpat)
#########################################################################################################################
	### <3> ### Norman AERPIMP Methods																					#
#########################################################################################################################
	def calculateScore(self,pattern,aaOcc): ### Calculates Information Score - Taken from AERPIMP
		'''
		Calculates Information Score - Taken from AERPIMP. Autor: Norman Davey.
		>> pattern:str = pattern
		>> aaOcc:dictionary of {aa:freq}
		<< information score
		'''
		try:	   
			count = 0	
			p1 = re.compile('\[[^\]]*\]')
			p2 = re.compile('(?<=])[^\[]*')

			data = p1.findall(pattern)
			patternTemp = "]" + pattern
			informationContent = 0
			for values in data:
				count += 1
				sum = 0.0
				for AA in values[1:-1]:
					if sum == 0.0:	
						sum = aaOcc[AA]
					else:
						sum += aaOcc[AA]

				informationContent +=  ((sum)*math.log(sum,2) - 4.1) /2  #/len(values[1:-1])

			count += pattern.count(".")
			data = p2.findall(patternTemp)
			
			stringTmp = ""
			for values in data:
				stringTmp += values
			stringTmp = stringTmp.replace(".","")
			#informationContent += 0.05*pattern.count(".")

			for AA in stringTmp:	
				try:
					informationContent +=  aaOcc[AA]*math.log((aaOcc[AA]),2) - 4.1 
				except:
					print pattern
			#if pattern == "D.[ILMV]Y..L..R....Y..L":
			#	sum = 0
			#	for values in aaOcc:
			#		try:
			#			print values,aaOcc[values],aaOcc[values]*math.log(aaOcc[values],2)
			#			sum += aaOcc[values]*math.log(aaOcc[values],2)
			##		except:
			#			pass
			#			print
			#	print 
			#	print informationContent
			#	print sum
			#print math.log(20,2),informationContent
			return	-informationContent - 1*pattern.count(".")
		except:
			self.log.errorLog('Problem with Norman\'s calculateScore Method.')
#########################################################################################################################
	def calculateInformationContent(self,pattern='',aafreq={},maxinfo=0.0,gap_pen=0.0):	### Calculates Information Content
		'''
		Calculates Information Content of a pattern given aa frequencies. Adapted from Norman Davey.
		>> pattern:str = motif pattern
		>> aafreq:dict = aa frequencies
		>> maxinfo:float = maximum information content score
		'''
		#!# NB. This needs neatening up and testing, I think. #!#
		try:
			### Setup ###
			aaOcc = aafreq	# ND variable
			if not maxinfo:
				maxinfo = rje_motif.maxInfo(aafreq)
			#gapPenalty = 0.5	# Not sure what the gap penalty is for. Something to do with wildcards
			gapPenalty = gap_pen
			count = 0
			p1 = re.compile('\[[^\]]*\]')
			p2 = re.compile('(?<=])[^\[]*')

			data = p1.findall(pattern)
			patternTemp = "]" + pattern
			informationContent = 0

			for values in data:
				count += 1
				sum = 0.0
				denominator = 0
				for AA in values[1:-1]:
					denominator += aaOcc[AA]

				for AA in values[1:-1]:
					if sum == 0.0:
						sum = -aaOcc[AA]/denominator*math.log((aaOcc[AA]/denominator),2)
					else:
						sum += -aaOcc[AA]/denominator*math.log((aaOcc[AA]/denominator),2)

				informationContent +=  (-sum - maxinfo)

			count += pattern.count(".")
			data = p2.findall(patternTemp)

			stringTmp = string.join(data,'')
			stringTmp = stringTmp.replace(".","")

			for AA in stringTmp:
				try:
					informationContent +=  -aaOcc[AA]/aaOcc[AA]*math.log((aaOcc[AA]/aaOcc[AA]),2) - maxinfo
				except:
					print pattern

			return  informationContent - gapPenalty*pattern.count(".")

		except:
			self.log.errorLog('Disaster during calculateInformationContent(%s).' % pattern)
			raise
#########################################################################################################################
	### <4> ### Run SLiMDisc on Files																					#
#########################################################################################################################
	def slimDisc(self):	### Runs SLiMDisc on batch of files
		'''Runs SLiMDisc on batch of files.'''
		try:
			### Setup ###
			if self.stat['MinSup'] > self.stat['SlimSupport'] and self.stat['SlimSupport'] > 1:
				self.stat['MinSup'] = self.stat['SlimSupport']
			if self.stat['MaxSup'] > 0  and self.stat['MaxSup'] < self.stat['SlimSupport'] and self.stat['SlimSupport'] > 1:
				self.stat['MaxSup'] = self.stat['SlimSupport']
			### Make File List ##
			_stage = 'Make File List'
			if self.info['SeqIn'].lower() not in ['','none']:
				if os.path.exists(self.info['SeqIn']):
					gfiles = [self.info['SeqIn']]
				else:
					self.log.errorLog('"seqin" file "%s" not found! No SLiMDisc analysis.' % self.info['SeqIn'],printerror=False)
					return False
			else:
				gfiles = rje.getFileList(callobj=self,filelist=self.list['SlimFiles'],subfolders=False,summary=False)
			self.log.printLog('#FILES','%s files identified for SLiMDisc analysis.' % rje.integerString(len(gfiles)))
			## Sort by size and filter by MinSup and MaxSup ###
			datasize = {}   # Dictionary for crude sorting of files by total AA content
			seqnum = {}		# Number of sequences in each file
			qry = {}		# Query sequence name (if any) for file
			tmpseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autofilter=F'])
			gx = 0
			while gx < len(gfiles):
				seqfilename = gfiles[gx]
				gx += 1
				seqfile = seqfilename[0:]
				tmpseq.seq = []
				tmpseq.loadSeqs(seqfile)
				## *** Special RemHub process *** ##
				checkhub = True
				for hubtype in ['rem','kept','no']:
					if seqfile.find('-%shub.fas' % hubtype) > 0:
						checkhub = False
				if self.stat['RemHub'] > 0.0 and checkhub:
					if rje.matchExp('(\S+)_PPI',seqfile):
						hub_acc = rje.matchExp('(\S+)_PPI',rje.baseFile(seqfile,strip_path=True))[0]
					else:
						hub_acc = rje.baseFile(seqfile,strip_path=True)
					hub_base = rje.matchExp('(\S+)%s' % hub_acc,seqfilename)[0]
					basefile = seqfile
					while rje.baseFile(basefile) != basefile:
						basefile = rje.baseFile(basefile)
					if tmpseq.querySeq(query=hub_acc):     ### Sets Hub as Query Sequence
						self.log.printLog('#HUB','Removing hub protein %s and >=%.1f%% ID from PPI dataset %s.' % (hub_acc,self.stat['RemHub'],seqfile))
						tmpseq.makeNR(text='Hub protein homologues',nrid=self.stat['RemHub'],blast=tmpseq.seqNum(),nrsim=0,nr_qry=tmpseq.obj['QuerySeq'])
						tmpseq.removeSeq(text='PPI Hub Protein (self-interactor)',seq=tmpseq.obj['QuerySeq'])
						tmpseq.obj['QuerySeq'] = None
						seqfile = '%s-remhub.fas' % basefile
						tmpseq.saveFasta(seqfile=seqfile)	### Saves sequences in fasta format
						keptfile = '%s-kepthub.fas' % basefile
						os.rename(seqfilename,keptfile)
						gfiles.append(keptfile)
					else:
						seqfile = '%s-nohub.fas' % basefile
						os.rename(seqfilename,seqfile)
						self.log.printLog('#HUB','Hub protein %s not in PPI dataset %s => %s.' % (hub_acc,seqfilename,seqfile))
						#X#print tmpseq.obj['QuerySeq']
				## Support Range ###				
				if tmpseq.seqNum() < self.stat['MinSup'] or (self.stat['MaxSup'] > 0 and tmpseq.seqNum() > self.stat['MaxSup']):
					self.log.printLog('#REJ','%s rejected: %s sequences = outside acceptable range of %d-%d.' % (seqfile,rje.integerString(tmpseq.seqNum()),self.stat['MinSup'],self.stat['MaxSup']))
					continue
				aasize = tmpseq.aaCount()
				self.log.printLog('#AA','%s = %s aa.' % (seqfile,rje.integerString(aasize)))
				while datasize.has_key(aasize):
					aasize += 1
				datasize[aasize] = seqfile
				seqnum[seqfile] = tmpseq.seqNum()
				## Query ##
				qry[seqfile] = None
				if self.opt['SlimQuery']:
					if rje.matchExp('qry_(\S+)\.',seqfilename):
						if tmpseq.querySeq(query=rje.matchExp('qry_(\S+)\.',seqfilename)[0]):     ### Sets Query Sequence if appropriate
							qry[seqfile] = tmpseq.obj['QuerySeq'].shortName()
			self.log.printLog('#INF','%s Datasets to process.' % rje.integerString(len(seqnum)))

			### Batch Output Mode ###
			batchout = None
			if self.info['BatchOut'].lower() not in ['','none']:
				batchout = self.info['BatchOut']
				if not self.opt['Append'] and os.path.exists(batchout):
					rje.backup(self,batchout)

			### Work through Files ###
			_stage = 'Work through files'
			for key in rje.sortKeys(datasize,revsort=self.opt['BigFirst']):
				seqfile = datasize[key]
				basefile = seqfile
				while rje.baseFile(basefile) != basefile:
					basefile = rje.baseFile(basefile)
				base = rje.baseFile(basefile,True)
				self.log.printLog('#DAT',seqfile,timeout=False)
				if not self.opt['UseRes']:
					slim_cmd = '-BT -TT'
				else:
					## Detect old files ##
					_stage = 'Detect old files'
					old_rank = '%s/%s.rank' % (basefile,base)
					self.log.printLog('#RES','Existing SLiMDisc Output?: %s' % (os.path.exists(old_rank)))
					old_b_list = glob.glob('%s/results/*.blastp' % basefile)
					old_t_file = '%s/%s.fasta.out' % (basefile,base)
					self.log.printLog('#RES','Existng TEIRESIAS Output?: %s' % (os.path.exists(old_t_file)))
					self.log.printLog('#RES','%s of %s BLAST files detected.' % (rje.integerString(len(old_b_list)),rje.integerString(seqnum[seqfile])))
					## TEIRESIAS ##
					if (os.path.exists(old_rank) or len(old_b_list) > 0) and os.path.exists(old_t_file):  # BLAST started: TEIRESIAS finished!
						slim_cmd = '-TF'
					else:
						slim_cmd = '-TT'
					## BLAST ##
					if len(old_b_list) != seqnum[seqfile]:	# Need BLAST
						slim_cmd += ' -BT'
					else:
						slim_cmd += ' -BF'
				## Query ##
				if self.opt['SlimQuery'] and qry[seqfile]:
					slim_cmd += ' -q %s' % qry[seqfile]
				## Ranks ##
				slim_cmd += ' -n %d' % self.stat['SlimRanks']
				## Support ##
				if self.stat['SlimSupport'] > 0 and self.stat['SlimSupport'] < 1:
					slim_cmd += ' -S %.1f' % self.stat['SlimSupport']
				elif self.stat['SlimSupport'] > 0:
					slim_cmd += ' -S %d' % self.stat['SlimSupport']
				## WallTime ##
				slim_cmd += ' -W %d' % self.stat['SlimWall']
				## MemSaver ##
				if self.opt['MemSaver']:
					slim_cmd += ' -X T'
				else:
					slim_cmd += ' -X F'
				## SlimOpt ##
				if self.info['SlimOpt']:
					slim_cmd += ' %s' % self.info['SlimOpt']
				## Perform SLiMDisc Run ##
				_stage = 'Peform SLiMDisc Run (%s)' % (seqfile)
				if batchout:
					BATCH = open(batchout,'a')
					BATCH.write('%s -i %s -Q0 %s\n' % (self.info['SlimCall'],seqfile,slim_cmd))
					BATCH.close()
				else:
					if self.stat['Verbose'] > 0:
						syscmd = 'python /home/richard/Python_Modules/slimdisc_V%s.py -i %s -Q2 %s' % (self.info['SlimVersion'],seqfile,slim_cmd)
					else:
						syscmd = 'python /home/richard/Python_Modules/slimdisc_V%s.py -i %s -Q0 %s' % (self.info['SlimVersion'],seqfile,slim_cmd)
					self.log.printLog('#SYS',syscmd)
					os.system(syscmd)
				if not batchout:
					new_rank = '%s/%s.rank' % (basefile,base)
					self.log.printLog('#RES','New rank result %s produced?: %s' % (new_rank,os.path.exists(new_rank)))

		except:
			self.log.errorLog('rje_pattern_discovery banjaxed in slimDisc() %s' % _stage,quitchoice=True)
#########################################################################################################################
## End of SECTION II: Class PatternDiscovery																			#
#########################################################################################################################

													### ~ ### ~ ###

#########################################################################################################################
### SECTION III:  Class Pattern:																						#
#########################################################################################################################
class Pattern(rje.RJE_Object):	   
	'''
	Pattern Class. Author: Rich Edwards (2005).

	Info:str
	- Pattern = Pattern RegExp
	
	Opt:boolean

	Stat:numeric
	- OccCount = Number of occurrences in all sequences
	- SeqCount = Number of sequences it occurs in
	- Info = Information Content
	- Length = Length of pattern
	- Fixed = Number of fixed postions in pattern
	- Wildcards = Number of wildcard positions in pattern

	Obj:RJE_Objects
	'''
	### Attributes
	occ = {}	# Occurrence dictionary {seq_index:[positions]}
#########################################################################################################################
	### <1> ### Class Initiation etc.: sets attributes																	#
#########################################################################################################################
	def _setAttributes(self):	### Sets Attributes of Object
		'''
		Sets Attributes of Object:
		- Info:str ['Pattern']
		- Opt:boolean []
		- Stats:float ['OccCount','SeqCount','Info','Length','Fixed','Wildcards']
		- Obj:RJE_Object []
		'''
		### <a> ### Basics 
		self.infolist = ['Pattern']
		self.optlist = []
		self.statlist = ['OccCount','SeqCount','Info','Length','Fixed','Wildcards']
		self.objlist = []
		### <b> ### Defaults
		self._setDefaults(info='.',opt=False,stat=0,obj=None)
		self.info['TeiresiasPath'] = 'c:/bioware/Teiresias/teiresias_char.exe'
		self.info['TeiresiasOpt'] = '-l3 -w10 -c1 -k2 -p -bc:\\bioware\\Teiresias\\equiv.txt'
		### <c> ### Other Attributes
		self.occ = {}
#########################################################################################################################
	### <2> ### Main Class Methods																						#
#########################################################################################################################
	def _makeLength(self):	### Makes length from pattern
		'''
		Makes length from pattern.
		'''
		add = 1
		plen = 0
		flen = 0
		for s in self.info['Pattern']:
			if s == '[':
				add = 0
			flen += add
			if s == ']':
				add = 1
			plen += add
		self.stat['Length'] = plen
		self.stat['Wildcards'] = string.count(self.info['Pattern'],'.')
		self.stat['Fixed'] = flen - self.stat['Wildcards']
#########################################################################################################################
### END OF SECTION III: Class Pattern																					#
#########################################################################################################################

													### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM																							#
#########################################################################################################################
def runMain():
	try:
		### <0>  ### Basic Setup of Program
		[info,out,mainlog,cmd_list] = setupProgram()
		
		### <1> ### Rest...
		patterns = PatternDiscovery(mainlog,cmd_list)
		patterns.run()

		### <X> ### End
	except SystemExit:
		return	# Fork exit etc.
	except KeyboardInterrupt:
		mainlog.errorLog("User terminated.\n")
	except:
		print "Unexpected error:", sys.exc_info()[0]
	mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#########################################################################################################################
if __name__ == "__main__":		### Call runMain 
	try:
		runMain()
	except:
		print 'Cataclysmic run error:', sys.exc_info()[0]
	sys.exit()
#########################################################################################################################
### END OF SECTION IV																									#
#########################################################################################################################
