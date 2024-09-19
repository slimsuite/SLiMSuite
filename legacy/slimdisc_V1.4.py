#!/usr/bin/env python

"""
Program:      SLiMDisc
Description:  Short, Linear Motif Discovery
Version:      1.4
Last Edit:    28/08/06
Copyright (C) 2006 Norman Davey - See website for GNU License Notice

Function:
	SLiMDisc provides a simple but effective method for identifying novel short linear motifs (SLiMs) within a set of
	biologically (functionally) related sequences, such as protein-protein interaction (PPI) datasets.

	First, the TEIRESIAS pattern discovery algorithm is used to identify common patterns in the dataset. These will
	generally be dominated by patterns arising in sequences with shared homology. The number of occurrences of the motif
	in the dataset is then normalised to reduce the impact of similarities arising through common evolutionary descent.
	This reduction is based on GABLAM treatment of BLAST local alignments.

	Motifs are then ranked using a score derived from the product of the normalised number of occurrences and the
	information content. The method has been shown to significantly out-perform methods that do not discount evolutionary
	relatedness, when applied to a subset of the ELM (Eukaryotic Linear Motif) database. This method provides a simple
	means to maximise the chance of finding motifs among unrelated sequences, without forcing the user to represent large
	protein families by a single arbitrarily chosen member.

Commandline:

	FILE PARAMETERS
	------------------------------------------------------------------------
	>filename
	i=[filename]	
	sets the input FASTA format file to be the file specified.[default=]
	
	>output
	o=[string]		
	output sets the output file prefix to be the string specified.[default=]
	
	>aminoAcidOccuranceFile 
	a=[filename]	
	aminoAcidOccuranceFile sets the output file for the amino acid 
	probabilities to be the file specified.[default=PeptideOccurances.txt]
	
	>batch	
	b=[T/F]	
	specifies whether the run is for a single file or a batch of files. 
	[default=F]                         
	
	>filetype	
	t=[filetype]	
	sets the batch filetype to be the filetype specified eg. 'fas'. 
	[default= ]    
	
	>Mask_path
	m=[filename]	
	sets the masking file to be the file specified. [default=mask.dat]    
	
	>ini_file
	d= [filename]	
	sets the ini file to be the file specified. [default=slimdisc.ini]              
	
	>memory_saver
	X=[T/F]	
	specifies whether the files created by the run are discarded or not.
	[default=F]   
	
	>overwrite
	O=[T/F]
	specifies if files already created should be over written or not
	[default=T]
	
	
	
	OUTPUT PARAMETERS
	------------------------------------------------------------------------
	>logfile
	l=[filename]	
	logfile sets the logfile to the file specified.[default=]
	
	>no_of_outputs
	n=[integer]		
	no_of_outputs sets the number of ranked results to be returned.
	[default=20]
	
	>view_matrix
	v=[T/F]			
	view_matrix prints the all X all sequence similarity matrix to the log 
	file.[default=F]
	
	>verbosity
	Q=[integer]
	sets the amount information written to the screen by Slimdisc
	[default=0]

	FILTERING PARAMETERS
	------------------------------------------------------------------------
	>use_evolution
	E=[T/F]		
	use_evolution specifies whether the evolutionary filtering and 
	evolutionary scoring are to be used.[default=T]
	
	>use_filtering
	A=[T/F]		
	use_filtering specifies whether the non-evolutionary filters are to be
	used. Non-evolutionary filters are complexity, surface accessibility and
	information content.[default=T]
	
	>information_cutOff
	y=[float]	
	information_cutOff sets the cutoff for information content. Any patterns
	with IC lower than this value will be filtered.[default=0.0]
	
	>complexity_cutOff
	z=[float]	
	complexity_cutOff sets the cutoff for complexity. Any patterns with
	complexity lower than this value will be filtered.[default=0.3]
	
	>surface_cutOff
	z=[float]	
	surface_cutOff sets the cutoff for surface accessibility. Any patterns 
	with surface accessibility lower than this value will be filtered. 
	[default=0.3]
	
	>self_hit	
	H=[T/F]	
	Specifies whether self hits are counted or not. [default=F]     
	
	>query_protein
	q=[protein_name]	
	Specifies a protein which motifs returned must be contained in[default=]                         
	
	>strict_inclusive_masking 
	M=[T/F]	
	When true, proteins not containing inclusive masking region are removed 
	completely else the whole protein is kept. [default=F]                   


	SCORING PARAMETERS
	------------------------------------------------------------------------
	>normalisation
	f=[integer]	
	flavour sets the type of normalisation to be.[default=1]
		Options:
			1=MST - Similarity based normalisation
			2=UHS - Unique domains based normalisation
			3=UP - Non-Overlapping protein clusters normalisation
	
	>normalisation_value
	C=[integer]	
	Specifies the normalisation value for stretching the MST branches
	[default=1]                   

	BLAST PARAMETERS
	------------------------------------------------------------------------
	>eVal 
	e=[float]	
	eVal sets the evalue to be used in the BLAST search.[default=1.0]
	
	>run_formatDB
	D=[T/F]		
	run_formatDB specifies if the input file is to be formatted to be used 
	as a BLAST database.[default=F]
	
	>run_BLAST
	B=[T/F]		
	run_BLAST specifies if the input file is to be BLASTed protein by 
	protein against the inputfiles formatted BLAST database.[default=F]
	
	>BLAST_results
	r=[string]	
	BLAST_results sets the output file prefix to be the string specified.
	[default=]

	TEIRESIAS PARAMETERS
	------------------------------------------------------------------------
	>run_TEIRESIAS
	T=[T/F]		
	run_TEIRESIAS specifies if TEIRESIAS is to be run on the input file.
	[default=F]
	
	>TEIRESIAS_output
	O=[filename] 
	TEIRESIAS_output is the output TEIRESIAS results file.[default=]
	
	>TEIRESIAS_input
	I=[filename] 
	TEIRESIAS_input is the TEIRESIAS input FASTA format file. [default=]
	
	>TEIRESIAS_patternLength
	P=[integer]	
	TEIRESIAS_patternLength is the maximum extent of an elementary pattern. 
	[default=10]
	
	>TEIRESIAS_fixedPositions
	F=[integer]	
	TEIRESIAS_fixedPositions is the number of nondot characters in the 
	pattern. [default=3]
	
	>TEIRESIAS_supportString
	S=[integer]	
	TEIRESIAS_supportString is the minimum allowed support for a pattern.
	[default=3]
	
	>TEIRESIAS_equiv
	R=[integer]	
	TEIRESIAS_equiv is the name of the file that contains the symbol 
	equivalences. [default=T]

	------------------------------------------------------------------------
"""


################################################################################
__doc__ = """
################################################################################

	     
 01010011 01101100 01101001 01101101 01000100 01101001 01110011 01100011
  ____    _   _               ____    _              
 / ___|  | | (_)  _ __ ___   |  _ \  (_)  ___    ___ 
 \___ \  | | | | | '_ ` _ \  | | | | | | / __|  / __|
  ___) | | | | | | | | | | | | |_| | | | \__ \ | (__ 
 |____/  |_| |_| |_| |_| |_| |____/  |_| |___/  \___|

 
  Author:     Norman Davey <norman.davey@ucd.ie>
  Program:    SLimDisc(Short LInear Motif DISCovery)
		
  Date:       08/02/2006
  Version:

  Description:  Application for finding short protein motifs using evolutionary 
  		information to remove patterns from closely related proteins

 01010011 01101100 01101001 01101101 01000100 01101001 01110011 01100011



Commandline:

	FILE PARAMETERS
	------------------------------------------------------------------------
	>filename
	i=[filename]	
	sets the input FASTA format file to be the file specified.[default=]
	
	>output
	o=[string]		
	output sets the output file prefix to be the string specified.[default=]
	
	>aminoAcidOccuranceFile 
	a=[filename]	
	aminoAcidOccuranceFile sets the output file for the amino acid 
	probabilities to be the file specified.[default=PeptideOccurances.txt]
	
	>batch	
	b=[T/F]	
	specifies whether the run is for a single file or a batch of files. 
	[default=F]                         
	
	>filetype	
	t=[filetype]	
	sets the batch filetype to be the filetype specified eg. 'fas'. 
	[default= ]    
	
	>Mask_path
	m=[filename]	
	sets the masking file to be the file specified. [default=mask.dat]    
	
	>ini_file
	d= [filename]	
	sets the ini file to be the file specified. [default=slimdisc.ini]              
	
	>memory_saver
	X=[T/F]	
	specifies whether the files created by the run are discarded or not.
	[default=F]   
	
	>overwrite
	O=[T/F]
	specifies if files already created should be over written or not
	[default=T]
	
	
	
	OUTPUT PARAMETERS
	------------------------------------------------------------------------
	>logfile
	l=[filename]	
	logfile sets the logfile to the file specified.[default=]
	
	>no_of_outputs
	n=[integer]		
	no_of_outputs sets the number of ranked results to be returned.
	[default=20]
	
	>view_matrix
	v=[T/F]			
	view_matrix prints the all X all sequence similarity matrix to the log 
	file.[default=F]
	
	>verbosity
	Q=[integer]
	sets the amount information written to the screen by Slimdisc
	[default=0]

	FILTERING PARAMETERS
	------------------------------------------------------------------------
	>use_evolution
	E=[T/F]		
	use_evolution specifies whether the evolutionary filtering and 
	evolutionary scoring are to be used.[default=T]
	
	>use_filtering
	A=[T/F]		
	use_filtering specifies whether the non-evolutionary filters are to be
	used. Non-evolutionary filters are complexity, surface accessibility and
	information content.[default=T]
	
	>information_cutOff
	y=[float]	
	information_cutOff sets the cutoff for information content. Any patterns
	with IC lower than this value will be filtered.[default=0.0]
	
	>complexity_cutOff
	z=[float]	
	complexity_cutOff sets the cutoff for complexity. Any patterns with
	complexity lower than this value will be filtered.[default=0.3]
	
	>surface_cutOff
	z=[float]	
	surface_cutOff sets the cutoff for surface accessibility. Any patterns 
	with surface accessibility lower than this value will be filtered. 
	[default=0.3]
	
	>self_hit	
	H=[T/F]	
	Specifies whether self hits are counted or not. [default=F]     
	
	>query_protein
	q=[protein_name]	
	Specifies a protein which motifs returned must be contained in[default=]                         
	
	>strict_inclusive_masking 
	M=[T/F]	
	When true, proteins not containing inclusive masking region are removed 
	completely else the whole protein is kept. [default=F]                   


	SCORING PARAMETERS
	------------------------------------------------------------------------
	>normalisation
	f=[integer]	
	flavour sets the type of normalisation to be.[default=1]
		Options:
			1=MST - Similarity based normalisation
			2=UHS - Unique domains based normalisation
			3=UP - Non-Overlapping protein clusters normalisation
	
	>normalisation_value
	C=[integer]	
	Specifies the normalisation value for stretching the MST branches
	[default=1]                   

	BLAST PARAMETERS
	------------------------------------------------------------------------
	>eVal 
	e=[float]	
	eVal sets the evalue to be used in the BLAST search.[default=1.0]
	
	>run_formatDB
	D=[T/F]		
	run_formatDB specifies if the input file is to be formatted to be used 
	as a BLAST database.[default=F]
	
	>run_BLAST
	B=[T/F]		
	run_BLAST specifies if the input file is to be BLASTed protein by 
	protein against the inputfiles formatted BLAST database.[default=F]
	
	>BLAST_results
	r=[string]	
	BLAST_results sets the output file prefix to be the string specified.
	[default=]

	TEIRESIAS PARAMETERS
	------------------------------------------------------------------------
	>run_TEIRESIAS
	T=[T/F]		
	run_TEIRESIAS specifies if TEIRESIAS is to be run on the input file.
	[default=F]
	
	>TEIRESIAS_output
	O=[filename] 
	TEIRESIAS_output is the output TEIRESIAS results file.[default=]
	
	>TEIRESIAS_input
	I=[filename] 
	TEIRESIAS_input is the TEIRESIAS input FASTA format file. [default=]
	
	>TEIRESIAS_patternLength
	P=[integer]	
	TEIRESIAS_patternLength is the maximum extent of an elementary pattern. 
	[default=10]
	
	>TEIRESIAS_fixedPositions
	F=[integer]	
	TEIRESIAS_fixedPositions is the number of nondot characters in the 
	pattern. [default=3]
	
	>TEIRESIAS_supportString
	S=[integer]	
	TEIRESIAS_supportString is the minimum allowed support for a pattern.
	[default=3]
	
	>TEIRESIAS_equiv
	R=[integer]	
	TEIRESIAS_equiv is the name of the file that contains the symbol 
	equivalences. [default=T]

	------------------------------------------------------------------------

################################################################################
"""
################################################################################
import sys,os,re,sets,time,math,traceback,getopt,random,tempfile,string,copy,shutil
from threading import Thread

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),"RJE_modules"))

import rje_uniprot,rje_blast,rje_seq
try:
	from SurfaceAccessibilityPlotter import SurfaceAccessibilityPlotter
	from HomologyPlotter import HomologyPlotter
	from DomainPlotter import DomainPlotter
	from DistributionPlotter import DistributionPlotter
except:
	warning_PIL_not_installed = """
	####################################################################################\n
	The PIL library does not seem to be installed no images or graphs will be generated.
		     To create images associated with the analysis go to 

		           ---------------------------------------
		           http://www.pythonware.com/products/pil/
		           ---------------------------------------
	
			    and download and install PIL library

	#####################################################################################\n
	"""

	if '-Q0' not in sys.argv:
		print warning_PIL_not_installed
		time.sleep(0)
		
	pass


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class fileManipulation:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  fileManipulation

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   readFile(options)
	   @param options

	   createTEIRESIASinputfromProteinList(proteinList,options)
	   @param proteinList
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass 
	#--------------------------------------------------------------------------------------------------------------------#
	def readFileOLD(self,options,filename):
		proteinList = {}
		
		try:
			fileIn = open(filename,"r")	
		except:
			if os.path.isfile(options["filename"]):
				print "Error opening file"
			else:
				errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nInput file does not exist.\nCheck the filepath and try again.\n" + "+"*70
				print errorfile_not_found

			sys.exit()

		sequence = ""
		name = ""
			
		proteinData = fileIn.read().strip()
		pattern1 = re.compile('^>')
		if len(pattern1.findall(proteinData)) == 0:
			print "Input not in a recognisable fasta format.\nIf input in UNIPROT format please use the -LT option"
			sys.exit()
			
		for lines in proteinData.split('\n'):
	
			if lines[0] == '>':
				if len(sequence) > 0:
					proteinList[name] = sequence.strip()
					sequence = ""
				try:
					name = lines.split("|")[1].split()[0]
				except:
					name = lines[1:].split()[0]
			else:
				sequence += lines.replace("\n","")

		proteinList[name] = sequence.strip()

		fileIn.close()
		return proteinList

	def readFile(self,options,filename):
		proteinList = {}
		
		if os.path.isfile(options["filename"]):
			pass
		else:
			errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nInput file does not exist.\nCheck the filepath and try again.\n" + "+"*70
			print errorfile_not_found
			sys.exit()

		seqs = rje_seq.SeqList(log=None,cmd_list=['v=-1','i=-1'])
		seqs.loadSeqs(seqfile=options["filename"],nodup=True)

		for seq in seqs.seq:
			proteinList[seq.info['AccNum']] = seq.info['Sequence']

		return proteinList

	def createTEIRESIASinputfromProteinList(self,proteinList,options):
		outString = ""
		for proteins in proteinList.keys():
			outString += ">" + proteins + " 1\n" + proteinList[proteins] + "\n"

		open(options["TEIRESIAS_input"],"w").write(outString)
	#----------------------------------------------END OF CLASS----------------------------------------------------------#


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class fileManipulationFull:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  fileManipulation

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   readFile(options)
	   @param options

	   createTEIRESIASinputfromProteinList(proteinList,options)
	   @param proteinList
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass 
	#--------------------------------------------------------------------------------------------------------------------#
	def readFile(self,options):
		try:
			fileIn = open(options["filename"],"r")	
		except:
			errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nInput file does not exist.\nCheck the filepath and try again.\n" + "+"*70
			print errorfile_not_found

			sys.exit()		

		proteinData = fileIn.read().strip()
		pattern1 = re.compile('^>')
		
		if len(pattern1.findall(proteinData)) != 0:
			print "Input not in a recognisable UNIPROT format.\nIf input in fasta format please use the -LF option"
			sys.exit()
			
		counter = 0
		mwcounter = 0
		features = {}
		proteinList = {}
		teiresiasInput = []
		offsets = []
		count = 0

		if options['memory_saver'] == 'F':
			try:
				plotter = DomainPlotter()
			except:
				pass
		else:
			pass
			#printLevel(1,options, "Skipping domain plotting")
		
		tempFile = open(options["output"].replace(".rank",".domains.html"),"w")
	
		tempFile.write("<html>")
		tempFile.write("<h2>Domains<br></h2>") 
		
		printLevel(1,options, options["Mask_path"])
		try:
			masksTemp = open(options["Mask_path"],"r").read().strip().split("\n")
			exmasks = []
			inmasks = []
			for mask in masksTemp:
				try:
					if mask.split()[1] == "e":
						exmasks.append(mask.split()[0])
					if mask.split()[1] == "i":
						inmasks.append(mask.split()[0])
				except:
					print mask
		except:
			print "Cannot find mask file"
			sys.exit()

		try:
			exception = open(os.path.join(os.path.dirname(options['Mask_path']),"exceptions.dat"),"r").read().split("\n")
		except:
			print "Cannot find exceptions file"
			open(os.path.join(os.path.dirname(options['Mask_path']),"exceptions.dat"),'w')
			sys.exit()

		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		printLevel(1,options,"Masking TEIRESIAS Input")
		printLevel(1,options,"Inclusive Masking:\t" + str(inmasks))
		printLevel(1,options,"Exclusive Masking:\t" + str(exmasks))

		printLevel(1,options,"exceptions: \t" + str(exception))
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")

		uniprotParser = rje_uniprot.UniProt(cmd_list=['v=-1'])
		uniprotParser.readUniProt(filename=options["filename"])
	
		try:
			maskingString = ''
			for entry in uniprotParser.list['Entry']:
				maskingString += '-'*80 + '\n'
				inclusiveMasking = 0
				entry_name = entry.info['Name']
				sequence = entry.obj['Sequence'].info['Sequence']
				unmasked = [0,len(sequence)]
				maskingString += 'Reading protein ' + entry_name + '\n'

				if entry_name not in proteinList:
					tempFile.write("<h3>" + entry_name + "</h3>\n")

					tempFile.write("<p><img src='" + os.path.abspath("./" + options["input_path"] + "/domains/" + entry_name + ".domain.gif'") + "></p><br>\n")
					proteinList[entry_name] = sequence

					features[entry_name] = {}

					for ft in entry.list['Feature']:

						ftlist = [ft['Start'],ft['End'],ft['Desc']]
						if ft['Type'] in features[entry_name]:
							features[entry_name][ft['Type']].append(list(ftlist))
						else:
							features[entry_name][ft['Type']] = [list(ftlist)]


					splits = []
					temp = sequence
					featureList = []
					
					tempFile.write("<FONT SIZE=2 face='courier'><table width='900' border='1' cellspacing='0' cellpadding='5'>" )
					
					proteinMaskInclusive = copy.deepcopy(list(sequence))
					proteinMaskExclusive = list(sequence)
					
					for x in range(0,len(sequence)):
						proteinMaskInclusive[x] = "x"

					
					maskingString += '\n'
					for inmask in inmasks:
						if inmask in features[entry_name]:
						

							maskingString += 'Inclusive masking ' + inmask + '\n'
							for iter in range(len(features[entry_name][inmask])):
								if features[entry_name][inmask][iter][2].replace('.','') in exception:
									maskingString += 'skipping ' + feature_data[2].replace('.','') + ' (Contained in exceptions file)' + '\n'

									pass
								else:	
									inclusiveMasking = 1
									maskingString +=  features[entry_name][inmask][iter][2].replace('.','') + ' from ' + str(features[entry_name][inmask][iter][0]) + ' to ' + str(features[entry_name][inmask][iter][1]) + '\n'
									for x in range(0,len(sequence)):
										if x >= features[entry_name][inmask][iter][0] -1  and x <= features[entry_name][inmask][iter][1] - 1:
											proteinMaskInclusive[x] = proteinList[entry_name][x]
											
					maskingString += '\n'	
					for featureType in exmasks:
						if featureType in features[entry_name]:
							maskingString += 'Exclusive masking ' + featureType + '\n'
							for feature_data in features[entry_name][featureType]:
								try:
									if feature_data[2].replace('.','') in exception:
										maskingString += 'skipping ' + feature_data[2].replace('.','') + ' (Contained in exceptions file)' + '\n'
	
										pass
									else:
										maskingString +=  feature_data[2].replace('.','') + ' from ' + str(feature_data[0]) + ' to ' + str(feature_data[1]) + '\n'
	
										featureList.append([featureType,feature_data[0],feature_data[1],feature_data[2]])
										tempFile.write("<tr><td align='left'><FONT SIZE=2 face='courier'>" + featureType + "</td><td align='left'><FONT SIZE=2 face='courier'>" + str(feature_data[0]) + "</font></td><td align='left'><FONT SIZE=2 face='courier'><FONT SIZE=2 face='courier'>" + str(feature_data[1]) + "</td><td align='left'><FONT SIZE=2 face='courier'>" + str(feature_data[2]) + "</td></tr>")
										
										
										for i in range(0,len(unmasked) - 1):
											if int(feature_data[0])-1 >= unmasked[i] and int(feature_data[1]) <= unmasked[i + 1]:
												unmasked.insert(i + 1,int(feature_data[0])-1) 
												unmasked.insert(i + 2,int(feature_data[1])) 
												for j in range(int(feature_data[0])-1,int(feature_data[1])):
													proteinMaskExclusive[j] = "x"
								except:
									print feature_data
					
					tempFile.write("</table><br>")
					
					if inclusiveMasking == 0:
						if options['strict_inclusive_masking'] == 'T':
							proteinMask = proteinMaskInclusive
						else:
							proteinMask = proteinMaskExclusive
					else:
						for iter in range(len(proteinMaskInclusive)):
							if proteinMaskExclusive[iter] == 'x':
								 proteinMaskInclusive[iter] = 'x'
							
						proteinMask = proteinMaskInclusive 
						
					maskedString = rje.join(proteinMask,"")
					pattern = re.compile("[^x]+")

					maskingString +=  "\n>" + entry_name + '\n'
					
					for i in range(0,len(maskedString)/60  +1): 
						tempFile.write(maskedString[60*i:60*(i+1)]  + "<br>")
						maskingString +=  maskedString[60*i:60*(i+1)] + '\n'

					iterMatch = pattern.finditer(maskedString)
					for match in iterMatch:
						offsets.append([entry_name,match.span()[0],count])
						teiresiasInput.append([entry_name,match.group()])
				
					if options['memory_saver'] == 'F':
						try:
							plotter.plotDomain(entry_name,featureList,len(sequence),options)
						except:
							maskingString +=  'Skipping plotting domains'
					else:
						pass
						#printLevel(1,options,  "Skipping domain plotting")
						
					count += 1
				else:
					record = iterator.next()
		except:
			errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nFile Format Invalid.\n" + "+"*70
			print errorfile_not_found
			raise
			sys.exit()


		teiresiasInputFile = open(options["TEIRESIAS_input"],"w")

		for i in range(len(teiresiasInput)):
			teiresiasInputFile.write( ">" + teiresiasInput[i][0] + " " + str(i) + "\n") 
			teiresiasInputFile.write( teiresiasInput[i][1]  + "\n") 
		
		
		printLevel(1,options, maskingString)
		sys.stderr.write(maskingString)
		
		return offsets,proteinList



	def createTEIRESIASinputfromProteinList(self,proteinList,options):
		outString = ""
		for proteins in proteinList.keys():
			outString += ">" + proteins + " 1\n" + proteinList[proteins] + "\n"

		open(options["TEIRESIAS_input"],"w").write(outString)
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class RunBLAST:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  RunBLAST

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   runBLAST( proteinData,options)
	   @param  proteinData
	   @param  options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################CREATING CLUSTERS#######################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		pass
	#--------------------------------------------------------------------------------------------------------------------#
	def runBLAST(self, proteinData,options):
		blast = BLAST()
		blastParser = parseBLAST()
		
		if options["run_formatDB"] == "T":
			printLevel(1,options,"Formatting BLAST database")
			sys.stderr.write(time.strftime('%X') + ":\n" + "Formatting BLAST database\n")

			blast.createDatabase(options["filename"],proteinData)
				
		blastHits = {}
		blastGlobalHits = {}
		blastGlobalAlignments = {}


		for protein in proteinData:
			open(options["input_path"] + "/query/" + protein + ".seq","w").write(">" + protein + "\n" + proteinData[protein])
			
			if options["run_BLAST"] == "T":
				printLevel(1,options,"Running BLAST searches")
				sys.stderr.write(time.strftime('%X') + ":\n" + "Running BLAST searches\n")
				printLevel(1,options,"Blasting " + protein + " against Database\r",)
				sys.stderr.write(time.strftime('%X') + ":\t" + "Blasting " + protein + " against Database\n")
				blast.queryDatabase(options["filename"],protein,proteinData,options["eVal"])

			printLevel(1,options,options["BLAST_results"]  +protein + ".blastp")
			blastHits[protein] = blastParser.parseBLASTFile(options["BLAST_results"]  + protein + ".blastp", options)
			[blastGlobalHits[protein],blastGlobalAlignments[protein]] = blast.globalAlignment(protein,proteinData)

		
		printLevel(1,options,"Finished BLAST searches\t\t\t\n")
		sys.stderr.write("\n")
			

		sys.stderr.write(time.strftime('%X') + ":\n" + "Finished BLAST searches\n")

		return [blastHits,blastGlobalHits,blastGlobalAlignments]
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class BLAST:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  BLAST

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   createDatabase(filename)
	   @param filename

	   queryDatabase(filename,protein,proteinData,eValue)
	   @param filename
	   @param protein
	   @param proteinData
	   @param eValue

	   collapseList(listTemp)
	   @param listTemp

	   globalAlignment(protein)
	   @param protein

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass
	#--------------------------------------------------------------------------------------------------------------------#	
	def createDatabase(self,filename,proteinData):
		temp = open(filename + ".temp","w")

		for protein in proteinData:
			temp.write(">" + protein + "\n" + proteinData[protein] + "\n")

		temp.close()

		filesuffix = os.path.basename(options["filename"]).split(".")[0]

		args = ["formatdb","-p","T","-i",filename + ".temp","-n",options["input_path"] + "/database/" +filesuffix + "db"]
		blastDir = options["BLAST_path"] + "formatdb"


		try:
			os.spawnv(os.P_WAIT,blastDir,args)
		except:
			errorBLAST_not_found ="ERROR : formatDB cannot be found at: " +  options["BLAST_path"] + ".\n" + "-"*66
			print errorBLAST_not_found
			sys.exit()

		os.remove(filename + ".temp")
	#--------------------------------------------------------------------------------------------------------------------#	
	def queryDatabase(self,filename,protein,proteinData,eValue):
		#open(options["input_path"] + "/query/" + protein + ".seq","w").write(">" + protein + "\n" + proteinData[protein])
		filesuffix = os.path.basename(options["filename"]).split(".")[0]
		
		args = ["blastall","-p","blastp","-d",options["input_path"].replace("\\\\","/") + "/database/" + filesuffix + "db","-i",options["input_path"] + "/query/"  +protein + ".seq","-e",str(float(eValue)),"-m","0","-o",options["input_path"] + "/results/" + protein + ".blastp","-F","F"]
		
 		blastDir = options["BLAST_path"] + "blastall"
 		#print blastDir
		try:
			
			if options['overwrite'] == 'T':
				#print args
				os.spawnv(os.P_WAIT,blastDir,args)		
			else:
				if os.path.exists(options["input_path"] + "/results/" + protein + ".blastp"):
					printLevel(1,options,'Skipping ' + protein + ', blast results file already exists.')
					pass
				else:
					os.spawnv(os.P_WAIT,blastDir,args)
		except:
			errorBLAST_not_found ="ERROR blastall cannot be found at: " +  options["BLAST_path"] + ".\n" + "-"*66
			print errorBLAST_not_found
			sys.exit()

	#--------------------------------------------------------------------------------------------------------------------#
	def collapseList(self,listTemp):
			stringTemp = ""
			for values in listTemp:
				stringTemp += values

			return stringTemp
	#--------------------------------------------------------------------------------------------------------------------#
	def globalAlignment(self,protein,proteinData):
		
		htmlFile = open( options["input_path"] + "/"+"alignments/" + protein + ".alignment.html","w")
		htmlFile.write("<html><PRE>")
		htmlFile.write("<h2>Alignments for " + protein + "<br></h2><h3>Cut-off : \t" + str(options["eVal"]) + "</h3><br>") 

		blast = rje_blast.BLASTRun(cmd_list=['i=-1', 'v=0'])

		blastFile = options["input_path"] + "/results/" + protein + ".blastp"
		
		blast.readBLAST(resfile=blastFile,clear=True)

		
		tempFile = open(options["output"].replace(".rank",".alignments.html"),"a")
		sequence = open(options["input_path"] + "/query/" + protein + ".seq","r").read().split("\n")[1]

		tempString = ""
		globalHits = {}
		globalAlignments = {}
		printBlast = 1

		tempString += "<FONT SIZE=2 face='courier'>Query: " + protein + "<br></FONT><FONT SIZE=1 color='#ff0000' face='courier'>" 
		tempString +=  "<a href='" + "alignments/" +  protein + ".alignment.html" + "'>Graphical View</a></font><br>"
		count = 0

		
		for hit in blast.search[0].hit:
			
			tempFile.write("<FONT SIZE=1 face='courier'>")
			
			tempString += "<BR>Hit: %-20s"%(hit.info['Name'])[1:21] + "<br>"  + "%-5s"%str(hit.stat['E-Value']) + "<br>" + proteinData[protein]  + "<br>"
			htmlFile.write("<h3>" + hit.info['Name']+ "</h3>\n")
			
			check = 0
			count += 1

			matchList = list("x"*len(sequence))
			offsets = []

			for aln in hit.aln:
				bit_score = aln.stat['BitScore']
				queryStart = aln.stat['QryStart']
				subjectLen = len(aln.info['SbjSeq'])
				subject = aln.info['SbjSeq']
				subjectStart = aln.stat['SbjStart']
				match =  aln.info['AlnSeq']
				query = aln.info['QrySeq']
				

				htmlFile.write("<FONT SIZE=1 face='courier'>")
				htmlFile.write("subject : match<br>score  : %s"%bit_score + "<br>e-value: %s"%bit_score + "<BR>")
				htmlFile.write("Query location: %s to %s" % ( queryStart,queryStart + subjectLen) + "<BR>")
				htmlFile.write("Target location: %s to %s" % ( subjectStart,subjectStart + subjectLen) + "<BR>")
				htmlFile.write("\n" + subject.replace(" ",".") + "<BR>")
				htmlFile.write("\n" + match.replace(" ",".") + "<BR>" + query.replace(" ",".")  + "<BR><BR></FONT>")
			 
				offsets.append([[queryStart,queryStart + subjectLen],[subjectStart,subjectStart + subjectLen]])
										
				temp1 = "."*(queryStart - 1)
				temp2 = "."*(queryStart - 1)
				gap = 0
				
				for x in range(subjectLen):
					if query[x] == "-":
						gap += 1
					else:
						temp1 += match[x]
						if match[x] != " ":
							if matchList[x + queryStart - 1 - gap] == "x":
								matchList[x + queryStart - 1 - gap] = match[x]
						else:
							if matchList[x + queryStart - 1 - gap] == "x":
								matchList[x + queryStart - 1 - gap] = "."
						temp2 += query[x]

				

			data = [[protein,hit.info['Name']],[len(proteinData[protein]),len(proteinData[hit.info['Name']])],offsets]
	 		
			try:
				HomologyPlotter().plotHomology(data,options)
			except:
				#print 'Skipping homology plot'
				pass
			
			htmlFile.write("<p><img src='" + os.path.abspath("./" + options["input_path"] + "./alignments/"  + protein + "_" + hit.info['Name'] + ".alignment.gif'") + "></p><br>\n")
		
			alignments = self.collapseList(matchList)
			mismatch = self.collapseList(matchList).count("x")
			partialmatch = self.collapseList(matchList).count("+")
			fullmatch = (len(self.collapseList(matchList))  - (self.collapseList(matchList).count(".") + self.collapseList(matchList).count("x") + self.collapseList(matchList).count("+")))
			globalHits[hit.info['Name']] = float(fullmatch)/len(sequence)
			
			globalAlignments[hit.info['Name']] = self.collapseList(matchList).replace("x",".")

			tempString +=  self.collapseList(matchList) + "\n"

			tempFile.write(tempString + "<br>")
			tempString = ""
			
		del blast
		tempFile.write("<br><br>------------------------------------------------------------------------<br><br>")
		tempFile.close()
		htmlFile.write("</PRE></html>")
		
		return [globalHits,globalAlignments]
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class parseBLAST:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  parseBLAST

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   parseBLASTFile(filename,options)
	   @param filename
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass 
	#--------------------------------------------------------------------------------------------------------------------#
	def parseBLASTFile(self,filename,options):
		
		try:
			fileIn = open(filename,"r")	
		except:
			errorBLAST_file_not_found = "\n\nERROR:\n" + "+"*70 
			errorBLAST_file_not_found += "\nBLAST output file does not exist.\n"
			errorBLAST_file_not_found +="Try running AERPIMP again using the -BT and -DT options.\n" + "+"*70
			print errorBLAST_file_not_found
			sys.exit()

		data = fileIn.read()

		pattern = re.compile("Sequences producing significant alignments:[^>]*")

		hitDict = {}
		try:
			for lines in pattern.findall(data)[0].split("\n")[2:-2]:
				hitDetails =lines.split()

				if hitDetails[2][0] == "e":
					eValue = "1" + hitDetails[2]
				else:
					eValue = hitDetails[2]

				if float(eValue) < options["eVal"]:
					hitDict[hitDetails[0]] = [hitDetails[1],float(eValue)]
			
		except Exception,e:
		
			errorBLAST_file_Long = "\n\nERROR:\n" + "+"*70 
			errorBLAST_file_Long += "\nBLAST output file does not contain any data.\n"
			if options["long_format"] == 'F':
				errorBLAST_file_Long +="Try running SLimdisc again using the -LT option.\n" + "+"*70
				
			print errorBLAST_file_Long
			sys.exit()
	
		
		del data
		del hitDetails
		fileIn.close()
		return hitDict
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Teiresias:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Teiresias

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   teiresias(options)
	   @param options
	#################################################################################################
	"""
	#################################################################################################	

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################CREATING PATTERNS#######################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		pass
	#--------------------------------------------------------------------------------------------------------------------#
	def teiresias(self,options):
		wallTime =  options["TEIRESIAS_walltime"] * 60
		sleepTime = 0.5

		if options["TEIRESIAS_walltime"] == 0:
			self.teiresiasNonThread(options)
		else:
			print time.ctime()
			print 'walltime == ' + str(wallTime) + ' sec'
			startTime = time.time()
			
			teiresiasThread = TeiresiasThread(options)
			teiresiasThread.start()
			
			
			while teiresiasThread.isAlive():
				print 'TEIRESIAS run time : %-1.3f'%(time.time() - startTime) + '\t(Walltime==' + str(options["TEIRESIAS_walltime"] * 60) +  ')\r',
				sys.stdout.flush()
				
				if (time.time() - startTime) >  wallTime:
					print '\n********WALLTIME of ' + str(options["TEIRESIAS_walltime"])  + ' minutes reached ********'
					
					print 'Killing TEIRESIAS ' + str(teiresiasThread.pid)
					
					if sys.platform[0:3] == 'win':
						os.system('taskkill /F /PID ' + str(teiresiasThread.pid))
					else:
						os.system('kill -KILL ' + str(teiresiasThread.pid))
					
					print '**' + '*' * len('******WALLTIME of ' + str(options["TEIRESIAS_walltime"])  + ' minutes reached ********')
					sys.exit()
					
				time.sleep(sleepTime)
		
	#--------------------------------------------------------------------------------------------------------------------#		
	def teiresiasNonThread(self,options):
		if options["TEIRESIAS_equiv"] == 'F':
			equiv = '-p'
		elif options["TEIRESIAS_equiv"] == 'T':
			if os.path.exists(options["calling_folder"] + "/Teiresias/equiv.txt"):
				
				equiv = "-b" + options["calling_folder"] + "/Teiresias/equiv.txt"
				
				printLevel(1,options, "%-49s"%"Teiresias/equiv.txt" + "\t" + "present")
			else:
				print os.environ 
				equiv = '-p'
				print "#"*50 +'\nTeiresias equiv file path does not exist\n' + "Teiresias/equiv.txt" + '\nRunning without ambiguity\n' + "#"*50
				
		elif len(options["TEIRESIAS_equiv"]) > 1:
			
			if os.path.exists(options["TEIRESIAS_equiv"]):
				equiv = '-b' + options["TEIRESIAS_equiv"] 
				printLevel(1,options, "%-49s"%equiv + "\t" + "present")
			else:
				equiv = '-p'
				print "#"*50 +'\nTeiresias equiv file path does not exist\n' + options["TEIRESIAS_equiv"]+ '\nRunning without ambiguity\n' + "#"*50
			
		tempFile = tempfile.mktemp(".txt") 
		

		if options["TEIRESIAS_local(IBM_path_length_bug_fix)"] == 'T':
			args = [options["TEIRESIAS_path"] + 'teiresias_char',"-v","-r","-i" + options["TEIRESIAS_input"],"-n2","-o"+tempFile,"-l"+str(options["TEIRESIAS_fixedPositions"]),"-w"+str(options["TEIRESIAS_patternLength"]),"-c1","-k"+str(options["TEIRESIAS_supportString"]),"-p",equiv]
		else:
			args = [options["TEIRESIAS_path"] + 'teiresias_char',"-v","-r","-i" + options["TEIRESIAS_input"],"-n2","-o"+options["TEIRESIAS_output"],"-l"+str(options["TEIRESIAS_fixedPositions"]),"-w"+str(options["TEIRESIAS_patternLength"]),"-c1","-k"+str(options["TEIRESIAS_supportString"]),"-p",equiv]
		
		try:
			if options['overwrite'] == 'T':
				printLevel(1,options, "Running TEIRESIAS")
				sys.stderr.write("Running TEIRESIAS\n")
				#print options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' ')
				#info = os.popen(options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' '))
				processId = os.spawnv(os.P_NOWAIT,args[0],args)
				self.pid  = processId
				#print self.pid
				#print args
				os.waitpid(processId,0)
				
				
				if options['TEIRESIAS_local(IBM_path_length_bug_fix)']  == 'T':
					try:
						os.remove(options["TEIRESIAS_output"])
					except:
						pass
					
					try:
						shutil.copyfile(tempFile,options["TEIRESIAS_output"])
						os.remove(tempFile)
					except:
						print 'Unable to find temporary TEIRESIAS files'
	
			else:
				if os.path.exists(options["TEIRESIAS_output"]):
					print 'TEIRESIAS output already exists'
					pass
				else:
					printLevel(1,options, "Overwriting TEIRESIAS output")
					sys.stderr.write("Running TEIRESIAS\n")
					processId = os.spawnv(os.P_NOWAIT,args[0],args)
					self.pid  = processId
					os.waitpid(processId,0)
	
					print options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' ')
	
					if options['TEIRESIAS_local(IBM_path_length_bug_fix)']  == 'T':
						print tempFile,options["TEIRESIAS_output"]
					
						try:
							os.remove(options["TEIRESIAS_output"])
						except:
							pass
							
						try:
							shutil.copyfile(tempFile,options["TEIRESIAS_output"])
							os.remove(tempFile)
						except:
							print 'Unable to find temporary TEIRESIAS files'
	
			#os.spawnv(os.P_WAIT,options["TEIRESIAS_path"] + "teiresias_char",args)
		except:
			errorBLAST_not_found ="ERROR : while running TEIRESIAS "
			print errorBLAST_not_found
	#----------------------------------------------END OF CLASS----------------------------------------------------------#


class TeiresiasThread(Thread):
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Teiresias

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   teiresias(options)
	   @param options
	#################################################################################################
	"""
	#################################################################################################	

	def __init__ (self,option):
		Thread.__init__(self)
		self.status = -1
     		self.pid = 0
     		self.options = option
     		
	#--------------------------------------------------------------------------------------------------------------------#
	def run(self):
		self.time = 20
		self.teiresias()
	#--------------------------------------------------------------------------------------------------------------------#
	def teiresias(self):
	
		if options["TEIRESIAS_equiv"] == 'F':
			equiv = '-p'
		elif options["TEIRESIAS_equiv"] == 'T':
			if os.path.exists("Teiresias/equiv.txt"):
				equiv = "-bTeiresias/equiv.txt"
				printLevel(1,options, "%-49s"%"Teiresias/equiv.txt" + "\t" + "present")
			else:
				equiv = '-p'
				print "#"*50 +'\nTeiresias equiv file path does not exist\n' + "Teiresias/equiv.txt" + '\nRunning without ambiguity\n' + "#"*50
			
		elif len(options["TEIRESIAS_equiv"]) > 1:
			
			if os.path.exists("Teiresias/equiv.txt"):
				equiv = '-b' + options["TEIRESIAS_equiv"] 
				printLevel(1,options, "%-49s"%equiv + "\t" + "present")
			else:
				equiv = '-p'
				print "#"*50 +'\nTeiresias equiv file path does not exist\n' + options["TEIRESIAS_equiv"] + '\nRunning without ambiguity\n' + "#"*50
			
		tempFile = tempfile.mktemp(".txt") 
		

		if options["TEIRESIAS_local(IBM_path_length_bug_fix)"] == 'T':
			args = [options["TEIRESIAS_path"] + 'teiresias_char',"-v","-r","-i" + options["TEIRESIAS_input"],"-n2","-o"+tempFile,"-l"+str(options["TEIRESIAS_fixedPositions"]),"-w"+str(options["TEIRESIAS_patternLength"]),"-c1","-k"+str(options["TEIRESIAS_supportString"]),"-p",equiv]
		else:
			args = [options["TEIRESIAS_path"] + 'teiresias_char',"-v","-r","-i" + options["TEIRESIAS_input"],"-n2","-o"+options["TEIRESIAS_output"],"-l"+str(options["TEIRESIAS_fixedPositions"]),"-w"+str(options["TEIRESIAS_patternLength"]),"-c1","-k"+str(options["TEIRESIAS_supportString"]),"-p",equiv]
		
		try:
			if self.options['overwrite'] == 'T':
				printLevel(1,self.options, "Running TEIRESIAS")
				sys.stderr.write("Running TEIRESIAS\n")
				print self.options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' ')
				#info = os.popen(self.options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' '))
				processId = os.spawnv(os.P_NOWAIT,args[0],args)
				self.pid  = processId
				#print self.pid
				#print args
				os.waitpid(processId,0)
				
				
				if self.options['TEIRESIAS_local(IBM_path_length_bug_fix)']  == 'T':
					try:
						os.remove(self.options["TEIRESIAS_output"])
					except:
						pass
					
					try:
						shutil.copyfile(tempFile,self.options["TEIRESIAS_output"])
						os.remove(tempFile)
					except:
						print 'Unable to find temporary TEIRESIAS files'
	
			else:
				if os.path.exists(self.options["TEIRESIAS_output"]):
					print 'TEIRESIAS output already exists'
					pass
				else:
					printLevel(1,self.options, "Overwriting TEIRESIAS output")
					sys.stderr.write("Running TEIRESIAS\n")
					processId = os.spawnv(os.P_NOWAIT,args[0],args)
					self.pid  = processId
					os.waitpid(processId,0)
	
					print self.options["TEIRESIAS_path"] + "teiresias_char " + rje.join(args[1:],' ')
	
					if self.options['TEIRESIAS_local(IBM_path_length_bug_fix)']  == 'T':
						print tempFile,self.options["TEIRESIAS_output"]
					
						try:
							os.remove(self.options["TEIRESIAS_output"])
						except:
							pass
							
						try:
							shutil.copyfile(tempFile,self.options["TEIRESIAS_output"])
							os.remove(tempFile)
						except:
							print 'Unable to find temporary TEIRESIAS files'
	
			#os.spawnv(os.P_WAIT,self.options["TEIRESIAS_path"] + "teiresias_char",args)
		except:
			errorBLAST_not_found ="ERROR : while running TEIRESIAS "
			print errorBLAST_not_found
	
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class ParseTeiresias:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  ParseTeiresias

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   parse(filename,options)
	   @param filename
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################PARSING PATTERNS########################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		self.dictPatterns = {}
		self.dictOccurances = {}
	#--------------------------------------------------------------------------------------------------------------------#
	def parse(self,filename,options):
		self.dictPatterns.clear()
		start = time.time()	
		pattern2 = re.compile('\[[^\]]*\]|.')

		scorer = Score()
		complexity = Complexity()

		aminoAcidOcc = scorer.loadaminoAcidOccurances(options["aminoAcidOccuranceFile"])

		try:
			input = open(filename,'r')
		except:
			errorTEIRESIAS_file_not_found = "\nERROR:\n" + "+"*70 + "\nTEIRESIAS output file does not exist.\n"
			
			if options["run_TEIRESIAS"] == 'F':
				errorTEIRESIAS_file_not_found += "Try running SLiMDisc again using the -TT option.\n" 
				
			errorTEIRESIAS_file_not_found += "+"*70
			
			print errorTEIRESIAS_file_not_found
			sys.exit()
			
		pattern1 = re.compile('\S+')
		count = 0
		counter = 0
		
		input.seek(0,2)
		length = input.tell()
		input.seek(0)

		printLevel(1,options,"\nProcessing TEIRESIAS output")
		printLevel(1,options,"%-10s"%"Input Size" + "\t%-10s"%"Finished" )
			
		while input.tell() < length:
			data = input.readlines(10000)
		
			for lines in data:
				counter += 1
				#print lines
				if counter%1000 == 0:	
					if options['verbosity'] > 1:
						print "%-10s"%(str(length/1000) + "kb") + "\t%-10s"%str(str(input.tell()/1000) + "kb") + "\r", 
					
				
				dictTemp = {}
				
				if lines[0] != "#":
					data = pattern1.findall(lines)
					
					number_of_times = data.pop(0)
					number_of_unique_times = data.pop(0)

					pattern = data.pop(0)

					self.dictOccurances[pattern] = number_of_times
		

					while len(data) > 1:
						seq = data.pop()
						offset = data.pop() 

						if offset in dictTemp:
							dictTemp[offset].append(seq)
						else:
							dictTemp[offset] = [seq]



						
					information_content = scorer.calculateScore(pattern,aminoAcidOcc)		
	
					patternParts = pattern2.findall(pattern)
#


					if options["use_filtering"] == "T":
						if information_content > float(options["information_cutOff"]) and len(patternParts) <= int(options["TEIRESIAS_patternLength"]) and len(patternParts) > 2 and complexity.complexity(pattern) > float(options["complexity_cutOff"]):
							count += 1
							self.dictPatterns[pattern] = [dictTemp,information_content]
					else:
						count += 1
						self.dictPatterns[pattern] = [dictTemp,information_content]
					
		if options['verbosity'] > 1:
			print "%-10s"%(str(length/1000) + "kb") + "\t%-10s"%str(str(input.tell()/1000) + "kb") + "\r",
					
		try:		
			outString = "\n" 
	
			outString +=  "Patterns processed            \t:" + str(counter - 5) + "\n"
			outString +=  "Patterns kept                 \t:" + str(count) + "\n"
			outString +=  "Percentage processed          \t:" "%2.2f"%(float(count)/(counter -5)) + "%\n"
			outString +=  "Time			     \t:" + "%2.2f seconds"%float(time.time() - start) + "\n"
			
			printLevel(1,options, "\n\nInput Statistics" + outString)
			sys.stderr.write("Input Statistics" + outString  + "\n")
		except:
			print 'Problem with parsing of TEIRESIAS input'
			
		
		del data
		input.close()
		return [self.dictPatterns,self.dictOccurances]
	#----------------------------------------------END OF CLASS----------------------------------------------------------#

class ParseTeiresiasFull:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  ParseTeiresias

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   parse(filename,options)
	   @param filename
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################PARSING PATTERNS########################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		self.dictPatterns = {}
		self.dictOccurances = {}
	#--------------------------------------------------------------------------------------------------------------------#
	def parse(self,filename,offsets,options):
		self.dictPatterns.clear()
		start = time.time()	
		pattern2 = re.compile('\[[^\]]*\]|.')

		scorer = Score()
		complexity = Complexity()

		aminoAcidOcc = scorer.loadaminoAcidOccurances(options["aminoAcidOccuranceFile"])

		try:
			input = open(filename,'r')
		except:
			errorTEIRESIAS_file_not_found = "\nERROR:\n" + "+"*70 + "\nTEIRESIAS output file does not exist.\nTry running AERPIMP again using the -TT option.\n" + "+"*70
			print errorTEIRESIAS_file_not_found
			sys.exit()

		pattern1 = re.compile('\S+')
		count = 0
		counter = 0
		
		input.seek(0,2)
		length = input.tell()
		input.seek(0)

		printLevel(1,options,"\nProcessing TEIRESIAS output")
		printLevel(1,options,"%-10s"%"Input Size" + "\t%-10s"%"Finished" )
			
		while input.tell() < length:
			data = input.readlines(10000)
		
			for lines in data:
				counter += 1

				if counter%1000 == 0:
					if options['verbosity'] > 1:
						print "%-10s"%(str(length/1000) + "kb") + "\t%-10s"%str(str(input.tell()/1000) + "kb") + "\r",
					
				
				dictTemp = {}
				
				if lines[0] != "#":
					
					data = pattern1.findall(lines)
					
					number_of_times = data.pop(0)
					number_of_unique_times = data.pop(0)

					pattern = data.pop(0)

					self.dictOccurances[pattern] = number_of_times
		

					while len(data) > 1:
						
						seq = data.pop()
						offset = data.pop() 
					
						if offsets[int(offset)][2] in dictTemp:
							dictTemp[offsets[int(offset)][2]].append(int(seq) + int(offsets[int(offset)][1]))
						else:
							dictTemp[offsets[int(offset)][2]] = [int(seq) + int(offsets[int(offset)][1])]

				
					information_content = scorer.calculateScore(pattern,aminoAcidOcc)		
	
					patternParts = pattern2.findall(pattern)

					if options["use_filtering"] == "T":
				
						if information_content > float(options["information_cutOff"]) and len(patternParts) <= int(options["TEIRESIAS_patternLength"]) and len(patternParts) > 2 and complexity.complexity(pattern) >= float(options["complexity_cutOff"]):
							count += 1
						
							self.dictPatterns[pattern] = [dictTemp,information_content]
					else:
						count += 1
						self.dictPatterns[pattern] = [dictTemp,information_content]
					

		printLevel(1,options,"%-10s"%(str(length/1000) + "kb") + "\t%-10s"%str(str(input.tell()/1000) + "kb") + "\r",)
					
			
		outString = "\n" 

		outString +=  "Patterns processed            \t:" + str(counter - 5) + "\n"
		outString +=  "Patterns kept                 \t:" + str(count) + "\n"
	
		try:	
			outString +=  "Percentage processed          \t:" "%2.2f"%(float(count)/(counter -5)) + "%\n"
		except:
			pass
			
		outString +=  "Time			     \t:" + "%2.2f seconds"%float(time.time() - start) + "\n"
			
		printLevel(1,options,"\n\nInput Statistics" + outString)
		sys.stderr.write("Input Statistics" + outString  + "\n")
		
		del data
		input.close()
		return [self.dictPatterns,self.dictOccurances]
	#----------------------------------------------END OF CLASS----------------------------------------------------------#


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class aaOccurance:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	   

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		self.aaDict = {'A':0,'R':0,'N':0,'D':0,'C':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'Q':0,'X':0,'B':0,'Z':0}
	#--------------------------------------------------------------------------------------------------------------------#

	def resetDict(self):
		self.aaDict = {'A':0,'R':0,'N':0,'D':0,'C':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'Q':0,'X':0,'B':0,'Z':0}
	#--------------------------------------------------------------------------------------------------------------------#
	def countAAbyAcc(self,options,proteinData):
		self.resetDict()
		count = 0

		for protein in proteinData.keys():
			sequence = proteinData[protein]
			count += len(sequence)

			try:
				for aminoAcids in sequence:
					self.aaDict[aminoAcids] += 1
			except:
				pass
		return count
	#--------------------------------------------------------------------------------------------------------------------#
	def createAAOccuranceFilefromProteinList(self,options,proteinData):
		count = self.countAAbyAcc(options,proteinData)
		self.printOutput(options,count)
	#--------------------------------------------------------------------------------------------------------------------#
	def printOutput(self,options,count):

		output = open(options["aminoAcidOccuranceFile"],'w')
		
		for aminoAcids in self.aaDict:
			try:
				output.write(aminoAcids + " %1.3f\n"%(float(self.aaDict[aminoAcids])/count))
			except:
				print 'Error calculating amino acid probabilities'
				sys.exit()
				
		output.close()
	#----------------------------------------------END OF CLASS----------------------------------------------------------#
		
		
#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Filereader:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Filereader

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   createTeiresiasIndexList(filename)
	   @param filename

	   readTeiresiasData(options)
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass
	#--------------------------------------------------------------------------------------------------------------------#
	def createTeiresiasIndexList(self,filename):
		indexList = {}
		input = open(filename,'r')

		count = 0
		for lines in input.readlines():
			if(lines[0] == ">"):
				
				acc = lines.split(" ")
				acc[0].replace("\r","")
				indexList[count] = acc[0][1:]
				count += 1

		input.close()
		return indexList
	#--------------------------------------------------------------------------------------------------------------------#
	def readTeiresiasData(self,offset,options):
		printLevel(1,options,"Loading TEIRESIAS file\r",)
		start = time.time()
		sys.stderr.write(time.strftime('%X') + ":\n" + "Loading TEIRESIAS file\n")
		
		teiresiasIndex = {}
			
		x = 1

		if options["long_format"] == "T":	
			teiresiasParser =ParseTeiresiasFull()
			[dictPatterns,dictOccurances] = teiresiasParser.parse(options["TEIRESIAS_output"],offset,options)
			for subSequence in offset:
				teiresiasIndex[subSequence[2]] = subSequence[0]
				
		else:
			teiresiasParser =ParseTeiresias()
			[dictPatterns,dictOccurances] = teiresiasParser.parse(options["TEIRESIAS_output"],options)
			teiresiasIndex = self.createTeiresiasIndexList(options["TEIRESIAS_input"])

		teiresiasData = {}

		printLevel(1,options,"Formatting Pattern Data")
		length = len(dictPatterns.keys())
		count = 0

		printLevel(1,options, "%-10s"%"Patterns" + "\t%-10s"%"Complete" + "\r",)

		for pattern in dictPatterns.keys():
			if count%1000 == 0: 
				if options['verbosity'] > 1:
						print "%-10s"%str(length) + "\t"  +"%-10s"%str(count) + "\r",
			
			count += 1
			
			proteinList = []
			offsetList = []


			for occurances in dictPatterns[pattern][0]:
				proteinList.append(teiresiasIndex[int(occurances)])
				offsetList.append(dictPatterns[pattern][0][occurances])

		
			if len(options["query_protein"]) > 1:
				if options["query_protein"] in proteinList:
					teiresiasData[pattern] = [proteinList,offsetList,dictPatterns[pattern][1]]
				else:
					pass

			else:
				teiresiasData[pattern] = [proteinList,offsetList,dictPatterns[pattern][1]]

			del dictPatterns[pattern]

		printLevel(1,options,"%-10s"%str(count) + "\t"  +"%-10s"%str(length))
		outString = "\n" 

		outString +=  "Patterns processed            \t:" + str(count) + "\n"
		outString +=  "Time			     \t:" + "%2.2f seconds"%float(time.time() - start) + "\n"
		
		printLevel(1,options,"\n\nFormatting Statistics" + outString  +"\n")
		sys.stderr.write("Formatting Statistics" + outString  + "\n")

		printLevel(1,options, "Finished loading TEIRESIAS file\n")
		sys.stderr.write(time.strftime('%X') + ":\n" + "Finished loading TEIRESIAS file\n")


		del teiresiasIndex
		del teiresiasParser
		del dictOccurances
		return teiresiasData
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class SuffixTree:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  SuffixTree

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   add_pattern(tree,pattern,proteinData,pos,count)
	   @param tree
	   @param pattern
	   @param proteinData
	   @param pos
	   @param count

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass
	#--------------------------------------------------------------------------------------------------------------------#
	def add_pattern(tree,pattern,proteinData,pos,count):
		if pattern[count] in tree[pos] and count != len(pattern) - 1:
			temptree = tree[pos]

		
			if pattern[count + 1] in tree[pos][pattern[count]]:
				pass
			else:
				tree[pos][pattern[count]][pattern[count + 1]] = {}

			tree[pos] = add_pattern(temptree,pattern,proteinData,pattern[count],count + 1)
		else:
			
			pass

		if count == len(pattern) - 1:
			if "index" in tree[pos][pattern[-1]]:
				for offsets in proteinData:	
					tree[pos][pattern[-1]]["index"][offsets]= proteinData[offsets]
			else:
				tree[pos][pattern[-1]]["index"] = proteinData
				

		return tree
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Score:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Score

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   collapseList(listTemp)
	   @param listTemp

	   calculateScore(pattern,aaOcc)
	   @param pattern
	   @param aaOcc

	   loadaminoAcidOccurances(aminoAcidOccuranceFile)
	   @param aminoAcidOccuranceFile

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		self.part1 = 0
		self.average = 0

		pass
	#--------------------------------------------------------------------------------------------------------------------#
	def collapseList(self,listTemp):
			stringTemp = ""
			for values in listTemp:
				stringTemp += values

			return stringTemp
	#--------------------------------------------------------------------------------------------------------------------#
	def calculateInformationContent(self,pattern,aaOcc):

		count = 0	
		p1 = re.compile('\[[^\]]*\]')
		p2 = re.compile('(?<=])[^\[]*')

		data = p1.findall(pattern)
		patternTemp = "]" + pattern
		informationContent = 0
		for values in data:
			#print values
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
				
			informationContent +=  (-sum - self.part1)

		count += pattern.count(".")
		data = p2.findall(patternTemp)
		
		stringTmp = self.collapseList(data)
		stringTmp = stringTmp.replace(".","")

		for AA in stringTmp:	
			try:
				informationContent +=  -aaOcc[AA]/aaOcc[AA]*math.log((aaOcc[AA]/aaOcc[AA]),2) - self.part1 
			except:
				print pattern

		return  informationContent - options["gap_weight"]*pattern.count(".")
	#--------------------------------------------------------------------------------------------------------------------#
	def loadaminoAcidOccurances(self,aminoAcidOccuranceFile):
		printLevel(1,options,"\nLoading Peptide Occurances\r",)
		sys.stderr.write("Load Peptide Occurances\n")
		aaOcc = open(aminoAcidOccuranceFile,"r")
		
		aminoAcidOccurance = {}

		for AA in aaOcc.readlines():
			AA = AA.strip()
			data = AA.split(" ")
			aminoAcidOccurance[data[0]] = float(data[1])

		printLevel(1,options,"Peptide Occurances Loaded\t")
		aaOcc.close()

		for aa in aminoAcidOccurance:
			if aminoAcidOccurance[aa] > 0:
				self.part1 +=  aminoAcidOccurance[aa]*math.log((aminoAcidOccurance[aa]),2)
			
		return aminoAcidOccurance
	#--------------------------------------------------------------------------------------------------------------------#
	def calculateScore(self,pattern,aaOcc):
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
				
			informationContent +=  (-sum - self.part1)/count

		count += pattern.count(".")
		data = p2.findall(patternTemp)
		
		stringTmp = self.collapseList(data)
		stringTmp = stringTmp.replace(".","")

		for AA in stringTmp:	
			try:
				informationContent +=  -aaOcc[AA]/aaOcc[AA]*math.log(1 + aaOcc[AA] - 0.05,2) - self.part1 
			except:
				pass

		return  (informationContent - options["gap_weight"]*pattern.count("."))
	#--------------------------------------------------------------------------------------------------------------------#
	def loadaminoAcidOccurances(self,aminoAcidOccuranceFile):
		printLevel(1,options, "\nLoading Peptide Occurances\r",)
		sys.stderr.write("Load Peptide Occurances\n")
		aaOcc = open(aminoAcidOccuranceFile,"r")
		
		aminoAcidOccurance = {}

		for AA in aaOcc.readlines():
			AA = AA.strip()
			data = AA.split(" ")
			aminoAcidOccurance[data[0]] = float(data[1])

		printLevel(1,options,"Peptide Occurances Loaded\t")
		aaOcc.close()

		for aa in aminoAcidOccurance:
			if aminoAcidOccurance[aa] > 0:
				self.part1 +=  aminoAcidOccurance[aa]*math.log((aminoAcidOccurance[aa]),2)
		
		self.average = sum(aminoAcidOccurance.values())/20

		return aminoAcidOccurance
	#----------------------------------------------END OF CLASS----------------------------------------------------------#


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class SurfaceProbability:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  SurfaceProbability

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	
	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		self.surface = ""
		self.janin = {'Z':0.00,'B':0.00,'X':0.00,'A':0.49,'R':0.95,'N':0.81,'D':0.78,'C':0.26,'Q':0.84,'E':0.84,'G':0.48,'H':0.66,'I':0.34,'L':0.40,'K':0.97,'M':0.48,'F':0.42,'P':0.75,'S':0.65,'T':0.70,'W':0.51,'Y':0.76,'V':0.36}
		self.average = sum(self.janin.values())/len(self.janin)
		
	#--------------------------------------------------------------------------------------------------------------------#
	def eminiMethod(self,sequence,pattern,portionPos,portionEnd,options):
		self.surface = ""
		scores = []
		count = 0
		countGreater = 0
		average = 0.0
		position = 0

		for i in range(portionPos,portionEnd):
			if pattern[i - portionPos] == '.':
				pass
			else:
				count += 1
				product = 0.0
				divisor = 0
				for j in range(1,7):
					if i + 4 - j < len(sequence) and i + 4 - j > 0:
						divisor += 1
						if product == 0.0:
							product = self.janin[sequence[i + 4 - j]]
						else:
							try:
								product *= self.janin[sequence[i + 4 - j]]*1.61
							except:
								print sequence
				
				if divisor != 6:
					score = product*self.average*pow(1.61,6 - (divisor))
				else:
					score = product

				scores.append(score)

				if  score > options["surface_cutOff"]:
					self.surface += sequence[i]
					countGreater += 1
				else:
					self.surface += "X"

				average += score
			
			position += 1
		
		percentageSurface = float(countGreater)/count
		averageSurfaceProb = average/count

		surfaceString = self.surface

		return {'averageSurfaceProb':averageSurfaceProb,'percentageSurface':percentageSurface,'surfaceString':surfaceString,"scores":scores}
	#--------------------------------------------------------------------------------------------------------------------#
	def checkPattern(self,sequence,pattern,offset,options):
		printLevel(1,options,sequence)
		for offsets in offset:
			eminiResult = self.eminiMethod(sequence,pattern,int(offsets),int(offsets) + len(pattern),options)
			
		eminiResult['pattern'] = pattern
			
	
		return [eminiResult["scores"],eminiResult['percentageSurface']]
	#--------------------------------------------------------------------------------------------------------------------#
	def checkSequence(self,sequence,options,protein):
		self.surface = ""
		scores = []
		count = 0
		countGreater = 0
		average = 0.0
		position = 0

		for i in range(0,len(sequence)):
				count += 1
				product = 0.0
				divisor = 0
				for j in range(1,7):
					if i + 4 - j < len(sequence) and i + 4 - j > 0:
						divisor += 1
						try:
							if product == 0.0:
								product = self.janin[sequence[i + 4 - j]]
							else:
								
								product *= self.janin[sequence[i + 4 - j]]*1.61
						except:
							print sequence,'<br>----',sequence[i + 4 - j]
							sys.exit()
				
				if divisor != 6:
					score = product*self.average*pow(1.61,6 - (divisor - 1))
				else:
					score = product

				scores.append(score)

				if  score > options["surface_cutOff"]:
					self.surface += sequence[i]
					countGreater += 1
				else:
					self.surface += "X"

				average += product
		
		try:
			SurfaceAccessibilityPlotter().plotSurfaceAccessibility(protein,scores,options)
		except:
			#print 'Skipping surface accessability plot'
			pass
			
		return self.surface
	#--------------------------------------------------------------------------------------------------------------------#
	def checkPosition(self,sequence,position,options):
		self.surface = ""
		i = position - 1

		for j in range(1,7):
			if i + 4 - j < len(sequence) and i + 4 - j > 0:
				if j == 1:
					product = self.janin[sequence[i + 4 - j]]
				else:
					product *= self.janin[sequence[i + 4 - j]]
					
		if  product*17.6 > options["surface_cutOff"]:
			self.surface += sequence[i]
		else:
			self.surface += "X"

		printLevel(1,options,self.surface)
		printLevel(1,options,sum(self.janin.values()))
	#--------------------------------------------------------------------------------------------------------------------#
	def checkPatterns(self,sequences,pattern,proteins,options):
		
		pattern1 = re.compile('\[[A-Z.]*\]|[A-Z.]')
		pattern =  pattern1.findall(pattern)

		results = {}
		for protein in proteins.keys():
				
			results[protein] = self.checkPattern(sequences[protein],pattern,proteins[protein],options)	
					
		length = len(pattern) - pattern.count(".")

		#######################
		#TODO ambiguous
		#######################

		accessibility = {}
		for i in range(0,length):
			if pattern[i] == ".":
				x =1
			else:
				accessibility[i] = []
				for protein in proteins.keys():
						if results[protein][0][i] > options["surface_cutOff"]:
							accessibility[i].append(1)
						else: 
							accessibility[i].append(0)


		accessibilityPercentage = []
		sumAverage = 0
		for vals in accessibility.keys():
			average = float(sum(accessibility[vals]))/len(accessibility[vals])
			accessibilityPercentage.append(average%1)
			sumAverage += average
		
		percentage_residues_accesible = float(sumAverage)/len(accessibility)

		correlation = float(accessibilityPercentage.count(0))/len(accessibility)
		return [percentage_residues_accesible,sum(accessibility[vals])]#correlation

	#----------------------------------------------END OF CLASS----------------------------------------------------------#


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Complexity:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Complexity

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   countChar(inp)
	   @param inp

	   collapseList(listTemp)
	   @param listTemp

	   factorial(inp)
	   @param inp

	   complexity(input)
	   @param input

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		results = {}
	#--------------------------------------------------------------------------------------------------------------------#
	def countChar(self,inp):
		result = {'A': 0,'R': 0,'N': 0,'D': 0,'C': 0,'Q': 0,'E': 0,'G': 0,'H': 0,'I': 0,'L': 0,'K': 0,'M': 0,'F': 0,'P': 0,'S': 0,'T': 0,'W': 0,'Y': 0,'V': 0, 'B': 0,'Z': 0,'length' : 0}

		for i in range(0,len(inp)):
			try:
				if inp[i] != "." or inp[i] != "[" or inp[i] != "]":
					result[inp[i]] += 1
					result["length"] += 1
			except:
				x = 0
	
		temp = len(inp)
	
		for i in result.keys():
			if i  == "length":
				x = 0
			else:
				result[i] = int((100*(float(result[i])/result["length"])))
		return result
	#--------------------------------------------------------------------------------------------------------------------#
	def collapseList(self,listTemp):
		stringTemp = ""
		for values in listTemp:
			stringTemp += values

		return stringTemp
	#--------------------------------------------------------------------------------------------------------------------#
	def factorial(self,inp):
		out = 1

		for i in range(1,inp + 1):
			out = i*out

		return out
	#--------------------------------------------------------------------------------------------------------------------#
	def complexity(self,input):
		try:
			p = re.compile('(?<=])[^\[]*')
			data = p.findall("]" + input)
			input = self.collapseList(data)

			acc = 1
			counter = []
			result = self.countChar(input)

			for i in result.keys():
				if i  == "length":
					x = 0
				else:
					acc = self.factorial(result[i]) * acc
					
			stringLength = 100

			x = self.factorial(stringLength)
			return 1/float(stringLength)*((math.log(float(x)/acc)/math.log(22)))
		except:
			return 0
	#----------------------------------------------END OF CLASS----------------------------------------------------------#


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Prim:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Prim

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   findEdges(proteinList,protein)
	   @param proteinList
	   @param protein

	   sortMinQueue(min_queue,highScoreList)
	   @param min_queue
	   @param highScoreList

	   prim(positions,blastData,proteinList)	   
	   @param positions
	   @param blastData
	   @param proteinList

	   scorePrim(blastData,proteinList)
	   @param blastData
	   @param proteinList

	#################################################################################################
	"""
	#################################################################################################
	
	#--------------------------------------------------------------------------------------------------------------------#
	def findEdges(self,proteinList, protein1):
		edgeList = []
		for protein2 in proteinList:
			if protein2 != protein1:
				edgeList.insert(0,protein2)
		return edgeList
	#--------------------------------------------------------------------------------------------------------------------#
	def sortMinQueue(self,min_queue,highScoreList):
		queueLength = len(min_queue)
		for i in range(0,queueLength):
			for j in range(0,queueLength):
				if highScoreList[min_queue[i]] < highScoreList[min_queue[j]]:
					tempProtein = min_queue[i]
					min_queue[i] = min_queue[j]
					min_queue[j] = tempProtein
		return min_queue
	#--------------------------------------------------------------------------------------------------------------------#
	def prim(self,positions,blastData,proteinList):
		INFINITY = pow(2,16)
		
		minimumSpanningTree =["empty"]*len(positions)

		highScoreList = {}
		for values in proteinList:
			highScoreList[values] = INFINITY

		min_queue=[]
		for protein in proteinList:
			min_queue.append(protein)

		highScoreList[proteinList[0]] = 0
		self.sortMinQueue(min_queue,highScoreList)    
		
		while len(min_queue) > 0:
			protein1 = min_queue[0]
			min_queue.remove(min_queue[0])
		
			edges = self.findEdges(positions, protein1)
			for protein2 in edges:
				try:
					
					dissimilarity =  math.pow(1 - max(blastData[protein1][protein2],blastData[protein2][protein1]),float(options['normalisation_value']))
				except:
					dissimilarity = 1

				if protein2 in min_queue and dissimilarity < highScoreList[protein2]:
					minimumSpanningTree[positions[protein2]] = protein1
					highScoreList[protein2] = dissimilarity
					min_queue = self.sortMinQueue(min_queue,highScoreList)


		return [minimumSpanningTree,highScoreList]
	#--------------------------------------------------------------------------------------------------------------------#
	def scorePrimBlast(self,blastData,proteinList):
		positions = {}
		count = 0
		for proteins in proteinList:
			positions[proteins] = count
			count += 1
			
		[minimumSpanningTree,key]  = self.prim(positions, blastData, proteinList)
		return 1 + sum(key.values())
	#--------------------------------------------------------------------------------------------------------------------#
	
	def scorePrimBinary(self,binaryDataIn,proteinList):
		positions = {}
		count = 0

		binaryData = copy.deepcopy(binaryDataIn)
		
		for key1 in binaryDataIn:
			for key2 in binaryDataIn:
				if key1 == key2:
					binaryData[key1][key2] = 1.0
				if binaryData[key1][key2] == 0:
					del binaryData[key1][key2]
				if binaryDataIn[key1][key2] == 1 or binaryDataIn[key2][key1] == 1:
					binaryData[key1][key2] = 1
					binaryData[key2][key1] = 1
					
		for proteins in proteinList:
			positions[proteins] = count
			count += 1
			
		[minimumSpanningTree,key]  = self.prim(positions, binaryData, proteinList)

		return 1 + sum(key.values())
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#

  
#----------------------------------------------------------------------------------------------------------------------------#
class rje_Score:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  rje_Score

	  Author:    
	  Rich Edwards

	  Created:  

	  Description:  
	  --
	  
	  Requirements:

	  Fixes:

	  Functions:   

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass

	def expect(self,pattern,aafreq,aanum,seqnum):    	
		### Returns the expected number of occurrences for that pattern
		'''
		Returns the expected number of occurrences for given pattern. Xs and .s both count as wildcards.
		>> aafreq:dictionary of AA frequencies {aa:freq}
		>> aanum:int = sum total of positions in dataset
		>> seqnum:int = number of different sequence fragments searched
		<< expdict:dictionary of {mismatch:expectation}
		'''
		### Setup ###
		terminal_constraint = False
		expvar = rje.replace(pattern,')',') ')   # String to analyses
		expvar = rje.replace(expvar,'[A-Z]','X')
		expvar = rje.replace(expvar,'.','X')
		expvar = rje.replace(expvar,']','] ')
		if expvar[0] == '^' or expvar[-1] == '$':
			terminal_constraint = True
		patlen = 0          # Length of pattern for calculating possible number of positions
		aafreq['X'] = 1.0   # Wildcards do not affect expectation (prob=1)
		prob_per_site = 1.0 # Probability of the pattern occurring at any given site.
		
		### Calculate prob_per_site ###    
		while expvar:       # Still some pattern to look at
		## Deal with spacers inserted for ease of pattern matching ##
			if expvar[0] in [' ','^','$']:
				expvar = expvar[1:]
				continue
		## Wildcard ##
			if expvar[:1] == 'X':   
				expvar = expvar[1:]
				patlen += 1
				continue
		## Update probability per site ##
			if expvar[0] == '[' and expvar.find(']') > 0:   # Choices
				csum = 0.0      # Summed frequency over choices
				for c in expvar[1:expvar.find(']')]:    # Region between []
					csum += aafreq[c]   # Prob = sum of all choices
				if csum < 1.0:
					prob_per_site *= csum
				expvar = expvar[expvar.find(']')+1:]
				patlen += 1
			elif expvar[0] == '(' and expvar.find(')') > 0: # Complex choice!
		    		csum = 0.0      # Sum prob of whole choice
		    		cvar = rje.split(expvar[1:expvar.find(')')],'|')     # List of different options
		    		for cv in cvar:
					cvexp = 1.0 # Probability of just one portion of choice
					while rje.matchExp('(\[(\S+)\])',cv):
			   		 	msum = 0.0  # sum for choice within portion!
			   			cvm = rje.matchExp('(\[(\S+)\])',cv)
			   			cv = rje.replace(cv,cvm[0],'',maxsplit=1)
			    			for m in cvm[1]:
							msum += aafreq[m]   # Prob = sum of all choices
			    			cvexp *= msum
					cv = rje.replace(cv,' ','')
					for c in cv:
			    			cvexp *= aafreq[c]
					csum += cvexp
		   		if csum < 1.0:
					prob_per_site *= csum
		    		expvar = expvar[expvar.find(')')+1:]
		   		patlen += float(len(rje.join(cvar,'')))/len(cvar)    # Add mean length of options
			elif expvar[0] not in aafreq.keys():    # Problem!
		    		print '! aafreq missing <%s> for pattern %s!' % (expvar[0],pattern)
		    		expvar = expvar[1:]
			else:   # Simple
		    		prob_per_site *= aafreq[expvar[0]]
		    		expvar = expvar[1:]
	
	   	 ### Convert to Expectation ###
		if terminal_constraint:
			num_sites = seqnum
		else:
			num_sites = aanum - (seqnum * (patlen - 1))
			
		return prob_per_site * num_sites
	#----------------------------------------------------------------------------------------------------------------------------#
	def occProb(self,observed,expected):     ### Returns the poisson probability of observed+ occurrences, given expected
		'''Returns the poisson probability of observed+ occurrences, given expected.'''
		prob = 0
		for x in range(0,observed):
			try:        #!# Fudge for OverflowError: long int too large to convert to float
				prob += (math.exp(-expected) * pow(expected,x) / self.factorial(x))
			except:
				break
		return 1 - prob
	    
	#----------------------------------------------------------------------------------------------------------------------------#
	def rje_Probability(self,observed,pattern,aminoAcidOcc,proteinData):
		expected = 0
		for protein in proteinData:
			aanum = len(proteinData[protein])
			seqnum = 1
			expectedProtein = self.expect(pattern,aminoAcidOcc,aanum,seqnum)
			expected += expectedProtein
	
		probScore =  self.occProb(observed,expected)
		probabilityScores = [expected,probScore]
		return probabilityScores
	#----------------------------------------------------------------------------------------------------------------------------#
	def factorial(self,m): ### Returns the factorial of the number m
		'''Returns the factorial of the number m.'''
		value = 1
		if m != 0:
			while m != 1:
				value = value*m
				m = m - 1
		return value
	#----------------------------------------------END OF CLASS----------------------------------------------------------#
	


#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Scorer:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Scorer

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   checkScore(scoreDict,score,data,options)
	   @param scoreDict
	   @param score
	   @param data
	   @param options

	   printGlobalMatrix(proteinData,blastHits,blastGlobalHits)
	   @param proteinData
	   @param blastHits
	   @param blastGlobalHits
	   printMatrix(proteinData,blastHits)
	   @param proteinData
	   @param blastHits

	   normalise_Occurances(pattern,teiresiasData,blastGlobalHits)
	   @param pattern
	   @param teiresiasData
	   @param blastGlobalHits

	   scorePatterns(teiresiasOutputData,topRanking,proteinData,options,blastData)
	   @param teiresiasOutputData
	   @param topRanking
	   @param proteinData
	   @param options
	   @param blastData

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################RANKING PATTERNS########################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		self.setList = {}
		self.distribution = {}

	#--------------------------------------------------------------------------------------------------------------------#
	def checkScore(self,scoreDict,score,data,fileHandle,options):
		if (int(score)/5)*5 in self.distribution:
			self.distribution[(int(score)/5)*5] += 1
		else:
			self.distribution[(int(score)/5)*5] = 1

		try:
			fileHandle.write(str(data[2]) + "\t" + str(data[3]) + "\n")
			
		except:
			pass

		
		if score > min(scoreDict.keys()):
			if len(scoreDict) > options["no_of_outputs"] - 1:
				del scoreDict[min(scoreDict.keys())]
			if score in scoreDict:
				scoreDict[score + float(random.randint(1,1000))/100000] = data
			else:
				scoreDict[score] = data

		return scoreDict
	#--------------------------------------------------------------------------------------------------------------------#
	def printGlobalMatrix(self,proteinData,blastHits,blastGlobalHits):
		outString = "\n"

		for protein1 in proteinData:
			outString += "%-10s"%protein1[0:10] + "\t"
			for protein2 in proteinData:
				if protein2 in blastHits[protein1].keys():
					outString += "  %2.2f"%float(blastGlobalHits[protein1][protein2]) + "\t"
				else:
					try:
						outString += " *%2.2f"%float(blastGlobalHits[protein1][protein2]) + "\t"
					except:
						outString += "  0.0\t" 
			outString += "\t\n"


		sys.stderr.write("\n" + outString)

		sys.stderr.write("\n")
	
	#--------------------------------------------------------------------------------------------------------------------#
	def printMatrix(self,proteinData,blastHits):

		outString = "\n"

		for protein1 in proteinData:
			outString += "%-10s"%protein1[0:10] + "\t"
			for protein2 in proteinData:
				if protein2 in blastHits[protein1].keys():
					outString += "   X\t"
				else:
					outString += "   -\t"
			outString += "\n"

		sys.stderr.write("\n" + outString)

	#--------------------------------------------------------------------------------------------------------------------#
	def normalise_Occurances(self,pattern,teiresiasData,blastGlobalHits,options):
		normalised_occurrance = 0

		primAlgorithm = Prim()
		#print blastGlobalHits
		#print teiresiasData
		normalised_occurrance = primAlgorithm.scorePrimBlast(blastGlobalHits,teiresiasData)
		del primAlgorithm	
		
		return normalised_occurrance
		
	#--------------------------------------------------------------------------------------------------------------------#
	def clusterBinary(self,clusterDict,tempDictIn):	
		tempDict = copy.deepcopy(tempDictIn)
		
		cluster = []
		count = 0
		
		for p1 in clusterDict.values():
			cluster.append(sets.Set(p1))
			clusterDict[p1[0]] = count
			count += 1
			
		for p1 in tempDict.keys():
			for p2 in tempDict.keys():
				if p1 != p2:
					if tempDict[p1][p2] == 1:
						temp = []
						temp = cluster[clusterDict[p1]].union(cluster[clusterDict[p2]])
						
						if temp not in cluster:
							cluster.append(temp)
							for protein in temp:
								clusterDict[protein] = len(cluster) - 1
 				else:
 					pass
 				
		countList = []

		for i in clusterDict.values():
			if i not in countList:
				countList.append(i)
		
		return len(countList)
	#--------------------------------------------------------------------------------------------------------------------#
		
	def normalise_Occurances_Nondistance(self,pattern,teiresiasData,blastGlobalAlignments,options):
		pattern1 = re.compile('\[[^\]]*\]|.')

		patternParts = pattern1.findall(pattern)
		
		tempDomainDict = {}
		clusterDomainDict = {}
		tempHomologyDict = {}
		clusterHomologyDict = {}
		offsetDict = {}

		querySet = sets.Set(teiresiasData[pattern][0])

		for i in range(len(teiresiasData[pattern][0])):
			offsetDict[teiresiasData[pattern][0][i]] = teiresiasData[pattern][1][i]

		count = 0
			
		if options['self_hit'] == "T":
			for proteins1 in offsetDict.keys():
				for x in range(len(offsetDict[proteins1])):
					tempDomainDict[proteins1  + str(x)] = {}
					tempHomologyDict[proteins1  + str(x)] = {}

					clusterDomainDict[proteins1 + str(x)] = [proteins1 + str(x)] 
					clusterHomologyDict[proteins1 + str(x)] = [proteins1 + str(x)] 
					

					for proteins2 in offsetDict.keys():
						for y in range(len(offsetDict[proteins2])):
							
							try:
								tempHomologyDict[proteins1 + str(x)][proteins2 + str(y)] = 1
								count = 0
								match = 1

								matched = 0
								mismatched = 0

								for values in patternParts:
									if len(values) == 1:
										if values == blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] or values == "."  or blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] == "+":
											matched += 1
											pass
										else:
											mismatched += 1	
							
									else:
										if blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] in values or blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] == "+":
											matched += 1
											pass
										else:
											mismatched += 1
										
									count += 1
									
								if float(mismatched)/matched > options["match/mismatch"]:	
									tempDomainDict[proteins1 + str(x)][proteins2 + str(y)] = 0
								else:
									tempDomainDict[proteins1 + str(x)][proteins2 + str(y)] = 1

								if proteins1 == proteins2:
									tempDomainDict[proteins1 + str(x)][proteins2 + str(y)] = 0
									tempHomologyDict[proteins1 + str(x)][proteins2 + str(y)] = 0
								
							except:
								tempDomainDict[proteins1 + str(x)][proteins2 + str(y)] = 0
								tempHomologyDict[proteins1 + str(x)][proteins2 + str(y)] = 0
				

		elif options['self_hit'] == "F":
			for proteins1 in offsetDict.keys():
				for x in range(len(offsetDict[proteins1])):
					tempDomainDict[proteins1 ] = {}
					tempHomologyDict[proteins1 ] = {}

					clusterDomainDict[proteins1] = [proteins1] 
					clusterHomologyDict[proteins1] = [proteins1] 
					
					for proteins2 in offsetDict.keys():
						for y in range(len(offsetDict[proteins2])):
							
							
							if proteins2 in blastGlobalAlignments[proteins1].keys(): 
								tempHomologyDict[proteins1][proteins2] = 1
								count = 0
								match = 1

								matched = 0
								mismatched = 0
								for values in patternParts:
									if  values != ".":
										if len(values) == 1:
											if  blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] != "+" or blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] != ".":
												if values == blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count]:
														matched += 1
														pass
												else:														
													mismatched += 1
											else:
												mismatched += 1	
								
										else:
											if blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] in values or blastGlobalAlignments[proteins1][proteins2][(int(offsetDict[proteins1][x])) + count] == "+":
												matched += 1
												pass
											else:
												mismatched += 1
											
									count += 1
								
								
								if matched > 0:
									if float(mismatched)/matched + 0.00000001 >  options["match/mismatch"]:	
										tempDomainDict[proteins1][proteins2] = 0
										break
									else:
										tempDomainDict[proteins1][proteins2] = 1
								else:
									tempDomainDict[proteins1][proteins2] = 0
									
								if proteins1 == proteins2:
									
									tempDomainDict[proteins1][proteins2] = 0
									tempHomologyDict[proteins1][proteins2] = 0
								
							else:
								tempDomainDict[proteins1][proteins2] = 0
								tempHomologyDict[proteins1][proteins2] = 0

			
		primAlgorithm = Prim()

		domainSupport = self.clusterBinary(clusterDomainDict,tempDomainDict)
		homologySupport = self.clusterBinary(clusterHomologyDict,tempHomologyDict)

		'''if 1:
			print '-h'*30
			for key1 in tempHomologyDict:
				print key1,'\t',
				for key2 in tempHomologyDict[key1]:
					print tempHomologyDict[key1][key2],
				print 
				
			print float(homologySupport)
			print '-'*30'''
		return [domainSupport,homologySupport]
	#--------------------------------------------------------------------------------------------------------------------#

	def normalise_Occurances_SlidingWindow_MST(self,pattern,teiresiasData,blastGlobalAlignments,options):
		pattern1 = re.compile('\[[^\]]*\]|.')
		
		patternParts = pattern1.findall(pattern)
		
		tempDict = {}
		offsetDict = {}
		clusterDict = {}
		
		querySet = sets.Set(teiresiasData[pattern][0])
		for i in range(len(teiresiasData[pattern][0])):
			offsetDict[teiresiasData[pattern][0][i]] = teiresiasData[pattern][1][i]

		print offsetDict
		windowSize = 30
		for proteins1 in offsetDict.keys():
			
			tempDict[proteins1] = {}
			clusterDict[proteins1] = [proteins1] 
			for proteins2 in offsetDict.keys():
				count = 0
				for offset in offsetDict[proteins1]:
					count += 1
					print proteins1 + '_' + str(count),'\t', proteins2,'\t', 
					try:
						if proteins2 in blastGlobalAlignments[proteins1].keys():
							start = offset - windowSize/2
							end = offset + len(pattern) + windowSize/2
							
							print blastGlobalAlignments[proteins1][proteins2][start:end],
							mismatch =  blastGlobalAlignments[proteins1][proteins2][start:start+windowSize/2].count('.') +  blastGlobalAlignments[proteins1][proteins2][end:end+windowSize/2].count('.')
							print start,end,
							print blastGlobalAlignments[proteins1][proteins2][start:start+windowSize/2],
							print blastGlobalAlignments[proteins1][proteins2][end:end+windowSize/2],
							print blastGlobalAlignments[proteins1][proteins2][start:start+windowSize/2].count('.'),
							print blastGlobalAlignments[proteins1][proteins2][end:end+windowSize/2].count('.'),
							#print len(blastGlobalAlignments[proteins1][proteins2][start:start+windowSize/2] + blastGlobalAlignments[proteins1][proteins2][end:end+windowSize/2])
							print "%-2.2f"%(float(mismatch)/len(blastGlobalAlignments[proteins1][proteins2][start:start+windowSize/2] + blastGlobalAlignments[proteins1][proteins2][end:end+windowSize/2]))
						else:
							print "No hit"
							
					except:
						print "failed"
					
	#--------------------------------------------------------------------------------------------------------------------#
	def scorePatterns(self,teiresiasOutputData,topRanking,proteinData,options,blastData):
		printLevel(1,options,"Ranking Patterns")
		
		[blastHits,blastGlobalHits,blastGlobalAlignments] = blastData

		tempString = ''
		for protein1 in blastGlobalAlignments:
			for protein2 in blastGlobalAlignments[protein1]:
				tempString += protein1 + "\t" + protein2 + "\t" + "%-3.2f"%(blastGlobalHits[protein1][protein2]*100) + "\t" + blastGlobalAlignments[protein1][protein2] + "\n"
		
		open(options["output"].replace(".rank", ".homology"),"w").write(tempString)
		
		[percentage_residues_accessible,count] = [0,0]
		
		surfaceProbabilityCheck = SurfaceProbability()
		scorer = Score()
		aminoAcidOcc = scorer.loadaminoAcidOccurances(options["aminoAcidOccuranceFile"])
		del scorer
	
		tempFile = open(options["output"].replace(".rank", ".surface.html"),"w")
		tempFile.write("<html>")
		tempFile.write("<h2>Surface Accessibility<br></h2><h3>Cut-off : \t" + str(options["surface_cutOff"]) + "</h3><br>") 
		
		for protein in proteinData.keys():
			tempFile.write("<h3>" + protein + "</h3>\n")
			tempFile.write("<p><img src='" +  os.path.abspath("./" + options["input_path"] + "/surface/"  + protein + ".gif'") + "></p><br>\n")
			stringTemp = surfaceProbabilityCheck.checkSequence(proteinData[protein],options,protein)

			tempFile.write("<FONT SIZE=1 color=#ff0000 face='courier'>")
			

			for i in range(0,(len(stringTemp)/100) + 1):
				if 60*(i+1) < len(stringTemp):
					stringOut = stringTemp[100*i:100*(i+1)].replace("X","<FONT color=#000000 >X</FONT>")
				
					tempFile.write(stringOut  + "<br>\n")
				else:				
					stringOut = stringTemp[100*i:100*(i+1)].replace("X","<FONT color=#000000 >X</FONT>")
				
					tempFile.write(stringOut  + "<br>\n")
					
			
			tempFile.write("</FONT>")
		
		tempFile.write("</html>")
		tempFile.close()

		fileHnd = open(options["output"].replace("rank","distribution"),"w")
		fileHnd.write("norm\tIC\n")

		sys.stderr.write(time.strftime('%X') + ":\n" + "Ranking Patterns\n\n")
		
		start =  time.time()
		counter = 0
		countOut = 0

		length = len(teiresiasOutputData.keys())

		if options["view_matrix"] == "T":
			self.printMatrix(proteinData,blastHits)
			self.printGlobalMatrix(proteinData,blastHits,blastGlobalHits)

		rje_Scorer = rje_Score()
		fileReader = fileManipulation()
		teiresiasInputProteins =  fileReader.readFile(options,options["TEIRESIAS_input"])
		del fileReader			
		
		if int(options["flavour"]) == 1 or int(options["flavour"]) == 2 or int(options["flavour"]) == 3:
	
			
			for pattern in teiresiasOutputData.keys():
				#self.normalise_Occurances_SlidingWindow_MST(pattern,teiresiasOutputData,blastGlobalAlignments,options)
				#sys.exit()
				"""pos1 = 0
				windowSize = 20
				
				for protein1 in teiresiasOutputData[pattern][0]:
					pos2 = 0
					for protein2 in teiresiasOutputData[pattern][0]:
						print protein1, '-', protein2,
						if protein2 in blastGlobalAlignments[protein1].keys():
							start1 = teiresiasOutputData[pattern][1][pos1][0] - windowSize/2
							start2 = teiresiasOutputData[pattern][1][pos2][0] - windowSize/2
							
							end1 = start1 + len(pattern) + windowSize/2
							end2 = start2 + len(pattern) + windowSize/2
							
							print '?'
							print blastGlobalAlignments[protein1][protein2]
							print blastGlobalAlignments[protein1][protein2][start1:end1],
							print blastGlobalAlignments[protein1][protein2][start2:end2]
						else:
							print 0
						pos2 += 1
					pos1 += 1
					
				sys.exit()"""
				
				if countOut%1000 == 0:
					if options['verbosity'] > 1:
						print "%-10s"%str(length) + "\t%-10s"%str(counter) + "\r",

				counter += 1
				similarity_Normalised_Occurances =  self.normalise_Occurances(pattern,teiresiasOutputData[pattern][0],blastGlobalHits,options)
				[domain_Normalised_Occurances,homology_Normalised_Occurances] = self.normalise_Occurances_Nondistance(pattern,teiresiasOutputData,blastGlobalAlignments,options)
				
				normalised_Supports = [similarity_Normalised_Occurances,domain_Normalised_Occurances,homology_Normalised_Occurances]
				
				if int(options["flavour"]) == 1:
					normalised_occurances = similarity_Normalised_Occurances 
				elif int(options["flavour"]) == 2:
					normalised_occurances = domain_Normalised_Occurances 
				elif int(options["flavour"]) == 3:
					normalised_occurances = homology_Normalised_Occurances
				else:
					print "Error occured with normalised_occurances"
					print int(options["flavour"])

					
				if homology_Normalised_Occurances > 1:
					countOut += 1
					observed = len(teiresiasOutputData[pattern][0])
					
					probabilityScores = rje_Scorer.rje_Probability(observed,pattern,aminoAcidOcc,teiresiasInputProteins)
					[expected,probScore] = probabilityScores
					
					information_content = teiresiasOutputData[pattern][2]
					
					if options["use_evolution"] == "T":	
						
						if int(options["ranking"]) == 1:
							score = information_content*(normalised_occurances)
						elif int(options["ranking"]) == 2:
							score = information_content*(normalised_occurances)*(observed/expected)
						elif int(options["ranking"]) == 10:
							score = information_content*observed
						else:
							errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nOption ranking = " + str(options["ranking"]) + " does not exist.\nPlease chosse {1-2] try again.\n" + "+"*70
							print errorfile_not_found
							sys.exit()
					else:
						
						normalised_occurances = len(teiresiasOutputData[pattern][0])
						score = information_content*normalised_occurances

					offsetDict = {}
					
					for i in range(len(teiresiasOutputData[pattern][0])):
						offsetDict[teiresiasOutputData[pattern][0][i]] = teiresiasOutputData[pattern][1][i]
										
					patternInfo = [pattern,teiresiasOutputData[pattern],normalised_occurances,information_content,percentage_residues_accessible,normalised_Supports,probabilityScores]#,correlation]
					topRanking = self.checkScore(topRanking,score,patternInfo,fileHnd,options)

				del teiresiasOutputData[pattern]		
		else:
			errorfile_not_found = "\nERROR:\n" + "+"*70 + "\nOption flavour = " + str(options["flavour"]) + " does not exist.\nPlease chosse {1-3] try again.\n" + "+"*70
			print errorfile_not_found

			sys.exit()

		printLevel(1,options,"%-10s"%str(length) + "\t%-10s"%str(counter))
				
		printLevel(1,options,"Finished Scoring Clusters")

		outString = "\n" 
		outString +=  "Patterns                      \t:" + str(counter) + "\n"
		outString +=  "Patterns in unrelated proteins\t:" + str(countOut) + "\n"
		
		try:
			outString +=  "Percentage processed          \t:" "%2.2f"%(float(countOut)/counter) + "%\n"
		except:
			pass
			
		outString +=  "Time			     \t:" + "%2.2f seconds"%float(time.time() - start) + "\n"
		
		printLevel(1,options,outString)
		sys.stderr.write("Search Statistics" + outString  + "\n")
		sys.stderr.write(time.strftime('%X') + ":\n" + "Finished Ranking Patterns\n")
		
		try:
			sortList = self.distribution.keys()
			numPatterns = sum(self.distribution.values())
			#print numPatterns
			sortList.sort()

			printLevel(1,options, "------------------------------------------------------------------")
			printLevel(1,options, "Score Distribution")
			for i in range(min(sortList),max(sortList) + 5,5):
				try:
					printLevel(1,options, i ,"\t",self.distribution[i],"\t","%-2.5f"%(float(self.distribution[i])/numPatterns),(self.distribution[i]/(numPatterns/100)+1)*"*")
				except:
					printLevel(1,options, 0,"\t","0.00000")
					pass
		except:
			pass

	
		return topRanking
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class FileWriter:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  FileWriter
	  
	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   printOutput(topRanking,options)
	   @param topRanking
	   @param options

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		sys.stderr.write("\n----------------------------------------------------------------\n")
		sys.stderr.write("########################CREATING OUTPUT#########################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		printLevel(1,options,"-----------------------------" + time.strftime('%X') + "-----------------------------")
		pass

	#--------------------------------------------------------------------------------------------------------------------#
	def htmlOutput(self,topRanking,	teiresiasData,options):
		'''HTML Output method.'''
		printLevel(1,options,os.path.abspath(options["filename"]))
		printLevel(1,options, os.path.abspath(options["aminoAcidOccuranceFile"]))

		if options['memory_saver'] == 'F':
			header = "<h2 align='left'>Data</h2><table width='900' border='1' cellspacing='0' cellpadding='5'>" 
			header += "<tr><td align='left'><a href='" + os.path.abspath(options["filename"]) + "'>Input File</a></td>"
			header += "<td align='left'><a href='" + os.path.abspath(options["output"]) + "'>Tab Delimited Output</a></td></tr>"
			header += "<tr><td align='left'><a href='" + os.path.abspath(options["logfile"])+ "'>Logfile</a></td>"
			header += "<td align='left'><a href='" + os.path.abspath(options["aminoAcidOccuranceFile"]) + "'>Amino Acid Occurance Probabilities</a></td></tr>"
			header += "<tr><td align='left'><a href='" +	os.path.abspath(options["TEIRESIAS_output"]) + "'>TEIRESIAS output</a></td>"
			header += "<td align='left'><a href='" + os.path.abspath(options["output"].replace(".rank",".alignments.html"))+ "'>Global BLAST Generated Alignments</a></td></tr>"
			header += "<tr><td align='left'><a href='" + os.path.abspath(options["output"].replace(".rank",".surface.html")) + "'>Emini Generate Surface Maps</a></td>"
			if options["long_format"] == "T":
				header += "<td align='left'><a href='" + os.path.abspath(options["output"].replace(".rank",".domains.html")) + "'>Domain View</a></td></tr>"
			header += "</table><br>"
		else:
			header = ""
		
		tableHeader = """<h2 align='left'>Results</h2>
		<table width="900" border="1" cellspacing="0" cellpadding="5">
		<align="left">
		<tr>
		<td align="left">Rank</td>
		<td align="left">Pattern</td>
		<td align="left">Score</td>
		<td align="left">Occ</td>
		<td align="left">Proteins</td>
		<td align="left">IC</td>
		<td align="left">N Occ</td>		
		<td align="left">Sim</td>
		</tr>
		</align>"""


		tableFooter = """<br><b>IC = Information Content<br>
		Occ = Occurances<br>
		N Occ = Normalised Occurances<br>
		Sim = Similarity<br>
		</html>"""

		row = """<tr>
		<td align="center"><b>RANK</b></td>
		<td align="left">MOTIF</td>
		<td align="left">SCORE</td>
		<td align="left">OCCURANCES</td>
		<td align="left">PROTEINS</td>
		<td align="left">IC</td>
		<td align="left">NORMALISED</td>
		<td align="left">SIMILARITY</td>

		</tr>
		"""

		html = open(options["filename"]+ ".output.html","w")
		html.write("<html><h1 align='center'>" + options["filename"] + " Results</h1>")
		
		html.write(header)
		html.write(tableHeader)
		temp = ""

		sort = topRanking.keys()
		sort.sort()

		sort.reverse()
		counter = 1
		
#		0 pattern
#		1 proteins
#		2 normalised_occurances
#		3 information_content
#		4 minExpected
#		5 maxExpected
#		6 percentage_residues_accessible
#		7 correlation
		##############################print dat.rank####################################
		proteinDict = {}
		motifDict = {}
		motifList = []
		
		if 0 in sort:
			sort.remove(0)
			
		for scores in sort:
			motifList.append(topRanking[scores][0])
			motifDict[topRanking[scores][0]] = {}
			for protein in topRanking[scores][1][0]:
				if protein in proteinDict:
					pass
				else:
					proteinDict[protein] = len(proteinDict)
					
			for i in range(0,len(topRanking[scores][1][0])):
				motifDict[topRanking[scores][0]][topRanking[scores][1][0][i]] = topRanking[scores][1][1][i]
		
		compressedOutput = ''
		
		compressedOutput += '#Proteins\n'
		for i in range(0,len(proteinDict)):
			for protein in proteinDict:
				if i == proteinDict[protein]:
					compressedOutput += str(proteinDict[protein]) + '\t' + protein  +'\n'
		
		compressedOutput += '#Motif\n'
		count = 0
		for motif in motifList:
			count += 1
			compressedOutput += str(count) + '\t' + motif + '\t'
			tempDict = {}
			for protein in motifDict[motif]:
				tempDict[proteinDict[protein]] = []
				for offset in motifDict[motif][protein]:
					tempDict[proteinDict[protein]].append(str(offset))
					#compressedOutput += str(proteinDict[protein]) + ':' + str(offset)  + ' '
			
			sortTemp = tempDict.keys()
			sortTemp.sort()
			
			for key in sortTemp:
				for offset in tempDict[key]:
					compressedOutput += str(key) + ':' + str(int(offset)  + 1)  + ' '

			compressedOutput += '\n'
		

		open(options['output'].replace('rank','dat.rank'),'w').write(compressedOutput)
		

		for scores in sort:
			#print scores
			try:	
				#print topRanking[scores][1][0]
				proteinString = "<TABLE border='0' table width='100%' cellspacing='0' cellpadding='2'>"
				for i in range(0,len(topRanking[scores][1][0])):
					offsetsTemp = topRanking[scores][1][1][i]
					offsetsTemp.sort()
					offsetString = ""

					for offset in offsetsTemp :
						
						offsetString += str(int(str(offset).replace("'",""))) + "-" +  str(int(str(offset).replace("'","")) + len(topRanking[scores][0]))  + " , "
					
					offsetString = offsetString
					proteinString += "<tr><td width='50%'><a href='http://us.expasy.org/cgi-bin/niceprot.pl?" + str(topRanking[scores][1][0][i]) + "'>" + str(topRanking[scores][1][0][i]) +  "</a></td><td>" + offsetString[:-2] + "<br></td></tr>"
				proteinString += "</TABLE>"

				temp = row.replace("RANK",str(counter))
				temp = temp.replace("MOTIF",str(topRanking[scores][0]))
				temp = temp.replace("SCORE","%-2.3f"%(scores))
				temp = temp.replace("OCCURANCES",str(len(topRanking[scores][1][0])))
				temp = temp.replace("PROTEINS",proteinString)
				temp = temp.replace("IC","%-2.3f"%(topRanking[scores][3]))
				temp = temp.replace("NORMALISED",str(topRanking[scores][2]))
				temp = temp.replace("SIMILARITY",str(topRanking[scores][2]/len(topRanking[scores][1][0])))
				html.write(temp)
				counter += 1
			except Exception,e:
				print e
				pass
		html.write("</table><br>")
		html.write(tableFooter)
		html.close()

	#--------------------------------------------------------------------------------------------------------------------#
	def printOutput(self,topRanking,stats,options):
		printLevel(1,options,"Formatting output")
		
		 
		outHnd = open(options["output"],"w")
		sort = topRanking.keys()
		sort.sort()
		sort.reverse()
		
		counter = 0
		printLevel(1,options,'')
		inputStatsString =  "------------------------------------------------------------------\n"

		sys.stderr.write(time.strftime('%X') + ":\n" + "Scores\n")
		
		outstring = "Rank\tScore\t%-20s"%"Pattern" + "\tOcc\t"  + "IC\t"+ "Norm\t" + "Sim\t" +"\tMST\tUHS\tUP\t\tExpect\tProb\n"
		outHnd.write( '#' + (stats['ProteinData'] + '\n' +'No. of motifs: '  + str(stats['NumMotifs']) ).replace('\n','\n#\t') + '\n' +"#------------------------------------------------------------------\n" + outstring)
		
		
		sys.stderr.write(outstring)
		printLevel(1,options,outstring)
		
		for scores in sort:	
			try:
				counter += 1
				
				if scores > 1000:
					eTemp = ("%-1.2e"%float(scores)).split("+")
					outstring = "(" + str(counter) + ")\t" + eTemp[0] + "+" + str(int(eTemp[1]))
				else:
					outstring = "(" + str(counter) + ")\t%2.3f"%scores 

				if len(str(topRanking[scores][0])) > 35:
					outstring += "SEQUENCE " 
				else:
					outstring += "SEQUENCE "

				outstring += "%-2s\t"%str(len(topRanking[scores][1][0])) 
				outstring += "%-2.3f\t"%float(topRanking[scores][3])
				outstring += "%-2.3f\t"%float(topRanking[scores][2])
				outstring += "%-2.3f\t"%(float(topRanking[scores][2])/len(topRanking[scores][1][0]))
				outstring += "\t"
				outstring += "%-2.3f\t"%float(topRanking[scores][-2][0])
				outstring += "%-2.3f\t"%float(topRanking[scores][-2][1])
				outstring += "%-2.3f\t"%float(topRanking[scores][-2][2])
				outstring += "\t"		
				outstring += "%-2.3f\t"%float(topRanking[scores][-1][0])
				outstring += "%-2.2g\t"%float(topRanking[scores][-1][1])
				
				outstring +="\n"
		
					
				outHnd.write(outrje.replace("SEQUENCE","\t%-20s "%str(topRanking[scores][0])) )
				sys.stderr.write(outrje.replace("SEQUENCE","\t%-20s "%str(topRanking[scores][0])) )
				
				if len(str(topRanking[scores][0])) > 20:
					printLevel(1,options, outrje.replace("SEQUENCE","\t%-20s"%(str(topRanking[scores][0])[0:17])  + "...  ")[:-1])
				else:
					printLevel(1,options, outrje.replace("SEQUENCE","\t%-20s"%str(topRanking[scores][0]))[:-1])
			
			except Exception,e:
				#print e
				#printLevel(1,options, scores)
				pass
				
		outHnd.close()
	#----------------------------------------------END OF CLASS----------------------------------------------------------#



#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class CommandLine:	
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  CommandLine

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   
	   isInt(i) 
	   @param i

	   isFloat(i) 
	   @param i

	   printOptions()

	   usage()

	   parseCommandLine()

	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		
		self.optionsVariables = {
				"filename":"-i",
				"output":"-o",
				"logfile":"-l",
				"batch":"-b",
				"run_TEIRESIAS":"-T",
				"run_formatDB":"-D",
				"run_BLAST":"-B",
				"BLAST_results":"-r",
				"self_hit":"-H",
				"aminoAcidOccuranceFile" : "-a",
				"no_of_outputs": "-n",
				"information_cutOff":"-y",
				"complexity_cutOff":"-z",
				"surface_cutOff":"-s",
				"ranking":'-K',
				"filetype":'-t',
				"view_matrix": "-v",
				"use_evolution":"-E",
				"use_filtering":"-A",	
				"long_format":"-L",
				"query_protein":"-q",
				"gap_weight":"J",
				"eVal" : "-e",
				"flavour":"-f",
				"match/mismatch":"-p",
				'overwrite':'-O',
				"TEIRESIAS_output":"",
				"TEIRESIAS_input":"-I",
				"TEIRESIAS_patternLength":"-P",
				"TEIRESIAS_fixedPositions": "-F",
				"TEIRESIAS_supportString":"-S",
				"TEIRESIAS_equiv":"-R",
				"TEIRESIAS_walltime":"-W",
				'normalisation_value':"-C",
				'TEIRESIAS_local(IBM_path_length_bug_fix)':"-G",
				"Mask_path":"-m",
				"BLAST_path":"",
				"TEIRESIAS_path":"",
				"input_path":"",
				'memory_saver':'-X',
				'ini_file':'-d',
				'verbosity':'-Q',
				'strict_inclusive_masking':'-M'
				}

		self.options = {
				"filename":"",
				"output":"",
				"logfile":"",
				"batch":"F",
				"run_TEIRESIAS":"F",
				"run_formatDB":"T",
				"run_BLAST":"T",
				"self_hit":"F",
				"BLAST_results":"./results/",
				"strict_inclusive_masking":"F",
				"aminoAcidOccuranceFile" : "PeptideOccurances.txt",
				"no_of_outputs": 100,
				"information_cutOff":8,
				"complexity_cutOff":-1,
				"surface_cutOff":0.0,
				"view_matrix": "F",
				"use_evolution":"T",
				"use_filtering":"T",
				"long_format":"F",
				"query_protein":"",
				"gap_weight":0.5,
				"eVal" : 0.001,
				"flavour":"1",
				"filetype":'',
				"ranking":'1',
				'overwrite':'T',
				'TEIRESIAS_local(IBM_path_length_bug_fix)':'F',
				"match/mismatch":0.2,
				"TEIRESIAS_output":"",
				"TEIRESIAS_input":"",
				"TEIRESIAS_patternLength":8,
				"TEIRESIAS_fixedPositions": 2,
				"TEIRESIAS_supportString":3,
				"TEIRESIAS_equiv":"F",
				"TEIRESIAS_walltime":0,
				"BLAST_path":"",
				'verbosity':'0',
				"Mask_path":"mask.dat",
				"TEIRESIAS_path":"",
				'normalisation_value':1,
				'memory_saver':'F',
				'ini_file':'slimdisc.ini',
				"input_path":""
				}
	#--------------------------------------------------------------------------------------------------------------------#
	def isInt(self,i): 
		try: 
			int(i) 
			return 1 
		except: 
			return 
	#--------------------------------------------------------------------------------------------------------------------#
	def isFloat(self,i): 
		try:
			float(i) 
			return 1 
		except: 
			return 
	#--------------------------------------------------------------------------------------------------------------------#
	def printOptions(self):
		outString = ""
		printLevel(1,self.options, "-----------------------------" + time.strftime('%X') + "-----------------------------")
		sys.stderr.write("----------------------------------------------------------------\n")
		sys.stderr.write("########################OPTIONS#################################\n")
		sys.stderr.write("----------------------------------------------------------------\n")
		
		sort = self.options.keys()
		sort.sort()
		for option in sort:
			outString += "%-30s"%option + "\t" + str(self.options[option]) + "\n"

		sys.stderr.write(outString)
		printLevel(1,self.options,outString)
	#--------------------------------------------------------------------------------------------------------------------#
	def usage(self):
		printLevel(1,self.options,"%-30s"%str("Options")  + "\t" + "%-10s"%"Prefix" + "\t" + "%-20s"%str("Defaults"))
		variable = "-a"
		printLevel(1,self.options, "------------------------------------------------------------------------------------------")
		sort = self.options.keys()
		sort.sort()
		for option in sort:
			printLevel(1,self.options,"%-30s"%option  + "\t" + "%-10s"%self.optionsVariables[option]  + "\t" + "%-20s"%str(self.options[option]))
		
		printLevel(1,self.options, "------------------------------------------------------------------------------------------")

 	#--------------------------------------------------------------------------------------------------------------------#
	def loadINI(self,dir):
		try:
			#print os.path.join(dir,self.options['ini_file'])
			
			iniFile = open(os.path.join(dir,self.options['ini_file']),"r").read()
			data = re.sub('\n+','\n',iniFile).split("\n")
			
			for value in data:
				if len(value) > 2:
					iniValues = value.split("=")

					if self.isInt(self.options[iniValues[0].strip()]):
						if iniValues[1].find(".") == -1:
							self.options[iniValues[0].strip()] = int(iniValues[1].strip())
						else:
							self.options[iniValues[0].strip()] = float(iniValues[1].strip())
					else:
						self.options[iniValues[0].strip()] = str(iniValues[1].strip())	
			
		except:
			if os.path.exists(self.options['ini_file']):
				error_iniFileError = '#'*40 + "\nINI file " + self.options['ini_file'] + " corrupt using default settings\n" + '#'*40 
			else:
				error_iniFileError = '#'*40 + "\nINI file " + self.options['ini_file'] + " missing using default setting\nCheck path\n" + '#'*40 
			print error_iniFileError
			time.sleep(2)

 	#--------------------------------------------------------------------------------------------------------------------#
	def parseCommandLine(self):
		if len(sys.argv) == 1:
			print '-'*23 + '\nPlease enter arguements\n' + '-'*23
			self.usage()
			sys.exit(2)
		else:
			if os.path.isfile('slimdisc.ini'):
				self.loadINI('')
			elif os.path.isfile(os.path.join(os.path.dirname(sys.argv[0]),'slimdisc.ini')): 
				self.loadINI(os.path.dirname(sys.argv[0]))

			try:
				sort = self.optionsVariables.values()
				sort.sort()
				opts, args = getopt.getopt(sys.argv[1:],"A:B:C:D:E:F:G:H:I:J:K:L:O:M:N:P:R:S:T:X:Q:W:a:b:c:d:e:f:i:l:m:n:o:p:q:r:s:t:v:w:x:y:z:h", self.options.keys())

			except getopt.GetoptError,e:
				print e
				sys.exit(2)
			
			output = None
		
			for o, a in opts:
				for option in self.optionsVariables.keys():
					if o ==self.optionsVariables[option]:
						if self.isInt(self.options[option]):
							if a.find(".") == -1:
								self.options[option] = int(a)
							else:
								self.options[option] = float(a)
						else:
							self.options[option] = a	
				
				if o in ("-h", "--help"):
					self.usage()
					
					print __doc__
					sys.exit(2)

			if self.options['ini_file'] != 'slimdisc.ini':
				self.loadINI('')
			
			self.options["calling_folder"] = os.path.dirname(sys.argv[0])
		
			return self.options
	#----------------------------------------------END OF CLASS----------------------------------------------------------#
		
		
#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class Expectation:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  Expectation

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to run Teiresias Algorithm

	  Requirements:

	  Fixes:

	  Functions:   


	#################################################################################################
	"""
	#################################################################################################

	def __init__(self,blastData,proteinData,aaOcc):
		self.totalLength = 0
		self.proteinData = proteinData
		self.aaOcc = aaOcc

		for proteins in proteinData.keys():	
			self.totalLength += len(proteinData[proteins])

	#--------------------------------------------------------------------------------------------------------------------#
	def collapseList(self,listTemp):
		stringTemp = ""
		for values in listTemp:
			stringTemp += values

		return stringTemp
	#--------------------------------------------------------------------------------------------------------------------#
	def calculateExpectationRange(self,pattern):
		maxProb = (((self.totalLength/len(self.proteinData.keys()) - (self.patternLength(pattern) + 1)))*len(self.proteinData.keys()))*(self.calculateScore(pattern))
		minProb = maxProb * (len(self.proteinData.keys())/self.normalised_occurrance)
		return [maxProb,minProb]
	#--------------------------------------------------------------------------------------------------------------------#
	def calculateScore(self,pattern):
		aaOcc = {"A":0.066,"C":0.017,"B":0.000,"E":0.062,"D":0.056,"G":0.089,"F":0.030,"I":0.049,"H":0.020,"K":0.049,"M":0.018,"L":0.107,"N":0.039,"Q":0.049,"P":0.056,"S":0.062,"R":0.063,"T":0.060,"W":0.013,"V":0.058,"Y":0.037,"X":0.000,"Z":0.000}
		
		count = 0	
		p1 = re.compile('\[[^\]]*\]')
		p2 = re.compile('(?<=])[^\[]*')

		data = p1.findall(pattern)
		pattern = "]" + pattern
		probability = 1.0
		for values in data:
			count += 1
			sum = 0.0
			for AA in values[1:-1]:
				if sum == 0.0:	
					sum = aaOcc[AA]
				else:
					sum += aaOcc[AA]
			
			probability *= sum

		count += pattern.count(".")
		data = p2.findall(pattern)
		
		stringTmp = self.collapseList(data)
		stringTmp = stringTmp.replace(".","")
		
		for AA in stringTmp:
			count += 1
			xx = -aaOcc[AA]*math.log(aaOcc[AA])
			probability *= aaOcc[AA]
		
		return probability
	#--------------------------------------------------------------------------------------------------------------------#
	def patternLength(self,pattern):
		p = re.compile('(?<=])[^\[]*')
		
		data = p.findall("]" + pattern)
		length = pattern.count("[") + len(self.collapseList(data))
		return length
	#----------------------------------------------END OF CLASS----------------------------------------------------------#
	
#------------------------------------------------------CLASS-----------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
class FileChecker:
	#################################################################################################
	"""
	#################################################################################################

	  Class:
	  FileChecker

	  Author:    
	  Norman Davey 

	  Created:       
	  08/03/2005

	  Description:  
	  Class used to check if all files are present

	  Requirements:

	  Fixes:

	  Functions:   


	#################################################################################################
	"""
	#################################################################################################

	def __init__(self):
		pass

	#--------------------------------------------------------------------------------------------------------------------#
	def checkFiles(self,options):
		printLevel(1,options, "-----------------------------" + time.strftime('%X') + "-----------------------------")
		printLevel(1,options, "Checking presence of Files")
	
		if os.path.exists(options["filename"]):
			printLevel(1,options, "%-49s"%options["filename"] + "\t" + "present")
		else:
			print "#"*50 +'\nInput file path does not exist\n' + options["filename"] + '\n' + "#"*50
			sys.exit()

		if os.path.exists(os.path.join(options["BLAST_path"],'blastall.exe')) or os.path.exists(os.path.join(options["BLAST_path"],'blastall')):
			printLevel(1,options,"%-49s"%"blastall" + "\t" + "present")
		else:
			print "#"*20 +"\nERROR: blastall not present at '" + options["BLAST_path"] + "' .Please check path\n"+ "#"*97  +"\n" 
			sys.exit()
		
		if os.path.exists(os.path.join(options["TEIRESIAS_path"],"teiresias_char.exe")) or os.path.exists(os.path.join(options["TEIRESIAS_path"],"teiresias_char")):
			printLevel(1,options, "%-49s"%"TEIRESIAS" + "\t" + "present")
		else:
			print "#"*97 +"\nERROR: TEIRESIAS not present at " + options["TEIRESIAS_path"]  +". Please check path\n"+ "#"*97  +"\n" 
			sys.exit()
			
	#----------------------------------------------END OF CLASS------------------


#------------------------------------------------------MAIN------------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
def printLevel(level,options,string):
	if int(options['verbosity']) < int(level):
		pass
	else:
		print string

def runMain(options):

	topRanking = {0:[]}	
	offsets = {}
	stats = {}
	
	x = 1

	if options["long_format"] == "T":
		fileObject = fileManipulationFull()
		offsets,proteinData = fileObject.readFile(options)
		if len(proteinData.keys()) < 3:
			print "Dataset contains only " + str(len(proteinData.keys())) + " proteins. Too small for analysis"
			sys.exit()
			
		aaOcc = aaOccurance()
		aaOcc.createAAOccuranceFilefromProteinList(options,proteinData)

	else:
		fileObject = fileManipulation()
		proteinData = fileObject.readFile(options,options["filename"])
		if len(proteinData.keys()) < 3:
			print "Dataset contains only " + str(len(proteinData.keys())) + " proteins. Too small for analysis"
			sys.exit()
			
		aaOcc = aaOccurance()
		aaOcc.createAAOccuranceFilefromProteinList(options,proteinData)
		fileObject.createTEIRESIASinputfromProteinList(proteinData,options)

	
	if float(options["TEIRESIAS_supportString"]) < 1:
		if int(len(proteinData)*float(options["TEIRESIAS_supportString"])) > 3:
			options["TEIRESIAS_supportString"] = int(len(proteinData)*float(options["TEIRESIAS_supportString"]))
		else:
			options["TEIRESIAS_supportString"] = 3

	
	if options["run_TEIRESIAS"] == "T":
		teiresiasRunner = Teiresias()
		teiresiasRunner.teiresias(options)
		del teiresiasRunner

	filereader = Filereader()	
	teiresiasData = filereader.readTeiresiasData(offsets,options)

	del filereader

	run = RunBLAST()
	blastOutput = run.runBLAST(proteinData,options)
	#print blastOutput
	del run

	primAlgorithm = Prim()
	normalised_occurrance = primAlgorithm.scorePrimBlast(blastOutput[1],proteinData.keys())

	stats['overall_MST'] = normalised_occurrance
	stats['NumMotifs'] = len(teiresiasData)
	
	inputStatsString =  "---------------------------Input stats----------------------------\n"
	inputStatsString += 'Input dataset:   ' + options['filename'] + '\n'
	inputStatsString += 'No. of proteins: ' + str(len(proteinData)) + '\n'
	inputStatsString += 'Overall MST:     ' + "%-2.5f"%stats['overall_MST'] + '\n\n'
	
	for protein in proteinData:
		inputStatsString += '-%-20s'%protein + '\t' + '%5d'%(len(proteinData[protein])) + ' a.a' + '\n'

	inputStatsString += '%-20s'%' ' + '\t' + '-----' + '\n'
	inputStatsString += '%-20s'%'No. of residues: ' + '\t' + '%5d'%(len(''.join(proteinData.values()))) + '\n'
	
	stats['ProteinData'] = inputStatsString
	
	printLevel(1,options, inputStatsString)
	sys.stderr.write(inputStatsString)
	
	stats['ProteinData'] = inputStatsString
	
	scorer = Scorer()
	
	topRanking = scorer.scorePatterns(teiresiasData,topRanking,proteinData,options,blastOutput)
	del scorer
	
	fileWriter = FileWriter()
	fileWriter.htmlOutput(topRanking,teiresiasData,options)
		
	fileWriter.printOutput(topRanking,stats,options)
	del fileWriter

	if options['memory_saver'] == 'F':
		try:
			distPlot = DistributionPlotter(options)
			distPlot.plotDistribution(options,10)
			del distPlot
		except:
			pass
			#printLevel(1,options, 'Skipping distribution plotting')
	else:
		pass
		#printLevel(1,options, 'Skipping distribution plotting')

	sys.stderr.write("----------------------------------------------------------------\n")
	printLevel(1,options, "-----------------------------" + time.strftime('%X') + "-----------------------------")


def setupFiles(options):
	
		#print "-----------------------------" + time.strftime('%X') + "-----------------------------"
		
		fileChecker = FileChecker()
		fileChecker.checkFiles(options)
		
		options["filename"] = os.path.abspath(options["filename"])
	
		if os.path.basename(options["filename"]).split(".")[0] in os.listdir(os.path.dirname(options["filename"])):
			pass
		else:
			printLevel(1,options,os.path.abspath(os.path.dirname(options["filename"]) + "/" + os.path.split(options["filename"])[-1].split(".")[0]))
			os.mkdir(os.path.abspath(os.path.dirname(options["filename"]) + "/" + os.path.split(options["filename"])[-1].split(".")[0]))
		
		
		path = os.path.dirname(options["filename"]) + "/" + os.path.basename(options["filename"]).split(".")[0]
		
		options["input_path"] = os.path.abspath(path)

		path = options["input_path"]

		if "logs" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/logs/")

		if "alignments" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/alignments/")

		if "domains" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/domains/")

		if "surface" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/surface/")

		if "query" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/query/")

		if "results" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/results/")

		if "database" in os.listdir(path):
			pass
		else:
			os.mkdir(path + "/database/")

		options["output"] = path + "/" + os.path.basename(options["filename"]).split(".")[0] + ".rank"
		options["logfile"]= path + "/logs/" + time.strftime('%d%b%Y_%X_').replace(":","-") + os.path.basename(options["filename"]).split(".")[0] +  '.log'
		options["TEIRESIAS_output"] = path + "/" + os.path.basename(options["filename"]).split(".")[0] + ".fasta.out"
		options["TEIRESIAS_input"] = path + "/" + os.path.basename(options["filename"]).split(".")[0] + ".fasta"
		options["aminoAcidOccuranceFile"] = path + "/" + "PeptideOccurances.txt"
		options["BLAST_results"] = path + "/results/"
		options["filename"] = os.path.abspath(options["filename"])

		if options['Mask_path'] == 'mask.dat':
			options['Mask_path'] = os.path.join(os.path.dirname(sys.argv[0]),'mask.dat')
	
def removeFiles(options):
	if options['memory_saver'] == 'T':
		for file in os.listdir(options["input_path"]):
			if os.path.isfile(os.path.join(options["input_path"],file)):
				if file.split('.')[-1] not in ['fasta','rank','out']:
					os.remove(os.path.join(options["input_path"],file))
				else:
					pass
					
			elif os.path.isdir(os.path.join(options["input_path"],file)):
				try:
					shutil.rmtree(os.path.join(options["input_path"],file))
				except:
					pass
			else:
				pass
#----------------------------------------------------RUN---------------------------------------------------------------------#
#
#----------------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
	'''
	Main
	'''
	try:
		__version__ = "2.0"
		
		start = time.time()
		commandline = CommandLine()
		options = commandline.parseCommandLine()	
	
		if options["batch"] == 'F':
			setupFiles(options)
			
			logfile = open(options["logfile"],"w")
			sys.stderr = logfile	
			
			commandline.printOptions()
			del commandline
			
			runMain(options)
			
			logfile.close()
			removeFiles(options)
		elif  options["batch"] == 'T':
			
			try:
				os.listdir(options["filename"])
			except:
				print 'Input directory path does not exist'
				sys.exit()
				
			fileList = os.listdir(options["filename"])
			file_path = copy.deepcopy(options["filename"])
			for file in fileList:
				print '\n\n'
				print '+'*15 +'Running ' + file + '+'*15
				
				if file.split('.')[-1] != options['filetype']:
					pass
				else:
					options["filename"] = file_path + '/' + file
					setupFiles(options)
					
					logfile = open(options["logfile"],"w")
					sys.stderr = logfile
					
					try:
						runMain(options)
					except Exception,e:
						print 'Error running ' + file
						print e
					logfile.close()
					removeFiles(options)
	
			del commandline

	
		printLevel(1,options,"Total Time Elapsed\n%2.4f"%(time.time() - start))
		printLevel(1,options,"------------------------------------------------------------------\n")

		sys.stderr.write("End Time            : " + str(time.strftime('%X')) + "\n")
		sys.stderr.write("Total Time Elapsed  : %2.4f"%(time.time() - start)  + "\n")
		sys.stderr.write("----------------------------------------------------------------\n")
	

	except SystemExit:
		pass

	except Exception,e:
		error_type = sys.exc_info()[0]
		error_value = sys.exc_info()[1]
		error_traceback = traceback.extract_tb(sys.exc_info()[2])
	
		sys.stderr.write("\n")
		sys.stderr.write('Error in routine\n')
		sys.stderr.write('Error Type       : ' + str(error_type) + '\n')
		sys.stderr.write('Error Value      : ' + str(error_value) + '\n')
		sys.stderr.write('File             : ' + str(error_traceback[-1][0]) + '\n')	
		sys.stderr.write('Method           : ' + str(error_traceback[-1][2]) + '\n')	
		sys.stderr.write('Line             : ' + str(error_traceback[-1][1]) + '\n')		
		sys.stderr.write('Error            : ' + str(error_traceback[-1][3]) + '\n')	
		
		try:
			logfile.close()
			printLevel(1,options, open(options["logfile"],"r").read(-1))
		except:
			pass
	except:
		raise
#----------------------------------------------------------------------------------------------------------------------------#
