#!/usr/bin/python
# -*- coding: utf-8 -*-

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_busco
Description:  BUSCO data manipulation library
Version:      0.1.2
Last Edit:    10/01/22
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module will read in a BUSCO v5 MetaEuk results file and generate the missing *.fna files needed by BUSCOMP. It
    is designed to be incorporated into BUSCOMP but can be run standalone.

Commandline:
    ### ~ Core input/output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [$SEQIN basefile]
    busco5run=PATH  : Path to BUSCO v5 run (e.g. run_SGDR64.2.1/run_saccharomycetes_odb10/) [./]
    outdir=PATH     : Output directory for reformatted/processed files, if different to [busco5run=PATH]
    ### ~ MetaEuk conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    metaeukfna=T/F  : Perform metaeuk nucleotide busco sequences extraction [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_sequence
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added recognition and parsing of transcriptome mode output, which lacks Start and End.
    # 0.1.1 - Fixed bug with sequence names containing pipe characters. (Why?!)
    # 0.1.2 - Fixed bug with odd MetaEuk runs with appended letter on BUSCO ID.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add to BUSCOMP
    # [ ] : Fix bug when pipe character in sequence names - will need to shift indexing for sequence mapping.
    # [ ] : Future version might want to find and correct overlapping gene errors like this (if v5.2.2 still makes them):
        #WARN	00:00:00	Cannot find sequence for 8205at4891 sgdXVI_YEAST__BK006949:657600-661049
        #8205at4891	Duplicated	sgdXVI_YEAST__BK006949	657600	661049	+	523.5	440	https://www.orthodb.org/v10?query=8205at4891	NADPH-dependent diflavin oxidoreductase 1
        #8205at4891	Duplicated	sgdXVI_YEAST__BK006949	659187	661049	+	526.9	440	https://www.orthodb.org/v10?query=8205at4891	NADPH-dependent diflavin oxidoreductase 1
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_BUSCO', '0.1.2', 'January 2022', '2021')
    description = 'BUSCO data manipulation library'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
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
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: BUSCO Class                                                                                             #
#########################################################################################################################
class BUSCO(rje_obj.RJE_Object):
    '''
    BUSCO Class. Author: Rich Edwards (2021).

    Str:str
    - Busco5Run=PATH  : Path to BUSCO v5 run (e.g. run_SGDR64.2.1/run_saccharomycetes_odb10/) [./]
    - OutDir=PATH     : Output directory for reformatted/processed files, if different to [busco5run=PATH]
    - SeqIn

    Bool:boolean
    - MetaEukFNA=T/F  : Perform metaeuk nucleotide busco sequences extraction [False]

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB
    - SeqIn
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Busco5Run','OutDir','SeqIn']
        self.boollist = ['MetaEukFNA']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['DB','SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'Busco5Run':'./','OutDir':''})
        self.setBool({'MetaEukFNA':False})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
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
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                self._cmdReadList(cmd,'path',['Busco5Run','OutDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['MetaEukFNA'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MetaEukFNA'):
                if self.getStrLC('Busco5Run') and rje.exists(self.getStr('Busco5Run')):
                    return self.fnaFromBUSCOv5()
                else:
                    self.warnLog('metaeukfna=T but no buscorun=PATH found')
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def seqinObj(self,summarise=False,expect=False): ### Returns the a SeqList object for the SeqIn file
        '''
        Returns the a SeqList object for the SeqIn file.
        :return: self.obj['SeqIn']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['SeqIn'] or not self.obj['SeqIn'].seqNum():
                if not self.getStrLC('SeqIn'):
                    if expect: raise IOError('No seqin=FILE given')
                    else:
                        self.obj['SeqIn'] = None
                        return None
                seqcmd = self.cmd_list
                if summarise: seqcmd = ['summarise=T']+self.cmd_list
                self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
            if not self.obj['SeqIn'].seqNum():
                if expect: raise IOError('No sequences loaded from seqin=FILE!')
                else: self.warnLog('No sequences loaded from seqin=FILE!')
        except:
            self.errorLog('{0} seqinObj() error'.format(self.prog())); raise
        return self.obj['SeqIn']
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('OutDir'): self.setStr({'OutDir':self.getStr('Busco5Run')})
            self.obj['SeqIn'] = self.seqinObj()
            ### ~ [2] Check BUSCO Path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # .
            # ├── busco_sequences
            # │   ├── fragmented_busco_sequences
            # │   ├── multi_copy_busco_sequences
            # │   └── single_copy_busco_sequences
            # ├── full_table.tsv
            # ├── hmmer_output
            # │   ├── initial_run_results
            # │   └── rerun_results
            # ├── metaeuk_output
            # │   ├── combined_pred_proteins.fas
            # │   ├── initial_results
            # │   ├── refseq_db_rerun.faa
            # │   └── rerun_results
            # ├── missing_busco_list.tsv
            # └── short_summary.txt



            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def fnaFromBUSCOv5(self):      ### Loads BUSCOv5 metaeuk fasta files and full table and generates *.fna files
        '''
        Loads BUSCOv5 metaeuk fasta files and full table and generates *.fna files.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            shortfile = '{0}short_summary.txt'.format(self.getStr('Busco5Run'))
            mode = 'genome'
            if rje.exists(shortfile):
                for fline in open(shortfile,'r').readlines():
                    if rje.matchExp('# BUSCO was run in mode: (\S+)',fline):
                        mode = rje.matchExp('# BUSCO was run in mode: (\S+)',fline)[0]
            else:
                self.warnLog('Failed to read in {0}'.format(shortfile))
            self.printLog('#MODE','BUSCO was run in mode: {0}'.format(mode))
            ## ~ [1a] Load BUSCO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            filename = '{0}full_table.tsv'.format(self.getStr('Busco5Run'))
            v5head = ["BuscoID", "Status", "SeqName", "Start", "End", "Strand", "Score", "Length", "OrthoDBURL", "Description"]
            if mode == 'transcriptome':
                v5head = ["BuscoID", "Status", "SeqName", "Score", "Length", "OrthoDBURL", "Description"]
                buscodb = db.addTable(filename, mainkeys=['BuscoID','SeqName'], headers=v5head, ignore=['#'], name='Full')
            else:
                buscodb = db.addTable(filename, mainkeys=['BuscoID','SeqName','Start','End'], headers=v5head, ignore=['#'], name='Full')
            buscodb.indexReport('Status')
            ## ~ [1b] Load fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Load sequences and make (BuscoID,SeqName,Start,End):(seqname,sequence) dictionary
            seqdict = {}
            files = glob.glob('{0}metaeuk_output/initial_results/*.codon.fas'.format(self.getStr('Busco5Run')))
            files += glob.glob('{0}metaeuk_output/rerun_results/*.codon.fas'.format(self.getStr('Busco5Run')))
            if not files:
                raise IOError('Cannot find {0}metaeuk_output/*_results/*.codon.fas!'.format(self.getStr('Busco5Run')))
            for seqfile in files:
                seqcmd = self.cmd_list+['seqin={0}'.format(seqfile),'seqmode=list','autoload=T','summarise=T']
                for (seqname,sequence) in rje_seqlist.SeqList(self.log,seqcmd).seqs():
                    # >2618at4891_0|sgdI_YEAST__BK006935|-|452|2.166e-129|1|89888|92047|92047[92047]:89888[89888]:2160[2160]
                    seqdata = seqname.split('|')
                    buscoid = seqdata[0].split('_')[0]
                    seqname = seqdata[1]
                    if mode == 'transcriptome':
                        seqdict[(buscoid,seqname)] = (seqname,sequence)
                    else:
                        i = 2
                        # Adjust for fucking pipes in sequence names
                        while seqdata[i] not in '+-': i += 1
                        seqname = '|'.join(seqdata[1:i])
                        outname = '{0}:{1}-{2}'.format(seqname, seqdata[i+4], seqdata[i+5])
                        seqdict[(buscoid,seqname,seqdata[i+4],seqdata[i+5])] = (outname,sequence)
            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # xref to table and output to subdirectories
            rate2dir = {'Complete':'single_copy','Duplicated':'multi_copy','Fragmented':'fragmented'}
            for rating in ['Complete','Duplicated','Fragmented']:
                if rating not in buscodb.index('Status'):
                    self.printLog('\r#RATE', 'No {0} sequences.'.format(rating))
                    continue
                rkeys = buscodb.index('Status')[rating]
                self.progLog('\r#RATE','{0} {1} sequences...'.format(rje.iLen(rkeys),rating))
                self.deBug(rkeys)
                ## ~ [2a] Compile sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                buscoseq = {}  # Dictionary of {BuscoID:[(seqname,sequence)]
                bx = 0; nx = 0
                for bkey in rkeys:
                    self.debug(bkey)
                    if bkey not in seqdict and rje.matchExp('^(\S+):\d+-\d+',bkey[1]):
                        if len(bkey) == 4:
                            bkey = (bkey[0],rje.matchExp('^(\S+):\d+-\d+',bkey[1])[0],bkey[2],bkey[3])
                        else:
                            bkey = (bkey[0], rje.matchExp('^(\S+):\d+-\d+', bkey[1])[0])
                    if bkey not in seqdict and rje.matchExp('^(\S+\d+)\D+$', bkey[0]):
                        if len(bkey) == 4:
                            bkey = (rje.matchExp('^(\S+\d+)\D+$', bkey[0])[0],bkey[1],bkey[2],bkey[3])
                        else:
                            bkey = (rje.matchExp('^(\S+\d+)\D+$', bkey[0])[0], bkey[1])
                    if bkey in seqdict:
                        if bkey[0] not in buscoseq: buscoseq[bkey[0]] = []
                        buscoseq[bkey[0]].append(seqdict[bkey])
                        bx += 1
                    else:
                        if len(bkey) == 4:
                            self.warnLog('Cannot find sequence for {0} {1}:{2}-{3}'.format(bkey[0],bkey[1],bkey[2],bkey[3]))
                        else:
                            self.warnLog('Cannot find sequence for {0} {1}'.format(bkey[0],bkey[1]))
                        nx += 1
                self.printLog('\r#RATE','{0} sequences extracted for {1} {2} BUSCO IDs; {3} missing.'.format(bx,len(buscoseq),rating,nx))
                ## ~ [2b] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #X# Check translation? Manual check shows -ve strand is already reverse complemented
                outdir = '{0}busco_sequences/{1}_busco_sequences/'.format(self.getStr('OutDir'),rate2dir[rating])
                rje.mkDir(self,outdir); fx = 0; faax = 0
                for buscoid in rje.sortKeys(buscoseq):
                    fasta = ''
                    for (seqname,sequence) in buscoseq[buscoid]:
                        fasta += '>{0}\n{1}\n'.format(seqname,sequence)
                    fnafile = '{0}{1}.fna'.format(outdir,buscoid)
                    open(fnafile,'w').write(fasta)
                    fx += 1
                    faafile = '{0}{1}.faa'.format(outdir, buscoid)
                    if not rje.exists(faafile):
                        fasta = ''
                        for (seqname, sequence) in buscoseq[buscoid]:
                            fasta += '>{0}\n{1}\n'.format(seqname, rje_sequence.dna2prot(sequence,warnobj=self))
                        open(faafile, 'w').write(fasta)
                        faax += 1
                self.printLog('#FNA','{0} {1}*.fna files output.'.format(rje.iStr(fx),outdir))
                self.printLog('#FAA','{0} missing {1}*.faa files output.'.format(rje.iStr(faax),outdir))
            return True
        except:
            self.errorLog('%s.fnaFromBUSCOv5 error' % self.prog())
            return False
#########################################################################################################################
### End of SECTION II: BUSCO Class                                                                                      #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: BUSCO(mainlog,cmd_list+['tuplekeys=T']).run()

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
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
