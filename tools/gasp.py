#!/usr/bin/python

# See below for name and description
# Copyright (C) 2004 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      GASP
Description:  Gapped Ancestral Sequence Prediction
Version:      2.0.0
Last Edit:    08/03/18
Citation:     Edwards & Shields (2004), BMC Bioinformatics 5(1):123. [PMID: 15350199]
Copyright (C) 2004  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is purely for running the GASP algorithm as a standalone program. All the main functionality is encoded in
    the modules listed below. The GASP algorithm itself is now encoded in the rje_ancseq module. GASP may also be run
    interactively from the rje_tree module.

Input/Output:
    seqin=FILE      : Fasta/ClustalW/Phylip Format sequence file [infile.fas].
    nsfin=FILE      : Newick Standard Format treefile [infile.nsf].
    indeltree=T/F   : Whether to output a text tree describing indel patterns [True]

GASP Model Options:
    fixpam=INT      : PAM distance fixed to X [0].
    rarecut=NUM     : Rare aa cut-off [0.05].
    fixup=T/F       : Fix AAs on way up (keep probabilities) [True].
    fixdown=T/F     : Fix AAs on initial pass down tree [False].
    ordered=T/F     : Order ancestral sequence output by node number [False].
    pamtree=T/F     : Calculate and output ancestral tree with PAM distances [True].
    desconly=T/F    : Limits ancestral AAs to those found in descendants [True].
    xpass=INT       : How many extra passes to make down & up tree after initial GASP [1].

"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj
import rje_ancseq, rje_seq, rje_tree   # rje_ancseq now sorts out rje_pam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Improved version with second pass.
    # 1.1 - Improved OO. Restriction to descendant AAs. (Good for BAD etc.)
    # 1.2 - No Out Object in Objects
    # 1.3 - Added more interactive load options
    # 1.4 - Minor tweaks to imports.
    # 2.0.0 - Upgraded to rje_obj framework for REST server.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Improved Interactive menus for module and major activities - currently use rje_tree.py. Add to rje_ancseq./py?
    # [ ] : Convert to newer module style.
    # [ ] : Take interactive running out of rje_tree and reinstate GASP as top-level module?
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('GASP', '2.0.0', 'March 2016', '2004')
    description = 'Gapped Ancestral Sequence Prediction'
    author = 'Dr Richard J. Edwards.'
    comments = ['Cite: Edwards & Shields (2004), BMC Bioinformatics 5(1):123. [PMID: 15350199]',
                 'Please report bugs to Richard.Edwards@UNSW.edu.au']
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class GASP(rje_obj.RJE_Object):
    '''
    GASP Class. Author: Rich Edwards (2018).

    Str:str
    - SeqIn=FILE      : Fasta/ClustalW/Phylip Format sequence file [infile.fas].
    - NSFIn=FILE      : Newick Standard Format treefile [infile.nsf].

    Bool:boolean
    - IndelTree=T/F  : Prints a text tree describing indel patterns.

    Int:integer

    Num:float

    File:file handles with matching str filenames

    List:list

    Dict:dictionary

    Obj:RJE_Objects
    - GASP = rje_ancseq.Gasp object
    - Tree - rje_tree.Tree object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['SeqIn','NSFIn']
        self.boollist = ['IndelTree']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'SeqIn':'infile.fas'})
        self.setBool({'IndelTree':True})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
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
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path
                self._cmdReadList(cmd,'file',['SeqIn','NSFIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['IndelTree'])  # True/False Booleans
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
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.gasp()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Read in Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = self.getStr('SeqIn')
            seqs = rje_seq.SeqList(log=self.log,cmd_list=['i=0']+self.cmd_list+['autofilter=F','autoload=F','seqin=None'])
            self.printLog('#SEQS','Loading sequences from %s' % seqfile)
            if not seqs.loadSeqs(seqfile=seqfile,seqtype='protein',aln=True):
                raise IOError('Cannot load from %s' % seqfile)
            seqfile = seqs.info['Name']
            basefile = rje.baseFile(seqfile)
            if not self.getStrLC('Basefile'): self.baseFile(basefile)
            self.printLog('#SEQ',"%s protein sequences read from %s\n" % (str(seqs.seqNum()),seqfile),1)
            #?# Add option to generate alignment?
            self.printLog('#SEQ',"Alignment = %s. (%d aa)\n" % (seqs.opt['Aligned'],seqs.seq[0].seqLen()),1)
            self.dict['Output']['seqin'] = seqfile

            ### ~ [1] Read in Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('NSFIn'):
                nsfin = self.getStr('NSFIn')
            else:
                nsfin = basefile + '.nsf'
            while not os.path.exists(nsfin):
                if self.i() >= 0:
                    nsfin = rje.choice(text='Input tree file "%s" not found. Input filename? (Blank to exit.)' % nsfin)
                    if nsfin == '':
                        raise KeyboardInterrupt
                else:
                    raise IOError('File %s not found. Cannot load tree!' % nsfin)
            self.dict['Output']['nsfin'] = nsfin
            self.cmd_list.append('nsfin=' + nsfin)
            self.printLog('#TREE','Loading tree from %s' % nsfin)
            self.obj['Tree'] = mytree = rje_tree.Tree(log=self.log,cmd_list=['root=yes']+self.cmd_list)
            mytree.mapSeq(seqlist=seqs)
            mytree.textTree()
            if mytree.opt['ReRooted']:
                mytree.saveTree(filename='%s.nsf' % basefile)
            return True     # Setup successful
        except KeyboardInterrupt: self.printLog('#CANCEL','User terminated.'); return False
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        seqin = Input sequence alignment [fas]
        nsfin = Input phylogenetic tree (Newick) [nsf]
        indeltree = [Optional] Text tree describing indel patterns. [txt]
        anc.fas = Ancestral sequence prediction alignment [fas]
        anc.txt = Text tree with ancestral nodes named and numbered [txt]
        anc.nsf = Newick version of tree with ancestral nodes named and numbered [nsf]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['seqin','nsfin','indeltree','anc.fas','anc.txt','anc.nsf']
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def gasp(self):      ### Main GASP Method, copied from GASP v1.4
        '''
        Main GASP Method, copied from GASP v1.4.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mytree = self.obj['Tree']
            indeltree = self.getBool('IndelTree')
            if self.baseFile() == 'infile': self.baseFile('gasp')
            ### ~ [2] GASP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not indeltree or not mytree.node[-1].obj['Sequence']:  # Perform GASP
                self.printLog('#SEQ','GASP: Gapped Ancestral Sequence Prediction')
                mygasp = rje_ancseq.Gasp(tree=mytree,ancfile='%s' % self.baseFile(),cmd_list=self.cmd_list,log=self.log)
                self.verbose(0,2,'%s' % mygasp.details(),1)
                if self.i() > 0:
                    if not rje.yesNo('Use these parameters?'):
                        mygasp.edit()
                mygasp.gasp()
                self.printLog('#GASP',"GASP run completed OK!")
                self.dict['Output']['anc.fas'] = '%s.anc.fas' % self.baseFile()
                self.dict['Output']['anc.nsf'] = '%s.anc.nsf' % self.baseFile()
                self.dict['Output']['anc.txt'] = '%s.anc.txt' % self.baseFile()
            ### ~ [3] InDel Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if indeltree:
                self.dict['Output']['indeltree'] = '%s.indel.txt' % self.baseFile()
                mytree.indelTree(filename='%s.indel.txt' % self.baseFile())
        except: self.errorLog('%s.gasp error' % self.prog())
#########################################################################################################################
### End of SECTION II: GASP Class                                                                                       #
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
    try: GASP(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################

