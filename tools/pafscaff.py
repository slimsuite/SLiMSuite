#!/usr/bin/python

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
Module:       PAFScaff
Description:  Pairwise mApping Format reference-based scaffold anchoring and super-scaffolding.
Version:      0.4.1
Last Edit:    25/08/20
Citation:     Field et al. (2020), GigaScience 9(4):giaa027. [PMID: 32236524]
GitHub:       https://github.com/slimsuite/pafscaff
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    PAFScaff is designed for mapping genome assembly scaffolds to a closely-related chromosome-level reference genome
    assembly. It uses (or runs) [Minimap2](https://github.com/lh3/minimap2) to perform an efficient (if rough) all-
    against-all mapping, then parses the output to assign assembly scaffolds to reference chromosomes.

    Mapping is based on minimap2-aligned assembly scaffold ("Query") coverage against the reference chromosomes.
    Scaffolds are "placed" on the reference scaffold with most coverage. Any scaffolds failing to map onto any
    chromosome are rated as "Unplaced". For each reference chromosome, PAFScaff then "anchors" placed assembly
    scaffolds starting with the longest assembly scaffold. Each placed scaffold is then assessed in order of
    decreasing scaffold length. Any scaffolds that do not overlap with already anchored scaffolds in terms of the
    Reference chromosome positions they map onto are also considered "Anchored". if `newprefix=X` is set, scaffolds
    are renamed with the Reference chromosome they match onto. The original scaffold name and mapping details are
    included in the description. Unplaced scaffolds are not renamed.

    Finally, Anchored scaffolds are super-scaffolded by inserting gaps of `NnNnNnNnNn` sequence between anchored
    scaffolds. The lengths of these gaps are determined by the space between the reference positions, modified by
    overhanging query scaffold regions (min. length 10). The alternating case of these gaps makes them easy to
    identify later.

    ## Output

    PAFScaff outputs renamed, sorted and reoriented scaffolds in fasta format, along with mapping details:

    * `*.anchored.fasta`, `*.placed.fasta` and `*.unplaced.fasta` contain the relevant subsets of assembly scaffolds,
    renamed and/or reverse-complemented if appropriate.
    * `*.scaffolds.fasta` contains the super-scaffolded anchored scaffolds.
    * `*.scaffolds.tdt` contains the details of the PAFScaff mapping of scaffolds to chromosomes.
    * `*.log` contains run details, including any warnings or errors encountered.

    **NOTE:** The precise ordering, orientation and naming of the output scaffolds depends on the settings for:
    `refprefix=X newprefix=X sorted=T/F revcomp=T/F`.

    For full documentation of the PAFScaff workflow, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    pafin=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping; or run minimap2 [minimap2]
    basefile=STR    : Base for file outputs [PAFIN basefile]
    seqin=FASFILE   : Input genome assembly to map/scaffold onto $REFERENCE (minimap2 $ASSEMBLY) []
    reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
    assembly=FASFILE: As seqin=FASFILE
    refprefix=X     : Reference chromosome prefix. If None, will use all $REFERENCE scaffolds [None]
    newprefix=X     : Assembly chromosome prefix. If None, will not rename $ASSEMBLY scaffolds [None]
    unplaced=X      : Unplaced scaffold prefix. If None, will not rename unplaced $ASSEMBLY scaffolds [None]
    sorted=X        : Criterion for $ASSEMBLY scaffold sorting (QryLen/Coverage/RefStart/None) [QryLen]
    minmap=PERC     : Minimum percentage mapping to a chromosome for assignment [0.0]
    minpurity=PERC  : Minimum percentage "purity" for assignment to Ref chromosome [50.0]
    revcomp=T/F     : Whether to reverse complement relevant scaffolds to maximise concordance [True]
    scaffold=T/F    : Whether to "anchor" non-overlapping scaffolds by Coverage and then scaffold [True]
    dochtml=T/F     : Generate HTML PAFScaff documentation (*.info.html) instead of main run [False]
    pagsat=T/F      : Whether to output sequence names in special PAGSAT-compatible format [False]
    newchr=X        : Prefix for short PAGSAT sequence identifiers [ctg]
    spcode=X        : Species code for renaming assembly sequences in PAGSAT mode [PAFSCAFF]
    ### ~ Mapping/Classification options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    mmsecnum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [0]
    mmpcut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]
    mapopt=CDICT    : Dictionary of additional minimap2 options to apply (caution: over-rides conflicting settings) []
    ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_paf, rje_seqlist, rje_sequence, rje_rmd
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version with basic documentation. Added scaffold=T/F as an option.
    # 0.2.0 - Added sorted=X : Criterion for $ASSEMBLY scaffold sorting (QryLen/Coverage/RefStart/None) [QryLen]
    # 0.2.1 - Add documentation and fixed setting of Minimap2 N and p.
    # 0.3.0 - Added pagsat=T/F : Whether to output sequence names in special PAGSAT-compatible format [False]
    # 0.4.0 - Added purity criteria for more stringent assignment to chromosomes.
    # 0.4.1 - Fixed some issues with ambiguous scaffold output.
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
    # [Y] : Add option to actually scaffold with defined gap length or Reference gap length (and min gap length)
    # [Y] : Fix gap lengths to incorporate Qry overhangs.
    # [Y] : Make anchoring an option: alternatively output all placed scaffolds in start position order
    # [Y] : Add PreZero to scaffold numbering.
    # [Y] : Add option to name scaffolds in Size order (scaffold=T), or Start Position (scaffold=F)
    # [ ] : Add option to mask ends of sequences prior to scaffolding using KAT self-kmers.
    # [ ] : Finish the dochtml information.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('PAFScaff', '0.4.1', 'August 2020', '2019')
    description = 'Pairwise mApping Format reference-based scaffold anchoring and super-scaffolding'
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
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class PAFScaff(rje_obj.RJE_Object):
    '''
    PAFScaff Class. Author: Rich Edwards (2015).

    Str:str
    - NewChr=X        : Prefix for short PAGSAT sequence identifiers [ctg]
    - NewPrefix=X     : Assembly chromosome prefix. If None, will not rename $ASSEMBLY scaffolds [None]
    - PAFIn=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping []
    - Reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
    - RefPrefix=X     : Reference chromosome prefix. If None, will used all $REFERENCE scaffolds [None]
    - SeqIn=FASFILE   : Input genome to identify variants in ($ASSEMBLY) []
    - Sorted=X        : Criterion for $ASSEMBLY scaffold sorting (QryLen/Coverage/RefStart/None) [QryLen]
    - SpCode=X        : Species code for renaming assembly sequences in PAGSAT mode [PAFSCAFF]
    - Unplaced=X      : Unplaced scaffold prefix. If None, will not rename unplaced $ASSEMBLY scaffolds [None]


    Bool:boolean
    - DocHTML=T/F     : Generate HTML PAFScaff documentation (*.info.html) instead of main run [False]
    - PAGSAT=T/F      : Whether to output sequence names in special PAGSAT-compatible format [False]
    - Scaffold=T/F    : Whether to "anchore" non-overlapping scaffolds by Coverage and then scaffold [True]
    - Sorted=T/F      : Whether to sort $ASSEMBLY scaffold outputs [True]
    - RevComp=T/F     : Whether to reverse complement relevant scaffolds to maximise concordance [True]

    Int:integer
    - MMSecNum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [0]

    Num:float
    - MinMap=PERC     : Minimum percentage mapping to a chromosome for assignment [0.0]
    - MinPurity=PERC  : Minimum percentage "purity" for assignment to Ref chromosome [50.0]
    - MMPCut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - Assembly
    - DB
    - PAF
    - Reference
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['NewChr','NewPrefix','PAFIn','Reference','RefPrefix','SeqIn','Sorted','SpCode','Unplaced']
        self.boollist = ['DocHTML','PAGSAT','RevComp','Scaffold','Sorted']
        self.intlist = ['MMSecNum']
        self.numlist = ['MinMap','MinPurity','MMPCut']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['Assembly','PAF','Reference']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'NewChr,':'ctg','PAFIn':'minimap2','Sorted':'QryLen','SpCode':'PAFSCAFF'})
        self.setBool({'DocHTML':False,'PAGSAT':False,'Scaffold':True,'Sorted':True,'RevComp':True})
        self.setInt({'MMSecNum':0})
        self.setNum({'MMPCut':0.0,'MinMap':0.0,'MinPurity':50.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['PAF'] = rje_paf.PAF(self.log,['basefile=pafscaff','mmsecnum=0','pafin=minimap2']+self.cmd_list)
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
                self._cmdRead(cmd,type='file',att='SeqIn',arg='Assembly')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['NewChr','NewPrefix','RefPrefix','SpCode','Unplaced'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['PAFIn','Reference','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML','PAGSAT','RevComp','Scaffold'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MMSecNum'])   # Integers
                self._cmdReadList(cmd,'float',['MMPCut']) # Floats
                self._cmdReadList(cmd,'perc',['MinMap','MinPurity']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ## Tidy commands ##
        if self.getStrLC('Sorted') in ['','none','f','false']:
            self.setBool({'Sorted':False})
            self.printLog('#SORT','sorted=%s: no $ASSEMBLY sorting' % self.getStr('Sorted'))
        elif self.getStrLC('Sorted') not in ['qrylen','coverage','refstart']:
            self.warnLog('sorted=%s not recognised: set sorted=QryLen' % self.getStr('Sorted'))
            self.setStr({'Sorted':'QryLen'})
            self.setBool({'Sorted':True})
        else: self.setBool({'Sorted':True})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main PAFScaff run method
        '''
        # PAFScaff: Pairwise mApping Format reference-based scaffold anchoring and super-scaffolding.

        PAFScaff is designed for mapping genome assembly scaffolds to a closely-related chromosome-level reference genome
        assembly. It uses (or runs) [Minimap2](https://github.com/lh3/minimap2) to perform an efficient (if rough) all-
        against-all mapping, then parses the output to assign assembly scaffolds to reference chromosomes.

        Mapping is based on minimap2-aligned assembly scaffold ("Query") coverage against the reference chromosomes.
        Scaffolds are "placed" on the reference scaffold with most coverage. Any scaffolds failing to map onto any
        chromosome are rated as "Unplaced". For each reference chromosome, PAFScaff then "anchors" placed assembly
        scaffolds starting with the longest assembly scaffold. Each placed scaffold is then assessed in order of
        decreasing scaffold length. Any scaffolds that do not overlap with already anchored scaffolds in terms of the
        Reference chromosome positions they map onto are also considered "Anchored". if `newprefix=X` is set, scaffolds
        are renamed with the Reference chromosome they match onto. The original scaffold name and mapping details are
        included in the description. Unplaced scaffolds are not renamed.

        Finally, Anchored scaffolds are super-scaffolded by inserting gaps of `NnNnNnNnNn` sequence between anchored
        scaffolds. The lengths of these gaps are determined by the space between the reference positions, modified by
        overhanging query scaffold regions (min. length 10). The alternating case of these gaps makes them easy to
        identify later.

        ## Output

        PAFScaff outputs renamed, sorted and reoriented scaffolds in fasta format, along with mapping details:

        * `*.anchored.fasta`, `*.placed.fasta` and `*.unplaced.fasta` contain the relevant subsets of assembly scaffolds,
        renamed and/or reverse-complemented if appropriate.
        * `*.scaffolds.fasta` contains the super-scaffolded anchored scaffolds.
        * `*.scaffolds.tdt` contains the details of the PAFScaff mapping of scaffolds to chromosomes.
        * `*.log` contains run details, including any warnings or errors encountered.

        **NOTE:** The precise ordering, orientation and naming of the output scaffolds depends on the settings for:
        `refprefix=X newprefix=X sorted=T/F revcomp=T/F`. See main documentation (below) for details.

        ---

        # Running PAFScaff

        PAFScaff is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/pafscaff.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [PAFScaff git repo](https://github.com/slimsuite/pafscaff), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [PAFScaff git repo](https://github.com/slimsuite/pafscaff)
        for running on example data.

        For mapping prior to parsing, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
        environment `$PATH` or given to PAFScaff with the `minimap2=PROG` setting.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        pafin=PAFFILE   : PAF generated from $REFERENCE $ASSEMBLY mapping; or run minimap2 [minimap2]
        basefile=STR    : Base for file outputs [PAFIN basefile]
        seqin=FASFILE   : Input genome assembly to map/scaffold onto $REFERENCE (minimap2 $ASSEMBLY) []
        reference=FILE  : Fasta (with accession numbers matching Locus IDs) ($REFERENCE) []
        assembly=FASFILE: As seqin=FASFILE
        refprefix=X     : Reference chromosome prefix. If None, will used all $REFERENCE scaffolds [None]
        newprefix=X     : Assembly chromosome prefix. If None, will not rename $ASSEMBLY scaffolds [None]
        unplaced=X      : Unplaced scaffold prefix. If None, will not rename unplaced $ASSEMBLY scaffolds [None]
        sorted=X        : Criterion for $ASSEMBLY scaffold sorting (QryLen/Coverage/RefStart/None) [QryLen]
        minmap=PERC     : Minimum percentage mapping to a chromosome for assignment [0.0]
        minpurity=PERC  : Minimum percentage "purity" for assignment to Ref chromosome [50.0]
        revcomp=T/F     : Whether to reverse complement relevant scaffolds to maximise concordance [True]
        scaffold=T/F    : Whether to "anchor" non-overlapping scaffolds by Coverage and then scaffold [True]
        dochtml=T/F     : Generate HTML PAFScaff documentation (*.info.html) instead of main run [False]
        pagsat=T/F      : Whether to output sequence names in special PAGSAT-compatible format [False]
        newchr=X        : Prefix for short PAGSAT sequence identifiers [ctg]
        spcode=X        : Species code for renaming assembly sequences in PAGSAT mode [PAFSCAFF]
        ### ~ Mapping/Classification options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        minimap2=PROG   : Full path to run minimap2 [minimap2]
        mmsecnum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [0]
        mmpcut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]
        mapopt=CDICT    : Dictionary of additional minimap2 options to apply (caution: over-rides conflicting settings) []
        ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ## PAFScaff input and output

        ### Reference genome naming

        _Details of reference genome naming requirements to be added. Chromosomes should start with the prefix given by_
        `refprefix=X`.

        ### Generating the PAF input file

        The PAF file is generated by running Minimap2 without secondary alignments. See `rje_paf.py` for details.

        ## Anchoring and mapping

        The PAF file is first parsed and the following fields extracted:

        * Qry QryLen QryStart QryEnd Strand Ref RefLen RefStart RefEnd Identity Length

        If the `refprefix=X` setting was provided, hits are then cleaned up to (a) strip the prefix from chromosome names
        and convert numerical chromosomes to numbers for proper sorting, and (b) drop any unplaced reference scaffolds.
        If no `refprefix` is given, no reference scaffold renaming/filtering will be performed.

        For each hit, `Coverage` is then calculated as `QryEnd`-`QryStart`+1. This `Coverage` stat is ultimately used for
        ranking the assignments of assembly scaffolds to reference chromosomes/scaffolds. For each Qry/Ref/Strand
        combination, the `Coverage`, `Identity` and `Length` statistics are summed, and a new field `N` added with the
        number of individual alignments. The top ranked reference strand is selected for each query, and two more fields
        added:

        * `Inv` = Query coverage on the opposite strand to the best reference strand
        * `Trans` = Summed query coverage for all other reference chromosomes/scaffolds.

        Together, these are used to establish the total percentage of query coverage that is scaffolded, versus mapping
        to a different reference sequence. These are converted into:

        * `Purity` = 100 * (`Coverage` + `Inv`) / (`Coverage` + `Inv` + `Trans`)

        v0.4.0 has introduced a couple of additional parameters than can be used to increase the stringency of any
        mapping. This is mainly for the purpose of reducing situtations where highly repetitive multi-mapping sequences
        are assigned to a single chromosome. By default, `minpurity=50.0`, meaning that at least half of the chromosome
        mapping should be to the main reference chromosome.

        * minmap=PERC   : Minimum percentage mapping to a chromosome for assignment [0.0]
        * minpurity=PERC: Minimum percentage "purity" for assignment to Ref chromosome [50.0]

        Once queries have been assigned to reference scaffolds, they are then ordered according to the reference
        scaffold and start position of the match to that scaffold. Queries are also sorted for the purposes of renaming,
        using the `sorted=X` statistic. By default, this is query length, but could be switched to `Coverage` or `RefStart`.

        NOTE: This ordering is quite crude and assumes the reference start position is part of the main mapping block. If
        the query has a fragment mapping to an upstream repeat sequence, for example, this ordering could be messed up.
        Future releases will check and try to fix this. Visualisation of the mapping is recommended.

        ### Anchoring and Scaffolding

        Scaffolds mapped to a reference chromosome are designated `Placed`. `Anchored` scaffolds are a special subset of
        `Placed` scaffolds, which can be combined into a super-scaffold. Any scaffolds failing to map onto a scaffold
        will be designated `Unplaced`.

        `Anchor` assignment is only performed if `scaffold=T`. Assembly scaffolds are processed in decreasing size and
        assigned to the `Anchor` set if either (a) no previous scaffold has been assigned to the mapped reference
        chromosome, or (b) the scaffold mapping positions on the reference chromosome do not overlap with any previous
        `Anchor` scaffold for that chromosome. As with the ordering (above), this scaffolding will get messed up if the
        earlier mapped scaffolds have fragments a great distance up or downstream of the main mapped region.

        If `scaffold=T`, the 'Anchor` sequences for each reference chromosome will be scaffolded by concatenating them
        in reference position order, reverse-complementing where the main mapping was to the negative strand. Gaps will
        be inserted between scaffolds, consisting of a stretch of `Nn` repeats that matches the size of the gap between
        reference positions (e.g. the end of the one match and the start of the next), with a minimum gap of 10 nt.

        NOTE: This explicitly expects chromosome-level reference scaffolds. It will not scaffold the reference using the
        assembly to maximise the overall scaffolding. It is purely for chromosome assignment and then additional
        scaffolding based on assumed synteny.

        ### Sequence naming

        If `newprefix=X` is given, the sequence name will be replaced with the prefix and reference chromosome identifier
        (parsed using `refprefix=X`). If multiple assembly sequences map to the same reference chromosome, these will be
        numbered `.1`, `.2` .. `.N`. (NOTE: these numbers will have zero prefixes to make them sort correctly, e.g. 10+
        sequences will start `.01`, 100+ sequences will start `.001` etc.

        Sequence descriptions take the form:

            SEQNAME len=SEQLEN COV% REF(STRAND) START:END; [INV% REF(INVSTRAND);] [OTHER% other;]

        where:

        * `SEQNAME` = original assembly sequence name.
        * `SEQLEN` = assembly sequence length.
        * `COV%` = percentage of query mapping onto reference strand.
        * `REF(STRAND)` = reference chromosome and strand. If mapped to the -ve strand, the description will have a `RevComp` prefix.
        * `START:END` = start and end positions on reference chromosome of query mapping.
        * `INV% REF(INVSTRAND)` = percentage mapped onto inverse strand of reference chromosome. (Omitted if 0%.)
        * `OTHER%` = percentage mapped onto other reference chromosome. (Omitted if 0%.)

        ### Output

        Assembly mapping an scaffolding details will be save to `$BASEFILE.scaffolds.tdt`.

        * _Details of output fields to follow._

        Sequences will be output to:

        * `$BASEFILE.anchored.fasta`
        * `$BASEFILE.placed.fasta`
        * `$BASEFILE.unplaced.fasta`

        If `scaffold=T`, the scaffolded assembly will be saved to `$BASEFILE.scaffolded.fasta` and the sequences
        summaried in  `$BASEFILE.scaffolded.tdt`.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.getBool('DocHTML'): return self.docHTML()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.parsePAF()
            if not self.scaffold(): raise ValueError('Scaffolding error')
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            paf = self.obj['PAF']
            if self.getBool('DocHTML'): paf.setInt({'Interactive':-1})
            paf.setup()
            self.baseFile(paf.baseFile())
            self.obj['Assembly'] = paf.obj['Assembly']
            self.obj['Reference'] = paf.obj['Reference']
            self.obj['DB'] = paf.obj['DB']
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def parsePAF(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            paf = self.obj['PAF']
            if paf.getStrLC('PAFIn') in ['minimap','minimap2']:
                paf.setStr({'PAFIn':'%s.paf' % paf.baseFile()})
                self.printLog('#PAFIN','Minimap2 PAF file set: %s' % paf.getStr('PAFIn'))
                pafopt = {'p':self.getNum('MMPCut'),'N':self.getInt('MMSecNum')}
                paf.dict['MapOpt'] = rje.combineDict(paf.dict['MapOpt'],pafopt,overwrite=True)
                paf.minimap2()
            paf.parsePAF()
            return True     # Setup successful
        except: self.errorLog('Problem during %s parsePAF.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def scaffold(self): ### Scaffold $ASSEMBLY using $REFERENCE mapping from PAF file.
        '''
        Scaffold $ASSEMBLY using $REFERENCE mapping from PAF file.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basename = self.baseFile()
            self.headLog('Assign scaffolds to reference chromosomes')
            #i# pafhead = string.split('Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length Quality')
            #i# pafaln = string.split('tp cm s1 s2 NM MD AS ms nn ts cg cs dv')
            pafdb = self.db('paf')
            pafdb.keepFields(string.split('# Qry QryLen QryStart QryEnd Strand Hit SbjLen SbjStart SbjEnd Identity Length'))
            for field in ['Hit','SbjLen','SbjStart','SbjEnd']: pafdb.renameField(field,'Ref%s' % field[3:])
            revcomp = self.getBool('RevComp')
            minmap = self.getNum('MinMap')
            minpurity = self.getNum('MinPurity')
            ## ~ [1a] Cleanup Reference chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('RefPrefix'):
                refpref = self.getStr('RefPrefix')
                conversion = {}
                for refchr in pafdb.index('Ref'):
                    if rje.matchExp('^%s(\d+)$' % refpref,refchr):
                        chr = rje.matchExp('^%s(\d+)$' % refpref,refchr)[0]
                        conversion[refchr] = string.atoi(chr)
                    elif rje.matchExp('^%s(\S+)' % refpref,refchr):
                        chr = rje.matchExp('^%s(\S+)' % refpref,refchr)[0]
                        conversion[refchr] = chr
                pafdb.dropEntriesDirect('Ref',conversion.keys(),inverse=True)
                for refchr in pafdb.index('Ref'):
                    for entry in pafdb.indexEntries('Ref',refchr): entry['Ref'] = conversion[refchr]

            ### ~ [2] Compress PAF hits and assign Assembly scaffolds to strands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pafdb.makeField('QryEnd-QryStart+1','Coverage')
            rankfield = 'Coverage'  # Could offer option to rank by Identity or Length
            rules = {'QryLen':'max','QryStart':'min','QryEnd':'max','RefLen':'max','RefStart':'min','RefEnd':'max','Identity':'sum','Length':'sum','Coverage':'sum','N':'sum'}
            pafdb.addField('N',evalue=1)
            pafdb.compress(['Qry','Ref','Strand'],rules)
            pafdb.dropField('#')
            ## ~ [2a] Rank Qry hits by Length to find dominant reference strand ~~~~~~~~~~~~~~~~~~~ ##
            pafdb.rankFieldByIndex('Qry',rankfield,newfield='Rank',rev=True,absolute=True,lowest=True,unique=True,warn=True)
            #i# Rank=1 should now be the best
            pafdb.addField('Inv',evalue=0)    # This is the total Length on the inverse strand
            pafdb.addField('Trans',evalue=0)  # This is the total Length on the non-dominant scaffolds
            qpaf = {}   # Dictionary of Qry:Rank 1 entry
            for entry in pafdb.entries():
                if entry['Rank'] == 1: qpaf[entry['Qry']] = entry
            for entry in pafdb.entries():
                if entry['Rank'] == 1: continue
                qentry = qpaf[entry['Qry']]
                if entry['Ref'] == qentry['Ref']: qentry['Inv'] += entry['Coverage']
                else: qentry['Trans'] += entry['Coverage']
            ## ~ [2b] Reduce Qry hits to dominant reference strand ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pafdb.dropEntriesDirect('Rank',[1],inverse=True)
            ## ~ [2c] Assign ambiguous ratings based on coverage and/or purity filters ~~~~~~~~~~~~ ##
            pafdb.addField('Purity',after='Trans')
            pafdb.addField('RefMap',after='RefEnd')
            for entry in pafdb.entries(sorted=True):
                htot = entry['Coverage'] + entry['Inv'] + entry['Trans']
                entry['Purity'] = 100.0 * (entry['Coverage'] + entry['Inv']) / htot
                hperc = 100.0 * entry['Coverage'] / entry['QryLen']
                if hperc < minmap or entry['Purity'] < minpurity:
                    entry['Ref'] = 'Ambig-%s' % entry['Ref']
                    entry['RefMap'] = 'Ambiguous'

            ### ~ [3] Report Scaffolding ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Rename %s scaffolds' % basename)
            newpref = None
            if self.getStrLC('NewPrefix'): newpref = self.getStr('NewPrefix')
            pafdb.newKey(['Ref','RefStart','RefEnd','Qry'])
            #?# Do we ever want to name the scaffolds by position rather than size? Add a naming option?
            rankfield = 'Qry'
            if self.getBool('Sorted'):
                for field in ['QryLen','Coverage','RefStart']:
                    if field.lower() == self.getStrLC('Sorted'):
                        rankfield = field
            #i# Add RefN, which counts the hits in rankfield order for each Reference Scaffold
            self.printLog('#REFN','Ranking by %s for Reference Chromosome subnumbering' % rankfield)
            pafdb.rankFieldByIndex('Ref',rankfield,newfield='RefN',rev=rankfield in ['QryLen','Coverage'],absolute=True,lowest=True,unique=True,warn=True)
            #x# rankfield = 'RefStart'
            #x# pafdb.rankFieldByIndex('Ref',rankfield,newfield='RefN',rev=False,absolute=True,lowest=True,unique=True,warn=False)
            pafdb.makeField('#Qry#','Scaffold',after='Qry')
            pafdb.addField('Description',after='Scaffold')
            for entry in pafdb.entries(sorted=True):
                invstrand = {'+':'-','-':'+'}[entry['Strand']]
                htot = entry['Coverage'] + entry['Inv'] + entry['Trans']
                hperc = 100.0 * entry['Coverage'] / entry['QryLen']
                iperc = 100.0 * entry['Inv'] / entry['QryLen']
                tperc = 100.0 * entry['Trans'] / entry['QryLen']
                cperc = 100.0 * htot / entry['QryLen']
                ref = entry['Ref']
                if entry['RefMap'] == 'Ambiguous':
                    ref = entry['Ref'][6:]
                    entry['Description'] = '%s len=%s; Ambiguous %.2f%% %s(%s) %s:%s;' % (entry['Qry'],rje_seqlist.dnaLen(entry['QryLen'],dp=0,sf=4),hperc,ref,entry['Strand'],rje.iStr(entry['RefStart']),rje.iStr(entry['RefEnd']))
                else:
                    entry['Description'] = '%s len=%s; %.2f%% %s(%s) %s:%s;' % (entry['Qry'],rje_seqlist.dnaLen(entry['QryLen'],dp=0,sf=4),hperc,entry['Ref'],entry['Strand'],rje.iStr(entry['RefStart']),rje.iStr(entry['RefEnd']))
                if iperc > 0: entry['Description'] += ' %.2f%% %s(%s);' % (iperc,ref,invstrand)
                if tperc > 0: entry['Description'] += ' %.2f%% other;' % (tperc)
                if revcomp and entry['Strand'] == '-': entry['Description'] = 'RevComp %s' % entry['Description']
                if newpref and entry['RefMap'] != 'Ambiguous':
                    entry['Qry'] = '%s%s' % (newpref,entry['Ref'])
                    refcount = len(pafdb.index('Ref',entry['Ref']))
                    if refcount > 1 or self.getBool('PAGSAT'): entry['Qry'] = '%s.%s' % (entry['Qry'],rje.preZero(entry['RefN'],refcount))
                    if self.getBool('PAGSAT'): entry['Qry'] = '%s%s.%s_%s__%s' % (self.getStr('NewChr'),entry['Ref'],rje.preZero(entry['RefN'],refcount),self.getStr('SpCode'),entry['Qry'])
                self.printLog('#CHRMAP','%s: %s' % (entry['Qry'],entry['Description']))
            ## ~ [3a] Classify into anchored, placed, and unplaced ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Anchored scaffolds are a special subset of placed scaffolds, which can be combined into a super-scaffold
            #i# NOTE: This explicitly expects chromosome-level reference scaffolds. It will not scaffold the reference
            #i# using the assembly to maximise the overall scaffolding. It is purely for chromosome assignment and then
            #i# additional scaffolding based on assumed synteny.
            anchored = {}  # Dictionary fof {RefChr:List of entries to scaffold}
            placed = {} # Dictionary fof {RefChr:List of entries that overlap scaffolded entries}
            #i# Assembly scaffolds are processed in decreasing size.
            #!# Consider making this a commandline option: QryLen/Coverage/Length/Identity
            if self.getBool('Scaffold'):
                scaffsort = 'QryLen'
                for entry in pafdb.sortedEntries(scaffsort,reverse=True):
                    if entry['RefMap'] == 'Ambiguous': continue
                    refchr = entry['Ref']
                    reflen = entry['RefLen']
                    if refchr not in anchored: anchored[refchr] = [{'RefStart':0,'RefEnd':0},{'RefStart':reflen+1,'RefEnd':reflen+1}]; placed[refchr] = []
                    i = 1
                    while len(anchored[refchr]) > i:
                        # Look whether entry can insert between i-0 and i
                        ref1 = anchored[refchr][i-1]
                        ref2 = anchored[refchr][i]
                        if entry['RefStart'] > ref1['RefEnd'] and ref2['RefStart'] > entry['RefEnd']:
                            entry['RefMap'] = 'Anchored'
                            anchored[refchr].insert(i,entry); break
                        # Look for overlap with i
                        if entry['RefStart'] <= ref1['RefEnd'] and ref1['RefStart'] <= entry['RefEnd']:
                            entry['RefMap'] = 'Placed'
                            placed[refchr].append(entry); break
                        # Move to next
                        i += 1
            #i# If scaffold=F, output all placed scaffolds in RefStart order.
            else:
                for entry in pafdb.sortedEntries('RefN',reverse=False):
                    if entry['RefMap'] == 'Ambiguous': continue
                    refchr = entry['Ref']
                    if refchr not in placed: placed[refchr] = []
                    placed[refchr].append(entry)
                    entry['RefMap'] = 'Placed'

            ### ~ [4] Output anchored, placed and unplaced assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Output %s scaffolds' % basename)
            pafdb.newKey(['Qry'])
            #!# Add option to update descriptions with more meaningful descriptions and species
            self.progLog('#OUT','Preparing fasta output...')
            assembly = self.obj['Assembly']
            seqdict = assembly.makeSeqNameDic('short')
            ax = 0
            if anchored:
                afile = '%s.anchored.fasta' % self.baseFile()
                rje.backup(self,afile)
                AFILE = open(afile,'w')
            pfile = '%s.placed.fasta' % self.baseFile(); px = 0
            rje.backup(self,pfile)
            PFILE = open(pfile,'w')
            ufile = '%s.unplaced.fasta' % self.baseFile(); ux = 0
            rje.backup(self,ufile)
            UFILE = open(ufile,'w')
            ## ~ [4a] Sorted output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Sorted'):
                self.progLog('#OUT','Sorted scaffold output... ')
                for entry in pafdb.entries(sorted=True):
                    if entry['RefMap'] == 'Ambiguous': continue
                    aseq = seqdict[entry['Scaffold']]
                    sequence = assembly.seqSequence(aseq)
                    if revcomp and entry['Strand'] == '-': sequence = rje_sequence.reverseComplement(sequence)
                    if entry['RefMap'] == 'Anchored':
                        AFILE.write('>%s %s\n%s\n' % (entry['Qry'],entry['Description'],sequence))
                        ax += 1
                    else:
                        PFILE.write('>%s %s\n%s\n' % (entry['Qry'],entry['Description'],sequence))
                        px += 1
            #i# Sort scaffolds by length for output
            scaffnum = {}   # Dictionary of {sname:scaffnumber}
            scaffseq = []   # Sequences for scaffold output
            if self.getStrLC('Unplaced') or self.getBool('PAGSAT'):
                scaffsort = []
                for sname in rje.sortKeys(seqdict):
                    if sname not in pafdb.index('Scaffold') or pafdb.indexEntries('Scaffold',sname)[0]['RefMap'] == 'Ambiguous':
                        aseq = assembly.getSeq(seqdict[sname],'tuple')
                        if self.getBool('Sorted'): scaffsort.append((-assembly.seqLen(aseq),sname))
                        else: scaffsort.append((sname,sname))
                scaffsort.sort()
                for (sorted,sname) in scaffsort:
                    ux += 1
                    scaffnum[sname] = ux
                    scaffseq.append(sname)
            #i# Output scaffolds
            aseq = assembly.nextSeq()
            if scaffnum: aseq = assembly.getSeq(seqdict[scaffseq.pop(0)],'tuple')
            while aseq:
                if self.getBool('Sorted'):
                    self.progLog('#OUT','Unplaced scaffold output...')
                else: self.progLog('#OUT','Unsorted scaffold output...')
                (aname,sequence) = aseq
                sname = string.split(aname)[0]
                ambentry = pafdb.data(sname)
                if sname not in pafdb.index('Scaffold') or pafdb.indexEntries('Scaffold',sname)[0]['RefMap'] == 'Ambiguous':
                    #UFILE.write(assembly.fasta[aseq])
                    #ux += 1
                    ux = scaffnum[sname]
                    if self.getStrLC('Unplaced') or self.getBool('PAGSAT'):
                        newname = '%s%s' % (self.getStr('Unplaced'),rje.preZero(ux,assembly.seqNum()))
                        if self.getBool('PAGSAT'):
                            newname = '%sUn.%s_%s__%s' %  (self.getStr('NewChr'),rje.preZero(ux,assembly.seqNum()),self.getStr('SpCode'),newname)
                        stype = 'contig'
                        if 'NNNNNNNNNN' in sequence.upper(): stype = 'scaffold'
                        UFILE.write('>%s %s len=%s; Unplaced %s\n%s\n' % (newname,aname,rje_seqlist.dnaLen(len(sequence),dp=0,sf=4),stype,sequence))
                        if ambentry:
                            ambentry['Qry'] = newname
                            pafdb.dict['Data'][newname] = pafdb.dict['Data'].pop(sname)
                        else: pafdb.addEntry({'Qry':newname,'Scaffold':aname,'RefMap':'Unplaced','QryLen':len(sequence)})
                    else:
                        UFILE.write('>%s\n%s\n' % (aname,sequence))
                        if not ambentry: pafdb.addEntry({'Qry':aname,'Scaffold':aname,'RefMap':'Unplaced','QryLen':len(sequence)})
                elif not self.getBool('Sorted'):
                    entry = pafdb.indexEntries('Scaffold',sname)[0]
                    if revcomp and entry['Strand'] == '-': sequence = rje_sequence.reverseComplement(sequence)
                    if entry['RefMap'] == 'Anchored':
                        AFILE.write('>%s %s\n%s\n' % (entry['Qry'],entry['Description'],sequence))
                        ax += 1
                    else:
                        PFILE.write('>%s %s\n%s\n' % (entry['Qry'],entry['Description'],sequence))
                        px += 1
                if scaffnum:
                    if scaffseq: aseq = assembly.getSeq(seqdict[scaffseq.pop(0)],'tuple')
                    else: aseq = None
                else: aseq = assembly.nextSeq()
            if anchored:
                AFILE.close()
                self.printLog('#OUT','%s anchored scaffolds output to: %s' % (rje.iStr(ax),afile))
            PFILE.close()
            self.printLog('#OUT','%s placed scaffolds output to: %s' % (rje.iStr(px),pfile))
            UFILE.close()
            self.printLog('#OUT','%s unplaced scaffolds output to: %s' % (rje.iStr(ux),ufile))
            ## ~ [4b] Add sequence summary of assigned sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if anchored: rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % afile,'seqmode=file','summarise','dna','autoload'])
            if px: rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % pfile,'seqmode=file','summarise','dna','autoload'])
            if ux: rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % ufile,'seqmode=file','summarise','dna','autoload'])
            ## ~ [4c] Save table, including unplaced scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pafdb.dropField('RefN')
            pafdb.saveToFile('%s.scaffolds.tdt' % self.baseFile())



            ### ~ [5] Scaffold assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add option to actually scaffold with defined gap length or Reference gap length (and min gap length)
            #i# Outputting scaffolded assembly to *.scaffolded.fasta not *.scaffolds.fasta to avoid tdt conflict when summarising
            if self.getBool('Scaffold'):
                mingap = 10
                sfile = '%s.scaffolded.fasta' % self.baseFile()
                rje.backup(self,sfile)
                SFILE = open(sfile,'w'); sx = 0; gx = 0; ngap = 0
                for refchr in rje.sortKeys(anchored):
                    refname = ['%s%s' % (newpref,refchr)]
                    refseq = ''
                    prevend = 0
                    prevqry = 0
                    for entry in anchored[refchr][1:-1]:
                        refname.append(string.split(entry['Description'],';')[0] + ';')
                        aseq = seqdict[entry['Scaffold']]
                        sequence = assembly.seqSequence(aseq)
                        if entry['Strand'] == '-': sequence = rje_sequence.reverseComplement(sequence)
                        if prevend:
                            gapx = entry['RefStart'] - prevend - 1 - prevqry
                            gapx = max(mingap,gapx)
                            gapseq = 'Nn' * (gapx / 2)
                            if rje.isOdd(gapx): gapseq += 'N'
                            if len(gapseq) != gapx: raise ValueError
                            refseq += gapseq; ngap += 1
                        refseq += sequence
                        prevend = entry['RefEnd']
                        if entry['Strand'] == '-':
                            prevqry = entry['QryStart'] - 1
                        else:
                            prevqry = entry['QryLen'] - entry['QryEnd']
                    refname = string.join(refname)
                    if len(anchored[refchr]) > 3:
                        refname += ' [PAFScaff reference-based scaffold]'
                        gx += 1
                    SFILE.write('>%s\n%s\n' % (refname,refseq))
                    sx += 1
                SFILE.close()
                self.printLog('#GAP','%s reference-scaffolded super-scaffolds generated: %s gap regions inserted' % (rje.iStr(gx),rje.iStr(ngap)))
                self.printLog('#OUT','%s scaffold sequences output to: %s' % (rje.iStr(sx),sfile))
                ## ~ [5a] Add sequence summary of scaffolded sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % sfile,'seqmode=db','summarise','dna','autoload'])



            return True
        except: self.errorLog('%s.scaffold error' % self.prog()); return False
#########################################################################################################################
    def docHTML(self):  ### Generate the PAFScaff Rmd and HTML documents.                                        # v0.1.0
        '''Generate the PAFScaff Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='PAFScaff Documentation',author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2019 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown PAFScaff documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: PAFScaff Class                                                                                   #
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
    try: PAFScaff(mainlog,cmd_list+['tuplekeys=T']).run()

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
