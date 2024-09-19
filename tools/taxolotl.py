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
Module:       taxolotl
Description:  Taxolotl genome assembly taxonomy summary and assessment tool
Version:      0.1.4
Last Edit:    09/05/23
Citation:     Tobias PA, Edwards RJ, et al. (preprint) bioRxiv 2022.07.29.502101 doi: 10.1101/2022.07.29.502101
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    Taxolotl combines the MMseqs2 `easy-taxonomy` with GFF parsing to perform taxonomic analysis of a genome assembly
    (and any subsets given by `taxsubsets=LIST`) using an annotated proteome. Taxonomic assignments are mapped onto
    genes as well as assembly scaffolds and (if `assembly=FILE` is given) contigs.

    See <https://slimsuite.github.io/taxolotl/> for details. General SLiMSuite run documentation can be found at
    <https://github.com/slimsuite/SLiMSuite>.

    Taxolotl is available as part of SLiMSuite, or via a standalone GitHub repo at
    <https://github.com/slimsuite/taxolotl>.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Protein annotation file to assess [annotation.faa]
    gffin=FILE      : Protein annotation GFF file [annotation.gff]
    cdsin=FILE      : Optional transcript annotation file for renaming and/or longest isoform extraction [annotation.fna]
    assembly=FILE   : Optional genome fasta file (required for some outputs) [None]
    basefile=X      : Prefix for output files [$SEQBASE]
    gffgene=X       : Label for GFF gene feature type ['gene']
    gffcds=X        : Label for GFF CDS feature type ['CDS']
    gffmrna=X       : Label for GFF mRNA feature type ['mRNA']
    taxlevels=LIST  : List of taxonomic levels to report (* for superkingdom and below) ['*']
    ### ~ Run mode options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    dochtml=T/F     : Generate HTML Taxolotl documentation (*.docs.html) instead of main run [False]
    ### ~ Taxolotl options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    taxdb=FILE      : MMseqs2 taxonomy database for taxonomy assignment [seqTaxDB]
    taxbase=X       : Output prefix for taxonomy output [$SEQBASE.$TAXADB]
    taxorfs=T/F     : Whether to generate ORFs from assembly if no seqin=FILE given [True]
    taxbyseq=T/F    : Whether to parse and generate taxonomy output for each assembly (GFF) sequence [True]
    taxbycontig=T/F : Whether to generate taxonomy output for each contig if the assembly is loaded [True]
    taxbyseqfull=T/F: Whether generate full easy taxonomy report outputs for each assembly (GFF) sequence [False]
    taxsubsets=FILELIST : Files (fasta/id) with sets of assembly input sequences (matching GFF) to summarise []
    taxwarnrank=X   : Taxonomic rank (and above) to warn when deviating for consensus [family]
    bestlineage=T/F : Whether to enforce a single lineage for best taxa ratings [True]
    mintaxnum=INT   : Minimum gene count in main dataset to keep taxon, else merge with higher level [2]
    ### ~ TabReport options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    tabreport=FILE  : Convert MMseqs2 report into taxonomy table with counts (if True use taxbase=X) [None]
    taxhigh=X       : Highest taxonomic level for tabreport [class]
    taxlow=X        : Lowest taxonomic level for tabreport [species]
    taxpart=T/F     : Whether to output entries with partial taxonomic levels to tabreport [False]
    ### ~ System options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    tmpdir=PATH     : Temporary directory path for running mmseqs2 [./tmp/]
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
import rje, rje_db, rje_gff, rje_obj, rje_rmd, rje_seqlist
import rje_kat, saaga   # Taxolotl class will inherit these classes
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added tabreport function.
    # 0.1.1 - Fix bug with contig output. Added seqname, start and end to contig summary.
    # 0.1.2 - Added input file checking.
    # 0.1.3 - Updated citation data.
    # 0.1.4 - Fixed a // bug for Python3.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [ ] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add SAAGA taxonomy.
    # [ ] : Add KAT.
    # [ ] : Add GC and CpGratio calculations.
    # [ ] : Add depth calculation. (SAMTools?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Taxolotl', '0.1.4', 'May 2023', '2021')
    description = 'Taxolotl genome assembly taxonomy summary and assessment tool'
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
### SECTION II: Taxolotl Class                                                                                          #
#########################################################################################################################
class Taxolotl(saaga.SAAGA, rje_kat.KAT):
    '''
    Taxolotl Class. Author: Rich Edwards (2021).

    Str:str
    - Assembly=FILE   : Optional genome fasta file (required for some outputs) [None]
    - CDSIn=FILE      : Optional transcript annotation file [annotation.fna]
    - GFFIn=FILE      : Protein annotation GFF file [annotation.gff]
    - GFFGene=X       : Label for GFF gene feature type ['gene']
    - GFFCDS=X        : Label for GFF CDS feature type ['CDS']
    - GFFmRNA=X       : Label for GFF mRNA feature type ['mRNA']
    - SeqIn = Protein annotation file to assess [annotation.faa]
    - TabReport=FILE  : Convert MMseqs2 report into taxonomy table with counts (if True use taxbase=X) [None]
    - TaxBase=X       : Output prefix for taxonomy output [$SEQBASE.$TAXADB]
    - TaxDB=FILE      : MMseqs2 taxonomy database for taxonomy assignment [None]
    - TaxHigh=X       : Highest taxonomic level for tabreport [class]
    - TaxLow=X        : Lowest taxonomic level for tabreport [species]
    - TaxWarnRank=X   : Taxonomic rank (and above) to warn when deviating for consensus [family]
    - TmpDir=PATH     : Temporary directory path for running mmseqs2 [./tmp/]

    Bool:boolean
    - BestLineage=T/F : Whether to enforce a single lineage for best taxa ratings [True]
    - DocHTML=T/F     : Generate HTML Taxolotl documentation (*.info.html) instead of main run [False]
    - TaxByContig=T/F : Whether to generate taxonomy output for each contig if the assembly is loaded [True]
    - TaxBySeq=T/F    : Whether to parse and generate taxonomy output for each assembly (GFF) sequence [True]
    - TaxBySeqFull=T/F: Whether generate full easy taxonomy report outputs for each assembly (GFF) sequence [False]
    - Taxonomy=T/F    : Summarise taxonomic assignments for contamination assessments [False]
    - TaxORFs=T/F     : Whether to generate ORFs from assembly if no seqin=FILE given [True]
    - TaxPart=T/F     : Whether to output entries with partial taxonomic levels to tabreport [False]

    Int:integer
    - MinTaxNum=INT   : Minimum gene count in main dataset to keep taxon, else merge with higher level [2]

    Num:float

    File:file handles with matching str filenames

    List:list
    - TaxLevels=LIST  : List of taxonomic levels to report (* for superkingdom and below) ['*']
    - TaxRanks : Actual taxonomic levels used in the taxonomy analysis (high to low)
    - TaxSubsets=FILELIST : Files (fasta/id) with sets of assembly input sequences (matching GFF) to summarise []

    Dict:dictionary
    - TaxLineage = taxtuple:list of (taxrank,taxid,taxname) in high to low order
    - TaxTuple = taxid:(taxrank,taxid,taxname)

    Obj:RJE_Objects
    - Assembly = SeqList: Optional genome assembly
    - CDSIn = SeqList: Optional transcript annotation file
    - DB = Database object
    - GFF = rje_gff.GFF
    - Rmd = RMarkdown control object
    - SeqIn = SeqList: Protein annotation file to assess [predicted_cds.fasta]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Assembly','CDSIn','GFFIn','GFFGene','GFFCDS','GFFmRNA','SeqIn','TabReport','TaxHigh','TaxLow','TaxBase','TaxDB','TaxWarnRank','TmpDir']
        self.boollist = ['BestLineage','DocHTML','TaxByContig','TaxBySeq','TaxBySeqFull','Taxonomy','TaxORFs','TaxPart']
        self.intlist = ['MinTaxNum']
        self.numlist = []
        self.filelist = []
        self.listlist = ['TaxLevels','TaxRanks','TaxSubsets']
        self.dictlist = ['TaxLineage','TaxTuple']
        self.objlist = ['Assembly','SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'CDSIn':'annotation.fna','GFFGene':'gene','GFFCDS':'CDS','GFFmRNA':'mRNA',
                     'GFFIn':'annotation.gff','SeqIn':'annotation.fasta',
                     'TaxWarnRank':'family','TaxHigh':'class','TaxLow':'species',
                     'TmpDir':'./tmp/'})
        self.setBool({'BestLineage':True,'DocHTML':False,'TaxByContig':True,'TaxBySeq':True,'TaxBySeqFull':False,'Taxonomy':True,'TaxORFs':True})
        self.setInt({'MinTaxNum':2,'TopHits':250})
        self.setNum({'MinGlobID':40.0})
        self.list['TaxLevels'] = ['kingdom','phylum','class','order','family','genus','species']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
        self.obj['GFF'] = rje_gff.GFF(self.log,['attfield=T']+self.cmd_list+['tuplekeys=T','warnfield='])
        self.obj['GFF'].obj['DB'] = self.obj['DB']
        self.obj['Rmd'] = rje_rmd.Rmd(self.log,self.cmd_list)
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
                self._katCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['GFFGene','GFFCDS','GFFmRNA','TaxHigh','TaxLow','TaxWarnRank'])   # Normal strings
                self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['Assembly','CDSIn','GFFIn','SeqIn','TabReport','TaxBase','TaxDB'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BestLineage','DocHTML','TaxByContig','TaxBySeq','TaxBySeqFull','TaxORFs','TaxPart'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinTaxNum'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'perc',['']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['TaxLevels'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['TaxSubsets']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if not self.baseFile(return_none=''): self.baseFile('taxolotl')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # Taxolotl: Taxolotl genome assembly taxonomy summary and assessment tool

        Taxolotl combines the MMseqs2 `easy-taxonomy` with GFF parsing to perform taxonomic analysis of a genome assembly
        (and any subsets given by `taxsubsets=LIST`) using an annotated proteome. Taxonomic assignments are mapped onto
        genes as well as assembly scaffolds and (if `assembly=FILE` is given) contigs.

        See <https://slimsuite.github.io/taxolotl/>  and the documentation below for details.
        General SLiMSuite run documentation can be found at <https://github.com/slimsuite/SLiMSuite>.

        Taxolotl is available as part of SLiMSuite, or via a standalone GitHub repo at
        <https://github.com/slimsuite/taxolotl>.

        ## Citing Taxolotl

        Taxolotl is currently unpublished. Please cite the GitHub page and this bioRxiv paper, which has an example of
        Taxolotl in action:

        * Tobias PA, Edwards RJ, Surana P, Mangelson H, Inácio V, do Céu Silva M, Várzea V, Park RF & Batista D.
        "A chromosome-level genome resource for studying virulence mechanisms and evolution of the coffee rust pathogen
        _Hemileia vastatrix_. bioRxiv 2022.07.29.502101 doi: [10.1101/2022.07.29.502101](https://doi.org/10.1101/2022.07.29.502101)

        ---

        # Running Taxolotl

        Taxolotl is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/taxolotl.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [Taxolotl git repo](https://github.com/slimsuite/taxolotl), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [Taxolotl git repo](https://github.com/slimsuite/taxolotl)
        for running on example data.

        [MMseqs2](https://github.com/soedinglab/MMseqs2) must be installed and either added to the environment `$PATH`.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Protein annotation file to assess [annotation.faa]
        gffin=FILE      : Protein annotation GFF file [annotation.gff]
        cdsin=FILE      : Optional transcript annotation file [annotation.fna]
        assembly=FILE   : Optional genome fasta file (required for some outputs) [None]
        basefile=X      : Prefix for output files [taxolotl]
        gffgene=X       : Label for GFF gene feature type ['gene']
        gffcds=X        : Label for GFF CDS feature type ['CDS']
        gffmrna=X       : Label for GFF mRNA feature type ['mRNA']
        gffdesc=X       : GFF output field label for annotated proteins (e.g. note, product) [product]
        taxlevels=LIST  : List of taxonomic levels to report (* for superkingdom and below) ['*']
        ### ~ Run mode options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        dochtml=T/F     : Generate HTML Taxolotl documentation (*.docs.html) instead of main run [False]
        ### ~ Taxonomy options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        taxdb=FILE      : MMseqs2 taxonomy database for taxonomy assignment [None]
        taxbase=X       : Output prefix for taxonomy output [$SEQBASE.$TAXADB]
        taxorfs=T/F     : Whether to generate ORFs from assembly if no seqin=FILE given [True]
        taxbyseq=T/F    : Whether to parse and generate taxonomy output for each assembly (GFF) sequence [True]
        taxbyseqfull=T/F: Whether generate full easy taxonomy report outputs for each assembly (GFF) sequence [False]
        taxsubsets=FILELIST : Files (fasta/id) with sets of assembly input sequences (matching GFF) to summarise []
        taxwarnrank=X   : Taxonomic rank (and above) to warn when deviating for consensus [family]
        bestlineage=T/F : Whether to enforce a single lineage for best taxa ratings [True]
        mintaxnum=INT   : Minimum gene count in main dataset to keep taxon, else merge with higher level [2]
        ### ~ TabReport options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        tabreport=FILE  : Convert MMseqs2 report into taxonomy table with counts (if True use taxbase=X) [None]
        taxhigh=X       : Highest taxonomic level for tabreport [class]
        taxlow=X        : Lowest taxonomic level for tabreport [species]
        taxpart=T/F     : Whether to output entries with partial taxonomic levels to tabreport [False]
        ### ~ System options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
        tmpdir=PATH     : Temporary directory path for running mmseqs2 [./tmp/]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---


        ## Taxolotl overview

        The first step is to run MMseqs2:

            mmseqs easy-taxonomy $PROTEOME $TAXDB $TAXBASE $TMPDIR

        Where `$PROTEOME` is the proteome provided with `seqin=FILE`, `$TAXDB` is a MMseqs2 taxonomic database
        (see below for creation), provided with `taxdb=FILE`, `$TAXBASE` is the `easy-taxonomy` output prefix, and
        `$TMPDIR` is the temporary directory (default `tmp`). If pre-existing results exist (`$TAXBASE._report` and
        `$TAXBASE_lca.tsv`) then these will be loaded, unless `force=T` is set. If MMseqs2 is not installed, pre-computed
        results *must* be provided. In principle, `report` and `lca.tsv` files generate by other tools should work as long
        as the format is the same.

        The core of Taxolotl is the MMSeqs2 "Lowest Common Ancestor" (LCA) assignment, in which each sequence is associated
        with the lowest unabmigious taxonomic rank possible. Where amibiguity exists, a sequence will be assigned to a
        higher level. Higher levels also receive all the taxonomic assignments of their daughter taxa, and so the sequence
        count for any given taxonomic group will always be equal or greater than its lower subdivisions. Conceptually,
        Taxolotl separates out the counts into `taxnum`, which are counts at that level or below, and `taxpure`, which are
        the numbers assigned specifically to that level. (i.e. `taxnum` will be the sum of `taxpure` for that taxonomic
        group and all lower divisions.) See the MMseqs2 documentation for more details.

        Taxolotl will first read in the `*_report` file to build its internal taxonomy tree for the samples. By default, mmseqs will report all possible taxonomic levels, and Taxolotl will retain the following:

            species, species subgroup, species group, subgenus, genus, subtribe, tribe, subfamily, family, superfamily, parvorder, infraorder, suborder, order, superorder, infraclass, subclass, class, superclass, subphylum, phylum, superphylum, subkingdom, kingdom, superkingdom

        This can be reduced further by specifying a subset of taxonomic levels of interest with `taxlevels=LIST`. Any missing levels, along with
        "no rank" or "clade" taxa (except `unclassified`, `root`, and `cellular organisms`), will be mapped to the next highest taxonomic level. Any MMseqs2 assignments to that level will be transferred to the higher level. Any taxa failing to meet the `mintaxnum=INT` threshold (default=2) will also be mapped onto higher levels.

        Next, the `*_lca.tsv` file is read and mapped onto the `gffin=FILE` GFF file to assign proteins to genes and
        sequences. The lowest-level hit for each gene will be kept, remapping to `taxlevels` as required. These
        collated ratings will be output to `*.lca_genes.tsv` and `*.lca_genes.gff` Gene ratings are then summed for each assembly sequence, and the dominant
        classification for each taxonomic level established for (a) each sequence, and (b) the whole dataset. Full
        collated ratings will be output to `*.taxolotl_report.tsv`. Ratings per sequence are output to `*.taxbyseq.tsv`. Dominant taxa are reported in the log file as `#BEST` entries.

        To flag contamination, each sequence is assessed against the dominant taxonomic rating at each taxonomic level.
        The percentage of genes matching each dominant rating is reported for each sequence in `*.taxolotl.tsv`
        along with the number of genes with a rating at that level, separated with a `|`. This will exclude any genes
        without ratings at that taxonomic level. A `:consensus:` entry will also report the overall values for the whole
        assembly.

        Any sequences that have a dominant taxonomic label deviating from the overall consensus at any ranking levels
        set by `taxwarnrank=X` (default family) or above will raise a contamination warning and be output in the log file with a `#BADTAX` rating. These sequences will have their dominant taxon and it's
        precentage appended to the consensus percentage, also separated by `|`. For example, `25.00|20|Chordata|50.00`
        would indicate that 25% of the 20 genes with ratings at that level matched the consensus, whilst the dominant
        classification was `Chordata` with 50% of 20 rated genes assigned to this category. Such sequences will also have `badtax` rating in the `rating` field of `*.taxolotl.tsv`. Sequences matching the dominant taxa will have a `goodtax` rating, whilst sequences without any genes mapped onto taxa by MMseqs2 will be rated `notax`.

        Good, Bad and missing sequence counts will be summarised in the log file in `#BEST`, `BADTAX`, and `#NOTAX` entries.
        Sequence subsets are output to `*.id` and `*.fasta` files, and summarised along with the full assembly in
        `*.seqsummary.tsv`. (Any ratings without sequences will not be output/summarised.) If `assembly=FILE` is provided,
        sequences without genes will also be summarised. Taxonomy ratings for these subsets are also output to
        `*.$RATING.taxolotl_report.tsv` files. Any sequence subsets provided by `taxsubsets=LIST` (see below) will also be
        summarised in `*.$SUBSET.taxolotl_report.tsv` files. It is recommended that all the MMseqs2 `_report` file is loaded
        with all the `*.taxolotl_report.tsv` for visualisation with [Pavian](https://github.com/fbreitwieser/pavian)
        (Breitwieser FP and Salzberg SL (2020) [Bioinformatics 36(4):1303-1304](https://doi.org/10.1093/bioinformatics/btz715))
        through its [Shiny App](https://fbreitwieser.shinyapps.io/pavian/).

        Finally, if `assembly=FILE` is provided (unless `taxbycontig=F`), contigs will be extracted by splitting scaffolds on `mingap=INT` (default 10) consecutive `N`s. Genes will be remapped onto contigs as with sequences, and taxonomic ratings output to `*.taxbyctg.tsv` and `*.ctgtaxolotl.tsv`. These are the contig equivalents of `*.taxbyseq.tsv` and `*.taxolotl.tsv`. Contigs without taxonomic ratings will be listed in the log file with `#BADTAX` entries, unless already reported as an assembly sequence.

        ## Main taxonomy outputs

        Outputs will be given a file prefix set by `taxbase=X`. By default, this will be `$SEQBASE.$TAXADB`, where
        `$SEQBASE` is the basename of `seqin=FILE` and `$TAXADB` is the taxonomy database set by `taxdb=FILE`.

        The main mmseqs `easy-taxonomy` output will generate:

        * `*_lca.tsv` = best assignments per protein sequence (protein, taxid, rank, taxname): required.
        * `*_report` = text summary of overall taxonomy that can be loaded by Pavian etc.: required.
        * `*_tophit_aln` = top database hits for each protein (not currently used): not required.
        * `*_tophit_report` = taxonomic classification of the top hit proteins: not required.

        In addition, Taxolotl will output:

        * `*.taxbyseq.tsv` = Rating counts for each taxonomic group by assembly sequence (scaffold).
        * `*.taxolotl_report.tsv` = Collated Kraken-style report file.
        * `*.lca_genes.tsv` = Best assignments (lowest taxonomic level) for each gene.
        * `*.lca_genes.gff` = GFF file with Taxolotl ratings for each gene.
        * `*.taxolotl.tsv` = Tab separated file with consensus taxonomic assignment at each taxonomic rank, and ratings per sequence.
        * `*.$SUBSET.id` = Sequence identifiers for assembly subsets based on Taxolotl ratings.
        * `*.$SUBSET.fasta` = Fasta files of assembly subsets based on Taxolotl ratings.
        * `*.seqsummary.tsv` = Summary statistics for assembly subset fasta files.
        * `*.taxbyctg.tsv` = Rating counts for each taxonomic group by assembly contig.
        * `*.ctgtaxolotl.tsv` = Taxolotl ratings by assembly contig.

        ### Taxonomy by sequence output

        If `taxbyseq=T` then an additional `*.taxbyseq.tsv` file will be produced, with the following fields:

        * `seqname` = assembly sequence name
        * `genenum` = number of genes parsed for that sequence
        * `protnum` = number of proteins parsed for that sequence
        * `rank` = taxonomic rank of rating
        * `genetax` = number of genes with assignment at that level
        * `taxid` = taxonomic label identifier number
        * `taxname` = taxonomic label name at that rank
        * `taxperc` = percentage assignment to this rank or lower
        * `taxnum` = number of genes assigned to this rank or lower
        * `taxpure` = number of genes assigned to this rank specifically

        ## Sequence subset analysis

        In addition to the main output for the whole proteome, any subsets given by `taxsubsets=LIST` will have their own `*.taxolotl_report.tsv` file, which can be visualised with Pavian. These must be lists of IDs that match the assembly sequence names in the GFF file. Subsets will be named after the subset file prefix, e.g. `assembly.suspect.id` would generate `*.assembly.suspect.taxolotl_report.tsv`.


        ## Generating a taxonomic database

        Please see the MMseqs2 documentation for generating a taxonomic database. To date, Taxolotl has been tested with taxonomy databases generated from NCBI nr, using BLAST+ and MMSeqs2 and the NCBI taxonomy dump (<https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>):

        ```
        blastdbcmd -db $NCBIPATH/nr -entry all > ncbinr.faa
        blastdbcmd -db $NCBIPATH/nr -entry all -outfmt "%a %T" > ncbinr.faa.taxidmapping

        mmseqs createdb ncbinr.faa ncbinr.faaDB
        mmseqs createtaxdb ncbinr.faaDB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file ncbinr.faa.taxidmapping
        mmseqs createindex ncbinr.faaDB tmp
        ```

        If the assembly itself is already in RefSeq, it is recommended that the taxa of the assembly is removed before running Taxolotl.


        ## Simple ORF mode

        If no proteins are given, ORFs will be generated by `SeqSuite` with default settings `minorf=100 rftran=6 terminorf=50 orfgaps=F`, i.e. ORFs of 100+ amino acids from all six reading frames, or 50+ amino acids if truncated at the end of a sequence. ORFs will not span assembly gaps, and any ambiguous (`X`) translations will be replaced with stop codons (`*`), unless `orfgaps=T` is set. Note that, due to introns, it is expected that these ORFs will often represent partial coding sequences, and many will be random junk translations.

        The idea of ORF mode is to provide a quick, crude impression of the taxonomic profile. However, for large assemblies it can be very slow to process.

        In ORF mode, each ORF is assumed to represent a different gene, although this may not be the case. Currently, `SeqSuite` will not generate a GFF file for the ORFs. As a result, the `taxbycontig` output is not available.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# The setup will parse the input sequences and GFF
            if self.getBool('DocHTML'): return self.docHTML()
            if self.getStrLC('TaxDB') and not self.setup(): self.printLog('#ABORT','Problem during setup: aborted'); return False
            elif not self.getStrLC('TaxDB') and not self.setTaxBase(): self.printLog('#ABORT','Problem during setTaxBase: aborted'); return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxbase = self.getStr('TaxBase')
            if self.getStrLC('TabReport') in ['t','true']:
                self.setStr({'TabReport':'{0}_report'.format(taxbase)})
            if self.getStrLC('TabReport'):
                self.tabReport()
            if self.getStrLC('TaxDB'):
                self.headLog('TAXOLOTL', line='=')
                if not self.taxonomy(): return False
            else:
                self.printLog('#TAXDB','No TaxDB=FILE set')
            #!# Add running and parsing of KAT and the GC content and CpGratio (and CAI?)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def runMode(self,checkmodes,allmodes=False):   # Check run modes and return True if any active or False if not
        '''
        Check run modes and return True if any active or False if not
        :param checkmodes: list of modes to check. Will print settings to log if none given
        :param allmodes: whether all modes given must be set to return True.
        :return: True/False
        '''
        if not checkmodes: return False
        for modestr in checkmodes:
            if modestr.lower() in ['taxonomy','taxolotl']: return True
        return False
#########################################################################################################################
    def pureTaxonomy(self): return True
#########################################################################################################################
    def setTaxBase(self):   ### Sets self.getStr('TaxBase')
        '''Sets self.getStr('TaxBase')'''
        if not self.getStrLC('Basefile'): self.baseFile(rje.baseFile(self.getStr('TabReport'), strip_path=True))
        if not self.getStrLC('TaxBase'): self.setStr({'TaxBase':self.getStr('TabReport')})
        self.printLog('#BASE', self.basefile())
        self.printLog('#TXBASE', self.getStr('TaxBase'))
        return True
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def tabReport(self):    ### Convert *_report file into taxonomy count report.
        '''
        Convert *_report file into taxonomy count report.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('TABREPORT', line='=')
            reportfile = self.getStr('TabReport')
            taxbase = self.getStr('TaxBase')
            taxranks = self.setTaxRanks()
            self.list['TaxRanks'] = taxranks = ['no rank'] + taxranks  #i# Add 'no rank' for root and 'cellular organisms'
            db = self.db()
            ## ~ [1a] Setup tabreport output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            try: ti = taxranks.index(self.getStrLC('TaxHigh'))
            except: raise ValueError('Problem finding TaxHigh "{0}" in taxonomic levels'.format(self.getStrLC('TaxHigh')))
            try: tj = taxranks.index(self.getStrLC('TaxLow'))
            except: raise ValueError('Problem finding TaxLow "{0}" in taxonomic levels'.format(self.getStrLC('TaxLow')))
            if ti > tj:
                self.warnLog('{0} > {1} => swapping TaxHigh and TaxLow'.format(self.getStrLC('TaxLow'),self.getStrLC('TaxHigh')))
                self.setStr({'TaxHigh':self.getStrLC('TaxLow'),'TaxLow':self.getStrLC('TaxHigh')})
                (ti,tj) = (tj,ti)
            self.debug(ti)
            self.debug(tj)
            taxhead = ['taxid'] + taxranks[ti:tj+1] + ['taxnum']
            self.debug(taxhead)
            tabdb = db.addEmptyTable('tabreport',taxhead,['taxid'])
            ### ~ [2] Parse easy-taxonomy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Parse easy-taxonomy', line='-')
            repdb = db.addTable(reportfile,mainkeys=['#'],datakeys='All',delimit='\t',headers=['taxperc','taxnum','taxpure','taxrank','taxid','taxname'],ignore=[''],name='report',expect=True)
            #i# Keep TaxID as strings
            repdb.dataFormat({'taxperc':'float','taxnum':'int','taxpure':'int'}) #,'taxid':'int'})
            taxtree = {}    # {(rank,id,name):{taxtree}}
            levels = []     # List of (rank,id,name) tuples corresponding to current level being parsed
            #X# taxmapping = {} # {(rank,id,name):(rank,id,name)} -> Will be used to reclassify lca table
            taxmapping = {} # {id:id} -> Will be used to reclassify lca table
            taxcounts = {}  # {id:taxnum}
            taxlineage = self.dict['TaxLineage'] = {} # {(rank,id,name):[list of higher levels]}
            taxtuple = self.dict['TaxTuple'] = {} # {id:(rank,id,name)}
            taxparent = {}  # {id: parent id}
            ## ~ [2a] Parse taxonomy tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mintaxnum = self.getInt('MinTaxNum') # Minimum gene count in main dataset to keep taxon, else merge with higher level [2]
            self.debug(mintaxnum)
            rx = 0.0; rtot = repdb.entryNum()
            for rentry in repdb.entries(sorted=True):
                self.progLog('\r#REPORT','Parsing taxonomy: {0:.2f}%'.format(rx/rtot)); rx += 100.0
                #self.bugPrint(rentry)
                #i# Set up root
                if rentry['taxname'] == 'unclassified':
                    taxtuple['0'] = ('no rank','0','unclassified')
                    taxcounts['0'] = rentry['taxnum']
                    taxlineage[('no rank','0','unclassified')] = []
                    continue
                if rentry['taxname'] == 'root':
                    levels = [('no rank','1','root')]
                    taxtree[('no rank','1','root')] = {}
                    taxtuple['1'] = ('no rank','1','root')
                    taxcounts['1'] = rentry['taxnum']
                    taxlineage[('no rank','1','root')] = []
                    continue
                #i# Parse levels
                i = len(rje.matchExp('^(\s+)',rentry['taxname'])[0]) / 2
                levels = levels[:i]
                #if rx < 1000: self.debug('Level {1}::{0}\n'.format(levels,i))
                mytree = taxtree
                for level in levels:
                    if level in mytree:
                        mytree = mytree[level]
                    elif level[0] in taxranks:
                        raise ValueError('{0} missing from {1}'.format(level,rje.sortKeys(mytree)))
                taxon = (rentry['taxrank'],rentry['taxid'],rentry['taxname'][i*2:])
                if taxon in mytree: raise ValueError
                #i# Map taxa based on absence from taxranks
                if taxon[0] not in taxranks and taxon[2] not in ['unclassified', 'root', 'cellular organisms']:
                    i = len(levels) - 1
                    while levels[i][0] not in taxranks: i -= 1
                    taxmapping[taxon[1]] = levels[i][1]
                    #i# Check for additional mapping due to low occurrence levels
                    while taxmapping[taxon[1]] in taxmapping:
                        taxmapping[taxon[1]] = taxmapping[taxmapping[taxon[1]]]
                    #self.debug(taxmapping)
                #!# Add additional mapping due to small numbers
                elif rentry['taxnum'] < mintaxnum:
                    #i# Just map up and then iteratively map
                    taxmapping[taxon[1]] = levels[-1][1]
                    while taxmapping[taxon[1]] in taxmapping:
                        taxmapping[taxon[1]] = taxmapping[taxmapping[taxon[1]]]
                    mytree[taxon] = {}
                else:
                    mytree[taxon] = {}
                    taxcounts[taxon[1]] = rentry['taxnum']
                taxparent[taxon[1]] = levels[-1][1]
                while taxparent[taxon[1]] in taxmapping:
                    taxparent[taxon[1]] = taxmapping[taxparent[taxon[1]]]
                taxlineage[taxon] = levels[0:]
                levels.append(taxon)
                taxtuple[taxon[1]] = taxon
                self.bugPrint('TaxTree::{0}\n'.format(taxtree))
                self.debug('Levels::{0}\n'.format(levels))
                #i# Add tabdb entry
                if rentry['taxnum'] < mintaxnum: continue
                if rentry['taxrank'] not in taxhead: continue
                if not self.getBool('TaxPart') and rentry['taxrank'] != self.getStrLC('TaxLow'): continue
                tentry = {'taxid':rentry['taxid'],'taxnum':rentry['taxnum']}
                for taxon in levels:
                    if taxon[0] in taxhead: tentry[taxon[0]] = taxon[2]
                tabdb.addEntry(tentry)
            self.printLog('\r#REPORT', 'Parsing taxonomy complete')
            ## ~ [2b] Save tabreport table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tabdb.saveToFile('{0}.tabreport.tsv'.format(taxbase))
            return True
        except: self.errorLog('%s.tabReport error' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: Taxolotl Class                                                                                   #
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
    try: Taxolotl(mainlog,cmd_list).run()

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
