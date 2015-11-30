#!/usr/bin/python

# rje_embl.py - RJE Module to Handle EMBL nucleotide Files
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <software@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_embl
Description:  RJE Module to Handle EMBL Nucleotide Files
Version:      0.1
Last Edit:    20/11/08
Copyright (C) 2008 Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains methods for handling EMBL Nucleotide files, primarily in other rje modules but also with some
    standalone functionality. It is built on rje_uniprot and a lot of the initial functionality is designed around
    converting EMBL files to be compatible with programs designed with UniProt files in mind.

    This module can be used to extract a list of entries from a larger database and/or to produce summary tables from
    flat files.

Input Options:
    embldat=FILE    : Name of EMBL file [None]
    emblpath=PATH   : Path to EMBL Datafile (will look here for DB Index file)
    dbindex=FILE    : Database index file [uniprot.index]
    extract=LIST    : Extract IDs/AccNums in list. LIST can be FILE or list of IDs/AccNums X,Y,.. []
    acclist=LIST    : As Extract.
    specdat=LIST    : Make a DAT file of the listed species from the index (over-rules extract=LIST) []
    
Output Options:
    makeindex=T/F   : Generate UniProt index files [False]
    makespec=T/F    : Generate species table [False]
    makefas=T/F     : Generate fasta files [False]
    datout=FILE     : Name of new (reduced) DAT file of extracted sequences [None]
    outformat=X     : Format for generated DAT file (uniprot/embl) [embl]
    splitout=PATH   : If path given, will split output into individual files per entry into PATH []
    tabout=FILE     : Table of extracted details [None]
    linkout=FILE    : Table of extracted Database links [None]
    ftout=FILE      : Outputs table of features into FILE [None]

Format Conversion Options:
    ucft=X          : Feature to add for UpperCase portions of sequence []
    lcft=X          : Feature to add for LowerCase portions of sequence []
    maskft=LIST     : List of Features to mask out []
    invmask=T/F     : Whether to invert the masking and only retain maskft features [False]
    caseft=LIST     : List of Features to make upper case with rest of sequence lower case []
    ftskip=LIST     : List of feature details to skip ['transl_table','translation']

General Options:
    append=T/F      : Append to results files rather than overwrite [False]
    memsaver=T/F    : Memsaver option to save memory usage - does not retain entries in object [False]    

Uses general modules: glob, os, re, string, sys, time
Uses RJE modules: rje, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_sequence, rje_uniprot
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation modified from rje_uniprot.
    # 0.1 - Added type storage and file splitting during conversion.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Lots of functionality to add! Look also to BioPython.
    # [ ] : Add counts of feature types for output
    # [ ] : Get basic version working
    # [ ] : Update processDAT() and rje_uniprot.processUniProt() methods
    # [ ] : Add output format options - Fasta/EMBL/UniProt
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_EMBL', '0.1', 'November 2008', '2008')
    description = 'RJE EMBL Nucleotide Parsing/Extraction Module'
    author = 'Dr Richard J. Edwards.'
    comments = []
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
### SECTION II: MODULE CONSTANTS                                                                                        #
#########################################################################################################################
### UniProt Parsing dictionary: Add crucial information to parse out here. Used by UniProtEntry.process()   ###
emblparse = {
    'ID' : string.join(['^(\S+);*','(\S.+);','(\S+);\s+(\d+)\s+BP'], '\s+'),       # ID, Extra Info, Length
    'AC' : '(\S+);',            # Primary Acc
    'DE' : '\s*(\S.+)',         # Description
    'OS' : '^(\S.+)\s*$',       # Species
    'RX' : 'PubMed=(\d+);',     # PubMed ID
    'RC' : 'TISSUE=(.+);',      # Tissue(s)                             ??? Not in EMBL ???
    'DR' : '^(\S+);\s+(\S.+)$', # Database links (Dbase,Details)
    'CC' : '^-!-\s+(\S.+):\s+(\S.+)$',  # Comments (Type, Details)  ??? Not in EMBL ???
    'FT' : '(\S+)\s+<*(\d+)\.\.(\d+)',    # Feature type, start, stop
    'FCOMP' : '(\S+)\s+complement\((\d+)\.\.(\d+)\)',    # Feature type, start, stop
    'FJOIN' : '(\S+)\s+join\((\d+)\..*\.(\d+)\)',    # Feature type, start, stop
    'FJCOM' : '(\S+)\s+join\(complement\((\d+)\..*\.(\d+)\)\)',    # Feature type, start, stop
    'FNOTE' : '/(\S+)=(\S.*)$'             # Feature details (Detail, Information)
    }
#########################################################################################################################
useful_data = ['ID','AC','DE','GN','OS','OC','RX','CC','DR','RC','KW','FT']     # Data to retain following parsing #
#########################################################################################################################
featurelist = ['LIPID','TRANSMEM','MOD_RES','DOMAIN']   #!# Features for function table. Add more! #!#
ftskip = ['transl_table','translation','codon_start','citation']
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: EMBL Class                                                                                             # 
#########################################################################################################################
class EMBL(rje_uniprot.UniProt):     
    '''
    EMBL Download Class. Author: Rich Edwards (2008).

    Info:str
    - Name = Name of EMBL File 
    - UniPath = Path to EMBL Datafile (will look here for DB Index file) [EMBL/]
    - DBIndex = Database index file [embl.index]
    - DATOut = Name of new (reduced) UniProt DAT file of extracted sequences [None]
    - TabOut = Name of table of extracted sequence details [None]
    - LinkOut = Table of extracted Database links [None]
    - FTOut = Outputs table of features into FILE [None]
    - OutFormat = Format for generated DAT file (uniprot/embl) [embl]
    - SplitOut = If path given, will split output into individual files per entry into PATH []
    - UCFT = Feature to add for UpperCase portions of sequence []
    - LCFT = Feature to add for LowerCase portions of sequence []
    
    Opt:boolean
    - MakeIndex = Generate EMBL index files [False]
    - MakeSpec = Generate species table [False]
    - MakeFas = Generate fasta files [False]

    Stat:numeric

    List:list
    - Entry = list of EMBLEntry objects
    - Extract = Extract AccNums/IDs in list. LIST can be FILE or list of AccNums X,Y,.. []
    - SpecDat = Make a DAT file of the listed species from the index []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#####################4####################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','UniPath','DBIndex','DATOut','TabOut','LinkOut','FTOut','UCFT','LCFT','OutFormat','SplitOut']
        self.optlist = ['MakeIndex','MakeSpec','MakeFas']
        self.statlist = []
        self.listlist = ['Extract','Entry','SpecDat']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'UniPath':rje.makePath('EMBL/'),'DBIndex':'embl.index','UCFT':'','LCFT':'','OutFormat':'embl'})
        self.opt['SpliceVar'] = False   # Needed for some inherited rje_uniprot functions
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdRead(cmd,'path',att='UniPath',arg='emblpath')
                self._cmdReadList(cmd,'file',['DBIndex','DATOut','TabOut','LinkOut','FTOut'])
                self._cmdReadList(cmd,'path',['SplitOut'])
                self._cmdReadList(cmd,'info',['UCFT','LCFT','OutFormat'])
                self._cmdReadList(cmd,'opt',['MakeIndex','MakeSpec','MakeFas'])
                self._cmdReadList(cmd,'list',['Extract','SpecDat'])
                self._cmdRead(cmd,type='file',att='Name',arg='embldat')
                self._cmdRead(cmd,type='list',att='Extract',arg='acclist')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def run(self):  ### Main Run Method if called direct from commandline. Returns self or None.
        '''Main Run Method if called direct from commandline. Returns True if no Errors, else False.'''
        try:### ~ [1] ~ Higher Level Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Make Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['MakeIndex'] or self.opt['MakeSpec'] or self.opt['MakeFas']:
                processDAT(self,makeindex=self.opt['MakeIndex'],makespec=self.opt['MakeSpec'],makefas=self.opt['MakeFas'])
                return
            ## ~ [1b] ~ Get Extract List from species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['SpecDat']: self.extractSpecies()
                
            ### ~ [2] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Check for input needs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['Name'].lower() in ['','none'] and not self.list['Extract']:   # No AccNum!
                self.errorLog('No input file or acclist given. Use "help" option for parameters.',printerror=False)
                return False         
            ## ~ [2b] ~ Setup Output Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for file in ['DATOut','TabOut','LinkOut','FTOut']:
                if self.info[file].lower() == 'none': self.info[file] = ''
                if self.info[file]: rje.backup(self,self.info[file])
            ## ~ [2c] ~ Extracted details & MemSaver ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (self.info['TabOut'] or self.info['LinkOut'] or self.info['FTOut']) and self.opt['MemSaver']:
                memtext = 'TabOut, LinkOut and FTOut will not function with MemSaver mode.'
                if self.stat['Interactive'] >= 0 and rje.yesNo('%s Switch Memsaver off?' % memtext):
                    self.opt['MemSaver'] = False
                    self.cmd_list += ['memsaver=F']
                else: self.printLog('#MEM',memtext)
            
            ### ~ [3] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Read UniProt File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.readUniProt(reformat=self.info['OutFormat'].lower() in ['uniprot'])
            ## ~ [3b] ~ Special Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.tableOutput() 
            if self.info['FTOut'] and not self.opt['MemSaver']: self.ftTable(self.info['FTOut'])
            return self
        except:
            self.errorLog('Fundamental error during EMBL.run()')
            return None
#########################################################################################################################
    ### <2> ### Reading DAT Entries                                                                                     #
#########################################################################################################################
    def _newEntryObject(self): return EMBLEntry(log=self.log,cmd_list=self.cmd_list)
#########################################################################################################################
    def addFromSeq(self,seq=None,sequence='',name='',data={},ft=[]):    ### Converts into UniProtEntry object 
        '''
        Converts into UniProtEntry object and adds to self.
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns entry if successful or None if fails
        '''
        ### Add case features ###
        addft = ft[0:]
        if self.info['UCFT'] and seq:
            for uc in seq.dict['Case']['Upper']:
                addft.append({'Type':self.info['UCFT'].upper(),'Desc':self.info['UCFT'],'Start':uc[0]+1,'End':uc[1]+1})
        if self.info['LCFT'] and seq:
            for lc in seq.dict['Case']['Lower']:
                addft.append({'Type':self.info['LCFT'].upper(),'Desc':self.info['LCFT'],'Start':lc[0]+1,'End':lc[1]+1})
        ### Make Entry ###        
        newentry = UniProtEntry(self.log,self.cmd_list)
        entry = newentry.uniProtFromSeq(seq,sequence,name,data,addft)
        if entry: self.list['Entry'].append(entry)
        return entry
#########################################################################################################################
    ### <3> ### UniProt Info Output
#########################################################################################################################
    def tableOutput(self):   ### Tabulated output of UniProt information
        '''Tabulated output of UniProt information. Divided into TabOut (UniProt summary) and LinkOut (database links)'''
        try:

            #!# Needs updating for EMBL #!#
            
            ### Database Links Table ###
            if self.info['LinkOut'].lower() not in ['','none']:
                self.linkOutput()
            if self.info['TabOut'].lower() in ['','none']:
                return

            ### Setup Groups and Columns ###
            headers = {'#1# Basic Details #':['#','AccNum','ID','Len','Description','Gene','Species'],
                       '#2# Function & Activity #':['Function','GO_MF','GO_BP','Activity','Interactions','Phenotype',
                                                    'Similarity'],
                       '#3# Expression, Location & Structure #':['Tissue','Cell_Loc','PDB','InterPro','Pfam','PROSITE',
                                                                 'Isoforms','Ensembl'],
                       '#4# References & Links #':['GeneCards','PubMed','Keywords','Comments']
                       }

            ### Open File and write header ###
            delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['TabOut']))
            TABOUT = open(self.info['TabOut'],'a')   # Already deleted if append=F
            TABOUT.write('# Generated by %s: %s\n' % (self.log.info['Name'],time.asctime(time.localtime(time.time()))))
            if self.info['Name'] != 'None':
                TABOUT.write('# Source: %s\n' % os.path.abspath(self.info['Name']))
            else:
                TABOUT.write('#Source: %s\n' % os.path.abspath(self.info['UniPath']+'*.dat'))
            TABOUT.write('# Seqnum: %s\n\n' % rje.integerString(self.entryNum()))
            ## Headers ##
            head1 = []
            head2 = []
            for h in rje.sortKeys(headers):
                head1.append(h)
                head1 += [''] * (len(headers[h])-1)
                head2 += headers[h]
            rje.writeDelimit(TABOUT,head1,delimit)
            rje.writeDelimit(TABOUT,head2,delimit)

            ### Write Data for Entries ###
            ex = 0
            for entry in self.list['Entry']:
                seq = entry.obj['Sequence']
                ex += 1
                data = []
                comments = rje.sortKeys(entry.dict['Comments'])
                for h in head2:     # Column headers:
                    #1# Basic Details #
                    if h == '#':
                        data.append(rje.preZero(ex,self.entryNum()))
                    elif h in ['AccNum','ID','Description']:
                        data.append(seq.info[h])
                    elif h == 'Len':
                        data.append('%d' % seq.aaLen())
                    elif h == 'Gene':
                        data.append(string.join([seq.info['Gene']] + entry.list['Synonyms'],'; '))
                    elif h == 'Species':
                        data.append('%s [%s]' % (seq.info['Species'],seq.info['SpecCode']))
                    #2# Function & Activity #
                    elif h == 'Function':
                        text = ''
                        for cc in ['FUNCTION','PATHWAY','DOMAIN']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'GO_MF':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; F:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'GO_BP':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; P:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'Activity':   #!# Join to Function #!#
                        text = ''
                        for cc in ['CATALYTIC ACTIVITY','COFACTOR']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Interactions':
                        text = ''
                        for cc in ['INTERACTION','ENZYME REGULATION','SUBUNIT','PTM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Phenotype':
                        text = ''
                        for cc in ['DISEASE','POLYMORPHISM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Similarity':
                        text = ''
                        for cc in ['SIMILARITY']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #3# Expression, Location & Structure #
                    elif h == 'Tissue':
                        text = ''
                        if entry.list['Tissues']:
                            text = 'TISSUES: %s' % string.join(entry.list['Tissues']+[''],'; ')
                        for cc in ['TISSUE SPECIFICITY','DEVELOPMENTAL STAGE','INDUCTION']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Cell_Loc':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; C:') < 0:
                                    go.remove(g)
                        cc = 'SUBCELLULAR LOCATION'
                        if entry.dict['Comments'].has_key(cc):
                            comments.remove(cc)
                            go = entry.dict['Comments'][cc] + go
                        data.append(string.join(go,' >> '))
                    elif h in ['PDB','InterPro','Pfam','PROSITE','Ensembl']:
                        if entry.dict['DBLinks'].has_key(h):
                            data.append(string.join(entry.dict['DBLinks'][h],' >> '))
                        else:
                            data.append('')
                    elif h == 'Isoforms':
                        text = ''
                        for cc in ['ALTERNATIVE PRODUCTS']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #4# References & Links
                    elif h == 'GeneCards':
                        if entry.dict['DBLinks'].has_key('HGNC'):
                            data.append(string.join(entry.dict['DBLinks']['HGNC'],' >> '))
                        else:
                            data.append('')
                    elif h in ['PubMed']:
                        if entry.list[h]:
                            data.append('http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=%s' % string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h in ['Keywords']:
                        if entry.list[h]:
                            data.append(string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h == 'Comments':
                        text = ''
                        for cc in comments:
                            if text:
                                text += ' '
                            text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                            if text[-1] != '.':
                                text += '.'
                        data.append(text)

                rje.writeDelimit(TABOUT,data,delimit)
                self.log.printLog('\r#OUT','UniProt summary output: %.1f%%' % (100.0 * ex / self.entryNum()),log=False,newline=False)
            self.log.printLog('\r#OUT','UniProt summary output for %s entries.' % rje.integerString(self.entryNum()))
            TABOUT.close()


        except:
            self.log.errorLog('Major error during UniProt.tableOutput()!',quitchoice=True)            
#########################################################################################################################
    def ftTable(self,outfile):  ### Outputs features into a table
        '''
        Outputs features into a table.
        >> outfile:str = Name of output file
        '''
        try:

            #!# Needs updating for EMBL #!#
            
            ### Setup ###
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            if self.opt['Append']:
                FT = open(outfile,'a')
            else:
                FT = open(outfile,'w')
                rje.writeDelimit(FT,['acc_num','feature','ft_start','ft_end','description'],delimit)

            ### Output ###
            (fx,ex) = (0,0.0)
            for entry in self.list['Entry']:
                ex += 100.0
                acc = entry.obj['Sequence'].info['AccNum']
                ## Make dictionary of {start:{end:[features]}}
                ft_dict = {}
                for ft in entry.list['Feature']:
                    ft_start = ft['Start']
                    if not ft_dict.has_key(ft_start):
                        ft_dict[ft_start] = {}
                    ft_end = ft['End']
                    if not ft_dict[ft_start].has_key(ft_end):
                        ft_dict[ft_start][ft_end] = []
                    ft_dict[ft_start][ft_end].append(ft)
                ## Sort and output ##
                for ft_start in rje.sortKeys(ft_dict):
                    for ft_end in rje.sortKeys(ft_dict[ft_start]):
                        for ft in ft_dict[ft_start][ft_end]:
                            outlist = [acc]
                            for fk in ['Type','Start','End','Desc']:
                                outlist.append('%s' % ft[fk])
                            rje.writeDelimit(FT,outlist,delimit)
                        fx += 1
                self.log.printLog('\r#FT','Feature output: %.1f%% (%s features)' % (ex/len(self.list['Entry']),rje.integerString(fx)),log=False,newline=False)

            ### End ###
            FT.close()
            self.log.printLog('\r#FT','Feature output complete: %s features, %s entries.' % (rje.integerString(fx),rje.integerString(len(self.list['Entry']))))
        except:
            self.log.errorLog('Program error during rje_uniprot.ftTable()',quitchoice=True)
#########################################################################################################################
    def saveEMBL(self,outfile,entries=[],append=False):    ### Saves self as a DAT file
        '''
        Saves self as a DAT file.
        >> outfile:str = Name of output file
        >> entries:list of entries (self.list['Entry'] if none given)
        >> append:boolean = whether to append file
        '''
        try:

            #!# Modify this to save in EMBL format #!#
            
            ### Setup Objects ###
            if not append: rje.backup(self,outfile)
            if not entries: entries = self.list['Entry'][0:]

            ### Output ###
            OUT = open(outfile,'a')                    
            for entry in entries[0:]:
                if not entry.dict['Data'] and not entry.uniProtFromSeq():
                    entries.remove(entry)
                    self.log.errorLog('Problem with %s (%s) - cannot output' % entry.info['Name'],printerror=False)
                    continue
                seq = entry.obj['Sequence']
                ## Standard info ##
                for key in ['ID','AC','DT','DE','GN','OS']:
                    if entry.dict['Data'].has_key(key):
                        for rest in entry.dict['Data'][key]:
                            OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## Other data, except Features and sequence ##
                for key in rje.sortKeys(entry.dict['Data']):
                    if key not in ['ID','AC','DT','DE','GN','OS','FT','SQ','SEQ','//']:
                        for rest in entry.dict['Data'][key]:
                            OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## Features ##
                entry.orderFT()
                for ftdict in entry.list['Feature']:
                    (p1,p2) = (ftdict['Start'],ftdict['End'])
                    ftxt = 'FT   %s' % ftdict['Type']
                    while len(ftxt) < 14 or ftxt[-1] != ' ': ftxt += ' '
                    ftxt += '%6s' % ('%d' % p1)
                    while len(ftxt) > 20 and ftxt[-(len('%d' % p1)+2):-len('%d' % p1)] == '  ': ftxt = ftxt[:-(len('%d' % p1)+1)] + ftxt[-len('%d' % p1):]
                    ftxt += '%7s' % ('%d' % p2)
                    while len(ftxt) > 27 and ftxt[-(len('%d' % p2)+2):-len('%d' % p2)] == '  ': ftxt = ftxt[:-(len('%d' % p2)+1)] + ftxt[-len('%d' % p2):]
                    ftxt += ' %s\n' % ftdict['Desc']
                    OUT.write(ftxt)
                ## Sequence/End ##
                OUT.write('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen(),rje_sequence.MWt(seq.info['Sequence'])))
                uniseq = seq.info['Sequence'][0:]
                while len(uniseq) > 0:
                    OUT.write('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                    uniseq = uniseq[60:]
                OUT.write('//\n')
            OUT.close()

            self.log.printLog('#OUT','UniProt format for %d entries saved to %s' % (len(entries),outfile))            
        except:
            self.log.errorLog('Major problem with UniProt.saveUniProt()')
#########################################################################################################################
    ### <4> ### Additional UniProt Tools                                                                                #
#########################################################################################################################
    #!# Currently sticking with rje_uniprot versions until problems arise:
    # - accNameSeq(self,acc_list=[],spec=None,justsequence=True):  ### Method to extract dictionaries of {acc:'ID__PrimaryAcc Desc'} & {acc:seq}
    # - accDict(self,acc_list=[],cleardata=True):      ### Method to extract dictionaries of {acc:UniProtEntry}
#########################################################################################################################
## End of SECTION III: EMBL Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: EMBLEntry Class                                                                                         # 
#########################################################################################################################
class EMBLEntry(rje_uniprot.UniProtEntry):     
    '''
    EMBL Entry Class. Author: Rich Edwards (2008).

    Info:str
    - Name = UniProt ID of Entry
    - Type = Preliminary or Standard
    - FullText = Full Text of Entry
    
    Opt:boolean
    - InvMask = Whether to invert the masking and only retain maskft features [False]    

    Stat:numeric
    - Length = Length of Sequence as annotated

    List:list
    - CaseFT = List of Features to make upper case with rest of sequence lower case []
    - Feature = List of feature dictionaries: [Type,Start,End,Desc]
    - FTSkip = List of feature details to skip ['transl_table','translation']
    - MaskFT = List of Features to mask out []
    - PubMed = List of PubMed IDs (as strings)
    - Keywords = List of UniProt Keywords
    - Tissues = List of UniProt Tissues
    - Synonyms = List of Gene synonyms
    
    Dict:dictionary
    - Data = Dictionary of lists of UniProt data (Keys are line headers ID/AC/CC etc.)
    - DB = Specific extractions from DR lines for use in other programs. {DB:[AccNum/ID]}
    - Comments = Dictionary of comments: {Type:List of Comments}
    - DBLinks = List of Database Link dictionaries {Dbase,List of Details} for dblinks output

    Obj:RJE_Objects
    - Sequence = rje_sequence.Sequence object
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','FullText']
        self.statlist = ['Length']
        self.optlist = ['CC2FT','InvMask','TMConvert']
        self.listlist = ['Feature','PubMed','Keywords','Tissues','Synonyms','MaskFT','CaseFT','FTSkip']
        self.dictlist = ['Data','Comments','DBLinks','DB']
        self.objlist = ['Sequence']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.obj['Sequence'] = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list)
        self.info['FullText'] = ''
        self.list['Feature'] = []   # List of features = {'Type':str,'Start':int,'End':int,'Desc':str}
        self.list['MaskFT'] = []
        self.list['FTSkip'] = ftskip[0:]
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'opt',['CC2FT','InvMask','TMConvert'])
                self._cmdReadList(cmd,'list',['CaseFT','MaskFT','FTSkip'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Attribute Processing
#########################################################################################################################
    def _uniParse(self,key):    ### Parses list of elements from self.dict['Data'] (and pops) 
        '''
        Parses list of elements from self.dict['Data'] (and pops).
        >> key:str = Key of UniProt entry type
        << List of matched elements or False if failure.
        '''
        try:### ~ [1] Check key and parse data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if key not in self.dict['Data'].keys(): return None
            return rje.matchExp(emblparse[key],self.dict['Data'][key][0])
        except: self.errorLog('EMBLEntry: Cataclysmic error during _uniParse(%s)!' % key)
        return None
#########################################################################################################################
    def process(self,logft=True,cleardata=True):  ### Extract Details from self.dict['Data'] to Sequence object
        '''
        Extract Details from self.dict['Data'] to Sequence object.
        >> logft:boolean = whether to write number of features to log file
        >> cleardata=Whether to clear self.dict['Data'] after processing (to save memory) [True]
        << True if OK, False if not.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Setup'
            seq_info = ['Name','Type','Description','Sequence','ID','AccNum','DBase','Gene','Species','SpecCode','Format']
            seqi = self.obj['Sequence'].info
            seqi['Type'] = 'DNA'    #!# RNA? #!#
            for key in ['DE']:
                if self.dict['Data'].has_key(key): self.dict['Data'][key] = [string.join(self.dict['Data'][key],' ')]
            for key in ['FH','RT','RL','RA','RN','RP']:
                try: self.dict['Data'].pop(key)
                except: pass
                                    
            ### ~ [1] ~ Basic Sequence Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            parse = self._uniParse('OS')
            if parse:
                seqi['Species'] = parse[0]
                seqi['SpecCode'] = rje_sequence.getSpecCode(seqi['Species'])
            else:
                seqi['Species'] = 'Unknown'
                seqi['SpecCode'] = 'UNK'
            ## ~ [1b] ~ Sequence ID (ID) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            parse = self._uniParse('ID')
            if parse: (self.info['Name'],self.info['TypeExtra'],self.info['Type'],self.stat['Length']) = parse
            if self.info['Name'][-1:] == ';': self.info['Name'] = self.info['Name'][:-1]
            self.stat['Length'] = string.atoi(self.stat['Length'])
            self.info['ID'] = seqi['ID'] = '%s_%s' % (self.info['Name'],seqi['SpecCode'])
            self.dict['Data']['ID'][0] = string.join([self.info['ID']]+string.split(self.dict['Data']['ID'][0])[1:])
            ## ~ [1c] ~ AccNum (AC) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            full_acc = string.split(string.join(self.dict['Data']['AC']))
            self.obj['Sequence'].list['Secondary ID'] = string.split(string.join(full_acc,''),';')[1:-1]
            parse = self._uniParse('AC')
            if parse: seqi['AccNum'] = parse[0]
            seqi['DBase'] = 'embl'
            seqi['Name'] = self.info['Name'] = '%s__%s' % (seqi['ID'],seqi['AccNum'])
            seqi['Format'] = 'gn_sp__acc'
            ## ~ [1d] ~ Description (DE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            parse = self._uniParse('DE')
            if parse: seqi['Description'] = parse[0]
            ## ~ [1e] ~ Sequence (SEQ) - initial descriptive line in ['SQ'] ~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('SEQ'): seqi['Sequence'] = self.dict['Data'].pop('SEQ')[0]
            seqi['Sequence'] = re.sub('\s+','',seqi['Sequence'])
            seqi['Sequence'] = re.sub('\d+','',seqi['Sequence']).upper()

            ### ~ [2] ~ Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dict['Data'].has_key('FT'):
                ftlist = []
                for ft in self.dict['Data']['FT']:
                    ftsplit = string.split(ft,'/')
                    ftlist.append(ftsplit[0])
                    for eft in ftsplit[1:]: ftlist.append('/%s' % eft)
                self.deBug(ftlist)
                extendft = False
                for ft in ftlist:   #x#self.dict['Data']['FT']:
                    parse = rje.matchExp(emblparse['FT'],ft)
                    parse_comp = rje.matchExp(emblparse['FCOMP'],ft)
                    if parse_comp: parse = parse_comp
                    parse_join = rje.matchExp(emblparse['FJOIN'],ft)
                    parse_jcom = rje.matchExp(emblparse['FJCOM'],ft)
                    parse_embl = rje.matchExp(emblparse['FNOTE'],ft)
                    parse_cntd = rje.matchExp('^\s*(\S.*)$',ft)
                    #x#self.deBug('%s\n%s\n%s\n%s' % (ft,parse,parse_embl,parse_cntd))
                    if parse:
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : '',
                            'EMBL' : {}
                            }
                        if parse_comp: ftdic['Desc'] = 'complement'
                        self.list['Feature'].append(ftdic)
                    elif parse_join:
                        ftdic = {
                            'Type' : parse_join[0],
                            'Start' : string.atoi(parse_join[1]),
                            'End' : string.atoi(parse_join[2]),
                            'Desc' : string.split(ft)[1],
                            'EMBL' : {}
                            }
                        if parse_comp: ftdic['Desc'] = 'complement'
                        self.list['Feature'].append(ftdic)
                    elif parse_jcom:
                        ftdic = {
                            'Type' : parse_jcom[0],
                            'Start' : string.atoi(parse_jcom[2]),
                            'End' : string.atoi(parse_jcom[1]),
                            'Desc' : string.split(ft)[1],
                            'EMBL' : {}
                            }
                        if parse_comp: ftdic['Desc'] = 'complement'
                        self.list['Feature'].append(ftdic)
                    else:
                        try:
                            ftdic = self.list['Feature'][-1]
                            if parse_embl:
                                (type,details) = parse_embl
                                if type.lower() in self.list['FTSkip']: extendft = False; continue
                                extendft = True
                                if ftdic['Desc']: ftdic['Desc'] = '%s ' % ftdic['Desc']
                                ftdic['Desc'] = re.sub('\s+',' ','%s/%s=%s' % (ftdic['Desc'],type,details))
                                if type in ftdic['EMBL']: ftdic['EMBL'][type] += details
                                else: ftdic['EMBL'][type] = details
                            elif not extendft: continue
                            elif parse_cntd:
                                ftdic['Desc'] = re.sub('\s+',' ','%s %s' % (ftdic['Desc'],parse_cntd[0]))
                            else: raise ValueError
                        except: self.errorLog('Cannot parse feature details from "%s"' % ft)

            ### ~ [3] ~ Other Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Tissues (RC) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('RC'):
                tissues = []
                self.list['Tissues'] = []
                for rc in self.dict['Data']['RC']:
                    parse = rje.matchExp(emblparse['RC'],rc)
                    if parse: tissues += string.split(parse[0],', ')
                for tissue in tissues:
                    if tissue[:4] == 'and ': self.list['Tissues'].append(tissue[4:])
                    else: self.list['Tissues'].append(tissue)
            ## ~ [3b] ~ Keywords (KW) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('KW'):
                keywords = string.join(self.dict['Data']['KW'])
                if keywords[-1:] == '.': keywords = keywords[:-1]    # Remove full stop
                self.list['Keywords'] = string.split(keywords,'; ')
            ## ~ [3c] ~ References (RX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('RX'):
                self.list['PubMed'] = []
                for rx in self.dict['Data']['RX']:
                    parse = rje.matchExp(emblparse['RX'],rx)
                    if parse: self.list['PubMed'].append(parse[0])
            ## ~ [3d] ~ Comments (CC) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('CC'):
                self.dict['Comments'] = {}  # Dictionary of comments: {Type:List of Comments}
                for cc in self.dict['Data']['CC']:
                    if cc.find('-----') == 0: break
                    csplit = string.split(cc[4:],': ')
                    ctype = csplit[0]
                    cdetail = string.join(csplit[1:],': ')
                    if self.dict['Comments'].has_key(ctype): self.dict['Comments'][ctype].append(cdetail)
                    else: self.dict['Comments'][ctype] = [cdetail]
            ## ~ [3e] ~ Database Links (DR) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dict['Data'].has_key('DR'):
                self.dict['DBLinks'] = {}   # Database Link dictionary {Dbase,List of Details}
                for dr in self.dict['Data']['DR']:
                    if rje.matchExp(emblparse['DR'],dr):
                        (ctype,cdetail) = rje.matchExp(emblparse['DR'],dr)
                        if self.dict['DBLinks'].has_key(ctype): self.dict['DBLinks'][ctype].append(cdetail)
                        else: self.dict['DBLinks'][ctype] = [cdetail]
                        self.specialDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']
                    elif rje.matchExp('^(Entrez Gene); (\S+)$',dr):
                        (ctype,cdetail) = rje.matchExp('^(Entrez Gene); (\S+)$',dr)
                        self.specialDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']

            ### ~ [4] ~ FT Masking/Case Change ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['MaskFT']: self.maskFT(self.list['MaskFT'],inverse=self.opt['InvMask'])
            if self.list['CaseFT']: self.caseFT(self.list['CaseFT'])

            ### ~ [5] ~ Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata: self.dict['Data'] = {}    # Save memory!
            if logft: self.printLog('#FT','%d features for %s.' % (len(self.list['Feature']),self.obj['Sequence'].info['AccNum']))
            return True
        except: self.errorLog('Cataclysmic error during EMBLEntry.process()!'); return False
#########################################################################################################################
    ### <3> ### UniProt Conversion and Saving                                                                           #
#########################################################################################################################
    def uniProtFromSeq(self,seq=None,sequence='',name='',data={},ft=[]): ### Converts into UniProtEntry object (self!)
        '''
        Converts into UniProtEntry object (self!).
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns self if successful or None if fails
        '''
        try:
            #!# Update to EMBL method #!#
            ### Setup ###
            if not seq and (sequence and name):
                seq = rje_sequence.Sequence(log=self.log)
                seq.info['Name'] = name
                seq.info['Sequence'] = sequence.upper()  #!# Change at some point to allow mixed case!
                seq.info['Type'] = 'Protein'
                seq.extractDetails()    #!#gnspacc=self.opt['GeneSpAcc'])
            if not seq:
                seq = self.obj['Sequence']
            if not seq:
                raise ValueError, 'No sequence information given'
            self.obj['Sequence'] = seq

            ### Update self ###
            for key in data:
                self.dict['Data'][key] = data[key]
            self.list['Feature'] += ft
            if seq.info['DBase'] != 'trembl':
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['ID'],seq.info['Type'],seq.aaLen())]
            elif seq.info['SpecCode'] not in ['None','UNK']:
                self.dict['Data']['ID'] = ['%s_%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['SpecCode'],seq.info['Type'],seq.aaLen())]
            else:
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['Type'],seq.aaLen())]
            if self.dict['Data'].has_key('AC'):
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']] + self.dict['Data']['AC']
            else:
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']]
            if seq.info['Description'].lower() != 'none':
                self.dict['Data']['DE'] = [seq.info['Description']]
            else:
                self.dict['Data']['DE'] = ['']
            if seq.info['Species'] not in ['None','Unknown',seq.info['SpecCode'],'']:
                if self.dict['Data'].has_key('OS'):
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']] + self.dict['Data']['OS']
                else:
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']]
            dt = string.split(time.ctime())
            if self.dict['Data'].has_key('DT'):
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])] + self.dict['Data']['DT']
            else:
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])]
            if self.dict['Data'].has_key('CC'):
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()] + self.dict['Data']['CC']
            else:
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()]

            #X#self.deBug(self.dict['Data'])
            if self.process(logft=False,cleardata=False):
                return self
            return None
        except:
            self.log.errorLog('UniProtEntry.uniProtFromSeq() has gone wrong.',quitchoice=True)
            return None
#########################################################################################################################
## End of SECTION III : EMBLEntry Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: General UniProt Methods                                                                                 #
#########################################################################################################################
def processDAT(callobj,makeindex=True,makespec=True,makefas=True):  ### Processes DAT file making index file and spectable as appropriate
    '''Processes DAT file making index file and spectable as appropriate.'''  
    try:### ~ [1] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #!-- indexfile = unipath + 'uniprot.index'
        #!# Add additonal basefile argument (and dbindex reading) to change 'uniprot' in method
        return rje_uniprot.processUniProt(callobj,makeindex,makespec,makefas)
    except: callobj.errorLog('UniProt/EMBL processing incompatability Error',printerror=True)
#########################################################################################################################
## End of SECTION IV : Generic Module Methods                                                                           #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MAIN PROGRAM                                                                                             #
#########################################################################################################################
def runMain():
    ### ~ [0] ~ Basic Setup of Program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
    ### ~ [1] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:EMBL(mainlog,cmd_list).run()
        
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
### END OF SECTION V                                                                                                    #
#########################################################################################################################

