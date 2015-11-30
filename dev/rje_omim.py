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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       RJE_OMIM
Description:  OMIM Parsing module
Version:      0.1
Last Edit:    13/08/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for parsing and storing OMIM gene/mutation data. In this capacitity it may be used by other modules.
    It also has a standalone capability specific to the OMIM-SLiMFinder grant application of Dr Richard Edwards. In this
    capacity, it should not be used without the permission of the author.

Commandline:
    omim=FILE       : Input file (OMIM text download) ['omim.txt']
    slimdir=PATH    : Directory where *.occ.csv files can be found []

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Incorporation of pingu and slimfinder etc. in standalone run mode.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Parsing mutations from omim.txt OMIM download.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_OMIM', '0.1', 'August 2007', '2007')
    description = 'OMIM Parsing Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please attain author\'s permission to run as standalone application',rje_zen.Zen().wisdom()]
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
### SECTION II: OMIM Class                                                                                              #
#########################################################################################################################
class OMIM(rje.RJE_Object):     
    '''
    OMIM Class. Author: Rich Edwards (2007).

    Info:str
    - Name = Input file (OMIM text download) ['omim.txt']
    - SlimDir = Directory where *.occ.csv files can be found []
    
    Opt:boolean

    Stat:numeric

    List:list

    Dict:dictionary
    - Fudge = Dictionary of {Gene:Fudge factor for mutations vs EnsLoci sequence}
    - Mutations = Dictionary of {Gene:{SubID:(Disease,Mutation)}}
    - Records = Dictionary of {Gene:[OMIM IDs]}

    Obj:RJE_Objects
    - Pingu = PINGU object storing PPI data and mapping genes to sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','SlimDir']
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = ['Mutations','Records']
        self.objlist = ['Pingu']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Name':'omim.txt'})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._cmdRead(cmd,type='file',att='Name',arg='omim')
                self._cmdReadList(cmd,'path',['SlimDir'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main OMIM Run Methods                                                                                   #
#########################################################################################################################
    def run(self):  ### Main Run Method
        '''Main Run Method.'''
        try:### ~ [1] Parse/Read Mutation data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Force'] or not self.loadMutations(): self.parseOMIM()

            ### ~ [2] Additional Pingu incorporation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Load PPI data using Pingu, map genes to sequences and check mutation residues #!#
            ## ~ [2a] Setup Pingu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            import pingu
            pcmd = self.cmd_list + ['fulloutput=F']
            ping = self.obj['Pingu'] = pingu.PINGU(self.log,pcmd)
            ping.run()
            ## ~ [2b] Read in EnsLoci sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not ping.obj['GeneCards']: return self.log.errorLog('Cannot map EnsLoci without GeneCards.', printerror=False)
            genecards = ping.obj['GeneCards'].dict['GeneCard']      # GeneCards dictionary
            ensloci = ping.getEnsLoci()     # EnsLoci SeqList object (ping.obj['EnsLoci'])
            seqdict = ensloci.seqNameDic()  
            if not seqdict: return self.log.errorLog('Failed to read in EnsLoci sequences.', printerror=False)
            ## ~ [2c] Calculate fudge factor for each gene ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Fudge'] = {}
            ensback = {}    # Dictionary of {EnsLoci name:OMIM gene}
            mutations = {}  # Reorganised dictionary of {gene:{pos:Mutation}}
            for gene in rje.sortKeys(self.dict['Mutations']):
                try: seq = seqdict[genecards[gene]['EnsLoci']]
                except:
                    self.log.printLog('#MAP','No EnsLoci protein mapped for %s' % gene)
                    continue
                mutations[gene] = {}
                ensback[genecards[gene]['EnsLoci']] = gene
                mutpos = {}     # Dictionary of {pos:AA} to map onto sequence
                for subid in rje.sortKeys(self.dict['Mutations'][gene]):                    
                    (disease,mutation) = self.dict['Mutations'][gene][subid]
                    (wild,pos,mut) = rje.matchExp('(\D\D\D)(\d+)(\D\D\D)',mutation)
                    mutpos[int(pos)] = rje_sequence.aa_3to1[wild.upper()]
                    mutations[gene][int(pos)] = self.dict['Mutations'][gene][subid]
                self.dict['Fudge'][seq] = seq.fudgeFactor(mutpos)
            self.deBug(self.dict['Fudge'])

            ### ~ [3] Cross-reference to SLiMFinder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            allslims = {}   # Full dictionary of SLiMFinder results matching OMIM genes
            slimomim = []   # List of (gene,pos) overlapping with SLiMs
            outfile = 'rje_omim.slimfinder.tdt'
            dataheaders = string.split('Dataset,Rank,Pattern,Hit,Pos,EndPos,SeqLen,Variant,Match,AbsChg,NetChg,PepSeq,PepDesign',',')
            headers = ['Gene','OMIM','SubID','Mutation','Disease'] + dataheaders
            rje.delimitedFileOutput(self,outfile,headers,delimit='\t',rje_backup=True)
            for file in glob.glob(self.info['SlimDir'] + '*.occ.csv'):      # Potential SLiM
                slimdata = rje.dataDict(self,file,['Pattern','Hit','Pos','Match'],dataheaders,delimit=',')
                for occ in slimdata:
                    if slimdata[occ]['Hit'] in ensback:     # OMIM gene - possible overlap
                        gene = ensback[slimdata[occ]['Hit']]
                        (start,end) = (int(slimdata[occ]['Pos']),int(slimdata[occ]['EndPos']))
                        if gene not in allslims: allslims[gene] = {}
                        allslims[gene][occ] = slimdata[occ]
                        for mpos in mutations[gene]:
                            if start <= (mpos + self.dict['Fudge'][seqdict[genecards[gene]['EnsLoci']]]) <= end:
                                self.log.printLog('#OMIMSLIM','%s %s %s (%d-%d) = %s' % (slimdata[occ]['Dataset'],slimdata[occ]['Hit'],slimdata[occ]['Pattern'],start,end,mutations[gene][mpos]))
                                slimdata[occ]['Gene'] = gene
                                slimdata[occ]['OMIM'] = string.join(self.dict['Records'][gene])
                                slimdata[occ]['Mutation'] = mutations[gene][mpos][1]
                                slimdata[occ]['Disease'] = mutations[gene][mpos][0]
                                rje.delimitedFileOutput(self,outfile,headers,'\t',slimdata[occ])
                                if (gene,mpos) not in slimomim: slimomim.append((gene,mpos))
            
            ### ~ [4] Calculate coverage of SLiMs for "significance" assessment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (inslim,resx,mutx) = (0,0,0)  # No. of residues in SLiMs, total residue count + no. mutations that may overlap
            for gene in mutations:      # These are just the genes that mapped to sequences
                mutx += len(mutations[gene])
                resx += seqdict[genecards[gene]['EnsLoci']].aaLen()
                if gene in allslims:    # Partially covered by SLiMs
                    res = [0] * seqdict[genecards[gene]['EnsLoci']].aaLen()
                    for occ in allslims[gene]:
                        (start,end) = (int(allslims[gene][occ]['Pos'])-1,int(allslims[gene][occ]['EndPos']))
                        res = res[:start] + [1] * (end-start) + res[end-1:]
                    self.deBug('%s %d (%d)' % (gene,sum(res),seqdict[genecards[gene]['EnsLoci']].aaLen()))
                    inslim += sum(res)
            self.log.printLog('#COV','SLiMs have %.1f%% coverage of OMIM gene sequences' % (100.0*inslim/resx))
            self.log.printLog('#MUT','%d mutations that could potentially occur in SLiMs' % mutx)
            self.log.printLog('#PROB','Probability of observed %d mutation overlap = %.4f' % (len(slimomim),rje.binomial(len(slimomim),mutx,float(inslim)/resx,callobj=self)))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <3> ### Main OMIM Parsing Methods                                                                               #
#########################################################################################################################
    def parseOMIM(self):    ### Main parsing method
        '''Main parsing method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Records'] = {}
            self.dict['Mutations'] = {}
            aas = string.split(string.join(rje_sequence.aa_code_3.values()).upper())
            oline = os.path.exists(self.info['Name'])
            (olen,ox,mx) = (len(open(self.info['Name'],'r').readlines()),0.0,0)
            OMIM = open(self.info['Name'],'r')

            ### ~ [2] Extract data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            record = gene = subid = disease = mutation = ''
            av = False      # Whether reading *FIELD* AV for mutation data
            while oline:
                oline = OMIM.readline()
                self.log.printLog('\r#OMIM','Processing OMIM: %.2f%% (%s genes)' % (ox/olen,rje.integerString(len(self.dict['Records']))),newline=False,log=False)
                ox += 100.0
                if not av and oline[:1] != '*': continue
                line = rje.chomp(oline)
                while line[-1:] == ' ': line = line[:-1]
                ## ~ [2a] New record ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line == '*RECORD*': (record,av) = ('',False)
                elif line == '*FIELD* NO':    # New record
                    record = rje.chomp(OMIM.readline())
                    gene = ''
                    ox += 100.0
                ## ~ [2b] Gene ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line == '*FIELD* TI':      # New gene
                    gene = string.split(rje.chomp(OMIM.readline()))[-1]
                    subid = ''
                    av = False
                    ox += 100.0
                ## ~ [2c] Mutations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line == '*FIELD* AV': av = True        # Start of mutation records
                elif av and rje.matchExp('^(\.\d+)',line):  # New subid mutation record
                    subid = rje.matchExp('^(\.\d+)',line)[0]
                    disease = rje.chomp(OMIM.readline())
                    ox += 100.0
                    try: mutation = rje.matchExp('^%s, (\D\D\D\d+\D\D\D)' % gene,rje.chomp(OMIM.readline()))[0]
                    except: continue    # No mutation or not coding change
                    ox += 100.0
                    subaa = rje.matchExp('(\D\D\D)\d+(\D\D\D)',mutation)
                    if subaa[0] not in aas or subaa[1] not in aas: continue
                    if gene not in self.dict['Records']: self.dict['Records'][gene] = [record]
                    if record not in self.dict['Records'][gene]: self.dict['Records'][gene] += [record]
                    if gene not in self.dict['Mutations']: self.dict['Mutations'][gene] = {}
                    mx += 1
                    self.dict['Mutations'][gene][subid] = (disease,mutation)
                        
            ### ~ [3] Finish & Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OMIM.close()
            self.log.printLog('\r#OMIM','Processing OMIM complete! (%s genes; %s mutations)' % (rje.integerString(len(self.dict['Records'])),rje.integerString(mx)))
            self.saveMutations()
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def saveMutations(self):    ### Outputs parsed mutations into a delimited file
        '''Outputs parsed mutations into a delimited file.'''
        try:### ~ [1] Setup output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = ['OMIM_ID','SubID','Gene','Pos','WildAA','MutAA','Disease']
            outfile = 'omim_mutations.tdt'
            rje.delimitedFileOutput(self,outfile,headers,'\t',rje_backup=True)

            ### ~ [2] Output mutations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gene in rje.sortKeys(self.dict['Mutations']):
                for subid in rje.sortKeys(self.dict['Mutations'][gene]):
                    (disease,mutation) = self.dict['Mutations'][gene][subid]
                    (wild,pos,mut) = rje.matchExp('(\D\D\D)(\d+)(\D\D\D)',mutation)
                    datadict = {'OMIM_ID':string.join(self.dict['Records'][gene],'; '),'SubID':subid,'Gene':gene,
                                'Pos':pos,'WildAA':wild,'MutAA':mut,'Disease':disease}
                    rje.delimitedFileOutput(self,outfile,headers,'\t',datadict)
            self.log.printLog('#OUT','OMIM Mutation output to %s complete' % outfile)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def loadMutations(self):    ### Inputs parsed mutations back into dictionaries
        '''Inputs parsed mutations back into dictionaries.'''
        try:### ~ [1] Setup input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Records'] = {}
            self.dict['Mutations'] = {}
            headers = ['OMIM_ID','SubID','Gene','Pos','WildAA','MutAA','Disease']
            infile = 'omim_mutations.tdt'
            if not os.path.exists(infile): return False
            datadict = rje.dataDict(self,infile,headers[:2],headers,'\t')
            mx = len(datadict)

            ### ~ [2] Process into dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dkey in datadict.keys()[0:]:
                data = datadict.pop(dkey)
                record = data['OMIM_ID']
                subid = data['SubID']
                gene = data['Gene']
                mutation = '%s%s%s' % (data['WildAA'],data['Pos'],data['MutAA'])
                disease = data['Disease']
                if gene not in self.dict['Records']: self.dict['Records'][gene] = [record]
                if record not in self.dict['Records'][gene]: self.dict['Records'][gene] += [record]
                if gene not in self.dict['Mutations']: self.dict['Mutations'][gene] = {}
                self.dict['Mutations'][gene][subid] = (disease,mutation)
            self.log.printLog('\r#OMIM','Loaded %s OMIM mutations (%s genes).' % (rje.integerString(mx),rje.integerString(len(self.dict['Records']))))
            return True
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            return False
#########################################################################################################################
### End of SECTION II: OMIM Class                                                                                       #
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
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
        
    ### Rest of Functionality... ###
    try: OMIM(mainlog,cmd_list).run()
        
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
