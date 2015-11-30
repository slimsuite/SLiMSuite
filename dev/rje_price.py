#!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_price
Description:  Experimental Price Equation Tool
Version:      0.0
Last Edit:    05/01/09
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    Current fitness scores:
    * 

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE  : Input sequence alignment in fasta format (over-rides batch=LIST) [None]
    batch=LIST  : List of alignment files to use as input [*.fas,*.fasta]
    query=X     : Identifier of query sequence [None]
    dna=T/F     : Whether sequences are DNA or protein [True]
    ### ~ PRICE EQUATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    fitness=X   : Fitness measurement [cons]
    phenotype=X : Phenotype measurement [cons]
    seqgroup=X  : Sequence grouping method [triplets]
    qrygaps=T/F : Whether to include gaps in the query sequence as positions to score [False]
    special=X   : Instigate special run, e.g. allbyall [None]
    normfit=X   : Normalise fitness to have mean of 1 [False]
    weighted=T/F: Weight the mean covariance by size of group [True]  
    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    resfile=FILE: Results file [price.tdt]
    append=T/F  : Append results file if already present [True]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : List here
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_PRICE', '0.0', 'August 2009', '2009')
    description = 'Experimental Price Equation Tool'
    author = 'Dr Richard J. Edwards & Dr Joel D. Parker.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### SECTION II: Price Class                                                                                             #
#########################################################################################################################
class Price(rje.RJE_Object):     
    '''
    Price Class. Author: Rich Edwards (2009).

    Info:str
    - Fitness = Fitness measurement [cons]
    - Phenotype = Phenotype measurement [cons]
    - ResFile = Results file [price.tdt]
    - SeqGroup = Sequence grouping method [triplets]
    - Special = Instigate special run, e.g. allbyall [None]

    Opt:boolean
    - NormFit = Normalise fitness to have mean of 1 [False]
    - QryGaps = Whether to include gaps in the query sequence as positions to score [False]
    - Weighted = Weight the mean covariance by size of group [False]

    Stat:numeric

    List:list
    - Batch = List of alignment files to use as input [*.fas,*.fasta]
    - Fitness = Fitness measurement vector (matches query sequence) 
    - Phenotype = Phenotype measurement (matches query sequence) 
    - SeqGroup = Sequence grouping method (matches query sequence)

    Dict:dictionary    

    Obj:RJE_Objects
    - SeqList = Sequence list object
    '''
#########################################################################################################################
    def qry(self): return self.obj['SeqList'].obj['QuerySeq']
    def seqs(self): return self.obj['SeqList'].seqs()
    def dna(self): return self.obj['SeqList'].dna()
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Fitness','Phenotype','ResFile','SeqGroup','Special']
        self.optlist = ['QryGaps','NormFit','Weighted']
        self.statlist = []
        self.listlist = ['Batch','Fitness','Phenotype','SeqGroup']
        self.dictlist = []
        self.objlist = ['SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Fitness':'cons','Phenotype':'cons','SeqGroup':'triplets','ResFile':'price.tdt'})
        self.setOpt({'Append':True,'Weighted':True})
        self.list['Batch'] = ['*.fas','*.fasta']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['SeqList'] = rje_seq.SeqList(self.log,['query=1']+self.cmd_list+['autoload=T'])
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdReadList(cmd,'info',['Fitness','Phenotype','SeqGroup','Special'])
                self._cmdReadList(cmd,'file',['ResFile'])
                self._cmdReadList(cmd,'opt',['QryGaps','NormFit','Weighted'])
                self._cmdReadList(cmd,'list',['Batch'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch: self.setupResults()
            ## ~ [1b] ~ Batch run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch and not self.obj['SeqList'].seqs():    ### Look for batch files and run for each
                batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                if not batchfiles: self.errorLog('No input files found!',printerror=False)
                else:
                    bx = 0
                    for infile in batchfiles:
                        bx += 1
                        self.printLog('#BATCH','Batch running %s' % infile)
                        bcmd = ['query=1']+self.cmd_list+['autoload=T','seqin=%s' % infile]
                        self.obj['SeqList'] = rje_seq.SeqList(self.log,bcmd)
                        self.run(batch=True)
                        self.opt['Append'] = True
                        self.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.opt['Win32'] and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return
            ## ~ [1c] ~ Special run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['Special'].lower() == 'allbyall':
                self.printLog('#RUN','Performing special "all-by-all" pairwise run')
                self.info['Special'] = ''
                for i in range(len(self.seqs())-1):
                    self.obj['SeqList'].obj['QuerySeq'] = self.seqs()[i]
                    for j in range(i+1,len(self.seqs())):
                        self.info['Fitness'] = self.info['Phenotype'] = '%d' % (j + 1)
                        self.run(batch=True)
                        self.opt['Append'] = True
                self.info['Special'] = 'allbyall'; return                
            ## ~ [1d] ~ General setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setup()
            ### ~ [2] ~ Price calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.fitness()
            self.phenotype()
            self.grouping()
            for vector in ['Fitness','Phenotype','SeqGroup']:
                if len(self.list[vector]) != self.qry().seqLen():
                    self.errorLog('%s vector length (%s) does not match %s sequence length (%s)' % (vector,len(self.list[vector]),self.qry().seqLen()),printerror=False)
                    raise ValueError
            results = self.price()
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            results['Dataset'] = rje.baseFile(self.obj['SeqList'].info['Name'],True)
            results['Query'] = self.qry().shortName()
            results['Fitness'] = self.info['Fmethod']
            results['Phenotype'] = self.info['Pmethod']
            results['SeqGroup'] = self.info['SeqGroup']
            rje.delimitedFileOutput(self,self.info['ResFile'],self.list['Headers'],datadict=results)
            self.printLog('#OUT','Results output to %s' % self.info['ResFile'])
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupResults(self):    ### Main results setup method.
        '''Main results setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Headers'] = ['Dataset','Query','Fitness','Phenotype','SeqGroup','CovP','CovB','CovW','Price','Ratio']
            rje.delimitedFileOutput(self,self.info['ResFile'],self.list['Headers'],rje_backup=True)
        except: self.errorLog('Problem during %s setupResults().' % self); raise
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList'] 
            seqlist._checkAln(aln=True,realign=True)
            if not seqlist.obj['QuerySeq']:
                seqlist.obj['QuerySeq'] = seqlist.seqs()[0]
                self.printLog('#QRY','No query sequence: will use %s' % seqlist.obj['QuerySeq'].shortName())
        except: self.errorLog('Problem during %s setup.' % self); raise
#########################################################################################################################
    ### <3> ### Price Equation Methods                                                                                  #
#########################################################################################################################
    def fitness(self):  ### Calculates fitness vector
        '''Calculates fitness vector.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            methodlist = ['cons','seqnumber']
            self.info['Fmethod'] = method = self.info['Fitness'].lower()
            if method not in methodlist:
                try:
                    method = string.atoi(method)
                    try: method = self.seqs()[method-1]
                    except: self.errorLog('Cannot use sequence "%s" for comparison!' % method); raise
                    self.info['Fmethod'] = method.shortName()
                except: 
                    self.errorLog('Fitness method "%s" not recognised!' % method,printerror=False)
                    self.errorLog('Check fitness=%s' % string.join(methodlist,'/'),printerror=False)
                    raise ValueError
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if method == 'cons': self.list['Fitness'] = self.posPercID()
            elif method in self.seqs(): self.list['Fitness'] = self.posPercID(comp=method)
            elif os.path.exists(method):
                self.list['Fitness'] = rje.listFromCommand(method,checkfile=True)
                self.printLog('#FIT','Vector of %s fitness values read from %s' % (len(self.list['Fitness']),method))
            return
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   
#########################################################################################################################
    def phenotype(self):  ### Calculates phenotype vector
        '''Calculates phenotype vector.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            methodlist = ['cons','seqnumber','hyd']
            self.info['Pmethod'] = method = self.info['Phenotype'].lower()
            if method not in methodlist:
                try:
                    method = string.atoi(method)
                    try: method = self.seqs()[method-1]
                    except: self.errorLog('Cannot use sequence "%s" for comparison!' % method); raise
                    self.info['Pmethod'] = method.shortName()
                except: 
                    self.errorLog('Phenotype method "%s" not recognised!' % method,printerror=False)
                    self.errorLog('Check phenotype=%s' % string.join(methodlist,'/'),printerror=False)
                    raise ValueError
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if method == 'cons': self.list['Phenotype'] = self.posPercID()
            elif method in self.seqs(): self.list['Phenotype'] = self.posPercID(comp=method)
            elif method == 'hyd': self.list['Phenotype'] = rje_sequence.eisenbergHydropathy(self.qry().info['Sequence'],returnlist=True)
            elif os.path.exists(method):
                self.list['Phenotype'] = rje.listFromCommand(method,checkfile=True)
                self.printLog('#PHEN','Vector of %s phenotype values read from %s' % (len(self.list['Phenotype']),method))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   
#########################################################################################################################
    def grouping(self):  ### Calculates grouping vector
        '''Calculates grouping vector.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            methodlist = ['triplets','codons','casechange','case','disorder']
            method = self.info['SeqGroup'].lower()
            if method not in methodlist:
                self.errorLog('SeqGroup method "%s" not recognised!' % method,printerror=False)
                self.errorLog('Check seqgroup=%s' % string.join(methodlist,'/'),printerror=False)
                raise ValueError
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if method == 'triplets': self.list['SeqGroup'] = self.triplets()
            elif method == 'codons': self.list['SeqGroup'] = self.codons()
            elif method == 'casechange': self.list['SeqGroup'] = self.caseChange()
            elif method == 'case': self.list['SeqGroup'] = self.case()
            elif method == 'disorder':
                if self.opt['QryGaps']: self.list['SeqGroup'] = self.qry().gappedDisorder()
                else: self.list['SeqGroup'] = self.qry().gappedDisorder(gap=None)
                for i in range(self.qry().seqLen()):
                    if self.list['SeqGroup'][i]:
                        if self.list['SeqGroup'][i] > self.qry().obj['Disorder'].stat['IUCut']: self.list['SeqGroup'][i] = 'Dis'
                        else: self.list['SeqGroup'][i] = 'Ord'
            elif os.path.exists(method):
                self.list['SeqGroup'] = rje.listFromCommand(method,checkfile=True)
                self.printLog('#GRP','Vector of %s group values read from %s' % (len(self.list['SeqGroup']),method))
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   
#########################################################################################################################
    def price(self):  ### Calculates price equation, using Fitness, Phenotype and SeqGroup vectors
        '''Calculates price equation, using Fitness, Phenotype and SeqGroup vectors.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pop = {'z':[],'w':[]}       # w = fitness, z = phenotype
            grp = {}                    # Each group will have its own w and z
            grpmean = {'z':[],'w':[]}   # Calculate means for each group
            grpcov = []                 # List of group covariances
            ### ~ [2] Populate data vectors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(self.list['SeqGroup'])
            self.deBug(self.list['Fitness'])
            self.deBug(self.list['Phenotype'])
            for i in range(len(self.list['SeqGroup'])):
                if not self.list['SeqGroup'][i]: continue
                g = self.list['SeqGroup'][i]
                w = self.list['Fitness'][i]
                z = self.list['Phenotype'][i]
                pop['z'].append(z); pop['w'].append(w)
                if g not in grp: grp[g] = {'z':[],'w':[]}
                grp[g]['z'].append(z); grp[g]['w'].append(w)
            ## ~ [2a] Normalise fitness? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['NormFit']:
                meanfit = float(rje.meansd(pop['w'])[0])
                for i in range(len(pop['w'])): pop['w'][i] = pop['w'][i] / meanfit
                for g in grp:
                    for i in range(len(grp[g]['w'])): grp[g]['w'][i] = grp[g]['w'][i] / meanfit
            ## ~ [2b] Group means ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            covw = 0.0                      # Mean covariance within groups
            for g in grp:
                grp[g]['cov'] = self.covariance(grp[g]['z'],grp[g]['w'])
                grpcov.append(grp[g]['cov'])
                if self.opt['Weighted']: covw += grp[g]['cov'] * len(grp[g]['w']) / len(pop['w'])
                else: covw += grp[g]['cov'] / len(grp)
                grpmean['z'].append(rje.meansd(grp[g]['z'])[0])
                grpmean['w'].append(rje.meansd(grp[g]['w'])[0])
            ### ~ [3] Calculate within and between group covariance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            covp = self.covariance(pop['z'],pop['w'])           # Covariance of whole population
            covb = self.covariance(grpmean['z'],grpmean['w'])   # Covariance between groups
            #x#covw = rje.meansd(grpcov)[0]                        # Mean covariance within groups
            price = covp / rje.meansd(pop['w'])[0]
            try: ratio = covb / covw
            except: ratio = -1
            ## ~ [3a] Perform checks of calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#CHECK','CovP = %s; (CovB + CovW) = %s' % (rje.expectString(covp),rje.expectString(covb+covw)))
            self.printLog('#PRICE','Price value = %s; CovB/CovW ratio = %s' % (rje.expectString(price),rje.expectString(ratio)))
            return {'CovP':rje.expectString(covp),'CovB':rje.expectString(covb),'CovW':rje.expectString(covw),'Price':rje.expectString(price),'Ratio':rje.expectString(ratio)}
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   
#########################################################################################################################
    def covariance(self,list1,list2):   ### Calculates the covariance of two lists and returns
        '''Calculates the covariance of two lists and returns.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            n = len(list1)
            if not n: self.errorLog('Lists for covariance are empty!',printerror=False); return 0.0
            if len(list2) != n: self.errorLog('Lists for covariance of different lengths!',printerror=False); raise ValueError
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return covariance(list1,list2)
        except: self.errorLog(rje_zen.Zen().wisdom()); raise   
#########################################################################################################################
    ### <4> ### Fitness/Phenotype Methods                                                                               #
#########################################################################################################################
    def posPercID(self,gaps=True,xval=0.0,default=1.0,comp=None):  ### Returns a list of absolute pecentage conservation across each position
        '''
        Returns a list of absolute pecentage conservation across each position.
        >> gaps:bool [True] = Whether to include gapped sequences in calculation [True]
        >> xval:num [0.0] = The value (0-1) to give undefined residues matching defined residues
        >> default:num [1.0] = Value to return if no homologues for position
        >> comp:Sequence object = sequence for pairwise comparison
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qry = self.qry()
            if comp: compseq = [comp]
            else: compseq = self.seqs()[0:]; compseq.remove(qry)
            poslist = [default] * qry.seqLen()    # List of percentage ID values
            xval = min(1.0,max(0.0,xval))
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for r in range(qry.seqLen()):
                q = qry.info['Sequence'].upper()[r]
                i = 0.0; n = 0
                for seq in compseq:
                    s = seq.info['Sequence'].upper()[r]
                    if s == q: i += 1; n += 1
                    elif 'X' in [s,q]: i += xval; n += 1
                    elif s == '-' and not gaps: continue
                    else: n += 1
                if n: poslist[r] = i / n
            return poslist
        except: self.errorLog(rje_zen.Zen().wisdom()); raise
#########################################################################################################################
    ### <5> ### Grouping Methods                                                                                        #
#########################################################################################################################
    def triplets(self): ### Returns grouping vector based on DNA triplets
        '''Returns grouping vector based on DNA triplets.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qry = self.qry()
            grplist = [0] * qry.seqLen()    # List of groups (0 = no group)
            gx = 0
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            trip = 0
            for r in range(qry.seqLen()):
                q = qry.info['Sequence'].upper()[r]
                if self.opt['QryGaps'] and q == '-': continue
                if not trip: gx += 1
                grplist[r] = gx
                if trip == 2: trip = 0
                else: trip += 1
            self.printLog('#GRP','%s triplet groups from %s%s' % (rje.integerString(gx),rje.integerString(qry.seqLen()),self.obj['SeqList'].units()))
            return grplist
        except: self.errorLog(rje_zen.Zen().wisdom()); raise           
#########################################################################################################################
    def codons(self): ### Returns grouping vector based on DNA codon positions (three groups)
        '''Returns grouping vector based on DNA codon positions (three groups).'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qry = self.qry()
            grplist = [0] * qry.seqLen()    # List of groups (0 = no group)
            ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            trip = 1
            for r in range(qry.seqLen()):
                q = qry.info['Sequence'].upper()[r]
                if self.opt['QryGaps'] and q == '-': continue
                grplist[r] = trip
                if trip == 3: trip = 1
                else: trip += 1
            self.printLog('#GRP','3 codon groups from %s%s' % (rje.integerString(qry.seqLen()),self.obj['SeqList'].units()))
            return grplist
        except: self.errorLog(rje_zen.Zen().wisdom()); raise           
#########################################################################################################################
    def caseChange(self):     ### Returns groupings based on Case boundaries of query
        '''Returns groupings based on Case boundaries of query.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qry = self.qry()
            self.deBug(qry.getSequence(case=True))
            grplist = ['UC'] * qry.seqLen()    # List of groups (None = no group)
            ### ~ [1] Map Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for (start,end) in self.qry().dict['Case']['Lower']:
                for i in range(start-1,end): grplist[i] = 'LC'
            caselist = grplist[0:]
            gx = 1
            for r in range(qry.seqLen()):
                q = qry.info['Sequence'].upper()[r]
                if not self.opt['QryGaps'] and q == '-': grplist[r] = 0
                elif r > 0 and caselist[r] != caselist[r-1]: gx += 1
                grplist[r] = gx
            self.printLog('#GRP','%s case groups from %s%s' % (rje.integerString(gx),rje.integerString(qry.seqLen()),self.obj['SeqList'].units()))
            self.deBug(grplist)
            return grplist
        except: self.errorLog(rje_zen.Zen().wisdom()); raise           
#########################################################################################################################
    def case(self):     ### Returns groupings based on Case boundaries of query
        '''Returns groupings based on Case boundaries of query.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qry = self.qry()
            self.deBug(qry.getSequence(case=True))
            grplist = ['UC'] * qry.seqLen()    # List of groups (None = no group)
            ### ~ [1] Map Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for (start,end) in self.qry().dict['Case']['Lower']:
                for i in range(start-1,end): grplist[i] = 'LC'
            if self.opt['QryGaps']: return grplist
            for r in range(qry.seqLen()):
                if qry.info['Sequence'].upper()[r] == '-': grplist[r] = None
            self.deBug(grplist)
            return grplist
        except: self.errorLog(rje_zen.Zen().wisdom()); raise           
#########################################################################################################################
### End of SECTION II: Price Class                                                                                      #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def covariance(list1,list2,biased=True):   ### Calculates the covariance of two lists and returns
    '''Calculates the covariance of two lists and returns.'''
    n = len(list1)
    x = list1; y = list2
    ### ~ [2] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sumprod = 0.0
    for i in range(n): sumprod += x[i] * y[i]
    cov = (sumprod - (1.0/n) * sum(x) * sum(y))
    if biased: return (1.0 / n) *  cov
    else: return (1.0 / (n-1)) *  cov
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
    try: Price(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
