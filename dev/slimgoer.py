#!/usr/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMGOER
Description:  Short Linear Motif GO Enrichment of RLC
Version:      0.1
Last Edit:    16/09/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module looks for enriched RLC scores for a given SLiM over its reversed/scrambled versions for GO categories.
    Initial versions of this module are simply for proof of principle analysis, coverting SLiMSearch results into a table
    for R analysis. Later versions will hopefully incorporate the SLiMSearch and statistical tests too.

Commandline:
    ### ~ INPUT/OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    resfile=FILE    : Basefile for SLiMSearch results CSV files []
    minocc=X        : Minimum number of occurrences for a given Motif/GO combo [5]
    minhom=X        : Minimum number of homologues for each protein [3]
    maxgenes=X      : Maximum number of genes for a GO term to apply to [2000]
    
See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_go, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_go, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Modifications to improve functionality, reduce run times and return better data.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add Mann-Whitney U tests to module itself / call R and parse results.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SLiMGOER', '0.1', 'September 2008', '2008')
    description = 'Short Linear Motif GO Enrichment of RLC'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
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
### SECTION II: GOER Class                                                                                               #
#########################################################################################################################
class GOER(rje.RJE_Object):     
    '''
    GOER Class. Author: Rich Edwards (2008).

    Info:str
    - ResFile = SLiMSearch results CSV file (full occ) []
    
    Opt:boolean

    Stat:numeric
    - MinHom = Minimum number of homologues for each protein [3]
    - MinOcc = Minimum number of occurrences for a given Motif/GO combo [5]
    - MaxGenes = Maximum number of genes for a GO term to apply to [2000]

    List:list

    Dict:dictionary
    - Occ = SLiMSearch results {motif:{type:{seq:[occdicts]}}}

    Obj:RJE_Objects
    - GO = rje_go.GO object containing GO Data
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['EnsGOPath','ResFile']
        self.optlist = []
        self.statlist = ['MinHom','MinOcc','MaxGenes']
        self.listlist = []
        self.dictlist = ['Occ']
        self.objlist = ['GO']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'MinOcc':5,'MaxGenes':2000})
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
                self._cmdReadList(cmd,'file',['ResFile'])
                self._cmdReadList(cmd,'int',['MinHom','MinOcc','MaxGenes'])  
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.readSLiMSearch()
            self.makeGOFile()
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup GO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['GO'] = rje_go.GO(self.log,self.cmd_list)
            self.obj['GO'].readGO()
            ### ~ [2] Setup EnsEMBL GO Linkage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['GO'].mapEnsGO()
            self.obj['GO'].makeGOGenes()    # Makes a dictionary of {GO:[Genes]}
            self.reduceGO()     # Reduce GO terms to those with enough sequences 
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### GO Methods                                                                                              # 
#########################################################################################################################
    def ensGO(self,gene):   ### Returns list of GO terms associated with gene (if any)
        '''Returns list of GO terms associated with gene (if any).'''
        if gene not in self.obj['GO'].dict['EnsGO']: return []
        return self.obj['GO'].dict['EnsGO'][gene][0:]
#########################################################################################################################
    def reduceGO(self):   ### Reduce GO terms to those with enough sequences
        '''Reduce GO terms to those with enough sequences.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            go = self.obj['GO']
            minocc = self.stat['MinOcc']
            maxocc = self.stat['MaxGenes']
            gokey = 'EnsGO'
            ### ~ [2] ~ Reduce GO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for id in rje.sortKeys(go.dict['GO']):
                try:
                    idgenes = go.getGOGenes(id,gokey)
                    self.deBug('%s: %s' % (id,idgenes))
                    if len(idgenes) < minocc or len(idgenes) > maxocc:   # Remove
                        go.dict['GO'].pop(id)
                        for gene in idgenes:
                            go.dict[gokey][gene].remove(id)
                            if not go.dict[gokey][gene]: go.dict[gokey].pop(gene)
                except: self.errorLog('GOER.reduceGO(%s) problem' % id)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### SLiMSearch Methods                                                                                      # 
#########################################################################################################################
    def readSLiMSearch(self):   ### Reads SLiMSearch results into data dictionary
        '''Reads SLiMSearch results into data dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sumfile = '%s.summary.csv' % self.info['ResFile']
            occfile = '%s.csv' % self.info['ResFile']
            if not os.path.exists(sumfile): return self.errorLog('No Summary file "%s"!' % sumfile,printerror=False)
            if not os.path.exists(occfile): return self.errorLog('No Occurrence file "%s"!' % occfile,printerror=False)
            ### ~ [2] Read Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            esum = rje.dataDict(self,sumfile,mainkeys=['Motif'],datakeys='All',getheaders=False)
            occmotifs = []      # List of motifs with enough occurrences
            for motif in rje.sortKeys(esum):
                if string.atoi(esum[motif]['N_Occ']) < self.stat['MinOcc']: continue
                occmotifs.append(motif)
            ### ~ [3] Read Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#MOTIF','%d motifs with N_Occ >= MinOcc (%d)' % (len(occmotifs),self.stat['MinOcc']))
            self.readSLiMSearchOcc(occmotifs)
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def readSLiMSearchOcc(self,motifs=[]):   ### Reads SLiMSearch results into data dictionary
        '''Reads SLiMSearch results into data dictionary.'''
        try:### ~ [1] Read ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not motifs: self.printLog('#OCC','Cannot process occurrences for No motifs!')
            occfile = '%s.csv' % self.info['ResFile']
            delimit = rje.delimitFromExt(filename=occfile)
            data = rje.dataDict(self,occfile,mainkeys=['Motif','Seq','Start_Pos','End_Pos'],datakeys=rje.split('Seq,Desc,Start_Pos,End_Pos,Cons,HomNum,GlobID,LocID,Hyd,SA',','))
            self.dict['Occ'] = {}
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mx,ox,otot) = (0,0.0,len(data))
            for occ in data:
                self.progLog('\r#OCC','Processing occurrences (%d motifs): %.2f%%' % (mx,ox/otot)); ox += 100.0
                #x#self.deBug('%s vs MinHom %d' % (data[occ],self.stat['MinHom']))
                if string.atoi(data[occ]['HomNum']) < self.stat['MinHom']: continue
                (motif,seq,start,end) = rje.split(occ,delimit)
                if motif not in motifs: continue
                try:
                    gene = rje.matchExp('gene:(\S+)\]',data[occ]['Desc'])[0]
                    self.deBug('%s:%s' % (gene,self.ensGO(gene)))
                    if not self.ensGO(gene): continue
                except: continue
                if motif[-3:] == 'rev': (motif,type) = (motif[:-4],'Rev')
                elif motif[-5:] == 'scram': (motif,type) = (motif[:-6],'Scr')
                else: type = 'ELM'
                if motif not in self.dict['Occ']: self.dict['Occ'][motif] = {}; mx += 1
                if type not in self.dict['Occ'][motif]: self.dict['Occ'][motif][type] = {}
                if gene not in self.dict['Occ'][motif][type]: self.dict['Occ'][motif][type][gene] = []
                self.dict['Occ'][motif][type][gene].append(data[occ])
            self.printLog('\r#OCC','Processed %s occurrences: %d motifs with GO-links' % (rje.integerString(otot),mx))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def makeGOFile(self):   ### Maps GO to sequences and outputs table for R analysis
        '''Maps GO to sequences and outputs table for R analysis.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.goer.tdt' % self.info['ResFile']
            headers = ['GOID','Motif','Type','Gene','Cons','HomNum','GlobID','LocID','Hyd','SA']
            rje.delimitedFileOutput(self,outfile,headers,rje_backup=True)
            ### ~ [2] ~ Work through dictionary and output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (mx,mtot) = (-100.0,len(self.dict['Occ']))
            for motif in rje.sortKeys(self.dict['Occ']):
                mx += 100.0; self.progLog('\r#OUT','Generating %s output: %.1f%% (%s|CheckSeq)         ' % (outfile,(mx/mtot),motif))
                ## ~ [2a] ~ Check MinOcc in terms of sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for type in rje.sortKeys(self.dict['Occ'][motif]):
                    if len(self.dict['Occ'][motif][type]) < self.stat['MinOcc']: self.dict['Occ'][motif].pop(type)
                if 'ELM' not in self.dict['Occ'][motif] or len(self.dict['Occ'][motif]) < 2: continue
                for type in self.dict['Occ'][motif]:
                    ## ~ [2b] ~ Map GO terms and check MinOcc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.progLog('\r#OUT','Generating %s output: %.1f%% (%s|Check%s) ' % (outfile,(mx/mtot),motif,type)); 
                    godict = {}     # Temp dictionary of {GOID:[Seqs]}
                    for gene in self.dict['Occ'][motif][type]:
                        for go in self.ensGO(gene):
                            if go not in godict: godict[go] = [gene]
                            else: godict[go].append(gene)
                    self.progLog('\r#OUT','Generating %s output: %.1f%% (%s|OccGO%s) ' % (outfile,(mx/mtot),motif,type)); 
                    for go in rje.sortKeys(godict):
                        if len(godict[go]) < self.stat['MinOcc']: godict.pop(go)
                    ## ~ [2c] ~ Output remaining GO terms occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.progLog('\r#OUT','Generating %s output: %.1f%% (%s|Output%s)' % (outfile,(mx/mtot),motif,type)); 
                    for go in rje.sortKeys(godict):
                        for gene in godict[go]:
                            for occdict in self.dict['Occ'][motif][type][gene]:
                                datadict = rje.combineDict({'GOID':'GO:%s' % go,'Motif':motif,'Type':type,'Gene':gene},occdict)
                                rje.delimitedFileOutput(self,outfile,headers,datadict=datadict)
                self.printLog('#OUT','Output for %s %s complete.' % (motif,rje.sortKeys(self.dict['Occ'][motif])),screen=False)
            self.printLog('\r#OUT','Generating %s output complete!         ' % (outfile))
        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: New Class                                                                                        #
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
    try:GOER(mainlog,cmd_list).run()

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
