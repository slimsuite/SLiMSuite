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
Module:       QSF_Analysis
Description:  Custom QSF analysis pipeline
Version:      0.0
Last Edit:    08/08/09
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_zen
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
    (program, version, last_edit, copyright) = ('QSF_Analysis', '0.0', 'August 2009', '2009')
    description = 'Custom QSF analysis pipeline'
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
### SECTION II: QSF_Analysis                                                                                            #
#########################################################################################################################
class QSF_Analysis(rje.RJE_Object):     
    '''
    QSF_Analysis Class. Author: Rich Edwards (2009).

    Info:str
    
    Opt:boolean

    Stat:numeric

    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
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
                ### Class Options ### 
                #self._cmdRead(cmd,type='info',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.tabulatePPIRegion()
            self.mapRegionsToSequences()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Pairwise PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ppipairwise = '/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/Pingu/pingu.pairwise.tdt'
            self.progLog('\r#PPI','Loading pairwise data...')
            pairwise = rje.dataDict(self,ppipairwise,['Hub','Spoke'],['Spoke','SpokeSeq','Evidence'])
            gene2seq = {}; seq2gene = {}
            fullppi = {}; px = 0.0; ptot = len(pairwise); ppix = 0
            for pair in rje.sortKeys(pairwise):
                self.progLog('\r#PPI','Processing full pairwise PPI: %.2f%%' % (px/ptot)); px += 100.0
                [hub,spoke] = string.split(pair,'\t')
                if spoke not in gene2seq:
                    sseq = pairwise[pair]['SpokeSeq']
                    gene2seq[spoke] = sseq; seq2gene[string.split(sseq,'__')[0]] = spoke
                if hub not in fullppi: fullppi[hub] = {}
                if spoke not in fullppi[hub]: fullppi[hub][spoke] = pairwise.pop(pair)['Evidence']; ppix += 1
            self.printLog('\r#PPI','Processed full pairwise PPI: %s genes; %s ppi.' % (rje.integerString(len(fullppi)),rje.integerString(ppix/2)))
            ### ~ [2] Filter complexes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodppifile = '/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/Pingu/hybrid.txt'
            goodppi = self.loadFromFile(goodppifile,chomplines=True)
            self.dict['PPI'] = {}
            px = 0.0; ptot = len(fullppi); fppix = ppix; ppix = 0
            for hub in fullppi:
                self.progLog('\r#PPI','Filtering complexes: %.2f%% (%s hubs; %s ppi)' % (px/ptot,rje.integerString(len(self.dict['PPI'])),rje.integerString(ppix))); px +=100.0
                self.dict['PPI'][hub] = []
                for spoke in fullppi[hub]:
                    goodspoke = False
                    for ptype in goodppi:
                        if rje.matchExp(':(%s)($|\|)' % ptype, fullppi[hub][spoke]): goodspoke = True; break
                    if goodspoke: self.dict['PPI'][hub].append(spoke); continue
                    goodspoke = True
                    for spoke2 in fullppi[hub]:
                        if spoke2 in [hub,spoke]: continue
                        if spoke2 in fullppi[spoke]: goodspoke = False; break
                    if goodspoke: self.dict['PPI'][hub].append(spoke)
                ppix += len(self.dict['PPI'][hub])
                if not self.dict['PPI'][hub]: self.dict['PPI'].pop(hub)
            self.printLog('\r#PPI','Filtered complexes: (%s -> %s hubs; %s -> %s ppi)' % (rje.integerString(len(fullppi)),rje.integerString(len(self.dict['PPI'])),rje.integerString(fppix/2),rje.integerString(ppix/2)))
            ### ~ [3] SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = '/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/EnsEMBL/ens_HUMAN.loci.fas'
            scmd = ['accnr=F','seqnr=F','seqin=%s' % seqfile] + self.cmd_list + ['autoload=T']
            seqlist = self.obj['SeqList'] = rje_seq.SeqList(self.log,scmd)
            self.dict['SeqObj'] = seqlist.seqNameDic('Max')
            self.dict['Gene2Seq'] = gene2seq; self.dict['Seq2Gene'] = seq2gene
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def tabulatePPIRegion(self):    ### Tabulates regions of known PPI from DAT file
        '''Tabulates regions of known PPI from DAT file.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tabfile = 'ppi_region.tdt'
            unifile = '/scratch/RJE_Filestore/SBSBINF/Databases/DBase_090505/UniFake/Human/ens_HUMAN.unifake.dat'
            if os.path.exists(tabfile) and not self.opt['Force']: return self.printLog('#REGTAB','%s found. (Force=F)' % tabfile)
            headers = ['Protein','Start','End','Interactor']
            rje.delimitedFileOutput(self,tabfile,headers,rje_backup=True)
            ### ~ [2] Extract and tabulate data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gcmd = "grep -P '(ID   |REGION)' %s | grep -P '(HUMAN|interact)' -i | grep REGION -B 1" % unifile
            self.printLog('#GREP',gcmd)
            prot = None; rx = 0; plist = []; ilist = []
            for gline in os.popen(gcmd).readlines():
                if rje.matchExp('ID   (\S+)',gline): prot = rje.matchExp('ID   (\S+)',gline)[0]
                if rje.matchExp('FT   REGION\s+(\d+)\s+(\d+).+nteract\S+ with (\S.+)',gline):
                    (rstart,rend,rint) = rje.matchExp('FT   REGION\s+(\d+)\s+(\d+).+nteract\S+ with (\S.+)',gline)
                    for ppi in string.split(rint):
                        if rje.matchExp('^([A-Z0-9][A-Z0-9]+)',ppi):
                            datadict = {'Protein':prot,'Start':rstart,'End':rend,'Interactor':rje.matchExp('^([A-Z0-9][A-Z0-9]+)',ppi)[0]}
                            rje.delimitedFileOutput(self,tabfile,headers,datadict=datadict); rx += 1
                            if prot not in plist: plist.append(prot)
                            if datadict['Interactor'] not in ilist: ilist.append(datadict['Interactor'])
                            self.progLog('\r#REGTAB','Tabulating regions: %s proteins; %s interactors; %s regions' % (rje.integerString(len(plist)),rje.integerString(len(ilist)), rje.integerString(rx)))
            self.printLog('\r#REGTAB','Tabulated regions (%s proteins; %s interactors; %s regions) => %s' % (rje.integerString(len(plist)),rje.integerString(len(ilist)),rje.integerString(rx),tabfile))
            return True
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def mapRegionsToSequences(self):    ### Maps tabulates PPI regions onto sequence datasets
        '''Maps tabulates PPI regions onto sequence datasets.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            minseq = 3
            outdir = 'RegPPI/'
            adddir = 'RegPPIAdd/'
            rje.mkDir(self,outdir)
            rje.mkDir(self,adddir)
            tabfile = 'ppi_region.tdt'
            region = rje.dataDict(self,tabfile,['Interactor','Protein'],['Start','End'],lists=True)
            ### ~ [2] Work through each pair in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            px = 0.0; ptot = len(region); fx = 0
            for pair in rje.sortKeys(region):
                self.progLog('\r#FAS','Generating fasta files: %.2f%%' % (px/ptot)); px += 100.0
                ## ~ [2a] Map sequences to PPI dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                [hub, spoke] = string.split(pair,'\t')
                try: qryseq = self.dict['SeqObj'][spoke]
                except: self.printLog('\n#QRY','Spoke gene "%s" missing from Sequence file' % spoke); continue
                try: spoke = self.dict['Seq2Gene'][spoke]
                except: self.printLog('\n#QRY','Spoke protein "%s" missing from PPI dictionary' % spoke); continue
                if hub not in self.dict['PPI']: self.printLog('\n#HUB','Hub gene "%s" missing from PPI dictionary' % hub); continue
                addspoke = spoke not in self.dict['PPI'][hub]
                if addspoke:
                    self.dict['PPI'][hub].append(spoke)
                    self.printLog('\n#PPI','Added spoke gene "%s" to hub "%s" interactome' % (spoke,hub))
                if len(self.dict['PPI'][hub]) < minseq: self.printLog('\n#HUB','Hub "%s" interactome too small (<%s spokes)' % (hub,minseq)); continue
                ## ~ [2b] Identify query sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                reglist = []
                for pos in region[pair]['Start'] + region[pair]['End']: reglist.append(string.atoi(pos))
                reglist.sort()
                qsequence = qryseq.info['Sequence'].lower()
                self.deBug(len(qsequence))
                self.deBug(qsequence)
                prelen = len(qsequence)
                while reglist:
                    self.deBug(reglist)
                    try: startx = reglist.pop(0) - 1; endx = reglist.pop(0)
                    except: self.errorLog('%s PPI Region problem: %s' % (pair,region[pair])); continue
                    self.deBug(qsequence[startx-1:endx+1].upper())
                    qsequence = qsequence[:startx] + qsequence[startx:endx].upper() + qsequence[endx:]
                self.deBug(qsequence)
                if len(qsequence) != prelen:
                    self.printLog('#FUCK','%s' % region[pair])
                    self.printLog('#FUCK',qryseq.info['Sequence'].lower())
                    self.printLog('#FUCK',qsequence)
                    raise ValueError
                ## ~ [2c] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if addspoke: outfile = '%s%s.%s.fas' % (adddir,hub,spoke); ox = 1
                else: outfile = '%s%s.%s.fas' % (outdir,hub,spoke); ox = 1
                open(outfile,'w').write('>%s\n%s\n' % (qryseq.info['Name'],qsequence))
                for spoke2 in self.dict['PPI'][hub]:
                    if spoke2 == spoke: continue
                    try:
                        sseq = self.dict['SeqObj'][self.dict['Gene2Seq'][spoke2]]
                        open(outfile,'a').write('>%s\n%s\n' % (sseq.info['Name'],sseq.info['Sequence']))
                        ox += 1
                    except: pass
                self.printLog('\n#FAS','%s sequences output to %s' % (ox,outfile))
                    
                
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: QSF_Analysis                                                                                     #
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
    try: QSF_Analysis(mainlog,cmd_list).run()

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
