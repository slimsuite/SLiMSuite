#!/usr/local/bin/python

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
Module:       SLiMPID
Description:  Short Linear Motif Protein Interaction Datasets
Version:      0.0
Last Edit:    04/06/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    To be added.

Commandline:
    ### ~ INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    *See RJE_SEQ and RJE_UNIPROT commandline options. Must have one or other loaded with sequences.*
    domft=LIST  : List of UniProt features to treat as Domains [PFAM,DOMAIN]
    ppi=FILE    : Pingu PPI File with Gene, EnsLoci & PPI headers []
    ddi=FILE    : iPFam domain-domain interaction file []
    fam=FILE    : File of GABLAM results for family-based methods []

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_uniprot, rje_zen
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
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('SLiMPID', '0.0', 'June 2008', '2008')
    description = 'Short Linear Motif Protein Interaction Datasets'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
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
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje_seq.__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje_uniprot.__doc__)
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
### SECTION II: SLiMPID Class                                                                                               #
#########################################################################################################################
class SLiMPID(rje.RJE_Object):     
    '''
    SLiMPID interaction dataset generating class. Author: Rich Edwards (2008).

    Info:str
    - Fam = File of GABLAM results for family-based methods []
    - PPI = Pingu PPI File with Gene, EnsLoci & PPI headers []
    - DDI = iPFam domain-domain interaction file []
    
    Opt:boolean

    Stat:numeric

    List:list
    - DomFT = List of UniProt features to treat as Domains [PFAM,DOMAIN]

    Dict:dictionary
    - DDI = Dictionary of {Domain:[Domain list]}
    - Domain = Dictionary of {Domain:[Seqs]}
    - Entry = Dictionary of {ShortName:UniProtEntry}
    - Gene = Dictionary of {ShortName:Gene}
    - PPI = Dictionary of {ShortName:[ShortName list]}
    - Seq = Dictionary of {ShortName:Sequence}

    Obj:RJE_Objects
    - SeqList = SeqList object containing sequence objects.
    - UniProt = UniProt object containing annotated entries.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['PPI','DDI','Fam']
        self.optlist = []
        self.statlist = []
        self.listlist = ['DomFT']
        self.dictlist = ['Domain','DDI','Gene','Entry','PPI','Seq']
        self.objlist = ['SeqList','UniProt']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.list['DomFT'] = ['PFAM','DOMAIN']
        ### Other Attributes ###
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
                self._cmdReadList(cmd,'file',['PPI','DDI','Fam'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def run(self):  ### Main run method controlling primary program flow.
        '''Main run method controlling primary program flow.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] Generate Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.ddi()  # Domain-domain interactions
            self.dpi()  # Domain-protein interactions
            self.fpi()  # Family-protein interactions
            self.ppi()  # Remaining protein-protein interactions
        except: self.errorLog('SLiMPID run() error')
#########################################################################################################################
    def setup(self):    ### Loads data into attributes.
        '''Loads data into attributes.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ UniProt Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            uniprot = self.obj['UniProt'] = rje_uniprot.UniProt(self.log,self.cmd_list)
            uniprot.readUniProt()
            if uniprot.entryNum() > 0:  ### UniProt data loaded. Populate seqlist and domain dictionary.
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F'])
                for entry in uniprot.list['Entry']:
                    seq = entry.obj['Sequence']
                    seqlist.seq.append(entry.obj['Sequence'])
                    name = seq.shortName()
                    self.dict['Entry'][name] = entry
                    self.dict['Seq'][name] = seq
                    for ft in entry.list['Feature']:
                        if ft['Type'] in self.list['DomFT']:
                            try:
                                dom = string.split(ft['Desc'])[0]
                                if dom not in self.dict['Domain']: self.dict['Domain'][dom] = []
                                if name not in self.dict['Domain'][dom]: self.dict['Domain'][dom].append(name)
                            except: self.errorLog('Trouble with %s feature %s' % (name,ft))
            ## ~ [1b] ~ SeqList only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                seqlist = rje_seq.SeqList(self.log,self.cmd_list)
                for seq in seqlist.seq:
                    name = seq.shortName()
                    self.dict['Entry'][name] = None
                    self.dict['Seq'][name] = seq
                    #!# Consider adding loading domains from a table #!#
            ## ~ [1c] ~ Add PPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['PPI']    # Dictionary of ShortName-centred 
            ppi = rje.dataDict(self,self.info['PPI'])
            for hub in ppi:
                if ppi[hub]['EnsLoci'] == '-': continue
                ens = ppi[hub]['EnsLoci']
                if ens not in self.dict['PPI']: self.dict['PPI'][ens] = []
                self.dict['Gene'][ens] = hub
                for gene in string.split(ppi[hub]['PPI'],','):
                    if ppi[gene]['EnsLoci'] == '-': continue
                    if ppi[gene]['EnsLoci'] not in self.dict['PPI'][ens]: self.dict['PPI'][ens].append(ppi[gene]['EnsLoci'])
            ## ~ [1d] ~ Add DDI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['DDI'] = {}
            if self.info['DDI'].lower() not in ['','none']:                    
                data = rje.dataDict(self,self.info['DDI'],mainkeys=['Name1'],datakeys=['Name2'],
                                    headers=['Pfam1','Pfam2','Name1','Name2','Acc1','Acc2','Code1','Code2'],lists=True)
                ## ~ Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                (dx,dtot) = (0.0,len(data))
                self.deBug(data)
                try: rje.sortKeys(data)
                except: self.errorLog('Fuck',quitchoice=True)
                for p1 in rje.sortKeys(data):
                    self.progLog('\r#DDI','Parsing DDI from iPFam: %.1f%%' % (dx/dtot))
                    if p1 not in self.dict['DDI']: self.dict['DDI'][p1] = []
                    for p2 in data[p1]['Name2']:
                        if p2 not in self.dict['DDI']: self.dict['DDI'][p2] = []
                        if p2 not in self.dict['DDI'][p1]: self.dict['DDI'][p1].append(p2)
                        if p1 not in self.dict['DDI'][p2]: self.dict['DDI'][p2].append(p1)
                self.printLog('\r#DDI','Parsing DDI from iPFam: %s domains' % (rje.integerString(dtot)))
            ## ~ [1e] ~ Family data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Fam'] = {}
            if self.info['Fam'].lower() not in ['','none']:                    
                data = rje.dataDict(self,self.info['Fam'],mainkeys=['Qry'],datakeys=['Hit'],lists=True)
                for qry in self.dict['Seq']:
                    self.dict['Fam'][qry] = []
                    if qry in data: self.dict['Fam'][qry] = data[qry]['Hit']
                    elif self.dict['Seq'][qry].info['AccNum'] in data: self.dict['Fam'][qry] = data[self.dict['Seq'][qry].info['AccNum']]['Hit']
                    if qry not in self.dict['Fam'][qry]: self.dict['Fam'][qry].append(qry)
        except: self.errorLog('Problem with SLiMPID.setup()',quitchoice=True)
#########################################################################################################################
    ### <3> ### Dataset Generation Methods                                                                              #
#########################################################################################################################
    def ddi(self):  ### Domain-domain interactions
        '''Domain-domain interactions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ddx = 0
            (dx,dtot) = (0.0,len(self.dict['DDI']))
            if not self.dict['DDI'] or not self.dict['Domain']: return
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dom in rje.sortKeys(self.dict['DDI']):
                self.progLog('\r#DDI','Screening domain-domain interactions: %.1f%%; %s removed' % ((dx/dtot),rje.integerString(ddx))); dx += 100
                if dom not in self.dict['Domain']: self.printLog('#DOM','No sequences with "%s" domains' % dom); continue
                for ddi in self.dict['DDI'][dom]:
                    if ddi not in self.dict['Domain']: continue
                    for hub in self.dict['Domain'][dom]:
                        if hub not in self.dict['PPI']: continue
                        for spoke in self.dict['PPI'][hub][0:]:
                            if spoke in self.dict['Domain'][ddi]: ddx+=1; self.dict['PPI'][hub].remove(spoke)
                    for hub in self.dict['Domain'][ddi]:
                        if hub not in self.dict['PPI']: continue
                        for spoke in self.dict['PPI'][hub][0:]:
                            if spoke in self.dict['Domain'][dom]: ddx+=1; self.dict['PPI'][hub].remove(spoke)
            self.printLog('\r#DDI','Screening domain-domain interactions complete: %s removed.' % (rje.integerString(ddx)))
            ### ~ [3] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hx = len(self.dict['PPI'])
            for hub in rje.sortKeys(self.dict['PPI']):
                if hub and self.dict['PPI'][hub]: continue
                self.dict['PPI'].pop(hub)
                self.printLog('#DDI','No %s interactions left after DDI removed' % hub,screen=False)
            self.printLog('#PPX','%s of %s PPI hubs remain after DDI removed' % (rje.integerString(len(self.dict['PPI'])),rje.integerString(hx)))
        except: self.errorLog('Problem with SLiMPID.ddi()',quitchoice=True)
#########################################################################################################################
    def dpi(self):  ### Domain-protein interactions
        '''Domain-protein interactions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Domain']: return
            outdir = 'SLiMPID_DPI'
            rje.mkDir(self,outdir)
            dpi = {}            # Dictionary of {domain:[interactors]}
            badname = []
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dom in rje.sortKeys(self.dict['Domain']):
                dpi[dom] = []
                for hub in self.dict['Domain'][dom]:
                    if hub in self.dict['PPI']: dpi[dom] += self.dict['PPI'][hub]      # Add with redundancy
                for spoke in dpi[dom][0:]:
                    if dpi[dom].count(spoke) == 1: dpi[dom].remove(spoke)   # Must have 2+ domain interactions
                for hub in self.dict['Domain'][dom]:
                    if hub not in self.dict['PPI']: continue
                    for spoke in self.dict['PPI'][hub][0:]:
                        if spoke in dpi[dom]:
                            self.dict['PPI'][hub].remove(spoke)
                            if spoke in self.dict['PPI'] and hub in self.dict['PPI'][spoke]: self.dict['PPI'][spoke].remove(hub)
                dpi[dom] = rje.sortUnique(dpi[dom],False,False)
                acc = []
                for name in dpi[dom]:
                    if not name: continue
                    if name in self.dict['Seq']: acc.append(self.dict['Seq'][name].info['AccNum'])
                    elif name not in badname: badname.append(name) 
                open('%s/%s.dpi.acc' % (outdir,dom),'w').write(string.join(acc,'\n'))
                self.printLog('#DPI','%s domain => %d interactors' % (dom,len(acc)))
            if badname:
                badname.sort()
                self.printLog('#BAD','%d "bad" protein names: %s' % (len(badname),string.join(badname,'; ')))
            ### ~ [3] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hx = len(self.dict['PPI'])
            for hub in rje.sortKeys(self.dict['PPI']):
                if hub and self.dict['PPI'][hub]: continue
                self.dict['PPI'].pop(hub)
                self.printLog('#DPI','No %s PPI left after DPI removed' % hub,screen=False)
            self.printLog('#PPX','%s of %s PPI hubs remain after DPI removed' % (rje.integerString(len(self.dict['PPI'])),rje.integerString(hx)))
        except: self.errorLog('Problem with SLiMPID.dpi()',quitchoice=True)
#########################################################################################################################
    def fpi(self):  ### Family-protein interactions
        '''Family-protein interactions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Domain']: return
            outdir = 'SLiMPID_FPI'
            rje.mkDir(self,outdir)
            fpi = {}            # Dictionary of {family:[interactors]}
            badname = []
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qry in rje.sortKeys(self.dict['PPI']):
                try:
                    fam = self.dict['Fam'][qry]
                    if len(fam) < 2: continue
                except: self.errorLog('Problem with "%s" protein family' % qry); continue
                fpi[qry] = []
                for hub in fam:
                    if hub not in self.dict['PPI']: continue
                    fpi[qry] += self.dict['PPI'][hub]      # Add with redundancy
                for spoke in fpi[qry][0:]:
                    if fpi[qry].count(spoke) == 1: fpi[qry].remove(spoke)   # Must have 2+ family interactions
                for hub in fam:
                    if hub not in self.dict['PPI']: continue
                    for spoke in self.dict['PPI'][hub][0:]:
                        if spoke in fpi[qry]:
                            self.dict['PPI'][hub].remove(spoke)
                            if spoke in self.dict['PPI'] and hub in self.dict['PPI'][spoke]: self.dict['PPI'][spoke].remove(hub)
                fpi[qry] = rje.sortUnique(fpi[qry],False,False)
                acc = []
                gene = self.dict['Gene'][qry]
                for name in fpi[qry]:
                    if not name: continue
                    if name in self.dict['Seq']: acc.append(self.dict['Seq'][name].info['AccNum'])
                    elif name not in badname: badname.append(name)                     
                open('%s/%s.fpi.acc' % (outdir,gene),'w').write(string.join(acc,'\n'))
                self.printLog('#FPI','%s family => %d interactors' % (gene,len(acc)))
            if badname:
                badname.sort()
                self.printLog('#BAD','%d "bad" protein names: %s' % (len(badname),string.join(badname,'; ')))
            ### ~ [3] Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hx = len(self.dict['PPI'])
            for hub in rje.sortKeys(self.dict['PPI']):
                if hub and self.dict['PPI'][hub]: continue
                self.dict['PPI'].pop(hub)
                self.printLog('#FPI','No %s PPI left after FPI removed' % hub)
            self.printLog('#PPX','%s of %s PPI hubs remain after FPI removed' % (rje.integerString(len(self.dict['PPI'])),rje.integerString(hx)))
        except: self.errorLog('Problem with SLiMPID.fpi()',quitchoice=True)
#########################################################################################################################
    def ppi(self):  ### Remaining protein-protein interactions
        '''Remaining protein-protein interactions.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['PPI']: return
            outdir = 'SLiMPID_PPI'
            rje.mkDir(self,outdir)
            badname = []
            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for hub in rje.sortKeys(self.dict['PPI']):
                gene = self.dict['Gene'][hub]
                acc = []
                for name in self.dict['PPI'][hub]:
                    if not name: continue
                    if name in self.dict['Seq']: acc.append(self.dict['Seq'][name].info['AccNum'])
                    elif name not in badname: badname.append(name)                     
                open('%s/%s.ppi.acc' % (outdir,gene),'w').write(string.join(acc,'\n'))
                self.printLog('#PPI','%s => %d individual interactors' % (gene,len(acc)))
            if badname:
                badname.sort()
                self.printLog('#BAD','%d "bad" protein names: %s' % (len(badname),string.join(badname,'; ')))
        except: self.errorLog('Problem with SLiMPID.setup()',quitchoice=True)
#########################################################################################################################
### End of SECTION II: SLiMPID  Class                                                                                   #
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
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:SLiMPID(mainlog,cmd_list).run()
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
