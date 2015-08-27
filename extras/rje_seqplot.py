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
Module:       rje_seqplot
Description:  Sequence plotting module
Version:      0.0
Last Edit:    16/06/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to work in conjunction with Norman's sequence plotting tool(s).

    See also rje_seq, rje_uniprot and rje_disorder options.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : File containing input sequences [None]
    occfile=FILE    : File containing motif occurrences to plot [None]
    plotstat=LIST   : List of stats to plot (cons/rel/dis/sa/hyd) [cons/rel/dis]
    plotft=LIST     : List of features to plot if input is UniProt [SIGNALP,TRANSMEMBRANE,DOMAIN,PFAM]
    plotre=LIST     : List of regular expressions to plot []
    plotwin=X       : Window for stats plot +/- [7]
    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    outfile=X       : Leader for output files [None]
    occonly=T/F     : Only output sequences with motif occurrences [False]
    rescale=T/F     : Whether to rescale plotstats (Truncate at 0.0 and normalise to max 1.0) [True]
    
Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_disorder, rje_seq, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_disorder, rje_seq, rje_slimcalc, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : combine=T/F     : Whether to combine results from different datasets into one file per seq [False]
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_SEQPLOT', '0.0', 'June 2008', '2008')
    description = 'Sequence plotting module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',
                'This program is designed to be used with the plotting webservers of Norman Davey']
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
            if rje.yesNo('Show rje_seq commandline options?'): out.verbose(-1,4,text=rje_seq.__doc__)
            if rje.yesNo('Show rje_disorder commandline options?'): out.verbose(-1,4,text=rje_disorder.__doc__)
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
### SECTION II: SeqPlotter Class                                                                                        #
#########################################################################################################################
class SeqPlotter(rje.RJE_Object):     
    '''
    SeqPlotter Class. Author: Rich Edwards (2008).

    Info:str
    - OccFile = File containing motif occurrences to plot [None]
    
    Opt:boolean
    - OccOnly = Only output sequences with motif occurrences [False]
    - ReScale = Whether to rescale plotstats (Truncate at 0.0 and normalise to max 1.0) [True]

    Stat:numeric
    - PlotWin = Window for stats plot +/- [7]
    
    List:list
    - PlotStat = List of stats to plot (cons/rel/dis/sa/hyd) [cons,rel,dis]
    - PlotFT = List of features to plot if input is UniProt [SIGNALP,TRANSMEMBRANE,DOMAIN,PFAM]
    - PlotRE = List of regular expressions to plot []
    
    Dict:dictionary
    - OccData = Data loaded from Occ dictionary [Dataset,Rank,Pattern,Seq,Start_Pos,End_Pos]

    Obj:RJE_Objects
    - SeqList = SeqList object containing sequences
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['OccFile']
        self.optlist = ['OccOnly','ReScale']
        self.statlist = ['PlotWin']
        self.listlist = ['PlotStat','PlotFT','PlotRE']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setStat({'PlotWin':7})
        self.setOpt({'ReScale':True})
        self.list['PlotStat'] = ['cons','rel','dis']
        self.list['PlotFT'] = ['SIGNALP','TRANSMEMBRANE','DOMAIN','PFAM']
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
                self._cmdReadList(cmd,'info',['OccFile'])
                self._cmdReadList(cmd,'opt',['OccOnly','ReScale'])
                self._cmdReadList(cmd,'stat',['PlotWin'])
                self._cmdReadList(cmd,'list',['PlotStat','PlotFT','PlotRE'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = ''
            elif self.info['Basefile'][-1] != '.': self.info['Basefile'] += '.'
            self.obj['SeqList'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T'])
            self.list['PlotFT'] = string.split(string.join(self.list['PlotFT']).upper())
            if self.info['OccFile'].lower() not in ['','none']:
                self.info['Delimit'] = rje.delimitFromExt(filename=self.info['OccFile'])
                self.dict['OccData'] = {}
                occdata = rje.dataDict(self,self.info['OccFile'],['Seq','Dataset','Pattern','Start_Pos','End_Pos'],['Seq','Dataset','Pattern','Start_Pos','End_Pos'])
                for key in rje.sortKeys(occdata):
                    seq = occdata[key].pop('Seq')
                    if seq not in self.dict['OccData']: self.dict['OccData'][seq] = {}
                    dataset = occdata[key].pop('Dataset')
                    if dataset not in self.dict['OccData'][seq]: self.dict['OccData'][seq][dataset] = []
                    self.dict['OccData'][seq][dataset].append(occdata[key])
                self.printLog('#OCC','Loaded data for %s occurrences in %s sequences' % (rje.integerString(len(occdata)),rje.integerString(len(self.dict['OccData']))))
                self.obj['SeqList'].autoFilter(['GoodSeq=%s' % string.join(rje.sortKeys(self.dict['OccData']),',')])
            ### ~ [2] Calculate Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['PlotStat'] = string.split(string.join(self.list['PlotStat']).lower())
            if 'cons' in self.list['PlotStat'] or 'rel' in self.list['PlotStat']: slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list)
            seqdict = self.obj['SeqList'].seqNameDic()
            for name in rje.sortKeys(seqdict):
                if self.opt['OccOnly'] and not name in self.dict['OccData']: continue
                seq = seqdict[name]
                sequence = seq.getSequence(gaps=False)
                seq.dict['PlotStat'] = {}
                if 'sa' in self.list['PlotStat']: seq.dict['PlotStat']['SA'] = rje_seq.surfaceAccessibility(sequence,returnlist=True)
                if 'hyd' in self.list['PlotStat']: seq.dict['PlotStat']['Hydropathy'] = rje_seq.eisenbergHydropathy(sequence,returnlist=True)
                if 'dis' in self.list['PlotStat']: seq.dict['PlotStat']['Disorder'] = seq.disorder(returnlist=True)
                if 'cons' in self.list['PlotStat'] or 'rel' in self.list['PlotStat']:
                    slimcalc.relConListFromSeq(seq,slimcalc.stat['RelConWin'],store=True)
                    try:
                        seq.dict['PlotStat']['Cons_Abs'] = seq.list.pop('Cons')
                        seq.dict['PlotStat']['Cons_Rel'] = seq.list.pop('RelCons')
                    except: self.printLog('#CONS','No conservation stats for %s' % name)
                self.printLog('#STAT','PlotStats calculated for %s' % name)
                for stat in seq.dict['PlotStat']:
                    if stat != 'Cons_Rel' and self.stat['PlotWin'] >= 0: seq.dict['PlotStat'][stat] = self.plotWin(seq.dict['PlotStat'][stat])
                    seq.dict['PlotStat'][stat] = self.convertStat(seq.dict['PlotStat'][stat])
                self.printLog('#STAT','PlotStats converted for %s' % name)                
            ### ~ [3] Output Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if name in self.dict['OccData']:
                    for dataset in self.dict['OccData'][name]:
                        ofile = '%s%s.%s.plot.txt' % (self.info['Basefile'],dataset,seq.info['AccNum'])
                        self.output(seq,ofile,self.dict['OccData'][name][dataset])
                else: self.output(seq,'%s%s.plot.txt' % (self.info['Basefile'],seq.info['AccNum']))
            return
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def plotWin(self,plotstat): ### Returns window conversion
        '''Returns window conversion.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            winstat = []
            ### ~ [2] ~ Convert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(plotstat)):
                win = plotstat[max(0,i-self.stat['PlotWin']):i+self.stat['PlotWin']+1]
                winstat.append(float(sum(win))/len(win))
            return winstat
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def convertStat(self,plotstat):     ### Returns text/scaled conversion
        '''Text/scaled conversion.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            convstat = []            
            for i in range(len(plotstat)):
                v = plotstat[i]
                if self.opt['ReScale']: v = max(0,v/max(plotstat))
                convstat.append('%.5f' % v)
            return convstat
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def output(self,seq,outfile,occdata=[]):   ### Output to file
        '''Output to file.'''
        try:### ~ [1] ~ Basic Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['OccOnly'] and not occdata: return
            odata = ['Name\t%s' % (seq.shortName()),
                     'Sequence\t%s' % (seq.getSequence(gaps=False)),
                     'Output\t%s' % (string.join(string.split(outfile,'.')[:-1],'.')),
                     'RE\t%s' % (string.join(self.list['PlotRE'],',')),
                     'TrueELMs\tY',
                     'Description\t%s' % (seq.info['Description']),
                     '',]
            ### ~ [2] ~ PlotStats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for plot in rje.sortKeys(seq.dict['PlotStat']): odata.append('Plot\t%s\t%s' % (plot,string.join(seq.dict['PlotStat'][plot],', ')))
            if seq.dict['PlotStat']: odata.append('')
            ### ~ [3] ~ PlotFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq.obj['Entry']:
                for ft in seq.obj['Entry'].list['Feature']:
                    if ft['Type'] in self.list['PlotFT']: odata.append('Region\t%s %s\t%s:%s' % (ft['Type'],ft['Desc'],ft['Start'],ft['End']))
                odata.append('')
            ### ~ [4] ~ MotifOcc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if occdata:
                for occ in occdata: odata.append('Motif\t%s\t%s:%s' % (occ['Pattern'],occ['Start_Pos'],occ['End_Pos']))
            ### ~ [5] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            open(outfile,'w').write(string.join(odata,'\n'))
            self.printLog('#PLOT','SeqPlot output saved as %s' % (outfile))
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: SeqPlotter Class                                                                                 #
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
    try:SeqPlotter(mainlog,cmd_list).run()
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
