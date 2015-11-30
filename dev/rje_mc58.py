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
Module:       rje_mc58
Description:  Custom analysis pipeline for MC58 data
Version:      0.0
Last Edit:    12/06/08
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    To be added.

Commandline:

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
import rje, rje_seq, rje_zen
import gablam
import rje_blast_V1 as rje_blast
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
    (program, version, last_edit, copyright) = ('RJE_MC58', '0.0', 'June 2008', '2008')
    description = 'Custom analysis pipeline for MC58 data'
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
### SECTION II: MC58 Class                                                                                               #
#########################################################################################################################
class MC58(rje.RJE_Object):     
    '''
    MC58 Class. Author: Rich Edwards (2008).

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
        ### Basics ###
        self.infolist = []
        self.optlist = []
        self.statlist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
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
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def addLoadObj(self,object_type=None,log=None,cmd_list=[]):  ### Returns a new object of the right type
        '''
        Returns a new object of the right type.
        >> object_type:str = Code from object used to identify correct object to create.
        >> log:Log object to feed new object
        >> cmd_list:List of commands for new object
        '''
        try:
            return None     #!# Edit this subroutine for the class #!#
        except:
            self.log.errorLog('Oh, the shame! Something has gone wrong in addLoadObj()')
            return None
#########################################################################################################################
    def selfLoadTidy(self,obj_dict):    ### Tidies up unusual data types, such as dictionaries
        '''
        Tidies up unusual data types, such as dictionaries.
        >> obj_dict:Dictionary of object codes and objects
        '''
        try:
            #!# Add class-specific code #!#
            return True                        
        except:
            self.log.errorLog('Oh dear! Major problem in selfLoadTidy()')
            return False
#########################################################################################################################
    def loadSpecial(self,line=''):  ### Loads a single piece of data 
        '''
        Loads a single piece of data.
        >> line:str = Special data to load
        '''
        try:
            ### Setup ###
            #!# Replace this method with the appropriate one for 'Special:::Data' input #!#
            self.log.errorLog('%s cannot handle Special input: %s.' % (self.info['Name'],line),False,False)
            return False
        except:
            self.log.errorLog('Major problem in loadSpecial(%s)' % line)
            return False
#########################################################################################################################
    ### <2> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] Reformat Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for fasta in glob.glob('*.fasta'):
                fas = fasta[:-2]
                if os.path.exists(fas): continue
                sx = 0
                for line in open(fasta,'r').readlines():
                    if line[:1] == '>':
                        try: (name,desc) = rje.matchExp('^>(\S+) (\S.+)$',line)
                        except: name = rje.matchExp('^>(\S+)',line)[0]
                        if len(string.split(name,'|')) == 3:
                            name = '6rf_NEIME__%s' % string.split(name,'|')[2]
                            open(fas,'a').write('>%s\n' % name)
                        elif len(string.split(name,'|')) == 5:
                            name = 'ref_NEIME__%s' % string.split(name,'|')[3]
                            open(fas,'a').write('>%s %s\n' % (name,desc))
                        else: print string.split(name,'|'); raise ValueError
                        self.progLog('\r#FAS','Processing %s: %s seqs' % (fas, rje.integerString(sx))); sx += 1
                    else: open(fas,'a').write(line)
                self.printLog('\r#FAS','Processed %s: %s seqs from %s' % (fas, rje.integerString(sx), fasta))
                rje_blast.BLASTRun(self.log,self.cmd_list).formatDB(fas,protein=True,force=True)
            ### ~ [2] Read in CSV Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rfhits = {}     # Dictionary of {hit:['File:hit_num']}
            acc = 'MC58_6RF_Hits.acc'; open(acc,'w')
            gfile = 'MC58_6RF_Hits.vs.MC58_1.hitsum.tdt'
            cx = 0
            for csv in glob.glob('MC58_6RF_CSV/*.CSV'):
                cx += 1
                file = os.path.basename(csv)[:-4]
                hits = False
                for line in open(csv,'r').readlines():
                    if line.find('prot_hit_num,prot_acc') == 0: hits = True
                    elif hits:
                        data = rje.readDelimit(line,',')
                        if len(data) < 2: continue
                        [num,name] = data[:2]
                        try: name = string.split(name,'|')[2]
                        except: continue
                        if name not in rfhits:
                            open(acc,'a').write('6rf_NEIME__%s\n' % name)
                            rfhits[name] = []
                        id = '%s:%s' % (file,num)
                        if id not in rfhits[name]: rfhits[name].append(id)
                        self.progLog('\r#CSV','Reading %d CSV files: %s 6RF Hits' % (cx,rje.integerString(len(rfhits))))
            self.printLog('\r#CSV','Read %d CSV files: %s 6RF Hits output to %s' % (cx,rje.integerString(len(rfhits)),acc))
            ### ~ [3] Extract sequences and perform GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.exists(gfile):
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % acc,'fasdb=MC58_6RF.fas','seqout=MC58_6RF_Hits.fas','autoload=T','accnr=F','seqnr=F'])
                seqlist.info['Name'] = 'MC58_6RF_Hits.fas'
                seqlist.saveFasta()
                gablam.GABLAM(self.log,self.cmd_list+['seqin=MC58_6RF_Hits.fas','searchdb=MC58_1.fas','qryacc=F']).gablam()
            ### ~ [4] Read in GABLAM and ID Hits without genomic homology ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdata = rje.dataDict(self,gfile,['Qry'],['HitNum'])
            zeros = []
            for hit in gdata:
                if string.atoi(gdata[hit]['HitNum']) == 0: zeros.append(hit)
            zeros = rje.sortUnique(zeros,False)
            open('6rf_zeros.acc','w').write(string.join(zeros,'\n'))
            self.printLog('#ZERO','%d 6RF hits with 0 BLAST hits to MC58_1' % len(zeros))
            ufile = 'MC58_6RF_Zeros.vs.embl_bacteria.hitsum.tdt'
            if not os.path.exists(ufile):
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=6rf_zeros.acc','fasdb=MC58_6RF.fas','seqout=MC58_6RF_Zeros.fas','autoload=T','accnr=F','seqnr=F'])
                seqlist.info['Name'] = 'MC58_6RF_Zeros.fas'
                seqlist.saveFasta()
                gablam.GABLAM(self.log,self.cmd_list+['seqin=MC58_6RF_Zeros.fas','searchdb=/scratch/Databases/NewDB/TaxaDB/embl_bacteria.fas','qryacc=F']).gablam()
            gdata = rje.dataDict(self,ufile,['Qry'],getheaders=True)
            fdata = rje.dataDict(self,string.replace(ufile,'hitsum','gablam'),['Qry'],['Hit'],lists=True)
            headers = gdata.pop('Headers')
            headers.insert(1,'Sample')
            headers.append('BestHit')
            rje.delimitedFileOutput(self,'MC58_6RF_Zeros.tdt',headers,rje_backup=True)
            for rf in rje.sortKeys(gdata):
                rfcut = string.split(rf,'__')[1]
                gdata[rf]['Sample'] = string.join(rfhits[rfcut],'; ')
                gdata[rf]['Qry'] = rfcut
                try: gdata[rf]['BestHit'] = fdata[rf]['Hit'][0]
                except: gdata[rf]['BestHit']  = '-'
                rje.delimitedFileOutput(self,'MC58_6RF_Zeros.tdt',headers,datadict=gdata[rf])
            
        except: self.errorLog(rje_zen.Zen().wisdom())
        self.printLog('#ZEN',rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: MC58 Class                                                                                       #
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
    try:MC58(mainlog,cmd_list).run()
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
