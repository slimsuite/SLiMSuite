#!/usr/local/bin/python

# rje_pam - Contains Objects for PAM matrices
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_pam
Description:  Contains Objects for PAM matrices
Version:      1.2
Last Edit:    16/11/10
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module handles functions associated with PAM matrices. A PAM1 matrix is read from the given input file and
    multiplied by itself to give PAM matrices corresponding to greater evolutionary distance. (PAM1 equates to one amino acid
    substitition per 100aa of sequence.) 

Commandline:
    pamfile=X   : Sets PAM1 input file [jones.pam]
    pammax=X    : Initial maximum PAM matrix to generate [100]
    pamcut=X    : Absolute maximum PAM matrix [1000]

Alternative PAM matrix commands:
    altpam=FILE : Alternative to PAM file input = matrix needing scaling by aafreq [None]
    seqin=FILE  : Sequence file from which to calculate AA freq for scaling [None]
    pamout=FILE : Name for rescaled PAM matrix output [*.pam named after altpam=FILE]

Classes:
    PamCtrl(rje.RJE_Object):
        - Controls a set of PAM matrices.
    PAM(pam,rawpamp,alpha):
        - Individual PAM matrix.
        >> pam:int = PAM distance
        >> rawpamp:dictionary of substitution probabilities
        >> alpha:list of amino acids (alphabet)

Uses general modules: os, string, sys, time
Uses RJE modules: rje
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
### History
# 0.0 - Initial Working Compilation.
# 0.1 - Reworking to better OO.
# 0.2 - No Out Object in Objects
# 1.0 - Better documentation to go with GASP V:1.2
# 1.1 - Tidied module up and improved text logging slightly
# 1.2 - Added WAG to PAM conversion
#########################################################################################################################
### Major Functionality to Add
# [ ] : Update Module in line with more recent structure. (Tidy up.)
# [ ] : Add fractional PAM distances. (Continuum.)
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_PAM'
        version = '1.2'
        last_edit = 'November 2010'
        description = 'PAM Matrix Module'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'):
                out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'):
                sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    ### <0> ### Objects setup
    try:
        ## <a> ## Initial Command Setup & Info
        cmd_list = sys.argv[1:]
        info = makeInfo()
        cmd_list = rje.getCmdList(cmd_list,info=info)      ### Load defaults from program.ini
        ## <b> ## Out object
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: PamCtrl Class: Controls a set of PAM matices                                                            #
#########################################################################################################################
class PamCtrl(rje.RJE_Object):     
    '''
    PamCtrl Class. Author: Rich Edwards (2005).
    Controls a set of PAM matrices.

    Info:str
    - Name = Name (filename) [jones.pam]
    - AltPam = Alternative to PAM file input = matrix needing scaling by aafreq [None]
    - SeqIn = Sequence file from which to calculate AA freq for scaling (*Must be FASTA*) [None]
    - PamOut = Name for rescaled PAM matrix output [*.pam named after altpam=FILE]
    
    Opt:boolean
    - X-Value = Whether X is included in matrix probabilities [True]
    - GapValue = Whether - is included in matrix probabilities [True]

    Stat:numeric
    - PamCut = Absolute Upper Limit for PAM [1000]
    - PamMax = Initial Limit for PAM [100]

    Obj:RJE_Objects
    - None

    Other:
    - matrix:list       # List of matrices (PAM Objects)
    - alphabet:list     # List of letters for use in PAM matrices (String)
    '''
    ## Other
    matrix = []       # List of matrices
    alphabet = []     # List of letters for use in PAM matrices
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','Type','AltPam','SeqIn','PamOut']
        - Stats:float ['PamCut','PamMax']
        - Opt:boolean ['X-Value','GapValue']
        - Obj:RJE_Object []
        - Other: matrix, alphabet
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','AltPam','SeqIn','PamOut']
        self.statlist = ['PamCut','PamMax']
        self.optlist = ['X-Value','GapValue']
        self.objlist = []
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=False,setdict=False)
        ### <b> ### Defaults
        self.info['Name'] = 'jones.pam'
        self.stat['PamCut'] = 1000
        self.stat['PamMax'] = 100
        ### <c> ### Other Attributes
        self.matrix = []
        self.alphabet = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='Name',arg='pamfile')
                self._cmdRead(cmd,type='int',att='PamCut')
                self._cmdRead(cmd,type='int',att='PamMax')
                self._cmdRead(cmd,type='info',att='AltPam')
                self._cmdRead(cmd,type='info',att='SeqIn')
                self._cmdRead(cmd,type='info',att='PamOut')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        self.buildPam()
        return
#########################################################################################################################
    ### <2> ### PAM Matrix Construction                                                                                 #
#########################################################################################################################
    def buildPam(self):     ### Builds PAM Matrix in memory
        '''Builds PAM matrix in memory.'''
        try:
            ### Check for Alternative PAM Matrix ###
            if self.info['AltPam'].lower() not in ['','none']:
                self.altPAM()
            
            self.verbose(0,3,"Reading PAM1 matrix from %s" % self.info['Name'],2)
            ### <a> ### Open file & Read Lines
            pamfiles = [self.info['Name'],rje.makePath(self.info['Path']) + self.info['Name'],rje.makePath(self.info['Path']) + rje.makePath('../libraries/') + self.info['Name']]
            self.info['Name'] = None
            for pfile in pamfiles:
                if rje.checkForFile(pfile):
                    file_lines = open(pfile, 'r').readlines()
                    self.info['Name'] = pfile
                    break
            if not self.info['Name']:
                for pfile in pamfiles: self.printLog('#ERR','File "%s" not found' % pfile)
                self.printLog('#ERR','No PAM file found!')
                raise ValueError

            ### <b> ### Read in alphabet
            self.verbose(0,3,file_lines[0],1)
            if file_lines[0].upper().find('X') >= 0:
                self.opt['X-Value'] = False
            if file_lines[0].find('-') >= 0:
                self.opt['GapValue'] = False
            self.alphabet = file_lines[0].split()

            ### <c> ### Make PAM0
            ## <i> ## Clear dics
            zeropamp = {}
            for r in self.alphabet:
                for c in self.alphabet:
                    zeropamp[r + c] = 0
                zeropamp[r + r] = 1
                if self.opt['X-Value']:
                    zeropamp['X' + r] = 1
                    zeropamp[r + 'X'] = 1
                if self.opt['GapValue']:
                    zeropamp['-' + r] = 1
                    zeropamp[r + '-'] = 1
            if self.opt['X-Value']:
                zeropamp['XX'] = 1
            if self.opt['GapValue']:
                zeropamp['--'] = 1
            if self.opt['X-Value'] and self.opt['GapValue']:
                zeropamp['-X'] = 1
                zeropamp['X-'] = 1
            ## <ii> ## New Matrix
            newmatrix = PAM(pam=0,rawpamp=zeropamp,alpha=self.alphabet)
            self.matrix.append(newmatrix)

            ## <d> ## Read in PAM1
            rawpamp = {}
            line = 1
            for r in self.alphabet:
                pamline = file_lines[line].split()
                if len(pamline) != (len(self.alphabet)+1):
                    self.log.errorLog("%s has wrong format! Does not match %s" % (pamline, self.alphabet),printerror=False,quitchoice=True)
                    raise
                for c in range(int(len(self.alphabet))):
                    prob = float(pamline[c+1])
                    rawpamp[r + self.alphabet[c]] = prob
                if self.opt['X-Value']:
                    rawpamp['X' + r] = 1
                    rawpamp[r + 'X'] = 1
                if self.opt['GapValue']:
                    rawpamp['-' + r] = 1
                    rawpamp[r + '-'] = 1
                line += 1
            if self.opt['X-Value']:
                rawpamp['XX'] = 1
            if self.opt['GapValue']:
                rawpamp['--'] = 1
            if self.opt['X-Value'] and self.opt['GapValue']:
                rawpamp['-X'] = 1
                rawpamp['X-'] = 1
            newmatrix = PAM(pam=1,rawpamp=rawpamp,alpha=self.alphabet)
            self.matrix.append(newmatrix)

            ## <e> ## Raise to pammax
            self.log.printLog('\r#PAM','Building PAM Matrices <= %d: ' % self.stat['PamMax'],log=False,newline=False)
            self.pamUp()
            self.log.printLog('\r#PAM','Building PAM Matrices <= %d: Complete.' % self.stat['PamMax'])
        except:
            self.log.errorLog('Fatal Error in PamCtrl.buildPam().')
            raise
#########################################################################################################################
    def pamUp(self):    ### Makes PAM matrix upto self.stat['PamMax']
        '''Makes PAM matrix upto self.stat['PamMax'].'''
        try:
            while len(self.matrix) <= self.stat['PamMax']:
                p = len(self.matrix)
                self.log.printLog('\r#PAM','Building PAM Matrices <= %d: %d' % (self.stat['PamMax'],p),log=False,newline=False)
                newpamp = {}
                for anc in self.alphabet:
                    for desc in self.alphabet:
                        newpamp[anc+desc] = 0
                        for mid in self.alphabet:
                            newpamp[anc+desc] += self.matrix[p-1].pamp[anc+mid] * self.matrix[1].pamp[mid+desc]
                    if self.opt['X-Value']:
                        newpamp[anc + 'X'] = 1
                        newpamp['X' + anc] = 1
                    if self.opt['GapValue']:
                        newpamp[anc + '-'] = 1
                        newpamp['-' + anc] = 1
                if self.opt['X-Value']:
                    newpamp['XX']=1
                if self.opt['GapValue']:
                    newpamp['--']=1
                if self.opt['X-Value'] and self.opt['GapValue']:
                    newpamp['X-']=1
                    newpamp['-X']=1
                newmatrix = PAM(p,newpamp,self.alphabet)
                self.matrix.append(newmatrix)
        except:
            self.log.errorLog('Fatal Error in PamCtrl.PamUp().')
            raise
#########################################################################################################################
    def altPAM(self):   ### Alternative PAM matrix construction
        '''Alternative PAM matrix construction.'''
        try:
            ### Setup ##
            wlines = self.loadFromFile(self.info['AltPam'])
            if not wlines:
                raise IOError
            aas = string.split(wlines[0].upper())
            codes = string.split(wlines[1])
            rawfreqs = string.split(wlines[2])
            freq = {}
            for i in range(len(rawfreqs)):
                freq[aas[i]] = string.atof(rawfreqs[i])
            prob = {}
            for r in range(3,22):
                subs = string.split(wlines[r])
                for i in range(len(subs)):
                    prob['%s%s' % (aas[i],aas[r-2])] = string.atof(subs[i])
                    prob['%s%s' % (aas[r-2],aas[i])] = string.atof(subs[i])

            ### Alternative freqs ###
            if self.info['SeqIn'].lower() not in ['','none'] and os.path.exists(self.info['SeqIn']):
                ## Clear freq ##
                freq = {}
                for a in aas:
                    freq[a] = 0.0
                ## Count freq ##
                slines = self.loadFromFile(self.info['SeqIn'])
                for line in slines:
                    if line[:1] == '>':
                        continue
                    for a in aas:
                        freq[a] += string.count(line.upper(),a)
                ## Convert to freq ##
                total = sum(freq.values())
                if total > 0:
                    for a in aas:
                        freq[a] = freq[a] / total
                self.log.printLog('#AA','Rescaling matrix based on %s aa from %s.' % (rje.integerString(total),self.info['SeqIn']))
            
            ### Calculate s ###
            s = 0.01
            step = 0.000001
            solve = True
            bests = 1.000000
            bestdif = -1
            while solve and s >= step:
                ## Scaler ##
                s = s - step
                self.log.printLog('\r#WAG','Considering s = %.6f; Best s = %.6f (Dif = %.6f)' % (s,bests,bestdif),log=False,newline=False)
                ## Self Subs ##
                newprobs = rje.scaledict(dict=prob,scale=s)
                toobig = False
                for a in aas:
                    newprobs['%s%s' % (a,a)] = 1.0
                    for key in prob.keys():
                        if key[0] == a:
                            newprobs['%s%s' % (a,a)] -= newprobs[key]
                            if newprobs['%s%s' % (a,a)] < 0.0:  # Overshot possibility
                                toobig = True
                                break
                    if toobig:
                        break
                if toobig:
                    continue
                #print 'PAM!!', 
                ## PAM1 ##
                dsum = 0.0
                for a in aas:
                    dsum += freq[a] * newprobs['%s%s' % (a,a)]
                dif = 0.99 - dsum
                if dif < 0:
                    dif = -dif
                if dif < bestdif or bestdif < 0:
                    bestdif = dif
                    bests = s

            ### Output best s ###
            self.log.printLog('\r#WAG','Considered all s <= 0.010000; Best s = %.6f (Dif = %.6f)' % (bests,bestdif))
            if self.info['PamOut'].lower() in ['','none']:
                self.info['PamOut'] = self.info['AltPam'] + '.pam'
            self.log.printLog('#PAM','Rescaled PAM matrix output to %s' % self.info['PamOut'])
            PAM = open(self.info['PamOut'],'w')
            rje.writeDelimit(PAM,aas,' ')
            newprobs = rje.scaledict(dict=prob,scale=bests)
            for a in aas:
                newprobs['%s%s' % (a,a)] = 1.0
                for key in prob.keys():
                    if key[0] == a:
                        newprobs['%s%s' % (a,a)] -= newprobs[key]
            for i in range(len(aas)):
                out = [codes[i]]
                a = aas[i]
                for b in aas:
                    out.append('%.6f' % newprobs['%s%s' % (a,b)])
                rje.writeDelimit(PAM,out,' ')
            PAM.close()
            self.info['Name'] = self.info['PamOut']

        except:
            self.log.errorLog('Major Error with PamCtrl.altPAM().',quitchoice=True)
#########################################################################################################################
    ### <3> ### PAM Distance Calculations                                                                               #
#########################################################################################################################
    def _getPamL(self,pam=0,ancseq='',descseq=''):     ### Returns probability for given PAM and given sequences
        '''
        Returns probability for given PAM and given sequences. Used in PamML.
        >> pam:int = PAM matrix
        >> ancseq:str = ancestral sequence
        >> descseq:str = descendant sequence
        << p:float = probability
        '''
        try:
            p = self.matrix[pam].pamL(ancseq,descseq)
            return p
        except:
            self.log.errorLog('Fatal Error in PamCtrl._getPamL().')
#########################################################################################################################
    def pamML(self,desc='pamML',ancseq='',descseq=''):    ### Calculates PAM distance between two sequences.
        '''
        Calculates PAM distance between two sequences.
        >> desc:str = description
        >> ancseq:str = ancestral sequence
        >> descseq:str = descendant sequence
        << pamml:int = PAM distance
        '''
        try:
            ### <a> ### Setup
            ## <i> ## Checks
            if len(ancseq) != len(descseq):
                self.log.errorLog("%s: Attempting PamML for different length sequences!" % desc,printerror=False)
                raise
            ## <ii> ## Check for zero length branch
            if ancseq == descseq or self._getPamL(0,ancseq,descseq) != 0:
                return 0
            ## <iii> ## Clear Variables
            pamlow = 0
            pamhigh = 0
            pammid = 0
            pamml = 0
            paml = [0] * len(self.matrix)
            # print self.matrix
            
            ### <b> ### Find first high, low and mid PAM
            pamjump = float(self.stat['PamMax'] - 1)/3
            mid1 = int(pamjump)
            mid2 = int(pamjump*2)
            paml[1] = self._getPamL(1,ancseq,descseq)   # matrix[1].PamL(ancseq,descseq)
            paml[mid1] = self._getPamL(mid1,ancseq,descseq)
            paml[mid2] = self._getPamL(mid2,ancseq,descseq)
            paml[self.stat['PamMax']] = self._getPamL(self.stat['PamMax'],ancseq,descseq)
            ## <i> ## 1 <= PamML <= mid1
            if paml[1] >= paml[mid1]:   
                #print "paml[1] >= paml[mid1]"
                pamlow = 1
                pamhigh = mid1
                while pammid == 0:
                    if (pamhigh - pamlow) <= 1: # PamML=1
                        pammid = pamlow
                        pamml = pammid
                    else:
                        mid = int((pamhigh+pamlow)/2)
                        paml[mid] = self._getPamL(mid,ancseq,descseq) # 1/2 between high and low
                        if paml[pamlow] >= paml[mid]:
                            pamhigh = mid   # Still 1 <= PamML <= mid
                        else:
                            pammid = mid    # Now pamlow <= ML <= pamhigh, highest=pammid
                # Now have pamlow <= pammid <= pamhigh
            ## <ii> ## 1 <= ML <= mid2
            elif paml[mid1] >= paml[mid2]:  
                #print "paml[mid1] >= paml[mid2]"
                pamlow = 1
                pammid = mid1
                pamhigh = mid2
                # Now have pamlow <= pammid <= pamhigh
            ## <iii> ## mid1 <= ML <= Max
            elif paml[mid2] >= paml[self.stat['PamMax']]:   
                #print "paml[mid2] >= paml[self.stat['PamMax']]"
                pamlow = mid1
                pammid = mid2
                pamhigh = self.stat['PamMax']
                # Now have pamlow <= pammid <= pamhigh
            ## <iv> ## Need a higher Max!
            else:   
                #print "Need a higher Max!", pamlow, pammid, pamhigh, mid1, mid2, self.stat['PamMax']
                pamlow = mid2
                pammid = self.stat['PamMax']
                while pamhigh == 0:
                    #self.deBug("Upping %d < %d < %d" % (pamlow, pammid, pamhigh))
                    self.stat['PamMax'] += mid1
                    if self.stat['PamMax'] > self.stat['PamCut']:
                        self.stat['PamMax'] = self.stat['PamCut']
                    self.pamUp()
                    while len(paml) <= self.stat['PamMax']:
                        paml.append(0)
                    paml[self.stat['PamMax']] = self._getPamL(self.stat['PamMax'],ancseq,descseq)
                    #print "Upped OK"
                    if paml[pammid] >= paml[self.stat['PamMax']]:
                        pamhigh = self.stat['PamMax']
                    elif self.stat['PamMax'] == self.stat['PamCut']:
                        pamlow = pammid
                        pamhigh = self.stat['PamCut']
                        if (pamhigh - pamlow) <= 1: # 'Found' it = pamcut
                            pammid=pamhigh
                            pamml=pammid
                        else:
                            mid = int((pamhigh + pamlow)/2)
                            paml[mid] = self._getPamL(mid,ancseq,descseq) # 1/2 between high and low
                            if paml[pamlow] >= paml[mid]:
                                pammid = mid    # pamlow <= ML <= pamhigh
                            else:
                                pamlow = mid    # Still mid <= ML <= pamcut
                    else:
                        pamlow = pammid
                        pammid = self.stat['PamMax']  # Need higher pammax still!

            ### <c> ### Should now have low, high and mid... pamlow < pammid < pamhigh
            while pamml==0: # PamML not yet found
                #self.DeBug("Still looking %d < %d < %d" % (pamlow, pammid, pamhigh))
                #print "\n%d >= %d & %d >= %d?" % (pamlow, (pammid-1), pammid ,(pamhigh-1))
                ## <i> ## Found Best PAM!
                if (pamlow >= (pammid-1)) & (pammid >= (pamhigh-1)):    
                    pamml = pammid
                    break
                ## <ii> ## Look between low and mid
                mid1 = int((pammid+pamlow)/2)
                paml[mid1] = self._getPamL(mid1,ancseq,descseq)
                if paml[mid1] > paml[pammid]:
                    pamhigh = pammid
                    pammid = mid1
                    continue
                ## <iii> ## Look between high and mid
                mid2 = int((pammid+pamhigh)/2)
                paml[mid2] = self._getPamL(mid2,ancseq,descseq)
                if paml[mid2] > paml[pammid]:   # Between high and mid
                    pamlow = pammid
                    pammid = mid2
                    continue
                ## <iv> ## Middle is highest
                else:   
                    pamlow = mid1
                    pamhigh = mid2

            ### <d> ### Last catch - not sure when this will apply: problematic (long?) protein sequences?
            if paml[pamml] == 0:
                pamml = self._pamMLSlow(desc,ancseq,descseq)

            ### <e> ### Return best PAM
            #self.DeBug("ML PAM = %d" % pamml)
            return pamml        
        except:
            self.log.errorLog('Fatal Error in PamCtrl.pamML().')
            raise
#########################################################################################################################
    def _pamMLSlow(self,desc,ancseq,descseq):   ### Creeps up through PAM to find best
        '''
        Creeps up through PAM to find best.
        >> desc:str = description
        >> ancseq:str = ancestral sequence
        >> descseq:str = descendant sequence
        << sp:int = PAM distance
        '''
        try:
            ### <a> ### Setup
            sp = 1  # PAM distance to try
            pam0 = 0
            pam1 = self._getPamL(1,ancseq,descseq)
            ### <b> ### Creep until next PAM p is less than current
            while pam1 >= pam0: 
                if sp >= self.stat['PamCut']:   # Max reached - return
                    return self.stat['PamCut']
                if sp == self.stat['PamMax']:
                    self.stat['PamMax'] += 1
                    self.pamUp()
                sp += 1
                pam0 = pam1
                pam1 = self._getPamL(sp,ancseq,descseq)
                if pam1 == 0 and self.stat['Verbose'] > 1:
                    self.log.printLog('#PAM',"%s PAM%d p=0\n" % (desc,sp),0)
            ### <c> ### Return best PAM distance (one before now)
            sp -= 1
            if self.stat['Verbose'] > 1:
                tmp = self._getPamL(sp-1,ancseq,descseq)
                self.log.printLog('#PAM',"%s Finished at %d (%e) vs < %d (%e) > %d (%e)\n" % (desc,sp,pam0,sp-1,tmp,sp+1,pam1),0)
            return sp
        except:
            self.log.errorLog('Fatal Error in PamCtrl._pamMLSlow().')
            raise
#########################################################################################################################
### End of SECTION II: PamCtrl                                                                                          #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: PAM Class: Individual PAM matrix                                                                       #
#########################################################################################################################
class PAM(object):
    ### <0> ### Define the Attributes of the Class
    PAM=0   # PAM value
    def setPam(self,x): self.PAM=x
    def getPam(self): return self.PAM
    pamp={} # Dictionary of probabilities
    def setPamp(self,x): self.pamp=x
    alphabet=[]     # List of letters for use in PAM matrix
#########################################################################################################################
    ### <1> ### Class Initiation: sets attributes
    def __init__(self,pam,rawpamp,alpha):
        '''
        Individual PAM matrix.
        >> pam:int = PAM distance
        >> rawpamp:dictionary of substitution probabilities
        >> alpha:list of amino acids (alphabet)
        '''
        self.setPam(pam)
        self.setPamp(rawpamp)
        self.alphabet=alpha
        self.pamAdjust()
#########################################################################################################################
    ### <2> ### General Class Methods
    def pamCheck(self):
        row={}
        col={}
        for r in self.alphabet:
            row[r]=0
            col[r]=0
        for r in self.alphabet:
            for c in self.alphabet:
                row[r]=row[r]+self.pamp[r+c]
                col[c]=col[c]+self.pamp[r+c]
        for r in self.alphabet:
            print "PAM%d %r: row=%f, col=%f" % (self.PAM, r, row[r], col[r])
#########################################################################################################################
    def pamAdjust(self):    # Makes sure rows total 1
        for r in self.alphabet:
            rowtotal=0
            for c in self.alphabet: rowtotal = rowtotal + self.pamp[r+c]
            for c in self.alphabet: self.pamp[r+c] = self.pamp[r+c] / rowtotal
#########################################################################################################################
    def pamL(self,ancseq,descseq):   # Calculates likelihood of PAM distance between two sequences
        #print "Paml for PAM%d" % self.PAM
        paml=0
        prevpaml=0
        prevsub="r=0"
        for r in range(len(ancseq)):
            sub=ancseq[r]+descseq[r]
            #print sub,
            #print ": p=%f (%e)" % (self.pamp[sub],paml)
            if r==0:
                paml=1e304*self.pamp[sub]  # First residue => Start with very big probability
                prevpaml=paml
            elif paml:
                paml*=self.pamp[sub]
                prevpaml=paml
            else:
                # print "Major problem with PamL, PAM %d residue %d (%s): no likelihood (%e)!" % (self.PAM,r,sub,paml)
                # print "Previous residue: L=%e (%s)" % (prevpaml,prevsub)
                return 0
            prevsub=sub
        #print "return %e" % paml
        return paml
#########################################################################################################################
### END OF SECTION III: PAM Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try:
        [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit:
        return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:
        mypam = PamCtrl(mainlog,cmd_list)
        #X#if mypam.info['AltPam'].lower() not in ['','none']:
        #X#    mypam.buildPam()
        
    ### End ###
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog('User terminated.')
    except:
        mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
