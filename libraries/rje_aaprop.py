#!/usr/local/bin/python

# rje_aaprop - AA Property Matrix Module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_aaprop
Description:  AA Property Matrix Module
Version:      0.2.0
Last Edit:    09/04/15
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Takes an amino acid property matrix file and reads into an AAPropMatrix object. Converts in an all by all property
    difference matrix. By default, gaps and Xs will be given null properties (None) unless part of input file.

Commandline:
    aaprop=FILE : Amino Acid property matrix file. [aaprop.txt]
    aagapdif=X  : Property difference given to amino acid vs gap comparisons [5]
    aanulldif=X : Property difference given to amino acid vs null values (e.g. X) [0.5]

Uses general modules: re, string, sys, time
Uses RJE modules: rje
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import re
import string
import os, sys
import time
#############################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#############################################################################################################################
### History
# 0.0 - Initial Working Compilation.
# 0.1 - Fixed some bugs and updated use of general methods
# 0.2.0 - Added $ and ^ to be recognised as gaps. Modified STDOUT.
#############################################################################################################################
### Major Functionality to Add
# [ ] Add a method to read from a property difference matrix file
#############################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_AAPROP'
        version = '0.2.0'
        last_edit = 'April 2015'
        description = 'Amino Acid Property Matrix Module'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#############################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
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
#############################################################################################################################
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
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#############################################################################################################################

#############################################################################################################################
### AAPropMatrix Class: 
#############################################################################################################################
class AAPropMatrix(rje.RJE_Object):     
    '''
    Amino Acid Property Matrix Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of file from which data is loaded
    
    Opt:boolean

    Stat:numeric
    - GapDif = Property difference given to amino acid vs gap comparisons [5]
    - NullDif = Property difference given to amino acid vs null values (e.g. X) [0.5]

    Obj:RJE_Objects
    '''
    ### Attributes
    prop = {}   # Dictionary of properties: keys [Property:str][AA:str]
    pdif = {}   # Dictionary of property differences: key['%s%s' % (aa1,aa2)]
    alphabet = []    # List of single letter codes for use in matrices
#############################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#############################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name']
        - Stats:float ['GapDif','NullDif']
        - Opt:boolean []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name']
        self.statlist = ['GapDif','NullDif']
        self.optlist = []
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.5,obj=None,setlist=False,setdict=False)
        self.info['Name'] = 'aaprop.txt'
        self.stat['GapDif'] = 5
        ### <c> ### Other Attributes
        self.prop = {}
        self.pdif = {}
        self.alphabet = []
#############################################################################################################################
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
                self._cmdRead(type='info',att='Name',arg='aaprop',cmd=cmd)
                self._cmdRead(type='stat',att='GapDif',arg='aagapdif',cmd=cmd)
                self._cmdRead(type='stat',att='NullDif',arg='aanulldif',cmd=cmd)
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        self.readAAProp()
#############################################################################################################################
    ### <2> ### File Input/Output
#############################################################################################################################
    def readAAProp(self,filename=None): ### Reads AA Property Matrix from file
        '''
        Reads AA Property Matrix from file.
        >> filename:str = Filename. If None, will use self.info['Name']
        '''
        try:
            ### <a> ### Load and read
            if filename:
                self.info['Name'] = filename
            else:
                filename = self.info['Name']
            readtxt = 'Reading AA Properties from %s...' % filename
            self.progLog('\r#AAPROP',readtxt)
            proplines = self.loadFromFile(filename,v=2)
            ### <b> ### Process
            self.alphabet = []
            self.prop = {}
            ## <i> ## Properties and alphabet
            for line in proplines:
                line = rje.chomp(line)
                if line.find('#') == 0: # Comment line
                    continue
                elif line.find('PROP') == 0:  # Header line - has amino acids
                    line = rje.matchExp('^\S+(\s.+)',line)[0]
                    while re.search('^\s+\S.*',line):
                        (aa,line) = rje.matchExp('^\s+(\S)(.*)',line)
                        self.alphabet.append(aa)
                    readtxt += ' ...%s' % string.join(self.alphabet)
                    self.progLog('\r#AAPROP',readtxt)
                elif re.search('^\S',line) and self.alphabet:   # Property line
                    (aaproperty,line) = rje.matchExp('^(\S+)(\s.+)',line)
                    readtxt += ' ...%s' % aaproperty
                    self.progLog('\r#AAPROP',readtxt)
                    self.prop[aaproperty] = {}
                    for aa in self.alphabet:
                        (p,line) = rje.matchExp('^\s+(\S)(.*)',line)
                        self.prop[aaproperty][aa] = p
                    #self.verbose(2,3,'...%s' % self.prop[property],0)
            readtxt += ' ...Done!'
            self.printLog('\r#AAPROP',readtxt)
        except IOError:
            self.log.errorLog('AA Property matrix file %s missing?' % self.info['Name'],True)
            raise
        except:
            self.log.errorLog('Major Problem reading AA Property matrix(%s)' % self.info['Name'],True)
            return
        add = []
        if 'X' not in self.alphabet:
            add.append('X')
        if '-' not in self.alphabet:
            add.append('-')
        if add:
            add = self.alphabet + add
            self.useAlphabet(alphabet=add)
        self.makePropDif()
#############################################################################################################################
    def saveAAProp(self,filename=None): ### Saves AA Property Matrix to file
        '''
        Saves AA Property Matrix to file.
        >> filename:str = Filename. If None, will use self.info['Name']
        '''
        try:
            ### <a> ### Open
            if filename:
                self.info['Name'] = filename
            else:
                filename = self.info['Name']
            self.verbose(0,3,'Writing AA Properties to %s...' % filename,0)
            AAPROP = open(filename,'w')
            ### <b> ### Write
            ## <i> ## Headers
            AAPROP.write('# File created by rje_aaprop.py\n')
            maxlen = 8
            for aaproperty in self.prop.keys():
                if len(aaproperty) > maxlen:
                    maxlen = len(aaproperty)
            header = 'PROPERTY'
            while len(header) < maxlen:
                header = '%s ' % header
            alph = self.alphabet[0:]
            alph.sort()
            for aa in alph:
                header = '%s %s' % (header,aa)
            AAPROP.write('%s \n' % header)
            ## <ii> ## Properties
            pkey = self.prop.keys()
            pkey.sort()
            for aaproperty in pkey:
                self.verbose(1,3,'...%s' % aaproperty,0)
                text = aaproperty
                while len(text) < maxlen:
                    text = '%s ' % text
                for aa in alph:
                    text = '%s %s' % (text,str(self.prop[aaproperty][aa]))
                AAPROP.write('%s \n' % text)
            AAPROP.close()
            self.verbose(0,2,'...Done!',1)
        except IOError:
            self.log.errorLog('Cannot create AA Property matrix file %s?!' % self.info['Name'])
            return
        except:
            self.log.errorLog('Major Problem saving AA Property matrix(%s)' % self.info['Name'])
            return
#############################################################################################################################
    def savePropDif(self,filename='aapdif.txt'): ### Saves AA Property Difference Matrix to file
        '''
        Saves AA Property Difference Matrix to file.
        >> filename:str = Filename. If None, will use self.info['Name']
        '''
        try:
            ### <a> ### Open
            if filename:
                self.info['Name'] = filename
            else:
                filename = self.info['Name']
            self.verbose(0,3,'Writing AA Property Differences to %s...' % filename,0)
            AAPROP = open(filename,'w')
            ### <b> ### Write
            ## <i> ## Headers
            AAPROP.write('# File created by rje_aaprop.py\n')
            header = 'AA'
            alph = self.alphabet[0:]
            alph.sort()
            for aa in alph:
                header = '%s %s' % (header,aa)
            AAPROP.write('%s \n' % header)
            ## <ii> ## Property Differences
            for a1 in alph:
                self.verbose(1,3,'...%s' % a1,0)
                text = '%s ' % a1
                for a2 in alph:
                    text = '%s %s' % (text,str(self.pdif['%s%s' % (a1,a2)]))
                AAPROP.write('%s \n' % text)
            AAPROP.close()
            self.verbose(0,2,'...Done!',1)
        except IOError:
            self.log.errorLog('Cannot create AA Property Difference matrix file %s?!' % self.info['Name'])
            return
        except:
            self.log.errorLog('Major Problem saving AA Property Differencematrix(%s)' % self.info['Name'])
            return
#############################################################################################################################
    ### <3> ### Property Matrix Maniuplations
#############################################################################################################################
    def makePropDif(self):  ### Converts the property matrix into a property difference matrix
        '''
        Converts the property matrix into a property difference matrix.
        '''
        try:
            ### <a> ### Setup
            self.pdif = {}
            ### <b> ### Convert
            for a1 in self.alphabet:
                for a2 in self.alphabet[self.alphabet.index(a1):]:
                    key1 = '%s%s' % (a1,a2)
                    key2 = '%s%s' % (a2,a1)
                    if a1 in '-^$' or a2 in '-^$':
                        self.pdif[key1] = self.stat['GapDif']
                    else:
                        self.pdif[key1] = 0.0
                        for propdic in self.prop.values():
                            if propdic[a1] == None or propdic[a2] == None:
                                self.pdif[key1] += self.stat['NullDif']
                            elif propdic[a1] == propdic[a2]:
                                self.pdif[key1] += 0
                            else:
                                self.pdif[key1] += 1
                        self.pdif[key1] = int(self.pdif[key1]+0.5)
                    self.pdif[key2] = self.pdif[key1]
        except:
            self.log.errorLog('Major Problem with makePropDif().',True)
#############################################################################################################################
    def useAlphabet(self,alphabet,missing=None,trim=False): ### Makes sure matrix is using supplied alphabet (no missing)
        '''
        Makes sure matrix is using supplied alphabet (no missing values to cause errors).
        >> alphabet:list = single-letter codes to ues.
        >> missing:num = values to give properties currently missing a given letter.
        >> trim:boolean = whether to delete parts of the property matrix that are not in given alphabet
        '''
        try:
            ### <a> ### Setup
            newaa = []      # List of new letters to add
            remaa = []      # List of letters to remove
            ### <b> ### Compare alphabets
            for letter in self.alphabet:
                if trim and letter not in alphabet:
                    remaa.append(letter)
                    self.alphabet.remove(letter)
            for letter in alphabet:
                if letter not in self.alphabet:
                    newaa.append(letter)
                    self.alphabet.append(letter)
            if len(newaa) == 0 and len(remaa) == 0:     # No changes
                return
            ### <c> ### Make changes
            ## <i> ## Update self.prop
            for propdic in self.prop.values():
                for aa in propdic.keys():
                    if trim and aa in remaa:
                        propdic.pop(aa)
                for aa in newaa:
                    propdic[aa] = missing
            ## <ii> ## Update self.pdif
            self.makePropDif()
        except:
            self.log.errorLog('Major Problem during useAlphabet().')
#############################################################################################################################
### End of AAPropmatrix Class
#############################################################################################################################
    
#############################################################################################################################
## End of SECTION II
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                           #
#############################################################################################################################

### Define module-specific methods here
    
#############################################################################################################################
### END OF SECTION III
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                               #
#############################################################################################################################
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
        aaprop = AAPropMatrix(log=mainlog, cmd_list=cmd_list)
        print 'Not for Standalone running.'        

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
### END OF SECTION IV
#############################################################################################################################
