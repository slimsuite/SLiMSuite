#!/usr/bin/python

# Generic Methods
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_examples
Description:  Example python script for some basic analysis
Version:      0.0
Last Edit:    22/01/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a stripped down template for methods only. This is for when a class has too many methods and becomes
    untidy. In this case, methods can be moved into a methods module and 'self' replaced with the relevant object.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module.

Uses general modules: copy, os, string, sys
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### Userimport rje
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: METHODS                                                                                                 #
#########################################################################################################################
def powX(x, power=2, stdout=True):   # Three parameters but stdout is now optional
    '''
    This method returns the value of x raised to a given power. If power is not
    given, x will be squared. Returns -1 in case of error.
    >> x:num = base value
    >> power:num [2] = power to raise x to
    >> stdout:bool [True] = whether to print the results to screen
    '''    # It is good to have descriptions of parameters in the document string.
    try:
        y = x ** power
        if stdout: print '%s to the power %s = %s' % (x, power, y)
        return y
    except:
        print 'powX error:', sys.exc_info()[0]
        print 'See correct usage below:'
        print powX.__doc__
        return -1
#########################################################################################################################
def loadWormSeq(seqfile='D:\\Teaching\\2008-08 Programming Course\\caeel_cdna.fas'):  ### Loads sequences from a fasta file and returns a dictionary
    '''Loads sequences from a fasta file and returns a dictionary.'''
    ### Load ###
    seqlines = open(seqfile,'r').readlines()[:20000]    # Temp truncation to speed up read time
    print '%d lines read from %s' % (len(seqlines),seqfile)
    ### Process ###
    seqdict = {}                    # Dictionary of {name:sequence}
    name = ''                       # Name of sequence currently being read
    seq = ''                        # DNA Sequence currently being read
    while seqlines:                 # Loop until seqlines has been fully processed
        line = seqlines.pop(0)              # Remove and return first item of seqlines
        line = string.replace(line,'\n','') # Remove line ending \n
        if line[:1] == '>':                 # Look for start of new sequence
            if name:                        # Must save old sequence if one has been read (i.e. if not first sequence)
                seqdict[name] = seq             # Add previous name and sequence to dictionary
                seq = ''                        # Clear sequence to be read in again for next entry
            name = line[1:]                 # Do not include the '>' in the name
        else:                               # No '>' indicates we are reading sequence
            seq = seq + line                # Add this part of the sequence to what (if anything) we have already read
    if name: seqdict[name] = seq    # The last sequence still needs to be added.
    print len(seqdict), 'sequences read into dictionary'
    ### Return ###
    return seqdict
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: 'MAIN' PROGRAM                                                                                         #
#########################################################################################################################
if __name__ == "__main__":  ### Print message to screen if called from commandline.
    try:
        print __doc__
        # Read in sequences
        # Read in table
        # Make NR
        # Calculate nt freqs
        # Calculate GC content - save for R analysis
        ## >> Upgrade to proper rje class-based approach << ##
        # Calculate dnt freqs
        # Assess for <> expectation
        # Save dnt freqs for analysis in R
        # Analyse GC content vs proximity
    except: print 'Run error:', sys.exc_info()[0]; sys.exit(1)
    sys.exit(0)
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
