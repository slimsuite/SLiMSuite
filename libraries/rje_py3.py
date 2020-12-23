#!/usr/bin/python

# Python 2.x Methods
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_py3
Description:  Python 3.x versions of core methods
Version:      0.0.0
Last Edit:    21/08/20
Copyright (C) 2020  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a stripped down collection of core methods that have diverged in terms of acceptable Python 2.7 vs
    Python 3.x syntax. See the corresponding rje_py2.py module for the Python 2.7 versions. The main rje.py will import
    one or other of these modules according to the version of python detected.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help').

Uses general modules: copy, os, string, sys
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation with printf() function from rje.py v4.22.5.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Test with updated rje.py.
    '''
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: METHODS                                                                                                 #
#########################################################################################################################
def printf(text,newline=True):    ### Python backwards-compatible print function
    '''Python backwards-compatible print function.'''
    if newline: print(text,flush=True)
    else: print(text, end=' ', flush=True)
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: 'MAIN' PROGRAM                                                                                         #
#########################################################################################################################
if __name__ == "__main__":  ### Print message to screen if called from commandline.
    try:
        print(__doc__)
    except:
        print('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
