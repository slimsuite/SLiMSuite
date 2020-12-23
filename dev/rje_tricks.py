#!/usr/bin/python

# Generic Methods
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_tricks
Description:  Random code from Python courses
Version:      0.0
Last Edit:    11/02/20
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
### User
import rje
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
# Add methods here
# os.listdir = files in directory
# spacy = module for dealing with language
# ? get getcharset used to pull out encoding and then read with encoding=encoding

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
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
