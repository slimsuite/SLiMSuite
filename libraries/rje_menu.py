#!/usr/local/bin/python

# RJE Menu Methods
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_menu
Description:  Generic Menu Methods Module
Version:      0.4.0
Last Edit:    18/01/16
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to contain generic menu methods for use with any RJE Object. At least, that's the plan...

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module.

Uses general modules: os, string, sys
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added infile and outfile
    # 0.2 - Added extra addcommand and printtext menu options.
    # 0.3 - Modified to work with new object types.
    # 0.4.0 - Changed handling of default for exiting menu loop. May affect behaviour of some existing menus.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add a "Flat" menu option that just lists options in each section.
    '''
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: METHODS                                                                                                 #
#########################################################################################################################
def menu(callobj,headtext='',menulist=[],choicetext='Please select:',changecase=True,default='',jointxt=' = '):   ### Main Menu method
    '''
    Main Menu method.
    >> callobj:Object for which attributes are to be read and altered. Also controls interactivity and log.
    >> headtext:str [''] = Introductory text for menu system.
    >> menulist:list [] = List of menu item tuples (edit code,description,optiontype,optionkey)
        - e.g. ('0','Sequence file','info','Name') would edit callobj.info['Name'])
        - If optiontype == 'return' then menu will return the value given in optionkey
        - If optiontype == '' then description will be printed as a breaker
        - If optiontype == 'infile' then callobj.info['Name'] would be changed using rje.getFileName(mustexist=True)
        - If optiontype == 'outfile' then callobj.info['Name'] would be changed using rje.getFileName(confirm=True)
        - If optiontype == 'showtext' then optionkey should contain text to be printed with verbose
        - If optiontype == 'addcmd' then commands can be added.
    >> choicetext:str ['Please select:'] = Text to display for choice option
    >> changecase:boolean [True] = change all choices and codes to upper text
    >> default:str [''] = What to return if nothing selected.
    >> jointxt:str [' = '] = What to join code and description with when listing options.
    << returns optionkey if appropriate, else True
    '''
    try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ## ~ [0a] Choice Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        choicedict = {}
        for (code,desc,vtype,key) in menulist:
            if not vtype: continue
            if changecase: choicedict[code.upper()] = (vtype,key)
            else: choicedict[code] = (vtype,key)
        ## ~ [0b] Setup Header Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        maxlen = 0
        for line in string.split(headtext,'\n'):
            if len(line) > maxlen: maxlen = len(line)
        headlist = ['#' * (maxlen + 10)]
        for line in string.split(headtext,'\n')[0:]:
            while len(line) < maxlen: line += ' '
            headlist.append('# #> %s <# #' % line)
        headlist.append(headlist[0])
        ### ~ [1] Main Menu Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while menulist:
            ## ~ [1a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mtxt = '\n%s' % string.join(headlist,'\n')
            while mtxt[-2:] != '\n\n': mtxt += '\n'
            for (code,desc,vtype,key) in menulist:
                if vtype and (code or desc):
                    if code and desc: mtxt += '%s%s%s' % (code,jointxt,desc)
                    elif code: mtxt += code
                    elif desc: mtxt += desc
                    if vtype in ['info','list','opt','stat','int','str','bool','num']: mtxt += ': %s' % callobj.getAttribute(vtype,key,default='#!#ERROR#!#')
                    elif vtype in ['infile','outfile']: mtxt += ': %s' % callobj.getAttribute('info',key,default='#!#ERROR#!#')
                else: mtxt += desc
                mtxt += '\n'
            ## ~ [1b] Give Choices ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            print mtxt
            while mtxt:
                try:## ~ Input user choice ~~~ ##
                    choice = rje.choice(choicetext,default=default)
                    if changecase: choice = choice.upper()
                    ## ~ Process user choice ~ ##
                    if choicedict.has_key(choice):
                        (vtype,key) = choicedict[choice]
                        if vtype in ['str','info']: callobj.setInfo({key:callobj._editChoice(key,callobj.getStr(key))})
                        if vtype in ['num','stat']: callobj.setStat({key:callobj._editChoice(key,callobj.getNum(key),numeric=True)})
                        if vtype == 'int': callobj.setStat({key:int(callobj._editChoice(key,callobj.getInt(key),numeric=True))})
                        if vtype in ['bool','opt']: callobj.setOpt({key: not callobj.getBool(key)})
                        if vtype == 'list': callobj.list[key] = string.split(callobj._editChoice(key,callobj.list[key]))
                        if vtype == 'infile': callobj.setInfo({key: rje.getFileName('%s File Name?' % key,callobj.getStr(key))})
                        if vtype == 'outfile': callobj.setInfo({key: rje.getFileName('%s File Name?' % key,callobj.getStr(key),mustexist=False,confirm=True)})
                        if vtype == 'showtext': callobj.verbose(-1,-1,key); break
                        if vtype == 'addcmd':
                            prevcmd = callobj.cmd_list
                            callobj.cmd_list = rje.inputCmds(out,choice)
                            callobj.printLog('#CMD','User Added commands: %s' % callobj.cmd_list)
                            callobj._cmdList()
                            callobj.cmd_list = prevcmd + callobj.cmd_list
                            break
                        if vtype in ['info','list','opt','stat','infile','outfile','str','bool','int','num']:
                            callobj.printLog('#%s' % vtype.upper(),'User edited %s parameter' % key); break
                        elif vtype == 'return': return key
                    print 'Choice "%s" not recognised!\n' % choice
                except KeyboardInterrupt:
                    if rje.yesNo('Terminate program?'): raise
                    if rje.yesNo('Exit menu and proceed?'): return default
                except: raise
        ### End ###
        return True
    except KeyboardInterrupt: raise
    except:
        if callobj: callobj.errorLog('Major disaster in rje_menu.menu()',quitchoice=True)
        else: raise
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: 'MAIN' PROGRAM                                                                                         #
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: print 'This module is not for standalone running.'
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
