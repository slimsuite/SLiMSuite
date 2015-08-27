#!/usr/local/bin/python

# Scoring and Ranking Methods for RJE Python Modules
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_scoring
Description:  Scoring and Ranking Methods for RJE Python Modules
Version:      0.0
Last Edit:    22/01/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains methods only for ranking, filtering and generating new scores from python dictionaries. At its
    conception, this is for unifying and clarifying the new scoring and filtering options used by PRESTO & SLiMPickings,
    though it is conceived that the methods will also be suitable for use in other/future programs.

    The general format of expected data is a list of column headers, on which data may be filtered/ranked etc. or
    combined to make new scores, and a dictionary containing the data for a given entry. The keys for the dictionary
    should match the headers in a *case-insensitive* fashion. (The keys and headers will not be changed but will match
    without using case, so do not have two case-sensitive variables, such as "A" and "a" unless they have the same
    values.) !NB! For some methods, the case should have been matched.

    Methods in this module will either return the input dictionary or list with additional elements (if calculating new
    scores) or take a list of data dictionaries and return a ranked or filtered list.

    Methods in this module:
    * setupStatFilter(callobj,statlist,filterlist) = Makes StatFilter dictionary from statlist and filterlist
    * statFilter(callobj,data,statfilter) = Filters data dictionary according to statfilter dictionary.
    * setupCustomScores(callobj,statlist,scorelist,scoredict) = Checks and returns Custom Scores and related lists
    
Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent modules.

Uses general modules: copy, os, string, sys
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, string, sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
#########################################################################################################################
### Major Functionality to Add
# [ ] : Basic scoring functions from SLiMPickings and PRESTO
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: FILTERING AND SCORING METHODS                                                                           #
#########################################################################################################################
    ### <1> ### Object scoring, ranking and filtering methods                                                           #
#########################################################################################################################
def setupStatFilter(callobj,statlist=[],filterlist=[]):  ### Makes StatFilter dictionary from statlist and filterlist
    '''
    Makes StatFilter dictionary from statlist and filterlist (from cmd_list) !!! Changes case of statfilter keys. !!!
    >> callobj:RJE_Object [None] = calling object for Error Messages etc.
    >> statlist:list of stats that are allowed for filtering. Generally column headers for output.
    >> filterlist:list of StatFilters read in from commandline consisting of StatOperatorValue 
    << statfilter:dictionary of StatFilter {Stat:(Operator,String,Numeric)}
    '''
    try:
        ## Setup dictionary ##
        statfilter = {}
        for filter in filterlist:
            ## Extract details ##
            match = rje.matchExp('^(\S*[A-Za-z0-9])(>|>=|=<|=>|<=|==|=|<|!=|<>)(-*[A-Za-z0-9]\S*)$',filter)
            if not match:
                callobj.log.errorLog('Filter "%s" not recognised.' % filter,printerror=False)
                continue
            (stat,op,cutoff) = match
            if op == '<>':
                op = '!='
            if op == '=':
                op = '=='
            if op in ['=>','=<']:
                op = rje.strReverse(op)
            if op not in ['=>','=<','!=','==','>','<']:
                callobj.log.errorLog('Filter "%s" operator "%s" not known!' % (filter,op),printerror=False)
                continue
            ## Check for numeric value ##
            try:
                numcut = float(cutoff)
            except:
                numcut = None
            ## Check stat ##
            if stat not in statlist:
                for h in statlist:
                    if h.lower() == stat.lower():
                        stat = h
                        break
            if stat not in statlist:
                callobj.log.errorLog('Stat "%s" in filter "%s" not found.' % (stat,filter),printerror=False)
                continue
            ## Update dictionary ##
            statfilter[stat] = (op,cutoff,numcut)
        ### Finish ###
        return statfilter
    except:
        callobj.log.errorLog('Error in rje_scoring.setupStatFilter()',quitchoice=True)
        return statfilter
#########################################################################################################################
def statFilter(callobj,data={},statfilter={},inverse=False,filtermissing=False):    ### Filters data according to statfilter, returning filtered data. 
    '''
    Filters data according to statfilter and returns filtered data. Filtering is *exclusive* based on statfilter.
    >> callobj:RJE_Object [None] = calling object for Error Messages etc.
    >> data:dictionary of data on which to filter {AnyKey:{stat:value}}
    >> statfilter:dictionary of StatFilter {Stat:(Operator,String,Numeric)}
    >> inverse:bool = Whether to reverse filter [False]
    >> filtermissing:bool [False] = whether to treat missing data as a "fail" and filter it [False]
    << data:dictionary of filtered data. This is a *the same* dictionary and must be pre-reassigned if original needed.
    '''
    try:
        ### ~ [1] ~ Filter patterns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for stat in statfilter.keys():      # {stat:(op,cutoff,numcut)}
            (op,strcut,numcut) = statfilter[stat]
            for key in data.keys()[0:]:
                ## ~ [1a] ~ Check for stat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not data[key].has_key(stat):
                    if filtermissing:
                        if not inverse: data.pop(key)
                    else: callobj.log.errorLog('Data for "%s" missing stat "%s"!' % (key,stat),printerror=False)
                    continue
                value = data[key][stat]
                ## ~ [1b] ~ Numeric? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                numeric = None
                if numcut != None:
                    try: numeric = float(data[key][stat])
                    except: numeric = None
                ## ~ [1c] ~ Evaluate and Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if numcut != None and numeric != None: (value,cutoff) = (numeric,numcut)
                else: cutoff = strcut
                #X#print stat, key, numcut, numeric, value, op, cutoff
                try:
                    popdata = False
                    if op == '==' and value == cutoff: popdata = True
                    elif op == '!=' and value != cutoff: popdata = True
                    elif op == '>=' and value >= cutoff: popdata = True
                    elif op == '>' and value > cutoff: popdata = True
                    elif op == '<=' and value <= cutoff: popdata = True
                    elif op == '<' and value < cutoff: popdata = True
                    if inverse: popdata = not popdata
                    if popdata: data.pop(key)
                except:
                    callobj.log.errorLog('Problem filtering %s by %s %s' % (key,stat,statfilter[stat]),printerror=False)
                    break
        ### Finish ###
        return data
    except:
        callobj.log.errorLog('Problem during rje_scoring.statFilter()',quitchoice=True)
        return data
#########################################################################################################################
def setupCustomScores(callobj,statlist=[],scorelist=[],scoredict={}):   ### Sets up Custom Scores using existing statlist
    '''
    Sets up Custom Scores using existing statlist.
    >> callobj:RJE_Object [None] = calling object for Error Messages etc.
    >> statlist:list of stats that are allowed for custom score. Generally column headers for output.
    >> scorelist:list of Custom Score Names in order they were read in (may use prev. scores)   
    >> scoredict:dictionary of Custom Scores: {Name:Formula}
    << (statlist,scorelist,scoredict):(list,list,dictionary) of acceptable Custom Scores ([Stats],[Names],{Name:Formula})
    '''
    try:
        ### Setup Custom Scores ###
        if not scorelist:
            scorelist = rje.sortKeys(scoredict)
        for new in scorelist[0:]:   # self.dict['NewScore'] keys() in order they were read in
            if new in statlist:
                callobj.log.errorLog('Score "%s" exists: custom score cannot be made.' % (new),printerror=False)
                scorelist.remove(new)
                scoredict.pop(new)
                continue
            if not rje.formula(callobj,formula=scoredict[new],varlist=statlist[0:],check=True,calculate=False):
                callobj.log.errorLog('Custom score "%s" cannot be made.' % (new),printerror=False)
                scorelist.remove(new)
                scoredict.pop(new)
                continue
            statlist.append(new)
        return (statlist,scorelist,scoredict)   ### Returns same things given ###
    except:
        callobj.log.errorLog('Problem during rje_scoring.setupCustomScores()',quitchoice=True)
        return scoredict
#########################################################################################################################
def rankObj(callobj,objlist,dkey,dlist=['stat'],case=False,default=None,rev=True,absolute=True,lowest=True,addstat='Rank',cutoff=0,convert=True): ### Ranks objects using numerical data
    '''
    Ranks objects using data in object and rje.getData()
    >> objlist:list of objects to be ranked, ordered and returned
    >> dkey:str = Key for dictionaries
    >> dlist:list [self.stat] = list of dictionaries to try after dict['Data']
    >> case:bool [False] = whether to match case for dkey
    >> default [None] = what value to give entry if no dictionary has key (if None, will not be returned in ranked list)
    >> rev:bool [True] = whether high values should be ranked number 1
    >> absolute:boolean [True] = return 1 to n, rather than 0 to 1
    >> lowest:boolean [True] = returns lowest rank rather mean rank in case of ties
    >> addstat:str ['Rank'] = key for callobj.stat to add rank to (will not add if None)
    >> cutoff:int [0] = only returns top X ranked motifs (if given)  (Can be float if absolute=False)
    >> convert:bool [True] = convert returned data into numeric
    << returns list of ranked objects
    '''
    try:
        ### NewRanks ###
        scores = []
        objlist = objlist[0:]
        for obj in objlist[0:]:
            score = obj.getData(dkey,dlist,case,str=False,default=default,dp=-1)
            if convert:
                try:
                    score = float(score)
                except:
                    score = default
            if score == None:
                objlist.remove(obj)
            else:
                scores.append(score)
        newranks = rje.rankList(scores,rev,absolute,lowest)

        ### Rerank ###
        rankdict = {}
        for i in range(len(objlist)):
            obj = objlist[i]
            r = newranks[i]
            if rankdict.has_key(r):
                rankdict[r].append(obj)
            else:
                rankdict[r] = [obj]
            if addstat:
                obj.stat[addstat] = r
                obj.dict['Data'][addstat] = r
        newlist = []
        for r in rje.sortKeys(rankdict):
            if cutoff <=0 or r <= cutoff:
                newlist += rankdict[r]

        ### Finish ###
        return newlist
    except:
        callobj.log.errorLog('Problem during rje_scoring.rankObj()',quitchoice=True)
        return objlist[0:]
#########################################################################################################################
def statFilterObj(callobj,objlist,statfilter={}):    ### Filters data according to statfilter, returning filtered data. 
    '''
    Filters data according to statfilter and returns filtered data. Filtering is *exclusive* based on statfilter.
    >> callobj:RJE_Object [None] = calling object for Error Messages etc.
    >> objlist:list of objects to be filtered
    >> statfilter:dictionary of StatFilter {Stat:(Operator,String,Numeric)}
    << objlist:list of filtered objects. This is a *the same* list and must be pre-reassigned if original needed.
    '''
    try:
        ### Filter patterns ###
        for stat in statfilter.keys():      # {stat:(op,cutoff,numcut)}
            (op,strcut,numcut) = statfilter[stat]
            for obj in objlist[0:]:
                ## Check for stat ##
                value = obj.getData(stat,str=False)
                if value == None:
                    callobj.log.errorLog('Object data for missing stat "%s"!' % (stat),printerror=False)
                    continue
                ## Numeric? ##
                numeric = None
                if numcut != None:
                    try:
                        numeric = float(value)
                    except:
                        numeric = None
                ## Evaluate and Filter ##
                if numcut != None and numeric != None:
                    (value,cutoff) = (numeric,numcut)
                else:
                    cutoff = strcut
                #X#print stat, key, numcut, numeric, value, op, cutoff
                try:
                    if op == '==' and value == cutoff:
                        objlist.remove(obj)
                    elif op == '!=' and value != cutoff:
                        objlist.remove(obj)
                    elif op == '>=' and value >= cutoff:
                        objlist.remove(obj)
                    elif op == '>' and value > cutoff:
                        objlist.remove(obj)
                    elif op == '<=' and value <= cutoff:
                        objlist.remove(obj)
                    elif op == '<' and value < cutoff:
                        objlist.remove(obj)
                except:
                    callobj.log.errorLog('Problem filtering by %s %s' % (stat,statfilter[stat]),printerror=False)
                    break
        ### Finish ###
        return objlist
    except:
        callobj.log.errorLog('Problem during rje_scoring.statFilter()',quitchoice=True)
        return data
#########################################################################################################################
    ### <2> ### Miscellaneous probability and ranking methods                                                           #
#########################################################################################################################
def adjustedProb(scorelist,reverse=False):  ### Returns the adjust probability value for each score
    '''
    Returns the adjust probability value for each score.
    >> scorelist:list of scores (low is bad)
    >> reverse:bool [False] = reverse so that low is good!
    '''
    ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if len(scorelist) == 1: return [0.5]
    looklist = scorelist[0:]            # Scorelist in order
    looklist.sort()             
    L = len(looklist)                   # Number of scores
    (minp,mult) = (1.0/L,(L-2.0)/L)     # Scaling factors
    adj = {}                            # Adjusted probabilities
    ### ~ [2] ~ Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    print looklist
    for i in range(L):
        s = looklist[i]
        if adj.has_key(s):
            print s, i, adj[s]
            continue
        j = i + 1
        while j < L and looklist[j] == s: j += 1    # j is now the first position after i where s(j) > s(i)
        if j < L:
            x = i / float(L - j)
            p = x / (1 + x)
        else:
            x = 1.0
            p = 1.0
        adj[s] = minp + mult * p
        print s, i, (L-j), x, p, adj[s]
    ### ~ [3] ~ Return List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rlist = []
    for s in scorelist:
        if reverse: rlist.append(adj[s])
        else: rlist.append(1.0 - adj[s])
    return rlist
#########################################################################################################################
### END OF SECTION II                                                                                                   #
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
