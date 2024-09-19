#!/usr/bin/python

# Modified SLiMFinder stats module
# Copyright (C) 2009 Norman E. Davey & Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       ned_rankbydistribution
Description:  Modified SLiMFinder stats module
Version:      1.2
Last Edit:    02/02/14
Copyright (C) 2009 Norman E. Davey & Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a stripped down template for methods only. This is for when a class has too many methods and becomes
    untidy. In this case, methods can be moved into a methods module and 'self' replaced with the relevant object.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module.

Uses general modules: re, copy, random, math, sys, time, os, pickle, sets, string, traceback
Uses RJE modules: rje_seq, rje_uniprot, rje, rje_blast, rje_slim
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import re, copy, random, math, sys, time, os, pickle, string, traceback
import rje, rje_slim, rje_slimcore
try:
   set
except NameError:
   from sets import Set as set
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Initial Compilation based on "25 september version". (Probably not!)
    # 1.1 - Minor bug fixes for SLiMFinder 4.5 SigV and SigPrime implementation.
    # 1.2 - Replaced depracated Set module.
    # 1.2.1 = Updated to work with newer module styles.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Fix any bugs with the new implementation.
    '''
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GENERAL METHODS                                                                                         #
#########################################################################################################################
def patternFromCode(slim): return rje_slim.patternFromCode(slim)  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
#########################################################################################################################
def slimPos(slim): return (rje.count(slim,'-') / 2) + 1  ### Returns the number of positions in a slim
#########################################################################################################################
def slimLen(slim): return len(patternFromCode(slim))    ### Returns length of slim
#########################################################################################################################
def slimDif(slim1,slim2): return rje_slimcore.slimDif(slim1,slim2)  ### Returns no. of dif. pos. between slim1 and slim2
#########################################################################################################################
def reformatMotif(slim):    ### Reformats motif from SLiMFinder "SLiM" format
    '''Reformats motif from SLiMFinder "SLiM" format.'''
    ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    slim = rje_slim.patternFromCode(slim)         # Reformatted motif
    positions = rje_slim.slimPosFromCode(slim)    # No. Defined positions
    gaps = rje_slim.slimLenFromCode(slim) - positions    # No. wildcard "gaps"
    return [slim,positions,gaps]
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: rankByDistribution Class                                                                               #
#########################################################################################################################
class rankByDistribution:
    '''Main class for new SLiMFinder stats. Modified from Norman's RankByDistribution.py module.'''
#########################################################################################################################
    ### <1> ### Setup Class                                                                                             #
#########################################################################################################################
    def __init__(self,options):
        self.options = options
        self.slimFinderObj = ''
#########################################################################################################################
    ### <2> ### CalculateRankings Methods                                                                               #
#########################################################################################################################
    def calculateRankings(self,SLiMFinderObj):    ### Main Class Ranking calculation method
        '''
        Main Class Ranking calculation method. This is run by SLiMFinder if SigPrime selected. If SigV only is True then
        calculateRankingsMeanFix() is used instead.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.slimFinderObj = SLiMFinderObj        # SLiMFinder Object
            #self.options["minimumSupport"] = 0        # ??? Should this always be 0? Was once SLiMFinderObj.stat['MinOcc']
            cutOff = {}                                # Dictionary of ?
            ### ~ [2] Create motif space ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SLiMFinderObj.deBug(self.options)
            ## ~ [2a] Setup motif maker ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #SLiMFinderObj.printLog("#NORM","Creating motif space")
            motifSpace = SLiMFinderObj.dict['Extremf.']
            motifMake = motifMaker()
            amb = {}
            if False:    #!# >>>> REDUNDANT CODE? <<<< #!#
                for a in SLiMFinderObj.list['Equiv']:
                    for upc in SLiMFinderObj.dict['AAFreq']:
                        summer = 0
                        for aa in a:
                            try: summer += SLiMFinderObj.dict['AAFreq'][upc][aa]
                            except: pass
                        SLiMFinderObj.dict['AAFreq'][upc][str(len(amb))] = summer
                    amb[str(len(amb))] = a
            ## ~ [2b] Create Motif List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            motifCounter = motifMake.createMotifsList(SLiMFinderObj)
            ### ~ [3] Estimating Percentage Each Support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            createDistribution = Distribution_creator(self.options,SLiMFinderObj)
            SLiMFinderObj.printLog("#PVAL","Calculating expected p-value distribution",log=False)
            binnedMeanOccProb = createDistribution.calculateExpectationSLiMFinder(motifCounter,SLiMFinderObj)
            ## ~ [3a] Calculating p-values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pValDict = createDistribution.calculateSigDistributions(binnedMeanOccProb,SLiMFinderObj)
            del binnedMeanOccProb
            del motifCounter
            for val in pValDict: cutOff[val] = max(pValDict[val])
            topRanking = self.processSlimFinderObject(SLiMFinderObj,cutOff)            
            pVal_sig = createDistribution.mapMotifsToSig(pValDict,topRanking)            
            ### ~ [4] Calculating Estimated Distribution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SLiMFinderObj.printLog("#RANK","Ranking motifs by true Significance",log=False)
            if self.options["fixes"] == "3":    # SigPrime only
                for i in range(self.options["minLength"],self.options["maxLength"] + 1):    #!# Formerly started at 3
                    for score in topRanking[i]:
                        for motif in topRanking[i][score]:
                            self.sigSlimDist(motif[0],pVal_sig[i][score])
            elif self.options["fixes"] == "4":  # SigPrime and SigV corrections
                for i in range(self.options["minLength"],self.options["maxLength"] + 1):    #!# Formerly started at 3
                    SLiMFinderObj.printLog("#NORM","Scoring motifs of length " +str(i))
                    count = 0
                    for score in topRanking[i]:
                        for motif in topRanking[i][score]:
                            count += 1
                            SLiMFinderObj.progLog("\r#NORM","Motif " + str(count) + " of " + str(len(topRanking[i])))
                            self.sigSlimDist(motif[0],pVal_sig[i][score])
                    SLiMFinderObj.printLog("\r#NORM","Motif " + str(count) + " of " + str(len(topRanking[i])))
            else: raise ValueError
        except: SLiMFinderObj.errorLog('Problem during rankByDistribution.calculateRankings')
        return self.slimFinderObj
#########################################################################################################################
    def calculateRankingsMeanFix(self,SLiMFinderObj):   ### Main Class Ranking calculation 
        '''Main Class Ranking calculation method when SigPrime correction is not used (SigV only).'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.slimFinderObj = SLiMFinderObj
            #self.options["minimumSupport"] = 0
            
            cutOff = {}
            for i in range(self.options["minLength"],self.options["maxLength"] + 1):    #!# Formerly started at 3
                cutOff[i] = 1
                
            topRanking= self.processSlimFinderObject(SLiMFinderObj,cutOff)
            createDistribution = Distribution_creator(self.options,SLiMFinderObj)
        
            SLiMFinderObj.printLog("#RANK","Ranking motifs by expectation",log=False)
            for i in range(3,self.options["maxLength"] + 1):
                ptxt = "Scoring motifs of length " +str(i)
                SLiMFinderObj.progLog("\r#RANK",ptxt)
            
                count = 0
                for score in topRanking[i]:
                    for motif in topRanking[i][score]:
                        count += 1
                        SLiMFinderObj.progLog("\r#RANK",'%s: Motif %s of %s' % (ptxt,rje.integerString(count),rje.integerString(len(topRanking[i]))))
                    
                        k = self.slimFinderObj.slimUP(motif[0])
                        
                        successProb = self.slimProbMeanFix(motif[0]).values()
                    
                        #self.sigSlimDist(motif[0],pVal_sig[i][score])
                        prob = createDistribution.meanFixProb(k,successProb)
                        
                        #self.slimFinderObj.dict['Slim'][motif[0]]['Sig'],self.slimFinderObj.dict['Slim'][motif[0]]['E']
                        
                        #print self.slimFinderObj.dict['Extremf.'][i] ,
                        e = self.slimFinderObj.dict['Extremf.'][i]*prob[0]
                        
                        try:
                            Sigv = rje.poisson(1,e,callobj=self.slimFinderObj)
                        except:
                            print(e)
                            Sigv = 1
                            
                        
                        self.sigSlimDist(motif[0],Sigv)
                SLiMFinderObj.printLog("\r#RANK",'%s: Motif %s of %s.' % (ptxt,rje.integerString(count),rje.integerString(len(topRanking[i]))))
            return self.slimFinderObj
        except: SLiMFinderObj.errorLog('Big Problem with calculateRankingsMeanFix()')
#########################################################################################################################
    def slimProbMeanFix(self,slim):
        try:
            #print self.slimFinderObj.list['UP']
            ###~Calculate prob of 1+ occ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            p1 = {}         # Dictionary of {upc:chance of 1+ occ in upc}
            ##~~Setup pattern and variable-lenght multiplier~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            poslist = []    # List of AA positions in SLiM
            wildlist = []   # List of wildcard lengths in SLiM
            wild = False    # Whether next part is a wildcard length
            mult = 1        # Variable-length multiplier
            for part in rje.split(slim,'-'):      # Split SLiM code in components
                ## Update lists ##
                if wild: wildlist.append(part)
                else: poslist.append(part)
                ## Calculate multiplier ##
                if wild:
                    (minx,maxx) = (self.slimFinderObj.getStat('MaxWild'),0)
                    for x in part:
                        minx = min(minx,int(x))
                        maxx = max(maxx,int(x))
                    mult *= (int(maxx) - int(minx) + 1)
                wild = not wild
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.slimFinderObj.list['UP']:
                if self.slimFinderObj.dict['AADimerFreq']: (k,N,p) = self.slimFinderObj.aaDp1(slim,upc)
                else:
                    ## Setup  parameters for binomial ##
                    N = self.slimFinderObj.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
                    p = 1.0                                 # Probability of SLiM at each position
                    k = 1                                   # Number of successful trials (occurrences)
                    if self.slimFinderObj.getBool('SeqOcc') and self.slimFinderObj.slimOccNum(slim,upc) > 1: k = self.slimFinderObj.slimOccNum(slim,upc)
                    ## Calculate p and N from AAFreq and DimFreq ##
                    
                    for pos in poslist:     # AA position
                        posfreq = 0.0
                        for aa in pos: 
                            #print aa,rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)
                            posfreq += rje.getFromDict(self.slimFinderObj.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        
                        p *= posfreq
                        
                    for dim in wildlist:    # DimerFreq
                        dimfreq = 0.0
                        for x in dim:
                            try: dimfreq += self.slimFinderObj.dict['DimFreq'][upc][int(x)]   # Options for wildcard length
                            except: pass
                        N *= (dimfreq / len(dim))       # Mutliply by mean dimer frequency
                    N *= mult       # Each length variant is effectively another position the SLiM could occur
                    if p > 1: p = 1.0   # Cannot in reality have p > 1!
                    ## Calculate binomial ##
                    
                    p1[upc] = rje.binomial(k,N,p,usepoisson=False,callobj=self.slimFinderObj)
            return p1
            
        except:
            self.slimFinderObj.errorLog('Error with slimProb(%s)' % slim)
            self.slimFinderObj.dict['Slim'][slim]['Prob'] = 1.0
#########################################################################################################################
    def sigSlimDist(self,slim,score):     ### Adds Score to dictionary and checks against significance cutoff
        '''Adds Score to dictionary and checks against significance cutoff.'''
        self.slimFinderObj.dict['Slim'][slim]["Normf."] = score
        if score <= self.slimFinderObj.getStat('ProbCut'): self.slimFinderObj.setStat({'NormSig': self.slimFinderObj.getStat('NormSig')+1}); return True
        return False
#########################################################################################################################
    def reformatMotif(self,motif):
        motifSplit = motif.split("-")
        slim = ""
        gaps = 0
        positions = 0        
        for i in range(0,len(motifSplit)):
            if i%2 == 0:
                positions += 1
                if len(motifSplit[i]) == 1: slim += motifSplit[i]                    
                else: slim += "[" + motifSplit[i] + "]"
            else:
                if len(motifSplit[i]) > 1:
                    slim += "{"
                    for aa in motifSplit[i]: slim += aa + "," 
                    slim = slim[:-1] + "}"
                else: slim += "."*int(motifSplit[i])                    
                gaps += int(motifSplit[i])        
        return [slim,positions,gaps]
#########################################################################################################################
    def processSlimFinderObject(self,SLiMFinderObj,cutOff):
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motifSLiMFinder = {}
            motifcounterSLiMFinder = {}
            topRanking = {}
            counter = 0
            sigtype = {'2':'SigV','3':'SigPrime','4':'SigPrimeV'}[self.options["fixes"]]
            for i in range(self.options["minLength"],self.options["maxLength"] + 1): topRanking[i] = {}
            ### ~ [1] Calculate True Distribution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for slim in SLiMFinderObj.dict["Slim"]:
                try:
                    SLiMFinderObj.progLog('\r#PROB','Calculating SLiM Probabilities %d of %d (%3.2f%%)' % (counter, SLiMFinderObj.slimNum(),100*float(counter)/SLiMFinderObj.slimNum()))                    
                    counter += 1
                    [motif,positions,gaps] =  self.reformatMotif(slim)
                    if positions not in topRanking: topRanking[positions] = {}  # Added due to options setup error
                    motif = slim
                    SLiMFinderObj.sigSlim(slim,rankfilter=False)
                    cum_bi_score = SLiMFinderObj.dict['Slim'][slim]['Prob']
                    if cum_bi_score < cutOff[positions]:
                        if len(topRanking[positions]) < 100:
                            if cum_bi_score in topRanking[positions]: topRanking[positions][cum_bi_score].append([motif,0])
                            else: topRanking[positions][cum_bi_score] = [[motif,0]]                            
                        elif max(topRanking[positions].keys()) > cum_bi_score:
                            del topRanking[positions][max(topRanking[positions].keys())]
                            if cum_bi_score in topRanking[positions]: topRanking[positions][cum_bi_score].append([motif,0])
                            else: topRanking[positions][cum_bi_score] = [[motif,0]]
                except: SLiMFinderObj.errorLog('Problem with processSlimFinderObject for %s' % slim)
            SLiMFinderObj.printLog('\r#PROB','Calculation of %s %s SLiM Probabilities complete.' % (SLiMFinderObj.slimNum(),sigtype))
            return topRanking        
        except: SLiMFinderObj.errorLog('Big Problem with processSlimFinderObject'); raise
#########################################################################################################################
    def createBLASTOBjectPointerDict(self,blastObject):     #!!! Cannot find called anywhere !!!#
        pointerDict = {}
        for object in blastObject:
            pointerDict[object.info["Name"]] = object
        return pointerDict
#########################################################################################################################
    def collapseList(self,listTemp):
        stringTemp = ""
        for values in listTemp:
            stringTemp += values
        return stringTemp
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: motifMaker Class                                                                                        #
#########################################################################################################################
class motifMaker:
#########################################################################################################################
    ### <1> ### Setup Class                                                                                             #
#########################################################################################################################
    def __init__(self):
        self.options = {}
#########################################################################################################################
    ### <2> ### Main MotifMaker methods                                                                                 #
#########################################################################################################################
    def createMotifsList(self,slimfinderObj):    ### Creates Motif Space
        '''Creates Motif Space.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.options["maxLength"] = slimfinderObj.getStat('SlimLen')
            aaOcc = ['A','R','N','D','C','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Q']
            motifCounter = {}
            ### ~ [2] Calculate motif space for each range of motif            
            for k in range(2,self.options["maxLength"] + 1):    #!# Changed to 2 for minic < 3 #!#
                c = range(1,k + 1)
                paths = self.diophantine(k,c)
                merlist = []
                motifs = {}
                for path in paths:
                    mul = 1
                    counter = self.convertPath(path)
                    comb = self.factorial(sum(counter))/self.path_product(counter)
                    [count,motifs] = self.buildMotifPath(path,counter,copy.deepcopy(aaOcc),len(counter) - 1,len(counter) - 1,0,[0]*(sum(path)),motifs,comb,aaOcc)
                motifCounter[k] = motifs
            return motifCounter
        except: slimfinderObj.errorLog('Problem during createMotifsList'); raise
#########################################################################################################################
    def buildMotif(self,aa_list,k,k_max,count,parent,motifList):
        if k == k_max: parent = [0]*(k_max + 1)
        for aa in aa_list:
            parent[k] = aa
            if k == 0:
                count += 1
                motifList.append("".join(parent))
            if k > 0: [count,motifList] = self.buildMotif(aa_list,k -1,k_max,count,parent,motifList)
        return [count,motifList]
#########################################################################################################################
    ### <3> ### Path methods                                                                                            #
#########################################################################################################################
    def diophantine(self,n,c):      ### Not sure what this does
        '''RJE: Not sure what this does.'''
        results = self.diophantine_rec(n,c,len(c) - 1,[],[],0,0)
        nodes = {}
        for path in results[1]:
            for node in path:
                nodes[node.keys()[0]] = node.values()[0]
        paths = []
        for path in results[1]:
            cur_bins = []
            start = path[len(path)-1].keys()[0]
            for i in range(1,nodes[start][0]): cur_bins.append(0)
            cur_bins.append(nodes[start][1])
            while nodes[start][3] > 0:
                start = nodes[start][3]
                cur_bins.append(nodes[start][1])
            paths.append(cur_bins)
        return paths
#########################################################################################################################
    def diophantine_rec(self,n,coef,l,path,results,count,parent):
        parent = count
        for r in range(n/coef[l] + 1):
            count += 1
            if l == 0:
                path.append({count:[l+1,n/coef[l],n-(n/coef[l]),parent]})
                break
            else:
                if n - coef[l]*r > 0:
                    path.append({count:[ l+1,r,n - coef[l]*r,parent]})
                    [path,results,count] = self.diophantine_rec(n - coef[l]*r,coef,l-1,path,results,count,parent)                    
                elif n - coef[l]*r == 0:                    
                    path.append({count:[l+1,r,n - coef[l]*r,parent]})
                    break
        if len(path) >0: results.append(path)
        path = []
        return [path,results,count]     # Why return an empty list? (Never used)
#########################################################################################################################
    def convertPath(self,path):
        counter = []
        for i in range(len(path) - 1, -1,-1): counter += [i +1]*path[i] 
        counter.reverse()
        return counter
#########################################################################################################################
    def buildMotifPath(self,path,counter,aa_list,k,k_max,count,parent,motifs,comb,aaOcc):
        if k == k_max: parent = [""]*(k_max + 1)
        if (len(path[1:]) -(path[1:].count(1) + path[1:].count(0)) > 0):
            if k < k_max:
                for i in range(0,k - 1): parent[i] = ""
                if counter[k] != counter[k+1]:
                    aa_list = copy.deepcopy(aaOcc)
                    for aa in parent[1:]:    
                        try: aa_list.remove(aa[0])
                        except: pass
            while len(aa_list) > 0:
                aa = aa_list[len(aa_list)-1]
                aa_list.remove(aa)
                parent[k] = aa*(counter[k])
                if k == 0:
                    count += 1    
                    motif = "".join(parent)
                    motifs[motif] = comb
                if k > 0:
                    temp = copy.deepcopy(aa_list)    
                    [count,motifs] = self.buildMotifPath(path,counter,temp,k -1,k_max,count,parent,motifs,comb,aaOcc)
        elif counter[k] == 1:
            while len(aa_list) > 0:
                aa = aa_list.pop()
                if aa in parent:
                    [count,motif] = self.buildMotifPath(path,counter,copy.deepcopy(aa_list),k -1,k_max,count,parent,motifs,comb,aaOcc)
                else:
                    parent[k] = aa*(counter[k])
                    if k == 0:
                        count += 1
                        motif = "".join(parent)
                        motifs[motif] = comb
                    if k > 0:
                        [count,motif] = self.buildMotifPath(path,counter,copy.deepcopy(aa_list),k -1,k_max,count,parent,motifs,comb,aaOcc)
        else:
            for i in range(0,len(aa_list)):    
                try: aa = aa_list[i]
                except: aa = aa_list[len(aa_list) - 1]                    
                if counter.count(counter[k]) > 1: aa_list.remove(aa)
                parent[k] = aa*(counter[k])                
                if k == 0:
                    count += 1
                    motif = "".join(parent)
                    motifs[motif] = comb                    
                if k > 0:
                    try:
                        temp = copy.deepcopy(aa_list)
                        temp.remove(aa)
                    except: pass                        
                    [count,motifs] = self.buildMotifPath(path,counter,temp,k -1,k_max,count,parent,motifs,comb,aaOcc)
        if len(aa_list) == 0 and k == -1:
            count += 1
            motif = "".join(parent)
            motifs[motif] = comb
        return [count,motifs]
#########################################################################################################################   
    def loopLength(self,n,depth,max_depth,aas,path,count,merList,combs):
        if depth <= max_depth:
            for i in range(n,-1,-1):
                temp = []
                temp += path + [i]
                [count,merList] =self.loopLength(i,depth+1,max_depth,aas,temp,count,merList,combs)
        else:
            test_path = [path.count(p) for p in set(path)]
            test_path.sort()                
            motif = "".join([aas[v] for v in path])
            merList[motif] = combs[str(test_path)]            
            count +=1            
        return [count,merList]
#########################################################################################################################
    ### <4> ### Maths methods                                                                                           #
#########################################################################################################################
    def permutations(self,n,r): #no replacement
        return self.factorial(n)/(self.factorial(n - r))
#########################################################################################################################
    def chose(self,n,k):
        return self.factorial(n)/(self.factorial(n - k)*self.factorial(k))
#########################################################################################################################
    def multichose(self,n,k):
        return self.factorial(n + k - 1)/(self.factorial(k)*self.factorial(n - 1))
#########################################################################################################################
    def factorial(self,m,s=0):
        value = 1
        if m != 0:
            while m != s:
                value = value*m
                m = m - 1
        return value
#########################################################################################################################
    def path_product(self,counter):
        prod =  1
        for val in counter:
            if val > 0: prod *= self.factorial(val)
        return prod
#########################################################################################################################
### END OF SECTION VI                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: Distribution_creator Class                                                                               #
#########################################################################################################################
class Distribution_creator:   
#########################################################################################################################
    ### <1> ### Setup Class                                                                                             #
#########################################################################################################################
    def __init__(self,options,SLiMFinderObj):
        self.options = options        
        self.startTime = time.time()
        self.SLiMFinderObj = SLiMFinderObj
#########################################################################################################################
    def walltime(self): return self.SLiMFinderObj.wallTime()
#########################################################################################################################
    ### <2> ### Main Methods                                                                                            #
#########################################################################################################################            
    def calculateExpectationSLiMFinder(self,merList,SLiMFinderObject):  ### Corrected expectation calculation
        '''Corrected expectation calculation.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            expectation_Dict = {}
            expectation_Fixed_Dict = {}
            distributionCumulative = {}
            expectationBins = {}
            min_expect = {}
            max_expect = {}
            range_expect = {}
            jumps_expect = {}
            posfreqDict = {}
            binnedDict = {}
            ### ~ [2] Calculate for each length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            for lengths in range(self.options["minLength"],self.options["maxLength"] + 1):            
                sortedMerlist = merList[lengths].keys()
                sortedMerlist.sort()
                expectation_Dict[lengths] = {}
                expectation_Fixed_Dict[lengths] = {}
                ptxt = 'Calculating %dmer motifs' % lengths
                counter = 0                
                for motif in sortedMerlist:                    
                    SLiMFinderObject.progLog("\r#PVAL",'%s: %s of %s (%.2f%%)' % (ptxt,rje.integerString(counter),rje.integerString(len(merList[lengths])),(100*float(counter)/len(merList[lengths]))))
                    sys.stdout.flush()
                    counter += 1                    
                    p1 = {}
                    for upc in SLiMFinderObject.list['UP']:
                        N = SLiMFinderObject.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
                        p = 1.0                                 # Probability of SLiM at each position
                        k = 1                                   # Number of successful trials (occurrences)                        
                        for aa in motif: p *= rje.getFromDict(SLiMFinderObject.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                       
                        if p > 1:
                            p = 1.0   
                                
                        p1[upc] = rje.binomial(k,N,p,usepoisson=False,callobj=SLiMFinderObject)
                            
                        expectation_Dict[lengths][motif] = sum(p1.values())/SLiMFinderObject.UPNum()
                        expectation_Fixed_Dict[lengths][motif] = p1
            
                        self.walltime()
                        
                sortMe = expectation_Dict[lengths].values()
                SLiMFinderObject.printLog("\r#PVAL",'Calculation of corrected expectation for %s %dmer motifs complete.' % (rje.integerString(len(merList[lengths])),lengths))
                sortMe.sort()
                
                val = min(expectation_Dict[lengths].values())
            
                try:
                    exponent = math.floor(math.log10(math.fabs(val)))
                except:
                    exponent = 32
                    
                temp = {}
                #print len(expectation_Dict[lengths])
                for motif in expectation_Dict[lengths]:
                    rounded = round(expectation_Dict[lengths][motif],abs(int(exponent))+1)
                    if rounded in temp:
                        temp[rounded][0] += merList[lengths][motif]
                        temp[rounded][1] += [motif]
                        temp[rounded][2] += [expectation_Dict[lengths][motif]]*merList[lengths][motif]
                        #temp[rounded][3] = expectation_Fixed_Dict[lengths][motif]
                    else:
                        temp[rounded] = [merList[lengths][motif],[motif],[expectation_Dict[lengths][motif]]*merList[lengths][motif],expectation_Fixed_Dict[lengths][motif].values()]
                
                    self.walltime()
                    
                binnedDict[lengths] = {}
            
                
                for entry in temp:
                    meanScore = sum(temp[entry][2])/temp[entry][0]
                    binnedDict[lengths][meanScore] = temp[entry]
                    
        
            return  binnedDict
        except: SLiMFinderObject.errorLog('Problem during calculateExpectationSLiMFinder()'); raise
#########################################################################################################################            
    def multList(self,listM,counts,length,SLiMFinderObject):
        prod = 1.0
        for x in listM:
            prod *= (1 - x)**(counts[x][0][2]*pow(SLiMFinderObject.getStat('MaxWild')-SLiMFinderObject.getStat('MinWild')+1,(length-1)))
        return prod
#########################################################################################################################
    def findIndex(self,p,sorter,index_start,index_end):
        if not sorter or p > sorter[index_end]: return [-1,-1]
        if p < sorter[index_start]: return [0,0]
        if index_end - index_start <= 1: return [index_start,index_end]
        else: 
            if sorter[index_start + (index_end- index_start )/2] > p:
                [index_1,index_2] = self.findIndex(p,sorter,index_start,index_start + ( index_end- index_start )/2)
            elif sorter[index_start + (index_end- index_start )/2] < p:
                [index_1,index_2] = self.findIndex(p,sorter,index_start + (index_end- index_start)/2,index_end)
            else:
                return [index_start + ( index_end- index_start )/2,index_start + (index_end- index_start)/2]
        return [index_1,index_2]
#########################################################################################################################
    def calculateSigDistributions(self,meanMotifProbData,SLiMFinderObject):
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            p_valueDict = {}    # Dictionary of {slimlength:
            shuffleMatrix = {}  # Dictionary of {slimlength:
            pVal_Sid_dict = {}  # Dictionary of {slimlength:
            ### ~ [1] Generate Shuffle Matrices ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for length in range(self.options["minLength"],self.options["maxLength"] + 1):
                ## ~ [1a] Setup sorted Mean Motif Probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #SLiMFinderObject.deBug('%s: meanMotifProbData[length]' % length)
                sortedMeanMotifProb = meanMotifProbData[length].keys()
                sortedMeanMotifProb.sort()
                shuffleMatrix[length] = {}
                p_valueDict[length] = {}
                pVal_Sid_dict[length] = {}
                
                count = 0
                for prob in sortedMeanMotifProb:
                    #print prob,meanMotifProbData[length][prob]
                    
                    count += 1
                    SLiMFinderObject.progLog("\r#SHUFF","Initialising Shuffle Matrix %.3f%%" %(float(count*100)/len(sortedMeanMotifProb)))
                    sys.stdout.flush()
                    
                    temp = [[],[],0]
                    #for k in range(0,SLiMFinderObject.UPNum() + 1):
                    try:
                        for k in range(SLiMFinderObject.UPNum(),0,-1):
                            #print "\n",k,prob,
                            if self.options["fixes"] == "3":
                                bi_vals = self.calculate_p(k,prob,SLiMFinderObject)
                            elif self.options["fixes"] == "4":
                                bi_vals = self.meanFixProb(k,copy.deepcopy(meanMotifProbData[length][prob][3]))
                            else: raise ValueError
                                
                            if bi_vals[0] < 0.1**length and bi_vals[0] > 0.1**16:
                                temp[0] = [bi_vals[0]] + temp[0]
                                temp[1] = [bi_vals[1]] + temp[1]
                                temp[2] = meanMotifProbData[length][prob][0]
                            
                        if temp[0][0] in shuffleMatrix[length]:
                            shuffleMatrix[length][temp[0][0]].append(temp)
                        else:
                            shuffleMatrix[length][temp[0][0]] = [temp]
                    except KeyboardInterrupt: raise
                    except:
                        #if SLiMFinderObject.opt['DeBug']: SLiMFinderObject.errorLog('ShuffleProblem')
                        pass
                    #count += meanMotifProbData[length][prob][0]
                    self.walltime()
                
                sorter = shuffleMatrix[length].keys()
                sorter.sort()
                
                summer = 0
                try: cum_p = sorter.pop(-1)        #!# Pop from empty list for Sigv #!#
                except: SLiMFinderObject.errorLog('Empty shuffleMatrix for length %s!' % length); continue
                temp = shuffleMatrix[length][cum_p]
                last = [temp[0][2],cum_p]
                    
                #mul = (1 - sorter.pop(0))**(shuffleMatrix[length][sorter[0]][0][2]*pow(SLiMFinderObject.stat['MaxWild']-SLiMFinderObject.stat['MinWild']+1,(length-1)))
                mul = self.multList(sorter,shuffleMatrix[length],length,SLiMFinderObject)
                #print "%1.16e"%mul
                count = 0
                while len(sorter) > 0:
                    
                    SLiMFinderObject.progLog("\r#ITER ","%.5f iterating towards 1: Length %s" % (float(mul),str(length)))
                    sys.stdout.flush()
                    
                    cum_p = sorter.pop(-1)
                    #print len(sorter),
                    #mul_old = mul
                    #mul = self.multList(sorter,shuffleMatrix[length],length,SLiMFinderObject)
                    #print mul, mul_old,cum_p,shuffleMatrix[length][cum_p][0][2]
                    #print mal,len(sorter),
                    #print last
                    try:
                        d = ((1 - cum_p)**(shuffleMatrix[length][cum_p][0][2]*pow(SLiMFinderObject.getStat('MaxWild')-SLiMFinderObject.getStat('MinWild')+1,(length-1))))
                        m = (1 - last[-1])**(last[0]*pow(SLiMFinderObject.getStat('MaxWild')-SLiMFinderObject.getStat('MinWild')+1,(length-1)))
                        #print  "%1.5e"%m,"%1.5e"%d,mul,self.multList(sorter,shuffleMatrix[length],length,SLiMFinderObject)
                        if count%500 == 0:
                            mul = self.multList(sorter,shuffleMatrix[length],length,SLiMFinderObject)
                        else:
                            mul = m*(mul/d)
                    except:
                        raise
                    
                    temp = shuffleMatrix[length][cum_p]
                    new_p = temp[0][0][0]
                    last = [temp[0][2],new_p]
                    
                    temp[0][0].pop(0)
                    temp[0][1].pop(0)
                    
                    count += 1#temp[0][2]
                    #mul *= (1 - cum_p)**(temp[0][2]*pow(SLiMFinderObject.stat['MaxWild']-SLiMFinderObject.stat['MinWild']+1,(length-1)))
                    #print mul,count
                    
                    try:
                        p_valueDict[length][cum_p] = 1 - mul
                    except:
                        pass
                    
                    
                    del shuffleMatrix[length][cum_p]
                    
                    
                    if len(temp[0][0]) > 0:
                        try:
                            shuffleMatrix[length][new_p ] = temp
                            [index,index2] = self.findIndex(new_p,sorter,0,len(sorter)-1)
                            
                            if index == -1:
                                sorter.append(new_p)
                            else:
                                sorter.insert(index ,new_p)
                        except:
                            print( temp)
                            SLiMFinderObject.errorLog('Problem during calculateSigDistributions()')
    
                    self.walltime()
                    
                            
                #print count        
                SLiMFinderObject.printLog("#\r#SHUFF","ShuffleMatrix created for " + str(length) + "mers" + " "*30,log=False)
            return p_valueDict
                    
        except: SLiMFinderObject.errorLog('Error during calculateSigDistributions()'); raise
#########################################################################################################################
    def mapMotifsToSig(self,p_valueDict,topRanking):
            pVal_Sid_dict = {}
            
            for length in topRanking:
                try:
                    pVal_Sid_dict[length] = {}
                    sorter = p_valueDict[length].keys()
                    sorter.sort()
                    
                    sortMotifs = topRanking[length].keys()
                    sortMotifs.sort()
                    
                    
                    for p in sortMotifs:
                        try:
                            index = self.findIndex(p,sorter,0,len(sorter)-1)
                            pVal_Sid_dict[length][p] = p_valueDict[length][sorter[index[0]]]
                        except:
                            pVal_Sid_dict[length][p] = 1
                except:
                    raise
                    pass
                    
            return pVal_Sid_dict        
    
    def chose(self,n,k):
        return self.factorial(n)/(self.factorial(n - k)*self.factorial(k))
#########################################################################################################################
    def factorial(self,m,s=0):
        value = 1
        if m != 0:
            while m != s:
                value = value*m
                m = m - 1
        return value        
#########################################################################################################################
    def buildMotifPath(self,aa_list,k,k_max,count,parent,motifs,aaOcc):
    
        if k == k_max:
            parent = [""]*(k_max + 1)
            
        while len(aa_list) > 0:
            aa = aa_list[len(aa_list)-1]
                
            aa_list.remove(aa)
            parent[k] = aa
            
            if k == 0:
                count += 1    
                motif = parent
                motifs[count] = copy.deepcopy(motif)
        
            if k > 0:
                temp = copy.deepcopy(aa_list)
            
                [count,motifs] = self.buildMotifPath(temp,k -1,k_max,count,parent,motifs,aaOcc)
    
        return [count,motifs]
#########################################################################################################################
    def probabilityPerm(self,list1,list2):  ### Product of successes of list1 and failures of list2 (I think!)
        '''Product of successes of list1 and failures of list2 (I think!).'''
        for x in list1: list2.remove(x)            
        prob = 1
        for i in list1: prob *= i
        for k in list2: prob *= 1 - k            
        return prob
#########################################################################################################################
    def meanFixProb(self,k,successProb):
        #print "-"*10
        #print k,len(successProb) ,successProb 
        
        [count,motifs] = self.buildMotifPath(copy.deepcopy(successProb),k - 1,k - 1,0,{},{},successProb)
        
        bi_p = 0
    
        for motif in motifs:
            bi_p += self.probabilityPerm(motifs[motif],copy.deepcopy(successProb))
            
        cum_p = bi_p
        
        for i in range(k+1,len(successProb) +1):
            #successProb = self.slimProbMeanFix(slim).values()
            #print k,len(successProb) ,successProb 
            [count,motifs] = self.buildMotifPath(copy.deepcopy(successProb),i - 1,i - 1,0,{},{},successProb)
    
            p = 0
            
            for motif in motifs:
                p += self.probabilityPerm(motifs[motif],copy.deepcopy(successProb))
            
            cum_p += p
            
        return [cum_p,bi_p]
#########################################################################################################################
    def calculate_p(self,k,expected,SLiMFinderObject):
        try:            
            n = SLiMFinderObject.UPNum()
            p = rje.binomial(k,n,expected,usepoisson=False,exact=True,callobj=SLiMFinderObject)
            cum_p = rje.binomial(k,n,expected,usepoisson=False,callobj=SLiMFinderObject)
            return [cum_p,p]
        except: SLiMFinderObject.errorLog('Error during calculate_p()'); raise
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION VI: 'MAIN' PROGRAM                                                                                          #
#########################################################################################################################
if __name__ == "__main__":  ### Print message to screen if called from commandline.
    try: print( __doc__)
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION VI                                                                                                   #
#########################################################################################################################

