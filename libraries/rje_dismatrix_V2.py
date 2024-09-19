#!/usr/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_dismatrix
Description:  Distance Matrix Module 
Version:      2.11.1
Last Edit:    12/05/23
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    DisMatrix Class. Stores distance matrix data and contains methods for extra calculations, such as MST. This module
    is primarily for use within other modules but can be used for slmple distance matrix conversions and UPGMA tree
    construction.    

Commandline:
    loadmatrix=FILE : Loads a matrix from FILE [None]
    symmetric=T/F   : Whether the matrix should be symmetrical (e.g. DisAB = DisBA) [False]
    outmatrix=X     : Type for output matrix - text / mysql / phylip / png
    nsf2nwk=T/F     : Whether to convert extension for Newick Standard Format from nsf to nwk (for MEGA) [False]

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, random, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Working version used in several programs.
    # 2.0 - Updated module inline with newer modules and layout. Incorporation of Minimum Spanning Tree.
    # 2.1 - Attempted in incorporate pylab heatmap generation from distance matrix.
    # 2.2 - Added loading of matrix from GABLAM-style database file.
    # 2.3 - Added UPGMA method.
    # 2.4 - Added UPGMA branch length return.
    # 2.5 - Added PNG output based on rje_slimhtml.
    # 2.6 - Added nsf2nwk=T/F - Whether to convert extension for Newick Standard Format from nsf to nwk (for MEGA)
    # 2.7 - Added reading of BLOSUM matrix and standalone functionality. Removed pylab heatmap from docstring.
    # 2.8 - Slightly modified clustering output.
    # 2.9 - Modified PNG output to use rje_tree code instead. Added Graph output using rje_ppi.
    # 2.10- Minor modifications for SLiMCore UPC.
    # 2.11.0 - Initial Python3 code conversion.
    # 2.11.1 - Fixed some Py3 keys errors.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Possibly add a mergeCase() method to combine mixed case versions of same keys.
    # [ ] : Finish heat map implementation at some point.
    # [ ] : Consider implementing simple NJ method. (To learn how!)
    # [ ] : Add XGMML output
    # [ ] : Check and replace old rje_dismatrix in all modules
    # [ ] : Check/Add NSF output from loaded dismatrix
    # [ ] : Upgrade to new rje_obj class structure.
    # [ ] : General improvement of docstrings etc.
    # [ ] : Add manual details for standalone functionality.
    # [ ] : Add rest of standalone function, including UPGMA.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_DISMATRIX', '2.11.1', 'May 2023', '2007')
    description = 'Distance Matrix Module'
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
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
        print('Problem during initial setup.')
        raise
#########################################################################################################################
mat_ext = {'none':'csv','text':'txt','mysql':'csv','csv':'csv','tdt':'tdt','phylip':'phy'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: DisMatrix Class                                                                                         #
#########################################################################################################################
class DisMatrix(rje.RJE_Object):     
    '''
    Sequence Distance Matrix Class. Author: Rich Edwards (2007). This class can handle *positive* distances only!

    Info:str
    - Name = Name of matrix
    - Type = Identifying type (e.g. 'PWAln ID'). Should be same as key in SeqList.obj
    - Description = Extra information if desired
    - OutMatrix = Type for output matrix - text / mysql / phylip
    
    Opt:boolean
    - nsf2nwk = Whether to convert extension for Newick Standard Format from nsf to nwk (for MEGA) [False]
    - Symmetric = Whether Seq1->Seq2 = Seq2->Seq1

    Stat:numeric

    List:list

    Dict:dictionary
    - Matrix = dictionary of Object pairs and their distance {Obj1:{Obj2:Dis}}

    Obj:RJE_Objects
    '''
    def objNum(self): return len(self.dict['Matrix'])
    def sortObj(self): return rje.sortKeys(self.dict['Matrix'])
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','Description','OutMatrix']
        self.optlist = ['Symmetric','nsf2nwk']
        self.statlist = []
        self.listlist = []
        self.dictlist = ['Matrix']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdRead(cmd,type='info',att='OutMatrix')
                self._cmdRead(cmd,type='file',att='Name',arg='loadmatrix')
                self._cmdReadList(cmd,'opt',['Symmetric','nsf2nwk'])  
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Basic Matrix Methods                                                                                    #
#########################################################################################################################
    def run(self):  ### Will try simple loading and saving of distance matrix
        '''Will try simple loading and saving of distance matrix.'''
        try:### ~ [1] ~ Load distance matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.loadMatrix()
            ### ~ [2] ~ Save reformatted matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: filext = mat_ext[self.info['OutMatrix'].lower()]
            except: filext = 'csv'
            outfile = '%s.%s' % (self.basefile(),filext)
            rje.backup(self,outfile)
            delimit = rje.delimitFromExt(filext,write=True)
            self.saveMatrix(filename=outfile,delimit=delimit,default=0.0)
        except: self.errorLog('Problem with dismatrix conversion run')            
#########################################################################################################################
    def addDis(self,obj1,obj2,dis):     ### Adds a distance to matrix
        '''
        Adds a distance to matrix.
        >> obj1 and obj2:key Objects
        >> dis:Float = 'Distance' measure
        '''
        if obj1 in self.dict['Matrix']: self.dict['Matrix'][obj1][obj2] = dis
        else: self.dict['Matrix'][obj1] = {obj2:dis}
#########################################################################################################################
    def getDis(self,obj1,obj2,default=None):     ### Returns distance from matrix or None if comparison not made.
        '''
        Returns distance from matrix or None if comparison not made.
        >> obj1 and obj2: Objects
        >> default:anything = distance returned if comparison not made.
        << Float = 'Distance' measure
        '''
        try: return self.dict['Matrix'][obj1][obj2]
        except: pass    # Cannot do that!
        if self.opt['Symmetric']:
            try: return self.dict['Matrix'][obj2][obj1]
            except: pass    # Cannot do that!
        return default
#########################################################################################################################
    def objName(self,obj):  ### Returns object name
        '''Returns object name.'''
        try: return obj.shortName()
        except:
            try: return obj.info['Name']
            except: return obj
#########################################################################################################################
    def rename(self,newnames={},missing='keep'):   ### Goes through matrix and renames objects using given dictionary
        '''
        Goes through matrix and renames objects using given dictionary.
        >> newnames:dict = mapping of existing names to new names
        >> missing:str = treatment of names missing as keys (keep/delete)
        '''
        try:### ~ [1] Go through self.dict['Matrix'] and replace everything ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            missing = missing.lower()
            newmatrix = {}
            for obj1 in list(self.dict['Matrix'].keys())[0:]:
                if obj1 not in newnames and missing[:3] in ['del','rem']: continue
                if obj1 in newnames: new1 = newnames[obj1]
                else: new1 = obj1
                newmatrix[new1] = {}
                for obj2 in self.dict['Matrix'][obj1]:
                    if obj2 not in newnames and missing[:3] in ['del','rem']: continue
                    if obj2 in newnames: new2 = newnames[obj2]
                    else: new2 = obj2
                    newmatrix[new1][new2] = self.dict['Matrix'][obj1][obj2]
            self.dict['Matrix'] = newmatrix
        except: self.log.errorLog('Problem during rje_dismatrix_V2.rename()')
#########################################################################################################################
    def remove(self,obj):   ### Removes object for self.dict['Matrix']
        '''Removes object for self.dict['Matrix'].'''
        if obj in self.dict['Matrix']: self.dict['Matrix'].pop(obj)
        for obj1 in self.dict['Matrix']:
            if obj in self.dict['Matrix'][obj1]: self.dict['Matrix'][obj1].pop(obj)
#########################################################################################################################
    ### <3> ### Symmetry Methods                                                                                        #
#########################################################################################################################
    def checkSymmetry(self,force=False):    ### Checks symmetry of matrix and forces if desired
        '''Checks symmetry of matrix and forces if desired.'''
        try:### ~ [1] Assess symmetry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sym = True
            for obj1 in self.dict['Matrix']:
                for obj2 in self.dict['Matrix'][obj1]:
                    if self.getDis(obj1,obj2) != self.getDis(obj1,obj2):
                        sym = False
                        break

            ### ~ [2] Force and/or update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.opt['Symmetric'] = sym
            if force: self.forceSymmetry()
            return self.opt['Symmetric']
        except: self.log.errorLog(rje_zen.Zen().wisdom())
        return False
#########################################################################################################################
    def forceSymmetry(self,method='min',missing=-1):   ### Compresses an Asymmetrical matrix into a symmetrical one
        '''
        Compresses an Asymmetrical matrix into a symmetrical one.
        >> method:str = method for combining distances: min/max/mean
        >> missing:float = Value to be given to missing distances. If < 0 will ignore these values.
        '''
        ### Setup ###
        if self.opt['Symmetric']: return
        self.opt['Symmetric'] = True
        if method not in ['min','max','mean']:
            self.log.errorLog('Force Symmetry method "%s" not recognised. Will use "min".' % method,printerror=False)
            method = 'min'
        ### Force symmetry ###
        newmatrix = DisMatrix(self.log,self.cmd_list)
        for obj1 in self.dict['Matrix']:
            for obj2 in self.dict['Matrix'][obj1]:
                if newmatrix.getDis(obj1,obj2,-1) >= 0: continue    # Already there
                dis1 = self.getDis(obj1,obj2,default=missing)
                dis2 = self.getDis(obj2,obj1,default=missing)
                if dis1 >=0 and dis2 >=0:   ## Combine data
                    if dis1 == dis2: newmatrix.addDis(obj1,obj2,dis1)
                    elif method == 'min': newmatrix.addDis(obj1,obj2,min(dis1,dis2))
                    elif method == 'max': newmatrix.addDis(obj1,obj2,max(dis1,dis2))
                    elif method == 'mean': newmatrix.addDis(obj1,obj2,sum([dis1,dis2])/2.0)
                elif dis1 >= 0 or dis2 >= 0: newmatrix.addDis(obj1,obj2,max(dis1,dis2))
                else: continue
                newmatrix.addDis(obj2,obj1,newmatrix.dict['Matrix'][obj1][obj2])
        ### Update ###
        self.dict['Matrix'] = newmatrix.dict['Matrix']
#########################################################################################################################
    def maxDis(self):   ### Returns maximum distance in matrix
        '''Returns maximum distance in matrix.'''
        maxdis = 0
        for obj1 in self.dict['Matrix']:
            for obj2 in self.dict['Matrix'][obj1]: maxdis = max(maxdis,self.getDis(obj1,obj2,-1))
        return maxdis
#########################################################################################################################
    def minDis(self,missing=1.0):   ### Returns minimum distance in matrix
        '''
        Returns maximum distance in matrix.
        >> missing:float [1.0] = Value to be assigned missing distances (and starting min dis)
        '''
        mindis = missing
        for obj1 in self.dict['Matrix']:
            for obj2 in self.dict['Matrix'][obj1]: mindis = max(mindis,self.getDis(obj1,obj2,missing))
        return mindis
#########################################################################################################################
    def minDisPair(self,missing=1.0):   ### Returns minimum distance in matrix
        '''
        Returns maximum distance in matrix.
        >> missing:float [1.0] = Value to be assigned missing distances (and starting min dis)
        >> returnpair:bool [False] = whether to return pair of objects with mindis rather than distance
        '''
        if len(self.dict['Matrix']) < 2: raise ValueError
        mindis = missing
        pair = self.sortObj()[:2]
        for obj1 in self.dict['Matrix']:
            for obj2 in self.dict['Matrix'][obj1]:
                if obj2 == obj1: continue
                if not pair or self.getDis(obj1,obj2,missing) < mindis:
                    pair = (obj1,obj2)
                    mindis = self.getDis(obj1,obj2,missing)
        return (pair[0],pair[1],mindis)
#########################################################################################################################
    ### <4> ### Minimum Spanning Tree (MST) Methods - adapted from Norman                                               #
#########################################################################################################################
    def MST(self,objkeys=[],normalisation=1.0): ### Calculate MST size for listed objects in matrix
        '''
        Calculate MST size for listed objects in matrix. Adapted from Norman Davey's SlimDiscApp.
        >> objkeys:list = list of object keys for MST calculation
        >> normalisation:float = MST weighting factor (similarities raised to this power)
        '''
        ### Setup ###
        self.forceSymmetry()
        if not objkeys: objkeys = rje.sortKeys(self.dict['Matrix'])
        min_span_tree = [None] * len(objkeys)
        maxdis = max(1.0,self.maxDis())
        high_scores = {}        #objkeys:0.0}
        for obj in objkeys: high_scores[obj] = maxdis
        high_scores[objkeys[0]] = 0.0
        min_queue = objkeys[0:]

        ### Object Loop ###
        while min_queue:
            obj1 = min_queue.pop(0)     # Obj with lowest high_score
            edges = self.findEdges(objkeys,obj1)
            for obj2 in edges:  # Works through other proteins backwards (why backwards?!)
                dis = pow(self.getDis(obj1,obj2,1.0),normalisation)
                if obj2 in min_queue and dis < high_scores[obj2]:
                    min_span_tree[objkeys.index(obj2)] = obj1
                    high_scores[obj2] = dis
                    min_queue = self.sortMinQueue(min_queue,high_scores)
        return 1 + sum(high_scores.values())

        ### Neighbours ###
        neighbours = {}
        for obj1 in objkeys:
            obj2 = min_span_tree[objkeys.index(obj1)]
            if obj2:
                sim = pow(self.getDis(obj1,obj2,1.0),normalisation)
                if obj1 in neighbours: neighbours[obj1].append([obj2,sim])
                else: neighbours[obj1] = [[obj2,sim]]
                if obj2 in neighbours: neighbours[obj2].append([obj1,sim])
                else: neighbours[obj2] = [[obj1,sim]]

        return [1+sum(high_scores.values()),neighbours]
#########################################################################################################################
    def sortMinQueue(self,min_queue,high_scores):   ### Returns min_queue list ordered by high_score
        '''Returns min_queue list ordered by high_score (smallest to highest).'''
        for i in range(len(min_queue)):
            for j in range(i):
                if high_scores[min_queue[i]] <= high_scores[min_queue[j]]:
                    (min_queue[i],min_queue[j]) = (min_queue[j],min_queue[i])
        return min_queue
#########################################################################################################################
    def findEdges(self,objkeys,obj):   ### Not sure why but it returns a reversed list of all other objects
        '''Not sure why but it returns a reversed list of all other objects.'''
        edges = objkeys[0:]
        edges.reverse()
        if obj in edges: edges.remove(obj)
        return edges
#########################################################################################################################
    ### <5> ### Matrix input/output methods                                                                             #
#########################################################################################################################
    def loadFromDataTable(self,filename=None,delimit=None,clear=True,key1='Qry',key2='Hit',distance='Qry_OrderedAlnID',normalise=100.0,inverse=True,checksym=False):  ### Loads distance matrix from e.g. GABLAM output
        '''
        Loads distance matrix from e.g. GABLAM output. (Defaults are for GABLAM)
        >> filename:str = Input file name. Will use self.info['Name'] if None.
        >> delimit:str = Text delimiter. Will ascertain from filename if None.
        >> clear:bool = Whether to clear the current matrix before loading the new one.
        >> key1:str ['Qry'] = column identifying obj1 key for matrix
        >> key2:str ['Hit'] = column identifying obj2 key for matrix
        >> distance:str ['Qry_OrderedAlnID'] = column storing distance measure to read
        >> normalise:float [100.0] = normalise distances to 0-1 scale by dividing by this number (if non-zero)
        >> inverse:bool [True] = whether to subtract value from 1 to get distance (reading a similarity measure)
        >> checksym:bool [False] = Whether to check symmetry of matrix. Will enforce and/or update self.opt['Symmetry']
        '''
        try:### ~ [1] Setup filename and options. Check for file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not filename: filename = self.info['Name']
            if not os.path.exists(filename):
                self.log.errorLog('Matrix data file "%s" not found' % filename,printerror=False)
                return False
            if not delimit: delimit = rje.delimitFromExt(filename=filename)
            if clear: self.dict['Matrix'] = {}

            ### ~ [2] Load distance matrix from filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = rje.dataDict(self,filename,mainkeys=[key1,key2],datakeys=[distance],delimit=delimit)
            (dx,dtot) = (0.0,float(len(data)))
            for pair in list(data.keys())[0:]:
                self.log.printLog('\r#DIS','Reading "%s" matrix from %s: %.1f%%' % (distance,filename,dx/dtot),newline=False,log=False)
                dx += 100
                (obj1,obj2) = rje.split(pair,delimit)
                dis = rje.atof(data.pop(pair)[distance])
                if normalise: dis /= normalise
                if inverse: dis = 1.0 - dis
                self.addDis(obj1,obj2,dis)
            self.log.printLog('\r#DIS','Reading "%s" matrix from %s: %s distances read.' % (distance,filename,rje.integerString(dtot)))
            
            ### ~ [3] Check symmetry of matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if checksym: self.checkSymmetry(force=self.opt['Symmetric'])
            return True
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            return False
#########################################################################################################################
    def loadMatrix(self,filename=None,delimit=None,clear=True,usecase=False,checksym=True,default=0.0):   ### Loads distance matrix 
        '''
        Loads distance matrix from a delimited file. The first column should contain unique identifiers, each of which
        should also have a column heading of its own. If self.opt['Symmetry'] is True then symmetry will be enforced.
        >> filename:str = Input file name. Will use self.info['Name'] if None.
        >> delimit:str = Text delimiter. Will ascertain from filename if None.
        >> clear:bool = Whether to clear the current matrix before loading the new one.
        >> usecase:bool = Enforces matching of case for rows and columns. Else with match any case column to rows.
        >> checksym:bool = Whether to check symmetry of matrix. Will enforce and/or update self.opt['Symmetry']
        >> default:float = Default values to be given to empty cells.
        '''
        try:### ~ [1] Setup filename and options. Check for file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not filename: filename = self.info['Name']
            if not os.path.exists(filename):
                self.log.errorLog('Matrix file "%s" not found' % filename,printerror=False)
                return False
            if not delimit: delimit = rje.delimitFromExt(filename=filename)
            if clear: self.dict['Matrix'] = {}
            ## ~ [1a] Check for BLOSUM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blosum = False
            MAT = open(filename,'r')
            line = MAT.readline()
            while line:
                if line[:1] != '#': break
                if 'BLOSUM' in rje.split(line): blosum = True; break
                line = MAT.readline()
            MAT.close()
            if blosum: return self.loadBLOSUM(filename)

            ### ~ [2] Load distance matrix from filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = rje.dataDict(self,filename,delimit=delimit,getheaders=True)
            ## ~ [2a] Setup matching dictionary for different case in column headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = data.pop('Headers')
            casemap = {}
            if not usecase:
                for h in headers:
                    if h not in data:   # Needs case mapping
                        for obj1 in data:
                            if obj1.lower() == h.lower():
                                casemap[h] = obj1
                                break
            ### ~ [2b] Update matrix from data dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for obj1 in data:
                for h in data[obj1]:
                    if h in casemap: obj2 = casemap[h]
                    else: obj2 = h
                    try: self.addDis(obj1,obj2,rje.atof(data[obj1][h]))
                    except: self.addDis(obj1,obj2,default)

            ### ~ [3] Check symmetry of matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if checksym: self.checkSymmetry(force=self.opt['Symmetric'])
            return True
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            return False
#########################################################################################################################
    def loadBLOSUM(self,filename):   ### Loads distance matrix from BLOSUM file
        '''Loads distance matrix from BLOSUM file.'''
        try:### ~ [1] Read data from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            MAT = open(filename,'r')
            line = MAT.readline()
            while line:
                if line[:1] != '#': aalist = rje.split(line); break
                line = MAT.readline()
            for a1 in aalist:
                data = rje.split(MAT.readline())
                for a2 in aalist:
                    score = rje.atoi(data[aalist.index(a2)])
                    self.addDis(a1,a2,score)
            MAT.close()
            ### ~ [2] Check symmetry of matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.checkSymmetry(force=True)
            return True
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            return False
#########################################################################################################################
    def saveMatrix(self,objkeys=None,filename='matrix.txt',delimit=',',format=None,log=True,default=1.0):   ### Saves matrix 
        '''
        Saves matrix in file. Uses None for missing values.
        >> objkeys:list of Objects that form keys of matrix in output order
        >> filename:str = output file name
        >> delimit:str = separator between columns
        >> format:str = 'text' [Default], 'mysql' = lower case header, 'phylip' = phylip
        >> log:Boolean = whether to print report to log [True]
        >> default:anything = distance to use when no distance present in dictionary
        '''
        try:### ~ [1] Setup formatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if objkeys == None: objkeys = self.sortObj()
            if not format: format = self.info['OutMatrix']
            if format == 'None' and self.opt['MySQL']: format = 'mysql'
            elif format == 'None': format = 'text'
                
            ### ~ [2] Output to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            MATRIX = open(filename, 'w')
            ## ~ [2a] Header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            header = ['SEQ']
            for obj in objkeys: header.append(self.objName(obj))
            if format == 'mysql': MATRIX.write('%s\n' % rje.join(header,delimit).lower())
            elif format == 'phylip':
                MATRIX.write('%d\n' % len(objkeys))
                delimit = ' '
            else: MATRIX.write('%s\n' % rje.join(header,delimit))
            ## ~ [2b] Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for obj1 in objkeys:
                if format == 'phylip':
                    try:
                        if len(self.objName(obj1)) < 10: sname = self.objName(obj1)
                        else: sname = '%d' % (objkeys.index(obj1)+1)
                    except: sname = obj1[:10]
                    while len(sname) < 10: sname = '%s ' % sname
                    matline = [sname]
                else: matline = [self.objName(obj1)]
                for obj2 in objkeys: matline.append(str(self.getDis(obj1,obj2,default)))
                MATRIX.write('%s\n' % rje.join(matline,delimit))
            ## ~ [2c] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            MATRIX.close()
            if log: self.log.printLog('#MAT','%s Distance matrix out to %s (%s format).' % (self.info['Name'],filename,format))
        except: self.log.errorLog(rje_zen.Zen().wisdom())  
#########################################################################################################################
    def saveCytoscape(self,basename=None,type='dis',cutoff=1.0,inverse=False,selflink=False):   ### Outputs matrix in Cytoscape format
        '''
        Outputs matrix in Cytoscape format.
        >> basename:str [None] = Basic name for files. Will use rje.baseFile(self.info['Name']) if missing.
        >> type:str ['dis'] = Type of interaction for *.sif file.
        >> cutoff:float [1.0] = Distances >= this will not be output as links
        >> inverse:bool [False] = If True, distances below cutoff will not be output
        >> selflink:bool [False] = whether to output self distances
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basename: basename = rje.baseFile(self.info['Name'])
            sif = '%s.sif' % basename   # Basic network
            eoa = '%s.eoa' % basename   # Edge attributes
            SIF = open(sif,'w')
            EDGE = open(eoa,'w')
            EDGE.write('Distance\n')

            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            objlist = rje.sortKeys(self.dict['Matrix'])
            for obj1 in objlist:
                for obj2 in objlist:
                    if self.opt['Symmetric'] and objlist.index(obj2) < objlist.index(obj1): continue
                    elif obj1 == obj2 and not selflink: continue
                    dis = self.getDis(obj1,obj2,cutoff)
                    if dis >= cutoff and not inverse: continue
                    elif inverse and dis <= cutoff: continue
                    SIF.write('%s %s %s\n' % (self.objName(obj1),type,self.objName(obj2)))
                    EDGE.write('%s (%s) %s = %s\n' % (self.objName(obj1),type,self.objName(obj2),dis))
            SIF.close()
            EDGE.close()
            self.log.printLog('#SIF','%s matrix output in Cytoscape format to %s and %s' % (basename,sif,eoa))
        except: self.log.errorLog('rje_dismatrix.saveCytoscape(%s) failure' % basename)            
#########################################################################################################################
    def treeName(self,name,nospace=False):     ### Reformats name to be OK in NSF file
        '''Reformats name to be OK in NSF file.'''
        rep = name[0:]
        rep = re.sub(',',' -',rep)
        if nospace: rep = re.sub('\s','_',rep)
        rep = re.sub('[\(\[]','{',rep)
        rep = re.sub('[\)\]]','}',rep)
        rep = re.sub('[:]','-',rep)
        rep = re.sub('[;]','.',rep)
        return rep
#########################################################################################################################
    def savePNG(self,tree,basefile=None,nsftree=None,savensf=True,bycluster=0,singletons=True,byclusters=[]):
        '''
        Saves Matrix as heatmap and tree to file using R code.
        >> tree:rje_tree.Tree object = Needs to be given to method to avoid circularity
        >> basefile:str [None] = Name for output file. Will use basefile by default.
        >> nsftree:str [None] = NSF Tree for output. Will make from distance matrix using UPGMA if None.
        >> savensf:bool [True] = Whether to save trees using tree.saveTrees() - will not save separate PNG if True.
        >> bycluster:int [0] = Number of sequences to split up according to clustering
        >> singletons:bool [True] = Whether to include singletons in main output
        >> byclusters:list [] = Pre-existing list of clusters to use for bycluster analysis []
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basefile: basefile = self.basefile()
            byclusters = byclusters[0:]
            ## ~ [0a] ByCluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if (bycluster > 0 or not singletons) and not nsftree:
                if byclusters:
                    singles = []
                    for clust in byclusters[0:]:
                        if len(clust) < 2: singles += clust; byclusters.remove(clust)
                    if singles: byclusters.insert(0,singles)
                else: byclusters = self.cluster(singletons=False)
                if singletons: self.savePNG(tree,basefile,nsftree,savensf)
                else:
                    nosing = []
                    for cluster in byclusters[1:]: nosing += cluster
                    cnsf = self.upgma(objkeys=nosing)
                    self.savePNG(tree,basefile,cnsf,savensf)
                if bycluster > 1:
                    newcluster = []
                    for cluster in byclusters[0:]:
                        #self.deBug('%s: %d < %d = %s' % (cluster,len(cluster),bycluster,len(cluster) < byclusters))
                        if len(cluster) < bycluster: newcluster += cluster; byclusters.remove(cluster)
                    if newcluster: byclusters.insert(0,newcluster)
                    cx = 0
                    for cluster in byclusters:
                        cbase = '%s.%d-%d' % (basefile,cx,len(cluster)); cx += 1
                        self.printLog('#PNG','Generate tree and PNG for %s' % cbase)
                        cnsf = self.upgma(objkeys=cluster)
                        self.savePNG(tree,cbase,cnsf,savensf)
            ## ~ [0b] NoCluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not nsftree: nsftree = self.upgma()
            #if savensf:
            #    if self.opt['nsf2nwk']: open('%s.nwk' % basefile,'w').write(nsftree)
            #    else: open('%s.nsf' % basefile,'w').write(nsftree)
            ### ~ [1] ~ Build tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tree.buildTree(nsftree,type='nsf',postprocess=False)
            if savensf:
                tree.info['Basefile'] = basefile
                tree.saveTrees()
                return
            ### ~ [2] ~ Old Output files for R code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists('%s.tree.csv' % basefile): os.unlink('%s.tree.csv' % basefile)
            tree.rTree('%s.tree.csv' % basefile,seqname='short')
            reduced = tree._vertOrder(internal=False,namelist=True)
            if os.path.exists('%s.heatmap.tdt' % basefile): os.unlink('%s.heatmap.tdt' % basefile)
            self.saveMatrix(reduced,basefile+'.heatmap.tdt',delimit='\t')
            ### ~ [4] ~ Call R to generate graphics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rcmd = '%s --no-restore --no-save --args "dispng" "%s"' % (self.getStr('RPath'),basefile)
            rslimjim = '%s../libraries/r/rje.r' % self.getStr('Path')
            rcmd += ' < "%s" > "%s.r.tmp.txt"' % (rslimjim,basefile)
            self.printLog('#RPNG',rcmd)
            problems = os.popen(rcmd).read()
            ## ~ [4a] ~ Clear up input files for R script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists('%s.png' % basefile) and not self.getBool('Test'): 
                for ext in ['heatmap.tdt','tree.csv','ppi.tdt','r.tmp.txt']:
                    if os.path.exists('%s.%s' % (basefile,ext)): os.unlink('%s.%s' % (basefile,ext))
            else: self.printLog('#RPNG','Generation of %s.png failed!' % basefile)
        except: self.errorLog('Major problem with %s.savePNG' % self)
#########################################################################################################################
    def saveClusters(self,clusters,basefile=None,upc=False,ctext=''):    ### Saves clusters to a file
        '''
        Saves clusters to a file.
        >> clusters:bool = List of lists containing clusters.append
        >> basefile:str [None] = basefile to use for output
        >> upc:bool [False] = whether to format as UPC file
        >> ctext:str [''] = Extra text to include in header (e.g. clustering parameters)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basefile: basefile = self.basefile()
            if upc: filename = '%s.upc' % basefile
            else: filename = '%s.clusters.txt' % basefile
            if ctext: ctext = '; %s' % ctext
            rje.backup(self,filename,appendable=False)
            missing = []
            outclusters = []
            for cluster in clusters:
                newcluster = []
                for obj in cluster: newcluster.append(rje.split('%s' % obj)[0])
                outclusters.append(newcluster)
                missing += newcluster
            ### ~ [1] ~ Output Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT = open(filename,'w')
            if upc:
                mst = self.MST()
                self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST\n' % (self.info['Name'],self.objNum(),len(clusters),mst))
                OUT.write('#%s# %d Seq; %d UPC; %.3f MST%s\n' % (basefile,self.objNum(),len(clusters),mst,ctext))
                rje.writeDelimit(OUT,['UP','N','MST','Seqs'],'\t')
            else:
                OUT.write('#Clustering of "%s": %d objects; %d clusters%s.\n' % (self.info['Name'],self.objNum(),len(clusters),ctext))
                OUT.write('Cluster\tN\tObjects\n')
            cx = 1
            for cluster in outclusters:
                if upc: rje.writeDelimit(OUT,[cx,len(cluster),rje.join(cluster)],'\t')
                else: OUT.write('%s\t%d\t%s.\n' % (rje.preZero(cx,len(clusters)),len(cluster),rje.join(cluster,'; ')))
                cx += 1
                for obj in cluster:
                    if obj in missing: missing.remove(obj)
            if missing:
                if upc: self.printLog('#ERR','Missing! %d: %s.' % (len(missing),rje.join(missing,'; ')))
                else: OUT.write('Missing!\t%d\t%s.' % (len(missing),rje.join(missing,'; ')))
            OUT.close()
        except: self.errorLog('Major problem with %s.saveClusters' % self)
#########################################################################################################################
    ### <6> ### Clustering methods                                                                                      #
#########################################################################################################################
    def upgma(self,sym='mean',nosim=1.0,log=True,objkeys=[],returnlen=False,checksym=True):  ### Generates a UPGMA tree (NSF) from distance matrix
        '''
        Generates a UPGMA tree (NSF) from distance matrix.
        >> sym:str ['mean'] = Symmetry method to employ
        >> nosim:float [1.0] = Distance to be given for missing distances (no similarity)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rename = {}
            backname = {}
            if not objkeys: objkeys = rje.sortKeys(self.dict['Matrix'])
            for obj in rje.sortKeys(self.dict['Matrix']):
                rename[obj] = self.treeName(self.objName(obj),nospace=True)
                backname[rename[obj]] = obj
                if obj in objkeys: objkeys.insert(objkeys.index(obj),rename[obj]); objkeys.remove(obj)
            self.rename(rename)
            oxx = len(objkeys)
            if checksym: self.checkSymmetry()   # UPGMA needs a symmetrical distance matrix
            self.forceSymmetry(method=sym)      # Force symmetry if not already a feature
            cx = 0                              # No. clades generated so far
            clades = {}                         # Dictionary of clade: [objects]
            depth = {}                          # Depth of clades
            totlen = 0.0
            nosim = max(nosim,self.maxDis())
            ## ~ [0a] ~ Setup UPGMA Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.deBug(objkeys); self.deBug(self.dict['Matrix'])
            upgma = DisMatrix(self.log,self.cmd_list)   # UPGMA dismatix object for calculations
            upgma.opt['Symmetric'] = True
            for obj1 in objkeys:
                for obj2 in self.dict['Matrix'][obj1]:
                    if obj2 in objkeys: upgma.addDis(obj1,obj2,self.getDis(obj1,obj2,nosim))
            cladex = upgma.objNum() - 1          # Max number of clades to generate
            if cladex <= 0:
                self.printLog('#UPGMA','Too few (%d) objects for UPGMA clustering!' % self.objNum())
                if returnlen: self.rename(backname); return 0.0
                if upgma.objNum(): self.rename(backname); return '%s;' % self.objName(list(self.dict['Matrix'].keys())[0])
                self.rename(backname); return ';'
            ### ~ [1] ~ UPGMA clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while upgma.objNum() > 1:                       # More than one cluster: keep clustering
                ## ~ [1a] ~ Determine next clade pair to cluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if log: self.progLog('\r#UPGMA','Clustering %d clades   |' % upgma.objNum())
                cluster = upgma.minDisPair(nosim)           # Returns (obj1,obj2,distance)
                #self.deBug(cluster); self.deBug(self.getDis(cluster[0],cluster[1])); self.deBug(upgma.getDis(cluster[0],cluster[1]));
                if log and self.v() > 1: self.progLog('\r#UPGMA','Clustering %s and %s   |' % (upgma.objName(cluster[0]),upgma.objName(cluster[1])))
                ## ~ [1b] ~ New clade details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                branchlen = [float(cluster[2])/2.0] * 2     # Ancestral branch length for each clade
                for i in range(2):
                    if cluster[i] in depth: branchlen[i] -= depth[cluster[i]]
                newclade = '(%s:%.4f,%s:%.4f)' % (upgma.objName(cluster[0]),branchlen[0],upgma.objName(cluster[1]),branchlen[1])
                totlen = totlen + sum(branchlen)
                #self.deBug(newclade)
                depth[newclade] = float(cluster[2]) / 2.0   # Store depth of new clade
                clades[newclade] = []                       # Store list of original objects in new clade
                for i in range(2):
                    if cluster[i] in clades: clades[newclade] += clades[cluster[i]]
                    else: clades[newclade].append(cluster[i])
                    upgma.remove(cluster[i])                # Remove clustered object from matrix
                ## ~ [1c] ~ Replace old clades with new clade ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newmatrix = {}                              # Dictionary of {upgma obj:new clade distance}
                for obj in list(upgma.dict['Matrix'].keys()):     # Take each object in turn
                    newdis = []                             # Going to take all-by-all mean
                    oldclade = [obj]
                    if obj in clades: oldclade = clades[obj]
                    for obj1 in oldclade:
                        for obj2 in clades[newclade]:
                            newdis.append(self.getDis(obj1,obj2,nosim))     # Use original distances
                    upgma.addDis(newclade,obj,rje.meanse(newdis)[0])        # Add mean to upgma for next cycle
            if log: self.printLog('\n#UPGMA','UPGMA Clustering of %d objects complete!' % oxx)
            ### ~ [2] ~ Return NSF tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug(newclade)
            self.rename(backname)
            if returnlen: return totlen
            return '%s;' % newclade
        except:
            self.log.errorLog('UPGMA Error')
            raise
#########################################################################################################################
    def wpgma(self,sym='mean',nosim=1.0,log=True,objkeys=[],returnlen=False,checksym=True):  ### Generates a WPGMA tree (NSF) from distance matrix
        '''
        Generates a WPGMA tree (NSF) from distance matrix.
        >> sym:str ['mean'] = Symmetry method to employ
        >> nosim:float [1.0] = Distance to be given for missing distances (no similarity)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rename = {}
            backname = {}
            if not objkeys: objkeys = rje.sortKeys(self.dict['Matrix'])
            for obj in rje.sortKeys(self.dict['Matrix']):
                rename[obj] = self.treeName(self.objName(obj),nospace=True)
                backname[rename[obj]] = obj
                if obj in objkeys: objkeys.insert(objkeys.index(obj),rename[obj]); objkeys.remove(obj)
            self.rename(rename)
            oxx = len(objkeys)
            if checksym: self.checkSymmetry()   # WPGMA needs a symmetrical distance matrix
            self.forceSymmetry(method=sym)      # Force symmetry if not already a feature
            cx = 0                              # No. clades generated so far
            clades = {}                         # Dictionary of clade: [objects]
            depth = {}                          # Depth of clades
            totlen = 0.0
            nosim = max(nosim,self.maxDis())
            ## ~ [0a] ~ Setup WPGMA Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.deBug(objkeys); self.deBug(self.dict['Matrix'])
            wpgma = DisMatrix(self.log,self.cmd_list)   # WPGMA dismatix object for calculations
            wpgma.opt['Symmetric'] = True
            for obj1 in objkeys:
                for obj2 in self.dict['Matrix'][obj1]:
                    if obj2 in objkeys: wpgma.addDis(obj1,obj2,self.getDis(obj1,obj2,nosim))
            cladex = wpgma.objNum() - 1          # Max number of clades to generate
            if cladex <= 0:
                self.printLog('#WPGMA','Too few (%d) objects for WPGMA clustering!' % self.objNum())
                if returnlen: self.rename(backname); return 0.0
                if wpgma.objNum(): self.rename(backname); return '%s;' % self.objName(list(self.dict['Matrix'].keys())[0])
                self.rename(backname); return ';'
            ### ~ [1] ~ WPGMA clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while wpgma.objNum() > 1:                       # More than one cluster: keep clustering
                ## ~ [1a] ~ Determine next clade pair to cluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if log: self.progLog('\r#WPGMA','Clustering %d clades   |' % wpgma.objNum())
                cluster = wpgma.minDisPair(nosim)           # Returns (obj1,obj2,distance)
                #self.deBug(cluster); self.deBug(self.getDis(cluster[0],cluster[1])); self.deBug(wpgma.getDis(cluster[0],cluster[1]));
                if log: self.progLog('\r#WPGMA','Clustering %s and %s   |' % (wpgma.objName(cluster[0]),wpgma.objName(cluster[1])))
                ## ~ [1b] ~ New clade details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                branchlen = [float(cluster[2])/2.0] * 2     # Ancestral branch length for each clade
                for i in range(2):
                    if cluster[i] in depth: branchlen[i] -= depth[cluster[i]]
                newclade = '(%s:%.4f,%s:%.4f)' % (wpgma.objName(cluster[0]),branchlen[0],wpgma.objName(cluster[1]),branchlen[1])
                totlen = totlen + sum(branchlen)
                #self.deBug(newclade)
                depth[newclade] = float(cluster[2]) / 2.0   # Store depth of new clade
                clades[newclade] = []                       # Store list of original objects in new clade
                for i in range(2):
                    if cluster[i] in clades: clades[newclade] += clades[cluster[i]]
                    else: clades[newclade].append(cluster[i])
                ## ~ [1c] ~ Replace old clades with new clade ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.deBug(cluster); self.deBug(wpgma.dict['Matrix'])
                for obj in list(wpgma.dict['Matrix'].keys()):     # Take each object in turn
                    if obj in cluster: continue             # This will be removed
                    # WPGMA is a simple mean, no need to go back to originals
                    self.deBug(obj);
                    wpgma.addDis(newclade,obj,(wpgma.getDis(obj,cluster[0])+wpgma.getDis(obj,cluster[1]))/2.0)      # Add mean to wpgma for next cycle
                for i in range(2): wpgma.remove(cluster[i])                # Remove clustered object from matrix
            if log: self.printLog('\n#WPGMA','WPGMA Clustering of %d objects complete!' % oxx)
            ### ~ [2] ~ Return NSF tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug(newclade)
            self.rename(backname)
            if returnlen: return totlen
            return '%s;' % newclade
        except:
            self.log.errorLog('WPGMA Error')
            raise
#########################################################################################################################
    def cluster(self,maxdis=0,singletons=True):  ### Returns a list of clusters objects
        '''
        Returns a list of clusters objects.
        >> maxdis:num [0] = Value representing maximum distance. (Will read from matrix if 0)
        >> singletons:bool [True] = Whether to include singleton clusters (True), combine into single cluster (False), or remove (None).
        '''
        try:### ~ [1] Setup  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Clusters'] = clusters = []
            clusdict = {}
            if not maxdis: maxdis = self.maxDis()
            self.progLog('\r#CLUST','Peparing clustering...')
            for obj in self.sortObj():
                clusdict[obj] = []
                for obj2 in self.sortObj():
                    if self.getDis(obj,obj2,maxdis) < maxdis: clusdict[obj].append(obj2)
                #self.debug('%s: %s' % (obj,clusdict[obj]))
            ### ~ [2] Make clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while clusdict:
                obj = rje.sortKeys(clusdict)[0]
                clusters.append([obj])
                cluslen = 0
                while cluslen < len(clusters[-1]):
                    self.progLog('\r#CLUST','%d clusters (%d keys remaining)' % (len(clusters),len(clusdict)))
                    cluslen = len(clusters[-1])
                    for obj in clusters[-1]:
                        if obj in clusdict:
                            for p in clusdict.pop(obj):
                                if p not in clusters[-1]: clusters[-1].append(p)
                clusters[-1].sort()
            self.printLog('\r#CLUST','%d clusters generated from %d keys (maxdis=%.2f)' % (len(clusters),self.objNum(),maxdis))
            ### ~ [3] Deal with singletons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not singletons:
                singles = []
                for clust in clusters[0:]:
                    if len(clust) == 1: singles += clust; clusters.remove(clust)
                if singletons == False:
                    clusters.insert(0,singles)
                self.printLog('#CLUST','%d singletons' % len(singles))
            ### ~ [4] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Clusters'] = clusters
            return clusters
        except: self.errorLog('%s matrix.cluster() error' % self.info['Name'])
        return []
#########################################################################################################################
    def makeGraph(self,cutoff=1.0,singletons=True):     ### Returns a nested dictionary for distances below cut-off
        '''
        Returns a nested dictionary for distances below cut-off.
        >> cutoff:num [1.0] = Value representing maximum distance for clustering. (Will cluster BELOW this value)
        >> singletons:bool [True] = Whether to include singleton clusters (True) or remove (False).
        '''
        try:### ~ [1] Setup  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            G = {}
            for obj in self.sortObj():
                gobj = rje.split('%s' % obj)[0]
                G[gobj] = {}
                for obj2 in self.sortObj():
                    gobj2 = rje.split('%s' % obj2)[0]
                    dis = self.getDis(obj,obj2,cutoff)
                    if dis < cutoff: G[gobj][gobj2] = dis
                if not G[gobj] and not singletons: G.pop(gobj)
            return G
        except: self.errorLog('%s matrix.makeGraph() error' % self.info['Name'])
        return {}
#########################################################################################################################
### End of SECTION II: DisMatrix Class                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def sortQueue(min_queue,high_scores):   ### Returns min_queue list ordered by high_score
    '''Returns min_queue list ordered by high_score (smallest to highest).'''
    for i in range(len(min_queue)):
        for j in range(i):
            if high_scores[min_queue[i]] <= high_scores[min_queue[j]]:
                (min_queue[i],min_queue[j]) = (min_queue[j],min_queue[i])
    return min_queue
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
        print('Unexpected error during program setup:', sys.exc_info()[0])
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: DisMatrix(mainlog,cmd_list).run()  #heatMap()
        
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
