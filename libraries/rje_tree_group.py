#!/usr/bin/python

# rje_tree_group - Contains all the Grouping Methods for rje_tree.py
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_tree_group
Description:  Contains all the Grouping Methods for rje_tree.py
Version:      1.2.1
Last Edit:    28/01/15
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is a stripped down template for methods only. This is for when a class has too many methods and becomes
    untidy. In this case, methods can be moved into a methods module and 'self' replaced with the relevant object. For
    this module, 'self' becomes '_tree'.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module: rje_tree.py

Uses general modules: copy, re, os, string, sys
Uses RJE modules: rje, rje_seq
Other modules needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, re, string, sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_seq
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 1.0 - Initial fully-functional module for use with rje_tree.py
# 1.1 - Modified some of the automated group naming
# 1.2 - Added
# 1.2.1 - Tweaked QryVar interactivity.
#########################################################################################################################
### Major Functionality to Add
# [Y] Change name of group when renaming genes in group
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: METHODS                                                                                                 #
#########################################################################################################################
    ### <A> ### Master Tree Subfamily/Groupings Methods
#########################################################################################################################
def treeGroup(_tree,callmenu=False):     ### Master Tree Grouping loop.
    '''
    Master Tree Grouping loop. This method will:
    1. Attempt to automatically group the sequences according to the current grouping option.
    2. Assess the grouping returned and enter manual mode if:
        (a) _tree.stat['Interactive'] > 1
    -OR-(b) _tree.stat['Interactive'] == 0 and grouping is bad or callmenu is True
    '''
    try:
        ### <1> ### Attempt Auto Grouping
        _stage = '<1> AutoGroup'
        _autoGroups(_tree,method=_tree.info['Grouping'])
        ### <2> ### Exit if appropriate
        _stage = '<2> Exit'
        if _tree.stat['Interactive'] < 0:
            return
        elif _checkGroups(_tree) and _tree.stat['Interactive'] < 1 and callmenu == False:
            return     
        ### <3> ### Give options
        _stage = '<3> Options'
        while _groupChoice(_tree):
            continue
    except:
        _tree.log.errorLog('Problem with rje_tree_group.treeGroup(%s)' % _stage,True)
#########################################################################################################################
def _autoGroups(_tree,method='man'):      ### Automatically groups sequences according to options
    '''
    Automatically goes into grouping routines.
    >> method:str = grouping method
        - man = manual grouping (unless i<0).
        - dup = duplication (all species unless groupspec specified).
        - qry = duplication with species of Query sequence (or Sequence 1) of treeseq
        - one = all sequences in one group
        - None = no group (case sensitive)
        - FILE = load groups from file
    '''
    try:
        ### Setup ###
        if not _tree.node:
            _clearGroups(_tree)
            return
        ### Automated Grouping ###
        if method == 'man':
            _clearGroups(_tree)
        elif method == 'dup':
            _dupGroup(_tree)
        elif method == 'qry':
            if _tree.obj['SeqList'].obj['QuerySeq']:
                _tree.info['GroupSpecies'] = _tree.obj['SeqList'].obj['QuerySeq'].info['SpecCode']
            else:
                _tree.info['GroupSpecies'] = _tree.obj['SeqList'].seq[0].info['SpecCode']
            _dupGroup(_tree)
        elif method == 'one':
            _tree.stat['MinFamNum'] = 1
            _tree.subfam = [_tree.node[-1]]
            _resetGroups(_tree)
            _tree.node[-1].info['CladeName'] = _tree.info['Name']
        elif method == 'None':
            _tree.stat['MinFamNum'] = 0
            _tree.opt['Orphans'] = True
            _clearGroups(_tree)
        else:                
            _loadGroups(_tree,filename=method)
        ### Check Group Names ###
        _checkGroupNames(_tree)
    except:
        _tree.log.errorLog('Major problem with rje_tree_group._autoGroups(%s).' % method)
#########################################################################################################################
def _checkGroups(_tree): ### Checks that Group selection does not break 'rules'
    '''
    Checks that Group selection does not break 'rules':
    - enough groups
    - enough sequences in each group
    - sufficient bootstrap support for each group
    - no variants in group if allowvar=False
    - no orphans if orphans=False
    << True if OK, False if rule(s) broken.
    '''
    try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        query = _tree.obj['SeqList'].obj['QuerySeq']
        ### ~ [1] ~ Group Number ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if len(_tree.subfam) < _tree.stat['MinFamNum']: return False
        ### ~ [2] ~ Group Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for node in _tree.subfam:
            ## ~ [2a] ~ Family Size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            clade = _tree._nodeClade(node)
            if len(clade) < _tree.stat['MinFamSize']: return False
            ## ~ [2b] ~ Variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not _tree.opt['AllowVar']:
                for node1 in clade:
                    for node2 in clade:
                        if node1 != node2 and node1.obj['Sequence'].sameSpec(node2.obj['Sequence']): # Same species (see rje_sequence)
                            if not query or not _tree.opt['QryVar'] or not node1.obj['Sequence'].sameSpec(query): return False
            ## ~ [2c] ~ BootCut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if node != _tree.node[-1]:
                ancbranch = node.ancBranch()
                if int(_tree.stat['Bootstraps'] * _tree.stat['BootCut']) > ancbranch.stat['Bootstrap']: return False
        ### ~ [3] Orphans ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not _tree.opt['Orphans'] and _orphanCount(_tree) > 0: return False
        return True
    except: _tree.errorLog('Fatal error in rje_tree_group._checkGroups()'); raise
#########################################################################################################################
    ### <B> ### Gubbins
#########################################################################################################################
def _checkGroupNames(_tree):    ### Automatically names groups after given gene unless already renamed
    '''Automatically names groups after given gene unless already renamed.'''
    try:
        for fam in _tree.subfam:
            ### Check existing name ###
            clade = _tree._nodeClade(fam)
            if clade[0].obj['Sequence'].shortName() not in fam.info['CladeName']:  ## Renamed already
                #X#print fam.info['CladeName'], 'renamed already!'
                continue
            ### Assemble Gene List ###
            genelist = []
            for node in clade:
                seq = node.obj['Sequence']
                if seq.info['Gene'] not in genelist:
                    genelist.append(seq.info['Gene'])
            ### Get candidate gene name ###
            geneocc = {}
            maxocc = 0
            bestgene = None
            for gene in genelist:
                geneocc[gene] = 0
                for node in clade:
                    seq = node.obj['Sequence']
                    if gene in [seq.info['Gene'],seq.info['Gene'].lower()]:
                        geneocc[gene] += 1
                if geneocc[gene] > maxocc:
                    maxocc = geneocc[gene]
                    bestgene = gene
            ### Assess best gene name ###
            if bestgene.lower() in ['ipi','ens','nvl','unknown','ref','hyp','nem','p','tr']:
                bestgene = None
                maxocc = 0
            #X#print 'Best gene = %s (%.1f%%) %s' % (bestgene,(100.0*maxocc/len(clade)),genelist)
            if maxocc > (0.9 * len(clade)):  ### 90%+
                fam.info['CladeName'] = bestgene
            elif maxocc > (0.5 * len(clade)) and (_tree.stat['Interactive'] < 0 or rje.yesNo('Use "%s" for new name of group "%s"?' % (bestgene,fam.info['CladeName']))):  ### 90%+
                fam.info['CladeName'] = bestgene
            elif _tree.stat['Interactive'] > 0:
                fam.info['CladeName'] = rje.choice('Name for group %s' % fam.info['CladeName'],default=fam.info['CladeName'])
            if fam.info['CladeName'].upper() in genelist:
                fam.info['CladeName'] = fam.info['CladeName'].upper()
    except:
        _tree.log.errorLog('Non-Fatal error in rje_tree_group._checkGroupNames()')
#########################################################################################################################
def _orphanCount(_tree):    ### Returns number of orphan sequences
    '''Returns number of orphan sequences.'''
    try:
        famnodes = []
        for fam in _tree.subfam:
            famnodes += _tree._nodeClade(fam)
        ox = 0
        for node in _tree.node[0:_tree.seqNum()]:
            if node not in famnodes:
                ox += 1
        return ox           
    except:
        _tree.log.errorLog('Problem with rje_tree_group._orphanCount()')
        return 0
#########################################################################################################################
def _purgeOrphans(_tree):   ### Removes orphan nodes
    '''Removes orphan nodes.'''
    try:
        famnodes = []
        for fam in _tree.subfam:
            famnodes += _tree._nodeClade(fam)
        ox = 0
        termini = _tree.node[0:_tree.seqNum()]
        for node in termini:
            if node not in famnodes:
                ox += 1
                _tree._prune(branch=node.ancBranch(),remtext='Orphan Sequence')
        _tree.log.printLog('#GRP','%d Orphan Sequences Removed.' % ox)
    except:
        _tree.log.errorLog('Problem with rje_tree_group._purgeOrphans()')
#########################################################################################################################
def _resetGroups(_tree):    ### Resets node compression to match current group selection
    '''Resets node compression to match current group selection.'''
    try:
        for node in _tree.node:
            if node in _tree.subfam:
                node.opt['Compress'] = True
            else:
                node.opt['Compress'] = False
    except:
        _tree.log.errorLog('Fatal Error with rje_tree_group._resetGroups().')
        raise
#########################################################################################################################
def _addGroup(_tree,node):  ### Adds group based at node, removing existing descendant groups
    '''
    Adds group based at node, removing existing descendant groups.
    >> node:Node Object
    '''
    try:
        ### <a> ## Add group
        _tree.subfam.append(node)
        ### <b> ### Check no descendant nodes are groups
        fams = _tree.subfam[0:]
        for othernode in fams:
            if _tree._isAnc(node,othernode):
                othernode.opt['Compress'] = False
                _tree.subfam.remove(othernode)
    except:
        _tree.log.errorLog('Major problem with rje_tree_group._addGroup()')
        raise
#########################################################################################################################
def _grpVarDel(_tree,delseq=[],kept=None,fam=None):     ### Deletes all variants in list
    '''
    Deletes all variants in list.
    >> delseq:list of Sequences to delete
    >> kept:Sequence Object of kept variant
    >> fam:Node Object that defines subfam
    '''
    try:
        text = 'Probable sequence variant vs %s' % kept.shortName()
        famnode = _tree._nodeClade(fam)
        for seq in delseq:
            for node in famnode:
                if node.obj['Sequence'] == seq:
                    famnode.remove(node)
                    _tree._prune(node.ancBranch(),remtext=text)
                    break
        if fam not in _tree.node:
            _tree.verbose(1,3,'Redefining subfamily node.',1)
            for f in range(len(_tree.subfam)):
                if _tree.subfam[f] == fam:
                   _tree.subfam[f] = _tree._bestCladeNode(famnode) 
    except:
        _tree.log.errorLog('Major problem during rje_tree_group._grpVarDel', True)
#########################################################################################################################
def _grpSeqSort(_tree,seqs=[],compseq=None):  ### Reorders seqs according to %ID, Gaps and Extra   
    '''
    Reorders seqs according to %ID, Gaps and Extra.
    >> seqs:list of sequences to reorder
    >> compseq:Sequence Object to compare seqlist to
    << reordered seqlist
    '''
    try:
        ### <0> ### Setup
        _stage = '<0> Setup'
        seqlist = _tree.obj['SeqList']
        seqsort = {}
        dislist = ['ID','Gaps','Extra']
        dblist = ['ipi', 'sprot', 'trembl', 'ens_known', 'ens_novel', 'ens_scan', 'None']
        for seq in seqs:
            if seq.info['DBase'] not in dblist:
                dblist.append(seq.info['DBase'])
        ### <1> ### Make dislist
        _stage = '<1> Make dislist'
        for seq in seqs:
            seqsort[seq] = ''
            for key in dislist: 
                if key == 'ID':
                    dis = rje.preZero(seq.seqLen() - int(seqlist.getDis(seq,compseq,'MSA %s' % key) * seq.aaLen() / 100),seq.seqLen())
                    seqsort[seq] = '%s%s%s' % (seqsort[seq], dis, key)
                else:
                    dis = seqlist.getDis(seq,compseq,'MSA %s' % key) / 100
                    seqsort[seq] = '%s%.5f%s' % (seqsort[seq], dis, key)
            seqsort[seq] = '%s%s' % (seqsort[seq], rje.preZero(dblist.index(seq.info['DBase']),len(dblist)))
        ### <2> ### Sort values of seqsort and then reform list from seqlist
        _stage = '<2> Sort and reform'
        sortlist = seqsort.values()
        sortlist.sort()
        seqorder = []
        for val in sortlist:
            for seq in seqs:
                if seqsort[seq] == val:
                    seqorder.append(seq)
                    seqs.remove(seq)
                    break
        if len(seqorder) != len(sortlist) or len(seqs) != 0:
            raise ValueError
        return seqorder
    except:
        _tree.log.errorLog('Major problem during rje_tree_group._grpSeqSort(%s):' % _stage)
        raise
#########################################################################################################################
def _reorderSeqToGroups(_tree): ### Reorders sequences according to groups #!# Add query? And _grpSeqSort()?
    '''
    Reorders sequences according to groups.
    '''
    try:
        _stage = '<0> Setup'
        seqlist = _tree.obj['SeqList']
        query = seqlist.obj['QuerySeq']
        neworder = []

        _stage = '<1> Query'
        qfam = _queryGroup(_tree)
        if qfam:
            _tree.subfam.remove(qfam)
            _tree.subfam = [qfam] + _tree.subfam       

        _stage = '<2> Families'
        for fam in _tree.subfam:
            famseq = []
            clade = _tree._nodeClade(fam)
            for node in clade:
                famseq.append(node.obj['Sequence'])
            if query:
                neworder += _grpSeqSort(_tree,seqs=famseq,compseq=query)
            else:
                neworder += famseq

        _stage = '<3> SeqList and Orphans'
        for seq in seqlist.seq:
            if seq not in neworder:
                neworder.append(seq)
        seqlist.seq = neworder
    except:
        _tree.log.errorLog('Major problem during rje_tree_group._reorderSeqToGroups(%s)' % _stage)
#########################################################################################################################
def _queryGroup(_tree): ### Returns Query subfam node or None if none
    '''
    Returns Query subfam node or None if none.
    '''
    try:
        query = _tree.obj['SeqList'].obj['QuerySeq']
        if query == None:
            return None
        for fam in _tree.subfam:
            clade = _tree._nodeClade(fam)
            for node in clade:
                if node.obj['Sequence'] == query:
                    return fam
        return None
    except:
        _tree.log.errorLog('Major problem during rje_tree_group._queryGroup()')
        return None
#########################################################################################################################
    ### <C> ### Interactive Grouping
#########################################################################################################################
def _groupChoice(_tree):    ### Gives manual options for grouping.
    '''
    Gives manual choices for grouping.
    << False if unchanged, True if changed
    '''
    try:
        if 'QryVar' not in _tree.opt: _tree.opt['QryVar'] = _tree.opt['AllowVar']
        ### <0> ### Summary
        _stage = '<0> Summary'
        _checkGroupNames(_tree)
        print '\n#*# Grouping Options #*#\n'
        print '[Bootstrap cut-off: %d (%f)]' % (int(_tree.stat['Bootstraps'] * _tree.stat['BootCut']), _tree.stat['BootCut'])
        print '[Min Family Size: %d]' % _tree.stat['MinFamSize']
        print '[Min Family Number: %d]' % _tree.stat['MinFamNum']
        print '[Group by Query Species: %s (%s)]' % (_tree.opt['QueryGroup'],_tree.info['GroupSpecies'])
        print '[Allow Sequence Variants: %s]' % _tree.opt['AllowVar']
        print '[Allow Query species Variants: %s]' % _tree.opt['QryVar']
        print '[Allow Orphan Sequences: %s]' % _tree.opt['Orphans']
        print '\n#*# Grouping Summary #*#\n',
        _sumGroups(_tree)
#        if _tree.stat['MinFamNum'] > 0:
#            print 'Grouping Needed (Min. %d Families).' % _tree.stat['MinFamNum']
#        else:
#            print 'No Grouping Needed.'

        ### <1> ### Options
        _stage = '<1> Options'
        print ' <K>eep current grouping. ',
        keepok = _checkGroups(_tree)
        if keepok:
            print '(Current Grouping OK.)'
        else:
            print '*** WARNING: Current Grouping breaks 1+ rules. ***'
        print ' --- '
        print " <O>ptions (Change Grouping 'rules')."
        print ' --- '
        print ' <L>oad Grouping.'
        print ' <D>uplication Grouping. ',
        if _tree.opt['QueryGroup']:
            print ' (%s Duplications)' % _tree.info['GroupSpecies']
        else:
            print '(All Duplications)'
        print ' <M>anual Grouping.'
        print ' <N>o Grouping (MinFamNum=0).'
        print ' <A>ll sequences in one group.'
        if len(_tree.subfam) > 0:
            print ' <E>dit Grouping (Manual)'
            print ' --- '
            print ' <R>eview Groups'
            print ' <G>roup in SeqList (Reorder)'
            print ' <P>urge Orphan Sequences (not in a group)'
            print ' <S>ave Groups'
        print ' --- \n <F>ull auto mode'
        print ' --- \n <Q>uit'
        # ! # Line up group options to save vertical space?
        
        ### <2> ### Make Choice
        _stage = '<2> Choice'
        choice = rje.choice('\nChoice for Grouping?',default='K').upper()
        if choice == 'K' and rje.yesNo('Keep Groups?'):   # Keep current grouping
            if _tree.opt['Orphans'] == False:
                if _orphanCount(_tree) > 0:
                    if rje.yesNo('%d Orphan sequences will be deleted. OK?' % _orphanCount(_tree)):
                        _purgeOrphans(_tree)
                    else:
                        return _groupChoice(_tree)
                return False
            if keepok:
                if _tree.stat['Interactive'] > 0 and rje.yesNo('Save groups?'):
                    grpfile = _tree.obj['SeqList'].info['Name']
                    if grpfile[-4] == '.':
                        grpfile = grpfile[:-4]
                    grpfile = rje.choice('Name of Groupfile?','%s.grp' % grpfile)
                    grpnames = rje.yesNo('Write Group Names?',default='N')
                    _saveGroups(_tree,filename=grpfile,groupnames=grpnames)
                return False
            elif rje.yesNo('Grouping rules broken. Are you sure?'):
                return False
            return _groupChoice(_tree)
        elif choice == 'O': # Change Rules
            _groupRules(_tree)
            return _groupChoice(_tree)
        elif choice == 'L': # Load Grouping
            grpfile = _tree.obj['SeqList'].info['Name']
            if grpfile[-4] == '.':
                grpfile = grpfile[:-4]
            grpfile = _tree._editChoice('Name of Groupfile','%s.grp' % grpfile)
            _loadGroups(_tree,filename=grpfile)
            return True
        elif choice == 'D': # Duplication Grouping
            _dupGroup(_tree)
            _tree.info['Grouping'] = 'dup'
            return True
        elif choice == 'M': # Manual Grouping
            _clearGroups(_tree)
            _tree.editTree(groupmode=True)
            _tree.info['Grouping'] = 'man'
            return True
        elif choice == 'N': # No Grouping
            _clearGroups(_tree)
            _tree.stat['MinFamNum'] = 0
            _tree.info['Grouping'] = 'None'
            return True
        elif choice == 'A': # One Group
            _tree.stat['MinFamNum'] = 1
            _tree.subfam = [_tree.node[-1]]
            _resetGroups(_tree)
            _tree.node[-1].info['CladeName'] = _tree.info['Name']
            _tree.info['Grouping'] = 'one'
            return True
        elif choice == 'E': # Manual Grouping
            _resetGroups(_tree)
            _tree.editTree(groupmode=True)
            return True
        elif choice == 'R': # Review Groups
            _reviewGroups(_tree)
            return _groupChoice(_tree)
        elif choice == 'G': # Reorder SeqList according to Grouping
            _reorderSeqToGroups(_tree)
            return True
        elif choice == 'P': # Purge Orphans
            _purgeOrphans(_tree)
            return True
        elif choice == 'S': # Save Grouping
            grpfile = _tree.obj['SeqList'].info['Name']
            if grpfile[-4] == '.':
                grpfile = grpfile[:-4]
            grpfile = _tree._editChoice('Name of Groupfile','%s.grp' % grpfile)
            grpnames = rje.yesNo('Write Group Names?')
            _saveGroups(_tree,filename=grpfile,groupnames=grpnames)
            return True
        elif choice.find('F') == 0 and rje.yesNo('Switch on Full-auto mode? (Cannot be undone)'): # i=-1
            _tree.setStat({'Interactive':-1})
            return False
        elif choice.find('Q') == 0 and rje.yesNo('Quit Grouping?'): # Quit
            return False
        else:   # Failed to choose
            return _groupChoice(_tree)
    except SystemExit:
        raise 
    except KeyboardInterrupt:
        raise
    except:
        _tree.log.errorLog('Major Problem with rje_tree_group._groupChoice(%s)' % _stage)
        return False
#########################################################################################################################
def _sumGroups(_tree):      ### Prints summary of Groups
    '''Prints summary of Groups.'''
    try:
        _tree.verbose(0,3,'Currently %d groups. (%d Orphans)' % (_tree.groupNum(),_tree._orphanCount()),1)
        gx = 0
        for fam in _tree.subfam:
            gx += 1
            clade = _tree._nodeClade(fam)
            _tree.verbose(0,3,'Group %d: %d seqs %s' % (gx,len(clade),_tree.nodeList(clade)),1)
        _tree.verbose(0,1,'',1)
    except KeyboardInterrupt:
        raise
    except:
        _tree.log.errorLog('Problem with rje_tree_group._sumGroups().')
#########################################################################################################################
def _groupRules(_tree):     ### Options to change grouping options
    '''Options to change grouping options.'''
    try:
        print '\nEdit Grouping Rules.\nEnter new values or leave Blank to retain.\n'
        for stat in ['BootCut','MinFamSize','MinFamNum']:
            _tree.stat[stat] = _tree._editChoice(stat,_tree.stat[stat],numeric=True)
        for opt in ['QueryGroup','Orphans','AllowVar','QryVar']:
            _tree.opt[opt] = _tree._editChoice(opt,_tree.opt[opt],boolean=True)
        if _tree.opt['QueryGroup']:
            if _tree.obj['SeqList'].obj['QuerySeq']:
                _tree.info['GroupSpecies'] = _tree.obj['SeqList'].obj['QuerySeq'].info['SpecCode']
            else:
                _tree.info['GroupSpecies'] = _tree.obj['SeqList'].seq[0].info['SpecCode']
        for info in ['GroupSpecies']:
            _tree.info[info] = _tree._editChoice(info,_tree.info[info])
    except:
        _tree.log.errorLog('Major problem with _groupRules().')
#########################################################################################################################
def _saveGroups(_tree,filename='rje_tree.grp',groupnames=True):   # Saves sequence names in Groups
    '''
    Saves sequence names in Groups.
    >> filename:str = group filename
    >> groupnames:boolean = whether to save groupnames
    '''
    try:
        GROUPS = open(filename, 'w')
        gx = 0
        sx = 0
        for fam in _tree.subfam:
            gx += 1
            GROUPS.write('Group %d:' % gx)
            if groupnames:
                GROUPS.write(' %s' % fam.info['CladeName'])
            GROUPS.write('\n')
            clade = _tree._nodeClade(fam)
            for node in clade:
                sx += 1
                GROUPS.write('%s\n' % node.obj['Sequence'].shortName())
            GROUPS.write('\n')
        GROUPS.close()
        _tree.log.printLog('#GRP','%d Groups (%d Sequences) output in %s.' % (gx,sx,filename))
    except:
        _tree.log.errorLog('Major problem with _saveGroups(%s).' % filename)
#########################################################################################################################
    ### <D> ### Grouping Methods
#########################################################################################################################
def _clearGroups(_tree):     ### Clears current group selection
    '''Clears current group selection.'''
    try:
        for node in _tree.node:
            node.opt['Compress'] = False
        _tree.subfam = []
    except:
        _tree.log.errorLog('Fatal Error with rje_tree_group._clearGroups().')
        raise
#########################################################################################################################
def _loadGroups(_tree,filename='rje_tree.grp'):   # Saves sequence names in Groups
    '''
    Saves sequence names in Groups.
    >> filename:str = group filename
    '''
    try:
        ### <a> ### Read in groups
        ## <i> ## File
        gtxt = 'Grouping sequences using %s' % filename
        _tree.log.printLog('#GRP','%s (Loading) ... 0.0%%' % gtxt,newline=False,log=False)
        try:
            GROUPS = open(filename, 'r')
            glines = GROUPS.readlines()
            GROUPS.close()
        except(IOError):
            _tree.log.errorLog("File %s not found" % rootfile)
            raise
        ## <ii> ## Group data
        groups = [] # List of name lists
        gnames = []
        lx = 0
        while lx < len(glines):
            _tree.log.printLog('\r#GRP','%s (Loading) ... %.1f%% (%d groups)' % (gtxt,100.0*lx/len(glines),len(groups)),newline=False,log=False)
            line = glines[lx]
            lx += 1
            rje.chomp(line)
            if re.search('^Group \d+:',line):  # New group
                newgroup = []
                if re.search('^Group \d+:\s+(\S.+)$',line):
                    gname = rje.matchExp('^Group \d+:\s+(\S.+)$',line)[0]
                else:
                    gname = 'None'
                while lx < len(glines):
                    line = glines[lx]
                    lx += 1
                    rje.chomp(line)
                    if re.search('^(\S+)',line):
                        newgroup.append(rje.matchExp('^(\S+)',line)[0])
                    else:
                        break
                groups.append(newgroup[0:])
                gnames.append(gname)
        _tree.log.printLog('\r#GRP','%s (Loading) ... 100.0%% (%d groups)' % (gtxt,len(groups)),log=False)

        ### <b> ### Map to nodes
        _clearGroups(_tree)
        gx = 0
        sx = 0
        for groupseq in groups:    # Each group is a list of sequence names
            _tree.log.printLog('\r#GRP','%s (Mapping) ... %.1f%% (%d seqs)' % (gtxt,100.0*groups.index(groupseq)/len(groups),sx),newline=False,log=False)
            ## <i> ## Convert list of names into list of nodes
            group = []
            for name in groupseq:
                for node in _tree.node[:_tree.stat['SeqNum']]:
                    if node.info['Name'].find(name) == 0:
                        group.append(node)
            ## <ii> ## Identify common ancestral node   
            gx += 1
            fam = _tree._bestCladeNode(group)
            if fam == None:
                _tree.log.errorLog('No sequences read in for Group %d (%s)' % (gx,gname[gx-1]))
                continue
            else:
                _tree.subfam.append(fam)
                fam.info['CladeName'] = gnames[gx-1]
                clade = _tree._nodeClade(fam)
                if fam.info['CladeName'] == 'None':    #!# Make this a method?
                    fam.info['CladeName'] = '%d Seqs %s' % (len(clade),_tree.nodeList(clade))
                sx += len(clade)
        _tree.log.printLog('\r#GRP','%s (Mapping) ... 100.0%%' % gtxt,log=False)
        _resetGroups(_tree)
        _tree.log.printLog('#GRP','%d Groups (%d Sequences) loaded from %s.' % (_tree.groupNum(),sx,filename))
        _sumGroups(_tree)
        _tree.textTree()
        _tree.info['Grouping'] = filename
    except:
        _tree.log.errorLog('Major problem with rje_tree_group._loadGroups(%s).' % filename)
#########################################################################################################################
def _dupGroup(_tree):     ### Duplication grouping
    '''Duplication grouping.'''
    try:
        ### <a> ## Find Duplications - maybe species specific
        if _tree.opt['QueryGroup'] and _tree.info['GroupSpecies'] != 'None':
            _tree.findDuplications(species=_tree.info['GroupSpecies'])
        else:
            _tree.findDuplications()
        ### <b> ### Identify duplications
        _clearGroups(_tree)
        ## <i> ## Work down tree and assign
        for node in _tree.node[:-1]:
            #print node.info['Name']
            if node.ancNode().opt['Duplication'] and node.ancBranch().stat['Bootstrap'] >= int(_tree.stat['Bootstraps'] * _tree.stat['BootCut']):
                clades = _tree._descClades(node)
                clades = clades[0] + clades[1]
                if len(clades) < _tree.stat['MinFamSize']: continue
                node.opt['Compress'] = True
                _tree.subfam.append(node)
                if node.info['CladeName'] == 'None':
                    node.info['CladeName'] = '%d Seqs %s' % (len(clades),_tree.nodeList(clades))
        ## <ii> ## Check no descendant nodes are duplications
        fams = _tree.subfam[0:]
        for node in fams:
            for othernode in fams:
                if node in _tree.subfam and _tree._isAnc(node,othernode):
                    node.opt['Compress'] = False
                    _tree.subfam.remove(node)
                    #break
        ## <iib> ## Lineage-specific duplications ##
        fams = _tree.subfam[0:]
        for node in _tree.node:
            if node in fams: continue
            elif node.opt['SpecDup']:
                clades = _tree._descClades(node)
                clades = clades[0] + clades[1]
                node.opt['Compress'] = True
                _tree.subfam.append(node)
                if node.info['CladeName'] == 'None': node.info['CladeName'] = '%d Seqs %s' % (len(clades),_tree.nodeList(clades))
                for othernode in _tree.subfam:
                    if node in _tree.subfam and _tree._isAnc(othernode,node):
                        node.opt['Compress'] = False
                        _tree.subfam.remove(node)
                        break
                for othernode in _tree.subfam[0:]:
                    if node.opt['Compress'] and _tree._isAnc(node,othernode):
                        othernode.opt['Compress'] = False
                        _tree.subfam.remove(othernode)
                        
        _tree.textTree()
        _tree.info['Grouping'] = 'dup'
        ## <iii> ## If no duplications, put all sequences in one group!
        if _tree.groupNum() == 0:
            node = _tree.node[-1]
            _tree.subfam.append(node)
            node.opt['Compress'] = True
            if node.info['CladeName'] == 'None':
                node.info['CladeName'] = 'All sequences.'
        return True
    except:
        _tree.log.errorLog('Major problem with rje_tree_group._dupGroup().')
        return False
#########################################################################################################################
    ### <E> ### Review Grouping 
#########################################################################################################################
def _reviewGroups(_tree,interactive=1):    ### Summarise, scan for variants (same species), edit group
    '''
    Summarise, scan for variants (same species), edit group.
    '''
    try:
        ### <0> ### Setup
        _stage = '<0> Setup'
        if _tree.groupNum() == 0:
            return
        masterspec = None
        seqlist = _tree.obj['SeqList']
        tmpdismatrix = rje_seq.DisMatrix(log=_tree.log)
        for key in ['MSA ID','MSA Gaps', 'MSA Extra']:
            if seqlist.obj[key] == None:
                seqlist.obj[key] = copy.deepcopy(tmpdismatrix)
                seqlist.obj[key].matrix = {}
        query = seqlist.obj['QuerySeq']
        qtxt = ''
        if query:
            qtxt = 'vs Query (%s)' % query.shortName()
        if _tree.info['GroupSpecies'] != 'None':
            masterspec = _tree.info['GroupSpecies']
        elif query:
            masterspec = query.info['SpecCode']
        _tree.log.printLog('#GRP','Reviewing %d Groups. (Master SpecCode = %s)' % (_tree.groupNum(),masterspec))

        ### <1> ### Loop through groups until satisfied
        _stage ='<1> Group Loop'
        g = 0
        while g < _tree.groupNum() and _tree.groupNum() > 0:  ### Continue review
            ## <a> ## Draw tree of group
            _stage ='<1a> Group Loop tree'
            fam = _tree.subfam[g]
            if fam == None or fam not in _tree.node:    # Been deleted
                _tree.log.errorLog('Subfam %d has been removed but clings to life in _tree.subfam. Destroying it and moving on.',printerror=False)
                _tree.subfam.pop(g)
                continue
            fam.opt['Compress'] = False
            if _tree.stat['Verbose'] > 0 or _tree.stat['Interactive'] >= 0:
                _tree.textTree(fromnode=fam,seqname='long',maxnamelen=0)
            fam.opt['Compress'] = True

            ## <b> ## Summary info = length + %ID, gaps & extra vs first (human) sequence
            _stage ='<1b> Summary Info'
            grpseq = _tree._nodeSeqs(_tree._nodeClade(fam))     # List of Sequence Objects in group
            seqspec = {}    # Dictionary of species counts
            for node in _tree._nodeClade(fam):
                seqspec[node.obj['Sequence'].info['SpecCode']] = 0
            for node in _tree._nodeClade(fam):
                seqspec[node.obj['Sequence'].info['SpecCode']] += 1

            ## <c> ## Select grpmaster = node of correct species (if any)
            ## .. within species (or all) pick node with smallest mean distance to all other nodes (make a MSA ID matrix if missing)
            _stage ='<1c> Group Master'
            masters = []
            mastervar = 0
            if masterspec:
                for seq in grpseq:
                    if seq.info['SpecCode'] == masterspec:
                        masters.append(seq)
                        mastervar += 1
            if len(masters) == 0:
                masters = grpseq
            if query:
                masters = _grpSeqSort(_tree,masters,query)
                grpmaster = masters[0]
            else:
                grpmaster = seqlist.bestMeanID(queries=masters,seqlist=grpseq,key='MSA ID')
            if mastervar == 1:
                masters = [grpmaster]
            _tree.verbose(0,4,'GroupMaster: Group %d, %s' % (g,fam.info['CladeName']),1)
            _tree.verbose(0,4,_groupDisSum(_tree,grpmaster,query,qtxt),2)

            ## <d> ## Calculate: %id, %gaps, %extra for each sequence versus grpmaster
            ## .. for each species, order nodes according to %gaps then %id then % extra            
            _stage ='<1d> Summary vs Group Master'
            grpseq = _tree._nodeSeqs(_tree._nodeClade(fam))     # List of Sequence Objects in group
            grpseq = _grpSeqSort(_tree,grpseq,grpmaster)
            for seq in grpseq[1:]:
                _tree.verbose(0,4,_groupDisSum(_tree,seq,grpmaster),1)
                
            ## <e> ## Look for variants in Master sequence! (vs Query)
            ## .. work through qryspec in order and accept as best sequence or remove 
            _stage ='<1e> Group Master Species Variants'
            for seq in masters[1:]: # Only deal with same species as GroupMaster
                if seq.info['SpecCode'] != grpmaster.info['SpecCode']:
                    masters.remove(seq)
            vardel = False; nextbreak = False
            while len(masters) > 1:
                _tree.verbose(0,3,'\n%d %s variants. ("Best" = %s) Choice:' % (len(masters),grpmaster.info['SpecCode'],masters[0].shortName()),1)
                ctext = ''
                if _tree.stat['Verbose'] > 0 or _tree.stat['Interactive'] >= 0:
                    for m in range(len(masters)):
                        print '\n<%d> %s' % (m, _groupDisSum(_tree,masters[m],query,qtxt)),
                        ctext += '<%d> %s; ' % (m, masters[m].shortName())
                    ctext += '<K> keep variants; <N>ext group:'
                if _tree.stat['Interactive'] < 1 and _tree.opt['QryVar'] and query and grpmaster.info['SpecCode'] == query.info['SpecCode']: choice = len(masters)
                elif _tree.stat['Interactive'] >= 0:
                    if (_tree.opt['AllowVar'] or _tree.opt['QryVar']) and query and grpmaster.info['SpecCode'] == query.info['SpecCode']: choice = rje.choice('\n %s?: ' % ctext,default='k').lower()
                    #choice = rje.getInt('\n %s?: ' % ctext,blank0=True)
                    else: choice = rje.choice('\n %s?: ' % ctext,default='0').lower()
                    if choice == 'k': choice = len(masters)
                    elif choice == 'n': choice = -1
                    else:
                        try: choice = string.atoi(choice)
                        except: continue
                elif _tree.opt['AllowVar']: choice = len(masters)
                elif _tree.opt['QryVar'] and query and grpmaster.info['SpecCode'] == query.info['SpecCode']: choice = len(masters)
                else: choice = 0

                if choice == 0:
                    masters.pop(0)
                    _grpVarDel(_tree,masters,grpmaster,fam)
                    vardel = True
                    break
                elif choice > 0 and choice < len(masters):
                    grpmaster = masters.pop(choice)
                    _grpVarDel(_tree,masters,grpmaster,fam)
                    vardel = True
                    break
                elif choice == len(masters):
                    if _tree.opt['AllowVar'] or (_tree.opt['QryVar'] and query and grpmaster.info['SpecCode'] == query.info['SpecCode']):
                        if _tree.stat['Interactive'] < 1 or rje.yesNo('Keep all variants?'): break
                    elif rje.yesNo('Keep all variants? (WARNING: AllowVar = False!)'):
                        _tree.opt['AllowVar'] = True
                        break
                elif choice < 0:
                    if g < (_tree.groupNum() - 1): g += 1; nextbreak = True
                    break
            if nextbreak: continue

            ## <f> ## Scan for Variants (i.e. same species)
            _stage ='<1f> Scan for Variants'
            seqspec.pop(grpmaster.info['SpecCode'])
            for spec in seqspec.keys():
                vseq = []
                if seqspec[spec] > 1:   # Variants
                    fam = _tree.subfam[g]
                    grpseq = _tree._nodeSeqs(_tree._nodeClade(fam))     # List of Sequence Objects in group
                    for seq in grpseq:
                        if seq.info['SpecCode'] == spec:
                            vseq.append(seq)
                while len(vseq) > 1:
                    vseq = _grpSeqSort(_tree,vseq,grpmaster)
                    _tree.verbose(0,3,'\n%d %s variants. ("Best" = %s) Choice:' % (len(vseq),spec,vseq[0].shortName()),1)
                    ctext = ''
                    if _tree.stat['Verbose'] > 0 or _tree.stat['Interactive'] >= 0:
                        for v in range(len(vseq)):
                            print '<%d> %s' % (v, _groupDisSum(_tree,vseq[v],grpmaster))
                            ctext += '<%d> %s; ' % (v, vseq[v].shortName())
                        ctext += '<K> keep variants; <N>ext group'

                    if _tree.stat['Interactive'] < 1 and _tree.opt['QryVar'] and query and vseq[0].info['SpecCode'] == query.info['SpecCode']: choice = len(vseq)
                    elif _tree.stat['Interactive'] >= 0:
                        #choice = rje.getInt('%s?: ' % ctext,blank0=True)
                        if (_tree.opt['AllowVar'] or _tree.opt['QryVar']) and query and vseq[0].info['SpecCode'] == query.info['SpecCode']: choice = rje.choice('\n %s?: ' % ctext,default='k').lower()
                        choice = rje.choice('\n %s?: ' % ctext,default='0').lower()
                        if choice == 'k': choice = len(vseq)
                        elif choice == 'n': choice = -1
                        else:
                            try: choice = string.atoi(choice)
                            except: continue
                    elif _tree.opt['AllowVar']: choice = len(vseq)
                    elif _tree.opt['QryVar'] and query and vseq[0].info['SpecCode'] == query.info['SpecCode']: choice = len(vseq)
                    else: choice = 0

                    if choice < 0:
                        if g < (_tree.groupNum() - 1): g += 1; nextbreak = True
                        break
                    elif choice < len(vseq):
                        var = vseq.pop(choice)
                        _grpVarDel(_tree,vseq,var,fam)
                        vardel = True
                        break
                    elif choice == len(vseq):
                        if _tree.opt['AllowVar'] or (_tree.opt['QryVar'] and query and vseq[0].info['SpecCode'] == query.info['SpecCode']):
                            if _tree.stat['Interactive'] < 1 or rje.yesNo('Keep all variants?'): break
                        elif rje.yesNo('Keep all variants? (WARNING: AllowVar = False!)'):
                            _tree.opt['AllowVar'] = True
                            break
            if vardel or nextbreak: continue
        
            ## <g> ## Option to Edit group/species - treeEdit/gene names/full sequence edit, Option to Keep or Delete group = proceed or prune tree
            _stage ='<1g> Group Edit Options'
            _tree.verbose(0,4,'\nSubfamily Options:',1)
            while g < _tree.groupNum():
                if g > 0:
                    _tree.verbose(0,4,'<P>revious, ',0)
                _tree.verbose(0,4,'<T>ree Edit, <G>ene Edit, <U>ngroup, <D>elete Group, Group <M>enu, ',0)
                if g < (_tree.groupNum() - 1):
                    _tree.verbose(0,4,'<Q>uit Review, <N>ext',0)
                    default = 'N'
                else:
                    _tree.verbose(0,4,'<F>inish Review',0)
                    default = 'F'
                if _tree.stat['Interactive'] < 0:
                    choice = default
                else:
                    choice = rje.choice(' ?:',default=default).upper()
                if choice not in ['P','T','G','U','D','Q','N','','F','M']:
                    continue
                ### Tree Edit
                if choice.find('T') == 0:
                    fam.opt['Compress'] = False
                    _tree.editTree(fromnode=fam)    #!# Check this works!
                    if fam in _tree.node:
                        fam.opt['Compress'] = True
                    else:
                        _tree.subfam.remove(fam)
                ### New Gene
                elif choice.find('G') == 0:
                    newgene = rje.choice('New Gene?:')
                    for node in _tree._nodeClade(fam):
                        seq = node.obj['Sequence']
                        #_tree.verbose(0,0,'%s => %s (gnspacc=%s)' % (seq.shortName(),seq.info['Format'],_tree.obj['SeqList'].opt['GeneSpAcc']),1)
                        seq.newGene(gene=newgene,keepsp=True,gnspacc=_tree.obj['SeqList'].opt['GeneSpAcc'])
                        node.info['Name'] = seq.info['Name']
                    node.info['CladeName'] = newgene
                ### UnGroup
                elif choice.find('U') == 0 and rje.yesNo('Remove this Group (Orphan sequences)?'):
                    _tree.subfam.remove(fam)
                ### Delete Group
                elif choice.find('D') == 0 and rje.yesNo('Delete this Group and all its sequences?'):
                    _tree.subfam.remove(fam)
                    _tree._prune(branch=fam.ancBranch(),remtext='Group Deleted')
                ### Next Group
                elif choice[0] == 'N' and g < (_tree.groupNum() - 1):
                    g += 1
                ### Previous Group
                elif choice.find('P') == 0:
                    g -= 1
                ### Finish
                elif (choice.find('Q') == 0  or (choice[0] == 'F' and g == (_tree.groupNum() - 1))) and (_tree.stat['Interactive'] < 0 or rje.yesNo('\nAccept Groups?')):
                    if _tree.stat['Interactive'] >= 0:
                        logtext = rje.choice('Text for log file? (Blank to skip): ')
                        if logtext != '':
                            _tree.log.printLog('#NB',logtext)
                    if _tree.stat['Interactive'] >= 0 and rje.yesNo('Save current sequence alignment and tree for backup?',default='N'):
                        fbase = ''
                        while fbase == '':
                            fbase = rje.choice('Base FILE name for output files, FILE.fas and FILE.nsf?: ')
                        _tree.obj['SeqList'].saveFasta(seqs=_tree.obj['SeqList'].seq[0:int(_tree.stat['SeqNum'])],seqfile='%s.fas' % fbase)
                        _tree.saveTree(filename='%s.nsf' % fbase)
                    else:
                        _tree.verbose(0,3,'No backups saved. ',0)
                    _tree.verbose(0,2,'Group Review Finished!',2)
                    return
                ### Group Menu
                elif choice.find('M') == 0 and rje.yesNo('Quit Review to Grouping Menu?'):
                    treeGroup(_tree,callmenu=True)
                    return _reviewGroups(_tree,interactive=interactive)
                break
            ## <h> ## Next Group or Finish
            _stage ='<1h> Next/finish groups'
            if _tree.groupNum() == 0:
                _tree.log.printLog('#GRP','All groups removed in Review Stage.')
            if g < 0:
                g = 0
            elif g >= _tree.groupNum():
                g = _tree.groupNum() - 1
    except:
        _tree.log.errorLog("Major problem during rje_tree_group._reviewGroups(%s) of %d groups." % (_stage,_tree.groupNum()), True)
#########################################################################################################################
def _groupDisSum(_tree,seq1,seq2,text=''):    ### Prints ID, Gaps and Extra Summary of seq1 vs seq2
    '''
    Prints MSA ID, MSA Gaps and MSA Extra Summary of seq1 vs seq2
    >> seq1:Sequence Object to be compared to...
    >> seq2:Sequence Object
    >> text:str  = vs text ['vs seq2.shortName()']
    << returns string
    '''
    if seq2 == None:
        return '%s:\t%d aa;' % (seq1.shortName(),seq1.aaLen())
    if text == '':
        text = 'vs %s' % seq2.shortName()
    seqlist = _tree.obj['SeqList']
    dis = (seqlist.getDis(seq1,seq2,'MSA ID'),seqlist.getDis(seq1,seq2,'MSA Gaps'),seqlist.getDis(seq1,seq2,'MSA Extra'))
    aas = (int(0.5+(seq1.aaLen()*dis[0]/100)),int(0.5+(seq1.aaLen()*dis[1]/100)),int(0.5+(seq1.aaLen()*dis[2]/100)))
    return '%s:\t%.2f%% ID (%d aa); %.2f%% Gaps (%d aa); %.2f%% Extra (%d aa); %s' % (seq1.shortName(),dis[0],aas[0],dis[1],aas[1],dis[2],aas[2],text)
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: 'MAIN' PROGRAM                                                                                         #
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        print 'This module is not for standalone running.'
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
