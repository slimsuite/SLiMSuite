#!/usr/local/bin/python

# BADASP - Burst After Duplication with Ancestral Sequence Prediction
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      BADASP
Description:  Burst After Duplication with Ancestral Sequence Prediction
Version:      1.3.1
Last Edit:    28/03/15
Citation:     Edwards & Shields (2005), Bioinformatics 21(22):4190-1. [PMID: 16159912]
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    BADASP implements the previously published Burst After Duplication (BAD) algorithm, plus two variants that have been used
    successfully in identifying functionally interesting sites in platelet signalling proteins and can identify Type I and
    Type II divergence. In addition, several other measures of functional specificity and conservation are calculated and
    output in plain text format for easy import into other applications. See Manual for details.

Commandline:
    # General Dataset Input/Output Options #
    seqin=FILE  : Loads sequences from FILE
    query=X     : Selects query sequence by name (or part of name, e.g. Accession Number)
    basefile=X  : Basic 'root' for all files X.* [By default will use 'root' of seqin=FILE if given or haq_AccNum if qblast]
    v=X         : Sets verbosity (-1 for silent) [0]
    i=X         : Sets interactivity (-1 for full auto) [0]
    log=FILE    : Redirect log to FILE [Default = calling_program.log or basefile.log]
    newlog=T/F  : Create new log file. [Default = False: append log file]
    rank=T/F    : Whether to output ranks as well as scores [True]
    append=FILE : Append results to FILE instead of standard output to *.badasp
    trimtrunc=T/F   : Whether to trim the leading and trailing gaps (within groups) -> change to X [False]
    winsize=X       : Window size for window scores

    # BADASP Statistics #
    funcspec=X1[,X2,..] : List of functional specificity methods to apply X1,X2,..,XN
        - BAD   = Burst After Duplication (2 Subfamilies)
        - BADX  = Burst After Duplication Extra (Query Subfam versus X subfams)
        - BADN  = Burst After Duplication vs N Subfams (2+ Subfams)
        - SSC   = Livingstone and Barton Score
        - PDAD  = Variant of Livingstone and Barton
        - ETA   = Evolutionary Trace Analysis (Basic)
        - ETAQ  = Evolutionary Trace Analysis (Quantitative)
        - all   = All of the above!
    seqcon=X1[,X2,..] : List of sequence conservation measures to apply X1,X2,..,XN
        - Info  = Information content
        - PCon  = Property Conservation (Absolute)
        - MPCon = Mean Property Conservation
        - QPCon = Mean Property Conservation with Query
        - all   = All of the above

    # Tree and Grouping Options #
    nsfin=FILE  : File name for Newick Standard Format tree
    root=X      : Rooting of tree (rje_tree.py), where X is:
        - mid = midpoint root tree. [Default]
        - ran = random branch.
        - ranwt = random branch, weighted by branch lengths.
        - man = always ask for rooting options (unless i<0).
        - FILE = with seqs in FILE as outgroup. (Any option other than above)
    bootcut=X   : cut-off percentage of tree bootstraps for grouping.
    mfs=X       : minimum family size [3]
    fam=X       : minimum number of families (If 0, no subfam grouping) [0]
    orphan=T/F  : Whether orphans sequences (not in subfam) allowed. [True]
    allowvar=T/F: Allow variants of same species within a group. [False]
    gnspacc=T/F : Convert sequences into gene_SPECIES__AccNum format wherever possible. [True] 
    groupspec=X : Species for duplication grouping [None]
    group=X     : Grouping of tree
        - man = manual grouping (unless i<0).
        - dup = duplication (all species unless groupspec specified).
        - qry = duplication with species of Query sequence (or Sequence 1) of treeseq
        - one = all sequences in one group
        - None = no group (case sensitive)
        - FILE = load groups from file

    # GASP ancestral sequence prediction options #
    useanc=FILE : Gives file of predicted ancestral sequences
    pamfile=FILE: Sets PAM1 input file [jones.pam]
    pammax=X    : Initial maximum PAM matrix to generate [100]
    pamcut=X    : Absolute maximum PAM matrix [1000]
    fixpam=X    : PAM distance fixed to X [100].
    rarecut=X   : Rare aa cut-off [0.05].
    fixup=T/F   : Fix AAs on way up (keep probabilities) [True].
    fixdown=T/F : Fix AAs on initial pass down tree [False].
    ordered=T/F : Order ancestral sequence output by node number [False].
    pamtree=T/F : Calculate and output ancestral tree with PAM distances [True].
    desconly=T/F: Limits ancestral AAs to those found in descendants [True].
    xpass=X     : How many extra passes to make down & up tree after initial GASP [1].

    # System Info Options #
    win32=T/F       : Run in Win32 Mode [False]

Please see help for rje_tree.py and rje_seq.py for additional options not covered here.    

Uses general modules: copy, string, sys, time
Uses RJE modules: rje, rje_aaprop, rje_ancseq, rje_conseq, rje_seq, rje_specificity, rje_tree
Additional modules required: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_tree_group, rje_uniprot
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import copy
import os
import string
import sys
import time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
#############################################################################################################################
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_aaprop
import rje_ancseq
import rje_conseq
import rje_seq
import rje_specificity
import rje_tree
#############################################################################################################################
### History
# 0.0 - Initial Development Compilation.
# 1.0 - Intial working version
# 1.1 - Added output options
# 1.2 - Neatened and updated
# 1.3 - Recognise only one subfamily and deal with accordingly.
# 1.3.1 - Modified to recognise *.nwk as well as *.nsf tree input.
#############################################################################################################################
### Major Functionality to Add
# [ ] Improve menu in general and choice of output in particular
# [y] Option to load ancestral sequences (in rje_tree?)
# [ ] Turn BADASP into an Object and neaten up.
#############################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'BADASP'
        version = '1.3.1'
        last_edit = 'March 2015'
        description = "Burst After Duplication with Ancestral Sequence Prediction"
        author = "Dr Richard J. Edwards."
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
            if rje.yesNo('Show sequence commandline options?'):
                out.verbose(-1,4,text=rje_seq.__doc__)
            if rje.yesNo('Show tree commandline options?'):
                out.verbose(-1,4,text=rje_tree.__doc__)
            if rje.yesNo('Show ancestral sequence prediction commandline options?'):
                out.verbose(-1,4,text=rje_ancseq.__doc__)
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
## End of SECTION II
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                           #
#############################################################################################################################
def badasp(out,mainlog,cmd_list,tree=None): ### Main BADASP Method
    '''
    Main BADASP Method. Automated run if interactive < 1
    <1> Load Sequences and Tree
    <2> Define Subfamilies
    <3> GASP Ancestral Sequence Prediction
    <4> Peform Functional Specificity and Sequence Conservation Calculations
    <5> Output Results
    '''
    try:    ### <0> ### Setup
        _seqfile = None
        _treefile = None
        append_file = None
        basefile = None
        for cmd in cmd_list:
            if cmd.find('seqin=') == 0:
                _seqfile = cmd[len('seqin='):]
                if _seqfile[-4] == '.':
                    _seqfile = _seqfile[:-4]
            if cmd.find('useanc=') == 0:
                _seqfile = cmd[len('useanc='):]
                if _seqfile[-8:] == '.anc.fas':
                    _seqfile = _seqfile[:-8]
            if cmd.find('nsfin=') == 0:
                _treefile = cmd[len('nsfin='):]
            if cmd.find('append=') == 0:
                append_file = cmd[len('append='):]
            if cmd.find('basefile=') == 0:
                basefile = cmd[len('basefile='):]
        if _seqfile and os.path.exists('%s.grp' % _seqfile):
            cmd_list.append('group=%s.grp' % _seqfile)
        if _seqfile and _treefile == None:
            if rje.checkForFile('%s.nwk' % _seqfile): _treefile = '%s.nwk' % _seqfile
            else: _treefile = '%s.nsf' % _seqfile
            out.verbose(0,2,'Looking for treefile %s.' % _treefile,1)
            if rje.checkForFile(_treefile):
                cmd_list.append('nsfin=%s' % _treefile)
                
        if tree == None:
            mainlog.verbose(0,1,'Tree: %s' % cmd_list,2)
            tree = rje_tree.Tree(log=mainlog,cmd_list=cmd_list)
            #tree._setupFromCmd()
        if tree.stat['MinFamNum'] < 2:
            tree.stat['MinFamNum'] = 2

            ### <1> ### Load Sequences and Tree
        while out.stat['Interactive'] > 0 or tree.obj['SeqList'] == None:
            tree = rje_tree.treeMenu(out,mainlog,['root=yes']+cmd_list,tree)
            if tree.obj['SeqList'] and tree.opt['Rooted']:
                break
            else:
                print '\n ** Must have loaded sequences and a rooted tree. ** \n'
                if out.stat['Interactive'] < 0 or rje.yesNo('Quit BADASP?',default='N'):
                    sys.exit()

        basename = tree.obj['SeqList'].info['Name']
        if basename[-4:] == '.fas':
            basename = basename[:-4]
        if basename[-4:] == '.anc':
            basename = basename[:-4]
        if basefile:
            basename = basefile
                                
    except SystemExit:
        raise
    except:
        mainlog.errorLog('Major Error in badasp loading sequences and tree',True)

    try:    ### <2> ### Define Subfamilies
        while out.stat['Interactive'] > 0 or tree.groupNum() < 2:
            tree.treeGroup(callmenu=True)
            if tree.groupNum() >= 2:
                break
            else:
                mainlog.errorLog('Must have at least two subfamilies for specificity analyses.',printerror=False)
                if out.stat['Interactive'] < 0 or rje.yesNo('Continue without specificity analyses?'):
                    cmd_list.append('funcspec=')
                    break
                elif rje.yesNo('Abort BADASP?'):
                    sys.exit()
    except SystemExit:
        raise
    except:
        mainlog.errorLog('Major Error in BADASP subfamilies',True)

    try:    ### <3> ### GASP Ancestral Sequence Prediction
        if tree.node[-1].obj['Sequence'] == None:   # No ancseq loaded
            while out.stat['Interactive'] > 0 and rje.yesNo('Use %s for output filenames?' % basename) == False:
                basename = rje.choice('FILEname (FILE.anc.fas, FILE.anc.nsf, FILE.txt)?: ', default=basename)
            mygasp = rje_ancseq.Gasp(tree=tree,ancfile=basename,cmd_list=cmd_list,log=mainlog)
            out.verbose(0,2,'%s' % mygasp.details(),1)
            if out.stat['Interactive'] > 0:
                if rje.yesNo('Use these parameters?') == False:
                    mygasp.edit()
            mygasp.gasp()
    except:
        mainlog.errorLog('Major Error in BADASP GASP',True)
        
    try:    ### <4> ### Peform Functional Specificity and Sequence Conservation Calculations
        _stage = '<4> Specificity/Conservation Analyses'
        aaprop = rje_aaprop.AAPropMatrix(log=mainlog,cmd_list=cmd_list)
        query = tree.obj['SeqList'].obj['QuerySeq']
        ## <a> ## Chosen Methods
        _stage = '<4a> Specificity/Conservation Analyses - Chosen Methods'
        funcspec = rje_specificity.methodlist   # ['BAD','BADN','BADX']
        seqcon = rje_conseq.methodlist # ['info']
        for cmd in cmd_list:
            if cmd.find('funcspec=') == 0:
                funcspec = cmd[9:].split(',')
            if cmd.find('seqcon=') == 0:
                seqcon = cmd[len('seqcon='):].split(',')
        if 'all' in funcspec:
            funcspec = rje_specificity.methodlist
        if 'all' in seqcon:
            seqcon = rje_conseq.methodlist
        for method in ['BADX','BADN','QPCon_Mean','QPCon_Abs','QPCon_Mean_All']:
            while method in funcspec and query == None:
                if rje.yesNo('Method %s needs query but none given. Drop %s from specificity methods?' % (method,method)):
                    funcspec.remove(method)
                    break
                for seq in tree.obj['SeqList'].seq:
                    if rje.yesNo('Method %s needs query but none given. Use sequence 1 (%s)?' % (method,seq.shortName()),default='N'):
                        query = seq
                        tree.obj['SeqList'].obj['Query'] = seq
                        break
            while method in seqcon and query == None:
                if rje.yesNo('Method %s needs query but none given. Drop %s from conservation methods?' % (method,method)):
                    seqcon.remove(method)
                    break
                for seq in tree.obj['SeqList'].seq:
                    if rje.yesNo('Method %s needs query but none given. Use sequence 1 (%s)?' % (method,seq.shortName()),default='N'):
                        query = seq
                        tree.obj['SeqList'].obj['Query'] = seq
                        break
                
        qname = query
        if query:
            qname = query.info['Name']
        out.verbose(0,3,'\nQuery = %s' % qname,2)

        ## <b> ## Spec Calculations
        _stage = '<4b> Specificity Calculations'
        specmatrix = rje_specificity.FuncSpec(log=mainlog,cmd_list=cmd_list,tree=tree,aaprop=aaprop)
        specmatrix.calcScore(query=query,methods=funcspec)

        ## <c> ## Conservation Calculations
        _stage = '<4c> Specificity/Conservation Analyses - Conservation Calculations'
        conseq = rje_conseq.SeqStat(log=mainlog,cmd_list=cmd_list,tree=tree,aaprop=aaprop)
        conseq.calcScore(query=query,methods=seqcon)   ### Sends appropriate seqlist to self.calcScore()

        ## <d> ## Special Case: QPCon vs All seqs
        _stage = '<4d> Specificity/Conservation Analyses - QPCon vs All'
        qpconall = []
        #if 'QPCon_Abs_All' in seqcon and query:
        #    qpconall.append('QPCon_Abs')
        if 'QPCon_Mean_All' in seqcon and query:
            qpconall.append('QPCon_Mean')
        for qp in qpconall:
            conseq.score['%s_All' % qp] = conseq.score[qp] 
            if conseq.alnwin.has_key(qp):
                conseq.alnwin['%s_All' % qp] = conseq.alnwin[qp] 
            if conseq.qrywin.has_key(qp):
                conseq.qrywin['%s_All' % qp] = conseq.qrywin[qp] 
            if conseq.rank.has_key(qp):
                conseq.rank['%s_All' % qp] = conseq.rank[qp] 
            if conseq.alnrankwin.has_key(qp):
                conseq.alnrankwin['%s_All' % qp] = conseq.alnrankwin[qp] 
            if conseq.qryrankwin.has_key(qp):
                conseq.qryrankwin['%s_All' % qp] = conseq.qryrankwin[qp] 
        
        _stage = '<4d> Specificity/Conservation Analyses - FamQP'
        famqp = []
        if 'QPCon_Mean' in seqcon:
            famqp.append('QPCon_Mean')
        if 'QPCon_Abs' in seqcon:
            famqp.append('QPCon_Abs')
        if len(famqp) > 0 and query:    #!# And subfam option?
            qseq = []
            for fam in tree.subfam:
                for node in tree._nodeClade(fam):
                    if query == node.obj['Sequence']:
                        for qnode in tree._nodeClade(fam):
                            qseq.append(qnode.obj['Sequence'])
            conseq.calcScore(query=query,seqlist=qseq,methods=famqp)   ### Sends appropriate seqlist to self.calcScore()

    except:
        mainlog.errorLog('Major Error in BADASP Specificity Analysis (%s):' % _stage,True)
        
    try:    ### <5> ### Full Output Results
        _stage = '<5> Full Output'
        # This output is in a tab- or comma-delimited file for easy manipulation or viewing with other programs.
        # (1) statistics for a given residue;
        # (2) statistics for a given window size across
        # - (a) the whole alignment,    (node=None)
        # - (b) the Query protein of interest (if given) and    (node=QueryNode)
        # - (c) the ancestral sequence of each subfamily;   (node=ancnode)
        # (3) Predicted ancestral sequences at
        # - (a) the root and
        # - (b) the ancestor of each subfamily.
        delimit = rje.getDelimit(cmd_list)

        ## <a> ## Setup        
        _stage = '<5a> Output - Setup'
        rankout = specmatrix.opt['Rank']
        #tree._regenerateSeqList(tree.obj['SeqList'],tree.node)
        root = tree.node[-1].obj['Sequence']      #!# At some point, make sure this is the most ancient duplication!
        out.verbose(0,3,'\nBADASP Results Output (%s.badasp) ...' % basename,0)

        ## <b> ## Header
        _stage = '<5b> Output - Header'
        _header = True
        if append_file:
            if rje.checkForFile(append_file):
               _header = False
            BADASP = open(append_file, 'a')
        else:
            BADASP = open('%s.badasp' % basename, 'w')
            BADASP.write("BADASP Output: %s\n" % (time.asctime(time.localtime(time.time()))))
            BADASP.write('%s\n\n' % cmd_list)        
        header = ['aln_pos','anc_aa'] # Aln Pos and AA
        alnlen = 0
        statlist = funcspec + seqcon
        _stage = '<5b-i> Output - Header Query'
        if query:
            header += ['qry_pos','qry_aa']  # Qry Pos and AA
        _stage = '<5b-ii> Output - Header Subfam'
        for f in range(len(tree.subfam)):
            header += ['fam%d_pos' % (f+1),'fam%d_aa' % (f+1)]   # Subfam Pos and AA
        for func in statlist:
            _stage = '<5b-iii> Output - Header %s' % func
            statobj = statObj(method=func,objlist=[specmatrix,conseq])
            fs = func.lower()
            alnlen = len(statobj.score[func])                
            header.append(fs)                   # Score
            if rankout:
                header.append('%s_rank' % fs)   # Rank
            if statobj.stat['WinSize'] > 1:
                header.append('%s_alnwin' % fs)     # Full align window
                if rankout:
                    header.append('%s_alnrankwin' % fs)   # Rank
                if query:
                    header.append('%s_qrywin' % fs) # Qry window
                    if rankout:
                        header.append('%s_qryrankwin' % fs)   # Rank
                if func in funcspec:
                    for f in range(len(tree.subfam)):
                        header.append('%s_fam%d_win' % (fs,f+1))  # Subfam windows
                        if rankout:
                            header.append('%s_fam%d_rankwin' % (fs,f+1))  # Subfam windows                    
        #if _header:
        BADASP.write('%s\n' % string.join(header, delimit))
        out.verbose(1,3,'%s...' % string.join(header, delimit),0)

        ## <c> ## Stats
        _stage = '<5c> Stats'
        qr = 0  # Qry pos
        fr = [0] * len(tree.subfam) # List of subfam positions
        aa = '' # Root aa
        qa = '' # Qry aa
        fa = [''] * len(tree.subfam) # List of subfam aas
        for r in range(alnlen):
            # <i> # Positions and aas
            _stage = '<5c-i> Output - Stats, positions & aas'
            aa = root.info['Sequence'][r]
            if query:
                qa = query.info['Sequence'][r]
                if qa != '-':
                    qr += 1
            for f in range(len(tree.subfam)):
                fa[f] = tree.subfam[f].obj['Sequence'].info['Sequence'][r]
                if fa[f] != '-':
                    fr[f] += 1

            # <ii> # Positions and AAs ii
            _stage = '<5c-ii> Output - Pos & AA ii'
            line = ['%d' % (r+1), aa]    # Aln Pos and AA
            if query:
                if qa == '-':               
                    line += ['-',qa]  # Qry Pos and AA
                else:               
                    line += ['%d' % qr,qa]  # Qry Pos and AA
            for f in range(len(tree.subfam)):
                if fa[f] == '-':
                    line += ['-',fa[f]]   # Subfam Pos and AA
                else:
                    line += ['%d' % fr[f],fa[f]]   # Subfam Pos and AA

            # <iii> # Stats                        
            _stage = '<5c-iii> Output - Stats'
            for func in statlist:
                statobj = statObj(method=func,objlist=[specmatrix,conseq])
                fs = func.lower()
                line.append(str(statobj.score[func][r]))   # Score
                if rankout:
                    line.append(str(statobj.rank[func][r]))   # Rank
                if specmatrix.stat['WinSize'] > 1:
                    line.append(str(statobj.alnwin[func][r]))     # Full align window
                    if rankout:
                        line.append(str(statobj.alnrankwin[func][r]))   # Rank
                    if query:
                        line.append(str(statobj.qrywin[func][r])) # Qry window
                        if rankout:
                            line.append(str(statobj.qryrankwin[func][r]))   # Rank
                    if func in funcspec:
                        for f in range(len(tree.subfam)):
                            line.append(str(statobj.famwin[func][tree.subfam[f]][r]))  # Subfam windows
                            if rankout:
                                line.append(str(statobj.famrankwin[func][tree.subfam[f]][r]))   # Subfam windows                    
            # <iv> # Writing
            _stage = '<5c-iv> Output - Writing'
            BADASP.write('%s\n' % string.join(line, delimit))
        BADASP.close()
        out.verbose(0,2,'Done!',2)

    except:
        mainlog.errorLog('Fatal Error in BADASP Full output (%s):' % _stage,True)
        BADASP.write('%s\n' % string.join(line, delimit))
        BADASP.close()
                        
    try:    ### <6> ### Partial Results Output 
        _stage = '<6> Partial Output'

        ## <a> ## Setup        
        _stage = '<6a> Output - Setup'
        # statlist & alnlen from above
        _part_append = False
        if out.stat['Interactive'] > 0 and rje.yesNo('Output additional, filtered results?',default='N'):
            partfile = rje.choice('Name for partial results file?:','%s.partial.badasp' % basename,confirm=True)
            if rje.checkForFile(partfile) and rje.yesNo('File %s exists. Append file without headers?' % partfile):
                _part_append = True
        else:
            return
        if rje.yesNo('Filter output columns?',default='N'):
            if rje.yesNo('Output query details (pos,aa & win)?') == False:
                query = None
            f = 1
            for fam in tree.subfam[0:]:
                if rje.yesNo('Output subfam %d (%s) details (pos,aa & win)?' % (f,fam.info['CladeName'])) == False:
                    tree.subfam.remove(fam)
                f += 1
            for func in statlist[0:]:
                if rje.yesNo('Output %s results?' % func) == False:
                    statlist.remove(func)

        alnout = [True] * alnlen
        if rje.yesNo('Filter Rows by Results VALUES?'):
            out.verbose(0,0,'Initial Defaults are minmum values. Accept intital default for no filtering of given Stat.',1)
            for stat in statlist:
                ### Filter by value? ###
                statobj = statObj(method=stat,objlist=[specmatrix,conseq])
                scores = statobj.score[stat][0:]
                scores.sort()
                cutoff = rje.getFloat('Min. value for %s?:' % stat,default='%f' % scores[0],confirm=True)
                for r in range(alnlen):
                    if statobj.score[stat][r] < cutoff:
                        alnout[r] = False
        if rankout and rje.yesNo('Filter Rows by Results RANKS?'):
            out.verbose(0,0,'Ranks range from 0 (low) to 1 (high).',1)
            for stat in statlist:
                ### Filter by Rank? ###
                statobj = statObj(method=stat,objlist=[specmatrix,conseq])
                cutoff = rje.getFloat('Min. rank for %s?:' % stat,default='0.0',confirm=True)
                for r in range(alnlen):
                    if statobj.rank[stat][r] < cutoff:
                        alnout[r] = False
                
        out.verbose(0,3,'\nBADASP Partial Results Output (%s) ...' % partfile,0)

        ## <b> ## Header
        _stage = '<6b> Partial Output - Header'
        if _part_append:
            BADASP = open(partfile, 'a')
        else:
            BADASP = open(partfile, 'w')
            BADASP.write("Partial BADASP Output: %s\n" % (time.asctime(time.localtime(time.time()))))
            BADASP.write('%s\n\n' % cmd_list)        
        header = ['aln_pos','anc_aa'] # Aln Pos and AA
        _stage = '<6b-i> Partial Output - Header Query'
        if query:
            header += ['qry_pos','qry_aa']  # Qry Pos and AA
        _stage = '<6b-ii> Partial Output - Header Subfam'
        for f in range(len(tree.subfam)):
            header += ['fam%d_pos' % (f+1),'fam%d_aa' % (f+1)]   # Subfam Pos and AA
        for func in statlist:
            _stage = '<6b-iii> Partial Output - Header %s' % func
            statobj = statObj(method=func,objlist=[specmatrix,conseq])
            fs = func.lower()
            header.append(fs)                   # Score
            if rankout:
                header.append('%s_rank' % fs)   # Rank
            if statobj.stat['WinSize'] > 1:
                header.append('%s_alnwin' % fs)     # Full align window
                if rankout:
                    header.append('%s_alnrankwin' % fs)   # Rank
                if query:
                    header.append('%s_qrywin' % fs) # Qry window
                    if rankout:
                        header.append('%s_qryrankwin' % fs)   # Rank
                if func in funcspec:
                    for f in range(len(tree.subfam)):
                        header.append('%s_fam%d_win' % (fs,f+1))  # Subfam windows
                        if rankout:
                            header.append('%s_fam%d_rankwin' % (fs,f+1))  # Subfam windows                    
        #if not _part_append:
        BADASP.write('%s\n' % string.join(header, delimit))
        out.verbose(1,3,'%s...' % string.join(header, delimit),0)

        ## <c> ## Stats
        _stage = '<6c> Stats'
        qr = 0  # Qry pos
        fr = [0] * len(tree.subfam) # List of subfam positions
        aa = '' # Root aa
        qa = '' # Qry aa
        fa = [''] * len(tree.subfam) # List of subfam aas
        for r in range(alnlen):
            if alnout[r] == False:
                continue
            # <i> # Positions and aas
            _stage = '<6c-i> Partial Output - Stats, positions & aas'
            aa = root.info['Sequence'][r]
            if query:
                qa = query.info['Sequence'][r]
                if qa != '-':
                    qr += 1
            for f in range(len(tree.subfam)):
                fa[f] = tree.subfam[f].obj['Sequence'].info['Sequence'][r]
                if fa[f] != '-':
                    fr[f] += 1

            # <ii> # Positions and AAs ii
            _stage = '<6c-ii> Partial Output - Pos & AA ii'
            line = ['%d' % (r+1), aa]    # Aln Pos and AA
            if query:
                if qa == '-':               
                    line += ['-',qa]  # Qry Pos and AA
                else:               
                    line += ['%d' % qr,qa]  # Qry Pos and AA
            for f in range(len(tree.subfam)):
                if fa[f] == '-':
                    line += ['-',fa[f]]   # Subfam Pos and AA
                else:
                    line += ['%d' % fr[f],fa[f]]   # Subfam Pos and AA

            # <iii> # Stats                        
            _stage = '<6c-iii> Partial Output - Stats'
            for func in statlist:
                statobj = statObj(method=func,objlist=[specmatrix,conseq])
                fs = func.lower()
                line.append(str(statobj.score[func][r]))   # Score
                if rankout:
                    line.append(str(statobj.rank[func][r]))   # Rank
                if specmatrix.stat['WinSize'] > 1:
                    line.append(str(statobj.alnwin[func][r]))     # Full align window
                    if rankout:
                        line.append(str(statobj.alnrankwin[func][r]))   # Rank
                    if query:
                        line.append(str(statobj.qrywin[func][r])) # Qry window
                        if rankout:
                            line.append(str(statobj.qryrankwin[func][r]))   # Rank
                    if func in funcspec:
                        for f in range(len(tree.subfam)):
                            line.append(str(statobj.famwin[func][tree.subfam[f]][r]))  # Subfam windows
                            if rankout:
                                line.append(str(statobj.famrankwin[func][tree.subfam[f]][r]))   # Subfam windows
            # <iv> # Writing
            _stage = '<6c-iv> Partial Output - Writing'
            BADASP.write('%s\n' % string.join(line, delimit))
        BADASP.close()
        out.verbose(0,2,'Done!',2)
        
    except:
        mainlog.errorLog('Fatal Error in BADASP Partial output (%s):' % _stage,True)
        BADASP.write('%s\n' % string.join(line, delimit))
        BADASP.close()
#############################################################################################################################
def statObj(method=None,objlist=[]):    ### Returns the appropriate Object for scores of 'method'
    '''
    Returns the appropriate Object for scores of 'method'.
    >> method:str = Method name
    >> objlist:list of rje_conseq.SeqStat type objects
    '''
    for statobj in objlist:
        if method in statobj.score.keys():
            return statobj
    return None
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
        badasp(out,mainlog,cmd_list)

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
