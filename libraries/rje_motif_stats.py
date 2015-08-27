#!/usr/local/bin/python

# Motif Statistics Methods Module
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_motif_stats
Description:  Motif Statistics Methods Module
Version:      1.0
Last Edit:    01/02/07
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the Alignment Conservation methods for motifs, as well as other calculations needing occurrence
    data. This module is designed to be used by the MotifList class, which contains the relevant commandline options.

Commandline:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module.

Uses general modules: copy, os, string, sys
Uses RJE modules: gopher_V2, rje, rje_blast, rje_disorder, rje_motif_V3, rje_seq, rje_sequence
Other modules needed: rje_seq modules
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, string, sys
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../legacy/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import gopher_V2, rje, rje_disorder, rje_motif_V3, rje_seq, rje_sequence
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Initial compilation based of rje_motif_cons V1.1.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Check full functionality with new modules.
    '''
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: General Occurrence Stats                                                                                #
#########################################################################################################################
def statLog(callobj,logtext,extra): ### Output to log
    '''Output to log.'''
    if len(logtext) == 2:
        callobj.log.printLog(logtext[0],logtext[1]+extra,log=False,newline=False)
#########################################################################################################################
def occStats(callobj,occlist,logtext=()):     ### Calculates general occurrence stats for occlist
    '''
    Calculates general occurrence stats for occlist.
    >> callobj:Object containing settings for stats generation (MotifList, generally).
    >> occlist:list of MotifOcc objects to calculate stats for (must all have same Seq)
    >> logtext:tuple of (logleader,logtext) to use as basis for log messages. (None if len != 2)
    '''
    try:
        ### Setup ###
        ## OccList ##
        if not occlist:
            return
        Seq = occlist[0].obj['Seq']
        ## Sequence information, inc. XPad ##
        Seq.deGap()
        xpaddb = callobj.getStat('XPadDB')
        sequence = 'X' * xpaddb + Seq.info['Sequence'] + 'X' * xpaddb
        seq_sa = rje_seq.surfaceAccessibility(sequence,returnlist=True)
        seq_eis = rje_seq.eisenbergHydropathy(sequence,returnlist=True)
        seq_dis = {}    # Dictionary of lists of scores for different disorder predictions
        seq_dom = {}    # Dictionary of lists of scores for different domain filters
        ## Disorder ##
        if callobj.opt.has_key('SlimIUP') and callobj.opt['SlimIUP']:
            callobj.opt['IUPred'] = True
        if callobj.opt.has_key('SlimFold') and callobj.opt['SlimFold']:
            callobj.opt['FoldIndex'] = True
        for o in ['IUPred','FoldIndex']:
            if callobj.opt[o]:
                callobj.log.stat['Verbose'] -= 2
                dcmd = callobj.cmd_list+['disorder=%s' % o.lower(),'sequence=%s' % sequence,'v=%d' % callobj.log.stat['Verbose']]
                disorder = rje_disorder.Disorder(log=callobj.log,cmd_list=dcmd)
                if o == 'FoldIndex' or disorder.stat['IUCut'] > 0.0:
                    disorder.flatten()   # Converts to 1/0 #
                callobj.log.stat['Verbose'] += 2
                seq_dis[o] = disorder.list['ResidueDisorder'][0:]
        seq_dom = seqDom(callobj,Seq,seq_dis)
        
        ### Stats ###
        ox = 0.0
        for Occ in occlist:
            statLog(callobj,logtext,' Stats %.1f%%' % (ox/len(occlist)))
            ox += 100.0
            r = int(Occ.stat['Pos'] - 1) + xpaddb
            if Occ.stat['Pos'] < 1:     # XPadDB
                r += 1
            match = Occ.getData('Match')
            ## SA, Hyd, Charge and disorder ##
            fundict = {'SA':seq_sa,'Hyd':seq_eis,'Chg':sequence}
            for fun in ['SA','Hyd','Chg','Dis']:
                ## Skip? ##
                if fun == 'Chg' and not callobj.getOpt('SlimChg'):
                    continue
                ## Set window ##
                win = int(callobj.getStat('Win%s' % fun))    # Will return 0 if missing
                if win >= 0:
                    w = r - win
                else:
                    w = r + win
                if w < 0:
                    w = 0
                ## Get appropriate Data region ##
                if fundict.has_key(fun):
                    if win >= 0:
                        winreg = fundict[fun][w:r+len(match)+win]     # Region includes occurrence
                    else:   
                        winreg = fundict[fun][w:r] + fundict[fun][r+len(match):r+len(match)-win]    # Flanks only
                    if fun == 'Chg':
                        Occ.setStat(rje_sequence.chargeDict(winreg))
                        #X#Occ.stat['CHG_ABS'] = chgdict['AbsChg']      #!# Changed, so update PRESTO Manual #!#
                        #X#Occ.stat['CHG_NET'] = chgdict['NetChg']
                        #X#Occ.stat['CHG_BAL'] = chgdict['BalChg']
                    else:                        
                        Occ.stat[fun] = 0   # Note that the stat is now Hyd not Hydropathy!
                        if winreg:
                            Occ.stat[fun] = sum(winreg) / len(winreg)
                elif fun == 'Dis':
                    for o in ['IUPred','FoldIndex']:    #!# Standardise this for SLiMPickings and PRESTO?! (IUP/FI?) #!#
                        if callobj.opt[o]:
                            if win >= 0:
                                winreg = seq_dis[o][w:r+len(match)+win]     # Region includes occurrence
                            else:
                                winreg = seq_dis[o][w:r] + seq_dis[o][r+len(match):r+len(match)-win]   # Flanks only
                            Occ.stat[o] = 0                              
                            if winreg:
                                Occ.stat[o] = sum(winreg) / len(winreg)
                            Occ.stat[o[:int(len(o)/2)]] = Occ.stat[o]

            ## Domain Filtering ##
            Occ.stat['MASK'] = max(seq_dom['MASK'][r:r+len(match)])
            Occ.stat['PROP'] = sum(seq_dom['MASK'][r:r+len(match)]) / len(seq_dom['MASK'][r:r+len(match)])
            Occ.stat['DIS'] = min(seq_dom['DIS'][r:r+len(match)])
            Occ.stat['COMB'] = sum(seq_dom['COMB'][r:r+len(match)]) / len(seq_dom['COMB'][r:r+len(match)])

            ## Peptide Design ##
            pepwin = rje.modulus(callobj.stat['WinSize'])
            w = r - pepwin
            if w < 0:
                w = 0
            Occ.info['PepSeq'] = sequence[w:r+len(match)+callobj.stat['WinSize']]     # Region for calculation
        statLog(callobj,logtext,' Stats %.1f%%' % (ox/len(occlist)))

        ### AlnCon ###  
        if callobj.opt['UseAln']:
            hitAlnCon(callobj,occlist)

    except:
        callobj.log.errorLog('Error in rje_motif_stats.occStats',quitchoice=True)  
#########################################################################################################################
def seqDom(callobj,seq,seq_dis):  ### Returns dictionary of lists of scores for different domain scores
    '''
    Returns dictionary of lists of scores for different domain scores.
    >> callobj:Object controlling attributes
    >> Seq:Sequence Object
    >> seq_dis:dictionary of lists of disorder scores
    << seq_dom:dicitonary of domain scores
    '''
    try:
        ### Setup ###
        seq_dom = {'MASK':[0.0] * seq.aaLen(),'DIS':[1.0] * seq.aaLen()}
        if seq_dis.has_key('IUPred'):
            seq_dom['COMB'] = seq_dis['IUPred'][0:]
        else:
            seq_dom['COMB'] = [0.0] * seq.aaLen()
        
        ### Domain Filter ###
        if not callobj.dict.has_key('DomFilter'):
            try:
                callobj.setupDomFilter()
            except:
                return seq_dom
        if callobj.dict.has_key('DomFilter') and callobj.dict['DomFilter'].has_key(seq.shortName()):
            ## MASK ##
            for dom in callobj.dict['DomFilter'][seq.shortName()]:
                (start,end) = dom
                for i in range(start-1,end):
                    seq_dom['MASK'][i] = float(end - start + 1)
                ## DIS ##
                if seq_dis.has_key('IUPred'):
                    (mean,se) = rje.meanse(seq_dis['IUPred'][start-1:end])
                    for i in range(start-1,end):
                        seq_dom['DIS'][i] = mean
                        seq_dom['COMB'][i] = mean
                    
        ### Finish ###
        return seq_dom
    except:
        callobj.log.errorLog('rje_motif_stats.seqDom() is not right!')
        return seq_dom
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: Conservation Methods                                                                                   #
#########################################################################################################################
def hitAlnCon(callobj,occlist,logtext=()):    ### Looks for alignment and, if appropriate, calculate conservation stats.
    '''
    Looks for alignment and, if appropriate, calculate conservation stats.
    
    Any homologues with masked (X) residues that coincide to non-wildcard positions of the motif occurrence will be
    ignored from conservation calculations. Gaps, however, shall be treated as divergence. The exception is that when the
    alngap=F option is used, 100% gapped regions of homologues are also ignored.

    This method deals with all the occurrences of all motifs for a single sequence and its alignment. Global alignment
    statistics are calculated first, then each occurrence for each motif is processed. Since version 1.1, subtaxa are
    treated the same as all taxa to reduce the coding: the default all taxa is now effectively an additional subtaxa set.
    
    >> callobj:Object containing settings for stats generation (MotifList, generally).
    >> occlist:list of MotifOcc objects to calculate stats for (must all have same Seq)
    >> logtext:tuple of (logleader,logtext) to use as basis for log messages. (None if len != 2)
    '''
    try:
        ### Setup ###
        if not occlist:
            return
        OccSeq = occlist[0].obj['Seq']

        ### Default Values (0.0) ###
        for Occ in occlist:
            Motif = Occ.obj['Motif']
            for taxa in ['ALL'] + rje.sortKeys(callobj.dict['ConsSpecLists']):
                Occ.stat['%s_HOM' % taxa] = 0
                Occ.stat['%s_GLOB_ID' % taxa] = 0.0
                Occ.stat['%s_LOC_ID' % taxa] = 0.0
                Occ.stat['%s_CONS' % taxa] = 1.0        #!# Default 1.0 if no homologues #!#
                if callobj.info['ConScore'] == 'all':
                    for method in ['ABS','POS','PROP']:
                        Occ.stat['%s_CONS_%s' % (taxa,method)] = 1.0    #!# Default 1.0 if no homologues #!#
            #X#callobj.deBug(Occ.stat)

        ### Identify and Load Alignment ###
        aln = loadOrthAln(callobj,OccSeq)
        if not aln:     # Alignment rejected during loadOrthAln. No conservation statistics here. #!# Check when GOPHER is called?! #!#
            return
        qry = aln.obj['QuerySeq']
        seqs = aln.seq[0:]      # This is a list of Sequence objects, minus the query, which is stored in qry
        seqs.remove(qry)
        callobj.log.printLog('#ALN','Comparing %d homologues for %s.' % (len(seqs),qry.shortName()),screen=False)   #!# Log? #!#

        ### Calculate Global Statistics ###
        ## Sequence-specific statistics ##
        seqid = {}      # Dictionary of global %ID vs query for each sequence (consweight<>0)
        seqwt = {}      # Dictionary of global %ID weighting for each sequence (consweight<>0)
        for seq in seqs:
            igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'],seq.info['Sequence'],nomatch=['X'])
            seqid[seq] = float(igedic['ID'][0]) / igedic['Len'][0]
            seqwt[seq] = seqid[seq] ** callobj.stat['ConsWeight']
        ## Taxa-specific Sequence lists ##
        conseqs = {'ALL':seqs[0:]}    # Dictionary of sequence lists {speclist:seqs}
        for taxa in rje.sortKeys(callobj.dict['ConsSpecLists']):
            conseqs[taxa] = []
            codelist = callobj.dict['ConsSpecLists'][taxa]
            for seq in seqs:
                if seq.info['SpecCode'] in codelist or seq.info['Species'] in codelist:
                    conseqs[taxa].append(seq)
        ## Taxa-specific Global ID ##
        globid = {}     # Dictionary of mean global %ID for each taxonomic grouping
        for taxa in conseqs.keys():
            globid[taxa] = 0.0
            for seq in conseqs[taxa]:
                globid[taxa] += seqid[seq]
            ## Average Global identity ##
            if len(conseqs[taxa]) > 0:
                globid[taxa] /= len(conseqs[taxa])

        ### Process Occurrences ###        
        for Occ in occlist:
            Motif = Occ.obj['Motif']

            ## Find position in aligment ##
            (start,end) = findOccPos(callobj,Occ,qry)      
            if (start,end) == (-1,-1): continue # Could not find!
                
            ## Homologues, Sequence Fragments and Local Identity ##
            hithom = {}     # Dictionary of {taxa:homologues} for this occurrence
            locid = {}      # Dictionary of {seq:local %ID} and {taxa:local %ID}
            alnfrag = {qry:qry.info['Sequence'][start:end]} # Partial alignments across occurrence
            seqfrag = {}    # Dictionary of {seq:degapped sequence frag across occurrence}
            for taxa in conseqs.keys():
                hithom[taxa] = []
                locid[taxa] = 0.0
            for seq in seqs:
                igedic = rje_seq.pwIDGapExtra(qry.info['Sequence'][start:end],seq.info['Sequence'][start:end],nomatch=['X'])
                locid[seq] = float(igedic['ID'][0]) / igedic['Len'][0]
                alnfrag[seq] = seq.info['Sequence'][start:end]
                seqfrag[seq] = rje_seq.deGap(seq.info['Sequence'][start:end])
                if (len(seqfrag[seq]) >= 1 or callobj.opt['AlnGap']) and len(seqfrag[seq]) > string.count(seqfrag[seq].upper(),'X'):
                    for taxa in conseqs.keys():
                        if seq in conseqs[taxa]:
                            hithom[taxa].append(seq)
                            locid[taxa] += locid[seq]
            for taxa in conseqs.keys():
                if hithom[taxa]: 
                    locid[taxa] /= len(hithom[taxa])

            ## Reduced alignment for Pos and Prop methods ##
            if callobj.info['ConScore'] in ['all','pos','prop']:
                red_aln = {}    # Dictionary of {seq:sequence to consider}
                for seq in [qry]+seqs:
                    red_aln[seq] = ''
                    m = 0   # Number of positions checked
                    for r in range(len(alnfrag[qry])):
                        a = alnfrag[qry][r]
                        if a == '-':    # Skip
                            continue
                        if alnfrag[seq][r] == '-':
                            red_aln[seq] += '-'
                        elif Occ.getData('Variant')[m] in ['X','*']:    # Not important for match
                            red_aln[seq] += 'X'
                        else:
                            red_aln[seq] += alnfrag[seq][r]
                        m += 1
                    if m != len(Occ.getData('Match')):
                        callobj.log.errorLog('Something wrong with Reduced alignment: %s' % red_aln.values(),printerror=False)
                        raise ValueError
                #X#callobj.deBug(red_aln)
                # >> red_aln sequence should now have an X at each wildcard position, otherwise the relevant residue from the alignment << #

            ## Conservation Scores ##
            for taxa in conseqs.keys():
                hitcons = {'ABS':0.0,'POS':0.0,'PROP':0.0}
                if callobj.info['ConScore'] in ['all','abs']:
                    hitcons['ABS'] = absCons(callobj,Occ,hithom[taxa],seqfrag,seqwt)
                if callobj.info['ConScore'] in ['all','pos']:                        
                    hitcons['POS'] = posCons(callobj,Occ,hithom[taxa],red_aln,seqwt,callobj.dict['PosMatrix'])
                if callobj.info['ConScore'] in ['all','prop']:
                    hitcons['PROP'] = posCons(callobj,Occ,hithom[taxa],red_aln,seqwt,callobj.dict['PropPosMatrix'])

            ## Update Stats ##
                Occ.stat['%s_HOM' % taxa] = len(hithom[taxa])
                Occ.stat['%s_GLOB_ID' % taxa] = globid[taxa]
                Occ.stat['%s_LOC_ID' % taxa] = locid[taxa]
                if callobj.info['ConScore'] == 'all':
                    Occ.stat['%s_CONS' % taxa] = sum(hitcons.values()) / len(hitcons)
                    for method in hitcons.keys():
                        Occ.stat['%s_CONS_%s' % (taxa,method)] = hitcons[method]
                else:
                    Occ.stat['%s_CONS' % taxa] = hitcons[callobj.info['ConScore'].upper()]

    except:
        callobj.log.errorLog('Error in rje_motif_stats.hitAlnCon()',quitchoice=True) 
        return
#########################################################################################################################
def absCons(callobj,Occ,hithom,seqfrag,seqwt):    ### Absolute conservation score.
    '''Absolute conservation score.'''
    try:
        ### Absolute matching of motif in corresponding homologous region ###
        Motif = Occ.obj['Motif']
        hitcon = {}     # Dictionary of {seq:conservation}
        for seq in hithom:
            hitcon[seq] = 0.0
            if callobj.opt['ConsAmb']:   # Search degenerate motif
                vlist = Motif.dict['Search'][0]     
            else:                       # Search with matched variant
                vlist = [Occ.getData('Variant')]
            for variant in vlist:
                searchvar = '(%s)' % string.replace(variant,'X','[A-Z]')
                if rje.matchExp(searchvar,seqfrag[seq]):
                    hitcon[seq] = 1.0
                    break
                    
        ### Weight by distance? ###
        return consWeight(callobj,hitcon,seqwt)
    except:
        callobj.log.errorLog('Error in rje_motif_cons.absCons()',quitchoice=True) 
        return
#########################################################################################################################
def consWeight(callobj,hitcon,seqwt):    ### Weights conservation and returns final conservation score
    '''
    Weights conservation and returns final conservation score.
    >> callobj:Calling object
    >> hitcon: raw dictionary of {seq:conservation}
    >> seqwt: weighting dictionary of {seq:weighting}
    << cons: conservation score
    '''
    try:
        ### Weight conservation of sequences ###
        cons = 0.0
        totweight = 0.0
        for seq in hitcon.keys():
            totweight += seqwt[seq]
            cons += hitcon[seq] * seqwt[seq]
            if hitcon[seq] > 1.0:
                callobj.log.errorLog('HitCon = %.3f > 1 (%s)!' % (hitcon[seq],seq.shortName()),printerror=False)
        if totweight:
            return cons / totweight
        return 0.0
    except:
        callobj.log.errorLog('Error in rje_motif_cons.absCons()',quitchoice=True) 
        return -1.0
#########################################################################################################################
def posCons(callobj,Occ,hithom,red_aln,seqwt,posmatrix):    ### Positional conservation score.
    '''Positional conservation score.'''
    try:
        ### Setup Lists for calculations ###
        infowt = []     # IC at each position for weighting
        ## Setup Information Content Weighting ##
        if callobj.opt['ConsInfo']:
            searchvar = Occ.info['SearchVar'][0:]
        else:
            searchvar = Occ.info['Variant'][0:]
        tvar = searchvar[0:]
        while searchvar:
            if searchvar.find('[A-Z]') == 0:
                el = 'X'
                searchvar = searchvar[5:]
            elif searchvar[0] == '[':
                el = searchvar[1:searchvar.find(']')]
                searchvar = searchvar[searchvar.find(']')+1:]
            else:
                el = searchvar[:1]
                searchvar = searchvar[1:]
            if el not in ['(',')']:
                if callobj.dict['ElementIC'].has_key(el):
                    infowt.append(callobj.dict['ElementIC'][el])
                else:
                    infowt.append(rje_motif_V3.elementIC(el))
        infosum = sum(infowt)
        #X#callobj.deBug('%s varlist = %s' % (hit.info['Name'],varlist))
        ## VarList ##
        varlist = []    # List of AAs at each position for checking (just use hit.info['Variant'] if not callobj.opt['ConsAmb']
        if callobj.opt['ConsAmb']:       
            searchvar = Occ.info['SearchVar'][0:]
            while searchvar:
                if searchvar.find('[A-Z]') == 0:
                    el = 'X'
                    searchvar = searchvar[5:]
                elif searchvar[0] == '[':
                    el = searchvar[1:searchvar.find(']')]
                    searchvar = searchvar[searchvar.find(']')+1:]
                else:
                    el = searchvar[:1]
                    searchvar = searchvar[1:]
                if el not in ['(',')']:
                    varlist.append(el)
        else:       # Search matched variant only
            varlist = Occ.info['Variant']   # Can use "for X in Y" for both lists and strings, so this is OK!
        #X#callobj.deBug('%s varlist = %s => %s' % (hit.info['Name'],varlist,infowt))

        ### Check ###
        if len(infowt) != len(Occ.info['Match']):
            callobj.log.errorLog('Search variant (%s) does not match length of match (%s)' % (tvar,Occ.info['Match']),printerror=False)
            return -1.0

        ### Positional conservation score ###
        hitcon = {}     # Dictionary of {seq:conservation}
        #X#print hithom
        #X#print red_aln
        for seq in hithom:
            hitcon[seq] = 0.0
            #X#print infowt, varlist, red_aln[seq]
            try:
                for p in range(len(varlist)):
                    a = red_aln[seq][p]
                    hitcon[seq] += bestScore(a,varlist[p],posmatrix) * infowt[p]
                hitcon[seq] /= infosum
            except:
                callobj.log.errorLog('Problem with %s:%s posCons(%s) aln "%s" pos %d' % (Occ.hit(),Occ.info['Match'],seq.shortName(),red_aln[seq],p),printerror=True,quitchoice=True)
                hitcon[seq] = 0.0

        ### Weight by distance? ###
        #X#callobj.deBug(hitcon)
        return consWeight(callobj,hitcon,seqwt)
    except:
        callobj.log.errorLog('Error in rje_motif_cons.posCons()',quitchoice=True) 
        return -1.0
#########################################################################################################################
def bestScore(aa,aalist,posmatrix):  ### Best score for compared AAs.
    '''
    >> aa:str = Amino acid to be compared
    >> aalist:str = List of aas to compare to
    >> posmatrix:dict of {'a1a2':score} to get score from
    '''
    if aa == 'X' or aalist == 'A-Z':
        return 0.0
    if aa in aalist:
        return 1.0
    best = 0.0
    for a in aalist:
        if posmatrix.has_key('%s%s' % (aa,a)) and posmatrix['%s%s' % (aa,a)] > best:
            best = posmatrix['%s%s' % (aa,a)]
        elif posmatrix.has_key('%s%s' % (a,aa)) and posmatrix['%s%s' % (a,aa)] > best:
            best = posmatrix['%s%s' % (a,aa)]
    return best
#########################################################################################################################
def loadOrthAln(callobj,seq,gopher=True):    ### Identifies file, loads and checks alignment.
    '''
    Identifies file, loads and checks alignment. If the identified file is not actually aligned, then RJE_SEQ will try to
    align the proteins using MUSCLE or ClustalW.
    >> callobj:Object containing settings for stats generation (MotifList, generally).
    >> seq:Sequence being analysed.
    >> gopher:bool [True] = whether to try to generate alignment with GOPHER if callobj.opt['Gopher']
    << aln = SeqList object containing alignment with queryseq
    '''
    try:
        ### Setup Attributes ###
        v = callobj.stat['Verbose']
        alndir = rje.makePath(callobj.info['AlnDir'])
        alnext = callobj.info['AlnExt']
        
        ### Identify File ###
        if alnext[0] != '.': alnext = '.%s' % alnext
        alnstart = [seq.info['AccNum'],seq.info['ID'],seq.shortName(),None]
        if v > 2: callobj.log.printLog('#PRESTO','%s' % callobj.opt)  #!# Old debugging? #!#
        if callobj.opt['Gopher'] and callobj.opt['FullForce']:
            if v > 0: callobj.log.printLog('#ALN','FullForce=T. Will call Gopher for %s regardless of existing files' % seq.shortName())
            alnstart = [None]
        for file in alnstart:
            if file:
                file = '%s%s%s' % (alndir,file,alnext)
                if rje.checkForFile(file): break  # File found
            else:
                #!# Sort out logging and see if Gopher can be used directly rather than just run() #!#
                ### Run GOPHER ###
                if gopher and callobj.opt['Gopher']:  #!# Add working version for PRESTO and SlimPickings #!#
                    callobj.deBug('Run GOPHER in %s' % callobj.info['GopherDir'])
                    mydir = os.getcwd()
                    os.chdir(callobj.info['GopherDir'])
                    callobj.log.printLog('\n#GOPHER','Running GOPHER on %s' % seq.shortName())
                    try:    #!# Add log.silent() method? #!#
                        gcmd = ['orthtree'] + callobj.cmd_list + ['gnspacc=T','i=-1']
                        solo_gopher = gopher_V2.GopherFork(log=callobj.log,cmd_list=gcmd)
                        solo_gopher.info['Name'] = seq.shortName()
                        solo_gopher.obj['Sequence'] = seq
                        solo_gopher.obj['BLAST'] = gopher_V2.Gopher(callobj.log,gcmd).setupBlast()  #!# Contemplate setting up Gopher in callobj #!#
                        solo_gopher.obj['BLAST'].log = callobj.log
                        solo_gopher.run('orthalign')    #X#gopher_V2.Gopher(callobj.log,gcmd).setMode())
                    except:
                        os.chdir(mydir)
                        callobj.log.errorLog('Problem with Gopher run!')
                        return None
                        
                    if not 'old_school':                            
                        inputseq = 'tmp%s.fas' % rje.randomString(8)
                        TMP = open(inputseq,'w')
                        TMP.write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
                        TMP.close()
                        gcmd = ['orthtree'] + callobj.cmd_list + ['gopher=%s' % inputseq, 'gnspacc=T','i=-1']
                        try:
                            mygopher = gopher_V2.Gopher(log=callobj.log,cmd_list=gcmd)
                            mygopher.run()
                        except:
                            os.chdir(mydir)
                            callobj.log.errorLog('Problem with Gopher run!',printerror=False)
                            return None
                        rje_blast.cleanupDB(callobj,dbfile=inputseq,deletesource=True)
                    os.chdir(mydir)
                if callobj.opt['Gopher']:  
                    file = '%s%s%s' % (alndir,seq.info['AccNum'],alnext)
                    if not os.path.exists(file):
                        file = None
                if not file:
                    callobj.log.printLog('#ALN','No alignment file found for %s in %s.' % (seq.shortName(),alndir),screen=False)
                    return None
        
        ### Load Alignment ###
        callobj.log.stat['Verbose'] = v - 1
        alncmd = ['seqin=None','query=%s' % seq.shortName(),'accnr=F','seqnr=F','autofilter=F','align=T','gnspacc=F'] 
        aln = rje_seq.SeqList(log=callobj.log,cmd_list=callobj.cmd_list+alncmd)
        #X#print file
        aln.loadSeqs(seqfile=file,seqtype='Protein',aln=True,nodup=None)
        callobj.log.stat['Verbose'] = v 
        ## Check Query ##
        qry = aln.obj['QuerySeq']
        if not qry:
            if aln.querySeq(query=seq.info['AccNum']):
                qry = aln.obj['QuerySeq']
            else:
                callobj.log.printLog('#ALN','Problem finding %s in %s.' % (seq.shortName(),file),screen=False)
                return None

        ### Check Alignment ###
        if aln.seqNum() < 2:
            callobj.log.printLog('#ALN','Not enough sequences for %s in %s.' % (seq.shortName(),file),screen=False)
            return None
        if aln._checkAln(aln=True,realign=True):
            return aln
        else:
            callobj.log.printLog('#ERR','%s not aligned!!!' % (file))
            return None       
    except:
        callobj.log.errorLog('Something bad has happened in rje_motif_stats.loadOrthAln()')
        callobj.log.stat['Verbose'] = v 
        return None
#########################################################################################################################
def findOccPos(callobj,Occ,qry,fudge=0):     ### Finds Motif Occurence in alignment
    '''
    Finds Motif Occurence in alignment.
    >> callobj = calling MotifList object
    >> Occ = MotifOcc object
    >> qry = query Sequence object from alignment file
    >> fudge = amount to try shifting match to find occurrence is non-matching sequence
    << (start,end) = start and end position in aligment to allow sequence[start:end]
    '''
    try:
        ### Find Hit in Alignment ###
        (start,end) = (-1,qry.seqLen())             # Start and end positins of match *in alignment*
        qpos = Occ.stat['Pos'] + fudge              # Starting position of hit (from 1->L)
        qmatch = Occ.getData('Match')
        qend = qpos + len(qmatch) - 1    # Ending position of hit (1->L)
        (r,a) = (0,0)                               # Counters for aln residues (r) and amino acid positions (a)
        while r < qry.seqLen():     # Keep looking
            #x#print r, a, qpos, qend, start, end
            if qry.info['Sequence'][r] != '-':      # Not a gap: increment a by 1
                a += 1
            if a == qpos and start < 0:             # Start of match (not yet r+=1 because pos is 1->L)
                start = r
            r += 1                                  # Move on through aligned sequences
            if a == qend:                           # End of match
                end = r
                break

        ### Assess whether hit is right! ###
        amatch = string.replace(qry.info['Sequence'][start:end],'-','')
        if amatch == qmatch:     # Everything is OK!
            return (start,end)
        ## Check whether already fudging! ##
        if fudge != 0:
            raise ValueError

        ### Something is wrong! Try to find real match! ###
        etxt = 'Alignment sequence (%s) does not match occurence (%s)' % (amatch,qmatch)
        #X#if string.replace(qry.info['Sequence'],'-','') == Occ.obj['Seq'].info['Sequence']:  # But sequence matches!
        #X#callobj.log.errorLog('Problem with %s pos %d. Sequences match but %s.' % (qry.shortName(),Occ.stat['Pos'],etxt),printerror=False)
        #X#return (-1,-1)
        ## Try to find match by moving start (using fudge) ##
        if callobj.stat['Interactive'] < 1 or rje.yesNo('%s. Try to find closest correct match?' % etxt):
            fudge = findFudge(string.replace(qry.info['Sequence'],'-',''),qmatch,Occ.stat['Pos']-1)
            if fudge:
                if callobj.stat['Interactive'] > 0: callobj.log.errorLog('%s in alignment differs from input: Fudged %s by %d aa!' % (qry.shortName(),qmatch,fudge),printerror=False)
                else: callobj.log.printLog('#ERR','%s in alignment differs from input: Fudged %s by %d aa!' % (qry.shortName(),qmatch,fudge),screen=False)
                return findOccPos(callobj,Occ,qry,fudge)
            callobj.log.errorLog('%s in alignment differs from input: Cannot find %s anywhere' % (qry.shortName(),qmatch),printerror=False)
            return (-1,-1)
        callobj.log.errorLog(etxt,printerror=False)
        return (-1,-1)
    except:
        callobj.log.errorLog('Something bad has happened in rje_motif_stats.findOccPos()')
        return (-1,-1)
#########################################################################################################################
def findFudge(qryseq,match,pos):    ### Returns fudge to closest match
    '''
    >> qryseq: degapped query sequence
    >> match:match sequence 
    >> pos:position match is meant to be from 0 to L
    << fudge:int = amount to move pos to find match in qryseq. 0 = not there!
    '''
    ### No match! ###
    if qryseq.find(match) < 0:
        return 0
    f = 1
    while (pos - f) >= 0  or (pos + f) < len(qryseq):
        if (pos - f) >= 0 and qryseq[(pos - f):].find(match) == 0:
            return -f
        elif (pos + f) < len(qryseq) and qryseq[(pos + f):].find(match) == 0:
            return f
        f += 1
    raise ValueError
#########################################################################################################################
### END OF SECTION II: Conservation Methods                                                                             #
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
