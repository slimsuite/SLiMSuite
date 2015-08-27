#!/usr/local/bin/python

# GASP - Gapped Ancestral Sequence Prediction
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
Program:      GASP
Description:  Gapped Ancestral Sequence Prediction
Version:      1.4
Last Edit:    16/07/13
Citation:     Edwards & Shields (2004), BMC Bioinformatics 5(1):123. [PMID: 15350199]
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is purely for running the GASP algorithm as a standalone program. All the main functionality is encoded in
    the modules listed below. The GASP algorithm itself is now encoded in the rje_ancseq module. GASP may also be run
    interactively from the rje_tree module.

Commandline: 
    seqin=FILE  : Fasta/ClustalW/Phylip Format sequence file [infile.fas].
    nsfin=FILE  : Newick Standard Format treefile [infile.nsf].

    help    : Triggers Help = list of commandline arguments.

    indeltree=FILE  : Prints a text tree describing indel patterns to FILE.

Uses general modules: os, re, sys, time
Uses RJE modules: rje, rje_ancseq, rje_pam, rje_seq, rje_tree
Additional modules needed: rje_blast, rje_dismatrix, rje_sequence, rje_tree_group, rje_uniprot
"""
#############################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#############################################################################################################################
import os, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
#############################################################################################################################
import rje
import rje_ancseq, rje_seq, rje_tree   # rje_ancseq now sorts out rje_pam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Improved version with second pass.
    # 1.1 - Improved OO. Restriction to descendant AAs. (Good for BAD etc.)
    # 1.2 - No Out Object in Objects
    # 1.3 - Added more interactive load options
    # 1.4 - Minor tweaks to imports.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] Improved Interactive menus for module and major activities - currently use rje_tree.py. Add to rje_ancseq./py?
    # [y] Read in from files and perform GASP
    # [ ] Convert to newer module style.
    # [ ] Take interactive running out of rje_tree and reinstate GASP as top-level module?
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        start_time=time.time()
        program='GASP'
        version="1.4"
        last_edit="July 13"
        description="Gapped Ancestral Sequence Prediction in proteins"
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
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print "\n\nHelp for " + info.program + " V:" + info.version + ": " + time.asctime(time.localtime(info.start_time)) + "\n"
            out.verbose(-1,-1,text=__doc__)
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except KeyboardInterrupt:
        raise
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
    except:
        print "Problem during initial setup."
        raise
#############################################################################################################################
### END OF SECTION I
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                           #
#############################################################################################################################

### Define module-specific methods here
    
#############################################################################################################################
### END OF SECTION III
#############################################################################################################################

                                                    ### ~ ### ~ ###

#############################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                               #
#############################################################################################################################
def runMain():
    try:
        ### <0>  ### Basic Setup of Program
        [info,out,mainlog,cmd_list] = setupProgram()
        
        ### <1> ### Load Data
        ##  <a>  ## Read in Sequences
        try:
            out.verbose(1,3,'Loading sequences...',0)
            seqfile = 'infile.fas'
            nsfin = None
            for cmd in cmd_list:
                if cmd.find('seqin=') == 0:
                    seqfile=cmd[len('seqin='):]
                if cmd.find('nsfin=') == 0:
                    nsfin = cmd[len('nsfin='):]
            basefile = seqfile
            extension = seqfile[-4:]
            if (extension == '.fas') or (extension == '.phy') or (extension == '.aln'):
                basefile = seqfile[:-4]
            seqs = rje_seq.SeqList(log=mainlog,cmd_list=['i=0']+cmd_list+['autofilter=F','autoload=F','seqin=None']) 
            out.verbose(1,3,"from %s" % seqfile,1)
            if not seqs.loadSeqs(seqfile=seqfile,seqtype='protein',aln=True):
                raise
            seqfile = seqs.info['Name']
            basefile = rje.baseFile(seqfile)
            mainlog.printLog('#SEQ',"%s protein sequences read from %s\n" % (str(seqs.seqNum()),seqfile),1)
            mainlog.printLog('#SEQ',"Alignment = %s. (%d aa)\n" % (seqs.opt['Aligned'],seqs.seq[0].seqLen()),1)
        except:
            mainlog.errorLog("Fatal run Exception during Sequence Input\n")
            raise
        ##  <b>  ## Read in Tree
        try:
            if not nsfin:
                nsfin = basefile + '.nsf'
            while not os.path.exists(nsfin):
                if out.stat['Interactive'] >= 0: 
                    nsfin = rje.choice(text='Input tree file "%s" not found. Input filename? (Blank to exit.)' % nsfin)
                    if nsfin == '':
                        raise KeyboardInterrupt
                else:
                    mainlog.log.errorLog('File %s not found. Cannot load tree!' % nsfin,printerror=False,quitchoice=True)
                    raise
            cmd_list.append('nsfin=' + nsfin)
            out.verbose(1,3,'Loading tree from %s...' % nsfin,1)
            mytree = rje_tree.Tree(log=mainlog,cmd_list=['root=yes']+cmd_list)
            mytree.mapSeq(seqlist=seqs)
            mytree.textTree()
            if mytree.opt['ReRooted']:
                mytree.saveTree(filename='%s.nsf' % basefile)
        except KeyboardInterrupt:
            mainlog.errorLog("User terminated.\n")
            raise
        except:
            mainlog.errorLog("Fatal run Exception during Tree Input\n")       
            raise

        ### <2> ### GASP
        try:
            ## <a> ## InDel Tree Setup
            indeltree = None
            for cmd in cmd_list:
                if cmd.find('indeltree=') == 0:
                    indeltree=cmd[len('indeltree='):]

            ## <b> ## GASP
            if indeltree == None or mytree.node[-1].obj['Sequence'] == None:  # Perform GASP
                out.verbose(0,2,'',3)
                mainlog.printLog('#SEQ','GASP: Gapped Ancestral Sequence Prediction',1)
                if basefile == 'infile':
                    basefile = 'gasp'
                mygasp = rje_ancseq.Gasp(tree=mytree,ancfile='%s' % basefile,cmd_list=cmd_list,log=mainlog)
                out.verbose(0,2,'%s' % mygasp.details(),1)
                if out.stat['Interactive'] > 0:
                    if rje.yesNo('Use these parameters?') == False:
                        mygasp.edit()
                mygasp.gasp()
                out.verbose(0,1,"\n\nGASP run completed OK!",2)

            ## <c> ## InDel Tree
            if indeltree:
                mytree.indelTree(filename=indeltree)

        except KeyboardInterrupt:
            mainlog.errorLog("User terminated.\n")
            raise
        except:
            mainlog.errorLog("Fatal run Exception during GASP\n")       
            raise

        ### <X> ### End
    except KeyboardInterrupt:
        mainlog.errorLog("User terminated.\n")
    except:
        print "Unexpected error:", sys.exc_info()[0]
    mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#############################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#############################################################################################################################
### END OF SECTION IV
#############################################################################################################################
